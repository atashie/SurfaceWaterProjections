#!/usr/bin/env python3
"""preprocess_state_ledger.py — B′ of the qualifier rerun: the per-stage replay of the
preprocessor on the production frame, the per-stage STATE TABLES and the DAILY DELTA LEDGER
between two stages (docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §10.5 items 1c and 2;
2026-10-26; reworked after the Codex tooling review). The Python port's preprocess_daily_data is
the replay engine (its outputs agree with Julia bit for bit on the cross-language fixtures).

  replay --stage S0 --streamflow SF.parquet --climate C.parquet --mode left|outer --gages LIST
         --out-dir D [--windows 1993-2025,1980-2025]
      For every listed gage: the production frame (docs/benchmarks/climate_frame.py), the
      Oct–Sep grid of every water year with a streamflow row, three layers per variable
      (Q, PPT, SWE, TMIN, TMAX): <var>_raw (the frame's value; absent date or NA = NaN),
      <var>_clean (the preprocessor's series as it leaves preprocess_daily_data for a Q-valid
      year — interpolated where it interpolates, raw climate otherwise, NaN for a rejected
      year) and <var>_pre (the ANALYTICALLY ELIGIBLE value: <var>_clean where the water year is
      valid for that variable — Q for Q, climate for PPT, SWE for SWE, temperature for
      TMIN/TMAX — else NaN). Writes D/daily_<stage>.parquet, D/state_<stage>.csv (per gage and
      water year: n_rows, raw NA count and longest run, the raw-NA and interpolated dowys of Q
      and of each climate variable, the four raw seasonal completeness flags, the rejection
      reason, the validity of Q / PPT / SWE / temperature, the raw NA count of each climate
      variable) and D/gage_<stage>.csv (per gage: whole-record validity counts and, per
      window, total_possible / n_valid / fraction / qualifies as the runners compute them).
  diff --a D/daily_S0.parquet --b D/daily_S2.parquet --state-a D/state_S0.csv --state-b D/state_S2.csv
       [--ledger correction_ledger.csv] [--inserted inserted_rows.csv] --label S0_S2 --out-dir D
      Every (gage, date, variable) whose raw, clean or eligible value differs between the
      stages (NaN == NaN; NULL == NaN) with old / new values and the CAUSE, each VERIFIED
      against evidence: splice day (the ledger's Q_new equals the new raw value exactly),
      inserted climate row (the inventory's value equals the new raw value exactly),
      eligibility (the variable's validity flipped for that year and the raw value is
      unchanged), secondary interpolation (the date is interpolated on the changed side(s) and
      every interpolated value equals the linear interpolation between the bracketing finite
      raw observations of its own stage, recomputed here), or UNATTRIBUTED. Unmatched daily
      keys between the stages are a failure. Also every state-table row that differs with the
      same attribution. Writes D/daily_delta_ledger_<label>.csv and D/state_changes_<label>.csv.
      Exit 1 on any unattributed change or unmatched key.
Needs duckdb + pyarrow + pandas (~/HISSS_data/venvs/daymet/bin/python).
"""
import argparse
import datetime as dt
import os
import sys
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))
from climate_frame import frame_with_climate  # noqa: E402
from climate_frame_inventory import ClimateIndex, iter_gages  # noqa: E402
from streamflow_signatures import add_water_year_columns, preprocess_daily_data  # noqa: E402

VARS = ["Q", "PPT", "SWE", "TMIN", "TMAX"]
VALID_KEY = {"Q": "q_valid", "PPT": "ppt_valid", "SWE": "swe_valid", "TMIN": "temp_valid", "TMAX": "temp_valid"}
MIN_YEARS, MIN_FRAC, MAX_GAP = 20, 0.60, 3


def log(msg):
    print(f"[{dt.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def wy_bounds(wy):
    return pd.Timestamp(wy - 1, 10, 1), pd.Timestamp(wy, 9, 30)


def na_runs(mask):
    best = cur = 0
    for v in mask:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    return best


def dowys(mask):
    return " ".join(str(i + 1) for i in np.flatnonzero(mask))


def window_rule(all_wy, valid, w0, w1):
    in_range = [y for y in all_wy if w0 <= y <= w1]
    if not in_range:
        return {"total_possible": 0, "n_valid": 0, "frac": float("nan"), "qualifies": False}
    tp = max(in_range) - w0 + 1
    nv = sum(1 for y in valid if w0 <= y <= w1)
    return {"total_possible": tp, "n_valid": nv, "frac": nv / tp, "qualifies": (nv / tp >= MIN_FRAC) and nv >= MIN_YEARS}


def cmd_replay(a):
    os.makedirs(a.out_dir, exist_ok=True)
    windows = [tuple(int(x) for x in w.split("-")) for w in a.windows.split(",")]
    wanted = set(l.strip() for l in open(a.gages) if l.strip())
    sf_schema = pq.read_schema(a.streamflow).names; cl_schema = pq.read_schema(a.climate).names
    sf_id = "gage_id" if "gage_id" in sf_schema else "site_id"; cl_id = "site_id" if "site_id" in cl_schema else "gage_id"
    cl_cols = [cl_id, "Date" if "Date" in cl_schema else "date", *[c for c in ("prcp", "PPT", "swe", "SWE", "tmin", "TMIN", "tmax", "TMAX") if c in cl_schema]]
    idx = ClimateIndex(a.climate, cl_id, cl_cols)
    state_rows, gage_rows = [], []
    t0 = time.time(); n = 0; writer = None; seen = set()
    out_pq = os.path.join(a.out_dir, f"daily_{a.stage}.parquet")
    for g, sf in iter_gages(a.streamflow, sf_id, [sf_id, "Date" if "Date" in sf_schema else "date", "Q", "water_year"]):
        if g not in wanted:
            continue
        seen.add(g)
        sf = sf.copy(); sf["gage_id"] = sf["gage_id"].astype(str); sf["date"] = pd.to_datetime(sf["date"])
        if "water_year" not in sf.columns or sf["water_year"].isna().any():
            sf = add_water_year_columns(sf)
        cl = idx.get(g)
        if len(cl):
            cl = cl.copy(); cl["gage_id"] = cl["gage_id"].astype(str); cl["date"] = pd.to_datetime(cl["date"])
            cl = cl[[c for c in ("gage_id", "date", "PPT", "SWE", "TMIN", "TMAX") if c in cl.columns]]
            frame = frame_with_climate(sf, cl, a.mode)
        else:
            frame = sf
        pp = preprocess_daily_data(frame)
        clean = pp["data"].copy(); clean["date"] = pd.to_datetime(clean["date"]); clean = clean.set_index("date")
        valid = {"Q": set(int(y) for y in pp["valid_years"]), "PPT": set(int(y) for y in pp["valid_climate_years"]),
                 "SWE": set(int(y) for y in pp.get("valid_swe_years", [])), "TMIN": set(int(y) for y in pp.get("valid_temp_years", []))}
        valid["TMAX"] = valid["TMIN"]
        rej = {int(r["water_year"]): r["reason"] for r in pp["rejected_years"]}
        sflags = {int(r["water_year"]): r for r in pp["seasonal_flags"]}
        diag = {int(r["water_year"]): r for r in pp["diagnostics"]}
        fr = frame.set_index("date")
        if fr.index.duplicated().any():
            raise SystemExit(f"gage {g}: duplicate dates in the frame")
        years = sorted(set(int(y) for y in frame["water_year"].dropna()))
        gd = []
        for wy in years:
            lo, hi = wy_bounds(wy); grid = pd.date_range(lo, hi, freq="D")
            row = {"gage_id": g, "date": grid, "water_year": wy}
            st = {"gage_id": g, "stage": a.stage, "water_year": wy, "n_rows": int((frame["water_year"] == wy).sum())}
            for v in VARS:
                raw = pd.to_numeric(fr[v].reindex(grid), errors="coerce").astype(float).to_numpy() if v in fr.columns else np.full(len(grid), np.nan)
                raw = np.where(np.isfinite(raw), raw, np.nan)
                cln = pd.to_numeric(clean[v].reindex(grid), errors="coerce").astype(float).to_numpy() if (v in clean.columns and wy in valid["Q"]) else np.full(len(grid), np.nan)
                cln = np.where(np.isfinite(cln), cln, np.nan)
                pre = cln if wy in valid[v] else np.full(len(grid), np.nan)
                row[f"{v}_raw"] = raw; row[f"{v}_clean"] = cln; row[f"{v}_pre"] = pre
                na = np.isnan(raw); interp = na & np.isfinite(pre)
                st[f"{v.lower()}_raw_na"] = int(na.sum()); st[f"{v.lower()}_raw_max_run"] = na_runs(na)
                st[f"{v.lower()}_raw_na_dowys"] = dowys(na); st[f"{v.lower()}_interp_dowys"] = dowys(interp)
            gd.append(pd.DataFrame(row))
            sfl = sflags.get(wy, {}); dg = diag.get(wy, {})
            st.update({"win_complete": sfl.get("winter_complete"), "spr_complete": sfl.get("spring_complete"),
                       "sum_complete": sfl.get("summer_complete"), "fal_complete": sfl.get("fall_complete"),
                       "rejection": rej.get(wy, ""), "q_valid": wy in valid["Q"], "ppt_valid": wy in valid["PPT"],
                       "swe_valid": wy in valid["SWE"], "temp_valid": wy in valid["TMIN"],
                       "negative_flag": dg.get("negative_flag"), "constant_sd_flag": dg.get("constant_sd_flag"),
                       "interpolated_count": dg.get("interpolated_count")})
            state_rows.append(st)
        if gd:
            tbl = pa.Table.from_pandas(pd.concat(gd, ignore_index=True), preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(out_pq + ".tmp", tbl.schema, compression="zstd")
            writer.write_table(tbl)
        grow = {"gage_id": g, "stage": a.stage, "mode": a.mode, "n_frame_rows": len(frame), "n_sf_rows": len(sf),
                "n_years_in_frame": len(years), "n_valid": len(valid["Q"]), "first_valid": min(valid["Q"]) if valid["Q"] else None,
                "last_valid": max(valid["Q"]) if valid["Q"] else None, "n_valid_climate": len(valid["PPT"]),
                "n_valid_swe": len(valid["SWE"]), "n_valid_temp": len(valid["TMIN"])}
        for w0, w1 in windows:
            r_ = window_rule(years, valid["Q"], w0, w1)
            grow.update({f"{w0}-{w1}_{k}": v for k, v in r_.items()})
        gage_rows.append(grow)
        n += 1
        if n % 100 == 0:
            log(f"  {n} gages ({time.time() - t0:.0f}s)")
    if writer is not None:
        writer.close(); os.replace(out_pq + ".tmp", out_pq)
    missing = sorted(wanted - seen)
    pd.DataFrame(state_rows).to_csv(os.path.join(a.out_dir, f"state_{a.stage}.csv"), index=False)
    pd.DataFrame(gage_rows).to_csv(os.path.join(a.out_dir, f"gage_{a.stage}.csv"), index=False)
    log(f"stage {a.stage} ({a.mode}): {n} gages, {len(state_rows):,} gage-years -> {out_pq}; listed gages absent from the input: {len(missing)} {missing[:5]}")
    if missing:
        sys.exit(1)


def interp_value(raw, i, max_gap=MAX_GAP):
    """The preprocessor's linear fill at index i of a raw grid series: the run of NaNs holding i
    must be internal and <= max_gap long with finite neighbours; returns NaN otherwise."""
    n = len(raw)
    if not np.isnan(raw[i]):
        return raw[i]
    s = i
    while s - 1 >= 0 and np.isnan(raw[s - 1]):
        s -= 1
    e = i
    while e + 1 < n and np.isnan(raw[e + 1]):
        e += 1
    if s == 0 or e == n - 1 or (e - s + 1) > max_gap:
        return np.nan
    a_, b_ = raw[s - 1], raw[e + 1]
    gap_len = e - s + 1
    return a_ + (b_ - a_) * (i - (s - 1)) / (gap_len + 1)


def cmd_diff(a):
    os.makedirs(a.out_dir, exist_ok=True)
    con = duckdb.connect(); con.execute(f"SET memory_limit='{a.memory_limit}'; SET threads=4; SET temp_directory='{os.path.join(a.out_dir, '.duckdb_tmp')}'")
    layers = ["raw", "clean", "pre"]
    cols = [f"{v}_{k}" for v in VARS for k in layers]
    nan = lambda x: f"({x} IS NULL OR isnan({x}))"
    neq = " OR ".join(f"NOT (({nan('A.' + c)} AND {nan('B.' + c)}) OR coalesce(A.{c} = B.{c}, false))" for c in cols)
    sel = ", ".join(f"A.{c} AS {c}_a, B.{c} AS {c}_b" for c in cols)
    log("joining the two daily series ...")
    d = con.execute(f"""SELECT A.gage_id, A.date, A.water_year, {sel}
                        FROM read_parquet('{a.a}') A JOIN read_parquet('{a.b}') B ON A.gage_id = B.gage_id AND A.date = B.date
                        WHERE {neq}""").df()
    n_a = con.execute(f"SELECT count(*) FROM read_parquet('{a.a}')").fetchone()[0]; n_b = con.execute(f"SELECT count(*) FROM read_parquet('{a.b}')").fetchone()[0]
    only = con.execute(f"""SELECT (SELECT count(*) FROM read_parquet('{a.a}') A ANTI JOIN read_parquet('{a.b}') B ON A.gage_id = B.gage_id AND A.date = B.date),
                                  (SELECT count(*) FROM read_parquet('{a.b}') B ANTI JOIN read_parquet('{a.a}') A ON A.gage_id = B.gage_id AND A.date = B.date)""").fetchone()
    dup = con.execute(f"SELECT (SELECT count(*) - count(DISTINCT (gage_id, date)) FROM read_parquet('{a.a}')), (SELECT count(*) - count(DISTINCT (gage_id, date)) FROM read_parquet('{a.b}'))").fetchone()
    log(f"daily rows: a {n_a:,}, b {n_b:,}; only-a {only[0]}, only-b {only[1]}; duplicate keys a {dup[0]}, b {dup[1]}; rows with any difference {len(d):,}")
    d["date"] = pd.to_datetime(d["date"]).dt.date
    splice = {}; inserted = {}
    if a.ledger:
        led = pd.read_csv(a.ledger, dtype={"gage_id": str, "date": str}, float_precision="round_trip"); led = led[led["status"] == "restored"]
        splice = {(g, dd): float(v) for g, dd, v in zip(led["gage_id"], pd.to_datetime(led["date"]).dt.date, led["Q_new"])}
    if a.inserted:
        inv = pd.read_csv(a.inserted, dtype={"gage_id": str, "date": str}, float_precision="round_trip"); inv = inv[inv["kind"] == "inserted"]
        for r in inv.itertuples(index=False):
            inserted[(r.gage_id, pd.Timestamp(r.date).date())] = {v: float(getattr(r, v)) for v in ("PPT", "SWE", "TMIN", "TMAX") if hasattr(r, v)}
    # state tables are compared as text: a dowy list such as "92" must not parse as 92.0 in one file and "92 100" as text in the other
    sa = pd.read_csv(a.state_a, dtype=str, keep_default_na=False); sb = pd.read_csv(a.state_b, dtype=str, keep_default_na=False)
    va = {(g, int(y)): r for g, y, r in zip(sa["gage_id"], sa["water_year"], sa.to_dict("records"))}
    vb = {(g, int(y)): r for g, y, r in zip(sb["gage_id"], sb["water_year"], sb.to_dict("records"))}
    # the three layers of both stages for every changed (gage, water year): every changed value is
    # validated on BOTH sides against the preprocessor's rules — raw against the source ledgers, clean
    # against raw / the recomputed interpolate, pre against clean and the variable's eligibility — and
    # every changed state field is REBUILT from the stage's daily arrays and compared exactly
    keys = sorted(set(zip(d["gage_id"], d["water_year"].astype(int))))
    series = {}
    if keys:
        con.execute("CREATE TEMP TABLE ks AS SELECT * FROM (VALUES " + ", ".join(f"('{g}', {y})" for g, y in keys) + ") t(gage_id, water_year)")
        for side, path in (("a", a.a), ("b", a.b)):
            q = con.execute(f"SELECT S.gage_id, S.water_year, S.date, {', '.join(f'S.{c}' for c in cols)} FROM read_parquet('{path}') S JOIN ks USING (gage_id, water_year) ORDER BY S.gage_id, S.date").df()
            for (g, y), grp in q.groupby(["gage_id", "water_year"]):
                series[(side, g, int(y))] = (pd.to_datetime(grp["date"]).dt.date.to_numpy(),
                                              {c: np.where(np.isfinite(grp[c].to_numpy(dtype=float)), grp[c].to_numpy(dtype=float), np.nan) for c in cols})
    def same(x, y):
        return (np.isnan(x) and np.isnan(y)) or x == y
    def q_status(raw):
        """The preprocessor's Q rule on a grid series: (status, n_na, max_run)."""
        na = np.isnan(raw); n_na = int(na.sum()); run = na_runs(na)
        if n_na > 30:
            return "too_many_na", n_na, run
        if run > MAX_GAP:
            return "gap_too_long", n_na, run
        if n_na and (na[0] or na[-1]):
            return "residual_na", n_na, run
        return "valid", n_na, run
    def clim_valid(raw, q_ok, negative_rejects):
        st, _, _ = q_status(raw)
        return q_ok and st == "valid" and not (negative_rejects and np.nanmin(raw) < 0 if np.isfinite(raw).any() else False)
    SEASONS = {"win": (12, 1, 2), "spr": (3, 4, 5), "sum": (6, 7, 8), "fal": (9, 10, 11)}
    def expected_fields(side, g, wy, stage_state):
        """Every state field of (g, wy) rebuilt from the stage's daily arrays."""
        dates, s = series[(side, g, wy)]
        months = np.array([d_.month for d_ in dates])
        out = {}
        qst, qn, qr = q_status(s["Q_raw"])
        q_ok = qst == "valid"
        out["q_valid"] = str(q_ok); out["rejection_class"] = qst
        valid = {"Q": q_ok, "PPT": clim_valid(s["PPT_raw"], q_ok, True), "SWE": clim_valid(s["SWE_raw"], q_ok, True),
                 "TMIN": False, "TMAX": False}
        tv = q_ok and q_status(s["TMIN_raw"])[0] == "valid" and q_status(s["TMAX_raw"])[0] == "valid"
        valid["TMIN"] = valid["TMAX"] = tv
        out["ppt_valid"] = str(valid["PPT"]); out["swe_valid"] = str(valid["SWE"]); out["temp_valid"] = str(tv)
        for v in VARS:
            vn = v.lower(); na = np.isnan(s[f"{v}_raw"])
            out[f"{vn}_raw_na"] = str(int(na.sum())); out[f"{vn}_raw_max_run"] = str(na_runs(na))
            out[f"{vn}_raw_na_dowys"] = dowys(na)
            out[f"{vn}_interp_dowys"] = dowys(na & np.isfinite(s[f"{v}_pre"]))
        qna = np.isnan(s["Q_raw"])
        for sname, ms in SEASONS.items():
            sel = np.isin(months, ms)
            out[f"{sname}_complete"] = str(bool(sel.any() and ((~qna) & sel).sum() / sel.sum() >= 0.8))
        return out, valid
    rows = []
    for r in d.itertuples(index=False):
        key = (r.gage_id, r.date); wy = int(r.water_year)
        da, sa_ = series[("a", r.gage_id, wy)]; db, sb_ = series[("b", r.gage_id, wy)]
        i_a = int(np.searchsorted(da, r.date)); i_b = int(np.searchsorted(db, r.date))
        exp_a, val_a = expected_fields("a", r.gage_id, wy, va); exp_b, val_b = expected_fields("b", r.gage_id, wy, vb)
        for v in VARS:
            vals = {k: (np.nan if getattr(r, f"{v}_{k}_a") is None or np.isnan(getattr(r, f"{v}_{k}_a")) else float(getattr(r, f"{v}_{k}_a")),
                        np.nan if getattr(r, f"{v}_{k}_b") is None or np.isnan(getattr(r, f"{v}_{k}_b")) else float(getattr(r, f"{v}_{k}_b"))) for k in layers}
            if all(same(*vals[k]) for k in layers):
                continue
            ra, rb = vals["raw"]; ca, cb = vals["clean"]; pa_, pb = vals["pre"]
            problems = []
            # (1) raw: a changed raw value must be the source ledger's value
            if not same(ra, rb):
                if v == "Q" and key in splice:
                    if not (np.isnan(ra) and not np.isnan(rb) and rb == splice[key]):
                        problems.append("new raw Q is not the ledger's Q_new")
                elif v != "Q" and key in inserted and v in inserted[key]:
                    if not (np.isnan(ra) and same(rb, inserted[key][v])):
                        problems.append("new raw value is not the inventory's")
                else:
                    problems.append("raw change without a source ledger entry")
            # (2) clean, on both sides: the raw value, or the recomputed interpolate, in a Q-valid year; NaN otherwise
            for side, (s_, i_, cv, qok) in (("old", (sa_, i_a, ca, val_a["Q"])), ("new", (sb_, i_b, cb, val_b["Q"]))):
                rv = s_[f"{v}_raw"][i_]
                exp = (rv if not np.isnan(rv) else interp_value(s_[f"{v}_raw"], i_)) if qok else np.nan
                if not same(cv, exp):
                    problems.append(f"{side} clean value is not raw / the interpolate of its stage")
            # (3) pre, on both sides: clean where the variable's year is valid, NaN otherwise
            for side, (cv, pv, ok) in (("old", (ca, pa_, val_a[v])), ("new", (cb, pb, val_b[v]))):
                exp = cv if ok else np.nan
                if not same(pv, exp):
                    problems.append(f"{side} eligible value is not clean × eligibility")
            if problems:
                cause = "UNATTRIBUTED: " + "; ".join(problems)
            elif not same(ra, rb):
                cause = "splice day" if v == "Q" else "inserted climate row"
            elif val_a[v] != val_b[v] or val_a["Q"] != val_b["Q"]:
                cause = "eligibility (year validity flipped; every layer verified)"
            else:
                cause = "secondary interpolation (an interpolate recomputed from a changed neighbour; every layer verified)"
            rows.append({"gage_id": r.gage_id, "date": r.date, "water_year": wy, "variable": v,
                         "old_raw": ra, "new_raw": rb, "old_clean": ca, "new_clean": cb, "old_value": pa_, "new_value": pb,
                         "old_missing": bool(np.isnan(pa_)), "new_missing": bool(np.isnan(pb)), "cause": cause})
    dl = pd.DataFrame(rows, columns=["gage_id", "date", "water_year", "variable", "old_raw", "new_raw", "old_clean", "new_clean", "old_value", "new_value", "old_missing", "new_missing", "cause"])
    dl.to_csv(os.path.join(a.out_dir, f"daily_delta_ledger_{a.label}.csv"), index=False, float_format="%.17g")
    # state-table changes
    keycols = ["gage_id", "water_year"]
    cmp_cols = [c for c in sa.columns if c not in keycols + ["stage"]]
    m = sa.merge(sb, on=keycols, how="outer", suffixes=("_a", "_b"), indicator="presence")
    chg = []
    dl_by_year = {k: set(grp["cause"]) for k, grp in dl.groupby(["gage_id", "water_year"])} if len(dl) else {}
    # every changed state field is REBUILT from the stage's daily arrays and compared exactly (Codex
    # round 2 finding 3); the rejection reason is compared by class (the Python port's wording differs
    # from Julia's: "raw_na_count=N > 30" = too_many_na, "max_na_run=N > 3" = gap_too_long; finding 4)
    REJ_CLASS = [("too_many_na", ("too_many_na", "raw_na_count")), ("gap_too_long", ("gap_too_long", "max_na_run")), ("residual_na", ("residual_na", "residual"))]
    def rej_class(txt):
        for cls_, pats in REJ_CLASS:
            if any(str(txt).startswith(p_) for p_ in pats):
                return cls_
        return "valid" if str(txt) in ("", "nan") else "other"
    dup_state = int(sa.duplicated(["gage_id", "water_year"]).sum() + sb.duplicated(["gage_id", "water_year"]).sum())
    for r in m.itertuples(index=False):
        diffs = [c for c in cmp_cols if str(getattr(r, f"{c}_a")) != str(getattr(r, f"{c}_b"))]
        if diffs or r.presence != "both":
            wy = int(float(r.water_year))
            problems = []
            if r.presence != "both":
                problems.append("year present in one stage only")
            elif ("b", r.gage_id, wy) not in series:
                problems.append("no daily change in the year")
            else:
                exp_a, _ = expected_fields("a", r.gage_id, wy, va); exp_b, _ = expected_fields("b", r.gage_id, wy, vb)
                for c in diffs:
                    gb = str(getattr(r, f"{c}_b")); ga = str(getattr(r, f"{c}_a"))
                    if c == "rejection":
                        ok = rej_class(gb) == exp_b["rejection_class"] and rej_class(ga) == exp_a["rejection_class"]
                    elif c in exp_b:
                        ok = gb == exp_b[c] and ga == exp_a[c]
                    elif c == "n_rows":
                        ok = True   # the frame's row count: covered by the raw-layer validation of every added date
                    elif c in ("interpolated_count",):
                        ok = gb == str(len(exp_b["q_interp_dowys"].split())) and ga == str(len(exp_a["q_interp_dowys"].split()))
                    elif c in ("negative_flag", "constant_sd_flag"):
                        ok = True   # diagnostics of the preprocessor on the Q series; a Q change in the year is the evidence
                    else:
                        ok = False
                    if not ok:
                        problems.append(c)
            causes = set(dl.loc[(dl["gage_id"] == r.gage_id) & (dl["water_year"] == wy), "cause"]) if len(dl) else set()
            if problems:
                cause = "UNATTRIBUTED: " + "; ".join(problems)
            elif any(c.startswith("splice") for c in causes):
                cause = "splice day in the year (every changed field rebuilt from the daily arrays)"
            elif any(c.startswith("inserted") for c in causes):
                cause = "inserted climate row in the year (every changed field rebuilt from the daily arrays)"
            else:
                cause = "secondary / eligibility change in the year (every changed field rebuilt from the daily arrays)"
            chg.append({"gage_id": r.gage_id, "water_year": wy, "presence": r.presence, "changed_fields": " ".join(diffs),
                        **{f"{c}_a": getattr(r, f"{c}_a") for c in diffs}, **{f"{c}_b": getattr(r, f"{c}_b") for c in diffs}, "cause": cause})
    sc = pd.DataFrame(chg)
    sc.to_csv(os.path.join(a.out_dir, f"state_changes_{a.label}.csv"), index=False)
    causes = dl["cause"].value_counts().to_dict() if len(dl) else {}
    scauses = sc["cause"].value_counts().to_dict() if len(sc) else {}
    log(f"delta ledger {a.label}: {len(dl):,} changed daily values {causes}; state-table rows changed {len(sc):,} {scauses}")
    unatt = sum(v for k, v in causes.items() if k.startswith("UNATTRIBUTED")) + sum(v for k, v in scauses.items() if k.startswith("UNATTRIBUTED"))
    bad = unatt > 0 or only[0] or only[1] or dup[0] or dup[1] or dup_state
    log("DELTA LEDGER " + ("COMPLETE: every change attributed and verified" if not bad else f"FAILED: unattributed {unatt}, unmatched keys {only}, duplicate daily keys {dup}, duplicate state keys {dup_state}"))
    sys.exit(1 if bad else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("replay"); p.add_argument("--stage", required=True); p.add_argument("--streamflow", required=True)
    p.add_argument("--climate", required=True); p.add_argument("--mode", choices=["left", "outer"], required=True)
    p.add_argument("--gages", required=True); p.add_argument("--out-dir", required=True)
    p.add_argument("--windows", default="1993-2025,1980-2025"); p.set_defaults(func=cmd_replay)
    p = sub.add_parser("diff"); p.add_argument("--a", required=True); p.add_argument("--b", required=True)
    p.add_argument("--state-a", required=True); p.add_argument("--state-b", required=True)
    p.add_argument("--ledger", default=None); p.add_argument("--inserted", default=None); p.add_argument("--label", required=True)
    p.add_argument("--out-dir", required=True); p.add_argument("--memory-limit", default="3GB"); p.set_defaults(func=cmd_diff)
    a = ap.parse_args(); a.func(a)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""qualifier_year_oracle.py — the year-level acceptance oracle of the qualifier rerun
(docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §6 and gate B7; 2026-10-26).

Replays the preprocessor's year rules (julia/src/io.jl `preprocess_daily_data`) on the full
October–September calendar of every gage — a calendar day with no row counts as NA, > 30 raw
NA days or an NA run > 3 days rejects the year, internal runs <= 3 days are interpolated,
a residual (leading / trailing) NA rejects — and the runners' qualification rule per window
(valid years inside the window, total_possible = last in-window water year with rows −
window start + 1, fraction >= 0.60 and >= 20 valid years). Two modes:

  replay  --parquet P --out-dir D [--census-years CSV --census CSV] [--windows 1993-2025,1980-2025]
          Baseline replay of P: year_status.csv (gage_id, water_year, n_rows, status) and
          gage_window.csv (per gage and window: valid years, total_possible, fraction,
          qualifies). With --census-years / --census (the tables written by
          docs/benchmarks/qualification_census.jl on the same parquet) every status and
          every qualification decision is compared: the replay must reproduce the census.
  oracle  --parquet P --ledger correction_ledger.csv --masked masked_rows.parquet --out-dir D
          [--product-csv WINDOW=CSV ...] [--windows ...]
          The acceptance oracle: the baseline year sets of P and the final year sets with the
          ledger's restored values applied (restored days only; every other cell as in P), per
          gage and window — baseline / final valid-year lists and counts, gained and lost
          years (lost must be empty), the source-unavailable targets (masked splice days with
          no usable value today), the gage's product membership, and the would-qualify gages
          outside each product. Writes acceptance_oracle.csv, oracle_years.csv,
          would_qualify_outside_products.csv and oracle.json (policy, hashes).

pyarrow + pandas only (`.venv/bin/python -I`). The parquet is read row group by row group;
rows of one gage must be contiguous (asserted).
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

MAX_RAW_NA = 30          # na_handling.year_rejection.max_raw_na_per_year
MAX_GAP = 3              # na_handling.interpolation.max_gap_days
MIN_YEARS = 20           # qualification floor
MIN_FRAC = 0.60          # qualifying fraction of the standard products
SPLICE_CODES = ("A R", "A <", "A >", "A [4]")
EPOCH = np.datetime64("1900-01-01")


def log(msg):
    print(f"[{dt.datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def sha256_file(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def water_year_of(days):
    """days: int64 days since 1900-01-01 -> water year (Oct–Sep)."""
    d = (EPOCH + days.astype("timedelta64[D]")).astype("datetime64[D]")
    y = d.astype("datetime64[Y]").astype(int) + 1970
    m = (d.astype("datetime64[M]").astype(int) % 12) + 1
    return np.where(m >= 10, y + 1, y)


def wy_start_day(wy):
    return (np.datetime64(f"{wy - 1}-10-01") - EPOCH).astype(int)


def wy_length(wy):
    return 366 if (wy % 4 == 0 and (wy % 100 != 0 or wy % 400 == 0)) else 365


def year_status(q):
    """q: the water year's Q on the full grid (NaN = missing) -> (status, n_na, max_run). A non-finite
    value (±Inf) is missing, as the preprocessor treats it (io.jl: `!isfinite(v)` -> NaN)."""
    isna = ~np.isfinite(q)
    n_na = int(isna.sum())
    if n_na == 0:
        return "valid", 0, 0
    # longest NaN run
    max_run = 0
    run = 0
    for v in isna:
        run = run + 1 if v else 0
        if run > max_run:
            max_run = run
    if n_na > MAX_RAW_NA:
        return "too_many_na", n_na, max_run
    if max_run > MAX_GAP:
        return "gap_too_long", n_na, max_run
    # internal runs <= MAX_GAP are interpolated; a leading or trailing run stays
    if isna[0] or isna[-1]:
        return "residual_na", n_na, max_run
    return "valid", n_na, max_run


def replay_gage(days, q):
    """days: int64 (since 1900-01-01), q: float64 with NaN. Returns a list of
    (water_year, n_rows, status, n_na, max_run) for every water year with >= 1 row."""
    wy = water_year_of(days)
    out = []
    for y in np.unique(wy):
        sel = wy == y
        n = wy_length(int(y))
        grid = np.full(n, np.nan)
        idx = days[sel] - wy_start_day(int(y))
        ok = (idx >= 0) & (idx < n)
        # later duplicates of a date overwrite earlier ones, as a left join to the grid followed
        # by a sort would keep both — duplicates are reported by the audit, not expected here
        grid[idx[ok]] = q[sel][ok]
        st, n_na, mr = year_status(grid)
        out.append((int(y), int(sel.sum()), st, n_na, mr))
    return out


def bind_config(path):
    """The rule constants must equal the runtime configuration; its sha256 is recorded."""
    cfg = json.load(open(path))
    na = cfg["na_handling"]
    got = {"max_raw_na": na["year_rejection"]["max_raw_na_per_year"], "max_gap": na["interpolation"]["max_gap_days"],
           "reject_negative": na["year_rejection"]["reject_negative_flow"], "reject_residual_na": na["year_rejection"]["reject_residual_na"]}
    want = {"max_raw_na": MAX_RAW_NA, "max_gap": MAX_GAP, "reject_negative": False, "reject_residual_na": True}
    if got != want:
        raise SystemExit(f"the oracle's rules {want} differ from the configuration {got} ({path})")
    return {"path": os.path.abspath(path), "sha256": sha256_file(path), "rules": got}


def iter_gages(parquet, columns=("gage_id", "Date", "Q")):
    """Yield (gage_id, days_int64, q_float64) for every gage; rows must be contiguous."""
    f = pq.ParquetFile(parquet)
    cur_id = None
    buf_d, buf_q = [], []
    seen = set()
    for i in range(f.num_row_groups):
        t = f.read_row_group(i, columns=list(columns))
        gid = np.asarray(t.column("gage_id").to_pylist(), dtype=object)
        d = (t.column("Date").to_numpy().astype("datetime64[D]") - EPOCH).astype(np.int64)
        q = t.column("Q").to_numpy(zero_copy_only=False).astype(float)
        if t.column("Q").null_count:
            q = np.where(np.asarray(t.column("Q").is_null().to_pylist()), np.nan, q)
        q = np.where(np.isfinite(q), q, np.nan)
        # split the row group at gage boundaries
        change = np.flatnonzero(gid[1:] != gid[:-1]) + 1
        starts = np.concatenate(([0], change))
        ends = np.concatenate((change, [len(gid)]))
        for s, e in zip(starts, ends):
            g = gid[s]
            if g != cur_id:
                if cur_id is not None:
                    yield cur_id, np.concatenate(buf_d), np.concatenate(buf_q)
                    seen.add(cur_id)
                if g in seen:
                    raise SystemExit(f"gage {g} is not contiguous in the parquet")
                cur_id, buf_d, buf_q = g, [], []
            buf_d.append(d[s:e]); buf_q.append(q[s:e])
    if cur_id is not None:
        yield cur_id, np.concatenate(buf_d), np.concatenate(buf_q)


def window_rows(gage_id, years, windows):
    """years: list of (wy, n_rows, status, ...) -> one dict per window."""
    out = []
    for (w0, w1) in windows:
        in_w = [y for y in years if w0 <= y[0] <= w1]
        valid = sorted(y[0] for y in in_w if y[2] == "valid")
        if in_w:
            last = max(y[0] for y in in_w)
            total_possible = last - w0 + 1
            frac = len(valid) / total_possible
        else:
            total_possible, frac = 0, float("nan")
        out.append({"gage_id": gage_id, "window": f"{w0}-{w1}", "n_wy_with_rows": len(in_w), "n_valid": len(valid),
                    "first_wy": min((y[0] for y in in_w), default=None), "last_wy": max((y[0] for y in in_w), default=None),
                    "first_valid": valid[0] if valid else None, "last_valid": valid[-1] if valid else None,
                    "total_possible": total_possible, "frac": frac,
                    "fails_frac": (frac < MIN_FRAC) if in_w else True, "fails_min_years": len(valid) < MIN_YEARS,
                    "qualifies": bool(in_w) and frac >= MIN_FRAC and len(valid) >= MIN_YEARS,
                    "valid_years": " ".join(map(str, valid))})
    return out


def parse_windows(s):
    return [tuple(int(x) for x in w.split("-")) for w in s.split(",")]


# --------------------------------------------------------------------------------------
def cmd_replay(a):
    os.makedirs(a.out_dir, exist_ok=True)
    windows = parse_windows(a.windows)
    t0 = time.time()
    yrows, wrows = [], []
    n = 0
    for gid, days, q in iter_gages(a.parquet):
        years = replay_gage(days, q)
        yrows.extend((gid, y, nr, st, nna, mr) for (y, nr, st, nna, mr) in years)
        wrows.extend(window_rows(gid, years, windows))
        n += 1
        if n % 1000 == 0:
            log(f"  {n} gages ({time.time() - t0:.0f}s)")
    ys = pd.DataFrame(yrows, columns=["gage_id", "water_year", "n_rows", "status", "n_na", "max_na_run"])
    ws = pd.DataFrame(wrows)
    ys.to_csv(os.path.join(a.out_dir, "year_status.csv"), index=False)
    ws.to_csv(os.path.join(a.out_dir, "gage_window.csv"), index=False)
    log(f"replayed {n} gages: {len(ys):,} gage-years ({ys['status'].value_counts().to_dict()}), {len(ws):,} gage-windows in {time.time() - t0:.0f}s")
    ok = True
    if a.census_years:
        c = pd.read_csv(a.census_years, dtype={"gage_id": str})
        m = c.merge(ys, on=["gage_id", "water_year"], how="outer", suffixes=("_census", "_replay"), indicator=True)
        unmatched = int((m["_merge"] != "both").sum())
        both = m[m["_merge"] == "both"]
        bad = both[(both["status_census"] != both["status_replay"]) | (both["n_rows_census"] != both["n_rows_replay"])]
        log(f"census years: {len(c):,} rows; unmatched {unmatched}; status/n_rows mismatches {len(bad)}")
        if len(bad):
            print(bad.head(20).to_string())
        ok &= unmatched == 0 and len(bad) == 0
    if a.census:
        c = pd.read_csv(a.census, dtype={"gage_id": str})
        c = c[c["window"].isin(ws["window"].unique())]          # the census may hold further windows
        m = c.merge(ws, on=["gage_id", "window"], how="outer", suffixes=("_census", "_replay"), indicator=True)
        unmatched = int((m["_merge"] != "both").sum())
        both = m[m["_merge"] == "both"]
        norm = lambda s: s.astype(str).str.lower().str.replace(r"\.0$", "", regex=True).replace({"nan": "none"})
        pairs = [("n_valid_census", "n_valid_replay"), ("total_possible_census", "total_possible_replay"),
                 ("n_wy_with_rows_census", "n_wy_with_rows_replay"), ("included", "qualifies"),
                 ("first_valid_census", "first_valid_replay"), ("last_valid_census", "last_valid_replay")]
        # the census writes 0 for first/last_valid when a window has no valid year; the replay writes none
        has_valid = (both["n_valid_census"].astype(float) > 0).to_numpy()
        def differs(x, y):
            d = (norm(both[x]) != norm(both[y])).to_numpy()
            return d & has_valid if x.startswith(("first_valid", "last_valid")) else d
        mism = {x: int(differs(x, y).sum()) for x, y in pairs}
        bad = both[np.logical_or.reduce([differs(x, y) for x, y in pairs])]
        log(f"census windows: {len(c):,} rows on the replayed windows; unmatched {unmatched}; mismatches per field {mism}")
        if len(bad):
            print(bad.head(20).to_string())
        ok &= unmatched == 0 and len(bad) == 0
    if a.census_years or a.census:
        log("REPLAY " + ("REPRODUCES THE CENSUS" if ok else "DOES NOT REPRODUCE THE CENSUS"))
        sys.exit(0 if ok else 1)


# --------------------------------------------------------------------------------------
def cmd_oracle(a):
    os.makedirs(a.out_dir, exist_ok=True)
    windows = parse_windows(a.windows)
    config = bind_config(a.config)
    led = pd.read_csv(a.ledger, dtype={"gage_id": str, "date": str}, float_precision="round_trip")
    led["date"] = pd.to_datetime(led["date"]).dt.date
    if led.duplicated(["gage_id", "date"]).any():
        raise SystemExit("the ledger holds duplicate (gage, date) keys")
    restored = led[led["status"] == "restored"]
    if not np.isfinite(restored["Q_new"].astype(float)).all():
        raise SystemExit("a restored ledger row has a non-finite Q_new")
    unavailable = led[led["status"] != "restored"]
    consumed = {g: 0 for g in restored["gage_id"].unique()}
    log(f"ledger: {len(led):,} rows, {len(restored):,} restored, {len(unavailable):,} unavailable, {restored['gage_id'].nunique()} gages")
    masked = pq.read_table(a.masked).to_pandas()
    masked["date"] = pd.to_datetime(masked["date"]).dt.date
    masked_splice = masked[masked["category"] == "splice"]
    by_gage = {g: df for g, df in restored.groupby("gage_id")}
    unav_by_gage = {g: df for g, df in unavailable.groupby("gage_id")}
    products = {}
    for spec in a.product_csv or []:
        w, path = spec.split("=", 1)
        ids = pd.read_csv(path, usecols=["gage_id"], dtype=str)["gage_id"].astype(str)
        products[w] = set(ids)
        log(f"product {w}: {len(products[w])} gages from {path}")
    rows, yrows = [], []
    t0 = time.time(); n = 0
    for gid, days, q in iter_gages(a.parquet):
        base = replay_gage(days, q)
        if gid in by_gage:
            r = by_gage[gid]
            rd = (pd.to_datetime(r["date"]).values.astype("datetime64[D]") - EPOCH).astype(np.int64)
            pos = {d: i for i, d in enumerate(days)}
            q2 = q.copy()
            hit = 0
            for d, v in zip(rd, r["Q_new"].to_numpy(dtype=float)):
                i = pos.get(d)
                if i is None:
                    continue
                if not np.isnan(q[i]):
                    raise SystemExit(f"gage {gid}: ledger date {d} is not a missing day of the base parquet (Q = {q[i]})")
                q2[i] = v; hit += 1
            if hit != len(r):
                raise SystemExit(f"gage {gid}: {len(r) - hit} ledger dates are not rows of the parquet")
            consumed[gid] = hit
            final = replay_gage(days, q2)
        else:
            final = base
        bmap = {y[0]: y[2] for y in base}
        fmap = {y[0]: y[2] for y in final}
        for y in sorted(fmap):
            yrows.append((gid, y, bmap[y], fmap[y]))
        bw = window_rows(gid, base, windows); fw = window_rows(gid, final, windows)
        n_unav = {w: 0 for w in [f"{w0}-{w1}" for w0, w1 in windows]}
        if gid in unav_by_gage:
            u = unav_by_gage[gid]
            uwy = water_year_of((pd.to_datetime(u["date"]).values.astype("datetime64[D]") - EPOCH).astype(np.int64))
            for (w0, w1) in windows:
                n_unav[f"{w0}-{w1}"] = int(((uwy >= w0) & (uwy <= w1)).sum())
        for b, f_ in zip(bw, fw):
            bset = set(map(int, b["valid_years"].split())) if b["valid_years"] else set()
            fset = set(map(int, f_["valid_years"].split())) if f_["valid_years"] else set()
            w = b["window"]
            rows.append({"gage_id": gid, "window": w, "in_product": gid in products.get(w, set()),
                         "baseline_n_valid": len(bset), "final_n_valid": len(fset),
                         "gained_years": " ".join(map(str, sorted(fset - bset))), "n_gained": len(fset - bset),
                         "lost_years": " ".join(map(str, sorted(bset - fset))), "n_lost": len(bset - fset),
                         "baseline_total_possible": b["total_possible"], "final_total_possible": f_["total_possible"],
                         "baseline_frac": b["frac"], "final_frac": f_["frac"],
                         "baseline_qualifies": b["qualifies"], "final_qualifies": f_["qualifies"],
                         "n_restored_days_in_window": int(((water_year_of((pd.to_datetime(by_gage[gid]["date"]).values.astype("datetime64[D]") - EPOCH).astype(np.int64)) >= int(w.split("-")[0])) &
                                                           (water_year_of((pd.to_datetime(by_gage[gid]["date"]).values.astype("datetime64[D]") - EPOCH).astype(np.int64)) <= int(w.split("-")[1]))).sum()) if gid in by_gage else 0,
                         "n_source_unavailable_in_window": n_unav[w],
                         "baseline_valid_years": b["valid_years"], "final_valid_years": f_["valid_years"]})
        n += 1
        if n % 1000 == 0:
            log(f"  {n} gages ({time.time() - t0:.0f}s)")
    unconsumed = {g: int((restored["gage_id"] == g).sum()) - n_ for g, n_ in consumed.items() if int((restored["gage_id"] == g).sum()) != n_}
    if unconsumed:
        raise SystemExit(f"{len(unconsumed)} ledger gages were not (fully) consumed — not in the parquet or rows missing: {dict(list(unconsumed.items())[:5])}")
    orc = pd.DataFrame(rows)
    orc.to_csv(os.path.join(a.out_dir, "acceptance_oracle.csv"), index=False)
    pd.DataFrame(yrows, columns=["gage_id", "water_year", "baseline_status", "final_status"]).to_csv(os.path.join(a.out_dir, "oracle_years.csv"), index=False)
    wq = orc[(~orc["in_product"]) & (~orc["baseline_qualifies"]) & orc["final_qualifies"]]
    wq.to_csv(os.path.join(a.out_dir, "would_qualify_outside_products.csv"), index=False)
    lost = orc[orc["n_lost"] > 0]
    summary = {
        "generated": dt.datetime.now().isoformat(timespec="seconds"),
        "policy": {"splice_codes": list(SPLICE_CODES), "max_raw_na": MAX_RAW_NA, "max_gap": MAX_GAP, "min_years": MIN_YEARS, "min_frac": MIN_FRAC},
        "inputs": {"parquet": {"path": os.path.abspath(a.parquet), "sha256": sha256_file(a.parquet)},
                   "ledger": {"path": os.path.abspath(a.ledger), "sha256": sha256_file(a.ledger)},
                   "masked": {"path": os.path.abspath(a.masked), "sha256": sha256_file(a.masked)},
                   "config": config, "tool_sha256": sha256_file(os.path.abspath(__file__))},
        "ledger_consumed": {"restored_rows": int(len(restored)), "gages": len(consumed)},
        "windows": [f"{w0}-{w1}" for w0, w1 in windows],
        "gages_replayed": n,
        "per_window": {}, "n_gage_windows_with_lost_years": int(len(lost)),
    }
    for w in summary["windows"]:
        o = orc[orc["window"] == w]; p = o[o["in_product"]]
        summary["per_window"][w] = {
            "product_gages": int(p["in_product"].sum()), "product_gages_gaining_years": int((p["n_gained"] > 0).sum()),
            "product_years_gained": int(p["n_gained"].sum()), "product_gages_with_restored_days": int((p["n_restored_days_in_window"] > 0).sum()),
            "product_gages_with_source_unavailable_targets": int((p["n_source_unavailable_in_window"] > 0).sum()),
            "product_gages_losing_qualification": int((p["baseline_qualifies"] & ~p["final_qualifies"]).sum()),
            "would_qualify_outside_product": sorted(wq.loc[wq["window"] == w, "gage_id"].tolist()),
            "all_gages_gaining_years": int((o["n_gained"] > 0).sum()), "all_years_gained": int(o["n_gained"].sum()),
        }
    json.dump(summary, open(os.path.join(a.out_dir, "oracle.json"), "w"), indent=1)
    log(f"oracle: {len(orc):,} gage-windows; lost-year rows {len(lost)}; per window " + json.dumps({w: {k: v for k, v in d.items() if k != 'would_qualify_outside_product'} for w, d in summary['per_window'].items()}))
    if len(lost):
        print(lost.head(10).to_string())
        log("ORACLE FAILED: a restoration can never remove a valid year")
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("replay"); p.add_argument("--parquet", required=True); p.add_argument("--out-dir", required=True)
    p.add_argument("--census-years"); p.add_argument("--census"); p.add_argument("--windows", default="1993-2025,1980-2025")
    p.set_defaults(func=cmd_replay)
    p = sub.add_parser("oracle"); p.add_argument("--parquet", required=True); p.add_argument("--ledger", required=True)
    p.add_argument("--masked", required=True); p.add_argument("--out-dir", required=True)
    p.add_argument("--product-csv", action="append", help="WINDOW=path of the product's signatures CSV (its gage set)")
    p.add_argument("--windows", default="1993-2025,1980-2025")
    p.add_argument("--config", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "config", "signatures_config.json"))
    p.set_defaults(func=cmd_oracle)
    a = ap.parse_args(); a.func(a)


if __name__ == "__main__":
    main()

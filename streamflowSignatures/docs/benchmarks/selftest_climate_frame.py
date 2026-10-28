#!/usr/bin/env python
"""selftest_climate_frame.py — the Python side of the climate-frame fixture self-test
(climate_frame.py; plan §10.2 tests (i)–(iv), 2026-10-26).

    python docs/benchmarks/selftest_climate_frame.py OUT_DIR [--rtol 1e-9] [--summary-rtol 1e-6]

Reads the fixture that selftest_climate_frame.jl wrote to OUT_DIR (bit-identical parquet
inputs), builds the frame in both modes with the Python helper, runs the Python
preprocessor and the production orchestrator call, repeats the fixture checks, and
compares the annual values and the summary values with canonical Julia per mode: annual
NaN patterns must agree and the values match within --rtol (the E7 annual tolerance; the
bit-identical count is reported), summaries within --summary-rtol with gate-state
agreement. Rank statistics (Mann-Kendall, Spearman, the Pettitt fields) of a base whose
annual series holds a tie or a near-tie (two values within 1e-9 relative) AND whose annual
values are not bit-identical in the two languages are reported, not gated: a last-bit
difference then reorders ranks — the named residual class of docs/CROSS_LANGUAGE_STATUS.md.
Exit 1 on any failed check.
"""
SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval",
            "_mean", "_median", "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
            "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
            "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"]
RANK_SUFFIXES = {s for s in SUFFIXES if s not in ("_senn_slp", "_linear_slp", "_mean", "_median")}


def split_key(k):
    for s in sorted(SUFFIXES, key=len, reverse=True):
        if k.endswith(s):
            return k[:-len(s)], s
    return k, None

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))
from climate_frame import climate_frame_mode, frame_with_climate  # noqa: E402
from streamflow_signatures import (AnnualCollector, calculate_all_signatures,  # noqa: E402
                                   preprocess_daily_data)
from streamflow_signatures.stats import generate_stats  # noqa: E402
from streamflow_signatures.config import (CP_END_WATER_YEAR, CP_MIN_SEGMENT_OBS,  # noqa: E402
                                          CP_MIN_TOTAL_OBS, CP_START_WATER_YEAR,
                                          MIN_VALUES_FOR_STATS, NA_DECADE_MIN_FRACTION,
                                          NA_TREND_MIN_FRACTION)

GID = "CF1"
n_pass = 0
n_fail = 0


def check(label, cond):
    global n_pass, n_fail
    if cond:
        n_pass += 1
        print(f"  PASS: {label}")
    else:
        n_fail += 1
        print(f"  FAIL: {label}")


# default mode (user decision 2026-10-28): "outer" when the variable is unset
_cf_saved = os.environ.pop("STREAMFLOW_CLIMATE_FRAME", None)
try:
    check('climate_frame_mode() is "outer" when STREAMFLOW_CLIMATE_FRAME is unset (default since 2026-10-28)',
          climate_frame_mode() == "outer")
finally:
    if _cf_saved is not None:
        os.environ["STREAMFLOW_CLIMATE_FRAME"] = _cf_saved


def same(a, b):
    a = float(a); b = float(b)
    return (math.isnan(a) and math.isnan(b)) or a == b


def close(a, b, atol=1e-9, rtol=1e-9):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def run_full(pp):
    """The production orchestrator call of run_python_benchmark.py (new path)."""
    data = pp["data"]
    cdata = data[data["water_year"].isin(set(pp["valid_climate_years"]))]
    sdata = data[data["water_year"].isin(set(pp.get("valid_swe_years", [])))]
    tdata = data[data["water_year"].isin(set(pp.get("valid_temp_years", [])))]
    cp = {"start_year": CP_START_WATER_YEAR, "end_year": CP_END_WATER_YEAR,
          "min_total_obs": CP_MIN_TOTAL_OBS, "min_segment_obs": CP_MIN_SEGMENT_OBS}
    coll = AnnualCollector()
    res = calculate_all_signatures(
        data, True, seasonal_flags=pp.get("seasonal_flags"),
        trend_completeness=NA_TREND_MIN_FRACTION, decade_completeness=NA_DECADE_MIN_FRACTION,
        climate_data=cdata, area_normalized=True, changepoint=cp,
        min_values_for_stats=MIN_VALUES_FOR_STATS, collector=coll, snow_data=sdata,
        snow_climate_years=pp["valid_climate_years"], gage_id=GID, temp_data=tdata)
    ann = pd.DataFrame({"signature": list(coll.signature), "water_year": [int(y) for y in coll.water_year],
                        "value": [float(v) for v in coll.value]})
    return ann, res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir")
    ap.add_argument("--rtol", type=float, default=1e-9, help="annual-value tolerance (relative)")
    ap.add_argument("--atol", type=float, default=1e-12)
    ap.add_argument("--summary-rtol", type=float, default=1e-6)
    a = ap.parse_args()
    out = Path(a.out_dir)

    sf = pd.read_parquet(out / "climate_frame_fixture_streamflow.parquet")
    cl = pd.read_parquet(out / "climate_frame_fixture_climate.parquet")
    for df in (sf, cl):
        df["gage_id"] = df["gage_id"].astype(str)
        df["date"] = pd.to_datetime(df["date"])
    ref = {m: json.load(open(out / f"climate_frame_fixture_summary_{m}.json")) for m in ("left", "outer")}
    ann_ref = {m: pd.read_parquet(out / f"climate_frame_fixture_annual_{m}.parquet") for m in ("left", "outer")}

    frames = {m: frame_with_climate(sf, cl, m) for m in ("left", "outer")}
    pp = {m: preprocess_daily_data(frames[m]) for m in ("left", "outer")}
    L, O = frames["left"], frames["outer"]

    print("\n(ii) frame composition")
    check("left frame = the streamflow rows", len(L) == len(sf) and set(L["date"]) == set(sf["date"]))
    added = sorted(set(O["date"]) - set(L["date"]))
    check(f"outer adds exactly the {ref['outer']['n_added_rows']} absent dates of the Julia reference, date for date",
          [d.strftime("%Y-%m-%d") for d in added] == ref["outer"]["added_dates"] and len(O) == ref["outer"]["n_frame_rows"])
    check("left adds nothing (the Julia reference lists no added date)", ref["left"]["added_dates"] == [])
    check("no added row in the absent water year 2005", 2005 not in set(O["water_year"]))
    check("no added row before the record start or after its end",
          O["date"].min() == sf["date"].min() and O["date"].max() == sf["date"].max())
    check("the water-year set is identical between modes", set(L["water_year"]) == set(O["water_year"]))
    check("Q is NaN on every added row", O["Q"].iloc[len(L):].isna().all())
    check("the joined rows of the outer frame equal the left frame", O.iloc[:len(L)][L.columns].equals(L))

    print("\n(ii) the runner's window denominator (run_python_benchmark.py: total_possible from the frame's water years)")
    def window_denominator(frame, p, start_wy, end_wy):
        all_wy = frame["water_year"].dropna().astype(int).unique()
        wy_in_range = [yr for yr in all_wy if start_wy <= yr <= end_wy]
        if not wy_in_range:
            return {"total_possible": 0, "n_valid": 0}
        return {"total_possible": int(max(wy_in_range) - start_wy + 1),
                "n_valid": int(sum(1 for yr in p["valid_years"] if start_wy <= yr <= end_wy))}
    for m in ("left", "outer"):
        got = {k: window_denominator(frames[m], pp[m], *map(int, k.split("-"))) for k in ref[m]["windows"]}
        check(f"{m}: total_possible and the valid-year count equal Julia for every window", got == ref[m]["windows"])

    print("\n(ii) Q-only preprocessor outputs identical between modes")
    pl, po = pp["left"], pp["outer"]
    check("valid_years identical", list(pl["valid_years"]) == list(po["valid_years"]))
    check("rejected_years identical, WY 2005 absent", pl["rejected_years"] == po["rejected_years"]
          and 2005 not in {r["water_year"] for r in pl["rejected_years"]})
    check("seasonal_flags identical", pl["seasonal_flags"] == po["seasonal_flags"])
    check("diagnostics identical", pl["diagnostics"] == po["diagnostics"])
    dl, do_ = pl["data"].reset_index(drop=True), po["data"].reset_index(drop=True)
    check("preprocessed Q identical (bitwise, NaN == NaN)", len(dl) == len(do_)
          and (dl["date"].to_numpy() == do_["date"].to_numpy()).all()
          and all(same(x, y) for x, y in zip(dl["Q"].to_numpy(), do_["Q"].to_numpy())))
    cp = {"start_year": CP_START_WATER_YEAR, "end_year": CP_END_WATER_YEAR,
          "min_total_obs": CP_MIN_TOTAL_OBS, "min_segment_obs": CP_MIN_SEGMENT_OBS}
    qonly = {m: calculate_all_signatures(pp[m]["data"], False, seasonal_flags=pp[m]["seasonal_flags"],
                                         trend_completeness=NA_TREND_MIN_FRACTION,
                                         decade_completeness=NA_DECADE_MIN_FRACTION, changepoint=cp,
                                         min_values_for_stats=MIN_VALUES_FOR_STATS, gage_id=GID)
             for m in ("left", "outer")}
    qkeys = sorted(qonly["left"])
    check(f"Q-only orchestrator: the same {len(qkeys)} keys in both modes (Julia: {ref['left']['checks']['n_qonly_keys']})",
          qkeys == sorted(qonly["outer"]) and len(qkeys) == ref["left"]["checks"]["n_qonly_keys"])
    check("Q-only orchestrator: every value bit-identical between modes",
          all(same(qonly["left"][k], qonly["outer"][k]) for k in qkeys))

    print("\n(i) + (iii) against the Julia reference numbers")
    anns = {}; sums = {}
    for m in ("left", "outer"):
        anns[m], sums[m] = run_full(pp[m])

    def ppt_at(m, d):
        df = pp[m]["data"]
        hit = df.loc[pd.to_datetime(df["date"]) == pd.Timestamp(d), "PPT"]
        return float(hit.iloc[0]) if len(hit) else float("nan")

    def prcp(m, wy):
        x = anns[m]
        hit = x.loc[(x["signature"] == "prcp_annual") & (x["water_year"] == wy), "value"]
        return float(hit.iloc[0]) if len(hit) else float("nan")

    for m in ("left", "outer"):
        c = ref[m]["checks"]
        check(f"{m}: preprocessed PPT on 2000-12-30 / 12-31 equals Julia",
              same(ppt_at(m, "2000-12-30"), c["ppt_2000_12_30"]) and same(ppt_at(m, "2000-12-31"), c["ppt_2000_12_31"]))
        check(f"{m}: prcp_annual WY 1995 equals Julia ({c['prcp_annual_1995']:.6f})", same(prcp(m, 1995), c["prcp_annual_1995"]))
        for key in ("valid_years", "valid_climate_years", "valid_swe_years", "valid_temp_years"):
            check(f"{m}: {key} equal Julia", [int(y) for y in pp[m].get(key, [])] == [int(y) for y in ref[m][key]])
    c = ref["outer"]["checks"]
    check("outer: prcp_annual WY 1995 = Daymet sum; left = Daymet sum − interpolation residual",
          close(prcp("outer", 1995), c["daymet_sum_1995"]) and close(prcp("left", 1995), c["daymet_sum_1995"] - c["interp_resid_1995"]))
    check("WY 2009 and 2013 climate-valid only in outer mode",
          sorted(set(pp["outer"]["valid_climate_years"]) - set(pp["left"]["valid_climate_years"])) == [2009, 2013])

    print("\n(iv) cross-language: annual values and summaries vs Julia, per mode")
    for m in ("left", "outer"):
        r = ann_ref[m].copy(); r["signature"] = r["signature"].astype(str); r["water_year"] = r["water_year"].astype(int)
        mm = r.merge(anns[m], on=["signature", "water_year"], how="outer", suffixes=("_jl", "_py"), indicator=True)
        unmatched = int((mm["_merge"] != "both").sum())
        both = mm[mm["_merge"] == "both"]
        nan_mis = int((both["value_jl"].isna() != both["value_py"].isna()).sum())
        ok = both[both["value_jl"].notna() & both["value_py"].notna()]
        vj, vp = ok["value_jl"].to_numpy(), ok["value_py"].to_numpy()
        n_bit = int((vj == vp).sum())
        rel = np.abs(vj - vp) / np.maximum(np.maximum(np.abs(vj), np.abs(vp)), 1e-300)
        within = np.abs(vj - vp) <= a.atol + a.rtol * np.maximum(np.abs(vj), np.abs(vp))
        worst_i = int(np.argmax(rel)) if len(rel) else None
        print(f"  {m}: {len(r)} Julia rows, {len(anns[m])} python rows, unmatched {unmatched}, NaN-pattern mismatches {nan_mis}, "
              f"bit-identical {n_bit}/{len(ok)}, worst rel {rel[worst_i]:.3e} at {ok['signature'].iloc[worst_i]} WY {ok['water_year'].iloc[worst_i]}"
              if worst_i is not None else f"  {m}: no comparable annual values")
        check(f"{m}: annual values match Julia (keys, NaN pattern, rel <= {a.rtol:g})",
              unmatched == 0 and nan_mis == 0 and bool(within.all()))
        # tie-sensitive bases: a tie or near-tie in Julia's annual series and annual values
        # that are not bit-identical in the port
        tie_sensitive = set()
        for sig, g in both.groupby("signature"):
            vj_s = np.sort(g["value_jl"].to_numpy())
            if len(vj_s) < 2:
                continue
            gaps = np.diff(vj_s) / np.maximum(np.abs(vj_s[1:]), 1e-300)
            a_, b_ = g["value_jl"].to_numpy(), g["value_py"].to_numpy()
            bit_identical = bool(((a_ == b_) | (np.isnan(a_) & np.isnan(b_))).all())
            if (gaps <= 1e-9).any() and not bit_identical:
                tie_sensitive.add(sig)
        jl = ref[m]["values"]; py = sums[m]
        missing = sorted(set(jl) - set(py)); extra = sorted(set(py) - set(jl))
        bad = []; diag = []; worst = 0.0; worst_key = None
        for k, v in jl.items():
            if k not in py:
                continue
            vp_ = py[k]
            base, suf = split_key(k)
            sensitive = suf in RANK_SUFFIXES and base in tie_sensitive
            jn = v == "NaN"; pn = isinstance(vp_, float) and math.isnan(vp_)
            if jn or pn:
                if jn != pn:
                    bad.append((k, v, vp_))     # a gate-state difference is never exempt (review 2026-10-26)
                continue
            v = float(v); vp_ = float(vp_)
            scale = max(abs(v), abs(vp_), 1e-300); rl = abs(v - vp_) / scale
            if abs(v - vp_) > a.atol + a.summary_rtol * scale:
                (diag if sensitive else bad).append((k, v, vp_))
            elif rl > worst:
                worst, worst_key = rl, k
        print(f"  {m}: {len(jl)} Julia summary values, python {len(py)}, missing {len(missing)}, extra {len(extra)}, "
              f"worst gated rel {worst:.3e} at {worst_key}, mismatches {len(bad)}; "
              f"tie-sensitive rank statistics reported, not gated: {len(diag)} on {sorted(tie_sensitive)}")
        for b in bad[:20]:
            print("    MISMATCH", b)
        for b in diag[:10]:
            print("    tie-sensitive", b)
        check(f"{m}: summary values match Julia (gate state always; rel <= {a.summary_rtol:g}; numeric rank statistics on tie-bearing series excepted)",
              not missing and not extra and not bad)
        # the exemption is justified only if the port's statistics code reproduces Julia's
        # rank statistics when fed JULIA's annual values (review 2026-10-26, finding 6)
        n_cmp = 0; n_skip = 0; bad2 = []
        for base in sorted(tie_sensitive):
            g = ann_ref[m][ann_ref[m]["signature"].astype(str) == base]
            df = pd.DataFrame({"water_year": g["water_year"].astype(int).to_numpy(), base: g["value"].to_numpy()})
            st = generate_stats(df, value_cols=[base], trend_completeness=NA_TREND_MIN_FRACTION,
                                decade_completeness=NA_DECADE_MIN_FRACTION, changepoint=cp,
                                min_values_for_stats=MIN_VALUES_FOR_STATS)
            for suf in sorted(RANK_SUFFIXES):
                k = base + suf
                if k not in jl or k not in st:
                    continue
                vj = jl[k]; vp_ = st[k]
                jn = vj == "NaN"; pn = isinstance(vp_, float) and math.isnan(vp_)
                if jn or pn:
                    n_skip += 1      # a family-specific gate (floor / completeness exemption) the direct call cannot reproduce
                    continue
                n_cmp += 1
                vj = float(vj); vp_ = float(vp_)
                if abs(vj - vp_) > 1e-12 + 1e-9 * max(abs(vj), abs(vp_)):
                    bad2.append((k, vj, vp_))
        for b in bad2[:10]:
            print("    STATS MISMATCH on Julia's values", b)
        check(f"{m}: the port's rank statistics reproduce Julia's on Julia's own annual values for the tie-sensitive bases "
              f"({n_cmp} compared, {n_skip} gated differently in the direct call)", not bad2)

    print(f"\n{n_pass} passed, {n_fail} failed")
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()

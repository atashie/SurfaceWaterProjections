#!/usr/bin/env python
"""Compare the Python climate metrics with canonical Julia on the shared fixture.

    python docs/benchmarks/crosscheck_climate_fixture.py OUT_DIR [--rtol 4e-14]

Reads OUT_DIR/climate_fixture_input.parquet (written by crosscheck_climate_fixture.jl) so
the input is bit-identical, runs the Python preprocessor + both climate families with
the production gates, requires the ANNUAL values to be bit-identical with
OUT_DIR/climate_fixture_julia_annual.parquet, and compares every one of the 39 x 16 summary
values with OUT_DIR/climate_fixture_julia.json. NaN must match NaN. Exit 1 on any mismatch.
"""
import argparse
import json
import math
import sys

import pandas as pd

from streamflow_signatures.io import add_water_year_columns, preprocess_daily_data
from streamflow_signatures.climate import calculate_precipitation_metrics, calculate_temperature_metrics
from streamflow_signatures.config import (NA_TREND_MIN_FRACTION, NA_DECADE_MIN_FRACTION,
                                          MIN_VALUES_FOR_STATS, CP_START_WATER_YEAR,
                                          CP_END_WATER_YEAR, CP_MIN_TOTAL_OBS, CP_MIN_SEGMENT_OBS)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir")
    ap.add_argument("--rtol", type=float, default=4e-14)
    ap.add_argument("--atol", type=float, default=1e-12,
                    help="absolute tolerance for values whose true result is ~0 (e.g. a Pettitt "
                         "delta of two equal means): |d| <= atol + rtol * scale")
    ap.add_argument("--mk-rtol", type=float, default=None,
                    help="looser relative tolerance for the Mann-Kendall-derived statistics "
                         "(_mk_rho, _mk_pval, _pettitt_pre/post_mk_pval); default = --rtol")
    a = ap.parse_args()
    mk_rtol = a.rtol if a.mk_rtol is None else a.mk_rtol
    raw = pd.read_parquet(f"{a.out_dir}/climate_fixture_input.parquet")  # exact doubles (CSV parsers are not)
    raw["gage_id"] = raw["gage_id"].astype(str)
    raw["date"] = pd.to_datetime(raw["date"])
    ref = json.load(open(f"{a.out_dir}/climate_fixture_julia.json"))
    df = add_water_year_columns(raw)
    pp = preprocess_daily_data(df)
    assert pp["valid_temp_years"] == ref["valid_temp_years"], (pp["valid_temp_years"], ref["valid_temp_years"])
    assert pp["valid_climate_years"] == ref["valid_climate_years"]
    data = pp["data"]
    cdata = data[data["water_year"].isin(set(pp["valid_climate_years"]))]
    tdata = data[data["water_year"].isin(set(pp["valid_temp_years"]))]
    cp = {"start_year": CP_START_WATER_YEAR, "end_year": CP_END_WATER_YEAR,
          "min_total_obs": CP_MIN_TOTAL_OBS, "min_segment_obs": CP_MIN_SEGMENT_OBS}
    kw = dict(trend_completeness=NA_TREND_MIN_FRACTION, decade_completeness=NA_DECADE_MIN_FRACTION,
              min_values_for_stats=MIN_VALUES_FOR_STATS, changepoint=cp)
    from streamflow_signatures import AnnualCollector
    coll = AnnualCollector()
    res = {}
    res.update(calculate_precipitation_metrics(cdata, collector=coll, **kw))
    res.update(calculate_temperature_metrics(tdata, collector=coll, **kw))
    # Annual values must be BIT-IDENTICAL (same inputs, same sequential arithmetic)
    ann_ref = pd.read_parquet(f"{a.out_dir}/climate_fixture_julia_annual.parquet")
    ann_ref["signature"] = ann_ref["signature"].astype(str); ann_ref["water_year"] = ann_ref["water_year"].astype(int)
    ann_py = pd.DataFrame({"signature": list(coll.signature), "water_year": [int(y) for y in coll.water_year],
                           "value": [float(v) for v in coll.value]})
    m = ann_ref.merge(ann_py, on=["signature", "water_year"], how="outer", suffixes=("_jl", "_py"), indicator=True)
    ann_missing = int((m["_merge"] != "both").sum())
    both = m[m["_merge"] == "both"]
    nan_mismatch = int((both["value_jl"].isna() != both["value_py"].isna()).sum())
    ok_both = both[both["value_jl"].notna() & both["value_py"].notna()]
    ann_diff = int((ok_both["value_jl"].to_numpy() != ok_both["value_py"].to_numpy()).sum())
    print(f"annual values: {len(ann_ref)} Julia rows, {len(ann_py)} python rows, unmatched {ann_missing}, "
          f"NaN-pattern mismatches {nan_mismatch}, non-identical values {ann_diff}")
    jl = ref["values"]
    missing = sorted(set(jl) - set(res)); extra = sorted(set(res) - set(jl))
    worst = 0.0; worst_key = None; bad = []
    for k, vj in jl.items():
        if k not in res:
            continue
        vp = res[k]
        vj_nan = vj == "NaN"
        vp_nan = isinstance(vp, float) and math.isnan(vp)
        if vj_nan or vp_nan:
            if vj_nan != vp_nan:
                bad.append((k, vj, vp, "nan-mismatch"))
            continue
        vj = float(vj); vp = float(vp)
        scale = max(abs(vj), abs(vp), 1e-300)
        rel = abs(vj - vp) / scale
        if rel > worst:
            worst, worst_key = rel, k
        tol = mk_rtol if k.endswith(("_mk_rho", "_mk_pval", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval")) else a.rtol
        if abs(vj - vp) > a.atol + tol * scale:
            bad.append((k, vj, vp, rel))
    print(f"compared {len(jl)} Julia values; python has {len(res)}; missing in python: {len(missing)}; extra: {len(extra)}")
    print(f"worst relative difference: {worst:.3e} at {worst_key}")
    for b in bad[:30]:
        print("  MISMATCH", b)
    ok = not missing and not extra and not bad and ann_missing == 0 and nan_mismatch == 0 and ann_diff == 0
    print("CROSSCHECK", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

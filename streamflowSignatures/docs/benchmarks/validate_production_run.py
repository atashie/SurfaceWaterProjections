"""Artifact-level validation gates for a production benchmark run.

Repo-resident successor of the (drive-lost) per-run validate_gates scripts used for
the July 2026 reruns. Gates: window bounds, optional column-set / gage-overlap
comparison against a reference run, snow sanity (incl. snow-implies-climate;
Daymet covers Canadian gages), annual-parquet volume/keys vs the timing JSON,
per-gage year-count consistency, and (when a stats floor is active) a floor-respect
check on non-exempt signatures.

Complements: validate_annual_values.py (parquet<->summary value consistency) and
audit_qualification.jl (independent window/fraction inclusion audit).

Usage:
  python validate_production_run.py --csv <signatures.csv> --annual <annual.parquet>
      --timing <timing.json> --start-year 1993 [--end-year 2025]
      [--reference-csv <other_run.csv>] [--floor 20]
"""
import argparse
import json
import sys

import numpy as np
import pandas as pd

FLOOR_EXEMPT = {"log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events",
                "concavity", "n_recession_events", "alpha_linear",
                "elasticity_rolling", "elasticity_annual"}

results = []
def gate(name, ok, detail):
    results.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--annual", required=True)
    ap.add_argument("--timing", required=True)
    ap.add_argument("--start-year", type=int, required=True)
    ap.add_argument("--end-year", type=int, default=2025,
                    help="max permissible water year (default 2025 = data end)")
    ap.add_argument("--reference-csv", default=None,
                    help="another run's CSV for column-set + overlap comparison")
    ap.add_argument("--floor", type=int, default=None,
                    help="stats floor in effect (checks non-exempt signatures)")
    args = ap.parse_args()

    prod = pd.read_csv(args.csv, dtype={"gage_id": str}, low_memory=False)
    print(f"run: {prod.shape}\n")

    # G1: window bounds
    g1 = (prod.start_water_year.min() >= args.start_year and
          prod.end_water_year.max() <= args.end_year)
    gate("1 window bounds", g1,
         f"start min={prod.start_water_year.min():.0f} (>= {args.start_year}), "
         f"end max={prod.end_water_year.max():.0f} (<= {args.end_year})")

    # G2/G5: reference comparison
    if args.reference_csv:
        ref = pd.read_csv(args.reference_csv, dtype={"gage_id": str}, low_memory=False)
        gate("2 column set == reference", list(prod.columns) == list(ref.columns),
             f"run cols={prod.shape[1]}, reference cols={ref.shape[1]}")
        sp, sr = set(prod.gage_id), set(ref.gage_id)
        m = prod.merge(ref[["gage_id", "latitude", "longitude"]], on="gage_id",
                       suffixes=("", "_r"))
        meta_ok = bool(np.isclose(m.latitude, m.latitude_r, rtol=1e-9).all() and
                       np.isclose(m.longitude, m.longitude_r, rtol=1e-9).all())
        gate("5 overlap + metadata consistency", meta_ok,
             f"run={len(sp)}, ref={len(sr)}, shared={len(sp & sr)}, "
             f"run-only={len(sp - sr)}, ref-only={len(sr - sp)}; lat/lon match={meta_ok}")

    # G3: snow sanity
    sm = pd.to_numeric(prod["swe_max_mean"], errors="coerce")
    ssm = pd.to_numeric(prod["ssm_mean"], errors="coerce")
    dowy = pd.to_numeric(prod["swe_max_dowy_mean"], errors="coerce")
    canadian = prod.gage_id.str.contains(r"[A-Za-z]", regex=True, na=False)
    an_false = prod["area_normalized"].astype(str).str.upper().eq("FALSE")
    sic = bool((sm.notna() <= (prod["annual_runoff_ratio_mean"].notna() | an_false)).all())
    g3 = ((sm.dropna() >= 0).all() and ssm.dropna().between(-1, 1).all()
          and dowy.dropna().between(1, 366).all() and sic)
    gate("3 snow sanity", g3,
         f"snow gages={sm.notna().sum()} (US={(sm.notna() & ~canadian).sum()}, "
         f"CAN={(sm.notna() & canadian).sum()}); ssm=[{ssm.min():.3f},{ssm.max():.3f}]; "
         f"snow_implies_climate={sic}")

    # G4: annual parquet volume/keys
    ann = pd.read_parquet(args.annual)
    with open(args.timing) as f:
        tj = json.load(f)
    # 100 = 90 (July 2026: 76 base + 14 snow) + 10 drought (2 metrics x 5 levels)
    g4 = (len(ann) == tj.get("n_annual_rows") and ann.signature.nunique() == 100
          and ann.water_year.min() >= args.start_year
          and ann.water_year.max() <= args.end_year
          and set(ann.gage_id.unique()) <= set(prod.gage_id)
          and not ann.duplicated(["gage_id", "signature", "water_year"]).any())
    gate("4 annual parquet volume/keys", g4,
         f"rows={len(ann)} (timing={tj.get('n_annual_rows')}), "
         f"sigs={ann.signature.nunique()}, wy=[{ann.water_year.min()},"
         f"{ann.water_year.max()}], dup="
         f"{ann.duplicated(['gage_id','signature','water_year']).sum()}")

    # G6: per-gage year counts
    span = prod.end_water_year - prod.start_water_year + 1
    g6 = bool((prod.num_water_years >= 20).all() and (prod.num_water_years <= span).all())
    gate("6 per-gage year-count consistency", g6,
         f"num_water_years=[{prod.num_water_years.min():.0f},"
         f"{prod.num_water_years.max():.0f}]")

    # G7: stats floor respected (non-exempt signatures)
    if args.floor:
        counts = (ann.dropna(subset=["value"])
                     .groupby(["gage_id", "signature"]).size().unstack(fill_value=0))
        bases = sorted({c[:-5] for c in prod.columns
                        if c.endswith("_mean") and "_pettitt_" not in c})
        violations = []
        for b in bases:
            if b in FLOOR_EXEMPT or b not in counts.columns:
                continue
            n = counts[b].reindex(prod["gage_id"]).fillna(0).to_numpy()
            has_mean = pd.to_numeric(prod[b + "_mean"], errors="coerce").notna().to_numpy()
            bad = has_mean & (n < args.floor)
            if bad.any():
                violations.append((b, int(bad.sum())))
        gate("7 stats floor respected", not violations,
             f"floor={args.floor}; non-exempt signatures with sub-floor non-NaN means: "
             f"{violations[:8] if violations else 'none'}")

    print("\n" + "=" * 60)
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print(f"GATES: {sum(1 for _, ok, _ in results if ok)}/{len(results)} PASS"
          + (f" — {n_fail} FAIL" if n_fail else " — ALL PASS"))
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()

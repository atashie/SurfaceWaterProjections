"""Cross-file consistency validation for the annual values export.

Verifies that the long-format annual parquet written by run_julia_benchmark.jl
contains exactly the series the summary statistics were computed from:
for every (gage, signature), the mean/median recomputed from the annual values
must match the `_mean` / `_median` columns in the summary CSV.

Consistency rules (see docs/plans/annual_values_export_plan.md section 4.2):
  1. Summary mean non-NaN  -> the parquet must hold >= 3 non-NaN rows for that
     (gage, signature) — generate_stats only computes a mean from >= min_rows
     (3) values, all of which were collected — AND the recomputed mean/median
     must match (rtol 1e-9). Fewer than 3 non-NaN rows behind a non-NaN mean
     is a coverage failure even if the values happen to average right.
  2. Summary mean NaN + no rows in parquet -> consistent (early empty_stats
     return, caller-side pruning to an empty frame, or a signature that never
     ran for this gage).
  3. Summary mean NaN + rows present with count(non-NaN) < 3 -> consistent
     (generate_stats min_rows guard NaNs the stats but values were computed).
  4. Summary mean NaN + rows present with count(non-NaN) >= 3 -> INCONSISTENT
     (generate_stats always computes mean/median when >= min_rows values
     survive NaN filtering — this flags a broken collector thread).

Coverage check: for every gage, the set of signatures with a non-NaN summary
mean must equal the set of signatures with >= 3 non-NaN rows in the parquet.

Zero-gage runs: when the benchmark qualifies zero gages it writes an empty
summary CSV and intentionally skips the annual parquet — that case exits 0.

Usage:
    python docs/benchmarks/validate_annual_values.py [--prefix julia]
    python docs/benchmarks/validate_annual_values.py --annual path.parquet --summary path.csv
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

RTOL = 1e-9
ATOL = 1e-12
MIN_ROWS = 3  # generate_stats min_rows


def load_summary_bases(summary_df):
    """Signature base names = columns ending '_mean', excluding Pettitt fields."""
    bases = []
    for col in summary_df.columns:
        if col.endswith("_mean") and "_pettitt_" not in col:
            bases.append(col[: -len("_mean")])
    return sorted(bases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default="julia",
                        help="Output prefix used by the benchmark (default: julia)")
    parser.add_argument("--annual", default=None,
                        help="Path to the annual values parquet (overrides --prefix)")
    parser.add_argument("--summary", default=None,
                        help="Path to the summary signatures CSV (overrides --prefix)")
    parser.add_argument("--max-report", type=int, default=25,
                        help="Max mismatches to print per category")
    parser.add_argument("--floor", type=int, default=None,
                        help="stats floor (min_values_for_stats): non-exempt "
                             "signatures with fewer non-NaN annual values are "
                             "EXPECTED to have NaN statistics")
    parser.add_argument("--floor-exempt", default=("log_a_pointcloud,log_a_events,"
                        "b_pointcloud,b_events,concavity,n_recession_events,"
                        "alpha_linear,elasticity_rolling,elasticity_annual"),
                        help="comma-separated signature bases exempt from --floor "
                             "(default: recession + elasticity families)")
    args = parser.parse_args()
    floor_exempt = set(s for s in args.floor_exempt.split(",") if s)

    annual_path = args.annual or os.path.join(SCRIPT_DIR, f"{args.prefix}_signatures_annual.parquet")
    summary_path = args.summary or os.path.join(SCRIPT_DIR, f"{args.prefix}_signatures.csv")

    if not os.path.exists(summary_path):
        print(f"ERROR: summary CSV not found: {summary_path}")
        return 1

    print(f"Summary: {summary_path}")
    print(f"Annual:  {annual_path}")

    summary = pd.read_csv(summary_path, dtype={"gage_id": str})

    if not os.path.exists(annual_path):
        if summary.shape[0] == 0:
            # Zero-gage benchmark run: empty summary CSV, annual parquet
            # intentionally skipped — nothing to validate.
            print("PASS (no-op): zero-gage run — empty summary CSV and no "
                  "annual parquet, as expected.")
            return 0
        print(f"ERROR: annual parquet not found: {annual_path}")
        print("(The summary CSV has rows, so the benchmark should have written "
              "the annual parquet unless annual_values.save=false.)")
        return 1
    annual = pd.read_parquet(annual_path)
    annual["gage_id"] = annual["gage_id"].astype(str)
    annual["signature"] = annual["signature"].astype(str)

    bases = load_summary_bases(summary)
    print(f"Gages in summary: {summary.shape[0]}, signature bases: {len(bases)}")
    print(f"Annual rows: {len(annual)}, gages: {annual['gage_id'].nunique()}, "
          f"signatures: {annual['signature'].nunique()}")

    # Duplicate key check
    dups = annual.duplicated(subset=["gage_id", "signature", "water_year"]).sum()
    if dups:
        print(f"FAIL: {dups} duplicate (gage_id, signature, water_year) rows")

    # Signatures in the parquet that the summary has never heard of
    unknown_sigs = sorted(set(annual["signature"].unique()) - set(bases))
    if unknown_sigs:
        print(f"FAIL: {len(unknown_sigs)} parquet signatures missing from summary: "
              f"{unknown_sigs[:10]}")

    # Recompute per-(gage, signature) stats from the annual values
    grouped = annual.groupby(["gage_id", "signature"])["value"].agg(
        recomputed_mean=lambda v: np.nanmean(v) if np.isfinite(v).any() else np.nan,
        recomputed_median=lambda v: np.nanmedian(v) if np.isfinite(v).any() else np.nan,
        n_nonnan=lambda v: int(np.isfinite(v).sum()),
    ).reset_index()

    # Long-format view of the summary means/medians
    mean_cols = {b: f"{b}_mean" for b in bases}
    median_cols = {b: f"{b}_median" for b in bases}
    records = []
    for b in bases:
        sub = pd.DataFrame({
            "gage_id": summary["gage_id"],
            "signature": b,
            "summary_mean": summary[mean_cols[b]],
            "summary_median": summary[median_cols[b]] if median_cols[b] in summary.columns else np.nan,
        })
        records.append(sub)
    summary_long = pd.concat(records, ignore_index=True)

    # Named indicator: itertuples() renames a leading-underscore "_merge" column to a
    # positional field, so row._merge would raise AttributeError.
    merged = summary_long.merge(grouped, on=["gage_id", "signature"], how="outer",
                                indicator="merge_side")

    n_checked = 0
    mismatch_value = []     # rule 1 violations
    inconsistent_nan = []   # rule 4 violations
    orphan_rows = []        # parquet rows for a (gage, sig) absent from summary
    coverage_miss = []      # summary non-NaN mean but no parquet rows

    for row in merged.itertuples(index=False):
        in_summary = row.merge_side in ("both", "left_only")
        in_annual = row.merge_side in ("both", "right_only")
        s_mean = row.summary_mean if in_summary else np.nan
        has_rows = in_annual
        n_nonnan = int(row.n_nonnan) if has_rows and pd.notna(row.n_nonnan) else 0

        if not in_summary:
            orphan_rows.append((row.gage_id, row.signature))
            continue

        # Effective minimum for THIS signature: the stats floor raises the bar for
        # non-exempt signatures (recession/elasticity keep the plain MIN_ROWS)
        eff_min = MIN_ROWS
        if args.floor is not None and row.signature not in floor_exempt:
            eff_min = max(MIN_ROWS, args.floor)

        if pd.notna(s_mean):
            n_checked += 1
            # generate_stats only computes a mean from >= eff_min non-NaN
            # values, all of which the collector captured — so anything less
            # here is missing rows (or an unmasked/unfloored statistic).
            if not has_rows or n_nonnan < eff_min:
                coverage_miss.append((row.gage_id, row.signature))
            else:
                ok_mean = np.isclose(row.recomputed_mean, s_mean, rtol=RTOL, atol=ATOL)
                ok_median = (pd.isna(row.summary_median) and not np.isfinite(row.recomputed_median)) or \
                    np.isclose(row.recomputed_median, row.summary_median, rtol=RTOL, atol=ATOL)
                if not (ok_mean and ok_median):
                    mismatch_value.append((row.gage_id, row.signature, s_mean,
                                           row.recomputed_mean, row.summary_median,
                                           row.recomputed_median))
        else:
            # Summary mean is NaN
            if has_rows and n_nonnan >= eff_min:
                inconsistent_nan.append((row.gage_id, row.signature, n_nonnan))
            # else: rules 2 & 3 — consistent (incl. floor-gated series)

    print()
    print("=" * 70)
    print("VALIDATION RESULTS")
    print("=" * 70)
    floor_note = (f" [stats floor {args.floor}, {len(floor_exempt)} exempt bases]"
                  if args.floor is not None else "")
    print(f"(gage, signature) pairs with non-NaN summary mean: {n_checked}{floor_note}")
    print(f"Value mismatches (rule 1): {len(mismatch_value)}")
    print(f"NaN-summary-but->=min-values (rule 4): {len(inconsistent_nan)}")
    print(f"Coverage misses (non-NaN mean, < min non-NaN annual rows): {len(coverage_miss)}")
    print(f"Orphan parquet groups (unknown to summary): {len(orphan_rows)}")
    print(f"Duplicate keys: {dups}")

    for label, items in [("VALUE MISMATCH", mismatch_value),
                         ("RULE-4 INCONSISTENCY", inconsistent_nan),
                         ("COVERAGE MISS", coverage_miss),
                         ("ORPHAN", orphan_rows)]:
        for item in items[: args.max_report]:
            print(f"  {label}: {item}")
        if len(items) > args.max_report:
            print(f"  ... and {len(items) - args.max_report} more {label} entries")

    failed = bool(mismatch_value or inconsistent_nan or coverage_miss or
                  orphan_rows or dups or unknown_sigs)
    print()
    print("FAIL" if failed else "PASS: annual values are consistent with the summary CSV")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

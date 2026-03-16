"""
Smoke Test - Runs signature extraction on a small subset of gages.

Usage:
    python tests/smoke_test.py

Requires production parquet files. Edit paths below if needed.
"""

import sys
import time
from pathlib import Path

# Add parent for imports when running as script
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd

from streamflow_signatures import (
    read_parquet,
    add_water_year_columns,
    filter_qualifying_years,
    calculate_all_signatures,
)

# ── Configuration ──────────────────────────────────────────────────────────────
STREAMFLOW_PATH = r"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet"
CLIMATE_PATH = r"D:\processedOuts_feb2026\daymet_1980_2023.parquet"

# 10 hardcoded test gages (same as R smoke test seed=42 selection for reproducibility)
TEST_GAGES = [
    "01011000", "01030500", "01057000", "01073000", "01094400",
    "01118300", "01137500", "01152500", "01169000", "01181000",
]

# Expected minimum signature columns (non-climate)
MIN_EXPECTED_COLS = 400

# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print("=" * 70)
    print("PYTHON SMOKE TEST - Streamflow Signatures")
    print("=" * 70)

    t_start = time.perf_counter()

    # Phase 1: Load streamflow data
    print("\nPhase 1: Loading streamflow data...")
    if not Path(STREAMFLOW_PATH).exists():
        print(f"ERROR: Streamflow file not found: {STREAMFLOW_PATH}")
        print("Edit STREAMFLOW_PATH in this file to point to your data.")
        return 1

    df = read_parquet(STREAMFLOW_PATH, gage_ids=TEST_GAGES)
    df = add_water_year_columns(df)
    print(f"  Loaded {len(df):,} rows for {df['gage_id'].nunique()} gages")

    # Phase 2: Load and merge climate data (optional)
    has_climate = False
    if Path(CLIMATE_PATH).exists():
        print("\nPhase 2: Loading climate data...")
        climate = read_parquet(CLIMATE_PATH, gage_ids=TEST_GAGES)
        climate = add_water_year_columns(climate)
        df = df.merge(
            climate[["gage_id", "date", "PPT"]],
            on=["gage_id", "date"],
            how="left",
        )
        has_climate = True
        print(f"  Merged climate data. PPT non-null: {df['PPT'].notna().sum():,}")
    else:
        print("\nPhase 2: Skipping climate data (file not found)")

    # Phase 3: Process each gage
    print(f"\nPhase 3: Processing {len(TEST_GAGES)} gages...")
    results = []
    n_passed = 0
    n_skipped = 0

    for gage_id in TEST_GAGES:
        gage_data = df[df["gage_id"] == gage_id]
        if gage_data.empty:
            print(f"  [{gage_id}] Not found in data, skipping")
            n_skipped += 1
            continue

        qual_years, qualifies = filter_qualifying_years(gage_data)
        if not qualifies:
            print(f"  [{gage_id}] Only {len(qual_years)} qualifying years, skipping")
            n_skipped += 1
            continue

        gage_data = gage_data[gage_data["water_year"].isin(qual_years)]
        sigs = calculate_all_signatures(gage_data, has_climate=has_climate)
        sigs["gage_id"] = gage_id

        n_cols = len(sigs) - 1  # exclude gage_id
        n_nonnull = sum(1 for v in sigs.values() if v is not None and not (isinstance(v, float) and np.isnan(v)))
        print(f"  [{gage_id}] {len(qual_years)} years, {n_cols} signature cols, {n_nonnull} non-null")

        results.append(sigs)
        n_passed += 1

    # Phase 4: Validation
    print(f"\n{'=' * 70}")
    print("VALIDATION")
    print("=" * 70)

    all_ok = True

    if n_passed == 0:
        print("FAIL: No gages processed")
        return 1

    results_df = pd.DataFrame(results)
    n_sig_cols = len(results_df.columns) - 1  # exclude gage_id

    # Check 1: Minimum column count
    if n_sig_cols >= MIN_EXPECTED_COLS:
        print(f"  [PASS] Column count: {n_sig_cols} >= {MIN_EXPECTED_COLS}")
    else:
        print(f"  [FAIL] Column count: {n_sig_cols} < {MIN_EXPECTED_COLS}")
        all_ok = False

    # Check 2: Key signatures present
    key_sigs = ["Qann_mean", "flashinessRB_mean", "BFI_Eckhardt_mean", "D50_day_mean"]
    for sig in key_sigs:
        if sig in results_df.columns:
            vals = results_df[sig].dropna()
            print(f"  [PASS] {sig}: {len(vals)}/{n_passed} non-NA, range [{vals.min():.4f}, {vals.max():.4f}]")
        else:
            print(f"  [FAIL] {sig}: missing from output")
            all_ok = False

    # Check 3: Range validation
    if "Qann_mean" in results_df.columns:
        qann = results_df["Qann_mean"].dropna()
        if (qann >= 0).all() and (qann <= 100000).all():
            print(f"  [PASS] Qann_mean range check (0-100000)")
        else:
            print(f"  [FAIL] Qann_mean out of range")
            all_ok = False

    if "BFI_Eckhardt_mean" in results_df.columns:
        bfi = results_df["BFI_Eckhardt_mean"].dropna()
        if (bfi >= 0).all() and (bfi <= 1).all():
            print(f"  [PASS] BFI_Eckhardt_mean range check (0-1)")
        else:
            print(f"  [FAIL] BFI_Eckhardt_mean out of range")
            all_ok = False

    # Check 4: Climate signatures if available
    if has_climate:
        climate_sigs = ["annual_runoff_ratio_mean", "elasticity_static", "avg_storage_mean"]
        for sig in climate_sigs:
            if sig in results_df.columns and results_df[sig].notna().any():
                print(f"  [PASS] Climate sig {sig} present with data")
            else:
                print(f"  [WARN] Climate sig {sig} missing or all-NA")

    elapsed = time.perf_counter() - t_start

    print(f"\n{'=' * 70}")
    print(f"Gages: {n_passed} passed, {n_skipped} skipped")
    print(f"Total time: {elapsed:.1f}s")

    if all_ok:
        print("\nSTATUS: SMOKE TEST PASSED")
        return 0
    else:
        print("\nSTATUS: SMOKE TEST FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())

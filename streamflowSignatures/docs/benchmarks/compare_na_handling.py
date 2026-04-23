"""
Compare old (legacy) vs new (preprocessed) Julia benchmark outputs.
Validates that changes from the NA handling standardization match expectations.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import io

# Force UTF-8 output on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

BENCHMARKS_DIR = Path(__file__).parent

def load_signatures(path):
    """Load a signatures CSV and set gage_id as index."""
    df = pd.read_csv(path, dtype={"gage_id": str})
    df = df.set_index("gage_id")
    return df

def identity_r2(x, y):
    """R² of the identity line (y = x), not a fitted line."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    x_clean, y_clean = x[mask], y[mask]
    ss_res = np.sum((y_clean - x_clean) ** 2)
    ss_tot = np.sum((x_clean - np.mean(x_clean)) ** 2)
    if ss_tot == 0:
        return 1.0 if ss_res == 0 else 0.0
    return 1.0 - ss_res / ss_tot

def main():
    legacy_path = BENCHMARKS_DIR / "julia_signatures_legacy.csv"
    new_path = BENCHMARKS_DIR / "julia_signatures.csv"

    if not legacy_path.exists():
        print(f"ERROR: Legacy output not found: {legacy_path}")
        return 1
    if not new_path.exists():
        print(f"ERROR: New output not found: {new_path}")
        return 1

    print("=" * 70)
    print("NA HANDLING VALIDATION: Legacy vs New Julia Outputs")
    print("=" * 70)

    legacy = load_signatures(legacy_path)
    new = load_signatures(new_path)

    # Also load golden outputs for comparison
    golden_path = BENCHMARKS_DIR.parent.parent / "golden-outputs" / "streamflow_signatures_full_10feb2026.csv"
    has_golden = golden_path.exists()
    if has_golden:
        golden = load_signatures(golden_path)
        print(f"\nGolden outputs loaded: {golden.shape}")

    # --- 1. Gage Population Changes ---
    print(f"\n{'─' * 50}")
    print("1. GAGE POPULATION")
    print(f"{'─' * 50}")
    print(f"  Legacy gages:  {len(legacy)}")
    print(f"  New gages:     {len(new)}")

    legacy_ids = set(legacy.index)
    new_ids = set(new.index)
    gained = new_ids - legacy_ids
    lost = legacy_ids - new_ids
    common = legacy_ids & new_ids

    print(f"  Common:        {len(common)}")
    print(f"  Gained (new preprocessor rescued): {len(gained)}")
    print(f"  Lost (new preprocessor rejected):  {len(lost)}")

    if len(gained) > 0 and len(gained) <= 20:
        print(f"  Gained IDs: {sorted(gained)}")
    elif len(gained) > 20:
        print(f"  Gained IDs (first 20): {sorted(gained)[:20]}")

    if len(lost) > 0 and len(lost) <= 20:
        print(f"  Lost IDs: {sorted(lost)}")
    elif len(lost) > 20:
        print(f"  Lost IDs (first 20): {sorted(lost)[:20]}")

    # --- 2. Column Comparison ---
    print(f"\n{'─' * 50}")
    print("2. COLUMN SCHEMA")
    print(f"{'─' * 50}")
    legacy_cols = set(legacy.columns)
    new_cols = set(new.columns)
    print(f"  Legacy columns: {len(legacy_cols)}")
    print(f"  New columns:    {len(new_cols)}")
    if new_cols - legacy_cols:
        print(f"  New columns added: {new_cols - legacy_cols}")
    if legacy_cols - new_cols:
        print(f"  Columns removed: {legacy_cols - new_cols}")

    # --- 3. Per-Column R² Comparison (common gages) ---
    print(f"\n{'─' * 50}")
    print("3. PER-COLUMN IDENTITY R² (common gages)")
    print(f"{'─' * 50}")

    # Identify signature columns (exclude metadata)
    meta_cols = {"latitude", "longitude", "basin_area", "gage_type", "area_normalized",
                 "num_water_years", "start_water_year", "end_water_year",
                 "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
                 "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
                 "HYDRO_DISTURB_INDX", "CLASS", "RHBN", "REGULATED",
                 "human_interference_class"}
    flag_cols = {c for c in legacy.columns if c.startswith("flagged_")}
    sig_cols = sorted(set(legacy.columns) & set(new.columns) - meta_cols - flag_cols)

    results = []
    for col in sig_cols:
        if legacy[col].dtype == "object" or new[col].dtype == "object":
            continue
        old_vals = legacy.loc[list(common), col].astype(float).values
        new_vals = new.loc[list(common), col].astype(float).values
        r2 = identity_r2(old_vals, new_vals)
        n_valid = np.sum(np.isfinite(old_vals) & np.isfinite(new_vals))
        n_old_na = np.sum(~np.isfinite(old_vals))
        n_new_na = np.sum(~np.isfinite(new_vals))
        results.append({
            "column": col,
            "r2": r2,
            "n_valid": int(n_valid),
            "old_na": int(n_old_na),
            "new_na": int(n_new_na),
            "na_diff": int(n_new_na) - int(n_old_na),
        })

    results_df = pd.DataFrame(results)

    # Summary buckets
    perfect = (results_df["r2"] >= 0.999).sum()
    good = ((results_df["r2"] >= 0.99) & (results_df["r2"] < 0.999)).sum()
    fair = ((results_df["r2"] >= 0.95) & (results_df["r2"] < 0.99)).sum()
    poor = (results_df["r2"] < 0.95).sum()
    na_r2 = results_df["r2"].isna().sum()

    print(f"\n  Total signature columns compared: {len(results_df)}")
    print(f"  Perfect (R² >= 0.999):  {perfect}")
    print(f"  Good (0.99 - 0.999):    {good}")
    print(f"  Fair (0.95 - 0.99):     {fair}")
    print(f"  Poor (R² < 0.95):       {poor}")
    print(f"  N/A (too few values):   {na_r2}")

    # --- 4. Changed Columns (focus on non-perfect) ---
    print(f"\n{'─' * 50}")
    print("4. COLUMNS WITH R² < 0.999 (changed by NA handling)")
    print(f"{'─' * 50}")

    changed = results_df[results_df["r2"] < 0.999].sort_values("r2")
    if len(changed) == 0:
        print("  All columns are perfect matches!")
    else:
        print(f"\n  {'Column':<45} {'R²':>8} {'N':>6} {'Old NA':>7} {'New NA':>7} {'NA Δ':>6}")
        print(f"  {'─'*45} {'─'*8} {'─'*6} {'─'*7} {'─'*7} {'─'*6}")
        for _, row in changed.iterrows():
            print(f"  {row['column']:<45} {row['r2']:>8.4f} {row['n_valid']:>6} {row['old_na']:>7} {row['new_na']:>7} {row['na_diff']:>+6}")

    # --- 5. Focus: fillna(0) Affected Signatures ---
    print(f"\n{'─' * 50}")
    print("5. FILLNA(0) AFFECTED SIGNATURES (timing, storage, qp_seasonality)")
    print(f"{'─' * 50}")

    fillna_prefixes = ["D5_day", "D10_day", "D20_day", "D30_day", "D40_day",
                       "D50_day", "D60_day", "D70_day", "D80_day", "D90_day",
                       "D95_day", "D25_to_D75", "Dmax",
                       "avg_storage", "qp_slope_sd", "qp_bimodality"]
    fillna_cols = [c for c in results_df["column"] if any(c.startswith(p) for p in fillna_prefixes)]

    if fillna_cols:
        fillna_df = results_df[results_df["column"].isin(fillna_cols)].sort_values("r2")
        print(f"\n  {'Column':<45} {'R²':>8} {'Old NA':>7} {'New NA':>7} {'NA Δ':>6}")
        print(f"  {'─'*45} {'─'*8} {'─'*7} {'─'*7} {'─'*6}")
        for _, row in fillna_df.iterrows():
            marker = " ***" if row["r2"] < 0.999 else ""
            print(f"  {row['column']:<45} {row['r2']:>8.4f} {row['old_na']:>7} {row['new_na']:>7} {row['na_diff']:>+6}{marker}")

    # --- 6. Num Water Years Comparison ---
    print(f"\n{'─' * 50}")
    print("6. NUM_WATER_YEARS COMPARISON (common gages)")
    print(f"{'─' * 50}")
    if "num_water_years" in legacy.columns and "num_water_years" in new.columns:
        old_nwy = legacy.loc[list(common), "num_water_years"].astype(float)
        new_nwy = new.loc[list(common), "num_water_years"].astype(float)
        diff = new_nwy - old_nwy
        print(f"  Gages with more years (rescued):   {(diff > 0).sum()}")
        print(f"  Gages with same years:             {(diff == 0).sum()}")
        print(f"  Gages with fewer years (rejected): {(diff < 0).sum()}")
        if (diff != 0).any():
            print(f"  Mean year change: {diff.mean():+.2f}")
            print(f"  Max gained: {diff.max():.0f}, Max lost: {diff.min():.0f}")

    # --- 7. Golden Output Comparison (if available) ---
    if has_golden:
        print(f"\n{'─' * 50}")
        print("7. NEW JULIA vs GOLDEN OUTPUTS (R historical, Feb 2026)")
        print(f"{'─' * 50}")

        golden_ids = set(golden.index)
        common_golden = new_ids & golden_ids
        print(f"  Common gages: {len(common_golden)}")

        golden_sig_cols = sorted(set(golden.columns) & set(new.columns) - meta_cols - flag_cols)
        golden_results = []
        for col in golden_sig_cols:
            if golden[col].dtype == "object" or new[col].dtype == "object":
                continue
            g_vals = golden.loc[list(common_golden), col].astype(float).values
            n_vals = new.loc[list(common_golden), col].astype(float).values
            r2 = identity_r2(g_vals, n_vals)
            golden_results.append({"column": col, "r2": r2})

        golden_df = pd.DataFrame(golden_results)
        g_perfect = (golden_df["r2"] >= 0.999).sum()
        g_good = ((golden_df["r2"] >= 0.99) & (golden_df["r2"] < 0.999)).sum()
        g_poor = (golden_df["r2"] < 0.99).sum()
        print(f"  Perfect (R² >= 0.999): {g_perfect}")
        print(f"  Good (0.99 - 0.999):   {g_good}")
        print(f"  Poor (R² < 0.99):      {g_poor}")

    # Save detailed results
    output_path = BENCHMARKS_DIR / "na_handling_comparison.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nDetailed results saved to: {output_path}")

    print(f"\n{'=' * 70}")
    print("VALIDATION COMPLETE")
    print(f"{'=' * 70}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

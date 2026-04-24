"""
Compare rpkg output against canonical Julia, Python, and R benchmark outputs.

Uses the same R² identity-line metric as compare_three_way.py.
"""

import sys
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import stats

BENCHMARK_DIR = Path(__file__).parent

META_COLS = {
    "gage_id", "gage_id_metadata", "latitude", "longitude", "basin_area",
    "basin_area_km2", "gage_type", "num_water_years", "start_year", "end_year",
    "start_water_year", "end_water_year", "country",
    "processing_status", "drainage_area_km2", "area_normalized",
    "human_interference_class", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
    "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX",
    "CLASS", "RHBN", "REGULATED",
}


def normalize_col(col):
    col = col.replace("-", "_")
    col = col.replace(".", "_")
    col = col.replace("FDC_all", "FDCall")
    col = col.replace("FDC_90th", "FDC90th")
    col = col.replace("FDC_mid", "FDCmid")
    return col


def is_sig_col(col):
    if col in META_COLS:
        return False
    if col.startswith("flagged_for_") or col.startswith("flagged_"):
        return False
    return True


def r2_identity(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    n = mask.sum()
    if n < 5:
        return np.nan, n
    xv, yv = x[mask], y[mask]
    ss_res = np.sum((yv - xv) ** 2)
    ss_tot = np.sum((yv - np.mean(yv)) ** 2)
    if ss_tot == 0:
        return (1.0 if ss_res == 0 else np.nan), n
    return 1.0 - ss_res / ss_tot, n


def spearman_corr(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    n = mask.sum()
    if n < 5:
        return np.nan, n
    rho, _ = stats.spearmanr(x[mask], y[mask])
    return rho, n


STAT_SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                 "_mk_rho", "_mk_pval", "_mean", "_median"]


def get_base_metric(col):
    for suf in STAT_SUFFIXES:
        if col.endswith(suf):
            return col[:-len(suf)]
    return col


def categorize_metric(base):
    if base.startswith(("Qann", "Qwin", "Qspr", "Qsum", "Qfal")):
        return "Flow Volumes"
    if base.startswith("Q") and len(base) > 1 and any(c.isdigit() for c in base[1:3]):
        return "Flow Percentiles"
    if base.startswith("FDC"):
        return "FDC"
    if base.startswith("BFI"):
        return "Baseflow"
    if base.startswith(("log_a", "b_", "concavity", "n_recession", "alpha_linear",
                        "recession_alpha")):
        return "Recession"
    if base.startswith(("n_high", "n_low", "dur_high", "dur_low", "TQmean", "Flow_Reversal")):
        return "Pulse Metrics"
    if base.startswith("flashiness"):
        return "Flashiness"
    if base.startswith("negative_ann"):
        return "Negative Days"
    if base.startswith("D") and ("_day" in base or base in ("D25_to_D75", "Dmax")):
        return "Flow Timing"
    if "runoff_ratio" in base:
        return "Runoff Ratios"
    if base.startswith("elasticity"):
        return "Elasticity"
    if base.startswith("qp_"):
        return "Q-P Seasonality"
    if "storage" in base:
        return "Average Storage"
    return "Other"


def compare_pair(df_a, df_b, label_a, label_b, common_gages):
    """Compare two DataFrames on common gages and signature columns."""
    a = df_a[df_a["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()
    b = df_b[df_b["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()

    # Normalize column names
    a_norm = {normalize_col(c): c for c in a.columns if is_sig_col(normalize_col(c))}
    b_norm = {normalize_col(c): c for c in b.columns if is_sig_col(normalize_col(c))}

    common_cols = sorted(set(a_norm.keys()) & set(b_norm.keys()))

    results = []
    for nc in common_cols:
        ac, bc = a_norm[nc], b_norm[nc]
        x = pd.to_numeric(a[ac], errors="coerce").values
        y = pd.to_numeric(b[bc], errors="coerce").values
        r2, n = r2_identity(x, y)
        rho, _ = spearman_corr(x, y)
        na_a = int((~np.isfinite(x)).sum())
        na_b = int((~np.isfinite(y)).sum())
        results.append({
            "column": nc,
            "r2": r2,
            "rho": rho,
            "n": n,
            "na_a": na_a,
            "na_b": na_b,
            "category": categorize_metric(get_base_metric(nc)),
        })

    return pd.DataFrame(results)


def print_summary(res_df, label):
    """Print summary statistics for a comparison."""
    valid = res_df["r2"].dropna()
    n_total = len(valid)
    if n_total == 0:
        print(f"  No valid comparisons for {label}")
        return

    perfect = (valid >= 0.999).sum()
    good = ((valid >= 0.99) & (valid < 0.999)).sum()
    poor = (valid < 0.99).sum()

    print(f"\n  {label}:")
    print(f"    Columns compared:  {n_total}")
    print(f"    Perfect (>=0.999):  {perfect}")
    print(f"    Good (0.99-0.999): {good}")
    print(f"    Poor (<0.99):      {poor}")
    print(f"    Mean R²:           {valid.mean():.6f}")
    print(f"    Median R²:         {valid.median():.6f}")
    print(f"    Min R²:            {valid.min():.6f}")

    if poor > 0:
        print(f"\n    Poor columns (<0.99 R²):")
        poor_df = res_df[res_df["r2"] < 0.99].sort_values("r2")
        for _, row in poor_df.iterrows():
            print(f"      {row['column']:40s}  R²={row['r2']:.4f}  rho={row['rho']:.4f}  "
                  f"n={row['n']:.0f}  cat={row['category']}")


def print_category_summary(res_df, label):
    """Print per-category summary."""
    print(f"\n  {label} — Per-category breakdown:")
    cats = res_df.groupby("category").agg(
        n_cols=("r2", "count"),
        mean_r2=("r2", "mean"),
        min_r2=("r2", "min"),
        n_poor=("r2", lambda x: (x < 0.99).sum()),
    ).sort_values("min_r2")

    print(f"    {'Category':20s} {'Cols':>5s} {'Mean R²':>9s} {'Min R²':>9s} {'Poor':>5s}")
    for cat, row in cats.iterrows():
        flag = " ***" if row["n_poor"] > 0 else ""
        print(f"    {cat:20s} {row['n_cols']:5.0f} {row['mean_r2']:9.6f} "
              f"{row['min_r2']:9.6f} {row['n_poor']:5.0f}{flag}")


def main():
    print("=" * 72)
    print("rpkg Cross-Language Comparison")
    print("=" * 72)

    # Load rpkg output
    rpkg_path = BENCHMARK_DIR / "rpkg_signatures.csv"
    if not rpkg_path.exists():
        print(f"ERROR: rpkg output not found: {rpkg_path}")
        sys.exit(1)

    rpkg = pd.read_csv(rpkg_path, low_memory=False)
    rpkg["gage_id"] = rpkg["gage_id"].astype(str).str.strip()
    print(f"\nrpkg:   {len(rpkg)} gages, {len(rpkg.columns)} columns")

    # Load other benchmarks
    pairs = []
    for name, filename in [("R (legacy)", "r_signatures.csv"),
                            ("Python", "python_signatures.csv"),
                            ("Julia", "julia_signatures.csv")]:
        path = BENCHMARK_DIR / filename
        if path.exists():
            df = pd.read_csv(path, low_memory=False)
            df["gage_id"] = df["gage_id"].astype(str).str.strip()
            print(f"{name:15s}: {len(df)} gages, {len(df.columns)} columns")
            pairs.append((name, df))
        else:
            print(f"{name:15s}: NOT FOUND ({path})")

    if not pairs:
        print("ERROR: No benchmark files found to compare against")
        sys.exit(1)

    # Compare rpkg against each
    all_results = {}
    for name, other_df in pairs:
        common = set(rpkg["gage_id"]) & set(other_df["gage_id"])
        print(f"\n{'=' * 72}")
        print(f"rpkg vs {name}: {len(common)} common gages")
        print(f"{'=' * 72}")

        if len(common) < 10:
            print(f"  Too few common gages ({len(common)}), skipping")
            continue

        res = compare_pair(rpkg, other_df, "rpkg", name, common)
        all_results[name] = res

        print_summary(res, f"rpkg vs {name}")
        print_category_summary(res, f"rpkg vs {name}")

    # Save detailed results
    if all_results:
        out_path = BENCHMARK_DIR / "rpkg_comparison.csv"
        frames = []
        for name, res in all_results.items():
            res = res.copy()
            res["comparison"] = f"rpkg_vs_{name.replace(' ', '_').replace('(', '').replace(')', '')}"
            frames.append(res)
        combined = pd.concat(frames, ignore_index=True)
        combined.to_csv(out_path, index=False)
        print(f"\nDetailed results saved to {out_path}")

    # Overall verdict
    print(f"\n{'=' * 72}")
    print("OVERALL VERDICT")
    print(f"{'=' * 72}")
    for name, res in all_results.items():
        valid = res["r2"].dropna()
        perfect = (valid >= 0.999).sum()
        poor = (valid < 0.99).sum()
        print(f"  rpkg vs {name:15s}: {perfect}/{len(valid)} perfect, "
              f"{poor} poor, min R²={valid.min():.4f}")


if __name__ == "__main__":
    main()

"""
Compare rpkg benchmark output against Golden Julia (canonical) reference.

Produces:
  - Console summary
  - docs/benchmarks/rpkg_vs_golden_julia_comparison.csv   (per-column detail)
  - docs/benchmarks/rpkg_vs_golden_julia_summary.md        (high-level report)

Usage:
    python docs/benchmarks/compare_rpkg_vs_golden_julia.py
"""

import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

# Paths overridable via ENV (port campaign, 2026-08-24): COMPARE_REFERENCE_CSV
# points at any Julia reference run; COMPARE_CANDIDATE_CSV at the port output;
# COMPARE_OUTPUT_DIR redirects both artifacts (one-folder convention).
GOLDEN_JULIA_PATH = Path(os.environ.get("COMPARE_REFERENCE_CSV",
    str(PROJECT_ROOT / "golden-outputs" / "streamflow_signatures_julia_apr2026.csv")))
RPKG_PATH = Path(os.environ.get("COMPARE_CANDIDATE_CSV", str(SCRIPT_DIR / "rpkg_signatures.csv")))
_OUT_DIR = Path(os.environ.get("COMPARE_OUTPUT_DIR", str(SCRIPT_DIR)))
OUTPUT_CSV = _OUT_DIR / "rpkg_vs_golden_julia_comparison.csv"
OUTPUT_MD = _OUT_DIR / "rpkg_vs_golden_julia_summary.md"

# ── Helpers ──────────────────────────────────────────────────────────

META_COLS = {
    "gage_id", "gage_id_metadata", "latitude", "longitude", "basin_area",
    "basin_area_km2", "gage_type", "num_water_years", "start_year", "end_year",
    "start_water_year", "end_water_year", "country",
    "processing_status", "drainage_area_km2", "area_normalized",
    "human_interference_class", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
    "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX",
    "CLASS", "RHBN", "REGULATED",
    "ice_affected_days_total",
}

QAQC_PREFIX = "flagged_for_"

STAT_SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                 "_mk_rho", "_mk_pval", "_mean", "_median"]


def is_signature_col(col):
    if col in META_COLS:
        return False
    if col.startswith(QAQC_PREFIX):
        return False
    return True


def normalize_col(col):
    col = col.replace("-", "_")
    col = col.replace(".", "_")
    col = col.replace("FDC_all", "FDCall")
    col = col.replace("FDC_90th", "FDC90th")
    col = col.replace("FDC_mid", "FDCmid")
    return col


def r2_identity(x, y):
    """R-squared of the identity line y = x."""
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
    rho, _ = sp_stats.spearmanr(x[mask], y[mask])
    return rho, n


def count_na_mismatch(x, y):
    x_na = ~np.isfinite(x)
    y_na = ~np.isfinite(y)
    return int((x_na & ~y_na).sum() + (~x_na & y_na).sum())


def get_base_metric(col):
    for suf in STAT_SUFFIXES:
        if col.endswith(suf):
            return col[:-len(suf)]
    return col


def get_stat_suffix(col):
    for suf in STAT_SUFFIXES:
        if col.endswith(suf):
            return suf.lstrip("_")
    return "scalar"


def categorize_metric(base):
    if base.startswith(("Qann", "Qwin", "Qspr", "Qsum", "Qfal")):
        return "Flow Volumes"
    if base.startswith("Q") and len(base) > 1 and any(c.isdigit() for c in base[1:4]):
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
    if base.startswith(("D1_", "D5_", "D10", "D20", "D25", "D30", "D40", "D50",
                        "D60", "D70", "D80", "D90", "D95", "D99", "Dmax")):
        return "Flow Timing"
    if "runoff_ratio" in base:
        return "Runoff Ratios"
    if "elasticity" in base:
        return "Elasticity"
    if base.startswith("qp_"):
        return "Q-P Seasonality"
    if base.startswith("avg_storage"):
        return "Storage"
    if base.startswith("negative"):
        return "Negative Flow"
    return "Other"


def classify_r2(r2_val):
    if np.isnan(r2_val):
        return "N/A"
    if r2_val >= 0.999:
        return "Perfect"
    if r2_val >= 0.99:
        return "Good"
    if r2_val >= 0.95:
        return "Poor"
    if r2_val >= 0.9:
        return "Low"
    if r2_val >= 0.5:
        return "Very Low"
    return "Extremely Low"


TIERS = [
    ("Perfect",       "R2 >= 0.999"),
    ("Good",          "0.99 <= R2 < 0.999"),
    ("Poor",          "0.95 <= R2 < 0.99"),
    ("Low",           "0.90 <= R2 < 0.95"),
    ("Very Low",      "0.50 <= R2 < 0.90"),
    ("Extremely Low", "R2 < 0.50"),
]


def tier_counts(r2_series):
    counts = {name: 0 for name, _ in TIERS}
    for v in r2_series:
        t = classify_r2(v)
        if t in counts:
            counts[t] += 1
    return counts


# ── Main comparison ───────────────────────────────────────────────────

def main():
    print("=" * 80)
    print("RPKG vs GOLDEN JULIA (CANONICAL) REFERENCE")
    print("=" * 80)

    if not GOLDEN_JULIA_PATH.exists():
        print(f"\nERROR: Golden Julia file not found: {GOLDEN_JULIA_PATH}")
        print("Run `julia docs/benchmarks/run_julia_benchmark.jl` and copy output to golden-outputs/")
        return 1
    if not RPKG_PATH.exists():
        print(f"\nERROR: rpkg benchmark output not found: {RPKG_PATH}")
        print("Run `Rscript docs/benchmarks/run_rpkg_benchmark.R` first")
        return 1

    print("\nLoading Golden Julia output...")
    jl_df = pd.read_csv(GOLDEN_JULIA_PATH, low_memory=False)
    jl_df["gage_id"] = jl_df["gage_id"].astype(str).str.strip()
    print(f"  {jl_df.shape[0]:,} gages x {jl_df.shape[1]:,} cols")

    print("Loading rpkg benchmark output...")
    rpkg_df = pd.read_csv(RPKG_PATH, low_memory=False)
    rpkg_df["gage_id"] = rpkg_df["gage_id"].astype(str).str.strip()
    print(f"  {rpkg_df.shape[0]:,} gages x {rpkg_df.shape[1]:,} cols")

    # Common gages
    jl_gages = set(jl_df["gage_id"])
    rpkg_gages = set(rpkg_df["gage_id"])
    common_gages = sorted(jl_gages & rpkg_gages)
    print(f"\n  Common gages: {len(common_gages):,}")
    print(f"  Only in Golden Julia: {len(jl_gages - rpkg_gages):,}")
    print(f"  Only in rpkg: {len(rpkg_gages - jl_gages):,}")

    jl = jl_df[jl_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()
    rp = rpkg_df[rpkg_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()

    # Column mapping
    jl_sig_cols = {c for c in jl.columns if is_signature_col(c)}
    rp_sig_cols = {c for c in rp.columns if is_signature_col(c)}

    jl_norm = {normalize_col(c): c for c in jl_sig_cols}
    rp_norm = {normalize_col(c): c for c in rp_sig_cols}

    common_norms = sorted(set(jl_norm) & set(rp_norm))
    jl_only_norms = sorted(set(jl_norm) - set(rp_norm))
    rp_only_norms = sorted(set(rp_norm) - set(jl_norm))

    print(f"\n  Common signature columns: {len(common_norms):,}")
    print(f"  Only in Golden Julia: {len(jl_only_norms)}")
    print(f"  Only in rpkg: {len(rp_only_norms)}")

    # ── Per-column comparison ──
    print(f"\n{'=' * 80}")
    print("COMPUTING PER-COLUMN METRICS...")
    print(f"{'=' * 80}")

    results = []
    for norm_col in common_norms:
        jl_col = jl_norm[norm_col]
        rp_col = rp_norm[norm_col]

        jl_vals = pd.to_numeric(jl[jl_col], errors="coerce").values
        rp_vals = pd.to_numeric(rp[rp_col], errors="coerce").values

        r2, n_valid = r2_identity(jl_vals, rp_vals)
        rho, _ = spearman_corr(jl_vals, rp_vals)
        na_mm = count_na_mismatch(jl_vals, rp_vals)

        jl_na = int((~np.isfinite(jl_vals)).sum())
        rp_na = int((~np.isfinite(rp_vals)).sum())

        mask = np.isfinite(jl_vals) & np.isfinite(rp_vals)
        if mask.sum() > 0:
            diff = rp_vals[mask] - jl_vals[mask]
            abs_diff = np.abs(diff)
            mad = np.mean(abs_diff)
            max_diff = np.max(abs_diff)
            median_diff = np.median(diff)
        else:
            mad = max_diff = median_diff = np.nan

        base = get_base_metric(norm_col)
        stat = get_stat_suffix(norm_col)
        cat = categorize_metric(base)

        results.append({
            "column": norm_col,
            "jl_col": jl_col,
            "rpkg_col": rp_col,
            "category": cat,
            "base_metric": base,
            "stat_type": stat,
            "r2_identity": r2,
            "spearman_rho": rho,
            "n_valid_pairs": n_valid,
            "na_mismatch": na_mm,
            "jl_na_count": jl_na,
            "rpkg_na_count": rp_na,
            "mean_abs_diff": mad,
            "max_abs_diff": max_diff,
            "median_diff": median_diff,
        })

    res = pd.DataFrame(results)
    if len(res) == 0:
        print("\nNo common signature columns found. Nothing to compare.")
        return 1
    res["agreement_tier"] = res["r2_identity"].apply(classify_r2)
    res.to_csv(OUTPUT_CSV, index=False)
    print(f"\nPer-column CSV saved to: {OUTPUT_CSV}")

    # ── Generate markdown summary ──
    generate_summary_report(res, common_norms, jl_only_norms, rp_only_norms,
                            len(common_gages), jl_df.shape[0], rpkg_df.shape[0],
                            jl_norm, rp_norm)

    # ── Console summary ──
    print_console_summary(res)

    return 0


def print_console_summary(res):
    valid = res.dropna(subset=["r2_identity"])

    print(f"\n{'=' * 80}")
    print("OVERALL ALIGNMENT: rpkg vs Golden Julia (Canonical)")
    print(f"{'=' * 80}")

    r2 = valid["r2_identity"]
    print(f"\n  Columns compared: {len(valid)}")
    print(f"  Mean R2:   {r2.mean():.6f}")
    print(f"  Median R2: {r2.median():.6f}")
    print(f"  Min R2:    {r2.min():.6f}")

    tc = tier_counts(r2)
    total = len(valid)
    print()
    for name, desc in TIERS:
        c = tc[name]
        print(f"  {name:<16} ({desc:<24}): {c:>4} ({100*c/total:.1f}%)")

    # By category
    hdr_tiers = "".join(f"{t[0]:>8}" for t in TIERS)
    print(f"\n  {'Category':<20} {'Cols':>5}{hdr_tiers} {'Mean R2':>10} {'Min R2':>10}")
    print(f"  {'-'*105}")
    for cat in sorted(res["category"].unique()):
        cat_df = res[res["category"] == cat].dropna(subset=["r2_identity"])
        n = len(cat_df)
        r2c = cat_df["r2_identity"]
        tc_cat = tier_counts(r2c)
        vals = "".join(f"{tc_cat[t[0]]:>8}" for t in TIERS)
        print(f"  {cat:<20} {n:>5}{vals} {r2c.mean():>10.6f} {r2c.min():>10.4f}")

    # Worst columns
    poor_df = res[res["r2_identity"] < 0.99].sort_values("r2_identity")
    if len(poor_df) > 0:
        print(f"\n  Worst columns (R2 < 0.99):")
        for _, row in poor_df.head(20).iterrows():
            r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
            tier = classify_r2(row['r2_identity'])
            print(f"    {row['column']:<45} R2={r2v}  [{tier}]  MAD={row['mean_abs_diff']:.4f}")


def generate_summary_report(res, common_norms, jl_only_norms, rp_only_norms,
                            n_common_gages, n_jl_total, n_rpkg_total,
                            jl_norm, rp_norm):
    valid = res.dropna(subset=["r2_identity"])
    r2 = valid["r2_identity"]

    with open(OUTPUT_MD, "w", encoding="utf-8") as f:
        f.write("# rpkg vs Golden Julia (Canonical): Comparison Report\n\n")
        f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("## Input Files\n\n")
        f.write("| Dataset | File | Gages | Columns |\n")
        f.write("|---------|------|-------|---------|\n")
        f.write(f"| Julia reference | `{GOLDEN_JULIA_PATH}` | {n_jl_total:,} | {len(jl_norm):,} |\n")
        f.write(f"| rpkg benchmark | `{RPKG_PATH}` | {n_rpkg_total:,} | — |\n")
        f.write(f"\n**Common gages**: {n_common_gages:,}\n")
        f.write(f"**Common signature columns**: {len(common_norms):,}\n\n")

        # High-level summary
        total = len(valid)
        tc = tier_counts(r2)

        f.write("## Agreement Summary\n\n")
        f.write("| Tier | Threshold | Count | % |\n")
        f.write("|------|-----------|-------|---|\n")
        for name, desc in TIERS:
            c = tc[name]
            f.write(f"| {name} | {desc} | {c} | {100*c/total:.1f}% |\n")
        f.write(f"| **Total** | | **{total}** | **100%** |\n\n")

        f.write(f"**Mean R2**: {r2.mean():.6f} | **Median R2**: {r2.median():.6f} | **Min R2**: {r2.min():.6f}\n\n")

        # By category
        f.write("## Agreement by Category\n\n")
        tier_hdrs = " | ".join(t[0] for t in TIERS)
        f.write(f"| Category | Cols | {tier_hdrs} | Mean R2 | Min R2 |\n")
        f.write(f"|----------|------|-" + "-|-".join("-" * len(t[0]) for t in TIERS) + "-|---------|--------|\n")
        for cat in sorted(res["category"].unique()):
            cat_df = res[res["category"] == cat].dropna(subset=["r2_identity"])
            n = len(cat_df)
            if n == 0:
                continue
            r2c = cat_df["r2_identity"]
            tc_cat = tier_counts(r2c)
            tier_vals = " | ".join(f"{tc_cat[t[0]]}" for t in TIERS)
            f.write(f"| {cat} | {n} | {tier_vals} | {r2c.mean():.6f} | {r2c.min():.4f} |\n")
        f.write("\n")

        # Worst columns
        poor_df = res[res["r2_identity"] < 0.99].sort_values("r2_identity")
        f.write(f"## Columns with R2 < 0.99 ({len(poor_df)} columns)\n\n")
        if len(poor_df) > 0:
            f.write("| Column | Category | Stat | R2 | Spearman | MAD | Julia NAs | rpkg NAs | NA Mismatch |\n")
            f.write("|--------|----------|------|----|----------|-----|-----------|----------|-------------|\n")
            for _, row in poor_df.iterrows():
                r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
                rho = f"{row['spearman_rho']:.4f}" if np.isfinite(row['spearman_rho']) else "N/A"
                mad = f"{row['mean_abs_diff']:.4f}" if np.isfinite(row['mean_abs_diff']) else "N/A"
                f.write(f"| `{row['column']}` | {row['category']} | {row['stat_type']} | "
                        f"{r2v} | {rho} | {mad} | {row['jl_na_count']} | {row['rpkg_na_count']} | {row['na_mismatch']} |\n")
        else:
            f.write("All columns have R2 >= 0.99.\n")
        f.write("\n")

        f.write("---\n")
        f.write("*Generated by `docs/benchmarks/compare_rpkg_vs_golden_julia.py`*\n")

    print(f"\nMarkdown report saved to: {OUTPUT_MD}")


if __name__ == "__main__":
    sys.exit(main())

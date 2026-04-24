"""
Compare updated Julia (post-Section 3) output against Golden R reference.

Produces:
  - Console summary
  - docs/benchmarks/julia_vs_golden_r_comparison.csv   (per-column detail)
  - docs/benchmarks/julia_vs_golden_r_summary.md        (high-level report)

Usage:
    python docs/benchmarks/compare_julia_vs_golden_r.py
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

GOLDEN_R_PATH = PROJECT_ROOT / "golden-outputs" / "streamflow_signatures_full_10feb2026.csv"
JULIA_PATH = SCRIPT_DIR / "julia_signatures.csv"

OUTPUT_CSV = SCRIPT_DIR / "julia_vs_golden_r_comparison.csv"
OUTPUT_MD = SCRIPT_DIR / "julia_vs_golden_r_summary.md"

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


# Columns renamed between golden R (Feb 2026) and Julia post-Section 3 (Apr 2026)
# Golden R has "elasticity_*", Julia has "elasticity_rolling_*"
RENAMED_COLS = {}
for suf in STAT_SUFFIXES:
    old = f"elasticity{suf}"
    new = f"elasticity_rolling{suf}"
    RENAMED_COLS[normalize_col(old)] = normalize_col(new)


# ── Main comparison ───────────────────────────────────────────────────

def main():
    print("=" * 80)
    print("JULIA (POST-SECTION 3) vs GOLDEN R REFERENCE")
    print("=" * 80)

    # Load data
    print("\nLoading Golden R output...")
    r_df = pd.read_csv(GOLDEN_R_PATH, low_memory=False)
    r_df["gage_id"] = r_df["gage_id"].astype(str).str.strip()
    print(f"  {r_df.shape[0]:,} gages x {r_df.shape[1]:,} cols")

    print("Loading Julia post-Section 3 output...")
    jl_df = pd.read_csv(JULIA_PATH, low_memory=False)
    jl_df["gage_id"] = jl_df["gage_id"].astype(str).str.strip()
    print(f"  {jl_df.shape[0]:,} gages x {jl_df.shape[1]:,} cols")

    # Common gages
    r_gages = set(r_df["gage_id"])
    jl_gages = set(jl_df["gage_id"])
    common_gages = sorted(r_gages & jl_gages)
    print(f"\n  Common gages: {len(common_gages):,}")
    print(f"  Only in Golden R: {len(r_gages - jl_gages):,}")
    print(f"  Only in Julia: {len(jl_gages - r_gages):,}")

    r = r_df[r_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()
    j = jl_df[jl_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()

    # Column mapping
    r_sig_cols = {c for c in r.columns if is_signature_col(c)}
    jl_sig_cols = {c for c in j.columns if is_signature_col(c)}

    r_norm = {normalize_col(c): c for c in r_sig_cols}
    jl_norm = {normalize_col(c): c for c in jl_sig_cols}

    # Handle renames: if a golden-R normalized name maps to a Julia renamed name
    r_norm_resolved = {}
    for norm, orig in r_norm.items():
        if norm in jl_norm:
            r_norm_resolved[norm] = orig
        elif norm in RENAMED_COLS and RENAMED_COLS[norm] in jl_norm:
            # Map old name to new name for matching
            r_norm_resolved[RENAMED_COLS[norm]] = orig
        else:
            r_norm_resolved[norm] = orig

    common_norms = sorted(set(r_norm_resolved) & set(jl_norm))
    r_only_norms = sorted(set(r_norm_resolved) - set(jl_norm))
    jl_only_norms = sorted(set(jl_norm) - set(r_norm_resolved))

    print(f"\n  Common signature columns: {len(common_norms):,}")
    print(f"  Only in Golden R: {len(r_only_norms)}")
    print(f"  Only in Julia (new): {len(jl_only_norms)}")

    if jl_only_norms:
        print(f"\n  Julia-only columns ({len(jl_only_norms)}):")
        for c in jl_only_norms:
            print(f"    + {jl_norm[c]}")

    if r_only_norms:
        print(f"\n  Golden R-only columns ({len(r_only_norms)}):")
        for c in r_only_norms[:10]:
            print(f"    - {r_norm_resolved[c]}")
        if len(r_only_norms) > 10:
            print(f"    ... and {len(r_only_norms) - 10} more")

    # ── Per-column comparison ──
    print(f"\n{'=' * 80}")
    print("COMPUTING PER-COLUMN METRICS...")
    print(f"{'=' * 80}")

    results = []
    for norm_col in common_norms:
        r_col = r_norm_resolved[norm_col]
        jl_col = jl_norm[norm_col]

        r_vals = pd.to_numeric(r[r_col], errors="coerce").values
        jl_vals = pd.to_numeric(j[jl_col], errors="coerce").values

        r2, n_valid = r2_identity(r_vals, jl_vals)
        rho, _ = spearman_corr(r_vals, jl_vals)
        na_mm = count_na_mismatch(r_vals, jl_vals)

        r_na = int((~np.isfinite(r_vals)).sum())
        jl_na = int((~np.isfinite(jl_vals)).sum())

        # Per-column distribution stats
        mask = np.isfinite(r_vals) & np.isfinite(jl_vals)
        if mask.sum() > 0:
            diff = jl_vals[mask] - r_vals[mask]
            abs_diff = np.abs(diff)
            mad = np.mean(abs_diff)
            max_diff = np.max(abs_diff)
            median_diff = np.median(diff)
            r_mean = np.mean(r_vals[mask])
            r_median = np.median(r_vals[mask])
            r_sd = np.std(r_vals[mask], ddof=1) if mask.sum() > 1 else 0.0
            jl_mean = np.mean(jl_vals[mask])
            jl_median = np.median(jl_vals[mask])
            jl_sd = np.std(jl_vals[mask], ddof=1) if mask.sum() > 1 else 0.0
        else:
            mad = max_diff = median_diff = np.nan
            r_mean = r_median = r_sd = np.nan
            jl_mean = jl_median = jl_sd = np.nan

        base = get_base_metric(norm_col)
        stat = get_stat_suffix(norm_col)
        cat = categorize_metric(base)

        results.append({
            "column": norm_col,
            "r_col": r_col,
            "jl_col": jl_col,
            "category": cat,
            "base_metric": base,
            "stat_type": stat,
            "r2_identity": r2,
            "spearman_rho": rho,
            "n_valid_pairs": n_valid,
            "na_mismatch": na_mm,
            "r_na_count": r_na,
            "jl_na_count": jl_na,
            "mean_abs_diff": mad,
            "max_abs_diff": max_diff,
            "median_diff": median_diff,
            "r_mean": r_mean,
            "r_median": r_median,
            "r_sd": r_sd,
            "jl_mean": jl_mean,
            "jl_median": jl_median,
            "jl_sd": jl_sd,
        })

    res = pd.DataFrame(results)
    res["agreement_tier"] = res["r2_identity"].apply(classify_r2)
    res.to_csv(OUTPUT_CSV, index=False)
    print(f"\nPer-column CSV saved to: {OUTPUT_CSV}")

    # ── Generate markdown summary ──
    generate_summary_report(res, r, j, r_norm_resolved, jl_norm,
                            common_norms, r_only_norms, jl_only_norms,
                            len(common_gages), r_df.shape[0], jl_df.shape[0])

    # ── Console summary ──
    print_console_summary(res)

    return 0


def classify_r2(r2_val):
    """Classify a single R2 value into agreement tiers."""
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


# Ordered list of tiers for consistent display
TIERS = [
    ("Perfect",       "R2 >= 0.999"),
    ("Good",          "0.99 <= R2 < 0.999"),
    ("Poor",          "0.95 <= R2 < 0.99"),
    ("Low",           "0.90 <= R2 < 0.95"),
    ("Very Low",      "0.50 <= R2 < 0.90"),
    ("Extremely Low", "R2 < 0.50"),
]


def tier_counts(r2_series):
    """Return dict of tier_name -> count for a Series of R2 values."""
    counts = {}
    for name, _ in TIERS:
        counts[name] = 0
    for v in r2_series:
        t = classify_r2(v)
        if t in counts:
            counts[t] += 1
    return counts


def print_console_summary(res):
    """Print concise summary to console."""
    valid = res.dropna(subset=["r2_identity"])

    print(f"\n{'=' * 80}")
    print("OVERALL ALIGNMENT: Julia (post-Section 3) vs Golden R")
    print(f"{'=' * 80}")

    r2 = valid["r2_identity"]
    print(f"\n  Columns compared: {len(valid)}")
    print(f"  Mean R2:   {r2.mean():.6f}")
    print(f"  Median R2: {r2.median():.6f}")
    print(f"  SD R2:     {r2.std():.6f}")
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
            print(f"    {row['column']:<45} R2={r2v}  [{tier}]  MAD={row['mean_abs_diff']:.4f}  max_diff={row['max_abs_diff']:.4f}")


def generate_summary_report(res, r_df, jl_df, r_norm, jl_norm,
                            common_norms, r_only_norms, jl_only_norms,
                            n_common_gages, n_r_total, n_jl_total):
    """Write comprehensive markdown report."""

    valid = res.dropna(subset=["r2_identity"])
    r2 = valid["r2_identity"]

    with open(OUTPUT_MD, "w", encoding="utf-8") as f:
        f.write("# Julia (Post-Section 3) vs Golden R Reference: Comparison Report\n\n")
        f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # ── Input files ──
        f.write("## Input Files\n\n")
        f.write("| Dataset | File | Gages | Columns | Date |\n")
        f.write("|---------|------|-------|---------|------|\n")
        r_ts = datetime.fromtimestamp(os.path.getmtime(GOLDEN_R_PATH)).strftime("%Y-%m-%d") if GOLDEN_R_PATH.exists() else "N/A"
        jl_ts = datetime.fromtimestamp(os.path.getmtime(JULIA_PATH)).strftime("%Y-%m-%d") if JULIA_PATH.exists() else "N/A"
        f.write(f"| Golden R | `golden-outputs/streamflow_signatures_full_10feb2026.csv` | {n_r_total:,} | {len(r_df.columns)+1} | {r_ts} |\n")
        f.write(f"| Julia (post-Section 3) | `docs/benchmarks/julia_signatures.csv` | {n_jl_total:,} | {len(jl_df.columns)+1} | {jl_ts} |\n")
        f.write(f"\n**Common gages**: {n_common_gages:,}\n")
        f.write(f"**Common signature columns**: {len(common_norms):,}\n\n")

        # ── Context ──
        f.write("## Context\n\n")
        f.write("The Golden R output (Feb 2026) predates several code changes:\n")
        f.write("1. **Section 3 Guidelines changes** (Apr 2026, Julia only): D1/D99 timing, recession algorithm fix, n_recession_events, elasticity rename/annual, runoff ratio flag, elasticity diagnostics\n")
        f.write("2. **Constant-SD flag-only** (Apr 2026): Years with constant monthly SD are flagged but no longer rejected -- adds ~64 extra qualifying gages\n")
        f.write("3. **Negative-Q conditional rejection** (Apr 2026): Q<0 rejection now config-driven (default: off)\n")
        f.write("4. **R canonical alignment fixes** (Apr 2026): Removed leftover min_Q filter, FDC negative Q guard, seasonal_flags passthrough\n")
        f.write("5. **Pre-existing cross-language divergences**: Recession OLS library diffs, FDC90th floating-point precision, Mann-Kendall implementation diffs\n\n")
        f.write("Differences reflect ALL of the above, not just Section 3 changes.\n\n")

        # ── High-level summary ──
        f.write("## High-Level Alignment Summary\n\n")
        total = len(valid)
        tc = tier_counts(r2)

        f.write("### Distribution Statistics\n\n")
        f.write("| Metric | Value |\n")
        f.write("|--------|-------|\n")
        f.write(f"| Columns compared | {total} |\n")
        f.write(f"| Mean R2 (identity) | {r2.mean():.6f} |\n")
        f.write(f"| Median R2 | {r2.median():.6f} |\n")
        f.write(f"| SD of R2 | {r2.std():.6f} |\n")
        f.write(f"| Min R2 | {r2.min():.4f} |\n")
        f.write("\n")

        f.write("### Agreement Tiers\n\n")
        f.write("| Tier | Threshold | Count | % |\n")
        f.write("|------|-----------|-------|---|\n")
        for name, desc in TIERS:
            c = tc[name]
            f.write(f"| {name} | {desc} | {c} | {100*c/total:.1f}% |\n")
        f.write(f"| **Total** | | **{total}** | **100%** |\n")
        f.write("\n")

        # ── Agreement by category ──
        f.write("## Agreement by Signature Category\n\n")
        tier_hdrs = " | ".join(t[0] for t in TIERS)
        f.write(f"| Category | Cols | {tier_hdrs} | Mean R2 | Median R2 | SD R2 | Min R2 |\n")
        f.write(f"|----------|------|-" + "-|-".join("-" * len(t[0]) for t in TIERS) + "-|---------|-----------|-------|--------|\n")
        for cat in sorted(res["category"].unique()):
            cat_df = res[res["category"] == cat].dropna(subset=["r2_identity"])
            n = len(cat_df)
            if n == 0:
                continue
            r2c = cat_df["r2_identity"]
            tc_cat = tier_counts(r2c)
            tier_vals = " | ".join(f"{tc_cat[t[0]]}" for t in TIERS)
            f.write(f"| {cat} | {n} | {tier_vals} | "
                    f"{r2c.mean():.6f} | {r2c.median():.6f} | {r2c.std():.4f} | {r2c.min():.4f} |\n")
        f.write("\n")

        # ── Agreement by stat type ──
        f.write("## Agreement by Statistic Type\n\n")
        tier_hdrs = " | ".join(t[0] for t in TIERS)
        f.write(f"| Stat Type | Cols | {tier_hdrs} | Mean R2 | Min R2 |\n")
        f.write(f"|-----------|------|-" + "-|-".join("-" * len(t[0]) for t in TIERS) + "-|---------|--------|\n")
        for stat in ["mean", "median", "senn_slp", "linear_slp", "spearman_rho",
                      "spearman_pval", "mk_rho", "mk_pval", "scalar"]:
            stat_df = valid[valid["stat_type"] == stat] if stat in valid["stat_type"].values else pd.DataFrame()
            if len(stat_df) == 0:
                continue
            r2s = stat_df["r2_identity"]
            tc_stat = tier_counts(r2s)
            tier_vals = " | ".join(f"{tc_stat[t[0]]}" for t in TIERS)
            f.write(f"| {stat} | {len(stat_df)} | {tier_vals} | "
                    f"{r2s.mean():.6f} | {r2s.min():.4f} |\n")
        f.write("\n")

        # ── Columns with poor alignment ──
        poor_df = res[res["r2_identity"] < 0.99].sort_values("r2_identity")
        f.write(f"## Columns with R2 < 0.99 ({len(poor_df)} columns)\n\n")
        if len(poor_df) > 0:
            f.write("| Column | Category | Stat | R2 | Spearman | MAD | Max Diff | R NAs | Jl NAs | NA Mismatch |\n")
            f.write("|--------|----------|------|----|----------|-----|----------|-------|--------|-------------|\n")
            for _, row in poor_df.iterrows():
                r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
                rho = f"{row['spearman_rho']:.4f}" if np.isfinite(row['spearman_rho']) else "N/A"
                mad = f"{row['mean_abs_diff']:.4f}" if np.isfinite(row['mean_abs_diff']) else "N/A"
                maxd = f"{row['max_abs_diff']:.4f}" if np.isfinite(row['max_abs_diff']) else "N/A"
                f.write(f"| `{row['column']}` | {row['category']} | {row['stat_type']} | "
                        f"{r2v} | {rho} | {mad} | {maxd} | "
                        f"{row['r_na_count']} | {row['jl_na_count']} | {row['na_mismatch']} |\n")
        else:
            f.write("All columns have R2 >= 0.99.\n")
        f.write("\n")

        # ── Substantial changes analysis ──
        f.write("## Substantially Changed Signatures\n\n")
        f.write("Columns where mean or median values shifted meaningfully (outside typical cross-language noise).\n\n")

        # Filter to mean/median stat types and look at relative change
        central_stats = res[res["stat_type"].isin(["mean", "median"])].copy()
        central_stats["rel_change"] = np.where(
            np.abs(central_stats["r_mean"]) > 1e-6,
            np.abs(central_stats["jl_mean"] - central_stats["r_mean"]) / np.abs(central_stats["r_mean"]),
            np.abs(central_stats["jl_mean"] - central_stats["r_mean"])
        )

        # Substantial = R2 < 0.99 AND rel_change > 5% for means/medians
        substantial = central_stats[
            (central_stats["r2_identity"] < 0.99) |
            (central_stats["rel_change"] > 0.05)
        ].sort_values("r2_identity")

        if len(substantial) > 0:
            f.write("| Column | Category | R2 | R Mean | Julia Mean | Rel Change | Max Diff |\n")
            f.write("|--------|----------|----|--------|------------|------------|----------|\n")
            for _, row in substantial.head(40).iterrows():
                r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
                rc = f"{row['rel_change']:.1%}" if np.isfinite(row['rel_change']) else "N/A"
                r_m = f"{row['r_mean']:.4f}" if np.isfinite(row['r_mean']) else "N/A"
                j_m = f"{row['jl_mean']:.4f}" if np.isfinite(row['jl_mean']) else "N/A"
                maxd = f"{row['max_abs_diff']:.4f}" if np.isfinite(row['max_abs_diff']) else "N/A"
                f.write(f"| `{row['column']}` | {row['category']} | {r2v} | {r_m} | {j_m} | {rc} | {maxd} |\n")
            if len(substantial) > 40:
                f.write(f"\n*... and {len(substantial) - 40} more*\n")
        else:
            f.write("No substantially changed columns found.\n")
        f.write("\n")

        # ── NA mismatch analysis ──
        f.write("## NA Mismatch Analysis\n\n")
        f.write("Columns where the number of NAs differs between R and Julia by >50 gages.\n\n")
        na_issues = res[res["na_mismatch"] > 50].sort_values("na_mismatch", ascending=False)
        if len(na_issues) > 0:
            f.write("| Column | Category | R NAs | Julia NAs | Mismatch | R2 |\n")
            f.write("|--------|----------|-------|-----------|----------|----|\n")
            for _, row in na_issues.head(30).iterrows():
                r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
                f.write(f"| `{row['column']}` | {row['category']} | {row['r_na_count']} | {row['jl_na_count']} | {row['na_mismatch']} | {r2v} |\n")
            if len(na_issues) > 30:
                f.write(f"\n*... and {len(na_issues) - 30} more*\n")
        else:
            f.write("No columns with >50 NA mismatches.\n")
        f.write("\n")

        # ── New Julia-only columns ──
        f.write("## New Julia-Only Columns (not in Golden R)\n\n")
        f.write(f"**{len(jl_only_norms)} columns** present in Julia but not in Golden R reference:\n\n")
        if jl_only_norms:
            # Group by category
            jl_new_by_cat = {}
            for nc in jl_only_norms:
                orig = jl_norm[nc]
                base = get_base_metric(nc)
                cat = categorize_metric(base)
                jl_new_by_cat.setdefault(cat, []).append(orig)
            for cat in sorted(jl_new_by_cat):
                cols = jl_new_by_cat[cat]
                f.write(f"### {cat} ({len(cols)} columns)\n")
                for c in sorted(cols):
                    f.write(f"- `{c}`\n")
                f.write("\n")
        f.write("\n")

        # ── Deep dive on worst ──
        worst = poor_df.head(10)
        if len(worst) > 0:
            f.write(f"## Deep Dive: Worst {len(worst)} Columns\n\n")
            for _, row in worst.iterrows():
                col = row["column"]
                r_col = row["r_col"]
                jl_col = row["jl_col"]

                r_vals = pd.to_numeric(r_df[r_col], errors="coerce")
                jl_vals = pd.to_numeric(jl_df[jl_col], errors="coerce")

                r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
                rho = f"{row['spearman_rho']:.4f}" if np.isfinite(row['spearman_rho']) else "N/A"

                f.write(f"### `{col}` ({row['category']})\n\n")
                f.write(f"- **Identity R2**: {r2v}\n")
                f.write(f"- **Spearman rho**: {rho}\n\n")
                f.write("| Stat | Golden R | Julia |\n")
                f.write("|------|----------|-------|\n")
                f.write(f"| Count | {r_vals.notna().sum()} | {jl_vals.notna().sum()} |\n")
                f.write(f"| Mean | {r_vals.mean():.6f} | {jl_vals.mean():.6f} |\n")
                f.write(f"| Median | {r_vals.median():.6f} | {jl_vals.median():.6f} |\n")
                f.write(f"| SD | {r_vals.std():.6f} | {jl_vals.std():.6f} |\n")
                f.write(f"| Min | {r_vals.min():.6f} | {jl_vals.min():.6f} |\n")
                f.write(f"| Max | {r_vals.max():.6f} | {jl_vals.max():.6f} |\n")
                f.write(f"| NAs | {r_vals.isna().sum()} | {jl_vals.isna().sum()} |\n")
                f.write("\n")

        # ── Final summary box ──
        f.write("## Summary\n\n")
        f.write("| Agreement Tier | Threshold | Count | % |\n")
        f.write("|----------------|-----------|-------|---|\n")
        for name, desc in TIERS:
            c = tc[name]
            f.write(f"| {name} | {desc} | {c} | {100*c/total:.1f}% |\n")
        f.write(f"| **Total compared** | | **{total}** | **100%** |\n")
        f.write(f"\nNew Julia-only columns (Section 3 additions): **{len(jl_only_norms)}**\n")
        f.write("\n---\n")
        f.write("*Generated by `docs/benchmarks/compare_julia_vs_golden_r.py`*\n")

    print(f"\nMarkdown report saved to: {OUTPUT_MD}")


if __name__ == "__main__":
    sys.exit(main())

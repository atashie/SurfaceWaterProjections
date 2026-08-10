"""
Compare experiment output against Julia baseline (julia_signatures.csv).

Produces:
  - Console summary
  - docs/benchmarks/{name}_vs_julia_comparison.csv   (per-column detail)
  - docs/benchmarks/{name}_vs_julia_summary.md        (high-level report)

Usage:
    python docs/benchmarks/compare_experiment_vs_julia.py startIn1993
    python docs/benchmarks/compare_experiment_vs_julia.py startIn1993_60pct
    python docs/benchmarks/compare_experiment_vs_julia.py startIn1993 --baseline path/to/baseline.csv --experiment path/to/experiment.csv
"""

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent

# ── Helpers ──────────────────────────────────────────────────────────

# SCOPE NOTE (2026-07-22, Codex finding): the four season_excluded_years_* per-gage
# diagnostics below are treated as METADATA here (excluded from the compared-column
# summary and tiers), but build_experiment_vs_julia_dashboard.py DOES visualize them
# as single-value targets — so the dashboard's target count exceeds this script's
# compared-column count by 4. Deliberate: they are window-dependent diagnostic
# counters, not signatures (an uncapped-vs-capped window comparison shifts them
# uniformly by ±1 per partial trailing year), and would distort signature tiers.
META_COLS = {
    "gage_id", "gage_id_metadata", "latitude", "longitude", "basin_area",
    "basin_area_km2", "gage_type", "num_water_years", "start_year", "end_year",
    "start_water_year", "end_water_year", "country",
    "processing_status", "drainage_area_km2", "area_normalized",
    "human_interference_class", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
    "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX",
    "CLASS", "RHBN", "REGULATED",
    "ice_affected_days_total",
    "season_excluded_years_winter", "season_excluded_years_spring",
    "season_excluded_years_summer", "season_excluded_years_fall",
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
    # Drought and Snow are matched first: both families are newer than the rules
    # below and would otherwise fall through to "Other" (e.g. b_* would claim
    # nothing here, but a future rename could collide).
    if base.startswith("drought_"):
        return "Streamflow Drought"
    if base.startswith(("swe_", "snow_", "melt_")) or base == "ssm":
        return "Snow"
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


def get_experiment_description(name):
    """Auto-detect experiment type and return description."""
    name_lower = name.lower()
    if "1993" in name_lower and "80pct" in name_lower:
        return ("Water years >= 1993 AND at least 80% of possible years "
                "(WY1993 to gage max) must have qualifying data.")
    if "1993" in name_lower and "60pct" in name_lower:
        return ("Water years >= 1993 AND at least 60% of possible years "
                "(WY1993 to gage max) must have qualifying data.")
    if "1993" in name_lower:
        return "Water years >= 1993 only (restricts analysis period from ~1980)."
    return f"Custom experiment: {name}"


# ── Main comparison ───────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compare experiment vs Julia baseline")
    parser.add_argument("experiment_name", help="Name of the experiment (e.g., startIn1993)")
    parser.add_argument("--baseline", default=None, help="Path to baseline CSV")
    parser.add_argument("--experiment", default=None, help="Path to experiment CSV")
    parser.add_argument("--output-dir", default=None,
                        help="Directory for the comparison CSV + summary MD. Default: the "
                             "experiment CSV's folder when --experiment is given, else "
                             "docs/benchmarks (legacy). Convention (July 2026): ALL "
                             "experiment artifacts live in the experiment's own folder.")
    args = parser.parse_args()

    name = args.experiment_name
    baseline_path = Path(args.baseline) if args.baseline else SCRIPT_DIR / "julia_signatures.csv"
    experiment_path = Path(args.experiment) if args.experiment else SCRIPT_DIR / f"{name}_signatures.csv"
    out_dir = Path(args.output_dir) if args.output_dir else (
        experiment_path.parent if args.experiment else SCRIPT_DIR)
    output_csv = out_dir / f"{name}_vs_julia_comparison.csv"
    output_md = out_dir / f"{name}_vs_julia_summary.md"

    print("=" * 80)
    print(f"EXPERIMENT '{name}' vs JULIA BASELINE")
    print("=" * 80)
    print(f"  Filter: {get_experiment_description(name)}")

    # Load data
    print("\nLoading Julia baseline...")
    if not baseline_path.exists():
        print(f"ERROR: Baseline not found: {baseline_path}")
        return 1
    base_df = pd.read_csv(baseline_path, low_memory=False)
    base_df["gage_id"] = base_df["gage_id"].astype(str).str.strip()
    print(f"  {base_df.shape[0]:,} gages x {base_df.shape[1]:,} cols")

    print(f"Loading experiment '{name}'...")
    if not experiment_path.exists():
        print(f"ERROR: Experiment not found: {experiment_path}")
        return 1
    exp_df = pd.read_csv(experiment_path, low_memory=False)
    exp_df["gage_id"] = exp_df["gage_id"].astype(str).str.strip()
    print(f"  {exp_df.shape[0]:,} gages x {exp_df.shape[1]:,} cols")

    # Gage analysis
    base_gages = set(base_df["gage_id"])
    exp_gages = set(exp_df["gage_id"])
    common_gages = sorted(base_gages & exp_gages)
    baseline_only = sorted(base_gages - exp_gages)
    experiment_only = sorted(exp_gages - base_gages)

    print(f"\n  Common gages: {len(common_gages):,}")
    print(f"  Baseline only (dropped): {len(baseline_only):,}")
    print(f"  Experiment only (new): {len(experiment_only):,}")

    # Gage diff breakdown by type and interference class
    gage_diff_report = _gage_diff_analysis(base_df, exp_df, baseline_only, experiment_only)

    # Subset to common gages
    b = base_df[base_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()
    e = exp_df[exp_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()

    # Years-per-gage comparison
    years_report = ""
    if "num_water_years" in b.columns and "num_water_years" in e.columns:
        b_years = pd.to_numeric(b["num_water_years"], errors="coerce")
        e_years = pd.to_numeric(e["num_water_years"], errors="coerce")
        years_report = (
            f"\n  Years per gage (common gages):\n"
            f"    Baseline:   mean={b_years.mean():.1f}, median={b_years.median():.1f}\n"
            f"    Experiment: mean={e_years.mean():.1f}, median={e_years.median():.1f}\n"
            f"    Mean diff:  {e_years.mean() - b_years.mean():.1f} years"
        )
        print(years_report)

    # Column matching (same code = same columns, no renames needed)
    base_sig_cols = {c for c in b.columns if is_signature_col(c)}
    exp_sig_cols = {c for c in e.columns if is_signature_col(c)}
    common_cols = sorted(base_sig_cols & exp_sig_cols)

    print(f"\n  Common signature columns: {len(common_cols):,}")
    if base_sig_cols - exp_sig_cols:
        print(f"  Baseline-only columns: {len(base_sig_cols - exp_sig_cols)}")
    if exp_sig_cols - base_sig_cols:
        print(f"  Experiment-only columns: {len(exp_sig_cols - base_sig_cols)}")

    if not common_gages:
        print("\nERROR: No common gages between baseline and experiment. Nothing to compare.")
        return 1
    if not common_cols:
        print("\nERROR: No common signature columns between baseline and experiment. Nothing to compare.")
        return 1

    # ── Per-column comparison ──
    print(f"\n{'=' * 80}")
    print("COMPUTING PER-COLUMN METRICS...")
    print(f"{'=' * 80}")

    results = []
    for col in common_cols:
        base_vals = pd.to_numeric(b[col], errors="coerce").values
        exp_vals = pd.to_numeric(e[col], errors="coerce").values

        r2, n_valid = r2_identity(base_vals, exp_vals)
        rho, _ = spearman_corr(base_vals, exp_vals)
        na_mm = count_na_mismatch(base_vals, exp_vals)

        base_na = int((~np.isfinite(base_vals)).sum())
        exp_na = int((~np.isfinite(exp_vals)).sum())

        mask = np.isfinite(base_vals) & np.isfinite(exp_vals)
        if mask.sum() > 0:
            diff = exp_vals[mask] - base_vals[mask]
            abs_diff = np.abs(diff)
            mad = np.mean(abs_diff)
            max_diff = np.max(abs_diff)
            median_diff = np.median(diff)
            base_mean = np.mean(base_vals[mask])
            exp_mean = np.mean(exp_vals[mask])
        else:
            mad = max_diff = median_diff = np.nan
            base_mean = exp_mean = np.nan

        base_metric = get_base_metric(col)
        stat = get_stat_suffix(col)
        cat = categorize_metric(base_metric)

        results.append({
            "column": col,
            "category": cat,
            "base_metric": base_metric,
            "stat_type": stat,
            "r2_identity": r2,
            "spearman_rho": rho,
            "n_valid_pairs": n_valid,
            "na_mismatch": na_mm,
            "baseline_na_count": base_na,
            "experiment_na_count": exp_na,
            "mean_abs_diff": mad,
            "max_abs_diff": max_diff,
            "median_diff": median_diff,
            "baseline_mean": base_mean,
            "experiment_mean": exp_mean,
        })

    res = pd.DataFrame(results)
    res["agreement_tier"] = res["r2_identity"].apply(classify_r2)
    res.to_csv(output_csv, index=False)
    print(f"\nPer-column CSV saved to: {output_csv}")

    # ── Generate markdown summary ──
    _generate_summary_report(
        res, name, output_md,
        baseline_path, experiment_path,
        base_df.shape[0], exp_df.shape[0],
        len(common_gages), len(baseline_only), len(experiment_only),
        gage_diff_report, years_report,
        b, e
    )

    # ── Console summary ──
    _print_console_summary(res, name)

    return 0


def _gage_diff_analysis(base_df, exp_df, baseline_only, experiment_only):
    """Analyze dropped/added gages by type and interference class."""
    lines = []

    if not baseline_only and not experiment_only:
        return "No gage differences.\n"

    for label, gage_list, source_df in [
        ("Dropped (baseline only)", baseline_only, base_df),
        ("Added (experiment only)", experiment_only, exp_df),
    ]:
        if not gage_list:
            continue
        subset = source_df[source_df["gage_id"].isin(gage_list)]
        lines.append(f"\n  {label}: {len(gage_list)} gages")

        if "gage_type" in subset.columns:
            type_counts = subset["gage_type"].value_counts()
            for gt, count in type_counts.items():
                lines.append(f"    {gt}: {count}")

        if "human_interference_class" in subset.columns:
            hic_counts = subset["human_interference_class"].value_counts()
            for hic, count in hic_counts.items():
                lines.append(f"    {hic}: {count}")

        # Show a few representative gages
        sample = gage_list[:5]
        lines.append(f"    Examples: {', '.join(sample)}")

    report = "\n".join(lines)
    print(report)
    return report


def _print_console_summary(res, name):
    valid = res.dropna(subset=["r2_identity"])

    print(f"\n{'=' * 80}")
    print(f"OVERALL ALIGNMENT: '{name}' vs Julia Baseline")
    print(f"{'=' * 80}")

    r2 = valid["r2_identity"]
    print(f"\n  Columns compared: {len(valid)}")
    print(f"  Mean R2:   {r2.mean():.6f}")
    print(f"  Median R2: {r2.median():.6f}")
    print(f"  Min R2:    {r2.min():.6f}")

    tc = tier_counts(r2)
    total = len(valid)
    print()
    for t_name, desc in TIERS:
        c = tc[t_name]
        print(f"  {t_name:<16} ({desc:<24}): {c:>4} ({100*c/total:.1f}%)")

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


def _generate_summary_report(res, name, output_md,
                              baseline_path, experiment_path,
                              n_base_total, n_exp_total,
                              n_common, n_dropped, n_added,
                              gage_diff_report, years_report,
                              b_df, e_df):
    valid = res.dropna(subset=["r2_identity"])
    r2 = valid["r2_identity"]

    with open(output_md, "w", encoding="utf-8") as f:
        f.write(f"# Experiment '{name}' vs Julia Baseline: Comparison Report\n\n")
        f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Context
        f.write("## Experiment Description\n\n")
        f.write(f"{get_experiment_description(name)}\n\n")

        # Input files
        f.write("## Input Files\n\n")
        f.write("| Dataset | File | Gages | Columns |\n")
        f.write("|---------|------|-------|---------|\n")
        f.write(f"| Julia Baseline | `{baseline_path.name}` | {n_base_total:,} | {len(b_df.columns)+1} |\n")
        f.write(f"| Experiment ({name}) | `{experiment_path.name}` | {n_exp_total:,} | {len(e_df.columns)+1} |\n")
        f.write(f"\n**Common gages**: {n_common:,}\n")
        f.write(f"**Dropped gages** (in baseline, not in experiment): {n_dropped:,}\n")
        f.write(f"**Added gages** (in experiment, not in baseline): {n_added:,}\n\n")

        # Gage diff
        if n_dropped > 0 or n_added > 0:
            f.write("### Gage Differences\n\n")
            f.write("```\n")
            f.write(gage_diff_report)
            f.write("\n```\n\n")

        # Years comparison
        if years_report:
            f.write("### Years Per Gage (Common Gages)\n\n")
            f.write("```\n")
            f.write(years_report.strip())
            f.write("\n```\n\n")

        # High-level summary
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
        for t_name, desc in TIERS:
            c = tc[t_name]
            f.write(f"| {t_name} | {desc} | {c} | {100*c/total:.1f}% |\n")
        f.write(f"| **Total** | | **{total}** | **100%** |\n")
        f.write("\n")

        # By category
        f.write("## Agreement by Signature Category\n\n")
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
            f.write(f"| {cat} | {n} | {tier_vals} | "
                    f"{r2c.mean():.6f} | {r2c.min():.4f} |\n")
        f.write("\n")

        # By stat type
        f.write("## Agreement by Statistic Type\n\n")
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

        # Worst columns
        poor_df = res[res["r2_identity"] < 0.99].sort_values("r2_identity")
        f.write(f"## Columns with R2 < 0.99 ({len(poor_df)} columns)\n\n")
        if len(poor_df) > 0:
            f.write("| Column | Category | Stat | R2 | Spearman | MAD | Max Diff | Baseline NAs | Experiment NAs | NA Mismatch |\n")
            f.write("|--------|----------|------|----|----------|-----|----------|--------------|----------------|-------------|\n")
            for _, row in poor_df.iterrows():
                r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
                rho = f"{row['spearman_rho']:.4f}" if np.isfinite(row['spearman_rho']) else "N/A"
                mad = f"{row['mean_abs_diff']:.4f}" if np.isfinite(row['mean_abs_diff']) else "N/A"
                maxd = f"{row['max_abs_diff']:.4f}" if np.isfinite(row['max_abs_diff']) else "N/A"
                f.write(f"| `{row['column']}` | {row['category']} | {row['stat_type']} | "
                        f"{r2v} | {rho} | {mad} | {maxd} | "
                        f"{row['baseline_na_count']} | {row['experiment_na_count']} | {row['na_mismatch']} |\n")
        else:
            f.write("All columns have R2 >= 0.99.\n")
        f.write("\n")

        # NA mismatch
        f.write("## NA Mismatch Analysis\n\n")
        f.write("Columns where the number of NAs differs by >50 gages.\n\n")
        na_issues = res[res["na_mismatch"] > 50].sort_values("na_mismatch", ascending=False)
        if len(na_issues) > 0:
            f.write("| Column | Category | Baseline NAs | Experiment NAs | Mismatch | R2 |\n")
            f.write("|--------|----------|--------------|----------------|----------|----|\n")
            for _, row in na_issues.head(30).iterrows():
                r2v = f"{row['r2_identity']:.4f}" if np.isfinite(row['r2_identity']) else "N/A"
                f.write(f"| `{row['column']}` | {row['category']} | {row['baseline_na_count']} | "
                        f"{row['experiment_na_count']} | {row['na_mismatch']} | {r2v} |\n")
        else:
            f.write("No columns with >50 NA mismatches.\n")
        f.write("\n")

        # Summary
        f.write("## Summary\n\n")
        f.write("| Agreement Tier | Threshold | Count | % |\n")
        f.write("|----------------|-----------|-------|---|\n")
        for t_name, desc in TIERS:
            c = tc[t_name]
            f.write(f"| {t_name} | {desc} | {c} | {100*c/total:.1f}% |\n")
        f.write(f"| **Total compared** | | **{total}** | **100%** |\n")
        f.write(f"\nGages dropped by experiment filter: **{n_dropped}**\n")
        f.write("\n---\n")
        f.write(f"*Generated by `docs/benchmarks/compare_experiment_vs_julia.py {name}`*\n")

    print(f"\nMarkdown report saved to: {output_md}")


if __name__ == "__main__":
    sys.exit(main())

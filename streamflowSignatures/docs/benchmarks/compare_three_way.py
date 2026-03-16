"""
Three-Way Cross-Language Comparison: R (canonical) vs Python vs Julia

Loads all three benchmark CSVs, identifies common gages and columns,
computes Spearman correlations, flags divergent metrics and gages.
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import stats

BENCHMARK_DIR = Path(__file__).parent

# ── Load data ──────────────────────────────────────────────────────────

def load_all():
    r_df = pd.read_csv(BENCHMARK_DIR / "r_signatures.csv", low_memory=False)
    py_df = pd.read_csv(BENCHMARK_DIR / "python_signatures.csv", low_memory=False)
    jl_df = pd.read_csv(BENCHMARK_DIR / "julia_signatures.csv", low_memory=False)

    # Standardize gage_id to string
    for df in (r_df, py_df, jl_df):
        df["gage_id"] = df["gage_id"].astype(str).str.strip()

    return r_df, py_df, jl_df


# ── Helpers ────────────────────────────────────────────────────────────

# Metadata / non-signature columns to exclude from comparison
META_COLS = {
    "gage_id", "gage_id_metadata", "latitude", "longitude", "basin_area",
    "basin_area_km2", "gage_type", "num_water_years", "start_year", "end_year",
    "start_water_year", "end_water_year", "country",
    "processing_status", "drainage_area_km2", "area_normalized",
    "human_interference_class", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
    "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX",
    "CLASS", "RHBN", "REGULATED",
}

# QA/QC flag columns
QAQC_PREFIX = "flagged_for_"

def is_signature_col(col):
    """Is this column a numeric signature (not metadata or QA/QC flag)?"""
    if col in META_COLS:
        return False
    if col.startswith(QAQC_PREFIX):
        return False
    return True

def normalize_col(col):
    """Normalize column name for cross-language matching.
    Handles Q95-Q10 vs Q95_Q10 and FDC_all vs FDCall naming."""
    col = col.replace("-", "_")
    col = col.replace(".", "_")  # R uses dots (Q95.Q10), Python/Julia use underscores
    col = col.replace("FDC_all", "FDCall")
    col = col.replace("FDC_90th", "FDC90th")
    col = col.replace("FDC_mid", "FDCmid")
    return col

def spearman_corr(x, y):
    """Spearman correlation on paired non-NA values. Returns (rho, pval, n)."""
    mask = np.isfinite(x) & np.isfinite(y)
    n = mask.sum()
    if n < 5:
        return np.nan, np.nan, n
    rho, pval = stats.spearmanr(x[mask], y[mask])
    return rho, pval, n

def count_na_mismatch(x, y):
    """Count how many values are NA in one but not the other."""
    x_na = ~np.isfinite(x)
    y_na = ~np.isfinite(y)
    return int((x_na & ~y_na).sum() + (~x_na & y_na).sum())


# ── Signature categorization ──────────────────────────────────────────

STAT_SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                 "_mk_rho", "_mk_pval", "_mean", "_median"]

def get_base_metric(col):
    """Extract the base metric name from a column with stat suffix."""
    for suf in STAT_SUFFIXES:
        if col.endswith(suf):
            return col[:-len(suf)]
    return col

def categorize_metric(base):
    """Assign a category to a base metric name."""
    if base.startswith(("Qann", "Qwin", "Qspr", "Qsum", "Qfal")):
        return "Flow Volumes"
    if base.startswith("Q") and any(c.isdigit() for c in base[1:3]):
        return "Flow Percentiles"
    if base.startswith("FDC"):
        return "FDC"
    if base.startswith("BFI"):
        return "Baseflow"
    if base.startswith(("log_a", "b_", "concavity")):
        return "Recession"
    if base.startswith(("n_high", "n_low", "dur_high", "dur_low", "TQmean", "Flow_Reversal")):
        return "Pulse Metrics"
    if base.startswith("flashiness"):
        return "Flashiness"
    if base.startswith(("D5", "D10", "D20", "D25", "D30", "D40", "D50",
                        "D60", "D70", "D80", "D90", "D95", "Dmax")):
        return "Flow Timing"
    if "runoff_ratio" in base:
        return "Runoff Ratios"
    if "elasticity" in base:
        return "Elasticity"
    if base.startswith("qp_"):
        return "Q-P Seasonality"
    if base.startswith("avg_storage"):
        return "Storage"
    return "Other"


def _file_timestamp(path):
    """Return modification time of a file as a string, or 'N/A'."""
    try:
        mtime = os.path.getmtime(path)
        return datetime.fromtimestamp(mtime).strftime("%Y-%m-%d %H:%M:%S")
    except OSError:
        return "N/A"


def _load_timing(name):
    """Load timing JSON for a language, returning dict or empty."""
    path = BENCHMARK_DIR / f"{name}_timing.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def generate_three_way_report(res_df, r_df, py_df, jl_df,
                              common_norms, r_norm, py_norm, jl_norm,
                              pj_only, rp_only, rj_only):
    """Write a comprehensive markdown report to benchmarks/comparison_report.md."""
    report_path = BENCHMARK_DIR / "comparison_report.md"
    n_common_gages = len(r_df)

    # Load timing
    py_timing = _load_timing("python")
    jl_timing = _load_timing("julia")

    # Categorize all columns
    for _, row in res_df.iterrows():
        base = get_base_metric(row["column"])
        res_df.loc[_, "category"] = categorize_metric(base)

    with open(report_path, "w") as f:
        # ── Header ──
        f.write("# Three-Way Cross-Language Comparison Report\n\n")
        f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Input Files\n\n")
        f.write("| File | Timestamp | Gages | Columns |\n")
        f.write("|------|-----------|-------|---------|\n")
        for name, csv_name in [("R", "r_signatures.csv"),
                                ("Python", "python_signatures.csv"),
                                ("Julia", "julia_signatures.csv")]:
            p = BENCHMARK_DIR / csv_name
            ts = _file_timestamp(p)
            df_temp = pd.read_csv(p, nrows=0)
            n_cols = len(df_temp.columns)
            n_gages_line = {"R": r_df, "Python": py_df, "Julia": jl_df}
            # Read gage count from full file via timing or csv
            if name == "R":
                ng = pd.read_csv(p, usecols=["gage_id"]).shape[0]
            else:
                t = _load_timing(name.lower())
                ng = t.get("n_gages_processed", "?")
            f.write(f"| {name} (`{csv_name}`) | {ts} | {ng} | {n_cols} |\n")

        f.write(f"\n**Common gages (all 3)**: {n_common_gages:,}\n")
        f.write(f"**Common signature columns**: {len(common_norms):,}\n\n")

        # ── Timing ──
        f.write("## Performance Comparison\n\n")
        f.write("| Metric | Python | Julia | Ratio (Py/Jl) |\n")
        f.write("|--------|--------|-------|---------------|\n")
        py_total = py_timing.get("total_seconds", 0)
        jl_total = jl_timing.get("total_seconds", 0)
        ratio = py_total / jl_total if jl_total > 0 else float("inf")
        f.write(f"| Total time | {py_total:.0f}s ({py_total/60:.1f} min) | {jl_total:.0f}s ({jl_total/60:.1f} min) | {ratio:.1f}x |\n")
        py_ng = py_timing.get("n_gages_processed", 0)
        jl_ng = jl_timing.get("n_gages_processed", 0)
        f.write(f"| Gages processed | {py_ng:,} | {jl_ng:,} | - |\n")
        py_rate = py_ng / py_total if py_total > 0 else 0
        jl_rate = jl_ng / jl_total if jl_total > 0 else 0
        f.write(f"| Processing rate | {py_rate:.2f}/s | {jl_rate:.2f}/s | {jl_rate/py_rate:.1f}x faster |\n")
        f.write(f"| R benchmark | ~1-2 hours (estimated) | - | - |\n")
        f.write("\n")

        # Phase breakdown
        if "phases" in py_timing and "phases" in jl_timing:
            f.write("### Phase Breakdown\n\n")
            f.write("| Phase | Python | Julia |\n")
            f.write("|-------|--------|-------|\n")
            all_phases = sorted(set(py_timing["phases"]) | set(jl_timing["phases"]))
            for phase in all_phases:
                py_p = py_timing["phases"].get(phase, 0)
                jl_p = jl_timing["phases"].get(phase, 0)
                f.write(f"| {phase} | {py_p:.0f}s | {jl_p:.0f}s |\n")
            f.write("\n")

        # ── Overall Spearman summary ──
        f.write("## Overall Spearman Correlation Summary\n\n")
        f.write("| Pair | Mean rho | Median rho | Min rho | Cols < 0.99 | Cols < 0.95 | Cols < 0.90 |\n")
        f.write("|------|----------|------------|---------|-------------|-------------|-------------|\n")
        valid = res_df.dropna(subset=["rp_rho", "rj_rho", "pj_rho"])
        for pair, label in [("rp_rho", "R vs Python"), ("rj_rho", "R vs Julia"), ("pj_rho", "Python vs Julia")]:
            vals = valid[pair]
            f.write(f"| {label} | {vals.mean():.6f} | {vals.median():.6f} | {vals.min():.4f} | "
                    f"{(vals < 0.99).sum()} | {(vals < 0.95).sum()} | {(vals < 0.90).sum()} |\n")
        f.write("\n")

        # ── Agreement by category ──
        f.write("## Agreement by Signature Category\n\n")
        f.write("| Category | Total Cols | Perfect (>=0.999) | Good (>=0.99) | Poor (<0.99) | Min rho |\n")
        f.write("|----------|-----------|-------------------|---------------|-------------|--------|\n")
        categories = sorted(res_df["category"].unique())
        for cat in categories:
            cat_df = res_df[res_df["category"] == cat]
            n_total = len(cat_df)
            n_perfect = (cat_df["min_rho"] >= 0.999).sum()
            n_good = ((cat_df["min_rho"] >= 0.99) & (cat_df["min_rho"] < 0.999)).sum()
            n_poor = (cat_df["min_rho"] < 0.99).sum()
            min_rho = cat_df["min_rho"].min()
            min_str = f"{min_rho:.4f}" if np.isfinite(min_rho) else "N/A"
            f.write(f"| {cat} | {n_total} | {n_perfect} | {n_good} | {n_poor} | {min_str} |\n")
        f.write("\n")

        # ── Columns with min rho < 0.99 ──
        poor = res_df[res_df["min_rho"] < 0.99].sort_values("min_rho")
        f.write(f"## Columns with min Spearman < 0.99 ({len(poor)} columns)\n\n")
        if len(poor) > 0:
            f.write("| Column | Category | R-Py | R-Jl | Py-Jl | R NA | Py NA | Jl NA |\n")
            f.write("|--------|----------|------|------|-------|------|-------|-------|\n")
            for _, row in poor.iterrows():
                rp = f"{row['rp_rho']:.4f}" if np.isfinite(row['rp_rho']) else "N/A"
                rj = f"{row['rj_rho']:.4f}" if np.isfinite(row['rj_rho']) else "N/A"
                pj = f"{row['pj_rho']:.4f}" if np.isfinite(row['pj_rho']) else "N/A"
                cat = row.get("category", "")
                f.write(f"| `{row['column']}` | {cat} | {rp} | {rj} | {pj} | "
                        f"{row['r_na_count']} | {row['py_na_count']} | {row['jl_na_count']} |\n")
        else:
            f.write("All columns have rho >= 0.99 across all pairs.\n")
        f.write("\n")

        # ── NA mismatch analysis ──
        f.write("## NA Mismatch Analysis\n\n")
        na_issues = res_df[(res_df["na_mismatch_rp"] > 100) |
                           (res_df["na_mismatch_rj"] > 100) |
                           (res_df["na_mismatch_pj"] > 100)]
        na_issues = na_issues.sort_values("na_mismatch_rp", ascending=False)
        f.write(f"Columns with >100 NA mismatches in any pair: **{len(na_issues)}**\n\n")
        if len(na_issues) > 0:
            f.write("| Column | NA R-Py | NA R-Jl | NA Py-Jl | R NA | Py NA | Jl NA |\n")
            f.write("|--------|---------|---------|----------|------|-------|-------|\n")
            for _, row in na_issues.head(50).iterrows():
                f.write(f"| `{row['column']}` | {row['na_mismatch_rp']} | {row['na_mismatch_rj']} | "
                        f"{row['na_mismatch_pj']} | {row['r_na_count']} | {row['py_na_count']} | "
                        f"{row['jl_na_count']} |\n")
            if len(na_issues) > 50:
                f.write(f"\n*... and {len(na_issues) - 50} more columns*\n")
        else:
            f.write("No columns with >100 NA mismatches.\n")
        f.write("\n")

        # ── Columns only in subset of languages ──
        f.write("## Column Coverage\n\n")
        r_only = sorted(set(r_norm) - set(py_norm) - set(jl_norm))
        py_only_cols = sorted(set(py_norm) - set(r_norm) - set(jl_norm))
        jl_only_cols = sorted(set(jl_norm) - set(r_norm) - set(py_norm))

        f.write(f"| Scope | Count |\n")
        f.write(f"|-------|-------|\n")
        f.write(f"| All 3 languages | {len(common_norms)} |\n")
        f.write(f"| R + Python only | {len(rp_only)} |\n")
        f.write(f"| R + Julia only | {len(rj_only)} |\n")
        f.write(f"| Python + Julia only | {len(pj_only)} |\n")
        f.write(f"| Only R | {len(r_only)} |\n")
        f.write(f"| Only Python | {len(py_only_cols)} |\n")
        f.write(f"| Only Julia | {len(jl_only_cols)} |\n")
        f.write("\n")

        if r_only:
            f.write(f"### Columns only in R ({len(r_only)})\n\n")
            for c in r_only:
                f.write(f"- `{r_norm[c]}`\n")
            f.write("\n")
        if py_only_cols:
            f.write(f"### Columns only in Python ({len(py_only_cols)})\n\n")
            for c in py_only_cols:
                f.write(f"- `{py_norm[c]}`\n")
            f.write("\n")
        if jl_only_cols:
            f.write(f"### Columns only in Julia ({len(jl_only_cols)})\n\n")
            for c in jl_only_cols:
                f.write(f"- `{jl_norm[c]}`\n")
            f.write("\n")

        # ── Deep dive on worst metrics ──
        worst = poor.head(10)
        if len(worst) > 0:
            f.write(f"## Deep Dive: Worst {len(worst)} Metrics\n\n")
            for _, row in worst.iterrows():
                col = row["column"]
                r_col = row["r_col"]
                py_col = row["py_col"]
                jl_col = row["jl_col"]

                r_vals = pd.to_numeric(r_df[r_col], errors="coerce")
                py_vals = pd.to_numeric(py_df[py_col], errors="coerce")
                jl_vals = pd.to_numeric(jl_df[jl_col], errors="coerce")

                rp = f"{row['rp_rho']:.4f}" if np.isfinite(row['rp_rho']) else "N/A"
                rj = f"{row['rj_rho']:.4f}" if np.isfinite(row['rj_rho']) else "N/A"
                pj = f"{row['pj_rho']:.4f}" if np.isfinite(row['pj_rho']) else "N/A"

                f.write(f"### `{col}`\n\n")
                f.write(f"**Spearman**: R-Py={rp}, R-Jl={rj}, Py-Jl={pj}\n\n")
                f.write("| Stat | R | Python | Julia |\n")
                f.write("|------|---|--------|-------|\n")
                f.write(f"| Min | {r_vals.min():.4f} | {py_vals.min():.4f} | {jl_vals.min():.4f} |\n")
                f.write(f"| Median | {r_vals.median():.4f} | {py_vals.median():.4f} | {jl_vals.median():.4f} |\n")
                f.write(f"| Max | {r_vals.max():.4f} | {py_vals.max():.4f} | {jl_vals.max():.4f} |\n")
                f.write(f"| NAs | {r_vals.isna().sum()} | {py_vals.isna().sum()} | {jl_vals.isna().sum()} |\n")
                f.write("\n")

        # ── Agreement summary ──
        perfect = res_df[res_df["min_rho"] >= 0.999]
        good = res_df[(res_df["min_rho"] >= 0.99) & (res_df["min_rho"] < 0.999)]
        f.write("## Summary\n\n")
        f.write(f"| Agreement Level | Count | % |\n")
        f.write(f"|-----------------|-------|---|\n")
        total = len(res_df)
        f.write(f"| Perfect (rho >= 0.999) | {len(perfect)} | {100*len(perfect)/total:.1f}% |\n")
        f.write(f"| Good (0.99 <= rho < 0.999) | {len(good)} | {100*len(good)/total:.1f}% |\n")
        f.write(f"| Poor (rho < 0.99) | {len(poor)} | {100*len(poor)/total:.1f}% |\n")
        f.write(f"| **Total compared** | **{total}** | **100%** |\n")
        f.write("\n---\n")
        f.write("*Generated by `docs/benchmarks/compare_three_way.py`*\n")

    print(f"\nMarkdown report saved to: {report_path}")


# ── Main comparison ───────────────────────────────────────────────────

def main():
    print("=" * 80)
    print("THREE-WAY CROSS-LANGUAGE COMPARISON: R vs Python vs Julia")
    print("=" * 80)

    # Load
    print("\nLoading CSVs...")
    r_df, py_df, jl_df = load_all()
    print(f"  R:      {r_df.shape[0]:,} gages, {r_df.shape[1]:,} columns")
    print(f"  Python: {py_df.shape[0]:,} gages, {py_df.shape[1]:,} columns")
    print(f"  Julia:  {jl_df.shape[0]:,} gages, {jl_df.shape[1]:,} columns")

    # ── Common gages ──
    r_gages = set(r_df["gage_id"])
    py_gages = set(py_df["gage_id"])
    jl_gages = set(jl_df["gage_id"])
    common_gages = r_gages & py_gages & jl_gages
    print(f"\n  Common gages (all 3): {len(common_gages):,}")
    print(f"  Only in R: {len(r_gages - py_gages - jl_gages):,}")
    print(f"  Only in Python/Julia (not R): {len(py_gages - r_gages):,}")

    # Filter to common gages
    r_df = r_df[r_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()
    py_df = py_df[py_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()
    jl_df = jl_df[jl_df["gage_id"].isin(common_gages)].set_index("gage_id").sort_index()

    # ── Column mapping ──
    r_sig_cols = {c for c in r_df.columns if is_signature_col(c)}
    py_sig_cols = {c for c in py_df.columns if is_signature_col(c)}
    jl_sig_cols = {c for c in jl_df.columns if is_signature_col(c)}

    # Build normalized name → original name mapping
    r_norm = {normalize_col(c): c for c in r_sig_cols}
    py_norm = {normalize_col(c): c for c in py_sig_cols}
    jl_norm = {normalize_col(c): c for c in jl_sig_cols}

    all_norms = set(r_norm) | set(py_norm) | set(jl_norm)
    common_norms = set(r_norm) & set(py_norm) & set(jl_norm)
    rp_only = set(r_norm) & set(py_norm) - set(jl_norm)
    rj_only = set(r_norm) & set(jl_norm) - set(py_norm)
    pj_only = set(py_norm) & set(jl_norm) - set(r_norm)

    print(f"\n  Signature columns across all 3: {len(common_norms):,}")
    print(f"  In R + Python but not Julia: {len(rp_only)}")
    print(f"  In R + Julia but not Python: {len(rj_only)}")
    print(f"  In Python + Julia but not R: {len(pj_only)}")
    print(f"  Only in R: {len(set(r_norm) - set(py_norm) - set(jl_norm))}")
    print(f"  Only in Python: {len(set(py_norm) - set(r_norm) - set(jl_norm))}")
    print(f"  Only in Julia: {len(set(jl_norm) - set(r_norm) - set(py_norm))}")

    # Print columns only in one
    r_only = sorted(set(r_norm) - set(py_norm) - set(jl_norm))
    py_only = sorted(set(py_norm) - set(r_norm) - set(jl_norm))
    jl_only_cols = sorted(set(jl_norm) - set(r_norm) - set(py_norm))

    if r_only:
        print(f"\n  Columns only in R ({len(r_only)}):")
        for c in r_only[:20]:
            print(f"    - {r_norm[c]}")
        if len(r_only) > 20:
            print(f"    ... and {len(r_only) - 20} more")

    if py_only:
        print(f"\n  Columns only in Python ({len(py_only)}):")
        for c in py_only[:20]:
            print(f"    - {py_norm[c]}")
        if len(py_only) > 20:
            print(f"    ... and {len(py_only) - 20} more")

    if jl_only_cols:
        print(f"\n  Columns only in Julia ({len(jl_only_cols)}):")
        for c in jl_only_cols[:20]:
            print(f"    - {jl_norm[c]}")
        if len(jl_only_cols) > 20:
            print(f"    ... and {len(jl_only_cols) - 20} more")

    # ── Three-way Spearman comparison on common columns ──
    print(f"\n{'=' * 80}")
    print("SPEARMAN CORRELATION COMPARISON (common columns, common gages)")
    print(f"{'=' * 80}")

    results = []
    for norm_col in sorted(common_norms):
        r_col = r_norm[norm_col]
        py_col = py_norm[norm_col]
        jl_col = jl_norm[norm_col]

        r_vals = pd.to_numeric(r_df[r_col], errors="coerce").values
        py_vals = pd.to_numeric(py_df[py_col], errors="coerce").values
        jl_vals = pd.to_numeric(jl_df[jl_col], errors="coerce").values

        rp_rho, rp_pval, rp_n = spearman_corr(r_vals, py_vals)
        rj_rho, rj_pval, rj_n = spearman_corr(r_vals, jl_vals)
        pj_rho, pj_pval, pj_n = spearman_corr(py_vals, jl_vals)

        r_na = int((~np.isfinite(r_vals)).sum())
        py_na = int((~np.isfinite(py_vals)).sum())
        jl_na = int((~np.isfinite(jl_vals)).sum())

        na_mismatch_rp = count_na_mismatch(r_vals, py_vals)
        na_mismatch_rj = count_na_mismatch(r_vals, jl_vals)
        na_mismatch_pj = count_na_mismatch(py_vals, jl_vals)

        results.append({
            "column": norm_col,
            "r_col": r_col,
            "py_col": py_col,
            "jl_col": jl_col,
            "rp_rho": rp_rho,
            "rj_rho": rj_rho,
            "pj_rho": pj_rho,
            "rp_n": rp_n,
            "rj_n": rj_n,
            "pj_n": pj_n,
            "r_na_count": r_na,
            "py_na_count": py_na,
            "jl_na_count": jl_na,
            "na_mismatch_rp": na_mismatch_rp,
            "na_mismatch_rj": na_mismatch_rj,
            "na_mismatch_pj": na_mismatch_pj,
            "min_rho": min(
                rp_rho if np.isfinite(rp_rho) else -999,
                rj_rho if np.isfinite(rj_rho) else -999,
                pj_rho if np.isfinite(pj_rho) else -999,
            ),
        })

    res_df = pd.DataFrame(results)

    # Save full results
    res_df.to_csv(BENCHMARK_DIR / "three_way_comparison.csv", index=False)
    print(f"\nFull results saved to benchmarks/three_way_comparison.csv")

    # ── Summary statistics ──
    valid = res_df.dropna(subset=["rp_rho", "rj_rho", "pj_rho"])
    print(f"\n  Columns with valid correlations across all 3 pairs: {len(valid)}")

    for pair, label in [("rp_rho", "R vs Python"), ("rj_rho", "R vs Julia"), ("pj_rho", "Python vs Julia")]:
        vals = valid[pair]
        print(f"\n  {label}:")
        print(f"    Mean rho:   {vals.mean():.6f}")
        print(f"    Median rho: {vals.median():.6f}")
        print(f"    Min rho:    {vals.min():.6f}")
        print(f"    Cols < 0.99:  {(vals < 0.99).sum()}")
        print(f"    Cols < 0.95:  {(vals < 0.95).sum()}")
        print(f"    Cols < 0.90:  {(vals < 0.90).sum()}")

    # ── Columns with poor agreement (< 0.99 in any pair) ──
    threshold = 0.99
    poor = res_df[res_df["min_rho"] < threshold].sort_values("min_rho")

    print(f"\n{'=' * 80}")
    print(f"COLUMNS WITH MIN SPEARMAN < {threshold} IN ANY PAIR ({len(poor)} columns)")
    print(f"{'=' * 80}")

    if len(poor) > 0:
        print(f"\n{'Column':<45} {'R-Py':>8} {'R-Jl':>8} {'Py-Jl':>8} {'R_NA':>6} {'Py_NA':>6} {'Jl_NA':>6}")
        print("-" * 97)
        for _, row in poor.iterrows():
            rp = f"{row['rp_rho']:.4f}" if np.isfinite(row['rp_rho']) else "  N/A "
            rj = f"{row['rj_rho']:.4f}" if np.isfinite(row['rj_rho']) else "  N/A "
            pj = f"{row['pj_rho']:.4f}" if np.isfinite(row['pj_rho']) else "  N/A "
            print(f"{row['column']:<45} {rp:>8} {rj:>8} {pj:>8} {row['r_na_count']:>6} {row['py_na_count']:>6} {row['jl_na_count']:>6}")

    # ── Categorize problems by signature group ──
    print(f"\n{'=' * 80}")
    print("PROBLEM SUMMARY BY SIGNATURE GROUP")
    print(f"{'=' * 80}")

    # Categorize all poor columns
    poor_by_group = {}
    for _, row in poor.iterrows():
        base = get_base_metric(row["column"])
        cat = categorize_metric(base)
        poor_by_group.setdefault(cat, []).append(row)

    for cat in sorted(poor_by_group.keys()):
        rows = poor_by_group[cat]
        print(f"\n  {cat} ({len(rows)} problematic columns):")
        for row in rows[:5]:
            rp = f"{row['rp_rho']:.4f}" if np.isfinite(row['rp_rho']) else "N/A"
            rj = f"{row['rj_rho']:.4f}" if np.isfinite(row['rj_rho']) else "N/A"
            pj = f"{row['pj_rho']:.4f}" if np.isfinite(row['pj_rho']) else "N/A"
            print(f"    {row['column']:<40} R-Py={rp}  R-Jl={rj}  Py-Jl={pj}")
        if len(rows) > 5:
            print(f"    ... and {len(rows) - 5} more")

    # ── Columns with perfect or near-perfect agreement ──
    perfect = res_df[res_df["min_rho"] >= 0.999]
    good = res_df[(res_df["min_rho"] >= threshold) & (res_df["min_rho"] < 0.999)]
    print(f"\n{'=' * 80}")
    print("AGREEMENT SUMMARY")
    print(f"{'=' * 80}")
    print(f"  Perfect (rho >= 0.999 all pairs): {len(perfect)} columns")
    print(f"  Good (0.99 <= rho < 0.999):       {len(good)} columns")
    print(f"  Poor (rho < 0.99 any pair):        {len(poor)} columns")
    print(f"  Total compared:                    {len(res_df)} columns")

    # ── NA pattern analysis ──
    print(f"\n{'=' * 80}")
    print("NA PATTERN ANALYSIS (columns with > 100 NA mismatches)")
    print(f"{'=' * 80}")

    na_issues = res_df[(res_df["na_mismatch_rp"] > 100) | (res_df["na_mismatch_rj"] > 100) | (res_df["na_mismatch_pj"] > 100)]
    na_issues = na_issues.sort_values("na_mismatch_rp", ascending=False)

    if len(na_issues) > 0:
        print(f"\n{'Column':<40} {'NA_RP':>8} {'NA_RJ':>8} {'NA_PJ':>8} {'R_NA':>6} {'Py_NA':>6} {'Jl_NA':>6}")
        print("-" * 86)
        for _, row in na_issues.head(40).iterrows():
            print(f"{row['column']:<40} {row['na_mismatch_rp']:>8} {row['na_mismatch_rj']:>8} {row['na_mismatch_pj']:>8} {row['r_na_count']:>6} {row['py_na_count']:>6} {row['jl_na_count']:>6}")
        if len(na_issues) > 40:
            print(f"... and {len(na_issues) - 40} more columns")
    else:
        print("  No columns with > 100 NA mismatches")

    # ── Pairwise columns (columns in 2 but not all 3) ──
    print(f"\n{'=' * 80}")
    print("PAIRWISE-ONLY COLUMNS")
    print(f"{'=' * 80}")

    for pair_cols, label, df_a, df_b, name_a, name_b, norm_a, norm_b in [
        (pj_only, "Python + Julia only (not in R)", py_df, jl_df, "Python", "Julia", py_norm, jl_norm),
        (rp_only, "R + Python only (not in Julia)", r_df, py_df, "R", "Python", r_norm, py_norm),
        (rj_only, "R + Julia only (not in Python)", r_df, jl_df, "R", "Julia", r_norm, jl_norm),
    ]:
        if pair_cols:
            print(f"\n  {label} ({len(pair_cols)} columns):")
            for nc in sorted(pair_cols)[:15]:
                col_a = norm_a[nc]
                col_b = norm_b[nc]
                va = pd.to_numeric(df_a[col_a], errors="coerce").values
                vb = pd.to_numeric(df_b[col_b], errors="coerce").values
                rho, _, n = spearman_corr(va, vb)
                rho_str = f"{rho:.4f}" if np.isfinite(rho) else "N/A"
                print(f"    {nc:<45} {name_a}-{name_b} rho={rho_str} (n={n})")
            if len(pair_cols) > 15:
                print(f"    ... and {len(pair_cols) - 15} more")

    # ── Deep dive: worst metrics with per-gage outlier analysis ──
    print(f"\n{'=' * 80}")
    print("DEEP DIVE: WORST 10 METRICS (per-gage outlier analysis)")
    print(f"{'=' * 80}")

    worst = poor.head(10)
    for _, row in worst.iterrows():
        col = row["column"]
        r_col = row["r_col"]
        py_col = row["py_col"]
        jl_col = row["jl_col"]

        r_vals = pd.to_numeric(r_df[r_col], errors="coerce")
        py_vals = pd.to_numeric(py_df[py_col], errors="coerce")
        jl_vals = pd.to_numeric(jl_df[jl_col], errors="coerce")

        print(f"\n  {col}:")
        print(f"    R range:  [{r_vals.min():.4f}, {r_vals.max():.4f}], median={r_vals.median():.4f}, NAs={r_vals.isna().sum()}")
        print(f"    Py range: [{py_vals.min():.4f}, {py_vals.max():.4f}], median={py_vals.median():.4f}, NAs={py_vals.isna().sum()}")
        print(f"    Jl range: [{jl_vals.min():.4f}, {jl_vals.max():.4f}], median={jl_vals.median():.4f}, NAs={jl_vals.isna().sum()}")

        # Find biggest discrepancies R vs Python
        diff_rp = (r_vals - py_vals).abs()
        top_rp = diff_rp.nlargest(5)
        if len(top_rp) > 0:
            print(f"    Top 5 R-Python discrepancies:")
            for gage_id in top_rp.index:
                rv = r_vals.get(gage_id, np.nan)
                pv = py_vals.get(gage_id, np.nan)
                jv = jl_vals.get(gage_id, np.nan)
                print(f"      gage={gage_id}: R={rv:.6f}, Py={pv:.6f}, Jl={jv:.6f}")

    # ── Generate markdown report ──
    generate_three_way_report(res_df, r_df, py_df, jl_df,
                              common_norms, r_norm, py_norm, jl_norm,
                              pj_only, rp_only, rj_only)

    print(f"\n{'=' * 80}")
    print("DONE")
    print(f"{'=' * 80}")

    return 0 if len(poor) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

"""
Cross-Language Validation Script

Compares Python signature outputs against R golden reference outputs.
Generates validation report with relative errors and pass/fail status.
"""

import sys
import os
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd

# Add project root for imports (docs/benchmarks/ -> project root)
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))

from streamflow_signatures import (
    read_parquet,
    add_water_year_columns,
    generate_stats,
    calculate_flow_vols_by_year,
    analyze_flashiness_trends,
    analyze_flow_timing_trends,
    analyze_fdc_trends,
    analyze_baseflow_indices,
    analyze_recession_parameters,
    calculate_pulse_metrics,
    analyze_Q_PPT_relationships,
    calculate_streamflow_elasticity,
    calculate_qp_seasonality,
    calculate_average_storage,
)


# Configuration
RELATIVE_TOLERANCE = 0.0001  # 0.01% relative error threshold
ABSOLUTE_TOLERANCE = 1e-10   # For values near zero


def load_golden_outputs(golden_path: Path) -> pd.DataFrame:
    """Load R-generated golden reference outputs."""
    if not golden_path.exists():
        raise FileNotFoundError(f"Golden outputs not found: {golden_path}")
    return pd.read_csv(golden_path)


def load_test_data(test_data_dir: Path) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """Load test streamflow and climate data."""
    streamflow_path = test_data_dir / "sample_streamflow.parquet"
    climate_path = test_data_dir / "sample_climate.parquet"

    if not streamflow_path.exists():
        raise FileNotFoundError(f"Test streamflow data not found: {streamflow_path}")

    streamflow = read_parquet(streamflow_path)
    streamflow = add_water_year_columns(streamflow, date_col="date")

    climate = None
    if climate_path.exists():
        climate = read_parquet(climate_path)
        climate = add_water_year_columns(climate, date_col="date")

    return streamflow, climate


def calculate_all_signatures(
    gage_data: pd.DataFrame,
    climate_data: Optional[pd.DataFrame] = None
) -> Dict[str, float]:
    """Calculate all signatures for a single gage."""
    results = {}

    # Simple signatures (no climate required)
    try:
        results.update(calculate_flow_vols_by_year(gage_data))
    except Exception as e:
        print(f"  Warning: flow_vols failed: {e}")

    try:
        results.update(analyze_flashiness_trends(gage_data))
    except Exception as e:
        print(f"  Warning: flashiness failed: {e}")

    try:
        results.update(analyze_flow_timing_trends(gage_data))
    except Exception as e:
        print(f"  Warning: timing failed: {e}")

    try:
        results.update(analyze_fdc_trends(gage_data))
    except Exception as e:
        print(f"  Warning: fdc failed: {e}")

    # Complex signatures
    try:
        results.update(analyze_baseflow_indices(gage_data))
    except Exception as e:
        print(f"  Warning: baseflow failed: {e}")

    try:
        results.update(analyze_recession_parameters(gage_data))
    except Exception as e:
        print(f"  Warning: recession failed: {e}")

    try:
        results.update(calculate_pulse_metrics(gage_data))
    except Exception as e:
        print(f"  Warning: pulses failed: {e}")

    # Climate-dependent signatures (if climate data available)
    if climate_data is not None and len(climate_data) > 0:
        # Merge climate data with streamflow
        merged = pd.merge(
            gage_data,
            climate_data[["gage_id", "date", "PPT"]],
            on=["gage_id", "date"],
            how="left"
        )

        if "PPT" in merged.columns and merged["PPT"].notna().sum() > 0:
            try:
                results.update(analyze_Q_PPT_relationships(merged))
            except Exception as e:
                print(f"  Warning: runoff_ratios failed: {e}")

            try:
                results.update(calculate_streamflow_elasticity(merged))
            except Exception as e:
                print(f"  Warning: elasticity failed: {e}")

            try:
                results.update(calculate_qp_seasonality(merged))
            except Exception as e:
                print(f"  Warning: qp_seasonality failed: {e}")

            try:
                results.update(calculate_average_storage(merged))
            except Exception as e:
                print(f"  Warning: storage failed: {e}")

    return results


def run_python_extraction(
    streamflow: pd.DataFrame,
    climate: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """Run Python signature extraction on all test gages."""
    gages = streamflow["gage_id"].unique()
    all_results = []

    print(f"\nProcessing {len(gages)} gages with Python...")

    for i, gage_id in enumerate(gages):
        print(f"  [{i+1}/{len(gages)}] Processing gage {gage_id}")

        gage_data = streamflow[streamflow["gage_id"] == gage_id].copy()

        # Get climate data for this gage if available
        gage_climate = None
        if climate is not None:
            gage_climate = climate[climate["gage_id"] == gage_id].copy()

        signatures = calculate_all_signatures(gage_data, gage_climate)
        signatures["gage_id"] = gage_id
        all_results.append(signatures)

    return pd.DataFrame(all_results)


def compute_relative_error(actual: float, expected: float) -> float:
    """Compute relative error between two values."""
    if pd.isna(actual) and pd.isna(expected):
        return 0.0  # Both NA is a match
    if pd.isna(actual) or pd.isna(expected):
        return np.inf  # One NA, one not is a mismatch
    if abs(expected) < ABSOLUTE_TOLERANCE:
        if abs(actual) < ABSOLUTE_TOLERANCE:
            return 0.0  # Both near zero
        return np.inf  # Expected zero, actual not
    return abs(actual - expected) / abs(expected)


def compare_outputs(
    python_results: pd.DataFrame,
    golden_results: pd.DataFrame,
    tolerance: float = RELATIVE_TOLERANCE
) -> Dict:
    """Compare Python outputs against R golden outputs."""

    # Find common columns (excluding gage_id)
    python_cols = set(python_results.columns) - {"gage_id"}
    golden_cols = set(golden_results.columns) - {"gage_id"}
    common_cols = python_cols & golden_cols

    missing_in_python = golden_cols - python_cols
    extra_in_python = python_cols - golden_cols

    print(f"\nColumn comparison:")
    print(f"  Common columns: {len(common_cols)}")
    print(f"  Missing in Python: {len(missing_in_python)}")
    print(f"  Extra in Python: {len(extra_in_python)}")

    if missing_in_python:
        print(f"  Missing columns: {sorted(missing_in_python)[:10]}...")

    # Merge on gage_id
    merged = pd.merge(
        python_results,
        golden_results,
        on="gage_id",
        suffixes=("_py", "_r")
    )

    # Calculate errors for each column
    column_results = {}
    all_errors = []

    for col in sorted(common_cols):
        py_col = f"{col}_py"
        r_col = f"{col}_r"

        if py_col not in merged.columns or r_col not in merged.columns:
            continue

        errors = []
        for _, row in merged.iterrows():
            error = compute_relative_error(row[py_col], row[r_col])
            errors.append(error)

        errors = np.array(errors)
        valid_errors = errors[np.isfinite(errors)]

        if len(valid_errors) > 0:
            max_error = np.max(valid_errors)
            mean_error = np.mean(valid_errors)
            passed = max_error <= tolerance
        else:
            max_error = np.nan
            mean_error = np.nan
            passed = len(errors) == 0

        column_results[col] = {
            "max_error": max_error,
            "mean_error": mean_error,
            "n_comparisons": len(errors),
            "n_valid": len(valid_errors),
            "n_both_na": sum(np.array(errors) == 0),
            "passed": passed
        }

        if len(valid_errors) > 0:
            all_errors.extend(valid_errors.tolist())

    # Summary statistics
    n_passed = sum(1 for r in column_results.values() if r["passed"])
    n_failed = sum(1 for r in column_results.values() if not r["passed"])

    summary = {
        "n_gages": len(merged),
        "n_columns_compared": len(common_cols),
        "n_columns_passed": n_passed,
        "n_columns_failed": n_failed,
        "missing_in_python": sorted(missing_in_python),
        "extra_in_python": sorted(extra_in_python),
        "overall_max_error": max(all_errors) if all_errors else np.nan,
        "overall_mean_error": np.mean(all_errors) if all_errors else np.nan,
        "tolerance": tolerance,
        "passed": n_failed == 0 and len(missing_in_python) == 0
    }

    return {
        "summary": summary,
        "column_results": column_results
    }


def generate_report(
    comparison: Dict,
    output_path: Path,
    python_results: pd.DataFrame,
    golden_results: pd.DataFrame
) -> None:
    """Generate validation report in markdown format."""
    summary = comparison["summary"]
    column_results = comparison["column_results"]

    with open(output_path, "w") as f:
        f.write("# Cross-Language Validation Report\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Overall status
        status = "PASSED" if summary["passed"] else "FAILED"
        status_emoji = "✅" if summary["passed"] else "❌"
        f.write(f"## Overall Status: {status_emoji} {status}\n\n")

        # Summary statistics
        f.write("## Summary\n\n")
        f.write(f"| Metric | Value |\n")
        f.write(f"|--------|-------|\n")
        f.write(f"| Gages compared | {summary['n_gages']} |\n")
        f.write(f"| Columns compared | {summary['n_columns_compared']} |\n")
        f.write(f"| Columns passed | {summary['n_columns_passed']} |\n")
        f.write(f"| Columns failed | {summary['n_columns_failed']} |\n")
        f.write(f"| Tolerance | {summary['tolerance']:.4%} |\n")
        f.write(f"| Overall max error | {summary['overall_max_error']:.6%} |\n")
        f.write(f"| Overall mean error | {summary['overall_mean_error']:.6%} |\n")
        f.write("\n")

        # Missing columns
        if summary["missing_in_python"]:
            f.write("## Missing Columns in Python\n\n")
            f.write("These columns are in R output but not in Python:\n\n")
            for col in summary["missing_in_python"]:
                f.write(f"- `{col}`\n")
            f.write("\n")

        # Failed columns
        failed_cols = {k: v for k, v in column_results.items() if not v["passed"]}
        if failed_cols:
            f.write("## Failed Columns\n\n")
            f.write("| Column | Max Error | Mean Error | N Valid |\n")
            f.write("|--------|-----------|------------|----------|\n")
            for col, result in sorted(failed_cols.items(), key=lambda x: -x[1]["max_error"]):
                f.write(f"| `{col}` | {result['max_error']:.6%} | {result['mean_error']:.6%} | {result['n_valid']} |\n")
            f.write("\n")

        # Passed columns (summary)
        passed_cols = {k: v for k, v in column_results.items() if v["passed"]}
        if passed_cols:
            f.write("## Passed Columns\n\n")
            f.write(f"All {len(passed_cols)} columns passed with max error < {summary['tolerance']:.4%}\n\n")

            # Show worst 10 passed columns
            sorted_passed = sorted(passed_cols.items(), key=lambda x: -x[1]["max_error"])[:10]
            if sorted_passed:
                f.write("Top 10 columns by error (still passing):\n\n")
                f.write("| Column | Max Error | Mean Error |\n")
                f.write("|--------|-----------|------------|\n")
                for col, result in sorted_passed:
                    f.write(f"| `{col}` | {result['max_error']:.6%} | {result['mean_error']:.6%} |\n")
            f.write("\n")

        # Extra columns (informational)
        if summary["extra_in_python"]:
            f.write("## Extra Columns in Python\n\n")
            f.write("These columns are in Python output but not in R golden output:\n\n")
            for col in summary["extra_in_python"]:
                f.write(f"- `{col}`\n")
            f.write("\n")

        f.write("---\n")
        f.write("*Generated by docs/benchmarks/compare_outputs.py*\n")

    print(f"\nReport saved to: {output_path}")


def main():
    """Main validation entry point."""
    print("=" * 60)
    print("CROSS-LANGUAGE VALIDATION")
    print("=" * 60)

    # Paths (docs/benchmarks/ -> project root)
    base_dir = Path(__file__).parent.parent.parent
    test_data_dir = base_dir / "test-data"
    golden_output_dir = base_dir / "golden-outputs"

    golden_path = golden_output_dir / "expected_signatures.csv"
    golden_climate_path = golden_output_dir / "expected_signatures_climate.csv"

    # Check if test data exists
    if not test_data_dir.exists():
        print(f"\nERROR: Test data directory not found: {test_data_dir}")
        print("Please run R/tests/generate_golden_outputs.R first to create test data.")
        sys.exit(1)

    if not golden_path.exists():
        print(f"\nERROR: Golden outputs not found: {golden_path}")
        print("Please run R/tests/generate_golden_outputs.R first to generate golden outputs.")
        sys.exit(1)

    # Load data
    print("\nLoading data...")
    streamflow, climate = load_test_data(test_data_dir)
    print(f"  Streamflow: {len(streamflow)} rows, {len(streamflow['gage_id'].unique())} gages")
    if climate is not None:
        print(f"  Climate: {len(climate)} rows")

    # Load golden outputs
    golden = load_golden_outputs(golden_path)
    print(f"  Golden outputs: {len(golden)} gages, {len(golden.columns)} columns")

    # Run Python extraction
    python_results = run_python_extraction(streamflow, climate)
    print(f"\nPython results: {len(python_results)} gages, {len(python_results.columns)} columns")

    # Save Python results for debugging
    python_output_path = validation_dir / "python_signatures.csv"
    python_results.to_csv(python_output_path, index=False)
    print(f"Python results saved to: {python_output_path}")

    # Compare outputs
    print("\nComparing outputs...")
    comparison = compare_outputs(python_results, golden, RELATIVE_TOLERANCE)

    # Generate report
    report_path = validation_dir / "validation_report.md"
    generate_report(comparison, report_path, python_results, golden)

    # Print summary
    summary = comparison["summary"]
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Status: {'PASSED' if summary['passed'] else 'FAILED'}")
    print(f"Columns: {summary['n_columns_passed']} passed, {summary['n_columns_failed']} failed")
    print(f"Max error: {summary['overall_max_error']:.6%}")
    print(f"Mean error: {summary['overall_mean_error']:.6%}")

    if not summary["passed"]:
        print("\nTo fix failing columns, check validation_report.md for details.")
        sys.exit(1)

    print("\nValidation complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

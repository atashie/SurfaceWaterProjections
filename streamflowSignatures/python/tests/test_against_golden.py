"""
Test Python signatures against R golden outputs.

This test suite validates that the Python implementation produces
results within tolerance of the canonical R implementation.

Run with: pytest tests/test_against_golden.py -v
"""

import sys
from pathlib import Path
import pytest
import numpy as np
import pandas as pd

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from streamflow_signatures import (
    read_parquet,
    add_water_year_columns,
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
RELATIVE_TOLERANCE = 0.0001  # 0.01% relative error
ABSOLUTE_TOLERANCE = 1e-10

# Paths
BASE_DIR = Path(__file__).parent.parent.parent
TEST_DATA_DIR = BASE_DIR / "test-data"
GOLDEN_OUTPUT_DIR = BASE_DIR / "golden-outputs"


def relative_error(actual: float, expected: float) -> float:
    """Calculate relative error between two values."""
    if pd.isna(actual) and pd.isna(expected):
        return 0.0
    if pd.isna(actual) or pd.isna(expected):
        return np.inf
    if abs(expected) < ABSOLUTE_TOLERANCE:
        if abs(actual) < ABSOLUTE_TOLERANCE:
            return 0.0
        return np.inf
    return abs(actual - expected) / abs(expected)


@pytest.fixture(scope="module")
def test_data():
    """Load test streamflow data."""
    streamflow_path = TEST_DATA_DIR / "sample_streamflow.parquet"
    if not streamflow_path.exists():
        pytest.skip(f"Test data not found: {streamflow_path}")

    df = read_parquet(streamflow_path)
    df = add_water_year_columns(df, date_col="date")
    return df


@pytest.fixture(scope="module")
def climate_data():
    """Load test climate data if available."""
    climate_path = TEST_DATA_DIR / "sample_climate.parquet"
    if not climate_path.exists():
        return None

    df = read_parquet(climate_path)
    df = add_water_year_columns(df, date_col="date")
    return df


@pytest.fixture(scope="module")
def golden_outputs():
    """Load R golden reference outputs."""
    golden_path = GOLDEN_OUTPUT_DIR / "expected_signatures.csv"
    if not golden_path.exists():
        pytest.skip(f"Golden outputs not found: {golden_path}")

    return pd.read_csv(golden_path)


@pytest.fixture(scope="module")
def golden_outputs_climate():
    """Load R golden reference outputs with climate signatures."""
    golden_path = GOLDEN_OUTPUT_DIR / "expected_signatures_climate.csv"
    if not golden_path.exists():
        return None

    return pd.read_csv(golden_path)


def get_gage_data(test_data, gage_id):
    """Extract data for a single gage."""
    return test_data[test_data["gage_id"] == gage_id].copy()


def get_gage_golden(golden_outputs, gage_id):
    """Get golden outputs for a single gage."""
    # Handle potential type mismatches in gage_id
    mask = golden_outputs["gage_id"].astype(str) == str(gage_id)
    if mask.sum() == 0:
        return None
    return golden_outputs[mask].iloc[0].to_dict()


class TestFlowVolumes:
    """Test flow volume signatures."""

    def test_flow_volumes_match_golden(self, test_data, golden_outputs):
        """Test that flow volume metrics match R outputs."""
        gages = test_data["gage_id"].unique()

        for gage_id in gages[:5]:  # Test subset for speed
            gage_data = get_gage_data(test_data, gage_id)
            golden = get_gage_golden(golden_outputs, gage_id)

            if golden is None:
                continue

            results = calculate_flow_vols_by_year(gage_data)

            # Check key metrics
            for metric in ["Qann_mean", "Qann_senn_slp", "Q50_mean", "Q95_mean"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestFlashiness:
    """Test flashiness signatures."""

    def test_flashiness_match_golden(self, test_data, golden_outputs):
        """Test that flashiness metrics match R outputs."""
        gages = test_data["gage_id"].unique()

        for gage_id in gages[:5]:
            gage_data = get_gage_data(test_data, gage_id)
            golden = get_gage_golden(golden_outputs, gage_id)

            if golden is None:
                continue

            results = analyze_flashiness_trends(gage_data)

            for metric in ["flashinessRB_mean", "flashinessRB_senn_slp"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestTiming:
    """Test flow timing signatures."""

    def test_timing_match_golden(self, test_data, golden_outputs):
        """Test that timing metrics match R outputs."""
        gages = test_data["gage_id"].unique()

        for gage_id in gages[:5]:
            gage_data = get_gage_data(test_data, gage_id)
            golden = get_gage_golden(golden_outputs, gage_id)

            if golden is None:
                continue

            results = analyze_flow_timing_trends(gage_data)

            for metric in ["D50_day_mean", "D50_day_senn_slp", "Dmax_mean"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestFDC:
    """Test flow duration curve signatures."""

    def test_fdc_match_golden(self, test_data, golden_outputs):
        """Test that FDC metrics match R outputs."""
        gages = test_data["gage_id"].unique()

        for gage_id in gages[:5]:
            gage_data = get_gage_data(test_data, gage_id)
            golden = get_gage_golden(golden_outputs, gage_id)

            if golden is None:
                continue

            results = analyze_fdc_trends(gage_data)

            for metric in ["FDCall_mean", "FDC90th_mean", "FDCmid_mean"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestBaseflow:
    """Test baseflow signatures."""

    def test_baseflow_match_golden(self, test_data, golden_outputs):
        """Test that baseflow indices match R outputs."""
        gages = test_data["gage_id"].unique()

        for gage_id in gages[:5]:
            gage_data = get_gage_data(test_data, gage_id)
            golden = get_gage_golden(golden_outputs, gage_id)

            if golden is None:
                continue

            results = analyze_baseflow_indices(gage_data)

            for metric in ["BFI_Eckhardt_mean", "BFI_LyneHollick_mean"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestRecession:
    """Test recession signatures."""

    def test_recession_match_golden(self, test_data, golden_outputs):
        """Test that recession metrics match R outputs."""
        gages = test_data["gage_id"].unique()

        for gage_id in gages[:5]:
            gage_data = get_gage_data(test_data, gage_id)
            golden = get_gage_golden(golden_outputs, gage_id)

            if golden is None:
                continue

            results = analyze_recession_parameters(gage_data)

            for metric in ["log_a_pointcloud_mean", "b_pointcloud_mean"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestPulses:
    """Test pulse signatures."""

    def test_pulses_match_golden(self, test_data, golden_outputs):
        """Test that pulse metrics match R outputs."""
        gages = test_data["gage_id"].unique()

        for gage_id in gages[:5]:
            gage_data = get_gage_data(test_data, gage_id)
            golden = get_gage_golden(golden_outputs, gage_id)

            if golden is None:
                continue

            results = calculate_pulse_metrics(gage_data)

            for metric in ["n_high_pulses_year_mean", "TQmean_mean"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestClimateSignatures:
    """Test climate-dependent signatures."""

    def test_runoff_ratios_match_golden(self, test_data, climate_data, golden_outputs_climate):
        """Test that runoff ratios match R outputs."""
        if climate_data is None or golden_outputs_climate is None:
            pytest.skip("Climate data or golden outputs not available")

        gages = test_data["gage_id"].unique()

        for gage_id in gages[:3]:
            gage_data = get_gage_data(test_data, gage_id)
            gage_climate = climate_data[climate_data["gage_id"] == gage_id].copy()
            golden = get_gage_golden(golden_outputs_climate, gage_id)

            if golden is None or len(gage_climate) == 0:
                continue

            merged = pd.merge(
                gage_data,
                gage_climate[["gage_id", "date", "PPT"]],
                on=["gage_id", "date"],
                how="left"
            )

            if "PPT" not in merged.columns or merged["PPT"].isna().all():
                continue

            results = analyze_Q_PPT_relationships(merged)

            for metric in ["annual_runoff_ratio_mean"]:
                if metric in golden and metric in results:
                    error = relative_error(results[metric], golden[metric])
                    assert error <= RELATIVE_TOLERANCE, \
                        f"Gage {gage_id}, {metric}: error {error:.6%}"


class TestAllSignatures:
    """Integration test for all signatures."""

    def test_all_signatures_exist(self, test_data, golden_outputs):
        """Test that Python produces all expected signature columns."""
        gage_id = test_data["gage_id"].iloc[0]
        gage_data = get_gage_data(test_data, gage_id)

        # Collect all signatures
        all_sigs = {}
        all_sigs.update(calculate_flow_vols_by_year(gage_data))
        all_sigs.update(analyze_flashiness_trends(gage_data))
        all_sigs.update(analyze_flow_timing_trends(gage_data))
        all_sigs.update(analyze_fdc_trends(gage_data))
        all_sigs.update(analyze_baseflow_indices(gage_data))
        all_sigs.update(analyze_recession_parameters(gage_data))
        all_sigs.update(calculate_pulse_metrics(gage_data))

        # Check that we have the expected signature suffixes
        expected_suffixes = [
            "_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
            "_mk_rho", "_mk_pval", "_mean", "_median"
        ]

        # Check Qann metrics as example
        for suffix in expected_suffixes:
            assert f"Qann{suffix}" in all_sigs, f"Missing Qann{suffix}"

        # Check total count
        golden_cols = set(golden_outputs.columns) - {"gage_id"}
        python_cols = set(all_sigs.keys())

        # Allow some tolerance for missing climate signatures
        min_expected = len(golden_cols) * 0.8  # At least 80% of columns
        assert len(python_cols) >= min_expected, \
            f"Python has {len(python_cols)} columns, expected at least {min_expected}"

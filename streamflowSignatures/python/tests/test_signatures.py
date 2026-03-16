"""
Unit tests for streamflow signatures using synthetic data.

These tests verify the Python implementation works correctly
without requiring R golden outputs.

Run with: pytest tests/test_signatures.py -v
"""

import numpy as np
import pandas as pd
import pytest
from datetime import datetime, timedelta

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from streamflow_signatures import (
    generate_stats,
    theil_sen_slope,
    mann_kendall_test,
    add_water_year_columns,
    calculate_flow_vols_by_year,
    analyze_flashiness_trends,
    analyze_flow_timing_trends,
    analyze_fdc_trends,
    analyze_baseflow_indices,
    eckhardt_filter,
    lyne_hollick_filter,
    analyze_recession_parameters,
    calculate_pulse_metrics,
)


# ============================================================================
# Fixtures for synthetic data
# ============================================================================

@pytest.fixture
def synthetic_streamflow():
    """Generate synthetic streamflow data for testing."""
    np.random.seed(42)

    # Generate 25 years of daily data
    start_date = datetime(1990, 10, 1)
    n_days = 365 * 25

    dates = [start_date + timedelta(days=i) for i in range(n_days)]

    # Create seasonal flow pattern with some noise
    doy = np.array([(d - datetime(d.year, 1, 1)).days + 1 for d in dates])

    # Base flow with seasonal variation (peak in spring)
    base_flow = 2.0 + 1.5 * np.sin(2 * np.pi * (doy - 100) / 365)
    noise = np.random.exponential(0.5, n_days)
    Q = np.maximum(0.01, base_flow + noise)

    df = pd.DataFrame({
        "gage_id": "TEST001",
        "date": dates,
        "Q": Q
    })

    df = add_water_year_columns(df, date_col="date")
    return df


@pytest.fixture
def synthetic_streamflow_with_climate():
    """Generate synthetic streamflow + climate data."""
    np.random.seed(42)

    start_date = datetime(1990, 10, 1)
    n_days = 365 * 25

    dates = [start_date + timedelta(days=i) for i in range(n_days)]
    doy = np.array([(d - datetime(d.year, 1, 1)).days + 1 for d in dates])

    # Flow
    base_flow = 2.0 + 1.5 * np.sin(2 * np.pi * (doy - 100) / 365)
    noise = np.random.exponential(0.5, n_days)
    Q = np.maximum(0.01, base_flow + noise)

    # Precipitation - correlated with flow
    PPT = np.maximum(0, Q * 2.5 + np.random.exponential(1.0, n_days))

    df = pd.DataFrame({
        "gage_id": "TEST001",
        "date": dates,
        "Q": Q,
        "PPT": PPT
    })

    df = add_water_year_columns(df, date_col="date")
    return df


# ============================================================================
# Stats module tests
# ============================================================================

class TestStats:
    """Test statistical functions."""

    def test_theil_sen_positive_trend(self):
        """Test Theil-Sen detects positive trend."""
        x = np.arange(20)
        y = 2 * x + 1 + np.random.normal(0, 0.5, 20)

        slope = theil_sen_slope(x, y)

        assert slope > 1.5, f"Expected slope ~2, got {slope}"
        assert slope < 2.5, f"Expected slope ~2, got {slope}"

    def test_theil_sen_no_trend(self):
        """Test Theil-Sen with no trend."""
        np.random.seed(42)
        x = np.arange(20)
        y = 5.0 + np.random.normal(0, 0.1, 20)

        slope = theil_sen_slope(x, y)

        assert abs(slope) < 0.1, f"Expected slope ~0, got {slope}"

    def test_mann_kendall_positive_trend(self):
        """Test Mann-Kendall detects positive trend."""
        x = np.arange(30)
        y = 2 * x + np.random.normal(0, 1, 30)

        tau, pval = mann_kendall_test(y)

        assert tau > 0.7, f"Expected strong positive tau, got {tau}"
        assert pval < 0.05, f"Expected significant p-value, got {pval}"

    def test_generate_stats_output_format(self):
        """Test generate_stats produces correct output structure."""
        data = pd.DataFrame({
            "water_year": range(2000, 2020),
            "metric1": np.random.random(20),
            "metric2": np.random.random(20)
        })

        stats = generate_stats(data, value_cols=["metric1", "metric2"])

        expected_suffixes = [
            "_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
            "_mk_rho", "_mk_pval", "_mean", "_median"
        ]

        for metric in ["metric1", "metric2"]:
            for suffix in expected_suffixes:
                key = f"{metric}{suffix}"
                assert key in stats, f"Missing {key}"
                assert isinstance(stats[key], (int, float, np.floating)), \
                    f"{key} is not numeric: {type(stats[key])}"


# ============================================================================
# Flow volumes tests
# ============================================================================

class TestFlowVolumes:
    """Test flow volume signatures."""

    def test_flow_volumes_output_structure(self, synthetic_streamflow):
        """Test flow volumes produces expected metrics."""
        results = calculate_flow_vols_by_year(synthetic_streamflow)

        # Check key metrics exist
        assert "Qann_mean" in results
        assert "Qann_senn_slp" in results
        assert "Q50_mean" in results
        assert "Q95_mean" in results

        # Check values are reasonable
        assert results["Qann_mean"] > 0, "Qann should be positive"
        assert results["Q50_mean"] > 0, "Q50 should be positive"

    def test_seasonal_volumes(self, synthetic_streamflow):
        """Test seasonal flow volumes."""
        results = calculate_flow_vols_by_year(synthetic_streamflow)

        for season in ["Qwin", "Qspr", "Qsum", "Qfal"]:
            assert f"{season}_mean" in results, f"Missing {season}_mean"
            assert results[f"{season}_mean"] > 0, f"{season} should be positive"


# ============================================================================
# Baseflow tests
# ============================================================================

class TestBaseflow:
    """Test baseflow signatures."""

    def test_eckhardt_filter_range(self):
        """Test Eckhardt filter produces values in [0, 1] ratio to Q."""
        np.random.seed(42)
        Q = np.maximum(0.01, np.random.exponential(2, 1000))

        baseflow = eckhardt_filter(Q)

        assert np.all(baseflow >= 0), "Baseflow should be non-negative"
        assert np.all(baseflow <= Q + 1e-10), "Baseflow should not exceed Q"

    def test_lyne_hollick_filter_range(self):
        """Test Lyne-Hollick filter produces values in valid range."""
        np.random.seed(42)
        Q = np.maximum(0.01, np.random.exponential(2, 1000))

        baseflow = lyne_hollick_filter(Q)

        assert np.all(baseflow >= 0), "Baseflow should be non-negative"
        assert np.all(baseflow <= Q + 1e-10), "Baseflow should not exceed Q"

    def test_baseflow_indices(self, synthetic_streamflow):
        """Test baseflow indices are in valid range [0, 1]."""
        results = analyze_baseflow_indices(synthetic_streamflow)

        for metric in ["BFI_Eckhardt_mean", "BFI_LyneHollick_mean"]:
            assert metric in results, f"Missing {metric}"
            if not np.isnan(results[metric]):
                assert 0 <= results[metric] <= 1, f"{metric} should be in [0, 1]"

    def test_baseflow_relationship(self, synthetic_streamflow):
        """Test BFI_Eckhardt < BFI_LyneHollick typically."""
        results = analyze_baseflow_indices(synthetic_streamflow)

        eck = results.get("BFI_Eckhardt_mean", np.nan)
        lh = results.get("BFI_LyneHollick_mean", np.nan)

        if not np.isnan(eck) and not np.isnan(lh):
            # Eckhardt typically produces lower values
            assert eck <= lh + 0.1, "BFI_Eckhardt should be <= BFI_LyneHollick"


# ============================================================================
# Flashiness tests
# ============================================================================

class TestFlashiness:
    """Test flashiness signatures."""

    def test_flashiness_output(self, synthetic_streamflow):
        """Test flashiness produces expected output."""
        results = analyze_flashiness_trends(synthetic_streamflow)

        assert "flashinessRB_mean" in results
        assert "flashinessRB_senn_slp" in results

        # Flashiness should be non-negative
        if not np.isnan(results["flashinessRB_mean"]):
            assert results["flashinessRB_mean"] >= 0

    def test_flashiness_constant_flow(self):
        """Test flashiness is zero for constant flow."""
        # Create constant flow data
        dates = pd.date_range("2000-10-01", periods=365*10, freq="D")
        df = pd.DataFrame({
            "gage_id": "CONST",
            "date": dates,
            "Q": 1.0  # Constant flow
        })
        df = add_water_year_columns(df, date_col="date")

        results = analyze_flashiness_trends(df)

        # Constant flow should have zero flashiness
        assert results["flashinessRB_mean"] == 0.0 or np.isnan(results["flashinessRB_mean"])


# ============================================================================
# Timing tests
# ============================================================================

class TestTiming:
    """Test flow timing signatures."""

    def test_timing_output(self, synthetic_streamflow):
        """Test timing produces expected metrics."""
        results = analyze_flow_timing_trends(synthetic_streamflow)

        # Check D-day metrics
        for pct in [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]:
            assert f"D{pct}_day_mean" in results, f"Missing D{pct}_day_mean"

        # Check Dmax
        assert "Dmax_mean" in results

    def test_timing_ordering(self, synthetic_streamflow):
        """Test D-day metrics are properly ordered."""
        results = analyze_flow_timing_trends(synthetic_streamflow)

        # D10 should be earlier than D50 should be earlier than D90
        d10 = results.get("D10_day_mean", np.nan)
        d50 = results.get("D50_day_mean", np.nan)
        d90 = results.get("D90_day_mean", np.nan)

        if not any(np.isnan([d10, d50, d90])):
            assert d10 < d50 < d90, "D-day metrics should be ordered"


# ============================================================================
# FDC tests
# ============================================================================

class TestFDC:
    """Test flow duration curve signatures."""

    def test_fdc_output(self, synthetic_streamflow):
        """Test FDC produces expected metrics."""
        results = analyze_fdc_trends(synthetic_streamflow)

        for metric in ["FDCall_mean", "FDC90th_mean", "FDCmid_mean"]:
            assert metric in results, f"Missing {metric}"


# ============================================================================
# Recession tests
# ============================================================================

class TestRecession:
    """Test recession signatures."""

    def test_recession_output(self, synthetic_streamflow):
        """Test recession produces expected metrics."""
        results = analyze_recession_parameters(synthetic_streamflow)

        # Check pointcloud metrics
        assert "log_a_pointcloud_mean" in results
        assert "b_pointcloud_mean" in results

    def test_recession_b_range(self, synthetic_streamflow):
        """Test recession exponent b is in reasonable range."""
        results = analyze_recession_parameters(synthetic_streamflow)

        b = results.get("b_pointcloud_mean", np.nan)
        if not np.isnan(b):
            # b typically between 0 and 3
            assert 0 < b < 5, f"b_pointcloud should be in reasonable range, got {b}"


# ============================================================================
# Pulse tests
# ============================================================================

class TestPulses:
    """Test pulse signatures."""

    def test_pulse_output(self, synthetic_streamflow):
        """Test pulses produces expected metrics."""
        results = calculate_pulse_metrics(synthetic_streamflow)

        assert "n_high_pulses_year_mean" in results
        assert "n_low_pulses_year_mean" in results
        assert "TQmean_mean" in results
        assert "Flow_Reversals_annual_mean" in results

    def test_tqmean_range(self, synthetic_streamflow):
        """Test TQmean is a valid percentage."""
        results = calculate_pulse_metrics(synthetic_streamflow)

        tqmean = results.get("TQmean_mean", np.nan)
        if not np.isnan(tqmean):
            assert 0 <= tqmean <= 100, f"TQmean should be in [0, 100], got {tqmean}"


# ============================================================================
# Integration tests
# ============================================================================

class TestIntegration:
    """Integration tests for multiple signatures."""

    def test_all_non_climate_signatures(self, synthetic_streamflow):
        """Test all non-climate signatures can be computed."""
        results = {}

        results.update(calculate_flow_vols_by_year(synthetic_streamflow))
        results.update(analyze_flashiness_trends(synthetic_streamflow))
        results.update(analyze_flow_timing_trends(synthetic_streamflow))
        results.update(analyze_fdc_trends(synthetic_streamflow))
        results.update(analyze_baseflow_indices(synthetic_streamflow))
        results.update(analyze_recession_parameters(synthetic_streamflow))
        results.update(calculate_pulse_metrics(synthetic_streamflow))

        # Should have many metrics
        assert len(results) > 400, f"Expected 400+ metrics, got {len(results)}"

        # Count non-NA values
        non_na = sum(1 for v in results.values() if not np.isnan(v))
        assert non_na > 300, f"Expected 300+ non-NA metrics, got {non_na}"

    def test_empty_data_handling(self):
        """Test signatures handle empty data gracefully."""
        empty_df = pd.DataFrame({
            "gage_id": [],
            "date": [],
            "Q": [],
            "water_year": [],
            "month": [],
            "dowy": []
        })

        # Should not raise exceptions
        result = calculate_flow_vols_by_year(empty_df)
        assert all(np.isnan(v) for v in result.values())

    def test_short_data_handling(self):
        """Test signatures handle short time series."""
        dates = pd.date_range("2000-10-01", periods=100, freq="D")
        short_df = pd.DataFrame({
            "gage_id": "SHORT",
            "date": dates,
            "Q": np.random.exponential(1, 100)
        })
        short_df = add_water_year_columns(short_df, date_col="date")

        # Should not raise exceptions (but may return NaN)
        result = calculate_flow_vols_by_year(short_df)
        assert isinstance(result, dict)

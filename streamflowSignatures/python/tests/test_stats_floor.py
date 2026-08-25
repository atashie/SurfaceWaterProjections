"""
Tests for the stats floor (min_values_for_stats) and force_skip_trends —
Phase 1b of the port campaign, porting julia/test/test_stats_floor.jl's
generate_stats-level scope plus direct force_skip_trends mechanism tests
(the plan's Codex review, F10: the mechanism must be tested when it lands,
not when its first consumer — snow — arrives in Phase 4).

Contract (mirrors julia/src/stats.jl):
- A metric with fewer non-NaN annual values than the floor emits NaN for ALL
  8 statistics AND its changepoint fields (same branch as min_rows).
- force_skip_trends suppresses the 6 trend stats through the SAME path as
  trend_completeness: mean/median still computed, changepoint unaffected.
- Recession and elasticity are exempt at the ORCHESTRATION layer: their
  function signatures do not even accept the kwarg (as in Julia).

Run with: pytest tests/test_stats_floor.py -v
"""

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures import generate_stats, STAT_SUFFIXES, CP_SUFFIXES
from streamflow_signatures.recession import analyze_recession_parameters
from streamflow_signatures.elasticity import calculate_streamflow_elasticity

CP_CONFIG = {"start_year": 1900, "end_year": 2100,
             "min_total_obs": 20, "min_segment_obs": 10}


def _series_df(n, start=1990, name="metric"):
    return pd.DataFrame({
        "water_year": list(range(start, start + n)),
        name: [10.0 + 0.5 * i + (i % 3) * 0.1 for i in range(n)],
    })


class TestStatsFloor:
    def test_below_floor_all_nan_including_cp(self):
        df = _series_df(19)
        result = generate_stats(df, value_cols=["metric"],
                                min_values_for_stats=20, changepoint=CP_CONFIG)
        for suffix in STAT_SUFFIXES + CP_SUFFIXES:
            assert np.isnan(result[f"metric{suffix}"]), f"metric{suffix} not NaN"

    def test_at_floor_stats_computed(self):
        df = _series_df(20)
        result = generate_stats(df, value_cols=["metric"],
                                min_values_for_stats=20, changepoint=CP_CONFIG)
        assert not np.isnan(result["metric_mean"])
        assert not np.isnan(result["metric_senn_slp"])
        assert not np.isnan(result["metric_pettitt_cp_year"])

    def test_no_floor_preserves_legacy_behavior(self):
        df = _series_df(5)
        with_none = generate_stats(df, value_cols=["metric"],
                                   min_values_for_stats=None)
        legacy = generate_stats(df, value_cols=["metric"])
        assert with_none == legacy
        assert not np.isnan(with_none["metric_mean"])

    def test_floor_counts_non_nan_values_only(self):
        """25 rows but only 19 non-NaN -> below a floor of 20."""
        df = _series_df(25)
        df.loc[df.index[:6], "metric"] = np.nan
        result = generate_stats(df, value_cols=["metric"],
                                min_values_for_stats=20)
        for suffix in STAT_SUFFIXES:
            assert np.isnan(result[f"metric{suffix}"])

    def test_clustered_series_regression(self):
        """The motivating 07292500 shape: 4 clustered non-NaN years in a long
        frame pass the own-span trend gate (4/4 = 100%) and the min_rows floor
        (4 >= 3) — only the stats floor stops the 4-point Theil-Sen slope."""
        years = list(range(1980, 2020))
        vals = [np.nan] * 40
        for i, v in zip(range(2, 6), [5.0, 6.0, 4.5, 7.0]):
            vals[i] = v
        df = pd.DataFrame({"water_year": years, "metric": vals})
        without_floor = generate_stats(df, value_cols=["metric"],
                                       trend_completeness=0.8,
                                       decade_completeness=0.8)
        assert not np.isnan(without_floor["metric_senn_slp"])  # the bug shape
        with_floor = generate_stats(df, value_cols=["metric"],
                                    trend_completeness=0.8,
                                    decade_completeness=0.8,
                                    min_values_for_stats=20)
        for suffix in STAT_SUFFIXES:
            assert np.isnan(with_floor[f"metric{suffix}"])

    def test_per_metric_independence(self):
        """One metric below the floor, one above — gated independently."""
        n = 25
        df = _series_df(n, name="dense")
        df["sparse"] = [df["dense"][i] if i < 10 else np.nan for i in range(n)]
        result = generate_stats(df, value_cols=["dense", "sparse"],
                                min_values_for_stats=20)
        assert not np.isnan(result["dense_mean"])
        assert np.isnan(result["sparse_mean"])

    def test_recession_and_elasticity_do_not_accept_floor(self):
        """Exemption is structural (matching Julia): the exempt functions'
        signatures do not accept the kwarg, so the orchestrator cannot leak it."""
        df = pd.DataFrame({
            "water_year": [2000] * 30, "Q": np.linspace(10, 1, 30),
            "dowy": list(range(1, 31)), "month": [10] * 30,
            "date": pd.date_range("1999-10-01", periods=30),
        })
        with pytest.raises(TypeError):
            analyze_recession_parameters(df, min_values_for_stats=20)
        with pytest.raises(TypeError):
            calculate_streamflow_elasticity(df, min_values_for_stats=20)


class TestForceSkipTrends:
    def test_gated_metric_trends_nan_mean_survives(self):
        df = _series_df(25)
        result = generate_stats(df, value_cols=["metric"],
                                force_skip_trends={"metric"},
                                changepoint=CP_CONFIG)
        for suffix in ["_senn_slp", "_linear_slp", "_spearman_rho",
                       "_spearman_pval", "_mk_rho", "_mk_pval"]:
            assert np.isnan(result[f"metric{suffix}"])
        assert not np.isnan(result["metric_mean"])
        assert not np.isnan(result["metric_median"])
        # Changepoint unaffected by the forced skip
        assert not np.isnan(result["metric_pettitt_cp_year"])

    def test_ungated_metric_unaffected(self):
        n = 25
        df = _series_df(n, name="gated")
        df["free"] = df["gated"] + 1.0
        result = generate_stats(df, value_cols=["gated", "free"],
                                force_skip_trends={"gated"})
        assert np.isnan(result["gated_senn_slp"])
        assert not np.isnan(result["free_senn_slp"])

    def test_empty_set_is_inert(self):
        df = _series_df(25)
        gated = generate_stats(df, value_cols=["metric"], force_skip_trends=set())
        plain = generate_stats(df, value_cols=["metric"])
        assert gated == plain

    def test_mean_median_match_ungated_values(self):
        """The forced-skip path must compute the same mean/median as ungated."""
        df = _series_df(25)
        gated = generate_stats(df, value_cols=["metric"],
                               force_skip_trends={"metric"})
        plain = generate_stats(df, value_cols=["metric"])
        assert gated["metric_mean"] == plain["metric_mean"]
        assert gated["metric_median"] == plain["metric_median"]

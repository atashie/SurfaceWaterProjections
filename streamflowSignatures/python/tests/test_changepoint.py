"""
Tests for the Pettitt changepoint port (Phase 1a of the port campaign).

Port of julia/test/test_changepoint.jl's pipeline-reachable scope (the BIC
model-selection testsets cover exported Julia-only API excluded from parity
scope), PLUS the deterministic exact-equality fixtures the plan's Codex review
required (F12): tied maxima, unsorted years, all-constant series, boundary
segment sizes, NaN handling. Every frozen expectation below was cross-checked
against the canonical Julia implementation on 2026-08-24 (identical cp_year on
all fixtures incl. the tie; p-values bit-identical; means within 2 ulp).

Run with: pytest tests/test_changepoint.py -v
"""

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures import generate_stats, CP_SUFFIXES
from streamflow_signatures.changepoint import pettitt_test, segment_differential_metrics

CP_CONFIG = {"start_year": 1900, "end_year": 2100,
             "min_total_obs": 20, "min_segment_obs": 10}


# ---------------------------------------------------------------------------
# Deterministic fixtures — exact equality, cross-checked against Julia
# ---------------------------------------------------------------------------

class TestPettittFixtures:
    def test_clean_step(self):
        """10 low + 15 high (distinct within-block values): cp at the boundary."""
        yrs = np.arange(1990, 2015).astype(float)
        vals = np.concatenate([10 + 0.01 * np.arange(10), 20 + 0.01 * np.arange(15)])
        r = pettitt_test(yrs, vals)
        assert r["cp_year"] == 2001.0  # frozen, Julia-cross-checked (t=12, 1-based)
        assert r["pval"] == pytest.approx(0.00025039952406474066, rel=1e-14)
        assert r["pre_mean"] == pytest.approx(11.705, rel=1e-12)
        assert r["post_mean"] == pytest.approx(20.08, rel=1e-12)
        d = segment_differential_metrics(yrs, vals, r["cp_year"])
        assert d["delta_mean"] == pytest.approx(8.375, rel=1e-12)
        assert d["pct_change"] == pytest.approx(71.5506193934216, rel=1e-12)

    def test_tied_maximum_takes_first(self):
        """Antisymmetric 'tent' series: |U| attains its max (30) at TWO candidate
        indices — the FIRST must win (Julia's strict `>` update rule)."""
        yrs = np.arange(1980, 2003).astype(float)  # n = 23
        vals = np.concatenate([np.arange(1, 13), np.arange(11, 0, -1)]).astype(float)
        # Prove the tie is real, then pin the first-max choice
        U = np.cumsum(np.sign(vals[:, None] - vals[None, :]).sum(axis=1))
        seg = np.abs(U[9:13])
        assert seg.tolist() == [30.0, 11.0, 11.0, 30.0]  # tie between t=10 and t=13
        r = pettitt_test(yrs, vals)
        assert r["cp_year"] == 1989.0  # t=10 (first max), NOT 1992 (t=13)
        assert r["pval"] == 1.0

    def test_constant_series_nonzero(self):
        yrs = np.arange(2000, 2025).astype(float)
        r = pettitt_test(yrs, np.full(25, 5.0))
        assert r["cp_year"] == 2009.0  # all |U|=0: first candidate (t=min_seg)
        assert r["pval"] == 1.0
        assert r["pre_mean"] == 5.0 and r["post_mean"] == 5.0
        d = segment_differential_metrics(yrs, np.full(25, 5.0), r["cp_year"])
        assert d["delta_mean"] == 0.0
        assert d["pct_change"] == 0.0

    def test_constant_series_zero_pct_nan(self):
        """pre_mean == 0 -> pct_change NaN (|pre_mean| <= 1e-10 guard)."""
        yrs = np.arange(2000, 2025).astype(float)
        vals = np.zeros(25)
        r = pettitt_test(yrs, vals)
        assert r["cp_year"] == 2009.0
        d = segment_differential_metrics(yrs, vals, r["cp_year"])
        assert d["delta_mean"] == 0.0
        assert np.isnan(d["pct_change"])

    def test_boundary_n20_single_candidate(self):
        """n == min_total_obs == 2*min_segment_obs: exactly one candidate (t=10)."""
        yrs = np.arange(1991, 2011).astype(float)
        vals = np.concatenate([np.full(10, 1.0) + 0.001 * np.arange(10),
                               np.full(10, 3.0) + 0.001 * np.arange(10)])
        r = pettitt_test(yrs, vals)
        assert r["cp_year"] == 2000.0
        assert r["pval"] == pytest.approx(0.0015809806462399323, rel=1e-14)

    def test_below_min_total_obs(self):
        yrs = np.arange(1991, 2010).astype(float)  # n = 19 < 20
        r = pettitt_test(yrs, np.arange(19).astype(float))
        assert all(np.isnan(v) for v in r.values())

    def test_nan_pairs_filtered(self):
        """NaNs removed before the test; result matches the frozen expectation."""
        yrs = np.arange(1985, 2013).astype(float)
        vals = np.concatenate([5 + 0.01 * np.arange(12), 15 + 0.01 * np.arange(16)])
        vals[[3, 7, 20, 25]] = np.nan
        r = pettitt_test(yrs, vals)
        assert r["cp_year"] == 1998.0
        assert r["pval"] == pytest.approx(0.0003537738044851332, rel=1e-14)

    def test_unsorted_input_invariance(self):
        """Shuffled rows must give the identical result to sorted input."""
        yrs = np.arange(1990, 2015).astype(float)
        vals = np.concatenate([10 + 0.01 * np.arange(10), 20 + 0.01 * np.arange(15)])
        rng = np.random.default_rng(7)
        perm = rng.permutation(len(yrs))
        r_sorted = pettitt_test(yrs, vals)
        r_shuffled = pettitt_test(yrs[perm], vals[perm])
        assert r_sorted == r_shuffled

    def test_pval_formula(self):
        """p = min(1, 2 exp(-6 K^2 / (n^3 + n^2))) — independent recomputation."""
        yrs = np.arange(1990, 2015).astype(float)
        vals = np.concatenate([10 + 0.01 * np.arange(10), 20 + 0.01 * np.arange(15)])
        U = np.cumsum(np.sign(vals[:, None] - vals[None, :]).sum(axis=1))
        n = len(vals)
        K = np.abs(U[9:n - 10]).max()
        expected = min(1.0, 2.0 * np.exp(-6.0 * K ** 2 / (n ** 3 + n ** 2)))
        assert pettitt_test(yrs, vals)["pval"] == pytest.approx(expected, rel=1e-15)

    def test_segment_min_three_values(self):
        """A cp_year splitting off <3 values on a side -> NaN differentials."""
        yrs = np.arange(2000, 2025).astype(float)
        vals = np.arange(25).astype(float)
        d = segment_differential_metrics(yrs, vals, 2001.0)  # 2 pre values
        assert all(np.isnan(v) for v in d.values())
        d = segment_differential_metrics(yrs, vals, np.nan)
        assert all(np.isnan(v) for v in d.values())


# ---------------------------------------------------------------------------
# generate_stats() integration — ports the Julia integration testset
# ---------------------------------------------------------------------------

def _step_df(n_pre=15, n_post=15, start=1985, lo=10.0, hi=20.0):
    years = list(range(start, start + n_pre + n_post))
    vals = ([lo + 0.01 * i for i in range(n_pre)]
            + [hi + 0.01 * i for i in range(n_post)])
    return pd.DataFrame({"water_year": years, "metric": vals})


class TestGenerateStatsIntegration:
    def test_pettitt_columns_appended(self):
        df = _step_df()
        result = generate_stats(df, value_cols=["metric"], changepoint=CP_CONFIG)
        for suffix in CP_SUFFIXES:
            assert f"metric{suffix}" in result
        assert not np.isnan(result["metric_pettitt_cp_year"])
        assert result["metric_pettitt_pval"] < 0.05
        assert result["metric_pettitt_delta_mean"] == pytest.approx(10.0, abs=2.0)

    def test_no_changepoint_no_columns(self):
        df = _step_df()
        result = generate_stats(df, value_cols=["metric"], changepoint=None)
        assert not any(k.startswith("metric_pettitt") for k in result)

    def test_independent_of_trend_gating(self):
        """A trend-gated metric (sparse series) still gets its Pettitt fields."""
        years = list(range(1980, 2020))
        vals = [float(i) if i % 2 == 0 else np.nan for i in range(40)]
        df = pd.DataFrame({"water_year": years, "metric": vals})
        result = generate_stats(df, value_cols=["metric"],
                                trend_completeness=0.9, decade_completeness=0.9,
                                changepoint=CP_CONFIG)
        assert np.isnan(result["metric_senn_slp"])         # trend gated
        assert not np.isnan(result["metric_mean"])          # mean survives
        assert not np.isnan(result["metric_pettitt_cp_year"])  # CP independent

    def test_below_min_rows_emits_nan_cp_keys(self):
        df = pd.DataFrame({"water_year": [2000, 2001], "metric": [1.0, 2.0]})
        result = generate_stats(df, value_cols=["metric"], min_rows=3,
                                changepoint=CP_CONFIG)
        for suffix in CP_SUFFIXES:
            assert np.isnan(result[f"metric{suffix}"])

    def test_window_filter(self):
        """Years outside [start_year, end_year] are excluded from the CP series;
        a step living entirely outside the window disappears."""
        # Step at 1975 (outside window starting 1980); inside-window data has a
        # deterministic repeating pattern (no trend, no step)
        years = list(range(1960, 2010))
        vals = [5.0 + 0.01 * (i % 3) if y < 1975 else 6.0 + 0.01 * (i % 3)
                for i, y in enumerate(years)]
        df = pd.DataFrame({"water_year": years, "metric": vals})
        cp = dict(CP_CONFIG, start_year=1980, end_year=2009)
        result = generate_stats(df, value_cols=["metric"], changepoint=cp)
        # Only the post-1980 (flat-ish trend) series is analyzed: no significant step
        assert result["metric_pettitt_pval"] > 0.05
        assert result["metric_pettitt_cp_year"] >= 1980

    def test_window_too_short_nan(self):
        """Window trims the series below min_total_obs -> NaN CP fields."""
        df = _step_df(n_pre=15, n_post=15, start=1985)
        cp = dict(CP_CONFIG, start_year=2000, end_year=2005)
        result = generate_stats(df, value_cols=["metric"], changepoint=cp)
        assert np.isnan(result["metric_pettitt_cp_year"])
        assert np.isnan(result["metric_pettitt_pval"])

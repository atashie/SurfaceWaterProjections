"""
Tests for the 14 snow metrics (Phase 4 of the port campaign; port of
julia/test/test_snow_metrics.jl + test_snow_record_decade_gate.jl scope).

Every expectation is either hand-derived from the definitions or was
cross-checked against canonical Julia on 2026-08-25 (126 values across the
three fixtures below: 0 mismatches, worst relative difference 3.5e-16).

Run with: pytest tests/test_snow_metrics.py -v
"""

from datetime import date

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures.snow import (
    calculate_snow_metrics, snow_spells, SNOW_METRICS,
)
from streamflow_signatures import calculate_all_signatures, AnnualCollector

CP_CONFIG = {"start_year": 1900, "end_year": 2100,
             "min_total_obs": 20, "min_segment_obs": 10}


def _year_rows(wy, swe_of_dowy, ppt=2.0):
    d0, d1 = date(wy - 1, 10, 1), date(wy, 9, 30)
    nd = (d1 - d0).days + 1
    return [(wy, i, float(swe_of_dowy(i, wy)), ppt) for i in range(1, nd + 1)]


def make_gage(swe_of_dowy, n_years=25, start=1996, ppt=2.0):
    rows = []
    for wy in range(start, start + n_years):
        rows += _year_rows(wy, swe_of_dowy, ppt)
    return pd.DataFrame(rows, columns=["water_year", "dowy", "SWE", "PPT"])


def triangle(i, wy):
    if i <= 90:
        return 0.0
    if i <= 180:
        return (i - 90) * 1.5
    return max(0.0, 135.0 - (i - 180) * 1.7) if i <= 260 else 0.0


def thaw(i, wy):
    if i <= 60:
        return 0.0
    if i <= 110:
        return (i - 60) * 2.0
    if i <= 130:
        return 0.0          # mid-winter thaw -> two spells
    if i <= 200:
        return (i - 130) * 2.5
    if i <= 270:
        return max(0.0, 175.0 - (i - 200) * 2.6)
    return 0.0


def vanishing(i, wy, last_snowy=2015):
    """Snowpack present every year through `last_snowy`, then gone. Over
    1996-2020 that is 20 snowy years (>= the 20-value stats floor, so the floor
    does NOT fire) while the record's LAST decade (2011-2020) is only 50%
    snowy — below the 80% decade requirement, so the record-anchored gate is
    what suppresses the trends. This isolates the gate from the floor."""
    if wy <= last_snowy:
        return triangle(i, wy)
    return 0.0 if i > 250 else 3.0


class TestSpellHelper:
    def test_maximal_runs(self):
        m = np.array([False, True, True, False, True, False, False, True])
        sp = snow_spells(m)
        assert [list(s) for s in sp] == [[1, 2], [4], [7]]

    def test_edges(self):
        assert snow_spells(np.array([True, True])) == [range(0, 2)]
        assert snow_spells(np.zeros(5, bool)) == []


class TestTriangleGage:
    """One clean snowpack per year: every metric is hand-derivable."""

    @pytest.fixture(scope="class")
    def r(self):
        return calculate_snow_metrics(make_gage(triangle), trend_completeness=0.6,
                                      decade_completeness=0.8, min_values_for_stats=20)

    def test_peak_magnitude_and_timing(self, r):
        # SWE rises 1.5/day from dowy 91; peak at dowy 180 = 90*1.5 = 135
        assert r["swe_max_mean"] == pytest.approx(135.0)
        assert r["swe_max_dowy_mean"] == pytest.approx(180.0)

    def test_snow_on_off_and_cover(self, r):
        # First day at/above the 10 mm threshold: 1.5*k >= 10 -> k = 7 -> dowy 97
        assert r["snow_on_dowy_mean"] == pytest.approx(97.0)
        # Melt-out: 135 - 1.7*k >= 10 -> k <= 73.5 -> k = 73 -> last snow day
        # 253, so snow_off = 254 (first snow-free day after the anchor spell)
        assert r["snow_off_dowy_mean"] == pytest.approx(254.0)
        assert r["snow_cover_days_mean"] == pytest.approx(254.0 - 97.0)

    def test_melt_season_and_rate(self, r):
        assert r["melt_season_days_mean"] == pytest.approx(254.0 - 180.0)
        assert r["melt_rate_mean"] == pytest.approx(135.0 / 74.0)

    def test_ssm_fully_seasonal(self, r):
        # A single spell of 158 days >= 60 -> fully seasonal
        assert r["ssm_mean"] == pytest.approx(1.0)

    def test_no_melt_before_peak(self, r):
        # Monotone rise to the peak: nothing melts before it
        assert r["melt_before_peak_mean"] == pytest.approx(0.0)
        assert r["melt_before_peak_pct_mean"] == pytest.approx(0.0)
        assert r["melt_before_peak_to_max_swe_mean"] == pytest.approx(0.0)

    def test_swe_apr1_leap_safe(self, r):
        # April 1 is dowy 183 (non-leap) / 184 (leap): SWE = 135 - 1.7*(dowy-180)
        # -> 129.9 or 128.2; the mean over 25 years sits between
        assert 128.0 < r["swe_apr1_mean"] < 130.0

    def test_swe_max_to_ppt(self, r):
        # PPT = 2 mm/day -> annual total 730/732; ratio = 135 / total
        assert r["swe_max_to_ppt_mean"] == pytest.approx(135.0 / 730.5, rel=2e-3)

    def test_all_metrics_present(self, r):
        for m in SNOW_METRICS:
            assert f"{m}_mean" in r and f"{m}_median" in r


class TestThawGage:
    """Mid-winter thaw: two spells, so the anchor spell and SSM both matter."""

    @pytest.fixture(scope="class")
    def r(self):
        return calculate_snow_metrics(make_gage(thaw), trend_completeness=0.6,
                                      decade_completeness=0.8, min_values_for_stats=20)

    def test_anchor_is_the_peak_spell(self, r):
        # Second spell peaks at dowy 200 (175 mm) > first spell's 100 mm
        assert r["swe_max_mean"] == pytest.approx(175.0)
        assert r["swe_max_dowy_mean"] == pytest.approx(200.0)
        # snow_on belongs to the SECOND spell (the anchor), not the first
        assert r["snow_on_dowy_mean"] > 130

    def test_melt_before_peak_counts_the_thaw(self, r):
        # The first spell fully melts (100 mm) before the annual peak
        assert r["melt_before_peak_mean"] == pytest.approx(100.0)
        assert r["melt_before_peak_to_max_swe_mean"] == pytest.approx(100.0 / 175.0)
        assert 0.0 < r["melt_before_peak_pct_mean"] < 100.0

    def test_ssm_pools_all_spells(self, r):
        # Spell 1 ~ 46 days (< 60, ephemeral), spell 2 ~ 130 days (seasonal)
        assert -1.0 <= r["ssm_mean"] <= 1.0
        assert r["ssm_mean"] < 1.0        # not fully seasonal
        assert r["ssm_mean"] > 0.0        # seasonal days dominate


class TestSnowFreeAndZeroPolicy:
    def test_magnitudes_zero_timing_nan(self):
        """Operationally snow-free years: magnitude metrics emit valid ZEROS,
        timing/melt/regime metrics emit NaN."""
        r = calculate_snow_metrics(make_gage(lambda i, wy: 5.0),  # always < 10 mm
                                   trend_completeness=0.6, decade_completeness=0.8)
        assert r["swe_max_mean"] == 0.0
        assert r["snow_cover_days_mean"] == 0.0
        assert r["swe_apr1_mean"] == 0.0
        assert r["swe_max_to_ppt_mean"] == 0.0
        for m in ("swe_max_dowy", "snow_on_dowy", "snow_off_dowy",
                  "melt_season_days", "melt_rate", "ssm", "melt_com_dowy"):
            assert np.isnan(r[f"{m}_mean"]), m

    def test_boundary_censoring(self):
        """An anchor spell touching Oct 1 / Sep 30 censors snow_on / snow_off."""
        r = calculate_snow_metrics(make_gage(lambda i, wy: 100.0),  # snow all year
                                   trend_completeness=0.6, decade_completeness=0.8)
        assert np.isnan(r["snow_on_dowy_mean"])
        assert np.isnan(r["snow_off_dowy_mean"])
        assert np.isnan(r["melt_season_days_mean"])
        assert r["snow_cover_days_mean"] in (365.0, pytest.approx(365.25, abs=0.5))


class TestGridGuard:
    def test_incomplete_year_left_nan(self):
        """A year that is not a complete 365/366-day grid is skipped entirely."""
        df = make_gage(triangle, n_years=25)
        df = df[~((df.water_year == 2000) & (df.dowy == 50))]     # punch a hole
        r = calculate_snow_metrics(df, trend_completeness=0.6, decade_completeness=0.8)
        c = AnnualCollector()
        calculate_snow_metrics(df, collector=c)
        rows = [(s, y, v) for s, y, v in zip(c.signature, c.water_year, c.value)
                if s == "swe_max" and y == 2000]
        assert len(rows) == 1 and np.isnan(rows[0][2])

    def test_nan_swe_year_left_nan(self):
        df = make_gage(triangle, n_years=25)
        df.loc[(df.water_year == 2001) & (df.dowy == 100), "SWE"] = np.nan
        c = AnnualCollector()
        calculate_snow_metrics(df, collector=c)
        rows = [v for s, y, v in zip(c.signature, c.water_year, c.value)
                if s == "swe_max" and y == 2001]
        assert len(rows) == 1 and np.isnan(rows[0])


class TestRecordAnchoredDecadeGate:
    def test_alternating_years_suppress_trends(self):
        """Snow vanishing mid-record -> the gated metrics fail the 80% LAST
        decade requirement and lose their 6 trend stats, keeping mean/median.
        20 snowy years keeps the stats floor from firing, so this isolates the
        record-anchored gate."""
        r = calculate_snow_metrics(make_gage(vanishing), trend_completeness=0.6,
                                   decade_completeness=0.8, min_values_for_stats=20,
                                   changepoint=CP_CONFIG)
        for m in ("swe_max_dowy", "ssm", "melt_com_dowy"):
            assert np.isnan(r[f"{m}_senn_slp"]), m
            assert not np.isnan(r[f"{m}_mean"]), m      # mean/median survive
        # magnitude metrics are EXEMPT: dense zero-including series keep trends
        assert not np.isnan(r["swe_max_mean"])
        assert "swe_max_pettitt_cp_year" in r           # changepoint unaffected

    def test_gate_inert_on_consistent_snow(self):
        r = calculate_snow_metrics(make_gage(triangle), trend_completeness=0.6,
                                   decade_completeness=0.8, min_values_for_stats=20)
        # every year snowy -> no suppression (constant series still give NaN
        # rank stats, so assert on the slope, which is defined)
        assert not np.isnan(r["swe_max_dowy_senn_slp"])


class TestOrchestratorIntegration:
    def test_requires_explicit_snow_data(self):
        """Snow metrics never run off an SWE column in the main gage frame."""
        gage = make_gage(triangle, n_years=25)
        gage["Q"] = 1.0
        gage["date"] = pd.NaT
        out = calculate_all_signatures(gage, has_climate=False, gage_id="T")
        assert not any(k.startswith("swe_max") for k in out)

        out2 = calculate_all_signatures(gage, has_climate=False, snow_data=gage,
                                        trend_completeness=0.6, decade_completeness=0.8,
                                        gage_id="T")
        assert "swe_max_mean" in out2 and not np.isnan(out2["swe_max_mean"])

    def test_zero_row_snow_data_emits_full_key_set(self):
        """A gage with SWE but no valid SWE year -> all 14 metrics present, NaN."""
        empty = pd.DataFrame({"water_year": pd.Series([], dtype="int64"),
                              "dowy": pd.Series([], dtype="int64"),
                              "SWE": pd.Series([], dtype="float64")})
        r = calculate_snow_metrics(empty, changepoint=CP_CONFIG)
        for m in SNOW_METRICS:
            assert np.isnan(r[f"{m}_mean"]), m
            assert f"{m}_pettitt_cp_year" in r, m

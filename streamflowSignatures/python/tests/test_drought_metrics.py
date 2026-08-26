"""
Tests for the 10 streamflow drought metrics + 5 threshold scalars (Phase 5 of
the port campaign; port of julia/test/test_drought_metrics.jl).

Cross-checked against canonical Julia on 2026-08-25: 105 values across three
fixtures (seasonal low-flow, intermittent/zero-threshold, and a gapped record)
— 0 mismatches, worst relative difference 6.0e-15.

Includes the boundary-attribution fixture Codex required of the Julia original
(MAJOR-1): a conservation check alone would pass under a uniform water-year
shift, so a low-flow block straddling Sep 30 -> Oct 1 is split explicitly.

Run with: pytest tests/test_drought_metrics.py -v
"""

from datetime import date, timedelta

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures.drought import (
    calculate_drought_metrics, weibull_quantile, smooth_daily_flow,
    DROUGHT_METRICS, DROUGHT_DURATION_METRICS, DROUGHT_DEFICIT_METRICS,
    DROUGHT_THRESHOLD_SCALARS,
)
from streamflow_signatures import AnnualCollector

CP_CONFIG = {"start_year": 1900, "end_year": 2100,
             "min_total_obs": 20, "min_segment_obs": 10}


def make_gage(q_of, n_years=25, start=1996, skip=None):
    rows = []
    for wy in range(start, start + n_years):
        d0, d1 = date(wy - 1, 10, 1), date(wy, 9, 30)
        nd = (d1 - d0).days + 1
        for i in range(1, nd + 1):
            if skip is not None and skip(i, wy):
                continue
            rows.append((pd.Timestamp(d0 + timedelta(days=i - 1)), wy, i, float(q_of(i, wy))))
    return pd.DataFrame(rows, columns=["date", "water_year", "dowy", "Q"])


class TestWeibullQuantile:
    def test_hand_derived_interpolation(self):
        """p_i = i/(n+1): for n=9, p=0.5 -> h=5 -> exactly the 5th value."""
        x = np.arange(1.0, 10.0)                    # 1..9, n = 9
        assert weibull_quantile(x, 0.5) == pytest.approx(5.0)
        # p = 0.25 -> h = 2.5 -> midpoint of the 2nd and 3rd values
        assert weibull_quantile(x, 0.25) == pytest.approx(2.5)

    def test_below_plotting_range_is_na_not_clamped(self):
        """Below 1/(n+1) the sample cannot resolve the probability: NaN, NOT the
        minimum (a project policy layered on Hyndman-Fan type 6)."""
        x = np.arange(1.0, 10.0)                    # smallest position = 0.1
        assert np.isnan(weibull_quantile(x, 0.05))
        assert weibull_quantile(x, 0.05, below_range="clamp") == pytest.approx(1.0)

    def test_above_range_returns_max(self):
        x = np.arange(1.0, 10.0)
        assert weibull_quantile(x, 0.99) == pytest.approx(9.0)

    def test_nans_dropped_and_empty(self):
        x = np.array([1.0, np.nan, 2.0, 3.0, np.nan])
        assert weibull_quantile(x, 0.5) == pytest.approx(2.0)
        assert np.isnan(weibull_quantile(np.array([np.nan]), 0.5))


class TestSmoothing:
    def _dates(self, n, start="2000-01-01"):
        return pd.date_range(start, periods=n, freq="D").to_numpy()

    def test_centered_window_and_edge_shrink(self):
        q = np.arange(1.0, 11.0)
        out = smooth_daily_flow(self._dates(10), q, window=7, alignment="center",
                                min_valid=1)
        # interior day 4 (0-based 3): mean of days 1..7 = 4
        assert out[3] == pytest.approx(4.0)
        # first day: window shrinks to days 1..4 -> mean 2.5
        assert out[0] == pytest.approx(2.5)

    def test_min_valid_yields_nan(self):
        q = np.arange(1.0, 11.0)
        out = smooth_daily_flow(self._dates(10), q, window=7, alignment="center",
                                min_valid=4)
        assert not np.isnan(out[3])
        # a 2-day run cannot reach min_valid=4
        short = smooth_daily_flow(self._dates(2), np.array([1.0, 2.0]), window=7,
                                  alignment="center", min_valid=4)
        assert np.isnan(short).all()

    def test_never_averages_across_a_gap(self):
        """Two runs separated by a missing day must not blend."""
        dates = np.array([np.datetime64("2000-01-01") + np.timedelta64(d, "D")
                          for d in [0, 1, 2, 3, 4, 6, 7, 8, 9, 10]])
        q = np.array([1.0] * 5 + [100.0] * 5)
        out = smooth_daily_flow(dates, q, window=7, alignment="center", min_valid=1)
        assert out[:5] == pytest.approx([1.0] * 5)      # untouched by the 100s
        assert out[5:] == pytest.approx([100.0] * 5)

    def test_duplicate_date_breaks_the_run(self):
        dates = np.array([np.datetime64("2000-01-01") + np.timedelta64(d, "D")
                          for d in [0, 1, 1, 2]])
        out = smooth_daily_flow(dates, np.array([1.0, 2.0, 3.0, 4.0]), window=3,
                                alignment="center", min_valid=1)
        assert np.isfinite(out).all()                   # runs are just shorter
        assert out[0] == pytest.approx(1.5)             # run = [1.0, 2.0]

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            smooth_daily_flow(self._dates(3), np.array([1.0, 2.0]))


class TestSeasonalGage:
    @pytest.fixture(scope="class")
    def r(self):
        gage = make_gage(lambda i, wy: 2.0 + 1.8*np.sin(2*np.pi*(i-30)/365)
                         + 0.3*((i*7919) % 11)/11)
        return calculate_drought_metrics(gage, trend_completeness=0.6,
                                         decade_completeness=0.8,
                                         min_values_for_stats=20)

    def test_all_metrics_and_scalars_present(self, r):
        for m in DROUGHT_METRICS:
            assert f"{m}_mean" in r
        for s in DROUGHT_THRESHOLD_SCALARS:
            assert s in r and np.isfinite(r[s])

    def test_level_monotonicity(self, r):
        """Duration and deficit are non-decreasing in the severity level by
        construction (a higher percentile threshold catches more days)."""
        durs = [r[f"{m}_mean"] for m in DROUGHT_DURATION_METRICS]
        defs = [r[f"{m}_mean"] for m in DROUGHT_DEFICIT_METRICS]
        assert durs == sorted(durs)
        assert defs == sorted(defs)

    def test_thresholds_ordered(self, r):
        thr = [r[s] for s in DROUGHT_THRESHOLD_SCALARS]
        assert thr == sorted(thr)

    def test_weak_construction_anchor(self, r):
        """A level-p threshold sits below ~p% of all days by construction, so
        duration_p10 averages ~36.5 d/yr. Near-circular — a sanity anchor only,
        as documented."""
        assert 30.0 < r["drought_duration_fixed_p10_mean"] < 43.0


class TestIntermittentGage:
    def test_zero_threshold_strict_comparison(self):
        """Where most days are zero flow the low thresholds are exactly 0, and
        the STRICT `<` then reports 0 duration/deficit — the documented
        interpretation trap, pinned here."""
        gage = make_gage(lambda i, wy: 0.0 if 150 <= i <= 300
                         else 1.0 + ((i*104729) % 7)/7)
        r = calculate_drought_metrics(gage, trend_completeness=0.6,
                                      decade_completeness=0.8)
        assert r["drought_threshold_fixed_p2"] == pytest.approx(0.0)
        assert r["drought_duration_fixed_p2_mean"] == pytest.approx(0.0)
        assert r["drought_deficit_fixed_p2_mean"] == pytest.approx(0.0)


class TestBoundaryAttribution:
    """Codex MAJOR-1 on the Julia original: a conservation check alone would
    pass under a UNIFORM water-year shift, so attribution must be pinned
    explicitly with a low-flow block straddling Sep 30 -> Oct 1."""

    def test_block_split_across_the_water_year_boundary(self):
        # High baseline everywhere except a 13-day low block spanning
        # 2005-09-25..2005-10-07 — i.e. 6 days in WY2005 and 7 days in WY2006.
        lo_start, lo_end = date(2005, 9, 25), date(2005, 10, 7)

        def q_of(i, wy):
            d = date(wy - 1, 10, 1) + timedelta(days=i - 1)
            return 0.0 if lo_start <= d <= lo_end else 10.0

        gage = make_gage(q_of, n_years=25)
        c = AnnualCollector()
        calculate_drought_metrics(gage, collector=c)
        dur = {y: v for s, y, v in zip(c.signature, c.water_year, c.value)
               if s == "drought_duration_fixed_p30"}
        # The block must be SPLIT, not attributed wholly to either year
        assert dur[2005] > 0 and dur[2006] > 0
        assert dur[2005] + dur[2006] >= 13          # smoothing widens the dip
        # ...and no other year sees drought (uniform-shift guard)
        others = [v for y, v in dur.items() if y not in (2005, 2006)]
        assert all(v == 0 for v in others)


class TestSchemaContract:
    def test_empty_frame_emits_full_key_set(self):
        empty = pd.DataFrame({"date": pd.Series([], dtype="datetime64[ns]"),
                              "water_year": pd.Series([], dtype="int64"),
                              "Q": pd.Series([], dtype="float64")})
        r = calculate_drought_metrics(empty, changepoint=CP_CONFIG)
        for m in DROUGHT_METRICS:
            assert np.isnan(r[f"{m}_mean"])
            assert f"{m}_pettitt_cp_year" in r
        for s in DROUGHT_THRESHOLD_SCALARS:
            assert np.isnan(r[s])

    def test_missing_columns_emit_full_key_set(self):
        r = calculate_drought_metrics(pd.DataFrame({"Q": [1.0, 2.0]}),
                                      changepoint=CP_CONFIG)
        for m in DROUGHT_METRICS:
            assert np.isnan(r[f"{m}_mean"])

    def test_short_record_gets_nan_thresholds(self):
        """Fewer than min_years_for_threshold qualifying years -> NaN thresholds
        (and hence NaN metrics) rather than an unstable low tail."""
        gage = make_gage(lambda i, wy: 1.0 + 0.5*np.sin(i/30), n_years=5)
        r = calculate_drought_metrics(gage)
        for s in DROUGHT_THRESHOLD_SCALARS:
            assert np.isnan(r[s])
        for m in DROUGHT_METRICS:
            assert np.isnan(r[f"{m}_mean"])

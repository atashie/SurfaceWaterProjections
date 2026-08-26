"""
Tests for the fixed-b=1 recession alpha convention (Phase 2 of the port
campaign; port of julia/test/test_recession_alpha_b1.jl).

Convention (July 2026, canonical): ALL alpha outputs (log_a_pointcloud,
log_a_events, the log_a_seasonality_* scalars) assume a linear reservoir —
log(a) = median(log(-dQ/dt) - log(Q)), no regression — while b_pointcloud /
b_events / concavity keep their FREE power-law fits. The quadratic-reservoir
gage is the discriminating case: the free fit recovers b = 2 exactly while the
b=1 alpha differs strongly from the free-fit intercept log(k).

Run with: pytest tests/test_recession_alpha_b1.py -v
"""

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures import add_water_year_columns
from streamflow_signatures.recession import (
    analyze_recession_parameters,
    identify_recession_events,
    fit_recession_event,
    fit_sinusoidal_model,
    _event_log_a_b1,
)


def _daily_frame(n_years, q_update, peak=10.0):
    dates = pd.date_range("1990-10-01", periods=365 * n_years, freq="D")
    Q = np.empty(len(dates))
    q = peak
    for i, d in enumerate(dates):
        if d.day == 1:
            q = peak
        else:
            q = q_update(q, d)
        Q[i] = max(q, 1e-6)
    df = pd.DataFrame({"gage_id": "TESTB1", "date": dates, "Q": Q})
    return add_water_year_columns(df, date_col="date")


def b1_linear_gage(n_years=25, alpha=0.95, peak=10.0):
    return _daily_frame(n_years, lambda q, d: q * alpha, peak)


def b1_quadratic_gage(n_years=25, k=0.02, peak=10.0):
    return _daily_frame(n_years, lambda q, d: q - k * q * q, peak)


def b1_seasonal_gage(n_years=25, peak=10.0):
    # alpha = 0.90 Oct-Mar, 0.98 Apr-Sep -> per-event log_a alternates between
    # log(0.10) and log(0.02), keyed to the event's mid-event day-of-water-year
    return _daily_frame(
        n_years, lambda q, d: q * (0.90 if d.month in (10, 11, 12, 1, 2, 3) else 0.98), peak)


def _recompute_event_log_a_dowy(df):
    """Independent per-event (b=1 log_a, mid-event dowy) recompute mirroring
    the pipeline's event acceptance (free-fit success)."""
    la, dw = [], []
    for _, ydf in df.sort_values(["water_year", "dowy"]).groupby("water_year", sort=False):
        Q = ydf["Q"].values
        dowy = ydf["dowy"].values
        for ev in identify_recession_events(Q):
            Q_event = Q[ev["indices"]]
            log_a, b = fit_recession_event(Q_event, remove_first_day=True)
            if np.isnan(log_a) or np.isnan(b):
                continue
            Q_work = Q_event[1:]
            dQ = -np.diff(Q_work)
            Qm = Q_work[:-1]
            valid = (Qm > 0) & (dQ > 0)
            if not valid.any():
                continue
            la.append(np.median(np.log(dQ[valid]) - np.log(Qm[valid])))
            # canonical mid-event index: FLOOR midpoint of the index range
            # (julia `div(start+end, 2)`), not the upper-middle element
            idx = ev["indices"]
            dw.append(dowy[(idx[0] + idx[-1]) // 2])
    return np.array(la), np.array(dw)


def _recompute_log_a_pc_b1_means(df):
    """Independent per-year pointcloud b=1 recompute; returns the mean of the
    per-year medians (what log_a_pointcloud_mean summarizes)."""
    per_year = []
    for _, ydf in df.sort_values(["water_year", "dowy"]).groupby("water_year", sort=False):
        Q = ydf["Q"].values
        pool_q, pool_dq = [], []
        for ev in identify_recession_events(Q):
            Q_event = Q[ev["indices"]]
            log_a, b = fit_recession_event(Q_event, remove_first_day=True)
            if np.isnan(log_a) or np.isnan(b):
                continue
            Q_sub = Q_event[1:]
            dQ = -np.diff(Q_event[1:])
            Qm = Q_sub[:-1]
            valid = (Qm > 0) & (dQ > 0)
            pool_q.extend(Qm[valid])
            pool_dq.extend(dQ[valid])
        if len(pool_q) > 10:
            per_year.append(np.median(np.log(pool_dq) - np.log(pool_q)))
    return float(np.mean(per_year)), len(per_year)


class TestHelper:
    def test_pure_exponential_exact(self):
        Q_exp = 10.0 * 0.9 ** np.arange(20)
        assert _event_log_a_b1(Q_exp) == pytest.approx(np.log(0.1), abs=1e-12)

    def test_degenerate_inputs_nan(self):
        assert np.isnan(_event_log_a_b1(np.array([5.0, 4.0])))       # too short
        assert np.isnan(_event_log_a_b1(np.array([1.0, 2.0, 3.0, 4.0])))  # increasing


class TestLinearReservoir:
    """b=1 convention coincides with the free fit on a true linear reservoir."""

    @pytest.fixture(scope="class")
    def r(self):
        return analyze_recession_parameters(b1_linear_gage(alpha=0.95))

    def test_alpha_and_b(self, r):
        assert r["log_a_pointcloud_mean"] == pytest.approx(np.log(0.05), abs=1e-9)
        assert r["log_a_events_mean"] == pytest.approx(np.log(0.05), abs=1e-9)
        assert r["b_pointcloud_mean"] == pytest.approx(1.0, abs=1e-6)
        assert r["b_events_mean"] == pytest.approx(1.0, abs=1e-6)
        assert abs(r["concavity_mean"]) < 1e-9

    def test_discrete_alpha_consistency(self, r):
        assert r["recession_alpha_point_cloud_linear_reservoir"] == pytest.approx(0.95, abs=1e-9)
        assert r["log_a_pointcloud_mean"] == pytest.approx(
            np.log(1 - r["recession_alpha_point_cloud_linear_reservoir"]), abs=1e-9)

    def test_event_counts_unchanged(self, r):
        # Monthly spikes -> exactly 12 events/year
        assert r["n_recession_events_mean"] == pytest.approx(12.0, abs=1e-9)


class TestQuadraticReservoir:
    """The discriminating case: free b stays 2, alpha is decoupled from b."""

    @pytest.fixture(scope="class")
    def setup(self):
        k = 0.02
        df = b1_quadratic_gage(k=k)
        return k, df, analyze_recession_parameters(df)

    def test_free_b_recovers_two(self, setup):
        _, _, r = setup
        assert r["b_pointcloud_mean"] == pytest.approx(2.0, abs=1e-6)
        assert r["b_events_mean"] == pytest.approx(2.0, abs=1e-6)

    def test_b1_alpha_matches_independent_recompute(self, setup):
        _, df, r = setup
        expected_mean, n_years = _recompute_log_a_pc_b1_means(df)
        assert n_years >= 20
        assert r["log_a_pointcloud_mean"] == pytest.approx(expected_mean, abs=1e-9)

    def test_differs_from_free_fit_intercept(self, setup):
        # The OLD behavior (fitted-b intercept) gives ~log(k); the b=1 value is
        # median(log(k*Q_j)) with Q_j in (Q_min, peak) — strictly different.
        k, df, r = setup
        q_min = df["Q"].min()
        for key in ("log_a_pointcloud_mean", "log_a_events_mean"):
            assert np.log(k * q_min) < r[key] < np.log(k * 10.0)
            assert abs(r[key] - np.log(k)) > 0.3

    def test_seasonality_still_produced(self, setup):
        _, _, r = setup
        assert not np.isnan(r["log_a_seasonality_amplitude_all"])
        assert not np.isnan(r["log_a_seasonality_minimum_all"])


class TestSeasonalityAlignment:
    """The seasonality sinusoid must consume the per-event b=1 values, exactly
    aligned with mid-event dowy — reproduced by an independent recompute."""

    def test_alignment(self):
        df = b1_seasonal_gage()
        r = analyze_recession_parameters(df)

        la, dw = _recompute_event_log_a_dowy(df)
        assert len(la) >= 250  # ~12 events/yr x 25 yr
        amp, min_doy = fit_sinusoidal_model(dw, la)

        assert r["log_a_seasonality_amplitude_all"] == pytest.approx(amp, rel=1e-9)
        assert r["log_a_seasonality_minimum_all"] == pytest.approx(min_doy, rel=1e-9)

        # Non-trivial signal: log_a alternates between log(0.10) and log(0.02)
        assert r["log_a_seasonality_amplitude_all"] > 0.4
        # Minimum (most negative log_a = slow-recession season) in Apr-Sep
        assert 200 < r["log_a_seasonality_minimum_all"] < 350

    def test_mid_event_index_is_floor_midpoint(self):
        """Regression for the 2026-08-25 off-by-one: the seasonality timing must
        use the FLOOR midpoint of each event's index range (canonical Julia
        `div(start+end, 2)`; rpkg's 1-based `ceiling(n/2)` is the same). The old
        `indices[len//2]` put EVEN-length events one day late."""
        df = b1_seasonal_gage()
        ydf = df[df.water_year == 1992].sort_values("dowy")
        evs = identify_recession_events(ydf["Q"].values)
        even = [ev for ev in evs if len(ev["indices"]) % 2 == 0]
        assert even, "fixture must contain even-length events to exercise this"
        for ev in even:
            idx = ev["indices"]
            canonical = (idx[0] + idx[-1]) // 2
            old_buggy = idx[len(idx) // 2]
            assert canonical == old_buggy - 1  # the defect this pins

        # End-to-end: the pipeline's seasonality must equal a recompute that
        # uses the canonical midpoint (already asserted in test_alignment, but
        # pinned here against the specific formula)
        r = analyze_recession_parameters(df)
        la, dw = _recompute_event_log_a_dowy(df)
        amp, min_doy = fit_sinusoidal_model(dw, la)
        assert r["log_a_seasonality_minimum_all"] == pytest.approx(min_doy, rel=1e-12)

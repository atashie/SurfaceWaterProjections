"""
Regression tests mandated by the Codex Phase-1 exit review (2026-08-25):

- config resolution: STREAMFLOW_CONFIG override honored; packaged fallback copy
  stays byte-identical to the canonical repo config (drift guard).
- generate_stats frame-level gate: short frames emit keys for requested-but-
  absent metrics (Julia's earliest nrow < min_rows return).
- orchestrator topology: changepoint reaches the exempt families (recession,
  elasticity) while the stats floor does not.
- qp_seasonality / storage retain their changepoint keys end-to-end (the
  explicit-copy defect fixed in Phase 1).
- flashiness negative-total guard (canonical <= 0) at the annual level.
- pettitt: NaN-in-years filtering; pct_change 1e-10 boundary; frozen segment
  MK p-values (Python-frozen; cross-language p-values are the documented
  library-difference class, NOT expected bit-equal vs Julia).

Run with: pytest tests/test_phase1_review_fixes.py -v
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures import generate_stats, CP_SUFFIXES, STAT_SUFFIXES, calculate_all_signatures
from streamflow_signatures.changepoint import pettitt_test, segment_differential_metrics
from streamflow_signatures.flashiness import analyze_flashiness_trends

REPO = Path(__file__).parent.parent.parent
CP_CONFIG = {"start_year": 1900, "end_year": 2100,
             "min_total_obs": 20, "min_segment_obs": 10}


class TestConfigResolution:
    def test_packaged_copy_matches_canonical(self):
        """The wheel-bundled config must never drift from the canonical file."""
        canonical = json.loads((REPO / "config" / "signatures_config.json").read_text())
        packaged = json.loads(
            (REPO / "python" / "streamflow_signatures" / "data" / "signatures_config.json").read_text())
        assert packaged == canonical

    def test_streamflow_config_env_override(self, tmp_path):
        """STREAMFLOW_CONFIG must actually be consumed (fresh interpreter)."""
        variant = json.loads((REPO / "config" / "signatures_config.json").read_text())
        variant["stats_floor"]["min_values_for_stats"] = 37
        vpath = tmp_path / "variant.json"
        vpath.write_text(json.dumps(variant))
        out = subprocess.run(
            [sys.executable, "-c",
             "from streamflow_signatures.config import MIN_VALUES_FOR_STATS; print(MIN_VALUES_FOR_STATS)"],
            capture_output=True, text=True,
            env={"STREAMFLOW_CONFIG": str(vpath), "PATH": "/usr/bin:/bin"},
            cwd=str(REPO / "python"))
        assert out.returncode == 0, out.stderr
        assert out.stdout.strip() == "37"

    def test_missing_override_fails_loudly(self, tmp_path):
        out = subprocess.run(
            [sys.executable, "-c", "import streamflow_signatures.config"],
            capture_output=True, text=True,
            env={"STREAMFLOW_CONFIG": str(tmp_path / "nope.json"), "PATH": "/usr/bin:/bin"},
            cwd=str(REPO / "python"))
        assert out.returncode != 0
        assert "STREAMFLOW_CONFIG" in out.stderr


class TestFrameLevelGate:
    def test_short_frame_emits_keys_for_absent_metrics(self):
        """Julia's earliest gate emits keys for every REQUESTED metric, even ones
        absent from the frame — the per-column path would silently skip them."""
        df = pd.DataFrame({"water_year": [2000, 2001], "present": [1.0, 2.0]})
        result = generate_stats(df, value_cols=["present", "absent"], min_rows=3,
                                changepoint=CP_CONFIG)
        for col in ("present", "absent"):
            for suffix in STAT_SUFFIXES + CP_SUFFIXES:
                assert np.isnan(result[f"{col}{suffix}"]), f"{col}{suffix}"

    def test_short_frame_no_changepoint_no_cp_keys(self):
        df = pd.DataFrame({"water_year": [2000], "m": [1.0]})
        result = generate_stats(df, value_cols=["m"], min_rows=3)
        assert set(result) == {f"m{s}" for s in STAT_SUFFIXES}


def _topology_gage(n_years=21):
    """Deterministic daily gage: monthly exponential recessions (>=25 events),
    positive PPT, 21 water years."""
    rows = []
    for wy in range(2000, 2000 + n_years):
        dates = pd.date_range(f"{wy-1}-10-01", f"{wy}-09-30", freq="D")
        for i, d in enumerate(dates):
            day_in_month = d.day
            q = 10.0 * np.exp(-0.15 * day_in_month) + 0.05
            rows.append((d, q, wy, d.month, i + 1))
    df = pd.DataFrame(rows, columns=["date", "Q", "water_year", "month", "dowy"])
    df["gage_id"] = "SYN"
    df["PPT"] = 2.0 + (df["water_year"] % 5) * 0.3
    return df


class TestOrchestratorTopology:
    @pytest.fixture(scope="class")
    def out(self):
        df = _topology_gage()
        return calculate_all_signatures(
            df, has_climate=True, changepoint=CP_CONFIG,
            min_values_for_stats=25,  # > 21 years: floors every NON-exempt family
            gage_id="SYN")

    def test_floor_hits_non_exempt_families(self, out):
        assert np.isnan(out["Qann_mean"])          # flow volumes floored
        assert np.isnan(out["flashinessRB_mean"])  # flashiness floored

    def test_floor_exempts_recession_and_elasticity(self, out):
        assert not np.isnan(out["log_a_pointcloud_mean"])     # recession exempt
        assert not np.isnan(out["elasticity_rolling_mean"])   # elasticity exempt

    def test_changepoint_reaches_exempt_families(self, out):
        # Key EXISTENCE proves the kwarg threading (values may be NaN when the
        # series is shorter than min_total_obs)
        assert "log_a_pointcloud_pettitt_cp_year" in out
        assert "n_recession_events_pettitt_cp_year" in out
        assert "elasticity_rolling_pettitt_cp_year" in out

    def test_qp_and_storage_retain_cp_keys(self, out):
        for base in ("qp_slope_sd", "qp_bimodality", "avg_storage"):
            for suffix in CP_SUFFIXES:
                assert f"{base}{suffix}" in out, f"{base}{suffix} missing"


class TestFlashinessNegativeTotalGuard:
    def test_negative_total_year_non_computable(self):
        rows = []
        for wy, base in ((2001, 1.0), (2002, -0.2), (2003, 1.5)):
            for i in range(60):
                rows.append((wy, base + 0.01 * (i % 5), i + 1))
        df = pd.DataFrame(rows, columns=["water_year", "Q", "dowy"])
        assert df[df.water_year == 2002].Q.sum() < 0
        out = analyze_flashiness_trends(df)
        # 2 computable years < min_rows=3 -> all stats NaN; the pre-fix guard
        # (== 0) computed a negative-denominator value for 2002 and returned
        # finite statistics from 3 "computable" years
        assert np.isnan(out["flashinessRB_mean"])


class TestPettittEdges:
    def test_nan_in_years_filtered(self):
        yrs = np.arange(1990, 2015).astype(float)
        vals = np.concatenate([10 + 0.01 * np.arange(10), 20 + 0.01 * np.arange(15)])
        yrs_nan = yrs.copy()
        # NaN two YEAR entries (not values); pairs must be dropped identically
        yrs_nan[[2, 20]] = np.nan
        r = pettitt_test(yrs_nan, vals)
        mask = ~np.isnan(yrs_nan)
        r_ref = pettitt_test(yrs[mask], vals[mask])
        assert r == r_ref

    def test_pct_change_guard_boundary(self):
        """|pre_mean| must exceed 1e-10 strictly (Julia: `> 1e-10`). The exact
        ==1e-10 razor edge is not robustly testable in FP (a mean of equal
        values can land a ulp either side, identically fragile in Julia), so
        pin the guard a decade either side."""
        years = np.arange(2000, 2030).astype(float)
        below = np.concatenate([np.full(15, 1e-11), np.full(15, 5.0)])
        d = segment_differential_metrics(years, below, 2014.0)
        assert np.isnan(d["pct_change"])
        assert d["delta_mean"] == pytest.approx(5.0 - 1e-11)
        above = np.concatenate([np.full(15, 1e-9), np.full(15, 5.0)])
        d2 = segment_differential_metrics(years, above, 2014.0)
        assert np.isfinite(d2["pct_change"])

    def test_frozen_segment_mk_pvals(self):
        """Frozen segment MK p-values. Since 2026-08-26 these are the CANONICAL
        values (Python now uses Julia/R's continuity-corrected formula), so the
        freeze is a cross-language pin, not merely a drift guard."""
        yrs = np.arange(1990, 2015).astype(float)
        vals = np.concatenate([10 + 0.5 * np.arange(10), 30 + 0.25 * np.arange(15)])
        r = pettitt_test(yrs, vals)
        d = segment_differential_metrics(yrs, vals, r["cp_year"])
        assert r["cp_year"] == 2001.0
        assert d["pre_mk_pval"] == pytest.approx(8.303107353668793e-06, rel=1e-12)
        assert d["post_mk_pval"] == pytest.approx(2.631276625475465e-06, rel=1e-12)
        assert d == segment_differential_metrics(yrs, vals, r["cp_year"])  # deterministic


class TestSparseFamilyRowEmission:
    """Phase-3 cross-language comparator finding (2026-08-25): flashiness and
    both BFI functions must OMIT non-computable years, not emit NaN rows —
    matching Julia's `total_Q <= 0 -> continue`. Also guards the coverage gap
    that let a NameError into `analyze_baseflow_indices_with_parameters`:
    the orchestrator's try/except LOGS family failures, so a broken family
    shows up as missing columns rather than a test failure — these call the
    functions DIRECTLY."""

    def _frame(self, bad_year=None):
        import numpy as np
        import pandas as pd
        df = pd.DataFrame({
            "water_year": np.repeat(np.arange(2000, 2026), 40),
            "dowy": np.tile(np.arange(1, 41), 26),
            "Q": np.tile(np.linspace(5.0, 1.0, 40), 26),
        })
        if bad_year is not None:
            df.loc[df.water_year == bad_year, "Q"] = -1.0   # negative total
        return df

    def test_both_bfi_functions_return_finite(self):
        from streamflow_signatures.baseflow import (
            analyze_baseflow_indices, analyze_baseflow_indices_with_parameters)
        r1 = analyze_baseflow_indices(self._frame())
        r2 = analyze_baseflow_indices_with_parameters(self._frame(), 0.95)
        assert np.isfinite(r1["BFI_Eckhardt_mean"])
        assert np.isfinite(r2["BFI_Eckhardt_param_mean"])

    def test_negative_total_year_omitted_not_nan_rowed(self):
        from streamflow_signatures.baseflow import analyze_baseflow_indices
        from streamflow_signatures import AnnualCollector
        clean = analyze_baseflow_indices(self._frame())
        withbad = analyze_baseflow_indices(self._frame(bad_year=2010))
        # An omitted year leaves the mean untouched; a NaN row would too, so
        # assert on the COLLECTED ROW SET, which is what actually differs.
        c = AnnualCollector()
        analyze_baseflow_indices(self._frame(bad_year=2010), collector=c)
        years = [y for s, y in zip(c.signature, c.water_year) if s == "BFI_Eckhardt"]
        assert 2010 not in years, "non-computable year must be omitted, not NaN-rowed"
        assert len(years) == 25
        assert clean["BFI_Eckhardt_mean"] == pytest.approx(withbad["BFI_Eckhardt_mean"])

    def test_flashiness_omits_non_computable_year(self):
        from streamflow_signatures.flashiness import analyze_flashiness_trends
        from streamflow_signatures import AnnualCollector
        c = AnnualCollector()
        analyze_flashiness_trends(self._frame(bad_year=2010), collector=c)
        years = [y for s, y in zip(c.signature, c.water_year) if s == "flashinessRB"]
        assert 2010 not in years
        assert len(years) == 25


class TestMannKendallCanonicalFormula:
    """MK p-value method aligned to canonical Julia/R (2026-08-26, option c).

    Julia is a faithful port of R's `Kendall::MannKendall` (verified identical
    to 1.000 at n = 10…33); scipy's `kendalltau` was the outlier — it selects on
    TIES (not sample size), returning the exact p-value when untied and an
    asymptotic branch with NO continuity correction when tied. Ties dominate
    here (day counts, pulse counts, drought durations), which shifted ~0.24 %
    of main-path and ~0.45 % of Pettitt-segment significance calls.

    Every expectation below is the value canonical Julia produces, verified
    BIT-EXACT (tau and p) on 2026-08-26.
    """

    JULIA = {
        "untied_n10_up":    ([1, 3, 2, 5, 4, 7, 6, 9, 8, 10],
                             0.8222222222222222, 0.0012821837418053317),
        "tied_ints_n10":    ([1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
                             0.9428090415820634, 0.00037057721471467353),
        "tied_ints_n14":    ([3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7],
                             0.48892055554204805, 0.020200365358799655),
        "tied_ints_n33":    ([(i % 6) for i in range(1, 34)],
                             0.02658136567731214, 0.8500761168897997),
        "tied_heavy_n20":   ([1.0 if i <= 10 else 2.0 for i in range(1, 21)],
                             0.7254762501100116, 0.00018267179110953435),
        "decreasing_n15":   ([15 - i for i in range(1, 16)],
                             -1.0, 2.6515787521219636e-07),
        "days_like_n25":    ([180, 181, 179, 182, 180, 178, 183, 181, 180, 179, 184,
                              182, 180, 177, 181, 183, 180, 182, 179, 181, 180, 178,
                              182, 181, 180],
                             -0.007188851546895894, 0.9809721615938882),
        "n4":               ([1, 2, 2, 3], 0.9128709291752769, 0.14856177489186861),
    }

    @pytest.mark.parametrize("name", sorted(JULIA))
    def test_matches_julia_bit_exact(self, name):
        from streamflow_signatures.stats import mann_kendall_test
        vals, tau_j, p_j = self.JULIA[name]
        tau, p = mann_kendall_test(np.asarray(vals, dtype=float))
        assert tau == tau_j, f"tau {tau!r} != julia {tau_j!r}"
        assert p == p_j, f"pval {p!r} != julia {p_j!r}"

    def test_continuity_correction_is_applied(self):
        """The p-value must use z = (S-1)/sqrt(var_S), not S/sqrt(var_S) —
        i.e. it must NOT equal scipy's uncorrected asymptotic value."""
        from scipy import stats as sp
        from streamflow_signatures.stats import mann_kendall_test
        y = np.array([1., 1., 2., 2., 3., 3., 4., 4., 5., 5., 6., 6.])
        _, p = mann_kendall_test(y)
        p_scipy_asym = sp.kendalltau(np.arange(len(y)), y, method="asymptotic").pvalue
        assert not np.isclose(p, p_scipy_asym, rtol=1e-6)
        assert p > p_scipy_asym          # the correction shrinks |z| -> larger p

    def test_degenerate_contracts(self):
        from streamflow_signatures.stats import mann_kendall_test
        tau, p = mann_kendall_test(np.full(12, 4.0))      # constant -> var_S <= 0
        assert np.isnan(tau) and np.isnan(p)
        tau, p = mann_kendall_test(np.array([1.0, 2.0]))  # n < 3
        assert np.isnan(tau) and np.isnan(p)
        tau, p = mann_kendall_test(np.array([1.0, 2.0, 3.0]))  # n == 3 -> p = 1.0
        assert p == 1.0
        # NaNs dropped before anything else
        t1, p1 = mann_kendall_test(np.array([1., np.nan, 3., 2., np.nan, 5., 4.]))
        t2, p2 = mann_kendall_test(np.array([1., 3., 2., 5., 4.]))
        assert (t1, p1) == (t2, p2)

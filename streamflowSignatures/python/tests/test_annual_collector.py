"""
Tests for the annual-values collector (Phase 3 of the port campaign; port of
julia/test/test_annual_collector.jl).

Contract (julia/src/stats.jl): the collector receives the annual series EXACTLY
as passed to generate_stats — before min_rows, the stats floor, and the trend
gates — and is read-only with respect to the statistics. NaN rows are KEPT
("year present, metric not computable"); rows whose YEAR is missing/non-finite
are skipped (they cannot be keyed).

Run with: pytest tests/test_annual_collector.py -v
"""

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures import (
    AnnualCollector, generate_stats, calculate_all_signatures, add_water_year_columns,
)

CP_CONFIG = {"start_year": 1900, "end_year": 2100,
             "min_total_obs": 20, "min_segment_obs": 10}


def _annual_df(n=25, start=1990):
    return pd.DataFrame({
        "water_year": range(start, start + n),
        "m1": [10.0 + 0.5 * i for i in range(n)],
        "m2": [100.0 - 0.25 * i for i in range(n)],
    })


def _synthetic_gage(n_years=25):
    dates = pd.date_range("1990-10-01", periods=365 * n_years, freq="D")
    doy = dates.dayofyear.values
    rng = np.random.default_rng(11)
    Q = np.maximum(0.05, 2.0 + 1.5 * np.sin(2 * np.pi * (doy - 100) / 365)
                   + rng.random(len(dates)) * 0.4)
    df = pd.DataFrame({"gage_id": "SYN", "date": dates, "Q": Q})
    df = add_water_year_columns(df, date_col="date")
    df["PPT"] = np.maximum(0.0, df["Q"] * 2.5 + rng.random(len(df)))
    return df


class TestReadOnly:
    def test_collector_has_zero_effect_on_statistics(self):
        df = _annual_df()
        without = generate_stats(df, value_cols=["m1", "m2"], changepoint=CP_CONFIG)
        c = AnnualCollector()
        with_coll = generate_stats(df, value_cols=["m1", "m2"], changepoint=CP_CONFIG,
                                   collector=c)
        assert without == with_coll
        assert len(c) == 50   # 2 metrics x 25 years

    def test_default_none_leaves_behavior_unchanged(self):
        df = _annual_df()
        assert generate_stats(df, value_cols=["m1"]) == \
               generate_stats(df, value_cols=["m1"], collector=None)

    def test_full_gage_collector_does_not_change_signatures(self):
        gage = _synthetic_gage()
        a = calculate_all_signatures(gage, has_climate=True, changepoint=CP_CONFIG,
                                     min_values_for_stats=20, gage_id="SYN")
        c = AnnualCollector()
        b = calculate_all_signatures(gage, has_climate=True, changepoint=CP_CONFIG,
                                     min_values_for_stats=20, collector=c, gage_id="SYN")
        assert set(a) == set(b)
        for k in a:
            va, vb = a[k], b[k]
            if isinstance(va, float) and np.isnan(va):
                assert np.isnan(vb), k
            else:
                assert va == vb, k
        assert len(c) > 0


class TestCollectionSemantics:
    def test_collected_before_min_rows_gate(self):
        """A frame below min_rows returns all-NaN stats but is still collected."""
        df = pd.DataFrame({"water_year": [2000, 2001], "m": [1.0, 2.0]})
        c = AnnualCollector()
        res = generate_stats(df, value_cols=["m"], min_rows=3, collector=c)
        assert np.isnan(res["m_mean"])
        assert len(c) == 2                      # collection precedes the gate
        assert c.value == [1.0, 2.0]

    def test_collected_before_stats_floor(self):
        df = _annual_df(n=19)
        c = AnnualCollector()
        res = generate_stats(df, value_cols=["m1"], min_values_for_stats=20, collector=c)
        assert np.isnan(res["m1_mean"])         # floored
        assert len([s for s in c.signature if s == "m1"]) == 19   # still collected

    def test_nan_values_kept(self):
        df = _annual_df(n=10)
        df.loc[df.index[2:5], "m1"] = np.nan
        c = AnnualCollector()
        generate_stats(df, value_cols=["m1"], collector=c)
        vals = [v for s, v in zip(c.signature, c.value) if s == "m1"]
        assert len(vals) == 10                   # rows kept, not dropped
        assert sum(np.isnan(v) for v in vals) == 3

    def test_unkeyable_years_skipped(self):
        df = pd.DataFrame({"water_year": [2000.0, np.nan, 2002.0, np.inf],
                           "m": [1.0, 2.0, 3.0, 4.0]})
        c = AnnualCollector()
        generate_stats(df, value_cols=["m"], collector=c)
        assert c.water_year == [2000, 2002]
        assert c.value == [1.0, 3.0]

    def test_explicit_year_col(self):
        df = pd.DataFrame({"yr": [2000, 2001, 2002], "m": [1.0, 2.0, 3.0]})
        c = AnnualCollector()
        generate_stats(df, value_cols=["m"], year_col="yr", collector=c)
        assert c.water_year == [2000, 2001, 2002]

    def test_water_year_dtype_is_int(self):
        c = AnnualCollector()
        generate_stats(_annual_df(n=5), value_cols=["m1"], collector=c)
        assert all(isinstance(y, int) for y in c.water_year)


class TestFullGageIntegrity:
    @pytest.fixture(scope="class")
    def collected(self):
        gage = _synthetic_gage()
        c = AnnualCollector()
        res = calculate_all_signatures(gage, has_climate=True, changepoint=CP_CONFIG,
                                       min_values_for_stats=20, collector=c, gage_id="SYN")
        df = pd.DataFrame({"signature": c.signature, "water_year": c.water_year,
                           "value": c.value})
        return res, df

    def test_no_duplicate_keys(self, collected):
        _, df = collected
        dups = df.duplicated(subset=["signature", "water_year"]).sum()
        assert dups == 0, f"{dups} duplicate (signature, water_year) rows"

    def test_signature_coverage(self, collected):
        """Dense families must all appear in the collected output."""
        _, df = collected
        sigs = set(df["signature"])
        for expected in ["Qann", "Q50", "flashinessRB", "D50_day", "FDCall",
                         "BFI_Eckhardt", "TQmean", "annual_runoff_ratio",
                         "n_recession_events"]:
            assert expected in sigs, expected

    def test_recomputed_mean_median_reproduce_summary(self, collected):
        """The collector's rows must reproduce the summary _mean/_median exactly
        for ungated dense signatures — the check validate_annual_values.py runs
        at benchmark scale."""
        res, df = collected
        checked = 0
        for sig, grp in df.groupby("signature"):
            vals = grp["value"].to_numpy()
            vals = vals[~np.isnan(vals)]
            key_mean, key_med = f"{sig}_mean", f"{sig}_median"
            if key_mean not in res or len(vals) < 20:
                continue
            if np.isnan(res[key_mean]):
                continue
            assert res[key_mean] == pytest.approx(float(np.mean(vals)), rel=1e-12), key_mean
            assert res[key_med] == pytest.approx(float(np.median(vals)), rel=1e-12), key_med
            checked += 1
        assert checked >= 20, f"only {checked} signatures cross-checked"

    def test_qann_rows_match_independent_sums(self, collected):
        """Collected Qann values must equal independently computed annual sums."""
        _, df = collected
        gage = _synthetic_gage()
        expected = gage.groupby("water_year")["Q"].sum()
        got = df[df["signature"] == "Qann"].set_index("water_year")["value"]
        assert len(got) == len(expected)
        for wy, v in expected.items():
            assert got.loc[wy] == pytest.approx(v, rel=1e-12)

    def test_elasticity_key_years(self, collected):
        """elasticity_rolling is keyed to the END year of each 11-yr window;
        elasticity_annual to the later year of each consecutive pair."""
        _, df = collected
        n_years = 25
        roll = df[df["signature"] == "elasticity_rolling"]
        ann = df[df["signature"] == "elasticity_annual"]
        if len(roll):
            assert len(roll) == n_years - 11 + 1
        if len(ann):
            assert len(ann) == n_years - 1

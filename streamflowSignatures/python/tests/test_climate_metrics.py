"""
Climate metrics (guidelines 3.9; decisions 2026-10-12) — closed-form unit tests.

Port of julia/test/test_climate_metrics.jl. Every expected value is derived by hand
from the definitions in docs/plans/2026-10-12-climate-metrics-plan.md, never from the
implementation. Annual values are read back through the AnnualCollector (collected
BEFORE any gating), so a single hand-built water year pins a definition.

Run with: pytest tests/test_climate_metrics.py -v
"""

from datetime import date, timedelta

import numpy as np
import pandas as pd
import pytest

from streamflow_signatures.climate import (
    calculate_precipitation_metrics, calculate_temperature_metrics,
    PRECIP_METRICS, TEMP_METRICS,
)
from streamflow_signatures.io import preprocess_daily_data, add_water_year_columns
from streamflow_signatures import calculate_all_signatures, AnnualCollector

STAT_KEYS = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
             "_mk_rho", "_mk_pval", "_mean", "_median"]
CP_CONFIG = {"start_year": 1980, "end_year": 2024, "min_total_obs": 20, "min_segment_obs": 10}


def climate_year_frame(wy, P=None, TN=None, TX=None, Q=None):
    """One complete water year (Oct 1 of wy-1 .. Sep 30 of wy) with optional climate."""
    d0, d1 = date(wy - 1, 10, 1), date(wy, 9, 30)
    n = (d1 - d0).days + 1
    dates = pd.to_datetime([d0 + timedelta(days=i) for i in range(n)])
    df = pd.DataFrame({"gage_id": ["TCLIM"] * n, "date": dates,
                       "Q": np.ones(n) if Q is None else np.asarray(Q, dtype=float)})
    df = add_water_year_columns(df)
    if P is not None:
        df["PPT"] = np.asarray(P, dtype=float)
    if TN is not None:
        df["TMIN"] = np.asarray(TN, dtype=float)
    if TX is not None:
        df["TMAX"] = np.asarray(TX, dtype=float)
    return df


def climate_multi_year(wys, dry_every=7):
    frames = []
    for wy in wys:
        base = climate_year_frame(wy)
        n = len(base)
        d = np.arange(1, n + 1)
        P = np.where(d % dry_every == 0, 0.0, base["month"].to_numpy(dtype=float) + 0.01 * d)
        TX = 20.0 + 15.0 * np.sin(2 * np.pi * (d - 100) / 365)
        TN = TX - 12.0
        frames.append(climate_year_frame(wy, P=P, TN=TN, TX=TX))
    return pd.concat(frames, ignore_index=True)


def annual_lookup(c, sig, wy):
    for s, y, v in zip(c.signature, c.water_year, c.value):
        if s == sig and int(y) == wy:
            return v
    return None


def collect_precip(df):
    c = AnnualCollector()
    r = calculate_precipitation_metrics(df, collector=c)
    return r, c


def collect_temp(df, **kwargs):
    c = AnnualCollector()
    r = calculate_temperature_metrics(df, collector=c, **kwargs)
    return r, c


def test_registry_names():
    assert len(PRECIP_METRICS) == 20
    assert len(TEMP_METRICS) == 19
    assert len(set(PRECIP_METRICS + TEMP_METRICS)) == 39
    for name in ["prcp_5d_max", "prcp_50_dowy", "prcp_ndays_abv_10mm",
                 "prcp_ndays_abv_25mm", "prcp_ndays_abv_50mm"]:
        assert name in PRECIP_METRICS
    for s in ["win", "spr", "sum", "fal"]:
        assert f"prcp_{s}" in PRECIP_METRICS
        assert f"temp_max_{s}" in TEMP_METRICS and f"temp_min_{s}" in TEMP_METRICS
        assert f"gdd_cum_{s}" in TEMP_METRICS
    assert "prcp_sumr" not in PRECIP_METRICS and "prcp_fall" not in PRECIP_METRICS


def test_single_storm_year():
    wy = 2001
    P = np.zeros(365); P[99] = 100.0            # dowy 100 = Jan 8 2001 (month 1)
    r, c = collect_precip(climate_year_frame(wy, P=P))
    v = lambda s: annual_lookup(c, s, wy)
    assert v("prcp_annual") == 100.0
    assert v("prcp_1d_max") == 100.0 and v("prcp_5d_max") == 100.0
    assert v("prcp_max_dowy") == 100.0 and v("prcp_50_dowy") == 100.0
    assert v("prcp_cm_dowy") == 100.0
    assert v("prcp_ndays_abv_10mm") == 1.0 and v("prcp_ndays_abv_25mm") == 1.0
    assert v("prcp_ndays_abv_50mm") == 1.0
    assert v("prcp_annual_sd") == pytest.approx(np.std(P, ddof=1), rel=1e-12)
    assert v("prcp_annual_cv") == pytest.approx(np.std(P, ddof=1) / P.mean(), rel=1e-12)
    assert v("prcp_gini") == pytest.approx(364 / 365, rel=1e-12)
    assert v("prcp_seasonality_index") == pytest.approx(11 / 6, rel=1e-12)
    assert v("prcp_intermit_annual") == 265.0
    assert v("prcp_intermit_annual_dowy_start") == 101.0
    assert v("prcp_intermit_annual_dowy_end") == 365.0
    assert v("prcp_win") == 100.0
    assert v("prcp_spr") == 0.0 and v("prcp_sum") == 0.0 and v("prcp_fal") == 0.0
    for m in PRECIP_METRICS:
        for s in STAT_KEYS:
            assert f"{m}{s}" in r
    assert len(r) == 20 * 8


def test_uniform_year():
    wy = 2001
    _, c = collect_precip(climate_year_frame(wy, P=np.full(365, 2.0)))
    v = lambda s: annual_lookup(c, s, wy)
    assert v("prcp_annual") == 730.0
    assert abs(v("prcp_gini")) < 1e-12
    mlen = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    X = 2.0 * np.array(mlen); R = X.sum()
    assert v("prcp_seasonality_index") == pytest.approx(np.abs(X - R / 12).sum() / R, rel=1e-12)
    assert v("prcp_50_dowy") == 183.0
    assert v("prcp_cm_dowy") == pytest.approx(183.0, rel=1e-12)
    assert v("prcp_5d_max") == 10.0
    assert v("prcp_max_dowy") == 1.0
    assert v("prcp_intermit_annual") == 0.0
    assert np.isnan(v("prcp_intermit_annual_dowy_start"))
    assert np.isnan(v("prcp_intermit_annual_dowy_end"))
    assert v("prcp_annual_sd") == 0.0 and v("prcp_annual_cv") == 0.0
    assert v("prcp_ndays_abv_10mm") == 0.0


def test_threshold_boundaries():
    wy = 2001
    P = np.full(365, 5.0)
    P[9] = 0.1; P[10] = 0.1000001; P[11] = 0.0; P[12] = 0.1
    P[49] = 10.0; P[50] = 10.0001; P[51] = 25.0; P[52] = 50.0; P[53] = 50.5
    _, c = collect_precip(climate_year_frame(wy, P=P))
    v = lambda s: annual_lookup(c, s, wy)
    assert v("prcp_ndays_abv_10mm") == 4.0
    assert v("prcp_ndays_abv_25mm") == 2.0
    assert v("prcp_ndays_abv_50mm") == 1.0
    assert v("prcp_intermit_annual") == 2.0
    assert v("prcp_intermit_annual_dowy_start") == 12.0
    assert v("prcp_intermit_annual_dowy_end") == 13.0


def test_ties_and_window():
    wy = 2001
    P = np.zeros(365)
    P[0:5] = [10.0, 0.0, 0.0, 0.0, 0.0]
    P[360:365] = [1.0, 2.0, 3.0, 4.0, 5.0]
    P[199] = 10.0
    _, c = collect_precip(climate_year_frame(wy, P=P))
    assert annual_lookup(c, "prcp_max_dowy", wy) == 1.0
    assert annual_lookup(c, "prcp_5d_max", wy) == 15.0
    P2 = np.zeros(365); P2[[0, 100, 200, 300]] = 1.0
    _, c2 = collect_precip(climate_year_frame(wy, P=P2))
    assert annual_lookup(c2, "prcp_intermit_annual", wy) == 99.0
    assert annual_lookup(c2, "prcp_intermit_annual_dowy_start", wy) == 2.0
    assert annual_lookup(c2, "prcp_intermit_annual_dowy_end", wy) == 100.0


def test_seasonal_totals_fall_convention():
    wy = 2001
    base = climate_year_frame(wy)
    P = base["month"].to_numpy(dtype=float)
    _, c = collect_precip(climate_year_frame(wy, P=P))
    v = lambda s: annual_lookup(c, s, wy)
    assert v("prcp_win") == 12 * 31 + 1 * 31 + 2 * 28
    assert v("prcp_spr") == 3 * 31 + 4 * 30 + 5 * 31
    assert v("prcp_sum") == 6 * 30 + 7 * 31 + 8 * 31
    assert v("prcp_fal") == 9 * 30 + 10 * 31 + 11 * 30
    assert v("prcp_annual") == v("prcp_win") + v("prcp_spr") + v("prcp_sum") + v("prcp_fal")


def test_grid_guard():
    _, c = collect_precip(climate_year_frame(2000, P=np.ones(366)))
    assert annual_lookup(c, "prcp_annual", 2000) == 366.0
    assert annual_lookup(c, "prcp_50_dowy", 2000) == 183.0
    short = climate_year_frame(2001, P=np.ones(365)).iloc[:364]
    _, c2 = collect_precip(short)
    assert np.isnan(annual_lookup(c2, "prcp_annual", 2001))
    P_nan = np.ones(365); P_nan[49] = np.nan
    _, c3 = collect_precip(climate_year_frame(2001, P=P_nan))
    assert np.isnan(annual_lookup(c3, "prcp_annual", 2001))
    assert np.isnan(annual_lookup(c3, "prcp_gini", 2001))


def test_temperature_metrics():
    wy = 2001
    base = climate_year_frame(wy)
    TX = base["month"].to_numpy(dtype=float)
    TN = TX - 30.0
    _, c = collect_temp(climate_year_frame(wy, TN=TN, TX=TX))
    v = lambda s: annual_lookup(c, s, wy)
    assert v("temp_max_annual") == 12.0 and v("temp_max_annual_dowy") == 62.0
    assert v("temp_min_annual") == 1.0 - 30.0 and v("temp_min_annual_dowy") == 93.0
    assert v("temp_max_win") == 12.0 and v("temp_max_spr") == 5.0
    assert v("temp_max_sum") == 8.0 and v("temp_max_fal") == 11.0
    assert v("temp_min_win") == -29.0 and v("temp_min_spr") == -27.0
    assert v("temp_min_sum") == -24.0 and v("temp_min_fal") == -21.0
    assert v("frost_first_dowy") == 1.0 and v("frost_last_dowy") == 273.0
    assert v("gdd_cum_annual") == 0.0

    TN2 = np.full(365, 5.0); TN2[19] = -1.0; TN2[249] = -0.5; TN2[299] = -2.0
    TX2 = np.full(365, 15.0)
    _, c2 = collect_temp(climate_year_frame(wy, TN=TN2, TX=TX2))
    assert annual_lookup(c2, "frost_first_dowy", wy) == 20.0
    assert annual_lookup(c2, "frost_last_dowy", wy) == 250.0
    TN3 = np.full(365, 5.0); TN3[272] = -1.0; TN3[273] = -1.0
    _, c3 = collect_temp(climate_year_frame(wy, TN=TN3, TX=TX2))
    assert annual_lookup(c3, "frost_first_dowy", wy) == 273.0
    assert annual_lookup(c3, "frost_last_dowy", wy) == 273.0
    _, c4 = collect_temp(climate_year_frame(wy, TN=np.zeros(365), TX=TX2))
    assert np.isnan(annual_lookup(c4, "frost_first_dowy", wy))
    assert np.isnan(annual_lookup(c4, "frost_last_dowy", wy))

    _, c5 = collect_temp(climate_year_frame(wy, TN=np.full(365, 10.0), TX=np.full(365, 30.0)))
    assert annual_lookup(c5, "gdd_cum_annual", wy) == 20.0 * 365
    assert annual_lookup(c5, "gdd_cum_win", wy) == 20.0 * (31 + 31 + 28)
    assert annual_lookup(c5, "gdd_cum_fal", wy) == 20.0 * (31 + 30 + 30)
    _, c6 = collect_temp(climate_year_frame(wy, TN=np.full(365, -20.0), TX=np.full(365, 10.0)))
    assert annual_lookup(c6, "gdd_cum_annual", wy) == 0.0
    d = np.arange(1, 366)
    _, c7 = collect_temp(climate_year_frame(wy, TN=np.where(d <= 180, -10.0, 0.0),
                                            TX=np.where(d <= 180, 0.0, 10.0)))
    assert annual_lookup(c7, "gdd_cum_annual", wy) == 5.0 * 185
    TXn = TX2.copy(); TXn[6] = np.nan
    _, c8 = collect_temp(climate_year_frame(wy, TN=TN2, TX=TXn))
    assert np.isnan(annual_lookup(c8, "temp_max_annual", wy))
    assert np.isnan(annual_lookup(c8, "frost_first_dowy", wy))


def test_record_anchored_gate():
    frames = []
    for k, wy in enumerate(range(2001, 2021), start=1):
        n = (date(wy, 9, 30) - date(wy - 1, 10, 1)).days + 1
        TN = np.full(n, 5.0)
        if k <= 10:
            TN[30 + k - 1] = -1.0
        frames.append(climate_year_frame(wy, TN=TN, TX=np.full(n, 15.0)))
    df = pd.concat(frames, ignore_index=True)
    r_gated = calculate_temperature_metrics(df, trend_completeness=0.6, decade_completeness=0.8)
    assert np.isnan(r_gated["frost_first_dowy_senn_slp"])
    assert np.isnan(r_gated["frost_first_dowy_mk_pval"])
    assert not np.isnan(r_gated["frost_first_dowy_mean"])
    assert r_gated["frost_first_dowy_mean"] == pytest.approx(np.mean(range(31, 41)), abs=1e-12)
    assert not np.isnan(r_gated["temp_max_annual_senn_slp"])
    r_open = calculate_temperature_metrics(df)
    assert r_open["frost_first_dowy_senn_slp"] == pytest.approx(1.0, abs=1e-12)


def test_multi_year_statistics_and_identities():
    df = climate_multi_year(range(1991, 2021))
    rp = calculate_precipitation_metrics(df, trend_completeness=0.6, decade_completeness=0.8,
                                         min_values_for_stats=20, changepoint=CP_CONFIG)
    rt = calculate_temperature_metrics(df, trend_completeness=0.6, decade_completeness=0.8,
                                       min_values_for_stats=20, changepoint=CP_CONFIG)
    assert len(rp) == 20 * 16 and len(rt) == 19 * 16
    assert rp["prcp_annual_mean"] > 0
    assert 0.0 <= rp["prcp_gini_mean"] <= 1.0
    assert 0.0 <= rp["prcp_seasonality_index_mean"] <= 11 / 6
    assert rt["gdd_cum_annual_mean"] > 0
    c = AnnualCollector()
    calculate_temperature_metrics(df, collector=c)
    for wy in range(1991, 2021):
        s = sum(annual_lookup(c, f"gdd_cum_{k}", wy) for k in ["win", "spr", "sum", "fal"])
        assert s == pytest.approx(annual_lookup(c, "gdd_cum_annual", wy), rel=1e-12)
    r_short = calculate_precipitation_metrics(climate_multi_year(range(2011, 2021)), min_values_for_stats=20)
    assert all(np.isnan(v) for v in r_short.values())


def test_preprocessor_temperature_validity():
    frames = []
    for wy in range(2001, 2005):
        n = (date(wy, 9, 30) - date(wy - 1, 10, 1)).days + 1
        P = np.full(n, 2.0); TN = np.full(n, -15.0); TX = np.full(n, 5.0)
        if wy == 2002:
            TN[99:101] = np.nan
        elif wy == 2003:
            TX[49:85] = np.nan
        elif wy == 2004:
            TN[199:203] = np.nan
        frames.append(climate_year_frame(wy, P=P, TN=TN, TX=TX))
    df = pd.concat(frames, ignore_index=True)
    pp = preprocess_daily_data(df)
    assert pp["valid_years"] == [2001, 2002, 2003, 2004]
    assert pp["valid_climate_years"] == [2001, 2002, 2003, 2004]
    assert pp["valid_temp_years"] == [2001, 2002]
    y2002 = pp["data"][pp["data"]["water_year"] == 2002]
    assert not y2002["TMIN"].isna().any()
    assert y2002["TMIN"].iloc[99] == -15.0
    pp0 = preprocess_daily_data(df.drop(columns=["TMIN", "TMAX"]))
    assert pp0["valid_temp_years"] == []
    assert pp0["valid_years"] == pp["valid_years"]
    assert pp0["valid_climate_years"] == pp["valid_climate_years"]


def test_orchestrator_contract():
    df = climate_multi_year(range(1991, 2016))
    pp = preprocess_daily_data(df)
    data = pp["data"]
    cdata = data[data["water_year"].isin(set(pp["valid_climate_years"]))]
    temp_data = data[data["water_year"].isin(set(pp["valid_temp_years"]))]

    # no explicit climate_data -> no precipitation keys (no legacy fallback)
    r_bare = calculate_all_signatures(data, True)
    assert "prcp_annual_mean" not in r_bare and "temp_max_annual_mean" not in r_bare

    r_no_temp = calculate_all_signatures(data, True, climate_data=cdata)
    assert "prcp_annual_mean" in r_no_temp
    assert "temp_max_annual_mean" not in r_no_temp

    r_temp = calculate_all_signatures(data, True, climate_data=cdata, temp_data=temp_data)
    assert "temp_max_annual_mean" in r_temp and not np.isnan(r_temp["temp_max_annual_mean"])
    for k, v in r_no_temp.items():
        v2 = r_temp[k]
        assert (isinstance(v, float) and np.isnan(v) and np.isnan(v2)) or v == v2

    r_empty = calculate_all_signatures(data, True, climate_data=cdata, temp_data=temp_data.iloc[0:0])
    assert np.isnan(r_empty["temp_max_annual_mean"])
    n_temp_keys = sum(1 for k in r_empty if k.startswith(("temp_", "frost_", "gdd_")))
    assert n_temp_keys == 19 * 8

    r_q = calculate_all_signatures(data.drop(columns=["PPT", "TMIN", "TMAX"]), False)
    assert "prcp_annual_mean" not in r_q
    r_unnorm = calculate_all_signatures(data, True, climate_data=cdata, area_normalized=False)
    assert "prcp_annual_mean" in r_unnorm
    assert "annual_runoff_ratio_mean" not in r_unnorm

    bad = df.copy(); bad["PPT"] = -1.0
    ppb = preprocess_daily_data(bad)
    assert ppb["valid_climate_years"] == []
    assert "prcp_annual_mean" not in calculate_all_signatures(ppb["data"], True)


def test_round2_parity_regressions():
    """Codex round 2: interpolation arithmetic, signed zeros, zero gap, label format."""
    from streamflow_signatures.io import _interpolate_internal_gaps
    from streamflow_signatures.climate import _wet_threshold_label
    # Julia's expression order: left + (right - left) * offset / span
    x = np.array([-5.0, np.nan, np.nan, np.nan, 1.7])
    f = _interpolate_internal_gaps(x, 3)
    assert f[3] == -5.0 + (1.7 - -5.0) * 3 / 4           # 0.025000000000000355 in Julia
    x2 = np.array([-0.30000000000000004, np.nan, np.nan, np.nan, 0.1])
    f2 = _interpolate_internal_gaps(x2, 3)
    assert f2[3] == -0.30000000000000004 + (0.1 - -0.30000000000000004) * 3 / 4   # exactly 0.0 in Julia
    assert f2[3] == 0.0
    # a run longer than the gap stays entirely NaN; max_gap 0 fills nothing, no crash
    assert np.isnan(_interpolate_internal_gaps(np.array([1.0, np.nan, np.nan, np.nan, np.nan, 2.0]), 3)[1:5]).all()
    assert np.isnan(_interpolate_internal_gaps(np.array([1.0, np.nan, 2.0]), 0)[1])
    # signed zero: first numerical tie wins, value is +0.0
    wy = 2001
    P = np.zeros(365); P[0] = -0.0; P[1] = 0.0
    _, c = collect_precip(climate_year_frame(wy, P=P))
    assert annual_lookup(c, "prcp_max_dowy", wy) == 1.0
    assert str(annual_lookup(c, "prcp_1d_max", wy)) == "0.0"
    TN = np.full(365, 5.0); TN[0] = 0.0; TN[1] = -0.0
    _, c2 = collect_temp(climate_year_frame(wy, TN=TN, TX=np.full(365, 15.0)))
    assert annual_lookup(c2, "temp_min_annual_dowy", wy) == 1.0
    # infinite inputs void the year in the family kernels
    Pi = np.ones(365); Pi[20] = np.inf
    _, c3 = collect_precip(climate_year_frame(wy, P=Pi))
    assert np.isnan(annual_lookup(c3, "prcp_annual", wy)) and np.isnan(annual_lookup(c3, "prcp_50_dowy", wy))
    # the shared label format
    assert _wet_threshold_label(10) == "10" and _wet_threshold_label(0.5) == "0p5"
    assert _wet_threshold_label(1e-5) == "1em05" and _wet_threshold_label(2147483648) == "2147483648"


def test_round3_regressions():
    """Round 2: crossing tolerance, non-finite inputs are missing, canonical fraction."""
    from streamflow_signatures.io import _interpolate_internal_gaps
    from streamflow_signatures.config import CLIMATE_CUM_FRACTION
    wy = 2001
    P = np.zeros(365); P[0:7] = 0.1; P[300:307] = 0.1
    _, c = collect_precip(climate_year_frame(wy, P=P))
    assert annual_lookup(c, "prcp_50_dowy", wy) == 7.0
    _, c2 = collect_precip(climate_year_frame(2000, P=np.full(366, 0.1)))
    assert annual_lookup(c2, "prcp_50_dowy", 2000) == 183.0
    P3 = np.zeros(365); P3[0] = 49.0; P3[1] = 51.0
    _, c3 = collect_precip(climate_year_frame(wy, P=P3))
    assert annual_lookup(c3, "prcp_50_dowy", wy) == 2.0
    Pi = np.full(365, 2.0); Pi[99] = np.inf
    ppi = preprocess_daily_data(climate_year_frame(wy, P=Pi))
    assert ppi["valid_climate_years"] == [wy]
    assert ppi["data"]["PPT"].iloc[99] == 2.0
    Pn = np.full(365, 2.0); Pn[99] = -np.inf
    assert preprocess_daily_data(climate_year_frame(wy, P=Pn))["valid_climate_years"] == [wy]
    Pm = np.full(365, 2.0); Pm[99:103] = np.inf
    assert preprocess_daily_data(climate_year_frame(wy, P=Pm))["valid_climate_years"] == []
    Qi = np.ones(365); Qi[0:31] = np.inf
    assert preprocess_daily_data(climate_year_frame(wy, P=np.ones(365), Q=Qi))["valid_years"] == []
    f = _interpolate_internal_gaps(np.array([0.0, np.nan, np.inf]), 3)
    assert np.isnan(f[1]) and np.isinf(f[2])
    with pytest.raises(ValueError):
        preprocess_daily_data(climate_year_frame(wy, P=np.ones(365)), config={"max_gap_temp": -1})
    assert CLIMATE_CUM_FRACTION == 0.5
    # observation contract: strings are missing in the kernels (Julia's coalesce_q)
    dfs = climate_year_frame(wy, P=np.full(365, 2.0))
    dfs["PPT"] = dfs["PPT"].astype(object); dfs.loc[dfs.index[0], "PPT"] = "2.0"
    _, c4 = collect_precip(dfs)
    assert np.isnan(annual_lookup(c4, "prcp_annual", wy))


def test_crossing_tolerance_boundary():
    """Decision 2026-10-12: the 50 % crossing tolerates a relative deficit of 1e-12 of the
    total (CLIMATE_CROSSING_REL_TOL), so a mathematical tie is not decided by rounding.
    Two storms of 1 - d and 1 + d mm (total 2 mm): d = 1e-12 lies inside the allowance
    (first storm crosses, day 7); d = 4e-12 lies outside (second storm crosses, day 301)."""
    from streamflow_signatures.climate import CLIMATE_CROSSING_REL_TOL
    assert CLIMATE_CROSSING_REL_TOL == 1e-12
    wy = 2001
    for d, expected in ((1e-12, 7.0), (4e-12, 301.0)):
        P = np.zeros(365); P[6] = 1.0 - d; P[300] = 1.0 + d
        _, c = collect_precip(climate_year_frame(wy, P=P))
        assert annual_lookup(c, "prcp_50_dowy", wy) == expected, d


def test_preprocessing_fidelity_review():
    """Sonnet review 2026-10-12 (F1, F3, F4): PPT / SWE are written back only when the year
    is valid, a residual NaN always rejects (independent of reject_residual_na), NumPy
    integer limits are accepted, and a numeric string is a missing observation."""
    wy = 2001
    base = climate_year_frame(wy, P=np.full(365, 2.0))
    base["SWE"] = 1.0
    # leading gap: cannot be filled -> residual -> rejected even with reject_residual_na False
    df = base.copy(); df.loc[0, "PPT"] = np.nan; df.loc[0, "SWE"] = np.nan
    pp = preprocess_daily_data(df, config={"reject_residual_na": False})
    assert pp["valid_climate_years"] == [] and pp["valid_swe_years"] == []
    assert np.isnan(pp["data"]["PPT"].iloc[0]) and np.isnan(pp["data"]["SWE"].iloc[0])
    # an invalid year keeps its raw column: 4-day internal gap stays unfilled
    df = base.copy(); df.loc[100:103, "PPT"] = np.nan
    pp = preprocess_daily_data(df)
    assert pp["valid_climate_years"] == [] and pp["data"]["PPT"].iloc[100:104].isna().all()
    # a valid year is filled with Julia's expression
    df = base.copy(); df.loc[100, "PPT"] = np.nan
    pp = preprocess_daily_data(df)
    assert pp["valid_climate_years"] == [wy] and pp["data"]["PPT"].iloc[100] == 2.0
    # NumPy integer limits are whole numbers >= 0
    pp = preprocess_daily_data(base.copy(), config={"max_raw_na_ppt": np.int64(30), "max_gap_swe": np.int32(3)})
    assert pp["valid_climate_years"] == [wy]
    with pytest.raises(ValueError):
        preprocess_daily_data(base.copy(), config={"max_raw_na_ppt": np.bool_(True)})
    # numeric strings are not observations (Julia's _real_vector contract) for PPT, SWE and Q
    df = base.copy(); df["PPT"] = df["PPT"].astype(object); df.loc[100, "PPT"] = "3"
    assert preprocess_daily_data(df)["valid_climate_years"] == [wy]          # one missing day, filled
    df = base.copy(); df["PPT"] = df["PPT"].astype(object); df.loc[100:103, "PPT"] = "3"
    assert preprocess_daily_data(df)["valid_climate_years"] == []            # 4-day gap of strings
    df = base.copy(); df["Q"] = df["Q"].astype(object); df.loc[100:103, "Q"] = "1"
    assert preprocess_daily_data(df)["valid_years"] == []


def test_review_regressions():
    wy = 2001
    P = np.zeros(365); P[0:50] = 0.1; P[300:350] = 0.1
    _, c = collect_precip(climate_year_frame(wy, P=P))
    assert annual_lookup(c, "prcp_50_dowy", wy) == 50.0
    s = 0.0
    for v in P.tolist():   # plain sequential sum (the builtin sum is compensated since 3.12)
        s += v
    assert annual_lookup(c, "prcp_annual", wy) == s
    trunc = climate_year_frame(2000, P=np.full(366, 2.0), TN=np.full(366, -1.0), TX=np.full(366, 5.0)).iloc[:365]
    _, c2 = collect_precip(trunc)
    assert np.isnan(annual_lookup(c2, "prcp_annual", 2000))
    _, c3 = collect_temp(trunc)
    assert np.isnan(annual_lookup(c3, "gdd_cum_annual", 2000))
    df_m = climate_year_frame(wy, P=np.full(365, 2.0))
    df_m.loc[df_m.index[0], "month"] = np.nan
    _, c4 = collect_precip(df_m)
    parts = sum(annual_lookup(c4, f"prcp_{s}", wy) for s in ["win", "spr", "sum", "fal"])
    assert parts == 730.0

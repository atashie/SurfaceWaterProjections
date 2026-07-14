"""
Tests for the area_normalized gate on Q-to-PPT signatures (July 2026).

Gages with no drainage area (HYDAT gap -- canals, dam outflows, channel splits)
carry Q in raw m3/s instead of mm/day. Q-to-PPT signatures (runoff ratios,
elasticity, Q-P seasonality, storage) mix Q and PPT units and are therefore
skipped when area_normalized=False. Q-only signatures must be unaffected.

Run with: pytest tests/test_area_normalized_gate.py -v
"""

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from streamflow_signatures import calculate_all_signatures, add_water_year_columns

# Prefixes of every Q-to-PPT (climate-dependent) output column
CLIMATE_SIGNATURE_PREFIXES = (
    "annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
    "summer_runoff_ratio", "fall_runoff_ratio", "runoff_ratio_high_count",
    "elasticity_static", "elasticity_rolling", "elasticity_annual",
    "elasticity_years_total", "elasticity_years_low_ppt",
    "qp_slope_sd", "qp_bimodality",
    "avg_storage",
)

# Representative Q-only columns that must survive the gate
QONLY_SENTINELS = (
    "Qann_mean", "Q50_mean", "BFI_Eckhardt_mean", "flashinessRB_mean",
    "D50_day_mean", "TQmean_mean", "FDCall_mean",
)


@pytest.fixture
def gage_with_climate():
    rng = np.random.default_rng(7)
    n_days = 365 * 25
    dates = pd.date_range("1990-10-01", periods=n_days, freq="D")
    doy = dates.dayofyear.values
    q = np.maximum(0.01, 2.0 + 1.5 * np.sin(2 * np.pi * (doy - 100) / 365)
                   + rng.random(n_days) * 0.5)
    df = pd.DataFrame({"gage_id": "TESTGATE1", "date": dates, "Q": q})
    df = add_water_year_columns(df)
    df["PPT"] = np.maximum(0.0, df["Q"] * 2.5 + rng.random(len(df)))
    return df


def _is_climate_key(key):
    return key.startswith(CLIMATE_SIGNATURE_PREFIXES)


def _same(a, b):
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    try:
        return bool(a == b)
    except Exception:
        return False


def test_default_computes_climate_signatures(gage_with_climate):
    r = calculate_all_signatures(gage_with_climate, True)
    assert any(_is_climate_key(k) for k in r)
    assert "annual_runoff_ratio_mean" in r
    assert "avg_storage_mean" in r
    for k in QONLY_SENTINELS:
        assert k in r


def test_gate_skips_all_climate_signatures(gage_with_climate):
    r = calculate_all_signatures(gage_with_climate, True, area_normalized=False)
    assert not [k for k in r if _is_climate_key(k)]
    for k in QONLY_SENTINELS:
        assert k in r


def test_gated_identical_to_no_climate(gage_with_climate):
    r_gated = calculate_all_signatures(gage_with_climate, True, area_normalized=False)
    r_noclimate = calculate_all_signatures(
        gage_with_climate.drop(columns="PPT"), False)
    assert set(r_gated) == set(r_noclimate)
    mismatches = [k for k in r_gated if not _same(r_gated[k], r_noclimate[k])]
    assert not mismatches


def test_gate_inert_when_true(gage_with_climate):
    r_default = calculate_all_signatures(gage_with_climate, True)
    r_explicit = calculate_all_signatures(gage_with_climate, True, area_normalized=True)
    assert set(r_default) == set(r_explicit)
    mismatches = [k for k in r_default if not _same(r_default[k], r_explicit[k])]
    assert not mismatches

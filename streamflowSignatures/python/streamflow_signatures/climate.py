"""
Climate metrics (guidelines §3.9, synced 2026-10-12; decisions: user, 2026-10-12).

Port of julia/src/climate.jl (canonical). Thirty-nine per-water-year metrics from
basin-mean Daymet daily precipitation (``prcp`` → PPT, mm/day) and air temperature
(``tmin`` → TMIN, ``tmax`` → TMAX, °C), in two function families:

- ``calculate_precipitation_metrics`` — 20 bases on PPT-valid years (the frame the
  Q-to-PPT families receive; pure climate, so it runs regardless of ``area_normalized``).
- ``calculate_temperature_metrics`` — 19 bases on an explicit ``temp_data`` frame
  filtered to temperature-valid years (``preprocess_daily_data``'s ``valid_temp_years``;
  the snow family's opt-in precedent — TMIN/TMAX inside the main gage frame are never
  used implicitly).

Conventions (design record: docs/plans/2026-10-12-climate-metrics-plan.md):
- Seasons are the config's ``season_definitions`` (DJF, MAM, JJA, SON) applied WITHIN
  the water year, the Qwin/Qspr/Qsum/Qfal convention: ``fal`` = Oct + Nov of the
  water-year start plus Sep of its end. Suffixes ``win``, ``spr``, ``sum``, ``fal``.
- dowy 1 = Oct 1. Ties take the first occurrence (maximum day, longest dry spell).
- A year is used only when its daily grid is complete (365/366 rows, dowy == 1..n, no
  NaN); anything else leaves the whole row NaN (defensive — callers pass qualified years).
- No Q-derived seasonal mask: these are climate series, not streamflow.
- A dry day has P <= ``dry_day_threshold_mm``; wet-day counts are STRICTLY above each
  threshold; the multi-day maximum is a window inside the water year; growing degree
  days are Σ max(0, (TMAX + TMIN)/2 − ``gdd_base_c``); frost is TMIN < ``frost_temp_c``.
- Gating: standard trend completeness and stats floor via ``generate_stats``; the
  record-anchored decade gate (snow precedent) on ``frost_first_dowy`` /
  ``frost_last_dowy``, anchored to the temperature-valid years.
- Config: the ``climate`` section (CLIMATE_*); an absent section disables the family in
  the orchestrator (no climate columns — the additivity-proof switch).
"""

import numbers
from calendar import monthrange
from datetime import date, timedelta
from typing import Dict, List, Optional, Set

import numpy as np
import pandas as pd

from .stats import generate_stats
from .config import (
    CLIMATE_DRY_DAY_MM,
    CLIMATE_WET_THRESHOLDS_MM,
    CLIMATE_MULTIDAY_WINDOW,
    CLIMATE_CUM_FRACTION,
    CLIMATE_FROST_TEMP_C,
    CLIMATE_FROST_LAST_END_MONTH,
    CLIMATE_GDD_BASE_C,
    CLIMATE_RECORD_DECADE_GATE,
    NA_SEASONAL_DEFINITIONS,
)

_SEASON_KEYS = ["win", "spr", "sum", "fal"]
_SEASON_FULL = {"win": "winter", "spr": "spring", "sum": "summer", "fal": "fall"}
_SEASON_DEFAULTS = {"winter": [12, 1, 2], "spring": [3, 4, 5],
                    "summer": [6, 7, 8], "fall": [9, 10, 11]}


def _season_months(key: str) -> Set[int]:
    full = _SEASON_FULL[key]
    return set(NA_SEASONAL_DEFINITIONS.get(full, _SEASON_DEFAULTS[full]))   # validated ints (config.py)


def _real_values(series) -> np.ndarray:
    """Observation contract shared with Julia's coalesce_q: a value is an observation when
    it is a real number (bool included, as Julia's Bool <: Real); strings, complex numbers
    and None are missing. Numeric dtypes pass straight through."""
    # fast path only for plain NumPy float / int / unsigned / bool storage; extension dtypes
    # (e.g. an Arrow decimal column, whose scalars are Decimal, not numbers.Real) and
    # object columns go through the scalar contract (review round 3)
    if isinstance(series.dtype, np.dtype) and series.dtype.kind in "fiub":
        return series.to_numpy(dtype=float)
    out = np.full(len(series), np.nan)
    for i, v in enumerate(series.tolist()):
        if isinstance(v, numbers.Real) and not isinstance(v, complex):
            out[i] = float(v)
    return out


def _wet_threshold_label(thr: float) -> str:
    """Column-name fragment for a wet-day threshold, identical in the three languages:
    an integral value prints as an integer ("10"), anything else with C's %g ("0.5",
    "1e-05"), then "." -> "p" and "-" -> "m" ("0p5", "1em05")."""
    thr = float(thr)
    s = str(int(thr)) if (thr.is_integer() and abs(thr) < 1e15) else ("%g" % thr)
    return s.replace(".", "p").replace("-", "m")


# Config-derived column names (the D{p}_day pattern): a changed window, fraction or
# threshold list changes the names, never silently the meaning behind a fixed name.
CLIMATE_MULTIDAY_MAX_COL = f"prcp_{int(CLIMATE_MULTIDAY_WINDOW)}d_max"
CLIMATE_CUM_DOWY_COL = f"prcp_{int(round(100 * float(CLIMATE_CUM_FRACTION)))}_dowy"
CLIMATE_WET_DAY_COLS = [f"prcp_ndays_abv_{_wet_threshold_label(t)}mm" for t in CLIMATE_WET_THRESHOLDS_MM]

PRECIP_METRICS: List[str] = (
    ["prcp_annual"]
    + [f"prcp_{s}" for s in _SEASON_KEYS]
    + ["prcp_1d_max", CLIMATE_MULTIDAY_MAX_COL, "prcp_max_dowy", CLIMATE_CUM_DOWY_COL, "prcp_cm_dowy"]
    + CLIMATE_WET_DAY_COLS
    + ["prcp_annual_sd", "prcp_annual_cv"]
    + ["prcp_intermit_annual", "prcp_intermit_annual_dowy_start", "prcp_intermit_annual_dowy_end"]
    + ["prcp_gini", "prcp_seasonality_index"]
)

if len(set(PRECIP_METRICS)) != len(PRECIP_METRICS):
    raise ValueError(f"climate: the configured window / thresholds produce duplicate precipitation column names: {PRECIP_METRICS}")

TEMP_METRICS: List[str] = (
    ["temp_max_annual", "temp_max_annual_dowy", "temp_min_annual", "temp_min_annual_dowy"]
    + [f"temp_max_{s}" for s in _SEASON_KEYS]
    + [f"temp_min_{s}" for s in _SEASON_KEYS]
    + ["frost_first_dowy", "frost_last_dowy"]
    + ["gdd_cum_annual"]
    + [f"gdd_cum_{s}" for s in _SEASON_KEYS]
)

# Timing metrics that are NaN in frost-free years: gated like the snow timing metrics
TEMP_RECORD_GATED_METRICS = ["frost_first_dowy", "frost_last_dowy"]

# Relative tolerance of the cumulative-fraction crossing (mirrors julia/src/climate.jl
# CLIMATE_CROSSING_REL_TOL): a mathematical tie must not be decided by the last bit
CLIMATE_CROSSING_REL_TOL = 1e-12


def _seqsum(x) -> float:
    """Plain left-to-right double accumulation, as an explicit loop. NOT the builtin
    ``sum`` — since Python 3.12 it uses compensated (Neumaier) summation for floats,
    which is more accurate but differs from Julia's loop and from ``_seqcumsum``.
    Every sum, mean and cumulative sum in this family uses this helper, so a
    cumulative sum's last element is exactly the total it is compared against
    (review finding 1, 2026-10-12) and the annual values are bit-identical with
    Julia and R."""
    s = 0.0
    for v in np.asarray(x, dtype=float).tolist():
        s += v
    return s


def _seqcumsum(x: np.ndarray) -> np.ndarray:
    out = np.empty(len(x), dtype=float)
    s = 0.0
    for i, v in enumerate(np.asarray(x, dtype=float).tolist()):
        s += v
        out[i] = s
    return out


def _water_year_days(yr: int) -> int:
    return (date(yr, 9, 30) - date(yr - 1, 10, 1)).days + 1


def _complete_year_grid(dowy: np.ndarray, n: int, yr: int) -> bool:
    """True when the year's rows are EXACTLY dowy 1..n and n is the number of days
    of water year yr (review finding 3: a leap water year truncated to 365 rows
    must not pass)."""
    if n != _water_year_days(yr):
        return False
    if np.isnan(dowy).any():
        return False
    return bool(np.array_equal(dowy, np.arange(1, n + 1, dtype=float)))


def _grid_months(yr: int, n: int) -> np.ndarray:
    """Calendar month of each day of water year yr on the complete grid, derived from
    the date rather than read from a month column (review finding 4)."""
    d0 = date(yr - 1, 10, 1)
    return np.array([(d0 + timedelta(days=i)).month for i in range(n)], dtype=int)


def _empty_family_stats(metrics: List[str], **kwargs) -> Dict[str, float]:
    """The full NaN key set of a family via generate_stats on an empty annual frame."""
    empty_annual = pd.DataFrame({"water_year": pd.Series([], dtype="int64")})
    for m in metrics:
        empty_annual[m] = pd.Series([], dtype="float64")
    return generate_stats(empty_annual, value_cols=metrics, year_col="water_year", **kwargs)


def _record_anchored_force_skip(annual: pd.DataFrame, anchor_col: str,
                                gated: List[str], decade_completeness: float) -> Optional[Set[str]]:
    """Record-anchored decade gate (the snow family's rule): the gated metrics need
    >= decade_completeness of the anchor years (rows where anchor_col is non-NaN) in
    BOTH the first and the last decade of the record to be non-NaN, else their six
    trend statistics are forced to NaN. Record span < 10 years -> no gate."""
    anchor_mask = ~annual[anchor_col].isna().to_numpy()
    anchor_years = annual.loc[anchor_mask, "water_year"].astype(int).to_numpy()
    if anchor_years.size == 0:
        return None
    rec_min, rec_max = int(anchor_years.min()), int(anchor_years.max())
    if rec_max - rec_min + 1 < 10:
        return None
    windows = ((rec_min, rec_min + 9), (rec_max - 9, rec_max))
    all_years = annual["water_year"].astype(int).to_numpy()
    skip: Set[str] = set()
    for m_name in gated:
        vals = annual[m_name].to_numpy(dtype=float)
        for lo, hi in windows:
            in_win = anchor_mask & (all_years >= lo) & (all_years <= hi)
            den = int(in_win.sum())
            if den == 0:
                continue
            num = int((~np.isnan(vals[in_win])).sum())
            if (num / den) < decade_completeness:
                skip.add(m_name)
                break
    return skip if skip else None


def _stats_kwargs(trend_completeness, decade_completeness, min_values_for_stats,
                  changepoint, collector):
    return dict(trend_completeness=trend_completeness,
                decade_completeness=decade_completeness,
                min_values_for_stats=min_values_for_stats,
                changepoint=changepoint, collector=collector)


def calculate_precipitation_metrics(
    df: pd.DataFrame,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """Twenty precipitation metrics with trend statistics from daily basin-mean
    precipitation (PPT, mm/day), per water year. See the module docstring and
    julia/src/climate.jl for the definitions; callers pass the PPT-valid frame."""
    kw = _stats_kwargs(trend_completeness, decade_completeness, min_values_for_stats,
                       changepoint, collector)
    required = ["PPT", "water_year", "dowy"]
    if any(c not in df.columns for c in required):
        return _empty_family_stats(PRECIP_METRICS, **kw)

    dry_thr = float(CLIMATE_DRY_DAY_MM)
    wet_thrs = [float(t) for t in CLIMATE_WET_THRESHOLDS_MM]
    w = int(CLIMATE_MULTIDAY_WINDOW)
    frac = float(CLIMATE_CUM_FRACTION)
    season_sets = {s: _season_months(s) for s in _SEASON_KEYS}

    years = list(pd.unique(df["water_year"]))
    annual = pd.DataFrame({"water_year": years})
    for m in PRECIP_METRICS:
        annual[m] = np.nan
    idx_of = {yr: i for i, yr in enumerate(years)}

    for yr, year_df in df.groupby("water_year", sort=False):
        yr_idx = idx_of.get(yr)
        if yr_idx is None or len(year_df) == 0:
            continue
        year_df = year_df.sort_values("dowy")

        # `+ 0.0` turns a signed -0.0 into +0.0 so extrema and ties do not depend on the
        # sign of zero; the guard requires FINITE values (±Inf is not data)
        P = _real_values(year_df["PPT"]) + 0.0
        dowy = pd.to_numeric(year_df["dowy"], errors="coerce").to_numpy(dtype=float)
        n = len(P)

        # Complete-grid + finite guard (snow precedent, year-aware): anything else -> NaN row
        if not _complete_year_grid(dowy, n, int(yr)) or not np.isfinite(P).all():
            continue
        months = _grid_months(int(yr), n)

        total = _seqsum(P)
        annual.at[yr_idx, "prcp_annual"] = total

        # Seasonal totals over the season's calendar months within the water year
        for s in _SEASON_KEYS:
            smask = np.isin(months, list(season_sets[s]))
            annual.at[yr_idx, f"prcp_{s}"] = _seqsum(P[smask])

        # 1-day maximum and its first day
        imax = int(np.argmax(P))
        annual.at[yr_idx, "prcp_1d_max"] = float(P[imax])
        annual.at[yr_idx, "prcp_max_dowy"] = float(dowy[imax])

        # Multi-day maximum: every window of w consecutive days inside the year
        if n >= w:
            best = -np.inf
            for i in range(0, n - w + 1):
                s_ = _seqsum(P[i:i + w])
                if s_ > best:
                    best = s_
            annual.at[yr_idx, CLIMATE_MULTIDAY_MAX_COL] = best

        # Timing of the cumulative total, centre of mass, Gini, seasonality index
        # (the cumulative sum and `total` share one accumulation order, so the last
        # cumulative value equals the total exactly; `frac * total` is exact for 0.5)
        if total > 0.0:
            cum = _seqcumsum(P)
            thr = frac * total - CLIMATE_CROSSING_REL_TOL * total
            reached = np.flatnonzero(cum >= thr)
            if reached.size:
                annual.at[yr_idx, CLIMATE_CUM_DOWY_COL] = float(dowy[reached[0]])
            annual.at[yr_idx, "prcp_cm_dowy"] = _seqsum(dowy * P) / total

            y = np.sort(P)
            i_idx = np.arange(1, n + 1, dtype=float)
            annual.at[yr_idx, "prcp_gini"] = _seqsum((2.0 * i_idx - n - 1.0) * y) / (n * total)

            X = np.array([_seqsum(P[months == m]) for m in range(1, 13)])
            R = _seqsum(X)
            if R > 0.0:
                annual.at[yr_idx, "prcp_seasonality_index"] = _seqsum(np.abs(X - R / 12.0)) / R

        # Wet-day counts, strictly above each threshold
        for col, t in zip(CLIMATE_WET_DAY_COLS, wet_thrs):
            annual.at[yr_idx, col] = float(int((P > t).sum()))

        # Within-year spread of the daily values (sample sd, n - 1; sequential sums)
        mu = total / n
        sd = float(np.sqrt(_seqsum((P - mu) ** 2) / (n - 1)))
        annual.at[yr_idx, "prcp_annual_sd"] = sd
        annual.at[yr_idx, "prcp_annual_cv"] = sd / mu if mu > 0.0 else np.nan

        # Longest dry spell (first on ties) and its first / last day
        dry = P <= dry_thr
        best_len = 0
        best_start = -1
        i = 0
        while i < n:
            if dry[i]:
                j = i
                while j < n - 1 and dry[j + 1]:
                    j += 1
                length = j - i + 1
                if length > best_len:
                    best_len = length
                    best_start = i
                i = j + 1
            else:
                i += 1
        annual.at[yr_idx, "prcp_intermit_annual"] = float(best_len)
        if best_len > 0:
            annual.at[yr_idx, "prcp_intermit_annual_dowy_start"] = float(dowy[best_start])
            annual.at[yr_idx, "prcp_intermit_annual_dowy_end"] = float(dowy[best_start + best_len - 1])

    return generate_stats(annual, value_cols=PRECIP_METRICS, year_col="water_year", **kw)


def calculate_temperature_metrics(
    df: pd.DataFrame,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """Nineteen temperature metrics with trend statistics from daily basin-mean
    minimum and maximum air temperature (TMIN, TMAX, °C), per water year. Callers
    pass the temperature-valid frame (``temp_data``); see julia/src/climate.jl."""
    kw = _stats_kwargs(trend_completeness, decade_completeness, min_values_for_stats,
                       changepoint, collector)
    required = ["TMIN", "TMAX", "water_year", "dowy"]
    if any(c not in df.columns for c in required):
        return _empty_family_stats(TEMP_METRICS, **kw)

    frost_t = float(CLIMATE_FROST_TEMP_C)
    end_month = int(CLIMATE_FROST_LAST_END_MONTH)
    base = float(CLIMATE_GDD_BASE_C)
    season_sets = {s: _season_months(s) for s in _SEASON_KEYS}

    years = list(pd.unique(df["water_year"]))
    annual = pd.DataFrame({"water_year": years})
    for m in TEMP_METRICS:
        annual[m] = np.nan
    idx_of = {yr: i for i, yr in enumerate(years)}

    for yr, year_df in df.groupby("water_year", sort=False):
        yr_idx = idx_of.get(yr)
        if yr_idx is None or len(year_df) == 0:
            continue
        year_df = year_df.sort_values("dowy")

        TN = _real_values(year_df["TMIN"]) + 0.0
        TX = _real_values(year_df["TMAX"]) + 0.0
        dowy = pd.to_numeric(year_df["dowy"], errors="coerce").to_numpy(dtype=float)
        n = len(TN)

        if not _complete_year_grid(dowy, n, int(yr)) or not np.isfinite(TN).all() or not np.isfinite(TX).all():
            continue
        months = _grid_months(int(yr), n)

        # Annual extremes and their first day
        imax = int(np.argmax(TX))
        imin = int(np.argmin(TN))
        annual.at[yr_idx, "temp_max_annual"] = float(TX[imax])
        annual.at[yr_idx, "temp_max_annual_dowy"] = float(dowy[imax])
        annual.at[yr_idx, "temp_min_annual"] = float(TN[imin])
        annual.at[yr_idx, "temp_min_annual_dowy"] = float(dowy[imin])

        # Growing degree days, clamped at the base temperature
        gdd = np.maximum(0.0, (TX + TN) / 2.0 - base)
        annual.at[yr_idx, "gdd_cum_annual"] = _seqsum(gdd)

        # Seasonal extremes and GDD over the season's months within the water year
        for s in _SEASON_KEYS:
            smask = np.isin(months, list(season_sets[s]))
            if smask.any():
                annual.at[yr_idx, f"temp_max_{s}"] = float(TX[smask].max())
                annual.at[yr_idx, f"temp_min_{s}"] = float(TN[smask].min())
                annual.at[yr_idx, f"gdd_cum_{s}"] = _seqsum(gdd[smask])

        # Frost: first day below the threshold; last such day up to the end of
        # frost_last_end_month (a month >= 10 belongs to the water year's first
        # calendar year). dowy == 1..n after the grid guard.
        frost = TN < frost_t
        first = np.flatnonzero(frost)
        if first.size:
            annual.at[yr_idx, "frost_first_dowy"] = float(dowy[first[0]])
        yr_i = int(yr)
        cal_year = yr_i - 1 if end_month >= 10 else yr_i
        cutoff_dowy = (date(cal_year, end_month, monthrange(cal_year, end_month)[1])
                       - date(yr_i - 1, 10, 1)).days + 1
        last = np.flatnonzero(frost & (dowy <= cutoff_dowy))
        if last.size:
            annual.at[yr_idx, "frost_last_dowy"] = float(dowy[last[-1]])

    force_skip = None
    if CLIMATE_RECORD_DECADE_GATE and decade_completeness is not None:
        force_skip = _record_anchored_force_skip(annual, "temp_max_annual",
                                                 TEMP_RECORD_GATED_METRICS, decade_completeness)

    return generate_stats(annual, value_cols=TEMP_METRICS, year_col="water_year",
                          force_skip_trends=force_skip, **kw)

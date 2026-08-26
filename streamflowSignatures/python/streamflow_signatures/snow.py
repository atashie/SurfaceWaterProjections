"""
Snow metrics signatures.

Fourteen per-water-year metrics computed from daily snow water equivalent (SWE,
mm; Daymet `swe`). Port of julia/src/snow.jl (canonical). ALL metrics operate on
a thresholded series ``SWE* = SWE if SWE >= SNOW_SWE_THRESHOLD_MM else 0`` —
days below the threshold are treated as snow-free for durations AND magnitudes
(domain decision 2026-07-14).

Definitions (full specification: docs/plans/snow_signatures_plan.md):
- Snow-day mask: SWE* > 0; spells are maximal runs of consecutive snow days.
- Anchor spell: the spell containing the annual maximum of SWE* (first-max tie
  rule). snow_on_dowy / snow_off_dowy are censored (NaN) when the anchor spell
  touches the first / last day of the water year.
- Melt increments m_t = max(0, SWE*_{t-1} - SWE*_t), attributed to day t;
  within-year only (no attribution across the Sep 30 -> Oct 1 boundary).
- SSM (snow seasonality metric; Hatchett 2021, Hydrology 8(1):32): spells of
  >= SNOW_SEASONAL_MIN_DAYS days are seasonal, shorter spells ephemeral;
  SSM = (seasonal_days - ephemeral_days) / total snow days, in [-1, +1].

Year qualification: callers pass data already filtered to SWE-valid years
(``preprocess_daily_data``'s ``valid_swe_years``). Defensively, any year whose
SWE series has NaNs or is not a complete Oct 1 - Sep 30 daily grid is left
entirely NaN (all 14 metrics).

References: Hatchett, B.J. (2021). Seasonal and Ephemeral Snowpacks of the
Conterminous United States. Hydrology 8(1), 32. Petersky & Harpold (2018), HESS.
"""

from datetime import date
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .stats import generate_stats
from .config import (
    SNOW_SWE_THRESHOLD_MM,
    SNOW_SEASONAL_MIN_DAYS,
    SNOW_MELT_COM_FRACTION,
    SNOW_MIN_ANNUAL_PPT_MM,
    SNOW_RECORD_DECADE_GATE,
)

SNOW_METRICS = [
    "swe_max",
    "swe_max_dowy",
    "snow_cover_days",
    "snow_on_dowy",
    "snow_off_dowy",
    "melt_season_days",
    "melt_rate",
    "ssm",
    "swe_apr1",
    "melt_before_peak",
    "melt_before_peak_pct",
    "melt_before_peak_to_max_swe",
    "melt_com_dowy",
    "swe_max_to_ppt",
]

# Threshold-dependent metrics (NaN for operationally snow-free years) gated by the
# record-anchored decade gate (July 2026): a trend requires computable values in
# >= decade_completeness of the SWE-valid years in BOTH the first and last decade
# of the gage's SWE record. Magnitude metrics (valid zeros) are NOT gated — their
# dense series carry the snow-decline signal legitimately.
SNOW_RECORD_GATED_METRICS = [
    "swe_max_dowy",
    "snow_on_dowy",
    "snow_off_dowy",
    "melt_season_days",
    "melt_rate",
    "ssm",
    "melt_before_peak",
    "melt_before_peak_pct",
    "melt_before_peak_to_max_swe",
    "melt_com_dowy",
]


def snow_spells(mask: np.ndarray) -> List[range]:
    """Maximal runs of consecutive True values, as 0-based ranges."""
    spells = []
    n = len(mask)
    i = 0
    while i < n:
        if mask[i]:
            j = i
            while j < n - 1 and mask[j + 1]:
                j += 1
            spells.append(range(i, j + 1))
            i = j + 1
        else:
            i += 1
    return spells


def calculate_snow_metrics(
    streamflow_data: pd.DataFrame,
    valid_climate_years: Optional[List[int]] = None,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate 14 snow metrics with trend statistics from daily SWE.

    Magnitude metrics (swe_max, snow_cover_days, swe_apr1, swe_max_to_ppt) emit
    valid ZEROS for operationally snow-free years; timing/melt/regime metrics
    emit NaN. Each metric produces 8 statistics via generate_stats() (+ 8
    changepoint fields when enabled).

    Callers must pass data filtered to SWE-valid years. In
    ``calculate_all_signatures`` this only runs when an explicit ``snow_data``
    frame is provided — an SWE column inside the main gage frame is never used
    implicitly.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily data with columns: water_year, dowy, SWE (+ optional PPT)
    valid_climate_years : list of int, optional
        Years qualified for PPT use (metric swe_max_to_ppt). When None (legacy
        path), the annual-PPT convention of the runoff ratios applies directly.

    Returns
    -------
    dict
        14 metrics x 8 statistics (+ 8 changepoint fields each when enabled)
    """
    df = streamflow_data

    # Defensive: callers gate on SWE presence; direct-API misuse still emits the
    # IDENTICAL key set (incl. changepoint fields when enabled) via the same
    # explicit-value_cols machinery as the 0-row case — schema contract.
    if "SWE" not in df.columns:
        empty_annual = pd.DataFrame({"water_year": pd.Series([], dtype="int64")})
        for m in SNOW_METRICS:
            empty_annual[m] = pd.Series([], dtype="float64")
        return generate_stats(
            empty_annual, value_cols=SNOW_METRICS, year_col="water_year",
            trend_completeness=trend_completeness,
            decade_completeness=decade_completeness,
            min_values_for_stats=min_values_for_stats,
            changepoint=changepoint, collector=collector)

    thr = float(SNOW_SWE_THRESHOLD_MM)
    seasonal_min = int(SNOW_SEASONAL_MIN_DAYS)
    com_frac = float(SNOW_MELT_COM_FRACTION)
    min_ppt = float(SNOW_MIN_ANNUAL_PPT_MM)

    has_ppt = "PPT" in df.columns
    climate_set = None if valid_climate_years is None else set(valid_climate_years)

    years = list(pd.unique(df["water_year"]))
    annual = pd.DataFrame({"water_year": years})
    for m in SNOW_METRICS:
        annual[m] = np.nan
    idx_of = {yr: i for i, yr in enumerate(years)}

    for yr, year_df in df.groupby("water_year", sort=False):
        yr_idx = idx_of.get(yr)
        if yr_idx is None or len(year_df) == 0:
            continue
        year_df = year_df.sort_values("dowy")

        swe_raw = pd.to_numeric(year_df["SWE"], errors="coerce").to_numpy(dtype=float)
        dowy = pd.to_numeric(year_df["dowy"], errors="coerce").to_numpy(dtype=float)
        n = len(swe_raw)

        # Complete-grid + no-NaN guard: metrics are only computed on full water
        # years whose days are EXACTLY the sequence 1..n (n in {365, 366}) with no
        # missing SWE. Endpoint checks alone would accept a year with a duplicated
        # interior day masking a missing one.
        if n not in (365, 366) or np.isnan(dowy).any() \
                or not np.array_equal(dowy, np.arange(1, n + 1, dtype=float)) \
                or np.isnan(swe_raw).any():
            continue

        # Thresholded series: sub-threshold days are snow-free everywhere
        swe_star = np.where(swe_raw >= thr, swe_raw, 0.0)
        mask = swe_star > 0.0
        n_snow = int(mask.sum())

        # --- Magnitude metrics (valid zeros for operationally snow-free years) ---
        swe_max_val = float(swe_star.max()) if n_snow > 0 else 0.0
        annual.at[yr_idx, "swe_max"] = swe_max_val
        annual.at[yr_idx, "snow_cover_days"] = float(n_snow)

        # April 1 SWE: calendar April 1 of water year yr (leap-safe; dowy == 1..n
        # after the grid guard, so this indexes directly)
        yr_i = int(yr)
        apr1_dowy = (date(yr_i, 4, 1) - date(yr_i - 1, 10, 1)).days + 1
        annual.at[yr_idx, "swe_apr1"] = float(swe_star[apr1_dowy - 1])

        # --- Timing / melt / regime metrics (need a qualifying snowpack) ---
        if n_snow > 0:
            peak_idx = int(np.argmax(swe_star))   # first occurrence on ties
            annual.at[yr_idx, "swe_max_dowy"] = float(dowy[peak_idx])

            spells = snow_spells(mask)

            # SSM pools ALL spells in the year
            seasonal_days = sum(len(sp) for sp in spells if len(sp) >= seasonal_min)
            ephemeral_days = sum(len(sp) for sp in spells if len(sp) < seasonal_min)
            annual.at[yr_idx, "ssm"] = (seasonal_days - ephemeral_days) / float(n_snow)

            # Anchor spell = the spell containing the peak
            anchor = next(sp for sp in spells if peak_idx in sp)
            a_first, a_last = anchor[0], anchor[-1]

            # Censoring: an anchor spell touching a water-year boundary -> NaN
            snow_on = np.nan if a_first == 0 else float(dowy[a_first])
            snow_off = np.nan if a_last == n - 1 else float(dowy[a_last] + 1)
            annual.at[yr_idx, "snow_on_dowy"] = snow_on
            annual.at[yr_idx, "snow_off_dowy"] = snow_off

            if not np.isnan(snow_off):
                melt_days = snow_off - float(dowy[peak_idx])
                annual.at[yr_idx, "melt_season_days"] = melt_days
                annual.at[yr_idx, "melt_rate"] = swe_max_val / melt_days

            # Melt increments, attributed to the day the drop lands (day t).
            # With the first-max tie rule, m at the peak day itself is 0.
            m = np.zeros(n)
            drops = swe_star[:-1] - swe_star[1:]
            m[1:] = np.where(drops > 0.0, drops, 0.0)
            total_melt = float(m.sum())
            melt_before = float(m[:peak_idx + 1].sum())
            annual.at[yr_idx, "melt_before_peak"] = melt_before
            annual.at[yr_idx, "melt_before_peak_to_max_swe"] = melt_before / swe_max_val
            if total_melt > 0.0:
                annual.at[yr_idx, "melt_before_peak_pct"] = 100.0 * melt_before / total_melt
                cum = np.cumsum(m)
                reached = np.flatnonzero(cum >= com_frac * total_melt)
                if reached.size:
                    annual.at[yr_idx, "melt_com_dowy"] = float(dowy[reached[0]])

        # --- swe_max_to_ppt: needs a PPT-qualified year ---
        if has_ppt and (climate_set is None or yr in climate_set):
            # Same annual-PPT convention as the runoff ratios: sum the year's
            # valid (non-NaN) PPT days, require the total above the minimum floor
            P = pd.to_numeric(year_df["PPT"], errors="coerce").to_numpy(dtype=float)
            p_valid = ~np.isnan(P)
            total_P = float(P[p_valid].sum()) if p_valid.any() else 0.0
            if total_P > min_ppt:
                annual.at[yr_idx, "swe_max_to_ppt"] = swe_max_val / total_P

    # Record-anchored decade gate (July 2026, user decision 2026-07-22): anchored
    # to the SWE-valid, grid-complete years (swe_max non-NaN — dense incl. zeros),
    # NOT the metric's own span, so clustered snowy years cannot slip through. The
    # threshold is the SAME decade_completeness the streamflow gate uses (linked
    # knob). Denominators count anchor years in each window — a year with no SWE
    # data must not count against snow presence. Record span < 10 -> gate skipped.
    force_skip = None
    if SNOW_RECORD_DECADE_GATE and decade_completeness is not None:
        anchor_mask = ~annual["swe_max"].isna().to_numpy()
        anchor_years = annual.loc[anchor_mask, "water_year"].astype(int).to_numpy()
        if anchor_years.size:
            rec_min, rec_max = int(anchor_years.min()), int(anchor_years.max())
            if rec_max - rec_min + 1 >= 10:
                windows = ((rec_min, rec_min + 9), (rec_max - 9, rec_max))
                all_years = annual["water_year"].astype(int).to_numpy()
                skip = set()
                for m_name in SNOW_RECORD_GATED_METRICS:
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
                if skip:
                    force_skip = skip

    return generate_stats(
        annual, value_cols=SNOW_METRICS, year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness,
        min_values_for_stats=min_values_for_stats,
        changepoint=changepoint, collector=collector,
        force_skip_trends=force_skip)

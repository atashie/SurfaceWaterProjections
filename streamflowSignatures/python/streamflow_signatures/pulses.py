"""
Pulse metrics signatures.

Calculates high and low flow pulse characteristics:
- Pulse counts and durations
- TQmean (percentage of days above mean)
- Flow reversals (direction changes)
"""

from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import (
    PULSES_MIN_DAYS,
    HIGH_PULSE_PERCENTILE,
    LOW_PULSE_PERCENTILE,
    FLOW_REVERSAL_THRESHOLD,
)


def identify_pulses(
    flow: np.ndarray,
    threshold: float,
    above: bool = True
) -> Dict:
    """
    Identify flow pulses (consecutive days above/below threshold).

    Parameters
    ----------
    flow : np.ndarray
        Daily discharge values
    threshold : float
        Flow threshold for pulse identification
    above : bool, default True
        If True, identify pulses above threshold; if False, below

    Returns
    -------
    dict
        Contains 'n_pulses', 'durations', and 'mean_duration'
    """
    if above:
        exceeds = flow > threshold
    else:
        exceeds = flow < threshold

    # Handle NAs by treating as False
    exceeds = np.where(np.isnan(flow), False, exceeds)

    # Find runs of consecutive True values
    # Pad with False at ends to detect edge pulses
    padded = np.concatenate([[False], exceeds, [False]])
    diffs = np.diff(padded.astype(int))

    # Starts are where diff = 1, ends are where diff = -1
    starts = np.where(diffs == 1)[0]
    ends = np.where(diffs == -1)[0]

    if len(starts) == 0:
        return {
            'n_pulses': 0,
            'durations': [],
            'mean_duration': np.nan
        }

    durations = ends - starts

    return {
        'n_pulses': len(durations),
        'durations': durations.tolist(),
        'mean_duration': np.mean(durations)
    }


def count_flow_reversals(
    flow: np.ndarray,
    threshold_pct: float = FLOW_REVERSAL_THRESHOLD
) -> int:
    """
    Count flow reversals (direction changes in flow).

    A reversal occurs when flow changes from increasing to decreasing
    or vice versa, with the change exceeding a threshold.

    Parameters
    ----------
    flow : np.ndarray
        Daily discharge values
    threshold_pct : float, default 0.02
        Minimum change as fraction of current flow (2%)

    Returns
    -------
    int
        Number of flow reversals
    """
    # Interpolate NA gaps to avoid false adjacencies
    if np.any(np.isnan(flow)):
        non_na_idx = np.where(~np.isnan(flow))[0]
        if len(non_na_idx) < 3:
            return 0
        flow_clean = np.interp(
            np.arange(len(flow)),
            non_na_idx,
            flow[non_na_idx]
        )
    else:
        flow_clean = flow.copy()

    n = len(flow_clean)
    if n < 3:
        return 0

    reversal_count = 0

    for i in range(1, n - 1):
        # Calculate changes
        prev_change = flow_clean[i] - flow_clean[i - 1]
        next_change = flow_clean[i + 1] - flow_clean[i]

        # Threshold based on current flow
        threshold = abs(threshold_pct * flow_clean[i])

        # Check for reversal with significant change
        if abs(next_change) > threshold:
            # Was increasing (or flat), now decreasing
            if prev_change >= 0 and next_change < 0:
                reversal_count += 1
            # Was decreasing (or flat), now increasing
            elif prev_change <= 0 and next_change > 0:
                reversal_count += 1

    return reversal_count


def calculate_pulse_metrics(
    streamflow_data: pd.DataFrame,
    min_nona_days: int = PULSES_MIN_DAYS,
) -> Dict[str, float]:
    """
    Calculate pulse metrics and flow reversal trends.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - dowy: int, day of water year
            - month: int, calendar month (1-12)
    min_nona_days : int, default 250
        Minimum non-NA days required per water year

    Returns
    -------
    dict
        Dictionary of signature statistics with keys for:
            - n_high_pulses_year: Count of high pulse events (>90th percentile)
            - n_low_pulses_year: Count of low pulse events (<10th percentile)
            - dur_high_pulses_year: Mean duration of high pulses
            - dur_low_pulses_year: Mean duration of low pulses
            - TQmean: Percentage of days above annual mean
            - Flow_Reversals_*: Annual and seasonal reversal counts

    Notes
    -----
    - High pulses: days > 90th percentile of that year's flow
    - Low pulses: days < 10th percentile of that year's flow
    - TQmean: Typically ~30-35% for natural streams
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "dowy", "month"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Calculate overall percentiles for entire period (from config)
    q90_all = df["Q"].quantile(HIGH_PULSE_PERCENTILE)
    q10_all = df["Q"].quantile(LOW_PULSE_PERCENTILE)

    years = df["water_year"].unique()

    # Initialize pulse metrics DataFrame
    pulse_metrics = pd.DataFrame({
        "water_year": years,
        "n_high_pulses_year": np.nan,
        "n_low_pulses_year": np.nan,
        "dur_high_pulses_year": np.nan,
        "dur_low_pulses_year": np.nan,
        # Period-of-record threshold metrics (matching R's *_all metrics)
        "n_high_pulses_all": np.nan,
        "n_low_pulses_all": np.nan,
        "dur_high_pulses_all": np.nan,
        "dur_low_pulses_all": np.nan,
        "TQmean": np.nan,
        "Flow_Reversals_annual": np.nan,
        "Flow_Reversals_winter": np.nan,
        "Flow_Reversals_spring": np.nan,
        "Flow_Reversals_summer": np.nan,
        "Flow_Reversals_fall": np.nan,
    })

    # Process each year
    for yr, year_data in df.groupby("water_year", sort=False):
        year_data = year_data.copy()  # needed because we mutate below (sort_values)

        # Skip years with insufficient data
        non_na_count = year_data["Q"].notna().sum()
        if non_na_count < min_nona_days:
            continue

        # Sort by day of water year
        year_data = year_data.sort_values("dowy")

        Q = year_data["Q"].values
        months = year_data["month"].values

        # Calculate year-specific thresholds (from config)
        q90_year = np.nanquantile(Q, HIGH_PULSE_PERCENTILE)
        q10_year = np.nanquantile(Q, LOW_PULSE_PERCENTILE)

        if np.isnan(q90_year) or np.isnan(q10_year):
            continue

        # Analyze pulses for year-specific thresholds
        high_pulses = identify_pulses(Q, q90_year, above=True)
        low_pulses = identify_pulses(Q, q10_year, above=False)

        # Analyze pulses for period-of-record thresholds (*_all metrics)
        high_pulses_all = identify_pulses(Q, q90_all, above=True)
        low_pulses_all = identify_pulses(Q, q10_all, above=False)

        # Calculate TQmean
        annual_mean = np.nanmean(Q)
        days_above_mean = np.nansum(Q > annual_mean)
        total_days = np.sum(~np.isnan(Q))
        tqmean_pct = (days_above_mean / total_days) * 100 if total_days > 0 else np.nan

        # Calculate annual flow reversals
        annual_reversals = count_flow_reversals(Q)

        # Calculate seasonal flow reversals
        winter_mask = np.isin(months, [12, 1, 2])
        spring_mask = np.isin(months, [3, 4, 5])
        summer_mask = np.isin(months, [6, 7, 8])
        fall_mask = np.isin(months, [9, 10, 11])

        winter_reversals = count_flow_reversals(Q[winter_mask]) if winter_mask.sum() >= 30 else np.nan
        spring_reversals = count_flow_reversals(Q[spring_mask]) if spring_mask.sum() >= 30 else np.nan
        summer_reversals = count_flow_reversals(Q[summer_mask]) if summer_mask.sum() >= 30 else np.nan
        fall_reversals = count_flow_reversals(Q[fall_mask]) if fall_mask.sum() >= 30 else np.nan

        # Store results using year-based lookup
        yr_mask = pulse_metrics["water_year"] == yr
        # Store results - year-specific thresholds
        pulse_metrics.loc[yr_mask, "n_high_pulses_year"] = high_pulses['n_pulses']
        pulse_metrics.loc[yr_mask, "n_low_pulses_year"] = low_pulses['n_pulses']
        pulse_metrics.loc[yr_mask, "dur_high_pulses_year"] = high_pulses['mean_duration']
        pulse_metrics.loc[yr_mask, "dur_low_pulses_year"] = low_pulses['mean_duration']
        # Store results - period-of-record thresholds (*_all metrics)
        pulse_metrics.loc[yr_mask, "n_high_pulses_all"] = high_pulses_all['n_pulses']
        pulse_metrics.loc[yr_mask, "n_low_pulses_all"] = low_pulses_all['n_pulses']
        pulse_metrics.loc[yr_mask, "dur_high_pulses_all"] = high_pulses_all['mean_duration']
        pulse_metrics.loc[yr_mask, "dur_low_pulses_all"] = low_pulses_all['mean_duration']
        pulse_metrics.loc[yr_mask, "TQmean"] = tqmean_pct
        pulse_metrics.loc[yr_mask, "Flow_Reversals_annual"] = annual_reversals
        pulse_metrics.loc[yr_mask, "Flow_Reversals_winter"] = winter_reversals
        pulse_metrics.loc[yr_mask, "Flow_Reversals_spring"] = spring_reversals
        pulse_metrics.loc[yr_mask, "Flow_Reversals_summer"] = summer_reversals
        pulse_metrics.loc[yr_mask, "Flow_Reversals_fall"] = fall_reversals

    # Get metric columns
    metric_cols = [c for c in pulse_metrics.columns if c != "water_year"]

    # Ensure numeric types
    for col in metric_cols:
        pulse_metrics[col] = pd.to_numeric(pulse_metrics[col], errors="coerce")

    # Generate statistics
    return generate_stats(
        pulse_metrics,
        value_cols=metric_cols,
        year_col="water_year"
    )

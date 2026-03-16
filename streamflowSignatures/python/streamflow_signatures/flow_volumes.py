"""
Flow volume signatures.

Calculates annual and seasonal flow totals and percentiles.
"""

from typing import Dict
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import FLOW_VOLUMES_MIN_DAYS, FLOW_PERCENTILES


def calculate_flow_vols_by_year(
    streamflow_data: pd.DataFrame,
    min_nona_days: int = FLOW_VOLUMES_MIN_DAYS,
) -> Dict[str, float]:
    """
    Calculate annual and seasonal flow volumes and percentiles.

    This function calculates 22 flow volume metrics for each water year,
    then applies generate_stats() to produce 8 statistics per metric
    (176 total output columns).

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - month: int, calendar month (1-12)
            - dowy: int, day of water year (1-366)
    min_nona_days : int, default 250
        Minimum number of non-NA days required per water year

    Returns
    -------
    dict
        Dictionary of signature statistics. For each of 22 metrics,
        8 statistics are produced (see generate_stats).

        Metrics:
            - Qann: Annual total flow (mm)
            - Qwin: Winter total (Dec-Feb) (mm)
            - Qspr: Spring total (Mar-May) (mm)
            - Qsum: Summer total (Jun-Aug) (mm)
            - Qfal: Fall total (Sep-Nov) (mm)
            - Q1, Q5, Q10, ..., Q99: Flow percentiles (mm/day)
            - Q95_Q10: High-low flow difference (mm/day)

    Notes
    -----
    - Seasons are defined as:
        - Winter: December, January, February
        - Spring: March, April, May
        - Summer: June, July, August
        - Fall: September, October, November
    - Water years with fewer than min_nona_days non-NA days are excluded
    - Percentiles are calculated per water year, then statistics are
      computed across years
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "month", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Filter to water years with sufficient data
    year_counts = df.groupby("water_year")["Q"].apply(lambda x: x.notna().sum())
    valid_years = year_counts[year_counts >= min_nona_days].index
    df = df[df["water_year"].isin(valid_years)]

    if len(df) == 0:
        # Return empty stats
        return generate_stats(pd.DataFrame({"water_year": []}), value_cols=[])

    # Calculate annual totals (sum of daily Q in mm/day -> total mm)
    annual_totals = df.groupby("water_year")["Q"].sum().reset_index()
    annual_totals.columns = ["water_year", "Qann"]

    # Calculate seasonal totals
    # Winter: Dec (12), Jan (1), Feb (2)
    winter = df[df["month"].isin([12, 1, 2])]
    winter_totals = winter.groupby("water_year")["Q"].sum().reset_index()
    winter_totals.columns = ["water_year", "Qwin"]

    # Spring: Mar (3), Apr (4), May (5)
    spring = df[df["month"].isin([3, 4, 5])]
    spring_totals = spring.groupby("water_year")["Q"].sum().reset_index()
    spring_totals.columns = ["water_year", "Qspr"]

    # Summer: Jun (6), Jul (7), Aug (8)
    summer = df[df["month"].isin([6, 7, 8])]
    summer_totals = summer.groupby("water_year")["Q"].sum().reset_index()
    summer_totals.columns = ["water_year", "Qsum"]

    # Fall: Sep (9), Oct (10), Nov (11)
    fall = df[df["month"].isin([9, 10, 11])]
    fall_totals = fall.groupby("water_year")["Q"].sum().reset_index()
    fall_totals.columns = ["water_year", "Qfal"]

    # Calculate percentiles by water year (from config)
    percentiles = FLOW_PERCENTILES
    percentile_dfs = {}

    for p in percentiles:
        pct = df.groupby("water_year")["Q"].quantile(p / 100).reset_index()
        pct.columns = ["water_year", f"Q{p}"]
        percentile_dfs[p] = pct

    # Merge all metrics into single DataFrame
    all_metrics = annual_totals.copy()

    # Merge seasonal totals
    for seasonal_df in [winter_totals, spring_totals, summer_totals, fall_totals]:
        if len(seasonal_df) > 0:
            all_metrics = all_metrics.merge(seasonal_df, on="water_year", how="left")

    # Merge percentiles
    for p in percentiles:
        if len(percentile_dfs[p]) > 0:
            all_metrics = all_metrics.merge(
                percentile_dfs[p], on="water_year", how="left"
            )

    # Calculate Q95_Q10 difference
    if "Q95" in all_metrics.columns and "Q10" in all_metrics.columns:
        all_metrics["Q95_Q10"] = all_metrics["Q95"] - all_metrics["Q10"]
    else:
        all_metrics["Q95_Q10"] = np.nan

    # Get metric columns (all except water_year)
    metric_cols = [c for c in all_metrics.columns if c != "water_year"]

    # Generate statistics for all metrics
    return generate_stats(all_metrics, value_cols=metric_cols, year_col="water_year")

"""
Flow volume signatures.

Calculates annual and seasonal flow totals and percentiles.
"""

from typing import Dict, Optional, List
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import FLOW_PERCENTILES


def calculate_flow_vols_by_year(
    streamflow_data: pd.DataFrame,
    seasonal_flags: Optional[List[Dict]] = None,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate annual and seasonal flow volumes and percentiles.

    This function calculates 21 flow volume metrics for each water year,
    then applies generate_stats() to produce 8 statistics per metric
    (168 total output columns).

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - month: int, calendar month (1-12)
            - dowy: int, day of water year (1-366)
    seasonal_flags : list of dict, optional
        Per-year seasonal completeness flags from preprocess_daily_data().
        Each dict has keys: water_year, winter_complete, spring_complete,
        summer_complete, fall_complete (bool). When provided, seasonal totals
        (Qwin, Qspr, Qsum, Qfal) are set to NaN for years where the
        corresponding season is flagged incomplete.

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
    - Percentiles are calculated per water year, then statistics are
      computed across years
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "month", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

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

    # Apply seasonal completeness flags: set incomplete seasons to NaN
    if seasonal_flags is not None:
        flags_df = pd.DataFrame(seasonal_flags)
        season_col_map = {
            "winter_complete": "Qwin",
            "spring_complete": "Qspr",
            "summer_complete": "Qsum",
            "fall_complete": "Qfal",
        }
        for flag_col, q_col in season_col_map.items():
            if flag_col in flags_df.columns and q_col in all_metrics.columns:
                # Find water years where the season is incomplete
                incomplete_wys = flags_df.loc[
                    ~flags_df[flag_col].astype(bool), "water_year"
                ].values
                if len(incomplete_wys) > 0:
                    mask = all_metrics["water_year"].isin(incomplete_wys)
                    # Need to convert column to float to support NaN
                    all_metrics[q_col] = all_metrics[q_col].astype(float)
                    all_metrics.loc[mask, q_col] = np.nan

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
    return generate_stats(
        all_metrics,
        value_cols=metric_cols,
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )

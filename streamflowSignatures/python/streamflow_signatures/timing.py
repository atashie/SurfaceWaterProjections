"""
Flow timing signatures.

Calculates when during the water year cumulative flow reaches
various percentile thresholds (center of mass timing).
"""

from typing import Dict, Optional
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import D_PERCENTILES


def analyze_flow_timing_trends(
    streamflow_data: pd.DataFrame,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate flow timing metrics trends.

    Determines the day of water year when cumulative flow reaches
    various percentile thresholds (1%, 5%, 10%, ..., 95%, 99%) —
    15 metrics: 13 D-day percentiles (config timing.d_percentiles)
    plus D25_to_D75 and Dmax.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - dowy: int, day of water year (1-366)

    Returns
    -------
    dict
        Dictionary of signature statistics with keys for each metric:
            - D1_day, D5_day, ..., D99_day: Day when cumulative
              flow reaches each percentile
            - D25_to_D75: Duration between 25% and 75% cumulative flow
            - Dmax: Day of maximum discharge

    Notes
    -----
    - Day 1 = October 1 (start of water year)
    - Day 365/366 = September 30 (end of water year)
    - Years with any remaining NAs after coercion are skipped
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Get unique years
    years = df["water_year"].unique()

    # Define percentiles to calculate (from config)
    percentiles = D_PERCENTILES

    # Initialize results DataFrame
    columns = ["water_year"] + [f"D{p}_day" for p in percentiles] + ["D25_to_D75", "Dmax"]
    timing_by_year = pd.DataFrame(index=range(len(years)), columns=columns)
    timing_by_year["water_year"] = years

    # Build a lookup from water_year to row index for fast assignment
    yr_to_idx = {yr: i for i, yr in enumerate(years)}

    # Calculate timing metrics for each year
    for yr, year_data in df.groupby("water_year", sort=False):

        idx = yr_to_idx[yr]

        # Sort by day of water year
        year_data = year_data.sort_values("dowy")

        # Skip years with any NAs in Q
        q_vals = year_data["Q"].values
        if np.any(np.isnan(q_vals)):
            continue

        # Calculate total annual flow
        total_flow = np.sum(q_vals)

        # Skip years with zero or negative total flow
        if total_flow <= 0:
            continue

        cum_flow = np.cumsum(q_vals)

        # Calculate cumulative flow as percentage of total
        cum_pct = (cum_flow / total_flow) * 100

        dowy_values = year_data["dowy"].values

        # Find day when cumulative flow reaches each percentile
        for p in percentiles:
            above_threshold = np.where(cum_pct >= p)[0]
            if len(above_threshold) > 0:
                timing_by_year.loc[idx, f"D{p}_day"] = dowy_values[above_threshold[0]]

        # Calculate D25_to_D75 (days between 25% and 75% cumulative flow)
        above_25 = np.where(cum_pct >= 25)[0]
        above_75 = np.where(cum_pct >= 75)[0]

        if len(above_25) > 0 and len(above_75) > 0:
            day_25 = dowy_values[above_25[0]]
            day_75 = dowy_values[above_75[0]]
            timing_by_year.loc[idx, "D25_to_D75"] = day_75 - day_25

        # Calculate Dmax (day of maximum discharge)
        q_values = year_data["Q"].values
        if not np.all(np.isnan(q_values)):
            max_idx = np.nanargmax(q_values)
            timing_by_year.loc[idx, "Dmax"] = dowy_values[max_idx]

    # Convert columns to numeric
    metric_cols = [f"D{p}_day" for p in percentiles] + ["D25_to_D75", "Dmax"]
    for col in metric_cols:
        timing_by_year[col] = pd.to_numeric(timing_by_year[col], errors="coerce")
    timing_by_year["water_year"] = pd.to_numeric(timing_by_year["water_year"])

    # Generate statistics
    return generate_stats(
        timing_by_year,
        value_cols=metric_cols,
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )

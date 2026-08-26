"""
Runoff ratio signatures.

Calculates annual and seasonal runoff ratios (Q/P).
Requires precipitation (PPT) data.
"""

from typing import Dict, List, Optional
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import RUNOFF_MIN_ANNUAL_PPT, RUNOFF_MIN_SEASONAL_PPT


def analyze_Q_PPT_relationships(
    streamflow_data: pd.DataFrame,
    seasonal_flags: Optional[List[Dict]] = None,
    min_annual_ppt: float = RUNOFF_MIN_ANNUAL_PPT,
    min_seasonal_ppt: float = RUNOFF_MIN_SEASONAL_PPT,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate runoff ratio (Q/P) trends.

    Runoff ratio is the fraction of precipitation that becomes streamflow.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow and climate data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - PPT: float, daily precipitation in mm/day
            - month: int, calendar month (1-12)
    min_annual_ppt : float, default 10
        Minimum annual precipitation (mm) for valid ratio
    min_seasonal_ppt : float, default 1
        Minimum seasonal precipitation (mm) for valid ratio

    Returns
    -------
    dict
        Dictionary of signature statistics with keys:
            - annual_runoff_ratio_*: Annual Q/P statistics
            - winter_runoff_ratio_*: Winter (Dec-Feb) Q/P statistics
            - spring_runoff_ratio_*: Spring (Mar-May) Q/P statistics
            - summer_runoff_ratio_*: Summer (Jun-Aug) Q/P statistics
            - fall_runoff_ratio_*: Fall (Sep-Nov) Q/P statistics

    Notes
    -----
    - Ratios >1 indicate more runoff than precipitation (e.g., snowmelt)
    - Ratios <1 indicate losses to ET, infiltration, etc.
    - Years/seasons with insufficient precipitation return NA
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "PPT", "month"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Drop rows where either Q or PPT is NaN before summing
    # Matches R's aggregate(cbind(Q, PPT) ~ water_year, na.action=na.omit)
    df_valid = df.dropna(subset=["Q", "PPT"])

    # Calculate annual totals
    annual_totals = df_valid.groupby("water_year").agg({
        "Q": "sum",
        "PPT": "sum"
    }).reset_index()

    # Calculate annual runoff ratio with minimum PPT threshold
    annual_totals["annual_runoff_ratio"] = np.where(
        annual_totals["PPT"] > min_annual_ppt,
        annual_totals["Q"] / annual_totals["PPT"],
        np.nan
    )

    # Calculate seasonal totals and ratios
    def calc_seasonal_ratio(months, season_name):
        """Calculate seasonal runoff ratio."""
        seasonal = df_valid[df_valid["month"].isin(months)].copy()
        if len(seasonal) == 0:
            return pd.DataFrame({
                "water_year": annual_totals["water_year"],
                f"{season_name}_runoff_ratio": np.nan
            })

        totals = seasonal.groupby("water_year").agg({
            "Q": "sum",
            "PPT": "sum"
        }).reset_index()

        totals[f"{season_name}_runoff_ratio"] = np.where(
            totals["PPT"] > min_seasonal_ppt,
            totals["Q"] / totals["PPT"],
            np.nan
        )

        return totals[["water_year", f"{season_name}_runoff_ratio"]]

    # Calculate seasonal ratios
    winter_ratios = calc_seasonal_ratio([12, 1, 2], "winter")
    spring_ratios = calc_seasonal_ratio([3, 4, 5], "spring")
    summer_ratios = calc_seasonal_ratio([6, 7, 8], "summer")
    fall_ratios = calc_seasonal_ratio([9, 10, 11], "fall")

    # Merge all ratios
    all_ratios = annual_totals[["water_year", "annual_runoff_ratio"]].copy()

    for seasonal_df in [winter_ratios, spring_ratios, summer_ratios, fall_ratios]:
        all_ratios = all_ratios.merge(seasonal_df, on="water_year", how="left")

    # Apply seasonal completeness flags: set incomplete seasons to NA
    if seasonal_flags is not None and len(seasonal_flags) > 0:
        ratio_season_map = {
            "winter_runoff_ratio": "winter_complete",
            "spring_runoff_ratio": "spring_complete",
            "summer_runoff_ratio": "summer_complete",
            "fall_runoff_ratio": "fall_complete",
        }
        flags_df = pd.DataFrame(seasonal_flags)
        for ratio_col, flag_col in ratio_season_map.items():
            if ratio_col in all_ratios.columns and flag_col in flags_df.columns:
                merged = all_ratios[["water_year"]].merge(
                    flags_df[["water_year", flag_col]], on="water_year", how="left"
                )
                incomplete_mask = merged[flag_col].eq(False)
                all_ratios.loc[incomplete_mask, ratio_col] = np.nan

    # Get metric columns
    metric_cols = [c for c in all_ratios.columns if c != "water_year"]

    # Generate statistics
    result = generate_stats(
        all_ratios,
        value_cols=metric_cols,
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )

    # Count years where annual runoff ratio > 2.0 (per-gage scalar diagnostic)
    arr = all_ratios["annual_runoff_ratio"].values
    n_high = int(np.nansum(arr > 2.0))
    result["runoff_ratio_high_count"] = float(n_high)

    return result

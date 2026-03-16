"""
Streamflow elasticity signatures.

Measures how sensitively streamflow responds to precipitation changes.
Requires precipitation (PPT) data.
"""

from typing import Dict, Optional
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import (
    ELASTICITY_WINDOW_YEARS,
    ELASTICITY_MIN_YEARS,
    ELASTICITY_MIN_DATA_COMPLETENESS,
    ELASTICITY_MIN_ANNUAL_PPT,
)


def calculate_streamflow_elasticity(
    streamflow_data: pd.DataFrame,
    rolling_window: int = ELASTICITY_WINDOW_YEARS,
    min_years: int = ELASTICITY_MIN_YEARS,
    min_completeness: float = ELASTICITY_MIN_DATA_COMPLETENESS,
    min_annual_ppt: float = ELASTICITY_MIN_ANNUAL_PPT,
) -> Dict[str, float]:
    """
    Calculate streamflow elasticity trends.

    Elasticity measures how streamflow changes in response to
    precipitation changes. E = (dQ/dP) / (Q_mean/P_mean)

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow and climate data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - PPT: float, daily precipitation in mm/day
    rolling_window : int, default 11
        Size of rolling window for trend calculation (years)
    min_years : int, default 15
        Minimum number of valid years required
    min_completeness : float, default 0.9
        Minimum fraction of valid days per year
    min_annual_ppt : float, default 10
        Minimum annual precipitation (mm) for valid year

    Returns
    -------
    dict
        Dictionary of signature statistics with keys:
            - elasticity_static: Overall median elasticity (single value)
            - elasticity_senn_slp, elasticity_linear_slp, etc.: Trend statistics

    Notes
    -----
    - E ~ 1.0: Proportional response (10% more P -> 10% more Q)
    - E > 1.0: Amplified response (more sensitive to P changes)
    - E < 1.0: Dampened response (less sensitive to P changes)

    Reference:
        Sawicz, K., et al. (2011). Catchment classification:
        empirical analysis of hydrologic similarity.
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "PPT"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Aggregate to annual totals with data quality check
    annual = df.groupby("water_year").agg({
        "Q": lambda x: x.sum(),
        "PPT": lambda x: x.sum(),
    }).reset_index()
    annual.columns = ["water_year", "Q_annual", "P_annual"]

    # Count valid days per year
    valid_counts = df.groupby("water_year").apply(
        lambda x: ((~x["Q"].isna()) & (~x["PPT"].isna())).sum()
    ).reset_index()
    valid_counts.columns = ["water_year", "n_valid_days"]
    annual = annual.merge(valid_counts, on="water_year")

    # Calculate expected days per water year (accounting for leap years)
    def expected_days(wy):
        # Water year ends in year wy, so check if wy is a leap year
        is_leap = (wy % 4 == 0 and wy % 100 != 0) or (wy % 400 == 0)
        return 366 if is_leap else 365

    annual["expected_days"] = annual["water_year"].apply(expected_days)

    # Filter: remove years with >10% missing data
    annual = annual[annual["n_valid_days"] >= min_completeness * annual["expected_days"]]

    # Remove years with zero or very low precipitation
    annual = annual[annual["P_annual"] > min_annual_ppt]

    # Return NAs if insufficient years
    if len(annual) < min_years:
        return {
            "elasticity_static": np.nan,
            "elasticity_senn_slp": np.nan,
            "elasticity_linear_slp": np.nan,
            "elasticity_spearman_rho": np.nan,
            "elasticity_spearman_pval": np.nan,
            "elasticity_mk_rho": np.nan,
            "elasticity_mk_pval": np.nan,
            "elasticity_mean": np.nan,
            "elasticity_median": np.nan,
        }

    # Calculate long-term means
    Q_mean = annual["Q_annual"].mean()
    P_mean = annual["P_annual"].mean()

    # Calculate annual elasticity values
    # E_i = (dQ_i/dP_i) / (Q_mean/P_mean)
    annual["dQ"] = annual["Q_annual"] - Q_mean
    annual["dP"] = annual["P_annual"] - P_mean

    # Avoid division by zero
    annual["elasticity"] = np.where(
        np.abs(annual["dP"]) > 0.1,
        (annual["dQ"] / annual["dP"]) / (Q_mean / P_mean),
        np.nan
    )

    # Static elasticity is the median
    elasticity_static = annual["elasticity"].median()

    # Calculate rolling window elasticity
    if rolling_window is not None and len(annual) >= rolling_window:
        annual = annual.sort_values("water_year")

        rolling_elasticity = []
        rolling_years = []

        for end_idx in range(rolling_window - 1, len(annual)):
            start_idx = end_idx - rolling_window + 1
            window = annual.iloc[start_idx:end_idx + 1]

            Q_mean_w = window["Q_annual"].mean()
            P_mean_w = window["P_annual"].mean()

            window = window.copy()
            window["dQ_w"] = window["Q_annual"] - Q_mean_w
            window["dP_w"] = window["P_annual"] - P_mean_w

            window["e_w"] = np.where(
                np.abs(window["dP_w"]) > 0.1,
                (window["dQ_w"] / window["dP_w"]) / (Q_mean_w / P_mean_w),
                np.nan
            )

            rolling_elasticity.append(window["e_w"].median())
            rolling_years.append(annual.iloc[end_idx]["water_year"])

        rolling_df = pd.DataFrame({
            "water_year": rolling_years,
            "elasticity_rolling": rolling_elasticity
        })

        # Calculate trend statistics on rolling elasticity
        trend_stats = generate_stats(
            rolling_df,
            value_cols=["elasticity_rolling"],
            year_col="water_year"
        )

        result = {
            "elasticity_static": elasticity_static,
            "elasticity_senn_slp": trend_stats.get("elasticity_rolling_senn_slp", np.nan),
            "elasticity_linear_slp": trend_stats.get("elasticity_rolling_linear_slp", np.nan),
            "elasticity_spearman_rho": trend_stats.get("elasticity_rolling_spearman_rho", np.nan),
            "elasticity_spearman_pval": trend_stats.get("elasticity_rolling_spearman_pval", np.nan),
            "elasticity_mk_rho": trend_stats.get("elasticity_rolling_mk_rho", np.nan),
            "elasticity_mk_pval": trend_stats.get("elasticity_rolling_mk_pval", np.nan),
            "elasticity_mean": trend_stats.get("elasticity_rolling_mean", np.nan),
            "elasticity_median": trend_stats.get("elasticity_rolling_median", np.nan),
        }
    else:
        # No rolling window - use annual values for stats
        annual_valid = annual[~annual["elasticity"].isna()]
        trend_stats = generate_stats(
            annual_valid,
            value_cols=["elasticity"],
            year_col="water_year"
        )

        result = {
            "elasticity_static": elasticity_static,
            "elasticity_senn_slp": trend_stats.get("elasticity_senn_slp", np.nan),
            "elasticity_linear_slp": trend_stats.get("elasticity_linear_slp", np.nan),
            "elasticity_spearman_rho": trend_stats.get("elasticity_spearman_rho", np.nan),
            "elasticity_spearman_pval": trend_stats.get("elasticity_spearman_pval", np.nan),
            "elasticity_mk_rho": trend_stats.get("elasticity_mk_rho", np.nan),
            "elasticity_mk_pval": trend_stats.get("elasticity_mk_pval", np.nan),
            "elasticity_mean": trend_stats.get("elasticity_mean", np.nan),
            "elasticity_median": trend_stats.get("elasticity_median", np.nan),
        }

    return result

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
    ELASTICITY_MIN_ANNUAL_PPT,
)

_SUFFIXES = [
    "_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
    "_mk_rho", "_mk_pval", "_mean", "_median",
]


def _empty_elasticity() -> Dict[str, float]:
    """Return NaN dict for all elasticity outputs."""
    result = {"elasticity_static": np.nan}
    for prefix in ("elasticity_rolling", "elasticity_annual"):
        for suffix in _SUFFIXES:
            result[f"{prefix}{suffix}"] = np.nan
    result["elasticity_years_total"] = np.nan
    result["elasticity_years_low_ppt"] = np.nan
    return result


def calculate_streamflow_elasticity(
    streamflow_data: pd.DataFrame,
    rolling_window: int = ELASTICITY_WINDOW_YEARS,
    min_years: int = ELASTICITY_MIN_YEARS,
    min_annual_ppt: float = ELASTICITY_MIN_ANNUAL_PPT,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate streamflow elasticity trends.

    Calculates:
    - elasticity_static: Overall catchment elasticity (single value)
    - elasticity_rolling: Rolling window (11-year) elasticity trends (8 statistics)
      Uses departure-from-mean method (Sawicz et al. 2011)
    - elasticity_annual: Year-over-year elasticity trends (8 statistics)
      Uses consecutive-year differences for year-to-year sensitivity
    - elasticity_years_total: Total years available (diagnostic scalar)
    - elasticity_years_low_ppt: Years excluded due to low PPT (diagnostic scalar)

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
    min_annual_ppt : float, default 10
        Minimum annual precipitation (mm) for valid year

    Returns
    -------
    dict
        Dictionary with elasticity_static, elasticity_rolling_*, elasticity_annual_*
        statistics, and diagnostic scalars.

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

    # Aggregate to annual totals
    annual = df.groupby("water_year").agg({
        "Q": lambda x: x.sum(),
        "PPT": lambda x: x.sum(),
    }).reset_index()
    annual.columns = ["water_year", "Q_annual", "P_annual"]

    # Track totals BEFORE filtering (Step 5: diagnostics)
    elasticity_years_total = float(len(annual))

    # Remove years with zero or very low precipitation
    annual = annual[annual["P_annual"] > min_annual_ppt]
    elasticity_years_low_ppt = elasticity_years_total - float(len(annual))

    # Return NAs if insufficient years
    if len(annual) < min_years:
        result = _empty_elasticity()
        result["elasticity_years_total"] = elasticity_years_total
        result["elasticity_years_low_ppt"] = elasticity_years_low_ppt
        return result

    # Calculate long-term means
    Q_mean = annual["Q_annual"].mean()
    P_mean = annual["P_annual"].mean()

    if P_mean <= 0 or Q_mean <= 0:
        result = _empty_elasticity()
        result["elasticity_years_total"] = elasticity_years_total
        result["elasticity_years_low_ppt"] = elasticity_years_low_ppt
        return result

    # Calculate annual elasticity values (departure-from-mean)
    annual = annual.copy()
    annual["dQ"] = annual["Q_annual"] - Q_mean
    annual["dP"] = annual["P_annual"] - P_mean

    annual["elasticity"] = np.where(
        np.abs(annual["dP"]) > 0.1,
        (annual["dQ"] / annual["dP"]) / (Q_mean / P_mean),
        np.nan
    )

    # Static elasticity is the median
    elasticity_static = annual["elasticity"].median()

    result = {
        "elasticity_static": elasticity_static,
        "elasticity_years_total": elasticity_years_total,
        "elasticity_years_low_ppt": elasticity_years_low_ppt,
    }

    # --- Rolling window elasticity (Sawicz et al. 2011) ---
    annual = annual.sort_values("water_year")

    if rolling_window is not None and len(annual) >= rolling_window + 3:
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

        if len(rolling_elasticity) >= 3:
            rolling_df = pd.DataFrame({
                "water_year": rolling_years,
                "elasticity_rolling": rolling_elasticity
            })

            rolling_stats = generate_stats(
                rolling_df,
                value_cols=["elasticity_rolling"],
                year_col="water_year",
                trend_completeness=trend_completeness,
                decade_completeness=decade_completeness, changepoint=changepoint, collector=collector,
            )
            result.update(rolling_stats)
        else:
            for suffix in _SUFFIXES:
                result[f"elasticity_rolling{suffix}"] = np.nan
    else:
        for suffix in _SUFFIXES:
            result[f"elasticity_rolling{suffix}"] = np.nan

    # --- Year-over-year elasticity (consecutive-year differences) ---
    # E_t = ((Q_t - Q_{t-1}) / (P_t - P_{t-1})) / (Q_mean / P_mean)
    Q_valid = annual["Q_annual"].values
    P_valid = annual["P_annual"].values
    years_valid = annual["water_year"].values

    yoy_elasticity = []
    yoy_years = []
    for i in range(1, len(Q_valid)):
        dQ = Q_valid[i] - Q_valid[i - 1]
        dP = P_valid[i] - P_valid[i - 1]
        if abs(dP) > 0.1:
            e = (dQ / dP) / (Q_mean / P_mean)
            yoy_elasticity.append(e)
            yoy_years.append(int(years_valid[i]))

    if len(yoy_elasticity) >= 3:
        yoy_df = pd.DataFrame({
            "water_year": yoy_years,
            "elasticity_annual": yoy_elasticity
        })
        yoy_stats = generate_stats(
            yoy_df,
            value_cols=["elasticity_annual"],
            year_col="water_year",
            trend_completeness=trend_completeness,
            decade_completeness=decade_completeness, changepoint=changepoint, collector=collector,
        )
        result.update(yoy_stats)
    else:
        for suffix in _SUFFIXES:
            result[f"elasticity_annual{suffix}"] = np.nan

    return result

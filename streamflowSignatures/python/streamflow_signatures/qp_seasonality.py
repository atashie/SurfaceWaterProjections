"""
Q-P seasonality signatures.

Quantifies seasonality in the cumulative streamflow vs. cumulative
precipitation relationship. Requires precipitation (PPT) data.
"""

from typing import Dict
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import (
    QP_SLOPE_WINDOW_DAYS,
    QP_MIN_YEARS,
    QP_MIN_DAYS_PER_YEAR,
    QP_MAX_NA_FRAC,
)


def calculate_qp_seasonality(
    streamflow_data: pd.DataFrame,
    slope_window_days: int = QP_SLOPE_WINDOW_DAYS,
    min_years: int = QP_MIN_YEARS,
    min_days_per_year: int = QP_MIN_DAYS_PER_YEAR,
    max_na_frac: float = QP_MAX_NA_FRAC,
) -> Dict[str, float]:
    """
    Calculate Q-P seasonality metrics.

    Quantifies how the relationship between cumulative streamflow and
    cumulative precipitation varies seasonally.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow and climate data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - PPT: float, daily precipitation in mm/day
            - month: int, calendar month (1-12)
            - dowy: int, day of water year
    slope_window_days : int, default 30
        Rolling window size for slope calculation
    min_years : int, default 10
        Minimum number of valid years required
    min_days_per_year : int, default 300
        Minimum days required per water year
    max_na_frac : float, default 0.1
        Maximum fraction of NA values allowed per year

    Returns
    -------
    dict
        Dictionary of signature statistics with keys:
            - qp_slope_sd_*: Standard deviation of monthly Q-P slopes
            - qp_bimodality_*: Bimodality coefficient of slope distribution

    Notes
    -----
    - qp_slope_sd: Higher values indicate stronger seasonal variation
    - qp_bimodality: Values > 0.555 suggest bimodal (seasonal) patterns

    Reference:
        Wrede, S., et al. (2015). Towards a common classification
        framework for hydrological models.
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "PPT", "month", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data
    years = df["water_year"].unique()

    # Return NAs if insufficient years
    if len(years) < min_years:
        return _empty_result()

    # Calculate annual metrics
    annual_metrics = []

    for yr, year_data in df.groupby("water_year", sort=False):
        year_data = year_data.copy()  # needed because we mutate below (fillna)

        if len(year_data) < min_days_per_year:
            continue

        # Check for too many NA values
        na_frac_Q = year_data["Q"].isna().sum() / len(year_data)
        na_frac_PPT = year_data["PPT"].isna().sum() / len(year_data)
        if na_frac_Q > max_na_frac or na_frac_PPT > max_na_frac:
            continue

        # Replace remaining NA with 0 for cumsum
        year_data["Q"] = year_data["Q"].fillna(0)
        year_data["PPT"] = year_data["PPT"].fillna(0)

        # Sort by day of water year
        year_data = year_data.sort_values("dowy")

        # Calculate cumulative Q and P
        cum_Q = year_data["Q"].cumsum().values
        cum_P = year_data["PPT"].cumsum().values
        months = year_data["month"].values

        n = len(year_data)
        if n < slope_window_days:
            continue

        # Calculate rolling slope of cum_Q vs cum_P
        slopes = []
        slope_months = []

        for end_idx in range(slope_window_days - 1, n):
            start_idx = end_idx - slope_window_days + 1

            delta_P = cum_P[end_idx] - cum_P[start_idx]
            delta_Q = cum_Q[end_idx] - cum_Q[start_idx]

            if abs(delta_P) < 0.01:
                slopes.append(np.nan)
            else:
                slopes.append(delta_Q / delta_P)

            # Associate with middle of window (match R: end_idx - floor(window/2))
            mid_idx = end_idx - slope_window_days // 2
            slope_months.append(months[mid_idx])

        slopes = np.array(slopes)
        slope_months = np.array(slope_months)

        # Calculate monthly mean slopes
        monthly_means = {}
        for m in range(1, 13):
            month_slopes = slopes[slope_months == m]
            valid_slopes = month_slopes[~np.isnan(month_slopes)]
            if len(valid_slopes) > 0:
                monthly_means[m] = np.mean(valid_slopes)

        if len(monthly_means) < 6:
            continue

        # Metric 1: SD of monthly slopes (ddof=1 to match R's sd())
        qp_slope_sd = np.std(list(monthly_means.values()), ddof=1)

        # Metric 2: Bimodality coefficient
        # B = (skewness^2 + 1) / kurtosis
        # Values > 5/9 (0.555) suggest bimodality
        slopes_clean = slopes[~np.isnan(slopes)]
        if len(slopes_clean) < 30:
            qp_bimodality = np.nan
        else:
            n_s = len(slopes_clean)
            m = np.mean(slopes_clean)
            s = np.std(slopes_clean, ddof=1)

            if s < 1e-10:
                qp_bimodality = np.nan
            else:
                skewness = np.sum((slopes_clean - m)**3) / (n_s * s**3)
                kurtosis = np.sum((slopes_clean - m)**4) / (n_s * s**4)

                if kurtosis < 1e-10:
                    qp_bimodality = np.nan
                else:
                    qp_bimodality = (skewness**2 + 1) / kurtosis

        annual_metrics.append({
            "water_year": yr,
            "qp_slope_sd": qp_slope_sd,
            "qp_bimodality": qp_bimodality,
        })

    if len(annual_metrics) < 5:
        return _empty_result()

    annual_df = pd.DataFrame(annual_metrics)

    # Calculate trend statistics
    sd_stats = generate_stats(
        annual_df,
        value_cols=["qp_slope_sd"],
        year_col="water_year"
    )

    bi_df = annual_df[~annual_df["qp_bimodality"].isna()]
    if len(bi_df) >= 3:
        bi_stats = generate_stats(
            bi_df,
            value_cols=["qp_bimodality"],
            year_col="water_year"
        )
    else:
        bi_stats = {}

    result = {
        "qp_slope_sd_senn_slp": sd_stats.get("qp_slope_sd_senn_slp", np.nan),
        "qp_slope_sd_linear_slp": sd_stats.get("qp_slope_sd_linear_slp", np.nan),
        "qp_slope_sd_spearman_rho": sd_stats.get("qp_slope_sd_spearman_rho", np.nan),
        "qp_slope_sd_spearman_pval": sd_stats.get("qp_slope_sd_spearman_pval", np.nan),
        "qp_slope_sd_mk_rho": sd_stats.get("qp_slope_sd_mk_rho", np.nan),
        "qp_slope_sd_mk_pval": sd_stats.get("qp_slope_sd_mk_pval", np.nan),
        "qp_slope_sd_mean": sd_stats.get("qp_slope_sd_mean", np.nan),
        "qp_slope_sd_median": sd_stats.get("qp_slope_sd_median", np.nan),
        "qp_bimodality_senn_slp": bi_stats.get("qp_bimodality_senn_slp", np.nan),
        "qp_bimodality_linear_slp": bi_stats.get("qp_bimodality_linear_slp", np.nan),
        "qp_bimodality_spearman_rho": bi_stats.get("qp_bimodality_spearman_rho", np.nan),
        "qp_bimodality_spearman_pval": bi_stats.get("qp_bimodality_spearman_pval", np.nan),
        "qp_bimodality_mk_rho": bi_stats.get("qp_bimodality_mk_rho", np.nan),
        "qp_bimodality_mk_pval": bi_stats.get("qp_bimodality_mk_pval", np.nan),
        "qp_bimodality_mean": bi_stats.get("qp_bimodality_mean", np.nan),
        "qp_bimodality_median": bi_stats.get("qp_bimodality_median", np.nan),
    }

    return result


def _empty_result() -> Dict[str, float]:
    """Return empty result with all NAs."""
    return {
        "qp_slope_sd_senn_slp": np.nan,
        "qp_slope_sd_linear_slp": np.nan,
        "qp_slope_sd_spearman_rho": np.nan,
        "qp_slope_sd_spearman_pval": np.nan,
        "qp_slope_sd_mk_rho": np.nan,
        "qp_slope_sd_mk_pval": np.nan,
        "qp_slope_sd_mean": np.nan,
        "qp_slope_sd_median": np.nan,
        "qp_bimodality_senn_slp": np.nan,
        "qp_bimodality_linear_slp": np.nan,
        "qp_bimodality_spearman_rho": np.nan,
        "qp_bimodality_spearman_pval": np.nan,
        "qp_bimodality_mk_rho": np.nan,
        "qp_bimodality_mk_pval": np.nan,
        "qp_bimodality_mean": np.nan,
        "qp_bimodality_median": np.nan,
    }

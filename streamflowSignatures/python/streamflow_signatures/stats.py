"""
Core statistical functions for streamflow signature analysis.

This module provides the statistics engine that produces 8 standard statistics
for each signature metric, matching the R implementation's generate_stats().
"""

from typing import Optional, Dict, List, Union
import numpy as np
import pandas as pd
from scipy import stats


def theil_sen_slope(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate the Theil-Sen slope estimator.

    This is a robust non-parametric estimator of linear trend, calculated as
    the median of all pairwise slopes between points.

    Parameters
    ----------
    x : np.ndarray
        Independent variable (typically years)
    y : np.ndarray
        Dependent variable (metric values)

    Returns
    -------
    float
        Theil-Sen slope estimate, or np.nan if calculation fails

    Notes
    -----
    Matches R's zyp::zyp.sen implementation.
    """
    # Remove NaN pairs
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) < 2:
        return np.nan

    try:
        result = stats.theilslopes(y_clean, x_clean)
        return result.slope
    except Exception:
        return np.nan


def mann_kendall_test(y: np.ndarray) -> tuple:
    """
    Perform the Mann-Kendall trend test.

    The Mann-Kendall test is a non-parametric test for monotonic trend.

    Parameters
    ----------
    y : np.ndarray
        Time series values (assumed to be in temporal order)

    Returns
    -------
    tuple
        (tau, p_value) where tau is Kendall's tau and p_value is the
        two-sided p-value for the test. Returns (np.nan, np.nan) if
        calculation fails.

    Notes
    -----
    Matches R's Kendall::MannKendall implementation.
    The Mann-Kendall test is equivalent to Kendall's tau correlation
    between the values and their time indices.
    """
    # Remove NaN values
    y_clean = y[~np.isnan(y)]

    if len(y_clean) < 3:
        return np.nan, np.nan

    try:
        # Mann-Kendall is Kendall's tau between values and time indices
        x = np.arange(len(y_clean))
        tau, p_value = stats.kendalltau(x, y_clean)
        return tau, p_value
    except Exception:
        return np.nan, np.nan


def spearman_correlation(x: np.ndarray, y: np.ndarray) -> tuple:
    """
    Calculate Spearman rank correlation and p-value.

    Parameters
    ----------
    x : np.ndarray
        First variable (typically years)
    y : np.ndarray
        Second variable (metric values)

    Returns
    -------
    tuple
        (rho, p_value) where rho is Spearman's correlation coefficient
        and p_value is the two-sided p-value. Returns (np.nan, np.nan)
        if calculation fails.
    """
    # Remove NaN pairs
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) < 3:
        return np.nan, np.nan

    try:
        rho, p_value = stats.spearmanr(x_clean, y_clean)
        return rho, p_value
    except Exception:
        return np.nan, np.nan


def linear_slope(x: np.ndarray, y: np.ndarray) -> float:
    """
    Calculate ordinary least squares linear regression slope.

    Parameters
    ----------
    x : np.ndarray
        Independent variable (typically years)
    y : np.ndarray
        Dependent variable (metric values)

    Returns
    -------
    float
        Linear regression slope, or np.nan if calculation fails
    """
    # Remove NaN pairs
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean = x[mask]
    y_clean = y[mask]

    if len(x_clean) < 2:
        return np.nan

    try:
        result = stats.linregress(x_clean, y_clean)
        return result.slope
    except Exception:
        return np.nan


def generate_stats(
    data: pd.DataFrame,
    value_cols: Optional[List[str]] = None,
    year_col: str = "water_year",
    min_rows: int = 3
) -> Dict[str, float]:
    """
    Generate 8 standard statistics for each metric column.

    This is the core statistics engine that produces standardized trend
    and summary statistics for each signature metric.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing annual metric values with a year column
    value_cols : list of str, optional
        Columns to analyze. If None, uses all numeric columns except year_col
    year_col : str, default "water_year"
        Name of the year column
    min_rows : int, default 3
        Minimum number of non-NA rows required for statistics

    Returns
    -------
    dict
        Dictionary with keys for each statistic. For each metric column,
        8 statistics are produced with suffixes:
            - _senn_slp: Theil-Sen slope
            - _linear_slp: Linear regression slope
            - _spearman_rho: Spearman correlation coefficient
            - _spearman_pval: Spearman p-value
            - _mk_rho: Mann-Kendall tau
            - _mk_pval: Mann-Kendall p-value
            - _mean: Arithmetic mean
            - _median: Median

    Examples
    --------
    >>> data = pd.DataFrame({
    ...     'water_year': [2000, 2001, 2002, 2003, 2004],
    ...     'Qann': [100, 110, 105, 120, 115]
    ... })
    >>> stats = generate_stats(data, value_cols=['Qann'])
    >>> print(stats['Qann_mean'])
    110.0
    """
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Input 'data' must be a DataFrame")

    if year_col not in data.columns:
        raise ValueError(f"Year column '{year_col}' not found in data")

    # Sort by year for trend calculations
    data = data.sort_values(year_col).copy()

    # Determine value columns
    if value_cols is None:
        numeric_cols = data.select_dtypes(include=[np.number]).columns.tolist()
        value_cols = [c for c in numeric_cols if c != year_col]

    results = {}

    for col in value_cols:
        if col not in data.columns:
            # Skip missing columns with warning
            continue

        # Extract non-NA data
        mask = ~data[col].isna()
        years = data.loc[mask, year_col].values.astype(float)
        values = data.loc[mask, col].values.astype(float)

        # Check minimum data requirement
        if len(values) < min_rows:
            results[f"{col}_senn_slp"] = np.nan
            results[f"{col}_linear_slp"] = np.nan
            results[f"{col}_spearman_rho"] = np.nan
            results[f"{col}_spearman_pval"] = np.nan
            results[f"{col}_mk_rho"] = np.nan
            results[f"{col}_mk_pval"] = np.nan
            results[f"{col}_mean"] = np.nan
            results[f"{col}_median"] = np.nan
            continue

        # Calculate Theil-Sen slope
        results[f"{col}_senn_slp"] = theil_sen_slope(years, values)

        # Calculate linear regression slope
        results[f"{col}_linear_slp"] = linear_slope(years, values)

        # Calculate Spearman correlation
        rho, pval = spearman_correlation(years, values)
        results[f"{col}_spearman_rho"] = rho
        results[f"{col}_spearman_pval"] = pval

        # Calculate Mann-Kendall test
        mk_tau, mk_pval = mann_kendall_test(values)
        results[f"{col}_mk_rho"] = mk_tau
        results[f"{col}_mk_pval"] = mk_pval

        # Calculate mean and median
        results[f"{col}_mean"] = np.mean(values)
        results[f"{col}_median"] = np.median(values)

    return results

"""
Flashiness signature.

Calculates the Richards-Baker flashiness index, which measures
how quickly streamflow rises and falls.
"""

from typing import Dict
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import FLASHINESS_MIN_DAYS, FLASHINESS_MAX_MISSING_FRAC


def analyze_flashiness_trends(
    streamflow_data: pd.DataFrame,
    min_days: int = FLASHINESS_MIN_DAYS,
    max_missing_frac: float = FLASHINESS_MAX_MISSING_FRAC,
) -> Dict[str, float]:
    """
    Calculate Richards-Baker flashiness index trends.

    The R-B flashiness index is the sum of absolute day-to-day flow
    changes divided by total flow. Higher values indicate more
    "flashy" (rapidly changing) streamflow.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - dowy: int, day of water year (optional, for sorting)
    min_days : int, default 30
        Minimum number of non-NA days required per water year
    max_missing_frac : float, default 0.2
        Maximum fraction of missing values allowed (beyond this, year is skipped)

    Returns
    -------
    dict
        Dictionary of signature statistics with keys:
            - flashinessRB_senn_slp, flashinessRB_linear_slp, etc.

    Notes
    -----
    R-B Index formula:
        sum(|Q[i] - Q[i-1]|) / sum(Q)

    Reference:
        Baker, D.B., et al. (2004). A new flashiness index: characteristics
        and applications to midwestern rivers and streams.
    """
    # Validate required columns
    required_cols = ["water_year", "Q"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Collect results as list of dicts (P13: avoid pre-allocation + boolean .loc)
    results_list = []

    # Calculate R-B index for each year
    for yr, year_data in df.groupby("water_year", sort=False):

        # Check minimum data requirement
        non_na_count = year_data["Q"].notna().sum()
        if non_na_count < min_days:
            results_list.append({"water_year": yr, "RB_index": np.nan})
            continue

        # Sort by day of water year if available
        if "dowy" in year_data.columns:
            year_data = year_data.sort_values("dowy")

        # Get Q values
        q_values = year_data["Q"].values.copy()

        # Check missing value fraction
        na_frac = np.isnan(q_values).sum() / len(q_values)
        if na_frac > max_missing_frac:
            results_list.append({"water_year": yr, "RB_index": np.nan})
            continue

        # Interpolate missing values if some are present
        if np.any(np.isnan(q_values)):
            # Linear interpolation
            indices = np.arange(len(q_values))
            valid_mask = ~np.isnan(q_values)
            if valid_mask.sum() >= 2:
                q_values = np.interp(
                    indices,
                    indices[valid_mask],
                    q_values[valid_mask]
                )

        # Calculate absolute day-to-day changes
        q_diff = np.abs(np.diff(q_values))

        # Calculate total flow
        total_q = np.nansum(q_values)

        # Guard against division by zero
        if total_q == 0:
            results_list.append({"water_year": yr, "RB_index": np.nan})
            continue

        # Calculate R-B index
        rb_index = np.nansum(q_diff) / total_q

        # Store result
        results_list.append({"water_year": yr, "RB_index": rb_index})

    flashiness_by_year = pd.DataFrame(results_list)

    # Generate statistics
    result = generate_stats(
        flashiness_by_year,
        value_cols=["RB_index"],
        year_col="water_year"
    )

    # Rename columns to match expected output
    renamed = {}
    for key, value in result.items():
        new_key = key.replace("RB_index", "flashinessRB")
        renamed[new_key] = value

    return renamed

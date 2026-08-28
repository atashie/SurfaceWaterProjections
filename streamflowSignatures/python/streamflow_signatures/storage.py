"""
Average storage signature.

Calculates catchment storage using simplified water balance.
Requires precipitation (PPT) data.
"""

from typing import Dict, Optional
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import STORAGE_MIN_YEARS


def calculate_average_storage(
    streamflow_data: pd.DataFrame,
    min_years: int = STORAGE_MIN_YEARS,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate average catchment storage trends.

    Uses simplified water balance (dS = P - Q) without ET estimation.
    Calculates storage amount corresponding to average discharge.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow and climate data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - PPT: float, daily precipitation in mm/day
            - dowy: int, day of water year
    min_years : int, default 10
        Minimum number of valid years required

    Returns
    -------
    dict
        Dictionary of signature statistics with keys:
            - avg_storage_*: Average storage statistics (mm)

    Notes
    -----
    This calculation ignores evapotranspiration (ET), using only P - Q
    for the water balance. This simplification may overestimate storage
    in watersheds with significant ET losses.

    Reference:
        Peters, N.E., & Aulenbach, B.T. (2011). Water storage at
        the Panola Mountain Research Watershed.
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "PPT", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data
    years = df["water_year"].unique()

    # Return NAs if insufficient years
    if len(years) < min_years:
        return _empty_result()

    # Calculate storage for each year
    annual_storage = []

    for yr, year_data in df.groupby("water_year", sort=False):
        year_data = year_data.copy()  # defensive copy of the groupby slice (no in-place mutation follows)

        # Skip years with any remaining NAs (data should be pre-interpolated)
        if year_data["Q"].isna().any() or year_data["PPT"].isna().any():
            continue

        # Sort by day of water year
        year_data = year_data.sort_values("dowy")

        # Water balance: dS = P - Q
        dS = year_data["PPT"].values - year_data["Q"].values

        # Cumulative storage (relative to start of year)
        S = np.cumsum(dS)

        # Calculate mean Q for this year
        mean_Q = year_data["Q"].mean()

        # Find storage at mean Q using interpolation
        Q_values = year_data["Q"].values

        # Sort by Q for interpolation
        sort_idx = np.argsort(Q_values)
        Q_sorted = Q_values[sort_idx]
        S_sorted = S[sort_idx]

        # Remove NAs
        valid_mask = ~np.isnan(Q_sorted) & ~np.isnan(S_sorted)
        Q_valid = Q_sorted[valid_mask]
        S_valid = S_sorted[valid_mask]

        if len(Q_valid) < 10:
            continue

        # Average S values at duplicate Q points to match R's approx(ties=mean)
        # Without this, np.interp picks the last S value at tied Q points,
        # introducing noise that corrupts trend statistics
        Q_unique, inverse_idx = np.unique(Q_valid, return_inverse=True)
        S_means = np.zeros(len(Q_unique))
        counts = np.zeros(len(Q_unique))
        for idx_i, grp in enumerate(inverse_idx):
            S_means[grp] += S_valid[idx_i]
            counts[grp] += 1
        S_means /= counts

        if len(Q_unique) < 10:
            continue

        # Linear interpolation to find S at mean_Q
        try:
            avg_storage = np.interp(mean_Q, Q_unique, S_means)
        except Exception:
            avg_storage = np.nan

        if not np.isnan(avg_storage):
            annual_storage.append({
                "water_year": yr,
                "avg_storage": avg_storage
            })

    if len(annual_storage) < 5:
        return _empty_result()

    annual_df = pd.DataFrame(annual_storage)

    # Calculate trend statistics
    stats = generate_stats(
        annual_df,
        value_cols=["avg_storage"],
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )

    # Return the full generate_stats output (merged wholesale, as Julia does):
    # the previous explicit 8-key copy silently DROPPED the changepoint fields
    # (found 2026-08-24 by the Phase-1 schema gate: avg_storage was one of 3
    # bases missing its 8 _pettitt_* columns).
    return dict(stats)


def _empty_result() -> Dict[str, float]:
    """Return empty result with all NAs."""
    return {
        "avg_storage_senn_slp": np.nan,
        "avg_storage_linear_slp": np.nan,
        "avg_storage_spearman_rho": np.nan,
        "avg_storage_spearman_pval": np.nan,
        "avg_storage_mk_rho": np.nan,
        "avg_storage_mk_pval": np.nan,
        "avg_storage_mean": np.nan,
        "avg_storage_median": np.nan,
    }

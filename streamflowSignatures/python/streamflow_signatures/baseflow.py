"""
Baseflow signatures.

Calculates baseflow indices using recursive digital filters:
- Eckhardt filter (BFI_Eckhardt)
- Lyne-Hollick filter (BFI_LyneHollick)
"""

from typing import Dict, Optional
import numpy as np
import pandas as pd
from .stats import generate_stats
from .config import (
    ECKHARDT_BFIMAX,
    ECKHARDT_ALPHA,
    LYNE_HOLLICK_ALPHA,
    LYNE_HOLLICK_PASSES,
)


def eckhardt_filter(Q: np.ndarray, BFImax: float = ECKHARDT_BFIMAX, a: float = ECKHARDT_ALPHA) -> np.ndarray:
    """
    Apply Eckhardt recursive digital filter for baseflow separation.

    Parameters
    ----------
    Q : np.ndarray
        Daily discharge values
    BFImax : float, default 0.8
        Maximum baseflow index (fraction of flow that can be baseflow)
    a : float, default 0.98
        Filter parameter (recession constant)

    Returns
    -------
    np.ndarray
        Baseflow time series

    References
    ----------
    Eckhardt, K. (2005). How to construct recursive digital filters
    for baseflow separation.
    """
    n = len(Q)
    baseflow = np.zeros(n)

    # Initialize (matching Julia: min(BFImax * Q[1], Q[1]) if valid and > 0, else 0)
    if not np.isnan(Q[0]) and Q[0] > 0:
        baseflow[0] = min(BFImax * Q[0], Q[0])
    else:
        baseflow[0] = 0.0

    for i in range(1, n):
        if np.isnan(Q[i]):
            # Forward-fill baseflow on NaN Q (matching Julia — no cascade)
            baseflow[i] = baseflow[i - 1]
            continue
        # Normal computation — baseflow[i-1] is always valid due to forward-fill
        numerator = (1 - BFImax) * a * baseflow[i - 1] + (1 - a) * BFImax * Q[i]
        denominator = 1 - a * BFImax
        baseflow[i] = numerator / denominator
        # Baseflow cannot exceed total flow
        baseflow[i] = min(baseflow[i], Q[i])
        # Baseflow cannot be negative
        baseflow[i] = max(baseflow[i], 0)

    return baseflow


def lyne_hollick_filter(Q: np.ndarray, alpha: float = LYNE_HOLLICK_ALPHA, passes: int = LYNE_HOLLICK_PASSES) -> np.ndarray:
    """
    Apply Lyne-Hollick recursive digital filter for baseflow separation.

    Uses forward-backward filtering with multiple passes.

    Parameters
    ----------
    Q : np.ndarray
        Daily discharge values
    alpha : float, default 0.925
        Filter parameter
    passes : int, default 2
        Number of forward-backward passes

    Returns
    -------
    np.ndarray
        Baseflow time series

    References
    ----------
    Lyne, V., & Hollick, M. (1979). Stochastic time-variable
    rainfall-runoff modelling.
    """
    n = len(Q)
    input_signal = Q.copy()

    for pass_num in range(passes):
        # Forward pass
        qf_forward = np.zeros(n)
        qf_forward[0] = 0

        for i in range(1, n):
            if (np.isnan(input_signal[i]) or np.isnan(input_signal[i-1]) or
                np.isnan(qf_forward[i-1])):
                qf_forward[i] = np.nan
            else:
                qf_forward[i] = (alpha * qf_forward[i-1] +
                                 ((1 + alpha) / 2) * (input_signal[i] - input_signal[i-1]))
                # Constrain quickflow to [0, input_signal]
                qf_forward[i] = max(0, min(qf_forward[i], input_signal[i]))

        # Backward pass
        qf_backward = np.zeros(n)
        qf_backward[n-1] = qf_forward[n-1]

        for i in range(n-2, -1, -1):
            if (np.isnan(qf_forward[i]) or np.isnan(qf_forward[i+1]) or
                np.isnan(qf_backward[i+1])):
                qf_backward[i] = np.nan
            else:
                qf_backward[i] = (alpha * qf_backward[i+1] +
                                  ((1 + alpha) / 2) * (qf_forward[i] - qf_forward[i+1]))
                # Constrain quickflow to [0, input_signal]
                qf_backward[i] = max(0, min(qf_backward[i], input_signal[i]))

        # Baseflow from this pass becomes input for next pass
        input_signal = input_signal - qf_backward
        input_signal[input_signal < 0] = 0
        input_signal[np.isnan(Q)] = np.nan

    # After all passes, input_signal holds the final baseflow
    return input_signal


def analyze_baseflow_indices(
    streamflow_data: pd.DataFrame,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate baseflow index trends using recursive digital filters.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day
            - dowy: int, day of water year (for sorting)

    Returns
    -------
    dict
        Dictionary of signature statistics with keys:
            - BFI_Eckhardt_*: Statistics for Eckhardt baseflow index
            - BFI_LyneHollick_*: Statistics for Lyne-Hollick baseflow index

    Notes
    -----
    BFI = sum(baseflow) / sum(total_flow) for each water year.
    Expected relationship: BFI_Eckhardt < BFI_LyneHollick
    """
    # Validate required columns
    required_cols = ["water_year", "Q", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Initialize results DataFrame
    # Rows are built only for COMPUTABLE years (total valid Q > 0), matching
    # Julia's `total_Q <= 0 -> continue`: a pre-allocated all-years frame would
    # export NaN rows Julia omits and would make the `< 3 rows` guard count
    # non-computable years (found 2026-08-25, Phase-3 cross-language comparator).
    bfi_rows = []

    # Process each year
    for yr, year_data in df.groupby("water_year", sort=False):
        year_data = year_data.copy()  # needed because we mutate below (sort_values)

        # Sort by day of water year
        year_data = year_data.sort_values("dowy")

        # Get streamflow values
        Q = year_data["Q"].values.astype(float)

        # Apply Eckhardt filter
        baseflow_eckhardt = eckhardt_filter(Q)

        # Apply Lyne-Hollick filter
        baseflow_lyne = lyne_hollick_filter(Q)

        # Calculate annual totals and BFI using paired masking
        # The recursive filters can propagate NaN from adjacent missing values
        # to positions where Q itself is valid. Using paired masking (only sum
        # where BOTH Q and baseflow are valid) empirically matches R's behavior
        # for both Eckhardt and Lyne-Hollick filters.

        # Year-level computability gate (Julia: total_Q <= 0 -> skip the year)
        valid_q = ~np.isnan(Q)
        if valid_q.sum() == 0 or np.sum(Q[valid_q]) <= 0:
            continue

        bfi_eckhardt = np.nan
        valid_eck = valid_q & ~np.isnan(baseflow_eckhardt)
        if valid_eck.sum() > 0 and np.sum(Q[valid_eck]) > 0:
            bfi_eckhardt = np.sum(baseflow_eckhardt[valid_eck]) / np.sum(Q[valid_eck])

        bfi_lyne = np.nan
        valid_lh = valid_q & ~np.isnan(baseflow_lyne)
        if valid_lh.sum() > 0 and np.sum(Q[valid_lh]) > 0:
            bfi_lyne = np.sum(baseflow_lyne[valid_lh]) / np.sum(Q[valid_lh])

        bfi_rows.append({"water_year": yr, "BFI_Eckhardt": bfi_eckhardt,
                         "BFI_LyneHollick": bfi_lyne})

    bfi_by_year = pd.DataFrame(bfi_rows, columns=["water_year", "BFI_Eckhardt",
                                                  "BFI_LyneHollick"])

    # Generate statistics
    return generate_stats(
        bfi_by_year,
        value_cols=["BFI_Eckhardt", "BFI_LyneHollick"],
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )


def analyze_baseflow_indices_with_parameters(
    streamflow_data: pd.DataFrame,
    alpha: float,
    BFImax: float = ECKHARDT_BFIMAX,
    passes: int = LYNE_HOLLICK_PASSES,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate baseflow indices using recession-derived filter parameters.

    Same as analyze_baseflow_indices() but uses recession-derived alpha
    instead of fixed defaults. BFImax remains fixed at 0.8.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns: water_year, Q, dowy.
    alpha : float
        Recession-derived filter parameter (from recession_alpha_point_cloud_linear_reservoir).
        Must be in (0, 1). NaN or out-of-range returns all NaN stats.
    BFImax : float, default 0.8
        Maximum baseflow index for Eckhardt filter.
    passes : int, default 2
        Number of passes for Lyne-Hollick filter.

    Returns
    -------
    dict
        Dictionary with BFI_Eckhardt_param_* and BFI_LyneHollick_param_* statistics.
    """
    metrics = ["BFI_Eckhardt_param", "BFI_LyneHollick_param"]

    # Validate alpha
    if np.isnan(alpha) or alpha <= 0.0 or alpha >= 1.0:
        result = {}
        for m in metrics:
            for suffix in ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                           "_mk_rho", "_mk_pval", "_mean", "_median"]:
                result[f"{m}{suffix}"] = np.nan
        return result

    # Validate required columns
    required_cols = ["water_year", "Q", "dowy"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Computable years only (see the fixed-parameter function above)
    bfi_rows = []

    for yr, year_data in df.groupby("water_year", sort=False):
        year_data = year_data.copy()
        year_data = year_data.sort_values("dowy")
        Q = year_data["Q"].values.astype(float)

        # Apply filters with recession-derived alpha
        bf_eck = eckhardt_filter(Q, BFImax=BFImax, a=alpha)
        bf_lh = lyne_hollick_filter(Q, alpha=alpha, passes=passes)

        # Year-level computability gate (Julia: total_Q <= 0 -> skip the year)
        valid_q = ~np.isnan(Q)
        if valid_q.sum() == 0 or np.sum(Q[valid_q]) <= 0:
            continue

        bfi_eck = np.nan
        valid_eck = valid_q & ~np.isnan(bf_eck)
        if valid_eck.sum() > 0 and np.sum(Q[valid_eck]) > 0:
            bfi_eck = float(np.clip(np.sum(bf_eck[valid_eck]) / np.sum(Q[valid_eck]), 0.0, 1.0))

        bfi_lh = np.nan
        valid_lh = valid_q & ~np.isnan(bf_lh)
        if valid_lh.sum() > 0 and np.sum(Q[valid_lh]) > 0:
            bfi_lh = float(np.clip(np.sum(bf_lh[valid_lh]) / np.sum(Q[valid_lh]), 0.0, 1.0))

        bfi_rows.append({"water_year": yr, "BFI_Eckhardt_param": bfi_eck,
                         "BFI_LyneHollick_param": bfi_lh})

    bfi_by_year = pd.DataFrame(bfi_rows, columns=["water_year"] + metrics)

    return generate_stats(
        bfi_by_year,
        value_cols=metrics,
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )

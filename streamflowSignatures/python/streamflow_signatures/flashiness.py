"""
Flashiness signature.

Calculates the Richards-Baker flashiness index, which measures
how quickly streamflow rises and falls.
"""

from typing import Dict, Optional
import numpy as np
import pandas as pd
from .stats import generate_stats


def analyze_flashiness_trends(
    streamflow_data: pd.DataFrame,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
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

        # Sort by day of water year if available
        if "dowy" in year_data.columns:
            year_data = year_data.sort_values("dowy")

        # Get Q values (preprocessor guarantees no NAs in valid years)
        q_values = year_data["Q"].values.copy()

        # Calculate absolute day-to-day changes
        q_diff = np.abs(np.diff(q_values))

        # Calculate total flow
        total_q = np.nansum(q_values)

        # Non-computable year: zero or NEGATIVE total flow (negative-Q days are
        # retained by default; a negative denominator is meaningless). `<= 0`
        # matches canonical Julia — the `== 0` guard let WY totals < 0 produce
        # values Julia declines (found 2026-08-24 at gage 02244440 WY1999,
        # total Q = -78.8, via the Phase-1 floor-transition gate).
        if total_q <= 0:
            # Non-computable year: OMIT the row entirely (Julia pushes only when
            # the value is non-NaN). Absent row and NaN are semantically
            # equivalent per the annual-values contract, but omitting keeps the
            # exported parquet structurally identical across languages and keeps
            # the `< 3 rows` early return counting COMPUTABLE years, as Julia
            # does (found 2026-08-25 by the Phase-3 cross-language comparator).
            continue

        # Calculate R-B index
        rb_index = np.nansum(q_diff) / total_q

        # Store result
        results_list.append({"water_year": yr, "flashinessRB": rb_index})

    flashiness_by_year = pd.DataFrame(results_list)

    # Generate statistics
    # Metric columns use the CANONICAL names directly (matching Julia, which
    # never renames): the annual-values collector captures the column name it
    # is given, so a placeholder + post-hoc rename would put the wrong
    # signature name in the annual parquet (found 2026-08-25, Phase 3).
    result = generate_stats(
        flashiness_by_year,
        value_cols=["flashinessRB"],
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )

    return result

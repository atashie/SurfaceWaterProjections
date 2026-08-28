"""
Flow Duration Curve signatures.

Calculates slopes of the flow duration curve at different
exceedance probability ranges.
"""

from typing import Dict, Optional
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from .stats import generate_stats


def analyze_fdc_trends(
    streamflow_data: pd.DataFrame,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    min_values_for_stats: Optional[int] = None,
    changepoint: Optional[dict] = None,
    collector: Optional[object] = None,
) -> Dict[str, float]:
    """
    Calculate Flow Duration Curve slope trends.

    The FDC shows the relationship between flow magnitude and
    exceedance probability. Slopes are calculated for different
    segments of the curve to characterize flow variability.

    Parameters
    ----------
    streamflow_data : pd.DataFrame
        Daily streamflow data with columns:
            - water_year: int, water year
            - Q: float, daily discharge in mm/day

    Returns
    -------
    dict
        Dictionary of signature statistics with keys for each metric:
            - FDCall: Slope of entire FDC (log-flow vs exceedance)
            - FDC90th: Slope for low flows (exceedance >= 0.9)
            - FDCmid: Slope for mid-range flows (0.2 <= exceedance <= 0.8)

    Notes
    -----
    - FDC is calculated as log10(Q) vs exceedance probability
    - Exceedance probability = rank / (n + 1)
    - More negative slopes indicate more variable flow regimes
    """
    # Validate required columns
    required_cols = ["water_year", "Q"]
    missing = [c for c in required_cols if c not in streamflow_data.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = streamflow_data

    # Collect results as list of dicts (P13: avoid pre-allocation + boolean .loc)
    results_list = []

    # Calculate FDC slopes for each year
    for yr, year_data in df.groupby("water_year", sort=False):

        row = {"water_year": yr, "FDCall": np.nan, "FDC90th": np.nan, "FDCmid": np.nan}

        # Remove NA and negative values from Q
        q_values = year_data["Q"].dropna().values
        q_values = q_values[q_values >= 0]

        # Sort flows in descending order
        sorted_flows = np.sort(q_values)[::-1]
        n = len(sorted_flows)

        # Calculate exceedance probabilities
        exceedance = np.arange(1, n + 1) / (n + 1)

        # Log-transform flow (add small constant to handle zeros)
        log_flow = np.log10(sorted_flows + 1e-10)

        # Calculate overall slope
        if n >= 10:
            try:
                result = scipy_stats.linregress(exceedance, log_flow)
                row["FDCall"] = result.slope
            except Exception:
                pass

            # Slope for low flows (90th percentile and above)
            low_flow_mask = exceedance >= 0.9
            if low_flow_mask.sum() >= 3:
                try:
                    result = scipy_stats.linregress(
                        exceedance[low_flow_mask],
                        log_flow[low_flow_mask]
                    )
                    row["FDC90th"] = result.slope
                except Exception:
                    pass

            # Slope for mid-range flows (20th to 80th percentile)
            mid_flow_mask = (exceedance >= 0.2) & (exceedance <= 0.8)
            if mid_flow_mask.sum() >= 3:
                try:
                    result = scipy_stats.linregress(
                        exceedance[mid_flow_mask],
                        log_flow[mid_flow_mask]
                    )
                    row["FDCmid"] = result.slope
                except Exception:
                    pass

        results_list.append(row)

    fdc_by_year = pd.DataFrame(results_list)

    # Generate statistics
    result = generate_stats(
        fdc_by_year,
        value_cols=["FDCall", "FDC90th", "FDCmid"],
        year_col="water_year",
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector,
    )

    # Metric columns use the CANONICAL names directly (matching Julia, which
    # never renames): the annual-values collector captures the column name it
    # is given, so a placeholder + post-hoc rename would put the wrong
    # signature name in the annual parquet (found 2026-08-25, Phase 3).
    return result

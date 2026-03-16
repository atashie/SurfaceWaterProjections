"""
Flow Duration Curve signatures.

Calculates slopes of the flow duration curve at different
exceedance probability ranges.
"""

from typing import Dict
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from .stats import generate_stats
from .config import FDC_MIN_DAYS


def analyze_fdc_trends(
    streamflow_data: pd.DataFrame,
    min_days: int = FDC_MIN_DAYS,
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
    min_days : int, default 30
        Minimum number of non-NA days required per water year

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

        row = {"water_year": yr, "slp_all": np.nan, "slp_90th": np.nan, "slp_mid": np.nan}

        # Check minimum data requirement
        if len(year_data) < min_days:
            results_list.append(row)
            continue

        # Remove NA and negative values from Q
        q_values = year_data["Q"].dropna().values
        q_values = q_values[q_values >= 0]

        if len(q_values) < min_days:
            results_list.append(row)
            continue

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
                row["slp_all"] = result.slope
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
                    row["slp_90th"] = result.slope
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
                    row["slp_mid"] = result.slope
                except Exception:
                    pass

        results_list.append(row)

    fdc_by_year = pd.DataFrame(results_list)

    # Generate statistics
    result = generate_stats(
        fdc_by_year,
        value_cols=["slp_all", "slp_90th", "slp_mid"],
        year_col="water_year"
    )

    # Rename columns to match expected output
    renamed = {}
    for key, value in result.items():
        new_key = key.replace("slp_all", "FDCall")
        new_key = new_key.replace("slp_90th", "FDC90th")
        new_key = new_key.replace("slp_mid", "FDCmid")
        renamed[new_key] = value

    return renamed

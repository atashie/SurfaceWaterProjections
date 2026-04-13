"""
High-level convenience function for calculating all streamflow signatures at once.
"""

from typing import Optional

import pandas as pd

from .flow_volumes import calculate_flow_vols_by_year
from .flashiness import analyze_flashiness_trends
from .timing import analyze_flow_timing_trends
from .fdc import analyze_fdc_trends
from .baseflow import analyze_baseflow_indices
from .recession import analyze_recession_parameters
from .pulses import calculate_pulse_metrics, calculate_negative_days
from .runoff_ratios import analyze_Q_PPT_relationships
from .elasticity import calculate_streamflow_elasticity
from .qp_seasonality import calculate_qp_seasonality
from .storage import calculate_average_storage


def calculate_all_signatures(
    gage_data: pd.DataFrame,
    has_climate: bool = False,
    seasonal_flags: list = None,
    climate_data: Optional[pd.DataFrame] = None,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
) -> dict:
    """Calculate all signatures for a single gage.

    Calls each signature function in sequence, catching exceptions per-signature
    so that a failure in one does not prevent others from being calculated.

    Parameters
    ----------
    gage_data : pd.DataFrame
        DataFrame for a single gage with columns: gage_id, date, Q, water_year,
        month, dowy. If has_climate is True, must also have PPT column.
    has_climate : bool, default False
        Whether climate data (PPT column) is available. If True but PPT column
        is missing, climate signatures are silently skipped.
    seasonal_flags : list of dict, optional
        Per-year seasonal completeness flags from preprocess_daily_data().
        Passed through to calculate_flow_vols_by_year() and
        analyze_Q_PPT_relationships() to NaN-out incomplete seasons.
    climate_data : pd.DataFrame, optional
        Separate climate DataFrame to use for climate signatures. If provided,
        this is used instead of gage_data for climate-dependent signatures.
        Must contain PPT column and be joinable with gage_data on water_year.
    trend_completeness : float, optional
        Minimum fraction of non-NaN values required in the time series for
        trend statistics to be computed. Forwarded to generate_stats().
    decade_completeness : float, optional
        Minimum fraction of decades that must contain data for trend statistics
        to be computed. Forwarded to generate_stats().

    Returns
    -------
    dict
        Dictionary mapping signature column names to values (floats or NaN).
    """
    results = {}

    # Common kwargs for trend completeness
    trend_kwargs = dict(
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness,
    )

    # Non-climate signatures
    try:
        results.update(calculate_flow_vols_by_year(
            gage_data, seasonal_flags=seasonal_flags, **trend_kwargs))
    except Exception:
        pass

    try:
        results.update(analyze_flashiness_trends(gage_data, **trend_kwargs))
    except Exception:
        pass

    try:
        results.update(analyze_flow_timing_trends(gage_data, **trend_kwargs))
    except Exception:
        pass

    try:
        results.update(analyze_fdc_trends(gage_data, **trend_kwargs))
    except Exception:
        pass

    try:
        results.update(analyze_baseflow_indices(gage_data, **trend_kwargs))
    except Exception:
        pass

    try:
        results.update(analyze_recession_parameters(gage_data))
    except Exception:
        pass

    try:
        results.update(calculate_pulse_metrics(gage_data, **trend_kwargs))
    except Exception:
        pass

    try:
        results.update(calculate_negative_days(gage_data, **trend_kwargs))
    except Exception:
        pass

    # Climate-dependent signatures
    # Use climate_data if provided, otherwise fall back to gage_data
    climate_df = climate_data if climate_data is not None else gage_data
    if has_climate and "PPT" in climate_df.columns:
        try:
            results.update(analyze_Q_PPT_relationships(
                climate_df, seasonal_flags=seasonal_flags, **trend_kwargs
            ))
        except Exception:
            pass

        try:
            results.update(calculate_streamflow_elasticity(climate_df))
        except Exception:
            pass

        try:
            results.update(calculate_qp_seasonality(climate_df, **trend_kwargs))
        except Exception:
            pass

        try:
            results.update(calculate_average_storage(climate_df, **trend_kwargs))
        except Exception:
            pass

    return results

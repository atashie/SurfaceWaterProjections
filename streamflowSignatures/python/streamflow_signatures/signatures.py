"""
High-level convenience function for calculating all streamflow signatures at once.
"""

import pandas as pd

from .flow_volumes import calculate_flow_vols_by_year
from .flashiness import analyze_flashiness_trends
from .timing import analyze_flow_timing_trends
from .fdc import analyze_fdc_trends
from .baseflow import analyze_baseflow_indices
from .recession import analyze_recession_parameters
from .pulses import calculate_pulse_metrics
from .runoff_ratios import analyze_Q_PPT_relationships
from .elasticity import calculate_streamflow_elasticity
from .qp_seasonality import calculate_qp_seasonality
from .storage import calculate_average_storage


def calculate_all_signatures(gage_data: pd.DataFrame, has_climate: bool = False) -> dict:
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

    Returns
    -------
    dict
        Dictionary mapping signature column names to values (floats or NaN).
    """
    results = {}

    # Non-climate signatures
    try:
        results.update(calculate_flow_vols_by_year(gage_data))
    except Exception:
        pass

    try:
        results.update(analyze_flashiness_trends(gage_data))
    except Exception:
        pass

    try:
        results.update(analyze_flow_timing_trends(gage_data))
    except Exception:
        pass

    try:
        results.update(analyze_fdc_trends(gage_data))
    except Exception:
        pass

    try:
        results.update(analyze_baseflow_indices(gage_data))
    except Exception:
        pass

    try:
        results.update(analyze_recession_parameters(gage_data))
    except Exception:
        pass

    try:
        results.update(calculate_pulse_metrics(gage_data))
    except Exception:
        pass

    # Climate-dependent signatures
    if has_climate and "PPT" in gage_data.columns:
        try:
            results.update(analyze_Q_PPT_relationships(gage_data))
        except Exception:
            pass

        try:
            results.update(calculate_streamflow_elasticity(gage_data))
        except Exception:
            pass

        try:
            results.update(calculate_qp_seasonality(gage_data))
        except Exception:
            pass

        try:
            results.update(calculate_average_storage(gage_data))
        except Exception:
            pass

    return results

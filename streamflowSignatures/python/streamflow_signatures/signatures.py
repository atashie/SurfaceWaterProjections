"""
High-level convenience function for calculating all streamflow signatures at once.
"""

from typing import Optional

import pandas as pd

from .flow_volumes import calculate_flow_vols_by_year
from .flashiness import analyze_flashiness_trends
from .timing import analyze_flow_timing_trends
from .fdc import analyze_fdc_trends
from .baseflow import analyze_baseflow_indices, analyze_baseflow_indices_with_parameters
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
    include_qa_flags: bool = False,
    area_normalized: bool = True,
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
    include_qa_flags : bool, default False
        If True, append 12 QA/QC flag columns from compute_qa_flags() to output.
    area_normalized : bool, default True
        Whether Q is area-normalized to mm/day. When False (gages with no
        drainage area -- Q left in raw m3/s), all Q-to-PPT signatures (runoff
        ratios, elasticity, Q-P seasonality, storage) are skipped because Q and
        PPT units don't match. Q-only signatures are unaffected.

    Returns
    -------
    dict
        Dictionary mapping signature column names to values (floats or NaN),
        plus 12 boolean flag columns if include_qa_flags is True.
    """
    results = {}

    # Season exclusion year counts (per-gage scalar diagnostics)
    if seasonal_flags is not None and len(seasonal_flags) > 0:
        import pandas as _pd
        _sf = _pd.DataFrame(seasonal_flags)
        for season, col in [("winter", "winter_complete"), ("spring", "spring_complete"),
                            ("summer", "summer_complete"), ("fall", "fall_complete")]:
            if col in _sf.columns:
                results[f"season_excluded_years_{season}"] = float(
                    (~_sf[col].astype(bool)).sum()
                )
            else:
                results[f"season_excluded_years_{season}"] = 0.0

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

    # Recession: inherently sparse (event-based), no trend completeness
    recession_alpha = float('nan')
    try:
        recession_results = analyze_recession_parameters(gage_data)
        # Extract scalar for parameterized BFI before merging
        recession_alpha = recession_results.get(
            "recession_alpha_point_cloud_linear_reservoir", float('nan'))
        results.update(recession_results)
    except Exception:
        pass

    # Parameterized BFI using recession-derived alpha (requires recession to have run first)
    # Uses trend completeness (same as fixed-parameter BFI)
    try:
        results.update(analyze_baseflow_indices_with_parameters(
            gage_data, recession_alpha, **trend_kwargs))
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
    # Use climate_data if provided, otherwise fall back to gage_data.
    # Gate: Q-to-PPT signatures are undefined when Q is not area-normalized
    # (raw m3/s vs PPT in mm -- units don't match, so Q/P ratios, dQ/dP, and
    # cumsum(P - Q) are all meaningless).
    if has_climate and not area_normalized:
        import logging as _logging
        _logging.info(
            "Skipping Q-to-PPT signatures (area_normalized=False -- "
            "Q in raw m3/s, PPT in mm)")
    climate_df = climate_data if climate_data is not None else gage_data
    if has_climate and area_normalized and "PPT" in climate_df.columns:
        try:
            results.update(analyze_Q_PPT_relationships(
                climate_df, seasonal_flags=seasonal_flags, **trend_kwargs
            ))
        except Exception:
            pass

        # Elasticity: rolling window produces fewer values than years, no trend completeness
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

    # QA/QC flags (optional)
    if include_qa_flags:
        import logging as _logging
        try:
            from .qa_qc import compute_qa_flags, get_flag_columns
            row_df = pd.DataFrame([results])
            flagged_df = compute_qa_flags(row_df)
            for col in get_flag_columns():
                if col in flagged_df.columns:
                    val = flagged_df[col].iloc[0]
                    results[col] = bool(val) if not pd.isna(val) else False
        except Exception as _e:
            _logging.warning("QA/QC flags failed: %s", _e)

    return results

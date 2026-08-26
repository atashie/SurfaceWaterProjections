"""
High-level convenience function for calculating all streamflow signatures at once.
"""

import logging
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
from .snow import calculate_snow_metrics
from .drought import calculate_drought_metrics
from .config import DROUGHT_ENABLED

logger = logging.getLogger(__name__)


def calculate_all_signatures(
    gage_data: pd.DataFrame,
    has_climate: bool = False,
    seasonal_flags: list = None,
    climate_data: Optional[pd.DataFrame] = None,
    trend_completeness: Optional[float] = None,
    decade_completeness: Optional[float] = None,
    include_qa_flags: bool = False,
    area_normalized: bool = True,
    changepoint: Optional[dict] = None,
    min_values_for_stats: Optional[int] = None,
    collector: Optional[object] = None,
    snow_data: Optional[pd.DataFrame] = None,
    snow_climate_years: Optional[list] = None,
    gage_id: Optional[str] = None,
) -> dict:
    """Calculate all signatures for a single gage.

    Calls each signature function in sequence, catching exceptions per-signature
    family so that a failure in one does not prevent others from being
    calculated. Every caught failure is LOGGED with the gage id and family name
    (mirroring the Julia orchestrator's @warn) — a silent failure would surface
    only as mysteriously missing columns for this gage.

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
    changepoint : dict, optional
        Changepoint (Pettitt) configuration forwarded to every signature
        family's generate_stats() call — including recession and elasticity
        (their trend-gate/stats-floor exemptions do NOT extend to changepoint,
        matching the Julia orchestrator). Keys: start_year, end_year,
        min_total_obs, min_segment_obs.
    collector : AnnualCollector, optional
        When supplied, every family's annual series is appended to it in long
        format (before gating). Threaded to ALL families including the
        floor-exempt recession and elasticity, matching the Julia orchestrator.
    min_values_for_stats : int, optional
        Stats floor: metrics with fewer non-NaN annual values emit NaN for all
        8 statistics + changepoint fields. NOT passed to recession/elasticity
        (inherently sparse; same exemption as the trend gates).
    snow_data : pd.DataFrame, optional
        EXPLICIT frame filtered to the gage's SWE-valid years. Snow metrics run
        only when this is supplied and carries an SWE column — an SWE column in
        `gage_data` is NEVER used implicitly (matching Julia; prevents
        SWE-invalid years leaking in). May legitimately have 0 rows.
    snow_climate_years : list of int, optional
        valid_climate_years for the PPT-dependent `swe_max_to_ppt` metric.
    gage_id : str, optional
        Gage identifier used only for warning context on per-family failures.

    Returns
    -------
    dict
        Dictionary mapping signature column names to values (floats or NaN),
        plus 12 boolean flag columns if include_qa_flags is True.
    """
    results = {}
    gid = gage_id if gage_id is not None else "unknown"

    def run_family(family, func, *args, **kwargs):
        """Run one signature family; log (never swallow silently) failures."""
        try:
            out = func(*args, **kwargs)
            results.update(out)
            return out
        except Exception as e:
            logger.warning("Gage %s: %s signatures failed: %s: %s",
                           gid, family, type(e).__name__, e)
            return None

    # Season exclusion year counts (per-gage scalar diagnostics)
    if seasonal_flags is not None and len(seasonal_flags) > 0:
        _sf = pd.DataFrame(seasonal_flags)
        for season, col in [("winter", "winter_complete"), ("spring", "spring_complete"),
                            ("summer", "summer_complete"), ("fall", "fall_complete")]:
            if col in _sf.columns:
                results[f"season_excluded_years_{season}"] = float(
                    (~_sf[col].astype(bool)).sum()
                )
            else:
                results[f"season_excluded_years_{season}"] = 0.0

    # Common kwargs: trend gates + stats floor + changepoint. NOTE the exemption
    # asymmetry (matching Julia's orchestrator): recession and elasticity are
    # exempt from the trend gates AND the stats floor (inherently sparse), but
    # changepoint IS passed to them explicitly below.
    trend_kwargs = dict(
        trend_completeness=trend_completeness,
        decade_completeness=decade_completeness,
        min_values_for_stats=min_values_for_stats,
        changepoint=changepoint,
        collector=collector,
    )

    # Non-climate signatures
    run_family("flow volumes", calculate_flow_vols_by_year,
               gage_data, seasonal_flags=seasonal_flags, **trend_kwargs)
    run_family("flashiness", analyze_flashiness_trends, gage_data, **trend_kwargs)
    run_family("flow timing", analyze_flow_timing_trends, gage_data, **trend_kwargs)
    run_family("FDC", analyze_fdc_trends, gage_data, **trend_kwargs)
    run_family("baseflow", analyze_baseflow_indices, gage_data, **trend_kwargs)

    # Recession: inherently sparse (event-based), no trend completeness —
    # changepoint still applies
    recession_alpha = float("nan")
    recession_results = run_family("recession", analyze_recession_parameters,
                                   gage_data, changepoint=changepoint, collector=collector)
    if recession_results is not None:
        recession_alpha = recession_results.get(
            "recession_alpha_point_cloud_linear_reservoir", float("nan"))

    # Parameterized BFI using recession-derived alpha (requires recession first).
    # Uses trend completeness (same as fixed-parameter BFI)
    run_family("parameterized baseflow", analyze_baseflow_indices_with_parameters,
               gage_data, recession_alpha, **trend_kwargs)

    run_family("pulse metrics", calculate_pulse_metrics, gage_data, **trend_kwargs)
    run_family("negative days", calculate_negative_days, gage_data, **trend_kwargs)

    # Streamflow drought (config-gated: absent `drought` section => family off)
    if DROUGHT_ENABLED:
        run_family("drought", calculate_drought_metrics, gage_data, **trend_kwargs)

    # Climate-dependent signatures
    # Use climate_data if provided, otherwise fall back to gage_data.
    # Gate: Q-to-PPT signatures are undefined when Q is not area-normalized
    # (raw m3/s vs PPT in mm -- units don't match, so Q/P ratios, dQ/dP, and
    # cumsum(P - Q) are all meaningless).
    if has_climate and not area_normalized:
        logger.info(
            "Gage %s: skipping Q-to-PPT signatures (area_normalized=False -- "
            "Q in raw m3/s, PPT in mm)", gid)
    climate_df = climate_data if climate_data is not None else gage_data
    if has_climate and area_normalized and "PPT" in climate_df.columns:
        run_family("runoff ratios", analyze_Q_PPT_relationships,
                   climate_df, seasonal_flags=seasonal_flags, **trend_kwargs)
        # Elasticity: rolling window produces fewer values than years, no trend
        # completeness — changepoint still applies
        run_family("elasticity", calculate_streamflow_elasticity,
                   climate_df, changepoint=changepoint, collector=collector)
        run_family("Q-P seasonality", calculate_qp_seasonality, climate_df, **trend_kwargs)
        run_family("storage", calculate_average_storage, climate_df, **trend_kwargs)

    # Snow metrics: explicit opt-in frame only (never derived from gage_data)
    if snow_data is not None and "SWE" in snow_data.columns:
        run_family("snow", calculate_snow_metrics, snow_data,
                   valid_climate_years=snow_climate_years, **trend_kwargs)

    # QA/QC flags (optional)
    if include_qa_flags:
        try:
            from .qa_qc import compute_qa_flags, get_flag_columns
            row_df = pd.DataFrame([results])
            flagged_df = compute_qa_flags(row_df)
            for col in get_flag_columns():
                if col in flagged_df.columns:
                    val = flagged_df[col].iloc[0]
                    results[col] = bool(val) if not pd.isna(val) else False
        except Exception as _e:
            logger.warning("Gage %s: QA/QC flags failed: %s", gid, _e)

    return results

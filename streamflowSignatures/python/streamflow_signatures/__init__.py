"""
Streamflow Signatures - Hydrological signature extraction from streamflow time series.

This package provides functions to calculate standardized hydrological signatures
from daily streamflow (and optionally climate) data. All signatures produce
8 statistics via generate_stats():

    - _senn_slp: Theil-Sen slope
    - _linear_slp: Linear regression slope
    - _spearman_rho: Spearman correlation coefficient
    - _spearman_pval: Spearman p-value
    - _mk_rho: Mann-Kendall tau
    - _mk_pval: Mann-Kendall p-value
    - _mean: Arithmetic mean
    - _median: Median
"""

__version__ = "0.1.0"

from .stats import generate_stats, theil_sen_slope, mann_kendall_test, empty_stats, STAT_SUFFIXES, CP_SUFFIXES, AnnualCollector
from .changepoint import pettitt_test, segment_differential_metrics
from .io import read_parquet, write_signatures, validate_schema, add_water_year_columns, filter_qualifying_years, preprocess_daily_data
from .signatures import calculate_all_signatures
from .flow_volumes import calculate_flow_vols_by_year
from .flashiness import analyze_flashiness_trends
from .timing import analyze_flow_timing_trends
from .fdc import analyze_fdc_trends
from .baseflow import analyze_baseflow_indices, analyze_baseflow_indices_with_parameters, eckhardt_filter, lyne_hollick_filter
from .recession import analyze_recession_parameters
from .pulses import calculate_pulse_metrics, calculate_negative_days
from .runoff_ratios import analyze_Q_PPT_relationships
from .elasticity import calculate_streamflow_elasticity
from .qp_seasonality import calculate_qp_seasonality
from .storage import calculate_average_storage
from .snow import calculate_snow_metrics, SNOW_METRICS
from .drought import calculate_drought_metrics, DROUGHT_METRICS, DROUGHT_THRESHOLD_SCALARS

__all__ = [
    # Stats
    "generate_stats",
    "theil_sen_slope",
    "mann_kendall_test",
    "empty_stats",
    "AnnualCollector",
    "STAT_SUFFIXES",
    "CP_SUFFIXES",
    # Changepoint (Pettitt)
    "pettitt_test",
    "segment_differential_metrics",
    # I/O
    "read_parquet",
    "write_signatures",
    "validate_schema",
    "add_water_year_columns",
    "filter_qualifying_years",
    "preprocess_daily_data",
    # High-level
    "calculate_all_signatures",
    # Signatures - Simple (no climate data required)
    "calculate_flow_vols_by_year",
    "analyze_flashiness_trends",
    "analyze_flow_timing_trends",
    "analyze_fdc_trends",
    # Signatures - Complex (no climate data required)
    "analyze_baseflow_indices",
    "analyze_baseflow_indices_with_parameters",
    "eckhardt_filter",
    "lyne_hollick_filter",
    "analyze_recession_parameters",
    "calculate_pulse_metrics",
    "calculate_negative_days",
    # Signatures - Climate-dependent (require PPT column)
    "analyze_Q_PPT_relationships",
    "calculate_streamflow_elasticity",
    "calculate_qp_seasonality",
    "calculate_average_storage",
    # Snow (Daymet SWE)
    "calculate_snow_metrics",
    "SNOW_METRICS",
    # Streamflow drought
    "calculate_drought_metrics",
    "DROUGHT_METRICS",
    "DROUGHT_THRESHOLD_SCALARS",
]

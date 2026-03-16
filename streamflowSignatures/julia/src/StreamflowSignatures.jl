"""
StreamflowSignatures - Hydrological signature extraction from streamflow time series.

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
module StreamflowSignatures

using DataFrames
using Dates
using Statistics
using StatsBase
using LinearAlgebra
using HypothesisTests
using JSON3

# Include configuration first
include("config.jl")

# Include submodules
include("stats.jl")
include("io.jl")
include("flow_volumes.jl")
include("flashiness.jl")
include("timing.jl")
include("fdc.jl")
include("baseflow.jl")
include("recession.jl")
include("pulses.jl")
include("runoff_ratios.jl")
include("elasticity.jl")
include("qp_seasonality.jl")
include("storage.jl")
include("qa_qc.jl")
include("metadata.jl")
include("signatures.jl")

# Export stats functions
export generate_stats, theil_sen_slope, mann_kendall_test

# Export I/O functions
export read_parquet, write_signatures, add_water_year_columns, filter_qualifying_years

# Export high-level functions
export calculate_all_signatures

# Export simple signatures
export calculate_flow_vols_by_year
export analyze_flashiness_trends
export analyze_flow_timing_trends
export analyze_fdc_trends

# Export complex signatures
export analyze_baseflow_indices, eckhardt_filter, lyne_hollick_filter
export analyze_recession_parameters
export calculate_pulse_metrics

# Export climate-dependent signatures
export analyze_Q_PPT_relationships
export calculate_streamflow_elasticity
export calculate_qp_seasonality
export calculate_average_storage

# Export QA/QC functions
export compute_qa_flags, get_flag_columns

# Export metadata functions
export load_gages_ii_interference, enrich_signatures_with_metadata

# Export config constants (used by benchmark runners)
export CFG_MIN_NUM_YEARS, CFG_MIN_FRAC_GOOD_DATA, CFG_MIN_Q_VALUE,
       CFG_MIN_DAYS_ABOVE_THRESHOLD, CFG_MIN_NONA_DAYS_ANNUAL

end # module

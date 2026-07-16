"""
Centralized configuration for streamflow signature calculations.

Loads configuration from the shared JSON config file and exports
parameters as constants for use across all signature modules.
"""

using JSON3

# Find config file relative to this module
const _CONFIG_DIR = joinpath(@__DIR__, "..", "..", "config")
const _DEFAULT_CONFIG_FILE = joinpath(_CONFIG_DIR, "signatures_config.json")
const _CONFIG_FILE = get(ENV, "STREAMFLOW_CONFIG", _DEFAULT_CONFIG_FILE)

function load_config()
    if !isfile(_CONFIG_FILE)
        error("Config file not found: $_CONFIG_FILE")
    end
    return JSON3.read(read(_CONFIG_FILE, String))
end

# Load config at module load time
const _config = load_config()

# =============================================================================
# Filtering Parameters
# =============================================================================
const CFG_MIN_NUM_YEARS = _config.filtering.min_num_years
const CFG_MIN_FRAC_GOOD_DATA = _config.filtering.min_frac_good_data
const CFG_MIN_Q_VALUE = _config.filtering.min_q_value
const CFG_MIN_DAYS_ABOVE_THRESHOLD = _config.filtering.min_days_above_threshold
const CFG_MIN_NONA_DAYS_ANNUAL = _config.filtering.min_nona_days_annual

# Experiment controls — JSON defaults with ENV overrides
const CFG_START_WATER_YEAR = let
    json_val = get(get(_config, "filtering", Dict()), "start_water_year", nothing)
    env_val = get(ENV, "STREAMFLOW_START_WATER_YEAR", "")
    v = if env_val != ""
        try parse(Int, env_val) catch e
            error("Invalid STREAMFLOW_START_WATER_YEAR='$env_val': must be an integer")
        end
    else
        json_val
    end
    v === nothing ? nothing : Int(v)
end

const CFG_MIN_QUALIFYING_DATA_FRACTION = let
    json_val = get(get(_config, "filtering", Dict()), "min_qualifying_data_fraction", nothing)
    env_val = get(ENV, "STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION", "")
    v = if env_val != ""
        try parse(Float64, env_val) catch e
            error("Invalid STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION='$env_val': must be a number")
        end
    else
        json_val
    end
    v === nothing ? nothing : Float64(v)
end

# =============================================================================
# Water Year
# =============================================================================
const CFG_WATER_YEAR_START_MONTH = _config.water_year.start_month

# =============================================================================
# Baseflow Parameters
# =============================================================================
const CFG_ECKHARDT_BFIMAX = _config.baseflow.eckhardt_bfimax
const CFG_ECKHARDT_ALPHA = _config.baseflow.eckhardt_alpha
const CFG_LYNE_HOLLICK_ALPHA = _config.baseflow.lyne_hollick_alpha
const CFG_LYNE_HOLLICK_PASSES = _config.baseflow.lyne_hollick_passes

# =============================================================================
# Recession Parameters
# =============================================================================
const CFG_RECESSION_MIN_LENGTH = _config.recession.min_length
const CFG_RECESSION_MIN_EVENTS = _config.recession.min_events

# =============================================================================
# Pulse Parameters
# =============================================================================
const CFG_HIGH_PULSE_PERCENTILE = _config.pulses.high_percentile
const CFG_LOW_PULSE_PERCENTILE = _config.pulses.low_percentile
const CFG_FLOW_REVERSAL_THRESHOLD = _config.pulses.flow_reversal_threshold

# =============================================================================
# Timing Parameters
# =============================================================================
const CFG_D_PERCENTILES = collect(_config.timing.d_percentiles)

# =============================================================================
# Elasticity Parameters
# =============================================================================
const CFG_ELASTICITY_WINDOW_YEARS = _config.elasticity.window_years
const CFG_ELASTICITY_MIN_YEARS = _config.elasticity.min_years
const CFG_ELASTICITY_MIN_ANNUAL_PPT = _config.elasticity.min_annual_ppt

# =============================================================================
# Q-P Seasonality Parameters
# =============================================================================
const CFG_QP_SLOPE_WINDOW_DAYS = _config.qp_seasonality.slope_window_days
const CFG_QP_MIN_YEARS = _config.qp_seasonality.min_years

# =============================================================================
# Runoff Ratio Parameters
# =============================================================================
const CFG_RUNOFF_MIN_ANNUAL_PPT = _config.runoff_ratios.min_annual_ppt
const CFG_RUNOFF_MIN_SEASONAL_PPT = _config.runoff_ratios.min_seasonal_ppt

# =============================================================================
# Flow Volumes Parameters
# =============================================================================
const CFG_FLOW_PERCENTILES = collect(_config.flow_volumes.percentiles)

# =============================================================================
# Storage Parameters
# =============================================================================
const CFG_STORAGE_MIN_YEARS = _config.storage.min_years

# =============================================================================
# Snow Parameters (absent-section defaults match the shipped config)
# =============================================================================
const _snow_config = get(_config, "snow", Dict())
const CFG_SNOW_SWE_THRESHOLD_MM = Float64(get(_snow_config, "swe_day_threshold_mm", 10.0))
const CFG_SNOW_SEASONAL_MIN_DAYS = Int(get(_snow_config, "seasonal_spell_min_days", 60))
const CFG_SNOW_MELT_COM_FRACTION = Float64(get(_snow_config, "melt_com_fraction", 0.5))
const CFG_SNOW_MIN_ANNUAL_PPT_MM = Float64(get(_snow_config, "min_annual_ppt_mm", 10))

# =============================================================================
# QA/QC Parameters
# =============================================================================
const CFG_QAQC_QANN_RANGE = Tuple(_config.qa_qc.qann_range)
const CFG_QAQC_BFI_RANGE = Tuple(_config.qa_qc.bfi_range)
const CFG_QAQC_FLASHINESS_RANGE = Tuple(_config.qa_qc.flashiness_range)
const CFG_QAQC_TQMEAN_RANGE = Tuple(_config.qa_qc.tqmean_range)
const CFG_QAQC_D50_RANGE = Tuple(_config.qa_qc.d50_range)
const CFG_QAQC_ELASTICITY_RANGE = Tuple(_config.qa_qc.elasticity_range)
const CFG_QAQC_RUNOFF_RATIO_RANGE = Tuple(_config.qa_qc.runoff_ratio_range)
const CFG_QAQC_SEASONAL_SUM_TOLERANCE = _config.qa_qc.seasonal_sum_tolerance
const CFG_QAQC_MAX_NA_FRACTION = _config.qa_qc.max_na_fraction

# =============================================================================
# Metadata Parameters
# =============================================================================
const CFG_INCLUDE_HUMAN_INTERFERENCE = _config.metadata.include_human_interference
const CFG_GAGES_II_DIR = _config.metadata.gages_ii_dir
const CFG_HYDAT_PATH = _config.metadata.hydat_path
const CFG_INTERFERENCE_COLUMNS = collect(_config.metadata.interference_columns)

# ==============================================================================
# NA HANDLING
# ==============================================================================
const _na_handling = get(_config, "na_handling", Dict())
const _na_interp = get(_na_handling, "interpolation", Dict())
const _na_reject = get(_na_handling, "year_rejection", Dict())
const _na_const_sd = get(_na_handling, "constant_sd_flag", Dict())
const _na_trend = get(_na_handling, "trend_completeness", Dict())
const _na_seasonal = get(_na_handling, "seasonal_completeness", Dict())
const _na_climate = get(_na_handling, "climate_na_policy", Dict())

const CFG_NA_MAX_GAP_DAYS = Int(get(_na_interp, "max_gap_days", 3))
const CFG_NA_INTERPOLATION_METHOD = get(_na_interp, "method", "linear")
const CFG_NA_INTERNAL_ONLY = get(_na_interp, "internal_only", true)
const CFG_NA_MAX_RAW_NA_PER_YEAR = Int(get(_na_reject, "max_raw_na_per_year", 30))
const CFG_NA_REJECT_NEGATIVE_FLOW = get(_na_reject, "reject_negative_flow", false)
const CFG_NA_REJECT_RESIDUAL_NA = get(_na_reject, "reject_residual_na", true)
const CFG_NA_CONSTANT_SD_ENABLED = get(_na_const_sd, "enabled", true)
const CFG_NA_CONSTANT_SD_MIN_DAYS = Int(get(_na_const_sd, "min_nonzero_days_per_month", 15))
const CFG_NA_CONSTANT_SD_MAX_UNIQUE = Int(get(_na_const_sd, "max_unique_values", 1))
const CFG_NA_TREND_MIN_FRACTION = Float64(get(_na_trend, "min_fraction", 0.80))
const CFG_NA_DECADE_MIN_FRACTION = Float64(get(_na_trend, "decade_min_fraction", 0.80))
const CFG_NA_SEASONAL_MIN_FRACTION = Float64(get(_na_seasonal, "min_fraction", 0.80))
const CFG_NA_SEASONAL_USE_RAW = get(_na_seasonal, "use_raw_observations", true)
const _na_season_defs = get(_na_seasonal, "season_definitions", Dict())
const CFG_NA_SEASON_MONTHS = Dict{String, Vector{Int}}(
    "win" => Int.(collect(get(_na_season_defs, "winter", [12, 1, 2]))),
    "spr" => Int.(collect(get(_na_season_defs, "spring", [3, 4, 5]))),
    "sum" => Int.(collect(get(_na_season_defs, "summer", [6, 7, 8]))),
    "fal" => Int.(collect(get(_na_season_defs, "fall", [9, 10, 11])))
)
const CFG_NA_MAX_RAW_NA_PPT = Int(get(_na_climate, "max_raw_na_per_year_ppt", 30))
const CFG_NA_MAX_GAP_PPT = Int(get(_na_climate, "max_interpolation_gap_ppt", 3))
const CFG_NA_REJECT_NEGATIVE_PPT = get(_na_climate, "reject_negative_ppt", true)
const _na_snow = get(_na_handling, "snow_na_policy", Dict())
const CFG_NA_MAX_RAW_NA_SWE = Int(get(_na_snow, "max_raw_na_per_year_swe", 30))
const CFG_NA_MAX_GAP_SWE = Int(get(_na_snow, "max_interpolation_gap_swe", 3))
const CFG_NA_REJECT_NEGATIVE_SWE = get(_na_snow, "reject_negative_swe", true)

# Changepoint analysis configuration
const _cp_config = get(_config, "changepoint", Dict())
const CFG_CHANGEPOINT_ENABLED = Bool(get(_cp_config, "enabled", false))
const CFG_CP_START_WATER_YEAR = Int(get(_cp_config, "start_water_year", 1980))
const CFG_CP_END_WATER_YEAR = Int(get(_cp_config, "end_water_year", 2024))
const CFG_CP_MIN_TOTAL_OBS = Int(get(_cp_config, "min_total_obs", 20))
const CFG_CP_MIN_SEGMENT_OBS = Int(get(_cp_config, "min_segment_obs", 10))

# Legacy filtering flag
const CFG_USE_LEGACY_FILTERING = get(get(_config, "filtering", Dict()), "use_legacy_filtering", true)

# Stats floor (July 2026): minimum non-NA annual values a metric needs before ANY
# of its 8 statistics (+ changepoint fields) are computed. `nothing` when the
# section is absent => no floor (legacy behavior). Exemptions (recession,
# elasticity) are enforced at the orchestration layer, like trend_completeness.
const _stats_floor = get(_config, "stats_floor", Dict())
const CFG_MIN_VALUES_FOR_STATS = let
    v = get(_stats_floor, "min_values_for_stats", nothing)
    v === nothing ? nothing : Int(v)
end

# Annual values export (per-year annual signature values, long-format parquet).
# Default FALSE when the section is absent: config is loaded once at module import,
# so a permissive fallback would silently enable a large new artifact for old configs.
const _annual_values = get(_config, "annual_values", Dict())
const CFG_SAVE_ANNUAL_VALUES = Bool(get(_annual_values, "save", false))

"""
    get_config()

Return the full configuration dictionary.
"""
function get_config()
    return _config
end

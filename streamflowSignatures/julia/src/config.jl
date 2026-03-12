"""
Centralized configuration for streamflow signature calculations.

Loads configuration from the shared JSON config file and exports
parameters as constants for use across all signature modules.
"""

using JSON3

# Find config file relative to this module
const _CONFIG_DIR = joinpath(@__DIR__, "..", "..", "config")
const _CONFIG_FILE = joinpath(_CONFIG_DIR, "signatures_config.json")

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
const CFG_BASEFLOW_MIN_DAYS = _config.baseflow.min_days
const CFG_BASEFLOW_MAX_MISSING_FRAC = _config.baseflow.max_missing_frac

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
const CFG_PULSES_MIN_DAYS = _config.pulses.min_days

# =============================================================================
# Timing Parameters
# =============================================================================
const CFG_TIMING_MIN_DAYS = _config.timing.min_days
const CFG_D_PERCENTILES = collect(_config.timing.d_percentiles)

# =============================================================================
# Elasticity Parameters
# =============================================================================
const CFG_ELASTICITY_WINDOW_YEARS = _config.elasticity.window_years
const CFG_ELASTICITY_MIN_YEARS = _config.elasticity.min_years
const CFG_ELASTICITY_MIN_DATA_COMPLETENESS = _config.elasticity.min_data_completeness
const CFG_ELASTICITY_MIN_ANNUAL_PPT = _config.elasticity.min_annual_ppt

# =============================================================================
# Q-P Seasonality Parameters
# =============================================================================
const CFG_QP_SLOPE_WINDOW_DAYS = _config.qp_seasonality.slope_window_days
const CFG_QP_MIN_YEARS = _config.qp_seasonality.min_years
const CFG_QP_MIN_DAYS_PER_YEAR = _config.qp_seasonality.min_days_per_year
const CFG_QP_MAX_NA_FRAC = _config.qp_seasonality.max_na_frac

# =============================================================================
# Runoff Ratio Parameters
# =============================================================================
const CFG_RUNOFF_MIN_ANNUAL_PPT = _config.runoff_ratios.min_annual_ppt
const CFG_RUNOFF_MIN_SEASONAL_PPT = _config.runoff_ratios.min_seasonal_ppt

# =============================================================================
# Flow Volumes Parameters
# =============================================================================
const CFG_FLOW_VOLUMES_MIN_DAYS = _config.flow_volumes.min_days
const CFG_FLOW_PERCENTILES = collect(_config.flow_volumes.percentiles)

# =============================================================================
# FDC Parameters
# =============================================================================
const CFG_FDC_MIN_DAYS = _config.fdc.min_days

# =============================================================================
# Flashiness Parameters
# =============================================================================
const CFG_FLASHINESS_MIN_DAYS = _config.flashiness.min_days
const CFG_FLASHINESS_MAX_MISSING_FRAC = _config.flashiness.max_missing_frac

# =============================================================================
# Storage Parameters
# =============================================================================
const CFG_STORAGE_MIN_YEARS = _config.storage.min_years
const CFG_STORAGE_MIN_DAYS_PER_YEAR = _config.storage.min_days_per_year
const CFG_STORAGE_MAX_NA_FRAC = _config.storage.max_na_frac

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

"""
    get_config()

Return the full configuration dictionary.
"""
function get_config()
    return _config
end

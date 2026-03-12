# Streamflow Signatures - JSON Config Loader
# Loads configuration from the shared signatures_config.json file
# ==============================================================================
# This provides an alternative to hardcoded values in config.R
# Usage: source("R/load_config.R") at the start of processing scripts
# ==============================================================================

library(jsonlite)

# Find config file relative to this script
# Works whether sourced from project root or from R/ directory
find_config_file <- function() {
  # Check multiple potential locations
  potential_paths <- c(
    "config/signatures_config.json",            # From project root
    "../config/signatures_config.json",         # From R/ directory
    file.path(dirname(sys.frame(1)$ofile), "..", "config", "signatures_config.json")  # Relative to this file
  )

  for (path in potential_paths) {
    if (file.exists(path)) {
      return(normalizePath(path))
    }
  }

  # If not found, try using the script path if available
  script_path <- tryCatch({
    dirname(sys.frame(1)$ofile)
  }, error = function(e) getwd())

  config_path <- file.path(script_path, "..", "config", "signatures_config.json")
  if (file.exists(config_path)) {
    return(normalizePath(config_path))
  }

  stop("Config file not found. Expected at config/signatures_config.json")
}

# Load configuration from JSON
load_signatures_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    config_path <- find_config_file()
  }

  if (!file.exists(config_path)) {
    stop(paste("Config file not found:", config_path))
  }

  message("Loading config from: ", config_path)
  config <- fromJSON(config_path)
  return(config)
}

# Load the config
CFG <- load_signatures_config()

# ==============================================================================
# Export config values as R variables (matching existing config.R structure)
# ==============================================================================

# Filtering Parameters
CFG_MIN_NUM_YEARS <- CFG$filtering$min_num_years
CFG_MIN_FRAC_GOOD_DATA <- CFG$filtering$min_frac_good_data
CFG_MIN_Q_VALUE <- CFG$filtering$min_q_value
CFG_MIN_DAYS_ABOVE_THRESHOLD <- CFG$filtering$min_days_above_threshold
CFG_MIN_NONA_DAYS_ANNUAL <- CFG$filtering$min_nona_days_annual

# Combined format matching existing MIN_Q_VALUE_AND_DAYS
CFG_MIN_Q_VALUE_AND_DAYS <- c(CFG_MIN_Q_VALUE, CFG_MIN_DAYS_ABOVE_THRESHOLD)

# Water Year
CFG_WATER_YEAR_START_MONTH <- CFG$water_year$start_month

# Baseflow Parameters
CFG_ECKHARDT_BFIMAX <- CFG$baseflow$eckhardt_bfimax
CFG_ECKHARDT_ALPHA <- CFG$baseflow$eckhardt_alpha
CFG_LYNE_HOLLICK_ALPHA <- CFG$baseflow$lyne_hollick_alpha
CFG_LYNE_HOLLICK_PASSES <- CFG$baseflow$lyne_hollick_passes
CFG_BASEFLOW_MIN_DAYS <- CFG$baseflow$min_days
CFG_BASEFLOW_MAX_MISSING_FRAC <- CFG$baseflow$max_missing_frac

# Recession Parameters
CFG_RECESSION_MIN_LENGTH <- CFG$recession$min_length
CFG_RECESSION_MIN_EVENTS <- CFG$recession$min_events

# Pulse Parameters
CFG_HIGH_PULSE_PERCENTILE <- CFG$pulses$high_percentile
CFG_LOW_PULSE_PERCENTILE <- CFG$pulses$low_percentile
CFG_FLOW_REVERSAL_THRESHOLD <- CFG$pulses$flow_reversal_threshold
CFG_PULSES_MIN_DAYS <- CFG$pulses$min_days

# Timing Parameters
CFG_TIMING_MIN_DAYS <- CFG$timing$min_days
CFG_D_PERCENTILES <- CFG$timing$d_percentiles

# Elasticity Parameters
CFG_ELASTICITY_WINDOW_YEARS <- CFG$elasticity$window_years
CFG_ELASTICITY_MIN_YEARS <- CFG$elasticity$min_years
CFG_ELASTICITY_MIN_DATA_COMPLETENESS <- CFG$elasticity$min_data_completeness
CFG_ELASTICITY_MIN_ANNUAL_PPT <- CFG$elasticity$min_annual_ppt

# Q-P Seasonality Parameters
CFG_QP_SLOPE_WINDOW_DAYS <- CFG$qp_seasonality$slope_window_days
CFG_QP_MIN_YEARS <- CFG$qp_seasonality$min_years
CFG_QP_MIN_DAYS_PER_YEAR <- CFG$qp_seasonality$min_days_per_year
CFG_QP_MAX_NA_FRAC <- CFG$qp_seasonality$max_na_frac

# Runoff Ratio Parameters
CFG_RUNOFF_MIN_ANNUAL_PPT <- CFG$runoff_ratios$min_annual_ppt
CFG_RUNOFF_MIN_SEASONAL_PPT <- CFG$runoff_ratios$min_seasonal_ppt

# Flow Volumes Parameters
CFG_FLOW_VOLUMES_MIN_DAYS <- CFG$flow_volumes$min_days
CFG_FLOW_PERCENTILES <- CFG$flow_volumes$percentiles

# FDC Parameters
CFG_FDC_MIN_DAYS <- CFG$fdc$min_days

# Flashiness Parameters
CFG_FLASHINESS_MIN_DAYS <- CFG$flashiness$min_days
CFG_FLASHINESS_MAX_MISSING_FRAC <- CFG$flashiness$max_missing_frac

# Storage Parameters
CFG_STORAGE_MIN_YEARS <- CFG$storage$min_years
CFG_STORAGE_MIN_DAYS_PER_YEAR <- CFG$storage$min_days_per_year
CFG_STORAGE_MAX_NA_FRAC <- CFG$storage$max_na_frac

# QA/QC Parameters
CFG_QAQC_QANN_RANGE <- CFG$qa_qc$qann_range
CFG_QAQC_BFI_RANGE <- CFG$qa_qc$bfi_range
CFG_QAQC_FLASHINESS_RANGE <- CFG$qa_qc$flashiness_range
CFG_QAQC_TQMEAN_RANGE <- CFG$qa_qc$tqmean_range
CFG_QAQC_D50_RANGE <- CFG$qa_qc$d50_range
CFG_QAQC_ELASTICITY_RANGE <- CFG$qa_qc$elasticity_range
CFG_QAQC_RUNOFF_RATIO_RANGE <- CFG$qa_qc$runoff_ratio_range
CFG_QAQC_SEASONAL_SUM_TOLERANCE <- CFG$qa_qc$seasonal_sum_tolerance
CFG_QAQC_MAX_NA_FRACTION <- CFG$qa_qc$max_na_fraction

# ==============================================================================
# Utility function to get the full config
# ==============================================================================

get_signatures_config <- function() {
  return(CFG)
}

message("Config loaded successfully. ", length(names(CFG)), " sections available.")

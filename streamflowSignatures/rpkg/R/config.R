# Package-internal environment for configuration values
pkg_env <- new.env(parent = emptyenv())

# %||% operator for NULL-coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

#' @importFrom jsonlite fromJSON
.onLoad <- function(libname, pkgname) {
  # Allow override via environment variable
  config_path <- Sys.getenv("STREAMFLOW_SIGNATURES_CONFIG", unset = "")
  if (config_path == "") {
    config_path <- system.file("config", "signatures_config.json",
                               package = pkgname)
  }

  if (!file.exists(config_path)) {
    stop("streamflowsignatures: config file not found at ", config_path)
  }

  cfg <- jsonlite::fromJSON(config_path)

  # Filtering thresholds
  pkg_env$min_num_years           <- cfg$filtering$min_num_years
  pkg_env$min_frac_good_data      <- cfg$filtering$min_frac_good_data
  pkg_env$min_q_value             <- cfg$filtering$min_q_value
  pkg_env$min_days_above_threshold <- cfg$filtering$min_days_above_threshold
  pkg_env$min_nona_days_annual    <- cfg$filtering$min_nona_days_annual

  # Water year
  pkg_env$water_year_start_month  <- cfg$water_year$start_month

  # Baseflow
  pkg_env$eckhardt_bfimax         <- cfg$baseflow$eckhardt_bfimax
  pkg_env$eckhardt_alpha          <- cfg$baseflow$eckhardt_alpha
  pkg_env$lyne_hollick_alpha      <- cfg$baseflow$lyne_hollick_alpha
  pkg_env$lyne_hollick_passes     <- cfg$baseflow$lyne_hollick_passes
  # Recession
  pkg_env$recession_min_length    <- cfg$recession$min_length
  pkg_env$recession_min_events    <- cfg$recession$min_events

  # Pulses
  pkg_env$high_pulse_percentile   <- cfg$pulses$high_percentile
  pkg_env$low_pulse_percentile    <- cfg$pulses$low_percentile
  pkg_env$flow_reversal_threshold <- cfg$pulses$flow_reversal_threshold

  # Timing
  pkg_env$d_percentiles           <- cfg$timing$d_percentiles

  # FDC
  pkg_env$fdc_flow_floor          <- 1e-10

  # Flow volumes
  pkg_env$flow_volumes_percentiles <- cfg$flow_volumes$percentiles

  # Elasticity
  pkg_env$elasticity_window_years <- cfg$elasticity$window_years
  pkg_env$elasticity_min_years    <- cfg$elasticity$min_years
  pkg_env$elasticity_min_annual_ppt <- cfg$elasticity$min_annual_ppt

  # QP seasonality
  pkg_env$qp_slope_window_days    <- cfg$qp_seasonality$slope_window_days
  pkg_env$qp_min_years            <- cfg$qp_seasonality$min_years
  # Runoff ratios
  pkg_env$runoff_min_annual_ppt   <- cfg$runoff_ratios$min_annual_ppt
  pkg_env$runoff_min_seasonal_ppt <- cfg$runoff_ratios$min_seasonal_ppt

  # Storage
  pkg_env$storage_min_years       <- cfg$storage$min_years
  # QA/QC
  pkg_env$qa_qann_range           <- cfg$qa_qc$qann_range
  pkg_env$qa_bfi_range            <- cfg$qa_qc$bfi_range
  pkg_env$qa_flashiness_range     <- cfg$qa_qc$flashiness_range
  pkg_env$qa_tqmean_range         <- cfg$qa_qc$tqmean_range
  pkg_env$qa_d50_range            <- cfg$qa_qc$d50_range
  pkg_env$qa_elasticity_range     <- cfg$qa_qc$elasticity_range
  pkg_env$qa_runoff_ratio_range   <- cfg$qa_qc$runoff_ratio_range
  pkg_env$qa_seasonal_sum_tolerance <- cfg$qa_qc$seasonal_sum_tolerance
  pkg_env$qa_max_na_fraction      <- cfg$qa_qc$max_na_fraction

  # Metadata
  pkg_env$metadata_gages_ii_dir   <- cfg$metadata$gages_ii_dir
  pkg_env$metadata_include_human_interference <- cfg$metadata$include_human_interference
  pkg_env$metadata_hydat_path     <- cfg$metadata$hydat_path
  pkg_env$metadata_interference_columns <- cfg$metadata$interference_columns

  # Stat suffixes
  pkg_env$stat_suffixes <- c(
    "_senn_slp", "_linear_slp",
    "_spearman_rho", "_spearman_pval",
    "_mk_rho", "_mk_pval",
    "_mean", "_median"
  )

  # NA handling (from na_handling section)
  nah <- cfg$na_handling
  if (!is.null(nah)) {
    pkg_env$na_max_gap_days            <- nah$interpolation$max_gap_days %||% 3L
    pkg_env$na_interpolation_method     <- nah$interpolation$method %||% "linear"
    pkg_env$na_internal_only            <- nah$interpolation$internal_only %||% TRUE
    pkg_env$na_max_raw_na_per_year      <- nah$year_rejection$max_raw_na_per_year %||% 30L
    pkg_env$na_reject_negative_flow     <- nah$year_rejection$reject_negative_flow %||% TRUE
    pkg_env$na_reject_residual_na       <- nah$year_rejection$reject_residual_na %||% TRUE
    pkg_env$na_constant_sd_enabled      <- nah$constant_sd_flag$enabled %||% TRUE
    pkg_env$na_constant_sd_min_days     <- nah$constant_sd_flag$min_nonzero_days_per_month %||% 15L
    pkg_env$na_constant_sd_max_unique   <- nah$constant_sd_flag$max_unique_values %||% 1L
    pkg_env$na_trend_min_fraction       <- nah$trend_completeness$min_fraction %||% 0.80
    pkg_env$na_decade_min_fraction      <- nah$trend_completeness$decade_min_fraction %||% 0.80
    pkg_env$na_seasonal_min_fraction    <- nah$seasonal_completeness$min_fraction %||% 0.80
    pkg_env$na_seasonal_use_raw         <- nah$seasonal_completeness$use_raw_observations %||% TRUE
    pkg_env$na_max_raw_na_ppt           <- nah$climate_na_policy$max_raw_na_per_year_ppt %||% 30L
    pkg_env$na_max_gap_ppt              <- nah$climate_na_policy$max_interpolation_gap_ppt %||% 3L
    pkg_env$na_reject_negative_ppt      <- nah$climate_na_policy$reject_negative_ppt %||% TRUE
  } else {
    # Defaults if na_handling section not present
    pkg_env$na_max_gap_days            <- 3L
    pkg_env$na_interpolation_method     <- "linear"
    pkg_env$na_internal_only            <- TRUE
    pkg_env$na_max_raw_na_per_year      <- 30L
    pkg_env$na_reject_negative_flow     <- TRUE
    pkg_env$na_reject_residual_na       <- TRUE
    pkg_env$na_constant_sd_enabled      <- TRUE
    pkg_env$na_constant_sd_min_days     <- 15L
    pkg_env$na_constant_sd_max_unique   <- 1L
    pkg_env$na_trend_min_fraction       <- 0.80
    pkg_env$na_decade_min_fraction      <- 0.80
    pkg_env$na_seasonal_min_fraction    <- 0.80
    pkg_env$na_seasonal_use_raw         <- TRUE
    pkg_env$na_max_raw_na_ppt           <- 30L
    pkg_env$na_max_gap_ppt              <- 3L
    pkg_env$na_reject_negative_ppt      <- TRUE
  }

  # Legacy filtering flag
  pkg_env$use_legacy_filtering <- cfg$filtering$use_legacy_filtering %||% TRUE

  invisible(NULL)
}

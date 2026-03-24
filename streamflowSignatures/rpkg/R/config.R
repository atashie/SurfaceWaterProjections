# Package-internal environment for configuration values
pkg_env <- new.env(parent = emptyenv())

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
  pkg_env$baseflow_min_days       <- cfg$baseflow$min_days
  pkg_env$baseflow_max_missing_frac <- cfg$baseflow$max_missing_frac

  # Recession
  pkg_env$recession_min_length    <- cfg$recession$min_length
  pkg_env$recession_min_events    <- cfg$recession$min_events

  # Pulses
  pkg_env$high_pulse_percentile   <- cfg$pulses$high_percentile
  pkg_env$low_pulse_percentile    <- cfg$pulses$low_percentile
  pkg_env$flow_reversal_threshold <- cfg$pulses$flow_reversal_threshold
  pkg_env$pulses_min_days         <- cfg$pulses$min_days

  # Timing
  pkg_env$timing_min_days         <- cfg$timing$min_days
  pkg_env$d_percentiles           <- cfg$timing$d_percentiles

  # FDC
  pkg_env$fdc_min_days            <- cfg$fdc$min_days
  pkg_env$fdc_flow_floor          <- 1e-10

  # Flashiness
  pkg_env$flashiness_min_days     <- cfg$flashiness$min_days
  pkg_env$flashiness_max_missing_frac <- cfg$flashiness$max_missing_frac

  # Flow volumes
  pkg_env$flow_volumes_min_days   <- cfg$flow_volumes$min_days
  pkg_env$flow_volumes_percentiles <- cfg$flow_volumes$percentiles

  # Elasticity
  pkg_env$elasticity_window_years <- cfg$elasticity$window_years
  pkg_env$elasticity_min_years    <- cfg$elasticity$min_years
  pkg_env$elasticity_min_data_completeness <- cfg$elasticity$min_data_completeness
  pkg_env$elasticity_min_annual_ppt <- cfg$elasticity$min_annual_ppt

  # QP seasonality
  pkg_env$qp_slope_window_days    <- cfg$qp_seasonality$slope_window_days
  pkg_env$qp_min_years            <- cfg$qp_seasonality$min_years
  pkg_env$qp_min_days_per_year    <- cfg$qp_seasonality$min_days_per_year
  pkg_env$qp_max_na_frac          <- cfg$qp_seasonality$max_na_frac

  # Runoff ratios
  pkg_env$runoff_min_annual_ppt   <- cfg$runoff_ratios$min_annual_ppt
  pkg_env$runoff_min_seasonal_ppt <- cfg$runoff_ratios$min_seasonal_ppt

  # Storage
  pkg_env$storage_min_years       <- cfg$storage$min_years
  pkg_env$storage_min_days_per_year <- cfg$storage$min_days_per_year
  pkg_env$storage_max_na_frac     <- cfg$storage$max_na_frac

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

  invisible(NULL)
}

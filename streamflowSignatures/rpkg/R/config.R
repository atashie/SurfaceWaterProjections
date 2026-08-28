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
    pkg_env$na_reject_negative_flow     <- nah$year_rejection$reject_negative_flow %||% FALSE  # key-missing fallback aligned with canonical (FALSE, 2026-08-24; was TRUE). NOTE: the section-absent branch below still defaults this to TRUE — a known deferred fix tracked in CHANGELOG [Unreleased].
    # SWE NA policy (mirrors PPT; fallbacks match julia/src/config.jl)
    pkg_env$na_max_raw_na_swe   <- nah$snow_na_policy$max_raw_na_per_year_swe %||% 30L
    pkg_env$na_max_gap_swe      <- nah$snow_na_policy$max_interpolation_gap_swe %||% 3L
    pkg_env$na_reject_negative_swe <- nah$snow_na_policy$reject_negative_swe %||% TRUE
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

  # ---- Changepoint (Pettitt); fallbacks mirror julia/src/config.jl
  cp <- cfg$changepoint
  pkg_env$changepoint_enabled   <- if (!is.null(cp$enabled)) isTRUE(cp$enabled) else FALSE
  pkg_env$cp_start_water_year   <- if (!is.null(cp$start_water_year)) cp$start_water_year else 1980L
  pkg_env$cp_end_water_year     <- if (!is.null(cp$end_water_year)) cp$end_water_year else 2024L
  pkg_env$cp_min_total_obs      <- if (!is.null(cp$min_total_obs)) cp$min_total_obs else 20L
  pkg_env$cp_min_segment_obs    <- if (!is.null(cp$min_segment_obs)) cp$min_segment_obs else 10L

  # ---- Stats floor: NULL when the section is absent (no floor; legacy behavior)
  pkg_env$min_values_for_stats <- if (!is.null(cfg$stats_floor$min_values_for_stats))
    as.integer(cfg$stats_floor$min_values_for_stats) else NULL

  # ---- Snow (fallbacks mirror julia/src/config.jl; absent section => gate off)
  sn <- cfg$snow
  pkg_env$snow_swe_threshold_mm   <- if (!is.null(sn$swe_day_threshold_mm)) sn$swe_day_threshold_mm else 10
  pkg_env$snow_seasonal_min_days  <- if (!is.null(sn$seasonal_spell_min_days)) sn$seasonal_spell_min_days else 60L
  pkg_env$snow_melt_com_fraction  <- if (!is.null(sn$melt_com_fraction)) sn$melt_com_fraction else 0.5
  pkg_env$snow_min_annual_ppt_mm  <- if (!is.null(sn$min_annual_ppt_mm)) sn$min_annual_ppt_mm else 10
  pkg_env$snow_record_decade_gate <- if (!is.null(sn$record_anchored_decade_gate)) isTRUE(sn$record_anchored_decade_gate) else FALSE

  # ---- Drought (absent section => family disabled)
  dr <- cfg$drought
  pkg_env$drought_enabled          <- !is.null(dr) && length(dr) > 0
  pkg_env$drought_smooth_window    <- if (!is.null(dr$smoothing_window_days)) dr$smoothing_window_days else 7L
  pkg_env$drought_smooth_alignment <- if (!is.null(dr$smoothing_alignment)) dr$smoothing_alignment else "center"
  pkg_env$drought_smooth_min_valid <- if (!is.null(dr$smoothing_min_valid_days)) dr$smoothing_min_valid_days else 4L
  pkg_env$drought_percentiles      <- if (!is.null(dr$threshold_percentiles)) as.integer(dr$threshold_percentiles) else c(2L,5L,10L,20L,30L)
  pkg_env$drought_below_range_policy <- if (!is.null(dr$below_plotting_range_policy)) dr$below_plotting_range_policy else "na"
  pkg_env$drought_min_years        <- if (!is.null(dr$min_years_for_threshold)) dr$min_years_for_threshold else 10L
  if (isTRUE(pkg_env$drought_enabled)) {
    if (!is.null(dr$threshold_method) && dr$threshold_method != "fixed")
      stop("drought.threshold_method must be 'fixed' (variable method not implemented)")
    if (!is.null(dr$plotting_position) && dr$plotting_position != "weibull")
      stop("drought.plotting_position must be 'weibull'")
    if (!pkg_env$drought_below_range_policy %in% c("na", "clamp"))
      stop("drought.below_plotting_range_policy must be 'na' or 'clamp'")
    if (pkg_env$drought_smooth_alignment == "center" && pkg_env$drought_smooth_window %% 2 == 0)
      stop("drought.smoothing_window_days must be ODD for centered alignment")
    if (pkg_env$drought_smooth_min_valid > pkg_env$drought_smooth_window)
      stop("drought.smoothing_min_valid_days must be <= smoothing_window_days")
    if (!pkg_env$drought_smooth_alignment %in% c("center", "trailing"))
      stop("drought.smoothing_alignment must be 'center' or 'trailing'")
    if (pkg_env$drought_smooth_window < 1L)
      stop("drought.smoothing_window_days must be >= 1")
    if (pkg_env$drought_smooth_min_valid < 1L)
      stop("drought.smoothing_min_valid_days must be >= 1")
    pc <- pkg_env$drought_percentiles
    if (length(pc) == 0L || any(pc <= 0) || any(pc >= 100))
      stop("drought.threshold_percentiles must be non-empty and all in (0, 100)")
    if (is.unsorted(pc)) stop("drought.threshold_percentiles must be sorted ascending")
    if (anyDuplicated(pc) > 0) stop("drought.threshold_percentiles must be unique")
    if (pkg_env$drought_min_years < 1L)
      stop("drought.min_years_for_threshold must be >= 1")
  }

  # ---- Annual values export
  pkg_env$save_annual_values <- if (!is.null(cfg$annual_values$save)) isTRUE(cfg$annual_values$save) else FALSE

  invisible(NULL)

}

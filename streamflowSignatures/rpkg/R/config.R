# Package-internal environment for configuration values
pkg_env <- new.env(parent = emptyenv())

# %||% operator for NULL-coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Strict scalar readers mirroring Julia's typed constants (review 2026-10-12): an integer
# setting must be a finite whole number (no truncation), a number finite, a flag a real
# logical. Used for the climate section and the temperature NA policy.
.cfg_int <- function(v, key, default) {
  if (is.null(v)) return(as.integer(default))
  if (is.logical(v) || !is.numeric(v) || length(v) != 1L || !is.finite(v) || v != floor(v))
    stop(sprintf("config %s: expected an integer, got %s", key, paste(format(v), collapse = ",")))
  as.integer(v)
}
.cfg_num <- function(v, key, default) {
  if (is.null(v)) return(as.numeric(default))
  if (is.logical(v) || !is.numeric(v) || length(v) != 1L || !is.finite(v))
    stop(sprintf("config %s: expected a finite number, got %s", key, paste(format(v), collapse = ",")))
  as.numeric(v)
}
.cfg_bool <- function(v, key, default) {
  if (is.null(v)) return(isTRUE(default))
  if (!is.logical(v) || length(v) != 1L || is.na(v))
    stop(sprintf("config %s: expected true/false, got %s", key, paste(format(v), collapse = ",")))
  v
}
# Observation contract shared with Julia's _real_vector (review round 4): a value is an
# observation when it is a real number (logical included, as Julia's Bool <: Real);
# character, complex, list and anything else are missing. Never coerces strings.
.real_values <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x)) return(as.numeric(x))
  if (is.list(x)) return(vapply(x, function(v) if (length(v) == 1L && (is.numeric(v) || is.logical(v))) as.numeric(v) else NA_real_, numeric(1)))
  rep(NA_real_, length(x))   # character, complex, factor, ...
}
# Season months: whole calendar months 1..12, no logicals (used by .onLoad and by
# preprocess_daily_data for explicit overrides — review round 3)
.cfg_months <- function(v, key, default) {
  if (is.null(v)) return(as.integer(default))
  if (is.logical(v) || !is.numeric(v) || length(v) == 0L || any(!is.finite(v)) || any(v != floor(v)) || any(v < 1 | v > 12))
    stop(sprintf("season_definitions.%s: expected calendar months 1..12", key))
  as.integer(v)
}

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
  # flagged_for_high_na denominator manifest (2026-09-04; mirrors julia/src/config.jl —
  # an absent section falls back to the same lists).
  hn <- cfg$qa_qc$high_na_denominator
  pkg_env$qa_high_na_suffixes <- if (is.null(hn$statistic_suffixes)) c(
    "_senn_slp",
    "_linear_slp",
    "_spearman_rho",
    "_spearman_pval",
    "_mk_rho",
    "_mk_pval",
    "_mean",
    "_median",
    "_pettitt_cp_year",
    "_pettitt_pval",
    "_pettitt_pre_mean",
    "_pettitt_post_mean",
    "_pettitt_delta_mean",
    "_pettitt_pct_change",
    "_pettitt_pre_mk_pval",
    "_pettitt_post_mk_pval") else as.character(hn$statistic_suffixes)
  pkg_env$qa_high_na_scalars <- if (is.null(hn$per_gage_scalars)) c(
    "elasticity_static",
    "log_a_seasonality_amplitude_all",
    "log_a_seasonality_amplitude_first_half",
    "log_a_seasonality_amplitude_last_half",
    "log_a_seasonality_minimum_all",
    "log_a_seasonality_minimum_first_half",
    "log_a_seasonality_minimum_last_half",
    "runoff_ratio_high_count",
    "elasticity_years_total",
    "elasticity_years_low_ppt",
    "ice_affected_days_total",
    "recession_alpha_point_cloud_linear_reservoir",
    "season_excluded_years_winter",
    "season_excluded_years_spring",
    "season_excluded_years_summer",
    "season_excluded_years_fall",
    "drought_threshold_fixed_p2",
    "drought_threshold_fixed_p5",
    "drought_threshold_fixed_p10",
    "drought_threshold_fixed_p20",
    "drought_threshold_fixed_p30") else as.character(hn$per_gage_scalars)

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
    pkg_env$na_max_raw_na_swe   <- .cfg_int(nah$snow_na_policy$max_raw_na_per_year_swe, "max_raw_na_per_year_swe", 30L)
    pkg_env$na_max_gap_swe      <- .cfg_int(nah$snow_na_policy$max_interpolation_gap_swe, "max_interpolation_gap_swe", 3L)
    pkg_env$na_reject_negative_swe <- .cfg_bool(nah$snow_na_policy$reject_negative_swe, "reject_negative_swe", TRUE)
    if (pkg_env$na_max_raw_na_swe < 0L || pkg_env$na_max_gap_swe < 0L) stop("na_handling.snow_na_policy limits must be >= 0")
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
    # Temperature (tmin/tmax) policy for the climate metrics (2026-10-12): the PPT/SWE
    # shape without a negative-value check (fallbacks mirror julia/src/config.jl)
    pkg_env$na_max_raw_na_temp          <- .cfg_int(nah$temperature_na_policy$max_raw_na_per_year_temp, "max_raw_na_per_year_temp", 30L)
    pkg_env$na_max_gap_temp             <- .cfg_int(nah$temperature_na_policy$max_interpolation_gap_temp, "max_interpolation_gap_temp", 3L)
    if (pkg_env$na_max_raw_na_temp < 0L || pkg_env$na_max_gap_temp < 0L)
      stop("na_handling.temperature_na_policy: max_raw_na_per_year_temp and max_interpolation_gap_temp must be >= 0")
    # season definitions (calendar months per season), as Julia's CFG_NA_SEASON_MONTHS;
    # whole months 1..12 only (review round 2: no truncation of 1.9 -> 1)
    sd <- nah$seasonal_completeness$season_definitions
    pkg_env$na_season_definitions <- list(
      winter = .cfg_months(sd$winter, "winter", c(12L, 1L, 2L)), spring = .cfg_months(sd$spring, "spring", c(3L, 4L, 5L)),
      summer = .cfg_months(sd$summer, "summer", c(6L, 7L, 8L)), fall = .cfg_months(sd$fall, "fall", c(9L, 10L, 11L)))
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
    pkg_env$na_max_raw_na_temp          <- 30L
    pkg_env$na_max_gap_temp             <- 3L
    pkg_env$na_season_definitions <- list(winter = c(12L, 1L, 2L), spring = c(3L, 4L, 5L),
                                          summer = c(6L, 7L, 8L), fall = c(9L, 10L, 11L))
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

  # ---- Climate (guidelines 3.9, 2026-10-12; fallbacks mirror julia/src/config.jl).
  # An ABSENT section disables the family entirely (no climate columns).
  # An explicit JSON null is not a value: jsonlite turns it into an absent key so the default
  # would silently apply, while Julia rejects it (Sonnet review 2026-10-13). Checked on the
  # unsimplified JSON for the climate section and every na_handling sub-section.
  .reject_json_nulls <- function(x, path) {
    if (!is.list(x) || is.null(names(x))) return(invisible(NULL))
    for (key in names(x)) {
      if (is.null(x[[key]])) stop(sprintf("%s.%s must not be null", path, key))
      if (is.list(x[[key]])) .reject_json_nulls(x[[key]], paste(path, key, sep = "."))
    }
  }
  raw_all <- tryCatch(jsonlite::fromJSON(config_path, simplifyVector = FALSE), error = function(e) NULL)
  if (!is.null(raw_all)) {
    .reject_json_nulls(raw_all$climate, "climate")
    .reject_json_nulls(raw_all$na_handling, "na_handling")
  }
  cl <- cfg$climate
  pkg_env$climate_enabled <- !is.null(cl) && length(cl) > 0
  pkg_env$climate_dry_day_mm <- .cfg_num(cl$dry_day_threshold_mm, "dry_day_threshold_mm", 0.1)
  wet <- cl$wet_day_thresholds_mm %||% c(10, 25, 50)
  if (is.logical(wet) || !is.numeric(wet) || any(!is.finite(wet))) stop("climate.wet_day_thresholds_mm must be a list of finite numbers")
  pkg_env$climate_wet_thresholds_mm <- as.numeric(wet)
  pkg_env$climate_multiday_window <- .cfg_int(cl$multi_day_max_window_days, "multi_day_max_window_days", 5L)
  cum_fraction_raw <- .cfg_num(cl$cumulative_fraction, "cumulative_fraction", 0.5)
  pkg_env$climate_cum_fraction <- cum_fraction_raw
  pkg_env$climate_frost_temp_c <- .cfg_num(cl$frost_temp_c, "frost_temp_c", 0)
  pkg_env$climate_frost_last_end_month <- .cfg_int(cl$frost_last_end_month, "frost_last_end_month", 6L)
  pkg_env$climate_gdd_base_c <- .cfg_num(cl$gdd_base_c, "gdd_base_c", 0)
  pkg_env$climate_record_decade_gate <- .cfg_bool(cl$record_anchored_decade_gate, "record_anchored_decade_gate", FALSE)
  if (isTRUE(pkg_env$climate_enabled)) {
    thr <- pkg_env$climate_wet_thresholds_mm
    if (pkg_env$climate_dry_day_mm < 0) stop("climate.dry_day_threshold_mm must be >= 0")
    if (length(thr) == 0 || any(thr <= 0) || is.unsorted(thr, strictly = TRUE))
      stop("climate.wet_day_thresholds_mm must be a non-empty, ascending, unique list of positive values")
    if (anyDuplicated(precip_metrics()) > 0)
      stop("climate: the configured window / thresholds produce duplicate precipitation column names")
    if (pkg_env$climate_multiday_window < 2L || pkg_env$climate_multiday_window > 366L)
      stop("climate.multi_day_max_window_days must lie in 2..366")
    if (!(cum_fraction_raw > 0 && cum_fraction_raw < 1))
      stop("climate.cumulative_fraction must lie strictly inside (0, 1)")
    pct <- round(100 * cum_fraction_raw, 9)
    if (abs(pct - round(pct)) > 0) stop("climate.cumulative_fraction must be a whole percentage (e.g. 0.5 -> prcp_50_dowy)")
    # the fraction USED is the canonical whole percentage the column name shows (mirrors Julia)
    pkg_env$climate_cum_fraction <- round(100 * cum_fraction_raw) / 100
    # JSON types before jsonlite's simplification: a Boolean inside the threshold list or a
    # one-element array for a scalar must fail (review round 2)
    raw_cfg <- tryCatch(jsonlite::fromJSON(config_path, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(raw_cfg) && !is.null(raw_cfg$climate)) {
      rc <- raw_cfg$climate
      for (key in c("dry_day_threshold_mm", "multi_day_max_window_days", "cumulative_fraction", "frost_temp_c",
                    "frost_last_end_month", "gdd_base_c", "record_anchored_decade_gate")) {
        if (!is.null(rc[[key]]) && is.list(rc[[key]])) stop(sprintf("climate.%s must be a scalar, not an array", key))
      }
      if (!is.null(rc$wet_day_thresholds_mm)) {
        for (x in rc$wet_day_thresholds_mm) if (is.logical(x) || !is.numeric(x)) stop("climate.wet_day_thresholds_mm must contain numbers only")
      }
      tp <- raw_cfg$na_handling$temperature_na_policy
      for (key in c("max_raw_na_per_year_temp", "max_interpolation_gap_temp")) {
        if (!is.null(tp[[key]]) && is.list(tp[[key]])) stop(sprintf("na_handling.temperature_na_policy.%s must be a scalar", key))
      }
      # season lists must hold numbers only (jsonlite turns [true, 2] into c(1, 2))
      sdr <- raw_cfg$na_handling$seasonal_completeness$season_definitions
      for (key in names(sdr)) for (x in sdr[[key]]) if (is.logical(x) || !is.numeric(x))
        stop(sprintf("na_handling.seasonal_completeness.season_definitions.%s must contain numbers only", key))
    }
    # the canonical percentage must stay strictly inside the open interval (review round 3)
    if (!(round(100 * pkg_env$climate_cum_fraction) >= 1 && round(100 * pkg_env$climate_cum_fraction) <= 99))
      stop("climate.cumulative_fraction must canonicalise to a whole percentage between 1 and 99")
    if (pkg_env$climate_frost_last_end_month < 1L || pkg_env$climate_frost_last_end_month > 12L)
      stop("climate.frost_last_end_month must be a calendar month 1-12")
  }

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

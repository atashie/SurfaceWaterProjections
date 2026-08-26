#' Calculate all streamflow signatures for a single gage
#'
#' Orchestrator function that calls all signature modules and merges results
#' into a single named list. Mirrors Julia's \code{calculate_all_signatures()}.
#'
#' @param streamflow_data A data.frame/data.table with columns: water_year, Q,
#'   month, dowy. For climate signatures, also needs PPT.
#' @param has_climate Logical. If TRUE and PPT column exists, compute climate
#'   signatures (runoff ratios, elasticity, Q-P seasonality, avg storage).
#' @param seasonal_flags Optional data.table from \code{preprocess_daily_data()}
#'   with per-year seasonal completeness flags. Passed to
#'   \code{calculate_flow_vols_by_year()} and \code{analyze_Q_PPT_relationships()}.
#' @param trend_completeness Numeric in (0,1] or NULL. Forwarded to
#'   \code{generate_stats()} in all signature functions.
#' @param decade_completeness Numeric in (0,1] or NULL. Forwarded to
#'   \code{generate_stats()} in all signature functions.
#' @param climate_data Optional data.frame/data.table with climate columns (PPT).
#'   If provided, used instead of \code{streamflow_data} for climate signatures.
#' @param include_qa_flags Logical. If TRUE, append 12 QA/QC flag columns from
#'   \code{compute_qa_flags()} to the output. Default FALSE.
#' @param area_normalized Logical. Whether Q is area-normalized to mm/day.
#'   When FALSE (gages with no drainage area — Q left in raw m3/s), all
#'   Q-to-PPT signatures (runoff ratios, elasticity, Q-P seasonality, storage)
#'   are skipped because Q and PPT units don't match. Q-only signatures are
#'   unaffected. Default TRUE.
#' @return Named list with ~1,653 elements for a full climate+SWE gage (8 stats
#'   + 8 Pettitt fields per signature base, plus per-gage scalars); fewer when
#'   climate, SWE or the drought family are unavailable/disabled,
#'   plus 12 flag columns if \code{include_qa_flags = TRUE}.
#' @export
calculate_all_signatures <- function(streamflow_data, has_climate = FALSE,
                                     seasonal_flags = NULL,
                                     trend_completeness = NULL,
                                     decade_completeness = NULL,
                                     climate_data = NULL,
                                     include_qa_flags = FALSE,
                                     area_normalized = TRUE,
                                     changepoint = NULL,
                                     min_values_for_stats = NULL,
                                     collector = NULL,
                                     snow_data = NULL,
                                     snow_climate_years = NULL,
                                     gage_id = NULL) {
  results <- list()

  # Season exclusion year counts (per-gage scalar diagnostics)
  if (!is.null(seasonal_flags) && nrow(seasonal_flags) > 0) {
    for (pair in list(
      list(season = "winter", col = "winter_complete"),
      list(season = "spring", col = "spring_complete"),
      list(season = "summer", col = "summer_complete"),
      list(season = "fall",   col = "fall_complete")
    )) {
      key <- paste0("season_excluded_years_", pair$season)
      if (pair$col %in% names(seasonal_flags)) {
        results[[key]] <- as.numeric(sum(!seasonal_flags[[pair$col]]))
      } else {
        results[[key]] <- 0
      }
    }
  }

  # Safe wrapper for calling signature functions
  gid <- if (is.null(gage_id)) "unknown" else as.character(gage_id)
  safe_call <- function(fn, label, ...) {
    tryCatch(fn(...), error = function(e) {
      warning(paste0("Gage ", gid, ": ", label, " failed: ", e$message))
      NULL
    })
  }

  # Flow volumes (pass seasonal_flags)
  out <- safe_call(calculate_flow_vols_by_year, "Flow volumes",
                   streamflow_data, seasonal_flags = seasonal_flags,
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness,
                   min_values_for_stats = min_values_for_stats,
                   changepoint = changepoint, collector = collector)
  if (!is.null(out)) results <- c(results, out)

  # Non-climate signatures (no seasonal_flags needed)
  specs <- list(
    list(fn = analyze_fdc_trends,           label = "FDC trends"),
    list(fn = analyze_flashiness_trends,    label = "Flashiness"),
    list(fn = analyze_flow_timing_trends,   label = "Flow timing"),
    list(fn = calculate_pulse_metrics,      label = "Pulse metrics"),
    list(fn = analyze_baseflow_indices,     label = "Baseflow indices")
  )

  for (spec in specs) {
    out <- safe_call(spec$fn, spec$label, streamflow_data,
                     trend_completeness = trend_completeness,
                     decade_completeness = decade_completeness,
                   min_values_for_stats = min_values_for_stats,
                   changepoint = changepoint, collector = collector)
    if (!is.null(out)) results <- c(results, out)
  }

  # Recession: event-based (inherently sparse), no trend completeness
  recession_alpha <- NaN
  out <- safe_call(analyze_recession_parameters, "Recession parameters",
                   streamflow_data, changepoint = changepoint, collector = collector)
  if (!is.null(out)) {
    # Extract scalar for parameterized BFI before merging
    if ("recession_alpha_point_cloud_linear_reservoir" %in% names(out)) {
      recession_alpha <- out[["recession_alpha_point_cloud_linear_reservoir"]]
    }
    results <- c(results, out)
  }

  # Parameterized BFI using recession-derived alpha (requires recession to have run first)
  # Uses trend completeness (same as fixed-parameter BFI)
  out <- safe_call(analyze_baseflow_indices_with_parameters, "Parameterized baseflow",
                   streamflow_data, alpha = recession_alpha,
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness,
                   min_values_for_stats = min_values_for_stats,
                   changepoint = changepoint, collector = collector)
  if (!is.null(out)) results <- c(results, out)

  # Negative days
  out <- safe_call(calculate_negative_days, "Negative days",
                   streamflow_data,
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness,
                   min_values_for_stats = min_values_for_stats,
                   changepoint = changepoint, collector = collector)
  if (!is.null(out)) results <- c(results, out)

  # Streamflow drought (config-gated: absent `drought` section => family off)
  if (isTRUE(pkg_env$drought_enabled)) {
    out <- safe_call(calculate_drought_metrics, "Drought",
                     streamflow_data,
                     trend_completeness = trend_completeness,
                     decade_completeness = decade_completeness,
                     min_values_for_stats = min_values_for_stats,
                     changepoint = changepoint, collector = collector)
    if (!is.null(out)) results <- c(results, out)
  }

  # Climate signatures
  # Gate: Q-to-PPT signatures are undefined when Q is not area-normalized
  # (raw m3/s vs PPT in mm — units don't match, so Q/P ratios, dQ/dP, and
  # cumsum(P - Q) are all meaningless).
  if (has_climate && !isTRUE(area_normalized)) {
    message("Skipping Q-to-PPT signatures (area_normalized = FALSE — Q in raw m3/s, PPT in mm)")
  }
  clim_input <- if (!is.null(climate_data)) climate_data else streamflow_data
  if (has_climate && isTRUE(area_normalized) && "PPT" %in% colnames(clim_input)) {
    # Q-PPT ratios (pass seasonal_flags)
    out <- safe_call(analyze_Q_PPT_relationships, "Q-PPT ratios",
                     clim_input, seasonal_flags = seasonal_flags,
                     trend_completeness = trend_completeness,
                     decade_completeness = decade_completeness,
                   min_values_for_stats = min_values_for_stats,
                   changepoint = changepoint, collector = collector)
    if (!is.null(out)) results <- c(results, out)

    # Elasticity: rolling-window based, no trend completeness
    out <- safe_call(calculate_streamflow_elasticity, "Elasticity",
                     clim_input, changepoint = changepoint, collector = collector)
    if (!is.null(out)) results <- c(results, out)

    climate_specs <- list(
      list(fn = calculate_qp_seasonality,         label = "Q-P seasonality"),
      list(fn = calculate_average_storage,        label = "Average storage")
    )

    for (spec in climate_specs) {
      out <- safe_call(spec$fn, spec$label, clim_input,
                       trend_completeness = trend_completeness,
                       decade_completeness = decade_completeness,
                   min_values_for_stats = min_values_for_stats,
                   changepoint = changepoint, collector = collector)
      if (!is.null(out)) results <- c(results, out)
    }
  }

  # Snow metrics: EXPLICIT opt-in frame only — an SWE column in streamflow_data
  # is never used implicitly (matching Julia; prevents SWE-invalid years leaking)
  if (!is.null(snow_data) && "SWE" %in% names(snow_data)) {
    out <- safe_call(calculate_snow_metrics, "Snow",
                     snow_data, valid_climate_years = snow_climate_years,
                     trend_completeness = trend_completeness,
                     decade_completeness = decade_completeness,
                     min_values_for_stats = min_values_for_stats,
                     changepoint = changepoint, collector = collector)
    if (!is.null(out)) results <- c(results, out)
  }

  # QA/QC flags (optional)
  if (isTRUE(include_qa_flags)) {
    qa_flags <- tryCatch(compute_qa_flags(results), error = function(e) {
      warning(paste("QA/QC flags failed:", e$message))
      NULL
    })
    if (!is.null(qa_flags)) results <- c(results, qa_flags)
  }

  results
}

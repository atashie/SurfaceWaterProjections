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
#' @return Named list with ~478 (no climate) or ~551 (with climate) elements.
#' @export
calculate_all_signatures <- function(streamflow_data, has_climate = FALSE,
                                     seasonal_flags = NULL,
                                     trend_completeness = NULL,
                                     decade_completeness = NULL,
                                     climate_data = NULL) {
  results <- list()

  # Safe wrapper for calling signature functions
  safe_call <- function(fn, label, ...) {
    tryCatch(fn(...), error = function(e) {
      warning(paste(label, "failed:", e$message))
      NULL
    })
  }

  # Flow volumes (pass seasonal_flags)
  out <- safe_call(calculate_flow_vols_by_year, "Flow volumes",
                   streamflow_data, seasonal_flags = seasonal_flags,
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness)
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
                     decade_completeness = decade_completeness)
    if (!is.null(out)) results <- c(results, out)
  }

  # Recession: event-based (inherently sparse), no trend completeness
  out <- safe_call(analyze_recession_parameters, "Recession parameters",
                   streamflow_data)
  if (!is.null(out)) results <- c(results, out)

  # Negative days
  out <- safe_call(calculate_negative_days, "Negative days",
                   streamflow_data,
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness)
  if (!is.null(out)) results <- c(results, out)

  # Climate signatures
  clim_input <- if (!is.null(climate_data)) climate_data else streamflow_data
  if (has_climate && "PPT" %in% colnames(clim_input)) {
    # Q-PPT ratios (pass seasonal_flags)
    out <- safe_call(analyze_Q_PPT_relationships, "Q-PPT ratios",
                     clim_input, seasonal_flags = seasonal_flags,
                     trend_completeness = trend_completeness,
                     decade_completeness = decade_completeness)
    if (!is.null(out)) results <- c(results, out)

    # Elasticity: rolling-window based, no trend completeness
    out <- safe_call(calculate_streamflow_elasticity, "Elasticity",
                     clim_input)
    if (!is.null(out)) results <- c(results, out)

    climate_specs <- list(
      list(fn = calculate_qp_seasonality,         label = "Q-P seasonality"),
      list(fn = calculate_average_storage,        label = "Average storage")
    )

    for (spec in climate_specs) {
      out <- safe_call(spec$fn, spec$label, clim_input,
                       trend_completeness = trend_completeness,
                       decade_completeness = decade_completeness)
      if (!is.null(out)) results <- c(results, out)
    }
  }

  results
}

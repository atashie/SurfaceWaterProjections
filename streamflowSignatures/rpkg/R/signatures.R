#' Calculate all streamflow signatures for a single gage
#'
#' Orchestrator function that calls all signature modules and merges results
#' into a single named list. Mirrors Julia's \code{calculate_all_signatures()}.
#'
#' @param streamflow_data A data.frame/data.table with columns: water_year, Q,
#'   month, dowy. For climate signatures, also needs PPT.
#' @param has_climate Logical. If TRUE and PPT column exists, compute climate
#'   signatures (runoff ratios, elasticity, Q-P seasonality, avg storage).
#' @return Named list with ~478 (no climate) or ~551 (with climate) elements.
#' @export
calculate_all_signatures <- function(streamflow_data, has_climate = FALSE) {
  results <- list()

  # Non-climate signatures
  safe_call <- function(fn, label) {
    tryCatch(fn(streamflow_data), error = function(e) {
      warning(paste(label, "failed:", e$message))
      NULL
    })
  }

  specs <- list(
    list(fn = calculate_flow_vols_by_year,  label = "Flow volumes"),
    list(fn = analyze_fdc_trends,           label = "FDC trends"),
    list(fn = analyze_flashiness_trends,    label = "Flashiness"),
    list(fn = analyze_flow_timing_trends,   label = "Flow timing"),
    list(fn = calculate_pulse_metrics,      label = "Pulse metrics"),
    list(fn = analyze_baseflow_indices,     label = "Baseflow indices"),
    list(fn = analyze_recession_parameters, label = "Recession parameters")
  )

  for (spec in specs) {
    out <- safe_call(spec$fn, spec$label)
    if (!is.null(out)) results <- c(results, out)
  }

  # Climate signatures
  if (has_climate && "PPT" %in% colnames(streamflow_data)) {
    climate_specs <- list(
      list(fn = analyze_Q_PPT_relationships,     label = "Q-PPT ratios"),
      list(fn = calculate_streamflow_elasticity,  label = "Elasticity"),
      list(fn = calculate_qp_seasonality,         label = "Q-P seasonality"),
      list(fn = calculate_average_storage,        label = "Average storage")
    )

    for (spec in climate_specs) {
      out <- safe_call(spec$fn, spec$label)
      if (!is.null(out)) results <- c(results, out)
    }
  }

  results
}

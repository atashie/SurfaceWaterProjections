#' Calculate Richards-Baker flashiness index trends
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list of 8 statistics for flashinessRB.
#' @export
analyze_flashiness_trends <- function(streamflow_data,
                                     trend_completeness = NULL,
                                     decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)

    # Metric columns use the CANONICAL names directly (matching Julia, which never
  # renames): the annual-values collector captures the column name it is given,
  # so a placeholder + post-hoc rename would put the wrong signature name in the
  # annual parquet (found 2026-08-25, Phase 3).
  flashiness_by_year <- data.frame(water_year = years, flashinessRB = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]

    if ("dowy" %in% colnames(year_data)) {
      year_data <- year_data[order(year_data$dowy), ]
    }

    q_values <- year_data$Q
    stopifnot(!any(is.na(q_values)))

    total_q <- sum(q_values)
    # <= 0 matches canonical Julia: zero or NEGATIVE total flow (negative-Q
    # days retained by default) is non-computable (2026-08-24 alignment fix)
    if (total_q <= 0) next

    rb_index <- sum(abs(diff(q_values))) / total_q
    flashiness_by_year$flashinessRB[flashiness_by_year$water_year == yr] <- rb_index
  }

  result <- generate_stats(flashiness_by_year, value_cols = "flashinessRB",
                           year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
  result
}

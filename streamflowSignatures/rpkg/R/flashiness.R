#' Calculate Richards-Baker flashiness index trends
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list of 8 statistics for flashinessRB.
#' @export
analyze_flashiness_trends <- function(streamflow_data,
                                     trend_completeness = NULL,
                                     decade_completeness = NULL) {
  required_cols <- c("water_year", "Q")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)

  flashiness_by_year <- data.frame(water_year = years, RB_index = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]

    if ("dowy" %in% colnames(year_data)) {
      year_data <- year_data[order(year_data$dowy), ]
    }

    q_values <- year_data$Q
    stopifnot(!any(is.na(q_values)))

    total_q <- sum(q_values)
    if (total_q == 0) next

    rb_index <- sum(abs(diff(q_values))) / total_q
    flashiness_by_year$RB_index[flashiness_by_year$water_year == yr] <- rb_index
  }

  result <- generate_stats(flashiness_by_year, value_cols = "RB_index",
                           year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness)
  names(result) <- gsub("RB_index", "flashinessRB", names(result))
  result
}

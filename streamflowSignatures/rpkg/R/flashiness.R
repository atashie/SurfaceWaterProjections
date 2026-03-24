#' Calculate Richards-Baker flashiness index trends
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list of 8 statistics for flashinessRB.
#' @export
analyze_flashiness_trends <- function(streamflow_data) {
  required_cols <- c("water_year", "Q")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)

  flashiness_by_year <- data.frame(water_year = years, RB_index = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]

    if (sum(!is.na(year_data$Q)) < pkg_env$flashiness_min_days) next

    if ("dowy" %in% colnames(year_data)) {
      year_data <- year_data[order(year_data$dowy), ]
    }

    q_values <- year_data$Q

    if (any(is.na(q_values))) {
      if (sum(is.na(q_values)) / length(q_values) > 0.2) next
      q_values <- approx(seq_along(q_values), q_values, seq_along(q_values), rule = 2)$y
    }

    total_q <- sum(q_values, na.rm = TRUE)
    if (total_q == 0) next

    rb_index <- sum(abs(diff(q_values)), na.rm = TRUE) / total_q
    flashiness_by_year$RB_index[flashiness_by_year$water_year == yr] <- rb_index
  }

  result <- generate_stats(flashiness_by_year, value_cols = "RB_index",
                           year_col = "water_year")
  names(result) <- gsub("RB_index", "flashinessRB", names(result))
  result
}

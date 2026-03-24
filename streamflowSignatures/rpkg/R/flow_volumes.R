#' Calculate flow volume signatures by water year
#'
#' Computes annual and seasonal flow totals plus 15 flow percentiles and
#' the Q95_Q10 ratio, then applies \code{generate_stats()} to produce 8
#' statistics per metric (22 metrics x 8 stats = 176 columns).
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, month, dowy.
#' @return Named list of statistics.
#' @export
calculate_flow_vols_by_year <- function(streamflow_data) {
  required_cols <- c("water_year", "Q", "month", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  min_days <- pkg_env$flow_volumes_min_days

  # Filter water years with insufficient data
  years <- unique(streamflow_data$water_year)
  valid_years <- vapply(years, function(yr) {
    sum(!is.na(streamflow_data$Q[streamflow_data$water_year == yr])) >= min_days
  }, logical(1))
  streamflow_data <- streamflow_data[streamflow_data$water_year %in% years[valid_years], ]

  if (nrow(streamflow_data) == 0) return(list())

  # Annual totals
  annual_totals <- aggregate(Q ~ water_year, data = streamflow_data, FUN = sum, na.rm = TRUE)
  names(annual_totals)[2] <- "Qann"

  # Seasonal totals
  season_months <- list(Qwin = c(12, 1, 2), Qspr = c(3, 4, 5),
                        Qsum = c(6, 7, 8), Qfal = c(9, 10, 11))
  all_metrics <- annual_totals
  for (nm in names(season_months)) {
    sub <- streamflow_data[streamflow_data$month %in% season_months[[nm]], ]
    if (nrow(sub) > 0) {
      seas <- aggregate(Q ~ water_year, data = sub, FUN = sum, na.rm = TRUE)
      names(seas)[2] <- nm
    } else {
      seas <- data.frame(water_year = numeric(), x = numeric())
      names(seas)[2] <- nm
    }
    all_metrics <- merge(all_metrics, seas, by = "water_year", all.x = TRUE)
  }

  # Percentiles
  pcts <- pkg_env$flow_volumes_percentiles
  for (p in pcts) {
    pct_df <- aggregate(Q ~ water_year, data = streamflow_data, FUN = function(x) {
      if (all(is.na(x))) return(NA_real_)
      quantile(x, probs = p / 100, na.rm = TRUE)
    })
    names(pct_df)[2] <- paste0("Q", p)
    all_metrics <- merge(all_metrics, pct_df, by = "water_year", all.x = TRUE)
  }

  # Q95_Q10 ratio
  if ("Q95" %in% names(all_metrics) && "Q10" %in% names(all_metrics)) {
    all_metrics[["Q95_Q10"]] <- all_metrics$Q95 - all_metrics$Q10
  } else {
    all_metrics[["Q95_Q10"]] <- NA_real_
  }

  metric_columns <- setdiff(names(all_metrics), "water_year")
  generate_stats(all_metrics, value_cols = metric_columns, year_col = "water_year")
}

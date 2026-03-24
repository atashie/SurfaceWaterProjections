#' Analyze flow timing trends
#'
#' Computes day of water year when cumulative flow reaches various percentiles,
#' the D25_to_D75 duration, and Dmax (day of peak flow).
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list of 8 statistics per timing metric (13 metrics).
#' @export
analyze_flow_timing_trends <- function(streamflow_data) {
  required_cols <- c("water_year", "Q", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)
  percentiles <- pkg_env$d_percentiles

  timing_by_year <- data.frame(water_year = years)
  for (p in percentiles) timing_by_year[[paste0("D", p, "_day")]] <- NA_real_
  timing_by_year$D25_to_D75 <- NA_real_
  timing_by_year$Dmax <- NA_real_

  min_days <- pkg_env$timing_min_days

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    if (nrow(year_data) < min_days) next

    year_data <- year_data[order(year_data$dowy), ]

    total_flow <- sum(year_data$Q, na.rm = TRUE)
    if (total_flow <= 0 || is.na(total_flow)) next

    q_for_cumsum <- year_data$Q
    q_for_cumsum[is.na(q_for_cumsum)] <- 0
    cum_pct <- (cumsum(q_for_cumsum) / total_flow) * 100

    idx <- which(timing_by_year$water_year == yr)

    for (p in percentiles) {
      above <- which(cum_pct >= p)
      if (length(above) > 0) {
        timing_by_year[idx, paste0("D", p, "_day")] <- year_data$dowy[above[1]]
      }
    }

    # D25_to_D75
    above_25 <- which(cum_pct >= 25)
    above_75 <- which(cum_pct >= 75)
    if (length(above_25) > 0 && length(above_75) > 0) {
      timing_by_year$D25_to_D75[idx] <- year_data$dowy[above_75[1]] - year_data$dowy[above_25[1]]
    }

    # Dmax
    max_idx <- which.max(year_data$Q)
    if (length(max_idx) > 0) {
      timing_by_year$Dmax[idx] <- year_data$dowy[max_idx]
    }
  }

  metric_columns <- c(paste0("D", percentiles, "_day"), "D25_to_D75", "Dmax")
  generate_stats(timing_by_year, value_cols = metric_columns, year_col = "water_year")
}

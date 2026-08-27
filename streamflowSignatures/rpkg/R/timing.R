#' Analyze flow timing trends
#'
#' Computes day of water year when cumulative flow reaches various percentiles,
#' the D25_to_D75 duration, and Dmax (day of peak flow).
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list of 8 statistics per timing metric (13 metrics).
#' @export
analyze_flow_timing_trends <- function(streamflow_data,
                                      trend_completeness = NULL,
                                      decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)
  percentiles <- pkg_env$d_percentiles

  # Zero-row input is legitimate and canonical Julia guards it. Without this,
  # data.frame(water_year = <empty>, X = NA_real_) raises "arguments imply
  # differing number of rows", which safe_call() swallows — silently dropping
  # the whole family for that gage (the class of the snow defect, 2026-08-26).
  # Route through generate_stats so the key set is identical to the normal path.
  if (nrow(streamflow_data) == 0L || length(years) == 0L) {
    .empty <- data.frame(water_year = integer(0))
    for (.m in c(paste0("D", percentiles, "_day"), "D25_to_D75", "Dmax")) .empty[[.m]] <- numeric(0)
    return(generate_stats(.empty, value_cols = c(paste0("D", percentiles, "_day"), "D25_to_D75", "Dmax"),
                          year_col = "water_year",
                          trend_completeness = trend_completeness,
                          decade_completeness = decade_completeness,
                          min_values_for_stats = min_values_for_stats,
                          changepoint = changepoint, collector = collector))
  }

  timing_by_year <- data.frame(water_year = years)
  for (p in percentiles) timing_by_year[[paste0("D", p, "_day")]] <- NA_real_
  timing_by_year$D25_to_D75 <- NA_real_
  timing_by_year$Dmax <- NA_real_

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]

    year_data <- year_data[order(year_data$dowy), ]

    total_flow <- sum(year_data$Q, na.rm = TRUE)
    if (total_flow <= 0 || is.na(total_flow)) next

    if (any(is.na(year_data$Q))) next
    cum_pct <- (cumsum(year_data$Q) / total_flow) * 100

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
  generate_stats(timing_by_year, value_cols = metric_columns, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
}

#' Calculate flow volume signatures by water year
#'
#' Computes annual and seasonal flow totals plus 15 flow percentiles and
#' the Q95_Q10 difference (Q95 - Q10), then applies \code{generate_stats()} to produce 8
#' statistics per metric (22 metrics x 8 stats = 176 columns).
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, month, dowy.
#' @param seasonal_flags Optional data.table from \code{preprocess_daily_data()}
#'   with columns water_year plus \code{<season>_complete} logicals. When
#'   supplied, seasonal totals for incomplete seasons are set to NA.
#' @return Named list of statistics.
#' @export
calculate_flow_vols_by_year <- function(streamflow_data, seasonal_flags = NULL,
                                       trend_completeness = NULL,
                                       decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q", "month", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

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

  # Apply seasonal_flags: set incomplete seasons to NA
  if (!is.null(seasonal_flags)) {
    sf <- data.table::as.data.table(seasonal_flags)
    season_col_map <- list(
      Qwin = "winter_complete",
      Qspr = "spring_complete",
      Qsum = "summer_complete",
      Qfal = "fall_complete"
    )
    for (metric_col in names(season_col_map)) {
      flag_col <- season_col_map[[metric_col]]
      if (metric_col %in% names(all_metrics) && flag_col %in% names(sf)) {
        for (i in seq_len(nrow(all_metrics))) {
          wy <- all_metrics$water_year[i]
          flag_row <- sf[sf$water_year == wy, ]
          if (nrow(flag_row) == 1 && !isTRUE(flag_row[[flag_col]])) {
            all_metrics[[metric_col]][i] <- NA_real_
          }
        }
      }
    }
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

  # Q95_Q10 difference
  if ("Q95" %in% names(all_metrics) && "Q10" %in% names(all_metrics)) {
    all_metrics[["Q95_Q10"]] <- all_metrics$Q95 - all_metrics$Q10
  } else {
    all_metrics[["Q95_Q10"]] <- NA_real_
  }

  metric_columns <- setdiff(names(all_metrics), "water_year")
  generate_stats(all_metrics, value_cols = metric_columns, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
}

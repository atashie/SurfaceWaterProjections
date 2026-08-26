#' Analyze Q-PPT relationships (runoff ratios)
#'
#' Computes annual and seasonal runoff ratios (Q/P) for each water year,
#' then applies generate_stats().
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, PPT, month.
#' @param seasonal_flags Optional data.table from \code{preprocess_daily_data()}
#'   with per-year seasonal completeness flags. When supplied, seasonal runoff
#'   ratios for incomplete seasons are set to NA.
#' @return Named list of 8 statistics per ratio (5 ratios = 40 values).
#' @export
analyze_Q_PPT_relationships <- function(streamflow_data, seasonal_flags = NULL,
                                       trend_completeness = NULL,
                                       decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q", "PPT", "month")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  min_annual_ppt   <- pkg_env$runoff_min_annual_ppt
  min_seasonal_ppt <- pkg_env$runoff_min_seasonal_ppt

  # Drop rows where either Q or PPT is NA (paired masking, matches canonical R na.action=na.omit)
  complete_data <- streamflow_data[!is.na(streamflow_data$Q) & !is.na(streamflow_data$PPT), ]

  # Annual totals
  annual <- aggregate(cbind(Q, PPT) ~ water_year, data = complete_data, FUN = sum)
  annual$annual_runoff_ratio <- ifelse(annual$PPT > min_annual_ppt, annual$Q / annual$PPT, NA_real_)

  # Seasonal
  season_defs <- list(
    winter = c(12, 1, 2), spring = c(3, 4, 5),
    summer = c(6, 7, 8), fall = c(9, 10, 11)
  )

  all_ratios <- data.frame(water_year = annual$water_year)
  all_ratios <- merge(all_ratios, annual[, c("water_year", "annual_runoff_ratio")],
                      by = "water_year", all.x = TRUE)

  for (nm in names(season_defs)) {
    sub <- complete_data[complete_data$month %in% season_defs[[nm]], ]
    if (nrow(sub) > 0) {
      seas <- aggregate(cbind(Q, PPT) ~ water_year, data = sub, FUN = sum)
      col_name <- paste0(nm, "_runoff_ratio")
      seas[[col_name]] <- ifelse(seas$PPT > min_seasonal_ppt, seas$Q / seas$PPT, NA_real_)
      all_ratios <- merge(all_ratios, seas[, c("water_year", col_name)],
                          by = "water_year", all.x = TRUE)
    }
  }

  # Apply seasonal_flags: set incomplete seasons to NA
  if (!is.null(seasonal_flags)) {
    sf <- data.table::as.data.table(seasonal_flags)
    season_col_map <- list(
      winter_runoff_ratio = "winter_complete",
      spring_runoff_ratio = "spring_complete",
      summer_runoff_ratio = "summer_complete",
      fall_runoff_ratio   = "fall_complete"
    )
    for (metric_col in names(season_col_map)) {
      flag_col <- season_col_map[[metric_col]]
      if (metric_col %in% names(all_ratios) && flag_col %in% names(sf)) {
        for (i in seq_len(nrow(all_ratios))) {
          wy <- all_ratios$water_year[i]
          flag_row <- sf[sf$water_year == wy, ]
          if (nrow(flag_row) == 1 && !isTRUE(flag_row[[flag_col]])) {
            all_ratios[[metric_col]][i] <- NA_real_
          }
        }
      }
    }
  }

  metric_cols <- setdiff(names(all_ratios), "water_year")
  result <- generate_stats(all_ratios, value_cols = metric_cols, year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)

  # Count years where annual runoff ratio > 2.0 (per-gage scalar diagnostic)
  n_high <- sum(all_ratios$annual_runoff_ratio > 2.0, na.rm = TRUE)
  result$runoff_ratio_high_count <- as.numeric(n_high)

  result
}

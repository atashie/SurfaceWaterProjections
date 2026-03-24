#' Analyze Q-PPT relationships (runoff ratios)
#'
#' Computes annual and seasonal runoff ratios (Q/P) for each water year,
#' then applies generate_stats().
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, PPT, month.
#' @return Named list of 8 statistics per ratio (5 ratios = 40 values).
#' @export
analyze_Q_PPT_relationships <- function(streamflow_data) {
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

  metric_cols <- setdiff(names(all_ratios), "water_year")
  generate_stats(all_ratios, value_cols = metric_cols, year_col = "water_year")
}

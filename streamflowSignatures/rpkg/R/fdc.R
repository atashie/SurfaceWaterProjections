#' Analyze flow duration curve trends
#'
#' Computes FDC slopes for the full range (FDCall), low-flow tail (FDC90th),
#' and mid-range (FDCmid) for each water year, then applies generate_stats().
#'
#' @param streamflow_data A data.frame with columns: water_year, Q.
#' @return Named list of 8 statistics per FDC metric (3 metrics = 24 values).
#' @export
analyze_fdc_trends <- function(streamflow_data) {
  required_cols <- c("water_year", "Q")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)
  if (length(years) == 0 || all(is.na(years))) {
    out <- c(empty_stats("FDCall"), empty_stats("FDC90th"), empty_stats("FDCmid"))
    return(out)
  }

  fdc_flow_floor <- pkg_env$fdc_flow_floor

  fdc_by_year <- data.frame(water_year = years, slp_all = NA_real_,
                            slp_90th = NA_real_, slp_mid = NA_real_)

  fdc_min_days <- pkg_env$fdc_min_days

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    if (nrow(year_data) < fdc_min_days) next

    Q_values <- year_data$Q[!is.na(year_data$Q)]
    Q_values <- Q_values[Q_values >= 0]
    if (length(Q_values) < fdc_min_days) next

    sorted_flows <- sort(Q_values, decreasing = TRUE)
    n <- length(sorted_flows)
    exceedance <- (1:n) / (n + 1)

    fdc <- data.frame(exceedance = exceedance,
                      log_flow = log10(sorted_flows + fdc_flow_floor))

    if (n < 10) next

    # Overall slope
    fit_all <- tryCatch(lm(log_flow ~ exceedance, data = fdc), error = function(e) NULL)
    if (!is.null(fit_all) && !is.na(coef(fit_all)[2]))
      fdc_by_year$slp_all[fdc_by_year$water_year == yr] <- coef(fit_all)[2]

    # 90th percentile and above (low flows)
    low_data <- fdc[fdc$exceedance >= 0.9, ]
    if (nrow(low_data) >= 3) {
      fit_low <- tryCatch(lm(log_flow ~ exceedance, data = low_data), error = function(e) NULL)
      if (!is.null(fit_low) && !is.na(coef(fit_low)[2]))
        fdc_by_year$slp_90th[fdc_by_year$water_year == yr] <- coef(fit_low)[2]
    }

    # Mid-range (20-80%)
    mid_data <- fdc[fdc$exceedance >= 0.2 & fdc$exceedance <= 0.8, ]
    if (nrow(mid_data) >= 3) {
      fit_mid <- tryCatch(lm(log_flow ~ exceedance, data = mid_data), error = function(e) NULL)
      if (!is.null(fit_mid) && !is.na(coef(fit_mid)[2]))
        fdc_by_year$slp_mid[fdc_by_year$water_year == yr] <- coef(fit_mid)[2]
    }
  }

  result <- generate_stats(fdc_by_year, value_cols = c("slp_all", "slp_90th", "slp_mid"),
                           year_col = "water_year")
  names(result) <- gsub("slp_all", "FDCall", names(result))
  names(result) <- gsub("slp_90th", "FDC90th", names(result))
  names(result) <- gsub("slp_mid", "FDCmid", names(result))
  result
}

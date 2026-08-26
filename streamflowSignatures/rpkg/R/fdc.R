#' Analyze flow duration curve trends
#'
#' Computes FDC slopes for the full range (FDCall), low-flow tail (FDC90th),
#' and mid-range (FDCmid) for each water year, then applies generate_stats().
#'
#' @param streamflow_data A data.frame with columns: water_year, Q.
#' @return Named list of 8 statistics per FDC metric (3 metrics = 24 values).
#' @export
analyze_fdc_trends <- function(streamflow_data,
                               trend_completeness = NULL,
                               decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  years <- unique(streamflow_data$water_year)
  if (length(years) == 0 || all(is.na(years))) {
    out <- c(empty_stats("FDCall"), empty_stats("FDC90th"), empty_stats("FDCmid"))
    return(out)
  }

  fdc_flow_floor <- pkg_env$fdc_flow_floor

  # Metric columns use the CANONICAL names directly (matching Julia, which never
  # renames): the annual-values collector captures the column name it is given,
  # so a placeholder + post-hoc rename would put the wrong signature name in the
  # annual parquet (found 2026-08-25, Phase 3).
  fdc_by_year <- data.frame(water_year = years, FDCall = NA_real_,
                            FDC90th = NA_real_, FDCmid = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]

    Q_values <- year_data$Q[!is.na(year_data$Q)]
    Q_values <- Q_values[Q_values >= 0]

    sorted_flows <- sort(Q_values, decreasing = TRUE)
    n <- length(sorted_flows)
    exceedance <- (1:n) / (n + 1)

    fdc <- data.frame(exceedance = exceedance,
                      log_flow = log10(sorted_flows + fdc_flow_floor))

    if (n < 10) next

    # Overall slope
    fit_all <- tryCatch(lm(log_flow ~ exceedance, data = fdc), error = function(e) NULL)
    if (!is.null(fit_all) && !is.na(coef(fit_all)[2]))
      fdc_by_year$FDCall[fdc_by_year$water_year == yr] <- coef(fit_all)[2]

    # 90th percentile and above (low flows)
    low_data <- fdc[fdc$exceedance >= 0.9, ]
    if (nrow(low_data) >= 3) {
      fit_low <- tryCatch(lm(log_flow ~ exceedance, data = low_data), error = function(e) NULL)
      if (!is.null(fit_low) && !is.na(coef(fit_low)[2]))
        fdc_by_year$FDC90th[fdc_by_year$water_year == yr] <- coef(fit_low)[2]
    }

    # Mid-range (20-80%)
    mid_data <- fdc[fdc$exceedance >= 0.2 & fdc$exceedance <= 0.8, ]
    if (nrow(mid_data) >= 3) {
      fit_mid <- tryCatch(lm(log_flow ~ exceedance, data = mid_data), error = function(e) NULL)
      if (!is.null(fit_mid) && !is.na(coef(fit_mid)[2]))
        fdc_by_year$FDCmid[fdc_by_year$water_year == yr] <- coef(fit_mid)[2]
    }
  }

  result <- generate_stats(fdc_by_year, value_cols = c("FDCall", "FDC90th", "FDCmid"),
                           year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
  result
}

#' Calculate Q-P seasonality metrics (Wrede et al. 2015)
#'
#' Computes qp_slope_sd (SD of monthly cumulative Q-P slopes) and
#' qp_bimodality (bimodality coefficient) for each water year.
#'
#' @param streamflow_data A data.table with columns: water_year, Q, PPT, month, dowy.
#' @param slope_window_days Window size for rolling slope (default from config).
#' @return Named list of 8 statistics per metric (2 x 8 = 16 values).
#' @export
#' @importFrom data.table as.data.table copy setorder data.table rbindlist
calculate_qp_seasonality <- function(streamflow_data, slope_window_days = NULL,
                                    trend_completeness = NULL,
                                    decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  if (is.null(slope_window_days)) slope_window_days <- pkg_env$qp_slope_window_days

  required_cols <- c("water_year", "Q", "PPT", "month", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  dt <- data.table::as.data.table(streamflow_data)
  years <- unique(dt$water_year)

  na_result <- c(empty_stats("qp_slope_sd"), empty_stats("qp_bimodality"))

  if (length(years) < pkg_env$qp_min_years) return(na_result)

  annual_list <- lapply(years, function(yr) {
    yd <- data.table::copy(dt[water_year == yr])

    if (any(is.na(yd$Q)) || any(is.na(yd$PPT))) return(NULL)
    data.table::setorder(yd, dowy)

    yd[, cum_Q := cumsum(Q)]
    yd[, cum_P := cumsum(PPT)]

    n <- nrow(yd)
    if (n < slope_window_days) return(NULL)

    slopes <- vapply(slope_window_days:n, function(end_idx) {
      si <- end_idx - slope_window_days + 1
      delta_P <- yd$cum_P[end_idx] - yd$cum_P[si]
      delta_Q <- yd$cum_Q[end_idx] - yd$cum_Q[si]
      if (abs(delta_P) < 0.01) NA_real_ else delta_Q / delta_P
    }, numeric(1))

    # Mid-point assignment matching Julia canonical
    mid_offsets <- slope_window_days:n - floor(slope_window_days / 2)
    slope_months <- yd$month[mid_offsets]

    slope_dt <- data.table::data.table(month = slope_months, slope = slopes)
    monthly <- slope_dt[!is.na(slope), .(mean_slope = mean(slope)), by = month]
    if (nrow(monthly) < 6) return(NULL)

    qp_slope_sd <- sd(monthly$mean_slope, na.rm = TRUE)

    slopes_clean <- slopes[!is.na(slopes)]
    if (length(slopes_clean) < 30) {
      qp_bimodality <- NA_real_
    } else {
      n_s <- length(slopes_clean)
      m <- mean(slopes_clean)
      s <- sd(slopes_clean)
      if (s < 1e-10) {
        qp_bimodality <- NA_real_
      } else {
        skewness <- sum((slopes_clean - m)^3) / (n_s * s^3)
        kurtosis <- sum((slopes_clean - m)^4) / (n_s * s^4)
        qp_bimodality <- if (kurtosis < 1e-10) NA_real_ else (skewness^2 + 1) / kurtosis
      }
    }

    data.table::data.table(water_year = yr, qp_slope_sd = qp_slope_sd,
                           qp_bimodality = qp_bimodality)
  })

  annual_metrics <- data.table::rbindlist(annual_list)
  if (is.null(annual_metrics) || nrow(annual_metrics) < 5) return(na_result)

  sd_stats <- generate_stats(annual_metrics, value_cols = "qp_slope_sd", year_col = "water_year",
                             trend_completeness = trend_completeness,
                             decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
  bi_df    <- annual_metrics[!is.na(qp_bimodality)]
  bi_stats <- if (nrow(bi_df) >= 3) {
    generate_stats(bi_df, value_cols = "qp_bimodality", year_col = "water_year",
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
  } else {
    empty_stats("qp_bimodality")
  }

  c(sd_stats, bi_stats)
}

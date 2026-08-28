#' Identify pulse events (runs above/below threshold)
#' @keywords internal
identify_pulses <- function(flow_vector, threshold, above = TRUE) {
  exceeds <- if (above) flow_vector > threshold else flow_vector < threshold
  exceeds[is.na(exceeds)] <- FALSE

  runs <- rle(exceeds)
  pulse_lengths <- runs$lengths[runs$values == TRUE]

  if (length(pulse_lengths) > 0) {
    list(n_pulses = length(pulse_lengths), mean_duration = mean(pulse_lengths))
  } else {
    list(n_pulses = 0L, mean_duration = NA_real_)
  }
}

#' Count flow reversals
#' @keywords internal
count_flow_reversals <- function(flow_vector, threshold_pct = NULL) {
  if (is.null(threshold_pct)) threshold_pct <- pkg_env$flow_reversal_threshold

  stopifnot(!any(is.na(flow_vector)))
  # No NA-cleaning happens here (the stopifnot above guarantees none);
  # the "_clean" name is historical.
  flow_clean <- flow_vector

  n <- length(flow_clean)
  if (n < 3) return(0L)

  reversal_count <- 0L
  for (i in 2:(n - 1)) {
    prev_change <- flow_clean[i] - flow_clean[i - 1]
    next_change <- flow_clean[i + 1] - flow_clean[i]
    threshold   <- abs(threshold_pct * flow_clean[i])

    if (abs(next_change) > threshold) {
      if (prev_change >= 0 && next_change < 0) reversal_count <- reversal_count + 1L
      else if (prev_change <= 0 && next_change > 0) reversal_count <- reversal_count + 1L
    }
  }
  reversal_count
}

#' Calculate pulse metrics
#'
#' Computes high/low pulse counts and durations (using year-specific and
#' period-of-record thresholds), TQmean, and annual/seasonal flow reversals.
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy, month.
#' @return Named list of 8 statistics per pulse metric (14 metrics).
#' @export
calculate_pulse_metrics <- function(streamflow_data,
                                   trend_completeness = NULL,
                                   decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q", "dowy", "month")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  q90_all <- quantile(streamflow_data$Q, probs = 0.90, na.rm = TRUE)
  q10_all <- quantile(streamflow_data$Q, probs = 0.10, na.rm = TRUE)

  years <- unique(streamflow_data$water_year)

  metric_names <- c("n_high_pulses_year", "n_low_pulses_year",
                    "n_high_pulses_all", "n_low_pulses_all",
                    "dur_high_pulses_year", "dur_low_pulses_year",
                    "dur_high_pulses_all", "dur_low_pulses_all",
                    "TQmean",
                    "Flow_Reversals_annual", "Flow_Reversals_winter",
                    "Flow_Reversals_spring", "Flow_Reversals_summer",
                    "Flow_Reversals_fall")

  pulse_metrics <- data.frame(water_year = years)
  for (nm in metric_names) pulse_metrics[[nm]] <- NA_real_

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]
    year_data <- year_data[order(year_data$dowy), ]

    q90_year <- quantile(year_data$Q, probs = 0.90, na.rm = TRUE)
    q10_year <- quantile(year_data$Q, probs = 0.10, na.rm = TRUE)
    if (is.na(q90_year) || is.na(q10_year)) next

    hp_yr  <- identify_pulses(year_data$Q, q90_year, above = TRUE)
    lp_yr  <- identify_pulses(year_data$Q, q10_year, above = FALSE)
    hp_all <- identify_pulses(year_data$Q, q90_all,  above = TRUE)
    lp_all <- identify_pulses(year_data$Q, q10_all,  above = FALSE)

    annual_mean <- mean(year_data$Q, na.rm = TRUE)
    days_above  <- sum(year_data$Q > annual_mean, na.rm = TRUE)
    total_valid <- sum(!is.na(year_data$Q))
    tqmean      <- (days_above / total_valid) * 100

    annual_rev <- count_flow_reversals(year_data$Q)

    season_months <- list(winter = c(12, 1, 2), spring = c(3, 4, 5),
                          summer = c(6, 7, 8), fall = c(9, 10, 11))
    seasonal_rev <- vapply(season_months, function(mos) {
      sub <- year_data[year_data$month %in% mos, ]
      if (nrow(sub) >= 30) count_flow_reversals(sub$Q) else NA_real_
    }, numeric(1))

    idx <- which(pulse_metrics$water_year == yr)
    pulse_metrics$n_high_pulses_year[idx]    <- hp_yr$n_pulses
    pulse_metrics$n_low_pulses_year[idx]     <- lp_yr$n_pulses
    pulse_metrics$n_high_pulses_all[idx]     <- hp_all$n_pulses
    pulse_metrics$n_low_pulses_all[idx]      <- lp_all$n_pulses
    pulse_metrics$dur_high_pulses_year[idx]  <- hp_yr$mean_duration
    pulse_metrics$dur_low_pulses_year[idx]   <- lp_yr$mean_duration
    pulse_metrics$dur_high_pulses_all[idx]   <- hp_all$mean_duration
    pulse_metrics$dur_low_pulses_all[idx]    <- lp_all$mean_duration
    pulse_metrics$TQmean[idx]                <- tqmean
    pulse_metrics$Flow_Reversals_annual[idx] <- annual_rev
    pulse_metrics$Flow_Reversals_winter[idx] <- seasonal_rev["winter"]
    pulse_metrics$Flow_Reversals_spring[idx] <- seasonal_rev["spring"]
    pulse_metrics$Flow_Reversals_summer[idx] <- seasonal_rev["summer"]
    pulse_metrics$Flow_Reversals_fall[idx]   <- seasonal_rev["fall"]
  }

  generate_stats(pulse_metrics, value_cols = metric_names, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
}

#' Count days with negative streamflow per water year
#'
#' Counts the number of days with Q < 0 in each water year, then produces
#' 8 trend statistics via \code{generate_stats()}.
#'
#' @param streamflow_data A data.frame with columns: water_year, Q.
#' @param trend_completeness Forwarded to generate_stats.
#' @param decade_completeness Forwarded to generate_stats.
#' @return Named list with 8 statistics for negative_ann.
#' @export
calculate_negative_days <- function(streamflow_data,
                                    trend_completeness = NULL,
                                    decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) {
    warning(paste("calculate_negative_days: Missing columns:", paste(missing, collapse = ", ")))
    return(empty_stats("negative_ann"))
  }

  annual <- aggregate(Q ~ water_year, data = streamflow_data,
                      FUN = function(q) sum(!is.na(q) & q < 0))
  names(annual)[2] <- "negative_ann"

  generate_stats(annual, value_cols = "negative_ann", year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
}

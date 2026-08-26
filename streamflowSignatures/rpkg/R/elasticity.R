#' Calculate streamflow elasticity
#'
#' Computes static elasticity, rolling-window elasticity trends (Sawicz et al.
#' 2011), and year-over-year elasticity trends using consecutive-year
#' differences.
#'
#' @param streamflow_data A data.table with columns: water_year, Q, PPT.
#' @param rolling_window Size of rolling window in years (default from config).
#' @param trend_completeness Forwarded to generate_stats.
#' @param decade_completeness Forwarded to generate_stats.
#' @return Named list with elasticity_static, elasticity_rolling_* (8 stats),
#'   elasticity_annual_* (8 stats), and diagnostic scalars.
#' @export
#' @importFrom data.table as.data.table setorder data.table
calculate_streamflow_elasticity <- function(streamflow_data,
                                            rolling_window = NULL,
                                            trend_completeness = NULL,
                                            decade_completeness = NULL, changepoint = NULL, collector = NULL) {
  if (is.null(rolling_window)) rolling_window <- pkg_env$elasticity_window_years

  required_cols <- c("water_year", "Q", "PPT")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  dt <- data.table::as.data.table(streamflow_data)

  annual <- dt[, .(Q_annual = sum(Q, na.rm = TRUE),
                   P_annual = sum(PPT, na.rm = TRUE)),
               by = water_year]

  min_ppt   <- pkg_env$elasticity_min_annual_ppt
  min_years <- pkg_env$elasticity_min_years

  # Track totals BEFORE filtering (diagnostics)
  elasticity_years_total <- as.numeric(nrow(annual))

  annual <- annual[P_annual > min_ppt]
  elasticity_years_low_ppt <- elasticity_years_total - as.numeric(nrow(annual))

  # Build NA result template
  na_result <- c(
    list(elasticity_static = NA_real_),
    empty_stats("elasticity_rolling"),
    empty_stats("elasticity_annual"),
    list(elasticity_years_total = elasticity_years_total,
         elasticity_years_low_ppt = elasticity_years_low_ppt)
  )

  if (nrow(annual) < min_years) return(na_result)

  Q_mean <- mean(annual$Q_annual, na.rm = TRUE)
  P_mean <- mean(annual$P_annual, na.rm = TRUE)

  if (P_mean <= 0 || Q_mean <= 0) return(na_result)

  annual[, dQ := Q_annual - Q_mean]
  annual[, dP := P_annual - P_mean]
  annual[, elasticity := ifelse(abs(dP) > 0.1, (dQ / dP) / (Q_mean / P_mean), NA_real_)]

  elasticity_static <- median(annual$elasticity, na.rm = TRUE)

  result <- list(
    elasticity_static = elasticity_static,
    elasticity_years_total = elasticity_years_total,
    elasticity_years_low_ppt = elasticity_years_low_ppt
  )

  # --- Rolling window elasticity (Sawicz et al. 2011) ---
  data.table::setorder(annual, water_year)
  n <- nrow(annual)

  if (n >= rolling_window + 3) {
    roll_years <- annual$water_year[rolling_window:n]
    roll_vals  <- vapply(rolling_window:n, function(end_idx) {
      start_idx <- end_idx - rolling_window + 1
      w <- annual[start_idx:end_idx]
      Qm <- mean(w$Q_annual, na.rm = TRUE)
      Pm <- mean(w$P_annual, na.rm = TRUE)
      w[, dQ_w := Q_annual - Qm]
      w[, dP_w := P_annual - Pm]
      w[, e_w := ifelse(abs(dP_w) > 0.1, (dQ_w / dP_w) / (Qm / Pm), NA_real_)]
      median(w$e_w, na.rm = TRUE)
    }, numeric(1))

    if (length(roll_vals) >= 3) {
      roll_dt <- data.table::data.table(water_year = roll_years,
                                        elasticity_rolling = roll_vals)
      ts <- generate_stats(roll_dt, value_cols = "elasticity_rolling",
                           year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness, changepoint = changepoint, collector = collector)
      result <- c(result, ts)
    } else {
      result <- c(result, empty_stats("elasticity_rolling"))
    }
  } else {
    result <- c(result, empty_stats("elasticity_rolling"))
  }

  # --- Year-over-year elasticity (consecutive-year differences) ---
  # E_t = ((Q_t - Q_{t-1}) / (P_t - P_{t-1})) / (Q_mean / P_mean)
  Q_valid <- annual$Q_annual
  P_valid <- annual$P_annual
  years_valid <- annual$water_year

  yoy_elasticity <- numeric(0)
  yoy_years <- numeric(0)
  for (i in 2:length(Q_valid)) {
    dQ <- Q_valid[i] - Q_valid[i - 1]
    dP <- P_valid[i] - P_valid[i - 1]
    if (abs(dP) > 0.1) {
      yoy_elasticity <- c(yoy_elasticity, (dQ / dP) / (Q_mean / P_mean))
      yoy_years <- c(yoy_years, years_valid[i])
    }
  }

  if (length(yoy_elasticity) >= 3) {
    yoy_dt <- data.table::data.table(water_year = yoy_years,
                                     elasticity_annual = yoy_elasticity)
    yoy_stats <- generate_stats(yoy_dt, value_cols = "elasticity_annual",
                                year_col = "water_year",
                                trend_completeness = trend_completeness,
                                decade_completeness = decade_completeness, changepoint = changepoint, collector = collector)
    result <- c(result, yoy_stats)
  } else {
    result <- c(result, empty_stats("elasticity_annual"))
  }

  result
}

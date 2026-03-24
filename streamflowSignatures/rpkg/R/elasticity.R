#' Calculate streamflow elasticity
#'
#' Computes static elasticity and rolling-window elasticity trend statistics.
#'
#' @param streamflow_data A data.table with columns: water_year, Q, PPT.
#' @param rolling_window Size of rolling window in years (default from config).
#' @return Named list with elasticity_static and 8 trend statistics.
#' @export
#' @importFrom data.table as.data.table setorder data.table
calculate_streamflow_elasticity <- function(streamflow_data,
                                            rolling_window = NULL) {
  if (is.null(rolling_window)) rolling_window <- pkg_env$elasticity_window_years

  required_cols <- c("water_year", "Q", "PPT")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  dt <- data.table::as.data.table(streamflow_data)

  annual <- dt[, .(Q_annual = sum(Q, na.rm = TRUE),
                   P_annual = sum(PPT, na.rm = TRUE),
                   n_valid  = sum(!is.na(Q) & !is.na(PPT))),
               by = water_year]

  annual[, expected_days := ifelse(
    ((water_year %% 4 == 0) & (water_year %% 100 != 0)) | (water_year %% 400 == 0),
    366L, 365L)]

  min_completeness <- pkg_env$elasticity_min_data_completeness
  min_ppt          <- pkg_env$elasticity_min_annual_ppt
  min_years        <- pkg_env$elasticity_min_years

  annual <- annual[n_valid >= min_completeness * expected_days]
  annual <- annual[P_annual > min_ppt]

  na_result <- list(
    elasticity_static = NA_real_,
    elasticity_senn_slp = NA_real_, elasticity_linear_slp = NA_real_,
    elasticity_spearman_rho = NA_real_, elasticity_spearman_pval = NA_real_,
    elasticity_mk_rho = NA_real_, elasticity_mk_pval = NA_real_,
    elasticity_mean = NA_real_, elasticity_median = NA_real_
  )

  if (nrow(annual) < min_years) return(na_result)

  Q_mean <- mean(annual$Q_annual, na.rm = TRUE)
  P_mean <- mean(annual$P_annual, na.rm = TRUE)

  annual[, dQ := Q_annual - Q_mean]
  annual[, dP := P_annual - P_mean]
  annual[, elasticity := ifelse(abs(dP) > 0.1, (dQ / dP) / (Q_mean / P_mean), NA_real_)]

  elasticity_static <- median(annual$elasticity, na.rm = TRUE)

  if (nrow(annual) < rolling_window) {
    na_result$elasticity_static <- elasticity_static
    return(na_result)
  }

  data.table::setorder(annual, water_year)
  n <- nrow(annual)

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

  roll_dt <- data.table::data.table(water_year = roll_years,
                                    elasticity_rolling = roll_vals)
  ts <- generate_stats(roll_dt, value_cols = "elasticity_rolling", year_col = "water_year")

  list(
    elasticity_static       = elasticity_static,
    elasticity_senn_slp     = ts$elasticity_rolling_senn_slp,
    elasticity_linear_slp   = ts$elasticity_rolling_linear_slp,
    elasticity_spearman_rho = ts$elasticity_rolling_spearman_rho,
    elasticity_spearman_pval = ts$elasticity_rolling_spearman_pval,
    elasticity_mk_rho       = ts$elasticity_rolling_mk_rho,
    elasticity_mk_pval      = ts$elasticity_rolling_mk_pval,
    elasticity_mean         = ts$elasticity_rolling_mean,
    elasticity_median       = ts$elasticity_rolling_median
  )
}

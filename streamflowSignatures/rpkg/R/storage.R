#' Calculate average storage metric (Peters & Aulenbach 2011)
#'
#' Uses simplified water balance dS = P - Q (no ET) to estimate mean annual
#' catchment storage in mm.
#'
#' @param streamflow_data A data.table with columns: water_year, Q, PPT, dowy.
#' @return Named list of 8 statistics for avg_storage.
#' @export
#' @importFrom data.table as.data.table copy setorder rbindlist
calculate_average_storage <- function(streamflow_data,
                                     trend_completeness = NULL,
                                     decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q", "PPT", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  dt <- data.table::as.data.table(streamflow_data)
  years <- unique(dt$water_year)
  na_result <- empty_stats("avg_storage")

  if (length(years) < pkg_env$storage_min_years) return(na_result)

  annual_list <- lapply(years, function(yr) {
    yd <- data.table::copy(dt[water_year == yr])

    if (any(is.na(yd$Q)) || any(is.na(yd$PPT))) return(NULL)
    data.table::setorder(yd, dowy)

    yd[, dS := PPT - Q]
    yd[, S := cumsum(dS)]

    mean_Q  <- mean(yd$Q, na.rm = TRUE)
    sq_data <- yd[!is.na(Q) & !is.na(S)]
    if (nrow(sq_data) < 10) return(NULL)

    # Julia applies a SECOND guard here: at least 10 DISTINCT Q values, not just
    # 10 valid rows (julia/src/storage.jl "if length(Q_unique) < 10"). rpkg had
    # only the row-count guard, so a year with e.g. 365 rows but 5 distinct flows
    # produced an avg_storage Julia declines to compute — 176 such gage-years,
    # caught 2026-08-26 by the cross-language annual gate.
    # (R's unique() uses == semantics, so -0.0 and +0.0 collapse; Julia's uses
    # isequal and keeps them distinct. That is the documented signed-zero
    # divergence, worth exactly 1 annual row in 18.9M — see CHANGELOG.)
    if (length(unique(sq_data$Q)) < 10) return(NULL)

    # The interpolation itself needs no de-duplication: R's approx() defaults to
    # ties = mean, which is precisely the behaviour Julia hand-rolls via S_means.

    data.table::setorder(sq_data, Q)
    avg_stor <- tryCatch(approx(sq_data$Q, sq_data$S, xout = mean_Q, rule = 2)$y,
                         error = function(e) NA_real_)

    data.table::data.table(water_year = yr, avg_storage = avg_stor)
  })

  annual_storage <- data.table::rbindlist(annual_list)
  if (is.null(annual_storage) || nrow(annual_storage) < 5) return(na_result)

  storage_df <- annual_storage[!is.na(avg_storage)]
  if (nrow(storage_df) < 3) return(na_result)

  generate_stats(storage_df, value_cols = "avg_storage", year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
}

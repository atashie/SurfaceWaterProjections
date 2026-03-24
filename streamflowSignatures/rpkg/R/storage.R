#' Calculate average storage metric (Peters & Aulenbach 2011)
#'
#' Uses simplified water balance dS = P - Q (no ET) to estimate mean annual
#' catchment storage in mm.
#'
#' @param streamflow_data A data.table with columns: water_year, Q, PPT, dowy.
#' @return Named list of 8 statistics for avg_storage.
#' @export
#' @importFrom data.table as.data.table copy setorder rbindlist
calculate_average_storage <- function(streamflow_data) {
  required_cols <- c("water_year", "Q", "PPT", "dowy")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  dt <- data.table::as.data.table(streamflow_data)
  years <- unique(dt$water_year)
  na_result <- empty_stats("avg_storage")

  if (length(years) < pkg_env$storage_min_years) return(na_result)

  annual_list <- lapply(years, function(yr) {
    yd <- data.table::copy(dt[water_year == yr])
    if (nrow(yd) < pkg_env$storage_min_days_per_year) return(NULL)

    na_frac_Q   <- sum(is.na(yd$Q)) / nrow(yd)
    na_frac_PPT <- sum(is.na(yd$PPT)) / nrow(yd)
    if (na_frac_Q > pkg_env$storage_max_na_frac || na_frac_PPT > pkg_env$storage_max_na_frac) return(NULL)

    yd[is.na(Q), Q := 0]
    yd[is.na(PPT), PPT := 0]
    data.table::setorder(yd, dowy)

    yd[, dS := PPT - Q]
    yd[, S := cumsum(dS)]

    mean_Q  <- mean(yd$Q, na.rm = TRUE)
    sq_data <- yd[!is.na(Q) & !is.na(S)]
    if (nrow(sq_data) < 10) return(NULL)

    data.table::setorder(sq_data, Q)
    avg_stor <- tryCatch(approx(sq_data$Q, sq_data$S, xout = mean_Q, rule = 2)$y,
                         error = function(e) NA_real_)

    data.table::data.table(water_year = yr, avg_storage = avg_stor)
  })

  annual_storage <- data.table::rbindlist(annual_list)
  if (is.null(annual_storage) || nrow(annual_storage) < 5) return(na_result)

  storage_df <- annual_storage[!is.na(avg_storage)]
  if (nrow(storage_df) < 3) return(na_result)

  generate_stats(storage_df, value_cols = "avg_storage", year_col = "water_year")
}

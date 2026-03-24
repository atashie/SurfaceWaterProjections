#' Read a parquet file into a data.table
#'
#' @param path Path to parquet file.
#' @return A data.table.
#' @export
read_parquet <- function(path) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is needed. Install with install.packages('arrow')")
  }
  dt <- data.table::as.data.table(arrow::read_parquet(path))

  # Normalize column names to match expected schema
  nmap <- c(Date = "date", site_id = "gage_id", prcp = "PPT",
            site_no = "gage_id", PRCP = "PPT")
  for (old in names(nmap)) {
    if (old %in% names(dt) && !nmap[[old]] %in% names(dt)) {
      data.table::setnames(dt, old, nmap[[old]])
    }
  }
  dt
}

#' Add water year columns to a data.table
#'
#' Adds \code{water_year}, \code{month}, and \code{dowy} (day of water year)
#' columns computed from the \code{date} column.
#'
#' @param dt A data.table with a \code{date} column (Date class).
#' @return The data.table, modified in place, with new columns added.
#' @export
#' @importFrom data.table :=
add_water_year_columns <- function(dt) {
  if (!"date" %in% names(dt)) stop("Column 'date' not found")
  dt <- data.table::as.data.table(dt)

  yr  <- as.integer(format(dt$date, "%Y"))
  mon <- as.integer(format(dt$date, "%m"))

  dt[, water_year := ifelse(mon >= 10L, yr + 1L, yr)]
  dt[, month := mon]

  wy_start <- as.Date(paste0(ifelse(mon >= 10L, yr, yr - 1L), "-10-01"))
  dt[, dowy := as.integer(date - wy_start + 1L)]

  dt
}

#' Filter to qualifying water years using three-stage per-year filtering
#'
#' Matches the canonical R filtering in process_signatures_from_parquet():
#' \enumerate{
#'   \item Minimum days with Q > threshold per water year
#'   \item Minimum fraction of non-NA days per water year (accounting for leap years)
#'   \item Minimum number of qualifying water years overall
#' }
#'
#' @param dt A data.table with \code{water_year} and \code{Q} columns.
#' @return A list with \code{qualifying_years} (integer vector) and
#'   \code{qualifies} (logical).
#' @export
filter_qualifying_years <- function(dt) {
  min_q_val   <- pkg_env$min_q_value

  min_days_q  <- pkg_env$min_days_above_threshold
  min_frac    <- pkg_env$min_frac_good_data
  min_years   <- pkg_env$min_num_years

  water_years_to_use <- integer(0)

  for (wy in unique(dt$water_year)) {
    wy_data <- dt[dt$water_year == wy, ]

    # Stage 1: minimum days above Q threshold
    n_above <- sum(wy_data$Q > min_q_val, na.rm = TRUE)
    if (n_above < min_days_q) next

    # Stage 2: minimum fraction of non-NA data (leap-year aware)
    is_leap <- ((wy %% 4 == 0) & (wy %% 100 != 0)) | (wy %% 400 == 0)
    expected_days <- ifelse(is_leap, 366L, 365L)
    min_required  <- floor(expected_days * min_frac)
    n_good <- sum(!is.na(wy_data$Q))
    if (n_good < min_required) next

    water_years_to_use <- c(water_years_to_use, wy)
  }

  list(
    qualifying_years = water_years_to_use,
    qualifies        = length(water_years_to_use) >= min_years
  )
}

#' Write signature results to CSV
#'
#' @param results A data.frame or data.table of signature results.
#' @param path Output CSV path.
#' @export
write_signatures <- function(results, path) {
  data.table::fwrite(data.table::as.data.table(results), path)
}

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
            site_no = "gage_id", PRCP = "PPT", swe = "SWE",
            tmin = "TMIN", tmax = "TMAX")  # Daymet temperature pair (2026-10-12)
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
#' LEGACY-path filter: mirrors the per-year filtering of the deprecated R shim's
#' process_signatures_from_parquet() (R/helperFunctions.R). The non-legacy
#' (default) path uses preprocess_daily_data() instead:
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

#' Linearly fill internal NA runs of length <= max_gap with Julia's exact expression
#'
#' \code{left + (right - left) * offset / span} evaluated in that order (offset and span
#' integers), so filled values are bit-identical with \code{julia/src/io.jl} (review
#' 2026-10-12: \code{stats::approx} divides first and rounds differently). Leading /
#' trailing runs and runs longer than max_gap stay NA; max_gap 0 fills nothing.
#' @keywords internal
#' @noRd
.fill_internal_gaps <- function(x, max_gap) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0L || !anyNA(x)) return(x)
  r <- rle(is.na(x))
  ends <- cumsum(r$lengths); starts <- c(1L, ends[-length(ends)] + 1L)
  for (ri in seq_along(r$lengths)) {
    if (!r$values[ri]) next
    gap_len <- r$lengths[ri]
    if (gap_len > max_gap) next
    s <- starts[ri]; e <- ends[ri]
    if (s > 1L && e < n && is.finite(x[s - 1L]) && is.finite(x[e + 1L])) {   # finite neighbours only
      left_val <- x[s - 1L]; right_val <- x[e + 1L]
      span <- (e + 1L) - (s - 1L)
      for (k in s:e) x[k] <- left_val + (right_val - left_val) * (k - (s - 1L)) / span
    }
  }
  x
}

#' Preprocess daily data with NA handling, interpolation, and diagnostics
#'
#' For each water year, normalizes to a complete daily grid (Oct 1 - Sep 30),
#' computes raw diagnostics (NA count, max NA run, seasonal completeness,
#' negative-Q flag, constant-SD flag), rejects years exceeding thresholds,
#' interpolates internal gaps, and checks post-interpolation residuals.
#' The same logic is applied to PPT if present.
#'
#' @param gage_flow A data.table with columns: Date (Date class), Q, water_year,
#'   month, dowy. Optionally PPT for climate data.
#' @param config A list mirroring the \code{na_handling} section of the JSON
#'   config. If NULL, defaults are built from the package-internal \code{pkg_env}.
#' @return A list with components:
#'   \describe{
#'     \item{data}{Cleaned data.table (only valid years, interpolated)}
#'     \item{valid_years}{Integer vector of accepted water years for Q}
#'     \item{valid_climate_years}{Integer vector of accepted water years for PPT}
#'     \item{rejected_years}{data.table with water_year and reason columns}
#'     \item{seasonal_flags}{data.table with per-year seasonal completeness}
#'     \item{diagnostics}{data.table with per-year NA diagnostics}
#'   }
#' @export
#' @importFrom data.table as.data.table copy setorder data.table rbindlist
preprocess_daily_data <- function(gage_flow, config = NULL) {
  # Build config from pkg_env defaults if not provided
  if (is.null(config)) {
    config <- list(
      interpolation = list(
        max_gap_days = pkg_env$na_max_gap_days,
        method = pkg_env$na_interpolation_method,
        internal_only = pkg_env$na_internal_only
      ),
      year_rejection = list(
        max_raw_na_per_year = pkg_env$na_max_raw_na_per_year,
        reject_negative_flow = pkg_env$na_reject_negative_flow,
        reject_residual_na = pkg_env$na_reject_residual_na
      ),
      constant_sd_flag = list(
        enabled = pkg_env$na_constant_sd_enabled,
        min_nonzero_days_per_month = pkg_env$na_constant_sd_min_days,
        max_unique_values = pkg_env$na_constant_sd_max_unique
      ),
      seasonal_completeness = list(
        min_fraction = pkg_env$na_seasonal_min_fraction,
        season_definitions = pkg_env$na_season_definitions   # configured, not hard-coded (review round 3)
      ),
      climate_na_policy = list(
        max_raw_na_per_year_ppt = pkg_env$na_max_raw_na_ppt,
        max_interpolation_gap_ppt = pkg_env$na_max_gap_ppt,
        reject_negative_ppt = pkg_env$na_reject_negative_ppt
      ),
      temperature_na_policy = list(
        max_raw_na_per_year_temp = pkg_env$na_max_raw_na_temp,
        max_interpolation_gap_temp = pkg_env$na_max_gap_temp
      ),
      snow_na_policy = list(
        max_raw_na_per_year_swe = pkg_env$na_max_raw_na_swe,
        max_interpolation_gap_swe = pkg_env$na_max_gap_swe,
        reject_negative_swe = pkg_env$na_reject_negative_swe
      )
    )
  }

  # Extract parameters from config — every limit a whole number >= 0, every flag a real
  # logical, also when supplied explicitly (review round 3; Julia types them as Int)
  max_gap_days       <- .cfg_int(config$interpolation$max_gap_days, "max_gap_days", 3L)
  max_raw_na         <- .cfg_int(config$year_rejection$max_raw_na_per_year, "max_raw_na_per_year", 30L)
  reject_negative    <- .cfg_bool(config$year_rejection$reject_negative_flow, "reject_negative_flow", TRUE)
  reject_residual    <- .cfg_bool(config$year_rejection$reject_residual_na, "reject_residual_na", TRUE)
  constant_sd_enabled   <- .cfg_bool(config$constant_sd_flag$enabled, "constant_sd_enabled", TRUE)
  constant_sd_min_days  <- .cfg_int(config$constant_sd_flag$min_nonzero_days_per_month, "min_nonzero_days_per_month", 15L)
  constant_sd_max_unique <- .cfg_int(config$constant_sd_flag$max_unique_values, "max_unique_values", 1L)
  seasonal_min_frac  <- .cfg_num(config$seasonal_completeness$min_fraction, "min_fraction", 0.80)
  if (max_gap_days < 0L || max_raw_na < 0L) stop("preprocess_daily_data: missing-day limits and interpolation gaps must be >= 0")
  season_defs <- if (!is.null(config$seasonal_completeness$season_definitions)) {
    sd_in <- config$seasonal_completeness$season_definitions
    list(winter = .cfg_months(sd_in$winter, "winter", c(12L, 1L, 2L)), spring = .cfg_months(sd_in$spring, "spring", c(3L, 4L, 5L)),
         summer = .cfg_months(sd_in$summer, "summer", c(6L, 7L, 8L)), fall = .cfg_months(sd_in$fall, "fall", c(9L, 10L, 11L)))
  } else if (!is.null(pkg_env$na_season_definitions)) {
    pkg_env$na_season_definitions
  } else {
    list(winter = c(12L, 1L, 2L), spring = c(3L, 4L, 5L),
         summer = c(6L, 7L, 8L), fall = c(9L, 10L, 11L))
  }

  # Climate NA policy
  has_ppt <- "PPT" %in% names(gage_flow)
  has_swe <- "SWE" %in% names(gage_flow)
  has_temp <- "TMIN" %in% names(gage_flow) && "TMAX" %in% names(gage_flow)
  max_raw_na_ppt     <- .cfg_int(config$climate_na_policy$max_raw_na_per_year_ppt, "max_raw_na_per_year_ppt", 30L)
  max_gap_ppt        <- .cfg_int(config$climate_na_policy$max_interpolation_gap_ppt, "max_interpolation_gap_ppt", 3L)
  reject_negative_ppt <- .cfg_bool(config$climate_na_policy$reject_negative_ppt, "reject_negative_ppt", TRUE)
  if (max_raw_na_ppt < 0L || max_gap_ppt < 0L) stop("climate_na_policy limits must be >= 0")
  # SWE policy: the supplied config first, the package defaults when absent (review round 4)
  max_raw_na_swe     <- .cfg_int(config$snow_na_policy$max_raw_na_per_year_swe, "max_raw_na_per_year_swe", pkg_env$na_max_raw_na_swe %||% 30L)
  max_gap_swe        <- .cfg_int(config$snow_na_policy$max_interpolation_gap_swe, "max_interpolation_gap_swe", pkg_env$na_max_gap_swe %||% 3L)
  reject_negative_swe <- .cfg_bool(config$snow_na_policy$reject_negative_swe, "reject_negative_swe", pkg_env$na_reject_negative_swe %||% TRUE)
  if (max_raw_na_swe < 0L || max_gap_swe < 0L) stop("snow_na_policy limits must be >= 0")
  # temperature policy: the supplied config first, the package defaults when absent —
  # validated like the package config (whole non-negative counts; review round 2)
  max_raw_na_temp    <- .cfg_int(config$temperature_na_policy$max_raw_na_per_year_temp, "max_raw_na_per_year_temp", pkg_env$na_max_raw_na_temp %||% 30L)
  max_gap_temp       <- .cfg_int(config$temperature_na_policy$max_interpolation_gap_temp, "max_interpolation_gap_temp", pkg_env$na_max_gap_temp %||% 3L)
  if (max_raw_na_temp < 0L || max_gap_temp < 0L) stop("temperature_na_policy limits must be >= 0")

  gage_flow <- data.table::copy(data.table::as.data.table(gage_flow))

  # Normalize Date column name
  if (!"Date" %in% names(gage_flow) && "date" %in% names(gage_flow)) {
    data.table::setnames(gage_flow, "date", "Date")
  }
  if (!"Date" %in% names(gage_flow)) stop("gage_flow must have a 'Date' or 'date' column")

  if (!inherits(gage_flow$Date, "Date")) {
    gage_flow[, Date := as.Date(Date)]
  }

  # Observation columns are numeric from here on (Julia's _real_vector contract: real
  # numbers incl. logicals are observations, anything else is missing). Converting the
  # whole frame once keeps the per-year rbind below from restoring a list or character
  # type that would break the downstream families (review round 5, 2026-10-13).
  for (col in intersect(c("Q", "PPT", "SWE", "TMIN", "TMAX"), names(gage_flow))) {
    if (!is.double(gage_flow[[col]])) data.table::set(gage_flow, j = col, value = .real_values(gage_flow[[col]]))
  }

  water_years <- sort(unique(gage_flow$water_year))

  # Initialize result containers
  rejected_list  <- list()
  diag_list      <- list()
  seasonal_list  <- list()
  valid_years    <- integer(0)
  valid_climate_years <- integer(0)
  valid_swe_years <- integer(0)
  valid_temp_years <- integer(0)

  for (wy in water_years) {
    # Step a: Normalize to complete daily grid
    wy_start <- as.Date(paste0(wy - 1L, "-10-01"))
    wy_end   <- as.Date(paste0(wy, "-09-30"))
    full_dates <- seq.Date(wy_start, wy_end, by = "day")

    wy_data <- gage_flow[water_year == wy]
    wy_data <- wy_data[!duplicated(Date)]
    data.table::setorder(wy_data, Date)

    grid <- data.table::data.table(Date = full_dates)
    wy_data <- merge(grid, wy_data, by = "Date", all.x = TRUE)
    wy_data[, water_year := wy]
    wy_data[, month := as.integer(format(Date, "%m"))]
    wy_data[, dowy := as.integer(Date - wy_start + 1L)]

    # Observations are real numbers; non-finite values (NaN, ±Inf), strings, complex
    # and other types are missing, for every series (review 2026-10-12; mirrors
    # julia/src/io.jl and _real_vector)
    wy_data$Q <- .real_values(wy_data$Q); wy_data$Q[!is.finite(wy_data$Q)] <- NA_real_
    if (has_ppt) { wy_data$PPT <- .real_values(wy_data$PPT); wy_data$PPT[!is.finite(wy_data$PPT)] <- NA_real_ }
    if (has_swe) { wy_data$SWE <- .real_values(wy_data$SWE); wy_data$SWE[!is.finite(wy_data$SWE)] <- NA_real_ }

    # Step b: Compute raw diagnostics (BEFORE interpolation)
    q_vals <- wy_data$Q
    raw_na_count <- sum(is.na(q_vals))

    # NA-cause tracking: USGS marks ice-affected days with an "Ice" qualifier in
    # the flag column. Mirrors julia/src/io.jl (regex [Ii]ce over the flags of
    # NA days); feeds the per-gage `ice_affected_days_total` diagnostic. Absent
    # a flag column this is simply 0.
    na_cause_ice <- 0L
    if (raw_na_count > 0 && "flag" %in% names(wy_data)) {
      na_flags <- wy_data$flag[is.na(q_vals)]
      na_cause_ice <- sum(!is.na(na_flags) & grepl("[Ii]ce", as.character(na_flags)))
    }

    raw_max_na_run <- 0L
    if (raw_na_count > 0) {
      rle_na <- rle(is.na(q_vals))
      na_runs <- rle_na$lengths[rle_na$values]
      if (length(na_runs) > 0) raw_max_na_run <- max(na_runs)
    }

    # Seasonal completeness from RAW observations
    season_flags <- list()
    for (sname in names(season_defs)) {
      s_months <- season_defs[[sname]]
      s_rows <- wy_data[month %in% s_months]
      n_total <- nrow(s_rows)
      n_valid <- if (n_total > 0) sum(!is.na(s_rows$Q)) else 0L
      frac <- if (n_total > 0) n_valid / n_total else 0
      season_flags[[paste0(sname, "_complete")]] <- frac >= seasonal_min_frac
      season_flags[[paste0(sname, "_frac")]] <- frac
    }
    season_flags$water_year <- wy

    # Negative Q check
    negative_flag <- any(q_vals < 0, na.rm = TRUE)

    # Constant-SD flag (sensor flatline detection)
    constant_sd_months <- 0L
    if (constant_sd_enabled) {
      for (m in unique(wy_data$month)) {
        m_vals <- wy_data[month == m & !is.na(Q) & Q > 0, Q]
        if (length(m_vals) >= constant_sd_min_days) {
          if (length(unique(m_vals)) <= constant_sd_max_unique) {
            constant_sd_months <- constant_sd_months + 1L
          }
        }
      }
    }

    # Step c: Year rejection (pre-interpolation)
    reject_reason <- NULL
    if (raw_na_count > max_raw_na) {
      reject_reason <- paste0("too many raw NAs (", raw_na_count, " > ", max_raw_na, ")")
    } else if (raw_max_na_run > max_gap_days) {
      reject_reason <- paste0("gap too large to interpolate (", raw_max_na_run, " > ", max_gap_days, " consecutive NAs)")
    } else if (reject_negative && negative_flag) {
      reject_reason <- "negative Q values"
    }

    if (!is.null(reject_reason)) {
      rejected_list[[length(rejected_list) + 1L]] <- data.table::data.table(
        water_year = wy, reason = reject_reason
      )
      diag_list[[length(diag_list) + 1L]] <- data.table::data.table(
        water_year = wy, raw_na_count = raw_na_count,
        raw_max_na_run = raw_max_na_run, interpolated_count = 0L,
        residual_na_count = raw_na_count, negative_flag = negative_flag,
        constant_sd_months = constant_sd_months, na_cause_ice = na_cause_ice
      )
      seasonal_list[[length(seasonal_list) + 1L]] <- data.table::as.data.table(season_flags)
      next
    }

    # Step d: Interpolate internal gaps <= max_gap_days
    interpolated_count <- 0L
    if (raw_na_count > 0) {
      n <- nrow(wy_data)
      na_mask <- is.na(wy_data$Q)

      rle_result <- rle(na_mask)
      run_ends <- cumsum(rle_result$lengths)
      run_starts <- c(1L, run_ends[-length(run_ends)] + 1L)

      for (ri in seq_along(rle_result$lengths)) {
        if (!rle_result$values[ri]) next
        gap_len <- rle_result$lengths[ri]
        if (gap_len > max_gap_days) next

        s <- run_starts[ri]
        e <- run_ends[ri]

        left_ok  <- (s > 1L) && !is.na(wy_data$Q[s - 1L])
        right_ok <- (e < n)  && !is.na(wy_data$Q[e + 1L])

        if (left_ok && right_ok) {
          # Julia's exact expression order (bit-identical fills across languages)
          left_val  <- wy_data$Q[s - 1L]
          right_val <- wy_data$Q[e + 1L]
          span <- (e + 1L) - (s - 1L)
          for (k in s:e) wy_data$Q[k] <- left_val + (right_val - left_val) * (k - (s - 1L)) / span
          interpolated_count <- interpolated_count + gap_len
        }
      }
    }

    # Step e: Post-interpolation residual check
    residual_na_count <- sum(is.na(wy_data$Q))
    if (reject_residual && residual_na_count > 0) {
      rejected_list[[length(rejected_list) + 1L]] <- data.table::data.table(
        water_year = wy, reason = paste0("residual boundary NAs (", residual_na_count, " remaining after interpolation)")
      )
      diag_list[[length(diag_list) + 1L]] <- data.table::data.table(
        water_year = wy, raw_na_count = raw_na_count,
        raw_max_na_run = raw_max_na_run,
        interpolated_count = interpolated_count,
        residual_na_count = residual_na_count,
        negative_flag = negative_flag,
        constant_sd_months = constant_sd_months, na_cause_ice = na_cause_ice
      )
      seasonal_list[[length(seasonal_list) + 1L]] <- data.table::as.data.table(season_flags)
      next
    }

    # Step f: Apply same logic to PPT if present
    ppt_valid <- TRUE
    if (has_ppt) {
      ppt_vals <- wy_data$PPT
      ppt_raw_na <- sum(is.na(ppt_vals))

      if (ppt_raw_na > max_raw_na_ppt) {
        ppt_valid <- FALSE
      } else if (ppt_raw_na > 0) {
        rle_ppt <- rle(is.na(ppt_vals))
        ppt_max_run <- max(rle_ppt$lengths[rle_ppt$values], 0L)
        if (ppt_max_run > max_gap_ppt) {
          ppt_valid <- FALSE
        } else {
          rle_result_ppt <- rle(is.na(wy_data$PPT))
          run_ends_ppt <- cumsum(rle_result_ppt$lengths)
          run_starts_ppt <- c(1L, run_ends_ppt[-length(run_ends_ppt)] + 1L)
          n_ppt <- nrow(wy_data)

          for (ri in seq_along(rle_result_ppt$lengths)) {
            if (!rle_result_ppt$values[ri]) next
            gap_len <- rle_result_ppt$lengths[ri]
            if (gap_len > max_gap_ppt) next

            s <- run_starts_ppt[ri]
            e <- run_ends_ppt[ri]
            left_ok  <- (s > 1L) && !is.na(wy_data$PPT[s - 1L])
            right_ok <- (e < n_ppt) && !is.na(wy_data$PPT[e + 1L])

            if (left_ok && right_ok) {
              left_val <- wy_data$PPT[s - 1L]; right_val <- wy_data$PPT[e + 1L]
              span <- (e + 1L) - (s - 1L)
              for (k in s:e) wy_data$PPT[k] <- left_val + (right_val - left_val) * (k - (s - 1L)) / span
            }
          }

          if (sum(is.na(wy_data$PPT)) > 0) ppt_valid <- FALSE
        }
      }

      if (ppt_valid && reject_negative_ppt && any(wy_data$PPT < 0, na.rm = TRUE)) {
        ppt_valid <- FALSE
      }
    }

    # Step f2: SWE handling (same rules as PPT, tracked separately)
    swe_valid <- TRUE
    if (has_swe) {
      swe_vals <- wy_data$SWE
      swe_raw_na <- sum(is.na(swe_vals))

      if (swe_raw_na > max_raw_na_swe) {
        swe_valid <- FALSE
      } else if (swe_raw_na > 0) {
        rle_swe <- rle(is.na(swe_vals))
        swe_max_run <- max(rle_swe$lengths[rle_swe$values], 0L)
        if (swe_max_run > max_gap_swe) {
          swe_valid <- FALSE
        } else {
          rle_res_swe <- rle(is.na(wy_data$SWE))
          run_ends_swe <- cumsum(rle_res_swe$lengths)
          run_starts_swe <- c(1L, run_ends_swe[-length(run_ends_swe)] + 1L)
          n_swe <- nrow(wy_data)
          for (ri in seq_along(rle_res_swe$lengths)) {
            if (!rle_res_swe$values[ri]) next
            if (rle_res_swe$lengths[ri] > max_gap_swe) next
            s2 <- run_starts_swe[ri]; e2 <- run_ends_swe[ri]
            left_ok  <- (s2 > 1L) && !is.na(wy_data$SWE[s2 - 1L])
            right_ok <- (e2 < n_swe) && !is.na(wy_data$SWE[e2 + 1L])
            if (left_ok && right_ok) {
              left_val <- wy_data$SWE[s2 - 1L]; right_val <- wy_data$SWE[e2 + 1L]
              span <- (e2 + 1L) - (s2 - 1L)
              for (k in s2:e2) wy_data$SWE[k] <- left_val + (right_val - left_val) * (k - (s2 - 1L)) / span
            }
          }
          if (sum(is.na(wy_data$SWE)) > 0) swe_valid <- FALSE
        }
      }
      if (swe_valid && isTRUE(reject_negative_swe) &&
          any(wy_data$SWE < 0, na.rm = TRUE)) {
        swe_valid <- FALSE
      }
    }

    # Step f3: TMIN/TMAX handling (2026-10-12): the SWE rules without the sign check
    # (temperatures are signed); non-finite values count as missing; the year is
    # temperature-valid only when BOTH series pass, and both are written back together
    temp_valid <- TRUE
    if (has_temp) {
      filled <- list()
      for (col in c("TMIN", "TMAX")) {
        tv <- .real_values(wy_data[[col]])
        tv[!is.finite(tv)] <- NA_real_
        if (sum(is.na(tv)) > max_raw_na_temp) {
          temp_valid <- FALSE
        } else if (sum(is.na(tv)) > 0) {
          tv <- .fill_internal_gaps(tv, max_gap_temp)
          if (sum(is.na(tv)) > 0) temp_valid <- FALSE
        }
        filled[[col]] <- tv
      }
      if (temp_valid) {
        wy_data$TMIN <- filled$TMIN
        wy_data$TMAX <- filled$TMAX
      }
    }

    # Step g: Build return objects for this year
    valid_years <- c(valid_years, wy)
    if (has_ppt && ppt_valid) {
      valid_climate_years <- c(valid_climate_years, wy)
    }
    if (has_swe && swe_valid) {
      valid_swe_years <- c(valid_swe_years, wy)
    }
    if (has_temp && temp_valid) {
      valid_temp_years <- c(valid_temp_years, wy)
    }

    diag_list[[length(diag_list) + 1L]] <- data.table::data.table(
      water_year = wy, raw_na_count = raw_na_count,
      raw_max_na_run = raw_max_na_run,
      interpolated_count = interpolated_count,
      residual_na_count = residual_na_count,
      negative_flag = negative_flag,
      constant_sd_months = constant_sd_months, na_cause_ice = na_cause_ice
    )
    seasonal_list[[length(seasonal_list) + 1L]] <- data.table::as.data.table(season_flags)

    # Update gage_flow with interpolated data for this year
    gage_flow <- gage_flow[water_year != wy]
    gage_flow <- rbind(gage_flow, wy_data, fill = TRUE)
  }

  # Sort the final data by Date
  data.table::setorder(gage_flow, Date)

  # Filter to valid years only
  cleaned_data <- gage_flow[water_year %in% valid_years]

  # Assemble return objects
  rejected_years <- if (length(rejected_list) > 0) {
    data.table::rbindlist(rejected_list)
  } else {
    data.table::data.table(water_year = integer(0), reason = character(0))
  }
  diagnostics <- if (length(diag_list) > 0) {
    data.table::rbindlist(diag_list)
  } else {
    data.table::data.table(
      water_year = integer(0), raw_na_count = integer(0),
      raw_max_na_run = integer(0), interpolated_count = integer(0),
      residual_na_count = integer(0), negative_flag = logical(0),
      constant_sd_months = integer(0), na_cause_ice = integer(0)
    )
  }
  seasonal_flags_dt <- if (length(seasonal_list) > 0) {
    data.table::rbindlist(seasonal_list, fill = TRUE)
  } else {
    data.table::data.table(
      water_year = integer(0), winter_complete = logical(0),
      spring_complete = logical(0), summer_complete = logical(0),
      fall_complete = logical(0)
    )
  }

  list(
    data = cleaned_data,
    valid_years = valid_years,
    valid_climate_years = valid_climate_years,
    valid_swe_years = valid_swe_years,
    valid_temp_years = valid_temp_years,
    rejected_years = rejected_years,
    seasonal_flags = seasonal_flags_dt,
    diagnostics = diagnostics
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

#' Climate metrics (guidelines 3.9, synced 2026-10-12; decisions: user, 2026-10-12)
#'
#' Port of julia/src/climate.jl (canonical). Thirty-nine per-water-year metrics from
#' basin-mean Daymet daily precipitation (prcp -> PPT, mm/day) and air temperature
#' (tmin -> TMIN, tmax -> TMAX, deg C), in two function families:
#' \code{calculate_precipitation_metrics} (20 bases on PPT-valid years; pure climate,
#' so it runs regardless of area_normalized) and \code{calculate_temperature_metrics}
#' (19 bases on an explicit temp_data frame of temperature-valid years, the snow
#' family's opt-in precedent).
#'
#' Conventions (docs/plans/2026-10-12-climate-metrics-plan.md): seasons DJF/MAM/JJA/SON
#' applied WITHIN the water year (fall = Oct + Nov of the water-year start plus Sep of
#' its end; suffixes win/spr/sum/fal); dowy 1 = Oct 1; ties take the first occurrence;
#' a year is used only on a complete grid (the water year's exact day count, dowy 1..n,
#' no NA); no Q-derived seasonal mask; a dry day has P <= dry_day_threshold_mm; wet-day
#' counts are strictly above each threshold; the multi-day maximum is a window inside
#' the water year; GDD = sum(max(0, (TMAX + TMIN)/2 - gdd_base_c)); frost is TMIN <
#' frost_temp_c. Sums use plain sequential double accumulation (\code{seq_sum}) so a
#' cumulative sum's last element equals the total it is compared against and the annual
#' values match Julia and Python bit for bit.
#' @name climate_metrics
NULL

.climate_season_keys <- c("win", "spr", "sum", "fal")
# Season months from the configured definitions (na_handling.seasonal_completeness
# .season_definitions, read in .onLoad), DJF/MAM/JJA/SON by default — as Julia's
# CFG_NA_SEASON_MONTHS (review 2026-10-12)
.climate_season_defs <- function() {
  sd <- pkg_env$na_season_definitions
  if (is.null(sd)) sd <- list(winter = c(12L, 1L, 2L), spring = c(3L, 4L, 5L), summer = c(6L, 7L, 8L), fall = c(9L, 10L, 11L))
  list(win = as.integer(sd$winter), spr = as.integer(sd$spring), sum = as.integer(sd$summer), fal = as.integer(sd$fall))
}

#' Plain left-to-right double accumulation (bit-identical with Julia/Python loops)
#' @keywords internal
seq_sum <- function(x) {
  s <- 0
  for (v in x) s <- s + v
  s
}

#' Sequential cumulative sum in double precision (R's cumsum accumulates in long double)
#' @keywords internal
seq_cumsum <- function(x) {
  out <- numeric(length(x))
  s <- 0
  for (i in seq_along(x)) { s <- s + x[i]; out[i] <- s }
  out
}

# Column-name fragment for a wet-day threshold, identical in the three languages: an
# integral value prints as an integer ("10"), anything else with C's %g ("0.5",
# "1e-05"), then "." -> "p" and "-" -> "m" ("0p5", "1em05")
.wet_threshold_label <- function(thr) {
  thr <- as.numeric(thr)
  s <- if (thr == floor(thr) && abs(thr) < 1e15) sprintf("%.0f", thr) else sprintf("%g", thr)
  gsub("-", "m", gsub(".", "p", s, fixed = TRUE), fixed = TRUE)
}

#' Precipitation metric names (20), derived from the config like Julia's registry
#' @export
precip_metrics <- function() {
  w <- as.integer(pkg_env$climate_multiday_window %||% 5L)
  frac <- as.numeric(pkg_env$climate_cum_fraction %||% 0.5)
  thrs <- as.numeric(pkg_env$climate_wet_thresholds_mm %||% c(10, 25, 50))
  c("prcp_annual", paste0("prcp_", .climate_season_keys),
    "prcp_1d_max", paste0("prcp_", w, "d_max"), "prcp_max_dowy",
    paste0("prcp_", as.integer(round(100 * frac)), "_dowy"), "prcp_cm_dowy",
    paste0("prcp_ndays_abv_", vapply(thrs, .wet_threshold_label, character(1)), "mm"),
    "prcp_annual_sd", "prcp_annual_cv",
    "prcp_intermit_annual", "prcp_intermit_annual_dowy_start", "prcp_intermit_annual_dowy_end",
    "prcp_gini", "prcp_seasonality_index")
}

#' Temperature metric names (19)
#' @export
temp_metrics <- function() {
  c("temp_max_annual", "temp_max_annual_dowy", "temp_min_annual", "temp_min_annual_dowy",
    paste0("temp_max_", .climate_season_keys), paste0("temp_min_", .climate_season_keys),
    "frost_first_dowy", "frost_last_dowy",
    "gdd_cum_annual", paste0("gdd_cum_", .climate_season_keys))
}

.temp_record_gated_metrics <- c("frost_first_dowy", "frost_last_dowy")

# Relative tolerance of the cumulative-fraction crossing (mirrors julia/src/climate.jl
# CLIMATE_CROSSING_REL_TOL): a mathematical tie must not be decided by the last bit
.climate_crossing_rel_tol <- 1e-12

.water_year_days <- function(yr) {
  as.integer(as.Date(paste0(yr, "-09-30")) - as.Date(paste0(yr - 1L, "-10-01"))) + 1L
}

# TRUE when the rows are EXACTLY dowy 1..n (compared as the numbers given — no
# truncation) and n is the water year's day count
.complete_year_grid <- function(dowy, n, yr) {
  if (n != .water_year_days(yr)) return(FALSE)
  if (any(!is.finite(dowy))) return(FALSE)
  isTRUE(all(dowy == seq_len(n)))
}

# Calendar month of each day of the complete grid, from the date (not a month column)
.grid_months <- function(yr, n) {
  as.integer(format(as.Date(paste0(yr - 1L, "-10-01")) + (seq_len(n) - 1L), "%m"))
}

.empty_family_stats <- function(metrics, trend_completeness, decade_completeness,
                                min_values_for_stats, changepoint, collector) {
  empty <- data.frame(water_year = integer(0))
  for (m in metrics) empty[[m]] <- numeric(0)
  generate_stats(empty, value_cols = metrics, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness,
                 min_values_for_stats = min_values_for_stats,
                 changepoint = changepoint, collector = collector)
}

# Record-anchored decade gate (the snow family's rule): gated metrics need
# >= decade_completeness of the anchor years in BOTH the first and last decade of
# the record to be non-NA, else their six trend statistics are forced to NA.
.record_anchored_force_skip <- function(annual, anchor_col, gated, decade_completeness) {
  anchor_mask <- !is.na(annual[[anchor_col]])
  anchor_years <- as.integer(annual$water_year[anchor_mask])
  if (length(anchor_years) == 0) return(NULL)
  rec_min <- min(anchor_years); rec_max <- max(anchor_years)
  if (rec_max - rec_min + 1 < 10) return(NULL)
  windows <- list(c(rec_min, rec_min + 9), c(rec_max - 9, rec_max))
  all_years <- as.integer(annual$water_year)
  skip <- character(0)
  for (mn in gated) {
    vals <- annual[[mn]]
    for (w in windows) {
      in_win <- anchor_mask & all_years >= w[1] & all_years <= w[2]
      den <- sum(in_win)
      if (den == 0) next
      num <- sum(!is.na(vals[in_win]))
      if ((num / den) < decade_completeness) { skip <- c(skip, mn); break }
    }
  }
  if (length(skip) > 0) skip else NULL
}

#' Calculate 20 precipitation metrics with trend statistics from daily basin-mean PPT
#'
#' Callers pass the PPT-qualified frame (\code{climate_data}, filtered to
#' \code{valid_climate_years}). A year whose grid is incomplete or holds an NA is left
#' NA in every metric. See julia/src/climate.jl for the definitions.
#'
#' @param df data.frame with water_year, dowy, PPT.
#' @param trend_completeness,decade_completeness,min_values_for_stats,changepoint,collector
#'   Forwarded to \code{generate_stats}.
#' @return Named list: 20 metrics x 8 statistics (+ 8 changepoint fields each).
#' @export
calculate_precipitation_metrics <- function(df, trend_completeness = NULL,
                                            decade_completeness = NULL,
                                            min_values_for_stats = NULL,
                                            changepoint = NULL, collector = NULL) {
  df <- as.data.frame(df)
  metrics <- precip_metrics()
  if (!all(c("PPT", "water_year", "dowy") %in% names(df))) {
    return(.empty_family_stats(metrics, trend_completeness, decade_completeness,
                               min_values_for_stats, changepoint, collector))
  }
  dry_thr <- as.numeric(pkg_env$climate_dry_day_mm %||% 0.1)
  wet_thrs <- as.numeric(pkg_env$climate_wet_thresholds_mm %||% c(10, 25, 50))
  w <- as.integer(pkg_env$climate_multiday_window %||% 5L)
  frac <- as.numeric(pkg_env$climate_cum_fraction %||% 0.5)
  col_multi <- paste0("prcp_", w, "d_max")
  col_cum <- paste0("prcp_", as.integer(round(100 * frac)), "_dowy")
  wet_cols <- paste0("prcp_ndays_abv_", vapply(wet_thrs, .wet_threshold_label, character(1)), "mm")

  years <- unique(df$water_year)
  annual <- data.frame(water_year = years)
  for (m in metrics) annual[[m]] <- rep(NA_real_, nrow(annual))  # 0-row safe (schema contract)

  for (k in seq_along(years)) {
    yr <- as.integer(years[k])
    yd <- df[df$water_year == yr, , drop = FALSE]
    if (nrow(yd) == 0) next
    yd <- yd[order(yd$dowy), , drop = FALSE]
    P <- .real_values(yd$PPT) + 0   # real observations only; + 0 turns a signed -0 into +0
    dowy <- as.numeric(yd$dowy)
    n <- length(P)
    if (!.complete_year_grid(dowy, n, yr) || any(!is.finite(P))) next
    months <- .grid_months(yr, n)
    season_defs <- .climate_season_defs()

    total <- seq_sum(P)
    annual$prcp_annual[k] <- total
    for (s in .climate_season_keys) {
      smask <- months %in% season_defs[[s]]
      annual[[paste0("prcp_", s)]][k] <- seq_sum(P[smask])
    }

    imax <- which.max(P)                      # first occurrence on ties
    annual$prcp_1d_max[k] <- P[imax]
    annual$prcp_max_dowy[k] <- dowy[imax]

    if (n >= w) {
      best <- -Inf
      for (i in seq_len(n - w + 1L)) {
        s_ <- seq_sum(P[i:(i + w - 1L)])
        if (s_ > best) best <- s_
      }
      annual[[col_multi]][k] <- best
    }

    if (total > 0) {
      cum <- seq_cumsum(P)
      thr <- frac * total - .climate_crossing_rel_tol * total
      reached <- which(cum >= thr)
      if (length(reached) > 0) annual[[col_cum]][k] <- dowy[reached[1]]
      annual$prcp_cm_dowy[k] <- seq_sum(dowy * P) / total
      y <- sort(P)
      annual$prcp_gini[k] <- seq_sum((2 * seq_len(n) - n - 1) * y) / (n * total)
      X <- vapply(1:12, function(m) seq_sum(P[months == m]), numeric(1))
      R <- seq_sum(X)
      if (R > 0) annual$prcp_seasonality_index[k] <- seq_sum(abs(X - R / 12)) / R
    }

    # counts of logicals are exact integers (<= 366): the sequential-sum contract covers
    # floating-point accumulation only (review round 2 note 6, exempted)
    for (j in seq_along(wet_thrs)) annual[[wet_cols[j]]][k] <- as.numeric(sum(P > wet_thrs[j]))

    mu <- total / n
    sd_ <- sqrt(seq_sum((P - mu)^2) / (n - 1))
    annual$prcp_annual_sd[k] <- sd_
    annual$prcp_annual_cv[k] <- if (mu > 0) sd_ / mu else NA_real_

    dry <- P <= dry_thr
    best_len <- 0L; best_start <- 0L
    i <- 1L
    while (i <= n) {
      if (dry[i]) {
        j <- i
        while (j < n && dry[j + 1L]) j <- j + 1L
        len <- j - i + 1L
        if (len > best_len) { best_len <- len; best_start <- i }
        i <- j + 1L
      } else {
        i <- i + 1L
      }
    }
    annual$prcp_intermit_annual[k] <- as.numeric(best_len)
    if (best_len > 0) {
      annual$prcp_intermit_annual_dowy_start[k] <- dowy[best_start]
      annual$prcp_intermit_annual_dowy_end[k] <- dowy[best_start + best_len - 1L]
    }
  }

  generate_stats(annual, value_cols = metrics, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness,
                 min_values_for_stats = min_values_for_stats,
                 changepoint = changepoint, collector = collector)
}

#' Calculate 19 temperature metrics with trend statistics from daily TMIN/TMAX
#'
#' Callers pass the temperature-valid frame (\code{temp_data}); in
#' \code{calculate_all_signatures} this runs ONLY when such a frame is provided.
#' The record-anchored decade gate applies to the two frost metrics when
#' \code{climate.record_anchored_decade_gate} is set.
#'
#' @param df data.frame with water_year, dowy, TMIN, TMAX.
#' @param trend_completeness,decade_completeness,min_values_for_stats,changepoint,collector
#'   Forwarded to \code{generate_stats}.
#' @return Named list: 19 metrics x 8 statistics (+ 8 changepoint fields each).
#' @export
calculate_temperature_metrics <- function(df, trend_completeness = NULL,
                                          decade_completeness = NULL,
                                          min_values_for_stats = NULL,
                                          changepoint = NULL, collector = NULL) {
  df <- as.data.frame(df)
  metrics <- temp_metrics()
  if (!all(c("TMIN", "TMAX", "water_year", "dowy") %in% names(df))) {
    return(.empty_family_stats(metrics, trend_completeness, decade_completeness,
                               min_values_for_stats, changepoint, collector))
  }
  frost_t <- as.numeric(pkg_env$climate_frost_temp_c %||% 0)
  end_month <- as.integer(pkg_env$climate_frost_last_end_month %||% 6L)
  base <- as.numeric(pkg_env$climate_gdd_base_c %||% 0)

  years <- unique(df$water_year)
  annual <- data.frame(water_year = years)
  for (m in metrics) annual[[m]] <- rep(NA_real_, nrow(annual))  # 0-row safe (schema contract)

  for (k in seq_along(years)) {
    yr <- as.integer(years[k])
    yd <- df[df$water_year == yr, , drop = FALSE]
    if (nrow(yd) == 0) next
    yd <- yd[order(yd$dowy), , drop = FALSE]
    TN <- .real_values(yd$TMIN) + 0; TX <- .real_values(yd$TMAX) + 0
    dowy <- as.numeric(yd$dowy)
    n <- length(TN)
    if (!.complete_year_grid(dowy, n, yr) || any(!is.finite(TN)) || any(!is.finite(TX))) next
    months <- .grid_months(yr, n)
    season_defs <- .climate_season_defs()

    imax <- which.max(TX); imin <- which.min(TN)
    annual$temp_max_annual[k] <- TX[imax]; annual$temp_max_annual_dowy[k] <- dowy[imax]
    annual$temp_min_annual[k] <- TN[imin]; annual$temp_min_annual_dowy[k] <- dowy[imin]

    gdd <- pmax(0, (TX + TN) / 2 - base)
    annual$gdd_cum_annual[k] <- seq_sum(gdd)
    for (s in .climate_season_keys) {
      smask <- months %in% season_defs[[s]]
      if (any(smask)) {
        annual[[paste0("temp_max_", s)]][k] <- max(TX[smask])
        annual[[paste0("temp_min_", s)]][k] <- min(TN[smask])
        annual[[paste0("gdd_cum_", s)]][k] <- seq_sum(gdd[smask])
      }
    }

    frost <- TN < frost_t
    f <- which(frost)
    if (length(f) > 0) annual$frost_first_dowy[k] <- dowy[f[1]]
    cal_year <- if (end_month >= 10L) yr - 1L else yr
    last_day <- as.Date(paste0(cal_year + (end_month == 12L), "-", sprintf("%02d", if (end_month == 12L) 1L else end_month + 1L), "-01")) - 1
    cutoff_dowy <- as.integer(last_day - as.Date(paste0(yr - 1L, "-10-01"))) + 1L
    l <- which(frost & dowy <= cutoff_dowy)
    if (length(l) > 0) annual$frost_last_dowy[k] <- dowy[l[length(l)]]
  }

  force_skip <- NULL
  if (isTRUE(pkg_env$climate_record_decade_gate) && !is.null(decade_completeness)) {
    force_skip <- .record_anchored_force_skip(annual, "temp_max_annual",
                                              .temp_record_gated_metrics, decade_completeness)
  }

  generate_stats(annual, value_cols = metrics, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness,
                 min_values_for_stats = min_values_for_stats,
                 changepoint = changepoint, collector = collector,
                 force_skip_trends = force_skip)
}

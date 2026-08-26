#' Snow metric names (14 per-water-year metrics)
#'
#' Port of julia/src/snow.jl (canonical). ALL metrics operate on a thresholded
#' series \code{SWE* = SWE if SWE >= swe_day_threshold_mm else 0} — days below
#' the threshold are treated as snow-free for durations AND magnitudes.
#' @export
SNOW_METRICS <- c(
  "swe_max", "swe_max_dowy", "snow_cover_days", "snow_on_dowy", "snow_off_dowy",
  "melt_season_days", "melt_rate", "ssm", "swe_apr1", "melt_before_peak",
  "melt_before_peak_pct", "melt_before_peak_to_max_swe", "melt_com_dowy",
  "swe_max_to_ppt"
)

# Threshold-dependent metrics gated by the record-anchored decade gate. The four
# magnitude metrics are NOT gated — their dense zero-including series carry the
# snow-decline signal legitimately.
SNOW_RECORD_GATED_METRICS <- c(
  "swe_max_dowy", "snow_on_dowy", "snow_off_dowy", "melt_season_days",
  "melt_rate", "ssm", "melt_before_peak", "melt_before_peak_pct",
  "melt_before_peak_to_max_swe", "melt_com_dowy"
)

#' Maximal runs of consecutive TRUE values
#'
#' @param mask Logical vector.
#' @return list of integer vectors (1-based index runs).
#' @keywords internal
snow_spells <- function(mask) {
  n <- length(mask)
  spells <- list()
  i <- 1L
  while (i <= n) {
    if (isTRUE(mask[i])) {
      j <- i
      while (j < n && isTRUE(mask[j + 1L])) j <- j + 1L
      spells[[length(spells) + 1L]] <- i:j
      i <- j + 1L
    } else {
      i <- i + 1L
    }
  }
  spells
}

#' Calculate 14 snow metrics with trend statistics from daily SWE
#'
#' Magnitude metrics (swe_max, snow_cover_days, swe_apr1, swe_max_to_ppt) emit
#' valid ZEROS for operationally snow-free years; timing/melt/regime metrics
#' emit NA. Callers must pass data already filtered to SWE-valid years
#' (\code{preprocess_daily_data}'s \code{valid_swe_years}); an SWE column in the
#' main gage frame is never used implicitly.
#'
#' @param streamflow_data data.frame with water_year, dowy, SWE (+ optional PPT).
#' @param valid_climate_years Years qualified for PPT use (swe_max_to_ppt), or NULL.
#' @param trend_completeness,decade_completeness,min_values_for_stats,changepoint,collector
#'   Forwarded to \code{generate_stats}.
#' @return Named list: 14 metrics x 8 statistics (+ 8 changepoint fields each).
#' @export
calculate_snow_metrics <- function(streamflow_data, valid_climate_years = NULL,
                                   trend_completeness = NULL,
                                   decade_completeness = NULL,
                                   min_values_for_stats = NULL,
                                   changepoint = NULL, collector = NULL) {
  df <- as.data.frame(streamflow_data)

  # Emit the IDENTICAL key set via the same explicit-value_cols machinery for
  # every degenerate input (schema contract).
  empty_snow_result <- function() {
    empty <- data.frame(water_year = integer(0))
    for (m in SNOW_METRICS) empty[[m]] <- numeric(0)
    generate_stats(empty, value_cols = SNOW_METRICS, year_col = "water_year",
                   trend_completeness = trend_completeness,
                   decade_completeness = decade_completeness,
                   min_values_for_stats = min_values_for_stats,
                   changepoint = changepoint, collector = collector)
  }

  # A 0-ROW frame is a legitimate input, not an error: the orchestrator passes an
  # explicit snow_data filtered to valid_swe_years, which is empty for a gage that
  # has SWE columns but no SWE-valid year (e.g. warm-climate gages). Without this
  # branch the code below builds a 0-row `annual` and then assigns a length-1
  # NA into it, raising "replacement has 1 row, data has 0" — swallowed by the
  # orchestrator's try/catch, so the gage silently lost all 224 snow columns.
  # Caught 2026-08-26 by check_signature_failures.py on 4 Florida gages.
  if (!"SWE" %in% names(df) || nrow(df) == 0L) {
    return(empty_snow_result())
  }

  thr          <- pkg_env$snow_swe_threshold_mm
  seasonal_min <- pkg_env$snow_seasonal_min_days
  com_frac     <- pkg_env$snow_melt_com_fraction
  min_ppt      <- pkg_env$snow_min_annual_ppt_mm

  has_ppt <- "PPT" %in% names(df)
  years <- unique(df$water_year)
  annual <- data.frame(water_year = years)
  for (m in SNOW_METRICS) annual[[m]] <- NA_real_

  for (k in seq_along(years)) {
    yr <- years[k]
    yd <- df[df$water_year == yr, , drop = FALSE]
    if (nrow(yd) == 0) next
    yd <- yd[order(yd$dowy), , drop = FALSE]

    swe_raw <- as.numeric(yd$SWE)
    dowy <- as.numeric(yd$dowy)
    n <- length(swe_raw)

    # Complete-grid + no-NA guard: full water years whose days are EXACTLY
    # 1..n (n in {365, 366}) with no missing SWE
    if (!(n == 365L || n == 366L) || any(is.na(dowy)) ||
        !identical(as.integer(dowy), seq_len(n)) || any(is.na(swe_raw))) next

    swe_star <- ifelse(swe_raw >= thr, swe_raw, 0)
    mask <- swe_star > 0
    n_snow <- sum(mask)

    # --- Magnitude metrics (valid zeros when operationally snow-free) ---
    swe_max_val <- if (n_snow > 0) max(swe_star) else 0
    annual$swe_max[k] <- swe_max_val
    annual$snow_cover_days[k] <- as.numeric(n_snow)

    # April 1 SWE (leap-safe; dowy == 1..n after the grid guard)
    apr1_dowy <- as.integer(as.Date(paste0(yr, "-04-01")) -
                            as.Date(paste0(yr - 1, "-10-01"))) + 1L
    annual$swe_apr1[k] <- swe_star[apr1_dowy]

    # --- Timing / melt / regime metrics ---
    if (n_snow > 0) {
      peak_idx <- which.max(swe_star)          # first occurrence on ties
      annual$swe_max_dowy[k] <- dowy[peak_idx]

      spells <- snow_spells(mask)
      lens <- vapply(spells, length, integer(1))
      seasonal_days  <- sum(lens[lens >= seasonal_min])
      ephemeral_days <- sum(lens[lens <  seasonal_min])
      annual$ssm[k] <- (seasonal_days - ephemeral_days) / n_snow

      anchor <- spells[[which(vapply(spells, function(sp) peak_idx %in% sp, logical(1)))[1]]]
      a_first <- anchor[1]; a_last <- anchor[length(anchor)]

      snow_on  <- if (a_first == 1L) NA_real_ else dowy[a_first]
      snow_off <- if (a_last == n)   NA_real_ else dowy[a_last] + 1
      annual$snow_on_dowy[k]  <- snow_on
      annual$snow_off_dowy[k] <- snow_off

      if (!is.na(snow_off)) {
        melt_days <- snow_off - dowy[peak_idx]
        annual$melt_season_days[k] <- melt_days
        annual$melt_rate[k] <- swe_max_val / melt_days
      }

      # Melt increments attributed to the day the drop lands (within-year only)
      m <- numeric(n)
      drops <- swe_star[-n] - swe_star[-1]
      m[-1] <- ifelse(drops > 0, drops, 0)
      total_melt <- sum(m)
      melt_before <- sum(m[seq_len(peak_idx)])
      annual$melt_before_peak[k] <- melt_before
      annual$melt_before_peak_to_max_swe[k] <- melt_before / swe_max_val
      if (total_melt > 0) {
        annual$melt_before_peak_pct[k] <- 100 * melt_before / total_melt
        cum <- cumsum(m)
        reached <- which(cum >= com_frac * total_melt)
        if (length(reached) > 0) annual$melt_com_dowy[k] <- dowy[reached[1]]
      }
    }

    # --- swe_max_to_ppt: needs a PPT-qualified year ---
    if (has_ppt && (is.null(valid_climate_years) || yr %in% valid_climate_years)) {
      P <- as.numeric(yd$PPT)
      pv <- !is.na(P)
      total_P <- if (any(pv)) sum(P[pv]) else 0
      if (total_P > min_ppt) annual$swe_max_to_ppt[k] <- swe_max_val / total_P
    }
  }

  # Record-anchored decade gate: anchored to the SWE-valid, grid-complete years
  # (swe_max non-NA — dense incl. zeros), NOT the metric's own span. Threshold is
  # the SAME decade_completeness knob the streamflow gate uses (linked by design).
  force_skip <- NULL
  if (isTRUE(pkg_env$snow_record_decade_gate) && !is.null(decade_completeness)) {
    anchor_mask <- !is.na(annual$swe_max)
    anchor_years <- as.integer(annual$water_year[anchor_mask])
    if (length(anchor_years) > 0) {
      rec_min <- min(anchor_years); rec_max <- max(anchor_years)
      if (rec_max - rec_min + 1 >= 10) {
        windows <- list(c(rec_min, rec_min + 9), c(rec_max - 9, rec_max))
        all_years <- as.integer(annual$water_year)
        skip <- character(0)
        for (mn in SNOW_RECORD_GATED_METRICS) {
          vals <- annual[[mn]]
          for (w in windows) {
            in_win <- anchor_mask & all_years >= w[1] & all_years <= w[2]
            den <- sum(in_win)
            if (den == 0) next
            num <- sum(!is.na(vals[in_win]))
            if ((num / den) < decade_completeness) { skip <- c(skip, mn); break }
          }
        }
        if (length(skip) > 0) force_skip <- skip
      }
    }
  }

  generate_stats(annual, value_cols = SNOW_METRICS, year_col = "water_year",
                 trend_completeness = trend_completeness,
                 decade_completeness = decade_completeness,
                 min_values_for_stats = min_values_for_stats,
                 changepoint = changepoint, collector = collector,
                 force_skip_trends = force_skip)
}

#' Weibull-plotting-position quantile (Hyndman-Fan definition 6)
#'
#' Quantile at non-exceedance probability \code{p} using \code{p_i = i/(n+1)}.
#' Port of \code{weibull_quantile} in julia/src/drought.jl. NAs are dropped.
#' \code{p} below the smallest plotting position \code{1/(n+1)} is NOT
#' interpolable: with \code{below_range = "na"} the result is NA (honest — the
#' sample cannot resolve that probability); \code{"clamp"} returns the minimum.
#'
#' @param x Numeric vector.
#' @param p Non-exceedance probability in (0, 1).
#' @param below_range "na" (default) or "clamp".
#' @return Numeric quantile, or NA.
#' @export
weibull_quantile <- function(x, p, below_range = "na") {
  vals <- sort(x[!is.na(x)])
  n <- length(vals)
  if (n == 0) return(NA_real_)
  h <- (n + 1) * p
  if (h < 1) return(if (below_range == "clamp") vals[1] else NA_real_)
  if (h >= n) return(vals[n])
  fl <- floor(h)
  vals[fl] + (h - fl) * (vals[fl + 1] - vals[fl])
}

#' Moving-average smoothing within maximal runs of consecutive dates
#'
#' Port of \code{smooth_daily_flow} in julia/src/drought.jl. The window never
#' averages across a temporal gap (rejected years, duplicated dates or missing
#' dates all break a run); it shrinks at run edges, and a day whose in-run
#' window holds fewer than \code{min_valid} non-NA values yields NA.
#' \code{dates} must be sorted ascending (callers sort first).
#'
#' @param dates Date vector (sorted ascending).
#' @param Q Numeric vector, same length as dates.
#' @param window,alignment,min_valid Smoothing parameters (default from config).
#' @return Numeric vector of smoothed values.
#' @export
smooth_daily_flow <- function(dates, Q, window = NULL, alignment = NULL,
                              min_valid = NULL) {
  if (is.null(window))    window    <- pkg_env$drought_smooth_window
  if (is.null(alignment)) alignment <- pkg_env$drought_smooth_alignment
  if (is.null(min_valid)) min_valid <- pkg_env$drought_smooth_min_valid

  Q <- as.numeric(Q)
  n <- length(Q)
  out <- rep(NA_real_, n)
  if (n == 0) return(out)
  if (length(dates) != n) {
    stop(sprintf("smooth_daily_flow: dates and Q length mismatch (%d vs %d)",
                 length(dates), n))
  }

  d <- suppressWarnings(as.Date(dates))
  ok <- !is.na(d)
  nums <- rep(0L, n)
  nums[ok] <- as.integer(d[ok])
  half <- (window - 1L) %/% 2L

  run_start <- 1L
  for (i in seq_len(n)) {
    ends_run <- (i == n) || !ok[i] || !ok[i + 1L] || (nums[i + 1L] != nums[i] + 1L)
    if (!ends_run) next
    # A row with a missing/unparseable date belongs to NO run
    run_end <- if (ok[i]) i else i - 1L
    if (run_end >= run_start) {
      for (j in run_start:run_end) {
        if (alignment == "center") {
          lo <- max(run_start, j - half); hi <- min(run_end, j + half)
        } else {
          lo <- max(run_start, j - window + 1L); hi <- j
        }
        # Accumulate SEQUENTIALLY in double precision, mirroring the loop in
        # julia/src/drought.jl exactly. R's mean() uses a long-double accumulator
        # plus a correction pass, so it returns a different last bit than Julia's
        # `s += v` — measured on gage 01589795: 4,513 of 10,227 smoothed values
        # differed, by at most 3.6e-15. Harmless in isolation, but the fixed
        # thresholds are percentiles OF this series, so when a threshold lands on
        # a flow plateau a last-bit shift flips every plateau day at once through
        # the strict `<` (that gage's WY2002 drought duration read 116 in Julia
        # vs 60 here). numpy's mean happens to match Julia's order for windows
        # this short, which is why the Python port never showed this.
        # This is also marginally FASTER: no per-day subvector allocation.
        s <- 0
        cnt <- 0L
        for (k in lo:hi) {
          vk <- Q[k]
          if (!is.na(vk)) {
            s <- s + vk
            cnt <- cnt + 1L
          }
        }
        out[j] <- if (cnt >= min_valid) s / cnt else NA_real_
      }
    }
    run_start <- i + 1L
  }
  out
}

#' Streamflow drought duration and deficit (Adelsperger et al., in review)
#'
#' Port of \code{calculate_drought_metrics} in julia/src/drought.jl. For each
#' configured percentile level p: \code{drought_duration_fixed_p{p}} (days per
#' water year below the fixed whole-record threshold) and
#' \code{drought_deficit_fixed_p{p}} (summed departures below it). Also emits
#' per-gage scalar diagnostics \code{drought_threshold_fixed_p{p}}.
#'
#' Both metrics are dense series whose ZEROS are meaningful, so the trend gate
#' and stats floor apply normally (this family is NOT exempt).
#'
#' @param streamflow_data data.frame with Q, water_year, date.
#' @param trend_completeness,decade_completeness,min_values_for_stats,changepoint,collector
#'   Forwarded to \code{generate_stats}.
#' @return Named list of statistics plus the threshold scalars.
#' @export
calculate_drought_metrics <- function(streamflow_data, trend_completeness = NULL,
                                      decade_completeness = NULL,
                                      min_values_for_stats = NULL,
                                      changepoint = NULL, collector = NULL) {
  pcts <- pkg_env$drought_percentiles
  dur_metrics <- paste0("drought_duration_fixed_p", pcts)
  def_metrics <- paste0("drought_deficit_fixed_p", pcts)
  metrics <- c(dur_metrics, def_metrics)
  scalars <- paste0("drought_threshold_fixed_p", pcts)

  all_nan_result <- function() {
    empty <- data.frame(water_year = integer(0))
    for (m in metrics) empty[[m]] <- numeric(0)
    res <- generate_stats(empty, value_cols = metrics, year_col = "water_year",
                          trend_completeness = trend_completeness,
                          decade_completeness = decade_completeness,
                          min_values_for_stats = min_values_for_stats,
                          changepoint = changepoint, collector = collector)
    for (s in scalars) res[[s]] <- NA_real_
    res
  }

  df <- as.data.frame(streamflow_data)
  # rpkg's preprocess_daily_data() renames `date` -> `Date` internally (io.R),
  # while Julia/Python keep lowercase throughout. Accept EITHER so the module
  # works on both a raw frame and a preprocessor output — without this, drought
  # silently returned an all-NA family for every real gage (caught 2026-08-26 by
  # an end-to-end smoke on real data; the unit tests built frames with `date`
  # directly and so never exercised the preprocessor -> drought path).
  if (!"date" %in% names(df) && "Date" %in% names(df)) {
    df$date <- df$Date
  }
  needed <- c("Q", "water_year", "date")
  if (length(setdiff(needed, names(df))) > 0) {
    warning(paste("calculate_drought_metrics: Missing columns:",
                  paste(setdiff(needed, names(df)), collapse = ", ")))
    return(all_nan_result())
  }
  if (nrow(df) == 0) return(all_nan_result())

  dfs <- df[order(df$date), , drop = FALSE]
  Q_smooth <- smooth_daily_flow(dfs$date, dfs$Q)

  wy <- suppressWarnings(as.numeric(dfs$water_year))
  years <- sort(unique(as.integer(wy[!is.na(wy)])))
  if (length(years) == 0) return(all_nan_result())

  # Fixed thresholds: whole-record percentiles of the smoothed series. Short
  # records get NA thresholds rather than an unstable low tail.
  thresholds <- rep(NA_real_, length(pcts))
  if (length(years) >= pkg_env$drought_min_years) {
    pool <- Q_smooth[!is.na(Q_smooth)]
    if (length(pool) > 0) {
      for (k in seq_along(pcts)) {
        thresholds[k] <- weibull_quantile(pool, pcts[k] / 100,
                                          below_range = pkg_env$drought_below_range_policy)
      }
    }
  }

  annual <- data.frame(water_year = years)
  for (m in metrics) annual[[m]] <- NA_real_

  for (i in seq_along(years)) {
    yr <- years[i]
    sel <- !is.na(wy) & wy == yr
    if (!any(sel)) next
    qs <- Q_smooth[sel]
    qs <- qs[!is.na(qs)]
    for (k in seq_along(pcts)) {
      thr <- thresholds[k]
      if (is.na(thr)) next
      below <- qs[qs < thr]                      # STRICT comparison
      annual[[dur_metrics[k]]][i] <- length(below)
      annual[[def_metrics[k]]][i] <- sum(thr - below)
    }
  }

  res <- generate_stats(annual, value_cols = metrics, year_col = "water_year",
                        trend_completeness = trend_completeness,
                        decade_completeness = decade_completeness,
                        min_values_for_stats = min_values_for_stats,
                        changepoint = changepoint, collector = collector)
  for (k in seq_along(scalars)) res[[scalars[k]]] <- thresholds[k]
  res
}

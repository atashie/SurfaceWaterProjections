#' Pettitt non-parametric changepoint test
#'
#' Port of the pipeline-reachable subset of \code{julia/src/changepoint.jl}
#' (canonical): \code{pettitt_test} and \code{segment_differential_metrics},
#' called from \code{generate_stats}'s changepoint block. The Julia module's
#' BIC model-selection path (\code{detect_changepoint}) is exported Julia-only
#' API that the production pipeline never calls — intentionally excluded from
#' parity scope (port plan 2026-08-24 §3).
#'
#' Semantics pinned to the canonical implementation:
#' \itemize{
#'   \item NaN/NA pairs filtered, then sorted by year (R and Julia are both
#'     1-based, so the index arithmetic translates directly).
#'   \item \code{U_t} via the recursion \code{U_t = U_\{t-1\} + V_t} with
#'     \code{V_t = sum_j sign(x_t - x_j)}.
#'   \item The changepoint is the FIRST \code{t} attaining the strict maximum of
#'     \code{|U_t|} over \code{t in [min_segment_obs, n - min_segment_obs]} —
#'     ties select the earliest candidate (Julia updates only on \code{>}).
#'   \item Asymptotic p-value \code{p = min(1, 2 exp(-6 K^2 / (n^3 + n^2)))}
#'     (Pettitt 1979).
#' }
#'
#' @param years Numeric vector of years.
#' @param values Numeric vector of annual metric values.
#' @param min_total_obs Minimum valid observations for the test to run.
#' @param min_segment_obs Minimum observations either side of a candidate.
#' @return Named list with \code{cp_year}, \code{pval}, \code{pre_mean},
#'   \code{post_mean} (all NA when there is insufficient data).
#' @export
pettitt_test <- function(years, values, min_total_obs = 20, min_segment_obs = 10) {
  na_result <- list(cp_year = NA_real_, pval = NA_real_,
                    pre_mean = NA_real_, post_mean = NA_real_)

  years <- as.numeric(years)
  values <- as.numeric(values)
  valid <- !is.na(years) & !is.na(values)
  yrs <- years[valid]
  vals <- values[valid]
  n <- length(yrs)

  if (n < min_total_obs) return(na_result)

  ord <- order(yrs)          # stable (order() uses a stable sort by default)
  yrs <- yrs[ord]
  vals <- vals[ord]

  # V_t = sum_j sign(x_t - x_j); U = cumsum(V). n <= ~46 for annual series.
  V <- rowSums(sign(outer(vals, vals, "-")))
  U <- cumsum(V)

  # Candidates t in [min_segment_obs, n - min_segment_obs]; which.max returns
  # the FIRST maximum, matching Julia's strict `>` update rule.
  lo <- min_segment_obs
  hi <- n - min_segment_obs
  if (hi < lo) return(na_result)
  seg <- abs(U[lo:hi])
  best_t <- lo + which.max(seg) - 1L
  K <- seg[which.max(seg)]

  pval <- 2 * exp(-6 * K^2 / (n^3 + n^2))
  pval <- min(pval, 1)

  list(cp_year   = yrs[best_t],
       pval      = pval,
       pre_mean  = mean(vals[seq_len(best_t)]),
       post_mean = mean(vals[(best_t + 1L):n]))
}


#' Differential metrics at a changepoint
#'
#' Computes delta_mean, pct_change, and the pre/post Mann-Kendall p-values at a
#' given changepoint year. Mirrors \code{segment_differential_metrics} in
#' \code{julia/src/changepoint.jl}, including the >= 3 values-per-segment guard
#' and the \code{abs(pre_mean) > 1e-10} percentage guard. Each language reuses
#' its own MK implementation (here \code{mann_kendall_test}, the wrapper
#' \code{generate_stats} uses, so constant segments give NA as in Julia).
#'
#' @param years Numeric vector of years.
#' @param values Numeric vector of values.
#' @param cp_year Changepoint year (NA yields an all-NA result).
#' @return Named list with \code{delta_mean}, \code{pct_change},
#'   \code{pre_mk_pval}, \code{post_mk_pval}.
#' @export
segment_differential_metrics <- function(years, values, cp_year) {
  na_result <- list(delta_mean = NA_real_, pct_change = NA_real_,
                    pre_mk_pval = NA_real_, post_mk_pval = NA_real_)

  years <- as.numeric(years)
  values <- as.numeric(values)
  if (is.na(cp_year) || length(years) == 0) return(na_result)

  pre_vals  <- values[years <= cp_year]
  post_vals <- values[years >  cp_year]
  if (length(pre_vals) < 3 || length(post_vals) < 3) return(na_result)

  pre_m  <- mean(pre_vals)
  post_m <- mean(post_vals)
  delta_mean <- post_m - pre_m
  pct_change <- if (abs(pre_m) > 1e-10) delta_mean / abs(pre_m) * 100 else NA_real_

  # Same NA contract as generate_stats (n < 3 or constant segment -> NA)
  mk_p <- function(v) mann_kendall_test(v)$pval

  list(delta_mean = delta_mean, pct_change = pct_change,
       pre_mk_pval = mk_p(pre_vals), post_mk_pval = mk_p(post_vals))
}


# Changepoint column suffixes (Pettitt test only) — mirrors CP_SUFFIXES in
# julia/src/stats.jl
CP_SUFFIXES <- c("_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
                 "_pettitt_post_mean", "_pettitt_delta_mean",
                 "_pettitt_pct_change", "_pettitt_pre_mk_pval",
                 "_pettitt_post_mk_pval")


#' Internal: NA changepoint entries for one metric
#' @param metric Base metric name.
#' @return Named list of 8 NA_real_ entries.
#' @keywords internal
empty_changepoint <- function(metric) {
  out <- as.list(rep(NA_real_, length(CP_SUFFIXES)))
  names(out) <- paste0(metric, CP_SUFFIXES)
  out
}


#' Internal: run the Pettitt block for one metric
#'
#' Mirrors \code{_run_changepoint_block!} in julia/src/stats.jl: filter to the
#' changepoint analysis window, run Pettitt (4 fields), then the differential
#' metrics (4 fields).
#'
#' @param col Metric name.
#' @param valid_years,valid_values The metric's non-NA annual series.
#' @param changepoint Named list with start_year, end_year, min_total_obs,
#'   min_segment_obs.
#' @return Named list of 8 changepoint entries.
#' @keywords internal
run_changepoint_block <- function(col, valid_years, valid_values, changepoint) {
  cp_start   <- if (!is.null(changepoint$start_year))      changepoint$start_year      else 1980
  cp_end     <- if (!is.null(changepoint$end_year))        changepoint$end_year        else 2024
  cp_min_obs <- if (!is.null(changepoint$min_total_obs))   changepoint$min_total_obs   else 20
  cp_min_seg <- if (!is.null(changepoint$min_segment_obs)) changepoint$min_segment_obs else 10

  in_win <- valid_years >= cp_start & valid_years <= cp_end
  cp_years <- valid_years[in_win]
  cp_vals  <- valid_values[in_win]

  pt <- pettitt_test(cp_years, cp_vals,
                     min_total_obs = cp_min_obs, min_segment_obs = cp_min_seg)
  df <- segment_differential_metrics(cp_years, cp_vals, pt$cp_year)

  out <- list(pt$cp_year, pt$pval, pt$pre_mean, pt$post_mean,
              df$delta_mean, df$pct_change, df$pre_mk_pval, df$post_mk_pval)
  names(out) <- paste0(col, CP_SUFFIXES)
  out
}

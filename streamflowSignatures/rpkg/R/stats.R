#' Generate 8 summary statistics for each metric column
#'
#' For each value column, produces Theil-Sen slope, linear slope, Spearman
#' correlation and p-value, Mann-Kendall tau and p-value, mean, and median.
#'
#' @param data A data.frame or data.table with a year column and value columns.
#' @param value_cols Character vector of column names to compute stats for.
#'   If NULL, uses all numeric columns except year_col.
#' @param year_col Name of the year column (default "water_year").
#' @param min_rows Minimum non-NA rows required for trend statistics.
#' @param trend_completeness Numeric in (0,1] or NULL. If non-NULL, trend
#'   statistics are set to NA when the fraction of non-NA values falls below
#'   this threshold. Mean and median are still computed.
#' @param decade_completeness Numeric in (0,1] or NULL. If non-NULL (and
#'   \code{trend_completeness} is also non-NULL), the first and last decades
#'   of the series must each have at least this fraction of non-NA values,
#'   otherwise trend statistics are set to NA.
#' @param changepoint Named list (start_year, end_year, min_total_obs,
#'   min_segment_obs) or NULL. When supplied, 8 additional \code{_pettitt_*}
#'   fields are emitted per metric. Changepoint analysis runs INDEPENDENTLY of
#'   the trend-completeness gates (matching Julia): a trend-gated metric still
#'   gets its Pettitt fields, while a metric below \code{min_rows} or the stats
#'   floor gets NA Pettitt fields.
#' @param min_values_for_stats Integer or NULL. Stats floor: a metric with
#'   fewer non-NA annual values emits NA for ALL 8 statistics AND its
#'   changepoint fields. Recession and elasticity are exempt at the
#'   orchestration layer (the argument is never passed to them).
#' @param force_skip_trends Character vector of metric names whose 6 trend
#'   statistics are suppressed through the SAME path as trend_completeness —
#'   mean/median and changepoint unaffected.
#' @param collector Optional \code{annual_collector()} environment. When
#'   supplied, each metric's annual series is appended in long format BEFORE
#'   any gating (export contract). Read-only with respect to the statistics.
#' @return Named list of statistics (metric_stat = value).
#' @export
generate_stats <- function(data, value_cols = NULL, year_col = "water_year",
                           min_rows = 3, trend_completeness = NULL,
                           decade_completeness = NULL, changepoint = NULL,
                           min_values_for_stats = NULL,
                           force_skip_trends = NULL, collector = NULL) {
  if (!is.data.frame(data)) stop("Input 'data' must be a data.frame or data.table")
  if (!year_col %in% colnames(data)) stop(paste("Year column", year_col, "not found"))

  # Resolve value columns BEFORE sorting/gating (mirrors Julia's hoisted
  # resolution — column types are unaffected by row order)
  if (is.null(value_cols)) {
    numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
    value_cols <- setdiff(numeric_cols, year_col)
  }

  results <- list()

  # Opt-in annual-values collection: capture each series exactly as passed,
  # BEFORE any gating below (matching Julia's placement).
  if (!is.null(collector)) {
    for (col in value_cols) {
      if (col %in% colnames(data)) {
        collect_annual(collector, data[[year_col]], data[[col]], col)
      }
    }
  }

  # Frame-level gate (mirrors Julia's earliest nrow < min_rows return): every
  # REQUESTED metric — including requested-but-absent columns — emits the 8 NA
  # statistics (+ NA changepoint fields when enabled).
  if (nrow(data) < min_rows) {
    for (col in value_cols) {
      results <- c(results, empty_stats(col))
      if (!is.null(changepoint)) results <- c(results, empty_changepoint(col))
    }
    return(results)
  }

  data <- data[order(data[[year_col]]), , drop = FALSE]

  for (col in value_cols) {
    if (!col %in% colnames(data)) {
      warning(paste("Column", col, "not found in data. Skipping."))
      next
    }

    years  <- data[[year_col]]
    values <- data[[col]]
    valid  <- !is.na(values)
    w_years  <- years[valid]
    w_values <- values[valid]

    # Stats floor gates through the same branch as min_rows: below the floor,
    # ALL 8 statistics AND the changepoint fields are NA (matching Julia).
    if (length(w_values) < min_rows ||
        (!is.null(min_values_for_stats) && length(w_values) < min_values_for_stats)) {
      results <- c(results, empty_stats(col))
      if (!is.null(changepoint)) results <- c(results, empty_changepoint(col))
      next
    }

    # Trend completeness check: uses metric's own non-NA year range for decades
    skip_trends <- FALSE
    # Externally-forced trend skip (e.g. the snow record-anchored decade gate):
    # suppresses the 6 trend stats; mean/median + changepoint unaffected.
    if (!is.null(force_skip_trends) && col %in% force_skip_trends) skip_trends <- TRUE
    if (!skip_trends && !is.null(trend_completeness)) {
      all_values <- data[[col]]
      all_years  <- data[[year_col]]
      valid_mask <- !is.na(all_values)
      valid_yrs  <- all_years[valid_mask]
      if (length(valid_yrs) > 0) {
        metric_min <- min(valid_yrs)
        metric_max <- max(valid_yrs)
        metric_span <- metric_max - metric_min + 1
        if (metric_span > 0 && (length(valid_yrs) / metric_span) < trend_completeness) skip_trends <- TRUE
        if (!skip_trends && !is.null(decade_completeness) && metric_span >= 10) {
          first_yrs <- valid_yrs[valid_yrs >= metric_min & valid_yrs <= (metric_min + 9)]
          first_exp <- min(10, metric_span)
          if (first_exp > 0 && (length(first_yrs) / first_exp) < decade_completeness) skip_trends <- TRUE
          if (!skip_trends) {
            last_start <- metric_max - 9
            last_yrs <- valid_yrs[valid_yrs >= last_start & valid_yrs <= metric_max]
            last_exp <- min(10, metric_max - max(last_start, metric_min) + 1)
            if (last_exp > 0 && (length(last_yrs) / last_exp) < decade_completeness) skip_trends <- TRUE
          }
        }
      }
    }
    if (skip_trends) {
      results[[paste0(col, "_senn_slp")]]      <- NA_real_
      results[[paste0(col, "_linear_slp")]]    <- NA_real_
      results[[paste0(col, "_spearman_rho")]]  <- NA_real_
      results[[paste0(col, "_spearman_pval")]] <- NA_real_
      results[[paste0(col, "_mk_rho")]]        <- NA_real_
      results[[paste0(col, "_mk_pval")]]       <- NA_real_
      results[[paste0(col, "_mean")]]          <- mean(w_values, na.rm = TRUE)
      results[[paste0(col, "_median")]]        <- median(w_values, na.rm = TRUE)
      # Changepoint runs INDEPENDENTLY of trend gating (matching Julia)
      if (!is.null(changepoint)) {
        results <- c(results, run_changepoint_block(col, w_years, w_values, changepoint))
      }
      next
    }

    # Theil-Sen slope
    sen_slope <- tryCatch({
      res <- zyp::zyp.sen(w_values ~ w_years)
      unname(res$coefficients[2])
    }, error = function(e) NA_real_)

    # Linear regression slope
    linear_slp <- tryCatch({
      fit <- lm(w_values ~ w_years)
      unname(coef(fit)[2])
    }, error = function(e) NA_real_)

    # Spearman correlation
    sp <- tryCatch({
      cor.test(w_years, w_values, method = "spearman")
    }, error = function(e) NULL)
    sp_rho  <- if (!is.null(sp)) unname(sp$estimate) else NA_real_
    sp_pval <- if (!is.null(sp)) sp$p.value          else NA_real_

    # Mann-Kendall
    mk <- tryCatch({
      Kendall::MannKendall(w_values)
    }, error = function(e) NULL)
    mk_tau  <- if (!is.null(mk)) as.numeric(mk$tau) else NA_real_
    mk_pval <- if (!is.null(mk)) as.numeric(mk$sl)  else NA_real_

    mean_val   <- mean(w_values, na.rm = TRUE)
    median_val <- median(w_values, na.rm = TRUE)

    results[[paste0(col, "_senn_slp")]]      <- sen_slope
    results[[paste0(col, "_linear_slp")]]    <- linear_slp
    results[[paste0(col, "_spearman_rho")]]  <- sp_rho
    results[[paste0(col, "_spearman_pval")]] <- sp_pval
    results[[paste0(col, "_mk_rho")]]        <- mk_tau
    results[[paste0(col, "_mk_pval")]]       <- mk_pval
    results[[paste0(col, "_mean")]]          <- mean_val
    results[[paste0(col, "_median")]]        <- median_val

    # Changepoint analysis (opt-in)
    if (!is.null(changepoint)) {
      results <- c(results, run_changepoint_block(col, w_years, w_values, changepoint))
    }
  }

  results
}

#' Return a named list of 8 NA statistics for a metric
#'
#' @param metric Base metric name.
#' @return Named list with 8 NA_real_ entries.
#' @export
empty_stats <- function(metric) {
  suffixes <- c("_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                "_mk_rho", "_mk_pval", "_mean", "_median")
  out <- as.list(rep(NA_real_, 8))
  names(out) <- paste0(metric, suffixes)
  out
}


# ---------------------------------------------------------------------------
# Annual-values collector (Phase 3 of the port campaign)
# ---------------------------------------------------------------------------

#' Create an annual-values collector
#'
#' Opt-in accumulator for per-year signature values in long format. Port of the
#' \code{AnnualCollector} struct in julia/src/stats.jl. When passed to
#' \code{generate_stats}, the annual series is appended EXACTLY as received —
#' before any min_rows, stats-floor or trend-completeness gating — so the
#' collected rows are precisely what the summary statistics were computed from.
#' Collection is read-only with respect to the statistics.
#'
#' Uses an environment holding lists of per-metric chunks, appended by index
#' assignment (\code{e$sig[[k]] <- ...}) and flattened once at drain time —
#' avoiding the per-row \code{c()} append that would be O(n^2) in R. One
#' instance per gage; the benchmark runner drains it with
#' \code{collector_drain()} and tags the rows with gage_id.
#'
#' @return An environment with fields \code{sig}, \code{wy}, \code{val}
#'   (lists of accumulated chunks) and \code{n} (chunk count).
#' @export
annual_collector <- function() {
  e <- new.env(parent = emptyenv())
  e$sig <- list(); e$wy <- list(); e$val <- list(); e$n <- 0L
  e
}

#' Drain a collector into a long-format data.frame
#'
#' @param collector An \code{annual_collector()} environment.
#' @return data.frame with columns signature, water_year, value (0 rows if empty).
#' @export
collector_drain <- function(collector) {
  if (is.null(collector) || collector$n == 0L) {
    return(data.frame(signature = character(0), water_year = integer(0),
                      value = numeric(0), stringsAsFactors = FALSE))
  }
  data.frame(
    signature  = unlist(collector$sig, use.names = FALSE),
    water_year = as.integer(unlist(collector$wy, use.names = FALSE)),
    value      = as.numeric(unlist(collector$val, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

#' Internal: append one metric's annual series to a collector
#'
#' Rows whose YEAR is missing/non-finite are skipped (they cannot be keyed);
#' NA values are KEPT — they record "year present, metric not computable".
#' @keywords internal
collect_annual <- function(collector, years, values, col) {
  keep <- !is.na(years) & is.finite(years)
  if (!any(keep)) return(invisible(NULL))
  k <- collector$n + 1L
  collector$sig[[k]] <- rep(col, sum(keep))
  collector$wy[[k]]  <- as.integer(round(years[keep]))
  collector$val[[k]] <- as.numeric(values[keep])
  collector$n <- k
  invisible(NULL)
}

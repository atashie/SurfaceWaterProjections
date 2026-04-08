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
#' @return Named list of statistics (metric_stat = value).
#' @export
generate_stats <- function(data, value_cols = NULL, year_col = "water_year",
                           min_rows = 3, trend_completeness = NULL,
                           decade_completeness = NULL) {
  if (!is.data.frame(data)) stop("Input 'data' must be a data.frame or data.table")
  if (!year_col %in% colnames(data)) stop(paste("Year column", year_col, "not found"))

  data <- data[order(data[[year_col]]), , drop = FALSE]

  if (is.null(value_cols)) {
    numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
    value_cols <- setdiff(numeric_cols, year_col)
  }

  results <- list()

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

    if (length(w_values) < min_rows) {
      results <- c(results, empty_stats(col))
      next
    }

    # Trend completeness check: uses metric's own non-NA year range for decades
    skip_trends <- FALSE
    if (!is.null(trend_completeness)) {
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

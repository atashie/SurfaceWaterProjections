#' Calculate Richards-Baker flashiness index trends
#'
#' @param streamflow_data A data.frame with columns: water_year, Q, dowy.
#' @return Named list of 8 statistics for flashinessRB.
#' @export
analyze_flashiness_trends <- function(streamflow_data,
                                     trend_completeness = NULL,
                                     decade_completeness = NULL, min_values_for_stats = NULL, changepoint = NULL, collector = NULL) {
  required_cols <- c("water_year", "Q")
  missing <- setdiff(required_cols, colnames(streamflow_data))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))

  # Zero-row input is legitimate (an orchestrator may pass an empty frame) and
  # canonical Julia guards it explicitly. Without this the code below builds a
  # 0-row frame and assigns a length-1 NA into it, raising "replacement has 1
  # row, data has 0" — which safe_call() swallows, silently dropping the whole
  # family for that gage (exactly the snow defect found 2026-08-26).
  if (nrow(streamflow_data) == 0L) return(empty_stats("flashinessRB"))

  years <- unique(streamflow_data$water_year)

    # Metric columns use the CANONICAL names directly (matching Julia, which never
  # renames): the annual-values collector captures the column name it is given,
  # so a placeholder + post-hoc rename would put the wrong signature name in the
  # annual parquet (found 2026-08-25, Phase 3).
  flashiness_by_year <- data.frame(water_year = years, flashinessRB = NA_real_)

  for (yr in years) {
    year_data <- streamflow_data[streamflow_data$water_year == yr, ]

    if ("dowy" %in% colnames(year_data)) {
      year_data <- year_data[order(year_data$dowy), ]
    }

    q_values <- year_data$Q
    stopifnot(!any(is.na(q_values)))

    total_q <- sum(q_values)
    # <= 0 matches canonical Julia: zero or NEGATIVE total flow (negative-Q
    # days retained by default) is non-computable (2026-08-24 alignment fix)
    if (total_q <= 0) next

    rb_index <- sum(abs(diff(q_values))) / total_q
    flashiness_by_year$flashinessRB[flashiness_by_year$water_year == yr] <- rb_index
  }

  # Julia APPENDS only computable years, so a non-computable year is ABSENT
  # from its frame. rpkg pre-allocates a row per year and assigns into it, so a
  # skipped year left an NA row behind — which the annual collector then exported
  # as a row Julia never emits (1,039 such rows across the sparse families,
  # caught 2026-08-26 by the cross-language annual gate). Drop them to match.
  flashiness_by_year <- flashiness_by_year[!is.na(flashiness_by_year$flashinessRB), , drop = FALSE]

  # Julia returns empty_stats outright below 3 computable years
  # (julia/src/flashiness.jl:100), so generate_stats is never reached and the
  # collector receives NO rows for that gage. Without this, rpkg would still emit
  # 1-2 annual rows Julia omits.
  if (nrow(flashiness_by_year) < 3) return(empty_stats("flashinessRB"))

  result <- generate_stats(flashiness_by_year, value_cols = "flashinessRB",
                           year_col = "water_year",
                           trend_completeness = trend_completeness,
                           decade_completeness = decade_completeness, min_values_for_stats = min_values_for_stats, changepoint = changepoint, collector = collector)
  result
}

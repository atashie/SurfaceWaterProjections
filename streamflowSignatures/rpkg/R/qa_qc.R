#' Columns that enter the \code{flagged_for_high_na} denominator
#'
#' Every name ending in one of the statistic suffixes
#' (\code{pkg_env$qa_high_na_suffixes}) plus the per-gage scalar registry
#' (\code{pkg_env$qa_high_na_scalars}), in input order; \code{flagged_*} columns
#' excluded. Shared definition with \code{julia/src/qa_qc.jl} and
#' \code{python/streamflow_signatures/qa_qc.py} (config
#' \code{qa_qc.high_na_denominator}).
#'
#' @param cols Character vector of column / key names.
#' @return Character vector (subset of \code{cols}).
#' @export
high_na_denominator_columns <- function(cols) {
  cols <- as.character(cols)
  is_flag <- startsWith(cols, "flagged_")
  is_scalar <- cols %in% pkg_env$qa_high_na_scalars
  is_stat <- vapply(cols, function(nm) any(endsWith(nm, pkg_env$qa_high_na_suffixes)),
                    logical(1), USE.NAMES = FALSE)
  cols[!is_flag & (is_scalar | is_stat)]
}

# Internal: the 11 range/consistency flags for ONE gage given as a named list.
.qa_row_flags <- function(sigs) {
  # Initialize all 12 flags to FALSE (always present, matching Python/Julia)
  flag_names <- c(
    "flagged_for_qann_range",
    "flagged_for_bfi_eckhardt_range",
    "flagged_for_bfi_lynehollick_range",
    "flagged_for_flashiness_range",
    "flagged_for_tqmean_range",
    "flagged_for_d50_range",
    "flagged_for_elasticity_range",
    "flagged_for_runoff_ratio_range",
    "flagged_for_seasonal_sum",
    "flagged_for_percentile_order",
    "flagged_for_timing_order",
    "flagged_for_high_na"
  )
  flags <- setNames(as.list(rep(FALSE, length(flag_names))), flag_names)

  get_val <- function(name) {
    v <- sigs[[name]]
    if (is.null(v)) NA_real_ else v
  }

  # 1. Qann range
  qann_mean <- get_val("Qann_mean")
  r <- pkg_env$qa_qann_range
  flags$flagged_for_qann_range <- !is.na(qann_mean) && (qann_mean < r[1] || qann_mean > r[2])

  # 2. BFI Eckhardt range
  bfi_e <- get_val("BFI_Eckhardt_mean")
  r <- pkg_env$qa_bfi_range
  flags$flagged_for_bfi_eckhardt_range <- !is.na(bfi_e) && (bfi_e < r[1] || bfi_e > r[2])

  # 3. BFI LyneHollick range
  bfi_l <- get_val("BFI_LyneHollick_mean")
  flags$flagged_for_bfi_lynehollick_range <- !is.na(bfi_l) && (bfi_l < r[1] || bfi_l > r[2])

  # 4. Flashiness range
  flash <- get_val("flashinessRB_mean")
  r <- pkg_env$qa_flashiness_range
  flags$flagged_for_flashiness_range <- !is.na(flash) && (flash < r[1] || flash > r[2])

  # 5. TQmean range
  tq <- get_val("TQmean_mean")
  r <- pkg_env$qa_tqmean_range
  flags$flagged_for_tqmean_range <- !is.na(tq) && (tq < r[1] || tq > r[2])

  # 6. D50 range
  d50 <- get_val("D50_day_mean")
  r <- pkg_env$qa_d50_range
  flags$flagged_for_d50_range <- !is.na(d50) && (d50 < r[1] || d50 > r[2])

  # 7. Elasticity range
  elast <- get_val("elasticity_static")
  r <- pkg_env$qa_elasticity_range
  flags$flagged_for_elasticity_range <- !is.na(elast) && (elast < r[1] || elast > r[2])

  # 8. Runoff ratio range
  rr <- get_val("annual_runoff_ratio_mean")
  r <- pkg_env$qa_runoff_ratio_range
  flags$flagged_for_runoff_ratio_range <- !is.na(rr) && (rr < r[1] || rr > r[2])

  # 9. Seasonal sum tolerance: |Qwin+Qspr+Qsum+Qfal - Qann| / Qann
  qann <- get_val("Qann_mean")
  qwin <- get_val("Qwin_mean")
  qspr <- get_val("Qspr_mean")
  qsum <- get_val("Qsum_mean")
  qfal <- get_val("Qfal_mean")
  if (!is.na(qann) && qann > 0 && !is.na(qwin) && !is.na(qspr) && !is.na(qsum) && !is.na(qfal)) {
    seasonal_sum <- qwin + qspr + qsum + qfal
    flags$flagged_for_seasonal_sum <- abs(seasonal_sum - qann) / qann > pkg_env$qa_seasonal_sum_tolerance
  }

  # 10. Percentile ordering: Q5 < Q25 < Q50 < Q75 < Q95
  q5  <- get_val("Q5_mean")
  q25 <- get_val("Q25_mean")
  q50 <- get_val("Q50_mean")
  q75 <- get_val("Q75_mean")
  q95 <- get_val("Q95_mean")
  pcts <- c(q5, q25, q50, q75, q95)
  if (all(!is.na(pcts))) {
    flags$flagged_for_percentile_order <- !all(diff(pcts) >= 0)
  }

  # 11. Timing order: D5 < D50 < D95
  d5  <- get_val("D5_day_mean")
  d50_t <- get_val("D50_day_mean")
  d95 <- get_val("D95_day_mean")
  if (!is.na(d5) && !is.na(d50_t) && !is.na(d95)) {
    flags$flagged_for_timing_order <- (d5 > d50_t) || (d50_t > d95)
  }

  flags
}

#' Compute QA/QC flags for signature results
#'
#' Checks ranges for common signatures and returns the 12 \code{flagged_for_*}
#' flags matching the Python/Julia convention.
#'
#' Two input shapes:
#' \itemize{
#'   \item a named \strong{list} (one gage, the output of
#'     \code{calculate_all_signatures}) — returns a named list of 12 logicals;
#'     \code{flagged_for_high_na} is evaluated over the signature keys the gage
#'     EMITTED (families it never computed are absent, not NA);
#'   \item a \strong{data.frame} (the assembled output table, one row per gage) —
#'     returns a data.frame of 12 logical columns; \code{flagged_for_high_na} is
#'     evaluated over the signature columns PRESENT in the table, so families
#'     NA-filled by the table union DO count. This is the production path (the
#'     benchmark runner) and matches the Julia and Python runners exactly.
#' }
#' The denominator is \code{high_na_denominator_columns(names)}: the 16 statistic
#' suffixes plus the per-gage scalar registry. Metadata, unregistered diagnostics
#' and flag columns never enter it (before 2026-09-04 rpkg counted only the
#' \code{_mean}/\code{_median} keys the gage emitted).
#'
#' @param sigs Named list (one gage) or data.frame (one row per gage).
#' @return Named list of 12 logicals, or a data.frame with 12 logical columns.
#' @export
compute_qa_flags <- function(sigs) {
  if (is.data.frame(sigs)) {
    df <- as.data.frame(sigs, stringsAsFactors = FALSE)
    rows <- lapply(seq_len(nrow(df)), function(i) .qa_row_flags(as.list(df[i, , drop = FALSE])))
    out <- as.data.frame(do.call(rbind, lapply(rows, function(r) unlist(r))),
                         stringsAsFactors = FALSE)
    out[] <- lapply(out, as.logical)
    if (nrow(df) == 0) {
      out <- as.data.frame(setNames(replicate(12, logical(0), simplify = FALSE),
                                    names(.qa_row_flags(list()))))
    }
    cols <- high_na_denominator_columns(names(df))
    if (length(cols) > 0) {
      na_frac <- rowMeans(is.na(df[, cols, drop = FALSE]))
      out$flagged_for_high_na <- na_frac > pkg_env$qa_max_na_fraction
    }
    rownames(out) <- NULL
    return(out)
  }

  flags <- .qa_row_flags(sigs)
  cols <- high_na_denominator_columns(names(sigs))
  if (length(cols) > 0) {
    vals <- unlist(lapply(cols, function(nm) { v <- sigs[[nm]]; if (is.null(v)) NA else v[1] }))
    flags$flagged_for_high_na <- sum(is.na(vals)) / length(vals) > pkg_env$qa_max_na_fraction
  }
  flags
}

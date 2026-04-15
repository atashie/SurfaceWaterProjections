#' Compute QA/QC flags for a named list of signature results
#'
#' Checks ranges for common signatures and returns a named logical list
#' of 12 flags matching the Python/Julia \code{flagged_for_*} convention.
#'
#' @param sigs Named list of signature values (output of calculate_all_signatures).
#' @return Named list of 12 logical flags (TRUE = potential issue).
#' @export
compute_qa_flags <- function(sigs) {
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

  # 12. NA fraction
  sig_vals <- unlist(sigs[grep("_(mean|median)$", names(sigs))])
  if (length(sig_vals) > 0) {
    flags$flagged_for_high_na <- sum(is.na(sig_vals)) / length(sig_vals) > pkg_env$qa_max_na_fraction
  }

  flags
}

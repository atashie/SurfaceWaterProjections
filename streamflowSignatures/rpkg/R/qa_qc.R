#' Compute QA/QC flags for a named list of signature results
#'
#' Checks ranges for common signatures and returns a named logical list
#' of flags.
#'
#' @param sigs Named list of signature values (output of calculate_all_signatures).
#' @return Named list of logical flags (TRUE = potential issue).
#' @export
compute_qa_flags <- function(sigs) {
  flags <- list()

  get_val <- function(name) {
    v <- sigs[[name]]
    if (is.null(v)) NA_real_ else v
  }

  # Qann range
  qann_mean <- get_val("Qann_mean")
  r <- pkg_env$qa_qann_range
  flags$qann_out_of_range <- !is.na(qann_mean) && (qann_mean < r[1] || qann_mean > r[2])

  # BFI range
  for (bfi in c("BFI_Eckhardt_mean", "BFI_LyneHollick_mean")) {
    v <- get_val(bfi)
    r <- pkg_env$qa_bfi_range
    flags[[paste0(bfi, "_out_of_range")]] <- !is.na(v) && (v < r[1] || v > r[2])
  }

  # BFI consistency: Eckhardt < LyneHollick
  bfi_e <- get_val("BFI_Eckhardt_mean")
  bfi_l <- get_val("BFI_LyneHollick_mean")
  flags$bfi_inconsistent <- !is.na(bfi_e) && !is.na(bfi_l) && bfi_e > bfi_l

  # Flashiness range
  flash <- get_val("flashinessRB_mean")
  r <- pkg_env$qa_flashiness_range
  flags$flashiness_out_of_range <- !is.na(flash) && (flash < r[1] || flash > r[2])

  # TQmean range
  tq <- get_val("TQmean_mean")
  r <- pkg_env$qa_tqmean_range
  flags$tqmean_out_of_range <- !is.na(tq) && (tq < r[1] || tq > r[2])

  # D50 range
  d50 <- get_val("D50_day_mean")
  r <- pkg_env$qa_d50_range
  flags$d50_out_of_range <- !is.na(d50) && (d50 < r[1] || d50 > r[2])

  # Elasticity range
  elast <- get_val("elasticity_static")
  r <- pkg_env$qa_elasticity_range
  flags$elasticity_out_of_range <- !is.na(elast) && (elast < r[1] || elast > r[2])

  # Runoff ratio range
  rr <- get_val("annual_runoff_ratio_mean")
  r <- pkg_env$qa_runoff_ratio_range
  flags$runoff_ratio_out_of_range <- !is.na(rr) && (rr < r[1] || rr > r[2])

  # Seasonal sum tolerance: |Qwin+Qspr+Qsum+Qfal - Qann| / Qann
  qann <- get_val("Qann_mean")
  qwin <- get_val("Qwin_mean")
  qspr <- get_val("Qspr_mean")
  qsum <- get_val("Qsum_mean")
  qfal <- get_val("Qfal_mean")
  if (!is.na(qann) && qann > 0 && !is.na(qwin) && !is.na(qspr) && !is.na(qsum) && !is.na(qfal)) {
    seasonal_sum <- qwin + qspr + qsum + qfal
    flags$seasonal_sum_mismatch <- abs(seasonal_sum - qann) / qann > pkg_env$qa_seasonal_sum_tolerance
  }

  # Percentile ordering: Q5 < Q25 < Q50 < Q75 < Q95
  q5  <- get_val("Q5_mean")
  q25 <- get_val("Q25_mean")
  q50 <- get_val("Q50_mean")
  q75 <- get_val("Q75_mean")
  q95 <- get_val("Q95_mean")
  pcts <- c(q5, q25, q50, q75, q95)
  if (all(!is.na(pcts))) {
    flags$percentile_order_violated <- !all(diff(pcts) >= 0)
  }

  # NA fraction
  sig_vals <- unlist(sigs[grep("_(mean|median)$", names(sigs))])
  if (length(sig_vals) > 0) {
    flags$high_na_fraction <- sum(is.na(sig_vals)) / length(sig_vals) > pkg_env$qa_max_na_fraction
  }

  flags
}

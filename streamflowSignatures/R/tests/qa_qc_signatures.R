################################################################################
# QA/QC Validation for Streamflow Signatures
#
# Tests for hydrological validity, physical constraints, and statistical sanity
# Output: data_out/qa_qc_report.txt
################################################################################

# Clear environment
rm(list = ls())

# Set working directory
main_dir <- getwd()
if (!file.exists(file.path(main_dir, "config.R"))) {
  main_dir <- dirname(main_dir)
  if (!file.exists(file.path(main_dir, "config.R"))) {
    main_dir <- dirname(main_dir)
  }
  setwd(main_dir)
}

# Load configuration and libraries
source("config.R")
library(data.table)

cat("========== STREAMFLOW SIGNATURES QA/QC VALIDATION ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Define paths
input_path <- "data_out/streamflow_signatures_full_JAN2026.csv"
report_path <- "data_out/qa_qc_report.txt"

# Check if input exists
if (!file.exists(input_path)) {
  stop("Signature file not found: ", input_path, "\n  Run run_full_processing.R first")
}

# Load data
cat("Loading signature data...\n")
signatures <- fread(input_path)
cat("  Loaded", nrow(signatures), "gages with", ncol(signatures), "columns\n\n")

# Initialize report
report <- c(
  "========================================",
  "STREAMFLOW SIGNATURES QA/QC REPORT",
  paste("Generated:", Sys.time()),
  paste("Input file:", input_path),
  paste("Gages:", nrow(signatures)),
  paste("Columns:", ncol(signatures)),
  "========================================",
  ""
)

# Helper function to add to report
add_to_report <- function(...) {
  msg <- paste(...)
  report <<- c(report, msg)
  cat(msg, "\n")
}

# Test counter
tests_passed <- 0
tests_failed <- 0
tests_warning <- 0

run_test <- function(test_name, condition, warning_only = FALSE) {
  if (condition) {
    add_to_report("[PASS]", test_name)
    tests_passed <<- tests_passed + 1
  } else if (warning_only) {
    add_to_report("[WARN]", test_name)
    tests_warning <<- tests_warning + 1
  } else {
    add_to_report("[FAIL]", test_name)
    tests_failed <<- tests_failed + 1
  }
}

################################################################################
# 1. RANGE VALIDATION TESTS
################################################################################
add_to_report("")
add_to_report("==================== RANGE VALIDATION ====================")

# Define expected ranges for key signatures
range_tests <- list(
  list(col = "Qann_mean", min = 0, max = 2000, name = "Annual mean runoff (mm/yr)"),
  list(col = "BFI_Eckhardt_mean", min = 0, max = 1, name = "BFI Eckhardt (fraction)"),
  list(col = "BFI_LyneHollick_mean", min = 0, max = 1, name = "BFI Lyne-Hollick (fraction)"),
  list(col = "flashinessRB_mean", min = 0, max = 2, name = "R-B Flashiness Index"),
  list(col = "TQmean_mean", min = 0, max = 100, name = "TQmean (% days above mean)"),
  list(col = "D50_day_mean", min = 1, max = 366, name = "D50 day of year"),
  list(col = "elasticity_static", min = 0.1, max = 5, name = "Streamflow elasticity"),
  list(col = "annual_runoff_ratio_mean", min = 0.01, max = 1.5, name = "Annual runoff ratio")
)

for (test in range_tests) {
  if (test$col %in% names(signatures)) {
    vals <- signatures[[test$col]]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      in_range <- sum(vals >= test$min & vals <= test$max)
      pct_in_range <- in_range / length(vals) * 100
      test_name <- sprintf("%s: %.1f%% in [%.2f, %.1f]",
                          test$name, pct_in_range, test$min, test$max)
      run_test(test_name, pct_in_range >= 95, warning_only = pct_in_range >= 80)

      # Report outliers
      n_low <- sum(vals < test$min)
      n_high <- sum(vals > test$max)
      if (n_low > 0 || n_high > 0) {
        add_to_report(sprintf("       Outliers: %d below min, %d above max", n_low, n_high))
      }
    } else {
      add_to_report("[SKIP]", test$name, "- no valid values")
    }
  } else {
    add_to_report("[SKIP]", test$name, "- column not found")
  }
}

################################################################################
# 2. CONSISTENCY TESTS
################################################################################
add_to_report("")
add_to_report("==================== CONSISTENCY TESTS ====================")

# Test 2.1: Seasonal sum check (Qwin + Qspr + Qsum + Qfal ~ Qann)
# NOTE: Qann and seasonal values are TOTALS (mm/year, mm/season), not daily means
# So seasonal totals should sum to annual total
if (all(c("Qann_mean", "Qwin_mean", "Qspr_mean", "Qsum_mean", "Qfal_mean") %in% names(signatures))) {
  seasonal_sum <- signatures$Qwin_mean + signatures$Qspr_mean +
                  signatures$Qsum_mean + signatures$Qfal_mean
  expected <- signatures$Qann_mean  # Seasonal totals should sum to annual total
  valid_rows <- !is.na(seasonal_sum) & !is.na(expected) & expected > 0
  if (sum(valid_rows) > 0) {
    ratio <- seasonal_sum[valid_rows] / expected[valid_rows]
    pct_close <- sum(ratio > 0.8 & ratio < 1.2) / sum(valid_rows) * 100
    run_test(sprintf("Seasonal sum ~ Qann: %.1f%% within 20%%", pct_close),
             pct_close >= 90, warning_only = pct_close >= 75)
  }
}

# Test 2.2: Percentile ordering (Q05 < Q25 < Q50 < Q75 < Q95)
percentile_cols <- c("Q5_mean", "Q25_mean", "Q50_mean", "Q75_mean", "Q95_mean")
if (all(percentile_cols %in% names(signatures))) {
  valid_rows <- complete.cases(signatures[, ..percentile_cols])
  if (sum(valid_rows) > 0) {
    ordered <- signatures[valid_rows, Q5_mean < Q25_mean &
                          Q25_mean < Q50_mean &
                          Q50_mean < Q75_mean &
                          Q75_mean < Q95_mean]
    pct_ordered <- sum(ordered) / sum(valid_rows) * 100
    run_test(sprintf("Flow percentile ordering: %.1f%% correctly ordered", pct_ordered),
             pct_ordered >= 99, warning_only = pct_ordered >= 95)
  }
}

# Test 2.3: BFI method consistency
if (all(c("BFI_Eckhardt_mean", "BFI_LyneHollick_mean") %in% names(signatures))) {
  valid_rows <- !is.na(signatures$BFI_Eckhardt_mean) & !is.na(signatures$BFI_LyneHollick_mean)
  if (sum(valid_rows) > 10) {
    bfi_cor <- cor(signatures$BFI_Eckhardt_mean[valid_rows],
                   signatures$BFI_LyneHollick_mean[valid_rows])
    run_test(sprintf("BFI method correlation: r = %.3f", bfi_cor),
             bfi_cor >= 0.7, warning_only = bfi_cor >= 0.5)
  }
}

# Test 2.4: Flow timing ordering (D05 < D25 < D50 < D75 < D95)
timing_cols <- c("D5_day_mean", "D25_to_D75_mean", "D50_day_mean", "D95_day_mean")
if ("D5_day_mean" %in% names(signatures) && "D50_day_mean" %in% names(signatures) &&
    "D95_day_mean" %in% names(signatures)) {
  valid_rows <- !is.na(signatures$D5_day_mean) & !is.na(signatures$D50_day_mean) &
                !is.na(signatures$D95_day_mean)
  if (sum(valid_rows) > 0) {
    ordered <- signatures[valid_rows, D5_day_mean < D50_day_mean &
                          D50_day_mean < D95_day_mean]
    pct_ordered <- sum(ordered) / sum(valid_rows) * 100
    run_test(sprintf("Flow timing ordering: %.1f%% correctly ordered", pct_ordered),
             pct_ordered >= 95, warning_only = pct_ordered >= 85)
  }
}

################################################################################
# 3. STATISTICAL SANITY TESTS
################################################################################
add_to_report("")
add_to_report("==================== STATISTICAL SANITY ====================")

# Test 3.1: Trend method agreement (Theil-Sen vs Linear regression)
# Find all _senn_slp and corresponding _linear_slp columns
senn_cols <- grep("_senn_slp$", names(signatures), value = TRUE)
if (length(senn_cols) > 5) {
  same_sign_count <- 0
  total_comparisons <- 0

  for (senn_col in head(senn_cols, 20)) {  # Sample first 20
    linear_col <- gsub("_senn_slp$", "_linear_slp", senn_col)
    if (linear_col %in% names(signatures)) {
      valid_rows <- !is.na(signatures[[senn_col]]) & !is.na(signatures[[linear_col]])
      if (sum(valid_rows) > 0) {
        same_sign <- sum(sign(signatures[[senn_col]][valid_rows]) ==
                         sign(signatures[[linear_col]][valid_rows]))
        same_sign_count <- same_sign_count + same_sign
        total_comparisons <- total_comparisons + sum(valid_rows)
      }
    }
  }

  if (total_comparisons > 0) {
    pct_same_sign <- same_sign_count / total_comparisons * 100
    run_test(sprintf("Theil-Sen vs Linear slope sign agreement: %.1f%%", pct_same_sign),
             pct_same_sign >= 80, warning_only = pct_same_sign >= 70)
  }
}

# Test 3.2: Spearman vs Mann-Kendall correlation
spearman_cols <- grep("_spearman_rho$", names(signatures), value = TRUE)
if (length(spearman_cols) > 5) {
  all_spearman <- c()
  all_mk <- c()

  for (sp_col in head(spearman_cols, 20)) {
    mk_col <- gsub("_spearman_rho$", "_mk_rho", sp_col)
    if (mk_col %in% names(signatures)) {
      valid_rows <- !is.na(signatures[[sp_col]]) & !is.na(signatures[[mk_col]])
      all_spearman <- c(all_spearman, signatures[[sp_col]][valid_rows])
      all_mk <- c(all_mk, signatures[[mk_col]][valid_rows])
    }
  }

  if (length(all_spearman) > 100) {
    rho_cor <- cor(all_spearman, all_mk)
    run_test(sprintf("Spearman vs Mann-Kendall tau correlation: r = %.3f", rho_cor),
             rho_cor >= 0.9, warning_only = rho_cor >= 0.8)
  }
}

# Test 3.3: NA fraction check
sig_cols <- grep("_(mean|median)$", names(signatures), value = TRUE)
if (length(sig_cols) > 0) {
  na_fractions <- sapply(signatures[, ..sig_cols], function(x) sum(is.na(x)) / length(x))
  overall_na <- mean(na_fractions) * 100
  run_test(sprintf("Overall NA fraction in signatures: %.1f%%", overall_na),
           overall_na < 20, warning_only = overall_na < 30)

  # Flag gages with high NA
  per_gage_na <- rowSums(is.na(signatures[, ..sig_cols])) / length(sig_cols) * 100
  high_na_gages <- sum(per_gage_na > 30)
  add_to_report(sprintf("       Gages with >30%% NA signatures: %d (%.1f%%)",
                        high_na_gages, high_na_gages / nrow(signatures) * 100))
}

################################################################################
# 4. CLIMATE SIGNATURE TESTS
################################################################################
add_to_report("")
add_to_report("==================== CLIMATE SIGNATURES ====================")

# Test 4.1: Elasticity distribution
if ("elasticity_static" %in% names(signatures)) {
  e_vals <- signatures$elasticity_static[!is.na(signatures$elasticity_static)]
  if (length(e_vals) > 0) {
    e_median <- median(e_vals)
    e_iqr <- IQR(e_vals)
    run_test(sprintf("Elasticity distribution: median=%.2f, IQR=%.2f", e_median, e_iqr),
             e_median > 0.5 && e_median < 3 && e_iqr < 3,
             warning_only = TRUE)
    add_to_report(sprintf("       Range: %.2f to %.2f", min(e_vals), max(e_vals)))
  }
}

# Test 4.2: Storage validity
if ("avg_storage_mean" %in% names(signatures)) {
  s_vals <- signatures$avg_storage_mean[!is.na(signatures$avg_storage_mean)]
  if (length(s_vals) > 0) {
    pct_positive <- sum(s_vals > 0) / length(s_vals) * 100
    pct_reasonable <- sum(s_vals > 0 & s_vals < 1000) / length(s_vals) * 100
    run_test(sprintf("Storage: %.1f%% positive, %.1f%% < 1000mm", pct_positive, pct_reasonable),
             pct_positive >= 80 && pct_reasonable >= 90,
             warning_only = pct_positive >= 60)
  }
}

# Test 4.3: Climate signature coverage
climate_cols <- c("elasticity_static", "qp_slope_sd_mean", "qp_bimodality_mean",
                  "avg_storage_mean", "annual_runoff_ratio_mean")
existing_climate <- intersect(climate_cols, names(signatures))
if (length(existing_climate) > 0) {
  coverages <- sapply(existing_climate, function(col) {
    sum(!is.na(signatures[[col]])) / nrow(signatures) * 100
  })
  avg_coverage <- mean(coverages)
  add_to_report(sprintf("Climate signature coverage: %.1f%% average", avg_coverage))
  for (i in seq_along(existing_climate)) {
    add_to_report(sprintf("       %s: %.1f%%", existing_climate[i], coverages[i]))
  }
}

################################################################################
# 5. HYDROLOGICAL RELATIONSHIP TESTS
################################################################################
add_to_report("")
add_to_report("==================== HYDROLOGICAL RELATIONSHIPS ====================")

# Test 5.1: BFI vs Flashiness (expect negative correlation)
if (all(c("BFI_Eckhardt_mean", "flashinessRB_mean") %in% names(signatures))) {
  valid_rows <- !is.na(signatures$BFI_Eckhardt_mean) & !is.na(signatures$flashinessRB_mean)
  if (sum(valid_rows) > 10) {
    bfi_flash_cor <- cor(signatures$BFI_Eckhardt_mean[valid_rows],
                         signatures$flashinessRB_mean[valid_rows])
    run_test(sprintf("BFI vs Flashiness correlation: r = %.3f (expect negative)", bfi_flash_cor),
             bfi_flash_cor < 0, warning_only = bfi_flash_cor < 0.2)
  }
}

# Test 5.2: Qann vs basin_area (expect positive correlation in log space)
if (all(c("Qann_mean", "basin_area") %in% names(signatures))) {
  valid_rows <- !is.na(signatures$Qann_mean) & !is.na(signatures$basin_area) &
                signatures$Qann_mean > 0 & signatures$basin_area > 0
  if (sum(valid_rows) > 10) {
    # Note: Qann in mm/day is area-normalized, so may not correlate directly
    # But Q (total) = Qann * area, so test makes sense for volumetric context
    add_to_report("[INFO] Qann is already area-normalized (mm/day) - correlation with area not expected")
  }
}

################################################################################
# 6. CREATE PER-GAGE FLAG COLUMNS
################################################################################
add_to_report("")
add_to_report("==================== FLAGGING INDIVIDUAL GAGES ====================")

# Initialize all flag columns as FALSE
signatures[, flagged_for_qann_range := FALSE]
signatures[, flagged_for_bfi_eckhardt_range := FALSE]
signatures[, flagged_for_bfi_lynehollick_range := FALSE]
signatures[, flagged_for_flashiness_range := FALSE]
signatures[, flagged_for_tqmean_range := FALSE]
signatures[, flagged_for_d50_range := FALSE]
signatures[, flagged_for_elasticity_range := FALSE]
signatures[, flagged_for_runoff_ratio_range := FALSE]
signatures[, flagged_for_seasonal_sum := FALSE]
signatures[, flagged_for_percentile_order := FALSE]
signatures[, flagged_for_timing_order := FALSE]
signatures[, flagged_for_high_na := FALSE]

# Apply range flags
if ("Qann_mean" %in% names(signatures)) {
  signatures[!is.na(Qann_mean) & (Qann_mean < 0 | Qann_mean > 2000),
             flagged_for_qann_range := TRUE]
}

if ("BFI_Eckhardt_mean" %in% names(signatures)) {
  signatures[!is.na(BFI_Eckhardt_mean) & (BFI_Eckhardt_mean < 0 | BFI_Eckhardt_mean > 1),
             flagged_for_bfi_eckhardt_range := TRUE]
}

if ("BFI_LyneHollick_mean" %in% names(signatures)) {
  signatures[!is.na(BFI_LyneHollick_mean) & (BFI_LyneHollick_mean < 0 | BFI_LyneHollick_mean > 1),
             flagged_for_bfi_lynehollick_range := TRUE]
}

if ("flashinessRB_mean" %in% names(signatures)) {
  signatures[!is.na(flashinessRB_mean) & (flashinessRB_mean < 0 | flashinessRB_mean > 2),
             flagged_for_flashiness_range := TRUE]
}

if ("TQmean_mean" %in% names(signatures)) {
  signatures[!is.na(TQmean_mean) & (TQmean_mean < 0 | TQmean_mean > 100),
             flagged_for_tqmean_range := TRUE]
}

if ("D50_day_mean" %in% names(signatures)) {
  signatures[!is.na(D50_day_mean) & (D50_day_mean < 1 | D50_day_mean > 366),
             flagged_for_d50_range := TRUE]
}

if ("elasticity_static" %in% names(signatures)) {
  signatures[!is.na(elasticity_static) & (elasticity_static < 0.1 | elasticity_static > 5),
             flagged_for_elasticity_range := TRUE]
}

if ("annual_runoff_ratio_mean" %in% names(signatures)) {
  signatures[!is.na(annual_runoff_ratio_mean) & (annual_runoff_ratio_mean < 0.01 | annual_runoff_ratio_mean > 1.5),
             flagged_for_runoff_ratio_range := TRUE]
}

# Seasonal sum flag
if (all(c("Qann_mean", "Qwin_mean", "Qspr_mean", "Qsum_mean", "Qfal_mean") %in% names(signatures))) {
  signatures[, seasonal_sum_temp := Qwin_mean + Qspr_mean + Qsum_mean + Qfal_mean]
  signatures[!is.na(seasonal_sum_temp) & !is.na(Qann_mean) & Qann_mean > 0,
             flagged_for_seasonal_sum := abs(seasonal_sum_temp / Qann_mean - 1) > 0.2]
  signatures[, seasonal_sum_temp := NULL]
}

# Percentile ordering flag
if (all(c("Q5_mean", "Q25_mean", "Q50_mean", "Q75_mean", "Q95_mean") %in% names(signatures))) {
  signatures[complete.cases(signatures[, .(Q5_mean, Q25_mean, Q50_mean, Q75_mean, Q95_mean)]),
             flagged_for_percentile_order := !(Q5_mean < Q25_mean & Q25_mean < Q50_mean &
                                                Q50_mean < Q75_mean & Q75_mean < Q95_mean)]
}

# Timing ordering flag
if (all(c("D5_day_mean", "D50_day_mean", "D95_day_mean") %in% names(signatures))) {
  signatures[!is.na(D5_day_mean) & !is.na(D50_day_mean) & !is.na(D95_day_mean),
             flagged_for_timing_order := !(D5_day_mean < D50_day_mean & D50_day_mean < D95_day_mean)]
}

# High NA flag
sig_cols <- grep("_(mean|median)$", names(signatures), value = TRUE)
sig_cols <- setdiff(sig_cols, grep("^flagged_", sig_cols, value = TRUE))  # Exclude flag columns
if (length(sig_cols) > 0) {
  signatures[, na_frac := rowSums(is.na(.SD)) / length(sig_cols), .SDcols = sig_cols]
  signatures[na_frac > 0.3, flagged_for_high_na := TRUE]
  signatures[, na_frac := NULL]
}

# Report flag counts
flag_cols <- grep("^flagged_for_", names(signatures), value = TRUE)
add_to_report(sprintf("Created %d flag columns", length(flag_cols)))
for (fc in flag_cols) {
  n_flagged <- sum(signatures[[fc]], na.rm = TRUE)
  if (n_flagged > 0) {
    add_to_report(sprintf("  %s: %d gages (%.1f%%)", fc, n_flagged, n_flagged/nrow(signatures)*100))
  }
}

# Save updated signatures with flags
output_path <- "data_out/streamflow_signatures_full_JAN2026.csv"
fwrite(signatures, output_path)
add_to_report(sprintf("\nUpdated signatures file saved with flag columns: %s", output_path))

################################################################################
# SUMMARY
################################################################################
add_to_report("")
add_to_report("==================== SUMMARY ====================")
add_to_report(sprintf("Tests passed: %d", tests_passed))
add_to_report(sprintf("Tests warning: %d", tests_warning))
add_to_report(sprintf("Tests failed: %d", tests_failed))
add_to_report(sprintf("Total tests: %d", tests_passed + tests_warning + tests_failed))
add_to_report("")

if (tests_failed == 0) {
  add_to_report("OVERALL STATUS: PASS (all critical tests passed)")
} else if (tests_failed <= 2) {
  add_to_report("OVERALL STATUS: WARNING (minor issues detected)")
} else {
  add_to_report("OVERALL STATUS: REVIEW NEEDED (multiple issues detected)")
}

# Save report
writeLines(report, report_path)
cat("\n\nReport saved to:", report_path, "\n")

cat("\n========== QA/QC COMPLETE ==========\n")

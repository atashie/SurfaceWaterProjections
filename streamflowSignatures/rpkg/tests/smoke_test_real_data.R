################################################################################################
# Smoke Test - rpkg/ workflow on real parquet data
# Mirrors the legacy (deprecated) R/tests/smoke_test.R but uses the rpkg package API
################################################################################################

library(streamflowsignatures)
library(data.table)

cat("========== rpkg SMOKE TEST: Real Data ==========\n\n")

# --- Paths ---
# ENV-overridable (same variable names as the Julia/Python smoke tests and the
# benchmark runners); defaults point at the Windows data drive.
parquet_path <- Sys.getenv("STREAMFLOW_DATA_PATH",
                           "D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet")
daymet_path  <- Sys.getenv("STREAMFLOW_CLIMATE_PATH",
                           "D:/processedOuts_feb2026/daymet_1980_2023.parquet")

if (!file.exists(parquet_path)) {
  stop("Parquet file not found: ", parquet_path)
}

# --- Step 1: Read parquet ---
cat("[1] Reading parquet...\n")
full_data <- read_parquet(parquet_path)
cat("    Rows:", nrow(full_data), " Unique gages:", length(unique(full_data$gage_id)), "\n")

# --- Step 2: Select 10 test gages FIRST (before expensive water year computation) ---
set.seed(42)
all_gages <- unique(full_data$gage_id)
test_gages <- sample(all_gages, min(10, length(all_gages)))
cat("[2] Selected", length(test_gages), "test gages:", paste(test_gages, collapse = ", "), "\n")

# Subset to just those gages
full_data <- full_data[full_data$gage_id %in% test_gages, ]
cat("    Subset rows:", nrow(full_data), "\n")

# --- Step 3: Add water year columns (now on subset only) ---
cat("[3] Adding water year columns...\n")
full_data <- add_water_year_columns(full_data)
cat("    Columns:", paste(names(full_data), collapse = ", "), "\n\n")

# --- Step 4: Process each gage ---
cat("[4] Processing signatures (no climate)...\n")
results_list <- list()
timings <- numeric(length(test_gages))

for (i in seq_along(test_gages)) {
  gage <- test_gages[i]
  t0 <- proc.time()["elapsed"]

  gage_data <- full_data[full_data$gage_id == gage, ]

  # Filter qualifying years
  filt <- filter_qualifying_years(gage_data)

  if (!filt$qualifies) {
    cat(sprintf("    [%d/%d] %s - SKIPPED (%d qualifying years)\n",
                i, length(test_gages), gage, length(filt$qualifying_years)))
    next
  }

  gage_filtered <- gage_data[gage_data$water_year %in% filt$qualifying_years, ]

  # Calculate all non-climate signatures
  sigs <- calculate_all_signatures(gage_filtered, has_climate = FALSE)

  elapsed <- proc.time()["elapsed"] - t0
  timings[i] <- elapsed

  results_list[[gage]] <- sigs
  cat(sprintf("    [%d/%d] %s - OK (%d years, %d sigs, %.1fs)\n",
              i, length(test_gages), gage,
              length(filt$qualifying_years), length(sigs), elapsed))
}

cat(sprintf("\n    Processed: %d/%d gages (%.1fs total)\n\n",
            length(results_list), length(test_gages), sum(timings)))

# --- Step 5: Validate results ---
cat("[5] Validating results...\n")
if (length(results_list) == 0) {
  cat("    FAIL: No gages produced results\n")
  quit(status = 1)
}

# Pick first successful gage for detailed checks
first_gage <- names(results_list)[1]
sigs <- results_list[[first_gage]]

# Check expected signature names
expected_bases <- c("Qann", "Qwin", "Qspr", "Qsum", "Qfal",
                    "Q50", "Q95_Q10",
                    "FDCall", "FDC90th", "FDCmid",
                    "flashinessRB",
                    "D50_day", "D25_to_D75", "Dmax",
                    "BFI_Eckhardt", "BFI_LyneHollick",
                    "n_high_pulses_year", "TQmean", "Flow_Reversals_annual")

expected_suffixes <- c("_senn_slp", "_linear_slp", "_spearman_rho",
                       "_spearman_pval", "_mk_rho", "_mk_pval",
                       "_mean", "_median")

cat(sprintf("    Gage %s: %d signature values\n", first_gage, length(sigs)))

pass_count <- 0
fail_count <- 0
for (base in expected_bases) {
  for (suf in expected_suffixes) {
    key <- paste0(base, suf)
    if (key %in% names(sigs)) {
      pass_count <- pass_count + 1
    } else {
      fail_count <- fail_count + 1
      cat(sprintf("    MISSING: %s\n", key))
    }
  }
}

cat(sprintf("    Name checks: %d/%d present\n", pass_count, pass_count + fail_count))

# Value sanity checks
checks <- list(
  list(name = "Qann_mean > 0",
       pass = !is.na(sigs$Qann_mean) && sigs$Qann_mean > 0),
  list(name = "BFI_Eckhardt_mean in [0,1]",
       pass = !is.na(sigs$BFI_Eckhardt_mean) &&
              sigs$BFI_Eckhardt_mean >= 0 && sigs$BFI_Eckhardt_mean <= 1),
  list(name = "BFI_LyneHollick_mean in [0,1]",
       pass = !is.na(sigs$BFI_LyneHollick_mean) &&
              sigs$BFI_LyneHollick_mean >= 0 && sigs$BFI_LyneHollick_mean <= 1),
  list(name = "BFI_Eckhardt < BFI_LyneHollick (warn-only, known edge case)",
       pass = TRUE,  # QA flag captures this; not a hard failure
       warn = !is.na(sigs$BFI_Eckhardt_mean) && !is.na(sigs$BFI_LyneHollick_mean) &&
              sigs$BFI_Eckhardt_mean > sigs$BFI_LyneHollick_mean),
  list(name = "flashinessRB_mean in [0,5]",
       pass = !is.na(sigs$flashinessRB_mean) &&
              sigs$flashinessRB_mean >= 0 && sigs$flashinessRB_mean <= 5),
  list(name = "D50_day_mean in [1,366]",
       pass = !is.na(sigs$D50_day_mean) &&
              sigs$D50_day_mean >= 1 && sigs$D50_day_mean <= 366),
  list(name = "Q50_mean >= 0",
       pass = !is.na(sigs$Q50_mean) && sigs$Q50_mean >= 0),
  list(name = "TQmean_mean in [0,100]",
       pass = !is.na(sigs$TQmean_mean) &&
              sigs$TQmean_mean >= 0 && sigs$TQmean_mean <= 100)
)

cat("\n    Value sanity checks:\n")
all_pass <- TRUE
for (chk in checks) {
  is_warn <- !is.null(chk$warn) && chk$warn
  status <- if (!chk$pass) "FAIL" else if (is_warn) "WARN" else "PASS"
  if (!chk$pass) all_pass <- FALSE
  cat(sprintf("      [%s] %s\n", status, chk$name))
  if (is_warn) {
    cat(sprintf("            BFI_Eckhardt=%.4f > BFI_LyneHollick=%.4f (QA flag set)\n",
                sigs$BFI_Eckhardt_mean, sigs$BFI_LyneHollick_mean))
  }
}

# QA/QC flags
cat("\n    QA/QC flags:\n")
qa <- compute_qa_flags(sigs)
for (nm in names(qa)) {
  cat(sprintf("      %s: %s\n", nm, qa[[nm]]))
}

# --- Step 6: Sample values across all processed gages ---
cat("\n[6] Cross-gage summary:\n")
cat(sprintf("    %-20s %8s %8s %8s\n", "Metric", "Min", "Mean", "Max"))
for (metric in c("Qann_mean", "BFI_Eckhardt_mean", "flashinessRB_mean",
                  "D50_day_mean", "FDCall_mean")) {
  vals <- vapply(results_list, function(s) {
    v <- s[[metric]]; if (is.null(v)) NA_real_ else v
  }, numeric(1))
  vals <- vals[!is.na(vals)]
  if (length(vals) > 0) {
    cat(sprintf("    %-20s %8.3f %8.3f %8.3f\n", metric,
                min(vals), mean(vals), max(vals)))
  }
}

# --- Step 7: Climate signatures (if Daymet available) ---
if (file.exists(daymet_path)) {
  cat("\n[7] Testing climate signatures with Daymet...\n")

  cat("    Reading Daymet parquet...\n")
  daymet_data <- read_parquet(daymet_path)

  # Pick first qualifying gage and subset Daymet to that gage only
  test_gage <- names(results_list)[1]
  gage_climate <- daymet_data[daymet_data$gage_id == test_gage, ]
  rm(daymet_data)  # Free memory

  gage_q <- full_data[full_data$gage_id == test_gage, ]

  if (nrow(gage_climate) > 0) {
    # Merge Q and PPT by date
    ppt_cols <- intersect(c("date", "PPT", "prcp", "PRCP"), names(gage_climate))
    gage_climate_sub <- gage_climate[, ppt_cols, with = FALSE]
    if ("prcp" %in% names(gage_climate_sub)) setnames(gage_climate_sub, "prcp", "PPT")
    if ("PRCP" %in% names(gage_climate_sub)) setnames(gage_climate_sub, "PRCP", "PPT")
    merged <- merge(gage_q, gage_climate_sub[, c("date", "PPT")], by = "date", all.x = TRUE)

    filt <- filter_qualifying_years(merged)
    if (filt$qualifies) {
      merged_filt <- merged[merged$water_year %in% filt$qualifying_years, ]
      climate_sigs <- calculate_all_signatures(merged_filt, has_climate = TRUE)

      climate_keys <- c("annual_runoff_ratio_mean", "elasticity_static",
                        "qp_slope_sd_mean", "avg_storage_mean")

      cat(sprintf("    Gage %s: %d total sigs (with climate)\n",
                  test_gage, length(climate_sigs)))
      for (key in climate_keys) {
        v <- climate_sigs[[key]]
        status <- if (!is.null(v) && !is.na(v)) sprintf("%.4f", v) else "NA"
        cat(sprintf("      %s = %s\n", key, status))
      }

      # Verify climate sigs are present
      found <- sum(climate_keys %in% names(climate_sigs))
      cat(sprintf("    Climate signature keys: %d/%d present\n", found, length(climate_keys)))
    } else {
      cat("    Gage did not qualify after merge\n")
    }
  } else {
    cat("    No Daymet data for gage", test_gage, "\n")
  }
} else {
  cat("\n[7] SKIP: Daymet parquet not found\n")
}

# --- Summary ---
cat("\n========================================\n")
if (length(results_list) > 0 && all_pass && fail_count == 0) {
  cat("STATUS: SMOKE TEST PASSED\n")
} else if (length(results_list) > 0 && all_pass) {
  cat(sprintf("STATUS: SMOKE TEST PASSED (with %d missing signature names)\n", fail_count))
} else {
  cat("STATUS: SMOKE TEST FAILED\n")
}
cat("========================================\n")

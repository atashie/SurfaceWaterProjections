################################################################################################
# Smoke Test - Uses existing parquet data with a subset of gages
################################################################################################

# Paths - auto-detect based on script location
main_dir <- getwd()
if (!file.exists(file.path(main_dir, "config.R"))) {
  # Try parent if running from subdirectory (e.g., R/)
  main_dir <- dirname(main_dir)
  if (!file.exists(file.path(main_dir, "config.R"))) {
    # Try grandparent (e.g., R/tests/)
    main_dir <- dirname(main_dir)
  }
}
metadata_dir <- "D:"

# Load config and functions
source(file.path(main_dir, "config.R"))
source(file.path(main_dir, "R", "helperFunctions.R"))

# Set log level
set_log_level("INFO")

# Input/output paths
parquet_path <- file.path(metadata_dir, "processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet")
metadata_path <- file.path(metadata_dir, "processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv")
output_path <- file.path(main_dir, "test_output/smoke_test_signatures.csv")

# Create output directory
dir.create(file.path(main_dir, "test_output"), showWarnings = FALSE)

cat("========== SMOKE TEST: Subset Signature Extraction ==========\n\n")

# Read parquet and select subset
log_info("Reading parquet file...")
full_data <- arrow::read_parquet(parquet_path)
full_data <- as.data.table(full_data)

log_info("Total rows in parquet:", nrow(full_data))
log_info("Total unique gages:", length(unique(full_data$gage_id)))

# Select 10 random gages for smoke test
set.seed(42)
all_gages <- unique(full_data$gage_id)
test_gages <- sample(all_gages, min(10, length(all_gages)))

log_info("Selected", length(test_gages), "gages for smoke test")
cat("Test gages:", paste(test_gages, collapse = ", "), "\n\n")

# Filter to subset
subset_data <- full_data[gage_id %in% test_gages]
log_info("Subset rows:", nrow(subset_data))

# Write subset to temporary parquet
temp_parquet <- file.path(main_dir, "test_output/subset_streamflow.parquet")
arrow::write_parquet(subset_data, temp_parquet)
log_info("Wrote subset parquet to:", temp_parquet)

# Run signature extraction on subset
log_info("Running process_signatures_from_parquet on subset...")

summary_output <- process_signatures_from_parquet(
  parquet_file_path = temp_parquet,
  metadata_file_path = metadata_path,
  output_file = output_path,
  min_Q_value_and_days = c(0.0001, 30),
  min_num_years = 20,
  min_frac_good_data = 0.95
)

# Summary
cat("\n========== SMOKE TEST RESULTS ==========\n")
if (nrow(summary_output) > 0) {
  cat("Successfully processed:", nrow(summary_output), "gages\n")
  cat("Output columns:", ncol(summary_output), "\n")
  cat("Sample signature values:\n")
  print(summary_output[1, .(gage_id, Qann_mean, Qann_senn_slp, flashinessRB_mean, BFI_Eckhardt_mean)])
  cat("\nSTATUS: SMOKE TEST PASSED\n")
} else {
  cat("STATUS: SMOKE TEST FAILED - No gages processed\n")
}

# ==============================================================================
# CLIMATE SIGNATURES TEST (Optional - requires Daymet parquet)
# ==============================================================================

daymet_parquet <- file.path(main_dir, DAYMET_PARQUET_PATH)

if (file.exists(daymet_parquet)) {
  cat("\n========== CLIMATE SIGNATURES SMOKE TEST ==========\n")

  output_path_climate <- file.path(main_dir, "test_output/smoke_test_climate_signatures.csv")

  log_info("Running process_signatures_from_parquet with Daymet integration...")

  summary_with_climate <- process_signatures_from_parquet(
    parquet_file_path = temp_parquet,
    metadata_file_path = metadata_path,
    output_file = output_path_climate,
    daymet_parquet_path = daymet_parquet,
    min_Q_value_and_days = c(0.0001, 30),
    min_num_years = 20,
    min_frac_good_data = 0.95
  )

  cat("\n========== CLIMATE SIGNATURES RESULTS ==========\n")

  # Check for climate signature columns
  climate_cols <- c("elasticity_static", "elasticity_mean",
                    "qp_slope_sd_mean", "qp_bimodality_mean",
                    "avg_storage_mean", "annual_runoff_ratio_mean")

  found_cols <- intersect(climate_cols, names(summary_with_climate))
  cat("Climate signature columns found:", length(found_cols), "/", length(climate_cols), "\n")

  for (col in found_cols) {
    non_na <- sum(!is.na(summary_with_climate[[col]]))
    cat(sprintf("  %s: %d/%d non-NA values\n",
                col, non_na, nrow(summary_with_climate)))
  }

  if (length(found_cols) > 0 && any(!is.na(summary_with_climate[[found_cols[1]]]))) {
    cat("\nSample climate signature values:\n")
    cols_to_show <- intersect(c("gage_id", found_cols), names(summary_with_climate))
    print(summary_with_climate[1, ..cols_to_show])
    cat("\nCLIMATE SIGNATURES STATUS: PASSED\n")
  } else {
    cat("\nCLIMATE SIGNATURES STATUS: No data populated (check Daymet-gage ID matching)\n")
  }
} else {
  cat("\n[SKIP] Climate signatures test - Daymet parquet not found\n")
  cat("       To enable: run convert_daymet_zip_to_parquet() first\n")
}

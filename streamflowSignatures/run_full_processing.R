################################################################################
# Full Streamflow Signature Processing Script
#
# Processes all available streamflow gages with climate data integration
# Output: data_out/streamflow_signatures_full_MAR2026.csv
#
# Expected runtime: 30-60 minutes for ~6,000 gages
#
# NOTE: This script expects to be run from the streamflowSignatures directory
# Paths are configured in config.R - edit PARQUET_DATA_DIR there for your system
################################################################################

# Clear environment
rm(list = ls())

# Load configuration and functions
cat("========== LOADING CONFIGURATION ==========\n")
source("config.R")
source("R/helperFunctions.R")

# Configure logging
set_log_level("INFO")
log_file_path <- "data_out/processing_log_MAR2026.txt"
set_log_file(log_file_path)

cat("\n========== FULL SIGNATURE EXTRACTION ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Define paths using config variables
parquet_path <- file.path(PARQUET_DATA_DIR, "combined_streamflow_data_09feb2026.parquet")
metadata_path <- file.path(PARQUET_DATA_DIR, "combined_watershed_metadata_09feb2026.csv")
daymet_path <- DAYMET_PARQUET_PATH
output_path <- "data_out/streamflow_signatures_full_MAR2026.csv"

# Verify input files exist
cat("Verifying input files...\n")
if (!file.exists(parquet_path)) {
  stop("Streamflow parquet not found: ", parquet_path,
       "\nEdit PARQUET_DATA_DIR in config.R or set STREAMFLOW_PARQUET_DIR environment variable")
}
if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path)
}
if (!file.exists(daymet_path)) {
  warning("Daymet parquet not found - climate signatures will be skipped")
  daymet_path <- NULL
}

# Check input data
cat("\nInput data summary:\n")
cat("  Streamflow parquet:", parquet_path, "\n")
cat("  Metadata file:", metadata_path, "\n")
cat("  Daymet parquet:", ifelse(is.null(daymet_path), "NOT AVAILABLE", daymet_path), "\n")
cat("  Output file:", output_path, "\n\n")

# Create output directory if needed
dir.create("data_out", showWarnings = FALSE)

# Run full processing
log_info("Starting full signature extraction...")
log_info("Parameters: MIN_NUM_YEARS =", MIN_NUM_YEARS,
         ", MIN_FRAC_GOOD_DATA =", MIN_FRAC_GOOD_DATA,
         ", MIN_Q_VALUE_AND_DAYS =", paste(MIN_Q_VALUE_AND_DAYS, collapse = ","))

start_time <- Sys.time()

summary_output <- process_signatures_from_parquet(
  parquet_file_path = parquet_path,
  metadata_file_path = metadata_path,
  output_file = output_path,
  daymet_parquet_path = daymet_path,
  min_Q_value_and_days = MIN_Q_VALUE_AND_DAYS,
  min_num_years = MIN_NUM_YEARS,
  min_frac_good_data = MIN_FRAC_GOOD_DATA
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

# Summary
cat("\n========== PROCESSING COMPLETE ==========\n")
cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Elapsed time:", round(as.numeric(elapsed), 1), "minutes\n\n")

if (nrow(summary_output) > 0) {
  cat("Results:\n")
  cat("  Gages processed:", nrow(summary_output), "\n")
  cat("  Columns in output:", ncol(summary_output), "\n")
  cat("  Output saved to:", output_path, "\n")
  cat("  Log saved to:", log_file_path, "\n")

  # Quick stats
  cat("\nQuick statistics:\n")
  if ("gage_type" %in% names(summary_output)) {
    cat("  By gage type:\n")
    print(table(summary_output$gage_type))
  }
  if ("num_water_years" %in% names(summary_output)) {
    cat("  Water years - min:", min(summary_output$num_water_years, na.rm = TRUE),
        ", median:", median(summary_output$num_water_years, na.rm = TRUE),
        ", max:", max(summary_output$num_water_years, na.rm = TRUE), "\n")
  }

  # Check climate signature coverage
  if ("elasticity_static" %in% names(summary_output)) {
    climate_coverage <- sum(!is.na(summary_output$elasticity_static)) / nrow(summary_output) * 100
    cat("  Climate signature coverage:", round(climate_coverage, 1), "%\n")
  }

  cat("\nSTATUS: SUCCESS\n")
} else {
  cat("STATUS: FAILED - No gages processed\n")
}

cat("\n========== DONE ==========\n")

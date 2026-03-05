################################################################################
# Restricted Water Year Streamflow Signature Processing Script
#
# Processes gages restricted to water years 1993-2022 with 20-year minimum coverage
# Output: D:/processedOuts_feb2026/streamflow_signatures_full_1993-to-2022-min-20yrs.csv
#
# Comparison baseline: D:/processedOuts_feb2026/streamflow_signatures_full_10feb2026.csv
#
# NOTE: This script expects to be run from the streamflowSignatures directory
################################################################################

# Clear environment
rm(list = ls())

# Load configuration and functions
cat("========== LOADING CONFIGURATION ==========\n")
source("config.R")
source("R/helperFunctions.R")

# Configure logging
set_log_level("INFO")
log_file_path <- "D:/processedOuts_feb2026/processing_log_1993-to-2022-min-20yrs_v2.txt"
set_log_file(log_file_path)

# ============== ANALYSIS PARAMETERS ==============
START_WATER_YEAR <- 1993
END_WATER_YEAR <- 2022
MIN_YEARS_REQUIRED <- 20  # Require at least 20 qualifying years in the range

cat("\n========== RESTRICTED SIGNATURE EXTRACTION ==========\n")
cat("Water year range:", START_WATER_YEAR, "to", END_WATER_YEAR, "\n")
cat("Minimum years required:", MIN_YEARS_REQUIRED, "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Define paths
# Note: Using FEB2026 parquet to match baseline - the older parquet has corrupted Q values
# for Canadian gages without basin area (99999 multiplier bug)
parquet_path <- "D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet"
metadata_path <- "D:/processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv"
daymet_path <- "D:/processedOuts_feb2026/daymet_1980_2023.parquet"
output_path <- "D:/processedOuts_feb2026/streamflow_signatures_full_1993-to-2022-min-20yrs_v2.csv"

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
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

# Run restricted processing
log_info("Starting restricted signature extraction...")
log_info("Parameters: START_WATER_YEAR =", START_WATER_YEAR,
         ", END_WATER_YEAR =", END_WATER_YEAR,
         ", MIN_YEARS_REQUIRED =", MIN_YEARS_REQUIRED)
log_info("Additional filters: MIN_FRAC_GOOD_DATA =", MIN_FRAC_GOOD_DATA,
         ", MIN_Q_VALUE_AND_DAYS =", paste(MIN_Q_VALUE_AND_DAYS, collapse = ","))

start_time <- Sys.time()

summary_output <- process_signatures_from_parquet(
  parquet_file_path = parquet_path,
  metadata_file_path = metadata_path,
  output_file = output_path,
  daymet_parquet_path = daymet_path,
  min_Q_value_and_days = MIN_Q_VALUE_AND_DAYS,
  min_num_years = MIN_YEARS_REQUIRED,
  min_frac_good_data = MIN_FRAC_GOOD_DATA,
  start_water_year = START_WATER_YEAR,
  end_water_year = END_WATER_YEAR
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
  if ("start_water_year" %in% names(summary_output)) {
    cat("  Actual start water year range:", min(summary_output$start_water_year, na.rm = TRUE),
        "to", max(summary_output$start_water_year, na.rm = TRUE), "\n")
  }
  if ("end_water_year" %in% names(summary_output)) {
    cat("  Actual end water year range:", min(summary_output$end_water_year, na.rm = TRUE),
        "to", max(summary_output$end_water_year, na.rm = TRUE), "\n")
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

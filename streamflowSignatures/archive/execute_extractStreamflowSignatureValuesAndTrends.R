######################################################################################## 
# Main execution:
# Process streamflow timeseries (from a pre-processed parquet file) to extract 
# streamflow signature average values and trends
#
# NOTE: This script expects to be run from the streamflowSignatures directory
# Paths are configured in config.R - edit PARQUET_DATA_DIR and METADATA_DATA_DIR there
#
# For most use cases, run_full_processing.R is the recommended entry point.

# Load configuration and helper functions
source("config.R")           # Centralized configuration
source("R/helperFunctions.R")  # All core functions

# Define paths using config variables
parquet_path <- file.path(PARQUET_DATA_DIR, "combined_streamflow_data.parquet")
metadata_path <- file.path(PARQUET_DATA_DIR, "combined_watershed_metadata.csv")
output_path <- file.path(PARQUET_DATA_DIR, "streamflowSignature_summaryData_OCT2025.csv")

# Process signatures from the concatenated parquet file
summary_output <- process_signatures_from_parquet(
  parquet_file_path = parquet_path,
  metadata_file_path = metadata_path,
  output_file = output_path,
  min_Q_value_and_days = MIN_Q_VALUE_AND_DAYS,
  min_num_years = MIN_NUM_YEARS,
  min_frac_good_data = MIN_FRAC_GOOD_DATA
)

# View summary statistics
if (nrow(summary_output) > 0) {
  cat("\nSummary of processed data:\n")
  cat("Total gages processed:", nrow(summary_output), "\n")
  cat("Columns in output:", ncol(summary_output), "\n")
  cat("Column names:", paste(names(summary_output)[1:min(10, ncol(summary_output))], collapse = ", "), "...\n")
}

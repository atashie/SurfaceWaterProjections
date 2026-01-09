######################################################################################## 
# Main execution:
# process streamflow timeseries (from a pre-processed parquet file) to extract streamflow signature average values and trends


# !!!!!! USER INPUT !!!!!!
# Set paths based on your file structure - i.e, where are your data / outputs?
metadata_dir <- "D:"
parquet_path <- file.path(metadata_dir, "combined_streamflow_output/combined_streamflow_data.parquet")
metadata_path <- file.path(metadata_dir, "combined_streamflow_output/combined_watershed_metadata.csv")
output_path <- file.path(metadata_dir, "combined_streamflow_output/streamflowSignature_summaryData_OCT2025.csv")

# !!!!!! USER INPUT !!!!!!
# Source helper functions - i.e., where is your code?
main_dir <- "C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"
source(file.path(main_dir, "config.R"))           # Centralized configuration
source(file.path(main_dir, "helperFunctions.R"))  # All core functions

# Process signatures from the concatenated parquet file
summary_output <- process_signatures_from_parquet(
  parquet_file_path = parquet_path,
  metadata_file_path = metadata_path,
  output_file = output_path,
  min_Q_value_and_days = c(0.0001, 30), # what is the minimum number of day that must exceed some minimum q (in mm) per WY in order to include a WY for analysis?
  min_num_years = 20, # what is the minimum number of "good" WYs in order to include a watershed for analysis?
  min_frac_good_data = 0.95 # what % of days in a WY must non-NA in order to include a WY for analysis?
)

# View summary statistics
if (nrow(summary_output) > 0) {
  cat("\nSummary of processed data:\n")
  cat("Total gages processed:", nrow(summary_output), "\n")
  cat("Columns in output:", ncol(summary_output), "\n")
  cat("Column names:", paste(names(summary_output)[1:min(10, ncol(summary_output))], collapse = ", "), "...\n")
}



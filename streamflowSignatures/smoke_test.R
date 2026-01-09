################################################################################################
# Smoke Test - Uses existing parquet data with a subset of gages
################################################################################################

# Paths
main_dir <- "C:/Users/18033/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"
metadata_dir <- "D:"

# Load config and functions
source(file.path(main_dir, "config.R"))
source(file.path(main_dir, "helperFunctions.R"))

# Set log level
set_log_level("INFO")

# Input/output paths
parquet_path <- file.path(metadata_dir, "combined_streamflow_output/combined_streamflow_data.parquet")
metadata_path <- file.path(metadata_dir, "combined_streamflow_output/combined_watershed_metadata.csv")
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
  print(summary_output[1, .(gage_id, Qann_mean, Qann_slp, flashinessRB_mean, BFI_Eckhardt_mean)])
  cat("\nSTATUS: SMOKE TEST PASSED\n")
} else {
  cat("STATUS: SMOKE TEST FAILED - No gages processed\n")
}

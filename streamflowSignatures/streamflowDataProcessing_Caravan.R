# Caravan Data Processing Script (Legacy)
# Processes Caravan data using the older process_caravan_gages function
#
# NOTE: This script expects to be run from the streamflowSignatures directory
# For new work, consider using caravan_to_annualized.R instead
# Paths are configured in config.R - edit CARAVAN_DATA_DIR there for your system

# Load configuration and helper functions
source("config.R")           # Centralized configuration (includes CARAVAN_DATA_DIR)
source("helperFunctions.R")

# Use Caravan directory from config (can be overridden via environment variable)
caravan_main_dir <- CARAVAN_DATA_DIR

# Define output file for a specific Caravan project
caravan_output_file_camels <- file.path(".", "summary_data_caravan_hysets.csv")

# User configuration
start_date_analysis <- as.Date("1970-01-01")
end_date_analysis <- as.Date("2024-12-31")
min_num_years_analysis <- 20   # minimum period of record (years) required for inclusion in analysis
min_Q_val_days <- c(0.0001, 30) # minimum flow (mm) required before setting to NA

# Example: Process the "hysets" data_project from Caravan
cat("\n========== PROCESSING CARAVAN HYSETS DATA ==========\n")
summary_camels_caravan <- process_caravan_gages(
  data_project_arg = "hysets", # e.g., "camels", "hysets", etc.
  caravan_base_dir = caravan_main_dir,
  min_num_years = min_num_years_analysis,
  start_date = start_date_analysis,
  end_date = end_date_analysis,
  min_Q_value_and_days = min_Q_val_days,
  output_file = caravan_output_file_camels
)

cat("\n========== CARAVAN ANALYSIS COMPLETE ==========\n")
if (exists("summary_camels_caravan") && nrow(summary_camels_caravan) > 0) {
  cat("Completed processing Caravan 'hysets' project. Final summary data has", nrow(summary_camels_caravan), "rows\n")
  cat("Results saved to:", caravan_output_file_camels, "\n")
}

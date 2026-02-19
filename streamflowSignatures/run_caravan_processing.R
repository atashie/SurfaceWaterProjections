# Caravan Data Annualization Script
# Processes Caravan NetCDF data into annualized signatures
#
# NOTE: This script expects to be run from the streamflowSignatures directory
# Paths are configured in config.R - edit CARAVAN_DATA_DIR there for your system

# Load configuration and helper functions
source("config.R")           # Centralized configuration (includes CARAVAN_DATA_DIR)
source("R/helperFunctions.R")

# Use Caravan directory from config (can be overridden via environment variable)
caravan_main_dir <- CARAVAN_DATA_DIR

# Define output file for a specific Caravan project
caravan_output_file_camels <- file.path(".", "summary_data_caravan_camels.csv")
caravan_input_files_camels <- file.path(caravan_main_dir, "timeseries/netcdf/camels")

# User configuration
start_date_analysis <- as.Date("1979-09-01")
end_date_analysis <- as.Date("2025-06-01")
min_num_years_analysis <- 30   # minimum period of record (years) required for inclusion in analysis
min_Q_val_days <- c(0.00001, 30) # minimum flow (mm) required before setting to NA

# Run processing
process_caravan_to_annual(caravan_directory = caravan_input_files_camels,
                          data_project = "camels",
                          min_num_years_data = 30,
                          start_date_filter = as.Date("1979-09-01"),
                          end_date_filter = as.Date("2025-06-01"),
                          output_dir = "annualized_caravan_data"
)

# Daymet ZIP to Parquet conversion script
#
# NOTE: This script expects to be run from the streamflowSignatures directory
# Run with: Rscript R/run_conversion.R

cat("Loading required packages...\n")
library(data.table)
library(arrow)
library(lubridate)

cat("Loading config and helper functions...\n")
source("config.R")
source("R/helperFunctions.R")

cat("Starting Daymet ZIP to Parquet conversion...\n")
cat("This may take 10-20 minutes for 4.3GB of data...\n\n")

convert_daymet_zip_to_parquet(
  daymet_zip_path = DAYMET_ZIP_PATH,
  output_parquet_path = DAYMET_PARQUET_PATH,
  years = DAYMET_START_YEAR:DAYMET_END_YEAR
)

cat("\n========================================\n")
cat("Conversion complete!\n")
cat("========================================\n")

# Verify output
df <- arrow::read_parquet(DAYMET_PARQUET_PATH)
cat("\nVerification:\n")
cat("  Total rows:", nrow(df), "\n")
cat("  Unique sites:", length(unique(df$site_id)), "\n")
cat("  Date range:", as.character(min(df$Date)), "to", as.character(max(df$Date)), "\n")
cat("  File size:", round(file.info(DAYMET_PARQUET_PATH)$size / 1024 / 1024, 1), "MB\n")

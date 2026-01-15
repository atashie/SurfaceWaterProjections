# Daymet ZIP to Parquet conversion script
# Run with: Rscript run_conversion.R

setwd("C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures")

cat("Loading required packages...\n")
library(data.table)
library(arrow)
library(lubridate)

cat("Loading config and helper functions...\n")
source("config.R")
source("helperFunctions.R")

cat("Starting Daymet ZIP to Parquet conversion...\n")
cat("This may take 10-20 minutes for 4.3GB of data...\n\n")

convert_daymet_zip_to_parquet(
  daymet_zip_path = "data_out/daymet_1980_2023.zip",
  output_parquet_path = "data_out/daymet_1980_2023.parquet",
  years = 1980:2023
)

cat("\n========================================\n")
cat("Conversion complete!\n")
cat("========================================\n")

# Verify output
df <- arrow::read_parquet("data_out/daymet_1980_2023.parquet")
cat("\nVerification:\n")
cat("  Total rows:", nrow(df), "\n")
cat("  Unique sites:", length(unique(df$site_id)), "\n")
cat("  Date range:", as.character(min(df$Date)), "to", as.character(max(df$Date)), "\n")
cat("  File size:", round(file.info("data_out/daymet_1980_2023.parquet")$size / 1024 / 1024, 1), "MB\n")

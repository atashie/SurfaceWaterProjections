################################################################################
# Regenerate Combined Watershed Metadata
#
# This script regenerates combined_watershed_metadata.csv from the individual
# processed directory metadata files, then enriches with human interference data.
################################################################################

cat("========== REGENERATING COMBINED METADATA ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

source("config.R")
source("R/helperFunctions.R")

input_directories <- c(
  "processed_streamflow_data_AK",
  "processed_streamflow_data_canada",
  "processed_streamflow_data_conus"
)

output_dir <- "combined_streamflow_output"

# Create output directory if needed
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Step 1: Combine metadata files (without parquet concatenation to avoid crash)
cat("Step 1: Combining metadata files...\n")
all_metadata <- data.table()

for (dir in input_directories) {
  metadata_file <- file.path(dir, "watershed_metadata.csv")
  if (file.exists(metadata_file)) {
    cat("  Reading:", metadata_file, "\n")
    metadata <- fread(metadata_file, colClasses = c("gage_id" = "character"))
    metadata$source_dir <- basename(dir)
    all_metadata <- rbind(all_metadata, metadata, fill = TRUE)
    cat("    Found", nrow(metadata), "gages\n")
  } else {
    cat("  WARNING: Not found:", metadata_file, "\n")
  }
}

# Save combined metadata
output_file <- file.path(output_dir, "combined_watershed_metadata.csv")
fwrite(all_metadata, output_file)
cat("\nCombined metadata saved to:", output_file, "\n")
cat("Total gages:", nrow(all_metadata), "\n")

# Print summary
cat("\nBy source directory:\n")
print(table(all_metadata$source_dir))

cat("\nBy processing status:\n")
print(table(all_metadata$processing_status))

# Step 2: Enrich with human interference metadata
cat("\n========== ENRICHING WITH HUMAN INTERFERENCE DATA ==========\n")

enriched_metadata <- enrich_metadata_with_interference(
  metadata_file_path = output_file,
  gages_ii_dir = GAGES_II_DIR
)

# Final summary
cat("\n========== FINAL SUMMARY ==========\n")
cat("Total gages:", nrow(enriched_metadata), "\n")
cat("\nBy gage type:\n")
print(table(enriched_metadata$gage_type))
cat("\nBy human_interference_class:\n")
print(table(enriched_metadata$human_interference_class, useNA = "ifany"))

cat("\n========== DONE ==========\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

################################################################################
# Metadata Enrichment Script
#
# Adds human interference metadata from GAGES-II (USGS) and Canadian HYDAT
# to the combined watershed metadata file.
#
# Run this BEFORE run_full_processing.R to ensure interference data is available
# in the final signature output.
#
# Usage:
#   Rscript run_enrich_metadata.R
#
# Or in R:
#   source("run_enrich_metadata.R")
#
# Prerequisites:
#   - combined_watershed_metadata.csv must exist in PARQUET_DATA_DIR
#   - GAGES-II metadata files must exist in GAGES_II_DIR (D:/gagesMetadata)
#   - tidyhydat package installed for Canadian metadata
################################################################################

cat("========== METADATA ENRICHMENT ==========\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Source configuration and helper functions
source("config.R")
source("R/helperFunctions.R")

# Define paths
gages_ii_path <- GAGES_II_DIR

# Find the most recent combined_watershed_metadata file
metadata_candidates <- list.files(
  PARQUET_DATA_DIR,
  pattern = "^combined_watershed_metadata.*\\.csv$",
  full.names = TRUE
)

if (length(metadata_candidates) == 0) {
  stop("No combined_watershed_metadata*.csv found in: ", PARQUET_DATA_DIR,
       "\n  Make sure PARQUET_DATA_DIR is set correctly in config.R",
       "\n  Current value: ", PARQUET_DATA_DIR)
}

# Use most recently modified file
metadata_path <- metadata_candidates[which.max(file.mtime(metadata_candidates))]

# Print configuration
cat("Configuration:\n")
cat("  Metadata file:", metadata_path, "\n")
cat("  GAGES-II directory:", gages_ii_path, "\n\n")

# Verify inputs exist (metadata_path already confirmed above)
if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path)
}

if (!dir.exists(gages_ii_path)) {
  stop("GAGES-II directory not found: ", gages_ii_path,
       "\n  Make sure GAGES_II_DIR is set correctly in config.R or as environment variable")
}

# Run enrichment
enriched_metadata <- enrich_metadata_with_interference(
  metadata_file_path = metadata_path,
  gages_ii_dir = gages_ii_path
)

# Print detailed summary
cat("\n========== ENRICHMENT COMPLETE ==========\n\n")

cat("Total gages:", nrow(enriched_metadata), "\n\n")

# By gage type
cat("By gage type:\n")
print(table(enriched_metadata$gage_type, useNA = "ifany"))

cat("\n")

# By interference class
cat("By human_interference_class:\n")
print(table(enriched_metadata$human_interference_class, useNA = "ifany"))

cat("\n")

# USGS-specific stats
usgs_gages <- enriched_metadata[gage_type == "USGS"]
if (nrow(usgs_gages) > 0) {
  cat("USGS gages with GAGES-II data:\n")
  cat("  Total USGS gages:", nrow(usgs_gages), "\n")
  cat("  With dam count (NDAMS_2009):", sum(!is.na(usgs_gages$NDAMS_2009)), "\n")
  cat("  With disturbance index:", sum(!is.na(usgs_gages$HYDRO_DISTURB_INDX)), "\n")
  cat("  With impervious surface:", sum(!is.na(usgs_gages$IMPNLCD06)), "\n")
  cat("  With developed area:", sum(!is.na(usgs_gages$DEVNLCD06)), "\n")
  cat("  Reference (CLASS=Ref):", sum(usgs_gages$CLASS == "Ref", na.rm = TRUE), "\n")
  cat("  Non-reference (CLASS=Non-ref):", sum(usgs_gages$CLASS == "Non-ref", na.rm = TRUE), "\n")
  cat("\n")
}

# Canadian-specific stats
canadian_gages <- enriched_metadata[gage_type %in% c("Canada", "CANADIAN")]
if (nrow(canadian_gages) > 0) {
  cat("Canadian gages with HYDAT data:\n")
  cat("  Total Canadian gages:", nrow(canadian_gages), "\n")
  cat("  With RHBN flag:", sum(!is.na(canadian_gages$RHBN)), "\n")
  cat("  RHBN = TRUE (reference):", sum(canadian_gages$RHBN == TRUE, na.rm = TRUE), "\n")
  cat("  With regulation status:", sum(!is.na(canadian_gages$REGULATED)), "\n")
  cat("  REGULATED = TRUE:", sum(canadian_gages$REGULATED == TRUE, na.rm = TRUE), "\n")
  cat("\n")
}

# Sample of enriched data
cat("Sample of enriched metadata (first 5 rows, key columns):\n")
sample_cols <- c("gage_id", "gage_type", "NDAMS_2009", "HYDRO_DISTURB_INDX",
                 "CLASS", "RHBN", "REGULATED", "human_interference_class")
available_sample_cols <- intersect(sample_cols, names(enriched_metadata))
print(head(enriched_metadata[, ..available_sample_cols], 5))

cat("\n========== DONE ==========\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nNext step: Run run_full_processing.R to generate signatures with interference metadata\n")

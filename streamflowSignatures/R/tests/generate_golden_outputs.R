################################################################################################
# Generate Golden Outputs for Cross-Language Validation
#
# This script:
# 1. Selects a representative subset of gages
# 2. Saves test data to test-data/
# 3. Runs R signature extraction
# 4. Saves golden outputs to golden-outputs/
################################################################################################

# Paths - auto-detect based on script location
main_dir <- getwd()
if (!file.exists(file.path(main_dir, "config.R"))) {
  main_dir <- dirname(main_dir)
  if (!file.exists(file.path(main_dir, "config.R"))) {
    main_dir <- dirname(main_dir)
  }
}

# Load config and functions
source(file.path(main_dir, "config.R"))
source(file.path(main_dir, "R", "helperFunctions.R"))

set_log_level("INFO")

cat("========== GOLDEN OUTPUT GENERATION ==========\n\n")

# ==============================================================================
# Configuration
# ==============================================================================

# Number of test gages (representative sample)
N_TEST_GAGES <- 15

# Paths
parquet_path <- "D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet"
metadata_path <- "D:/processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv"
daymet_path <- "D:/processedOuts_feb2026/daymet_1980_2023.parquet"

# Output paths
test_data_dir <- file.path(main_dir, "test-data")
golden_output_dir <- file.path(main_dir, "golden-outputs")

# Create directories
dir.create(test_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(golden_output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Step 1: Select Representative Test Gages
# ==============================================================================

log_info("Reading full parquet to select test gages...")

if (!file.exists(parquet_path)) {
  stop("Parquet file not found: ", parquet_path,
       "\nPlease update parquet_path in this script.")
}

full_data <- arrow::read_parquet(parquet_path)
full_data <- as.data.table(full_data)

log_info("Total gages in parquet:", length(unique(full_data$gage_id)))

# Get gage statistics for selection
gage_stats <- full_data[, .(
  n_years = length(unique(water_year)),
  mean_Q = mean(Q, na.rm = TRUE),
  cv_Q = sd(Q, na.rm = TRUE) / mean(Q, na.rm = TRUE),
  na_frac = sum(is.na(Q)) / .N
), by = gage_id]

# Filter to gages with good data (20+ years, <5% NA)
good_gages <- gage_stats[n_years >= 20 & na_frac < 0.05]
log_info("Gages meeting criteria:", nrow(good_gages))

# Select diverse sample: low/medium/high flow variability
set.seed(42)  # Reproducibility

# Stratified sampling by coefficient of variation
good_gages[, cv_group := cut(cv_Q, quantile(cv_Q, c(0, 0.33, 0.67, 1), na.rm = TRUE),
                              labels = c("low", "medium", "high"), include.lowest = TRUE)]

test_gages <- c()
for (grp in c("low", "medium", "high")) {
  grp_gages <- good_gages[cv_group == grp]$gage_id
  n_select <- min(5, length(grp_gages))
  test_gages <- c(test_gages, sample(grp_gages, n_select))
}

# Ensure we have the target number
if (length(test_gages) < N_TEST_GAGES) {
  remaining <- setdiff(good_gages$gage_id, test_gages)
  n_more <- min(N_TEST_GAGES - length(test_gages), length(remaining))
  test_gages <- c(test_gages, sample(remaining, n_more))
}

test_gages <- head(test_gages, N_TEST_GAGES)

cat("\nSelected test gages:\n")
print(test_gages)
cat("\n")

# ==============================================================================
# Step 2: Extract and Save Test Data
# ==============================================================================

log_info("Extracting test data for", length(test_gages), "gages...")

# Filter streamflow data
test_streamflow <- full_data[gage_id %in% test_gages]
log_info("Test streamflow rows:", nrow(test_streamflow))

# Save to parquet
test_parquet_path <- file.path(test_data_dir, "sample_streamflow.parquet")
arrow::write_parquet(test_streamflow, test_parquet_path)
log_info("Saved streamflow to:", test_parquet_path)

# Extract and save metadata
if (file.exists(metadata_path)) {
  metadata <- fread(metadata_path)
  test_metadata <- metadata[gage_id %in% test_gages |
                            as.character(gage_id) %in% as.character(test_gages)]

  # Handle leading zeros in gage IDs
  if (nrow(test_metadata) == 0) {
    test_metadata <- metadata[sapply(gage_id, function(x) {
      any(grepl(paste0("^0*", x, "$"), test_gages) |
          grepl(paste0("^0*", x, "$"), as.character(test_gages)))
    })]
  }

  test_metadata_path <- file.path(test_data_dir, "sample_metadata.csv")
  fwrite(test_metadata, test_metadata_path)
  log_info("Saved metadata to:", test_metadata_path)
} else {
  log_warn("Metadata file not found:", metadata_path)
  test_metadata_path <- NULL
}

# Extract and save climate data if available
if (file.exists(daymet_path)) {
  log_info("Extracting Daymet climate data...")
  daymet_data <- arrow::read_parquet(daymet_path)
  daymet_data <- as.data.table(daymet_data)

  # Filter to test gages
  test_climate <- daymet_data[gage_id %in% test_gages |
                               as.character(gage_id) %in% as.character(test_gages)]

  if (nrow(test_climate) > 0) {
    test_climate_path <- file.path(test_data_dir, "sample_climate.parquet")
    arrow::write_parquet(test_climate, test_climate_path)
    log_info("Saved climate data to:", test_climate_path)
  } else {
    log_warn("No matching climate data found for test gages")
  }
} else {
  log_info("Daymet parquet not found - skipping climate data")
}

# ==============================================================================
# Step 3: Generate Golden Outputs (without climate)
# ==============================================================================

log_info("\nRunning R signature extraction (without climate)...")

golden_output_path <- file.path(golden_output_dir, "expected_signatures.csv")

summary_output <- process_signatures_from_parquet(
  parquet_file_path = test_parquet_path,
  metadata_file_path = if (!is.null(test_metadata_path)) test_metadata_path else metadata_path,
  output_file = golden_output_path,
  min_Q_value_and_days = c(0.0001, 30),
  min_num_years = 20,
  min_frac_good_data = 0.95
)

log_info("Golden output saved to:", golden_output_path)
log_info("Gages processed:", nrow(summary_output))
log_info("Signature columns:", ncol(summary_output))

# ==============================================================================
# Step 4: Generate Golden Outputs (with climate if available)
# ==============================================================================

test_climate_path <- file.path(test_data_dir, "sample_climate.parquet")

if (file.exists(test_climate_path)) {
  log_info("\nRunning R signature extraction (with climate)...")

  golden_climate_path <- file.path(golden_output_dir, "expected_signatures_climate.csv")

  summary_climate <- process_signatures_from_parquet(
    parquet_file_path = test_parquet_path,
    metadata_file_path = if (!is.null(test_metadata_path)) test_metadata_path else metadata_path,
    output_file = golden_climate_path,
    daymet_parquet_path = test_climate_path,
    min_Q_value_and_days = c(0.0001, 30),
    min_num_years = 20,
    min_frac_good_data = 0.95
  )

  log_info("Golden output (with climate) saved to:", golden_climate_path)

  # Check climate columns
  climate_cols <- c("elasticity_static", "annual_runoff_ratio_mean",
                    "qp_slope_sd_mean", "avg_storage_mean")
  found <- intersect(climate_cols, names(summary_climate))
  log_info("Climate signature columns found:", length(found))
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n========== GOLDEN OUTPUT GENERATION COMPLETE ==========\n\n")
cat("Test data:\n")
cat("  Streamflow:", test_parquet_path, "\n")
if (!is.null(test_metadata_path)) cat("  Metadata:", test_metadata_path, "\n")
if (file.exists(test_climate_path)) cat("  Climate:", test_climate_path, "\n")
cat("\nGolden outputs:\n")
cat("  Signatures:", golden_output_path, "\n")
if (file.exists(test_climate_path)) {
  cat("  With climate:", golden_climate_path, "\n")
}
cat("\nNext steps:\n")
cat("  1. Run Python: python docs/benchmarks/compare_outputs.py\n")
cat("  2. Check validation_report.md for results\n")

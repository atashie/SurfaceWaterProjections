################################################################################################
# Regression Verification — Compare fresh output against golden reference
#
# Purpose: Verify that code quality changes (R1-R7) produce identical signature values.
# Runs on 20 gages from golden output, re-computes signatures, and checks for deviation.
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
set_log_level("WARN")  # Quiet logging for verification

cat("========== REGRESSION VERIFICATION ==========\n\n")

# Paths
parquet_path <- "D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet"
metadata_path <- "D:/processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv"
golden_path <- file.path(main_dir, "golden-outputs/streamflow_signatures_full_10feb2026.csv")
daymet_path <- "D:/processedOuts_feb2026/daymet_1980_2023.parquet"
output_path <- file.path(main_dir, "test_output/regression_verify.csv")

# Verify all files exist
for (p in c(parquet_path, metadata_path, golden_path, daymet_path)) {
  if (!file.exists(p)) stop("Missing required file: ", p)
}
dir.create(file.path(main_dir, "test_output"), showWarnings = FALSE)

# Load golden output
cat("Loading golden reference output...\n")
golden <- fread(golden_path, colClasses = list(character = "gage_id"))
cat("  Golden output: ", nrow(golden), " gages, ", ncol(golden), " columns\n")

# Select 20 gages from golden output (deterministic seed)
set.seed(123)
test_gages <- sample(golden$gage_id, min(20, nrow(golden)))
cat("  Selected ", length(test_gages), " gages for regression check\n\n")

# Write subset parquet for these gages
cat("Extracting subset from parquet...\n")
full_data <- arrow::read_parquet(parquet_path)
full_data <- as.data.table(full_data)
subset_data <- full_data[gage_id %in% test_gages]

temp_parquet <- file.path(main_dir, "test_output/regression_subset.parquet")
arrow::write_parquet(subset_data, temp_parquet)
rm(full_data)  # free memory
cat("  Subset: ", nrow(subset_data), " rows for ", length(unique(subset_data$gage_id)), " gages\n\n")

# Run signature extraction with climate data (full pipeline)
cat("Running process_signatures_from_parquet (with Daymet)...\n")
fresh_output <- process_signatures_from_parquet(
  parquet_file_path = temp_parquet,
  metadata_file_path = metadata_path,
  output_file = output_path,
  daymet_parquet_path = daymet_path,
  min_Q_value_and_days = MIN_Q_VALUE_AND_DAYS,
  min_num_years = MIN_NUM_YEARS,
  min_frac_good_data = MIN_FRAC_GOOD_DATA
)
cat("\n  Fresh output: ", nrow(fresh_output), " gages, ", ncol(fresh_output), " columns\n\n")

# Compare fresh vs golden
cat("========== COMPARING FRESH vs GOLDEN ==========\n\n")

# Find common gages
fresh_output[, gage_id := as.character(gage_id)]
golden[, gage_id := as.character(gage_id)]
common_gages <- intersect(fresh_output$gage_id, golden$gage_id)
cat("Common gages: ", length(common_gages), "\n")

if (length(common_gages) == 0) {
  cat("FAIL: No common gages found!\n")
  cat("Fresh gages: ", paste(head(fresh_output$gage_id, 5), collapse = ", "), "\n")
  cat("Golden gages: ", paste(head(golden$gage_id, 5), collapse = ", "), "\n")
  quit(status = 1)
}

# Subset to common gages, sorted
fresh_sub <- fresh_output[gage_id %in% common_gages]
golden_sub <- golden[gage_id %in% common_gages]
setkey(fresh_sub, gage_id)
setkey(golden_sub, gage_id)

# Find common signature columns (exclude metadata)
meta_cols <- c("gage_id", "gage_id_metadata", "latitude", "longitude", "basin_area",
               "gage_type", "num_water_years", "start_water_year", "end_water_year",
               "analysis_start_water_year", "analysis_end_water_year",
               "num_years", "start_year", "end_year",
               "num_upstream_basins", "area_normalized",
               "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009", "IMPNLCD06",
               "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX", "CLASS",
               "RHBN", "REGULATED", "human_interference_class", "country")
sig_cols_fresh <- setdiff(names(fresh_sub), meta_cols)
sig_cols_golden <- setdiff(names(golden_sub), meta_cols)
common_sig_cols <- intersect(sig_cols_fresh, sig_cols_golden)
cat("Common signature columns: ", length(common_sig_cols), "\n\n")

# Compare each column
n_perfect <- 0
n_close <- 0
n_different <- 0
different_cols <- character()

for (col in common_sig_cols) {
  fresh_vals <- as.numeric(fresh_sub[[col]])
  golden_vals <- as.numeric(golden_sub[[col]])

  # Both all-NA: perfect match
  if (all(is.na(fresh_vals)) && all(is.na(golden_vals))) {
    n_perfect <- n_perfect + 1
    next
  }

  # NA pattern mismatch
  if (!identical(is.na(fresh_vals), is.na(golden_vals))) {
    n_different <- n_different + 1
    different_cols <- c(different_cols, col)
    na_diff <- sum(is.na(fresh_vals) != is.na(golden_vals))
    cat(sprintf("  DIFF [NA pattern]: %-40s  %d/%d NAs differ\n", col, na_diff, length(fresh_vals)))
    next
  }

  # Numeric comparison on non-NA values
  valid <- !is.na(fresh_vals) & !is.na(golden_vals)
  if (sum(valid) == 0) {
    n_perfect <- n_perfect + 1
    next
  }

  f <- fresh_vals[valid]
  g <- golden_vals[valid]

  # Check exact match
  if (identical(f, g)) {
    n_perfect <- n_perfect + 1
    next
  }

  # Check relative tolerance (1e-10)
  max_abs_diff <- max(abs(f - g))
  max_rel_diff <- max(abs(f - g) / pmax(abs(g), 1e-15))

  if (max_abs_diff < 1e-10) {
    n_perfect <- n_perfect + 1  # Floating point noise
  } else if (max_rel_diff < 1e-6) {
    n_close <- n_close + 1
    cat(sprintf("  CLOSE: %-40s  max_rel_diff=%.2e\n", col, max_rel_diff))
  } else {
    n_different <- n_different + 1
    different_cols <- c(different_cols, col)
    cat(sprintf("  DIFF:  %-40s  max_abs_diff=%.6g  max_rel_diff=%.2e\n", col, max_abs_diff, max_rel_diff))
  }
}

cat("\n========== SUMMARY ==========\n")
cat(sprintf("Perfect match:  %d / %d columns\n", n_perfect, length(common_sig_cols)))
cat(sprintf("Close match:    %d / %d columns  (rel_diff < 1e-6)\n", n_close, length(common_sig_cols)))
cat(sprintf("Different:      %d / %d columns\n", n_different, length(common_sig_cols)))

if (n_different == 0) {
  cat("\nSTATUS: PASS — Zero signature regressions detected.\n")
} else {
  cat("\nSTATUS: FAIL — ", n_different, " columns differ from golden reference.\n")
  cat("Different columns: ", paste(different_cols, collapse = ", "), "\n")
}

# Cleanup temp files
unlink(temp_parquet)
unlink(output_path)

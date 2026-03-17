#!/usr/bin/env Rscript
# ==============================================================================
# R Benchmark - Streamflow Signatures
# ==============================================================================
# This script runs the full R signature extraction workflow and:
# 1. Uses the standard config.R configuration
# 2. Renames Q95.Q10 columns to Q95_Q10 (underscore convention)
# 3. Outputs to docs/benchmarks/r_signatures.csv
# 4. Saves timing info to docs/benchmarks/r_timing.json
# ==============================================================================

cat("======================================================================\n")
cat("R BENCHMARK - Streamflow Signatures\n")
cat("======================================================================\n")

start_time <- Sys.time()
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# Initialize timing
timing <- list(
  start_time = format(start_time, "%Y-%m-%dT%H:%M:%S"),
  phases = list(),
  language = "R"
)

# Set working directory to project root
# Use commandArgs to get script path (works with Rscript)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
} else {
  # Fallback: assume we're in the benchmarks directory
  script_dir <- getwd()
}
project_root <- normalizePath(file.path(script_dir, "../.."))
setwd(project_root)
cat("Working directory:", getwd(), "\n\n")

# Phase 1: Load configuration
cat("Phase 1: Loading configuration...\n")
phase_start <- Sys.time()

source("config.R")
source("R/helperFunctions.R")

timing$phases$load_config <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat("  Config loaded in", round(timing$phases$load_config, 1), "s\n")

# Phase 2: Set up paths
cat("\nPhase 2: Setting up paths...\n")
phase_start <- Sys.time()

parquet_path <- "D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet"
metadata_path <- "D:/processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv"
daymet_path <- "D:/processedOuts_feb2026/daymet_1980_2023.parquet"
output_path <- "docs/benchmarks/r_signatures.csv"

# Verify files exist
if (!file.exists(parquet_path)) {
  stop("Streamflow parquet not found: ", parquet_path)
}
if (!file.exists(metadata_path)) {
  stop("Metadata file not found: ", metadata_path)
}
if (!file.exists(daymet_path)) {
  warning("Climate data not found, running without climate signatures: ", daymet_path)
  daymet_path <- NULL
}

timing$phases$setup <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat("  Paths verified in", round(timing$phases$setup, 1), "s\n")

# Phase 3: Process signatures
cat("\nPhase 3: Processing signatures...\n")
cat("  This may take 30-60 minutes...\n")
phase_start <- Sys.time()

summary_output <- process_signatures_from_parquet(
  parquet_file_path = parquet_path,
  metadata_file_path = metadata_path,
  output_file = output_path,
  daymet_parquet_path = daymet_path,
  min_Q_value_and_days = MIN_Q_VALUE_AND_DAYS,
  min_num_years = MIN_NUM_YEARS,
  min_frac_good_data = MIN_FRAC_GOOD_DATA
)

timing$phases$process_signatures <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat("  Signatures processed in", round(timing$phases$process_signatures, 1), "s (",
    round(timing$phases$process_signatures/60, 1), " min)\n")

# Phase 4: Rename columns (Q95.Q10 -> Q95_Q10)
cat("\nPhase 4: Renaming columns (Q95.Q10 -> Q95_Q10)...\n")
phase_start <- Sys.time()

old_names <- names(summary_output)
new_names <- gsub("Q95\\.Q10", "Q95_Q10", old_names)
n_renamed <- sum(old_names != new_names)

if (n_renamed > 0) {
  setnames(summary_output, old_names, new_names)
  cat("  Renamed", n_renamed, "columns\n")
  
  # Re-save with new column names
  fwrite(summary_output, output_path)
  cat("  Re-saved to", output_path, "\n")
} else {
  cat("  No columns needed renaming\n")
}

timing$phases$rename_save <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))

# Finalize timing
end_time <- Sys.time()
timing$end_time <- format(end_time, "%Y-%m-%dT%H:%M:%S")
timing$total_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
timing$n_gages_processed <- nrow(summary_output)
timing$n_columns <- ncol(summary_output)
timing$n_signature_columns <- sum(sapply(new_names, function(n) {
  any(sapply(c("_mean", "_median", "_senn_slp", "_linear_slp", 
               "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval",
               "elasticity_static", "log_a_seasonality"), 
             function(s) grepl(s, n, fixed = TRUE)))
}))
timing$n_metadata_columns <- ncol(summary_output) - timing$n_signature_columns - 
  sum(grepl("^flagged_", new_names))
timing$n_qaqc_flags <- sum(grepl("^flagged_", new_names))

# Save timing JSON
timing_path <- "docs/benchmarks/r_timing.json"
jsonlite::write_json(timing, timing_path, auto_unbox = TRUE, pretty = TRUE)
cat("  Timing saved to", timing_path, "\n")

cat("\n======================================================================\n")
cat("BENCHMARK COMPLETE\n")
cat("======================================================================\n")
cat("Total time:", round(timing$total_seconds, 1), "s (", 
    round(timing$total_seconds/60, 1), " min)\n")
cat("Gages processed:", timing$n_gages_processed, "\n")
cat("Total columns:", timing$n_columns, "\n")
cat("  Signature columns:", timing$n_signature_columns, "\n")
cat("  Metadata columns:", timing$n_metadata_columns, "\n")
cat("  QA/QC flag columns:", timing$n_qaqc_flags, "\n")
cat("Rate:", round(timing$n_gages_processed / timing$total_seconds, 2), "gages/s\n")

#!/usr/bin/env Rscript
# ==============================================================================
# rpkg Benchmark - Streamflow Signatures
# ==============================================================================
# Runs the rpkg package (streamflowsignatures) on all qualifying gages with
# climate data, producing output in the same format as run_r_benchmark.R,
# run_python_benchmark.py, and run_julia_benchmark.jl.
#
# Output: docs/benchmarks/rpkg_signatures.csv
# ==============================================================================

cat("======================================================================\n")
cat("rpkg BENCHMARK - Streamflow Signatures (streamflowsignatures package)\n")
cat("======================================================================\n")

start_time <- Sys.time()
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

timing <- list(
  start_time = format(start_time, "%Y-%m-%dT%H:%M:%S"),
  phases = list(),
  language = "rpkg"
)

# Set working directory to project root
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}
project_root <- normalizePath(file.path(script_dir, "../.."))
setwd(project_root)
cat("Working directory:", getwd(), "\n\n")

# Phase 1: Load rpkg
cat("Phase 1: Loading rpkg...\n")
phase_start <- Sys.time()

library(devtools)
load_all("rpkg")
library(data.table)

timing$phases$load_config <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat("  rpkg loaded in", round(timing$phases$load_config, 1), "s\n")

# Phase 2: Load data
cat("\nPhase 2: Loading data...\n")
phase_start <- Sys.time()

parquet_path <- "D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet"
daymet_path  <- "D:/processedOuts_feb2026/daymet_1980_2023.parquet"
metadata_path <- "D:/processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv"
output_path  <- "docs/benchmarks/rpkg_signatures.csv"

if (!file.exists(parquet_path)) stop("Parquet not found: ", parquet_path)

cat("  Reading streamflow parquet...\n")
streamflow <- read_parquet(parquet_path)
cat("    ", nrow(streamflow), "rows,", length(unique(streamflow$gage_id)), "gages\n")

cat("  Reading Daymet parquet...\n")
has_daymet <- file.exists(daymet_path)
if (has_daymet) {
  daymet <- read_parquet(daymet_path)
  daymet_gages <- unique(daymet$gage_id)
  cat("    ", nrow(daymet), "rows,", length(daymet_gages), "gages with climate data\n")
} else {
  warning("Daymet not found, running without climate")
  daymet <- NULL
  daymet_gages <- character(0)
}

cat("  Reading metadata...\n")
metadata <- fread(metadata_path)
cat("    ", nrow(metadata), "gages in metadata\n")

timing$phases$load_data <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat("  Data loaded in", round(timing$phases$load_data, 1), "s\n")

# Phase 3: Per-gage filtering and signature extraction
use_legacy <- pkg_env$use_legacy_filtering
cat(sprintf("\nPhase 3: Per-gage filtering and signature extraction (legacy=%s)...\n",
            as.character(use_legacy)))
phase_start <- Sys.time()

all_gages <- unique(streamflow$gage_id)
n_gages <- length(all_gages)
cat("  Total gages to attempt:", n_gages, "\n")

# Pre-index streamflow and daymet by gage_id for fast lookup
streamflow_idx <- split(seq_len(nrow(streamflow)), streamflow$gage_id)

if (has_daymet) {
  daymet_idx <- split(seq_len(nrow(daymet)), daymet$gage_id)
} else {
  daymet_idx <- list()
}

results_list <- vector("list", n_gages)
names(results_list) <- all_gages
processed <- 0
skipped <- 0
errored <- 0

# Cache preprocessing results for non-legacy path
preprocess_cache <- list()

# Progress tracking
report_interval <- 500

for (i in seq_along(all_gages)) {
  gage <- all_gages[i]

  tryCatch({
    # Extract gage data
    rows <- streamflow_idx[[gage]]
    gage_data <- streamflow[rows, ]
    gage_data <- add_water_year_columns(gage_data)

    # Merge climate data if available (needed before preprocessing)
    has_climate <- FALSE
    if (gage %in% names(daymet_idx)) {
      daymet_rows <- daymet_idx[[gage]]
      gage_climate <- daymet[daymet_rows, ]

      # Get PPT column (may be named PPT, prcp, or PRCP)
      ppt_col <- intersect(c("PPT", "prcp", "PRCP"), names(gage_climate))
      if (length(ppt_col) > 0) {
        climate_merge <- gage_climate[, c("date", ppt_col[1]), with = FALSE]
        if (ppt_col[1] != "PPT") setnames(climate_merge, ppt_col[1], "PPT")
        gage_data <- merge(gage_data, climate_merge, by = "date", all.x = TRUE)
        has_climate <- TRUE
      }
    }

    if (use_legacy) {
      # Legacy path: filter qualifying years
      filt <- filter_qualifying_years(gage_data)
      if (!filt$qualifies) {
        skipped <- skipped + 1
        next
      }

      gage_filtered <- gage_data[gage_data$water_year %in% filt$qualifying_years, ]

      # Calculate all signatures
      sigs <- calculate_all_signatures(gage_filtered, has_climate = has_climate)

      # Add metadata
      sigs$gage_id <- gage
      sigs$start_water_year <- min(filt$qualifying_years)
      sigs$end_water_year <- max(filt$qualifying_years)
      sigs$num_water_years <- length(filt$qualifying_years)
    } else {
      # New path: preprocess and use trend completeness
      pp <- preprocess_daily_data(gage_data)
      if (length(pp$valid_years) < pkg_env$min_num_years) {
        skipped <- skipped + 1
        next
      }

      pp_data <- pp$data
      gage_has_climate <- has_climate && length(pp$valid_climate_years) > 0

      # Filter to valid_climate_years for climate signatures
      climate_data_subset <- NULL
      if (gage_has_climate) {
        climate_data_subset <- pp_data[pp_data$water_year %in% pp$valid_climate_years, ]
      }

      sigs <- calculate_all_signatures(
        pp_data, has_climate = gage_has_climate,
        seasonal_flags = pp$seasonal_flags,
        trend_completeness = pkg_env$na_trend_min_fraction,
        decade_completeness = pkg_env$na_decade_min_fraction,
        climate_data = climate_data_subset
      )

      # Add metadata
      sigs$gage_id <- gage
      sigs$start_water_year <- min(pp$valid_years)
      sigs$end_water_year <- max(pp$valid_years)
      sigs$num_water_years <- length(pp$valid_years)
    }

    # Match metadata row
    gage_stripped <- sub("^0+", "", gage)
    meta_row <- metadata[metadata$gage_id == gage | as.character(metadata$gage_id) == gage_stripped, ]
    if (nrow(meta_row) > 0) {
      meta_row <- meta_row[1, ]
      for (col in c("latitude", "longitude", "basin_area", "gage_type", "area_normalized")) {
        if (col %in% names(meta_row)) sigs[[col]] <- meta_row[[col]]
      }
    }

    # QA/QC flags
    qa <- compute_qa_flags(sigs)
    for (nm in names(qa)) {
      sigs[[paste0("flagged_for_", nm)]] <- qa[[nm]]
    }

    results_list[[gage]] <- sigs
    processed <- processed + 1

  }, error = function(e) {
    errored <<- errored + 1
    if (errored <= 10) {
      cat("  ERROR on gage", gage, ":", e$message, "\n")
    }
  })

  # Progress report
  if (i %% report_interval == 0 || i == n_gages) {
    elapsed <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
    rate <- i / elapsed
    eta <- (n_gages - i) / rate
    cat(sprintf("  [%d/%d] processed=%d skipped=%d errors=%d | %.1f gages/s | ETA %.0fs\n",
                i, n_gages, processed, skipped, errored, rate, eta))
  }
}

timing$phases$process_signatures <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat(sprintf("  Signatures processed in %.1fs (%.1f min)\n",
            timing$phases$process_signatures, timing$phases$process_signatures / 60))

# Phase 4: Assemble output
cat("\nPhase 4: Assembling output...\n")
phase_start <- Sys.time()

# Filter to non-NULL results
results_list <- results_list[!vapply(results_list, is.null, logical(1))]
cat("  Gages with results:", length(results_list), "\n")

# Get union of all column names
all_cols <- unique(unlist(lapply(results_list, names)))

# Build data.table row by row
out_dt <- rbindlist(lapply(results_list, function(sigs) {
  # Ensure all columns present
  row <- vector("list", length(all_cols))
  names(row) <- all_cols
  for (col in all_cols) {
    v <- sigs[[col]]
    row[[col]] <- if (is.null(v)) NA else v
  }
  as.data.table(row)
}), fill = TRUE)

# Order columns: metadata first, then signatures (sorted), then flags
meta_cols <- intersect(c("gage_id", "latitude", "longitude", "basin_area", "gage_type",
                         "num_water_years", "start_water_year", "end_water_year",
                         "area_normalized"), names(out_dt))
flag_cols <- sort(grep("^flagged_", names(out_dt), value = TRUE))
sig_cols <- sort(setdiff(names(out_dt), c(meta_cols, flag_cols)))

col_order <- c(meta_cols, sig_cols, flag_cols)
col_order <- intersect(col_order, names(out_dt))
out_dt <- out_dt[, col_order, with = FALSE]

# Save
fwrite(out_dt, output_path)
cat("  Saved to", output_path, "\n")
cat("  Rows:", nrow(out_dt), " Cols:", ncol(out_dt), "\n")

timing$phases$assemble <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))

# Finalize timing
end_time <- Sys.time()
timing$end_time <- format(end_time, "%Y-%m-%dT%H:%M:%S")
timing$total_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
timing$n_gages_processed <- nrow(out_dt)
timing$n_columns <- ncol(out_dt)
timing$n_signature_columns <- length(sig_cols)
timing$n_metadata_columns <- length(meta_cols)
timing$n_qaqc_flags <- length(flag_cols)

# Save timing
timing_path <- "docs/benchmarks/rpkg_timing.json"
jsonlite::write_json(timing, timing_path, auto_unbox = TRUE, pretty = TRUE)
cat("  Timing saved to", timing_path, "\n")

cat("\n======================================================================\n")
cat("BENCHMARK COMPLETE\n")
cat("======================================================================\n")
cat("Total time:", round(timing$total_seconds, 1), "s (",
    round(timing$total_seconds / 60, 1), " min)\n")
cat("Gages processed:", timing$n_gages_processed, "\n")
cat("Gages skipped:", skipped, "\n")
cat("Gages errored:", errored, "\n")
cat("Total columns:", timing$n_columns, "\n")
cat("  Signature columns:", timing$n_signature_columns, "\n")
cat("  Metadata columns:", timing$n_metadata_columns, "\n")
cat("  QA/QC flag columns:", timing$n_qaqc_flags, "\n")
cat("Rate:", round(timing$n_gages_processed / timing$total_seconds, 2), "gages/s\n")

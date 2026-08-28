#!/usr/bin/env Rscript

# Print warnings as they occur. R defaults to deferring them to the end of
# the script and truncating at 50, which hides per-gage signature-family
# failures and leaves check_signature_failures.py unable to certify the log.
options(warn = 1)
# ==============================================================================
# rpkg Benchmark - Streamflow Signatures
# ==============================================================================
# Runs the rpkg package (streamflowsignatures) on all qualifying gages with
# climate data, producing output in the same format as run_python_benchmark.py
# and run_julia_benchmark.jl (the canonical runner).
#
# Rebuilt 2026-08-24 for the port campaign (Phase 0 —
# docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md; Codex
# findings F8/F9/F14): ENV overrides mirroring the Julia runner, runtime
# window + qualifying-fraction gates, column-projected climate read, loads the
# INSTALLED package by default (devtools::load_all conceals installed-package
# failures — opt back in with STREAMFLOW_RPKG_LOADALL=1), legacy-filtering
# fail-fast, a provenance block in the timing JSON, and a nonzero exit status
# when any gage errors (silent whole-gage drops are a gate failure).
#
# Output: {STREAMFLOW_OUTPUT_DIR}/{STREAMFLOW_OUTPUT_PREFIX}_signatures.csv
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

env_or <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (v == "") default else v
}
env_int <- function(name) {
  v <- Sys.getenv(name, unset = "")
  if (v == "") NULL else as.integer(v)
}
env_num <- function(name) {
  v <- Sys.getenv(name, unset = "")
  if (v == "") NULL else as.numeric(v)
}

# Experiment controls (read at runtime; ENV names mirror the Julia runner)
start_water_year    <- env_int("STREAMFLOW_START_WATER_YEAR")
end_water_year      <- env_int("STREAMFLOW_END_WATER_YEAR")
min_qualifying_frac <- env_num("STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION")

parquet_path  <- env_or("STREAMFLOW_DATA_PATH", "D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet")
daymet_path   <- env_or("STREAMFLOW_CLIMATE_PATH", "D:/processedOuts_feb2026/daymet_1980_2023.parquet")
metadata_path <- env_or("STREAMFLOW_METADATA_PATH", "D:/processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv")
output_dir    <- env_or("STREAMFLOW_OUTPUT_DIR", "docs/benchmarks")
output_prefix <- env_or("STREAMFLOW_OUTPUT_PREFIX", "rpkg")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_dir, paste0(output_prefix, "_signatures.csv"))
timing_path <- file.path(output_dir, paste0(output_prefix, "_timing.json"))

cat("Experiment: OUTPUT_PREFIX=", output_prefix,
    " START_WY=", ifelse(is.null(start_water_year), "NULL", start_water_year),
    " END_WY=", ifelse(is.null(end_water_year), "NULL", end_water_year),
    " MIN_QUAL_FRAC=", ifelse(is.null(min_qualifying_frac), "NULL", min_qualifying_frac),
    "\nOutput dir: ", output_dir, "\n\n", sep = "")

# Phase 1: Load rpkg — INSTALLED package by default. devtools::load_all() can
# conceal installed-package failures (missing NAMESPACE exports, stale bundled
# config); only use it for local iteration.
cat("Phase 1: Loading rpkg...\n")
phase_start <- Sys.time()

if (Sys.getenv("STREAMFLOW_RPKG_LOADALL") == "1") {
  cat("  STREAMFLOW_RPKG_LOADALL=1 -> devtools::load_all (dev mode)\n")
  library(devtools)
  load_all("rpkg")
} else {
  library(streamflowsignatures)
}
library(data.table)

# pkg_env is internal (not exported): reach it via ::: so the runner works with
# the INSTALLED package, not only under devtools::load_all()
pkg_env <- streamflowsignatures:::pkg_env

timing$phases$load_config <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat("  rpkg loaded in", round(timing$phases$load_config, 1), "s\n")

# Legacy filtering bypasses the modern preprocessor plumbing — a
# reference/port-campaign run on that path compares apples to oranges.
use_legacy <- pkg_env$use_legacy_filtering
if (isTRUE(use_legacy) && Sys.getenv("STREAMFLOW_ALLOW_LEGACY") != "1") {
  stop("config has na_handling.use_legacy_filtering=true. This runner requires ",
       "the modern preprocessing path; set STREAMFLOW_ALLOW_LEGACY=1 only for ",
       "deliberate legacy-compat experiments.")
}

# Phase 2: Load data
cat("\nPhase 2: Loading data...\n")
phase_start <- Sys.time()

if (!file.exists(parquet_path)) stop("Parquet not found: ", parquet_path)

cat("  Reading streamflow parquet...\n")
# Memory: read ONLY the columns this runner uses. The parquet stores 9
# (gage_id, Date, Q, year, month, doy, water_year, dowy, flag) but the loop
# calls add_water_year_columns() right after subsetting, which re-derives
# water_year/month/dowy from date; flag/doy/year are never read. Loading all
# 9 for 111.6M rows put ~4 GB of dead weight resident and drove this 16 GB
# machine into swap (measured 2026-08-26: 163 MB/s page-ins, CPU at 54%).
# `flag` is deliberately NOT loaded. rpkg's preprocessor can count ice-affected
# NA days from it, but canonical Julia emits ice_affected_days_total = 0 for
# ALL 6,678 gages — including 01011000 / 01029200 / 15129500, which provably
# have 'P Ice' rows with NULL Q. Loading flag here would make rpkg emit
# non-zero values and DIVERGE from canonical. Julia is canonical: match its
# output, and fix the canonical side first (logged in CHANGELOG).
streamflow <- data.table::as.data.table(arrow::read_parquet(
  parquet_path, col_select = c("gage_id", "Date", "Q")))
data.table::setnames(streamflow, "Date", "date")
cat("    ", nrow(streamflow), "rows,", length(unique(streamflow$gage_id)), "gages\n")

cat("  Reading Daymet parquet (column-projected)...\n")
has_daymet <- file.exists(daymet_path)
if (has_daymet) {
  # Memory: read only the columns this runner uses (the full 8-column ~98M-row
  # frame is several GB of dead weight) — PPT and SWE for the climate/snow
  # signatures. Raw names cover pre-/post-normalization variants.
  sch_names <- arrow::open_dataset(daymet_path)$schema$names
  raw_wanted <- c("gage_id", "site_id", "site_no", "date", "Date",
                  "PPT", "prcp", "PRCP", "SWE", "swe")
  sel <- intersect(raw_wanted, sch_names)
  daymet <- data.table::as.data.table(arrow::read_parquet(daymet_path, col_select = sel))
  # Same normalization map as rpkg::read_parquet
  nmap <- c(Date = "date", site_id = "gage_id", prcp = "PPT",
            site_no = "gage_id", PRCP = "PPT", swe = "SWE")
  for (old in names(nmap)) {
    if (old %in% names(daymet) && !nmap[[old]] %in% names(daymet)) {
      data.table::setnames(daymet, old, nmap[[old]])
    }
  }
  daymet_gages <- unique(daymet$gage_id)
  cat("    ", nrow(daymet), "rows,", length(daymet_gages),
      "gages with climate data (columns:", paste(names(daymet), collapse = ","), ")\n")
} else {
  warning("Daymet not found, running without climate")
  daymet <- NULL
  daymet_gages <- character(0)
}

cat("  Reading metadata...\n")
metadata <- fread(metadata_path)
cat("    ", nrow(metadata), "gages in metadata\n")

# Area-normalization lookup (Q-to-PPT gate): gages with no drainage area carry
# Q in raw m3/s — Q-to-PPT signatures (runoff ratios, elasticity, Q-P
# seasonality, storage) are skipped for them because Q and PPT units don't
# match. Keys are leading-zero-stripped for numeric ids (metadata stores USGS
# ids without leading zeros). Only an explicit FALSE gates; NA defaults TRUE.
normalize_join_id <- function(id) {
  id <- as.character(id)
  is_num <- grepl("^[0-9]+$", id)
  stripped <- sub("^0+", "", id)
  ifelse(is_num, ifelse(stripped == "", "0", stripped), id)
}
area_norm_lookup <- logical(0)
if ("area_normalized" %in% names(metadata)) {
  an <- metadata$area_normalized
  # Only an explicit false-like value gates (logical FALSE, numeric 0,
  # "FALSE"/"0" strings — harmonized with the Julia/Python runners); NA
  # defaults to normalized
  explicit_false <- if (is.logical(an)) {
    !is.na(an) & !an
  } else if (is.numeric(an)) {
    !is.na(an) & an == 0
  } else {
    !is.na(an) & toupper(trimws(as.character(an))) %in% c("FALSE", "0")
  }
  area_norm_lookup <- stats::setNames(!explicit_false, normalize_join_id(metadata$gage_id))
  cat("    Area-normalization lookup:", length(area_norm_lookup), "gages,",
      sum(!area_norm_lookup), "not normalized (Q-to-PPT signatures will be skipped)\n")
}

timing$phases$load_data <- as.numeric(difftime(Sys.time(), phase_start, units = "secs"))
cat("  Data loaded in", round(timing$phases$load_data, 1), "s\n")

# Phase 3: Per-gage filtering and signature extraction (single pass)
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

# Changepoint configuration (mirrors the Julia/Python runners)
cp_config <- NULL
if (isTRUE(pkg_env$changepoint_enabled)) {
  cp_config <- list(start_year = pkg_env$cp_start_water_year,
                    end_year = pkg_env$cp_end_water_year,
                    min_total_obs = pkg_env$cp_min_total_obs,
                    min_segment_obs = pkg_env$cp_min_segment_obs)
}
save_annual <- isTRUE(pkg_env$save_annual_values)
cat("Changepoint: enabled=", isTRUE(pkg_env$changepoint_enabled),
    " WY=", pkg_env$cp_start_water_year, "-", pkg_env$cp_end_water_year, "\n", sep = "")
cat("Stats floor: min_values_for_stats=",
    if (is.null(pkg_env$min_values_for_stats)) "NULL" else pkg_env$min_values_for_stats,
    " (recession/elasticity exempt)\n", sep = "")
cat("Annual values: save=", save_annual, "\n", sep = "")
annual_chunks <- list()

results_list <- vector("list", n_gages)
names(results_list) <- all_gages
processed <- 0
skipped <- 0
errored <- 0
error_gages <- character(0)

# Progress tracking
report_interval <- 50

for (i in seq_along(all_gages)) {
  gage <- all_gages[i]

  tryCatch({
    # Extract gage data
    rows <- streamflow_idx[[gage]]
    gage_data <- streamflow[rows, ]
    gage_data <- add_water_year_columns(gage_data)

    # Q-to-PPT gate: default TRUE when gage absent from lookup
    gan <- area_norm_lookup[normalize_join_id(gage)]
    gage_area_normalized <- is.na(gan) || gan

    # Merge climate data if available (needed before preprocessing).
    # PPT and SWE are both carried through; the preprocessor derives
    # valid_swe_years from the SWE column.
    has_climate <- FALSE
    if (gage %in% names(daymet_idx)) {
      daymet_rows <- daymet_idx[[gage]]
      gage_climate <- daymet[daymet_rows, ]
      clim_cols <- intersect(c("date", "PPT", "SWE"), names(gage_climate))
      if ("PPT" %in% clim_cols) {
        climate_merge <- gage_climate[, clim_cols, with = FALSE]
        gage_data <- merge(gage_data, climate_merge, by = "date", all.x = TRUE)
        has_climate <- TRUE
      }
    }

    if (use_legacy) {
      # Legacy path (deprecated compat mode; blocked above unless explicitly
      # allowed): filter qualifying years
      filt <- filter_qualifying_years(gage_data)
      if (!filt$qualifies) {
        skipped <- skipped + 1
        next
      }

      gage_filtered <- gage_data[gage_data$water_year %in% filt$qualifying_years, ]

      sigs <- calculate_all_signatures(gage_filtered, has_climate = has_climate,
                                       area_normalized = gage_area_normalized)

      sigs$gage_id <- gage
      sigs$start_water_year <- min(filt$qualifying_years)
      sigs$end_water_year <- max(filt$qualifying_years)
      sigs$num_water_years <- length(filt$qualifying_years)
    } else {
      # New path: preprocess, then window + fraction gates (semantics mirror the
      # Julia runner: preprocess the FULL record, then filter the RESULT — never
      # pre-filter raw rows; interpolation and rejection must see the same
      # series in every language)
      pp <- preprocess_daily_data(gage_data)

      if (!is.null(start_water_year) || !is.null(end_water_year)) {
        wy_lo <- if (is.null(start_water_year)) -Inf else start_water_year
        wy_hi <- if (is.null(end_water_year)) Inf else end_water_year
        pp$data <- pp$data[pp$data$water_year >= wy_lo & pp$data$water_year <= wy_hi, ]
        pp$valid_years <- pp$valid_years[pp$valid_years >= wy_lo & pp$valid_years <= wy_hi]
        pp$valid_climate_years <- pp$valid_climate_years[
          pp$valid_climate_years >= wy_lo & pp$valid_climate_years <= wy_hi]
        if (!is.null(pp$valid_swe_years)) {  # NULL only when the preprocessor saw no SWE
          pp$valid_swe_years <- pp$valid_swe_years[
            pp$valid_swe_years >= wy_lo & pp$valid_swe_years <= wy_hi]
        }
        # seasonal_flags / diagnostics are data.tables in rpkg (not row lists)
        pp$seasonal_flags <- pp$seasonal_flags[
          pp$seasonal_flags$water_year >= wy_lo & pp$seasonal_flags$water_year <= wy_hi, ]
        pp$diagnostics <- pp$diagnostics[
          pp$diagnostics$water_year >= wy_lo & pp$diagnostics$water_year <= wy_hi, ]
      }

      # Min qualifying data fraction gate. Denominator = per-gage possible years
      # within the [start, end] window: window start .. the gage's last
      # in-window year (window-capped) — identical to the Julia runner.
      if (!is.null(min_qualifying_frac)) {
        all_wy <- unique(gage_data$water_year)
        all_wy <- all_wy[!is.na(all_wy)]
        start_yr <- if (!is.null(start_water_year)) start_water_year else min(all_wy)
        end_yr <- if (!is.null(end_water_year)) end_water_year else Inf
        wy_in_range <- all_wy[all_wy >= start_yr & all_wy <= end_yr]
        if (length(wy_in_range) > 0) {
          total_possible <- max(wy_in_range) - start_yr + 1
          frac <- length(pp$valid_years) / total_possible
          if (frac < min_qualifying_frac) {
            skipped <- skipped + 1
            next
          }
        }
      }

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

      # Snow data: EXPLICIT opt-in frame filtered to SWE-valid years. Possibly
      # 0 rows; NULL when the gage has no SWE column at all. No implicit
      # fallback inside the orchestrator (mirrors the Julia/Python runners).
      snow_data <- NULL
      if ("SWE" %in% names(pp_data)) {
        snow_data <- pp_data[pp_data$water_year %in% pp$valid_swe_years, , drop = FALSE]
      }

      gage_collector <- if (save_annual) annual_collector() else NULL

      sigs <- calculate_all_signatures(
        pp_data, has_climate = gage_has_climate,
        seasonal_flags = pp$seasonal_flags,
        trend_completeness = pkg_env$na_trend_min_fraction,
        decade_completeness = pkg_env$na_decade_min_fraction,
        min_values_for_stats = pkg_env$min_values_for_stats,
        climate_data = climate_data_subset,
        area_normalized = gage_area_normalized,
        changepoint = cp_config,
        collector = gage_collector,
        snow_data = snow_data,
        snow_climate_years = pp$valid_climate_years,
        gage_id = gage
      )

      # Drain this gage's annual values
      if (!is.null(gage_collector)) {
        d <- collector_drain(gage_collector)
        if (nrow(d) > 0) {
          d$gage_id <- gage
          annual_chunks[[length(annual_chunks) + 1L]] <- d
        }
      }

      sigs$gage_id <- gage
      sigs$start_water_year <- min(pp$valid_years)
      sigs$end_water_year <- max(pp$valid_years)
      sigs$num_water_years <- length(pp$valid_years)

      # Per-gage ice diagnostic (mirrors run_julia_benchmark.jl)
      sigs$ice_affected_days_total <- if (!is.null(pp$diagnostics) &&
                                          "na_cause_ice" %in% names(pp$diagnostics)) {
        sum(pp$diagnostics$na_cause_ice, na.rm = TRUE)
      } else 0
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

    # QA/QC flags. compute_qa_flags ALREADY returns names prefixed with
    # "flagged_for_" — the previous `paste0("flagged_for_", nm)` produced
    # `flagged_for_flagged_for_*` for all 12 flags. Pre-existing; invisible to
    # the intersection-based comparison scripts, caught 2026-08-25 by the
    # strict schema gate (check_schema_equality.py) on the rpkg baseline.
    qa <- compute_qa_flags(sigs)
    for (nm in names(qa)) {
      sigs[[nm]] <- qa[[nm]]
    }

    results_list[[gage]] <- sigs
    processed <- processed + 1

  }, error = function(e) {
    errored <<- errored + 1
    error_gages <<- c(error_gages, gage)
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

# Phase 4a: write the annual values parquet (before the metadata merge, so a
# metadata failure cannot lose it). Long format, deterministically sorted.
if (save_annual && length(annual_chunks) > 0) {
  cat("\nWriting annual values parquet...\n")
  annual_dt <- data.table::rbindlist(annual_chunks)
  data.table::setcolorder(annual_dt, c("gage_id", "signature", "water_year", "value"))
  data.table::setorder(annual_dt, gage_id, signature, water_year)
  annual_path <- file.path(output_dir, paste0(output_prefix, "_signatures_annual.parquet"))
  arrow::write_parquet(annual_dt, annual_path, compression = "zstd")
  cat("  Saved", nrow(annual_dt), "rows (",
      length(unique(annual_dt$signature)), "signatures ) to", annual_path, "\n")
  timing$n_annual_rows <- nrow(annual_dt)
  rm(annual_dt); invisible(gc())
}

# Phase 4b: interference metadata (GAGES-II + Canadian HYDAT), mirroring the
# Julia runner's Phase 5. Closes a pre-existing rpkg schema gap (2026-08-24):
# the rpkg output previously lacked all 11 interference columns.
gages_dir_env <- Sys.getenv("STREAMFLOW_GAGES_II_DIR", unset = "")
gages_ii <- tryCatch({
  if (gages_dir_env == "") load_gages_ii_interference() else load_gages_ii_interference(gages_dir_env)
}, error = function(e) { cat("  Warning: GAGES-II load failed:", e$message, "\n"); NULL })
if (!is.null(gages_ii) && nrow(gages_ii) > 0 && "STAID" %in% names(gages_ii)) {
  gages_ii <- data.table::as.data.table(gages_ii)
  interference_cols <- c("NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
                         "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
                         "HYDRO_DISTURB_INDX", "CLASS")
  avail_int_cols <- intersect(interference_cols, names(gages_ii))
  gsub_dt <- gages_ii[, c("STAID", avail_int_cols), with = FALSE]
  data.table::setnames(gsub_dt, "STAID", "gage_id")
  gsub_dt[, gage_id := as.character(gage_id)]
  out_dt[, gage_id := as.character(gage_id)]
  out_dt <- merge(out_dt, gsub_dt, by = "gage_id", all.x = TRUE, sort = FALSE)
  cat("  Matched", sum(!is.na(out_dt[[avail_int_cols[1]]])), "gages with GAGES-II data\n")
  cls <- trimws(as.character(out_dt$CLASS))
  out_dt[, human_interference_class := ifelse(!is.na(cls) & cls == "Ref", "reference",
                                       ifelse(!is.na(cls) & cls == "Non-ref", "non-reference",
                                              "unknown"))]
  # Canadian HYDAT interference from the pre-exported CSV (same file the Julia
  # runner reads); hic derivation mirrors julia/src/metadata.jl (RHBN true ->
  # reference, false -> non-reference, else unknown)
  out_dt[, RHBN := NA]
  out_dt[, REGULATED := NA]
  hydat_csv <- Sys.getenv("STREAMFLOW_HYDAT_PATH",
                          unset = file.path(project_root, "metadata", "canadian_hydat_interference.csv"))
  if (file.exists(hydat_csv)) {
    can <- fread(hydat_csv)
    if ("STATION_NUMBER" %in% names(can)) {
      data.table::setnames(can, "STATION_NUMBER", "gage_id")
      can[, gage_id := as.character(gage_id)]
      can <- unique(can, by = "gage_id", fromLast = TRUE)
      m <- match(out_dt$gage_id, can$gage_id)
      hit <- !is.na(m)
      out_dt$RHBN[hit] <- can$RHBN[m[hit]]
      out_dt$REGULATED[hit] <- can$REGULATED[m[hit]]
      rh <- can$RHBN[m[hit]]
      out_dt$human_interference_class[hit] <-
        ifelse(!is.na(rh) & rh, "reference",
        ifelse(!is.na(rh) & !rh, "non-reference", "unknown"))
      cat("  Matched", sum(hit), "Canadian gages with HYDAT metadata\n")
    }
  } else {
    cat("  Warning: Canadian HYDAT CSV not found:", hydat_csv, "\n")
  }
} else {
  cat("  Warning: GAGES-II data not available\n")
}

# Order columns: metadata first, then signatures (sorted), then flags
meta_cols <- intersect(c("gage_id", "latitude", "longitude", "basin_area", "gage_type",
                         "num_water_years", "start_water_year", "end_water_year",
                         "area_normalized",
                         "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
                         "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
                         "HYDRO_DISTURB_INDX", "CLASS", "RHBN", "REGULATED",
                         "human_interference_class"), names(out_dt))
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
timing$n_gages_skipped <- skipped
timing$n_gages_errored <- errored
timing$error_gages <- error_gages
timing$n_columns <- ncol(out_dt)
timing$n_signature_columns <- length(sig_cols)
timing$n_metadata_columns <- length(meta_cols)
timing$n_qaqc_flags <- length(flag_cols)

# Provenance (mirrors the Julia runner's block; Codex MAJOR-4 / port-plan §4)
file_info <- function(path, hash_big) {
  if (!file.exists(path)) return(list(path = path, exists = FALSE))
  info <- list(path = normalizePath(path), exists = TRUE,
               bytes = file.size(path),
               mtime = format(file.mtime(path), "%Y-%m-%dT%H:%M:%S"))
  if (hash_big || file.size(path) < 50e6) {
    info$sha256 <- digest::digest(path, algo = "sha256", file = TRUE)
  }
  info
}
hash_big <- Sys.getenv("STREAMFLOW_HASH_INPUTS") == "1"
git_rev <- tryCatch(system2("git", c("-C", shQuote(project_root), "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = FALSE)[1],
                    error = function(e) "unavailable")
git_dirty <- tryCatch(length(system2("git", c("-C", shQuote(project_root), "status", "--porcelain"),
                                     stdout = TRUE, stderr = FALSE)) > 0,
                      error = function(e) NA)
config_file <- Sys.getenv("STREAMFLOW_SIGNATURES_CONFIG",
                          unset = system.file("config", "signatures_config.json",
                                              package = "streamflowsignatures"))
env_keys <- c("STREAMFLOW_START_WATER_YEAR", "STREAMFLOW_END_WATER_YEAR",
              "STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION", "STREAMFLOW_SIGNATURES_CONFIG",
              "STREAMFLOW_OUTPUT_PREFIX", "STREAMFLOW_ALLOW_LEGACY",
              "STREAMFLOW_RPKG_LOADALL")
env_overrides <- Filter(function(x) x != "", stats::setNames(Sys.getenv(env_keys), env_keys))
pkg_names <- c("streamflowsignatures", "data.table", "arrow", "zyp", "Kendall")
pkg_versions <- lapply(stats::setNames(pkg_names, pkg_names), function(p)
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) "unavailable"))
timing$provenance <- list(
  streamflow = file_info(parquet_path, hash_big),
  climate = file_info(daymet_path, hash_big),
  metadata = file_info(metadata_path, hash_big),
  config = file_info(config_file, hash_big),
  git_revision = git_rev,
  git_working_tree_dirty = git_dirty,
  r_version = as.character(getRversion()),
  package_versions = pkg_versions,
  hostname = Sys.info()[["nodename"]],
  env_overrides = as.list(env_overrides)
)

# Save timing
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

# Zero-error gate (Codex F9): any errored gage means silently missing rows —
# fail loudly so a validation pipeline can't mistake this run for clean.
if (errored > 0) {
  cat("\nERROR:", errored, "gages errored (first gages: ",
      paste(utils::head(error_gages, 10), collapse = ", "),
      ") — exiting nonzero.\n")
  quit(status = 1)
}

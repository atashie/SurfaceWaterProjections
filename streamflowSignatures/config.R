# Streamflow Signatures - Configuration File
# Centralized parameters for data processing and signature extraction
# Source this file before running processing scripts
# ==============================================================================
# USER-CONFIGURABLE PATHS (edit these for your system, or set environment vars)
# ==============================================================================
# These paths point to external data directories. Set via environment variables
# or edit the defaults below to match your local setup.
#
# Environment variables (optional, override defaults):
#   STREAMFLOW_PARQUET_DIR - Directory containing processed parquet files
#   STREAMFLOW_CARAVAN_DIR - Directory containing Caravan NetCDF data
#   STREAMFLOW_METADATA_DIR - Directory containing gage metadata files
#   STREAMFLOW_DOWNLOADS_DIR - Directory for downloaded files

PARQUET_DATA_DIR <- Sys.getenv("STREAMFLOW_PARQUET_DIR",
                               unset = "D:/processedOuts_feb2026")

CARAVAN_DATA_DIR <- Sys.getenv("STREAMFLOW_CARAVAN_DIR",
                               unset = "D:/Caravan-nc/usr/local/google/home/kratzert/Data/Caravan-Jan25-nc")

METADATA_DATA_DIR <- Sys.getenv("STREAMFLOW_METADATA_DIR",
                                unset = "D:/")

DOWNLOADS_DIR <- Sys.getenv("STREAMFLOW_DOWNLOADS_DIR",
                            unset = file.path(Sys.getenv("USERPROFILE", "~"), "Downloads"))


# ==============================================================================
# DATA QUALITY THRESHOLDS
# ==============================================================================

# Minimum flow threshold and required days
# First value: minimum Q value in mm/day (values below are set to NA)
# Second value: minimum days per water year that must exceed the threshold
MIN_Q_VALUE_AND_DAYS <- c(0.0001, 30)

# Minimum number of "good" water years required for a gage to be processed
MIN_NUM_YEARS <- 20

# Minimum fraction of non-NA days per water year to count as "good"
MIN_FRAC_GOOD_DATA <- 0.95

# Minimum non-NA days per year for annual calculations (floor(365 * 0.68))
MIN_NONA_DAYS_ANNUAL <- 250

# ==============================================================================
# DATE RANGES
# ==============================================================================

# Default analysis period
DEFAULT_START_DATE <- as.Date("1979-10-01")  # Start of water year 1980
DEFAULT_END_DATE <- as.Date("2025-10-01")    # Current

# Alternative shorter period (modern data)
MODERN_START_DATE <- as.Date("1973-01-01")

# ==============================================================================
# WATER YEAR CONFIGURATION
# ==============================================================================
# Northern Hemisphere: Water year starts October 1
# Southern Hemisphere: Would start April 1 (not currently implemented)

WATER_YEAR_START_MONTH <- 10  # October

# ==============================================================================
# BASEFLOW FILTER PARAMETERS
# ==============================================================================

# Eckhardt filter parameters
ECKHARDT_BFIMAX <- 0.8   # Maximum baseflow index (range: 0.25-0.80)
ECKHARDT_ALPHA <- 0.98   # Recession constant (range: 0.90-0.99)

# Lyne-Hollick filter parameters
LYNE_HOLLICK_ALPHA <- 0.925  # Filter parameter
LYNE_HOLLICK_PASSES <- 2     # Number of filter passes (forward and backward)

# ==============================================================================
# RECESSION ANALYSIS PARAMETERS
# ==============================================================================

# Minimum consecutive days for recession event
RECESSION_MIN_DAYS <- 5

# Minimum number of recession events for analysis
RECESSION_MIN_EVENTS <- 25

# ==============================================================================
# FDC (FLOW DURATION CURVE) PARAMETERS
# ==============================================================================

# Small constant added to flow values before log10 transform to handle zeros
FDC_FLOW_FLOOR <- 1e-10

# ==============================================================================
# PULSE ANALYSIS PARAMETERS
# ==============================================================================

# Percentile thresholds for pulse identification
HIGH_PULSE_PERCENTILE <- 0.90  # 90th percentile
LOW_PULSE_PERCENTILE <- 0.10   # 10th percentile

# Minimum threshold for flow reversal detection (fraction of current flow)
FLOW_REVERSAL_THRESHOLD <- 0.02  # 2%

# ==============================================================================
# USGS DATA QUALITY CODES
# ==============================================================================
# Codes that indicate acceptable data quality from USGS NWIS

USGS_ACCEPTED_CODES <- c("A", "A e", "P", "P e")
# A = Approved, P = Provisional, e = estimated

# ==============================================================================
# PARQUET PROCESSING
# ==============================================================================

# Number of gages per chunk when saving to parquet
PARQUET_CHUNK_SIZE <- 1000

# Number of gages per batch during signature processing
PROCESSING_BATCH_SIZE <- 500

# Progress reporting interval (number of gages)
PROGRESS_INTERVAL <- 500

# ==============================================================================
# DAYMET CLIMATE DATA CONFIGURATION
# ==============================================================================

# Path to Daymet ZIP file (source data)
DAYMET_ZIP_PATH <- "data_out/daymet_1980_2023.zip"

# Path to converted parquet file (processed data with Date column)
DAYMET_PARQUET_PATH <- "D:/processedOuts_feb2026/daymet_1980_2023.parquet"

# Daymet date range
DAYMET_START_YEAR <- 1980
DAYMET_END_YEAR <- 2023

# ==============================================================================
# HUMAN INTERFERENCE METADATA (GAGES-II)
# ==============================================================================

# Path to GAGES-II metadata directory
GAGES_II_DIR <- Sys.getenv("GAGES_II_DIR", unset = "D:/gagesMetadata")

# GAGES-II files to load (CONUS)
GAGES_II_FILES_CONUS <- list(
  hydromod_dams = "conterm_hydromod_dams.txt",
  pop_infrastr = "conterm_pop_infrastr.txt",
  hydromod_other = "conterm_hydromod_other.txt",
  bas_classif = "conterm_bas_classif.txt",
  lc06_basin = "conterm_lc06_basin.txt"
)

# GAGES-II files to load (Alaska/Hawaii/Puerto Rico)
GAGES_II_FILES_AKHIPR <- list(
  hydromod_dams = "AKHIPR_hydromod_dams.txt",
  pop_infrastr = "AKHIPR_pop_infrastr.txt",
  hydromod_other = "AKHIPR_hydromod_other.txt",
  bas_classif = "AKHIPR_bas_classif.txt"
  # Note: AKHIPR does not have lc06_basin file
)

# Columns to extract from each GAGES-II file type
GAGES_II_COLUMNS <- list(
  hydromod_dams = c("STAID", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009"),
  pop_infrastr = c("STAID", "IMPNLCD06"),
  hydromod_other = c("STAID", "FRESHW_WITHDRAWAL"),
  bas_classif = c("STAID", "CLASS", "HYDRO_DISTURB_INDX"),
  lc06_basin = c("STAID", "DEVNLCD06")
)

# Missing value sentinel in GAGES-II files
GAGES_II_MISSING_VALUE <- -999

# Human interference metadata columns (for schema validation)
EXPECTED_INTERFERENCE_COLS <- c(
  "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
  "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
  "HYDRO_DISTURB_INDX", "CLASS",
  "RHBN", "REGULATED",
  "human_interference_class"
)

# ==============================================================================
# CLIMATE SIGNATURE PARAMETERS
# ==============================================================================

# Streamflow elasticity parameters (Sawicz et al. 2011)
ELASTICITY_WINDOW_YEARS <- 11      # Rolling window size for elasticity trends
MIN_YEARS_ELASTICITY <- 15         # Minimum years required for elasticity calculation
ELASTICITY_MIN_DATA_COMPLETENESS <- 0.9  # Minimum fraction of valid days per year (90%)
ELASTICITY_MIN_ANNUAL_PPT <- 10    # Minimum annual precipitation in mm

# Q-P seasonality parameters (Wrede et al. 2015)
QP_SLOPE_WINDOW_DAYS <- 30         # Rolling window for cumulative Q-P slope calculation
QP_MIN_YEARS <- 10                 # Minimum years required for Q-P seasonality calculation

# Average storage parameters (Peters & Aulenbach 2011)
# Note: Simplified water balance uses P - Q only (no ET estimation)

# ==============================================================================
# OUTPUT CONFIGURATION
# ==============================================================================

# Expected output columns for CSV schema validation
# (This helps ensure backward compatibility)
EXPECTED_METADATA_COLS <- c(
  "gage_id", "gage_id_metadata", "latitude", "longitude", "basin_area",
  "gage_type", "area_normalized",
  "num_water_years", "start_water_year", "end_water_year",
  # Human interference metadata (optional - may be NA for gages without GAGES-II data)
  "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
  "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
  "HYDRO_DISTURB_INDX", "CLASS", "RHBN", "REGULATED",
  "human_interference_class"
)

# Signature column suffixes (8 statistics per metric)
STAT_SUFFIXES <- c(
  "_senn_slp",      # Theil-Sen slope (robust non-parametric trend)
  "_linear_slp",    # Linear regression slope (parametric trend)
  "_spearman_rho",  # Spearman's rank correlation coefficient
  "_spearman_pval", # Spearman's p-value for trend significance
  "_mk_rho",        # Mann-Kendall tau (non-parametric trend correlation)
  "_mk_pval",       # Mann-Kendall p-value
  "_mean",          # Arithmetic mean across water years
  "_median"         # Median across water years
)

# ==============================================================================
# LOGGING SYSTEM
# ==============================================================================

# Log levels (lower number = more verbose)
LOG_LEVELS <- list(
  DEBUG = 10,
  INFO = 20,
  WARN = 30,
  ERROR = 40,
  NONE = 100
)

# Current log level (set to INFO by default, DEBUG for development)
LOG_LEVEL <- LOG_LEVELS$INFO

# Optional log file path (NULL = console only)
LOG_FILE <- NULL

# Logging functions
log_message <- function(level, level_name, ..., context = NULL) {
  if (level >= LOG_LEVEL) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    msg_parts <- list(...)
    msg <- paste(sapply(msg_parts, as.character), collapse = " ")

    # Add context if provided
    if (!is.null(context)) {
      msg <- paste0("[", context, "] ", msg)
    }

    formatted_msg <- sprintf("[%s] %s: %s", timestamp, level_name, msg)

    # Output to console
    if (level_name == "ERROR") {
      message(formatted_msg)  # Use message() for errors (goes to stderr)
    } else {
      cat(formatted_msg, "\n")
    }

    # Output to file if configured
    if (!is.null(LOG_FILE)) {
      tryCatch({
        cat(formatted_msg, "\n", file = LOG_FILE, append = TRUE)
      }, error = function(e) {
        # Silently fail if log file write fails
      })
    }
  }
  invisible(NULL)
}

log_debug <- function(..., context = NULL) {
  log_message(LOG_LEVELS$DEBUG, "DEBUG", ..., context = context)
}

log_info <- function(..., context = NULL) {
  log_message(LOG_LEVELS$INFO, "INFO", ..., context = context)
}

log_warn <- function(..., context = NULL) {
  log_message(LOG_LEVELS$WARN, "WARN", ..., context = context)
}

log_error <- function(..., context = NULL) {
  log_message(LOG_LEVELS$ERROR, "ERROR", ..., context = context)
}

# Set log level helper
set_log_level <- function(level) {
  if (is.character(level)) {
    level <- toupper(level)
    if (level %in% names(LOG_LEVELS)) {
      LOG_LEVEL <<- LOG_LEVELS[[level]]
      log_info("Log level set to", level)
    } else {
      log_warn("Unknown log level:", level, ". Valid levels:",
               paste(names(LOG_LEVELS), collapse = ", "))
    }
  } else if (is.numeric(level)) {
    LOG_LEVEL <<- level
  }
  invisible(LOG_LEVEL)
}

# Set log file helper
set_log_file <- function(file_path) {
  if (is.null(file_path)) {
    LOG_FILE <<- NULL
    log_info("File logging disabled")
  } else {
    # Create directory if needed
    dir_path <- dirname(file_path)
    if (!dir.exists(dir_path) && dir_path != ".") {
      dir.create(dir_path, recursive = TRUE)
    }
    LOG_FILE <<- file_path
    log_info("Logging to file:", file_path)
  }
  invisible(LOG_FILE)
}

# ==============================================================================
# INPUT VALIDATION FUNCTIONS
# ==============================================================================

# Validate file exists and is readable
validate_file_exists <- function(file_path, param_name = "file",
                                  required_ext = NULL, context = NULL) {
  if (is.null(file_path) || !is.character(file_path) || nchar(file_path) == 0) {
    log_error(param_name, "path is NULL or empty", context = context)
    stop(paste0(param_name, " path must be a non-empty string"))
  }

  if (!file.exists(file_path)) {
    log_error(param_name, "not found:", file_path, context = context)
    stop(paste0(param_name, " not found: ", file_path))
  }

  if (!is.null(required_ext)) {
    file_ext <- tolower(tools::file_ext(file_path))
    valid_exts <- tolower(required_ext)
    if (!file_ext %in% valid_exts) {
      log_error(param_name, "has invalid extension:", file_ext,
                ". Expected:", paste(required_ext, collapse = ", "), context = context)
      stop(paste0(param_name, " must have extension: ",
                  paste(required_ext, collapse = " or ")))
    }
  }

  log_debug("Validated file:", file_path, context = context)
  invisible(TRUE)
}

# Validate directory exists (optionally create)
validate_directory <- function(dir_path, param_name = "directory",
                                create = FALSE, context = NULL) {
  if (is.null(dir_path) || !is.character(dir_path) || nchar(dir_path) == 0) {
    log_error(param_name, "path is NULL or empty", context = context)
    stop(paste0(param_name, " path must be a non-empty string"))
  }

  if (!dir.exists(dir_path)) {
    if (create) {
      tryCatch({
        dir.create(dir_path, recursive = TRUE)
        log_info("Created directory:", dir_path, context = context)
      }, error = function(e) {
        log_error("Failed to create", param_name, ":", e$message, context = context)
        stop(paste0("Failed to create ", param_name, ": ", e$message))
      })
    } else {
      log_error(param_name, "not found:", dir_path, context = context)
      stop(paste0(param_name, " not found: ", dir_path))
    }
  }

  log_debug("Validated directory:", dir_path, context = context)
  invisible(TRUE)
}

# Validate numeric parameter within range
validate_numeric <- function(value, param_name, min_val = NULL, max_val = NULL,
                              allow_null = FALSE, context = NULL) {
  if (is.null(value)) {
    if (allow_null) return(invisible(TRUE))
    log_error(param_name, "is NULL but required", context = context)
    stop(paste0(param_name, " cannot be NULL"))
  }

  if (!is.numeric(value)) {
    log_error(param_name, "must be numeric, got:", class(value), context = context)
    stop(paste0(param_name, " must be numeric"))
  }

  if (!is.null(min_val) && any(value < min_val)) {
    log_error(param_name, "=", value, "is below minimum:", min_val, context = context)
    stop(paste0(param_name, " must be >= ", min_val))
  }

  if (!is.null(max_val) && any(value > max_val)) {
    log_error(param_name, "=", value, "is above maximum:", max_val, context = context)
    stop(paste0(param_name, " must be <= ", max_val))
  }

  log_debug("Validated", param_name, "=", paste(value, collapse = ", "), context = context)
  invisible(TRUE)
}

# Validate data frame/data.table has required columns
validate_columns <- function(df, required_cols, df_name = "data", context = NULL) {
  if (is.null(df)) {
    log_error(df_name, "is NULL", context = context)
    stop(paste0(df_name, " cannot be NULL"))
  }

  if (!is.data.frame(df)) {
    log_error(df_name, "must be a data.frame, got:", class(df), context = context)
    stop(paste0(df_name, " must be a data.frame or data.table"))
  }

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    log_error(df_name, "missing required columns:",
              paste(missing_cols, collapse = ", "), context = context)
    stop(paste0(df_name, " missing required columns: ",
                paste(missing_cols, collapse = ", ")))
  }

  log_debug("Validated", df_name, "has columns:",
            paste(required_cols, collapse = ", "), context = context)
  invisible(TRUE)
}

# Validate date parameter
validate_date <- function(value, param_name, allow_null = FALSE, context = NULL) {
  if (is.null(value)) {
    if (allow_null) return(invisible(TRUE))
    log_error(param_name, "is NULL but required", context = context)
    stop(paste0(param_name, " cannot be NULL"))
  }

  if (!inherits(value, "Date")) {
    # Try to convert
    tryCatch({
      value <- as.Date(value)
    }, error = function(e) {
      log_error(param_name, "cannot be converted to Date:", value, context = context)
      stop(paste0(param_name, " must be a valid Date"))
    })
  }

  log_debug("Validated date", param_name, "=", as.character(value), context = context)
  invisible(TRUE)
}

# Validate gage type
validate_gage_type <- function(gage_type, context = NULL) {
  valid_types <- c("USGS", "Canada", "CANADIAN", "Caravan")
  if (!gage_type %in% valid_types) {
    log_error("Invalid gage_type:", gage_type,
              ". Valid types:", paste(valid_types, collapse = ", "), context = context)
    stop(paste0("gage_type must be one of: ", paste(valid_types, collapse = ", ")))
  }
  log_debug("Validated gage_type:", gage_type, context = context)
  invisible(TRUE)
}

# ==============================================================================
# OUTPUT SCHEMA VALIDATION
# ==============================================================================
#
# REQUIREMENT: Every signature produces exactly 8 statistics via generate_stats():
#   _senn_slp     = Theil-Sen slope (robust non-parametric trend)
#   _linear_slp   = Linear regression slope (parametric trend)
#   _spearman_rho = Spearman's rank correlation coefficient
#   _spearman_pval= Spearman's p-value for trend significance
#   _mk_rho       = Mann-Kendall tau (non-parametric trend correlation)
#   _mk_pval      = Mann-Kendall p-value for trend significance
#   _mean         = Arithmetic mean across water years
#   _median       = Median across water years
#
# EXCEPTIONS (explicitly listed below):
#   - elasticity_static: Single value, not a time series
#   - Recession seasonality columns: Derived from model fitting
#
# To add a new signature:
#   1. Add base name to EXPECTED_SIGNATURE_BASES
#   2. Implement function that returns generate_stats() output
#   3. Run smoke test to verify schema validation
# ==============================================================================

# Define expected signature columns (base names without suffixes)
# Updated to match actual output from helperFunctions.R
EXPECTED_SIGNATURE_BASES <- c(
  # Flow volumes (seasonal and annual)
  "Qann", "Qwin", "Qspr", "Qsum", "Qfal",
  # Flow percentiles (Q1 through Q99)
  "Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40", "Q50",
  "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99",
  "Q95-Q10",  # Ratio of high to low flows
  # Flow Duration Curve metrics
  "FDC90th", "FDCall", "FDCmid",
  # Baseflow indices
  "BFI_Eckhardt", "BFI_LyneHollick",
  # Recession parameters
  "log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity",
  # Pulse metrics (yearly and overall)
  "n_high_pulses_year", "n_low_pulses_year", "n_high_pulses_all", "n_low_pulses_all",
  "dur_high_pulses_year", "dur_low_pulses_year", "dur_high_pulses_all", "dur_low_pulses_all",
  "TQmean",
  # Flow reversals (annual and seasonal)
  "Flow_Reversals_annual", "Flow_Reversals_winter", "Flow_Reversals_spring",
  "Flow_Reversals_summer", "Flow_Reversals_fall",
  # Flashiness
  "flashinessRB",
  # Flow timing (day of water year for cumulative flow percentiles)
  "D5_day", "D10_day", "D20_day", "D30_day", "D40_day", "D50_day",
  "D60_day", "D70_day", "D80_day", "D90_day", "D95_day",
  "D25_to_D75", "Dmax",
  # Q-PPT runoff ratios (climate-dependent)
  "annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
  "summer_runoff_ratio", "fall_runoff_ratio",
  # Streamflow elasticity (Sawicz et al. 2011)
  "elasticity",
  # Q-P seasonality (Wrede et al. 2015)
  "qp_slope_sd", "qp_bimodality",
  # Average storage (Peters & Aulenbach 2011)
  "avg_storage"
)

# Static elasticity column (doesn't follow standard suffix pattern)
EXPECTED_ELASTICITY_STATIC <- "elasticity_static"

# Recession seasonality columns (not standard suffix pattern)
EXPECTED_RECESSION_SEASONALITY <- c(
  "log_a_seasonality_amplitude_all", "log_a_seasonality_amplitude_first_half",
  "log_a_seasonality_amplitude_last_half", "log_a_seasonality_minimum_all",
  "log_a_seasonality_minimum_first_half", "log_a_seasonality_minimum_last_half"
)

# Validate output CSV schema
validate_output_schema <- function(output_df, strict = FALSE, context = NULL) {
  log_info("Validating output schema...", context = context)

  # Check metadata columns
  missing_meta <- setdiff(EXPECTED_METADATA_COLS, names(output_df))
  if (length(missing_meta) > 0) {
    log_warn("Missing metadata columns:", paste(missing_meta, collapse = ", "),
             context = context)
    if (strict) {
      stop(paste0("Output missing required metadata columns: ",
                  paste(missing_meta, collapse = ", ")))
    }
  }

  # Check for expected signature columns (with suffixes)
  expected_sig_cols <- unlist(lapply(EXPECTED_SIGNATURE_BASES, function(base) {
    paste0(base, STAT_SUFFIXES)
  }))

  # Add recession seasonality columns (don't follow standard suffix pattern)
  expected_sig_cols <- c(expected_sig_cols, EXPECTED_RECESSION_SEASONALITY)

  present_sig_cols <- intersect(expected_sig_cols, names(output_df))
  missing_sig_cols <- setdiff(expected_sig_cols, names(output_df))

  log_info("Found", length(present_sig_cols), "of", length(expected_sig_cols),
           "expected signature columns", context = context)

  if (length(missing_sig_cols) > 0) {
    log_debug("Missing signature columns:",
              paste(head(missing_sig_cols, 10), collapse = ", "),
              if(length(missing_sig_cols) > 10) paste0("... and ",
                 length(missing_sig_cols) - 10, " more"),
              context = context)
  }

  # Check for unexpected columns (potential schema drift)
  all_expected <- c(EXPECTED_METADATA_COLS, expected_sig_cols)
  unexpected_cols <- setdiff(names(output_df), all_expected)
  if (length(unexpected_cols) > 0) {
    log_info("Found", length(unexpected_cols), "additional columns not in base schema",
             context = context)
    log_debug("Additional columns:", paste(head(unexpected_cols, 10), collapse = ", "),
              context = context)
  }

  # Return validation summary
  validation_result <- list(
    valid = length(missing_meta) == 0,
    n_metadata_cols = length(intersect(EXPECTED_METADATA_COLS, names(output_df))),
    n_signature_cols = length(present_sig_cols),
    n_missing_meta = length(missing_meta),
    n_missing_sig = length(missing_sig_cols),
    n_unexpected = length(unexpected_cols),
    missing_meta_cols = missing_meta,
    unexpected_cols = unexpected_cols
  )

  log_info("Schema validation complete. Valid:", validation_result$valid,
           context = context)

  return(validation_result)
}

# Validate a single gage's output before adding to summary
validate_gage_output <- function(gage_row, gage_id, context = NULL) {
  # Check for all-NA signature values (indicates calculation failure)
  sig_cols <- grep("_(senn_slp|linear_slp|spearman_rho|spearman_pval|mk_rho|mk_pval|mean|median)$", names(gage_row), value = TRUE)

  if (length(sig_cols) == 0) {
    log_warn("No signature columns found for gage:", gage_id, context = context)
    return(FALSE)
  }

  # Count NA values
  na_count <- sum(is.na(gage_row[, ..sig_cols]))
  total_count <- length(sig_cols)
  na_fraction <- na_count / total_count

  if (na_fraction > 0.5) {
    log_warn("Gage", gage_id, "has", round(na_fraction * 100, 1),
             "% NA values in signatures", context = context)
  }

  # Check for infinite values
  numeric_cols <- names(gage_row)[sapply(gage_row, is.numeric)]
  inf_count <- sum(sapply(gage_row[, ..numeric_cols], function(x) sum(is.infinite(x))))
  if (inf_count > 0) {
    log_warn("Gage", gage_id, "has", inf_count, "infinite values", context = context)
  }

  log_debug("Validated output for gage:", gage_id,
            "- NA:", na_count, "/", total_count, context = context)

  return(TRUE)
}

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Function to print current configuration
print_config <- function() {
  cat("Streamflow Signatures Configuration\n")
  cat("====================================\n")
  cat("\nData Quality:\n")
  cat("  MIN_Q_VALUE_AND_DAYS:", MIN_Q_VALUE_AND_DAYS, "\n")
  cat("  MIN_NUM_YEARS:", MIN_NUM_YEARS, "\n")
  cat("  MIN_FRAC_GOOD_DATA:", MIN_FRAC_GOOD_DATA, "\n")
  cat("\nDate Range:\n")
  cat("  DEFAULT_START_DATE:", as.character(DEFAULT_START_DATE), "\n")
  cat("  DEFAULT_END_DATE:", as.character(DEFAULT_END_DATE), "\n")
  cat("\nBaseflow Parameters:\n")
  cat("  ECKHARDT_BFIMAX:", ECKHARDT_BFIMAX, "\n")
  cat("  ECKHARDT_ALPHA:", ECKHARDT_ALPHA, "\n")
  cat("  LYNE_HOLLICK_ALPHA:", LYNE_HOLLICK_ALPHA, "\n")
  cat("\nHuman Interference Metadata:\n")
  cat("  GAGES_II_DIR:", GAGES_II_DIR, "\n")
  cat("  GAGES_II_DIR exists:", dir.exists(GAGES_II_DIR), "\n")
  cat("\nLogging:\n")
  cat("  LOG_LEVEL:", names(LOG_LEVELS)[LOG_LEVELS == LOG_LEVEL], "\n")
  cat("  LOG_FILE:", ifelse(is.null(LOG_FILE), "None (console only)", LOG_FILE), "\n")
}

# Test file for climate-dependent streamflow signatures
# Run with: source("R/tests/test_climate_signatures.R")

# Load dependencies
library(data.table)
library(arrow)

# Get the project root directory
test_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
if (is.null(test_dir) || test_dir == "") {
  test_dir <- getwd()
}
project_dir <- dirname(test_dir)
# Navigate up from R/tests/ to project root
if (basename(project_dir) == "tests" || basename(project_dir) == "R") {
  project_dir <- dirname(project_dir)
}
if (basename(project_dir) == "R") {
  project_dir <- dirname(project_dir)
}

# Source config and helper functions
source(file.path(project_dir, "config.R"))
source(file.path(project_dir, "R", "helperFunctions.R"))

# ==============================================================================
# TEST UTILITIES
# ==============================================================================

test_count <- 0
pass_count <- 0
fail_count <- 0

run_test <- function(test_name, test_expr) {
  test_count <<- test_count + 1
  cat(sprintf("\n[TEST %d] %s\n", test_count, test_name))

  result <- tryCatch({
    eval(test_expr)
    pass_count <<- pass_count + 1
    cat("  PASS\n")
    TRUE
  }, error = function(e) {
    fail_count <<- fail_count + 1
    cat(sprintf("  FAIL: %s\n", e$message))
    FALSE
  })

  return(result)
}

assert_true <- function(condition, message = "Condition should be TRUE") {
  if (!isTRUE(condition)) {
    stop(message)
  }
}

assert_equal <- function(actual, expected, message = NULL) {
  if (!identical(actual, expected)) {
    msg <- if (is.null(message)) {
      sprintf("Expected %s but got %s", deparse(expected), deparse(actual))
    } else {
      message
    }
    stop(msg)
  }
}

assert_columns_exist <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing columns: %s", paste(missing, collapse = ", ")))
  }
}

# ==============================================================================
# TEST: convert_daymet_zip_to_parquet
# ==============================================================================

test_convert_daymet_zip_to_parquet <- function() {
  cat("\n========== Testing convert_daymet_zip_to_parquet ==========\n")

  # Test 1: Function exists
  run_test("convert_daymet_zip_to_parquet function exists", {
    assert_true(exists("convert_daymet_zip_to_parquet"),
                "Function convert_daymet_zip_to_parquet should exist")
  })

  # Test 2: Output has correct columns
  run_test("Output parquet has required columns", {
    # Use a small subset for testing (just 1980)
    temp_parquet <- tempfile(fileext = ".parquet")
    zip_path <- file.path(project_dir, DAYMET_ZIP_PATH)

    if (!file.exists(zip_path)) {
      stop("Daymet ZIP file not found - skipping test")
    }

    # Convert just one year for speed
    result <- convert_daymet_zip_to_parquet(
      daymet_zip_path = zip_path,
      output_parquet_path = temp_parquet,
      years = 1980,
      max_sites = 5  # Limit to 5 sites for fast testing
    )

    # Read back and check
    df <- read_parquet(temp_parquet)
    expected_cols <- c("site_id", "Date", "prcp", "tmin", "tmax", "swe", "vp", "srad")
    assert_columns_exist(df, expected_cols)

    # Cleanup
    unlink(temp_parquet)
  })

  # Test 3: Date column is properly reconstructed
  run_test("Date column is properly reconstructed", {
    temp_parquet <- tempfile(fileext = ".parquet")
    zip_path <- file.path(project_dir, DAYMET_ZIP_PATH)

    if (!file.exists(zip_path)) {
      stop("Daymet ZIP file not found - skipping test")
    }

    result <- convert_daymet_zip_to_parquet(
      daymet_zip_path = zip_path,
      output_parquet_path = temp_parquet,
      years = 1980,
      max_sites = 1
    )

    df <- as.data.table(read_parquet(temp_parquet))

    # Check that we have 365 or 366 days for the site
    n_days <- nrow(df)
    assert_true(n_days %in% c(365, 366),
                sprintf("Expected 365 or 366 days, got %d", n_days))

    # Check date range
    min_date <- min(df$Date)
    max_date <- max(df$Date)
    assert_true(format(min_date, "%Y-%m-%d") == "1980-01-01",
                sprintf("Expected min date 1980-01-01, got %s", min_date))
    assert_true(format(max_date, "%Y-%m-%d") == "1980-12-31",
                sprintf("Expected max date 1980-12-31, got %s", max_date))

    # Cleanup
    unlink(temp_parquet)
  })

  # Test 4: No duplicate site_id + Date combinations
  run_test("No duplicate site_id + Date combinations", {
    temp_parquet <- tempfile(fileext = ".parquet")
    zip_path <- file.path(project_dir, DAYMET_ZIP_PATH)

    if (!file.exists(zip_path)) {
      stop("Daymet ZIP file not found - skipping test")
    }

    result <- convert_daymet_zip_to_parquet(
      daymet_zip_path = zip_path,
      output_parquet_path = temp_parquet,
      years = 1980,
      max_sites = 3
    )

    df <- as.data.table(read_parquet(temp_parquet))

    # Check for duplicates
    n_rows <- nrow(df)
    n_unique <- nrow(unique(df[, .(site_id, Date)]))
    assert_equal(n_rows, n_unique,
                 sprintf("Found %d duplicate site_id + Date combinations", n_rows - n_unique))

    # Cleanup
    unlink(temp_parquet)
  })
}

# ==============================================================================
# TEST: load_daymet_for_gage
# ==============================================================================

test_load_daymet_for_gage <- function() {
  cat("\n========== Testing load_daymet_for_gage ==========\n")

  # Test 1: Function exists
  run_test("load_daymet_for_gage function exists", {
    assert_true(exists("load_daymet_for_gage"),
                "Function load_daymet_for_gage should exist")
  })

  # Test 2: Returns correct structure
  run_test("Returns data.table with correct columns", {
    parquet_path <- file.path(project_dir, DAYMET_PARQUET_PATH)

    if (!file.exists(parquet_path)) {
      stop("Daymet parquet file not found - run conversion first")
    }

    result <- load_daymet_for_gage(
      gage_id = "01011000",
      daymet_parquet_path = parquet_path,
      start_year = 1980,
      end_year = 1981
    )

    assert_true(is.data.table(result), "Result should be a data.table")
    expected_cols <- c("Date", "prcp", "tmin", "tmax", "swe", "vp", "srad")
    assert_columns_exist(result, expected_cols)
  })

  # Test 3: Filters by gage correctly
  run_test("Filters by gage_id correctly", {
    parquet_path <- file.path(project_dir, DAYMET_PARQUET_PATH)

    if (!file.exists(parquet_path)) {
      stop("Daymet parquet file not found - run conversion first")
    }

    result <- load_daymet_for_gage(
      gage_id = "01011000",
      daymet_parquet_path = parquet_path,
      start_year = 1980,
      end_year = 1980
    )

    # Should have ~365 days for one year
    assert_true(nrow(result) %in% c(365, 366),
                sprintf("Expected ~365 rows for one year, got %d", nrow(result)))
  })
}

# ==============================================================================
# TEST: calculate_streamflow_elasticity
# ==============================================================================

test_calculate_streamflow_elasticity <- function() {
  cat("\n========== Testing calculate_streamflow_elasticity ==========\n")

  # Test 1: Function exists
  run_test("calculate_streamflow_elasticity function exists", {
    assert_true(exists("calculate_streamflow_elasticity"),
                "Function calculate_streamflow_elasticity should exist")
  })

  # Test 2: Known elasticity with synthetic data
  run_test("Calculates known elasticity correctly", {
    # Create synthetic data where Q = 0.5 * P (elasticity should be ~1.0)
    years <- 1980:2010
    set.seed(42)
    P_annual <- 1000 + rnorm(length(years), 0, 100)  # ~1000 mm/year
    Q_annual <- 0.5 * P_annual + rnorm(length(years), 0, 10)  # ~50% runoff ratio

    # Create daily data for each year
    daily_data <- rbindlist(lapply(seq_along(years), function(i) {
      yr <- years[i]
      n_days <- ifelse(yr %% 4 == 0, 366, 365)
      data.table(
        water_year = yr,
        Q = Q_annual[i] / n_days,  # Distribute evenly
        prcp = P_annual[i] / n_days,
        month = rep(1:12, each = ceiling(n_days/12))[1:n_days],
        dowy = 1:n_days
      )
    }))

    result <- calculate_streamflow_elasticity(daily_data)

    # Check output columns exist
    assert_columns_exist(result, c("elasticity_static", "elasticity_mean", "elasticity_median"))

    # Elasticity should be close to 1.0 for proportional Q-P relationship
    assert_true(abs(result$elasticity_static - 1.0) < 0.5,
                sprintf("Expected elasticity near 1.0, got %.2f", result$elasticity_static))
  })
}

# ==============================================================================
# TEST: calculate_qp_seasonality
# ==============================================================================

test_calculate_qp_seasonality <- function() {
  cat("\n========== Testing calculate_qp_seasonality ==========\n")

  # Test 1: Function exists
  run_test("calculate_qp_seasonality function exists", {
    assert_true(exists("calculate_qp_seasonality"),
                "Function calculate_qp_seasonality should exist")
  })

  # Test 2: Output has expected columns
  run_test("Output has expected columns", {
    # Create synthetic data with some seasonality
    years <- 1980:2000
    daily_data <- rbindlist(lapply(years, function(yr) {
      n_days <- ifelse(yr %% 4 == 0, 366, 365)
      months <- rep(1:12, each = ceiling(n_days/12))[1:n_days]

      # Seasonal pattern: higher Q/P ratio in winter
      seasonal_factor <- ifelse(months %in% c(12, 1, 2), 0.8, 0.3)
      prcp <- abs(rnorm(n_days, 5, 2))
      Q <- prcp * seasonal_factor + rnorm(n_days, 0, 0.1)
      Q[Q < 0] <- 0.001

      data.table(
        water_year = yr,
        Q = Q,
        prcp = prcp,
        month = months,
        dowy = 1:n_days
      )
    }))

    result <- calculate_qp_seasonality(daily_data)

    expected_cols <- c("qp_slope_sd_mean", "qp_slope_sd_median",
                       "qp_bimodality_mean", "qp_bimodality_median")
    assert_columns_exist(result, expected_cols)
  })
}

# ==============================================================================
# TEST: calculate_average_storage
# ==============================================================================

test_calculate_average_storage <- function() {
  cat("\n========== Testing calculate_average_storage ==========\n")

  # Test 1: Function exists
  run_test("calculate_average_storage function exists", {
    assert_true(exists("calculate_average_storage"),
                "Function calculate_average_storage should exist")
  })

  # Test 2: Output has expected columns
  run_test("Output has expected columns", {
    # Create synthetic data
    years <- 1980:2000
    daily_data <- rbindlist(lapply(years, function(yr) {
      n_days <- ifelse(yr %% 4 == 0, 366, 365)
      prcp <- abs(rnorm(n_days, 3, 2))
      Q <- abs(rnorm(n_days, 2, 1))

      data.table(
        water_year = yr,
        Q = Q,
        prcp = prcp,
        month = rep(1:12, each = ceiling(n_days/12))[1:n_days],
        dowy = 1:n_days
      )
    }))

    result <- calculate_average_storage(daily_data)

    expected_cols <- c("avg_storage_mean", "avg_storage_median",
                       "avg_storage_slp", "avg_storage_rho", "avg_storage_pval")
    assert_columns_exist(result, expected_cols)
  })
}

# ==============================================================================
# RUN ALL TESTS
# ==============================================================================

run_all_tests <- function() {
  cat("\n")
  cat("================================================================\n")
  cat("         CLIMATE SIGNATURES TEST SUITE\n")
  cat("================================================================\n")

  # Reset counters
  test_count <<- 0
  pass_count <<- 0
  fail_count <<- 0

  # Run test suites
  tryCatch(test_convert_daymet_zip_to_parquet(), error = function(e) {
    cat(sprintf("\nTest suite error: %s\n", e$message))
  })

  tryCatch(test_load_daymet_for_gage(), error = function(e) {
    cat(sprintf("\nTest suite error: %s\n", e$message))
  })

  tryCatch(test_calculate_streamflow_elasticity(), error = function(e) {
    cat(sprintf("\nTest suite error: %s\n", e$message))
  })

  tryCatch(test_calculate_qp_seasonality(), error = function(e) {
    cat(sprintf("\nTest suite error: %s\n", e$message))
  })

  tryCatch(test_calculate_average_storage(), error = function(e) {
    cat(sprintf("\nTest suite error: %s\n", e$message))
  })

  # Summary
  cat("\n")
  cat("================================================================\n")
  cat(sprintf("         RESULTS: %d/%d tests passed\n", pass_count, test_count))
  cat("================================================================\n")

  if (fail_count > 0) {
    cat(sprintf("\n  WARNING: %d tests failed!\n", fail_count))
  } else {
    cat("\n  All tests passed!\n")
  }

  return(fail_count == 0)
}

# Run tests if this file is sourced directly
if (interactive()) {
  cat("\nTo run all tests, call: run_all_tests()\n")
} else {
  run_all_tests()
}

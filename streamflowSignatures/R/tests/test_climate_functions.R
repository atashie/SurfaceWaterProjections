# Quick validation of climate signature functions with synthetic data
#
# NOTE: This script expects to be run from the streamflowSignatures directory
# Run with: Rscript R/tests/test_climate_functions.R

cat("Loading packages and functions...\n")
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(lubridate)
})

# Auto-detect project root
main_dir <- getwd()
if (!file.exists(file.path(main_dir, "config.R"))) {
  main_dir <- dirname(main_dir)
  if (!file.exists(file.path(main_dir, "config.R"))) {
    main_dir <- dirname(main_dir)
  }
}

source(file.path(main_dir, "config.R"))
source(file.path(main_dir, "R", "helperFunctions.R"))

cat("\n========== TESTING CLIMATE SIGNATURE FUNCTIONS ==========\n\n")

# Test 1: load_daymet_for_gage
cat("[TEST 1] load_daymet_for_gage\n")
parquet_path <- DAYMET_PARQUET_PATH

if (file.exists(parquet_path)) {
  daymet_data <- load_daymet_for_gage(
    gage_id = "01011000",
    daymet_parquet_path = parquet_path,
    start_year = 1980,
    end_year = 1985
  )
  cat(sprintf("  Loaded %d rows for gage 01011000 (1980-1985)\n", nrow(daymet_data)))
  cat(sprintf("  Columns: %s\n", paste(names(daymet_data), collapse = ", ")))
  cat("  PASS\n\n")
} else {
  cat("  SKIP - parquet file not found at:", parquet_path, "\n\n")
}

# Test 2: calculate_streamflow_elasticity with synthetic data
cat("[TEST 2] calculate_streamflow_elasticity\n")
set.seed(42)
years <- 1980:2010
P_annual <- 1000 + rnorm(length(years), 0, 100)
Q_annual <- 0.5 * P_annual + rnorm(length(years), 0, 10)

synth_data <- rbindlist(lapply(seq_along(years), function(i) {
  yr <- years[i]
  n_days <- ifelse(yr %% 4 == 0, 366, 365)
  data.table(
    water_year = yr,
    Q = Q_annual[i] / n_days,
    PPT = P_annual[i] / n_days,
    month = rep(1:12, each = ceiling(n_days/12))[1:n_days],
    dowy = 1:n_days
  )
}))

elasticity_result <- calculate_streamflow_elasticity(synth_data)
cat(sprintf("  elasticity_static: %.3f\n", elasticity_result$elasticity_static))
cat(sprintf("  elasticity_mean: %.3f\n", elasticity_result$elasticity_mean))
if (abs(elasticity_result$elasticity_static - 1.0) < 0.5) {
  cat("  PASS (elasticity near 1.0 as expected for proportional Q-P)\n\n")
} else {
  cat("  WARNING - elasticity not near expected value\n\n")
}

# Test 3: calculate_qp_seasonality with synthetic seasonal data
cat("[TEST 3] calculate_qp_seasonality\n")
years <- 1980:2000
synth_seasonal <- rbindlist(lapply(years, function(yr) {
  n_days <- ifelse(yr %% 4 == 0, 366, 365)
  months <- rep(1:12, each = ceiling(n_days/12))[1:n_days]
  # Higher Q/P ratio in winter months
  seasonal_factor <- ifelse(months %in% c(12, 1, 2), 0.8, 0.3)
  PPT <- abs(rnorm(n_days, 5, 2))
  Q <- PPT * seasonal_factor + abs(rnorm(n_days, 0, 0.1))
  data.table(
    water_year = yr,
    Q = Q,
    PPT = PPT,
    month = months,
    dowy = 1:n_days
  )
}))

seasonality_result <- calculate_qp_seasonality(synth_seasonal)
cat(sprintf("  qp_slope_sd_mean: %.4f\n", seasonality_result$qp_slope_sd_mean))
cat(sprintf("  qp_bimodality_mean: %.4f\n", seasonality_result$qp_bimodality_mean))
cat("  PASS\n\n")

# Test 4: calculate_average_storage
cat("[TEST 4] calculate_average_storage\n")
years <- 1980:2000
synth_storage <- rbindlist(lapply(years, function(yr) {
  n_days <- ifelse(yr %% 4 == 0, 366, 365)
  PPT <- abs(rnorm(n_days, 3, 2))
  Q <- abs(rnorm(n_days, 2, 1))
  data.table(
    water_year = yr,
    Q = Q,
    PPT = PPT,
    month = rep(1:12, each = ceiling(n_days/12))[1:n_days],
    dowy = 1:n_days
  )
}))

storage_result <- calculate_average_storage(synth_storage)
cat(sprintf("  avg_storage_mean: %.2f mm\n", storage_result$avg_storage_mean))
cat(sprintf("  avg_storage_slp: %.4f\n", storage_result$avg_storage_slp))
cat("  PASS\n\n")

# Test 5: Daymet parquet statistics
cat("[TEST 5] Daymet parquet verification\n")
if (file.exists(parquet_path)) {
  df <- arrow::read_parquet(parquet_path)
  cat(sprintf("  Total rows: %d\n", nrow(df)))
  cat(sprintf("  Unique sites: %d\n", length(unique(df$site_id))))
  cat(sprintf("  Date range: %s to %s\n", min(df$Date), max(df$Date)))
  cat(sprintf("  File size: %.1f MB\n", file.info(parquet_path)$size / 1024 / 1024))
  cat("  PASS\n\n")
} else {
  cat("  SKIP - parquet file not found\n\n")
}

cat("========== ALL TESTS COMPLETED ==========\n")

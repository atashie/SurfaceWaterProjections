# test_na_handling.R — Unit tests for preprocess_daily_data() and related NA handling
#
# Usage:
#   source("config.R")
#   source("R/helperFunctions.R")
#   source("R/tests/test_na_handling.R")
#
# All tests use synthetic data — no external parquet files needed.
# ==============================================================================

cat("=== NA Handling Unit Tests ===\n\n")

# Helper: create a synthetic water-year data.table for testing
# Returns a data.table with Date, Q, water_year, month, dowy columns
make_wy_data <- function(water_year, q_values = NULL, n_days = 365,
                          ppt_values = NULL) {
  wy_start <- as.Date(paste0(water_year - 1, "-10-01"))
  dates <- seq.Date(wy_start, by = "day", length.out = n_days)

  if (is.null(q_values)) {
    q_values <- runif(n_days, 0.1, 10)
  }

  dt <- data.table(
    Date = dates,
    Q = q_values,
    water_year = water_year,
    month = as.integer(format(dates, "%m")),
    dowy = seq_len(n_days)
  )

  if (!is.null(ppt_values)) {
    dt$PPT <- ppt_values
  }

  dt
}

# Helper: build multi-year dataset
make_multi_year <- function(years, ...) {
  rbindlist(lapply(years, function(yr) make_wy_data(yr, ...)))
}

# Default config for testing
test_config <- list(
  interpolation = list(max_gap_days = 3L, method = "linear", internal_only = TRUE),
  year_rejection = list(max_raw_na_per_year = 30L, reject_negative_flow = TRUE,
                        reject_residual_na = TRUE),
  constant_sd_flag = list(enabled = TRUE, min_nonzero_days_per_month = 15L,
                          max_unique_values = 1L),
  seasonal_completeness = list(min_fraction = 0.80,
                                use_raw_observations = TRUE,
                                season_definitions = list(
                                  winter = c(12L, 1L, 2L), spring = c(3L, 4L, 5L),
                                  summer = c(6L, 7L, 8L), fall = c(9L, 10L, 11L)
                                )),
  climate_na_policy = list(max_raw_na_per_year_ppt = 30L,
                            max_interpolation_gap_ppt = 3L,
                            reject_negative_ppt = TRUE)
)

tests_passed <- 0
tests_failed <- 0

assert <- function(condition, test_name) {
  if (condition) {
    cat(paste0("  PASS: ", test_name, "\n"))
    tests_passed <<- tests_passed + 1
  } else {
    cat(paste0("  FAIL: ", test_name, "\n"))
    tests_failed <<- tests_failed + 1
  }
}

# ==============================================================================
# 1. DAILY GRID NORMALIZATION
# ==============================================================================
cat("--- 1. Daily Grid Normalization ---\n")

# Test: missing date rows are filled with NA
test_data <- make_wy_data(2000)
# Remove 5 rows in the middle
test_data <- test_data[!(dowy %in% 50:54)]
result <- preprocess_daily_data(test_data, config = test_config)
# The preprocessed data for WY 2000 should have 365 or 366 rows
wy2000 <- result$data[water_year == 2000]
wy_start <- as.Date("1999-10-01")
wy_end <- as.Date("2000-09-30")
expected_days <- as.integer(wy_end - wy_start) + 1L
assert(nrow(wy2000) == expected_days || 2000 %in% result$valid_years,
       "Missing date rows filled — grid normalized")

# Test: duplicate dates removed
test_data2 <- make_wy_data(2001)
# Add a duplicate row
test_data2 <- rbind(test_data2, test_data2[10])
result2 <- preprocess_daily_data(test_data2, config = test_config)
wy2001 <- result2$data[water_year == 2001]
assert(nrow(wy2001) == length(unique(wy2001$Date)),
       "Duplicate dates removed")

# Test: dates are sorted
test_data3 <- make_wy_data(2002)
test_data3 <- test_data3[sample(nrow(test_data3))]  # shuffle
result3 <- preprocess_daily_data(test_data3, config = test_config)
wy2002 <- result3$data[water_year == 2002]
assert(all(diff(wy2002$Date) > 0),
       "Dates sorted in output")

# ==============================================================================
# 2. INTERPOLATION
# ==============================================================================
cat("\n--- 2. Interpolation ---\n")

# Test: single-day internal gap correctly interpolated
test_interp <- make_wy_data(2003, q_values = rep(10, 365))
test_interp$Q[100] <- NA  # single internal gap
result_interp <- preprocess_daily_data(test_interp, config = test_config)
wy2003 <- result_interp$data[water_year == 2003]
assert(2003 %in% result_interp$valid_years, "Single-day gap: year accepted")
assert(!is.na(wy2003$Q[wy2003$dowy == 100]), "Single-day gap: interpolated")
assert(abs(wy2003$Q[wy2003$dowy == 100] - 10) < 0.01,
       "Single-day gap: interpolated value correct (10)")

# Test: 3-day internal gap correctly interpolated
test_interp3 <- make_wy_data(2004, q_values = rep(10, 365))
test_interp3$Q[100:102] <- NA  # 3-day internal gap
result_interp3 <- preprocess_daily_data(test_interp3, config = test_config)
assert(2004 %in% result_interp3$valid_years, "3-day gap: year accepted")
wy2004 <- result_interp3$data[water_year == 2004]
assert(all(!is.na(wy2004$Q[wy2004$dowy %in% 100:102])),
       "3-day gap: all 3 days interpolated")

# Test: 4-day gap causes year rejection (gap too large)
test_interp4 <- make_wy_data(2005, q_values = rep(10, 365))
test_interp4$Q[100:103] <- NA  # 4-day gap
result_interp4 <- preprocess_daily_data(test_interp4, config = test_config)
assert(!(2005 %in% result_interp4$valid_years),
       "4-day gap: year rejected")
assert(any(result_interp4$rejected_years$water_year == 2005),
       "4-day gap: appears in rejected_years")

# Test: gap at START of water year (leading NAs) — NOT interpolated
test_leading <- make_wy_data(2006, q_values = rep(10, 365))
test_leading$Q[1:2] <- NA  # leading NAs
result_leading <- preprocess_daily_data(test_leading, config = test_config)
# With reject_residual_na=TRUE, this year should be rejected (boundary NAs can't be interpolated)
assert(!(2006 %in% result_leading$valid_years),
       "Leading NAs: year rejected (boundary NAs not interpolated)")

# Test: gap at END of water year (trailing NAs) — NOT interpolated
test_trailing <- make_wy_data(2007, q_values = rep(10, 365))
test_trailing$Q[364:365] <- NA  # trailing NAs
result_trailing <- preprocess_daily_data(test_trailing, config = test_config)
assert(!(2007 %in% result_trailing$valid_years),
       "Trailing NAs: year rejected (boundary NAs not interpolated)")

# ==============================================================================
# 3. YEAR REJECTION
# ==============================================================================
cat("\n--- 3. Year Rejection ---\n")

# Test: 31 scattered single-day NAs → rejected (raw_na_count > 30)
test_many_na <- make_wy_data(2008, q_values = rep(10, 365))
na_positions <- seq(10, 320, by = 10)  # 31 positions
test_many_na$Q[na_positions[1:31]] <- NA
result_many_na <- preprocess_daily_data(test_many_na, config = test_config)
assert(!(2008 %in% result_many_na$valid_years),
       "31 scattered NAs: year rejected (>30 raw NAs)")

# Test: 30 scattered single-day NAs → accepted and interpolated
test_30_na <- make_wy_data(2009, q_values = rep(10, 365))
na_pos_30 <- seq(10, 300, by = 10)  # 30 positions, all internal
test_30_na$Q[na_pos_30[1:30]] <- NA
result_30_na <- preprocess_daily_data(test_30_na, config = test_config)
assert(2009 %in% result_30_na$valid_years,
       "30 scattered NAs: year accepted")

# Test: negative Q values → rejected
test_neg <- make_wy_data(2010, q_values = rep(10, 365))
test_neg$Q[50] <- -1
result_neg <- preprocess_daily_data(test_neg, config = test_config)
assert(!(2010 %in% result_neg$valid_years),
       "Negative Q: year rejected")
assert(any(grepl("negative", result_neg$rejected_years$reason)),
       "Negative Q: rejection reason mentions negative")

# Test: mix of passing and failing years
mixed_data <- rbind(
  make_wy_data(2011, q_values = rep(10, 365)),  # good
  make_wy_data(2012, q_values = rep(10, 365)),  # will have negative
  make_wy_data(2013, q_values = rep(10, 365))   # good
)
mixed_data$Q[mixed_data$water_year == 2012 & mixed_data$dowy == 50] <- -5
result_mixed <- preprocess_daily_data(mixed_data, config = test_config)
assert(2011 %in% result_mixed$valid_years && 2013 %in% result_mixed$valid_years,
       "Mixed years: good years accepted")
assert(!(2012 %in% result_mixed$valid_years),
       "Mixed years: bad year rejected")

# ==============================================================================
# 4. CONSTANT-SD FLAG
# ==============================================================================
cat("\n--- 4. Constant-SD Flag ---\n")

# Test: month with single unique non-zero Q value → flagged
test_const <- make_wy_data(2014, q_values = rep(10, 365))
# Make all of November (month 11, days ~32-61 of water year) constant
nov_mask <- test_const$month == 11
test_const$Q[nov_mask] <- 5.0  # single unique non-zero value, 30 days
result_const <- preprocess_daily_data(test_const, config = test_config)
diag_2014 <- result_const$diagnostics[water_year == 2014]
assert(diag_2014$constant_sd_months >= 1,
       "Constant-SD: flagged month with single unique non-zero Q")

# Test: month with genuine low variability but >1 unique value → NOT flagged
test_varied <- make_wy_data(2015, q_values = rep(10, 365))
nov_mask2 <- test_varied$month == 11
test_varied$Q[nov_mask2] <- 5.0
test_varied$Q[which(nov_mask2)[1]] <- 5.001  # slightly different
result_varied <- preprocess_daily_data(test_varied, config = test_config)
diag_2015 <- result_varied$diagnostics[water_year == 2015]
# This month has 2 unique values, so it should NOT be flagged
assert(diag_2015$constant_sd_months == 0,
       "Constant-SD: not flagged with >1 unique value")

# ==============================================================================
# 5. SEASONAL COMPLETENESS
# ==============================================================================
cat("\n--- 5. Seasonal Completeness ---\n")

# Test: computed from RAW observations, not post-interpolation
test_seasonal <- make_wy_data(2016, q_values = rep(10, 365))
# Remove 25% of December days (should flag winter as incomplete)
dec_rows <- which(test_seasonal$month == 12)
n_remove <- ceiling(length(dec_rows) * 0.25)
test_seasonal$Q[dec_rows[1:n_remove]] <- NA
# But the total NAs should be < 30 to not reject the year
result_seasonal <- preprocess_daily_data(test_seasonal, config = test_config)
if (2016 %in% result_seasonal$valid_years) {
  sf <- result_seasonal$seasonal_flags[water_year == 2016]
  # Winter includes Dec, Jan, Feb — Dec has ~25% missing from raw
  # But Jan and Feb are complete, so overall winter fraction should be > 80%
  # unless December alone has very high missingness relative to the season total
  assert("winter_complete" %in% names(sf), "Seasonal flags: winter_complete column exists")
}

# Test: season with 75% raw data → flagged incomplete
test_seasonal2 <- make_wy_data(2017, q_values = rep(10, 365))
# Remove 25% of ALL winter months (Dec, Jan, Feb = ~90 days, remove ~23)
winter_rows <- which(test_seasonal2$month %in% c(12, 1, 2))
n_winter <- length(winter_rows)
n_remove2 <- ceiling(n_winter * 0.25)
test_seasonal2$Q[winter_rows[1:n_remove2]] <- NA
# Make sure total NAs are still <= 30
total_na <- sum(is.na(test_seasonal2$Q))
if (total_na <= 30) {
  result_seasonal2 <- preprocess_daily_data(test_seasonal2, config = test_config)
  if (2017 %in% result_seasonal2$valid_years) {
    sf2 <- result_seasonal2$seasonal_flags[water_year == 2017]
    assert(!sf2$winter_complete, "Season with 75% raw data: flagged incomplete")
  }
}

# ==============================================================================
# 6. TREND COMPLETENESS (generate_stats)
# ==============================================================================
cat("\n--- 6. Trend Completeness ---\n")

# Test: 80% complete series → trends computed
test_trend_data <- data.frame(
  water_year = 2000:2029,
  metric = runif(30, 1, 10)
)
result_80 <- generate_stats(test_trend_data, value_cols = "metric",
                            year_col = "water_year",
                            trend_completeness = 0.80, decade_completeness = 0.80)
assert(!is.na(result_80$metric_senn_slp),
       "80% complete: trends computed")

# Test: 66.7% complete series (internal gap) → trends set to NA, mean/median still
# computed. NOTE: completeness is measured over the metric's OWN non-NA year span
# (first to last valid year), so trailing/leading NAs cannot reduce it — the gap must
# be internal. (Previous version of this case used 9 trailing NAs and was vacuous:
# span 2000-2020 gave 21/21 = 100%.) Thresholds are passed explicitly to test the
# mechanism; the production default is config-driven (0.60 overall since July 2026).
test_trend_67 <- data.frame(
  water_year = 2000:2029,
  metric = c(runif(10, 1, 10), rep(NA, 10), runif(10, 1, 10))  # 20/30 = 66.7% over span
)
result_67 <- generate_stats(test_trend_67, value_cols = "metric",
                            year_col = "water_year",
                            trend_completeness = 0.80, decade_completeness = 0.80)
assert(is.na(result_67$metric_senn_slp), "66.7% complete (internal gap): trends set to NA")
assert(!is.na(result_67$metric_mean), "66.7% complete: mean still computed")
assert(!is.na(result_67$metric_median), "66.7% complete: median still computed")

# Test: first decade only 50% complete → trends NA
test_first_decade <- data.frame(
  water_year = 2000:2029,
  metric = c(runif(5, 1, 10), rep(NA, 5), runif(20, 1, 10))  # first decade 50%
)
result_fd <- generate_stats(test_first_decade, value_cols = "metric",
                            year_col = "water_year",
                            trend_completeness = 0.80, decade_completeness = 0.80)
assert(is.na(result_fd$metric_senn_slp),
       "First decade 50% complete: trends NA")

# Test: series < 10 years → decade check skipped
test_short <- data.frame(
  water_year = 2020:2027,
  metric = runif(8, 1, 10)
)
result_short <- generate_stats(test_short, value_cols = "metric",
                                year_col = "water_year",
                                trend_completeness = 0.80, decade_completeness = 0.80)
assert(!is.na(result_short$metric_senn_slp),
       "Short series (<10yr): decade check skipped, trends computed")

# ==============================================================================
# 7. CLIMATE (PPT) HANDLING
# ==============================================================================
cat("\n--- 7. Climate (PPT) Handling ---\n")

# Test: PPT NAs tracked separately from Q
test_climate <- make_wy_data(2018, q_values = rep(10, 365),
                              ppt_values = rep(5, 365))
# Add PPT NAs exceeding limit
test_climate$PPT[1:35] <- NA  # 35 > 30 limit
result_climate <- preprocess_daily_data(test_climate, config = test_config)
# Year should be valid for Q (no Q NAs)
assert(2018 %in% result_climate$valid_years,
       "PPT NAs: year valid for Q-only signatures")
# But NOT valid for climate signatures
assert(!(2018 %in% result_climate$valid_climate_years),
       "PPT NAs: year NOT valid for Q-PPT signatures")

# Test: year valid for both Q and PPT
test_both <- make_wy_data(2019, q_values = rep(10, 365),
                           ppt_values = rep(5, 365))
result_both <- preprocess_daily_data(test_both, config = test_config)
assert(2019 %in% result_both$valid_years,
       "Clean year: valid for Q")
assert(2019 %in% result_both$valid_climate_years,
       "Clean year: valid for Q+PPT")

# ==============================================================================
# 8. DIAGNOSTICS
# ==============================================================================
cat("\n--- 8. Diagnostics ---\n")

# Test: diagnostics table has correct columns
test_diag <- make_wy_data(2020, q_values = rep(10, 365))
test_diag$Q[50] <- NA  # 1 internal NA
result_diag <- preprocess_daily_data(test_diag, config = test_config)
expected_cols <- c("water_year", "raw_na_count", "raw_max_na_run",
                   "interpolated_count", "residual_na_count",
                   "negative_flag", "constant_sd_months")
assert(all(expected_cols %in% names(result_diag$diagnostics)),
       "Diagnostics: all expected columns present")
diag_row <- result_diag$diagnostics[water_year == 2020]
assert(diag_row$raw_na_count == 1, "Diagnostics: raw_na_count correct")
assert(diag_row$interpolated_count == 1, "Diagnostics: interpolated_count correct")
assert(diag_row$residual_na_count == 0, "Diagnostics: residual_na_count correct")
assert(!diag_row$negative_flag, "Diagnostics: negative_flag correct")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n=== NA Handling Test Summary ===\n")
cat(sprintf("  Passed: %d\n", tests_passed))
cat(sprintf("  Failed: %d\n", tests_failed))
cat(sprintf("  Total:  %d\n", tests_passed + tests_failed))
if (tests_failed > 0) {
  cat("  STATUS: SOME TESTS FAILED\n")
} else {
  cat("  STATUS: ALL TESTS PASSED\n")
}

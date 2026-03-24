test_that("calculate_all_signatures works on synthetic data", {
  # Generate 25 years of synthetic daily streamflow
  set.seed(123)
  dates <- seq(as.Date("1979-10-01"), as.Date("2004-09-30"), by = "day")
  n <- length(dates)

  # Seasonal pattern with noise
  day_of_year <- as.numeric(format(dates, "%j"))
  seasonal <- sin(2 * pi * day_of_year / 365) * 3 + 5
  Q <- pmax(seasonal + rnorm(n, sd = 1), 0.001)

  df <- data.frame(date = dates, Q = Q)
  df <- add_water_year_columns(df)

  # Non-climate signatures
  result <- calculate_all_signatures(df, has_climate = FALSE)

  expect_true(is.list(result))
  expect_true(length(result) > 100)  # Should have many signature values

  # Check key signatures exist
  expect_true("Qann_mean" %in% names(result))
  expect_true("BFI_Eckhardt_mean" %in% names(result))
  expect_true("flashinessRB_mean" %in% names(result))
  expect_true("D50_day_mean" %in% names(result))
  expect_true("FDCall_mean" %in% names(result))

  # Values should be reasonable
  expect_true(result$Qann_mean > 0)
  expect_true(result$BFI_Eckhardt_mean >= 0 && result$BFI_Eckhardt_mean <= 1)
})

test_that("calculate_all_signatures works with climate data", {
  set.seed(456)
  dates <- seq(as.Date("1979-10-01"), as.Date("2004-09-30"), by = "day")
  n <- length(dates)

  day_of_year <- as.numeric(format(dates, "%j"))
  seasonal <- sin(2 * pi * day_of_year / 365) * 3 + 5
  Q   <- pmax(seasonal + rnorm(n, sd = 1), 0.001)
  PPT <- pmax(rnorm(n, mean = 3, sd = 2), 0)

  df <- data.frame(date = dates, Q = Q, PPT = PPT)
  df <- add_water_year_columns(df)

  result <- calculate_all_signatures(df, has_climate = TRUE)

  expect_true(is.list(result))
  expect_true(length(result) > 200)

  # Climate-specific signatures
  expect_true("annual_runoff_ratio_mean" %in% names(result))
  expect_true("elasticity_static" %in% names(result))
  expect_true("qp_slope_sd_mean" %in% names(result))
  expect_true("avg_storage_mean" %in% names(result))
})

test_that("filter_qualifying_years works correctly", {
  # Create data with some bad years
  set.seed(789)
  dates <- seq(as.Date("1979-10-01"), as.Date("2004-09-30"), by = "day")
  df <- data.frame(date = dates, Q = abs(rnorm(length(dates), mean = 5, sd = 2)))
  df <- add_water_year_columns(df)

  result <- filter_qualifying_years(df)

  expect_true(is.list(result))
  expect_true("qualifying_years" %in% names(result))
  expect_true("qualifies" %in% names(result))
  expect_true(result$qualifies)
  expect_true(length(result$qualifying_years) >= 20)
})

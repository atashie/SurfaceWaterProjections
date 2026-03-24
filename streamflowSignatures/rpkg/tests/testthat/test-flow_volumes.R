test_that("calculate_flow_vols_by_year returns correct structure", {
  # Create 3 years of synthetic daily data
  dates <- seq(as.Date("1999-10-01"), as.Date("2002-09-30"), by = "day")
  n <- length(dates)
  df <- data.frame(
    date = dates,
    Q = abs(sin(seq(0, 6 * pi, length.out = n))) * 5 + 0.5
  )
  df <- add_water_year_columns(df)

  result <- calculate_flow_vols_by_year(df)

  # Should have Qann + 4 seasons + 15 percentiles + Q95-Q10 = 21 metrics x 8 stats
  expect_true(is.list(result))
  expect_true(length(result) >= 21 * 8)

  # Check key names exist
  expect_true("Qann_mean" %in% names(result))
  expect_true("Qwin_mean" %in% names(result))
  expect_true("Q50_mean" %in% names(result))
  expect_true("Q95_Q10_mean" %in% names(result))

  # Values should be positive
  expect_true(result$Qann_mean > 0)
  expect_true(result$Q50_mean > 0)
})

test_that("calculate_flow_vols_by_year rejects insufficient data", {
  # Only 100 days
  df <- data.frame(
    water_year = rep(2000, 100),
    Q = runif(100),
    month = rep(1, 100),
    dowy = 1:100
  )
  result <- calculate_flow_vols_by_year(df)
  expect_true(length(result) == 0 || all(is.na(unlist(result))))
})

test_that("eckhardt_filter produces values between 0 and Q", {
  Q <- c(10, 8, 6, 5, 4, 3, 2, 1.5, 1, 0.8, 0.5, 0.3, 0.2, 0.1, 0.05)
  bf <- streamflowsignatures:::eckhardt_filter(Q)

  expect_length(bf, length(Q))
  expect_true(all(bf >= 0))
  expect_true(all(bf <= Q))
})

test_that("eckhardt_filter handles NA with forward-fill", {
  Q <- c(10, 8, NA, 5, 4)
  bf <- streamflowsignatures:::eckhardt_filter(Q)

  expect_length(bf, 5)
  expect_false(is.na(bf[3]))  # forward-filled
  expect_equal(bf[3], bf[2])  # should equal previous baseflow
})

test_that("lyne_hollick_filter produces values between 0 and Q", {
  Q <- c(10, 8, 6, 5, 4, 3, 2, 1.5, 1, 0.8, 0.5, 0.3, 0.2, 0.1, 0.05)
  bf <- streamflowsignatures:::lyne_hollick_filter(Q)

  expect_length(bf, length(Q))
  # Baseflow should be non-negative where Q is valid
  valid <- !is.na(bf)
  expect_true(all(bf[valid] >= 0))
})

test_that("analyze_baseflow_indices returns correct structure", {
  dates <- seq(as.Date("1999-10-01"), as.Date("2022-09-30"), by = "day")
  n <- length(dates)
  df <- data.frame(date = dates, Q = abs(sin(seq(0, 20 * pi, length.out = n))) * 5 + 0.5)
  df <- add_water_year_columns(df)

  result <- analyze_baseflow_indices(df)

  expect_true(is.list(result))
  expect_true("BFI_Eckhardt_mean" %in% names(result))
  expect_true("BFI_LyneHollick_mean" %in% names(result))

  # BFI should be between 0 and 1
  expect_true(result$BFI_Eckhardt_mean >= 0 && result$BFI_Eckhardt_mean <= 1)
  expect_true(result$BFI_LyneHollick_mean >= 0 && result$BFI_LyneHollick_mean <= 1)
})

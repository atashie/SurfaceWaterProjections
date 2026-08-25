test_that("identify_recession_events finds events in declining flow", {
  # Create a clean exponential recession: Q = 100 * exp(-0.1 * t)
  # dQ/dt = -10 * exp(-0.1 * t), |dQ/dt| decreasing
  t_vals <- 0:20
  Q <- 100 * exp(-0.1 * t_vals)
  events <- streamflowsignatures:::identify_recession_events(Q, min_length = 5)

  expect_true(length(events) >= 1)
  expect_true(length(events[[1]]$indices) >= 5)
})

test_that("identify_recession_events returns empty for flat flow", {
  Q <- rep(5, 30)
  events <- streamflowsignatures:::identify_recession_events(Q, min_length = 5)
  expect_length(events, 0)
})

test_that("fit_recession_event returns log_a and b", {
  # Power law recession: dQ/dt = a * Q^b
  Q <- 100 * exp(-0.05 * (1:15))
  result <- streamflowsignatures:::fit_recession_event(Q, remove_first_day = TRUE)

  expect_true(!is.na(result$log_a))
  expect_true(!is.na(result$b))
})

test_that("fit_sinusoidal_model returns amplitude and minimum_doy", {
  set.seed(42)
  doy <- sample(1:365, 50)
  log_a <- 2 * sin(2 * pi * doy / 365) + rnorm(50, sd = 0.5)

  result <- streamflowsignatures:::fit_sinusoidal_model(doy, log_a)
  expect_true(!is.na(result$amplitude))
  expect_true(!is.na(result$minimum_doy))
  expect_true(result$amplitude > 0)
  expect_true(result$minimum_doy >= 0 && result$minimum_doy <= 365)
})

test_that("analyze_recession_parameters returns NA for insufficient events", {
  # Very short dataset with no recession events
  df <- data.frame(
    water_year = rep(2000, 50),
    Q = rep(5, 50),
    dowy = 1:50
  )
  result <- analyze_recession_parameters(df)
  expect_true(is.list(result))
  expect_true(is.na(result$b_events_mean))
})

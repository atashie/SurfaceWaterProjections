test_that("package exports match expected signature count", {
  # Verify calculate_all_signatures exists and is callable

  expect_true(is.function(calculate_all_signatures))
  expect_true(is.function(generate_stats))
  expect_true(is.function(filter_qualifying_years))
  expect_true(is.function(analyze_fdc_trends))
  expect_true(is.function(analyze_baseflow_indices))
  expect_true(is.function(analyze_recession_parameters))
})

test_that("config loads correctly on package load", {
  # Access internal config to verify it loaded
  env <- streamflowsignatures:::pkg_env

  expect_equal(env$min_num_years, 20)
  expect_equal(env$min_frac_good_data, 0.95)
  expect_equal(env$eckhardt_bfimax, 0.8)
  expect_equal(env$eckhardt_alpha, 0.98)
  expect_equal(env$lyne_hollick_alpha, 0.925)
  expect_equal(env$recession_min_events, 25)
  expect_equal(env$recession_min_length, 5)
  expect_equal(env$elasticity_window_years, 11)
  expect_length(env$d_percentiles, 13)  # D1..D99, 13 since April 2026 (was 11)
  expect_length(env$stat_suffixes, 8)
})

test_that("empty_stats matches stat_suffixes from config", {
  result <- empty_stats("test")
  env <- streamflowsignatures:::pkg_env
  expected_names <- paste0("test", env$stat_suffixes)
  expect_equal(sort(names(result)), sort(expected_names))
})

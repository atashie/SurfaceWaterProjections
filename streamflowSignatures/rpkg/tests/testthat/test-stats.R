test_that("generate_stats produces 8 named values per column", {
  df <- data.frame(water_year = 2000:2020, val = seq(10, 30, length.out = 21))
  result <- generate_stats(df, value_cols = "val", year_col = "water_year")

  expect_true(is.list(result))
  expect_length(result, 8)

  expected_names <- paste0("val", c("_senn_slp", "_linear_slp", "_spearman_rho",
                                    "_spearman_pval", "_mk_rho", "_mk_pval",
                                    "_mean", "_median"))
  expect_equal(sort(names(result)), sort(expected_names))
})

test_that("generate_stats returns all NA for insufficient data", {
  df <- data.frame(water_year = 2000:2001, val = c(1, 2))
  result <- generate_stats(df, value_cols = "val", year_col = "water_year")

  expect_true(all(is.na(unlist(result))))
})

test_that("generate_stats handles NA values in data", {
  df <- data.frame(water_year = 2000:2024,
                   val = c(NA, NA, seq(1, 23)))
  result <- generate_stats(df, value_cols = "val", year_col = "water_year")

  expect_false(is.na(result$val_mean))
  expect_equal(result$val_mean, mean(1:23))
})

test_that("generate_stats produces correct mean and median", {
  df <- data.frame(water_year = 2000:2009, val = 1:10)
  result <- generate_stats(df, value_cols = "val", year_col = "water_year")

  expect_equal(result$val_mean, 5.5)
  expect_equal(result$val_median, 5.5)
})

test_that("generate_stats detects positive trend", {
  df <- data.frame(water_year = 2000:2029, val = seq(1, 30))
  result <- generate_stats(df, value_cols = "val", year_col = "water_year")

  expect_true(result$val_senn_slp > 0)
  expect_true(result$val_linear_slp > 0)
  expect_true(result$val_spearman_rho > 0)
  expect_true(result$val_mk_rho > 0)
})

test_that("empty_stats returns 8 named NAs", {
  result <- empty_stats("test_metric")
  expect_length(result, 8)
  expect_true(all(is.na(unlist(result))))
  expect_true(all(grepl("^test_metric_", names(result))))
})

test_that("generate_stats handles multiple columns", {
  df <- data.frame(water_year = 2000:2019, a = 1:20, b = 20:1)
  result <- generate_stats(df, value_cols = c("a", "b"), year_col = "water_year")

  expect_length(result, 16)
  expect_true(result$a_spearman_rho > 0)  # increasing
  expect_true(result$b_spearman_rho < 0)  # decreasing
})

# Mann-Kendall NA contract parity with Julia/Python (2026-09-04).
# Julia's mann_kendall_test returns (NaN, NaN) for n < 3 or when the tie-corrected
# variance of S is 0 (i.e. a constant series). Kendall::MannKendall instead returns
# tau = 1, sl = 1 with an IFAULT = 12 warning on a constant series.

test_that("constant series gives NA tau and p-value", {
  r <- mann_kendall_test(rep(0, 25))
  expect_true(is.na(r$tau)); expect_true(is.na(r$pval))
  r <- mann_kendall_test(c(3, 3, 3, NA, 3))
  expect_true(is.na(r$tau)); expect_true(is.na(r$pval))
})

test_that("fewer than 3 values gives NA", {
  r <- mann_kendall_test(c(1, 2))
  expect_true(is.na(r$tau)); expect_true(is.na(r$pval))
  r <- mann_kendall_test(c(NA, 1, NA, 2))
  expect_true(is.na(r$tau)); expect_true(is.na(r$pval))
})

test_that("a non-constant series still returns Kendall's values", {
  v <- c(1, 3, 2, 5, 4, 6, 8, 7, 9, 10)
  r <- mann_kendall_test(v)
  k <- Kendall::MannKendall(v)
  expect_equal(r$tau, as.numeric(k$tau))
  expect_equal(r$pval, as.numeric(k$sl))
  expect_true(r$tau > 0.8 && r$pval < 0.01)
})

test_that("generate_stats emits NA MK fields for an all-zero metric (negative_ann shape)", {
  df <- data.frame(water_year = 1990:2019, negative_ann = rep(0, 30))
  res <- generate_stats(df, value_cols = "negative_ann", year_col = "water_year")
  expect_true(is.na(res$negative_ann_mk_rho))
  expect_true(is.na(res$negative_ann_mk_pval))
  expect_equal(res$negative_ann_mean, 0)
  expect_equal(res$negative_ann_median, 0)
})

test_that("segment_differential_metrics gives NA segment p-values on constant segments", {
  years <- 1990:2019
  values <- c(rep(0, 15), rep(1, 15))          # step change, both segments constant
  r <- segment_differential_metrics(years, values, 2004)
  expect_equal(r$delta_mean, 1)
  expect_true(is.na(r$pre_mk_pval)); expect_true(is.na(r$post_mk_pval))
  # a trending segment still gets a p-value
  values2 <- c(rep(0, 15), 1:15)
  r2 <- segment_differential_metrics(years, values2, 2004)
  expect_true(is.na(r2$pre_mk_pval)); expect_false(is.na(r2$post_mk_pval))
})

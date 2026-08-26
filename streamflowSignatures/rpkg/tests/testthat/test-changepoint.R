# Tests for the Pettitt changepoint port (Phase 1a of the port campaign).
#
# The frozen expectations below are the SAME fixture values cross-checked
# against canonical Julia on 2026-08-24 (identical cp_year on every fixture
# including the tie; bit-identical p-values), so this file pins R to the same
# cross-language contract the Python suite pins.

CP_CONFIG <- list(start_year = 1900, end_year = 2100,
                  min_total_obs = 20, min_segment_obs = 10)

step_series <- function() {
  list(years = as.numeric(1990:2014),
       values = c(10 + 0.01 * (0:9), 20 + 0.01 * (0:14)))
}

test_that("clean step: frozen Julia-cross-checked values", {
  s <- step_series()
  r <- pettitt_test(s$years, s$values)
  expect_equal(r$cp_year, 2001)
  expect_equal(r$pval, 0.00025039952406474066, tolerance = 1e-12)
  expect_equal(r$pre_mean, 11.705, tolerance = 1e-12)
  expect_equal(r$post_mean, 20.08, tolerance = 1e-12)

  d <- segment_differential_metrics(s$years, s$values, r$cp_year)
  expect_equal(d$delta_mean, 8.375, tolerance = 1e-12)
  expect_equal(d$pct_change, 71.5506193934216, tolerance = 1e-10)
})

test_that("tied maximum selects the FIRST candidate", {
  # Antisymmetric 'tent': |U| attains its max (30) at TWO candidates
  years <- as.numeric(1980:2002)          # n = 23
  values <- as.numeric(c(1:12, 11:1))
  U <- cumsum(rowSums(sign(outer(values, values, "-"))))
  expect_equal(abs(U[10:13]), c(30, 11, 11, 30))   # the tie is real
  r <- pettitt_test(years, values)
  expect_equal(r$cp_year, 1989)           # first max (t=10), not 1992 (t=13)
  expect_equal(r$pval, 1)
})

test_that("constant series: zero U, first candidate, guards", {
  years <- as.numeric(2000:2024)
  r5 <- pettitt_test(years, rep(5, 25))
  expect_equal(r5$cp_year, 2009)
  expect_equal(r5$pval, 1)
  d5 <- segment_differential_metrics(years, rep(5, 25), r5$cp_year)
  expect_equal(d5$delta_mean, 0)
  expect_equal(d5$pct_change, 0)

  # pre_mean == 0 -> pct_change NA (abs(pre_mean) > 1e-10 guard)
  d0 <- segment_differential_metrics(years, rep(0, 25),
                                     pettitt_test(years, rep(0, 25))$cp_year)
  expect_equal(d0$delta_mean, 0)
  expect_true(is.na(d0$pct_change))
})

test_that("boundary observation counts", {
  # n == min_total_obs == 2 * min_segment_obs: exactly one candidate
  y20 <- as.numeric(1991:2010)
  v20 <- c(rep(1, 10) + 0.001 * (0:9), rep(3, 10) + 0.001 * (0:9))
  r <- pettitt_test(y20, v20)
  expect_equal(r$cp_year, 2000)
  expect_equal(r$pval, 0.0015809806462399323, tolerance = 1e-12)

  # n = 19 < min_total_obs -> all NA
  r19 <- pettitt_test(as.numeric(1991:2009), as.numeric(0:18))
  expect_true(all(vapply(r19, is.na, logical(1))))
})

test_that("NA pairs filtered and input order does not matter", {
  s <- step_series()
  v <- s$values; v[c(4, 8, 21, 26)] <- NA     # R 1-based == Python idx 3,7,20,25
  years28 <- as.numeric(1985:2012)
  vals28 <- c(5 + 0.01 * (0:11), 15 + 0.01 * (0:15))
  vals28[c(4, 8, 21, 26)] <- NA
  r <- pettitt_test(years28, vals28)
  expect_equal(r$cp_year, 1998)
  expect_equal(r$pval, 0.0003537738044851332, tolerance = 1e-12)

  # shuffled input must give an identical result
  set.seed(11)
  perm <- sample(seq_along(s$years))
  expect_equal(pettitt_test(s$years[perm], s$values[perm]), pettitt_test(s$years, s$values))
})

test_that("segment guard: <3 values per side -> all NA", {
  years <- as.numeric(2000:2024)
  vals <- as.numeric(0:24)
  d <- segment_differential_metrics(years, vals, 2001)   # 2 pre values
  expect_true(all(vapply(d, is.na, logical(1))))
  expect_true(all(vapply(segment_differential_metrics(years, vals, NA), is.na, logical(1))))
})

# ---- generate_stats integration ------------------------------------------

step_df <- function(n_pre = 15, n_post = 15, start = 1985) {
  data.frame(water_year = start:(start + n_pre + n_post - 1),
             metric = c(10 + 0.01 * seq_len(n_pre), 20 + 0.01 * seq_len(n_post)))
}

test_that("Pettitt columns appended, absent without config", {
  r <- generate_stats(step_df(), value_cols = "metric", changepoint = CP_CONFIG)
  for (sfx in c("_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
                "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
                "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval")) {
    expect_true(paste0("metric", sfx) %in% names(r), info = sfx)
  }
  expect_false(is.na(r[["metric_pettitt_cp_year"]]))
  expect_lt(r[["metric_pettitt_pval"]], 0.05)

  r_off <- generate_stats(step_df(), value_cols = "metric")
  expect_length(grep("_pettitt_", names(r_off)), 0)
})

test_that("changepoint is independent of trend gating", {
  vals <- ifelse(seq_len(40) %% 2 == 0, as.numeric(seq_len(40)), NA_real_)
  df <- data.frame(water_year = 1980:2019, metric = vals)
  r <- generate_stats(df, value_cols = "metric",
                      trend_completeness = 0.9, decade_completeness = 0.9,
                      changepoint = CP_CONFIG)
  expect_true(is.na(r[["metric_senn_slp"]]))       # trend gated
  expect_false(is.na(r[["metric_mean"]]))          # mean survives
  expect_false(is.na(r[["metric_pettitt_cp_year"]]))  # CP independent
})

test_that("below min_rows emits NA changepoint keys", {
  df <- data.frame(water_year = c(2000, 2001), metric = c(1, 2))
  r <- generate_stats(df, value_cols = "metric", min_rows = 3, changepoint = CP_CONFIG)
  expect_true(is.na(r[["metric_pettitt_cp_year"]]))
  expect_true(is.na(r[["metric_pettitt_post_mk_pval"]]))
})

test_that("frame-level gate emits keys for requested-but-absent metrics", {
  df <- data.frame(water_year = c(2000, 2001), present = c(1, 2))
  r <- generate_stats(df, value_cols = c("present", "absent"), min_rows = 3,
                      changepoint = CP_CONFIG)
  expect_true(all(c("present_mean", "absent_mean",
                    "present_pettitt_pval", "absent_pettitt_pval") %in% names(r)))
  expect_true(all(vapply(r, function(x) is.na(x), logical(1))))
})

test_that("changepoint window filter is applied", {
  df <- step_df(n_pre = 15, n_post = 15, start = 1985)
  cp <- modifyList(CP_CONFIG, list(start_year = 2000, end_year = 2005))
  r <- generate_stats(df, value_cols = "metric", changepoint = cp)
  expect_true(is.na(r[["metric_pettitt_cp_year"]]))   # window too short
})

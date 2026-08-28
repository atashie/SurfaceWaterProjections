# Tests for the stats floor (min_values_for_stats) and force_skip_trends —
# Phase 1b of the port campaign (mirrors python/tests/test_stats_floor.py and
# julia/test/test_stats_floor.jl).
#
# Contract (julia/src/stats.jl equivalent): a metric with fewer non-NA annual
# values than the floor emits NA for ALL 8 statistics AND its changepoint
# fields; force_skip_trends suppresses the 6 trend stats only (mean/median and
# changepoint survive); recession and elasticity are exempt at the
# ORCHESTRATION layer — their signatures do not accept the argument.

CP_CFG <- list(start_year = 1900, end_year = 2100,
               min_total_obs = 20, min_segment_obs = 10)

STAT_SFX <- c("_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
              "_mk_rho", "_mk_pval", "_mean", "_median")
CP_SFX <- c("_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
            "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
            "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval")

series_df <- function(n, start = 1990, name = "metric") {
  df <- data.frame(water_year = start:(start + n - 1),
                   x = 10 + 0.5 * seq_len(n) + (seq_len(n) %% 3) * 0.1)
  names(df)[2] <- name
  df
}

test_that("below the floor: all 8 stats AND changepoint fields are NA", {
  r <- suppressWarnings(generate_stats(series_df(19), value_cols = "metric",
                                       min_values_for_stats = 20, changepoint = CP_CFG))
  for (sfx in c(STAT_SFX, CP_SFX)) {
    expect_true(is.na(r[[paste0("metric", sfx)]]), info = sfx)
  }
})

test_that("at the floor: statistics are computed", {
  r <- suppressWarnings(generate_stats(series_df(20), value_cols = "metric",
                                       min_values_for_stats = 20, changepoint = CP_CFG))
  expect_false(is.na(r[["metric_mean"]]))
  expect_false(is.na(r[["metric_senn_slp"]]))
  expect_false(is.na(r[["metric_pettitt_cp_year"]]))
})

test_that("no floor preserves legacy behavior", {
  a <- suppressWarnings(generate_stats(series_df(5), value_cols = "metric",
                                       min_values_for_stats = NULL))
  b <- suppressWarnings(generate_stats(series_df(5), value_cols = "metric"))
  expect_equal(a, b)
  expect_false(is.na(a[["metric_mean"]]))
})

test_that("floor counts non-NA values, not rows", {
  df <- series_df(25)
  df$metric[1:6] <- NA          # 19 non-NA
  r <- suppressWarnings(generate_stats(df, value_cols = "metric",
                                       min_values_for_stats = 20))
  for (sfx in STAT_SFX) expect_true(is.na(r[[paste0("metric", sfx)]]), info = sfx)
})

test_that("clustered-series regression (the 07292500 shape)", {
  vals <- rep(NA_real_, 40)
  vals[3:6] <- c(5, 6, 4.5, 7)   # 4 clustered non-NA years
  df <- data.frame(water_year = 1980:2019, metric = vals)
  without <- suppressWarnings(generate_stats(df, value_cols = "metric",
                                             trend_completeness = 0.8,
                                             decade_completeness = 0.8))
  expect_false(is.na(without[["metric_senn_slp"]]))   # the bug shape
  with_floor <- suppressWarnings(generate_stats(df, value_cols = "metric",
                                                trend_completeness = 0.8,
                                                decade_completeness = 0.8,
                                                min_values_for_stats = 20))
  for (sfx in STAT_SFX) expect_true(is.na(with_floor[[paste0("metric", sfx)]]), info = sfx)
})

test_that("per-metric independence", {
  df <- series_df(25, name = "dense")
  df$sparse <- ifelse(seq_len(25) <= 10, df$dense, NA_real_)
  r <- suppressWarnings(generate_stats(df, value_cols = c("dense", "sparse"),
                                       min_values_for_stats = 20))
  expect_false(is.na(r[["dense_mean"]]))
  expect_true(is.na(r[["sparse_mean"]]))
})

test_that("recession and elasticity do not accept the floor argument", {
  # Structural exemption, matching Julia: the orchestrator cannot leak it
  expect_false("min_values_for_stats" %in% names(formals(analyze_recession_parameters)))
  expect_false("min_values_for_stats" %in% names(formals(calculate_streamflow_elasticity)))
  # ...but they DO accept changepoint
  expect_true("changepoint" %in% names(formals(analyze_recession_parameters)))
  expect_true("changepoint" %in% names(formals(calculate_streamflow_elasticity)))
})

test_that("force_skip_trends suppresses trends only", {
  r <- suppressWarnings(generate_stats(series_df(25), value_cols = "metric",
                                       force_skip_trends = "metric",
                                       changepoint = CP_CFG))
  for (sfx in c("_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                "_mk_rho", "_mk_pval")) {
    expect_true(is.na(r[[paste0("metric", sfx)]]), info = sfx)
  }
  expect_false(is.na(r[["metric_mean"]]))
  expect_false(is.na(r[["metric_median"]]))
  expect_false(is.na(r[["metric_pettitt_cp_year"]]))   # CP unaffected
})

test_that("force_skip_trends is per-metric and inert when empty", {
  df <- series_df(25, name = "gated")
  df$free <- df$gated + 1
  r <- suppressWarnings(generate_stats(df, value_cols = c("gated", "free"),
                                       force_skip_trends = "gated"))
  expect_true(is.na(r[["gated_senn_slp"]]))
  expect_false(is.na(r[["free_senn_slp"]]))

  plain <- suppressWarnings(generate_stats(series_df(25), value_cols = "metric"))
  none <- suppressWarnings(generate_stats(series_df(25), value_cols = "metric",
                                          force_skip_trends = character(0)))
  expect_equal(none, plain)
})

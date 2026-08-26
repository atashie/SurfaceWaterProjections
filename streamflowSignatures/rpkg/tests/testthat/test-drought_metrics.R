# Streamflow drought (Phase 5). Cross-checked against canonical Julia on
# 2026-08-26: 105 values over three fixtures (seasonal, intermittent, gapped),
# 0 mismatches, worst relative difference 3.6e-14.

mk_dr <- function(f, n_years = 25, skip = NULL) {
  rows <- list()
  for (wy in 1996:(1995 + n_years)) {
    d0 <- as.Date(paste0(wy - 1, "-10-01")); d1 <- as.Date(paste0(wy, "-09-30"))
    nd <- as.integer(d1 - d0) + 1L
    keep <- rep(TRUE, nd)
    if (!is.null(skip)) for (i in 1:nd) keep[i] <- !skip(i, wy)
    df <- data.frame(date = d0 + (0:(nd - 1L)), water_year = wy, dowy = 1:nd,
                     Q = vapply(1:nd, function(i) f(i, wy), numeric(1)))
    rows[[length(rows) + 1L]] <- df[keep, , drop = FALSE]
  }
  do.call(rbind, rows)
}
PCTS <- c(2, 5, 10, 20, 30)
DUR <- paste0("drought_duration_fixed_p", PCTS)
DEF <- paste0("drought_deficit_fixed_p", PCTS)
THR <- paste0("drought_threshold_fixed_p", PCTS)

test_that("weibull_quantile: hand-derived interpolation and range policy", {
  x <- 1:9                                    # n = 9, positions i/10
  expect_equal(weibull_quantile(x, 0.5), 5)   # h = 5 -> the 5th value
  expect_equal(weibull_quantile(x, 0.25), 2.5)
  expect_true(is.na(weibull_quantile(x, 0.05)))               # below 1/(n+1)
  expect_equal(weibull_quantile(x, 0.05, below_range = "clamp"), 1)
  expect_equal(weibull_quantile(x, 0.99), 9)                  # above n/(n+1)
  expect_equal(weibull_quantile(c(1, NA, 2, 3, NA), 0.5), 2)  # NAs dropped
})

test_that("smoothing shrinks at edges and NEVER averages across a gap", {
  d <- as.Date("2000-01-01") + 0:9
  out <- smooth_daily_flow(d, 1:10, window = 7, alignment = "center", min_valid = 1)
  expect_equal(out[4], 4)      # interior: mean of days 1..7
  expect_equal(out[1], 2.5)    # edge: window shrinks to days 1..4

  dg <- as.Date("2000-01-01") + c(0, 1, 2, 3, 4, 6, 7, 8, 9, 10)   # 1-day hole
  og <- smooth_daily_flow(dg, c(rep(1, 5), rep(100, 5)), window = 7,
                          alignment = "center", min_valid = 1)
  expect_equal(og[1:5], rep(1, 5))      # untouched by the 100s across the gap
  expect_equal(og[6:10], rep(100, 5))
})

test_that("smoothing: min_valid and duplicate dates", {
  d2 <- as.Date("2000-01-01") + 0:1
  expect_true(all(is.na(smooth_daily_flow(d2, c(1, 2), window = 7,
                                          alignment = "center", min_valid = 4))))
  ddup <- as.Date("2000-01-01") + c(0, 1, 1, 2)   # duplicate breaks the run
  od <- smooth_daily_flow(ddup, c(1, 2, 3, 4), window = 3,
                          alignment = "center", min_valid = 1)
  expect_equal(od[1], 1.5)
  expect_false(any(is.na(od)))
})

test_that("seasonal gage: all metrics, monotone levels, ordered thresholds", {
  r <- calculate_drought_metrics(
    mk_dr(function(i, wy) 2 + 1.8*sin(2*pi*(i-30)/365) + 0.3*((i*7919) %% 11)/11),
    trend_completeness = 0.6, decade_completeness = 0.8, min_values_for_stats = 20)
  for (m in c(DUR, DEF)) expect_true(paste0(m, "_mean") %in% names(r), info = m)
  durs <- vapply(DUR, function(m) r[[paste0(m, "_mean")]], numeric(1))
  defs <- vapply(DEF, function(m) r[[paste0(m, "_mean")]], numeric(1))
  thr  <- vapply(THR, function(s) r[[s]], numeric(1))
  expect_false(is.unsorted(durs))     # non-decreasing in severity level
  expect_false(is.unsorted(defs))
  expect_false(is.unsorted(thr))
  # weak construction anchor: a level-p threshold sits below ~p% of days
  expect_gt(r$drought_duration_fixed_p10_mean, 30)
  expect_lt(r$drought_duration_fixed_p10_mean, 43)
})

test_that("intermittent gage: zero thresholds + STRICT < give zero duration", {
  r <- calculate_drought_metrics(
    mk_dr(function(i, wy) if (i >= 150 && i <= 300) 0 else 1 + ((i*104729) %% 7)/7),
    trend_completeness = 0.6, decade_completeness = 0.8)
  expect_equal(r$drought_threshold_fixed_p2, 0)
  expect_equal(r$drought_duration_fixed_p2_mean, 0)
  expect_equal(r$drought_deficit_fixed_p2_mean, 0)
})

test_that("a low-flow block straddling Sep 30 -> Oct 1 is SPLIT across years", {
  # Codex MAJOR-1 on the Julia original: conservation alone would pass under a
  # uniform water-year shift, so attribution must be pinned explicitly.
  lo_start <- as.Date("2005-09-25"); lo_end <- as.Date("2005-10-07")
  f <- function(i, wy) {
    d <- as.Date(paste0(wy - 1, "-10-01")) + (i - 1L)
    if (d >= lo_start && d <= lo_end) 0 else 10
  }
  co <- annual_collector()
  calculate_drought_metrics(mk_dr(f), collector = co)
  d <- collector_drain(co)
  dur <- d[d$signature == "drought_duration_fixed_p30", ]
  v05 <- dur$value[dur$water_year == 2005]; v06 <- dur$value[dur$water_year == 2006]
  expect_gt(v05, 0); expect_gt(v06, 0)                    # SPLIT, not lumped
  expect_gte(v05 + v06, 13)
  others <- dur$value[!dur$water_year %in% c(2005, 2006)]
  expect_true(all(others == 0))                            # uniform-shift guard
})

test_that("schema contract on empty / missing-column / short records", {
  cp <- list(start_year=1900, end_year=2100, min_total_obs=20, min_segment_obs=10)
  empty <- data.frame(date = as.Date(character(0)), water_year = integer(0), Q = numeric(0))
  r <- calculate_drought_metrics(empty, changepoint = cp)
  for (m in c(DUR, DEF)) {
    expect_true(is.na(r[[paste0(m, "_mean")]]), info = m)
    expect_true(paste0(m, "_pettitt_cp_year") %in% names(r), info = m)
  }
  for (s in THR) expect_true(is.na(r[[s]]), info = s)

  r2 <- suppressWarnings(calculate_drought_metrics(data.frame(Q = c(1, 2)), changepoint = cp))
  for (m in c(DUR, DEF)) expect_true(is.na(r2[[paste0(m, "_mean")]]), info = m)

  # fewer than min_years_for_threshold -> NA thresholds, hence NA metrics
  r3 <- calculate_drought_metrics(mk_dr(function(i, wy) 1 + 0.5*sin(i/30), n_years = 5))
  for (s in THR) expect_true(is.na(r3[[s]]), info = s)
})

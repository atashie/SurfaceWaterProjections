# Annual-values collector (Phase 3). Contract mirrors julia/src/stats.jl:
# the collector receives the series EXACTLY as passed to generate_stats —
# before min_rows, the stats floor and the trend gates — and is read-only
# with respect to the statistics. NA values are KEPT; unkeyable years skipped.

CP <- list(start_year = 1900, end_year = 2100, min_total_obs = 20, min_segment_obs = 10)

adf <- function(n = 25, start = 1990) {
  data.frame(water_year = start:(start + n - 1),
             m1 = 10 + 0.5 * seq_len(n),
             m2 = 100 - 0.25 * seq_len(n))
}

test_that("collector is read-only with respect to the statistics", {
  df <- adf()
  a <- generate_stats(df, value_cols = c("m1", "m2"), changepoint = CP)
  co <- annual_collector()
  b <- generate_stats(df, value_cols = c("m1", "m2"), changepoint = CP, collector = co)
  expect_identical(a, b)
  expect_equal(nrow(collector_drain(co)), 50)   # 2 metrics x 25 years
})

test_that("collected values equal the input series", {
  df <- adf()
  co <- annual_collector()
  generate_stats(df, value_cols = c("m1", "m2"), collector = co)
  d <- collector_drain(co)
  expect_equal(d$value[d$signature == "m1"], df$m1)
  expect_equal(d$water_year[d$signature == "m1"], df$water_year)
  expect_true(is.integer(d$water_year))
  expect_equal(sum(duplicated(d[, c("signature", "water_year")])), 0)
})

test_that("collection happens BEFORE the min_rows gate and the stats floor", {
  short <- data.frame(water_year = c(2000, 2001), m = c(1, 2))
  co <- annual_collector()
  r <- generate_stats(short, value_cols = "m", min_rows = 3, collector = co)
  expect_true(is.na(r$m_mean))
  expect_equal(nrow(collector_drain(co)), 2)     # collected anyway

  df <- adf(19)
  co2 <- annual_collector()
  r2 <- generate_stats(df, value_cols = "m1", min_values_for_stats = 20, collector = co2)
  expect_true(is.na(r2$m1_mean))
  expect_equal(sum(collector_drain(co2)$signature == "m1"), 19)
})

test_that("NA values are kept; unkeyable years skipped", {
  df <- adf(10); df$m1[3:5] <- NA
  co <- annual_collector()
  generate_stats(df, value_cols = "m1", collector = co)
  d <- collector_drain(co)
  expect_equal(nrow(d), 10)
  expect_equal(sum(is.na(d$value)), 3)

  df2 <- data.frame(water_year = c(2000, NA, 2002), m = c(1, 2, 3))
  co2 <- annual_collector()
  generate_stats(df2, value_cols = "m", collector = co2)
  d2 <- collector_drain(co2)
  expect_equal(d2$water_year, c(2000L, 2002L))
  expect_equal(d2$value, c(1, 3))
})

test_that("empty collector drains to a 0-row frame with the right columns", {
  d <- collector_drain(annual_collector())
  expect_equal(nrow(d), 0)
  expect_equal(names(d), c("signature", "water_year", "value"))
})

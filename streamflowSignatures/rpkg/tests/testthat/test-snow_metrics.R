# Snow metrics (Phase 4). Expectations cross-checked against canonical Julia
# on 2026-08-26: 126 values over the three fixtures, 0 mismatches, worst
# relative difference 3.9e-16.

mk_snow <- function(f, n_years = 25, ppt = 2) {
  rows <- list()
  for (wy in 1996:(1995 + n_years)) {
    d0 <- as.Date(paste0(wy - 1, "-10-01")); d1 <- as.Date(paste0(wy, "-09-30"))
    nd <- as.integer(d1 - d0) + 1L
    rows[[length(rows) + 1L]] <- data.frame(
      water_year = wy, dowy = 1:nd,
      SWE = vapply(1:nd, function(i) f(i, wy), numeric(1)), PPT = ppt)
  }
  do.call(rbind, rows)
}
triangle <- function(i, wy) if (i <= 90) 0 else if (i <= 180) (i-90)*1.5 else
  if (i <= 260) max(0, 135 - (i-180)*1.7) else 0
thaw <- function(i, wy) if (i <= 60) 0 else if (i <= 110) (i-60)*2 else
  if (i <= 130) 0 else if (i <= 200) (i-130)*2.5 else
  if (i <= 270) max(0, 175-(i-200)*2.6) else 0
vanishing <- function(i, wy) if (wy <= 2015) triangle(i, wy) else if (i > 250) 0 else 3

test_that("snow_spells finds maximal runs", {
  sp <- streamflowsignatures:::snow_spells(c(FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE))
  expect_equal(lapply(sp, as.integer), list(2:3, 5L, 8L))
  expect_equal(length(streamflowsignatures:::snow_spells(rep(FALSE, 5))), 0)
})

test_that("triangle gage: every metric is hand-derivable", {
  r <- calculate_snow_metrics(mk_snow(triangle), trend_completeness = 0.6,
                              decade_completeness = 0.8, min_values_for_stats = 20)
  expect_equal(r$swe_max_mean, 135)
  expect_equal(r$swe_max_dowy_mean, 180)
  expect_equal(r$snow_on_dowy_mean, 97)      # 1.5k >= 10 -> k = 7
  expect_equal(r$snow_off_dowy_mean, 254)    # 135 - 1.7k >= 10 -> k <= 73.5
  expect_equal(r$snow_cover_days_mean, 254 - 97)
  expect_equal(r$melt_season_days_mean, 254 - 180)
  expect_equal(r$melt_rate_mean, 135 / 74)
  expect_equal(r$ssm_mean, 1)                # single spell >= 60 days
  expect_equal(r$melt_before_peak_mean, 0)   # monotone rise to the peak
})

test_that("thaw gage: anchor spell is the PEAK spell; SSM pools all spells", {
  r <- calculate_snow_metrics(mk_snow(thaw), trend_completeness = 0.6,
                              decade_completeness = 0.8, min_values_for_stats = 20)
  expect_equal(r$swe_max_mean, 175)
  expect_equal(r$swe_max_dowy_mean, 200)
  expect_gt(r$snow_on_dowy_mean, 130)                 # second spell, not the first
  expect_equal(r$melt_before_peak_mean, 100)          # first spell melts entirely
  expect_equal(r$melt_before_peak_to_max_swe_mean, 100/175)
  expect_lt(r$ssm_mean, 1); expect_gt(r$ssm_mean, 0)
})

test_that("snow-free years: magnitudes are ZERO, timing metrics NA", {
  r <- calculate_snow_metrics(mk_snow(function(i, wy) 5),   # always below threshold
                              trend_completeness = 0.6, decade_completeness = 0.8)
  expect_equal(r$swe_max_mean, 0)
  expect_equal(r$snow_cover_days_mean, 0)
  expect_equal(r$swe_apr1_mean, 0)
  for (m in c("swe_max_dowy","snow_on_dowy","snow_off_dowy","melt_season_days",
              "melt_rate","ssm","melt_com_dowy")) {
    expect_true(is.na(r[[paste0(m, "_mean")]]), info = m)
  }
})

test_that("boundary censoring when the anchor spell touches a year edge", {
  r <- calculate_snow_metrics(mk_snow(function(i, wy) 100),
                              trend_completeness = 0.6, decade_completeness = 0.8)
  expect_true(is.na(r$snow_on_dowy_mean))
  expect_true(is.na(r$snow_off_dowy_mean))
  expect_true(is.na(r$melt_season_days_mean))
})

test_that("incomplete / NA years are skipped entirely", {
  df <- mk_snow(triangle)
  df <- df[!(df$water_year == 2000 & df$dowy == 50), ]     # punch a hole
  co <- annual_collector()
  calculate_snow_metrics(df, collector = co)
  d <- collector_drain(co)
  expect_true(is.na(d$value[d$signature == "swe_max" & d$water_year == 2000]))
})

test_that("record-anchored decade gate suppresses trends, keeps mean/median", {
  r <- calculate_snow_metrics(mk_snow(vanishing), trend_completeness = 0.6,
                              decade_completeness = 0.8, min_values_for_stats = 20)
  for (m in c("swe_max_dowy", "ssm", "melt_com_dowy")) {
    expect_true(is.na(r[[paste0(m, "_senn_slp")]]), info = m)
    expect_false(is.na(r[[paste0(m, "_mean")]]), info = m)
  }
  expect_false(is.na(r$swe_max_mean))    # magnitude metrics are EXEMPT
})

test_that("no SWE column still emits the full 14-metric key set", {
  r <- calculate_snow_metrics(data.frame(water_year = integer(0), dowy = integer(0)),
                              changepoint = list(start_year=1900, end_year=2100,
                                                 min_total_obs=20, min_segment_obs=10))
  for (m in SNOW_METRICS) {
    expect_true(paste0(m, "_mean") %in% names(r), info = m)
    expect_true(paste0(m, "_pettitt_cp_year") %in% names(r), info = m)
  }
})

test_that("a ZERO-ROW snow_data frame returns the full key set, not an error", {
  # Regression for 2026-08-26: the orchestrator passes an explicit snow_data
  # filtered to valid_swe_years, which is legitimately EMPTY for a gage with SWE
  # columns but no SWE-valid year. snow.R then built a 0-row `annual` and
  # assigned a length-1 NA into it -> "replacement has 1 row, data has 0", which
  # the orchestrator's try/catch swallowed into a warning, silently costing that
  # gage all 224 snow columns. Found on 4 Florida gages by the run-log gate; no
  # unit test covered a 0-row frame with the SWE column present.
  d <- data.table::data.table(
    gage_id = "T", date = seq(as.Date("1992-10-01"), as.Date("1995-09-30"), by = "day"))
  d[, Q := 1]; d[, SWE := 0]
  d <- add_water_year_columns(d)

  r0 <- calculate_snow_metrics(d[0], valid_climate_years = integer(0))
  rN <- calculate_snow_metrics(d, valid_climate_years = unique(d$water_year))

  expect_identical(sort(names(r0)), sort(names(rN)))   # identical schema
  expect_true(all(is.na(unlist(r0))))                  # and all-NA values
  expect_gt(length(r0), 0)

  # ...and through the orchestrator, the way the runner calls it
  s <- calculate_all_signatures(d, has_climate = FALSE, snow_data = d[0],
                                snow_climate_years = integer(0), gage_id = "T")
  expect_true(any(grepl("^swe_max", names(s))))
})

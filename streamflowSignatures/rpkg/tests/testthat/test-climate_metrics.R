# Climate metrics (guidelines 3.9; decisions 2026-10-12). Port of
# julia/test/test_climate_metrics.jl: every expectation is derived by hand from the
# definitions in docs/plans/2026-10-12-climate-metrics-plan.md. Annual values are read
# back through the annual collector (collected BEFORE any gating).

climate_year_frame <- function(wy, P = NULL, TN = NULL, TX = NULL, Q = NULL) {
  d0 <- as.Date(paste0(wy - 1, "-10-01")); d1 <- as.Date(paste0(wy, "-09-30"))
  dates <- seq.Date(d0, d1, by = "day")
  n <- length(dates)
  df <- data.frame(gage_id = "TCLIM", date = dates, Q = if (is.null(Q)) rep(1, n) else Q)
  df <- add_water_year_columns(df)
  if (!is.null(P)) df$PPT <- as.numeric(P)
  if (!is.null(TN)) df$TMIN <- as.numeric(TN)
  if (!is.null(TX)) df$TMAX <- as.numeric(TX)
  as.data.frame(df)
}

annual_lookup <- function(c, sig, wy) {
  d <- collector_drain(c)
  v <- d$value[d$signature == sig & d$water_year == wy]
  if (length(v) == 0) NULL else v[1]
}

collect_precip <- function(df) {
  c <- annual_collector()
  r <- calculate_precipitation_metrics(df, collector = c)
  list(r = r, c = c)
}
collect_temp <- function(df, ...) {
  c <- annual_collector()
  r <- calculate_temperature_metrics(df, collector = c, ...)
  list(r = r, c = c)
}

test_that("registry: 20 precipitation + 19 temperature bases with the code's season suffixes", {
  expect_equal(length(precip_metrics()), 20)
  expect_equal(length(temp_metrics()), 19)
  expect_true(all(c("prcp_5d_max", "prcp_50_dowy", "prcp_ndays_abv_10mm") %in% precip_metrics()))
  for (s in c("win", "spr", "sum", "fal")) {
    expect_true(paste0("prcp_", s) %in% precip_metrics())
    expect_true(paste0("gdd_cum_", s) %in% temp_metrics())
  }
  expect_false("prcp_sumr" %in% precip_metrics())
})

test_that("single-storm year: timing, concentration and spells by hand", {
  wy <- 2001
  P <- rep(0, 365); P[100] <- 100
  o <- collect_precip(climate_year_frame(wy, P = P))
  v <- function(s) annual_lookup(o$c, s, wy)
  expect_equal(v("prcp_annual"), 100)
  expect_equal(v("prcp_1d_max"), 100); expect_equal(v("prcp_5d_max"), 100)
  expect_equal(v("prcp_max_dowy"), 100); expect_equal(v("prcp_50_dowy"), 100)
  expect_equal(v("prcp_cm_dowy"), 100)
  expect_equal(v("prcp_ndays_abv_10mm"), 1); expect_equal(v("prcp_ndays_abv_50mm"), 1)
  expect_equal(v("prcp_annual_sd"), sd(P), tolerance = 1e-12)
  expect_equal(v("prcp_gini"), 364 / 365, tolerance = 1e-12)
  expect_equal(v("prcp_seasonality_index"), 11 / 6, tolerance = 1e-12)
  expect_equal(v("prcp_intermit_annual"), 265)
  expect_equal(v("prcp_intermit_annual_dowy_start"), 101)
  expect_equal(v("prcp_intermit_annual_dowy_end"), 365)
  expect_equal(v("prcp_win"), 100); expect_equal(v("prcp_fal"), 0)
  expect_equal(length(o$r), 20 * 8)
})

test_that("uniform year, thresholds, ties and the fall season convention", {
  wy <- 2001
  o <- collect_precip(climate_year_frame(wy, P = rep(2, 365)))
  expect_equal(annual_lookup(o$c, "prcp_50_dowy", wy), 183)
  expect_equal(annual_lookup(o$c, "prcp_cm_dowy", wy), 183, tolerance = 1e-12)
  expect_equal(annual_lookup(o$c, "prcp_intermit_annual", wy), 0)
  expect_true(is.na(annual_lookup(o$c, "prcp_intermit_annual_dowy_start", wy)))
  expect_equal(annual_lookup(o$c, "prcp_gini", wy), 0, tolerance = 1e-12)

  P <- rep(5, 365); P[10] <- 0.1; P[11] <- 0.1000001; P[12] <- 0; P[13] <- 0.1
  P[50] <- 10; P[51] <- 10.0001; P[52] <- 25; P[53] <- 50; P[54] <- 50.5
  o2 <- collect_precip(climate_year_frame(wy, P = P))
  expect_equal(annual_lookup(o2$c, "prcp_ndays_abv_10mm", wy), 4)
  expect_equal(annual_lookup(o2$c, "prcp_ndays_abv_25mm", wy), 2)
  expect_equal(annual_lookup(o2$c, "prcp_ndays_abv_50mm", wy), 1)
  expect_equal(annual_lookup(o2$c, "prcp_intermit_annual", wy), 2)
  expect_equal(annual_lookup(o2$c, "prcp_intermit_annual_dowy_start", wy), 12)

  base <- climate_year_frame(wy)
  o3 <- collect_precip(climate_year_frame(wy, P = base$month))
  expect_equal(annual_lookup(o3$c, "prcp_win", wy), 12 * 31 + 31 + 2 * 28)
  expect_equal(annual_lookup(o3$c, "prcp_fal", wy), 9 * 30 + 10 * 31 + 11 * 30)

  P4 <- rep(0, 365); P4[c(1, 101, 201, 301)] <- 1
  o4 <- collect_precip(climate_year_frame(wy, P = P4))
  expect_equal(annual_lookup(o4$c, "prcp_intermit_annual", wy), 99)
  expect_equal(annual_lookup(o4$c, "prcp_intermit_annual_dowy_start", wy), 2)
})

test_that("crossing-tolerance boundary (decision 2026-10-12: 1e-12 of the total)", {
  # two storms of 1 - d and 1 + d mm (total 2 mm): d = 1e-12 lies inside the allowance
  # (the first storm crosses, day 7); d = 4e-12 lies outside (the second, day 301)
  expect_equal(streamflowsignatures:::.climate_crossing_rel_tol, 1e-12)
  wy <- 2001
  for (case in list(c(1e-12, 7), c(4e-12, 301))) {
    P <- rep(0, 365); P[7] <- 1 - case[1]; P[301] <- 1 + case[1]
    expect_equal(annual_lookup(collect_precip(climate_year_frame(wy, P = P))$c, "prcp_50_dowy", wy), case[2])
  }
})

test_that("preprocessing fidelity (review 2026-10-12): observations are real numbers; SWE policy from the config", {
  d0 <- as.Date("2000-10-01"); dates <- seq.Date(d0, as.Date("2001-09-30"), by = "day")
  base <- as.data.frame(add_water_year_columns(data.frame(gage_id = "T", date = dates, Q = 1, PPT = 2, SWE = 1,
                                                          TMIN = 1, TMAX = 10, stringsAsFactors = FALSE)))
  # strings are missing observations (Julia's _real_vector): a character column is entirely
  # missing; in a list column one string fills, four reject
  df <- base; df$PPT <- as.character(df$PPT)
  expect_equal(length(preprocess_daily_data(df)$valid_climate_years), 0)
  df <- base; df$PPT <- as.list(df$PPT); df$PPT[[101]] <- "3"
  expect_equal(preprocess_daily_data(df)$valid_climate_years, 2001)
  df <- base; df$PPT <- as.list(df$PPT); df$PPT[101:104] <- list("3", "3", "3", "3")
  expect_equal(length(preprocess_daily_data(df)$valid_climate_years), 0)
  df <- base; df$TMIN <- as.character(df$TMIN); df$TMAX <- as.character(df$TMAX)
  expect_equal(length(preprocess_daily_data(df)$valid_temp_years), 0)
  df <- base; df$TMIN <- complex(real = df$TMIN, imaginary = 0); df$TMAX <- complex(real = df$TMAX, imaginary = 0)
  expect_equal(length(preprocess_daily_data(df)$valid_temp_years), 0)
  # the SWE policy honours the supplied config (limits typed, flag a real logical)
  df <- base; df$SWE[101] <- NA
  expect_equal(preprocess_daily_data(df)$valid_swe_years, 2001)
  expect_equal(length(preprocess_daily_data(df, config = list(snow_na_policy = list(
    max_raw_na_per_year_swe = 0, max_interpolation_gap_swe = 0)))$valid_swe_years), 0)
  expect_error(preprocess_daily_data(df, config = list(snow_na_policy = list(max_raw_na_per_year_swe = 3.5))))
  expect_error(preprocess_daily_data(df, config = list(snow_na_policy = list(max_raw_na_per_year_swe = TRUE))))
  expect_error(preprocess_daily_data(df, config = list(snow_na_policy = list(reject_negative_swe = 1))))
  df <- base; df$SWE[101] <- -1
  expect_equal(length(preprocess_daily_data(df)$valid_swe_years), 0)
  expect_equal(preprocess_daily_data(df, config = list(snow_na_policy = list(reject_negative_swe = FALSE)))$valid_swe_years, 2001)
})

test_that("list-typed observation columns come out numeric and feed the existing families (review round 5)", {
  d <- as.data.frame(add_water_year_columns(data.frame(
    gage_id = "T", date = seq(as.Date("2000-10-01"), as.Date("2003-09-30"), by = "day"), Q = 1, PPT = 1,
    stringsAsFactors = FALSE)))
  d$PPT <- as.list(d$PPT); d$PPT[[100]] <- "3"
  d$Q <- as.list(d$Q)
  p <- preprocess_daily_data(d)
  expect_equal(p$valid_years, c(2001, 2002, 2003))
  expect_equal(p$valid_climate_years, c(2001, 2002, 2003))
  expect_true(is.double(p$data$PPT) && is.double(p$data$Q))
  expect_equal(p$data$PPT[100], 1)
  r <- calculate_all_signatures(p$data, TRUE, climate_data = p$data, area_normalized = TRUE)
  expect_equal(r$Qann_mean, 365)
  expect_equal(r$annual_runoff_ratio_mean, 1)
  expect_equal(r$prcp_annual_mean, 365)
})

test_that("round-2 regressions: crossing tolerance, non-finite inputs are missing", {
  wy <- 2001
  P <- rep(0, 365); P[1:7] <- 0.1; P[301:307] <- 0.1
  expect_equal(annual_lookup(collect_precip(climate_year_frame(wy, P = P))$c, "prcp_50_dowy", wy), 7)
  expect_equal(annual_lookup(collect_precip(climate_year_frame(2000, P = rep(0.1, 366)))$c, "prcp_50_dowy", 2000), 183)
  P3 <- rep(0, 365); P3[1] <- 49; P3[2] <- 51
  expect_equal(annual_lookup(collect_precip(climate_year_frame(wy, P = P3))$c, "prcp_50_dowy", wy), 2)
  Pi <- rep(2, 365); Pi[100] <- Inf
  ppi <- preprocess_daily_data(climate_year_frame(wy, P = Pi))
  expect_equal(ppi$valid_climate_years, wy)
  expect_equal(as.numeric(ppi$data$PPT[100]), 2)
  Pn <- rep(2, 365); Pn[100] <- -Inf
  expect_equal(preprocess_daily_data(climate_year_frame(wy, P = Pn))$valid_climate_years, wy)
  Pm <- rep(2, 365); Pm[100:103] <- Inf
  expect_equal(length(preprocess_daily_data(climate_year_frame(wy, P = Pm))$valid_climate_years), 0)
  Qi <- rep(1, 365); Qi[1:31] <- Inf
  expect_equal(length(preprocess_daily_data(climate_year_frame(wy, P = rep(1, 365), Q = Qi))$valid_years), 0)
  f <- streamflowsignatures:::.fill_internal_gaps(c(0, NA, Inf), 3)
  expect_true(is.na(f[2]) && is.infinite(f[3]))
  expect_error(preprocess_daily_data(climate_year_frame(wy, P = rep(1, 365)),
                                     config = list(temperature_na_policy = list(max_interpolation_gap_temp = 3.5))))
  expect_equal(streamflowsignatures:::pkg_env$climate_cum_fraction, 0.5)
})

test_that("review regressions: equal blocks, truncated leap year, month column ignored", {
  wy <- 2001
  P <- rep(0, 365); P[1:50] <- 0.1; P[301:350] <- 0.1
  o <- collect_precip(climate_year_frame(wy, P = P))
  expect_equal(annual_lookup(o$c, "prcp_50_dowy", wy), 50)
  trunc <- climate_year_frame(2000, P = rep(2, 366), TN = rep(-1, 366), TX = rep(5, 366))[1:365, ]
  expect_true(is.na(annual_lookup(collect_precip(trunc)$c, "prcp_annual", 2000)))
  expect_true(is.na(annual_lookup(collect_temp(trunc)$c, "gdd_cum_annual", 2000)))
  df_m <- climate_year_frame(wy, P = rep(2, 365)); df_m$month[1] <- NA
  o4 <- collect_precip(df_m)
  parts <- sum(sapply(c("win", "spr", "sum", "fal"), function(s) annual_lookup(o4$c, paste0("prcp_", s), wy)))
  expect_equal(parts, 730)
})

test_that("temperature: extremes, seasons, frost window, GDD clamp", {
  wy <- 2001
  base <- climate_year_frame(wy)
  TX <- as.numeric(base$month); TN <- TX - 30
  o <- collect_temp(climate_year_frame(wy, TN = TN, TX = TX))
  v <- function(s) annual_lookup(o$c, s, wy)
  expect_equal(v("temp_max_annual"), 12); expect_equal(v("temp_max_annual_dowy"), 62)
  expect_equal(v("temp_min_annual"), -29); expect_equal(v("temp_min_annual_dowy"), 93)
  expect_equal(v("temp_max_fal"), 11); expect_equal(v("temp_min_sum"), -24)
  expect_equal(v("frost_first_dowy"), 1); expect_equal(v("frost_last_dowy"), 273)
  expect_equal(v("gdd_cum_annual"), 0)

  TN2 <- rep(5, 365); TN2[20] <- -1; TN2[250] <- -0.5; TN2[300] <- -2
  o2 <- collect_temp(climate_year_frame(wy, TN = TN2, TX = rep(15, 365)))
  expect_equal(annual_lookup(o2$c, "frost_first_dowy", wy), 20)
  expect_equal(annual_lookup(o2$c, "frost_last_dowy", wy), 250)
  TN3 <- rep(5, 365); TN3[273] <- -1; TN3[274] <- -1
  o3 <- collect_temp(climate_year_frame(wy, TN = TN3, TX = rep(15, 365)))
  expect_equal(annual_lookup(o3$c, "frost_last_dowy", wy), 273)
  o4 <- collect_temp(climate_year_frame(wy, TN = rep(0, 365), TX = rep(15, 365)))
  expect_true(is.na(annual_lookup(o4$c, "frost_first_dowy", wy)))

  o5 <- collect_temp(climate_year_frame(wy, TN = rep(10, 365), TX = rep(30, 365)))
  expect_equal(annual_lookup(o5$c, "gdd_cum_annual", wy), 20 * 365)
  expect_equal(annual_lookup(o5$c, "gdd_cum_win", wy), 20 * (31 + 31 + 28))
  o6 <- collect_temp(climate_year_frame(wy, TN = rep(-20, 365), TX = rep(10, 365)))
  expect_equal(annual_lookup(o6$c, "gdd_cum_annual", wy), 0)
})

test_that("record-anchored decade gate on the frost metrics", {
  frames <- list()
  for (k in 1:20) {
    wy <- 2000 + k
    n <- as.integer(as.Date(paste0(wy, "-09-30")) - as.Date(paste0(wy - 1, "-10-01"))) + 1L
    TN <- rep(5, n); if (k <= 10) TN[30 + k] <- -1
    frames[[k]] <- climate_year_frame(wy, TN = TN, TX = rep(15, n))
  }
  df <- do.call(rbind, frames)
  r <- calculate_temperature_metrics(df, trend_completeness = 0.6, decade_completeness = 0.8)
  expect_true(is.na(r$frost_first_dowy_senn_slp))
  expect_equal(r$frost_first_dowy_mean, mean(31:40))
  expect_false(is.na(r$temp_max_annual_senn_slp))
  r_open <- calculate_temperature_metrics(df)
  expect_equal(r_open$frost_first_dowy_senn_slp, 1, tolerance = 1e-12)
})

test_that("preprocessor: temperature validity mirrors the SWE policy without a sign check", {
  frames <- list()
  for (wy in 2001:2004) {
    n <- as.integer(as.Date(paste0(wy, "-09-30")) - as.Date(paste0(wy - 1, "-10-01"))) + 1L
    P <- rep(2, n); TN <- rep(-15, n); TX <- rep(5, n)
    if (wy == 2002) TN[100:101] <- NA
    if (wy == 2003) TX[50:85] <- NA
    if (wy == 2004) TN[200:203] <- NA
    frames[[length(frames) + 1]] <- climate_year_frame(wy, P = P, TN = TN, TX = TX)
  }
  df <- do.call(rbind, frames)
  pp <- preprocess_daily_data(df)
  expect_equal(pp$valid_years, 2001:2004)
  expect_equal(pp$valid_climate_years, 2001:2004)
  expect_equal(pp$valid_temp_years, 2001:2002)
  y2002 <- pp$data[pp$data$water_year == 2002, ]
  expect_false(any(is.na(y2002$TMIN)))
  pp0 <- preprocess_daily_data(df[, setdiff(names(df), c("TMIN", "TMAX"))])
  expect_equal(pp0$valid_temp_years, integer(0))
})

test_that("orchestrator: explicit PPT-qualified and temperature frames only", {
  frames <- list()
  for (wy in 1991:2015) {
    base <- climate_year_frame(wy)
    n <- nrow(base); d <- seq_len(n)
    P <- ifelse(d %% 7 == 0, 0, base$month + 0.01 * d)
    TX <- 20 + 15 * sin(2 * pi * (d - 100) / 365); TN <- TX - 12
    frames[[length(frames) + 1]] <- climate_year_frame(wy, P = P, TN = TN, TX = TX)
  }
  df <- do.call(rbind, frames)
  pp <- preprocess_daily_data(df)
  data <- as.data.frame(pp$data)
  cdata <- data[data$water_year %in% pp$valid_climate_years, ]
  tdata <- data[data$water_year %in% pp$valid_temp_years, ]
  r_bare <- calculate_all_signatures(data, has_climate = TRUE)
  expect_false("prcp_annual_mean" %in% names(r_bare))
  expect_false("temp_max_annual_mean" %in% names(r_bare))
  r1 <- calculate_all_signatures(data, has_climate = TRUE, climate_data = cdata)
  expect_true("prcp_annual_mean" %in% names(r1))
  expect_false("temp_max_annual_mean" %in% names(r1))
  r2 <- calculate_all_signatures(data, has_climate = TRUE, climate_data = cdata, temp_data = tdata)
  expect_true("temp_max_annual_mean" %in% names(r2))
  expect_false(is.na(r2$temp_max_annual_mean))
  r3 <- calculate_all_signatures(data, has_climate = TRUE, climate_data = cdata, temp_data = tdata[0, ])
  expect_true(is.na(r3$temp_max_annual_mean))
  r4 <- calculate_all_signatures(data, has_climate = TRUE, climate_data = cdata, area_normalized = FALSE)
  expect_true("prcp_annual_mean" %in% names(r4))
  expect_false("annual_runoff_ratio_mean" %in% names(r4))
})

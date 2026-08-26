# Tests for the fixed-b=1 recession alpha convention (Phase 2 of the port
# campaign; mirrors python/tests/test_recession_alpha_b1.py and
# julia/test/test_recession_alpha_b1.jl).
#
# Convention (July 2026, canonical): ALL alpha outputs (log_a_pointcloud,
# log_a_events, log_a_seasonality_*) assume a linear reservoir —
# log(a) = median(log(-dQ/dt) - log(Q)), no regression — while b_pointcloud /
# b_events / concavity keep their FREE power-law fits. The quadratic-reservoir
# gage is the discriminating case: the free fit recovers b = 2 exactly while
# the b=1 alpha differs strongly from the free-fit intercept log(k).

daily_frame <- function(n_years, q_update, peak = 10) {
  dates <- seq(as.Date("1990-10-01"), by = "day", length.out = 365 * n_years)
  Q <- numeric(length(dates))
  q <- peak
  dom <- as.integer(format(dates, "%d"))
  mon <- as.integer(format(dates, "%m"))
  for (i in seq_along(dates)) {
    q <- if (dom[i] == 1L) peak else q_update(q, mon[i])
    Q[i] <- max(q, 1e-6)
  }
  add_water_year_columns(data.frame(gage_id = "TESTB1", date = dates, Q = Q))
}

linear_gage <- function(n_years = 25, alpha = 0.95) {
  daily_frame(n_years, function(q, m) q * alpha)
}
quadratic_gage <- function(n_years = 25, k = 0.02) {
  daily_frame(n_years, function(q, m) q - k * q^2)
}
seasonal_gage <- function(n_years = 25) {
  daily_frame(n_years, function(q, m) q * if (m %in% c(10, 11, 12, 1, 2, 3)) 0.90 else 0.98)
}

test_that("event_log_a_b1 helper: exact on a pure exponential", {
  Q_exp <- 10 * 0.9^(0:19)
  expect_equal(streamflowsignatures:::event_log_a_b1(Q_exp), log(0.1), tolerance = 1e-12)
  expect_true(is.na(streamflowsignatures:::event_log_a_b1(c(5, 4))))        # too short
  expect_true(is.na(streamflowsignatures:::event_log_a_b1(c(1, 2, 3, 4))))  # increasing
})

test_that("linear reservoir: b=1 convention coincides with the free fit", {
  r <- suppressWarnings(analyze_recession_parameters(linear_gage(alpha = 0.95)))
  expect_equal(r[["log_a_pointcloud_mean"]], log(0.05), tolerance = 1e-9)
  expect_equal(r[["log_a_events_mean"]], log(0.05), tolerance = 1e-9)
  expect_equal(r[["b_pointcloud_mean"]], 1, tolerance = 1e-6)
  expect_equal(r[["b_events_mean"]], 1, tolerance = 1e-6)
  expect_lt(abs(r[["concavity_mean"]]), 1e-9)
  expect_equal(r[["recession_alpha_point_cloud_linear_reservoir"]], 0.95, tolerance = 1e-9)
  # log_a = log(1 - alpha) consistency with the discrete alpha
  expect_equal(r[["log_a_pointcloud_mean"]],
               log(1 - r[["recession_alpha_point_cloud_linear_reservoir"]]), tolerance = 1e-9)
  # event counts unchanged by the convention (monthly spikes -> 12/yr)
  expect_equal(r[["n_recession_events_mean"]], 12, tolerance = 1e-9)
})

test_that("quadratic reservoir: b stays free at 2, alpha is decoupled", {
  k <- 0.02
  df <- quadratic_gage(k = k)
  r <- suppressWarnings(analyze_recession_parameters(df))

  # FREE fits must still recover b = 2 (-dQ = k*Q^2 exactly by construction)
  expect_equal(r[["b_pointcloud_mean"]], 2, tolerance = 1e-6)
  expect_equal(r[["b_events_mean"]], 2, tolerance = 1e-6)

  # The OLD behavior (fitted-b intercept) would give ~log(k); the b=1 value is
  # median(log(k*Q_j)) with Q_j in (Q_min, peak) — strictly different.
  q_min <- min(df$Q)
  for (key in c("log_a_pointcloud_mean", "log_a_events_mean")) {
    expect_gt(r[[key]], log(k * q_min))
    expect_lt(r[[key]], log(k * 10))
    expect_gt(abs(r[[key]] - log(k)), 0.3)
  }
  expect_false(is.na(r[["log_a_seasonality_amplitude_all"]]))
})

test_that("seasonality consumes the per-event b=1 values at the canonical mid-day", {
  # Seasonally alternating alpha => strong seasonal log_a signal. The pipeline's
  # seasonality must reproduce a sinusoid fit to independently recomputed
  # (per-event b=1 log_a, mid-event dowy) pairs — a misalignment would change it.
  df <- seasonal_gage()
  r <- suppressWarnings(analyze_recession_parameters(df))

  la <- numeric(0); dw <- numeric(0)
  for (yr in unique(df$water_year)) {
    yd <- df[df$water_year == yr, ]
    yd <- yd[order(yd$dowy), ]
    for (ev in streamflowsignatures:::identify_recession_events(yd$Q)) {
      Qe <- yd$Q[ev$indices]
      p <- streamflowsignatures:::fit_recession_event(Qe, remove_first_day = TRUE)
      if (is.na(p$log_a) || is.na(p$b)) next
      v <- streamflowsignatures:::event_log_a_b1(Qe)
      if (is.na(v)) next
      la <- c(la, v)
      # canonical mid-event index (1-based ceiling(n/2) == Julia's floor midpoint)
      dw <- c(dw, yd$dowy[ev$indices[ceiling(length(ev$indices) / 2)]])
    }
  }
  expect_gte(length(la), 250)   # ~12 events/yr x 25 yr
  fit <- streamflowsignatures:::fit_sinusoidal_model(dw, la)

  expect_equal(r[["log_a_seasonality_amplitude_all"]], fit$amplitude, tolerance = 1e-9)
  expect_equal(r[["log_a_seasonality_minimum_all"]], fit$minimum_doy, tolerance = 1e-9)

  # Signal must be non-trivial: log_a alternates between log(0.10) and log(0.02)
  expect_gt(r[["log_a_seasonality_amplitude_all"]], 0.4)
  # Minimum (most negative log_a = slow-recession season) in Apr-Sep
  expect_gt(r[["log_a_seasonality_minimum_all"]], 200)
  expect_lt(r[["log_a_seasonality_minimum_all"]], 350)
})

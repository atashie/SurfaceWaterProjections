# Tests for the area_normalized gate on Q-to-PPT signatures (July 2026).
#
# Gages with no drainage area (HYDAT gap — canals, dam outflows, channel splits)
# carry Q in raw m3/s instead of mm/day. Q-to-PPT signatures (runoff ratios,
# elasticity, Q-P seasonality, storage) mix Q and PPT units and are therefore
# skipped when area_normalized = FALSE. Q-only signatures must be unaffected.

make_gate_gage <- function() {
  set.seed(7)
  dates <- seq(as.Date("1990-10-01"), by = "day", length.out = 365 * 25)
  day_of_year <- as.numeric(format(dates, "%j"))
  Q <- pmax(0.01, 2.0 + 1.5 * sin(2 * pi * (day_of_year - 100) / 365) +
              runif(length(dates)) * 0.5)
  df <- data.frame(date = dates, Q = Q)
  df <- add_water_year_columns(df)
  df$PPT <- pmax(0, df$Q * 2.5 + runif(nrow(df)))
  df
}

clim_prefixes <- c(
  "annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
  "summer_runoff_ratio", "fall_runoff_ratio", "runoff_ratio_high_count",
  "elasticity_static", "elasticity_rolling", "elasticity_annual",
  "elasticity_years_total", "elasticity_years_low_ppt",
  "qp_slope_sd", "qp_bimodality", "avg_storage"
)
is_clim_key <- function(k) any(startsWith(k, clim_prefixes))

qonly_sentinels <- c("Qann_mean", "Q50_mean", "BFI_Eckhardt_mean",
                     "flashinessRB_mean", "D50_day_mean", "TQmean_mean",
                     "FDCall_mean")

test_that("area_normalized = FALSE skips ALL Q-to-PPT signatures", {
  df <- make_gate_gage()

  r_on <- suppressWarnings(suppressMessages(
    calculate_all_signatures(df, has_climate = TRUE)))
  r_off <- suppressWarnings(suppressMessages(
    calculate_all_signatures(df, has_climate = TRUE, area_normalized = FALSE)))

  expect_gt(length(Filter(is_clim_key, names(r_on))), 0)
  expect_length(Filter(is_clim_key, names(r_off)), 0)

  # Q-only signatures unaffected
  for (k in qonly_sentinels) {
    expect_true(k %in% names(r_off), info = k)
  }
})

test_that("gated result is identical to a gage without climate data", {
  df <- make_gate_gage()

  r_gated <- suppressWarnings(suppressMessages(
    calculate_all_signatures(df, has_climate = TRUE, area_normalized = FALSE)))
  # add_water_year_columns returns a data.table: an expression j needs
  # with = FALSE to select columns (a bare expression yields the character
  # vector itself under current data.table — broke this test on 1.16+)
  df_noclimate <- df[, setdiff(names(df), "PPT"), with = FALSE]
  r_noclimate <- suppressWarnings(suppressMessages(
    calculate_all_signatures(df_noclimate, has_climate = FALSE)))

  expect_identical(sort(names(r_gated)), sort(names(r_noclimate)))
  for (k in names(r_gated)) {
    a <- r_gated[[k]]
    b <- r_noclimate[[k]]
    same <- isTRUE(all.equal(a, b)) ||
      (is.numeric(a) && is.numeric(b) && is.na(a) && is.na(b))
    expect_true(same, info = k)
  }
})

test_that("gate is inert when area_normalized = TRUE", {
  df <- make_gate_gage()

  r_default <- suppressWarnings(suppressMessages(
    calculate_all_signatures(df, has_climate = TRUE)))
  r_explicit <- suppressWarnings(suppressMessages(
    calculate_all_signatures(df, has_climate = TRUE, area_normalized = TRUE)))

  expect_identical(sort(names(r_default)), sort(names(r_explicit)))
})

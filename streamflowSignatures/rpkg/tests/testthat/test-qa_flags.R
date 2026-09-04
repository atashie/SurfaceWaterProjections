# flagged_for_high_na denominator — regression for the 2026-09-04 finding.
# Mirrors julia/test/test_qa_high_na.jl and python/tests/test_qa_flags.py.

test_that("manifest is loaded from config", {
  expect_length(streamflowsignatures:::pkg_env$qa_high_na_suffixes, 16)
  expect_length(streamflowsignatures:::pkg_env$qa_high_na_scalars, 21)
})

test_that("denominator selects signature columns by name only", {
  hdr <- c("gage_id", "latitude", "NDAMS_2009", "area_normalized", "gage_type",
           "Qann_mean", "Qann_median", "Q5_pettitt_cp_year", "elasticity_static",
           "swe_max_mean", "flagged_for_qann_range", "num_water_years")
  expect_equal(high_na_denominator_columns(hdr),
               c("Qann_mean", "Qann_median", "Q5_pettitt_cp_year", "elasticity_static", "swe_max_mean"))
  bases <- paste0("base", 1:100)
  full <- c(paste0(rep(bases, each = 16), streamflowsignatures:::pkg_env$qa_high_na_suffixes),
            streamflowsignatures:::pkg_env$qa_high_na_scalars, hdr)
  expect_length(high_na_denominator_columns(full), 1621 + 5)
})

.frame <- function() {
  data.frame(
    gage_id = c("g1", "g2", "g3", "g4"),
    Qann_mean = c(1, 1, 1, NA),
    Qann_median = c(1, NaN, 1, NA),                 # NaN counts as NA
    Q5_pettitt_cp_year = c(2000, NA, NA, NA),
    elasticity_static = c(1.2, NA, NA, NA),         # registered scalar
    swe_max_mean = c(NA, NA, 3, NA),                # NA-filled family counts
    latitude = c(40, 41, 42, 43),                   # metadata: excluded
    NDAMS_2009 = c(NA, NA, NA, 3),
    MAJ_DDENS_2009 = c(NA, NA, NA, 1),
    num_water_years = c(20L, 21L, 22L, 23L),
    area_normalized = c(TRUE, TRUE, TRUE, FALSE),
    gage_type = c("Canada", "Canada", "USGS", "USGS"),
    stringsAsFactors = FALSE
  )
}

test_that("table path: NA fraction over present signature columns (production semantics)", {
  df <- .frame()
  out <- compute_qa_flags(df)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), 4)
  expect_setequal(names(out), c("flagged_for_qann_range", "flagged_for_bfi_eckhardt_range",
    "flagged_for_bfi_lynehollick_range", "flagged_for_flashiness_range", "flagged_for_tqmean_range",
    "flagged_for_d50_range", "flagged_for_elasticity_range", "flagged_for_runoff_ratio_range",
    "flagged_for_seasonal_sum", "flagged_for_percentile_order", "flagged_for_timing_order",
    "flagged_for_high_na"))
  # 1/5, 4/5, 3/5, 5/5
  expect_equal(out$flagged_for_high_na, c(FALSE, TRUE, TRUE, TRUE))
  # data.table input behaves identically
  expect_equal(compute_qa_flags(data.table::as.data.table(df))$flagged_for_high_na,
               c(FALSE, TRUE, TRUE, TRUE))
  # the 11 range flags agree with the per-gage list path row by row
  for (i in seq_len(nrow(df))) {
    row <- compute_qa_flags(as.list(df[i, , drop = FALSE]))
    for (nm in setdiff(names(out), "flagged_for_high_na")) {
      expect_identical(unname(out[[nm]][i]), unname(row[[nm]]), info = paste(nm, i))
    }
  }
})

test_that("list path: per-gage semantics over EMITTED keys", {
  sigs <- list(gage_id = "g", Qann_mean = 1, Qann_median = NaN, Q5_pettitt_cp_year = NA_real_,
               latitude = 40, NDAMS_2009 = NA_real_)   # 2 of 3 signature keys NA
  expect_true(compute_qa_flags(sigs)$flagged_for_high_na)
  sigs2 <- list(gage_id = "g", Qann_mean = 1, Qann_median = 1, Q5_pettitt_cp_year = 2000,
                NDAMS_2009 = NA_real_, MAJ_DDENS_2009 = NA_real_)   # metadata NA must not flag
  expect_false(compute_qa_flags(sigs2)$flagged_for_high_na)
})

test_that("metadata-only input never flags", {
  df <- data.frame(gage_id = c("a", "b"), latitude = c(1, 2), NDAMS_2009 = c(NA, NA))
  expect_length(high_na_denominator_columns(names(df)), 0)
  expect_equal(compute_qa_flags(df)$flagged_for_high_na, c(FALSE, FALSE))
})

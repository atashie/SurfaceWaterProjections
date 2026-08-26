# Regression: the flashiness guard must treat a water year with NEGATIVE total
# flow as non-computable (canonical Julia: total_Q <= 0), not just == 0.
# Pre-fix (== 0), a negative-total year produced a negative-denominator R-B
# value that Julia declines — found 2026-08-24 by the port campaign's Phase-1
# floor-transition gate at gage 02244440 (WY1999 total = -78.8 mm).

test_that("negative-total water year is non-computable for flashiness", {
  rows <- do.call(rbind, lapply(list(
    list(wy = 2001, base = 1.0),
    list(wy = 2002, base = -0.2),   # negative total
    list(wy = 2003, base = 1.5)
  ), function(sp) {
    data.frame(water_year = sp$wy,
               Q = sp$base + 0.01 * ((seq_len(60) - 1) %% 5),
               dowy = seq_len(60))
  }))
  expect_lt(sum(rows$Q[rows$water_year == 2002]), 0)

  out <- suppressWarnings(analyze_flashiness_trends(rows))
  # Only 2 computable years remain (< min_rows 3) -> all stats NaN. The old
  # == 0 guard kept 3 "computable" years and returned finite statistics.
  expect_true(is.na(out[["flashinessRB_mean"]]))
})

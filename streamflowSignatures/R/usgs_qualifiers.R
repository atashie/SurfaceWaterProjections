# usgs_qualifiers.R — the USGS NWIS daily-value acceptance rule of the ingestion step
# (decision 1d of docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md, approved by the
# user 2026-10-26).
#
# A daily value is KEPT when all three hold:
#   (1) its qualification string starts with the status letter A (approved) or
#       P (provisional) — the first token of the string;
#   (2) the value is a finite number;
#   (3) the value is not the service's no-data sentinel (-999999, the `noDataValue` of the
#       NWIS daily-value service).
# Remark codes after the status letter (e estimated, R revised, < and > censored at a
# reporting limit, [4] "no matching qualifier", Ice, Bkw, ZFL, Rat, Dis, Eqp, Fld, Ssn,
# Mnt, ***) NEVER mask a value at storage: the full code string is kept verbatim in the
# `flag` column and downstream analysis decides. Before 2026-10-26 the rule kept only the
# exact strings "A", "A e", "P", "P e" and set 193,542 approved days (A R, A <, A >, A [4])
# on 987 compiled gages to NA (docs/plans/2026-09-04-filtering-stage-census.md §10.1).
#
# Pure base-R functions with no package dependency, unit-tested by
# R/tests/test_qualifier_codes.R (`Rscript R/tests/test_qualifier_codes.R` from the repo
# root). The Python splice tool docs/benchmarks/splice_nwis_qualifiers.py implements the
# same rule to patch the February 2026 parquet; keep the two in step.

USGS_NWIS_NO_DATA_VALUE <- -999999

USGS_STATUS_LETTERS <- c("A", "P")

#' Status letter of a NWIS qualification string: its first token ("" when absent).
#' Tokens are separated by whitespace, ':' or ',' (dataRetrieval serves "A e";
#' the WaterML JSON serves a qualifier array that the fetch tool joins with spaces).
usgs_qualifier_status <- function(flag) {
  f <- as.character(flag)
  f[is.na(f)] <- ""
  f <- trimws(f)
  sub("[[:space:]:,].*$", "", f)
}

#' TRUE where the value is accepted under the rule above (vectorised; same length as
#' `flag` and `value`). Non-numeric values are not finite and are never kept.
usgs_qualifier_keep <- function(flag, value, no_data_value = USGS_NWIS_NO_DATA_VALUE) {
  if (length(flag) != length(value)) {
    stop("usgs_qualifier_keep: flag and value must have the same length")
  }
  v <- suppressWarnings(as.numeric(value))
  status_ok <- usgs_qualifier_status(flag) %in% USGS_STATUS_LETTERS
  value_ok <- is.finite(v) & !(v == no_data_value)
  value_ok[is.na(value_ok)] <- FALSE
  status_ok & value_ok
}

#' The numeric value vector with every non-accepted entry set to NA_real_; accepted
#' entries are returned unchanged (as.numeric of the input).
usgs_qualifier_mask <- function(flag, value, no_data_value = USGS_NWIS_NO_DATA_VALUE) {
  v <- suppressWarnings(as.numeric(value))
  v[!usgs_qualifier_keep(flag, value, no_data_value)] <- NA_real_
  v
}

#' Per-code tally for the input audit ("flag x is.na(Q)"): one row per distinct code with
#' the number of days kept, masked by the rule although a finite value was served, and
#' served without a usable value (NA, non-numeric or the no-data sentinel).
usgs_qualifier_tally <- function(flag, value, no_data_value = USGS_NWIS_NO_DATA_VALUE) {
  f <- as.character(flag)
  f[is.na(f)] <- "<NA>"
  v <- suppressWarnings(as.numeric(value))
  has_value <- is.finite(v) & !(v == no_data_value)
  has_value[is.na(has_value)] <- FALSE
  kept <- usgs_qualifier_keep(flag, value, no_data_value)
  codes <- sort(unique(f))
  out <- data.frame(
    flag = codes,
    days = vapply(codes, function(cd) sum(f == cd), integer(1)),
    kept = vapply(codes, function(cd) sum(kept[f == cd]), integer(1)),
    masked_with_value = vapply(codes, function(cd) sum(has_value[f == cd] & !kept[f == cd]), integer(1)),
    no_value = vapply(codes, function(cd) sum(!has_value[f == cd]), integer(1)),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

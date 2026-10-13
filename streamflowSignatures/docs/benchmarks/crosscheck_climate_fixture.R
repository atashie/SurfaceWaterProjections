# Compare the rpkg climate metrics with canonical Julia on the shared fixture.
#   Rscript docs/benchmarks/crosscheck_climate_fixture.R OUT_DIR [rtol] [atol] [mk_rtol]
# The annual values must be bit-identical; summary statistics within rtol/atol. The
# Mann-Kendall-derived ones (_mk_rho, _mk_pval, _pettitt_*_mk_pval) take mk_rtol, which
# defaults to rtol since rpkg's Mann-Kendall became a double-precision port of Julia's
# formula (2026-10-12; the earlier Kendall::MannKendall backend needed ~1e-6).
# Reads OUT_DIR/climate_fixture_input.csv (written by crosscheck_climate_fixture.jl)
# and compares every one of the 39 x 16 values with OUT_DIR/climate_fixture_julia.json:
# |d| <= atol + rtol * scale; NaN must match NaN. Exit status 1 on any mismatch.
suppressPackageStartupMessages(library(streamflowsignatures))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
rtol <- if (length(args) >= 2) as.numeric(args[2]) else 4e-14
atol <- if (length(args) >= 3) as.numeric(args[3]) else 1e-12
mk_rtol <- if (length(args) >= 4) as.numeric(args[4]) else rtol
raw <- as.data.frame(arrow::read_parquet(file.path(out_dir, "climate_fixture_input.parquet")))  # exact doubles
raw$gage_id <- as.character(raw$gage_id)
raw$date <- as.Date(raw$date)
ref <- jsonlite::fromJSON(file.path(out_dir, "climate_fixture_julia.json"))
pkg_env <- streamflowsignatures:::pkg_env
df <- add_water_year_columns(raw)
pp <- preprocess_daily_data(df)
stopifnot(identical(as.integer(pp$valid_temp_years), as.integer(ref$valid_temp_years)))
stopifnot(identical(as.integer(pp$valid_climate_years), as.integer(ref$valid_climate_years)))
data <- as.data.frame(pp$data)
cdata <- data[data$water_year %in% pp$valid_climate_years, ]
tdata <- data[data$water_year %in% pp$valid_temp_years, ]
cp <- list(start_year = pkg_env$cp_start_water_year, end_year = pkg_env$cp_end_water_year,
           min_total_obs = pkg_env$cp_min_total_obs, min_segment_obs = pkg_env$cp_min_segment_obs)
coll <- annual_collector()
res <- c(calculate_precipitation_metrics(cdata, trend_completeness = pkg_env$na_trend_min_fraction,
                                         decade_completeness = pkg_env$na_decade_min_fraction,
                                         min_values_for_stats = pkg_env$min_values_for_stats, changepoint = cp,
                                         collector = coll),
         calculate_temperature_metrics(tdata, trend_completeness = pkg_env$na_trend_min_fraction,
                                       decade_completeness = pkg_env$na_decade_min_fraction,
                                       min_values_for_stats = pkg_env$min_values_for_stats, changepoint = cp,
                                       collector = coll))
ann_ref <- as.data.frame(arrow::read_parquet(file.path(out_dir, "climate_fixture_julia_annual.parquet")))
ann_ref$signature <- as.character(ann_ref$signature); ann_ref$water_year <- as.integer(ann_ref$water_year)
ann_r <- as.data.frame(collector_drain(coll))
m <- merge(ann_ref, ann_r, by = c("signature", "water_year"), all = TRUE, suffixes = c("_jl", "_r"))
unmatched <- sum(is.na(m$value_jl) & is.na(m$value_r) & FALSE)  # placeholder for clarity
missing_rows <- nrow(m) - nrow(ann_ref)
nan_mismatch <- sum(is.na(m$value_jl) != is.na(m$value_r))
both <- m[!is.na(m$value_jl) & !is.na(m$value_r), ]
ann_diff <- sum(both$value_jl != both$value_r)
cat(sprintf("annual values: %d Julia rows, %d rpkg rows, unmatched %d, NaN-pattern mismatches %d, non-identical values %d\n",
            nrow(ann_ref), nrow(ann_r), missing_rows + (nrow(m) - nrow(ann_r)), nan_mismatch, ann_diff))
jl <- ref$values
missing <- setdiff(names(jl), names(res)); extra <- setdiff(names(res), names(jl))
worst <- 0; worst_key <- NA; bad <- character(0)
for (k in names(jl)) {
  if (!k %in% names(res)) next
  vj <- jl[[k]]; vr <- res[[k]]
  vj_nan <- identical(vj, "NaN"); vr_nan <- is.na(vr)
  if (vj_nan || vr_nan) { if (vj_nan != vr_nan) bad <- c(bad, sprintf("%s nan-mismatch julia=%s r=%s", k, vj, vr)); next }
  vj <- as.numeric(vj); vr <- as.numeric(vr)
  scale <- max(abs(vj), abs(vr), 1e-300); rel <- abs(vj - vr) / scale
  if (rel > worst) { worst <- rel; worst_key <- k }
  tol <- if (grepl("(_mk_rho|_mk_pval|_pettitt_pre_mk_pval|_pettitt_post_mk_pval)$", k)) mk_rtol else rtol
  if (abs(vj - vr) > atol + tol * scale) bad <- c(bad, sprintf("%s julia=%.17g r=%.17g rel=%.3e", k, vj, vr, rel))
}
cat(sprintf("compared %d Julia values; rpkg has %d; missing in rpkg: %d; extra: %d\n", length(jl), length(res), length(missing), length(extra)))
cat(sprintf("worst relative difference: %.3e at %s\n", worst, worst_key))
for (b in head(bad, 30)) cat("  MISMATCH", b, "\n")
ok <- length(missing) == 0 && length(extra) == 0 && length(bad) == 0 && nan_mismatch == 0 && ann_diff == 0 && nrow(m) == nrow(ann_ref)
cat("CROSSCHECK", if (ok) "PASS" else "FAIL", "\n")
quit(status = if (ok) 0 else 1)

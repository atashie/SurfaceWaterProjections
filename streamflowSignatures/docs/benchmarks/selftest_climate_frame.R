# selftest_climate_frame.R — the rpkg side of the climate-frame fixture self-test
# (climate_frame.R; plan §10.2 tests (i)–(iv), 2026-10-26).
#   Rscript docs/benchmarks/selftest_climate_frame.R OUT_DIR [rtol] [summary_rtol]
# Reads the fixture that selftest_climate_frame.jl wrote to OUT_DIR, builds the frame in both
# modes with the R helper (the rpkg runner's column set: gage_id, date, Q), runs the rpkg
# preprocessor and the production orchestrator call, repeats the fixture checks and compares
# the annual values and the summary values with canonical Julia per mode (annual NaN
# patterns agree and values within rtol — bit-identical count reported; summaries within
# summary_rtol with gate-state agreement). Rank statistics (Mann-Kendall, Spearman, the
# Pettitt fields) of a base whose annual series holds a tie or a near-tie (two values within
# 1e-9 relative) AND whose annual values are not bit-identical in the two languages are
# reported, not gated: a last-bit difference then reorders ranks — the named residual class
# of docs/CROSS_LANGUAGE_STATUS.md. Exit status 1 on any failed check.
suppressPackageStartupMessages({ library(streamflowsignatures); library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
rtol <- if (length(args) >= 2) as.numeric(args[2]) else 1e-9
summary_rtol <- if (length(args) >= 3) as.numeric(args[3]) else 1e-6
atol <- 1e-12
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
source(file.path(script_dir, "climate_frame.R"))
pkg_env <- streamflowsignatures:::pkg_env
GID <- "CF1"
n_pass <- 0L; n_fail <- 0L
check <- function(label, cond) {
  if (isTRUE(cond)) { n_pass <<- n_pass + 1L; cat("  PASS:", label, "\n") }
  else { n_fail <<- n_fail + 1L; cat("  FAIL:", label, "\n") }
}
# default mode (user decision 2026-10-28): "outer" when the variable is unset
.cf_saved <- Sys.getenv("STREAMFLOW_CLIMATE_FRAME", unset = NA)
Sys.unsetenv("STREAMFLOW_CLIMATE_FRAME")
check("climate_frame_mode() is \"outer\" when STREAMFLOW_CLIMATE_FRAME is unset (default since 2026-10-28)",
      identical(climate_frame_mode(), "outer"))
if (!is.na(.cf_saved)) Sys.setenv(STREAMFLOW_CLIMATE_FRAME = .cf_saved)
same <- function(a, b) (is.na(a) && is.na(b)) || (!is.na(a) && !is.na(b) && a == b)
SUFFIXES <- c("_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval",
              "_mean", "_median", "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
              "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
              "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval")
RANK_SUFFIXES <- setdiff(SUFFIXES, c("_senn_slp", "_linear_slp", "_mean", "_median"))
split_key <- function(k) {
  for (s in SUFFIXES[order(-nchar(SUFFIXES))]) if (endsWith(k, s)) return(list(base = substr(k, 1, nchar(k) - nchar(s)), suffix = s))
  list(base = k, suffix = NA_character_)
}
same_vec <- function(a, b) length(a) == length(b) && all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
close_to <- function(a, b, at = 1e-9, rt = 1e-9) abs(a - b) <= at + rt * max(abs(a), abs(b))

sf <- as.data.table(arrow::read_parquet(file.path(out, "climate_frame_fixture_streamflow.parquet")))
cl <- as.data.table(arrow::read_parquet(file.path(out, "climate_frame_fixture_climate.parquet")))
sf[, gage_id := as.character(gage_id)]; sf[, date := as.Date(date)]
cl[, gage_id := as.character(gage_id)]; cl[, date := as.Date(date)]
# the rpkg runner reads gage_id, Date, Q and re-derives the water-year columns
sf <- add_water_year_columns(sf[, .(gage_id, date, Q)])
climate_merge <- cl[, .(date, PPT, SWE, TMIN, TMAX)]
ref <- list(left = jsonlite::fromJSON(file.path(out, "climate_frame_fixture_summary_left.json")),
            outer = jsonlite::fromJSON(file.path(out, "climate_frame_fixture_summary_outer.json")))
ann_ref <- list(left = as.data.frame(arrow::read_parquet(file.path(out, "climate_frame_fixture_annual_left.parquet"))),
                outer = as.data.frame(arrow::read_parquet(file.path(out, "climate_frame_fixture_annual_outer.parquet"))))

frames <- list(left = frame_with_climate(copy(sf), climate_merge, "left"),
               outer = frame_with_climate(copy(sf), climate_merge, "outer"))
pp <- list(left = preprocess_daily_data(copy(frames$left)), outer = preprocess_daily_data(copy(frames$outer)))
L <- frames$left; O <- frames$outer

cat("\n(ii) frame composition\n")
check("left frame = the streamflow rows", nrow(L) == nrow(sf) && setequal(L$date, sf$date))
added <- sort(setdiff(O$date, L$date))
check(sprintf("outer adds exactly the %d absent dates of the Julia reference, date for date", ref$outer$n_added_rows),
      identical(format(as.Date(added)), as.character(ref$outer$added_dates)) && nrow(O) == ref$outer$n_frame_rows)
check("left adds nothing (the Julia reference lists no added date)", length(ref$left$added_dates) == 0)
# the full parquet schema (gage_id, date, Q, year, doy, flag + the water-year columns): the
# added rows carry year and doy as well (review 2026-10-26 finding 7)
sf_full <- as.data.table(arrow::read_parquet(file.path(out, "climate_frame_fixture_streamflow.parquet")))
sf_full[, gage_id := as.character(gage_id)]; sf_full[, date := as.Date(date)]
o_full <- frame_with_climate(copy(sf_full), climate_merge, "outer")
added_rows <- o_full[date %in% as.Date(added)]
check("full schema: added rows carry year and doy from the date, Q and flag NA, the same column classes",
      nrow(added_rows) == length(added) && all(added_rows$year == as.integer(format(added_rows$date, "%Y"))) &&
      all(added_rows$doy == as.integer(format(added_rows$date, "%j"))) && all(is.na(added_rows$Q)) && all(is.na(added_rows$flag)) &&
      identical(sapply(o_full, class), sapply(merge(copy(sf_full), climate_merge, by = "date", all.x = TRUE), class)))
check("full schema: the rows on streamflow dates equal the former merge expression",
      isTRUE(all.equal(as.data.frame(o_full[!date %in% as.Date(added)]), as.data.frame(merge(copy(sf_full), climate_merge, by = "date", all.x = TRUE)), check.attributes = FALSE)))
check("a POSIXct streamflow date column is refused in both modes",
      all(vapply(c("left", "outer"), function(m) inherits(tryCatch(frame_with_climate(copy(sf)[, date := as.POSIXct(date)], climate_merge, m), error = function(e) e), "error"), logical(1))))
cl_dup <- rbind(climate_merge, climate_merge[date == as.Date(added[1])])
check("an identical duplicate climate row on an absent date is added once",
      sum(frame_with_climate(copy(sf), cl_dup, "outer")$date == as.Date(added[1])) == 1)
cl_conf <- copy(cl_dup); cl_conf[.N, PPT := PPT + 1]
check("a conflicting duplicate climate row on an absent date is an error",
      inherits(tryCatch(frame_with_climate(copy(sf), cl_conf, "outer"), error = function(e) e), "error"))
check("no added row in the absent water year 2005", !(2005 %in% O$water_year))
check("no added row before the record start or after its end", min(O$date) == min(sf$date) && max(O$date) == max(sf$date))
check("the water-year set is identical between modes", setequal(L$water_year, O$water_year))
check("Q is NA on every added row", all(is.na(O$Q[O$date %in% as.Date(added)])))
check("the rows on streamflow dates of the outer frame equal the left frame",
      isTRUE(all.equal(as.data.frame(O[!date %in% as.Date(added), names(L), with = FALSE]), as.data.frame(L), check.attributes = FALSE)))

cat("\n(ii) the runner's window denominator (run_rpkg_benchmark.R: total_possible from the frame's water years)\n")
window_denominator <- function(frame, p, start_wy, end_wy) {
  all_wy <- unique(frame$water_year); all_wy <- all_wy[!is.na(all_wy)]
  wy_in_range <- all_wy[all_wy >= start_wy & all_wy <= end_wy]
  if (length(wy_in_range) == 0) return(list(total_possible = 0L, n_valid = 0L))
  list(total_possible = as.integer(max(wy_in_range) - start_wy + 1), n_valid = as.integer(sum(p$valid_years >= start_wy & p$valid_years <= end_wy)))
}
for (m in c("left", "outer")) {
  ok <- TRUE
  for (k in names(ref[[m]]$windows)) {
    b <- as.integer(strsplit(k, "-")[[1]]); got <- window_denominator(frames[[m]], pp[[m]], b[1], b[2]); want <- ref[[m]]$windows[[k]]
    if (got$total_possible != want$total_possible || got$n_valid != want$n_valid) ok <- FALSE
  }
  check(sprintf("%s: total_possible and the valid-year count equal Julia for every window", m), ok)
}

cat("\n(ii) Q-only preprocessor outputs identical between modes\n")
pl <- pp$left; po <- pp$outer
check("valid_years identical", identical(as.integer(pl$valid_years), as.integer(po$valid_years)))
check("rejected_years identical, WY 2005 absent", isTRUE(all.equal(pl$rejected_years, po$rejected_years)) &&
      !(2005 %in% as.data.frame(pl$rejected_years)$water_year))
check("seasonal_flags identical", isTRUE(all.equal(pl$seasonal_flags, po$seasonal_flags)))
check("diagnostics identical", isTRUE(all.equal(pl$diagnostics, po$diagnostics)))
dl <- as.data.frame(pl$data); do_ <- as.data.frame(po$data)
dcol <- if ("Date" %in% names(dl)) "Date" else "date"
check("preprocessed Q identical (bitwise, NA == NA)", nrow(dl) == nrow(do_) && identical(dl[[dcol]], do_[[dcol]]) && same_vec(dl$Q, do_$Q))
cp <- list(start_year = pkg_env$cp_start_water_year, end_year = pkg_env$cp_end_water_year,
           min_total_obs = pkg_env$cp_min_total_obs, min_segment_obs = pkg_env$cp_min_segment_obs)
qonly <- lapply(pp, function(p) calculate_all_signatures(
  p$data, has_climate = FALSE, seasonal_flags = p$seasonal_flags,
  trend_completeness = pkg_env$na_trend_min_fraction, decade_completeness = pkg_env$na_decade_min_fraction,
  min_values_for_stats = pkg_env$min_values_for_stats, changepoint = cp, gage_id = GID))
qkeys <- sort(names(qonly$left))
check(sprintf("Q-only orchestrator: the same %d keys in both modes (Julia: %d)", length(qkeys), ref$left$checks$n_qonly_keys),
      identical(qkeys, sort(names(qonly$outer))) && length(qkeys) == ref$left$checks$n_qonly_keys)
check("Q-only orchestrator: every value bit-identical between modes",
      all(vapply(qkeys, function(k) same(as.numeric(qonly$left[[k]]), as.numeric(qonly$outer[[k]])), logical(1))))

run_full <- function(p) {   # the production orchestrator call of run_rpkg_benchmark.R (new path)
  data <- p$data
  cdata <- data[data$water_year %in% p$valid_climate_years, , drop = FALSE]
  sdata <- data[data$water_year %in% p$valid_swe_years, , drop = FALSE]
  tdata <- data[data$water_year %in% p$valid_temp_years, , drop = FALSE]
  coll <- annual_collector()
  res <- calculate_all_signatures(
    data, has_climate = TRUE, seasonal_flags = p$seasonal_flags,
    trend_completeness = pkg_env$na_trend_min_fraction, decade_completeness = pkg_env$na_decade_min_fraction,
    min_values_for_stats = pkg_env$min_values_for_stats, climate_data = cdata, area_normalized = TRUE,
    changepoint = cp, collector = coll, snow_data = sdata, snow_climate_years = p$valid_climate_years,
    gage_id = GID, temp_data = tdata)
  list(ann = as.data.frame(collector_drain(coll)), res = res)
}
full <- lapply(pp, run_full)

cat("\n(i) + (iii) against the Julia reference numbers\n")
ppt_at <- function(m, d) { df <- as.data.frame(pp[[m]]$data); v <- df$PPT[as.Date(df[[dcol]]) == as.Date(d)]; if (length(v)) v[1] else NA_real_ }
prcp <- function(m, wy) { x <- full[[m]]$ann; v <- x$value[x$signature == "prcp_annual" & x$water_year == wy]; if (length(v)) v[1] else NA_real_ }
num_or_na <- function(v) if (is.character(v)) NA_real_ else as.numeric(v)
for (m in c("left", "outer")) {
  ck <- ref[[m]]$checks
  check(sprintf("%s: preprocessed PPT on 2000-12-30 / 12-31 equals Julia", m),
        same(ppt_at(m, "2000-12-30"), num_or_na(ck$ppt_2000_12_30)) && same(ppt_at(m, "2000-12-31"), num_or_na(ck$ppt_2000_12_31)))
  check(sprintf("%s: prcp_annual WY 1995 equals Julia (%.6f)", m, num_or_na(ck$prcp_annual_1995)), same(prcp(m, 1995), num_or_na(ck$prcp_annual_1995)))
  for (key in c("valid_years", "valid_climate_years", "valid_swe_years", "valid_temp_years"))
    check(sprintf("%s: %s equal Julia", m, key), identical(as.integer(pp[[m]][[key]]), as.integer(ref[[m]][[key]])))
}
ck <- ref$outer$checks
check("outer: prcp_annual WY 1995 = Daymet sum; left = Daymet sum - interpolation residual",
      close_to(prcp("outer", 1995), ck$daymet_sum_1995) && close_to(prcp("left", 1995), ck$daymet_sum_1995 - ck$interp_resid_1995))
check("WY 2009 and 2013 climate-valid only in outer mode",
      identical(as.integer(sort(setdiff(pp$outer$valid_climate_years, pp$left$valid_climate_years))), c(2009L, 2013L)))

cat("\n(iv) cross-language: annual values and summaries vs Julia, per mode\n")
for (m in c("left", "outer")) {
  r <- ann_ref[[m]]; r$signature <- as.character(r$signature); r$water_year <- as.integer(r$water_year)
  x <- full[[m]]$ann; x$signature <- as.character(x$signature); x$water_year <- as.integer(x$water_year)
  mm <- merge(r, x, by = c("signature", "water_year"), all = TRUE, suffixes = c("_jl", "_r"))
  unmatched <- (nrow(mm) - nrow(r)) + (nrow(mm) - nrow(x))
  nan_mis <- sum(is.na(mm$value_jl) != is.na(mm$value_r))
  both <- mm[!is.na(mm$value_jl) & !is.na(mm$value_r), ]
  n_bit <- sum(both$value_jl == both$value_r)
  scale <- pmax(abs(both$value_jl), abs(both$value_r), 1e-300)
  rel <- abs(both$value_jl - both$value_r) / scale
  within <- abs(both$value_jl - both$value_r) <= atol + rtol * scale
  wi <- which.max(rel)
  cat(sprintf("  %s: %d Julia rows, %d rpkg rows, unmatched %d, NaN-pattern mismatches %d, bit-identical %d/%d, worst rel %.3e at %s WY %d\n",
              m, nrow(r), nrow(x), unmatched, nan_mis, n_bit, nrow(both), rel[wi], both$signature[wi], both$water_year[wi]))
  check(sprintf("%s: annual values match Julia (keys, NaN pattern, rel <= %g)", m, rtol), unmatched == 0 && nan_mis == 0 && all(within))
  # tie-sensitive bases: a tie or near-tie in Julia's annual series and annual values that
  # are not bit-identical in the port
  tie_sensitive <- character(0)
  for (sig in unique(both$signature)) {
    g <- both[both$signature == sig, ]
    if (nrow(g) < 2) next
    vs <- sort(g$value_jl); gaps <- diff(vs) / pmax(abs(vs[-1]), 1e-300)
    if (any(gaps <= 1e-9) && !same_vec(g$value_jl, g$value_r)) tie_sensitive <- c(tie_sensitive, sig)
  }
  jl <- ref[[m]]$values; rs <- full[[m]]$res
  missing <- setdiff(names(jl), names(rs)); extra <- setdiff(names(rs), names(jl))
  bad <- character(0); diag <- character(0); worst <- 0; worst_key <- NA
  for (k in names(jl)) {
    if (!k %in% names(rs)) next
    vj <- jl[[k]]; vr <- as.numeric(rs[[k]])
    sk <- split_key(k); sensitive <- !is.na(sk$suffix) && sk$suffix %in% RANK_SUFFIXES && sk$base %in% tie_sensitive
    jn <- identical(vj, "NaN"); rn <- is.na(vr)
    if (jn || rn) {
      if (jn != rn) bad <- c(bad, sprintf("%s nan-mismatch julia=%s r=%s", k, vj, vr))   # a gate-state difference is never exempt
      next
    }
    vj <- as.numeric(vj); sc <- max(abs(vj), abs(vr), 1e-300); rl <- abs(vj - vr) / sc
    if (abs(vj - vr) > atol + summary_rtol * sc) {
      msg <- sprintf("%s julia=%.17g r=%.17g rel=%.3e", k, vj, vr, rl); if (sensitive) diag <- c(diag, msg) else bad <- c(bad, msg)
    } else if (rl > worst) { worst <- rl; worst_key <- k }
  }
  cat(sprintf("  %s: %d Julia summary values, rpkg %d, missing %d, extra %d, worst gated rel %.3e at %s, mismatches %d; tie-sensitive rank statistics reported, not gated: %d on %s\n",
              m, length(jl), length(rs), length(missing), length(extra), worst, worst_key, length(bad), length(diag), paste(tie_sensitive, collapse = ", ")))
  for (b in head(bad, 20)) cat("    MISMATCH", b, "\n")
  for (b in head(diag, 10)) cat("    tie-sensitive", b, "\n")
  check(sprintf("%s: summary values match Julia (gate state always; rel <= %g; numeric rank statistics on tie-bearing series excepted)", m, summary_rtol),
        length(missing) == 0 && length(extra) == 0 && length(bad) == 0)
  # the exemption is justified only if the port's statistics code reproduces Julia's rank
  # statistics when fed JULIA's annual values (review 2026-10-26, finding 6)
  n_cmp <- 0L; n_skip <- 0L; bad2 <- character(0)
  for (base in sort(tie_sensitive)) {
    g <- r[r$signature == base, ]
    df <- data.frame(water_year = as.integer(g$water_year)); df[[base]] <- as.numeric(g$value)
    st <- generate_stats(df, value_cols = base, trend_completeness = pkg_env$na_trend_min_fraction,
                         decade_completeness = pkg_env$na_decade_min_fraction, changepoint = cp,
                         min_values_for_stats = pkg_env$min_values_for_stats)
    for (suf in sort(RANK_SUFFIXES)) {
      k <- paste0(base, suf)
      if (!k %in% names(jl) || !k %in% names(st)) next
      vj <- jl[[k]]; vr <- as.numeric(st[[k]])
      jn <- identical(vj, "NaN"); rn <- is.na(vr)
      if (jn || rn) { n_skip <- n_skip + 1L; next }
      n_cmp <- n_cmp + 1L; vj <- as.numeric(vj)
      if (abs(vj - vr) > 1e-12 + 1e-9 * max(abs(vj), abs(vr))) bad2 <- c(bad2, sprintf("%s julia=%.17g r=%.17g", k, vj, vr))
    }
  }
  for (b in head(bad2, 10)) cat("    STATS MISMATCH on Julia's values", b, "\n")
  check(sprintf("%s: the port's rank statistics reproduce Julia's on Julia's own annual values for the tie-sensitive bases (%d compared, %d gated differently in the direct call)", m, n_cmp, n_skip),
        length(bad2) == 0)
}
cat(sprintf("\n%d passed, %d failed\n", n_pass, n_fail))
quit(status = if (n_fail == 0) 0 else 1)

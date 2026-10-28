# climate_frame.R — the production climate-frame helper of the rpkg runner
# (Phase A′ of docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md §10.2; decisions
# 12–13, user 2026-10-26). Mirrors climate_frame.jl (canonical Julia runner) and
# climate_frame.py; self-test: selftest_climate_frame.R. Needs data.table and the rpkg's
# exported add_water_year_columns() (library(streamflowsignatures) or load_all first).
#
# The runner merges each gage's Daymet rows onto that gage's streamflow ROWS before the
# package builds the water-year grid, so a date with no streamflow row reached
# preprocess_daily_data() without climate and was interpolated (CHANGELOG Known Issue
# "climate rides on the streamflow dates", 2026-10-16). frame_with_climate() builds the
# frame under the runtime switch STREAMFLOW_CLIMATE_FRAME:
#   left   today's merge(all.x = TRUE), bit for bit — the attribution baseline;
#   outer  that merge PLUS the gage's Daymet rows on dates with no streamflow row,
#          restricted to the streamflow record span AND to water years that hold at least
#          one streamflow row (whole interior water years without streamflow must not
#          enter the preprocessor's year loop — Codex round 3).
# Both frames must carry a Date column `date` (asserted). Added rows carry gage_id, date,
# water_year / month / dowy (add_water_year_columns; year and doy from the date when the
# frame has them), the climate columns and NA for Q and every other column; the result is
# ordered by date, as data.table's merge() orders the left frame. Default mode `outer` since
# 2026-10-28 (user decision; the version-3 products are built with it); `left` reproduces every
# product before the qualifier rerun and the S0 / S1 attribution baselines.

CLIMATE_FRAME_MODES <- c("left", "outer")

climate_frame_mode <- function() {
  v <- tolower(trimws(Sys.getenv("STREAMFLOW_CLIMATE_FRAME", unset = "outer")))   # default outer since 2026-10-28
  if (!v %in% CLIMATE_FRAME_MODES) {
    stop("STREAMFLOW_CLIMATE_FRAME must be one of ", paste(CLIMATE_FRAME_MODES, collapse = ", "),
         "; got \"", v, "\"")
  }
  v
}

# gage_data: one gage's streamflow rows (data.table: gage_id, date, Q, water_year, month,
# dowy, ...); climate_merge: that gage's climate rows (date, PPT[, SWE, TMIN, TMAX]);
# mode: "left" or "outer".
frame_with_climate <- function(gage_data, climate_merge, mode) {
  # both frames must carry a Date column `date`: a POSIXct or character column would match
  # no climate row in the merge and silently yield a climate-free frame (review 2026-10-26)
  for (nm in c("gage_data", "climate_merge")) {
    df <- get(nm)
    if (!"date" %in% names(df)) stop("frame_with_climate: the ", nm, " frame has no date column")
    if (!inherits(df$date, "Date")) stop("frame_with_climate: the ", nm, " frame's date column is ",
                                         paste(class(df$date), collapse = "/"), ", expected Date")
  }
  if (mode == "left") {
    return(merge(gage_data, climate_merge, by = "date", all.x = TRUE))
  }
  if (mode != "outer") stop("frame_with_climate: unknown mode \"", mode, "\"")
  joined <- merge(gage_data, climate_merge, by = "date", all.x = TRUE)
  if (nrow(gage_data) == 0 || nrow(climate_merge) == 0) return(joined)
  if (!"water_year" %in% names(gage_data)) {
    stop("frame_with_climate: the streamflow frame has no water_year column")
  }

  sf_dates <- as.Date(gage_data$date)
  lo <- min(sf_dates)
  hi <- max(sf_dates)
  wy_set <- unique(gage_data$water_year[!is.na(gage_data$water_year)])
  c_dates <- as.Date(climate_merge$date)
  in_span <- c_dates >= lo & c_dates <= hi
  if (!any(in_span)) return(joined)

  cand <- data.table::copy(data.table::as.data.table(climate_merge)[in_span])
  # duplicate climate dates: identical duplicate rows collapse to one, conflicting ones are
  # an error — an appended duplicate would put two rows on one grid day
  if (anyDuplicated(cand$date) > 0) {
    cand <- unique(cand)
    if (anyDuplicated(cand$date) > 0) {
      stop("frame_with_climate: conflicting duplicate climate rows for gage ", gage_data$gage_id[1],
           " on ", format(cand$date[duplicated(cand$date)][1]))
    }
  }
  cand <- add_water_year_columns(cand)          # water_year, month, dowy (in place)
  keep <- cand$water_year %in% wy_set & !(cand$date %in% sf_dates)
  if (!any(keep)) return(joined)
  add <- cand[keep]

  add[, gage_id := gage_data$gage_id[1]]
  # the other temporal columns of the frame, when it has them (year, doy as the parquet
  # stores them), so the frame is uniform across the three languages
  same_type <- function(x, template) if (is.double(template)) as.double(x) else if (is.integer(template)) as.integer(x) else x
  if ("year" %in% names(joined)) add[, year := same_type(as.integer(format(date, "%Y")), joined$year)]
  if ("doy" %in% names(joined)) add[, doy := same_type(as.integer(format(date, "%j")), joined$doy)]
  for (cn in setdiff(names(joined), names(add))) add[, (cn) := NA]
  add <- add[, names(joined), with = FALSE]
  out <- rbind(joined, add)
  data.table::setorder(out, date)
  out
}

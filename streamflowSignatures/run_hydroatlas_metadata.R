################################################################################
# run_hydroatlas_metadata.R
#
# Build a per-gage WATERSHED hydro-geophysical metadata file from HydroATLAS /
# BasinATLAS v10, sitting alongside the streamflow signatures (keyed by gage_id).
#
# Pipeline: gage point -> outlet HydroBASINS unit (vectorized spatial join) ->
# upstream member set (cached RDS + NEXT_DOWN BFS) -> hybrid per-attribute
# aggregation (prefer HydroATLAS _u; area-weight _s where no _u) -> per gage.
#
# Usage:
#   Rscript run_hydroatlas_metadata.R                 # full run (success gages) -> data_out/
#   Rscript run_hydroatlas_metadata.R --subset 40     # dry-run, 40 spread gages -> test_output/
#   Rscript run_hydroatlas_metadata.R --status all    # all gages with valid coords
#   Rscript run_hydroatlas_metadata.R --metadata <path/to/combined_watershed_metadata.csv>
#
# Prerequisites: HYDROATLAS_GPKG (config.R), upstream_hydrobasins.rds (repo root),
#   a combined_watershed_metadata*.csv (auto-located: PARQUET_DATA_DIR,
#   combined_streamflow_output/, data_out/ — newest by mtime).
################################################################################

t0 <- Sys.time()
suppressMessages({ library(data.table); library(sf); library(arrow) })
source("config.R")
source("R/aggregate_hydroatlas_metadata.R")

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && i[1] < length(args)) args[i[1] + 1L] else default
}
subset_n <- suppressWarnings(as.integer(get_arg("--subset", NA)))
status_filter <- get_arg("--status", "success")
metadata_override <- get_arg("--metadata", NA)

is_dry  <- !is.na(subset_n)
out_dir <- if (is_dry) "test_output" else "data_out"
date_tag <- tolower(format(Sys.Date(), "%d%b%Y"))
if (is_dry) date_tag <- paste0(date_tag, "_dryrun", subset_n)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("========== HYDROATLAS WATERSHED METADATA ==========\n")
cat("Start:", format(t0, "%Y-%m-%d %H:%M:%S"),
    "| mode:", if (is_dry) paste0("DRY-RUN(", subset_n, ")") else "FULL",
    "| status:", status_filter, "| out:", out_dir, "\n")

# ---- locate metadata ----
if (!is.na(metadata_override)) {
  metadata_path <- metadata_override
} else {
  search_dirs <- c(PARQUET_DATA_DIR, "combined_streamflow_output", "data_out")
  cands <- unlist(lapply(search_dirs, function(d)
    if (dir.exists(d)) list.files(d, "^combined_watershed_metadata.*\\.csv$", full.names = TRUE)
    else character(0)))
  if (!length(cands)) stop("No combined_watershed_metadata*.csv found in: ",
                           paste(search_dirs, collapse = ", "))
  # Prefer the date embedded in the filename (e.g. 19feb2026); fall back to mtime
  # (copy timestamps on D: can be misleading).
  fdate <- vapply(basename(cands), function(b) {
    g <- regmatches(b, regexpr("[0-9]{1,2}[a-z]{3}[0-9]{4}", b, ignore.case = TRUE))
    if (length(g)) suppressWarnings(as.numeric(as.Date(g, format = "%d%b%Y"))) else NA_real_
  }, numeric(1))
  metadata_path <- if (any(!is.na(fdate))) cands[which.max(fdate)] else cands[which.max(file.mtime(cands))]
}
cat("Metadata:", metadata_path, "\n")
stopifnot(file.exists(metadata_path), file.exists(HYDROATLAS_GPKG))

md <- fread(metadata_path, colClasses = c(gage_id = "character"))
if (status_filter != "all" && "processing_status" %in% names(md)) {
  md <- md[processing_status == status_filter]
}
md <- md[!is.na(latitude) & !is.na(longitude)]
if (is_dry && nrow(md) > subset_n) md <- md[round(seq(1, .N, length.out = subset_n))]  # spread for variety
cat("Gages to process:", nrow(md), "\n")

# ---- load HydroATLAS gpkg once (peak memory: ~1.5-2 GB) ----
cat("Loading gpkg (167k polygons; this is the slow step)...\n")
basinAt <- sf::st_read(HYDROATLAS_GPKG, quiet = TRUE)

cat("Mapping gages to outlet basins (single st_join)...\n")
gage2outlet <- map_gages_to_outlets(md, basinAt)

# strip geometry -> attribute data.table, then free the heavy sf object
HB_strip <- basinAt; sf::st_geometry(HB_strip) <- NULL
HB_dt <- data.table::as.data.table(HB_strip)
rm(HB_strip, basinAt); gc()

rule_table <- classify_hydroatlas_attributes(names(HB_dt))
cat("Attribute rules:",
    paste(sprintf("%s=%d", names(table(rule_table$rule)), as.integer(table(rule_table$rule))),
          collapse = "  "), "\n")

outlet_ids <- gage2outlet[grepl("^[0-9]+$", Downstream_HB_ID), unique(Downstream_HB_ID)]
cat(sprintf("Resolved outlets: %d unique (of %d gages; %d unresolved)\n",
            length(outlet_ids), nrow(gage2outlet),
            gage2outlet[Downstream_HB_ID == "NO_BASIN_FOUND", .N]))

# ---- upstream member sets (cache first, BFS the rest) ----
cat("Building upstream member sets...\n")
cached <- if (file.exists("upstream_hydrobasins.rds")) readRDS("upstream_hydrobasins.rds") else list()
upstream_members <- build_upstream_members(outlet_ids, HB_dt, cached)
member_union <- unique(unlist(upstream_members, use.names = FALSE))
cat(sprintf("  outlets: %d | unique member basins: %d | from cache: %d\n",
            length(upstream_members), length(member_union), sum(outlet_ids %in% names(cached))))

# ---- attribute table (subset to members; cache parquet) ----
attr_cache <- file.path(out_dir, sprintf("hydroatlas_member_attrs_%s.parquet", date_tag))
attr_dt <- load_hydroatlas_attributes(HB_dt, member_union, rule_table$attribute, cache_path = attr_cache)
rm(HB_dt); gc()
cat(sprintf("Attribute table: %d basins x %d cols  (cached: %s)\n",
            nrow(attr_dt), ncol(attr_dt), attr_cache))

# ---- aggregate (per unique outlet) + assemble (per gage) ----
cat("Aggregating per outlet...\n")
outlet_agg <- aggregate_hydroatlas_per_outlet(attr_dt, upstream_members, rule_table, outlet_ids)
gage_dt <- assemble_gage_hydroatlas(md, gage2outlet, outlet_agg)

# ---- validate + write ----
validate_hydroatlas_metadata(gage_dt, rule_table, attr_dt, upstream_members)
write_hydroatlas_metadata(gage_dt, rule_table, out_dir, date_tag)

cat(sprintf("\nDONE in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

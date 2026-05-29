#!/usr/bin/env Rscript
# =============================================================================
# assemble_explorer.R  (PHASE 2: points + self-contained HTML)
# Reads the HydroATLAS metadata + streamflow signatures + watershed borders,
# builds the per-gage point data (color variables from both datasets), and
# injects everything into explorer_template.html -> data_out/streamflow_explorer.html
# Run AFTER build_explorer.R has produced data_out/watershed_borders.geojson.
# =============================================================================
suppressMessages({ library(data.table); library(jsonlite) })
source("R/aggregate_hydroatlas_metadata.R")   # canonical_gage_id()
t0 <- Sys.time()
OUT_HTML <- "data_out/streamflow_explorer.html"
TPL <- "explorer_template.html"
BORDERS_GEOJSON <- "data_out/watershed_borders.geojson"

ha <- fread("data_out/watershed_hydroatlas_metadata_29may2026.csv",
            colClasses = c(gage_id = "character", Downstream_HB_ID = "character"))
sig_f <- tail(sort(list.files("data_out", "^streamflow_signatures_full.*\\.csv$", full.names = TRUE)), 1)
sig <- fread(sig_f, colClasses = c(gage_id = "character"))
cat(sprintf("signatures: %s (%d x %d)\n", basename(sig_f), nrow(sig), ncol(sig)))

ha[, ck := canonical_gage_id(gage_id)]
sig[, ck := canonical_gage_id(gage_id)]

pts <- ha[!is.na(longitude) & !is.na(latitude),
          .(gage_id, lon = round(longitude, 5), lat = round(latitude, 5),
            hb = Downstream_HB_ID, gage_type, ck)]

# ---- variable selection ----
ha_excl <- c("gage_id", "latitude", "longitude", "gage_type", "Downstream_HB_ID", "ck")
ha_vars <- setdiff(names(ha), ha_excl)
ha_vars <- ha_vars[vapply(ha[, ..ha_vars], is.numeric, logical(1))]

sig_num  <- names(sig)[vapply(sig, is.numeric, logical(1))]
sig_vars <- grep("(_mean|_senn_slp|_mk_pval)$", sig_num, value = TRUE)
sig_vars <- c(sig_vars, intersect(c("elasticity_static",
                                    "recession_alpha_point_cloud_linear_reservoir"), sig_num))

# merge HA vars onto pts (all present); merge signature vars by canonical id
pts <- merge(pts, ha[, c("ck", ha_vars), with = FALSE], by = "ck", all.x = TRUE, sort = FALSE)
sig_sub <- unique(sig[, c("ck", sig_vars), with = FALSE], by = "ck")
pts <- merge(pts, sig_sub, by = "ck", all.x = TRUE, sort = FALSE)

is_cat <- function(v) grepl("_cl_smj$|_id_smj$", v)
mk_meta <- function(v, ds) {
  x <- pts[[v]]
  if (is_cat(v)) {
    list(key = v, label = v, ds = ds, type = "categorical", cats = sort(unique(x[!is.na(x)])))
  } else {
    q <- as.numeric(quantile(x, c(0.02, 0.98), na.rm = TRUE))
    if (!all(is.finite(q)) || q[1] == q[2]) q <- range(x, na.rm = TRUE)
    if (!all(is.finite(q))) q <- c(0, 1)
    list(key = v, label = v, ds = ds, type = "continuous", domain = signif(q, 4))
  }
}
meta <- c(
  list(list(key = "gage_type", label = "gage_type", ds = "HydroATLAS",
            type = "categorical", cats = c("USGS", "Canada"))),
  lapply(ha_vars, mk_meta, ds = "HydroATLAS"),
  lapply(sig_vars, mk_meta, ds = "Signatures")
)

# value columns (rounded), column-oriented
valcols <- list(gage_type = pts$gage_type)
for (v in c(ha_vars, sig_vars)) {
  x <- pts[[v]]
  valcols[[v]] <- if (is.numeric(x) && !is_cat(v)) signif(x, 5) else x
}

POINTS <- list(n = nrow(pts), gage_id = pts$gage_id, lon = pts$lon, lat = pts$lat,
               hb = pts$hb, v = valcols)
points_json <- toJSON(POINTS, na = "null", auto_unbox = TRUE, digits = 6)
meta_json   <- toJSON(meta, na = "null", auto_unbox = TRUE)
cat(sprintf("points: %d gages, %d color vars | POINTS JSON %.1f MB | meta %.0f KB  (%.1f min)\n",
            nrow(pts), length(meta), nchar(points_json) / 1e6, nchar(meta_json) / 1e3,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

if (!file.exists(BORDERS_GEOJSON)) {
  cat("Borders not ready (", BORDERS_GEOJSON, "). Re-run after build_explorer.R finishes.\n", sep = "")
  quit(save = "no")
}
borders_json <- paste(readLines(BORDERS_GEOJSON, warn = FALSE), collapse = "")
cat(sprintf("borders GeoJSON: %.1f MB\n", nchar(borders_json) / 1e6))

tpl <- paste(readLines(TPL, warn = FALSE), collapse = "\n")
inject <- function(tpl, ph, val) { p <- strsplit(tpl, ph, fixed = TRUE)[[1]]
  paste0(p[1], val, if (length(p) > 1) paste(p[-1], collapse = ph) else "") }
tpl <- inject(tpl, "__POINTS__", points_json)
tpl <- inject(tpl, "__VARS__", meta_json)
tpl <- inject(tpl, "__BORDERS__", borders_json)
writeLines(tpl, OUT_HTML, useBytes = TRUE)
cat(sprintf("wrote %s : %.1f MB  (%.1f min)\n", OUT_HTML, file.size(OUT_HTML) / 1e6,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

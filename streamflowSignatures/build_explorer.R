#!/usr/bin/env Rscript
# =============================================================================
# build_explorer.R  (PHASE 1: watershed border outlines)
# Dissolved, simplified, border-only watershed outlines (one per unique gage
# outlet) -> compact GeoJSON for the static HTML explorer.
#   Rscript build_explorer.R --n 60     # test on 60 outlets
#   Rscript build_explorer.R            # all outlets (run in background; ~minutes)
# Output: data_out/watershed_borders.geojson  (keyed by `outlet` = Downstream_HB_ID)
# =============================================================================
suppressMessages({ library(sf); library(data.table) })
options(scipen = 999); sf::sf_use_s2(FALSE)
source("R/aggregate_hydroatlas_metadata.R")   # build_upstream_members(), as_hybas_chr()
t0 <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(f, d = NA) { i <- which(args == f); if (length(i) && i[1] < length(args)) args[i[1]+1] else d }
n_lim <- suppressWarnings(as.integer(getarg("--n", NA)))
tol   <- as.numeric(getarg("--tol", "0.01"))             # ~1 km base simplification
maxv  <- suppressWarnings(as.integer(getarg("--maxv", "250")))  # per-outline vertex cap (big basins coarsened)
out_geojson <- file.path("data_out", "watershed_borders.geojson")

# Drop interior rings (sliver holes from the per-basin pre-simplify) and tiny detached
# sliver parts, so a dissolved watershed renders as a clean border with no interior lines.
.ring_area <- function(m) { x <- m[, 1]; y <- m[, 2]; n <- length(x); abs(sum(x[-n]*y[-1] - x[-1]*y[-n]))/2 }
remove_holes <- function(g) {
  geom <- sf::st_geometry(g)[[1]]
  if (inherits(geom, "GEOMETRYCOLLECTION")) {                # union sometimes yields polygon + line debris
    polys <- Filter(function(x) inherits(x, c("POLYGON", "MULTIPOLYGON")), geom)
    if (!length(polys)) return(g)
    geom <- if (length(polys) == 1) polys[[1]]
            else sf::st_geometry(sf::st_union(sf::st_sfc(polys, crs = 4326)))[[1]]
  }
  if (inherits(geom, "POLYGON")) return(sf::st_sfc(sf::st_polygon(geom[1]), crs = 4326))
  if (inherits(geom, "MULTIPOLYGON")) {
    ext <- lapply(geom, function(p) p[1])                  # exterior ring of each part only
    if (length(ext) > 1) { a <- vapply(ext, function(r) .ring_area(r[[1]]), numeric(1))
                           ext <- ext[a >= max(a) * 0.01] }    # drop parts < 1% of the largest
    return(sf::st_sfc(sf::st_multipolygon(ext), crs = 4326))
  }
  g
}

ha <- fread("data_out/watershed_hydroatlas_metadata_29may2026.csv",
            colClasses = c(gage_id = "character", Downstream_HB_ID = "character"))
outlets <- unique(ha[grepl("^[0-9]+$", Downstream_HB_ID), Downstream_HB_ID])
if (!is.na(n_lim)) outlets <- outlets[round(seq(1, length(outlets), length.out = min(n_lim, length(outlets))))]
cat(sprintf("outlets to build: %d (tol=%.3f deg)\n", length(outlets), tol))

b <- sf::st_read("D:/geospatial_derivedData/basinAt_NorAm_polys.gpkg", quiet = TRUE)
HB_dt <- data.table(HYBAS_ID = b$HYBAS_ID, NEXT_DOWN = b$NEXT_DOWN)
cached <- if (file.exists("upstream_hydrobasins.rds")) readRDS("upstream_hydrobasins.rds") else list()
members_list <- build_upstream_members(outlets, HB_dt, cached)   # cache + NEXT_DOWN BFS for the rest
cat(sprintf("member sets ready for %d outlets (%.1f min)\n", length(members_list),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

member_union <- unique(unlist(members_list, use.names = FALSE)); member_union <- member_union[!is.na(member_union)]
hb_all <- as.character(b$HYBAS_ID)
sel <- which(hb_all %in% member_union)
# Use RAW basin geometry (NOT pre-simplified): adjacent basins keep identical shared edges, so the
# union dissolves cleanly -> no sliver holes/spikes. Simplify happens AFTER the union, on the clean outline.
bsub_geom <- sf::st_geometry(b)[sel]
bmap <- data.table(hybas_chr = hb_all[sel], gidx = seq_along(sel)); setkey(bmap, hybas_chr)
cat(sprintf("using %d raw needed basins (%.1f min)\n", length(sel),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

geoms <- vector("list", length(outlets)); ids <- character(length(outlets)); miss <- 0L
for (i in seq_along(outlets)) {
  mem <- members_list[[outlets[i]]]
  gid <- if (!is.null(mem)) bmap[.(mem), gidx, nomatch = NULL] else integer(0)
  if (!length(gid)) { miss <- miss + 1L; next }
  u <- tryCatch(suppressWarnings(sf::st_union(bsub_geom[gid])), error = function(e) NULL)
  if (is.null(u) || length(u) == 0 || sf::st_is_empty(u)[1]) { miss <- miss + 1L; next }
  u <- tryCatch(remove_holes(u), error = function(e) u)    # drop interior sliver holes/parts -> clean fill
  # simplify the clean (hole-free) outline; preserveTopology=TRUE avoids self-intersections
  g <- suppressWarnings(sf::st_simplify(u, dTolerance = tol, preserveTopology = TRUE))
  np <- nrow(sf::st_coordinates(g)); t2 <- tol; it2 <- 0L
  while (np > maxv && it2 < 6L) {
    t2 <- t2 * 1.8
    g <- suppressWarnings(sf::st_simplify(g, dTolerance = t2, preserveTopology = TRUE))
    np <- nrow(sf::st_coordinates(g)); it2 <- it2 + 1L
  }
  if (length(g) == 0 || sf::st_is_empty(g)[1]) { miss <- miss + 1L; next }
  geoms[[i]] <- sf::st_geometry(g)[[1]]; ids[i] <- outlets[i]
  if (i %% 1000 == 0) cat(sprintf("  %d/%d  (%.1f min)\n", i, length(outlets),
                                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
keep <- which(nzchar(ids))
borders <- sf::st_sf(outlet = ids[keep], geometry = sf::st_sfc(geoms[keep], crs = 4326))
cat(sprintf("built %d outlines; %d missing\n", nrow(borders), miss))

if (!dir.exists("data_out")) dir.create("data_out", recursive = TRUE)
suppressWarnings(sf::st_write(borders, out_geojson, driver = "GeoJSON", delete_dsn = TRUE,
                              quiet = TRUE, layer_options = c("COORDINATE_PRECISION=4", "RFC7946=YES")))
# Final hole-strip on the ROUNDED geometries: the topology-preserving simplifier can pinch a
# concavity into a hole, which errors remove_holes on full-precision in-memory geometry but
# strips cleanly once GDAL rounds coords to 4 dp. So re-read, strip, and re-write.
ww <- sf::st_read(out_geojson, quiet = TRUE)
ggv <- suppressWarnings(sf::st_make_valid(sf::st_geometry(ww)))   # fix simplifier self-intersections
ww_clean <- lapply(ggv, function(gi)
  sf::st_geometry(tryCatch(remove_holes(sf::st_sfc(gi, crs = 4326)),
                           error = function(e) sf::st_sfc(gi, crs = 4326)))[[1]])
sf::st_geometry(ww) <- sf::st_sfc(ww_clean, crs = 4326)
suppressWarnings(sf::st_write(ww, out_geojson, driver = "GeoJSON", delete_dsn = TRUE,
                              quiet = TRUE, layer_options = c("COORDINATE_PRECISION=4", "RFC7946=YES")))
cat(sprintf("wrote %s : %.2f MB  (total %.1f min)\n", out_geojson, file.size(out_geojson)/1e6,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

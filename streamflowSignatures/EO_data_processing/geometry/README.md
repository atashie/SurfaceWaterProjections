# Watershed geometry build (EO universe)

Builds one watershed polygon per QA-passing gage (universe = **8,018** = metadata
`processing_status==success` 8,014 + 4 signature-bearing `processing` USGS gages),
**official-primary / HydroBasins-fallback**, keyed by zero-padded `gage_id` (+ `canon_id`).
**Delivered geometry = 7,964** after excluding **54 basins with drainage >100,000 km²**
(53 huge Canadian/HB rivers + `05KH009`, pruned 27 Jun — see QA).

**Deliverable** (on `s3://climate-ai-data-science-shiny-app-data/streamflow/`):
`watershed_polygons_26jun2026.{gpkg,parquet}` (EPSG:4326) + `_qa.csv`.

## Sources
| Source | Gages | Where | Key | CRS |
|---|---|---|---|---|
| GAGES-II (US official) | 6,164 | ScienceBase DOI 10.5066/P96CPHOT `boundaries_shapefiles_by_aggeco.zip` | `GAGE_ID` | NAD83 Albers |
| ECCC MDA_ADP (CAN official) | 1,823 | `collaboration.cmc.ec.gc.ca/.../HydrometricNetworkBasinPolygons/gpkg/MDA_ADP_{01..11}.gpkg` layer `DrainageBasin_BassinDeDrainage` | `StationNum` | ESRI:102001 |
| HydroBasins (fallback) | 29 | `basinAt_NorAm_polys.gpkg` (S3) lev12, BFS over `NEXT_DOWN` from `Downstream_HB_ID` or lat/lon→outlet | `HYBAS_ID` | EPSG:4326 |

## Pipeline (run order)
1. `coverage.py` — intersect passing gages vs each official source's station list.
2. `geometry_build_us.py` / `geometry_build_ca.py` — per-source official polygons → 4326 (`official/us_gagesII_4326.gpkg`, `ca_mda_4326.gpkg`).
3. `geometry_build_hb_fallback.py` — delineate the 31 residual from `basinAt` (`hb_fallback_4326.gpkg`).
4. `reconcile.py` — leading-zero reconciliation → canonical universe (8,018); canonical key rule.
5. `build_v3.py` — **canonical builder**: merge sources, dual `gage_id`/`canon_id`, metadata attrs,
   **adaptive simplification** (~200 m; keep full-res where it would shift area >2% — 141 basins),
   area from delivered geometry (EPSG:6933 equal-area), flags (`geom_simplified`, `area_flag`,
   `low_confidence`), assertions (no canon/gage_id collisions). Writes the deliverable.

## Output columns
`gage_id` (zero-padded; joins signatures), `canon_id` (zero-stripped; joins metadata),
`watershed_geom_source` {gagesii,wsc_eccc,hydrobasins}, `geom_area_km2`, `basin_area`,
`area_rel_diff`, `area_flag` (>50%), `geom_simplified`, `low_confidence`, `latitude`,
`longitude`, `gage_type`, `geometry`.

## QA (v3, updated 27 Jun — Codex-reviewed)
**7,964 delivered** of the 8,018 universe — 54 basins with drainage >100,000 km² excluded
(53 huge rivers + `05KH009`, an HB fallback that mis-delineated a 328,000 km² basin as
200 km² and slipped the geom-area filter; pruned 27 Jun). 0 dups/invalid/empty.
Area vs reported `basin_area`: median ~0.001, 97.8% <10% — **US is a consistency check
only** (`basin_area` is GAGES-II-derived, so not independent validation). `low_confidence`
= 29 HB fallbacks + area outliers, incl. 3 official CAN polygons (`04HA001`/`06FD001`/
`06FB001`) whose reported `basin_area` ≫ official polygon area (suspect metadata; kept).
Known limitation: the 15 lat/lon-derived HB fallbacks are approximate (gage may sit
mid-basin) — all `low_confidence`.

> Notes: paths in these scripts were session-scratch absolute paths — adjust before re-running.
> HydroBasins point-derived fallbacks (15 of the 31, no `Downstream_HB_ID`) are approximate
> (gage may sit mid-basin) — all carry `low_confidence=True`. To promote to a proper package,
> consolidate into `eo_processing/geometry.py` with the data-download steps inlined.

# EO_data_processing — per-watershed MODIS (LAI, LULC), Annual NLCD and Daymet products

Metadata/ingestion products (Python), NOT cross-language signatures — nothing here is
ported to Julia or R. The MODIS and NLCD pipelines are COMPLETE and their outputs are staged
as HydroShare Resource 5; do not rebuild unless asked. Design, pipeline notes and the dated
build history: `README.md` (MODIS LAI/LULC + the 7,964-watershed geometry layer),
`README_NLCD.md` (CONUS Annual NLCD 1985–2025) and `README_DAYMET.md` (the reprocessed
Daymet climate input, 8,017 basins, 1980–2025; tools in `daymet/`, run plan in
`docs/plans/2026-09-29-daymet-reprocessing-action-plan.md`).

Rules that still bite:
- **Ids.** Outputs carry `gage_id` (zero-padded) and `canon_id` (zero-stripped). The
  build's `zfill(8)` fallback left 44 nine-to-ten-digit USGS ids UN-padded in the
  delivered EO tables (9 in the boundary layer). User decision 2026-09-04: the files
  stay as delivered; every join uses the zero-stripped form on BOTH sides and never
  re-pads. Any future build must carry the agency id from the streamflow parquet rather
  than re-derive it.
- **S3 is gone** (access lost 2026-08-24). Every `s3://` path in the READMEs is a
  historical delivery record; the project Google Drive folder is the backup, and the
  geometry layer was rebuilt locally 2026-08-25 (every June target matched exactly).
- **Universe.** 8,014 `processing_status == success` gages (+ 4 → 8,018 canonical
  union); delivered geometry = 7,964 after excluding 54 basins > 100,000 km²; NLCD =
  6,119 CONUS gages (45 Alaska gages excluded, never published as zeros).
- **Extraction discipline.** Reproject polygons to each raster's native CRS read from
  `src.crs` (never hardcoded, never EPSG:5070 for NLCD); never resample categorical
  rasters; mask raw fill BEFORE scaling; reindex class-% columns to the committed legend
  manifests (`eo_processing/lulc_legends.csv`, `nlcd_legends.csv`); unknown codes hard-fail.
- **Run bookkeeping.** Per-(tile, period) checkpoints with skip-if-exists; delete empty
  checkpoint parquets before resuming; pop static `AWS_*` env vars in worker processes.
- **Daymet.** Weights are coverage × TRUE cell area over the UNSIMPLIFIED polygons (that
  reproduces the co-authors' series; coverage-only or the simplified layer do not); fill is
  excluded per day, never NaN-propagated. The Earthdata token lives in
  `~/.config/earthdata/edl_token` (mode 600), never in the repo, docs or memory. The user
  writes it from their own terminal, since a token pasted into a session persists in the
  Claude transcripts. The stream hands it to curl on stdin. `--workers 6` (the default)
  suits the 16 GB laptop with the 8,017-basin weights. Run
  `daymet/selftest_daymet_tools.py` after changing the tools.

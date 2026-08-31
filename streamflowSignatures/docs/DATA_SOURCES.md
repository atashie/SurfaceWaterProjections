# Data Sources

All external data sources used across the project — signature extraction, watershed
metadata, and Earth Observation products. Each row lists the data class, what it is
used for in this project, the data provider, and where the original data is accessible.

> Compiled 2026-07-14 from the workflow docs (`CLAUDE.md`, `docs/DEVELOPMENT.md`,
> `EO_data_processing/README.md`, `EO_data_processing/geometry/README.md`), `config.R`,
> and the ingestion / EO pipeline code.

## Primary data sources

| # | Data class | What it's used for | Provider | Original data access |
|---|---|---|---|---|
| 1 | **US daily streamflow** (param 00060, cfs; incl. site metadata / drainage area) | Core signal for all 100+ signatures; ingested via `run_ingest_usgs_hydat.R` (`dataRetrieval::readNWISdv`) → parquet | U.S. Geological Survey (NWIS) | https://waterdata.usgs.gov / web services at https://waterservices.usgs.gov (R package `dataRetrieval`) |
| 2 | **Canadian daily streamflow** (m³/s) + station regulation metadata (RHBN, REGULATED via `hy_stn_regulation`) | Same signature pipeline; regulation flags feed `human_interference_class` (exported to `metadata/canadian_hydat_interference.csv`). Caveat: HYDAT publishes no drainage area for some stations (mostly canals/dam outflows/channel splits) — those gages stay in raw m³/s with `area_normalized = FALSE` (see DEVELOPMENT.md → Canadian HYDAT) | Water Survey of Canada / Environment and Climate Change Canada (HYDAT database) | https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/ (HYDAT SQLite archive; R package `tidyhydat`) |
| 3 | **Caravan** — bundled daily streamflow + climate forcing, NetCDF (CAMELS, HYSETS, etc.) | Alternative processing pipeline (`run_caravan_processing.R`); shorter records but has co-located climate | Kratzert et al. (community dataset, Google Research) | https://github.com/kratzert/Caravan → Zenodo archive (Kratzert et al. 2023, *Scientific Data*) |
| 4 | **Daymet daily climate** — PPT, tmin/tmax, SWE, vapor pressure, solar radiation, 1980–2023, 1 km, basin-aggregated to 6,087 sites | Climate-dependent signatures: runoff ratios, elasticity, Q-P seasonality, storage, and all 14 snow metrics. Source of record is the **44 annual CSVs** (`daymet_1980_2023/`, 10.08 GiB, columns `site_id, month, year, prcp, tmin, tmax, swe, vp, srad` — **no day column**); the pipeline consumes a parquet built from them (97,757,220 rows, 3.76 GiB), joined at runtime. **Two traps**: Daymet publishes a **365-day calendar** — leap years keep Feb 29 and DROP Dec 31 — so the series has a one-day hole each leap year; and the canonical `daymet_1980_2023.parquet` on the drive is TRUNCATED, so use `daymet_1980_2023_rebuilt_10aug2026.parquet` (rebuild tool: `docs/benchmarks/convert_daymet_csvs_to_parquet.py`). See DEVELOPMENT.md → Active Parquet Files | ORNL DAAC (NASA); basin aggregation by co-authors via `gdptools` | https://daymet.ornl.gov (Daymet V4, ORNL DAAC) |
| 5 | **US watershed boundaries** — GAGES-II curated basin polygons (`boundaries_shapefiles_by_aggeco.zip`) | Primary US geometry for EO zonal statistics (6,164 basins; 100% coverage, so the NLDI/NHDPlus/StreamStats backstops were never needed) | USGS | ScienceBase DOI [10.5066/P96CPHOT](https://doi.org/10.5066/P96CPHOT) |
| 6 | **Canadian watershed boundaries** — National Hydrometric Network Basin Polygons (MDA_ADP GeoPackages) | Primary Canadian geometry for EO zonal statistics (1,823 basins) | ECCC / Water Survey of Canada | https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/HydrometricNetworkBasinPolygons/gpkg/ (pin June-2024 version) |
| 7 | **HydroBASINS level-12 polygons + topology** (`basinAt_NorAm_polys.gpkg`, 167k basins, `NEXT_DOWN`) | Fallback watershed delineation (BFS upstream), explorer watershed borders, and the spatial units underlying the HydroATLAS aggregation | HydroSHEDS project (WWF; Lehner & Grill 2013) | https://www.hydrosheds.org/products/hydrobasins |
| 8 | **HydroATLAS / BasinATLAS v10** — 281 hydro-geophysical attributes (climate, terrain, soils, land cover, anthropogenic) | Per-gage upstream-watershed metadata product (`run_hydroatlas_metadata.R`, ~211 output columns) | HydroSHEDS project (Linke et al. 2019) | https://www.hydrosheds.org/hydroatlas (BasinATLAS v10) |
| 9 | **MODIS LAI** — MCD15A3H v061, 4-day, 500 m, 2002–2024 | Per-watershed monthly LAI panel (mean + spatial heterogeneity + QA), 7,964 basins × 270 months | NASA LP DAAC (Terra+Aqua MODIS) | https://lpdaac.usgs.gov/products/mcd15a3hv061/ via Earthdata; operationally read from the ClimateAI catalog mirror (`s3://climateai-geospatial-data/data/MCD15A3H/`) with LP DAAC backfill |
| 10 | **MODIS LULC** — MCD12Q1 v061, annual, 500 m, 2001–2024, all 8 classification bands | Per-watershed annual land-cover composition (% per class, 102 class columns) | NASA LP DAAC | https://lpdaac.usgs.gov/products/mcd12q1v061/ via Earthdata CMR search + authenticated HTTPS download |
| 11 | **GAGES-II attribute tables** (`conterm_hydromod_dams.txt`, `conterm_bas_classif.txt`, …) | US human-interference metadata: dam counts/storage, impervious/developed %, withdrawals, disturbance index, Ref/Non-ref classification → `human_interference_class` | USGS (Falcone 2011, GAGES-II) | Same GAGES-II ScienceBase release as the boundaries — DOI [10.5066/P96CPHOT](https://doi.org/10.5066/P96CPHOT) (local copy at `D:/gagesMetadata`) |

### Governance input (methodology, not data)

The **collaborative signature guidelines Google Doc** — domain-expert plain-English
signature definitions and QA/QC requirements, auto-synced each session to
`docs/SIGNATURE_GUIDELINES.md`. It is the source of truth for *how* the data above
gets processed.
URL: https://docs.google.com/document/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub

## Deliberately excluded

Boundary cases that are not counted as project data sources, with reasons:

- **Planned but not yet used**: ERA5/PRISM climate fetching and MODIS ET are on the
  roadmap (`CHANGELOG.md` → `[Unreleased]`) but nothing consumes them today.
- **Not origins**: the ClimateAI S3 buckets (`climate-ai-data-science-shiny-app-data`,
  `climateai-geospatial-data`) hosted derived project outputs and a MODIS mirror —
  access routes, not providers. **Access to these buckets was lost on 2026-08-24**;
  most delivered outputs were backed up to the project Google Drive first (link in
  CHANGELOG → August 2026).
- **Viz-only third parties**: Leaflet/Plotly CDN and CARTO basemap tiles used by the
  HTML explorers and dashboards; they render maps but contribute no analysis data.
- **US geometry backstops never exercised**: USGS NLDI
  (`https://api.water.usgs.gov/nldi/`), Streamgage NHDPlus V1 Basins
  (DOI [10.5066/P9EZY2M6](https://doi.org/10.5066/P9EZY2M6)), and StreamStats were
  designed as fallbacks for US watershed boundaries, but GAGES-II covered 6,164/6,164
  gages so they were never used.

## Naming note: HydroBASINS vs HydroATLAS

These are distinct HydroSHEDS products and easy to conflate: **HydroBASINS** supplies
the level-12 polygon units and `NEXT_DOWN` topology (rows 7); **HydroATLAS/BasinATLAS**
supplies the 281 attributes layered onto those units (row 8). This project uses both.

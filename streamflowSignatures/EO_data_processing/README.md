# EO Data Processing

Ingestion and processing of global MODIS Earth Observation (EO) products — Land
Use / Land Cover (**LULC**, MCD12Q1) and Leaf Area Index (**LAI**, MCD15A3H) —
aggregated per watershed into per-product temporal records that sit alongside the
streamflow signatures and join by `gage_id`.

This is a **metadata / ingestion product** (like the HydroATLAS watershed metadata),
**not** a cross-language signature — so it is implemented once and is **not** ported
to Julia/R. See "Implementation language" below for why this one is **Python** rather
than R.

> **Sibling product:** a CONUS-only, higher-resolution **Annual NLCD** land-cover + fractional
> impervious product (USGS/MRLC, 30 m, 1985–2025) is built by the same routine, adapted — see
> [`README_NLCD.md`](README_NLCD.md) (`eo_processing/nlcd_pipeline.py`, `nlcd_finalize.py`,
> `viz/build_nlcd_explorer.py`). It complements this continental MODIS LULC (500 m, N. America).

> **⚠ ACCESS CHANGE (24 Aug 2026)**: the project **no longer has access to the S3 buckets**
> referenced throughout this README (`climate-ai-data-science-shiny-app-data`, the
> ClimateAI catalog). Most delivered artifacts were backed up beforehand to the project
> Google Drive folder
> (https://drive.google.com/drive/u/0/folders/1DVuq4nC5j_Y01sBaDj9cwjbv7S8sndjj).
> Treat every S3 URL below as a historical delivery record, and verify the Drive
> inventory before relying on any file. See CHANGELOG → August 2026.

> **Status (30 Jun 2026)**: ALL THREE DATA LAYERS DONE + on S3 — geometry (7,964 watersheds),
> LAI (monthly, 270-month panel; **`good_coverage_frac` QA column added → `..._30jun2026.parquet`**),
> and LULC (annual, 24-yr, Codex-GO). **BOTH QA/QC explorers built + on S3**:
> `watershed_modis_lai_explorer.html` (builder `viz/build_lai_explorer.py`) and
> `watershed_modis_lulc_explorer.html` (builder `viz/build_lulc_explorer.py`). The LULC explorer
> maps 27 variables (10 IGBP-derived summary + all 17 individual IGBP classes; non-overlapping —
> static `pct_*` roll-ups dropped) and draws per-gage stacked-area class composition 2001–2024 with
> a band switcher across all 8 MODIS schemes. **EO pipeline + viz COMPLETE.**
>
> **Status (29 Jun 2026)**: ALL THREE LAYERS DONE + on S3 — geometry (7,964 watersheds),
> LAI (monthly, 270-month panel), and LULC (annual, 24-yr, Codex-GO). Pipeline COMPLETE.
>
> **Status (27 Jun 2026)**: Geometry layer DONE + on S3 (7,964 watersheds, §4/§11).
> **LAI pipeline DONE + backfilled** — 270-month full panel (7,964×270 = 2,150,280 rows;
> far-N + 2024-11 backfilled from LP DAAC; only 17 urban basins NA).
> **Delivered to S3** (`watershed_modis_lai_monthly_29jun2026.parquet` + `_na_basins_` +
> `_dictionary_`); `partial_month` flag added. **LAI COMPLETE.**
> **LULC pipeline DONE + Codex-GO + DELIVERED to S3 (29 Jun)** — full 24-yr run (2001–2024),
> finalized + reviewed (gpt-5.5, GO; 2 latent code hazards hardened, output byte-identical).
> **191,136 rows = 7,964 gages × 24 yr**, all 8 bands sum to exactly 100 across all years, 0 dup
> keys. On S3 `s3://climate-ai-data-science-shiny-app-data/streamflow/`:
> `watershed_modis_lulc_annual_29jun2026.{parquet,csv}` + `_dictionary.csv` + `_granule_manifest_29jun2026.parquet`.
> **Both LAI + LULC products COMPLETE + on S3.**
> **Resume after a crash → §11 task 5** (re-run the same command; checkpoints persist on
> `/home`). Durable cross-session notes: memory `eo-data-access-paths`.

---

## 1. Goal & deliverables

**Universe = most-inclusive QA-passing set = 8,018 gages** (6,164 US + 1,854 Canadian) —
the canonical union of metadata `processing_status == success` (8,014) and all signature
CSVs (adds 4 `processing`-status USGS gages that carry signatures). From
`combined_watershed_metadata_09feb2026.csv` (16,994 rows); `no_data` 8,787 etc. excluded.

> **⚠ CORRECTION (2026-09-04)**: the `zfill(8)` rule below cannot restore the leading zero
> of USGS site numbers that are 9–10 digits long (`021650905`, `0212433550`): after
> stripping they are 8–9 digits and `zfill(8)` leaves them unchanged. The delivered EO
> tables therefore store **44** such ids un-padded in `gage_id` (and the rebuilt boundary
> layer 9). **Decision (user, 2026-09-04): the delivered files stay byte-identical; users
> join on the zero-stripped form (`canon_id`, or strip leading zeros on both sides) — never
> re-pad.** Documented in the root README → "Reading and Joining the Data" and in the
> HydroShare resource READMEs. Any future EO build should carry the agency id from the
> streamflow parquet instead of re-deriving it.

**Canonical `gage_id` key (IMPORTANT):** the metadata + the watershed-geometry source
files store gage_ids **zero-stripped** (`1011000`); the signature CSVs store them
**zero-padded** (`01011000`). Reconciled rule: numeric → strip leading zeros for the
canonical match, **store the zero-padded form** (US `zfill(8)`, or the authoritative
padded id from the signature files; 15-digit lat-long ids unchanged); Canadian HYDAT →
uppercase as-is. Outputs carry **`gage_id`** (zero-padded, joins to the signatures) and
**`canon_id`** (zero-stripped, joins to the metadata). Two flat tables (CSV **and**
parquet, mirroring the signature / HydroATLAS outputs), each with a data dictionary:

### 1a. LULC table — annual (`watershed_modis_lulc_{date}.{csv,parquet}`)
One row per **(gage_id, year)**, year ∈ 2001–2024.

| Column group | Columns |
|---|---|
| Keys | `gage_id`, `year` |
| Class coverage | `% coverage per class, for every classification band` — `{band}_c{code}_pct` (e.g. `lc_type1_c12_pct` = IGBP cropland %). All bands: **LC_Type1** (IGBP, 17), **LC_Type2** (UMD, 16), **LC_Type3** (LAI, 11), **LC_Type4** (BGC, 9), **LC_Type5** (PFT, 12), **LC_Prop1/2/3** (FAO-LCCS, sparse codes). ~200 class columns total. |
| QA | `n_modis_pixels` (effective coverage-weighted pixel count in basin), `valid_coverage_frac` (valid / total intersecting pixels, fill=255 excluded from the % denominator), `watershed_geom_source` |

Per row, the per-band class %s sum to ~100 over valid (non-fill) pixels.

### 1b. LAI table — monthly (`watershed_modis_lai_{date}.{csv,parquet}`)
One row per **(gage_id, year, month)**.

| Column group | Columns |
|---|---|
| Keys | `gage_id`, `year`, `month` |
| Monthly LAI stats | `lai_min, lai_q05, lai_q25, lai_q50, lai_q75, lai_q95, lai_max, lai_mean, lai_std` (LAI units; raw DN masked **before** the 0.1 scale — see §6; valid 0–10) |
| QA (implemented) | `n_eff_pixels` (usable-LAI pixel count), `n_basin_pixels` (land-pixel denominator), `valid_pixel_frac` = n_eff/n_basin (spatial % valid), `mean_obs_valid_frac` (fraction of land observations that were clear/usable — the discriminating quality metric), `frac_cloud` (cloud∪cirrus∪shadow), `frac_snow_ice`, `frac_aerosol`, `frac_backup_algorithm` (empirical vs main-RT retrievals), `n_composites`, `low_confidence`. Flags decoded from `FparExtra_QC` (B04) + `SCF_QC` (B03); all fractions ∈ [0,1] (Codex-verified). |

> **Note (exactextract semantics)**: `lai_mean / lai_std / lai_q*` are **coverage-weighted**
> across the basin's pixels; `lai_min / lai_max` are **unweighted extrema** over any
> intersecting non-fill pixel, so a single boundary sliver can drive them. They are
> emitted as documented unweighted extrema; downstream prefer q05/q95 for robust range.

**Monthly aggregation semantics — confirmed.** A **two-stage
spatial-distribution-of-monthly-mean**: (A) per pixel, mean of the month's QA-good
4-day composites → a monthly-mean LAI field; (B) the 9 stats describe the spatial
distribution across the watershed's pixels of that monthly-mean field. So `lai_mean` =
**basin mean monthly LAI**, and `lai_std`/quantiles = **spatial heterogeneity** — both
wanted (per user).

**Composite → month assignment**: MCD15A3H composites are 4-day periods stamped by
start date and **do not align to calendar months** (and cross year/leap boundaries).
Stage A weights each composite's contribution to a calendar month by its **day-overlap**
with that month (a composite spanning month-end contributes fractionally to both),
making monthly values reproducible and leap-year-safe. The assignment rule is config'd
and recorded.

---

## 2. Architecture & data flow

```
combined_watershed_metadata_09feb2026.csv   (gage list, gage_id key)
         │
         ▼
 Watershed geometry layer  (§4)   [universe: 8,014 QA-passing gages]
   PRIMARY:  GAGES-II (US) + WSC/ECCC (CAN) official basin polygons
   FALLBACK: HydroBasins delineation (basinAt + Downstream_HB_ID, Python)
             + watershed_geom_source provenance
         │
         ├───────────────────────────────┐
         ▼                                ▼
   LULC source (§5)                   LAI source (§6)
   MCD12Q1 v061  ── earthaccess ──    MCD15A3H ── GEOLATELY catalog ──
   LP DAAC (HTTPS, any region)        ClimateAI catalog (us-east-2)
   (HDF4, annual, ~10 GB)             (reprocessed per-band COGs, 4-day)
         │                                │
         └──────────────┬─────────────────┘
                        ▼
        Shared zonal-stats engine  (§7)   [runs from us-east-2]
        reproject polygons → each raster's native CRS (never resample categorical);
        tile-batched exactextract (coverage-weighted); mask raw fill before scaling
                        │
         ┌──────────────┴─────────────────┐
         ▼                                 ▼
   annual per-class %  (LULC)        monthly LAI stats (LAI)
         │                                 │
         ▼                                 ▼
   data_out/  CSV + parquet + data dictionary   (keyed gage_id, leading-zero-safe)
```

---

## 3. Region — runs from us-east-2 (no us-west-2 instance needed)

Resolved. The `wildfireModeling` pipeline runs this exact stack from **us-east-2**
(reference: `wildfireModeling/data-collection/v1-derived-features/MODIS_ACCESS_NOTES.md`):

- **LAI (catalog MCD15A3H)**: the ClimateAI catalog stores MCD15A3H as **reprocessed
  per-band COGs in its own bucket** (`s3://climateai-geospatial-data/data/MCD15A3H/`),
  read from **us-east-2** — the DuckDB catalog `SECRET` hardcodes `REGION 'us-east-2'`.
  The us-west-2 `403` lock applies to **NASA-origin** S3 (e.g. astraea-opendata
  MCD43A4), **not** to these ClimateAI copies. So LAI works directly from our instance.
- **LULC (MCD12Q1)**: not in the catalog → LP DAAC. Direct S3 would `403` from
  us-east-2, but **`earthaccess` HTTPS works from any region**; at ~10 GB the
  cross-region egress is cheap (chosen over provisioning us-west-2).

Net: **no us-west-2 provisioning** — both products reachable from the existing
us-east-2 instance. Operational caveat: **materialize AWS creds into env vars BEFORE
importing geolately** (DuckDB httpfs reads them at import), and refresh them
periodically — SageMaker creds rotate ~30 min (see §12).

---

## 4. Watershed geometry layer

- **Gage list / key**: `combined_watershed_metadata_09feb2026.csv` (S3), `gage_id`
  (zero-padded, leading-zero-safe canonical matching as in HydroATLAS). **Universe =
  the 8,014 `processing_status == success` gages** (6,160 USGS-like, 1,854 Canadian-like).
- **Source priority (per user): official polygons for ALL locations; HydroBasins only
  as a secondary fallback where official data is missing.**
  1. **PRIMARY — official basin polygons:**
     - **US (USGS)** — layered official sources (final ordering = §9 decision):
       - **GAGES-II** curated basin boundaries — ScienceBase DOI `10.5066/P96CPHOT`,
         `boundaries_shapefiles_by_aggeco.zip`; key **`GAGE_ID`**; CRS NAD83 Albers;
         9,322 basins; **fixed ~2011 set → misses gages activated after ~2011**.
       - **USGS NLDI** on-demand per site:
         `https://api.water.usgs.gov/nldi/linked-data/nwissite/USGS-<siteno>/basin`
         → GeoJSON polygon, **EPSG:4326**, ~0.3 s/call, broad NWIS coverage (active +
         historical), **CONUS only** (NHDPlusV2 — no AK/HI/PR), catchment-snapped
         (coarser at the outlet).
       - Backstops: **Streamgage NHDPlus V1 Basins** (DOI `10.5066/P9EZY2M6`, 19,031,
         CONUS) and **StreamStats** delineate (AK / true point delineation).
     - **Canada (ECCC/WSC)** — **National Hydrometric Network Basin Polygons** (MDA_ADP),
       `https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/HydrometricNetworkBasinPolygons/gpkg/`
       — 11 GeoPackages by Major Drainage Area (= first 2 digits of station no.); layer
       **`DrainageBasin_BassinDeDrainage`**; key **`StationNum`** (clean 7-char HYDAT);
       CRS **ESRI:102001 (Canada Albers)** — reproject; ~98% active / 93% discontinued;
       OGL-Canada; ~2.8 GB; **pin the June-2024 version**. (No live delineation API.)
  2. **FALLBACK — HydroBasins delineation** (only where no official polygon): BFS upstream
     over `NEXT_DOWN` in `basinAt_NorAm_polys.gpkg` (167,665 lev12 basins, valid GeoPackage)
     from each gage's **`Downstream_HB_ID`** (a column in the metadata), union + dissolve
     members. Pure-Python, in us-east-2, no R — reproduces the project's R explorer logic.
  3. Record **`watershed_geom_source ∈ {gagesii, wsc_eccc, hydrobasins}`** per gage;
     gages with neither official nor HB geometry are dropped + reported.
- **Coverage of the 8,014 passing gages — MEASURED**: official polygons cover **7,983
  (99.6%)** directly — **US GAGES-II 6,160/6,160 = 100%** (NLDI fill not needed), **Canada
  ECCC MDA_ADP 1,823/1,854 = 98.3%**. Only **31** Canadian gages fall to fallback: 16 have
  HydroBasins via `Downstream_HB_ID`; **15** have no geometry from any source → derive
  outlet from lat/lon (point-in-`basinAt`) then BFS, or drop + report. (HydroBasins
  fallback availability across all passing: 6,042/8,014 have a numeric `Downstream_HB_ID`.)
- **⚠️ `unified_watersheds_simplified.gpkg` (S3) is unusable** — despite the name/extension
  it is a corrupted 108k-line CSV, not a polygon layer (GDAL cannot open it). Do not use.
- **Simplification**: **topology-preserving** simplify to ~250 m tolerance (≤ half a
  500 m pixel) with a vertex cap (reuse the explorer's clean-dissolve approach). Cheap at
  500 m for large basins, but it can materially distort **narrow / small / coastal /
  island / headwater** polygons — so **validate simplified-vs-unsimplified zonal output
  on a stratified sample** (by area), apply simplification only above a size threshold,
  and **record the tolerance + a geometry version** in outputs.
- **Area QA**: `watershed_area_rel_diff` = |Σ coverage-weighted pixel-area − reported
  `basin_area`| / `basin_area` to flag mis-snapped / divergent geometries (as in
  HydroATLAS). MODIS sinusoidal is **equal-area on its sphere**, so coverage fractions in
  the raster CRS are the right basis for pixel area; document the area convention so QA
  doesn't misread datum/ellipsoid differences vs. the reported drainage area as geometry
  errors.

---

## 5. LULC source — MCD12Q1 v061 (via LP DAAC `earthaccess`)

- **Not in the ClimateAI catalog** (catalog LULC = IO_LULC / WORLDCOVER / USDA_CDL
  only) → fetch from **NASA LP DAAC Earthdata Cloud**: `short_name="MCD12Q1",
  version="061"` via `earthaccess` (Earthdata Login / `.netrc`). Direct S3 would `403`
  from us-east-2, but **earthaccess HTTPS works from any region**; ~10 GB → cheap egress.
- **Coverage**: annual, 500 m, sinusoidal, stamped Jan 1. LP DAAC lists the temporal
  extent as **2001–present**; the latest v061 year confirmed is **2024** (use 2001–latest,
  not a hardcoded end). LP DAAC flags reduced land-cover confidence for **2021+** years
  (degraded inputs) — carry that caveat in the dictionary. ~44 tiles cover CONUS + Alaska
  + Canada; ~10 GB for the full N. America archive — trivial.
- **Reproducibility**: pin **exact CMR granule IDs + production timestamps** in an input
  manifest, not just "year" — LP DAAC reprocesses, so reruns could otherwise silently
  pick different granules.
- **Format**: HDF-EOS2 / **HDF4**. ⚠️ The env's pip-GDAL has **no HDF4 driver** and
  `/vsicurl/`+Bearer fails on the LP DAAC S3 redirect → read SDS via **`pyhdf`** (proven in the
  LAI backfill): `SD(path).select("LC_Type1")[:]`; build the sinusoidal transform from tile h,v.
- **Bands (all stable, identical legend every year — confirmed)**: LC_Type1–5 +
  LC_Prop1–3 (+ assessment/QC/LW, not used for %). **Caution**: LC_Prop1/2/3 use
  **sparse, non-contiguous class codes** — enumerate from the v061 User Guide, don't
  assume `0..N`. **Fill = 255**, masked (excluded from the % denominator).
- **Per-class %**: exactextract `frac` op → coverage-weighted class area fraction →
  ×100. The full class-code list per band is a **committed manifest** (`legends.py` /
  dictionary), *not* discovered from observed pixels — reindex each band to that fixed
  list so missing (zero-count) classes become explicit `0.0` columns. This is what keeps
  the ~200-column schema stable across years/gages.

---

## 6. LAI source — MCD15A3H (ClimateAI catalog via GEOLATELY)

Proven reference (us-east-2): `wildfireModeling/.../ingest_mcd15a3h.py`, `modis_io.py`,
`MODIS_ACCESS_NOTES.md`. Adapt these directly.

- **Stored as reprocessed per-band COGs** (NOT HDF): `B01=Fpar, B02=Lai_500m,
  B03=FparLai_QC, B04=FparExtra_QC, B05=FparStdDev, B06=LaiStdDev`. We read **`B02`**
  (filter `URI.str.endswith("_B02.tif")`) + B03/B04 for QC. **No HDF4 needed for LAI**
  (only LULC needs it). 500 m, 4-day. Each composite date may span **multiple tiles** →
  group returned URIs by `START_DATE` + tile.
- **Coverage**: combined Terra+Aqua starts **~2002-07-04** → **start in 2002, no
  MOD15A2H splice** (confirmed). Runs to present.
- **Access** (us-east-2 catalog):
  ```python
  # materialize AWS creds into env BEFORE importing geolately (DuckDB httpfs reads them)
  from climateai.geolately.db_utils import Geospatial_DB
  db = Geospatial_DB()
  gdf = gpd.GeoDataFrame({"name":[aoi]}, geometry=[box(lon0,lat0,lon1,lat1)], crs="EPSG:4326")
  df = db.query_asset_catalog(product="MCD15A3H", start_date=s, end_date=e, query_gdf=gdf)
  lai_uris = df[df.URI.str.endswith("_B02.tif")]   # raster index → COG URIs; we aggregate
  ```
  Install GEOLATELY from private GitHub (`GITHUB_TOKEN` PAT) or its ECR docker image.
  COG CRS is read per-granule (`src.crs`) — the §7 engine reprojects polygons to it.
- **Fill masking (order matters — confirmed in reference)**: `Lai_500m` raw 0–100 valid;
  **mask raw ≥ 101 → NaN BEFORE** the ×0.1 scale (flags 249–255 sit above valid range;
  the reference masks all ≥101 conservatively), then defensively cap scaled LAI ≤ 10.
- **QC filtering (explicit bit rules; logic duplicated inline in `lai_pipeline.py` + `lai_backfill_lpdaac.py` — not yet a shared module/config)**: decode `FparLai_QC`
  (B03) — `MODLAND_QC` (bit 0) == 0, `SCF_QC` (bits 5–7) ∈ {0,1} (main RT algorithm)
  for the strict default; `FparExtra_QC` (B04: snow/ice, aerosol, cirrus, cloud-shadow)
  handling in `legends.py`, emitted in run metadata. Strict vs. backup-algorithm is a
  pilot/sensitivity decision (§9).
- **Monthly aggregation**: two-stage (§1b) — per-pixel day-weighted monthly mean, then
  **basin monthly mean + spatial heterogeneity (std, quantiles)** across basin pixels.

---

## 7. Shared zonal-stats engine (runs from us-east-2)

One engine, both products:

1. **Reproject polygons → the raster's native CRS** (read per-granule from `src.crs`):
   MODIS **sinusoidal** (`+proj=sinu +a=6371007.181 …`) for the LULC HDF granules;
   whatever the **catalog COG** reports for LAI — **never resample the categorical
   raster**. (For continuous LAI, reprojecting polygons is still preferred to avoid
   resampling artifacts.)
2. **Tile-batched**: for each MODIS tile, run exactextract over **all polygons
   intersecting that tile at once** (vectorized). **Multi-tile basins** (those crossing
   a sinusoidal tile boundary — common for large/Alaska/dateline basins) need care:
   - **Means & class fractions** combine exactly from per-tile coverage-weighted **sums
     and weight totals** (Σwx, Σw, Σwx²) → pooled mean/std/% with no error.
   - **Quantiles and min/max CANNOT be area-weighted from tile-level quantiles** (that
     is mathematically wrong). For multi-tile basins, **mosaic the basin's tiles to a
     per-basin VRT** (sinusoidal, lossless) and run a single exactextract, *or* pool the
     per-tile (value, coverage-weight) pixel samples and quantile once. The pipeline
     detects tile-crossing basins and routes them down this pooled path.
3. **LULC**: `frac` (per-class fraction) + coverage → % + QA counts.
4. **LAI**: Stage A builds per-pixel monthly-mean rasters per tile (QA-filtered,
   day-weighted composites per §1b); Stage B runs exactextract for `mean / stdev /
   quantile(.05,.25,.5,.75,.95)` (**coverage-weighted**) plus `min / max` (**unweighted
   extrema** — documented as such; optionally restrict extrema to pixels with coverage ≥
   a threshold to suppress slivers).
5. **Fill handling**: mask **raw** fill (255 for LULC; DN > 100 for LAI, **before**
   scaling) so NoData never pollutes weights/denominators.
6. **Small basins** (< ~1–2 pixels): exactextract partial coverage still yields a
   value from fractional pixels; emit `n_eff_pixels` and flag low-confidence rows
   (no upsampling — it adds no information).

---

## 8. Compute, cost & runtime (rough; confirm via pilot)

| | LULC (MCD12Q1) | LAI (MCD15A3H) |
|---|---|---|
| Files | 24 yr × ~44 tiles ≈ 1k tile-reads (HDF) | ~92 composites/yr × ~22 yr × ~44 tiles × 3 bands (B02/B03/B04 COGs) ≈ **~250k** reads |
| Volume / egress | ~10 GB **HTTPS egress** (cheap, out-of-region) | **in-region** (catalog us-east-2) → free |
| Compute | minutes–low hours | **the pole** — hours, heavily parallel per-tile; pilot to measure |
| Fallback | n/a | Terra-only **MOD15A2H** (8-day) ≈ halves the reads |

The pilot on a few tiles/months prices MCD15A3H and decides whether the Terra-only
fallback is needed. Note each composite read pulls **multiple bands** (`Lai_500m`,
`FparLai_QC`, and likely `FparExtra_QC`), so the read count above is a floor — the LAI
run must materialize **restartable per-(tile, month) intermediate artifacts** (Stage A
monthly-mean rasters / partial tables), not stream end-to-end.

---

## 9. Decisions

**Resolved (user, June 2026):**
1. **Language = Python** (GEOLATELY is Python-only; EO stack is Python-native). Unlike
   HydroATLAS (R), this product is Python.
2. **LAI monthly semantics** = two-stage (§1b): basin mean monthly LAI + spatial
   heterogeneity (std, quantiles).
3. **LULC width** = single wide table, **all bands** (LC_Type1–5 + LC_Prop1–3, ~200
   cols) + data dictionary.
4. **LAI start = 2002** (MCD15A3H combined), no MOD15A2H splice.
5. **Region** = run from us-east-2 (§3); no us-west-2 instance.
6. **Universe** = QA-passing gages (`processing_status == success`, 8,014), not all 16,994.
7. **Geometry priority** = **official polygons primary for all locations** (GAGES-II +
   WSC/ECCC); HydroBasins delineation is fallback only. `unified_watersheds_simplified.gpkg`
   is corrupt/unusable (Task 1 finding).

8. **Official-polygon sources + US ordering** — **US = GAGES-II first → NLDI fill →
   NHDPlus/StreamStats backstop; Canada = ECCC MDA_ADP**; HydroBasins final fallback
   (datasets/URLs in §4).

**Remaining (resolve during build / pilot):**
9. **LAI QC strictness** — strict main-algorithm (default) vs. allow backup; pilot.
10. **catalog COG CRS** — confirm at pilot.
11. **✅ Universe reconciliation — RESOLVED**: the 2,081 raw overlap was 100% a
    leading-zero artifact (canonical match → 5,707; golden ⊂ metadata-success). Most-inclusive
    universe = **8,018** (canonical union). Canonical key rule + dual `gage_id`/`canon_id`
    in §1. Geometry rebuilt to cover all 8,018.

---

## 10. Proposed module layout

```
EO_data_processing/
├── README.md                 # this plan / design doc
├── config/
│   └── eo_config.yaml        # paths, product params, QA thresholds, tile list
├── eo_processing/            # Python package
│   ├── geometry.py           # load watersheds; official+fallback merge; provenance; simplify
│   ├── catalog_lai.py        # GEOLATELY query for MCD15A3H
│   ├── fetch_lulc.py         # earthaccess search/open for MCD12Q1
│   ├── zonal.py              # shared engine: reproject + tile-batched exactextract
│   ├── legends.py            # MODIS class codes/names; LAI QC bit decoding
│   ├── lulc_pipeline.py      # → annual per-class % table
│   ├── lai_pipeline.py       # → monthly LAI stats table
│   ├── validate.py           # range / sum-to-100 / monotonic-quantile / area checks
│   └── io.py                 # CSV+parquet writers; data dictionary generator
├── notebooks/                # pilots (LAI compute pricing, small-basin checks)
├── viz/                      # self-contained QA/QC HTML explorers (Leaflet + Canvas)
│   ├── build_lai_explorer.py   # → lai_explorer.html  (S3: watershed_modis_lai_explorer.html)
│   └── build_lulc_explorer.py  # → lulc_explorer.html (S3: watershed_modis_lulc_explorer.html)
└── tests/
```
Outputs land in `data_out/` (gitignored), alongside the signature outputs.
(NOTE: the `eo_processing/` listing above is the original *proposal*; the as-built modules are
`lai_pipeline.py`, `lulc_pipeline.py`, `lai_backfill_lpdaac.py`, `lai_finalize.py`,
`lulc_finalize.py`, and `lulc_legends.csv` — see §11 for what each does.)

---

## 11. Task sequence (milestones)

1. **✅ DONE — Official polygons (PRIMARY geometry)** — GAGES-II (US, 6,160/6,160) +
   ECCC MDA_ADP (Canada, 1,823/1,854) → 99.6% of passing gages.
2. **✅ DONE — HydroBasins fallback + unified layer** — residual delineated (via
   `Downstream_HB_ID` or lat/lon→outlet BFS over `basinAt`) + 4 inclusive-universe USGS
   gages from GAGES-II. Universe = 8,018; **delivered geometry = 7,964** after excluding
   **54 basins with drainage >100,000 km²** (per user — 53 correctly-delineated huge
   Canadian/HB rivers + `05KH009`, an HB fallback that mis-delineated a 328,000 km² basin
   as 200 km² and slipped the geom-area filter; pruned 27 Jun). Unified per-gage layer,
   zero-padded `gage_id` + `canon_id`, `watershed_geom_source`. **Codex-reviewed (27 Jun)**:
   adaptive simplification (~200 m; 141 small/narrow basins kept full-res where simplify
   would shift area >2%); `geom_area_km2` recomputed from the **delivered** geometry
   (EPSG:6933); `low_confidence` flag (29 HydroBasins fallbacks + area outliers, incl. 3
   official CAN polygons `04HA001`/`06FD001`/`06FB001` whose reported `basin_area` ≫ official
   polygon — suspect metadata, kept + flagged). 0 dups/invalid. **Caveat**: US
   `area_rel_diff` is a consistency check, not independent validation — `basin_area` is
   GAGES-II-derived. **Deliverable on S3** `s3://climate-ai-data-science-shiny-app-data/streamflow/`:
   `watershed_polygons_26jun2026.{gpkg,parquet}` (EPSG:4326, 42/36 MB; sources: gagesii
   6,164 / wsc_eccc 1,771 / hydrobasins 29; cols incl. `geom_simplified`, `low_confidence`,
   `area_flag`) + `watershed_polygons_26jun2026_qa.csv`. Build scripts: `EO_data_processing/geometry/`.
3. **✅ DONE — Shared engine** — tile-batched, day-weighted monthly mean, VRT-mosaic +
   coverage-weighted exactextract (multi-tile safe). Lives inside `lai_pipeline.py`.
4. **✅ DONE — LULC pipeline + Codex-GO + on S3** — MCD12Q1 → annual all-band % table (29 Jun).
   `eo_processing/lulc_pipeline.py`: per-tile all-band HDF (LP DAAC **CMR search + curl
   `-H "Authorization: Bearer <~/.edl_token>"` + `pyhdf`** — pip-GDAL lacks HDF4 & `/vsicurl/`+Bearer
   fails) → uint8 sinusoidal GeoTIFF → VRT-mosaic per band (explicit nearest+nodata, lossless for
   aligned categorical) → `exact_extract(band, basins, ["unique","frac"])` → reindex to committed
   manifest (`eo_processing/lulc_legends.csv`, 102 class cols). **completeness guard** (write a year only
   if all cache-tiles downloaded *and* a tile counts only if all 8 bands decoded; else skip→re-run fills)
   + unknown-code hard-fail + search retry/backoff handle the Earthdata account-lock. Cache 28/30 tiles
   (2 ocean tiles = 0 granules). `eo_processing/lulc_finalize.py` adds
   `watershed_geom_source`/`low_pixel_support`/`low_confidence` + dictionary (idempotent; merge
   `validate="many_to_one"`). **Codex-reviewed 29 Jun (gpt-5.5) = GO** — no CRITICAL; output verified
   clean (all 8 bands sum to 100 within 1e-13, manifest matches v061 spec); 2 MAJOR latent code hazards
   hardened (output byte-identical). **Deliverable on S3** `s3://climate-ai-data-science-shiny-app-data/streamflow/`:
   `watershed_modis_lulc_annual_29jun2026.{parquet,csv}` + `_dictionary.csv` +
   `_granule_manifest_29jun2026.parquet` — **191,136 rows = 7,964 gages × 24 yr**, 109 cols, 0 dup keys,
   n_modis_pixels 100% populated. ~24 yr × ~28 tiles ≈ ~670 downloads.
5. **✅ DONE — LAI pipeline** — `EO_data_processing/eo_processing/lai_pipeline.py`, enhanced
   QA, Codex-GO. **Run complete + QA'd 28 Jun**: 2.14M rows, 7,964 gages, **269/270 months**,
   LAI 0.10–5.80 (0 out-of-range), monotonic quantiles, QA fracs ∈ [0,1], seasonality OK.
   **Codex results-review 28 Jun: values/geography/panel check out; NO-GO-as-is pending fixes.**
   The 109 always-NA basins = **92 no-catalog far-N** (catalog lacks those tiles → backfillable
   from LP DAAC) + **17 urban CONUS** (MCD15A3H Lai fill code 250 = built-up; MODIS retrieves no
   LAI over urban — permanent NA, NOT a bug, NOT backfillable; polygons confirmed correctly
   located). All 109 flagged in `out_qa/lai_always_na_basins.csv` (na_reason + backfillable flag).
   **Remaining**: (b) ✅ **DONE — backfilled + merged** via `eo_processing/lai_backfill_lpdaac.py`
   (LP DAAC: CMR search + curl-Bearer download + pyhdf; lock-aware downloader, workers≤2, never
   concurrent — after an Earthdata rate-limit incident). far-N (92 basins) + 2024-11 → **270-month
   full panel, 2,150,280 rows = 7,964×270, 0 dup keys**; only the 17 urban basins remain NA;
   pre-backfill backed up `_prebackfill.parquet`. **Codex-reviewed GO 29 Jun** (merge integrity +
   cross-source LP-DAAC/catalog consistency verified); `lai_always_na_basins.csv` regenerated to the
   17 urban. **✅ S3 upload DONE 29 Jun** (`watershed_modis_lai_monthly_29jun2026.parquet` +
   `_na_basins_29jun2026.csv` + `_dictionary_29jun2026.csv`). (c) ✅ **late-2024 partial flag DONE**
   — `partial_month` column (n_composites<6: 2024-10/12 + 2022-04/10). Winter NA = legit snow.
   Final product rebuildable via `eo_processing/lai_finalize.py` (merge + partial_month + NA-regen + dictionary).
   Command (re-run / resume):
   ```
   /home/sagemaker-user/.conda/envs/geo/bin/python \
     EO_data_processing/eo_processing/lai_pipeline.py --years 2002-2024 --workers 5
   ```
   Defaults: geometry `data_out/eo_lai/geometry_7964.gpkg`, outdir `data_out/eo_lai/out_qa`
   (persistent `/home`). **Resume after a crash: re-run the SAME command** — `process_month`
   skips done `out_qa/lai_YYYY_MM.parquet`, `build_uri_cache` reuses `out_qa/_uri_cache.parquet`
   (full range). **Do NOT delete `_uri_cache.parquet` on resume** (only if changing year range —
   it returns any existing cache blindly). Final: `out_qa/watershed_modis_lai_monthly.{parquet,csv}`.
   **Creds (critical, fixed 28 Jun)**: workers pop `AWS_*` env vars so boto3 uses the
   auto-refreshing `container-role` — without it, static env creds expire ~30 min in and
   silently empty every month after the first ~21 (symptom: a long tail of `: empty` despite a
   full-range cache). **If re-running after such a failure, delete the empty `lai_*.parquet`
   first** (skip-if-exists would otherwise keep them) — empties have only a `gage_id` column.
   **Post-run (DONE)**: far-N (92) + 2024-11 backfilled via `lai_backfill_lpdaac.py`, merged +
   flagged via `lai_finalize.py`, uploaded to S3 — see (b).
6. **✅ DONE — Outputs, data dictionary, validation** — both LAI (monthly) and LULC (annual)
   tables + dictionaries written, QA'd, Codex-reviewed (GO), and delivered to S3. gage_id join to
   the signatures confirmed (leading-zero-safe `gage_id`/`canon_id`). **EO pipeline COMPLETE.**
7. **✅ DONE — QA/QC HTML explorers (both products) + on S3** — self-contained Leaflet maps for
   manual human review (`viz/build_lai_explorer.py`, `viz/build_lulc_explorer.py`; outputs in
   `data_out/eo_lai/`, on S3 as `watershed_modis_lai_explorer.html` / `watershed_modis_lulc_explorer.html`).
   LAI: 7,964 points by 17 QA/summary vars + per-gage monthly time series (quantile bands +
   good-coverage strip). LULC: 27 map variables (10 IGBP-derived summary + 17 individual IGBP classes;
   overlapping static `pct_*` roll-ups intentionally dropped) + per-gage stacked-area composition
   2001–2024 with a band switcher across all 8 schemes; glossary carries the MODIS caveats (SE-US
   woody-savanna labeling; 2021+ reduced confidence). Both headless-validated (puppeteer/Chrome).
   **EO pipeline + both explorers COMPLETE.**

Environment (runs from **us-east-2**): Python, **`pyhdf`** for HDF4 reads (env pip-GDAL has NO HDF4 driver),
`earthaccess rioxarray rasterio exactextract geopandas pyproj shapely pandas pyarrow
numpy`, GEOLATELY (private GitHub PAT / ECR image), Earthdata Login `.netrc`, ClimateAI
AWS SSO. The `wildfireModeling` conda env `geo` already bundles
rasterio/geopandas/s3fs/`climateai.geolately`/pyarrow — reuse or mirror it (add
`earthaccess` + `exactextract`).

---

## 12. Operational controls & reproducibility

- **Credentials** — materialize AWS creds into env vars **before importing geolately**
  (DuckDB httpfs reads them at import); call a `refresh_aws_env()` at run start and
  periodically — SageMaker creds rotate ~30 min (pattern in `wildfireModeling/.../lib.py`).
- **Input manifest** — record exact CMR granule IDs + production timestamps (LULC) and
  catalog URIs + `REGION` (LAI) per run, so reruns are deterministic and reprocessed
  granules don't silently change outputs.
- **Checkpointing / restart** — per-(tile, year) (LULC) and per-(tile, month) (LAI)
  intermediate artifacts; idempotent writes so a failed long LAI run resumes without
  recompute. Retry policy for transient S3/HDF read errors.
- **Versioned config + env lock** — `eo_config.yaml` (products, QC bit rules, tile list,
  QA thresholds, simplification tolerance, geometry version) under version control; a
  pinned environment lockfile.
- **Validation gates (`validate.py`)** — LAI ∈ [0,10]; monotonic quantiles
  (q05≤q25≤…≤q95); per-band LULC class %s sum ≈ 100 over valid pixels; expected
  (gage × year) / (gage × year × month) row completeness; `watershed_area_rel_diff`
  within threshold; 100% leading-zero-safe `gage_id` join to the signatures. Emit a QA
  report per run.

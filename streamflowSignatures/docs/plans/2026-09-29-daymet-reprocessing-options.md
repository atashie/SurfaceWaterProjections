# Daymet climate input — reprocessing options and feasibility plan

**Date**: 2026-09-29 · **Status**: options review + test plan (no data pulled — the
polygons are not on this machine) · **Author**: Claude Code session for A. Tashie

**Action plan (how to run it on the dedicated machine)**: `2026-09-29-daymet-reprocessing-action-plan.md`.

**Update 2026-09-29 (later)**: option A was run for calendar 2023 on the M5 MacBook (T1, most
of T0) and passed — the job is download-bound (≈ 50 h for 1980–2025 on this link), the
co-authors' values are reproduced once cells are weighted by true area over unsimplified
polygons, and T2 (AWS) / T3 (Earth Engine) are no longer needed for feasibility. User
decisions the same day: full-resolution polygons, all six variables, the user's EDL token
for 2025, and the 53 basins > 100,000 km² if RAM allows. Details: action plan §0.

**Purpose.** Replace the stale climate input (`daymet_1980_2023` → 6,087 basins, calendar
1980–2023) with a fresh Daymet extraction that (1) reaches the latest published year and
(2) covers **all 7,964 watersheds** of the delivered boundary layer, closing the hole that
leaves ~1,160 (WY 1993–2025) / ~610 (WY 1980–2025) product gages without climate, snow
and precipitation signatures.

Everything marked **[V]** was verified on a live page on 2026-09-29 by the four research
sub-agents (URLs inline or in §8); **[U]** is unverified. Repo facts come from the docs
and the local metadata CSV.

---

## 0. Bottom line

1. **Daymet V4 R1 now ends at calendar 2025** (mosaics released 2026-05-22) [V]. Calendar
   2024 and 2025 are available; **calendar 2026 will not exist before ~spring 2027**, so
   the achievable target is **WY 1980–2025** — which is exactly where both streamflow
   products end. Nothing can supply WY 2026 climate today.
2. **The access landscape changed under us.** ORNL's THREDDS server, its NetCDF Subset
   Service and the 2°×2° tiles are gone (2025) [V]; every convenience library the
   community relied on (`pydaymet`, `daymetr`, `climateR`, the gdptools `daymet4` catalog
   entry) is broken for gridded pulls and none has shipped a fix [V]. What remains at
   scale: **whole-file annual mosaics** (ORNL over HTTPS with an Earthdata Login, or S3
   in `us-west-2`; **plus a login-free NCAR mirror through 2024 with Globus**) and
   **Google Earth Engine** (through 2025-12-31) [V].
3. **No published product substitutes for a recompute** — nothing carries basin-averaged
   Daymet with SWE for more than ~1,000 US basins, none on our polygons, none keyed to
   HYDAT ids at scale [V]. A full recompute over our own polygons is unavoidable.
4. **The coverage hole is not a basin-size problem.** Only 58 of the 8,014 usable gages
   exceed the 85,000 km² Daymet cap, yet 2,049 lack Daymet; the 6,087-basin candidate list
   simply predates the 7,964-polygon geometry layer (the manuscript's new §2.1.3 says so).
   Reprocessing on the 7,964 layer closes it.
5. **Recommended path**: our own zonal statistics over the whole-file mosaics with
   `exactextract` (coverage weights computed once for the 7,964 polygons), in two
   tranches — **prcp + swe first** (≈ 0.57 TB of input; the two variables every climate
   and snow signature needs), the four secondary variables later — validated by
   reproducing the co-authors' 2023 values on the 6,087 overlapping basins. Run the
   feasibility tests in §4 first (≈ 1 week of calendar time, mostly waiting on transfers);
   T1 alone decides whether a laptop can do tranche 1 or whether a `us-west-2` instance /
   Earth Engine is needed.

---

## 1. Current state and what "missing" actually means

| Fact | Value | Source |
|---|---|---|
| Current input | 44 annual CSVs `site_id, month, year, prcp, tmin, tmax, swe, vp, srad` (no day column), CY 1980–2023 → `daymet_1980_2023_rebuilt_10aug2026.parquet`, 97,757,220 rows | DATA_SOURCES.md, DEVELOPMENT.md |
| Basins covered | 6,087 sites (122 never compiled) → **5,965 of 8,014** usable gages; **5,517 of 6,678** (product #1) and **5,638 of 6,250** (product #2) | census 2026-09-04 |
| Aggregation | co-authors' USGS `gdptools` (`agg_gen`), area-weighted fractional overlap, basins < 85,000 km² | manuscript §2.1.3 |
| Signatures that need it | PPT: runoff ratios (+`runoff_ratio_high_count`), elasticity (+2 diagnostics), Q-P seasonality, `avg_storage`; SWE: the 14 snow metrics. tmin/tmax/vp/srad are used by no signature today (future ET work; HydroShare Resource 3 table only) | SIGNATURES.md |
| Product gages without climate | #1: **1,161** (1,006 US / 155 CAN); #2: **612** (494 US / 118 CAN) | census §2 minus §1 |
| Gages > 85,000 km² among the 8,014 | **58** (52 CAN > 100k, 6 CAN 85–100k, 0 US) | local metadata, 2026-09-29 |
| Basins > 100,000 km² excluded from the boundary layer | 54 → the 7,964 layer | census |

So the size cap explains at most 58 of the 2,049 missing gages. The rest were never in
the co-authors' candidate list. Consequence: the gap audit (T0) must run on the actual
site list, but the fix is the same either way — compute for all 7,964 polygons. Whether
to *keep* the 85,000 km² cap as a downstream filter is a separate decision (§6).

Two inherited facts a reprocess must preserve or consciously change:
- **365-day calendar**: Daymet drops Dec 31 in leap years [V]. The preprocessor absorbs
  the one-day hole; keep the source calendar faithful and document it (as now).
- **Snow metrics only see SWE-valid years**; 32 / 3 Daymet gages in the products have no
  SWE-valid year in window — a data property, not a coverage gap.

---

## 2. Daymet today (verified 2026-09-29)

### 2.1 Product status
- **Latest year: calendar 2025**, released **2026-05-22**; the V4 R1 daily record is 1,176
  granules / 3.205 TB (7 variables × NA 46 yr + HI 46 yr + PR 76 yr) [V]
  (https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-daymet-daily-v4r1-2129-4.1;
  release table in the user guide
  https://data.ornldaac.earthdata.nasa.gov/public/daymet/Daymet_Daily_V4R1/comp/Daymet_Daily_V4R1.pdf).
- **Cadence**: 2022 → 2023-03-01; 2023 → 2024-04-23; 2024 → 2025-09-12 / 10-01 (slipped:
  ECCC input-feed problems); 2025 → 2026-05-22. Nominally March–May of year N+1 [V]. The
  monthly-latency product ended 2023-03 and is access-restricted; a Single Pixel request
  for Jan 2026 returns nothing [V]. **Plan on no 2026 data before ~spring 2027.**
- **R1 vs V4**: only the 2020 and 2021 files changed (high-latitude fix); "files outside
  of 2020 and 2021 have not changed" [V]. Our pull reaches 2023, so it is R1 already.
- **Canadian input change**: for 2024 and 2025 the GHCNd `CA0` stations were not used;
  ECCC/CCCS stations were substituted and bias-corrected [V]. Treat 2024–2025 as a
  station-network seam for Canadian trend endpoints (document in run notes and the
  Resource 3 README).
- No Daymet V5 or discontinuation notice found [V, absence]. `daymet.ornl.gov/news`
  returned HTTP 500 [U].

### 2.2 Access routes (what still works)

| Route | Status | Scale fit | Notes |
|---|---|---|---|
| **NA daily mosaics, one NetCDF-4 per variable-year** (`daymet_v4_daily_na_<var>_<yyyy>.nc`) | **Live** [V] | The only bulk route | Sizes (2025 / 1980, MB): prcp 4,360 / 3,689 · swe 7,646 / 8,775 · tmin 14,019 / 11,405 · tmax 13,735 / 11,645 · vp 17,023 / 21,251 · srad 15,852 / 18,648 · dayl 1,484 / 1,290 [V]. **prcp+swe ≈ 12.5 GB/yr ≈ 0.57 TB for 1980–2025; all six ≈ 73 GB/yr ≈ 3.3–3.5 TB.** Internal chunking undocumented [U] — measure in T1. Grid 7,814 × 8,075, LCC. |
| **NCAR GDEX d682806** "convenience copy of version 4", NA/HI **1980–2024**, 3.45 TB, CC BY 4.0, yearly updates | **Live, anonymous** [V] | Bulk; **Globus, HTTPS, THREDDS/OPeNDAP + subsetting** | Same file names; `fileServer/.../daymet_v4_daily_na_prcp_2024.nc` → 200, 4,471,051,114 B, modified 2026-02-19 [V]. **No 2025 yet.** Global attribute reads "Version 4.0" in both 2021 and 2024 files, so R1-ness of its 2020/2021 files is [U] → verify by size/checksum against ORNL's CMR record (T1). https://gdex.ucar.edu/datasets/d682806/ |
| ORNL HTTPS with Earthdata Login (`.netrc` + curl/wget, or `earthaccess`) | Live [V] | Bulk, from anywhere | `https://data.ornldaac.earthdata.nasa.gov/protected/daymet/Daymet_Daily_V4R1/data/<file>`; CMR short_name `Daymet_Daily_V4R1_2129`, concept-id `C2532426483-ORNL_CLOUD` [V]. No published rate limit [V, absence]. Needed for the 2025 files. |
| **ORNL S3 direct** `s3://ornl-cumulus-prod-protected/daymet/Daymet_Daily_V4R1/data/` | Live [V] | Bulk, in-region only | `us-west-2`; 1-hour temporary credentials from `https://data.ornldaac.earthdata.nasa.gov/s3credentials`; read-only; same-region compute required [V]. Zero egress from a `us-west-2` instance. |
| Google Earth Engine `NASA/ORNL/DAYMET_V4` | Live, **through 2025-12-31**, R1 2020/2021 [V] | Whole job server-side | Coverage-weighted reducers, weights quantised to 1/256 [V]. Licensing and quotas: §3 D. |
| THREDDS / NCSS / OPeNDAP aggregate `daymet-v4-agg/na.ncml` / 2°×2° tiles | **Retired** (THREDDS "permanently down" Aug 2025; tiles "no longer available" Jan 2026) [V] | — | `daymet.ornl.gov/getdata` still advertises them (stale) [V]. |
| Earthdata Hyrax OPeNDAP (per granule) | Live; index-range slicing only, per file, EDL login [V] | Poor | Hundreds of requests per variable-year; 2024 users saw 2 kb/s–5 Mb/s and timeouts [V]. |
| NASA Harmony subsetting | "should not be used" for Daymet (ORNL, Aug 2025) [V] | — | `grid_mapping` errors; > 1 year fails. |
| AppEEARS (`DAYMET.004`, 1950–2025) | Live [V] | Unknown for 8k polygons [U] | Area requests return clipped grids, not means; limits undocumented. |
| Single Pixel Extraction API | Live, serves 2025 [V] | Points only | Spot-checks: `https://daymet.ornl.gov/single-pixel/api/data?lat=&lon=&vars=&start=&end=`. |
| Microsoft Planetary Computer Zarr (`daymet-daily-na`, westeurope) | Live but **ends 2020-12-30**, pre-R1; SAS token required [V] | Stale | Chunks [365, 284, 584]. Not useful. |
| Zarr / Kerchunk / VirtualiZarr references | None published (pangeo-forge recipe WIP 2022, abandoned) [V] | — | Only a NASA Openscapes Kerchunk tutorial (shows the 1-h credential expiry) [V]. |

### 2.3 Libraries and reference code
- **pydaymet** 0.19.4 (2025-05): `get_bygeom` used NCSS → HTTP 401 since Sept 2025; issue
  #72 open, maintainer "looking for suitable alternatives"; `get_bystac` = Planetary
  Computer (2020 only); Alaska never supported [V].
- **daymetr** 1.7.1: `download_daymet_ncss` (6 GB cap) dead; README says migrate to
  `{appeears}`; single-pixel still works [V].
- **climateR** `getDaymet()`: OPeNDAP → "Access denied" (issue #113, open); the
  September-2026 climateR-catalog release has **zero** Daymet rows [V]. **zonal**
  (`weighting_grid` once, `execute_zonal` over annual files) remains a sound pattern [V].
- **gdptools** 0.4.0 (2026-09-02): very active — `exactextract` engine (0.3.4), its own
  zonal engines deprecated in favour of it (0.3.11), memory-safe chunked mode,
  `estimate_memory_gb()`, `UserCatData` for local NetCDF/Zarr/OPeNDAP; parallel with
  batching recommended for 5,000+ polygons; USGS ran it over ~110,000 NHM HRUs with cached
  weights on SLURM (128 GB tasks); no published timing for thousands of polygons × decades
  of daily data [V]. Its `daymet4` ClimateR-catalog entry was a THREDDS source — gone [V].
- **exactextract** (Python 0.3.0): multi-band rasters in one call (`band_1_mean, …`),
  `coverage`/`cell_id` ops for a one-time weight matrix; maintainer's benchmark: ~55 % of
  time is coverage computation with one band, ~5 % with ten — **stack many days per
  call** [V] (https://github.com/isciences/exactextract/issues/183). Already used in
  this repo (MODIS, NLCD).
- **CAMELS-SPAT** processed Daymet for 1,426 US+CA basins the whole-file way: download
  the annual NA files (10 threads), subset per basin with `datatool`, area-average with
  **EASYMORE** (remap weights once; no swe) [V]
  (https://github.com/CH-Earth/camels_spat, `7_forcing_data/6a_download_daymet.py`,
  `6d_daymet_to_basins.py`). EASYMORE: 500,000 subbasins × 7 variables in 1.2 s per time
  step once the remap file exists [V].
- ORNL's own notebooks (`ornldaac/daymet`, `daymet-python-opendap-xarray`) cover bbox
  subsetting via Hyrax/earthaccess (365 requests per variable-year) — no polygon
  averaging [V]. Earth Engine's official `reduceRegions` example is Daymet V3 × HUC06;
  geemap #1188 reports per-image `reduceRegion` exports "very slow" for ~3,000 polygons
  and `toBands()` + zonal stats as the fix [V].
- **Existing basin-averaged products** (all [V]): HYSETS (no Daymet; ERA5-Land, NRCan,
  Livneh, SCDNA), CAMELS-SPAT (≤ 1,426, no swe, MERIT polygons), CAMELS (671, 1980–2014),
  Caravan / CAMELSH (ERA5-Land), MacroSheds (169 research sites), NHM-PRMS by HRU (P/T
  only), MACH (1,014 US, NHDPlus polygons, to 2023), River & Floodplain (505),
  Kovacek/Borealis (BC region, CC BY-SA), BASINGRID (9,067 GAGES-II basins but 64×64
  raster patches, 1985+). None substitutes.

---

## 3. Options

Common to A–D: the polygon layer is the delivered 7,964-basin layer (Resource 4,
EPSG:4326; Drive backup — S3 is gone), reprojected to Daymet's LCC grid
(`+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m`)
so the raster is never resampled; output = one area-weighted daily value per basin per
variable, 365-day calendar preserved.

### A. Whole-file mosaics → own zonal statistics on our hardware (GDEX + ORNL)
- **Source**: 1980–2024 from **NCAR GDEX via Globus** (anonymous, resumable, built for
  multi-TB transfers) or HTTPS; 2025 (two files for tranche 1, six for tranche 2) from
  ORNL over HTTPS with an Earthdata Login. Verify GDEX's 2020/2021 files against ORNL's
  CMR sizes/checksums once (R1-ness) [U until T1].
- **Compute**: for each variable-year file either (a) `exact_extract(year_stack, basins,
  ["mean"])` with the 365 days as bands, or (b) a sparse pixel×polygon coverage matrix
  built once (`coverage` + `cell_id`) and `W @ day_slice` per day — the CAMELS-SPAT /
  EASYMORE / zonal pattern; nested basins cost nothing extra in (b) because each pixel is
  read once. A full NA day slice is 63 M px ≈ 250 MB uncompressed; tranche 1 is
  365 × 46 × 2 ≈ 8.4 TB of decompressed reads — hours to a couple of days of I/O on an
  M1 or the 16 GB laptop, memory bounded by chunked reads.
- **Data**: tranche 1 prcp+swe ≈ 0.57 TB; tranche 2 the other four ≈ 2.8 TB. Stream or
  stage on an internal/USB-C SSD — **never the exFAT thumb drive** (three silent
  truncations already).
- **Cost**: $0 (a 2–4 TB SSD if tranche 2 is done locally). **Effort**: ~2 days of code
  (reuse the repo's exactextract pipelines and the structural checks of
  `convert_daymet_csvs_to_parquet.py`). **Risk**: transfer time; laptop tied up for
  days; polygon-version mismatch with the co-authors' run (T1 validation).
- **Verdict**: the default for tranche 1; T1 decides.

### A-addendum — measured streaming budget (2026-09-29, from this Windows laptop)

Probed remotely on the NCAR mirror (HDF5 headers via byte-range reads; h5py) and the
mirror's catalog; decompression timed locally; download timed with curl.

| Quantity | Measured |
|---|---|
| File layout **1980–2019** (`prcp`; same pattern expected for all vars) | float32 (365, 8075, 7814), chunks **(1, 1000, 1000)**, gzip-4, no shuffle → 26,280 chunks/file, 4 MB uncompressed each |
| File layout **2020–2025** | chunks **(10, 300, 300)**, gzip-4, shuffle → 26,973 chunks/file, 3.6 MB each. The layout switch at 2020 is the R1 remake, so the mirror's 2020/2021 files are R1 (supports, does not prove, byte-identity with ORNL — still size/md5-check once) |
| Uncompressed size per variable-year | 63.1 M cells × 365 × 4 B = **92 GB** (never load whole) |
| Compressed totals 1980–2024 (mirror catalog, 315 NA files) | prcp 0.174 TB · swe 0.390 · tmin 0.521 · tmax 0.533 · vp 0.938 · srad 0.825 · dayl 0.059 → **six vars 3.38 TB, prcp+swe 0.564 TB**; per year 72–77 GB (six) / ~12.5 GB (prcp+swe); 2025 from ORNL adds ≈ 72.6 / 12 GB |
| Download throughput, mirror → this machine | **18.3 MB/s** single stream (300 MB ranges, twice); 4 parallel streams ≈ 25 MB/s aggregate → the link, not the server, is the cap here |
| Decompression on this laptop | **150–170 MB/s uncompressed per core** (Python zlib + unshuffle; prcp ratio 8×, vp 3.7× on land tiles) → ≈ 9–10 min per variable-year per core reading every tile |
| Basin pixel budget (metadata, 7,889 basins ≤ 100k km²) | Σ area = **25.6 M km²** ≈ nnz of the coverage matrix (≈ 205 MB as CSR int32 + float32); 662 basins > 10k km², 36 > 50k, max 93,900 |

**Per-year cycle (download → process → parquet → delete), six variables:**
- Download ≈ 50–70 min (72–77 GB at 18–25 MB/s). This is the bottleneck.
- Compute ≈ 6–15 min on 4 worker processes: read only the spatial tiles that intersect any
  basin (est. 30–50 % of the grid), chunk-aligned (`read_direct_chunk` / chunk-boundary
  slices), accumulate `W_tile @ chunk` into a 7,964 × 365 array per variable; the sparse
  matvec itself is ≈ 30–50 ms per day (≈ 2 min/yr for six vars). Conservatively 30 min.
- Pipelined (download year N+1 while processing N) the run is download-bound:
  **≈ 2–2.5 days for six variables, ≈ 9–10 h for prcp+swe**, unattended.
- Output ≈ 17.4 M values/yr → **~50–70 MB parquet per year, ~3 GB for 1980–2025**.

**RAM (peak, chunk-aligned design):** coverage matrix ~0.2–0.3 GB (+ tile index) shared or
sliced per worker; in-flight chunks 4 MB each (tens per worker); accumulators 23 MB per
variable (140 MB for six); Python/numpy/h5py baseline ~0.3 GB per process → **≈ 2–3 GB
with 4 workers, ≈ 4–5 GB with 8**. Known spike traps: a whole-day slice is 252 MB and a
10-day slab (2020+ chunking) 2.5 GB per process — acceptable but avoidable; `var[:]` on a
file is 92 GB — never. The exactextract-multi-band alternative (GDAL netCDF windows per
polygon) peaks at ~1–2 GB per process but re-decompresses shared chunks for nested basins,
so it is slower, not lighter.

**Disk:** one year of raw input 72–77 GB (six) / 12.5 GB (prcp+swe); pipelined ≈ 2× that
(≈ 155 GB / 25 GB) plus ~3 GB of outputs and ~0.3 GB of cached weights. Variable-at-a-time
download→process→delete bounds it at ≈ 2 × the largest file (vp ≈ 21 GB) ≈ 45 GB.
**Internal or USB-C SSD only** — the exFAT thumb drive has already silently truncated
three files.

### B. Same computation, in AWS `us-west-2` next to the ORNL bucket
- **What**: a spot instance (8 vCPU / 32–64 GB, 1–2 TB gp3 scratch) reads the mosaics via
  `earthaccess` S3 credentials (renew hourly), runs the same code, ships back only the
  results (7,964 × 16,800 × 6 ≈ 0.8 G values ≈ 3–6 GB parquet).
- **Cost**: order of $20–60 of spot instance time for all six variables plus ~$1 egress
  [U — check current pricing]. **Effort**: as A plus account/instance setup (the
  project's previous AWS access is gone; a new account or a co-author's). **Risk**:
  credential expiry mid-run (re-auth loop), spot interruption (checkpoint per
  variable-year).
- **Verdict**: the route for tranche 2 (2.8 TB) and the fallback for tranche 1 if T1 shows
  the transfer is too slow.

### C. gdptools by the USGS co-author (reproduces the existing product)
- **What**: the co-author who ran the 2023 aggregation re-runs it: `UserCatData` on the
  downloaded NA files, `WeightGen` once on the 7,964 layer, `AggGen` with the
  `exactextract` engine, parallel + batching, `agg_data_2_parquet`. gdptools 0.4.0 is
  current and is the tool the manuscript names.
- **Cost**: their time and compute. **Effort**: an email plus a shared polygon file.
  **Risk**: availability before the Nov 9 submission; we cannot inspect the run.
- **Verdict**: ask now regardless — even if we compute in-house, their script, polygon
  set and runtime are the ground truth for T1's validation.

### D. Google Earth Engine `reduceRegions` on `NASA/ORNL/DAYMET_V4`
- **What**: server-side daily basin means, exported per year as tables; no download, no
  local compute; covers Alaska/Canada and 2025; would make the annual refresh trivial.
- **Licensing** [V]: noncommercial tiers are for academic staff, individuals for
  noncommercial purposes, nonprofits; Google states "if you work at a private company,
  you need to configure a paid (commercial) Earth Engine account" and there is "no
  official grace period". The employee-of-a-company-on-an-academic-project case is not
  addressed anywhere → **run under an academic co-author's registered noncommercial
  Cloud project**, or pay (Basic $500/month; batch $0.40/EECU-h on the usage-fee plan).
- **Quotas** [V]: Community tier 150 EECU-h/month, Contributor 1,000 (billing account on
  file, not charged); ~2 concurrent batch tasks, 3,000 queued, 10-day task lifetime.
  `maxPixels` default 1e7, `tileScale` for memory. Per-image `reduceRegion` exports are
  slow; `toBands()` per year + one `reduceRegions` is the documented fix.
- **Precision**: coverage weights quantised to 1/256 — immaterial for basins of hundreds
  of pixels; check the 1,199 US gages ≤ 100 km².
- **Verdict**: strong second option; cheap to test (T3).

### E. Incremental append only (2024–2025 for the existing 6,087 basins)
≈ 145 GB for six variables, ~25 GB for prcp+swe; needs the co-authors' polygons to be
consistent. **Does not close the coverage hole** and mixes two aggregation runs in one
series. Only a stopgap if the full recompute cannot land before submission.

### F. Not viable (for the record)
pydaymet / daymetr / climateR gridded pulls (THREDDS gone); Planetary Computer (2020);
Hyrax / Harmony / AppEEARS / Single Pixel at 7,964 polygons × 46 years; any existing
basin-averaged product (§2.3).

### Comparison

| | Currency | Input moved | $ | Engineering | Reproduces co-authors' run? | Main risk |
|---|---|---|---|---|---|---|
| A GDEX+ORNL mosaics, local | 2025 | 0.57 TB (tranche 1) / 3.4 TB | 0 (+SSD) | ~2 d | validate on 6,087 | transfer time, disk |
| B ORNL S3, `us-west-2` | 2025 | 0 (in-region) | ~$20–60 [U] | ~2 d + AWS | same code as A | account, creds, spot |
| C gdptools (co-author) | 2025 | their side | 0 | an email | yes (same tool) | availability |
| D Earth Engine | 2025 | 0 | 0 if eligible, else $500/mo | 1–2 d | approximates | licence, 1/256 weights |
| E append 2024–25 | 2025 | 25–145 GB | 0 | 1 d | partial | leaves the hole |

---

## 4. Feasibility tests (each ≤ 1 day of work; T0 first, T1–T4 in parallel)

### T0 — prerequisites and gap audit (½ day, no bulk downloads)
1. Earthdata Login account; `.netrc` on the Mac; confirm a 1 MB range request from
   `…/protected/daymet/Daymet_Daily_V4R1/data/daymet_v4_daily_na_dayl_2025.nc` succeeds.
   Install Globus Connect Personal; confirm the GDEX endpoint lists `d682806`.
2. Record CMR granule ids, sizes and `updated` dates for every ORNL file we will use, and
   GDEX sizes/mtimes for the 1980–2024 files (provenance manifest, as the LULC pipeline
   does). Sizes that differ between GDEX and CMR for the same file = stop and check.
3. Retrieve the boundary layer and the Daymet site list: the Resource 4 layer (Drive
   backup, 7,964 rows) and `SELECT DISTINCT site_id FROM
   daymet_1980_2023_rebuilt_10aug2026.parquet` (6,087). Join on the zero-stripped id
   (never re-pad — 9 boundary ids are mis-padded).
4. **Gap audit**: classify the polygons without Daymet (≈ 1,999) by country,
   `watershed_geom_source`, area bin, and product membership; list the 122 never-compiled
   Daymet sites. Output: one CSV + a paragraph for the manuscript's §2.1.3 (which now
   claims the gap is a timing artefact — this audit proves or corrects it).
5. Reproject the 7,964 polygons to Daymet LCC; validity check; pixel count per basin
   (`geom_area_km2` ≈ pixels); flag basins < 4 pixels (weighting sensitivity).
- **Pass**: logins work; join reproduces 5,965; audit CSV written.

### T1 — Option A on one year (1 day; decides tranche 1)
1. Transfer `daymet_v4_daily_na_prcp_2023.nc` (≈ 4.4 GB) and `…_swe_2023.nc` (≈ 8 GB)
   from GDEX (Globus and HTTPS — **record MB/s for each**) and the same two 2025 files
   from ORNL HTTPS (record MB/s). Extrapolate 0.57 TB and 3.4 TB.
2. Verify GDEX = ORNL: sizes of the 2020 and 2021 prcp files vs CMR; md5 of one file from
   both sources if sizes match.
3. Inspect chunking (`ncdump -hs` / xarray `encoding`), CRS and grid metadata.
4. Run both zonal strategies on all 7,964 polygons: (a) `exact_extract` with 365 bands;
   (b) coverage matrix once (`coverage` + `cell_id`) then `W @ slice` per day. Record
   wall-clock, peak RSS, and which is faster on this file layout.
5. **Validate against the co-authors' 2023 values** for the 6,087 overlapping basins:
   per-basin identity R² and max |Δ| for daily prcp and swe; annual totals; count basins
   with R² < 0.999 and inspect them (expected culprits: different polygon version, basins
   < 4 pixels, coastal slivers). This is the acceptance test for *any* option.
6. Spot-check five basins against the Single Pixel API at their centroids (sanity).
- **Pass**: ≥ 95 % of shared basins at R² ≥ 0.999 on daily prcp with explainable
  residuals; per-variable-year wall-clock × 92 (tranche 1) fits in ≤ 3 days on the Mac;
  projected transfer ≤ 1 week. Transfer the only failure → B.

### T2 — Option B on the same year (½ day + AWS access)
1. `us-west-2` spot instance (8 vCPU / 32 GB / 500 GB gp3), same environment, S3
   credentials via `earthaccess`; time (a) `aws s3 cp` of one file to local disk and (b)
   `s3fs` + `h5netcdf` direct reads.
2. Re-run T1 step 4 there; extrapolate cost = instance-hours × price for 92 (tranche 1)
   and 276 (all six) variable-years, plus results egress.
- **Pass**: ≤ 30 min per variable-year end-to-end and total ≤ ~$100 [U pricing]; results
  bit-identical to T1 (same code).

### T3 — Option D on Earth Engine (1 day)
1. Confirm eligibility in writing; register (or borrow) an academic co-author's
   noncommercial Cloud project.
2. Upload the 7,964 polygons as an asset; for 2023: `ImageCollection('NASA/ORNL/DAYMET_V4')
   .filterDate(...).select(['prcp','swe']).toBands()` → one `reduceRegions(mean,
   scale=1000, crs=<Daymet LCC>)` → export table. Also 200 basins for 2020–2025 and all
   7,964 for one month, to measure queue and EECU-h per year.
3. Compare with T1 values (identity R² per basin; quantify the 1/256 effect on the
   smallest basins) and with the co-authors' 2023 values.
- **Pass**: a year-task completes in ≤ a few hours within quota; agreement with T1 within
  the quantisation noise; licence question answered.

### T4 — Option C: ask the co-author (an email; answers within the week)
Ask for: the gdptools script and version used, the polygon file and its date, the
candidate list (why 6,087 / who was excluded), runtime and machine, willingness to re-run
1980–2025 on the 7,964 layer, and whether they can run the T1 year so the two
implementations can be diffed.
- **Pass**: script + polygons in hand (T1 validation becomes exact), or a committed
  re-run date before mid-October.

### Decision (target: go/no-go by ~2026-10-10)
- T1 pass → **A for tranche 1 immediately** (prcp+swe, 1980–2025, all 7,964), then B or A
  for tranche 2.
- T1 fails only on transfer speed → **B**; fails on agreement → investigate polygons
  first (T4), then decide.
- T3 passes and licensing is clean → GEE becomes the annual-refresh mechanism (2026 data
  in 2027) even if A/B produce this year's deliverable.

---

## 5. Output contract (so nothing downstream changes)

- **Parquet** `daymet_1980_2025_<date>.parquet`: `site_id` (zero-padded agency id as a
  STRING — the join key the runner uses), `Date`, `prcp`, `tmin`, `tmax`, `swe`, `vp`,
  `srad` (+ `dayl` if wanted), one row per site-day, **365 rows per site-year** (Dec 31
  absent in leap years, as in the source). Tranche 1 may ship as a prcp/swe-only file;
  the runner keeps only `PPT`/`SWE` (`run_julia_benchmark.jl` intersects
  `gage_id, date, PPT, SWE`).
- **Provenance JSON** beside it: granule ids + CMR `updated` stamps / GDEX mtimes, polygon
  layer file + md5, CRS, weighting method (coverage-weighted mean, fill masked), software
  versions, machine, wall-clock; plus a per-basin QA table (pixel count, coverage
  fraction, `low_pixel_support`) mirroring the NLCD product.
- **Validation before use** (both already tooled):
  1. Structural checks as in `convert_daymet_csvs_to_parquet.py` (365 rows per site-year,
     no duplicate (site, Date), site set identical across years).
  2. Replay: run the WY 1993–2025 standard config with the new climate input and
     `docs/benchmarks/check_additivity.jl` against product #1. Expect the 6,087-basin
     climate signatures to reproduce for CY ≤ 2023 **only if** the polygons match the
     co-authors'; otherwise bounded drift that the T1 report must explain. New gages
     populate; WY 2024–2025 climate years appear.
- This reprocess is the "next rerun of any portion of the data" that regenerates
  `flagged_for_high_na` (CHANGELOG → Known Issues) — plan the product regeneration and the
  HydroShare re-staging together.

---

## 6. Decisions to take (user / co-authors)

1. **Variables**: prcp + swe now, the other four later (recommended), or all six in one
   pass (needs B, or a week of transfers and a 4 TB SSD).
2. **Size cap**: compute for all 7,964 polygons (cheap) and (a) keep the 85,000 km² rule
   as a documented downstream filter, or (b) drop it — the boundary layer already stops at
   100,000 km² and only 6 basins sit between the thresholds. The manuscript's HUC8
   justification argues for (a); a uniform product argues for (b).
3. **Polygons**: the rebuilt 2026-08-25 layer (Resource 4, the published one) —
   recommended — versus whatever the co-authors used (unknown until T4).
4. **Schedule vs. the Nov 9 submission**: a reprocess implies regenerating both standard
   products (~30 min each in Julia), re-staging Resources 1–3, and rewriting §2.1.3 / §3
   and three READMEs. Feasible in October if T1 passes in the first week; otherwise ship
   v1 as is and reprocess for a v1.1.
5. **Earth Engine**: whether an academic co-author hosts the project, or skip GEE.

---

## 7. Risks and gotchas

- **Transfer**: 0.57–3.4 TB; Globus from GDEX is the robust path; ORNL HTTPS has no
  published rate limit but unknown throughput; never on the exFAT drive.
- **GDEX R1-ness** of 2020/2021 unverified (attributes say 4.0 either way) — size/md5
  check in T1, or take those two years from ORNL.
- **S3 credentials expire hourly**; the Kerchunk tutorial itself shows the failure mode.
- **Memory on the 16 GB laptop**: never load a full NA year (92 GB uncompressed); read
  day slices or basin windows; the coverage-matrix path bounds memory by the union pixel
  count (tens of millions), not by the grid.
- **Polygon version drift** between our layer and the co-authors' run — the most likely
  cause of validation residuals; also the 9 mis-padded ids in the boundary layer (join on
  stripped ids).
- **2024–2025 Canadian station seam** (§2.1) — annotate in run notes, Resource 3 README,
  and the manuscript if 2024–2025 Canadian trends are shown.
- **Leap-year calendar** — keep; the preprocessor handles it; do not "fix" Dec 31.
- **Units**: prcp mm/day, swe kg/m² (= mm), tmin/tmax °C, vp Pa, srad W/m² (daylight
  average) — record in the dictionary.
- **Manuscript**: §2.1.3, the abstract's "6,041", §3 Resource 3 and Table 1 all change
  if this lands before submission.

---

## 8. Research trail

Four sub-agent reports (official routes; cloud mirrors; libraries and workflows; existing
products), 2026-09-29. Key URLs:
- Earthdata catalog record https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-daymet-daily-v4r1-2129-4.1;
  CMR granules `collection_concept_id=C2532426483-ORNL_CLOUD`; S3 credentials README
  https://data.ornldaac.earthdata.nasa.gov/s3credentialsREADME; AWS open-data registry
  https://registry.opendata.aws/nasa-daymet/; user guide PDF (release table, Canadian
  station note, calendar) https://data.ornldaac.earthdata.nasa.gov/public/daymet/Daymet_Daily_V4R1/comp/Daymet_Daily_V4R1.pdf.
- THREDDS retirement https://forum.earthdata.nasa.gov/viewtopic.php?t=7091 and t=7585;
  2025 release t=7638; 2024 release t=7277; Harmony warning t=7158; cloud OPeNDAP
  https://www.opendap.org/accessing-daymet-data-from-nasas-archives/.
- NCAR GDEX d682806 https://gdex.ucar.edu/datasets/d682806/ (data access
  https://gdex.ucar.edu/datasets/d682806/dataaccess/; THREDDS
  https://tds.gdex.ucar.edu/thredds/catalog/files/d682806/catalog.html).
- Google Earth Engine catalog https://developers.google.com/earth-engine/datasets/catalog/NASA_ORNL_DAYMET_V4;
  terms https://earthengine.google.com/terms/; noncommercial eligibility
  https://earthengine.google.com/noncommercial/; commercial transition
  https://developers.google.com/earth-engine/guides/transition_to_commercial; quotas
  https://developers.google.com/earth-engine/guides/noncommercial_tiers; pricing
  https://cloud.google.com/earth-engine/pricing; geemap discussion #1188.
- Planetary Computer STAC https://planetarycomputer.microsoft.com/api/stac/v1/collections/daymet-daily-na.
- pydaymet #72, #73; climateR #113; climateR-catalogs September-2026 release; gdptools
  HISTORY https://code.usgs.gov/wma/nhgf/toolsteam/gdptools/-/raw/develop/HISTORY.md;
  exactextract #183; CAMELS-SPAT https://github.com/CH-Earth/camels_spat; EASYMORE
  https://github.com/ShervanGharari/EASYMORE; NHM spatial targets
  https://github.com/rmcd-mscb/nhf-spatial-targets.
- Products checked: HYSETS OSF `rpc3w` (2023 update notes), CAMELS-SPAT HESS 2025,
  CAMELS DASH, Caravan Zenodo 14673536, CAMELSH Zenodo 16729675, MacroSheds, ScienceBase
  NHM-PRMS items, Borealis 10.5683/SP3/65FXAS, MACH Zenodo 18686475, BASINGRID Zenodo
  19358585.

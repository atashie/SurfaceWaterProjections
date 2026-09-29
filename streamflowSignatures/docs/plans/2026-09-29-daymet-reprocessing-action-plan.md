# Daymet reprocessing — action plan for the dedicated machine

**Date**: 2026-09-29 · **Status**: single-year gate (calendar 2023) PASSED on the M5 MacBook
the same day — read §0 first; it supersedes §§1–3 and parts of §§5 and 9 where they differ.
§§1–10 below are the original plan (numbers marked **M** measured, **E** estimated).

**Background and option comparison**: `2026-09-29-daymet-reprocessing-options.md`
(same folder). This document is the "how", that one is the "why".

---

## 0. Single-year test, 2026-09-29 — results and what changes

**User decisions (2026-09-29).** Reprocess every watershed that has a polygon (the 7,964-basin
boundary layer), entirely from scratch — no append to the stale series (option E is out).
Check the new values against the stale co-author product. Before any multi-year run,
process ONE year, time it, and compare it with the stale product.

**What ran.** Calendar 2023 (last year of the stale product; same chunk layout as 2020–2025),
all six variables, NCAR mirror → `~/Downloads/daymet_raw/`, processed on the 16 GB M5 MacBook
(10 cores); run folder `/Volumes/Untitled/daymet_processed_sep2026/` (every file md5-verified
by re-reading from the drive; staging copy `~/Downloads/daymet_work/`; its `RUN_NOTES.md` has
the full tables). Tools (written in `docs/benchmarks/daymet/`, moved the same day to
`EO_data_processing/daymet/`): `daymet_probe`, `_weights`, `_aggregate`, `_crosscheck`,
`_pixelcheck`, `_validate`, `_assemble`, `_stream` (the §5.5 orchestrator). An extra
old-layout + leap-year probe ran on prcp 1980.

**Measured** (replaces the **E** cells of §3):

| Step | 2023 |
|---|---|
| Download, six files, 72.2 GB | 66 min (3,972 s); 15.0–20.9 MB/s single stream, 18.2 MB/s overall (4 parallel streams ≈ the same total — the link is the cap) |
| Checksums | SHA-256 of all six = NASA CMR's → mirror ≡ ORNL byte for byte. CMR publishes a SHA-256 for every granule; all 270 mirror files 1980–2024 have exactly the CMR size (2025 is ORNL-only) |
| Weights (once) | 25.6 M (basin, cell) weights over 9.4 M distinct cells: 4 s on the published layer, 80 s on full-resolution polygons |
| Aggregation, one variable-year, 8 workers, 2020+ layout (10,300,300) | 13–21 s wall (103 s for all six); 214 of 729 tiles hold basins (29 %) → 28 GB decompressed; peak RSS 3.7 GB parent + ≤ 0.6 GB per worker |
| Same, 1980–2019 layout (1,1000,1000) — prcp 1980 | 32 s; 29 of 72 tiles (40 %) → 41.8 GB decompressed; peak RSS 3.6 GB parent + 0.75 GB per worker (≈ 9.6 GB with 8 workers — use `--workers 6` if the laptop is also in use) |
| Output | 8–24 MB zstd parquet per variable-year (7,964 × 365 rows); 110 MB per assembled year → ≈ 5 GB for 1980–2025 |

→ Compute ≈ 1.5 min per year of six variables. The job is purely download-bound:
1980–2025 = 3.455 TB ≈ **50 h** at 19 MB/s (prcp + swe 0.57 TB ≈ 8.5 h). The "dedicated
machine" is unnecessary: this laptop runs it with ≈ 65 GB of free disk (download → verify →
aggregate → delete, one file ahead).

**Correctness, independent of the stale product.** (1) exactextract's own `weighted_mean`
on five days per variable, all 7,964 basins, matches the chunk-aligned aggregation to
≤ 5e-11 absolute, ≤ 1e-13 relative (`daymet_crosscheck.py`). (2) ORNL's Single Pixel API
at six points (TN, WA, CO, QC, FL, YT), all six 2023 variables and prcp 1980: our LCC x/y equal ORNL's to 0.00 m and the
pixel series match to float32 rounding while the neighbouring pixel differs by up to 18 mm
(`daymet_pixelcheck.py`). (3) Output is bit-identical across reruns and worker counts
(partial sums accumulated in a fixed order).

**Agreement with the stale product** (5,969 shared basins; table in `RUN_NOTES.md`) forced
three method findings:
1. **Weights must be TRUE AREAS, not coverage fractions.** The co-authors' gdptools values
   are an area-weighted mean. Coverage fractions of the conformal LCC grid (areal scale
   0.91–1.21 over the basins) over-weight cells where k < 1 and move annual totals of large
   northern basins by up to 0.3 %. With coverage × true cell area and full-resolution
   polygons, all six 2023 variables reproduce the stale values: R² ≥ 0.99999996 in every
   basin (swe: 0.99979 in a basin whose mean SWE is 4e-8 mm); prcp annual totals within
   ±0.0011 % (1st–99th percentile).
2. **The published Resource 4 layer is 200 m-simplified (7,823 of 7,964 basins), and that is
   the whole remaining difference.** Published layer: prcp annual totals within ±0.1 %
   (p01–p99), up to 1.2 % in small high-relief basins; swe annual sums within ±2.5 %
   (p01–p99), far more in trace-snow basins. A full-resolution rebuild — the same rebuild
   script with the simplification step disabled (`polygons_fullres/` in the run folder,
   1.3 GB) — removes it.
3. **Fill handling differs.** The stale product is NaN for every basin that touches even one
   fill cell: 4 St. Johns River basins in Florida (02234500, 02236000, 02236125, 02244040;
   99.95–99.98 % valid) are NaN in all 44 years, as is 01372058 (not a usable gage, no
   polygon). The four sit in the delivered products (4 in #1, 3 in #2) with every climate
   signature NaN, so the "signatures + climate" counts 5,517 / 5,638 overstate by 4 / 3. The
   new masked mean renormalises over valid cells, so they gain climate.

**Coverage.** Every one of the 5,965 usable gages that had Daymet keeps it (all have
polygons). Product #1 goes from 5,513 gages with a non-NaN climate series to 6,634
(+1,121), product #2 from 5,635 to 6,205 (+570). Still without climate: the 54 usable basins > 100,000 km² that
the boundary layer excludes (44 / 45 product gages); their source polygons exist.

**Unknowns of §9 now.** U1 closed (sizes above; SHA-256 per file at download). U2
15–21 MB/s here. U3 absent in prcp 1980 (all 26,280 chunks stored); `daymet_stream.py`
probes every file and stops if an unstored chunk lies under a basin. U4 29 % (2020+ layout),
40 % (1980–2019). U5 validity constant in time for every basin in 2023. U6 answered
(finding 3). U7 answered (finding 2). U8 measured. U9 no date shift (lag-0 r = 0.99999997 vs
0.15 at ±1 day); leap years verified — 1980 has Feb 29 and no Dec 31, and prcp 1980
reproduces the stale values (annual totals within ±0.001 %, p01–p99). **U10 open: no
Earthdata Login on this Mac; the six 2025 files need one.** U12 measured.

**DECISIONS (user, 2026-09-29, after the gate):**
- **D1 polygons: FULL RESOLUTION.** The Resource 3/4 READMEs will say the climate uses the
  unsimplified source polygons, of which Resource 4 is a 200 m simplification. The
  full-resolution layer comes from the now-committed
  `EO_data_processing/geometry/rebuild_watershed_polygons.py --no-simplify --include-large`
  (its defaults reproduce Resource 4: same 7,964 basins, ids, order and vertices; areas to
  ~1e-11, a PROJ-build difference).
- **D2 variables: ALL SIX in one pass** (≈ 50 h of downloads).
- **D3 Earthdata Login: the user supplied an EDL bearer token** (account arik.tash@gmail.com;
  the user states access is free and not rate-limited). The token expires ≈ 2026-10-25, so
  the run must finish before then; fetch the six 2025 files first. Store it at run time
  outside the repo (e.g. `~/.config/earthdata/edl_token`, mode 600) and never in docs,
  memory or git. DONE: `daymet_stream.py` sends the token through a mode-600 header file
  (never on a command line; curl drops it on the redirect to ORNL's CloudFront store) and
  downloads the ORNL-only 2025 files first; a 1 MB range request with it returned HTTP 206.
- **D4 large basins: INCLUDE the 53 correctly delineated basins > 100,000 km² if RAM
  allows.** 05KH009 stays out: its HydroBASINS fallback drew a 328,000 km² river as 200 km².
  The 53 (all Canadian, heavily nested; Mackenzie 10LC014 = 1.68 M km²) sum to 16.6 M km²
  of reported area against 26.1 M km² for the 7,964. MEASURED → INCLUDED: the 8,017-basin
  weights hold 45.8 M entries over 11.2 M cells (build 118 s, peak ≈ 8 GB once); after a
  planning fix (int32 keys, slices; output byte-identical) aggregation takes 33 s on the
  1980–2019 layout with `--workers 6`, peak 4.3 GB parent + < 1 GB per worker, and touches
  no new tiles. These basins have no stale series to validate against and no polygon in
  Resource 4, so the deposit documents them separately.
- **LOCATIONS (user, 2026-09-29):** code in `EO_data_processing/` (`daymet/`, `viz/`,
  `geometry/`); the run folder on the internal SSD, `~/HISSS_data/daymet-processed-29sep2026/`,
  then md5-verified copies to `/Volumes/Untitled/daymet-processed-29sep2026/`; commit and
  merge when ready.

Then: `caffeinate -i python daymet_stream.py --years 1980-2025 --workers 6 …` (keeps the
laptop awake for the ~2-day run) → `daymet_assemble.py` → `daymet_validate.py --years
1980-2023` → verified copy to the drive → the Phase 1 replay (§6) against product #1.

---

## 1. Objective and defaults

Produce a new basin-averaged Daymet daily climate input for the signature pipeline:

| | Current input | Target |
|---|---|---|
| Basins | 6,087 (5,965 of the 8,014 usable gages) | **all 7,964** watersheds of the delivered boundary layer (Resource 4) |
| Years | calendar 1980–2023 | **calendar 1980–2025** (= WY 1980–2025; 2026 does not exist yet) |
| Variables | prcp, tmin, tmax, swe, vp, srad | same six; `dayl` optional (off by default) |
| Source | Daymet V4 (R1) via co-authors' gdptools | Daymet V4 R1 annual NA mosaics: **NCAR GDEX mirror for 1980–2024**, ORNL for 2025 |
| Method | gdptools area-weighted mean | **coverage-weighted mean** (exactextract weights once, chunk-aligned accumulation) |
| Output | parquet `site_id, Date, prcp, tmin, tmax, swe, vp, srad` | same schema + provenance + per-basin QA |

**Defaults chosen for this plan** (change them consciously, they are not user decisions
yet — see options doc §6): process variables in the order **prcp, swe** (tranche 1,
usable on its own), then **tmin, tmax, vp, srad** (tranche 2); compute for all 7,964
polygons regardless of size (the 85,000 km² rule becomes a downstream filter if the
co-authors want to keep it); use the published Resource 4 polygons; keep Daymet's
365-day calendar exactly.

---

## 2. Measured facts the plan relies on (2026-09-29)

All measured from this Windows laptop against the NCAR mirror
(`https://tds.gdex.ucar.edu/thredds/fileServer/files/d682806/<file>`), which serves the
files anonymously with `Accept-Ranges: bytes`. Probe scripts:
`docs/benchmarks/daymet/probe_*.py` are NOT written yet — the one-off scripts used lived
in the session scratchpad; §5 step 0 re-creates the probe as a permanent tool.

| Fact | Value | How |
|---|---|---|
| File naming | `daymet_v4_daily_na_<var>_<yyyy>.nc`, one per variable-year; mirror holds 315 NA files = 7 vars × 1980–2024; ORNL holds 1980–2025 | mirror catalog; CMR (`C2532426483-ORNL_CLOUD`) **M** |
| Grid | float32 (time 365, y 8,075, x 7,814) = 63,098,050 cells; 1 km; `x` −4,560,250 … 3,252,750 m, `y` 4,984,000 … −3,090,000 m (cell centres, y descending) | file header **M** |
| CRS (`lambert_conformal_conic`) | central meridian −100, latitude of origin 42.5, standard parallels 25 / 60, false E/N 0, semi-major 6,378,137, inverse flattening 298.257223563 (WGS84) | file header **M** |
| Time | `days since 1950-01-01 00:00:00`, calendar `standard`, values at noon (e.g. 27028.5 = 2024-01-01), `yearday` 1…365, "24-hour day based on local time"; **365 values every year, leap years drop Dec 31** | file header **M**; calendar rule from the user guide |
| Fill | `_FillValue = missing_value = −9999` on every variable; HDF5 dataset fill value is **0.0 in 1980–2019 files** and −9999 in 2020+ files (matters only for unstored chunks — §9 U3) | file header **M** |
| Chunking 1980–2019 | (1, 1000, 1000), gzip level 4, no shuffle → 365 × 9 × 8 = 26,280 chunks, 4 MB uncompressed each | **M** on prcp 1980/1990/2000/2010/2018/2019 |
| Chunking 2020–2025 | (10, 300, 300), gzip 4, shuffle → 37 × 27 × 27 = 26,973 chunks, 3.6 MB each; written with netCDF 4.9.2 / HDF5 1.14.3 (the R1 remake) | **M** on prcp 2020–2024, swe 2024, vp 2024 |
| Compressed sizes 1980–2024 (mirror) | prcp 0.174 TB · swe 0.390 · tmin 0.521 · tmax 0.533 · vp 0.938 · srad 0.825 · dayl 0.059; **six vars 3.38 TB, prcp+swe 0.564 TB; 72–77 GB per year (six), ≈ 12.5 GB (prcp+swe)** | catalog sum **M** |
| 2025 sizes (ORNL CMR, MB) | prcp 4,360 · swe 7,646 · tmin 14,019 · tmax 13,735 · vp 17,023 · srad 15,852 · dayl 1,484 | CMR **M** |
| Uncompressed per variable-year | 92 GB | arithmetic |
| Download, mirror → this laptop | 18.3 MB/s single stream; ≈ 25 MB/s aggregate with 4 streams | curl 300 MB ranges **M** |
| Decompression on this laptop | 150–170 MB/s of uncompressed output per core (Python zlib + unshuffle) | 9 land chunks each of prcp/vp 2024 **M** |
| Basin pixel budget | 7,889 basins ≤ 100,000 km² in the metadata sum to 25.6 M km² ≈ 25.6 M weight entries; 662 basins > 10,000 km², 36 > 50,000, largest 93,900 | metadata CSV **M** |
| Sanity of the mirror | 2020/2021 files carry the 2020+ layout while 2019 carries the old one → consistent with R1 remakes; byte-identity with ORNL **not** verified | **M** layout, identity **U** |

---

## 3. Resource budget (per year, six variables) — E unless noted

| Resource | Budget | Basis |
|---|---|---|
| Download per year | 50–70 min at 18–25 MB/s (**M** rate on this laptop's link; the dedicated machine's link is **U**) | 72–77 GB **M** |
| Compute per year | 6–15 min on 4 workers; **plan 30 min** | decompression rate **M** × fraction of tiles touching basins (**E** 30–50 %) |
| Wall-clock, pipelined | **≈ 2–2.5 days for six variables; ≈ 9–10 h for prcp+swe** | download-bound |
| RAM peak | **2–3 GB with 4 workers, 4–5 GB with 8** (weights ≈ 0.2–0.3 GB, chunks 4 MB each in flight, accumulators 23 MB per variable, ≈ 0.3 GB Python baseline per process) | arithmetic, not measured end-to-end |
| Disk, pipelined | ≈ 155 GB free (two years of raw input) + 3 GB outputs + 0.3 GB weights; ≈ 25 GB for prcp+swe; ≈ 45 GB if variable-at-a-time | sizes **M** |
| Output size | 7,964 × 365 × 6 ≈ 17.4 M values/yr → ≈ 50–70 MB parquet/yr, ≈ 3 GB total | arithmetic |
| Cores | 4 is enough (download-bound); 8 shortens the compute phase, irrelevant when pipelined | — |
| Network total | 3.4 TB (six) / 0.58 TB (prcp+swe) | **M** |

The exFAT thumb drive must not be used for any of this (three silent truncations on
record). Use an internal or USB-C SSD.

---

## 4. Inputs to bring to the machine

1. **Boundary layer** — the published Resource 4 file (`hisss_watershed_boundaries.gpkg` /
   `.parquet`, 7,964 rows, EPSG:4326, columns incl. `gage_id`, `canon_id`,
   `watershed_geom_source`, `geom_area_km2`, `low_confidence`). Sources: the HydroShare
   staging folder on the Mac (`~/Downloads/Signatures/resource4_*`) or the Google Drive
   backup. Record its md5.
2. **Old climate parquet** `daymet_1980_2023_rebuilt_10aug2026.parquet` (3.76 GiB; on the
   thumb drive `D:/processedOuts_feb2026/` — copy it, verify size 4,040,997,608 B and the
   `PAR1` footer). Needed only for validation (§6).
3. **Metadata CSV** `golden-outputs/combined_watershed_metadata_09feb2026.csv` (in the
   repo) for ids, areas, country.
4. **Product #1 CSV + config** (`processedOuts_drought_28jul2026`) for the replay in §6
   step 4 — optional but the strongest test.
5. **Earthdata Login** (free; needed for the 2025 files and the identity check). Put it in
   `~/.netrc` (`machine urs.earthdata.nasa.gov login … password …`, mode 600).
6. Optional **Globus Connect Personal** if the mirror's HTTPS is slow from that machine.

Software (Python ≥ 3.11): `numpy scipy h5py fsspec aiohttp requests pyarrow pandas
duckdb geopandas shapely pyproj rasterio exactextract>=0.3 tqdm` (+ `earthaccess` for the
ORNL files, `xarray` optional). Pin the versions you install in the run's provenance.
The Julia environment of this repo is needed only for the §6 replay.

---

## 5. Tools to write (all under `docs/benchmarks/daymet/`, outputs in the run folder)

Repo convention: tools live in `docs/benchmarks/`, every run artifact lives in ONE
experiment folder on the data drive, e.g. `daymet_1980_2025_<ddmonyyyy>/`
(CLAUDE.md Critical Constraint #5). Each tool below is small; the contracts matter more
than the code.

### 5.0 `daymet_probe.py` — re-verify the facts of §2 (½ h)
Opens a file remotely via `fsspec` + `h5py` (byte ranges) and prints shape, chunks,
compression, fill, CRS attrs, time attrs, and `get_num_chunks()` for one old and one new
file. Run it first on the dedicated machine; if anything differs from §2, stop and
update this plan.

### 5.1 `daymet_manifest.py` — the file list with provenance (1 h)
- For 1980–2024: parse the mirror catalog (`…/thredds/catalog/files/d682806/catalog.html`)
  → url, size, last-modified for each `daymet_v4_daily_na_<var>_<yyyy>.nc`.
- For 1980–2025: query CMR (`collection_concept_id=C2532426483-ORNL_CLOUD`, page through
  `granules.json`) → granule id, size, checksum if present, `updated`.
- Cross-check mirror size == CMR size for every 1980–2024 file; list mismatches
  (expect none — **U**). Write `manifest.json`. Decide per file: source = mirror if sizes
  match, else ORNL.

### 5.2 `daymet_fetch.py` — resumable download of one file (2 h)
`--file <name> --source mirror|ornl --dest <dir>`: HTTP range-resume into `<name>.part`,
verify final size against the manifest (and md5 when CMR provides one), atomic rename,
log MB/s. ORNL needs the Earthdata cookie/redirect dance — simplest is
`earthaccess.login()` + `earthaccess.download()`, or `curl -n -b -c -L`. Retries with
backoff. Never delete anything here.

### 5.3 `daymet_weights.py` — coverage matrix, once (3 h)
1. Read polygons; drop nothing; reproject to the Daymet CRS built from the header attrs
   (§2; `pyproj.CRS.from_cf(...)` or the proj string
   `+proj=lcc +lon_0=-100 +lat_0=42.5 +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +ellps=WGS84 +units=m`).
2. Build the reference grid: transform from `x[0]−500, y[0]+500` with 1,000 m pixels,
   7,814 × 8,075 (assert against the file's x/y vectors).
3. `exact_extract(reference_raster, polygons, ["cell_id", "coverage"])` per polygon (a
   dummy raster with the grid geometry suffices — exactextract only needs the grid to
   compute coverage) → for each polygon, the flat cell ids and coverage fractions.
4. Assemble CSR `W` (7,964 × 63,098,050; nnz ≈ 25.6 M): row = basin, column = flat cell
   index, value = coverage fraction (NOT yet normalised — normalisation happens per day
   against valid pixels, §5.4). Save `weights.npz` + `basin_qa.csv` (`gage_id`,
   `n_cells`, `coverage_sum`, `n_cells_lt_half`, `low_pixel_support` = coverage_sum < 4).
5. Tile index: for each spatial tile of BOTH layouts (1000×1000 and 300×300 tilings),
   the list of basins with any cell in the tile and the sub-CSR for that tile. Save.
   Report the fraction of tiles touched by any basin (this is estimate E in §3 — record
   the measured value).
- Check: Σ coverage per basin ≈ `geom_area_km2` (1 km² pixels) within a few %.

### 5.4 `daymet_aggregate.py` — one variable-year → per-basin daily means (1 day incl. tests)
Contract: `--file <path> --weights <dir> --out <parquet>`; deterministic; resumable at
the file level (writes `<out>.done` on success).

Core loop (chunk-aligned; works for both layouts):
```
open h5py dataset d; ct, cy, cx = d.chunks; ndays = 365
num = zeros(n_basins, ndays); den = zeros(n_basins, ndays)
for each spatial tile (iy, ix) that touches a basin:            # from the tile index
    Wt = tile sub-CSR (rows = basins touching, cols = cells in tile order)
    for each time chunk t0 in range(0, ndays, ct):
        block = d[t0:t0+ct, iy*cy:(iy+1)*cy, ix*cx:(ix+1)*cx]   # exactly one HDF5 chunk
        vals = block.reshape(ct_actual, -1).T                   # cells × days
        valid = vals != -9999 (and not NaN)
        num[rows, t0:t0+ct] += Wt @ where(valid, vals, 0)
        den[rows, t0:t0+ct] += Wt @ valid.astype(float32)
mean = num / den   (NaN where den == 0)
```
- Parallelism: multiprocessing over tiles (each worker opens its own `h5py.File`; on
  Windows use `spawn`). 4–8 workers.
- The old layout stores one day per chunk: the same loop with `ct = 1` (365 reads per
  tile); the new layout gives 10 days per read. Both are one full pass over the touched
  tiles.
- Output parquet: `site_id` (string, zero-padded as in the boundary layer), `Date`
  (from `time` → `1950-01-01 + floor(days)`), `<var>` float32; 365 rows per basin;
  plus `<var>_den` (sum of valid weights) kept in a sidecar for QA, not in the product.
- Tests (synthetic file, 20 basins, both chunk layouts): constant field → mean equals
  constant; a basin straddling a tile boundary; a basin over fill; a chunk with mixed
  fill; results identical with 1 vs 8 workers.

### 5.5 `daymet_stream.py` — the orchestrator (½ day)
For `year in 1980..2025`, for `var in [prcp, swe, tmin, tmax, vp, srad]`:
fetch (unless present and size-verified) → aggregate → verify the output (365 rows ×
7,964 basins, no NaN except basins with `den == 0`, value ranges: prcp ≥ 0, swe ≥ 0,
−70 < t < 60, vp > 0, srad ≥ 0) → write `.done` → **delete the raw file** → append a
line to `stream_log.csv` (file, source, bytes, MB/s, seconds to aggregate, peak RSS via
`psutil`). Prefetch the next file in a background thread while aggregating (bounded to
one file ahead; that is the 2-year disk budget). Fully resumable: a variable-year with a
`.done` is skipped; a `.part` is resumed.
Order of the loops: variable-major within tranche 1 (`prcp` 1980…2025, then `swe`
1980…2025) so that tranche 1 is complete and usable as early as possible; then the four
others.

### 5.6 `daymet_assemble.py` — the deliverable (½ day)
Join the per-variable-year parquets into `daymet_1980_2025_<date>.parquet`
(`site_id, Date, prcp, tmin, tmax, swe, vp, srad`; one row per site-day; **365 rows per
site-year**; sorted by site, Date). Structural checks copied from
`docs/benchmarks/convert_daymet_csvs_to_parquet.py` (per-site-year row counts, no
duplicate keys, site set identical across years). Write `provenance.json` (manifest
subset actually used, polygon file + md5, weights md5, software versions, machine,
per-file MB/s and seconds, git commit of the tools) and `daymet_basin_qa.csv`. A
prcp+swe-only assembly after tranche 1 is allowed and expected.

### 5.7 `daymet_validate.py` — against the old input (½ day)
For the 6,087 shared basins (join on zero-stripped ids), CY 1980–2023, per variable:
per-basin identity R² and max |Δ| of daily values, annual totals ratio; histogram of R²;
list of basins with R² < 0.999 joined to `basin_qa.csv` (pixel count,
`watershed_geom_source`, coastal?) so residuals can be explained. Plus five centroids
checked against the Single Pixel API (`https://daymet.ornl.gov/single-pixel/api/data`).

---

## 6. Runbook (phases, checkpoints, go/no-go)

**Phase 0 — setup and one-year proof (½–1 day, before any long run)**
1. Bring inputs (§4); create the run folder; `pip install`; record versions.
2. `daymet_probe.py` on prcp 2000 and prcp/swe 2024 → §2 confirmed (else stop).
3. `daymet_manifest.py` → mismatches = 0 (else investigate; fall back to ORNL for
   mismatching files).
4. `daymet_weights.py` → `weights.npz`; Σcoverage vs area check; touched-tile fraction
   recorded.
5. Fetch prcp 2023 + swe 2023 from the mirror (record MB/s — **this is the number that
   sets the schedule**) and prcp 2025 from ORNL (record MB/s and that EDL works).
6. `daymet_aggregate.py` on the three files with 4 workers; record wall-clock and peak
   RSS.
7. `daymet_validate.py` for 2023 against the old parquet.
- **Go** if: ≥ 95 % of shared basins have daily-prcp R² ≥ 0.999 with the rest explained
  (small/coastal/different polygon), aggregate ≤ 30 min per variable-year, download rate
  makes tranche 1 ≤ 1 day. **No-go**: agreement fails broadly → the polygons or the
  weighting differ from the co-authors' run; stop and consult (see options doc T4 —
  ask the USGS co-author for their script/polygons) before spending days of transfer.

**Phase 1 — tranche 1 (prcp, swe, 1980–2025; ≈ 10 h unattended)**
`daymet_stream.py --vars prcp,swe`. Then `daymet_assemble.py` (prcp+swe file),
`daymet_validate.py` (all years), and the **replay**: run the WY 1993–2025 standard
config with `STREAMFLOW_CLIMATE_PATH` pointing at the new file and
`docs/benchmarks/check_additivity.jl` against product #1. Expected: climate/snow
signatures reproduce for the 6,087 basins where the daily series agree; new basins
populate; WY 2024–2025 climate years appear. Ship this file to the signature work even if
tranche 2 is still running.

**Phase 2 — tranche 2 (tmin, tmax, vp, srad; ≈ 2 days unattended)**
`daymet_stream.py --vars tmin,tmax,vp,srad`, then the full assembly + validation. `dayl`
only if asked (1.3 GB/yr).

**Phase 3 — downstream (separate sessions)**
Regenerate both standard products with the new input (this is the rerun that also
regenerates `flagged_for_high_na`); update the HydroShare Resource 3 Daymet table and
READMEs (site count, years, method, the Canadian station seam), the dictionary, the
manuscript §2.1.3 / §3 / abstract counts, `docs/DATA_SOURCES.md`, DEVELOPMENT.md's
Active Parquet Files, and the claude-skill. Log everything in CHANGELOG.

---

## 7. Acceptance criteria for the new input

1. Structure: 7,964 sites × 46 years × 365 rows, no duplicate (site, Date), the site set
   identical across years, Dec 31 absent in leap years.
2. Coverage: every basin has `den > 0` on every day for prcp (fill-only basins would be
   a polygon error — expect zero of them); `basin_qa.csv` lists any basin with
   coverage_sum < 4 pixels.
3. Agreement with the old input on the 6,087 shared basins: report the R² distribution;
   the acceptance threshold is the Phase 0 one (≥ 95 % at ≥ 0.999 for daily prcp) with a
   written explanation of the rest. This is agreement with the co-authors' *run*, not a
   proof of correctness — the synthetic tests in §5.4 are the correctness check.
4. Replay: `check_additivity.jl` reports 0 columns added/dropped and an identical gage
   set for the shared basins' streamflow-only columns (they must not move at all), and the
   climate-dependent columns move only where the daily series differ.
5. Provenance JSON complete (manifest, md5s, versions, timings).

---

## 8. Deliverables and where they go

| Artifact | Location |
|---|---|
| Tools | `docs/benchmarks/daymet/*.py` (+ tests) — committed |
| Run folder | `daymet_1980_2025_<date>/` on the data drive: `manifest.json`, `weights.npz`, `basin_qa.csv`, per-variable-year parquets, `stream_log.csv`, `provenance.json`, `validation_*.csv/md`, the assembled parquet(s) — NOT in the repo |
| The input for the pipeline | `daymet_1980_2025_<date>.parquet` (and the prcp+swe interim file) — referenced by `STREAMFLOW_CLIMATE_PATH`; add to DEVELOPMENT.md → Active Parquet Files with size and footer check |
| Docs | this plan updated with the measured numbers (fill the **E** cells), CHANGELOG entry, DATA_SOURCES.md row 4 |

---

## 9. Unknowns — honest list

- **U1 Mirror vs ORNL byte-identity.** The mirror's 2020/2021 files have the R1 layout
  (good sign) but no checksum comparison has been made. Step 5.1 compares sizes; if CMR
  exposes md5, compare it for at least 2020, 2021 and one old year. If it doesn't, take
  2020–2021 from ORNL directly (≈ 25 GB for six variables) to be safe.
- **U2 Download rate on the dedicated machine** — 18–25 MB/s was this laptop on this
  network; Globus and a different link are untested. The whole schedule scales with it.
  Also observed: the mirror dropped a long-lived HTTP connection once during the
  probing (connection reset after several minutes of range reads), so §5.2's
  range-resume is not optional — expect a few resets per terabyte.
- **U3 Unstored chunks in the old files.** The 1980–2019 files have an HDF5 fill value of
  0.0 while `_FillValue` is −9999. If any all-ocean chunk is *omitted* rather than
  written, h5py returns 0.0 for it — a valid-looking zero, not fill — and a coastal
  basin touching that chunk would be biased. A `get_num_chunks()` count against the
  26,280 possible chunks was attempted remotely and FAILED after ~5 min — the mirror
  reset the connection mid-transfer (`ConnectionResetError 10054`, "received 1,351,328
  of 8,391,744 bytes"). **Do this check locally in Phase 0**; if chunks are
  missing, treat unstored chunks as fill (h5py `read_direct_chunk` raises for a missing
  chunk, or compare against the 2020+ files' explicit −9999 mask).
- **U4 Fraction of tiles touched by basins** (30–50 % assumed) — sets compute time; the
  tile index reports it.
- **U5 Fill mask stability** — whether the −9999 mask is identical across variables and
  days (ocean only) or whether swe/prcp carry additional fill. The per-day `den`
  handles either case; check `den` variability in QA.
- **U6 What the co-authors' gdptools run did at partial coverage and fill** —
  renormalisation over valid pixels is assumed here; differences would show up at
  coastal basins in validation.
- **U7 Polygon versions** — the co-authors' 6,087-basin polygons may differ from the
  published Resource 4 layer; residuals in validation may be polygon, not method.
- **U8 Peak RAM** is computed (§3), not measured end-to-end; measure in Phase 0 step 6
  and fill in.
- **U9 Time semantics** — Daymet days are "24-hour day based on local time"; the old
  input was rebuilt positionally from CSVs with the same convention, so no shift is
  expected; the 2023 validation will confirm (a one-day shift gives a very visible R²
  drop).
- **U10 ORNL HTTPS throughput and EDL from a script** — only a 1 MB range was planned;
  the 2025 files are ≈ 73 GB.
- **U11 Canadian 2024–2025 station seam** (ECCC/CCCS replaced GHCNd CA0 stations) — a
  data property to document, not something this run can fix.
- **U12 exactextract on a 63 M-cell reference grid** — building coverage for 7,964
  polygons is expected to take minutes, but the largest basins (≈ 94,000 cells) and the
  `cell_id` output size (25.6 M rows) have not been timed.

---

## 10. Do-not list

- Do not read a whole variable (`d[:]`, 92 GB) or a whole day-stack into memory; stay
  chunk-aligned.
- Do not resample the raster or reproject it; reproject the polygons.
- Do not "fix" the missing Dec 31 in leap years; the preprocessor expects the hole.
- Do not re-pad ids; join on the zero-stripped form (9 boundary ids are mis-padded).
- Do not put any of this on the exFAT thumb drive.
- Do not delete a raw file before its output passed the checks and `.done` exists.
- Do not run the streamflow products against the new input without the §6 replay.

# NLCD Annual — per-watershed land-cover processing plan

Per-watershed **Annual NLCD** (USGS/MRLC Collection 1) land cover + fractional impervious
surface, aggregated over each gage's upstream watershed into an annual per-class-% record
that sits alongside the streamflow signatures and the MODIS EO products, joining by
leading-zero-safe `gage_id` / `canon_id`.

This is a **CONUS-only, US** companion to the continental MODIS LULC product
(`README.md`). It is a **metadata / ingestion product** (like HydroATLAS and the MODIS EO
layers) — **Python**, not a cross-language signature, not ported to Julia/R.

> **✅ QA COMPLETE + STAGED FOR HYDROSHARE (25 Aug 2026)**: the finalized product was
> retrieved from the Drive backup, re-verified against every documented invariant
> (250,879 rows / 6,119 gages / C1V2 / class sums exact), and staged for deposit as
> HydroShare Resource 5 (`~/Downloads/Signatures/resource5_eo_products/`, deposit names
> `hisss_nlcd_annual.parquet` etc.). **The human-QA pass via `nlcd_explorer.html` is
> DONE — user signed off 25 Aug 2026** — closing the last publication gate. Note the
> original `nlcd_out_of_footprint.csv` was not in the backup; it was exactly
> reconstructed as `hisss_nlcd_excluded_gages.csv` (45 gages, all verified Alaska).
> See CHANGELOG → August 2026.

> **⚠ ACCESS CHANGE (24 Aug 2026)**: S3 access is LOST (see `README.md` and CHANGELOG →
> August 2026). The pending `nlcd_finalize.py --upload` S3 publish is defunct, and the
> `temp_lulc_conus/` staging (which lived in the lost bucket) is moot. The finalized
> product **is confirmed saved in the project Google Drive backup** (folder
> `1DVuq4nC5j_Y01sBaDj9cwjbv7S8sndjj`; user-verified 24 Aug 2026) — retrieve it from
> there. ~~The human-QA pass via `nlcd_explorer.html` is still owed before publication
> (now targeting HydroShare, not S3).~~ *(Done 25 Aug — see above.)*

> **Status (24 Jul 2026) — BUILD COMPLETE (local); pending human QA + S3 publish.** Extraction
> pipeline, finalize, and QA/QC explorer are all built and validated. Finalized CONUS product =
> **250,879 rows (6,119 gages × 41 yr, 1985–2025)**; 45 Alaska gages excluded. **Remaining:**
> human QA via `nlcd_explorer.html` → `nlcd_finalize.py --upload` (publish to S3) → delete the
> `temp_lulc_conus/` staging (~93 GB) → CHANGELOG. The dated history below traces how we got here.
>
> **Status (22 Jul 2026)**: PLAN — design-reviewed + **access/method pilot PASSED + micro-pilot
> PASSED + source mosaics staged**. Legend manifest (`nlcd_legends.csv`) built. Source data
> staged to `streamflow/temp_lulc_conus/` (82 mosaics, 93.4 GB, C1V2; delete after final QA/QC).
> `nlcd_pipeline.py` written + **Codex adversarial review (22 Jul) = FIX-FIRST → all blockers +
> majors + minors fixed** and re-smoke-tested clean (20 basins, sum100 dev ≤8e-13, imp_oob=0).
> Fixes: publish-integrity gates (per-year `_validate_year`; final table written ONLY if every
> requested year has a valid checkpoint, else exit nonzero — no silent partial publish), atomic
> checkpoint/final writes (tmp+os.replace), download validated vs S3 ContentLength, both rasters
> schema-checked, C1 window rounded+intersected, `valid_coverage_frac` hard-fails >1.1 instead of
> silent clip, final concat scoped to requested years, merges `validate="one_to_one"`, columns
> reindexed to a committed order, `assert`→explicit exceptions. **Next: launch the ~10 h full run
> (awaiting go-ahead) → `nlcd_finalize.py` + S3 → `build_nlcd_explorer.py` → delete temp_lulc_conus.**
> Derived product not on S3 yet.
>
> **Full run COMPLETE + validated (24 Jul):** `watershed_nlcd_annual.parquet` — 252,724 rows =
> 6,164 gages × 41 yr; class %s exactly 100 on covered rows; impervious 0–68.8%; 45 all-fill
> Alaska gages isolated. **Codex reviewed results + finalize/explorer plans = GO-with-fixes.**
> Reconciliation (24 Jul) folded into the finalize/explorer specs below:
> - **Never publish out-of-footprint zeros as real** → finalize DROPS the 45 Alaska gages
>   (identified by `max(n_nlcd_pixels)==0`) to `nlcd_out_of_footprint.csv` (reason "Alaska —
>   outside CONUS Annual NLCD footprint"; QA-confirmed 45/45 AK, 0 HI/PR) + an all-gage QA
>   companion; CONUS table = 6,119 × 41 = 250,879 rows.
> - **Separate, documented QA flags** (not a single opaque union): `geom_low_confidence`,
>   `low_pixel_support` (<100 eff px; QA-measured flags 0), `partial_coverage`
>   (`valid_coverage_frac`<0.99; flags 3), `low_confidence` = documented union. Coverage is
>   gage-level stable (0 gages vary across years) and near-complete.
> - **Join integrity**: gage_id str, assert geometry unique, canon_id equality check, merge
>   validated one-to-one; **provenance** metadata JSON + `nlcd_collection='C1V2'` col +
>   `valid_area_km2`; dictionary states units/valid-pixel-conditional %/impervious semantics.
> - **Explorer**: smoothed endpoint-difference deltas mean(1985–87) vs mean(2023–25) LABELLED as
>   endpoint-diff (not trend) with the annual-update seam disclosed; developed% and impervious%
>   shown separately (non-additive; QA-confirmed corr 0.963); a largest-one-year-change / shrub↔
>   grass artifact view (QA-confirmed −0.71 swap correlation, 77–90 pp swings); partial-coverage/
>   low-support badges; excluded-45 disclosed; normalized Shannon; compact 0.1% encoding.
> - **HOLD S3 publish until human QA via the explorer** (Codex advice); build both locally first.
>
> **Micro-pilot (22 Jul, LOCAL 2023 mosaics, ALL 6,164 US `gage_type=='USGS'` basins):**
> Per-year both-layer extraction = **58.8 min** (land cover 32.4 min @315 ms/basin; impervious
> 26.4 min @257 ms/basin). **Full 41-yr run ≈ 40 h single-thread.** This instance = **4 cores /
> 30 GB**, so ~4 per-year workers → **~10 h here** (3.34 GB/worker × 4 ≈ 13 GB RAM, 4×2.4 GB
> disk — both fit); a larger instance scales down ~linearly with cores. Per-year download ~13 s
> (2.4 GB in-region). Peak RSS **3.34 GB/worker**. ✅ Sum-to-100 **exactly 100.000 across all 6,119 in-CONUS basins**,
> **0 unknown codes**; valid-coverage filter caught **45 AK/HI/out-of-CONUS** USGS basins.
> ✅ **C1 impervious out-of-range = 0 basins across ALL 6,164** (overall max 100.0) → V2 appears
> clean; keep the [0,100] clip defensively + assert on a couple more years during the run.
> Local extraction is **~8× faster than `/vsis3/` streaming** — confirms download/stage-then-extract.
>
> **Pilot results (22 Jul, live vs `s3://usgs-landcover` C1V2, year 2023, 36+22-basin samples):**
> ✅ `/vsis3/` requester-pays access works; keys/versions/6-products/41-yr (1985–2025) confirmed.
> ✅ CRS = `AEA WGS84` (Albers/WGS84, **no EPSG authority** — confirms M1, not 5070); reproject
>   basins to granule WKT works. ✅ nodata=250 set on both COGs. ✅ **Land cover sums to exactly
>   100.000** for US `gage_type=='USGS'` CONUS basins (18/18), **0 unknown class codes**. ✅
>   **CONUS filter = `gage_type=='USGS'`(6,164) + valid-coverage**: cleanly separates US (sum 100)
>   from Canada (sum 0, all-fill) — a lon/lat bbox wrongly kept 652 Canadian basins. ✅ Peak RSS
>   **2.44 GB** even with a 91,268 km² basin — memory is a non-issue. ⚠ **Impervious C1 bug did
>   NOT reproduce in V2-2023** (max=100.0, no >100/underflow in-sample) — keep the [0,100] clip
>   defensively (cheap; guide documents it) but it may be fixed in V2 or rare; validate across the
>   full run. ⚠ **`/vsis3/` streaming is too slow at scale** (~4.7 s/basin → ~330 h single-thread
>   for the full run) → **download-once-per-year local extraction is mandatory** (removes per-basin
>   network latency); run a local-mosaic micro-pilot to get the true per-year wall-clock.
>
> **Review fixes incorporated (22 Jul 2026):** (C1) Fractional Impervious has a *documented*
> bug — non-truncated regression predictions produce values >100 and UINT8 underflow
> (wrapped 251–255) that are **not** the 250 fill; treat valid = [0,100], mask/clip everything
> else, and validate the **raw pixel distribution** (not just the mean). (M1) CRS is **WGS84
> Albers**, not EPSG:5070 (5070 = legacy NAD83 NLCD) — reproject basins to the granule WKT
> read from `src.crs`, never hardcode. (M2) Latest coverage is **1985–2025** (C1V2), not
> 1985–2023 (that was C1V0) — 41 years. (M3) Add change-detection caveats (shrub↔grass, Great
> Lakes over-water artifacts) + smoothed deltas + optional Confidence-layer QA.

---

## 1. Goal & deliverables

**Universe** = the CONUS-US subset of the existing 7,964-watershed geometry layer
(`geometry_7964.gpkg`). NLCD Annual is **CONUS only**, so we take the `gagesii` US basins
and **drop / flag Alaska, Hawaii, PR and any basin outside the CONUS NLCD extent**;
Canadian (`wsc_eccc`) basins are out of scope. Expected ≈ 6,000 US watersheds (confirm the
non-CONUS drop count at build; GAGES-II includes AK/HI).

Two flat tables (CSV **and** parquet, mirroring the signatures / HydroATLAS / MODIS
outputs) plus a data dictionary, and a self-contained QA/QC HTML explorer.

### 1a. NLCD table — annual (`watershed_nlcd_annual_{date}.{csv,parquet}`)
One row per **(gage_id, year)**, year ∈ **1985–2025** (41 years, C1V2 — see §9 for the V0/V2
version pin). ~6,000 gages × 41 yr ≈ **~246,000 rows**, **25 columns**.

| Column group | Columns |
|---|---|
| Keys | `gage_id` (zero-padded, joins signatures), `canon_id` (zero-stripped, joins metadata), `year` |
| Land-cover class % | `nlcd_c{code}_pct` for the 16 CONUS classes (codes 11–95) — coverage-weighted % of valid (non-fill) 30 m pixels; per row sums to ~100 |
| Impervious | `impervious_mean_pct` — coverage-weighted **mean** fractional-impervious-surface % over **[0,100]-valid** pixels (basin-mean imperviousness; urbanization signal). **⚠ Mask fill 250 AND out-of-range/underflow pixels (>100, 251–255) — see §6/C1.** |
| QA | `n_nlcd_pixels` (coverage-weighted valid 30 m pixel count in basin), `valid_coverage_frac` (valid / total intersecting pixels; fill excluded from the % denominator — **net-new vs MODIS, needs a fill-inclusive denominator, see §7**), `watershed_geom_source`, `low_pixel_support`, `low_confidence` |

**Class legend (committed manifest `nlcd_legends.csv`, CONUS, 16 classes):**
11 Open Water · 12 Perennial Ice/Snow · 21 Developed Open Space · 22 Developed Low ·
23 Developed Medium · 24 Developed High · 31 Barren Land · 41 Deciduous Forest ·
42 Evergreen Forest · 43 Mixed Forest · 52 Shrub/Scrub · 71 Grassland/Herbaceous ·
81 Pasture/Hay · 82 Cultivated Crops · 90 Woody Wetlands · 95 Emergent Herbaceous Wetlands.
Fill / NoData = **250** (masked, excluded from the % denominator). Class %s are reindexed
to this committed list so missing classes become explicit `0.0` — stable schema across
years/gages (same discipline as `lulc_legends.csv`). **Unknown codes hard-fail** (a stray
AK-only code like 51 would otherwise bias the band total low).

---

## 2. Architecture & data flow

```
geometry_7964.gpkg  (existing)  ── filter to CONUS-US basins by LOCATION ──►  ~6,000 basins
        │                                                                    (EPSG:4326)
        ▼
NLCD source (§5): Annual NLCD Collection 1 (C1V2, 1985–2025)  ── boto3 / GDAL /vsis3/ ──
  s3://usgs-landcover  (us-west-2, REQUESTER-PAYS)   [Land Cover + Fractional Impervious]
        │  download the ~41 annual CONUS mosaics per layer once (COG GeoTIFF, uint8)
        ▼
Shared zonal-stats engine (§7)  [runs from us-east-2]
  reproject basins → granule WKT (WGS84 Albers, per src.crs; never resample categorical);
  windowed coverage-weighted exactextract; mask fill 250 (+ C1 out-of-range) before the % / mean
        │
        ▼
annual per-class % + basin-mean impervious   (one row per gage × year)
        │
        ▼
data_out/eo_nlcd/  CSV + parquet + data dictionary  (keyed gage_id, leading-zero-safe)
        │
        ├─►  S3  s3://climate-ai-data-science-shiny-app-data/streamflow/
        └─►  build_nlcd_explorer.py  →  self-contained Leaflet QA/QC HTML  →  S3
```

**Contrast with MODIS (what gets simpler / harder):**
- **Simpler**: COG GeoTIFF → GDAL-native, no `pyhdf`, no HDF4 driver, no Earthdata/CMR
  Bearer token. Single classification scheme (16 cols vs 102). WGS84 Albers is equal-area.
  ~41 mosaic files/layer vs ~670 MODIS tile downloads. Downloading annual CONUS mosaics
  once removes the per-tile VRT-mosaic step entirely.
- **Harder**: **30 m vs 500 m ≈ 278× more pixels per unit area** → zonal extraction is the
  compute pole (MODIS LULC was minutes–low hours; this is hours–day, pilot to measure).
  Must read **windowed** per basin, not whole-CONUS-in-memory. Requester-pays + cross-region
  egress (us-west-2 → us-east-2) — trivial in $ for ~20–40 GB but set `RequestPayer=requester`.

---

## 3. Region & credentials — runs from us-east-2

- Data lives in **`s3://usgs-landcover` (us-west-2, requester-pays)**. From our us-east-2
  instance: set `AWS_REQUEST_PAYER=requester` (and `RequestPayer='requester'` for boto3 /
  `--request-payer requester` for the CLI). No Earthdata login needed — public USGS bucket,
  standard AWS creds.
- **Process one year at a time: download that year's two CONUS mosaics → extract → delete**
  (each ~1.4 + ~0.9 GB, so peak disk ~5 GB, cumulative egress ~90–100 GB / ~$2). Keeps disk
  bounded and the run restartable. Exact keys (verified): `.../c1/v2/cu/mosaic/
  Annual_NLCD_{LndCov,FctImp}_{year}_CU_C1V2.tif`. Alternative: `/vsis3/...` windowed reads
  with `AWS_REQUEST_PAYER=requester` (no local storage) — viable since these are COGs, but
  ~6,000 basins collectively cover most of CONUS, so a whole-mosaic read per year is simpler
  and reads comparable bytes. A `public_inventory.csv` in each `cu/mosaic/` lists every object
  (+ ETag) for manifest pinning.
- Same creds discipline as the MODIS pipelines: materialize/refresh AWS creds (SageMaker
  creds rotate ~30 min); for a ProcessPool, pop static `AWS_*` env vars in workers so boto3
  resolves the auto-refreshing container role (the fix in `lai_pipeline.py`).
- **Temp staging (user, 22 Jul):** the source mosaics are copied **once** from usgs-landcover
  (requester-pays, us-west-2) into our own bucket at
  **`s3://climate-ai-data-science-shiny-app-data/streamflow/temp_lulc_conus/`** — which is
  **us-east-2, same region as compute**, so the requester-pays/cross-region egress is paid
  once and all pipeline runs/reruns read in-region (free, low-latency, either `/vsis3/` or a
  cheap in-region download). **These staged mosaics are temporary — DELETE `temp_lulc_conus/`
  once the final NLCD product is QA/QC'd** (task 7). The derived per-watershed table + explorer
  go to the permanent `streamflow/` prefix; only the bulky raw mosaics are transient.

---

## 4. Watershed geometry — reuse, filtered to CONUS

- **Reuse `geometry_7964.gpkg`** (already built, Codex-reviewed; carries `gage_id`,
  `canon_id`, `watershed_geom_source`, `low_confidence`, `latitude`, `longitude`).
- **Filter to US by `gage_type`, then to CONUS by actual NLCD coverage — NOT by a lon/lat
  bbox and NOT by `watershed_geom_source`.** Pilot finding (22 Jul): a CONUS lon/lat box
  (−125..−66.5, 24.5..49.5) wrongly keeps **652 `wsc_eccc` (Canadian) basins** whose centroids
  straddle the southern border, and keying "US" to `gagesii` would drop the ~19 US CONUS
  HydroBasins-fallback basins. **[Correction 2026-08-24: that "~19 US HB-fallback" claim is
  inconsistent with the delivered layer's own accounting (README.md §11: sources gagesii
  6,164 / wsc_eccc 1,771 / hydrobasins 29 vs gage_type USGS 6,164 / Canada 1,800 — gagesii
  == USGS exactly, and 1,771 + 29 = 1,800, so all 29 HB fallbacks are Canadian). **Settled
  2026-08-25 by the geometry rebuild** (the June layer was lost with S3; the rebuild matched
  every recorded target exactly): gagesii = 6,164 = every US gage, all 29 HB fallbacks
  Canadian. The gage_type filter remains the right choice either way.]** Correct filter: **`gage_type == 'USGS'`** (clean split —
  `gage_type` ∈ {`USGS` 6,164, `Canada` 1,800}) to get US basins including US HydroBasins-
  fallbacks, **then drop AK/HI/PR** — robustly via **valid NLCD coverage** (a basin outside the
  CONUS footprint extracts ~all fill 250 → `valid_coverage_frac ≈ 0`), backstopped by a CONUS
  lon/lat sanity check. Record the dropped set (`nlcd_out_of_conus.csv`) with reason + count.
  Canada (`gage_type=='Canada'`/`wsc_eccc`) is out of scope. (Universe is US `gage_type=='USGS'`
  ≈ 6,164 minus AK/HI; confirm the CONUS count at build.)
- Reproject basins to **the granule's native CRS read from `src.crs`** (WGS84 Albers, §5) —
  do not hardcode EPSG:5070. Our geometry is EPSG:4326/WGS84, so 4326 → WGS84-Albers is
  datum-consistent (no spurious NAD83 shift).
- Carry `watershed_geom_source` + `low_confidence` through to the output as in
  `lulc_finalize.py`; add `low_pixel_support` (n_nlcd_pixels below a threshold — at 30 m a
  basin has ~1,111 pixels/km², so **do NOT inherit the MODIS `LOWPX=5.0`** — it would never
  fire; set to e.g. a few hundred pixels ≈ <~0.3 km², confirm at pilot).

---

## 5. NLCD source — Annual NLCD Collection 1 (USGS/MRLC)

- **Product**: Annual NLCD Collection 1. Versioning uses a `C{C}V{V}` token in the filename
  and an S3 `c1/v{V}/` path: **V0** (Oct 2024) = 1985–**2023**; **V1** (Jun 2025) = 1985–2024;
  **V2** (Jun 2026, current) = 1985–**2025**. USGS's public "1.0 / 1.2" shorthand ≈ V0 / V2.
  CONUS, **30 m**, modified Anderson Level II. Six-product suite; we ingest **Land Cover** +
  **Fractional Impervious Surface** (§9). **Pin one version explicitly** (§9 → C1V2, 41 yr).
- **Format**: Cloud-Optimized GeoTIFF. Land Cover = uint8, 16 classes, fill 250. Fractional
  Impervious = uint8, valid **[0,100] %**, fill 250 — **plus documented out-of-range/underflow
  pixels (>100, wrapped 251–255)** that must be masked (§6/C1).
- **Projection**: **Albers Equal Area on the WGS84 datum** (Landsat C2 ARD CONUS Albers) —
  Albers params match 5070's (std parallels 29.5/45.5, origin −96/23) **but the datum is WGS84,
  not NAD83, so there is NO standard EPSG code** (≠ EPSG:5070, which is the legacy NLCD CRS).
  **Read the WKT per-granule from `src.crs` and reproject basins to it; never force 5070.**
  Confirm the exact WKT at pilot.
- **S3 layout** (confirm exact keys by listing at pilot): mosaics under
  `s3://usgs-landcover/annual-nlcd/c1/v{V}/cu/mosaic/` (CONUS = `cu`), ARD tiles under
  `.../cu/tile/h{xx}v{yy}/`. Example key: `.../cu/mosaic/Annual_NLCD_FctImp_1985_CU_C1V0.tif`.
  Prefer the **CONUS mosaic per year per layer**.
- **Reproducibility**: pin the exact `c1/v{V}` path + `C1V{V}` token + object ETags in an
  input manifest per run. Note two discontinuities to record: (a) re-running under a newer
  version **reprocesses historical years too** (V0→V2 can change 1985-era values); (b) the
  most-recent 1–2 years of any version come from the rule-based ICC "annual update" path, not
  the full baseline model — a methodological seam at the series tail.

---

## 6. Impervious layer semantics

`impervious_mean_pct` = coverage-weighted **mean** of the Fractional Impervious Surface
raster over the basin's valid pixels — a continuous basin-mean imperviousness, the natural
hydrologic companion to the categorical land cover (captures sub-pixel urbanization that the
16-class scheme lumps into Developed bins).

**⚠ C1 — masking (must-do):** the User Guide §5 documents that regression predictions were
*not truncated*, so the raster contains values **>100 and UINT8 underflow (wrapped 251–255)**
that are **NOT** the 250 fill. Masking only `==250` leaves a "255%" pixel in the mean and
silently biases low-impervious rural basins high — and a mean can still land in [0,100] while
being wrong, so the mean-range gate won't catch it. **Rule: valid = raw ∈ [0,100]; mask/clip
everything else ({250} ∪ >100 ∪ 251–255) BEFORE the coverage-weighted mean.** Validate the
**raw per-granule pixel-value distribution** at pilot (assert no valid pixels outside
[0,100]∪{250}). Same discipline the MODIS LAI pipeline used (mask raw ≥101 before scaling).
(Optional future add: `impervious_p90` for spatial spread — not in v1.)

---

## 7. Shared zonal-stats engine (30 m specifics)

Adapt the MODIS engine; key differences for 30 m:
1. **Reproject basins → the granule WKT read per-granule from `src.crs`** (WGS84 Albers, §5;
   NOT EPSG:5070). Never resample the categorical land-cover raster.
2. **Windowed reads per basin** — point `exactextract` at the local mosaic (or `/vsis3/`) as
   a **lazy dataset/file, never a pre-loaded array**, so it reads only each basin's window and
   the ~15 GB-uncompressed CONUS mosaic is never materialized. Memory is driven by the basin's
   **bounding-box window, not its area** — a sprawling basin's bbox can be several× its area
   (~0.5–1 GB), so confirm peak RSS on the largest/most-sprawling basin at pilot.
3. **Land cover**: `exact_extract(lc, basins, ["unique","frac","count"])` → per-class
   coverage-weighted fraction ×100 (+ `n_nlcd_pixels` from count). Reindex to the 16-class
   manifest; **unknown code → hard-fail** (catches AK-only 51/72–74 leaking through the CONUS
   filter). Mask fill 250.
4. **Impervious**: `exact_extract(imperv, basins, ["mean"])` coverage-weighted, over
   **[0,100]-clipped** pixels (§6/C1) → `impervious_mean_pct`.
5. **`valid_coverage_frac` is net-new** (the shipped MODIS code emits only `n_modis_pixels`,
   never a valid fraction). It needs a **fill-inclusive denominator**: either a second
   exactextract pass with nodata unset, or `polygon_area / 900 m²`. Specify the mechanism in
   code (or drop the column) — don't assume it copies over.
6. **Confirm the COG advertises `nodata = 250`.** exactextract's `unique/frac/count/mean` only
   exclude fill if the dataset tags nodata; if the FIS COG omits it (plausible given the
   underflow bug), set it explicitly (as the MODIS code did via `srcNodata`/`VRTNodata`), else
   fill pollutes the denominator/mean. Verify at pilot.
7. **Checkpointing**: per-year parquet (`nlcd_{year}.parquet`), skip-if-exists, so the long
   run resumes after a crash — same pattern as `lulc_pipeline.py`. A year is written only if
   both layers extracted cleanly for all basins.
8. **No multi-tile mosaic math needed** if using the CONUS mosaic (single raster). If ARD
   tiles are used instead, per-tile coverage-weighted **sums** combine exactly for both mean
   and class-fraction (no quantiles here), so tiles are fine too; the single mosaic is simpler.

---

## 8. Compute, cost & runtime (rough; confirm via pilot)

| | NLCD (Annual C1V2) |
|---|---|
| Files | ~41 yr × 2 layers = **82 CONUS mosaic COGs** (verified in-bucket 22 Jul) |
| Volume / egress | **measured in-bucket**: LndCov ~1.4 GB/yr, FctImp ~0.8–1.0 GB/yr → **~90–100 GB** cumulative for both layers (not the earlier 20–40 GB guess; COGs are ~10× compressed but the CONUS grid is ~14.5 B pixels). Requester-pays us-west-2→us-east-2 egress ≈ **~$2**. **Bound peak disk with a per-year download→extract→delete loop (~5 GB resident), not an 82-file bulk download.** ⚠ The Confidence layer is **~7.7 GB/yr** (~315 GB for 41 yr) — prohibitive; use only sampled if at all (M3). |
| Compute | **the pole** — 30 m exactextract over ~6,000 basins × 41 yr × 2 layers; **hours–day**, heavily parallel per-year; **pilot 2–3 years to price** before committing the full run |
| Rows out | ~6,000 gages × 41 yr ≈ **~246,000 rows**, 25 cols |

**Pilot checklist (deliverables):** wall-clock per year; **peak memory on the largest/most-
sprawling basin** (bbox-window, not area); **confirm the granule WKT** (WGS84 Albers, not
5070); **confirm the COG nodata tag = 250**; exact `c1/v{V}` S3 keys; **raw impervious pixel
distribution** (no valid pixels outside [0,100]∪{250} → validates the C1 mask); class %s sum
to 100 over valid pixels; impervious mean ∈ [0,100] after clipping.

---

## 9. Decisions

**Locked (user, 22 Jul 2026):**
1. **Layers** = **Land Cover + Fractional Impervious Surface** (not the full six-product
   suite; Change/DOY use fill 9999 and change is derivable from the annual series).
2. **Period** = **full record to current = 1985–2025 (41 yr), pinned to C1V2.** User asked for
   "1985 to current"; current is now 2025 (V2, Jun 2026), not 2023. ⚠ **CONFIRM the version
   pin:** C1V2 gives the best-available estimates + 2024–25 but reprocesses history and moves
   with each release; **C1V0 (1985–2023)** is the alternative if a frozen, maximally
   reproducible baseline is preferred. Legend + pilot are identical either way.
3. **Deliverables** = **full parity** — pipeline + table (parquet/CSV + dictionary on S3) +
   self-contained Leaflet QA/QC explorer.
4. **Language** = Python (mirrors the MODIS EO stack).
5. **Universe** = CONUS-US basins by **location** (not geom source); drop/flag AK/HI/non-CONUS;
   Canada out of scope.
6. **Access** = stage source mosaics once into `streamflow/temp_lulc_conus/` (our us-east-2
   bucket, in-region with compute), extract from there; **delete the temp staging after final
   QA/QC** (task 7). Avoids repeated requester-pays/cross-region pulls on reruns.

**Locked build requirements from the design review (22 Jul 2026):**
7. **Impervious mask** = clip to [0,100], mask {250}∪>100∪251–255; validate raw distribution (C1).
8. **CRS** = reproject basins to the granule WKT (WGS84 Albers), never hardcode EPSG:5070 (M1).
9. **Change caveats** = add shrub↔grass / Great-Lakes-over-water caveats to dict + explorer;
   report smoothed (mean-first-3 vs mean-last-3) deltas; optionally ingest Confidence as QA (M3).

**Resolve at pilot:**
10. Confirm the granule WKT, exact `c1/v{V}` S3 keys, and the COG nodata tag.
11. `low_pixel_support` threshold at 30 m (NOT 5); non-CONUS drop count; `valid_coverage_frac`
    denominator mechanism.
12. Whether to add `impervious_p90` (spatial spread) — defer unless wanted.

---

## 10. Module layout (direct analogs of the LULC modules)

```
EO_data_processing/
├── README_NLCD.md              # this plan
├── eo_processing/
│   ├── nlcd_legends.csv         # 16 CONUS classes: code, class_name, output_column, hex_color
│   ├── nlcd_pipeline.py         # download CONUS mosaics + windowed exactextract → annual per-class % + impervious mean
│   └── nlcd_finalize.py         # + geom source / low_confidence / low_pixel_support, dictionary, date-stamped parquet+csv, S3 upload
└── viz/
    └── build_nlcd_explorer.py   # self-contained Leaflet HTML QA/QC explorer (NLCD palette)
```
Outputs → `data_out/eo_nlcd/` (gitignored); products → `s3://climate-ai-data-science-shiny-app-data/streamflow/`.

---

## 11. Task sequence (milestones)

1. **`nlcd_legends.csv`** — commit the 16 CONUS classes (code, name, `nlcd_c{code}_pct`,
   official hex color for the explorer). Version-independent → build now.
2. **Access + pilot** — list `s3://usgs-landcover` (requester-pays) to confirm `c1/v{V}` keys
   + pin the version; download 2–3 land-cover + impervious CONUS mosaics; **confirm the granule
   WKT (WGS84 Albers) + COG nodata tag + raw impervious pixel distribution**; run extraction on
   a stratified basin sample; price wall-clock/**peak memory on the largest/most-sprawling
   basin**; validate sum-to-100 + impervious mean ∈ [0,100] after the C1 clip.
3. **`nlcd_pipeline.py`** — full 1985–2025 (C1V2) run, per-year checkpoints, both layers,
   **C1 impervious clip**, **basins reprojected to granule WKT**, unknown-code hard-fail,
   completeness guard (write a year only if both layers extracted for all basins).
4. **✅ DONE (local, 24 Jul) — `nlcd_finalize.py`** — drops 45 Alaska out-of-footprint gages to
   `nlcd_out_of_footprint` (never published as real zeros) + all-gage QA companion; separate flags
   (geom_low_confidence / low_pixel_support<100px / partial_coverage<0.99 / low_confidence union);
   `nlcd_collection='C1V2'` + `valid_area_km2` + metadata JSON; dictionary with M3 + non-additive
   caveats; date-stamped parquet/CSV; join `validate="many_to_one"` + canon_id equality assert.
   **S3 upload opt-in (`--upload`), HELD for QA.** Output: 250,879 rows = 6,119 CONUS gages × 41 yr.
5. **Validation gates** — class %s sum ≈ 100 over valid pixels; impervious mean ∈ [0,100] AND
   no raw pixel outside [0,100]∪{250}; expected (gage × year) completeness; 100% leading-zero-
   safe `gage_id` join to signatures; non-CONUS drops reported; **Great-Lakes/coastal basins
   QA-checked** for over-water artifacts. Emit a QA report.
6. **✅ DONE (local, 24 Jul) — `build_nlcd_explorer.py`** → `nlcd_explorer.html` (12.6 MB): 6,119
   points, 32 map vars (dominant class/group, normalized Shannon, smoothed endpoint-diff deltas
   [labeled, seam-disclosed], impervious mean, annual_volatility / shrub_grass_swap artifact views,
   per-class %, QA flags); click a gage → stacked-area composition 1985–2025 + dashed impervious
   line, 2024–25 update-seam shaded, per-year hover; partial-coverage badges + "only fully-covered"
   filter; glossary with M3 + non-additive caveats + excluded-45 disclosure. Structural + JS-syntax
   validated (browser render = the human-QA step; no Chrome on this box).
7. **Cleanup + CHANGELOG + memory** — **after final QA/QC of the NLCD product, DELETE
   `s3://climate-ai-data-science-shiny-app-data/streamflow/temp_lulc_conus/`** (the staged raw
   mosaics, ~90–100 GB — transient by design). Document the new product in CHANGELOG; note it
   complements (does not replace) the CONUS+Canada MODIS LULC.

Environment (us-east-2): reuse the `geo` conda env (rasterio/geopandas/exactextract/pyarrow/
boto3). **No `pyhdf`, no earthaccess needed** — GDAL reads NLCD COGs natively.

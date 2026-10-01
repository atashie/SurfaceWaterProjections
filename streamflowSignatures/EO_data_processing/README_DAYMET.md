# Daymet basin-mean daily climate (HISSS climate input, reprocessed 2026)

Daily basin means of six Daymet V4 R1 variables for **8,017 watersheds**, calendar
**1980–2025**, computed in this repo from the gridded Daymet mosaics. It is meant to replace
the earlier co-author product (gdptools, 6,087 basins, 1980–2023) as the climate input of
the signature pipeline. **As of 2026-10-01 neither delivered product uses it.** Whether to
rerun them on it is an open user decision (`docs/STATUS.md`). The user judged the raw-data
agreement below sufficient and dropped a signature-level replay (2026-10-01). Tools: `daymet/` (see
`daymet/README.md`); polygons: `geometry/`.

## Source
- **Files.** Daymet V4 R1 daily, North America, 1 km (ORNL DAAC, doi:10.3334/ORNLDAAC/2129;
  NASA CMR collection C2532426483-ORNL_CLOUD), one NetCDF-4 file per variable and year. The
  2026-09-29 run took 269 files from ORNL (Earthdata Login). The other 7 came from the NCAR
  GDEX mirror (dataset d682806): prcp 1980 and the six 2023 files, downloaded for the gate.
- **Checksums.** Every file's SHA-256 matched NASA CMR's. A CMR re-query on 2026-10-01 lists
  the same size and SHA-256 for all 276. The mirror's 270 files for 1980–2024 have exactly
  CMR's sizes (2025 is ORNL-only), but only the 7 used were checksummed.
- **Calendar.** 365 days per year: Daymet drops Dec 31 in leap years, and the output keeps
  that calendar.
- **Temperature floor.** Daymet's tmin/tmax cell values bottom out at −60 °C in most years.
  No basin mean reaches that (lowest −58.4 °C).
- **Canadian station change.** From the 2024 annual update, Canadian inputs come from
  ECCC/CCCS stations, bias-corrected, instead of the GHCNd `CA0` stations (ORNL
  documentation; not confirmed for 2025). That is a station-network seam to keep in mind for
  Canadian trend endpoints. The processing adds none:
  - ORNL's own cell values reproduce the file in 2024–2025 (below);
  - the Canada/US prcp ratio drifts from about 2015 with no step at 2024.

## Polygons
- **The basins.** The 7,964 basins of the delivered boundary layer (HydroShare Resource 4)
  at **full resolution** (Resource 4 is their 200 m simplification), plus the **53 Canadian
  basins > 100,000 km²** that Resource 4 leaves out. 05KH009 stays excluded: its HydroBASINS
  fallback polygon is wrong.
- **Build.** `geometry/rebuild_watershed_polygons.py --no-simplify --include-large` (its
  defaults reproduce Resource 4). On 2026-10-01 the layer rebuilt byte-identically from the
  recorded inputs; their md5s are in the run's
  `polygons/*.provenance_rebuild_2026-10-01.json`.
- **Validity.** All 8,017 polygons are valid and non-empty. Their ids match the streamflow
  parquet exactly.

## Method
- **Grid.** Polygons are reprojected to the Daymet Lambert Conformal Conic grid; the raster
  is never resampled.
  - The CRS comes from the CF attributes of the file the weights were built on (prcp 2023).
  - Every file's CRS parameters, grid corners, spacing and calendar are asserted equal to
    the Daymet NA grid. That is enforced from 2026-10-01. In the 2026-09-29 run the 276
    probe records carry the identical CRS.
- **Weights.** exactextract gives each polygon's coverage fraction of every 1 km cell. Each
  weight is coverage × the cell's true area: the grid is conformal, so a projected 1 km²
  cell covers 0.83–1.10 km² on the ellipsoid.
- **Daily value.** Σ w·v / Σ w over the cells whose value is not −9999 that day.
- **Reproducibility.** Every HDF5 chunk is decompressed once. With a fixed weights file the
  result is bit-reproducible across reruns and worker counts. A different weights file (more
  basins) changes the summation order and so the last bits (≤ 2e-13 relative).

| Column | Unit |
|---|---|
| `site_id` | gage id as spelled in the streamflow parquet (the join key) |
| `Date` | day |
| `prcp` | mm/day |
| `tmin`, `tmax` | °C |
| `swe` | kg/m² (= mm) |
| `vp` | Pa |
| `srad` | W/m², daylight average |

Rows are ordered (year, site_id, Date), so a site is not contiguous across years.

## Validation (single-year gate, calendar 2023 + prcp 1980)
- **exactextract's own weighted mean,** recomputed on five days per variable for all 7,964
  basins of the gate layer, agrees with the chunked aggregation to ≤ 5e-11 absolute
  (≤ 3e-13 relative).
  - The check shares the polygons, the CRS and the cell-area function with the code under
    test, so it verifies the chunked accumulation.
  - The 53 basins added later were not cross-checked.
- **ORNL's Single Pixel API at six points:** projected coordinates identical to 0.00 m, and
  pixel series identical to float32 rounding.
- **Against the co-author product,** on the 5,965 basins it has data for (5,969 shared, 4 of
  them all-NaN there):
  - per-basin daily R² ≥ 0.99999995 for prcp, tmin, tmax, vp and srad;
  - swe ≥ 0.99979 (a basin with 4e-8 mm mean SWE);
  - prcp annual totals within ±0.0011 % (1st–99th percentile);
  - no date shift.
  prcp 1980 (the 1980–2019 file layout, a leap year) agrees equally.
- **How the co-author product was made, as far as its numbers show** (inferred from the
  match, not confirmed by the co-authors):
  - It behaves as if computed on the unsimplified polygons with true-area weights.
  - It behaves as NaN for any basin touching a fill cell. 02234500, 02236000, 02236125 and
    02244040 are NaN on every day. So is 01372058, which has no polygon; cause unknown.
    This product averages the valid cells instead.
  - The alternatives do not reproduce it. Coverage-only weights would shift large northern
    basins by up to 0.3 %. The simplified Resource 4 polygons shift prcp totals by about
    0.1 % (p01–p99) and swe by about 2.5 %.

## The file and its full-run validation (run 2026-09-29/30; reviewed 2026-10-01)
**The file.** `daymet_1980_2025_29sep2026.parquet` in `~/HISSS_data/daymet-processed-29sep2026/`.
- 134,605,430 rows = 8,017 sites × 46 years × 365 days, no NaN in any variable.
- 5,134,193,404 bytes, md5 `c059076088cd8abd963c64094468e90f`.
- An md5-verified copy is in `/Volumes/Untitled/daymet-processed-29sep2026/`; the
  `RUN_NOTES.md` there has the full tables.
- Its `.provenance.json` records: the inputs and code commits; the source granules' CMR
  checksums; the weights and polygon md5s; the verification below.

**What the agreement means.** Both products read the same Daymet V4 R1 cells over the same
official GAGES-II and ECCC outlines. Agreement therefore shows that this toolchain
reproduces the co-authors' aggregation; it is not an independent test of Daymet's accuracy.
The weighting and the polygon version were chosen in the 2023 gate because they matched, so
1980–2022 are out of sample for that choice.

**Against the co-author product, 1980–2023,** on the 5,965 basins it has data for (262,460
basin-years per variable), daily values:
- **Largest daily differences:**
  - 0.015 mm (prcp), 0.003 °C (tmin, tmax), 0.26 Pa (vp), 0.09 W/m² (srad), 0.26 mm (swe).
  - The tmin, tmax, vp and srad maxima all come from one 2.9 km², 8-cell basin (08GA061).
    Without it: 0.0016 °C, 0.13 Pa, 0.08 W/m².
- **prcp annual totals** within ±0.0011 % (p01–p99), ±0.012 % at worst.
- **Identity R² per basin-year** is ≥ 0.9999999 for prcp, tmin, tmax, vp and srad. The
  seasonal cycle dominates it, so the absolute differences above are the informative
  numbers.
- **swe:** R² ≥ 0.999 in all but 2 basin-years, both trace snow (mean < 1e-6 mm). 25,105 swe
  basin-years are all-zero in both products, identically.
- **No date offset.** Every basin-year joins on all 365 dates with this agreement. A ±1-day
  lag test never preferred a shift, but on its own it has little power.

`daymet_validate.py` compared the per-variable-year files. The assembled parquet was then
verified equal to them in all 807,632,580 values (`daymet_assemble.py --provenance-only`,
2026-10-01).

**Not compared with any earlier product:**
- 2,048 of the 8,017 basins, among them all 53 basins > 100,000 km², all 30
  HydroBASINS-fallback polygons and 47 of the 57 low-confidence polygons;
- calendar 2024–2025, for every basin;
- the 7 basins that touch Daymet fill cells (≤ 0.14 % of their weight).

For these, the evidence is ORNL's Single Pixel API. The basin means of 8 small basins (4
USGS, 4 Canadian; 3–6 cells), recomputed from ORNL's own cell series, match the file in
1980, 1996, 2012, 2019, 2020, 2023, 2024 and 2025:
- to ≤ 1.1e-4 for vp;
- to ≤ 3e-5 for the other variables.
The check is `daymet/daymet_outputcheck.py`; its result is in the run's
`validation/outputcheck_single_pixel_2026-10-01.json`. The 2026-10-01 reviews' own API
checks of further small basins agree. The large basins rest on the weight checks: Σ coverage
× true area equals the polygon area to 1e-7 for all 53.

**Explorer.** `validation/daymet_record_explorer_2026-10-01.html` in the run folder is
self-contained (10.9 MB). It shows these metrics as distributions and maps every basin. For
30 basins it plots the original and new daily series side by side: the 10 least matching,
10 random and 10 without an original series. `viz/build_daymet_record_explorer.py` builds
it.

**Coverage.**
- Every gage of both delivered products is in this file (#1 6,678, #2 6,250).
- The stale file has 5,517 / 5,638 of them; 4 / 3 of those are all-NaN (St. Johns River, FL).
- A rerun on this file would give a usable climate series to 1,165 (#1) / 615 (#2) gages:
  the 1,161 / 612 absent from the stale file plus the 4 / 3 all-NaN ones.
- The stale file's other 118 sites belong to neither product.

## Caveats
- **Low-confidence polygons, kept and flagged** (user decision 2026-10-01). 57 polygons are
  `low_confidence`: 30 HydroBASINS fallbacks and 28 whose area is > 50 % off the metadata area
  (one is both). 12 of them are in product #1 and 13 in #2.
  - Example: 06FD001, 28,997 km² against a reported 289,000.
  - 27 of the HydroBASINS polygons have no metadata area to check against.
  - The flags are in `daymet_basin_flags.csv`, next to the climate file. It has one row per
    flagged basin (101, the 45 small basins included): the reason, both areas, and product
    membership. `daymet/daymet_basin_flags.py` builds it.
  - At a rerun the table goes with the products as a companion file; the signature CSV's
    column contract is unchanged.
- **Small basins.** 45 basins cover fewer than 4 grid cells (`low_pixel_support`; 35 in
  #1, 21 in #2).
- **No streamflow.** 4 basins have no streamflow record (01591000, 01591400, 01591610,
  01591700); they keep their boundary id.
- **Memory.** The file is 37 % larger than the stale input. The Julia runner reads all eight
  columns before keeping four, so expect about 8.6 GB resident (extrapolated from one row
  group). On the 16 GB laptop, run it in gage batches with
  `docs/benchmarks/run_batched_julia.py`, which writes 4-column batch inputs. In the
  2026-10-01 control run (stale input), each of 6 batches (~1,336 gages) peaked at
  4.8–5.5 GB RSS.

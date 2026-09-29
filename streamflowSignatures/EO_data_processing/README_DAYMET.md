# Daymet basin-mean daily climate (HISSS climate input, reprocessed 2026)

Daily basin means of six Daymet V4 R1 variables for **8,017 watersheds**, calendar
**1980–2025**, computed in this repo from the gridded Daymet mosaics. It replaces the earlier
co-author product (gdptools, 6,087 basins, 1980–2023) as the climate input of the signature
pipeline. Tools: `daymet/` (see `daymet/README.md`); polygons: `geometry/`.

## Source
- Daymet V4 R1 daily, North America, 1 km (ORNL DAAC, doi:10.3334/ORNLDAAC/2129): one
  NetCDF-4 file per variable and year. 1980–2024 from the NCAR GDEX mirror (dataset
  d682806), 2025 from ORNL (Earthdata Login). Every file's SHA-256 matches NASA CMR's.
- 365 days per year: Daymet drops Dec 31 in leap years; the output keeps that calendar.
- Canadian inputs change in 2024–2025 (ECCC/CCCS stations replaced the GHCNd `CA0` stations,
  bias-corrected) — a station-network seam to keep in mind for Canadian trend endpoints.

## Polygons
The 7,964 basins of the delivered boundary layer (HydroShare Resource 4) at **full
resolution** (Resource 4 is their 200 m simplification), plus the **53 Canadian basins
> 100,000 km²** that Resource 4 leaves out. 05KH009 stays excluded: its HydroBASINS
fallback polygon is wrong. Build: `geometry/rebuild_watershed_polygons.py --no-simplify
--include-large` (its defaults reproduce Resource 4).

## Method
Polygons are reprojected to the Daymet Lambert Conformal Conic grid (read from each file's
CF attributes); the raster is never resampled. exactextract gives each polygon's coverage
fraction of every 1 km cell; each weight is coverage × the cell's true area (the grid is
conformal, so a projected 1 km² cell covers 0.83–1.10 km² on the ellipsoid). The daily
basin value is Σ w·v / Σ w over the cells whose value is not −9999 that day. Every HDF5
chunk is decompressed once; the result is bit-reproducible across reruns and worker counts.

| Column | Unit |
|---|---|
| `site_id` | gage id as spelled in the streamflow parquet (the join key) |
| `Date` | day |
| `prcp` | mm/day |
| `tmin`, `tmax` | °C |
| `swe` | kg/m² (= mm) |
| `vp` | Pa |
| `srad` | W/m², daylight average |

## Validation (single-year gate, calendar 2023 + prcp 1980)
- exactextract's own weighted mean, recomputed independently on five days per variable for
  every basin, agrees to ≤ 1e-13 relative.
- ORNL's Single Pixel API at six points: projected coordinates identical to 0.00 m, pixel
  series identical to float32 rounding.
- Against the co-author product on the 5,969 shared basins: per-basin daily R² ≥ 0.99999996
  for prcp, tmin, tmax, vp and srad; swe ≥ 0.99979 (a basin with 4e-8 mm mean SWE); prcp
  annual totals within ±0.0011 % (1st–99th percentile); no date shift. prcp 1980 (the
  1980–2019 file layout, a leap year) agrees equally.
- Differences from the co-author product: it returned NaN for any basin touching a single
  fill cell, so 02234500, 02236000, 02236125, 02244040 and 01372058 were NaN on every day;
  this product averages the valid cells. Its basins were the unsimplified polygons, and
  its weights true areas: coverage-only weights would shift large northern basins by up to
  0.3 %, the simplified Resource 4 polygons prcp totals by ~0.1 % and swe by ~2.5 %.

Full-run results (1980–2023 comparison, timings): added when the run completes.

## Caveats
- 45 basins cover fewer than 4 grid cells (`low_pixel_support` in the basin table).
- 4 basins have no streamflow record (01591000, 01591400, 01591610, 01591700); they keep
  their boundary id.
- The 53 large basins have no earlier product to compare against and no polygon in
  Resource 4.

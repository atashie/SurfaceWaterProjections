# Daymet basin aggregation tools

Daily basin means of Daymet V4 R1 (prcp, tmin, tmax, swe, vp, srad) for the HISSS
watershed polygons, computed from the annual North America mosaics. Product description,
method and validation: `../README_DAYMET.md`. Run artifacts go in the run folder, never
here.

## Environment
```
uv venv --python 3.12 pyenv
uv pip install --python pyenv/bin/python -r requirements.txt
```
Files come from ORNL (Earthdata Login) when auth is available — about three times faster
than the anonymous NCAR GDEX mirror from the M5 laptop — else from the mirror, which holds
1980–2024 byte-identically; every file is checked against NASA CMR's SHA-256. Calendar 2025
exists only at ORNL. Auth: a bearer token in `~/.config/earthdata/edl_token` (mode 600;
never in the repo) or `~/.netrc` (`machine urs.earthdata.nasa.gov login … password …`).

## Tools, in run order
| Tool | Does |
|---|---|
| `../geometry/rebuild_watershed_polygons.py` | the polygon layer; `--no-simplify --include-large` gives the 8,017-basin full-resolution layer the climate uses (defaults reproduce Resource 4) |
| `daymet_weights.py` | once per polygon layer: weights (basin, cell, coverage, area_km2) + per-basin table with the streamflow `site_id` spelling |
| `daymet_stream.py` | per variable-year: fetch (resumable; ORNL-only files first) → checksum vs CMR → probe → aggregate → delete the raw file; `--dry-run [--check-mirror]` first |
| `daymet_probe.py` | layout facts; fails if an unstored chunk lies under a basin |
| `daymet_aggregate.py` | one file → `<var>_<year>.parquet` + QA + timing; area weights; bit-reproducible |
| `daymet_assemble.py` | per-variable-years → the pipeline's climate parquet (`site_id, Date, prcp, tmin, tmax, swe, vp, srad`) |
| `daymet_validate.py` | per-basin comparison with the stale co-author product (1980–2023) |
| `daymet_crosscheck.py` | independent check: exactextract's own weighted mean on a few days |
| `daymet_pixelcheck.py` | georeferencing check against ORNL's Single Pixel API |
| `copy_verify.py` | copy run outputs to the exFAT drive and re-read every file from the device (md5 + size) |
| `../viz/build_daymet_comparison_dashboard.py` | HTML dashboard: fresh vs original values and daily series |

## Resources (16 GB M5 laptop, measured 2026-09-29)
A full 1980–2025 run downloads 3.455 TB (~50 h at the mirror's ~19 MB/s). The 2026-09-29
run took 25.9 h from ORNL: median 50.5 MB/s per file, but 10 % of the files crawled at
3–30 MB/s on one slow connection (6.8 h lost) — curl's stall floor (200 KB/s for 60 s)
does not catch that; a higher floor that drops and resumes the connection would. Then
assemble 1 min, validate 7 min, drive copy 27 min (13 GB). Aggregation takes 15–56 s per
variable-year and hides behind the downloads. With the 8,017-basin weights (45.8 M
entries) use `--workers 6`: peak 5.0 GB parent + 1.5 GB per worker. The one-time weights
build peaks at ≈ 8 GB. Start long runs under `caffeinate -i`.

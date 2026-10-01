# Daymet basin aggregation tools

Daily basin means of Daymet V4 R1 (prcp, tmin, tmax, swe, vp, srad) for the HISSS
watershed polygons, computed from the annual North America mosaics. Product description,
method and validation: `../README_DAYMET.md`. Run artifacts go in the run folder, never
here.

## Environment
```
uv venv --python 3.12 pyenv
uv pip install --python pyenv/bin/python -r requirements.txt
pyenv/bin/python selftest_daymet_tools.py      # 23 checks, ~15 s, no network
```

**Download sources.** Files come from ORNL (Earthdata Login) when auth is available,
about three times faster than the anonymous NCAR GDEX mirror from the M5 laptop.
Otherwise they come from the mirror, whose 1980–2024 files have exactly CMR's sizes.
Calendar 2025 exists only at ORNL. Every file is checked against NASA CMR's SHA-256.

**Auth.** One of:
- a bearer token in `~/.config/earthdata/edl_token` (mode 600). `daymet_stream.py` hands it
  to curl on stdin, never on a command line or in a file. Write the file from your own
  terminal: a token pasted into a Claude session stays in its transcripts. Revoke it when
  the run is done.
- `~/.netrc` (`machine urs.earthdata.nasa.gov login … password …`).

## Tools, in run order
| Tool | Does |
|---|---|
| `../geometry/rebuild_watershed_polygons.py` | the polygon layer; `--no-simplify --include-large` gives the 8,017-basin full-resolution layer the climate uses (defaults reproduce Resource 4); writes an input-hash `.provenance.json` |
| `daymet_weights.py` | once per polygon layer: weights (basin, cell, coverage, area_km2) + per-basin table with the streamflow `site_id` spelling |
| `daymet_stream.py` | per variable-year: fetch (resumable; ORNL-only files first) → checksum vs CMR → probe → aggregate → delete the raw file; `--dry-run [--check-mirror]` first. Re-queries CMR on every start and stops on a changed record; refuses `.done` files built on other weights; drops and resumes a connection below `--min-speed-mbps` (20) at once; HTTP 401/403/404 end a source; terminates curl on any exit |
| `daymet_probe.py` | layout facts; fails if an unstored chunk lies under a basin, or (through `daymet_common.read_grid`) if the grid, CRS or calendar differs |
| `daymet_aggregate.py` | one file → `<var>_<year>.parquet` + QA + timing (with the verified source SHA-256, weights md5, commit); area weights; bit-reproducible for a fixed weights file; an out-of-range cell value is fatal |
| `daymet_assemble.py` | per-variable-years → the pipeline's climate parquet (`site_id, Date, prcp, tmin, tmax, swe, vp, srad`, ordered year, site_id, Date); re-reads it and verifies every value against the inputs; `--provenance-only` upgrades an existing run's sidecar |
| `daymet_validate.py` | per-basin comparison with the stale co-author product (1980–2023) |
| `daymet_outputcheck.py` | the assembled file against ORNL's Single Pixel API for a few small basins (any year; no raw files) |
| `daymet_basin_flags.py` | `daymet_basin_flags.csv`, the companion table of questionable basins (HydroBASINS fallback, area mismatch, no area reference, < 4 cells) with product membership |
| `../viz/build_daymet_record_explorer.py` | self-contained HTML explorer of the whole record: the agreement metrics as distributions, a basin map, and the original and new daily series overlaid (with an optional anomaly axis, new − original, and zoom on dates and values) for 3 × N embedded basins (least matching, random, no original). It checks the embedded series against the validation table |
| `daymet_crosscheck.py` | independent check: exactextract's own weighted mean on a few days (needs a raw file) |
| `daymet_pixelcheck.py` | georeferencing check of a raw file against ORNL's Single Pixel API |
| `copy_verify.py` | copy run outputs to the exFAT drive and re-read every file from the device (md5 + size); re-copying replaces the manifest line |
| `../viz/build_daymet_comparison_dashboard.py` | HTML dashboard: fresh vs original values and daily series |
| `selftest_daymet_tools.py` | synthetic-file and local-server tests of the above (grid asserts, aggregation, assembly verification, manifest, download path) |

After a run: `daymet_assemble.py` (it verifies), `daymet_validate.py` for the years the stale
product covers, then `daymet_outputcheck.py --years <a spread incl. the newest>`,
`../viz/build_daymet_record_explorer.py` and `copy_verify.py`.

## Resources (16 GB M5 laptop, measured 2026-09-29)
**Downloads.** A full 1980–2025 run downloads 3.455 TB (~50 h at the mirror's ~19 MB/s).
The 2026-09-29 run took 25.9 h from ORNL, median 50.5 MB/s per file. But 26 of the 269
files crawled at 3–30 MB/s, mostly one at a time between full-speed files, and lost 6.8 h.
The old 200 KB/s stall floor never caught that. The default floor is now 20 MB/s for
30 s, with an immediate resume, which would have recovered most of the loss.

**Processing.**
- Afterwards: assemble 1 min (plus ~40 s for its verification), validate 7 min, drive copy
  27 min (13 GB).
- Aggregation takes 15–56 s per variable-year and hides behind the downloads.
- With the 8,017-basin weights (45.8 M entries) use `--workers 6`: peak 5.0 GB parent +
  1.5 GB per worker. The one-time weights build peaks at ≈ 8 GB.
- Start long runs under `caffeinate -i`.

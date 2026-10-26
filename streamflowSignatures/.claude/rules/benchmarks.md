---
paths:
  - "docs/benchmarks/**"
  - "julia/src/config.jl"
  - "julia/src/io.jl"
---

# Benchmarks, experiments and standard products

- `docs/benchmarks/` holds TOOLS and the long-lived cross-language reference CSVs only.
  Run artifacts (signatures CSV, annual parquet, timing JSON, log, explorer, dashboards,
  comparison CSVs) go in the run's own folder, one folder per experiment (user
  convention, 2026-07-22).
- Two delivered standard products, both 1,653 columns at 60 % qualifying fraction:
  WY 1993–2025 (`processedOuts_drought_28jul2026`, 6,678 gages) and WY 1980–2025
  (`processedOuts_1980_2025_11aug2026`, 6,250). Neither is a subset of the other, and
  record-dependent signatures (drought thresholds, `*_all` pulses, elasticity,
  parameterized BFI) are never compared across them nor re-aggregated onto another
  window. Never rewrite a delivered product without an explicit user decision
  (`docs/STATUS.md` records the `flagged_for_high_na` decision).
- `STREAMFLOW_CONFIG` is read at PRECOMPILE time; the window, fraction, path and output
  overrides are read at runtime inside `main()`. A config-variant result obtained
  without a cache purge, a probe, or an observed expected delta is untrustworthy.
- Climate input since 2026-10-12 is the reprocessed `daymet_1980_2025_29sep2026.parquet`
  (8,017 basins, six variables; the runner keeps prcp, swe, tmin, tmax); the 2026
  products used `daymet_1980_2023_rebuilt_10aug2026.parquet`; the canonical-named
  `daymet_1980_2023.parquet` is truncated. Verify parquet byte sizes and `PAR1` footers
  against the last timing JSON's `provenance` block before any long run; prefer the
  verified local copies under `~/HISSS_data/` to the exFAT drive.
- `run_batched_julia.py merge` sorts the annual parquet in DuckDB under `--memory-limit`
  (default 1GB): the 46-year product's 34.7 M rows need `--memory-limit 3GB` (2026-10-13).
  Never trust a chain's echoed exit codes — `echo "[$(date)] exit=$?"` reports the date
  substitution's status; read `merge.log` and the `validation_*.txt` files.
- Gates vs diagnostics: `check_schema_equality.py`, `check_annual_parquet_equality.py`,
  `check_summary_equality.py` (summary CSV cell for cell, gate state always; tie-bearing
  rank statistics and propagated annual residuals reported by name, 2026-10-26),
  `check_signature_failures.py` and `check_additivity.jl` are GATES (non-zero exit,
  waivers named on the command line). The `compare_*` scripts intersect columns and
  exit 0 on a missing family — diagnostics only. `validate_annual_values.py` is
  self-referential (necessary, not sufficient).
- Qualifier rerun tool chain (2026-10-26; plan `docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md`):
  `splice_nwis_qualifiers.py` (list → fetch → splice → control → audit), `qualifier_year_oracle.py`
  (replay / oracle), `climate_frame_inventory.py` and `preprocess_state_ledger.py` (B′ ledgers and
  state tables), `assess_rerun_changes.py` (the change census, duckdb venv), `build_stage_cohorts.py`
  (S1 / S1′ manifests), `docs/qualifier_rerun_registry.csv` (the dependency registry; docs/benchmarks/*.csv is git-ignored). Production stages
  run through `run_batched_julia.py` with `--gage-list`, `--climate-frame`, `--qa-schema`, `--stage`
  and `--require-clean-tree`; a stage's provenance is in its merged timing JSON.
- Dashboards and explorers group by FUNCTION FAMILY (the 17 computing functions); the
  nine signature CATEGORIES come from `docs/signature_categories.csv`. Keep every
  `SIGNATURE_GROUPS` list complete when a family is added — the snow family was silently
  absent from every dashboard before 2026-08-10 (all four lists completed 2026-10-12).
- Basemaps: keyless tiles only. Every Leaflet builder calls `hisssBasemap()` (Esri gray canvas →
  OpenStreetMap fallback → layers control); CARTO's tiles watermark keyless requests since
  2026-10 — never reintroduce them. Check a rebuilt HTML with `check_html_basemap.py` (tiles
  loaded, attribution, screenshot) before staging it.
- Explorer sidecars: `build_signature_explorer.py --annual` writes the per-signature time series to
  `<out-stem>_annual/` AND embeds that folder name in the page (`D.meta.annual.dir`). Renaming the
  HTML for a deposit means renaming the folder and rewriting the embedded name together
  (`stage_hydroshare_climate_v2.py` does; the 2026-10-14 staging defect did not) — then load a
  sidecar headlessly to prove it. A deposit ships the HTML with its complete `_annual/` folder.
- Climate frame (2026-10-26): the three runners build each gage's frame through
  `docs/benchmarks/climate_frame.{jl,py,R}` under `STREAMFLOW_CLIMATE_FRAME = left | outer`
  (`left` = the streamflow dates only, bit-identical to every product before the qualifier
  rerun; `outer` = + the gage's Daymet rows on absent streamflow dates, inside the record
  span and the gage's water years only). Keep the three helpers in step; their fixture
  self-tests are `selftest_climate_frame.{jl,py,R}` (Julia writes the fixture and the
  reference, the ports read it) — run all three after touching a runner's join.
- Procedure, env vars and the post-run checklist: `/run-benchmark`.

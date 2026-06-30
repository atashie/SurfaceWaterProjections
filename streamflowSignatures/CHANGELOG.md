# Changelog

All notable changes to the Streamflow Signatures project.
For full historical detail (Dec 2025 – April 2026), see [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

> **Convention** — the most recent month(s) appear here as a condensed summary with pointers; full per-change detail lives in the archive. When a month's work is complete, condense it here and move the detail to the archive. File-level change lists belong in `git log`, not here; analysis and benchmark tables belong in the canonical docs (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`) and are linked rather than re-hosted.

## [Unreleased]

### Planned
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)
- Port data ingestion utilities to Julia (long-term — currently R-only via dataRetrieval/tidyhydat)
- Generate Julia golden outputs (624 cols, 7,313 gages)
- BFImax estimation via Collischonn & Fan (2013) backward filter — would give BFI_Eckhardt_param per-gage BFImax instead of fixed 0.8, improving discriminating power (currently range [0.47–0.80] due to BFImax saturation)

### Guidelines Document TODOs
<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

---

## [June 2026]

### New: per-watershed MODIS Earth Observation products (LAI + LULC, non-signature)
Two new per-gage EO products aggregated over each gage's upstream watershed from MODIS Collection 061, sitting alongside the signatures and joining by leading-zero-safe `gage_id`/`canon_id`. Python (not a cross-language signature — like HydroATLAS, this is metadata/ingestion). 7,964 watersheds. New subproject `EO_data_processing/` (its own README, pipelines, manifests). All artifacts on S3 `s3://climate-ai-data-science-shiny-app-data/streamflow/`. Full detail: [EO_data_processing/README.md](EO_data_processing/README.md).

- **Watershed geometry layer** (prerequisite) — official basin polygons primary (US GAGES-II + Canada ECCC MDA_ADP), HydroBasins delineation fallback; 7,964 basins (8,018 universe − 54 >100,000 km²). `watershed_polygons_26jun2026.{gpkg,parquet}`. Codex-reviewed.
- **LAI (MCD15A3H, monthly, 2002–2024)** — 270-month panel, 2,150,280 rows (7,964 × 270). Two-stage spatial-distribution-of-monthly-mean (day-weighted per-pixel monthly mean → coverage-weighted basin mean + spatial heterogeneity stats), QA fractions from FparLai_QC/FparExtra_QC. Catalog (GEOLATELY/us-east-2) + LP DAAC backfill for far-N + 2024-11; 17 urban basins permanently NA (MODIS fill code 250). Codex-GO.
- **LULC (MCD12Q1, annual, 2001–2024)** — 191,136 rows (7,964 × 24), 109 cols: per-class % coverage for all 8 classification bands (LC_Type1–5 + LC_Prop1–3, 102 class columns) via coverage-weighted exactextract over VRT-mosaicked sinusoidal tiles. All 8 bands sum to 100 within 1e-13; manifest matches v061 spec. **Codex-reviewed 29 Jun (gpt-5.5) = GO**; 2 latent code hazards (tile-success accounting, unknown-code drop) hardened with the delivered output byte-identical.
- **Cross-cutting**: leading-zero canonical key reconciliation (metadata zero-stripped vs signatures zero-padded); LP DAAC HDF4 read via `pyhdf` (pip-GDAL lacks the driver); Earthdata account-lock handled by completeness-guard + retry/backoff (self-healing, resumable).
- **LAI `good_coverage_frac`** (30 Jun) — per basin-month continuous QA flag = valid pixel-obs / (all basin pixels × all expected composite dates); generalizes the binary `partial_month`. Added to `watershed_modis_lai_monthly_30jun2026.parquet` (data otherwise unchanged, hash-verified).
- **LAI QA/QC explorer** (30 Jun) — self-contained Leaflet HTML (`EO_data_processing/viz/build_lai_explorer.py` → `watershed_modis_lai_explorer.html`) for manual review: 7,964 points colored by any of 17 QA/summary variables, click a gage for its full monthly LAI time series (quantile bands + per-month good-coverage strip), map/TS hovers, glossary, processing notes. Mirrors `streamflow_explorer.html`; LULC counterpart next.

## [May 2026]

### New: per-gage watershed HydroATLAS metadata (non-signature product)
A standalone metadata file pairing each gage with the hydro-geophysical character of its **entire upstream watershed**, from HydroATLAS / BasinATLAS v10 (281 source attributes → ~211 output columns: climate, hydrology, terrain, land cover, soils/geology, anthropogenic). Sits alongside the signatures, keyed by leading-zero-safe `gage_id`: `data_out/watershed_hydroatlas_metadata_{date}.{csv,parquet}` + a generated data dictionary. 8,014 gages, ~13 MB.

- **Hybrid aggregation** over each gage's delineated upstream basin set (`upstream_hydrobasins.rds`): HydroATLAS upstream-accumulated `_u` / pour-point `_p` values pass through from the outlet basin (91 attrs — the delineated set equals HydroATLAS's own upstream extent, so `_u` is the rigorous watershed value); `_s`-only attributes are SUB_AREA-weighted (92 area-weighted means incl. monthly climate, with the HydroATLAS `-9999` NoData sentinel masked; elevation as spatial min/max; 11 categoricals as area-weighted majority — glc/pnv via argmax of upstream `_u` fractions, wet + the other 8 via area-weighted mode of the source class). Percentage / areal-extent fields are always area-weighted mean (never max). A `watershed_area_rel_diff` column ( |summed member SUB_AREA − gage `basin_area`| / `basin_area`, gage-reported drainage as truth) flags gages where the delineated watershed area diverges from the reported drainage (mis-snapped outlets).
- New `R/aggregate_hydroatlas_metadata.R` + `run_hydroatlas_metadata.R`; `HYDROATLAS_GPKG` added to `config.R`. R-only (metadata/ingestion) — not a cross-language signature, so no Julia/Python port.
- Validated: member SUB_AREA sums to outlet UP_AREA (median rel.diff 1e-4); area-weighted `_s` reproduces HydroATLAS `_u` (elevation cor 0.9994); 100% join to the signatures file (5,707/5,707). See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Watershed HydroATLAS Metadata.
- **Independent code review** (in-session) confirmed the classifier, NA-aware weighting, BFS delineation, and leading-zero joins; correctness fixes applied — `wet_cl_smj` switched from `_u`-fraction argmax to area-weighted mode (its fraction basis omits classes 10/11), `-9999` NoData masked in weighted means (`sgr_dk_sav`), and a malformed legacy-cache member id sanitized.

### New: static HTML data explorer
A self-contained, double-clickable `data_out/streamflow_explorer.html` (Leaflet + Canvas) for exploring the gages: 8,014 points colored by any of **412 variables** (dataset switcher between HydroATLAS watershed metadata and streamflow signatures); click a gage to draw its **watershed boundary** (border-only). Generators: `build_explorer.R` (unions the RAW lev12 basins per outlet for a clean dissolve — no interior sliver lines — then strips holes, simplifies, and caps at ~250 vertices → `watershed_borders.geojson` ~9.5 MB), `assemble_explorer.R` (joins points + injects into `explorer_template.html`). Full-resolution borders are ~111 MB (too heavy to inline); the clean-dissolve build is ~28 MB total. Border click-through and clean rendering verified headlessly (puppeteer-core + Chrome). See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Static HTML Explorer.

---

## [April 2026]

A large release: new signature families, the canonical-language transition, centralized NA handling, and cross-language alignment to 99.5% agreement. Output grew from 551 to **1,264 columns** (656 base/metadata + 608 Pettitt changepoint) across 7,313 gages. **Full per-change detail is in [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md) → [April 2026].**

### New signatures
- **Pettitt changepoint detection** applied to all 76 time-series signatures (8 fields each → 608 columns: cp year, p-value, pre/post means, delta, pct change, pre/post MK p-values). ~13.4% significance rate (2.7× null); signal-robustness analysis in [docs/SIGNATURES.md](docs/SIGNATURES.md) → Changepoint Detection.
- **Recession-parameterized baseflow** — `alpha_linear`, `BFI_Eckhardt_param`, `BFI_LyneHollick_param`, plus the per-gage `recession_alpha_point_cloud_linear_reservoir` scalar (Collischonn & Fan 2013). +25 columns.
- **Guidelines Section 3 additions** — D1/D99 flow timing, `n_recession_events`, `elasticity_annual`, `negative_ann`, `runoff_ratio_high_count`, 4 seasonal exclusion-year diagnostics, and 2 elasticity diagnostics.

### Canonical language transition
- **Julia replaced R as the canonical implementation** (~40× faster, modular). `rpkg/` is the production R port; monolithic `R/helperFunctions.R` is deprecated (ingestion utilities only). See [CLAUDE.md](CLAUDE.md).

### Centralized NA handling
- New `preprocess_daily_data()` runs once per gage before signatures: daily-grid normalization, ≤3-day interpolation, year rejection, seasonal completeness flags. Removed all per-signature min-days thresholds and `fillna(0)` calls; added the 80% trend-completeness gate to `generate_stats()`. Negative-Q rejection is config-driven (off by default); constant-SD is a flag only. See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → NA Handling Architecture.

### Cross-language alignment
- Synced Section 3 plus the recession-algorithm rewrite to Python and rpkg; fixed Julia's Mann-Kendall tau-b, a `leftjoin` row-order bug, and several R-canonical divergences. The three synced implementations (Julia, Python, rpkg) agree at 99.5% R² (615/623 columns Perfect). Benchmark tables: [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Cross-Language Benchmarks.

### Tooling & audits
- Julia sensitivity-experiment framework (WY≥1993 ± qualifying-fraction filters), a Julia-vs-golden comparison pipeline, and interactive HTML validation dashboards (`docs/benchmarks/`).
- Three audits: Guidelines-vs-implementation, NA-handling (Codex), and a Julia adversarial review — all findings since addressed or documented.

---

## [March 2026]

### rpkg: Proper R Package
New R package (`rpkg/`) rewriting the monolithic `R/helperFunctions.R` as a proper installable package. 16 R source files, 6 test files, 27 man pages. R CMD check: 0 errors, 0 warnings. See [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md) for full details.

### Cross-Language Alignment (Rounds 1–6, complete)
Six rounds of alignment work brought R, Python, and Julia implementations to production-ready agreement. Final status (March 16-17, 2026):

| Pair | Mean R² | Median R² | Cols < 0.99 |
|------|---------|-----------|-------------|
| R vs Python | 0.9968 | 1.0000 | 17 |
| R vs Julia | 0.9965 | 1.0000 | 19 |
| Python vs Julia | 0.9997 | 1.0000 | 3 |

- 475 perfect (R²>=0.999), 56 good (0.99-0.999), 20 poor (<0.99)
- 531 of 551 columns (96.4%) have R² >= 0.99 across all 3 pairs
- Golden output regression: 551/551 perfect match against Feb 2026 golden outputs
- 20 remaining poor columns are all trend statistics (slopes, p-values) — underlying means/medians are near-identical
- Full round-by-round details: [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md), [docs/CROSS_LANGUAGE_STATUS.md](docs/CROSS_LANGUAGE_STATUS.md)

### Code Quality Improvements (Rounds 1–2)
Cross-language code review covering ~9,000+ lines. Key fixes: O(n²) metadata lookup → O(1), Python vectorized water year calc (~50x speedup), Julia DataFrame pre-allocation, Eckhardt BFI forward-fill aligned across all 3 languages. Full details in `docs/CODE_REVIEW.md` and [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

### Benchmark Timing (March 16-17, 2026)

| Implementation | Time | Gages | Rate |
|---------------|------|-------|------|
| Julia | 9.2 min | 7,369 | 13.4/s |
| Python | 78.9 min | 7,369 | 1.56/s |
| rpkg | 114 min | 7,369 | 1.08/s |
| R (canonical) | 874 min* | 5,707 | 0.11/s |

*R ran concurrently with Python/Julia — timing inflated by I/O contention.

---

## Version History Notes

This project uses date-based versioning (MONTH YEAR) rather than semantic versioning, reflecting its nature as a research tool with continuous development.

### Output File Naming Convention
Output files include date stamps: `streamflow_signatures_full_JAN2026.csv`

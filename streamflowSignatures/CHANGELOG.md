# Changelog

All notable changes to the Streamflow Signatures project.
For full historical detail (Dec 2025 – April 2026), see [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

> **Convention** — the most recent month(s) appear here as a condensed summary with pointers; full per-change detail lives in the archive. When a month's work is complete, condense it here and move the detail to the archive. File-level change lists belong in `git log`, not here; analysis and benchmark tables belong in the canonical docs (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`) and are linked rather than re-hosted.

## [Unreleased]

### Planned
- Port annual-values collector, b=1 recession alpha, AND snow metrics to Python and
  rpkg (after the Julia benchmark re-run validates all three)
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)
- Port data ingestion utilities to Julia (long-term — currently R-only via dataRetrieval/tidyhydat)
- Generate Julia golden outputs (624 cols, 7,313 gages)
- BFImax estimation via Collischonn & Fan (2013) backward filter — would give BFI_Eckhardt_param per-gage BFImax instead of fixed 0.8, improving discriminating power (currently range [0.47–0.80] due to BFImax saturation)

### Known Issues (discovered 2026-07-14, Codex review — not yet fixed)
- **MEDIUM (documented limitation, by design) — 37 Canadian gages in the signature
  output carry raw m³/s units (`area_normalized = FALSE`)**: HYDAT publishes NO
  drainage area (neither `DRAINAGE_AREA_GROSS` nor `DRAINAGE_AREA_EFFECT` — verified
  directly against `Hydat.sqlite3`) for 73 successfully-processed Canadian stations;
  37 survive the 20-year filter into the Julia canonical output (7,313 gages).
  Station names show most are not natural watersheds: ~40 irrigation/diversion canals
  + ditches, ~15 dam/powerhouse outflows, several huge-river channel splits and lake
  outlets (St. Lawrence, Mackenzie, Nelson, Lake of the Woods); only ~8 look like
  natural streams. 62/73 are `REGULATED = TRUE`.
  **DECISION (user, July 2026): keep these gages with raw m³/s flow — NO area
  backfill** (HydroBasins `UP_AREA` was assessed on 1,383 validation gages: accurate
  only for main-stem dam outflows, wrong for canals and channel splits, unusable
  <100 km²). Q-to-PPT signatures are now structurally gated for these gages (see the
  July 2026 entry below). **Remaining limitation**: unit-carrying Q-only signatures
  (Q volumes, percentiles, Q95_Q10, log_a) stay in m³/s for these 37 rows —
  incomparable with mm/day gages (Qann_mean up to 3.18M for the St. Lawrence).
  **Flag gap**: `flagged_for_qann_range` catches only 27/37 — 10 small canals/creeks
  land inside [0, 2000] unflagged; downstream users must filter on
  `area_normalized == TRUE` before any cross-gage comparison of unit-carrying
  signatures. See docs/DEVELOPMENT.md → Canadian HYDAT → Missing drainage areas.
- **MEDIUM — seasonal runoff ratios ignore seasonal completeness flags**:
  `julia/src/runoff_ratios.jl:113-118` looks up flags named
  `winter_complete/spring_complete/summer_complete/fall_complete`, but the preprocessor
  emits `win_complete/spr_complete/sum_complete/fal_complete` (`julia/src/io.jl:462`;
  cf. correct usage in `flow_volumes.jl:136-137`). The existence check at
  `runoff_ratios.jl:120` silently fails, so seasonal runoff ratios are never NaN'd for
  incomplete seasons — the guidelines' 80% seasonal completeness rule is not applied to
  this category. Python and rpkg likely mirror the bug (synced ports — verify). Fix
  changes outputs → schedule with the next behavior-changing release (Julia first, then
  ports, benchmark re-run + golden comparison).

### Guidelines Document TODOs
<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

---

## [July 2026]

### New: Q-to-PPT unit gate for un-normalized gages (Julia + Python + rpkg)
Gages with no drainage area in HYDAT keep flow in raw m³/s (decision: no backfill).
`calculate_all_signatures()` gains an `area_normalized = true` kwarg in all three
languages: when `false`, ALL Q-to-PPT signatures — runoff ratios (+
`runoff_ratio_high_count`), elasticity (+ diagnostics), Q-P seasonality, and
`avg_storage` — are skipped, because Q (m³/s) and PPT (mm) units don't match, making
Q/P, dQ/dP, and cumsum(P − Q) meaningless. Q-only signatures (incl. Q##/D##
percentiles, recession, BFI) are unaffected. All three benchmark runners look up
`area_normalized` from the metadata CSV (leading-zero-safe join; only an explicit
FALSE gates, missing defaults to normalized) and pass it per gage.

- **Intentional output change at the next benchmark**: 16 of the 37 un-normalized
  gages currently carry mixed-unit climate signatures in the April 2026 golden
  (e.g. runoff ratios computed from m³/s Q against mm PPT) — those become NA.
- Elasticity and `qp_bimodality` are technically scale-invariant in Q, but are gated
  with the rest: mixed-unit inputs are conceptually invalid, and the family's
  internal PPT thresholds and P − Q terms are not scale-invariant.
- Tests (all three languages verify: 84 climate keys normally → 0 gated; gated
  output IDENTICAL to a no-climate run; gate inert when true):
  `julia/test/test_area_normalized_gate.jl` (also covers the AnnualCollector path),
  `python/tests/test_area_normalized_gate.py`,
  `rpkg/tests/testthat/test-area_normalized_gate.R`.
- **Codex adversarial review (2026-07-14)**: 2 MEDIUM findings, both fixed — (1) the
  Julia/Python runners hard-failed on metadata files predating the `area_normalized`
  column (now header-guarded, matching rpkg's existing guard); (2) explicit-FALSE
  parsing diverged across languages for 0/1-serialized booleans (now harmonized in
  all three runners: Bool false / numeric 0 / "FALSE"/"0" strings gate; anything
  else, incl. missing, defaults to normalized). Fix verification surfaced a third
  latent bug the review missed: the Julia runner's join-id helper was typed
  `::String` but CSV.jl yields InlineStrings — Phase 2b would have crashed on the
  production metadata (loosened to `AbstractString`). Runner parsing was then
  functionally verified in all 3 languages against TRUE/FALSE, 0/1, and
  missing-column metadata variants, plus the production CSV (16,994 gages parsed,
  exactly the expected 1,601 gated).

### New: snow metrics signature family (14 metrics, Daymet SWE — Julia)
Fourteen per-water-year snow metrics from daily SWE, through the standard
`generate_stats()` path (8 stats + 8 Pettitt fields each → **+224 columns**; 76 → 90
time-series signatures; annual values exported automatically via the collector). New
module `julia/src/snow.jl` (`calculate_snow_metrics`); full spec + review record:
`docs/plans/snow_signatures_plan.md`. Requested by the user with four in-session domain
decisions (2026-07-14).

- **Metrics**: `swe_max`, `swe_max_dowy`, `snow_cover_days`, `snow_on_dowy`,
  `snow_off_dowy`, `melt_season_days`, `melt_rate`, `ssm` (Hatchett 2021, Hydrology
  8(1):32), `swe_apr1`, `melt_before_peak`, `melt_before_peak_pct`,
  `melt_before_peak_to_max_swe`, `melt_com_dowy`, `swe_max_to_ppt`.
- **Domain decisions**: SWE ≥ 10 mm day threshold applied to durations AND magnitudes
  (thresholded series SWE*; sub-threshold years are operationally snow-free);
  snow-on/off anchored to the spell containing the annual peak (boundary-censored →
  NaN); melt rate = max SWE ÷ melt-season days; SSM seasonal spells ≥ 60 continuous
  days. Config: `signatures_config.json` → `snow` + `na_handling.snow_na_policy`.
- **Plumbing**: Daymet `swe` was previously dropped at the climate join — now
  normalized (`swe`→`SWE`) and carried through; `preprocess_daily_data()` gained a
  per-year SWE policy (mirrors PPT: >30 NAs, ≤3-day interpolation, negative rejection)
  and a new `valid_swe_years` return field. Snow metrics run ONLY on an explicitly
  passed `snow_data` frame filtered to those years — NO implicit fallback to the gage
  frame (plan review finding: prevents SWE-invalid years leaking in). Canadian gages
  (no Daymet) → all snow columns NA.
- **Known existing-column change (accepted)**: `flagged_for_high_na` counts the 224 new
  NA fields in its denominator, so some no-SWE gages will flip at the re-run
  (April 2026 Pettitt precedent; whitelisted + quantified; semantics pinned by test).
- **Two-round Codex plan review**: initial NO-GO (2 CRITICAL: orchestrator fallback
  leak; the false "no golden divergence" claim) → 7 findings + 2 delta-round residuals
  incorporated → GO. Record: plan doc §12.
- **Tests**: `julia/test/test_snow_metrics.jl` (~330 assertions: exact hand-derived
  gages — triangle, mid-winter thaw, spell arithmetic, threshold splitting, boundary
  censoring, tie rule, leap-year April 1, PPT gating incl. legacy runoff-ratio parity,
  preprocessor `valid_swe_years`, no-fallback + 0-row schema contract, SWE-invalid-year
  leak regression, collector/zero-warn invariants, QA-flag semantics pin);
  `smoke_test.jl` now joins SWE and validates snow presence/ranges.
- Julia only; Python/rpkg ports deferred (three-feature port queue). Benchmark re-run
  pending — bundled with annual values + b=1 alpha.

### Docs: data-source inventory
New `docs/DATA_SOURCES.md` — table of all 11 external data sources (data class /
project use / provider / original access): USGS NWIS, HYDAT, Caravan, Daymet,
GAGES-II boundaries + attribute tables, ECCC MDA_ADP basin polygons, HydroBASINS,
HydroATLAS/BasinATLAS v10, MODIS MCD15A3H (LAI) and MCD12Q1 (LULC). Also records
the guidelines Google Doc as a governance input and the deliberately excluded
boundary cases (planned ERA5/PRISM, S3 mirrors, viz-only CDNs, unused geometry
backstops).

### Changed: recession alpha now assumes a linear reservoir (b = 1) — Julia
All alpha outputs (`log_a_pointcloud`, `log_a_events`, the 6 `log_a_seasonality_*`
scalars) are now computed with the recession exponent **fixed at b = 1** across all
locations and periods: `log(a) = median(log(-dQ/dt) - log(Q))`, no regression. The
exponent analyses are unchanged — `b_pointcloud`, `b_events`, and `concavity` keep
their free power-law fits, as do event identification, the 25-event gate, and
`alpha_linear` / `recession_alpha_point_cloud_linear_reservoir` (already b=1).

- **Rationale** (domain decision): log(a) is the intercept of a regression whose slope
  is b, so free-fit alpha estimates are convolved with b — trends/seasonality in alpha
  partly reflected changes in b. Fixing b = 1 decouples them.
- Column names unchanged (frozen CSV contract) — methodology change behind existing
  columns, like the April 2026 recession rewrite. At the next benchmark, log_a columns
  will diverge from golden intentionally; b/concavity/alpha_linear must not.
- Consequence: under the forward-difference discretization each b=1 point equals
  `log(1 − Q_{i+1}/Q_i)`, so per-year `log_a_pointcloud` APPROXIMATES
  `log(1 − alpha_linear)`. The relation is approximate, not a same-sample identity:
  `alpha_linear` also includes events whose free power-law fit failed (the point cloud
  only pools fit-success events), and medians only commute exactly with the monotone
  transform at odd pair counts.
- New tests (`julia/test/test_recession_alpha_b1.jl`): exact helper values on a pure
  exponential; a quadratic-reservoir gage (−dQ = k·Q² exactly) proving b stays 2 under
  the free fit while alpha matches an independent b=1 recompute (and is far from the
  old free-fit intercept); linear-reservoir continuity (b=1 gage values unchanged).
- Julia only; Python/rpkg sync deferred (bundled with the annual-values port).
  Design note: `docs/plans/recession_alpha_linear_reservoir_plan.md`.

### New: annual values export (per-year signature values, Julia)
The per-year annualized metric values — previously computed inside every signature
function and discarded after `generate_stats()` collapsed them into the 8 summary
statistics — are now saved as one long-format parquet alongside the summary CSV:
`{prefix}_signatures_annual.parquet` with columns `gage_id, signature, water_year, value`
(~20M rows expected for 7,313 gages × ~76 signatures, ~170 MB zstd).

- **Mechanism**: opt-in `AnnualCollector` threaded through `calculate_all_signatures()`
  and all 13 signature functions into `generate_stats()` (mirrors the `changepoint`
  kwarg plumbing). Collection happens BEFORE min_rows / trend-completeness gating and is
  read-only w.r.t. the statistics — the summary CSV contract is untouched (verified by a
  with/without-collector equality test).
- **Semantics**: the export records the exact series `generate_stats()` consumed. Dense
  signatures carry NaN placeholder rows ("year qualified, metric not computable");
  caller-pruned signatures (`qp_bimodality`, `n_recession_events`, flashiness, storage,
  baseflow, elasticity) omit those rows — absent row ≡ NaN for consumers.
  `elasticity_rolling` is keyed to the END year of each 11-yr window; `elasticity_annual`
  to the later year of each consecutive pair.
- **Config**: `config/signatures_config.json` → `annual_values.save` (shipped `true`;
  defaults to `false` when the section is absent). Constant `CFG_SAVE_ANNUAL_VALUES`.
- **Validation**: new tests in `julia/test/test_annual_collector.jl` (stat-identity,
  zero-warnings, signature coverage, NaN/missing-year guards, no duplicate keys) +
  `docs/benchmarks/validate_annual_values.py` (recomputes mean/median from the parquet
  and cross-checks the summary CSV; runs after the next benchmark).
- **Pre-flight**: Parquet2 pinned (`0.2`, resolved v0.2.27) after a write→read round-trip
  check (zstd/snappy/uncompressed, NaN-safe, empty-frame-safe; 1M rows = 8.5 MB, 0.2s).
- **Two Codex reviews**: plan-stage (NO-GO → 5 amendments incorporated pre-implementation)
  and code-stage (core implementation confirmed correct; 3 harness findings fixed:
  validator ≥3-non-NaN-rows coverage rule, zero-gage no-op exit, and a deterministic
  exponential-recession test gage covering the recession + parameterized-BFI collector
  paths the noisy synthetic gage never exercised). Full records:
  `docs/plans/annual_values_export_plan.md` §9–§10.
- Julia only — Python/rpkg ports deferred.
- **Benchmark re-run pending** — bundled with an upcoming workflow feature.

## [June 2026]

### New: per-watershed MODIS Earth Observation products (LAI + LULC, non-signature)
Two new per-gage EO products aggregated over each gage's upstream watershed from MODIS Collection 061, sitting alongside the signatures and joining by leading-zero-safe `gage_id`/`canon_id`. Python (not a cross-language signature — like HydroATLAS, this is metadata/ingestion). 7,964 watersheds. New subproject `EO_data_processing/` (its own README, pipelines, manifests). All artifacts on S3 `s3://climate-ai-data-science-shiny-app-data/streamflow/`. Full detail: [EO_data_processing/README.md](EO_data_processing/README.md).

- **Watershed geometry layer** (prerequisite) — official basin polygons primary (US GAGES-II + Canada ECCC MDA_ADP), HydroBasins delineation fallback; 7,964 basins (8,018 universe − 54 >100,000 km²). `watershed_polygons_26jun2026.{gpkg,parquet}`. Codex-reviewed.
- **LAI (MCD15A3H, monthly, 2002–2024)** — 270-month panel, 2,150,280 rows (7,964 × 270). Two-stage spatial-distribution-of-monthly-mean (day-weighted per-pixel monthly mean → coverage-weighted basin mean + spatial heterogeneity stats), QA fractions from FparLai_QC/FparExtra_QC. Catalog (GEOLATELY/us-east-2) + LP DAAC backfill for far-N + 2024-11; 17 urban basins permanently NA (MODIS fill code 250). Codex-GO.
- **LULC (MCD12Q1, annual, 2001–2024)** — 191,136 rows (7,964 × 24), 109 cols: per-class % coverage for all 8 classification bands (LC_Type1–5 + LC_Prop1–3, 102 class columns) via coverage-weighted exactextract over VRT-mosaicked sinusoidal tiles. All 8 bands sum to 100 within 1e-13; manifest matches v061 spec. **Codex-reviewed 29 Jun (gpt-5.5) = GO**; 2 latent code hazards (tile-success accounting, unknown-code drop) hardened with the delivered output byte-identical.
- **Cross-cutting**: leading-zero canonical key reconciliation (metadata zero-stripped vs signatures zero-padded); LP DAAC HDF4 read via `pyhdf` (pip-GDAL lacks the driver); Earthdata account-lock handled by completeness-guard + retry/backoff (self-healing, resumable).
- **LAI `good_coverage_frac`** (30 Jun) — per basin-month continuous QA flag = valid pixel-obs / (all basin pixels × all expected composite dates); generalizes the binary `partial_month`. Added to `watershed_modis_lai_monthly_30jun2026.parquet` (data otherwise unchanged, hash-verified).
- **LAI QA/QC explorer** (30 Jun) — self-contained Leaflet HTML (`EO_data_processing/viz/build_lai_explorer.py` → `watershed_modis_lai_explorer.html`) for manual review: 7,964 points colored by any of 17 QA/summary variables, click a gage for its full monthly LAI time series (quantile bands + per-month good-coverage strip), map/TS hovers, glossary, processing notes. Mirrors `streamflow_explorer.html`.
- **LULC QA/QC explorer** (30 Jun) — self-contained Leaflet HTML (`EO_data_processing/viz/build_lulc_explorer.py` → `watershed_modis_lulc_explorer.html`, ~49 MB, on S3). 7,964 points colored by any of **27 map variables** in two non-overlapping groups: 10 IGBP-derived summary features (`class_diversity`, `delta_forest/urban/cropland` [robust mean(2001–03) vs mean(2022–24) endpoints], `dominant_class`, `dominant_group` [lump-robust aggregate dominance — fixes IGBP's class-splitting bias], `dominant_change`, `n_modis_pixels`, geom source / low-confidence) + all 17 individual IGBP classes (static `pct_*` roll-ups dropped as redundant with the per-class maps). Click a gage → stacked-area land-cover composition 2001–2024 with a **band switcher across all 8 MODIS schemes** (IGBP/UMD/LAI-fPAR/BGC/PFT/FAO-LCCS×3), official MODIS IGBP palette, per-year class-breakdown hover (0.1 % resolution). Glossary carries the MODIS interpretation caveats (SE-US woody-savanna behavior; 2021+ reduced confidence). Headless-validated (puppeteer/Chrome). **EO pipeline + viz COMPLETE.**

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

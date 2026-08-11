# Streamflow Signatures Project

## Project Context

This project has two primary goals and two secondary goals:

1. **Data Processing** — Ingest raw streamflow data (USGS, HYDAT, Caravan), clean/filter/collate metadata, standardize outputs.
2. **Signature Extraction** — Extract 100+ hydrological signatures under strict guardrails. Domain experts update methodology via plain-English instructions in `docs/SIGNATURE_GUIDELINES.md`; code implements those definitions.
3. **Visualization** (secondary) — Shiny dashboard for broader audience exploration.
4. **Cross-Language Implementations** (secondary) — Julia is canonical; Python and rpkg produce near-identical results for community sharing (future: publishable packages/libraries).

## Multi-Language Architecture

This project provides identical signature calculations in Julia (canonical), Python, and R:

| Directory | Language | Status | Description |
|-----------|----------|--------|-------------|
| `julia/` | Julia | **Canonical** | Reference implementation - all changes start here |
| `python/` | Python | Active | Port of Julia signatures |
| `rpkg/` | R | Active | R port of Julia signatures |
| `R/` | R | Deprecated | Legacy shim (still functional for ingestion) |

**Change Workflow**: Julia is canonical. Changes are made in Julia first, then propagated to Python and rpkg. Golden outputs from Julia (April 2026) validate other implementations. Historical note: R was canonical through March 2026; Guidelines Section 3 changes (April 2026) were implemented Julia-first for faster iteration (~10 min benchmark vs hours for R), which drove the transition. Python and rpkg synced April 14-15.

**Canadian HYDAT Metadata**: RHBN and REGULATED status for Canadian gages is pre-exported to `metadata/canadian_hydat_interference.csv` (via `R/export_hydat_metadata.R` using tidyhydat). Julia reads this CSV directly; R uses tidyhydat at runtime.

## Canonical Code

**Julia canonical source:**
- `julia/src/` - All canonical signature modules (17 files under `StreamflowSignatures.jl`)
- `config/signatures_config.json` - Cross-language configuration (source of truth)
- `config.R` - R-side configuration, logging, validation (still active for R ingestion scripts)

**R port:**
- `rpkg/` - Proper R package mirroring Julia structure (active, production-ready)

**Deprecated:**
- `R/helperFunctions.R` - Legacy monolithic R implementation (deprecated; use rpkg instead)
- Other `helperFunctions*.R` files in `archive/` are deprecated.

## Key Entry Points

- `julia/src/StreamflowSignatures.jl` - PRIMARY: Julia package module (canonical)
- `docs/benchmarks/run_julia_benchmark.jl` - Full signature extraction benchmark (~27 min, 7,313 gages)
- `run_ingest_usgs_hydat.R` - Raw USGS/HYDAT data ingestion to parquet (R, still active)
- `run_full_processing.R` - Legacy: Full R signature extraction with climate (still functional for ingestion)
- `run_caravan_processing.R` - Caravan data processing (lower priority)
- `streamflowAndClimateVisualizationApp/app.R` - Shiny dashboard

## Critical Constraints

1. **CSV Output Format**: MUST remain unchanged - downstream tools depend on exact column names
2. **Water Year**: Oct 1 - Sep 30
3. **Flow Units**: mm/day (converted from cfs/m3s)
4. **Minimum Data Requirements**:
   - 20+ water years per gage
   - Year qualification via `preprocess_daily_data()`: rejects years with >30 NAs, >3-day gaps; negative Q rejection conditional on `reject_negative_flow` config (default: false)
   - No per-signature min_days thresholds — preprocessor is single source of truth
   - No additional min_Q_value_and_days filter in non-legacy path (removed April 2026 — was causing R to diverge from Python/Julia by excluding low-flow years)
5. **Experiment outputs live together** (user convention, 2026-07-22): ALL
   artifacts of a production/experiment run — signatures CSV, annual parquet,
   timing JSON, run log, signature explorer (+ its `_annual/` sidecar folder),
   and every comparison dashboard/CSV/summary generated for that run — belong in
   that experiment's OWN output folder, one unique folder per experiment (e.g.
   `processedOuts_drought_28jul2026` for the current WY 1993–2025 standard
   product). Do NOT leave
   run-specific artifacts in `docs/benchmarks/` — the repo keeps only the tools
   and the long-lived cross-language reference CSVs (April golden/experiment
   files).
6. **NA Handling** (April 2026): `preprocess_daily_data()` runs ONCE per gage BEFORE signatures.
   - Interpolates internal gaps <= 3 days; rejects years with >30 raw NAs, >3-day gaps
   - Negative Q rejection is config-driven (`reject_negative_flow: false` by default); `negative_ann` signature counts Q<0 days instead
   - Constant-SD is a QA flag only (never causes year rejection)
   - Config: `config/signatures_config.json` → `na_handling` section
   - `use_legacy_filtering: false` — new preprocessing is the default
   - Do NOT use fillna(0) in signature functions — the preprocessor handles NAs centrally
   - Do NOT add per-year min_days or max_na_frac checks in signature functions

## Signature Statistics Rule

**Every signature MUST produce exactly 8 statistics using `generate_stats()`:**

| Suffix | Statistic |
|--------|-----------|
| `_senn_slp` | Theil-Sen slope |
| `_linear_slp` | Linear regression slope |
| `_spearman_rho` | Spearman correlation |
| `_spearman_pval` | Spearman p-value |
| `_mk_rho` | Mann-Kendall tau |
| `_mk_pval` | Mann-Kendall p-value |
| `_mean` | Arithmetic mean |
| `_median` | Median |

**Exceptions** (documented in `config.R`):
- `elasticity_static` - single value, not time series
- `log_a_seasonality_amplitude`, `log_a_seasonality_minimum` - recession seasonality
- `runoff_ratio_high_count` - per-gage scalar (count of years with ratio > 2.0)
- `elasticity_years_total`, `elasticity_years_low_ppt` - per-gage diagnostics
- `ice_affected_days_total` - per-gage diagnostic
- `recession_alpha_point_cloud_linear_reservoir` - per-gage scalar (whole-record median Q_{i+1}/Q_i)

## Code Status

| File | Status | Notes |
|------|--------|-------|
| `julia/src/*.jl` | **CANONICAL** | All canonical signature modules (17 files) |
| `docs/benchmarks/run_julia_benchmark.jl` | **PRIMARY** | Full benchmark entry point |
| `config/signatures_config.json` | **ACTIVE** | Cross-language configuration |
| `python/streamflow_signatures/` | **ACTIVE** | Python port |
| `rpkg/` | **ACTIVE** | R port (proper package) |
| `config.R` | **ACTIVE** | R-side config (ingestion scripts) |
| `run_ingest_usgs_hydat.R` | **ACTIVE** | Raw data ingestion (R) |
| `R/helperFunctions.R` | **DEPRECATED** | Legacy shim - use rpkg instead |
| `run_full_processing.R` | **LEGACY** | Still functional for R ingestion |
| `R/tests/smoke_test.R` | **ACTIVE** | Quick R validation |
| `R/tests/qa_qc_signatures.R` | **ACTIVE** | Output validation |
| `archive/*` | **DO NOT USE** | Reference only |

## Adding New Signatures

1. Create function in appropriate `julia/src/*.jl` module returning annual values
2. Call `generate_stats()` to produce 8 statistics
3. Register in `julia/src/signatures.jl` orchestration function
4. Add base name to `EXPECTED_SIGNATURE_BASES` in `config.R` (per-gage scalars that
   don't follow the 8-stat pattern need their own `EXPECTED_*` constant, wired into
   `validate_output_schema`)
5. **Register in the test/validation registries** — easy to miss, and each one fails
   loudly only after the fact:
   - `EXPECTED_DENSE_SIGNATURES` in `julia/test/test_annual_collector.jl` — asserts
     **set equality** of collected annual series, so any new dense signature fails it
     until listed
   - the signature-count gate in `docs/benchmarks/validate_production_run.py`
     (`ann.signature.nunique() == N`)
   - the **expected total summary-column count** in the docs — count both the
     8-stat + 8-Pettitt fields AND any non-8-stat scalars (getting this wrong is easy:
     the drought family shipped documented as +160 when it is +165)
6. Run the Julia unit suite (`julia --project=julia julia/test/runtests.jl`), then the
   benchmark (`docs/benchmarks/run_julia_benchmark.jl`, ~27 min) to verify
7. **Prove additivity, don't assume it** — diff the benchmark output against the previous
   canonical run: every pre-existing column must be unchanged (only `flagged_for_high_na`
   is expected to shift, since its denominator counts all fields), AND the new columns
   must be populated. The orchestrator's per-signature `try/catch` turns an unexpected
   failure into silently missing columns, so a green unit suite proves nothing here.
   Smoke tests should assert new values are FINITE, not merely that the keys exist.
8. Port to Python (`python/streamflow_signatures/`) and rpkg (`rpkg/R/`)

## Current Standard Products (as of 2026-08-11)

Two delivered products, both @ 60% qualifying fraction, both **1,653 columns** (incl. the
drought family), each in its own folder with explorer + comparison dashboards + retained
validation reports:

| # | Window | Folder | Gages | Annual parquet |
|---|---|---|---|---|
| 1 | WY 1993–2025 | `processedOuts_drought_28jul2026` | 6,678 | 18,898,406 rows / 100 signatures |
| 2 | WY 1980–2025 | `processedOuts_1980_2025_11aug2026` | 6,250 | 24,366,487 rows / 100 signatures |

Both supersede the 22 Jul folders (1,488 columns, no drought). **Neither is a subset of the
other** and **record-dependent signatures (drought thresholds, `*_all` pulses, elasticity,
parameterized BFI) must never be compared across them — nor re-aggregated from the annual
values onto a different window, since their thresholds/record means come from the run's own
window.** Each Julia run also writes the per-signature annual values parquet alongside the
summary CSV (when `annual_values.save` is on — the shipped default; not implemented in the
Python/rpkg ports) — see DEVELOPMENT.md → Annual Values Export.

⚠️ **Climate input**: the canonical `daymet_1980_2023.parquet` is TRUNCATED; use
`daymet_1980_2023_rebuilt_10aug2026.parquet` (product #1 predates the rebuild, product #2
uses it; difference bounded at ≤ 3.4e-13). See DEVELOPMENT.md → Active Parquet Files.

## References

- **@docs/DEVELOPMENT.md** - Architecture, file structure, common tasks, workflows
- **@docs/SIGNATURES.md** - Detailed signature documentation (11 categories)
- **@CHANGELOG.md** - Current work, roadmap (historical: `docs/CHANGELOG_ARCHIVE.md`)
- **@docs/SIGNATURE_GUIDELINES.md** - Collaborative guidelines from hydrology colleagues (auto-synced)
- **@docs/MANUSCRIPT_DRAFT.md** - HISSS manuscript draft snapshot (auto-synced; reconciliation review)
- **@EO_data_processing/README.md** - Earth Observation (MODIS LULC & LAI) per-watershed ingestion and processing
- **@EO_data_processing/README_NLCD.md** - Per-watershed Annual NLCD (CONUS land cover + impervious, 30m, 1985-2025) — MODIS LULC sibling product

## Session-Start Workflow: Document Sync & Reconciliation

Two collaborative Google Docs govern this project and are synced **at the start of each session**. Fetch technique for both: WebFetch paraphrases through a small model, so for faithful sync/diffing download the published page's raw HTML (curl) and extract the text from the `<div id="contents">` block.

### A. Signature Guidelines (ground truth for methodology)

The guidelines document is a core design feature: domain experts write plain-English signature definitions and QA/QC requirements, and those are translated into code. The current doc (created July 2026) supplanted the previous version at the same publish URL.

1. **Fetch fresh content** from the Google Doc:
   ```
   URL: https://docs.google.com/document/d/e/2PACX-1vQnt7OCPm19vnWF4yynXL9JTzTvq9CrGoEaDv7yFSngLoFsypiWsx6fZLKWwaO5YQ/pub
   ```

2. **Save to `docs/SIGNATURE_GUIDELINES.md`** (overwrite previous content; update the `Last synced` date in the header)

3. **Compare with previous version** to identify changes:
   - New signature definitions or requirements
   - Updated QA/QC flags or thresholds
   - New function requirements or parameters
   - Comments or suggestions from colleagues

4. **Add new TODOs to `CHANGELOG.md`** under `[Unreleased]` → `### Guidelines Document TODOs`

5. **Present changes to user**:
   > "Guidelines document has X new/changed items. Would you like to review and implement?"

6. **Implementation workflow**: For each suggestion:
   - Create todo item
   - Implement the change in Julia first (`julia/src/`)
   - Run Julia benchmark (`docs/benchmarks/run_julia_benchmark.jl`, ~27 min) to verify
   - Port to Python and rpkg
   - Mark todo complete

### B. Manuscript Draft (HISSS paper — reconciliation review)

The first manuscript leveraging this workflow (Scientific Data, "HISSS", submission target Nov 9 2026) describes the methods implemented here. Its claims must stay consistent with both the code and the repo documentation.

1. **Fetch fresh content** from the Google Doc:
   ```
   URL: https://docs.google.com/document/d/e/2PACX-1vS7j4FRp7SEwlXoBUVA8NA7cj_I0XzyS0u58r3bl8SOz4BfpZPrdPJge4RMcFocnX8Gnllkc1M-CTJ3/pub
   ```

2. **Compare with the saved snapshot** `docs/MANUSCRIPT_DRAFT.md` (diff the extracted text against the snapshot body below its header). If unchanged, report "manuscript unchanged" and stop here.

3. **If edited**: overwrite the snapshot (update `Last synced`), then run a **reconciliation review** of every changed methods claim in three directions:
   - **Manuscript vs code** — is the described method correctly and precisely implemented (Julia canonical `julia/src/` + `config/signatures_config.json`)?
   - **Manuscript vs repo docs** — is it correctly and precisely documented (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`, `docs/SIGNATURE_GUIDELINES.md`)?
   - **Direction of fix** — if the code/docs are right and the manuscript is wrong, the fix belongs in the Google Doc (flag for the user to relay to co-authors — the doc cannot be edited from here); if the manuscript states the agreed methodology and the code/docs lag, it becomes an implementation TODO.

4. **Log all resulting updates (implemented and planned)** in `CHANGELOG.md` under `[Unreleased]` → `### Manuscript Reconciliation Log`, as dated entries. Move entries to the dated release section once resolved.

5. **Present findings to user** with the discrepancies grouped by direction of fix.

## Changelog Maintenance

**CHANGELOG.md must be kept updated consistently:**

1. **Document all code changes** - Every bug fix, feature, or modification
2. **Track guidelines implementation** - When implementing suggestions from `docs/SIGNATURE_GUIDELINES.md`
3. **Use date-based versioning** - Format: `[Month Year]` (e.g., `[March 2026]`)
4. **Severity labels** - Use HIGH/MEDIUM/LOW for bug fixes
5. **New suggestions** - Add under `[Unreleased]` → `### Guidelines Document TODOs`
6. **Manuscript reconciliation** - Log implemented and planned updates from manuscript syncs under `[Unreleased]` → `### Manuscript Reconciliation Log` (dated entries)

## Claude Skill Maintenance

**Update `claude-skill/streamflow-signatures.md` whenever:**

1. **User feedback** - Recurring questions, confusion points, or feature requests
2. **Novel findings** - New understanding of signatures, edge cases, or best practices
3. **Workflow updates** - Changes to processing pipelines, data formats, or validation
4. **Methodology changes** - Updated formulas, parameters, or statistical approaches
5. **Cross-language updates** - When Python/Julia implementations are added or modified

The skill helps users interpret outputs, understand methodology, and troubleshoot issues.

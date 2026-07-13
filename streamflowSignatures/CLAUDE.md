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
5. **NA Handling** (April 2026): `preprocess_daily_data()` runs ONCE per gage BEFORE signatures.
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
4. Add base name to `EXPECTED_SIGNATURE_BASES` in `config.R`
5. Run Julia benchmark (`docs/benchmarks/run_julia_benchmark.jl`, ~27 min) to verify
6. Port to Python (`python/streamflow_signatures/`) and rpkg (`rpkg/R/`)

## References

- **@docs/DEVELOPMENT.md** - Architecture, file structure, common tasks, workflows
- **@docs/SIGNATURES.md** - Detailed signature documentation (11 categories)
- **@CHANGELOG.md** - Current work, roadmap (historical: `docs/CHANGELOG_ARCHIVE.md`)
- **@docs/SIGNATURE_GUIDELINES.md** - Collaborative guidelines from hydrology colleagues (auto-synced)
- **@EO_data_processing/README.md** - Earth Observation (MODIS LULC & LAI) per-watershed ingestion and processing

## Session-Start Workflow: Guidelines Sync

The guidelines document is a core design feature: domain experts write plain-English signature definitions and QA/QC requirements, and those are translated into code. **At the start of each session**, sync the collaborative guidelines document:

1. **Fetch fresh content** from the Google Doc:
   ```
   URL: https://docs.google.com/document/u/1/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub
   ```

2. **Save to `docs/SIGNATURE_GUIDELINES.md`** (overwrite previous content)

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

## Changelog Maintenance

**CHANGELOG.md must be kept updated consistently:**

1. **Document all code changes** - Every bug fix, feature, or modification
2. **Track guidelines implementation** - When implementing suggestions from `docs/SIGNATURE_GUIDELINES.md`
3. **Use date-based versioning** - Format: `[Month Year]` (e.g., `[March 2026]`)
4. **Severity labels** - Use HIGH/MEDIUM/LOW for bug fixes
5. **New suggestions** - Add under `[Unreleased]` → `### Guidelines Document TODOs`

## Claude Skill Maintenance

**Update `claude-skill/streamflow-signatures.md` whenever:**

1. **User feedback** - Recurring questions, confusion points, or feature requests
2. **Novel findings** - New understanding of signatures, edge cases, or best practices
3. **Workflow updates** - Changes to processing pipelines, data formats, or validation
4. **Methodology changes** - Updated formulas, parameters, or statistical approaches
5. **Cross-language updates** - When Python/Julia implementations are added or modified

The skill helps users interpret outputs, understand methodology, and troubleshoot issues.

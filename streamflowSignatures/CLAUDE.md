# Streamflow Signatures Project

## Project Context

This project has two primary goals and two secondary goals:

1. **Data Processing** — Ingest raw streamflow data (USGS, HYDAT, Caravan), clean/filter/collate metadata, standardize outputs.
2. **Signature Extraction** — Extract 100+ hydrological signatures under strict guardrails. Domain experts update methodology via plain-English instructions in `docs/SIGNATURE_GUIDELINES.md`; code implements those definitions.
3. **Visualization** (secondary) — Shiny dashboard for broader audience exploration.
4. **Cross-Language Implementations** (secondary) — R is canonical; Python and Julia produce near-identical results for community sharing (future: publishable packages/libraries).

## Multi-Language Architecture

This project provides identical signature calculations in R (canonical), Python, and Julia:

| Directory | Language | Status | Description |
|-----------|----------|--------|-------------|
| `R/` | R | **Canonical** | Reference implementation - all changes start here |
| `python/` | Python | Active | Port of R signatures |
| `julia/` | Julia | Active | Port of R signatures |

**Change Workflow**: R is canonical. Changes are made in R first, then propagated to Python/Julia. Golden outputs from R validate other implementations.

## Canonical Code

**Always source in this order:**
1. `config.R` - Configuration, logging, validation (source FIRST)
2. `R/helperFunctions.R` - All core functions (45+ functions)

Other `helperFunctions*.R` files in `archive/` are deprecated.

## Key Entry Points

- `run_full_processing.R` - PRIMARY: Full signature extraction with climate
- `run_ingest_usgs_hydat.R` - Raw USGS/HYDAT data ingestion to parquet
- `run_caravan_processing.R` - Caravan data processing (lower priority)
- `streamflowAndClimateVisualizationApp/app.R` - Shiny dashboard

## Critical Constraints

1. **CSV Output Format**: MUST remain unchanged - downstream tools depend on exact column names
2. **Water Year**: Oct 1 - Sep 30
3. **Flow Units**: mm/day (converted from cfs/m3s)
4. **Minimum Data Requirements**:
   - 20+ water years per gage
   - 95% non-NA days per water year
   - 30+ days above minimum flow threshold

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

## Code Status

| File | Status | Notes |
|------|--------|-------|
| `config.R` | **ACTIVE** | Source first |
| `R/helperFunctions.R` | **CANONICAL** | All core functions |
| `run_full_processing.R` | **PRIMARY** | Main entry point |
| `run_ingest_usgs_hydat.R` | **ACTIVE** | Raw data ingestion |
| `R/tests/smoke_test.R` | **ACTIVE** | Quick validation |
| `R/tests/qa_qc_signatures.R` | **ACTIVE** | Output validation |
| `archive/*` | **DO NOT USE** | Reference only |

## Adding New Signatures

1. Create function in `R/helperFunctions.R` returning annual values
2. Call `generate_stats()` to produce 8 statistics
3. Add base name to `EXPECTED_SIGNATURE_BASES` in `config.R`
4. Run `R/tests/smoke_test.R` to verify

## References

- **@docs/DEVELOPMENT.md** - Architecture, file structure, common tasks, workflows
- **@docs/SIGNATURES.md** - Detailed signature documentation (10 categories)
- **@CHANGELOG.md** - Bug fixes, known issues, roadmap
- **@docs/SIGNATURE_GUIDELINES.md** - Collaborative guidelines from hydrology colleagues (auto-synced)

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
   - Implement the change
   - Run `R/tests/smoke_test.R` to verify
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

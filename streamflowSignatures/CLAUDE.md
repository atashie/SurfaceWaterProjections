# Streamflow Signatures Project

R-based hydrological signature extraction from USGS, HYDAT, and Caravan data. Calculates many metrics with standardized trend analysis for all metrics.

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
| `tests/smoke_test.R` | **ACTIVE** | Quick validation |
| `tests/qa_qc_signatures.R` | **ACTIVE** | Output validation |
| `archive/*` | **DO NOT USE** | Reference only |
| `orphans/*` | **DO NOT USE** | Deprecated files |

## Adding New Signatures

1. Create function in `R/helperFunctions.R` returning annual values
2. Call `generate_stats()` to produce 8 statistics
3. Add base name to `EXPECTED_SIGNATURE_BASES` in `config.R`
4. Run `tests/smoke_test.R` to verify

## References

- **@DEVELOPMENT.md** - Architecture, file structure, common tasks, workflows
- **@SIGNATURES.md** - Detailed signature documentation (10 categories)
- **@CHANGELOG.md** - Bug fixes, known issues, roadmap
- **@SIGNATURE_GUIDELINES.md** - Collaborative guidelines from hydrology colleagues (auto-synced)

## Session-Start Workflow: Guidelines Sync

**At the start of each session**, sync the collaborative guidelines document:

1. **Fetch fresh content** from the Google Doc:
   ```
   URL: https://docs.google.com/document/u/1/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub
   ```

2. **Save to `SIGNATURE_GUIDELINES.md`** (overwrite previous content)

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
   - Run `tests/smoke_test.R` to verify
   - Mark todo complete

## Changelog Maintenance

**CHANGELOG.md must be kept updated consistently:**

1. **Document all code changes** - Every bug fix, feature, or modification
2. **Track guidelines implementation** - When implementing suggestions from `SIGNATURE_GUIDELINES.md`
3. **Use date-based versioning** - Format: `[Month Year]` (e.g., `[March 2026]`)
4. **Severity labels** - Use HIGH/MEDIUM/LOW for bug fixes
5. **New suggestions** - Add under `[Unreleased]` → `### Guidelines Document TODOs`

# Streamflow Signatures Project

R-based hydrological signature extraction from USGS, HYDAT, and Caravan data. Calculates many metrics with standardized trend analysis for all metrics.

## Canonical Code

**Always source in this order:**
1. `config.R` - Configuration, logging, validation (source FIRST)
2. `helperFunctions.R` - All core functions (45+ functions)

Other `helperFunctions*.R` files are deprecated.

## Key Entry Points

- `run_full_processing.R` - PRIMARY: Full signature extraction with climate
- `caravan_to_annualized.R` - Caravan data processing - no longer a priority, as we preference processing "live" data via api from USGS and HYDAT
- `streamflowVisualizationApp/app.R` - Shiny dashboard

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
| `helperFunctions.R` | **CANONICAL** | All core functions |
| `run_full_processing.R` | **PRIMARY** | Main entry point |
| `tests/smoke_test.R` | **ACTIVE** | Quick validation |
| `tests/qa_qc_signatures.R` | **ACTIVE** | Output validation |
| `streamflowDataProcessing*.R` | **LEGACY** | Use sparingly |
| `archive/*` | **DO NOT USE** | Reference only |

## Adding New Signatures

1. Create function in `helperFunctions.R` returning annual values
2. Call `generate_stats()` to produce 8 statistics
3. Add base name to `EXPECTED_SIGNATURE_BASES` in `config.R`
4. Run `tests/smoke_test.R` to verify

## References

- **@DEVELOPMENT.md** - Architecture, file structure, common tasks, workflows
- **@SIGNATURES.md** - Detailed signature documentation (10 categories)
- **@CHANGELOG.md** - Bug fixes, known issues, roadmap
- **@.claude/rules/validation.md** - Logging and validation API reference

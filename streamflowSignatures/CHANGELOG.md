# Changelog

All notable changes to the Streamflow Signatures project.
For historical entries (Dec 2025 – March 2026), see [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

## [Unreleased]

### Planned
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)

### Julia leftjoin Row-Order Bug Fix (April 2026)
**HIGH**: `preprocess_daily_data()` in `julia/src/io.jl` produced incorrect `valid_climate_years` because `DataFrames.jl`'s `leftjoin` does not preserve row order. After joining the daily grid with year data, rows were reordered so that missing-PPT dates (e.g., Dec 31 at dowy=92) appeared at boundary positions (position 365 = Sep 30). The interpolation logic uses positional indices, so 1-day internal gaps were misidentified as boundary gaps and could not be interpolated — causing those years to be excluded from `valid_climate_years`.

**Fix**: Added `sort!(joined, :dowy)` after the leftjoin (1 line). For gage 01011000, `valid_climate_years` went from 32 → 43 (matching Python). Across the full benchmark, 0 entirely-NA climate columns (was 48 before fix).

### Interactive Comparison Dashboard (April 2026)
New `docs/benchmarks/build_comparison_dashboard.py` generates an interactive HTML dashboard (`signature_comparison.html`) comparing Julia benchmark output against golden R reference. Features: two synced Leaflet maps (original filter vs new filter), Plotly scatterplot with identity line, click-to-zoom from scatterplot to maps, R² summary. Covers 159 signature columns across 5,687 common gages.

### Remove Per-Signature Min-Days Thresholds (April 2026)
Removed all per-year min_days/max_na_frac/min_data_completeness checks from signature functions across all 4 codebases. The preprocessor (`preprocess_daily_data()`) is now the single source of truth for year qualification.

**Root cause**: Two independent year-qualification systems were fighting each other — per-signature min_days thresholds rejected preprocessor-approved years, creating artificial NaN years that then failed the 80% trend completeness check, causing massive NA inflation (e.g., elasticity dropped from 5,683 to 19 gages).

**Removed from `config/signatures_config.json`** (13 per-year thresholds):
- `baseflow.min_days`, `baseflow.max_missing_frac`
- `fdc.min_days`, `flashiness.min_days`, `flashiness.max_missing_frac`
- `flow_volumes.min_days`, `timing.min_days`, `pulses.min_days`
- `elasticity.min_data_completeness`
- `qp_seasonality.min_days_per_year`, `qp_seasonality.max_na_frac`
- `storage.min_days_per_year`, `storage.max_na_frac`

**Retained** (gage-level mathematical minimums, not per-year data quality):
- `recession.min_length=5`, `recession.min_events=25`
- `elasticity.min_years=15`, `elasticity.min_annual_ppt=10`
- `qp_seasonality.min_years=10`, `storage.min_years=10`

**Set `use_legacy_filtering: false`** permanently.

**Exempted recession and elasticity from trend_completeness** across all 4 codebases (recession is event-based/inherently sparse, elasticity uses rolling windows producing fewer values than input years).

**Updated benchmark runners** (Python, rpkg) to use `preprocess_daily_data()` path with trend_completeness/decade_completeness, matching Julia benchmark runner.

### Julia Adversarial Review (April 2026, Codex)
6 HIGH, 1 MEDIUM, 1 LOW findings. No critical bugs. Fixes deferred.

- **HIGH** FDC/flashiness min-days threshold 250 vs R's 30 (`fdc.jl`, `flashiness.jl`, `config.jl`) — silently drops years R keeps
- **HIGH** `Q95_Q10` column name vs R's `Q95-Q10` (`flow_volumes.jl`) — schema mismatch
- **HIGH** `generate_stats()` returns empty Dict on short data instead of NaN row (`stats.jl`) — breaks schema
- **HIGH** Runoff ratios aggregate Q/PPT jointly not independently (`runoff_ratios.jl`) — systematic bias when missingness differs
- **HIGH** Mann-Kendall tau denominator not tie-adjusted (`stats.jl`) — likely explains some of the 20 poor benchmark columns
- **MEDIUM** All-NA season yields NaN vs R's 0 (`flow_volumes.jl`)
- **LOW** Arrow IPC output mislabeled as parquet (`io.jl`)

### Standardized NA Handling (April 2026)
Centralized pre-processing for missing data, implemented identically across R canonical, rpkg, Python, and Julia. Driven by Guidelines Document "NA Handling" section.

**New `preprocess_daily_data()` function** (all 4 codebases):
- Daily grid normalization: one row per day per water year, sorted, unique dates
- Linear interpolation of internal gaps <= 3 consecutive days (no extrapolation)
- Year rejection: >30 raw NAs, gaps >3 days, negative Q values
- Residual boundary NAs (leading/trailing) cause year rejection
- Constant-SD QA flag: detects sensor flatline (flag only, not rejection)
- Seasonal completeness flags from RAW observations (pre-interpolation)
- Separate climate (PPT) NA policy with independent tracking

**Config changes** (`config/signatures_config.json`):
- New `na_handling` section with interpolation, year_rejection, constant_sd_flag, trend_completeness, seasonal_completeness, climate_na_policy
- `use_legacy_filtering: true` flag for backward-compatible staged migration

**fillna(0) removal** (3 functions, all 4 codebases):
- `analyze_flow_timing_trends`: no longer replaces NA Q with 0 before cumsum
- `calculate_average_storage`: skips years with residual NAs instead of zero-filling
- `calculate_qp_seasonality`: skips years with residual NAs instead of zero-filling

**Trend completeness** in `generate_stats()`:
- New `trend_completeness` and `decade_completeness` parameters
- Requires >=80% non-NA values overall and in first/last decade for trend stats
- Series <10 years: decade check skipped
- If incomplete: 6 trend stats set to NA, mean/median still computed

**Seasonal completeness** in flow volumes and runoff ratios:
- `calculate_flow_vols_by_year()` and `analyze_Q_PPT_relationships()` accept `seasonal_flags`
- Incomplete seasons (RAW observation fraction <80%) set to NA per year

**Tests**: `R/tests/test_na_handling.R` with 30+ test cases covering all edge cases

### Guidelines Document TODOs
<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

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

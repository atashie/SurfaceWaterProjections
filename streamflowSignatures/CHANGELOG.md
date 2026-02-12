# Changelog

All notable changes to the Streamflow Signatures project.

## [Unreleased]

### Planned
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)

---

## [February 2026]

### Added

#### Signature Cross-Correlation Plot (Means vs Trends)
- **New Shiny dashboard section**: Interactive scatter plot showing how pairwise signature relationships in spatial patterns (means across gages) compare to relationships in temporal trends (Theil-Sen slopes across gages)
- **Architecture**: Pre-computed offline to avoid startup delay; CSV loaded from S3 (~4,556 pairs)
- **Features**:
  - Each point = one ordered pair (A → B); x = Theil-Sen slope of z-scored means, y = Theil-Sen slope of z-scored trends
  - 13-category color coding (Flow Volume, Flow Percentiles, FDC, Baseflow, Recession, Pulse Metrics, Flow Reversals, Flashiness, Flow Timing, Runoff Ratios, Elasticity, Q-P Seasonality, Average Storage)
  - Category checkbox filter to show/hide groups
  - 1:1 reference line showing where spatial and temporal co-variation are equal
  - Click any point to reveal gage-level detail panel with two scatter plots (A_mean vs B_mean, A_senn_slp vs B_senn_slp) including Theil-Sen trend lines and Spearman statistics
- **New files**:
  - `precompute_cross_signature_analysis.R` — Offline script computing all pairwise slopes across 68 metrics × 5,691 gages, uploads CSV to S3
  - `streamflowAndClimateVisualizationApp/cross_signature_analysis.csv` — Pre-computed output (4,556 rows)
- **Modified files**:
  - `streamflowAndClimateVisualizationApp/app.R` — Added category constants, UI section, and server logic
  - `streamflowAndClimateVisualizationApp/helperFunctions.R` — Added `fast_theil_sen_slope()` utility for batch Theil-Sen computation

### Fixed

### Fixed

#### Visualization App Gage ID Mismatch (HIGH)
- **Issue**: Scatter plot incorrectly filtered out 3,647 valid gages as "failed QC"
- **Root Cause**: App compared `gage_id` (with leading zeros, e.g., `01011000`) to metadata `gage_id` (without leading zeros, e.g., `1011000`)
- **Additional Issue**: App redundantly applied `processing_status == "success"` filter already enforced during pre-processing
- **Fix**: Removed `goodGages` filter from scatter plot; only `flagged_for_qann_range` is now applied
- **Location**: `streamflowAndClimateVisualizationApp/app.R`

#### Code Review Bug Fixes (Feb 2026)

**HIGH severity:**
- **H1**: `calculate_flow_vols_by_year()` — Renamed misleading variables from `annual_means`/`*_means` to `annual_totals`/`*_totals`; updated SIGNATURES.md to say "total" with units "mm" instead of "mean" with "mm/day"
- **H2**: `calculate_flow_vols_by_year()` — Added 250-day minimum data filter per water year; previously years with as few as 1 day could pass through
- **H3**: `analyze_recession_parameters()` — Fixed early-return schema: was producing 5 columns with wrong suffixes (`slp`, `rho`, `pval`), now produces correct 8 columns matching `STAT_SUFFIXES`
- **H4**: `analyze_flow_timing_trends()` — `cumsum()` now handles NAs by replacing with 0 before accumulation; previously a single NA corrupted all D-day metrics for that year
- **H5**: Canadian unit conversion — Replaced sentinel value `99999` with `NA` when basin area is missing; previously produced absurd Q values (~100,000+ mm/day)
- **H5 follow-up**: Missing basin area — The H5 fix (`conversion = NA`) caused all Q values to become NA, silently dropping gages at the qualifying-years gate. Now sets `conversion = 1` to keep raw units (cfs for USGS, m3/s for Canadian) and adds `area_normalized` flag to `watershed_metadata.csv`. Also added the same NA guard for USGS gages (previously unprotected).

**MEDIUM severity:**
- **M1**: `lyne_hollick_filter()` — Fixed multi-pass to use previous pass output as input; the forward pass was always referencing original Q regardless of pass number, making `passes=2` effectively single-pass
- **M2**: `analyze_recession_parameters()` — Point cloud `log_a` now uses `b_pointcloud` instead of `median_b` from per-event fits
- **M4**: `analyze_flashiness_trends()` — Added division-by-zero guard when all flow values are zero
- **M5**: `count_flow_reversals()` — NA gaps are now interpolated instead of removed, preventing false adjacencies between non-adjacent days
- **M6**: `calculate_qp_seasonality()` and `calculate_average_storage()` — Use `copy()` to prevent modifying parent data.table by reference via `:=`
- **M7**: `calculate_qp_seasonality()` — Slopes now assigned to middle of rolling window instead of end
- **M8**: `analyze_Q_PPT_relationships()` — Raised near-zero PPT threshold from 0.001mm to 10mm (annual) and 1mm (seasonal) to prevent extreme runoff ratios

**LOW severity:**
- **L1**: `generate_stats()` — Now explicitly sorts by year before trend calculations
- **L3**: `analyze_baseflow_indices()` — Removed unused `doy` from required columns (only `dowy` is used)
- **L5**: `process_signatures_from_parquet()` — Non-climate error handlers now log warnings instead of silently swallowing errors
- **L7**: `analyze_flashiness_trends()` — Minimum-day check now counts non-NA values instead of total rows
- **L8**: `process_signatures_from_parquet()` — Renamed loop variable from `gage_id` to `current_gage_id` to avoid data.table column name shadowing

### Added
- Comprehensive workflow review document (`WORKFLOW_REVIEW.md`)
- Documentation of gage ID format differences (`gage_id` vs `gage_id_metadata` columns)
- `area_normalized` column in `watershed_metadata.csv` — `TRUE` when basin area is available and Q is in mm/day, `FALSE` when basin area is missing and Q is in raw units

### Changed
- SIGNATURES.md updated to match actual code: column names, metric counts, units, PPT thresholds, D-day naming (D5 not D05), recession seasonality variants

---

## [January 2026]

### Fixed

#### Metadata Lookup Bug (HIGH)
- **Issue**: Data.table scoping caused all gages to receive the first row's metadata
- **Root Cause**: `metadata_lookup[gage_id == gage_id]` compared column to itself (always TRUE)
- **Fix**: Renamed `find_metadata()` parameter to `target_gage_id`
- **Location**: `helperFunctions.R` line ~3425

#### Canadian Basin Area Missing (HIGH)
- **Issue**: Basin area was hardcoded as `NA` for Canadian stations
- **Impact**: Canadian gages had no drainage area in output
- **Fix**: Now fetches `DRAINAGE_AREA_GROSS` from `tidyhydat::hy_stations()`
- **Locations**:
  - Metadata creation functions (for new processing)
  - `process_signatures_from_parquet()` runtime fallback (for existing metadata)

### Added
- Centralized `config.R` for all configuration parameters
- Structured logging system with levels (DEBUG, INFO, WARN, ERROR)
- Input validation functions (`validate_file_exists`, `validate_numeric`, etc.)
- Output schema validation (`validate_output_schema`, `validate_gage_output`)
- QA/QC scripts (`qa_qc_signatures.R`, `visualize_qa_qc.R`)
- Smoke test for quick validation (`smoke_test.R`)
- Climate function tests with synthetic data (`test_climate_functions.R`)

### Changed
- Consolidated all helper function variants into canonical `helperFunctions.R`
- Moved deprecated helper files to `archive/` directory

### Completed Milestones
- Consolidate active helper file variants into canonical `helperFunctions.R`
- Add centralized `config.R` for parameters
- Implement structured logging
- Fix metadata lookup bug
- Fix Canadian basin_area bug

---

## [December 2025]

### Added
- Initial climate signatures: elasticity, Q-P seasonality, average storage
- Daymet climate data integration
- Caravan processing pipeline (`caravan_to_annualized.R`)

---

## Version History Notes

This project uses date-based versioning (MONTH YEAR) rather than semantic versioning, reflecting its nature as a research tool with continuous development.

### Output File Naming Convention
Output files include date stamps: `streamflow_signatures_full_JAN2026.csv`

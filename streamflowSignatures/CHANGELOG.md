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

### Fixed

#### Visualization App Gage ID Mismatch (HIGH)
- **Issue**: Scatter plot incorrectly filtered out 3,647 valid gages as "failed QC"
- **Root Cause**: App compared `gage_id` (with leading zeros, e.g., `01011000`) to metadata `gage_id` (without leading zeros, e.g., `1011000`)
- **Additional Issue**: App redundantly applied `processing_status == "success"` filter already enforced during pre-processing
- **Fix**: Removed `goodGages` filter from scatter plot; only `flagged_for_qann_range` is now applied
- **Location**: `streamflowAndClimateVisualizationApp/app.R`

### Added
- Comprehensive workflow review document (`WORKFLOW_REVIEW.md`)
- Documentation of gage ID format differences (`gage_id` vs `gage_id_metadata` columns)

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

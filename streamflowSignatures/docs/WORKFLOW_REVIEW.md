# Streamflow Signatures Workflow Review - Comprehensive Analysis

**Date:** February 2026
**Focus:** Data quality, data quality flags, and potential anomalies in outputs

---

## Executive Summary

This document provides a comprehensive review of the entire streamflow signatures processing workflow, from raw data extraction through visualization. Seven specialized analyses were conducted covering:

1. Streamflow data extraction (USGS/HYDAT)
2. Weather data processing (Daymet CSV to Parquet)
3. Data concatenation/joining
4. Metadata generation
5. Signature calculations
6. Visualization app
7. Data quality flagging

**Key Findings:**
- Two HIGH-priority bugs were fixed in January 2026 (metadata lookup, Canadian basin area)
- 11 data quality flags are generated per gage
- Multiple data quality thresholds enforce minimum standards (20 years, 95% completeness)
- Several gaps identified in anomaly detection and cross-gage validation

---

## 1. Streamflow Data Extraction (USGS & HYDAT)

### Data Sources & APIs

| Source | API | Parameter | Raw Units |
|--------|-----|-----------|-----------|
| USGS | `dataRetrieval::readNWISdv()` | 00060 (Discharge) | cfs |
| HYDAT | `tidyhydat::hy_daily()` | Flow | m³/s |

### Quality Code Filtering

**USGS Accepted Codes:** `"A", "A e", "P", "P e"`
- A = Approved, P = Provisional, e = Estimated
- All other codes → data masked as `NA`

**Location:** `helperFunctions.R` lines 657-658

### Unit Conversions

| Source | Formula | Result |
|--------|---------|--------|
| USGS | `cfs × 86400 / (km² × 3280.84³) × 10⁶` | mm/day |
| HYDAT | `m³/s × 86400 × 10⁹ / (km² × 10¹²)` | mm/day |

**Location:** `helperFunctions.R` lines 663, 699

### Data Quality Thresholds (config.R)

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `MIN_Q_VALUE_AND_DAYS` | `c(0.0001, 30)` | Min 30 days > 0.0001 mm/day per water year |
| `MIN_NUM_YEARS` | 20 | Minimum qualifying water years |
| `MIN_FRAC_GOOD_DATA` | 0.95 | 95% non-NA days per water year |

### Known Issues (FIXED)

1. **Metadata Lookup Bug** - All gages received first row's metadata due to data.table scoping collision (`gage_id == gage_id`). Fixed by renaming parameter to `target_gage_id`.

2. **Canadian Basin Area** - Was hardcoded as NA. Now dynamically fetched from `tidyhydat::hy_stations()`.

### Potential Anomalies

- **Zero drainage area:** Causes conversion factor = 99999 (implicit data masking)
- **Leading zeros inconsistency:** Multi-format lookup implemented (original, stripped, padded)
- **Regulated Canadian stations:** Filtered out via `REGULATED != TRUE`

---

## 2. Weather Data Processing (Daymet)

### Conversion Pipeline

**Input:** `daymet_1980_2023.zip` (CSV files per year)
**Output:** `daymet_1980_2023.parquet` (Snappy compression)
**Function:** `convert_daymet_zip_to_parquet()` at `helperFunctions.R` lines 2686-2820

### Variables Stored

| Variable | Full Name | Units |
|----------|-----------|-------|
| prcp | Precipitation | mm/day |
| tmin | Min Temperature | °C |
| tmax | Max Temperature | °C |
| swe | Snow Water Equivalent | mm |
| vp | Vapor Pressure | Pa |
| srad | Solar Radiation | W/m² |

### Data Quality Checks

1. **Day count validation:** Warns if sites have ≠365/366 days per year
2. **Duplicate detection:** Warns if duplicate (site_id, Date) combinations found
3. **Coverage reporting:** Warns if <95% of streamflow dates have matching climate data

### Potential Anomalies

- **Day reconstruction assumption:** Days reconstructed sequentially within groups - could be incorrect if source data has gaps
- **NA replacement with zero:** Small gaps in Q or PPT replaced with 0 (not interpolation)
- **Temperature variables unused:** tmin, tmax, vp, srad available but not used in signatures

---

## 3. Data Concatenation (Streamflow + Climate)

### Join Mechanism

**Function:** `integrate_daymet_with_streamflow()` at lines 2880-2931
**Join Type:** LEFT OUTER JOIN on `Date` column
**Key:** Single-column join on exact date match

### Join Sequence

```
1. Filter Daymet parquet by gage_id (site_id match)
2. Filter by date range from streamflow data
3. merge(streamflow, daymet, by="Date", all.x=TRUE)
4. Rename prcp → PPT for downstream compatibility
```

### Missing Data Handling

- **All streamflow records preserved** (left outer join)
- **Missing climate → NA** (no imputation)
- **Coverage tracking:** Warns if <95% dates matched
- **Graceful degradation:** If Daymet unavailable, non-climate signatures still calculated

### Caravan Difference

Caravan NetCDF files bundle Q + climate data together - **no join required**. Variables aligned by array index position.

---

## 4. Metadata Generation

### Fields Tracked

| Field | Source (USGS) | Source (HYDAT) | Source (Caravan) |
|-------|---------------|----------------|------------------|
| gage_id | STAID | STATION_NUMBER | Filename parsing |
| latitude | LAT_GAGE | LATITUDE | NA |
| longitude | LNG_GAGE | LONGITUDE | NA |
| basin_area | DRAIN_SQKM | tidyhydat API | NA |
| gage_type | "USGS" | "Canada" | "Caravan" |

### Metadata Lookup Strategy

**Function:** `find_metadata()` at lines 4164-4193

1. Exact match on gage_id
2. Padded match (add up to 4 leading zeros)
3. Stripped match (remove leading zeros)
4. Return NULL if no match

### Data Quality Indicators

- **Unmatched gages:** Tracked and saved to `*_unmatched_gages.txt`
- **Missing basin_area:** Canadian gages attempt runtime fetch from tidyhydat
- **Caravan limitation:** All spatial metadata is NA (no workaround)

### Important: Gage ID Formats

The signatures output file contains TWO gage_id columns:

| Column | Format | Example | Use For |
|--------|--------|---------|---------|
| `gage_id` | With leading zeros | `01011000` | Display, parquet lookups |
| `gage_id_metadata` | Without leading zeros | `1011000` | Metadata file matching |

**Critical:** The metadata file (`combined_watershed_metadata.csv`) stores gage IDs WITHOUT leading zeros. When joining signatures to metadata, use `gage_id_metadata`, not `gage_id`.

### Processing Status Field

The `processing_status` field in metadata indicates pre-processing outcome:

| Status | Meaning | Count (Feb 2026) |
|--------|---------|------------------|
| `success` | Passed all pre-processing QC | 6,046 |
| `no_data` | No streamflow data retrieved from API | 10,733 |
| `insufficient_years` | < 20 qualifying water years | 81 |
| `processing` | Incomplete (interrupted) | 134 |

**Important:** The `processing_status` filter is applied ONLY during pre-processing to determine which gages to fetch and process. It should NOT be used to filter gages in downstream applications (visualization app, analysis scripts). All gages present in the signatures output file have already passed pre-processing QC.

---

## 5. Signature Calculations

### Orchestrator Function

`process_signatures_from_parquet()` at lines 4050-4472

### The 8-Statistic Rule

Every signature metric produces 8 statistics via `generate_stats()`:

| Suffix | Statistic | Method |
|--------|-----------|--------|
| `_senn_slp` | Theil-Sen slope | `zyp::zyp.sen` |
| `_linear_slp` | Linear regression slope | `lm()` |
| `_spearman_rho` | Spearman correlation | `cor.test` |
| `_spearman_pval` | Spearman p-value | `cor.test` |
| `_mk_rho` | Mann-Kendall tau | `Kendall::MannKendall` |
| `_mk_pval` | Mann-Kendall p-value | `Kendall::MannKendall` |
| `_mean` | Arithmetic mean | `mean()` |
| `_median` | Median | `median()` |

### Signature Categories

| Category | Function | Requires Climate | Output Columns |
|----------|----------|------------------|----------------|
| Flow Volumes | `calculate_flow_vols_by_year()` | No | ~200 |
| Baseflow | `analyze_baseflow_indices()` | No | 16 |
| Recession | `analyze_recession_parameters()` | No | 46 |
| Pulse Metrics | `calculate_pulse_metrics()` | No | 120 |
| Flashiness | `analyze_flashiness_trends()` | No | 8 |
| Flow Timing | `analyze_flow_timing_trends()` | No | 104 |
| FDC | `analyze_fdc_trends_from_streamflow()` | No | 24 |
| Q-PPT Ratios | `analyze_Q_PPT_relationships()` | Yes | 40 |
| Elasticity | `calculate_streamflow_elasticity()` | Yes | 9 |
| Q-P Seasonality | `calculate_qp_seasonality()` | Yes | 16 |
| Average Storage | `calculate_average_storage()` | Yes | 8 |

### Exceptions to 8-Stat Rule

| Column | Type | Reason |
|--------|------|--------|
| `elasticity_static` | Single value | Overall catchment elasticity |
| `log_a_seasonality_amplitude_*` | Single values (3) | Recession seasonality |
| `log_a_seasonality_minimum_*` | Single values (3) | Recession seasonality |

### Per-Signature Data Requirements

| Signature | Min Days/Year | Other Requirements |
|-----------|---------------|-------------------|
| Baseflow | 250 | ≤20% missing |
| Flashiness | 30 | Sorted by dowy |
| Flow Timing | 300 | Total annual flow > 0 |
| FDC | 30 | - |
| Elasticity | - | 90% complete, P > 10mm/yr, ≥15 years |
| Q-P Seasonality | 300 | ≤10% NA in Q or PPT |
| Storage | 300 | ≤10% NA, ≥10 S-Q pairs |

### Known Limitations

1. **Average Storage ignores ET:** Uses P - Q only, may overestimate storage
2. **Zero precipitation:** Handled with thresholds (0.001-0.1 mm depending on metric)
3. **Sort order dependencies:** Several signatures require dowy sort

---

## 6. Data Quality Flagging System

### Flag Definitions

| Flag | Metric | Valid Range | Threshold |
|------|--------|-------------|-----------|
| `flagged_for_qann_range` | Qann_mean | [0, 2000] mm/yr | 95% in range |
| `flagged_for_bfi_eckhardt_range` | BFI_Eckhardt_mean | [0, 1] | 95% in range |
| `flagged_for_bfi_lynehollick_range` | BFI_LyneHollick_mean | [0, 1] | 95% in range |
| `flagged_for_flashiness_range` | flashinessRB_mean | [0, 2] | 95% in range |
| `flagged_for_tqmean_range` | TQmean_mean | [0, 100]% | 95% in range |
| `flagged_for_d50_range` | D50_day_mean | [1, 366] | 95% in range |
| `flagged_for_elasticity_range` | elasticity_static | [0.1, 5] | 95% in range |
| `flagged_for_runoff_ratio_range` | RR_ann_mean | [0.01, 1.5] | 95% in range |
| `flagged_for_seasonal_sum` | Seasonal Q sum | [0.8, 1.2] × Qann | Ratio check |
| `flagged_for_percentile_order` | Q percentiles | Q05 < Q25 < Q50 < Q75 < Q95 | 99% ordered |
| `flagged_for_high_na` | All signatures | >30% NA | Per-gage |

### Consistency Checks

| Check | Expected Relationship | Pass Threshold |
|-------|----------------------|----------------|
| BFI methods | BFI_Eckhardt < BFI_LyneHollick | Correlation ≥ 0.7 |
| Flow timing | D05 < D50 < D95 | 95% ordered |
| BFI vs Flashiness | Negative correlation | r < 0 |
| Slopes agreement | Theil-Sen same sign as Linear | 80% concordance |
| Rank tests | Spearman ~ Mann-Kendall | r ≥ 0.9 |

### Statistical Sanity Tests

- Overall NA fraction: <20% = PASS, <30% = WARNING
- Elasticity distribution: median [0.5, 3], IQR < 3
- Storage validity: 80% positive, 90% < 1000mm

**Location:** `R/tests/qa_qc_signatures.R` lines 77-430

---

## 7. Visualization App

### Data Sources (S3)

- Metadata: `streamflow/combined_watershed_metadata.csv`
- Signatures: `streamflow/streamflow_signatures_full_JAN2026.csv`
- Boundaries: `unified_watershedBoundaries_simplified.gpkg`
- Streamflow: `streamflow/combined_streamflow_data.parquet`
- Daymet: `streamflow/daymet_1980_2023.parquet`

### Filtering Philosophy

**Important:** The visualization app does NOT filter by `processing_status` or `goodGages`. All gages present in the signatures file are displayed because:

1. Pre-processing QC already filtered gages before signature calculation
2. Any gage in the signatures file has passed: metadata matching, 20+ water years, 95% data completeness
3. Additional filtering in the app would incorrectly exclude valid gages

**Scatter Plot Filtering:** Only `flagged_for_qann_range` is applied to remove gages with Qann values outside [0, 2000] mm/year.

### Quality Display Features

1. **Gap visualization:** Red rectangles for data gaps >1 day
2. **NA markers:** Orange "X" at dates with missing values
3. **Flag filtering:** Scatter plots remove only `flagged_for_qann_range` gages
4. **Significance highlighting:** P-value < 0.05 gets thick black outline

### Weather Overlays

6 Daymet variables available as optional overlays:
- Precipitation (inverted bars)
- Min/Max Temperature
- Snow Water Equivalent
- Vapor Pressure
- Solar Radiation

---

## 8. Identified Gaps & Recommendations

### Issues Fixed (February 2026)

| Issue | Description | Resolution |
|-------|-------------|------------|
| **Gage ID mismatch in app** | App compared `signature_data$gage_id` (with leading zeros) to `metadata$gage_id` (without), causing 3,647 gages to appear as "failed QC" | Removed `goodGages` filter from scatter plot; all gages in signatures file already passed pre-processing QC |
| **Redundant QC filtering** | App re-applied `processing_status == "success"` filter that was already enforced during pre-processing | Filter removed; signatures file is the authoritative source of valid gages |

### Current Gaps

| Gap | Description | Priority |
|-----|-------------|----------|
| No seasonal anomaly detection | Flags are annual only | Low |
| Recession seasonality not validated | Custom columns could contain NaN/Inf | Medium |
| No data traceability | Flags don't link to specific years/months | Low |
| ET not included in storage | Simplified P-Q balance | High (documented) |
| Caravan metadata missing | All spatial data is NA | Low |

### Recommended Improvements

1. **Add month-of-year statistics** to detect systematic seasonal biases
2. **Validate recession seasonality columns** for NaN/Inf values
3. **Enhanced climate signature validation:**
   - Check elasticity sign (should be positive)
   - Validate runoff ratio vs precipitation trends
   - Flag when Q-P window has <90% completeness
4. **Add data lineage tracking** to link flags to source data issues
5. **Implement regime classification** (natural vs regulated, perennial vs ephemeral)

---

## 9. Key File Reference

| Component | File | Key Lines |
|-----------|------|-----------|
| Configuration | `config.R` | 1-612 |
| Core functions | `helperFunctions.R` | 1-4472 |
| Main processing | `run_full_processing.R` | 1-96 |
| QA/QC flagging | `R/tests/qa_qc_signatures.R` | 77-430 |
| QA/QC visualization | `R/tests/visualize_qa_qc.R` | 70-650 |
| Smoke test | `R/tests/smoke_test.R` | 1-130 |
| Visualization app | `streamflowAndClimateVisualizationApp/app.R` | 1-1066 |

---

## 10. Summary Statistics

- **Total output columns:** ~595-605 per gage (550+ signatures + metadata + flags)
- **Flag columns:** 11 per gage
- **Signature categories:** 11
- **Climate-dependent signatures:** 4 categories (Q-PPT, Elasticity, Q-P Seasonality, Storage)
- **Data quality thresholds:** 20 years, 95% completeness, 30 days above minimum flow
- **QA/QC tests:** ~10 statistical tests, 16 diagnostic plot types

---

*This review was generated by exploring the codebase with 7 specialized analysis agents, focusing on data quality, flagging mechanisms, and potential anomalies in the processing workflow.*

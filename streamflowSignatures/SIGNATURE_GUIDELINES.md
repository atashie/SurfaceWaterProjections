# Summary Documentation for Streamflow Signatures (2026-02)

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/u/1/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub)
> **Last synced**: 2026-03-05

## Overview
This document provides comprehensive technical guidance for calculating and analyzing streamflow signatures—quantitative metrics describing watershed hydrology characteristics. The material covers calculation methodologies, statistical approaches, data quality protocols, and utility functions for hydrologic analysis.

## 1. Streamflow Signature Categories

### Flow Volume Metrics
Annual and seasonal mean flows (Qann, Qwin, Qspr, Qsum, Qfal) and percentile-based discharge values (Q10, Q50, Q90).

### Baseflow Indices
"BFI_Eckhardt" and "BFI_LyneHollick" quantify groundwater contributions using digital filtering with parameter constraints (maximum 0.8, alpha 0.98).

### Recession Analysis
Characterizes streamflow decline between precipitation events. Median log(a) and exponent b values describe recession behavior, with "concavity" measuring temporal variations across event phases.

### Pulse Metrics
High and low-flow pulse frequency, duration, and flow reversal counts using percentile thresholds (90th and 10th).

### Flashiness & Timing
Richards-Baker index measures flow variability. Cumulative flow timing (D10, D50, D90) indicates seasonal discharge patterns.

## 2. Primary Processing Functions

**process_gages_rawData()**: "Processes raw gage data (USGS/Canadian) to calculate streamflow metrics for individual gages and save results to a file."

**process_caravan_gages()**: Handles NetCDF timeseries with coincident climate data, managing dataset redundancy and truncated records.

**generate_streamflow_dt()**: Standardizes discharge data into consistent tabular format with derived temporal columns.

**calculate_flow_vols_by_year()**: "Seasonal periods are defined as: Winter: December-February, Spring: March-May, Summer: June-August, Fall: September-November."

## 3. Advanced Analysis Functions

**analyze_baseflow_indices()**: Requires minimum 250 valid days annually; excludes years exceeding 20% missing data.

**analyze_recession_parameters()**: Identifies events with ≥5 consecutive monotonic decreases; excludes watersheds with <25 discrete observations.

**calculate_streamflow_elasticity()**: "Measures how sensitively streamflow responds to precipitation changes. Values ~1.0 indicate proportional response; >1 indicates amplified response."

**calculate_qp_seasonality()**: Quantifies cumulative streamflow-precipitation relationships using 30-day rolling slopes, detecting bimodal patterns.

**calculate_average_storage()**: "Daily water balance: dS = P - Q" (simplified without evapotranspiration accounting).

## 4. Statistical Metrics

| Metric | Method | Robustness |
|--------|--------|-----------|
| slp | Theil-Sen | Non-parametric trend estimation |
| rho | Spearman rank | Monotonic correlation |
| mean/median | Central tendency | Arithmetic/robust measures |

## 5. Data Quality Protocols

Automated flags identify problematic records:
- Annual discharge outside 0-2000 mm range
- Baseflow indices beyond 0-1 bounds
- Missing data >30% of signatures
- Percentile ordering violations (Q5<Q25<Q50<Q75<Q95)

Minimum data thresholds: 250 days for baseflow/pulses, 300 days for timing metrics, 15 years for elasticity calculations.

## 6. Utility Functions

**integrate_daymet_with_streamflow()**: Joins precipitation/temperature data to discharge records from pre-processed parquet archive.

**convert_daymet_zip_to_parquet()**: One-time conversion of 44 annual CSV files into optimized storage (~70% compression).

## 7. Advanced Metrics

**Runoff Ratios**: Annual and seasonal Q/P relationships.

**Q-P Seasonality**: "Standard deviation of monthly cumulative Q-P slopes. Higher values indicate stronger seasonal variation."

**Storage Estimation**: Water balance-derived catchment storage in millimeters.

## 8. Key References

- Peters & Aulenbach (2011): Panola Mountain water storage methodology
- Sawicz et al. (2011): Catchment classification via hydrologic similarity
- Wrede et al. (2015): Perceptual model development for Q-P dynamics

## 9. Implementation Notes

Data filtering requires handling frozen days (USGS flags), missing values, and non-zero flow thresholds. Minimum valid years typically range 10-15 depending on metric complexity. Seasonal aggregations demand ≥68% data availability per period.

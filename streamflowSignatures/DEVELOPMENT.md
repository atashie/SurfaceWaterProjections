# Development Guide

This guide covers development workflows, architecture, and common tasks for the Streamflow Signatures project.

## Architecture

### Data Flow

```
Data Sources                    Processing                     Output
────────────                    ──────────                     ──────
USGS (dataRetrieval)  ──┐
  (streamflow only)     ├──> Parquet Storage ──┐
                        │                       │
Canadian HYDAT  ────────┘                       ├──> Signature ──> CSV Summary
  (streamflow only)                             │    Extraction     (550+ columns)
                                                │
Daymet (climate data)  ──> Parquet Storage ─────┘
  (PPT, temp, SWE)         (joined at runtime)

Alternative Pipeline:
Caravan NetCDF ──────────────> Direct Processing ──> Caravan Output
  (bundled Q + climate)        (annualized CSVs)
```

## File Structure

```
streamflowSignatures/
├── CLAUDE.md                    # Claude Code instructions (rules only)
├── README.md                    # User documentation
├── DEVELOPMENT.md               # This file
├── SIGNATURES.md                # Detailed signature documentation
├── CHANGELOG.md                 # Bug fixes, version history
│
│ ## Core Configuration & Functions
├── config.R                     # Centralized configuration parameters
├── helperFunctions.R            # CANONICAL - All core functions (45+ functions)
│
│ ## Main Entry Points
├── run_full_processing.R        # PRIMARY - Full signature extraction with climate data
├── execute_extractStreamflowSignatureValuesAndTrends.R  # Alternative entry point
├── caravan_to_annualized.R      # Caravan data processing
│
│ ## Data Processing Scripts (Legacy)
├── streamflowDataProcessing.R                          # Raw USGS/HYDAT processing
├── streamflowDataProcessing_Caravan.R                  # Caravan processing variant
├── streamflowDataProcessing_USGS-and-Hydat.R          # Reference gage processing
├── streamflowDataProcessing_USGS-and-HYDAT_fullTimeseries.R  # Full timeseries variant
├── streamflowSignatures_wrapperForPreprocessedParquet.R      # Wrapper for pre-processed data
│
│ ## Testing & QA/QC
├── smoke_test.R                 # Quick validation on subset (10 gages)
├── test_climate_functions.R     # Climate function tests with synthetic data
├── qa_qc_signatures.R           # Output validation and QA/QC checks
├── visualize_qa_qc.R            # QA/QC visualization plots
├── test.R                       # Development/exploratory tests (legacy)
├── tests/                       # Unit test directory
│   └── test_climate_signatures.R
│
│ ## Utilities
├── run_conversion.R             # Daymet ZIP to Parquet conversion
│
│ ## Shiny Visualization App
├── streamflowVisualizationApp/
│   ├── app.R                    # Main Shiny application
│   └── helperFunctions.R        # App-specific utilities (S3, legends)
│
│ ## Data & Metadata
├── metadata/                    # Basin and gage metadata (42 files)
├── data_out/                    # Processed outputs (gitignored)
├── test_output/                 # Test outputs (gitignored)
│
└── archive/                     # Archived files (DO NOT USE)
    ├── helperFunction.R
    ├── helperFunctions_sept2025.R
    ├── helperFunctions_extractStreamflowSignatureValuesAndTrends.R
    ├── helperFunctions_processRawStreamflowToParquet.R
    ├── helperWrapperFunctions.R
    └── test_code.txt
```

## Common Tasks

### Run Signature Extraction

```r
source("config.R")              # Load configuration
source("helperFunctions.R")     # Load all functions

summary_output <- process_signatures_from_parquet(
  parquet_file_path = "path/to/streamflow.parquet",
  metadata_file_path = "path/to/metadata.csv",
  output_file = "path/to/output.csv",
  min_Q_value_and_days = MIN_Q_VALUE_AND_DAYS,  # From config.R
  min_num_years = MIN_NUM_YEARS,                 # From config.R
  min_frac_good_data = MIN_FRAC_GOOD_DATA        # From config.R
)
```

### Run Full Processing Pipeline (Recommended)

The easiest way to run a complete signature extraction with climate data:

```bash
# From the streamflowSignatures directory:
Rscript run_full_processing.R
```

This script:
- Loads configuration from `config.R`
- Reads parquet data from `PARQUET_DATA_DIR` (configured in config.R)
- Integrates Daymet climate data if available
- Outputs to `data_out/streamflow_signatures_full_JAN2026.csv`
- Logs progress to `data_out/processing_log_JAN2026.txt`

**Prerequisites:** Edit `config.R` to set `PARQUET_DATA_DIR` to your data location.

### Validate Output Quality

After processing, run QA/QC validation:

```r
source("config.R")
source("helperFunctions.R")
source("qa_qc_signatures.R")
```

Or run the visualization script for diagnostic plots:

```r
source("visualize_qa_qc.R")
# Outputs to data_out/qa_plots/
```

QA/QC checks include:
- Range validation (e.g., BFI in [0,1])
- Baseflow consistency (BFI_Eckhardt < BFI_LyneHollick)
- Elasticity constraints
- Correlation checks between related metrics

### Process Caravan Data

```r
source("helperFunctions.R")

process_caravan_to_annual(
  caravan_directory = "path/to/caravan/netcdf",
  data_project = "camels",
  min_num_years_data = 30,
  start_date_filter = as.Date("1979-09-01"),
  end_date_filter = as.Date("2025-06-01"),
  output_dir = "annualized_caravan_data"
)
```

## Adding a New Signature

### Step-by-Step Process

1. **Create calculation function** in `helperFunctions.R`
   - Function should accept daily data and return annual values
   - Use `data.table` for the return value with a `water_year` column

2. **Apply the 8-statistic rule** using `generate_stats()`:
   ```r
   # Your function should return annual values, then call:
   stats <- generate_stats(annual_data, value_cols = "metric_name", year_col = "water_year")
   # This produces 8 columns: metric_senn_slp, metric_linear_slp,
   # metric_spearman_rho, metric_spearman_pval, metric_mk_rho,
   # metric_mk_pval, metric_mean, metric_median
   ```

3. **Add call to `process_signatures_from_parquet()`** in the signature extraction section

4. **Register the signature** in `config.R`:
   - Add base name to `EXPECTED_SIGNATURE_BASES`

5. **Test with smoke test**:
   ```r
   source("smoke_test.R")
   # Verify schema validation passes
   ```

### Output Column Naming Convention

All signatures follow the pattern: `{metric}_{stat}`

Examples:
- `Qann_senn_slp` - Theil-Sen slope for annual mean flow
- `BFI_Eckhardt_mk_pval` - Mann-Kendall p-value for Eckhardt baseflow index

## Testing

### Quick Validation (Smoke Test)

```r
source("smoke_test.R")
# Runs on 10 gages, validates output schema
```

### Climate Function Tests

```r
source("test_climate_functions.R")
# Tests climate signatures with synthetic data
```

### Unit Tests

```r
# Run all tests in tests/ directory
testthat::test_dir("tests/")
```

## Data Sources

### USGS (via dataRetrieval)
- Parameter code: `00060` (Discharge, cubic feet per second)
- Quality codes accepted: `A`, `A e`, `P`, `P e`
- Conversion: cfs -> mm/day using drainage area

### Canadian HYDAT (via tidyhydat)
- Parameter: Flow (m3/s)
- Excludes regulated stations
- Conversion: m3/s -> mm/day using drainage area

### Caravan
- NetCDF format with daily streamflow + climate variables
- Includes: PPT, SWE, temperature
- Trade-off: Shorter records (ends ~2018-2020) but has climate data

## Dependencies

```r
# Core data handling
library(data.table)
library(arrow)        # Parquet I/O
library(lubridate)    # Date handling

# Data retrieval
library(dataRetrieval) # USGS NWIS API
library(tidyhydat)     # Canadian HYDAT database

# Statistics
library(zyp)          # Theil-Sen slope estimation
library(Kendall)      # Mann-Kendall trend test
library(mblm)         # Alternative Theil-Sen

# Spatial (for basin delineation)
library(sf)
library(terra)

# Visualization app
library(shiny)
library(leaflet)
library(plotly)
library(aws.s3)       # S3 data storage
```

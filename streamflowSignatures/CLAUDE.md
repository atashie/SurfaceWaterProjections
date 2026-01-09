# Streamflow Signatures Project

## Overview

R-based system for extracting hydrological signatures from USGS, Canadian (HYDAT), and Caravan datasets. Calculates 100+ metrics with trend analysis for studying streamflow characteristics and their relationship to climate.

**Primary Research Goal**: Analyze synchrony between streamflow and climate patterns.
**Secondary Goal**: Create comprehensive dataset for broader research community.

## Architecture

### Data Flow
```
Raw Data Sources                Processing                    Output
─────────────────              ────────────                  ──────
USGS (dataRetrieval)  ─┐
                       ├──> Parquet Storage ──> Signature ──> CSV Summary
Canadian HYDAT        ─┤         │              Extraction     (100+ metrics)
                       │         │                   │
Caravan NetCDF ───────┴─────────┴───────────────────┘
  (includes climate)
```

### Canonical Code
- **`helperFunctions.R`** - Primary source file containing all core functions
- Other `helperFunctions*.R` files are deprecated variants (see Code Status section)

### Key Entry Points
- `execute_extractStreamflowSignatureValuesAndTrends.R` - Main signature extraction workflow
- `caravan_to_annualized.R` - Caravan data processing (includes climate variables)
- `streamflowVisualizationApp/app.R` - Shiny dashboard with AWS S3 integration

## Critical Constraints

1. **CSV Output Format**: MUST remain unchanged - downstream tools depend on exact column names and structure
2. **Water Year Convention**: Oct 1 - Sep 30 (Northern Hemisphere standard)
3. **Flow Units**: mm/day after conversion from raw API units (cfs for USGS, m³/s for Canadian)
4. **Minimum Data Requirements**:
   - 20+ water years of valid data per gage
   - 95% non-NA days per water year
   - 30+ days above minimum flow threshold per water year

## Signature Categories

### 1. Flow Volumes (`calculate_flow_vols_by_year`)
- **Qann**: Annual mean streamflow
- **Qwin, Qspr, Qsum, Qfal**: Seasonal means (Dec-Feb, Mar-May, Jun-Aug, Sep-Nov)
- **Q05-Q95**: Flow percentiles (Q05, Q10, Q20, Q25, Q30, Q40, Q50, Q60, Q70, Q75, Q80, Q90, Q95)

### 2. Baseflow (`analyze_baseflow_indices`)
- **BFI_Eckhardt**: Baseflow index using Eckhardt recursive digital filter (BFImax=0.8, a=0.98)
- **BFI_LyneHollick**: Baseflow index using Lyne-Hollick filter (alpha=0.925, 2 passes)

### 3. Recession (`analyze_recession_parameters`)
- **log_a_pointcloud, log_a_events**: Recession rate parameter from dQ/dt = a*Q^b
- **b_pointcloud, b_events**: Recession exponent
- **concavity**: Difference in b between first and second halves of recession
- **log_a_seasonality_amplitude, log_a_seasonality_minimum**: Seasonal variation in recession

### 4. Pulse Metrics (`calculate_pulse_metrics`)
- **n_high_pulses_year, n_low_pulses_year**: Pulse counts (>90th, <10th percentile)
- **dur_high_pulses_year, dur_low_pulses_year**: Mean pulse duration
- **TQmean**: Percentage of days with flow above annual mean
- **Flow_Reversals**: Direction changes in flow (annual and seasonal)

### 5. Flashiness (`analyze_flashiness_trends`)
- **R-B Index**: Richards-Baker flashiness index (sum of absolute day-to-day changes / total flow)

### 6. Flow Timing (`analyze_flow_timing_trends`)
- **D05_day through D95_day**: Day of water year when cumulative flow reaches percentile
- **D25_to_D75**: Duration between 25% and 75% cumulative flow
- **Dmax**: Day of maximum flow

### 7. Q-PPT Relationships (`analyze_Q_PPT_relationships`) - *Requires Climate Data*
- Runoff ratios by season (streamflow / precipitation)
- Currently needs PPT column in input data (available in Caravan, not in raw USGS/HYDAT)

## Statistics Per Signature

Each signature has 5 associated statistics:
- **_slp**: Theil-Sen slope (robust non-parametric trend)
- **_rho**: Spearman's rank correlation coefficient
- **_pval**: P-value for trend significance
- **_mean**: Arithmetic mean across all water years
- **_median**: Median value across all water years

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

## File Structure

```
streamflowSignatures/
├── CLAUDE.md                    # This file
├── README.md                    # User documentation
├── config.R                     # Centralized configuration parameters
├── helperFunctions.R            # CANONICAL - All core functions (45 functions)
├── execute_extractStreamflowSignatureValuesAndTrends.R  # Main entry point
├── caravan_to_annualized.R      # Caravan processing
├── streamflowDataProcessing*.R  # Data ingestion scripts
├── streamflowVisualizationApp/  # Shiny dashboard
│   ├── app.R
│   └── helperFunctions.R        # App-specific utilities (S3, legends)
├── metadata/                    # Basin and gage metadata (41 files)
│   ├── conterm_*.txt           # CONUS reference gage metadata
│   ├── AKHIPR_*.txt            # Alaska reference gage metadata
│   └── Canadian_gages_goodones.csv
├── data_out/                    # Processed outputs
└── archive/                     # Deprecated helper files (DO NOT USE)
    ├── helperFunction.R
    ├── helperFunctions_sept2025.R
    ├── helperFunctions_extractStreamflowSignatureValuesAndTrends.R
    ├── helperFunctions_processRawStreamflowToParquet.R
    └── helperWrapperFunctions.R
```

## Code Status

| File | Status | Notes |
|------|--------|-------|
| `config.R` | **ACTIVE** | Centralized configuration - source before helperFunctions.R |
| `helperFunctions.R` | **CANONICAL** | All 45 core functions - use this for all development |
| `archive/*` | **ARCHIVED** | DO NOT USE - kept for reference only |

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

### Add New Signature
1. Create calculation function in `helperFunctions.R`
2. Add call to `process_signatures_from_parquet()` signature extraction section
3. Ensure output columns follow naming convention: `{metric}_{stat}` (e.g., `Qann_slp`)
4. Test with small dataset before full run

## Data Sources

### USGS (via dataRetrieval)
- Parameter code: `00060` (Discharge, cubic feet per second)
- Quality codes accepted: `A`, `A e`, `P`, `P e`
- Conversion: cfs → mm/day using drainage area

### Canadian HYDAT (via tidyhydat)
- Parameter: Flow (m³/s)
- Excludes regulated stations
- Conversion: m³/s → mm/day using drainage area

### Caravan
- NetCDF format with daily streamflow + climate variables
- Includes: PPT, SWE, temperature
- Trade-off: Shorter records (ends ~2018-2020) but has climate data

## Logging System

The project includes a structured logging system in `config.R`:

### Log Levels
- `DEBUG` (10): Verbose debugging info (per-gage, per-year details)
- `INFO` (20): Normal operation messages (default)
- `WARN` (30): Warnings that don't stop processing
- `ERROR` (40): Errors that may stop processing
- `NONE` (100): Disable all logging

### Usage
```r
source("config.R")

# Set log level
set_log_level("DEBUG")  # Show all messages
set_log_level("WARN")   # Only warnings and errors

# Enable file logging
set_log_file("logs/processing.log")

# Log messages with context
log_info("Processing started", context = "my_function")
log_warn("Missing data", context = "gage:01234567")
log_error("File not found", context = "process_signatures")
```

## Input Validation Functions

Available in `config.R` after sourcing:

| Function | Purpose |
|----------|---------|
| `validate_file_exists(path, name, ext)` | Check file exists with optional extension |
| `validate_directory(path, name, create)` | Check/create directory |
| `validate_numeric(value, name, min, max)` | Numeric range validation |
| `validate_columns(df, cols, name)` | Check data frame has required columns |
| `validate_date(value, name)` | Date parameter validation |
| `validate_gage_type(type)` | Validate gage type enum |

## Output Validation

```r
# Validate output CSV schema
result <- validate_output_schema(output_df, strict = FALSE)
# Returns list with: valid, n_metadata_cols, n_signature_cols, missing columns

# Validate single gage output
validate_gage_output(gage_row, gage_id)
# Checks for NA/Inf values, warns if >50% NA
```

## Known Issues Fixed (Jan 2026)

### Metadata Lookup Bug
- **Issue**: Data.table scoping caused all gages to receive the first row's metadata
- **Root Cause**: `metadata_lookup[gage_id == gage_id]` compared column to itself (always TRUE)
- **Fix**: Renamed `find_metadata()` parameter to `target_gage_id`
- **Location**: `helperFunctions.R` line ~3425

### Canadian Basin Area Missing
- **Issue**: Basin area was hardcoded as `NA` for Canadian stations
- **Impact**: Canadian gages had no drainage area in output
- **Fix**: Now fetches `DRAINAGE_AREA_GROSS` from `tidyhydat::hy_stations()`
- **Locations**:
  - Metadata creation functions (for new processing)
  - `process_signatures_from_parquet()` runtime fallback (for existing metadata)

## Future Development

### Climate Integration
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)

### Code Organization
- ~~Consolidate active helper file variants into canonical `helperFunctions.R`~~ DONE
- ~~Add centralized `config.R` for parameters~~ DONE
- ~~Implement structured logging~~ DONE
- ~~Fix metadata lookup bug~~ DONE
- ~~Fix Canadian basin_area bug~~ DONE
- Add unit tests

## Contact

Part of the SurfaceWaterProjections repository.

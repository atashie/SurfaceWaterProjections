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

### Design Principles

1. **Plain-English Guardrails** — Domain experts define signature methodology in `SIGNATURE_GUIDELINES.md` (auto-synced from a shared Google Doc). Code implements those definitions. This separates hydrological expertise from implementation.

2. **R Canonical, Others Follow** — All changes are made in R first, validated, then propagated to Python/Julia. Golden outputs from R serve as the reference for cross-language validation.

3. **Strict Output Schema** — The CSV output format (column names, ordering) is a contract. Downstream tools depend on exact column names. Every signature produces exactly 8 statistics via `generate_stats()`.

4. **Per-Year Quality Filtering** — Each water year is independently evaluated against data completeness thresholds (95% non-NA days, 30+ days above minimum flow). Years that fail are excluded; gages need 20+ qualifying years.

## File Structure

```
streamflowSignatures/
├── README.md                    # User entry point
├── CHANGELOG.md                 # Bug fixes, version history
├── CLAUDE.md                    # Claude Code instructions
├── .gitignore                   # Git ignore patterns
│
│ ## R Workflow (unchanged at root for current users)
├── config.R                     # Centralized configuration parameters
├── run_full_processing.R        # PRIMARY - Full signature extraction with climate data
├── run_ingest_usgs_hydat.R      # Raw USGS/HYDAT data ingestion to parquet
├── run_caravan_processing.R     # Caravan data processing
├── run_restricted_processing.R  # Restricted processing
│
├── R/                           # R canonical implementation + tests
│   ├── helperFunctions.R        # CANONICAL - All core functions (45+ functions)
│   ├── load_config.R            # Config loader
│   ├── run_conversion.R         # Daymet ZIP to Parquet conversion
│   ├── run_enrich_metadata.R    # Human interference metadata enrichment
│   ├── run_regenerate_metadata.R # Regenerate combined metadata
│   ├── precompute_cross_signature_analysis.R  # Offline computation
│   └── tests/                   # R test suite
│       ├── smoke_test.R         # Quick validation on subset (10 gages)
│       ├── smoke_test_reorganization.R
│       ├── qa_qc_signatures.R   # Output validation and QA/QC checks
│       ├── visualize_qa_qc.R    # QA/QC visualization plots
│       ├── test_climate_functions.R # Climate function tests
│       ├── test_climate_signatures.R
│       ├── generate_golden_outputs.R
│       └── verify_no_regression.R  # Golden output regression check
│
├── python/                      # Python package (production-ready)
│   ├── README.md
│   ├── pyproject.toml
│   ├── streamflow_signatures/   # 17 modules
│   └── tests/                   # Python tests
│
├── julia/                       # Julia package (production-ready)
│   ├── README.md
│   ├── Project.toml
│   ├── src/                     # 17 modules
│   └── test/                    # Julia tests
│
├── config/                      # Cross-language configuration
│   └── signatures_config.json
│
├── golden-outputs/              # R reference outputs for validation (Feb 2026, pre-Round 6)
│   ├── README.md
│   ├── streamflow_signatures_full_10feb2026.csv
│   └── combined_watershed_metadata_09feb2026.csv
│   # Note: Golden outputs pre-date Round 6 code quality changes.
│   # Round 6 produced identical correlations, so these remain valid.
│   # Refresh after next methodology change.
│
├── docs/                        # Extended documentation
│   ├── DEVELOPMENT.md           # This file
│   ├── SIGNATURES.md            # Detailed signature documentation
│   ├── SIGNATURE_GUIDELINES.md  # Collaborative guidelines (auto-synced)
│   ├── WORKFLOW_REVIEW.md       # Workflow review
│   ├── CROSS_LANGUAGE_STATUS.md # Cross-language alignment detail
│   ├── CODE_REVIEW.md          # Cross-language code review findings
│   ├── benchmarks/              # Benchmark runners and results
│   │   ├── run_python_benchmark.py
│   │   ├── run_julia_benchmark.jl
│   │   ├── run_r_benchmark.R
│   │   ├── compare_three_way.py
│   │   ├── compare_outputs.py
│   │   ├── comparison_report.md
│   │   └── diagnostics/        # Archived diagnostic scripts
│   └── plans/                   # Planning notes
│
├── claude-skill/                # Claude AI skill
│   └── streamflow-signatures.md
│
├── streamflowAndClimateVisualizationApp/  # Shiny dashboard
│   ├── app.R                    # Main Shiny application
│   └── helperFunctions.R        # App-specific utilities
│
├── metadata/                    # Basin and gage metadata (42 files)
├── archive/                     # Archived/deprecated files (DO NOT USE)
│
├── data_out/                    # Processed outputs (gitignored)
├── test_output/                 # Test outputs (gitignored)
└── logs/                        # Processing logs (gitignored)
```

## Common Tasks

### Run Signature Extraction

```r
source("config.R")              # Load configuration
source("R/helperFunctions.R")   # Load all functions

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
source("R/helperFunctions.R")
source("R/tests/qa_qc_signatures.R")
```

Or run the visualization script for diagnostic plots:

```r
source("R/tests/visualize_qa_qc.R")
# Outputs to data_out/qa_plots/
```

QA/QC checks include:
- Range validation (e.g., BFI in [0,1])
- Baseflow consistency (BFI_Eckhardt < BFI_LyneHollick)
- Elasticity constraints
- Correlation checks between related metrics

### Process Caravan Data

```r
source("R/helperFunctions.R")

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

1. **Create calculation function** in `R/helperFunctions.R`
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
   source("R/tests/smoke_test.R")
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
source("R/tests/smoke_test.R")
# Runs on 10 gages, validates output schema
```

### Climate Function Tests

```r
source("R/tests/test_climate_functions.R")
# Tests climate signatures with synthetic data
```

### Unit Tests

```r
# Run all tests in R/tests/ directory
testthat::test_dir("R/tests/")
```

## Cross-Language Benchmarks

Python and Julia implementations are validated against R using the R² of the identity line (y = x) across 5,707 common gages and 551 signature columns. This measures whether values are identical (not just correlated). Spearman rank correlation is reported as a secondary diagnostic.

### Current Status (March 2026, Post-Round 6)

| Metric | Julia | Python | R |
|--------|-------|--------|---|
| Total Time | 9.2 min | 78.9 min | 874 min* |
| Processing Rate | 13.4/s | 1.56/s | 0.11/s* |

*March 16-17, 2026 re-run. R ran concurrently with Python/Julia — timing inflated by I/O contention. Previous solo R runs: ~1-2 hours.

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 |
|------|---------|-----------|--------|-------------|
| R vs Python | 0.9968 | 1.0000 | 0.6745 | 17 |
| R vs Julia | 0.9965 | 1.0000 | 0.6745 | 19 |
| Python vs Julia | 0.9997 | 1.0000 | 0.9771 | 3 |

475 perfect (R²>=0.999), 56 good (0.99-0.999), 20 poor (<0.99). 531 of 551 columns (96.4%) have R² >= 0.99 across all 3 pairs. All 3 languages are production-ready. The 20 poor columns are trend statistics (slopes, p-values) sensitive to small numerical differences — the underlying signature means/medians are near-identical.

#### Alignment Progress (Spearman rho, cols < 0.99 — historical metric through Round 6)

| Pair | Round 0 | Round 2 | Round 3 | Round 4 | Round 5 | Round 6 | Improvement |
|------|---------|---------|---------|---------|---------|---------|-------------|
| R vs Python | 323 | 21 | 7 | 6 | **4** | **4** | 98.8% reduction |
| R vs Julia | 321 | 49 | 5 | 4 | **4** | **4** | 98.8% reduction |
| Python vs Julia | 73 | 30 | 3 | 3 | **0** | **0** | 100% reduction |

Note: Rounds 0-6 used Spearman rank correlation. Post-Round 6, the primary metric switched to identity R² (see table above), which is stricter and reveals 20 poor columns vs Spearman's 4.

#### Known Remaining Divergences (20 columns with R² < 0.99)

Under the identity R² metric, 20 columns fall below 0.99. All are trend statistics (slopes, p-values) where small numerical differences amplify:
- 4 recession pointcloud p-values: Irreducible OLS library differences (R's `lm()` QR rank-checking vs Python/Julia SVD)
- 7 FDC90th trend stats: R has 1 extra NA gage; Python/Julia have 35 extra NAs from stricter filtering
- 2 BFI_LyneHollick p-values, 2 elasticity p-values, 2 avg_storage p-values: Small filter differences amplified through trend fitting
- 2 n_low_pulses trend stats, 1 flashiness trend stat: Edge-case gage differences

Python and Julia agree near-perfectly on all columns (min Py-Jl R² = 0.977, only 3 cols < 0.99).

### Running Benchmarks

```bash
# All scripts use __file__/@__DIR__ for paths, so they work from any directory.
# Run from project root for consistency:

# R benchmark (~1-5 hours)
Rscript docs/benchmarks/run_r_benchmark.R

# Python benchmark (~70-130 min)
python docs/benchmarks/run_python_benchmark.py

# Julia benchmark (~10 min)
julia docs/benchmarks/run_julia_benchmark.jl

# Three-way comparison (R vs Python vs Julia)
python docs/benchmarks/compare_three_way.py
```

### Benchmark Files

| File | Description |
|------|-------------|
| `docs/benchmarks/run_python_benchmark.py` | Python full signature extraction |
| `docs/benchmarks/run_julia_benchmark.jl` | Julia full signature extraction |
| `docs/benchmarks/compare_three_way.py` | **PRIMARY** — Three-way comparison using R² of identity line (y=x) |
| `docs/benchmarks/comparison_report.md` | Generated comparison report |

For implementation details, alignment history, and known divergences, see [`CROSS_LANGUAGE_STATUS.md`](CROSS_LANGUAGE_STATUS.md).

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

## Parquet Data Files

### Active Parquet Files (Use These)

| File | Location | Created | Description |
|------|----------|---------|-------------|
| `combined_streamflow_data_09feb2026.parquet` | `D:/processedOuts_feb2026/` | Feb 2026 | Current streamflow parquet with bug fixes |
| `combined_watershed_metadata_09feb2026.csv` | `D:/processedOuts_feb2026/` | Feb 2026 | Corresponding metadata |
| `daymet_1980_2023.parquet` | `D:/processedOuts_feb2026/` | Feb 2026 | Climate data (PPT, temp, SWE) |

### Deprecated Parquet Files (DO NOT USE)

| File | Location | Issue |
|------|----------|-------|
| `combined_streamflow_data.parquet` | `D:/combined_streamflow_output/` | **CORRUPTED** - Contains 99999 multiplier bug for Canadian gages without basin area |
| `combined_watershed_metadata.csv` | `D:/combined_streamflow_output/` | Outdated metadata |
| `streamflowSignature_summaryData_OCT2025.csv` | `D:/combined_streamflow_output/` | Generated from corrupted parquet |

### The 99999 Bug

The October 2025 parquet was created with buggy code that applied `conversion = 99999` for Canadian gages without basin area, instead of keeping values in raw units. This resulted in Q values ~100,000x too high.

**Example (gage 08ND025):**
- Raw HYDAT: 785.6 m³/s
- Corrupted parquet: 78,557,297 "mm/day" (785.6 × 99999)
- Fixed parquet: 785.6 (raw m³/s, flagged as `area_normalized = FALSE`)

See CHANGELOG.md entry for "H5 follow-up" for full details of the fix.

## Human Interference Metadata

Watershed metadata is automatically enriched with human interference indicators when `concatenate_with_metadata()` is called during data processing.

### Data Sources

**USGS Gages (GAGES-II):**
- Location: `D:/gagesMetadata/` (configured via `GAGES_II_DIR` in config.R)
- Files: `conterm_hydromod_dams.txt`, `conterm_bas_classif.txt`, etc.
- Columns extracted: NDAMS_2009, MAJ_DDENS_2009, STOR_NID_2009, IMPNLCD06, DEVNLCD06, FRESHW_WITHDRAWAL, HYDRO_DISTURB_INDX, CLASS

**Canadian Gages (HYDAT via tidyhydat):**
- RHBN: Reference Hydrometric Basin Network designation
- REGULATED: Station regulation status from `hy_stn_regulation()`

### Unified Classification

The `human_interference_class` column provides a unified classification:
- **reference**: USGS gages with CLASS="Ref" or Canadian gages with RHBN=TRUE
- **non-reference**: USGS gages with CLASS="Non-ref" or Canadian gages with RHBN=FALSE
- **unknown**: Gages without classification data

### Manual Enrichment

To re-enrich existing metadata (one-time use):

```r
source("config.R")
source("R/helperFunctions.R")
source("R/run_enrich_metadata.R")
```

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

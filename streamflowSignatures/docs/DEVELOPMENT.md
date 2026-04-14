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
  (streamflow only)                             │    Extraction     (560+ columns)
                                                │
Daymet (climate data)  ──> Parquet Storage ─────┘
  (PPT, temp, SWE)         (joined at runtime)

Alternative Pipeline:
Caravan NetCDF ──────────────> Direct Processing ──> Caravan Output
  (bundled Q + climate)        (annualized CSVs)
```

### Design Principles

1. **Plain-English Guardrails** — Domain experts define signature methodology in `SIGNATURE_GUIDELINES.md` (auto-synced from a shared Google Doc). Code implements those definitions. This separates hydrological expertise from implementation. Note: Notion was evaluated as an alternative hosting platform (March 2026) but rejected — Notion's JS-rendered published pages are unreliable for automated fetching, and Google Docs' static HTML + permanent URLs are better suited to the auto-sync workflow.

2. **R Canonical, Others Follow** — All changes are made in R first, validated, then propagated to Python/Julia. Golden outputs from R serve as the reference for cross-language validation.

3. **Strict Output Schema** — The CSV output format (column names, ordering) is a contract. Downstream tools depend on exact column names. Every signature produces exactly 8 statistics via `generate_stats()`.

4. **Per-Year Quality Filtering** — Each water year is independently evaluated by `preprocess_daily_data()` against data completeness thresholds (>30 raw NAs, >3-day gaps, residual boundary NAs). Negative Q rejection is config-driven (`reject_negative_flow: false` by default). Constant-SD is a QA flag only. Years that fail are excluded; gages need 20+ qualifying years.

5. **Centralized NA Handling** — Missing data is handled once per gage via `preprocess_daily_data()` before any signatures are computed. This replaces ad-hoc per-signature NA handling (fillna(0), forward-fill, etc.) with a standardized pipeline: daily grid normalization → interpolation (<=3 day internal gaps) → year rejection → residual check. Configuration is in `config/signatures_config.json` under `na_handling`.

### NA Handling Architecture

```
Raw gage data (may have NAs, gaps, duplicates)
    │
    ▼
preprocess_daily_data()           ◄── Called ONCE per gage, BEFORE all signatures
    │
    ├── 1. Daily grid normalization (one row/day, sorted, unique)
    ├── 2. Raw diagnostics (NA count, max run, seasonal completeness, negative check)
    ├── 3. Year rejection (>30 raw NAs, >3-day gaps, negative Q if config-enabled)
    ├── 4. Interpolation (internal gaps <=3 days, linear, no extrapolation)
    ├── 5. Residual check (boundary NAs → reject year)
    ├── 6. PPT handling (same rules, tracked separately)
    │
    └── Returns: cleaned data, valid_years, valid_climate_years,
                 seasonal_flags, diagnostics, rejected_years
    │
    ▼
Signature functions receive clean data
    ├── Flow volumes, timing, baseflow, recession, etc.
    ├── Climate signatures use valid_climate_years subset
    └── Seasonal signatures respect seasonal_flags (incomplete → NA)
```

**Key design decisions:**
- `config/signatures_config.json` → `na_handling` section is the single source of truth
- `use_legacy_filtering: false` — new preprocessing is the default
- Negative Q rejection is conditional on `reject_negative_flow` (default: false); `negative_ann` signature counts Q<0 days instead
- Constant-SD detection is a QA flag, not a year rejection criterion
- `ice_affected_days_total` per-gage output aggregates ice-related NA days from diagnostics
- Seasonal completeness is computed from RAW observations (pre-interpolation)
- `generate_stats()` has optional `trend_completeness` / `decade_completeness` params

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
│   ├── export_hydat_metadata.R  # One-time: export HYDAT RHBN/REGULATED to CSV
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
├── rpkg/                        # R package (production-ready, mirrors Julia/Python structure)
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   ├── README.md
│   ├── R/                       # 17 modules
│   └── tests/                   # testthat tests + real-data smoke test
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
├── golden-outputs/              # R reference outputs for validation (Feb 2026)
│   ├── README.md
│   ├── streamflow_signatures_full_10feb2026.csv
│   └── combined_watershed_metadata_09feb2026.csv
│   # Note: Golden outputs pre-date April 2026 changes (trend_completeness,
│   # Section 3 signatures, recession fix). Refresh planned.
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
│   │   ├── compare_julia_vs_golden_r.py  # 6-tier Julia vs Golden R comparison
│   │   ├── compare_outputs.py
│   │   ├── build_julia_vs_golden_r_dashboard.py  # Interactive HTML dashboard
│   │   ├── build_section3_dashboard.py   # Section 3 pre/post dashboard
│   │   ├── comparison_report.md
│   │   ├── julia_vs_golden_r_summary.md  # Generated comparison report
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
├── metadata/                    # Basin and gage metadata (43 files)
│   └── canadian_hydat_interference.csv  # RHBN + REGULATED for 8,012 Canadian stations
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

### Current Status (April 2026)

| Metric | rpkg | Julia | Python | R (canonical) |
|--------|------|-------|--------|---------------|
| Total Time | 114 min | ~10 min | 78.9 min | 874 min* |
| Gages Processed | 7,369 | 7,313 | 7,369 | 5,707 |
| Processing Rate | 1.08/s | 13.4/s | 1.56/s | 0.11/s* |

*March 16-17, 2026 re-run. R canonical ran concurrently with Python/Julia — timing inflated by I/O contention. Previous solo R runs: ~1-2 hours.

#### Three-Way Identity R² (April 2026, post R canonical fixes + Julia tau-b fix)

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 |
|------|---------|-----------|--------|-------------|
| R vs Python | 0.9990 | 1.0000 | 0.7621 | 4 |
| R vs Julia | 0.9991 | 1.0000 | 0.7621 | 4 |
| Python vs Julia | 0.9999 | 1.0000 | 0.9979 | 0 |

542 perfect (R²>=0.999), 5 good (0.99-0.999), 4 poor (<0.99). 547 of 551 columns (99.3%) have R² >= 0.99 across all 3 pairs. Python-Julia agreement is perfect (0 columns below 0.99). Remaining 4 poor columns are irreducible: 2 recession Spearman p-values (exact permutation test vs t-approximation for small n), 2 FDC90th p-values (28-gage NA mismatch from floating-point precision).

**R canonical fixes that drove this improvement** (from March 2026's 20 poor → 4 poor):
1. Removed leftover min_Q filter (affected ~150 gages' year populations)
2. FDC negative Q filter (`Q >= 0` guard matching Python/Julia)
3. seasonal_flags passthrough to flow volumes and Q-PPT functions

#### rpkg vs Other Implementations (March 2026)

| Pair | Perfect (>=0.999) | Good (0.99-0.999) | Poor (<0.99) |
|------|-------------------|-------------------|-------------|
| rpkg vs Python | **527** | 18 | 6 |
| rpkg vs Julia | **515** | 28 | 8 |
| rpkg vs R (canonical) | **490** | 51 | 10 |

rpkg aligns more closely with Python/Julia than canonical R does, because rpkg uses config-driven parameters consistently. The 10 poor columns vs canonical R are explained by 4 intentional design decisions documented in [rpkg/README.md](../rpkg/README.md#intentional-differences-from-canonical-r) and [CROSS_LANGUAGE_STATUS.md](CROSS_LANGUAGE_STATUS.md#rpkg-intentional-design-decisions).

#### Known Remaining Divergences (4 columns, April 2026)

All 4 poor columns are irreducible library-level differences:
- 2 recession pointcloud p-values: OLS library differences (R's `lm()` QR rank-checking vs Python/Julia SVD)
- 2 FDC90th p-values: 28-gage NA mismatch from floating-point precision in near-zero regression

Python and Julia agree perfectly: 0 columns below 0.99 (min Py-Jl R² = 0.998).

#### Julia Post-Section 3 vs Golden R (Feb 2026)

Julia's April 2026 output includes Guidelines Section 3 changes (new signatures, recession algorithm fix, trend_completeness) that the Feb 2026 Golden R reference predates. A separate comparison pipeline (`compare_julia_vs_golden_r.py`) uses 6-tier R² classification. See `docs/benchmarks/julia_vs_golden_r_summary.md` for the full report. Key divergence drivers:

| Root Cause | Cols Affected | Type |
|------------|--------------|------|
| trend_completeness (80% gate) | 220 | Temporal gap — Golden R predates feature |
| Recession algorithm rewrite | 46 | Intentional — R/Python sync pending |
| Elasticity operator bug (fixed) | 9 | Bug fix applied |
| dur_low_pulses_all NAs | 6 | Under investigation |

#### Alignment Progress (Spearman rho, cols < 0.99 — historical metric through Round 6)

| Pair | Round 0 | Round 2 | Round 3 | Round 4 | Round 5 | Round 6 | Improvement |
|------|---------|---------|---------|---------|---------|---------|-------------|
| R vs Python | 323 | 21 | 7 | 6 | **4** | **4** | 98.8% reduction |
| R vs Julia | 321 | 49 | 5 | 4 | **4** | **4** | 98.8% reduction |
| Python vs Julia | 73 | 30 | 3 | 3 | **0** | **0** | 100% reduction |

Note: Rounds 0-6 used Spearman rank correlation. Post-Round 6, the primary metric switched to identity R² (see tables above), which is stricter.

### Running Benchmarks

```bash
# All scripts use __file__/@__DIR__ for paths, so they work from any directory.
# Run from project root for consistency:

# R canonical benchmark (~1-5 hours)
Rscript docs/benchmarks/run_r_benchmark.R

# rpkg benchmark (~2 hours)
Rscript docs/benchmarks/run_rpkg_benchmark.R

# Python benchmark (~70-130 min)
python docs/benchmarks/run_python_benchmark.py

# Julia benchmark (~10 min)
julia docs/benchmarks/run_julia_benchmark.jl

# Three-way comparison (R canonical vs Python vs Julia)
python docs/benchmarks/compare_three_way.py

# rpkg comparison (rpkg vs R canonical, Python, Julia)
python docs/benchmarks/compare_rpkg.py
```

### Benchmark Files

| File | Description |
|------|-------------|
| `docs/benchmarks/run_r_benchmark.R` | R canonical full signature extraction |
| `docs/benchmarks/run_rpkg_benchmark.R` | rpkg package full signature extraction |
| `docs/benchmarks/run_python_benchmark.py` | Python full signature extraction |
| `docs/benchmarks/run_julia_benchmark.jl` | Julia full signature extraction |
| `docs/benchmarks/compare_three_way.py` | **PRIMARY** — Three-way comparison using R² of identity line (y=x) |
| `docs/benchmarks/compare_rpkg.py` | rpkg vs all other implementations |
| `docs/benchmarks/compare_julia_vs_golden_r.py` | Julia post-Section 3 vs Golden R (Feb 2026) — 6-tier R² classification |
| `docs/benchmarks/build_comparison_dashboard.py` | Generates interactive HTML dashboard (maps + scatterplot) |
| `docs/benchmarks/build_julia_vs_golden_r_dashboard.py` | Interactive HTML: Julia vs Golden R with dual maps + scatterplot |
| `docs/benchmarks/build_section3_dashboard.py` | Section 3 pre/post comparison dashboard |
| `docs/benchmarks/comparison_report.md` | Generated comparison report |
| `docs/benchmarks/julia_vs_golden_r_summary.md` | Generated Julia vs Golden R detailed report |

For implementation details, alignment history, and known divergences, see [`CROSS_LANGUAGE_STATUS.md`](CROSS_LANGUAGE_STATUS.md).

## Data Sources

### USGS (via dataRetrieval)
- Parameter code: `00060` (Discharge, cubic feet per second)
- Quality codes accepted: `A`, `A e`, `P`, `P e`
- Conversion: cfs -> mm/day using drainage area

### Canadian HYDAT (via tidyhydat)
- Parameter: Flow (m3/s)
- Includes both regulated and unregulated stations (REGULATED flag tracked in metadata)
- Conversion: m3/s -> mm/day using drainage area
- Interference metadata (RHBN, REGULATED) exported to `metadata/canadian_hydat_interference.csv` for cross-language use

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

See docs/CHANGELOG_ARCHIVE.md entry for "H5 follow-up" for full details of the fix.

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

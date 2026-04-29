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
  (streamflow only)                             │    Extraction     (1264 columns)
                                                │
Daymet (climate data)  ──> Parquet Storage ─────┘
  (PPT, temp, SWE)         (joined at runtime)

Alternative Pipeline:
Caravan NetCDF ──────────────> Direct Processing ──> Caravan Output
  (bundled Q + climate)        (annualized CSVs)
```

### Design Principles

1. **Plain-English Guardrails** — Domain experts define signature methodology in `SIGNATURE_GUIDELINES.md` (auto-synced from a shared Google Doc). Code implements those definitions. This separates hydrological expertise from implementation. Note: Notion was evaluated as an alternative hosting platform (March 2026) but rejected — Notion's JS-rendered published pages are unreliable for automated fetching, and Google Docs' static HTML + permanent URLs are better suited to the auto-sync workflow.

2. **Julia Canonical, Others Follow** — All changes are made in Julia first, validated via benchmark (~27 min), then propagated to Python and rpkg. Golden outputs from Julia serve as the reference for cross-language validation.

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
├── R/                           # Legacy R (deprecated for signatures — ingestion utilities still active)
│   ├── helperFunctions.R        # DEPRECATED — Legacy shim (ingestion utilities only)
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
├── rpkg/                        # R port (production-ready, mirrors Julia/Python structure)
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
├── julia/                       # Julia canonical implementation (production-ready)
│   ├── README.md
│   ├── Project.toml
│   ├── src/                     # 17 modules
│   └── test/                    # Julia tests
│
├── config/                      # Cross-language configuration
│   └── signatures_config.json
│
├── golden-outputs/              # Julia canonical (April 2026) + R historical (Feb 2026)
│   ├── README.md
│   ├── streamflow_signatures_full_10feb2026.csv
│   └── combined_watershed_metadata_09feb2026.csv
│   # Julia is the canonical golden output (April 2026). R golden outputs (Feb 2026)
│   # are historical — they pre-date trend_completeness, Section 3 signatures,
│   # and recession fix.
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
│   │   ├── compare_experiment_vs_julia.py  # Parameterized experiment comparison
│   │   ├── build_experiment_vs_julia_dashboard.py  # Experiment dashboard builder
│   │   ├── build_new_vs_golden_julia_dashboard.py  # New benchmark vs golden Julia validation
│   │   ├── run_julia_benchmark_startIn1993.jl  # Experiment: WY >= 1993
│   │   ├── run_julia_benchmark_startIn1993_60pct.jl  # Experiment: WY >= 1993 + 60% filter
│   │   ├── run_julia_benchmark_startIn1993_80pct.jl  # Experiment: WY >= 1993 + 80% filter
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

**Julia (canonical)**:

```julia
using StreamflowSignatures
# See docs/benchmarks/run_julia_benchmark.jl for full pipeline
```

**R (via rpkg)**:

```r
library(streamflowsignatures)
# See docs/benchmarks/run_rpkg_benchmark.R for full pipeline
```

### Run Full Processing Pipeline

The fastest way to run a complete signature extraction:

```bash
# Julia (canonical, ~27 min for 7,313 gages):
julia docs/benchmarks/run_julia_benchmark.jl

# R legacy (deprecated — use Julia or rpkg instead):
Rscript run_full_processing.R
```

**Prerequisites:** Data paths are configured in `config/signatures_config.json` and `config.R`.

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

### Step-by-Step Process (Julia-First Workflow)

1. **Create calculation function** in the appropriate `julia/src/` module file
   - Function should accept daily data and return annual values as a DataFrame
   - Use the existing module structure (e.g., `flow_volumes.jl`, `pulses.jl`, `recession.jl`)

2. **Apply the 8-statistic rule** using `generate_stats()`:
   ```julia
   # Your function should return annual values, then call:
   stats = generate_stats(annual_df, [:metric_name], :water_year)
   # This produces 8 columns: metric_senn_slp, metric_linear_slp,
   # metric_spearman_rho, metric_spearman_pval, metric_mk_rho,
   # metric_mk_pval, metric_mean, metric_median
   ```

3. **Add call in `julia/src/signatures.jl`** in the `calculate_all_signatures()` function

4. **Register the signature** in `config/signatures_config.json` and `config.R` (for R tests):
   - Add base name to `EXPECTED_SIGNATURE_BASES` in `config.R`

5. **Run Julia benchmark** (~27 min) to validate:
   ```bash
   julia docs/benchmarks/run_julia_benchmark.jl
   ```

6. **Port to Python and rpkg**:
   - Python: Add function in the appropriate `python/streamflow_signatures/` module, call from `signatures.py`
   - rpkg: Add function in the appropriate `rpkg/R/` module, call from `signatures.R`, export in `NAMESPACE`

7. **Run cross-language comparison** to verify alignment:
   ```bash
   python docs/benchmarks/compare_three_way.py
   python docs/benchmarks/compare_rpkg.py
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

Python, rpkg, and R canonical implementations are validated against Julia golden outputs using the R² of the identity line (y = x). This measures whether implementations produce identical values (not just correlated). Spearman rank correlation is reported as a secondary diagnostic.

### Current Status (April 28, 2026 — Pettitt changepoint only)

Julia now produces 1,264 columns (656 base/metadata + 608 changepoint = 76 sigs × 8 Pettitt fields) across 7,313 gages. Python and rpkg remain at 624 (changepoint port pending). R canonical still has the old recession algorithm (46 poor columns, sync pending).

| Metric | rpkg | Julia | Python | R (canonical) |
|--------|------|-------|--------|---------------|
| Total Time | 244.6 min* | ~15 min | 150.4 min | 874 min** |
| Gages Processed | 7,313 | 7,313 | 7,313 | 5,707 |
| Signature Columns | 624 | 1,264 | 624 | 551 |
| Processing Rate | 0.5/s* | ~8/s | 0.81/s | 0.11/s** |

*rpkg April 27 re-run had I/O contention (concurrent with Python benchmark). Previous solo run: ~114 min.
**R canonical March 16-17 re-run, also with I/O contention. Previous solo R runs: ~1-2 hours.

**Pettitt changepoint signal summary**: 13.4% overall significance rate (2.7x null expectation). Strongest signal in Flashiness (19.4%), Flow Percentiles (18.4%), Baseflow (15.1%); weakest in Flow Timing (3.7%, below null). Effective independence ~17/76 signatures (77% redundancy). After BH-FDR correction, ~3.5% of evaluations survive. See `docs/SIGNATURES.md` → Changepoint Detection → Signal Robustness for full analysis.

#### Synced Implementations (rpkg, Python, Julia — 624 columns, April 27, 2026)

| Pair | Perfect (>=0.999) | Good (0.99-0.999) | Poor (<0.99) | Min R² |
|------|-------------------|-------------------|-------------|--------|
| Python vs Julia | **615** | 5 | 3 | 0.984 |
| rpkg vs Julia | **614** | 4 | 5 | 0.969 |

Python-Julia: 620 of 623 columns (99.5%) have R² >= 0.99. All 25 new recession-parameterized BFI columns are Perfect.

rpkg-Julia: 618 of 623 columns (99.2%) have R² >= 0.99. 24 of 25 new recession-parameterized BFI columns are Perfect (1 poor: `alpha_linear_spearman_pval` — Spearman p-value library difference).

#### R Canonical vs Synced Implementations (April 14, 2026)

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 |
|------|---------|-----------|--------|-------------|
| R vs Python | 0.6018 | 0.9999 | -49.72 | 49 |
| R vs Julia | 0.6022 | 0.9999 | -49.72 | 49 |
| rpkg vs R | 0.6791 | 0.9999 | -50.04 | 46 |

All 46-49 poor columns are recession metrics — R canonical still has the old look-ahead algorithm. Non-recession categories show near-perfect agreement (0 poor columns).

#### Known Remaining Divergences (7 columns max, April 2026)

All poor columns across the 3 synced implementations are irreducible library-level differences:
- 2 recession pointcloud p-values: Spearman p-value calculation differences (exact permutation vs t-approximation for small n)
- 2 FDC90th p-values: 28-gage NA mismatch from floating-point precision in near-zero regression
- 3 recession seasonality minimum: Sinusoidal fit sensitivity (Python-Julia only, R²=0.984-0.989)

#### Julia Post-Section 3 vs Golden R (Feb 2026)

Julia's April 2026 output includes Guidelines Section 3 changes (new signatures, recession algorithm fix, trend_completeness) that the Feb 2026 Golden R reference predates. A separate comparison pipeline (`compare_julia_vs_golden_r.py`) uses 6-tier R² classification. See `docs/benchmarks/julia_vs_golden_r_summary.md` for the full report. Key divergence drivers:

| Root Cause | Cols Affected | Type |
|------------|--------------|------|
| trend_completeness (80% gate) | 220 | Temporal gap — Golden R predates feature |
| Recession algorithm rewrite | 46 | Intentional — R canonical sync pending |
| Elasticity operator bug (fixed) | 9 | Bug fix applied |
| dur_low_pulses_all NAs | 6 | Under investigation |

#### Alignment Progress (historical — Spearman rho through Round 6, then identity R²)

| Pair | Round 0 | Round 6 | Post-tau-b (Apr) | Post-Section 3 (Apr 15) |
|------|---------|---------|-----------------|------------------------|
| R vs Python | 323 | **4** | 4 | 49* |
| R vs Julia | 321 | **4** | 4 | 49* |
| Python vs Julia | 73 | **0** | 0 | **3** |
| rpkg vs Python | — | — | 6 | **7** |
| rpkg vs Julia | — | — | 8 | **4** |

*R canonical recession divergence (46 cols) is intentional — old algorithm, sync pending. Non-recession: 3 poor cols (FDC90th p-values + dur_low_pulses).

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

# Primary comparisons: Python/rpkg vs Julia golden (canonical)
python docs/benchmarks/compare_python_vs_golden_julia.py
python docs/benchmarks/compare_rpkg_vs_golden_julia.py

# Julia golden dashboards
python docs/benchmarks/build_julia_golden_dashboard.py python
python docs/benchmarks/build_julia_golden_dashboard.py rpkg

# Legacy three-way comparison (historical — uses R as baseline)
python docs/benchmarks/compare_three_way.py
python docs/benchmarks/compare_rpkg.py

# Sensitivity experiments (~10-20 min each)
julia docs/benchmarks/run_julia_benchmark_startIn1993.jl
julia docs/benchmarks/run_julia_benchmark_startIn1993_60pct.jl
julia docs/benchmarks/run_julia_benchmark_startIn1993_80pct.jl

# Experiment comparisons + dashboards
python docs/benchmarks/compare_experiment_vs_julia.py startIn1993
python docs/benchmarks/compare_experiment_vs_julia.py startIn1993_60pct
python docs/benchmarks/compare_experiment_vs_julia.py startIn1993_80pct
python docs/benchmarks/build_experiment_vs_julia_dashboard.py startIn1993
python docs/benchmarks/build_experiment_vs_julia_dashboard.py startIn1993_60pct
python docs/benchmarks/build_experiment_vs_julia_dashboard.py startIn1993_80pct
```

### Benchmark Files

| File | Description |
|------|-------------|
| `docs/benchmarks/run_r_benchmark.R` | R canonical full signature extraction |
| `docs/benchmarks/run_rpkg_benchmark.R` | rpkg package full signature extraction |
| `docs/benchmarks/run_python_benchmark.py` | Python full signature extraction |
| `docs/benchmarks/run_julia_benchmark.jl` | **CANONICAL** — Julia full signature extraction (reference output) |
| `docs/benchmarks/compare_python_vs_golden_julia.py` | **PRIMARY** — Python vs Julia golden, 6-tier R² classification |
| `docs/benchmarks/compare_rpkg_vs_golden_julia.py` | **PRIMARY** — rpkg vs Julia golden, 6-tier R² classification |
| `docs/benchmarks/build_julia_golden_dashboard.py` | Interactive HTML dashboard: Python or rpkg vs Julia golden |
| `docs/benchmarks/compare_three_way.py` | Legacy — Three-way comparison (historical, uses R as baseline) |
| `docs/benchmarks/compare_rpkg.py` | Legacy — rpkg vs all other implementations |
| `docs/benchmarks/compare_julia_vs_golden_r.py` | Historical — Julia vs Golden R (Feb 2026) — 6-tier R² classification |
| `docs/benchmarks/build_comparison_dashboard.py` | Historical — Interactive HTML dashboard (maps + scatterplot) |
| `docs/benchmarks/build_julia_vs_golden_r_dashboard.py` | Historical — Julia vs Golden R with dual maps + scatterplot |
| `docs/benchmarks/build_new_vs_golden_julia_dashboard.py` | New benchmark vs golden Julia validation dashboard |
| `docs/benchmarks/build_section3_dashboard.py` | Section 3 pre/post comparison dashboard |
| `docs/benchmarks/compare_experiment_vs_julia.py` | Parameterized experiment vs Julia baseline comparison |
| `docs/benchmarks/build_experiment_vs_julia_dashboard.py` | Parameterized experiment dashboard (HTML) |
| `docs/benchmarks/run_julia_benchmark_startIn1993.jl` | Experiment: WY >= 1993 wrapper |
| `docs/benchmarks/run_julia_benchmark_startIn1993_60pct.jl` | Experiment: WY >= 1993 + 60% qualifying fraction wrapper |
| `docs/benchmarks/run_julia_benchmark_startIn1993_80pct.jl` | Experiment: WY >= 1993 + 80% qualifying fraction wrapper |
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

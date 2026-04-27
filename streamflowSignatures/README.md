# Streamflow Signatures

A framework for extracting hydrological signatures from daily streamflow data, with implementations in Julia (canonical), Python, and R. Processes data from USGS, Canadian HYDAT, and Caravan datasets to calculate 100+ metrics characterizing flow regimes, trends, and variability.

## Project Goals

**Primary:**
1. **Data Processing** — Ingest raw streamflow data (USGS, HYDAT, Caravan), clean and filter by configurable quality thresholds, collate metadata, and produce standardized parquet/CSV outputs.
2. **Signature Extraction** — Calculate 100+ hydrological signatures under strict guardrails. Methodology is updatable by domain experts via a [plain-English guidelines document](https://docs.google.com/document/u/1/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub) (`docs/SIGNATURE_GUIDELINES.md`) — hydrologists define what to calculate, and code implements those definitions.

**Secondary:**
3. **Visualization** — Interactive Shiny dashboard for exploring signatures, trends, and cross-signature relationships across thousands of gages.
4. **Cross-Language Implementations** — Python and R ports produce near-identical results to the Julia canonical implementation (99.5% of columns have R² >= 0.99). Goal: publishable packages/libraries for community use.

## Features

- **Multi-source data ingestion**: USGS NWIS, Canadian HYDAT, Caravan/HYSETS datasets
- **Comprehensive signature suite**: Flow volumes, baseflow indices, recession parameters, pulse metrics, flashiness, and flow timing
- **Robust trend analysis**: Multiple trend methods (Theil-Sen, linear regression, Mann-Kendall, Spearman) for each signature
- **Quality control**: Configurable filters for data completeness and minimum record length
- **Interactive visualization**: Shiny dashboard with maps and time series plots

## Quick Start

### Prerequisites

**Julia** (canonical, 1.9+):
```julia
cd("julia/")
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

**Python** (3.9+):
```bash
cd python/
pip install -e .
# Or: pip install numpy pandas pyarrow scipy
```

**R** (legacy):
```r
install.packages(c(
  "data.table", "arrow", "lubridate",
  "dataRetrieval", "tidyhydat",
  "zyp", "Kendall", "sf", "terra"
))

# For visualization app
install.packages(c("shiny", "leaflet", "plotly", "aws.s3"))
```

### Julia: Extract Signatures (Canonical)

```julia
using StreamflowSignatures, DataFrames

df = read_parquet("path/to/streamflow.parquet")
df = add_water_year_columns(df)  # auto-detects "Date" or "date" column
gage_data = df[df.gage_id .== "01011000", :]

# Calculate signatures (each returns a Dict{String, Float64})
all_signatures = merge(
    calculate_flow_vols_by_year(gage_data),
    analyze_flashiness_trends(gage_data),
    analyze_flow_timing_trends(gage_data),
    analyze_fdc_trends(gage_data),
    analyze_baseflow_indices(gage_data),
    analyze_recession_parameters(gage_data),
    calculate_pulse_metrics(gage_data),
)
```

**Note:** Examples use `"path/to/streamflow.parquet"` as a placeholder. Your parquet needs columns `gage_id`, `Date` (or `date`), and `Q` (mm/day). To list available gages: `unique(df.gage_id)` (Julia) or `df["gage_id"].unique()` (Python). See `docs/DEVELOPMENT.md` → "Parquet Data Files" for data locations.

See [`julia/README.md`](julia/README.md) for full API details, input format, and climate-dependent signatures.

### Python: Extract Signatures

```python
import pandas as pd
from streamflow_signatures import (
    add_water_year_columns,
    calculate_flow_vols_by_year,
    analyze_flashiness_trends,
    analyze_flow_timing_trends,
    analyze_fdc_trends,
    analyze_baseflow_indices,
    analyze_recession_parameters,
    calculate_pulse_metrics,
    # Climate-dependent (require PPT column):
    # analyze_Q_PPT_relationships, calculate_streamflow_elasticity,
    # calculate_qp_seasonality, calculate_average_storage,
)

df = pd.read_parquet("path/to/streamflow.parquet")
df = add_water_year_columns(df)  # auto-detects "Date" or "date" column
gage_data = df[df["gage_id"] == "01011000"]

# Calculate signatures (each returns a dict of metric: value)
all_signatures = {
    **calculate_flow_vols_by_year(gage_data),
    **analyze_flashiness_trends(gage_data),
    **analyze_flow_timing_trends(gage_data),
    **analyze_fdc_trends(gage_data),
    **analyze_baseflow_indices(gage_data),
    **analyze_recession_parameters(gage_data),
    **calculate_pulse_metrics(gage_data),
}
```

See [`python/README.md`](python/README.md) for full API details, input format, and climate-dependent signatures.

### R (Legacy): Extract Signatures from Parquet

```r
# Source config and helper functions
source("config.R")
source("R/helperFunctions.R")

# Process signatures
summary_output <- process_signatures_from_parquet(
  parquet_file_path = "path/to/streamflow.parquet",
  metadata_file_path = "path/to/metadata.csv",
  output_file = "path/to/output.csv",
  min_Q_value_and_days = c(0.0001, 30),  # Min flow threshold
  min_num_years = 20,                      # Min years of data
  min_frac_good_data = 0.95               # Min data completeness
)
```

### R: Process Caravan Data

```r
source("config.R")
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

## Signature Categories

| Category | Metrics | Description |
|----------|---------|-------------|
| **Flow Volumes** | Qann, Qwin, Qspr, Qsum, Qfal, Q1-Q99 | Annual/seasonal totals and percentiles |
| **FDC** | FDCall, FDC90th, FDCmid | Flow duration curve slopes |
| **Baseflow** | BFI_Eckhardt, BFI_LyneHollick | Groundwater contribution indices |
| **Recession** | log_a, b, concavity, seasonality | Recession curve parameters |
| **Pulse** | n_pulses, dur_pulses, TQmean | High/low flow event characteristics |
| **Flow Reversals** | annual, seasonal | Direction change frequency |
| **Flashiness** | flashinessRB | Richards-Baker flashiness index |
| **Flow Timing** | D5-D95_day, D25_to_D75, Dmax | Cumulative flow timing |
| **Runoff Ratios** | annual, seasonal (requires climate) | Q/P ratios by season |
| **Elasticity** | elasticity, elasticity_static (requires climate) | Streamflow sensitivity to precipitation |
| **Q-P Seasonality** | qp_slope_sd, qp_bimodality (requires climate) | Precipitation-runoff relationship |
| **Average Storage** | avg_storage (requires climate) | Mean catchment storage |

Each signature includes 8 statistics: `_senn_slp` (Theil-Sen trend), `_linear_slp` (linear trend), `_spearman_rho` (correlation), `_spearman_pval` (significance), `_mk_rho` (Mann-Kendall tau), `_mk_pval` (Mann-Kendall significance), `_mean`, `_median`.

**Note:** Climate-dependent signatures require Daymet climate data to be integrated.

## Cross-Language Alignment

Julia is the canonical implementation. Python and R are validated ports producing near-identical results:

| Metric | Julia | Python | rpkg (R port) | R (legacy) |
|--------|-------|--------|---------------|------------|
| Approx. processing time (7,313 gages) | ~10 min | ~70 min | ~2 hrs | ~1-2 hrs |
| Signature columns | 624 | 624 | 624 | 551 |
| Signature cols with R² >= 0.99 vs Julia | -- | 591/594 (99.5%)* | 590/594 (99.3%)* | 502/551 (91.1%) |

*Pre-BFI-port comparison on 594 shared columns. Post-port (624 columns) validation pending.

All implementations share configuration via `config/signatures_config.json`. The few remaining columns below 0.99 are irreducible library-level differences (Spearman p-value calculations, floating-point precision in near-zero regression). R legacy still uses the old recession algorithm (46 divergent columns).

See [`docs/CROSS_LANGUAGE_STATUS.md`](docs/CROSS_LANGUAGE_STATUS.md) for full alignment history and methodology.

## Output Format

The output CSV contains one row per gage with columns:

| Column Type | Examples | Description |
|-------------|----------|-------------|
| Metadata | gage_id, latitude, longitude, basin_area | Gage identification |
| Record Info | num_water_years, start_water_year, end_water_year | Data coverage |
| Human Interference | NDAMS_2009, HYDRO_DISTURB_INDX, CLASS, human_interference_class | Watershed disturbance indicators |
| Signatures | Qann_senn_slp, Qann_linear_slp, Qann_spearman_rho, ... | Signature statistics (8 per metric) |

## Human Interference Metadata

Watershed metadata is automatically enriched with human interference indicators:

### USGS Gages (from GAGES-II)
| Column | Description |
|--------|-------------|
| NDAMS_2009 | Number of dams upstream |
| MAJ_DDENS_2009 | Major dam density |
| STOR_NID_2009 | Total dam storage |
| IMPNLCD06 | Impervious surface (%) |
| DEVNLCD06 | Developed area (%) |
| FRESHW_WITHDRAWAL | Freshwater withdrawal |
| HYDRO_DISTURB_INDX | Composite disturbance index (0-20) |
| CLASS | Reference classification (Ref/Non-ref) |

### Canadian Gages (from HYDAT)
| Column | Description |
|--------|-------------|
| RHBN | Reference Hydrometric Basin Network flag |
| REGULATED | Station regulation status |

### Unified Classification
| Column | Values | Description |
|--------|--------|-------------|
| human_interference_class | reference, non-reference, unknown | Unified classification across data sources |

## File Structure

```
streamflowSignatures/
├── julia/                         # Julia implementation (canonical)
├── python/                        # Python implementation (port)
├── rpkg/                          # R package (port)
├── config/signatures_config.json  # Shared cross-language config (thresholds, parameters)
├── golden-outputs/                # Reference outputs for validation
├── docs/                          # Extended documentation & benchmarks
├── metadata/                      # Basin/gage metadata
├── config.R                       # R legacy configuration
├── R/helperFunctions.R            # R legacy core functions (deprecated)
├── R/tests/                       # R legacy test suite
├── run_full_processing.R          # R legacy entry point
├── run_ingest_usgs_hydat.R        # Raw data ingestion (R)
├── run_caravan_processing.R       # Caravan processing (R)
├── streamflowAndClimateVisualizationApp/  # Shiny dashboard
├── logs/                          # Processing logs (gitignored)
└── data_out/                      # Output files (gitignored)
```

## Data Sources

### USGS
Data retrieved via `dataRetrieval::readNWISdv()` with parameter code 00060 (discharge).

### Canadian HYDAT
Data retrieved via `tidyhydat::hy_daily()`. Requires local HYDAT database (downloaded automatically by tidyhydat).

### Caravan
NetCDF files containing daily streamflow plus climate variables (precipitation, snow water equivalent, temperature).

## Visualization App

Launch the Shiny dashboard:

```r
shiny::runApp("streamflowAndClimateVisualizationApp")
```

Features:
- Interactive map with signature trend visualization (2x3 grid)
- Hydrograph with weather overlay (precipitation, SWE, temperature)
- Scatter plots with custom metric selection
- Cross-signature correlation analysis (spatial vs temporal patterns)
- Multi-gage comparison

## Configuration

Key parameters in processing functions:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_num_years` | 20 | Minimum water years required |
| `min_frac_good_data` | 0.95 | Minimum % non-NA days per year |
| `min_Q_value_and_days` | c(0.0001, 30) | Min flow (mm) and days above |

Water year convention: October 1 - September 30 (Northern Hemisphere)

## Logging and Validation

The project includes structured logging and input validation via `config.R`:

```r
source("config.R")

# Set log level (DEBUG, INFO, WARN, ERROR, NONE)
set_log_level("DEBUG")

# Enable file logging
set_log_file("logs/processing.log")
```

See `config.R` for logging functions and validation utilities.

## Development

See the following documentation:
- `docs/DEVELOPMENT.md` - Architecture, workflows, and common tasks
- `docs/SIGNATURES.md` - Detailed signature documentation
- `CHANGELOG.md` - Bug fixes and version history
- `CLAUDE.md` - Claude Code instructions (rules and constraints)
- `python/README.md` - Python implementation quickstart
- `julia/README.md` - Julia implementation quickstart
- `docs/CROSS_LANGUAGE_STATUS.md` - Cross-language alignment detail

### Adding New Signatures

1. Create calculation function in `julia/src/` (canonical implementation)
2. Add function call in the signature orchestration (`julia/src/signatures.jl`)
3. Follow naming convention: `{metric}_{stat}` (e.g., `NewMetric_senn_slp`, `NewMetric_mk_pval`)
4. Propagate to Python and rpkg
5. Update documentation

## References

- **Baseflow filters**: Eckhardt (2005), Lyne & Hollick (1979)
- **Trend analysis**: Sen (1968) - Theil-Sen estimator, Mann & Kendall (1945) - Mann-Kendall test
- **Flashiness**: Baker et al. (2004) - Richards-Baker Index
- **Data sources**: USGS NWIS, Water Survey of Canada HYDAT, Caravan dataset


## Contact

For questions or collaboration, please open an issue in the repository.

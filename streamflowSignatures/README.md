# Streamflow Signatures

An R-based framework for extracting hydrological signatures from daily streamflow data. Processes data from USGS, Canadian HYDAT, and Caravan datasets to calculate 100+ metrics characterizing flow regimes, trends, and variability.

## Features

- **Multi-source data ingestion**: USGS NWIS, Canadian HYDAT, Caravan/HYSETS datasets
- **Comprehensive signature suite**: Flow volumes, baseflow indices, recession parameters, pulse metrics, flashiness, and flow timing
- **Robust trend analysis**: Theil-Sen slopes and Spearman correlations for each signature
- **Quality control**: Configurable filters for data completeness and minimum record length
- **Interactive visualization**: Shiny dashboard with maps and time series plots

## Quick Start

### Prerequisites

```r
# Install required packages
install.packages(c(
  "data.table", "arrow", "lubridate",
  "dataRetrieval", "tidyhydat",
  "zyp", "sf", "terra"
))

# For visualization app
install.packages(c("shiny", "leaflet", "plotly", "aws.s3"))
```

### Extract Signatures from Parquet

```r
# Source helper functions
source("helperFunctions.R")

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

## Signature Categories

| Category | Metrics | Description |
|----------|---------|-------------|
| **Flow Volumes** | Qann, Qwin, Qspr, Qsum, Qfal, Q05-Q95 | Annual/seasonal means and percentiles |
| **Baseflow** | BFI_Eckhardt, BFI_LyneHollick | Groundwater contribution indices |
| **Recession** | log_a, b, concavity, seasonality | Recession curve parameters |
| **Pulse** | n_pulses, dur_pulses, TQmean, reversals | High/low flow event characteristics |
| **Flashiness** | R-B Index | Richards-Baker flashiness index |
| **Flow Timing** | D05-D95_day, D25_to_D75, Dmax | Cumulative flow timing |

Each signature includes 5 statistics: `_slp` (trend), `_rho` (correlation), `_pval` (significance), `_mean`, `_median`.

## Output Format

The output CSV contains one row per gage with columns:

| Column Type | Examples | Description |
|-------------|----------|-------------|
| Metadata | gage_id, latitude, longitude, basin_area | Gage identification |
| Record Info | num_water_years, start_water_year, end_water_year | Data coverage |
| Signatures | Qann_slp, Qann_rho, Qann_pval, Qann_mean, Qann_median | Signature statistics |

## File Structure

```
streamflowSignatures/
├── helperFunctions.R              # Core functions (canonical)
├── execute_*_.R                   # Entry point scripts
├── caravan_to_annualized.R        # Caravan processing
├── streamflowVisualizationApp/    # Shiny dashboard
├── metadata/                      # Basin/gage metadata
└── data_out/                      # Output files
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
shiny::runApp("streamflowVisualizationApp")
```

Features:
- Interactive map of gage locations
- Time series plots (linear and log scale)
- Annual pattern visualization
- Multi-gage comparison

## Configuration

Key parameters in processing functions:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_num_years` | 20 | Minimum water years required |
| `min_frac_good_data` | 0.95 | Minimum % non-NA days per year |
| `min_Q_value_and_days` | c(0.0001, 30) | Min flow (mm) and days above |

Water year convention: October 1 - September 30 (Northern Hemisphere)

## Development

See `CLAUDE.md` for detailed architecture documentation, code conventions, and development guidelines.

### Adding New Signatures

1. Create calculation function in `helperFunctions.R`
2. Add function call in `process_signatures_from_parquet()`
3. Follow naming convention: `{metric}_{stat}` (e.g., `NewMetric_slp`)
4. Update documentation

## References

- **Baseflow filters**: Eckhardt (2005), Lyne & Hollick (1979)
- **Trend analysis**: Sen (1968) - Theil-Sen estimator
- **Flashiness**: Baker et al. (2004) - Richards-Baker Index
- **Data sources**: USGS NWIS, Water Survey of Canada HYDAT, Caravan dataset


## Contact

For questions or collaboration, please open an issue in the repository.

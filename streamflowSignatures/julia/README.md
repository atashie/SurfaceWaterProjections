# StreamflowSignatures.jl

Julia implementation of hydrological signature extraction from daily streamflow time series.

## Installation

```julia
# From the julia/ directory
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# Or add to your project
Pkg.add(url="https://github.com/your-org/streamflow-signatures", subdir="julia")
```

## Quick Start

```julia
using StreamflowSignatures
using DataFrames

# Load your streamflow data
df = read_parquet("path/to/streamflow.parquet")

# Add water year columns if not present (auto-detects "Date" or "date" column)
df = add_water_year_columns(df)

# Calculate signatures for a single gage
gage_data = df[df.gage_id .== "01011000", :]

# Calculate flow volume signatures (22 metrics x 8 statistics)
flow_vols = calculate_flow_vols_by_year(gage_data)

# Calculate flashiness (1 metric x 8 statistics)
flashiness = analyze_flashiness_trends(gage_data)

# Calculate timing signatures (13 metrics x 8 statistics)
timing = analyze_flow_timing_trends(gage_data)

# Calculate FDC signatures (3 metrics x 8 statistics)
fdc = analyze_fdc_trends(gage_data)

# Calculate baseflow indices (2 metrics x 8 statistics)
baseflow = analyze_baseflow_indices(gage_data)

# Calculate recession parameters (5 metrics x 8 statistics + 6 seasonality)
recession = analyze_recession_parameters(gage_data)

# Calculate pulse metrics (10 metrics x 8 statistics)
pulses = calculate_pulse_metrics(gage_data)

# Climate-dependent signatures (require PPT column in data)
# runoff_ratios = analyze_Q_PPT_relationships(gage_data)
# elasticity = calculate_streamflow_elasticity(gage_data)
# qp_season = calculate_qp_seasonality(gage_data)
# storage = calculate_average_storage(gage_data)

# Combine all signatures
all_signatures = merge(
    flow_vols, flashiness, timing, fdc,
    baseflow, recession, pulses
)
```

## Input Data Format

Streamflow data must be a DataFrame with these columns:

| Column | Type | Description |
|--------|------|-------------|
| `gage_id` | String | Unique gage identifier |
| `date` | Date | Observation date |
| `Q` | Float64 | Daily discharge in mm/day |
| `water_year` | Int | Water year (Oct 1 - Sep 30) |
| `month` | Int | Calendar month (1-12) |
| `dowy` | Int | Day of water year (1-366) |

Use `add_water_year_columns()` to add the temporal columns if you only have `date`.

## Output Statistics

Each signature metric produces 8 statistics:

| Suffix | Statistic | Description |
|--------|-----------|-------------|
| `_senn_slp` | Theil-Sen slope | Robust non-parametric trend |
| `_linear_slp` | Linear slope | OLS regression trend |
| `_spearman_rho` | Spearman rho | Monotonic correlation |
| `_spearman_pval` | Spearman p-value | Significance of correlation |
| `_mk_rho` | Mann-Kendall tau | Monotonic trend statistic |
| `_mk_pval` | Mann-Kendall p-value | Significance of trend |
| `_mean` | Mean | Arithmetic mean |
| `_median` | Median | Central value |

## Available Signatures

### Simple Signatures

- **Flow Volumes** (`calculate_flow_vols_by_year`): 22 metrics
  - Qann, Qwin, Qspr, Qsum, Qfal (seasonal totals)
  - Q1, Q5, Q10, ..., Q99 (percentiles)
  - Q95_Q10 (high-low difference)

- **Flashiness** (`analyze_flashiness_trends`): 1 metric
  - Richards-Baker flashiness index

- **Flow Timing** (`analyze_flow_timing_trends`): 13 metrics
  - D5_day, D10_day, ..., D95_day (cumulative flow timing)
  - D25_to_D75 (duration of middle 50% of flow)
  - Dmax (day of maximum discharge)

- **Flow Duration Curve** (`analyze_fdc_trends`): 3 metrics
  - FDC_all (overall slope)
  - FDC_90th (low flow slope)
  - FDC_mid (mid-range slope)

### Complex Signatures

- **Baseflow Indices** (`analyze_baseflow_indices`): 2 metrics
  - BFI_Eckhardt (Eckhardt recursive digital filter)
  - BFI_LyneHollick (Lyne-Hollick filter with 2 passes)

- **Recession Parameters** (`analyze_recession_parameters`): 5 metrics + 6 seasonality
  - log_a_pointcloud, log_a_events (recession rate parameter)
  - b_pointcloud, b_events (recession exponent)
  - concavity (curvature of recession)
  - log_a_seasonality_* (sinusoidal seasonality of recession)

- **Pulse Metrics** (`calculate_pulse_metrics`): 10 metrics
  - n_high_pulses_year, n_low_pulses_year (pulse counts)
  - dur_high_pulses_year, dur_low_pulses_year (pulse durations)
  - TQmean (percentage of days above mean)
  - Flow_Reversals_annual, _winter, _spring, _summer, _fall

### Climate-Dependent Signatures (require PPT column)

- **Runoff Ratios** (`analyze_Q_PPT_relationships`): 5 metrics
  - annual_runoff_ratio, winter/spring/summer/fall_runoff_ratio

- **Streamflow Elasticity** (`calculate_streamflow_elasticity`): 1 static + 8 trend
  - elasticity_static (single value)
  - elasticity rolling window trend statistics

- **Q-P Seasonality** (`calculate_qp_seasonality`): 2 metrics
  - qp_slope_sd (seasonal variation in Q-P relationship)
  - qp_bimodality (bimodality coefficient, >0.555 suggests seasonal)

- **Average Storage** (`calculate_average_storage`): 1 metric
  - avg_storage (catchment storage at mean discharge, mm)

## Running Tests

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

## Cross-Language Validation

This Julia implementation is validated against the canonical R implementation.
Golden reference outputs are generated by R and stored in `../golden-outputs/`.

## License

MIT

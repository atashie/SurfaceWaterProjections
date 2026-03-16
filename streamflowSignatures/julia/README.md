# StreamflowSignatures.jl

Julia implementation of hydrological signature extraction from daily streamflow time series.

## Installation

```julia
# From the julia/ directory
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# If you see a warning about manifest resolution, run:
Pkg.resolve()
```

## Quick Start

```julia
using StreamflowSignatures
using DataFrames

# Load data (auto-normalizes columns: Date->date, site_id->gage_id, prcp->PPT)
df = read_parquet("path/to/streamflow.parquet")

# Note: if your gage_id column loads as integers, convert to strings:
# df.gage_id = string.(df.gage_id)
# (The production parquet already stores gage_id as String)

df = add_water_year_columns(df)

# Optional: merge climate data for climate-dependent signatures
climate = read_parquet("path/to/daymet.parquet")
climate = add_water_year_columns(climate)
df = leftjoin(df, climate[:, [:gage_id, :date, :PPT]], on=[:gage_id, :date])

# Process a single gage
gage_data = df[df.gage_id .== "01011000", :]

# Filter to qualifying water years (3-stage filter matching R):
#   1. Min 30 days with Q > 0.0001 mm/day
#   2. Min 95% non-NA days per water year
#   3. Min 20 qualifying water years total
# Returns: (Vector{Int} of qualifying year numbers, Bool whether gage qualifies)
qual_years, qualifies = filter_qualifying_years(gage_data)
if qualifies
    qual_set = Set(qual_years)
    gage_data = gage_data[in.(gage_data.water_year, Ref(qual_set)), :]
    results = calculate_all_signatures(gage_data, "PPT" in names(gage_data))
    # results is a Dict with ~551 keys (with climate) or ~478 keys (without climate)
    # e.g. "Qann_mean" => 573.7, "Qann_senn_slp" => -0.3, ...
end
```

## Data Setup

Production data paths (edit at the top of scripts). Paths shown are Windows; use forward slashes on Linux/macOS.

| File | Default Path |
|------|-------------|
| Streamflow | `D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet` |
| Climate (Daymet) | `D:\processedOuts_feb2026\daymet_1980_2023.parquet` |
| Metadata | `D:\processedOuts_feb2026\combined_watershed_metadata_09feb2026.csv` |

These paths are specific to the development machine. Adjust for your data location.

## Smoke Test

```bash
cd julia
julia --project=. test/smoke_test.jl
```

Runs 10 hardcoded gages through the full pipeline with validation checks.
Expected runtime: ~60-120 seconds (includes JIT compilation and parquet I/O on first run).
You should see `STATUS: SMOKE TEST PASSED` if everything works.

If your data files are in a different location, edit `STREAMFLOW_PATH` and
`CLIMATE_PATH` at the top of `test/smoke_test.jl`.

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
The `read_parquet()` function auto-normalizes common column name variants (`Date`->`date`, `site_id`->`gage_id`, `prcp`->`PPT`).

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
  - FDCall (overall slope)
  - FDC90th (low flow slope)
  - FDCmid (mid-range slope)

### Complex Signatures

- **Baseflow Indices** (`analyze_baseflow_indices`): 2 metrics
  - BFI_Eckhardt (Eckhardt recursive digital filter)
  - BFI_LyneHollick (Lyne-Hollick filter with 2 passes)

- **Recession Parameters** (`analyze_recession_parameters`): 5 metrics + 6 seasonality
  - log_a_pointcloud, log_a_events (recession rate parameter)
  - b_pointcloud, b_events (recession exponent)
  - concavity (curvature of recession)
  - log_a_seasonality_* (sinusoidal seasonality of recession)

- **Pulse Metrics** (`calculate_pulse_metrics`): 14 metrics
  - n_high_pulses_year, n_high_pulses_all, n_low_pulses_year, n_low_pulses_all (pulse counts)
  - dur_high_pulses_year, dur_high_pulses_all, dur_low_pulses_year, dur_low_pulses_all (pulse durations)
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

## Individual Signature Functions

For fine-grained control, call individual signature functions instead of `calculate_all_signatures()`:

```julia
flow_vols = calculate_flow_vols_by_year(gage_data)  # returns Dict
baseflow = analyze_baseflow_indices(gage_data)       # returns Dict
```

## Running Tests

```julia
# From the julia/ directory
using Pkg
Pkg.activate(".")
Pkg.test()
```

## Production Reference

For full-scale benchmark processing (all 7,000+ gages), see `docs/benchmarks/run_julia_benchmark.jl`.

## Cross-Language Validation

This Julia implementation is validated against the canonical R implementation.
Golden reference outputs from R are stored in `../golden-outputs/`.

## License

MIT

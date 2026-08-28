# StreamflowSignatures.jl

> **Julia is the canonical implementation (April 2026). All signature development starts here.**

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

# Filter to qualifying water years — LEGACY path (3-stage filter matching R):
#   1. Min 30 days with Q > 0.0001 mm/day
#   2. Min 95% non-NA days per water year
#   3. Min 20 qualifying water years total
# Returns: (Vector{Int} of qualifying year numbers, Bool whether gage qualifies)
# NOTE: production runs use preprocess_daily_data() (exported) for year
# qualification and NA handling instead (config na_handling;
# use_legacy_filtering: false).
qual_years, qualifies = filter_qualifying_years(gage_data)
if qualifies
    qual_set = Set(qual_years)
    gage_data = gage_data[in.(gage_data.water_year, Ref(qual_set)), :]
    results = calculate_all_signatures(gage_data, "PPT" in names(gage_data))
    # results is a Dict of signature statistics. The full production pipeline
    # emits the 1,653-column product for a climate+SWE gage — this bare call
    # emits the 8-stat keys only: the Pettitt changepoint fields require the
    # changepoint=... keyword and the 14 snow metrics require snow_data=...
    # (see docs/benchmarks/run_julia_benchmark.jl for the full wiring)
    # e.g. "Qann_mean" => 573.7, "Qann_senn_slp" => -0.3, ...
end
```

## Input Data

Streamflow input is a table of `gage_id`, `date`, and daily `Q` in mm/day (see
"Input Data Format" below for the full schema). Climate-dependent signatures
additionally need `PPT` — and the snow metrics `SWE` — merged in from a climate
table keyed by gage and date (e.g. basin-aggregated Daymet). The benchmark
scripts read paths from environment variables (`STREAMFLOW_DATA_PATH`,
`STREAMFLOW_CLIMATE_PATH`, `STREAMFLOW_METADATA_PATH`). The published HISSS
data resources (daily streamflow, gage metadata, basin-aggregated Daymet
climate, and the signature products) are distributed through the HISSS
HydroShare collection:
<https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/>.

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

When a `changepoint` configuration is passed to `calculate_all_signatures()`
(as the production pipeline does), each signature additionally produces 8
Pettitt changepoint fields:

| Suffix | Description |
|--------|-------------|
| `_pettitt_cp_year` | Most likely changepoint year |
| `_pettitt_pval` | Asymptotic p-value |
| `_pettitt_pre_mean` | Mean before the changepoint |
| `_pettitt_post_mean` | Mean after the changepoint |
| `_pettitt_delta_mean` | Post minus pre mean |
| `_pettitt_pct_change` | Percent change at the changepoint |
| `_pettitt_pre_mk_pval` | Mann-Kendall p-value, pre-changepoint segment |
| `_pettitt_post_mk_pval` | Mann-Kendall p-value, post-changepoint segment |

An opt-in `AnnualCollector` (pass `collector=...`) records the per-year annual
values behind every signature's statistics; the benchmark runner drains it into
a long-format parquet with columns `gage_id, signature, water_year, value`,
written alongside the summary CSV (same schema in Python and rpkg).

## Available Signatures

### Simple Signatures

- **Flow Volumes** (`calculate_flow_vols_by_year`): 21 metrics
  - Qann, Qwin, Qspr, Qsum, Qfal (annual + seasonal totals)
  - Q1, Q5, Q10, ..., Q99 (15 percentiles)
  - Q95_Q10 (high-low difference)

- **Flashiness** (`analyze_flashiness_trends`): 1 metric
  - Richards-Baker flashiness index

- **Flow Timing** (`analyze_flow_timing_trends`): 15 metrics
  - D1_day, D5_day, ..., D99_day (13 cumulative flow timing days, per config `timing.d_percentiles`)
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

- **Recession-Parameterized Baseflow** (`analyze_baseflow_indices_with_parameters`): 2 metrics + 1 scalar
  - BFI_Eckhardt_param (Eckhardt filter with recession-derived alpha)
  - BFI_LyneHollick_param (Lyne-Hollick filter with recession-derived alpha)
  - recession_alpha_point_cloud_linear_reservoir (per-gage scalar)

- **Recession Parameters** (`analyze_recession_parameters`): 7 metrics + 6 seasonality
  - log_a_pointcloud, log_a_events (recession rate parameter; b=1 linear-reservoir convention since July 2026)
  - b_pointcloud, b_events (recession exponent, free fit)
  - concavity (curvature of recession)
  - n_recession_events (recession event count per year)
  - alpha_linear (discrete recession constant, b=1 assumption)
  - log_a_seasonality_* (sinusoidal seasonality of recession; 6 scalars)

- **Pulse Metrics** (`calculate_pulse_metrics`): 14 metrics
  - n_high_pulses_year, n_high_pulses_all, n_low_pulses_year, n_low_pulses_all (pulse counts)
  - dur_high_pulses_year, dur_high_pulses_all, dur_low_pulses_year, dur_low_pulses_all (pulse durations)
  - TQmean (percentage of days above mean)
  - Flow_Reversals_annual, _winter, _spring, _summer, _fall

- **Negative Flow Days** (`calculate_negative_days`): 1 metric
  - negative_ann (count of days with Q < 0 per year)

- **Streamflow Drought** (`calculate_drought_metrics`): 10 metrics + 5 threshold scalars
  - drought_duration_fixed_p{2,5,10,20,30} (days below the fixed threshold, 7-day-smoothed Q)
  - drought_deficit_fixed_p{2,5,10,20,30} (summed departures below the threshold)
  - drought_threshold_fixed_p{2,5,10,20,30} (per-gage threshold scalars)
  - No climate data needed; thresholds are whole-record percentiles

### Climate-Dependent Signatures (require PPT column)

- **Runoff Ratios** (`analyze_Q_PPT_relationships`): 5 metrics
  - annual_runoff_ratio, winter/spring/summer/fall_runoff_ratio

- **Streamflow Elasticity** (`calculate_streamflow_elasticity`): 2 annual metrics + 1 static + 2 diagnostics
  - elasticity_rolling (11-year rolling window, 8 stats)
  - elasticity_annual (consecutive-year differences, 8 stats)
  - elasticity_static (single value)
  - elasticity_years_total, elasticity_years_low_ppt (per-gage diagnostics)

- **Q-P Seasonality** (`calculate_qp_seasonality`): 2 metrics
  - qp_slope_sd (seasonal variation in Q-P relationship)
  - qp_bimodality (bimodality coefficient, >0.555 suggests seasonal)

- **Average Storage** (`calculate_average_storage`): 1 metric
  - avg_storage (catchment storage at mean discharge, mm)

### Snow Signatures (require Daymet SWE)

- **Snow Metrics** (`calculate_snow_metrics`): 14 metrics
  - swe_max, swe_max_dowy, snow_cover_days, snow_on_dowy, snow_off_dowy
  - melt_season_days, melt_rate, ssm, swe_apr1
  - melt_before_peak, melt_before_peak_pct, melt_before_peak_to_max_swe
  - melt_com_dowy, swe_max_to_ppt
  - Run only on an explicitly passed, SWE-valid-year-filtered `snow_data`
    frame (`calculate_all_signatures(...; snow_data=...)`); an SWE column in
    the main gage frame is never used implicitly

### Preprocessing

- **`preprocess_daily_data`** (exported): the production year-qualification and
  NA-handling preprocessor, run once per gage BEFORE any signature function —
  daily-grid normalization, <=3-day gap interpolation, year rejection, seasonal
  completeness flags, and separate PPT/SWE-valid year tracking (config
  `na_handling`).

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

Julia is the canonical implementation. Python and rpkg are validated against Julia output.
Historical golden reference outputs from R (Feb 2026) are stored in `../golden-outputs/`.

## License

MIT

# streamflowsignatures

> **Note**: R port of the Julia canonical implementation. Mirrors the Julia and Python package structure.

R package for extracting hydrological signatures from daily streamflow data.

## Installation

```r
# From the project root:
devtools::install("rpkg")

# Or load for development:
devtools::load_all("rpkg")
```

## Quick Start

```r
library(streamflowsignatures)

# Read streamflow data
df <- read_parquet("path/to/streamflow.parquet")
df <- add_water_year_columns(df)

# Filter to qualifying water years for one gage
gage_data <- df[df$gage_id == "01011000", ]
filt <- filter_qualifying_years(gage_data)

if (filt$qualifies) {
  gage_filtered <- gage_data[gage_data$water_year %in% filt$qualifying_years, ]

  # Calculate all signatures (set has_climate = TRUE if PPT column exists)
  sigs <- calculate_all_signatures(gage_filtered, has_climate = "PPT" %in% names(df))

  # sigs is a named list with ~478-551 values
  sigs$Qann_mean
  sigs$BFI_Eckhardt_mean
  sigs$elasticity_static
}
```

## Individual Signature Functions

Each function accepts a data.frame with `water_year`, `Q`, and other required
columns, and returns a named list of statistics:

| Function | Metrics | Requires Climate |
|----------|---------|-----------------|
| `calculate_flow_vols_by_year()` | 21 (Qann, seasonal, percentiles, Q95_Q10) | No |
| `analyze_fdc_trends()` | 3 (FDCall, FDC90th, FDCmid) | No |
| `analyze_flashiness_trends()` | 1 (flashinessRB) | No |
| `analyze_flow_timing_trends()` | 13 (D-days, D25_to_D75, Dmax) | No |
| `calculate_pulse_metrics()` | 14 (pulses, TQmean, reversals) | No |
| `analyze_baseflow_indices()` | 2 (BFI_Eckhardt, BFI_LyneHollick) | No |
| `analyze_recession_parameters()` | 11 (recession + seasonality) | No |
| `analyze_Q_PPT_relationships()` | 5 (runoff ratios) | Yes |
| `calculate_streamflow_elasticity()` | 1 static + 8 trend | Yes |
| `calculate_qp_seasonality()` | 2 (qp_slope_sd, qp_bimodality) | Yes |
| `calculate_average_storage()` | 1 (avg_storage) | Yes |

## Cross-Language Compatibility

This package produces near-identical results to the Python (`streamflow_signatures`)
and Julia (`StreamflowSignatures`) packages, all sharing
`config/signatures_config.json` for parameter values.

### Benchmark Results (April 2026 — post Section 3 sync)

Full benchmark: 7,313 gages, 594 signature columns. Julia canonical is the reference.

| Pair | Perfect (>=0.999) | Good (0.99-0.999) | Poor (<0.99) |
|------|-------------------|-------------------|-------------|
| rpkg vs Julia | **586** | 4 | 4 |
| rpkg vs Python | **582** | 5 | 7 |
| rpkg vs R (monolithic) | **490** | — | 46 |

rpkg aligns closely with the Julia canonical implementation.
Poor columns are irreducible library-level differences (Spearman p-values, sinusoidal fit sensitivity).

### Intentional Differences from Monolithic R

rpkg incorporates four design decisions that improve cross-language alignment
but create small divergences from the monolithic `R/helperFunctions.R` script.
Users migrating from the monolithic R workflow should be aware of these:

#### 1. FDC `min_days` = 250 (config) vs 30 (monolithic R hardcoded)

rpkg uses `fdc.min_days = 250` from `config/signatures_config.json`, matching
Julia/Python. Monolithic R hardcodes `30`. This rejects water years with 30-249
valid days from FDC90th computation, eliminating noisy slopes from data-sparse
years. **Affects**: FDC90th trend statistics (4 columns).

To match monolithic R behavior, set `fdc.min_days = 30` in the config JSON.

#### 2. BFI_LyneHollick: paired masking vs independent sums

rpkg computes BFI using paired masking: only positions where both Q and
baseflow are non-NA contribute to the ratio. Monolithic R sums numerator and
denominator independently with `na.rm=TRUE`, which includes Q at positions
where the Lyne-Hollick filter propagated NA to the baseflow. **Affects**:
BFI_LyneHollick trend p-values (2 columns).

The paired masking approach matches Julia/Python and avoids a subtle denominator
mismatch when NAs propagate through the forward-backward filter passes.

#### 3. Elasticity and Avg Storage: NA Q row handling

Monolithic R's `process_signatures_from_parquet()` removes all NA-Q rows
**before** merging climate data. rpkg passes the full dataset to signature
functions, which handle NAs internally (matching Julia/Python). This means:
- PPT on NA-Q days is included in elasticity's `P_annual` (slightly higher)
- Storage's water balance includes days where Q is replaced with 0

**Affects**: elasticity and avg_storage trend p-values (4 columns).

To match monolithic R, pre-filter your data: `df <- df[!is.na(df$Q), ]` before
calling `calculate_all_signatures()`.

## Testing

```r
devtools::test("rpkg")
```

### Real-Data Smoke Test

```r
# Requires parquet data on disk
devtools::load_all("rpkg")
source("rpkg/tests/smoke_test_real_data.R")
```

### Full Benchmark

```r
Rscript docs/benchmarks/run_rpkg_benchmark.R
```

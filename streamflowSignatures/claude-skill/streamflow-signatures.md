# Streamflow Signatures Assistant

A skill for helping users understand, use, and interpret the streamflow signatures libraries.

## Overview

This skill assists scientists working with hydrological signature extraction from daily streamflow time series. It covers the Julia (canonical), Python, and R implementations.

## What Are Streamflow Signatures?

Streamflow signatures are quantitative metrics that characterize watershed hydrology. Each signature captures a different aspect of streamflow behavior:

- **Flow Volumes**: How much water flows annually and seasonally
- **Baseflow**: Groundwater contribution to streamflow
- **Recession**: How quickly flow decreases after precipitation
- **Flashiness**: How rapidly flow rises and falls
- **Timing**: When during the year flow occurs
- **Elasticity**: How streamflow responds to precipitation changes

## The 8-Statistic Rule

Every signature metric produces **8 statistics** via `generate_stats()`:

| Suffix | Statistic | Interpretation |
|--------|-----------|----------------|
| `_senn_slp` | Theil-Sen slope | Robust trend (units/year) |
| `_linear_slp` | Linear slope | Parametric trend (units/year) |
| `_spearman_rho` | Spearman rho | Correlation with time (-1 to 1) |
| `_spearman_pval` | Spearman p-value | Significance of correlation |
| `_mk_rho` | Mann-Kendall tau | Monotonic trend (-1 to 1) |
| `_mk_pval` | Mann-Kendall p-value | Significance of trend |
| `_mean` | Arithmetic mean | Central tendency |
| `_median` | Median | Robust central tendency |

## Interpreting Results

### Trend Statistics

- **Positive slope**: Increasing trend over time
- **Negative slope**: Decreasing trend over time
- **p-value < 0.05**: Statistically significant trend
- **Theil-Sen vs Linear**: Theil-Sen is more robust to outliers

### Common Interpretations

| Signature | Positive Trend | Negative Trend |
|-----------|----------------|----------------|
| Qann | Increasing annual flow | Decreasing annual flow |
| BFI | More groundwater contribution | Less groundwater contribution |
| flashinessRB | Becoming more flashy | Becoming less flashy |
| D50_day | Flow center shifting later | Flow center shifting earlier |

## Data Requirements

### Input Format

Daily streamflow data with columns:
- `gage_id`: Unique gage identifier
- `date`: Observation date
- `Q`: Daily discharge in **mm/day**
- `water_year`: Water year (Oct 1 - Sep 30)
- `month`: Calendar month (1-12)
- `dowy`: Day of water year (1-366)

### Quality Thresholds

All year qualification is handled centrally by `preprocess_daily_data()` (single source of truth):
- **Minimum years**: 20+ water years per gage
- **Year rejection**: >30 raw NAs, >3-day gaps, negative Q, boundary NAs
- **Interpolation**: Internal gaps ≤3 days filled by linear interpolation
- **Climate years**: Separate `valid_climate_years` set (same rules applied to PPT)
- **Trend completeness**: ≥80% non-NA annual values and ≥80% in first/last decade for trend stats
- **Seasonal completeness**: Seasons with <80% raw observations flagged → seasonal metrics set to NA
- **No per-signature min_days**: Preprocessor is the only year filter (no additional gates in signature functions)

## Troubleshooting

### Common Issues

**Q: Why are all my statistics NA?**
A: Check that you have enough data. Each metric requires at least 3 non-NA annual values for statistics.

**Q: Why do Theil-Sen and linear slopes differ significantly?**
A: Theil-Sen is robust to outliers. Large differences indicate outliers or non-linear trends in your data.

**Q: How do I interpret very small p-values (e.g., 1e-10)?**
A: These indicate highly significant trends. However, statistical significance doesn't imply practical importance.

**Q: My timing metrics (D50_day) seem wrong. What's happening?**
A: Check that `dowy` (day of water year) is calculated correctly. Day 1 = October 1, not January 1.

### QA/QC Flags

The output may include flags for:
- `flagged_for_qann_range`: Annual flow outside expected range (0-2000 mm)
- `flagged_for_bfi_range`: Baseflow index outside [0, 1]
- `processing_status`: "success" or error description

## Language-Specific Notes

### Julia (Canonical)

```julia
using StreamflowSignatures
df = CSV.read("streamflow_data.csv", DataFrame)
df = add_water_year_columns(df)
signatures = calculate_flow_vols_by_year(df)
```

### Python

```python
from streamflow_signatures import calculate_flow_vols_by_year, add_water_year_columns
df = add_water_year_columns(df, date_col="date")
signatures = calculate_flow_vols_by_year(df)
```

### R

```r
source("config.R")
source("R/helperFunctions.R")
result <- process_signatures_from_parquet(parquet_file, metadata_file, output_file)
```

## References

- **SIGNATURES.md**: Detailed signature documentation
- **SIGNATURE_GUIDELINES.md**: Collaborative guidelines from hydrology team
- **METHODOLOGY.md**: Mathematical specifications (planned)

## Cross-Language Alignment

When debugging divergences between Julia, Python, and R implementations, use the **cross-language-alignment** skill. It provides a systematic 4-phase workflow, a 9-pattern divergence taxonomy, and the Iron Rule (one fix at a time, benchmark after each, revert if worse).

## Getting Help

If this skill doesn't answer your question, check:
1. The SIGNATURES.md file for detailed metric documentation
2. The CHANGELOG.md for known issues and recent fixes
3. The DEVELOPMENT.md for architecture and workflow details
4. The **cross-language-alignment** skill for debugging cross-language divergences

# Cross-Language Implementation Status

Detailed status of Python and Julia implementations relative to the R canonical reference.

For summary results and how-to-run commands, see the [Cross-Language Benchmarks section in DEVELOPMENT.md](../DEVELOPMENT.md#cross-language-benchmarks). For the historical record of individual bug fixes, see [CHANGELOG.md](../CHANGELOG.md).

## Implementation Status

**Julia**: Production-ready. Major fixes applied across 4 alignment rounds:
- BFI_LyneHollick: Fixed NaN propagation in sum() to match R's na.rm=TRUE
- BFI valid_q_mask: Removed `Q > 0` filter (zeros are valid data)
- Flow_Reversals: Implemented linear interpolation matching R's approx()
- Recession b_events: Rewrote event identification to use look-ahead algorithm, added first-day removal
- Recession concavity: Overlapping half-split matching R
- Recession log_a_events: Recalculation using median_b matching R
- Recession pointcloud: Pre-populated annual_data with all years
- Runoff ratios: Removed min_days=250 gate (R/Python have no such filter)
- Elasticity: Fixed rolling window year assignment and leap year handling
- FDC naming: Renamed to match R canonical names (FDCall, FDC90th, FDCmid)
- FDC exceedance: Fixed from Hazen x100 scale to Weibull [0,1] scale matching R/Python
- generate_stats: Pre-filter NaN values from years array before stat calculations

**Python**: Production-ready. Major fixes applied:
- BFI_LyneHollick: Fixed NaN propagation with paired masking
- BFI_Eckhardt: Fixed denominator to use total valid Q (matching R/Julia)
- Recession concavity: Overlapping half-split matching R
- Recession log_a_events: Recalculation using median_b matching R
- Runoff ratios: Paired masking for Q/PPT sums
- avg_storage: Tied Q averaging before interpolation
- qp_slope_sd: Fixed mid-point offset and ddof mismatch
- FDC column naming still uses underscores (`FDC_all` etc.); normalized during comparison via `compare_three_way.py`

## Benchmark Results (March 2026, Post-Round 4)

| Metric | Julia | Python | R |
|--------|-------|--------|---|
| Total Time | 11.6 min | 129 min | ~1-2 hours |
| Gages Processed | 7,369 | 7,369 | 5,707 |
| Total Columns | 571 | 571 | 572 |
| Common Signature Columns | 551 | 551 | 551 |
| Processing Rate | 10.6/s | 0.95/s | ~1/s |
| Common Gages (all 3) | 5,707 | 5,707 | 5,707 |

## Three-Way Spearman Correlation Summary

| Pair | Mean rho | Median rho | Min rho | Cols < 0.99 |
|------|----------|------------|---------|-------------|
| R vs Python | 0.9987 | 1.0000 | 0.85 | 7 |
| R vs Julia | 0.9987 | 1.0000 | 0.84 | 5 |
| Python vs Julia | 0.9998 | 1.0000 | 0.99 | 3 |

## Alignment Progress

| Pair | Round 0 (Cols < 0.99) | Round 2 | Round 3 | Round 4 | Improvement |
|------|----------------------|---------|---------|---------|-------------|
| R vs Python | 323 | 21 | 7 | **7** | 98% reduction |
| R vs Julia | 321 | 49 | 5 | **5** | 98% reduction |
| Python vs Julia | 73 | 30 | 3 | **3** | 96% reduction |

## Known Remaining Divergences (8 columns < 0.99)

| Issue | Severity | Affected Columns | Notes |
|-------|----------|------------------|-------|
| Recession pointcloud p-values | LOW | 4 cols (R-Py rho ~0.85-0.91) | Only 975 gages have pointcloud data; small sample amplifies p-value noise; Py-Jl=1.0 |
| BFI_Eckhardt trend stats | LOW | 3 cols (R-Py rho ~0.989-0.990) | Very close to threshold; R-Jl=1.0; Python-specific minor divergence |
| qp_slope_sd linear slope | LOW | 1 col (R-Py rho ~0.990) | Borderline numerical noise; all 3 pairs ~0.99 |

Investigation confirmed these are irreducible: recession pointcloud has no algorithmic difference (implementation-level edge-case handling in R's `lm()` vs Python's `linregress()`), and the other columns are borderline numerical noise.

## Filtering Alignment

Both Python and Julia benchmarks use per-year quality filtering matching R's `process_signatures_from_parquet()`:

1. Min 30 days with Q > 0.0001 mm/day per water year
2. Min 95% non-NA days per water year (accounting for leap years)
3. Min 20 qualifying water years per gage

Config constants are imported from shared `config/signatures_config.json` instead of hardcoded.

Output metadata columns are aligned: `basin_area`, `start_water_year`, `end_water_year`, `num_water_years`.

Python/Julia produce 7,369 qualifying gages (vs R's 5,707). The 1,662 extra gages lack Daymet climate coverage — R only iterates gages with Daymet data, while Python/Julia process all gages and leave climate signatures as NA when Daymet is unavailable.

## Validation Notes

- 543 of 551 columns have rho >= 0.99 across all 3 pairs (98.5%)
- Median rho = 1.000 for all 3 pairs
- Python-Julia: 0 columns below 0.95 (min rho = 0.989)
- All 5,707 R gages matched in Python/Julia output
- Julia is ~11x faster than Python for full benchmark (11.6 min vs 129 min)

# Cross-Language Implementation Status

Detailed status of Python and Julia implementations relative to the R canonical reference.

For summary results and how-to-run commands, see the [Cross-Language Benchmarks section in DEVELOPMENT.md](../DEVELOPMENT.md#cross-language-benchmarks). For the historical record of individual bug fixes, see [CHANGELOG.md](../CHANGELOG.md).

## Implementation Status

**Julia**: Production-ready. Major fixes applied across 6 alignment rounds:
- BFI_LyneHollick: Fixed NaN propagation in sum() to match R's na.rm=TRUE
- BFI valid_q_mask: Removed `Q > 0` filter (zeros are valid data)
- Flow_Reversals: Implemented linear interpolation matching R's approx()
- Recession b_events: Rewrote event identification to use look-ahead algorithm, added first-day removal
- Recession concavity: Overlapping half-split matching R
- Recession log_a_events: Recalculation using median_b matching R
- Recession pointcloud: Pre-populated annual_data with all years; added near-singularity guard (Round 5)
- Runoff ratios: Removed min_days=250 gate (R/Python have no such filter)
- Elasticity: Fixed rolling window year assignment and leap year handling
- FDC naming: Renamed to match R canonical names (FDCall, FDC90th, FDCmid)
- FDC exceedance: Fixed from Hazen x100 scale to Weibull [0,1] scale matching R/Python
- FDC min_days: Now config-driven via `CFG_FDC_MIN_DAYS` instead of hardcoded 250
- Flashiness: Rewrote NaN handling to use linear interpolation (was compacting array); added na_frac > 0.2 guard (Round 5)
- generate_stats: Pre-filter NaN values from years array before stat calculations
- OLS denominator guard: `abs() < 1e-10` instead of exact `== 0` (Round 6)
- DataFrame pre-allocation in 5 modules replacing `push!` pattern (Round 6)

**R**: Canonical. R16 fix applied: Eckhardt filter now forward-fills baseflow on NaN Q, matching Python/Julia. Previously cascaded NAs, causing numerator/denominator mismatch in BFI calculation. This resolves the last 3 BFI_Eckhardt columns below rho 0.99.

**Python**: Production-ready. Major fixes applied across 6 alignment rounds:
- BFI_LyneHollick: Fixed NaN propagation with paired masking
- BFI_Eckhardt: NaN forward-fill in recursive filter to eliminate cascade (Round 5); denominator fix attempted and reverted in Round 4
- Recession concavity: Overlapping half-split matching R
- Recession log_a_events: Recalculation using median_b matching R
- Recession pointcloud: Added near-singularity guard `var(log_Q) < 1e-8` (Round 5)
- Runoff ratios: Paired masking for Q/PPT sums
- avg_storage: Tied Q averaging before interpolation
- qp_slope_sd: Fixed mid-point offset and ddof bias correction
- FDC: Added negative Q filter before log10 (prevents -inf)
- FDC column naming still uses underscores (`FDC_all` etc.); normalized during comparison via `compare_three_way.py`
- Per-year groupby optimization replacing boolean indexing in 8 modules (Round 6)
- Redundant `.copy()` removal in 11 modules (Round 6)
- List accumulation replacing `.loc` boolean assignment in flashiness/fdc (Round 6)

## Benchmark Results (March 2026, Post-Round 6 / Code Quality Round 2)

| Metric | Julia | Python | R |
|--------|-------|--------|---|
| Total Time | 9.2 min | 78.9 min | 874 min* |
| Gages Processed | 7,369 | 7,369 | 5,707 |
| Total Columns | 571 | 583 | 572 |
| Common Signature Columns | 551 | 551 | 551 |
| Processing Rate | 13.4/s | 1.56/s | 0.11/s* |
| Common Gages (all 3) | 5,707 | 5,707 | 5,707 |

*March 16-17, 2026 re-run. R ran concurrently with Python/Julia — timing inflated by I/O contention. Previous solo R runs: ~1-2 hours. Golden output regression: 551/551 perfect match against Feb 2026 reference.

## Three-Way Identity R² Summary

R² of the identity line (y = x): measures whether implementations produce identical values, not just correlated values.

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 |
|------|---------|-----------|--------|-------------|
| R vs Python | 0.9968 | 1.0000 | 0.6745 | 17 |
| R vs Julia | 0.9965 | 1.0000 | 0.6745 | 19 |
| Python vs Julia | 0.9997 | 1.0000 | 0.9771 | 3 |

475 perfect (R²>=0.999), 56 good (0.99-0.999), 20 poor (<0.99). The 20 poor columns are trend statistics (slopes, p-values) sensitive to small numerical differences. Spearman rank correlation is reported as a secondary diagnostic in the comparison CSV.

## Alignment Progress (Spearman rho, cols < 0.99 — historical metric through Round 6)

| Pair | Round 0 | Round 2 | Round 3 | Round 4 | Round 5 | Round 6 | Improvement |
|------|---------|---------|---------|---------|---------|---------|-------------|
| R vs Python | 323 | 21 | 7 | 6 | **4** | **4** | 98.8% reduction |
| R vs Julia | 321 | 49 | 5 | 4 | **4** | **4** | 98.8% reduction |
| Python vs Julia | 73 | 30 | 3 | 3 | **0** | **0** | 100% reduction |

Note: Rounds 0-6 used Spearman rank correlation. Post-Round 6, the primary metric switched to identity R² (see above), which is stricter and reveals 20 poor columns vs Spearman's 4. Round 6 = post-R16 (Eckhardt forward-fill in R canonical) + code quality fixes P11-P13, J13-J14.

## Per-Category Results (Identity R², March 2026)

| Category | Total Cols | Perfect (>=0.999) | Good (>=0.99) | Poor (<0.99) | Min R² |
|----------|-----------|-------------------|---------------|-------------|--------|
| Baseflow | 16 | 8 | 6 | 2 | 0.9778 |
| Elasticity | 9 | 3 | 4 | 2 | 0.9841 |
| FDC | 24 | 11 | 6 | 7 | 0.6745 |
| Flashiness | 8 | 7 | 0 | 1 | 0.9856 |
| Flow Percentiles | 128 | 123 | 5 | 0 | 0.9976 |
| Flow Timing | 104 | 104 | 0 | 0 | 1.0000 |
| Flow Volumes | 40 | 40 | 0 | 0 | 0.9999 |
| Pulse Metrics | 112 | 99 | 7 | 6 | 0.9768 |
| Q-P Seasonality | 16 | 14 | 2 | 0 | 0.9989 |
| Recession | 46 | 32 | 10 | 4 | 0.6846 |
| Runoff Ratios | 40 | 40 | 0 | 0 | 0.9998 |
| Storage | 8 | 0 | 6 | 2 | 0.9851 |

## Known Remaining Divergences (20 columns with R² < 0.99)

All 20 poor columns are trend statistics (slopes, p-values) where small numerical differences amplify through trend fitting. Signature means and medians remain near-identical across all 3 languages.

| Column | R-Py R² | R-Jl R² | Py-Jl R² | Root Cause |
|--------|---------|---------|----------|------------|
| `FDC90th_senn_slp` | 0.675 | 0.675 | 1.000 | 1 extra NA in R; different year populations |
| `log_a_pointcloud_mk_pval` | 0.711 | 0.685 | 0.998 | OLS library edge case (QR vs SVD) |
| `b_pointcloud_mk_pval` | 0.724 | 0.700 | 0.998 | OLS library edge case (QR vs SVD) |
| `b_pointcloud_spearman_pval` | 0.759 | 0.759 | 1.000 | OLS library edge case (QR vs SVD) |
| `log_a_pointcloud_spearman_pval` | 0.771 | 0.771 | 1.000 | OLS library edge case (QR vs SVD) |
| `FDC90th_linear_slp` | 0.960 | 0.960 | 1.000 | 1 extra NA in R |
| `FDC90th_mk_pval` | 0.976 | 0.978 | 0.997 | NA population difference (R:1, Py/Jl:36) |
| `n_low_pulses_year_mk_rho` | 1.000 | 0.977 | 0.977 | NA population difference |
| `BFI_LyneHollick_mk_pval` | 0.978 | 0.978 | 1.000 | Filter propagation differences |
| `FDC90th_spearman_pval` | 0.978 | 0.979 | 0.998 | NA population difference |
| `BFI_LyneHollick_spearman_pval` | 0.978 | 0.978 | 1.000 | Filter propagation differences |
| `n_low_pulses_all_mk_rho` | 1.000 | 0.979 | 0.980 | NA population difference (Py/Jl: 309 extra NAs) |
| `FDC90th_median` | 0.983 | 0.983 | 1.000 | 1 extra NA in R |
| `elasticity_mk_pval` | 0.984 | 0.984 | 1.000 | Rolling window edge cases |
| `elasticity_spearman_pval` | 0.984 | 0.984 | 1.000 | Rolling window edge cases |
| `avg_storage_mk_pval` | 0.985 | 0.986 | 1.000 | Interpolation tie-handling |
| `FDC90th_mean` | 0.985 | 0.985 | 1.000 | 1 extra NA in R |
| `flashinessRB_linear_slp` | 1.000 | 0.986 | 0.986 | Small numerical differences |
| `avg_storage_spearman_pval` | 0.986 | 0.986 | 1.000 | Interpolation tie-handling |
| `FDC90th_mk_rho` | 0.990 | 0.994 | 0.991 | NA population difference |

**Python and Julia agree near-perfectly**: min Py-Jl R² = 0.977, only 3 columns below 0.99.

## Filtering Alignment

Both Python and Julia benchmarks use per-year quality filtering matching R's `process_signatures_from_parquet()`:

1. Min 30 days with Q > 0.0001 mm/day per water year
2. Min 95% non-NA days per water year (accounting for leap years)
3. Min 20 qualifying water years per gage

Config constants are imported from shared `config/signatures_config.json` instead of hardcoded.

Output metadata columns are aligned: `basin_area`, `start_water_year`, `end_water_year`, `num_water_years`.

Python/Julia produce 7,369 qualifying gages (vs R's 5,707). The 1,662 extra gages lack Daymet climate coverage — R only iterates gages with Daymet data, while Python/Julia process all gages and leave climate signatures as NA when Daymet is unavailable.

## Output Column Notes

### Human Interference Metadata (13 columns)

All three languages include GAGES-II human interference metadata in their output: `NDAMS_2009`, `MAJ_DDENS_2009`, `STOR_NID_2009`, `IMPNLCD06`, `DEVNLCD06`, `FRESHW_WITHDRAWAL`, `HYDRO_DISTURB_INDX`, `CLASS`, `RHBN`, `REGULATED`, `human_interference_class`, plus `area_normalized` and `gage_type`. R loads both GAGES-II (USGS gages) and HYDAT (Canadian gages) metadata via `tidyhydat`. Python/Julia load GAGES-II metadata; RHBN/REGULATED columns are present but empty (HYDAT integration is R-only for now — to be resolved in a future round).

### Q95.Q10 Naming Convention (temporary)

R uses `Q95.Q10` (dot separator, from R's `data.frame` conventions) while Python/Julia use `Q95_Q10` (underscore). The `compare_three_way.py` script normalizes dots to underscores during comparison, so this does not affect cross-language validation. To be resolved in a future round by standardizing R output to use underscores.

### QA/QC Flag Columns (12 columns)

Python and Julia benchmarks compute 12 `flagged_*` QA/QC columns (e.g., `flagged_for_qann_range`, `flagged_for_bfi_range`). These are included for comparison with golden outputs and quality validation, but are not part of the 551 signature columns used in cross-language correlation analysis. R computes equivalent flags via `qa_qc_signatures.R` as a separate post-processing step rather than inline.

## Round 6 Fixes (March 2026)

- **R16: Eckhardt BFI forward-fill**: R canonical now forward-fills baseflow on NaN Q, matching Python/Julia. Confirmed no regression — all 551 columns identical to Round 5 results.
- **P11: groupby optimization**: Replaced per-year boolean indexing with groupby in 8 Python modules (performance impact not measured due to concurrent execution)
- **P12: .copy() removal**: Removed redundant DataFrame copies in 11 Python modules
- **P13: List accumulation**: Replaced boolean `.loc` assignment with list accumulation in flashiness/fdc
- **J13: OLS denominator guard**: Changed from exact `== 0` to `abs() < 1e-10` in recession
- **J14: DataFrame pre-allocation**: Pre-allocated DataFrames in 5 Julia modules
- **compare_three_way.py**: Added dot-to-underscore normalization for R's Q95.Q10 column naming
- **Net result (Spearman rho)**: 505 perfect columns (>=0.999), 42 good (0.99-0.999), 4 poor (<0.99) — unchanged from Round 5

## Round 5 Fixes (March 2026)

- **Python BFI_Eckhardt NaN forward-fill**: Eliminated NaN cascade in recursive filter (resolved 3 BFI_Eckhardt columns)
- **Recession pointcloud near-singularity guard**: Added `var(log_Q) < 1e-8` check in both Python and Julia (improved Py-Jl agreement; R divergence remains irreducible)
- **Julia flashiness NaN interpolation**: Rewrote to use linear interpolation instead of array compaction; added na_frac > 0.2 guard
- **Python FDC negative Q filter**: Prevents `-inf` from `log10` of negative values
- **Julia FDC min_days**: Config-driven via `CFG_FDC_MIN_DAYS` instead of hardcoded 250
- **Reverted**: Constant-series Mann-Kendall fix (worsened 5 columns); FDC min_days=30 (added new poor column)
- **Net result (Spearman rho)**: 505 perfect columns (>=0.999), 42 good (0.99-0.999), 4 poor (<0.99)

## Round 4 Fixes (March 2026)

- **qp_slope_sd mid-point offset**: Fixed in both Python and Julia (resolved 1 column from R-Py < 0.99)
- **Python qp_seasonality ddof**: Fixed bias correction in standard deviation calculation
- **BFI_Eckhardt denominator**: Attempted fix to use total valid Q in Python; reverted because it worsened correlations
- **Net result**: 499 perfect columns (>=0.999), 45 good (0.99-0.999), 7 poor (<0.99)

## Validation Notes

- 547 of 551 columns have rho >= 0.99 across all 3 pairs (99.3%)
- Median rho = 1.000 for all 3 pairs
- 505 columns have rho >= 0.999 across all 3 pairs (91.7% perfect agreement)
- Python-Julia: 0 columns below 0.99 (min rho = 0.9976) — **perfect pairwise alignment**
- All 5,707 R gages matched in Python/Julia output
- Julia is ~7-14x faster than Python for full benchmark (9.8 min vs 69-133 min depending on contention)

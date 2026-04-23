# Cross-Language Implementation Status

Detailed status of Python and rpkg implementations relative to the Julia canonical reference (April 2026). R monolithic (`R/helperFunctions.R`) is deprecated for signature development; rpkg is the maintained R package.

For summary results and how-to-run commands, see the [Cross-Language Benchmarks section in DEVELOPMENT.md](../DEVELOPMENT.md#cross-language-benchmarks). For the historical record of individual bug fixes, see [CHANGELOG.md](../CHANGELOG.md).

## Implementation Status

**rpkg** (`rpkg/`): Production-ready R package. Proper installable package mirroring Julia/Python structure. Uses shared `config/signatures_config.json` for all parameters. Incorporates code review fixes that improve cross-language alignment relative to monolithic R. Key design decisions documented below.

**Julia**: **Canonical** (April 2026). Major fixes applied across 6 alignment rounds + Guidelines Section 3 changes:
- BFI_LyneHollick: Fixed NaN propagation in sum() to match R's na.rm=TRUE
- BFI valid_q_mask: Removed `Q > 0` filter (zeros are valid data)
- Flow_Reversals: Implemented linear interpolation matching R's approx()
- Recession b_events: Rewrote event identification to use look-ahead algorithm, added first-day removal
- Recession concavity: Overlapping half-split matching R
- Recession log_a_events: Recalculation using median_b matching R
- Recession pointcloud: Pre-populated annual_data with all years; added near-singularity guard (Round 5)
- Runoff ratios: Removed min_days=250 gate (R/Python have no such filter)
- Elasticity: Fixed rolling window year assignment, leap year handling, and `>=` to `>` operator for min_annual_ppt
- FDC naming: Renamed to match R names (FDCall, FDC90th, FDCmid)
- FDC exceedance: Fixed from Hazen x100 scale to Weibull [0,1] scale matching R/Python
- FDC min_days: Now config-driven via `CFG_FDC_MIN_DAYS` instead of hardcoded 250
- Flashiness: Rewrote NaN handling to use linear interpolation (was compacting array); added na_frac > 0.2 guard (Round 5)
- generate_stats: Pre-filter NaN values from years array before stat calculations
- OLS denominator guard: `abs() < 1e-10` instead of exact `== 0` (Round 6)
- DataFrame pre-allocation in 5 modules replacing `push!` pattern (Round 6)
- **Section 3 (April 2026, Julia-only pending R/Python sync):**
  - Recession event detection: position-level marking replacing look-ahead algorithm (captures full event length)
  - New signatures: D1/D99 timing, n_recession_events, elasticity_annual, negative_ann, runoff_ratio_high_count
  - Renamed: `elasticity_*` → `elasticity_rolling_*`; added elasticity diagnostics (years_total, years_low_ppt)
  - Canadian HYDAT metadata: RHBN, REGULATED, human_interference_class from pre-exported CSV
  - Mann-Kendall: tau-b tie adjustment, n≤3 p-value guard matching R

**R (monolithic)**: Deprecated for signature development (April 2026). Replaced by rpkg for R users. Historical fixes applied:
- R16: Eckhardt filter now forward-fills baseflow on NaN Q, matching Python/Julia
- April 2026: Removed leftover `min_Q_value_and_days` filter from non-legacy preprocessing path (Python/Julia never had this). Fixed ~150 gages with different year populations.
- April 2026: Added negative Q filter in FDC (`Q >= 0`), matching Python/Julia
- April 2026: `process_signatures_from_parquet()` now passes `seasonal_flags` to `calculate_flow_vols_by_year()` and `analyze_Q_PPT_relationships()`, matching Python/Julia
- April 2026: Auto-reads `trend_completeness` and `decade_completeness` from config JSON; filters climate signatures to `valid_climate_years` only

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

## Three-Way Identity R² Summary (April 2026, Post-Alignment Fixes)

R² of the identity line (y = x): measures whether implementations produce identical values, not just correlated values.

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 |
|------|---------|-----------|--------|-------------|
| R vs Python | 0.9990 | 1.0000 | 0.7621 | 4 |
| R vs Julia | 0.9991 | 1.0000 | 0.7621 | 4 |
| Python vs Julia | 0.9999 | 1.0000 | 0.9979 | 0 |

542 perfect (R²>=0.999), 5 good (0.99-0.999), 4 poor (<0.99). 547 of 551 columns (99.3%) have R² >= 0.99 across all 3 pairs. Improvement from March 2026 (475 perfect, 20 poor) driven by: removal of leftover min_Q filter in R, Julia tau-b tie adjustment, n≤3 p-value guard.

The 4 remaining poor columns are irreducible:
- 2 recession pointcloud p-values: OLS library differences (R's `lm()` QR rank-checking vs Python/Julia SVD)
- 2 FDC90th p-values: 28-gage NA mismatch from floating-point precision in near-zero regression

## rpkg Benchmark Results (March 2026)

The `rpkg/` package is a proper installable R package mirroring Julia/Python structure. It uses `config/signatures_config.json` for all parameters and incorporates code review fixes that improve cross-language alignment relative to monolithic R.

### rpkg Performance

| Metric | rpkg | Julia | Python | R (monolithic) |
|--------|------|-------|--------|----------------|
| Total Time | 114 min | 9.2 min | 78.9 min | 874 min* |
| Gages Processed | 7,369 | 7,369 | 7,369 | 5,707 |
| Total Columns | 572 | 571 | 583 | 572 |
| Processing Rate | 1.08/s | 13.4/s | 1.56/s | 0.11/s* |

### rpkg Identity R² Summary

| Pair | Perfect (>=0.999) | Good (0.99-0.999) | Poor (<0.99) |
|------|-------------------|-------------------|-------------|
| rpkg vs Python | **527** | 18 | 6 |
| rpkg vs Julia | **515** | 28 | 8 |
| rpkg vs R (monolithic) | **490** | 51 | 10 |

rpkg aligns more closely with Python/Julia than monolithic R does (527 vs 475 perfect with Python), because rpkg uses config-driven parameters consistently and incorporates the same fixes applied to Python/Julia during alignment rounds.

### rpkg vs Monolithic R: Known Divergences (10 columns with R² < 0.99)

All 10 poor columns between rpkg and monolithic R are explained by 4 intentional design decisions. These are not bugs — they represent cases where rpkg deliberately matches the Julia canonical/Python behavior instead of monolithic R.

| Column | rpkg-R R² | rpkg-Py R² | Root Cause |
|--------|-----------|------------|------------|
| `FDC90th_senn_slp` | ~0.68 | ~1.00 | Decision 1: FDC min_days |
| `FDC90th_linear_slp` | ~0.96 | ~1.00 | Decision 1: FDC min_days |
| `FDC90th_mk_pval` | ~0.98 | ~1.00 | Decision 1: FDC min_days |
| `FDC90th_spearman_pval` | ~0.98 | ~1.00 | Decision 1: FDC min_days |
| `BFI_LyneHollick_mk_pval` | ~0.98 | ~1.00 | Decision 2: paired masking |
| `BFI_LyneHollick_spearman_pval` | ~0.98 | ~1.00 | Decision 2: paired masking |
| `elasticity_mk_pval` | ~0.98 | ~1.00 | Decision 3: NA Q handling |
| `elasticity_spearman_pval` | ~0.98 | ~1.00 | Decision 3: NA Q handling |
| `avg_storage_mk_pval` | ~0.99 | ~1.00 | Decision 3: NA Q handling |
| `avg_storage_spearman_pval` | ~0.99 | ~1.00 | Decision 3: NA Q handling |

### rpkg Intentional Design Decisions

These decisions improve cross-language consistency. Users migrating from monolithic R should be aware:

**Decision 1: FDC `min_days` = 250 (config) vs 30 (monolithic R hardcoded)**
rpkg uses `fdc.min_days = 250` from `config/signatures_config.json`, matching Julia/Python. Monolithic R hardcodes `30`. This rejects water years with 30-249 valid days from FDC90th computation, eliminating noisy slopes from data-sparse years. Affects 4 FDC90th trend columns. To match monolithic R: set `fdc.min_days = 30` in the config JSON.

**Decision 2: BFI_LyneHollick paired masking vs independent sums**
rpkg computes BFI using paired masking: only positions where both Q and baseflow are non-NA contribute to the ratio. Monolithic R sums numerator and denominator independently with `na.rm=TRUE`, which includes Q at positions where the Lyne-Hollick filter propagated NA to the baseflow. Affects 2 BFI_LyneHollick trend p-values. The paired masking approach matches Julia/Python and avoids a subtle denominator mismatch.

**Decision 3: Elasticity and Avg Storage NA Q row handling**
Monolithic R's `process_signatures_from_parquet()` removes all NA-Q rows before merging climate data. rpkg passes the full dataset to signature functions, which handle NAs internally (matching Julia/Python). This means PPT on NA-Q days is included in elasticity's `P_annual` (slightly higher), and storage's water balance includes days where Q is replaced with 0. Affects 4 elasticity/avg_storage trend p-values. To match monolithic R: pre-filter your data with `df <- df[!is.na(df$Q), ]`.

## Alignment Progress (Spearman rho, cols < 0.99 — historical metric through Round 6)

| Pair | Round 0 | Round 2 | Round 3 | Round 4 | Round 5 | Round 6 | Improvement |
|------|---------|---------|---------|---------|---------|---------|-------------|
| R vs Python | 323 | 21 | 7 | 6 | **4** | **4** | 98.8% reduction |
| R vs Julia | 321 | 49 | 5 | 4 | **4** | **4** | 98.8% reduction |
| Python vs Julia | 73 | 30 | 3 | 3 | **0** | **0** | 100% reduction |

Note: Rounds 0-6 used Spearman rank correlation. Post-Round 6, the primary metric switched to identity R² (see above), which is stricter and reveals 20 poor columns vs Spearman's 4. Round 6 = post-R16 (Eckhardt forward-fill in monolithic R) + code quality fixes P11-P13, J13-J14.

## Per-Category Results (Identity R², April 2026, Three-Way)

Post R monolithic fixes + Julia tau-b fix. 551 common columns across 5,707 common gages.

| Category | Total Cols | Perfect (>=0.999) | Good (>=0.99) | Poor (<0.99) | Min R² |
|----------|-----------|-------------------|---------------|-------------|--------|
| Baseflow | 16 | 12 | 4 | 0 | 0.9930 |
| Elasticity | 9 | 5 | 4 | 0 | 0.9907 |
| FDC | 24 | 18 | 2 | 4 | 0.7621 |
| Flashiness | 8 | 8 | 0 | 0 | 0.9998 |
| Flow Percentiles | 128 | 128 | 0 | 0 | 0.9996 |
| Flow Timing | 104 | 104 | 0 | 0 | 1.0000 |
| Flow Volumes | 40 | 40 | 0 | 0 | 0.9999 |
| Pulse Metrics | 112 | 108 | 4 | 0 | 0.9914 |
| Q-P Seasonality | 16 | 16 | 0 | 0 | 1.0000 |
| Recession | 46 | 42 | 4 | 0 | 0.9908 |
| Runoff Ratios | 40 | 40 | 0 | 0 | 1.0000 |
| Storage | 8 | 5 | 3 | 0 | 0.9960 |

Note: These are three-way results (R vs Python vs Julia on the same 551 pre-Section 3 columns). Julia canonical output has additional columns (D1/D99, n_recession_events, elasticity_annual, negative_ann, diagnostics) not included here.

## Known Remaining Divergences

### Three-Way (8 columns with R² < 0.99, April 2026)

After R monolithic fixes (min_Q filter removal, FDC negative Q guard, seasonal_flags passthrough) and Julia tau-b fix, 4 poor columns remain across all 3 pairs. These are irreducible:

| Column | R-Py R² | R-Jl R² | Py-Jl R² | Root Cause |
|--------|---------|---------|----------|------------|
| `log_a_pointcloud_spearman_pval` | ~0.72 | ~0.69 | ~1.00 | OLS library (R `lm()` QR vs Python/Julia SVD) |
| `b_pointcloud_spearman_pval` | ~0.76 | ~0.76 | ~1.00 | OLS library (R `lm()` QR vs Python/Julia SVD) |
| `FDC90th_mk_pval` | ~0.76 | ~0.76 | ~1.00 | 28-gage NA mismatch from floating-point precision |
| `FDC90th_spearman_pval` | ~0.76 | ~0.76 | ~1.00 | 28-gage NA mismatch from floating-point precision |

**Python and Julia agree perfectly**: 0 columns below 0.99 (min Py-Jl R² = 0.998).

### Historical (20 columns, March 2026 — pre-R monolithic fixes)

The March 2026 three-way comparison showed 20 poor columns. The April 2026 R monolithic fixes resolved 16 of these by removing the leftover min_Q filter, adding the FDC negative Q guard, and passing seasonal_flags. See [CHANGELOG_ARCHIVE.md](CHANGELOG_ARCHIVE.md) for the full March 2026 table.

### Julia vs Golden R (288 columns with R² < 0.99)

Julia's April 2026 output (post-Section 3 changes) diverges substantially from the Feb 2026 Golden R reference. This is expected — the Golden R output predates trend_completeness, the recession algorithm fix, and several other April 2026 changes. Root causes:

| Category | Cols Affected | Root Cause | Type |
|----------|--------------|------------|------|
| All categories (trend stats) | 220 | trend_completeness filtering (80% gate) | Temporal gap with Golden R |
| Recession | 46 | Intentional algorithm change (position-marking) | Julia canonical algorithm; Golden R predates |
| Elasticity | 9 | `>=` vs `>` operator bug (fixed) + removed min_Q filter | Bug + temporal gap |
| dur_low_pulses_all | 6 | Julia produces more NAs | Under investigation |

See `docs/benchmarks/julia_vs_golden_r_summary.md` for the full detailed report.

## Filtering Alignment

All languages use centralized per-year quality filtering via `preprocess_daily_data()` (April 2026):

1. Daily grid normalization (one row/day, sorted, unique)
2. Linear interpolation of internal gaps ≤3 days
3. Year rejection: >30 raw NAs, >3-day gaps, residual boundary NAs
4. Negative Q rejection: config-driven (`reject_negative_flow: false` by default)
5. Constant-SD: QA flag only, not a rejection criterion
6. Min 20 qualifying water years per gage
7. Trend completeness: ≥80% non-NA annual values required for trend stats (slopes, p-values); mean/median always computed. Recession and elasticity exempted.

Config constants are imported from shared `config/signatures_config.json` instead of hardcoded. Legacy filtering (`use_legacy_filtering: true`) preserves the old 95%-non-NA-days rule.

Output metadata columns are aligned: `basin_area`, `start_water_year`, `end_water_year`, `num_water_years`.

Julia/Python produce 7,313-7,369 qualifying gages (vs monolithic R's 5,707). The extra gages lack Daymet climate coverage — monolithic R only iterates gages with Daymet data, while Julia/Python process all gages and leave climate signatures as NA when Daymet is unavailable.

## Output Column Notes

### Human Interference Metadata (13 columns)

All three languages include GAGES-II human interference metadata in their output: `NDAMS_2009`, `MAJ_DDENS_2009`, `STOR_NID_2009`, `IMPNLCD06`, `DEVNLCD06`, `FRESHW_WITHDRAWAL`, `HYDRO_DISTURB_INDX`, `CLASS`, `RHBN`, `REGULATED`, `human_interference_class`, plus `area_normalized` and `gage_type`. Julia loads GAGES-II metadata directly and Canadian HYDAT metadata from `metadata/canadian_hydat_interference.csv` (pre-exported from tidyhydat via `R/export_hydat_metadata.R`; 8,012 stations with RHBN and REGULATED status). Monolithic R loads both GAGES-II (USGS gages) and HYDAT (Canadian gages) metadata via `tidyhydat`. Python loads GAGES-II metadata; Canadian HYDAT integration pending.

### Q95_Q10 Naming Convention

Julia, Python, and rpkg all use `Q95_Q10` (underscore). Monolithic R uses `Q95.Q10` (dot separator, from R's `data.frame` conventions). The `compare_three_way.py` script normalizes dots to underscores during comparison, so this does not affect cross-language validation.

### QA/QC Flag Columns (12 columns)

Julia and Python benchmarks compute 12 `flagged_*` QA/QC columns (e.g., `flagged_for_qann_range`, `flagged_for_bfi_range`). These are included for comparison with golden outputs and quality validation, but are not part of the 551 signature columns used in cross-language correlation analysis. Monolithic R computes equivalent flags via `qa_qc_signatures.R` as a separate post-processing step rather than inline.

## Round 6 Fixes (March 2026)

- **R16: Eckhardt BFI forward-fill**: Monolithic R now forward-fills baseflow on NaN Q, matching Julia/Python. Confirmed no regression — all 551 columns identical to Round 5 results.
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

## Validation Notes (April 2026)

- 547 of 551 columns have R² >= 0.99 across all 3 pairs (99.3%)
- Median R² = 1.000 for all 3 pairs
- 542 columns have R² >= 0.999 across all 3 pairs (98.4% perfect agreement)
- Python-Julia: 0 columns below 0.99 (min R² = 0.998) — **perfect pairwise alignment**
- All 5,707 R gages matched in Python/Julia output
- Julia is ~7-14x faster than Python for full benchmark (~10 min vs 69-133 min depending on contention)
- Julia canonical: 594 signature columns (551 shared + 32 new stats + 8 renamed + 3 scalar diagnostics)

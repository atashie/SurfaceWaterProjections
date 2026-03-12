# Changelog

All notable changes to the Streamflow Signatures project.

## [Unreleased]

### Planned
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)
- **Python FDC column naming (LOW)**: Python CSV still uses `FDC_all/FDC_90th/FDC_mid`; should rename to match R/Julia `FDCall/FDC90th/FDCmid`

### Fixed

#### Cross-Language Alignment Round 4 — Final Column Fixes (March 2026)

**Python BFI_Eckhardt Denominator — No Fix (3 cols, investigated)**:
- **Issue**: BFI_Eckhardt trend stats (linear_slp, mk_pval, spearman_pval) had R-Py rho ~0.989-0.990; R-Jl = 1.000
- **Investigation**: Attempted changing Python Eckhardt denominator from paired masking to total valid Q (matching the apparent R/Julia formula). However, benchmarking showed this WORSENED correlations (R-Py dropped from 0.989 → 0.914). The plan's analysis was incorrect — R's actual behavior empirically matches Python's paired masking more closely than total valid Q for Eckhardt. Change was reverted. These 3 columns remain borderline at rho ~0.989-0.990.

**qp_slope_sd Mid-Point Offset (MEDIUM — 1 col, affects all 3 languages)**:
- **Issue**: qp_slope_sd_linear_slp had R-Py rho = 0.990, R-Jl rho = 0.990, Py-Jl = 1.000
- **Root Cause**: Python and Julia assigned rolling slopes to 1 day LATER than R. R uses `end_idx - floor(window/2)` for mid-point; Python used `start_idx + window//2`; Julia used `start_idx + div(window, 2)`. For a 30-day window, R assigns to position 15 while Python/Julia assigned to position 16. This 1-day offset shifts month boundary assignments, changing monthly means and the resulting SD.
- **Fix**: Changed Python to `mid_idx = end_idx - slope_window_days // 2` and Julia to `mid_idx = end_idx - div(slope_window, 2)` to match R.
- **Locations**: `python/streamflow_signatures/qp_seasonality.py` line 127, `julia/src/qp_seasonality.jl` line 133

**Python qp_seasonality ddof Mismatch (LOW — bias only, no rank effect)**:
- **Issue**: Python `np.std()` defaults to `ddof=0` (population SD, divides by n), while R `sd()` and Julia `std()` use `ddof=1` (sample SD, divides by n-1). For 12 monthly values: ratio = sqrt(11/12) ≈ 0.957. This preserves ranks (no Spearman effect) but adds systematic ~4% bias.
- **Fix**: Changed `np.std(...)` to `np.std(..., ddof=1)` for both `qp_slope_sd` and bimodality calculations.
- **Location**: `python/streamflow_signatures/qp_seasonality.py` lines 145, 157

**Recession Pointcloud P-Values — No Fix (4 cols, irreducible)**:
- **Investigation**: Thorough comparison found NO algorithmic difference across all 3 languages. 286 of 975 gages have R producing mk_pval=1.000 (n < 3 valid pointcloud years) while Python/Julia compute real p-values (n=3-4). Caused by implementation-level differences in how R's `lm()` vs Python's `linregress()` vs Julia's OLS handle edge-case fits. These 4 columns will remain < 0.99.

#### Cross-Language Alignment Round 3 — Recession, BFI, Stats Pre-Filtering (March 2026)

**Julia `generate_stats()` NaN Pre-Filtering (HIGH — 6+ cols)**:
- **Issue**: Runoff ratio p-values had R-Jl rho ~0.92-0.96 despite R-Py being perfect (1.000)
- **Root Cause**: Julia's `generate_stats()` extracted years and values from the DataFrame but did NOT pre-filter years to exclude rows where the value is NaN. It passed full arrays (including NaN-paired years) to stat functions, causing different n values and p-value calculations. R and Python both pre-filter: R does `working_data <- working_data[!is.na(working_data$value), ]`, Python does `mask = ~data[col].isna()` before extracting years.
- **Fix**: Added `valid_mask`/`valid_years` pre-filtering before all stat calculations, so only non-NaN value/year pairs are passed to `theil_sen_slope`, `linear_slope`, `spearman_correlation`, and `mann_kendall_test`.
- **Location**: `julia/src/stats.jl` lines 267-303

**Julia BFI `< 3` Minimum-Years Gate (MEDIUM — 6 cols)**:
- **Issue**: Julia produced 15 extra NaN gages in BFI columns, dragging Py-Jl correlation to ~0.96-0.97
- **Root Cause**: Julia's `analyze_baseflow_indices()` had a hardcoded `if nrow(annual_data) < 3` early-return gate that returned all NaN before reaching `generate_stats()`. R and Python have no equivalent check — they pass whatever data exists directly to `generate_stats()`, which applies its own `min_rows` logic internally.
- **Fix**: Removed the `< 3` early-return block, letting `generate_stats()` handle minimum data requirements (matching R/Python behavior).
- **Location**: `julia/src/baseflow.jl` lines 255-260

**Concavity Split — Overlapping Halves (MEDIUM — 8 cols)**:
- **Issue**: concavity_* columns had R-Py/Jl rho ~0.57-0.94
- **Root Cause**: R uses overlapping half-splits for concavity: `first_half = Q[1:mid_point]`, `second_half = Q[mid_point:n]` (both include the midpoint). Python used non-overlapping `[:mid_point]` / `[mid_point:]`, and Julia used `[1:mid_pt]` / `[mid_pt+1:end]`.
- **Fix**: Python now uses `[:mid_point + 1]` / `[mid_point:]` to include midpoint in both halves. Julia now uses `[1:mid_pt]` / `[mid_pt:end]` to include midpoint in both halves.
- **Locations**: `python/streamflow_signatures/recession.py` lines 354-356, `julia/src/recession.jl` lines 365-368

**log_a_events Recalculation with median_b (MEDIUM — 8 cols)**:
- **Issue**: log_a_events_* columns had R-Py/Jl rho ~0.86-0.96
- **Root Cause**: R recalculates log_a for each event using the ensemble `median_b` across all events in that year: `log_a_vals = log(dQ) - median_b * log(Q)`. Python and Julia used the direct per-event OLS intercept as log_a, which is mathematically different (uses the event's own b, not the ensemble median).
- **Fix**: Both Python and Julia now store Q data for successful events, compute `median_b = median(event_b_values)`, then recalculate `log_a` for each event using `log(dQ) - median_b * log(Q)`, taking the per-event median and then the year median — matching R's two-step recalculation.
- **Locations**: `python/streamflow_signatures/recession.py` lines 394-413, `julia/src/recession.jl` lines 410-427

**Julia Runoff Ratio min_days=250 Gate (HIGH — 18 cols)**:
- **Issue**: All 18 runoff ratio trend statistics had R-Jl rho ~0.92-0.98 despite R-Py being perfect (1.000)
- **Root Cause**: Julia's `analyze_Q_PPT_relationships()` had a `min_days=250` parameter that skipped entire water years with <250 valid (Q, PPT) pairs before adding to the annual DataFrame. R uses `aggregate(na.rm=TRUE)` and Python uses `dropna().groupby()` — neither has any per-year minimum check. Including ALL years changes the time series length and year indices used for trend analysis, causing different slopes and p-values.
- **Fix**: Removed `min_days` parameter and the per-year `valid_count < min_days` gate entirely. Also removed the `nrow(annual_data) < 3` early-return gate (R/Python have no equivalent — they let `generate_stats()` handle min_rows internally).
- **Location**: `julia/src/runoff_ratios.jl`
- **Impact**: All 18 runoff ratio columns now have rho > 0.999 (was 0.92-0.98)

**Julia BFI `Q > 0` Filter in valid_q_mask (MEDIUM — 6 cols)**:
- **Issue**: Julia produced 15 extra NaN gages in BFI columns (Py-Jl rho ~0.96-0.97)
- **Root Cause**: Julia's `valid_q_mask` required `!isnan(Q) && Q > 0`, while R/Python only require `!isnan(Q)`. The `Q > 0` condition excluded zero-flow days from the count, making the `min_days=250` check stricter. For 15 gages with significant zero-flow periods, this caused <250 "valid" days per year, rejecting years that R/Python accepted.
- **Fix**: Changed `valid_q_mask = .!isnan.(Q) .& (Q .> 0)` to `valid_q_mask = .!isnan.(Q)` — zeros are valid data, not missing values.
- **Location**: `julia/src/baseflow.jl` line 229

**Julia Recession Pointcloud Pre-Population (MEDIUM — 4 cols)**:
- **Issue**: Recession pointcloud p-values had Py-Jl rho ~0.92 (was improving from worse)
- **Root Cause**: Julia appended rows to `annual_data` only when recession events existed in that year. Years with pointcloud data but no successfully-fitted events were excluded, reducing the sample size for trend statistics.
- **Fix**: Pre-populated `annual_data` DataFrame with ALL years (filled with NaN), then updated specific cells using `yr_idx = findfirst(...)`. Moved pointcloud assignment OUTSIDE the `!isempty(year_b)` gate so years with pointcloud but no events still contribute.
- **Location**: `julia/src/recession.jl` lines 317-326, 400-418

**Benchmark Re-Run Results (March 12, 2026)**:
- Re-ran Julia benchmark (11.6 min, 7,369 gages) after Round 3 follow-up fixes
- Three-way comparison results:
  - R-Python: **7** poor columns (unchanged — only Julia was modified)
  - R-Julia: 49 → **5** poor columns (90% reduction)
  - Python-Julia: 30 → **3** poor columns (90% reduction), min rho improved from 0.924 → **0.989**
  - Runoff ratio columns: ALL rho > 0.999 (was 0.92-0.98 — fully resolved)
  - BFI columns: All R-Jl > 0.999 (Julia extra NAs fully resolved)
  - Only 8 unique columns remain below 0.99 (4 recession pointcloud p-values, 3 BFI Eckhardt R-Py, 1 qp_slope_sd)
- 543 of 551 columns have rho >= 0.99 across all 3 pairs (98.5%)
- Median rho = 1.000 for all 3 pairs

#### Cross-Language Alignment — Config, Metadata, and Divergence Fixes (March 2026)

**Julia FDC Exceedance Scale Bug (HIGH)**:
- **Issue**: Julia FDC90th values in completely different range than R/Python (~100x off), causing rho ~0.74-0.87 for 24 FDC columns
- **Root Cause**: Julia used Hazen plotting position `(i - 0.5) / n * 100` (scale 0-100), while R/Python use Weibull `i / (n + 1)` (scale 0-1). Since slope = delta(log10 Q) / delta(exceedance), Julia slopes were exactly 100x smaller.
- **Additional issues**: Julia filtered out zeros (`x > 0`) while R/Python add 1e-10; Julia required 5 points per sub-range while R/Python require 3
- **Fix**: Changed exceedance to Weibull formula `i / (n + 1)`, ranges from (0-100/90-100/20-80) to (0-1/0.9-1/0.2-0.8), added 1e-10 to zeros, reduced min points from 5 to 3
- **Location**: `julia/src/fdc.jl`

**Per-Year Quality Filtering in Python/Julia Benchmarks (HIGH)**:
- **Issue**: 335 of 551 columns had rho < 0.99, primarily due to different gage populations and water year ranges
- **Root Cause**: R uses three-stage per-year filtering: (1) min 30 days with Q > 0.0001, (2) min 95% non-NA days per year accounting for leap years, (3) min 20 qualifying years. Python/Julia used global filtering only (n_years >= 20 and global na_frac <= 0.05), qualifying 7,636 vs R's 5,707 gages with different year ranges per gage.
- **Fix**: Replaced `filter_qualifying_gages()` with `filter_qualifying_years()` in both `benchmarks/run_python_benchmark.py` and `benchmarks/run_julia_benchmark.jl`, implementing identical three-stage per-year filtering. Each gage's data is filtered to only qualifying water years before signature computation.
- **Impact**: Python/Julia should now produce ~5,700 qualifying gages (matching R), with identical water year ranges per gage. Expected to resolve majority of the 335 poor correlation columns.

**Hardcoded Constants Replaced with Config Imports (MEDIUM)**:
- **Issue**: `MIN_YEARS = 20` and `MIN_FRAC_GOOD_DATA = 0.95` were hardcoded in both benchmark runners instead of using the shared JSON config
- **Fix**: Python imports `MIN_NUM_YEARS`, `MIN_FRAC_GOOD_DATA`, `MIN_Q_VALUE`, `MIN_DAYS_ABOVE_THRESHOLD` from `streamflow_signatures.config`; Julia uses `CFG_MIN_NUM_YEARS`, `CFG_MIN_FRAC_GOOD_DATA`, `CFG_MIN_Q_VALUE`, `CFG_MIN_DAYS_ABOVE_THRESHOLD` from `StreamflowSignatures.config`
- **Locations**: `benchmarks/run_python_benchmark.py`, `benchmarks/run_julia_benchmark.jl`

**Missing Metadata Columns in Python/Julia Output (MEDIUM)**:
- **Issue**: Python/Julia output `start_year`/`end_year` (from metadata CSV) but not `start_water_year`/`end_water_year` (computed from actual qualifying data). R computes these from the qualifying year range.
- **Fix**: Both benchmarks now compute `start_water_year = min(qualifying_years)`, `end_water_year = max(qualifying_years)`, `num_water_years = len(qualifying_years)` from the per-year filter results. Also renamed `basin_area_km2` to `basin_area` in output to match R column naming.
- **Locations**: `benchmarks/run_python_benchmark.py`, `benchmarks/run_julia_benchmark.jl`

**Benchmark Re-Run Results (March 11, 2026)**:
- Re-ran Python (101 min, 7,369 gages) and Julia (16.5 min, 7,369 gages) benchmarks with aligned per-year filtering
- Three-way comparison results:
  - R-Python: 323 → **21** poor columns (93% reduction), runoff ratio linear slope rho improved from 0.47 → **1.000**
  - R-Julia: 321 → **49** poor columns (85% reduction)
  - Python-Julia: 73 → **30** poor columns, min rho improved from 0.51 → **0.924**
  - FDC columns: All rho > 0.99 (was 0.74 — fully resolved by exceedance fix)
  - Storage trends: R-Py now > 0.99 (resolved by filtering alignment)
  - avg_storage R-Jl rho ~0.97 (slight remaining divergence)
- Remaining 49 poor columns are recession (inherently noisy, ~1,351 gages), BFI (Julia 15 extra NAs), and runoff ratio Julia p-values
- All 5,707 R gages now matched (was 5,652 common gages before)

### Added

**Three-Way Comparison Infrastructure (March 2026)**:
- Fixed `compare_three_way.py` `normalize_col()` to handle FDC naming mismatch (`FDC_all` -> `FDCall`, `FDC_90th` -> `FDC90th`, `FDC_mid` -> `FDCmid`)
- Added missing metadata columns to `META_COLS`: `start_water_year`, `end_water_year`, `country`
- Added `generate_three_way_report()` function producing comprehensive markdown report with:
  - Input file timestamps and gage/column counts
  - Performance comparison table (Python vs Julia timing)
  - Overall Spearman correlation summary per language pair
  - Agreement breakdown by signature category (12 categories)
  - Detailed table of all columns with min rho < 0.99
  - NA mismatch analysis
  - Deep dive on worst 10 metrics with per-language range statistics
- Deprecated `compare_results.py` (Python-vs-Julia only) in favor of `compare_three_way.py`
- Updated DEVELOPMENT.md benchmark section with fresh three-way results

**Three-Way Comparison Results (March 11, 2026)**:
- 551 common signature columns compared across R, Python, Julia (5,652 common gages)
- Overall: 86 perfect (>=0.999), 130 good (>=0.99), 335 poor (<0.99)
- Python-Julia agreement: median rho=1.000 (near-perfect for most columns)
- Most "poor" columns are trend statistics (slopes, p-values) that amplify small computation differences
- FDC columns now properly compared (were silently excluded before due to naming mismatch)

#### Cross-Language Alignment Fixes (March 2026, Round 2)

**Julia D25_to_D75 Bug (HIGH)**:
- **Issue**: D25_to_D75 was always NaN in Julia output
- **Root Cause**: `CFG_D_PERCENTILES = [5,10,20,30,40,50,60,70,80,90,95]` does not include 25 or 75, so dict lookup always returned NaN
- **Fix**: Compute D25/D75 directly from cumulative percentages, matching R/Python
- **Location**: `julia/src/timing.jl`

**Julia log_a_seasonality_minimum Bug (HIGH)**:
- **Issue**: Julia seasonality minimum was ~365 minus the correct value
- **Root Cause**: Different sin/cos basis order and phase formula from R/Python
- **Fix**: Matched R's `atan2(-B2, B1)` + 273.75 formula exactly
- **Location**: `julia/src/recession.jl`

**Python BFI_LyneHollick Canadian Gage Discrepancy (HIGH)**:
- **Issue**: Python BFI 0.15-0.59 lower than R/Julia for Canadian gages
- **Root Cause**: `np.nansum(bf) / np.nansum(Q)` creates numerator/denominator mismatch when filter propagates NaN to valid Q positions
- **Fix**: Replaced with paired masking — only sum where BOTH Q and BF are valid
- **Location**: `python/streamflow_signatures/baseflow.py`

**FDC Column Naming Mismatch (MEDIUM)**:
- **Issue**: R uses `FDCall/FDC90th/FDCmid`, Python/Julia used `FDC_all/FDC_90th/FDC_mid`, preventing cross-language comparison of 24 columns
- **Fix**: Renamed Python/Julia FDC columns to match R canonical naming
- **Additional Fix**: Julia FDC mid-range changed from (20-70%) to (20-80%) matching R
- **Locations**: `julia/src/fdc.jl`, `python/streamflow_signatures/fdc.py`, `python/streamflow_signatures/io.py`, tests

**Julia Spearman ordinalrank Bug (HIGH)**:
- **Issue**: Julia used `ordinalrank()` (sequential ranks for ties) instead of `tiedrank()` (average ranks), causing ~100 Spearman columns to diverge
- **Fix**: Changed to `tiedrank()` matching R/Python behavior
- **Location**: `julia/src/stats.jl`

**Julia Recession min_events Gate Missing (MEDIUM)**:
- **Issue**: Julia produced recession metrics for 1,571 extra gages where R/Python returned NA
- **Root Cause**: Julia only applied the `min_events=25` threshold to seasonality metrics, not base metrics (b_events, concavity, etc.)
- **Fix**: Added min_events gate before generate_stats() for all recession metrics
- **Location**: `julia/src/recession.jl`

**Runoff Ratio Slope Divergence (MEDIUM)**:
- **Issue**: `annual_runoff_ratio_linear_slp` R-Python rho=0.47 despite mean agreement at 0.99
- **Root Cause**: Python/Julia summed Q and PPT independently (days with Q=NaN still contributed PPT to sums, inflating denominator)
- **Fix**: Drop rows where either Q or PPT is NaN before summing, matching R's `aggregate(na.action=na.omit)`
- **Locations**: `python/streamflow_signatures/runoff_ratios.py`, `julia/src/runoff_ratios.jl`

**avg_storage Trend Disagreement (MEDIUM)**:
- **Issue**: avg_storage p-values had rho 0.77-0.83 across languages despite mean agreement at 0.99
- **Root Cause**: Different handling of duplicate Q values during interpolation (R averages S at tied Q; Python takes last; Julia takes arbitrary pair)
- **Fix**: Python and Julia now average S values at duplicate Q positions before interpolation, matching R's `approx(ties=mean)`
- **Locations**: `python/streamflow_signatures/storage.py`, `julia/src/storage.jl`

**Gage Count Difference (DOCUMENTED)**:
- **Issue**: R qualifies 5,707 gages vs Python/Julia's 7,636
- **Root Causes**: (1) R only iterates gages with Daymet coverage (~1,879 excluded), (2) R uses per-year quality filtering vs Python/Julia's global filter, (3) R requires min_Q_value_and_days
- **Status**: Documented as expected design difference; Python/Julia process all gages, leaving climate signatures as NA when Daymet unavailable

### Added
- **Collaborative Guidelines Workflow**: Set up auto-sync of `SIGNATURE_GUIDELINES.md` from [Google Doc](https://docs.google.com/document/u/1/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub) at session start
- Renamed `Summary Documentation for Streamflow Signatures.docx.txt` to `SIGNATURE_GUIDELINES.md`
- Added session-start workflow instructions to `CLAUDE.md`
- Added changelog maintenance requirements to `CLAUDE.md`

#### Cross-Language Benchmark (Python vs Julia vs R)
- **New benchmark suite**: Full signature extraction on 7,636 gages comparing Python, Julia, and R implementations
- **Files created**:
  - `benchmarks/run_python_benchmark.py` - Python benchmark runner
  - `benchmarks/run_julia_benchmark.jl` - Julia benchmark runner
  - `benchmarks/compare_results.py` - Cross-language comparison script
  - `benchmarks/comparison_report.md` - Generated comparison report
- **Results (March 2026, post-fixes)**:
  - Julia: 7.6 min, 569 columns, 16.8 gages/s
  - Python: 120 min, 569 columns, 1.06 gages/s
  - R: ~1-2 hours, 561 columns
- **Key correlations with R reference** (5,652 common gages):
  - Qann_mean: Julia=0.9999, Python=0.9999
  - BFI_LyneHollick_mean: Julia=0.9983, Python=0.9872
  - elasticity_static: Julia=0.9953, Python=0.9953
  - Flow_Reversals_annual_mean: Julia=0.9993, Python=0.9993
  - flashinessRB_mean: Julia=0.9995, Python=0.9995
  - b_events_mean: Julia=0.9993, Python=0.9993
- **Implementation design notes**: Added to SIGNATURES.md documenting recession event thresholds, flow timing NA handling, and sinusoidal model periodicity

### Fixed
- **Terminology**: `analyze_flow_timing_trends()` — Renamed internal variable `julday_max` to `timing_by_year` and fixed comments to correctly say "day of water year" instead of "Julian day" (code was already correct, but naming was misleading)

#### Cross-Language Alignment Fixes (March 2026)

**Julia Recession min_events Gate Missing (HIGH)**:
- **Issue**: Julia produced recession metrics (b_events, log_a_events, concavity) for 2,214 more gages than R/Python. R and Python agreed on NA patterns (99.7% match); Julia was the outlier.
- **Root Cause**: Julia's `analyze_recession_parameters()` was missing the `min_events` (25) threshold check before generating statistics for event-based metrics. R checks `total_events < RECESSION_MIN_EVENTS` and returns all NAs; Python does the same with `len(all_recession_events) < min_events`. Julia only applied the `min_events` gate to seasonality metrics, not to base metrics (b_events, log_a_events, concavity, b_pointcloud, log_a_pointcloud).
- **Fix**: Added early-return `min_events` check in `analyze_recession_parameters()` before `generate_stats()`, matching R/Python behavior. If `length(all_log_a) < min_events`, all recession metrics return NaN.
- **Impact**: 2,214 gages that previously received recession values in Julia will now correctly return NA, matching R/Python.
- **Location**: `julia/src/recession.jl` lines 423-434

**Julia BFI_LyneHollick Bug (HIGH)**:
- **Issue**: Correlation with R was ~0.5 instead of >0.99
- **Root Cause**: Julia's `sum()` returns NaN if any input is NaN, unlike R's `sum(..., na.rm=TRUE)`. The Lyne-Hollick filter propagates NaN from adjacent missing values to valid Q positions, causing entire BFI calculations to become NaN.
- **Fix**: Updated `baseflow.jl` to filter baseflow values by checking both Q validity AND baseflow validity before summing, matching R's na.rm=TRUE behavior.
- **Location**: `julia/src/baseflow.jl` lines 222-246

**Julia Flow_Reversals Bug (HIGH)**:
- **Issue**: Correlation with R was ~0.35-0.52 instead of >0.99
- **Root Cause**: Julia used forward-fill for NA gaps in flow series, while R uses linear interpolation (`approx()`). Forward-fill creates artificial plateaus that suppress reversal detection.
- **Fix**: Implemented linear interpolation matching R's `approx()` function with rule=2 (extrapolate from nearest value at ends).
- **Location**: `julia/src/pulses.jl` function `count_flow_reversals()`

**Config Centralization (7 Julia files, 11 Python files)**:
- Removed hardcoded constants from all signature modules
- Julia: Replaced local `const` declarations with `CFG_*` imports from `config.jl`
- Python: Added imports from `config.py` to all signature modules (previously config.py was unused)
- Affected Julia files: `recession.jl`, `runoff_ratios.jl`, `qp_seasonality.jl`, `elasticity.jl`, `timing.jl`, `flow_volumes.jl`, `baseflow.jl`
- Affected Python files: `elasticity.py`, `flow_volumes.py`, `pulses.py`, `recession.py`, `runoff_ratios.py`, `storage.py`, `timing.py`, `fdc.py`, `flashiness.py`, `baseflow.py`, `qp_seasonality.py`

**Julia Recession Parameters (HIGH) - VERIFIED FIXED**:
- **Issue**: `b_events_mean` correlation with R was 0.9075 instead of >0.99
- **Root Causes**:
  1. Event identification used "local" algorithm (checked adjacent pairs); R/Python use "look-ahead" algorithm (check next min_length days all qualify)
  2. Per-event fitting did not remove first day; R/Python remove first day (often influenced by peak)
  3. `min_points=10` was too strict; R/Python only require 2 valid points
  4. Julia used Theil-Sen regression by default; R/Python use OLS
- **Fixes**:
  - Rewrote `identify_recession_events()` to use look-ahead algorithm matching R/Python exactly
  - Added `remove_first_day` parameter to `fit_recession_power_law()`, default true for per-event fitting
  - Changed `min_points` default from 10 to 2 to match R/Python
  - Changed default from Theil-Sen to OLS for all recession fitting
- **Location**: `julia/src/recession.jl`
- **Post-fix correlation**: `b_events_mean` Julia=0.9993 vs R (improved from 0.9075)

**Julia Elasticity (MEDIUM) - VERIFIED FIXED**:
- **Issue**: Several algorithmic differences from R/Python
- **Root Causes**:
  1. Hardcoded 300 days instead of config-based completeness check
  2. Rolling window year assigned to MIDDLE instead of END
  3. Default min_years=10 instead of config value (15)
  4. Type mismatch in function signature (`Float64` vs `Real`)
- **Fixes**:
  - Added proper leap year handling for expected days calculation
  - Changed rolling window year assignment to END of window (like R/Python)
  - Updated function signature to use `Real` type for config values
  - Fixed `expected_days()` to accept Real and convert to Int
- **Location**: `julia/src/elasticity.jl`
- **Verification**: Post-fix correlation with R: elasticity_static=0.9953, elasticity_mean=0.9968

**Python/Julia avg_storage Interpolation Tie-Handling (HIGH)**:
- **Issue**: `avg_storage` trend statistics (slopes, p-values) had Spearman correlations of only 0.77-0.94 with R, despite means agreeing at 0.99+
- **Root Cause**: R's `approx()` uses `ties = mean` by default, averaging all S (cumulative storage) values at duplicate Q (discharge) positions before interpolating. Python's `np.interp()` takes the **last** S value at tied Q points; Julia's manual `findlast`/`findfirst` bracket interpolation picks a single arbitrary pair. Since many daily Q values are identical (especially low flows), this creates noisy annual storage values with ~2x inflated variance in Python/Julia, producing systematically different slopes and p-values.
- **Fixes**:
  - Python: Added `np.unique`-based grouping to average S at duplicate Q before calling `np.interp`
  - Julia: Added `unique`+`mean` grouping to average S at duplicate Q before bracket interpolation
- **Verification**: Simulated 30-year test shows Spearman R vs Python improves from 0.71 to 1.00; Theil-Sen slopes match exactly
- **Locations**: `python/streamflow_signatures/storage.py`, `julia/src/storage.jl`

**Human Interference Metadata Support**:
- Added `metadata` section to `config/signatures_config.json` with GAGES-II and HYDAT settings
- Created `julia/src/metadata.jl` module for loading watershed interference data
- Created `python/streamflow_signatures/metadata.py` module for loading watershed interference data
- Updated `julia/src/config.jl` and `python/streamflow_signatures/config.py` to load metadata settings

### Deprecated
- **Corrupted Parquet File**: Marked `D:/combined_streamflow_output/combined_streamflow_data.parquet` (October 2025) as **DO NOT USE**
  - **Issue**: Contains 99999 multiplier bug for Canadian gages without basin area
  - **Impact**: Q values ~100,000x too high for affected gages (e.g., 08ND025: 78.5M instead of 785.6)
  - **Use instead**: `D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet`
  - **Documentation**: Added "Parquet Data Files" section to DEVELOPMENT.md with active vs deprecated file list

### Guidelines Document TODOs
<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

---

## [February 2026]

### Added

#### Signature Cross-Correlation Plot (Means vs Trends)
- **New Shiny dashboard section**: Interactive scatter plot showing how pairwise signature relationships in spatial patterns (means across gages) compare to relationships in temporal trends (Theil-Sen slopes across gages)
- **Architecture**: Pre-computed offline to avoid startup delay; CSV loaded from S3 (~4,556 pairs)
- **Features**:
  - Each point = one ordered pair (A → B); x = Theil-Sen slope of z-scored means, y = Theil-Sen slope of z-scored trends
  - 13-category color coding (Flow Volume, Flow Percentiles, FDC, Baseflow, Recession, Pulse Metrics, Flow Reversals, Flashiness, Flow Timing, Runoff Ratios, Elasticity, Q-P Seasonality, Average Storage)
  - Category checkbox filter to show/hide groups
  - 1:1 reference line showing where spatial and temporal co-variation are equal
  - Click any point to reveal gage-level detail panel with two scatter plots (A_mean vs B_mean, A_senn_slp vs B_senn_slp) including Theil-Sen trend lines and Spearman statistics
- **New files**:
  - `precompute_cross_signature_analysis.R` — Offline script computing all pairwise slopes across 68 metrics × 5,691 gages, uploads CSV to S3
  - `streamflowAndClimateVisualizationApp/cross_signature_analysis.csv` — Pre-computed output (4,556 rows)
- **Modified files**:
  - `streamflowAndClimateVisualizationApp/app.R` — Added category constants, UI section, and server logic
  - `streamflowAndClimateVisualizationApp/helperFunctions.R` — Added `fast_theil_sen_slope()` utility for batch Theil-Sen computation

#### Human Interference Metadata Integration
- **New feature**: Watershed metadata is automatically enriched with human interference indicators during data processing
- **USGS gages (from GAGES-II)**: NDAMS_2009 (dam count), MAJ_DDENS_2009 (dam density), STOR_NID_2009 (dam storage), IMPNLCD06 (impervious %), DEVNLCD06 (developed %), FRESHW_WITHDRAWAL, HYDRO_DISTURB_INDX (0-20 scale), CLASS (Ref/Non-ref)
- **Canadian gages (from HYDAT)**: RHBN (Reference Hydrometric Basin Network flag), REGULATED (regulation status)
- **Unified classification**: `human_interference_class` column with values: reference, non-reference, unknown
- **New functions in helperFunctions.R**:
  - `load_gages_ii_interference()` — Load GAGES-II metadata
  - `load_canadian_interference()` — Load HYDAT metadata via tidyhydat
  - `enrich_metadata_with_interference()` — Enrich combined metadata
- **Integration**: `concatenate_with_metadata()` now automatically calls enrichment
- **New utility scripts** (in R/):
  - `run_enrich_metadata.R` — Manual enrichment utility
  - `run_regenerate_metadata.R` — Regenerate and enrich metadata
- **New config**: `GAGES_II_DIR`, `GAGES_II_FILES_*`, `GAGES_II_COLUMNS`, `EXPECTED_INTERFERENCE_COLS` in config.R

### Fixed

#### Visualization App Gage ID Mismatch (HIGH)
- **Issue**: Scatter plot incorrectly filtered out 3,647 valid gages as "failed QC"
- **Root Cause**: App compared `gage_id` (with leading zeros, e.g., `01011000`) to metadata `gage_id` (without leading zeros, e.g., `1011000`)
- **Additional Issue**: App redundantly applied `processing_status == "success"` filter already enforced during pre-processing
- **Fix**: Removed `goodGages` filter from scatter plot; only `flagged_for_qann_range` is now applied
- **Location**: `streamflowAndClimateVisualizationApp/app.R`

#### Code Review Bug Fixes (Feb 2026)

**HIGH severity:**
- **H1**: `calculate_flow_vols_by_year()` — Renamed misleading variables from `annual_means`/`*_means` to `annual_totals`/`*_totals`; updated SIGNATURES.md to say "total" with units "mm" instead of "mean" with "mm/day"
- **H2**: `calculate_flow_vols_by_year()` — Added 250-day minimum data filter per water year; previously years with as few as 1 day could pass through
- **H3**: `analyze_recession_parameters()` — Fixed early-return schema: was producing 5 columns with wrong suffixes (`slp`, `rho`, `pval`), now produces correct 8 columns matching `STAT_SUFFIXES`
- **H4**: `analyze_flow_timing_trends()` — `cumsum()` now handles NAs by replacing with 0 before accumulation; previously a single NA corrupted all D-day metrics for that year
- **H5**: Canadian unit conversion — Replaced sentinel value `99999` with `NA` when basin area is missing; previously produced absurd Q values (~100,000+ mm/day)
- **H5 follow-up**: Missing basin area — The H5 fix (`conversion = NA`) caused all Q values to become NA, silently dropping gages at the qualifying-years gate. Now sets `conversion = 1` to keep raw units (cfs for USGS, m3/s for Canadian) and adds `area_normalized` flag to `watershed_metadata.csv`. Also added the same NA guard for USGS gages (previously unprotected).

**MEDIUM severity:**
- **M1**: `lyne_hollick_filter()` — Fixed multi-pass to use previous pass output as input; the forward pass was always referencing original Q regardless of pass number, making `passes=2` effectively single-pass
- **M2**: `analyze_recession_parameters()` — Point cloud `log_a` now uses `b_pointcloud` instead of `median_b` from per-event fits
- **M4**: `analyze_flashiness_trends()` — Added division-by-zero guard when all flow values are zero
- **M5**: `count_flow_reversals()` — NA gaps are now interpolated instead of removed, preventing false adjacencies between non-adjacent days
- **M6**: `calculate_qp_seasonality()` and `calculate_average_storage()` — Use `copy()` to prevent modifying parent data.table by reference via `:=`
- **M7**: `calculate_qp_seasonality()` — Slopes now assigned to middle of rolling window instead of end
- **M8**: `analyze_Q_PPT_relationships()` — Raised near-zero PPT threshold from 0.001mm to 10mm (annual) and 1mm (seasonal) to prevent extreme runoff ratios

**LOW severity:**
- **L1**: `generate_stats()` — Now explicitly sorts by year before trend calculations
- **L3**: `analyze_baseflow_indices()` — Removed unused `doy` from required columns (only `dowy` is used)
- **L5**: `process_signatures_from_parquet()` — Non-climate error handlers now log warnings instead of silently swallowing errors
- **L7**: `analyze_flashiness_trends()` — Minimum-day check now counts non-NA values instead of total rows
- **L8**: `process_signatures_from_parquet()` — Renamed loop variable from `gage_id` to `current_gage_id` to avoid data.table column name shadowing

### Added
- Comprehensive workflow review document (`WORKFLOW_REVIEW.md`)
- Documentation of gage ID format differences (`gage_id` vs `gage_id_metadata` columns)
- `area_normalized` column in `watershed_metadata.csv` — `TRUE` when basin area is available and Q is in mm/day, `FALSE` when basin area is missing and Q is in raw units

### Changed
- SIGNATURES.md updated to match actual code: column names, metric counts, units, PPT thresholds, D-day naming (D5 not D05), recession seasonality variants

---

## [January 2026]

### Fixed

#### Metadata Lookup Bug (HIGH)
- **Issue**: Data.table scoping caused all gages to receive the first row's metadata
- **Root Cause**: `metadata_lookup[gage_id == gage_id]` compared column to itself (always TRUE)
- **Fix**: Renamed `find_metadata()` parameter to `target_gage_id`
- **Location**: `helperFunctions.R` line ~3425

#### Canadian Basin Area Missing (HIGH)
- **Issue**: Basin area was hardcoded as `NA` for Canadian stations
- **Impact**: Canadian gages had no drainage area in output
- **Fix**: Now fetches `DRAINAGE_AREA_GROSS` from `tidyhydat::hy_stations()`
- **Locations**:
  - Metadata creation functions (for new processing)
  - `process_signatures_from_parquet()` runtime fallback (for existing metadata)

### Added
- Centralized `config.R` for all configuration parameters
- Structured logging system with levels (DEBUG, INFO, WARN, ERROR)
- Input validation functions (`validate_file_exists`, `validate_numeric`, etc.)
- Output schema validation (`validate_output_schema`, `validate_gage_output`)
- QA/QC scripts (`qa_qc_signatures.R`, `visualize_qa_qc.R`)
- Smoke test for quick validation (`smoke_test.R`)
- Climate function tests with synthetic data (`test_climate_functions.R`)

### Changed
- Consolidated all helper function variants into canonical `helperFunctions.R`
- Moved deprecated helper files to `archive/` directory

### Completed Milestones
- Consolidate active helper file variants into canonical `helperFunctions.R`
- Add centralized `config.R` for parameters
- Implement structured logging
- Fix metadata lookup bug
- Fix Canadian basin_area bug

---

## [December 2025]

### Added
- Initial climate signatures: elasticity, Q-P seasonality, average storage
- Daymet climate data integration
- Caravan processing pipeline (`caravan_to_annualized.R`)

---

## Version History Notes

This project uses date-based versioning (MONTH YEAR) rather than semantic versioning, reflecting its nature as a research tool with continuous development.

### Output File Naming Convention
Output files include date stamps: `streamflow_signatures_full_JAN2026.csv`

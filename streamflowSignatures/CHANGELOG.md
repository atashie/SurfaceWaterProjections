# Changelog

All notable changes to the Streamflow Signatures project.
For historical entries (Dec 2025 – March 2026), see [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

## [Unreleased]

### Planned
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)
- Recession-informed BFI (Collischonn & Fan 2013) — deferred from Section 3 review
- Sync R/Python with Julia Section 3 changes (D1/D99, recession fix, n_recession_events, elasticity rename/annual, runoff ratio flag, diagnostics)

### Guidelines Section 3 Alignment — Julia-First (April 2026)

Six changes implementing Guidelines Document Section 3 function review. Julia-first; R/Python sync pending.

**Step 1 (LOW): Add D1 and D99 flow timing percentiles**
Added `1` and `99` to `d_percentiles` in `config/signatures_config.json`. Julia's `timing.jl` reads config dynamically — also fixed hardcoded DataFrame pre-allocation to use config-driven column creation. 16 new columns (`D1_day_*`, `D99_day_*`). Caveat: D1 may be near-constant for flashy streams; D99 may always be ~364-365.

**Step 2 (MEDIUM): Fix recession event detection algorithm**
Replaced look-ahead algorithm in `julia/src/recession.jl` `identify_recession_events()` with position-level marking. The old algorithm checked if the NEXT `min_length` days met criteria, which truncated events by up to `min_length` days from their true end. New algorithm marks each position where both criteria hold (Q decreasing AND |dQ/dt| decreasing), then finds contiguous runs >= `min_length`. Preserves inclusive-inclusive `(start, end)` tuple semantics. R/Python still have old algorithm — will diverge on recession metrics until synced.

**Step 3 (LOW): Add n_recession_events signature**
New `n_recession_events` column in `analyze_recession_parameters()` counting recession events per water year. Computed INDEPENDENTLY of the `min_events=25` gate — event counts are most informative for gages with few events. 8 new columns. Added to `EXPECTED_SIGNATURE_BASES` in `config.R`.

**Step 4 (LOW): Add runoff_ratio_high_count flag**
New per-gage scalar in `analyze_Q_PPT_relationships()` counting years where annual runoff ratio > 2.0. Registered as schema exception `EXPECTED_RUNOFF_RATIO_HIGH` in `config.R`.

**Step 5 (MEDIUM): Rename elasticity → elasticity_rolling + add elasticity_annual**
- Renamed rolling-window metric from `elasticity` → `elasticity_rolling` (8 cols renamed). Uses Sawicz et al. (2011) departure-from-mean method.
- Added `elasticity_annual`: year-over-year consecutive-difference method. Formula: `E_t = ((Q_t - Q_{t-1}) / (P_t - P_{t-1})) / (Q_mean / P_mean)`. 8 new columns. Breaking rename: downstream tools, R, Python need updating.

**Step 6 (LOW): Add elasticity data quality diagnostics**
New per-gage scalars `elasticity_years_total` and `elasticity_years_low_ppt` in `calculate_streamflow_elasticity()`. Registered as schema exception `EXPECTED_ELASTICITY_DIAGNOSTICS` in `config.R`.

**Net impact**: +32 new stat columns, 8 renamed, 3 scalar diagnostics = 594 total columns (was 551).

### Guidelines Alignment — Negative Q, Constant-SD, Ice Tracking (April 2026)

Three changes aligning code with user decisions on guidelines interpretation. All 3 codebases (R, Python, Julia) updated.

**Change 1 (MEDIUM): Constant-SD → flag only, never reject**
Removed constant-SD year rejection from `preprocess_daily_data()` in R, Python, and Julia. The per-month uniqueness check remains as a QA diagnostic in all 3 languages (`constant_sd_months`/`constant_sd_flag` in diagnostics). Years previously rejected for constant SD are now retained. This reverses NA Audit Fix 1 from earlier in April 2026 — user decided even 30+ consecutive identical days should only be flagged, not rejected.

**Change 2 (MEDIUM): Make negative-Q rejection config-driven + add `negative_ann` signature**
- Negative-Q year rejection is now conditional on `reject_negative_flow` in config (currently `false`). The `reject_negative` variable was previously dead code in all 3 languages — now properly wired into the rejection chain.
- Hardcoded fallback defaults changed from `true`/`TRUE` to `false`/`FALSE` in all 3 languages (R line 824, Python `config.py` line 135, Julia `config.jl` line 129) so behavior is consistent whether or not the config key is present.
- Added `calculate_negative_days()` function (R: `helperFunctions.R`, Python: `pulses.py`, Julia: `pulses.jl`): counts days with Q < 0 per water year, processed through `generate_stats()` for 8 statistics (`negative_ann_*`).
- Registered in signature orchestration (`process_signatures_from_parquet()` in R, `signatures.py`, `signatures.jl`).
- Added `"negative_ann"` to `EXPECTED_SIGNATURE_BASES` in `config.R`.

**Change 3 (LOW): Per-gage ice-affected day total in output**
- R: Added `ice_affected_days_total` column to `gage_row` in `process_signatures_from_parquet()`, summing `na_cause_ice` from preprocessor diagnostics.
- Python/Julia: Added same aggregation in benchmark runners (`run_python_benchmark.py`, `run_julia_benchmark.jl`).

**Change 4 (LOW): Julia season definitions from config**
Julia's `preprocess_daily_data()` previously hardcoded season month definitions. Now reads from `config/signatures_config.json` via new `CFG_NA_SEASON_MONTHS` constant in `julia/src/config.jl`, matching R and Python behavior.

**Change 5 (LOW): Julia `generate_stats()` short-data NaN fill**
Julia's `generate_stats()` returned an empty Dict when `nrow(df) < min_rows`, breaking the schema contract (callers expect all 8 stats per column). Now returns NaN for all 8 stats per column, matching R and Python.

### NA Handling Audit Fixes (April 2026)

Six fixes from Codex NA handling audit aligning code with `SIGNATURE_GUIDELINES.md` requirements. All 4 codebases (R, rpkg, Python, Julia) updated where applicable.

**Fix 1 (MEDIUM): Constant-SD year rejection**
Changed constant-SD from QA flag only to year rejection criterion, matching guidelines requirement that years with constant monthly standard deviations during non-zero flow be excluded. Added rejection logic after negative-flow check in `preprocess_daily_data()` across R, Python, and Julia.

**Fix 2 (LOW): NA-cause ice tracking**
Added per-year tracking of ice-affected NA days (from USGS qualifier flags) vs other NA causes. New `na_cause_ice` and `na_cause_other` fields in diagnostics output. Implemented in R (`R/helperFunctions.R`), Python (`python/streamflow_signatures/io.py`), and Julia (`julia/src/io.jl`).

**Fix 3 (LOW): seasonal_flags passthrough to runoff ratios**
Python and Julia `analyze_Q_PPT_relationships()` now accept and apply `seasonal_flags` to NaN-out incomplete-season runoff ratios, matching R canonical. Updated callers in `signatures.py`/`signatures.jl`.

**Fix 5 (LOW): Constant-SD detection filters Q > 0**
Python's constant-SD monthly check now filters to `Q > 0` before testing uniqueness, matching R/Julia behavior. Prevents zero-flow days from masking truly constant sensor readings.

**Fix 6 (LOW): Remove dead ad-hoc NA interpolation**
Removed legacy per-signature NA interpolation code from flashiness and pulse functions across Python and Julia. This code was dead — `preprocess_daily_data()` guarantees zero NAs in valid years — but was misleading. R's flashiness interpolation (4 lines) and pulse flow-reversal interpolation also cleaned up.

**Fix 8 (LOW): Guidelines seasonal completeness threshold**
Updated `docs/SIGNATURE_GUIDELINES.md` line 88 from "≥68%" to "≥80%" to match the actual config value in `config/signatures_config.json`.

### Mann-Kendall Tau-b Fix and n≤3 P-value Guard (April 2026)

Two fixes to statistical functions improve cross-language alignment from 8 poor columns to 4.

**Fix 1 (HIGH): Julia Mann-Kendall tau-b tie adjustment**
Julia's `mann_kendall_test()` in `julia/src/stats.jl` was computing tau-a (S / n_pairs) instead of tau-b (S / sqrt((n_pairs - T_y) * n_pairs)), where T_y is the number of tied pairs in y. R's `Kendall::MannKendall` and Python's `scipy.stats.kendalltau` both compute tau-b. This caused divergence in any column with tied values (pulse metrics, percentiles with zero-flow years).

**Fix 2 (LOW): n≤3 p-value guard in Python and Julia**
R's `Kendall::MannKendall` fails for n≤3 with IFAULT=12 and returns p=1.0. Python's `scipy.stats.kendalltau` returns correct exact p-values for small n (e.g., p=0.333 for n=3). Added `if n <= 3: p_value = 1.0` in both Python (`python/streamflow_signatures/stats.py`) and Julia (`julia/src/stats.jl`) to match R's behavior.

**Reverted: Constant-input MK fix** — Attempted changing (NaN, NaN) → (0.0, 1.0) for constant series in both Python and Julia. This caused severe regressions (n_low_pulses_all_mk_rho dropped from 0.98 to 0.47, Q1-Q30 mk_rho all regressed) because Julia/Python had constant annual values (all zeros) where R had tiny non-zero values from floating-point precision. Returning tau=0 instead of NaN created worse mismatches than NaN exclusion. Reverted after one benchmark cycle per Iron Rule.

**Investigation confirmed non-issues:**
- Recession OLS: Already identical across all 3 languages (R²=1.000000 for annual b/log_a values, max diff ≈ 5e-15)
- FDC90th 28-gage NA mismatch: R's `lm()` produces tiny non-zero slopes (≈1e-15) for near-zero regressions while Python's `linregress()` produces exactly 0.0. Documented as irreducible.

**Updated results (April 2026, post tau-b fix):**

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 |
|------|---------|-----------|--------|-------------|
| R vs Python | 0.9990 | 1.0000 | 0.7621 | 4 |
| R vs Julia | 0.9991 | 1.0000 | 0.7621 | 4 |
| Python vs Julia | 0.9999 | 1.0000 | 0.9979 | 0 |

542 perfect (R²>=0.999), 5 good (0.99-0.999), 4 poor (<0.99). 547 of 551 columns (99.3%) have R² >= 0.99. Remaining 4 poor columns are irreducible: 2 recession Spearman p-values (exact vs t-approximation for small n), 2 FDC90th p-values (28-gage NA mismatch from floating-point precision).

### R Canonical Alignment Fixes (April 2026)

Three fixes to `R/helperFunctions.R` align R canonical's `process_signatures_from_parquet()` with Python/Julia benchmark behavior. Driven by diagnostic investigation after re-running benchmarks showed regression (41 poor cols vs previous 20).

**Fix 1 (HIGH): Removed leftover min_Q_value_and_days filter from non-legacy path**
The non-legacy preprocessing path retained a min_Q filter loop (~10 lines at line 4908) that required 30+ days above minimum flow threshold per water year. This filter pre-dated `preprocess_daily_data()` and was never part of Python/Julia benchmarks. It caused two systematic divergences:
- **Climate year exclusion**: Low-flow years (e.g., 1984, 2011, 2018 for gage 08145000) were removed from `valid_climate_years` via `intersect(valid_climate_years, water_years_to_use)`, changing elasticity from 1.398 → 0.488 for that gage
- **Flow year exclusion**: R had fewer qualifying years than Python/Julia (e.g., 38 vs 46 for gage 11274500), changing FDCmid and other signatures
- Affected ~150 gages (2.5%) with measurable differences across all signature categories

**Fix 2 (LOW): Added negative Q filter in FDC calculation**
`analyze_fdc_trends_from_streamflow()` line 1525: changed `Q[!is.na(Q)]` to `Q[!is.na(Q) & Q >= 0]`, matching Python/Julia. Safety net — preprocessor already rejects years with negative Q, but 81 gages have some negative Q rows.

**Fix 3 (MEDIUM code quality, zero current impact): Pass seasonal_flags to signature functions**
`calculate_flow_vols_by_year()` and `analyze_Q_PPT_relationships()` now receive `seasonal_flags` from the preprocessor, matching Python/Julia. No current impact (all gages have complete seasons), but ensures correctness if incomplete-season gages appear.

**Diagnostic investigation confirmed non-issues:**
- `daymet_id_for_gage` mapping: correct (Daymet parquet uses "08145000" directly)
- Climate merge on Date: produces identical PPT values to Python
- `generate_stats()` parity: all stats match except p-value library differences (expected)
- Preprocessor output: identical between benchmark path and direct path

**Cross-language audit (Codex)** found 3 remaining low-impact divergences:
- Julia hardcodes season month definitions instead of reading from JSON
- Julia relies on residual-NA rejection for PPT max-gap instead of explicit check
- Config APIs differ structurally (flat keys vs nested JSON) but numeric values match

### Julia leftjoin Row-Order Bug Fix (April 2026)
**HIGH**: `preprocess_daily_data()` in `julia/src/io.jl` produced incorrect `valid_climate_years` because `DataFrames.jl`'s `leftjoin` does not preserve row order. After joining the daily grid with year data, rows were reordered so that missing-PPT dates (e.g., Dec 31 at dowy=92) appeared at boundary positions (position 365 = Sep 30). The interpolation logic uses positional indices, so 1-day internal gaps were misidentified as boundary gaps and could not be interpolated — causing those years to be excluded from `valid_climate_years`.

**Fix**: Added `sort!(joined, :dowy)` after the leftjoin (1 line). For gage 01011000, `valid_climate_years` went from 32 → 43 (matching Python). Across the full benchmark, 0 entirely-NA climate columns (was 48 before fix).

### Interactive Comparison Dashboard (April 2026)
New `docs/benchmarks/build_comparison_dashboard.py` generates an interactive HTML dashboard (`signature_comparison.html`) comparing Julia benchmark output against golden R reference. Features: two synced Leaflet maps (original filter vs new filter), Plotly scatterplot with identity line, click-to-zoom from scatterplot to maps, R² summary. Covers 159 signature columns across 5,687 common gages.

### Remove Per-Signature Min-Days Thresholds (April 2026)
Removed all per-year min_days/max_na_frac/min_data_completeness checks from signature functions across all 4 codebases. The preprocessor (`preprocess_daily_data()`) is now the single source of truth for year qualification.

**Root cause**: Two independent year-qualification systems were fighting each other — per-signature min_days thresholds rejected preprocessor-approved years, creating artificial NaN years that then failed the 80% trend completeness check, causing massive NA inflation (e.g., elasticity dropped from 5,683 to 19 gages).

**Removed from `config/signatures_config.json`** (13 per-year thresholds):
- `baseflow.min_days`, `baseflow.max_missing_frac`
- `fdc.min_days`, `flashiness.min_days`, `flashiness.max_missing_frac`
- `flow_volumes.min_days`, `timing.min_days`, `pulses.min_days`
- `elasticity.min_data_completeness`
- `qp_seasonality.min_days_per_year`, `qp_seasonality.max_na_frac`
- `storage.min_days_per_year`, `storage.max_na_frac`

**Retained** (gage-level mathematical minimums, not per-year data quality):
- `recession.min_length=5`, `recession.min_events=25`
- `elasticity.min_years=15`, `elasticity.min_annual_ppt=10`
- `qp_seasonality.min_years=10`, `storage.min_years=10`

**Set `use_legacy_filtering: false`** permanently.

**Exempted recession and elasticity from trend_completeness** across all 4 codebases (recession is event-based/inherently sparse, elasticity uses rolling windows producing fewer values than input years).

**Updated benchmark runners** (Python, rpkg) to use `preprocess_daily_data()` path with trend_completeness/decade_completeness, matching Julia benchmark runner.

### Julia Adversarial Review (April 2026, Codex)
6 HIGH, 1 MEDIUM, 1 LOW findings. No critical bugs. Fixes deferred.

- **HIGH** FDC/flashiness min-days threshold 250 vs R's 30 (`fdc.jl`, `flashiness.jl`, `config.jl`) — silently drops years R keeps
- **HIGH** `Q95_Q10` column name vs R's `Q95-Q10` (`flow_volumes.jl`) — schema mismatch
- **HIGH** `generate_stats()` returns empty Dict on short data instead of NaN row (`stats.jl`) — breaks schema
- **HIGH** Runoff ratios aggregate Q/PPT jointly not independently (`runoff_ratios.jl`) — systematic bias when missingness differs
- **HIGH** Mann-Kendall tau denominator not tie-adjusted (`stats.jl`) — likely explains some of the 20 poor benchmark columns
- **MEDIUM** All-NA season yields NaN vs R's 0 (`flow_volumes.jl`)
- **LOW** Arrow IPC output mislabeled as parquet (`io.jl`)

### Standardized NA Handling (April 2026)
Centralized pre-processing for missing data, implemented identically across R canonical, rpkg, Python, and Julia. Driven by Guidelines Document "NA Handling" section.

**New `preprocess_daily_data()` function** (all 4 codebases):
- Daily grid normalization: one row per day per water year, sorted, unique dates
- Linear interpolation of internal gaps <= 3 consecutive days (no extrapolation)
- Year rejection: >30 raw NAs, gaps >3 days, negative Q values
- Residual boundary NAs (leading/trailing) cause year rejection
- Constant-SD QA flag: detects sensor flatline (flag only, not rejection)
- Seasonal completeness flags from RAW observations (pre-interpolation)
- Separate climate (PPT) NA policy with independent tracking

**Config changes** (`config/signatures_config.json`):
- New `na_handling` section with interpolation, year_rejection, constant_sd_flag, trend_completeness, seasonal_completeness, climate_na_policy
- `use_legacy_filtering: true` flag for backward-compatible staged migration

**fillna(0) removal** (3 functions, all 4 codebases):
- `analyze_flow_timing_trends`: no longer replaces NA Q with 0 before cumsum
- `calculate_average_storage`: skips years with residual NAs instead of zero-filling
- `calculate_qp_seasonality`: skips years with residual NAs instead of zero-filling

**Trend completeness** in `generate_stats()`:
- New `trend_completeness` and `decade_completeness` parameters
- Requires >=80% non-NA values overall and in first/last decade for trend stats
- Series <10 years: decade check skipped
- If incomplete: 6 trend stats set to NA, mean/median still computed

**Seasonal completeness** in flow volumes and runoff ratios:
- `calculate_flow_vols_by_year()` and `analyze_Q_PPT_relationships()` accept `seasonal_flags`
- Incomplete seasons (RAW observation fraction <80%) set to NA per year

**Tests**: `R/tests/test_na_handling.R` with 30+ test cases covering all edge cases

### Guidelines Document TODOs
<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

---

## [March 2026]

### rpkg: Proper R Package
New R package (`rpkg/`) rewriting the monolithic `R/helperFunctions.R` as a proper installable package. 16 R source files, 6 test files, 27 man pages. R CMD check: 0 errors, 0 warnings. See [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md) for full details.

### Cross-Language Alignment (Rounds 1–6, complete)
Six rounds of alignment work brought R, Python, and Julia implementations to production-ready agreement. Final status (March 16-17, 2026):

| Pair | Mean R² | Median R² | Cols < 0.99 |
|------|---------|-----------|-------------|
| R vs Python | 0.9968 | 1.0000 | 17 |
| R vs Julia | 0.9965 | 1.0000 | 19 |
| Python vs Julia | 0.9997 | 1.0000 | 3 |

- 475 perfect (R²>=0.999), 56 good (0.99-0.999), 20 poor (<0.99)
- 531 of 551 columns (96.4%) have R² >= 0.99 across all 3 pairs
- Golden output regression: 551/551 perfect match against Feb 2026 golden outputs
- 20 remaining poor columns are all trend statistics (slopes, p-values) — underlying means/medians are near-identical
- Full round-by-round details: [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md), [docs/CROSS_LANGUAGE_STATUS.md](docs/CROSS_LANGUAGE_STATUS.md)

### Code Quality Improvements (Rounds 1–2)
Cross-language code review covering ~9,000+ lines. Key fixes: O(n²) metadata lookup → O(1), Python vectorized water year calc (~50x speedup), Julia DataFrame pre-allocation, Eckhardt BFI forward-fill aligned across all 3 languages. Full details in `docs/CODE_REVIEW.md` and [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

### Benchmark Timing (March 16-17, 2026)

| Implementation | Time | Gages | Rate |
|---------------|------|-------|------|
| Julia | 9.2 min | 7,369 | 13.4/s |
| Python | 78.9 min | 7,369 | 1.56/s |
| rpkg | 114 min | 7,369 | 1.08/s |
| R (canonical) | 874 min* | 5,707 | 0.11/s |

*R ran concurrently with Python/Julia — timing inflated by I/O contention.

---

## Version History Notes

This project uses date-based versioning (MONTH YEAR) rather than semantic versioning, reflecting its nature as a research tool with continuous development.

### Output File Naming Convention
Output files include date stamps: `streamflow_signatures_full_JAN2026.csv`

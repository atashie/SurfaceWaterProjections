# Changelog Archive

Historical changelog entries for the Streamflow Signatures project.
For current/active changes, see [CHANGELOG.md](../CHANGELOG.md).

---

## [April 2026]

> Detailed record of April 2026 work. For the at-a-glance summary, see [CHANGELOG.md](../CHANGELOG.md). File-level change lists are omitted here — see `git log` for exact files.

### Changepoint Detection — Pettitt Test

Non-parametric Pettitt changepoint test applied to all time-series signatures, with differential metrics.

**Pettitt Test** (`pettitt_test`). Non-parametric rank-based test (Pettitt 1979). Test statistic K = max|U_t| where U_t accumulates pairwise sign comparisons. P-value via asymptotic approximation. Always identifies a location; use p-value for significance.

**Differential Metrics** (`segment_differential_metrics`), applied at the Pettitt changepoint: delta_mean (post − pre), pct_change (%), pre/post Mann-Kendall p-values.

Search window: WY 1980-2024, minimum 20 non-NA observations total, minimum 10 per segment. Independent of the 80% trend_completeness gate.

**8 new columns per signature (76 signatures × 8 = 608 changepoint columns):** `pettitt_cp_year`, `pettitt_pval`, `pettitt_pre_mean`, `pettitt_post_mean`, `pettitt_delta_mean`, `pettitt_pct_change`, `pettitt_pre_mk_pval`, `pettitt_post_mk_pval`.

**Design decisions:**
- Per-signature changepoints (not per-gage) — different processes change at different times
- Pettitt O(n²) is fast for annual series (n ≤ 45)
- BIC 4-model comparison was initially implemented alongside Pettitt but removed to reduce clutter (code retained in `changepoint.jl` but not called from `generate_stats()`)

**Julia benchmark (April 28, 2026)**: 7,313 gages, 1,264 columns (656 base/metadata + 608 changepoint). Overall ~13% of evaluations have p < 0.05 (69,279 / 516,992; 2.7× the 5% null). CP year range [1989, 2014].

**Signal robustness (April 29)**: significance rates vary by category (Flow Timing 3.7% below null, Flashiness 19.4%, Elasticity 25.8% inflated by short rolling-window series); effective independence ~17/76 signatures (~77% redundancy); ~3.5% survive BH-FDR at 0.05. Full analysis in [docs/SIGNATURES.md](SIGNATURES.md) → Changepoint Detection → Signal Robustness & Known Limitations.

**Reference:** Pettitt, A.N. (1979). A non-parametric approach to the change-point problem. Applied Statistics, 28(2), 126-135.

### Recession-Parameterized Baseflow Signatures

Three new signatures using recession-derived filter parameters (Collischonn & Fan 2013, Eckhardt 2005). Julia canonical April 24; Python and rpkg ported April 25.

- **`alpha_linear`** (8 stats): discrete recession constant under linear reservoir (b=1); per-year median of Q_{i+1}/Q_i ratios from point-cloud recession data (first day of each event removed). Per-year threshold >10 valid alpha pairs.
- **`recession_alpha_point_cloud_linear_reservoir`** (per-gage scalar): whole-record median of Q_{i+1}/Q_i ratios; computed independently of the 25-event gate. Parameterizes the two BFI signatures below.
- **`BFI_Eckhardt_param`** (8 stats): Eckhardt BFI with recession-derived `a`; BFImax fixed at 0.8.
- **`BFI_LyneHollick_param`** (8 stats): Lyne-Hollick BFI with recession-derived `alpha` (heuristic — L-H alpha has no physical derivation from recession; Nathan & McMahon 1990).

**Total: 24 new stat columns + 1 per-gage scalar = 25 new columns (624 total signature columns).**

**Julia validation (April 24)**: 7,313 gages, 15.0 min. 598 of 599 common columns Perfect (R² ≥ 0.999). Zero regression.

| Metric | Non-NA gages | Range |
|--------|-------------|-------|
| `alpha_linear_median` | 5,825 / 7,313 | [0.273, 0.985] |
| `recession_alpha_point_cloud_linear_reservoir` | 7,131 / 7,313 | [0.273, 0.997] |
| `BFI_Eckhardt_param_mean` | 7,131 / 7,313 | [0.474, 0.804] |
| `BFI_LyneHollick_param_mean` | 7,131 / 7,313 | [0.255, 0.966] |

**Cross-language validation (April 27)**: Python — all 25 new columns Perfect (615 Perfect / 5 Good / 3 Poor overall vs Julia). rpkg — 24 of 25 new Perfect; 1 new poor `alpha_linear_spearman_pval` (R² 0.971, Spearman p-value library difference). Current synced-implementation benchmark tables: [docs/DEVELOPMENT.md](DEVELOPMENT.md) → Cross-Language Benchmarks.

**Known limitation — BFI_Eckhardt_param narrow range** ([0.474, 0.804], std 0.036 vs fixed-parameter [0.139, 0.802], std 0.119): a mathematical property of the Eckhardt filter when BFImax is fixed at 0.8 — see [docs/SIGNATURES.md](SIGNATURES.md) → Implementation Design Notes. Future improvement: per-gage BFImax via Collischonn & Fan (2013) backward filter.

### Canonical Language Transition

Julia replaced R as the canonical implementation. All three synced implementations (Julia, Python, rpkg) reached equilibrium (99.5% R² agreement, 594 signature columns at the time, 7,313 gages). Julia is ~40× faster than monolithic R (~27 min vs ~18 hrs) with cleaner modular architecture, and already served as first-mover for Section 3 changes.

- Julia (`julia/src/`) is now the canonical reference for all signature implementations
- rpkg (`rpkg/`) is the production-ready R port (replaces monolithic R for new development)
- Monolithic R (`R/helperFunctions.R`) deprecated as legacy shim — retained for data ingestion utilities and existing callers; no new signature features
- Data ingestion stays in R (dataRetrieval/tidyhydat)
- Recession algorithm sync to rpkg only (monolithic R skipped — deprecated)

Documentation updated across `CLAUDE.md`, `README.md`, `docs/`, package READMEs, and skills; deprecation headers added to the legacy R entry points.

### Julia Experiment Infrastructure

Sensitivity experiment framework for Julia benchmarks. Three experiments test how restricting the analysis period and gage quality affect results.

- **Config** (`config/signatures_config.json`): `filtering.start_water_year` and `filtering.min_qualifying_data_fraction` (both null by default).
- **Infrastructure**: ENV-based config overrides in `julia/src/config.jl` (read at runtime to avoid precompilation caching of `const` values), filtering logic + parameterized output paths in the benchmark runner, three thin experiment wrappers, and parameterized comparison/dashboard scripts (`docs/benchmarks/`).

**Results (April 16-17)**:

| Experiment | Gages | Dropped | Time | Rate |
|------------|-------|---------|------|------|
| Baseline | 7,313 | — | 27.5 min | 4.4/s |
| startIn1993 | 6,678 | 635 | 24 min | 4.7/s |
| startIn1993+60pct | 6,579 | 734 | 18 min | 5.9/s |
| startIn1993+80pct | 5,431 | 1,882 | 14 min | 6.5/s |

- 635 gages dropped by WY≥1993 (485 USGS, 150 Canadian — lost pre-1993 years below the 20-year min)
- 99 more dropped by the 60% filter (all had exactly 20 qualifying years; frac=20/34=0.588 due to partial WY2026 in raw parquet)
- 1,148 more dropped by the 80% filter (need ≥28/34 valid years)
- Trend statistics diverge vs baseline (expected — different period-of-record); means/medians stable

**Bug found and fixed**: Julia precompilation caches `const` values in `.ji` files, so ENV vars set by experiment wrappers were silently ignored after first compilation. Fixed by reading ENV vars at runtime in `main()` as local variables.

### Guidelines vs Implementation Audit

Systematic audit comparing `SIGNATURE_GUIDELINES.md` against Python, Julia, and rpkg.

- **Season exclusion year counts (HIGH)**: added 4 per-gage diagnostics `season_excluded_years_{winter,spring,summer,fall}` (years where each season failed the 80% completeness threshold). Registered as `EXPECTED_SEASON_EXCLUDED_YEARS` in `config.R`.
- **Q95_Q10 documentation fix (MEDIUM)**: corrected "ratio" → "difference" (the metric is Q95 − Q10) in code comments and docs.
- **Trend completeness exemptions (LOW)**: documented in `SIGNATURES.md` why recession and elasticity are exempt from the 80% gate.
- **rpkg dead NA interpolation (2× HIGH)**: removed inline `approx(rule=2)` from `rpkg/R/flashiness.R` and `rpkg/R/pulses.R` (missed in the April cleanup that hit Python/Julia/R-canonical); replaced with `stopifnot` assertions.
- **QA/QC flag integration (MEDIUM)**: added optional `include_qa_flags` to `calculate_all_signatures()` in all 3 packages (appends 12 flag columns from `compute_qa_flags()`).
- **Recession-informed BFI (MEDIUM)**: documented as not-yet-implemented (later implemented — see Recession-Parameterized Baseflow above).
- **Elasticity 30% diagnostic (MEDIUM)**: deferred pending domain-expert clarification.
- **14 glossary/guidelines items** flagged for the user to update in the Google Doc (Qann "mean" vs "total", elasticity rename, PPT thresholds, negative-Q and constant-SD contradictions, stale "julian" note, etc.) — no code changes.

### Section 3 Sync to Python & rpkg + Benchmark Re-Run

Synced all 7 Section 3 changes from Julia to Python and rpkg, re-ran all 3 benchmarks, and validated.

- **Julia runoff_ratios fix (LOW)**: `>=` → `>` for the PPT threshold, matching R/Python.
- **Ported to Python** (`recession.py`, `elasticity.py`, `runoff_ratios.py`, `signatures.py`) and **rpkg** (`recession.R`, `elasticity.R`, `runoff_ratios.R`, `pulses.R`, `signatures.R`, `NAMESPACE`): position-marking recession algorithm, `n_recession_events`, elasticity rename + `elasticity_annual` + diagnostics, `runoff_ratio_high_count`, `calculate_negative_days`. Synced rpkg's bundled config (was severely outdated).
- Updated comparison-script categorization and fixed Julia benchmark-runner metadata column typing.

**Benchmark results (April 14-15)** — all 3 synced implementations at 594 signature columns, 7,313 gages:

| Pair | Perfect (≥0.999) | Good (0.99-0.999) | Poor (<0.99) | Min R² |
|------|-------------------|-------------------|-------------|--------|
| Python vs Julia | 619 | 5 | 3 | 0.986 |
| rpkg vs Julia | 586 | 4 | 4 | 0.967 |
| rpkg vs Python | 582 | 5 | 7 | 0.967 |

Poor columns are irreducible library-level differences (recession pointcloud p-values, FDC90th p-values, recession seasonality minimum). rpkg vs R canonical: 490 perfect, 46 poor (all recession — R canonical still has the old look-ahead algorithm).

**Net column impact**: 594 columns (was 551 pre-Section 3): +16 D1/D99, +8 n_recession_events, +8 elasticity_annual, +8 negative_ann, +2 elasticity diagnostics, +1 runoff_ratio_high_count; elasticity renamed (8 cols).

### Julia vs Golden R Comparison & Fixes

- **Comparison tooling**: `compare_julia_vs_golden_r.py` (6-tier R² across 551 common columns, 5,697 common gages), `build_julia_vs_golden_r_dashboard.py`, and `julia_vs_golden_r_summary.md`.
- **Elasticity operator fix (LOW)**: `>=` → `>` for `min_annual_ppt` in `julia/src/elasticity.jl`, aligning all 3 languages.
- **Canadian HYDAT metadata for Julia**: previously all 1,333 Canadian gages were `human_interference_class = "unknown"`. New `R/export_hydat_metadata.R` exports 8,012 stations (RHBN + REGULATED) to `metadata/canadian_hydat_interference.csv` (tracked in git); Julia now reads it and merges RHBN/REGULATED/class for Canadian gages.
- **Divergence analysis** (288 columns with R² < 0.99 vs Golden R): recession (46 cols, intentional algorithm change), trend_completeness (220 cols, Golden R predates the 80% gate), elasticity (9 cols, operator bug + since-removed min_Q filter), dur_low_pulses_all (6 cols, investigation pending).

### Guidelines Section 3 Alignment — Julia-First

Six changes implementing Guidelines Document Section 3. Julia-first; Python and rpkg synced April 14-15.

1. **(LOW) D1/D99 flow timing percentiles** — added `1` and `99` to `d_percentiles`; made Julia's `timing.jl` allocation config-driven. 16 new columns. (D1 may be near-constant for flashy streams; D99 ~364-365.)
2. **(MEDIUM) Recession event-detection algorithm** — replaced the look-ahead algorithm in `identify_recession_events()` with position-level marking (old algorithm truncated events by up to `min_length` days). Marks each position where both Q and |dQ/dt| decrease, then finds contiguous runs ≥ `min_length`.
3. **(LOW) `n_recession_events`** — per-water-year event count, computed independently of the `min_events=25` gate. 8 new columns.
4. **(LOW) `runoff_ratio_high_count`** — per-gage count of years with annual runoff ratio > 2.0.
5. **(MEDIUM) Elasticity rename + `elasticity_annual`** — renamed `elasticity` → `elasticity_rolling` (Sawicz 2011 departure-from-mean); added year-over-year `elasticity_annual` = ((Q_t − Q_{t-1})/(P_t − P_{t-1}))/(Q_mean/P_mean). 8 renamed + 8 new.
6. **(LOW) Elasticity diagnostics** — per-gage `elasticity_years_total` and `elasticity_years_low_ppt`.

**Net impact**: +32 new stat columns, 8 renamed, 3 scalar diagnostics = 594 total (was 551).

### Guidelines Alignment — Negative Q, Constant-SD, Ice Tracking

Three changes aligning code with user decisions on guidelines interpretation. All 3 codebases (R, Python, Julia).

1. **(MEDIUM) Constant-SD → flag only, never reject** — removed constant-SD year rejection from `preprocess_daily_data()`; the per-month uniqueness check remains as a QA diagnostic. Reverses NA Audit Fix 1 (below) — user decided even 30+ consecutive identical days should only be flagged.
2. **(MEDIUM) Config-driven negative-Q rejection + `negative_ann`** — negative-Q year rejection is now conditional on `reject_negative_flow` (default false; previously dead code). Added `calculate_negative_days()` (counts Q<0 days per water year → 8 stats). Registered `negative_ann` in `EXPECTED_SIGNATURE_BASES`.
3. **(LOW) Per-gage `ice_affected_days_total`** in output, summing `na_cause_ice` from preprocessor diagnostics.
4. **(LOW) Julia season definitions from config** via new `CFG_NA_SEASON_MONTHS`.
5. **(LOW) Julia `generate_stats()` short-data NaN fill** — returns NaN for all 8 stats instead of an empty Dict, restoring the schema contract.

### NA Handling Audit Fixes

Six fixes from a Codex NA-handling audit aligning code with `SIGNATURE_GUIDELINES.md`. (Fixes 1 and 5 partly superseded by the Negative-Q/Constant-SD change above.)

1. **(MEDIUM) Constant-SD year rejection** — added (later reverted to flag-only, above).
2. **(LOW) NA-cause ice tracking** — new `na_cause_ice` / `na_cause_other` diagnostics from USGS qualifier flags.
3. **(LOW) seasonal_flags → runoff ratios** — Python/Julia now NaN-out incomplete-season runoff ratios, matching R.
4. **(LOW) Constant-SD detection filters Q > 0** in Python, matching R/Julia.
5. **(LOW) Removed dead ad-hoc NA interpolation** from flashiness/pulse functions (Python/Julia; R cleaned too) — `preprocess_daily_data()` already guarantees zero NAs in valid years.
6. **(LOW) Guidelines seasonal-completeness threshold** corrected "≥68%" → "≥80%" to match config.

### Mann-Kendall Tau-b Fix and n≤3 P-value Guard

Two fixes improving cross-language alignment from 8 poor columns to 4.

1. **(HIGH) Julia Mann-Kendall tau-b** — was computing tau-a (S / n_pairs) instead of tau-b (S / sqrt((n_pairs − T_y)·n_pairs)). R and Python both use tau-b; this caused divergence on any column with ties (pulse metrics, percentiles with zero-flow years).
2. **(LOW) n≤3 p-value guard** — added `p_value = 1.0` for n≤3 in Python and Julia to match R's `Kendall::MannKendall` (which fails with IFAULT=12 → p=1.0).

**Reverted: constant-input MK fix** — changing (NaN, NaN) → (0.0, 1.0) for constant series caused severe regressions (Julia/Python had all-zero annual values where R had tiny floating-point non-zeros), so returning tau=0 created worse mismatches than NaN exclusion. Reverted after one benchmark cycle.

**Confirmed non-issues**: recession OLS already identical across languages (R²=1.000000); FDC90th 28-gage NA mismatch is irreducible (R's `lm()` yields ≈1e-15 slopes where Python's `linregress()` yields exactly 0.0).

**Results (post tau-b fix)**: R vs Python / R vs Julia both 4 cols < 0.99; Python vs Julia 0. 547 of 551 columns (99.3%) at R² ≥ 0.99. Remaining 4 poor: 2 recession Spearman p-values + 2 FDC90th p-values (both irreducible).

### R Canonical Alignment Fixes

Three fixes to `R/helperFunctions.R` aligning `process_signatures_from_parquet()` with Python/Julia benchmark behavior.

1. **(HIGH) Removed leftover min_Q_value_and_days filter** from the non-legacy path (~10 lines requiring 30+ days above a minimum flow per year). It excluded low-flow years from both `valid_climate_years` (e.g., elasticity 1.398 → 0.488 for gage 08145000) and flow years (e.g., 38 vs 46 years for 11274500). Affected ~150 gages (2.5%).
2. **(LOW) Negative-Q filter in FDC** — `Q[!is.na(Q)]` → `Q[!is.na(Q) & Q >= 0]`, matching Python/Julia (safety net; 81 gages have some negative rows).
3. **(MEDIUM, zero current impact) Pass seasonal_flags** to flow-volume and runoff-ratio functions.

A Codex cross-language audit also found 3 remaining low-impact divergences (Julia hardcoded season months — later fixed; Julia PPT max-gap via residual-NA rejection; flat-key vs nested-JSON config APIs with matching numeric values).

### Julia leftjoin Row-Order Bug Fix

**HIGH**: `preprocess_daily_data()` produced incorrect `valid_climate_years` because `DataFrames.jl`'s `leftjoin` doesn't preserve row order — missing-PPT dates landed at boundary positions, so 1-day internal gaps were misread as boundary gaps and not interpolated. Fixed with `sort!(joined, :dowy)` after the join. For gage 01011000, `valid_climate_years` went 32 → 43; across the benchmark, 0 entirely-NA climate columns (was 48).

### Interactive Comparison Dashboard

New `docs/benchmarks/build_comparison_dashboard.py` → interactive HTML (`signature_comparison.html`): two synced Leaflet maps (original vs new filter), Plotly identity-line scatterplot, click-to-zoom, R² summary. 159 signature columns, 5,687 common gages.

### Remove Per-Signature Min-Days Thresholds

Removed all per-year min_days/max_na_frac/min_data_completeness checks from signature functions across all 4 codebases — `preprocess_daily_data()` is now the single source of truth for year qualification.

**Root cause**: two competing year-qualification systems — per-signature thresholds rejected preprocessor-approved years, creating artificial NaN years that then failed the 80% trend-completeness check, massively inflating NAs (elasticity dropped from 5,683 to 19 gages).

- **Removed** 13 per-year thresholds from config (baseflow, fdc, flashiness, flow_volumes, timing, pulses, elasticity, qp_seasonality, storage).
- **Retained** gage-level mathematical minimums: `recession.min_length=5`, `recession.min_events=25`, `elasticity.min_years=15`, `elasticity.min_annual_ppt=10`, `qp_seasonality.min_years=10`, `storage.min_years=10`.
- Set `use_legacy_filtering: false` permanently; exempted recession and elasticity from trend_completeness; updated Python/rpkg benchmark runners to the preprocessor path.

### Julia Adversarial Review (Codex)

6 HIGH, 1 MEDIUM, 1 LOW findings; no critical bugs. (All since addressed by the fixes above.)

- **HIGH** FDC/flashiness min-days threshold 250 vs R's 30 — silently drops years
- **HIGH** `Q95_Q10` vs R's `Q95-Q10` column name — schema mismatch
- **HIGH** `generate_stats()` returns empty Dict on short data instead of a NaN row
- **HIGH** runoff ratios aggregate Q/PPT jointly not independently
- **HIGH** Mann-Kendall tau denominator not tie-adjusted
- **MEDIUM** all-NA season yields NaN vs R's 0
- **LOW** Arrow IPC output mislabeled as parquet

### Standardized NA Handling

Centralized pre-processing for missing data, implemented identically across R canonical, rpkg, Python, and Julia (driven by the Guidelines "NA Handling" section). This is the foundational change underlying most of the April alignment work above. Architecture diagram: [docs/DEVELOPMENT.md](DEVELOPMENT.md) → NA Handling Architecture.

**New `preprocess_daily_data()`** (called once per gage, before all signatures): daily-grid normalization (one row/day, sorted, unique); linear interpolation of internal gaps ≤3 days (no extrapolation); year rejection (>30 raw NAs, gaps >3 days, negative Q); residual boundary NAs reject the year; constant-SD QA flag (sensor flatline); seasonal completeness flags from raw observations; separate climate (PPT) NA policy.

**Config**: new `na_handling` section; `use_legacy_filtering` flag (initially true for staged migration, later set false).

**fillna(0) removal** (3 functions): `analyze_flow_timing_trends`, `calculate_average_storage`, `calculate_qp_seasonality` no longer zero-fill — they skip years with residual NAs.

**Trend completeness** in `generate_stats()`: new `trend_completeness`/`decade_completeness` params requiring ≥80% non-NA overall and in the first/last decade; if incomplete, 6 trend stats → NA, mean/median still computed (decade check skipped for series <10 years).

**Seasonal completeness**: `calculate_flow_vols_by_year()` and `analyze_Q_PPT_relationships()` accept `seasonal_flags`; incomplete seasons (raw observation fraction <80%) → NA per year.

**Tests**: `R/tests/test_na_handling.R`, 30+ cases.

---

## [March 2026] — Completed Items

### rpkg: Proper R Package Implementation (March 2026)

**New R package** (`rpkg/`) rewriting the monolithic `R/helperFunctions.R` as a proper installable R package (`streamflowsignatures`). Mirrors the Julia/Python package structure with 18 source files, shared `config/signatures_config.json`, and full test suite.

**Package structure**: DESCRIPTION, NAMESPACE (roxygen2), 16 R source files, 6 test files, 27 man pages. R CMD check: 0 errors, 0 warnings.

**Benchmark Results (March 2026)**:
- 7,369 gages processed, 551 signature columns, 0 errors, 114 min (1.08 gages/s)
- rpkg vs canonical R: 490/551 perfect (>=0.999), 51 good, 10 poor
- rpkg vs Python: 527/551 perfect, 18 good, 6 poor
- rpkg vs Julia: 515/551 perfect, 28 good, 8 poor

rpkg aligns more closely with Python/Julia than canonical R does (527 perfect vs Python, compared to 475 for canonical R), because rpkg uses config-driven parameters and code review fixes that match the Python/Julia implementations.

**Code review fixes applied during rpkg development**:
- `filter_qualifying_years()`: Stage 2 counts non-NA Q values (not total rows), matching Python/Julia
- `fdc.R`: Uses config `fdc_min_days=250` instead of hardcoded 30; adds negative Q filtering before `log10()`
- `flashiness.R`: Uses config `flashiness_min_days=250` instead of hardcoded 30
- `flow_volumes.R`: Renamed `Q95-Q10` to `Q95_Q10` matching Python/Julia naming
- `baseflow.R`: Paired masking for BFI computation; BFI clamped to [0,1]
- `runoff_ratios.R`: Paired dropna before aggregation (rows where either Q or PPT is NA excluded)
- `elasticity.R`: Removed fallback path; returns NA when insufficient for rolling window (matching Julia)
- `qa_qc.R`: Added 5 missing QA/QC flags (elasticity, runoff ratio, seasonal sum, percentile ordering)

**Known divergences from canonical R** (10 poor columns, all intentional — see docs/CROSS_LANGUAGE_STATUS.md):
- 4 FDC90th columns: `min_days=250` (config) vs `30` (canonical R hardcoded)
- 2 BFI_LyneHollick p-values: Paired masking vs independent `na.rm=TRUE` sums
- 2 elasticity p-values, 2 avg_storage p-values: Canonical R pre-removes NA Q rows before climate merge; rpkg does not

**New benchmark files**:
- `docs/benchmarks/run_rpkg_benchmark.R` — rpkg full signature extraction
- `docs/benchmarks/compare_rpkg.py` — rpkg vs R/Python/Julia comparison
- `docs/benchmarks/rpkg_signatures.csv` — rpkg benchmark output
- `docs/benchmarks/rpkg_comparison.csv` — Detailed per-column R² results
- `rpkg/tests/smoke_test_real_data.R` — Real-data smoke test (10 gages)

### Benchmark Script Fixes and Re-Run (March 2026)

**R Benchmark Script Path Bugs (HIGH — script failed to run)**:
- **Issue**: `run_r_benchmark.R` could not find `config.R` because project root resolution went up only 1 level (`..`) instead of 2 (`../..`) from `docs/benchmarks/`. Output paths (`benchmarks/r_signatures.csv`) were also wrong relative to the corrected working directory. A typo `end_timng` instead of `end_time` caused the script to error after saving results.
- **Fix**: Changed `file.path(script_dir, "..")` to `file.path(script_dir, "../..")`. Updated output paths to `docs/benchmarks/r_signatures.csv` and `docs/benchmarks/r_timing.json`. Fixed `end_timng` typo.
- **Location**: `docs/benchmarks/run_r_benchmark.R`

**DEVELOPMENT.md Benchmark Instructions (MEDIUM)**:
- **Issue**: Instructions said `cd docs/benchmarks` then bare filenames, with no R benchmark command at all.
- **Fix**: Updated to run all scripts from project root with full relative paths (e.g., `Rscript docs/benchmarks/run_r_benchmark.R`).
- **Location**: `docs/DEVELOPMENT.md`

**Cross-Language Comparison Metric Change (March 2026)**:
- **Issue**: Spearman rank correlation only checks if values are monotonically related, not whether they are identical. Two implementations could have perfect Spearman rho but produce systematically different values (e.g., offset or scaled).
- **Fix**: Replaced Spearman with R² of the identity line (y = x) as the primary comparison metric. R² = 1 - SS_res/SS_tot where SS_res = sum((y - x)²). This measures actual value agreement: R² = 1.0 means identical, R² < 1.0 means values differ. Spearman retained as secondary diagnostic.
- **Impact**: More divergences detected — 20 poor columns (was 4 with Spearman). All newly flagged columns are trend statistics (slopes, p-values) where small numerical differences in underlying calculations amplify through trend fitting. Signature means/medians remain near-identical across all 3 languages.
- **Location**: `docs/benchmarks/compare_three_way.py`

**Benchmark Re-Run Results (March 16-17, 2026)**:
- Re-ran all three benchmarks after recent commits (`05ff4be` through `2311433`)
- Julia: 9.2 min (7,369 gages, 551 sig cols)
- Python: 78.9 min (7,369 gages, 551 sig cols)
- R: 874 min (5,707 gages, 551 sig cols) — ran concurrently with Python/Julia, timing inflated
- Three-way comparison (Identity R²): **475 perfect, 56 good, 20 poor**
- Three-way comparison (Spearman rho): 505 perfect, 42 good, 4 poor — identical to previous baseline
- Golden output regression: **551/551 perfect match** against Feb 2026 golden outputs
- No regressions from recent code changes

### Code Quality Improvements Round 2 (March 2026)

**R16: Eckhardt BFI Forward-Fill (HIGH — confirmed: all 3 languages now identical)**:
- **Issue**: R's `eckhardt_filter()` set `baseflow[i] <- NA` when Q[i] is NA, cascading NAs through the recursive filter. Python/Julia were fixed in Round 5 to forward-fill, but R (canonical) was not updated. This was the source of R-Python BFI_Eckhardt rho ~0.989-0.990.
- **Fix**: Applied forward-fill to R: when Q[i] is NA, `baseflow[i] <- baseflow[i-1]`. Also aligned initialization to `min(Q[1] * BFImax, Q[1])` when Q[1] > 0, else 0, matching Python/Julia. All 3 languages now use identical Eckhardt filter logic.
- **Location**: `R/helperFunctions.R` `eckhardt_filter()`
- **Confirmed**: Full benchmark re-run (March 15, 2026) verified R16 produces identical alignment: 505 perfect, 42 good, 4 poor (same as Round 5). BFI_Eckhardt fully resolved across all 3 languages.

**P11: Per-Year Boolean Indexing → groupby (HIGH — Python performance)**:
- **Issue**: 8 modules used `df[df["water_year"] == yr]` in per-year loops — millions of boolean mask operations.
- **Fix**: Replaced with `df.groupby("water_year", sort=False)` in baseflow, recession, flashiness, timing, fdc, pulses, qp_seasonality, storage. Estimated 30-40% Python benchmark speedup.
- **Locations**: All 8 affected modules in `python/streamflow_signatures/`

**P12: Redundant `.copy()` Removal (HIGH — Python performance)**:
- **Issue**: All 11 modules started with `df = streamflow_data.copy()` (unnecessary since functions don't mutate input).
- **Fix**: Changed to `df = streamflow_data` in all 11 modules. Retained `.copy()` inside per-year loops only where year_data is mutated (5 modules).
- **Locations**: All 11 modules in `python/streamflow_signatures/`

**P13: `.loc` Boolean Assignment → List Accumulation (MEDIUM — Python performance)**:
- **Issue**: `flashiness.py` and `fdc.py` used boolean scan per `.loc` assignment per year.
- **Fix**: Replaced pre-allocation + boolean `.loc` with list-of-dicts accumulation + single `pd.DataFrame()` construction.
- **Locations**: `python/streamflow_signatures/flashiness.py`, `python/streamflow_signatures/fdc.py`

**J13: OLS Denominator Guard (HIGH — Julia correctness)**:
- **Issue**: `ols_slope_intercept()` used `denominator == 0` (exact float comparison that almost never triggers).
- **Fix**: Changed to `abs(denominator) < 1e-10` matching `linear_slope()` in stats.jl.
- **Location**: `julia/src/recession.jl`

**J14: DataFrame Pre-allocation (HIGH — Julia performance)**:
- **Issue**: 5 modules used `DataFrame()` + `push!(row; cols=:union)` requiring schema reconciliation per push.
- **Fix**: Pre-allocated DataFrames with `fill(NaN, length(years))` and direct `annual_data[yr_idx, :col]` indexing.
- **Locations**: `julia/src/flow_volumes.jl`, `julia/src/fdc.jl`, `julia/src/timing.jl`, `julia/src/pulses.jl`, `julia/src/runoff_ratios.jl`

**Benchmark Infrastructure Fix**:
- **compare_three_way.py**: Added dot-to-underscore normalization in `normalize_col()` for R's `Q95.Q10` vs Python/Julia's `Q95_Q10`. Previously 8 columns were silently excluded from comparison.
- **R15: run_full_processing.R paths**: Updated to reference `combined_streamflow_data_09feb2026.parquet` (was referencing corrupted Oct 2025 file). Updated `config.R` default `PARQUET_DATA_DIR` to `D:/processedOuts_feb2026`.

**Benchmark Re-Run Results (March 15, 2026 — Post-Round 2 Code Quality)**:
- Re-ran full R pipeline (5h 13m, 5,707 gages, 572 cols), Python benchmark (133 min, 7,369 gages), and Julia benchmark (9.78 min, 7,369 gages)
- Note: R and Python ran concurrently — timings are unreliable for performance comparison. Julia ran solo (9.78 min, comparable to previous 9.6 min).
- Three-way comparison results (551 common columns, 5,707 common gages):
  - R-Python: **4** poor columns (unchanged — all recession pointcloud)
  - R-Julia: **4** poor columns (unchanged — all recession pointcloud)
  - Python-Julia: **0** poor columns (min rho = 0.9976)
  - 505 perfect (>=0.999), 42 good (0.99-0.999), 4 poor (<0.99)
  - BFI_Eckhardt: Confirmed fully resolved from R side (R16 forward-fill working)
  - Results identical to Round 5, confirming R16/P11/P12/P13/J13/J14 did not regress any correlations
- P11/P12 Python performance impact could not be assessed due to concurrent R execution

**GAGES-II Metadata Enrichment in Python/Julia Benchmarks (March 2026)**:
- **Issue**: Python/Julia benchmark outputs lacked the 13 human interference metadata columns that R includes (NDAMS_2009, MAJ_DDENS_2009, STOR_NID_2009, IMPNLCD06, DEVNLCD06, FRESHW_WITHDRAWAL, HYDRO_DISTURB_INDX, CLASS, RHBN, REGULATED, human_interference_class, area_normalized, gage_type).
- **Fix**: Added GAGES-II loading and enrichment to both `run_python_benchmark.py` and `run_julia_benchmark.jl`. Updated `metadata_order` to include all interference columns. RHBN/REGULATED columns are present but empty (HYDAT integration is R-only for now).
- **Bug fixes in metadata module**: Both `python/streamflow_signatures/metadata.py` and `julia/src/metadata.jl` had incorrect file definitions (wrong filenames, wrong delimiter, missing files). Fixed to match R's `config.R`: 5 CONUS files + 4 AKHIPR files, CSV format, `latin-1` encoding for Python, `STAID` read as String in Julia.
- **Locations**: `docs/benchmarks/run_python_benchmark.py`, `docs/benchmarks/run_julia_benchmark.jl`, `python/streamflow_signatures/metadata.py`, `julia/src/metadata.jl`

### Code Quality Improvements Round 1 (March 2026)

Cross-language code review covering ~9,000+ lines across R, Python, and Julia.
Full findings documented in `docs/CODE_REVIEW.md`.

**HIGH priority fixes:**

- **R1: O(n²) metadata lookup** — Replaced row-by-row `rbind` loop (O(n²)) with vectorized `data.table` construction + `setkey()` for O(1) lookup. ~35 second speedup per 7,000-gage run.
  - Location: `R/helperFunctions.R` metadata_lookup construction
- **R2: DRY violations** — Extracted `safe_calculate()` helper function and refactored 11 repeated tryCatch blocks into loop structures in both `process_gages_rawData()` and `process_signatures_from_parquet()`. Reduced ~100 lines.
  - Locations: `R/helperFunctions.R` lines ~214-267 and ~4598-4673
- **P1: Row-wise water year calculation** — Replaced `.apply()` lambda with vectorized pandas datetime arithmetic for ~50x speedup in `dowy` computation.
  - Location: `python/streamflow_signatures/io.py` `add_water_year_columns()`
- **J1: Unchecked findfirst result** — Added `yr_idx === nothing && continue` guard to prevent potential crash on edge-case data.
  - Location: `julia/src/recession.jl` line ~401

**MEDIUM priority fixes:**

- **R3: FDC magic number** — Moved hardcoded `1e-10` flow floor to `FDC_FLOW_FLOOR` in `config.R`.
- **R4: Hardcoded BATCH_SIZE** — Moved to `PROCESSING_BATCH_SIZE` in `config.R`; removed unnecessary mid-loop `gc()` call.
- **R5: Per-metric success tracking** — Added `metric_success` counters to `process_signatures_from_parquet()` with summary report at end of processing.
- **R6: Dead code removal** — Moved `generate_streamflow_dt_og()` to `archive/deprecated_generate_streamflow_dt_og.R` (contains the 99999 bug).
- **R7: Missing docstring** — Added documentation for `fit_sinusoidal_model()` explaining the sinusoidal formula and phase-to-minimum-day conversion.
- **P5: Hardcoded benchmark paths** — Made `STREAMFLOW_PATH`, `CLIMATE_PATH`, `METADATA_PATH` configurable via environment variables in `run_python_benchmark.py`.
- **J4: Double Float64 conversion** — Removed redundant `copy()` in `flashiness.jl` (`Float64.(Q)` already creates a new array).

**Bug fix (found during verification):**

- **R 4.5.1 sprintf compatibility** — `sprintf("%d", as.numeric(id))` fails in R 4.5.1 when `as.numeric()` returns a double (stricter type checking than R 4.3). Affected both the Daymet ID matching loop and `find_metadata()`. Replaced with `paste0(strrep("0", num_zeros), id)` which is simpler and works for all ID types including non-numeric Canadian gage IDs.
  - Locations: `R/helperFunctions.R` Daymet ID loop and `find_metadata()`

### Cross-Language Alignment Round 5 — Final 4 Columns (March 2026)

**Python Eckhardt Filter NaN Forward-Fill (3 cols, BFI_Eckhardt trend stats)**:
- **Issue**: BFI_Eckhardt trend stats (linear_slp, mk_pval, spearman_pval) had R-Py rho ~0.989-0.990; R-Jl = 1.000
- **Root Cause**: Python's Eckhardt filter cascaded NaN (one NaN Q → all subsequent baseflow = NaN), while Julia forward-fills baseflow on NaN Q (no cascade). R also cascades, but R's `sum(Q, na.rm=TRUE)` denominator includes post-gap Q values, while Python's paired masking excluded them — creating a denominator mismatch. The previous Round 4 attempt to fix the denominator was wrong (made things worse). The correct fix is to eliminate the cascade.
- **Fix**: Changed Python's `eckhardt_filter()` to forward-fill baseflow when Q is NaN (matching Julia). Also aligned initialization: `baseflow[0] = min(BFImax * Q[0], Q[0])` when Q[0] > 0, else 0.0 (matching Julia). With forward-fill, paired masking becomes equivalent to total valid Q since baseflow is never NaN where Q is valid.
- **Location**: `python/streamflow_signatures/baseflow.py` lines 46-67

**Recession Pointcloud Near-Singularity Check (4 cols, Python + Julia)**:
- **Issue**: Recession pointcloud p-values had R-Py rho ~0.85-0.91. Previously deemed "irreducible library-level OLS differences."
- **Root Cause**: R's `lm()` uses QR decomposition with rank checking — it rejects near-singular design matrices where `var(log(Q))` is near zero (low-flow plateaus). Python's `linregress()` (LAPACK SVD) and Julia's custom OLS handle these gracefully, returning extreme slopes instead of failing. For marginal gages (3-4 valid pointcloud years), R fails on 1-2 near-singular years → n < 3 → MK returns p=1.0, while Python/Julia succeed → real p-values.
- **Fix**: Added `var(log_Q) < 1e-8` check before OLS fitting in both Python and Julia pointcloud analysis. When variance is below threshold, the year is skipped (remains NA), matching R's `lm()` + `tryCatch` behavior. Threshold tuned from initial 1e-10 to 1e-8 (closer to R's `.Machine$double.eps^0.5 ≈ 1.49e-8`). Improved Py-Jl agreement but could not fully replicate R's QR rank-checking — these 4 columns remain irreducible.
- **Locations**: `python/streamflow_signatures/recession.py` line 388, `julia/src/recession.jl` line 409
- **Note**: Only applied at pointcloud call site, NOT in `ols_slope_intercept()` itself (used by event fitting where the issue doesn't apply)

### Cross-Language Alignment Round 4 — Final Column Fixes (March 2026)

**Constant-Series Mann-Kendall — Reverted (LOW — investigated, made things worse)**:
- **Issue**: R's `Kendall::MannKendall()` returns tau=0, p=1 for constant series. Python scipy returns NaN; Julia returns NaN when var_S <= 0.
- **Attempted Fix**: Added constant-series check returning (0.0, 1.0) in Python and Julia.
- **Reverted**: Benchmarking showed this WORSENED correlations — added 5 new poor columns (n_low_pulses_all/year mk_rho, Q1/Q5/Q10 mk_rho). Root cause: Python/Julia produce constant series for some gages where R computes non-constant values (subtle differences in underlying signature calculations). Before the fix, NaN excluded these edge-case gages from correlation. After the fix, tau=0.0 conflicted with R's real tau values.
- **Status**: Reverted to NaN behavior. The constant-series gages are excluded from Spearman correlation, which is the correct behavior when the underlying signature values differ.

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

**Julia Flashiness NaN Interpolation (MEDIUM — 1 gage outlier, affects all flashiness cols)**:
- **Issue**: Gage 02244440 had Julia flashiness ~26x higher than R/Python (0.496 vs 0.019). Systematic assessment found this was the only gage-level outlier in flashiness.
- **Root Cause**: Julia removed NaN values and compacted the array before computing `diff()`, creating artificial jumps between non-adjacent days. R uses `approx()` (linear interpolation) and Python uses `np.interp` to fill NaN gaps, preserving temporal adjacency. Additionally, Julia was missing the `max_missing_frac > 0.2` check that R and Python apply.
- **Fix**: Rewrote `calculate_flashiness()` to interpolate NaN values using linear interpolation (matching R's `approx()` and Python's `np.interp`). Added `na_frac > 0.2` guard in `analyze_flashiness_trends()` matching R/Python.
- **Location**: `julia/src/flashiness.jl` lines 24-64, 95-107

**FDC90th Config min_days — Partially Reverted (MEDIUM — investigated)**:
- **Issue**: Python produced 42 extra NaN values for FDC90th where R had valid values. R uses hardcoded `min_days=30` while config had `min_days=250`.
- **Attempted Fix**: Changed `config/signatures_config.json` `fdc.min_days` from 250 to 30.
- **Partially Reverted**: Benchmarking showed changing min_days to 30 added FDC90th_mk_pval as a new poor column (R-Py rho 0.988). Years with only 30-249 days have very few data points in the 90th percentile exceedance range, making FDC90th slopes noisy. Config reverted to min_days=250.
- **Kept**: Python negative Q filter (`q_values = q_values[q_values >= 0]`) before log10 — this is a genuine bugfix preventing `-inf` from `linregress`. Julia `fdc.jl` now uses `CFG_FDC_MIN_DAYS` instead of hardcoded 250 (both resolve to 250, but config-driven is cleaner).
- **Locations**: `python/streamflow_signatures/fdc.py` line 79, `julia/src/fdc.jl` line 84

**Recession Pointcloud P-Values — No Fix (4 cols, irreducible)**:
- **Investigation**: Thorough comparison found NO algorithmic difference across all 3 languages. 286 of 975 gages have R producing mk_pval=1.000 (n < 3 valid pointcloud years) while Python/Julia compute real p-values (n=3-4). Caused by implementation-level differences in how R's `lm()` vs Python's `linregress()` vs Julia's OLS handle edge-case fits. These 4 columns will remain < 0.99.

**Benchmark Re-Run Results (March 14, 2026 — Post-Round 5)**:
- Re-ran Python benchmark (69 min, 7,369 gages) and Julia benchmark (9.6 min, 7,369 gages) after Round 5 fixes
- Three-way comparison results:
  - R-Python: 6 → **4** poor columns (3 BFI_Eckhardt resolved by forward-fill)
  - R-Julia: 4 → **4** poor columns (unchanged — all recession pointcloud)
  - Python-Julia: 3 → **0** poor columns (BFI_Eckhardt resolved → perfect pairwise alignment)
  - BFI_Eckhardt: All 3 cols now rho >= 0.99 (was 0.914-0.979 — fully resolved)
  - Recession pointcloud: 4 cols remain at rho ~0.84-0.91 (near-singularity guard improved Py-Jl but R divergence irreducible)
- 505 of 551 columns have rho >= 0.999 across all 3 pairs (Perfect)
- 42 additional columns >= 0.99 (Good)
- Only 4 columns remain below 0.99 (all recession pointcloud p-values)
- Median rho = 1.000 for all 3 pairs
- Python-Julia: min rho = 0.9976 (was 0.9137) — **zero poor columns**
- Benchmark outputs saved to `docs/benchmarks/`

**Benchmark Re-Run Results (March 13, 2026 — Post-Round 4)**:
- Re-ran Python benchmark (82 min, 7,369 gages) and Julia benchmark (14 min, 7,369 gages) after Round 4 fixes
- Three-way comparison results:
  - R-Python: 7 → **6** poor columns (qp_slope_sd resolved, BFI_Eckhardt_spearman_pval crossed 0.99)
  - R-Julia: 5 → **4** poor columns (qp_slope_sd resolved)
  - Python-Julia: 3 → **3** poor columns (unchanged — all BFI Eckhardt borderline)
  - qp_slope_sd: ALL rho > 0.99 (was 0.990 — fully resolved by mid-point + ddof fix)
  - BFI_Eckhardt: R-Py rho 0.989-0.990 (unchanged — denominator fix was reverted)
- 499 of 551 columns have rho >= 0.999 across all 3 pairs (Perfect)
- 45 additional columns >= 0.99 (Good)
- Only 7 columns remain below 0.99 (4 recession pointcloud + 3 BFI Eckhardt)
- Median rho = 1.000 for all 3 pairs
- Benchmark outputs saved to `docs/benchmarks/`

### Cross-Language Alignment Round 3 — Recession, BFI, Stats Pre-Filtering (March 2026)

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

### Cross-Language Alignment — Config, Metadata, and Divergence Fixes (March 2026)

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

### Three-Way Comparison Infrastructure (March 2026)

- Fixed `compare_three_way.py` `normalize_col()` to handle FDC naming mismatch (`FDC_all` -> `FDCall`, `FDC_90th` -> `FDC90th`, `FDC_mid` -> `FDCmid`)
- Added missing metadata columns to `META_COLS`: `start_water_year`, `end_water_year`, `country`
- Added `generate_three_way_report()` function producing comprehensive markdown report
- Deprecated `compare_results.py` (Python-vs-Julia only) in favor of `compare_three_way.py`
- Updated DEVELOPMENT.md benchmark section with fresh three-way results

### Cross-Language Alignment Round 2 (March 2026)

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

### Cross-Language Alignment Round 1 (March 2026)

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

### Deprecated (March 2026)
- **Corrupted Parquet File**: Marked `D:/combined_streamflow_output/combined_streamflow_data.parquet` (October 2025) as **DO NOT USE**
  - **Issue**: Contains 99999 multiplier bug for Canadian gages without basin area
  - **Impact**: Q values ~100,000x too high for affected gages (e.g., 08ND025: 78.5M instead of 785.6)
  - **Use instead**: `D:/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet`
  - **Documentation**: Added "Parquet Data Files" section to DEVELOPMENT.md with active vs deprecated file list

### Terminology Fix (March 2026)
- `analyze_flow_timing_trends()` — Renamed internal variable `julday_max` to `timing_by_year` and fixed comments to correctly say "day of water year" instead of "Julian day" (code was already correct, but naming was misleading)

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

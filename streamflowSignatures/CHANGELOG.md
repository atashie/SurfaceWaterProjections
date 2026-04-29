# Changelog

All notable changes to the Streamflow Signatures project.
For historical entries (Dec 2025 – March 2026), see [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

## [Unreleased]

### Planned
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)
- Port data ingestion utilities to Julia (long-term — currently R-only via dataRetrieval/tidyhydat)
- Generate Julia golden outputs (624 cols, 7,313 gages)
- BFImax estimation via Collischonn & Fan (2013) backward filter — would give BFI_Eckhardt_param per-gage BFImax instead of fixed 0.8, improving discriminating power (currently range [0.47–0.80] due to BFImax saturation)

### Changepoint Detection — Pettitt Test (April 2026)

Non-parametric Pettitt changepoint test applied to all time-series signatures, with differential metrics.

**Pettitt Test** (`pettitt_test`). Non-parametric rank-based test (Pettitt 1979). Test statistic K = max|U_t| where U_t accumulates pairwise sign comparisons. P-value via asymptotic approximation. Always identifies a location; use p-value for significance.

**Differential Metrics** (`segment_differential_metrics`). Applied at the Pettitt changepoint:
- delta_mean (post − pre), pct_change (%), pre/post Mann-Kendall p-values

Changepoint search window: WY 1980-2024, minimum 20 non-NA observations total, minimum 10 per segment. Independent of the 80% trend_completeness gate.

**8 new columns per signature (76 signatures × 8 = 608 changepoint columns):**

| Columns | Description |
|---------|-------------|
| `pettitt_cp_year`, `pettitt_pval`, `pettitt_pre_mean`, `pettitt_post_mean`, `pettitt_delta_mean`, `pettitt_pct_change`, `pettitt_pre_mk_pval`, `pettitt_post_mk_pval` | Non-parametric test + differential |

**Design decisions:**
- Per-signature changepoints (not per-gage) — different processes change at different times
- Pettitt O(n²) is fast for annual series (n ≤ 45)
- BIC 4-model comparison was initially implemented alongside Pettitt but removed to reduce clutter (code retained in `changepoint.jl` but not called from `generate_stats()`)

**Files modified:**
- `julia/src/changepoint.jl`: `pettitt_test()` + `segment_differential_metrics()` (+ unused BIC engine)
- `julia/src/stats.jl`: `generate_stats()` with `changepoint` kwarg, `_run_changepoint_block!()` helper, `CP_SUFFIXES` constant (8)
- `julia/src/signatures.jl`: Passes changepoint config to all signature functions
- `julia/src/StreamflowSignatures.jl`: Exports for `pettitt_test`, `segment_differential_metrics`
- `julia/src/config.jl`: `CFG_CP_*` constants
- `config/signatures_config.json`: `changepoint` section
- All 13 signature function files: `changepoint` kwarg propagation
- `docs/benchmarks/run_julia_benchmark.jl`: Changepoint config NamedTuple
- `julia/test/test_changepoint.jl`: 78 tests (39 unit + 39 integration)
- `config.R`: `CHANGEPOINT_SUFFIXES` (8 entries)
- `docs/benchmarks/build_changepoint_dashboard.py`: Pettitt-only dashboard

**Julia benchmark (April 28, 2026)**: 7,313 gages, 1,264 columns (656 base/metadata + 608 changepoint). 76 signatures × 8 Pettitt fields.

**Pettitt test results:**
- Overall: ~13% of evaluations have p < 0.05 (69,279 / 516,992); 2.7x the 5% null expectation
- CP year range: [1989, 2014] — respects 10-year minimum segments within WY 1980-2024

**Signal robustness analysis (April 29, 2026):**
- Significance rates vary by category: Flow Timing 3.7% (below null — stationary), Flashiness 19.4%, Baseflow 15.1%, Flow Percentiles 18.4%, Elasticity 25.8% (inflated by short rolling-window series)
- Effective independence ~17/76 signatures (77% redundancy) — flow percentiles are highly correlated with Qann
- After BH-FDR correction at 0.05: ~3.5% survive — a core of robust detections
- CP year clustering: ~50% of significant detections in 1998–2006; partly center-of-window bias, partly real hydroclimatic signal
- elasticity_rolling at ~49% significance is a methodological artifact from short series (n=15–25), not a hydrological finding
- pct_change discriminates well: median |pct_change| ~40% (sig) vs ~15% (non-sig), 2.65x ratio
- Pre/post MK p-values largely uninformative (92% show no within-segment trend)
- Asymptotic p-value is conservative for small n: significance rate 4.3% for n=20–25 (below 5% null)

**References:**
- Pettitt, A.N. (1979). A non-parametric approach to the change-point problem. Applied Statistics, 28(2), 126-135.

### Recession-Parameterized Baseflow Signatures (April 2026)

Three new signatures using recession-derived filter parameters (Collischonn & Fan 2013, Eckhardt 2005). Implemented in all three languages (Julia canonical April 24, Python and rpkg ported April 25).

**Signature 1: `alpha_linear`** (8 stats via generate_stats)
- Discrete recession constant under linear reservoir assumption (b=1)
- Per-year computation: median of Q_{i+1}/Q_i ratios from point-cloud recession data
- First day of each recession event removed (storm peak influence)
- Alpha collection independent of power-law fit success — depends only on raw Q pairs
- Per-year threshold: >10 valid alpha pairs; same 25-event gate as other recession metrics for trend stats

**Scalar: `recession_alpha_point_cloud_linear_reservoir`** (per-gage diagnostic)
- Whole-record median of Q_{i+1}/Q_i ratios across all recession events
- Computed independently of 25-event gate (only requires >10 alpha pairs)
- Used to parameterize BFI_Eckhardt_param and BFI_LyneHollick_param

**Signature 2: `BFI_Eckhardt_param`** (8 stats)
- Eckhardt BFI with recession-derived `a = recession_alpha_point_cloud_linear_reservoir`
- BFImax fixed at 0.8 (not estimated from backward filtering)

**Signature 3: `BFI_LyneHollick_param`** (8 stats)
- Lyne-Hollick BFI with recession-derived `alpha = recession_alpha_point_cloud_linear_reservoir`
- Heuristic parameterization — L-H alpha has no physical derivation from recession (Nathan & McMahon 1990)

**Total: 24 new stat columns + 1 per-gage scalar = 25 new columns (624 total signature columns)**

**Julia benchmark validation (April 24, 2026)**: 7,313 gages, 15.0 min (8.1 gages/s). 598 of 599 common columns Perfect (R² >= 0.999), 1 N/A (<10 paired values). Zero regression confirmed.

| Metric | Non-NA gages | Range |
|--------|-------------|-------|
| `alpha_linear_median` | 5,825 / 7,313 | [0.273, 0.985] |
| `recession_alpha_point_cloud_linear_reservoir` | 7,131 / 7,313 | [0.273, 0.997] |
| `BFI_Eckhardt_param_mean` | 7,131 / 7,313 | [0.474, 0.804] |
| `BFI_LyneHollick_param_mean` | 7,131 / 7,313 | [0.255, 0.966] |

**Cross-language validation (April 27, 2026)**: Python and rpkg benchmarks re-run with 624 columns, compared against Julia golden.

| Pair | Perfect (>=0.999) | Good (0.99-0.999) | Poor (<0.99) | Min R² |
|------|-------------------|-------------------|-------------|--------|
| Python vs Julia | **615** | 5 | 3 | 0.984 |
| rpkg vs Julia | **614** | 4 | 5 | 0.969 |

- Python: all 25 new columns Perfect. 3 poor are existing recession seasonality minimum (sinusoidal fit sensitivity). 150.4 min, 0.81 gages/s.
- rpkg: 24 of 25 new columns Perfect. 1 new poor: `alpha_linear_spearman_pval` (R² = 0.971, Spearman p-value library difference — R `cor.test` vs Julia t-approximation on tied annual values). 4 pre-existing poor unchanged. 244.6 min (I/O contention).
- Zero regression on old columns confirmed: excluding 25 new columns, Python has 590 Perfect / 3 Poor (same 3 as before), rpkg has 590 Perfect / 4 Poor (same 4 as before).

**Validation dashboards**: `docs/benchmarks/build_new_vs_golden_julia_dashboard.py` (Julia new vs golden), `docs/benchmarks/build_julia_golden_dashboard.py python` (Python vs Julia golden), `docs/benchmarks/build_julia_golden_dashboard.py rpkg` (rpkg vs Julia golden).

**Known limitation — BFI_Eckhardt_param narrow range**: BFI_Eckhardt_param_mean has range [0.474, 0.804] with std=0.036, much narrower than fixed-parameter BFI_Eckhardt_mean (range [0.139, 0.802], std=0.119). This is a mathematical property of the Eckhardt filter: with BFImax=0.8 fixed, the `a` parameter only controls dynamics, while BFImax controls the steady-state level. Since 92% of gages have recession alpha < 0.95, the filter saturates at BFImax. BFI_LyneHollick_param is unaffected (range [0.255, 0.966], std=0.117). Future improvement: estimate BFImax per gage via Collischonn & Fan (2013) backward filter. See `docs/SIGNATURES.md` for details.

**Files modified:**
- `julia/src/recession.jl`: Added alpha_linear per-year computation and whole-record scalar to `analyze_recession_parameters()`
- `julia/src/baseflow.jl`: Added `analyze_baseflow_indices_with_parameters()` function
- `julia/src/signatures.jl`: Wired recession scalar extraction → parameterized BFI call
- `julia/src/StreamflowSignatures.jl`: Added export
- `python/streamflow_signatures/recession.py`: Ported alpha_linear + scalar (April 25)
- `python/streamflow_signatures/baseflow.py`: Ported `analyze_baseflow_indices_with_parameters()` (April 25)
- `python/streamflow_signatures/signatures.py`: Wired orchestration (April 25)
- `python/streamflow_signatures/__init__.py`: Added export (April 25)
- `rpkg/R/recession.R`: Ported alpha_linear + scalar (April 25)
- `rpkg/R/baseflow.R`: Ported `analyze_baseflow_indices_with_parameters()` (April 25)
- `rpkg/R/signatures.R`: Wired orchestration (April 25)
- `rpkg/NAMESPACE`: Added export (April 25)
- `config.R`: Registered new signature bases and scalar exception
- `docs/SIGNATURES.md`: Documented new metrics in Baseflow and Recession sections
- `docs/DEVELOPMENT.md`: Updated column counts (594 → 624), file structure, benchmark files table

### Canonical Language Transition (April 2026)

Julia replaced R as the canonical implementation. All three synced implementations (Julia, Python, rpkg) are in equilibrium (99.5% R² agreement, 594 signature columns, 7,313 gages). Julia is ~40x faster than monolithic R (~27 min vs ~18 hrs), has cleaner modular architecture (18 files vs monolithic helperFunctions.R), and already served as first-mover for Section 3 changes.

**Changes:**
- Julia (`julia/src/`) is now the canonical reference for all signature implementations
- rpkg (`rpkg/`) is the production-ready R port (replaces monolithic R for new development)
- Monolithic R (`R/helperFunctions.R`) deprecated as legacy shim — retained for data ingestion utilities and existing callers, but no new signature features will be added
- Data ingestion stays in R (dataRetrieval/tidyhydat — separate concern from signatures)
- Julia golden outputs to be generated alongside historical R golden outputs (Feb 2026)
- New Julia-primary comparison scripts for validating Python and rpkg against Julia golden outputs
- Recession algorithm sync to rpkg only (monolithic R skipped — deprecated)

**Documentation updated:**
- `CLAUDE.md`: Canonical references, multi-language table, change workflow, code status
- `docs/DEVELOPMENT.md`: Design principles, adding new signatures workflow, benchmark framing
- `README.md`: Quick-start, language table
- `docs/CROSS_LANGUAGE_STATUS.md`: Reframed baseline
- `docs/SIGNATURES.md`: Reference updates
- Package READMEs (`julia/`, `python/`, `rpkg/`): Canonical/port status notes
- `docs/WORKFLOW_REVIEW.md`, `docs/CODE_REVIEW.md`: Historical banners
- Skills (`claude-skill/cross-language-alignment.md`, `claude-skill/streamflow-signatures.md`): Julia-canonical guardrails
- Deprecation headers on `R/helperFunctions.R`, `run_full_processing.R`, `run_caravan_processing.R`, `run_restricted_processing.R`

### Julia Experiment Infrastructure (April 2026)

Sensitivity experiment framework for Julia benchmarks. Three experiments test how restricting the analysis period and gage quality affect results.

**Configuration** (`config/signatures_config.json`):
- Added `filtering.start_water_year` (null default) — restrict analysis to water years >= specified year
- Added `filtering.min_qualifying_data_fraction` (null default) — require minimum fraction of possible years with qualifying data

**Infrastructure**:
- `julia/src/config.jl`: ENV-based overrides (`STREAMFLOW_START_WATER_YEAR`, `STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION`, `STREAMFLOW_CONFIG`, `STREAMFLOW_OUTPUT_PREFIX`) with JSON fallback
- `docs/benchmarks/run_julia_benchmark.jl`: Added filtering logic, parameterized output paths, zero-gage guard. Experiment controls read from ENV at runtime (not module constants) to avoid Julia precompilation caching of `const` values
- `docs/benchmarks/run_julia_benchmark_startIn1993.jl`: Thin wrapper — WY >= 1993
- `docs/benchmarks/run_julia_benchmark_startIn1993_60pct.jl`: Thin wrapper — WY >= 1993 + 60% qualifying fraction
- `docs/benchmarks/run_julia_benchmark_startIn1993_80pct.jl`: Thin wrapper — WY >= 1993 + 80% qualifying fraction
- `docs/benchmarks/compare_experiment_vs_julia.py`: Parameterized comparison (CSV + markdown report) with gage diff analysis, years-per-gage stats, 6-tier R² classification
- `docs/benchmarks/build_experiment_vs_julia_dashboard.py`: Parameterized interactive HTML dashboard (dual Leaflet maps + Plotly scatterplot)

**Results (April 16-17, 2026)**:

| Experiment | Gages | Dropped | Time | Rate |
|------------|-------|---------|------|------|
| Baseline | 7,313 | — | 27.5 min | 4.4/s |
| startIn1993 | 6,678 | 635 | 24 min | 4.7/s |
| startIn1993+60pct | 6,579 | 734 | 18 min | 5.9/s |
| startIn1993+80pct | 5,431 | 1,882 | 14 min | 6.5/s |

- 635 gages dropped by WY>=1993 restriction (485 USGS, 150 Canadian — lost pre-1993 years below 20-year min)
- 99 additional gages dropped by 60% filter (all had exactly 20 qualifying years; `frac=20/34=0.588` due to partial WY2026 in raw parquet)
- 1,148 additional gages dropped by 80% filter (1,553 USGS + 329 Canadian total vs baseline). With denominator=34 (WY1993–2026), gages need >= 28 valid years (28/34=0.824 passes; 27/34=0.794 fails)
- Trend statistics diverge significantly vs baseline (expected — different period-of-record); means/medians stable

**Bug found and fixed**: Julia precompilation caches `const` values in `.ji` files. ENV vars set by experiment wrappers were silently ignored after first compilation. Fixed by reading ENV vars at runtime in `main()` as local variables.

### Guidelines vs Implementation Audit (April 2026)

Systematic audit comparing `SIGNATURE_GUIDELINES.md` against Python, Julia, and rpkg codebases.

**Season exclusion year counts (HIGH — Finding #1)**: Added 4 per-gage scalar diagnostics: `season_excluded_years_winter`, `season_excluded_years_spring`, `season_excluded_years_summer`, `season_excluded_years_fall`. Counts of years where each season failed the 80% completeness threshold. Guidelines explicitly requested this but it was never implemented. Added to Python (`signatures.py`), Julia (`signatures.jl`), rpkg (`signatures.R`). Registered in `config.R` as `EXPECTED_SEASON_EXCLUDED_YEARS`.

**Q95_Q10 documentation fix (MEDIUM — Finding #5)**: Fixed "ratio" → "difference" in code comments and documentation. The metric computes Q95 - Q10 (subtraction), not Q95/Q10 (division). Fixed in: `config.R`, `SIGNATURES.md`, `julia/src/flow_volumes.jl`, `rpkg/R/flow_volumes.R`, `rpkg/man/calculate_flow_vols_by_year.Rd`. Python already correct.

**Trend completeness exemptions documented (LOW — Finding #12)**: Added documentation in `SIGNATURES.md` explaining why recession and elasticity are exempt from the 80% trend completeness gate.

**Guidelines items requiring Google Doc update** (no code changes — user to update):
- Finding #4: Runoff PPT thresholds (10mm annual, 1mm seasonal) stricter than guidelines' zero-only spec
- Finding #7: Qann glossary says "mean" but code computes "total" (sum of daily mm/day)
- Finding #8: `elasticity` → `elasticity_rolling` rename not reflected in glossary
- Finding #9: Negative Q — guidelines say both "remove" and "count"; code resolves by config-driven rejection (off by default) + `negative_ann` count
- Finding #10: Constant-SD — guidelines say both "remove" and "flag"; code flags only (per user decision)
- Finding #11: "Need to fix: julian" note is stale — all codebases use `D*_day` names

**User decisions (documented in `SIGNATURES.md`)**:
- Finding #2: Pulse `*_year` variants (per-year percentiles) — **kept and documented**. Complement the guidelines-specified `*_all` (period-of-record) variants.
- Finding #3: Seasonal flow reversals (`Flow_Reversals_winter/spring/summer/fall`) — **kept and documented**. Complement the annual count.
- Finding #6: Q95_Q10 column name — **underscore (`Q95_Q10`) is the standard**. R canonical still outputs `Q95-Q10` (hyphen); will be updated when R is synced.

**rpkg flashiness dead ad-hoc NA interpolation (HIGH — Audit Finding #1)**: Removed inline `approx(rule=2)` interpolation from `rpkg/R/flashiness.R:26-28`. This was removed from Python, Julia, and R canonical in April 2026 (CHANGELOG Fix 6) but rpkg was missed. `rule=2` does constant boundary extrapolation, violating centralized NA design. Replaced with `stopifnot(!any(is.na(q_values)))` defensive assertion.

**rpkg flow reversals dead ad-hoc NA interpolation (HIGH — Audit Finding #2)**: Removed inline `approx(rule=2)` interpolation from `rpkg/R/pulses.R` `count_flow_reversals()` (lines 22-28). Python, Julia, and R canonical are all clean. Replaced with `stopifnot(!any(is.na(flow_vector)))` + direct assignment.

**QA/QC flag integration into calculate_all_signatures() (MEDIUM — Audit Finding #3)**: Added optional `include_qa_flags` parameter to `calculate_all_signatures()` in all 3 packages (rpkg, Python, Julia). When TRUE/True/true, appends 12 QA/QC flag columns from `compute_qa_flags()` to output. Previously QA flags were only called from benchmark runners, not the public API.

**Recession-informed BFI formally deferred (MEDIUM — Audit Finding #4)**: Guidelines describe `analyze_baseflow_indices_with_parameters()` (Collischonn & Fan 2013). Not implemented in any codebase. Added "Not Yet Implemented" note in `SIGNATURES.md`.

**Elasticity 30% diagnostic deferred (MEDIUM — Audit Finding #5)**: Guidelines request counterfactual "what if <30% missing" diagnostic. Pending domain expert clarification. Documented in `SIGNATURES.md`.

**Guidelines items requiring Google Doc update** (second round — no code changes):
- Audit Finding #6: Qann glossary says "mean" but code computes "total" (sum of daily mm/day); same for seasonal
- Audit Finding #7: `elasticity` → `elasticity_rolling` rename not in glossary
- Audit Finding #8: Missing `n_recession_events` in recession glossary
- Audit Finding #9: Missing `elasticity_annual` + "elastacity" typo in glossary
- Audit Finding #10: PPT threshold "zero" → "10mm annual / 1mm seasonal"
- Audit Finding #11: Negative Q contradictory ("remove" vs "count") — code resolves via config
- Audit Finding #12: Constant-SD contradictory ("remove" vs "flag") — code flags only
- Audit Finding #13: Stale "Need to fix: julian" note — all codebases use `D*_day`
- Audit Finding #14: `elasticity_static` output description ambiguous

### Section 3 Sync to Python & rpkg + Benchmark Re-Run (April 2026)

Synced all 7 Section 3 changes from Julia to Python and rpkg. Re-ran all 3 benchmarks and validated.

**Julia runoff_ratios fix (LOW)**: Changed `>=` to `>` for PPT threshold in `julia/src/runoff_ratios.jl` (lines 90, 105). Julia was using `>=` for `CFG_RUNOFF_MIN_ANNUAL_PPT` and `CFG_RUNOFF_MIN_SEASONAL_PPT`; R and Python use strict `>`. Two-character fix.

**Python Section 3 sync** (4 files):
- `python/streamflow_signatures/recession.py`: Replaced look-ahead algorithm with position-marking; added `n_recession_events` per-year count (independent of min_events gate)
- `python/streamflow_signatures/elasticity.py`: Renamed `elasticity` → `elasticity_rolling`; added `elasticity_annual` (year-over-year); added `elasticity_years_total` and `elasticity_years_low_ppt` diagnostics
- `python/streamflow_signatures/runoff_ratios.py`: Added `runoff_ratio_high_count` per-gage scalar
- `python/streamflow_signatures/signatures.py`: Added `calculate_negative_days` call

**rpkg Section 3 sync** (6 files):
- `rpkg/R/recession.R`: Replaced look-ahead algorithm with position-marking; added `n_recession_events`
- `rpkg/R/elasticity.R`: Added `elasticity_annual`; added `elasticity_years_total` and `elasticity_years_low_ppt`; fixed NA template
- `rpkg/R/runoff_ratios.R`: Added `runoff_ratio_high_count` per-gage scalar
- `rpkg/R/pulses.R`: Added `calculate_negative_days()` function
- `rpkg/R/signatures.R`: Added `calculate_negative_days` call
- `rpkg/NAMESPACE`: Exported `calculate_negative_days`
- `rpkg/inst/config/signatures_config.json`: Synced with project-level config (was severely outdated — missing `use_legacy_filtering: false`, `na_handling` section, D1/D99 percentiles)

**Comparison script updates**:
- `docs/benchmarks/compare_three_way.py`: Updated `categorize_metric()` for new Section 3 columns (n_recession, negative_ann, D1/D99)
- `docs/benchmarks/compare_rpkg.py`: Same categorization updates

**Julia benchmark runner fixes**:
- `julia/src/StreamflowSignatures.jl`: Added `load_canadian_interference` to module exports
- `docs/benchmarks/run_julia_benchmark.jl`: Fixed `Vector{Missing}` → `Vector{Union{Missing, Bool}}` for RHBN/REGULATED columns

**Benchmark results (April 14-15, 2026)** — All 3 synced implementations (rpkg, Python, Julia) now produce 594 signature columns across 7,313 gages:

| Pair | Perfect (>=0.999) | Good (0.99-0.999) | Poor (<0.99) | Min R² |
|------|-------------------|-------------------|-------------|--------|
| Python vs Julia | 619 | 5 | 3 | 0.986 |
| rpkg vs Julia | 586 | 4 | 4 | 0.967 |
| rpkg vs Python | 582 | 5 | 7 | 0.967 |

Poor columns are all irreducible library-level differences: recession pointcloud p-values (2), FDC90th p-values (2), recession seasonality minimum (3). Python-Julia agreement: 0 columns below 0.95.

rpkg vs R canonical: 490 perfect, 46 poor (all recession — R canonical still has old look-ahead algorithm).

**Net column impact**: 594 signature columns per gage (was 551 pre-Section 3). +16 D1/D99, +8 n_recession_events, +8 elasticity_annual, +8 negative_ann, +2 elasticity diagnostics, +1 runoff_ratio_high_count = +43 new; elasticity renamed (8 cols).

### Julia vs Golden R Comparison & Fixes (April 2026)

**Julia vs Golden R comparison tooling** — New comparison pipeline for validating Julia post-Section 3 output against Feb 2026 Golden R reference:
- `docs/benchmarks/compare_julia_vs_golden_r.py`: 6-tier R² comparison across 551 common columns and 5,697 common gages
- `docs/benchmarks/build_julia_vs_golden_r_dashboard.py`: Interactive HTML dashboard with dual Leaflet maps and Plotly scatterplot
- `docs/benchmarks/julia_vs_golden_r_summary.md`: Detailed markdown report with per-category and per-stat breakdowns

**Elasticity operator fix (LOW)**: Changed `>=` to `>` for `min_annual_ppt` filter in `julia/src/elasticity.jl`. R and Python use strict `>` (excluding years where `P_annual == 10mm`); Julia was using `>=` (including them). Single-character fix aligns all 3 languages.

**Canadian HYDAT metadata for Julia** — Previously all 1,333 Canadian gages had `human_interference_class = "unknown"` with empty RHBN/REGULATED columns in Julia output:
- New `R/export_hydat_metadata.R`: One-time script to export 8,012 Canadian stations (RHBN + REGULATED) from tidyhydat to CSV
- Exported to `metadata/canadian_hydat_interference.csv` (tracked in git for cross-language use)
- Updated `config/signatures_config.json`: `hydat_path` from `null` to `"metadata/canadian_hydat_interference.csv"`
- Updated `julia/src/metadata.jl`: Relative path resolution against project root in `load_canadian_interference()`
- Updated `docs/benchmarks/run_julia_benchmark.jl`: Calls `load_canadian_interference()` and merges RHBN/REGULATED/human_interference_class for Canadian gages

**Julia vs Golden R divergence analysis** — Root-cause investigation of 288 divergent columns (R² < 0.99):
- **Recession (46 cols)**: Intentional algorithm change (position-marking vs look-ahead). Expected; R/Python sync pending.
- **trend_completeness (220 cols)**: Julia applies 80% non-NA gate for trend stats; Golden R (Feb 2026) predates this feature. Affects slopes/p-values across all categories; means/medians unaffected. Expected; resolves when Golden R is regenerated.
- **Elasticity (9 cols)**: `>=` vs `>` operator bug (fixed above) + Golden R's since-removed min_Q filter.
- **dur_low_pulses_all (6 cols)**: Julia produces 3,469 more NAs than R — investigation pending.

### Guidelines Section 3 Alignment — Julia-First (April 2026)

Six changes implementing Guidelines Document Section 3 function review. Julia-first; Python and rpkg synced April 14-15, 2026.

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

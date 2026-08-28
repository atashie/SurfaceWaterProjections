# Cross-Language Implementation Status

> **STATUS (August 2026): the August 2026 port campaign implemented all six
> previously Julia-only features (Pettitt changepoint fields, the 20-value stats
> floor, the annual-values collector, the b=1 recession alpha, 14 snow metrics,
> 10 drought metrics) in BOTH Python and rpkg, and BOTH are now VALIDATED against
> canonical Julia at full scale — Python 2026-08-26, rpkg 2026-08-27.**
> Everything below the "Historical record" divider describes the APRIL 2026
> state, when the ports reproduced only a 623-signature-column subset. It is
> retained for provenance; do not read it as current.
> Campaign plan and phase-by-phase record:
> `docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md` and
> CHANGELOG → August 2026.

## Current status (August 2026)

"Implemented" = the feature exists in that language and is covered by its unit
suite. "Validated" is a separate, stronger claim — a full-scale benchmark diffed
against canonical Julia under the strict schema + value gates — and is recorded
per language below.

| | Julia | Python | rpkg |
|---|---|---|---|
| Role | **Canonical** | Full port, **validated** | Full port, **validated** |
| Summary columns (WY 1993–2025 @ 60 %) | 1,653 | **1,653 (gated)** | **1,653 (gated)** |
| Pettitt changepoint fields | yes | yes | yes |
| 20-value stats floor | yes | yes | yes |
| Annual-values collector + parquet | yes | yes | yes |
| b = 1 recession alpha | yes | yes | yes |
| Snow metrics (14) | yes | yes | yes |
| Drought metrics (10 + 5 scalars) | yes | yes | yes |

### Python vs Julia — feature-complete validation (2026-08-26)

Strict schema gate **PASS** with no waivers: identical column set (1,653) and
gage set (6,678). Over 1,620 shared signature columns:

| Tier | Count |
|---|---|
| Perfect (R² ≥ 0.999) | **1,615 (99.7 %)** |
| Good (0.99–0.999) | 5 |
| Poor / Low / Very Low / Extremely Low | **0 / 0 / 0 / 0** |

Mean R² **0.999988**, min R² **0.997935**. Annual parquet: identical
100-signature sets, 0 duplicate keys, **0 NA-pattern mismatches**, 18,898,405 of
18,898,406 rows shared.

**No divergence class remains.** The two residuals are characterised, not
merely tolerated:
1. **Pettitt changepoint ties.** `cp_year` agrees exactly on 597,505/597,527
   cells (99.9963 %); where a rank tie flips the changepoint year, that base's
   segment split moves and its segment p-values follow. This accounts for all
   5 Good columns (3 bases: FDC90th, Q90, Qsum) and, **at a 1e-6 tolerance**,
   for 267 of 18,217,552 finite annual pairs (0.0015 %) — every one a discrete
   threshold-crossing metric (`D*_day`, `D25_to_D75`, `TQmean`, one
   `dur_low_pulses_all`), where a last-bit tie moves a whole day.

   *Always quote the tolerance with these counts.* At 1e-9 the count is 465:
   the extra 198 rows are `qp_slope_sd` (116, max 1.3e-07) and
   `elasticity_annual` (82, max 5.5e-08), continuous quantities differing at
   ordinary cross-library floating-point noise rather than tie flips. The two
   figures are one result read at two thresholds.
2. **Signed-zero `unique` semantics** — 1 annual row in 18.9 M. Julia's
   `unique` uses `isequal`, under which `-0.0` and `+0.0` are DISTINCT values;
   numpy and R use `==`, under which they are equal. So a series containing
   both a negative and a positive zero yields one more distinct value in Julia
   than in the ports.

   **Direction of fix is Julia-side** (the canonical implementation is the
   outlier here, and `==` is the semantics the hydrology intends — a negative
   zero flow is not a different flow from zero). **Decision (user, 2026-08-26):
   note it as a future cleanup; it does NOT block the port campaign and does
   NOT justify regenerating any delivered product** — the measured cost is one
   annual row in 18.9 million. Apply the `==` semantics when this code is next
   touched for another reason, and re-measure rather than assuming the count
   stays at one. Tracked in CHANGELOG → `[Unreleased]` → Planned.

### rpkg vs Julia — feature-complete validation (2026-08-27)

131.9 min, **6,678 gages, 0 errored, 1,653 columns**. All four acceptance gates pass:

| Gate | Result |
|---|---|
| Strict schema equality | **PASS, no waivers** — identical column set (1,653) and gage set (6,678) |
| Swallowed family failures | **PASS** — 0 per-gage failures across a 997,194-line log |
| Cross-language annual parquet | **PASS** with two named residuals (below) |
| Identity R² (diagnostic) | **1,601 Perfect (98.8 %) / 10 Good / 9 Poor / 0 below 0.95**; mean 0.999843, min 0.971051 |

Annual parquet: identical 100-signature sets, 0 duplicate keys, **0 NA-pattern
mismatches**, 18,898,405 of 18,898,406 rows shared, and **518 of 18,217,552 finite
pairs (0.0028 %) over 1e-6**.

**The 19 non-Perfect columns are entirely the pre-existing irreducible classes** —
11 FDC90th + 3 FDCmid + 2 Q90 (near-zero-tail OLS on `log10(Q + 1e-10)`, and
Pettitt tie flips downstream of it) and 3 recession Spearman p-values (exact
permutation vs t-approximation for small n). No drought, storage, snow, BFI or
changepoint family appears.

Two residuals are waived by name at gate time rather than hidden:
1. **1 annual row in 18.9 M** — the signed-zero `unique` divergence in the storage
   distinct-value guard (see the accepted-divergence note in CHANGELOG).
2. **518 annual values** — 291 discrete threshold-crossing metrics (`D*_day`,
   `D25_to_D75`, `TQmean`, `dur_low_pulses_*`), where a last-bit tie moves a whole
   day; and ~199 recession-family rows traceable to **3 gages of 6,507** whose
   `recession_alpha_point_cloud_linear_reservoir` differs (max 4.9e-3) through
   tie-sensitive recession-event identification.

#### What the gates caught that nothing else could

The first rpkg benchmark passed its unit suite and still failed three of four
gates. Every finding was a real defect:

- `calculate_snow_metrics()` raised on a legitimately **zero-row** `snow_data`
  frame; the orchestrator swallowed it and 4 gages silently lost all 224 snow
  columns. The same guard was missing in flashiness, timing and both baseflow
  functions.
- `load_gages_ii_interference()` **intersected** the CONUS and AKHIPR column sets,
  dropping 4 metadata columns for the ~9,000 CONUS gages that have them.
- flashiness, both BFIs and storage **pre-allocated a row per year**, exporting
  1,039 annual rows Julia omits; storage also lacked Julia's ≥10-**distinct**-Q guard.
- **`smooth_daily_flow` used R's `mean()`.** Its long-double accumulator plus
  correction pass returns a different last bit than Julia's sequential `s += v`
  loop — 4,513 of 10,227 smoothed values differed by ≤ 3.6e-15. Harmless alone,
  but the drought thresholds are percentiles OF that series, so a threshold landing
  on a flow plateau flipped every plateau day at once through the strict `<`
  (gage 01589795 WY2002: Julia 116 days, rpkg 60). An explicit sequential double
  sum makes the series **bit-identical across all 10,227 days** — and is
  marginally faster. This alone was ~11,578 of the 12,096 annual mismatches.
  numpy's `mean` matches Julia's order for windows this short, which is why the
  Python port never showed it.

### Convention decisions taken during the campaign

- **Mann-Kendall p-value method.** scipy's `kendalltau` selects on TIES (not
  sample size): exact when untied, asymptotic WITHOUT continuity correction
  when tied. Julia and R both use the continuity-corrected normal
  approximation, and `Kendall::MannKendall` reproduces the Julia formula
  exactly — so Python was the sole outlier. Python was changed to the canonical
  formula (2026-08-26), which cut significance-call disagreement from 0.24 % to
  0.0009 % (main path) and 0.45 % to 0.0009 % (Pettitt segments).
- **Sparse-family row emission.** flashiness and both BFI families OMIT
  non-computable years rather than emitting NaN rows, matching Julia's
  `total_Q <= 0 -> continue`.
- **Canonical metric names throughout.** flashiness and FDC now use
  `flashinessRB` / `FDCall` / `FDC90th` / `FDCmid` internally in all three
  languages; the ports' old placeholder-then-rename pattern would have written
  the wrong signature names into the annual parquet.

### Defects found by the campaign's gates

The strict schema gate, the floor-transition check and the cross-language
annual comparator each caught real defects that the older intersection-based
comparison could not see: the ports' flashiness `== 0` guard (vs canonical
`<= 0`), Python's `qp_seasonality`/`storage` silently dropping their new
changepoint keys, placeholder metric names in BOTH ports, NaN-vs-absent row
emission, rpkg's 12 double-prefixed QA-flag columns, and a mask-tool float
parser. Full list: CHANGELOG → August 2026.

Two later defects came from a different kind of check — running the actual
pipeline end to end on real data rather than testing modules in isolation:

- **rpkg's benchmark runner did not pass the new arguments** (no `changepoint`,
  `min_values_for_stats`, `collector`, or `snow_data`; SWE discarded at the
  climate join), and rpkg's config carried no `changepoint` or `stats_floor`
  section at all. A 19-hour run was discarded because it provably could not have
  produced the ported columns.
- **rpkg's drought family was silently all-NA on every real gage.** rpkg's
  `preprocess_daily_data()` renames `date` → `Date` internally, while the ported
  drought module checks the Julia/Python lowercase name; the missing-column
  branch returned an all-NA family behind a warning. The unit tests could not
  catch it because they build frames with `date` directly and never traverse
  the preprocessor → drought path. Fixed by accepting either spelling, plus a
  new end-to-end test that goes through `preprocess_daily_data()` and the
  orchestrator.

The shared lesson: **module-level unit tests cannot prove a port is wired in.**
Both defects presented as clean, warning-only degradation — exactly the failure
mode the campaign's per-gage schema assertions are meant to convert into a hard
error.

---

# Historical record (April 2026 and earlier)

Everything below predates the August 2026 port campaign.

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
- **Section 3 (April 2026, synced to Python/rpkg April 14-15):**
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

**Python**: Production-ready. 624 columns (all Julia-equivalent signatures including recession-parameterized BFI, ported April 25). Major fixes applied across 6 alignment rounds:
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

Note: These are three-way results (R vs Python vs Julia on the same 551 pre-Section 3 columns). All three synced implementations (Julia, Python, rpkg) now produce 624 columns including Section 3 additions (D1/D99, n_recession_events, elasticity_annual, negative_ann, diagnostics) and recession-parameterized BFI (alpha_linear, BFI_Eckhardt_param, BFI_LyneHollick_param + scalar). Monolithic R is not synced with these additions.

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
7. Trend completeness: ≥60% non-NA annual values overall (lowered from 80% in July 2026) + ≥80% in first/last decade required for trend stats (slopes, p-values); mean/median always computed. Recession and elasticity exempted.

Config constants are imported from shared `config/signatures_config.json` instead of hardcoded. Legacy filtering (`use_legacy_filtering: true`) preserves the old 95%-non-NA-days rule.

Output metadata columns are aligned: `basin_area`, `start_water_year`, `end_water_year`, `num_water_years`.

Julia/Python produce 7,313-7,369 qualifying gages (vs monolithic R's 5,707). The extra gages lack Daymet climate coverage — monolithic R only iterates gages with Daymet data, while Julia/Python process all gages and leave climate signatures as NA when Daymet is unavailable.

## Output Column Notes

### Human Interference Metadata (13 columns)

All three languages include GAGES-II human interference metadata in their output: `NDAMS_2009`, `MAJ_DDENS_2009`, `STOR_NID_2009`, `IMPNLCD06`, `DEVNLCD06`, `FRESHW_WITHDRAWAL`, `HYDRO_DISTURB_INDX`, `CLASS`, `RHBN`, `REGULATED`, `human_interference_class`, plus `area_normalized` and `gage_type`. Julia loads GAGES-II metadata directly and Canadian HYDAT metadata from `metadata/canadian_hydat_interference.csv` (pre-exported from tidyhydat via `R/export_hydat_metadata.R`; 8,012 stations with RHBN and REGULATED status). Monolithic R loads both GAGES-II (USGS gages) and HYDAT (Canadian gages) metadata via `tidyhydat`. Python and rpkg load GAGES-II metadata directly and fill RHBN/REGULATED from `metadata/canadian_hydat_interference.csv` (added in the August 2026 runner rebuild — the earlier "Canadian HYDAT integration pending" state is superseded).

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

- All three synced implementations (Julia, Python, rpkg) produce 624 signature columns (551 shared three-way + 43 Section 3 + 25 recession-parameterized BFI + 5 scalars/diagnostics). Recession-parameterized BFI ported to Python and rpkg April 25, 2026.
- **Python vs Julia (624 cols, April 27)**: 615 Perfect, 5 Good, 3 Poor. Min R² = 0.984. All 25 new recession-parameterized BFI columns are Perfect (R² = 1.000).
- **rpkg vs Julia (624 cols, April 27)**: 614 Perfect, 4 Good, 5 Poor. Min R² = 0.969. 24 of 25 new recession-parameterized BFI columns are Perfect (1 poor: `alpha_linear_spearman_pval` — Spearman p-value library difference, R² = 0.971).
- Python-Julia: 3 columns below 0.99 (all recession seasonality minimum — irreducible sinusoidal fit sensitivity)
- rpkg-Julia: 5 columns below 0.99 (3 recession Spearman p-values + 2 FDC90th p-values — irreducible library differences)
- All 7,313 gages matched across all three implementations
- Julia is ~10x faster than Python (~15 min vs ~150 min), ~16x faster than rpkg (~15 min vs ~245 min with contention)
- Median R² = 1.000 for both Python-Julia and rpkg-Julia

# Changelog

All notable changes to the Streamflow Signatures project.
For full historical detail (Dec 2025 – April 2026), see [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

> **Convention** — the most recent month(s) appear here as a condensed summary with pointers; full per-change detail lives in the archive. When a month's work is complete, condense it here and move the detail to the archive. File-level change lists belong in `git log`, not here; analysis and benchmark tables belong in the canonical docs (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`) and are linked rather than re-hosted.

## [Unreleased]

### Planned
- Run `smoke_test.jl` (10-gage subset) for the drought family on a machine with the
  `D:\processedOuts_feb2026\` inputs — the unit suite is green, the smoke test is not
  yet exercised
- Port annual-values collector, b=1 recession alpha, snow metrics, AND drought metrics
  to Python and rpkg (after the Julia benchmark re-run validates all four)
- Regenerate both standard products (WY 1993–2025 and WY 1980–2025 @ 60 %) to carry the
  drought family — decide whether to batch with the ports
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)
- Port data ingestion utilities to Julia (long-term — currently R-only via dataRetrieval/tidyhydat)
- Generate Julia golden outputs (624 cols, 7,313 gages)
- BFImax estimation via Collischonn & Fan (2013) backward filter — would give BFI_Eckhardt_param per-gage BFImax instead of fixed 0.8, improving discriminating power (currently range [0.47–0.80] due to BFImax saturation)

### Known Issues (discovered 2026-07-14, Codex review — not yet fixed)
- **LOW (raw/legacy path only, discovered 2026-07-27 by the drought smoke run) —
  `calculate_negative_days` crashes on `missing` Q and silently drops its 8 columns.**
  `julia/src/pulses.jl:323` applies `:Q => (q -> sum(x -> !isnan(x) && x < 0, q))`
  straight to a `Union{Missing,Float64}` column; `!isnan(missing)` is `missing`, so
  `missing && …` throws `TypeError: non-boolean (Missing) used in boolean context`. The
  orchestrator's per-signature `try/catch` turns this into a warning, so the gage just
  loses all 8 `negative_ann` columns (smoke gage 01073000 emitted 808 signature columns
  vs 816 for the other nine). Every other signature routes Q through `coalesce_q`; this
  one does not. **Production output is unaffected** while `use_legacy_filtering: false`,
  because `preprocess_daily_data()` emits Float64 Q with NaN rather than `missing` — it
  only bites callers passing raw frames (smoke test, direct API use). Fix is one line
  (`coalesce_q(df.Q)` before the group-by, matching the rest of the codebase); left
  unapplied because it is outside the drought work's scope. Python/rpkg ports likely
  mirror the pattern — verify when touched.
- **LOW (legacy shim only, discovered 2026-07-21) — 6 pre-existing failures in the
  legacy R NA-handling test suite** (`R/tests/test_na_handling.R`, run per its
  documented usage after sourcing `config.R` + `R/helperFunctions.R`): grid
  normalization ("Missing date rows filled"), "3-day gap: year accepted",
  constant-SD flag, and raw/residual NA diagnostic counts. Verified present at
  HEAD before the 2026-07-21 trend-gate work (which only touched — and fixed —
  the vacuous trend-completeness case). Indicates drift between the deprecated
  `R/helperFunctions.R` shim and the evolved test expectations; rpkg is the
  active R implementation and its testthat suite is unaffected. Triage when the
  legacy shim is next touched, or retire the legacy suite with it.
- **MEDIUM (documented limitation, by design) — 37 Canadian gages in the signature
  output carry raw m³/s units (`area_normalized = FALSE`)**: HYDAT publishes NO
  drainage area (neither `DRAINAGE_AREA_GROSS` nor `DRAINAGE_AREA_EFFECT` — verified
  directly against `Hydat.sqlite3`) for 73 successfully-processed Canadian stations;
  37 survive the 20-year filter into the Julia canonical output (7,313 gages).
  Station names show most are not natural watersheds: ~40 irrigation/diversion canals
  + ditches, ~15 dam/powerhouse outflows, several huge-river channel splits and lake
  outlets (St. Lawrence, Mackenzie, Nelson, Lake of the Woods); only ~8 look like
  natural streams. 62/73 are `REGULATED = TRUE`.
  **DECISION (user, July 2026): keep these gages with raw m³/s flow — NO area
  backfill** (HydroBasins `UP_AREA` was assessed on 1,383 validation gages: accurate
  only for main-stem dam outflows, wrong for canals and channel splits, unusable
  <100 km²). Q-to-PPT signatures are now structurally gated for these gages (see the
  July 2026 entry below). **Remaining limitation**: unit-carrying Q-only signatures
  (Q volumes, percentiles, Q95_Q10, log_a) stay in m³/s for these 37 rows —
  incomparable with mm/day gages (Qann_mean up to 3.18M for the St. Lawrence).
  **Flag gap**: `flagged_for_qann_range` catches only 27/37 — 10 small canals/creeks
  land inside [0, 2000] unflagged; downstream users must filter on
  `area_normalized == TRUE` before any cross-gage comparison of unit-carrying
  signatures. See docs/DEVELOPMENT.md → Canadian HYDAT → Missing drainage areas.
- **MEDIUM — seasonal runoff ratios ignore seasonal completeness flags**:
  `julia/src/runoff_ratios.jl:113-118` looks up flags named
  `winter_complete/spring_complete/summer_complete/fall_complete`, but the preprocessor
  emits `win_complete/spr_complete/sum_complete/fal_complete` (`julia/src/io.jl:462`;
  cf. correct usage in `flow_volumes.jl:136-137`). The existence check at
  `runoff_ratios.jl:120` silently fails, so seasonal runoff ratios are never NaN'd for
  incomplete seasons — the guidelines' 80% seasonal completeness rule is not applied to
  this category. Python and rpkg likely mirror the bug (synced ports — verify). Fix
  changes outputs → schedule with the next behavior-changing release (Julia first, then
  ports, benchmark re-run + golden comparison).

### Guidelines Document TODOs
<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

Synced 2026-07-21 — first doc revision since 2026-04-15. Behavior-changing items:

- [x] **Trend completeness gate: overall series 80% → 60%** (source: NA Handling) —
  **IMPLEMENTED 2026-07-21** (see July 2026 entry below). Config-only change picked
  up by all languages; takes numerical effect at the next benchmark/production run.
- [x] **Decade gate 60% vs 80% — RESOLVED (user, 2026-07-22): 80% first/last
  decade + 60% overall CONFIRMED.** Matches the shipped config
  (`decade_min_fraction: 0.80`, `min_fraction: 0.60`), the guidelines doc, and
  manuscript §2.2.3 — no code change needed. The linked snow record-anchored
  gate therefore also runs at 0.80 decades.
- [ ] **Recession fit-quality flag: R² < 0.8** (source: analyze_recession_parameters).
  "To control the quality of fitted parameters, calculate recession fits and create a
  flag for any R2 < 0.8." Not currently implemented. Needs design: per-event vs
  per-gage flag, and note that b=1 alpha fits are medians (no regression R²) — the
  free-fit b regressions are the natural target. Clarify scope with domain experts.
- [ ] **avg_storage omitted from major analyses (4/23/26)** (source:
  calculate_average_storage). Doc header now says "OMITTING VARIABLE FROM MAJOR
  ANALYSES"; extensive redesign notes added (3 options incl. GLEAM ET water balance;
  Erin's dormant-season recession storage method needing PET, dormant-season
  definition, initial-storage assumption; open questions on snow). Decide: gate/flag
  `avg_storage` in outputs vs leave computed + documented as excluded downstream.
- [ ] **NA-handling wording conflict — item (i) "flag not remove"** (source: NA
  Handling). New sentence: "Items i, iii, and iv are set in the config to flag (not
  remove) by default." Items iii (negative Q) and iv (constant SD) match the current
  config-driven flag-only behavior, but item i (>3 consecutive days of NAs) currently
  REJECTS the year in `preprocess_daily_data()` (not config-toggleable). Clarify with
  domain experts whether item i should become a config-driven flag.

Documentation-only / already-implemented items from this sync (no action):
- FDC glossary section added (FDCall/FDC90th/FDCmid — already implemented).
- Recession b free-fit + a via b=1 linear aquifer now in the doc — matches the July
  2026 Julia implementation.
- Q-P seasonality: "backward-looking" 30-day rolling window confirmed — matches the
  implementation (trailing window `[end-29, end]`, mid-window month attribution);
  new caveat that the rolling approach may differ from Wrede et al. (2015).
- Elasticity: annual (t1−t0) second calculation marked "achieved as of 4/16"; the
  <30%-missing counterfactual marked "added as documentation not filter" — resolves
  the previously-pending 30% diagnostic item as documentation-only.
- Runoff ratio NA rule now states the implemented thresholds (annual PPT < 10mm,
  seasonal < 1mm → NA).
- Flow volume glossary corrected to "total" (was "mean") streamflow — matches
  implementation; flow-timing "labeled as julian" fix-note removed (resolved).

### Manuscript Reconciliation Log

Session-start reconciliation of the HISSS manuscript draft (Scientific Data,
submission target Nov 9 2026) against code + repo docs. Snapshot:
`docs/MANUSCRIPT_DRAFT.md`; workflow: CLAUDE.md → Session-Start Workflow → B.

**2026-07-21 — baseline established + initial reconciliation pass.** Snapshot
created and the manuscript added to the session-start review. Findings from the
first pass (manuscript §refs; direction of fix in brackets):

- **§2.2.3 asserts the 60% overall trend gate** — consistent with the July
  guidelines doc, but code still enforces 80% (`trend_completeness.min_fraction`).
  [Fix in code — same item as Guidelines TODO #1 above; the manuscript and
  guidelines agree, the code lags.]
- **§2.2.3 repeats "items i, iii, and iv … flag (not remove) by default"** —
  item i (>3 consecutive NA days) currently REJECTS the year in the preprocessor,
  not config-toggleable. [Clarify with domain experts — same as Guidelines TODO
  #4; whichever way it resolves, manuscript, guidelines, and code must end up
  aligned.]
- **§2.2.2 says statistics were computed for "a subset from 1993-2025
  wateryears"** — the July 2026 production run used an explicit WY 1993–2022
  window (user-directed end cap; see
  `docs/plans/production_rerun_1993_2022_60pct_plan.md`). [Flag for authors —
  or confirm a future 1993–2025 run is intended.]
  **[CORRECTION + RESOLVED 2026-07-22]**: the July 14 run was actually
  WY ≥ 1993 UNCAPPED — its plan's §10 execution record superseded the 1993–2022
  title (this log's original claim mis-read the plan). Moot either way: the
  **WY 1993–2025 standard run executed 2026-07-22 (Codex results-review GO)**,
  so the manuscript's "1993-2025" wording is now accurate as written.
- **§2.1.2 says flow was normalized "by watershed area as defined in the
  watershed boundary shapefile"** — implementation normalizes by the published
  drainage areas (GAGES-II metadata for USGS; HYDAT STATIONS
  `DRAINAGE_AREA_GROSS` for Canada), not shapefile polygon areas; and 37
  Canadian gages carry raw m³/s (`area_normalized = FALSE`, no published area).
  [Flag for authors — wording fix + the un-normalized-gage caveat belongs in
  Usage Notes/limitations.]
- **§2.2.2 metric-family list** duplicates "baseflow" and omits runoff ratios,
  Q-P seasonality, snow metrics, and pulse/timing details; stats list ("slope,
  Spearman's Rho, p-value, mean and median") omits the second slope (linear),
  Mann-Kendall tau/p-value, and the Pettitt changepoint fields. [Flag for
  authors — incomplete enumeration of the delivered columns.]
- **§2.1.3 says Daymet was aggregated to 6,041 basins** via gdptools/agg_gen —
  repo docs record ~6,087 sites in `daymet_1980_2023.parquet`, and the July
  rerun found snow values for 5,622 gages (4,533 US + 1,089 Canadian). [Verify
  counts with authors — the gdptools aggregation is the co-authors' upstream
  step, not in this repo.]
- **§1 describes the dataset as "at the HUC-8 scale"** — products are per-gage
  watersheds (GAGES-II / ECCC polygons + HydroBasins fallback), not HUC-8
  units. [Flag for authors.]
- Placeholders to watch on future syncs: methods summary paragraph, input-data
  table, §2.2.1 (empty), "n=xx", metadata file "(name)".

**2026-07-21 — user decisions on the initial findings.**

- **§2.2.3 (60% trend gate) → code updated.** `trend_completeness.min_fraction`
  0.80 → 0.60 implemented same day (see July 2026 entry). Manuscript, guidelines,
  and code now agree; takes numerical effect at the next run.
- **§2.2.2 (1993–2025 window) → outputs will catch up to the manuscript.** Two
  "standard" products planned: WY 1993–2025 @ 60% (first), WY 1980–2025 @ 60%
  (second, later) — see Planned. No manuscript edit needed for the window; §2.2.2
  should ultimately describe both standard outputs once they exist.
- **§2.1.2 (area normalization wording) → manuscript edit** (code is correct).
- Remaining findings → manuscript edits, queued below for Arik to apply in the
  Google Doc after this documentation pass + Codex review.

**Manuscript edits queued (to be applied in the Google Doc by Arik):**

1. **§2.1.2** — replace "normalized volumetric measurements by watershed area as
   defined in the watershed boundary shapefile (Sect. 2.1.1)" with the actual
   method: normalization uses agency-published drainage areas (GAGES-II metadata
   for USGS gages; HYDAT STATIONS `DRAINAGE_AREA_GROSS` for Canadian gages), not
   polygon areas from the boundary shapefiles. Add (here or in Usage Notes): 37
   Canadian gages have no published drainage area and are retained in raw m³/s
   with `area_normalized = FALSE`; their Q-to-PPT signatures are structurally NA,
   and downstream users must filter on `area_normalized == TRUE` before
   cross-gage comparison of unit-carrying signatures.
2. **§2.2.2** — fix the metric-family list: "baseflow" appears twice; runoff
   ratios, Q-P seasonality, and snow metrics are missing (and storage is
   intentionally omitted per the 4/23/26 guidelines decision — worth stating).
3. **§2.2.2** — complete the statistics list: per metric the output is 8
   statistics (Theil-Sen slope, linear-regression slope, Spearman's rho +
   p-value, Mann-Kendall tau + p-value, mean, median) plus 8 Pettitt changepoint
   fields — not just "slope, Spearman's Rho, p-value, mean and median".
4. **§2.2.2** — once the standard runs exist, describe the two standard output
   windows (user, 2026-07-22): WY 1993–2025, and "entire period of record"
   operationalized as WY 1980–2025 — both @ 60% qualifying fraction. The current
   sentence ("full period of record and a subset from 1993-2025") is close but
   should state the 1980 start of the entire-record product explicitly.
5. **§2.1.3** — verify "6,041 basins": the Daymet parquet used by this pipeline
   carries ~6,087 sites (repo docs); reconcile which count the gdptools
   aggregation produced.
6. **§1** — "at the HUC-8 scale" should be "at the gaged-watershed scale"
   (GAGES-II / WSC basin polygons, HydroBasins fallback), not HUC-8 units.
7. **§2.2.3** — "Items i, iii, and iv are set in the config to flag (not remove)
   by default" is accurate for iii (negative Q) and iv (constant SD) but NOT for
   i (>3 consecutive NA days), which unconditionally rejects the year in the
   preprocessor. Pending the domain-expert clarification (Guidelines TODO above),
   either the sentence drops item i or the code makes item i config-driven.
8. **§2.2.2 (new, 2026-07-27)** — add the **streamflow drought** family to the
   metric-family list: duration + deficit at five fixed percentile thresholds
   (2/5/10/20/30 %, U.S. Drought Monitor D4–D0), computed on 7-day-smoothed flow
   with the unbiased Weibull plotting position, citing Adelsperger et al. (in
   review) and Laaha et al. (2017). Two deviations from the source method must be
   stated: only the FIXED (whole-record) thresholds are implemented (the variable
   day-of-year thresholds are too uncertain at the low levels with 20–46 years of
   record, 2 % falling below the smallest plotting position outright), and aggregation
   is by WATER year rather than the paper's climate year
   (Apr–Mar), so droughts crossing Oct 1 are split across two annual values. Also
   worth a Usage Notes line: drought values are record-dependent (thresholds come
   from each product's own window) and `drought_deficit_*` is unit-carrying.

---

## [July 2026]

### New: streamflow drought signature family (10 metrics, Adelsperger et al. — Julia)
Per-water-year **drought duration** (days below threshold) and **drought deficit**
(summed departures below threshold, mm) at five fixed severity levels — 2/5/10/20/30 %
magnitude percentiles, mirroring the operational U.S. Drought Monitor D4→D0 ladder —
through the standard `generate_stats()` path (8 stats + 8 Pettitt fields each = 160
columns, plus 5 per-gage `drought_threshold_fixed_p{n}` scalar diagnostics →
**+165 columns, 1,488 → 1,653**; 90 → 100 time-series signature bases in the annual
parquet, where the scalars do not appear; annual values exported via the collector).
New module `julia/src/drought.jl` (`calculate_drought_metrics`); full spec + decision
record: `docs/plans/2026-07-27-drought-signatures-plan.md`. Requested by the user
(2026-07-27) from Adelsperger et al. (in review), "A novel severity-based approach for
assessing streamflow drought characteristics and drivers", with four in-session scope
decisions.

- **Method**: 7-day **centered** smoothing of daily Q applied to the CONTINUOUS
  date-indexed series (no artificial Oct 1 discontinuity), within maximal runs of
  consecutive dates so the window never averages across a rejected-year gap (shrinking
  at run edges; < 4 valid days → NaN). Thresholds are whole-record percentiles of the
  smoothed values via the **unbiased Weibull plotting position** `i/(n+1)`
  (Hyndman-Fan definition 6 = `quantile(x, p; alpha=0, beta=0)`) — deliberately NOT the
  type-7 default used elsewhere, which differs most in exactly this low tail.
  Comparison is strict (`<`); ≥ 10 qualifying years required for a threshold.
- **Scope decisions (user, 2026-07-27)**: **fixed thresholds only** — the paper's
  variable day-of-year method is NOT implemented because its per-day sample is one value
  per year, so at 20–46 years of record the low levels carry very large sampling
  uncertainty and the 2 % level falls below the smallest Weibull plotting position
  `1/(n+1)` outright (2 % needs ≥ 49 years, 5 % needs ≥ 19); type 6 would clamp to the
  sample minimum there, and this project's `below_plotting_range_policy: "na"` refuses to
  — an estimability policy layered on type 6, not a property of it; all five levels;
  **water-year** aggregation (documented deviation from the paper's Apr–Mar climate
  year, which splits droughts crossing Oct 1 — day-level indicators and record totals
  are identical either way, annual series and trends are not); NaN rather than
  clamp-to-minimum below the plotting range (unreachable for the fixed method, kept as
  a defensive invariant); explicit `_fixed_` infix in column names so a future variable
  method can be added non-breakingly.
- **Record-dependent**: thresholds come from the run's own window, matching the paper's
  per-analysis-period thresholds — so drought values join `*_all` pulses, elasticity,
  and the parameterized BFI on the "valid within a product, never compared across the
  WY 1993–2025 and WY 1980–2025 products" list. `drought_deficit_*` is also
  unit-carrying (mm only where `area_normalized = TRUE`); the durations are
  scale-invariant and valid everywhere.
- **Gates**: trend completeness and the 20-value stats floor both apply (NOT exempt).
  No snow-style record-anchored gate is needed — a drought-free year emits a valid 0,
  which is precisely the signal a trend should see. **NO existing column changes at all**:
  `flagged_for_high_na` was expected to shift (its denominator counts all signature
  fields — the April 2026 Pettitt / July 2026 snow precedent), but direct per-gage
  recomputation found **zero crossings** (1,224 flagged before and after). Do not
  generalize this to other windows — the closest gage sits only ~6.1e-5 from the 30 %
  threshold, so the outcome is window-specific.
- **Tests**: `julia/test/test_drought_metrics.jl` — hand-derived Weibull quantiles and
  the below-plotting-range guard; exact smoothing values incl. the never-blend-across-a-
  gap and duplicate-date cases; and two record invariants that pin the threshold
  definition AND the day → water-year attribution simultaneously (summed duration ==
  count of pooled values below the threshold == `floor((n+1)p)` for a distinct-valued
  sample; summed deficit == the pooled departure sum). Plus level monotonicity, wet-year
  valid zeros, the intermittent zero-threshold/strict-`<` case, threshold and stats
  floors, schema/collector/zero-warning contracts, and orchestrator + annual-export
  wiring. `smoke_test.jl` gained presence, range, and monotonicity checks.
- **Suite green**: `julia/test/runtests.jl` → **2,042 assertions, 0 failures** (drought
  file: 670). Julia 1.12.6 was installed locally (juliaup) to run it. The expectations
  had first been derived against a purpose-built pure-Python mirror of `drought.jl`;
  the Julia run then confirmed them. **One real registration gap the suite caught**:
  `EXPECTED_DENSE_SIGNATURES` in `test_annual_collector.jl` pins the exact set of dense
  annual series and must list every new dense signature — the 10 drought bases were
  added (that test asserts set EQUALITY, so a new family fails it until registered; now
  called out in the CLAUDE.md / DEVELOPMENT.md "Adding New Signatures" checklists).
- **Smoke test PASSED on real data** (10 gages × 45–46 WYs, feb2026 inputs): all 10
  bases present, ranges valid, level monotonicity holds. Mean durations came out at
  7.296 / 36.499 / 109.542 d/yr for p2 / p10 / p30 against the `p × 365.25` construction
  target of 7.30 / 36.52 / 109.58 — a **weak, near-circular consistency check** (the
  threshold is a percentile of the series it is then counted against), NOT a proof that
  the smoothing or plotting position is right; it catches gross mismatches only.
  `smoke_test.jl` gained ENV path overrides (`STREAMFLOW_DATA_PATH` /
  `STREAMFLOW_CLIMATE_PATH`, the same names the benchmark runner already uses) so it can
  run wherever the inputs are mounted instead of only the Windows `D:\` drive.
- **Codex adversarial review (2026-07-27, codex-cli 0.145.0, read-only): GO-WITH-FIXES.**
  No CRITICAL; 6 MAJOR + 5 MINOR, all resolved — full table in the plan doc §16. The
  substantive ones: (MAJOR-1) the tests' record invariant proved **conservation, not
  attribution** — a uniform water-year shift would have passed it, so a hand-derived
  boundary fixture (a low-flow block straddling Sep 30 → Oct 1 splitting exactly 6 days /
  27.0 mm into WY 2000 and 7 days / 36.0 mm into WY 2001) plus an independent per-year
  recompute were added; (MAJOR-2) the `p × 365.25` check is near-circular and its
  description as "threshold correctness confirmed" was withdrawn; (MAJOR-3) the
  redundancy ratio was statistically invalid (`mean(n) × mean(d) ≠ mean(n × d)`) and the
  non-redundancy claim is withdrawn pending an annual-series comparison at benchmark
  scale; (MAJOR-5) the column count was wrong — **+165 → 1,653**, not +160 → 1,648;
  (MAJOR-6) "unresolvable" overstated the variable-threshold problem. Two code fixes:
  rows with unparseable dates are no longer smoothed with their neighbours, and numeric
  config values (even "centered" windows, `min_valid > window`, out-of-range or unsorted
  percentiles) now fail fast. Post-fix suite: **2,700 assertions, 0 failures** (drought
  670 → 1,328).
- **Benchmark run (2026-07-28, WY 1993–2025 @ 60 %, thumbdrive inputs on the M1):**
  6.05 min, **6,678 gages × 1,653 columns**, annual parquet 18,898,406 rows / 100
  signatures. Gage set identical to standard product #1, confirming the family changes
  no gage qualification. Outputs in their own folder per the one-folder convention:
  `/Volumes/Untitled/processedOuts_drought_28jul2026/` (wrapper
  `docs/benchmarks/run_julia_benchmark_drought_1993_2025_60pct.jl`).
- **REDUNDANCY MEASURED — `drought_duration_fixed_p10` is largely redundant** with the
  existing pulse pair (new tool `docs/benchmarks/analyze_drought_redundancy.jl`, run on
  the ANNUAL series, where `n × mean-duration` reconstructs that year's sub-threshold day
  count — mathematically exact, subject to floating representation). Over 200,834
  gage-years: **r = 0.979**, ρ = 0.971; and judged against the interannual signal rather
  than the series mean (the quantity a trend consumes), **within-gage median r = 0.994**
  (p10 0.971) with disagreement ≈ **11.7 % of each gage's own duration SD** (p90 25.3 %).
  Only p10 collapses: the other four levels measure at r = 0.712 / 0.902 / 0.846 / 0.731
  (p2/p5/p20/p30) against the same pulse pair. Redundancy is an aggregate statement —
  32.5 % of gage-years agree exactly, max disagreement 318 days; on intermittent gages
  99.87 % agree exactly (13 of 9,652 differ, max 8 days). The family's non-redundant
  content is therefore (a) `drought_deficit_*`, the only magnitude-weighted low-flow
  measure in the output, and (b) the four non-p10 severity levels. SIGNATURES.md and the
  claude-skill steer users accordingly. **DECISION (user, 2026-07-28): `drought_duration_
  fixed_p10` is KEPT** — cross-family redundancy is not grounds for removal (err toward
  abundance; the same quantity via two documented methods has independent value, and the
  severity ladder stays complete at all five USDM rungs). The overlap is a documentation
  caveat, not a defect: don't present it and the pulse pair as independent evidence.
- **✅ ADDITIVITY GATE: PASS** (closes Codex MAJOR-4). The rigorous test is a same-machine,
  same-Julia baseline with the drought family disabled (config copy without the `drought`
  section via `STREAMFLOW_CONFIG`), diffed with EXACT equality against the drought-enabled
  run: 165 columns added (expected 165), no column dropped, gage sets identical, **all
  1,487 shared columns bitwise unchanged**, every added column populated. A
  **cross-machine** diff against the delivered product #1 canNOT serve as this gate: 431
  columns differ as pure FP noise (≤ 3.2e-06) and 66 differ materially on 1–38 of 6,678
  rows — all of them discretely FP-sensitive statistics (`TQmean*`, where one flipped day
  is 100/365 = 0.274 pp; `FDC90th*`, OLS on `log10(Q + 1e-10)` in the near-zero tail,
  already documented here as FP-fragile; and rank/p-value/Pettitt fields where a tie
  flips). None attributable to the drought family.
- **New validation tooling**: `docs/benchmarks/check_additivity.jl` — proves a run ADDED
  columns without changing pre-existing ones (column set, gage set, per-gage value
  identity with NaN==NaN, added-column population, `--allow-shift`, `--expect-added`,
  and a `--tol` mode for cross-machine comparisons). Required by the new CLAUDE.md /
  DEVELOPMENT.md checklist step.
- **Fixed (HIGH, config plumbing): GAGES-II directory is now resolved at RUNTIME, not
  precompile time.** Symptom: a benchmark silently emits 1,642 instead of 1,653 columns,
  dropping the 11 GAGES-II/HYDAT human-interference columns with only a warning, because
  `metadata.gages_ii_dir` still points at the Windows `D:/gagesMetadata`. First fix
  attempt added an ENV-derived **constant** `CFG_GAGES_II_DIR` — which does not work: a
  const is baked at precompile time, so a wrapper setting
  `ENV["STREAMFLOW_GAGES_II_DIR"]` at runtime is silently ignored (the same
  precompilation trap documented in `julia/src/config.jl` and DEVELOPMENT.md; it bit
  three times in one day). Now a **function** `gages_ii_dir()` reads ENV at call time and
  is the default argument of `load_gages_ii_interference`; `CFG_GAGES_II_DIR` remains as
  the JSON-only value for compatibility. Verified by probe: in a session where the ENV var
  is set, the constant still returns `D:/gagesMetadata` while the accessor returns the
  mounted path. **Lesson: any ENV override consumed through a `const` is unreliable** —
  this is why the runner already re-reads the window/fraction ENV vars at runtime.
- **Codex review #2 — RESULTS/ANALYSES (2026-07-28): GO-WITH-FIXES, product approved for
  promotion** once the doc errors were corrected; all 4 MAJOR + 7 MINOR resolved (table in
  plan §17). Codex independently re-verified the config equivalence, gage sets, and **all
  16,890,066 shared annual values** (zero mismatches). The findings that mattered:
  (MAJOR-1) my claim of "exactly zero disagreement on intermittent gage-years" was
  **FALSE** — r = 0.981 alongside ρ = 0.492 cannot coexist with identical series; the true
  figure is 99.87 % exact-equal, 13 of 9,652 differing, max 8 days, and the analysis script
  now reports exact-equal %, nonzero count and MAX |diff| so quantiles can't hide it again;
  (MAJOR-2) the redundancy argument scaled a median difference against the series *mean*
  when a trend consumes interannual variation — re-measured within-gage (median r 0.994,
  RMSE/SD 0.117), conclusion survives; (MAJOR-3) `check_additivity.jl`'s population gate
  passed a column with ONE finite value of 6,678, so it could not detect the failure mode
  it existed for — now `--min-finite-frac` (default 0.5) with a by-suffix breakdown, re-run
  at 0.75 = PASS; (MAJOR-4) runs lacked provenance.
- **New: provenance block in the benchmark timing JSON** (`run_julia_benchmark.jl`) —
  resolved paths + size/mtime for every input, SHA-256 for files < 50 MB (config,
  metadata; `STREAMFLOW_HASH_INPUTS=1` also hashes the multi-GB parquets), git revision +
  dirty flag, Julia version, hostname, and the experiment ENV overrides. Additive to the
  timing JSON only — no signature output changes. **The 28 Jul product predates it**; a
  6-min re-run would capture it.
- Also documented from this round: `flagged_for_high_na` stability is window-specific (the
  closest gage sits 6.1e-5 from the 30 % threshold — do not generalize); the ~20 % NaN rate
  in some drought columns is two separate mechanisms (495 gages / 7.4 % lose trend stats to
  the completeness gate — the *same* set that loses `Qann_senn_slp`, so not drought-specific
  — plus ~300 losing rank stats to constant series); and an audit confirming no past
  config-variant result is invalidated by the precompilation gotcha (DEVELOPMENT.md).
- Julia only; Python/rpkg ports deferred (now a four-feature port queue).
  `validate_production_run.py` signature-count gate updated 90 → 100.

### New: per-watershed Annual NLCD product (CONUS land cover + impervious, non-signature)
A CONUS-only US companion to the continental MODIS LULC product, built on **USGS/MRLC Annual NLCD
Collection 1 (C1V2, 30 m, 1985–2025)**: per-watershed annual land-cover % (16 classes) + basin-mean
fractional impervious surface for **6,119 CONUS gages × 41 years = 250,879 rows**, joining to the
signatures by leading-zero-safe `gage_id`/`canon_id`. Python (metadata/ingestion — like HydroATLAS
and the MODIS EO layers, not a cross-language signature). New modules under `EO_data_processing/`;
full design + review record: [EO_data_processing/README_NLCD.md](EO_data_processing/README_NLCD.md).

- **Source/access** — Annual NLCD COGs on `s3://usgs-landcover` (us-west-2, requester-pays); native
  CRS is **WGS84-Albers (no EPSG code — read from `src.crs`, not the legacy EPSG:5070)**; fill=250.
  Source mosaics staged once to `s3://…/streamflow/temp_lulc_conus/` (82 files, 93.4 GB, us-east-2
  in-region; **temporary — delete after QA**). Fractional-impervious carries a documented
  out-of-range/underflow bug — guarded by a [0,100] clip (0 occurrences in C1V2).
- **Pipeline** (`eo_processing/nlcd_pipeline.py`) — per-year-parallel download→extract→delete
  (scratch on `/tmp`, NOT the 2.6 GB home volume), coverage-weighted exactextract (basins reprojected
  to the granule CRS; never resample categorical); publish-integrity gates (per-year validation,
  atomic writes, final table written only if every requested year passes, else nonzero exit — no
  silent partial publish). ~10 h for 41 yr on 4 cores. Class %s sum to exactly 100 on covered rows.
- **Finalize** (`eo_processing/nlcd_finalize.py`) — drops the 45 out-of-footprint **Alaska** gages
  (never published as real zeros) to a traceable exclusion CSV + all-gage QA companion; separate
  documented QA flags (`geom_low_confidence`/`low_pixel_support`/`partial_coverage`/`low_confidence`);
  provenance (`nlcd_collection=C1V2`, `valid_area_km2`, metadata JSON); data dictionary with
  model-noise / annual-update-seam / developed-vs-impervious-non-additive caveats. **S3 upload
  opt-in (`--upload`), held for QA.**
- **QA/QC explorer** (`viz/build_nlcd_explorer.py` → `nlcd_explorer.html`, 12.6 MB) — self-contained
  Leaflet map, 6,119 points × 32 variables; click a gage → stacked-area composition 1985–2025 +
  dashed impervious line, 2024–25 update-seam shaded; smoothed endpoint-diff deltas (labeled),
  `annual_volatility`/`shrub_grass_swap` artifact views, partial-coverage badges, excluded-45
  disclosure, normalized Shannon.
- **Codex adversarial reviews** — pipeline (FIX-FIRST → all publish-integrity blockers fixed) and
  results+plans (GO-with-fixes → all reconciled). Validation: dev%↔impervious% corr 0.963; documented
  shrub↔grass annual artifact confirmed (−0.71 swap correlation, 77–90 pp single-year swings).
- **Known gotcha (this SageMaker box):** `/home/sagemaker-user` has only ~2.6 GB free — heavy temp
  I/O must go to `/tmp`; the pipeline stages downloads there via `--scratch` + a disk preflight.
- Python only (CONUS/US); complements — does not replace — the continental MODIS LULC. **Remaining:
  human QA → S3 publish → delete `temp_lulc_conus/`.**

### New: STANDARD OUTPUT #2 — "entire period of record" WY 1980–2025 @ 60% (2026-07-22, Codex results-review GO, zero findings)
Second standard product (user decision: entire record, operationalized as
WY 1980–2025). Wrapper `docs/benchmarks/run_julia_benchmark_prod_1980_2025_60pct.jl`;
plan + full gate log: `docs/plans/production_run_1980_2025_60pct_plan.md`. Outputs
(own experiment folder per the one-folder convention):
`D:/processedOuts_1980_2025_22jul2026/streamflow_1980_2025_60pct_22jul2026_*` +
signature explorer (61.6 MB, 1,456 mapped variables) + comparison dashboard/CSV/
summary vs the April full canonical (written into the run folder via the new
`--output-dir` defaults).

- **25.1 min, 6,250 gages × 1,488 columns; annual parquet 21.82M rows.** All
  gates PASS/attributed: column delta vs canonical = exactly the 224 snow
  columns; all 6,250 gages ⊆ the canonical 7,313; independent end-cap audit
  PASS (both denominator regimes verified — gage-capped for early-ending,
  1980-anchored 46-year for late starters); annual values 525,329 pairs with 0
  mismatches/0 dup keys; QA flags clean (0 percentile/timing/seasonal-sum
  violations).
- **Cleanest cross-validation to date**: vs the April FULL canonical,
  1,189/1,227 columns Perfect at min R² = 1.0000 in every non-recession
  category INCLUDING Pettitt fields; all 38 divergent columns are the
  intentional b=1 log_a family (0 NA mismatches). Root cause of the
  cleanliness: the source parquet begins at 1980, so the canonical's uncapped
  window ≈ this run's window — run #2 is effectively the canonical analysis,
  properly capped + fraction-filtered + July features.
- **Run #1 ↔ #2 relationship** (by design of the window-start-anchored 60%
  denominator): 5,771 shared gages; 479 only-in-#2 (1980–1992 starts whose
  records ended before #1's window); 907 only-in-#1 (late starters at
  0.435–0.587 < 0.60 over the 46-year anchor — every one verified). Neither
  standard is a subset of the other. `area_normalized = FALSE`: 28 rows here /
  32 in #1 / 37 full-record (window-dependent).
- Codex results review: **GO with zero findings**; all probes confirmed
  independently (incl. the `season_excluded_years_*` −1-only pattern at 4,404
  gages vs the uncapped canonical).
- **Cross-product consistency verified (2026-07-23)**: run #1's annual parquet
  is contained in run #2's for all 5,771 shared gages (0 missing rows; the
  1.39M absent rows are exactly the 907 run-#1-only gages), with overlap-year
  values bit-identical for within-year-computable signatures. Record-dependent
  signatures (`*_all` pulses, elasticity, parameterized BFI) differ between
  windows by design — do not mix them across products. Interpretation guidance
  added to the claude-skill ("The Two Standard Output Products"), incl. why the
  longer window has fewer gages (window-start-anchored 60% denominator: 60% ×
  46 ≈ 28-year effective floor vs the 20-year floor binding in #1).

### New: STANDARD OUTPUT #1 — production run WY 1993–2025 @ 60% (2026-07-22, Codex results-review GO)
First of the two standard products (user decision; HISSS manuscript §2.2.2).
**First production exercise of the `STREAMFLOW_END_WATER_YEAR` end cap** — the
July 14 run was WY ≥ 1993 UNCAPPED (correcting the 2026-07-21 reconciliation-log
claim of a 1993–2022 cap). Wrapper
`docs/benchmarks/run_julia_benchmark_prod_1993_2025_60pct.jl`; plan + full gate
log: `docs/plans/production_run_1993_2025_60pct_plan.md`. Outputs:
`D:/processedOuts_22jul2026/streamflow_1993_2025_60pct_22jul2026_{signatures.csv,
signatures_annual.parquet, timing.json}`.

- **27.6 min, 6,678 gages × 1,488 columns; annual parquet 16.89M rows / 90
  signatures.** Gage set = July's 6,579 + 99: the cap bounds the qualifying
  denominator at 33 (partial-WY2026 parquet rows inflated uncapped denominators),
  so inclusion can only grow — all 99 verified at exactly 20/33 = 0.606.
- First run carrying the **60% overall trend gate** (gained 45 / lost 0 trend
  gages vs April on shared gages, all in the [0.60, 0.80) completeness band) and
  the **snow record-anchored decade gate** (876 snow_on_dowy / 645 melt_com_dowy
  gages trend-suppressed with means + Pettitt intact; exempt swe_max only 333).
- **Gates all PASS/attributed**: independent end-cap window audit
  (`audit_qualification.jl`, extended with an end-year arg) PASS on stratified
  edge gages; annual-values consistency PASS (545,266 pairs, 0 mismatches, 0
  floor violations, 0 dup keys); production gates 6/7 with the strict
  column-equality "FAIL" fully attributed (reference = April pre-Pettitt/pre-snow
  experiment; delta = exactly 608 Pettitt + 224 snow columns — the July products
  were lost to the flash-drive rollback, so April is the only surviving
  reference). Codex results review **GO** with 2 attributed MINORs (32-of-37
  un-normalized gages qualify under this window — the 37 count is full-record;
  2 gages lost BFI_Eckhardt mean/median to the 20-value stats floor at 17
  annual values). Dense shared means/medians EXACT vs April (max abs diff 0.0).
- **QA/QC dashboards** (same day, all in the run folder): signature explorer
  (rebuilt with the Statistic picker extended 8 → **16 stats** — the 8 Pettitt
  fields added to `build_signature_explorer.py`, so all 1,456
  signature-statistic combinations are reviewable); primary QA comparison vs the
  same-window April experiment (597/619 Perfect at min R² 1.0000, all 22
  divergent columns in the comparator's signature scope = the intentional b=1
  log_a family, 0 NA mismatches); window-sensitivity comparison vs the April
  full-record canonical.
- **Codex dashboard/results review (initial NO-GO → resolved, GO to commit)**:
  (MAJOR) the comparison DASHBOARD visualizes 4 `season_excluded_years_*`
  diagnostics the comparator excludes as metadata — their divergence vs April
  verified benign (uniformly −1 at 5,038 gages: the uncapped April run counted
  phantom partial-WY2026 seasons as excluded; the capped run correctly doesn't),
  and the scope difference is now documented in both tools' source. (MINORs
  fixed) both comparison scripts gained `--output-dir` defaulting to the
  experiment CSV's folder (one-folder convention now self-enforcing);
  `audit_qualification.jl` edge-sample slices clamped (`first(v, n)`); `temp/`
  gitignored. Explorer HTML size (64.6 MB) noted as a watch item.
- **NEW CONVENTION (user)**: all outputs of an experiment — run artifacts,
  explorer + sidecars, comparison dashboards/CSVs/summaries — live together in
  the experiment's own unique folder (this run: `D:/processedOuts_22jul2026`);
  `docs/benchmarks/` keeps only tools + long-lived reference CSVs. Recorded in
  CLAUDE.md → Critical Constraints.

### New: record-anchored decade gate for threshold-dependent snow metrics (Julia, 2026-07-22)
The 10 timing/melt/regime snow metrics (`swe_max_dowy`, `snow_on_dowy`,
`snow_off_dowy`, `melt_season_days`, `melt_rate`, `ssm`, `melt_before_peak` +
`_pct` + `_to_max_swe`, `melt_com_dowy`) now additionally require ≥
`decade_min_fraction` of the SWE-valid years in BOTH the first and last decade of
the gage's SWE record to be computable (snowy) before trend statistics are
computed. Requested by the user (2026-07-22) with four in-session decisions:
threshold **LINKED to the streamflow decade gate** (same
`na_handling.trend_completeness.decade_min_fraction` knob — no new constant, one
config line governs both); anchored to the **SWE-valid record** (rows where
`swe_max` is non-NaN — dense incl. zeros), not the metric's own span; scope =
the 10 threshold-dependent metrics only (the 4 magnitude metrics emit valid zeros
whose dense series legitimately carry snow decline); failure NaNs the **6 trend
stats only** (mean/median + Pettitt fields survive, matching streamflow trend-gate
semantics). Motivation: the own-span trend gate and the 20-value stats floor both
pass a gage whose snowy years are clustered (e.g. snow present 1981–2000, absent
2001–2010 → own span 100% complete), yielding trends conditioned on snow-present
years. Mechanism: new opt-in `force_skip_trends` kwarg on `generate_stats()`
reusing the existing skip-trends path; `calculate_snow_metrics()` computes the
gate and passes failing metrics. Denominators count SWE-valid years per window
(missing Daymet years don't count against snow presence); record span < 10 →
gate skipped; config flag `snow.record_anchored_decade_gate` (shipped `true`;
absent → disabled); collector/annual parquet unaffected (collection precedes
gating). Tests: `julia/test/test_snow_record_decade_gate.jl` (96 assertions:
kwarg mechanics, vanishing/appearing snow, valid-year-denominator pin at the
dc=0.7 separation point, linked-threshold, short-record skip, collector
invariance, magnitude exemption); full suite green (1,302). Plan:
`docs/plans/2026-07-22-snow-record-anchored-decade-gate.md`. Julia only (snow
metrics are Julia-only; gate ports with the queued snow port — note rpkg's
bundled config has no `snow` section yet, so nothing to mirror there). Output
impact at next benchmark: gated gages flip 6 trend stats → NaN per gated metric;
`flagged_for_high_na` denominators shift (accepted precedent).
**Codex adversarial review (2026-07-22): GO** — no CRITICAL/MAJOR. Confirmed the
`stats.jl` restructure is byte-identical legacy behavior across the full branch
matrix (collector/min_rows/stats-floor/trend-gate/changepoint ordering), the
anchor-set and window arithmetic (incl. span-exactly-10 and `7/10 < 0.7 ==
false` boundary), the linked-knob production path (config → benchmark →
orchestrator → snow gate), the denominator-separation test's validity, and docs
consistency. One MINOR fixed: stale rpkg-mirror step in the plan doc (marked as
executed no-op).

### Changed: trend-completeness overall gate lowered 80% → 60% (all languages, 2026-07-21)
`config/signatures_config.json` (and rpkg's bundled copy `rpkg/inst/config/`) →
`trend_completeness.min_fraction` 0.80 → 0.60; `decade_min_fraction` stays 0.80.
Source: the July 2026 guidelines doc revision AND HISSS manuscript §2.2.3, which
both state "at least 60% of the entire annual streamflow metric time series must be
complete" (first/last-decade requirements remain 80%). Config-only change: Julia
(`CFG_NA_TREND_MIN_FRACTION`), Python (`NA_TREND_MIN_FRACTION`), rpkg
(`pkg_env$na_trend_min_fraction`), and legacy R (`NA_TREND_MIN_FRACTION`) all read
the JSON at load — verified loaded values 0.60/0.80 in Julia and Python. Language-side
fallbacks stay 0.80 (legacy behavior when the config section is absent, matching the
stats-floor convention). Effect: trend statistics (slopes, correlations, p-values)
populate for gages with 60–80% complete annual series that were previously NaN'd;
mean/median/Pettitt fields and the recession/elasticity exemptions are unaffected.
Takes numerical effect at the next benchmark; the planned WY 1993–2025 and
WY 1980–2025 standard runs (see Planned) include it. Docs updated: SIGNATURES.md,
CROSS_LANGUAGE_STATUS.md, claude-skill.

**Codex adversarial review (2026-07-21)**: confirmed config pickup in all four
consumers, the safety of the added JSON `comment` key, and every Manuscript
Reconciliation Log claim (gap rejection unconditional vs config-driven negative-Q/
constant-SD flags; drainage-area normalization; gate affects trend stats only;
recession/elasticity exemption). Two fixes from its findings, both verified by
execution: (1) the legacy R trend-completeness test's "70% complete" case was
VACUOUS — it used 9 trailing NAs, but completeness is measured over the metric's
own non-NA year span, so the case could never trigger (pre-existing test bug,
its assert fails when run; rewritten as an internal-gap 66.7% case that does
trigger); (2) stale comment in `R/helperFunctions.R` claiming Python/Julia
"always pass 0.80" (they pass config-derived values). Julia mechanism tests that
pass explicit 0.8 thresholds (`test_changepoint.jl`, `test_snow_metrics.jl`) are
intentionally parameterized and unchanged. Also verified by execution: legacy R
loads 0.60/0.80 from the shared JSON.

Second (full) Codex review pass: no CRITICAL; 2 MAJOR fixed — the two non-archive
plan docs (`production_rerun_1993_2022_60pct_plan.md`, `snow_signatures_plan.md`)
still described the 80% gate as current and could mislead future reruns; both now
carry dated update notes (historical statements preserved as accurate for the runs
as executed). Notes from its MINOR findings: an **already-installed rpkg build
stays at 0.80 until reinstalled** — the package loads its bundled
`inst/config/signatures_config.json` via `system.file()`, so reinstall rpkg before
any future rpkg benchmark; historical 80% mentions in DEVELOPMENT.md /
CROSS_LANGUAGE_STATUS.md comparison tables are intentionally unchanged (they
describe April 2026 runs). Rerun of the legacy R NA-handling test suite: the
rewritten trend-completeness cases all PASS; 6 unrelated pre-existing failures
surfaced (see Known Issues).

### New: stats floor — min 20 annual values before ANY statistics (Julia)
`min_values_for_stats = 20` (config `stats_floor`): a metric with fewer non-NA annual
values emits NaN for ALL 8 statistics AND its Pettitt fields. Requested by the user
after gage 07292500's `snow_on_dowy` carried a Theil-Sen slope built on 4 clustered
years (1982–85) — the 80% trend gate measures the metric's OWN year span (4/4 = 100%),
the decade check skips spans <10, and `generate_stats` admitted ≥3 values. The floor
harmonizes with Pettitt's existing `min_total_obs=20`. **Recession and elasticity are
exempt** (inherently sparse — same exemption as trend_completeness; enforced at the
orchestration layer: the kwarg is never passed to them). Collector/annual-parquet
export is unaffected (collection precedes gating). Threaded through `generate_stats`
+ orchestrator + all non-exempt signature functions; absent config section = no floor.
Tests: `julia/test/test_stats_floor.jl` (98 assertions incl. the 07292500-shaped
regression and exemption pins). Codex review: NO-GO (2 findings: post-hoc mask left QA
flags stale → new `refresh_qa_flags.jl` recomputes all 12 via canonical
`compute_qa_flags`; legacy-path scope documented) → delta-verified GO. Post-hoc
tooling for already-delivered CSVs: `docs/benchmarks/apply_stats_floor_mask.py` (+
floor-aware `validate_annual_values.py --floor/--floor-exempt`). Ports deferred
(four-feature queue).

### New: production rerun harness + validation tooling + signature explorer
Supporting the July production reruns (WY1993+/60% and WY1980+/60%, both validated
9-gates-PASS with Codex results-review GO; plan + full execution log:
`docs/plans/production_rerun_1993_2022_60pct_plan.md`):
- **Benchmark runner**: ENV overrides for output dir (`STREAMFLOW_OUTPUT_DIR`) and
  input paths (`STREAMFLOW_DATA_PATH/_CLIMATE_PATH/_METADATA_PATH`); end-water-year
  window support (`STREAMFLOW_END_WATER_YEAR`, non-legacy path); memory patches for
  the 16 GB machine (climate frame trimmed to used columns before the water-year
  copy; raw frames released after Phase 3; per-gage preprocess-cache eviction) after
  the WY1980 run thrashed at 20.8 GB commit. Production wrappers
  `run_julia_benchmark_prod_{1993,1980}_60pct.jl`.
- **Repo-resident validation** (after on-drive copies were lost to a flash-drive
  rollback): `validate_production_run.py` (window/columns/snow/parquet/year-count/
  floor gates), `audit_qualification.jl` (independent window+fraction inclusion
  audit), floor-aware `validate_annual_values.py` (+ fixed a latent itertuples
  `_merge` crash found on its first real-data run).
- **Signature explorer** `build_signature_explorer.py`: self-contained Leaflet map of
  any run's CSV — 90 bases × 8 stats + scalars via two pickers, robust p5–p95 color
  scaling, viridis-colored distribution histogram, per-gage click panel, and annual
  time-series plots (Sen trend line + Pettitt marker) lazily loaded from
  per-signature sidecar files built from the annual parquet.
- **Finding — Daymet covers Canadian gages** (~1,100 Canadian gages carry snow
  values); earlier "US-only Daymet" doc claims corrected in SIGNATURES.md §12 and the
  claude-skill.

### New: Q-to-PPT unit gate for un-normalized gages (Julia + Python + rpkg)
Gages with no drainage area in HYDAT keep flow in raw m³/s (decision: no backfill).
`calculate_all_signatures()` gains an `area_normalized = true` kwarg in all three
languages: when `false`, ALL Q-to-PPT signatures — runoff ratios (+
`runoff_ratio_high_count`), elasticity (+ diagnostics), Q-P seasonality, and
`avg_storage` — are skipped, because Q (m³/s) and PPT (mm) units don't match, making
Q/P, dQ/dP, and cumsum(P − Q) meaningless. Q-only signatures (incl. Q##/D##
percentiles, recession, BFI) are unaffected. All three benchmark runners look up
`area_normalized` from the metadata CSV (leading-zero-safe join; only an explicit
FALSE gates, missing defaults to normalized) and pass it per gage.

- **Intentional output change at the next benchmark**: 16 of the 37 un-normalized
  gages currently carry mixed-unit climate signatures in the April 2026 golden
  (e.g. runoff ratios computed from m³/s Q against mm PPT) — those become NA.
- Elasticity and `qp_bimodality` are technically scale-invariant in Q, but are gated
  with the rest: mixed-unit inputs are conceptually invalid, and the family's
  internal PPT thresholds and P − Q terms are not scale-invariant.
- Tests (all three languages verify: 84 climate keys normally → 0 gated; gated
  output IDENTICAL to a no-climate run; gate inert when true):
  `julia/test/test_area_normalized_gate.jl` (also covers the AnnualCollector path),
  `python/tests/test_area_normalized_gate.py`,
  `rpkg/tests/testthat/test-area_normalized_gate.R`.
- **Codex adversarial review (2026-07-14)**: 2 MEDIUM findings, both fixed — (1) the
  Julia/Python runners hard-failed on metadata files predating the `area_normalized`
  column (now header-guarded, matching rpkg's existing guard); (2) explicit-FALSE
  parsing diverged across languages for 0/1-serialized booleans (now harmonized in
  all three runners: Bool false / numeric 0 / "FALSE"/"0" strings gate; anything
  else, incl. missing, defaults to normalized). Fix verification surfaced a third
  latent bug the review missed: the Julia runner's join-id helper was typed
  `::String` but CSV.jl yields InlineStrings — Phase 2b would have crashed on the
  production metadata (loosened to `AbstractString`). Runner parsing was then
  functionally verified in all 3 languages against TRUE/FALSE, 0/1, and
  missing-column metadata variants, plus the production CSV (16,994 gages parsed,
  exactly the expected 1,601 gated).

### New: snow metrics signature family (14 metrics, Daymet SWE — Julia)
Fourteen per-water-year snow metrics from daily SWE, through the standard
`generate_stats()` path (8 stats + 8 Pettitt fields each → **+224 columns**; 76 → 90
time-series signatures; annual values exported automatically via the collector). New
module `julia/src/snow.jl` (`calculate_snow_metrics`); full spec + review record:
`docs/plans/snow_signatures_plan.md`. Requested by the user with four in-session domain
decisions (2026-07-14).

- **Metrics**: `swe_max`, `swe_max_dowy`, `snow_cover_days`, `snow_on_dowy`,
  `snow_off_dowy`, `melt_season_days`, `melt_rate`, `ssm` (Hatchett 2021, Hydrology
  8(1):32), `swe_apr1`, `melt_before_peak`, `melt_before_peak_pct`,
  `melt_before_peak_to_max_swe`, `melt_com_dowy`, `swe_max_to_ppt`.
- **Domain decisions**: SWE ≥ 10 mm day threshold applied to durations AND magnitudes
  (thresholded series SWE*; sub-threshold years are operationally snow-free);
  snow-on/off anchored to the spell containing the annual peak (boundary-censored →
  NaN); melt rate = max SWE ÷ melt-season days; SSM seasonal spells ≥ 60 continuous
  days. Config: `signatures_config.json` → `snow` + `na_handling.snow_na_policy`.
- **Plumbing**: Daymet `swe` was previously dropped at the climate join — now
  normalized (`swe`→`SWE`) and carried through; `preprocess_daily_data()` gained a
  per-year SWE policy (mirrors PPT: >30 NAs, ≤3-day interpolation, negative rejection)
  and a new `valid_swe_years` return field. Snow metrics run ONLY on an explicitly
  passed `snow_data` frame filtered to those years — NO implicit fallback to the gage
  frame (plan review finding: prevents SWE-invalid years leaking in). Canadian gages
  (no Daymet) → all snow columns NA.
- **Known existing-column change (accepted)**: `flagged_for_high_na` counts the 224 new
  NA fields in its denominator, so some no-SWE gages will flip at the re-run
  (April 2026 Pettitt precedent; whitelisted + quantified; semantics pinned by test).
- **Two-round Codex plan review**: initial NO-GO (2 CRITICAL: orchestrator fallback
  leak; the false "no golden divergence" claim) → 7 findings + 2 delta-round residuals
  incorporated → GO. Record: plan doc §12.
- **Tests**: `julia/test/test_snow_metrics.jl` (~330 assertions: exact hand-derived
  gages — triangle, mid-winter thaw, spell arithmetic, threshold splitting, boundary
  censoring, tie rule, leap-year April 1, PPT gating incl. legacy runoff-ratio parity,
  preprocessor `valid_swe_years`, no-fallback + 0-row schema contract, SWE-invalid-year
  leak regression, collector/zero-warn invariants, QA-flag semantics pin);
  `smoke_test.jl` now joins SWE and validates snow presence/ranges.
- Julia only; Python/rpkg ports deferred (three-feature port queue). Benchmark re-run
  pending — bundled with annual values + b=1 alpha.

### Docs: data-source inventory
New `docs/DATA_SOURCES.md` — table of all 11 external data sources (data class /
project use / provider / original access): USGS NWIS, HYDAT, Caravan, Daymet,
GAGES-II boundaries + attribute tables, ECCC MDA_ADP basin polygons, HydroBASINS,
HydroATLAS/BasinATLAS v10, MODIS MCD15A3H (LAI) and MCD12Q1 (LULC). Also records
the guidelines Google Doc as a governance input and the deliberately excluded
boundary cases (planned ERA5/PRISM, S3 mirrors, viz-only CDNs, unused geometry
backstops).

### Changed: recession alpha now assumes a linear reservoir (b = 1) — Julia
All alpha outputs (`log_a_pointcloud`, `log_a_events`, the 6 `log_a_seasonality_*`
scalars) are now computed with the recession exponent **fixed at b = 1** across all
locations and periods: `log(a) = median(log(-dQ/dt) - log(Q))`, no regression. The
exponent analyses are unchanged — `b_pointcloud`, `b_events`, and `concavity` keep
their free power-law fits, as do event identification, the 25-event gate, and
`alpha_linear` / `recession_alpha_point_cloud_linear_reservoir` (already b=1).

- **Rationale** (domain decision): log(a) is the intercept of a regression whose slope
  is b, so free-fit alpha estimates are convolved with b — trends/seasonality in alpha
  partly reflected changes in b. Fixing b = 1 decouples them.
- Column names unchanged (frozen CSV contract) — methodology change behind existing
  columns, like the April 2026 recession rewrite. At the next benchmark, log_a columns
  will diverge from golden intentionally; b/concavity/alpha_linear must not.
- Consequence: under the forward-difference discretization each b=1 point equals
  `log(1 − Q_{i+1}/Q_i)`, so per-year `log_a_pointcloud` APPROXIMATES
  `log(1 − alpha_linear)`. The relation is approximate, not a same-sample identity:
  `alpha_linear` also includes events whose free power-law fit failed (the point cloud
  only pools fit-success events), and medians only commute exactly with the monotone
  transform at odd pair counts.
- New tests (`julia/test/test_recession_alpha_b1.jl`): exact helper values on a pure
  exponential; a quadratic-reservoir gage (−dQ = k·Q² exactly) proving b stays 2 under
  the free fit while alpha matches an independent b=1 recompute (and is far from the
  old free-fit intercept); linear-reservoir continuity (b=1 gage values unchanged).
- Julia only; Python/rpkg sync deferred (bundled with the annual-values port).
  Design note: `docs/plans/recession_alpha_linear_reservoir_plan.md`.

### New: annual values export (per-year signature values, Julia)
The per-year annualized metric values — previously computed inside every signature
function and discarded after `generate_stats()` collapsed them into the 8 summary
statistics — are now saved as one long-format parquet alongside the summary CSV:
`{prefix}_signatures_annual.parquet` with columns `gage_id, signature, water_year, value`
(~20M rows expected for 7,313 gages × ~76 signatures, ~170 MB zstd).

- **Mechanism**: opt-in `AnnualCollector` threaded through `calculate_all_signatures()`
  and all 13 signature functions into `generate_stats()` (mirrors the `changepoint`
  kwarg plumbing). Collection happens BEFORE min_rows / trend-completeness gating and is
  read-only w.r.t. the statistics — the summary CSV contract is untouched (verified by a
  with/without-collector equality test).
- **Semantics**: the export records the exact series `generate_stats()` consumed. Dense
  signatures carry NaN placeholder rows ("year qualified, metric not computable");
  caller-pruned signatures (`qp_bimodality`, `n_recession_events`, flashiness, storage,
  baseflow, elasticity) omit those rows — absent row ≡ NaN for consumers.
  `elasticity_rolling` is keyed to the END year of each 11-yr window; `elasticity_annual`
  to the later year of each consecutive pair.
- **Config**: `config/signatures_config.json` → `annual_values.save` (shipped `true`;
  defaults to `false` when the section is absent). Constant `CFG_SAVE_ANNUAL_VALUES`.
- **Validation**: new tests in `julia/test/test_annual_collector.jl` (stat-identity,
  zero-warnings, signature coverage, NaN/missing-year guards, no duplicate keys) +
  `docs/benchmarks/validate_annual_values.py` (recomputes mean/median from the parquet
  and cross-checks the summary CSV; runs after the next benchmark).
- **Pre-flight**: Parquet2 pinned (`0.2`, resolved v0.2.27) after a write→read round-trip
  check (zstd/snappy/uncompressed, NaN-safe, empty-frame-safe; 1M rows = 8.5 MB, 0.2s).
- **Two Codex reviews**: plan-stage (NO-GO → 5 amendments incorporated pre-implementation)
  and code-stage (core implementation confirmed correct; 3 harness findings fixed:
  validator ≥3-non-NaN-rows coverage rule, zero-gage no-op exit, and a deterministic
  exponential-recession test gage covering the recession + parameterized-BFI collector
  paths the noisy synthetic gage never exercised). Full records:
  `docs/plans/annual_values_export_plan.md` §9–§10.
- Julia only — Python/rpkg ports deferred.
- **Benchmark re-run pending** — bundled with an upcoming workflow feature.

## [June 2026]

### New: per-watershed MODIS Earth Observation products (LAI + LULC, non-signature)
Two new per-gage EO products aggregated over each gage's upstream watershed from MODIS Collection 061, sitting alongside the signatures and joining by leading-zero-safe `gage_id`/`canon_id`. Python (not a cross-language signature — like HydroATLAS, this is metadata/ingestion). 7,964 watersheds. New subproject `EO_data_processing/` (its own README, pipelines, manifests). All artifacts on S3 `s3://climate-ai-data-science-shiny-app-data/streamflow/`. Full detail: [EO_data_processing/README.md](EO_data_processing/README.md).

- **Watershed geometry layer** (prerequisite) — official basin polygons primary (US GAGES-II + Canada ECCC MDA_ADP), HydroBasins delineation fallback; 7,964 basins (8,018 universe − 54 >100,000 km²). `watershed_polygons_26jun2026.{gpkg,parquet}`. Codex-reviewed.
- **LAI (MCD15A3H, monthly, 2002–2024)** — 270-month panel, 2,150,280 rows (7,964 × 270). Two-stage spatial-distribution-of-monthly-mean (day-weighted per-pixel monthly mean → coverage-weighted basin mean + spatial heterogeneity stats), QA fractions from FparLai_QC/FparExtra_QC. Catalog (GEOLATELY/us-east-2) + LP DAAC backfill for far-N + 2024-11; 17 urban basins permanently NA (MODIS fill code 250). Codex-GO.
- **LULC (MCD12Q1, annual, 2001–2024)** — 191,136 rows (7,964 × 24), 109 cols: per-class % coverage for all 8 classification bands (LC_Type1–5 + LC_Prop1–3, 102 class columns) via coverage-weighted exactextract over VRT-mosaicked sinusoidal tiles. All 8 bands sum to 100 within 1e-13; manifest matches v061 spec. **Codex-reviewed 29 Jun (gpt-5.5) = GO**; 2 latent code hazards (tile-success accounting, unknown-code drop) hardened with the delivered output byte-identical.
- **Cross-cutting**: leading-zero canonical key reconciliation (metadata zero-stripped vs signatures zero-padded); LP DAAC HDF4 read via `pyhdf` (pip-GDAL lacks the driver); Earthdata account-lock handled by completeness-guard + retry/backoff (self-healing, resumable).
- **LAI `good_coverage_frac`** (30 Jun) — per basin-month continuous QA flag = valid pixel-obs / (all basin pixels × all expected composite dates); generalizes the binary `partial_month`. Added to `watershed_modis_lai_monthly_30jun2026.parquet` (data otherwise unchanged, hash-verified).
- **LAI QA/QC explorer** (30 Jun) — self-contained Leaflet HTML (`EO_data_processing/viz/build_lai_explorer.py` → `watershed_modis_lai_explorer.html`) for manual review: 7,964 points colored by any of 17 QA/summary variables, click a gage for its full monthly LAI time series (quantile bands + per-month good-coverage strip), map/TS hovers, glossary, processing notes. Mirrors `streamflow_explorer.html`.
- **LULC QA/QC explorer** (30 Jun) — self-contained Leaflet HTML (`EO_data_processing/viz/build_lulc_explorer.py` → `watershed_modis_lulc_explorer.html`, ~49 MB, on S3). 7,964 points colored by any of **27 map variables** in two non-overlapping groups: 10 IGBP-derived summary features (`class_diversity`, `delta_forest/urban/cropland` [robust mean(2001–03) vs mean(2022–24) endpoints], `dominant_class`, `dominant_group` [lump-robust aggregate dominance — fixes IGBP's class-splitting bias], `dominant_change`, `n_modis_pixels`, geom source / low-confidence) + all 17 individual IGBP classes (static `pct_*` roll-ups dropped as redundant with the per-class maps). Click a gage → stacked-area land-cover composition 2001–2024 with a **band switcher across all 8 MODIS schemes** (IGBP/UMD/LAI-fPAR/BGC/PFT/FAO-LCCS×3), official MODIS IGBP palette, per-year class-breakdown hover (0.1 % resolution). Glossary carries the MODIS interpretation caveats (SE-US woody-savanna behavior; 2021+ reduced confidence). Headless-validated (puppeteer/Chrome). **EO pipeline + viz COMPLETE.**

## [May 2026]

### New: per-gage watershed HydroATLAS metadata (non-signature product)
A standalone metadata file pairing each gage with the hydro-geophysical character of its **entire upstream watershed**, from HydroATLAS / BasinATLAS v10 (281 source attributes → ~211 output columns: climate, hydrology, terrain, land cover, soils/geology, anthropogenic). Sits alongside the signatures, keyed by leading-zero-safe `gage_id`: `data_out/watershed_hydroatlas_metadata_{date}.{csv,parquet}` + a generated data dictionary. 8,014 gages, ~13 MB.

- **Hybrid aggregation** over each gage's delineated upstream basin set (`upstream_hydrobasins.rds`): HydroATLAS upstream-accumulated `_u` / pour-point `_p` values pass through from the outlet basin (91 attrs — the delineated set equals HydroATLAS's own upstream extent, so `_u` is the rigorous watershed value); `_s`-only attributes are SUB_AREA-weighted (92 area-weighted means incl. monthly climate, with the HydroATLAS `-9999` NoData sentinel masked; elevation as spatial min/max; 11 categoricals as area-weighted majority — glc/pnv via argmax of upstream `_u` fractions, wet + the other 8 via area-weighted mode of the source class). Percentage / areal-extent fields are always area-weighted mean (never max). A `watershed_area_rel_diff` column ( |summed member SUB_AREA − gage `basin_area`| / `basin_area`, gage-reported drainage as truth) flags gages where the delineated watershed area diverges from the reported drainage (mis-snapped outlets).
- New `R/aggregate_hydroatlas_metadata.R` + `run_hydroatlas_metadata.R`; `HYDROATLAS_GPKG` added to `config.R`. R-only (metadata/ingestion) — not a cross-language signature, so no Julia/Python port.
- Validated: member SUB_AREA sums to outlet UP_AREA (median rel.diff 1e-4); area-weighted `_s` reproduces HydroATLAS `_u` (elevation cor 0.9994); 100% join to the signatures file (5,707/5,707). See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Watershed HydroATLAS Metadata.
- **Independent code review** (in-session) confirmed the classifier, NA-aware weighting, BFS delineation, and leading-zero joins; correctness fixes applied — `wet_cl_smj` switched from `_u`-fraction argmax to area-weighted mode (its fraction basis omits classes 10/11), `-9999` NoData masked in weighted means (`sgr_dk_sav`), and a malformed legacy-cache member id sanitized.

### New: static HTML data explorer
A self-contained, double-clickable `data_out/streamflow_explorer.html` (Leaflet + Canvas) for exploring the gages: 8,014 points colored by any of **412 variables** (dataset switcher between HydroATLAS watershed metadata and streamflow signatures); click a gage to draw its **watershed boundary** (border-only). Generators: `build_explorer.R` (unions the RAW lev12 basins per outlet for a clean dissolve — no interior sliver lines — then strips holes, simplifies, and caps at ~250 vertices → `watershed_borders.geojson` ~9.5 MB), `assemble_explorer.R` (joins points + injects into `explorer_template.html`). Full-resolution borders are ~111 MB (too heavy to inline); the clean-dissolve build is ~28 MB total. Border click-through and clean rendering verified headlessly (puppeteer-core + Chrome). See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Static HTML Explorer.

---

## [April 2026]

A large release: new signature families, the canonical-language transition, centralized NA handling, and cross-language alignment to 99.5% agreement. Output grew from 551 to **1,264 columns** (656 base/metadata + 608 Pettitt changepoint) across 7,313 gages. **Full per-change detail is in [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md) → [April 2026].**

### New signatures
- **Pettitt changepoint detection** applied to all 76 time-series signatures (8 fields each → 608 columns: cp year, p-value, pre/post means, delta, pct change, pre/post MK p-values). ~13.4% significance rate (2.7× null); signal-robustness analysis in [docs/SIGNATURES.md](docs/SIGNATURES.md) → Changepoint Detection.
- **Recession-parameterized baseflow** — `alpha_linear`, `BFI_Eckhardt_param`, `BFI_LyneHollick_param`, plus the per-gage `recession_alpha_point_cloud_linear_reservoir` scalar (Collischonn & Fan 2013). +25 columns.
- **Guidelines Section 3 additions** — D1/D99 flow timing, `n_recession_events`, `elasticity_annual`, `negative_ann`, `runoff_ratio_high_count`, 4 seasonal exclusion-year diagnostics, and 2 elasticity diagnostics.

### Canonical language transition
- **Julia replaced R as the canonical implementation** (~40× faster, modular). `rpkg/` is the production R port; monolithic `R/helperFunctions.R` is deprecated (ingestion utilities only). See [CLAUDE.md](CLAUDE.md).

### Centralized NA handling
- New `preprocess_daily_data()` runs once per gage before signatures: daily-grid normalization, ≤3-day interpolation, year rejection, seasonal completeness flags. Removed all per-signature min-days thresholds and `fillna(0)` calls; added the 80% trend-completeness gate to `generate_stats()`. Negative-Q rejection is config-driven (off by default); constant-SD is a flag only. See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → NA Handling Architecture.

### Cross-language alignment
- Synced Section 3 plus the recession-algorithm rewrite to Python and rpkg; fixed Julia's Mann-Kendall tau-b, a `leftjoin` row-order bug, and several R-canonical divergences. The three synced implementations (Julia, Python, rpkg) agree at 99.5% R² (615/623 columns Perfect). Benchmark tables: [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Cross-Language Benchmarks.

### Tooling & audits
- Julia sensitivity-experiment framework (WY≥1993 ± qualifying-fraction filters), a Julia-vs-golden comparison pipeline, and interactive HTML validation dashboards (`docs/benchmarks/`).
- Three audits: Guidelines-vs-implementation, NA-handling (Codex), and a Julia adversarial review — all findings since addressed or documented.

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

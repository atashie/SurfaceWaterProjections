# Codex plan review, round 4 — plan version 5 (delta review after round 3; 2026-10-26; gpt-6-astra, reasoning high, read-only)

Brief: delta review of v5 — §10.2 (the span-AND-water-year-set restriction; whether it keeps every Q-only output unchanged), §10.4 (one tree with the `STREAMFLOW_CLIMATE_FRAME` switch; S1 on the splice-affected gages + control sample; the composed S1′ batch), §10.5 (ledger of changed preprocessed daily values; dependency closure by output kind; gates 5a–5f; magnitude census), §10.8 (coupled Q–P metrics), the consistency of §0 / §4 / §5 / §8 / §9 with §10, and decisions 12–17. Verdict and findings as returned:

**GO WITH CHANGES.** Keep the defaults of decisions 12–17, but correct the stage construction, dependency registry and acceptance gates below before implementation.

Disposition of the eight round-3 findings:

1. **RESOLVED — water-year support.** The span-plus-original-year-set restriction preserves the production Q grid and year loop. My Julia fixture preserved valid/rejected years, diagnostics, seasonal flags, 624 Q-only result fields and window denominators.
2. **RESOLVED — production versus legacy joins.** §10.2 now targets the production join, leaves legacy behavior explicit and requires helper failures to surface.
3. **PARTLY RESOLVED — daily ledger and numerical gates.** Secondary interpolation changes are covered, but preprocessed-value differences cannot replace raw provenance, and the degree-day delta rule remains incorrect.
4. **PARTLY RESOLVED — coupled metrics.** Direction gates and the small-residual claim are removed; §10.8 still overstates locality, storm exclusion and what exposure measures establish.
5. **PARTLY RESOLVED — attribution stages.** One tree and a reported overlap residual address the main design errors. Subset identity, missing stage coverage and contradictory Phase-D instructions remain.
6. **PARTLY RESOLVED — dependency closure.** Drought, pulses and annual elasticity are addressed, but raw seasonal masks, collector gates, several scalars and temperature anchors remain insufficiently specified.
7. **RESOLVED — forecast corrections.** I reproduced both forecast CSV sections row for row from the supplied census and product cohorts, including totals and decade counts; restored-year gaps also reproduce as 26/59 and 15/26 years/dates.
8. **PARTLY RESOLVED — census and release integration.** The statistical definitions and deposit inventory address most concerns; I reconfirmed 695 p-value columns per product. Validation coverage, stale instructions and scheduling arithmetic remain unresolved.

On the specific Q-only question: **yes, the restriction is sufficient for the current production preprocessing contract**, assuming unique date keys and consistent date-derived water years. Inserting missing-Q rows inside an existing partial year produces exactly the same Q grid as absent rows did. All-NA years retain their existing rejection and seasonal flags. Years outside the window remain outside it; windows without rows retain the same empty denominator branch. `years`, `total_possible` and `season_excluded_years_*` therefore remain unchanged. `flagged_for_high_na` is a mixed-output flag, not a Q-only invariant.

New findings, most severe first:

1. **HIGH — §10.5 items 1–3: the preprocessed ledger loses causally necessary information.**

   A restored observation can equal its previous interpolate exactly. In my Julia fixture, restoring one such winter observation preserved every preprocessed Q value and the valid-year set, but changed winter completeness **false→true** under the configured 80% threshold. Seasonal outputs and exclusions can consequently change without a daily-value delta or year-set change. The implementation computes completeness from [raw Q](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/io.jl:326).

   **Change:** retain the raw splice ledger and inserted-date inventory alongside the preprocessed ledger. Persist per-stage raw missingness/interpolation masks, seasonal flags and diagnostics, including rejected years. Derive affected membership and gate 5a from all these inputs. Preprocessed missingness also cannot identify every interpolated-Q exposure requested by item 7.

2. **HIGH — §10.4/5e: subset execution does not yet establish the stated identities.**

   The runner unions result keys across gages before computing QA; `flagged_for_high_na` uses that union as its denominator. My fixture changed this flag **false→true** for identical existing gage values merely by adding structurally missing signature columns. The [batch harness explicitly documents this dependency](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:11). Matching headers within a subset run does not prove equality with the full-run header.

   Moreover, S1 and S1′ do not contain “every gage without a splice change.” Those universal comparisons are unavailable as written.

   **Change:** pin or verify the full-product signature schema before QA for every stage, preserving v2’s denominator. Define full-cohort S1 as executed affected/control rows plus explicitly sourced S0 rows elsewhere. Restrict measured S1′ identities to its executed cohort; distinguish inferred equality elsewhere. Record parent-cohort and subset manifests separately, and record stage-specific input hashes rather than implying those hashes match.

   Compose the interaction batch from realized support too: I found restored-year gap candidates outside the baseline join set on **18/13 additional gages**, respectively. The fixed 147/122 overlap counts do not cover these candidates.

3. **HIGH — §10.5 item 3: the registry is not complete at the output level.**

   Recession events **do not span water years** here: event detection runs inside each year. However, the record-wide minimum-event gate precedes annual collection. My fixture added a 25th fitted event/year and caused `log_a_events` to appear for **24 unchanged earlier years**. See the [gate before collection](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/recession.jl:515).

   Q–P seasonality likewise has whole-record year/usable-result floors before collection. Snow’s record-anchored gate affects trend statistics, not annual values, means or changepoint fields; temperature has an analogous frost gate anchored to temperature-valid years. Ordinary trend completeness uses each metric’s finite-value span, so adding years can enable or disable it.

   **Change:** enumerate dependencies separately for numerical values, annual-key emission, summaries and scalar fields. Include recession seasonality/alpha, drought thresholds, static elasticity and its counters, and `runoff_ratio_high_count`. Preserve drought’s smoothed-record and pulse `*_all` dependencies already identified. Classify old-year keys newly emitted by family gates separately from newly Q/climate-eligible years. Do not permit arbitrary shared-year snow annual changes merely because an anchor changed.

4. **HIGH — §10.5 gates 5b–5d: numerical validation still has an incorrect formula and uncovered outputs.**

   `gdd_cum_*` sums **thresholded daily degree-days**, not raw temperature deltas. My fixture changed TMIN and TMAX from −2°C to +2°C on one day: annual GDD changed **+2**, versus **+4** for mean temperature and **+8** for summed TMIN/TMAX changes. This follows the [implemented formula](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/climate.jl:476).

   **Change:** sum
   `max(0,Tmean_new−base) − max(0,Tmean_old−base)`
   over the relevant dates, reconstructing both temperature components.

   Gate 5b covers totals only on shared valid years; 5c’s “every other” wording leaves newly eligible totals without an explicit numerical check. Gate 5d counts gains but does not validate their values. Require absolute recomputation for **all gained values**, including totals.

   Finally, E7’s annual-value comparison plus `check_schema_equality.py` does not validate changed scalar values or all summary statistics: the latter reads only headers and gage IDs. Add numerical checks for these outputs and explicit gate-state comparisons, with defined tolerances and existing port residual allowances.

5. **HIGH — §§0, 4 and 5: conflicting executable instructions remain.**

   §0 still limits affected changes to oracle-named years and promises general shared-year climate equality. E5’s U excludes both fixes, while A still includes only Q-spliced gages; join-only gages consequently fall into neither class. Phase D still names the pre-A′ tree and “batch 1,” contradicting §10.4. E4 still requires `qualified = selected = emitted`, contradicting the frozen-cohort policy.

   **Change:** rewrite these passages around the final contract: A is the complement of U; climate equality is S0→S1 on shared eligible years; all stages use the switched tree and composed cohort; `emitted = selected ⊆ qualified`. Resolve whether baseline drift stops execution or changes the reference—§10.4 requires S0=v2, while D0 permits replacing v2 after explanation.

6. **MEDIUM — §10.8: several scientific claims remain overstated.**

   “On the affected gage-years only” conflicts with record-normalized elasticity and family gates. Old P was interpolated on absent-row days, not every interpolated-Q day. Missing rows do not establish exclusion of storms whose runoff Q retained. Interpolated-day counts and associated rainfall quantify **exposure**, not the unknown runoff-error residual.

   **Change:** state these narrower claims explicitly. Limit any ratio identity to matched annual/seasonal eligibility and denominators; a ratio of summary means need not equal a mean runoff ratio. Retain the production complete-Q-series statement and decision 17.

7. **MEDIUM — §8: November 9 is not supported by the stated sequential budget.**

   The table sums to **11.5 person-days**, not 11. There are only **10 weekdays before November 9**, or 11 including submission day, starting October 26. “Starts this week” therefore does not support completion by that date for one person.

   S1’s stated 60–90 minutes across both windows plus S1′ exceeds the table’s 1–1.5 hours. The download allowance also needs inclusion or explicit overlap in the unattended total. Registry completion, state capture and expanded numerical checks make the hands-on estimates optimistic.

   **Change:** provide a dated critical path with explicit parallel staffing/weekend assumptions or revise the completion claim. Frozen gage counts do not guarantee unchanged manuscript signature statistics.

The gates are feasible below 3 GB using gage-partitioned replay and column-wise accumulation. **They are not fully computable from the currently specified artifacts:** raw masks/state, complete registry coverage and stage/schema provenance are missing. Gates 5a and 5d are legitimate coverage/accounting gates, but cannot substitute for numerical validation; 5b can reject correct GDD results, 5e overstates available comparisons, and 5f must distinguish full-product cohort equality from subset membership.

Verified independently:

- Production preprocessing, runner, QA and relevant family/collector/gate code.
- Julia fixtures for year-support invariance, raw-mask changes, subset QA, GDD and recession emission.
- Forecast row equality, totals, decades, restored-year gap candidates, p-value counts and schedule arithmetic.
- Read-only work; no internet, file changes or production runs.

Taken on trust:

- Full-domain census qualification labels and prior-round climate-value/eligibility replay results.
- Source retrieval/restoration assumptions and runtime estimates.
- Future three-language helper parity and implementation of the proposed artifacts and gates.

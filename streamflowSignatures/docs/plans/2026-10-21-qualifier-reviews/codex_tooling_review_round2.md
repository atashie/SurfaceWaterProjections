# Codex tooling review, round 2 (2026-10-26, gpt-6-astra, high reasoning, read-only) — after the round-1 and Sonnet round-2 fixes, with the full Phase B outputs

Brief: `scratchpad/codex/review_brief_B_C_tools_v2.md` of the session plus its round-2 addendum. Verdict: Phase B outputs GO WITH CHANGES (data independently reconfirmed), Phase D on this tooling NO-GO; dispositions in the plan's §11 progress log. Verbatim report:

196,625
**NO-GO for Phase D. The current Phase B data passed my independent checks, but several downstream gates still accept corrupted results.** No files were modified; corruption tests used in-memory copies.

The completed full audit now reports PASS. Independently, I verified:

- All **193,542 restorations** exactly equal accepted cached values × metadata conversion.
- All **111,624,189 rows**, **909 row groups**, schema metadata, row order and non-Q columns are preserved. Changed Q cells exactly match the ledger.
- The compile rule reproduces original and patched metadata for all **987 gages**.
- Patched replay matches **322,563 census year statuses** and **16,028 window decisions**, including the oracle’s exact final year sets.
- Nonzero-only B2 comparisons pass for every restored gage; minimum matching share is **66.984%**.
- The subset annual E7 comparison passes with the documented single `D40_day` residual.

1. **BLOCKER — An annual residual exempts arbitrary summary corruption.**  
   [check_summary_equality.py:125](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/check_summary_equality.py:125)

   Every summary statistic belonging to a `(gage, base)` with an annual residual is removed from numerical gating. Using the real subset artefacts, I changed `02160390 / D40_day_mean` from **142.34375 to 265.34375** in memory. The existing one-day annual residual remained unchanged; the summary gate still returned **PASS, exit 0**.

   **Fix:** Recompute the affected summaries from each annual series and validate the actual propagated difference. A named annual residual must not authorize unlimited differences in all its summary statistics.

2. **BLOCKER — The delta checker does not validate the eligible values it exports.**  
   [preprocess_state_ledger.py:249](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/preprocess_state_ledger.py:249)

   Raw-change branches validate `raw` against the source ledger but do not verify `clean` or `pre`. Secondary-interpolation and eligibility branches check `clean`, while exporting unchecked `pre` as `new_value`.

   I added **123** to `PPT_pre` at `08166140 / 1984-12-31`, preserving its correct raw and clean values. The checker reported **“every change attributed and verified,” exit 0**.

   **Fix:** Independently validate every layer: raw against source evidence; clean against preprocessing; pre against clean plus the variable’s eligibility. Perform these checks regardless of the attributed cause.

3. **MAJOR — State-field “reproduction” remains a permission check.**  
   [preprocess_state_ledger.py:344](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/preprocess_state_ledger.py:344)

   Missingness masks, longest runs, interpolation masks, seasonal flags and climate validity are generally accepted when a related daily change exists. Their resulting values are not reconstructed. I changed `ppt_raw_max_run` to **999** for `08166140 / WY 1985`; the checker still returned **COMPLETE, exit 0**.

   **Fix:** Rebuild masks, runs, seasonal fractions and validity from each stage’s daily arrays, then compare exact expected fields. Also reject duplicate state-table keys.

4. **MAJOR — Legitimate splice changes fail because rejection vocabularies differ.**  
   [preprocess_state_ledger.py:352](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/preprocess_state_ledger.py:352)

   `q_rule()` produces Julia-style `too_many_na` / `gap_too_long`, but replay uses Python’s `raw_na_count=…` / `max_na_run=…` reasons.

   Replaying the actual four restorations for `01021470` gives WY 2019 rejection changing from `raw_na_count=153 > 30` to `raw_na_count=149 > 30`. All four daily changes validate, yet the state comparison returns **UNATTRIBUTED: rejection, exit 1**.

   **Fix:** Compare structured rejection codes and numerical diagnostics, using one explicit mapping between implementations.

5. **MAJOR — E7 hashes bind files, but not validation coverage.**  
   [assess_rerun_changes.py:546](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:546), [check_summary_equality.py:93](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/check_summary_equality.py:93)

   A summary verdict can hash the complete reference CSV while comparing an arbitrary `--gages` subset. Selecting a nonexistent gage against the supplied corrupted candidate produced **PASS with zero gages and zero compared cells**. The census would accept that verdict’s matching reference hash.

   Conversely, annual E7 requires identical complete key sets, while census binding requires the full assessed file’s hash. Thus a separately materialized affected-gage reference cannot certify the planned full-product census.

   **Fix:** Record and verify exact compared gages/keys, window, exclusions, tolerances and both input identities. Require complete coverage of the census’s affected cells. Support a verifiable projection of the full product for subset E7.

6. **MAJOR — The summary gate accepts finite/infinite disagreement.**  
   [check_summary_equality.py:119](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/check_summary_equality.py:119)

   For finite versus infinity, both the difference and relative tolerance become infinity; `inf > inf` is false. A one-cell comparison of **1.0 versus +Inf** returned **PASS**, without any waiver.

   **Fix:** Handle infinities explicitly, as the annual checker already does: require identical infinity signs, or reject nonfinite outputs where prohibited. Apply numeric tolerance only to finite pairs.

7. **MAJOR — B2’s new per-gage scale proof counts zeros.**  
   [splice_nwis_qualifiers.py:639](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:639), [splice_nwis_qualifiers.py:809](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:809)

   Despite the comment, `n_compared` includes zero flows, which cannot establish the conversion factor. Doubling `08330600`’s area from **335.6 to 671.2** still yields **10,513 / 11,947 = 87.997%** equality and passes the 50% gate. Sixteen restored gages have at least half their stored values equal to zero.

   The audit also trusts the control CSV without binding it to the current metadata/cache.

   **Fix:** Require nonzero comparison evidence per gage and recompute it—or hash-bind the control inputs. My independent nonzero-only check confirms the **current** outputs pass this stronger test.

8. **MAJOR — B8 can silently treat an omitted gage as source-unavailable.**  
   [splice_nwis_qualifiers.py:922](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:922)

   Missing entries in `records` become `ok=False`, rather than an evidence failure. B1 only checks the supplied universe; B3 checks category totals without failing on `no_record`.

   Executing B8 with the five-gage test’s targets, but omitting `02338660` from the ledger and opened-source universe, returned **PASS**: **1,391 targets, only 394 considered eligible**, with the omitted 997 restorations unnoticed.

   **Fix:** Derive required source-gage coverage independently from the base targets and masked audit rows. Missing source evidence must fail completeness, not establish ineligibility.

9. **MAJOR — The fingerprint still misses supported metadata inputs.**  
   [run_batched_julia.py:254](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:254)

   GAGES-II content is hashed only when supplied through `--gages-ii-dir`; Julia also resolves it from `STREAMFLOW_GAGES_II_DIR` or configuration. With the environment route, my fingerprint probe never called the directory hasher and remained unchanged across simulated content changes.

   HYDAT fingerprinting similarly uses an environment/default path, whereas Julia reads `metadata.hydat_path` from configuration.

   **Fix:** Resolve paths exactly as the runner does, then hash the resolved dependencies. Check the execution snapshot again after each batch, including the last.

10. **MAJOR — Stale or absent annual output can still satisfy the runner.**  
    [run_batched_julia.py:389](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:389), [run_batched_julia.py:477](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:477)

    Freshness checks cover CSV and timing only. Fresh CSV/timing plus an older annual parquet returns `outputs_ok=True`; `--require-annual` checks existence only. Resume skips before enforcing `--require-annual`, so a previously successful run without annual output can also be skipped.

    Merge does not recheck the recorded output hashes.

    **Fix:** Require fresh annual output when requested, enforce that requirement before resume, and validate hashes during merge. Prefer isolated attempt outputs published only after validation.

11. **MAJOR — Emission explanations use different counters from the library gates.**  
    [assess_rerun_changes.py:237](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:237)

    Recession attribution sums `n_recession_events`, which counts detected events; the actual minimum-event gate uses **successfully fitted events**. Q–P attribution counts finite emitted `qp_bimodality` years, whereas the library separately checks input-year count, usable annual rows and finite bimodality rows. The registry’s collector-floor cases have no corresponding explanation branch.

    **Fix:** Export the actual per-stage family gate counters and decisions. Validate key gains/losses against those states, including collector floors. Output presence is not independent evidence of the gate that produced it.

12. **MAJOR — GDD gate 5b explicitly passes a contradictory result.**  
    [assess_rerun_changes.py:523](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:523)

    With no temperature delta, this branch unconditionally reports PASS and merely prints the number of unequal GDD cells. Executing it with an annual GDD change **100 → 223** produced:

    `PASS … every GDD cell must be unchanged: 1 non-equal GDD cells`

    The temperature-change branch also permits zero comparisons without checking expected coverage.

    **Fix:** Gate on equality of shared eligible cells when expected delta is zero, and verify complete expected comparison coverage. Handle eligibility gains separately.

13. **MAJOR — Stage identities and denominator freeze can pass with missing evidence.**  
    [assess_rerun_changes.py:380](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:380), [assess_rerun_changes.py:558](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:558)

    `total_possible` comparisons silently skip gages absent from either table. Gate 5e intersects the expected-identical list with `common`, so an entirely missing expected cohort passes with zero comparisons.

    The tool also does not implement the four-stage interaction residual or inferred-cell labelling required by §10.4.

    **Fix:** Require every expected gage’s evidence, validate stage manifests, and explicitly implement the stage comparisons, interaction residual and measured/inferred labels.

14. **MINOR — B4 plausibility uses the patched range, making the diagnostic self-referential.**  
    [splice_nwis_qualifiers.py:852](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:852)

    Restored extremes expand the range used to assess themselves. The report flags **one** range violation; applying the tool’s same range formula to the **base** series identifies **105**. Neighbours are selected by row position rather than calendar-day distance, and ratio quantiles come only from restored-day ratios.

    **Fix:** Use the February range and an explicit baseline neighbour-ratio distribution, with calendar-date neighbours and sign-correct range expansion. Regenerate the diagnostic disposition.

The core splice acceptance and metadata rules are sound on the inspected data. Absent records/dates, invalid numeric values and sentinels remain unavailable; ambiguous primary series, ambiguous method blocks, duplicate source dates, site mismatches and incorrect units are refused. The three observed multi-method sites resolve to their unqualified blocks. Bounds remain values, negatives remain permitted, and flags remain unchanged.

The oracle’s calendar, leap-year, boundary-gap and qualification rules matched the supplied census. Its restored ledger rows must exist and target missing Q. Remaining hardening: reject duplicate **base** dates instead of overwriting them, and reconcile unavailable ledger entries too. Census agreement validates year decisions; it does not validate interpolation magnitudes or source correctness.

The census’s revised SQL handles NULL/NaN equality correctly, and its **INCOMPLETE** verdict is an improvement. The supplied `census2` remains INCOMPLETE because 5d, 5e and 5f were not run. Its state columns cover the requested categories, but findings 2–4 prevent trusting their attribution. The pinned-header blocks are value-preserving when all signature columns already exist; `--require-clean-tree` correctly refuses the present dirty tree. I did not rerun the filesystem-writing batch self-test.

**Phase B outputs: GO WITH CHANGES** — retain the verified data; repair and rerun the affected audit checks and diagnostics.  
**Starting Phase D on this tooling: NO-GO.**

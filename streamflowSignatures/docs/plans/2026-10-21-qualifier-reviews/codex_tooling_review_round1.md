# Codex tooling review, round 1 (2026-10-26, gpt-6-astra, high reasoning, read-only) — Phase B / C / B′ / E8 tooling before the full Phase B chain

Brief: `scratchpad/codex/review_brief_B_C_tools.md` of the session (the tools, the 40-gage and 5-gage artefacts, six questions). Verdict NO-GO; dispositions in the plan's §11 progress log. Verbatim report:

145,947
Read-only review completed. I ran the full oracle replay, the 5-gage audit, and adversarial checks with writes intercepted in memory. No artefacts were modified.

1. **BLOCKER — B4 rejects correct restorations after CSV parsing.**  
   [splice_nwis_qualifiers.py:642](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:642), [line 710](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:710).

   The supplied ledger produces **1,389 mismatches among 1,391 restorations**. Reading it with `float_precision="round_trip"` produces **zero**. Default pandas parsing changes the last bits of `conversion`/`Q_new`, which B4 then compares exactly. The maximum difference against independently recomputed conversion was approximately `1.78e-15`.

   **Fix:** preserve float round trips, recompute conversion from metadata, and compare against the actual patched Q.

2. **BLOCKER — B6 incorrectly requires every updated metadata row to change.**  
   [splice_nwis_qualifiers.py:751](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:751).

   The supplied splice restores five gages, but only one needs different compile-year totals. B6 therefore fails: **one changed row versus five restored gages**. Updating a row correctly can leave its values unchanged.

   Moreover, B6 neither recomputes those totals nor verifies that every parquet gage joins metadata; line 750 leaves the latter unimplemented. The splice never clears `error_message`, and B6’s allowed-column list would reject that required change.

   **Fix:** compare each restored gage against independently recomputed expected metadata; require unrelated rows to remain identical; verify complete normalized joins; clear and permit `error_message`.

3. **BLOCKER — B4/B5/B8 never establish source → ledger → output value identity.**  
   [splice_nwis_qualifiers.py:704](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:704), [line 729](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:729), [line 755](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:755).

   B4 checks cache values against ledger values using the **ledger’s conversion**. B5 checks only where Q changed. B8 checks membership in the restored-key set. Consequently, an arbitrary finite output value on every ledger key can satisfy all three checks once finding 1 is fixed.

   B4 passes vacuously with no restored rows. B8 can pass with no eligible values. B1 also trusts report rows without reopening ordinary single-primary caches, so a stale report weakens the protection against unavailable caches.

   **Fix:** perform the required keyed source/metadata/base/ledger/output join, checking unique keys, actual base flags and missingness, expected conversion, exact restored values, and unchanged off-target cells.

4. **MAJOR — Splice does not enforce its declared acceptance rule.**  
   [splice_nwis_qualifiers.py:478](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:478).

   Any finite parsed value becomes `restored`, regardless of today’s status. A numeric value with no qualifier or a non-A/P status is therefore restored although `R/usgs_qualifiers.R` rejects it. Units are recorded but never required to be `ft3/s`; source site identity is not checked. Target eligibility comes from the supplied masked file’s `category`, without revalidating the base row’s flag.

   Duplicate source dates overwrite earlier values in `load_today`; a fresh B1 report detects duplicates later, but splice itself proceeds.

   **Fix:** centralize the A/P + finite + sentinel predicate, validate site/statistic/unit and actual base eligibility, and refuse duplicate or ambiguous source records before producing output.

5. **MAJOR — Fetch resume ignores the requested window.**  
   [splice_nwis_qualifiers.py:287](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:287).

   Cache existence plus `ok` is sufficient to skip retrieval. Reusing the September-30 cache with the corrected October-1 command silently retains the shorter responses. Changing the default does not invalidate them. The separate cache directory being used operationally avoids this particular occurrence, but the tool does not enforce it.

   **Fix:** bind cached responses to the complete request URL/parameters and a response hash. Refetch on mismatch; verify request-window provenance during audit.

6. **MAJOR — The oracle trusts incomplete or invalid ledgers, despite matching the baseline census.**  
   [qualifier_year_oracle.py:245](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/qualifier_year_oracle.py:245), [line 256](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/qualifier_year_oracle.py:256).

   Missing dates are rejected **only for gages encountered in the parquet**. An in-memory ledger containing one restored row for a nonexistent gage completed successfully and silently applied nothing. `masked_splice` is computed but unused; existing ledger targets need not be missing, eligible, unique, or finite. Lost years are printed without causing failure.

   The replay also treats ±Inf as observations because `year_status` uses `np.isnan`; Julia normalizes nonfinite values to missing. Four consecutive infinities yielded oracle **valid**, versus preprocessor **gap rejection**. Hard-coded rules are not checked against or hashed with the runtime configuration.

   **Fix:** require complete ledger consumption and target validation, normalize nonfinite values, bind configuration, and fail on lost years. Baseline census agreement cannot establish these unexercised cases.

7. **BLOCKER — Census explanations are permissions inferred from labels, without the required evidence checks.**  
   [assess_rerun_changes.py:196](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:196).

   Emission changes are accepted merely because the registry names a recession/Q–P gate; before/after gate counts are never recomputed. Inserted-row years broadly permit climate changes without checking variable-specific daily changes. Record-dependent labels permit distant changes without verifying the relevant record inputs. Unknown signatures can fall through into the coupled branch.

   **Demonstrated false pass:** on the supplied subset, I added **123** to `15052000`, WY 2020, `temp_max_annual`, and **123** to its `prcp_annual_mean`. The census returned **ALL GATES PASS**, including 5b.

   **Fix:** require registry coverage, attach concrete ledger/state/gate evidence to explanations, and keep numerical validation separate and mandatory.

8. **MAJOR — Classes and “valid-year identity” do not implement the specified state-based contract.**  
   [assess_rerun_changes.py:177](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:177), [line 267](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:267).

   There is no state-table input. A-J membership additionally requires inserted rows to intersect years represented in annual output, rather than implementing the plan’s ledger/state union. Q-valid years are inferred from **any emitted annual signature**, conflating eligibility with collection/emission.

   The oracle gate checks that expected gains are present and `num_water_years` matches. It does not compare complete final year sets or start/end years. A gain explicitly labelled **“NOT in the oracle”** is nevertheless treated as explained.

   **Fix:** derive classes from all three ledgers plus state changes; compare actual Q-valid sets and derived metadata exactly against the oracle.

9. **BLOCKER — “ALL GATES PASS” does not mean §10.5 gates 5b–5g passed.**  
   [assess_rerun_changes.py:342](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:342).

   Implemented 5b checks only **changed, finite `prcp_annual` cells**. It correctly reconciled the supplied 44 changes, but omits seasonal precipitation and every GDD total. It never checks unchanged annual cells against nonzero ledger sums; missing values are replaced with zero without establishing shared eligibility.

   Remaining gaps: **5c** absolute recomputation is absent; **5d** eligibility/emission attribution is unverified; **5e** executed-stage identities, inferred-cell labels and interaction residuals are absent; **5f** lacks denominator and exact valid-year identity checks; **5g** has no required E7 result linkage. Omitting the delta ledger still permits an overall PASS.

   **Fix:** represent required gates as PASS/FAIL/NOT RUN and refuse overall acceptance when incomplete. Implement scoped totals/GDD checks and require independently validated recomputation and stage evidence.

10. **MAJOR — Delta attribution can certify corruption as interpolation.**  
    [preprocess_state_ledger.py:197](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/preprocess_state_ledger.py:197).

    “Secondary interpolation” requires only a nearby raw change. It does not establish that the changed date was interpolated or that its interpolation endpoints explain the value.

    **Demonstrated:** a finite observed Q remained raw `1`, but its preprocessed value changed from `1` to `999` beside a splice. The tool labelled it secondary interpolation and exited **0**. Direct splice/inserted attribution likewise checks key membership, not value agreement.

    The daily comparison uses an inner join; unmatched keys are counted but do not fail acceptance.

    **Fix:** verify interpolation masks, actual bracketing observations and recomputed values; validate direct changes against source ledgers; fail on unexpected missing or duplicate daily keys.

11. **MAJOR — State outputs are insufficient for the promised reconstruction.**  
    [preprocess_state_ledger.py:116](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/preprocess_state_ledger.py:116), [line 125](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/preprocess_state_ledger.py:125).

    State tables contain Q masks/runs, seasonal flags and validity, but lack corresponding climate interpolation masks and rejection diagnostics. Gage tables contain whole-record counts, without window-specific `total_possible`, fraction or qualification decisions.

    The advertised “NaN where invalid for that variable” representation is also false: `_pre` extraction checks only Q validity. The supplied S0 daily file has **5,480 finite `PPT_pre` rows in PPT-invalid years**, reflecting retained raw climate columns in the library output.

    **Fix:** distinguish raw, cleaned and analytically eligible values explicitly; export the missing variable-specific states and window qualification metadata.

12. **MAJOR — Resume fingerprints and output checks permit stale batches.**  
    [run_batched_julia.py:194](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:194), [line 405](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:405), [line 342](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:342).

    The fingerprint omits GAGES-II/HYDAT file contents. An environment-supplied `STREAMFLOW_QA_SCHEMA` is hashed only when supplied through `--qa-schema`. Dirty source changes collapse to the same revision plus `dirty=True`, so editing an included runner/helper can leave the fingerprint unchanged.

    `outputs_ok` checks filenames and, for fresh runs, CSV/timing mtimes—not output identity. Annual output is optional and its freshness is unchecked.

    `--require-clean-tree` correctly rejects an initially dirty/unknown checkout, including untracked files. It checks once before the batch loop and does not protect external inputs or subsequent tree changes.

    **Fix:** hash all resolved dependencies and outputs; prohibit dirty-tree resume or hash its contents; require expected annual output; revalidate the execution snapshot before each batch.

13. **MAJOR — Qualification audit still rejects legitimate frozen-cohort exclusions.**  
    [audit_qualification.jl:73](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/audit_qualification.jl:73), [line 107](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/audit_qualification.jl:107).

    The original sample check still requires `qualified == emitted`. A sampled gage that qualifies but is intentionally outside the selection fails before the new cohort contract can help. Conversely, omitting the optional census skips the complete `selected ⊆ qualified` check.

    **Fix:** in cohort mode, independently check qualification, then compare emission against `qualified ∩ selected`. Require complete qualification evidence for every selected gage.

14. **MAJOR — The newly added E7 summary checker permits forbidden gate-state waivers.**  
    [check_summary_equality.py:99](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/check_summary_equality.py:99).

    This file appeared during review. `--waive-column` skips missingness-state differences, contrary to the plan’s unconditional gate-state requirement. A test with reference `Q_mean=1`, candidate `Q_mean=NaN`, and a column waiver returned **SUMMARY GATE PASS**, while reporting one gate-state difference.

    Its automatic tie exemption is also broad: any reference near-tie plus any annual difference exempts the pair’s designated statistics, without requiring the annual differences to satisfy port tolerances.

    **Fix:** gate missingness independently of waivers and restrict tie exceptions to validated, explicitly bounded residuals.

15. **MINOR — Several promised diagnostics and reporting checks are missing or altered.**  
    [splice_nwis_qualifiers.py:606](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:606), [line 676](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:676), [assess_rerun_changes.py:295](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:295).

    B2 adds an unapproved aggregate 95% threshold and does not require control coverage for every relevant gage. The control’s rule-effect count is restricted to February-missing rows and lacks year-level results. B3 verifies aggregate target counts, but collapses sentinel/non-numeric cases into `no_value` and drops today’s status from the written aggregation. B4’s required plausibility diagnostics/dispositions are absent.

    The magnitude census is partial: annual per-signature magnitude summaries, standardized per-gage effects, regional/seasonal views, interpolation exposure and population-mean comparisons are not completed.

    **Fix:** implement these explicitly; distinguish diagnostics, incomplete reporting and acceptance gates.

16. **NOTE — Several important parts did verify correctly.**  
    [qualifier_year_oracle.py:102](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/qualifier_year_oracle.py:102), [splice_nwis_qualifiers.py:500](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:500), [run_julia_benchmark.jl:613](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark.jl:613).

    The full replay reproduced **322,563 year statuses and 16,028 window decisions with zero mismatches**. Leap-year length, October–December assignment, boundary rejection and negative-Q acceptance checked correctly for finite inputs.

    The supplied patched parquet passed the full B5 scan: identical schema, row groups, row count, row order and non-Q columns; exactly **1,391 changed Q cells**, all on ledger keys. Conversion and compile-year arithmetic match R. Missing records/dates/non-numeric values remain unavailable; the three named multi-method caches resolve to the unqualified block.

    DuckDB’s implemented NULL/NaN comparison behaves correctly. The QA blocks add nothing when all columns exist; the supplied left versus left-with-QA runs had **zero changed cells** across all 40 gages. Merge’s pinned-header comparison is strict, though optional.

**NO-GO for Phase B’s full chain and Phase D on the current tooling.** The immediate B4/B6 false failures must be fixed, and the audit, resume and census false-pass paths must be closed before production acceptance is credible.

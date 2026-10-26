# Codex review of Phases A and A′ (2026-10-26, gpt-6-astra, high reasoning, read-only)

Brief: `scratchpad/codex/review_brief_A_Aprime.md` of the session. Verdict NO-GO; every finding applied the same day (CHANGELOG [October 2026]). Verbatim report:

101,203
**NO-GO for tagging the current tree.** The ordinary fixtures pass, but adversarial probes expose defects and a test exemption that accepts a real gate-state error.

Verification: Julia’s **26 checks passed** when evaluated in memory with artifact writes omitted; Python’s **32 checks passed** against the existing fixture. R could not start because the read-only sandbox prevents creating `R_TempDir`. Independently comparing the saved 40-gage CSVs confirmed **0/2,276 changed columns in left**, and **599 in outer**. I changed no files.

1. **MAJOR — Duplicate climate dates can change Q-only results.**  
   [climate_frame.jl:69](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/climate_frame.jl:69), [climate_frame.py:68](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/climate_frame.py:68)

   Every candidate climate row survives when its date is absent from streamflow; duplicate candidates are appended independently. Julia and Python subsequently duplicate that day in the preprocessor’s daily grid.

   **Evidence:** Duplicating the fixture’s climate row for 1994-11-14 changes Python’s WY 1995 raw missing-Q count from **5 in left to 6 in outer**. A separate Julia probe changes a 365-day cleaned year to **366 rows**, with missing/interpolated Q counts increasing from 1 to 2. R instead drops duplicate dates in its preprocessor, creating another cross-language difference. The real-data subset contains no duplicate keys, so its successful comparison does not cover this.

   **Fix:** Establish a duplicate-date policy before appending candidates: collapse identical duplicates and reject conflicting observations, or fail explicitly on duplicates. Test duplicates on both existing and added dates. Preserve the historical `left` behavior.

2. **MAJOR — Outer joins and absence checks use different date representations.**  
   [climate_frame.jl:54](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/climate_frame.jl:54), [climate_frame.py:50](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/climate_frame.py:50), [climate_frame.R:41](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/climate_frame.R:41)

   The initial join uses unmodified keys; the later span and membership checks convert them to calendar dates.

   **Evidence:** With Julia streamflow dates `DateTime(2001-01-01)` and `DateTime(2001-01-03)`, and climate `Date` rows for January 1–3, outer produces PPT values **missing, missing, 20**. It recognizes January 1 and 3 as existing dates and therefore never recovers their climate. Python likewise leaves climate missing when streamflow timestamps are noon and climate timestamps midnight. R follows the same inconsistent normalization sequence, though I could not execute it.

   **Fix:** For `outer`, normalize copies of both join keys consistently before joining and checking membership. Define whether non-midnight/time-zone-bearing timestamps are supported or rejected. Keep `left` unchanged for attribution, and test mixed Date/DateTime inputs explicitly.

3. **MAJOR — The tie exemption can hide arbitrary summary errors, including gate-state divergence.**  
   [selftest_climate_frame.py:200](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/selftest_climate_frame.py:200), [selftest_climate_frame.R:146](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/selftest_climate_frame.R:146)

   Once a base series qualifies as “tie-sensitive,” both numeric discrepancies and finite-versus-NaN discrepancies in its selected summaries become diagnostics rather than failures. The test does not establish that changed ranks actually explain those discrepancies. Python additionally uses ordinary equality when detecting annual differences, so matching NaNs can falsely establish “not bit-identical.”

   **Evidence:** I replaced Python’s finite `Q95_Q10_mk_pval` with NaN in memory, leaving annual values unchanged. The self-test still reported **32 passed, 0 failed; exit 0**.

   **Fix:** Use NaN-aware annual equality. Always gate output presence and eligibility, with narrowly documented exceptions for mathematically explained undefined statistics. Validate permitted rank residuals by recomputing summaries from a common annual series or proving the actual rank change; do not exempt every selected summary of a qualifying base.

4. **MAJOR — Required window/denominator coverage is absent.**  
   [selftest_climate_frame.jl:57](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/selftest_climate_frame.jl:57), [selftest_climate_frame.jl:80](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/selftest_climate_frame.jl:80)

   The fixture truncates the **record**, but never applies a runner analysis window ending inside that record. None of the three self-tests evaluates the runner’s `total_possible` or qualifying-fraction gate. Thus §10.2(ii)’s window-ending-in-an-interior-hole requirement remains untested.

   Other gaps: exact historical left-frame comparison—including types/order—input non-mutation assertions, duplicate dates, mixed date representations, empty frames, boundary-only years, and explicit missing climate values. Python/R check the **number** of added dates against Julia, rather than their exact identities.

   **Fix:** Add the actual window/fraction calculation with the window ending inside the absent WY 2005; assert denominator and inclusion decisions. Add focused helper fixtures for the listed cases and compare `left` directly with the former join expression.

5. **MINOR — Resolved-mode provenance is incomplete at the batch boundary.**  
   [run_batched_julia.py:89](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:89), [run_batched_julia.py:204](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:204), [run_batched_julia.py:405](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:405)

   All three runners correctly record resolved `provenance.climate_frame`, including default `left` when unset. Explicit environment overrides also reach the shared `env_overrides` comparison, so ordinary explicit left-versus-outer batches are protected.

   However, `SHARED_PROVENANCE` omits the resolved field, and `batches.json` records no mode. Merge copies batch 1’s resolved value without checking the others. This falls short of §10.2’s explicit provenance contract.

   **Fix:** Require and compare resolved `climate_frame` during merge. Record the resolved mode in the stage’s batch manifest, including when unset; account explicitly for input batches reused across stages with different modes.

6. **MINOR — R’s legacy branch no longer necessarily keeps its left join.**  
   [run_rpkg_benchmark.R:286](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_rpkg_benchmark.R:286)

   The helper executes **before** `if (use_legacy)`. Consequently, an authorized legacy run with `STREAMFLOW_CLIMATE_FRAME=outer` receives added rows. Julia and Python retain their separate legacy joins. This directly violates the stated legacy exception.

   **Fix:** Pass `"left"` when `use_legacy`, or retain a separate historical merge for that branch. Record the effective mode accurately.

7. **MINOR — R leaves existing `year` and `doy` columns missing on added rows.**  
   [climate_frame.R:62](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/climate_frame.R:62)

   Julia and Python derive these columns when present. R derives only `water_year`, `month`, and `dowy`, then fills all remaining columns with NA. Its self-test removes `year`, `doy`, and `flag` before calling the helper, hiding this difference.

   **Impact:** The current R runner selects only `gage_id`, `Date`, and `Q`, so this does not presently change its production signatures. It nevertheless violates the requested uniform temporal-column behavior.

   **Fix:** Populate `year` and `doy` from the added dates when those columns exist, and test the full input schema.

8. **NOTE — Left-mode parity and the narrowed Julia lookup look correct, within the intended exception.**  
   [climate_frame.jl:50](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/climate_frame.jl:50), [run_julia_benchmark.jl:230](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark.jl:230)

   Each left branch invokes the exact former join expression, preserving its ordering, column ordering and dtype behavior. Python’s exact frame comparison passed; Julia’s present/absent tuple-key lookup behaved correctly.

   The former catch swallowed missing-column and join errors as well as absent-gage lookups; those errors now propagate, as required. I found no converse case where the new lookup hides an error the old code exposed.

   For compatible, unique daily keys, the span **and original-WY-set** restriction is correct: partial years remain bounded, boundary rows establish membership, missing climate values are retained, and empty frames return before extrema calculations. Julia/Python input non-mutation probes passed; R copies candidates before its in-place temporal update. Missing Q is handled downstream, and the preserved WY set protects `all_wy`. Library directories are unchanged.

   **Fix:** Capture these properties in regression tests, especially the lookup-error distinction and full frame equality.

9. **NOTE — No Phase A acceptance-rule defect found for the specified inputs; integration coverage is missing.**  
   [usgs_qualifiers.R:30](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/R/usgs_qualifiers.R:30), [helperFunctions.R:652](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/R/helperFunctions.R:652)

   Inspection confirms whitespace/colon/comma token separation, arbitrary trailing remarks, rejection of lowercase and NA status, numeric-string parsing, and rejection of string `"-999999"`, NaN and infinities. The mask changes only `Q_rawUnits`; the original `flag` column and conversion formula remain unchanged.

   Sourcing works for the documented Windows ingestion invocation from the repository root. It is **cwd-dependent**, not independently robust when sourcing the helper by absolute path elsewhere. Existing ingestion sources already require that cwd, so this is not a new Windows-specific failure.

   **Fix:** Add comma/string-sentinel cases and a stubbed `generate_streamflow_dt` integration test asserting exact flag preservation and converted Q. If standalone helper sourcing is supported, resolve the qualifier file relative to the sourced script.

**NO-GO for committing Phases A and A′ as the tagged rerun tree until the outer-frame defects and acceptance-test gaps are addressed.**

# Codex plan review, round 1 (2026-10-23) — model gpt-6-astra, reasoning high, read-only

Brief: the plan of 2026-10-21 (`2026-10-21-qualifier-mask-fix-rerun-plan.md`, version 1), its oracle CSV, the census §10, the ingestion code, the Julia preprocessor, the harness and the local data. Verdict and findings as returned:

**NO-GO as written.** The proposed oracle is demonstrably wrong, the provisional-value premise is contradicted by the cached responses, a qualifying candidate is omitted, and appending candidate metadata would create duplicate keys. Correct those findings, separate qualifier policy from retrieval-vintage changes, and repair the acceptance and staging procedures before approving production regeneration. All analyses below were read-only and offline.

1. **BLOCKER — The oracle reproduces faulty arithmetic, not the canonical qualification rules.**

   Using `.venv/bin/python -B -c …`, I independently combined `qualification_census.csv` with `AR_gage_wy_budget.csv`: **6,707 / 6,278**, identical admission lists and CSV counts.

   However, [min_flow_rule_footprint.py:85](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/min_flow_rule_footprint.py:85) counts existing rows only. [io.jl:281](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/io.jl:281) constructs the complete October–September calendar, including leap days. It rejects **raw** NA counts above 30 **before interpolation**, then rejects residual boundary NAs.

   My batch-read, full-calendar replay matched every relevant baseline census status. Assuming all targeted approved days become finite:

   | Restoration policy; existing compiled gages only | WY 1993–2025 | WY 1980–2025 |
   |---|---:|---:|
   | Plan’s `A R` budget | 6,707 | 6,278 |
   | Calendar-correct `A R` | **6,707** | **6,276** |
   | Calendar-correct four approved codes | **6,709** | **6,279** |

   The budget has **50 false-positive and 33 missed restored years**: 1,197 recoverable `A R` years across 601 gages, not 1,214/613. Product #1 replaces `02313700` with `02312722`; product #2 loses predicted additions `03250322` and `07076517`. The window-start denominator is correctly reproduced; it is not necessarily 33/46 for early-ending records.

   Original-oracle zero-slack additions:

   - Short window: `01196561,01454700,01480638,01589290,02313700,02317620,03292475,03508050,04180988,072632962,07325850,08397600,12210700`.
   - Long window: `10308800,01349810,02146530,02153200,02483500,03028500,06131200,06935770,07076517,07363400,09301500,09419800,10243260,10251300,11180900,11261500`.

   **Change:** rebuild the oracle on full calendars and actual recoverable values; include exact year sets, candidates and every accepted code. “Every deviation explained” must not override an allegedly exact gate.

2. **BLOCKER — “Provisional-remark days carry no value today” is false even in the samples.**

   I joined cached NWIS values to **original parquet gage/date keys**, rather than searching current qualifier strings. Sampled formerly masked days with numeric values now include:

   | Original flag | Numeric days now |
   |---|---:|
   | `P Rat` | 703 |
   | `P Bkw` | 199 |
   | `P Ice` | 136 |
   | `P ZFL` | 86 |
   | `P [4]` | 74 |
   | `P <` / `P >` / `P Ssn` | 66 / 4 / 5 |

   Total: **1,273**. Examples: `03282120`’s `P Rat` days now carry `P e`/`P`; `08286500`’s `P Ice` days now carry `A`/`A e`. Some other sampled dates are absent from today’s response, which is different from an explicit no-data value.

   **Change:** audit every targeted historical gage/date against the cached current response, recording value, current status, units, absence and sentinel separately. Current availability cannot establish February availability. Neither a zero-value assumption nor the present oracle survives this finding.

3. **BLOCKER — The candidate list and metadata operation are wrong.**

   Independent cached-JSON replay gives `11129800` **20 qualifying calendar years** without masking, versus 18 masked. The finding CSV misclassifies it as “zero flow.”

   Running `~/.juliaup/bin/julia --startup-file=no --compiled-modules=existing --project=julia -e …` through canonical `preprocess_daily_data` confirmed:

   - Four listed candidates: **22, 25, 26, 25** valid shorter-window years; all qualify for product #1, none for #2.
   - `11129800`: **29 / 33** valid years; qualifies for both.

   All five lack a site in the 8,017-site Daymet parquet, verified by batch-reading `site_id`.

   The staged metadata already contains **16,994 rows**, including all five as `insufficient_years`. Appending rows creates duplicate normalized IDs in the runner’s metadata join.

   **Change:** identify candidates by counterfactual gates, not `removed_by`; incorporate these admissions into the oracle; **update existing metadata rows**, enforce unique normalized IDs, and clear stale failure messages.

4. **MAJOR — The footprint confuses approved-code gages with the broader masked population.**

   A complete `iter_batches(batch_size=131072)` tally over **111,624,189 rows** confirms every approved-code table row: `A R` 181,461/755; `A <` 11,444/139; `A [4]` 419/101; `A >` 218/30. Their per-product and outside-both counts also match individually.

   Their **union is 987 gages**, not 1,186: **908 / 757** in the products and **63** outside both, not 1,088/898/81. Provisional remarks total **15,257 days across 263 gages**. There are also **12,103 US rows with null flags**. No additional non-null US flag category was found.

   **Change:** distinguish historical masked days, distinct gages, current recoverable values and causal mask losses. Correct CHANGELOG/census’s stale 15,460 provisional-day figure too. The 193,542 count establishes approved-code NAs, not agency-value availability across every affected date.

5. **MAJOR — Climate attribution and its proposed tolerance are incorrect.**

   The mask changes Q but preserves rows. The runner joins climate onto those rows at [run_julia_benchmark.jl:220](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark.jl:220). Thus restoring Q does **not** restore Daymet join dates.

   The actual interaction is whole-year selection: `io.jl:444` skips rejected Q years before climate processing. Newly valid Q years can therefore add climate years. Existing absent-date climate biases remain.

   **Change:** require exact agreement for year-local, climate-only annual values on shared valid years, including years with restored Q. Report newly admitted climate years separately from changed Q-dependent metrics and record-dependent summaries. Drop “few percent on a handful” and the claim that this splice shrinks the absent-date defect. Keeping decision 5 separate is reasonable only with this corrected attribution.

6. **MAJOR — Accepting every numeric A/P value bundles distinct scientific policies.**

   Cached qualifier tables explicitly define `<`/`>` as bounds; `[4]` means no matching qualifier, not established measurement quality. Treating a reporting limit as an exact observation affects zero-flow frequencies, low-flow quantiles, drought thresholds and recessions. Approval does not remove censoring.

   Existing code already accepts estimated `A e`/`P e`; HYDAT symbols are retained without masking (`R/helperFunctions.R:662–687`). That supports retaining estimates with provenance, but does not justify interpreting every remark identically.

   **Change:** separate storage policy from analytical treatment. I would restore revised values and retain agency estimates, while requiring an explicit bound-substitution policy and sensitivity assessment for censored values. Unknown remarks need documented handling. Require finite values, recognized status and sentinel rejection. B4’s blanket prohibition on new negatives also conflicts with “every numeric value” and configured negative-flow acceptance.

7. **MAJOR — Splice-only cannot provide the claimed pure attribution.**

   October values may reflect revisions unavailable in February. Preserving old flags while inserting newer values can also misstate their status. Two retrieval dates in a README are insufficient provenance.

   A missed third option is **re-ingest complete records only for affected USGS gages**, keeping unaffected USGS/HYDAT records fixed. This broadens class-A changes but avoids mixing retrieval vintages within those gages.

   **Change:** compare that option explicitly. A splice remains publishable as a *patched historical snapshot* if it deposits the immutable base, raw responses and a gage/date correction ledger containing old/current qualifiers, values, conversion and retrieval details. Do not describe it as reconstructed February agency data. A current-source old-rule/new-rule control would separate rule effects from vintage effects.

8. **MAJOR — Input gates can reject legitimate source changes and miss wrong restored values.**

   B2 already fails on cached examples: `08286500` has **711/16,113 = 4.413%** differing overlapping values; `07362000` has **192/11,460 = 1.675%**. These exceed 0.1%. Calling every mismatch a revision is not verification.

   B5 restricts *where* changes occur but does not prove the restored values are correct. B4’s range/neighbor checks can accept a wrong splice or reject legitimate revised extremes.

   **Change:** independently verify every restored value against its exact cached site/date/statistic/unit record and conversion. Separate conversion checks from source-difference reporting. Enforce unique daily keys, target coverage, metadata uniqueness and candidate calendar semantics. Make plausibility checks diagnostics with explicit dispositions. Replace “approved NA must be zero” with “no eligible finite source value was discarded”; approved source no-data may legitimately remain NA.

9. **MAJOR — Class A is incomplete, and its comparison is mostly diagnostic.**

   The 552/447 counts describe gages gaining years, not all affected gages. Historical approved-code masks occur inside the windows at **899 / 757 existing product gages**; `A R` alone affects **683 / 562**. Filling previously interpolated gaps can change signatures without gaining a year.

   **Change:** derive classes from actual changed input keys; make them exhaustive and disjoint. Gate exact valid-year identities, unchanged metadata and unaffected annual keys. Use an explicit dependency registry for year-local versus record-dependent signatures. R² reporting alone cannot reject a broken affected-gage result. The 0.766 population fraction for new climate columns is not a defensible minimum for every family at newly admitted, sometimes climate-free gages.

10. **MAJOR — Production and provenance gates need implementation changes before launch.**

   Verified code conflicts:

   - `run_batched_julia.py:453–457` rejects additions with `--reference` unless `--allow-missing` is supplied.
   - `validate_production_run.py:69–75` gates overlapping coordinates, **not gage-set equality**.
   - `run_rpkg_acceptance_gates.sh:23–26,129–132` hard-codes the August reference and **1,653 columns**.
   - Batch resume checks command/environment, not changed code/input content (`run_batched_julia.py:318`).
   - Merge permits all batches to share `git_working_tree_dirty=true`.

   Both v2 timing files actually record dirty trees; the saved status lists untracked `julia/src/climate.jl`, which ordinary saved diffs do not preserve.

   **Change:** parameterize the harness, enforce the exact admission set, bind resumes to content hashes, and fail production on dirty/unknown provenance. Resolve the v2 source baseline before treating every discrepancy as a splice error.

11. **MINOR — Class-U exact equality is feasible, but “bitwise” needs definition.**

   Adding candidates changes sorted-rank batch assignments (`run_batched_julia.py:203`). That alone does not change per-gage calculations. Equal batch headers and equal v2/v3 headers stabilize the high-NA denominator; annual merge explicitly sorts by `(gage_id,signature,water_year)`.

   **Change:** compare keyed cell values, not parquet bytes, row offsets or row-group layout. Specify NaN, null and signed-zero semantics; existing `compare_climate_rerun.py:89–103` is value equality, not literal bit equality. Pin runtime/configuration and establish a baseline control. I found no inevitable numerical class-U change caused solely by reassignment under those conditions.

12. **MAJOR — Staging unchanged would publish false statements.**

   [stage_hydroshare_climate_v2.py:185](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/stage_hydroshare_climate_v2.py:185) hard-codes “same streamflow input,” “gage set is unchanged,” and “Q-only columns identical.” Lines 145–148 claim every gage has climate. Lines 85–94 copy old climate-review filenames and silently skip missing files.

   **Change:** adapt staging for the qualifier release, require its audit/oracle/comparison artifacts, and validate rendered claims. Distinguish internal build IDs from public release versions in one manifest. Include correction provenance in R3 and update R4/R5 coverage relationships even if their data remain unchanged.

13. **MINOR — The documentation inventory is incomplete; blanket replacement would corrupt history.**

   `rg` for the four counts and qualifier rule found omitted current documentation or procedures:

   - `docs/SIGNATURE_GUIDELINES.md:26,540`; `docs/WORKFLOW_REVIEW.md:41`.
   - `docs/CROSS_LANGUAGE_STATUS.md`; `rpkg/README.md`; `python/README.md`; `.claude/rules/signatures-code.md`.
   - `docs/plans/2026-09-04-dataset-join-guidance.md`, including R4/R5 intersections and “attributes for every product gage.”
   - `EO_data_processing/{CLAUDE.md,README.md,README_DAYMET.md,geometry/README.md}`: clarify the frozen geometry universe versus expanded streamflow coverage.
   - `docs/filtering_min_flow_rule_2026-10-21.csv`, the filtering-stage census CSV, reconciliation logs and the manuscript’s Figure 3 counts/underlying figures.
   - `summarize_qualification_census.py:45`, staging/harness text and wrapper descriptions.

   **Change:** create a version-aware update checklist. Preserve dated benchmark results, archives, earlier plans and unchanged **7,964-polygon** counts; add supersession pointers instead of rewriting every occurrence.

14. **MAJOR — The schedule excludes substantial implementation and review work.**

   Phase B budgets 1.5 hands-on hours for a new downloader, splice engine, compiler and seven audits. The ten-minute oracle estimate is contradicted by the required calendar and candidate corrections. Existing census/audit scripts load the full streamflow table; annual validators load whole parquets, and the comparison loads both.

   **Change:** put compiler validation, oracle correction, bounded-memory QA, tool adaptation and an input-review checkpoint **before** production. Cached service responses establish neither future availability nor a complete OGC fallback; that needs pagination, series, units and qualifier coverage checks. Keep production and rpkg sequential, budget reruns after review fixes, and treat two–three days as an unvalidated estimate.

**Verified independently**

- Full US flag/NA tally, distinct-gage unions and product intersections.
- Oracle arithmetic; full-calendar counterfactual qualification and baseline-status agreement.
- Cached historical-date/current-value matches and conversion-overlap discrepancies.
- Candidate calendar gates; canonical Julia qualification; existing metadata rows; Daymet absence.
- Runner, validator, staging behavior and repository documentation occurrences.

**Taken on trust**

- February retrieval provenance and byte identity of the separate compiled-file copies.
- Availability of historical raw agency values beyond supplied caches.
- Future NWIS/OGC behavior, retirement timing and download duration.
- Polygon absence for the candidates, historical review verdicts and prospective production runtimes.

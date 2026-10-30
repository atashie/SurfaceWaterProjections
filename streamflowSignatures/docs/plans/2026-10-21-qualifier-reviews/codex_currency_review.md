# Codex review — code and documentation currency after the default flip (2026-10-29, gpt-6-astra, read-only sandbox, reasoning high)

Brief: `review_brief_currency.md` (same folder) with a whole-repository assignment. Verdict GO WITH CHANGES, no BLOCKER; every finding applied 2026-10-29 (dispositions in `../2026-10-21-qualifier-mask-fix-rerun-plan.md` §11, row "currency review"). Its sandbox had no writable temp directory and no network, so the R suites, the batched self-test and the mirror diff were blocked there (all run by the Sonnet reviewers and by the main session). Verbatim report:

| Rank | Finding and evidence | One-line fix |
|---|---|---|
| **MAJOR** | **Batch provenance can record the wrong effective frame.** [run_batched_julia.py:465](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:465) records `outer (default)` when the environment is unset. A historical wrapper subsequently selects `left`, e.g. [run_julia_benchmark_drought_1993_2025_60pct.jl:46](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark_drought_1993_2025_60pct.jl:46). Thus the batch fingerprint can disagree with the child’s timing provenance. | Record requested and effective modes separately; obtain the effective mode from the child’s provenance. |
| **MAJOR** | **Current documentation still says the climate fix is on hold.** [docs/SIGNATURES.md:779](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/SIGNATURES.md:779), [docs/DEVELOPMENT.md:82](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/DEVELOPMENT.md:82), and [claude-skill/streamflow-signatures.md:125](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/claude-skill/streamflow-signatures.md:125) describe the old join as current. [run_batched_julia.py:7](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:7) likewise describes an unconditional left join. `docs/SIGNATURES.md` was not changed in the review range. | Mark the defect as historical and document `outer` as the default since 2026-10-28. |
| **MAJOR** | **Current Daymet status contradicts the products.** [CHANGELOG.md:34](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/CHANGELOG.md:34) says the rerun decision is pending and “No product uses it yet”; [docs/DATA_SOURCES.md:18](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/DATA_SOURCES.md:18) says the reprocessed input is “NOT yet used by either product.” Both product provenance records identify that input. | Replace these pending-status statements with the completed October reruns and outstanding deposit status. |
| **MINOR** | **Optional qualifier integration test contains contradictory assertions.** [R/tests/test_qualifier_codes.R:160](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/R/tests/test_qualifier_codes.R:160) inserts `−999999` on an `A R` row; line 170 requires every approved-remark row to be non-NA, while line 171 correctly requires the sentinel to be NA. This is a static finding; R execution was blocked. | Exclude unusable values from the approved-remark “kept” assertion; retain the sentinel assertion. |
| **MINOR** | **The legacy R runner is incorrectly presented as canonical.** [docs/DEVELOPMENT.md:739](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/DEVELOPMENT.md:739) calls `run_r_benchmark.R` canonical. That runner calls the legacy implementation, whose join remains `all.x=TRUE` at [R/helperFunctions.R:4851](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/R/helperFunctions.R:4851); it does not use the climate-frame switch. | Label this runner deprecated and direct current R workflows to `run_rpkg_benchmark.R`. |
| **MINOR** | **Public evidence references point to excluded material.** [README.md:181](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/README.md:181), [docs/CROSS_LANGUAGE_STATUS.md:35](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/CROSS_LANGUAGE_STATUS.md:35), and both staged run notes at line 79 refer readers to `docs/plans/`, excluded by [publish_to_hisss.sh:29](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/publish_to_hisss.sh:29). The staged notes explicitly say the plan is in the HISSS repository. | Publish a concise verification record outside `docs/plans/` and redirect public references. |
| **MINOR** | **Internal version naming reaches public documents.** [README.md:129](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/README.md:129), line 181, and [python/README.md:263](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/README.md:263) use “version 3”; both staged run-note titles expose “internal build v3.” The expected publish set contains 23 such Markdown lines across nine files. | Use public “version 2” in reader-facing prose; retain internal build identifiers only where needed for provenance. |
| **MINOR** | **Reproduction examples are incomplete or obsolete.** Both `D/resource{1,2}_signatures_*/hisss_run_notes_*.md:42` omit required `--out-dir`; executing the shown arguments produced an argparse error. [.claude/skills/run-benchmark/SKILL.md:63](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/.claude/skills/run-benchmark/SKILL.md:63) still recommends standalone product-wrapper launches, which now hit the guard and do not enforce the frozen cohort. | Supply a complete split/run/merge example with output paths and the frozen selection manifest. |
| **MINOR** | **Python README conflates two exception counts.** [python/README.md:266](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/README.md:266) says 8/7 gate-state exceptions are “of” the 103 waived cells. Both staged summary-gate logs, line 7, report **103 finite-value waivers plus separately 8/7 gate-state exceptions**, consistent with `docs/CROSS_LANGUAGE_STATUS.md:59`. | Change the wording to “103 finite-value cells and, separately, 8/7 gate-state cells.” |
| **MINOR** | **README contains stale current schema descriptions.** [README.md:72](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/README.md:72) calls the full summary 1,653 columns despite the current 2,277-column schema; line 114 says “eight” categories immediately after listing nine. | Update the current-product description to 2,277 columns and nine categories; retain dated historical comparisons. |
| **NOTE** | **Verification remains incomplete because of environment restrictions.** R could not create its temporary directory; batch self-tests could not create scratch files; mirror access failed. Python’s two setup errors also arose from temporary-directory restrictions. | Repeat the blocked checks in an authorized writable scratch environment with public network access. |

**Verdict: GO WITH CHANGES.** I found no demonstrated corruption of the reviewed products. The default flip itself is implemented, but documentation currency and effective-frame provenance need correction. This is not an unconditional test or public-mirror certification.

No repository, run-folder, staged-deposit or worktree files were edited. HEAD was `713fdbfc655dffab8489d0922fb31340811e20a4`. Contrary to the clean-tree premise, one untracked file existed throughout: `docs/plans/2026-10-21-qualifier-reviews/review_brief_currency.md`.

**A — Stale-text sweep**

I searched source, documentation, hidden rules/skills and staged text for the four-code rule, `keep_codes`, `USGS_ACCEPTED_CODES`, streamflow-date climate joins, left-default wording, `left (unset)`, both old row counts, version-3 wording and deposit claims.

The material distinctions are:

- **Old qualifier rule:** `docs/WORKFLOW_REVIEW.md:41` belongs to an explicitly historical February review. `docs/SIGNATURE_GUIDELINES.md:540` remains stale, but is a dated external-document snapshot; `docs/reconciliation/guidelines_todos.md:19` explicitly tracks its correction. These are not evidence that ingestion still uses the old rule.
- **Removed constant:** references remain in `config.R:122`, the changelog and plan history, but as descriptions of its removal. I found no executable reference.
- **Historical analysis tools:** `min_flow_rule_footprint.py:3` and `min_flow_rule_replay.py:34` describe the four-code rule as belonging to `generate_streamflow_dt`. Their historical replay purpose is legitimate; those comments should explicitly say “February 2026.”
- **Left default:** the actual helper defaults are corrected. No remaining literal `left (unset)` was found in the scanned current source. The substantive stale join descriptions are in the findings table.
- **Old annual counts:** 27,602,161 / 34,681,331 occur in dated October 12/13 build descriptions and old-versus-new comparisons. Those historical numbers should remain, with clear superseded-build labeling.
- **Deposit wording:** current product-status passages say the corrected products are not deposited. However, `D/resource3_input_data/hisss_abstract_resource3.txt:1` says “deposited twice,” and `docs/reconciliation/manuscript_log.md:17` says “the deposited products restore…”. These read as publication-ready wording, not evidence of an upload. “Included as two snapshots” and “staged products” would avoid ambiguity before deposit. References to the earlier HydroShare products are a separate, historical claim.

**B — Code and documentation agreement**

The ingestion call chain is intact:

`run_ingest_usgs_hydat.R:4–6` sources configuration, helpers and qualifier functions; its USGS calls at lines 93 and 105 enter `process_gages_rawToRaw`; `R/helperFunctions.R:386` calls `generate_streamflow_dt`; line 650 applies `usgs_qualifier_mask` before unit conversion.

`R/usgs_qualifiers.R:30–47` accepts status token `A` or `P`, requires a finite numeric value and rejects `−999999`. Remark codes do not mask values. The original flag survives selection and output at `R/helperFunctions.R:640` and 431–434. Gage admission remains a separate gate; accepted gages retain all dates in the requested window.

`splice_nwis_qualifiers.py:357–369` implements the corresponding numeric/status acceptance rule. Its **restoration scope is deliberately narrower**: the four approved remark codes, rather than retrospective restoration of every provisional remark.

The production runners resolve and use the helpers as follows:

| Entry point | Resolution | Application |
|---|---|---|
| Julia | `climate_frame.jl:36`; runner line 51 | Runner line 233 |
| Python | `climate_frame.py:38`; runner line 171 | Runner line 330 |
| rpkg | `climate_frame.R:27`; runner line 75 | Runner lines 288–289 |

All three normalize case/whitespace, default to `outer`, accept `left` explicitly and reject unsupported modes. Outer additions are restricted to the streamflow record span and water years containing streamflow rows. The deprecated legacy branches retain left joins.

I independently counted **12 Julia wrappers: 10 `left`, 2 `outer`**. The two outer wrappers are `run_julia_benchmark_prod_{1993,1980}_2025_60pct_climate.jl`. Every other `run_julia_benchmark_*.jl` wrapper selects `left` through `get!`, preserving an explicit environment override.

Both climate wrappers use the corrected input and `streamflow_<window>_60pct_climate_v3_26oct2026` prefix. Both guards fired during read-only evaluation. Their precise protection is: **unset output directory plus an existing, nonempty default product directory**. They are not an unconditional prohibition on explicitly supplying that directory.

The batch runner’s merged output takes the child’s actual frame from timing provenance (`run_batched_julia.py:725`); its earlier batch fingerprint has the unresolved-default defect described above.

The documented rerun CLI flags exist, including `--gage-list`, `--climate-frame`, `--qa-schema`, `--stage`, `--require-clean-tree`, `--gages`, `--allow-diff-signature`, `--waive-column`, `--allow-gate-state-column`, `--accept-waiver` and `--floor`. The demonstrated CLI problem is a **missing required argument**, not an invented flag.

`git diff 54f0483..HEAD -- julia/src python/streamflow_signatures rpkg` was empty. The libraries do not resolve `STREAMFLOW_CLIMATE_FRAME`; that policy belongs to the runners.

**C — Tests**

| Check | Observed result |
|---|---|
| Julia unit suite | Exit 0; **3,184 passes across 13 testsets**, summed from the emitted counts below |
| Python suite | `175 passed, 9 skipped, 145 warnings, 2 errors in 7.80s` |
| Installed rpkg testthat suite | Blocked before tests: `Fatal error: creating temporary file for '-e' failed` |
| R qualifier suite | Blocked before tests: `Fatal error: cannot create 'R_TempDir'` |
| Julia climate self-test | `41 passed, 0 failed` — read-only adaptation described below |
| Python climate self-test | `38 passed, 0 failed` |
| R climate self-test | Blocked before tests: `Fatal error: cannot create 'R_TempDir'` |
| Batch-runner self-test, daymet venv | Blocked at `tempfile.mkdtemp`: `No usable temporary directory found` |

Julia’s emitted pass/total counts were:

```text
StreamflowSignatures.jl                                           78 / 78
Annual Collector                                                639 / 639
Recession alpha with fixed b = 1                                  60 / 60
Area-normalized gate on Q-to-PPT signatures                        27 / 27
Snow metrics                                                    377 / 377
Stats floor (min_values_for_stats)                                98 / 98
Snow record-anchored decade gate                                  96 / 96
Drought signatures                                             1328 / 1328
Changepoint Detection                                             39 / 39
generate_stats() changepoint integration                          27 / 27
Seasonal runoff-ratio completeness masking                        32 / 32
flagged_for_high_na denominator                                   14 / 14
Climate metrics                                                 369 / 369
```

Python ran with bytecode/cache writes disabled and in-memory output capture. Its two errors were the `tmp_path` setups in `test_phase1_review_fixes.py:46` and 61, not failed numerical assertions.

For the Julia climate test, I evaluated an in-memory copy with redundant `Pkg.activate` and fixture filesystem writes suppressed; **all assertions remained unchanged**. Python read an existing fixture without modifying it. Thus 41/38 are reproduced, but fresh fixture-file generation was not verified.

The installed R package is present, version 0.1.0, built 2026-10-28. Its suite was **blocked**, not skipped because the package was absent.

I cannot certify the documented **43 / 28 / 37** runtime results. Static inspection shows 37 unconditional qualifier checks plus five optional integration checks; “37” does not include integration coverage. That optional branch also contains the sentinel assertion defect above.

**D — Independently checked artifact numbers**

Annual counts came from PyArrow metadata, without reading the annual tables in full. Both staged copies have the same counts and sizes.

| Recomputed quantity | WY 1993–2025 | WY 1980–2025 |
|---|---:|---:|
| Summary rows / unique gages | 6,678 / 6,678 | 6,250 / 6,250 |
| Summary columns | 2,277 | 2,277 |
| Annual rows | **27,766,578** | **34,813,881** |
| Annual parquet bytes | 145,973,271 | 181,665,874 |
| Annual metadata year bounds | 1993–2025 | 1980–2025 |
| Would-qualify list rows | 31 | 29 |
| Class-U gages | 5,069 | 4,894 |
| Restored ledger days inside window | 144,724 | 119,933 |
| Gages holding those restored days | 899 | 757 |
| Listed gained years / affected gages | 1,177 / 619 | 950 / 502 |
| Listed lost years | 0 | 0 |
| `flagged_for_high_na` TRUE | 471 | 345 |
| Raw-unit gages | 32 | 28 |

The summary gage order equals both the frozen selection manifest and the previous-build summary order. The gained-year entries match the corresponding oracle entries. These are independently counted artifact entries; I did not rerun the entire qualification computation.

Additional checks:

- Correction ledger: **193,542 rows, 987 gages, every row marked restored**. Codes: 181,461 `A R`; 11,444 `A <`; 218 `A >`; 419 `A [4]`.
- Both summary dictionaries contain **2,277 unique entries**, exactly matching summary column names.
- All **70 staging-manifest entries** passed: 26 input entries and 22 per product. File hashes and sizes matched; explorer directories matched recorded file counts and aggregate bytes.
- Corrected parquet, ledger and corrected metadata sizes/SHA-256 values matched `R/input/provenance.json`.
- R3 metadata confirms **111,624,189 rows** in each streamflow snapshot and **134,605,430 Daymet rows**. Both candidate metadata tables contain **16,994 rows**.
- Both timing records identify clean commit `54f0483b52fbb1bc743df26f5ee2753a43575be1`, `outer`, and execution on **2026-10-27**. The annotated baseline tag peels to that commit. Commits `161c56e` and `713fdbf` are dated **2026-10-28**.

These support the headline product and provenance claims. They do **not** establish that every historical number, diagnostic statistic, external date or filename mentioned anywhere in the long documents has been independently recomputed. The numerical inconsistency I identified in the updated validation prose is the 103-versus-8/7 exception wording.

**E — Public mirror**

I could not clone or independently diff the live mirror: the session offered no writable scratch location, local network access failed DNS resolution, and browser fetches of GitHub, raw content and the tree API also failed.

Applying `publish_to_hisss.sh` exclusions locally yields an expected **410-file publish set**. Within the scanned public source/document formats, I found **149 `docs/plans/` reference lines across 61 files**. Many are code-formatted references rather than clickable hyperlinks, but they still direct public readers to unavailable evidence. I found no private-repository URL matching `SurfaceWaterProjections` or `streamflowSignatures` in that scanned publish set.

Therefore, **mirror equality with HEAD and the claimed publication from `713fdbf` remain unverified**.

**F — Claims verified and outstanding**

- **Claim 1:** ingestion wiring and acceptance logic verified by source inspection; ledger counts verified. R execution remains blocked, and its optional integration test needs correction.
- **Claim 2:** production-helper defaults and runner wiring verified; libraries unchanged. The legacy R runner is an explicit exception.
- **Claim 3:** wrapper counts, inputs, prefixes, parsing and guards verified; Julia/Python climate assertions reproduced. Batch effective-frame recording is defective; R/batch test counts remain unverified.
- **Claim 4:** headline rows, cohorts, restored-day counts, manifests and input hashes verified. No upload was performed; remote HydroShare state was not independently checked.
- **Claim 5:** **not fully satisfied**. Several current-state passages remain stale, and public references/naming need cleanup.

The remaining release work is documentation correction, effective-frame provenance correction, and completion of the blocked tests and mirror comparison.
tokens used
296,861
| Rank | Finding and evidence | One-line fix |
|---|---|---|
| **MAJOR** | **Batch provenance can record the wrong effective frame.** [run_batched_julia.py:465](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:465) records `outer (default)` when the environment is unset. A historical wrapper subsequently selects `left`, e.g. [run_julia_benchmark_drought_1993_2025_60pct.jl:46](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark_drought_1993_2025_60pct.jl:46). Thus the batch fingerprint can disagree with the child’s timing provenance. | Record requested and effective modes separately; obtain the effective mode from the child’s provenance. |
| **MAJOR** | **Current documentation still says the climate fix is on hold.** [docs/SIGNATURES.md:779](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/SIGNATURES.md:779), [docs/DEVELOPMENT.md:82](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/DEVELOPMENT.md:82), and [claude-skill/streamflow-signatures.md:125](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/claude-skill/streamflow-signatures.md:125) describe the old join as current. [run_batched_julia.py:7](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_batched_julia.py:7) likewise describes an unconditional left join. `docs/SIGNATURES.md` was not changed in the review range. | Mark the defect as historical and document `outer` as the default since 2026-10-28. |
| **MAJOR** | **Current Daymet status contradicts the products.** [CHANGELOG.md:34](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/CHANGELOG.md:34) says the rerun decision is pending and “No product uses it yet”; [docs/DATA_SOURCES.md:18](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/DATA_SOURCES.md:18) says the reprocessed input is “NOT yet used by either product.” Both product provenance records identify that input. | Replace these pending-status statements with the completed October reruns and outstanding deposit status. |
| **MINOR** | **Optional qualifier integration test contains contradictory assertions.** [R/tests/test_qualifier_codes.R:160](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/R/tests/test_qualifier_codes.R:160) inserts `−999999` on an `A R` row; line 170 requires every approved-remark row to be non-NA, while line 171 correctly requires the sentinel to be NA. This is a static finding; R execution was blocked. | Exclude unusable values from the approved-remark “kept” assertion; retain the sentinel assertion. |
| **MINOR** | **The legacy R runner is incorrectly presented as canonical.** [docs/DEVELOPMENT.md:739](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/DEVELOPMENT.md:739) calls `run_r_benchmark.R` canonical. That runner calls the legacy implementation, whose join remains `all.x=TRUE` at [R/helperFunctions.R:4851](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/R/helperFunctions.R:4851); it does not use the climate-frame switch. | Label this runner deprecated and direct current R workflows to `run_rpkg_benchmark.R`. |
| **MINOR** | **Public evidence references point to excluded material.** [README.md:181](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/README.md:181), [docs/CROSS_LANGUAGE_STATUS.md:35](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/CROSS_LANGUAGE_STATUS.md:35), and both staged run notes at line 79 refer readers to `docs/plans/`, excluded by [publish_to_hisss.sh:29](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/publish_to_hisss.sh:29). The staged notes explicitly say the plan is in the HISSS repository. | Publish a concise verification record outside `docs/plans/` and redirect public references. |
| **MINOR** | **Internal version naming reaches public documents.** [README.md:129](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/README.md:129), line 181, and [python/README.md:263](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/README.md:263) use “version 3”; both staged run-note titles expose “internal build v3.” The expected publish set contains 23 such Markdown lines across nine files. | Use public “version 2” in reader-facing prose; retain internal build identifiers only where needed for provenance. |
| **MINOR** | **Reproduction examples are incomplete or obsolete.** Both `D/resource{1,2}_signatures_*/hisss_run_notes_*.md:42` omit required `--out-dir`; executing the shown arguments produced an argparse error. [.claude/skills/run-benchmark/SKILL.md:63](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/.claude/skills/run-benchmark/SKILL.md:63) still recommends standalone product-wrapper launches, which now hit the guard and do not enforce the frozen cohort. | Supply a complete split/run/merge example with output paths and the frozen selection manifest. |
| **MINOR** | **Python README conflates two exception counts.** [python/README.md:266](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/python/README.md:266) says 8/7 gate-state exceptions are “of” the 103 waived cells. Both staged summary-gate logs, line 7, report **103 finite-value waivers plus separately 8/7 gate-state exceptions**, consistent with `docs/CROSS_LANGUAGE_STATUS.md:59`. | Change the wording to “103 finite-value cells and, separately, 8/7 gate-state cells.” |
| **MINOR** | **README contains stale current schema descriptions.** [README.md:72](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/README.md:72) calls the full summary 1,653 columns despite the current 2,277-column schema; line 114 says “eight” categories immediately after listing nine. | Update the current-product description to 2,277 columns and nine categories; retain dated historical comparisons. |
| **NOTE** | **Verification remains incomplete because of environment restrictions.** R could not create its temporary directory; batch self-tests could not create scratch files; mirror access failed. Python’s two setup errors also arose from temporary-directory restrictions. | Repeat the blocked checks in an authorized writable scratch environment with public network access. |

**Verdict: GO WITH CHANGES.** I found no demonstrated corruption of the reviewed products. The default flip itself is implemented, but documentation currency and effective-frame provenance need correction. This is not an unconditional test or public-mirror certification.

No repository, run-folder, staged-deposit or worktree files were edited. HEAD was `713fdbfc655dffab8489d0922fb31340811e20a4`. Contrary to the clean-tree premise, one untracked file existed throughout: `docs/plans/2026-10-21-qualifier-reviews/review_brief_currency.md`.

**A — Stale-text sweep**

I searched source, documentation, hidden rules/skills and staged text for the four-code rule, `keep_codes`, `USGS_ACCEPTED_CODES`, streamflow-date climate joins, left-default wording, `left (unset)`, both old row counts, version-3 wording and deposit claims.

The material distinctions are:

- **Old qualifier rule:** `docs/WORKFLOW_REVIEW.md:41` belongs to an explicitly historical February review. `docs/SIGNATURE_GUIDELINES.md:540` remains stale, but is a dated external-document snapshot; `docs/reconciliation/guidelines_todos.md:19` explicitly tracks its correction. These are not evidence that ingestion still uses the old rule.
- **Removed constant:** references remain in `config.R:122`, the changelog and plan history, but as descriptions of its removal. I found no executable reference.
- **Historical analysis tools:** `min_flow_rule_footprint.py:3` and `min_flow_rule_replay.py:34` describe the four-code rule as belonging to `generate_streamflow_dt`. Their historical replay purpose is legitimate; those comments should explicitly say “February 2026.”
- **Left default:** the actual helper defaults are corrected. No remaining literal `left (unset)` was found in the scanned current source. The substantive stale join descriptions are in the findings table.
- **Old annual counts:** 27,602,161 / 34,681,331 occur in dated October 12/13 build descriptions and old-versus-new comparisons. Those historical numbers should remain, with clear superseded-build labeling.
- **Deposit wording:** current product-status passages say the corrected products are not deposited. However, `D/resource3_input_data/hisss_abstract_resource3.txt:1` says “deposited twice,” and `docs/reconciliation/manuscript_log.md:17` says “the deposited products restore…”. These read as publication-ready wording, not evidence of an upload. “Included as two snapshots” and “staged products” would avoid ambiguity before deposit. References to the earlier HydroShare products are a separate, historical claim.

**B — Code and documentation agreement**

The ingestion call chain is intact:

`run_ingest_usgs_hydat.R:4–6` sources configuration, helpers and qualifier functions; its USGS calls at lines 93 and 105 enter `process_gages_rawToRaw`; `R/helperFunctions.R:386` calls `generate_streamflow_dt`; line 650 applies `usgs_qualifier_mask` before unit conversion.

`R/usgs_qualifiers.R:30–47` accepts status token `A` or `P`, requires a finite numeric value and rejects `−999999`. Remark codes do not mask values. The original flag survives selection and output at `R/helperFunctions.R:640` and 431–434. Gage admission remains a separate gate; accepted gages retain all dates in the requested window.

`splice_nwis_qualifiers.py:357–369` implements the corresponding numeric/status acceptance rule. Its **restoration scope is deliberately narrower**: the four approved remark codes, rather than retrospective restoration of every provisional remark.

The production runners resolve and use the helpers as follows:

| Entry point | Resolution | Application |
|---|---|---|
| Julia | `climate_frame.jl:36`; runner line 51 | Runner line 233 |
| Python | `climate_frame.py:38`; runner line 171 | Runner line 330 |
| rpkg | `climate_frame.R:27`; runner line 75 | Runner lines 288–289 |

All three normalize case/whitespace, default to `outer`, accept `left` explicitly and reject unsupported modes. Outer additions are restricted to the streamflow record span and water years containing streamflow rows. The deprecated legacy branches retain left joins.

I independently counted **12 Julia wrappers: 10 `left`, 2 `outer`**. The two outer wrappers are `run_julia_benchmark_prod_{1993,1980}_2025_60pct_climate.jl`. Every other `run_julia_benchmark_*.jl` wrapper selects `left` through `get!`, preserving an explicit environment override.

Both climate wrappers use the corrected input and `streamflow_<window>_60pct_climate_v3_26oct2026` prefix. Both guards fired during read-only evaluation. Their precise protection is: **unset output directory plus an existing, nonempty default product directory**. They are not an unconditional prohibition on explicitly supplying that directory.

The batch runner’s merged output takes the child’s actual frame from timing provenance (`run_batched_julia.py:725`); its earlier batch fingerprint has the unresolved-default defect described above.

The documented rerun CLI flags exist, including `--gage-list`, `--climate-frame`, `--qa-schema`, `--stage`, `--require-clean-tree`, `--gages`, `--allow-diff-signature`, `--waive-column`, `--allow-gate-state-column`, `--accept-waiver` and `--floor`. The demonstrated CLI problem is a **missing required argument**, not an invented flag.

`git diff 54f0483..HEAD -- julia/src python/streamflow_signatures rpkg` was empty. The libraries do not resolve `STREAMFLOW_CLIMATE_FRAME`; that policy belongs to the runners.

**C — Tests**

| Check | Observed result |
|---|---|
| Julia unit suite | Exit 0; **3,184 passes across 13 testsets**, summed from the emitted counts below |
| Python suite | `175 passed, 9 skipped, 145 warnings, 2 errors in 7.80s` |
| Installed rpkg testthat suite | Blocked before tests: `Fatal error: creating temporary file for '-e' failed` |
| R qualifier suite | Blocked before tests: `Fatal error: cannot create 'R_TempDir'` |
| Julia climate self-test | `41 passed, 0 failed` — read-only adaptation described below |
| Python climate self-test | `38 passed, 0 failed` |
| R climate self-test | Blocked before tests: `Fatal error: cannot create 'R_TempDir'` |
| Batch-runner self-test, daymet venv | Blocked at `tempfile.mkdtemp`: `No usable temporary directory found` |

Julia’s emitted pass/total counts were:

```text
StreamflowSignatures.jl                                           78 / 78
Annual Collector                                                639 / 639
Recession alpha with fixed b = 1                                  60 / 60
Area-normalized gate on Q-to-PPT signatures                        27 / 27
Snow metrics                                                    377 / 377
Stats floor (min_values_for_stats)                                98 / 98
Snow record-anchored decade gate                                  96 / 96
Drought signatures                                             1328 / 1328
Changepoint Detection                                             39 / 39
generate_stats() changepoint integration                          27 / 27
Seasonal runoff-ratio completeness masking                        32 / 32
flagged_for_high_na denominator                                   14 / 14
Climate metrics                                                 369 / 369
```

Python ran with bytecode/cache writes disabled and in-memory output capture. Its two errors were the `tmp_path` setups in `test_phase1_review_fixes.py:46` and 61, not failed numerical assertions.

For the Julia climate test, I evaluated an in-memory copy with redundant `Pkg.activate` and fixture filesystem writes suppressed; **all assertions remained unchanged**. Python read an existing fixture without modifying it. Thus 41/38 are reproduced, but fresh fixture-file generation was not verified.

The installed R package is present, version 0.1.0, built 2026-10-28. Its suite was **blocked**, not skipped because the package was absent.

I cannot certify the documented **43 / 28 / 37** runtime results. Static inspection shows 37 unconditional qualifier checks plus five optional integration checks; “37” does not include integration coverage. That optional branch also contains the sentinel assertion defect above.

**D — Independently checked artifact numbers**

Annual counts came from PyArrow metadata, without reading the annual tables in full. Both staged copies have the same counts and sizes.

| Recomputed quantity | WY 1993–2025 | WY 1980–2025 |
|---|---:|---:|
| Summary rows / unique gages | 6,678 / 6,678 | 6,250 / 6,250 |
| Summary columns | 2,277 | 2,277 |
| Annual rows | **27,766,578** | **34,813,881** |
| Annual parquet bytes | 145,973,271 | 181,665,874 |
| Annual metadata year bounds | 1993–2025 | 1980–2025 |
| Would-qualify list rows | 31 | 29 |
| Class-U gages | 5,069 | 4,894 |
| Restored ledger days inside window | 144,724 | 119,933 |
| Gages holding those restored days | 899 | 757 |
| Listed gained years / affected gages | 1,177 / 619 | 950 / 502 |
| Listed lost years | 0 | 0 |
| `flagged_for_high_na` TRUE | 471 | 345 |
| Raw-unit gages | 32 | 28 |

The summary gage order equals both the frozen selection manifest and the previous-build summary order. The gained-year entries match the corresponding oracle entries. These are independently counted artifact entries; I did not rerun the entire qualification computation.

Additional checks:

- Correction ledger: **193,542 rows, 987 gages, every row marked restored**. Codes: 181,461 `A R`; 11,444 `A <`; 218 `A >`; 419 `A [4]`.
- Both summary dictionaries contain **2,277 unique entries**, exactly matching summary column names.
- All **70 staging-manifest entries** passed: 26 input entries and 22 per product. File hashes and sizes matched; explorer directories matched recorded file counts and aggregate bytes.
- Corrected parquet, ledger and corrected metadata sizes/SHA-256 values matched `R/input/provenance.json`.
- R3 metadata confirms **111,624,189 rows** in each streamflow snapshot and **134,605,430 Daymet rows**. Both candidate metadata tables contain **16,994 rows**.
- Both timing records identify clean commit `54f0483b52fbb1bc743df26f5ee2753a43575be1`, `outer`, and execution on **2026-10-27**. The annotated baseline tag peels to that commit. Commits `161c56e` and `713fdbf` are dated **2026-10-28**.

These support the headline product and provenance claims. They do **not** establish that every historical number, diagnostic statistic, external date or filename mentioned anywhere in the long documents has been independently recomputed. The numerical inconsistency I identified in the updated validation prose is the 103-versus-8/7 exception wording.

**E — Public mirror**

I could not clone or independently diff the live mirror: the session offered no writable scratch location, local network access failed DNS resolution, and browser fetches of GitHub, raw content and the tree API also failed.

Applying `publish_to_hisss.sh` exclusions locally yields an expected **410-file publish set**. Within the scanned public source/document formats, I found **149 `docs/plans/` reference lines across 61 files**. Many are code-formatted references rather than clickable hyperlinks, but they still direct public readers to unavailable evidence. I found no private-repository URL matching `SurfaceWaterProjections` or `streamflowSignatures` in that scanned publish set.

Therefore, **mirror equality with HEAD and the claimed publication from `713fdbf` remain unverified**.

**F — Claims verified and outstanding**

- **Claim 1:** ingestion wiring and acceptance logic verified by source inspection; ledger counts verified. R execution remains blocked, and its optional integration test needs correction.
- **Claim 2:** production-helper defaults and runner wiring verified; libraries unchanged. The legacy R runner is an explicit exception.
- **Claim 3:** wrapper counts, inputs, prefixes, parsing and guards verified; Julia/Python climate assertions reproduced. Batch effective-frame recording is defective; R/batch test counts remain unverified.
- **Claim 4:** headline rows, cohorts, restored-day counts, manifests and input hashes verified. No upload was performed; remote HydroShare state was not independently checked.
- **Claim 5:** **not fully satisfied**. Several current-state passages remain stale, and public references/naming need cleanup.

The remaining release work is documentation correction, effective-frame provenance correction, and completion of the blocked tests and mirror comparison.

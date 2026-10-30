# Sonnet review — code and documentation currency after the default flip, area R (2026-10-29, blind, read-only)

Brief: `review_brief_currency.md` (same folder). Verdict GO WITH CHANGES; every finding applied 2026-10-29 (dispositions in `../2026-10-21-qualifier-mask-fix-rerun-plan.md` §11, row "currency review"). Verbatim report:

# Currency review, area R (ingestion path + rpkg), HEAD 713fdbf, read-only

## Verdict: GO WITH CHANGES
No BLOCKER. The removed `USGS_ACCEPTED_CODES` is read nowhere; `source("config.R")` succeeds. The ingestion rule is applied at exactly one point and matches the docs. rpkg is byte-unchanged since 54f0483. All R-side tests and self-tests pass. The changes are doc wording: one undocumented code-versus-data divergence (P-remark days), a public README that leaks "version 3" and internal paths, and a few unmarked stale present-tense statements.

## Findings (ranked)
| # | Sev | file:line | Evidence | One-line fix |
|---|---|---|---|---|
| 1 | MAJOR | docs/DEVELOPMENT.md:782-794, docs/DATA_SOURCES.md:15 vs deposit resource3_input_data/hisss_readme_inputs.md:83 | The docs say the rule is "remark codes never mask" and the patched parquet is "built from it". The new R rule (R/usgs_qualifiers.R:41-48) keeps a `P <remark>` day that has a finite value. The patched parquet and the products leave the 15,257 provisional-remark days (263 gages) NA (splice only `SPLICE_CODES = A R, A <, A >, A [4]`, splice_nwis_qualifiers.py:63; deposit README says they "stay NA"). Re-ingesting with the fixed R code would therefore NOT reproduce the patched parquet. Only CHANGELOG.md:150 and :443 and the deposit README say so. | Add one sentence to DEVELOPMENT.md (USGS block) and DATA_SOURCES.md row 1: "the 15,257 `P`+remark days stay NA in the patched parquet and products; a fresh ingestion with `usgs_qualifiers.R` would keep those with a value". |
| 2 | MAJOR | README.md:179, :181, :129 (public README, in the mirror) | Line 179 describes "Version 2 (built 2026-10-12/13 ... staged for deposit, not yet published)" with private `~/HISSS_data/...` paths. Line 181 then describes "Version 3 ... to be published as 'version 2'". The public text has two "Version 2" entries (the second is named "Version 3" in its heading) and an internal-name leak. Both entries link `docs/plans/...`, which is excluded from the mirror, so those links are dead for public readers. | Fold 179 into 181 (version-2 folders were never deposited) and drop the "Version 3" and internal-path wording from the public README, or state the naming once. |
| 3 | MINOR | docs/benchmarks/min_flow_rule_footprint.py:3 | Docstring says "generate_streamflow_dt keeps only A, A e, P, P e" in the present tense. The tool is a dated 2026-10-21 forensic of the Feb pipeline, and the same file is header-dated. | Prefix "(rule until 2026-10-26)". |
| 4 | MINOR | docs/benchmarks/min_flow_rule_replay.py:34 | `KEEP_CODES = {"A","A e","P","P e"}  # USGS mask in generate_streamflow_dt`. The comment reads as current. It is the Feb-rule replay and is correct for that purpose. | Change the comment to "the Feb 2026 mask (replay)". |
| 5 | MINOR | docs/WORKFLOW_REVIEW.md:41-43 | "USGS Accepted Codes: A, A e, P, P e; All other codes -> masked as NA". The file header says it is a Feb 2026 document written when R was canonical, so it is historical, but this section carries no date of its own. | Add "(rule until 2026-10-26; see R/usgs_qualifiers.R)". |
| 6 | MINOR | docs/SIGNATURE_GUIDELINES.md:540 | Snapshot of the co-authors' Google Doc still says "kept only for qualifiers A, A e, P, P e (others set to NA)". Ground truth owned by the user. Relay item is already open and marked DUE in docs/reconciliation/guidelines_todos.md:20. | Relay to the user (cannot edit from here); no repo change. |
| 7 | MINOR | R/tests/test_qualifier_codes.R:147-175 | The 5 integration checks (section 7, `generate_streamflow_dt`) are SKIPPED here because `lubridate` is not installed. The run reports 37 passed, which excludes them (42 where lubridate exists). Fixture text for the documented count is unclear. | State "37 (+5 integration checks when lubridate is installed)" wherever 37 is quoted; DEVELOPMENT.md:794 quotes no count. |
| 8 | NOTE | R/tests/test_qualifier_codes.R | Untested rule branches: tab/newline whitespace in `trimws`; a factor `flag` (as.character path); a glued code with no separator ("Ae", "A<"), which is masked because the status token is "Ae"; non-numeric `value` through `usgs_qualifier_mask`, which `keep` covers but `mask` does not directly. Every named case in the brief is covered (see table below). | Optional: add 3 checks. |
| 9 | NOTE | docs/benchmarks/run_r_benchmark.R | It does NOT use `climate_frame.R` and has no `STREAMFLOW_CLIMATE_FRAME` handling. It is the legacy monolithic runner (`process_signatures_from_parquet`, R/helperFunctions.R:4467, with its own climate join and hard-coded D:/ paths). Brief claim 2 ("run_r_benchmark.R / run_rpkg_benchmark.R join through climate_frame") is true for `run_rpkg_benchmark.R` only. DEVELOPMENT.md:739 still describes it as "R canonical full signature extraction" (canonical is Julia). | Say in DEVELOPMENT.md that the legacy R runner has no frame switch (always its own join) and is not canonical. |
| 10 | NOTE | docs/benchmarks/run_rpkg_acceptance_gates.sh:30-37 | Defaults are the August port-campaign reference (`DIR=/Volumes/Untitled/...portref_24aug2026`, `PREFIX=rpkg_full_26aug2026`, `EXPECT_COLS=1653`). It is parameterised by env (header says so), so it is not stale. It has no frame handling; frame provenance is in the runner's timing.json. A caller must set EXPECT_COLS=2277 for the version-3 products. | None needed; the header says it. |
| 11 | NOTE | rpkg/README.md:141-143 | "Monolithic R's process_signatures_from_parquet() removes all NA-Q rows before merging climate data" describes the legacy path, not the runners' join. Current for the legacy path. Nothing in rpkg/README.md, rpkg/R roxygen or rpkg/man mentions the qualifier rule or the climate join. | None. |

## Call-chain trace (what the ingestion keeps)
1. `run_ingest_usgs_hydat.R:4-6`: `source("config.R")`, `source("R/helperFunctions.R")`, `source("R/usgs_qualifiers.R")`. R/helperFunctions.R:34-39 also auto-sources `R/usgs_qualifiers.R` when it is found from the cwd (the repo-root convention of every `run_*.R`).
2. `run_ingest_usgs_hydat.R:79,93,105` call `process_gages_rawToRaw(...)`. Three runs: Canada, AK and CONUS (the CONUS and AK calls use gage_type "USGS").
3. `process_gages_rawToRaw` (R/helperFunctions.R:303) calls `generate_streamflow_dt(current_gage, gage_type, ...)` at :386. The other callers are `process_gages_rawData` (:151) and `process_gages` (:4367), all the same function.
4. `generate_streamflow_dt` (R/helperFunctions.R:596) for USGS: `readNWISdv` (:612), column pick `Flow` / `Flow_cd` or `X_00060_00003(_cd)` (:622-639), rename to `Date, Q_rawUnits, flag` (:640).
5. R/helperFunctions.R:647-650: stop if `usgs_qualifier_mask` is missing, then `streamy$Q_rawUnits <- usgs_qualifier_mask(streamy$flag, streamy$Q_rawUnits)`. This is the ONLY keep/mask point for a USGS value.
6. Rule (R/usgs_qualifiers.R:39-49, 52-56): keep iff `usgs_qualifier_status(flag)` is in `c("A","P")` (first token after trimws, split on whitespace, ':' or ','; case-sensitive; NA becomes "") AND the value is finite AND it is not -999999. A value that is not kept becomes `NA_real_`. `flag` is stored verbatim (the `flag` column is never edited).
7. Then cfs to mm/day conversion (:652-662). Canada (:676+) is never masked; `Symbol` goes to `flag`.
8. Year rejection stays in `preprocess_daily_data()` and is not touched here. The per-calendar-year compile floor `MIN_Q_VALUE_AND_DAYS` is a count, not an edit (config.R:36-44).

`config.R` carries only a comment about the removed constant (config.R:118-123). `Rscript -e 'source("config.R")'` printed `config ok`; `exists("USGS_ACCEPTED_CODES")` is FALSE. No `.R`, `.sh`, `.py` or `.jl` file reads the name (grep of the repo excluding worktree, plans and .venv).

## Grep classification (USGS_ACCEPTED_CODES, keep_codes, "A e", P e, accepted codes, qualif)
| Hit | Class |
|---|---|
| config.R:44, 118-123 (rule pointer; names the removed constant in a comment) | current-and-correct |
| R/usgs_qualifiers.R:1-23, 28; R/helperFunctions.R:30-39, 642-650; run_ingest_usgs_hydat.R:6 | current-and-correct (the old four strings appear only as dated history) |
| R/tests/test_qualifier_codes.R:33, 50, 90, 125, 159 | current (test fixtures; :125 is the explicit old-rule comparator) |
| docs/DATA_SOURCES.md:15; docs/DEVELOPMENT.md:782-794 | current-and-correct (new rule plus dated history), gap in finding 1 |
| docs/STATUS.md:28 ("the R ingestion's `keep_codes` drops ...") | historical-and-dated (found 2026-10-21, same entry says FIXED 2026-10-26) |
| docs/benchmarks/splice_nwis_qualifiers.py:7, 39, 62-63, 83, 592 | current-and-correct (`KEPT_CODES_FEB` is the dated Feb rule, SPLICE_CODES the targets) |
| docs/benchmarks/min_flow_rule_replay.py:34; min_flow_rule_footprint.py:3 | historical, undated in text (MINOR 3, 4) |
| docs/WORKFLOW_REVIEW.md:41 | historical (file header says Feb 2026), MINOR 5 |
| docs/SIGNATURE_GUIDELINES.md:540 | STALE but user-owned snapshot; relay item open (guidelines_todos.md:20, 255-256) |
| docs/reconciliation/guidelines_todos.md:20, 255; docs/filtering_min_flow_rule_2026-10-21.csv (rows "A\|A e", "A R days set to NA") | historical-and-dated |
| CHANGELOG.md / changelog-old.md / CHANGELOG_ARCHIVE.md / docs/plans/* | not classified line by line, history by nature |
| claude-skill/streamflow-signatures.md:217, 271-279 | current (the text says the rule was fixed 2026-10-26) |
| README.md:129, 179, 181 | current except the public wording, MAJOR 2 |
| rpkg/** (README, R, man, tests) | no statement about the qualifier rule; `rpkg/R/io.R:286` is the Ice note and is correct |
| EO_data_processing/, .claude/rules, .claude/skills | no stale qualifier-rule hit (`benchmarks.md:41-49` and `run-benchmark/SKILL.md:47` describe the rerun correctly) |
| Deposit hisss_readme_inputs.md:74-83, hisss_input_data_dictionary.csv:10 | current. `flag` is "the February 2026 retrieval's code" kept verbatim, including restored days (old code stays). Sentence on which days are kept: status A or P plus finite value, 193,446 approved and 96 provisional restored; P-remark days stay NA |

Other stale-text sweeps: no "left (unset)" and no "left as default" statement in the R-side files. The `left` default text appears only in the Julia wrapper comments, which are correct ("pinned so this wrapper keeps reproducing"). The 27,602,161 and 34,681,331 counts appear as version-2 rows (STATUS.md:12-13, README.md:179, DEVELOPMENT.md:606), labelled version 2 and historical. No claim of a HydroShare deposit is made as done.

## Rule branch to test coverage (R/usgs_qualifiers.R)
| Case | Covered (test line) |
|---|---|
| A R, A <, A >, A [4], A e, P, P e, A R e, P Ice etc. | yes (:32-47) |
| Bare remark "e" | yes (:79-80) |
| -999999 sentinel (numeric and string), custom sentinel | yes (:56-61, 88) |
| NA value, NaN, Inf, "Ice" as value | yes (:54, 62-66) |
| NA / empty flag | yes (:75-78) |
| Lower-case, leading / trailing spaces, ':' and ',' separators | yes (:81-91) |
| Tally and the old-versus-new comparison | yes (sections 5, 6) |
| Integration through `generate_streamflow_dt` | present but SKIPPED here (no lubridate) |
| Tab / newline, factor input, glued "Ae" | no (NOTE 8) |

## Test outputs (verbatim)
- `Rscript R/tests/test_qualifier_codes.R`: last lines `  SKIP: lubridate not installed (run on the ingestion machine)` then `37 passed, 0 failed`. The documented 37 is the real count here (5 integration checks not run).
- `selftest_climate_frame.jl` (fixture written to scratch .../review_currency/r/fixture): `41 passed, 0 failed`.
- `Rscript docs/benchmarks/selftest_climate_frame.R <that fixture>`: `43 passed, 0 failed` (exit 0, 43 PASS lines, 0 FAIL). Includes the default check `climate_frame_mode() is "outer" when STREAMFLOW_CLIMATE_FRAME is unset`. The documented 41 / 43 are real.
- rpkg testthat, installed package 0.1.0 (run with `testthat::test_dir("rpkg/tests/testthat", package="streamflowsignatures")`, the repo `testthat.R` convention): `TOTAL tests 98 expectations 1218 failed 0 skipped 0 errors 0 warnings 6`. The warnings are "collapsing to unique 'x' values" in test-smoke.R:46. Package already installed; I did not install anything.
- `git diff 54f0483..HEAD --stat -- rpkg/` returned empty (no change). `rpkg/inst/config/signatures_config.json` is byte-identical to `config/signatures_config.json` (cmp).
- Not run (outside my area): 38 (Python), 28, Julia and pytest suites.
- `git status` after my runs: only the pre-existing untracked brief file. Nothing was written to the repo.

## Climate frame in the R runners
- `climate_frame.R:28-35`: `STREAMFLOW_CLIMATE_FRAME` unset gives "outer" (default since 2026-10-28); lower-cased and trimmed; any value other than left/outer errors. Matches the header comments (:20-22).
- `run_rpkg_benchmark.R:75`: `climate_frame <- climate_frame_mode()`, printed at :91. At :285-289 `frame_with_climate(gage_data, climate_merge, if (isTRUE(use_legacy)) "left" else climate_frame)`; the deprecated legacy path is pinned to left, commented. SWE: `clim_cols <- intersect(c("date","PPT","SWE","TMIN","TMAX"), names(gage_climate))` (:280) goes into the same merge, so SWE, TMIN and TMAX ride the same frame; `has_climate` means PPT is present.
- Provenance recorded: `timing$provenance$climate_frame` (:657), `STREAMFLOW_CLIMATE_FRAME` in the `env_overrides` list (:645), plus streamflow / climate / metadata / config hashes and info, git revision and dirty flag, package versions.
- `run_r_benchmark.R`: no switch (finding 9).

## Verified vs could not
Verified myself:
- The rule text, the single keep/mask point and the call chain.
- The config.R source and the absence of the removed constant.
- The test run (37), the Julia self-test (41) and the R self-test (43), the rpkg testthat run, and the rpkg diff (empty) and config identity.
- The runner default resolution and provenance.
- The grep classification.
- The deposit text on `flag` and the kept-day sentence.

Could not verify:
- The 5 integration checks of the R unit test (lubridate absent; I did not install it).
- The Python and Julia suites, `selftest_run_batched_julia.py` (28) and the Python self-test (38).
- Whether `readNWISdv` still returns these exact qualifier strings (no network use, by the hard rules).
- Row counts, sizes and manifests of the deposit and run folder (area D).
- The mirror diff (area E). README dead-link status is inferred from the "excluded plans/" note in the brief, not from a mirror clone.

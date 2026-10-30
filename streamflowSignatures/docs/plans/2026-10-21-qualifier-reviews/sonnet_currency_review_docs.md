# Sonnet review — code and documentation currency after the default flip, area DOCS (2026-10-29, blind, read-only)

Brief: `review_brief_currency.md` (same folder). Verdict GO WITH CHANGES; every finding applied 2026-10-29 (dispositions in `../2026-10-21-qualifier-mask-fix-rerun-plan.md` §11, row "currency review"). Verbatim report:

# Review: documentation, numbers and public mirror (blind, read-only; HEAD 713fdbf)

## Verdict: GO WITH CHANGES (no blocker)
Every number I could recompute matches the artefacts, and the mirror equals HEAD minus the exclusions. The defects are stale prose: the climate-frame "fix on hold" Known Issue is still written as current in five places, STATUS/CHANGELOG still carry closed Daymet items, and the public README exposes internal "version 3" naming and dead pointers.

## Findings (ranked)
| Rank | File:line | Evidence | Fix |
|---|---|---|---|
| MAJOR | docs/SIGNATURES.md:332-333, 540-541, 779-788; docs/DEVELOPMENT.md:82-86; claude-skill/streamflow-signatures.md:125-130 | "KNOWN ISSUE (2026-10-16, fix on hold)", "runners join Daymet onto the streamflow ROWS", "documented, not fixed until discussed" presented as current. Code default is `outer` (climate_frame.{jl,py,R}); products built with it. | Rewrite: fixed in the version-3 products, `outer` default since 2026-10-28; earlier products carry it. |
| MAJOR | docs/STATUS.md:23; CHANGELOG [Unreleased] Planned entry "Daymet climate input: product-rerun decision is PENDING (NEXT)" | "No product uses it yet", "NEXT: rerun" contradict the STATUS table (v2/v3 use `daymet_1980_2025_29sep2026.parquet`). Closed item still listed. | Delete or reduce to "used by v2/v3; open: co-author confirmation of polygons/weights/fill". |
| MAJOR | README.md:179, 181 (public front page) | Two "Version 2" labels (undeposited v2 and "Version 3 ... published as version 2"), internal build name, `~/HISSS_data`, `~/Downloads` paths, `docs/plans/...` pointers (excluded from mirror). | Public text: one "version 2 (October 2026, not yet published)" paragraph; move internal naming to CHANGELOG. |
| MINOR | docs/STATUS.md:30, claude-skill:211 | `flagged_for_high_na` 560 / 401 is the v2 value. v3 CSVs: 471 (#1) / 345 (#2) (recomputed; deposit run notes say 471). 471/345 appear in no repo doc. | Record 471 / 345 for the v3 products. |
| MINOR | CHANGELOG Known Issues USGS entry (~l.136-190); STATUS.md:28 | "Fix = the rule in `keep_codes` ... — NOT applied" followed by "Phase A-D DONE"; present tense "keep_codes drops". | Reword to past tense / "applied". |
| MINOR | .claude/skills/run-benchmark/SKILL.md:63-64 | Wrapper comments say "product #1 VERSION 2"; wrappers are now the version-3 ones. | Say version 3 (public "version 2"). |
| MINOR | docs/STATUS.md:22, 12-13, 25 | Climate-metrics item still says "NEXT (user): review staged folders, decide the deposit" (v2 superseded); rows 1 v2 / 2 v2 kept; Earthdata token expiry (2026-10-25) has passed. | Condense v2 rows/item; drop token line once revoked. |
| MINOR | Staged deposit (run notes l.8, l.79; readme l.36) | Public text names undeposited "version-2 build ... climate_12oct2026", cites `docs/plans/...§11` (not on mirror); provenance/validation files carry `/Users/arik/...` paths. | Accept as provenance, or strip the plan pointer. |
| NOTE | docs/DEVELOPMENT.md:769 | "the version-2 prefix and inputs until 2026-10-28" is ambiguous; wrappers now pin the v3 prefix (verified l.49/50). | Clarify. |
| NOTE | Brief claim 5 | docs/plans/2026-10-12-coauthor-relay.md and the HydroShare plan are unchanged in 54f0483..HEAD; the relay has no qualifier/frame item (guidelines_todos.md:19-26 holds it, marked DUE). `run_r_benchmark.R` does not exist (only run_rpkg_benchmark.R). | Add a relay line if wanted. |
| NOTE | docs/benchmarks/min_flow_rule_footprint.py:3, min_flow_rule_replay.py:34, SIGNATURE_GUIDELINES.md:540, WORKFLOW_REVIEW.md:41 | Old four-string rule: census tools, read-only guideline snapshot (relay open), dated Feb-2026 review. | Historical; fine. |
| NOTE | Mirror ships docs/STATUS.md, reconciliation/*, EO_data_processing/CLAUDE.md | They cite excluded `.claude/rules`, `docs/plans`. Pre-existing design. | None. |

## Numbers table
| Value | Documents | Artefact | Match |
|---|---|---|---|
| 6,678 / 6,250 gages | README, STATUS, results, readmes | summary CSV rows (header excluded) 6,678 / 6,250; selection manifests same | yes |
| 2,277 columns | all | CSV header 2,277 both | yes |
| 27,766,578 / 34,813,881 rows | README, STATUS, results, DEV, skill | parquet metadata 27,766,578 / 34,813,881 (run folder and staged copies) | yes |
| v2 27,602,161 / 34,681,331 | STATUS, README, CHANGELOG, DEV | census md: old rows 27,602,161 | yes, historical context only |
| 193,542 days / 987 gages | CHANGELOG, DEV, DATA_SOURCES, skill, results | ledger 193,542 rows, 987 gages; provenance.json restored 193,542 | yes |
| by code A R 181,461 (755 gages), A < 11,444 (139), A > 218 (30), A [4] 419 (101); 12,081 non-A R | STATUS, CHANGELOG | ledger group-by | yes |
| 144,724 / 119,933 in-window days, 899 / 757 gages | results, CHANGELOG, skill, readmes | ledger joined to selection manifests, window dates | yes |
| 908 gages (#1), 757 (#2) full calendar | STATUS, CHANGELOG | 144,885 days / 908; 119,935 / 757 | yes |
| 1,177 / 950 water years, 619 / 502 gages, 0 lost | results, README, skill, readmes | classes CSVs: 1,177 in 619 gages; 950 in 502; lost 0 | yes |
| 5,069 / 4,894 untouched | README, results | classes: U 5,069 / 4,894 | yes |
| 874 / 734 gages with inserted climate days | CHANGELOG, readme | A-J + A-QJ = 874; 599+135 = 734 | yes |
| 31 / 29 would-qualify | many | would_qualify CSVs 31 / 29 data rows; R3 file 60 | yes |
| 1,609 / 1,356 changed gages | results, python README | affected manifests (lines) 1,609 / 1,356 | yes |
| 387,055 / 360,065 non-equal cells | results | census md | yes |
| 117 / 188 values beyond 1e-9 | results, CROSS_LANG, python README | port-check files | yes |
| rpkg v3: 2,227 Perfect / 11 Good / 6 Poor, R2 0.999922, 323 of 26,939,257, 27,766,577 rows | README, CROSS_LANG, results | gate_annual_waived, summary md, parquet | yes |
| flagged_for_high_na 560 / 401 | STATUS, skill, CHANGELOG | v2 CSV 560 / 401; v3 CSV 471 / 345 | v2 yes; v3 not documented |
| 3,241 / 2,221 inserted days | CHANGELOG, skill, log | not recomputed | not verified |
| staged sizes 646 MB / 714 MB / 6.3 GB | brief | du -sh: 646M / 714M / 6.3G | yes |
| staging manifests | | all 22 / 22 / 26 entries: size match, none missing; sha256 checked for the two summary CSVs (match staged and run copies) | yes |
| 111.6 M observations | abstract | patched parquet 111,624,189 rows | yes |
| tag baseline = 54f0483 | results, run notes | `git rev-list` 54f0483b5...; tag on origin | yes |
| mirror "@ 713fdbf (2026-10-28)" | brief | mirror commit d9432d3 message | yes |
| self-test counts 41/38/43/28/37 | CHANGELOG:535-536, STATUS | not re-run (scope C) | not verified |

## Stale-text table
| Phrase | Hit | Class |
|---|---|---|
| `USGS_ACCEPTED_CODES` | config.R:122 (removal note only) | fine (no reader) |
| `keep_codes` | STATUS.md:28, CHANGELOG ~l.160-170 | present-tense wording of a fixed defect: MINOR |
| `A, A e, P, P e` | R/usgs_qualifiers.R:15, R/helperFunctions.R:645, DEV:786, DATA_SOURCES:15, splice tool, tests, census scripts | historical and dated: fine |
| same | SIGNATURE_GUIDELINES.md:540 | read-only snapshot; relay open (guidelines_todos:19) |
| same | WORKFLOW_REVIEW.md:41 | Feb 2026 historical note: fine |
| "left (unset)", default `left` | none found | clean |
| "joins climate onto streamflow dates/rows" as current | SIGNATURES:332,540,779; DEV:82; skill:125 | current-state, must fix (MAJOR) |
| 27,602,161 / 34,681,331 | STATUS rows 1v2/2v2, README:179, DEV:606, CHANGELOG:974,982 | labelled version-2: fine |
| "version 3"/"v3" in public text | README:181, DEV (3), DATA_SOURCES:15, CROSS_LANG (2), python/README:263, CHANGELOG (6) | README must be fixed; others tolerable |
| claims of deposit | none says deposited; all say not yet | clean |

## STATUS.md discipline
49 lines (< 60). One pointer per item. Closed-but-listed: the Daymet reprocess item (l.23), climate-metrics NEXT (l.22), Earthdata token (l.24, date past). Missing: v3 `flagged_for_high_na` 471 / 345. The default-flip decision (2026-10-28) is recorded.

## Mirror diff
Clone of CZ-Sync/HISSS (depth 1): HEAD d9432d3, "Sync from source repo @ 713fdbf (2026-10-28)" (no "+uncommitted" suffix; the source tree was clean). `git ls-files` minus the script's exclusion regexes gives 410 files; the mirror has the same 410; `cmp` found 0 content differences. Source commit is recorded in the mirror commit message.

## Dead-link / local-path list (mirror-shipped docs)
No link to the private repo or code-sandbox except historical mentions (changelog-old.md:280,336; guidelines_todos.md:349). python/README.md and rpkg/README.md have none.
| File:line | Target | Presentation |
|---|---|---|
| README.md:179 | docs/plans/2026-10-12-climate-reviews/, ~/HISSS_data, ~/Downloads | public pointer: fix |
| README.md:181 | docs/plans/2026-10-21-qualifier-rerun-results.md, ~/HISSS_data | public pointer: fix |
| docs/CROSS_LANGUAGE_STATUS.md:12,35,99 / 73,112 | docs/plans/... / ~/HISSS_data | internal record pointers |
| docs/DATA_SOURCES.md:18 / 15 | docs/plans/... / ~/HISSS_data | internal notes |
| docs/DEVELOPMENT.md:86,134,771,792 / 790,888,890 | docs/plans/... / ~/HISSS_data | internal notes |
| docs/SIGNATURES.md:788,802 | docs/plans/... | internal notes |
| CHANGELOG.md (22 plans, 5 local, 4 Downloads); changelog-old.md (27 / 11); reconciliation/guidelines_todos.md (9), manuscript_log.md (16); STATUS.md (5 / 3 / 2) | docs/plans, local paths, `.claude/rules` | internal logs, pre-existing |
| EO_data_processing READMEs / CLAUDE.md | docs/plans, local paths | pre-existing |
Code comments in ~60 files also cite docs/plans (pre-existing; harmless).

## Verified vs not
Verified: all rows marked "yes" above; code default `outer` in the three helpers and run_batched_julia.py help text; 12 Julia wrappers = 10 `left` + 2 `outer`; product-folder guard present in both climate wrappers; ingestion chain (run_ingest_usgs_hydat.R sources R/usgs_qualifiers.R; helperFunctions.R applies `usgs_qualifier_mask`); no other reader of the removed constant; staged dictionaries/categories identical across the three copies.
Not verified: unit suites and self-test counts (scope C, not run); 3,241 / 2,221 inserted-day counts; sha256 of staged files beyond the two CSVs; HydroShare live state; plan §11 line by line; documented CLI flags beyond `--climate-frame`.

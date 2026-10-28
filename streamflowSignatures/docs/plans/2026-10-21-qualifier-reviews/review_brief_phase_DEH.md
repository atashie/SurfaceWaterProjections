# Review brief — qualifier-mask + climate-frame rerun, Phases D / E / H (2026-10-28), blind, read-only

Repository: /Users/arik/github/SurfaceWaterProjections/streamflowSignatures (main tree, uncommitted edits since the tag
`qualifier-rerun-baseline-2026-10-26` = commit 54f0483). Run folder: /Users/arik/HISSS_data/qualifier-rerun-26oct2026 (R).
Staged deposit folders: /Users/arik/Downloads/Signatures_v3_28oct2026 (local only). Do NOT edit anything; do not run a
production job; reading, small Python/duckdb checks and recomputations are welcome (python: $R/../venvs/daymet/bin/python has
duckdb + pandas + pyarrow; the repo .venv has pandas/pyarrow).

## What was done (claims under review)
1. Phase D executed from a CLEAN worktree of the tag (R/worktree/streamflowSignatures, never edited): S0 (Feb input, left)
   bit-identical to the version-2 products on both windows; version-3 products R/processedOuts_{1993,1980}_2025_climate_v3_26oct2026
   (patched input, outer frame, pinned v2 schema, frozen selection manifests R/manifests/selection_<W>.txt);
   S1 / S1′ stages (R/S1_<WU>, R/S1prime_<WU>); B′ replays + daily delta ledgers (R/bprime_<WU>); E7 Python port on the
   affected gages (R/e7_<WU>); censuses (R/census_<WU>); attribution (R/attribution_<WU>). Driver scripts: the scratchpad
   phase_d1.sh, phase_d1b.sh, phase_d2.sh, phase_d2b.sh, phase_d2c.sh, phase_e6.sh, phase_f.sh (same folder as this brief);
   logs R/phase_d1.log, R/phase_d2.log, R/phase_e6.log.
2. Verification tools patched POST-TAG in the main tree (git diff vs the tag): docs/benchmarks/assess_rerun_changes.py
   (Qann fallback for non-replayed gages, record-dependent key emission on old years, scoped total_possible freeze, waiver key
   allow_gate_state_column), docs/benchmarks/preprocess_state_ledger.py (check 2: a climate NA day in a year the VARIABLE is
   invalid in stays NaN), docs/benchmarks/check_summary_equality.py (--allow-gate-state-column: gate-state cells of tie pairs
   accepted only when both sides are reproduced by recomputation), docs/benchmarks/check_annual_parquet_equality.py (--gages),
   docs/benchmarks/stage_attribution.py (new), docs/benchmarks/stage_hydroshare_qualifier_v3.py (new),
   docs/benchmarks/reconstruct_qualifier_patched_streamflow.py (new).
3. Results claimed (plan §11 rows dated 2026-10-27/28 in docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md; CHANGELOG
   [October 2026] Phase D paragraph; docs/STATUS.md): census verdict PASS both windows (0 unexplained of 387,055 / 360,065
   non-equal annual cells; 5,069 / 4,894 untouched gages identical), E7 PASS with named allowances (117 / 188 annual values in
   19 signatures; summary 0 gated, 8 / 7 gate-state cells verified by recomputation), stage identities PASS, E1–E4, E6 PASS.
4. Docs updated for version 3: README.md, docs/STATUS.md, CHANGELOG.md (Known Issues bullets + October entry),
   docs/DEVELOPMENT.md, docs/DATA_SOURCES.md, docs/CROSS_LANGUAGE_STATUS.md, claude-skill/streamflow-signatures.md,
   .claude/rules/benchmarks.md, python/README.md, docs/reconciliation/{guidelines_todos,manuscript_log}.md.
5. Staging (local): three resource folders + check manifests; README / run notes / dictionary text generated from the files.

## Questions (answer each with evidence; name file:line)
A. Is any claim in plan §11 / CHANGELOG / STATUS / README NOT supported by the artefacts in R? Recompute at least two
   headline numbers independently (e.g. the 5,069 untouched gages bit-identical at both layers; the 1,177 / 950 restored
   water years vs R/input/oracle/acceptance_oracle.csv; the 117 / 188 E7 residuals' classes).
B. Do the post-tag tool patches weaken a gate? In particular: (i) does the record-dependent "key emitted / withdrawn on an
   old year" explanation in assess_rerun_changes.py hide a real defect class? (ii) does --allow-gate-state-column accept a
   cell it should not (what exactly is recomputed, and from which side's series)? (iii) does the Qann fallback for
   non-replayed gages make 5f weaker than the oracle check it replaces? (iv) the scoped total_possible freeze.
C. Stage attribution: are the identities meaningful on the executed cohorts (control 200, join-only 197 / 199 — 3 / 1 of the
   sampled join-only gages fall inside the interaction cohort and are excluded)? Is the reported interaction residual
   (elasticity family) plausibly non-additive by construction, or a sign of something else?
D. Staging: open the staged README / run notes / dictionary / abstract and check every number and every sentence against
   the files beside them and the run folder; is anything misleading for a data user (e.g. the public name "version 2", the
   flag column semantics, the two streamflow tables, the would-qualify lists)? Does the reconstruction recipe's contract hold?
E. Anything that must be fixed BEFORE a commit / deposit, vs. can wait. Verdict: GO / GO WITH CHANGES / NO-GO, with findings
   ranked BLOCKER / MAJOR / MINOR / NOTE.

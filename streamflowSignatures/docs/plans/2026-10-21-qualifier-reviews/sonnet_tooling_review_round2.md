# Sonnet 5.5 tooling review, round 2 (2026-10-26, blind, read-only) — Phase B / C / B′ / E8 / E7 tooling after the Codex round-1 fixes

Brief: `scratchpad/codex/review_brief_B_C_tools_v2.md` of the session (the tools, the plan sections, the 40-gage and 5-gage artefacts, six questions, the round-2 context). Verdict GO WITH CHANGES; every finding 1–7, 9 and 12 applied the same day, 8 partly (`normpath`), 10–11 noted (plan §11). Verbatim report:

## What I ran
- **Oracle replay.** I ran `qualifier_year_oracle.py replay` on the February parquet against the Julia census from 2026-10-26. It reproduces all 322,563 year statuses and all 16,028 window rows with 0 mismatches.
- **Ledger re-derivation.** I recomputed every one of the 193,542 ledger `Q_new` values from the NWIS cache with my own code (unqualified block, `value × 86400/(area × 3280.84³) × 1e6`). There are 0 mismatches, and the values are exact.
- **Multi-method sites.** For 01467087, 021473428 and 06294000, all masked days lie in block 0 (the unqualified one). For 06294000 that block equals the February values (0 differences), while the "Adjusted for Agency Canal" block differs on 1,828 days. Choosing the unqualified block is correct.
- **Metadata rule.** The compile rule (`compile_year_stats`) applied to the unpatched February parquet reproduces the February `num_years`, `start_year`, `end_year` and `num_days` for all 987 splice gages. The rule is faithful. The tool never asserts this itself (finding 8).
- **Control run.** I ran `control` on all 987 gages into scratch: 13,786,734 February values compared, 99.770 % equal at 1e-9. The 97 splice-code days that `n_restored_by_rule` omits are exactly the days whose qualifier today is plain `P` (96) or `A` (1).
- **Batch runner and census.** The 28/28 runner self-test passes, and the pinned-header subset run is gate-equal to product #1 v2. I re-ran the census on the 40-gage subset with a fake oracle built from the state tables. The 5d/5f code path executes and passes. Both had been NOT RUN in the earlier census.
- **Chain state.** The chain was still writing the patched parquet when I looked: 193,542 target days, all 193,542 restored, 0 unavailable.

## Findings

1. **MAJOR: the B8 gate takes about 2 hours** (`splice_nwis_qualifiers.py:868-870`). The line `if g not in set(led["gage_id"]) and a.universe_scope:` builds a 193,542-element set on every loop iteration. I timed one build at 38 ms, times 193,542 base targets. The audit will look hung and may be killed.
   - Fix: hoist `led_gages = set(led["gage_id"])` above the loop, and test `a.universe_scope` first.
   - B5's per-row Python loop over about 111 M rows (roughly 10 min) is tolerable.

2. **MAJOR: gate 5c/5g trusts a bare verdict string** (`assess_rerun_changes.py:533`). It passes when `j.get("verdict") == "PASS"` and checks nothing else.
   - It does not check that the JSON's `ref`/`port` are the annual parquet under assessment (`--new-annual`) and its Python twin.
   - It does not check hashes, or that the compared gage set matches. The subset `e7_annual.json` compares `out_outer` against `out_outer_py`.
   - It prints `allow_value_diff=1` and `waived=1` but never fails on them. A hand-written `{"verdict":"PASS"}`, or the E7 JSON of another run, passes 5c.
   - This matters because 5c is the only real value validator for most climate bases. The 5a "0 unexplained" check is a closure test. The reviewer's +123 injection passes 5a (1,811 cells, 0 unexplained, in `inject/census`) and is caught only because E7 fails.
   - Fix: have `check_annual_parquet_equality.py --out` and `check_summary_equality.py --out` record sha256 of both inputs plus the gage count, window and waiver list. Make `assess` require `sha256(ref) == sha256(--new-annual)` (and the CSV for 5g), and fail on non-empty waivers unless an explicit `--accept-waivers N` is given.

3. **MAJOR: the B2 conversion gate is aggregate only** (`splice_nwis_qualifiers.py:723`). It passes at ≥ 95 % overall equality and gates no single gage.
   - In my control run, 03358000 has 33 % of its values differing (max relative difference 31.2). It has 4,265 February dates that are absent today.
   - 02310700 has a maximum relative difference of 11.8.
   - A gage with a wrong `basin_area` would be 100 % off and still pass the aggregate, because B4 re-derives with the same `conversion_for()` and is circular.
   - I verified that the restored days on 03358000 and 02310700 have neighbours (±2 days) equal to February (ratio 1.0).
   - Fix: add a per-gage gate. Require each restored gage's compared values equal at 1e-9, or at least its ±3-day neighbours of every restored day equal. Report and justify the exceptions.

4. **MAJOR: state-table attribution is coarse** (`preprocess_state_ledger.py:305-319`).
   - Any changed field in a water year that contains a splice or inserted daily delta is attributed to "splice day in the year" or "inserted climate row in the year". It does not matter whether the field can be caused by that day.
   - A corrupted `win_complete`, `rejection` or `q_valid` in an affected year would be accepted.
   - Also, a daily "eligibility" cause (`:261`, with a dead `if False` branch) only needs some validity flip in that year, not a root cause.
   - Fix: recompute the expected state change from the daily deltas and compare it. The raw-NA count deltas must equal the number of changed raw days. A seasonal flag may change only in the season of a restored day. `q_valid` and `rejection` must match the oracle or a recomputed rule.

5. **MINOR: B1 does not pin the request window** (`splice_nwis_qualifiers.py:707-714`). It requires one common window across gages but not `1979-10-01..2025-10-01`. The superseded cache `nwis_cache_to_2025-09-30/` would pass B1.
   - The real cache in `input/nwis_cache` is correct: all 1,185 URLs have `endDT=2025-10-01`.
   - Fix: add `--expect-window`, defaulting to the compile window, and gate on it.

6. **MINOR: 5b can pass with nothing compared** (`assess_rerun_changes.py`, 5b block). `n_cmp == 0` is a PASS even when the delta ledger holds PPT rows. The "no changed temperature date" GDD branch also passes by default.
   - Fix: require `n_cmp > 0` whenever PPT or TMIN/TMAX deltas exist.

7. **MINOR: `--qa-schema` silently masks a missing metadata layer** (`run_julia_benchmark.jl:617`, `run_python_benchmark.py:642`, and the R twin).
   - `added` covers every pinned column that is absent, including metadata columns such as `NDAMS_2009`, `CLASS` and `RHBN`. If GAGES-II or HYDAT fails to load, the header still matches and the columns are all missing.
   - Fix: restrict `added` to signature columns, or error when a column in the `METADATA` list is absent. The block is a verified no-op when every column is present.

8. **MINOR: Julia config baked at precompile time is not detected** (`run_batched_julia.py:232-262`).
   - The fingerprint hashes the config file, but `julia/src/config.jl` says the constants are fixed when the package is precompiled. A stale compiled cache with an older config would run under a fingerprint that claims the current config.
   - Fix: have the Julia wrapper write the config sha256 it was compiled with (or `CFG_*` values) into the timing JSON, and compare that in `merge`.
   - The HYDAT path is hashed from env `STREAMFLOW_HYDAT_PATH` or a default, but Julia reads `CFG_HYDAT_PATH` from the config.
   - `repo = os.path.dirname(a.project)` breaks with a trailing slash; use `normpath`.

9. **MINOR: the audit doesn't assert a few things it could** (`splice_nwis_qualifiers.py`).
   - B6 recomputes the new metadata but never asserts that the rule applied to the base series equals the February metadata. I checked this for all 987 gages: 0 differences. Add it as a gate.
   - B3 passes on count equality with the list stage alone, not on restoration. B8 covers that.
   - The `--expect-splice-days` and `--expect-provisional-days` defaults are hard-coded for this universe.
   - The oracle's `masked` table is loaded but unused.
   - B7 (the oracle's final sets against the replay of the patched parquet and the Julia census on it) is not a single command. Add a one-line equality check.

10. **NOTE: provenance of 97 restored days.** 96 days flagged `A R` and 1 flagged `A <` in February now carry status `P` or `A`. They are restored under the A/P rule with the stale February flag kept in `flag`. The ledger records `qualifiers_today`. Mention this in the data notes. Also, 554 restored values are negative, and the config does not reject negative flow.

11. **NOTE: census closure and coverage** (`assess_rerun_changes.py`).
    - Explanations are year-level, not season-level. For example, a winter-only base changing in a year whose restored day is in summer is explained.
    - The registry's `record_dependent` label is hand-set and unproven. An over-labelled base hides changes.
    - `neighbour_years` is computed but never used (Dec 31 cases are covered by the delta ledger and 5b).
    - Gate 5d's "join eligibility gains = validity flips" half is not a gate. It only appears as a per-cell explanation, and that explanation needs the state tables.
    - The inventory excludes Q-invalid-year rows from classes; I saw no problem with that.
    - Gate 5c would also have caught the E7 annual allowance of 1 value difference only because the injection added a second one.

12. **NOTE: run-hygiene details.**
    - `git status --porcelain` is repo-wide, so the parent `.DS_Store` makes the tree dirty. Run from a clean HEAD worktree.
    - A changing untracked file mid-run stops the run (`dirty_tree_sha256`).
    - `merge` does not compare `fingerprint.qa_schema_sha256` with the `--qa-schema` file it is given.
    - The pinned-header check compares raw lines, so CRLF versus LF would give a misleading error.
    - `ChunkedArray.equals` treats NaN as unequal in B5's non-Q column check. The base has no NaN in the non-Q columns (and no non-null NaN in Q), so it is fine today.

## Direct answers to the six questions
1. **Splice.** The rule is faithful to the R code: status A/P, finite, not the sentinel, with `Q_new = value_today × conversion`. It acts only on masked splice-code days, and the base-NaN and flag assertions are enforced at write time. The schema, row groups and row order are preserved (the base has no non-null NaN). Metadata is updated in place with the compile's rule, which reproduces the February metadata. Absent date, non-numeric value, several series, several unqualified blocks, duplicate date, wrong unit and wrong site are each refused or marked unavailable (tested with synthetic caches). A duplicate date anywhere in a gage's record voids the whole gage, which is conservative.
2. **Oracle.** The preprocessor and qualification rules are reproduced exactly, including leap years and Oct–Dec water years. The ledger is applied only to existing NaN days, and every ledger gage must be consumed. The gaps are items 9 and 11.
3. **Gates.** B1 doesn't pin the window; B2 is aggregate only; B3 is count-only; B4 is circular; B5 and B8 are sound. B8's runtime is finding 1.
4. **Census.** NaN and NULL are handled correctly in duckdb. Class definitions match §10.5. Value validity rests on 5c, which is the weakest link (finding 2). Gates 5c–5g are present but only 5b is computed locally; the rest are given by file.
5. **Delta ledger.** Daily attribution is verified (ledger and inventory values exact, interpolates recomputed). State attribution is coarse (finding 4). The state tables contain the raw mask, interpolation mask, seasonal flags, rejection and validity (for Q only on the reason), which is enough once the checks in finding 4 are added.
6. **Batch runner.** The fingerprint is good except for the baked config, the HYDAT source and the `qa_schema` equality (findings 8 and 12). `--require-clean-tree` does what it says. The `STREAMFLOW_QA_SCHEMA` block is a no-op when every column is present, but see finding 7.

## Verdict
**GO WITH CHANGES.** Accept Phase B's data outputs once the audit runs: the splice values, ledger, metadata rule and oracle are independently verified; the patched parquet is sound as far as could be checked. Finding 1 must be fixed before the audit runs. Before Phase D, fix findings 2, 3 and 4 (2 is the main one, because 5c carries the value validation); 5, 6 and 7 are quick and should go in at the same time. Findings 8–12 are not blockers.

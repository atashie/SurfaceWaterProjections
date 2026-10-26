# Fix the USGS qualifier mask — clean up both products for their existing gages and QA/QC them against the version-2 runs

**Date**: 2026-10-21; **version 6 of the plan, 2026-10-26** — v3 (Codex rounds 1 and 2 folded
in) plus **§10: the climate-frame fix folded into this run** (user with the research team,
2026-10-26; reverses decision 5), rewritten after Codex round 3 (NO-GO on v4) and corrected
after Codex round 4 (GO WITH CHANGES on v5: the preprocessed-value ledger loses raw-mask
changes; the subset stages cannot prove universal identities and `flagged_for_high_na` depends
on the executed cohort's columns; the dependency registry must be enumerated per output kind,
family emission gates included; the GDD gate formula; conflicting instructions between §0 / §4 /
§5 and §10; §10.8 overstated; the schedule does not reach 2026-11-09 for one person — all
verified, dispositions in §9).
**Status**: decisions 1–4 and 6–11 DECIDED 2026-10-26; decision 5 reversed; §10's decisions
12–19 proposed with defaults (12–17 endorsed by the reviewer; 18–19 new in v6). Nothing has
been run.
**Companions**: `2026-10-21-qualifier-fix-expected-years.csv` (splice forecast) and
`2026-10-26-outer-join-forecast.csv` (join forecast, §10.3).
**Trigger**: the `A R` finding of the filtering census (`2026-09-04-filtering-stage-census.md`
§10.1; CHANGELOG Known Issues; STATUS). **Companion**: `2026-10-21-qualifier-fix-expected-years.csv`
— the per-gage, per-window year-level oracle the rerun must reproduce (§6).

**USER DECISION (2026-10-23): no new watershed enters the analysis.** The corrected products
keep exactly the version-2 gage sets (6,678 for WY 1993–2025, 6,250 for WY 1980–2025); the
fix only cleans up values and water years of gages already in them. Gages that would newly
qualify under the corrected input (29 / 26 for the `A R` rule, 31 / 29 if the other approved
codes are restored too; §6) stay out and are listed in the release notes. The four compile-floor
candidates of census §10 stay out. Manuscript and product counts do not change.

## 0. Summary

The R ingestion (`R/helperFunctions.R::generate_streamflow_dt`) keeps a USGS daily value only
when its NWIS qualifier string is exactly `A`, `A e`, `P` or `P e`. Every other code sets the day
to NA while the code itself is kept in the parquet's `flag` column. Four of those codes mark
**approved** values: `A R` (records revised), `A <` / `A >` (censored at a reporting limit) and
`A [4]` (the service's "no matching qualifier" placeholder). **193,542 approved days on 987 of
the 8,014 compiled gages are NA for us** (908 of them in product #1 v2, 757 in #2 v2). Treating
them as finite would make 1,370 water years valid that the preprocessor now rejects (1,197 of
them from `A R` alone), 1,177 / 950 of those inside the two products.

Two defects are fixed in one version: the qualifier mask (§1–§9) and the climate frame
(§10). Route (details §3–§5, §10): fix the rule in R; re-fetch the affected gages; **splice** the restored
days into the February 2026 parquet as a patched snapshot with a correction ledger, the raw
responses and a rule-vs-vintage control; rerun each product as a **frozen cohort** — the full patched input with a
deposited per-product selection manifest (`--gage-list`), so that a reader can reconstruct the
run exactly; gate the result against version 2 in two exhaustive classes — gages unaffected by
BOTH fixes cell-for-cell equal, every other gage (the complement: splice-affected, join-affected
or both) changed only where the dependency closure of §10.5 permits, every changed cell
reconciled to the raw splice ledger, the inserted-row inventory or the per-stage state tables,
and climate-only annual values exactly equal between the baseline and the splice stage
(S0 → S1) on every shared eligible year. Harness, validator,
rpkg-gate and staging scripts need the changes in §4 before the run. Realistic effort: about
twelve and a half person-days hands-on plus ~23 h of unattended machine time, plus a 5–7 h
rerun contingency (§8), the production batches on a quiet machine; for one person working
sequentially the 2026-11-09 submission date is NOT reachable (§8, decision 19).

## 1. The defect and its footprint (verified on the full parquet, 111,624,189 rows)

| Masked code | Days NA | Gages | in #1 v2 | in #2 v2 | Outside both | NWIS meaning (the service's code table) |
|---|---|---|---|---|---|---|
| `A R` | 181,461 | 755 | 684 | 562 | 58 | Approved; "Records for these data have been revised" |
| `A <` | 11,444 | 139 | 132 | 114 | 4 | Approved; "Actual value is known to be less than reported value" |
| `A [4]` | 419 | 101 | 98 | 86 | 3 | Approved; `[4]` = "Returned when there is no matching qualifier" |
| `A >` | 218 | 30 | 28 | 22 | 2 | Approved; "known to be greater than reported value" |
| **union, approved codes** | **193,542** | **987** | **908** | **757** | **63** | |
| `P` + remark (`Dis`, `Rat`, `Ice`, `***`, `Bkw`, `Eqp`, `Fld`, `[4]`, `<`, `ZFL`, `Ssn`, `Mnt`, `>`) | 15,257 | 263 (198 of them carry no approved code: 12,277 days) | | | | Provisional + a remark; what the agency served in February 2026 is unknowable from our files |
| **all masked gages** | | **1,185** | | | | the audit-fetch universe (B3) |
| US rows with no code and no value | 12,103 | | | | | not part of the defect |

What the agency serves **today** for the masked days (cached responses, 17 sample gages, joined
on the parquet's gage/date keys — the method of the first draft, which searched today's code
strings, was wrong): every sampled `A R` / `A <` / `A >` / `A [4]` day carries a value under the
same code; of the provisional-remark days, 1,273 of the sampled ones now carry a value, mostly
re-coded `A`, `A e` or `P e` after approval (`P Rat` 703, `P Bkw` 199, `P Ice` 136, `P ZFL` 86,
`P [4]` 74), while `P ***`, `P Dis`, `P Eqp`, `P Fld`, `P Mnt` days carry none. So "the
provisional-remark mask changed nothing" is withdrawn; the audit (B3) records, per masked day,
value / status / absence / sentinel today, and decision 1e says what to do with them.

Revisions cluster late: 92,794 `A R` days in WY 2020–2025 and 77,532 in 2010–2019, so record
ends and the last-decade trend-completeness gate are hit hardest.

Year-level effect, replayed on the full October–September calendar with the preprocessor's
rules (`julia/src/io.jl`: missing calendar days count as NA; > 30 NA or a run > 3 rejects;
boundary NAs reject after interpolation; the replay reproduces every baseline census status):

| | `A R` only | four approved codes |
|---|---|---|
| Water years made valid, all compiled gages | 1,197 in 601 gages | 1,370 in 676 gages |
| inside product #1 v2 | 1,016 in 550 gages | 1,177 in 619 gages |
| inside product #2 v2 | 810 in 444 gages | 950 in 502 gages |
| Product gages with ≥ 1 restored day in the window (#1 / #2) | 683 / 562 | **899 / 757** |
| Already-valid water years that get real values where NA or interpolated (#1 / #2) | | 500 / 403 |

The first draft's budget (`AR_gage_wy_budget.csv`) counted NAs among existing rows only, so it
had 50 false positives and 33 misses; it is superseded by the full-calendar replay.

Not affected: Canadian rows (HYDAT symbols are never masked), the EO products, the Daymet basin
series. **Restoring a day's value does not change the climate join**: the runner left-joins
climate onto the streamflow rows (`run_julia_benchmark.jl` ~l. 220) and the masked rows exist
with Q = NA, so their Daymet values were always joined; the "climate rides on the streamflow
dates" defect (absent rows) is untouched by this fix. The interaction is at the year level:
climate years are evaluated only for streamflow-valid years, so newly valid Q years add climate
years, and restored values replace interpolated or NA days in already-valid years.

## 2. Our workflow, reviewed

| Stage | Where / how | Carries the fix? |
|---|---|---|
| Ingestion | `run_ingest_usgs_hydat.R` → `process_gages_rawToRaw` → `generate_streamflow_dt`; R + dataRetrieval + tidyhydat on the **Windows** machine, February 2026 | The rule lives here (one line). Not runnable on this Mac; raw retrievals were **not retained**; no test covers the mask; `config.R` documents the opposite ("values below are set to NA"). |
| Compiled inputs | `combined_streamflow_data_09feb2026.parquet` + metadata CSV (16,994 rows; **7,046 USGS ids stored without leading zeros**, the parquet padded — every join must be leading-zero-safe, as the runner's is); verified local copies under `~/HISSS_data/daymet-replay-01oct2026/inputs/` and the staged `resource3_input_data/` | Size and footer checks only — the content was never audited against the source. A `flag × is.na(Q)` tally would have shown it. |
| Qualification oracle | `qualification_census.jl` (101 s; reproduces both v2 gage sets exactly) | Baseline for §6; the counterfactual needs the full-calendar replay (this plan's `fullcal_replay`, to become a tool). |
| Production run | Julia canonical via `run_batched_julia.py` split/run/merge (8 batches), one folder per run | Yes, with the changes in §4 (reference merge, resume keyed on content, dirty-tree refusal). **Caveat**: both v2 products record `git_working_tree_dirty = true` at `c6d63dc` with the climate family uncommitted (untracked `julia/src/climate.jl` is not in the saved diff) — v2 cannot be tied to a tree, so class U (§5) doubles as the proof that the committed tree reproduces it. |
| Gates | `validate_production_run.py`, `validate_annual_values.py`, `check_signature_failures.py`, `audit_qualification.jl`, `check_additivity.jl` | `validate_production_run` gate 5 checks overlap + lat/lon, not set equality; `check_additivity` requires identical sets (true here). |
| Output review vs predecessor | `compare_climate_rerun.py` (2026-10-12) | Its classes are wrong for this change (Q-only columns *must* move for affected gages). A sibling tool (§4 E5). |
| Ports | rpkg reference re-baselined on #1 v2; `run_rpkg_acceptance_gates.sh` hard-codes the August reference and 1,653 columns | No signature code changes → rpkg reference re-baseline only, after the gate script is parameterized. |
| Deposit | `stage_hydroshare_climate_v2.py` hard-codes "same streamflow input", "gage set unchanged", "Q-only columns identical", "every gage has climate", and silently skips missing review files | Needs a qualifier-release mode whose claims are templated from the audit and comparison outputs. |
| Docs + close-out | CHANGELOG / STATUS / census / skill / README / DATA_SOURCES / manuscript relay; commit → Sonnet 5.5 + Codex reviewers by area → fix → commit → push → `publish_to_hisss.sh` | Yes. |

Beyond the one line, this defect asks for: (1) an input-audit gate in `/run-benchmark` — the
`flag × is.na(Q)` tally recorded with every run, with the rule "no eligible finite source value
was discarded"; (2) raw NWIS responses retained with the run that used them; (3) a Mac-side,
tested compile step (the splice tool is that port for the restored days); (4) the acceptance rule
documented in `docs/DATA_SOURCES.md`, the input dictionary and a Methods sentence.

## 3. Decisions (numbered; recommended default first — **every default APPROVED by the user, 2026-10-26**)

1. **Which masked days to restore, and how to treat them.**
   - 1a `A R` (revised, approved) and estimates (`e`): restore. Default yes.
   - 1b `A <` / `A >` (censored bounds): **storage** — restore the reported bound with the code
     kept in `flag`; **analysis** — use the bound as the value (default). Evidence for the choice:
     a **numerical sensitivity run** on the 169 gages that hold such days (Julia, one small batch,
     bound-as-value vs NA): qualification differences and per-signature differences, reported in
     the audit next to the exposure table (days, gages, share of each gage's zero/low-flow days).
     Alternative: keep them NA (removes whole low-flow years from intermittent streams, e.g.
     06400875's 1,498 days).
   - 1c `A [4]` (no matching qualifier, approved): restore. Default yes.
   - 1d The ingestion rule going forward: "status letter ∈ {A, P} ∧ finite value ∧ not the
     service's no-data sentinel"; remark codes never mask at **storage**; a remark category not in
     §1 is logged and reviewed before it is used analytically. For THIS patch the restored
     categories are exactly 1a–1c. Default yes.
   - 1e Provisional-remark days (15,257; many numeric today under a new status): default **out of
     scope for this patch** — an explicit scope exception, not a claim that their February state
     differs in kind from the approved days' (February availability is unknowable for both; for
     approved days it is near-certain, for provisional-remark days it is not); the ledger records
     how many are numeric today. Alternative: restore when numeric today.
2. **Scope of the data change.** Default: **splice-only, as a patched snapshot** — re-fetch the
   987 affected gages, replace Q on the masked-code days only; the immutable February parquet, the
   raw responses and a gage/date **correction ledger** (old flag, old Q, value today, status
   today, units, conversion, retrieval time) are deposited with it; a **control of three
   datasets** — the immutable base, the October fetch under the old rule, the October fetch under
   the new rule — isolates the rule effect **at October vintage** and reports, per gage, how many
   unmasked days differ from February (not applied; §5 B2 — the cached samples show 0–4.4 % of a
   gage's days); it cannot establish what the February-masked dates held in February, and the
   README says so. The README calls the file "February 2026 retrieval with the qualifier-mask
   correction applied from an October 2026 retrieval", never reconstructed February data. Alternatives: (B) re-ingest complete records for the 987 affected
   USGS gages only (no mixed vintage inside a gage, mixed across gages, wider class-A changes);
   (C) a full re-ingest of all 16,994 candidates (after submission).
3. **New gages.** RESOLVED by the user, 2026-10-23: none. Mechanism (round-2 change): the full
   patched input plus a **per-product selection manifest** (the v2 gage list, deposited with its
   hash, the v2 reference hash and the exact command) applied by a first-class `--gage-list` in
   `run_batched_julia.py split`; the auditor distinguishes **qualified**, **selected** and
   **emitted** sets (emitted = selected, every selected gage qualified); the README describes the
   products as a frozen cohort, not "all qualifying gages", and lists the would-qualify gages
   of §6 recomputed from the realized restorations. The candidates 02201230, 07069000, 07074420, 07188653
   and 11129800 (also a mask casualty: 20 qualifying calendar years without the mask; its row in
   `docs/filtering_min_flow_rule_2026-10-21.csv` now says so) stay out.
4. **US days with no code and no value in February** (12,103): leave alone; counted in the
   fetch report.
5. **Fold in the streamflow-frame climate fix** (outer join; ON HOLD since 2026-10-16)?
   **REVERSED 2026-10-26 (user, after discussion with the research team): YES — folded into
   this run.** Design, attribution and assessment in §10; its own decisions are 12–16.
6. **Naming.** Run folders `~/HISSS_data/qualifier-rerun-<date>/processedOuts_<window>_climate_v3_<date>`;
   internally "version 3"; the deposit stays "version 2 (October 2026)" because the 2026-10-12
   build was never uploaded. Alternative: call the deposit version 3 and say why.
7. **Metadata CSV.** Update the existing rows of the spliced gages in place (`num_years`,
   `start_year`, `end_year`, `num_days` recomputed with the compile rule; `error_message`
   cleared); never append; the join key is the leading-zero-normalized id and must stay unique.
8. **rpkg reference re-baseline** on product #1 v3 (8 batches, ~3 h of R, unattended), after the
   gate script is parameterized. Default yes. Python: schema gate on one batch.
9. **Reviews.** Sonnet 5.5 + Codex, blind, by area — (a) the splice tool, ledger and input audit
   against the NWIS code table; (b) the honesty of the v3-vs-v2 comparison and the oracle; (c)
   docs, counts and the release notes; after the production batches, never alongside. Codex
   round 2 on this plan before any code is written (§9).
10. **Deposit.** Stage locally only (R1, R2, R3 updated); upload after the user's review. The
    deposited 2026 products (1,653 columns) are superseded, not rewritten.
11. **Timing.** Now — the version-2 deposit is still pending and the known issue is HIGH; the
    chain fits before 2026-11-09 only if the tooling week (§8) starts immediately.

## 4. The plan

### Phase A — the rule, its tests and its documentation (½ day)
- `R/helperFunctions.R::generate_streamflow_dt`: the decision-1 rule; `config.R`: correct the
  `MIN_Q_VALUE_AND_DAYS` comment; `R/tests/test_qualifier_codes.R` on a toy frame covering every
  code of §1 with and without a value and the no-data sentinel.
- `docs/DEVELOPMENT.md` (Data Sources line), `docs/DATA_SOURCES.md` row 1, the R3 input
  dictionary `flag` row, census §10.1 (pointer), `claude-skill/streamflow-signatures.md`.

### Phase B — tooling (1–2 days), then fetch, splice, audit (½ day, mostly download)
New `docs/benchmarks/splice_nwis_qualifiers.py`, every artefact under
`~/HISSS_data/qualifier-rerun-<date>/input/`:
- `list` — two universes from the parquet (US rows, `flag` with status A/P outside the kept
  set, Q NA): the **audit-fetch universe** (all 1,185 masked gages, needed by B3) and the
  **splice universe** (decision 1; 987 gages for 1a–1c); per-code day counts; key =
  leading-zero-normalized id.
- `fetch` — legacy NWIS DV JSON per gage (1979-10-01 … 2025-10-01, `siteStatus=all`), gzip-cached,
  resumable, retried, ≥ 1 s apart; `fetch_report.csv` (series count, days, codes seen, no-data
  days, absent dates). OGC fallback only after a 20-gage check that values, dates, series and
  codes map one-to-one, with pagination and units covered.
- `splice` — `Q_new = value × conversion` on the masked-code days that carry a finite value
  today (`conversion` exactly as R: `86400 / (basin_area × 3280.84³) × 1e6`, area from the
  metadata CSV; no un-normalized USGS gage exists); all other cells untouched; output parquet with
  identical schema and row order; `correction_ledger.csv`; `provenance.json` (sha256, `PAR1`,
  sizes, retrieval times); metadata rows updated (decision 7).
- `control` — old rule vs new rule on the October fetch for the affected gages: days and years
  that differ, split into "restored by the rule" and "revised by the agency since February".
- `audit` — gates B1–B8 (§5) → `input_audit.md` + CSVs.
- `qualification_census.jl` on the new parquet, and the full-calendar replay promoted to a tool
  (`docs/benchmarks/qualifier_year_oracle.py`) that builds the **acceptance oracle from the base
  parquet + the eligible cached source records** (realized restorations, not the forecast):
  per gage and window the baseline and final valid-year sets and counts, source-unavailable
  targets, the selected policy and the configuration hashes (§6). The forecast CSV of this
  session is kept as the forecast.
- **Checkpoint**: the audit and the oracle are reviewed (user + Codex) before Phase D.

### Phase A′ — the climate-frame fix in the three runners (§10.2; 1.5 person-days incl. blind reviews)

### Phase C — harness changes (1 day; verified gaps from the round-1 review)
- `validate_production_run.py`: `--expect-same-gages` (exact set equality + metadata equality),
  not overlap.
- `run_batched_julia.py`: `split --gage-list` (selection manifest, hashed into `batches.json`);
  `merge --reference` with identical sets needs no `--allow-missing`; resume keyed on the full
  source revision (tree hash), the Julia version + `Manifest.toml` hash, the config, metadata and
  manifest hashes and the actual batch-input hashes, with output identity re-checked; **a
  production run refuses a dirty or unknown tree** (no `--allow-dirty` for products).
- `audit_qualification.jl`: qualified / selected / emitted sets (round-2 finding 2; the
  frozen-cohort contract is `emitted = selected ⊆ qualified`).
- **Pinned signature schema for QA** (round-4 finding 2): the runner builds
  `flagged_for_high_na` over the union of the result keys of the gages it executed
  (`run_julia_benchmark.jl` l. 466–473; `run_batched_julia.py` l. 11–14 documents it), so a
  subset run whose cohort lacks a family (no snow gage, say) gets a smaller denominator and
  the flag flips for gages whose values did not change. Every stage, every subset and the E7
  port run take the full-product header (v2's 2,277 columns) as the QA denominator: a
  `STREAMFLOW_QA_SCHEMA=<v2 header>` override read by the three runners (missing columns
  added as NaN before QA); `merge` verifies every batch header against it, not only across
  batches.
- **Stage provenance**: each stage records its parent-cohort manifest (the product's selection
  manifest), its executed-subset manifest, its OWN batch-input hashes (the stage's inputs, not
  the product's) and the `STREAMFLOW_CLIMATE_FRAME` value — `batches.json` and the timing JSON
  carry all four.
- `validate_annual_values.py`, `check_annual_parquet_equality.py` and the new comparison tool
  read the annual parquets by row group / gage partition (27.6 M and 34.7 M rows; the current
  validators materialize whole tables).
- `run_rpkg_acceptance_gates.sh`: reference, candidate and log paths, expected column count
  (2,277) and the named residual allowances as arguments.
- `stage_hydroshare_climate_v2.py`: a qualifier-release mode — every README claim templated from
  the audit, oracle and comparison outputs; missing review files are an error; one manifest
  distinguishing the internal build id from the public release label.
- `compare_qualifier_rerun.py` (§5 E5) and the dependency registry (year-local vs
  record-dependent signatures, from `docs/signature_categories.csv` + the family code).

### Phase D — production (≈ 2–3 h per product on a quiet machine; unattended)
- D0 = **S0, baseline reproduction**: the tagged tree (A, A′, C committed) with
  `STREAMFLOW_CLIMATE_FRAME=left` on the ORIGINAL input, both windows, all batches (≈ 5–7 h
  unattended), compared cell-for-cell with v2 — gives a tree-tied v2 baseline (v2 itself was
  built from a dirty tree) before the input changes. **S0 ≠ v2 stops Phase D** (decision 18):
  the drift is explained and resolved to equality; replacing v2 by the reproduction as the
  comparison reference is a user decision recorded in the release notes, never a default.
- D-S1 = **stage 1**: the spliced input with `STREAMFLOW_CLIMATE_FRAME=left`, the same tagged
  tree, on the splice-affected gages + the 200-gage control sample, both windows — isolates the
  splice effect (§10.4).
- D-S1′ = **interaction cohort**: the original input with `STREAMFLOW_CLIMATE_FRAME=outer`, the
  same tree, on the composed cohort of §10.4 (the realized overlap gages incl. those whose
  restored years hold absent dates, the eligibility gains, the Dec 31 neighbours, ≥ 200
  join-only gages), both windows.
- D1 (= **stage 2**, the product) clean tree at a tagged commit, `STREAMFLOW_HASH_INPUTS=1`, the two version-2 wrappers with
  `STREAMFLOW_DATA_PATH` / `STREAMFLOW_METADATA_PATH` at the spliced files, `split --gage-list`
  with the product's selection manifest, climate `daymet_1980_2025_29sep2026.parquet` (verify 5,134,193,404 B,
  `PAR1`, md5 `c0590760…`); `split` (8 batches) → `run` (check `vm.swapusage` first; raise
  `--max-minutes` if shared) → `merge --reference <v2 CSV>` (`--memory-limit 3GB` for the 46-year
  window); read `merge.log`. Product #1 v3, its gates, then #2 v3.

### Phase E — gates, comparison, artefacts (½ day)
E1 `validate_production_run.py --reference-csv <v2> --expect-same-gages`; E2
`validate_annual_values.py`; E3 `check_signature_failures.py` on every batch log; E4
`audit_qualification.jl` (`emitted = selected ⊆ qualified`, the frozen cohort); E5 `compare_qualifier_rerun.py --new
v3 --old v2 --oracle <csv> --ledger <csv> --registry <csv>` (§5); E6 explorers +
`check_html_basemap.py`, a v3-vs-v2 dashboard, all into the run folders; `flagged_for_high_na`
comes from the run; E8 the **change census** of §10.5 on every stage pair (`assess_rerun_changes.py`); E7
**independent numerical check of the affected results**: the Python port (with its own A′
join and the pinned QA schema) run on the affected gages of product #1 (≈ 1,600 gages,
~30 min): `check_annual_parquet_equality.py` on the annual values, `check_schema_equality.py`
on header + gage set — it reads nothing else (l. 29) and proves nothing about values (round-4
finding 4) — AND a numerical gate on the subset's summary CSV: every summary statistic (16
per base), every per-gage scalar and every QA flag compared cell for cell at the port
tolerances (annual values 1e-9 relative, summaries 1e-6, the named residuals of
`docs/CROSS_LANGUAGE_STATUS.md` as explicit waivers) with a gate-state comparison (a cell NaN
in one language and finite in the other is a failure, not a residual);
`compare_python_vs_golden_julia.py` is promoted to a gate for this (today it intersects
columns and exits 0 on a missing family). The diagnostics of E5 (vi–vii) are accepted only
with this gate PASS and a written disposition.

### Phase F — rpkg reference re-baseline (≈ 3 h of R; decision 8) with the parameterized gates.

### Phase G — adversarial reviews (1 day; decision 9), fixes, re-gates; a rerun if a fix touches the input.

### Phase H — staging, docs, close-out (1 day)
R1/R2 via the staging tool's new mode (selection manifests + hashes + the exact commands in
the run notes); R3: the **immutable February base**, the patched parquet, the ledger, the raw
responses (gzip), the three-dataset control report, the final acceptance oracle, the selection
manifests and an executable reconstruction recipe, updated metadata, README, abstract,
dictionary — the staging check verifies every one of them is present. Docs (§7); every number
the manuscript quotes diffed against v3 (the freeze keeps the counts; a signature statistic may
move) and the changed ones relayed (decision 19). Commit →
reviews (G) → commit → push → `publish_to_hisss.sh` from a clean checkout.

## 5. Acceptance criteria

**Input (Phase B audit; all must pass before Phase D):**
- B1 Fetch completeness: every listed gage fetched; exactly one primary `00060 / 00003` series;
  daily keys unique; window coverage reported; failures named.
- B2 Conversion check against the cached records and the metadata: for every stored non-NA day,
  `value_today × conversion` vs `Q` in normalized units — the share of days equal to 1e-9
  relative must be the overwhelming majority and the equal days prove area and conversion (the
  February raw values were not retained, so "identical in raw units" is not claimable);
  the differing days are the source-difference report (0–4.4 % per gage in the samples;
  08286500 711 of 16,113, 01381800 227 of 9,743) — reported per gage, NOT applied (decision 2),
  no threshold.
- B3 Restoration accounting per code over the **audit-fetch universe** (1,185 gages): masked
  days with a finite value today / sentinel / absent; totals against §1; the provisional-remark
  days tabulated the same way (decision 1e).
- B4 Restored values verified against the exact cached record (site, date, statistic, unit) and
  the ledger; plausibility checks (within the gage's compiled range widened by 10 %; ratio to
  the ±3-day neighbours inside the gage's own p0.1–p99.9) are **diagnostics with a written
  disposition**, not gates; negatives are allowed where the config allows them.
- B5 Cell-level diff of the whole parquet: schema, dtypes, row count, per-gage row counts, row
  order, every non-Q column identical; `Q` differs only on ledger rows; `PAR1`; sha256 recorded.
- B6 Metadata: rows updated in place, normalized ids unique, every parquet gage joins.
- B7 Oracle: the census on the new parquet gives, per product gage, exactly the acceptance
  oracle's final valid-year sets and counts (§6); zero drops; the would-qualify gages outside
  the products equal the recomputed lists. Deviations are failures to explain **and resolve**,
  not to waive.
- B8 "No eligible finite source value discarded": a daily **source / ledger / output join** —
  every US day with status A/P and a finite value in the cached response, in the splice
  universe, has that value in the output; every other cell is unchanged.

**Run (Phase E):**
- E1–E4 as in `/run-benchmark`: gates PASS, gage set **equal** to v2, annual validation PASS, 0
  swallowed failures, audit = runner.
- E5 two exhaustive, disjoint classes derived from the realized ledgers and state tables of
  §10.5 (a gage is affected when, inside the window, it has ≥ 1 raw splice-ledger day whose Q
  changed, ≥ 1 inserted climate row, or ≥ 1 changed per-stage state entry — raw mask, seasonal
  flag, rejection reason), compared as **keyed cells** (gage × column; gage × signature ×
  water year) with NaN == NaN, null == null, signed zeros equal:
  - **U — unaffected by BOTH fixes** (forecast 5,069 / 4,893 gages; final membership from the
    realized ledgers and state tables): every summary cell and every annual-parquet row equal
    to v2 (via S0). Any difference FAILS (D0 tells whether the tree or an input is the cause).
  - **A — the complement of U** (forecast 1,609 / 1,357), reported in three sub-classes: A-Q
    splice-only, A-J join-only, A-QJ both (the overlap, §10.4). The criteria below apply to
    the sub-classes with a splice change (A-Q, A-QJ; 899 / 757 if every targeted day
    restores); A-J gages are gated by the §10.5 census alone (every change permitted by the
    join closure, precipitation totals by summed deltas, everything else by recomputation):
    (i) valid-year identity per gage = v2's years ∪
    the oracle's restored years, exactly; (ii) **immutable site attributes** (latitude,
    longitude, basin_area, area_normalized, gage_type, the 11 human-interference columns)
    unchanged, while the **derived qualification metadata** (`num_water_years`,
    `start_water_year`, `end_water_year` — the runner computes them from the qualifying years)
    equal the oracle; (iii) on water years with no
    ledger day, no state-table change and unchanged validity, **year-local** signatures equal
    cell for cell unless a family emission gate flipped for the gage (§10.5 item 3; those keys
    are classed "gate-emitted old year" and validated by recomputation, never as deltas);
    **record-dependent** families (drought thresholds, `*_all` pulses, `elasticity_rolling`,
    parameterized BFI, `avg_storage`, per the registry) reported with R² and max |Δ| on those
    years; (iv) **climate-only** annual values (the precipitation and temperature families and the
    SWE bases, `swe_max_to_ppt` excluded — it is coupled) equal on every shared valid year
    including years with restored Q, **between S0 and S1** (the splice stage; §10.4) — an
    inequality there means the splice moved the climate frame; (v) newly valid years counted per gage = oracle;
    (vi) summary statistics move by construction — reported per column with R², max |Δ| and
    its gage, gained / lost cells, `flagged_for_high_na` flips; (vii) Q-dependent
    climate-coupled signatures (runoff ratios, elasticity, Q–P seasonality, storage, snow) on
    years with restored days reported separately from (iii). (vi) and (vii) are diagnostics:
    accepted only with E7 PASS and a written disposition; the family registry names the snow
    bases explicitly (several carry no `swe_` prefix).
- E6 explorers load their sidecars headlessly; basemap PASS.
- E8 the change census of §10.5 with its gates (unexplained changes = 0; totals equal the
  ledger sums, every other base equals recomputation; stage endpoint identities).
- F rpkg: schema gate (2,277 / 6,678), failure scan, annual gate at 1e-6 with the named residuals.

## 6. Oracle (year level; gage sets frozen)

Two artefacts. (1) The **forecast** of this session,
`docs/plans/2026-10-21-qualifier-fix-expected-years.csv` — one row per affected product gage and
window: masked days per approved code inside the window, water years touched, of which already
valid, and the water years that become valid under `A R` only and under the four approved codes
(with the year lists; 1,656 rows, every one reconciled by Codex against the replay), assuming
every targeted day restores. (2) The **acceptance oracle**, generated in Phase B from the base
parquet + the eligible cached source records for the decided rule: baseline and final valid-year
sets and counts per gage and window, source-unavailable targets, policy and configuration hashes
— frozen for B7 and E5.

Gages that WOULD qualify under the corrected input and are kept out (decision 3; for the
release notes) — `A R` rule, WY 1993–2025 (29): 01196561, 01454700, 01480638, 01527500,
01542500, 01581752, 01585200, 01589290, 021473428, 02312722, 02317620, 02369600, 03250322,
03292475, 03508050, 04180988, 04250750, 06917630, 07056700, 07076517, 072632962, 07325850,
07360200, 07362000, 07363200, 08330600, 08397600, 10308800, 12210700; WY 1980–2025 (26):
01160350, 01349150, 01349810, 01585200, 0208735012, 02087359, 02146530, 02153200, 02319800,
02320000, 02423630, 02483500, 03028500, 06131200, 06935770, 07363400, 09301500, 09419800,
10157500, 10243260, 10251300, 10308800, 11109000, 11180900, 11261500, 13140800. With the four
approved codes add 02274005 and 12202310 (short window) and 07010030, 07299890, 09400568 (long
window). Thirteen / fifteen of them sit exactly at the gate edge (20 years, or fraction 0.606 /
0.609). These are forecast lists; the release notes carry the lists recomputed from the realized
restorations.

## 7. Documentation inventory (version-aware; dated history is never rewritten)

Rule and input description: `docs/DEVELOPMENT.md`, `docs/DATA_SOURCES.md`, R3 README /
abstract / dictionary, `README.md`, `claude-skill/streamflow-signatures.md`,
`.claude/skills/run-benchmark/SKILL.md` (input-audit gate), `.claude/rules/benchmarks.md`.
Known issue closure and counts of restored years: CHANGELOG (Fixed HIGH + the user decision),
STATUS (product table v3 rows, issue line removed), census §10 / §10.1 and
`docs/filtering_min_flow_rule_2026-10-21.csv` (five mask casualties, not four),
`docs/filtering_stage_census_2026-09-04.csv` (pointer only), `docs/SIGNATURES.md` caveats that
cite NA years, `docs/CROSS_LANGUAGE_STATUS.md`, `rpkg/README.md`, `python/README.md`
(reference product names). Unchanged by the freeze and only to be cross-referenced:
`docs/plans/2026-09-04-dataset-join-guidance.md`, `EO_data_processing/{CLAUDE.md, README.md,
README_DAYMET.md, geometry/README.md}` (7,964 polygons, 8,014 compiled), the manuscript counts.
Co-author relays (the Google Docs cannot be edited here): the guidelines doc states the old
mask ("USGS values are kept only for qualifiers A, A e, P, P e", `docs/SIGNATURE_GUIDELINES.md`
l. 540 of the snapshot) — a doc-side item in `docs/reconciliation/guidelines_todos.md`; the
manuscript needs a Methods sentence on qualifier handling, the restored-year counts and the
frozen-cohort statement (the funnel sentence is already relayed).

## 8. Effort

| Phase | Hands-on | Unattended | Depends on |
|---|---|---|---|
| A rule, tests, docs | 0.5 | — | decisions 1 |
| B tooling + fetch (1,185 gages) + splice + audit + sensitivity run + oracle + checkpoint | 2 | 0.5 d download | A |
| C harness (`--gage-list`, resume fingerprint, dirty refusal), auditor, partitioned validators, rpkg gate, staging mode, comparison tool | 1.5 | — | — (parallel with B only with a second person) |
| D baseline reproduction + two products | 0.5 | 10–14 h quiet machine | B checkpoint, C |
| E gates, comparison, E7 port check, explorers | 0.5 | 1–2 h | D |
| A′ climate-frame fix: shared frame helper + `STREAMFLOW_CLIMATE_FRAME` switch in Julia, Python, rpkg; fixtures (interior years, Dec 31 neighbours, eligibility flips, window inside the span); cross-language fixture; blind reviews | 2 | — | — |
| D-S1 stage 1 on the splice-affected gages + control sample, both windows (30–45 min each); D-S1′ on the composed cohort (≈ 20 min per window) | 0.3 | 1.5–2.5 h | B checkpoint, C, A′ |
| B′ daily delta ledger + per-stage year-set tables (library replay on the affected gages) | 0.5 | 1 h | B, A′ |
| E8 change census tool + reports + dashboard | 1.5 | 1 h | D |
| Round-4 additions: raw masks + per-stage state tables (B′); registry enumerated per output kind incl. the family emission gates (E8); summary-CSV numerical gate + gate-state comparison (E7); pinned QA schema + stage provenance (C) | 1 | — | B′, C, E |
| F rpkg re-baseline | 0.2 | 3 h | E |
| G reviews + fixes (+ a rerun if the input changes) | 1 | (rerun: 5–7 h) | E, machine quiet |
| H staging, docs, relay, commit, publish | 1 | — | G |

**The table adds to 12.5 person-days hands-on (v5 said "about 11" over a table that summed
to 11.5) plus ~23 h of unattended machine time (D 10–14 h, D-S1 / S1′ 1.5–2.5 h, B′ 1 h, E
1–2 h, E8 1 h, F 3 h; B's half-day download overlaps with C), plus a 5–7 h rerun contingency
if a review fix touches an input or the join** (v3 said 7 / 20 before the climate-frame fix;
v4's 10 / 28 ran stage 1 in full). **Dated critical path, one person, sequential, from
Monday 2026-10-26**: A + A′ Oct 26–28 · B tooling Oct 28–30 (fetch overnight Oct 29) · C
Nov 2–3 · B splice / audit / oracle + checkpoint Nov 3–4 · D0 overnight Nov 4 · D-S1 / S1′
Nov 5 · products overnight Nov 5 and Nov 6 · B′ + E + E8 Nov 9–11 · F Nov 11 · G Nov 12 · H
Nov 13 → **earliest completion 2026-11-12/13 with no rerun, about 2026-11-16 with one**. There
are 10 weekdays between 2026-10-26 and the 2026-11-09 submission (11 counting the submission
day): **one person cannot finish by 2026-11-09**. With a second person (B ∥ C and A′, saving
≈ 3.5 days) and the machine running every night and the weekend of Oct 31, the earliest
completion is Friday 2026-11-06 with zero slack and no rerun. The frozen gage sets keep every
COUNT in the manuscript unchanged; a signature STATISTIC the manuscript quotes (a mean, a
median, a share of significant trends) can move with the restored years and the corrected P,
so Phase H diffs every manuscript number against v3 and relays the changed ones (decision 19).
The first draft's "two to three" did not budget the tooling; v2's "four to six" did not add
the table up; v5's "about 11" rounded a 11.5 down and claimed the submission date.

## 9. Review history

**Codex round 1 (2026-10-23, gpt-6-astra, high reasoning, read-only) — NO-GO as written; 14
findings.** Dispositions after independent verification (all numbers re-derived here):

| # | Finding | Disposition |
|---|---|---|
| 1 | Oracle counted NAs among existing rows, not the full calendar (50 false positives, 33 misses) | ACCEPTED; full-calendar replay reproduces every census status; year-level oracle rebuilt (§1, §6). |
| 2 | "Provisional-remark days carry no value today" false — 1,273 numeric in the samples | ACCEPTED; claim withdrawn; decision 1e added; B3 records every masked day's state today. |
| 3 | Candidate list misses 11129800; appending metadata rows would duplicate keys | ACCEPTED (label fixed; rows are updated in place); the candidate question itself is MOOT — user decision 2026-10-23, no new gages. |
| 4 | 1,186 gages mixed approved codes with the provisional population; union is 987 (908 / 757 / 63); provisional 15,257 days / 263 gages (CHANGELOG said 15,460) | ACCEPTED; corrected everywhere. |
| 5 | Restoring Q does not restore Daymet join dates; attribution and the "few percent" tolerance wrong | ACCEPTED; §1 and E5 (iv) rewritten; climate-only annual values must be exactly equal. |
| 6 | "Every numeric A/P value" bundles censoring with approval; B4 negatives clause contradicts the config | ACCEPTED in part: decision 1 split into 1a–1e with a storage/analysis distinction and a sensitivity table for `<` / `>`; negatives follow the config. The default still restores the bounds as values — the alternative is recorded. |
| 7 | Splice-only is not pure attribution; third option (affected gages re-ingested) missing; needs base + raw + ledger | ACCEPTED; decision 2 restated as a patched snapshot with ledger, raw responses and the rule-vs-vintage control; alternative B listed. |
| 8 | B2's 0.1 % gate fails on real gages (4.4 %); B4/B5 do not prove restored values | ACCEPTED; B2 split into conversion identity on identical values + an unthresholded difference report; B4 diagnostics with dispositions; restored values verified against the cached record. |
| 9 | Class A = gages gaining years (552 / 447) is incomplete; affected = 899 / 757 | ACCEPTED; classes derived from the ledger, exhaustive and disjoint; registry of year-local vs record-dependent signatures. |
| 10 | Harness: `--allow-missing`, overlap-not-equality gate, hard-coded rpkg gates, resume on cmd/env, dirty trees permitted; v2 built dirty | ACCEPTED; Phase C; D0 control batch; v2 provenance caveat recorded. |
| 11 | "Bitwise" undefined; batch reassignment harmless | ACCEPTED; keyed-cell equality with NaN / null / signed-zero semantics. |
| 12 | Staging would publish false hard-coded claims | ACCEPTED; qualifier-release mode in Phase C. |
| 13 | Doc inventory incomplete; blanket replacement would corrupt dated history | ACCEPTED; §7 rewritten, version-aware. |
| 14 | Schedule excluded tooling and review loops | ACCEPTED; §8 is now four to six days with a checkpoint before production. |

**Codex round 2 (2026-10-23, delta review of v2) — GO WITH CHANGES; 8 new findings; round-1
items 3, 4, 5, 11 RESOLVED, the rest PARTLY.** Independent re-derivation by the reviewer: all
1,656 oracle rows, the §1 totals, the §6 lists, 259 sampled gage-year replays and a Julia
experiment confirming E5 (iv) (53 shared climate annual values, zero differences, after a
restored Q year). Dispositions in this version:

| # | Finding | Disposition |
|---|---|---|
| 1 | E5 (ii) "metadata unchanged" rejects the intended result — the runner derives `num_water_years` / `start` / `end` from the qualifying years | ACCEPTED (verified `run_julia_benchmark.jl` l. 404–406); immutable attributes vs derived qualification metadata split; affected = actual Q change; diagnostics need E7 + disposition. |
| 2 | The freeze needs a deposited selection contract; the auditor equates qualification with inclusion | ACCEPTED (verified `audit_qualification.jl` l. 69–80); `--gage-list` + manifests + qualified / selected / emitted; "frozen cohort" wording. |
| 3 | The oracle is a forecast, not a source-backed gate; CSV lacks final year sets | ACCEPTED; forecast kept, acceptance oracle generated from base + cached records (§6). |
| 4 | Fetch universe (987) cannot support the B3 audit (198 provisional-only gages, 12,277 days); B2 invokes unavailable raw values; B8 not computable from a tally | ACCEPTED (verified: 198 / 12,277 / 1,185); audit-fetch vs splice universes; B2 and B8 restated; the control's limit stated. |
| 5 | Decisions 1b–1e blur policy and evidence | ACCEPTED; 1b gets a numerical sensitivity run; 1d scoped; 1e an explicit scope exception. |
| 6 | Resume fingerprint and D0 too weak; dirty trees | ACCEPTED; full fingerprint, production refuses dirty trees, D0 = full baseline reproduction of v2. |
| 7 | Validators materialize whole tables; rpkg gate paths; Phase H inventory | ACCEPTED; partitioned readers, all paths parameterized, inventory completed. |
| 8 | Effort table adds to 6.2–6.7 person-days; guidelines l. 540 and the CHANGELOG known-issue bullet still state the old numbers | ACCEPTED; §8 restated in person-days; guidelines relay added to §7; CHANGELOG bullet corrected 2026-10-26. |

**Decisions taken 2026-10-26** (user): the defaults of 1a–1e, 2, 6, 7, 8, 9, 10 and 11, with
3 (no new gages) and 4 (no-code days untouched) as recorded; decision 5 reversed the same day
(research team): the climate-frame fix is folded in (§10). **Codex round 3 (2026-10-26, plan v4) — NO-GO as written; 8 findings, all verified here
(the interior-year counts 3,562 / 5,116 across 759 / 894 gages with 560 / 677 in class U; the
418 up / 501 down split of the October-16 deviations; the 695 `_pval` columns; the second
Julia join inside the legacy branch; `season_excluded_years_*` counted over every year of the
loop):**

| # | Finding | Disposition |
|---|---|---|
| 1 | The span-only outer join adds interior water years with no streamflow (3,562 / 5,116 gage-years) — they reach the year loop, are rejected, and their seasonal flags change the Q-only `season_excluded_years_*` scalars; `total_possible` can move when a window ends inside the span | ACCEPTED; the frame helper restricts added rows to the original WATER-YEAR SET, not only the date span (§10.2); fixtures for both cases. |
| 2 | The runner's second join is in the deprecated legacy branch, not production; the broad `catch` would hide helper failures | ACCEPTED; production = the Phase-2 join only; legacy branch left as is and documented; the `catch` narrowed. |
| 3 | The ledger must cover secondary changes (the Daymet-missing Dec 31 interpolate moves when a neighbour is restored; 3 per window) and nonlinear bases cannot be gated by summed deltas | ACCEPTED; the ledger is every changed preprocessed daily value; totals gated by sums, everything else by recomputation (§10.5). |
| 4 | §10.8: keep full-frame P, but drop the direction gate and the "residual is small" claim (P moves both ways: 418 up / 501 down; storage is not monotone in annual ΔP; family-specific eligibility) | ACCEPTED; §10.8 rewritten; decision 17 default kept. |
| 5 | Stage trees underspecified; the interaction check mis-defined; a cheaper rigorous S1 exists | ACCEPTED; one tree with a fingerprinted `STREAMFLOW_CLIMATE_FRAME` switch; S1 on the splice-affected gages + a control sample; S1′ on a selected batch; endpoint identities on non-overlap gages, the overlap residual reported (§10.4). |
| 6 | "0 unexplained changes" needs an output-level dependency closure; per-stage year sets; §0 / E5 contradictions | ACCEPTED; §10.5 (closure by output kind; "permitted to change" ≠ validated); E5 (iv) now S0→S1 only, `swe_max_to_ppt` excluded; U = unaffected by both fixes. |
| 7 | Forecast table: long-window decades mislabeled; inserted rows 3,183 / 2,195; two eligibility gains explain 1,174 vs 1,172; restored years with absent dates (26 / 59 and 15 / 26) | ACCEPTED; §10.3 corrected. |
| 8 | Census details (tolerances at zero, R² definition, 695 `_pval` columns of 5 types, paper-facing views, Q-only population means may shift from the splice, deposit inventory, headline effort, decision numbering) | ACCEPTED; §10.5, §10.7, §0, §8 corrected. |

**Codex round 4 (2026-10-26, delta review of v5) — GO WITH CHANGES; 7 new findings; round-3
items 1, 2, 7 RESOLVED, 3–6 and 8 PARTLY.** The reviewer reproduced both forecast CSV sections
row for row, the restored-year gaps (26 / 59 and 15 / 26), the 695 `_pval` columns and the
schedule arithmetic, and ran Julia fixtures for the year-support invariance — confirming that
the span-AND-year-set restriction leaves `years`, `total_possible` and
`season_excluded_years_*` unchanged, so the Q-only question of round 3 is settled —, the
raw-mask flip, the subset QA denominator, GDD and the recession emission gate. All seven
verified here: `io.jl` l. 326–338, `run_julia_benchmark.jl` l. 466–473 with
`run_batched_julia.py` l. 11–14, `recession.jl` l. 515, `qp_seasonality.jl` l. 65, `snow.jl`
l. 285–298, `climate.jl` l. 185–209 / 476 / 509, `check_schema_equality.py` l. 29; the extra
overlap gages reproduced from the census (24 / 15 gages carry restored years with absent dates,
18 / 13 of them outside the join forecast); the §8 table summed to 11.5; 10 weekdays from
2026-10-26 to the submission.

| # | Finding | Disposition |
|---|---|---|
| 1 | The preprocessed-value ledger loses causal information: a restored observation equal to its interpolate changes no preprocessed value yet flips a RAW seasonal completeness flag (80 % of raw days per season, `io.jl` l. 326); preprocessed missingness cannot identify interpolated days | ACCEPTED; §10.5 items 1–2: three ledgers (raw splice, inserted rows, preprocessed deltas) plus per-stage state tables (raw masks, interpolation mask, seasonal flags, rejection reasons, rejected years included); affected membership and gate 5a derived from all of them. |
| 2 | Subset stages cannot prove the stated identities: `flagged_for_high_na`'s denominator is the union of the EXECUTED gages' result keys (a cohort lacking a family flips the flag on unchanged gages); S1 / S1′ do not contain "every gage without a splice change"; stage input hashes differ; the interaction cohort misses the 18 / 13 gages whose restored years hold absent dates | ACCEPTED; Phase C pins the full-product header as the QA denominator for every stage, subset and port run; §10.4 gates identities on executed gages only, defines the full-cohort S1 as executed rows ∪ S0 rows with INFERRED labels, records parent / subset manifests and per-stage hashes; the S1′ cohort is composed from the realized support (165 / 135 overlap gages). |
| 3 | The registry is incomplete at the output level: the recession minimum-event gate (record-wide, before collection) makes old years' keys appear; Q–P seasonality record floors; the snow and frost anchored gates touch the trend statistics only; trend completeness uses each metric's own span | ACCEPTED; §10.5 item 3 enumerated per output kind (annual values, key emission, summaries, scalars, flags); "gate-emitted old year" class; an anchor change never permits a shared-year snow / frost annual value to move. |
| 4 | Gate 5b's GDD formula is wrong (`gdd_cum_*` sums thresholded daily degree-days); 5c leaves the totals of newly eligible years unchecked; 5d counts gains without validating them; E7 + `check_schema_equality.py` (header + gage ids only) validate no scalar or summary value | ACCEPTED; 5b per-day thresholded formula; 5c covers every changed or gained value incl. totals; 5d validates gains under 5c; E7 gains a summary-CSV numerical gate with gate-state comparison and named residuals (gate 5g). |
| 5 | Conflicting instructions: §0 (oracle-named years, universal climate equality), E5 (join-only gages in neither class), Phase D (pre-A′ tree, "batch 1"), E4 (`qualified = selected = emitted`), D0 vs §10.4 on baseline drift | ACCEPTED; §0, Phase D, E4, E5 rewritten around the final contract (A = complement of U with sub-classes A-Q / A-J / A-QJ; `emitted = selected ⊆ qualified`; one switched tree; composed cohort); decision 18: drift stops execution. |
| 6 | §10.8 overstates: "affected gage-years only" (record-normalised families, emission gates), "old P interpolated on every interpolated-Q day" (absent-row days only), storm exclusion, exposure ≠ residual, a ratio of means ≠ a mean ratio | ACCEPTED; §10.8 rewritten with the narrower claims; the complete-Q-series statement and decision 17 kept. |
| 7 | §8: the table sums to 11.5, not 11; S1 + S1′ exceed the stated 1–1.5 h; the download is outside the unattended total; only 10 weekdays before 2026-11-09; frozen gage counts do not freeze the manuscript's signature statistics | ACCEPTED; §8 restated at 12.5 person-days with a dated critical path (one person: 2026-11-12/13 at the earliest); decision 19; Phase H diffs every manuscript number. |

The reviewer's own gate assessment: every gate is computable below 3 GB by gage-partitioned
replay and column-wise accumulation once the raw masks / state, the complete registry and the
stage / schema provenance exist — which v6 specifies. Verbatim report:
`2026-10-21-qualifier-reviews/codex_plan_review_round4.md`.
Implementation order: A and A′ → B (with its checkpoint) and C → D0 → D-S1 and D-S1′ →
D (stage 2) → B′ (ledgers, state tables) → E (+ E8) → F → G → H.

## 10. The climate-frame fix, folded in (user with the research team, 2026-10-26; v5 after Codex round 3)

### 10.1 The defect (CHANGELOG Known Issues, MEDIUM, 2026-10-16; `2026-10-16-precipitation-review.md` §5)

The production runners join each gage's Daymet series onto that gage's streamflow ROWS before
the library builds the water-year grid (`run_julia_benchmark.jl` l. 220 — the second join at
l. 332 belongs to the deprecated legacy branch and is not production; `run_python_benchmark.py`
l. 319; `run_rpkg_benchmark.R` l. 276; all left joins). A date with no streamflow row reaches
`preprocess_daily_data` without PPT / TMIN / TMAX / SWE, the grid makes it NaN in every series
and the ≤ 3-day interpolation fills it — climate that exists in the input is replaced by an
interpolate. The library is correct for the frame it receives; only the runners' frame
construction is wrong. Reach: the precipitation and temperature families, the snow family, and
every Q-to-PPT signature, in every product and in both ports.

### 10.2 The fix (decisions 12–13)

One shared helper per language, used by the production join, selected by a runtime switch
`STREAMFLOW_CLIMATE_FRAME = left | outer` that is recorded in the timing JSON provenance and in
`batches.json` (so every stage of §10.4 runs on ONE tree): `frame_with_climate(gage_df,
gage_climate, mode)` =
1. candidate climate rows = the gage's Daymet rows whose date lies inside the streamflow record
   span AND whose water year is in the gage's ORIGINAL set of water years with streamflow rows
   (the year-set restriction is what Codex round 3 showed to be necessary: the February
   records hold 3,562 / 5,116 whole interior water years with no streamflow at all in the
   two windows, on 759 / 894 product gages; a span-only join would add them to the
   preprocessor's year loop, where they are rejected but their seasonal flags still change
   the Q-only scalars `season_excluded_years_*`, and a window ending inside such a hole would
   move `total_possible`);
2. outer join on `date` (Julia `outerjoin`, pandas `how="outer"`, data.table `all = TRUE`);
3. added rows get `gage_id`, `date`, `water_year` (and the other temporal columns via the
   language's `add_water_year_columns`, so the frame is uniform), `Q = missing`,
   `flag = missing`. The production year loop takes only `date`, `Q` and the climate columns
   from the frame and derives `month` / `dowy` from the grid (`io.jl` l. 297), so nothing else
   is read; the existing ice diagnostic that never fires (input `flag` dropped before it) is
   left exactly as it is.
4. `mode = left` reproduces today's frame bit for bit (the switch exists for attribution, not
   as a product option).
Also: the Julia join's broad `catch` (l. 218–223) is narrowed to "climate not available for
this gage", so a helper failure surfaces instead of masquerading as missing climate. The
legacy branch keeps its left join (deprecated, not production; documented).

What the fix does NOT touch: rows that exist with Q NA (the masked days of §1 — their climate
was always joined), the Daymet 365-day calendar (Dec 31 of leap years stays absent in Daymet,
NaN on the grid and interpolated, documented in `docs/SIGNATURES.md` §14 — but its interpolate
CHANGES when a neighbouring Dec 30 / Jan 1 was an absent streamflow day and is now filled,
§10.3), the library, the gage sets.

Tests and reviews (Julia first, then the ports; "a change is done only when all three
agree"): (i) a fixture that drops streamflow rows on wet days — `prcp_annual` equals the
Daymet water-year sum (with the documented Dec 31 rule), every Q-only signature and scalar
bit-identical to the left-frame run; (ii) a fixture with an absent interior water year and
one with the window ending inside the record span — valid years, rejected years, seasonal
flags, diagnostics, `season_excluded_years_*` and `total_possible` identical between modes;
(iii) the Dec 31-neighbour fixture (the interpolate changes by the expected amount) and the
eligibility-flip fixtures (a 30-absent-day year plus Dec 31 that was climate-invalid becomes
valid; a 4-day climate gap closes); (iv) the cross-language fixture of 2026-10-12 extended
with these cases: annual parquet equality across the three runners in both modes; (v) blind
Codex + Sonnet reviews per language (decision 9 extended to A′); (vi) the smoke test and one
batch in both modes (`check_additivity.jl` against v2: zero changes in `left` mode).

### 10.3 Footprint forecast (census of 2026-10-26; Q-valid gage-years in the v2 products with ≥ 1 absent streamflow date; reviewer-reproduced row for row)

| | WY 1993–2025 | WY 1980–2025 |
|---|---|---|
| Q-valid gage-years in the product | 200,834 | 254,195 |
| with ≥ 1 absent date (join-affected) | 1,174 (0.58 %) | 938 (0.37 %) |
| gages (all USGS; no Canadian gage has absent dates) | 857 | 722 |
| absent dates / of which Daymet-missing Dec 31 / climate rows inserted | 3,186 / 3 / 3,183 | 2,198 / 3 / 2,195 |
| 1 / 2 / 3 absent dates in the year | 625 / 211 / 125 | 526 / 159 / 102 |
| by decade 1980s / 1990s / 2000s / 2010s / 2020s | — / 37 / 90 / 556 / 491 | 14 / 47 / 59 / 419 / 399 |
| join-affected ∩ splice-affected gages (a Q-valid year with absent dates AND a splice day) | 147 | 122 |
| gages whose splice-RESTORED years hold absent dates (realized overlap; 18 / 13 of them are outside the join forecast because their only absent dates sit in years v2 rejected) | 24 (18 new) | 15 (13 new) |
| interaction cohort for S1′ (union; final membership from B′) | 165 | 135 |
| affected by either fix | 1,609 | 1,357 |
| **unaffected by both (class U, forecast)** | **5,069** | **4,893** |
| splice-restored years that also hold absent dates (reviewer-derived, reproduced from the census; confirmed in B′) | 26 years / 59 dates / 24 gages | 15 / 26 / 15 |
| whole interior water years with no streamflow (must NOT enter the frame) | 3,562 (759 gages) | 5,116 (894) |

`2026-10-26-outer-join-forecast.csv` lists every join-affected gage-year (window, gage, water
year, streamflow rows, absent dates). Reconciliation with the 2026-10-16 check (1,172
gage-years, 919 deviating `prcp_annual` values): the two extra gage-years are ELIGIBILITY
gains, not bookkeeping — 0214269560 WY 1997 (30 absent dates plus the Daymet-missing Dec 31
exceeded the climate NA limit) and 02484000 WY 2017 (a 4-day climate gap) become climate-valid
under the new frame; product #2 gains one year in each climate family. Of the 919 deviating
values, 418 rise and 501 FALL under the correction (an interpolate between two wet neighbours
over-states a dry absent day), so there is no universal direction. Three existing-date climate
values per window change through the Dec 31 interpolate (e.g. 07285500, 2016-12-31, PPT
0.333 → 2.797 mm) — the ledger must carry them (§10.5). Worst gages by absent dates: 02146285
(373 over 22 years; `prcp_annual_mean` 892 → ≈ 1,182 mm before any splice-restored year enters
its mean), 02329342 (136 / 23), 0214297160 (81 / 15), 12137800 (62 / 10), 11363930 (52 / 13).

### 10.4 Attribution: stages on one tree (decision 14)

| Stage | Input | `STREAMFLOW_CLIMATE_FRAME` | Gages | Isolates |
|---|---|---|---|---|
| S0 = D0 | February parquet | left | all (both windows) | the tree-tied baseline (must equal v2) |
| S1 | spliced parquet | left | the splice-affected gages + a 200-gage control sample of unaffected gages | the splice effect = S1 − S0; the control sample must equal S0 exactly |
| S1′ | February parquet | outer | one deliberately composed cohort: every overlap gage of the REALIZED support (the 147 / 122 join ∩ splice gages plus the 18 / 13 gages whose restored years hold absent dates = 165 / 135; final membership from B′), the eligibility-gain and Dec 31-neighbour cases, and ≥ 200 join-only gages | the join effect alone; on the executed join-only gages S1′ must equal S2 exactly |
| S2 = v3 (the product) | spliced parquet | outer | all (both windows) | the product |

All four stages run on the same clean, tagged tree (Phases A, A′ and C committed), with the
same config and metadata, the pinned QA schema (Phase C), and per stage its parent-cohort
manifest, its executed-subset manifest and its OWN batch-input hashes recorded (the hashes
differ by stage and are never implied equal). **Measured identities are gated on executed
gages only**: S0 = v2 on every gage (a full run); S1 = S0 exactly on the 200 executed control
gages, and S1 ≠ S0 only where the splice closure permits on the executed affected gages;
S1′ = S2 exactly on the executed join-only gages of the cohort; S0 = S1 = S2 on class U (S0 and
S2 are full runs — U is measured, not inferred). Where a comparison needs S1 on a gage S1 did
not execute, the **full-cohort S1 is defined as the executed S1 rows ∪ the S0 rows of every
other gage**, and every such cell is labelled INFERRED (equal to S0 on the control-sample
evidence, not by execution) in the census. On the interaction cohort the residual
(S2 − S1) − (S1′ − S0) is REPORTED per cell with an absolute-plus-relative tolerance scaled to
the cell's magnitude, with gained / lost cells listed separately — it is not required to be
zero (the two effects interact through eligibility and interpolation neighbours). Cost: S1 ≈
30–45 min per window (1–1.5 h for both), S1′ ≈ 20 min per window, instead of v4's full second
product pair.

### 10.5 Systematic assessment of every change against the previous run (B′ + E8, `docs/benchmarks/assess_rerun_changes.py`)

Artefacts, all keyed, all deposited:
1. **Three input ledgers, none subsuming another** (round-4 finding 1): (a) the **raw splice
   ledger** of Phase B (every day whose stored Q changed: old / new value, code, source); (b)
   the **inserted-row inventory** of the join (every climate row added, with its Daymet values
   and the Daymet-missing dates among them); (c) the **preprocessed daily delta ledger** (B′):
   for every affected gage, `preprocess_daily_data` is replayed on the S0 frame and on the S2
   frame (and on the S1 / S1′ frames for the interaction cohort); every changed PREPROCESSED
   daily value of Q, PPT, SWE, TMIN, TMAX is written with old / new value, old / new
   missingness and the cause (splice day, inserted climate row, secondary interpolation change
   such as the Dec 31 neighbour). (c) cannot replace (a): a restored observation can equal its
   previous interpolate exactly, leaving every preprocessed value and the valid-year set
   unchanged while the RAW seasonal completeness (`io.jl` l. 326–338 counts raw non-NaN days
   per season against the 80 % threshold) flips a season from incomplete to complete — the
   seasonal outputs and `season_excluded_years_*` then move with no daily-value delta at all.
2. **Per-stage state tables** (B′; every gage of the executed cohort, rejected years included):
   per gage, stage and water year the raw NA count, the longest raw NA run, the raw
   missingness mask and the interpolation mask (which days of the preprocessed series are
   interpolates — the only source for item 7's interpolated-Q exposure; preprocessed
   missingness cannot identify them), the four raw seasonal completeness flags, the rejection
   reason, the Q-, PPT-, SWE- and temperature-validity, and the derived qualification metadata.
   Affected membership (E5) and gate 5a are derived from items 1 AND 2 together.
3. **Dependency closure, enumerated per OUTPUT KIND** (the registry; round-4 finding 3 — v5
   stopped at the annual values):
   - **annual values** (numerical): year-local signatures depend on their own year's daily
     values, masks and eligibility; record-dependent families (drought thresholds, `*_all`
     pulses, `elasticity_rolling` AND `elasticity_annual` — normalised by record-wide means,
     differencing adjacent RETAINED years —, `elasticity_static` and its counters,
     `runoff_ratio_high_count`, parameterized BFI, `avg_storage`, the recession seasonality
     bases and `recession_alpha_point_cloud_linear_reservoir`) depend on every year of the
     gage;
   - **annual-key emission** (which gage × year × signature cells exist in the parquet): the
     family gates that run over the whole record BEFORE collection — the recession
     minimum-event gate (`recession.jl` l. 515: fewer than `min_events` fitted events over
     the record → every base metric empty; a 25th event in ONE year makes `log_a_events` and
     its siblings appear for every earlier, unchanged year), the Q–P seasonality floor
     (`qp_seasonality.jl` l. 65: fewer than `min_years` = 10 usable years → nothing emitted;
     plus its usable-result floors), the ≥ 3-row collector floors. Keys of OLD years newly
     emitted (or withdrawn) through such a gate form the class **"gate-emitted old year"**,
     listed per gage with the gate and its before / after count, validated by recomputation —
     they are neither "newly eligible years" (gate 5d) nor numerical deltas;
   - **summary statistics** (16 per base): depend on the gage's whole annual series and its
     year set; the 20-value stats floor and the 60 % / 80 % trend and decade gates use each
     metric's OWN finite-value span, so a gained year can switch them on or off; the **snow
     record-anchored decade gate** (`snow.jl` l. 285–298) and the **frost gate anchored to the
     temperature-valid years** (`climate.jl` l. 185–209, 509) force-skip the TREND statistics
     only — an anchor change permits trend-statistic changes and NEVER a shared-year annual
     value, mean, median or Pettitt field of a snow or frost base;
   - **per-gage scalars**: `season_excluded_years_*` (every year in the loop, from the RAW
     seasonal flags), `num_water_years`, `start_` / `end_water_year`,
     `ice_affected_days_total`, `recession_alpha_point_cloud_linear_reservoir`,
     `elasticity_static`, the drought thresholds;
   - **QA flags**, incl. `flagged_for_high_na` (from the final fields over the PINNED schema,
     Phase C).
   "Permitted to change" = the closure of the three ledgers and the state-table changes under
   these rules, per output kind. It is a permission, not a validation.
4. **Observed changes**: keyed cell diffs of the summary CSV and of the annual parquet between
   each stage pair, NaN == NaN, null == null, signed zeros equal, computed by gage partition in
   bounded memory.
5. **Gates**: (a) **0 unexplained changes** — every changed cell (value, key emission, summary,
   scalar, flag) is permitted by item 3 and traced to a ledger or state-table entry; an
   accounting gate, not a numerical one; (b) for the linear totals (`prcp_annual` and the four
   seasonal precipitation totals) on shared eligible years, the annual delta equals the scoped
   sum of the ledger's daily PPT deltas (abs 1e-9 + rel 1e-9); for `gdd_cum_*` the delta
   equals Σ over the changed dates of `max(0, Tmean_new − base) − max(0, Tmean_old − base)`,
   `Tmean = (TMAX + TMIN) / 2` reconstructed from both components on each side (`climate.jl`
   l. 476 sums THRESHOLDED daily degree-days: a day moved from −2 °C to +2 °C changes the
   annual GDD by 2, not by the 4 of the mean or the 8 of the summed component deltas); (c)
   **every** other changed or gained annual value — the totals of NEWLY eligible years
   included — equals an absolute recomputation of the library on the reconstructed daily
   series (E7's port parity on the affected gages is that recomputation in a second language;
   extended with long-window-only and fixture cases); (d) newly eligible years are counted as
   gains and their values validated under (c), never as numerical deltas against an invented
   old value; gate-emitted old years (item 3) likewise; (e) the stage identities of §10.4 on
   the EXECUTED cohorts, inferred cells labelled; (f) the freeze (gage sets, `total_possible`
   and the full-product header identical to v2's) and the valid-year identity = oracle; (g)
   changed summary statistics, per-gage scalars and QA flags are gated numerically by E7's
   summary-CSV comparison in the second language, with gate-state agreement — the summary
   layer is validated, not only the annual layer.
6. **Magnitude census per column**: cells changed / gained / lost, max |Δ| with its gage, p50 /
   p99 |Δ| (absolute, and relative where the old value is non-zero; absolute only at zero),
   identity R² (1 − Σ(new − old)² / Σ(old − mean)², reported as max |Δ| only for constant
   columns); for the **695 `_pval` columns (5 suffix types)** the number of gages crossing
   0.05 and 0.01 in either direction; for `_pettitt_cp_year` the number moved and the shift
   distribution; for `_dowy` bases the days moved; for the 12 QA flags and
   `flagged_for_high_na` the flips; for `_senn_slp` / `_linear_slp` the sign changes.
7. **Per-gage and paper-facing views**: the gages with the most changed cells and the largest
   standardized |Δ| (per column, relative to its cross-gage sd), each with its cause; same-year
   value corrections separated from year-composition changes; regional (HUC2 / province) and
   seasonal distributions of the climate deltas; per affected year, the interpolated-Q
   exposure from the interpolation mask of item 2 (interpolated Q days and the Daymet P on
   them — exposure to the interpolation residual, not its size); product-wide median and mean of
   every signature's `_mean` before vs after — Q-only population means may shift slightly
   through the splice and are reported, not assumed unchanged.
8. **Outputs**: `change_census_<pair>_<window>.md` + CSVs + a comparison dashboard per stage
   pair in the run folders; STATUS / CHANGELOG / README quote them; the daily delta ledger,
   the stage manifests, the per-stage year-set tables and the gate-state census join the R3
   deposit inventory (Phase H).

### 10.6 Decisions (defaults first)

12. **Fix placement**: a runner-level shared helper behind the `STREAMFLOW_CLIMATE_FRAME`
    switch, production join only; the library untouched. Alternative: pass the climate frame
    into `preprocess_daily_data` as its own argument (a cleaner contract for a later release).
13. **Restriction of the inserted climate rows**: to the record span AND the original water-
    year set (proven by the fixtures of §10.2). Alternative: span only (rejected: interior
    years, Q-only scalars, denominators).
14. **Stages**: S1 on the splice-affected gages + control sample; S1′ on one composed batch;
    all on one tree. Alternative: full S1 products (5–7 h, no extra rigour).
15. **Daymet Dec 31 of leap years**: unchanged (absent in Daymet, interpolated on the grid,
    documented); the ledger carries its secondary changes. Alternative: fill it from the
    Dec 30 / Jan 1 mean in the runner (a further change to the climate values, not requested).
16. **Python runner**: validated by the cross-language fixture in both modes and E7 on the
    affected gages; no full Python product. rpkg: the reference re-baseline (F) covers its
    runner change.
17. **The coupled Q–P metrics** use the SAME full-frame P as the climate-only metrics (§10.8).
    Alternative — a second, left-frame P for the coupled families only (not recommended).
18. **Baseline drift** (round-4 finding 5): S0 ≠ v2 STOPS Phase D; the drift is explained and
    resolved to equality before any stage runs. Replacing v2 by the reproduction as the
    comparison reference is a user decision recorded in the release notes, never a default.
    Alternative: continue with the reproduction as the reference after a written explanation
    (v5's wording; rejected — it would let an unexplained tree change ride into the product).
19. **Schedule** (round-4 finding 7): one person, sequential, from 2026-10-26, earliest
    completion 2026-11-12/13 without a rerun — AFTER the 2026-11-09 submission; the submission
    goes out on the version-2 products and counts, and Phase H diffs every manuscript number
    against v3 and relays the changed ones. Alternative: a second person on B ∥ C and A′ plus
    night and weekend machine time → 2026-11-06 at the earliest, zero slack, no rerun.

### 10.7 Documentation when it lands
CHANGELOG "Fixed (MEDIUM)" and the Known Issue closed; STATUS line removed; `docs/SIGNATURES.md`
§8 / §12 / §14 caveats rewritten; `docs/DEVELOPMENT.md` NA-pipeline paragraph (l. 83);
`claude-skill/streamflow-signatures.md` l. 126; relay item 1 of `2026-10-12-coauthor-relay.md`
closed; R1 / R2 READMEs; a Methods sentence for the manuscript. The climate support is stated
as **Q-valid water years within the streamflow record** — the fix does not extend any climate
metric to years the streamflow record rejects (the 2026-10-12 decision 13 stands; the user
confirmed this intent on 2026-10-26).

### 10.8 Does the fix change the Q–P metrics, and is that sound? (user question 2026-10-26)

It changes them — directly on the affected gage-years, and indirectly on other years of the
same gage through the record-normalised families (`elasticity_annual`, `elasticity_static`,
the drought thresholds, the `*_all` pulses) and the family emission gates (§10.5 item 3) — and
the situation the question worries about ("a whole year's P against a partial year's Q")
cannot arise in the production path:
- A Q–P signature is only ever computed for a **Q-valid** water year (≤ 30 missing days, no
  run > 3), whose missing days are linearly interpolated for Q; the PPT-dependent families
  additionally receive only PPT-valid years (`run_julia_benchmark.jl` l. 348). The Q side is
  always a complete series with at most 30 interpolates. Years with a real portion of Q
  missing are rejected and get no signature of any family.
- On the interpolated days that were ABSENT ROWS the OLD frame interpolated P too; on the
  interpolated days served as Q = NA rows (the 193,542 masked days of §1, for instance) it
  already kept the Daymet P against an interpolated Q. The NEW frame keeps the Daymet P and
  the interpolated Q on both kinds (the standing NA policy, unchanged); the old frame made
  the treatment of a missing day depend on whether NWIS served an NA row or no row — an
  arbitrary distinction the fix removes.
- Family-specific eligibility still applies and is in the closure (§10.5): runoff ratios use
  paired finite days, precipitation floors and the raw-Q seasonal completeness flags;
  elasticity filters low-P years, differences adjacent RETAINED years and normalises by
  record-wide means, so distant years can move; Q–P seasonality uses rolling precipitation
  thresholds and minimum monthly coverage; `avg_storage` needs ≥ 10 unique Q values;
  `swe_max_to_ppt` needs SWE and PPT eligibility and sits under the snow record-anchored gate.
- No direction is asserted or gated: the restored P is smaller than the old interpolate on
  501 of the 919 deviating gage-years and larger on 418; `avg_storage` interpolates the
  cumulative balance at mean discharge, so a positive annual ΔP can lower it. Every coupled
  change is validated by recomputation (gate 5c), not by sign. Nor is the residual claimed
  small, or measured: an absent row does not say whether a storm fell that day, and Q is
  interpolated across it either way, so interpolated Q over a short flashy event can miss
  most of its runoff. That residual is the existing NA policy's, shared with every other
  ≤ 3-day gap, and its size is unknown; what §10.5 item 7 reports per affected year
  (interpolated Q days and the Daymet P on them) quantifies the EXPOSURE to it, not the
  residual itself.
- Keeping a separate, left-frame P for the coupled families (alternative of decision 17)
  would make the P behind a `runoff_ratio` row differ from the `prcp_*` value of the same row
  for no hydrological reason. What the single frame buys is narrow and stated as such: a
  ratio identity such as `Qfal / prcp_fal = fall_runoff_ratio` holds only where both sides
  share the annual or seasonal eligibility and denominator; a ratio of two summary MEANS is
  not a mean runoff ratio and is never gated. The production Q side stays a complete series
  with ≤ 30 interpolates per valid year; decision 17 stands.


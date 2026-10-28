# Fix the USGS qualifier mask — clean up both products for their existing gages and QA/QC them against the version-2 runs

> **EXECUTED 2026-10-26 … 2026-10-28.** The version-3 products are built and verified (Phases A–H through local
> staging); the consolidated fixes / final code / outputs / results are in `2026-10-21-qualifier-rerun-results.md`,
> the dated log with every defect and review disposition is §11, the reviews are in `2026-10-21-qualifier-reviews/`.
> Open: the user's commit / push / publish, the deposit, the climate-frame default (STATUS).

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


## 11. Progress log (implementation started 2026-10-26; decisions 12–19 taken as their defaults when the user said to proceed)

| When | Phase | Done / found | Deviation from the plan text |
|---|---|---|---|
| 2026-10-26 | A | Rule in `R/usgs_qualifiers.R`, applied by `generate_streamflow_dt`; `R/tests/test_qualifier_codes.R` 37 checks (integration check runs where lubridate is installed); docs updated. | The R3 input-dictionary `flag` row waits for the restaging (Phase H) — the staged file describes the February parquet until then. |
| 2026-10-26 | A′ | `climate_frame.{jl,py,R}` + `STREAMFLOW_CLIMATE_FRAME`; self-tests Julia 40 / Python 37 / R 42; batch 1 of #1 in `left` bitwise = v2 (`check_additivity.jl` PASS, 824 gages); 40-gage subset in `outer` moves 599 columns, all climate-dependent. Reviews: Sonnet GO WITH CHANGES, Codex NO-GO → all findings applied (duplicate climate dates, date-type assertions, gate state never exempt, window denominator tested, rpkg legacy branch, R `year`/`doy`). | Default mode stays `left` until the products land (then `outer`, Phase H). |
| 2026-10-26 | C | `run_batched_julia.py` (`--gage-list`, fingerprint, `--require-clean-tree`, `--qa-schema`, stage provenance; self-test 28/28); `STREAMFLOW_QA_SCHEMA` in the three runners (subset header = the 2,277-column product header, `flagged_for_high_na` = v2); `validate_production_run.py --expect-same-gages`; `audit_qualification.jl` cohort contract; `run_rpkg_acceptance_gates.sh` parameterized; registry `docs/qualifier_rerun_registry.csv`. | Whole-table readers of `validate_annual_values.py` / `check_annual_parquet_equality.py` kept (they validated 34.7 M rows on this machine 2026-10-13); the staging tool's qualifier-release mode moves to Phase H. |
| 2026-10-26 | B | `splice_nwis_qualifiers.py list` reproduces §1 exactly (193,542 / 987; 15,257 / 263 provisional; 1,185 masked gages; 12,103 no-code rows). `fetch`: 1,184 / 1,185 on the first pass (0158397967 HTTP 503, retried on the second pass). `qualifier_year_oracle.py replay` reproduces every one of the 322,563 census year statuses and every window decision. `splice` / `control` / `audit` tested on 5 gages: 1,391 restored (1,365 `A R`, 26 `A <`), B2 identity 99.965 %, metadata rule reproduces the February values of the unchanged gages. | (1) The February compile window ends 2025-10-01 inclusive (config.R `DEFAULT_END_DATE`), the first fetch ended 2025-09-30: 49 masked rows (4 splice targets, all WY 2026, outside both products) had no source day — the whole universe is re-fetched with `--end 2025-10-01` for a uniform deposited cache. (2) Three sites (06294000, 01467087, 021473428) hold two METHOD BLOCKS inside one 00060/00003 series; the unqualified block (dataRetrieval's `Flow`) is the one the compile stored and the only one read; B1 reports them. |
| 2026-10-26 | B′ / E8 | `climate_frame_inventory.py`: the realized inserted-row inventory equals the join forecast on the product gage-years — 1,174 / 857 / 3,183 (#1) and 937 / 721 / 2,195 (#2; the forecast's 938th gage-year, 02068500 WY 1985, has only the Daymet-missing 1984-12-31 absent, which cannot be inserted). `preprocess_state_ledger.py` (replay + diff → daily delta ledger with causes, state-table changes) and `assess_rerun_changes.py` (classes, keyed annual diff, explanations, magnitude census, gates incl. 5b) pass on the 40-gage subset (1,810 changed annual cells, 0 unexplained; 44 `prcp_annual` deltas = summed daily deltas). `build_stage_cohorts.py` composes S1 / S1′ from the realized ledgers (on the realized join alone it finds J = 857, the two eligibility-gain gages 0214269560 and 02484000 and six Dec 31-neighbour gages of #1). `check_summary_equality.py` is the E7 summary gate (cell for cell, gate state always, tie-bearing rank statistics reported). | Over the whole record the inventory holds 3.3 M inserted rows on 3,419 gages (mostly in rejected or partial years); the census and the cohorts use the rows in Q-valid years only. E5 and E8 are ONE tool (`assess_rerun_changes.py` run on the product vs v2 and on every stage pair); `compare_qualifier_rerun.py` is not a separate script. |
| 2026-10-26 | Codex tooling review | Blind review of the B / C / B′ / E8 tooling (`docs/plans/2026-10-21-qualifier-reviews/codex_tooling_review_round1.md`): NO-GO, 5 BLOCKER + 9 MAJOR + 1 MINOR, every one verified and applied: the audit's B4 round-trip and B6 "every row must change" false failures; the keyed source → ledger → output identity (B5 now checks every ledger key's base NaN and splice flag and the patched value, B4 re-derives every value from the cache and the metadata area, B8 recomputes the base targets from the parquet, B1 reopens every cache); the splice enforces the status-letter rule, the unqualified METHOD BLOCK, the site code, the unit and duplicate dates; the fetch resume is bound to the request URL (window included); the oracle treats ±Inf as missing, binds the rule constants to the config (sha256 recorded), consumes the whole ledger, validates every target against the base parquet and FAILS on a lost year; the census reports every gate of §10.5 as PASS / FAIL / NOT RUN with a verdict that is PASS only when every required gate passed (INCOMPLETE otherwise), backs each explanation with evidence (the delta ledger's changed variables, the state tables' validity flips, the recession / Q–P gate counts crossing their thresholds), gates 5b over EVERY shared eligible precipitation and GDD cell (unchanged cells included), takes the E7 verdict JSONs for 5c / 5g, an expected-identical list for 5e, the oracle's exact sets for 5d / 5f, and adds the population-mean, standardized, regional and interpolated-Q-exposure views; the delta ledger verifies every cause against the ledger / inventory values and recomputes every interpolate from its own stage's raw series, exports per-variable masks, three value layers (raw / clean / eligible) and the window qualification per gage, and fails on unmatched keys; the batch fingerprint covers the GAGES-II directory, the HYDAT file, the resolved QA schema, the climate frame and the dirty-tree content, resumes only with recorded output hashes on a clean tree and re-takes the snapshot before every batch; the cohort audit compares emission with qualified ∩ selected and requires the census; the summary gate never waives a gate-state difference and exempts tie-bearing rank statistics only within the annual tolerance. The reviewer's injection (+123 into an annual value and a summary mean of an inserted-row year) now FAILS the census through 5c / 5g (E7) — permissions are no longer validations. | The 40-gage and 5-gage tests pass end to end on the reworked tools; the first full chain run on the old tooling was stopped and rerun. |
| 2026-10-26 | Sonnet tooling review (round 2) | Blind review of the reworked tooling (`docs/plans/2026-10-21-qualifier-reviews/sonnet_tooling_review_round2.md`): GO WITH CHANGES, 4 MAJOR + 5 MINOR + 3 NOTE; the reviewer independently re-derived all 193,542 ledger values from the cache (0 mismatches), confirmed the unqualified method block is the one the compile stored (06294000: 0 differences vs February, the adjusted block 1,828), and that the compile rule reproduces the February metadata of all 987 splice gages. Applied: the B8 loop that would have taken ~2 h (a set rebuilt per iteration); the E7 verdict JSONs now carry the sha256 of their inputs and the census binds 5c / 5g to the assessed files and accepts waivers only by name (`--accept-waiver`); a per-gage B2 gate (the ±3-day neighbours of every restored day equal February at 1e-9, so a wrong basin area cannot hide in the aggregate); the state-table attribution reproduces every changed field from the daily deltas (raw-NA counts, masks, seasonal flags only in a season with a changed Q day, validity from the preprocessor's rule on the new mask); B1 pins the request window; 5b requires comparisons; the pinned schema adds signature columns only and errors on an absent metadata column; `merge` checks the schema the batches ran with; B6 asserts the rule reproduces the February metadata on the base series (batch-runner self-test 28 / 28 on a quiet machine; one timing check wobbled once while three Julia runs shared the CPU). Noted, not changed: 97 restored days carry status `P` (96) or `A` (1) today under a stale February flag (the ledger records `qualifiers_today`); 554 restored values are negative (the config keeps negative flow). | The decision-1b sensitivity run (`bounds_sensitivity.py`) was added as a Phase B deliverable the tooling had not covered. |
| 2026-10-26 | B, decision 1b | `bounds_sensitivity.py` on the 163 gages that hold `A <` / `A >` days (11,662 days; the patched input, Julia, outer frame, pinned schema, window 1993–2025): 156 gages qualify with the bounds as values, 155 with them as NA — the one that qualifies only with the bounds (12202310) is already on the §6 kept-out list; `num_water_years` differs on 79 of the 155 gages that qualify under both policies (the first report said 87: it counted NaN ≠ NaN and the newly qualifying gage); on those gages 22,325 annual cells change and 21,572 are gained (15,345 of the changes in record-dependent bases, 1,503 Q-only cells on restored-day years, 162 recession keys newly emitted). `~/HISSS_data/qualifier-rerun-26oct2026/sensitivity_1b/` (`bounds_sensitivity.md`, `qualification_sensitivity.csv`, `sensitivity_census/`). | A diagnostic of the policy, not a gate: the census's oracle and E7 gates are deliberately NOT RUN for it. |
| 2026-10-26 | B, full chain | `~/HISSS_data/qualifier-rerun-26oct2026/input/`: fetch 1,185 / 1,185 (window 1979-10-01…2025-10-01, 1.2 GB of verbatim responses); splice: 193,542 of 193,542 targets restored (A R 181,461, A < 11,444, A [4] 419, A > 218), 0 unavailable, 987 gages, patched parquet `combined_streamflow_data_09feb2026_qualifier_patched.parquet` (sha256 2d44c1af…; base 65fb5563…), metadata rows of the 987 gages updated in place (88 with new compile-year totals), `correction_ledger.csv`, `provenance.json`; control: 13,786,734 stored values compared with today, 99.770 % identical at 1e-9, 31,676 revised by the agency since February (66 gages above 1 %, 03358000 at 33 % — reported, never applied), 193,713 days restored by the rule on the October record (193,445 on splice codes; 97 splice-code days carry status `P` or `A` today; 1,991 February-NA days a kept code carries today and 8,129 provisional-remark days accepted today stay OUT by decisions 1e / 4). Audit: B1, B3, B4, B5 (193,542 Q cells changed, all on ledger keys, none off), B6, B8 PASS; B2 passed its aggregate (99.770 %) and the per-gage identity (every gage ≥ 67 % of its values identical; the first zero-tolerance neighbour gate flagged 02423397 and 03146500, whose neighbours the agency revised — re-operationalised as the per-gage share; the rerun at 16:01 reports ALL GATES PASS). Oracle (`oracle/`): product #1 619 gages gain 1,177 years, 899 with restored days; #2 502 / 950 / 757; all compiled gages 676 / 1,370; 0 gages lose qualification, 0 source-unavailable targets; would-qualify outside the products 31 / 29 = the §6 lists with the four codes; config bound (max_raw_na 30, gap 3, residual NA rejects, negatives kept). B7: the Julia census on the patched parquet = the oracle's final sets on 16,028 gage-windows and 322,563 gage-years (0 mismatches), and the Python replay reproduces that census (0 mismatches). | Every count equals the forecast of §1 / §6 / §10.3. |
| 2026-10-26 | Codex tooling review (round 2) | `docs/plans/2026-10-21-qualifier-reviews/codex_tooling_review_round2.md`: Phase B OUTPUTS GO WITH CHANGES — the reviewer independently reconfirmed all 193,542 restorations, all 111,624,189 rows / 909 row groups / schema / order, the 987 metadata rows, the 322,563 statuses and 16,028 window decisions, the oracle's exact final sets and a non-zero-only per-gage conversion identity (minimum 66.984 %) — and the VERIFICATION TOOLING NO-GO for Phase D (2 blockers, 11 majors, 1 minor). Applied the same evening: the summary gate no longer exempts a residual pair's cells but gates each against the statistics recomputed from its own side's annual series (blocker 1); the delta diff validates every layer on both sides — raw against the source ledgers, clean against raw / the recomputed interpolate in a Q-valid year, eligible against clean × validity — and rebuilds every changed state field from the stage's daily arrays (blockers 2–3), with the Python port's rejection wording mapped to Julia's classes (4); E7 verdict JSONs carry their compared-gage list and the census requires it to cover every affected gage (5); infinities are a gate state (6); the per-gage B2 share counts non-zero values only and the control is hash-bound to the metadata and cache (7); a base-target gage without an opened source record fails B8 coverage (8); the fingerprint hashes the GAGES-II directory and HYDAT file as the runner resolves them (env, else config) and the snapshot is re-taken after the last batch (9); a stale annual parquet counts as missing, `--require-annual` applies before resume and `merge` re-checks the recorded output hashes (10); the GDD gate requires equality when no temperature day changed and full coverage otherwise (12); the total_possible and 5e gates fail on absent evidence (13); B4's plausibility uses the February range, calendar-day neighbours and the base series' own ratio distribution (14). Not changed: the emission explanations stay heuristic (the library's gate counters are not exported — decision 12 keeps the library untouched) and are validated by E7's exact key-set gate (11); the interaction residual and inferred-cell labels of §10.4 await the user's decision on keeping S1 / S1′ (13). | Re-tested on the 40-gage and 5-gage sets; the full control + audit rerun with the round-2 gates reports ALL GATES PASS (16:23; non-zero per-gage B2 minimum 66.98 % at 03358000, every base-target gage's source record opened); batch-runner self-test 28 / 28. |
| 2026-10-26 | B CHECKPOINT | DECISION (user): approved; S1 / S1′ KEPT. Baseline tag `qualifier-rerun-baseline-2026-10-26`; Phase D runs from a clean `git worktree` of the tag so the main tree stays editable (the runner's per-batch snapshot would otherwise see every doc edit as a dirty tree). Preview of product #1 (patched input, outer, pinned schema; stage `S2preview`, dirty tree) + internal review summary for the user's data review — not a product, not documented beyond this row. | Measured batch speed 835 gages ≈ 2 min (7.5 gages/s): S0 for both windows is ~1 h, not the 5–7 h estimated in §0. |
| 2026-10-26 | Codex final scope review | `docs/plans/2026-10-21-qualifier-reviews/codex_final_scope_review.md` — the user's question "did we change nothing but the two facets, and does every output change follow from them physically?": the INPUT patch is within scope (the reviewer scanned both parquets in full — 111,624,189 rows, 909 row groups, 193,542 Q cells changed from NA and nothing else, 4,470 negative zeros preserved; every restored value re-derived from the cache; metadata changes confined to `num_years` (88), `num_days` (88), `end_year` (2); Daymet md5, the 16 reused batch files, config, GAGES-II / HYDAT, both selection manifests verified); the CODE change is within scope for the production path (library untouched; `left` the former join; batch 1 byte-identical at both layers; the smoke's 599 changed columns all climate-dependent, 1,810 annual cells, totals = summed daily deltas); the UNCONDITIONAL claim is not established until the products run (S2 absent). Verdict: NO-GO for the claim as written today, GO WITH CHANGES to proceed with S2 and the full-product validation. Reporting corrections applied: (1) 96 restored values at 12344000 (2024-05-15…08-18) are PROVISIONAL (`P`) today under their February `A R` flag — the audit now states that the restored days are selected by their February approved flags and valued from the October record under its status today, and names them; (2) the sensitivity report's "87 shared gages differ in `num_water_years`" counted NaN ≠ NaN and the newly qualifying gage — among the 155 gages qualifying under both policies it is 79; (3) the census and the review tables now separate NUMERICAL differences (within the annual tolerance, e.g. a drought deficit moving in its 15th digit) from SUBSTANTIVE ones. | The cross-language threshold / tie residuals are named, never counted as hydrological effects. |
| 2026-10-26 | E7 dry run | The Python runner on the 40-gage subset in `outer` mode with the pinned schema (1.64 gages/s → ~16 min for 1,600 gages) vs the Julia run: annual values identical on 162,373 keys except ONE last-bit threshold tie (`D40_day`, 02160390 WY 2001, 167 vs 168 — the documented `D*_day` residual class); the summary gate passes with that residual's propagation (10 cells) and 66 tie-sensitive rank statistics reported, and ONE named column waiver: `ice_affected_days_total`, which the Python port counts from the `flag` column (23 at 01011000) while canonical Julia emits 0 for every gage (Known Issue; rpkg deliberately matches Julia). | The E7 waivers are therefore: the `D*_day` / `D25_to_D75` / `TQmean` / `dur_low_pulses_*` threshold-tie class on the annual gate (`--allow-value-diff`, counts quoted), the propagated-residual and tie-bearing rank cells on the summary gate (reported by the tool), and `ice_affected_days_total` by name. |
| 2026-10-27 | D, S0 = v2 gate | Both windows from the clean worktree of the tag (`S0_1993_2025/`, `S0_1980_2025/`): summary CSVs byte-identical to v2, annual gate PASS (26,782,608 / 33,651,562 finite pairs, max diff 0), summary gate PASS (2,276 columns bitwise), additivity PASS — decision 18 satisfied (2026-10-26 21:45 / 2026-10-27 04:36). |
| 2026-10-27 | D, S2 products | #1 v3 `processedOuts_1993_2025_climate_v3_26oct2026/` (6,678 × 2,277; annual rows 27,766,578 = v2 + 164,417; merged 08:47; E1 6/6 PASS) and #2 v3 `processedOuts_1980_2025_climate_v3_26oct2026/` (6,250 × 2,277; annual rows 34,813,881 = v2 + 132,550; merged 13:36; E1 6/6 PASS), both `outer`, patched input, pinned v2 schema, stage S2 with the selection manifests. Overnight batches ran 15–25× slower than daytime (full swap / idle power management; `caffeinate -i` attached). Script 1 was stopped at 12:06 during batch 4 of #2 for the tool fixes below and resumed as script 1b at 13:24 (batches 1–3 reused by fingerprint; the interrupted b4 folder set aside in `logs_interrupted/`). |
| 2026-10-27 | B′ / E8 on #1 (full scale) | Three verification defects, none in the products or the product code (the worktree is unchanged; the fixes are post-tag in the main tree): (1) `bprime/inserted_rows_all.csv` had been concatenated through `pd.read_csv` with the default float parser, which is not round-trip — 1,277,602 of 3,352,085 rows 1–2 ulp off → 53,924 "new raw value is not the inventory's" daily deltas; rebuilt by text concatenation of the exact per-batch files (0 mismatches against them; the old file kept beside it). (2) `preprocess_state_ledger.py` diff expected the interpolate on a climate NA day of a year the variable is INVALID in (> 30 NAs); the preprocessor leaves those days NaN and keeps the finite days raw (0214269560 WY1997, 02484000 WY2017; 128 cells) — check 2 now uses the variable's own validity. (3) `assess_rerun_changes.py`: 5f and the total_possible freeze failed on the 5,069 non-replayed gages (B′ covers the affected manifest only) → Qann annual keys for them, the freeze scoped to the replayed gages with every other gage required to be class U; 116 keys emitted / withdrawn on OLD years of record-dependent bases (`BFI_*_param` 92, `elasticity_*` 24; `collector_min_rows` gate) now explained as record-dependent (validated by 5c) or restored-day years. After the fixes: delta ledger COMPLETE — 2,702,639 changed daily values (eligibility 2,012,871, inserted climate row 544,884, splice day 144,857, secondary interpolation 27), 3,903 state rows changed; census PASS on every gate that runs before E7 (0 unexplained of 387,055 non-equal annual cells; class U identical at both layers; 5b ×2, 5d, 5f PASS; 5c / 5g / 5e NOT RUN). Internal review artefacts for the user (not for the repository): `review_1993_2025/review_summary.md` + CSVs and `v3_1993_2025_vs_julia_dashboard.html`. |
| 2026-10-27 | B′ / E8 on #2 | With the corrected inventory and tools (script 1b, 13:37–14:01): replay of the 1,356 affected gages, delta ledger COMPLETE — 2,217,736 changed daily values (eligibility 1,620,831, inserted climate row 476,952, splice day 119,935, secondary interpolation 18), census PASS on every gate before E7 (0 unexplained of 360,065 non-equal annual cells; annual rows 34,681,331 → 34,813,881: gained 132,551, lost 1, changed 227,268; 4,894 class U gages identical at both layers; 5b ×2, 5d, 5f PASS). Review artefacts `review_1980_2025/` + `v3_1980_2025_vs_julia_dashboard.html` (not for the repository). Fractions touched (valid v3 gage-years and their days): #1 — restored days 142,693 = 0.19 % of daily Q in 1,677 gage-years (0.83 %; 1,177 gained + 500 existing), inserted climate days 3,241 = 0.0044 % in 1,199 gage-years (0.59 %), either 1.41 %; #2 — 116,705 = 0.13 % in 1,353 gage-years (0.53 %), 2,221 = 0.0024 % in 952 (0.37 %), either 0.90 %. Script 2 launched 14:05. |
| 2026-10-27 | D, S1 / S1′ + E2–E4 (script 2, 14:05–14:47; 2b 14:50–15:30) | Both attribution stages from the clean worktree on both windows: S1 (spliced input, `left`; 1,099 / 957 gages = splice-affected + 200 control) and S1′ (February input, `outer`; 375 / 348 gages = interaction cohort 178 / 149 + the join-only sample), replayed by B′. Stage censuses: the 200 control gages identical at both layers (5e PASS); the join-only sample identical on 197 / 199 gages — 3 / 1 of the 200 sampled lie inside the interaction cohort and are excluded from the identity, not failed; 5d / 5f PASS; 5a NOT RUN by design (no stage delta ledger; 0 unexplained so far). E2 PASS on both products with `--floor 20` (without it 15,458 / 11,697 RULE-4 rows, the floor v2 was validated with); E3 clean on all 16 batch logs; E4 PASS (frozen cohort emitted exactly, 6,678 / 6,250; would-qualify kept out 31 / 29). Script-2 defects, tooling only: the E7 climate subset read `gage_id` / `date` instead of Daymet's `site_id` / `Date` — fixed with an assertion that the subset carries every affected gage, because the Python runner silently runs WITHOUT climate when the file is absent; `duckdb` is absent from the repo venv (summary gate and attribution now run under the daymet venv); E2 had been launched without the floor. |
| 2026-10-27 | E7 both products | The Python port on the 1,609 / 1,356 affected gages (subset of the S2 inputs incl. climate, `outer`, pinned schema). Annual: key sets identical (6,725,235 / 7,778,988 shared keys, 0 one-sided); 117 / 188 of 6,478,680 / 7,501,467 finite pairs over 1e-9, every one in the documented cross-language classes — threshold ties on day-of-year and pulse counts (`D1_day`…`D99_day`, `D25_to_D75`, `TQmean`, `dur_low_pulses_all`, `n_low_pulses_all`; max 4.0 / 5.6) and float noise (`elasticity_annual`, `qp_slope_sd`); gate PASS with the 19 signatures NAMED (`--allow-diff-signature`; 17 of them carry a residual in #1, all 19 in #2) and the count bounded at the observed value. Summary: 0 gated mismatches of 3,338,102 / 2,858,935 numeric cells; 103 / 103 finite-value cells waived by column (`ice_affected_days_total`, and the FDC90th / FDCmid rank + Pettitt statistics of zero-heavy series: 15896000 all-zero in Julia vs 4.8e-29 in Python; 02304500 / 02313230 / 401733105392404 Pettitt change-point tie flips) and, separately, 8 / 7 gate-state (NaN vs finite) cells of those columns accepted only because the statistic recomputed from each side's own annual series reproduces that side (`--allow-gate-state-column`, the one named exception to "gate state is never waived", new 2026-10-27); tie-sensitive rank statistics verified by recomputation 1,510 / 2,503, propagated annual residuals verified 472 / 775, signed-zero ties reported 15 / 18. |
| 2026-10-27/28 | Attribution + final censuses (script 2c, 18:58–22:01) | Identities on executed gages PASS both windows — S1 = S0 on the 200 control gages (455,200 summary cells; 833,690 / 1,101,466 annual rows) and S1′ = S2 on the 197 / 199 join-only gages (448,372 / 452,924 cells; 813,718 / 1,136,764 rows), every cohort gage executed — after `stage_attribution.py` learnt NaN-aware string equality (pandas 3 keeps `RHBN` / `REGULATED` / `CLASS` NaN ≠ NaN: 354 / 394 false unequal cells with identical values). INFERRED accounting: the full-cohort S1 takes 5,579 / 5,293 gages from S0. Interaction residual (S2 − S1) − (S1′ − S0) on the 178 / 149 interaction gages REPORTED, not gated: 2,650 of 684,034 / 2,779 of 800,735 finite annual cells outside atol 1e-9 + rtol 1e-6 (656,930 / 773,641 of them on gages executed in all four stages, where every outside-tolerance cell lies), max 21.2 / 20.2, almost all `elasticity_annual` (2,329 / 2,533) and `elasticity_rolling` (292 / 227) — elasticity is built from year-over-year changes of Q and P normalised by record-wide means under a 0.1 threshold, so a Q change (splice) and a P change (frame) on the same gage multiply rather than add; every outside-tolerance cell sits on an overlap or restored-with-inserted gage (gained / lost keys are excluded from the residual and listed separately); the remainder ≤ 0.003 except `avg_storage` 1.3 on 3 gages. Realised interaction cohort 178 / 149 = overlap 147 / 122 + restored-with-inserted 24 / 22 + eligibility 2 / 1 + Dec 31 5 / 4 (the design's 165 / 135 and 18 / 13 were forecasts); 171 / 144 of them are executed in all four stages — the 7 / 5 eligibility-gain and Dec 31 gages are not in the S1 cohort (S1 = splice-affected + control), have no restored day inside the window (checked against the ledger; a coverage gate, Codex 2026-10-28) and carry S1 inferred = S0, labelled `s1_inferred` in the residual tables. Keys absent from a stage 34,487 / 29,318 (gained years `s1,s2` 34,361 / 29,255; climate-eligibility gains `s1p,s2` 126 / 63); keys present in all four stages but NaN in all four 26,344 / 30,313; NaN-state differences 69 / 57 — key presence and value finiteness are separate masks since the Codex review (the earlier 88,004 / 86,782 'absent' count had conflated them). The two STAGE censuses (S1 vs S0, S1′ vs S2) end INCOMPLETE by design — no stage delta ledger, no E7 on a stage — and gate only class U, 5d, 5e (the control / join-only identities) and 5f; the PASS verdict is the v3-vs-v2 census's. Final v3-vs-v2 census with the E7 verdicts bound by sha256 to the assessed files, every waiver accepted by exact name (19 signatures, `allow_value_diff`, the ice column, the FDC regex) and 5e on the 5,069 / 4,894 gages outside the affected manifests (`manifests/expected_identical_<W>.txt` = selection − affected): every gate PASS, **verdict PASS on both windows** (`census_<WU>/change_census_v3_vs_v2_<W>.md`; the #1 census took 51 min on a swapping machine). Phase D complete. |
| 2026-10-28 | E6 | Signature explorers of both version-3 products built into the run folders (`signature_explorer_<WU>_60pct_climate_v3.html` + `_annual/` sidecars; 100.4 MB, 6,678 gages, 2,245 mapped variables for #1) and the headless basemap check PASS on both explorers (60 / 60 Esri tiles) and both v3-vs-v2 dashboards (24 / 24; attribution present, no failed request). `flagged_for_high_na` comes from the runs (pinned schema). |
| 2026-10-28 | H staging (local) | `docs/benchmarks/stage_hydroshare_qualifier_v3.py` (new; `signatures` / `inputs` / `check`) stages `~/Downloads/Signatures_v3_28oct2026/` (public name "version 2"): R1 / R2 from the v3 run folders — deposit names of 2026, the explorer with its embedded sidecar name rewritten, the frozen selection manifest, the would-qualify list (31 / 29), the validation record (E1, E2, E4, the change census + gates + classes, the E7 port check), a README whose counts and the three new conventions bullets (corrected input, climate frame, frozen gage set) are read from the files, and run notes with the input hashes and the batch runner's recorded commands; R3 = the immutable February compile beside the corrected table (sha256 bound to the control provenance), the ledger, the raw NWIS responses (tar.gz, 1,185 sites), the control report (+ years, provenance), the acceptance oracle, the audit, the selection manifests, the executable reconstruction recipe (`reconstruct_qualifier_patched_streamflow.py`, row-group-wise; run on the staged files: 111,624,189 rows compared, 193,542 ledger days applied, 0 row groups differing — after its ledger read was switched to `float_precision="round_trip"`, pandas' default parser having put 457 row groups 1–2 ulp off), the metadata updated in place with the February copy and the change list, README / abstract / dictionary extended. `check`: every required file present, PAR1 footers, the explorer's sidecar name — 21 / 21, 21 / 21, 26 / 26 PASS, sha256 manifests written. Nothing uploaded. |
| 2026-10-28 | G, round 1 (Sonnet 5.5 ×2, blind) | Evidence + tooling review: GO WITH CHANGES, no BLOCKER / MAJOR (`docs/plans/2026-10-21-qualifier-reviews/sonnet_phase_DEH_review_evidence_tooling.md`): the 387,055 non-equal cells, the 5,069 / 4,894 untouched gages, the 1,177 / 950 restored years and the 117 / 188 E7 residuals reproduced independently; the post-tag patches found compensated by the class-U identity, 5e and the E7 key-set identity. Applied: the class-U identity S2 = S0 that `stage_attribution.py` documented but never ran now runs on `--class-u` (selection − affected): PASS both windows (5,069 gages: 11,537,044 summary cells, 21,041,343 annual rows; 4,894: 11,138,744 / 27,034,893; verdict INCOMPLETE when the list is absent); the stage censuses' INCOMPLETE-by-design scope, the elasticity residual's real mechanism (record-normalised year-over-year form; gained keys excluded) and the realised cohort composition recorded above; the "validated by 5c" labels of withdrawn / emitted keys now name the E7 key-set identity. Staging + docs review: GO WITH CHANGES, no BLOCKER (`sonnet_phase_DEH_review_staging_docs.md`): every count of the staged READMEs, run notes, abstracts and repo docs reproduced (incl. `flagged_for_high_na` 471 / 345); eight MAJOR documentation items applied the same day — the R2 README's unshipped run log (the eight batch logs now ship concatenated), a `?` fallback in the run notes (the E7 counts are now read from the gate text and the JSON, failing when absent), the 2026 dictionary's stale `flagged_for_high_na` row, the explorer titles (rebuilt under the public name; the `v3_vs_v2` validation labels explained in the README), the port-check wording (run notes, `python/README.md`), the R3 README's NA-row count (394,210 corrected / 587,752 compile, read from the tables), the dictionary rows of the new tables (wide `metadata_changes`, 88 of 987 rows change, the conversion factor 2.44657 / area), and the stale 29 / 26 would-qualify count (31 / 29 under the four codes) in STATUS and the manuscript log; `check` now verifies that every file a README or run-notes names is in the folder. |
| 2026-10-28 | G, round 1 (Codex, blind, read-only) | GO WITH CHANGES, no BLOCKER (`codex_phase_DEH_review.md`): independently reproduced the annual rows, the 1,177 / 950 restored and 0 lost years against the oracle (0 of 6,678 / 6,250 gages off), the untouched gages (21,041,343 / 27,034,893 annual cells, 0 differences; summary rows byte-equal), the 117 / 188 E7 residuals by family, all 8 / 7 gate-state exceptions by recomputation, the S0 = v2 hashes, the clean-tree fingerprints of all 16 batches, the recipe on four row groups; explained the single withdrawn key (05527500 `elasticity_annual` WY 2018: ΔP 23.24 → 0.077 mm, below the 0.1 threshold). Four MAJOR, applied the same day: (1) the interaction cohort is not fully executed in S1 (the 7 / 5 eligibility-gain and Dec 31 gages) and `stage_attribution.py` silently intersected → a coverage gate (admissible only without a restored day in the window) and explicit `s1_inferred` = S0 labels; (2) absent keys conflated with present-but-NaN values → separate presence / finiteness masks (numbers above); (3) staged dictionaries: the stale `flagged_for_high_na` row, `ice_affected_days_total` described as a real count though structurally 0, the R3 `gage_id` row advising 8-digit padding (66 compiled gages unmatched that way; strip on both sides) → rows rewritten with counts read from the tables; (4) the R3 manifest's recipe checksum was stale → manifests regenerated after the final edits and verified. MINOR applied: the E7 counts in the run notes; the A-Q label (899 / 757 = a restored day inside the window; 891 / 749 in valid years); R3 NA-row counts 394,210 corrected / 587,752 compile; the one compiled gage without forcing (05KH009) is in neither product; abstract 8,013 of 8,014 basins; the 96 restored values that are provisional today (12344000, 15 May – 18 Aug 2024) stated; waived-cell and gate-state counts separated; #2's 18 secondary-interpolation records; E6 tile counts; 17 / 19 named signatures carry a residual in #1. NOTE kept in the public texts: the census and port PASS mean agreement under the named allowances, not unconditional equality. Unverified by the review (recorded): a repeat of the source-cache restoration, every daily attribution mechanism, the complete Daymet checksum, every explorer interaction. |
| 2026-10-28 | F rpkg re-baseline on #1 v3 | `rpkg` reinstalled from the clean worktree of the tag (08:0x); the worktree's `run_rpkg_benchmark.R` on the 8 S2 input batches of #1 (patched streamflow + metadata, Daymet 2025, `outer`, pinned v2 schema, `STREAMFLOW_HASH_INPUTS=1`): 08:06–09:33, ~11 min per batch, 0 gages errored; merged in the Julia gage order by `merge_rpkg_batches.py` (scratchpad; byte-preserving CSV rows, annual sorted by gage rank) → `rpkg_reference_1993_2025_28oct2026/` (6,678 × 2,277; 27,766,577 annual rows = Julia − 1). Gates: schema PASS; failure scan PASS ×8; annual at 1e-6 WITHOUT waivers FAIL as expected (`gate_annual.txt`): 1 reference-only key (`avg_storage`, signed zero) + 323 of 26,939,257 finite pairs — D*_day 244, `D25_to_D75` 59, `TQmean` 18, `dur_low_pulses_all` / `_year` 1 + 1 (max 40 days, 08158840 `D90_day` 2006); the gate PASS with those 18 signatures named and the counts bounded (`gate_annual_waived.txt`); at 1e-9 also `qp_slope_sd` 146 (≤ 4.4e-11 relative). The 39 climate signatures bit-identical in every annual cell; 53 signatures carry last-bit noise ≤ 1e-9. R² diagnostic 2,227 / 2,244 Perfect, 11 Good, 6 Poor (FDC90th), mean 0.999922 (v2: 2,228 / 10 / 6, 0.999925). `docs/CROSS_LANGUAGE_STATUS.md`, `README.md` port table. |

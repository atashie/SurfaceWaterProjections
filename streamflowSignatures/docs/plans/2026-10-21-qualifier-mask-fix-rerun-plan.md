# Fix the USGS qualifier mask — clean up both products for their existing gages and QA/QC them against the version-2 runs

**Date**: 2026-10-21; **version 3 of the plan, 2026-10-26** after Codex's round-1 review
(NO-GO on v1) and round-2 delta review (GO WITH CHANGES on v2; both in
`2026-10-21-qualifier-reviews/`; dispositions in §9) and the user's decision of 2026-10-23.
**Status**: **DECIDED 2026-10-26** — the user approved the defaults of decisions 1, 2 and
5–11 (3 and 4 as recorded); nothing below has been run yet; implementation starts in a fresh
session, Phase A first.
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

Route (details §3–§5): fix the rule in R; re-fetch the affected gages; **splice** the restored
days into the February 2026 parquet as a patched snapshot with a correction ledger, the raw
responses and a rule-vs-vintage control; rerun each product as a **frozen cohort** — the full patched input with a
deposited per-product selection manifest (`--gage-list`), so that a reader can reconstruct the
run exactly; gate the result against version 2 in two exhaustive classes — unaffected gages
cell-for-cell equal, affected gages (those whose Q actually changed inside the window) changed
only on the years the acceptance oracle names, with climate-only annual values exactly equal on
every shared valid year. Harness, validator,
rpkg-gate and staging scripts need the changes in §4 before the run. Realistic effort: about seven
person-days hands-on plus ~20 h of unattended machine time (§8), the production batches on a
quiet machine.

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
5. **Fold in the streamflow-frame climate fix** (outer join; ON HOLD since 2026-10-16)? Default
   **no** — it is a different defect (absent rows), and this fix does not touch it.
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

### Phase C — harness changes (1 day; verified gaps from the round-1 review)
- `validate_production_run.py`: `--expect-same-gages` (exact set equality + metadata equality),
  not overlap.
- `run_batched_julia.py`: `split --gage-list` (selection manifest, hashed into `batches.json`);
  `merge --reference` with identical sets needs no `--allow-missing`; resume keyed on the full
  source revision (tree hash), the Julia version + `Manifest.toml` hash, the config, metadata and
  manifest hashes and the actual batch-input hashes, with output identity re-checked; **a
  production run refuses a dirty or unknown tree** (no `--allow-dirty` for products).
- `audit_qualification.jl`: qualified / selected / emitted sets (round-2 finding 2).
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
- D0 **baseline reproduction**: the committed tree on the ORIGINAL input, both windows, all
  batches (≈ 5–7 h unattended), compared cell-for-cell with v2 — gives a tree-tied v2 baseline
  (v2 itself was built from a dirty tree) before the input changes; if it does not reproduce v2,
  the comparison reference becomes this reproduction and the drift is explained first.
- D1 clean tree at a tagged commit, `STREAMFLOW_HASH_INPUTS=1`, the two version-2 wrappers with
  `STREAMFLOW_DATA_PATH` / `STREAMFLOW_METADATA_PATH` at the spliced files, `split --gage-list`
  with the product's selection manifest, climate `daymet_1980_2025_29sep2026.parquet` (verify 5,134,193,404 B,
  `PAR1`, md5 `c0590760…`); `split` (8 batches) → `run` (check `vm.swapusage` first; raise
  `--max-minutes` if shared) → `merge --reference <v2 CSV>` (`--memory-limit 3GB` for the 46-year
  window); read `merge.log`. Product #1 v3, its gates, then #2 v3.

### Phase E — gates, comparison, artefacts (½ day)
E1 `validate_production_run.py --reference-csv <v2> --expect-same-gages`; E2
`validate_annual_values.py`; E3 `check_signature_failures.py` on every batch log; E4
`audit_qualification.jl` (qualified = selected = emitted); E5 `compare_qualifier_rerun.py --new
v3 --old v2 --oracle <csv> --ledger <csv> --registry <csv>` (§5); E6 explorers +
`check_html_basemap.py`, a v3-vs-v2 dashboard, all into the run folders; `flagged_for_high_na`
comes from the run; E7 **independent numerical check of the affected results**: the Python port
run on the affected gages of product #1 (≈ 900 gages, ~15 min) and
`check_annual_parquet_equality.py` / `check_schema_equality.py` against the Julia v3 subset —
the diagnostics of E5 (vi–vii) are accepted only with this check PASS and a written
disposition.

### Phase F — rpkg reference re-baseline (≈ 3 h of R; decision 8) with the parameterized gates.

### Phase G — adversarial reviews (1 day; decision 9), fixes, re-gates; a rerun if a fix touches the input.

### Phase H — staging, docs, close-out (1 day)
R1/R2 via the staging tool's new mode (selection manifests + hashes + the exact commands in
the run notes); R3: the **immutable February base**, the patched parquet, the ledger, the raw
responses (gzip), the three-dataset control report, the final acceptance oracle, the selection
manifests and an executable reconstruction recipe, updated metadata, README, abstract,
dictionary — the staging check verifies every one of them is present. Docs (§7). Commit →
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
- E5 two exhaustive, disjoint classes derived from the ledger (gages with ≥ 1 ledger row inside
  the window = affected; all others = unaffected), compared as **keyed cells** (gage × column;
  gage × signature × water year) with NaN == NaN, null == null, signed zeros equal:
  - **U — unaffected** (≈ 5,779 / 5,493 gages): every summary cell and every annual-parquet row
    equal to v2. Any difference FAILS (D0 tells whether the tree or the splice is the cause).
  - **A — affected** (gages with ≥ 1 ledger row whose Q actually changed inside the window;
    899 / 757 if every targeted day restores): (i) valid-year identity per gage = v2's years ∪
    the oracle's restored years, exactly; (ii) **immutable site attributes** (latitude,
    longitude, basin_area, area_normalized, gage_type, the 11 human-interference columns)
    unchanged, while the **derived qualification metadata** (`num_water_years`,
    `start_water_year`, `end_water_year` — the runner computes them from the qualifying years)
    equal the oracle; (iii) on water years with no
    ledger day and unchanged validity, **year-local** signatures equal cell for cell;
    **record-dependent** families (drought thresholds, `*_all` pulses, `elasticity_rolling`,
    parameterized BFI, `avg_storage`, per the registry) reported with R² and max |Δ| on those
    years; (iv) **climate-only** annual values (`prcp_*`, `temp_*`, `gdd_*`, `frost_*`,
    `swe_*` bases) equal on every shared valid year including years with restored Q — an
    inequality means the climate join moved; (v) newly valid years counted per gage = oracle;
    (vi) summary statistics move by construction — reported per column with R², max |Δ| and
    its gage, gained / lost cells, `flagged_for_high_na` flips; (vii) Q-dependent
    climate-coupled signatures (runoff ratios, elasticity, Q–P seasonality, storage, snow) on
    years with restored days reported separately from (iii). (vi) and (vii) are diagnostics:
    accepted only with E7 PASS and a written disposition; the family registry names the snow
    bases explicitly (several carry no `swe_` prefix).
- E6 explorers load their sidecars headlessly; basemap PASS.
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
| F rpkg re-baseline | 0.2 | 3 h | E |
| G reviews + fixes (+ a rerun if the input changes) | 1 | (rerun: 5–7 h) | E, machine quiet |
| H staging, docs, relay, commit, publish | 1 | — | G |

**About 7 person-days hands-on plus ~20 h of unattended machine time**; for one person working
sequentially, seven to eight elapsed working days (B and C overlap only with a second person).
The first draft's "two to three" did not budget the tooling; v2's "four to six" did not add
the table up.

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

No further plan review round is planned; the next reviews are of the code (decision 9).
**Decisions taken 2026-10-26** (user): the defaults of 1a–1e, 2, 5, 6, 7, 8, 9, 10 and 11, with
3 (no new gages) and 4 (no-code days untouched) as recorded. Implementation order: A → B (with
its checkpoint) and C → D0 → D → E → F → G → H.

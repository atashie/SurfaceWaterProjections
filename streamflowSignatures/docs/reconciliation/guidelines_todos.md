# Guidelines document — sync log and open TODOs

Dated entries, newest first, written by `/sync-docs` after each sync of the signature
guidelines Google Doc (`docs/SIGNATURE_GUIDELINES.md` is the snapshot). Open items carry
`[ ]`; tick them when applied in the doc or implemented in code, and move closed dated
blocks to `changelog-old.md` when they are fully resolved. Moved here VERBATIM from
`CHANGELOG.md → [Unreleased] → Guidelines Document TODOs` on 2026-09-29.

---

**Synced 2026-09-29 — 20 wording regions changed since 2026-09-10; structure unchanged
(277 vs 270 paragraphs). Snapshot overwritten.** Landed from earlier queues: Part 4
"missing years excluded" and "handles flagged data" (both 2026-09-04 items, now checked
below), Part 5 `flagged_for_high_na` redefinition (1) + the cross-language qualifier (3),
and Part 5 now says the flags "can be recomputed from an output CSV without rerunning any
signature" (true — `refresh_qa_flags.jl`, `recompute_high_na_flag.py`). Verified against
code: the §2.2.2-style constant-flow rule is not in the doc, but Part 3.8 snow rewording,
the b = 1 recession paragraph, and the 2 % reversal note all still match. New doc-side
items (code is right):
- [ ] **Part 3.3 parameterized-BFI validity interval changed from (0, 1) to [0, 1] — wrong.**
  `analyze_baseflow_indices_with_parameters` returns NA when `alpha <= 0 || alpha >= 1`
  (`julia/src/baseflow.jl:283`), so the OPEN interval was correct; revert to (0, 1).
- [ ] Part 1.2 now cites "Hecht et al. 2024" for the four QA conditions; no such entry in
  Part 6 References — add the reference.
- [ ] Part 3.3 and 3.3-recession now use the name `recession_alpha_point_cloud` twice; the
  shipped column is `recession_alpha_point_cloud_linear_reservoir` (carried forward).
- [ ] Part 5 (2): the released-product high_na known-issue sentence still not landed.
- [ ] Carried forward unchanged: "(see 3.12)" → 3.1; `avg_storage` block absent;
  `ice_affected_days_total` absent; elasticity "11 consecutive qualifying observations";
  Pettitt window paragraph in Part 2. Nit: "occurs ." (stray space) in `swe_max_dowy`.

**Synced 2026-09-10 (pm) — the doc was RESTRUCTURED to the colleague's 8 categories; full
category crosswalk run (user request).** Part 3 now has eight modules (Flow Volume, Flow
Duration, Storage, Flashiness, Drought, Flow Timing, Precipitation Streamflow, Snow) with
the former 13 modules nested as function blocks; the legacy tail and START/END markers
are gone. Snapshot overwritten. Manuscript unchanged. Live HydroShare is private and
Chrome was not connected, so the deposit side was checked on the STAGED folder. Record
and per-output CSV: `docs/plans/2026-09-10-signature-category-review.md` (afternoon
section) + `2026-09-10-signature-category-crosswalk-pm.csv`. Open items, by product:
- [ ] **Manuscript** still says 14 categories in the §2 preamble and §2.2.1 (and points to
  the dictionary for details); the §4 "sites and signature categories" figure scheme is
  unknown. Needs the 8-category wording once the vocabulary is settled.
- [ ] **HydroShare R1/R2 (staged)**: dictionary `category` + `hisss_signature_categories.csv`
  carry the 9-class grouping (7 bases differ from the doc: flashinessRB, six "Snow
  Timing" bases); the 21 scalars have no group; README lede lists the 14 families; the
  README file-table row says "nine-class … finer methodological categories are a
  separate taxonomy"; the shipped explorer / validation summary / dashboard embed the
  repo's 16- and 15-group schemes.
- [ ] **Guidelines doc internal**: `avg_storage` / `calculate_average_storage()` absent
  (sheet + dictionary + manuscript all carry storage); TQmean in 3.4 Flashiness vs Flow
  Volume in the sheet and dictionary (the one annual base where doc and sheet disagree);
  Part 1.2 "(see 3.12)" → 3.1; `ice_affected_days_total` has no module; names
  "recession_alpha_point_cloud" / "Negative_ann" vs the shipped columns
  `recession_alpha_point_cloud_linear_reservoir` / `negative_ann`.
- [x] **Repo docs on the public mirror** — DONE 2026-09-10 (later): the user ruled that
  the eight are the categories; README.md, SIGNATURES.md, CLAUDE.md, claude-skill, the
  schematic and the explorer builder are aligned (see `[September 2026]`). Manuscript and
  guidelines-doc edits catalogued in `docs/plans/2026-09-10-manuscript-category-edits.md`;
  HydroShare docs deliberately untouched (Planned → HydroShare documentation updates).

**Synced 2026-09-10 — 19 wording regions + two Part 4 corrections landed; signature-category
review (user request).** Snapshot overwritten (header lists the edits; the legacy tail no
longer carries the page JavaScript the 2026-09-04 rebuild appended). Applied in the doc
since 2026-09-04: Part 4 entry point → `process_gages_rawToRaw()` (new + legacy parts)
and the compilation gate stated correctly — the two Part 4 items below are now checked;
flashiness / flow-timing / snow definitions reworded, all still consistent with the code.
Part 5 unchanged (the three `flagged_for_high_na` edits still pending). New nit: Part 1.2
says Negative_ann is "(see 3.12)" — the module is 3.11 in the pasted doc (3.12 is Snow);
the pointer is right only under the draft's numbering that had Storage at 3.11.

*Signature categories — three schemes are in circulation and they are not the same axis*
(record + 121-row matrix: `docs/plans/2026-09-10-signature-category-review.md` and
`docs/plans/2026-09-10-signature-category-matrix.csv`):
- **14 methodological categories** — manuscript §2 preamble + §2.2.1, README.md,
  SIGNATURES.md summary table, the workflow schematic, and the code's 14 family functions
  (baseflow's two functions counted once). Every 14-list agrees.
- **13 guidelines modules** (3.1–3.13) — the same list minus Storage; `avg_storage` is
  documented only in the legacy tail. That is the whole 13-vs-14 discrepancy (queued
  2026-09-04, still open).
- **8 analysis groups** (colleague's sheet `1sFVNL9bp…`, 100 annual bases, no scalars) —
  a coarser interpretive grouping, NOT a re-cut of the 14: Baseflow + Recession + Runoff
  ratios + Storage → "Storage"; Pulses + Reversals + Elasticity + flashinessRB →
  "Flashiness"; TQmean + negative_ann → "Flow Volume"; FDC → "Flow Duration"; Q-P
  seasonality → "Precipitation Streamflow"; Snow; Drought; Flow Timing. It is a revision
  of the **9-class exploratory grouping shipped in the staged HydroShare
  `hisss_signature_categories.csv` AND as the `category` column of
  `hisss_data_dictionary.csv`**: 7 of 100 rows differ — flashinessRB moves Storage →
  Flashiness (the staged placement is an error) and the 6 "Snow Timing" bases fold into
  "Snow".
- **The inconsistency that bites**: the manuscript says "14 categories … details in
  hisss_data_dictionary.csv", but the dictionary's only `category` column carries the
  9-class grouping — a reader finds 9 categories that do not match the 14 named in §2.2.1.
  The Resource 1/2 READMEs do say the CSV is a "coarse nine-class exploratory grouping"
  and the paper's categories are "a separate taxonomy", but nothing ships the 14.
- **Proposed solution (pending user decision)**: keep BOTH axes under distinct names —
  *signature category* (14; methods, guidelines modules, code families) and *analysis
  group* (the colleague's 8; figures and exploratory dashboards) — and (1) restore a
  Storage module to guidelines Part 3 so it reads 14 (renumbering 3.11–3.13 → 3.12–3.14,
  which also repairs the "(see 3.12)" pointer); (2) adopt the colleague's sheet over the
  staged 9-class file, extended to the 21 scalars by inheriting each scalar's family
  group; (3) give the dictionary two columns (`category` = 14, `analysis_group` = 8) and
  ship both in the categories CSV; (4) add one bridging sentence to the manuscript where
  the §4 figure uses the 8 groups; (5) align the repo tooling afterwards — dashboard
  `SIGNATURE_GROUPS` (negative_ann sits in "Pulses"; 21 scalars unselectable), the
  explorer/compare categorizers (split "Flow percentiles" from "Flow volumes", scalars
  fall to "Other"), the claude-skill overview (lists 6 categories), SIGNATURES.md (no
  numbered Negative Flow Days section). Judgment calls to put to the colleague, not
  errors: elasticity under "Flashiness" while runoff ratios sit in "Storage" and Q-P
  seasonality in "Precipitation Streamflow" (three P–Q families in three groups); TQmean
  and negative_ann under "Flow Volume".

**2026-09-04 (later) — user-requested accuracy review of Parts 4 and 5 (live doc
re-fetched; three small edits since the morning sync mirrored into the snapshot).**
Findings, by direction of fix:

*Part 4 (doc-side — the code is the legacy R ingestion and behaves as read):*
- [x] **Wrong entry point named.** APPLIED in the doc 2026-09-10. The production compilation (`run_ingest_usgs_hydat.R`)
  calls `process_gages_rawToRaw()`, which writes the daily parquet + the metadata table
  (`processing_status`, `area_normalized`) and computes NO signatures. The documented
  `process_gages_rawData()` is the legacy end-to-end R path (retrieval + R signatures +
  HydroBasins lookup) that did not produce HISSS. Part 4 should describe `rawToRaw` (or
  say the R utilities only compile inputs; signatures come from the Julia/Python/R library).
- [x] **The `min_Q_value_and_days` legacy note is wrong for the compilation stage.** APPLIED in the doc 2026-09-10. The
  gate is STILL applied there: a gage is `success` only if ≥ 20 calendar years each hold
  ≥ 30 days with Q > 0.0001 (`MIN_Q_VALUE_AND_DAYS = c(0.0001, 30)`, `MIN_NUM_YEARS = 20`);
  80 gages carry `insufficient_years` ("Only N years with valid data", max 19) in the
  production metadata. What April 2026 removed is the per-year application inside the
  SIGNATURE path. The note should say so, and the manuscript's "fewer than 20 total years
  containing valid daily observations" glosses the ≥ 30-day/0.0001 rule (relay).
- [x] APPLIED in the doc by 2026-09-29. "Missing or invalid years are excluded from analysis" — not at ingestion:
  `rawToRaw` keeps every day in the window once the gage qualifies; year exclusion is the
  preprocessor's job (Part 1.2). True only of the legacy `rawData` path.
- [x] APPLIED in the doc by 2026-09-29. `generate_streamflow_dt()` "handles flagged data": the qualifier mask (keep A, A e,
  P, P e) is applied to USGS only; HYDAT `Symbol` values are carried in `flag` but never
  mask Q. The "required years" check is a row-count proxy (`nrow > 365 × min_num_years`).
- [ ] `process_caravan_gages()` "handles redundancy between CAMELS and HYSETS" —
  overstated: it only skips a (watershed, project) pair already in the output file on
  resume; the same USGS gage present in both projects is processed twice. Also
  `run_caravan_processing.R` calls `process_caravan_to_annual()` (a 328-valid-day annual
  aggregator), not `process_caravan_gages()`. Caravan is not part of HISSS — consider
  saying so or dropping the module.
- [ ] `integrate_daymet_with_streamflow()` is not on the production path either — the
  Julia runner joins climate itself (`gage_id` + `date`, keeping PPT and SWE only). The
  R function's behaviour is otherwise as described (site_id match, join on Date, < 95 %
  coverage warning, prcp → PPT). The path "data_out/daymet_1980_2023.parquet" is a
  parameter; config.R points at `D:/processedOuts_feb2026/…`, and the products used the
  `_rebuilt_10aug2026` file (the canonical file is truncated).
- Dropping `convert_daymet_zip_to_parquet()` is fine, but the doc now says nothing about
  how the climate parquet was built; the current input came from the Python
  `docs/benchmarks/convert_daymet_csvs_to_parquet.py` (365-day Daymet calendar).

*Part 5 (thresholds all match `config/signatures_config.json`; 11 of 12 flags verified
identical across languages on the reference set):*
- [ ] **`flagged_for_high_na`** — code fixed the same day (one shared definition; see
  `[September 2026]`); three doc-side edits queued for the Google Doc: (1) redefine the
  flag as ">30 % of a gage's signature-derived columns (every column except the gage
  metadata and the flag columns) are NA"; (2) append the known-issue sentence for the
  released products (column computed over the 16 metadata columns; marks every
  Canadian gage; regenerated at the next data rerun); (3) qualify "identical across the
  Julia, Python, and R implementations" with a pointer to that known issue. Optional
  nit: the `flagged_for_elasticity_range` parenthetical should read "the only range
  check not on a _mean column".
- Part 4: the six concise replacement edits are in
  `docs/plans/2026-09-04-guidelines-part4-edits.md` (delivered to the user in-session).
- [x] Everything else verified against code: ranges, the BFI defensive note (per-year
  `clamp(bfi, 0, 1)` in `baseflow.jl`), `elasticity_static` as the one non-`_mean` input,
  ties allowed in both order checks, NA never triggers a flag, seasonal-sum tolerance 0.2.
  One negligible edge: Python flags `seasonal_sum` when `Qann_mean == 0` and the seasonal
  sum is nonzero (inf ratio) where Julia/rpkg return false — unreachable in practice.
- Note: the legacy `R/tests/qa_qc_signatures.R` uses STRICT ordering (ties flag); the
  doc's "R implementation" must mean rpkg, which matches Julia/Python.

**Synced 2026-09-04 — the doc now carries the 2026-08-31 RESTRUCTURE (pasted between
`START/END: NEW DOCUMENT` banners) followed by the previous document verbatim under a
`START: ORIGINAL DOCUMENT` banner; header link now https://github.com/CZ-Sync/HISSS.**
Snapshot overwritten. Substantive deviations from the delivered draft, **reconciliation
review PENDING** (surfaced during the 2026-09-04 dataset-join session, not yet
adjudicated): the entire 3.11 Storage module (`avg_storage`) was DROPPED from the new
part (modules renumbered; avg_storage survives only in the ORIGINAL part); Part 2's
Pettitt evaluation-window sentence (WY 1980–2024, ≥ 20 / ≥ 10 obs, WY 2025 excluded)
dropped; the "T-1 moving window" baseflow bullet and the elasticity decision note
dropped; `recession_alpha_point_cloud_linear_reservoir` renamed to
"recession_alpha_point_cloud" (no longer marked per-gage scalar) — the CSV column name
is unchanged, so this is a doc-side naming discrepancy; snow Purpose rewritten (tie rule
spelled out; `swe_apr1` lost "leap-year safe"; unclosed parenthesis in melt_season_days);
Yilmaz et al. (2008) added; typo "Method and parameters→".
Adjudication queued (2026-09-04, from the census session's own diff of the pasted text vs
the draft — 48 differing regions, most of them wording):
- [ ] **Storage module dropped** — the product still ships `avg_storage` (+16 columns)
  and the manuscript (synced the same day) now lists "catchment storage" among its
  **14** categories, so the guidelines (13 modules) and manuscript disagree on the family
  count. [Doc-side: restore the module (recommended) or state that storage is shipped but
  undocumented there.]
- [ ] **Pettitt evaluation-window paragraph dropped** — a product-level fact already in
  both HydroShare READMEs and all 800 Pettitt dictionary rows; Part 2 should carry it.
- [ ] Minor wording drift (all still consistent with the code as read): Lyne-Hollick
  "a = 0.925" (conventionally alpha); drought thresholds described as "calculated
  per-water-year" — they are whole-record (fixed) thresholds applied per water year.
- The legacy document trailing the new one should eventually be deleted from the Google
  Doc so the auto-sync diff stops carrying ~4,000 words of superseded text. [Doc-side.]

<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

**Synced 2026-08-31 — the guidelines doc MOVED to its current publish URL (declared ground
truth by the user) and was RESTRUCTURED into one self-contained module per signature
family; most of the doc-side fixes queued that day were APPLIED the same day** (flow-volume
units, pulse thresholds, NA-handling four-condition wording, recession b=1 sentence, storage /
snow / drought sections, high_na wording, references, the legacy 8-stat table, typos). The
applied-items checklist, the Codex review of that sync (GO-WITH-FIXES, 4 MAJOR + 4 MINOR),
and the **user decision that every external-facing repo pointer goes to
https://github.com/CZ-Sync/HISSS** are recorded in changelog-old.md → Guidelines Document
TODOs. Still open from that sync (doc-side; the code is correct):
- [ ] **Elasticity rolling window wording** (Codex finding): the "11-year window"
  is 11 consecutive QUALIFYING observations (`elasticity.jl` indexes the
  PPT-filtered valid series), usually but not necessarily 11 consecutive water
  years — same nuance already documented for `elasticity_annual`'s adjacent
  qualifying years. In the restructure draft (§3.9); not yet in the doc.
- [ ] Recommended additions: the Pettitt window WY 1980–2024 (+ ≥20 obs / ≥10 per
  segment ⇒ WY2025 excluded from all Pettitt fields); the 20-value stats floor;
  recession/elasticity trend-gate exemptions; BFI_*_param variants in the Baseflow
  glossary; Q95_Q10; references for Eckhardt 2005 / Lyne & Hollick 1979 /
  Baker et al. 2004 / Pettitt 1979. ALL included in the restructure draft
  (Parts 1.3 / 2, §§3.3–3.5, Part 6); land in the doc when it is pasted.
- [ ] **Header repo links** (new, user decision 2026-08-31): the doc's "github
  repo here or here" links point at the private repo and CZ-Sync/code-sandbox;
  they must point at https://github.com/CZ-Sync/HISSS. In the restructure draft
  (header line); the audit-sheet link is preserved.

*Code-side items the review surfaced (already tracked; nothing new):*
- Recession R² < 0.8 fit-quality flag — still requested by the doc, still
  unimplemented (existing TODO below).
- Ice-affected day count — requested by the doc; `ice_affected_days_total` ships
  but is structurally always 0 (Known Issue, MEDIUM, 2026-08-26).

Synced 2026-07-21 — first doc revision since 2026-04-15; still-open behavior-changing items
(the two implemented items — the 80 % → 60 % overall trend gate and the confirmed 80 %
decade gate — and the documentation-only list are in changelog-old.md → Guidelines
Document TODOs):

- [ ] **Recession fit-quality flag: R² < 0.8** (source: analyze_recession_parameters).
  "To control the quality of fitted parameters, calculate recession fits and create a
  flag for any R2 < 0.8." Not currently implemented. Needs design: per-event vs
  per-gage flag, and note that b=1 alpha fits are medians (no regression R²) — the
  free-fit b regressions are the natural target. Clarify scope with domain experts.
- [ ] **avg_storage omitted from major analyses (4/23/26)** (source:
  calculate_average_storage). Doc header now says "OMITTING VARIABLE FROM MAJOR
  ANALYSES"; extensive redesign notes added (3 options incl. GLEAM ET water balance;
  Erin's dormant-season recession storage method needing PET, dormant-season
  definition, initial-storage assumption; open questions on snow). Decide: gate/flag
  `avg_storage` in outputs vs leave computed + documented as excluded downstream.
- [ ] **NA-handling wording conflict — item (i) "flag not remove"** (source: NA
  Handling). New sentence: "Items i, iii, and iv are set in the config to flag (not
  remove) by default." Items iii (negative Q) and iv (constant SD) match the current
  config-driven flag-only behavior, but item i (>3 consecutive days of NAs) currently
  REJECTS the year in `preprocess_daily_data()` (not config-toggleable). Clarify with
  domain experts whether item i should become a config-driven flag.

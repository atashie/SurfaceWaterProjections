# Changelog

All notable changes to the Streamflow Signatures project.

This file holds only the CURRENT state: `[Unreleased]` (open plans, live known issues,
open guidelines/manuscript items), the most recent month in full, and condensed
summaries of the two months before it. Everything older — the full August and July 2026
entries, June–March 2026, and the resolved/superseded `[Unreleased]` items — lives in
[changelog-old.md](changelog-old.md) (verbatim, newest first); Dec 2025 – April 2026
detail is in [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

> **Convention** — keep this file short: it is loaded into context at the start of every
> session (CLAUDE.md → `@CHANGELOG.md`). When a month closes, condense it here to a
> headline-per-change summary with pointers and move its full text to `changelog-old.md`
> (verbatim, newest first); prune `[Unreleased]` to the items still open, moving resolved
> entries and superseded dated log entries to `changelog-old.md` as well. File-level change
> lists belong in `git log`, not here; analysis and benchmark tables belong in the canonical
> docs (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`) and are linked rather than re-hosted.

## [Unreleased]

### Planned
- **At the NEXT rerun of any portion of the data: regenerate `flagged_for_high_na` in
  both standard products** (user decision 2026-09-04 — the delivered column is a
  cataloged known issue, not rewritten in place). Either the full benchmark (the code
  is fixed) or `docs/benchmarks/recompute_high_na_flag.py --write` on the existing CSVs
  produces the corrected column (product #1: 1,224 → 791; #2: 1,243 → 598); then update
  the HydroShare READMEs/dictionary row and the run-notes counts, and re-run the rpkg
  benchmark so its reference CSV picks up the constant-series Mann-Kendall fix.
- **Canonical cleanup (LOW priority, non-blocking): make Julia's storage
  `unique` use `==` semantics rather than `isequal`.** `julia/src/storage.jl`
  builds `Q_unique = unique(Q_valid)` and gates the year on
  `length(Q_unique) < 10`. Julia's `unique` compares with `isequal`, under which
  **−0.0 and +0.0 are DISTINCT**, so a year containing both signed zeros counts
  one extra "unique" value; numpy (`==`) and R (`unique`, also `==`) do not.
  Measured impact across the whole WY 1993–2025 product: **exactly 1 gage-year
  in 18.9 M** (gage 08134000, `avg_storage`, WY2017 — 365 days, 9 distinct Q
  values plus both signed zeros). The port is arguably the more correct side
  here: for a *numerical* distinctness threshold, −0.0 == +0.0. **Direction of
  fix: change JULIA to `==` semantics** (e.g. `unique(x -> x + 0.0, Q_valid)`
  or an explicit tolerance-free numeric dedup), not the ports. **Explicitly NOT
  blocking**: it does not justify re-running any benchmark or delaying the port
  campaign (user decision, 2026-08-26); fold it into the next behavior-changing
  Julia release and let it land in the next product regeneration.
- **Port campaign COMPLETE (2026-08-27)** — both ports validated at full scale against
  canonical Julia (1,653 columns / 6,678 gages; see `[August 2026]` below, full record in
  changelog-old.md). **Four rpkg-side fixes identified 2026-08-26 remain DEFERRED** because
  they change package source and need an `R CMD INSTALL` between benchmark runs (none
  invalidates a delivered run — each was checked against the run it would have touched):
  1. `rpkg/R/config.R` reads `STREAMFLOW_SIGNATURES_CONFIG`, not the canonical
     `STREAMFLOW_CONFIG` — an experiment launched uniformly with `STREAMFLOW_CONFIG=variant.json`
     would use the variant in Julia/Python and the installed default in rpkg (same hazard
     class as the Julia precompilation gotcha).
  2. Non-canonical fallbacks when the whole `na_handling` section is absent
     (`reject_negative_flow` would default TRUE; the three SWE keys would be NULL) — latent
     only, the bundled config carries the section byte-identical to canonical.
  3. Emit a per-run identifier into the log header AND a MANIFEST written after every other
     artifact (all three runners, rpkg first), so a log / timing JSON / CSV set is bound to one
     execution rather than matched by OUTPUT_PREFIX and counts.
  4. `run_rpkg_benchmark.R` merges SWE only inside the `if ("PPT" %in% clim_cols)` branch; a
     climate parquet carrying SWE but no PPT would silently drop every snow key (Python merges
     the slice independently).
  The three runner-side fixes (`options(warn = 1)`, `report_interval` 500 → 50, column
  projection of the streamflow parquet) are APPLIED.
- Restore a readable copy of the ORIGINAL `daymet_1980_2023.parquet` if one exists
  elsewhere (**the Google Drive backup** — the Shiny app formerly read this file from
  S3, so a copy existed there and may have been backed up before S3 access was lost on
  2026-08-24; else another machine or the Windows `D:` drive). **Product #2 is built on the
  CSV-rebuilt input; product #1 predates the rebuild and used the original** (each run's
  `timing.json` → `provenance` records which). An original copy would let the ≤ 3.4e-13
  replay residual be attributed rather than merely bounded, would let the truncated file
  be deleted, and would make product #1 reproducible at all — today it is not, because its
  climate input no longer exists in readable form.
- **Require a clean working tree for a standard-product run** (or retain the diff).
  Both delivered products recorded `git_working_tree_dirty = true` in their provenance
  (#1 at `b7e8988`, #2 at `0487bbd`) and no patch or source-tree hash was kept, so neither
  product can be tied byte-for-byte to a reproducible source tree — the committed code
  strongly appears to be what ran, but that cannot be *proved* from what was retained.
  The multi-GB inputs also carry size/mtime but no hash unless `STREAMFLOW_HASH_INPUTS=1`.
- **Harden the annual-values export against a silent skip.** `CFG_SAVE_ANNUAL_VALUES`
  defaults to **false when the `annual_values` config section is absent**, and it is a
  `const` baked at precompile time — so a hand-made config variant that drops the section
  would produce no annual parquet, announced only by a quiet `Annual values: save=false`
  line in the log header. Both standard products are unaffected (both logged `save=true`
  and both parquets validate), but the default should arguably be `true`, or the runner
  should warn loudly. Same hazard class as the `STREAMFLOW_CONFIG` precompilation gotcha.
- Consider making `check_additivity.jl` report a cross-machine mode explicitly — its
  shared-value gate cannot pass across machines because rank statistics flip on last-bit
  ties (measured 2026-08-11: annual series agree to 5.7e-14 while `FDC90th_spearman_pval`
  moves 0.81). Today that shows up as a bare FAIL that a reader must interpret.
- Add unit tests for core functions
- Complete `analyze_Q_PPT_relationships()` for raw data pipeline
- Add ERA5/PRISM data fetching for USGS/HYDAT gages
- Implement synchrony metrics (cross-correlation, lag analysis)
- Port data ingestion utilities to Julia (long-term — currently R-only via dataRetrieval/tidyhydat)
- Generate Julia golden outputs (624 cols, 7,313 gages)
- BFImax estimation via Collischonn & Fan (2013) backward filter — would give BFI_Eckhardt_param per-gage BFImax instead of fixed 0.8, improving discriminating power (currently range [0.47–0.80] due to BFImax saturation)

### Known Issues
- **KNOWN ISSUE in BOTH delivered products — `flagged_for_high_na` (one column).** In the
  shipped CSVs it was computed over the 16 numeric metadata columns only (the runner built
  every signature column as `Vector{Any}`, which failed the numeric-eltype filter), so it is
  TRUE for every Canadian gage + USGS gages lacking GAGES-II attributes (1,224 in product #1,
  1,243 in #2) and says nothing about signature completeness — do not use it to screen
  gages. Code fixed 2026-09-04 (one shared definition in Julia/Python/rpkg — see
  `[September 2026]`). **DECISION (user, 2026-09-04): the products are NOT rewritten; the
  column is regenerated the next time any portion of the data is rerun** (a benchmark, or
  `docs/benchmarks/recompute_high_na_flag.py --write` on the existing CSVs: #1 1,224 → 791,
  #2 1,243 → 598). Cataloged in README.md, SIGNATURES.md, CLAUDE.md, the claude-skill, the
  HydroShare Resource 1/2 READMEs and the `hisss_data_dictionary.csv` row. Full root-cause
  record (per-language denominators, how the gates missed it): changelog-old.md → Known
  Issues — resolved.
- **rpkg constant-series Mann-Kendall — RESOLVED in code 2026-09-04, rpkg benchmark rerun
  PENDING.** rpkg emitted tau = 1, p = 1 where canonical emits NaN for zero-variance annual
  series (`negative_ann` is all-zero at 6,621 of 6,678 gages; also `swe_apr1` and the drought
  p2/p5 levels — 59,191 NA-pattern mismatches the port-campaign gates missed because they
  score finite pairs only). Fixed by the exported `mann_kendall_test()` wrapper in
  `rpkg/R/stats.R` (used by `generate_stats()` and the Pettitt segment p-values; parity test
  added). The 2026-08-27 rpkg reference CSV still carries the old values until the rpkg
  benchmark is re-run. No delivered product is affected (they are Julia-built).
- **HIGH (staged HydroShare deposit, discovered 2026-09-04) — mis-padded `gage_id` in
  Resource 5 (44 ids in all three EO tables) and the Resource 4 boundary layer (9 ids).**
  9–10-digit USGS site numbers with a leading zero (e.g. `011055566`) are stored
  re-padded to 8 digits (`11055566`) because the geometry build's fallback
  `c.zfill(8)` (`EO_data_processing/geometry/build_v3.py:25`;
  `~/Downloads/geometry_rebuild/rebuild_watershed_polygons.py:137`) is a no-op on a
  stripped id that already has 8–9 digits. `canon_id` is correct everywhere. A direct
  `gage_id` join therefore silently loses **35** WY 1993–2025 and **9** WY 1980–2025
  product gages from LAI/LULC/NLCD (and 2 from the boundaries in #2), and the R5
  README's coverage counts (6,599 / 5,419 / 6,196 / 4,998) are understated by exactly
  those amounts (true: 6,634 / 5,454 / 6,205 / 5,007). Independently found the same day by
  the filtering-stage census (`docs/plans/2026-09-04-filtering-stage-census.md` §5.1).
  **DECISION (user, 2026-09-04): the staged data files stay as delivered; the join rule is
  DOCUMENTED instead** — root README → "Reading and Joining the Data" (mirrored to
  CZ-Sync/HISSS), EO README §1 correction, and the Resource 3/4/5 READMEs, the R4
  abstract and the dictionaries now all say: read ids as text, join on the zero-stripped
  form on BOTH sides, never re-pad (ids are collision-free under stripping — 17,166
  stations, 0 collisions). The R5 README coverage counts were corrected to the canonical
  6,634 / 6,205 (MODIS) and 5,454 / 5,007 (NLCD). Related, also fixed in the R3 README:
  its `zfill(8) if len <= 8` recipe for the zero-stripped metadata ids loses 57 / 29
  product gages (66 of the 8,014 compiled gages) for the same reason — replaced with
  strip-leading-zeros-on-both-sides. The draft §5.1.3 join text in
  `docs/plans/2026-09-04-dataset-join-guidance.md` §4 was revised to the
  strip-on-both-sides wording the same day; manuscript §3's key-fields sentence
  ("All tables join on gage_id, the zero-padded …") needs the same qualification —
  replacement text delivered 2026-09-04 (plan doc §4b), pending paste.
- **ACCEPTED DIVERGENCE (rpkg vs Julia, 1 annual row in 18.9M) — the storage
  distinct-value guard and signed zeros.** `julia/src/storage.jl` gates a year on
  `length(unique(Q_valid)) < 10`, and Julia's `unique` uses `isequal`, under which
  `-0.0` and `+0.0` are DISTINCT. rpkg's `unique` (like numpy's) uses `==`, so a year
  holding both signed zeros counts one fewer distinct value and can be skipped where
  Julia keeps it. **A pre-run review recommended replicating Julia's signed-zero
  distinctness in rpkg so the strict annual gate passes. Deliberately NOT done**: the
  user already ruled (2026-08-26) that the ports are the more correct side here and
  that JULIA should move to `==` semantics, so reproducing the quirk downstream would
  propagate a behaviour already agreed to be wrong, in two languages instead of one.
  The cost is bounded and measured at exactly **1 annual row in 18,898,406**. It is
  recorded at gate time with an explicit, named waiver
  (`--allow-key-diff 1 --allow-diff-signature avg_storage`) rather than hidden, and
  closes for good when the canonical `unique` is fixed.
- **MEDIUM (canonical Julia, discovered 2026-08-26) — `ice_affected_days_total` is
  structurally always 0.** The WY 1993–2025 reference run emits `0.0` for **all 6,678
  gages**, yet the streamflow parquet contains **3,693 rows flagged `'P Ice'` across 90
  gages, every one with a NULL Q** — exactly the days the diagnostic is meant to count.
  Three of those gages (01011000, 01029200, 15129500) are in the product and still
  report 0. The Julia runner does load the `flag` column (it reads the full parquet),
  so the guard in `julia/src/io.jl:365` (`na_count > 0 && "flag" in names(joined)`)
  should fire; the cause is somewhere between the daily-grid normalization and that
  check, and is not yet pinned down. **Consequence**: the metric conveys no information
  in any delivered product, and the guidelines' request for "a total count of the number
  of days that are ice affected" per site is unmet.
  **Deliberately NOT worked around in rpkg** (2026-08-26): rpkg's preprocessor now
  tracks `na_cause_ice` correctly, but the benchmark runner does not load `flag`, so
  rpkg emits 0 and MATCHES canonical. Making rpkg "more correct" than Julia would break
  cross-language parity — fix the canonical side first, then enable it in both.
- **RESOLVED 2026-08-11 (workaround in place; the corrupt file is still on the drive) — the
  canonical `daymet_1980_2023.parquet` is TRUNCATED on the exFAT thumbdrive** (1,261,436,928
  of the 4,125,630,653 bytes recorded in the 28 Jul provenance block, mtime unchanged, no
  `PAR1` footer; the provenance block was what caught it). It is left in place at the user's
  direction; **every climate/SWE run must set `STREAMFLOW_CLIMATE_PATH` to
  `daymet_1980_2023_rebuilt_10aug2026.parquet`** (rebuilt from the 44 annual CSVs; product #1
  replays against it to ≤ 3.4e-13). The WY 1980–2025 wrapper defaults to it. Discovery and
  recovery record: changelog-old.md → Known Issues — resolved, and → August 2026.
- **LOW (raw/legacy path only, discovered 2026-07-27 by the drought smoke run) —
  `calculate_negative_days` crashes on `missing` Q and silently drops its 8 columns.**
  `julia/src/pulses.jl:323` applies `:Q => (q -> sum(x -> !isnan(x) && x < 0, q))`
  straight to a `Union{Missing,Float64}` column; `!isnan(missing)` is `missing`, so
  `missing && …` throws `TypeError: non-boolean (Missing) used in boolean context`. The
  orchestrator's per-signature `try/catch` turns this into a warning, so the gage just
  loses all 8 `negative_ann` columns (smoke gage 01073000 emitted 808 signature columns
  vs 816 for the other nine). Every other signature routes Q through `coalesce_q`; this
  one does not. **Production output is unaffected** while `use_legacy_filtering: false`,
  because `preprocess_daily_data()` emits Float64 Q with NaN rather than `missing` — it
  only bites callers passing raw frames (smoke test, direct API use). Fix is one line
  (`coalesce_q(df.Q)` before the group-by, matching the rest of the codebase); left
  unapplied because it is outside the drought work's scope. **VERIFIED 2026-08-24: the
  ports do NOT mirror the pattern** — Python's `(g["Q"] < 0).sum()` is NA-safe by pandas
  comparison semantics, and rpkg guards explicitly with `!is.na(q)`. One nuance to pin
  when touched: rpkg's `aggregate` formula interface drops an all-NA water year entirely,
  where Julia's groupby retains the group — add a parity test with the eventual fix.
- **LOW (legacy shim only, discovered 2026-07-21) — 6 pre-existing failures in the
  legacy R NA-handling test suite** (`R/tests/test_na_handling.R`, run per its
  documented usage after sourcing `config.R` + `R/helperFunctions.R`): grid
  normalization ("Missing date rows filled"), "3-day gap: year accepted",
  constant-SD flag, and raw/residual NA diagnostic counts. Verified present at
  HEAD before the 2026-07-21 trend-gate work (which only touched — and fixed —
  the vacuous trend-completeness case). Indicates drift between the deprecated
  `R/helperFunctions.R` shim and the evolved test expectations; rpkg is the
  active R implementation and its testthat suite is unaffected. Triage when the
  legacy shim is next touched, or retire the legacy suite with it.
- **MEDIUM (documented limitation, by design) — 37 Canadian gages in the signature
  output carry raw m³/s units (`area_normalized = FALSE`)**: HYDAT publishes NO
  drainage area (neither `DRAINAGE_AREA_GROSS` nor `DRAINAGE_AREA_EFFECT` — verified
  directly against `Hydat.sqlite3`) for 73 successfully-processed Canadian stations;
  37 survive the 20-year filter into the Julia canonical output (7,313 gages).
  Station names show most are not natural watersheds: ~40 irrigation/diversion canals
  + ditches, ~15 dam/powerhouse outflows, several huge-river channel splits and lake
  outlets (St. Lawrence, Mackenzie, Nelson, Lake of the Woods); only ~8 look like
  natural streams. 62/73 are `REGULATED = TRUE`.
  **DECISION (user, July 2026): keep these gages with raw m³/s flow — NO area
  backfill** (HydroBasins `UP_AREA` was assessed on 1,383 validation gages: accurate
  only for main-stem dam outflows, wrong for canals and channel splits, unusable
  <100 km²). Q-to-PPT signatures are now structurally gated for these gages (see the
  July 2026 entry below). **Remaining limitation**: unit-carrying Q-only signatures
  (Q volumes, percentiles, Q95_Q10, log_a) stay in m³/s for these 37 rows —
  incomparable with mm/day gages (Qann_mean up to 3.18M for the St. Lawrence).
  **Flag gap**: `flagged_for_qann_range` catches only 27/37 — 10 small canals/creeks
  land inside [0, 2000] unflagged; downstream users must filter on
  `area_normalized == TRUE` before any cross-gage comparison of unit-carrying
  signatures. See docs/DEVELOPMENT.md → Canadian HYDAT → Missing drainage areas.
- **RESOLVED 2026-08-24 — seasonal runoff ratios ignored the seasonal completeness flags
  in Julia** (`runoff_ratios.jl` looked up `winter_complete` while the preprocessor emits
  `win_complete`; Julia-only — Python and rpkg were internally consistent). Fixed Julia-first
  with a regression test; verified inert for the delivered WY 1993–2025 product (zero masking
  events in that window). Record: changelog-old.md → August 2026.

### Guidelines Document TODOs
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
- [ ] "Missing or invalid years are excluded from analysis" — not at ingestion:
  `rawToRaw` keeps every day in the window once the gage qualifies; year exclusion is the
  preprocessor's job (Part 1.2). True only of the legacy `rawData` path.
- [ ] `generate_streamflow_dt()` "handles flagged data": the qualifier mask (keep A, A e,
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

### Manuscript Reconciliation Log

Session-start reconciliation of the HISSS manuscript draft (Scientific Data,
submission target Nov 9 2026) against code + repo docs. Snapshot:
`docs/MANUSCRIPT_DRAFT.md`; workflow: CLAUDE.md → Session-Start Workflow → B.

**Earlier passes are in changelog-old.md → Manuscript Reconciliation Log** — the
2026-07-21 baseline with the eight queued manuscript edits, the 2026-08-24 second pass and
§2.1.4 / §3 drafting support, the 2026-08-28 third pass with its numbered open-item list
1–12, and the 2026-09-01 fourth pass with the paste-ready correction blocks 1a–7
(`docs/plans/2026-09-01-manuscript-correction-blocks.md`) and the §2 workflow schematic
(`docs/plans/dataset_workflow_schematic.md`). The queued-edit and item numbers cited below
refer to those lists. One thread carried forward from 2026-09-01: the committed schematic
fails GitHub's rich rendering because of a GitHub-side mermaid bundle crash (verified
against GitHub's own mermaid README, not this file) — re-check after GitHub ships a fixed
bundle; until then the 3x PNG beside the doc is the review copy.

**2026-09-10 — sync; no methods change.** Since the 2026-09-04 second sync: §3
Resources 1–2 now reads "21 signatures and related outputs that do not carry statistics"
(was "per-gage scalar diagnostics"; consistent with §2.2.1's "21 stand-alone diagnostic
metrics"); **§5.1.2 "How to join HISSS with other datasets" landed** (the
`docs/plans/2026-09-04-dataset-join-guidance.md` §4 paragraph with co-author edits) and the
record-dependent note moved to §5.1.3; the Acknowledgements double period is fixed. Snapshot
overwritten. Relay (manuscript-side): §5.1.2's example id "0103500" has seven digits (USGS
site numbers have eight or more — e.g. 01013500); "may be read as a character string" should
be "must"; an unclosed parenthesis "(most US site numbers begin with a zero … (see resource
READMEs for details)."; §5.1.3 "should not be compared *against* the WY 1993–2025 and
WY 1980–2025 products" → "*between*"; new §2.1.1 typo "7, 964". Still open from earlier
passes: avg_storage sentence (item 6), "Claude Code 0.145.0" (item 7), references +
"Linke et al. 2013" (item 5), CC-BY/DOI labelling (item 11), §2.1.2 end date,
"(e.g., CAMELS-Chem,", ">6,000stream". **Category consistency**: §2.2.1's 14 categories
match the code and README exactly; the guidelines doc has 13 modules and the shipped
dictionary a 9-class grouping — see Guidelines Document TODOs (2026-09-10) for the
proposed two-axis resolution; if the §4 "sites and signature categories" figure uses the
colleague's 8 groups, §2.2.1 or §4 needs one sentence saying the 14 categories are
aggregated into 8 analysis groups for presentation.

**2026-09-04 — filtering-stage watershed census + fifth reconciliation pass (user
request: count the watersheds at every stage of the workflow schematic and confirm the
manuscript matches).** Full record: `docs/plans/2026-09-04-filtering-stage-census.md`;
new tools `docs/benchmarks/qualification_census.jl` (runs the canonical preprocessor on
ALL 8,014 gages and re-implements the runner's inclusion gates — reproduces both products'
gage sets EXACTLY, 0 mismatches) and `summarize_qualification_census.py`. The manuscript
was revised by the co-authors since 2026-09-01 (snapshot overwritten): correction blocks
1a/1b/1c/2a/2b/3/4/5/6 are APPLIED (121 signatures / 14 categories incl. drought,
agency-published drainage areas, Daymet coverage 5,965 / 5,517 / 5,638, duplicate BFI
paragraph deleted, avg_storage listed, new Usage Notes 5.1.1 + 5.1.2, HISSS repo link,
"eight" MODIS schemes); block 7 (Acknowledgements) is NOT — "Claude Code 0.145.0" still
appears twice. **Every checkable count in the published text is CORRECT** (16,994 =
9,154 + 7,840; 8,980 excluded; 8,014 = 6,160 + 1,854; 111,624,189; 73 / 32 / 28;
6,678 / 6,250; 54 / 7,964; 100 % / 98.3 %; 6,087 / 5,965 / 5,517 / 5,638; 6,119;
1,653; 18,898,406 / 24,366,487; 2,150,280; 191,136; 250,879). Measured stage counts not
in the text, for the figure/Data Overview: 33,732 of 317,182 gage-years (10.6 %)
rejected by the preprocessor over WY 1980–2025 (31,091 >30 NA days, 1,803 gap >3 d, 838
boundary NA); within the products only 3.2 % / 2.9 % of gage-years are rejected;
exclusions 1,336 (924 both gates, 412 only <20 yr, 0 only <60 %) and 1,764 (1,063 only
<60 % with ≥20 valid years, 674 both, 27 only <20 yr); trend statistics survive for
6,183 / 5,851 gages on the dense signatures. Findings, by direction of fix:
- **DEPOSIT DEFECT (staging, fix before upload)**: `gage_id` is NOT zero-padded for 9
  gages in Resource 4's boundary layer (+ its QA CSV) and 44 gages in each Resource 5
  table (`zfill(8)` cannot restore the leading zero of 9–10-digit USGS ids), so the
  manuscript's and the READMEs' "join directly on gage_id" claim fails for those gages
  and the R5 README's coverage counts (6,599 / 6,196 MODIS, 5,419 / 4,998 NLCD — also
  quoted in this file's 2026-08-25 R5 entry) undercount; canonical values are 6,634 /
  6,205 and 5,454 / 5,007. `canon_id`, HydroATLAS and Resource 3 are clean. User decision same day: keep the
  files, document the strip-on-both-sides join rule in both repos' README and the
  resource READMEs (done — see Known Issues).
- Manuscript (relay): §2.1.2's "8,980 … fewer than 20 total years" bucket includes 113
  USGS gages whose retrieval never completed (`processing` status; 4 of them carry
  GAGES-II polygons and sit in the EO layers) — add "or could not be retrieved"; the
  year-rejection rule list omits the third rule (≤3-day NA run touching Oct 1 / Sep 30
  is not interpolated ⇒ `residual_na`, 838 gage-years); §3's "211 watershed-scale
  attributes" = 211 columns (198 attributes + 13 keys/diagnostics); new typos "wtih",
  "diagnostics.We", plus the surviving "s(e.g.,", "MOIDS", "fitted checksums …" garble,
  "Linke et al., 2013", "[tbd]" citations and the missing references (McMillan,
  Hatchett, Petersky & Harpold, Adelsperger, Laaha, Peters & Aulenbach, HYSETS, Caravan,
  CAMELS-SPAT, DeCicco, Albers, gdptools, Annual NLCD).
- Code/docs-side: none. The schematic is not contradicted by any measured count.

**2026-09-04 — co-author revision synced; dataset-join guidance drafted (§5 stub).**
The published draft changed since 2026-09-01: §2 preamble and §2.1.2 now state the
agency-published drainage areas as the normalization source (queued #1 core error
RESOLVED; the §5 caveat also landed as 5.1.1); counts corrected to **121 signatures /
14 categories incl. drought** and the "(name)" placeholder resolved to
`hisss_data_dictionary.csv` (queued #8 + open items 1a/1b RESOLVED; the drought
*methods* paragraph 1c is still absent); §2.1.3 Daymet coverage numbers added
(6,087 candidates / 5,965 of 8,014 / 5,517 and 5,638 — matches the repo); the 47 %
sentence replaced by the explicit "8,980 of 16,994 excluded" (item 3 RESOLVED); §2.2.1
/ §2.2.2 rewritten and the duplicate BFI-statistics paragraph deleted (item 4
RESOLVED); §2.1.4 "all five" → "all eight" (the 2026-09-01 finding RESOLVED); §3
mechanical typos fixed (`wy{window}`, WGS84, teh/ot, "[tbd])"); `<LINK TO REPO>` filled
with https://github.com/CZ-Sync/HISSS; **§5 Usage Notes drafted** — 5.1.1 unnormalized
gages and 5.1.2 record-dependent signatures (both verbatim from block 2b), with the
dataset-join item still a one-line stub. Still open: avg_storage sentence (item 6),
Acknowledgements "Claude Code 0.145.0" (item 7), missing references + "Linke et al.
2013" (item 5), CC-BY/DOI labelling (item 11), §2.1.2 end date, "(e.g., CAMELS-Chem,"
unclosed paren, ">6,000stream".

*Dataset-join guidance (user request, same day)*: CAMELS-Chem, MacroSheds, and the
Daymet-VPD product were reviewed against their own documentation and files, the
identifier conventions of HYSETS / CAMELS / CAMELS-SPAT / Caravan / CANOPEX / GAGES-II
/ WQP / NLDI→COMID / HydroBASINS were verified, and the actual join behaviour of the
five staged HydroShare resources was audited. Measured overlaps: CAMELS-Chem 507 / 516
(WY 1993–2025) and 515 / 516 (WY 1980–2025), all GAGES-II Ref; MacroSheds stream
gauges 146 / 224 inside ≥ 1 product watershed but almost always *nested* (median
smallest containing polygon 799 km²), only 11 co-located within 500 m. The audit found
two deposit defects (Known Issues, HIGH: mis-padded `gage_id` in R5/R4; R3 README
recipe). Verified facts, numbers, and the first-draft §5.1.3 paragraph:
`docs/plans/2026-09-04-dataset-join-guidance.md`.

**2026-09-04 (second sync, afternoon) — user edits; manuscript otherwise unchanged.**
Two typo fixes only (a missing space after "diagnostics." in the §2 preamble; "MOIDS" →
"MODIS" in §2.1.4); no methods claim changed, nothing to reconcile. The §5 "How to join
with other datasets" item is still the one-line stub — the paste-ready §5.1.3 text
(revised to the strip-leading-zeros-on-both-sides join rule after the user's decision
to leave the staged deposit files as delivered) is in
`docs/plans/2026-09-04-dataset-join-guidance.md` §4. The snapshot rebuild also removed
~170 lines of Google's page JavaScript that the morning rebuild had appended below the
body (the extractor did not skip a `<script>` block inside the contents div; the
paragraph-level diff was unaffected because the JS lines had no counterpart on either
side).
Later the same day the user asked for a replacement for §3's "Conventions and key
fields" join sentence (their interim edit — "All tables join on a zero-stripped gage_id,
the zero-padded agency station identifier …" — contradicts itself); two variants were
delivered and recorded in `docs/plans/2026-09-04-dataset-join-guidance.md` §4b. A
six-minute re-poll of the published copy showed no change, so the snapshot stands and
the next sync will confirm what landed.

---

## [September 2026]

### Fixed (HIGH): `flagged_for_high_na` now has ONE definition in all three languages (2026-09-04)
Found during the guidelines Parts 4–5 accuracy review and independently confirmed by a
Codex adversarial review (8/8 findings CONFIRMED, GO-WITH-FIXES) before the user
approved the fix. The flag meant three different things — Julia (the delivered
products) counted NA over the **16 numeric metadata columns only** (every signature
column arrived as `Vector{Any}` from the runner and failed the numeric-eltype filter,
so 1,224 = every Canadian gage + 41 USGS gages missing GAGES-II attributes), Python
counted every non-flag column incl. string metadata, and rpkg counted only the
`_mean`/`_median` keys each gage happened to emit. The other 11 flags agreed
cell-for-cell across the three languages.

- **Definition (config `qa_qc.high_na_denominator`, one manifest shared by
  Julia/Python/rpkg):** the fraction of the SIGNATURE columns present in the assembled
  output table — every column ending in one of the 16 statistic suffixes plus the 21
  registered per-gage scalars (1,621 of the 1,653 product columns) — whose value is
  NA or NaN, flagged when strictly greater than `max_na_fraction` (0.30). Metadata,
  unregistered diagnostics and flag columns never enter the denominator; families
  NA-filled by the table union DO count (that structural NA — no Daymet, no SWE — is
  exactly what the guidelines say the flag surfaces). Column selection is by NAME and
  the NA test is value-level, so untyped columns can no longer drop out.
- **Code**: `julia/src/qa_qc.jl` (+ `high_na_denominator_columns`, exported),
  `python/streamflow_signatures/qa_qc.py`, `rpkg/R/qa_qc.R` (now also accepts the
  assembled data.frame; the per-gage list path keeps emitted-key semantics, matching
  Julia's `include_qa_flags` path), config constants in all three languages with
  identical fallbacks, the manifest synced into the two bundled config copies, and
  `run_rpkg_benchmark.R` now computes the flags on the assembled table like the Julia
  and Python runners. **Tests** (`julia/test/test_qa_high_na.jl`,
  `python/tests/test_qa_flags.py`, `rpkg/tests/testthat/test-qa_flags.R`) go through
  the runner's untyped-column shape, pin metadata exclusion, NaN-counts-as-NA, the
  1,621-of-1,653 selection, and (Python) that the manifest equals the schema
  registries. Suites: Julia green (new file 14/14), Python 160 passed, rpkg 1,078
  passed against the installed package.
- **Measured on the WY 1993–2025 reference set** (new rule applied to each
  language's own CSV): Julia 1,224 → **791**, Python 787 → 791, rpkg 90 → 771 — Julia
  and Python now agree on all 6,678 gages. The 20 rpkg differences are a SEPARATE
  pre-existing divergence surfaced by the new rule — rpkg emitted tau = 1, p = 1 for
  constant annual series where canonical emits NaN — **also fixed the same day**
  (`mann_kendall_test()` wrapper in rpkg with the canonical NA contract, used by
  `generate_stats()` and the Pettitt segment p-values; parity test added; rpkg suite
  1,098 passed). The rpkg reference CSV still needs a benchmark re-run to reflect it.
- **Delivered products (dry run of `docs/benchmarks/recompute_high_na_flag.py`, a
  byte-preserving text-level rewrite of that one column):** WY 1993–2025 product #1
  1,224 → 791 (1,011 Canadian gages un-flagged, 611 USGS gages newly flagged, 33 USGS
  un-flagged); WY 1980–2025 product #2 1,243 → 598. **Deliberately NOT applied
  (user decision, 2026-09-04)**: the delivered CSVs and their staged HydroShare copies
  keep the column as shipped; it is cataloged as a known issue in the product docs
  and the dictionary row, and will be regenerated at the next rerun of any portion of
  the data. `recompute_high_na_flag.py` stays ready for that moment.
- Tooling: `refresh_qa_flags.jl` gained `--dry-run` (it re-serializes every value, so
  it is the canonical cross-check, not the rewrite tool).

## [August 2026]

Condensed summary — **the full text of every entry is in [changelog-old.md](changelog-old.md) → [August 2026]**.
The month closed the Python/rpkg port campaign, staged the entire HydroShare deposit,
published the code to the public mirror, and regenerated standard product #2 after the
climate input was found truncated.

- **2026-08-28 — pre-publication audit + PUBLIC CODE RELEASE** at
  https://github.com/CZ-Sync/HISSS (snapshot mirror via `publish_to_hisss.sh`; MIT license
  everywhere by user decision). Five parallel audit passes: zero sensitive-content blockers;
  stale port-status claims corrected in every doc; all three package READMEs brought to the
  validated-parity state; ~29 legacy files removed; Python `validate_schema()` reworked to
  the full 1,653-column schema. Codex review GO-WITH-FIXES (4 MAJOR + 6 MINOR, all fixed).
- **2026-08-27 — rpkg VALIDATED at full scale; the port campaign is COMPLETE.** 1,653
  columns × 6,678 gages, 0 errors, all four gates green (strict schema, swallowed-failure
  log scan, annual-parquet equality, identity-R² 1,601 Perfect / 10 Good / 9 Poor / 0 below
  0.95). The first rpkg run passed its unit suite yet failed three gates on real defects —
  zero-row crashes in snow, a CONUS↔AKHIPR metadata `intersect()`, pre-allocated annual
  frames, and `smooth_daily_flow` using R's `mean()` (last-bit differences flipping whole
  drought plateaus through the strict `<`; a sequential sum made it bit-identical). Column
  projection of the streamflow parquet cut the run from ~22 h projected to ~2 h. Accepted
  and named at gate time: the signed-zero `unique` divergence (1 annual row in 18.9 M).
- **2026-08-26 — PYTHON validated at scale** (1,615 Perfect / 5 Good / 0 below 0.99, mean
  R² 0.999988, 18,898,405 of 18,898,406 annual rows shared). Mann-Kendall p-value convention
  settled: scipy omits the continuity correction under ties, so Python was moved to the
  canonical Julia/R formula (user decision; bit-exact on 12 fixtures). Two rpkg wiring
  defects unit tests could not see: the runner passed none of the ported kwargs (a 19-hour
  run discarded) and the drought family returned all-NA on every real gage (`date` vs
  `Date`). New gate `check_signature_failures.py` scans run logs for swallowed family
  exceptions and refuses R's truncated-warning banner.
- **2026-08-24 → 26 — port campaign Phases 0–5**
  (`docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md`): both runners rebuilt
  with Julia-mirrored ENV overrides and bounded-memory climate handling; Pettitt + stats
  floor, b=1 recession alpha, annual collector, snow and drought ported to Python and rpkg
  with fixture cross-checks to ≤ 4e-14. Defects found on the way: flashiness guarded `== 0`
  not `<= 0` in both ports; Python `qp_seasonality`/`storage` dropped the Pettitt fields;
  `apply_stats_floor_mask.py` used pandas' inexact float parser; Python's mid-event
  day-of-water-year was one day late for even-length events; flashiness/FDC used
  placeholder metric names that would have mislabeled four annual-parquet signatures;
  rpkg's 12 QA-flag columns were double-prefixed. `check_schema_equality.py` added as the
  strict column/gage-set gate (the `compare_*` scripts are diagnostics only).
- **2026-08-25 — HydroShare COLLECTION created**
  (https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/; HydroShare is
  ground truth once a resource is uploaded and verified — user decision) and **all five
  resources staged + adversarially reviewed**: R1–2 signatures (READMEs + a 1,653-row data
  dictionary; Pettitt fields are NA outside the WY 1980–2024 changepoint window, so WY 2025
  is excluded from them), R3 inputs (`hisss_gage_metadata.csv` stores US ids ZERO-STRIPPED —
  documented, not re-padded), R4 geometry + HydroATLAS (the 7,964-basin layer had existed
  only on S3 and was **rebuilt the same day**, every June target matched exactly; the
  flash-drive GAGES-II zip was exFAT casualty #3), R5 EO products (LAI/LULC/NLCD
  md5-verified; NLCD human QA signed off). Product #1's explorer was found silently
  truncated on the thumbdrive (exFAT casualty #2) and restored from the Drive backup.
- **2026-08-24 — seasonal runoff-ratio completeness masking was dead code in Julia**
  (looked up `winter_complete`, the preprocessor emits `win_complete`); Julia-only, fixed
  with a regression test, inert for the delivered WY 1993–2025 product (zero masking
  events; reference run `processedOuts_portref_24aug2026`).
- **2026-08-24 — S3 access LOST**; the project Google Drive folder is the off-site backup
  (inventory unverified). The Shiny app is non-functional as deployed; final dashboard scope
  decided (signatures, annual values, raw series, MODIS/NLCD tables; a performant rebuild,
  not a port) and **app work deferred**.
- **2026-08-11 — STANDARD OUTPUT #2 regenerated with the drought family**
  (`processedOuts_1980_2025_11aug2026`, 6,250 × 1,653, gage set identical to the July build;
  same-machine additivity PASS — 165 columns added, 1,487 shared columns bitwise unchanged).
  Cross-machine diffs were partitioned with controls: the Windows → M1 machine change alone
  reproduces the failure pattern in discretely FP-sensitive statistics (rank stats, `TQmean`,
  Pettitt locations, `FDC90th`) while the annual series agree to ≤ 1.8e-12.
- **2026-08-11 — Daymet climate parquet REBUILT** from the 44 annual CSVs with the new
  `docs/benchmarks/convert_daymet_csvs_to_parquet.py` (stricter than the R original;
  `site_id` kept as a string; correctly-rounded float parsing). **Daymet publishes a 365-day
  calendar** (leap years drop Dec 31). Validated by replaying product #1: 0 columns
  added/dropped, identical gage set, ≤ 3.4e-13 on 98 climate-derived columns.
- **2026-08-10 — STANDARD OUTPUT #1 promoted with the drought family**
  (`processedOuts_drought_28jul2026`, 6,678 × 1,653; explorer + comparison dashboards in the
  run folder). Dashboard tooling gap fixed: `SIGNATURE_GROUPS` had never included the snow
  family, so every pre-2026-08-10 dashboard silently omitted its 14 bases.

## [July 2026]

Condensed summary — **the full text of every entry is in [changelog-old.md](changelog-old.md) → [July 2026]**.
Julia gained the annual-values export, the b=1 recession alpha, the snow and drought
signature families, the 20-value stats floor and the record-anchored snow gate; both
standard products were produced; the Annual NLCD product was built.

- **Streamflow drought family (10 metrics + 5 threshold scalars, +165 → 1,653 columns;
  Julia, 2026-07-27/28)** after Adelsperger et al. (in review): 7-day centered smoothing
  within contiguous date runs, whole-record Weibull (type-6) thresholds at the five USDM
  levels, strict `<`. Scope decisions (user): fixed thresholds only, water-year aggregation,
  NaN below the plotting range, `_fixed_` infix. Measured: `drought_duration_fixed_p10` is
  largely redundant with the low-pulse pair (within-gage median r = 0.994) — **kept by user
  decision**. Same-machine additivity PASS (1,487 shared columns bitwise unchanged) via the
  new `check_additivity.jl`; the benchmark timing JSON gained a **provenance block**; the
  GAGES-II directory is now resolved at runtime (`gages_ii_dir()`) after the
  precompile-constant trap bit three times in one day. Two Codex reviews, GO-WITH-FIXES each.
- **Annual NLCD per-watershed product** (CONUS, 30 m, 1985–2025; 6,119 gages × 41 years =
  250,879 rows; 16 classes + impervious) — `EO_data_processing/README_NLCD.md`.
- **STANDARD OUTPUT #2 — WY 1980–2025 @ 60 % (2026-07-22)**: 6,250 gages × 1,488 columns,
  Codex results review GO with zero findings; neither standard product is a subset of the
  other (window-start-anchored 60 % denominator).
- **STANDARD OUTPUT #1 — WY 1993–2025 @ 60 % (2026-07-22)**: 6,678 × 1,488; first production
  use of the `STREAMFLOW_END_WATER_YEAR` cap; explorer extended to all 16 statistics.
  **New convention (user)**: every artifact of a run lives in that run's own folder.
- **Record-anchored decade gate for the 10 timing/melt/regime snow metrics (2026-07-22)** —
  linked to the streamflow `decade_min_fraction` knob; NaNs the 6 trend stats only. Codex GO.
- **Trend-completeness overall gate 80 % → 60 % (all languages, 2026-07-21)** — config-only,
  per the guidelines doc and manuscript §2.2.3; the decade gate stays 80 %. Reinstall rpkg
  before its benchmarks (it bundles the config).
- **Stats floor — 20 annual values before ANY statistics (Julia)**; recession and elasticity
  exempt; post-hoc `apply_stats_floor_mask.py` + `refresh_qa_flags.jl`.
- **Production rerun harness**: ENV overrides for input paths / output dir / end year, memory
  patches for the 16 GB machine, `validate_production_run.py`, `audit_qualification.jl`, the
  signature explorer (`build_signature_explorer.py`). Finding: Daymet covers Canadian gages.
- **Q-to-PPT unit gate for un-normalized gages (Julia + Python + rpkg)** —
  `area_normalized = false` skips runoff ratios, elasticity, Q-P seasonality and storage; the
  runners read the flag from metadata (leading-zero-safe). Codex: 2 MEDIUM fixed, plus a
  latent InlineStrings typing bug caught during verification.
- **Snow metrics family (14 metrics, Daymet SWE; Julia)** — SWE ≥ 10 mm threshold,
  anchor-spell timing, SSM (Hatchett 2021); preprocessor `valid_swe_years`; runs only on an
  explicit `snow_data` frame.
- **Recession alpha assumes a linear reservoir (b = 1; Julia)** —
  `log(a) = median(log(-dQ/dt) - log(Q))`; b and concavity keep their free fits; column
  names unchanged.
- **Annual values export (Julia)** — opt-in `AnnualCollector` →
  `{prefix}_signatures_annual.parquet` (`gage_id, signature, water_year, value`); config
  `annual_values.save`.
- Docs: `docs/DATA_SOURCES.md` — inventory of the 11 external data sources.

## [June 2026] – [March 2026]

Condensed summaries of these months — the MODIS LAI/LULC EO products (June); the
HydroATLAS watershed metadata + static HTML explorer (May); the April 2026 release (Pettitt
changepoints, recession-parameterized BFI, Section 3 signatures, the Julia-canonical
transition, centralized NA handling, cross-language alignment); the rpkg package + alignment
rounds (March) — are in [changelog-old.md](changelog-old.md). Full per-change detail for
Dec 2025 – April 2026 is in [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

---

## Version History Notes

This project uses date-based versioning (MONTH YEAR) rather than semantic versioning, reflecting its nature as a research tool with continuous development.

### Output File Naming Convention
Output files include date stamps: `streamflow_signatures_full_JAN2026.csv`

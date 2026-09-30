# Changelog

All notable changes to the Streamflow Signatures project.

This file holds only the CURRENT state: `[Unreleased]` (open plans, live known issues),
the most recent month in full, and condensed
summaries of the two months before it. Everything older — the full August and July 2026
entries, June–March 2026, and the resolved/superseded `[Unreleased]` items — lives in
[changelog-old.md](changelog-old.md) (verbatim, newest first); Dec 2025 – April 2026
detail is in [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md). Open guidelines and
manuscript items live in `docs/reconciliation/`; the one-line project status that Claude
loads every session is `docs/STATUS.md`.

> **Convention** — keep this file short. Since 2026-09-29 it is NOT auto-loaded; the one-liners
> Claude reads every session live in `docs/STATUS.md` (rule `.claude/rules/changelog.md`).
> When a month closes, condense it here to a
> headline-per-change summary with pointers and move its full text to `changelog-old.md`
> (verbatim, newest first); prune `[Unreleased]` to the items still open, moving resolved
> entries and superseded dated log entries to `changelog-old.md` as well. File-level change
> lists belong in `git log`, not here; analysis and benchmark tables belong in the canonical
> docs (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`) and are linked rather than re-hosted.

## [Unreleased]

### Planned
- **Daymet climate input reprocess (options review 2026-09-29; single-year gate PASSED the
  same day; D1–D4 and the locations DECIDED 2026-09-29; the 1980–2025 file BUILT and
  validated 2026-09-30 — see `[September 2026]`; NEXT: the Phase 1 replay against
  product #1, then the user's decision on rerunning the products; action plan §0).**
  Daymet V4 R1 now ends at calendar 2025 (released 2026-05-22; no 2026 before ~spring
  2027), so WY 1980–2025 climate is achievable — matching the products. The 6,087-basin
  input predates the 7,964-polygon layer (basin size explains ≤ 58 of the 2,049 gages
  without Daymet); a recompute over all 7,964 polygons is the only way to close the hole
  (no published product substitutes). ORNL THREDDS/NCSS/tiles are gone and pydaymet /
  daymetr / climateR are broken for gridded pulls; viable routes are the annual NA
  mosaics (ORNL HTTPS/S3; NCAR GDEX mirror through 2024 via Globus) with our own
  exactextract zonal statistics, gdptools by the USGS co-author, or Google Earth Engine
  (through 2025; licensing question for a private-company author). Recommended: prcp +
  swe first (≈ 0.57 TB), validated by reproducing the co-authors' 2023 values on the
  6,087 shared basins; four ≤ 1-day feasibility tests (T0–T4) and a go/no-go by
  ~2026-10-10 are laid out in `docs/plans/2026-09-29-daymet-reprocessing-options.md`.
  This rerun would also regenerate `flagged_for_high_na` (below).
  Measured the same day (plan §3 A-addendum): files are gzip-4 NetCDF-4 chunked
  (1,1000,1000) for 1980–2019 and (10,300,300) for 2020+; six variables total 3.38 TB
  for 1980–2024 (72–77 GB/yr), prcp+swe 0.564 TB; the NCAR mirror delivers 18–25 MB/s to
  the Windows laptop, decompression runs 150–170 MB/s per core → a download-bound
  year-by-year stream needs ≈ 2–3 GB RAM, ≈ 155 GB (six) / 25 GB (prcp+swe) of SSD, and
  ≈ 2–2.5 days (six) / ≈ 9 h (prcp+swe) of wall-clock.
  **Action plan for the dedicated machine (tools to write, runbook, acceptance criteria,
  twelve named unknowns): `docs/plans/2026-09-29-daymet-reprocessing-action-plan.md`.**
- **HydroShare documentation updates PENDING (not applied 2026-09-10 — user decision to
  leave the staged files untouched this session).** (H1) Category terminology: regenerate
  `hisss_data_dictionary.csv` `category` + `hisss_signature_categories.csv` (root/R1/R2)
  from `docs/signature_categories.csv` (8 categories incl. the scalars; today they carry
  the 9-class grouping); R1/R2 README lede (14-family list → the eight) and the file-table
  row ("nine-class exploratory grouping … separate taxonomy"); rebuild the explorer HTML
  with the updated builder; relabel or annotate the validation summary/dashboard groupings
  as function families. (H2) The `flagged_for_high_na` known-issue text in the R1/R2
  READMEs + dictionary row, rewritten when the column is regenerated (already planned
  above). The user noted a further aspect of the HydroShare docs also needs updating —
  record it here when specified. (H3, added 2026-09-10 evening) R3 README + input
  dictionary: state the HYDAT release used (2025-10-14; retrieval 2026-02-07) and that
  Canadian records end 2024-12-31, so Canadian gages have no WY 2025 in either product;
  R1/R2 READMEs could carry the same one-liner. Details: `docs/plans/2026-09-10-manuscript-category-edits.md` §C.
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
- **LOW (both delivered products + staged Resource 3, found 2026-09-29) — five Daymet sites
  are NaN on every day of 1980–2023, all six variables**: 02234500, 02236000, 02236125,
  02244040 (St. Johns River, FL) and 01372058. The co-authors' aggregation returned NaN for
  any basin touching a single fill cell (these basins are 99.95–99.98 % valid). The four
  Florida gages are in product #1 (all four) and #2 (three) with every climate signature
  NaN, so "gages with climate" counts overstate by 4 / 3. The reprocessed input
  (2026-09-30, masked mean) has complete series for the four (01372058 has no polygon and
  is in neither product); the products change only when rerun on it — until then note it
  wherever those counts are quoted.
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

### Guidelines Document TODOs and Manuscript Reconciliation Log
Moved 2026-09-29 to `docs/reconciliation/guidelines_todos.md` and
`docs/reconciliation/manuscript_log.md` (dated entries, newest first, maintained by the
`/sync-docs` skill). This file no longer carries them.

---

## [September 2026]

### Added: reprocessed Daymet climate input, 8,017 basins, 1980–2025 (run 2026-09-29/30)
`daymet_1980_2025_29sep2026.parquet` (run folder `~/HISSS_data/daymet-processed-29sep2026/`,
md5-verified copy `/Volumes/Untitled/daymet-processed-29sep2026/`): 134,605,430 rows
(8,017 × 46 × 365), no NaN in any variable, md5 `c059076088cd8abd963c64094468e90f`.
NOT yet the input of any delivered product.
- Against the stale product 1980–2023 (5,969 shared basins, 262,460 basin-years per
  variable): R² ≥ 0.9999999 in every basin-year for prcp, tmin, tmax, vp and srad; swe
  ≥ 0.999 in all but 2 trace-snow basin-years; prcp annual totals within ±0.0011 %
  (p01–p99); no date offset. Every gage of both products is covered (#1 6,678, #2 6,250;
  the stale file had 5,513 / 5,635 non-NaN).
- Run: 25.9 h of downloads (3.379 TB from ORNL) — 6.8 h over the estimate, lost to slow
  single connections that curl's stall floor does not catch — then 35 min to assemble,
  validate and copy. Tables: the run folder's `RUN_NOTES.md`; `EO_data_processing/README_DAYMET.md`;
  action plan §0.

### Added: Daymet reprocessing toolchain; single-year gate on calendar 2023 PASSED (2026-09-29)
**DECISION (user, 2026-09-29)**: reprocess Daymet for every watershed that has a polygon
(the 7,964-basin layer), entirely from scratch — no append to the stale series; check the
result against the stale co-author product; process, time and compare ONE year before any
multi-year run.
- New tools in `EO_data_processing/daymet/` (moved there from `docs/benchmarks/daymet/` the
  same day; user decision): `daymet_probe`, `_weights`, `_aggregate`, `_crosscheck`,
  `_pixelcheck`, `_validate`, `_assemble`, `_stream`, `copy_verify`; product README
  `EO_data_processing/README_DAYMET.md`. Method: exactextract
  coverage × true cell area on the unresampled Daymet LCC grid, chunk-aligned accumulation,
  mean over the cells valid that day; every raw file's SHA-256 checked against NASA CMR.
- 2023, six variables (72.2 GB from the NCAR mirror, byte-identical to ORNL): download
  66 min at 15–21 MB/s; aggregation 13–21 s per variable-year on the 16 GB M5 laptop (32 s
  for the 1980–2019 chunk layout) → 1980–2025 is download-bound, ≈ 50 h (prcp + swe
  ≈ 8.5 h). Exactness: exactextract's own weighted mean agrees to ≤ 5e-11 in all basins;
  ORNL's Single Pixel API agrees at six points; output is bit-identical across reruns.
  Probe on prcp 1980 (old layout, leap year): all chunks stored, Feb 29 kept / Dec 31
  dropped as in the stale data, stale values reproduced.
- Against the stale product (5,969 shared basins): with true-area weights and
  full-resolution polygons, prcp reproduces it (annual totals within ±0.0011 %, p01–p99).
  Coverage-only weights shift large northern basins by up to 0.3 %; the published
  200 m-simplified Resource 4 polygons leave ≈ ±0.1 % (prcp) and ±2.5 % (swe) at p01–p99.
  Found the five all-NaN stale sites (Known Issues).
- Results: `docs/plans/2026-09-29-daymet-reprocessing-action-plan.md` §0. Run folder:
  `/Volumes/Untitled/daymet_processed_sep2026/` (`RUN_NOTES.md`).
- Comparison dashboard: `EO_data_processing/viz/build_daymet_comparison_dashboard.py` (+ its
  HTML template) — one year, one or more polygon versions: summary tiles, original-vs-fresh
  scatter, difference map, daily series for curated basins, least-agreeing table. The 2023
  page is `daymet_2023_comparison_dashboard.html` in the gate run folder.

**DECISIONS (user, 2026-09-29, after the gate)** for the 1980–2025 run: (D1) full-resolution
polygons; (D2) all six variables in one pass; (D3) the user's Earthdata Login bearer token
for the ORNL-only 2025 files (expires ≈ 2026-10-25; kept outside the repo, docs and memory);
(D4) add the 53 correctly delineated basins > 100,000 km² if RAM allows (05KH009 stays out,
wrong polygon). **LOCATIONS (user, 2026-09-29):** code in `EO_data_processing/`; run folder
on the internal SSD `~/HISSS_data/daymet-processed-29sep2026/` with md5-verified copies to
the exFAT drive; commit and merge when ready. Implemented the same day: the polygon rebuild
script is now committed (`EO_data_processing/geometry/rebuild_watershed_polygons.py`; defaults
reproduce Resource 4, `--no-simplify --include-large` gives the 8,017-basin layer); bearer-token
downloads with the ORNL-only 2025 files first; an aggregator planning fix (int32 keys, slices;
byte-identical output) that keeps the 8,017-basin runs at 4.3 GB parent + < 1 GB per worker
with `--workers 6`; downloads prefer ORNL (48–54 MB/s here vs ~18 MB/s from the mirror), so
the 1980–2025 run takes ≈ 19 h instead of ≈ 50 h. The run started 2026-09-29 20:23 UTC.

### Changed: Claude Code instruction files restructured for the context budget (2026-09-29)
The eight files auto-loaded at every session start (CLAUDE.md plus seven `@`-imports:
DEVELOPMENT.md, SIGNATURES.md, CHANGELOG.md, both Google-Doc snapshots, both EO READMEs)
totalled ~354 KB ≈ 88k tokens. Official guidance: keep CLAUDE.md under 200 lines, `@`-imports
load eagerly, procedures belong in skills, directory-scoped constraints in path-scoped rules.
Restructured per the user-approved plan the same day:
- `CLAUDE.md` rewritten to ~100 lines of always-true rules; its only import is the new
  `docs/STATUS.md` (one-liners: products, pending decisions, live issues, deferred fixes,
  sync dates; ≤ 60 lines). Startup context is now ≈ 5k tokens.
- `.claude/` is TRACKED (only `settings.local.json` ignored) and excluded from the HISSS
  mirror. Path-scoped rules: `signatures-code.md`, `benchmarks.md`, `changelog.md`. Skills
  (description at startup, body on demand): `/sync-docs` with `sync_google_docs.py` — a
  stdlib fetch + paragraph diff of both Google Docs, verified to reproduce today's snapshots
  exactly, so an unchanged sync costs no context — `/add-signature`, `/run-benchmark`, and
  `/cross-language-alignment` (moved from `claude-skill/`). `claude-skill/streamflow-signatures.md`
  stays as the user-facing interpretation skill.
- `[Unreleased] → Guidelines Document TODOs` and `→ Manuscript Reconciliation Log` (461 lines)
  moved verbatim to `docs/reconciliation/`; CHANGELOG.md itself is no longer auto-loaded.
- `EO_data_processing/CLAUDE.md` (nested, loads only when working there) carries the id-join
  and S3-loss rules; the stacked dated status banners of both EO READMEs moved verbatim to a
  "Build history" section at the end of each file.
- Nothing was deleted; `docs/DEVELOPMENT.md` and `docs/SIGNATURES.md` are unchanged and are
  read on demand via the pointers in CLAUDE.md and the rules.

### Changed: signature categories are the co-authors' EIGHT; repo docs aligned (2026-09-10)
**DECISION (user, 2026-09-10)**: the colleague's 8-category sheet is the reference
grouping — Flow Volume, Flow Duration, Storage, Flashiness, Drought, Flow Timing,
Precipitation Streamflow, Snow. Vocabulary from here on: *category* = one of the eight
(manuscript, guidelines doc, HydroShare dictionary); *function family* = the computing
function (15 functions of `calculate_all_signatures()`), a sub-level. The 21 scalars
inherit their function's category (`season_excluded_years_*` → Flow Volume;
`ice_affected_days_total` stays outside as a preprocessing diagnostic). TQmean → Flow
Volume and avg_storage → Storage follow the sheet.
- **New canonical mapping** `docs/signature_categories.csv` (121 rows: signature, kind,
  category, function, module, function_family, category_source) — the source for
  regenerating the HydroShare dictionary/categories CSV later.
- **Repo docs aligned (applied)**: README.md "Signature Categories" table rebuilt as
  8 categories × function families (+ line 237 pointer); docs/SIGNATURES.md gained a
  "Signature categories" map in the Overview, its Summary Table carries a Category column,
  and the Pettitt signal table is labeled by function family; CLAUDE.md pointer;
  claude-skill overview lists the eight; `docs/plans/dataset_workflow_schematic.md` +
  PNG re-rendered ("in 8 categories"; playwright + chromium installed into `.venv`,
  mermaid byproducts git-ignored); CROSS_LANGUAGE_STATUS "Per-Category Results" note
  (those tables are by function family). `build_signature_explorer.py` `category_of`
  now reads the canonical CSV (rule-based fallback kept), so the next explorer build
  shows the eight; the comparison scripts' function-family groupings are unchanged
  (diagnostics).
- **Not editable from here — catalogued** in
  `docs/plans/2026-09-10-manuscript-category-edits.md`: 8 manuscript locations (§2
  preamble, §2.2.1, §3, §4 figure, §5.1.1, §5.1.3, §1) and 8 residual guidelines-doc
  edits (TQmean to 3.1, add the avg_storage block to 3.3, "(see 3.12)" → 3.1,
  ice_affected_days_total, two column names, module/family wording, Part 5 high_na).
- **HydroShare documents NOT changed this session (user decision)** — pending list in
  the same file (§C) and under Planned.

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

# Changelog

All notable changes to the Streamflow Signatures project.

This file holds only the CURRENT state: `[Unreleased]` (open plans, live known issues),
the most recent month in full, and condensed
summaries of the two months before it. Everything older — the full September, August and July 2026
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
- **Daymet climate input: the product-rerun decision is PENDING (NEXT).** The reprocessed
  input `daymet_1980_2025_29sep2026.parquet` (8,017 basins, calendar 1980–2025, no NaN) was
  built 2026-09-30 and adversarially reviewed 2026-10-01. On the 5,965 basins with stale
  data it reproduces the co-authors' series. The user accepted it on raw-data equivalence
  and dropped the signature replay (DECISION, 2026-10-01); the questionable polygons are
  kept and flagged in `daymet_basin_flags.csv`. No product uses it yet. A rerun on it would
  give 1,165 (#1) / 615 (#2) more gages usable climate, complete the four all-NaN Florida
  gages (Known Issues) and regenerate `flagged_for_high_na` (below). On the 16 GB laptop,
  run Julia in gage batches with `docs/benchmarks/run_batched_julia.py`. Record:
  `[October 2026]`, `[September 2026]` and the action plan §0
  (`docs/plans/2026-09-29-daymet-reprocessing-action-plan.md`); the superseded planning
  text of this entry is in `changelog-old.md`.
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
  wherever those counts are quoted (usable: 5,513 of #1's gages, 5,635 of #2's).
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

## [October 2026]

### Fixed (HIGH): second review of the session's Daymet work — batched runner, explorer, tools (2026-10-05)
Four Sonnet 5.5 reviewers (2026-10-01) re-checked every update since 2026-09-29: the Daymet
tools and flags, the batched Julia runner, the record explorer, and the docs. They found **no
wrong value** in the climate file, the flags table or the explorer's numbers, all recomputed
independently. Fixed:
- **Batched runner (HIGH for a rerun; no product was affected).**
  - `merge` took whatever batch folders existed, so a missing or failed batch gave an
    incomplete product with exit 0. It now requires every batch of `batches.json`: done, one
    CSV and timing JSON each, and an annual parquet in all or none.
  - It also checks each CSV's gages against the batch assignment, and the CSV and annual row
    counts against the batch's timing JSON. It refuses batches from different code, config or
    inputs.
  - Rows follow `--reference`; a gage on one side only is an error unless `--allow-missing`.
    Without a reference they are sorted by gage id.
  - `split` records each source's sha256, size and PAR1 footer, which the merged provenance
    reports.
  - `run` checks that every input exists. Julia runs in its own process group, killed on low
    memory, a failed memory sample, `--max-minutes` or SIGTERM.
  - The annual merge sorts in gage chunks: peak 1.3 GB instead of 2.8 GB.
  - Re-merging the 2026-10-01 control batches reproduces its CSV and annual parquet byte for
    byte. New `selftest_run_batched_julia.py` (28 checks).
- **Record explorer.**
  - Holding ↑ could hang the tab. Value zoom now stops at 10 quanta, and the tick loops are
    capped.
  - Steps like 2.5 were labelled with rounded values; 58,814 fuzzed ticks now parse back
    exactly.
  - "Exactly 0" was false for 1 − R², where R² rounds to 1.
  - The map legend states its 2nd–98th percentile clipping. Map zoom survives page-height
    changes, and tables no longer widen the page on phones.
  - Wording: "2,052 basins without an original series" (2,048 absent, 4 NaN), cells touched
    vs cells of area, lag counts of basin-years that have r, and annual-table digits.
  - Edge cases: no anomaly axis without overlap; day ticks for short windows; 14-day windows at
    the record end.
  - Accessibility: − / + value-zoom buttons (touch), screen-reader announcements, coverage
    shown by shape too, and the map no longer an empty tab stop.
  - Builder: NaN guards on the embedded series, and `<title>` moved into `<head>`.
- **Daymet tools.**
  - `copy_verify.py` writes `<file>.new`, flushes it with F_FULLFSYNC and replaces the drive
    file and its manifest line only after verifying. Before, a truncated source could replace
    a good copy. A file whose line records other content now needs `--replace`.
  - `daymet_assemble.py` verifies before replacing the output; before, a failed verification
    left an unverified file in place of a good one. It requires every `.done` and one weights
    md5, and stores NaN as NaN.
  - `--workers` defaults to 6 (8 needs ~16.7 GB). The NCAR mirror has its own 10 MB/s floor
    and is queried only when used.
  - `git_state` records a git failure as unknown, not clean.
  - `daymet_basin_flags.py` stops when any basin lacks a QA row. It still regenerates the
    delivered table byte for byte.
  - Tool self-test: 25 checks.
- **Docs.**
  - 2,052 basins are not compared (not 2,048).
  - 45 basins have less than 4 cells of area; only 5 touch fewer than 4.
  - 46 ECCC outlines also lack a metadata area.
  - The file has 38 % more rows, but 27 % more bytes.
  - 8,013 of the 8,017 ids match the streamflow parquet.
  - The memory caveat now covers the new input's batches.
  - The plans had stale replay, token-handling and precision statements.
- **Deferred (LOW):**
  - `daymet_stream.py`'s CMR re-query has no retry or offline path, and it does not notice a
    granule that vanished from CMR.
  - In dark mode the explorer's blue and orange are close in luminance, so zoomed-out edges fade
    in greyscale. Its group markers reuse the series colours.
- **Left to the user:**
  - CLAUDE.md and `/add-signature` say to publish the mirror after every merge, while STATUS
    records the 2026-10-01 deferral.
  - `EO_data_processing/CLAUDE.md` ships to the public mirror and mentions the token pasted on
    2026-09-29.
  - The STATUS line "37 Canadian gages … 27/37" describes the canonical run; the products hold
    32 / 28 such gages, of which 22 / 25 are flagged.

### Added: Daymet record explorer, basin flags, batched Julia runner; signature replay dropped (2026-10-01)
**Decisions (user, 2026-10-01):**
- Keep the questionable polygons and FLAG them. The flags travel as a companion table; the
  signature CSV's column contract is unchanged.
- Protect memory on the 16 GB laptop, which another session shares.
- Keep `main` current; the HISSS mirror will be force-pushed later.
- Later the same day: STOP the signature replay against product #1. In the user's words: "if
  the raw data are sufficiently equivalent then that is all we need to know".

What was built:
- **Basin flags.** `EO_data_processing/daymet/daymet_basin_flags.py` writes
  `daymet_basin_flags.csv` next to the climate file: 101 flagged basins.
  - 30 HydroBASINS fallbacks, 28 area mismatches (> 50 %), 27 HydroBASINS polygons with no
    metadata area (46 ECCC outlines lack one too and are not flagged), 45 with less than 4
    Daymet cells of area (coverage sum; wording corrected 2026-10-05).
  - 57 of them are low-confidence: 12 in product #1, 13 in #2.
- **Batched runner.** `docs/benchmarks/run_batched_julia.py` splits, runs and merges a
  Julia run in gage batches.
  - A memory guard stops Julia when available memory falls below `--min-avail-gb`.
  - `merge` refuses batches whose headers differ (the `flagged_for_high_na` denominator
    would differ) and orders the gages like a reference CSV.
  - The replay's control arm ran as 6 batches of ~1,336 gages: code b5f4c13 (clean
    tree), the stale climate input, product #1's config. Each batch took 121–177 s and
    peaked at 4.8–5.5 GB RSS.
  - Merged, the batches give 6,678 gages × 1,653 columns and 18,898,406 annual rows,
    product #1's shape. No batch log has a family failure.
  - The control arm was never compared with product #1, because the replay stopped
    first. Step (a) finished 1 of its 6 batches and step (b) never started (folder
    `~/HISSS_data/daymet-replay-01oct2026/`).
- **Record explorer.** `EO_data_processing/viz/build_daymet_record_explorer.py` writes
  `validation/daymet_record_explorer_2026-10-01.html` in the run folder. The page is
  10.9 MB, has a verified drive copy, and makes no network request. It shows:
  - the 1980–2023 agreement metrics as per-variable distributions over 262,460
    basin-years: largest daily |Δ|, RMSE, mean |Δ|, mean Δ, 1 − R² and annual totals, with
    quantile and lag tables;
  - a map of all 8,017 basins, coloured by coverage or by per-basin record RMSE / largest
    daily |Δ|;
  - original and new daily series overlaid on one chart for 30 embedded basins: 10 least
    matching, 10 random and 10 without an original series (seed 20261001). A switch adds the
    anomaly, new − original, on its own right-hand axis.
  - Later the same day, at the user's request: the new line thins as the view zooms out.
    From 1.5 days per pixel, the new draws only the edges of its min-max band, so the
    original's filled band shows through. The chart also zooms on values as well as dates:
    drag a band or a box, drag along either axis, or press ↑ ↓.
  The builder re-derives every embedded basin-year's largest |Δ|, which equals the
  validation table. It also decodes every embedded block back against the source series.

### Fixed (MEDIUM): Daymet input — adversarial review, doc corrections, toolchain hardened (2026-10-01)
Three independent reviews of the 2026-09-30 file (ingest; processing and storage; honesty of
the comparison with the stale product) re-ran checks read-only and found **no wrong value**
in it. Evidence:
- ORNL's Single Pixel API reproduces the basin means of small basins to float32 rounding:
  all six variables, both chunk layouts, leap years, and 2024–2025. The archived check is
  `daymet_outputcheck.py` (8 basins × 8 years, ≤ 1.1e-4).
- All 807,632,580 values of the parquet equal the per-variable-year files.
- All 276 source granules still match NASA CMR.
- The polygon layer rebuilds byte-identically.

What they found, and what changed:
- **Docs corrected (overstated or wrong claims; no value affected).**
  - README_DAYMET said the file already "replaces" the stale input. It is the input of
    neither product.
  - The comparison covers the 5,965 basins with stale data, not 5,969 (5,969 × 44 ≠
    262,460).
  - A rerun gives usable climate to 1,165 / 615 product gages: 1,161 / 612 absent from the
    stale file plus the 4 / 3 all-NaN ones. The earlier wording was "1,161 / 612 that had
    none".
  - The STATUS line had dropped its qualifiers (swe; p01–p99).
  - The gate's R² floor was 0.99999995, not 0.99999996.
  - The mirror's "byte-identical" rests on sizes for 263 of its files.
  - The CRS was not asserted per file.
  - Bit-reproducibility holds only for a fixed weights file.
  - The crosscheck covered the 7,964 gate basins, not every basin.
  - Added what was not compared at all: 2,052 basins (2,048 absent from the stale file,
    including all 53 large and all 30 HydroBASINS polygons, and 4 NaN there; count corrected
    2026-10-05); 2024–2025; the 7 basins that touch fill cells.
  - Added that agreement shows reproduction of the co-authors' aggregation of the same
    Daymet cells, not accuracy.
- **Toolchain hardened (the existing file is unaffected).**
  - `daymet_common.read_grid` asserts the CRS parameters, the last cell centre and the
    calendar.
  - `daymet_stream.py`:
    - hands the Earthdata token to curl on stdin;
    - terminates curl and joins the downloader on every exit (a fatal error used to leave
      a token header file and an orphan curl);
    - ends a source at once on HTTP 401/403/404;
    - resumes connections below 20 MB/s at once (the 200 KB/s floor lost 6.8 h);
    - re-queries CMR on every start, stops on a changed record, and records granule
      concept and revision ids;
    - refuses `.done` files built on other weights;
    - guards disk space and renames a checksum-failing file to `.bad`.
  - `daymet_aggregate.py` records the verified source SHA-256, the weights md5 and the
    commit, and an out-of-range value is fatal.
  - `daymet_assemble.py` verifies the written file value by value, documents its (year,
    site_id, Date) order and writes a self-contained provenance sidecar; `--provenance-only`
    upgrades an existing run's sidecar.
  - New `daymet_outputcheck.py`: a Single Pixel API check of the output that needs no raw
    files.
  - `copy_verify.py` replaces manifest lines instead of appending.
  - The polygon builder writes input hashes.
  - Transitive packages are pinned.
  - New `selftest_daymet_tools.py` (synthetic file plus local HTTP server): 23/23 checks.
- **Run folder updated; the data are unchanged.**
  - The sidecar now carries the verification and the run's provenance.
  - Added `validation/outputcheck_single_pixel_2026-10-01.json`,
    `validation/cmr_requery_2026-10-01.json` and
    `polygons/…provenance_rebuild_2026-10-01.json`.
  - RUN_NOTES corrected; drive copy re-verified.
- **Open before a rerun** (action plan, review block; the first three items were settled
  later the same day, see the entry above):
  - The Julia runner loads all eight columns (~8.6 GB for this file, extrapolated); give it
    a 4-column copy on the 16 GB laptop.
  - Decide on the 57 low-confidence polygons (12 in #1, 13 in #2; 27 HydroBASINS polygons
    cannot be area-checked).
  - Replay in two steps: shared basins 1980–2023, then the full file.
- **Security (user action):** the Earthdata token pasted into chat on 2026-09-29 persists in
  the local Claude transcripts. Revoke it at Earthdata Login; the run no longer needs it.

## [September 2026]

Condensed summary — **the full text of every entry is in [changelog-old.md](changelog-old.md) → [September 2026]**.
The month rebuilt the climate input from the Daymet mosaics, restructured the Claude Code
instruction files, adopted the co-authors' eight signature categories and gave
`flagged_for_high_na` one definition in all three languages.

- **2026-09-30 — reprocessed Daymet climate input BUILT**: `daymet_1980_2025_29sep2026.parquet`
  (8,017 basins, calendar 1980–2025, no NaN). It reproduces the stale co-author series on
  the 5,965 basins with stale data (1980–2023) and is not yet the input of any product.
  Reviewed 2026-10-01 (`[October 2026]`).
- **2026-09-29 — Daymet reprocessing toolchain** (`EO_data_processing/daymet/`) and the
  calendar-2023 gate. **DECISIONS (user)**: every polygon, from scratch; full-resolution
  polygons; all six variables; the user's Earthdata token for 2025; the 53 basins
  > 100,000 km²; code in `EO_data_processing/`, run on the SSD with md5-verified drive
  copies.
- **2026-09-29 — Claude Code instruction files restructured**: CLAUDE.md ~100 lines +
  `docs/STATUS.md`, path-scoped rules, skills including `/sync-docs`, reconciliation logs
  moved to `docs/reconciliation/`.
- **2026-09-10 — signature categories are the co-authors' EIGHT (DECISION, user)**. The
  canonical map is `docs/signature_categories.csv`; repo docs aligned. The manuscript,
  guidelines and HydroShare edits are catalogued in
  `docs/plans/2026-09-10-manuscript-category-edits.md`.
- **2026-09-04 — Fixed (HIGH): `flagged_for_high_na` has one definition in all three
  languages**: signature columns only (1,224 → 791 on product #1). rpkg's constant-series
  Mann-Kendall was fixed too. The delivered products are NOT rewritten (DECISION, user).

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

## [July 2026] – [March 2026]

Condensed summaries of these months — the Julia annual-values export, b = 1 recession
alpha, snow and drought families, 20-value stats floor, Q-to-PPT unit gate and 60 % trend
gate, both standard products and the Annual NLCD product (July);
the MODIS LAI/LULC EO products (June); the
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

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
- **Climate on Daymet's own dates — the fix for the streamflow-frame Known Issue (found
  2026-10-16) — ON HOLD pending discussion with the research team (user decision,
  2026-10-16).** Scope when approved: outer join in the three runners, Julia first; gate:
  every Q-only column bit-identical (`check_additivity.jl --expect-added 0`); then rerun
  both version-2 products (the 1,172 gage-years with absent streamflow dates change in the
  PPT/TMIN/TMAX/SWE-derived families). Design: `docs/plans/2026-10-16-precipitation-review.md` §5.
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
- **Compiled USGS records drop every day whose NWIS qualifier is `A R` (approved, "records
  revised") — HIGH, every product (found 2026-10-21).** `generate_streamflow_dt`
  (`R/helperFunctions.R`) keeps a value only when the qualifier string is exactly `A`, `A e`,
  `P` or `P e`; `R` = "Records for these data have been revised" marks an APPROVED value. In
  `combined_streamflow_data_09feb2026.parquet`: 181,461 `A R` days on 755 gages are NA, plus
  the other approved remark codes `A <` (11,444 days / 139 gages), `A >` (218 / 30) and `A [4]`
  (419 / 101) — 193,542 approved days on 987 gages (908 in product #1 v2, 757 in #2 v2).
  Full-calendar replay of the preprocessor's rules (2026-10-23, Codex-verified): 1,370 water
  years in 676 gages become valid once those days are finite (1,197 / 601 from `A R` alone);
  1,177 of them in #1 v2 (619 gages) and 950 in #2 v2 (502 gages); 899 / 757 product gages hold
  such a day inside their window; revisions cluster in 2013–2025, so record ends and last-decade
  completeness are hit hardest; 5 candidates fell below the 20-year compilation floor because of
  it. The provisional remark codes (`P Ice`, `P Dis`, `P Rat`, `P ***`, `P Bkw`, `P Eqp`, `P ZFL`,
  …; 15,257 days on 263 gages) are masked too; many of those days carry a value at the agency
  today under a new status, their February state is unknowable. Fix = the rule in `keep_codes`,
  a re-fetch and splice, a rerun of both products with the gage sets FROZEN (user decision
  2026-10-23) — NOT applied; plan v3 `docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md`,
  **DECISIONS (user, 2026-10-26): every default of plan v3 approved** (restore `R`, `<`, `>`,
  `[4]` with the censored bounds kept as values; provisional-remark days out of scope; splice as a
  patched snapshot with ledger, raw responses and the three-dataset control; gage sets frozen;
  the climate outer-join fix stays out; internal v3 / public "version 2"; metadata updated in
  place; rpkg re-baseline; code reviews by area; local staging; start now) — implementation
  pending. Measurement: census §10.1
  (`docs/plans/2026-09-04-filtering-stage-census.md`); tools
  `docs/benchmarks/min_flow_rule_footprint.py`, `min_flow_rule_replay.py`.
- **MEDIUM (every product, found 2026-10-16 by the Pacific-coast precipitation review) —
  climate rides on the streamflow DATES: a day absent from a gage's streamflow record (no
  row, as opposed to a row with NA) loses its Daymet value as well, and the preprocessor
  then linearly interpolates PPT (and TMIN/TMAX/SWE) across the hole as if the climate were
  missing.** All three runners left-join each gage's climate onto its own streamflow rows
  before `preprocess_daily_data()` builds the date grid
  (`docs/benchmarks/run_julia_benchmark.jl:220` and `:332`, `run_python_benchmark.py:319`,
  `run_rpkg_benchmark.R:276`), so the library never sees the Daymet values of those days.
  Gages drop out preferentially during storms, so the interpolate is biased low. Measured on
  the WY 1993–2025 version-2 annual parquet against the Daymet input: 919 of 200,832
  `prcp_annual` gage-years (679 gages) differ from the true water-year sum, every one with
  ≥ 1 absent streamflow date (0 counterexamples; 253 further gage-years with absent dates
  happen to match); |diff| median 6 mm, p90 54, p99 327, max −606 mm (02146285 WY 2024: 22
  absent dates in 1–2-day runs carried 713 mm, 32 mm/day against 2 mm/day on the present
  days; product 943 mm, Daymet 1,549 mm). Per-gage `prcp_annual_mean` bias exceeds 1 % at 8
  gages (02146285 −24.5 %, 0214297160 −6.1 %, 0214269560 −3.4 %, 02329342 −3.3 %, 02336030
  −2.4 %, 11363930 −1.9 %, 11153650 −1.8 %, 03250190 −1.1 %); the 70 affected Pacific-coast
  gages are all within 1.9 %, so this does NOT explain the coastal 3–4.5 m totals (which the
  review found consistent with or below PRISM basin means and station normals). The same
  frame feeds the Q-to-PPT family (runoff ratios, elasticity, Q-P seasonality, avg_storage),
  the snow family (`swe_max_to_ppt`, SWE itself) and the temperature metrics, so by code
  reading the 2026 deliveries and the WY 1980–2025 version-2 product carry the same
  behaviour (not measured there). Not a library defect: `preprocess_daily_data()` is correct
  for the frame it is given. **Proposed fix**: give the preprocessor the climate on Daymet's
  own dates — an outer join on `date` in the three runners (Q becomes NA on the added rows,
  which the grid normalisation already produces for absent dates, so every Q-only column
  must stay bit-identical; `check_additivity.jl` is the gate) — Julia first, then the
  ports, then a rerun of both version-2 products. **DECISION (user, 2026-10-16): documented
  for future work; the fix is NOT implemented until it has been discussed with the research
  team.** Record and evidence: `docs/plans/2026-10-16-precipitation-review.md`;
  `~/HISSS_data/climate-rerun-12oct2026/precip_review_16oct2026/` (scripts +
  `full_domain_prcp_check.parquet`, one row per product gage-year).
- **Explorers and dashboards built before 2026-10-13 show CARTO's "API KEY REQUIRED"
  watermark instead of a basemap** (CARTO withdrew keyless tiles; points and panels still
  work). Affected deliveries: the 2026 explorers and validation dashboards of HydroShare
  R1 / R2 (superseded by the version-2 deposit) and the three EO explorers of R5
  (`hisss_modis_lai_explorer.html`, `hisss_modis_lulc_explorer.html`,
  `hisss_nlcd_explorer.html`) — rebuilding and redepositing R5's explorers is a user
  decision. Builders fixed 2026-10-13 (October 2026 → Fixed).
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

### Reviewed: what the compilation's flow rule removed; the `A R` qualifier mask found (2026-10-21)

User question: how many of the 16,994 candidate gages did "at least 20 years with > 30 days
over 0.0001 mm/day" remove, and how many because of the flow threshold itself? Re-fetched the
80 `insufficient_years` candidates (NWIS + ECCC GeoMet) and replayed the R compilation gate
under four rule variants (76 of 80 recorded year counts reproduce exactly), and scanned the
8,014 compiled gages. Result: the threshold proper removed **72** gages — 67 zero-flow /
diverted channels (the threshold acts as "> 0") and 5 at the 0.0001 level; 8,787 of the 8,980
excluded failed the 7,300-row length pre-check before any value was examined, 113 were never
retrieved, 4 were lost to the qualifier mask, 3 to length, 1 now passes. Among the compiled
gages the rule left 2,269 gage-years (0.7 %) uncounted, none of which affected admission.
Side finding, HIGH → Known Issues above: the `A R` (approved, revised) qualifier mask. Full
record `docs/plans/2026-09-04-filtering-stage-census.md` §10; per-gage table
`docs/filtering_min_flow_rule_2026-10-21.csv`; tools `docs/benchmarks/min_flow_rule_*.py`;
the manuscript's §2 funnel sentence wording is relayed in `manuscript_log.md` (2026-10-21).
Follow-up (user request 2026-10-21, revised 2026-10-23): a plan to fix the mask and regenerate
both products with a QA/QC against the version-2 runs —
`docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md` (version 2). **DECISION (user,
2026-10-23): no new watershed enters the analysis** — the corrected products keep the version-2
gage sets; gages that would newly qualify (29 / 26) stay out and are listed. Codex round 1
(`docs/plans/2026-10-21-qualifier-reviews/codex_plan_review_round1.md`) returned NO-GO on the
first draft: its oracle counted NAs among existing rows rather than on the full water-year
calendar (50 false positives, 33 misses — rebuilt: 1,370 restorable water years in 676 gages for
the four approved codes, 1,197 / 601 for `A R` alone; 899 / 757 product gages hold a masked
approved day inside their window); the claim that provisional-remark days carry no value today
was wrong (1,273 numeric in the samples, re-coded after approval); restoring Q does not change
the climate join (the rows exist); the harness, validator, rpkg gate and staging script carry
hard-coded assumptions that must be parameterized first; the schedule is four to six days, not
two to three. All 14 findings verified and folded in; the stale "15,460 provisional days" above
and in census §10.1 corrected to 15,257; `docs/filtering_min_flow_rule_2026-10-21.csv` now
marks 11129800 as a mask casualty too. Year-level oracle (forecast):
`docs/plans/2026-10-21-qualifier-fix-expected-years.csv`. Codex round 2 (delta review of plan
v2, `codex_plan_review_round2.md`): GO WITH CHANGES — the reviewer re-derived all 1,656 oracle
rows, the totals and the would-qualify lists, and confirmed by a Julia experiment that
climate-only annual values stay exactly equal on shared years; its eight changes (the runner
derives `num_water_years` / `start` / `end` from the qualifying years, so "metadata unchanged"
must exclude them; a deposited per-product selection manifest with a `--gage-list` in the
splitter and an auditor that distinguishes qualified / selected / emitted; an acceptance oracle
built from base + cached source records; an audit-fetch universe of all 1,185 masked gages; a
numerical sensitivity run for censored values; a full resume fingerprint and a dirty-tree
refusal for products; a full baseline reproduction of v2; partitioned validators; about seven
person-days) are folded into plan v3 (2026-10-26). **DECISIONS (user, 2026-10-26): all plan-v3
defaults approved** (recorded in the Known Issues bullet above and in plan §3 / §9). **Later the
same day (user, after discussion with the research team): decision 5 REVERSED — the
climate-frame fix (Known Issues, MEDIUM, 2026-10-16) is folded into this run**; plan v4 adds
§10: a runner-level, span-restricted outer join (an unrestricted one would move every gage's
qualification denominator and drop gages, breaking the freeze), tests in the three languages,
a three-stage attribution (S0 baseline, S1 spliced input + old join, S2 = the product, plus a
one-batch interaction check), and a change census that must reconcile every changed cell to
one of the two input ledgers (0 unexplained changes, per-column magnitude and p-value-crossing
tables, the eight > 1 % precipitation-bias gages re-checked). Forecast
`docs/plans/2026-10-26-outer-join-forecast.csv`: 1,174 / 938 Q-valid gage-years with absent
days in 857 / 722 product gages, 147 / 122 of them also splice-affected; 5,069 / 4,893 gages
untouched by either fix. Codex round 3 on v4: NO-GO — a span-only outer join would insert
3,562 / 5,116 whole interior water years that have no streamflow (759 / 894 product gages),
which the preprocessor rejects but whose seasonal flags move the Q-only
`season_excluded_years_*` scalars; the inserted rows must be restricted to the gage's original
water-year set (fixtures pin it); the ledger must carry secondary interpolation changes (the
Daymet-missing Dec 31 moves when a neighbour is filled; 3 per window) and nonlinear bases are
gated by recomputation, not summed deltas; the restored P falls on 501 and rises on 418 of the
919 deviating gage-years, so no direction is gated; the stages run on one tree behind a
fingerprinted `STREAMFLOW_CLIMATE_FRAME` switch, stage 1 only on the splice-affected gages;
695 `_pval` columns of five types. All eight findings verified and folded into plan v5
(`codex_plan_review_round3.md`). Codex round 4 on v5: GO WITH CHANGES, decisions 12–17
endorsed, the Q-only invariance of the year-set restriction confirmed by the reviewer's own
Julia fixture, 7 findings, all verified and folded into plan v6
(`codex_plan_review_round4.md`): a ledger of changed PREPROCESSED values loses raw-mask changes
(a restored day equal to its interpolate flips a raw seasonal completeness flag, `io.jl`
l. 326, with no value delta) → three input ledgers plus per-stage state tables (raw masks,
interpolation mask, seasonal flags, rejection reasons); `flagged_for_high_na` is computed over
the EXECUTED gages' columns (`run_julia_benchmark.jl` l. 466), so every subset stage pins the
full-product header as its denominator; identities gated on executed gages only, the S1′
cohort composed from the realized support (165 / 135 overlap gages — 24 / 15 gages carry
restored years with absent dates, 18 / 13 of them outside the join forecast); the registry
enumerated per output kind (the record-wide recession minimum-event gate, `recession.jl`
l. 515, and the Q–P seasonality floor make OLD years' keys appear; the snow / frost anchored
gates touch the trend statistics only); the GDD gate uses the thresholded per-day formula
(`climate.jl` l. 476); §0 / Phase D / E4 / E5 aligned with §10 (A = complement of U;
`emitted = selected ⊆ qualified`); §10.8 narrowed (exposure ≠ residual); the effort table sums
to 12.5 person-days + ~23 h unattended, and one person cannot finish before the 2026-11-09
submission (decisions 18 baseline drift and 19 schedule proposed). Nothing run yet.

### Reviewed: Pacific-coast annual precipitation totals; one defect found (2026-10-16)

User question: basin totals of ~4 m on the Pacific coast looked high. Findings
(`docs/plans/2026-10-16-precipitation-review.md`): (1) only 3 of 6,678 gages in the
WY 1993–2025 version-2 product exceed 4,000 mm and 58 exceed 3,000 mm — Coast Mountains,
windward Olympics/Cascades and SE Alaska headwaters; the explorer's robust colour scale
clips at the domain p95 (1,764 mm), which paints every wet coastal basin the top colour.
(2) Against GAGES-II PRISM 1971–2000 basin means the product's median ratio is 0.99–1.03
west of the Cascades (1.07 for the rest of CONUS) and 0.59–0.87 in the wettest Olympic
basins; station normals (Forks 3,046, Hoh 3,298, Tofino 3,271, Coquitlam Lake 3,616,
Ketchikan 3,798 mm; Quinault basin estimate 4,724 mm) bracket or exceed our values; in
BC / SE Alaska 43 gages have mean runoff > mean precipitation (Daymet low in windward
mountains, plus glacier melt). (3) `daymet_aggregate.py` is a coverage × true-area weighted
MEAN, not a sum, validated against exactextract and the co-authors' gdptools series — no
defect. (4) `calculate_precipitation_metrics` reproduces the Daymet water-year sum to
1e-6 mm for 199,913 of 200,832 gage-years; the Dec 31 fill in 364-day Daymet water years is
the documented behaviour. (5) The remaining 919 gage-years expose the streamflow-frame
defect recorded under Known Issues (MEDIUM; fix deferred — user decision 2026-10-16).
Docs: caveats in `docs/SIGNATURES.md` §8 / §12 / §14, `docs/DEVELOPMENT.md` (NA handling),
the claude-skill and the co-author relay list; STATUS one-liner. No code changed.

### Fixed (MEDIUM): the staged version-2 explorers could not find their `_annual/` sidecar (2026-10-14)

Answering "which sidecar files must be uploaded with the HTMLs" showed that both staged
explorers in `~/Downloads/Signatures_v2_12oct2026/` embed the RUN folder's sidecar name
(`signature_explorer_1993_2025_60pct_climate_annual`): `build_signature_explorer.py` derives
the folder name from the output stem and writes it into the page (`D.meta.annual.dir`, the
only thing the page loads at view time besides Leaflet and the tiles), and
`stage_hydroshare_climate_v2.py` renamed the folder to
`hisss_signature_explorer_wy1993_2025_annual` without rewriting the page. The map and the
signature values (embedded) worked; every per-gage time-series panel would have said "annual
sidecar not found" — confirmed headlessly (a sidecar `<script>` load fails under the embedded
name, succeeds under the staged one). The 2026 staged explorers embed the deposit name and are
not affected; the run-folder explorers keep their own names and never were. Fix: the stage
script now rewrites the embedded name right after renaming the folder (`rename_sidecar_dir()`,
exactly one match or abort), and the two staged pages were patched in place (that one string;
now 100,193,900 and 95,036,409 bytes, md5 `45d257ee…` / `16877dda…`); both load `Qann.js` from
the embedded name. Deposit rule: upload each explorer HTML together with its complete
`_annual/` folder (139 `.js` files; 178 MB for WY 1993–2025, 226 MB for WY 1980–2025) at the
same level, names unchanged — the page also needs internet for Leaflet (unpkg) and the basemap
tiles. Rule added to `.claude/rules/benchmarks.md`.

### Fixed (MEDIUM): explorer and dashboard basemaps — CARTO's tiles now require an API key (2026-10-13)

The user's review of the version-2 explorers found no basemap, only CARTO's "API KEY
REQUIRED" watermark: `basemaps.cartocdn.com` now answers every keyless request with that
2 KB tile (HTTP 200, etag `wm-…`; CARTO's offer is a free key at 5 M tiles/month for
non-commercial use), so every HTML the twelve Leaflet builders ever produced lost its basemap
(`explorer_template.html`, the three EO explorers, `build_signature_explorer.py`, the seven
`build_*_dashboard.py`). All twelve now call one `hisssBasemap(map, variant, maxZoom,
position)`: Esri's keyless World Light / Dark Gray canvas (base + reference labels, native
zoom 16) first, an automatic switch to OpenStreetMap when three tiles fail before any loads,
and a Leaflet layers control for a manual switch (bottom-right in the signature explorer,
whose top-right corner is the gage panel; top-right elsewhere). Esri's legacy raster
basemaps are in mature support (no updates since 2021) with retirement phases announced for
March 2028 and December 2029 — the fallback and the switch exist for that day. Verified
headlessly with Playwright (new `docs/benchmarks/check_html_basemap.py`): light and dark
test pages load 60 / 60 Esri tiles, a page whose Esri host does not resolve switches to
OpenStreetMap (30 / 30 tiles, "OpenStreetMap" selected); the two version-2 explorers and
comparison dashboards were rebuilt (each explorer differs from its 08:54 build in the 34-line
basemap block only; 60 / 60 and 24 / 24 tiles loaded; the 139-file `_annual/` sidecars
byte-identical) and the explorers re-staged in
`~/Downloads/Signatures_v2_12oct2026/` (the `_annual/` sidecars unchanged; the staged READMEs
now say the explorer needs internet for Leaflet and the tiles). KNOWN ISSUE below: every HTML
built before this date shows the watermark — the deposited 2026 explorers and validation
dashboards of HydroShare R1 / R2 (replaced by the version-2 deposit) and the three EO
explorers of R5 (rebuild on request). Docs: `docs/DEVELOPMENT.md`, `docs/DATA_SOURCES.md`,
`.claude/rules/benchmarks.md`.

### Added: Climate metrics family — guidelines §3.9, 39 annual signatures (2026-10-12, Julia first)
**Trigger.** The guidelines doc gained §3.9 "Climate Metrics" (synced 2026-10-12,
`docs/reconciliation/guidelines_todos.md`) and the reprocessed Daymet input
(`daymet_1980_2025_29sep2026.parquet`, every gage, calendar 1980–2025) is ready. The user
asked (2026-10-12) for: the code, blind adversarial reviews, a rerun of both standard
products on the new input, output reviews, explorers and a staged deposit. Design record
and run plan: `docs/plans/2026-10-12-climate-metrics-plan.md`; relay list for the
co-authors: `docs/plans/2026-10-12-coauthor-relay.md`.

**DECISIONS (user, 2026-10-12)** — shipped column names stay (the doc's `_annual` /
`_xSeason` convention applies to new metrics only; doc-side revert relayed); new season
suffixes are the code's `win/spr/sum/fal`; `precip_` → `prcp_`, `GDD_` → `gdd_`; the §3.8
snowfall lines are deleted from the doc (not implemented); GDD = Σ max(0, (tmax+tmin)/2 −
0 °C); a dry day has prcp ≤ 0.1 mm; `prcp_annual_sd/_cv` are the within-year spread of the
daily values; ties take the first occurrence; the 5-day window stays inside the water year;
frost = basin-mean tmin < 0 °C, last frost within Oct 1 – Jun 30; Gini over all days (zeros
included), Walsh–Lawler index on the water year's 12 months; year qualification and gating
follow the snow-family precedent (climate years evaluated for streamflow-valid years;
record-anchored decade gate on the two frost dates); both products are rerun on the new
Daymet input as new versions (the delivered ones are not rewritten), Julia first, the
ports validated before publication; `flagged_for_high_na` regenerated and the rpkg
reference rerun fold in; HydroShare staged only.

**What changed (Julia canonical).**
- `julia/src/climate.jl` (new): `calculate_precipitation_metrics` — 20 bases on PPT-valid
  years (`prcp_annual`, `prcp_win/spr/sum/fal`, `prcp_1d_max`, `prcp_5d_max`,
  `prcp_max_dowy`, `prcp_50_dowy`, `prcp_cm_dowy`, `prcp_ndays_abv_10/25/50mm`,
  `prcp_annual_sd`, `prcp_annual_cv`, `prcp_intermit_annual` + `_dowy_start` / `_dowy_end`,
  `prcp_gini`, `prcp_seasonality_index`); `calculate_temperature_metrics` — 19 bases on an
  explicit `temp_data` frame of temperature-valid years (`temp_max_annual` + `_dowy`,
  `temp_min_annual` + `_dowy`, `temp_max_<s>`, `temp_min_<s>`, `frost_first_dowy`,
  `frost_last_dowy`, `gdd_cum_annual`, `gdd_cum_<s>`). Every base carries the 16 statistics;
  no per-gage scalar. **39 × 16 = 624 new columns → 2,277 per product; 139 annual
  signatures.** Category **Climate** (ninth), function families Precipitation and
  Temperature (`docs/signature_categories.csv`).
- `io.jl`: `read_parquet` maps `tmin`/`tmax` → `TMIN`/`TMAX`; `preprocess_daily_data` gains
  the temperature block (≤ 30 raw NaN, gaps ≤ 3 days filled, no sign check; both series must
  pass) and returns `valid_temp_years`. `signatures.jl`: `temp_data` kwarg; the precipitation
  family runs on the PPT-valid frame regardless of `area_normalized` (pure climate — the snow
  precedent; the gate test's "identical to no-climate" invariant now excludes `prcp_*`).
- Config: `climate` section + `na_handling.temperature_na_policy` (bundled copies
  byte-identical, sha256 a9e8ed31…). Absent `climate` section ⇒ no climate columns (the
  additivity switch). Column names derive from the config (`prcp_5d_max`, `prcp_50_dowy`,
  `prcp_ndays_abv_<t>mm`), the `D{p}_day` pattern.
- Runner, smoke test and `run_batched_julia.py` keep six climate columns (`+ tmin, tmax`);
  the runner records `climate_columns` in the timing JSON. Registries: `config.R`
  `EXPECTED_SIGNATURE_BASES`, `EXPECTED_DENSE_SIGNATURES`, `validate_production_run.py`
  (`--expect-signatures`, default 139), the four dashboards' `SIGNATURE_GROUPS` (snow and
  drought added to the three lists that still lacked them).
- Tests: `julia/test/test_climate_metrics.jl` — closed-form years (single-storm Gini
  (n−1)/n and SI 11/6, uniform year, threshold boundaries, tie rules, the non-contiguous fall
  season, leap-year grid, frost window incl. the Jun 30 cut-off, GDD clamp, the decade gate,
  the preprocessor's temperature validity, the opt-in orchestrator contract). Unit suite green
  after the gate-test update; smoke test on 10 real gages with the new Daymet file: all 39
  bases present, ranges plausible.
- **Same-machine check (2026-10-12):** the 2026-10-01 control run of product #1's
  configuration on this laptop reproduces product #1 with every Q-only column bitwise
  identical; the 98 columns that differ are all climate-dependent, ≤ 3.4e-13 (the rebuilt
  input). So the reruns' Q-only columns must match the delivered products exactly.

**Blind review round 1 (Codex, Julia; `docs/plans/2026-10-12-climate-reviews/`) —
NO-GO, seven findings, all resolved the same day:**
1. MAJOR `prcp_50_dowy` compared a `cumsum` against a differently accumulated `sum`
   (two equal 5 mm blocks separated by a dry interval returned day 301, not 50). FIX:
   every sum, mean, sd and cumulative sum in the family is a plain sequential double
   accumulation (`_seqsum` / `_seqcumsum`), so the cumulative sum's last element equals
   the total and the crossing test is `cum ≥ frac × total`; the same loop in Python (the
   builtin `sum` is compensated since 3.12 and must not be used) and R (whose `sum` /
   `cumsum` accumulate in long double) — which also makes the annual values bit-identical
   across the three languages.
2. MAJOR the precipitation family inherited the legacy fallback to `gage_data` when no
   `climate_data` is passed, bypassing PPT qualification. FIX: it runs ONLY on the
   explicit PPT-qualified frame (the snow contract); tests, smoke test and the gate test
   pass the frame explicitly.
3. MAJOR the complete-grid guard accepted a leap water year truncated to 365 rows. FIX:
   the guard takes the water year and requires its exact day count.
4. MINOR a missing `month` value silently dropped days from every season. FIX: months are
   derived from the date on the grid; the `month` column is no longer read.
5. MINOR a fractional percentage would be mislabelled (`cumulative_fraction` 0.504 →
   `prcp_50_dowy`). FIX: config validation requires a whole percentage.
6. MINOR window 1 duplicated `prcp_1d_max`, window > 366 was always NaN. FIX: window
   validated to 2..366.
7. NOTE ±Inf counted as a valid interpolation neighbour. FIX: non-finite temperatures are
   treated as missing.
Regression tests for 1, 3 and 4 in all three languages. Round-2 additivity on batch 1
after the fixes: PASS again (1,652 shared columns bitwise, +624 populated, 0 failures).

**Ports (same day):** `python/streamflow_signatures/climate.py` and `rpkg/R/climate.R`
module for module (loader, preprocessor temperature block, orchestrator, runners,
schema registries, tests). Cross-language audit tool `docs/benchmarks/crosscheck_climate_fixture.{jl,py,R}`:
a deterministic 30-year fixture exchanged as Parquet (CSV float parsers in pandas and R
are not round-trip exact); **1,170 annual values bit-identical in Julia, Python and R;
all 624 summary values within 4e-14 / 1e-12**, except rpkg's Mann-Kendall statistics
at ~6e-8 relative (the `Kendall` package's single-precision Fortran — the documented
pre-existing residual; tolerance 1e-6 for those four suffixes). Unit suites green in
all three languages (Julia 355 climate assertions, Python 173 tests, rpkg 81).

**Blind review round 1 of the ports (Codex, Python and rpkg) — both NO-GO; all findings
resolved the same day, several of them pre-existing port-wide residuals:**
- BLOCKER (both ports) interpolation arithmetic: pandas `interpolate` and R's
  `stats::approx` round differently from Julia's `left + (right − left) × offset / span`,
  so a filled day (Daymet's absent Dec 31 in every leap year) broke bit-identity and could
  flip a frost date. FIX: one explicit gap filler with Julia's expression order in each
  port, used for Q, PPT, SWE and TMIN/TMAX alike (the ports' existing families move
  toward Julia by last bits on interpolated days).
- BLOCKER (rpkg, library-wide) Spearman p-values came from `cor.test`'s exact / AS89
  method; Julia and Python use the t-approximation. FIX: `exact = FALSE`.
- BLOCKER (rpkg, library-wide) `Kendall::MannKendall` computes in single precision
  (distinct doubles collapsed into ties, small p-values rounded to 0, tau off by ~1e-7).
  FIX: `mann_kendall_test` is now a double-precision port of Julia's formula; the climate
  fixture's Mann-Kendall statistics agree to 4e-14 (were 6e-8). Both library-wide fixes
  land in the rpkg reference rerun already scheduled (user decision 2026-10-12).
- MAJOR (both) signed zeros: −0.0 / +0.0 extrema and ties differed by language. FIX:
  `+ 0.0` normalisation before the kernels in all three languages.
- MAJOR (Python) `limit = 0` crashed pandas when a temperature gap setting is 0 — gone
  with the explicit filler; negative policy values now rejected in all three.
- MAJOR (both) config coercion: `int()` / `as.integer()` truncated 5.9 → 5, string or
  numeric "booleans" were accepted, Python's JSON reader accepted NaN / Infinity. FIX:
  strict typed readers (integrality, finiteness, real booleans) for the climate and
  temperature keys; NaN / Infinity JSON tokens rejected.
- MAJOR (rpkg) the complete-grid guard truncated fractional dowy values; explicit
  temperature-policy overrides were ignored; the hard-coded DJF/MAM/JJA/SON ignored a
  configured `season_definitions`; the runner dropped temperature-only climate input.
  All fixed (exact numeric comparison; policy from the supplied config; definitions read
  in `.onLoad`; climate columns merged independently of PPT).
- MINOR (both) threshold labels diverged for scientific notation / large values and could
  collide. FIX: one shared formatter (integer, else C `%g`, `.`→`p`, `-`→`m`) plus a
  uniqueness check in all three.
- NOTE (all three) ±Inf passed the annual guards. FIX: the kernels require finite values.
Python 174 tests green; rpkg full suite green (the Mann-Kendall test now tolerates the
documented 1e-7 Fortran residual of the old backend). Three-way cross-check after the
fixes: annual values bit-identical, every summary value within 4e-14 / 1e-12 in both
ports.

**Blind review round 2 (Codex, all three languages; Sonnet, Julia) — Codex NO-GO on a
residual list, Sonnet GO. Resolved the same day:**
- MAJOR (shared) the 50 % crossing at a mathematical tie: two equal 7-day rainfall blocks
  gave sequential prefix 0.7 against half-total 0.7000000000000001, moving the date by
  294 days; a constant-0.1 mm leap year gave day 184 for 183. DECISION (implementer
  2026-10-12, confirmed by the user 2026-10-13): the crossing is
  `cumulative ≥ fraction × total − 1e-12 × total` in all three languages
  (`CLIMATE_CROSSING_REL_TOL`); 1e-12 of an annual total is physically nothing, and the
  closest real margin in the 360,585 production water years is 1.3e-9 (Sonnet's scan), so
  no production value depends on it. `D50_day` keeps its strict comparison.
- MAJOR (shared) non-finite precipitation was an "observation" in preprocessing (an
  isolated +Inf qualified the year and then voided every metric; −Inf tripped the
  negative-value rejection). FIX: NaN, ±Inf (and in Julia any non-Real) are MISSING for Q,
  PPT, SWE, TMIN and TMAX alike — counted, sign-checked and interpolated as such; the
  helpers require finite neighbours. No effect on real data (the inputs hold no Inf).
- MINOR (Julia) config readers accepted `true` as 1 and 0/1 as a flag; (shared) the
  validated-but-not-used cumulative fraction (0.500000000001 under the `prcp_50_dowy`
  name); season months and explicit policy overrides not validated (R overrides accepted
  3.5 / Inf / TRUE); Python's kernels coerced numeric strings while its preprocessor did
  not. FIX: strict typed readers in Julia too (Bool excluded from numbers, real Bool for
  flags, list elements checked), the fraction canonicalised to the percentage shown, whole
  months 1..12 enforced in all three, overrides validated, one `numbers.Real` observation
  contract across Python's entry points. Residual (documented): R's jsonlite simplifies a
  one-element array to a scalar before validation (not in the shipped config).
- Sonnet (Julia, GO): 2.1 M real annual values recomputed independently with 0
  mismatches; MINOR documentation — the interpolated Dec 31 belongs to the water year
  AFTER a leap calendar year and touches every metric of that year (mean 0.25 %, max 2.1 %
  of the annual total) — `docs/SIGNATURES.md` §14 corrected; `Printf` compat relaxed to
  "1" (the generated "1.11.0" would have refused Julia 1.9/1.10).
Regression tests for all of it in the three languages; suites green (Julia 355 climate
assertions + the full suite, Python 174, rpkg 1,179); cross-check unchanged (annual values
bit-identical, summaries within 4e-14 / 1e-12, including rpkg's Mann-Kendall).

**Blind review round 3 (Codex) — Julia GO; ports NO-GO on custom-input boundaries only.
Fixed the same day:** the canonical cumulative percentage must lie in 1..99 (an input of
1e-14 canonicalised to 0 % and emitted `prcp_0_dowy`); Julia's kernels take a Real-only
observation conversion (`coalesce_q` let a complex value through); Python validates all
eight preprocessing limits and the five policy flags when overridden and routes Arrow
decimal / object columns through the `numbers.Real` scalar contract; rpkg validates every
resolved limit and flag with the typed readers, checks season lists in the raw JSON,
validates explicit season overrides, and its default preprocessing config now carries the
configured season definitions (it fell back to hard-coded DJF for a custom config — real,
shipped config unaffected); exact integer counts (sums of logicals) are documented as
exempt from the sequential-sum convention. Suites green (Julia 366 climate assertions,
Python 175, rpkg 1,190); cross-check unchanged.

**Blind review round 4 (Codex, verification of round 3) — Julia GO, Python GO, rpkg NO-GO
on two MAJORs, fixed 2026-10-12/13:** (1) R still coerced inadmissible observations —
`as.numeric()` parsed numeric strings and dropped imaginary parts, so a string day gave
`prcp_annual = 367` and all-string temperatures a qualifying year with GDD 2007.5 where
Julia gives missing / no valid year. FIX: a package-level `.real_values()` mirroring Julia's
`_real_vector` (numeric and logical → value; character, complex, factor and other types →
NA; list columns element-wise) is now the only observation conversion for Q, PPT, SWE, TMIN
and TMAX in preprocessing and in both climate kernels. (2) The SWE policy ignored
`config$snow_na_policy` overrides and `.onLoad` read its three settings untyped (an
override of 0/0 still filled a gap and qualified the year; JSON `reject_negative_swe = 1`
silently qualified negative SWE). FIX: the SWE limits and flag are resolved from the config
with the typed readers (`.cfg_int` ≥ 0, `.cfg_bool`), default to the package values, and
are used throughout the SWE block; the default preprocessing config carries
`snow_na_policy`. Reinstalled; full rpkg suite 1,205 passed; cross-check unchanged (annual
values bit-identical, summaries ≤ 4e-14 / 1e-12). Codex round 5 (rpkg verification): its interim
log confirms both round-4 reproductions and the round-3 regressions pass in R and Julia
("inadmissible observations produce missing metrics, the SWE 0/0 override leaves the gap
missing and excludes the year, invalid SWE limits and flags raise errors"); the process then
sat idle for 3.5 h under the overnight memory pressure (below) and was stopped before
writing its report. Rerun 2026-10-13 on the quiet machine (report in 4 min): both round-4
MAJORs and every round-3 regression PASS; one new MAJOR — a list-typed input column came
out of `preprocess_daily_data` as a list again (each cleaned year was re-bound into the
ORIGINAL frame, whose column type won), so qualification succeeded and the runoff-ratio,
elasticity and flow-volume families then failed on the list column (579 fields missing for
a list Q). FIX: the observation columns present (Q, PPT, SWE, TMIN, TMAX) are converted with
`.real_values()` once on the whole input frame before the yearly loop; end-to-end
regression test added (list PPT with a string day → double columns, `Qann_mean` 365,
`annual_runoff_ratio_mean` 1, `prcp_annual_mean` 365); reinstalled, full suite 1,212
passed, cross-check unchanged. **Codex round 6 (2026-10-13): GO** — the round-5
reproduction passes for every observation column (1,128 result fields and all 378 annual
rows equal the numeric control), the round-4 and round-3 reproductions pass, 11 regenerated
Julia fixtures are byte-identical (9,256 finite annual values), 114 climate assertions.

**Sonnet 5.5 blind review of the rpkg port (resumed on the round-5 build) — GO:** 38,047
synthetic annual cells with 0 bit differences against Julia (edge cases, interpolation
stress, tie-rich quantised rain, ±Inf/NaN, a 400-year direct-call fuzz), summaries within
5.2e-14; seven real Daymet gages through both runners 11,232 annual rows bit-identical, 12 QA
flags equal; an exact-rational reference matches 13,595 cells. The three defects it found
during the review were the ones fixed in rounds 1–2 (strict crossing rule, `approx`
rounding, single-precision Mann-Kendall / exact Spearman). Three MINORs fixed 2026-10-13:
the roxygen block of `preprocess_daily_data` had been glued to the new gap-filling helper
(`document()` would have dropped the export) — helper moved above the block with `@noRd`,
roxygen regenerated (NAMESPACE now carries the four climate exports, Rd pages added); an
explicit JSON `null` in the `climate` / `na_handling` sections was accepted by R (jsonlite
turns it into "absent") where Julia rejects it — rejected now on the unsimplified JSON;
the cross-check fixture's header and `mk_rtol = 1e-6` default were stale (Mann-Kendall is
double precision now; default = rtol). Also: `Kendall` dropped from Imports (unused since
2026-10-12), stale "irreducible Spearman" / "1e-7 residual" statements corrected. NOTEs
recorded, no code change: Gini of a constant series can be −5.9e-18 (floating noise, never
in real data); the crossing allowance is 1e-12 × total, so it flips at 2e-12 mm for a 2 mm
year (as designed); the R runner needs ~1.6 GB more memory than Julia's and orders columns
by locale (the gates compare sets); the parity runs were all aarch64 (long double = double
there), so an x86-64 run of the ports remains undone; `STREAMFLOW_SIGNATURES_CONFIG` vs
`STREAMFLOW_CONFIG` stays on the deferred list. The rpkg reference re-baseline (decision
18) ran on the round-5 build, batched like the Julia products; these later edits are
non-numeric and were re-verified by the suite and cross-check after reinstalling.

**rpkg reference re-baselined on product #1 v2 (decision 18; 2026-10-13).** The rpkg port
ran product #1's configuration (WY 1993–2025, 60 %, climate family, reprocessed Daymet)
in the same 8 gage batches as Julia — 166 min of R, 6,678 gages, 0 errors — merged in the
Julia gage order (`~/HISSS_data/climate-rerun-12oct2026/rpkg_reference_1993_2025_13oct2026`,
2,277 columns, 27,602,160 annual rows). Gates against the Julia product: strict schema
equality PASS (identical column and gage sets); swallowed-failure scan PASS on all eight
logs; annual parquet equality PASS at the August tolerance (1e-6) with the residuals named
on the command line (`gate_annual_waived.txt`): 1 reference-only row (`avg_storage`, the
signed-zero `unique` divergence, unchanged since August) and 319 of 26,782,607 finite pairs
(0.0012 %), every one a discrete threshold crossing moved by a last-bit tie (`D*_day` 290,
`D25_to_D75` 59, `TQmean` 18, `dur_low_pulses_*` 2); NA patterns identical. At 1e-9 the
only additional class is `qp_slope_sd` (145 rows, ≤ 4.4e-11 relative on values in the
thousands). **The 39 climate signatures are bit-identical in every one of their annual
cells**, and the ~199 recession-family residual rows of the August run are gone. Identity R² (diagnostic): 2,228 of 2,244 columns Perfect (99.3 %; August 98.8 %), 10 Good,
6 Poor (all FDC90th), 0 below 0.95; mean 0.999925. The 2026-08-27 rpkg reference (1,653
columns, single-precision Mann-Kendall) is superseded; `docs/CROSS_LANGUAGE_STATUS.md`
carries the new baseline. After the hygiene edits the reinstalled package re-verified: full
suite 1,212 passed, both cross-check fixtures PASS (R now at the default 4e-14 for the
Mann-Kendall fields too), an explicit JSON null in `climate` or `na_handling` refuses to
load, the shipped config loads, the five exports resolve.

**Sonnet 5.5 blind review of the Python port — GO (final snapshot; 593 differential
case-runs against Julia, 410,156 annual rows and 368,560 summary values with 0 mismatches,
worst summary difference 2e-12 relative; real Daymet data on 65 gages and the whole batch
b1 (968,666 finite annual pairs) bit-identical; an independent naive reference of 20,054
annual values agrees; existing families bitwise unchanged with or without TMIN/TMAX).**
Four MINORs fixed 2026-10-13: PPT / SWE interpolation was written back for invalid years
(Julia writes only when valid) and their residual-NaN rejection depended on
`reject_residual_na` (Julia's is unconditional) — the two blocks now mirror `io.jl`
exactly; NumPy integer limits (`np.int64`) were rejected by the round-3 validators
(`numbers.Real` now, booleans still excluded); the climate functions and registries are
exported from the package. Q, PPT and SWE preprocessing now also go through the
`numbers.Real` observation contract (numeric strings were coerced by NumPy). NOTE items:
the 1e-12 crossing tolerance is an implementer decision (round 2) not in the guidelines —
recorded as such and confirmed by the user 2026-10-13, and its boundary is now pinned in all three
suites (storms of 1 ∓ δ mm: δ = 1e-12 crosses on day 7, 4e-12 on day 301); the Python
runner needs ~40 % more memory than Julia's and has no batching; an all-zero
precipitation year gives `prcp_max_dowy = 1` in every language (none in the input: 376,799
gage-years, minimum annual total 1.9 mm) — `docs/SIGNATURES.md` §14. Suites: Python 177,
Julia 369 climate assertions + the full suite, rpkg climate file 107. The Sonnet rpkg review
stalled (stream watchdog, no report) and is resumed after the reruns.

**FINAL production reruns (round-3 Julia; 2026-10-12/13).** Product #1 (WY 1993–2025):
8 batches in 2 h 8 min (129–2,088 s each while the reviewers shared the laptop), merge
against the delivered product, 0 per-gage family failures, gates 5/5, annual validation
PASS (894,846 non-NaN (gage, signature) means, 0 mismatches), output review PASS on every
gate with the pre-final numbers above — and **bitwise identical to the pre-final run**
(2,277 columns × 6,678 gages; all 27,602,161 annual rows), so the round-3 and round-4 fixes
changed no production value, as designed. Product #2 (WY 1980–2025): batches 1–2 finished
(2,629 s, 4,236 s), batch 3 was stopped by the runner's 120-minute cap at 0.1 gages/s —
the laptop's 8 GB swap was 88 % used by other applications and the 2.5 GB Julia process
crawled; resumed 2026-10-13 08:19 from the two finished batches (same command and
environment, verified by the runner) with the cap at 480 min: the remaining six batches took
151–909 s on the quiet machine (all eight: 2,629 / 4,236 / 909 / 258 / 151 / 155 / 167 /
159 s). The annual-parquet merge then failed at the runner's default 1 GB DuckDB budget
(34,681,331 rows; `--memory-limit 3GB` completed it in 30 s — the default is too small for
the 46-year product). Product #2 result: 0 per-gage family failures, gates 5/5, annual
validation PASS (842,401 non-NaN (gage, signature) means, 0 mismatches), output review
PASS on every gate: 1,252 Q-only + metadata columns identical on all 6,250 gages; the 24
continuous climate-dependent signatures R² ≥ 0.999 on shared gage-years (runoff ratios
≤ 2.0e-4 absolute on 223,715 shared gage-years; elasticity_annual 0.9994); discrete snow
timings changed in ≤ 0.05 % of gage-years; the old `_mean`s reproduced from the new annual
values on the old year sets (R² ≥ 0.999); 624 new columns populated (min finite fraction
0.766) and within the physical ranges; seasonal identity holds.

**Version-2 products, coverage and flags (from the CSVs).** Every gage of product #1
(6,678) and 6,223 of the 6,250 of product #2 carry climate summary statistics (the other 27 have a
Daymet series but exactly 20 qualifying years from WY 1980, and WY 1980 has no climate
value because the input begins 1980-01-01 — 19 climate years fall below the 20-value floor,
so their climate summaries are NaN while the annual values are in the parquet) (the 2026 products: 5,513 / 5,635 usable);
`annual_runoff_ratio` is present for 6,646 / 6,195 (the un-normalized gages excluded);
the frost-date trend statistics survive the record-anchored gate on 5,995 / 5,671 gages.
`flagged_for_high_na` is TRUE for 560 (#1) / 401 (#2) gages — regenerated by the run
(the delivered columns say 1,224 / 1,243 and are wrong; the 2026-09-04 recomputation on the
old columns gave 791 / 598; fewer now because every gage carries the climate families).
`area_normalized = FALSE` for 32 / 28 gages (gage sets unchanged; the 37-vs-32/28
discrepancy noted 2026-10-05 stays open) and `flagged_for_qann_range` for 159 / 151.
Folders: `~/HISSS_data/climate-rerun-12oct2026/processedOuts_{1993,1980}_2025_climate_12oct2026`
(signatures CSV, annual parquet, timing JSON with provenance, batch logs, validation and
review files, explorers and dashboards); nothing deposited (decision 20).

**Output-review design and the pre-final dry run (2026-10-12).**
`docs/benchmarks/compare_climate_rerun.py` classifies every column of a rerun against
the superseded product: Q-only + metadata must be identical per gage; climate-dependent
columns are gated on SHARED GAGE-YEARS (decision 16) — continuous signatures R² ≥ 0.999
and 99th-percentile relative difference ≤ 1e-3 on the annual values, discrete snow
timings reported as the fraction of gage-years changed, elasticity reported only — plus
the old product's `_mean` recomputed from the new annual values on the old year set
(R² ≥ 0.999); the summary-level trend statistics move for every gage because two climate
years were added and are reported, not gated; new columns must be populated and pass
physical sniff ranges. The pre-final product #1 (round-2 code) passed every gate against
`processedOuts_drought_28jul2026`: 1,252 Q-only + metadata columns identical on all 6,678
gages; runoff ratios within 1.8e-4 absolute on 163,422 shared gage-years; discrete snow
timings changed in ≤ 0.05 % of gage-years; old means reproduced with R² 1.000
(elasticity_annual 0.9993); `flagged_for_high_na` regenerated (1,518 gages change — the
delivered column was wrong, Known Issues).


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
- **DECISIONS (user, 2026-10-05) on the three items the review left open:**
  - Publish the mirror at every merge as the general preference; skip only for a good one-off
    reason (the 2026-10-01 deferral was one). Republished 2026-10-05 from a clean checkout of
    HEAD, since other sessions had uncommitted edits.
  - The sentence about the pasted token is trimmed from `EO_data_processing/CLAUDE.md`, which
    ships to the mirror.
  - The STATUS line "37 Canadian gages … 27/37" describes the canonical run; the products hold
    32 / 28 such gages, of which 22 / 25 are flagged. To be corrected in an upcoming update.

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

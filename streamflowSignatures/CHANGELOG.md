# Changelog

All notable changes to the Streamflow Signatures project.
For full historical detail (Dec 2025 – April 2026), see [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

> **Convention** — the most recent month(s) appear here as a condensed summary with pointers; full per-change detail lives in the archive. When a month's work is complete, condense it here and move the detail to the archive. File-level change lists belong in `git log`, not here; analysis and benchmark tables belong in the canonical docs (`docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`) and are linked rather than re-hosted.

## [Unreleased]

### Planned
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
- **Port campaign COMPLETE (2026-08-27) — both ports validated at full scale.**
  Retained below for the follow-ups it still carries. All six formerly Julia-only features (Pettitt changepoint fields, the
  20-value stats floor, the annual-values collector, the b=1 recession alpha, the 14
  snow metrics, the 10 drought metrics) are now implemented in BOTH ports; the
  "623-column April subset" description no longer applies to either. **Python is
  validated at full scale** against canonical Julia (1,653 columns / 6,678 gages,
  1,615 of 1,620 shared columns Perfect, zero below 0.99 — August 2026 entry).
  **rpkg is now validated too** (2026-08-27, all four gates green — see the August 2026
  milestone). Remaining follow-ups:
    - **Seven rpkg-side fixes identified 2026-08-26. Three are APPLIED (the runner
      ones — the benchmark they were queued behind was stopped externally, which
      made them safe to apply); four remain DEFERRED** because they change rpkg
      package source and would require an `R CMD INSTALL` that must not happen while
      a benchmark runs. None invalidates any run — each was checked against the
      active one:
      1. ✅ **APPLIED** — `options(warn = 1)` in `run_rpkg_benchmark.R`, so R prints each per-gage
         family-failure warning immediately instead of deferring and truncating at 50.
         Until then `check_signature_failures.py` correctly REFUSES to certify an rpkg
         log ending in that banner, and rpkg's assurance rests on the schema +
         annual-parquet gates.
      2. **rpkg reads `STREAMFLOW_SIGNATURES_CONFIG`, not the canonical
         `STREAMFLOW_CONFIG`** (`rpkg/R/config.R`). An experiment launched uniformly
         with `STREAMFLOW_CONFIG=variant.json` would silently use the variant in
         Julia/Python and the *installed default* in rpkg — the same class of hazard
         as the Julia precompilation gotcha. **Basis for calling it inert here:** the
         2026-08-26 launch command sets no config environment variable of either
         name, so both resolve empty and rpkg falls through to its bundled config.
         Confirm against that run's `timing.json` provenance block once written —
         the in-flight log does not record the selected config path.
      3. **Non-canonical fallbacks when the whole `na_handling` section is absent**
         (`rpkg/R/config.R`): rpkg would default `reject_negative_flow` to TRUE where
         Julia/Python default FALSE, and would leave the three SWE keys `NULL`,
         producing length-zero control-flow errors in `io.R`. **Latent only** — this
         branch runs only when the whole `na_handling` section is absent, and the
         bundled config's `na_handling` was diffed against canonical as
         byte-identical (no key differing, added or missing), with the installed
         package resolving `na_reject_negative_flow = FALSE` and all three SWE keys
         populated. That establishes the fallback is unreachable for the bundled
         config; confirm the run actually used it from its provenance block.
      4. ✅ **APPLIED** — **`report_interval` 500 → 50.** 500 was far too coarse to monitor an R run,
         which is slow enough that the first marker can be over an hour away —
         during the 2026-08-26 run no progress line appeared for 65 minutes.
         Drop it to ~50. **Correction to an earlier diagnosis in this entry: R
         does NOT block-buffer `cat()` here.** That was tested and disproven —
         `cat()` output appeared within seconds while the writing script still
         slept, on the local SSD *and* on the exFAT volume the log lives on. The
         silence was real: the loop genuinely had not reached iteration 500.
         (The Kendall `tauk2` warnings do go to stdout — 31 on stdout, 0 on
         stderr for a 3-gage smoke — so they are a usable liveness signal, but
         their per-gage rate varies far too much to serve as a progress count.)
      5. ✅ **APPLIED — the single biggest win of the day (~15x).** The runner read
         the 111.6M-row streamflow parquet with NO column
         projection (`run_rpkg_benchmark.R:118` → `read_parquet(path)`, which
         takes no `col_select`), pulling all 9 columns — roughly 6–8 GB resident —
         on top of ~3 GB of climate. On a 16 GB machine this drives the system
         into swap: measured mid-run at 20.3 of 21.5 GB swap used, 163 MB/s
         sustained page-ins, 6.8 GB in the compressor, and the R process's RSS
         squeezed to 0.41 GB despite holding ~210M rows. The run still progresses
         but loses roughly a quarter of its wall-clock to paging (CPU time ran at
         ~67–75 % of elapsed). Project to `gage_id, Date, Q` and let
         `add_water_year_columns()` derive the rest — the same bounded-memory
         treatment `run_python_benchmark.py` already received.
      6. **Emit a per-run identifier into the log header AND a MANIFEST written
         last** (all three runners, rpkg first). *Corrected 2026-08-26 after
         review*: an earlier version of this item proposed the log + timing JSON
         alone, which does NOT bind either to the CSV — an aborted rerun leaves a
         new CSV beside an old log and old JSON that still share an identifier, so
         they would validate each other while describing different artifacts. The
         identifier must instead live in a manifest written **after every other
         artifact**, recording the run ID plus each artifact's name and row/column
         counts; its presence then implies the whole set came from that run.
         Today a log can be tied to artifacts only by OUTPUT_PREFIX and by counts
         (all three runners, rpkg first). Today a run log can be tied to an artifact
         set only by OUTPUT_PREFIX and by the footer counts, and both describe a
         run's *configuration* and its output *shape* rather than the execution —
         two runs of the same config match on all of them. A UUID (or the run's
         start timestamp) written to both places, and required to match by
         `check_signature_failures.py`, makes the binding exact and unforgeable by
         copying or touching a file. Until then the gate's binding is documented as
         making accidental mispairing hard, not as proof of provenance.
      7. **SWE is merged only inside an `if ("PPT" %in% clim_cols)` branch**
         (`run_rpkg_benchmark.R`), so a climate parquet carrying SWE but no PPT would
         silently drop every snow key; Python merges the slice independently. Inert
         here — the run's own log records `columns: gage_id,date,PPT,SWE`.
    - Julia canonical cleanup (non-blocking, no rerun needed): `unique` uses `isequal`,
      so −0.0 and +0.0 are distinct values; numpy and R use `==`. Cost measured at
      **1 annual row in 18.9 M**. Julia should switch to `==` semantics when this code
      is next touched.
  Campaign plan, phase records, and the full Codex review table:
  `docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md`. The campaign was
  elevated in priority by the 2026-08-24 manuscript sync, whose draft claims every
  signature is computed by the cross-language Julia/Python/R library — a claim that is
  now true for Python and will be true for rpkg once its benchmark is gated.
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

### Known Issues (discovered 2026-07-14, Codex review — not yet fixed)
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
- **RESOLVED 2026-08-11 (workaround in place; the corrupt file is still on the drive) —
  the truncated climate parquet.** The user restored the 44 annual Daymet CSVs to
  `/Volumes/Untitled/daymet_1980_2023/` (10 GB) and the parquet was rebuilt from them
  with the new `docs/benchmarks/convert_daymet_csvs_to_parquet.py` →
  `processedOuts_feb2026/daymet_1980_2023_rebuilt_10aug2026.parquet` (97,757,220 rows ×
  6,087 sites, 3.76 GB). **The canonical `daymet_1980_2023.parquet` path still holds the
  truncated file** (left untouched at the user's direction), so runs must point
  `STREAMFLOW_CLIMATE_PATH` at the rebuilt file — the WY 1980–2025 wrapper now defaults
  to it. Details of the rebuild, its validation, and the Daymet 365-day calendar finding
  are in the August 2026 entry below. Original issue, kept for the record:
- **BLOCKING (data integrity, discovered 2026-08-10) — the climate input
  `daymet_1980_2023.parquet` is TRUNCATED on the thumbdrive.** The file is
  **1,261,436,928 bytes** with no `PAR1` footer, against the **4,125,630,653 bytes**
  recorded in the 28 Jul drought run's own provenance block — with an unchanged mtime
  (`2026-07-16T14:43:10`), so nothing in the filesystem metadata signals the loss. Julia
  fails at Phase 2 with `ArgumentError: invalid parquet: final bytes are "@…", expect
  "PAR1"`. **No other copy exists on the drive** (no Daymet ZIP, no annual CSVs, no
  second parquet) and none is present locally. Consequence: **no run that needs
  precipitation or SWE can be reproduced** — this blocks standard product #2's drought
  regeneration and any future re-run of #1. The streamflow parquet
  (899,935,349 bytes) and metadata CSV (2,164,196 bytes) are intact: sizes match
  provenance and both have valid footers. Recovery options, in order of preference:
  (1) another copy (S3 `s3://climate-ai-data-science-shiny-app-data/streamflow/`,
  another machine, or the original Windows `D:` drive); (2) rebuild from the Daymet ZIP
  of 44 annual CSVs via `convert_daymet_zip_to_parquet()`; (3) re-download from ORNL
  DAAC and re-run the co-authors' gdptools basin aggregation. The drive is **exFAT**
  (`/dev/disk4s1`, and carries a `FOUND.000` from an earlier fsck), which is a plausible
  cause — the delivered products living on it should be checksummed and backed up.
  **This is the first thing the July 2026 provenance block caught**: without the
  recorded input sizes the truncation would have looked like an unexplained Julia error.
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
- **RESOLVED 2026-08-24 (fixed in Julia — see August 2026 entry; verified inert for
  the delivered WY 1993–2025 product, zero masking events in that window) —
  seasonal runoff ratios ignored seasonal completeness flags**:
  `julia/src/runoff_ratios.jl:113-118` looks up flags named
  `winter_complete/spring_complete/summer_complete/fall_complete`, but the preprocessor
  emits `win_complete/spr_complete/sum_complete/fal_complete` (`julia/src/io.jl:462`;
  cf. correct usage in `flow_volumes.jl:136-137`). The existence check at
  `runoff_ratios.jl:120` silently fails, so seasonal runoff ratios are never NaN'd for
  incomplete seasons — the guidelines' 80% seasonal completeness rule is not applied to
  this category. **VERIFIED 2026-08-24 (port-gap survey): the bug is JULIA-ONLY.**
  Python (`io.py` emits full names; `runoff_ratios.py:120-123` looks them up) and rpkg
  (`io.R:418-420` / `runoff_ratios.R:54-57`) are internally consistent — their masking
  works. Fixing Julia therefore moves it TOWARD the ports. **FIXED same day** as Phase 0
  of the port campaign (`docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md`),
  with a regression test and a fresh reference benchmark attributing the delta (zero
  masking events in the WY 1993–2025 window — full record in the August 2026 entry);
  delivered standard products are unaffected.

### Guidelines Document TODOs
<!-- New suggestions from hydrology colleagues will be tracked here -->
<!-- Format: - [ ] Description (source: section name in guidelines doc) -->

**Update 2026-08-31 (later same day) — most queued fixes APPLIED in the doc; full
RESTRUCTURE DRAFT delivered; HISSS external-link decision.** A second sync confirms the
user applied in the Google Doc: the NA-handling four-condition wording, flow-volume
totals in mm, the pulse headline mis-pairing + low-pulse copy-paste fix, the recession
b-free / a-fixed-at-1 sentence (partially — the word "a" was dropped in editing), the
Storage / Snow / Drought sections, the high_na "numeric output columns" wording, five
new references, removal of the legacy Part 3 8-stat table, the "elasticity_annual"
typo, and the title's helperFunctions.R reference. The remaining residuals (Qxx
percentiles say "mm", should be mm/day; the low-pulse zero-flow note; the dropped "a";
drought notes truncated with a stray "ß" and missing four caveats + the ≥10-year
threshold rule; every "recommended addition") are ALL folded into a **full restructured
draft of the doc**, built on a co-author suggestion to integrate the Glossary and
Functions sections — one self-contained module per signature family (Metrics → Function
→ Method and parameters → Requirements and decisions; Parts 1–6, modules 3.1–3.14
matching the product's 14 categories, `name: definition` lines kept machine-diffable
for the auto-sync). Delivered as a paste-ready artifact; full text + change list
recorded at `docs/plans/2026-08-31-guidelines-restructure-draft.md`. Snapshot
re-synced (its header lists applied fixes and residuals). **DECISION (user,
2026-08-31): all external-facing repo pointers — including the guidelines doc's header
links — go to the public golden-standard repo https://github.com/CZ-Sync/HISSS**,
replacing the private-repo and CZ-Sync/code-sandbox links (repo-internal docs were
already clean from the 2026-08-28 audit; re-verified — zero stale references in the
tree; the doc-header swap ships with the restructure draft). Checkbox statuses below
updated accordingly.

**Synced 2026-08-31 — the guidelines doc MOVED to a NEW publish URL** (declared ground
truth by the user; snapshot overwritten, sync URL repointed in CLAUDE.md / README.md /
docs/DATA_SOURCES.md / docs/SIGNATURES.md). Major revision incorporating the
2026-08-30/31 review session (baseflow filter parameters + per-year application,
pulse/reversal definitions, elasticity equations, 16-statistic glossary incl. Pettitt
fields, rewritten QA-flags section; the storage section was REMOVED). Full three-way
alignment review run same day. Findings, by direction of fix:

*Doc-side fixes queued for Arik to apply in the Google Doc (code is correct):*
- [x] **Flow-volume units** — PARTIALLY APPLIED 2026-08-31: the five totals now say
  mm, but the Qxx percentiles were overcorrected to "mm" too (should be mm/day).
  Residual fixed in the restructure draft (§3.1).
- [x] **Pulse glossary mis-pairs thresholds** — APPLIED 2026-08-31 (headline
  definitions no longer attach period-of-record thresholds to `_year`).
- [x] **Low-pulse sub-bullet errors** — PARTIALLY APPLIED 2026-08-31: the
  high-pulse copy-paste is fixed, but "for intermittent streams, 0 flow days are
  excluded" survives. Residual fixed in the restructure draft (§3.5).
- [x] **Recession REGRESSION** — PARTIALLY APPLIED 2026-08-31: the sentence now
  states b = free fit and the b = 1 constraint for a, but dropped the word "a"
  ("parameter similarly determined") and still reads as a line-of-best-fit for a
  (it is a median). Residual fixed in the restructure draft (§3.4).
- [x] **NA-handling "items i, iii, iv flag (not remove)"** — APPLIED 2026-08-31
  (the corrected four-condition paragraph, incl. the raw-count-before-
  interpolation nuance).
- [x] **Leftover Part 3 8-stat table** — APPLIED 2026-08-31 (removed).
- [x] **avg_storage reinstated** — APPLIED 2026-08-31 (with the eligibility
  qualification and OMITTED-FROM-MAJOR-ANALYSES note).
- [x] **Snow and drought sections** — APPLIED 2026-08-31 (all 14 snow metrics +
  the drought family), except the drought notes pasted truncated (stray "ß";
  missing the severity-ladder, intermittent-stream, p10-redundancy, and units
  caveats and the ≥10-year threshold rule). Residual fixed in the restructure
  draft (§3.14).
- [x] **high_na flag wording** — APPLIED 2026-08-31 ("numeric output columns";
  the fuller denominator description remains available in the restructure draft).
- [ ] **Elasticity rolling window wording** (Codex finding): the "11-year window"
  is 11 consecutive QUALIFYING observations (`elasticity.jl` indexes the
  PPT-filtered valid series), usually but not necessarily 11 consecutive water
  years — same nuance already documented for `elasticity_annual`'s adjacent
  qualifying years. In the restructure draft (§3.9); not yet in the doc.
- [x] Minor — MOSTLY APPLIED 2026-08-31: Section 2 "ß" removed (a new stray "ß"
  landed at the end of the drought notes), "elastacity_annual" fixed, title no
  longer says "(in helperFunctions.R)". Still open (all in the restructure
  draft): "Only included changes" grammar; the ~92%/80% annual-completeness
  wording; the `min_Q_value_and_days` legacy note.
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

**Codex adversarial review of this sync + the proposed doc text (2026-08-31,
codex exec read-only): GO-WITH-FIXES — 4 MAJOR + 4 MINOR, all addressed same
day.** Codex independently verified clean: the URL swaps, the flow-volume-units /
pulse-threshold / zero-flow / recession-b=1 / NA-rejection / Pettitt-window
findings, the baseflow and MK glossary claims, and all 14 snow definitions. The
corrections it forced: (1) proposed NA-handling text wrongly said interpolated
≤3-day gaps "do not count against the year" — the raw NA count is taken BEFORE
interpolation (`io.jl` step b precedes step d), so interpolated days DO count
toward the 30-day ceiling; (2) drought deficit/threshold units are mm and mm/day
only for `area_normalized = TRUE` gages (raw m³/s units otherwise — the family is
NOT gated on area_normalized); (3) the drought severity-ladder monotonicity must
be stated as "non-decreasing from p2 to p30" (ambiguous "level" wording read
backwards against the D0→D4 framing); (4) the high_na denominator wording (new
queued item above); plus the ~92% annual-rule nuance, the 11-qualifying-
observation rolling window (queued above), an avg_storage eligibility
qualification (requires PPT + area-normalized flow + ≥10 candidate years / ≥5
annual values — NA otherwise), and a mechanical snapshot-fidelity check
(performed: token ratio 0.988 vs the published doc's extracted text; every
difference is markdown mechanics or one of three deliberate trivial typo
normalizations — no content invented or omitted). Corrected paste-ready doc
blocks delivered in-session.

Synced 2026-07-21 — first doc revision since 2026-04-15. Behavior-changing items:

- [x] **Trend completeness gate: overall series 80% → 60%** (source: NA Handling) —
  **IMPLEMENTED 2026-07-21** (see July 2026 entry below). Config-only change picked
  up by all languages; takes numerical effect at the next benchmark/production run.
- [x] **Decade gate 60% vs 80% — RESOLVED (user, 2026-07-22): 80% first/last
  decade + 60% overall CONFIRMED.** Matches the shipped config
  (`decade_min_fraction: 0.80`, `min_fraction: 0.60`), the guidelines doc, and
  manuscript §2.2.3 — no code change needed. The linked snow record-anchored
  gate therefore also runs at 0.80 decades.
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

Documentation-only / already-implemented items from this sync (no action):
- FDC glossary section added (FDCall/FDC90th/FDCmid — already implemented).
- Recession b free-fit + a via b=1 linear aquifer now in the doc — matches the July
  2026 Julia implementation.
- Q-P seasonality: "backward-looking" 30-day rolling window confirmed — matches the
  implementation (trailing window `[end-29, end]`, mid-window month attribution);
  new caveat that the rolling approach may differ from Wrede et al. (2015).
- Elasticity: annual (t1−t0) second calculation marked "achieved as of 4/16"; the
  <30%-missing counterfactual marked "added as documentation not filter" — resolves
  the previously-pending 30% diagnostic item as documentation-only.
- Runoff ratio NA rule now states the implemented thresholds (annual PPT < 10mm,
  seasonal < 1mm → NA).
- Flow volume glossary corrected to "total" (was "mean") streamflow — matches
  implementation; flow-timing "labeled as julian" fix-note removed (resolved).

### Manuscript Reconciliation Log

Session-start reconciliation of the HISSS manuscript draft (Scientific Data,
submission target Nov 9 2026) against code + repo docs. Snapshot:
`docs/MANUSCRIPT_DRAFT.md`; workflow: CLAUDE.md → Session-Start Workflow → B.

**2026-07-21 — baseline established + initial reconciliation pass.** Snapshot
created and the manuscript added to the session-start review. Findings from the
first pass (manuscript §refs; direction of fix in brackets):

- **§2.2.3 asserts the 60% overall trend gate** — consistent with the July
  guidelines doc, but code still enforces 80% (`trend_completeness.min_fraction`).
  [Fix in code — same item as Guidelines TODO #1 above; the manuscript and
  guidelines agree, the code lags.]
- **§2.2.3 repeats "items i, iii, and iv … flag (not remove) by default"** —
  item i (>3 consecutive NA days) currently REJECTS the year in the preprocessor,
  not config-toggleable. [Clarify with domain experts — same as Guidelines TODO
  #4; whichever way it resolves, manuscript, guidelines, and code must end up
  aligned.]
- **§2.2.2 says statistics were computed for "a subset from 1993-2025
  wateryears"** — the July 2026 production run used an explicit WY 1993–2022
  window (user-directed end cap; see
  `docs/plans/production_rerun_1993_2022_60pct_plan.md`). [Flag for authors —
  or confirm a future 1993–2025 run is intended.]
  **[CORRECTION + RESOLVED 2026-07-22]**: the July 14 run was actually
  WY ≥ 1993 UNCAPPED — its plan's §10 execution record superseded the 1993–2022
  title (this log's original claim mis-read the plan). Moot either way: the
  **WY 1993–2025 standard run executed 2026-07-22 (Codex results-review GO)**,
  so the manuscript's "1993-2025" wording is now accurate as written.
- **§2.1.2 says flow was normalized "by watershed area as defined in the
  watershed boundary shapefile"** — implementation normalizes by the published
  drainage areas (GAGES-II metadata for USGS; HYDAT STATIONS
  `DRAINAGE_AREA_GROSS` for Canada), not shapefile polygon areas; and 37
  Canadian gages carry raw m³/s (`area_normalized = FALSE`, no published area).
  [Flag for authors — wording fix + the un-normalized-gage caveat belongs in
  Usage Notes/limitations.]
- **§2.2.2 metric-family list** duplicates "baseflow" and omits runoff ratios,
  Q-P seasonality, snow metrics, and pulse/timing details; stats list ("slope,
  Spearman's Rho, p-value, mean and median") omits the second slope (linear),
  Mann-Kendall tau/p-value, and the Pettitt changepoint fields. [Flag for
  authors — incomplete enumeration of the delivered columns.]
- **§2.1.3 says Daymet was aggregated to 6,041 basins** via gdptools/agg_gen —
  repo docs record ~6,087 sites in `daymet_1980_2023.parquet`, and the July
  rerun found snow values for 5,622 gages (4,533 US + 1,089 Canadian). [Verify
  counts with authors — the gdptools aggregation is the co-authors' upstream
  step, not in this repo.]
- **§1 describes the dataset as "at the HUC-8 scale"** — products are per-gage
  watersheds (GAGES-II / ECCC polygons + HydroBasins fallback), not HUC-8
  units. [Flag for authors.]
- Placeholders to watch on future syncs: methods summary paragraph, input-data
  table, §2.2.1 (empty), "n=xx", metadata file "(name)".

**2026-07-21 — user decisions on the initial findings.**

- **§2.2.3 (60% trend gate) → code updated.** `trend_completeness.min_fraction`
  0.80 → 0.60 implemented same day (see July 2026 entry). Manuscript, guidelines,
  and code now agree; takes numerical effect at the next run.
- **§2.2.2 (1993–2025 window) → outputs will catch up to the manuscript.** Two
  "standard" products planned: WY 1993–2025 @ 60% (first), WY 1980–2025 @ 60%
  (second, later) — see Planned. No manuscript edit needed for the window; §2.2.2
  should ultimately describe both standard outputs once they exist.
- **§2.1.2 (area normalization wording) → manuscript edit** (code is correct).
- Remaining findings → manuscript edits, queued below for Arik to apply in the
  Google Doc after this documentation pass + Codex review.

**Manuscript edits queued (to be applied in the Google Doc by Arik):**

1. **§2.1.2** — replace "normalized volumetric measurements by watershed area as
   defined in the watershed boundary shapefile (Sect. 2.1.1)" with the actual
   method: normalization uses agency-published drainage areas (GAGES-II metadata
   for USGS gages; HYDAT STATIONS `DRAINAGE_AREA_GROSS` for Canadian gages), not
   polygon areas from the boundary shapefiles. Add (here or in Usage Notes): 37
   Canadian gages have no published drainage area and are retained in raw m³/s
   with `area_normalized = FALSE`; their Q-to-PPT signatures are structurally NA,
   and downstream users must filter on `area_normalized == TRUE` before
   cross-gage comparison of unit-carrying signatures.
2. **§2.2.2** — fix the metric-family list: "baseflow" appears twice; runoff
   ratios, Q-P seasonality, and snow metrics are missing (and storage is
   intentionally omitted per the 4/23/26 guidelines decision — worth stating).
3. **§2.2.2** — complete the statistics list: per metric the output is 8
   statistics (Theil-Sen slope, linear-regression slope, Spearman's rho +
   p-value, Mann-Kendall tau + p-value, mean, median) plus 8 Pettitt changepoint
   fields — not just "slope, Spearman's Rho, p-value, mean and median".
4. **§2.2.2** — once the standard runs exist, describe the two standard output
   windows (user, 2026-07-22): WY 1993–2025, and "entire period of record"
   operationalized as WY 1980–2025 — both @ 60% qualifying fraction. The current
   sentence ("full period of record and a subset from 1993-2025") is close but
   should state the 1980 start of the entire-record product explicitly.
5. **§2.1.3** — verify "6,041 basins": the Daymet parquet used by this pipeline
   carries ~6,087 sites (repo docs); reconcile which count the gdptools
   aggregation produced.
6. **§1** — "at the HUC-8 scale" should be "at the gaged-watershed scale"
   (GAGES-II / WSC basin polygons, HydroBasins fallback), not HUC-8 units.
7. **§2.2.3** — "Items i, iii, and iv are set in the config to flag (not remove)
   by default" is accurate for iii (negative Q) and iv (constant SD) but NOT for
   i (>3 consecutive NA days), which unconditionally rejects the year in the
   preprocessor. Pending the domain-expert clarification (Guidelines TODO above),
   either the sentence drops item i or the code makes item i config-driven.
8. **§2.2.2 (new, 2026-07-27)** — add the **streamflow drought** family to the
   metric-family list: duration + deficit at five fixed percentile thresholds
   (2/5/10/20/30 %, U.S. Drought Monitor D4–D0), computed on 7-day-smoothed flow
   with the unbiased Weibull plotting position, citing Adelsperger et al. (in
   review) and Laaha et al. (2017). Two deviations from the source method must be
   stated: only the FIXED (whole-record) thresholds are implemented (the variable
   day-of-year thresholds are too uncertain at the low levels with 20–46 years of
   record, 2 % falling below the smallest plotting position outright), and aggregation
   is by WATER year rather than the paper's climate year
   (Apr–Mar), so droughts crossing Oct 1 are split across two annual values. Also
   worth a Usage Notes line: drought values are record-dependent (thresholds come
   from each product's own window) and `drought_deficit_*` is unit-carrying.

**2026-08-24 — major co-author revision synced; second reconciliation pass.** The
draft was substantially rewritten since 2026-07-21: new §1 framing (related-work
comparison vs HYSETS/Caravan/CAMELS-SPAT), the §2 methods summary paragraph filled
in (16,994 candidate gages → 8,014 usable; 111.6M observations), §2.1.1/2.1.2
rewritten with concrete QA numbers, a new §2.1.4 HydroATLAS paragraph, and §2.2
restructured. Snapshot overwritten.

*Queued-edit status* (numbers refer to the 2026-07-21 queue above):
- **APPLIED**: #2 (metric-family list — baseflow deduped, runoff ratios/Q-P
  seasonality/snow added, storage omission stated), #3 (stats list in §1: both
  slopes, Spearman + MK, Pettitt — but see the leftover-sentence finding below),
  #4 (both standard windows with gage counts 6,678 / 6,250), #6 (HUC-8 → gaged
  watersheds), #7 (flag-not-remove sentence now exactly matches the code: removal
  = >3-day gap + >30 NAs; flag-only = negative Q + constant SD).
- **PARTIAL**: #1 (area normalization) — "all but 73 gages" now stated, but the
  text still points at §2.1.1 (boundary polygons) as the area source; the actual
  normalization uses agency-published drainage areas (GAGES-II / HYDAT STATIONS),
  stated nowhere, and the Usage-Notes caveat (filter `area_normalized == TRUE`)
  is still missing.
- **OPEN**: #5 (§2.1.3 "6,041 basins" vs repo's ~6,087 sites), #8 (drought family
  absent — see counts finding below).
  - **#5 repo-side RESOLVED by direct count (2026-08-24)**: the delivered Daymet
    parquet contains exactly **6,087** distinct basins (97,757,220 rows; DuckDB
    `COUNT(DISTINCT site_id)` on the rebuilt file), and all 6,087 match metadata
    gages. §2.1.3 should say 6,087 unless the gdptools co-authors can explain
    "6,041" (possibly a stale count from an earlier aggregation run). Also measured
    for §2.1.3/Usage Notes — Daymet coverage is a strict subset of each gage set:
    5,965 of the 8,014 compiled gages, **5,517 of product #1's 6,678**, and
    **5,638 of product #2's 6,250** have Daymet series; climate/snow-dependent
    signatures are structurally NA for the remainder. (6,250/6,678 are signature-
    product gage counts and are unrelated to the Daymet aggregation.)

*New findings (direction of fix in brackets):*
- **§1 + §2 now claim every signature is computed by an open cross-language
  library (Julia, Python, and R)** — currently FALSE: Python/rpkg reproduce the
  April 623-signature-column subset (no Pettitt, no stats floor, no annual-values
  collector, free-fit alpha instead of b=1, no snow, no drought). [Fix in code —
  this is the six-feature port queue in Planned; port initiative planned
  2026-08-24. Fallback if the port slips past submission: soften the claim.]
- **Signature counts predate the drought family**: "106 signatures … 90 annually
  resolved + 16 per-gage scalars" and "13 categories" describe the 1,488-column
  pre-drought state; delivered products carry 100 annually resolved + 21 scalars
  (121) across 14 categories. [Manuscript — same item as queued edit #8.]
- **§2.1.2 conflates the decade trend gate with gage inclusion**: "calculated
  signatures for gages only if they … included 80% of the possible WYs in the
  first and final decade of the window" — the 80%-decade rule is a per-signature
  TREND gate (§2.2.2 states it correctly), not a gage-inclusion criterion.
  [Manuscript.]
- **Leftover duplicate paragraph in §2.2.1**: the old "For the Eckhardt filter,
  we used both default parameters … including slope, Spearman's Rho, p-value,
  mean and median" sentence survives right after the new (correct) BFI paragraph
  and understates the 8 stats + Pettitt. [Manuscript — delete/merge.]
- **§2.1.1 vs §2.1.3 inconsistency**: §2.1.1 says the 7,964-polygon layer was
  used "to aggregate gridded climate data and remotely sensed LULC and LAI";
  §2.1.3 says Daymet was aggregated by the co-authors' gdptools run to 6,041
  basins (a separate, earlier boundary set). Only LULC/LAI/NLCD use the 7,964
  layer. [Manuscript.]
- **Closed**: §2.2.2's trend gates now match the code exactly — 20-value stats
  floor, 60% overall, 80% first/last decades, recession + elasticity exemptions.
  The 2026-07-21 60%-gate finding is fully resolved.
- §2.1.4's aggregation description omits that 91 `_u`/`_p` attributes pass
  through from the outlet basin (HydroATLAS's own upstream accumulation); only
  `_s`-only attributes are re-aggregated. [Manuscript nuance — optional.]
- Minor, to relay: title "a nalyzing"; "indented"→"intended"; "Corack"→"Corak";
  unclosed paren "(e.g., CAMELS-Chem,"; missing references for McMillan 2021,
  Hatchett 2021, Petersky & Harpold 2018, Arsenault et al. 2020 (HYSETS),
  Kratzert et al. 2023 (Caravan), Knoben et al. 2025 (CAMELS-SPAT); §2.1.2's
  retrieval end "30 September 2025" while the Feb 2026 parquet carries partial
  WY-2026 rows. [Manuscript.]

**2026-08-24 — drafting support for blank/revised sections (user-directed).**
§2.1.4 (MODIS land cover/LAI) revised text reviewed against the delivered EO
products — verified accurate (8 schemes, 102 class columns, LAI algorithm/range),
minor fixes suggested (percent vs fraction; v061 stated for both products;
reference-list update to the v061 dataset citations), plus a ready-to-paste Annual
NLCD paragraph (30 m, 1985–2025, CONUS-only, 6,119 watersheds, impervious surface,
~280× pixel density vs MODIS). **§3 Data Record drafted in full** around the planned
five-resource HydroShare structure (DOIs intentionally TBD) — delivered in-session
for Arik to paste into the Google Doc; structure follows
`docs/plans/2026-08-24-hydroshare-publication-and-dashboard-plan.md`. New references
flagged for the list: HYDAT database, MCD12Q1 v061, MCD15A3H v061, Annual NLCD
Collection 1.

**2026-08-28 — major revision synced (publication-audit session); third reconciliation
pass.** The draft gained a full author list (Kaiser, Tashie, Lowman, Jennings, Gorski,
Murray), a rewritten §1, a corrected §2, and a fully drafted §3 Data Record. Snapshot
overwritten.

*Resolved since 2026-08-24*: queued #5 (§2.1.3 now says **6,087** Daymet basins —
matches the repo's direct count); the §1/§2 cross-language library claim is now TRUE
(both ports validated at full scale 2026-08-26/27 — no softening needed); §2.1.2 no
longer conflates the decade trend gate with gage inclusion; §2.1.1 vs §2.1.3
aggregation-boundary inconsistency resolved; typos fixed ("a nalyzing", "indented",
"Hammon"); **§3's checkable numbers are all CORRECT** (1,653 = 20 metadata + 100×16 +
21 scalars + 12 flags; annual rows 18,898,406 / 24,366,487; 111,624,189 obs / 8,014
gages; Daymet 6,087 / 97,757,220 / CY 1980–2023; 7,964 polygons; 211 HydroATLAS attrs;
LAI/LULC/NLCD counts; the 16 statistic suffixes; the HydroShare collection id).

*Still open, manuscript-side (relay to co-authors)*:
1. **Queued #1 core error remains**: §2 preamble still says watershed area "extracted
   for each gage [from the boundary polygons] … used to area-normalize discharge" and
   §2.1.2 still points at §2.1.1; truth = agency-published drainage areas (GAGES-II;
   HYDAT `DRAINAGE_AREA_GROSS`). §5 Usage-Notes caveat (filter
   `area_normalized == TRUE`) still missing.
2. **Queued #8 still open**: §2 preamble + §2.2.1 still say "106 signatures … 90
   annually resolved + 16 … 13 categories" and omit drought; truth = **121 (100
   annually resolved + 21 scalars) across 14 categories**; drought methods paragraph
   (queued edit #8 text) still needed.
3. §2.2.1 dropped the avg_storage sentence entirely (it's still in the product —
   reinstate "computed but omitted from major analyses, no ET term").
4. §2.2.1 leftover duplicate BFI paragraph still present.
5. Missing references unchanged (McMillan 2021; Hatchett 2021; Petersky & Harpold 2018;
   HYSETS/Caravan/CAMELS-SPAT; DeCicco; Albers; HYDAT; gdptools; Annual NLCD;
   correct MCD15A3H v061; Adelsperger) + §1 cites nonexistent "Linke et al. 2013".
6. Carried typos ("Corack", unclosed "(e.g., CAMELS-Chem,"); placeholders (`<LINK TO
   REPO>` — **now fillable: https://github.com/CZ-Sync/HISSS**; metadata file "(name)";
   schematic; input-data table; §4/§5 stubs); §2.1.2 end date "30 September 2025" vs
   partial-WY2026 rows.

*New issues in this revision (manuscript-side)*:
7. §2.1.2 "47% reduction in valid gages" — arithmetic wrong: 8,980/16,994 excluded =
   **53% reduction** (47% retained).
8. §2.1.3 "6,087 basins that both met both the aforementioned data quality standards" —
   overstated: 122 of the 6,087 were never compiled; only 5,965 of the 8,014 processed
   gages have Daymet series (also fix "both met both").
9. §3 Resources 1–2: garbled sentence ("fitted checksums…"); filename typo
   "hiss_signatures_wy}window}.csv"; the gage_id join claim needs the Resource-3
   zero-stripped-ids nuance (documented join recipe).
10. Acknowledgements: "Claude Code 0.145.0" is wrong — 0.145.0 is the **codex-cli**
    version from this repo's adversarial-review records, not a Claude Code version;
    duplicated AI sentences; consider crediting the second review tool (Codex CLI)
    explicitly.
11. §3 "released under CC-BY 4.0" — no license decision recorded in-repo; confirm it
    matches the HydroShare resources. §3 labels the collection URL as a "DOI" — the
    DOI form is https://doi.org/10.4211/hs.f702201faa5d46069a5ee83ffa4c9768.
12. Small typos: ">6,000stream gages", "applied.The", "MOIDS", "WSG84", duplicate
    "daily", "teh geometry layer", "joins ot legacy", "and_pre_mk_pval",
    "nlcd_c{code}_pct}", stray "[tbd])".

*Code/docs-side*: none new — every checkable changed claim matches code, or the error
is manuscript-side. One repo-side stale comment found on the way and **fixed same
day**: `config/signatures_config.json` → `annual_values.comment` no longer says
"Julia only; Python/rpkg ports deferred".

---

## [August 2026]

### Milestone: pre-publication audit + PUBLIC CODE RELEASE at github.com/CZ-Sync/HISSS (2026-08-28)
Full audit of code + documentation ahead of the manuscript submission, followed by the
first publication of the project to the public **CZ-Sync/HISSS** repository (created by
the user; snapshot-mirror model). Five parallel audit passes (docs accuracy, package
docs/metadata, legacy inventory, code hygiene + sensitive-content sweep, manuscript
reconciliation), all findings fixed same day:

- **Sensitive-content sweep: ZERO blockers** — no credentials, tokens, or unintended
  personal emails anywhere in the 376 tracked files (all credential handling is
  env-var/`.netrc` indirection; the maintainer contact added to rpkg/DESCRIPTION the
  same day is the one deliberate address). Remaining private-path/infra references classified and dispositioned
  (personal Windows-profile paths removed with the deleted April logs and the
  `prod_1980_60pct` wrapper fix; Drive folder ID published as-is per user decision).
- **Stale port-status claims corrected everywhere** — README.md (rewritten in full as
  the public HISSS README: current 1,653-column parity numbers, complete 14-category
  signature table incl. snow + drought, 16-statistics convention, standard products +
  HydroShare link, MPL-2.0 LICENSE added), CLAUDE.md ("rpkg not yet gated" contradiction,
  annual-values "not implemented in ports"), DEVELOPMENT.md (1,264-column diagram,
  "17 modules" ×4), CROSS_LANGUAGE_STATUS.md ("Canadian HYDAT integration pending"),
  claude-skill (wrong negative-Q rejection claim; stale 3-value stats rule → the
  20-value floor; legacy-shim R example → rpkg), SIGNATURES.md summary-table count
  errors (22→21 flow volumes, pulse double-count of negative_ann, recession scalar
  accounting), `config/signatures_config.json` `annual_values.comment` ("Julia only;
  ports deferred") — synced into the wheel-bundled and rpkg-bundled config copies too.
- **All three package READMEs updated to the validated-parity state** + metadata fixes:
  rpkg README's self-contradictory validation section; DESCRIPTION 60+→100+, real
  maintainer email, stats/utils Imports; pyproject `tests`-package install bug + wrong
  homepage URL (→ CZ-Sync/HISSS); `calculate_negative_days` now exported from Python
  (matching Julia/rpkg); Julia README's ~551-key claim and missing six families;
  unused `HypothesisTests` dropped from Project.toml. Quick Starts now label
  `filter_qualifying_years` as the LEGACY path (production = `preprocess_daily_data`)
  and state what the bare `calculate_all_signatures` call does/doesn't emit.
- **Code hygiene** (all suites green after: Julia 2,798 assertions / 0 fail; Python
  143 pass; rpkg 1,021 pass): Python's stale `EXPECTED_SIGNATURE_BASES` (5 phantom
  Flow_Reversals bases + nonexistent `elasticity`) replaced with the current 100-base
  registry; vestigial rename indirection removed from fdc/flashiness; dead
  pre-allocation lines removed; MK comments de-scipy'd; `changepoint.jl` header now
  says Pettitt is production (BIC retained-but-uncalled); `io.jl` empty-diagnostics
  fallback gained the 2 missing na_cause columns; elasticity docstring min_years
  10→15; outlived "Phase-4 pending" runner comments removed; dashboard embedded-JS
  `SINGLE_VALUE` sets synced to their Python lists (newer scalars were unselectable);
  compare-script reports now record the resolved candidate path; smoke tests gained
  ENV path overrides in Python/rpkg.
- **Legacy removal** (user-approved list, ~29 files / ~11.6 MB, each checked for
  inbound references — one stale pointer survived and was caught by the post-audit
  Codex review, see below): `archive/` (18 files), the three stray April benchmark logs,
  `compare_restricted_vs_baseline.R`, bit-rotted `compare_outputs.py`, root+data_out
  `upstream_hydrobasins.RData` stubs, pre-2026 `data_out/summary_data.csv`, duplicate
  `combined_watershed_metadata_09feb2026.csv`, header-only `all_us-can-al-grdc.csv`,
  `logs.gitkeep`. References updated in the same commit.
- **Public mirror tooling**: `publish_to_hisss.sh` mirrors exactly the git-tracked
  files minus a fixed exclusion list (Claude tooling, `docs/MANUSCRIPT_DRAFT.md`,
  `docs/plans/`, the defunct Shiny app, the two >50 MB golden CSVs — user-approved) to
  HISSS `main`, stamping the source revision. Documented in DEVELOPMENT.md → Common
  Tasks and CLAUDE.md → Public Mirror; run after every merge to main.
- **Manuscript synced (major co-author revision) + third reconciliation pass** — see
  `[Unreleased]` → Manuscript Reconciliation Log (2026-08-28): 6 prior findings
  resolved by the revision (incl. §3 Data Record fully drafted with every checkable
  number correct), 12 items queued for co-authors (area-normalization wording persists;
  signature counts still 106/13-categories vs the true 121/14; a 47%-vs-53% arithmetic
  slip; "Claude Code 0.145.0" is actually the codex-cli version; `<LINK TO REPO>` can
  now be filled with https://github.com/CZ-Sync/HISSS).
- **Codex adversarial review of the full audit diff (2026-08-30, codex-cli 0.145.0,
  read-only): GO-WITH-FIXES — 4 MAJOR + 6 MINOR, all fixed same day.** It independently
  verified the science-facing claims clean (the 100-base Python registry exact against
  the Julia registries; no behavior regression in any code cleanup; README counts and
  parity numbers sound). The substantive catches: (MAJOR) Python's `validate_schema()`
  was still incompatible with the real product — reworked to model the full 1,653-column
  schema (20 production metadata columns derived from the runner + 100 bases × 16
  suffixes + 21 scalars + 12 flags, strict/non-strict modes, 13 new tests asserting the
  arithmetic; suite 156 passed); (MAJOR) the public README's R Quick Start used
  `arrow::read_parquet`, which keeps `Date` and breaks rpkg's `add_water_year_columns`
  — switched to the package reader; (MAJOR) root MPL-2.0 vs package MIT license
  conflict — **user decision: MIT everywhere** (root LICENSE replaced; python/ gains a
  LICENSE copy; rpkg was already MIT); (MAJOR) `publish_to_hisss.sh` treated ANY clone
  failure as an empty repo — now discriminates via `git ls-remote --exit-code` and
  aborts loudly on real failures, plus a deliberate empty-publish-set guard (MINOR).
  Also fixed: a surviving pointer to the deleted `compare_outputs.py` in
  `R/tests/generate_golden_outputs.R` (the one reference the deletion sweep missed),
  two overclaims in this very changelog entry (softened above), the README's
  collector-required implication and a broken anchor.

### Milestone: rpkg VALIDATED at full scale — the port campaign is COMPLETE (2026-08-27)
Both ports now reproduce canonical Julia's **1,653 columns for the same 6,678 gages** on
the WY 1993–2025 @ 60 % standard configuration. rpkg's run: 131.9 min, **0 gages
errored**, all four acceptance gates green — strict schema equality with **no waivers**,
zero swallowed family failures across a 997,194-line log, cross-language annual-parquet
equality (two named residuals), and an identity-R² diagnostic of **1,601 Perfect
(98.8 %) / 10 Good / 9 Poor / 0 below 0.95**, mean R² 0.999843.

**The 19 non-Perfect columns are entirely pre-existing irreducible classes** — 11
FDC90th + 3 FDCmid + 2 Q90 (near-zero-tail OLS on `log10(Q + 1e-10)` and the Pettitt
tie flips downstream of it) and 3 recession Spearman p-values (exact permutation vs
t-approximation at small n). No drought, storage, snow, BFI or changepoint family
appears. Annual parquet: 0 NA-pattern mismatches, 18,898,405 of 18,898,406 rows shared,
**518 of 18,217,552 finite pairs (0.0028 %) over 1e-6**.

**The first rpkg benchmark passed its unit suite and still failed three of four gates,
and every finding was a real defect** (full detail in docs/CROSS_LANGUAGE_STATUS.md):
zero-row crashes in snow (and latent in flashiness/timing/baseflow) that the
orchestrator swallowed, costing 4 gages all 224 snow columns; a CONUS↔AKHIPR
`intersect()` that dropped 4 metadata columns for ~9,000 gages; pre-allocated annual
frames exporting 1,039 rows Julia omits, plus a missing ≥10-**distinct**-Q storage
guard; and the big one —

- **`smooth_daily_flow` used R's `mean()`.** Its long-double accumulator plus correction
  pass returns a different last bit than Julia's sequential `s += v`: 4,513 of 10,227
  smoothed values differed by ≤ 3.6e-15. Harmless in isolation, except the drought
  thresholds are percentiles OF that series, so a threshold landing on a flow plateau
  flips every plateau day at once through the strict `<` (gage 01589795 WY2002 read 116
  days in Julia and 60 in rpkg). An explicit sequential double sum makes the smoothed
  series **bit-identical across all 10,227 days**, and is marginally faster. This was
  ~11,578 of the 12,096 annual mismatches. **numpy's `mean` matches Julia's summation
  order for windows this short, which is why the Python port never exhibited it** — the
  clean Python drought result was the clue that bit-parity was achievable, and filing
  this under the documented tie-sensitivity class was the wrong initial call.

**Performance**: the runner read all 9 columns of the 111.6M-row streamflow parquet,
driving a 16 GB machine into swap (measured mid-run: 20.3/21.5 GB swap, 163 MB/s
page-ins, RSS squeezed to 0.41 GB, CPU at 54 %, 0.1 gages/s → 22 h projected).
Column-projecting to `gage_id/Date/Q` cut the table to a measured **2.50 GB** and the
run to ~2 h — a ~15x throughput gain, verified behaviour-neutral by identical
processed/skipped tallies at every checkpoint.

**Accepted, not fixed**: Julia's `unique` uses `isequal`, so `-0.0` and `+0.0` are
distinct in the storage distinct-value guard while rpkg's (like numpy's) uses `==`. A
pre-run review recommended replicating the quirk so the strict gate passes; declined —
the user ruled the ports are the correct side and Julia should move to `==`, so
reproducing it downstream would propagate an agreed-wrong behaviour into two languages.
Cost: **1 annual row in 18,898,406**, waived by name at gate time.

### Milestone: PYTHON reaches full feature parity with canonical Julia — VALIDATED at scale (2026-08-26)
The August port campaign's Python track is complete and gated. Python now produces the
**same 1,653 columns for the same 6,678 gages** as canonical Julia on the WY 1993–2025
@ 60 % standard configuration, carrying all six formerly Julia-only features (Pettitt
changepoint fields, the 20-value stats floor, the annual-values collector + parquet, the
b=1 recession alpha, 14 snow metrics, 10 drought metrics + 5 scalars).

- **Strict schema gate PASS with no waivers** — identical column set and gage set. Over
  the 1,620 shared signature columns: **1,615 Perfect (R² ≥ 0.999, 99.7 %)**, 5 Good,
  and **zero** in any lower tier. Mean R² 0.999988, min 0.997935. The annual parquet
  matches on the 100-signature set with 0 duplicate keys, **0 NA-pattern mismatches**,
  and 18,898,405 of 18,898,406 rows shared.
- **No divergence class remains** — both residuals are characterised rather than merely
  tolerated. (1) **Pettitt ties**: `cp_year` agrees on 597,505/597,527 cells (99.9963 %);
  where a rank tie flips the changepoint the segment split moves with it, which accounts
  for all 5 Good columns (FDC90th, Q90, Qsum) and, **at a 1e-6 tolerance, 267 of
  18,217,552 finite annual pairs (0.0015 %) — every one a discrete threshold-crossing
  metric** (`D*_day`, `D25_to_D75`, `TQmean`, one `dur_low_pulses_all`), where a
  last-bit tie moves a whole day. (2) **Signed-zero `unique` semantics**: 1 annual row
  in 18.9 M (an `avg_storage` row present only in the reference) — Julia's `unique` uses
  `isequal` so −0.0 ≠ +0.0, while numpy and R use `==`. Logged as a non-blocking
  canonical cleanup; **no rerun required** (user decision, 2026-08-26).

  **State the tolerance when quoting these counts.** Tightening to 1e-9 adds 198 more
  rows — `qp_slope_sd` (116, max 1.3e-07) and `elasticity_annual` (82, max 5.5e-08) —
  which are genuinely continuous quantities differing at ordinary floating-point noise
  between statistics libraries, not tie flips. 465 @ 1e-9 and 267 @ 1e-6 are the same
  result read at two thresholds, not a discrepancy.
- **Mann-Kendall p-value convention settled.** scipy's `kendalltau` selects on TIES, not
  sample size — exact when untied, asymptotic **without** continuity correction when
  tied. Julia and R both use the continuity-corrected normal approximation, and
  `Kendall::MannKendall` reproduces the Julia formula exactly, so **Python was the sole
  outlier**. Python was moved to the canonical formula (user decision, option C), cutting
  significance-call disagreement from 0.24 % → 0.0009 % on the main path and 0.45 % →
  0.0009 % on the Pettitt segments.
- **rpkg is code-complete and unit-tested (1,008 assertions against the INSTALLED
  package) but NOT yet validated at scale** — its benchmark is running. Do not cite rpkg
  parity numbers until that run passes the same three gates.

### Fixed: two wiring defects that unit tests structurally could not catch (2026-08-26)
Both were in rpkg, both presented as clean warning-only degradation, and both were found
by running the actual pipeline end to end on real data rather than testing modules in
isolation:

- **The rpkg benchmark runner did not pass the ported arguments at all** — no
  `changepoint`, `min_values_for_stats`, `collector`, or `snow_data`, and it discarded
  SWE at the climate join; rpkg's bundled config carried no `changepoint` or
  `stats_floor` section whatsoever. A 19-hour run in flight was **killed and discarded**,
  since it provably could not have produced the ported columns. Runner and config fixed;
  a 3-gage smoke now confirms the full key set (1,620 signature keys, 800 Pettitt, 224
  snow, 165 drought, 100-signature collector).
- **rpkg's drought family returned an all-NA family for every real gage.** rpkg's
  `preprocess_daily_data()` renames `date` → `Date` internally (io.R) while the ported
  drought module checks the Julia/Python lowercase name, so the missing-column branch
  fired behind a warning on every gage. The unit tests could not catch it: they build
  frames with `date` directly and never traverse the preprocessor → drought path. Fixed
  to accept either spelling, plus a new end-to-end test that goes through
  `preprocess_daily_data()` **and** the orchestrator. Post-fix the smoke gives
  `drought_duration_fixed_p10_mean` = 36.52 days/yr, matching the p × 365.25
  construction target.

**Verified NOT a defect** on the way: per-gage key-set variation is legitimate and
canonical. A gage that fails the 25-event recession gate emits 64 fewer keys (the
Pettitt fields of 8 recession/parameterized-BFI bases) because Julia's own
`empty_stats` returns the 8 statistics ONLY — no changepoint keys. Python's and rpkg's
`empty_stats` mirror that exactly, and the writers' union fills NA. A strict per-gage
key-set equality assertion would therefore false-positive.

### New: `check_signature_failures.py` — gate a run log for swallowed family exceptions (2026-08-26)
`calculate_all_signatures()` wraps every signature family in a try/catch in all three
languages, so an exception degrades to a warning and that gage simply loses those
columns — which reappear as ordinary NA in the rectangular output CSV, indistinguishable
from a legitimately non-computable value. Neither the unit suites nor the column-set gate
can see it. The new tool scans a run log for the per-gage failure lines all three runners
emit and exits nonzero on any unwaived failure (`--allow-family` waives explicitly and
still reports what was ignored). Validated against the completed Python run (**PASS, 0
failures**) and synthetic positive controls.

**It also refuses to certify a log that ends in R's deferred-warning truncation banner**
("There were 50 or more warnings"), because R defers warnings and truncates at 50 — a
pass cannot be substantiated from such a log. `run_rpkg_benchmark.R` needs
`options(warn = 1)` so each failure prints immediately; that must be applied **between**
runs, never mid-run, so it is queued behind the benchmark currently executing (tracked in
`[Unreleased]` → Planned). Until then rpkg's assurance rests on the schema, identity-R²,
and annual-parquet gates, which would light up on any family failing systematically.

### Milestone: HydroShare COLLECTION created; all five HISSS resources staged + reviewed (2026-08-25)
The HydroShare collection for the HISSS deposit now exists (verified live in-browser):
**"Hydrologic Information Signatures and Summary Statistics"** —
https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/ — which fixes
the collection DOI that activates at publication:
**https://doi.org/10.4211/hs.f702201faa5d46069a5ee83ffa4c9768**. Recorded in the plan
doc (`docs/plans/2026-08-24-hydroshare-publication-and-dashboard-plan.md`), project
memory, and appended to all five staged resource READMEs. Staging status: **all five
resources are deposit-ready** in `~/Downloads/Signatures/resource{1..5}_*/` — deposit
file names, per-resource READMEs + data dictionaries + abstracts, every document
Codex-adversarially reviewed with findings fixed (R1–R2 signatures, R4
geometry+attributes incl. the same-day boundary rebuild in one session; R3 inputs +
R5 EO products in a parallel session, R5 carrying the user's NLCD QA sign-off).
Remaining before submission: create the 5 member resources on HydroShare (ids → DOI
strings for manuscript §3), upload + integrity-check, coverage metadata, co-author
profiles/author order, reviewer private links; publish at acceptance per plan D4.
**DECIDED (user, 2026-08-25): HydroShare is GROUND TRUTH going forward** — once a
resource is uploaded and integrity-verified there, the flash drive and Google Drive
are convenient backups only; on any discrepancy between copies, HydroShare wins.
This makes the upload-time checksum verification step load-bearing.

### New: HydroShare Resource 5 (vegetation + land cover) staged + verified (2026-08-25)
`~/Downloads/Signatures/resource5_eo_products/`: the three EO products staged from the
Drive backup (`~/Downloads/LAI-LULCC/`) under deposit names — `hisss_modis_lai_monthly.parquet`
(195 MB), `hisss_modis_lulc_annual.parquet` (63 MB), `hisss_nlcd_annual.parquet` (31 MB) +
the three explorers (`hisss_modis_lai_explorer.html` 47 MB, `hisss_modis_lulc_explorer.html`
47 MB, `hisss_nlcd_explorer.html` 12 MB) — every copy md5-verified against its backup
source, and every parquet re-verified against the documented invariants: LAI 2,150,280 =
7,964 × 270 (2002-07…2024-12) with `good_coverage_frac` + `partial_month` (i.e. the 30 Jun
final), 0 dup keys, lai_mean ∈ [0.10, 5.80], exactly 17 always-NA urban basins; LULC
191,136 = 7,964 × 24 (2001–2024), IGBP sums 100 ± 1e-13, 0 dups; NLCD is the **finalized**
build — 250,879 = 6,119 × 41 (1985–2025), all finalize markers present (`nlcd_collection`
= C1V2, `valid_area_km2`, 4 QA flags), class sums 100 ± 1.3e-12, impervious 0–68.8 %.
All three explorers end cleanly (no exFAT-style truncation) and are self-contained except
the Leaflet CDN. New docs: `hisss_readme_eo.md`, merged `hisss_eo_data_dictionary.csv`
(173 rows: all 168 columns of the three parquets — incl. a note row for LAI's vestigial
pandas `__index_level_0__` column, kept for byte-identity — + NLCD's 4 `CAVEAT_*` rows +
the 5 columns of the exclusions CSV), abstract, and rename-manifest rows. Two backup
gaps handled: the original `nlcd_out_of_footprint.csv` was not in the Drive backup, so
`hisss_nlcd_excluded_gages.csv` was **exactly reconstructed** (EO-universe US gages absent
from the NLCD table = 45 gages, all verified in Alaska, lat 55.2–70.3° N); the MODIS
granule manifest was also absent — documented in the README as repo-level provenance
rather than shipped. Coverage measured for the README: MODIS covers 6,599/6,678 of
product #1's gages and 6,196/6,250 of #2's; NLCD covers 5,419 and 4,998. Deposit choice:
parquet-only for the data tables (the backup's CSV twins are exact re-serializations;
the 695 MB LAI CSV in particular is redundant). **User sign-off 2026-08-25 on the full
R5 staging — including the NLCD human QA via `hisss_nlcd_explorer.html`** (the pass owed
since 24 Jul; that closes the last NLCD publication gate). R5 is deposit-ready.

### Lost → rebuilding: the watershed geometry product; plus exFAT casualty #3 (2026-08-25)
`watershed_polygons_26jun2026.{gpkg,parquet}` + `_qa.csv` (the 7,964-basin delivered
layer) existed **only on S3 and was not captured in the Google Drive backup** — the
first permanent casualty of the S3 access loss (user-confirmed absent from Drive;
searched the flash drive exhaustively; SageMaker-side `geometry_7964.gpkg` also
unreachable). **Rebuild launched same day** from the committed pipeline
(`EO_data_processing/geometry/`, consolidated into
`~/Downloads/geometry_rebuild/rebuild_watershed_polygons.py`) with: a **fresh USGS
GAGES-II download** (204,809,860 B, zip integrity verified — the flash-drive zip is
**truncated to exactly 18 MiB, exFAT casualty #3**, and its extracted folder is also
partial at 107 MB vs 1.4 GB real; the repo docs' "18 MiB ZIP" note was recording the
already-truncated size); **fresh ECCC MDA_ADP** (the server serves a single version,
Last-Modified 2024-09-05 — slightly newer than the documented June-2024 pin;
the flash `ca/` copy was always a stub, 7 of 11 files zero-byte); and the flash drive's
`basinAt_NorAm_polys.gpkg` (verified: 167,665 features, EPSG:4326). The **54-basin
exclusion** (geometry > 100,000 km² + `05KH009`) was an uncommitted post-step in the
June delivery and is now implemented explicitly in the rebuild script. Validation
targets from the June QA record: 7,964 features; gagesii 6,164 / wsc_eccc 1,771 /
hydrobasins 29; 0 dups/invalid/empty; ~141 full-res-kept basins. **Caveat for R4 and
the EO products**: LAI/LULC/NLCD were computed on the ORIGINAL geometry; the rebuild
is expected to be equivalent but not bit-identical (ECCC source version differs), and
the R4 README will carry explicit rebuild provenance. Separately, the **HydroATLAS
attributes product was recovered from the Drive backup** (`HydroATLAS/` folder:
8,014 × 211 + its 211-row dictionary) and is staged + verified into Resource 4.

**Rebuild COMPLETE same day — every recorded June target matched EXACTLY**: 7,964
features; gagesii 6,164 / wsc_eccc 1,771 / hydrobasins 29; 0 dups/invalid/empty;
**141 full-res-kept basins (identical to June)**; area QA median 0.0010 with 97.8%
< 10%; 41.7 MB gpkg. Canadian matching produced 1,823 even on the Sep-2024 ECCC
version (per-district wanted/got identical in effect), the 31-gage residual resolved
31/31 outlets (15 by point-in-basin), and the 54-basin exclusion decomposed exactly
as predicted (52 wsc_eccc + 2 hydrobasins). This also **settles the geometry-source
question definitively: gagesii == 6,164 == every US gage; all 29 HB fallbacks are
Canadian** — confirming the manuscript's "100% of US gages from GAGES-II" claim and
the README_NLCD correction. Set relationship measured: attributes ∩ boundaries =
7,960; 4 geometry-only (the inclusive-universe extras); 54 attributes-only
(size-excluded). Staged as `hisss_watershed_boundaries.{gpkg,parquet}` +
`_qa.csv` in `resource4_geometry_attributes/`; rebuild workspace + script retained
at `~/Downloads/geometry_rebuild/`.

**Resource 4 docs Codex-reviewed same day: GO-WITH-FIXES (4 MAJOR + 4 MINOR) → all
fixed → delta review 8/8 RESOLVED + 2 residuals it exposed → fixed per its
prescribed wording → done.** Substantive catches: the "same committed pipeline"
provenance claim overstated (the 54-basin exclusion was newly implemented; ECCC
source is the Sep-2024 serving); "joins to every HISSS table" overstated coverage
(overlap measured: 7,960 shared / 4 geometry-only / 54 attributes-only); and the
inherited HydroATLAS dictionary carried two `x?` placeholder rows — resolved against
the BasinATLAS v10 catalog PDF, which also corrected a wrong unit (`ero_kh_uav` is
kg/hectare/yr, not kg/m²/yr; `hdi_ix_sav` is index × 1000) — plus 11 ambiguous
categorical-aggregation rows now stating the actual per-variable method (glc/pnv
argmax vs SUB_AREA-weighted mode). Codex independently verified all counts (7,964;
6,164/1,771/29; 141 full-res; 56 low_confidence; 8,014 × 211; dictionary set
equality; the overlap sets).

### Fixed: SECOND silent exFAT truncation — product #1 explorer restored from the Drive backup (2026-08-25)
The HydroShare staging review of Resources 1–2 (`~/Downloads/Signatures/`, pulled from
the Google Drive backup) found the thumbdrive copy of
`processedOuts_drought_28jul2026/signature_explorer_1993_2025_60pct_drought.html`
**silently truncated to 30,146,560 bytes** (ends mid-data, no `</html>`) against the
71,782,193 bytes recorded in RUN_NOTES/CHANGELOG — same failure signature as the
Daymet parquet on 2026-08-10 (size shrunk, mtime preserved). The Drive-backup copy has
the full recorded size and a clean ending; the thumbdrive file was **restored from it
and md5-verified**. Everything else checked bit-identical between the Drive backup and
the thumbdrive: both products' signatures CSVs, annual parquets, timing JSONs, #2's
production log, and all 200 explorer sidecar files (md5, file-by-file). The thumbdrive
has now silently truncated **two** files — treat it as failing; never hold a sole copy
there. Staging gaps fixed the other direction: RUN_NOTES, validation reports, and
vs-julia comparison artifacts were absent from the Drive-backup staging folders and
were copied in from the thumbdrive. Remaining for deposit readiness:
`signature_categories_260717.csv` is stale (90 pre-drought bases, 0 drought rows),
deposit renames pending co-author sign-off on the §3 structure, per-resource
README/dictionary not yet assembled.

### New: HydroShare Resource 3 (harmonized inputs) staged + documented, subagent-reviewed (2026-08-25)
`~/Downloads/Signatures/resource3_input_data/`: the four input files staged under deposit
names (`hisss_streamflow_daily.parquet`, `hisss_gage_metadata.csv`,
`hisss_daymet_basin_daily.parquet`, `hisss_hydat_interference.csv`), each md5-identical
to its canonical source (thumbdrive sizes + `PAR1` footers verified first). Built
`hisss_readme_inputs.md`, `hisss_input_data_dictionary.csv` (all 34 columns across the
4 tables, values profiled from the data), and the resource abstract. **The user caught a
framing error the docs shipped with**: "1,601 un-normalized Canadian gages" is the
all-candidates metadata count — only **73** are compiled and present in the streamflow
table (1,523 of the rest are `no_data`); all three docs corrected to use each count in
its proper universe. **Codex stalled twice** (hung before creating a session file —
tell: no new `~/.codex/sessions/.../rollout-*.jsonl` + near-zero CPU; two runs completed
fine earlier the same morning), so per user direction the adversarial review ran as a
**Claude subagent** instead: verdict **GO-WITH-FIXES**, 1 CRITICAL + 3 MAJOR + 5 MINOR,
all fixed. The CRITICAL: **`hisss_gage_metadata.csv` stores US gage ids ZERO-STRIPPED**
(7,046 seven-char ids; a naive join to the zero-padded parquets matches only 3,118 of
8,014 gages) — the docs had claimed zero-padded. Resolved by documenting (join recipe +
warning + corrected dictionary row) rather than re-padding the data, preserving the
resource's byte-identical-inputs provenance. Other fixes: Daymet coverage is 5,965 of
8,014 compiled gages (the 6,087 sites include 122 never-compiled candidates); the
`flag` column's real vocabulary (Canadian rows carry only B/E/A/NULL — no D/R; bare "A"
is agency-dependent; NULL ≠ "no qualifier" for US rows); ~588k NA-discharge rows exist
(missing days are NOT solely absent rows); rebuild-provenance wording tightened to the
CHANGELOG's precise claims; REGULATED has 371 empties; STN_REGULATION table credited;
1,528 = no-data + insufficient-years; float-typed calendar helpers noted. The reviewer
empirically verified all counts/extents, the 365-day Daymet calendar across all 267,828
site-years, and cross-document consistency.

### New: HydroShare deposit docs for Resources 1–2 — README + full data dictionary, Codex-reviewed (2026-08-25)
Staging (`~/Downloads/Signatures/`) renamed to the deposit convention
(`hisss_signatures_wy{window}.csv`, `hisss_signatures_annual_*.parquet`,
`hisss_signature_explorer_*` — explorer HTMLs' internal sidecar-folder references
patched to the renamed `_annual/` folders; `hisss_rename_manifest.csv` maps old → new).
Built per-resource `hisss_readme_wy{window}.md` and a generated
**`hisss_data_dictionary.csv` covering all 1,653 columns** (generator asserts full
header coverage; per-statistic units derived from per-base units), plus
`hisss_signature_categories.csv` updated to 100 bases (+10 drought, new "Drought"
class; unexplained `question` column dropped). **Codex adversarial review
(codex exec, read-only): NO-GO → all 10 findings fixed → delta review 9/10 RESOLVED,
1 PARTIAL → fixed → done.** The CRITICAL: docs claimed a Pettitt changepoint year is
"always reported" — in fact Pettitt fields are NA without ≥ 20 valid values and ≥ 10
per segment **within the WY 1980–2024 changepoint window** (config
`changepoint.end_water_year = 2024`), i.e. WY 2025 is excluded from all Pettitt
fields — a product-level fact now documented in both READMEs and all 800 Pettitt
dictionary rows. Other MAJORs: snow signatures are NOT suppressed for
`area_normalized = FALSE` gages (only Q-to-PPT is); flow totals for un-normalized
gages are m³/s·days (×86,400 → m³), stated for Qann/seasonal totals AND drought
deficits; "14 categories" claim conflicted with the shipped nine-class exploratory
grouping. Codex verified counts (6,678/6,250 rows, 1,653 cols, 100×16 + 21 + 12
schema, parquet row counts vs provenance) and the gate thresholds. Copies of
dictionary + categories are identical (md5) at staging root and in both resource
folders.

### Fixed (MEDIUM → resolved): seasonal runoff-ratio completeness masking was dead code in Julia (2026-08-24)
`julia/src/runoff_ratios.jl` looked up seasonal completeness flags by FULL names
(`winter_complete`, …) while the preprocessor emits ABBREVIATED names
(`win_complete`/`spr_complete`/`sum_complete`/`fal_complete`) — the existence check
failed silently, so the guidelines' 80% seasonal rule was never applied to seasonal
runoff ratios (Known Issue since 2026-07-14). Fixed to the abbreviated names (the
same lookup `flow_volumes.jl` always used), Julia-first, as Phase 0 of the port
campaign. Key facts established on the way:
- **The bug was JULIA-ONLY** — Python and rpkg emit AND look up full names
  consistently, so their masking always worked. The fix converges Julia toward the
  ports (2026-08-24 port-gap survey; contrary to the Known Issue's "ports likely
  mirror" assumption, now corrected).
- **New regression test** `julia/test/test_runoff_ratio_seasonal_mask.jl` (32
  assertions) asserts the mask actually FIRES with preprocessor-shaped flags — the
  bug survived precisely because no test did. Also wired the previously-orphaned
  `test_changepoint.jl` into `runtests.jl` (Codex review find: a green suite never
  ran it). Suite: 2,798 assertions, 0 failures.
- **The fix is behaviorally INERT for the delivered WY 1993–2025 product**: a
  fresh same-machine reference run (`processedOuts_portref_24aug2026`, 6,678 ×
  1,653, 10.9 min) diffed against standard product #1 shows ZERO masking events —
  no qualifying year in that window has a season under 80% completeness (the
  preprocessor's >3-day-gap year rejection makes that nearly impossible). The
  entire observed delta is the documented climate-input residual (product #1
  predates the daymet rebuild): exactly the replay's 98 columns at 1e-18…6.3e-13,
  max 2.56e-13 over 729 of 18.9M annual rows, all climate/SWE-derived. Evidence +
  interpretation retained in the run folder (`RUN_NOTES.md`,
  `validation_additivity_vs_product1.txt`, `validation_annual_rr_delta.txt`).
  Delivered products are NOT superseded by the reference run.
- Windows where seasons CAN fail completeness (other configs, legacy-path users)
  will see seasonal runoff ratios become NaN where flow volumes already did — the
  two families are now consistent.

### Port campaign: rpkg Phases 3–5 implemented + fixture-verified (2026-08-26)
Written and validated WHILE the rpkg Phase-1/2 benchmark ran — source edits and
`devtools::load_all` (which loads into a separate session) never touch the
installed library the running benchmark uses; only `R CMD INSTALL` would.
- **Phase 3 (collector)**: `annual_collector()` / `collector_drain()` /
  `collect_annual()` in `stats.R`, using an environment with pre-grown chunk
  lists rather than `c()` append (O(n²) in R); collection precedes every gate;
  threaded through all 13 module functions and the orchestrator. Contract
  verified: stats identical with/without a collector, values equal the input
  series, integer years, 0 duplicate keys, and rows still collected when the
  metric is below `min_rows` AND below the stats floor.
- **Phase 4 (snow)**: `rpkg/R/snow.R` (14 metrics) + preprocessor SWE support
  in `io.R` (`swe`→`SWE`, per-year policy mirroring PPT, new `valid_swe_years`)
  + record-anchored decade gate + the orchestrator's explicit-`snow_data` rule.
  **Cross-checked vs canonical Julia: 126 values over 3 fixtures, 0 mismatches,
  worst relative difference 3.9e-16.**
- **Phase 5 (drought)**: `rpkg/R/drought.R` (`weibull_quantile`, gap-aware
  `smooth_daily_flow`, duration/deficit at the 5 levels with strict `<`,
  threshold scalars) + config fail-fast + `drought_enabled` gate.
  **Cross-checked vs canonical Julia: 105 values over 3 fixtures (seasonal,
  intermittent, gapped), 0 mismatches, worst relative difference 3.6e-14.**
- **No MK-method change is needed for rpkg** — `Kendall::MannKendall` already
  reproduces the Julia formula exactly, so rpkg was canonical on that point
  throughout; Python was the sole outlier.
- rpkg suite: **1,003 passing / 0 failures** (was 760), including new
  `test-annual_collector.R`, `test-snow_metrics.R` and `test-drought_metrics.R`
  — the last carrying the same Sep 30 → Oct 1 boundary-attribution fixture
  Codex required of the Julia original.

### Port campaign: PYTHON PORT FEATURE-COMPLETE + VALIDATED (2026-08-25)
The Phase-4/5 benchmark (snow + drought, 83 min) produced **6,678 gages ×
1,653 columns — the Julia reference count exactly** — and passed the STRICT
schema gate with NO waivers (column set and gage set both equal to the
reference). Value comparison over 1,620 shared signature columns: mean R²
0.999071, 1,443 Perfect / 121 Good / 55 Poor / 1 Low, and **the Extremely-Low
and Very-Low tiers are EMPTY**. **All 56 columns below R² 0.99 are Pettitt
segment MK p-values** — the single documented methodology difference (see the
Phase-2 entry); zero non-Pettitt columns diverge. The annual parquets share
**18,898,405 of 18,898,406 rows** with identical 100-signature sets, 0 duplicate
keys and **0 NA-pattern mismatches**; 267 of 18.2 M finite pairs differ by
>1e-6 (0.0015 %), all discrete threshold-crossing metrics (the FP-tie class).
- **The one non-shared row is fully explained**: gage 08134000 / avg_storage /
  WY2017 has 365 days but only 9 distinct Q values — including BOTH −0.0 and
  +0.0. Julia's `unique` uses `isequal` (−0.0 ≠ +0.0 → 10 unique → passes its
  `< 10` guard); numpy's uses `==` (→ 9 unique → skips the year). A pure IEEE
  signed-zero semantics difference between the two languages' `unique`,
  1 row in 18.9 M. Like the MK p-value method it is a canonical-side
  convention choice, so it was NOT changed unilaterally — logged for Phase 6.
- rpkg Phase-0 baseline also complete: 516.9 min, 6,678 gages, **0 errored**,
  installed-package testthat **760/760**, and its **gage set is identical to the
  reference** — the rebuilt R runner reproduces Julia's qualification exactly.
  **Found + fixed (pre-existing): all 12 rpkg QA-flag columns were
  DOUBLE-prefixed** (`flagged_for_flagged_for_*`) because the runner
  re-prefixed names `compute_qa_flags` already prefixes — invisible to the
  intersection-based comparison scripts, caught by the strict schema gate.

### Port campaign Phase 5 (2026-08-25, Python drought — implemented)
All 10 drought metrics + 5 threshold scalars ported
(`python/streamflow_signatures/drought.py`): `weibull_quantile` (Hyndman-Fan
type 6 with the project's `below_plotting_range_policy="na"` refusal to
extrapolate), `smooth_daily_flow` (7-day centered, applied within maximal runs
of consecutive dates so it never averages across a gap), duration/deficit at the
five USDM levels with the STRICT `<` comparison, config fail-fast validation
mirroring Julia's, and the `DROUGHT_ENABLED` orchestrator gate (absent config
section => family off). **Cross-language fixture check: 105 values across three
fixtures — seasonal low-flow, intermittent (zero thresholds + strict `<`), and a
GAPPED record that proves the smoothing never blends across a temporal break —
0 mismatches, worst relative difference 6.0e-15.** Tests include the
boundary-attribution fixture Codex demanded of the Julia original (a low-flow
block straddling Sep 30 → Oct 1 must SPLIT across two water years; a
conservation check alone would pass under a uniform year shift).

**Python is now FEATURE-COMPLETE for all six queued features** (Pettitt, stats
floor, annual collector, b=1 alpha, snow, drought). Suite: 130/130.

### Port campaign Phase 4 (2026-08-25, Python snow — implemented)
All 14 snow metrics ported (`python/streamflow_signatures/snow.py`), plus the
preprocessor SWE plumbing (`swe`→`SWE` normalization, per-year SWE policy
mirroring PPT, new `valid_swe_years` return key), the record-anchored decade
gate via the Phase-1b `force_skip_trends` mechanism, orchestrator
`snow_data`/`snow_climate_years` kwargs (explicit frame only — an SWE column in
the main gage frame is NEVER used implicitly, matching Julia), and runner SWE
plumbing. **Cross-language fixture check: 126 values across three discriminating
synthetic gages (clean triangle snowpack, mid-winter thaw producing two spells,
and a vanishing snowpack that trips the decade gate) — 0 mismatches, worst
relative difference 3.5e-16.** Python suite 112/112 (21 new snow tests). Two
hand-derived test expectations were off by one against the code; the code was
correct (it matched Julia exactly), the expectations were corrected — and the
gate fixture was rebuilt so it isolates the record-anchored gate from the
20-value stats floor rather than conflating them.

### Port campaign Phase 3 findings (2026-08-25, Python — in progress)
Annual-values collector ported to Python (`AnnualCollector` + `_collect_annual`,
threaded through all 13 signature functions and the orchestrator — including the
floor-exempt recession/elasticity, as Julia does; collection happens BEFORE every
gate, so the exported rows are exactly what the statistics were computed from).
Runner drains per-gage collectors into a long-format zstd parquet
(`{prefix}_signatures_annual.parquet`, `gage_id/signature/water_year/value`),
written before the metadata merge so a metadata failure cannot lose it.
- **Fixed (schema, BOTH ports): flashiness and FDC used placeholder metric names
  internally (`RB_index`, `slp_all`/`slp_90th`/`slp_mid`) and renamed the stat
  keys afterwards, where canonical Julia uses the final names
  (`flashinessRB`, `FDCall`/`FDC90th`/`FDCmid`) throughout.** Invisible in the
  summary CSV (the rename produced identical column names) but the collector
  captures the column name it is GIVEN — so the annual parquet would have
  carried four wrong signature names in both ports. Found by the collector's
  signature-coverage test. Both ports now use the canonical names directly and
  the rename steps are gone; summary output is unchanged by construction.
  (Exactly the class Codex F6 predicted a self-referential validator would miss.)
Python suite 91/91. Phase-3 benchmark running.

### Port campaign Phase 2 findings (2026-08-25, Python — in progress)
b=1 recession alpha ported to Python (`log_a_pointcloud` / `log_a_events` /
seasonality now use `median(log(-dQ/dt) − log(Q))` with b fixed at 1; `b_*`,
`concavity`, `alpha_linear` keep their free fits). Cross-checked by running
three synthetic gages (linear reservoir, quadratic reservoir, seasonally
alternating alpha) through BOTH languages: **the b=1 values are bit-identical**
(log_a pointcloud/events mean+median exact; `recession_alpha_point_cloud_
linear_reservoir` exact; free-fit `b` within 1 ulp; event counts exact).
- **Fixed (alignment, pre-existing): Python's mid-event day-of-water-year used
  the UPPER-middle index (`indices[len//2]`) where canonical Julia uses the
  FLOOR midpoint (`div(start+end, 2)`)** — so every EVEN-length recession event
  was timed one day late (8 of 12 events/year on the monthly-recession
  fixture). It feeds only the seasonality sinusoid, and was invisible while the
  ports fed that sinusoid free-fit log_a values; the b=1 alignment exposed it.
  rpkg's 1-based `ceiling(n/2)` is algebraically identical to Julia and was
  already correct — **Python only**. After the fix the seasonality outputs
  converge from ~1e-3 relative to **0.0 / 1.7e-15** on the signal-bearing
  fixtures (the zero-amplitude linear-reservoir gage has no seasonal signal, so
  its fitted phase is meaningless in either language). Regression test pins the
  formula and the end-to-end equality.
Python suite 77/77. **Phase-2 exit gates PASS** (89.7 min benchmark, 6,678 ×
1,264): Gate A surgical — exactly 38 log_a-family columns moved vs Phase 1 and
**0 other columns changed bitwise**; Gate B convergence — all 38 now Perfect
(R² ≥ 0.999) against the Julia reference, against −3.13 … −0.28 before.
Overall 6-tier: mean R² 0.999203, **the Extremely-Low and Very-Low tiers are
now EMPTY** (Phase 1 had 39 columns below 0.50); 1,112 Perfect / 85 Good /
33 Poor / 1 Low of 1,231 shared columns; 5 NA-pattern mismatches in total.

**Last remaining divergence class isolated: Mann-Kendall p-value METHOD (full
mechanism established 2026-08-26 — see the diagnosis below, which SUPERSEDES
the first-pass account).** Every sub-0.99 column in the feature-complete port
is an MK p-value; taus, slopes, means, Spearman p-values and changepoint
locations are all unaffected (`_mk_rho` 100/100 Perfect, `_spearman_pval`
100/100 Perfect, `_mk_pval` 56 of 300 below 0.99).

Two distinct mechanisms, the first dominant:
1. **Continuity correction on TIED series.** With ties, scipy's `kendalltau`
   uses its asymptotic branch with `z = S/√var_S` and NO continuity
   correction; Julia (and R) use `z = (S∓1)/√var_S`. Ties are pervasive here —
   day counts (`D*_day`, `TQmean`), pulse counts, drought durations,
   `snow_cover_days`. Measured on integer-valued series: median |Δp| = 0.058
   at n = 10, 0.021 at n = 20, 0.006 at n = 46; significance-flip rate at
   α = 0.05 falls from 1.2 % to 0.1 % over the same range.
2. **Exact vs normal on UNTIED series.** With no ties scipy uses the EXACT
   Kendall distribution regardless of n (verified: scipy 1.18 selects on ties,
   NOT on sample size — correcting the first-pass claim). This effect is small:
   median |Δp| ≈ 0.002–0.004 and **zero** significance flips in 4,000 trials
   per n.

**R is NOT a third variant — `Kendall::MannKendall` reproduces the Julia
formula to 1.000 at every n tested.** Julia is a faithful port of R's
convention, so Julia and rpkg agree and **Python is the sole outlier**.

Real-world impact on the WY 1993–2025 product (Julia vs Python, per cell):
main-path `_mk_pval` 1,298 significance flips in 551,218 cells (**0.24 %**);
Pettitt segment p-values 5,256 in 1,161,770 (**0.45 %**); Spearman p-values
2 in 551,219 (0.0004 %, i.e. not affected).

Options were: (a) accept and document; (b) standardize all three on the exact
method (changes Julia AND rpkg, i.e. the delivered products, and exact is
unavailable under ties anyway — where most of the divergence lives);
(c) change PYTHON only to the continuity-corrected normal formula.

**DECIDED (user, 2026-08-26): option (c) — IMPLEMENTED same day.** Python's
`mann_kendall_test` no longer calls `scipy.stats.kendalltau`; it computes S,
the tie-corrected variance, tau-b and the continuity-corrected normal p-value
directly, mirroring `julia/src/stats.jl` line for line. tau is now computed
here too (algebraically identical to scipy's tau-b for a tie-free time index,
but this removes the last-ulp difference). Verified **BIT-EXACT against Julia
on 12 fixtures** spanning untied/tied, short/long, monotone/constant/reversing,
negative S, single-tie and the n≤3 and constant-series NaN contracts — both
tau and p-value, worst relative difference 0.000e+00. Rationale for the
direction: R's `Kendall::MannKendall` reproduces the Julia formula exactly, so
canonical is corroborated by the convention's origin, and the continuity
correction is the standard treatment for a discrete statistic — scipy's
omission of it under ties is the anomaly. Regression tests pin all 12 fixtures
plus an explicit assertion that the result differs from scipy's uncorrected
asymptotic value. Python suite 143/143.

**VALIDATED (86.6 min benchmark, 6,678 × 1,653, strict schema gate PASS):**
mean R² over 1,620 shared columns rose 0.999071 → **0.999988**, min R²
0.932816 → **0.997935**, and the tier table went 1,443/121/55/1 →
**1,615 Perfect (99.7 %) / 5 Good / 0 Poor / 0 Low** — *no column below
R² 0.99 and no tier below "Good" populated*. Significance flips at α = 0.05
collapsed from 1,298 → **5** on the main path (0.2355 % → 0.0009 %) and
5,256 → **11** on Pettitt segments (0.4524 % → 0.0009 %). The 5 residual Good
columns are Pettitt fields on 3 bases only and are *downstream of changepoint
tie flips*, not the MK formula: `cp_year` agrees exactly on 597,505/597,527
cells (99.9963 %), and a flipped changepoint year moves the segment split.
The annual parquet is **byte-identical** to the pre-fix run
(`DataFrame.equals` True), confirming the change touches statistics only.
**The Python port now has no remaining divergence class.**

### Port campaign Phase 1 findings (2026-08-24, in progress)
Phase 1 (Pettitt + stats floor → Python) surfaced four genuine defects, each caught
by a specific Codex-mandated gate — none would have shown in the old
intersection-only comparisons:
- **Fixed (alignment, MEDIUM): Python and rpkg flashiness guarded `total_Q == 0`
  where canonical Julia guards `<= 0`** — with negative-Q days retained by default,
  a water year with NEGATIVE total flow produced a negative-denominator R-B value
  in the ports that Julia declines. Verified: exactly 4 in-window gage-years across
  3 Florida gages (02236900 ×2, 02244440 WY1999 at −78.8 mm, 02313700). Found by
  the Phase-1 floor-transition EXACT-equality gate (the extra year pushed the
  Python count to 20, dodging the stats floor Julia applied at 19). Both ports
  fixed to `<= 0`.
- **Fixed: Python `qp_seasonality`/`storage` silently dropped the new changepoint
  fields** — they rebuilt their return dicts with explicit 8-key copies instead of
  merging the `generate_stats` output (Julia merges wholesale). Found by the
  schema gate: exactly 3 bases × 8 = 24 missing `_pettitt_*` columns. Now merged
  wholesale.
- **Fixed (tooling): `apply_stats_floor_mask.py` read CSVs with pandas' fast
  (not correctly-rounded) float parser** — rewriting perturbed last bits of
  unmasked values (~1 ulp; the same defect class as the Daymet converter fix),
  which broke exact-equality gates built on its output. Now
  `float_precision="round_trip"`.
- **rpkg testthat suite is green against the INSTALLED package (634/634)** after
  fixing three environment-fragile test defects the `devtools::load_all()` habit
  had concealed (Codex F8 vindicated): a data.table j-expression column-drop that
  returned a character vector under current data.table (test-area_normalized_gate),
  seven unexported-internal calls (`eckhardt_filter` etc. — now `:::`-prefixed),
  and a stale `d_percentiles` length-11 expectation (13 since April 2026). The R
  runner also gained `pkg_env <- streamflowsignatures:::pkg_env` (internal, not
  visible under `library()`).
Interim gate results (pre-fix run): gage sets identical; Pettitt cp_year exact
agreement vs the Julia reference **99.995%** excluding the two known b≠1 log_a
bases (23 single-cell tie flips across rank-sensitive columns, 0 NA-pattern
mismatches); floor transition exact after the mask-tool fix except the
flashiness finding above. Full-gate rerun with the fixes in flight.

### Port campaign Phase 0 (2026-08-24, in progress)
Per `docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md`: Python uv
env (py3.12, 26 unit tests green); rpkg bundled config synced from canonical (29
missing keys; `na_reject_negative_flow` fallback aligned to `false`); **both port
benchmark runners rebuilt** (`run_python_benchmark.py`, `run_rpkg_benchmark.R`):
Julia-mirrored ENV overrides, runtime WY-window + qualifying-fraction gates
(semantics copied from the Julia runner), legacy-filtering fail-fast,
bounded-memory climate handling (column-projected reads, per-gage joins, cache
eviction — the old Python runner merged the full ~98M-row climate frame wholesale),
provenance blocks, zero-error exit gate (R). Two more pre-existing schema gaps
closed at the runner level: rpkg output gains the 11 GAGES-II/HYDAT interference
columns + `human_interference_class` (previously absent entirely), and Python now
fills RHBN/REGULATED from `metadata/canadian_hydat_interference.csv` (the "not
available from Python" comment was stale — the CSV exists for exactly this).
Still open for rpkg: `na_cause_ice` preprocessor tracking behind
`ice_affected_days_total` (Phase 1). New acceptance-gate tool
`docs/benchmarks/check_schema_equality.py` (strict column-set + gage-set equality
with explicit, logged waivers — the compare_* scripts intersect columns and exit 0
on missing schemas, so they are diagnostics, not gates).

### Changed: S3 access LOST — Google Drive is now the off-site backup (2026-08-24)
Access to the project S3 buckets (`climate-ai-data-science-shiny-app-data`, and the
ClimateAI catalog environment generally) has ended. Most delivered artifacts were backed
up beforehand to the project Google Drive folder:
https://drive.google.com/drive/u/0/folders/1DVuq4nC5j_Y01sBaDj9cwjbv7S8sndjj —
"most", so the backup inventory still needs verifying against the delivery manifests.
Consequences:
- Every S3 URL in the docs (EO READMEs, `docs/DEVELOPMENT.md`, `docs/DATA_SOURCES.md`)
  is now a **historical delivery record**, not a live location — pull from the Drive
  backup instead.
- The **NLCD product's pending S3 publish is defunct**: it was awaiting human QA on the
  SageMaker box when access ended. **Resolved same day — the user confirmed the
  finalized NLCD product IS in the Drive backup** (no re-run needed; human QA still
  owed; publication now targets HydroShare). The `temp_lulc_conus/` staging (~93 GB,
  transient by design) is moot.
- The Shiny app read its data live from S3 → **non-functional as deployed**.
- The possible S3 copy of the ORIGINAL `daymet_1980_2023.parquet` (the app read that
  path) is recoverable only if it made the Drive backup — first check of the
  publication plan's staging pass (see the amended Planned item above).
- HydroShare publication planning is unaffected:
  `docs/plans/2026-08-24-hydroshare-publication-and-dashboard-plan.md`.

### Decided: final dashboard scope + performance bar; app work deferred (user, 2026-08-24)
The original Shiny app (`streamflowAndClimateVisualizationApp/`) is judged **very low
performance** — any final published app should be planned as a performant rebuild, not a
port. Required scope for the final app: (i) streamflow signatures and metrics,
(ii) annualized values, (iii) raw time series values, and (iv) the per-watershed MODIS
LAI/LULC and Annual NLCD annual/monthly values and metrics — aggregated tables only,
**no raw LAI or LULC**. Context: CUAHSI has offered hosting (their shinyApps server can
outlive the grant; the app must be Shiny; sunset expectations required) — options,
sunset terms, and sequencing live in the plan doc. **App work is deferred for now**,
decoupled from the Nov 9 submission track.

### Promoted: STANDARD OUTPUT #2 regenerated with the drought family (2026-08-11)
WY 1980–2025 @ 60 % rebuilt on the recovered climate input, superseding
`processedOuts_1980_2025_22jul2026`. Delivered folder:
`/Volumes/Untitled/processedOuts_1980_2025_11aug2026/` (wrapper
`run_julia_benchmark_prod_1980_2025_60pct_drought.jl`; `RUN_NOTES.md` carries the
inventory and gate log). **7.78 min, 6,250 gages × 1,653 columns**; annual parquet
24,366,487 rows / 100 signatures. Both standard products now carry drought, so the
asymmetry noted in the 2026-08-10 entry is closed.

- **Gates**: `validate_production_run.py` 6/7 with the same mis-applied gate-2 waiver
  (1,653 vs the 1,488-column July reference; delta verified as exactly +165 `drought_*`,
  0 dropped). Gage set **identical to the July build** (6,250 = 6,250, 0 either-only,
  lat/lon match). `validate_annual_values.py` **PASS** — 587,829 pairs, 0 mismatches, 0
  coverage misses, 0 orphans, 0 duplicate keys. Year counts [20, 46]; stats floor
  respected; snow 5,635 gages (4,546 US + 1,089 CAN). Comparison vs the July build:
  **1,447/1,451 Perfect**, 4 Good, min R² 0.9977 (FDC), mean R² 0.999992.
- **`check_additivity.jl` FAILS its shared-value check against the July build** — 74
  columns differ materially even at `--tol 1e-6` (worst: `FDC90th_spearman_pval` 0.81,
  `FDC90th_mk_pval` 0.80, `Qwin_*` rank stats ~0.05), plus 437 within tolerance
  (≤ 4.67e-06). That comparison confounds **three** changes at once — machine (Windows →
  M1), climate input (original → rebuilt), and the drought family — and several failing
  columns are **Q-only**, which no climate change can explain. Rather than assert "it's
  the machine", the three were **partitioned with controls**:

  | Comparison | What varies | Result |
  |---|---|---|
  | #1 replay: M1+original vs M1+rebuilt climate | **input only** | 98 cols, ≤ **3.4e-13**, no rank flips |
  | #2 vs its drought-disabled twin (same machine, same input) | **drought code only** | all 1,487 shared cols **bitwise** unchanged |
  | #1: M1+original vs Windows+original | **machine** | **66 cols material**, `FDC90th_spearman_rho` 0.53, `Qwin`/`Qfal` rank stats — *same signature families and order of magnitude as #2's 74* |

  The machine control reproduces the failure signature with **no rebuilt parquet on either
  side**. Two honest limits on that control: it is 66 vs 74 columns and the shared maxima
  are not equal (`FDC90th_spearman_pval` 0.48 there vs 0.81 in #2), so it accounts for the
  *pattern*, not for every column one-for-one; and the July Windows run predates the
  provenance block, so "original climate input on both sides" rests on the run's dated
  record rather than on machine-readable provenance.
  **Mechanism, at the annual level**: comparing the two machines' annual parquets for the
  failing Q-only bases over 200,834 gage-years gives max |diff| **5.68e-14** (`FDC90th`)
  and **1.82e-12** (`Qwin`), with **no rows present on only one side** and no NaN-pattern
  mismatches (independently reproduced by Codex to the same digits). The annual series are
  therefore numerically identical, and every material summary difference sits in a
  **discretely FP-sensitive** statistic — i.e. one where a last-bit change crosses a
  threshold rather than perturbing a value:
  - **rank statistics** (53 of the 67 columns measured this way) — Spearman/Mann-Kendall
    tau and p-values, where a tie flips;
  - **`TQmean`** — a COUNT of days above the annual mean, so one flipped day moves the
    year by 100/365 = 0.274 pp; this propagates to `TQmean_mean` (max 0.149),
    `_median` (0.137), `_linear_slp` (0.026), `_senn_slp` (0.018);
  - **Pettitt fields**, where the changepoint LOCATION jumps discretely
    (`FDC90th_pettitt_cp_year` by 5 years, moving `pct_change` by 26.8);
  - **`FDC90th`** itself — OLS on `log10(Q + 1e-10)` in the near-zero tail, already
    documented as FP-fragile.
  An earlier version of this entry said the divergence appeared "only in rank-based
  summary statistics". **That was wrong** (Codex delta review, 2026-08-11): **14 of the
  material columns are not rank statistics.** The accurate statement is the discreteness
  one above, which matches the 28 Jul cross-machine finding. (Column-count convention:
  `check_additivity.jl` flags **66** material columns — the authoritative figure, retained
  in `validation_machine_control_from_product1.txt`; an independent recompute using a plain
  absolute > 1e-6 rule over finite pairs gives 67. The 14 non-rank columns are the same set
  either way.)
- **✅ ADDITIVITY GATE: PASS** — rather than infer additivity from the WY 1993–2025 window,
  a drought-disabled twin of THIS run was produced on this machine with this climate input
  (`STREAMFLOW_CONFIG` variant with the `drought` section removed; precompile cache purged
  first and restored afterwards, with the `CFG_DROUGHT_ENABLED` probe confirming
  `false` before and `true` after — see DEVELOPMENT.md → config-variant gotcha). Diffed
  with exact equality: **all 1,487 shared columns bitwise unchanged**, gage sets identical
  (6,250), exactly 165 columns added, every added column populated at ≥ 0.75 finite
  fraction. `check_additivity.jl` is a same-machine instrument and cannot pass across
  machines — this is the comparison it was designed for, and product #2 now has the same
  class of proof product #1 has. Evidence retained beside the product as
  `validation_additivity_samemachine.txt`.

### Recovered: Daymet climate parquet rebuilt from the annual CSVs (2026-08-11)
After the canonical `daymet_1980_2023.parquet` was found truncated (Known Issues), the
user restored the **44 annual Daymet CSVs** (1980–2023, 10 GB) to
`/Volumes/Untitled/daymet_1980_2023/`, and the parquet was rebuilt from them with a new
tool, unblocking every climate/SWE-dependent run.

- **New tool `docs/benchmarks/convert_daymet_csvs_to_parquet.py`** — a Python
  reimplementation of the R `convert_daymet_zip_to_parquet()` (there is no R on this Mac)
  using the same positional reconstruction. **Not an exact mirror** — it is stricter (R
  only *warns* on unexpected per-site day counts and never checks per-month counts, input
  columns, year values, invalid dates, or duplicates; those are fatal here) and it keeps
  `site_id` as a string where `fread()` may infer numeric. The CSVs carry
  `site_id, month, year` but **no day column**, so `day` = row order within
  (site_id, year, month) and `Date` = year-month-day. Reads a DIRECTORY rather than a ZIP
  and streams year by year because ~98M rows do not fit in 16 GB; each year is sorted by
  (site_id, Date) and appended, and pyarrow splits by size into **~132 row groups, not one
  per year**. Row order is not semantically meaningful downstream (the runner groups
  climate by gage_id and `preprocess_daily_data()` re-sorts onto a normalized daily grid).
  Output: **97,757,220 rows × 6,087 sites, 4,040,997,608 bytes** — matching the ~98M rows
  / 6,087 sites the guidelines describe for the original.
- **`site_id` is read as a STRING and never coerced.** The streamflow parquet keys on
  zero-padded ids (verified: 6,636 of 8,014 gage ids begin with `0`), so numeric coercion
  would have silently broken the climate join for most of the domain.
- **DATA FINDING — Daymet publishes a 365-day calendar.** Verified across all 44 CSVs:
  leap years keep Feb 29 and instead **drop Dec 31**, so December has 30 rows and every
  year is exactly 365 rows per site. The positional reconstruction therefore yields
  Dec 1–30 in leap years and no Dec 31 row at all — a one-day hole in the climate series
  every leap year, absorbed by the preprocessor as an internal gap ≤ 3 days. This is a
  property of the source reproduced faithfully, **not new**: the R version wrote the same
  rows (its `expected_days` check assumed 366 and merely logged a warning). The converter
  asserts the per-(site, month) day counts under this convention. **That is a structural
  cardinality check, not a proof of chronology**: it catches a month with the wrong number
  of rows (a site short a day in February would otherwise shift every subsequent date
  silently) but **cannot detect a permutation of days *within* a month**. The chronology is
  corroborated end-to-end instead, by the product replay below — a within-month permutation
  could not survive that replay.
- **Validated end-to-end by reproducing a delivered product**, the only check available
  since the original binary is unreadable: the WY 1993–2025 standard config was re-run
  against the rebuilt input and diffed with `check_additivity.jl` — **0 columns added, 0
  dropped, identical 6,678-gage set**, and 98 of 1,653 columns differing by at most
  **3.4e-13** (≤ 69 of 6,678 rows each; most ≤ 1e-15), all climate/SWE-derived. For scale,
  that is ~7 orders of magnitude tighter than the cross-machine FP differences this
  project already accepts between Windows and M1 runs of the same config (≤ 3.2e-06 on
  431 columns).
- **Root cause of the residual, and a fixed defect**: the first rebuild used pandas'
  default CSV float parser, which is ~1 ULP off on ~10 % of these values, and 137 columns
  drifted up to 3.2e-07. The converter now uses `float_precision="round_trip"`, verified
  **bit-exact against a correctly-rounded `strtod` parse on 1,200,000 values (0
  mismatches)**. The remaining 3.4e-13 is **consistent with a last-bit
  ingestion/serialization difference** — R's `fread` uses a fast parser that is not
  guaranteed correctly rounded — but that cause is **not established**: the original
  parquet is unreadable, so it cannot be excluded that the restored CSVs differ from the
  ones that built it, or that some other ingestion detail (blank/NA handling,
  serialization) contributes. What IS established is that the rebuilt values are the
  exactly rounded representation of the restored source CSVs, and that the residual is
  numerically negligible for product use.
- **The truncated file is left in place** at the canonical path (user's direction),
  renamed nothing; runs must set `STREAMFLOW_CLIMATE_PATH` to the rebuilt file. The
  WY 1980–2025 wrapper defaults to it and the choice is recorded in each run's provenance
  block, so a future reader can tell which input any product was built from by its
  recorded size (rebuilt 4,040,997,608 vs original 4,125,630,653).

### Promoted: STANDARD OUTPUT #1 now carries the drought family (2026-08-10)
The 28 Jul drought validation run is **promoted to standard product #1**, superseding
`processedOuts_22jul2026`. Same window (WY 1993–2025 @ 60 %), same config, same 6,678
gages; the difference is the +165-column streamflow drought family (1,488 → 1,653).
Delivered folder: `/Volumes/Untitled/processedOuts_drought_28jul2026` (`RUN_NOTES.md`
carries the full artifact inventory and gate log). The drought family itself was
committed the same day (`53878be`).

- **Artifacts completed to the standard-product convention** (all in the run's own
  folder, per the one-folder rule): signature explorer
  `signature_explorer_1993_2025_60pct_drought.html` (71.8 MB; 1,621 mapped variables =
  100 bases × 16 stats + 21 scalars) + its `_annual/` sidecar folder (100 files, all 10
  drought series present); comparison CSV + summary MD + dashboard
  (`prod_1993_2025_60pct_drought_vs_julia_*`, 84.4 MB) against the superseded 22 Jul
  product.
- **Gates**: `validate_production_run.py` **6/7**, with gate 2 ("column set == reference")
  FAILING. That is a **documented waiver, not a PASS**: the validator was handed a
  reference whose column set was *expected* to differ (1,653 vs 1,488), so the gate is
  mis-applied here rather than satisfied. What compensates is the direct check — the delta
  is **exactly 165 columns, all `drought_*`, 0 dropped** — plus the separate same-machine
  additivity test.
  `validate_annual_values.py` **PASS** (612,046 pairs, 0 mismatches, 0 coverage misses,
  0 orphans, 0 dup keys). `audit_qualification.jl` **PASS** on 11 stratified edge gages
  (near-threshold 20/33 = 0.606 inclusions and 4 exclusions, independent recompute).
- **Comparison vs the superseded product**: 1,451 shared columns, **1,446 Perfect
  (R² ≥ 0.999), 5 Good, 0 below 0.99**, min R² 0.9977 (FDC). ⚠ Cross-machine (Windows →
  M1), so the residual is FP noise in discretely FP-sensitive statistics — this is a
  sanity check, **not** the additivity proof. The additivity proof remains the
  same-machine `check_additivity.jl` diff from 28 Jul (all 1,487 shared columns bitwise
  unchanged).
- **Tooling gap fixed (affects older dashboards too)**: `SIGNATURE_GROUPS` in
  `build_experiment_vs_julia_dashboard.py` gated which signatures become selectable
  targets, and had **never been updated for the snow family** — so every dashboard built
  before 2026-08-10, including standard product #1's and #2's July dashboards, silently
  omitted all 14 snow bases. Added **Snow** and **Streamflow Drought** groups + the 5
  drought threshold scalars to `SINGLE_VALUE_SIGS` (selectable targets for this
  comparison 623 → 820: +112 from the 14 snow bases × 8 stats, +80 from the 10 drought
  bases × 8 stats, +5 drought scalars).
  `compare_experiment_vs_julia.py` and `build_signature_explorer.py` categorize
  generically, so they only needed the two new category labels (drought would otherwise
  have shown as "Other") and the 5 scalars registered in the explorer's `SCALARS`.
- **Standard product #2 (WY 1980–2025) NOT regenerated — blocked.** The run was
  launched the same day and failed in Phase 2: the climate parquet is truncated on the
  drive (see Known Issues). Wrapper
  `run_julia_benchmark_prod_1980_2025_60pct_drought.jl` is committed and ready. Until it
  runs, the two standard products are **asymmetric** — only #1 carries drought — which
  is now stated in the claude-skill's standard-products section.

---

## [July 2026]

### New: streamflow drought signature family (10 metrics, Adelsperger et al. — Julia)
Per-water-year **drought duration** (days below threshold) and **drought deficit**
(summed departures below threshold, mm) at five fixed severity levels — 2/5/10/20/30 %
magnitude percentiles, mirroring the operational U.S. Drought Monitor D4→D0 ladder —
through the standard `generate_stats()` path (8 stats + 8 Pettitt fields each = 160
columns, plus 5 per-gage `drought_threshold_fixed_p{n}` scalar diagnostics →
**+165 columns, 1,488 → 1,653**; 90 → 100 time-series signature bases in the annual
parquet, where the scalars do not appear; annual values exported via the collector).
New module `julia/src/drought.jl` (`calculate_drought_metrics`); full spec + decision
record: `docs/plans/2026-07-27-drought-signatures-plan.md`. Requested by the user
(2026-07-27) from Adelsperger et al. (in review), "A novel severity-based approach for
assessing streamflow drought characteristics and drivers", with four in-session scope
decisions.

- **Method**: 7-day **centered** smoothing of daily Q applied to the CONTINUOUS
  date-indexed series (no artificial Oct 1 discontinuity), within maximal runs of
  consecutive dates so the window never averages across a rejected-year gap (shrinking
  at run edges; < 4 valid days → NaN). Thresholds are whole-record percentiles of the
  smoothed values via the **unbiased Weibull plotting position** `i/(n+1)`
  (Hyndman-Fan definition 6 = `quantile(x, p; alpha=0, beta=0)`) — deliberately NOT the
  type-7 default used elsewhere, which differs most in exactly this low tail.
  Comparison is strict (`<`); ≥ 10 qualifying years required for a threshold.
- **Scope decisions (user, 2026-07-27)**: **fixed thresholds only** — the paper's
  variable day-of-year method is NOT implemented because its per-day sample is one value
  per year, so at 20–46 years of record the low levels carry very large sampling
  uncertainty and the 2 % level falls below the smallest Weibull plotting position
  `1/(n+1)` outright (2 % needs ≥ 49 years, 5 % needs ≥ 19); type 6 would clamp to the
  sample minimum there, and this project's `below_plotting_range_policy: "na"` refuses to
  — an estimability policy layered on type 6, not a property of it; all five levels;
  **water-year** aggregation (documented deviation from the paper's Apr–Mar climate
  year, which splits droughts crossing Oct 1 — day-level indicators and record totals
  are identical either way, annual series and trends are not); NaN rather than
  clamp-to-minimum below the plotting range (unreachable for the fixed method, kept as
  a defensive invariant); explicit `_fixed_` infix in column names so a future variable
  method can be added non-breakingly.
- **Record-dependent**: thresholds come from the run's own window, matching the paper's
  per-analysis-period thresholds — so drought values join `*_all` pulses, elasticity,
  and the parameterized BFI on the "valid within a product, never compared across the
  WY 1993–2025 and WY 1980–2025 products" list. `drought_deficit_*` is also
  unit-carrying (mm only where `area_normalized = TRUE`); the durations are
  scale-invariant and valid everywhere.
- **Gates**: trend completeness and the 20-value stats floor both apply (NOT exempt).
  No snow-style record-anchored gate is needed — a drought-free year emits a valid 0,
  which is precisely the signal a trend should see. **NO existing column changes at all**:
  `flagged_for_high_na` was expected to shift (its denominator counts all signature
  fields — the April 2026 Pettitt / July 2026 snow precedent), but direct per-gage
  recomputation found **zero crossings** (1,224 flagged before and after). Do not
  generalize this to other windows — the closest gage sits only ~6.1e-5 from the 30 %
  threshold, so the outcome is window-specific.
- **Tests**: `julia/test/test_drought_metrics.jl` — hand-derived Weibull quantiles and
  the below-plotting-range guard; exact smoothing values incl. the never-blend-across-a-
  gap and duplicate-date cases; and two record invariants that pin the threshold
  definition AND the day → water-year attribution simultaneously (summed duration ==
  count of pooled values below the threshold == `floor((n+1)p)` for a distinct-valued
  sample; summed deficit == the pooled departure sum). Plus level monotonicity, wet-year
  valid zeros, the intermittent zero-threshold/strict-`<` case, threshold and stats
  floors, schema/collector/zero-warning contracts, and orchestrator + annual-export
  wiring. `smoke_test.jl` gained presence, range, and monotonicity checks.
- **Suite green**: `julia/test/runtests.jl` → **2,042 assertions, 0 failures** (drought
  file: 670). Julia 1.12.6 was installed locally (juliaup) to run it. The expectations
  had first been derived against a purpose-built pure-Python mirror of `drought.jl`;
  the Julia run then confirmed them. **One real registration gap the suite caught**:
  `EXPECTED_DENSE_SIGNATURES` in `test_annual_collector.jl` pins the exact set of dense
  annual series and must list every new dense signature — the 10 drought bases were
  added (that test asserts set EQUALITY, so a new family fails it until registered; now
  called out in the CLAUDE.md / DEVELOPMENT.md "Adding New Signatures" checklists).
- **Smoke test PASSED on real data** (10 gages × 45–46 WYs, feb2026 inputs): all 10
  bases present, ranges valid, level monotonicity holds. Mean durations came out at
  7.296 / 36.499 / 109.542 d/yr for p2 / p10 / p30 against the `p × 365.25` construction
  target of 7.30 / 36.52 / 109.58 — a **weak, near-circular consistency check** (the
  threshold is a percentile of the series it is then counted against), NOT a proof that
  the smoothing or plotting position is right; it catches gross mismatches only.
  `smoke_test.jl` gained ENV path overrides (`STREAMFLOW_DATA_PATH` /
  `STREAMFLOW_CLIMATE_PATH`, the same names the benchmark runner already uses) so it can
  run wherever the inputs are mounted instead of only the Windows `D:\` drive.
- **Codex adversarial review (2026-07-27, codex-cli 0.145.0, read-only): GO-WITH-FIXES.**
  No CRITICAL; 6 MAJOR + 5 MINOR, all resolved — full table in the plan doc §16. The
  substantive ones: (MAJOR-1) the tests' record invariant proved **conservation, not
  attribution** — a uniform water-year shift would have passed it, so a hand-derived
  boundary fixture (a low-flow block straddling Sep 30 → Oct 1 splitting exactly 6 days /
  27.0 mm into WY 2000 and 7 days / 36.0 mm into WY 2001) plus an independent per-year
  recompute were added; (MAJOR-2) the `p × 365.25` check is near-circular and its
  description as "threshold correctness confirmed" was withdrawn; (MAJOR-3) the
  redundancy ratio was statistically invalid (`mean(n) × mean(d) ≠ mean(n × d)`) and the
  non-redundancy claim is withdrawn pending an annual-series comparison at benchmark
  scale; (MAJOR-5) the column count was wrong — **+165 → 1,653**, not +160 → 1,648;
  (MAJOR-6) "unresolvable" overstated the variable-threshold problem. Two code fixes:
  rows with unparseable dates are no longer smoothed with their neighbours, and numeric
  config values (even "centered" windows, `min_valid > window`, out-of-range or unsorted
  percentiles) now fail fast. Post-fix suite: **2,700 assertions, 0 failures** (drought
  670 → 1,328).
- **Benchmark run (2026-07-28, WY 1993–2025 @ 60 %, thumbdrive inputs on the M1):**
  6.05 min, **6,678 gages × 1,653 columns**, annual parquet 18,898,406 rows / 100
  signatures. Gage set identical to standard product #1, confirming the family changes
  no gage qualification. Outputs in their own folder per the one-folder convention:
  `/Volumes/Untitled/processedOuts_drought_28jul2026/` (wrapper
  `docs/benchmarks/run_julia_benchmark_drought_1993_2025_60pct.jl`).
- **REDUNDANCY MEASURED — `drought_duration_fixed_p10` is largely redundant** with the
  existing pulse pair (new tool `docs/benchmarks/analyze_drought_redundancy.jl`, run on
  the ANNUAL series, where `n × mean-duration` reconstructs that year's sub-threshold day
  count — mathematically exact, subject to floating representation). Over 200,834
  gage-years: **r = 0.979**, ρ = 0.971; and judged against the interannual signal rather
  than the series mean (the quantity a trend consumes), **within-gage median r = 0.994**
  (p10 0.971) with disagreement ≈ **11.7 % of each gage's own duration SD** (p90 25.3 %).
  Only p10 collapses: the other four levels measure at r = 0.712 / 0.902 / 0.846 / 0.731
  (p2/p5/p20/p30) against the same pulse pair. Redundancy is an aggregate statement —
  32.5 % of gage-years agree exactly, max disagreement 318 days; on intermittent gages
  99.87 % agree exactly (13 of 9,652 differ, max 8 days). The family's non-redundant
  content is therefore (a) `drought_deficit_*`, the only magnitude-weighted low-flow
  measure in the output, and (b) the four non-p10 severity levels. SIGNATURES.md and the
  claude-skill steer users accordingly. **DECISION (user, 2026-07-28): `drought_duration_
  fixed_p10` is KEPT** — cross-family redundancy is not grounds for removal (err toward
  abundance; the same quantity via two documented methods has independent value, and the
  severity ladder stays complete at all five USDM rungs). The overlap is a documentation
  caveat, not a defect: don't present it and the pulse pair as independent evidence.
- **✅ ADDITIVITY GATE: PASS** (closes Codex MAJOR-4). The rigorous test is a same-machine,
  same-Julia baseline with the drought family disabled (config copy without the `drought`
  section via `STREAMFLOW_CONFIG`), diffed with EXACT equality against the drought-enabled
  run: 165 columns added (expected 165), no column dropped, gage sets identical, **all
  1,487 shared columns bitwise unchanged**, every added column populated. A
  **cross-machine** diff against the delivered product #1 canNOT serve as this gate: 431
  columns differ as pure FP noise (≤ 3.2e-06) and 66 differ materially on 1–38 of 6,678
  rows — all of them discretely FP-sensitive statistics (`TQmean*`, where one flipped day
  is 100/365 = 0.274 pp; `FDC90th*`, OLS on `log10(Q + 1e-10)` in the near-zero tail,
  already documented here as FP-fragile; and rank/p-value/Pettitt fields where a tie
  flips). None attributable to the drought family.
- **New validation tooling**: `docs/benchmarks/check_additivity.jl` — proves a run ADDED
  columns without changing pre-existing ones (column set, gage set, per-gage value
  identity with NaN==NaN, added-column population, `--allow-shift`, `--expect-added`,
  and a `--tol` mode for cross-machine comparisons). Required by the new CLAUDE.md /
  DEVELOPMENT.md checklist step.
- **Fixed (HIGH, config plumbing): GAGES-II directory is now resolved at RUNTIME, not
  precompile time.** Symptom: a benchmark silently emits 1,642 instead of 1,653 columns,
  dropping the 11 GAGES-II/HYDAT human-interference columns with only a warning, because
  `metadata.gages_ii_dir` still points at the Windows `D:/gagesMetadata`. First fix
  attempt added an ENV-derived **constant** `CFG_GAGES_II_DIR` — which does not work: a
  const is baked at precompile time, so a wrapper setting
  `ENV["STREAMFLOW_GAGES_II_DIR"]` at runtime is silently ignored (the same
  precompilation trap documented in `julia/src/config.jl` and DEVELOPMENT.md; it bit
  three times in one day). Now a **function** `gages_ii_dir()` reads ENV at call time and
  is the default argument of `load_gages_ii_interference`; `CFG_GAGES_II_DIR` remains as
  the JSON-only value for compatibility. Verified by probe: in a session where the ENV var
  is set, the constant still returns `D:/gagesMetadata` while the accessor returns the
  mounted path. **Lesson: any ENV override consumed through a `const` is unreliable** —
  this is why the runner already re-reads the window/fraction ENV vars at runtime.
- **Codex review #2 — RESULTS/ANALYSES (2026-07-28): GO-WITH-FIXES, product approved for
  promotion** once the doc errors were corrected; all 4 MAJOR + 7 MINOR resolved (table in
  plan §17). Codex independently re-verified the config equivalence, gage sets, and **all
  16,890,066 shared annual values** (zero mismatches). The findings that mattered:
  (MAJOR-1) my claim of "exactly zero disagreement on intermittent gage-years" was
  **FALSE** — r = 0.981 alongside ρ = 0.492 cannot coexist with identical series; the true
  figure is 99.87 % exact-equal, 13 of 9,652 differing, max 8 days, and the analysis script
  now reports exact-equal %, nonzero count and MAX |diff| so quantiles can't hide it again;
  (MAJOR-2) the redundancy argument scaled a median difference against the series *mean*
  when a trend consumes interannual variation — re-measured within-gage (median r 0.994,
  RMSE/SD 0.117), conclusion survives; (MAJOR-3) `check_additivity.jl`'s population gate
  passed a column with ONE finite value of 6,678, so it could not detect the failure mode
  it existed for — now `--min-finite-frac` (default 0.5) with a by-suffix breakdown, re-run
  at 0.75 = PASS; (MAJOR-4) runs lacked provenance.
- **New: provenance block in the benchmark timing JSON** (`run_julia_benchmark.jl`) —
  resolved paths + size/mtime for every input, SHA-256 for files < 50 MB (config,
  metadata; `STREAMFLOW_HASH_INPUTS=1` also hashes the multi-GB parquets), git revision +
  dirty flag, Julia version, hostname, and the experiment ENV overrides. Additive to the
  timing JSON only — no signature output changes. **The 28 Jul product predates it**; a
  6-min re-run would capture it.
- Also documented from this round: `flagged_for_high_na` stability is window-specific (the
  closest gage sits 6.1e-5 from the 30 % threshold — do not generalize); the ~20 % NaN rate
  in some drought columns is two separate mechanisms (495 gages / 7.4 % lose trend stats to
  the completeness gate — the *same* set that loses `Qann_senn_slp`, so not drought-specific
  — plus ~300 losing rank stats to constant series); and an audit confirming no past
  config-variant result is invalidated by the precompilation gotcha (DEVELOPMENT.md).
- Julia only; Python/rpkg ports deferred (now a four-feature port queue).
  `validate_production_run.py` signature-count gate updated 90 → 100.

### New: per-watershed Annual NLCD product (CONUS land cover + impervious, non-signature)
A CONUS-only US companion to the continental MODIS LULC product, built on **USGS/MRLC Annual NLCD
Collection 1 (C1V2, 30 m, 1985–2025)**: per-watershed annual land-cover % (16 classes) + basin-mean
fractional impervious surface for **6,119 CONUS gages × 41 years = 250,879 rows**, joining to the
signatures by leading-zero-safe `gage_id`/`canon_id`. Python (metadata/ingestion — like HydroATLAS
and the MODIS EO layers, not a cross-language signature). New modules under `EO_data_processing/`;
full design + review record: [EO_data_processing/README_NLCD.md](EO_data_processing/README_NLCD.md).

- **Source/access** — Annual NLCD COGs on `s3://usgs-landcover` (us-west-2, requester-pays); native
  CRS is **WGS84-Albers (no EPSG code — read from `src.crs`, not the legacy EPSG:5070)**; fill=250.
  Source mosaics staged once to `s3://…/streamflow/temp_lulc_conus/` (82 files, 93.4 GB, us-east-2
  in-region; **temporary — delete after QA**). Fractional-impervious carries a documented
  out-of-range/underflow bug — guarded by a [0,100] clip (0 occurrences in C1V2).
- **Pipeline** (`eo_processing/nlcd_pipeline.py`) — per-year-parallel download→extract→delete
  (scratch on `/tmp`, NOT the 2.6 GB home volume), coverage-weighted exactextract (basins reprojected
  to the granule CRS; never resample categorical); publish-integrity gates (per-year validation,
  atomic writes, final table written only if every requested year passes, else nonzero exit — no
  silent partial publish). ~10 h for 41 yr on 4 cores. Class %s sum to exactly 100 on covered rows.
- **Finalize** (`eo_processing/nlcd_finalize.py`) — drops the 45 out-of-footprint **Alaska** gages
  (never published as real zeros) to a traceable exclusion CSV + all-gage QA companion; separate
  documented QA flags (`geom_low_confidence`/`low_pixel_support`/`partial_coverage`/`low_confidence`);
  provenance (`nlcd_collection=C1V2`, `valid_area_km2`, metadata JSON); data dictionary with
  model-noise / annual-update-seam / developed-vs-impervious-non-additive caveats. **S3 upload
  opt-in (`--upload`), held for QA.**
- **QA/QC explorer** (`viz/build_nlcd_explorer.py` → `nlcd_explorer.html`, 12.6 MB) — self-contained
  Leaflet map, 6,119 points × 32 variables; click a gage → stacked-area composition 1985–2025 +
  dashed impervious line, 2024–25 update-seam shaded; smoothed endpoint-diff deltas (labeled),
  `annual_volatility`/`shrub_grass_swap` artifact views, partial-coverage badges, excluded-45
  disclosure, normalized Shannon.
- **Codex adversarial reviews** — pipeline (FIX-FIRST → all publish-integrity blockers fixed) and
  results+plans (GO-with-fixes → all reconciled). Validation: dev%↔impervious% corr 0.963; documented
  shrub↔grass annual artifact confirmed (−0.71 swap correlation, 77–90 pp single-year swings).
- **Known gotcha (this SageMaker box):** `/home/sagemaker-user` has only ~2.6 GB free — heavy temp
  I/O must go to `/tmp`; the pipeline stages downloads there via `--scratch` + a disk preflight.
- Python only (CONUS/US); complements — does not replace — the continental MODIS LULC. **Remaining:
  human QA → S3 publish → delete `temp_lulc_conus/`.**

### New: STANDARD OUTPUT #2 — "entire period of record" WY 1980–2025 @ 60% (2026-07-22, Codex results-review GO, zero findings)
Second standard product (user decision: entire record, operationalized as
WY 1980–2025). Wrapper `docs/benchmarks/run_julia_benchmark_prod_1980_2025_60pct.jl`;
plan + full gate log: `docs/plans/production_run_1980_2025_60pct_plan.md`. Outputs
(own experiment folder per the one-folder convention):
`D:/processedOuts_1980_2025_22jul2026/streamflow_1980_2025_60pct_22jul2026_*` +
signature explorer (61.6 MB, 1,456 mapped variables) + comparison dashboard/CSV/
summary vs the April full canonical (written into the run folder via the new
`--output-dir` defaults).

- **25.1 min, 6,250 gages × 1,488 columns; annual parquet 21.82M rows.** All
  gates PASS/attributed: column delta vs canonical = exactly the 224 snow
  columns; all 6,250 gages ⊆ the canonical 7,313; independent end-cap audit
  PASS (both denominator regimes verified — gage-capped for early-ending,
  1980-anchored 46-year for late starters); annual values 525,329 pairs with 0
  mismatches/0 dup keys; QA flags clean (0 percentile/timing/seasonal-sum
  violations).
- **Cleanest cross-validation to date**: vs the April FULL canonical,
  1,189/1,227 columns Perfect at min R² = 1.0000 in every non-recession
  category INCLUDING Pettitt fields; all 38 divergent columns are the
  intentional b=1 log_a family (0 NA mismatches). Root cause of the
  cleanliness: the source parquet begins at 1980, so the canonical's uncapped
  window ≈ this run's window — run #2 is effectively the canonical analysis,
  properly capped + fraction-filtered + July features.
- **Run #1 ↔ #2 relationship** (by design of the window-start-anchored 60%
  denominator): 5,771 shared gages; 479 only-in-#2 (1980–1992 starts whose
  records ended before #1's window); 907 only-in-#1 (late starters at
  0.435–0.587 < 0.60 over the 46-year anchor — every one verified). Neither
  standard is a subset of the other. `area_normalized = FALSE`: 28 rows here /
  32 in #1 / 37 full-record (window-dependent).
- Codex results review: **GO with zero findings**; all probes confirmed
  independently (incl. the `season_excluded_years_*` −1-only pattern at 4,404
  gages vs the uncapped canonical).
- **Cross-product consistency verified (2026-07-23)**: run #1's annual parquet
  is contained in run #2's for all 5,771 shared gages (0 missing rows; the
  1.39M absent rows are exactly the 907 run-#1-only gages), with overlap-year
  values bit-identical for within-year-computable signatures. Record-dependent
  signatures (`*_all` pulses, elasticity, parameterized BFI) differ between
  windows by design — do not mix them across products. Interpretation guidance
  added to the claude-skill ("The Two Standard Output Products"), incl. why the
  longer window has fewer gages (window-start-anchored 60% denominator: 60% ×
  46 ≈ 28-year effective floor vs the 20-year floor binding in #1).

### New: STANDARD OUTPUT #1 — production run WY 1993–2025 @ 60% (2026-07-22, Codex results-review GO)
First of the two standard products (user decision; HISSS manuscript §2.2.2).
**First production exercise of the `STREAMFLOW_END_WATER_YEAR` end cap** — the
July 14 run was WY ≥ 1993 UNCAPPED (correcting the 2026-07-21 reconciliation-log
claim of a 1993–2022 cap). Wrapper
`docs/benchmarks/run_julia_benchmark_prod_1993_2025_60pct.jl`; plan + full gate
log: `docs/plans/production_run_1993_2025_60pct_plan.md`. Outputs:
`D:/processedOuts_22jul2026/streamflow_1993_2025_60pct_22jul2026_{signatures.csv,
signatures_annual.parquet, timing.json}`.

- **27.6 min, 6,678 gages × 1,488 columns; annual parquet 16.89M rows / 90
  signatures.** Gage set = July's 6,579 + 99: the cap bounds the qualifying
  denominator at 33 (partial-WY2026 parquet rows inflated uncapped denominators),
  so inclusion can only grow — all 99 verified at exactly 20/33 = 0.606.
- First run carrying the **60% overall trend gate** (gained 45 / lost 0 trend
  gages vs April on shared gages, all in the [0.60, 0.80) completeness band) and
  the **snow record-anchored decade gate** (876 snow_on_dowy / 645 melt_com_dowy
  gages trend-suppressed with means + Pettitt intact; exempt swe_max only 333).
- **Gates all PASS/attributed**: independent end-cap window audit
  (`audit_qualification.jl`, extended with an end-year arg) PASS on stratified
  edge gages; annual-values consistency PASS (545,266 pairs, 0 mismatches, 0
  floor violations, 0 dup keys); production gates 6/7 with the strict
  column-equality "FAIL" fully attributed (reference = April pre-Pettitt/pre-snow
  experiment; delta = exactly 608 Pettitt + 224 snow columns — the July products
  were lost to the flash-drive rollback, so April is the only surviving
  reference). Codex results review **GO** with 2 attributed MINORs (32-of-37
  un-normalized gages qualify under this window — the 37 count is full-record;
  2 gages lost BFI_Eckhardt mean/median to the 20-value stats floor at 17
  annual values). Dense shared means/medians EXACT vs April (max abs diff 0.0).
- **QA/QC dashboards** (same day, all in the run folder): signature explorer
  (rebuilt with the Statistic picker extended 8 → **16 stats** — the 8 Pettitt
  fields added to `build_signature_explorer.py`, so all 1,456
  signature-statistic combinations are reviewable); primary QA comparison vs the
  same-window April experiment (597/619 Perfect at min R² 1.0000, all 22
  divergent columns in the comparator's signature scope = the intentional b=1
  log_a family, 0 NA mismatches); window-sensitivity comparison vs the April
  full-record canonical.
- **Codex dashboard/results review (initial NO-GO → resolved, GO to commit)**:
  (MAJOR) the comparison DASHBOARD visualizes 4 `season_excluded_years_*`
  diagnostics the comparator excludes as metadata — their divergence vs April
  verified benign (uniformly −1 at 5,038 gages: the uncapped April run counted
  phantom partial-WY2026 seasons as excluded; the capped run correctly doesn't),
  and the scope difference is now documented in both tools' source. (MINORs
  fixed) both comparison scripts gained `--output-dir` defaulting to the
  experiment CSV's folder (one-folder convention now self-enforcing);
  `audit_qualification.jl` edge-sample slices clamped (`first(v, n)`); `temp/`
  gitignored. Explorer HTML size (64.6 MB) noted as a watch item.
- **NEW CONVENTION (user)**: all outputs of an experiment — run artifacts,
  explorer + sidecars, comparison dashboards/CSVs/summaries — live together in
  the experiment's own unique folder (this run: `D:/processedOuts_22jul2026`);
  `docs/benchmarks/` keeps only tools + long-lived reference CSVs. Recorded in
  CLAUDE.md → Critical Constraints.

### New: record-anchored decade gate for threshold-dependent snow metrics (Julia, 2026-07-22)
The 10 timing/melt/regime snow metrics (`swe_max_dowy`, `snow_on_dowy`,
`snow_off_dowy`, `melt_season_days`, `melt_rate`, `ssm`, `melt_before_peak` +
`_pct` + `_to_max_swe`, `melt_com_dowy`) now additionally require ≥
`decade_min_fraction` of the SWE-valid years in BOTH the first and last decade of
the gage's SWE record to be computable (snowy) before trend statistics are
computed. Requested by the user (2026-07-22) with four in-session decisions:
threshold **LINKED to the streamflow decade gate** (same
`na_handling.trend_completeness.decade_min_fraction` knob — no new constant, one
config line governs both); anchored to the **SWE-valid record** (rows where
`swe_max` is non-NaN — dense incl. zeros), not the metric's own span; scope =
the 10 threshold-dependent metrics only (the 4 magnitude metrics emit valid zeros
whose dense series legitimately carry snow decline); failure NaNs the **6 trend
stats only** (mean/median + Pettitt fields survive, matching streamflow trend-gate
semantics). Motivation: the own-span trend gate and the 20-value stats floor both
pass a gage whose snowy years are clustered (e.g. snow present 1981–2000, absent
2001–2010 → own span 100% complete), yielding trends conditioned on snow-present
years. Mechanism: new opt-in `force_skip_trends` kwarg on `generate_stats()`
reusing the existing skip-trends path; `calculate_snow_metrics()` computes the
gate and passes failing metrics. Denominators count SWE-valid years per window
(missing Daymet years don't count against snow presence); record span < 10 →
gate skipped; config flag `snow.record_anchored_decade_gate` (shipped `true`;
absent → disabled); collector/annual parquet unaffected (collection precedes
gating). Tests: `julia/test/test_snow_record_decade_gate.jl` (96 assertions:
kwarg mechanics, vanishing/appearing snow, valid-year-denominator pin at the
dc=0.7 separation point, linked-threshold, short-record skip, collector
invariance, magnitude exemption); full suite green (1,302). Plan:
`docs/plans/2026-07-22-snow-record-anchored-decade-gate.md`. Julia only (snow
metrics are Julia-only; gate ports with the queued snow port — note rpkg's
bundled config has no `snow` section yet, so nothing to mirror there). Output
impact at next benchmark: gated gages flip 6 trend stats → NaN per gated metric;
`flagged_for_high_na` denominators shift (accepted precedent).
**Codex adversarial review (2026-07-22): GO** — no CRITICAL/MAJOR. Confirmed the
`stats.jl` restructure is byte-identical legacy behavior across the full branch
matrix (collector/min_rows/stats-floor/trend-gate/changepoint ordering), the
anchor-set and window arithmetic (incl. span-exactly-10 and `7/10 < 0.7 ==
false` boundary), the linked-knob production path (config → benchmark →
orchestrator → snow gate), the denominator-separation test's validity, and docs
consistency. One MINOR fixed: stale rpkg-mirror step in the plan doc (marked as
executed no-op).

### Changed: trend-completeness overall gate lowered 80% → 60% (all languages, 2026-07-21)
`config/signatures_config.json` (and rpkg's bundled copy `rpkg/inst/config/`) →
`trend_completeness.min_fraction` 0.80 → 0.60; `decade_min_fraction` stays 0.80.
Source: the July 2026 guidelines doc revision AND HISSS manuscript §2.2.3, which
both state "at least 60% of the entire annual streamflow metric time series must be
complete" (first/last-decade requirements remain 80%). Config-only change: Julia
(`CFG_NA_TREND_MIN_FRACTION`), Python (`NA_TREND_MIN_FRACTION`), rpkg
(`pkg_env$na_trend_min_fraction`), and legacy R (`NA_TREND_MIN_FRACTION`) all read
the JSON at load — verified loaded values 0.60/0.80 in Julia and Python. Language-side
fallbacks stay 0.80 (legacy behavior when the config section is absent, matching the
stats-floor convention). Effect: trend statistics (slopes, correlations, p-values)
populate for gages with 60–80% complete annual series that were previously NaN'd;
mean/median/Pettitt fields and the recession/elasticity exemptions are unaffected.
Takes numerical effect at the next benchmark; the planned WY 1993–2025 and
WY 1980–2025 standard runs (see Planned) include it. Docs updated: SIGNATURES.md,
CROSS_LANGUAGE_STATUS.md, claude-skill.

**Codex adversarial review (2026-07-21)**: confirmed config pickup in all four
consumers, the safety of the added JSON `comment` key, and every Manuscript
Reconciliation Log claim (gap rejection unconditional vs config-driven negative-Q/
constant-SD flags; drainage-area normalization; gate affects trend stats only;
recession/elasticity exemption). Two fixes from its findings, both verified by
execution: (1) the legacy R trend-completeness test's "70% complete" case was
VACUOUS — it used 9 trailing NAs, but completeness is measured over the metric's
own non-NA year span, so the case could never trigger (pre-existing test bug,
its assert fails when run; rewritten as an internal-gap 66.7% case that does
trigger); (2) stale comment in `R/helperFunctions.R` claiming Python/Julia
"always pass 0.80" (they pass config-derived values). Julia mechanism tests that
pass explicit 0.8 thresholds (`test_changepoint.jl`, `test_snow_metrics.jl`) are
intentionally parameterized and unchanged. Also verified by execution: legacy R
loads 0.60/0.80 from the shared JSON.

Second (full) Codex review pass: no CRITICAL; 2 MAJOR fixed — the two non-archive
plan docs (`production_rerun_1993_2022_60pct_plan.md`, `snow_signatures_plan.md`)
still described the 80% gate as current and could mislead future reruns; both now
carry dated update notes (historical statements preserved as accurate for the runs
as executed). Notes from its MINOR findings: an **already-installed rpkg build
stays at 0.80 until reinstalled** — the package loads its bundled
`inst/config/signatures_config.json` via `system.file()`, so reinstall rpkg before
any future rpkg benchmark; historical 80% mentions in DEVELOPMENT.md /
CROSS_LANGUAGE_STATUS.md comparison tables are intentionally unchanged (they
describe April 2026 runs). Rerun of the legacy R NA-handling test suite: the
rewritten trend-completeness cases all PASS; 6 unrelated pre-existing failures
surfaced (see Known Issues).

### New: stats floor — min 20 annual values before ANY statistics (Julia)
`min_values_for_stats = 20` (config `stats_floor`): a metric with fewer non-NA annual
values emits NaN for ALL 8 statistics AND its Pettitt fields. Requested by the user
after gage 07292500's `snow_on_dowy` carried a Theil-Sen slope built on 4 clustered
years (1982–85) — the 80% trend gate measures the metric's OWN year span (4/4 = 100%),
the decade check skips spans <10, and `generate_stats` admitted ≥3 values. The floor
harmonizes with Pettitt's existing `min_total_obs=20`. **Recession and elasticity are
exempt** (inherently sparse — same exemption as trend_completeness; enforced at the
orchestration layer: the kwarg is never passed to them). Collector/annual-parquet
export is unaffected (collection precedes gating). Threaded through `generate_stats`
+ orchestrator + all non-exempt signature functions; absent config section = no floor.
Tests: `julia/test/test_stats_floor.jl` (98 assertions incl. the 07292500-shaped
regression and exemption pins). Codex review: NO-GO (2 findings: post-hoc mask left QA
flags stale → new `refresh_qa_flags.jl` recomputes all 12 via canonical
`compute_qa_flags`; legacy-path scope documented) → delta-verified GO. Post-hoc
tooling for already-delivered CSVs: `docs/benchmarks/apply_stats_floor_mask.py` (+
floor-aware `validate_annual_values.py --floor/--floor-exempt`). Ports deferred
(four-feature queue).

### New: production rerun harness + validation tooling + signature explorer
Supporting the July production reruns (WY1993+/60% and WY1980+/60%, both validated
9-gates-PASS with Codex results-review GO; plan + full execution log:
`docs/plans/production_rerun_1993_2022_60pct_plan.md`):
- **Benchmark runner**: ENV overrides for output dir (`STREAMFLOW_OUTPUT_DIR`) and
  input paths (`STREAMFLOW_DATA_PATH/_CLIMATE_PATH/_METADATA_PATH`); end-water-year
  window support (`STREAMFLOW_END_WATER_YEAR`, non-legacy path); memory patches for
  the 16 GB machine (climate frame trimmed to used columns before the water-year
  copy; raw frames released after Phase 3; per-gage preprocess-cache eviction) after
  the WY1980 run thrashed at 20.8 GB commit. Production wrappers
  `run_julia_benchmark_prod_{1993,1980}_60pct.jl`.
- **Repo-resident validation** (after on-drive copies were lost to a flash-drive
  rollback): `validate_production_run.py` (window/columns/snow/parquet/year-count/
  floor gates), `audit_qualification.jl` (independent window+fraction inclusion
  audit), floor-aware `validate_annual_values.py` (+ fixed a latent itertuples
  `_merge` crash found on its first real-data run).
- **Signature explorer** `build_signature_explorer.py`: self-contained Leaflet map of
  any run's CSV — 90 bases × 8 stats + scalars via two pickers, robust p5–p95 color
  scaling, viridis-colored distribution histogram, per-gage click panel, and annual
  time-series plots (Sen trend line + Pettitt marker) lazily loaded from
  per-signature sidecar files built from the annual parquet.
- **Finding — Daymet covers Canadian gages** (~1,100 Canadian gages carry snow
  values); earlier "US-only Daymet" doc claims corrected in SIGNATURES.md §12 and the
  claude-skill.

### New: Q-to-PPT unit gate for un-normalized gages (Julia + Python + rpkg)
Gages with no drainage area in HYDAT keep flow in raw m³/s (decision: no backfill).
`calculate_all_signatures()` gains an `area_normalized = true` kwarg in all three
languages: when `false`, ALL Q-to-PPT signatures — runoff ratios (+
`runoff_ratio_high_count`), elasticity (+ diagnostics), Q-P seasonality, and
`avg_storage` — are skipped, because Q (m³/s) and PPT (mm) units don't match, making
Q/P, dQ/dP, and cumsum(P − Q) meaningless. Q-only signatures (incl. Q##/D##
percentiles, recession, BFI) are unaffected. All three benchmark runners look up
`area_normalized` from the metadata CSV (leading-zero-safe join; only an explicit
FALSE gates, missing defaults to normalized) and pass it per gage.

- **Intentional output change at the next benchmark**: 16 of the 37 un-normalized
  gages currently carry mixed-unit climate signatures in the April 2026 golden
  (e.g. runoff ratios computed from m³/s Q against mm PPT) — those become NA.
- Elasticity and `qp_bimodality` are technically scale-invariant in Q, but are gated
  with the rest: mixed-unit inputs are conceptually invalid, and the family's
  internal PPT thresholds and P − Q terms are not scale-invariant.
- Tests (all three languages verify: 84 climate keys normally → 0 gated; gated
  output IDENTICAL to a no-climate run; gate inert when true):
  `julia/test/test_area_normalized_gate.jl` (also covers the AnnualCollector path),
  `python/tests/test_area_normalized_gate.py`,
  `rpkg/tests/testthat/test-area_normalized_gate.R`.
- **Codex adversarial review (2026-07-14)**: 2 MEDIUM findings, both fixed — (1) the
  Julia/Python runners hard-failed on metadata files predating the `area_normalized`
  column (now header-guarded, matching rpkg's existing guard); (2) explicit-FALSE
  parsing diverged across languages for 0/1-serialized booleans (now harmonized in
  all three runners: Bool false / numeric 0 / "FALSE"/"0" strings gate; anything
  else, incl. missing, defaults to normalized). Fix verification surfaced a third
  latent bug the review missed: the Julia runner's join-id helper was typed
  `::String` but CSV.jl yields InlineStrings — Phase 2b would have crashed on the
  production metadata (loosened to `AbstractString`). Runner parsing was then
  functionally verified in all 3 languages against TRUE/FALSE, 0/1, and
  missing-column metadata variants, plus the production CSV (16,994 gages parsed,
  exactly the expected 1,601 gated).

### New: snow metrics signature family (14 metrics, Daymet SWE — Julia)
Fourteen per-water-year snow metrics from daily SWE, through the standard
`generate_stats()` path (8 stats + 8 Pettitt fields each → **+224 columns**; 76 → 90
time-series signatures; annual values exported automatically via the collector). New
module `julia/src/snow.jl` (`calculate_snow_metrics`); full spec + review record:
`docs/plans/snow_signatures_plan.md`. Requested by the user with four in-session domain
decisions (2026-07-14).

- **Metrics**: `swe_max`, `swe_max_dowy`, `snow_cover_days`, `snow_on_dowy`,
  `snow_off_dowy`, `melt_season_days`, `melt_rate`, `ssm` (Hatchett 2021, Hydrology
  8(1):32), `swe_apr1`, `melt_before_peak`, `melt_before_peak_pct`,
  `melt_before_peak_to_max_swe`, `melt_com_dowy`, `swe_max_to_ppt`.
- **Domain decisions**: SWE ≥ 10 mm day threshold applied to durations AND magnitudes
  (thresholded series SWE*; sub-threshold years are operationally snow-free);
  snow-on/off anchored to the spell containing the annual peak (boundary-censored →
  NaN); melt rate = max SWE ÷ melt-season days; SSM seasonal spells ≥ 60 continuous
  days. Config: `signatures_config.json` → `snow` + `na_handling.snow_na_policy`.
- **Plumbing**: Daymet `swe` was previously dropped at the climate join — now
  normalized (`swe`→`SWE`) and carried through; `preprocess_daily_data()` gained a
  per-year SWE policy (mirrors PPT: >30 NAs, ≤3-day interpolation, negative rejection)
  and a new `valid_swe_years` return field. Snow metrics run ONLY on an explicitly
  passed `snow_data` frame filtered to those years — NO implicit fallback to the gage
  frame (plan review finding: prevents SWE-invalid years leaking in). Canadian gages
  (no Daymet) → all snow columns NA.
- **Known existing-column change (accepted)**: `flagged_for_high_na` counts the 224 new
  NA fields in its denominator, so some no-SWE gages will flip at the re-run
  (April 2026 Pettitt precedent; whitelisted + quantified; semantics pinned by test).
- **Two-round Codex plan review**: initial NO-GO (2 CRITICAL: orchestrator fallback
  leak; the false "no golden divergence" claim) → 7 findings + 2 delta-round residuals
  incorporated → GO. Record: plan doc §12.
- **Tests**: `julia/test/test_snow_metrics.jl` (~330 assertions: exact hand-derived
  gages — triangle, mid-winter thaw, spell arithmetic, threshold splitting, boundary
  censoring, tie rule, leap-year April 1, PPT gating incl. legacy runoff-ratio parity,
  preprocessor `valid_swe_years`, no-fallback + 0-row schema contract, SWE-invalid-year
  leak regression, collector/zero-warn invariants, QA-flag semantics pin);
  `smoke_test.jl` now joins SWE and validates snow presence/ranges.
- Julia only; Python/rpkg ports deferred (three-feature port queue). Benchmark re-run
  pending — bundled with annual values + b=1 alpha.

### Docs: data-source inventory
New `docs/DATA_SOURCES.md` — table of all 11 external data sources (data class /
project use / provider / original access): USGS NWIS, HYDAT, Caravan, Daymet,
GAGES-II boundaries + attribute tables, ECCC MDA_ADP basin polygons, HydroBASINS,
HydroATLAS/BasinATLAS v10, MODIS MCD15A3H (LAI) and MCD12Q1 (LULC). Also records
the guidelines Google Doc as a governance input and the deliberately excluded
boundary cases (planned ERA5/PRISM, S3 mirrors, viz-only CDNs, unused geometry
backstops).

### Changed: recession alpha now assumes a linear reservoir (b = 1) — Julia
All alpha outputs (`log_a_pointcloud`, `log_a_events`, the 6 `log_a_seasonality_*`
scalars) are now computed with the recession exponent **fixed at b = 1** across all
locations and periods: `log(a) = median(log(-dQ/dt) - log(Q))`, no regression. The
exponent analyses are unchanged — `b_pointcloud`, `b_events`, and `concavity` keep
their free power-law fits, as do event identification, the 25-event gate, and
`alpha_linear` / `recession_alpha_point_cloud_linear_reservoir` (already b=1).

- **Rationale** (domain decision): log(a) is the intercept of a regression whose slope
  is b, so free-fit alpha estimates are convolved with b — trends/seasonality in alpha
  partly reflected changes in b. Fixing b = 1 decouples them.
- Column names unchanged (frozen CSV contract) — methodology change behind existing
  columns, like the April 2026 recession rewrite. At the next benchmark, log_a columns
  will diverge from golden intentionally; b/concavity/alpha_linear must not.
- Consequence: under the forward-difference discretization each b=1 point equals
  `log(1 − Q_{i+1}/Q_i)`, so per-year `log_a_pointcloud` APPROXIMATES
  `log(1 − alpha_linear)`. The relation is approximate, not a same-sample identity:
  `alpha_linear` also includes events whose free power-law fit failed (the point cloud
  only pools fit-success events), and medians only commute exactly with the monotone
  transform at odd pair counts.
- New tests (`julia/test/test_recession_alpha_b1.jl`): exact helper values on a pure
  exponential; a quadratic-reservoir gage (−dQ = k·Q² exactly) proving b stays 2 under
  the free fit while alpha matches an independent b=1 recompute (and is far from the
  old free-fit intercept); linear-reservoir continuity (b=1 gage values unchanged).
- Julia only; Python/rpkg sync deferred (bundled with the annual-values port).
  Design note: `docs/plans/recession_alpha_linear_reservoir_plan.md`.

### New: annual values export (per-year signature values, Julia)
The per-year annualized metric values — previously computed inside every signature
function and discarded after `generate_stats()` collapsed them into the 8 summary
statistics — are now saved as one long-format parquet alongside the summary CSV:
`{prefix}_signatures_annual.parquet` with columns `gage_id, signature, water_year, value`
(~20M rows expected for 7,313 gages × ~76 signatures, ~170 MB zstd).

- **Mechanism**: opt-in `AnnualCollector` threaded through `calculate_all_signatures()`
  and all 13 signature functions into `generate_stats()` (mirrors the `changepoint`
  kwarg plumbing). Collection happens BEFORE min_rows / trend-completeness gating and is
  read-only w.r.t. the statistics — the summary CSV contract is untouched (verified by a
  with/without-collector equality test).
- **Semantics**: the export records the exact series `generate_stats()` consumed. Dense
  signatures carry NaN placeholder rows ("year qualified, metric not computable");
  caller-pruned signatures (`qp_bimodality`, `n_recession_events`, flashiness, storage,
  baseflow, elasticity) omit those rows — absent row ≡ NaN for consumers.
  `elasticity_rolling` is keyed to the END year of each 11-yr window; `elasticity_annual`
  to the later year of each consecutive pair.
- **Config**: `config/signatures_config.json` → `annual_values.save` (shipped `true`;
  defaults to `false` when the section is absent). Constant `CFG_SAVE_ANNUAL_VALUES`.
- **Validation**: new tests in `julia/test/test_annual_collector.jl` (stat-identity,
  zero-warnings, signature coverage, NaN/missing-year guards, no duplicate keys) +
  `docs/benchmarks/validate_annual_values.py` (recomputes mean/median from the parquet
  and cross-checks the summary CSV; runs after the next benchmark).
- **Pre-flight**: Parquet2 pinned (`0.2`, resolved v0.2.27) after a write→read round-trip
  check (zstd/snappy/uncompressed, NaN-safe, empty-frame-safe; 1M rows = 8.5 MB, 0.2s).
- **Two Codex reviews**: plan-stage (NO-GO → 5 amendments incorporated pre-implementation)
  and code-stage (core implementation confirmed correct; 3 harness findings fixed:
  validator ≥3-non-NaN-rows coverage rule, zero-gage no-op exit, and a deterministic
  exponential-recession test gage covering the recession + parameterized-BFI collector
  paths the noisy synthetic gage never exercised). Full records:
  `docs/plans/annual_values_export_plan.md` §9–§10.
- Julia only — Python/rpkg ports deferred.
- **Benchmark re-run pending** — bundled with an upcoming workflow feature.

## [June 2026]

### New: per-watershed MODIS Earth Observation products (LAI + LULC, non-signature)
Two new per-gage EO products aggregated over each gage's upstream watershed from MODIS Collection 061, sitting alongside the signatures and joining by leading-zero-safe `gage_id`/`canon_id`. Python (not a cross-language signature — like HydroATLAS, this is metadata/ingestion). 7,964 watersheds. New subproject `EO_data_processing/` (its own README, pipelines, manifests). All artifacts on S3 `s3://climate-ai-data-science-shiny-app-data/streamflow/`. Full detail: [EO_data_processing/README.md](EO_data_processing/README.md).

- **Watershed geometry layer** (prerequisite) — official basin polygons primary (US GAGES-II + Canada ECCC MDA_ADP), HydroBasins delineation fallback; 7,964 basins (8,018 universe − 54 >100,000 km²). `watershed_polygons_26jun2026.{gpkg,parquet}`. Codex-reviewed.
- **LAI (MCD15A3H, monthly, 2002–2024)** — 270-month panel, 2,150,280 rows (7,964 × 270). Two-stage spatial-distribution-of-monthly-mean (day-weighted per-pixel monthly mean → coverage-weighted basin mean + spatial heterogeneity stats), QA fractions from FparLai_QC/FparExtra_QC. Catalog (GEOLATELY/us-east-2) + LP DAAC backfill for far-N + 2024-11; 17 urban basins permanently NA (MODIS fill code 250). Codex-GO.
- **LULC (MCD12Q1, annual, 2001–2024)** — 191,136 rows (7,964 × 24), 109 cols: per-class % coverage for all 8 classification bands (LC_Type1–5 + LC_Prop1–3, 102 class columns) via coverage-weighted exactextract over VRT-mosaicked sinusoidal tiles. All 8 bands sum to 100 within 1e-13; manifest matches v061 spec. **Codex-reviewed 29 Jun (gpt-5.5) = GO**; 2 latent code hazards (tile-success accounting, unknown-code drop) hardened with the delivered output byte-identical.
- **Cross-cutting**: leading-zero canonical key reconciliation (metadata zero-stripped vs signatures zero-padded); LP DAAC HDF4 read via `pyhdf` (pip-GDAL lacks the driver); Earthdata account-lock handled by completeness-guard + retry/backoff (self-healing, resumable).
- **LAI `good_coverage_frac`** (30 Jun) — per basin-month continuous QA flag = valid pixel-obs / (all basin pixels × all expected composite dates); generalizes the binary `partial_month`. Added to `watershed_modis_lai_monthly_30jun2026.parquet` (data otherwise unchanged, hash-verified).
- **LAI QA/QC explorer** (30 Jun) — self-contained Leaflet HTML (`EO_data_processing/viz/build_lai_explorer.py` → `watershed_modis_lai_explorer.html`) for manual review: 7,964 points colored by any of 17 QA/summary variables, click a gage for its full monthly LAI time series (quantile bands + per-month good-coverage strip), map/TS hovers, glossary, processing notes. Mirrors `streamflow_explorer.html`.
- **LULC QA/QC explorer** (30 Jun) — self-contained Leaflet HTML (`EO_data_processing/viz/build_lulc_explorer.py` → `watershed_modis_lulc_explorer.html`, ~49 MB, on S3). 7,964 points colored by any of **27 map variables** in two non-overlapping groups: 10 IGBP-derived summary features (`class_diversity`, `delta_forest/urban/cropland` [robust mean(2001–03) vs mean(2022–24) endpoints], `dominant_class`, `dominant_group` [lump-robust aggregate dominance — fixes IGBP's class-splitting bias], `dominant_change`, `n_modis_pixels`, geom source / low-confidence) + all 17 individual IGBP classes (static `pct_*` roll-ups dropped as redundant with the per-class maps). Click a gage → stacked-area land-cover composition 2001–2024 with a **band switcher across all 8 MODIS schemes** (IGBP/UMD/LAI-fPAR/BGC/PFT/FAO-LCCS×3), official MODIS IGBP palette, per-year class-breakdown hover (0.1 % resolution). Glossary carries the MODIS interpretation caveats (SE-US woody-savanna behavior; 2021+ reduced confidence). Headless-validated (puppeteer/Chrome). **EO pipeline + viz COMPLETE.**

## [May 2026]

### New: per-gage watershed HydroATLAS metadata (non-signature product)
A standalone metadata file pairing each gage with the hydro-geophysical character of its **entire upstream watershed**, from HydroATLAS / BasinATLAS v10 (281 source attributes → ~211 output columns: climate, hydrology, terrain, land cover, soils/geology, anthropogenic). Sits alongside the signatures, keyed by leading-zero-safe `gage_id`: `data_out/watershed_hydroatlas_metadata_{date}.{csv,parquet}` + a generated data dictionary. 8,014 gages, ~13 MB.

- **Hybrid aggregation** over each gage's delineated upstream basin set (`upstream_hydrobasins.rds`): HydroATLAS upstream-accumulated `_u` / pour-point `_p` values pass through from the outlet basin (91 attrs — the delineated set equals HydroATLAS's own upstream extent, so `_u` is the rigorous watershed value); `_s`-only attributes are SUB_AREA-weighted (92 area-weighted means incl. monthly climate, with the HydroATLAS `-9999` NoData sentinel masked; elevation as spatial min/max; 11 categoricals as area-weighted majority — glc/pnv via argmax of upstream `_u` fractions, wet + the other 8 via area-weighted mode of the source class). Percentage / areal-extent fields are always area-weighted mean (never max). A `watershed_area_rel_diff` column ( |summed member SUB_AREA − gage `basin_area`| / `basin_area`, gage-reported drainage as truth) flags gages where the delineated watershed area diverges from the reported drainage (mis-snapped outlets).
- New `R/aggregate_hydroatlas_metadata.R` + `run_hydroatlas_metadata.R`; `HYDROATLAS_GPKG` added to `config.R`. R-only (metadata/ingestion) — not a cross-language signature, so no Julia/Python port.
- Validated: member SUB_AREA sums to outlet UP_AREA (median rel.diff 1e-4); area-weighted `_s` reproduces HydroATLAS `_u` (elevation cor 0.9994); 100% join to the signatures file (5,707/5,707). See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Watershed HydroATLAS Metadata.
- **Independent code review** (in-session) confirmed the classifier, NA-aware weighting, BFS delineation, and leading-zero joins; correctness fixes applied — `wet_cl_smj` switched from `_u`-fraction argmax to area-weighted mode (its fraction basis omits classes 10/11), `-9999` NoData masked in weighted means (`sgr_dk_sav`), and a malformed legacy-cache member id sanitized.

### New: static HTML data explorer
A self-contained, double-clickable `data_out/streamflow_explorer.html` (Leaflet + Canvas) for exploring the gages: 8,014 points colored by any of **412 variables** (dataset switcher between HydroATLAS watershed metadata and streamflow signatures); click a gage to draw its **watershed boundary** (border-only). Generators: `build_explorer.R` (unions the RAW lev12 basins per outlet for a clean dissolve — no interior sliver lines — then strips holes, simplifies, and caps at ~250 vertices → `watershed_borders.geojson` ~9.5 MB), `assemble_explorer.R` (joins points + injects into `explorer_template.html`). Full-resolution borders are ~111 MB (too heavy to inline); the clean-dissolve build is ~28 MB total. Border click-through and clean rendering verified headlessly (puppeteer-core + Chrome). See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Static HTML Explorer.

---

## [April 2026]

A large release: new signature families, the canonical-language transition, centralized NA handling, and cross-language alignment to 99.5% agreement. Output grew from 551 to **1,264 columns** (656 base/metadata + 608 Pettitt changepoint) across 7,313 gages. **Full per-change detail is in [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md) → [April 2026].**

### New signatures
- **Pettitt changepoint detection** applied to all 76 time-series signatures (8 fields each → 608 columns: cp year, p-value, pre/post means, delta, pct change, pre/post MK p-values). ~13.4% significance rate (2.7× null); signal-robustness analysis in [docs/SIGNATURES.md](docs/SIGNATURES.md) → Changepoint Detection.
- **Recession-parameterized baseflow** — `alpha_linear`, `BFI_Eckhardt_param`, `BFI_LyneHollick_param`, plus the per-gage `recession_alpha_point_cloud_linear_reservoir` scalar (Collischonn & Fan 2013). +25 columns.
- **Guidelines Section 3 additions** — D1/D99 flow timing, `n_recession_events`, `elasticity_annual`, `negative_ann`, `runoff_ratio_high_count`, 4 seasonal exclusion-year diagnostics, and 2 elasticity diagnostics.

### Canonical language transition
- **Julia replaced R as the canonical implementation** (~40× faster, modular). `rpkg/` is the production R port; monolithic `R/helperFunctions.R` is deprecated (ingestion utilities only). See [CLAUDE.md](CLAUDE.md).

### Centralized NA handling
- New `preprocess_daily_data()` runs once per gage before signatures: daily-grid normalization, ≤3-day interpolation, year rejection, seasonal completeness flags. Removed all per-signature min-days thresholds and `fillna(0)` calls; added the 80% trend-completeness gate to `generate_stats()`. Negative-Q rejection is config-driven (off by default); constant-SD is a flag only. See [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → NA Handling Architecture.

### Cross-language alignment
- Synced Section 3 plus the recession-algorithm rewrite to Python and rpkg; fixed Julia's Mann-Kendall tau-b, a `leftjoin` row-order bug, and several R-canonical divergences. The three synced implementations (Julia, Python, rpkg) agree at 99.5% R² (615/623 columns Perfect). Benchmark tables: [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) → Cross-Language Benchmarks.

### Tooling & audits
- Julia sensitivity-experiment framework (WY≥1993 ± qualifying-fraction filters), a Julia-vs-golden comparison pipeline, and interactive HTML validation dashboards (`docs/benchmarks/`).
- Three audits: Guidelines-vs-implementation, NA-handling (Codex), and a Julia adversarial review — all findings since addressed or documented.

---

## [March 2026]

### rpkg: Proper R Package
New R package (`rpkg/`) rewriting the monolithic `R/helperFunctions.R` as a proper installable package. 16 R source files, 6 test files, 27 man pages. R CMD check: 0 errors, 0 warnings. See [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md) for full details.

### Cross-Language Alignment (Rounds 1–6, complete)
Six rounds of alignment work brought R, Python, and Julia implementations to production-ready agreement. Final status (March 16-17, 2026):

| Pair | Mean R² | Median R² | Cols < 0.99 |
|------|---------|-----------|-------------|
| R vs Python | 0.9968 | 1.0000 | 17 |
| R vs Julia | 0.9965 | 1.0000 | 19 |
| Python vs Julia | 0.9997 | 1.0000 | 3 |

- 475 perfect (R²>=0.999), 56 good (0.99-0.999), 20 poor (<0.99)
- 531 of 551 columns (96.4%) have R² >= 0.99 across all 3 pairs
- Golden output regression: 551/551 perfect match against Feb 2026 golden outputs
- 20 remaining poor columns are all trend statistics (slopes, p-values) — underlying means/medians are near-identical
- Full round-by-round details: [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md), [docs/CROSS_LANGUAGE_STATUS.md](docs/CROSS_LANGUAGE_STATUS.md)

### Code Quality Improvements (Rounds 1–2)
Cross-language code review covering ~9,000+ lines. Key fixes: O(n²) metadata lookup → O(1), Python vectorized water year calc (~50x speedup), Julia DataFrame pre-allocation, Eckhardt BFI forward-fill aligned across all 3 languages. Full details in `docs/CODE_REVIEW.md` and [docs/CHANGELOG_ARCHIVE.md](docs/CHANGELOG_ARCHIVE.md).

### Benchmark Timing (March 16-17, 2026)

| Implementation | Time | Gages | Rate |
|---------------|------|-------|------|
| Julia | 9.2 min | 7,369 | 13.4/s |
| Python | 78.9 min | 7,369 | 1.56/s |
| rpkg | 114 min | 7,369 | 1.08/s |
| R (canonical) | 874 min* | 5,707 | 0.11/s |

*R ran concurrently with Python/Julia — timing inflated by I/O contention.

---

## Version History Notes

This project uses date-based versioning (MONTH YEAR) rather than semantic versioning, reflecting its nature as a research tool with continuous development.

### Output File Naming Convention
Output files include date stamps: `streamflow_signatures_full_JAN2026.csv`

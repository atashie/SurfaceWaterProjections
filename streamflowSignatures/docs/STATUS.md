# Project status — loaded every session (keep under 60 lines)

One line per item, each with a pointer; edit the line when the state changes, delete it when the item closes
(rule `.claude/rules/changelog.md`). Full text: `CHANGELOG.md` (open items + current month), `changelog-old.md`.

## Standard products (60 % qualifying fraction, drought family included; the 2026 deliveries have 1,653 columns, the version-2 reruns 2,277)

| # | Window | Folder | Gages | Annual rows | Climate input |
|---|---|---|---|---|---|
| 1 | WY 1993–2025 | `processedOuts_drought_28jul2026` | 6,678 | 18,898,406 | original Daymet parquet (since truncated on the drive) |
| 2 | WY 1980–2025 | `processedOuts_1980_2025_11aug2026` | 6,250 | 24,366,487 | `daymet_1980_2023_rebuilt_10aug2026.parquet` |
| 1 v2 | WY 1993–2025 | `~/HISSS_data/climate-rerun-12oct2026/processedOuts_1993_2025_climate_12oct2026` (2,277 columns; built 2026-10-12, validated, NOT yet deposited) | 6,678 | 27,602,161 | `daymet_1980_2025_29sep2026.parquet` |
| 2 v2 | WY 1980–2025 | `~/HISSS_data/climate-rerun-12oct2026/processedOuts_1980_2025_climate_12oct2026` (2,277 columns; built 2026-10-13, validated, NOT yet deposited) | 6,250 | 34,681,331 | `daymet_1980_2025_29sep2026.parquet` |

Neither is a subset of the other; record-dependent signatures are never compared across them. HydroShare deposit
(5 resources, staged): collection `f702201faa5d46069a5ee83ffa4c9768`. Public code mirror https://github.com/CZ-Sync/HISSS: publish after every merge to main (user preference, reaffirmed 2026-10-05; skip only for a stated one-off reason), from a clean checkout of HEAD while other sessions have uncommitted edits.

## Pending user decisions
- **Climate-metrics update (user request 2026-10-12) — BUILT AND VERIFIED, awaiting the user's review and deposit decision.** Decisions taken (21 answers 2026-10-12 + the crossing tolerance, confirmed 2026-10-13; `docs/plans/2026-10-12-climate-metrics-plan.md` §1). DONE: Julia family (39 bases, 624 columns → 2,277), Python and rpkg ports (blind reviews: Codex 4/4/6 rounds + Sonnet ×3, all GO), both version-2 products (table above: gates, annual validation, output review PASS; #1 bitwise identical to the pre-final run; `flagged_for_high_na` regenerated 560 / 401), explorers + comparison dashboards in the run folders, rpkg reference re-baselined on product #1 v2 (climate family bit-identical; residuals named), staged HydroShare folders R1/R2/R3 in `~/Downloads/Signatures_v2_12oct2026/` (local only). NEXT (user): review the staged folders, decide the deposit; relay list `docs/plans/2026-10-12-coauthor-relay.md`. CHANGELOG [October 2026].
- Daymet reprocess (decisions 2026-09-29: every polygon from scratch, full-res, all 6 variables, + 53 basins > 100k km²) — BUILT 2026-09-30, reviewed 2026-10-01: `daymet_1980_2025_29sep2026.parquet` (8,017 sites, 1980–2025) reproduces the co-authors' series on the 5,965 basins with original data; new basins and 2024–25 are spot-checked only. No product uses it yet (a rerun adds usable climate for 1,165 #1 / 615 #2 gages). Decisions 2026-10-01: questionable polygons FLAGGED, signature replay DROPPED, memory protected (Julia in gage batches, `docs/benchmarks/run_batched_julia.py`). NEXT: the rerun is now planned as part of the climate-metrics update below (user, 2026-10-12); open: co-author confirmation of polygons/weights/fill. `EO_data_processing/README_DAYMET.md`, action plan §0.
- Revoke the 2026-09-29 Earthdata token at urs.earthdata.nasa.gov (valid to 2026-10-25; it persists in local Claude transcripts; the run no longer needs it) — the user is handling it (2026-10-01).
- HydroShare doc updates H1–H3 (category terminology, high_na caveat, HYDAT release note) deliberately NOT applied 2026-09-10 — `docs/plans/2026-09-10-manuscript-category-edits.md` §C.

## Live known issues (full text: CHANGELOG → Known Issues)
- `flagged_for_high_na` is wrong in BOTH 2026 deliveries (metadata-only denominator; TRUE for every Canadian gage). Code fixed 2026-09-04; those products NOT rewritten (user decision). The version-2 reruns carry the regenerated flag (560 / 401 TRUE); the HydroShare READMEs / dictionary row are updated with the version-2 deposit.
- Stale Daymet input (2026 deliveries only): 5 sites NaN on every day 1980–2023, so 4 FL gages in #1 (3 in #2) have every climate signature NaN and "with climate" counts overstate by 4/3 (usable: 5,513 / 5,635). Resolved in the version-2 reruns (every gage of #1 and 6,223 of 6,250 in #2 carry climate summaries; the 27 others have 19 climate years, below the 20-value floor).
- HTML explorers/dashboards built before 2026-10-13 show CARTO's "API key required" watermark instead of a basemap (keyless tiles withdrawn): the deposited 2026 R1/R2 explorers + validation dashboards (superseded by the v2 deposit) and R5's three EO explorers (rebuild = user decision). Builders fixed (Esri gray canvas, OSM fallback, layer switch); v2 explorers/dashboards rebuilt and re-staged — CHANGELOG [October 2026].
- `ice_affected_days_total` is structurally 0 for every gage (Julia; cause not yet pinned; rpkg deliberately matches).
- 37 Canadian gages carry raw m³/s (`area_normalized = FALSE`, no HYDAT drainage area; user decision: no backfill); their Q-to-PPT signatures are NA by design; `flagged_for_qann_range` catches only 27/37 — downstream must filter on `area_normalized`. These counts are the canonical run's; the products hold 32 (#1) / 28 (#2), 22 / 25 flagged — to be corrected in an upcoming update (user, 2026-10-05).
- Storage year gate: Julia `unique` treats −0.0 ≠ +0.0 (1 annual row in 18.9 M); agreed fix is Julia → `==`, low priority, named waiver at gate time.
- Canonical `daymet_1980_2023.parquet` is truncated on the exFAT drive — always use the `_rebuilt_10aug2026` file; a readable original may exist in the Drive backup (would make product #1 reproducible).
- Staged HydroShare R4/R5 tables store 44 (+9) USGS ids un-padded; files kept as delivered, join rule documented (strip leading zeros on both sides).
- Legacy path only: `calculate_negative_days` crashes on `missing` Q (one-line fix pending); 6 pre-existing failures in `R/tests/test_na_handling.R`.

## Deferred fixes (none invalidates a delivered run)
- Daymet review 2026-10-01, LOW items deferred: no CMR retry/offline path and no vanished-granule check (`daymet_stream.py`); the explorer's dark-mode greyscale contrast and its group colours — CHANGELOG [October 2026].
- rpkg: reads `STREAMFLOW_SIGNATURES_CONFIG` not `STREAMFLOW_CONFIG`; non-canonical fallbacks when `na_handling` is absent; SWE merged only inside the PPT branch of the runner; per-run identifier + MANIFEST wanted in all three runners.
- Standard-product provenance: require a clean tree (both products logged `git_working_tree_dirty = true`) and `STREAMFLOW_HASH_INPUTS=1`; `annual_values.save` silently defaults to false when the config section is absent; `check_additivity.jl` needs an explicit cross-machine mode.
- Long-standing backlog: BFImax backward filter (Collischonn & Fan 2013); recession R² < 0.8 fit-quality flag; avg_storage "omitted from major analyses" decision; NA-handling item (i) flag-vs-reject wording; synchrony metrics; Julia ingestion port.

## Doc sync state (`/sync-docs`)
- Guidelines doc: last synced 2026-10-12 — NEW §3.9 Climate Metrics (39) + 5 snowfall metrics + the `_annual/_xSeason` naming convention, none implemented (see Pending); doc-side items in `docs/reconciliation/guidelines_todos.md` (top: duplicate "1.1" numbering; Gini/seasonality equations pasted as images; GDD written without base or clamp; parameterized-BFI interval must read (0, 1); Hecht 2024 reference missing; `recession_alpha_point_cloud` name; avg_storage / ice_affected_days_total blocks absent).
- Manuscript: last synced 2026-10-12 — relay list in `docs/reconciliation/manuscript_log.md` (top: §2.1.1 now defines seasons as Jan–Mar / Apr–Jun / Jul–Sep / Oct–Dec, CONTRADICTING the code's DJF/MAM/JJA/SON; §2 says nine families, §2.2.1 eight, both 121 signatures; §2.1.3 still 1980–2023 + gdptools; the seasonal-80 % and 12-flag sentences were dropped though implemented; carried: gap-rule wording; in-window 20-year floor; 85,000 vs 100,000 km²; cross-references; citations).
- Backups: S3 access lost 2026-08-24; the project Google Drive folder `1DVuq4nC5j_Y01sBaDj9cwjbv7S8sndjj` is the off-site backup (inventory unverified).

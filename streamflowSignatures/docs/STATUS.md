# Project status — loaded every session (keep under 60 lines)

One line per item, each with a pointer; edit the line when the state changes, delete it when the item closes
(rule `.claude/rules/changelog.md`). Full text: `CHANGELOG.md` (open items + current month), `changelog-old.md`.

## Delivered standard products (both 1,653 columns, 60 % qualifying fraction, drought family included)

| # | Window | Folder | Gages | Annual rows | Climate input |
|---|---|---|---|---|---|
| 1 | WY 1993–2025 | `processedOuts_drought_28jul2026` | 6,678 | 18,898,406 | original Daymet parquet (since truncated on the drive) |
| 2 | WY 1980–2025 | `processedOuts_1980_2025_11aug2026` | 6,250 | 24,366,487 | `daymet_1980_2023_rebuilt_10aug2026.parquet` |

Neither is a subset of the other; record-dependent signatures are never compared across them. HydroShare deposit
(5 resources, staged): collection `f702201faa5d46069a5ee83ffa4c9768`. Public code mirror https://github.com/CZ-Sync/HISSS.

## Pending user decisions
- Daymet reprocess to calendar 2025 over all 7,964 polygons — options + action plan in `docs/plans/2026-09-29-daymet-reprocessing-*.md`; feasibility tests T0–T4, go/no-go ~2026-10-10.
- HydroShare doc updates H1–H3 (category terminology, high_na caveat, HYDAT release note) deliberately NOT applied 2026-09-10 — `docs/plans/2026-09-10-manuscript-category-edits.md` §C.

## Live known issues (full text: CHANGELOG → Known Issues)
- `flagged_for_high_na` is wrong in BOTH delivered products (metadata-only denominator; TRUE for every Canadian gage). Code fixed 2026-09-04; products NOT rewritten (user decision) — regenerate at the next rerun of any data (`docs/benchmarks/recompute_high_na_flag.py --write`: #1 1,224 → 791, #2 1,243 → 598), then update the HydroShare READMEs/dictionary row.
- rpkg constant-series Mann-Kendall fixed in code 2026-09-04; the rpkg reference CSV still needs a benchmark rerun.
- `ice_affected_days_total` is structurally 0 for every gage (Julia; cause not yet pinned; rpkg deliberately matches).
- 37 Canadian gages carry raw m³/s (`area_normalized = FALSE`, no HYDAT drainage area; user decision: no backfill); their Q-to-PPT signatures are NA by design; `flagged_for_qann_range` catches only 27/37 — downstream must filter on `area_normalized`.
- Storage year gate: Julia `unique` treats −0.0 ≠ +0.0 (1 annual row in 18.9 M); agreed fix is Julia → `==`, low priority, named waiver at gate time.
- Canonical `daymet_1980_2023.parquet` is truncated on the exFAT drive — always use the `_rebuilt_10aug2026` file; a readable original may exist in the Drive backup (would make product #1 reproducible).
- Staged HydroShare R4/R5 tables store 44 (+9) USGS ids un-padded; files kept as delivered, join rule documented (strip leading zeros on both sides).
- Legacy path only: `calculate_negative_days` crashes on `missing` Q (one-line fix pending); 6 pre-existing failures in `R/tests/test_na_handling.R`.

## Deferred fixes (none invalidates a delivered run)
- rpkg: reads `STREAMFLOW_SIGNATURES_CONFIG` not `STREAMFLOW_CONFIG`; non-canonical fallbacks when `na_handling` is absent; SWE merged only inside the PPT branch of the runner; per-run identifier + MANIFEST wanted in all three runners.
- Standard-product provenance: require a clean tree (both products logged `git_working_tree_dirty = true`) and `STREAMFLOW_HASH_INPUTS=1`; `annual_values.save` silently defaults to false when the config section is absent; `check_additivity.jl` needs an explicit cross-machine mode.
- Long-standing backlog: BFImax backward filter (Collischonn & Fan 2013); recession R² < 0.8 fit-quality flag; avg_storage "omitted from major analyses" decision; NA-handling item (i) flag-vs-reject wording; synchrony metrics; Julia ingestion port.

## Doc sync state (`/sync-docs`)
- Guidelines doc: last synced 2026-09-29 — open doc-side items in `docs/reconciliation/guidelines_todos.md` (top: parameterized-BFI interval must read (0, 1) not [0, 1]; Hecht 2024 reference missing; `recession_alpha_point_cloud` vs the shipped `_linear_reservoir` name; avg_storage and ice_affected_days_total blocks absent).
- Manuscript: last synced 2026-09-29 (major co-author revision) — relay list in `docs/reconciliation/manuscript_log.md` (gap-rule wording; 20-year floor is in-window; 85,000 vs 100,000 km²; 6,041 → 6,087; cross-references after the §2.1.1/2.1.2 swap; citations).
- Backups: S3 access lost 2026-08-24; the project Google Drive folder `1DVuq4nC5j_Y01sBaDj9cwjbv7S8sndjj` is the off-site backup (inventory unverified).

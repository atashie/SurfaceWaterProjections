# Plan: Production Standard Run #2 — WY 1980–2025 ("entire period of record"), 60% Qualifying Fraction → D:/processedOuts_1980_2025_22jul2026

**Date**: 2026-07-22. **Status**: executing. §4 is the execution log.

The second of the two standard timeframes (user, 2026-07-21/22): "entire period of
record", operationalized as WY 1980–2025. Companion to standard run #1
(WY 1993–2025, `docs/plans/production_run_1993_2025_60pct_plan.md`, Codex GO).

## 1. Configuration

| Item | Setting | Basis |
|---|---|---|
| Water-year window | **WY 1980–2025 inclusive** (`STREAMFLOW_START_WATER_YEAR=1980`, `STREAMFLOW_END_WATER_YEAR=2025`) | User decision. End-cap path production-validated by run #1's independent audit. |
| Qualifying fraction | 0.60; denominator = 1980 .. gage's last in-window year (max 46 possible years → 60% = 27.6, so the fraction binds harder than the 20-year floor for full-window gages — opposite of run #1) | Same inclusion rule as run #1. |
| Feature set | Committed `97e0be0` — identical code to run #1 (July features + 60% overall trend gate + snow record-anchored decade gate + dashboard/tool fixes) | No code changes between the two standard runs. |
| Outputs | `D:/processedOuts_1980_2025_22jul2026/streamflow_1980_2025_60pct_22jul2026_*` | One unique folder per experiment (CLAUDE.md Critical Constraint #5). |
| Wrapper | `docs/benchmarks/run_julia_benchmark_prod_1980_2025_60pct.jl` | ENV-before-`using` pattern. |
| Memory | 16 GB machine; the runner's committed memory patches (climate-frame trim, raw-frame release, cache eviction) | The July 14 uncapped WY1980 run thrashed at 20.8 GB commit BEFORE these patches; they are now in the committed runner. Watch the log if runtime exceeds ~1 h. |

## 2. Expected characteristics & references

- More gages than run #1 is NOT expected a priori (1980 start admits long-record
  gages already in #1; the binding 60% × 46-year denominator excludes gages with
  large pre-1993 gaps that #1 admitted). The July 14 uncapped WY1980+ run's gage
  count is unknown here (outputs lost to the flash-drive rollback).
- References: `julia_signatures.csv` (April full-record canonical, 7,313 gages —
  gage-superset for the audit; column delta vs new run = exactly the 224 snow
  columns since April full already has Pettitt). The same-window April
  experiment does NOT exist for 1980 — loose comparisons only.
- Snow/climate families still end at WY 2023 (Daymet); WY 1980 is partial in
  Daymet → SWE-invalid (documented boundary behavior).
- `season_excluded_years_*` diagnostics: no phantom-2026 artifact this time (the
  cap excludes WY 2026 rows), but vs the April FULL-record baseline they will
  legitimately differ for pre-1980 records (window truncation).

## 3. Validation gates (same protocol as run #1)

1. `validate_production_run.py --csv <prod> --annual <prod parquet> --timing <json>
   --start-year 1980 --end-year 2025 --reference-csv docs/benchmarks/julia_signatures.csv
   --floor 20` (column gate expected to "fail" attributably: +224 snow columns).
2. `audit_qualification.jl <prod> docs/benchmarks/julia_signatures.csv 1980
   <parquet> 20 0.60 2025`.
3. `validate_annual_values.py --annual <parquet> --summary <csv> --floor 20`.
4. Spot checks: window bounds; num_water_years ∈ [20, 46]; overlap vs run #1's
   6,678 (expect large shared core); snow-gate + trend-gate behavior consistent.
5. Codex results review before declaring standard output #2.

## 4. Execution log

- 2026-07-22: wrapper + plan authored (post-97e0be0); output dir created; run
  launched.
- Run complete in **25.1 min**: **6,250 gages** × 1,488 columns; annual parquet
  21,824,537 rows / 90 signatures / wy [1980, 2025]. Log clean (expected
  area-normalized skips only). Memory a non-issue (38.5 GB free RAM pre-launch —
  the "16 GB machine" note in the July plan is outdated).
- **Gate results — ALL PASS (or attributed)**:
  - `validate_production_run.py`: 6/7 PASS; gate 2 column-equality FAIL fully
    attributed — vs the April full canonical (1,264 cols) run-only = EXACTLY the
    224 snow columns, ref-only = ∅. Overlap: all 6,250 gages ⊆ the canonical
    7,313 (run-only = 0). Snow sanity: 5,635 snow gages (4,546 US + 1,089 CAN).
    Floor respected; year counts [20, 46]; window bounds [1980, 2025].
  - `audit_qualification.jl` (end-cap arg): **AUDIT PASS** — early-ending gages
    show gage-capped denominators (22–30); late starters correctly excluded by
    the 1980-anchored 46-year denominator (0.457–0.565 < 0.60). The clamped
    edge-sampling (97e0be0 fix) handled this run's thinner strata.
  - `validate_annual_values.py --floor 20`: **PASS** — 525,329 pairs, 0 value
    mismatches, 0 floor violations, 0 orphans, 0 duplicate keys.
- **Comparison vs April FULL canonical** (in the run folder, via the new
  `--output-dir` defaults): **1,189/1,227 Perfect, min R² = 1.0000 in every
  non-recession category INCLUDING Pettitt fields; all 38 divergent columns are
  the log_a family** (16 + 16 stats/Pettitt + 6 seasonality scalars) = the
  intentional b=1 alpha change. Why so clean vs "full record": the source
  parquet begins at 1980, so the canonical's uncapped window ≈ this run's
  window for shared gages — run #2 is effectively the canonical analysis,
  properly capped + fraction-filtered + July features.
- **Run #1 ↔ #2 relationship** (product property, by design of the
  window-start-anchored denominator): 5,771 shared; 479 only-in-#2 (start
  1980–1992, 20–32 years — records ended before #1's window); 907 only-in-#1
  (late starters at e.g. 26/46 = 0.565 < 0.60 under the 1980 anchor). Neither
  standard is a subset of the other. `area_normalized = FALSE`: 28 rows here
  (32 in #1, 37 full-record — window-dependent).
- **Explorer built**: `signature_explorer_1980_2025_60pct.html` (61.6 MB, 1,456
  mapped variables, 90 sidecars 1980–2025).
- **Codex results review (2026-07-22): GO — zero findings** (no CRITICAL/MAJOR/
  MINOR on the outputs; sandbox blocked re-running the Julia audit, replaced by
  independent arithmetic verification from the CSVs). Confirmed: all 38
  divergent columns are log_a (na_mismatch = 0 on all 38); compared-column
  accounting fully resolved (1,227 = 1,264 common − 37 metadata/QA/diagnostic
  exclusions; min valid pairs 12); only-in-#2 inclusion arithmetic verified
  (fracs 0.6087–1.0 over the 1980..end denominator) and ALL 907 only-in-#1
  exclusions verified (fracs 0.4348–0.5870 < 0.60 over the 46-year anchor);
  `season_excluded_years_*` −1 pattern re-verified (4,404 gages, never +1); log
  clean; timing/CSV/parquet exactly consistent; 0 duplicate keys; QA flags all
  boolean with 0 percentile/timing/seasonal-sum violations; snow sanity + ssm
  ∈ [−1, 1] exact; explorer integrity (6,250 gages, 90 sidecars, 1980–2025).
- **VERDICT: standard output #2 DECLARED** —
  `D:/processedOuts_1980_2025_22jul2026/streamflow_1980_2025_60pct_22jul2026_*`.
- **2026-07-23 — cross-product annual-data consistency check (user request)**:
  verified run #1's annual parquet is CONTAINED in run #2's for shared gages —
  of run #1's 16,890,066 rows, the 1,394,200 absent from run #2 belong exactly
  and exclusively to the 907 run-#1-only gages (0 missing rows across all 5,771
  shared gages). Overlap-year VALUES are bit-identical for within-year-computable
  signatures (sampled Qann, Q50, D50_day, BFI_Eckhardt, flashinessRB,
  snow_on_dowy, annual_runoff_ratio: 0 differences over ~1.2M keys);
  record-dependent signatures differ by design (`n_high_pulses_all` 50,348
  differing years — period-of-record thresholds; `elasticity_annual` 150,885 —
  record-mean normalization; `BFI_Eckhardt_param` 154,411 — whole-record alpha).
  Interpretation guidance added to the claude-skill ("The Two Standard Output
  Products").

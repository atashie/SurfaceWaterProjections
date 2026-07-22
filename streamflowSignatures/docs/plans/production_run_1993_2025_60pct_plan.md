# Plan: Production Standard Run #1 — WY 1993–2025, 60% Qualifying Fraction → D:/processedOuts_22jul2026

**Date**: 2026-07-22. **Status**: executing. §5 is the execution log.

The first of the two "standard output" timeframes (user, 2026-07-21/22): WY 1993–2025,
and later "entire period of record" = WY 1980–2025. Matches HISSS manuscript §2.2.2's
1993–2025 subset.

## 1. Configuration

| Item | Setting | Basis |
|---|---|---|
| Water-year window | **WY 1993–2025 inclusive** — `STREAMFLOW_START_WATER_YEAR=1993`, `STREAMFLOW_END_WATER_YEAR=2025` | User decision 2026-07-21 (standard output #1). **First production use of the end-cap ENV** — the July 14 run executed WITHOUT an end cap (its plan §10; the 1993–2022 window in that plan's title was superseded mid-execution). CORRECTION to the 2026-07-21 reconciliation-log claim that the July run used "an explicit WY 1993–2022 window" — it did not; it was WY ≥ 1993 uncapped. |
| Qualifying fraction | `STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION=0.60` | Same gage-inclusion rule as April/July. Denominator = window start .. gage's last in-window year (capped at 2025 → max 33). Since the Feb 2026 parquet carries partial WY 2026 rows, the July no-cap denominators could reach 34 while WY 2026 never qualifies — so the capped run's gage set should be a **superset** of the July set (6,579). |
| Minimum years | `CFG_MIN_NUM_YEARS = 20` (config, unchanged) | 33 possible years × 0.60 = 19.8 → the 20-year floor binds for full-window gages. |
| Feature set | Committed `20b8c52`: all four July features + **60% overall trend gate** (2026-07-21) + **snow record-anchored decade gate** (2026-07-22, Codex GO) | Both new gates are in `config/signatures_config.json` (`min_fraction: 0.60`, `decade_min_fraction: 0.80`, `snow.record_anchored_decade_gate: true`). 80% decades + 60% overall confirmed by user 2026-07-22. |
| Inputs | Feb 2026 parquets on `D:` (verified present pre-launch: streamflow 0.84 GB, Daymet 3.84 GB, metadata CSV) | Defaults in `run_julia_benchmark.jl`. |
| Outputs | `D:/processedOuts_22jul2026/streamflow_1993_2025_60pct_22jul2026_{signatures.csv, signatures_annual.parquet, timing.json}` | 559 GB free on D:. |
| Wrapper | `docs/benchmarks/run_julia_benchmark_prod_1993_2025_60pct.jl` | ENV-before-`using` pattern. |

## 2. Expected divergences vs available baselines

The July 14 products were lost to the flash-drive rollback (`processedOuts_14jul2026`
is gone — the incident that motivated the repo-resident validators). Surviving
references (repo, `docs/benchmarks/`):

- **`startIn1993_60pct_signatures.csv`** (April 24; same start + fraction, NO end
  cap, pre-July features) — primary reference. Expected: column delta = exactly the
  +224 snow columns; gage overlap ≈ superset (cap shrinks denominators); on shared
  gages, non-recession/non-climate means/medians should agree closely (log_a columns
  intentionally diverged with b=1 alpha; snow columns new; Q-to-PPT gated for the 37
  un-normalized gages; trend stats differ by design via the 60% overall gate and the
  snow record gate — the latter only affects snow columns absent from April anyway).
- **`julia_signatures.csv`** (April 28 full run, no window) — superset gage list for
  the audit's excluded-candidate sampling.

## 3. Validation gates (repo-resident tooling)

1. `validate_production_run.py --csv <prod> --annual <prod parquet> --timing <json>
   --start-year 1993 --end-year 2025 --reference-csv startIn1993_60pct_signatures.csv
   --floor 20` — window bounds, column/overlap comparison, snow sanity, annual-parquet
   keys/volume, year-count consistency, floor respect.
2. `audit_qualification.jl <prod> julia_signatures.csv 1993 <parquet> 20 0.60 2025` —
   independent window+fraction reimplementation on stratified edge gages. **Extended
   2026-07-22 with the end-cap 7th argument** (caps valid years AND the denominator,
   mirroring the runner) — critical gate: the end-cap path has never run in production.
3. `validate_annual_values.py <parquet> <csv> --floor 20` — parquet↔summary
   mean/median consistency with floor awareness.
4. Spot checks: gage-count vs July's 6,579 (expect ≥, small delta); trend-stat
   population increase on 60–80%-complete series (new overall gate); snow trend NaNs
   at record-gate-fired gages.
5. Codex results review before declaring the product standard.

## 4. Out of scope

- Standard output #2 ("entire period of record", WY 1980–2025) — subsequent run.
- Python/rpkg ports (queued), commits during the run (user must request).

## 5. Execution log

- 2026-07-22: pre-flight OK (inputs verified, 559 GB free, output dir created);
  audit end-cap extension + wrapper + this plan authored; run launched.
- Run complete in **27.6 min**: **6,678 gages** (July's 6,579 + 99), 1,488 columns
  (1,456 signatures + 20 metadata + 12 QA), annual parquet 16,890,066 rows /
  90 signatures / 82.7 MB. Log clean (expected area-normalized gate messages only).
- **Gate results — ALL PASS (or attributed)**:
  - `validate_production_run.py`: 6/7 PASS. Gate 2 (strict column equality) FAILED
    against the only surviving reference (April `startIn1993_60pct`, 656 cols) —
    delta verified EXACTLY attributable: run-only = 832 = 608 non-snow Pettitt +
    224 snow (112 stats + 112 snow-Pettitt); ref-only = ∅; unattributed = ∅.
    Overlap: shared 6,579, run-only 99 (superset as predicted), ref-only 0. Snow
    sanity: 5,485 snow gages (4,457 US + 1,028 CAN), ssm ∈ [−1, 1]. Floor
    respected. Window bounds [1993, 2025]; year counts [20, 33].
  - `audit_qualification.jl` (with the new end-cap arg): **AUDIT PASS**, 11 edge
    gages — capped denominator (33 = 2025−1993+1) independently confirmed; all
    inclusion/exclusion decisions and per-gage year columns match.
  - `validate_annual_values.py --floor 20`: **PASS** — 545,266 pairs, 0 value
    mismatches, 0 floor violations, 0 orphans, 0 duplicate keys.
  - Spot checks: 60% overall gate gained 45 / lost 0 trend gages (Qann, Q50) on
    shared gages; snow record gate visible (snow_on_dowy 876 and melt_com_dowy 645
    trend-suppressed with means intact vs 333 for exempt swe_max); Pettitt
    non-NaN at all 876 suppressed snow_on_dowy gages.
- **Codex results review (2026-07-22): GO** — no CRITICAL/MAJOR. Independently
  re-verified: all 99 run-only gages legitimately included (num_water_years=20,
  20/33 = 0.606, no floor violations); shared-gage dense means/medians EXACT vs
  April (Qann/TQmean max abs diff 0.0); the 45 gained trends all sit in the
  expected [0.636, 0.793] completeness band (≥0.60, <0.80, ≥20 values); timing
  JSON ↔ CSV ↔ parquet internally consistent; 0 duplicate keys; column order
  preserved; snow-gate counts reproduced (876/645/333, Pettitt non-null at all
  876); log clean. Two MINOR, both attributed:
  1. **32 (not 37) `area_normalized=FALSE` rows** — NOT stale docs: the 37 count
     is for the full-record canonical output (7,313 gages); under this 1993–2025
     @ 60% window only 32 of those gages qualify (April's identically-started
     reference also has 32). All 81 Q-to-PPT columns NA for all 32 ✓. (Only 14
     gate log-messages appear because the message fires only for gages that have
     climate data joined; the other 18 have no Daymet rows.)
  2. **2 shared gages (06762500, 08HD020) lost `BFI_Eckhardt_mean/median`** vs
     April — each has only 17 non-null annual BFI values, below the 20-value
     stats floor (which April predates). Floor working as designed.
- **VERDICT: standard output #1 DECLARED** —
  `D:/processedOuts_22jul2026/streamflow_1993_2025_60pct_22jul2026_*`.
- **QA/QC dashboards built (2026-07-22)**, all in the run folder:
  - `signature_explorer_1993_2025_60pct.html` (+ `_annual/` sidecars, 90 files) —
    per-gage map explorer. Rebuilt same day with the Statistic picker extended
    from 8 to **16 statistics** (the 8 Pettitt changepoint fields added) so every
    signature-statistic combination is reviewable: 90 bases × 16 + 16 scalars =
    1,456 mapped variables.
  - `prod_1993_2025_60pct_vsApril60pct_vs_julia_{dashboard.html, comparison.csv,
    summary.md}` — PRIMARY QA vs the same-window April experiment: 597/619
    Perfect (min R² = 1.0000 in every non-recession category); all 22 divergent
    columns **within the comparator's signature scope** are the log_a family =
    the intentional b=1 alpha change (0 NA mismatches; b/concavity/alpha_linear
    pins Perfect). SCOPE QUALIFICATION (Codex finding, 2026-07-22): the DASHBOARD
    additionally visualizes 4 targets the comparator excludes as metadata — the
    `season_excluded_years_*` per-gage diagnostics — which diverge at 5,038
    gages. Fully attributed: the difference is EXACTLY −1 at every differing
    gage, for all four seasons, concentrated on gages with data into partial
    WY 2026 — the uncapped April run counted each phantom 2026 incomplete season
    as an excluded year; the capped run correctly does not. The new run's counts
    are the correct ones for a 1993–2025 product.
  - `prod_1993_2025_60pct_vs_julia_{dashboard.html, comparison.csv, summary.md}` —
    window-sensitivity view vs the April full-record canonical (broad divergence
    EXPECTED; concentrated in Pettitt pct_change/cp_year + window-dependent
    diagnostics).
- **NEW CONVENTION (user, 2026-07-22)**: ALL outputs of an experiment live in the
  experiment's own folder (this run: `D:/processedOuts_22jul2026`) — recorded in
  CLAUDE.md Critical Constraints; the six comparison artifacts were moved out of
  `docs/benchmarks/` accordingly.
- **Codex dashboard/results review (2026-07-22): initial NO-GO → all findings
  resolved.** Confirmed: explorer 16-stat extension correct (90 bases, 0 phantom
  pettitt bases, 1,456 cols verified in the built HTML; suffix-driven JS paths
  sound), artifact move complete in both directions, CLAUDE.md renumbering
  coherent. Findings + dispositions: (MAJOR) dashboard-vs-comparator scope
  mismatch on the 4 `season_excluded_years_*` diagnostics — divergence verified
  benign (uniform −1 WY2026 window artifact, above), scope now documented in BOTH
  tools' source + qualified here and in CHANGELOG; (MINOR) comparison tools
  defaulted outputs to `docs/benchmarks/` — both scripts gained `--output-dir`
  defaulting to the experiment CSV's folder, making the one-folder convention
  self-enforcing; (MINOR) `audit_qualification.jl` edge-sample slices could
  BoundsError on small runs — clamped with `first(v, n)`; (MINOR) untracked
  `temp/` — added to `.gitignore`. Explorer HTML size (64.6 MB, above the prior
  ~12–36 MB envelope) flagged as a watch item: consider splitting the map payload
  by category/stat if it grows further (sidecar approach for annual series
  already in place).

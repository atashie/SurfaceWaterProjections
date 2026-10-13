# Climate metrics (guidelines §3.9) — design record and run plan

**Date**: 2026-10-12 · **Status**: decisions taken (user, 2026-10-12); implementation in
progress. Julia first, then the reviews, the reruns, the ports.
**Trigger**: the guidelines doc gained §3.9 "Climate Metrics" (synced 2026-10-12,
`docs/reconciliation/guidelines_todos.md`); the reprocessed Daymet input
(`daymet_1980_2025_29sep2026.parquet`, 8,017 basins, calendar 1980–2025) is ready and no
product uses it yet (`EO_data_processing/README_DAYMET.md`).

## 1. Decisions (user, 2026-10-12)

| # | Question | Decision |
|---|---|---|
| 1 | Rename shipped columns to the doc's `_annual` / `_xSeason` convention? | NO. Shipped names stay; the convention applies to new metrics only; relay a doc-side revert; a rename manifest only if the co-authors insist. |
| 2 | Season suffixes for new metrics | the Q seasonal totals' suffixes as in the code: `win`, `spr`, `sum`, `fal` (not the doc's `sumr`, `fall`). |
| 3 | Prefix normalisation | yes: `precip_intermit_*` → `prcp_intermit_*`; `GDD_cum_*` → `gdd_cum_*`; relay. |
| 4, 6, 11 | Snowfall metrics (§3.8 insert) | REMOVED from the doc by the user. Not implemented. |
| 5 | Growing degree days | base temperature 0 °C, clamped: Σ max(0, (Tmax+Tmin)/2 − 0); the user updated the doc. |
| 7 | Dry day | prcp ≤ 0.1 mm (not the ETCCDI 1 mm). |
| 8 | `prcp_annual_sd` / `_cv` | within-year spread of the daily values → annual signatures with the 16 statistics. |
| 9 | 5-day maximum, ties, dry spells | window inside the water year; ties → first occurrence; the first longest dry spell. |
| 10 | Frost | basin-mean tmin < 0 °C; first frost searched from Oct 1; last frost within Oct 1 – Jun 30; NA when no day qualifies. |
| 12 | Gini, seasonality index | Gini over all days of the water year, zeros included; Walsh–Lawler index on the 12 calendar months of the water year. |
| 13 | Year qualification | the snow-family precedent: climate years are evaluated only for streamflow-valid years (`preprocess_daily_data`), precipitation metrics on PPT-valid years, temperature metrics on temperature-valid years (new policy, same thresholds, no negative-value check). |
| 14 | Gating of NaN-prone timing metrics | the snow timing precedent: the record-anchored decade gate for `frost_first_dowy` and `frost_last_dowy`, anchored to the temperature-valid years. |
| 15 | Rerun | both windows (WY 1993–2025 @ 60 %, WY 1980–2025 @ 60 %), same config otherwise, as new product versions superseding #1 and #2; the delivered products are not rewritten. |
| 16 | "Slightly altered" acceptance | per-column R² ≥ 0.999 on shared gage-years, every maximum reported and explained. |
| 17 | Ports | Julia products first; Python and rpkg validated against them before publication. |
| 18 | Fold-ins | regenerate `flagged_for_high_na`; rerun the rpkg reference benchmark (constant-series Mann-Kendall fix). |
| 19 | Reviews | Codex + Sonnet 5.5 per language, blind (guidelines text + code, no implementation narrative), plus one numeric cross-language audit. |
| 20 | HydroShare | stage locally only; no upload before the user's review. |
| 21 | Relay list | yes, one line per item, incl. the manuscript season error. |
| 22 | Crossing tolerance (review round 2, implementer decision 2026-10-12, confirmed by the user 2026-10-13) | `prcp_50_dowy` crosses when `cumulative ≥ fraction × total − 1e-12 × total`; 1e-12 of an annual total is physically nothing and the closest production margin is 1.3e-9; boundary pinned by tests in all three languages. |

Open with the co-authors (relay): the published guidelines page still carries the snowfall
lines and the GDD formula without a base (the user's edits had not propagated at the
2026-10-12 re-sync); the Gini and seasonality equations are images.

## 2. The 39 metrics — definitions as implemented

Category **Climate** (ninth; the manuscript's §2 preamble already lists it). Two function
families in `julia/src/climate.jl`: `calculate_precipitation_metrics` (20 bases, PPT-valid
years) and `calculate_temperature_metrics` (19 bases, temperature-valid years). Every base
gets the 16 statistics; no per-gage scalar. Columns: 39 × 16 = **624** → the summary CSV
grows from 1,653 to **2,277**; annual signatures 100 → **139**.

Conventions shared by both families: a year is used only when its daily grid is complete
(365/366 rows, dowy 1..n, no NaN after preprocessing); seasons are the calendar-month sets
of `na_handling.seasonal_completeness.season_definitions` (DJF, MAM, JJA, SON) applied
WITHIN the water year, so `fal` = Oct + Nov of the water-year start plus Sep of its end,
exactly as `Qfal` and `fall_runoff_ratio`; dowy 1 = Oct 1; ties take the first occurrence;
no Q-derived seasonal mask is applied (these are pure climate series).

### Precipitation (daily P, mm; PPT-valid years)

| Base | Definition | Units | Year NaN when |
|---|---|---|---|
| `prcp_annual` | Σ P | mm | — |
| `prcp_win`, `prcp_spr`, `prcp_sum`, `prcp_fal` | Σ P over the season's months | mm | — |
| `prcp_1d_max` | max P | mm | — |
| `prcp_5d_max` | max over i of Σ P[i..i+4], window inside the water year | mm | — |
| `prcp_max_dowy` | dowy of the first max P | day | — |
| `prcp_50_dowy` | first dowy with cumulative P ≥ 0.5 Σ P | day | Σ P ≤ 0 |
| `prcp_cm_dowy` | Σ (dowy · P) / Σ P | day | Σ P ≤ 0 |
| `prcp_ndays_abv_10mm`, `_25mm`, `_50mm` | count(P > 10 / 25 / 50 mm) | days | — |
| `prcp_annual_sd` | sample sd (n − 1) of the daily P | mm | — |
| `prcp_annual_cv` | sd / mean of the daily P | – | mean ≤ 0 |
| `prcp_intermit_annual` | longest run of dry days (P ≤ 0.1 mm); 0 when every day is wet | days | — |
| `prcp_intermit_annual_dowy_start` / `_dowy_end` | first / last dowy of the first longest dry run | day | no dry day |
| `prcp_gini` | G = Σᵢ (2i − n − 1) y₍ᵢ₎ / (n Σ y), y ascending, all n days incl. zeros | – | Σ P ≤ 0 |
| `prcp_seasonality_index` | SI = (1/R) Σₘ₌₁¹² \|Xₘ − R/12\|, Xₘ monthly totals, R = Σ Xₘ (0 … 1.83) | – | R ≤ 0 |

### Temperature (daily tmin, tmax, °C; temperature-valid years)

| Base | Definition | Units | Year NaN when |
|---|---|---|---|
| `temp_max_annual`, `temp_max_annual_dowy` | max of daily tmax; dowy of its first occurrence | °C, day | — |
| `temp_min_annual`, `temp_min_annual_dowy` | min of daily tmin; dowy of its first occurrence | °C, day | — |
| `temp_max_win/spr/sum/fal` | max of daily tmax over the season's months | °C | — |
| `temp_min_win/spr/sum/fal` | min of daily tmin over the season's months | °C | — |
| `frost_first_dowy` | first dowy with tmin < 0 °C | day | no frost day (gated, §4) |
| `frost_last_dowy` | last dowy with tmin < 0 °C among days from Oct 1 to Jun 30 | day | no such day (gated, §4) |
| `gdd_cum_annual` | Σ max(0, (tmax + tmin)/2 − 0) | °C·day | — |
| `gdd_cum_win/spr/sum/fal` | the same over the season's months | °C·day | — |

## 3. Data-flow changes (Julia canonical)

1. `read_parquet` normalises `tmin` → `TMIN`, `tmax` → `TMAX` beside `prcp` → `PPT`,
   `swe` → `SWE`.
2. `preprocess_daily_data` gains a temperature block (f3) mirroring the SWE block: a year is
   temperature-valid when BOTH tmin and tmax have ≤ `max_raw_na_per_year_temp` raw NaNs,
   every internal gap ≤ `max_interpolation_gap_temp` days is linearly filled, and no NaN
   remains; no negative-value check (temperatures are signed). Returns `valid_temp_years`.
   The block runs only for years that passed the streamflow checks (decision 13). Existing
   outputs are byte-identical when TMIN/TMAX are absent.
3. `calculate_all_signatures` gains `temp_data` (a frame filtered to temperature-valid years,
   explicit opt-in like `snow_data`; a TMIN/TMAX column inside `gage_data` is never used
   implicitly). Precipitation metrics run on the PPT-valid `climate_data` frame whenever it
   carries PPT, independently of `area_normalized` (no Q is involved). Both families sit in
   the per-family `try/catch` and are config-gated: an absent `climate` section emits no
   climate columns (the drought precedent), which is also the additivity-proof switch.
4. Config: `climate` section (`dry_day_threshold_mm` 0.1, `wet_day_thresholds_mm`
   [10, 25, 50], `multi_day_max_window_days` 5, `cumulative_fraction` 0.5, `frost_temp_c`
   0.0, `frost_last_end_month` 6, `gdd_base_c` 0.0, `record_anchored_decade_gate` true) and
   `na_handling.temperature_na_policy` (`max_raw_na_per_year_temp` 30,
   `max_interpolation_gap_temp` 3). Bundled copies byte-identical. Julia bakes at precompile
   time — purge the cache.
5. Runner (`run_julia_benchmark.jl`), smoke test, and `run_batched_julia.py` keep six
   climate columns (`site_id, Date, prcp, swe, tmin, tmax`); the runner builds `temp_data`
   from `valid_temp_years`. Python and rpkg runners likewise.
6. Registries (`/add-signature` §2): `EXPECTED_SIGNATURE_BASES` (`config.R`),
   `EXPECTED_DENSE_SIGNATURES` (Julia collector test), the `nunique() == 100` gate in
   `validate_production_run.py` (→ 139), `docs/signature_categories.csv` (+39 rows),
   Python `validate_schema()`, rpkg schema tests, the documented column total.

## 4. Gates

Standard: trend completeness 60 % / 80 % first-and-last decade, the 20-value stats floor,
Pettitt fields. Record-anchored decade gate (snow precedent, `snow.jl`): anchor years =
temperature-valid, grid-complete years (`temp_max_annual` non-NaN); `frost_first_dowy` and
`frost_last_dowy` need ≥ `decade_min_fraction` of the anchor years in BOTH the first and the
last decade of the record to be non-NaN, else their six trend statistics are NaN
(mean/median + Pettitt still computed); record span < 10 years → gate skipped.

## 5. Verification and review plan

1. Unit tests (`julia/test/test_climate_metrics.jl`): hand-built years with closed-form
   answers (single-storm year → Gini (n−1)/n, SI 1.833…, 50 % and centre-of-mass dates on the
   storm day; dry-spell runs with the ≤ 0.1 mm rule and ties; thresholds exactly at 10 mm not
   counted; 5-day window at the year end; leap-year grid; the non-contiguous fall season;
   GDD clamp with sub-zero days; frost dates incl. a frost-free year and the Jun 30 cut-off;
   the record decade gate; key-set contract with and without the family; byte-identical
   existing outputs when TMIN/TMAX are absent).
2. Smoke test on real gages, then a one-batch benchmark; `check_signature_failures.py` on
   the log; `check_additivity.jl` NEW (climate on) vs TWIN (climate section removed), same
   input and gages: every shared column bitwise identical, exactly 624 added, populated.
3. Blind adversarial reviews (decision 19): each reviewer receives §3.9 of the guidelines,
   the decision table above, and the code of one language; it must derive expected values
   independently and report discrepancies with reproductions. Findings are resolved and
   re-reviewed before the reruns.
4. Ports: Python and rpkg module-for-module; cross-language fixtures ≤ 4e-14;
   `check_schema_equality.py`, `check_annual_parquet_equality.py`, rpkg acceptance gates.

### Verification record (2026-10-12)

| Step | Result |
|---|---|
| Julia unit suite | green; `test_climate_metrics.jl` 355 assertions (closed-form years + review regressions) |
| Smoke test, 10 real gages, reprocessed Daymet subset | PASS; all 39 bases present, ranges plausible, seasonal means add up |
| Additivity, batch 1 (824 gages), climate on vs config section removed | PASS twice (before and after the review fixes): 1,652 shared columns bitwise, +624, all populated ≥ 0.79, 0 family failures |
| Same-machine control (2026-10-01 run of product #1's config vs product #1) | Q-only columns bitwise identical; 98 climate-dependent columns ≤ 3.4e-13 |
| Blind review, Julia, Codex round 1 | NO-GO → 7 findings fixed (CHANGELOG [October 2026]) |
| Blind review, Python, Codex round 1 | NO-GO → 7 findings fixed (exact gap filler, signed zeros, strict config, shared labels, finite guards) |
| Blind review, rpkg, Codex round 1 | NO-GO → 10 findings fixed, incl. the library-wide Spearman (exact = FALSE) and double-precision Mann-Kendall rewrites |
| Three-way cross-check after round 1 | 1,170 annual values bit-identical; all 624 summary values within 4e-14 / 1e-12 in BOTH ports (R Mann-Kendall residual gone) |
| Unit suites after round 1 | Julia green (355 climate assertions); Python 174; rpkg 1,179 (full suite) |
| Codex round 2 (all three languages) | NO-GO on a residual list (crossing tie, non-finite policy, config strictness, Python observation contract, R overrides / seasons) → all fixed; regression tests added; suites green (Julia 366 climate assertions, Python 174, rpkg 1,190); cross-check unchanged |
| Sonnet blind review, Julia | GO; 2.1 M real annual values recomputed independently, 0 mismatches; documentation caveat (Dec 31 of leap calendar years) corrected |
| Codex round 3 | Julia GO (one MINOR: canonical percentage range, fixed); ports NO-GO on custom-input boundaries only (override validation, Arrow-decimal dtype, R season validation / preprocessing default, endpoint canonicalisation) → all fixed; suites green (Julia 366 climate assertions, Python 175, rpkg 1,190); cross-check unchanged |
| Codex round 4 (verification) | Julia GO; Python GO; rpkg NO-GO on two MAJORs (`as.numeric` coercion of strings / complex observations; SWE policy overrides ignored and untyped) → fixed (`.real_values()` contract, SWE policy from the config via the typed readers); rpkg climate tests 107; cross-check unchanged |
| Sonnet blind review, Python | GO; 593 differential case-runs vs Julia (410,156 annual rows, 368,560 summary values, 0 mismatches), real data bit-identical; 4 MINORs fixed (write-back only when valid, unconditional residual rejection, NumPy integer limits, exports); crossing-tolerance boundary pinned in all three suites (Python 177, Julia 369 climate assertions) |
| Codex round 5 (rpkg verification) | round-4 MAJORs and round-3 regressions PASS; one new MAJOR (a list-typed input column re-bound as a list, breaking the downstream families) → fixed (one up-front `.real_values()` conversion of the observation columns), regression test added, rpkg suite 1,212, cross-check unchanged; Codex round 6: GO (no new findings) |
| Sonnet blind review, rpkg (resumed on the round-5 build) | GO; 38,047 synthetic annual cells 0 bit differences, real-gage runner parity, exact-rational reference; 3 MINORs fixed (roxygen glue, JSON nulls, stale fixture header), NOTEs recorded |
| rpkg reference re-baseline (decision 18) | DONE 2026-10-13 on the round-5 build: 8 batches (166 min of R), 6,678 gages, 2,277 columns; schema gate PASS, failure scans PASS, annual gate PASS at 1e-6 with named residuals (1 signed-zero `avg_storage` row + 319 discrete tie flips of 26.8 M pairs); the 39 climate signatures bit-identical in every annual cell |
| Output-review tool on the PRE-FINAL product #1 (round-2 code) vs the delivered product | PASS on every gate: 1,252 Q-only + metadata columns identical per gage; annual values on shared gage-years R² ≥ 0.999 and p99 relative ≤ 1e-3 for all 24 continuous climate-dependent signatures (runoff ratios ≤ 1.8e-4 absolute); discrete snow timings changed in ≤ 0.05 % of gage-years (threshold crossings); the old `_mean`s reproduced from the new annual values on the old year sets (R² 1.000, elasticity_annual 0.9993); 624 new columns populated (min finite fraction 0.785) and within physical ranges; seasonal means add up. Summary-level trend statistics move for every gage because two climate years were added — expected and reported, not gated. |
| FINAL production rerun, product #1 (round-3 Julia) | DONE 2026-10-12: 0 family failures, gates 5/5, annual validation PASS, output review PASS; bitwise identical to the pre-final run (2,277 × 6,678; 27,602,161 annual rows) — `~/HISSS_data/climate-rerun-12oct2026/processedOuts_1993_2025_climate_12oct2026` |
| FINAL production rerun, product #2 | DONE 2026-10-13 (batch 3 had been stopped by the 120-min cap under swap pressure; resumed from the two finished batches; merge needed `--memory-limit 3GB`): 0 family failures, gates 5/5, annual validation PASS, output review PASS (1,252 Q-only + metadata columns identical on 6,250 gages; 24 continuous signatures R² ≥ 0.999 on 223,715 shared gage-years; discrete ≤ 0.05 %; new columns populated ≥ 0.766, in range) — `processedOuts_1980_2025_climate_12oct2026`, 34,681,331 annual rows |
| Explorer basemap fix (user review 2026-10-13) | CARTO tiles now watermark keyless requests ("API KEY REQUIRED"); the twelve Leaflet builders moved to `hisssBasemap()` (Esri gray canvas, OpenStreetMap fallback, layers control); headless check PASS (`check_html_basemap.py`); both v2 explorers + comparison dashboards rebuilt and the explorers re-staged |

## 6. Reruns (decision 15) and output reviews

Inputs: streamflow `combined_streamflow_data_09feb2026.parquet`, metadata
`combined_watershed_metadata_09feb2026.csv`, climate `daymet_1980_2025_29sep2026.parquet`
(verify size 5,134,193,404 and `PAR1`), GAGES-II dir. Batched (`run_batched_julia.py`, 8
batches), clean tree, `STREAMFLOW_HASH_INPUTS=1`, one folder per run:
`processedOuts_1993_2025_climate_<date>` and `processedOuts_1980_2025_climate_<date>`.

Expected versus the delivered products, by column class:
- **Q-only columns**: identical values — same streamflow, same qualification (gage sets
  6,678 / 6,250 unchanged). The same-machine question is settled by comparing the
  2026-10-01 control run (product #1's configuration on this laptop) with product #1.
- **Climate-dependent existing columns** (runoff ratios, elasticity, Q-P seasonality,
  storage, snow): three separable effects — (i) input values on the 5,965 basins the stale
  file covered (prcp totals within ±0.0011 %); (ii) two added climate years, WY 2024–2025,
  for every gage; (iii) 1,165 (#1) / 615 (#2) gages gaining a climate series. Review
  (`docs/benchmarks/compare_climate_rerun.py`, gates per decision 16 on shared GAGE-YEARS):
  (2a) annual values on shared gage-years — continuous signatures R² ≥ 0.999 and 99th-
  percentile relative difference ≤ 1e-3; discrete ones (snow days / dowy / ssm) reported
  as the fraction of gage-years changed; elasticity reported only (unstable by
  construction); (2b) the old product's `_mean` recomputed from the NEW annual values on
  the OLD year set, R² ≥ 0.999 per signature. The summary-level statistics are expected to
  move for every gage (effect ii) and are reported per column, not gated.
- **New columns**: sniff tests — physical ranges (prcp_annual vs the HydroATLAS climate
  attributes; temperature extremes vs latitude/elevation; frost dates vs `snow_on_dowy`;
  Gini and SI in range; dry spells longer in arid basins), internal identities
  (`fall_runoff_ratio × prcp_fal ≈ Qfal`, `swe_max_to_ppt × prcp_annual ≈ swe_max`,
  `prcp_50_dowy ≤ prcp_cm_dowy` ordering tendencies), spatial coherence maps in the explorer.
- `flagged_for_high_na` regenerated by the run itself (decision 18); the rpkg reference
  benchmark rerun on the new reference.

## 7. Explorers, deposit, docs

Explorer per product (`build_signature_explorer.py` + `_annual/` sidecar) and the
product-vs-product comparison dashboards into the run folders; staged HydroShare folders
(R1/R2 new versions, R3 climate input + basin flags, dictionary + category tables) in
`~/Downloads/Signatures`, not uploaded. Docs: `docs/SIGNATURES.md` §14, README, user skill,
DEVELOPMENT.md, `signature_categories.csv`, CHANGELOG, STATUS, run-benchmark skill, CLAUDE.md
column counts; relay list for the co-authors.

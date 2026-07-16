# Plan: Production Rerun — WY 1993–2022, 60% Qualifying Fraction → D:/processedOuts_14jul2026

**Status**: EXECUTED + VALIDATED (2026-07-15). Codex plan review NO-GO → 6 findings →
GO (§9); Codex RESULTS review: **GO, "well validated"** (3 LOW findings, addressed).
> **PROVENANCE NOTE — the executed window differs from this plan's title/§1/§5**: after
> two user course-corrections during execution, the run used the EXACT April
> `startIn1993_60pct` config — **WY ≥ 1993, NO end cap, 60% fraction** — not the
> 1993–2022 window planned below. §10 is the authoritative record of the executed
> configuration; the end-cap harness feature described in §2 was built, reviewed, and
> retained but UNUSED in this run. Validation packet:
> `D:/processedOuts_14jul2026/validation/RESULTS.md` (9/9 gates PASS).
**Requested by**: user, 2026-07-14: "rerun our analysis to capture our updated signatures
and the annualized signature parquet... for the period of 1993 to present with a
requirement of 60% completion (not 80%). Save all outputs to D:/processedOuts_14jul2026."
Corrected mid-planning: "match our previous 60% run exactly (apart from the additional
metrics, and the change to recession analysis alpha calculations) which I believe ran
from 1993 to 2022 wy, not 1993 to present wy."

## 1. Interpretation lock

| Item | Setting | Basis |
|---|---|---|
| Water-year window | **WY 1993–2022 inclusive** | User correction. Historical note for the record: the April 2026 JULIA experiment `startIn1993_60pct` ran "WY1993 to gage max" (no end cap — verified from its summary report); the 1993–2022 window matches the older R-era run (`streamflow_signatures_full_1993-to-2022-min-20yrs.csv`). The new run uses the explicit 1993–2022 window as directed. |
| "60% completion" | `min_qualifying_data_fraction = 0.60` — the gage-INCLUSION rule from the April experiment: ≥60% of the gage's possible water years within the window must pass year qualification | Matches the `startIn1993_60pct` experiment the user references; "(not 80%)" contrasts with the `_80pct` sibling. The per-metric 80% TREND-completeness gate (`na_handling.trend_completeness`) is a different knob and stays at 0.80, matching the previous 60% run. |
| Minimum years | `CFG_MIN_NUM_YEARS = 20` unchanged | Same as the previous run. Note the interplay: 30 possible years (1993–2022) × 60% = 18, so the 20-year floor (66.7%) is the binding constraint for full-window gages; the 60% rule additionally binds gages whose records end before 2022 (per-gage denominator = gage's max in-window year − 1993 + 1). Identical semantics to April, window-capped. |
| Feature set | The full uncommitted working tree: annual-values export (parquet), b=1 recession alpha, snow metrics (14), area-normalized gate | All four are user-approved July work; tests green (1,108 assertions), all Codex-GO. |
| Explicitly NOT included | Runoff-ratio seasonal-flag bug fix (CHANGELOG Known Issues) | User has not opted in; outputs will carry the documented bug, same as every prior run. |

**"Match exactly" divergence contract** (vs a baseline run of the PRE-July code on the
identical window/filters): the user named two expected differences — (a) the additional
metrics (14 snow bases → +224 columns; plus the new annual-values parquet artifact),
(b) the recession alpha change (log_a_pointcloud, log_a_events — 16 fields each — and
the 6 log_a_seasonality_* scalars). **Two further side effects of already-approved July
work must be surfaced because they also change outputs**:
- `flagged_for_high_na` — the 224 snow columns enter its NA-fraction denominator; some
  no-SWE (Canadian) gages flip. Accepted policy (snow plan §6, Pettitt precedent);
  flips will be counted and reported.
- The 37 `area_normalized = FALSE` Canadian gages now have the ENTIRE Q-to-PPT family
  structurally NA instead of unit-invalid garbage (user decision, July 2026). The
  whitelist must enumerate the full family, not just 8-stat time series (Codex
  finding 1): all runoff-ratio, elasticity, qp-seasonality, and avg_storage stat +
  Pettitt columns, PLUS the scalars `elasticity_static`, `elasticity_years_total`,
  `elasticity_years_low_ppt`, `runoff_ratio_high_count`, PLUS the QA flags these NAs
  feed: `flagged_for_elasticity_range`, `flagged_for_runoff_ratio_range` (both can
  flip when the underlying value becomes NA), and `flagged_for_high_na`. Divergent
  ROWS for this family must be exactly the 37-gage set.

Everything else shared between baseline and production must agree at R² = 1 /
identical values, and **both runs must select the IDENTICAL gage set** (no July feature
touches year qualification — SWE handling only adds `valid_swe_years`).

**Comparison-tooling caveat** (Codex finding 1): `compare_experiment_vs_julia.py`
EXCLUDES metadata and QA-flag columns from its shared-column comparison — it covers the
signature columns only. QA flags and metadata therefore get their own dedicated exact
check (§5 gate 4b); relying on the compare script alone would silently skip them.

## 2. Code changes required (small, additive, benchmark-harness only)

1. **`docs/benchmarks/run_julia_benchmark.jl`**:
   - `const OUTPUT_DIR = get(ENV, "STREAMFLOW_OUTPUT_DIR", @__DIR__)` (currently
     hardcoded `@__DIR__` — outputs must reach D:).
   - **End-water-year support** (does not exist today): read
     `STREAMFLOW_END_WATER_YEAR` at runtime alongside the start-year read; extend the
     existing post-preprocess filter block (data/valid_years/valid_climate_years/
     valid_swe_years/seasonal_flags/diagnostics masks gain `.<= end`), and cap the
     qualifying-fraction window (`wy_in_range = start <= yr <= end`; per-gage
     `total_possible = max(in-range years) − start + 1` as today, now range-capped).
     Print the window in the config banner.
   - **Precondition (Codex finding 4)**: the start/fraction filters — and therefore the
     new end cap — exist ONLY in the non-legacy preprocess branch; the legacy branch has
     no window logic at all. This rerun is valid because `use_legacy_filtering=false` in
     BOTH the working tree and b7014c2 (verified) — asserted at pre-flight (§4.0) and
     stated here explicitly: the end-cap patch is intentionally non-legacy-only.
2. **New production wrapper** `docs/benchmarks/run_julia_benchmark_prod_1993_2022_60pct.jl`
   (pattern of the existing experiment wrappers — ENV before `using`):
   ```julia
   ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
   ENV["STREAMFLOW_END_WATER_YEAR"] = "2022"
   ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
   ENV["STREAMFLOW_OUTPUT_PREFIX"] = "streamflow_1993_2022_60pct_14jul2026"
   ENV["STREAMFLOW_OUTPUT_DIR"] = raw"D:/processedOuts_14jul2026"
   include(joinpath(@__DIR__, "run_julia_benchmark.jl"))
   ```
3. No package (`julia/src/`) changes — the harness edits don't touch signature math, so
   the 1,108-assertion suite is unaffected (it will be re-run once anyway as a
   pre-flight formality).

## 3. Baseline strategy (the "match exactly" proof)

**Problem**: the April Julia outputs (full run + all experiments) are no longer on disk
(gitignored `docs/benchmarks/*.csv`, cleaned) — only the summary .md reports survive.
And no historical run used the 1993–2022@60% window anyway. So the baseline is
**generated fresh from the pre-July code**:

1. `git worktree add <scratchpad>/baseline_head b7014c2` — HEAD is docs-only above the
   April-era code; all July feature code is uncommitted (verified: 43 modified/untracked
   paths in the working tree, b7014c2 touches one docs file). The worktree is an
   isolated checkout; the working tree is never touched.
2. Copy the CURRENT `julia/Manifest.toml` into the worktree (gitignored, so the worktree
   has none) → identical package versions for both runs; `Pkg.instantiate()` to verify.
   **FAIL CLOSED** (Codex finding 3): if the copied Manifest does not instantiate
   cleanly against HEAD's Project.toml, ABORT the baseline (investigate; do not fall
   back to a fresh resolve — silent version drift would invalidate the "identical
   dependencies" premise). Record `Pkg.status` output for both runs in validation/.
   (Julia precompile caches are keyed per-project/environment, so the worktree and main
   tree do not share compiled StreamflowSignatures state; each process also reads its
   OWN checkout's `config/signatures_config.json` via `config.jl`'s `@__DIR__`-relative
   default — no config leakage.)
3. Apply the SAME end-water-year patch (§2.1) to the worktree's copy of
   `run_julia_benchmark.jl` (HEAD lacks it too). Harness-only, ~10 lines, identical
   logic — recorded as a known asymmetry (the baseline runs patched-HEAD, not pure
   HEAD), and independently audited by §5 gate 9 so a bug in the shared patch cannot
   hide behind matching outputs (Codex finding 2).
4. **Dedicated baseline wrapper FILE inside the worktree** (Codex finding 6) —
   `run_julia_benchmark_baseline_1993_2022_60pct.jl`, same ENV-before-`using` pattern
   (start 1993, end 2022, frac 0.60, prefix `baseline_head_1993_2022_60pct`). Both runs
   launch as separate `julia <wrapper>.jl` processes; NO shell-session ENV mutation, so
   nothing can leak between runs. Baseline outputs stay inside the worktree's
   docs/benchmarks (HEAD has no OUTPUT_DIR override and doesn't need one).
5. Baseline output: ~1,264 columns (656 base/metadata + 608 Pettitt), no snow, no
   annual parquet (that feature is uncommitted — expected).

**Comparison**: `python docs/benchmarks/compare_experiment_vs_julia.py
prod14jul_vs_head --baseline <worktree CSV> --experiment <D: CSV>` (the script accepts
explicit paths — verified). Verdict gates in §5.

**Serialization**: the two runs load ~8–12 GB of parquet each on a 16 GB machine —
baseline and production runs execute SEQUENTIALLY, never concurrently (I/O + memory
contention corrupted timings in past concurrent runs).

## 4. Execution sequence

0. **Pre-flight**: `D:/processedOuts_14jul2026` exists and is EMPTY (verified); D: has
   556 GB free (verified); source parquets present (verified paths in benchmark
   constants); no other Julia processes running; **assert `use_legacy_filtering=false`
   in BOTH the working-tree config and the worktree's checkout config** (Codex
   finding 4 precondition).
1. Apply §2 edits; re-run the Julia test suite (formality — harness edits only).
2. Build the worktree baseline (§3.1–3.3).
3. **Run 1 — baseline** (background, monitored; expect ~8–15 min): verify completion
   banner, gage count N_base, column count, output CSV row count = N_base.
4. **Run 2 — production** (background, monitored; expect ~10–18 min, snow adds a few
   minutes of compute + the annual parquet write): outputs land in
   `D:/processedOuts_14jul2026/`.
5. **Validation** (§5).
6. **Cleanup**: copy the baseline CSV + comparison report into
   `D:/processedOuts_14jul2026/validation/`; `git worktree remove` the baseline
   checkout; revert NOTHING in the main tree (the §2 edits are permanent additive
   features of the benchmark).
7. Report to user: gage count, column count, runtimes, every validation gate's verdict,
   flagged_for_high_na flip count, snow coverage summary, deliverables listing.

## 5. Validation gates (all must pass before reporting success)

| # | Gate | Expectation |
|---|---|---|
| 1 | Window honored | Production CSV: min(start_water_year) ≥ 1993 AND max(end_water_year) ≤ 2022; annual parquet water_year ∈ [1993, 2022] |
| 2 | Identical gage set | Baseline and production CSVs: same gage_id set, EXACT (∆ = 0 gages) |
| 3 | Column arithmetic | Production CSV = baseline CSV columns + the 224 snow columns, EXACTLY (the `area_normalized` metadata column already exists at b7014c2 — the July gate is behavioral, adding no CSV columns; the annual parquet is a new ARTIFACT, not a CSV delta). Name-level set diff must equal the snow-column list |
| 4 | Shared signature-column agreement | Via the compare script: all shared signature columns R² = 1 / identical EXCEPT the whitelist: `log_a_pointcloud*`, `log_a_events*` (16 fields each), `log_a_seasonality_*` (6 scalars), and the FULL Q-to-PPT family (runoff ratios, elasticity incl. `elasticity_static`/`elasticity_years_total`/`elasticity_years_low_ppt`/`runoff_ratio_high_count`, qp seasonality, avg_storage) — whose divergent rows must be exactly the 37 area_normalized=FALSE gages |
| 4b | QA-flag + metadata exact equality (dedicated check — the compare script EXCLUDES these columns) | Row-by-row equality of all `flagged_*` and metadata columns between baseline and production, with differences permitted ONLY in `flagged_for_high_na` (no-SWE gages), `flagged_for_elasticity_range` and `flagged_for_runoff_ratio_range` (restricted to the 37 gated gages); every differing row enumerated and attributed |
| 5 | Annual parquet ↔ summary consistency | `validate_annual_values.py --annual <D:>... --summary <D:>...` exits 0 (means/medians reproduce; no duplicate keys; coverage rules) |
| 6 | Snow sanity | swe_max_mean ≥ 0 everywhere; ssm_mean ∈ [−1, 1]; DOWY means ∈ [1, 366]; snow columns non-NA only for Daymet-covered (US) gages; snow coverage count reported |
| 7 | flagged_for_high_na accounting | Flip count vs baseline reported; flips restricted to gages with structurally-NA families (no-SWE / no-climate) |
| 8 | Annual parquet volume | Rows ≈ N_gages × ~90 signatures × ≤30 yr (order 15–20 M); zstd size ~100–200 MB; row count matches timing JSON `n_annual_rows` |
| 9 | **Independent qualification audit** (Codex finding 2 — the end-cap patch is shared by both runs, so matching outputs alone cannot prove the window math) | A standalone audit script INDEPENDENTLY reimplements the ~6 lines of window/fraction math (from `preprocess_daily_data` outputs, NOT the benchmark's code path) for a stratified sample of edge gages — records starting before 1993, ending before 2022, and sitting near the 20-year and 60% thresholds, drawn from both included and excluded sets — and emits an audit table (`all_wy`, capped `wy_in_range`, `total_possible`, `n_valid`, fraction, decision) that must match the run's actual include/exclude decisions. Plus two global assertions on the production CSV: every included gage has `num_water_years ≥ 20` AND recomputed fraction ≥ 0.60 |

## 6. Deliverables (all in `D:/processedOuts_14jul2026/`)

| File | Content |
|---|---|
| `streamflow_1993_2022_60pct_14jul2026_signatures.csv` | Summary signatures + metadata + QA flags (~N × ~1,490 columns) |
| `streamflow_1993_2022_60pct_14jul2026_signatures_annual.parquet` | Long-format per-year values (gage_id, signature, water_year, value) |
| `streamflow_1993_2022_60pct_14jul2026_timing.json` | Run timing + n_annual_rows |
| `validation/` | Baseline CSV (patched-HEAD, same window), comparison CSV + summary md, validate_annual_values log, gate-by-gate verdict notes |

## 7. Roadblocks identified & mitigations

| Roadblock | Mitigation |
|---|---|
| Benchmark has NO end-water-year mechanism (previous "1993" runs were start-only) | §2.1 adds it, mirrored into the baseline worktree; window gate §5.1 proves it worked in both runs |
| OUTPUT_DIR hardcoded to docs/benchmarks | §2.1 ENV override |
| April baseline outputs deleted from disk | Fresh baseline from a HEAD worktree (§3) — better anyway, since no historical run used this exact window |
| July features are uncommitted → any git operation touching the working tree risks them | Worktree checkout is fully isolated; NO stash/checkout/reset in the main tree at any point |
| Worktree package environment could resolve different versions | Copy the current Manifest.toml in; verify instantiate; record versions |
| 16 GB RAM vs ~8–12 GB per run | Strictly sequential runs; no concurrent Julia/R processes |
| Daymet ends calendar 2023 | Irrelevant under the 2022 cap — every window year has full climate + SWE coverage (cleaner than a to-present run; no partial-tail climate years) |
| Changepoint config window 1980–2024 | No-op vs the 1993–2022 series (window is a filter, series already inside it); cp columns compute with 30-yr series (min_total_obs=20 satisfied for dense metrics) |
| Gage-count expectation | April's 6,579 was 1993→gage-max (uncapped) — NOT directly comparable; the anchor is baseline == production EXACT (§5.2). Rough expectation: ≤ 6,579 |
| `elasticity_rolling` short series (30 yr → ≤20 rolling values) | Same as the April experiment; known Pettitt inflation caveat already documented |
| Known runoff-ratio seasonal-flag bug persists in outputs | Documented (CHANGELOG Known Issues); explicitly out of scope per user's silence — surfaced again at report time |
| validate/compare scripts assume docs/benchmarks paths | Both have explicit path overrides (`--annual/--summary`; `--baseline/--experiment`) — verified |
| Windows paths in ENV | Wrapper uses raw string; Julia handles forward slashes on Windows |
| Timing artifacts of two long runs | Both runs backgrounded with completion notification; no polling loops against the Julia processes |

## 8. Out of scope

- Promoting these outputs to `golden-outputs/` (user call, later).
- Any commits (user must request).
- Python/rpkg ports (three-feature queue, post-validation).
- Trend-completeness gate change (stays 0.80 — flagged in §1 in case the user meant
  this knob by "60% completion"; one-line change if so).
- Runoff-ratio seasonal-flag fix.

## 9. Codex review record (2026-07-14)

Initial verdict **NO-GO-as-is** (job task-mrl7a152-2q1kvz). All 6 findings incorporated:

| # | Severity | Finding | Disposition |
|---|---|---|---|
| 1 | HIGH | Divergence whitelist incomplete (Q-to-PPT scalars `elasticity_static`/`elasticity_years_*`/`runoff_ratio_high_count`; QA flags `flagged_for_elasticity_range`/`flagged_for_runoff_ratio_range`) AND the compare script EXCLUDES QA/metadata columns entirely | Whitelist fully enumerated (§1); new dedicated gate 4b: row-by-row QA-flag + metadata equality with attributed exceptions |
| 2 | HIGH | Both runs share the new end-cap patch — a bug in the window/fraction math would pass gates 2–4 undetected | New gate 9: standalone audit script independently reimplements the window math from `preprocess_daily_data` outputs for stratified edge gages (pre-1993 starts, pre-2022 ends, near-threshold) + global `num_water_years ≥ 20` / fraction ≥ 0.60 assertions on the production CSV |
| 3 | MEDIUM | "Fallback to fresh instantiate" would silently accept dependency drift in the baseline | FAIL CLOSED: abort if the copied Manifest doesn't instantiate; `Pkg.status` recorded for both runs (§3.2) |
| 4 | MEDIUM | End-cap patch is non-legacy-branch-only; plan didn't state the `use_legacy_filtering=false` precondition | Stated in §2.1; pre-flight assertion for BOTH checkouts (§4.0) |
| 5 | MEDIUM | Gate 3 wrongly expected an `area_normalized` column delta — the column already exists at b7014c2; the July gate is behavioral | Gate 3 corrected: CSV column delta = the 224 snow columns EXACTLY |
| 6 | LOW | Baseline launched via shell-ENV mutation risks leakage into the next run | Dedicated baseline wrapper FILE in the worktree; both runs as separate `julia <wrapper>.jl` processes (§3.4) |

Review also confirmed: per-checkout config isolation (each process reads its own
checkout's signatures_config.json via `@__DIR__`), per-project precompile caches (no
cross-tree contamination), and all claimed script path-overrides.

**Delta-verification** (job task-mrl7j50v-hf7nin): all 6 RESOLVED, "New problem: none."
**Final verdict: GO.** Execution proceeds per §4.

## 10. Execution log & deviations (2026-07-14/15)

1. **Feature commit landed mid-execution**: the user committed the four July features
   (`419ab84`), simplifying the tree — only the rerun harness remained uncommitted.
2. **Data-drive incident**: `D:` is a flash drive; the device holding
   `processedOuts_feb2026` (and the April outputs) was swapped/unplugged mid-session —
   the baseline worktree and output dir "vanished" (rebuilt), and the first baseline
   launch died on missing inputs. The user restored the three input files.
3. **OOM + memory patch**: the second baseline attempt hit OutOfMemoryError in Phase 2
   (`add_water_year_columns` copying the 8-column ~98M-row climate frame with the
   111M-row streamflow frame resident; 6.7/15.7 GB free). Fix applied IDENTICALLY to
   both trees: `select!(climate, intersect(["gage_id","date","PPT","SWE"], names))`
   before the water-year copy — result-identical (dropped columns are never read by
   any signature path in this runner).
4. **USER DECISION — fresh baseline run cancelled** (paused mid-load, then dropped):
   the April outputs turned out to SURVIVE on the reconnected drive —
   `startIn1993_60pct_signatures.csv` (April 60% experiment: same start/fraction, no
   2022 cap, pre-July features) and `julia_signatures_wPettitt-trends.csv` (April full
   run) — and the user judged these, plus
   `streamflow_signatures_full_1993-to-2022-min-20yrs.csv`, sufficient as baselines.
   Consequence for §5: gates 2/3/4/4b lose their EXACT-match form (no identically-
   windowed baseline exists; the April experiment's window runs past 2022, so all
   values legitimately differ) and are replaced by LOOSE comparisons vs the April
   files (column-set arithmetic vs the April full run; high-correlation sanity on
   shared gages/columns; snow/log_a expected-divergent). Gates 1, 5, 6, 7 (semantics
   via unit test), 8, and 9 (independent window audit) apply unchanged and carry the
   rigor. The features themselves remain validated by the 1,108-assertion suite + two
   Codex GO cycles.
5. **FINAL WINDOW (user, after two course-corrections)**: the first production launch
   (1993–2022) was killed when the user asked to match the April experiment more
   closely; a proposed 2024 cap was then corrected on evidence (the April run's own
   output shows 73% of gages ending at WY 2025 — no cap ever existed). Final config =
   **EXACT April window: WY ≥ 1993, NO end cap, 60% fraction** (wrapper
   `run_julia_benchmark_prod_1993_60pct.jl`, prefix `streamflow_1993_60pct_14jul2026`).
   Consequence: `startIn1993_60pct_signatures.csv` (April) is once again a STRICT
   baseline — identical window + filters, pre-July code — so the §5 exact-match gates
   (2/3/4/4b) are restored in full against that file, with the §1 divergence whitelist
   unchanged. The end-cap harness feature remains in the benchmark (unused this run;
   covered by the gate-9 audit only for the start/fraction math actually exercised).
   Expected gage count: 6,579 (April parity). Production log:
   `D:/processedOuts_14jul2026/validation/production_run.log`.

# Plan: Export Per-Year Annual Signature Values (Long-Format Parquet)

**Status**: IMPLEMENTED (2026-07-14, Julia). Plan was Codex-reviewed (initial verdict
NO-GO-as-is; all required amendments incorporated, see §9) before implementation.
Benchmark re-run pending — bundled with one additional workflow feature.
**Scope decision (user)**: long-format parquet only (no wide table, no CSV); Julia only
(Python/rpkg ports deferred).
**Date**: 2026-07-14

---

## 1. Problem

The pipeline persists (a) raw daily streamflow (parquet, R ingestion) and (b) the final
per-gage signature summary (`julia_signatures.csv`, 1,264 columns = 8 stats + 8 Pettitt
fields per signature). The intermediate product — **annualized metric values per water
year per gage** (e.g., Qann by year, BFI_Eckhardt by year, D50_day by year) — is computed
inside every signature function as a local `annual_data` DataFrame, consumed by
`generate_stats()`, and discarded. ~76 time-series signatures × 7,313 gages × ~20–45
years of per-year values are lost on every run.

Saving them enables: per-year QA/QC, custom trend windows, changepoint re-analysis,
year-level cross-gage analysis, and time-series views in the explorer/Shiny app — all
without re-running the ~15–27 min extraction.

## 2. Key architectural fact (verified; scope statement per Codex review)

**Every annual series that feeds summary statistics flows through `generate_stats()`**
(`julia/src/stats.jl:249`) as a DataFrame with a `water_year` column plus one or more
value columns. This includes the two irregular series: `elasticity_rolling` (keyed to
the END year of each 11-yr window, `elasticity.jl:157-190`) and `elasticity_annual`
(keyed to the later year of each consecutive-year pair, `elasticity.jl:196-216`).
Signature base names are globally unique across all call sites (verified — they merge
into one flat result Dict), so a single long table has no collisions.

**Export semantics (precise)**: the collector records **the exact series handed to
`generate_stats()`** — i.e., exactly what the summary statistics were computed from.
This is the property the §4.2 validation exploits. It is NOT "every annual intermediate
ever computed": some callers prune rows before the call, and NaN-year representation
varies by signature:

- **Dense signatures** (frame pre-allocated with one row per water year; NaN
  placeholders survive into the export): flow volumes (`flow_volumes.jl:73`), FDC
  (`fdc.jl:115`), timing (`timing.jl:57`), pulses (`pulses.jl:220`), runoff ratios
  (`runoff_ratios.jl:59`), recession base metrics (`recession.jl:303,494-498`).
- **Caller-pruned / sparse signatures** (rows for non-computable years are absent, not
  NaN): `qp_bimodality` (NaN rows filtered into `bi_df` before the call,
  `qp_seasonality.jl:188-190`), `n_recession_events` (NaN years filtered into
  `events_df`, `recession.jl:465-468`), and signatures whose frames are built by
  appending only computed years (flashiness `flashiness.jl:78`, storage
  `storage.jl:153`, baseflow `baseflow.jl:192,307`, elasticity rolling/yoy
  `elasticity.jl:178-216`).

Consumers must therefore treat **absent row** and **NaN value** as equivalent
("not computable that year"); the schema note in §3.5 documents this. Extending the
pruned callers to emit dense NaN rows is possible later but is explicitly NOT part of
this change (it would alter what generate_stats receives — risk to the frozen contract).

What does NOT flow through `generate_stats()` — and is correctly out of scope:

- Per-gage scalars (documented 8-stat exceptions): `elasticity_static`, recession
  seasonality (6 values), `runoff_ratio_high_count`, `elasticity_years_*`,
  `ice_affected_days_total`, `recession_alpha_point_cloud_linear_reservoir`,
  `season_excluded_years_*`. These are not annual series.
- Early-return paths that call `empty_stats()` without calling `generate_stats()`
  (e.g., `flow_volumes.jl:151-160` when `nrow(annual_data) < 3`; elasticity when
  `< 3` valid values). These gages/signatures produce all-NaN summary stats today and
  will simply have **no rows** in the annual file. With the 20-year minimum this is
  rare in practice; the validation script treats "no rows + NaN summary" as consistent.

## 3. Design

### 3.1 Collector type (`julia/src/stats.jl`)

```julia
mutable struct AnnualCollector
    signature::Vector{String}
    water_year::Vector{Int32}
    value::Vector{Float64}
end
AnnualCollector() = AnnualCollector(String[], Int32[], Float64[])
```

A lightweight append-only accumulator. One instance per gage (created in the benchmark
loop), drained into global vectors tagged with `gage_id` after each gage.

### 3.2 Hook inside `generate_stats()`

Add keyword `collector::Union{Nothing, AnnualCollector} = nothing`. Behavior:

1. **Hoist the `value_cols` resolution** to the top of the function (it is currently
   duplicated in the `nrow < min_rows` early-return branch at `stats.jl:262-268` and the
   main branch at `stats.jl:289-297`; both compute "all Real columns except year_col").
   The hoist must **literally preserve** the current predicate (`eltype(...) <: Real`,
   exclude `year_col`) and must resolve against the **unsorted** input `df` — resolving
   after `sort(df, year_col)` or broadening the eltype predicate (e.g., to admit
   `Union{Missing, Float64}`) would change behavior for hidden callers. Note hidden
   callers exist: tests exercise explicit `year_col` kwargs
   (`julia/test/test_changepoint.jl:162,201,220,238,250`) — the collector must key on
   the passed `year_col`, never a hardcoded `water_year`.
2. **Collect before any gating inside generate_stats**: immediately after resolving
   `value_cols` — i.e., *before* the `nrow < min_rows` early return, *before* the
   per-column `length(valid_values) < min_rows` NaN-out, and *before* the 80%
   trend-completeness gate. NaN values that reach `generate_stats` are collected (e.g.,
   seasonal flow-volume metrics NaN'd by `seasonal_flags` — "year qualified, metric not
   computable"). Caller-side pruning upstream of `generate_stats` is out of the hook's
   reach by design (see §2 export semantics).
3. Collection is **read-only** with respect to the stats computation. With
   `collector=nothing` (the default), behavior is byte-identical to today — no change
   to the summary CSV contract, no change for Python/rpkg parity, no change for any
   existing caller (tests, qa_qc).

```julia
if collector !== nothing
    for col in value_cols
        col in names(df) || continue
        for row in eachrow(df)
            # Year validation FIRST — skip rows that cannot be keyed.
            # (Int32(round(NaN)) would throw InexactError; never convert blind.)
            y = row[year_col]
            (y === missing || (y isa AbstractFloat && !isfinite(y))) && continue
            v = row[col]
            push!(collector.signature, col)
            push!(collector.water_year, Int32(round(Float64(y))))
            push!(collector.value, v === missing ? NaN : Float64(v))
        end
    end
end
```

(Implementation will use vectorized appends rather than `eachrow` for speed; sketch
above is for semantics. The `missing` handling on values is defensive — current
callers never pass `missing` into value columns, since `Float64.(df_sorted[!, col])`
at `stats.jl:304` would already error, but the early-return branch never converts, so
cheap insurance.)

**3.2a Year handling**: `water_year` arrives as Int (most callers) or numeric-
convertible. Rows whose year is `missing` or non-finite are **explicitly skipped**
(cannot be keyed) rather than coerced — coercion would throw. No such rows exist
today; the skip rule is a guard, and the unit test covers it with a synthetic frame.

**Collection order**: rows are collected in the order the caller's DataFrame provides
(not `df_sorted`); determinism comes from a final global sort (§3.4).

### 3.3 Threading the kwarg (mirrors the `changepoint` plumbing exactly)

`calculate_all_signatures()` (`julia/src/signatures.jl`) gains
`collector::Union{Nothing, AnnualCollector}=nothing` and forwards it to every signature
call. Each of the 13 signature functions gains the same kwarg and forwards it to its
`generate_stats()` call(s):

| Module | Function | generate_stats calls |
|---|---|---|
| flow_volumes.jl | calculate_flow_vols_by_year | 1 |
| flashiness.jl | analyze_flashiness_trends | 1 |
| timing.jl | analyze_flow_timing_trends | 1 |
| fdc.jl | analyze_fdc_trends | 1 |
| baseflow.jl | analyze_baseflow_indices | 1 |
| baseflow.jl | analyze_baseflow_indices_with_parameters | 1 |
| recession.jl | analyze_recession_parameters | 1 (annual metrics; seasonality scalars bypass) |
| pulses.jl | calculate_pulse_metrics | 1 |
| pulses.jl (or module hosting it) | calculate_negative_days | 1 |
| runoff_ratios.jl | analyze_Q_PPT_relationships | 1 |
| elasticity.jl | calculate_streamflow_elasticity | 2 (rolling + yoy) |
| qp_seasonality.jl | calculate_qp_seasonality | 1 |
| storage.jl | calculate_average_storage | 1 |

This works identically for the legacy (`use_legacy_filtering: true`) and new
preprocessing paths, since the hook is downstream of both.

### 3.4 Benchmark integration (`docs/benchmarks/run_julia_benchmark.jl`)

In the Phase 4 loop:

```julia
collector = CFG_SAVE_ANNUAL_VALUES ? AnnualCollector() : nothing
signatures = calculate_all_signatures(...; collector=collector)
if collector !== nothing
    n = length(collector.value)
    append!(annual_gage_id, fill(gage_id_str, n))   # n refs to one String — cheap
    append!(annual_signature, collector.signature)
    append!(annual_water_year, collector.water_year)
    append!(annual_value, collector.value)
end
```

After Phase 4 (new sub-phase, timed in the timing JSON):

```julia
annual_df = DataFrame(gage_id=annual_gage_id, signature=annual_signature,
                      water_year=annual_water_year, value=annual_value)
sort!(annual_df, [:gage_id, :signature, :water_year])   # deterministic output
Parquet2.writefile(joinpath(OUTPUT_DIR, "$(OUTPUT_PREFIX)_signatures_annual.parquet"),
                   annual_df; compression_codec=:zstd)
```

- **Parquet2 write path is currently UNVERIFIED** (Codex finding): `Parquet2` is a
  declared dependency (`julia/Project.toml:16`) but has **no `[compat]` entry** and no
  checked-in `Manifest.toml`, and the codebase only uses it for reads
  (`io.jl:30`). Implementation step 0 is therefore: add a `[compat]` pin for Parquet2,
  and run a standalone **write→read round-trip smoke check** (4-col typed DataFrame,
  incl. an empty frame) against the resolved version to confirm the `writefile`
  signature and codec kwarg (`:zstd` preferred, `:snappy` fallback, uncompressed
  acceptable — file stays modest either way) BEFORE integrating into the benchmark.
- The annual parquet respects `OUTPUT_PREFIX`, so sensitivity experiments
  (`startIn1993` etc.) automatically produce their own annual files, consistent with
  their filtered inputs.
- Written **before** Phase 5 (metadata merge) so a metadata failure can't lose it.
- **Zero-gage guard** (`run_julia_benchmark.jl:279-284`): when zero gages qualify, the
  runner writes an empty summary CSV and returns early — the annual parquet is
  **skipped** in that branch (log a note). The validation script treats a missing
  annual file alongside an empty summary CSV as consistent.
- Legacy-path note: by Phase 4, legacy runs are already filtered to `qual_years`
  (`run_julia_benchmark.jl:218-235`) and non-legacy climate signatures receive
  `climate_data` filtered to `valid_climate_years` (`:237-258`), so collected series
  are automatically consistent with each path's inputs — no extra handling needed.

### 3.5 Schema

Long format, 4 columns:

| Column | Type | Notes |
|---|---|---|
| `gage_id` | String | Zero-padded original ID (same as summary CSV `gage_id`) |
| `signature` | String | Signature base name (e.g., `Qann`, `BFI_Eckhardt`, `D50_day`) — matches summary column prefixes |
| `water_year` | Int32 | Water year (Oct 1 – Sep 30). Semantics for irregular series: `elasticity_rolling` = END year of the 11-yr window; `elasticity_annual` = later year of the consecutive pair |
| `value` | Float64 | Annual metric value. **NaN and absent-row are equivalent** ("not computable that year") — dense signatures emit NaN placeholders, caller-pruned signatures omit the row entirely (see §2 export semantics) |

Scale estimate: 7,313 gages × ~76 signatures × ~20–45 years ≈ **20M rows**, roughly
100–250 MB as compressed parquet. Memory during accumulation is dominated by 3 × 20M
primitive/pointer vectors (~500 MB peak) — trivial next to the 111M-row input parquet
already held in memory.

### 3.6 Config switch

`config/signatures_config.json` gains:

```json
"annual_values": {
    "save": true,
    "comment": "Write per-year annual signature values (long parquet) alongside the summary CSV. Julia only; Python/rpkg ports deferred."
}
```

`julia/src/config.jl` reads it into `CFG_SAVE_ANNUAL_VALUES` — **default `false` when
the section is absent** (Codex finding: config is loaded once at module import into
constants, `config.jl:15-23`, with no runtime override pattern; a `true` fallback would
silently enable a large new artifact for any old/external config). The repo's
`signatures_config.json` ships with the section present and `save: true`, so the
canonical pipeline gets the feature explicitly.

`CFG_SAVE_ANNUAL_VALUES` must also be added to the **export list** in
`julia/src/StreamflowSignatures.jl:84-91` — the benchmark uses bare `CFG_*` names
(`run_julia_benchmark.jl:43-55`), which only work for exported constants.

## 4. Validation

### 4.1 Unit + regression tests (`julia/test/`)

Synthetic gage (~25 years, known values):
1. Run `calculate_all_signatures` twice — with and without a collector — and assert the
   returned Dicts are **identical** (collector has zero effect on stats).
2. Assert collected `Qann` rows equal independently computed annual sums; assert row
   counts per signature; assert NaN rows appear for a deliberately incomplete season
   (flow volumes — a dense signature; runoff ratios cannot be used for this until the
   seasonal-flag name bug in §9/finding 2 is fixed).
3. Assert the collector's year-skip guard: a synthetic frame with a missing/NaN year
   row collects the remaining rows and does not throw.
4. Assert no duplicate (signature, water_year) pairs per gage.

**Anti-silent-failure guard** (Codex finding: `calculate_all_signatures` wraps every
signature in `try/catch` that only warns, `signatures.jl:52-157` — a botched kwarg
thread would silently drop a signature's summary AND annual rows):
5. Run the synthetic gage under `Test.@test_logs` (or a logger capture) and assert
   **zero warnings** are emitted with the collector enabled.
6. **Whole-pipeline signature-coverage regression**: assert the with-collector run
   produces the exact same set of result-Dict keys as the without-collector run, and
   that the collector contains rows for the full expected set of signature base names
   (a committed list, ~76 names) for a fully-populated synthetic gage with climate.

### 4.2 Cross-file consistency script (`docs/benchmarks/validate_annual_values.py`)

After the next benchmark run: read the annual parquet + `julia_signatures.csv`; for
every (gage, signature), recompute mean/median over non-NaN annual values and compare
to the `_mean`/`_median` summary columns (rtol 1e-9). Rules:

- Summary mean NaN + no rows in parquet → consistent. This arises from THREE paths:
  (a) early `empty_stats()` returns that never call `generate_stats`, (b) caller-side
  pruning to an empty/too-small frame (e.g., `bi_df`/`events_df` with <3 rows), and
  (c) the whole-frame `nrow < min_rows` early return — rows ARE collected in case (c),
  so specifically: NaN summary + `count(non-NaN rows) < 3` is consistent.
- Summary mean NaN + rows present with `count(non-NaN) < 3` → consistent
  (per-column `min_rows` guard at `stats.jl:313` NaNs the stats but values were computed).
- Otherwise mean/median must match. Note the summary computes mean/median even when the
  trend-completeness gate zeroes out the 6 trend stats (`stats.jl:366-381`), so
  gated signatures are still fully checkable.
- **Coverage check** (complements the anti-silent-failure guard in §4.1): for every
  gage, the set of signatures with a non-NaN `_mean` in the summary CSV must equal the
  set of signatures having ≥1 non-NaN row in the parquet. Any asymmetry flags a broken
  collector thread or a silently-swallowed signature exception.
- **Full-output regression**: the next benchmark run's summary CSV must be
  value-identical to a pre-change reference run (the collector must not perturb any
  statistic) — compare against the current golden Julia output.

This is a strong end-to-end proof that the file contains exactly the series the
statistics were computed from.

## 5. Documentation updates

- `CHANGELOG.md` → `[Unreleased]`: feature entry + "Port annual-values collector to
  Python/rpkg" added under Planned (deferred per user).
- `docs/DEVELOPMENT.md`: data-flow diagram + file-structure tree + a short "Annual
  Values Export" subsection (schema table, year-key semantics, validation script).
- `docs/SIGNATURES.md`: one paragraph under Overview pointing at the annual parquet and
  the NaN semantics.
- `claude-skill/streamflow-signatures.md`: mention the new artifact (skill-maintenance
  rule in CLAUDE.md).

## 6. Explicitly out of scope / deferred

- **No wide table, no CSV** (user decision — long parquet only).
- **Python/rpkg ports deferred** (user decision). The kwarg-with-nothing-default design
  means cross-language numerical parity of the summary outputs is unaffected meanwhile.
- **No change to the summary CSV** (frozen contract) and no change to any computed
  statistic.
- **No benchmark re-run yet** — implementation will be validated together with one
  additional upcoming workflow feature in a single re-run (user pause point).
- Per-year values for scalar-exception signatures (e.g., the per-year E_i values behind
  `elasticity_static`) — not collected; they never pass through `generate_stats`. Can be
  added later if domain experts want them.

## 7. Implementation order

0. **Parquet2 pre-flight**: add `[compat]` pin for Parquet2 in `julia/Project.toml`;
   standalone write→read round-trip smoke check (typed 4-col frame + empty frame) to
   confirm the `writefile` API and codec kwarg for the resolved version.
1. `stats.jl`: `AnnualCollector` + hoist `value_cols` resolution (literal predicate,
   unsorted df) + collection hook (year-skip guard); export `AnnualCollector` from
   `StreamflowSignatures.jl`.
2. Thread `collector` kwarg through the 13 signature functions + `calculate_all_signatures`.
3. `config/signatures_config.json` (section present, `save: true`) + `config.jl`
   (`CFG_SAVE_ANNUAL_VALUES`, default **false** when absent) + add the constant to the
   `StreamflowSignatures.jl` export list.
4. `run_julia_benchmark.jl`: per-gage collector, global accumulators, sorted parquet
   write (skipped in the zero-gage branch), timing entry.
5. Tests (§4.1, incl. zero-warn and coverage regression); run `julia/test` suite.
6. Validation script (§4.2) — committed now, executed at the next benchmark run
   (which must also pass the golden-output regression).
7. Docs (§5).

## 8. Risks & mitigations

| Risk | Mitigation |
|---|---|
| Accidental behavior change while hoisting `value_cols` resolution | Literal-predicate hoist on unsorted df; unit test asserts with/without-collector Dict equality; golden-output regression at next benchmark run |
| Missing/NaN year or `missing` value rows breaking collector conversion (`Int32(round(NaN))` throws) | Explicit year-skip guard before conversion; `missing` value → NaN; unit test covers both |
| Silent collector breakage swallowed by per-signature `try/catch` (`signatures.jl:52-157`) | Zero-warn assertion + signature-coverage regression (§4.1) + summary/parquet coverage cross-check (§4.2) |
| Parquet2 write API/codec unverified (no compat pin, no Manifest, reads-only usage today) | Step 0 pre-flight: compat pin + write/read round-trip before integration; fall back to snappy/uncompressed |
| Double-counting if a signature name ever appeared in two generate_stats calls | Names verified globally unique across all call sites (Codex finding 9); unit test asserts no duplicate (gage, signature, water_year) keys |
| Stale/external configs silently enabling the new artifact | `CFG_SAVE_ANNUAL_VALUES` defaults to false when the config section is absent |
| Memory growth on 20M rows | ~500 MB peak, measured against existing 111M-row input load; acceptable |
| Sensitivity experiments overwriting the production annual file | `OUTPUT_PREFIX` naming already isolates them |

## 9. Codex review record (2026-07-14, gpt-5.x via codex-rescue — initial verdict NO-GO-as-is)

All five required amendments are incorporated above. Findings and dispositions:

| # | Severity | Finding | Disposition |
|---|---|---|---|
| 1 | CRITICAL | Choke-point claim overstated: callers prune before `generate_stats` (`qp_seasonality.jl:188-190`, `recession.jl:465-468`; sparse frames in flashiness/storage/elasticity/baseflow) | **Verified & accepted** — §2 rewritten with precise export semantics (export = exact series generate_stats consumed; absent row ≡ NaN); §3.5 schema note updated |
| 2 | MAJOR | **Pre-existing bug found**: `runoff_ratios.jl:113-118` looks up seasonal flags as `winter_complete/spring_complete/summer_complete/fall_complete` but the preprocessor emits `win_complete/spr_complete/sum_complete/fal_complete` (`io.jl:462`, cf. `flow_volumes.jl:136-137`) — the existence check at `runoff_ratios.jl:120` silently fails, so **incomplete seasons are never NaN'd out of seasonal runoff ratios**. Likely mirrored in Python/rpkg (synced ports). | **Verified — independent of this feature.** Logged in CHANGELOG as a discovered bug. Fixing it changes outputs → bundle with the next behavior-changing release + cross-language sync; do NOT fix silently inside this additive feature. §4.1 avoids runoff ratios as the NaN test case meanwhile |
| 3 | MAJOR | `CFG_SAVE_ANNUAL_VALUES` must be exported (`StreamflowSignatures.jl:84-91`) for bare use in the benchmark | **Accepted** — §3.6 + step 3 |
| 4 | MAJOR | `save: true` on absent config section is unsafe (config loaded once at import; no runtime override) | **Accepted** — default false when absent; repo config ships the section with true |
| 5 | MAJOR | Parquet2 write path unverified (no compat pin, no Manifest, reads-only today) | **Accepted** — step 0 pre-flight added |
| 6 | MAJOR | Year-coercion sketch throws on NaN instead of skipping | **Accepted** — §3.2 sketch fixed (guard before conversion) |
| 7 | MAJOR | Per-signature `try/catch` can silently swallow collector bugs | **Accepted** — zero-warn assertion + coverage regressions added (§4.1/§4.2) |
| 8 | MINOR | Hoist safe only if literal predicate on unsorted df; tests use explicit `year_col` | **Accepted** — §3.2 item 1 tightened; collector keys on passed `year_col` |
| 9 | MINOR | No signature-name collisions or alternate year columns found (independent verification) | Confirms design assumption |
| 10 | MINOR | "Summary NaN + no rows" also arises from caller-side pruning | **Accepted** — §4.2 rules enumerate all three paths |
| 11 | MINOR | Zero-gage branch (`run_julia_benchmark.jl:279-284`) behavior unspecified | **Accepted** — skip annual write + log (§3.4) |

## 10. Codex code review (implementation, 2026-07-14 — initial verdict NO-GO-as-is, fixes applied)

The implemented code was reviewed by a second Codex pass. The **core implementation was
confirmed correct**: hoisted `value_cols` resolution preserves the frozen contract;
`_collect_annual!` guards hold (missing/non-finite years skipped, missing values → NaN,
no df mutation); benchmark wiring correct on both legacy and non-legacy paths (one
collector per gage, one drain, sorted pre-metadata write, zero-gage skip); config default
+ exports correct; Parquet2 0.2.27 `writefile(...; compression_codec=:zstd)` verified
against the installed package source. All findings were in the test/validation harness:

| # | Severity | Finding | Fix applied |
|---|---|---|---|
| 1 | MAJOR | Validator only flagged `n_nonnan == 0` behind a non-NaN summary mean; 1–2 rows could slip through if mean/median happened to match | Coverage rule tightened to `n_nonnan >= MIN_ROWS (3)` — generate_stats only computes a mean from ≥3 values, all collected |
| 2 | MAJOR | Tests asserted signature presence, not completeness; parameterized-BFI collector path untested (the noisy synthetic gage produces ZERO recession events, so recession + param-BFI paths never collected) | Added `assert_collector_completeness` (recomputed mean/median must reproduce summary for every collected signature, ≥3 rows); exact-set coverage (`collected == EXPECTED_DENSE_SIGNATURES`); deterministic per-signature row counts; new deterministic exponential-recession gage (`Q_{t+1}=0.95·Q_t`) exercising all recession metrics + both param BFIs (25 rows each, alpha scalar exactly 0.95) with frozen-contract + zero-warn checks on that gage too |
| 3 | MINOR | Validator hard-failed the legitimate zero-gage path (missing parquet → exit 1) | Missing parquet + empty summary CSV → PASS (no-op) exit 0; missing parquet + populated summary remains a failure |

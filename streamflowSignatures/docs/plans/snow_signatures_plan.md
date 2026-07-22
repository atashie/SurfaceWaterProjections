# Plan: Snow Metrics Signature Family (Daymet SWE)

**Status**: IMPLEMENTED (2026-07-14, Julia). Codex plan review 2026-07-14: initial
NO-GO-as-is → all 7 findings + 2 delta-round residuals incorporated → **GO** (§12).
Full test suite green post-implementation (78 + 569 + 60 + 24 + 329 = 1,060 assertions).
Open items resolved at implementation: Daymet `swe` units confirmed **mm** by sample
read (first row group: max 561.7, mean 53.5, Maine gages); per-gage spatial support
(basin-average vs centroid) confirmed NOT recorded in-repo → documented as an
interpretation caveat in SIGNATURES.md §12. Benchmark re-run pending — bundled with
annual values + b=1 alpha.
**Requested by**: user (domain decision), 2026-07-14, with four clarifying decisions resolved
in-session (threshold, spell anchoring, melt rate, gating — see §1).
**Scope**: Julia canonical only; Python/rpkg ports deferred (join the existing July port
queue: annual-values collector + b=1 recession alpha). Additive columns only — the frozen
CSV contract is untouched (no renames/removals; the VALUES of one existing QA flag column,
`flagged_for_high_na`, will shift for some gages — accepted policy, see §6).
**Reference**: Hatchett, B.J. (2021). Seasonal and Ephemeral Snowpacks of the Conterminous
United States. *Hydrology* 8(1), 32. (SSM definition, §2.2; 60-day seasonal threshold
traces to Petersky & Harpold 2018.)

## 0. Summary

Add a 14-metric snow signature family computed from daily SWE, per water year per gage,
flowing through the standard `generate_stats()` path: 8 statistics + 8 Pettitt fields per
metric (**224 new columns**), plus automatic per-year rows in the annual-values parquet.

**Data source**: Daymet `swe` (kg/m² = mm water equivalent), already present in
`daymet_1980_2023.parquet` but currently dropped at the climate join
(`run_julia_benchmark.jl:150,239` select only `[:gage_id, :date, :PPT]`). Canadian gages
have no Daymet → all-NA snow columns (same as runoff ratios/elasticity/storage today).
Daymet SWE is a model product with documented biases (esp. mountain terrain) — carried as
a documentation caveat; upgrading to UA 4-km SWE / SNODAS is a separate future ingestion
project.

## 1. Resolved design decisions (user, 2026-07-14)

| # | Decision | Choice |
|---|----------|--------|
| 1 | Snow-day threshold | **SWE ≥ 10 mm**. Values < 10 mm are ignored for ALL duration AND magnitude calculations. Config-exposed. |
| 2 | Snow-on/off anchor spell | **The continuous spell containing the annual max SWE.** |
| 3 | Melt rate | **max SWE ÷ melt-season days** (net linear ablation between peak and snow-off). |
| 4 | Timing-metric gate | Subsumed by decision 1: any year with ≥ 1 qualifying snow day has max SWE ≥ 10 mm by construction. SSM's seasonal component additionally requires ≥ 60-day continuous spells (paper rule). |

**The thresholded series.** All metrics are computed from
`SWE*_t = (SWE_t >= THR) ? SWE_t : 0.0` with `THR = 10.0` mm. The snow-day mask is
`SWE* > 0`. This makes "values < 10 mm are ignored" uniform across durations, magnitudes,
spells, and melt increments. Consequences (intentional, documented):

- A year whose SWE peaks at 8 mm is operationally snow-free: `swe_max = 0`,
  `snow_cover_days = 0`, `swe_apr1 = 0` (valid zeros); all timing/melt metrics NaN.
- `swe_apr1` reports 0 when April 1 SWE is 2 mm (sub-threshold censoring at 10 mm).
- A sub-10 mm day between two deep spells **splits** them (two spells, not one).
- Melt increments are differences of SWE*, so the final melt-out step (crossing from
  ≥ 10 mm to < 10 mm) attributes the full remaining SWE to the crossing day, and melt
  below 10 mm is not counted.
- Our SSM therefore diverges from Hatchett's literal `SWE > 0` rule; the config knob
  allows a paper-faithful sensitivity run (threshold → 0).

## 2. Metric specifications

Per water year, on the complete daily grid the preprocessor guarantees for valid years.
Let `mask = SWE* .> 0`; `spells` = maximal runs of consecutive true days; melt increments
`m_t = max(0, SWE*_{t-1} − SWE*_t)` for t = 2..n (within-year only; melt across the
Sep 30 → Oct 1 boundary is not attributed — documented); `total_melt = Σ m_t`.
**Day-attribution pin (implementation-ready)**: `m_t` is attributed to day `t` — the day
the drop lands. `melt_before_peak` sums `m_t` for `t ≤ swe_max_dowy` (with the first-max
tie rule, `m_t` on the peak day itself is necessarily 0). `melt_com_dowy` is the FIRST
day `t` where `cumsum(m)[t] ≥ melt_com_fraction × total_melt`.

| # | Base name | Definition | No-snow year | Units |
|---|-----------|------------|--------------|-------|
| 1 | `swe_max` | `maximum(SWE*)` | **0.0 (valid)** | mm |
| 2 | `swe_max_dowy` | DOWY of first occurrence of the max (ties → first) | NaN | day |
| 3 | `snow_cover_days` | `count(mask)` | **0.0 (valid)** | days |
| 4 | `snow_on_dowy` | First DOWY of the anchor spell (spell containing `swe_max_dowy`). **Censored → NaN** if that spell starts on DOWY 1 (carryover snowpack predates Oct 1). | NaN | day |
| 5 | `snow_off_dowy` | First snow-free DOWY after the anchor spell (= last spell day + 1). **Censored → NaN** if the spell runs through the last day of the WY (snowpack persists past Sep 30). | NaN | day |
| 6 | `melt_season_days` | `snow_off_dowy − swe_max_dowy` (≥ 1 by construction) | NaN | days |
| 7 | `melt_rate` | `swe_max / melt_season_days` | NaN | mm/day |
| 8 | `ssm` | Spells ≥ 60 days are *seasonal*, < 60 days *ephemeral* (both config). `(seasonal_days − ephemeral_days) / snow_cover_days`, over ALL spells in the WY. NaN when `snow_cover_days == 0`. Range [−1, +1]. | NaN | – |
| 9 | `swe_apr1` | `SWE*` on **calendar April 1** (`Date(yr, 4, 1)`; leap-safe — DOWY 183 or 184) | **0.0 (valid)** | mm |
| 10 | `melt_before_peak` | `Σ m_t for t ≤ swe_max_dowy` (whole-year increments — early ephemeral melt-outs and mid-winter thaws count) | NaN | mm |
| 11 | `melt_before_peak_pct` | `100 × melt_before_peak / total_melt`; NaN if `total_melt ≤ 0` | NaN | % |
| 12 | `melt_before_peak_to_max_swe` | `melt_before_peak / swe_max` (denominator ≥ 10 mm whenever defined → bounded) | NaN | – |
| 13 | `melt_com_dowy` | First DOWY where cumulative melt ≥ 50% (config) of `total_melt`; NaN if `total_melt ≤ 0` | NaN | day |
| 14 | `swe_max_to_ppt` | `swe_max / Σ PPT(WY)`. Requires the year ∈ `valid_climate_years` AND `Σ PPT ≥ 10 mm` (mirrors the runoff-ratio annual floor); else NaN. | **0.0 (valid)** when PPT qualifies | – |

Notes:
- An anchor spell spanning the ENTIRE water year censors BOTH ends: `snow_on_dowy`,
  `snow_off_dowy`, `melt_season_days`, `melt_rate` all NaN; magnitude metrics unaffected.
- Metrics 6/7 are anchored to the max-SWE spell; metrics 8/10/11/13 pool ALL spells /
  whole-year increments; metric 10's window is "before the peak DOWY" regardless of spell.
- Melt-out day for a spell ending on day `e` contributes `m_{e+1} = SWE*_e` (drop to 0).
- Zero-vs-NaN policy: magnitude metrics (1, 3, 9, 14) emit **valid zeros** for
  operationally snow-free years → dense series at snow-free gages (constant-zero trends,
  same behavior class as D1_day); timing/melt/regime metrics (2, 4–8, 10–13) emit NaN →
  the existing trend-completeness gate (NOT exempted; 80% overall at time of writing —
  lowered to 60% on 2026-07-21 per revised guidelines, decade gates remain 80%) correctly
  suppresses trends at marginal-snow gages.
- Daily SWE differencing understates melt when snowfall and melt co-occur within a day —
  inherent to the product; documented as a limitation.

## 3. Data plumbing

### 3a. `read_parquet` normalization (`julia/src/io.jl`)
Add `swe → SWE` rename alongside the existing `prcp → PPT` (guard: only if `SWE` absent).

### 3b. Preprocessor (`julia/src/io.jl` → `preprocess_daily_data`)
Mirror the per-year PPT block (lines ~380–442) for SWE:
- `has_swe = "SWE" in names(gage_flow)`; include `:SWE` in the grid left-join select.
- Per year: NaN count ≤ `max_raw_na_swe` (default 30); reject if any negative SWE
  (`reject_negative_swe`, default true — negative SWE is a data error, and Daymet never
  emits it); linearly interpolate internal gaps ≤ `max_gap_swe` (default 3 — defensible
  for a smooth state variable); residual NaN → year invalid for snow.
- `swe_ok` years: write interpolated `SWE` back; push to a new **`valid_swe_years`** vector.
- Return NamedTuple gains `valid_swe_years` (additive field). Empty-frame branch gains
  the `SWE Float64[]` column when `has_swe`.
- **Gotcha**: the benchmark's `start_water_year` experiment filter reconstructs the
  NamedTuple field-by-field (`run_julia_benchmark.jl:164-171`) — must add
  `valid_swe_years` filtering there or experiments silently pass unfiltered years.
- In practice Daymet is gridded and complete, so `valid_swe_years ≈ valid_years ∩ years
  with Daymet coverage`; the parquet spans calendar 1980–2023 → WY 1980 and WY 2024 are
  partial → auto-invalid. Same boundary behavior as existing climate signatures.

### 3c. Benchmark (`docs/benchmarks/run_julia_benchmark.jl`)
- After climate load: `has_swe = has_climate && ("SWE" in names(climate))`. (The full
  parquet is already loaded with all columns — zero additional load cost.)
- Both climate joins (non-legacy line ~150, legacy line ~239): build the select list
  dynamically — `[:gage_id, :date, :PPT]` plus `:SWE` when `has_swe`.
- Non-legacy call site: whenever `"SWE" in names(gage_data)` (i.e., the gage had Daymet),
  pass `snow_data = gage_data[in.(gage_data.water_year, Ref(Set(pp.valid_swe_years))), :]`
  — **possibly a 0-row frame** when no year passes the SWE checks (the snow module then
  emits the full NaN-stat key set, keeping per-gage schema stable — §4.1). Gages with no
  SWE column at all (Canadian) pass `snow_data = nothing` → snow keys absent → `missing`
  fill at assembly. Also pass `snow_climate_years=pp.valid_climate_years`.
- Legacy path: pass `snow_data = gage_data` EXPLICITLY when the SWE column is present
  (unfiltered years — deprecated compat mode, mirroring how legacy passes unfiltered
  climate), with `snow_climate_years = nothing`.
- **No `nothing → gage_data` fallback anywhere** (Codex finding 1): `snow_data === nothing`
  means "do not run snow", full stop. This is what prevents Q-valid/SWE-invalid years from
  leaking into snow metrics via `pp.data`, which carries an SWE column (including NaN-y
  invalid years) once the preprocessor handles SWE.

### 3d. Orchestrator (`julia/src/signatures.jl`)
New kwargs on `calculate_all_signatures`: `snow_data::Union{Nothing,DataFrame}=nothing`,
`snow_climate_years::Union{Nothing,Vector{Int}}=nothing`. New block (standalone — NOT
inside `if has_climate`; gated on the caller passing `snow_data` explicitly. This
deliberately does NOT mirror the `cdata` `nothing → gage_data` fallback — Codex finding 1):

```julia
if snow_data !== nothing && "SWE" in names(snow_data)
    try
        merge!(results, calculate_snow_metrics(snow_data;
            valid_climate_years=snow_climate_years,
            trend_completeness=tc, decade_completeness=dc,
            changepoint=cp, collector=coll))
    catch e
        @warn "calculate_snow_metrics failed for gage $gage_id" exception=(e, catch_backtrace())
    end
end
```

Direct-API callers must opt in by passing `snow_data` (documented in the docstring); an
SWE column sitting inside `gage_data` is never used implicitly.

Gages without SWE simply lack the 224 keys; the results-assembly key-union
(`run_julia_benchmark.jl:339-350`) fills `missing` — identical to how Canadian gages get
NA climate signatures today.

## 4. New module: `julia/src/snow.jl`

`calculate_snow_metrics(df::DataFrame; valid_climate_years=nothing, trend_completeness,
decade_completeness, changepoint, collector)`:

1. Requires columns `date, water_year, dowy, SWE` (+ optional `PPT`). **Schema contract
   (explicit — Codex finding 5)**: the module ALWAYS constructs the annual frame with all
   14 typed metric columns — even for a 0-row `snow_data` — and calls `generate_stats`
   with the explicit 14-column `value_cols` list, so every gage that reaches this function
   emits the identical 224-key set (NaN-filled where not computable). Gages that never
   reach it (no SWE data at all) lack the keys and get `missing` at assembly. Two tiers,
   both schema-stable, both stated.
2. Per water year present in `df`: build `SWE*`, the mask, spells (single pass), the 14
   metrics per §2. Metric 14: when `valid_climate_years !== nothing`, require membership
   AND the ≥ 10 mm annual-PPT floor; when `nothing` (legacy fallback), mirror the
   runoff-ratio convention (`runoff_ratios.jl`): sum PPT over the year's non-NaN PPT days
   and apply the same ≥ 10 mm floor — NOT a stricter fully-non-NaN requirement (Codex
   finding 4).
3. Assemble a **dense** annual DataFrame (one row per water year in `df`, NaN where not
   computable) → one `generate_stats` call over the 14 value columns → dense annual-export
   semantics (NaN placeholder rows), like flow volumes. Years Q-valid but SWE-invalid are
   absent from `snow_data` → omitted rows (≡ NaN for consumers), like climate signatures.
4. `seasonal_flags` not used (no seasonal sub-metrics). Trend/decade completeness gates
   APPLY (snow metrics are dense per-year values — no recession/elasticity-style exemption).
5. Register: `include("snow.jl")` + `export calculate_snow_metrics` in
   `StreamflowSignatures.jl` (18th module).

## 5. Configuration

`config/signatures_config.json`:

```json
"snow": {
  "swe_day_threshold_mm": 10.0,
  "seasonal_spell_min_days": 60,
  "melt_com_fraction": 0.5,
  "min_annual_ppt_mm": 10.0,
  "comment": "Snow metrics from Daymet SWE. Days with SWE < threshold are treated as snow-free everywhere (magnitudes AND durations, per domain decision 2026-07-14). seasonal_spell_min_days per Hatchett (2021) Hydrology 8(1):32 / Petersky & Harpold (2018). Set threshold near 0 to reproduce the paper's literal SWE>0 rule."
}
```

`na_handling` additions — as a **nested block mirroring the existing `climate_na_policy`
group** (the Julia loader reads nested groups via `get(..., Dict())`; the originally
proposed flat keys would have been a new, easy-to-miswire pattern — Codex finding 3):

```json
"na_handling": { ..., "snow_na_policy": { "max_raw_na": 30, "max_gap_days": 3, "reject_negative": true } }
```

Key names inside the block mirror whatever `climate_na_policy` actually uses — confirm at
implementation and stay symmetric.

`julia/src/config.jl`: `CFG_SNOW_SWE_THRESHOLD_MM`, `CFG_SNOW_SEASONAL_MIN_DAYS`,
`CFG_SNOW_MELT_COM_FRACTION`, `CFG_SNOW_MIN_ANNUAL_PPT_MM`, `CFG_NA_MAX_RAW_NA_SWE`,
`CFG_NA_MAX_GAP_SWE`, `CFG_NA_REJECT_NEGATIVE_SWE` — all with absent-section defaults
matching the values above (function must work when the `snow` section is missing).

`config.R`: append the 14 base names to `EXPECTED_SIGNATURE_BASES`.

## 6. Statistics / export integration (all automatic via `generate_stats`)

- 8 standard statistics + 8 Pettitt changepoint fields per metric → 14 × 16 = **224
  columns**; time-series signature count 76 → **90**; total columns 1,264 → **1,488**.
- Annual-values parquet: 14 new signatures collected per SWE-valid gage-year
  (≈ +2.5 M rows for ~6,000 Daymet gages × ~30 yr); `validate_annual_values.py` covers
  them with no changes (dense NaN-placeholder semantics).
- Changepoint window (1980–2024) vs SWE coverage (WY 1981–2023): edge years NaN — the
  existing `min_total_obs=20` handles it.
- **Known existing-column VALUE change (accepted policy — Codex finding 2)**:
  `flagged_for_high_na` computes its NA fraction over ALL signature columns
  (`julia/src/qa_qc.jl:183-189`, threshold `CFG_QAQC_MAX_NA_FRACTION`). The 224 new snow
  columns enter that denominator, and no-SWE gages (Canadian) carry them as NA — so some
  gages' flag values WILL flip at the bundled re-run. Precedent: the April 2026 Pettitt
  addition (+608 columns) already reshaped this denominator identically. Policy: accept +
  document; whitelist `flagged_for_high_na` as expected-divergent in the new-vs-golden
  comparison (alongside the log_a columns); QUANTIFY the number of flipped gages in the
  re-run report. (Rejected alternative — availability-aware denominators — would
  retroactively change the flag's existing semantics for the climate columns too.)

## 7. Tests (`julia/test/test_snow_metrics.jl`, + runtests include)

Deterministic synthetic gages with hand-derivable or independently-recomputed expectations:

1. **Triangular snowpack gage** (linear accumulation Nov→Mar, linear melt to May, 20+
   yrs): exact `swe_max`, `swe_max_dowy`, threshold-crossing `snow_on/off_dowy`,
   `melt_season_days`, `melt_rate = swe_max/length`, `swe_apr1` from geometry,
   `melt_before_peak == 0`, `pct == 0`, `ratio == 0`, `ssm == +1`, `melt_com_dowy`
   against an independent in-test recompute (threshold clipping makes hand formulas
   fiddly — recompute helper is authoritative, plus a few hand-derived anchors).
2. **Mid-winter-thaw gage** (two peaks, deep thaw between): exact `melt_before_peak > 0`,
   `pct`, `ratio`, earlier `melt_com_dowy`; anchor spell = the spell containing the
   LARGER peak.
3. **Spell-arithmetic gage** (one 90-day + two 15-day spells): `ssm == (90−30)/120 == 0.5`
   exactly; `snow_cover_days == 120`; ephemeral-only variant → `ssm == −1`.
4. **Sub-threshold censoring gage** (peak 8 mm): `swe_max == 0.0`, `snow_cover_days == 0`,
   `swe_apr1 == 0.0`, all timing/melt/ssm NaN.
5. **Spell-splitting**: a single 9.9 mm day inside a deep spell splits it into two spells
   (SSM arithmetic changes accordingly); a 10.0 mm day does not.
6. **Censoring**: SWE ≥ 10 on Oct 1 → `snow_on_dowy` NaN while max/off/melt compute;
   SWE ≥ 10 through Sep 30 → `snow_off/melt_season/melt_rate` NaN while `snow_on`/peak
   compute.
7. **Tie rule**: two days at the max → first DOWY.
8. **Leap year**: `swe_apr1` picks `Date(yr,4,1)` (DOWY 184 in WY 2004), value exact.
9. **Metric 14**: exact ratio under constant PPT; a PPT-invalid year (∉
   `valid_climate_years`) → NaN despite being SWE-valid; `ΣPPT < 10 mm` → NaN.
10. **Preprocessor**: >30 SWE NaNs → year ∉ `valid_swe_years` (but still ∈ `valid_years`);
    ≤3-day gap interpolated (exact linear values); negative SWE rejects; `valid_swe_years ⊆
    valid_years`; NamedTuple field present and correctly filtered by `start_water_year`
    reconstruction.
11. **Harness invariants** (patterns from the annual-collector suite): with/without
    collector stat-identity (`isequal`); zero-warnings through `calculate_all_signatures`
    with `snow_data`; exactly 16 keys per base (8 stats + 8 Pettitt); collector rows for
    all 14 signatures with expected counts and NaN placeholders.
12. **No-fallback gate** (Codex findings 1/5): `calculate_all_signatures` on a gage frame
    CONTAINING an SWE column but with `snow_data=nothing` emits ZERO snow keys (pins the
    orchestrator against any `nothing → gage_data` leak); a 0-row `snow_data` emits the
    full 224-key NaN set.
13. **Filtered-years leak test**: a gage with one SWE-corrupt year (> 30 SWE NaNs) — that
    year's SWE values must not influence any snow statistic (equality against the same
    gage with the year hard-removed).
14. **Legacy metric-14 consistency** (finding 4): with `valid_climate_years=nothing`, a
    year with scattered PPT NaNs reproduces the runoff-ratio convention (non-NaN-day sum
    + ≥ 10 mm floor), not a stricter all-non-NaN rule.
15. **Real-data harness exercises snow** (finding 6; revised at delta-verify):
    `julia/test/smoke_test.jl` is THE real-data snow harness — add `:SWE` to its climate
    join select (line ~94) and assert presence + sane ranges of the 224 keys for a snowy
    USGS gage. `julia/test/test_against_golden.jl` is explicitly NOT extended: it loads
    climate but never joins it into the gage frames (it exercises no climate signatures
    today), and golden outputs contain no snow columns to compare against — extending it
    would change that harness's purpose. Documented decision, not an oversight.
16. **Anchor-spell corner regressions** (finding 7): peak on the spell's LAST day
    (melt season = 1 day, `melt_rate == swe_max`); anchor spell spanning the full WY
    (both ends censored → timing NaN, magnitudes valid); `m_t` day-attribution
    off-by-one pin via a hand-built 5-day melt sequence with exact `melt_before_peak`
    and `melt_com_dowy`.
17. **QA-flag semantics pin** (finding 2; corrected at delta-verify): `compute_qa_flags`
    on a synthetic results table where the 224 snow columns are PRESENT (as they are
    after cross-gage key-union at assembly) and hold `missing`/NaN for a no-SWE gage row
    — assert they DO enter the NA-fraction denominator and can flip the flag. Columns
    absent from the table never enter the denominator (`qa_qc.jl` builds its column list
    from present columns), so the originally-worded missing-columns test would have
    validated the wrong semantics.

## 8. Docs

- `docs/SIGNATURES.md`: new "12. Snow Metrics" section (before the Summary Table): source
  + caveats (Daymet model SWE; Canadian gages NA; 10 mm operational threshold; daily-
  differencing melt limitation), metric table, references (Hatchett 2021; Petersky &
  Harpold 2018; Musselman et al. 2017 for melt-rate context); Summary Table row; update
  column counts.
- `CHANGELOG.md` [July 2026]: new-signature entry (decision log, column arithmetic,
  Julia-only + port deferral); note definitions originated from the user in-session,
  pending formal addition to the guidelines Google Doc by the domain team.
- `docs/DEVELOPMENT.md`: data-flow diagram (+SWE), NA-handling section (+`valid_swe_years`),
  benchmark notes.
- `claude-skill/streamflow-signatures.md`: snow family + interpretation notes.
- `CLAUDE.md`: signature stats table unchanged; add snow to relevant counts if listed.

## 9. Validation & rollout

1. Implement (Julia) → full test suite green (existing 78 + 569 + 60 + new snow suite).
2. Codex review of the implementation (house pattern), fixes, delta-verify.
3. **Single bundled benchmark re-run** validates all three July features together
   (annual export + b=1 alpha + snow). Expected-divergence whitelist vs golden: (a) the
   log_a / recession-seasonality columns (b=1 alpha change — intentional), (b)
   `flagged_for_high_na` (denominator now includes the 224 snow columns — §6 policy;
   quantify the flips). Every other column shared with golden must NOT diverge. Run
   `validate_annual_values.py`; new-vs-golden dashboard for the shared columns.
4. Ports to Python/rpkg: deferred, three-feature queue.

## 10. Open items / risks

| Item | Handling |
|------|----------|
| Daymet SWE provenance (basin-averaged vs centroid point per gage) | Verify from the R conversion lineage (`R/run_conversion.R` / metadata docs) during implementation; document either way in SIGNATURES.md. Interpretation caveat only. |
| Daymet `swe` units (kg/m² = mm) | Confirm on a sample read at implementation; assert plausible range in QA. |
| Daymet SWE model bias (esp. mountains, deep-pack underestimation) | Documentation caveat; future upgrade path = UA SWE / SNODAS ingestion (separate project). |
| WY 1980 / WY 2024 partial Daymet coverage | Auto-rejected by `valid_swe_years`; matches existing climate-signature behavior. |
| Melt from daily ΔSWE understates true melt under same-day snowfall+melt | Documented limitation (inherent to daily SWE). |
| Sub-threshold censoring step at 10 mm creates a discontinuity in magnitude metrics at marginal gages | Intentional per domain decision; config knob enables sensitivity runs. |
| `n_modules` grows to 18; collector/changepoint threading must match the other 13 call sites | Covered by zero-warnings + collector-coverage tests. |
| `flagged_for_high_na` values shift as 224 NA-heavy columns enter its denominator (esp. Canadian gages) | Accepted policy (§6, Pettitt precedent); whitelisted + quantified at the re-run; unit test pins the semantics (§7.17). |

## 11. Out of scope

- New SWE data ingestion (UA SWE, SNODAS, ERA5-Land) and Canadian SWE coverage.
- Caravan-pipeline snow metrics (Caravan bundles ERA5 SWE — possible later, R-side).
- Snow-fraction-of-precipitation metrics (would need temperature partitioning — future).
- Python/rpkg ports (deferred with the other two July features).

## 12. Codex plan-review record (2026-07-14)

Initial verdict: **NO-GO-as-is** — "one real control-flow bug in the non-legacy path and
one unaccounted-for existing-output regression." All 7 findings incorporated into this
document (job task-mrl3hrwe-fejgq1; a first review attempt, task-mrkyox9c-u8jwxn, hung
after 47 min of healthy reading and was cancelled/relaunched):

| # | Severity | Finding | Disposition |
|---|----------|---------|-------------|
| 1 | CRITICAL | `snow_data=nothing` + a `cdata`-style fallback to `gage_data` would leak Q-valid/SWE-invalid years into snow metrics (pp.data carries an SWE column once the preprocessor handles SWE) | Fallback REMOVED: snow runs iff `snow_data !== nothing`; benchmark passes a possibly-empty filtered frame whenever the gage has an SWE column; legacy passes `gage_data` explicitly (§3c/§3d). Regression tests §7.12/§7.13 |
| 2 | CRITICAL | "No golden divergence possible" was false: `flagged_for_high_na` (qa_qc.jl:183-189) counts ALL signature columns, so 224 new NA fields flip the flag for some no-SWE gages | Policy: accept (April 2026 Pettitt +608-column precedent), whitelist in the new-vs-golden comparison, quantify flips at the re-run, pin semantics with a unit test (§0/§6/§9/§7.17). qa_qc.jl claim independently re-verified in-session |
| 3 | MAJOR | Proposed flat `na_handling` keys didn't match the loader's nested-group pattern (`interpolation`/`year_rejection`/`climate_na_policy`) | Nested `na_handling.snow_na_policy` block mirroring `climate_na_policy` key style (§5) |
| 4 | MAJOR | Legacy metric-14 fallback ("fully non-NaN PPT") was STRICTER than the runoff-ratio convention it claimed to mirror (runoff ratios sum valid days + floor) | Legacy fallback now mirrors runoff_ratios.jl: non-NaN-day PPT sum + ≥ 10 mm floor (§4.2); consistency test §7.14 |
| 5 | MAJOR | Schema contract internally inconsistent between §3d ("may lack keys") and §4 ("full NaN key set") | Explicit two-tier contract: SWE column present → ALWAYS all 224 keys via explicit `value_cols` (even 0-row input); no SWE data at all → keys absent, `missing` fill (§4.1) |
| 6 | MAJOR | Test plan missed the four likeliest failures: fallback leak, QA-flag regression, legacy-14 strictness, and real-data harnesses (smoke_test.jl:94, test_against_golden.jl) still dropping SWE at their own joins | Tests §7.12–§7.17 added; smoke_test.jl + test_against_golden.jl SWE-join updates brought into scope. smoke_test.jl claim independently re-verified in-session |
| 7 | MINOR | Melt indexing (`m_t` day attribution) and the full-WY anchor-spell case not code-ready | Day-attribution pin + full-WY both-ends-censored rule added (§2); corner regressions §7.16 |

Review's verified-claims summary: read_parquet normalization set, preprocessor NamedTuple
+ benchmark field-by-field reconstruction, both climate joins dropping SWE, the
orchestrator climate-fallback pattern, the generate_stats/AnnualCollector choke point,
and config absent-section tolerance (an unknown top-level `snow` section cannot break
existing consumers that read specific keys) — all confirmed against code. Not verifiable
in-repo: Daymet parquet swe contents/units/provenance (external file — §10 open item).

**Delta-verification round 1**: findings 1–5 and 7 RESOLVED; finding 6 partially — the
review corrected §7.15 (`test_against_golden.jl` never joins climate at all, so "add SWE
to its join" was the wrong fix; smoke_test.jl is now the sole real-data snow harness,
with test_against_golden.jl explicitly out of scope) — and flagged one NEW problem: §7.17
originally tested `compute_qa_flags` with snow columns ABSENT, which would validate the
wrong semantics (absent columns never enter the denominator); corrected to
present-columns-holding-missing. Both §7.15 and §7.17 amended accordingly.

**Delta-verification round 2**: §7.15 RESOLVED, §7.17 RESOLVED, §12 record accurate, no
new issues. **Final verdict: GO** (job task-mrl491hg-9q4067, 2026-07-14).

## 13. Codex implementation-review record (2026-07-14)

Post-implementation review (same Codex session, job task-mrl6adep-flsef1). Plan
findings 1/3/4/7 verified CONFIRMED-IN-CODE (no-fallback control flow traced in
signatures.jl + benchmark; nested config + absent-section defaults; metric-14 legacy
runoff-ratio parity; melt attribution/tie/censoring). Initial verdict
**NO-GO-as-is** — 3 findings, all fixed:

| # | Severity | Finding | Fix applied |
|---|---|---|---|
| 1 | MAJOR | The defensive no-SWE-column path used `empty_stats()` (8 keys/metric) — with changepoint enabled it emitted 112 keys instead of the contractual 224 (plan finding 5 violated on that one path; suite only tested it without changepoint) | Early return rewritten to build a 0-row 14-column frame and call `generate_stats` with explicit `value_cols` — identical machinery to the 0-row case, so the key set is contract-identical on every path. Regression test added (no-SWE + changepoint → 224 NaN keys) |
| 2 | MAJOR | Complete-grid guard checked only endpoints (`n>=365`, `dowy[1]==1`, `dowy[end]==n`) — a duplicated interior day masking a missing one passes and computes wrong metrics (direct-API/malformed-legacy inputs only; preprocessed path unaffected) | Guard now requires the EXACT sequence: `n ∈ {365,366}`, no missing dowy, `dowy[i] == i ∀ i`. Regression test added (duplicate day 150 + drop day 151 → all-NaN year; intact control computes) |
| 3 | MINOR | smoke_test.jl exercised only the legacy-style unfiltered snow path — the canonical `valid_swe_years` plumbing had no real-data coverage | smoke_test now preprocesses per gage and passes the `valid_swe_years`-filtered frame (benchmark-identical provenance); rest of the harness stays legacy-style |

**Delta-verification** (job task-mrl6nmcr-yn17su): all three RESOLVED, "NEW problem
introduced: none found." **Final verdict: GO.** Post-fix suite: snow 377/377; full
package suite green (78 + 569 + 60 + 24 + 377 = 1,108 assertions).

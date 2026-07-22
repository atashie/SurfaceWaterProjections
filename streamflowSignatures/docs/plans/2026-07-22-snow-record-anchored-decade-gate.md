# Snow Record-Anchored Decade Gate — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Suppress trend statistics for the 10 threshold-dependent snow metrics at gages whose computable (snowy) years are absent from the first or last decade of the gage's SWE-valid record — closing the clustered-snowy-years hole that the metric-own-span trend gate and the 20-value stats floor leave open.

**Architecture:** A new opt-in `force_skip_trends` kwarg on `generate_stats()` reuses the existing skip-trends path (6 trend stats → NaN; mean/median computed; Pettitt fields still run; collector untouched). `calculate_snow_metrics()` computes the gate per metric — anchored to the gage's SWE-valid, grid-complete years (rows where `swe_max` is non-NaN) — and passes the failing metrics in that set. The threshold is **the same `decade_completeness` value the orchestrator already passes** (`CFG_NA_DECADE_MIN_FRACTION`, currently 0.80), so the SWE gate and the streamflow decade gate are structurally linked: one config knob governs both.

**Tech Stack:** Julia only (snow metrics are Julia-only; the gate ports to Python/rpkg with the queued snow port). Config: `config/signatures_config.json` + `rpkg/inst/config/signatures_config.json`.

---

## User decisions (2026-07-22, this session)

1. **Threshold**: "match the streamflow decade gate … in any case, they should be linked" → implement as linked to `decade_min_fraction` (no new threshold constant). ⚠️ **OPEN QUESTION for domain experts**: the user twice referred to the decade gate as 60%, but config, guidelines doc, and manuscript §2.2.3 all say **80%** (yesterday's change lowered only the overall gate to 60%). Because the gate is linked, flipping decades to 60% later is a one-line config change (`decade_min_fraction: 0.80 → 0.60`) — but the guidelines doc and manuscript must be updated to match. Logged in CHANGELOG (Task 6).
   **→ RESOLVED same-day (user, 2026-07-22): 80% first/last decade + 60% overall confirmed.** Shipped config already matches; no change needed.
2. **Anchor**: the gage's SWE-valid record (first/last decade of SWE-valid years), NOT the metric's own non-NA span and NOT the analysis window.
3. **Scope**: only the 10 threshold-dependent timing/melt/regime metrics: `swe_max_dowy`, `snow_on_dowy`, `snow_off_dowy`, `melt_season_days`, `melt_rate`, `ssm`, `melt_before_peak`, `melt_before_peak_pct`, `melt_before_peak_to_max_swe`, `melt_com_dowy`. The 4 magnitude metrics (`swe_max`, `snow_cover_days`, `swe_apr1`, `swe_max_to_ppt`) emit valid zeros and are NOT gated.
4. **On failure**: match streamflow trend-gate semantics — NaN the 6 trend statistics only; mean/median and Pettitt changepoint fields still computed.

## Design details

- **Anchor-year set**: `annual_data` rows where `swe_max` is non-NaN. `swe_max` is dense over grid-complete SWE-valid years (0.0 for snow-free years), so this equals "years where snow metrics were computable at all" and excludes both SWE-invalid years (absent from `snow_data`) and defensive grid-guard NaN rows.
- **Windows**: first decade `[rec_min, rec_min+9]`, last decade `[rec_max−9, rec_max]` where `rec_min`/`rec_max` are the anchor-set extremes. Record span < 10 → gate skipped entirely (mirrors the existing decade-check skip for short series).
- **Fractions**: per metric, per window: `non-NaN metric years in window ÷ anchor years in window`. Denominator is anchor years (not calendar 10) — a year with no SWE data must not count against snow presence. Either window below `decade_completeness` → metric gated.
- **Activation**: requires BOTH the config flag `snow.record_anchored_decade_gate` (absent → `false`, legacy) AND a non-`nothing` `decade_completeness` kwarg (mirrors the opt-in trend gate). The orchestrator already passes `dc` to `calculate_snow_metrics` (`julia/src/signatures.jl:181`), so production activates automatically once the flag ships `true`.
- **Additivity**: the existing own-span trend gate, decade checks, and the 20-value stats floor all still apply — this gate only ADDS suppressions ("additionally require", per user).
- **Output impact**: at the next benchmark, gated gages flip 6 trend stats → NaN per gated metric. `flagged_for_high_na` denominators shift accordingly (same accepted precedent as Pettitt/snow launches). Annual-values parquet is unaffected (collection precedes gating).

---

### Task 1: Config flag (both JSONs + Julia constant)

**Files:**
- Modify: `config/signatures_config.json` (snow section, ~line 98)
- Modify: `rpkg/inst/config/signatures_config.json` (same section)
- Modify: `julia/src/config.jl:128` (after `CFG_SNOW_MIN_ANNUAL_PPT_MM`)

**Step 1: Add the flag to `config/signatures_config.json`** — inside the `"snow"` object, after `"min_annual_ppt_mm": 10,`:

```json
    "record_anchored_decade_gate": true,
```

and extend the section's `"comment"` with: `" record_anchored_decade_gate (July 2026): threshold-dependent snow metrics require >= decade_min_fraction (na_handling.trend_completeness — LINKED to the streamflow decade gate) of the first and last decade of the gage's SWE-valid record to have computable (snowy) values before trend statistics are computed; absent key => disabled."`

**Step 2: Mirror the identical edit in `rpkg/inst/config/signatures_config.json`.**

> **[EXECUTION NOTE 2026-07-22: Step 2 was a NO-OP.]** rpkg's bundled config has no
> `snow` section at all (snow metrics are not yet ported to rpkg), so there was
> nothing to mirror. The flag ships to rpkg's config together with the queued snow
> port. Codex review (2026-07-22, GO verdict) flagged the stale file list here.

**Step 3: Add the constant in `julia/src/config.jl`** after line 128:

```julia
const CFG_SNOW_RECORD_DECADE_GATE = Bool(get(_snow_config, "record_anchored_decade_gate", false))
```

Check `julia/src/StreamflowSignatures.jl` — the other `CFG_SNOW_*` constants' export treatment; match it (if they are not exported, do not export this one; tests reference it as `StreamflowSignatures.CFG_SNOW_RECORD_DECADE_GATE`).

**Step 4: Verify** — `julia --project=julia -e "using StreamflowSignatures; println(StreamflowSignatures.CFG_SNOW_RECORD_DECADE_GATE)"` → `true`.

Do NOT commit (project rule: commits only on user request).

---

### Task 2: `force_skip_trends` kwarg on `generate_stats` (TDD)

**Files:**
- Create: `julia/test/test_snow_record_decade_gate.jl`
- Modify: `julia/src/stats.jl:321-459`
- Modify: `julia/test/runtests.jl:288` (add include AFTER `test_stats_floor.jl`)

**Step 1: Write the failing unit tests** — create `julia/test/test_snow_record_decade_gate.jl`:

```julia
# Tests for the snow record-anchored decade gate (plan:
# docs/plans/2026-07-22-snow-record-anchored-decade-gate.md).
# Part 1: generate_stats force_skip_trends mechanics.
# Part 2: calculate_snow_metrics record-anchored gate behavior.

using Test
using DataFrames
using Dates
using Statistics
using StreamflowSignatures

const RDG_CP = (start_year=1980, end_year=2024, min_total_obs=20, min_segment_obs=10)

"Daily SWE gage: snowy years get a triangle profile (peak 242 mm), others all-zero."
function rdg_gage(years; snowy::Function, ppt=2.0)
    tri(t) = 32 <= t <= 152 ? 2.0 * (t - 31) : (t > 152 ? max(0.0, 242.0 - 3.0 * (t - 152)) : 0.0)
    dfs = DataFrame[]
    for yr in years
        dates = collect(Date(yr - 1, 10, 1):Day(1):Date(yr, 9, 30))
        n = length(dates)
        swe = snowy(yr) ? Float64[tri(t) for t in 1:n] : zeros(Float64, n)
        push!(dfs, DataFrame(gage_id="rdgtest", date=dates, Q=fill(1.0, n),
                             water_year=fill(yr, n), month=month.(dates),
                             dowy=collect(1:n), PPT=fill(Float64(ppt), n), SWE=swe))
    end
    return vcat(dfs...)
end

@testset "Snow record-anchored decade gate" begin

    @testset "generate_stats force_skip_trends mechanics" begin
        df = DataFrame(water_year=2000:2024, a=collect(1.0:25.0), b=collect(2.0:2.0:50.0))

        # Baseline: no kwarg vs empty-set-equivalent (nothing) — identical
        base = generate_stats(df; value_cols=["a", "b"], changepoint=RDG_CP)
        same = generate_stats(df; value_cols=["a", "b"], changepoint=RDG_CP,
                              force_skip_trends=nothing)
        @test all(isequal(base[k], same[k]) for k in keys(base))

        # Skip "a" only: 6 trend stats NaN; mean/median computed; Pettitt still runs
        skipped = generate_stats(df; value_cols=["a", "b"], changepoint=RDG_CP,
                                 force_skip_trends=Set(["a"]))
        for sfx in ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                    "_mk_rho", "_mk_pval"]
            @test isnan(skipped["a$(sfx)"])
            @test !isnan(skipped["b$(sfx)"])   # untouched column unaffected
        end
        @test skipped["a_mean"] ≈ mean(1.0:25.0)
        @test skipped["a_median"] ≈ median(1.0:25.0)
        @test !isnan(skipped["a_pettitt_pval"])          # changepoint independent
        @test isequal(skipped["b_mean"], base["b_mean"]) # b fully identical
        @test isequal(skipped["b_pettitt_pval"], base["b_pettitt_pval"])

        # Forced skip works WITHOUT trend_completeness passed (orthogonal paths)
        solo = generate_stats(df; value_cols=["a"], force_skip_trends=Set(["a"]))
        @test isnan(solo["a_senn_slp"]) && !isnan(solo["a_mean"])
    end

end
```

**Step 2: Register the file** — in `julia/test/runtests.jl`, after `include("test_stats_floor.jl")` add:

```julia
include("test_snow_record_decade_gate.jl")
```

**Step 3: Run to verify failure** —
`julia --project=julia -e "using Pkg; Pkg.test()"` (or targeted: `julia --project=julia test/runtests.jl` from `julia/`; match how the suite is normally run).
Expected: FAIL — `generate_stats` has no `force_skip_trends` kwarg (MethodError/UndefKeywordError).

**Step 4: Implement in `julia/src/stats.jl`** — three coordinated edits:

(a) Add the kwarg to the signature (line 321-331):

```julia
    min_values_for_stats::Union{Nothing, Int} = nothing,
    force_skip_trends::Union{Nothing, Set{String}} = nothing
```

(b) Restructure the gating block (lines 403-459). Replace:

```julia
        compute_trends = true
        if trend_completeness !== nothing
```

with:

```julia
        compute_trends = true

        # Externally-forced trend skip (July 2026): caller-computed gates (e.g. the
        # snow record-anchored decade gate) suppress the 6 trend stats through the
        # SAME path as trend_completeness — mean/median and changepoint unaffected.
        if force_skip_trends !== nothing && col in force_skip_trends
            compute_trends = false
        end

        if compute_trends && trend_completeness !== nothing
```

(c) Move the `if !compute_trends ... continue end` block (lines 443-458) OUT of the `trend_completeness !== nothing` branch so it executes for forced skips too — i.e., de-indent it to sit after the (now-guarded) trend_completeness block, at the same level as `compute_trends = true`. The block's contents are unchanged. (When `force_skip_trends === nothing` and `trend_completeness === nothing`, `compute_trends` stays `true` and the block is dead — byte-identical legacy behavior.)

(d) Update the docstring (Parameters section, ~line 311): add

```
force_skip_trends : Set{String} or nothing
    Columns whose 6 trend statistics are forced to NaN regardless of their own
    completeness (mean/median and changepoint fields still computed). Used by the
    snow record-anchored decade gate (July 2026).
```

**Step 5: Run the new test file** — expected: Part 1 PASSES. Also run `julia/test/test_stats_floor.jl`, `test_annual_collector.jl`, `test_changepoint.jl`, `test_snow_metrics.jl` — expected: all PASS (no behavior change with the kwarg absent).

---

### Task 3: Record-anchored gate in `calculate_snow_metrics` (TDD)

**Files:**
- Modify: `julia/test/test_snow_record_decade_gate.jl` (append Part 2)
- Modify: `julia/src/snow.jl:31-46` (new constant) and `:257-263` (gate + pass-through)

**Step 1: Append the failing behavior tests** to the `@testset` in `test_snow_record_decade_gate.jl` (inside the outer testset, after Part 1):

```julia
    # Metrics gated by the record-anchored decade gate (the 10 threshold-dependent)
    GATED = ["swe_max_dowy", "snow_on_dowy", "snow_off_dowy", "melt_season_days",
             "melt_rate", "ssm", "melt_before_peak", "melt_before_peak_pct",
             "melt_before_peak_to_max_swe", "melt_com_dowy"]

    @testset "gate fires when snow vanishes from the last decade" begin
        # 30 SWE-valid years 1981-2010; snowy 1981-2000 (20 -> passes stats floor),
        # snow-free 2001-2010. Own-span gate (1981-2000, 100% complete) and floor
        # both pass -> without the new gate every trend WOULD be computed.
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2000)
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                   min_values_for_stats=20, changepoint=RDG_CP)
        for m in GATED
            @test isnan(r["$(m)_senn_slp"])
            @test isnan(r["$(m)_mk_pval"])
        end
        @test !isnan(r["snow_on_dowy_mean"])          # mean survives (over snowy years)
        @test r["snow_on_dowy_mean"] ≈ 32.0            # triangle onset, hand-derived
        @test !isnan(r["snow_on_dowy_pettitt_pval"])  # Pettitt survives (20 obs)
        # Magnitude metrics are NOT gated (dense zeros carry the decline signal)
        for m in ["swe_max", "snow_cover_days", "swe_apr1"]
            @test !isnan(r["$(m)_senn_slp"])
        end
    end

    @testset "control: gate off when decade_completeness not passed" begin
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2000)
        r = calculate_snow_metrics(df; trend_completeness=0.6,
                                   min_values_for_stats=20, changepoint=RDG_CP)
        @test !isnan(r["snow_on_dowy_senn_slp"])  # pre-gate behavior: trend computed
    end

    @testset "gate fires on first-decade absence" begin
        df = rdg_gage(1981:2010; snowy = yr -> yr >= 1991)
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                   min_values_for_stats=20)
        @test isnan(r["snow_on_dowy_senn_slp"])
    end

    @testset "denominator is SWE-valid years, not calendar years" begin
        # Last decade [2001,2010]: years 2004-2006 have NO SWE data (absent rows,
        # i.e. SWE-invalid). Of the 7 valid years, 6 are snowy (2010 snow-free):
        # 6/7 = 0.857 >= 0.8 -> PASSES. A calendar-10 denominator (6/10 = 0.6)
        # would wrongly fire. First decade fully snowy.
        yrs = [collect(1981:2003); collect(2007:2010)]
        df = rdg_gage(yrs; snowy = yr -> yr != 2010)
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                   min_values_for_stats=20)
        @test !isnan(r["snow_on_dowy_senn_slp"])
    end

    @testset "linked threshold: gate obeys the passed decade_completeness" begin
        # Last decade: 7 of 10 snowy -> frac 0.7. dc=0.6 passes; dc=0.8 fires.
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2000 || yr in (2001, 2003, 2005, 2007, 2009, 2002, 2004))
        r60 = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.6,
                                     min_values_for_stats=20)
        r80 = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                     min_values_for_stats=20)
        @test !isnan(r60["snow_on_dowy_senn_slp"])
        @test isnan(r80["snow_on_dowy_senn_slp"])
    end

    @testset "short record (<10-year span): gate skipped" begin
        df = rdg_gage(2001:2008; snowy = yr -> yr <= 2004)  # span 8
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8)
        # 4 snowy years: existing min_rows=3 admits trends; gate must NOT fire
        @test !isnan(r["snow_on_dowy_senn_slp"])
    end

    @testset "collector unaffected by the gate" begin
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2000)
        c_on = AnnualCollector("rdgtest")
        c_off = AnnualCollector("rdgtest")
        calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                               min_values_for_stats=20, collector=c_on)
        calculate_snow_metrics(df; trend_completeness=0.6,
                               min_values_for_stats=20, collector=c_off)
        @test isequal(c_on.signature, c_off.signature)
        @test isequal(c_on.water_year, c_off.water_year)
        @test all(isequal.(c_on.value, c_off.value))
    end

    @testset "config flag is shipped true" begin
        @test StreamflowSignatures.CFG_SNOW_RECORD_DECADE_GATE === true
    end
```

Notes for the executor: (1) check `AnnualCollector`'s constructor signature in `julia/src/stats.jl` before using `AnnualCollector("rdgtest")` — match how `test_annual_collector.jl` builds one; (2) in the linked-threshold gage, the snowy-set expression lists 7 years of 2001–2010 — verify count = 7; (3) the `melt_com_dowy`/`melt_before_peak_pct` values are NaN in all-snow-free years by construction, which is exactly what the gate counts.

**Step 2: Run to verify failure** — Part 2 fails: gated metrics still carry computed trends (kwarg not yet wired in snow.jl).

**Step 3: Implement in `julia/src/snow.jl`** — two edits:

(a) After the `SNOW_METRICS` constant (line 46), add:

```julia
# Threshold-dependent metrics (NaN for operationally snow-free years) gated by the
# record-anchored decade gate (July 2026): a trend requires computable values in
# >= decade_completeness of the SWE-valid years in BOTH the first and last decade
# of the gage's SWE record. Magnitude metrics (valid zeros) are NOT gated — their
# dense series carry the snow-decline signal legitimately.
const SNOW_RECORD_GATED_METRICS = [
    "swe_max_dowy",
    "snow_on_dowy",
    "snow_off_dowy",
    "melt_season_days",
    "melt_rate",
    "ssm",
    "melt_before_peak",
    "melt_before_peak_pct",
    "melt_before_peak_to_max_swe",
    "melt_com_dowy"
]
```

(b) Replace the closing `generate_stats` call (lines 259-263) with:

```julia
    # Record-anchored decade gate (July 2026, user decision 2026-07-22): anchored to
    # the SWE-valid, grid-complete years (swe_max non-NaN — dense incl. zeros), NOT
    # the metric's own span, so clustered snowy years cannot slip through. Threshold
    # is the SAME decade_completeness the streamflow gate uses (linked knob).
    force_skip = nothing
    if CFG_SNOW_RECORD_DECADE_GATE && decade_completeness !== nothing
        anchor_mask = .!isnan.(annual_data.swe_max)
        anchor_years = Int.(annual_data.water_year[anchor_mask])
        if !isempty(anchor_years)
            rec_min, rec_max = extrema(anchor_years)
            if rec_max - rec_min + 1 >= 10
                windows = ((rec_min, rec_min + 9), (rec_max - 9, rec_max))
                skip = Set{String}()
                for m in SNOW_RECORD_GATED_METRICS
                    vals = annual_data[!, Symbol(m)]
                    for (lo, hi) in windows
                        den = 0
                        num = 0
                        for i in 1:nrow(annual_data)
                            anchor_mask[i] || continue
                            yr = Int(annual_data.water_year[i])
                            (lo <= yr <= hi) || continue
                            den += 1
                            isnan(vals[i]) || (num += 1)
                        end
                        if den > 0 && (num / den) < decade_completeness
                            push!(skip, m)
                            break
                        end
                    end
                end
                if !isempty(skip)
                    force_skip = skip
                end
            end
        end
    end

    stats = generate_stats(annual_data; value_cols=SNOW_METRICS,
                           trend_completeness=trend_completeness,
                           decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats,
                           changepoint=changepoint, collector=collector,
                           force_skip_trends=force_skip)
    merge!(result, stats)
    return result
```

(The 0-row/no-SWE defensive branch at lines 122-132 needs no gate — nothing to anchor; leave unchanged.)

**Step 4: Run the new test file** — expected: ALL PASS.

**Step 5: Run the adjacent suites** — `test_snow_metrics.jl` (existing snow tests pass `decade_completeness=0.8` in three places — the gate now activates there; verify none of those synthetic gages trip it: they are all-years-snowy or explicitly boundary designs. If any existing assertion breaks, STOP and reassess the gage design rather than weakening the assertion), `test_stats_floor.jl`, `test_annual_collector.jl`, `test_changepoint.jl`, then the full `runtests.jl`. Expected: all PASS.

---

### Task 4: Documentation

**Files:**
- Modify: `docs/SIGNATURES.md` (§12 Snow Metrics → Method + Data Quality & Caveats)
- Modify: `claude-skill/streamflow-signatures.md` (snow/completeness bullets)
- Modify: `CHANGELOG.md` (`[July 2026]` new entry + `[Unreleased]` open question)

**Step 1: SIGNATURES.md** — in §12 Method, add a step 6:

> 6. **Record-anchored decade gate** (July 2026): the 10 timing/melt/regime metrics
> additionally require ≥ `decade_min_fraction` (linked to the streamflow decade
> gate; 0.80 as shipped) of the SWE-valid years in BOTH the first and last decade
> of the gage's SWE record to be computable (snowy). Failing metrics have their 6
> trend statistics set to NaN; mean/median and Pettitt fields are still computed.
> Anchored to the SWE record — not the metric's own span — so gages whose snow
> disappears (or appears) mid-record do not carry trends conditioned on
> snow-present years. Magnitude metrics are exempt (their dense zero-including
> series carry the decline signal legitimately). Config:
> `snow.record_anchored_decade_gate`.

And in Data Quality & Caveats, extend the zero-vs-NA bullet with a sentence pointing at the gate.

**Step 2: claude-skill** — after the Trend completeness bullet (line ~118), add:

> - **Snow record-anchored decade gate**: the 10 threshold-dependent snow metrics also require ≥ decade_min_fraction of SWE-valid years in the record's first AND last decade to be snowy/computable, else trend stats are NaN (linked to the streamflow decade knob)

**Step 3: CHANGELOG** — new entry at the top of `[July 2026]` (title: "New: record-anchored decade gate for threshold-dependent snow metrics (Julia)") covering: motivation (clustered snowy years pass the own-span gate + floor), the four user decisions, linked-knob design, anchor/denominator semantics, failure semantics, config flag, tests, output impact (trend NaNs at next run; `flagged_for_high_na` shift), ports deferred with the snow port. PLUS under `[Unreleased]` → Guidelines Document TODOs:

> - [ ] **Decade gate 60% vs 80% — needs domain-expert decision.** User (2026-07-22)
>   twice described the first/last-decade gate as 60%; config, guidelines doc, and
>   manuscript §2.2.3 all say 80% (only the overall gate moved to 60%). The snow
>   record-anchored gate is LINKED to the same knob, so a decision for 60% is a
>   one-line change (`decade_min_fraction: 0.80 → 0.60`) + guidelines/manuscript
>   edits. Until decided, both gates run at 0.80.

**Step 4: config comments** — already handled in Task 1.

---

### Task 5: Full verification + Codex review

**Step 1:** Full Julia suite: run `julia/test/runtests.jl` end-to-end. Expected: all pass, zero new warnings.

**Step 2:** Config sanity in all consumers unaffected (flag lives under `snow`, which only Julia reads today): `python -c "from streamflow_signatures.config import NA_TREND_MIN_FRACTION"` still imports clean.

**Step 3:** Offer the user a Codex adversarial review of the diff (project convention for signature-methodology changes) before any benchmark rerun. The gate lands in the planned WY 1993–2025 standard run.

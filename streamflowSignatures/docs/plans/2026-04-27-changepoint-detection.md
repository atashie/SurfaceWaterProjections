# Changepoint Detection Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add 4-model BIC changepoint detection to all time-series streamflow signatures, producing 7 new per-gage columns per signature (~490 new columns).

**Architecture:** A new `changepoint.jl` module implements the core 4-model BIC comparison (flat/trend/step/piecewise-continuous). Integration into `generate_stats()` via a `changepoint` kwarg propagated through all signature functions. Changepoint analysis runs on annual values within WY 1980-2024, independent of the existing trend_completeness gate.

**Tech Stack:** Julia (DataFrames, Statistics, LinearAlgebra), existing `theil_sen_slope()` from stats.jl

---

## Task 1: Create `julia/src/changepoint.jl` — Core Module

**Files:**
- Create: `julia/src/changepoint.jl`

**Step 1: Write the changepoint detection module**

```julia
"""
Changepoint detection via 4-model BIC comparison.

Fits four competing models to an annual time series and selects the best via BIC:
1. Flat:      y = μ + ε                                    (k=2)
2. Trend:     y = α + βt + ε                               (k=3)
3. Step:      y = μ₁(t≤τ) + μ₂(t>τ) + ε                   (k=4)
4. Piecewise: y = α+β₁t (t≤τ), α+β₁τ+β₂(t-τ) (t>τ) + ε  (k=5, continuous)
"""

"""
    detect_changepoint(years, values; min_total_obs=20, min_segment_obs=10)

Run 4-model BIC comparison on annual time series.

# Arguments
- `years::AbstractVector{<:Real}`: Water years (may contain NaN)
- `values::AbstractVector{<:Real}`: Annual metric values (may contain NaN)
- `min_total_obs::Int`: Minimum non-NA observations required (default: 20)
- `min_segment_obs::Int`: Minimum non-NA observations per segment (default: 10)

# Returns
NamedTuple with fields:
- `best_model::Int`: 1=flat, 2=trend, 3=step, 4=piecewise
- `cp_year::Float64`: Changepoint water year (NaN if flat/trend wins)
- `delta_bic::Float64`: BIC(best_no_cp) - BIC(best_cp); positive = CP preferred
- `pre_slope::Float64`: Theil-Sen slope before CP (NaN if flat/trend)
- `post_slope::Float64`: Theil-Sen slope after CP (NaN if flat/trend)
- `pre_mean::Float64`: Mean before CP (NaN if flat/trend)
- `post_mean::Float64`: Mean after CP (NaN if flat/trend)
"""
function detect_changepoint(
    years::AbstractVector{<:Real},
    values::AbstractVector{<:Real};
    min_total_obs::Int = 20,
    min_segment_obs::Int = 10
)
    na_result = (
        best_model = NaN,
        cp_year = NaN,
        delta_bic = NaN,
        pre_slope = NaN,
        post_slope = NaN,
        pre_mean = NaN,
        post_mean = NaN
    )

    # Filter NaN pairs
    valid_mask = .!isnan.(years) .& .!isnan.(values)
    yrs = Float64.(years[valid_mask])
    vals = Float64.(values[valid_mask])
    n = length(yrs)

    if n < min_total_obs
        return na_result
    end

    # Sort by year
    perm = sortperm(yrs)
    yrs = yrs[perm]
    vals = vals[perm]

    # --- Model 1: Flat (k=2) ---
    bic_flat = _bic_flat(vals, n)

    # --- Model 2: Trend (k=3) ---
    bic_trend = _bic_trend(yrs, vals, n)

    # --- Models 3 & 4: Search over candidate breakpoints ---
    # Candidate τ: each unique year where both segments have ≥ min_segment_obs
    unique_yrs = sort(unique(yrs))
    best_bic_step = Inf
    best_tau_step = NaN
    best_bic_piecewise = Inf
    best_tau_piecewise = NaN

    for τ in unique_yrs
        n_pre = count(y -> y <= τ, yrs)
        n_post = n - n_pre
        if n_pre < min_segment_obs || n_post < min_segment_obs
            continue
        end

        pre_idx = findall(y -> y <= τ, yrs)
        post_idx = findall(y -> y > τ, yrs)

        # Model 3: Step at τ (k=4)
        bic_s = _bic_step(vals, pre_idx, post_idx, n)
        if bic_s < best_bic_step
            best_bic_step = bic_s
            best_tau_step = τ
        end

        # Model 4: Piecewise continuous at τ (k=5)
        bic_p = _bic_piecewise(yrs, vals, τ, n)
        if bic_p < best_bic_piecewise
            best_bic_piecewise = bic_p
            best_tau_piecewise = τ
        end
    end

    # If no valid breakpoint found, only flat/trend compete
    no_cp_candidates = isinf(best_bic_step) && isinf(best_bic_piecewise)

    # --- Select best model ---
    best_no_cp = min(bic_flat, bic_trend)
    best_cp = min(best_bic_step, best_bic_piecewise)

    if no_cp_candidates
        # No valid breakpoint — flat or trend wins
        best_model = bic_flat <= bic_trend ? 1 : 2
        delta_bic = NaN  # No CP model to compare against
        return (
            best_model = Float64(best_model),
            cp_year = NaN,
            delta_bic = NaN,
            pre_slope = NaN,
            post_slope = NaN,
            pre_mean = NaN,
            post_mean = NaN
        )
    end

    # Compare all 4 models
    all_bics = [bic_flat, bic_trend, best_bic_step, best_bic_piecewise]
    best_model = argmin(all_bics)
    delta_bic = best_no_cp - best_cp  # positive = CP preferred

    if best_model <= 2
        # Flat or Trend won — no changepoint
        return (
            best_model = Float64(best_model),
            cp_year = NaN,
            delta_bic = delta_bic,
            pre_slope = NaN,
            post_slope = NaN,
            pre_mean = NaN,
            post_mean = NaN
        )
    end

    # A changepoint model won — compute segment statistics
    τ_best = best_model == 3 ? best_tau_step : best_tau_piecewise
    pre_mask = yrs .<= τ_best
    post_mask = yrs .> τ_best

    pre_yrs = yrs[pre_mask]
    pre_vals = vals[pre_mask]
    post_yrs = yrs[post_mask]
    post_vals = vals[post_mask]

    # Theil-Sen slopes on each segment
    pre_slope, _ = theil_sen_slope(pre_yrs, pre_vals)
    post_slope, _ = theil_sen_slope(post_yrs, post_vals)

    return (
        best_model = Float64(best_model),
        cp_year = τ_best,
        delta_bic = delta_bic,
        pre_slope = pre_slope,
        post_slope = post_slope,
        pre_mean = mean(pre_vals),
        post_mean = mean(post_vals)
    )
end


# --- Internal BIC helpers ---

"""Compute BIC for flat model: y = μ + ε. k=2 (μ, σ²)."""
function _bic_flat(vals::Vector{Float64}, n::Int)
    μ = mean(vals)
    rss = sum((vals .- μ).^2)
    return _bic(rss, n, 2)
end

"""Compute BIC for linear trend model: y = α + βt + ε. k=3 (α, β, σ²)."""
function _bic_trend(yrs::Vector{Float64}, vals::Vector{Float64}, n::Int)
    rss = _ols_rss(yrs, vals)
    return _bic(rss, n, 3)
end

"""Compute BIC for step model at τ: y = μ₁(t≤τ) + μ₂(t>τ) + ε. k=4 (μ₁, μ₂, τ, σ²)."""
function _bic_step(vals::Vector{Float64}, pre_idx::Vector{Int}, post_idx::Vector{Int}, n::Int)
    μ₁ = mean(vals[pre_idx])
    μ₂ = mean(vals[post_idx])
    rss = sum((vals[pre_idx] .- μ₁).^2) + sum((vals[post_idx] .- μ₂).^2)
    return _bic(rss, n, 4)
end

"""
Compute BIC for continuous piecewise model at τ. k=5 (α, β₁, β₂, τ, σ²).

Model: y = α + β₁t (t≤τ), y = α + β₁τ + β₂(t-τ) (t>τ)
Rewritten as OLS: y = α + β₁t + γ·max(0, t-τ) where γ = β₂-β₁
"""
function _bic_piecewise(yrs::Vector{Float64}, vals::Vector{Float64}, τ::Float64, n::Int)
    # Hinge variable: max(0, t - τ)
    hinge = max.(0.0, yrs .- τ)

    # OLS with design matrix [1, t, hinge]
    X = hcat(ones(n), yrs, hinge)
    # Solve normal equations: β = (X'X) \ (X'y)
    XtX = X' * X
    Xty = X' * vals

    # Check for singularity
    if abs(det(XtX)) < 1e-20
        return Inf
    end

    β = XtX \ Xty
    predicted = X * β
    rss = sum((vals .- predicted).^2)
    return _bic(rss, n, 5)
end

"""Compute BIC = n·log(RSS/n) + k·log(n). Guards against RSS ≤ 0."""
function _bic(rss::Float64, n::Int, k::Int)
    if rss <= 0.0 || n <= 0
        return Inf
    end
    return n * log(rss / n) + k * log(n)
end

"""Compute OLS residual sum of squares for simple linear regression y = α + βx."""
function _ols_rss(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    x_mean = mean(x)
    y_mean = mean(y)
    Sxy = sum((x .- x_mean) .* (y .- y_mean))
    Sxx = sum((x .- x_mean).^2)
    if abs(Sxx) < 1e-20
        # Constant x — RSS = total SS
        return sum((y .- y_mean).^2)
    end
    β = Sxy / Sxx
    α = y_mean - β * x_mean
    predicted = α .+ β .* x
    return sum((y .- predicted).^2)
end
```

**Step 2: Verify file was created**

Run: `julia -e 'include("julia/src/changepoint.jl"); println("OK")'` from project root.
Expected: compilation errors (theil_sen_slope not available standalone) — that's fine, confirms syntax is valid when we check after integration.

---

## Task 2: Add Config Entries

**Files:**
- Modify: `config/signatures_config.json`
- Modify: `julia/src/config.jl`

**Step 1: Add changepoint section to signatures_config.json**

Add after the existing `"qa_qc"` section (before the closing `}`):

```json
  "changepoint": {
    "enabled": true,
    "start_water_year": 1980,
    "end_water_year": 2024,
    "min_total_obs": 20,
    "min_segment_obs": 10
  }
```

**Step 2: Add config constants to julia/src/config.jl**

Add after the existing NA handling constants (after the `CFG_NA_SEASON_MONTHS` block):

```julia
# Changepoint analysis configuration
const _cp_config = get(_config, "changepoint", Dict())
const CFG_CHANGEPOINT_ENABLED = Bool(get(_cp_config, "enabled", false))
const CFG_CP_START_WATER_YEAR = Int(get(_cp_config, "start_water_year", 1980))
const CFG_CP_END_WATER_YEAR = Int(get(_cp_config, "end_water_year", 2024))
const CFG_CP_MIN_TOTAL_OBS = Int(get(_cp_config, "min_total_obs", 20))
const CFG_CP_MIN_SEGMENT_OBS = Int(get(_cp_config, "min_segment_obs", 10))
```

---

## Task 3: Integrate into `generate_stats()` in `julia/src/stats.jl`

**Files:**
- Modify: `julia/src/stats.jl` (lines 243-396)

**Step 1: Add changepoint kwarg to generate_stats signature**

Add `changepoint` parameter to the function signature (after `decade_completeness`):

```julia
function generate_stats(
    df::DataFrame;
    value_cols::Union{Nothing, Vector{String}, Vector{Symbol}} = nothing,
    year_col::Union{String, Symbol} = "water_year",
    min_rows::Int = 3,
    trend_completeness::Union{Nothing, Float64} = nothing,
    decade_completeness::Union{Nothing, Float64} = nothing,
    changepoint::Union{Nothing, NamedTuple} = nothing
)
```

**Step 2: Add changepoint detection block inside the per-column loop**

Insert AFTER the mean/median computation (after line 392 `result["$(col)_median"] = median(valid_values)`) and BEFORE the `end` of the for loop:

```julia
        # --- Changepoint analysis (opt-in) ---
        if changepoint !== nothing
            cp_start = get(changepoint, :start_year, 1980)
            cp_end = get(changepoint, :end_year, 2024)
            cp_min_obs = get(changepoint, :min_total_obs, 20)
            cp_min_seg = get(changepoint, :min_segment_obs, 10)

            # Filter to analysis window
            cp_mask = valid_years .>= cp_start .&& valid_years .<= cp_end
            cp_years = valid_years[cp_mask]
            cp_vals = valid_values[cp_mask]

            cp_result = detect_changepoint(cp_years, cp_vals;
                min_total_obs=cp_min_obs, min_segment_obs=cp_min_seg)

            result["$(col)_best_model"] = cp_result.best_model
            result["$(col)_cp_year"] = cp_result.cp_year
            result["$(col)_cp_delta_bic"] = cp_result.delta_bic
            result["$(col)_pre_slope"] = cp_result.pre_slope
            result["$(col)_post_slope"] = cp_result.post_slope
            result["$(col)_pre_mean"] = cp_result.pre_mean
            result["$(col)_post_mean"] = cp_result.post_mean
        end
```

**Step 3: Handle the min_rows early-return path**

In the `nrow(df) < min_rows` block (lines 253-268), add NaN changepoint columns when changepoint analysis is enabled:

After the existing `result["$(col)$(suffix)"] = NaN` loop, add:
```julia
            if changepoint !== nothing
                for cp_suffix in ["_best_model", "_cp_year", "_cp_delta_bic",
                                  "_pre_slope", "_post_slope", "_pre_mean", "_post_mean"]
                    result["$(col)$(cp_suffix)"] = NaN
                end
            end
```

**Step 4: Handle the per-column min_rows check**

In the `length(valid_values) < min_rows` block (lines 301-306), add the same NaN changepoint fill:

After `result["$(col)$(suffix)"] = NaN`, add:
```julia
            if changepoint !== nothing
                for cp_suffix in ["_best_model", "_cp_year", "_cp_delta_bic",
                                  "_pre_slope", "_post_slope", "_pre_mean", "_post_mean"]
                    result["$(col)$(cp_suffix)"] = NaN
                end
            end
```

---

## Task 4: Propagate `changepoint` kwarg through all signature functions

**Files to modify** (13 files — same mechanical change in each):

For each file, add `changepoint::Union{Nothing, NamedTuple} = nothing` to the function signature, and pass `changepoint=changepoint` to the `generate_stats()` call.

### Pattern for each file:

**Before:**
```julia
function some_signature(data; trend_completeness=nothing, decade_completeness=nothing)
    # ... build annual DataFrame ...
    return generate_stats(annual_df; trend_completeness=trend_completeness, decade_completeness=decade_completeness)
end
```

**After:**
```julia
function some_signature(data; trend_completeness=nothing, decade_completeness=nothing, changepoint=nothing)
    # ... build annual DataFrame ...
    return generate_stats(annual_df; trend_completeness=trend_completeness, decade_completeness=decade_completeness, changepoint=changepoint)
end
```

### Files and functions:

1. `julia/src/flow_volumes.jl` — `calculate_flow_vols_by_year()`: add kwarg, pass to generate_stats
2. `julia/src/flashiness.jl` — `analyze_flashiness_trends()`: add kwarg, pass to generate_stats
3. `julia/src/timing.jl` — `analyze_flow_timing_trends()`: add kwarg, pass to generate_stats
4. `julia/src/fdc.jl` — `analyze_fdc_trends()` (may call internal helper): add kwarg, pass through
5. `julia/src/baseflow.jl` — `analyze_baseflow_indices()`: add kwarg, pass to generate_stats
6. `julia/src/baseflow.jl` — `analyze_baseflow_indices_with_parameters()`: add kwarg, pass to generate_stats
7. `julia/src/recession.jl` — `analyze_recession_parameters()`: add kwarg, pass to generate_stats. **Note:** recession currently has NO trend_completeness — add changepoint independently
8. `julia/src/pulses.jl` — `calculate_pulse_metrics()`: add kwarg, pass to generate_stats
9. `julia/src/pulses.jl` — `calculate_negative_days()`: add kwarg, pass to generate_stats
10. `julia/src/runoff_ratios.jl` — `analyze_Q_PPT_relationships()`: add kwarg, pass to generate_stats
11. `julia/src/elasticity.jl` — `calculate_streamflow_elasticity()`: add kwarg, pass to generate_stats. **Note:** elasticity has NO trend_completeness — add changepoint independently
12. `julia/src/qp_seasonality.jl` — `calculate_qp_seasonality()`: add kwarg, pass to generate_stats
13. `julia/src/storage.jl` — `calculate_average_storage()`: add kwarg, pass to generate_stats

**Important:** Read each file first to find the exact function signature and generate_stats() call location before editing. Some functions have multiple generate_stats() calls (e.g., elasticity calls it for rolling and annual separately).

---

## Task 5: Wire through `calculate_all_signatures()` and module exports

**Files:**
- Modify: `julia/src/signatures.jl`
- Modify: `julia/src/StreamflowSignatures.jl`

**Step 1: Add changepoint kwarg to calculate_all_signatures**

In `julia/src/signatures.jl`, add to function signature (line 28, after `include_qa_flags`):

```julia
    include_qa_flags::Bool=false,
    changepoint::Union{Nothing, NamedTuple}=nothing
```

**Step 2: Pass changepoint to all signature function calls**

Add `changepoint=changepoint` to every `merge!(results, ...)` call in signatures.jl. For example:

```julia
# Line 51 — flow volumes
merge!(results, calculate_flow_vols_by_year(gage_data; seasonal_flags=seasonal_flags, trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 57 — flashiness
merge!(results, analyze_flashiness_trends(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 63 — timing
merge!(results, analyze_flow_timing_trends(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 69 — FDC
merge!(results, analyze_fdc_trends(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 75 — baseflow fixed
merge!(results, analyze_baseflow_indices(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 83 — recession (no trend_completeness, but YES changepoint)
recession_results = analyze_recession_parameters(gage_data; changepoint=changepoint)

# Line 94 — baseflow parameterized
merge!(results, analyze_baseflow_indices_with_parameters(gage_data, recession_alpha; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 101 — pulses
merge!(results, calculate_pulse_metrics(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 107 — negative days
merge!(results, calculate_negative_days(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 118 — Q-PPT relationships
merge!(results, analyze_Q_PPT_relationships(cdata; seasonal_flags=seasonal_flags, trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 125 — elasticity (no trend_completeness, but YES changepoint)
merge!(results, calculate_streamflow_elasticity(cdata; changepoint=changepoint))

# Line 131 — Q-P seasonality
merge!(results, calculate_qp_seasonality(cdata; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))

# Line 137 — storage
merge!(results, calculate_average_storage(cdata; trend_completeness=tc, decade_completeness=dc, changepoint=changepoint))
```

**Step 3: Include changepoint.jl and add exports in StreamflowSignatures.jl**

In `julia/src/StreamflowSignatures.jl`:

Add include AFTER `stats.jl` (line 31) and BEFORE `io.jl`:
```julia
include("changepoint.jl")
```

Add export:
```julia
export detect_changepoint
```

Add config exports (with existing CFG exports):
```julia
export CFG_CHANGEPOINT_ENABLED, CFG_CP_START_WATER_YEAR, CFG_CP_END_WATER_YEAR,
       CFG_CP_MIN_TOTAL_OBS, CFG_CP_MIN_SEGMENT_OBS
```

---

## Task 6: Update Benchmark Runner

**Files:**
- Modify: `docs/benchmarks/run_julia_benchmark.jl`

**Step 1: Build changepoint config NamedTuple**

Add after the existing config loading (in the `main()` function, after the `use_legacy` assignment), roughly around line 50:

```julia
    # Changepoint analysis configuration
    cp_config = if CFG_CHANGEPOINT_ENABLED
        (
            start_year = CFG_CP_START_WATER_YEAR,
            end_year = CFG_CP_END_WATER_YEAR,
            min_total_obs = CFG_CP_MIN_TOTAL_OBS,
            min_segment_obs = CFG_CP_MIN_SEGMENT_OBS
        )
    else
        nothing
    end
```

**Step 2: Pass changepoint to calculate_all_signatures**

In the non-legacy path (around line 237-244), add `changepoint=cp_config`:

```julia
            signatures = calculate_all_signatures(
                gage_data, gage_has_climate;
                gage_id=gage_id_str,
                seasonal_flags=pp.seasonal_flags,
                trend_completeness=CFG_NA_TREND_MIN_FRACTION,
                decade_completeness=CFG_NA_DECADE_MIN_FRACTION,
                climate_data=climate_data,
                changepoint=cp_config
            )
```

In the legacy path (around line 222), also pass it:
```julia
            signatures = calculate_all_signatures(gage_data, has_climate; gage_id=gage_id_str, changepoint=cp_config)
```

**No other benchmark runner changes needed** — the new changepoint columns flow through the existing Dict → DataFrame → CSV pipeline automatically.

---

## Task 7: Test with Synthetic Data

**Files:**
- Create: `julia/test/test_changepoint.jl`

```julia
using Test
using Statistics

# Load the module
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using StreamflowSignatures

@testset "Changepoint Detection" begin

    @testset "Flat series — should select model 1" begin
        years = Float64.(1980:2024)
        values = fill(5.0, 45) .+ randn(45) .* 0.1
        result = detect_changepoint(years, values)
        @test result.best_model == 1.0  # Flat
        @test isnan(result.cp_year)
        @test isnan(result.pre_slope)
        @test isnan(result.post_slope)
    end

    @testset "Strong trend — should select model 2" begin
        years = Float64.(1980:2024)
        values = collect(1.0:45.0) .+ randn(45) .* 0.5
        result = detect_changepoint(years, values)
        @test result.best_model == 2.0  # Trend
        @test isnan(result.cp_year)
    end

    @testset "Clear step change — should select model 3" begin
        years = Float64.(1980:2024)
        values = vcat(fill(10.0, 22) .+ randn(22) .* 0.3,
                      fill(20.0, 23) .+ randn(23) .* 0.3)
        result = detect_changepoint(years, values)
        @test result.best_model == 3.0  # Step
        @test !isnan(result.cp_year)
        @test 1999.0 <= result.cp_year <= 2003.0  # Should be near 2001
        @test result.delta_bic > 0  # CP model preferred
        @test result.pre_mean < result.post_mean
        @test abs(result.pre_mean - 10.0) < 1.0
        @test abs(result.post_mean - 20.0) < 1.0
    end

    @testset "Trend reversal — should select model 4" begin
        years = Float64.(1980:2024)
        n = length(years)
        mid = 22  # WY 2001
        values = vcat(
            collect(1.0:Float64(mid)),            # increasing
            collect(Float64(mid):-1.0:Float64(mid-n+mid+1))  # decreasing
        )
        # Simpler: increasing then decreasing
        values = Float64[i <= mid ? Float64(i) : Float64(2*mid - i) for i in 1:n]
        values .+= randn(n) .* 0.3
        result = detect_changepoint(years, values)
        @test result.best_model == 4.0  # Piecewise
        @test !isnan(result.cp_year)
        @test result.pre_slope > 0  # Increasing before
        @test result.post_slope < 0  # Decreasing after
    end

    @testset "Too few observations — returns NaN" begin
        years = Float64.(1980:1994)  # 15 years, below min_total_obs=20
        values = randn(15)
        result = detect_changepoint(years, values)
        @test isnan(result.best_model)
        @test isnan(result.cp_year)
    end

    @testset "NaN handling — filters NaN pairs" begin
        years = Float64.(1980:2024)
        values = fill(5.0, 45) .+ randn(45) .* 0.1
        values[5:10] .= NaN  # 6 NaN values, 39 remaining (> 20)
        result = detect_changepoint(years, values)
        @test !isnan(result.best_model)  # Should still work
    end

    @testset "Delta BIC sign convention" begin
        # Flat series: delta_bic should be negative or NaN (no CP preferred)
        years = Float64.(1980:2024)
        values = fill(5.0, 45) .+ randn(45) .* 0.01
        result = detect_changepoint(years, values)
        if !isnan(result.delta_bic)
            @test result.delta_bic <= 0  # No CP preferred
        end

        # Step series: delta_bic should be positive
        values_step = vcat(fill(0.0, 22), fill(10.0, 23))
        result_step = detect_changepoint(years, values_step)
        @test result_step.delta_bic > 0  # CP preferred
    end

    @testset "Segment means for step model" begin
        years = Float64.(1980:2024)
        values = vcat(fill(5.0, 22), fill(15.0, 23))
        result = detect_changepoint(years, values)
        @test result.best_model == 3.0
        @test abs(result.pre_mean - 5.0) < 0.01
        @test abs(result.post_mean - 15.0) < 0.01
    end

    @testset "Min segment obs respected" begin
        years = Float64.(1980:2024)
        # Changepoint at year 5 — pre-segment too small with min_segment_obs=10
        values = vcat(fill(0.0, 5), fill(10.0, 40))
        result = detect_changepoint(years, values; min_segment_obs=10)
        # CP should NOT be at year 5 (only 5 pre-obs)
        if result.best_model >= 3.0
            pre_count = count(y -> y <= result.cp_year, years)
            @test pre_count >= 10
        end
    end
end
```

**Step 2: Run tests**

Run: `cd julia && julia --project=. test/test_changepoint.jl`
Expected: All tests pass.

---

## Task 8: Run Full Benchmark

**Step 1: Run the Julia benchmark**

```bash
cd /c/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures
julia docs/benchmarks/run_julia_benchmark.jl
```

Expected: ~15-20 min (existing ~15 min + ~5 min for changepoint analysis). Output CSV should have ~1,100+ columns (624 existing + ~490 changepoint columns).

**Step 2: Verify output**

Quick sanity check in Julia:
```julia
using CSV, DataFrames
df = CSV.read("docs/benchmarks/julia_signatures.csv", DataFrame)
println("Columns: $(ncol(df))")
println("Rows: $(nrow(df))")

# Check changepoint columns exist
cp_cols = [c for c in names(df) if contains(c, "_best_model")]
println("Changepoint signature count: $(length(cp_cols))")

# Check best_model distribution
for col in cp_cols[1:5]
    vals = df[!, col] |> skipmissing |> collect
    println("$col: flat=$(count(==(1.0), vals)) trend=$(count(==(2.0), vals)) step=$(count(==(3.0), vals)) pw=$(count(==(4.0), vals)) NA=$(count(ismissing, df[!, col]))")
end

# Spot-check a well-known signature
qann_cp = df[!, "Qann_cp_year"] |> skipmissing |> collect
println("Qann changepoint years — min: $(minimum(qann_cp)), max: $(maximum(qann_cp)), median: $(median(qann_cp))")
```

---

## Task 9: Update Documentation

**Files:**
- Modify: `CHANGELOG.md`
- Modify: `docs/SIGNATURES.md`
- Modify: `docs/DEVELOPMENT.md`
- Modify: `config.R` (register new column name patterns)

**Step 1: Add CHANGELOG entry under [Unreleased]**

Add a new section:
```markdown
### Changepoint Detection — 4-Model BIC Comparison (April 2026)

New changepoint analysis applied to all time-series signatures. For each annual metric, fits 4 competing models via BIC:
1. **Flat**: constant mean (k=2)
2. **Trend**: linear trend (k=3) — the implicit assumption of existing MK/Theil-Sen statistics
3. **Step**: abrupt level shift at unknown breakpoint τ (k=4)
4. **Piecewise continuous**: two joined trend lines meeting at τ (k=5)

Best model selected by BIC. Changepoint search window: WY 1980-2024, minimum 20 non-NA observations total, minimum 10 per segment.

**7 new columns per signature:**
| Column | Description |
|--------|-------------|
| `{sig}_best_model` | Winning model: 1=flat, 2=trend, 3=step, 4=piecewise |
| `{sig}_cp_year` | Changepoint water year (NaN if flat/trend wins) |
| `{sig}_cp_delta_bic` | BIC evidence for changepoint: positive = CP preferred |
| `{sig}_pre_slope` | Theil-Sen slope before changepoint |
| `{sig}_post_slope` | Theil-Sen slope after changepoint |
| `{sig}_pre_mean` | Mean before changepoint |
| `{sig}_post_mean` | Mean after changepoint |

**ΔBIC interpretation** (Kass & Raftery 1995): |ΔBIC| > 2 positive evidence, > 6 strong, > 10 very strong. Sign: positive = changepoint preferred, negative = no-changepoint preferred.

**Files modified:**
- `julia/src/changepoint.jl`: New module with 4-model BIC comparison
- `julia/src/stats.jl`: `generate_stats()` gains `changepoint` kwarg
- `julia/src/signatures.jl`: Passes changepoint config to all signature functions
- `julia/src/StreamflowSignatures.jl`: Include + exports
- `julia/src/config.jl`: New `CFG_CP_*` constants
- `config/signatures_config.json`: New `changepoint` section
- All 13 signature function files: `changepoint` kwarg propagation
- `docs/benchmarks/run_julia_benchmark.jl`: Passes changepoint config

**References:**
- Beaulieu & Killick (2018). Distinguishing trends and shifts from memory in climate data. J. Climate.
- Kass & Raftery (1995). Bayes factors. JASA.
- Villarini et al. (2009). On the stationarity of annual flood peaks. WRR.
```

**Step 2: Add to docs/SIGNATURES.md**

Add a new section "Changepoint Analysis" after the Summary Table, documenting the 4 models, output columns, and interpretation guidelines.

**Step 3: Update column count in docs/DEVELOPMENT.md**

Update references to "624 columns" to reflect the new total.

---

## Implementation Notes

### Column count estimate

Signatures producing annual time series (get changepoint columns):
- Flow volumes: 22 metrics × 7 = 154
- Flashiness: 1 × 7 = 7
- Flow timing: 15 × 7 = 105
- FDC: 3 × 7 = 21
- Baseflow fixed: 2 × 7 = 14
- Baseflow param: 2 × 7 = 14
- Recession: 8 × 7 = 56 (log_a_pointcloud, log_a_events, b_pointcloud, b_events, concavity, n_recession_events, alpha_linear + seasonality amplitude metrics that go through generate_stats)
- Pulses: 15 × 7 = 105
- Negative days: 1 × 7 = 7
- Q-PPT relationships: 5 × 7 = 35
- Elasticity: 2 × 7 = 14 (elasticity_rolling, elasticity_annual)
- Q-P seasonality: 2 × 7 = 14
- Average storage: 1 × 7 = 7

**Total estimate: ~553 new columns** (varies based on exact metrics going through generate_stats)

Signatures that are per-gage scalars (NO changepoint columns):
- elasticity_static, elasticity_years_total, elasticity_years_low_ppt
- recession_alpha_point_cloud_linear_reservoir
- runoff_ratio_high_count
- season_excluded_years_* (4 columns)
- ice_affected_days_total
- log_a_seasonality_minimum_* (3), log_a_seasonality_amplitude_* (3) — these are per-gage scalars, not annual time series

### Recession seasonality special case

The 6 recession seasonality metrics (log_a_seasonality_amplitude_all/first_half/last_half, log_a_seasonality_minimum_all/first_half/last_half) are per-gage scalars computed from the whole record — they do NOT go through generate_stats(). They should NOT get changepoint columns. Verify this by confirming they are added directly to the results Dict in recession.jl, not via generate_stats().

### Performance expectation

The changepoint detection is O(n² × T) per metric per gage, where n ≤ 80 (annual values) and T ≤ 60 (candidate breakpoints). For 7,313 gages × ~70 metrics, total is ~500K calls to detect_changepoint(), each doing ~60 OLS fits on ≤80 points. Expected: 3-8 minutes additional runtime.

### Testing checklist

After full benchmark:
- [ ] Total column count is reasonable (~1,100-1,200)
- [ ] `_best_model` values are 1-4 (no out-of-range)
- [ ] `_cp_year` values are within WY 1990-2014 range (10-year buffer from edges)
- [ ] `_cp_year` is NaN when `_best_model` is 1 or 2
- [ ] `_delta_bic` is positive when `_best_model` is 3 or 4
- [ ] `_delta_bic` is negative when `_best_model` is 1 or 2
- [ ] `_pre_slope` and `_post_slope` are NaN when `_best_model` is 1 or 2
- [ ] Gages with <20 valid years in WY 1980-2024 have all-NaN changepoint columns
- [ ] Recession metrics have many NaN changepoint columns (expected — sparse data)
- [ ] `Qann_best_model` distribution is plausible (mix of all 4 types)

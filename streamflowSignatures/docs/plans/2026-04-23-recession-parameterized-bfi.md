# Recession-Parameterized Baseflow Signatures Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add 3 new signatures — a linear-reservoir recession constant (`alpha_linear`), and two parameterized baseflow indices (`BFI_Eckhardt_param`, `BFI_LyneHollick_param`) — totaling 25 new output columns.

**Architecture:** Compute the discrete recession constant `alpha = Q_{i+1}/Q_i` under the linear reservoir assumption (b=1) using existing recession event infrastructure. Per-year point-cloud medians produce the `alpha_linear` signature (8 stats). A whole-record median produces `recession_alpha_point_cloud_linear_reservoir` (per-gage scalar) that parameterizes both Eckhardt and Lyne-Hollick filters for the new BFI signatures (8 stats each). All new signatures are additive — no existing calculations are modified.

**Tech Stack:** Julia (canonical implementation in `julia/src/`). Python and rpkg ports are out of scope for this plan.

**Key Design Decisions (from Codex review):**

1. **Discrete-time estimation** (Codex finding #1): Compute `alpha_i = Q_{i+1}/Q_i` directly from recession pairs. This avoids the continuous/discrete mapping ambiguity that would arise from using `-dQ/dt` and then converting. The discrete alpha IS the Eckhardt filter parameter — no transformation needed.

2. **First day removal**: Consistent with existing point-cloud fitting, remove the first day of each recession event (storm peak influence) before computing alpha values.

3. **Validity constraints**: Require `Q_i > 0` and `0 < alpha_i < 1` for each pair. The `alpha_i < 1` constraint is guaranteed during recession (Q is decreasing), but enforce it explicitly as a guard.

4. **Median aggregation** (Codex finding #3): Use median for both per-year and whole-record aggregation. Robust to outlier years and regime mixing.

5. **Minimum data gates**: Same 25-event minimum as existing recession metrics. Per-year point cloud requires >10 valid alpha values (matching existing `year_pc_Q` threshold at `recession.jl:393` which uses strict `> 10`).

6. **Alpha collection is independent of power-law fit** (Codex review finding): Alpha values (`Q_{i+1}/Q_i`) depend only on raw Q pairs, not on whether the power-law regression succeeded. Alpha collection must happen outside the `if !isnan(log_a) && !isnan(b)` gate to avoid unnecessarily narrowing the data pool.

7. **NaN fallbacks on all early-return paths** (Codex BLOCKER): `analyze_recession_parameters()` has 3 early-return paths (missing columns, empty data, <25 events). All must populate `recession_alpha_point_cloud_linear_reservoir = NaN` and `alpha_linear` empty_stats to prevent missing keys in the output dict.

8. **BFImax fixed at 0.8** for parameterized Eckhardt. Only `a` is parameterized.

9. **L-H parameterization is heuristic** (Codex finding #4): The Lyne-Hollick `alpha` parameter has no physical derivation from recession analysis. Using the recession-derived alpha is an exploratory heuristic, not a physically calibrated parameter.

10. **Trend completeness**: `alpha_linear` is part of recession analysis and inherently sparse — exempt from trend completeness gate (same as existing recession metrics). Parameterized BFI signatures DO apply trend completeness (same as existing BFI).

---

## Task 1: Add `alpha_linear` and `recession_alpha_point_cloud_linear_reservoir` to recession analysis

**Files:**
- Modify: `julia/src/recession.jl` (lines 254-518, `analyze_recession_parameters()`)

**Context:** The existing `analyze_recession_parameters()` already iterates through recession events per year, collecting point-cloud Q/dQdt pairs. We add a parallel collection of `alpha_i = Q_{i+1}/Q_i` values, compute per-year medians for the `alpha_linear` signature, and a whole-record median for the scalar.

**Step 1: Add `alpha_linear` to `base_metrics` list, DataFrame pre-allocation, and ALL early-return paths (lines 258-288, 300-308, 447-458)**

In the `base_metrics` list at line 258, add `"alpha_linear"`:
```julia
base_metrics = ["log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity", "alpha_linear"]
```

**CRITICAL (Codex BLOCKER fix):** The function has 3 early-return paths that must populate the new keys. Add to EACH early-return block:

Early-return 1 (missing columns, lines 268-277) — add after the existing `for m in seasonality_metrics` block:
```julia
        result["recession_alpha_point_cloud_linear_reservoir"] = NaN
```

Early-return 2 (empty data, lines 279-288) — same addition:
```julia
        result["recession_alpha_point_cloud_linear_reservoir"] = NaN
```

Early-return 3 (min_events gate, lines 449-458) — add BEFORE the `return result`:
```julia
        # Scalar still computed from available data (not gated by min_events)
        # — handled in Step 5 below. But alpha_linear stats must be in base_metrics loop.
```
Note: The min_events early-return path at lines 449-458 iterates `base_metrics` to fill empty_stats. Since `"alpha_linear"` is now in `base_metrics`, it will automatically get NaN-filled. The scalar is computed BEFORE this gate (see Step 5).

In the `annual_data` DataFrame at line 300, add the column:
```julia
annual_data = DataFrame(
    water_year = years,
    log_a_events = fill(NaN, length(years)),
    b_events = fill(NaN, length(years)),
    log_a_pointcloud = fill(NaN, length(years)),
    b_pointcloud = fill(NaN, length(years)),
    concavity = fill(NaN, length(years)),
    n_recession_events = fill(NaN, length(years)),
    alpha_linear = fill(NaN, length(years))
)
```

**Step 2: Add whole-record alpha collector (after line 294)**

Add a collector vector alongside the existing `all_log_a`, `all_b`, etc.:
```julia
all_alpha_linear = Float64[]  # Whole-record alpha values for scalar
```

**Step 3: Compute alpha values within per-event loop — OUTSIDE the power-law fit gate**

**CRITICAL (Codex fix):** Alpha computation (`Q_{i+1}/Q_i`) depends only on raw Q values, NOT on whether the power-law regression succeeded. Place this OUTSIDE the `if !isnan(log_a) && !isnan(b)` block (which starts at line 351), directly inside the `for (start_idx, end_idx) in events` loop. This ensures alpha values are collected even when the power-law fit fails (e.g., too few points, singular data).

Insert BEFORE line 351 (before `if !isnan(log_a) && !isnan(b)`), right after `log_a, b = fit_recession_power_law(...)` at line 349:
```julia
            # Compute discrete recession constant alpha = Q_{i+1}/Q_i (b=1 linear reservoir)
            # INDEPENDENT of power-law fit success — alpha depends only on raw Q pairs
            # Remove first day (storm peak) consistent with point-cloud fitting
            if length(Q_event) > 2  # Need at least 3 days (2 after removing first)
                Q_alpha = Q_event[2:end]  # Remove first day
                for j in 1:(length(Q_alpha) - 1)
                    if Q_alpha[j] > 0 && !isnan(Q_alpha[j]) && !isnan(Q_alpha[j + 1])
                        a_i = Q_alpha[j + 1] / Q_alpha[j]
                        if 0 < a_i < 1  # Must be valid recession (decreasing Q)
                            push!(year_alpha_linear, a_i)
                            push!(all_alpha_linear, a_i)
                        end
                    end
                end
            end
```

Note: `year_alpha_linear` must be initialized as `Float64[]` at the start of the per-year loop (alongside `year_log_a`, `year_b`, etc. at ~line 335):
```julia
year_alpha_linear = Float64[]
```

**Step 4: Store per-year alpha_linear (after the events loop, ~line 392)**

After the existing `yr_idx` lookup (line 389), add:
```julia
# Per-year alpha_linear (point cloud with b=1)
if length(year_alpha_linear) > 10
    annual_data[yr_idx, :alpha_linear] = median(year_alpha_linear)
end
```

**Step 5: Add whole-record scalar to output BEFORE the min_events gate (~line 447)**

**CRITICAL (Codex BLOCKER fix):** The scalar must be computed BEFORE the min_events early-return at lines 449-458. The scalar depends only on having enough alpha values (>10), not on having >=25 successfully-fitted power-law events. This allows gages with, say, 15 recession events (below the 25-event gate for other recession metrics) to still produce a valid parameterized BFI if they have enough alpha pairs.

Insert BEFORE the `if length(all_log_a) < min_events` check (line 449):
```julia
    # Whole-record recession alpha (scalar — independent of min_events gate)
    if length(all_alpha_linear) > 10
        result["recession_alpha_point_cloud_linear_reservoir"] = median(all_alpha_linear)
    else
        result["recession_alpha_point_cloud_linear_reservoir"] = NaN
    end
```

**Step 6: Verify `alpha_linear` gets stats via existing loop**

The existing loop at lines 461-474 iterates over `base_metrics` and calls `generate_stats()` for each. Since we added `"alpha_linear"` to `base_metrics` and the column to `annual_data`, it will automatically get 8 statistics. No additional code needed.

**Step 7: Run Julia tests**

```bash
cd "C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"
julia -e 'using Pkg; Pkg.activate("julia"); Pkg.test()'
```

Expected: All existing tests pass. New `alpha_linear_*` keys appear in recession output.

**Step 8: Commit**

```bash
git add julia/src/recession.jl
git commit -m "feat: add alpha_linear recession constant (b=1 linear reservoir) and whole-record scalar"
```

---

## Task 2: Add `analyze_baseflow_indices_with_parameters()` to baseflow module

**Files:**
- Modify: `julia/src/baseflow.jl` (append after line 244)

**Context:** The existing `eckhardt_filter()` and `lyne_hollick_filter()` already accept custom parameters via keyword arguments. The new function follows the same pattern as `analyze_baseflow_indices()` but uses a caller-supplied `alpha` for both filters.

**Step 1: Add the new function (append to baseflow.jl after line 244)**

```julia
"""
    analyze_baseflow_indices_with_parameters(df::DataFrame, alpha::Float64; BFImax=0.8, passes=2) -> Dict

Calculate baseflow index signatures using recession-derived filter parameters.

Uses the recession-derived discrete alpha (from linear reservoir assumption, b=1)
to parameterize both Eckhardt and Lyne-Hollick filters. BFImax remains fixed at 0.8.
The Lyne-Hollick parameterization is heuristic — the L-H alpha parameter has no
physical derivation from recession analysis.

Produces 2 metrics × 8 statistics = 16 columns:
- BFI_Eckhardt_param_*
- BFI_LyneHollick_param_*

# Arguments
- `df::DataFrame`: Daily streamflow data with columns: water_year, Q, dowy
- `alpha::Float64`: Recession-derived filter parameter (from recession_alpha_point_cloud_linear_reservoir)
- `BFImax::Float64`: Maximum baseflow index for Eckhardt filter (default 0.8)
- `passes::Int`: Number of passes for Lyne-Hollick filter (default 2)
"""
function analyze_baseflow_indices_with_parameters(
    df::DataFrame,
    alpha::Float64;
    BFImax::Float64=CFG_ECKHARDT_BFIMAX,
    passes::Int=CFG_LYNE_HOLLICK_PASSES,
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing
)
    result = Dict{String, Float64}()
    metrics = ["BFI_Eckhardt_param", "BFI_LyneHollick_param"]

    # Validate alpha is in valid range
    if isnan(alpha) || alpha <= 0.0 || alpha >= 1.0
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year", "dowy"])
    if !valid
        @warn "analyze_baseflow_indices_with_parameters: Missing columns: $missing_cols"
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    if nrow(df) == 0
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    Q_clean = coalesce_q(df.Q)
    years = unique(df.water_year)
    annual_data = DataFrame()

    for yr in years
        year_mask = coalesce.(df.water_year .== yr, false)
        if !any(year_mask)
            continue
        end
        dowy_clean = [coalesce(d, 999) for d in df.dowy[year_mask]]
        year_Q = Q_clean[year_mask]

        sort_idx = sortperm(dowy_clean)
        Q = year_Q[sort_idx]

        # Apply filters with recession-derived alpha
        bf_eck = eckhardt_filter(Q; BFImax=BFImax, a=alpha)
        bf_lh = lyne_hollick_filter(Q; alpha=alpha, passes=passes)

        # Calculate BFI (same logic as analyze_baseflow_indices)
        valid_q_mask = .!isnan.(Q)
        total_Q = sum(Q[valid_q_mask])
        if total_Q <= 0
            continue
        end

        valid_eck_mask = valid_q_mask .& .!isnan.(bf_eck)
        valid_lh_mask = valid_q_mask .& .!isnan.(bf_lh)

        bfi_eck = sum(valid_eck_mask) > 0 ? sum(bf_eck[valid_eck_mask]) / total_Q : NaN
        bfi_lh = sum(valid_lh_mask) > 0 ? sum(bf_lh[valid_lh_mask]) / total_Q : NaN

        bfi_eck = clamp(bfi_eck, 0.0, 1.0)
        bfi_lh = clamp(bfi_lh, 0.0, 1.0)

        push!(annual_data, (water_year=yr, BFI_Eckhardt_param=bfi_eck, BFI_LyneHollick_param=bfi_lh))
    end

    result = generate_stats(annual_data; value_cols=metrics,
        trend_completeness=trend_completeness, decade_completeness=decade_completeness)

    return result
end
```

**Step 2: Run Julia tests**

```bash
julia -e 'using Pkg; Pkg.activate("julia"); Pkg.test()'
```

Expected: Tests pass. New function is callable.

**Step 3: Commit**

```bash
git add julia/src/baseflow.jl
git commit -m "feat: add analyze_baseflow_indices_with_parameters() for recession-parameterized BFI"
```

---

## Task 3: Wire new signatures into orchestration

**Files:**
- Modify: `julia/src/signatures.jl` (lines 80-85, after recession call)
- Modify: `julia/src/StreamflowSignatures.jl` (line 64, exports)

**Context:** The orchestration function `calculate_all_signatures()` calls each signature function sequentially. The parameterized BFI depends on the recession scalar, so it must run AFTER recession analysis. We extract the scalar from recession results and pass it to the new baseflow function.

**Step 1: Modify recession call to capture scalar and pass to parameterized BFI (signatures.jl, lines 80-85)**

Replace the existing recession try-catch block:
```julia
    # Recession: inherently sparse (event-based, many years with no events) — no trend completeness
    try
        merge!(results, analyze_recession_parameters(gage_data))
    catch e
        @warn "analyze_recession_parameters failed for gage $gage_id" exception=(e, catch_backtrace())
    end
```

With:
```julia
    # Recession: inherently sparse (event-based, many years with no events) — no trend completeness
    local recession_alpha = NaN
    try
        recession_results = analyze_recession_parameters(gage_data)
        # Extract scalar for parameterized BFI before merging
        recession_alpha = get(recession_results, "recession_alpha_point_cloud_linear_reservoir", NaN)
        merge!(results, recession_results)
    catch e
        @warn "analyze_recession_parameters failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    # Parameterized BFI using recession-derived alpha (requires recession to have run first)
    # Uses trend completeness (same as fixed-parameter BFI)
    try
        merge!(results, analyze_baseflow_indices_with_parameters(gage_data, recession_alpha;
            trend_completeness=tc, decade_completeness=dc))
    catch e
        @warn "analyze_baseflow_indices_with_parameters failed for gage $gage_id" exception=(e, catch_backtrace())
    end
```

**Step 2: Add export to StreamflowSignatures.jl (line 64)**

Change:
```julia
export analyze_baseflow_indices, eckhardt_filter, lyne_hollick_filter
```
To:
```julia
export analyze_baseflow_indices, analyze_baseflow_indices_with_parameters, eckhardt_filter, lyne_hollick_filter
```

**Step 3: Run Julia tests**

```bash
julia -e 'using Pkg; Pkg.activate("julia"); Pkg.test()'
```

Expected: All tests pass. `calculate_all_signatures()` now returns 25 new columns.

**Step 4: Commit**

```bash
git add julia/src/signatures.jl julia/src/StreamflowSignatures.jl
git commit -m "feat: wire recession-parameterized BFI into signature orchestration"
```

---

## Task 4: Register new signatures in config.R

**Files:**
- Modify: `config.R` (lines 571-574 and 611-615)

**Context:** `config.R` maintains the `EXPECTED_SIGNATURE_BASES` list for schema validation. New 8-stat signatures must be added to the list. The per-gage scalar must be added as a non-standard exception.

**Step 1: Add to EXPECTED_SIGNATURE_BASES (config.R, ~line 573)**

After the existing recession entries:
```r
  "log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity",
  "n_recession_events",
```

Add:
```r
  "log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity",
  "n_recession_events", "alpha_linear",
```

After the existing baseflow entries (line 571):
```r
  "BFI_Eckhardt", "BFI_LyneHollick",
```

Add:
```r
  "BFI_Eckhardt", "BFI_LyneHollick",
  "BFI_Eckhardt_param", "BFI_LyneHollick_param",
```

**Step 2: Add scalar exception (config.R, after line 615)**

After `EXPECTED_ELASTICITY_DIAGNOSTICS`, add:
```r
# Recession-derived alpha for parameterized BFI (per-gage scalar, not 8-stat pattern)
EXPECTED_RECESSION_ALPHA_SCALAR <- "recession_alpha_point_cloud_linear_reservoir"
```

And update the validation assembly (~line 644) to include it:
```r
  expected_sig_cols <- c(expected_sig_cols, EXPECTED_RECESSION_SEASONALITY,
                         EXPECTED_RUNOFF_RATIO_HIGH, EXPECTED_ELASTICITY_DIAGNOSTICS,
                         EXPECTED_SEASON_EXCLUDED_YEARS, EXPECTED_RECESSION_ALPHA_SCALAR)
```

**Step 3: Commit**

```bash
git add config.R
git commit -m "feat: register alpha_linear, BFI_Eckhardt_param, BFI_LyneHollick_param in config.R"
```

---

## Task 5: Run full Julia benchmark and validate

**Files:**
- None modified (validation only)

**Context:** The benchmark processes 7,313 gages (~27 min). This validates that the new signatures produce sensible values across all gages and don't break existing output.

**Step 1: Run benchmark**

```bash
cd "C:/Users/arikt/Documents/GitHub/SurfaceWaterProjections/streamflowSignatures"
julia docs/benchmarks/run_julia_benchmark.jl
```

Expected: ~27 min runtime, 619 signature columns (was 594).

**Step 2: Spot-check output**

After benchmark completes, verify in Julia REPL:
```julia
using CSV, DataFrames
df = CSV.read("path/to/output.csv", DataFrame)

# Check new columns exist
@assert "alpha_linear_mean" in names(df)
@assert "alpha_linear_median" in names(df)
@assert "BFI_Eckhardt_param_mean" in names(df)
@assert "BFI_LyneHollick_param_mean" in names(df)
@assert "recession_alpha_point_cloud_linear_reservoir" in names(df)

# Check value ranges
alpha_vals = df.recession_alpha_point_cloud_linear_reservoir
println("Alpha scalar: $(sum(.!ismissing.(alpha_vals) .& .!isnan.(alpha_vals))) non-NA values")
println("Alpha range: $(minimum(skipmissing(filter(x -> !isnan(x), collect(skipmissing(alpha_vals)))))) to $(maximum(skipmissing(filter(x -> !isnan(x), collect(skipmissing(alpha_vals))))))")
# Expected: values between ~0.90 and ~0.995, most near 0.95-0.99

bfi_eck = df.BFI_Eckhardt_param_mean
println("BFI_Eckhardt_param_mean range: $(minimum(skipmissing(filter(x -> !isnan(x), collect(skipmissing(bfi_eck)))))) to $(maximum(skipmissing(filter(x -> !isnan(x), collect(skipmissing(bfi_eck))))))")
# Expected: values between 0 and 1

# Compare with fixed-parameter BFI
println("Correlation with fixed-param Eckhardt: ", cor(collect(skipmissing(df.BFI_Eckhardt_mean)), collect(skipmissing(bfi_eck))))
```

**Step 3: Commit (if output file is tracked)**

No code commit needed here — this is validation only.

---

## Task 6: Update documentation

**Files:**
- Modify: `docs/SIGNATURES.md` (Section 2 Baseflow, Section 3 Recession)
- Modify: `CHANGELOG.md` (add entry under [Unreleased])

**Step 1: Add to SIGNATURES.md — Recession section (Section 3)**

Add `alpha_linear` to the recession metrics table:
```markdown
| **alpha_linear** | Discrete recession constant under linear reservoir (b=1); median of Q_{i+1}/Q_i ratios from point cloud |
```

Add scalar to the notes:
```markdown
- `recession_alpha_point_cloud_linear_reservoir` is a per-gage scalar (whole-record median of alpha ratios). Used to parameterize BFI_Eckhardt_param and BFI_LyneHollick_param.
```

**Step 2: Add to SIGNATURES.md — Baseflow section (Section 2)**

Add new metrics:
```markdown
| **BFI_Eckhardt_param** | Baseflow index using Eckhardt filter with recession-derived alpha | a=recession_alpha, BFImax=0.8 |
| **BFI_LyneHollick_param** | Baseflow index using Lyne-Hollick filter with recession-derived alpha (heuristic) | alpha=recession_alpha, 2 passes |
```

Add note:
```markdown
### Recession-Parameterized BFI

The parameterized BFI signatures use a recession-derived discrete filter constant instead of fixed defaults. The filter constant `alpha` is estimated as the median of `Q_{i+1}/Q_i` ratios across all recession events in the entire record (stored as `recession_alpha_point_cloud_linear_reservoir`).

- **Eckhardt**: `alpha` directly replaces the fixed `a=0.98`. BFImax remains fixed at 0.8.
- **Lyne-Hollick**: `alpha` replaces the fixed `alpha=0.925`. This is a heuristic — the L-H filter parameter has no physical derivation from recession analysis (Nathan & McMahon, 1990).
- Gages with <25 recession events produce NA for all parameterized BFI values (same gate as other recession metrics).

References:
- Eckhardt, K. (2005). How to construct recursive digital filters for baseflow separation.
- Collischonn, W., & Fan, F.M. (2013). Defining parameters for Eckhardt's digital baseflow filter.
```

**Step 3: Add CHANGELOG.md entry**

Under `[Unreleased]`, add:
```markdown
### Recession-Parameterized Baseflow Signatures (April 2026)

Three new signatures using recession-derived filter parameters (Collischonn & Fan 2013, Eckhardt 2005).

**Signature 1: `alpha_linear`** (8 stats via generate_stats)
- Discrete recession constant under linear reservoir assumption (b=1)
- Per-year computation: median of Q_{i+1}/Q_i ratios from point-cloud recession data
- First day of each recession event removed (storm peak influence)
- Minimum >10 valid alpha values per year (matching existing point-cloud threshold)
- Same 25-event gate as other recession metrics for alpha_linear trend stats
- Whole-record scalar computed independently of 25-event gate (only requires >10 alpha pairs)

**Scalar: `recession_alpha_point_cloud_linear_reservoir`** (per-gage diagnostic)
- Whole-record median of Q_{i+1}/Q_i ratios across all recession events
- Used to parameterize BFI_Eckhardt_param and BFI_LyneHollick_param

**Signature 2: `BFI_Eckhardt_param`** (8 stats)
- Eckhardt BFI with recession-derived `a = recession_alpha_point_cloud_linear_reservoir`
- BFImax fixed at 0.8 (not estimated from backward filtering)

**Signature 3: `BFI_LyneHollick_param`** (8 stats)
- Lyne-Hollick BFI with recession-derived `alpha = recession_alpha_point_cloud_linear_reservoir`
- Heuristic parameterization — L-H alpha has no physical derivation from recession (Nathan & McMahon 1990)

**Total: 24 new stat columns + 1 per-gage scalar = 25 new columns (619 total)**

Implementation: Julia only (canonical). Python and rpkg ports deferred.
```

**Step 4: Commit**

```bash
git add docs/SIGNATURES.md CHANGELOG.md
git commit -m "docs: add recession-parameterized BFI signatures to SIGNATURES.md and CHANGELOG"
```

---

## Summary

| Task | Files | Columns Added | Description |
|------|-------|---------------|-------------|
| 1 | `julia/src/recession.jl` | 8 (`alpha_linear_*`) + 1 scalar | Linear recession constant per-year + whole-record |
| 2 | `julia/src/baseflow.jl` | 16 (`BFI_*_param_*`) | Parameterized BFI function |
| 3 | `julia/src/signatures.jl`, `StreamflowSignatures.jl` | — | Wire into orchestration |
| 4 | `config.R` | — | Register for schema validation |
| 5 | — | — | Benchmark validation |
| 6 | `docs/SIGNATURES.md`, `CHANGELOG.md` | — | Documentation |

**Total: 25 new columns (619 total, was 594)**

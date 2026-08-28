"""
Changepoint detection.

`pettitt_test` (below) is the PRODUCTION changepoint method: it is what
`generate_stats()` applies to every signature's annual series to produce the
8 `_pettitt_*` output fields.

The 4-model BIC comparison (`detect_changepoint` and its helpers) is retained
by design but is NOT called by the pipeline. It was built to distinguish step
changes from trend breaks — something the Pettitt test cannot do — and is kept
for potential future use:
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
- `best_model::Float64`: 1=flat, 2=trend, 3=step, 4=piecewise (NaN if insufficient data)
- `cp_year::Float64`: Changepoint water year (NaN if flat/trend wins or insufficient data)
- `delta_bic::Float64`: BIC(best_no_cp) - BIC(best_cp); positive = CP preferred
- `pre_slope::Float64`: Theil-Sen slope before CP (NaN if flat/trend wins)
- `post_slope::Float64`: Theil-Sen slope after CP (NaN if flat/trend wins)
- `pre_mean::Float64`: Mean before CP (NaN if flat/trend wins)
- `post_mean::Float64`: Mean after CP (NaN if flat/trend wins)
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

    best_no_cp = min(bic_flat, bic_trend)
    best_cp = min(best_bic_step, best_bic_piecewise)

    if no_cp_candidates
        best_model = bic_flat <= bic_trend ? 1 : 2
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
    delta_bic = best_no_cp - best_cp

    if best_model <= 2
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
    hinge = max.(0.0, yrs .- τ)
    X = hcat(ones(n), yrs, hinge)
    XtX = X' * X
    Xty = X' * vals

    if abs(det(XtX)) < 1e-20
        return Inf
    end

    β = XtX \ Xty
    predicted = X * β
    rss = sum((vals .- predicted).^2)
    return _bic(rss, n, 5)
end

"""
Compute BIC = n·log(RSS/n) + k·log(n).

RSS=0 (perfect fit) yields -Inf (best possible BIC).
RSS<0 should not occur; returns Inf as a guard.
"""
function _bic(rss::Float64, n::Int, k::Int)
    if n <= 0
        return Inf
    end
    if rss < 0.0
        return Inf  # Should not happen; defensive guard
    end
    if rss == 0.0
        return -Inf  # Perfect fit
    end
    return n * log(rss / n) + k * log(n)
end

"""
    pettitt_test(years, values; min_total_obs=20, min_segment_obs=10)

Run the Pettitt non-parametric changepoint test on an annual time series.

The test statistic K_T = max|U_t| where U_t = Σ_{i≤t, j>t} sgn(x_i - x_j).
P-value uses the asymptotic approximation: p ≈ 2·exp(-6K²/(T³+T²)).

# Returns
NamedTuple with fields:
- `cp_year::Float64`: Most likely changepoint year (NaN if insufficient data)
- `pval::Float64`: Asymptotic p-value (NaN if insufficient data)
- `pre_mean::Float64`: Mean before changepoint
- `post_mean::Float64`: Mean after changepoint
"""
function pettitt_test(
    years::AbstractVector{<:Real},
    values::AbstractVector{<:Real};
    min_total_obs::Int = 20,
    min_segment_obs::Int = 10
)
    na_result = (cp_year = NaN, pval = NaN, pre_mean = NaN, post_mean = NaN)

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

    # Compute U_t via recursive relation: U_t = U_{t-1} + V_t
    # where V_t = Σ_{j=1}^{n} sgn(x_t - x_j)
    U = zeros(n)
    for t in 1:n
        V_t = 0.0
        for j in 1:n
            if vals[t] > vals[j]
                V_t += 1.0
            elseif vals[t] < vals[j]
                V_t -= 1.0
            end
        end
        U[t] = (t == 1 ? 0.0 : U[t-1]) + V_t
    end

    # Find τ* = argmax |U_t|, respecting min_segment_obs
    best_abs_U = -1.0
    best_t = 0
    for t in min_segment_obs:(n - min_segment_obs)
        abs_U = abs(U[t])
        if abs_U > best_abs_U
            best_abs_U = abs_U
            best_t = t
        end
    end

    if best_t == 0
        return na_result
    end

    K = best_abs_U
    # Asymptotic p-value (Pettitt 1979)
    pval = 2.0 * exp(-6.0 * K^2 / (Float64(n)^3 + Float64(n)^2))
    pval = min(pval, 1.0)

    pre_vals = vals[1:best_t]
    post_vals = vals[(best_t+1):n]

    return (
        cp_year = yrs[best_t],
        pval = pval,
        pre_mean = mean(pre_vals),
        post_mean = mean(post_vals)
    )
end

"""
    segment_differential_metrics(years, values, cp_year)

Compute differential metrics for a changepoint: delta_mean, pct_change,
pre/post Mann-Kendall p-values. Returns NaN tuple if cp_year is NaN.
"""
function segment_differential_metrics(
    years::Vector{Float64},
    values::Vector{Float64},
    cp_year::Float64
)
    na_result = (delta_mean = NaN, pct_change = NaN, pre_mk_pval = NaN, post_mk_pval = NaN)

    if isnan(cp_year) || length(years) == 0
        return na_result
    end

    pre_mask = years .<= cp_year
    post_mask = years .> cp_year

    pre_yrs = years[pre_mask]
    pre_vals = values[pre_mask]
    post_yrs = years[post_mask]
    post_vals = values[post_mask]

    if length(pre_vals) < 3 || length(post_vals) < 3
        return na_result
    end

    pre_m = mean(pre_vals)
    post_m = mean(post_vals)
    delta_mean = post_m - pre_m
    pct_change = abs(pre_m) > 1e-10 ? delta_mean / abs(pre_m) * 100.0 : NaN

    # Mann-Kendall p-values for each segment (MK test takes values only, sorted by time)
    _, pre_mk_pval = mann_kendall_test(pre_vals)
    _, post_mk_pval = mann_kendall_test(post_vals)

    return (delta_mean = delta_mean, pct_change = pct_change,
            pre_mk_pval = pre_mk_pval, post_mk_pval = post_mk_pval)
end

"""Compute OLS residual sum of squares for simple linear regression y = α + βx."""
function _ols_rss(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    x_mean = mean(x)
    y_mean = mean(y)
    Sxy = sum((x .- x_mean) .* (y .- y_mean))
    Sxx = sum((x .- x_mean).^2)
    if abs(Sxx) < 1e-20
        return sum((y .- y_mean).^2)
    end
    β = Sxy / Sxx
    α = y_mean - β * x_mean
    predicted = α .+ β .* x
    return sum((y .- predicted).^2)
end

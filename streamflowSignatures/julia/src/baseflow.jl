"""
Baseflow signatures.

Digital filter methods for baseflow separation and BFI calculation.
"""

# Configuration loaded from config.jl:
# CFG_ECKHARDT_BFIMAX, CFG_ECKHARDT_ALPHA, CFG_LYNE_HOLLICK_ALPHA,
# CFG_LYNE_HOLLICK_PASSES

"""
    eckhardt_filter(Q::Vector; BFImax=0.8, a=0.98) -> Vector

Apply Eckhardt recursive digital filter for baseflow separation.

Parameters
----------
Q : Vector
    Daily discharge values
BFImax : Float64
    Maximum baseflow index (default 0.8)
a : Float64
    Filter parameter (default 0.98)

Returns
-------
Vector
    Baseflow component
"""
function eckhardt_filter(Q::AbstractVector{<:Real}; BFImax::Float64=CFG_ECKHARDT_BFIMAX, a::Float64=CFG_ECKHARDT_ALPHA)
    n = length(Q)
    if n == 0
        return Float64[]
    end

    baseflow = zeros(n)

    # Initialize
    if !isnan(Q[1]) && Q[1] > 0
        baseflow[1] = min(BFImax * Q[1], Q[1])
    end

    # Recursive filter - PERFORMANCE FIX: Added @inbounds
    @inbounds for i in 2:n
        if isnan(Q[i])
            baseflow[i] = baseflow[i-1]
            continue
        end

        # Eckhardt equation
        numerator = (1 - BFImax) * a * baseflow[i-1] + (1 - a) * BFImax * Q[i]
        denominator = 1 - a * BFImax

        bf = numerator / denominator

        # Ensure baseflow <= Q and >= 0
        baseflow[i] = max(0.0, min(bf, Q[i]))
    end

    return baseflow
end


"""
    lyne_hollick_filter(Q::Vector; alpha=0.925, passes=2) -> Vector

Apply Lyne-Hollick digital filter for baseflow separation.

Parameters
----------
Q : Vector
    Daily discharge values
alpha : Float64
    Filter parameter (default 0.925)
passes : Int
    Number of forward-backward passes (default 2)

Returns
-------
Vector
    Baseflow component
"""
function lyne_hollick_filter(Q::AbstractVector{<:Real}; alpha::Float64=CFG_LYNE_HOLLICK_ALPHA, passes::Int=CFG_LYNE_HOLLICK_PASSES)
    n = length(Q)
    if n == 0
        return Float64[]
    end

    # Start with the input signal (will become baseflow after passes)
    input_signal = Float64.(Q)

    for pass_num in 1:passes
        # Forward pass: compute quickflow going forward
        qf_forward = zeros(n)
        qf_forward[1] = 0.0

        @inbounds for i in 2:n
            if isnan(input_signal[i]) || isnan(input_signal[i-1]) || isnan(qf_forward[i-1])
                qf_forward[i] = NaN
            else
                qf_forward[i] = alpha * qf_forward[i-1] +
                                ((1 + alpha) / 2) * (input_signal[i] - input_signal[i-1])
                # Constrain quickflow to [0, input_signal]
                qf_forward[i] = max(0.0, min(qf_forward[i], input_signal[i]))
            end
        end

        # Backward pass: refine quickflow going backward
        qf_backward = zeros(n)
        qf_backward[n] = qf_forward[n]

        @inbounds for i in (n-1):-1:1
            if isnan(qf_forward[i]) || isnan(qf_forward[i+1]) || isnan(qf_backward[i+1])
                qf_backward[i] = NaN
            else
                qf_backward[i] = alpha * qf_backward[i+1] +
                                 ((1 + alpha) / 2) * (qf_forward[i] - qf_forward[i+1])
                # Constrain quickflow to [0, input_signal]
                qf_backward[i] = max(0.0, min(qf_backward[i], input_signal[i]))
            end
        end

        # Baseflow from this pass becomes input for next pass
        @inbounds for i in 1:n
            input_signal[i] = input_signal[i] - qf_backward[i]
            if input_signal[i] < 0
                input_signal[i] = 0.0
            end
            if isnan(Q[i])
                input_signal[i] = NaN
            end
        end
    end

    # After all passes, input_signal holds the final baseflow
    return input_signal
end


"""
    analyze_baseflow_indices(df::DataFrame) -> Dict

Calculate baseflow index signatures with trend statistics.

Calculates 2 metrics:
- BFI_Eckhardt: Baseflow index using Eckhardt filter
- BFI_LyneHollick: Baseflow index using Lyne-Hollick filter

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, dowy
    (preprocessor guarantees all years are valid)

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (2 metrics × 8 stats = 16 values)
"""
function analyze_baseflow_indices(
    df::DataFrame;
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing,
    changepoint::Union{Nothing, NamedTuple}=nothing,
    collector::Union{Nothing, AnnualCollector}=nothing
)
    result = Dict{String, Float64}()
    metrics = ["BFI_Eckhardt", "BFI_LyneHollick"]

    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year", "dowy"])
    if !valid
        @warn "analyze_baseflow_indices: Missing columns: $missing_cols"
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

    # Convert Q to Float64, handling missing values
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

        # Sort by day of water year
        sort_idx = sortperm(dowy_clean)
        Q = year_Q[sort_idx]

        # Calculate baseflow using both filters
        bf_eck = eckhardt_filter(Q)
        bf_lh = lyne_hollick_filter(Q)

        # Calculate BFI (sum of baseflow / sum of total flow)
        # CRITICAL: Must filter by BOTH Q and baseflow being non-NaN
        # The Lyne-Hollick filter can propagate NaN from adjacent missing values
        # even to positions where Q itself is valid. R uses na.rm=TRUE which
        # automatically handles this, but Julia's sum() returns NaN if any input is NaN.

        # For total Q: only need Q to be valid (matching R/Python — zeros allowed)
        valid_q_mask = .!isnan.(Q)

        total_Q = sum(Q[valid_q_mask])
        if total_Q <= 0
            continue
        end

        # For baseflow sums: filter by BOTH Q valid AND baseflow valid
        # This matches R's na.rm=TRUE behavior for sum()
        valid_eck_mask = valid_q_mask .& .!isnan.(bf_eck)
        valid_lh_mask = valid_q_mask .& .!isnan.(bf_lh)

        # Calculate BFI using only valid values (matching R's na.rm=TRUE)
        bfi_eck = sum(valid_eck_mask) > 0 ? sum(bf_eck[valid_eck_mask]) / total_Q : NaN
        bfi_lh = sum(valid_lh_mask) > 0 ? sum(bf_lh[valid_lh_mask]) / total_Q : NaN

        # Ensure BFI in [0, 1]
        bfi_eck = clamp(bfi_eck, 0.0, 1.0)
        bfi_lh = clamp(bfi_lh, 0.0, 1.0)

        push!(annual_data, (water_year=yr, BFI_Eckhardt=bfi_eck, BFI_LyneHollick=bfi_lh))
    end

    # No early-return gate here — let generate_stats() handle min_rows internally
    # (matching R/Python, which pass whatever data exists to generate_stats)
    result = generate_stats(annual_data; value_cols=metrics, trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)

    return result
end


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
    decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing,
    changepoint::Union{Nothing, NamedTuple}=nothing,
    collector::Union{Nothing, AnnualCollector}=nothing
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
        trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)

    return result
end
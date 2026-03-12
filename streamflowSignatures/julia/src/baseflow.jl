"""
Baseflow signatures.

Digital filter methods for baseflow separation and BFI calculation.
"""

# Configuration loaded from config.jl:
# CFG_ECKHARDT_BFIMAX, CFG_ECKHARDT_ALPHA, CFG_LYNE_HOLLICK_ALPHA,
# CFG_LYNE_HOLLICK_PASSES, CFG_BASEFLOW_MIN_DAYS, CFG_BASEFLOW_MAX_MISSING_FRAC

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
    analyze_baseflow_indices(df::DataFrame; min_days=250, max_na_frac=0.2) -> Dict

Calculate baseflow index signatures with trend statistics.

Calculates 2 metrics:
- BFI_Eckhardt: Baseflow index using Eckhardt filter
- BFI_LyneHollick: Baseflow index using Lyne-Hollick filter

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, dowy
min_days : Int
    Minimum days per year for valid calculation
max_na_frac : Float64
    Maximum NA fraction allowed per year

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (2 metrics × 8 stats = 16 values)
"""
function analyze_baseflow_indices(
    df::DataFrame;
    min_days::Int=CFG_BASEFLOW_MIN_DAYS,
    max_na_frac::Float64=CFG_BASEFLOW_MAX_MISSING_FRAC
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

        if length(year_Q) < min_days
            continue
        end

        # Check NA fraction
        na_frac = sum(isnan.(year_Q)) / length(year_Q)
        if na_frac > max_na_frac
            continue
        end

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
        if sum(valid_q_mask) < min_days
            continue
        end

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
    result = generate_stats(annual_data; value_cols=metrics)

    return result
end

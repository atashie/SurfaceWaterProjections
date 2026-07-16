"""
Flashiness signature.

Richards-Baker flashiness index measuring flow variability.
"""

"""
    calculate_flashiness(Q::Vector) -> Float64

Calculate Richards-Baker flashiness index.

R-B Index = sum(|Q_i - Q_{i-1}|) / sum(Q_i)

Parameters
----------
Q : Vector
    Daily discharge values

Returns
-------
Float64
    Flashiness index (0 = constant flow, higher = more variable)
"""
function calculate_flashiness(Q::AbstractVector{<:Real})
    valid_mask = .!isnan.(Q)
    n_valid = sum(valid_mask)
    if n_valid < 2
        return NaN
    end

    # Preprocessor guarantees no NaNs in valid years
    Q_interp = Float64.(Q)

    total_Q = sum(Q_interp)
    if total_Q <= 0
        return NaN
    end

    sum_changes = sum(abs.(diff(Q_interp)))
    return sum_changes / total_Q
end


"""
    analyze_flashiness_trends(df::DataFrame) -> Dict

Calculate flashiness trends over time.

Calculates annual Richards-Baker flashiness index and trends.

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, dowy
    (preprocessor guarantees all years are valid)

Returns
-------
Dict{String, Float64}
    Dictionary with flashinessRB statistics (8 values)
"""
function analyze_flashiness_trends(df::DataFrame; trend_completeness::Union{Nothing, Float64}=nothing, decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing, changepoint::Union{Nothing, NamedTuple}=nothing, collector::Union{Nothing, AnnualCollector}=nothing)
    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year", "dowy"])
    if !valid
        @warn "analyze_flashiness_trends: Missing columns: $missing_cols"
        return empty_stats("flashinessRB")
    end

    if nrow(df) == 0
        return empty_stats("flashinessRB")
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
        year_Q_sorted = year_Q[sort_idx]

        # Calculate flashiness
        flashiness = calculate_flashiness(year_Q_sorted)

        if !isnan(flashiness)
            push!(annual_data, (water_year=yr, flashinessRB=flashiness))
        end
    end

    if nrow(annual_data) < 3
        return empty_stats("flashinessRB")
    end

    return generate_stats(annual_data; value_cols=["flashinessRB"], trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)
end

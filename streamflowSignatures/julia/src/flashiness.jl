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
    # Remove NaN values but keep track of indices for adjacency
    valid_mask = .!isnan.(Q)
    if sum(valid_mask) < 2
        return NaN
    end

    # Get valid values
    Q_valid = Q[valid_mask]
    total_Q = sum(Q_valid)

    if total_Q <= 0
        return NaN
    end

    # Sum of absolute day-to-day changes
    sum_changes = sum(abs.(diff(Q_valid)))

    return sum_changes / total_Q
end


"""
    analyze_flashiness_trends(df::DataFrame; min_days=250) -> Dict

Calculate flashiness trends over time.

Calculates annual Richards-Baker flashiness index and trends.

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, dowy
min_days : Int
    Minimum non-NA days per year for valid calculation

Returns
-------
Dict{String, Float64}
    Dictionary with flashinessRB statistics (8 values)
"""
function analyze_flashiness_trends(df::DataFrame; min_days::Int=250)
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

        # Check minimum non-NA data requirement
        n_valid = sum(.!isnan.(year_Q_sorted))
        if n_valid < min_days
            continue
        end

        # Calculate flashiness
        flashiness = calculate_flashiness(year_Q_sorted)

        if !isnan(flashiness)
            push!(annual_data, (water_year=yr, flashinessRB=flashiness))
        end
    end

    if nrow(annual_data) < 3
        return empty_stats("flashinessRB")
    end

    return generate_stats(annual_data; value_cols=["flashinessRB"])
end

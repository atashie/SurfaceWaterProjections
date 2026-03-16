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

    # Interpolate NaN values (matching R's approx() and Python's np.interp)
    # to preserve temporal adjacency — removing NaNs and compacting
    # creates artificial jumps between non-adjacent days
    Q_interp = Float64.(Q)  # Already creates a new array; no copy needed
    if n_valid < length(Q)
        valid_indices = findall(valid_mask)
        valid_values = Q_interp[valid_indices]
        for i in 1:length(Q_interp)
            if isnan(Q_interp[i])
                # Find surrounding valid indices for linear interpolation
                left = findlast(j -> j < i, valid_indices)
                right = findfirst(j -> j > i, valid_indices)
                if left === nothing && right !== nothing
                    Q_interp[i] = valid_values[right]  # extrapolate from nearest
                elseif right === nothing && left !== nothing
                    Q_interp[i] = valid_values[left]    # extrapolate from nearest
                elseif left !== nothing && right !== nothing
                    # Linear interpolation
                    li, ri = valid_indices[left], valid_indices[right]
                    frac = (i - li) / (ri - li)
                    Q_interp[i] = valid_values[left] + frac * (valid_values[right] - valid_values[left])
                end
            end
        end
    end

    total_Q = sum(Q_interp)
    if total_Q <= 0
        return NaN
    end

    sum_changes = sum(abs.(diff(Q_interp)))
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

        # Skip if too many missing values (>20%), matching R/Python
        na_frac = 1.0 - n_valid / length(year_Q_sorted)
        if na_frac > 0.2
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

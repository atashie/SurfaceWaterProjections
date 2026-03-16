"""
Flow timing signatures.

Cumulative flow timing metrics (D-day metrics) and peak timing.
"""

"""
    analyze_flow_timing_trends(df::DataFrame; min_days=300) -> Dict

Calculate flow timing signatures with trend statistics.

Calculates 13 metrics:
- D5_day through D95_day: Days when cumulative flow reaches X%
- D25_to_D75: Duration of middle 50% of flow
- Dmax: Day of maximum flow

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, dowy
min_days : Int
    Minimum days per year for valid calculation

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (13 metrics × 8 stats = 104 values)
"""
function analyze_flow_timing_trends(df::DataFrame; min_days::Int=300)
    result = Dict{String, Float64}()

    # D-day percentiles from config
    d_percentiles = CFG_D_PERCENTILES
    metrics = ["D$(p)_day" for p in d_percentiles]
    push!(metrics, "D25_to_D75", "Dmax")

    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year", "dowy"])
    if !valid
        @warn "analyze_flow_timing_trends: Missing columns: $missing_cols"
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

    years = unique(df.water_year)
    annual_data = DataFrame(
        water_year = years,
        D5_day = fill(NaN, length(years)),
        D10_day = fill(NaN, length(years)),
        D20_day = fill(NaN, length(years)),
        D30_day = fill(NaN, length(years)),
        D40_day = fill(NaN, length(years)),
        D50_day = fill(NaN, length(years)),
        D60_day = fill(NaN, length(years)),
        D70_day = fill(NaN, length(years)),
        D80_day = fill(NaN, length(years)),
        D90_day = fill(NaN, length(years)),
        D95_day = fill(NaN, length(years)),
        D25_to_D75 = fill(NaN, length(years)),
        Dmax = fill(NaN, length(years))
    )

    for yr in years
        year_mask = coalesce.(df.water_year .== yr, false)
        if !any(year_mask)
            continue
        end
        year_df = df[year_mask, :]

        if nrow(year_df) < min_days
            continue
        end

        yr_idx = findfirst(annual_data.water_year .== yr)
        yr_idx === nothing && continue

        # Sort by day of water year (handle Missing dowy)
        dowy_clean = [coalesce(d, 999) for d in year_df.dowy]
        perm = sortperm(dowy_clean)
        dowy_sorted = dowy_clean[perm]

        # Handle NaN/Missing in Q - replace with 0 for cumsum
        Q_values = coalesce_q(year_df.Q)[perm]
        Q_values[isnan.(Q_values)] .= 0.0

        # Cumulative flow
        cum_Q = cumsum(Q_values)
        total_Q = cum_Q[end]

        if total_Q <= 0
            continue
        end

        # Normalize to percentage
        cum_pct = cum_Q ./ total_Q .* 100

        # Find D-day for each percentile
        for p in d_percentiles
            # Find first day where cumulative flow >= p%
            idx = findfirst(cum_pct .>= p)
            if idx !== nothing
                annual_data[yr_idx, Symbol("D$(p)_day")] = dowy_sorted[idx]
            end
        end

        # D25 to D75 duration — compute directly from cumulative percentages
        # (NOT from D_PERCENTILES, which may not include 25 and 75)
        idx_25 = findfirst(cum_pct .>= 25)
        idx_75 = findfirst(cum_pct .>= 75)
        if idx_25 !== nothing && idx_75 !== nothing
            annual_data[yr_idx, :D25_to_D75] = dowy_sorted[idx_75] - dowy_sorted[idx_25]
        end

        # Day of maximum flow
        max_idx = argmax(Q_values)
        annual_data[yr_idx, :Dmax] = dowy_sorted[max_idx]
    end

    if nrow(annual_data) < 3
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    # Generate statistics for all metrics
    value_cols = filter(x -> x != "water_year", names(annual_data))
    result = generate_stats(annual_data; value_cols=value_cols)

    return result
end

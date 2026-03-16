"""
Runoff ratio signatures.

Q/P ratios (annual and seasonal). Requires precipitation data.
"""

# PPT thresholds loaded from config.jl (CFG_RUNOFF_MIN_ANNUAL_PPT, CFG_RUNOFF_MIN_SEASONAL_PPT)


"""
    analyze_Q_PPT_relationships(df::DataFrame) -> Dict

Calculate runoff ratio signatures with trend statistics.

Calculates 5 metrics:
- annual_runoff_ratio: Annual Q/P
- winter_runoff_ratio, spring_runoff_ratio, summer_runoff_ratio, fall_runoff_ratio

Each metric produces 8 statistics via generate_stats().
All years with any valid (Q, PPT) data are included (matching R/Python — no min_days filter).

Parameters
----------
df : DataFrame
    Daily streamflow and climate data with columns: water_year, Q, PPT, month

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (5 metrics × 8 stats = 40 values)
"""
function analyze_Q_PPT_relationships(df::DataFrame)
    result = Dict{String, Float64}()

    metrics = [
        "annual_runoff_ratio",
        "winter_runoff_ratio", "spring_runoff_ratio",
        "summer_runoff_ratio", "fall_runoff_ratio"
    ]

    # Check for PPT column
    if !("PPT" in names(df))
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
        annual_runoff_ratio = fill(NaN, length(years)),
        winter_runoff_ratio = fill(NaN, length(years)),
        spring_runoff_ratio = fill(NaN, length(years)),
        summer_runoff_ratio = fill(NaN, length(years)),
        fall_runoff_ratio = fill(NaN, length(years))
    )

    for yr in years
        year_mask = coalesce.(df.water_year .== yr, false)
        if !any(year_mask)
            continue
        end
        year_df = df[year_mask, :]

        # Convert Q and PPT to clean vectors
        Q_clean = coalesce_q(year_df.Q)
        PPT_clean = coalesce_q(year_df.PPT)

        yr_idx = findfirst(annual_data.water_year .== yr)
        yr_idx === nothing && continue

        # No min_days check — include ALL years with any valid data (matching R/Python)

        # Annual ratio - sum only days where BOTH Q and PPT are valid
        # Matches R's aggregate(cbind(Q, PPT), na.action=na.omit)
        valid_mask = .!isnan.(Q_clean) .& .!isnan.(PPT_clean)
        total_Q = sum(valid_mask) > 0 ? sum(Q_clean[valid_mask]) : 0.0
        total_P = sum(valid_mask) > 0 ? sum(PPT_clean[valid_mask]) : 0.0

        if total_P >= CFG_RUNOFF_MIN_ANNUAL_PPT
            annual_data[yr_idx, :annual_runoff_ratio] = total_Q / total_P
        end

        # Seasonal ratios
        for (season, months) in SEASONS
            month_mask = [coalesce(m, 0) in months for m in year_df.month]
            season_Q = Q_clean[month_mask]
            season_PPT = PPT_clean[month_mask]

            # No minimum-days check for seasonal ratios (matching R/Python)
            s_valid = .!isnan.(season_Q) .& .!isnan.(season_PPT)
            Q_season = sum(s_valid) > 0 ? sum(season_Q[s_valid]) : 0.0
            P_season = sum(s_valid) > 0 ? sum(season_PPT[s_valid]) : 0.0

            if P_season >= CFG_RUNOFF_MIN_SEASONAL_PPT
                annual_data[yr_idx, Symbol("$(season)_runoff_ratio")] = Q_season / P_season
            end
        end
    end

    # No early-return gate — let generate_stats() handle min_rows internally
    # (matching R/Python, which pass whatever data exists to generate_stats)
    result = generate_stats(annual_data; value_cols=metrics)

    return result
end

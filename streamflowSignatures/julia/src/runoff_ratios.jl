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
function analyze_Q_PPT_relationships(df::DataFrame; seasonal_flags::Union{Nothing, DataFrame}=nothing, trend_completeness::Union{Nothing, Float64}=nothing, decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing, changepoint::Union{Nothing, NamedTuple}=nothing, collector::Union{Nothing, AnnualCollector}=nothing)
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
        result["runoff_ratio_high_count"] = NaN
        return result
    end

    if nrow(df) == 0
        for m in metrics
            merge!(result, empty_stats(m))
        end
        result["runoff_ratio_high_count"] = NaN
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

        if total_P > CFG_RUNOFF_MIN_ANNUAL_PPT
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

            if P_season > CFG_RUNOFF_MIN_SEASONAL_PPT
                annual_data[yr_idx, Symbol("$(season)_runoff_ratio")] = Q_season / P_season
            end
        end
    end

    # Apply seasonal completeness flags: set incomplete seasons to NA
    if seasonal_flags !== nothing && nrow(seasonal_flags) > 0
        # Flag names must match the preprocessor's abbreviated seasonal columns
        # (win_/spr_/sum_/fal_complete — see io.jl seasonal_rows and the identical
        # lookup in flow_volumes.jl). Full-name lookups made this block dead code
        # until 2026-08-24 (CHANGELOG Known Issues).
        ratio_flag_map = Dict(
            :winter_runoff_ratio => :win_complete,
            :spring_runoff_ratio => :spr_complete,
            :summer_runoff_ratio => :sum_complete,
            :fall_runoff_ratio => :fal_complete
        )
        for (ratio_col, flag_col) in ratio_flag_map
            if String(ratio_col) in names(annual_data) && String(flag_col) in names(seasonal_flags)
                for i in 1:nrow(annual_data)
                    yr = annual_data.water_year[i]
                    flag_rows = filter(row -> row.water_year == yr, eachrow(seasonal_flags))
                    for fr in flag_rows
                        if haskey(fr, flag_col) && fr[flag_col] == false
                            annual_data[i, ratio_col] = NaN
                        end
                    end
                end
            end
        end
    end

    # No early-return gate — let generate_stats() handle min_rows internally
    # (matching R/Python, which pass whatever data exists to generate_stats)
    result = generate_stats(annual_data; value_cols=metrics, trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)

    # Count years where annual runoff ratio > 2.0 (per-gage scalar diagnostic)
    n_high_rr = count(x -> !isnan(x) && x > 2.0, annual_data.annual_runoff_ratio)
    result["runoff_ratio_high_count"] = Float64(n_high_rr)

    return result
end

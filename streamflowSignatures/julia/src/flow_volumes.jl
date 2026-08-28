"""
Flow volume signatures.

Annual and seasonal flow volumes plus percentile-based metrics.
"""

"""
    calculate_flow_vols_by_year(df::DataFrame) -> Dict

Calculate flow volume signatures with trend statistics.

Calculates 21 metrics:
- Qann: Annual total flow
- Qwin, Qspr, Qsum, Qfal: Seasonal totals
- Q1, Q5, Q10, Q20, Q25, Q30, Q40, Q50, Q60, Q70, Q75, Q80, Q90, Q95, Q99: Percentiles
- Q95_Q10: High-low flow difference (Q95 - Q10)

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, month, Q
    (preprocessor guarantees all years are valid)

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (21 metrics × 8 stats = 168 values)
"""
function calculate_flow_vols_by_year(df::DataFrame; seasonal_flags::Union{Nothing, DataFrame}=nothing, trend_completeness::Union{Nothing, Float64}=nothing, decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing, changepoint::Union{Nothing, NamedTuple}=nothing, collector::Union{Nothing, AnnualCollector}=nothing)
    result = Dict{String, Float64}()

    # Define metrics for empty results
    metrics = ["Qann", "Qwin", "Qspr", "Qsum", "Qfal",
               "Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40",
               "Q50", "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99",
               "Q95_Q10"]

    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year", "month"])
    if !valid
        @warn "calculate_flow_vols_by_year: Missing columns: $missing_cols"
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

    # Season definitions (by calendar month)
    seasons = Dict(
        "win" => [12, 1, 2],
        "spr" => [3, 4, 5],
        "sum" => [6, 7, 8],
        "fal" => [9, 10, 11]
    )

    # Percentiles from config
    percentiles = CFG_FLOW_PERCENTILES

    # Convert Q to Float64, handling missing values
    Q_clean = coalesce_q(df.Q)

    # Group by water year — pre-allocate DataFrame
    years = unique(df.water_year)
    annual_data = DataFrame(
        water_year = years,
        Qann = fill(NaN, length(years)),
        Qwin = fill(NaN, length(years)),
        Qspr = fill(NaN, length(years)),
        Qsum = fill(NaN, length(years)),
        Qfal = fill(NaN, length(years)),
        Q1 = fill(NaN, length(years)),
        Q5 = fill(NaN, length(years)),
        Q10 = fill(NaN, length(years)),
        Q20 = fill(NaN, length(years)),
        Q25 = fill(NaN, length(years)),
        Q30 = fill(NaN, length(years)),
        Q40 = fill(NaN, length(years)),
        Q50 = fill(NaN, length(years)),
        Q60 = fill(NaN, length(years)),
        Q70 = fill(NaN, length(years)),
        Q75 = fill(NaN, length(years)),
        Q80 = fill(NaN, length(years)),
        Q90 = fill(NaN, length(years)),
        Q95 = fill(NaN, length(years)),
        Q99 = fill(NaN, length(years)),
        Q95_Q10 = fill(NaN, length(years))
    )

    for yr in years
        year_mask = coalesce.(df.water_year .== yr, false)
        if !any(year_mask)
            continue
        end
        year_Q = Q_clean[year_mask]
        year_month = Int[coalesce(m, 0) for m in df.month[year_mask]]

        yr_idx = findfirst(annual_data.water_year .== yr)
        yr_idx === nothing && continue

        # Get valid (non-NaN) Q values for this year
        Q_valid = filter(!isnan, year_Q)

        # Annual total (sum of daily mm/day values = mm)
        annual_data[yr_idx, :Qann] = sum(Q_valid)

        # Seasonal totals
        for (season, season_months) in seasons
            season_mask = in.(year_month, Ref(Set(season_months)))
            season_Q = year_Q[season_mask]
            season_valid = filter(!isnan, season_Q)
            if length(season_valid) > 0
                annual_data[yr_idx, Symbol("Q$season")] = sum(season_valid)
            end
        end

        # Percentiles (based on daily values)
        for p in percentiles
            annual_data[yr_idx, Symbol("Q$p")] = quantile(Q_valid, p / 100)
        end

        # Q95-Q10 difference
        annual_data[yr_idx, :Q95_Q10] = annual_data[yr_idx, :Q95] - annual_data[yr_idx, :Q10]
    end

    # Apply seasonal_flags: set incomplete seasons to NaN
    if seasonal_flags !== nothing && nrow(seasonal_flags) > 0
        season_col_map = Dict("Qwin" => :win_complete, "Qspr" => :spr_complete,
                              "Qsum" => :sum_complete, "Qfal" => :fal_complete)
        for (qcol, flag_col) in season_col_map
            if Symbol(qcol) in propertynames(annual_data) && flag_col in propertynames(seasonal_flags)
                for i in 1:nrow(annual_data)
                    wy = annual_data[i, :water_year]
                    flag_row = findfirst(seasonal_flags.water_year .== wy)
                    if flag_row !== nothing && !seasonal_flags[flag_row, flag_col]
                        annual_data[i, Symbol(qcol)] = NaN
                    end
                end
            end
        end
    end

    if nrow(annual_data) < 3
        metrics = ["Qann", "Qwin", "Qspr", "Qsum", "Qfal",
                   "Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40",
                   "Q50", "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99",
                   "Q95_Q10"]
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    # Generate statistics for all metrics
    value_cols = filter(x -> x != "water_year", names(annual_data))
    result = generate_stats(annual_data; value_cols=value_cols, trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)

    return result
end

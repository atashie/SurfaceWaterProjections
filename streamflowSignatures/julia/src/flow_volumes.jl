"""
Flow volume signatures.

Annual and seasonal flow volumes plus percentile-based metrics.
"""

"""
    calculate_flow_vols_by_year(df::DataFrame; min_days=250) -> Dict

Calculate flow volume signatures with trend statistics.

Calculates 22 metrics:
- Qann: Annual total flow
- Qwin, Qspr, Qsum, Qfal: Seasonal totals
- Q1, Q5, Q10, Q20, Q25, Q30, Q40, Q50, Q60, Q70, Q75, Q80, Q90, Q95, Q99: Percentiles
- Q95_Q10: High-low flow ratio

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, month, Q
min_days : Int
    Minimum days per year for valid calculation

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (22 metrics × 8 stats = 176 values)
"""
function calculate_flow_vols_by_year(df::DataFrame; min_days::Int=250)
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

        # Check minimum data requirement - count non-NaN values
        n_valid = sum(.!isnan.(year_Q))
        if n_valid < min_days
            continue
        end

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
    result = generate_stats(annual_data; value_cols=value_cols)

    return result
end

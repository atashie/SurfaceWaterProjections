"""
Flow timing signatures.

Cumulative flow timing metrics (D-day metrics) and peak timing.
"""

"""
    analyze_flow_timing_trends(df::DataFrame) -> Dict

Calculate flow timing signatures with trend statistics.

Calculates 15 metrics:
- D1_day through D99_day: Days when cumulative flow reaches X% (13 D-day
  percentiles from config `timing.d_percentiles`)
- D25_to_D75: Duration of middle 50% of flow
- Dmax: Day of maximum flow

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, dowy
    (preprocessor guarantees all years are valid)

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (15 metrics × 8 stats = 120 values)
"""
function analyze_flow_timing_trends(df::DataFrame; trend_completeness::Union{Nothing, Float64}=nothing, decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing, changepoint::Union{Nothing, NamedTuple}=nothing, collector::Union{Nothing, AnnualCollector}=nothing)
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
    # Build DataFrame dynamically from config percentiles
    annual_data = DataFrame(water_year = years)
    for p in d_percentiles
        annual_data[!, Symbol("D$(p)_day")] = fill(NaN, length(years))
    end
    annual_data[!, :D25_to_D75] = fill(NaN, length(years))
    annual_data[!, :Dmax] = fill(NaN, length(years))

    for yr in years
        year_mask = coalesce.(df.water_year .== yr, false)
        if !any(year_mask)
            continue
        end
        year_df = df[year_mask, :]

        yr_idx = findfirst(annual_data.water_year .== yr)
        yr_idx === nothing && continue

        # Sort by day of water year (handle Missing dowy)
        dowy_clean = [coalesce(d, 999) for d in year_df.dowy]
        perm = sortperm(dowy_clean)
        dowy_sorted = dowy_clean[perm]

        # Handle NaN/Missing in Q - skip year if any NaNs remain
        # (preprocess_daily_data guarantees zero NAs for valid years)
        Q_values = coalesce_q(year_df.Q)[perm]
        n_na = count(isnan, Q_values)
        if n_na > 0
            continue
        end

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
    result = generate_stats(annual_data; value_cols=value_cols, trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)

    return result
end

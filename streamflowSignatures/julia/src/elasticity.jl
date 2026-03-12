"""
Streamflow elasticity signatures.

Sensitivity of streamflow to precipitation changes.
"""

# Configuration loaded from config.jl:
# CFG_ELASTICITY_WINDOW_YEARS, CFG_ELASTICITY_MIN_YEARS,
# CFG_ELASTICITY_MIN_DATA_COMPLETENESS, CFG_ELASTICITY_MIN_ANNUAL_PPT


"""
    calculate_streamflow_elasticity(df::DataFrame; min_years=10, window=11) -> Dict

Calculate streamflow elasticity signatures.

Calculates:
- elasticity_static: Overall catchment elasticity (single value)
- elasticity: Rolling window elasticity trends (8 statistics)

Elasticity E = (dQ/dP) / (Q_mean/P_mean)
- E ~ 1.0: Proportional response
- E > 1.0: Amplified response
- E < 1.0: Dampened response

Parameters
----------
df : DataFrame
    Daily streamflow and climate data with columns: water_year, Q, PPT
min_years : Int
    Minimum years required for calculation
window : Int
    Rolling window size for trend calculation

Returns
-------
Dict{String, Float64}
    Dictionary with elasticity_static and elasticity statistics
"""
function calculate_streamflow_elasticity(
    df::DataFrame;
    min_years::Int=CFG_ELASTICITY_MIN_YEARS,
    window::Int=CFG_ELASTICITY_WINDOW_YEARS,
    min_completeness::Real=CFG_ELASTICITY_MIN_DATA_COMPLETENESS,
    min_annual_ppt::Real=CFG_ELASTICITY_MIN_ANNUAL_PPT
)
    result = Dict{String, Float64}()

    # Check for PPT column
    if !("PPT" in names(df))
        result["elasticity_static"] = NaN
        merge!(result, empty_stats("elasticity"))
        return result
    end

    if nrow(df) == 0
        result["elasticity_static"] = NaN
        merge!(result, empty_stats("elasticity"))
        return result
    end

    # Aggregate to annual totals
    annual_data = combine(
        groupby(df, :water_year),
        :Q => (x -> sum(skipmissing(x); init=0.0)) => :Q_annual,
        :PPT => (x -> sum(skipmissing(x); init=0.0)) => :P_annual,
        [:Q, :PPT] => ((q, p) -> count(i -> !ismissing(q[i]) && !isnan(q[i]) && !ismissing(p[i]) && !isnan(p[i]), 1:length(q))) => :n_valid_days
    )

    # Calculate expected days per water year (accounting for leap years, like R/Python)
    function expected_days(wy::Real)
        # Convert to Int (water_year might be Float64 in DataFrame)
        wy_int = Int(wy)
        # Water year ends in year wy, so check if wy is a leap year
        is_leap = (wy_int % 4 == 0 && wy_int % 100 != 0) || (wy_int % 400 == 0)
        return is_leap ? 366 : 365
    end

    # Filter years with sufficient data (>= min_completeness of expected days)
    valid_years = [
        coalesce(annual_data.n_valid_days[i] >= min_completeness * expected_days(annual_data.water_year[i]), false)
        for i in 1:nrow(annual_data)
    ]
    annual_data = annual_data[valid_years, :]

    if nrow(annual_data) < min_years
        result["elasticity_static"] = NaN
        merge!(result, empty_stats("elasticity"))
        return result
    end

    # Sort by year
    sort!(annual_data, :water_year)

    Q = annual_data.Q_annual
    P = annual_data.P_annual

    # Filter out years with near-zero P (using config value, like R/Python)
    valid_mask = P .>= min_annual_ppt
    Q_valid = Q[valid_mask]
    P_valid = P[valid_mask]

    if length(Q_valid) < min_years
        result["elasticity_static"] = NaN
        merge!(result, empty_stats("elasticity"))
        return result
    end

    # Calculate static elasticity using per-year anomaly method (matching R/Python)
    Q_mean = mean(Q_valid)
    P_mean = mean(P_valid)

    if P_mean <= 0 || Q_mean <= 0
        result["elasticity_static"] = NaN
        merge!(result, empty_stats("elasticity"))
        return result
    end

    # FIX: Per-year anomaly method (matching R/Python)
    # E_i = (dQ_i / dP_i) / (Q_mean / P_mean)
    # where dQ_i = Q_i - Q_mean, dP_i = P_i - P_mean
    elasticity_values = Float64[]
    min_dP = 0.1  # Avoid division by very small P changes

    for i in 1:length(Q_valid)
        dQ = Q_valid[i] - Q_mean
        dP = P_valid[i] - P_mean

        if abs(dP) > min_dP
            e = (dQ / dP) / (Q_mean / P_mean)
            push!(elasticity_values, e)
        end
    end

    if isempty(elasticity_values)
        result["elasticity_static"] = NaN
    else
        elasticity_static = median(elasticity_values)
        result["elasticity_static"] = elasticity_static
    end

    # Rolling window elasticity using same per-year anomaly method
    n = length(Q_valid)
    if n < window + 3
        merge!(result, empty_stats("elasticity"))
        return result
    end

    years_valid = annual_data.water_year[valid_mask]
    rolling_elasticity = Float64[]
    rolling_years = Int[]  # Will convert Float64 years to Int when pushing

    # FIX: Use END of window for year assignment (matching R/Python)
    # R: water_year = annual$water_year[(rolling_window):n]
    # Python: rolling_years.append(annual.iloc[end_idx]["water_year"])
    for end_idx in window:n
        start_idx = end_idx - window + 1
        Q_window = Q_valid[start_idx:end_idx]
        P_window = P_valid[start_idx:end_idx]

        Q_mean_w = mean(Q_window)
        P_mean_w = mean(P_window)

        if P_mean_w > 0 && Q_mean_w > 0
            # Calculate per-year elasticities within window
            window_elasticities = Float64[]
            for j in 1:length(Q_window)
                dQ = Q_window[j] - Q_mean_w
                dP = P_window[j] - P_mean_w

                if abs(dP) > min_dP
                    e = (dQ / dP) / (Q_mean_w / P_mean_w)
                    push!(window_elasticities, e)
                end
            end

            if !isempty(window_elasticities)
                push!(rolling_elasticity, median(window_elasticities))
                push!(rolling_years, Int(years_valid[end_idx]))  # Use END of window, convert to Int
            end
        end
    end

    if length(rolling_elasticity) >= 3
        rolling_df = DataFrame(
            water_year = rolling_years,
            elasticity = rolling_elasticity
        )
        merge!(result, generate_stats(rolling_df; value_cols=["elasticity"]))
    else
        merge!(result, empty_stats("elasticity"))
    end

    return result
end

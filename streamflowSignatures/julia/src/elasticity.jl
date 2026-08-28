"""
Streamflow elasticity signatures.

Sensitivity of streamflow to precipitation changes.
"""

# Configuration loaded from config.jl:
# CFG_ELASTICITY_WINDOW_YEARS, CFG_ELASTICITY_MIN_YEARS,
# CFG_ELASTICITY_MIN_ANNUAL_PPT


"""
    calculate_streamflow_elasticity(df::DataFrame; min_years=15, window=11) -> Dict

Calculate streamflow elasticity signatures.

Calculates:
- elasticity_static: Overall catchment elasticity (single value)
- elasticity_rolling: Rolling window (11-year) elasticity trends (8 statistics)
  Uses departure-from-mean method (Sawicz et al. 2011)
- elasticity_annual: Year-over-year elasticity trends (8 statistics)
  Uses consecutive-year differences for year-to-year sensitivity
- elasticity_years_total: Total years available for this gage (diagnostic scalar)
- elasticity_years_low_ppt: Years excluded due to low PPT (diagnostic scalar)

Elasticity E = (dQ/dP) / (Q_mean/P_mean)
- E ~ 1.0: Proportional response
- E > 1.0: Amplified response
- E < 1.0: Dampened response

Parameters
----------
df : DataFrame
    Daily streamflow and climate data with columns: water_year, Q, PPT
min_years : Int
    Minimum years required for calculation (default CFG_ELASTICITY_MIN_YEARS
    from config `elasticity.min_years`; 15 as shipped)
window : Int
    Rolling window size for trend calculation

Returns
-------
Dict{String, Float64}
    Dictionary with elasticity_static, elasticity_rolling, elasticity_annual statistics,
    and diagnostic scalars
"""
function calculate_streamflow_elasticity(
    df::DataFrame;
    min_years::Int=CFG_ELASTICITY_MIN_YEARS,
    window::Int=CFG_ELASTICITY_WINDOW_YEARS,
    min_annual_ppt::Real=CFG_ELASTICITY_MIN_ANNUAL_PPT,
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing,
    changepoint::Union{Nothing, NamedTuple}=nothing,
    collector::Union{Nothing, AnnualCollector}=nothing
)
    result = Dict{String, Float64}()

    # Helper to fill all empty outputs
    function fill_empty_elasticity!(r)
        r["elasticity_static"] = NaN
        merge!(r, empty_stats("elasticity_rolling"))
        merge!(r, empty_stats("elasticity_annual"))
        r["elasticity_years_total"] = NaN
        r["elasticity_years_low_ppt"] = NaN
        return r
    end

    # Check for PPT column
    if !("PPT" in names(df))
        fill_empty_elasticity!(result)
        return result
    end

    if nrow(df) == 0
        fill_empty_elasticity!(result)
        return result
    end

    # Aggregate to annual totals (preprocessor guarantees all years are valid)
    annual_data = combine(
        groupby(df, :water_year),
        :Q => (x -> sum(skipmissing(x); init=0.0)) => :Q_annual,
        :PPT => (x -> sum(skipmissing(x); init=0.0)) => :P_annual
    )

    if nrow(annual_data) < min_years
        fill_empty_elasticity!(result)
        return result
    end

    # Sort by year
    sort!(annual_data, :water_year)

    Q = annual_data.Q_annual
    P = annual_data.P_annual

    # Filter out years with near-zero P (using config value, like R/Python)
    # Must use strict > to match R/Python (not >=)
    valid_mask = P .> min_annual_ppt
    Q_valid = Q[valid_mask]
    P_valid = P[valid_mask]

    # Step 6: Elasticity data quality diagnostics (per-gage scalars)
    result["elasticity_years_total"] = Float64(nrow(annual_data))
    result["elasticity_years_low_ppt"] = Float64(sum(.!valid_mask))

    if length(Q_valid) < min_years
        result["elasticity_static"] = NaN
        merge!(result, empty_stats("elasticity_rolling"))
        merge!(result, empty_stats("elasticity_annual"))
        return result
    end

    # Calculate static elasticity using per-year anomaly method (matching R/Python)
    Q_mean = mean(Q_valid)
    P_mean = mean(P_valid)

    if P_mean <= 0 || Q_mean <= 0
        result["elasticity_static"] = NaN
        merge!(result, empty_stats("elasticity_rolling"))
        merge!(result, empty_stats("elasticity_annual"))
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

    # Rolling window elasticity using same per-year anomaly method (Sawicz et al. 2011)
    n = length(Q_valid)
    years_valid = annual_data.water_year[valid_mask]

    if n < window + 3
        merge!(result, empty_stats("elasticity_rolling"))
    else
        rolling_elasticity = Float64[]
        rolling_years = Int[]

        # Use END of window for year assignment (matching R/Python)
        for end_idx in window:n
            start_idx = end_idx - window + 1
            Q_window = Q_valid[start_idx:end_idx]
            P_window = P_valid[start_idx:end_idx]

            Q_mean_w = mean(Q_window)
            P_mean_w = mean(P_window)

            if P_mean_w > 0 && Q_mean_w > 0
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
                    push!(rolling_years, Int(years_valid[end_idx]))
                end
            end
        end

        if length(rolling_elasticity) >= 3
            rolling_df = DataFrame(
                water_year = rolling_years,
                elasticity_rolling = rolling_elasticity
            )
            merge!(result, generate_stats(rolling_df; value_cols=["elasticity_rolling"], trend_completeness=trend_completeness, decade_completeness=decade_completeness, changepoint=changepoint, collector=collector))
        else
            merge!(result, empty_stats("elasticity_rolling"))
        end
    end

    # Year-over-year elasticity (NOT Sawicz — uses consecutive-year differences)
    # E_t = ((Q_t - Q_{t-1}) / (P_t - P_{t-1})) / (Q_mean / P_mean)
    # Measures year-to-year sensitivity rather than deviation-from-mean sensitivity
    yoy_elasticity = Float64[]
    yoy_years = Int[]
    for i in 2:length(Q_valid)
        dQ = Q_valid[i] - Q_valid[i - 1]
        dP = P_valid[i] - P_valid[i - 1]
        if abs(dP) > min_dP
            e = (dQ / dP) / (Q_mean / P_mean)
            push!(yoy_elasticity, e)
            push!(yoy_years, Int(years_valid[i]))
        end
    end

    if length(yoy_elasticity) >= 3
        yoy_df = DataFrame(
            water_year = yoy_years,
            elasticity_annual = yoy_elasticity
        )
        merge!(result, generate_stats(yoy_df; value_cols=["elasticity_annual"], trend_completeness=trend_completeness, decade_completeness=decade_completeness, changepoint=changepoint, collector=collector))
    else
        merge!(result, empty_stats("elasticity_annual"))
    end

    return result
end

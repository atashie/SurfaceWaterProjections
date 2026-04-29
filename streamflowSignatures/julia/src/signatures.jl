"""
High-level convenience function for calculating all streamflow signatures at once.
"""

"""
    calculate_all_signatures(gage_data::DataFrame, has_climate::Bool=false; gage_id::String="unknown") -> Dict{String, Any}

Calculate all signatures for a single gage.

Calls each signature function in sequence, catching exceptions per-signature
so that a failure in one does not prevent others from being calculated.

# Arguments
- `gage_data::DataFrame`: DataFrame for a single gage with columns: gage_id, date, Q, water_year, month, dowy.
- `has_climate::Bool`: Whether climate data (PPT column) is available.
- `gage_id::String`: Gage identifier for warning messages.

# Returns
- `Dict{String, Any}`: Dictionary mapping signature column names to values.
"""
function calculate_all_signatures(
    gage_data::DataFrame, has_climate::Bool=false;
    gage_id::String="unknown",
    seasonal_flags::Union{Nothing, DataFrame}=nothing,
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing,
    climate_data::Union{Nothing, DataFrame}=nothing,
    include_qa_flags::Bool=false,
    changepoint::Union{Nothing, NamedTuple}=nothing
)::Dict{String, Any}
    results = Dict{String, Any}()

    # Season exclusion year counts (per-gage scalar diagnostics)
    if seasonal_flags !== nothing && nrow(seasonal_flags) > 0
        for (season, col) in [("winter", :win_complete), ("spring", :spr_complete),
                               ("summer", :sum_complete), ("fall", :fal_complete)]
            if hasproperty(seasonal_flags, col)
                results["season_excluded_years_$(season)"] = Float64(sum(.!seasonal_flags[!, col]))
            else
                results["season_excluded_years_$(season)"] = 0.0
            end
        end
    end

    tc = trend_completeness
    dc = decade_completeness
    cp = changepoint

    # Non-climate signatures — apply trend completeness to signatures that
    # produce one value per year. Skip for inherently sparse signatures
    # (recession is event-based, elasticity uses rolling windows).
    try
        merge!(results, calculate_flow_vols_by_year(gage_data; seasonal_flags=seasonal_flags, trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "calculate_flow_vols_by_year failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    try
        merge!(results, analyze_flashiness_trends(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "analyze_flashiness_trends failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    try
        merge!(results, analyze_flow_timing_trends(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "analyze_flow_timing_trends failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    try
        merge!(results, analyze_fdc_trends(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "analyze_fdc_trends failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    try
        merge!(results, analyze_baseflow_indices(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "analyze_baseflow_indices failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    # Recession: inherently sparse (event-based, many years with no events) — no trend completeness
    local recession_alpha = NaN
    try
        recession_results = analyze_recession_parameters(gage_data; changepoint=cp)
        # Extract scalar for parameterized BFI before merging
        recession_alpha = get(recession_results, "recession_alpha_point_cloud_linear_reservoir", NaN)
        merge!(results, recession_results)
    catch e
        @warn "analyze_recession_parameters failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    # Parameterized BFI using recession-derived alpha (requires recession to have run first)
    # Uses trend completeness (same as fixed-parameter BFI)
    try
        merge!(results, analyze_baseflow_indices_with_parameters(gage_data, recession_alpha;
            trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "analyze_baseflow_indices_with_parameters failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    try
        merge!(results, calculate_pulse_metrics(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "calculate_pulse_metrics failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    try
        merge!(results, calculate_negative_days(gage_data; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
    catch e
        @warn "calculate_negative_days failed for gage $gage_id" exception=(e, catch_backtrace())
    end

    # Climate-dependent signatures — use climate_data if provided (filtered to valid climate years),
    # otherwise fall back to gage_data (legacy path)
    if has_climate
        cdata = climate_data !== nothing ? climate_data : gage_data
        if "PPT" in names(cdata)
            try
                merge!(results, analyze_Q_PPT_relationships(cdata; seasonal_flags=seasonal_flags, trend_completeness=tc, decade_completeness=dc, changepoint=cp))
            catch e
                @warn "analyze_Q_PPT_relationships failed for gage $gage_id" exception=(e, catch_backtrace())
            end

            # Elasticity: rolling window produces fewer values than years — no trend completeness
            try
                merge!(results, calculate_streamflow_elasticity(cdata; changepoint=cp))
            catch e
                @warn "calculate_streamflow_elasticity failed for gage $gage_id" exception=(e, catch_backtrace())
            end

            try
                merge!(results, calculate_qp_seasonality(cdata; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
            catch e
                @warn "calculate_qp_seasonality failed for gage $gage_id" exception=(e, catch_backtrace())
            end

            try
                merge!(results, calculate_average_storage(cdata; trend_completeness=tc, decade_completeness=dc, changepoint=cp))
            catch e
                @warn "calculate_average_storage failed for gage $gage_id" exception=(e, catch_backtrace())
            end
        end
    end

    # QA/QC flags (optional)
    if include_qa_flags
        try
            row_df = DataFrame(Dict(k => [v] for (k, v) in results))
            flagged_df = compute_qa_flags(row_df)
            for col in get_flag_columns()
                if col in names(flagged_df)
                    results[col] = flagged_df[1, col]
                end
            end
        catch e
            @warn "QA/QC flags failed for gage $gage_id" exception=(e, catch_backtrace())
        end
    end

    return results
end

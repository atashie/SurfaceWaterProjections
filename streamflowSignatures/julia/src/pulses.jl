"""
Pulse metrics signatures.

High/low pulse counts, durations, TQmean, and flow reversals.
"""

# Season definitions
const SEASONS = Dict(
    "winter" => [12, 1, 2],
    "spring" => [3, 4, 5],
    "summer" => [6, 7, 8],
    "fall" => [9, 10, 11]
)


"""
    count_pulses(Q::Vector, threshold::Real, above::Bool) -> (count, mean_duration)

Count pulse events and calculate mean duration.

Parameters
----------
Q : Vector
    Daily discharge values
threshold : Real
    Threshold value for pulse detection
above : Bool
    If true, count pulses above threshold; if false, below

Returns
-------
Tuple{Int, Float64}
    (number of pulses, mean duration in days)
"""
function count_pulses(Q::AbstractVector{<:Real}, threshold::Real, above::Bool)
    n = length(Q)
    if n == 0
        return (0, NaN)
    end

    # Determine which days are in pulse state
    if above
        in_pulse = Q .> threshold
    else
        in_pulse = Q .< threshold
    end

    # Handle NaN - not in pulse
    in_pulse[isnan.(Q)] .= false

    # Count pulse events and durations
    pulse_count = 0
    durations = Int[]
    current_duration = 0

    for i in 1:n
        if in_pulse[i]
            current_duration += 1
        else
            if current_duration > 0
                pulse_count += 1
                push!(durations, current_duration)
                current_duration = 0
            end
        end
    end

    # Handle pulse at end of series
    if current_duration > 0
        pulse_count += 1
        push!(durations, current_duration)
    end

    mean_duration = isempty(durations) ? NaN : mean(durations)

    return (pulse_count, mean_duration)
end


"""
    count_flow_reversals(Q::Vector; threshold_pct=0.02) -> Int

Count direction changes in flow that exceed a magnitude threshold.

A reversal occurs when flow changes from increasing to decreasing or vice versa,
and the magnitude of the change exceeds threshold_pct of the current flow value.

Parameters
----------
Q : Vector
    Daily discharge values
threshold_pct : Float64
    Minimum magnitude threshold as fraction of current flow (default 0.02 = 2%)
"""
function count_flow_reversals(Q::AbstractVector{<:Real}; threshold_pct::Float64=0.02)
    # Preprocessor guarantees no NaNs in valid years
    n = length(Q)
    if n < 3
        return 0
    end

    Q_interp = Float64.(Q)

    reversals = 0

    # Match R/Python implementation: check prev_change and next_change
    # Only count reversal if |next_change| > threshold
    @inbounds for i in 2:(n-1)
        prev_change = Q_interp[i] - Q_interp[i-1]
        next_change = Q_interp[i+1] - Q_interp[i]
        threshold = abs(threshold_pct * Q_interp[i])

        if abs(next_change) > threshold
            # Peak (rising to falling)
            if prev_change >= 0 && next_change < 0
                reversals += 1
            # Trough (falling to rising)
            elseif prev_change <= 0 && next_change > 0
                reversals += 1
            end
        end
    end

    return reversals
end


"""
    calculate_tqmean(Q::Vector) -> Float64

Calculate percentage of days with flow above mean.
"""
function calculate_tqmean(Q::AbstractVector{<:Real})
    Q_valid = filter(!isnan, Q)
    if isempty(Q_valid)
        return NaN
    end

    mean_Q = mean(Q_valid)
    days_above = sum(Q_valid .> mean_Q)

    return days_above / length(Q_valid) * 100
end


"""
    calculate_pulse_metrics(df::DataFrame) -> Dict

Calculate pulse metric signatures with trend statistics.

Calculates 14 metrics:
- n_high_pulses_year, n_low_pulses_year: Pulse counts (year-specific thresholds)
- dur_high_pulses_year, dur_low_pulses_year: Mean pulse durations (year-specific)
- n_high_pulses_all, n_low_pulses_all: Pulse counts (period-of-record thresholds)
- dur_high_pulses_all, dur_low_pulses_all: Mean pulse durations (period-of-record)
- TQmean: Percentage of days above mean
- Flow_Reversals_annual: Total annual reversals
- Flow_Reversals_winter/spring/summer/fall: Seasonal reversals

The *_year metrics use year-specific Q90/Q10 thresholds, while *_all metrics
use Q90/Q10 calculated over the entire period of record.

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, month, dowy
    (preprocessor guarantees all years are valid)

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics
"""
function calculate_pulse_metrics(df::DataFrame; trend_completeness::Union{Nothing, Float64}=nothing, decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing, changepoint::Union{Nothing, NamedTuple}=nothing, collector::Union{Nothing, AnnualCollector}=nothing)
    result = Dict{String, Float64}()

    metrics = [
        "n_high_pulses_year", "n_low_pulses_year",
        "dur_high_pulses_year", "dur_low_pulses_year",
        # Period-of-record threshold metrics (matching R's *_all metrics)
        "n_high_pulses_all", "n_low_pulses_all",
        "dur_high_pulses_all", "dur_low_pulses_all",
        "TQmean", "Flow_Reversals_annual",
        "Flow_Reversals_winter", "Flow_Reversals_spring",
        "Flow_Reversals_summer", "Flow_Reversals_fall"
    ]

    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year", "month", "dowy"])
    if !valid
        @warn "calculate_pulse_metrics: Missing columns: $missing_cols"
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

    # Calculate period-of-record thresholds for *_all metrics
    Q_all = coalesce_q(df.Q)
    Q_all_valid = filter(!isnan, Q_all)
    if length(Q_all_valid) < 30
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end
    q90_all = quantile(Q_all_valid, 0.90)
    q10_all = quantile(Q_all_valid, 0.10)

    years = unique(df.water_year)
    annual_data = DataFrame(
        water_year = years,
        n_high_pulses_year = fill(NaN, length(years)),
        n_low_pulses_year = fill(NaN, length(years)),
        dur_high_pulses_year = fill(NaN, length(years)),
        dur_low_pulses_year = fill(NaN, length(years)),
        n_high_pulses_all = fill(NaN, length(years)),
        n_low_pulses_all = fill(NaN, length(years)),
        dur_high_pulses_all = fill(NaN, length(years)),
        dur_low_pulses_all = fill(NaN, length(years)),
        TQmean = fill(NaN, length(years)),
        Flow_Reversals_annual = fill(NaN, length(years)),
        Flow_Reversals_winter = fill(NaN, length(years)),
        Flow_Reversals_spring = fill(NaN, length(years)),
        Flow_Reversals_summer = fill(NaN, length(years)),
        Flow_Reversals_fall = fill(NaN, length(years))
    )

    for yr in years
        year_mask = df.water_year .== yr
        if !any(skipmissing(year_mask))
            continue
        end
        year_df = df[coalesce.(year_mask, false), :]
        year_df = sort(year_df, :dowy)

        # Convert Q to clean Float64 vector
        Q_year = coalesce_q(year_df.Q)

        # FIX: Calculate year-specific thresholds (matching R/Python)
        Q_valid = filter(!isnan, Q_year)
        if length(Q_valid) < 30
            continue
        end
        high_threshold = quantile(Q_valid, 0.90)
        low_threshold = quantile(Q_valid, 0.10)

        yr_idx = findfirst(annual_data.water_year .== yr)
        yr_idx === nothing && continue

        # High and low pulses using year-specific thresholds
        n_high, dur_high = count_pulses(Q_year, high_threshold, true)
        n_low, dur_low = count_pulses(Q_year, low_threshold, false)

        annual_data[yr_idx, :n_high_pulses_year] = n_high
        annual_data[yr_idx, :n_low_pulses_year] = n_low
        annual_data[yr_idx, :dur_high_pulses_year] = dur_high
        annual_data[yr_idx, :dur_low_pulses_year] = dur_low

        # High and low pulses using period-of-record thresholds (*_all metrics)
        n_high_all, dur_high_all = count_pulses(Q_year, q90_all, true)
        n_low_all, dur_low_all = count_pulses(Q_year, q10_all, false)

        annual_data[yr_idx, :n_high_pulses_all] = n_high_all
        annual_data[yr_idx, :n_low_pulses_all] = n_low_all
        annual_data[yr_idx, :dur_high_pulses_all] = dur_high_all
        annual_data[yr_idx, :dur_low_pulses_all] = dur_low_all

        # TQmean
        annual_data[yr_idx, :TQmean] = calculate_tqmean(Q_year)

        # Flow reversals - annual
        annual_data[yr_idx, :Flow_Reversals_annual] = count_flow_reversals(Q_year)

        # Seasonal reversals
        for (season, months) in SEASONS
            month_mask = [coalesce(m, 0) in months for m in year_df.month]
            if sum(month_mask) >= 30
                season_Q = Q_year[month_mask]
                annual_data[yr_idx, Symbol("Flow_Reversals_$season")] = count_flow_reversals(season_Q)
            end
        end
    end

    if nrow(annual_data) < 3
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    # Generate statistics for all metrics
    result = generate_stats(annual_data; value_cols=metrics, trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)

    return result
end


"""
    calculate_negative_days(df; trend_completeness=nothing, decade_completeness=nothing)

Count days with negative streamflow (Q < 0) per water year.
"""
function calculate_negative_days(df::DataFrame; trend_completeness::Union{Nothing, Float64}=nothing, decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing, changepoint::Union{Nothing, NamedTuple}=nothing, collector::Union{Nothing, AnnualCollector}=nothing)
    result = Dict{String, Float64}()

    valid, missing_cols = validate_columns(df, ["Q", "water_year"])
    if !valid
        @warn "calculate_negative_days: Missing columns: $missing_cols"
        merge!(result, empty_stats("negative_ann"))
        return result
    end

    annual_data = combine(groupby(df, :water_year),
        :Q => (q -> sum(x -> !isnan(x) && x < 0, q)) => :negative_ann
    )

    result = generate_stats(annual_data; value_cols=["negative_ann"],
        trend_completeness=trend_completeness, decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector)

    return result
end

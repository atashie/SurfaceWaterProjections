"""
Recession analysis signatures.

Power-law recession curve fitting and seasonality analysis.
"""

# Configuration loaded from config.jl (CFG_RECESSION_MIN_LENGTH, CFG_RECESSION_MIN_EVENTS)

"""
    ols_slope_intercept(x::Vector, y::Vector) -> (slope, intercept)

Simple OLS linear regression (like scipy.stats.linregress).
"""
function ols_slope_intercept(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    if n < 2
        return (NaN, NaN)
    end
    x_mean = mean(x)
    y_mean = mean(y)

    numerator = sum((x .- x_mean) .* (y .- y_mean))
    denominator = sum((x .- x_mean).^2)

    if abs(denominator) < 1e-10
        return (NaN, NaN)
    end

    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    return (slope, intercept)
end


"""
    identify_recession_events(Q::Vector, dowy::Vector; min_length=5) -> Vector{Tuple}

Identify recession events where both Q and |dQ/dt| are decreasing.

Uses "look-ahead" algorithm matching R/Python implementation:
At each position i, checks if the NEXT min_length consecutive days ALL meet criteria:
1. Q[i+j+1] < Q[i+j] (flow strictly decreasing)
2. |dQ/dt[i+j+1]| < |dQ/dt[i+j]| (rate becoming less negative)

Returns vector of (start_idx, end_idx) tuples.
"""
function identify_recession_events(
    Q::AbstractVector{<:Real},
    dowy::AbstractVector{<:Integer};
    min_length::Int=CFG_RECESSION_MIN_LENGTH
)
    n = length(Q)
    events = Tuple{Int, Int}[]

    if n < min_length + 1
        return events
    end

    # Calculate dQ/dt (forward difference) with NA at end
    dQdt = vcat(diff(Q), NaN)  # dQdt[i] = Q[i+1] - Q[i]

    in_recession = false
    start_idx = 1

    # Loop from 1 to n-min_length (matching R/Python)
    @inbounds for i in 1:(n - min_length)
        # Check if we have valid data at position i
        if isnan(Q[i]) || isnan(dQdt[i])
            if in_recession
                # End current recession if we hit NA
                if (i - start_idx) >= min_length
                    push!(events, (start_idx, i - 1))
                end
                in_recession = false
            end
            continue
        end

        # Check for monotonic decrease in both Q and |dQ/dt| for next min_length days
        if i < n - 1
            is_recession = true

            # Check next min_length consecutive decreases (requires min_length+1 points)
            for j in 0:(min_length - 1)
                idx = i + j
                if idx + 1 > n || isnan(Q[idx]) || isnan(Q[idx + 1]) ||
                   isnan(dQdt[idx]) || isnan(dQdt[idx + 1])
                    is_recession = false
                    break
                end

                # Check if Q is decreasing
                if Q[idx + 1] >= Q[idx]
                    is_recession = false
                    break
                end

                # Check if |dQ/dt| is decreasing (becoming less negative)
                if abs(dQdt[idx + 1]) >= abs(dQdt[idx])
                    is_recession = false
                    break
                end
            end

            if is_recession && !in_recession
                # Start new recession
                in_recession = true
                start_idx = i
            elseif !is_recession && in_recession
                # End current recession
                if (i - start_idx) >= min_length
                    push!(events, (start_idx, i - 1))
                end
                in_recession = false
            end
        end
    end

    # Check if we ended in a recession
    if in_recession && (n - start_idx) >= min_length
        push!(events, (start_idx, n))
    end

    return events
end


"""
    fit_recession_power_law(Q::Vector; min_points=2, use_ols=true, remove_first_day=false) -> (log_a, b)

Fit power-law recession model: dQ/dt = a * Q^b

Uses point cloud method with OLS (default, matches R/Python) or Theil-Sen regression.

Parameters
----------
Q : Vector
    Discharge values for the recession event
min_points : Int, default 2
    Minimum number of valid points required (R/Python use 2)
use_ols : Bool, default true
    Use OLS regression (true, matches R/Python) or Theil-Sen (false)
remove_first_day : Bool, default false
    Whether to remove first day (often influenced by peak)
"""
function fit_recession_power_law(
    Q::AbstractVector{<:Real};
    min_points::Int=2,
    use_ols::Bool=true,
    remove_first_day::Bool=false
)
    # Remove first day if requested (for per-event fitting, like R/Python)
    Q_work = remove_first_day && length(Q) > 1 ? Q[2:end] : Q

    if length(Q_work) < 3
        return (NaN, NaN)
    end

    # Calculate -dQ/dt
    dQdt = -diff(Q_work)

    # Use Q values at start of each interval (like R/Python)
    Q_mid = Q_work[1:end-1]

    # Filter valid points
    valid_mask = (Q_mid .> 0) .& (dQdt .> 0) .& .!isnan.(Q_mid) .& .!isnan.(dQdt)
    Q_valid = Q_mid[valid_mask]
    dQdt_valid = dQdt[valid_mask]

    if length(Q_valid) < min_points
        return (NaN, NaN)
    end

    # Log transform (natural log to match R/Python)
    log_Q = log.(Q_valid)
    log_dQdt = log.(dQdt_valid)

    # Fit: log(dQ/dt) = log(a) + b * log(Q)
    if use_ols
        # OLS regression (matches Python's scipy.stats.linregress)
        slope, intercept = ols_slope_intercept(log_Q, log_dQdt)
    else
        # Theil-Sen for robustness
        slope, intercept = theil_sen_slope(log_Q, log_dQdt)
    end

    return (intercept, slope)  # log_a, b
end


"""
    fit_recession_seasonality(log_a_values::Vector, dowy_values::Vector) -> Dict

Fit sinusoidal seasonality model to recession parameter.

Returns amplitude and day of minimum.
"""
function fit_recession_seasonality(
    log_a_values::AbstractVector{<:Real},
    dowy_values::AbstractVector{<:Integer}
)
    # Remove NaN
    valid_mask = .!isnan.(log_a_values)
    log_a = log_a_values[valid_mask]
    dowy = dowy_values[valid_mask]

    if length(log_a) < 10
        return Dict("amplitude" => NaN, "minimum_day" => NaN)
    end

    # Convert dowy to radians (match R's 365, not 365.25)
    omega = 2 * pi / 365
    theta = omega .* dowy

    # Fit: log_a = B1*sin(theta) + B2*cos(theta) + C
    # Using least squares (same basis order as R/Python: sin, cos, intercept)
    n = length(log_a)
    X = hcat(sin.(theta), cos.(theta), ones(n))
    y = log_a

    # Solve normal equations (with error handling matching R's tryCatch / Python's try/except)
    local B1, B2
    try
        coeffs = X \ y
        B1, B2 = coeffs[1], coeffs[2]  # sin coeff, cos coeff
    catch
        return Dict("amplitude" => NaN, "minimum_day" => NaN)
    end

    # Calculate amplitude
    amplitude = sqrt(B1^2 + B2^2)

    # Calculate phase (matching R's atan2(-B2, B1))
    phase_rad = atan(-B2, B1)
    phase_days = phase_rad * 365 / (2 * pi)

    # Ensure phase is between 0 and 365
    if phase_days < 0
        phase_days += 365
    end

    # Minimum occurs at phase + 273.75 days (3/4 of a cycle, matching R)
    minimum_day = phase_days + 273.75
    if minimum_day > 365
        minimum_day -= 365
    end

    return Dict("amplitude" => amplitude, "minimum_day" => minimum_day)
end


"""
    analyze_recession_parameters(df::DataFrame; min_events=25) -> Dict

Calculate recession curve signature trends.

Calculates metrics:
- log_a_pointcloud, log_a_events: Recession rate parameter
- b_pointcloud, b_events: Recession exponent
- concavity: Difference in b between recession halves
- log_a_seasonality_*: Sinusoidal seasonality parameters

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q, dowy
min_events : Int
    Minimum recession events required

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics
"""
function analyze_recession_parameters(df::DataFrame; min_events::Int=CFG_RECESSION_MIN_EVENTS)
    result = Dict{String, Float64}()

    # Define all expected output keys
    base_metrics = ["log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity"]
    seasonality_metrics = [
        "log_a_seasonality_amplitude_all", "log_a_seasonality_minimum_all",
        "log_a_seasonality_amplitude_first_half", "log_a_seasonality_minimum_first_half",
        "log_a_seasonality_amplitude_last_half", "log_a_seasonality_minimum_last_half"
    ]

    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year", "dowy"])
    if !valid
        @warn "analyze_recession_parameters: Missing columns: $missing_cols"
        for m in base_metrics
            merge!(result, empty_stats(m))
        end
        for m in seasonality_metrics
            result[m] = NaN
        end
        return result
    end

    if nrow(df) == 0
        for m in base_metrics
            merge!(result, empty_stats(m))
        end
        for m in seasonality_metrics
            result[m] = NaN
        end
        return result
    end

    # Collect all recession events across years (for seasonality analysis)
    all_log_a = Float64[]
    all_b = Float64[]
    all_dowy = Int[]
    all_event_years = Int[]  # Track water year for each event (for proper half-split)

    years = unique(df.water_year)

    # Pre-populate annual_data with ALL years (matching R/Python architecture)
    # This ensures years with pointcloud data but no events still contribute to trends
    annual_data = DataFrame(
        water_year = years,
        log_a_events = fill(NaN, length(years)),
        b_events = fill(NaN, length(years)),
        log_a_pointcloud = fill(NaN, length(years)),
        b_pointcloud = fill(NaN, length(years)),
        concavity = fill(NaN, length(years))
    )

    for yr in years
        year_mask = df.water_year .== yr
        if !any(skipmissing(year_mask))
            continue
        end
        year_df = df[coalesce.(year_mask, false), :]
        year_df = sort(year_df, :dowy)

        # Convert Q and dowy to proper types for helper functions
        Q_clean = coalesce_q(year_df.Q)
        dowy_clean = Int[coalesce(d, 0) for d in year_df.dowy]

        # Identify recession events
        events = identify_recession_events(Q_clean, dowy_clean)

        if isempty(events)
            continue
        end

        year_log_a = Float64[]
        year_b = Float64[]
        year_concavity = Float64[]
        successful_events_Q = Vector{Vector{Float64}}()  # Store Q data for log_a recalc

        # Collect point cloud data for this year (like R/Python)
        year_pc_Q = Float64[]
        year_pc_dQdt = Float64[]

        for (start_idx, end_idx) in events
            Q_event = Q_clean[start_idx:end_idx]
            dowy_mid = dowy_clean[div(start_idx + end_idx, 2)]

            # Per-event fitting with remove_first_day=true (matching R/Python)
            log_a, b = fit_recession_power_law(Q_event; remove_first_day=true)

            if !isnan(log_a) && !isnan(b)
                push!(year_log_a, log_a)
                push!(year_b, b)
                push!(successful_events_Q, copy(Q_event))
                push!(all_log_a, log_a)
                push!(all_b, b)
                push!(all_dowy, dowy_mid)
                push!(all_event_years, yr)

                # Calculate concavity for this event by splitting into halves (like R/Python)
                n_event = length(Q_event)
                if n_event >= 6  # Need minimum points for two half-fits
                    mid_pt = div(n_event, 2)
                    # Overlapping split matching R: first=[1:mid_pt], second=[mid_pt:end]
                    _, b_first = fit_recession_power_law(Q_event[1:mid_pt]; min_points=2, use_ols=true)
                    _, b_second = fit_recession_power_law(Q_event[mid_pt:end]; min_points=2, use_ols=true)
                    if !isnan(b_first) && !isnan(b_second)
                        push!(year_concavity, b_second - b_first)
                    end
                end

                # Add to point cloud data (removing first day, like R/Python)
                if length(Q_event) > 1
                    Q_subset = Q_event[2:end]  # Remove first day
                    dQ_subset = -diff(Q_event[2:end])
                    Q_for_pc = Q_subset[1:end-1]  # Q values for dQ/dt pairs

                    for j in 1:length(dQ_subset)
                        if Q_for_pc[j] > 0 && dQ_subset[j] > 0
                            push!(year_pc_Q, Q_for_pc[j])
                            push!(year_pc_dQdt, dQ_subset[j])
                        end
                    end
                end
            end
        end

        # Find row index for this year
        yr_idx = findfirst(annual_data.water_year .== yr)
        yr_idx === nothing && continue

        # Per-year point cloud analysis — assign independently of events (matching R/Python)
        if length(year_pc_Q) > 10
            log_Q = log.(year_pc_Q)
            log_dQdt = log.(year_pc_dQdt)

            # Skip near-singular data (matching R's lm() QR rank check)
            # R's tolerance is .Machine$double.eps^0.5 ≈ 1.49e-8
            if var(log_Q) >= 1e-8
                # Fit using OLS (lm() in R, linregress in Python)
                b_pc, _ = ols_slope_intercept(log_Q, log_dQdt)

                if !isnan(b_pc)
                    log_a_values_pc = log_dQdt .- b_pc .* log_Q
                    log_a_pc = median(log_a_values_pc)

                    annual_data[yr_idx, :log_a_pointcloud] = log_a_pc
                    annual_data[yr_idx, :b_pointcloud] = b_pc
                end
            end
        end

        # Event-based metrics
        if !isempty(year_b)
            # Recalculate log_a_events using median_b (matching R)
            median_b_val = median(year_b)
            log_a_events_recalc = Float64[]
            for Q_ev in successful_events_Q
                if length(Q_ev) > 1
                    Q_subset = Q_ev[2:end]  # Remove first day
                    dQ_subset = -diff(Q_ev[2:end])
                    Q_for_calc = Q_subset[1:end-1]
                    valid_mask = (Q_for_calc .> 0) .& (dQ_subset .> 0)
                    if any(valid_mask)
                        log_a_vals = log.(dQ_subset[valid_mask]) .- median_b_val .* log.(Q_for_calc[valid_mask])
                        push!(log_a_events_recalc, median(log_a_vals))
                    end
                end
            end

            annual_data[yr_idx, :log_a_events] = isempty(log_a_events_recalc) ? NaN : median(log_a_events_recalc)
            annual_data[yr_idx, :b_events] = median(year_b)
            annual_data[yr_idx, :concavity] = isempty(year_concavity) ? NaN : mean(year_concavity)
        end
    end

    # Check minimum events requirement (matching R/Python)
    # total events = length(all_log_a) = number of successfully-fitted events
    if length(all_log_a) < min_events
        # Not enough recession events — return all NAs
        for m in base_metrics
            merge!(result, empty_stats(m))
        end
        for m in seasonality_metrics
            result[m] = NaN
        end
        return result
    end

    # Generate statistics for annual metrics
    if nrow(annual_data) >= 3
        for m in base_metrics
            if m in names(annual_data)
                stats = generate_stats(annual_data; value_cols=[m])
                merge!(result, stats)
            else
                merge!(result, empty_stats(m))
            end
        end
    else
        for m in base_metrics
            merge!(result, empty_stats(m))
        end
    end

    # Seasonality analysis (static values, not trends)
    if length(all_log_a) >= min_events
        # All data
        seasonality_all = fit_recession_seasonality(all_log_a, all_dowy)
        result["log_a_seasonality_amplitude_all"] = seasonality_all["amplitude"]
        result["log_a_seasonality_minimum_all"] = seasonality_all["minimum_day"]

        # Split into first and last half of years (like R/Python)
        median_water_year = median(unique(all_event_years))
        first_half_mask = all_event_years .<= median_water_year
        last_half_mask = all_event_years .> median_water_year

        # First half
        if sum(first_half_mask) >= 10
            seasonality_first = fit_recession_seasonality(
                all_log_a[first_half_mask], all_dowy[first_half_mask]
            )
            result["log_a_seasonality_amplitude_first_half"] = seasonality_first["amplitude"]
            result["log_a_seasonality_minimum_first_half"] = seasonality_first["minimum_day"]
        else
            result["log_a_seasonality_amplitude_first_half"] = NaN
            result["log_a_seasonality_minimum_first_half"] = NaN
        end

        # Last half
        if sum(last_half_mask) >= 10
            seasonality_last = fit_recession_seasonality(
                all_log_a[last_half_mask], all_dowy[last_half_mask]
            )
            result["log_a_seasonality_amplitude_last_half"] = seasonality_last["amplitude"]
            result["log_a_seasonality_minimum_last_half"] = seasonality_last["minimum_day"]
        else
            result["log_a_seasonality_amplitude_last_half"] = NaN
            result["log_a_seasonality_minimum_last_half"] = NaN
        end
    else
        for m in seasonality_metrics
            result[m] = NaN
        end
    end

    return result
end

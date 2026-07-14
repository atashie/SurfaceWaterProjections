"""
Recession analysis signatures.

Power-law recession curve fitting and seasonality analysis.

Parameter conventions (July 2026): the recession exponent b and concavity are
estimated with FREE power-law fits (dQ/dt = -a*Q^b), but all alpha outputs
(log_a_pointcloud, log_a_events, log_a seasonality, alpha_linear, and the
recession_alpha scalar) assume a LINEAR reservoir (b fixed at 1) across all
locations and periods. Rationale: log(a) is the intercept of a regression whose
slope is b, so free-fit alpha estimates are convolved with b; fixing b = 1
decouples them.
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

Uses position-level marking: marks each position where both recession criteria hold,
then finds contiguous runs of valid positions >= min_length. This captures the full
extent of each recession event, unlike the previous look-ahead algorithm which could
truncate events by up to min_length days from their true end.

Returns vector of (start_idx, end_idx) INCLUSIVE-INCLUSIVE tuples, where
Q[start_idx:end_idx] gives all days in the recession event.
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

    # Calculate dQ/dt (forward difference) with NaN at end
    dQdt = vcat(diff(Q), NaN)  # dQdt[i] = Q[i+1] - Q[i]

    # Step 1: Mark each position where both recession criteria hold
    # Position i is "valid" if Q[i+1] < Q[i] AND |dQdt[i+1]| < |dQdt[i]|
    is_valid = falses(n)
    @inbounds for i in 1:(n - 1)
        if !isnan(Q[i]) && !isnan(Q[i + 1]) && !isnan(dQdt[i]) && !isnan(dQdt[i + 1])
            if Q[i + 1] < Q[i] && abs(dQdt[i + 1]) < abs(dQdt[i])
                is_valid[i] = true
            end
        end
    end

    # Step 2: Find connected runs of valid positions >= min_length
    # A run of k valid transitions at positions start to start+k-1 covers
    # days start to start+k (k+1 days). Return INCLUSIVE-INCLUSIVE tuples.
    run_start = 0
    @inbounds for i in 1:n
        if is_valid[i]
            if run_start == 0
                run_start = i
            end
        else
            if run_start > 0
                run_length = i - run_start  # number of valid transitions
                if run_length >= min_length
                    # Last valid position is i-1; event spans start through i
                    # (i = first day after last valid transition, included in Q slice)
                    push!(events, (run_start, i))
                end
                run_start = 0
            end
        end
    end
    # Handle run ending at data boundary
    if run_start > 0
        run_length = n - run_start
        if run_length >= min_length
            push!(events, (run_start, n))
        end
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
    event_log_a_b1(Q_event::Vector) -> Float64

Per-event recession rate parameter log(a) under a fixed linear-reservoir
assumption (b = 1): with -dQ/dt = a*Q, log(a) = log(-dQ/dt) - log(Q), so the
event value is the median of that quantity over the event's valid pairs.

Fixing b = 1 removes the slope-intercept convolution of the free power-law fit
(log(a) is the intercept of a regression whose slope is b, so free-fit alpha
estimates absorb any variation in b). b itself is still reported from the free
fit — see fit_recession_power_law.

Conventions match the power-law fitting path: first day removed (storm peak),
pairs require Q > 0 and -dQ > 0, natural log. Returns NaN if no valid pairs.
"""
function event_log_a_b1(Q_event::AbstractVector{<:Real})
    if length(Q_event) < 3
        return NaN
    end
    Q_work = Q_event[2:end]        # Remove first day (storm peak)
    dQdt = -diff(Q_work)
    Q_mid = Q_work[1:end-1]        # Q at start of each interval

    valid_mask = (Q_mid .> 0) .& (dQdt .> 0) .& .!isnan.(Q_mid) .& .!isnan.(dQdt)
    if !any(valid_mask)
        return NaN
    end
    return median(log.(dQdt[valid_mask]) .- log.(Q_mid[valid_mask]))
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
- log_a_pointcloud, log_a_events: Recession rate parameter under a FIXED
  linear-reservoir assumption (b = 1): log(a) = median(log(-dQ/dt) - log(Q))
- b_pointcloud, b_events: Recession exponent (free power-law fit)
- concavity: Difference in b between recession halves (free fits)
- log_a_seasonality_*: Sinusoidal seasonality of the per-event b=1 log_a values

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
function analyze_recession_parameters(df::DataFrame; min_events::Int=CFG_RECESSION_MIN_EVENTS, trend_completeness::Union{Nothing, Float64}=nothing, decade_completeness::Union{Nothing, Float64}=nothing, changepoint::Union{Nothing, NamedTuple}=nothing, collector::Union{Nothing, AnnualCollector}=nothing)
    result = Dict{String, Float64}()

    # Define all expected output keys
    base_metrics = ["log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events", "concavity", "alpha_linear"]
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
        merge!(result, empty_stats("n_recession_events"))
        for m in seasonality_metrics
            result[m] = NaN
        end
        result["recession_alpha_point_cloud_linear_reservoir"] = NaN
        return result
    end

    if nrow(df) == 0
        for m in base_metrics
            merge!(result, empty_stats(m))
        end
        merge!(result, empty_stats("n_recession_events"))
        for m in seasonality_metrics
            result[m] = NaN
        end
        result["recession_alpha_point_cloud_linear_reservoir"] = NaN
        return result
    end

    # Collect all recession events across years (for seasonality analysis)
    all_log_a = Float64[]
    all_b = Float64[]
    all_dowy = Int[]
    all_event_years = Int[]  # Track water year for each event (for proper half-split)
    all_alpha_linear = Float64[]  # Whole-record alpha values for scalar

    years = unique(df.water_year)

    # Pre-populate annual_data with ALL years (matching R/Python architecture)
    # This ensures years with pointcloud data but no events still contribute to trends
    annual_data = DataFrame(
        water_year = years,
        log_a_events = fill(NaN, length(years)),
        b_events = fill(NaN, length(years)),
        log_a_pointcloud = fill(NaN, length(years)),
        b_pointcloud = fill(NaN, length(years)),
        concavity = fill(NaN, length(years)),
        n_recession_events = fill(NaN, length(years)),
        alpha_linear = fill(NaN, length(years))
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

        # Store event count for this year (INDEPENDENT of min_events gate)
        yr_idx_count = findfirst(annual_data.water_year .== yr)
        if yr_idx_count !== nothing
            annual_data[yr_idx_count, :n_recession_events] = Float64(length(events))
        end

        if isempty(events)
            continue
        end

        year_log_a = Float64[]     # per-event log_a under fixed b = 1 (linear reservoir)
        year_b = Float64[]         # per-event b from the free power-law fit
        year_concavity = Float64[]
        year_alpha_linear = Float64[]

        # Collect point cloud data for this year (like R/Python)
        year_pc_Q = Float64[]
        year_pc_dQdt = Float64[]

        for (start_idx, end_idx) in events
            Q_event = Q_clean[start_idx:end_idx]
            dowy_mid = dowy_clean[div(start_idx + end_idx, 2)]

            # Per-event FREE power-law fit (remove_first_day=true) — used for b only.
            # Fit success also gates event acceptance (unchanged event counts).
            log_a, b = fit_recession_power_law(Q_event; remove_first_day=true)

            # Per-event alpha under fixed b = 1 (linear reservoir). log(a) is the
            # intercept of the free fit and is convolved with its slope b; fixing
            # b = 1 across all locations and periods decouples the two. Same
            # valid-pair definition as the free fit, so non-NaN whenever the fit
            # succeeds.
            log_a_b1 = event_log_a_b1(Q_event)

            # Compute discrete recession constant alpha = Q_{i+1}/Q_i (b=1 linear reservoir)
            # INDEPENDENT of power-law fit — depends only on raw Q pairs
            # Remove first day (storm peak) consistent with point-cloud fitting
            if length(Q_event) > 2  # Need at least 3 days (2 after removing first)
                Q_alpha = Q_event[2:end]  # Remove first day
                for j in 1:(length(Q_alpha) - 1)
                    if Q_alpha[j] > 0 && !isnan(Q_alpha[j]) && !isnan(Q_alpha[j + 1])
                        a_i = Q_alpha[j + 1] / Q_alpha[j]
                        if 0 < a_i < 1  # Must be valid recession (decreasing Q)
                            push!(year_alpha_linear, a_i)
                            push!(all_alpha_linear, a_i)
                        end
                    end
                end
            end

            if !isnan(log_a) && !isnan(b)
                push!(year_log_a, log_a_b1)
                push!(year_b, b)
                push!(all_log_a, log_a_b1)
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

        # Per-year alpha_linear (point cloud with b=1, linear reservoir)
        if length(year_alpha_linear) > 10
            annual_data[yr_idx, :alpha_linear] = median(year_alpha_linear)
        end

        # Per-year point cloud analysis — assign independently of events (matching R/Python)
        if length(year_pc_Q) > 10
            log_Q = log.(year_pc_Q)
            log_dQdt = log.(year_pc_dQdt)

            # log_a under fixed b = 1 (linear reservoir): log(a) = log(-dQ/dt) - log(Q).
            # No regression is involved, so this needs no singularity gate and is
            # independent of the free b fit below.
            annual_data[yr_idx, :log_a_pointcloud] = median(log_dQdt .- log_Q)

            # b from the FREE point-cloud fit (unchanged).
            # Skip near-singular data (matching R's lm() QR rank check)
            # R's tolerance is .Machine$double.eps^0.5 ≈ 1.49e-8
            if var(log_Q) >= 1e-8
                # Fit using OLS (lm() in R, linregress in Python)
                b_pc, _ = ols_slope_intercept(log_Q, log_dQdt)

                if !isnan(b_pc)
                    annual_data[yr_idx, :b_pointcloud] = b_pc
                end
            end
        end

        # Event-based metrics. log_a_events uses the per-event fixed-b=1 values
        # (year_log_a) — the previous median-b recalculation is gone: alpha is
        # decoupled from b by assuming a linear reservoir everywhere.
        if !isempty(year_b)
            year_log_a_valid = filter(!isnan, year_log_a)
            annual_data[yr_idx, :log_a_events] = isempty(year_log_a_valid) ? NaN : median(year_log_a_valid)
            annual_data[yr_idx, :b_events] = median(year_b)
            annual_data[yr_idx, :concavity] = isempty(year_concavity) ? NaN : mean(year_concavity)
        end
    end

    # n_recession_events stats computed INDEPENDENTLY of min_events gate
    # Event counts are meaningful even for gages with few events
    events_df = annual_data[.!isnan.(annual_data.n_recession_events), [:water_year, :n_recession_events]]
    if nrow(events_df) >= 3
        merge!(result, generate_stats(events_df; value_cols=["n_recession_events"],
            trend_completeness=trend_completeness, decade_completeness=decade_completeness, changepoint=changepoint, collector=collector))
    else
        merge!(result, empty_stats("n_recession_events"))
    end

    # Whole-record recession alpha (scalar — independent of min_events gate)
    if length(all_alpha_linear) > 10
        result["recession_alpha_point_cloud_linear_reservoir"] = median(all_alpha_linear)
    else
        result["recession_alpha_point_cloud_linear_reservoir"] = NaN
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
                stats = generate_stats(annual_data; value_cols=[m], trend_completeness=trend_completeness, decade_completeness=decade_completeness, changepoint=changepoint, collector=collector)
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

"""
Flow duration curve signatures.

Slopes of the flow duration curve at different exceedance ranges.
"""

"""
    calculate_fdc_slope(Q::Vector, exceedance_range::Tuple) -> Float64

Calculate slope of FDC in log-space for a given exceedance range.

Parameters
----------
Q : Vector
    Daily discharge values
exceedance_range : Tuple{Float64, Float64}
    (min_exceedance, max_exceedance) probabilities in [0, 1]

Returns
-------
Float64
    Slope of FDC in log-space
"""
function calculate_fdc_slope(Q::AbstractVector{<:Real}, exceedance_range::Tuple{<:Real, <:Real})
    # Remove NaN and negative values (keep zeros, matching R/Python)
    Q_valid = filter(x -> !isnan(x) && x >= 0, Q)
    n = length(Q_valid)

    if n < 10
        return NaN
    end

    # Sort descending (highest flow first)
    Q_sorted = sort(Q_valid, rev=true)

    # Calculate exceedance probabilities (Weibull plotting position, matching R)
    exceedance = [i / (n + 1) for i in 1:n]

    # Filter to range
    min_exc, max_exc = exceedance_range
    mask = (exceedance .>= min_exc) .& (exceedance .<= max_exc)

    exc_range = exceedance[mask]
    Q_range = Q_sorted[mask]

    if length(exc_range) < 3
        return NaN
    end

    # Log-transform Q (add small constant to handle zeros, matching R/Python)
    log_Q = log10.(Q_range .+ 1e-10)

    # Linear regression slope
    slope = linear_slope(exc_range, log_Q)

    return slope
end


"""
    analyze_fdc_trends(df::DataFrame; min_days=CFG_FDC_MIN_DAYS) -> Dict

Calculate flow duration curve signature trends.

Calculates 3 FDC metrics:
- FDCall: Slope across full range (0-100%)
- FDC90th: Slope at low flows (90-100%)
- FDCmid: Slope at mid-range (20-80%)

Each metric produces 8 statistics via generate_stats().

Parameters
----------
df : DataFrame
    Daily streamflow data with columns: water_year, Q
min_days : Int
    Minimum days per year for valid calculation

Returns
-------
Dict{String, Float64}
    Dictionary of signature statistics (3 metrics × 8 stats = 24 values)
"""
function analyze_fdc_trends(df::DataFrame; min_days::Int=CFG_FDC_MIN_DAYS)
    result = Dict{String, Float64}()
    metrics = ["FDCall", "FDC90th", "FDCmid"]

    # Column validation
    valid, missing_cols = validate_columns(df, ["Q", "water_year"])
    if !valid
        @warn "analyze_fdc_trends: Missing columns: $missing_cols"
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

    # Convert Q to Float64, handling missing values
    Q_clean = coalesce_q(df.Q)

    # FDC ranges (exceedance probabilities, matching R/Python)
    fdc_ranges = Dict(
        "FDCall" => (0.0, 1.0),
        "FDC90th" => (0.9, 1.0),
        "FDCmid" => (0.2, 0.8)
    )

    years = unique(df.water_year)
    annual_data = DataFrame(
        water_year = years,
        FDCall = fill(NaN, length(years)),
        FDC90th = fill(NaN, length(years)),
        FDCmid = fill(NaN, length(years))
    )

    for yr in years
        year_mask = coalesce.(df.water_year .== yr, false)
        if !any(year_mask)
            continue
        end
        year_Q = Q_clean[year_mask]

        # Check minimum data
        n_valid = sum(.!isnan.(year_Q))
        if n_valid < min_days
            continue
        end

        yr_idx = findfirst(annual_data.water_year .== yr)
        yr_idx === nothing && continue

        for (metric, range) in fdc_ranges
            annual_data[yr_idx, Symbol(metric)] = calculate_fdc_slope(year_Q, range)
        end
    end

    if nrow(annual_data) < 3
        for m in metrics
            merge!(result, empty_stats(m))
        end
        return result
    end

    # Generate statistics for all metrics
    result = generate_stats(annual_data; value_cols=metrics)

    return result
end

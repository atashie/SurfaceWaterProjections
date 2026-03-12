"""
Statistical functions for streamflow signature analysis.

Core functions for trend analysis and summary statistics.
"""

using Distributions: Normal, TDist, cdf

# Stat suffixes matching R implementation
const STAT_SUFFIXES = [
    "_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
    "_mk_rho", "_mk_pval", "_mean", "_median"
]

"""
    theil_sen_slope(x::Vector, y::Vector) -> (slope, intercept)

Calculate Theil-Sen robust slope estimator.

Returns median of all pairwise slopes and corresponding intercept.
"""
function theil_sen_slope(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    if n != length(y)
        error("x and y must have same length")
    end
    if n < 2
        return (NaN, NaN)
    end

    # Remove NaN pairs
    valid_mask = .!isnan.(x) .& .!isnan.(y)
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]
    n_valid = length(x_valid)

    if n_valid < 2
        return (NaN, NaN)
    end

    # PERFORMANCE FIX: Preallocate slopes array instead of using push!()
    # Maximum possible pairs = n_valid * (n_valid - 1) / 2
    max_pairs = div(n_valid * (n_valid - 1), 2)
    slopes = Vector{Float64}(undef, max_pairs)
    n_slopes = 0

    @inbounds for i in 1:n_valid
        for j in (i+1):n_valid
            dx = x_valid[j] - x_valid[i]
            if abs(dx) > 1e-10
                n_slopes += 1
                slopes[n_slopes] = (y_valid[j] - y_valid[i]) / dx
            end
        end
    end

    if n_slopes == 0
        return (NaN, NaN)
    end

    # Resize to actual number of valid slopes
    resize!(slopes, n_slopes)

    slope = median(slopes)
    intercept = median(y_valid .- slope .* x_valid)

    return (slope, intercept)
end


"""
    mann_kendall_test(y::Vector) -> (tau, pvalue)

Perform Mann-Kendall trend test.

Returns Kendall's tau and two-sided p-value.
"""
function mann_kendall_test(y::AbstractVector{<:Real})
    # Remove NaN values
    y_valid = filter(!isnan, y)
    n = length(y_valid)

    if n < 3
        return (NaN, NaN)
    end

    # Calculate S statistic - PERFORMANCE FIX: Added @inbounds
    S = 0
    @inbounds for i in 1:n
        for j in (i+1):n
            diff = y_valid[j] - y_valid[i]
            if diff > 0
                S += 1
            elseif diff < 0
                S -= 1
            end
        end
    end

    # Calculate variance
    # Account for ties
    unique_vals = unique(y_valid)
    tie_groups = [count(==(v), y_valid) for v in unique_vals]
    tie_correction = sum((t * (t - 1) * (2t + 5) for t in tie_groups if t > 1); init=0)

    var_S = (n * (n - 1) * (2n + 5) - tie_correction) / 18

    if var_S <= 0
        return (NaN, NaN)
    end

    # Calculate tau
    n_pairs = n * (n - 1) / 2
    tau = S / n_pairs

    # Calculate z-score and p-value
    if S > 0
        z = (S - 1) / sqrt(var_S)
    elseif S < 0
        z = (S + 1) / sqrt(var_S)
    else
        z = 0.0
    end

    # Two-sided p-value using normal approximation
    pvalue = 2 * (1 - cdf(Normal(), abs(z)))

    return (tau, pvalue)
end



"""
    spearman_correlation(x::Vector, y::Vector) -> (rho, pvalue)

Calculate Spearman rank correlation coefficient and p-value.
"""
function spearman_correlation(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    # Remove NaN pairs
    valid_mask = .!isnan.(x) .& .!isnan.(y)
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]
    n = length(x_valid)

    if n < 3
        return (NaN, NaN)
    end

    # Calculate ranks (average/fractional ranks for ties, matching R and Python)
    rank_x = tiedrank(x_valid)
    rank_y = tiedrank(y_valid)

    # Spearman rho
    rho = cor(rank_x, rank_y)

    if isnan(rho)
        return (NaN, NaN)
    end

    # t-statistic for significance
    if abs(rho) >= 1.0
        pvalue = 0.0
    else
        t_stat = rho * sqrt((n - 2) / (1 - rho^2))
        # Two-sided p-value from t-distribution
        pvalue = 2 * (1 - cdf(TDist(n - 2), abs(t_stat)))
    end

    return (rho, pvalue)
end



"""
    linear_slope(x::Vector, y::Vector) -> slope

Calculate ordinary least squares linear regression slope.
"""
function linear_slope(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    # Remove NaN pairs
    valid_mask = .!isnan.(x) .& .!isnan.(y)
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]
    n = length(x_valid)

    if n < 2
        return NaN
    end

    x_mean = mean(x_valid)
    y_mean = mean(y_valid)

    numerator = sum((x_valid .- x_mean) .* (y_valid .- y_mean))
    denominator = sum((x_valid .- x_mean).^2)

    if abs(denominator) < 1e-10
        return NaN
    end

    return numerator / denominator
end


"""
    generate_stats(df::DataFrame; value_cols=nothing, year_col="water_year", min_rows=3) -> Dict

Generate 8 summary statistics for each metric column.

For each column, calculates:
- Theil-Sen slope (robust trend)
- Linear slope (OLS trend)
- Spearman rho and p-value (correlation with time)
- Mann-Kendall tau and p-value (trend test)
- Mean and median

Parameters
----------
df : DataFrame
    Data with year column and value columns
value_cols : Vector{String} or nothing
    Columns to analyze (default: all numeric except year_col)
year_col : String
    Name of year column (default: "water_year")
min_rows : Int
    Minimum rows required for statistics

Returns
-------
Dict{String, Float64}
    Dictionary with keys like "metric_senn_slp", "metric_mean", etc.
"""
function generate_stats(
    df::DataFrame;
    value_cols::Union{Nothing, Vector{String}, Vector{Symbol}} = nothing,
    year_col::Union{String, Symbol} = "water_year",
    min_rows::Int = 3
)
    result = Dict{String, Float64}()

    if nrow(df) < min_rows
        # Return empty dict - caller should handle
        return result
    end

    # Sort by year
    df_sorted = sort(df, year_col)

    # Get year values as numeric
    years = Float64.(df_sorted[!, year_col])

    # Determine columns to process
    if value_cols === nothing
        # All numeric columns except year
        numeric_cols = [name for name in names(df_sorted)
                       if eltype(df_sorted[!, name]) <: Real &&
                          String(name) != String(year_col)]
        value_cols = String.(numeric_cols)
    else
        value_cols = String.(value_cols)
    end

    for col in value_cols
        if !(col in names(df_sorted))
            continue
        end

        values = Float64.(df_sorted[!, col])

        # Pre-filter years and values to exclude NaN rows (matching R/Python)
        # R: working_data <- working_data[!is.na(working_data$value), ]
        # Python: mask = ~data[col].isna(); years = data.loc[mask, ...]
        valid_mask = .!isnan.(values)
        valid_values = values[valid_mask]
        valid_years = years[valid_mask]

        if length(valid_values) < min_rows
            # Set all stats to NaN
            for suffix in STAT_SUFFIXES
                result["$(col)$(suffix)"] = NaN
            end
            continue
        end

        # Theil-Sen slope (use pre-filtered years)
        slope, _ = theil_sen_slope(valid_years, valid_values)
        result["$(col)_senn_slp"] = slope

        # Linear slope (use pre-filtered years)
        result["$(col)_linear_slp"] = linear_slope(valid_years, valid_values)

        # Spearman correlation (use pre-filtered years)
        rho, pval = spearman_correlation(valid_years, valid_values)
        result["$(col)_spearman_rho"] = rho
        result["$(col)_spearman_pval"] = pval

        # Mann-Kendall test (uses only values, already pre-filtered)
        tau, mk_pval = mann_kendall_test(valid_values)
        result["$(col)_mk_rho"] = tau
        result["$(col)_mk_pval"] = mk_pval

        # Mean and median
        result["$(col)_mean"] = mean(valid_values)
        result["$(col)_median"] = median(valid_values)
    end

    return result
end


"""
    empty_stats(metric::String) -> Dict

Return dict with NaN values for all 8 statistics for a metric.
"""
function empty_stats(metric::String)
    return Dict(
        "$(metric)_senn_slp" => NaN,
        "$(metric)_linear_slp" => NaN,
        "$(metric)_spearman_rho" => NaN,
        "$(metric)_spearman_pval" => NaN,
        "$(metric)_mk_rho" => NaN,
        "$(metric)_mk_pval" => NaN,
        "$(metric)_mean" => NaN,
        "$(metric)_median" => NaN
    )
end


"""
    validate_columns(df::DataFrame, required_cols::Vector{String}) -> (Bool, Vector{String})

Check if DataFrame has required columns.

Returns tuple of (all_present, missing_cols).
"""
function validate_columns(df::DataFrame, required_cols::Vector{String})
    df_cols = String.(names(df))
    missing_cols = [c for c in required_cols if !(c in df_cols)]
    return (isempty(missing_cols), missing_cols)
end


"""
    count_valid_q(Q::AbstractVector) -> Int

Count non-missing, non-NaN values in a Q vector.
Handles both Missing and NaN types.
"""
function count_valid_q(Q::AbstractVector)
    count = 0
    for val in Q
        if !ismissing(val) && (!(val isa Number) || !isnan(val))
            count += 1
        end
    end
    return count
end


"""
    coalesce_q(Q::AbstractVector) -> Vector{Float64}

Convert Q vector to Float64, replacing missing/NaN with NaN.
"""
function coalesce_q(Q::AbstractVector)
    result = Vector{Float64}(undef, length(Q))
    # PERFORMANCE FIX: Added @inbounds
    @inbounds for (i, val) in enumerate(Q)
        if ismissing(val)
            result[i] = NaN
        elseif val isa Number
            result[i] = Float64(val)
        else
            result[i] = NaN
        end
    end
    return result
end

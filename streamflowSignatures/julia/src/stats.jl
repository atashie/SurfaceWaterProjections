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

# Changepoint column suffixes (Pettitt test only)
const CP_SUFFIXES = [
    "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean", "_pettitt_post_mean",
    "_pettitt_delta_mean", "_pettitt_pct_change", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"
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

    # Calculate tau-b. All three languages agree here: R's Kendall::MannKendall
    # computes the same tie-corrected tau, and Python implements this formula
    # directly since 2026-08-26 (no scipy.stats.kendalltau).
    # For Mann-Kendall, x = time indices (no ties in x), so:
    #   tau_b = S / sqrt((n_pairs - T_y) * n_pairs)
    # where T_y = number of tied pairs in y
    n_pairs = n * (n - 1) / 2
    T_y = sum(t * (t - 1) / 2 for t in tie_groups; init=0.0)
    denom = sqrt((n_pairs - T_y) * n_pairs)
    tau = denom > 0.0 ? S / denom : NaN

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

    # The continuity-corrected normal approximation above is the shared
    # convention across all three languages: R's Kendall::MannKendall
    # reproduces this formula, and Python implements it directly (2026-08-26,
    # no scipy). For n<=3 R's MannKendall fails (IFAULT=12, returns p=1.0);
    # match that behavior.
    if n <= 3
        pvalue = 1.0
    end

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
    AnnualCollector()

Opt-in accumulator for per-year annual signature values in long format.

When passed to `generate_stats()` via the `collector` keyword, the exact annual
series handed to `generate_stats()` is appended as (signature, water_year, value)
triplets — before any min_rows or trend-completeness gating — so the collected
rows are exactly what the summary statistics were computed from. Collection is
read-only with respect to the statistics.

One instance is used per gage; the benchmark runner drains it into a global
long-format table tagged with gage_id.
"""
mutable struct AnnualCollector
    signature::Vector{String}
    water_year::Vector{Int32}
    value::Vector{Float64}
end

AnnualCollector() = AnnualCollector(String[], Int32[], Float64[])


"""
    _collect_annual!(collector, df, value_cols, year_col)

Append the annual series exactly as passed to `generate_stats()` (long format).
Rows whose year is `missing` or non-finite are skipped (they cannot be keyed);
`missing` values become NaN. NaN values are kept — they record "year present,
metric not computable".
"""
function _collect_annual!(
    collector::AnnualCollector,
    df::DataFrame,
    value_cols::Vector{String},
    year_col::Union{String, Symbol}
)
    years_raw = df[!, year_col]
    n = length(years_raw)

    # Pre-compute keyable years; skip missing / non-finite (Int32(round(NaN)) would throw)
    keyed = Vector{Int32}(undef, n)
    keep = falses(n)
    @inbounds for i in 1:n
        y = years_raw[i]
        y === missing && continue
        yf = Float64(y)
        isfinite(yf) || continue
        keyed[i] = Int32(round(yf))
        keep[i] = true
    end

    df_names = names(df)
    for col in value_cols
        col in df_names || continue
        vals = df[!, col]
        @inbounds for i in 1:n
            keep[i] || continue
            v = vals[i]
            push!(collector.signature, col)
            push!(collector.water_year, keyed[i])
            push!(collector.value, v === missing ? NaN : Float64(v))
        end
    end
    return nothing
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
min_values_for_stats : Int or nothing
    Optional floor (July 2026): a metric with fewer non-NaN annual values than this
    emits NaN for all 8 statistics and changepoint fields. Applied by the
    orchestrator to all families except recession and elasticity.
force_skip_trends : Set{String} or nothing
    Columns whose 6 trend statistics are forced to NaN regardless of their own
    completeness (mean/median and changepoint fields still computed). Used by
    caller-computed gates, e.g. the snow record-anchored decade gate (July 2026).

Returns
-------
Dict{String, Float64}
    Dictionary with keys like "metric_senn_slp", "metric_mean", etc.
"""
function generate_stats(
    df::DataFrame;
    value_cols::Union{Nothing, Vector{String}, Vector{Symbol}} = nothing,
    year_col::Union{String, Symbol} = "water_year",
    min_rows::Int = 3,
    trend_completeness::Union{Nothing, Float64} = nothing,
    decade_completeness::Union{Nothing, Float64} = nothing,
    changepoint::Union{Nothing, NamedTuple} = nothing,
    collector::Union{Nothing, AnnualCollector} = nothing,
    min_values_for_stats::Union{Nothing, Int} = nothing,
    force_skip_trends::Union{Nothing, Set{String}} = nothing
)
    result = Dict{String, Float64}()

    # Resolve value columns once, against the unsorted input df (previously this
    # identical predicate was duplicated in both branches below; eltypes are
    # unaffected by sorting, so hoisting preserves semantics exactly)
    if value_cols === nothing
        value_cols = String[name for name in names(df)
                            if eltype(df[!, name]) <: Real &&
                               String(name) != String(year_col)]
    else
        value_cols = String.(value_cols)
    end

    # Opt-in annual-values collection: capture the series exactly as passed,
    # BEFORE any gating below. Read-only with respect to the statistics.
    if collector !== nothing
        _collect_annual!(collector, df, value_cols, year_col)
    end

    if nrow(df) < min_rows
        for col in value_cols
            for suffix in STAT_SUFFIXES
                result["$(col)$(suffix)"] = NaN
            end
            if changepoint !== nothing
                for cp_suffix in CP_SUFFIXES
                    result["$(col)$(cp_suffix)"] = NaN
                end
            end
        end
        return result
    end

    # Sort by year
    df_sorted = sort(df, year_col)

    # Get year values as numeric
    years = Float64.(df_sorted[!, year_col])

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

        # Stats floor (July 2026): with min_values_for_stats set, a metric below the
        # floor emits NaN for ALL 8 statistics AND the changepoint fields — a 4-point
        # clustered series must not carry a Theil-Sen slope. Collection into the
        # annual-values collector happened above, BEFORE this gate (export contract).
        if length(valid_values) < min_rows ||
           (min_values_for_stats !== nothing && length(valid_values) < min_values_for_stats)
            # Set all stats to NaN
            for suffix in STAT_SUFFIXES
                result["$(col)$(suffix)"] = NaN
            end
            if changepoint !== nothing
                for cp_suffix in CP_SUFFIXES
                    result["$(col)$(cp_suffix)"] = NaN
                end
            end
            continue
        end

        # --- Per-column trend completeness check (opt-in) ---
        # Only active when caller explicitly passes trend_completeness.
        # Uses the metric's own non-NA year range for decade boundaries,
        # so climate signatures use climate years and recession uses recession years.
        compute_trends = true

        # Externally-forced trend skip (July 2026): caller-computed gates (e.g. the
        # snow record-anchored decade gate) suppress the 6 trend stats through the
        # SAME path as trend_completeness — mean/median and changepoint unaffected.
        if force_skip_trends !== nothing && col in force_skip_trends
            compute_trends = false
        end

        if compute_trends && trend_completeness !== nothing
            tc_frac = trend_completeness
            metric_yr_min = Int(minimum(valid_years))
            metric_yr_max = Int(maximum(valid_years))
            metric_span = metric_yr_max - metric_yr_min + 1
            metric_count = length(valid_years)

            # Overall completeness: non-NA values / span of metric's year range
            if metric_span > 0 && (metric_count / metric_span) < tc_frac
                compute_trends = false
            end

            # First/last decade check (only if span >= 10 years)
            if compute_trends && decade_completeness !== nothing && metric_span >= 10
                dc_frac = decade_completeness
                metric_yr_ints = Int.(valid_years)

                # First decade: metric_yr_min to metric_yr_min + 9
                first_expected = min(10, metric_span)
                first_actual = count(y -> metric_yr_min <= y <= metric_yr_min + 9, metric_yr_ints)
                if first_expected > 0 && (first_actual / first_expected) < dc_frac
                    compute_trends = false
                end

                # Last decade: metric_yr_max - 9 to metric_yr_max
                if compute_trends
                    last_start = metric_yr_max - 9
                    last_expected = min(10, metric_yr_max - max(last_start, metric_yr_min) + 1)
                    last_actual = count(y -> last_start <= y <= metric_yr_max, metric_yr_ints)
                    if last_expected > 0 && (last_actual / last_expected) < dc_frac
                        compute_trends = false
                    end
                end
            end

        end

        if !compute_trends
            # Insufficient coverage or forced skip: set 6 trend stats to NaN, still
            # compute mean/median. Changepoint runs independently of trend gating.
            result["$(col)_senn_slp"] = NaN
            result["$(col)_linear_slp"] = NaN
            result["$(col)_spearman_rho"] = NaN
            result["$(col)_spearman_pval"] = NaN
            result["$(col)_mk_rho"] = NaN
            result["$(col)_mk_pval"] = NaN
            result["$(col)_mean"] = mean(valid_values)
            result["$(col)_median"] = length(valid_values) > 0 ? median(valid_values) : NaN
            if changepoint !== nothing
                _run_changepoint_block!(result, col, valid_years, valid_values, changepoint)
            end
            continue
        end

        if compute_trends
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
        else
            # Insufficient temporal coverage for trend statistics
            result["$(col)_senn_slp"] = NaN
            result["$(col)_linear_slp"] = NaN
            result["$(col)_spearman_rho"] = NaN
            result["$(col)_spearman_pval"] = NaN
            result["$(col)_mk_rho"] = NaN
            result["$(col)_mk_pval"] = NaN
        end

        # Mean and median are always computed
        result["$(col)_mean"] = mean(valid_values)
        result["$(col)_median"] = median(valid_values)

        # --- Changepoint analysis (opt-in) ---
        if changepoint !== nothing
            _run_changepoint_block!(result, col, valid_years, valid_values, changepoint)
        end
    end

    return result
end


"""
Internal helper: run Pettitt test and differential metrics.
Writes 8 keys into `result` for the given column.
"""
function _run_changepoint_block!(
    result::Dict{String, Float64},
    col::String,
    valid_years::Vector{Float64},
    valid_values::Vector{Float64},
    changepoint::NamedTuple
)
    cp_start = get(changepoint, :start_year, 1980)
    cp_end = get(changepoint, :end_year, 2024)
    cp_min_obs = get(changepoint, :min_total_obs, 20)
    cp_min_seg = get(changepoint, :min_segment_obs, 10)

    # Filter to analysis window
    cp_mask = valid_years .>= cp_start .&& valid_years .<= cp_end
    cp_years = valid_years[cp_mask]
    cp_vals = valid_values[cp_mask]

    # --- Pettitt test (4 base columns) ---
    pettitt = pettitt_test(cp_years, cp_vals;
        min_total_obs=cp_min_obs, min_segment_obs=cp_min_seg)

    result["$(col)_pettitt_cp_year"] = pettitt.cp_year
    result["$(col)_pettitt_pval"] = pettitt.pval
    result["$(col)_pettitt_pre_mean"] = pettitt.pre_mean
    result["$(col)_pettitt_post_mean"] = pettitt.post_mean

    # --- Pettitt differential metrics (4 columns) ---
    pettitt_diff = segment_differential_metrics(cp_years, cp_vals, pettitt.cp_year)
    result["$(col)_pettitt_delta_mean"] = pettitt_diff.delta_mean
    result["$(col)_pettitt_pct_change"] = pettitt_diff.pct_change
    result["$(col)_pettitt_pre_mk_pval"] = pettitt_diff.pre_mk_pval
    result["$(col)_pettitt_post_mk_pval"] = pettitt_diff.post_mk_pval
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

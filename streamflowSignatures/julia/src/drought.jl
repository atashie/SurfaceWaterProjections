"""
Streamflow drought signatures.

Per-water-year drought duration and deficit against fixed percentile thresholds,
after Adelsperger et al. (in review), "A novel severity-based approach for assessing
streamflow drought characteristics and drivers". Full specification:
`docs/plans/2026-07-27-drought-signatures-plan.md`.

Method (per gage):
1. Daily Q is smoothed with a 7-day CENTERED moving average over the CONTINUOUS
   date-indexed series — not per water year, so the Oct 1 boundary introduces no
   artificial discontinuity. The window never blends across a calendar gap: it is
   applied within maximal runs of consecutive daily dates and shrinks at run edges
   (a day whose in-run window holds fewer than `smoothing_min_valid_days` values is
   NaN). Because the preprocessor hands us complete daily grids for qualifying years,
   only runs shorter than the minimum produce NaN in practice.
2. Thresholds are MAGNITUDE (non-exceedance) percentiles of the smoothed values
   pooled over the whole record — the FIXED method of the paper — computed with the
   unbiased Weibull plotting position `p_i = i/(n+1)` (Hyndman & Fan definition 6).
   NOTE this is NOT Julia's default quantile (type 7), and the two differ most in
   exactly the low tail these thresholds occupy. Levels 2/5/10/20/30 % correspond to
   U.S. Drought Monitor classes D4/D3/D2/D1/D0; the paper's convention (and this
   repo's `Q{n}` columns) means the 10 % flow is the flow exceeded 90 % of the time.
3. For each water year and level:
   - `drought_duration_fixed_p{n}` = count of days with `Q_smooth < Q_thr` (days)
   - `drought_deficit_fixed_p{n}`  = `sum(Q_thr - Q_smooth)` over those days (mm)
   Comparison is STRICT (`<`), so on a gage whose threshold is exactly 0 (heavily
   intermittent) both metrics are a valid 0 rather than counting zero-flow days.

Conventions and deliberate deviations:
- Aggregation is by WATER year (Oct 1 - Sep 30), per repo convention; the paper uses
  the climate year (Apr 1 - Mar 31), so droughts spanning Oct 1 are split across two
  annual values here. Day-level indicators — and therefore record totals — are
  identical under either grouping; the annual series and trend statistics are not.
- Thresholds come from the run's own window, matching the paper's per-analysis-period
  thresholds. Both metrics are therefore RECORD-DEPENDENT: comparable within a
  product, never across products with different windows.
- The paper's VARIABLE (day-of-year) threshold method is intentionally not
  implemented: its per-day sample is one value per year, so at 20-46 years of record
  the 2 % level falls below the smallest Weibull plotting position (see plan §2.5).
- Both metrics are dense annual series whose ZEROS are meaningful (a wet year has no
  drought days), so no snow-style record-anchored trend gate is needed.
- `drought_deficit_*` is unit-carrying (mm only when Q is area-normalized); the
  duration metrics are scale-invariant.

References: Adelsperger et al. (in review). Laaha, G. et al. (2017) — unbiased
plotting positions for low-flow frequency analysis.
"""

# Metric names are derived from the configured levels so config and columns cannot
# drift apart. Durations first, then deficits (level order as configured).
const DROUGHT_DURATION_METRICS = ["drought_duration_fixed_p$(p)" for p in CFG_DROUGHT_PERCENTILES]
const DROUGHT_DEFICIT_METRICS = ["drought_deficit_fixed_p$(p)" for p in CFG_DROUGHT_PERCENTILES]
const DROUGHT_METRICS = vcat(DROUGHT_DURATION_METRICS, DROUGHT_DEFICIT_METRICS)

# Per-gage threshold values (mm/day), emitted as scalar diagnostics — not 8-stat
# metrics — so the thresholds behind the metrics are auditable in the delivered CSV
# (and gages whose low-level threshold is exactly 0 are immediately identifiable).
const DROUGHT_THRESHOLD_SCALARS = ["drought_threshold_fixed_p$(p)" for p in CFG_DROUGHT_PERCENTILES]


"""
    weibull_quantile(x, p; below_range="na") -> Float64

Quantile at non-exceedance probability `p` using the unbiased Weibull plotting
position `p_i = i/(n+1)` (Hyndman & Fan definition 6; equivalently
`quantile(x, p; alpha=0, beta=0)`).

NaNs are dropped. `p` below the smallest plotting position `1/(n+1)` is NOT
interpolable: with `below_range="na"` the result is NaN (honest — the sample cannot
resolve that probability), with `"clamp"` it is the sample minimum (Julia's default
behavior). Above `n/(n+1)` the sample maximum is returned.
"""
function weibull_quantile(x::AbstractVector{<:Real}, p::Real; below_range::String="na")
    vals = Float64[v for v in x if !isnan(v)]
    isempty(vals) && return NaN
    sort!(vals)
    n = length(vals)

    h = (n + 1) * float(p)
    if h < 1.0
        return below_range == "clamp" ? vals[1] : NaN
    elseif h >= n
        return vals[n]
    end

    fl = floor(Int, h)
    frac = h - fl
    return vals[fl] + frac * (vals[fl + 1] - vals[fl])
end


"""
    _day_numbers(dates) -> (Vector{Int}, BitVector)

Day numbers (`Dates.value`) plus a validity mask. Unparseable/missing dates are
marked invalid and never treated as contiguous with a neighbour.
"""
function _day_numbers(dates::AbstractVector)
    n = length(dates)
    nums = zeros(Int, n)
    ok = trues(n)
    for i in 1:n
        d = dates[i]
        if ismissing(d)
            ok[i] = false
        elseif d isa Date
            # Common path: avoid relying on a Date(::Date) constructor existing
            nums[i] = Dates.value(d)
        else
            try
                nums[i] = Dates.value(Date(d))   # DateTime / ISO string / convertible
            catch
                ok[i] = false
            end
        end
    end
    return nums, ok
end


"""
    smooth_daily_flow(dates, Q; window, alignment, min_valid) -> Vector{Float64}

Moving-average smoothing of a daily series, applied within maximal runs of
consecutive calendar dates so the window never averages across a temporal gap
(rejected years, duplicated dates, or missing dates all break a run).

`alignment="center"` uses `±(window-1)÷2` days; `"trailing"` uses the `window` days
ending on the day itself. The window shrinks at run edges; a day whose in-run window
contains fewer than `min_valid` non-NaN values yields NaN.

`dates` must be sorted ascending (callers sort first).
"""
function smooth_daily_flow(
    dates::AbstractVector, Q::AbstractVector{<:Real};
    window::Int=CFG_DROUGHT_SMOOTH_WINDOW,
    alignment::String=CFG_DROUGHT_SMOOTH_ALIGNMENT,
    min_valid::Int=CFG_DROUGHT_SMOOTH_MIN_VALID
)
    n = length(Q)
    out = fill(NaN, n)
    n == 0 && return out
    length(dates) == n || throw(ArgumentError("smooth_daily_flow: dates and Q length mismatch ($(length(dates)) vs $n)"))

    nums, ok = _day_numbers(dates)
    half = (window - 1) ÷ 2

    run_start = 1
    for i in 1:n
        # A run ends at i when i is the last row, either endpoint is invalid, or the
        # next date is not exactly one day later (gap OR duplicate).
        ends_run = (i == n) || !ok[i] || !ok[i + 1] || (nums[i + 1] != nums[i] + 1)
        if ends_run
            # A row with a missing/unparseable date belongs to NO run: it must not be
            # averaged with its neighbours (it has no known position in time), so the
            # run stops before it and its own output stays NaN.
            run_end = ok[i] ? i : i - 1
            @inbounds for j in run_start:run_end
                lo, hi = if alignment == "center"
                    (max(run_start, j - half), min(run_end, j + half))
                else  # trailing
                    (max(run_start, j - window + 1), j)
                end
                s = 0.0
                c = 0
                for k in lo:hi
                    v = Q[k]
                    if !isnan(v)
                        s += v
                        c += 1
                    end
                end
                out[j] = c >= min_valid ? s / c : NaN
            end
            run_start = i + 1
        end
    end

    return out
end


"""
    calculate_drought_metrics(df::DataFrame; kwargs...) -> Dict

Calculate streamflow drought duration and deficit with trend statistics.

For each configured percentile level `p` (default 2, 5, 10, 20, 30):
- `drought_duration_fixed_p{p}` — days per water year below the fixed threshold (days)
- `drought_deficit_fixed_p{p}`  — summed departures below the threshold (mm)

Each produces 8 statistics via `generate_stats()` (+ 8 Pettitt changepoint fields when
enabled). Additionally emits per-gage scalar diagnostics
`drought_threshold_fixed_p{p}` (mm/day) — the threshold values themselves.

Both metrics are dense with meaningful zeros, so the trend-completeness gate and the
stats floor apply normally (this family is NOT exempt).

Parameters
----------
df : DataFrame
    Daily data with columns: date, Q, water_year (preprocessor output; callers pass
    data already filtered to qualifying years).

!!! note "min_years_for_threshold counts DISTINCT water_year labels"
    The threshold floor counts distinct `water_year` values present in `df`, not
    complete years — it assumes the caller has already applied year qualification
    (`preprocess_daily_data`, as the benchmark runners do). A direct caller passing
    ten one-day years would clear the floor with ten observations.

Returns
-------
Dict{String, Float64}
    2 x n_levels metrics x 8 statistics (+ changepoint fields) + n_levels scalars
"""
function calculate_drought_metrics(
    df::DataFrame;
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing,
    min_values_for_stats::Union{Nothing, Int}=nothing,
    changepoint::Union{Nothing, NamedTuple}=nothing,
    collector::Union{Nothing, AnnualCollector}=nothing
)
    result = Dict{String, Float64}()

    # Every early-return path emits the IDENTICAL key set (including changepoint
    # fields when enabled) via the same explicit-value_cols machinery — schema contract.
    function _all_nan_result()
        empty_annual = DataFrame(water_year = Int[])
        for m in DROUGHT_METRICS
            empty_annual[!, Symbol(m)] = Float64[]
        end
        merge!(result, generate_stats(empty_annual; value_cols=DROUGHT_METRICS,
                                     trend_completeness=trend_completeness,
                                     decade_completeness=decade_completeness,
                                     min_values_for_stats=min_values_for_stats,
                                     changepoint=changepoint, collector=collector))
        for s in DROUGHT_THRESHOLD_SCALARS
            result[s] = NaN
        end
        return result
    end

    valid, missing_cols = validate_columns(df, ["Q", "water_year", "date"])
    if !valid
        @warn "calculate_drought_metrics: Missing columns: $missing_cols"
        return _all_nan_result()
    end
    if nrow(df) == 0
        return _all_nan_result()
    end

    dfs = sort(df, :date)
    Q_smooth = smooth_daily_flow(dfs.date, coalesce_q(dfs.Q))

    years = sort(unique(Int[y for y in dfs.water_year if !ismissing(y)]))
    if isempty(years)
        return _all_nan_result()
    end

    # Fixed thresholds: whole-record percentiles of the smoothed series. Short records
    # get NaN thresholds (and hence NaN metrics) rather than an unstable low tail.
    n_levels = length(CFG_DROUGHT_PERCENTILES)
    thresholds = fill(NaN, n_levels)
    if length(years) >= CFG_DROUGHT_MIN_YEARS
        pool = filter(!isnan, Q_smooth)
        if !isempty(pool)
            for (k, p) in enumerate(CFG_DROUGHT_PERCENTILES)
                thresholds[k] = weibull_quantile(pool, p / 100;
                                                 below_range=CFG_DROUGHT_BELOW_RANGE_POLICY)
            end
        end
    end

    annual_data = DataFrame(water_year = years)
    for m in DROUGHT_METRICS
        annual_data[!, Symbol(m)] = fill(NaN, length(years))
    end

    for (yr_idx, yr) in enumerate(years)
        year_mask = coalesce.(dfs.water_year .== yr, false)
        any(year_mask) || continue
        qs = Q_smooth[year_mask]

        for k in 1:n_levels
            thr = thresholds[k]
            isnan(thr) && continue

            dur = 0
            deficit = 0.0
            @inbounds for v in qs
                if !isnan(v) && v < thr
                    dur += 1
                    deficit += thr - v
                end
            end

            annual_data[yr_idx, Symbol(DROUGHT_DURATION_METRICS[k])] = Float64(dur)
            annual_data[yr_idx, Symbol(DROUGHT_DEFICIT_METRICS[k])] = deficit
        end
    end

    merge!(result, generate_stats(annual_data; value_cols=DROUGHT_METRICS,
                                 trend_completeness=trend_completeness,
                                 decade_completeness=decade_completeness,
                                 min_values_for_stats=min_values_for_stats,
                                 changepoint=changepoint, collector=collector))

    for (k, s) in enumerate(DROUGHT_THRESHOLD_SCALARS)
        result[s] = thresholds[k]
    end

    return result
end

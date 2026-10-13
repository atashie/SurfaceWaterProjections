"""
Climate metrics (guidelines §3.9, synced 2026-10-12; decisions: user, 2026-10-12).

Thirty-nine per-water-year metrics from basin-mean Daymet daily precipitation
(`prcp` → PPT, mm/day) and air temperature (`tmin` → TMIN, `tmax` → TMAX, °C), in
two function families:

- `calculate_precipitation_metrics` — 20 bases on PPT-valid years (the frame the
  Q-to-PPT families receive; pure climate, so it runs regardless of `area_normalized`).
- `calculate_temperature_metrics` — 19 bases on an explicit `temp_data` frame filtered
  to temperature-valid years (`preprocess_daily_data`'s `valid_temp_years`; the snow
  family's opt-in precedent — TMIN/TMAX inside the main gage frame are never used
  implicitly).

Conventions (design record: docs/plans/2026-10-12-climate-metrics-plan.md):
- Seasons are `CFG_NA_SEASON_MONTHS` (DJF, MAM, JJA, SON) applied WITHIN the water year,
  the Qwin/Qspr/Qsum/Qfal convention: `fal` = Oct + Nov of the water-year start plus Sep
  of its end. Suffixes `win`, `spr`, `sum`, `fal`.
- dowy 1 = Oct 1. Ties take the first occurrence (maximum day, longest dry spell).
- A year is used only when its daily grid is complete (365/366 rows, dowy == 1..n, no
  NaN); anything else leaves the whole row NaN (defensive — callers pass qualified years).
- No Q-derived seasonal mask: these are climate series, not streamflow.
- A dry day has P <= `dry_day_threshold_mm`; wet-day counts are STRICTLY above each
  threshold; the multi-day maximum is a window inside the water year; growing degree
  days are Σ max(0, (TMAX + TMIN)/2 − `gdd_base_c`); frost is TMIN < `frost_temp_c`.
- Gating: standard trend completeness and stats floor via `generate_stats`; the
  record-anchored decade gate (snow precedent) on `frost_first_dowy` / `frost_last_dowy`,
  anchored to the temperature-valid years.
- Config: the `climate` section (CFG_CLIMATE_*); an absent section disables the family
  in the orchestrator (no climate columns — the additivity-proof switch).
"""

const _CLIMATE_SEASON_KEYS = ["win", "spr", "sum", "fal"]

# Column-name fragment for a wet-day threshold, identical in the three languages
# (review round 2): an integral value prints as an integer ("10"), anything else with
# C's %g ("0.5", "1e-05"), then "." -> "p" and "-" -> "m" ("0p5", "1em05").
function _wet_threshold_label(thr::Real)
    s = (isinteger(thr) && abs(thr) < 1e15) ? string(Int(thr)) : @sprintf("%g", Float64(thr))
    return replace(replace(s, "." => "p"), "-" => "m")
end

# Config-derived column names (the D{p}_day pattern): a changed window, fraction or
# threshold list changes the names, never silently the meaning behind a fixed name.
const CLIMATE_MULTIDAY_MAX_COL = "prcp_$(CFG_CLIMATE_MULTIDAY_WINDOW)d_max"
const CLIMATE_CUM_DOWY_COL = "prcp_$(round(Int, 100 * CFG_CLIMATE_CUM_FRACTION))_dowy"
const CLIMATE_WET_DAY_COLS = ["prcp_ndays_abv_$(_wet_threshold_label(t))mm" for t in CFG_CLIMATE_WET_THRESHOLDS_MM]

const PRECIP_METRICS = vcat(
    ["prcp_annual"],
    ["prcp_$(s)" for s in _CLIMATE_SEASON_KEYS],
    ["prcp_1d_max", CLIMATE_MULTIDAY_MAX_COL, "prcp_max_dowy", CLIMATE_CUM_DOWY_COL, "prcp_cm_dowy"],
    CLIMATE_WET_DAY_COLS,
    ["prcp_annual_sd", "prcp_annual_cv"],
    ["prcp_intermit_annual", "prcp_intermit_annual_dowy_start", "prcp_intermit_annual_dowy_end"],
    ["prcp_gini", "prcp_seasonality_index"],
)

if length(unique(PRECIP_METRICS)) != length(PRECIP_METRICS)
    error("climate: the configured window / thresholds produce duplicate precipitation column names: $(PRECIP_METRICS)")
end

const TEMP_METRICS = vcat(
    ["temp_max_annual", "temp_max_annual_dowy", "temp_min_annual", "temp_min_annual_dowy"],
    ["temp_max_$(s)" for s in _CLIMATE_SEASON_KEYS],
    ["temp_min_$(s)" for s in _CLIMATE_SEASON_KEYS],
    ["frost_first_dowy", "frost_last_dowy"],
    ["gdd_cum_annual"],
    ["gdd_cum_$(s)" for s in _CLIMATE_SEASON_KEYS],
)

# Timing metrics that are NaN in frost-free years: gated like the snow timing metrics
const TEMP_RECORD_GATED_METRICS = ["frost_first_dowy", "frost_last_dowy"]

# Relative tolerance of the cumulative-fraction crossing (review round 2, 2026-10-12):
# two equal rainfall blocks are mathematically tied at 50 % but their sequential
# Float64 sums differ in the last bits (0.7 vs 0.7000000000000001), which moved the
# crossing by hundreds of days. The crossing is therefore "cumulative >= fraction x total
# minus 1e-12 x total". 1e-12 of an annual total is physically nothing; the closest real
# margin in the 360,585 production water years is 1.3e-9 (Sonnet review), so no production
# value is affected. D50_day keeps the strict comparison (its own precedent).
const CLIMATE_CROSSING_REL_TOL = 1e-12


"""
    _seqsum(x) -> Float64

Plain left-to-right Float64 accumulation. Every sum, mean and cumulative sum in this
family uses it, so (a) a cumulative sum's last element is exactly the total it is
compared against (review finding 1, 2026-10-12: `cumsum` and a pairwise `sum` of the
same series round differently and can miss an exact crossing), and (b) the annual
values are bit-identical across Julia, Python and R, whose library sums differ in
accumulation order.
"""
function _seqsum(x)
    s = 0.0
    for v in x
        s += v
    end
    return s
end

function _seqcumsum(x::AbstractVector{Float64})
    out = Vector{Float64}(undef, length(x))
    s = 0.0
    for i in eachindex(x)
        s += x[i]
        out[i] = s
    end
    return out
end

"""
    _water_year_days(yr) -> Int

Number of days in water year `yr` (Oct 1 of yr-1 .. Sep 30 of yr): 366 when Feb 29 of
`yr` exists, else 365.
"""
_water_year_days(yr::Integer) = Dates.value(Date(yr, 9, 30) - Date(yr - 1, 10, 1)) + 1

"""
    _real_vector(x) -> Vector{Float64}

Observation contract of the climate kernels (review round 3): a value is an observation
when it is a `Real` (Bool included); `missing`, strings, complex numbers and anything
else are NaN. (`coalesce_q` accepts any `Number`, so a complex value with a zero
imaginary part slipped through and one with a non-zero part raised.)
"""
function _real_vector(x::AbstractVector)
    out = Vector{Float64}(undef, length(x))
    for (i, v) in enumerate(x)
        out[i] = (!ismissing(v) && v isa Real) ? Float64(v) : NaN
    end
    return out
end


"""
    _complete_year_grid(dowy, n, yr) -> Bool

True when the year's rows are EXACTLY dowy 1..n and n is the number of days of water
year `yr` (review finding 3: accepting any n in {365, 366} let a leap water year
truncated to 365 rows pass). Endpoint checks alone would accept a year with a
duplicated interior day masking a missing one.
"""
function _complete_year_grid(dowy::AbstractVector, n::Int, yr::Integer)
    n == _water_year_days(yr) || return false
    any(ismissing, dowy) && return false
    return all(i -> dowy[i] == i, 1:n)
end

"""
    _grid_months(yr, n) -> Vector{Int}

Calendar month of each day of water year `yr` on the complete grid (dowy 1..n),
derived from the date rather than read from a `month` column (review finding 4: a
missing or wrong month value silently dropped days from every season).
"""
_grid_months(yr::Integer, n::Int) = [month(Date(yr - 1, 10, 1) + Day(i - 1)) for i in 1:n]


"""
    _empty_family_stats(metrics; kwargs...) -> Dict

The full NaN key set of a family (16 fields per metric when changepoint is enabled),
emitted through generate_stats on an empty annual frame — the schema contract for
missing inputs and direct-API misuse (snow precedent).
"""
function _empty_family_stats(metrics::Vector{String};
        trend_completeness=nothing, decade_completeness=nothing,
        min_values_for_stats=nothing, changepoint=nothing, collector=nothing)
    empty_annual = DataFrame(water_year = Int[])
    for m in metrics
        empty_annual[!, Symbol(m)] = Float64[]
    end
    return generate_stats(empty_annual; value_cols=metrics,
                          trend_completeness=trend_completeness,
                          decade_completeness=decade_completeness,
                          min_values_for_stats=min_values_for_stats,
                          changepoint=changepoint, collector=collector)
end


"""
    _record_anchored_force_skip(annual_data, anchor_col, gated, decade_completeness)

Record-anchored decade gate (the snow family's rule, `snow.jl`, July 2026): the gated
metrics need >= `decade_completeness` of the anchor years (rows where `anchor_col` is
non-NaN — the family's valid, grid-complete years) in BOTH the first and the last decade
of the record to be non-NaN, else their six trend statistics are forced to NaN
(mean/median and Pettitt fields still computed). Record span < 10 years → no gate.
Returns the Set of metric names to skip, or `nothing`.
"""
function _record_anchored_force_skip(annual_data::DataFrame, anchor_col::Symbol,
                                     gated::Vector{String}, decade_completeness::Float64)
    anchor_mask = .!isnan.(annual_data[!, anchor_col])
    anchor_years = Int.(annual_data.water_year[anchor_mask])
    isempty(anchor_years) && return nothing
    rec_min, rec_max = extrema(anchor_years)
    (rec_max - rec_min + 1 >= 10) || return nothing
    windows = ((rec_min, rec_min + 9), (rec_max - 9, rec_max))
    skip = Set{String}()
    for m in gated
        vals = annual_data[!, Symbol(m)]
        for (lo, hi) in windows
            den = 0
            num = 0
            for i in 1:nrow(annual_data)
                anchor_mask[i] || continue
                yr = Int(annual_data.water_year[i])
                (lo <= yr <= hi) || continue
                den += 1
                isnan(vals[i]) || (num += 1)
            end
            if den > 0 && (num / den) < decade_completeness
                push!(skip, m)
                break
            end
        end
    end
    return isempty(skip) ? nothing : skip
end


"""
    calculate_precipitation_metrics(df::DataFrame; kwargs...) -> Dict

Twenty precipitation metrics with trend statistics from daily basin-mean precipitation
(PPT, mm/day), per water year:

- `prcp_annual`, `prcp_win/spr/sum/fal` — totals (mm) over the water year / the season's
  calendar months within it
- `prcp_1d_max` (mm), `prcp_<w>d_max` (max total of `w` consecutive days inside the water
  year, mm), `prcp_max_dowy` (first day of the 1-day maximum)
- `prcp_<f>_dowy` — first day whose cumulative precipitation reaches fraction `f` of the
  water-year total; `prcp_cm_dowy` — centre of mass Σ(dowy·P)/ΣP (both NaN when ΣP <= 0)
- `prcp_ndays_abv_<t>mm` — days with P > t, for each configured threshold
- `prcp_annual_sd`, `prcp_annual_cv` — sample sd (n − 1) of the daily values and sd/mean
  (cv NaN when the mean is not > 0)
- `prcp_intermit_annual` — longest run of dry days (P <= dry threshold; 0 when every day
  is wet); `prcp_intermit_annual_dowy_start` / `_dowy_end` — first and last day of the
  FIRST longest run (NaN when there is no dry day)
- `prcp_gini` — Gini index of the daily values over all n days, zeros included:
  G = Σᵢ (2i − n − 1)·y₍ᵢ₎ / (n·Σy) with y sorted ascending (0 uniform … (n−1)/n one wet
  day; NaN when ΣP <= 0)
- `prcp_seasonality_index` — Walsh & Lawler (1981): SI = (1/R)·Σₘ₌₁¹² |Xₘ − R/12| over
  the water year's 12 calendar-month totals Xₘ, R = ΣXₘ (0 … 1.83; NaN when R <= 0)

Callers pass the PPT-valid frame (`climate_data`). A year whose grid is incomplete or
holds a NaN is left NaN in every metric. Each metric produces the 8 statistics (+ 8
Pettitt fields when enabled) via `generate_stats()`.
"""
function calculate_precipitation_metrics(
    df::DataFrame;
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing,
    min_values_for_stats::Union{Nothing, Int}=nothing,
    changepoint::Union{Nothing, NamedTuple}=nothing,
    collector::Union{Nothing, AnnualCollector}=nothing
)
    result = Dict{String, Float64}()

    valid, missing_cols = validate_columns(df, ["PPT", "water_year", "dowy"])
    if !valid
        @warn "calculate_precipitation_metrics: Missing columns: $missing_cols"
        merge!(result, _empty_family_stats(PRECIP_METRICS;
            trend_completeness=trend_completeness, decade_completeness=decade_completeness,
            min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector))
        return result
    end

    dry_thr = Float64(CFG_CLIMATE_DRY_DAY_MM)
    wet_thrs = CFG_CLIMATE_WET_THRESHOLDS_MM
    w = Int(CFG_CLIMATE_MULTIDAY_WINDOW)
    frac = Float64(CFG_CLIMATE_CUM_FRACTION)
    season_months = CFG_NA_SEASON_MONTHS
    col_multi = Symbol(CLIMATE_MULTIDAY_MAX_COL)
    col_cum = Symbol(CLIMATE_CUM_DOWY_COL)
    wet_cols = Symbol.(CLIMATE_WET_DAY_COLS)

    years = unique(df.water_year)
    n_years = length(years)
    annual_data = DataFrame(water_year = years)
    for m in PRECIP_METRICS
        annual_data[!, Symbol(m)] = fill(NaN, n_years)
    end

    for (yr_idx, yr) in enumerate(years)
        year_mask = coalesce.(df.water_year .== yr, false)
        year_df = df[year_mask, :]
        nrow(year_df) == 0 && continue
        sort!(year_df, :dowy)

        # `.+ 0.0` turns a signed -0.0 into +0.0 so extrema and ties do not depend on the
        # sign of zero (review round 2); the guard requires FINITE values (±Inf is not data)
        P = _real_vector(year_df.PPT) .+ 0.0
        dowy = year_df.dowy
        n = length(P)

        # Complete-grid + finite guard (snow precedent, year-aware): anything else -> NaN row
        (_complete_year_grid(dowy, n, Int(yr)) && all(isfinite, P)) || continue
        months = _grid_months(Int(yr), n)

        total = _seqsum(P)
        annual_data[yr_idx, :prcp_annual] = total

        # Seasonal totals over the season's calendar months within the water year
        for s in _CLIMATE_SEASON_KEYS
            smask = in.(months, Ref(Set(season_months[s])))
            annual_data[yr_idx, Symbol("prcp_$(s)")] = _seqsum(P[smask])
        end

        # 1-day maximum and its first day
        imax = argmax(P)
        annual_data[yr_idx, :prcp_1d_max] = P[imax]
        annual_data[yr_idx, :prcp_max_dowy] = Float64(dowy[imax])

        # Multi-day maximum: every window of w consecutive days inside the year
        if n >= w
            best = -Inf
            for i in 1:(n - w + 1)
                s = _seqsum(view(P, i:(i + w - 1)))
                if s > best
                    best = s
                end
            end
            annual_data[yr_idx, col_multi] = best
        end

        # Timing of the cumulative total, centre of mass, Gini, seasonality index.
        # The cumulative sum and `total` share one accumulation order (cum[end] == total);
        # the crossing test carries CLIMATE_CROSSING_REL_TOL so a mathematical tie is not
        # decided by the last bit of the sums.
        if total > 0.0
            cum = _seqcumsum(P)
            thr = frac * total - CLIMATE_CROSSING_REL_TOL * total
            idx = findfirst(c -> c >= thr, cum)
            if idx !== nothing
                annual_data[yr_idx, col_cum] = Float64(dowy[idx])
            end
            annual_data[yr_idx, :prcp_cm_dowy] = _seqsum(Float64.(dowy) .* P) / total

            y = sort(P)
            annual_data[yr_idx, :prcp_gini] = _seqsum((2 * i - n - 1) * y[i] for i in 1:n) / (n * total)

            X = [_seqsum(P[months .== m]) for m in 1:12]
            R = _seqsum(X)
            if R > 0.0
                annual_data[yr_idx, :prcp_seasonality_index] = _seqsum(abs.(X .- R / 12)) / R
            end
        end

        # Wet-day counts, strictly above each threshold
        for (k, t) in enumerate(wet_thrs)
            annual_data[yr_idx, wet_cols[k]] = Float64(count(p -> p > t, P))
        end

        # Within-year spread of the daily values (sample sd, n - 1; sequential sums)
        mu = total / n
        sd = sqrt(_seqsum((p - mu)^2 for p in P) / (n - 1))
        annual_data[yr_idx, :prcp_annual_sd] = sd
        annual_data[yr_idx, :prcp_annual_cv] = mu > 0.0 ? sd / mu : NaN

        # Longest dry spell (first on ties) and its first / last day
        dry = P .<= dry_thr
        best_len = 0
        best_start = 0
        i = 1
        while i <= n
            if dry[i]
                j = i
                while j < n && dry[j + 1]
                    j += 1
                end
                len = j - i + 1
                if len > best_len
                    best_len = len
                    best_start = i
                end
                i = j + 1
            else
                i += 1
            end
        end
        annual_data[yr_idx, :prcp_intermit_annual] = Float64(best_len)
        if best_len > 0
            annual_data[yr_idx, :prcp_intermit_annual_dowy_start] = Float64(dowy[best_start])
            annual_data[yr_idx, :prcp_intermit_annual_dowy_end] = Float64(dowy[best_start + best_len - 1])
        end
    end

    merge!(result, generate_stats(annual_data; value_cols=PRECIP_METRICS,
                                  trend_completeness=trend_completeness,
                                  decade_completeness=decade_completeness,
                                  min_values_for_stats=min_values_for_stats,
                                  changepoint=changepoint, collector=collector))
    return result
end


"""
    calculate_temperature_metrics(df::DataFrame; kwargs...) -> Dict

Nineteen temperature metrics with trend statistics from daily basin-mean minimum and
maximum air temperature (TMIN, TMAX, °C), per water year:

- `temp_max_annual`, `temp_max_annual_dowy` — maximum of the daily TMAX and the first
  day it occurs; `temp_min_annual`, `temp_min_annual_dowy` — minimum of the daily TMIN
  and its first day
- `temp_max_win/spr/sum/fal`, `temp_min_win/spr/sum/fal` — the same extremes over the
  season's calendar months within the water year
- `frost_first_dowy` — first day with TMIN < `frost_temp_c` (NaN when none);
  `frost_last_dowy` — last such day among the days from Oct 1 through the end of
  `frost_last_end_month` (NaN when none)
- `gdd_cum_annual`, `gdd_cum_win/spr/sum/fal` — growing degree days
  Σ max(0, (TMAX + TMIN)/2 − `gdd_base_c`) over the water year / the season (°C·day)

Callers pass the temperature-valid frame (`temp_data`; in `calculate_all_signatures`
this function runs ONLY when such a frame is provided). A year whose grid is incomplete
or holds a NaN in either series is left NaN in every metric. The record-anchored decade
gate applies to the two frost metrics when `climate.record_anchored_decade_gate` is set.
"""
function calculate_temperature_metrics(
    df::DataFrame;
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing,
    min_values_for_stats::Union{Nothing, Int}=nothing,
    changepoint::Union{Nothing, NamedTuple}=nothing,
    collector::Union{Nothing, AnnualCollector}=nothing
)
    result = Dict{String, Float64}()

    valid, missing_cols = validate_columns(df, ["TMIN", "TMAX", "water_year", "dowy"])
    if !valid
        @warn "calculate_temperature_metrics: Missing columns: $missing_cols"
        merge!(result, _empty_family_stats(TEMP_METRICS;
            trend_completeness=trend_completeness, decade_completeness=decade_completeness,
            min_values_for_stats=min_values_for_stats, changepoint=changepoint, collector=collector))
        return result
    end

    frost_t = Float64(CFG_CLIMATE_FROST_TEMP_C)
    end_month = Int(CFG_CLIMATE_FROST_LAST_END_MONTH)
    base = Float64(CFG_CLIMATE_GDD_BASE_C)
    season_months = CFG_NA_SEASON_MONTHS

    years = unique(df.water_year)
    n_years = length(years)
    annual_data = DataFrame(water_year = years)
    for m in TEMP_METRICS
        annual_data[!, Symbol(m)] = fill(NaN, n_years)
    end

    for (yr_idx, yr) in enumerate(years)
        year_mask = coalesce.(df.water_year .== yr, false)
        year_df = df[year_mask, :]
        nrow(year_df) == 0 && continue
        sort!(year_df, :dowy)

        TN = _real_vector(year_df.TMIN) .+ 0.0   # signed zero normalised, as for P
        TX = _real_vector(year_df.TMAX) .+ 0.0
        dowy = year_df.dowy
        n = length(TN)

        (_complete_year_grid(dowy, n, Int(yr)) && all(isfinite, TN) && all(isfinite, TX)) || continue
        months = _grid_months(Int(yr), n)

        # Annual extremes and their first day
        imax = argmax(TX)
        imin = argmin(TN)
        annual_data[yr_idx, :temp_max_annual] = TX[imax]
        annual_data[yr_idx, :temp_max_annual_dowy] = Float64(dowy[imax])
        annual_data[yr_idx, :temp_min_annual] = TN[imin]
        annual_data[yr_idx, :temp_min_annual_dowy] = Float64(dowy[imin])

        # Growing degree days, clamped at the base temperature
        gdd = max.(0.0, (TX .+ TN) ./ 2 .- base)
        annual_data[yr_idx, :gdd_cum_annual] = _seqsum(gdd)

        # Seasonal extremes and GDD over the season's months within the water year
        for s in _CLIMATE_SEASON_KEYS
            smask = in.(months, Ref(Set(season_months[s])))
            if any(smask)
                annual_data[yr_idx, Symbol("temp_max_$(s)")] = maximum(TX[smask])
                annual_data[yr_idx, Symbol("temp_min_$(s)")] = minimum(TN[smask])
                annual_data[yr_idx, Symbol("gdd_cum_$(s)")] = _seqsum(gdd[smask])
            end
        end

        # Frost: first day below the threshold; last such day up to the end of
        # frost_last_end_month (a month >= 10 belongs to the water year's first
        # calendar year). dowy == 1..n after the grid guard, so the cut-off is an index.
        frost = TN .< frost_t
        f = findfirst(frost)
        if f !== nothing
            annual_data[yr_idx, :frost_first_dowy] = Float64(dowy[f])
        end
        cal_year = end_month >= 10 ? yr - 1 : yr
        cutoff_dowy = Dates.value(Date(cal_year, end_month, Dates.daysinmonth(cal_year, end_month)) -
                                  Date(yr - 1, 10, 1)) + 1
        last_idx = findlast(i -> frost[i] && dowy[i] <= cutoff_dowy, 1:n)
        if last_idx !== nothing
            annual_data[yr_idx, :frost_last_dowy] = Float64(dowy[last_idx])
        end
    end

    force_skip = nothing
    if CFG_CLIMATE_RECORD_DECADE_GATE && decade_completeness !== nothing
        force_skip = _record_anchored_force_skip(annual_data, :temp_max_annual,
                                                 TEMP_RECORD_GATED_METRICS, decade_completeness)
    end

    merge!(result, generate_stats(annual_data; value_cols=TEMP_METRICS,
                                  trend_completeness=trend_completeness,
                                  decade_completeness=decade_completeness,
                                  min_values_for_stats=min_values_for_stats,
                                  changepoint=changepoint, collector=collector,
                                  force_skip_trends=force_skip))
    return result
end

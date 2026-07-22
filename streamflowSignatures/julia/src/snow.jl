"""
Snow metrics signatures.

Fourteen per-water-year metrics computed from daily snow water equivalent (SWE,
mm; Daymet `swe`). ALL metrics operate on a thresholded series
`SWE* = (SWE >= CFG_SNOW_SWE_THRESHOLD_MM) ? SWE : 0.0` — days below the threshold
are treated as snow-free for durations AND magnitudes (domain decision 2026-07-14).

Definitions (full specification: docs/plans/snow_signatures_plan.md):
- Snow-day mask: SWE* > 0; spells are maximal runs of consecutive snow days.
- Anchor spell: the spell containing the annual maximum of SWE* (first-max tie
  rule). snow_on_dowy / snow_off_dowy are censored (NaN) when the anchor spell
  touches the first / last day of the water year.
- Melt increments m_t = max(0, SWE*_{t-1} - SWE*_t), attributed to day t;
  within-year only (no attribution across the Sep 30 -> Oct 1 boundary).
- SSM (snow seasonality metric; Hatchett 2021, Hydrology 8(1):32): spells of
  >= CFG_SNOW_SEASONAL_MIN_DAYS days are seasonal, shorter spells ephemeral;
  SSM = (seasonal_days - ephemeral_days) / total snow days, in [-1, +1].

Year qualification: callers pass data already filtered to SWE-valid years
(`preprocess_daily_data`'s `valid_swe_years`). Defensively, any year whose SWE
series has NaNs or is not a complete Oct 1 - Sep 30 daily grid is left entirely
NaN (all 14 metrics).

References: Hatchett, B.J. (2021). Seasonal and Ephemeral Snowpacks of the
Conterminous United States. Hydrology 8(1), 32. Petersky & Harpold (2018), HESS.
"""

# Snow parameters loaded from config.jl (CFG_SNOW_*)

const SNOW_METRICS = [
    "swe_max",
    "swe_max_dowy",
    "snow_cover_days",
    "snow_on_dowy",
    "snow_off_dowy",
    "melt_season_days",
    "melt_rate",
    "ssm",
    "swe_apr1",
    "melt_before_peak",
    "melt_before_peak_pct",
    "melt_before_peak_to_max_swe",
    "melt_com_dowy",
    "swe_max_to_ppt"
]

# Threshold-dependent metrics (NaN for operationally snow-free years) gated by the
# record-anchored decade gate (July 2026): a trend requires computable values in
# >= decade_completeness of the SWE-valid years in BOTH the first and last decade
# of the gage's SWE record. Magnitude metrics (valid zeros) are NOT gated — their
# dense series carry the snow-decline signal legitimately.
const SNOW_RECORD_GATED_METRICS = [
    "swe_max_dowy",
    "snow_on_dowy",
    "snow_off_dowy",
    "melt_season_days",
    "melt_rate",
    "ssm",
    "melt_before_peak",
    "melt_before_peak_pct",
    "melt_before_peak_to_max_swe",
    "melt_com_dowy"
]


"""
    snow_spells(mask::AbstractVector{Bool}) -> Vector{UnitRange{Int}}

Maximal runs of consecutive `true` values, as index ranges.
"""
function snow_spells(mask::AbstractVector{Bool})
    spells = UnitRange{Int}[]
    i = 1
    n = length(mask)
    while i <= n
        if mask[i]
            j = i
            while j < n && mask[j + 1]
                j += 1
            end
            push!(spells, i:j)
            i = j + 1
        else
            i += 1
        end
    end
    return spells
end


"""
    calculate_snow_metrics(df::DataFrame; kwargs...) -> Dict

Calculate 14 snow metrics with trend statistics from daily SWE.

Metrics (per water year; see module docstring for conventions):
- swe_max (mm), swe_max_dowy, snow_cover_days, swe_apr1 (mm)
- snow_on_dowy, snow_off_dowy, melt_season_days, melt_rate (mm/day)
- ssm, melt_before_peak (mm), melt_before_peak_pct (%),
  melt_before_peak_to_max_swe, melt_com_dowy
- swe_max_to_ppt (requires PPT; year must be climate-qualified)

Magnitude metrics (swe_max, snow_cover_days, swe_apr1, swe_max_to_ppt) emit
valid ZEROS for operationally snow-free years; timing/melt/regime metrics emit
NaN. Each metric produces 8 statistics via generate_stats() (+ changepoint
fields when enabled).

Callers must pass data filtered to SWE-valid years. In `calculate_all_signatures`
this function only runs when an explicit `snow_data` frame is provided — an SWE
column inside the main gage frame is never used implicitly.

Parameters
----------
df : DataFrame
    Daily data with columns: water_year, dowy, SWE (+ optional PPT)
valid_climate_years : Vector{Int} or nothing
    Years qualified for PPT use (metric swe_max_to_ppt). When `nothing`
    (legacy path), the annual-PPT convention of the runoff ratios applies
    directly: sum of the year's non-NaN PPT days, subject to the minimum floor.

Returns
-------
Dict{String, Float64}
    14 metrics x 8 statistics (+ 8 changepoint fields each when enabled)
"""
function calculate_snow_metrics(
    df::DataFrame;
    valid_climate_years::Union{Nothing, Vector{Int}}=nothing,
    trend_completeness::Union{Nothing, Float64}=nothing,
    decade_completeness::Union{Nothing, Float64}=nothing, min_values_for_stats::Union{Nothing, Int}=nothing,
    changepoint::Union{Nothing, NamedTuple}=nothing,
    collector::Union{Nothing, AnnualCollector}=nothing
)
    result = Dict{String, Float64}()

    # Defensive: callers gate on SWE presence; direct-API misuse still emits the
    # IDENTICAL key set (incl. changepoint fields when enabled) via the same
    # explicit-value_cols machinery as the 0-row case — schema contract
    if !("SWE" in names(df))
        empty_annual = DataFrame(water_year = Int[])
        for m in SNOW_METRICS
            empty_annual[!, Symbol(m)] = Float64[]
        end
        merge!(result, generate_stats(empty_annual; value_cols=SNOW_METRICS,
                                      trend_completeness=trend_completeness,
                                      decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats,
                                      changepoint=changepoint, collector=collector))
        return result
    end

    thr = Float64(CFG_SNOW_SWE_THRESHOLD_MM)
    seasonal_min = Int(CFG_SNOW_SEASONAL_MIN_DAYS)
    com_frac = Float64(CFG_SNOW_MELT_COM_FRACTION)
    min_ppt = Float64(CFG_SNOW_MIN_ANNUAL_PPT_MM)

    has_ppt = "PPT" in names(df)
    climate_set = valid_climate_years === nothing ? nothing : Set(valid_climate_years)

    years = unique(df.water_year)
    n_years = length(years)
    annual_data = DataFrame(water_year = years)
    for m in SNOW_METRICS
        annual_data[!, Symbol(m)] = fill(NaN, n_years)
    end

    for (yr_idx, yr) in enumerate(years)
        year_mask = coalesce.(df.water_year .== yr, false)
        year_df = df[year_mask, :]
        if nrow(year_df) == 0
            continue
        end
        sort!(year_df, :dowy)

        swe_raw = coalesce_q(year_df.SWE)
        dowy = year_df.dowy
        n = length(swe_raw)

        # Complete-grid + no-NaN guard: metrics are only computed on full water
        # years whose days are EXACTLY the sequence 1..n (n in {365, 366}) with no
        # missing SWE. Endpoint checks alone would accept a year with a duplicated
        # interior day masking a missing one. Callers normally guarantee this via
        # valid_swe_years; anything else -> NaN row.
        if !(n == 365 || n == 366) || any(ismissing, dowy) ||
           !all(i -> dowy[i] == i, 1:n) || any(isnan, swe_raw)
            continue
        end

        # Thresholded series: sub-threshold days are snow-free everywhere
        swe_star = [s >= thr ? s : 0.0 for s in swe_raw]
        mask = swe_star .> 0.0
        n_snow = count(mask)

        # --- Magnitude metrics (valid zeros for operationally snow-free years) ---
        swe_max_val = n_snow > 0 ? maximum(swe_star) : 0.0
        annual_data[yr_idx, :swe_max] = swe_max_val
        annual_data[yr_idx, :snow_cover_days] = Float64(n_snow)

        # April 1 SWE: calendar April 1 of water year yr (leap-safe via Dates
        # arithmetic; dowy == 1..n after the grid guard, so this indexes directly)
        apr1_dowy = Dates.value(Date(yr, 4, 1) - Date(yr - 1, 10, 1)) + 1
        annual_data[yr_idx, :swe_apr1] = swe_star[apr1_dowy]

        # --- Timing / melt / regime metrics (need a qualifying snowpack) ---
        if n_snow > 0
            peak_idx = argmax(swe_star)  # first occurrence on ties
            annual_data[yr_idx, :swe_max_dowy] = Float64(dowy[peak_idx])

            spells = snow_spells(mask)

            # SSM pools ALL spells in the year
            seasonal_days = 0
            ephemeral_days = 0
            for sp in spells
                if length(sp) >= seasonal_min
                    seasonal_days += length(sp)
                else
                    ephemeral_days += length(sp)
                end
            end
            annual_data[yr_idx, :ssm] = (seasonal_days - ephemeral_days) / Float64(n_snow)

            # Anchor spell = the spell containing the peak
            anchor = spells[findfirst(sp -> peak_idx in sp, spells)]

            # Censoring: an anchor spell touching a water-year boundary -> NaN
            snow_on = first(anchor) == 1 ? NaN : Float64(dowy[first(anchor)])
            snow_off = last(anchor) == n ? NaN : Float64(dowy[last(anchor)] + 1)
            annual_data[yr_idx, :snow_on_dowy] = snow_on
            annual_data[yr_idx, :snow_off_dowy] = snow_off

            if !isnan(snow_off)
                melt_days = snow_off - Float64(dowy[peak_idx])
                annual_data[yr_idx, :melt_season_days] = melt_days
                annual_data[yr_idx, :melt_rate] = swe_max_val / melt_days
            end

            # Melt increments, attributed to the day the drop lands (day t).
            # With the first-max tie rule, m at the peak day itself is 0.
            m = zeros(Float64, n)
            for t in 2:n
                d = swe_star[t - 1] - swe_star[t]
                if d > 0.0
                    m[t] = d
                end
            end
            total_melt = sum(m)
            melt_before = sum(view(m, 1:peak_idx))
            annual_data[yr_idx, :melt_before_peak] = melt_before
            annual_data[yr_idx, :melt_before_peak_to_max_swe] = melt_before / swe_max_val
            if total_melt > 0.0
                annual_data[yr_idx, :melt_before_peak_pct] = 100.0 * melt_before / total_melt
                cum = 0.0
                for t in 1:n
                    cum += m[t]
                    if cum >= com_frac * total_melt
                        annual_data[yr_idx, :melt_com_dowy] = Float64(dowy[t])
                        break
                    end
                end
            end
        end

        # --- swe_max_to_ppt: needs a PPT-qualified year ---
        if has_ppt && (climate_set === nothing || yr in climate_set)
            # Same annual-PPT convention as the runoff ratios: sum the year's
            # valid (non-NaN) PPT days, require the total above the minimum floor
            P = coalesce_q(year_df.PPT)
            p_valid = .!isnan.(P)
            total_P = any(p_valid) ? sum(P[p_valid]) : 0.0
            if total_P > min_ppt
                annual_data[yr_idx, :swe_max_to_ppt] = swe_max_val / total_P
            end
        end
    end

    # Record-anchored decade gate (July 2026, user decision 2026-07-22): anchored to
    # the SWE-valid, grid-complete years (swe_max non-NaN — dense incl. zeros), NOT
    # the metric's own span, so clustered snowy years cannot slip through. The
    # threshold is the SAME decade_completeness the streamflow gate uses (linked
    # knob). Denominators count anchor years in each window — a year with no SWE
    # data must not count against snow presence. Record span < 10 -> gate skipped.
    force_skip = nothing
    if CFG_SNOW_RECORD_DECADE_GATE && decade_completeness !== nothing
        anchor_mask = .!isnan.(annual_data.swe_max)
        anchor_years = Int.(annual_data.water_year[anchor_mask])
        if !isempty(anchor_years)
            rec_min, rec_max = extrema(anchor_years)
            if rec_max - rec_min + 1 >= 10
                windows = ((rec_min, rec_min + 9), (rec_max - 9, rec_max))
                skip = Set{String}()
                for m in SNOW_RECORD_GATED_METRICS
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
                if !isempty(skip)
                    force_skip = skip
                end
            end
        end
    end

    stats = generate_stats(annual_data; value_cols=SNOW_METRICS,
                           trend_completeness=trend_completeness,
                           decade_completeness=decade_completeness, min_values_for_stats=min_values_for_stats,
                           changepoint=changepoint, collector=collector,
                           force_skip_trends=force_skip)
    merge!(result, stats)
    return result
end

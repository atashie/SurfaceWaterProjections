"""
I/O functions for streamflow signature analysis.

Functions for reading/writing data and data validation.
"""

using Arrow
using CSV
using DataFrames
using Dates
using Parquet2

# Required columns for streamflow data
const REQUIRED_COLS = ["gage_id", "date", "Q"]
const TEMPORAL_COLS = ["water_year", "month", "dowy"]


"""
    read_parquet(path::String; normalize_columns::Bool=true) -> DataFrame

Read a parquet file into a DataFrame.

If `normalize_columns` is true (default), auto-renames common column variants
to standard names: Date->date, site_id->gage_id, prcp->PPT, swe->SWE,
tmin->TMIN, tmax->TMAX (the climate metrics, 2026-10-12).
"""
function read_parquet(path::String; normalize_columns::Bool=true)
    if !isfile(path)
        error("File not found: $path")
    end
    df = DataFrame(Parquet2.Dataset(path))

    if normalize_columns
        if "Date" in names(df) && !("date" in names(df))
            rename!(df, :Date => :date)
        end
        if "site_id" in names(df) && !("gage_id" in names(df))
            rename!(df, :site_id => :gage_id)
        end
        if "prcp" in names(df) && !("PPT" in names(df))
            rename!(df, :prcp => :PPT)
        end
        if "swe" in names(df) && !("SWE" in names(df))
            rename!(df, :swe => :SWE)
        end
        if "tmin" in names(df) && !("TMIN" in names(df))
            rename!(df, :tmin => :TMIN)
        end
        if "tmax" in names(df) && !("TMAX" in names(df))
            rename!(df, :tmax => :TMAX)
        end
    end

    return df
end


"""
    write_signatures(df::DataFrame, path::String; format=:csv)

Write signature results to file.

Supports CSV and Parquet formats.
"""
function write_signatures(df::DataFrame, path::String; format::Symbol=:csv)
    if format == :csv
        CSV.write(path, df)
    elseif format == :parquet
        Arrow.write(path, df)
    else
        error("Unsupported format: $format. Use :csv or :parquet")
    end
end


"""
    validate_schema(df::DataFrame; require_climate=false) -> (valid, missing_cols)

Validate that DataFrame has required columns.

Parameters
----------
df : DataFrame
    Input data to validate
require_climate : Bool
    Whether to require climate columns (PPT)

Returns
-------
Tuple{Bool, Vector{String}}
    (is_valid, list of missing columns)
"""
function validate_schema(df::DataFrame; require_climate::Bool=false)
    required = copy(REQUIRED_COLS)
    if require_climate
        push!(required, "PPT")
    end

    col_names = String.(names(df))
    missing_cols = [col for col in required if !(col in col_names)]

    return (isempty(missing_cols), missing_cols)
end


"""
    add_water_year_columns(df::DataFrame; date_col="date") -> DataFrame

Add water year temporal columns to DataFrame.

Adds columns:
- water_year: Water year (Oct 1 - Sep 30)
- month: Calendar month (1-12)
- dowy: Day of water year (1-366)

Parameters
----------
df : DataFrame
    Input data with date column
date_col : String
    Name of date column

Returns
-------
DataFrame
    Copy of input with added temporal columns
"""
function add_water_year_columns(df::DataFrame; date_col::String="date")
    result = copy(df)

    # Auto-detect date column if specified name not found
    if !(date_col in names(result))
        variants = filter(c -> lowercase(c) == lowercase(date_col), names(result))
        if length(variants) == 1
            date_col = variants[1]
        else
            error("Date column '$date_col' not found. Available columns: $(names(result))")
        end
    end

    dates = result[!, date_col]

    # Convert to Date if needed
    if !(eltype(dates) <: Date)
        dates = Date.(dates)
    end

    # Water year: Oct-Dec is next year's water year
    water_years = [month(d) >= 10 ? year(d) + 1 : year(d) for d in dates]

    # Month
    months = month.(dates)

    # Day of water year
    # Oct 1 = day 1, Sep 30 = day 365/366
    function calc_dowy(d::Date)
        wy = month(d) >= 10 ? year(d) + 1 : year(d)
        wy_start = Date(wy - 1, 10, 1)
        return Dates.value(d - wy_start) + 1
    end

    dowy = [calc_dowy(d) for d in dates]

    result[!, :water_year] = water_years
    result[!, :month] = months
    result[!, :dowy] = dowy

    return result
end


"""
    _interpolate_internal_gaps!(x::Vector{Float64}, max_gap::Int)

Linearly fill every INTERNAL run of NaNs of length <= max_gap whose two neighbours are
finite (the rule the PPT and SWE blocks apply inline below); leading/trailing runs and
longer runs are left as NaN. Used by the temperature block (2026-10-12).
"""
function _interpolate_internal_gaps!(x::Vector{Float64}, max_gap::Int)
    i = 1
    n = length(x)
    while i <= n
        if isnan(x[i])
            gap_start = i
            gap_end = i
            while gap_end < n && isnan(x[gap_end + 1])
                gap_end += 1
            end
            gap_len = gap_end - gap_start + 1
            if gap_start > 1 && gap_end < n && gap_len <= max_gap
                left_val = x[gap_start - 1]
                right_val = x[gap_end + 1]
                if isfinite(left_val) && isfinite(right_val)   # finite neighbours only
                    span = (gap_end + 1) - (gap_start - 1)
                    for k in gap_start:gap_end
                        x[k] = left_val + (right_val - left_val) * (k - (gap_start - 1)) / span
                    end
                end
            end
            i = gap_end + 1
        else
            i += 1
        end
    end
    return x
end


"""
    preprocess_daily_data(gage_flow::DataFrame; kwargs...) -> NamedTuple

Preprocess daily streamflow data with NA handling, gap interpolation,
and per-year quality diagnostics.

Accepts a DataFrame with columns: date, Q, water_year, month, dowy (and optional
PPT, SWE and/or the TMIN + TMAX pair).
Returns a NamedTuple with fields:
- data: DataFrame with cleaned/interpolated daily data
- valid_years: Vector{Int} of years passing all quality checks
- valid_climate_years: Vector{Int} of years also passing PPT checks
- valid_swe_years: Vector{Int} of years also passing SWE checks
- valid_temp_years: Vector{Int} of years also passing the TMIN and TMAX checks
  (2026-10-12; same thresholds, no negative-value check; evaluated, like PPT and
  SWE, only for years that passed the streamflow checks)
- rejected_years: DataFrame with columns (water_year, reason)
- seasonal_flags: DataFrame with per-year seasonal completeness flags
- diagnostics: DataFrame with per-year raw diagnostics
"""
function preprocess_daily_data(
    gage_flow::DataFrame;
    max_raw_na::Int = CFG_NA_MAX_RAW_NA_PER_YEAR,
    max_gap_days::Int = CFG_NA_MAX_GAP_DAYS,
    reject_negative::Bool = CFG_NA_REJECT_NEGATIVE_FLOW,
    reject_residual_na::Bool = CFG_NA_REJECT_RESIDUAL_NA,
    seasonal_min_fraction::Float64 = Float64(CFG_NA_SEASONAL_MIN_FRACTION),
    max_raw_na_ppt::Int = CFG_NA_MAX_RAW_NA_PPT,
    max_gap_ppt::Int = CFG_NA_MAX_GAP_PPT,
    reject_negative_ppt::Bool = CFG_NA_REJECT_NEGATIVE_PPT,
    max_raw_na_swe::Int = CFG_NA_MAX_RAW_NA_SWE,
    max_gap_swe::Int = CFG_NA_MAX_GAP_SWE,
    reject_negative_swe::Bool = CFG_NA_REJECT_NEGATIVE_SWE,
    max_raw_na_temp::Int = CFG_NA_MAX_RAW_NA_TEMP,
    max_gap_temp::Int = CFG_NA_MAX_GAP_TEMP,
    constant_sd_enabled::Bool = CFG_NA_CONSTANT_SD_ENABLED,
    constant_sd_min_days::Int = CFG_NA_CONSTANT_SD_MIN_DAYS,
    constant_sd_max_unique::Int = CFG_NA_CONSTANT_SD_MAX_UNIQUE
)
    # Policy values must be whole non-negative counts also when passed directly
    (max_raw_na >= 0 && max_gap_days >= 0 && max_raw_na_ppt >= 0 && max_gap_ppt >= 0 &&
     max_raw_na_swe >= 0 && max_gap_swe >= 0 && max_raw_na_temp >= 0 && max_gap_temp >= 0) ||
        throw(ArgumentError("preprocess_daily_data: missing-day limits and interpolation gaps must be >= 0"))

    has_ppt = "PPT" in names(gage_flow)
    has_swe = "SWE" in names(gage_flow)
    has_temp = "TMIN" in names(gage_flow) && "TMAX" in names(gage_flow)
    years = sort(unique(gage_flow.water_year))

    valid_years = Int[]
    valid_climate_years = Int[]
    valid_swe_years = Int[]
    valid_temp_years = Int[]
    rejected_rows = NamedTuple{(:water_year, :reason), Tuple{Int, String}}[]
    seasonal_rows = NamedTuple{(:water_year, :win_complete, :spr_complete, :sum_complete, :fal_complete), Tuple{Int, Bool, Bool, Bool, Bool}}[]
    diag_rows = NamedTuple{(:water_year, :na_count, :max_na_run, :interpolated_count, :negative_flag, :constant_sd_flag, :na_cause_ice, :na_cause_other), Tuple{Int, Int, Int, Int, Bool, Bool, Int, Int}}[]

    # Season month mapping (from config)
    season_months = CFG_NA_SEASON_MONTHS

    # Collect all cleaned year DataFrames
    cleaned_dfs = DataFrame[]

    for yr in years
        yr_mask = coalesce.(gage_flow.water_year .== yr, false)
        yr_df = copy(gage_flow[yr_mask, :])

        if nrow(yr_df) == 0
            push!(rejected_rows, (water_year=yr, reason="no_data"))
            continue
        end

        # --- (a) Normalize to complete daily grid (Oct 1 - Sep 30) ---
        wy_start = Date(yr - 1, 10, 1)
        wy_end = Date(yr, 9, 30)
        full_dates = collect(wy_start:Day(1):wy_end)
        expected_days = length(full_dates)

        # Build a lookup from date -> row index in yr_df
        grid = DataFrame(date = full_dates)
        grid[!, :water_year] .= yr
        grid[!, :month] = month.(full_dates)
        grid[!, :dowy] = [Dates.value(d - wy_start) + 1 for d in full_dates]

        # Left-join to get Q (and PPT) aligned to full grid
        if eltype(yr_df.date) != Date
            yr_df.date = Date.(yr_df.date)
        end
        joined = leftjoin(grid, select(yr_df, :date, :Q, (has_ppt ? [:PPT] : Symbol[])..., (has_swe ? [:SWE] : Symbol[])..., (has_temp ? [:TMIN, :TMAX] : Symbol[])...); on=:date)
        sort!(joined, :dowy)

        # Convert Q to Float64 with NaN for missing. A non-finite value (NaN, ±Inf) or a
        # non-Real entry is missing (review 2026-10-12: ±Inf is not an observation)
        Q = Vector{Float64}(undef, nrow(joined))
        for i in 1:nrow(joined)
            v = joined[i, :Q]
            if ismissing(v) || !(v isa Real) || !isfinite(v)
                Q[i] = NaN
            else
                Q[i] = Float64(v)
            end
        end

        # --- (b) Raw diagnostics ---
        na_count = count(isnan, Q)
        # Max consecutive NaN run
        max_na_run = 0
        current_run = 0
        for i in 1:length(Q)
            if isnan(Q[i])
                current_run += 1
                max_na_run = max(max_na_run, current_run)
            else
                current_run = 0
            end
        end

        # Seasonal completeness from RAW data
        raw_months = joined.month
        win_obs = count(i -> raw_months[i] in season_months["win"] && !isnan(Q[i]), 1:length(Q))
        spr_obs = count(i -> raw_months[i] in season_months["spr"] && !isnan(Q[i]), 1:length(Q))
        sum_obs = count(i -> raw_months[i] in season_months["sum"] && !isnan(Q[i]), 1:length(Q))
        fal_obs = count(i -> raw_months[i] in season_months["fal"] && !isnan(Q[i]), 1:length(Q))
        win_expected = count(m -> m in season_months["win"], raw_months)
        spr_expected = count(m -> m in season_months["spr"], raw_months)
        sum_expected = count(m -> m in season_months["sum"], raw_months)
        fal_expected = count(m -> m in season_months["fal"], raw_months)
        win_complete = win_expected > 0 && (win_obs / win_expected) >= seasonal_min_fraction
        spr_complete = spr_expected > 0 && (spr_obs / spr_expected) >= seasonal_min_fraction
        sum_complete = sum_expected > 0 && (sum_obs / sum_expected) >= seasonal_min_fraction
        fal_complete = fal_expected > 0 && (fal_obs / fal_expected) >= seasonal_min_fraction

        # Negative Q flag
        negative_flag = any(q -> !isnan(q) && q < 0.0, Q)

        # Constant-SD flag: check if any month has suspiciously constant flow
        constant_sd_flag = false
        if constant_sd_enabled
            for m in 1:12
                m_mask = raw_months .== m
                m_vals = Q[m_mask]
                m_nonzero = filter(v -> !isnan(v) && v > 0.0, m_vals)
                if length(m_nonzero) >= constant_sd_min_days
                    if length(unique(m_nonzero)) <= constant_sd_max_unique
                        constant_sd_flag = true
                        break
                    end
                end
            end
        end

        # --- (c) Reject years ---
        rejected = false
        if na_count > max_raw_na
            push!(rejected_rows, (water_year=yr, reason="too_many_na ($na_count > $max_raw_na)"))
            rejected = true
        end
        if !rejected && max_na_run > max_gap_days
            push!(rejected_rows, (water_year=yr, reason="gap_too_long ($max_na_run > $max_gap_days)"))
            rejected = true
        end
        if !rejected && reject_negative && negative_flag
            push!(rejected_rows, (water_year=yr, reason="negative_flow"))
            rejected = true
        end
        # constant_sd_flag is tracked in diagnostics but no longer causes rejection

        interpolated_count = 0

        if !rejected
            # --- (d) Interpolate internal gaps <= max_gap_days ---
            i = 1
            n = length(Q)
            while i <= n
                if isnan(Q[i])
                    # Find end of gap
                    gap_start = i
                    gap_end = i
                    while gap_end < n && isnan(Q[gap_end + 1])
                        gap_end += 1
                    end
                    gap_len = gap_end - gap_start + 1

                    # Check boundaries for internal gap
                    if gap_start > 1 && gap_end < n && gap_len <= max_gap_days
                        left_val = Q[gap_start - 1]
                        right_val = Q[gap_end + 1]
                        if !isnan(left_val) && !isnan(right_val)
                            span = (gap_end + 1) - (gap_start - 1)
                            for k in gap_start:gap_end
                                Q[k] = left_val + (right_val - left_val) * (k - (gap_start - 1)) / span
                                interpolated_count += 1
                            end
                        end
                    end

                    i = gap_end + 1
                else
                    i += 1
                end
            end

            # --- (e) Post-interpolation: reject if residual NaNs ---
            residual_na = count(isnan, Q)
            if reject_residual_na && residual_na > 0
                push!(rejected_rows, (water_year=yr, reason="residual_na ($residual_na remaining)"))
                rejected = true
            end
        end

        # NA-cause tracking (ice-affected days from USGS qualifier flags)
        # KNOWN ISSUE (CHANGELOG → Known Issues, 2026-08-26): this guard can
        # never fire — `joined` is built above via select(yr_df, :date, :Q,
        # [:PPT], [:SWE]), which drops the `flag` column, so
        # `"flag" in names(joined)` is always false and
        # `ice_affected_days_total` is structurally 0 for every gage.
        # Do NOT "fix" this here without cross-language coordination: rpkg
        # deliberately matches this behavior for parity.
        na_cause_ice = 0
        na_cause_other = 0
        if na_count > 0 && "flag" in names(joined)
            na_mask = isnan.(Q)
            na_flags = joined.flag[na_mask]
            for f in na_flags
                if !ismissing(f) && occursin(r"[Ii]ce", string(f))
                    na_cause_ice += 1
                end
            end
            na_cause_other = na_count - na_cause_ice
        end

        # Record diagnostics
        push!(diag_rows, (water_year=yr, na_count=na_count, max_na_run=max_na_run, interpolated_count=interpolated_count, negative_flag=negative_flag, constant_sd_flag=constant_sd_flag, na_cause_ice=na_cause_ice, na_cause_other=na_cause_other))
        push!(seasonal_rows, (water_year=yr, win_complete=win_complete, spr_complete=spr_complete, sum_complete=sum_complete, fal_complete=fal_complete))

        if rejected
            continue
        end

        # Write cleaned Q back
        joined[!, :Q] = Q

        # --- (f) Handle PPT if present ---
        ppt_ok = true
        if has_ppt
            P = Vector{Float64}(undef, nrow(joined))
            for i in 1:nrow(joined)
                v = joined[i, :PPT]
                if ismissing(v) || !(v isa Real) || !isfinite(v)   # non-finite = missing
                    P[i] = NaN
                else
                    P[i] = Float64(v)
                end
            end

            ppt_na_count = count(isnan, P)
            if ppt_na_count > max_raw_na_ppt
                ppt_ok = false
            end

            if reject_negative_ppt && any(p -> !isnan(p) && p < 0.0, P)
                ppt_ok = false
            end

            if ppt_ok
                # Interpolate PPT gaps
                i = 1
                n = length(P)
                while i <= n
                    if isnan(P[i])
                        gap_start = i
                        gap_end = i
                        while gap_end < n && isnan(P[gap_end + 1])
                            gap_end += 1
                        end
                        gap_len = gap_end - gap_start + 1
                        if gap_start > 1 && gap_end < n && gap_len <= max_gap_ppt
                            left_val = P[gap_start - 1]
                            right_val = P[gap_end + 1]
                            if !isnan(left_val) && !isnan(right_val)
                                span = (gap_end + 1) - (gap_start - 1)
                                for k in gap_start:gap_end
                                    P[k] = left_val + (right_val - left_val) * (k - (gap_start - 1)) / span
                                end
                            end
                        end
                        i = gap_end + 1
                    else
                        i += 1
                    end
                end

                ppt_residual = count(isnan, P)
                if ppt_residual > 0
                    ppt_ok = false
                end
            end

            if ppt_ok
                joined[!, :PPT] = P
            end
        end

        # --- (f2) Handle SWE if present (mirrors the PPT policy) ---
        swe_ok = true
        if has_swe
            S = Vector{Float64}(undef, nrow(joined))
            for i in 1:nrow(joined)
                v = joined[i, :SWE]
                if ismissing(v) || !(v isa Real) || !isfinite(v)   # non-finite = missing
                    S[i] = NaN
                else
                    S[i] = Float64(v)
                end
            end

            swe_na_count = count(isnan, S)
            if swe_na_count > max_raw_na_swe
                swe_ok = false
            end

            if reject_negative_swe && any(s -> !isnan(s) && s < 0.0, S)
                swe_ok = false
            end

            if swe_ok
                # Interpolate SWE gaps (internal, <= max_gap_swe days)
                i = 1
                n = length(S)
                while i <= n
                    if isnan(S[i])
                        gap_start = i
                        gap_end = i
                        while gap_end < n && isnan(S[gap_end + 1])
                            gap_end += 1
                        end
                        gap_len = gap_end - gap_start + 1
                        if gap_start > 1 && gap_end < n && gap_len <= max_gap_swe
                            left_val = S[gap_start - 1]
                            right_val = S[gap_end + 1]
                            if !isnan(left_val) && !isnan(right_val)
                                span = (gap_end + 1) - (gap_start - 1)
                                for k in gap_start:gap_end
                                    S[k] = left_val + (right_val - left_val) * (k - (gap_start - 1)) / span
                                end
                            end
                        end
                        i = gap_end + 1
                    else
                        i += 1
                    end
                end

                swe_residual = count(isnan, S)
                if swe_residual > 0
                    swe_ok = false
                end
            end

            if swe_ok
                joined[!, :SWE] = S
            end
        end

        # --- (f3) Handle TMIN/TMAX if present (2026-10-12): the SWE policy without
        # the negative-value check (temperatures are signed). The year is
        # temperature-valid only when BOTH series pass; both are written back
        # together so a frame never carries one filled and one raw series.
        temp_ok = true
        if has_temp
            filled = Dict{Symbol, Vector{Float64}}()
            for col in (:TMIN, :TMAX)
                T = Vector{Float64}(undef, nrow(joined))
                for i in 1:nrow(joined)
                    v = joined[i, col]
                    # non-finite values (NaN, ±Inf, missing) are all "missing" here —
                    # an Inf neighbour must not be interpolated from (review note 7)
                    T[i] = (ismissing(v) || !(v isa Real) || !isfinite(v)) ? NaN : Float64(v)
                end
                if count(isnan, T) > max_raw_na_temp
                    temp_ok = false
                else
                    _interpolate_internal_gaps!(T, max_gap_temp)
                    if count(isnan, T) > 0
                        temp_ok = false
                    end
                end
                filled[col] = T
            end
            if temp_ok
                joined[!, :TMIN] = filled[:TMIN]
                joined[!, :TMAX] = filled[:TMAX]
            end
        end

        push!(cleaned_dfs, joined)
        push!(valid_years, yr)
        if has_ppt && ppt_ok
            push!(valid_climate_years, yr)
        end
        if has_swe && swe_ok
            push!(valid_swe_years, yr)
        end
        if has_temp && temp_ok
            push!(valid_temp_years, yr)
        end
    end

    # --- (g) Build output ---
    if isempty(cleaned_dfs)
        out_data = DataFrame(date=Date[], Q=Float64[], water_year=Int[], month=Int[], dowy=Int[])
        if has_ppt
            out_data[!, :PPT] = Float64[]
        end
        if has_swe
            out_data[!, :SWE] = Float64[]
        end
        if has_temp
            out_data[!, :TMIN] = Float64[]
            out_data[!, :TMAX] = Float64[]
        end
    else
        out_data = vcat(cleaned_dfs...)
    end

    rejected_df = if isempty(rejected_rows)
        DataFrame(water_year=Int[], reason=String[])
    else
        DataFrame(rejected_rows)
    end

    seasonal_df = if isempty(seasonal_rows)
        DataFrame(water_year=Int[], win_complete=Bool[], spr_complete=Bool[], sum_complete=Bool[], fal_complete=Bool[])
    else
        DataFrame(seasonal_rows)
    end

    diagnostics_df = if isempty(diag_rows)
        DataFrame(water_year=Int[], na_count=Int[], max_na_run=Int[], interpolated_count=Int[], negative_flag=Bool[], constant_sd_flag=Bool[], na_cause_ice=Int[], na_cause_other=Int[])
    else
        DataFrame(diag_rows)
    end

    return (
        data = out_data,
        valid_years = valid_years,
        valid_climate_years = valid_climate_years,
        valid_swe_years = valid_swe_years,
        valid_temp_years = valid_temp_years,
        rejected_years = rejected_df,
        seasonal_flags = seasonal_df,
        diagnostics = diagnostics_df
    )
end


"""
    filter_qualifying_years(gage_data::DataFrame; kwargs...) -> (Vector{Int}, Bool)

Filter water years per-gage matching R's process_signatures_from_parquet().

Three-stage per-year filtering:
1. Per water year, check at least min_days_above days with Q > min_q_value
2. Per water year, check data completeness (>= min_frac_good of expected days)
3. Gage qualifies if >= min_num_years pass both sub-checks

Returns (qualifying_years, gage_qualifies) tuple.
"""
function filter_qualifying_years(gage_data::DataFrame;
        min_q_value::Real=CFG_MIN_Q_VALUE,
        min_days_above::Int=CFG_MIN_DAYS_ABOVE_THRESHOLD,
        min_frac_good::Real=CFG_MIN_FRAC_GOOD_DATA,
        min_num_years::Int=CFG_MIN_NUM_YEARS)
    qualifying = Int[]
    for wy in unique(gage_data.water_year)
        yr_data = gage_data[gage_data.water_year .== wy, :]
        q = yr_data.Q

        # Count non-NA values (handle both missing and NaN)
        n_nona = sum(x -> !ismissing(x) && (x isa Number ? !isnan(x) : true), q)

        # Sub-check 1: days above threshold
        n_above = sum(x -> !ismissing(x) && (x isa Number ? (!isnan(x) && x > min_q_value) : false), q)
        n_above < min_days_above && continue

        # Sub-check 2: data completeness (accounting for leap years)
        # Water year Y spans Oct 1 (Y-1) to Sep 30 (Y)
        expected_days = isleapyear(wy) ? 366 : 365
        min_good = floor(Int, expected_days * min_frac_good)
        n_nona < min_good && continue

        push!(qualifying, wy)
    end
    return qualifying, length(qualifying) >= min_num_years
end

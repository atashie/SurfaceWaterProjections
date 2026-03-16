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
to standard names: Date->date, site_id->gage_id, prcp->PPT.
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

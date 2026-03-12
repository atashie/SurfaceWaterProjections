"""
I/O functions for streamflow signature analysis.

Functions for reading/writing data and data validation.
"""

using Arrow
using CSV
using DataFrames
using Dates

# Required columns for streamflow data
const REQUIRED_COLS = ["gage_id", "date", "Q"]
const TEMPORAL_COLS = ["water_year", "month", "dowy"]


"""
    read_parquet(path::String) -> DataFrame

Read a parquet file into a DataFrame.
"""
function read_parquet(path::String)
    if !isfile(path)
        error("File not found: $path")
    end
    return DataFrame(Arrow.Table(path))
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

    if !(date_col in names(result))
        error("Date column '$date_col' not found")
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
    filter_valid_years(df::DataFrame; min_days=250, max_na_frac=0.1) -> DataFrame

Filter to water years with sufficient data.

Parameters
----------
df : DataFrame
    Input data with water_year and Q columns
min_days : Int
    Minimum days per water year
max_na_frac : Float64
    Maximum fraction of NA values allowed

Returns
-------
DataFrame
    Filtered data with only valid years
"""
function filter_valid_years(
    df::DataFrame;
    min_days::Int=250,
    max_na_frac::Float64=0.1
)
    # Group by water year and calculate stats
    year_stats = combine(
        groupby(df, :water_year),
        :Q => length => :n_days,
        :Q => (x -> sum(ismissing.(x) .| isnan.(x)) / length(x)) => :na_frac
    )

    # Filter to valid years
    valid_years = year_stats[
        (year_stats.n_days .>= min_days) .& (year_stats.na_frac .<= max_na_frac),
        :water_year
    ]

    return df[in.(df.water_year, Ref(Set(valid_years))), :]
end

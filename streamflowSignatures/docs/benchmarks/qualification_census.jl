# Qualification census — gage-inclusion and year-rejection accounting for the two
# standard windows (WY 1993–2025, WY 1980–2025) plus the full record, over ALL gages
# in the streamflow parquet (not only the ones that qualified).
#
# Runs the canonical preprocessor (`preprocess_daily_data`) on every gage and
# re-implements the benchmark runner's inclusion math (window-capped valid years,
# frac >= 0.60 of the window-start-anchored possible years, >= 20 valid years).
# Written 2026-09-04 for the manuscript filtering-stage census
# (docs/plans/2026-09-04-filtering-stage-census.md). Reproduced both standard
# products' gage sets exactly (6,678 / 6,250; zero mismatches).
#
# Usage:
#   STREAMFLOW_DATA_PATH=/path/to/combined_streamflow_data.parquet \
#     julia --project=julia docs/benchmarks/qualification_census.jl <out.csv>
# Writes <out.csv> (one row per gage × window: valid/rejected-year counts by reason,
# fraction, gate outcomes) and <out>_years.csv (one row per gage × water year with
# the preprocessor's status). Summarize with summarize_qualification_census.py.
# Runtime ~3 min on the M1 (8,014 gages).
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "julia"))
using StreamflowSignatures, DataFrames, CSV, Dates

const SF = get(ENV, "STREAMFLOW_DATA_PATH", "/Volumes/Untitled/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet")
const OUT = length(ARGS) >= 1 ? ARGS[1] : "qualification_census.csv"
const WINDOWS = [(1993, 2025), (1980, 2025), (typemin(Int), typemax(Int))]
const MIN_YEARS = 20
const MIN_FRAC = 0.60

println("Loading parquet..."); flush(stdout)
t0 = time()
sf = read_parquet(SF)
sf.gage_id = string.(sf.gage_id)
println("  rows=$(nrow(sf)) cols=$(names(sf)) in $(round(time()-t0, digits=1))s"); flush(stdout)

reason_key(r) = startswith(r, "no_data") ? "no_data" :
                startswith(r, "too_many_na") ? "too_many_na" :
                startswith(r, "gap_too_long") ? "gap_too_long" :
                startswith(r, "negative_flow") ? "negative_flow" :
                startswith(r, "residual_na") ? "residual_na" : "other"

rows = NamedTuple[]
yearrows = NamedTuple[]   # per gage-year status for the whole record
g = groupby(sf, :gage_id)
n = length(g); i = 0
for (key, gv) in pairs(g)
    global i += 1
    i % 500 == 0 && (println("  gage $i/$n  ($(round(time()-t0, digits=0))s)"); flush(stdout))
    gdf = DataFrame(gv)
    res = preprocess_daily_data(gdf)
    all_wy = sort(unique(Int.(round.(gdf.water_year))))
    valid = Set(res.valid_years)
    rej = Dict(Int(r.water_year) => reason_key(String(r.reason)) for r in eachrow(res.rejected_years))
    for yr in all_wy
        nq = count(!ismissing, gdf.Q[Int.(round.(gdf.water_year)) .== yr])
        push!(yearrows, (gage_id=key.gage_id, water_year=yr, n_rows=sum(Int.(round.(gdf.water_year)) .== yr),
                         status = yr in valid ? "valid" : get(rej, yr, "unlisted")))
    end
    for (lo, hi) in WINDOWS
        wname = lo == typemin(Int) ? "full" : "$(lo)-$(hi)"
        wy_in = filter(y -> lo <= y <= hi, all_wy)
        vin = filter(y -> lo <= y <= hi, res.valid_years)
        rin = [get(rej, y, "") for y in wy_in if !(y in valid)]
        cnt(k) = count(==(k), rin)
        if isempty(wy_in)
            total_possible = 0; frac = NaN
        else
            start_yr = lo == typemin(Int) ? minimum(all_wy) : lo
            total_possible = maximum(wy_in) - start_yr + 1
            frac = length(vin) / total_possible
        end
        fails_frac = !isempty(wy_in) && frac < MIN_FRAC
        fails_min = length(vin) < MIN_YEARS
        included = !isempty(wy_in) && !fails_frac && !fails_min
        push!(rows, (gage_id=key.gage_id, window=wname, n_wy_with_rows=length(wy_in), n_valid=length(vin),
                     n_rejected=length(wy_in)-length(vin),
                     rej_no_data=cnt("no_data"), rej_too_many_na=cnt("too_many_na"), rej_gap=cnt("gap_too_long"),
                     rej_negative=cnt("negative_flow"), rej_residual=cnt("residual_na"), rej_unlisted=cnt(""),
                     first_wy = isempty(wy_in) ? 0 : minimum(wy_in), last_wy = isempty(wy_in) ? 0 : maximum(wy_in),
                     first_valid = isempty(vin) ? 0 : minimum(vin), last_valid = isempty(vin) ? 0 : maximum(vin),
                     total_possible=total_possible, frac=frac, fails_frac=fails_frac, fails_min_years=fails_min, included=included))
    end
end
CSV.write(OUT, DataFrame(rows))
CSV.write(replace(OUT, ".csv" => "_years.csv"), DataFrame(yearrows))
println("Wrote $(length(rows)) gage-window rows and $(length(yearrows)) gage-year rows in $(round(time()-t0, digits=0))s"); flush(stdout)

# Independent qualification audit for a production benchmark run (gate 9).
# Repo-resident successor of the (drive-lost) per-run audit scripts, July 2026.
#
# Reimplements the benchmark's window/fraction inclusion math (START_YEAR <= wy <=
# END_YEAR (optional cap, July 2026), frac >= MIN_FRAC of possible in-window years,
# >= MIN_YEARS valid years) from preprocess_daily_data outputs — NOT via the
# benchmark's code path — for a stratified sample of edge gages (near-threshold
# included, early-ending, late-starting, and excluded-by-filter), and cross-checks
# the run's actual include/exclude decisions and per-gage year columns.
#
# Usage:
#   julia audit_qualification.jl <run_signatures.csv> <reference_superset.csv> \
#         <start_year> [streamflow_parquet] [min_years] [min_frac] [end_year]
# The reference superset (any run/file whose gage list is a superset, e.g. a
# broader-window run) supplies candidate EXCLUDED gages. end_year (7th arg, added
# 2026-07-22 for the WY 1993-2025 standard run) caps both the valid-year window and
# the possible-years denominator, mirroring run_julia_benchmark.jl; omitted = no cap.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "julia"))
using StreamflowSignatures
using DataFrames, CSV, Dates

const PROD = ARGS[1]
const REFERENCE = ARGS[2]
const START_YEAR = parse(Int, ARGS[3])
const SF = length(ARGS) >= 4 ? ARGS[4] :
    raw"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet"
const MIN_YEARS = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 20
const MIN_FRAC = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 0.60
const END_YEAR = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : typemax(Int)

prod = CSV.read(PROD, DataFrame;
                select=["gage_id", "num_water_years", "start_water_year", "end_water_year"],
                types=Dict("gage_id" => String))
ref = CSV.read(REFERENCE, DataFrame; select=["gage_id"], types=Dict("gage_id" => String))
prod_set = Set(prod.gage_id)
dropped = sort(collect(setdiff(Set(ref.gage_id), prod_set)))

psort = sort(prod, :num_water_years)
edge_included = vcat(
    first(psort.gage_id, 4),                                                     # near the year floor
    first(psort.gage_id[psort.end_water_year .< maximum(psort.end_water_year)], 3),  # early-ending
    first(psort.gage_id[psort.start_water_year .> START_YEAR], 3),               # late-starting
)
edge_excluded = dropped[1:min(4, length(dropped))]
sample_ids = unique(vcat(edge_included, edge_excluded))

println("Loading streamflow parquet...")
sf = read_parquet(SF)
if !("water_year" in names(sf))
    sf = add_water_year_columns(sf)
end
sf.gage_id = string.(sf.gage_id)

println("\ngage_id | in_prod | n_valid_win | total_possible | frac | independent_decision | verdict")
fails = 0
for g in sample_ids
    gdf = sf[sf.gage_id .== g, :]
    if isempty(gdf)
        println("$g | NOT IN PARQUET — cannot audit")
        global fails += 1
        continue
    end
    pp = preprocess_daily_data(gdf)
    # --- independent reimplementation of the inclusion rule ---
    valid_win = [y for y in pp.valid_years if START_YEAR <= y <= END_YEAR]
    raw_in_range = [y for y in unique(gdf.water_year) if START_YEAR <= y <= END_YEAR]
    total_possible = isempty(raw_in_range) ? 0 : maximum(raw_in_range) - START_YEAR + 1
    frac = total_possible == 0 ? 0.0 : length(valid_win) / total_possible
    decide = (total_possible > 0) && (frac >= MIN_FRAC) && (length(valid_win) >= MIN_YEARS)
    in_prod = g in prod_set
    match = decide == in_prod
    detail = ""
    if in_prod
        row = prod[prod.gage_id .== g, :][1, :]
        col_ok = (row.num_water_years == length(valid_win)) &&
                 (row.start_water_year == minimum(valid_win)) &&
                 (row.end_water_year == maximum(valid_win))
        match = match && col_ok
        detail = " | cols n=$(row.num_water_years)/$(length(valid_win)) wy=$(row.start_water_year)-$(row.end_water_year) vs $(minimum(valid_win))-$(maximum(valid_win))"
    end
    if !match
        global fails += 1
    end
    println("$g | $in_prod | $(length(valid_win)) | $total_possible | $(round(frac, digits=3)) | $decide | $(match ? "OK" : "MISMATCH")$detail")
end

println(fails == 0 ? "\nAUDIT PASS ($(length(sample_ids)) gages)" : "\nAUDIT FAIL ($fails mismatches)")
exit(fails == 0 ? 0 : 1)

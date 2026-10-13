# Cross-language fixture for the climate metrics (2026-10-12): writes a deterministic
# 30-year daily input (dates, Q, PPT, TMIN, TMAX) and canonical Julia's 39 x 16 summary
# values for it, so the Python and rpkg ports can be checked against bit-identical input.
#
#     julia --project=julia docs/benchmarks/crosscheck_climate_fixture.jl OUT_DIR
#
# Writes OUT_DIR/climate_fixture_input.parquet (+ .csv for reading),
# OUT_DIR/climate_fixture_julia.json and OUT_DIR/climate_fixture_julia_annual.csv. Then:
#     python docs/benchmarks/crosscheck_climate_fixture.py OUT_DIR        # Python vs Julia
#     Rscript docs/benchmarks/crosscheck_climate_fixture.R OUT_DIR        # rpkg vs Julia
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "julia"))
using StreamflowSignatures, DataFrames, Dates, CSV, JSON, Parquet2

out = length(ARGS) >= 1 ? ARGS[1] : @__DIR__
mkpath(out)

frames = DataFrame[]
for wy in 1991:2020
    dates = collect(Date(wy - 1, 10, 1):Day(1):Date(wy, 9, 30))
    n = length(dates)
    d = collect(1:n)
    mon = month.(dates)
    # deterministic, no transcendental functions (bit-identical across languages)
    P = [(k % 7 == 0) ? 0.0 : Float64(mon[k]) + 0.01 * k + 0.001 * ((k * 37) % 11) for k in d]
    TX = [20.0 - 0.15 * abs(k - 100) + 0.02 * ((k * 13) % 17) for k in d]
    TN = TX .- 12.0 .- 0.01 .* [(k * 7) % 5 for k in d]
    Q = [1.0 + 0.001 * k for k in d]
    push!(frames, DataFrame(gage_id = fill("FIX1", n), date = dates, Q = Q, PPT = P, TMIN = TN, TMAX = TX))
end
raw = vcat(frames...)
# Parquet carries the doubles exactly; the CSV is for reading. (R's read.csv and pandas'
# default converter parse "10.014" one ulp off, which broke bit-identity in 2026-10-12 tests.)
Parquet2.writefile(joinpath(out, "climate_fixture_input.parquet"), raw)
CSV.write(joinpath(out, "climate_fixture_input.csv"), raw)

df = add_water_year_columns(raw)
pp = preprocess_daily_data(df)
cdata = pp.data[in.(pp.data.water_year, Ref(Set(pp.valid_climate_years))), :]
tdata = pp.data[in.(pp.data.water_year, Ref(Set(pp.valid_temp_years))), :]
cp = (start_year = CFG_CP_START_WATER_YEAR, end_year = CFG_CP_END_WATER_YEAR,
      min_total_obs = CFG_CP_MIN_TOTAL_OBS, min_segment_obs = CFG_CP_MIN_SEGMENT_OBS)
kw = (trend_completeness = CFG_NA_TREND_MIN_FRACTION, decade_completeness = CFG_NA_DECADE_MIN_FRACTION,
      min_values_for_stats = CFG_MIN_VALUES_FOR_STATS, changepoint = cp)
coll = AnnualCollector()
res = Dict{String, Any}()
merge!(res, calculate_precipitation_metrics(cdata; kw..., collector=coll))
merge!(res, calculate_temperature_metrics(tdata; kw..., collector=coll))
# annual values (signature, water_year, value) — the bit-identity check across languages
annual = DataFrame(signature = String.(coll.signature), water_year = Int64.(coll.water_year), value = coll.value)
Parquet2.writefile(joinpath(out, "climate_fixture_julia_annual.parquet"), annual)  # exact doubles
CSV.write(joinpath(out, "climate_fixture_julia_annual.csv"), annual)              # for reading
# NaN is not JSON: encode as the string "NaN"
enc = Dict(k => (isnan(v) ? "NaN" : v) for (k, v) in res)
open(joinpath(out, "climate_fixture_julia.json"), "w") do f
    JSON.print(f, Dict("n_values" => length(enc), "valid_temp_years" => pp.valid_temp_years,
                       "valid_climate_years" => pp.valid_climate_years, "values" => enc), 1)
end
println("wrote $(length(enc)) Julia summary values and $(nrow(annual)) annual values to $out")

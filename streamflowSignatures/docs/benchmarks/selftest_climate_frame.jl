# selftest_climate_frame.jl — fixture self-test of the climate-frame helper (climate_frame.jl;
# plan §10.2 tests (i)–(iv), 2026-10-26) and the Julia reference for the Python and rpkg
# self-tests (selftest_climate_frame.py / .R read the parquets this writes).
#
#     julia --project=julia docs/benchmarks/selftest_climate_frame.jl OUT_DIR
#
# One synthetic gage CF1, water years 1991–2020, Q / PPT / SWE / TMIN / TMAX deterministic
# (no transcendental functions — bit-identical across languages). Daymet-style climate on
# the 365-day calendar (Dec 31 of leap years absent). Streamflow rows absent on:
#   WY 1995  five storm days (dowy 45, 135, 243 isolated; 198–199 a 2-day gap)     -> (i)
#   WY 2005  the whole water year (an interior year with no streamflow at all)     -> (ii)
#   2000-12-30, the neighbour of the Daymet-missing 2000-12-31                      -> (iii)
#   WY 2009  thirty isolated days + the Daymet-missing 2008-12-31 (31 climate NaN)   -> (iii)
#   WY 2013  2012-12-29, 12-30 and 2013-01-01 around the Daymet-missing 2012-12-31   -> (iii)
#   the record starts 1990-10-15 and ends 2020-03-15 (window inside the span)        -> (ii)
# Checks: the outer frame adds exactly the 39 absent dates inside the span and the gage's
# water years (never WY 2005, never outside the record); every Q-only preprocessor output
# and every Q-only signature is bit-identical between modes; prcp_annual in outer mode
# equals the Daymet water-year sum (Dec 31 rule) and in left mode differs by the
# interpolation residual of the absent days; the Dec 31 interpolate moves by the expected
# amount; WY 2009 and 2013 are climate-, temperature- and SWE-valid only in outer mode.
# Writes OUT_DIR/climate_frame_fixture_streamflow.parquet, _climate.parquet, and per mode
# climate_frame_fixture_annual_<mode>.parquet + climate_frame_fixture_summary_<mode>.json.
# Exit status 1 on any failed check.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "julia"))
using StreamflowSignatures, DataFrames, Dates, CSV, JSON, Parquet2
include(joinpath(@__DIR__, "climate_frame.jl"))

out = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "climate_frame_fixture")
mkpath(out)
n_pass = Ref(0); n_fail = Ref(0)
function check(label, cond)
    if cond; n_pass[] += 1; println("  PASS: $label") else n_fail[] += 1; println("  FAIL: $label") end
end
same(a, b) = (isnan(a) && isnan(b)) || a === b || a == b
close_to(a, b; atol=1e-9, rtol=1e-9) = abs(a - b) <= atol + rtol * max(abs(a), abs(b))

# ---------------- fixture ----------------
const GID = "CF1"
wy_of(d) = month(d) >= 10 ? year(d) + 1 : year(d)
dowy_of(d) = Dates.value(d - Date(wy_of(d) - 1, 10, 1)) + 1
# dry every 7th day, a 5 mm storm every 9th day (so an absent storm day is NOT reproduced by
# the linear interpolate of its neighbours), otherwise a slowly varying base scaled by a
# year factor (so annual P differs between years and the Q–P families are well conditioned)
pval(d) = (k = dowy_of(d); w = wy_of(d); (k % 7 == 0) ? 0.0 :
           (Float64(month(d)) + 0.01 * k + 0.001 * ((k * 37) % 11)) * (1.0 + 0.05 * ((w * 5) % 11)) + (k % 9 == 0 ? 5.0 : 0.0))
txval(d) = (k = dowy_of(d); 20.0 - 0.15 * abs(k - 100) + 0.02 * ((k * 13) % 17))
tnval(d) = (k = dowy_of(d); txval(d) - 12.0 - 0.01 * ((k * 7) % 5))
sweval(d) = (k = dowy_of(d); max(0.0, 40.0 - 0.4 * abs(k - 120)))
# Q varies between years and jumps on storm days, so no annual series is constant across
# years and no year's mean sits exactly on an observed value (rank statistics and the
# TQmean threshold would otherwise amplify last-bit differences between languages)
qval(d) = (k = dowy_of(d); w = wy_of(d); 1.0 + 0.001 * k + 0.1 * ((w * 7) % 13) + 0.017 * (w - 1990) +
           (k % 9 == 0 ? 0.5 * (1 + 0.01 * ((w * 3) % 7)) : 0.0) + 0.013 * ((k * 5) % 7))

all_dates = collect(Date(1990, 10, 1):Day(1):Date(2020, 9, 30))
record_start, record_end = Date(1990, 10, 15), Date(2020, 3, 15)
clim_dates = filter(d -> !(month(d) == 12 && day(d) == 31 && isleapyear(d)), all_dates)
wyd(wy, k) = Date(wy - 1, 10, 1) + Day(k - 1)
absent = Set{Date}()
for k in (45, 135, 243, 198, 199); push!(absent, wyd(1995, k)); end
for d in Date(2004, 10, 1):Day(1):Date(2005, 9, 30); push!(absent, d); end
push!(absent, Date(2000, 12, 30))
for k in 10:10:300; push!(absent, wyd(2009, k)); end
for d in (Date(2012, 12, 29), Date(2012, 12, 30), Date(2013, 1, 1)); push!(absent, d); end
sf_dates = filter(d -> record_start <= d <= record_end && !(d in absent), all_dates)
expected_added = filter(d -> record_start <= d <= record_end && d in absent && wy_of(d) != 2005 && d in Set(clim_dates), all_dates)

streamflow = DataFrame(gage_id = fill(GID, length(sf_dates)), date = sf_dates, Q = qval.(sf_dates),
                       year = year.(sf_dates), doy = dayofyear.(sf_dates), flag = fill("A", length(sf_dates)))
streamflow = add_water_year_columns(streamflow)
climate = DataFrame(gage_id = fill(GID, length(clim_dates)), date = clim_dates, PPT = pval.(clim_dates),
                    SWE = sweval.(clim_dates), TMIN = tnval.(clim_dates), TMAX = txval.(clim_dates))
Parquet2.writefile(joinpath(out, "climate_frame_fixture_streamflow.parquet"), streamflow)
Parquet2.writefile(joinpath(out, "climate_frame_fixture_climate.parquet"), climate)
println("fixture: $(nrow(streamflow)) streamflow rows, $(nrow(climate)) climate rows, $(length(expected_added)) absent dates inside the span and the gage's water years")

# ---------------- frames and preprocessing ----------------
frames = Dict(m => frame_with_climate(streamflow, climate, m) for m in ("left", "outer"))
pp = Dict(m => preprocess_daily_data(frames[m]) for m in ("left", "outer"))
L, O = frames["left"], frames["outer"]

println("\n(ii) frame composition")
check("left frame = the streamflow rows", nrow(L) == nrow(streamflow) && Set(L.date) == Set(streamflow.date))
added = sort(setdiff(Date.(O.date), Date.(L.date)))
check("outer adds exactly the $(length(expected_added)) absent dates inside the span and the gage's water years", added == expected_added)
check("no added row in the absent water year 2005", !(2005 in Set(O.water_year)))
check("no added row before the record start or after its end", all(record_start .<= Date.(O.date) .<= record_end))
check("the water-year set is identical between modes", Set(L.water_year) == Set(O.water_year))
check("Q and flag are missing on every added row", all(ismissing, O.Q[nrow(L)+1:end]) && all(ismissing, O.flag[nrow(L)+1:end]))
check("the joined rows of the outer frame equal the left frame", isequal(O[1:nrow(L), names(L)], L))
check("added rows carry the climate values of the input", all(isequal(pval(Date(r.date)), r.PPT) for r in eachrow(O[nrow(L)+1:end, :])))

println("\n(ii) helper edge cases (review 2026-10-26)")
let sf2 = streamflow[1:400, :], cl2 = climate[1:600, :]
    dup = vcat(sf2, sf2[200:200, :])                      # a duplicated streamflow date
    o_dup = frame_with_climate(dup, cl2, "outer"); l_dup = frame_with_climate(dup, cl2, "left")
    check("a duplicated streamflow date is kept as the left join keeps it and never added again",
          nrow(l_dup) == nrow(dup) && nrow(o_dup) - nrow(l_dup) == length(setdiff(Date.(cl2.date[(Date.(cl2.date) .>= minimum(sf2.date)) .& (Date.(cl2.date) .<= maximum(sf2.date))]), Date.(dup.date))) &&
          count(==(dup.date[200]), o_dup.date) == 2)
    check("zero climate rows: both modes return the streamflow rows with the climate columns missing",
          nrow(frame_with_climate(sf2, cl2[1:0, :], "outer")) == nrow(sf2) && all(ismissing, frame_with_climate(sf2, cl2[1:0, :], "outer").PPT))
    cl_nan = copy(cl2); cl_nan.PPT = allowmissing(cl_nan.PPT)
    i3 = findfirst(==(sf2.date[10]), cl_nan.date); i5 = findfirst(==(sf2.date[20]), cl_nan.date)
    cl_nan.PPT[i3] = missing; cl_nan.TMIN[i5] = NaN
    o_nan = frame_with_climate(sf2, cl_nan, "outer")
    check("missing / NaN climate values pass through unchanged (the preprocessor treats them)",
          ismissing(o_nan.PPT[findfirst(==(sf2.date[10]), o_nan.date)]) && isnan(o_nan.TMIN[findfirst(==(sf2.date[20]), o_nan.date)]))
    sf_dt = copy(sf2); sf_dt.date = DateTime.(sf_dt.date)
    check("a DateTime streamflow date column is refused in both modes",
          all(m -> (try; frame_with_climate(sf_dt, cl2, m); false; catch e; e isa ErrorException; end), ("left", "outer")))
    check("an unknown mode is refused", (try; frame_with_climate(sf2, cl2, "inner"); false; catch e; e isa ErrorException; end))
    check("inputs are not mutated", sf2 == streamflow[1:400, :] && cl2 == climate[1:600, :])
    check("left mode IS the former join expression (isequal, column order and types included)",
          isequal(frame_with_climate(sf2, cl2, "left"), leftjoin(sf2, cl2, on=[:gage_id, :date])) &&
          eltype.(eachcol(frame_with_climate(sf2, cl2, "left"))) == eltype.(eachcol(leftjoin(sf2, cl2, on=[:gage_id, :date]))))
    check("an empty streamflow frame returns the (empty) join in both modes",
          nrow(frame_with_climate(sf2[1:0, :], cl2, "outer")) == 0 && nrow(frame_with_climate(sf2[1:0, :], cl2, "left")) == 0)
    bnd = sf2[[1, 366], :]                                   # one row each at the two ends of one water year
    o_bnd = frame_with_climate(bnd, cl2, "outer")
    n_between = count(d -> bnd.date[1] <= d <= bnd.date[2] && wy_of(d) in Set(bnd.water_year), cl2.date) - 2
    check("a boundary-only year: every climate date strictly between the two rows is added ($(n_between)), nothing outside",
          nrow(o_bnd) == 2 + n_between && all(bnd.date[1] .<= o_bnd.date .<= bnd.date[2]))
    cl_dup = vcat(cl2, cl2[20:20, :])                        # an identical duplicate climate row on an absent date? make it absent first
    sf_hole = sf2[Not(findfirst(==(cl2.date[20]), sf2.date)), :]
    o_dup2 = frame_with_climate(sf_hole, cl_dup, "outer")
    check("an identical duplicate climate row on an absent date is added once",
          count(==(cl2.date[20]), o_dup2.date) == 1 && nrow(o_dup2) == nrow(frame_with_climate(sf_hole, cl2, "outer")))
    cl_conf = copy(cl_dup); cl_conf.PPT[end] = cl_conf.PPT[end] + 1.0
    check("a conflicting duplicate climate row on an absent date is an error",
          (try; frame_with_climate(sf_hole, cl_conf, "outer"); false; catch e; e isa ErrorException; end))
end

println("\n(ii) the runner's window denominator (run_julia_benchmark.jl: total_possible from unique(gage_df.water_year))")
function window_denominator(frame, p, start_wy, end_wy)
    all_wy = unique(frame.water_year)
    wy_in_range = filter(yr -> start_wy <= yr <= end_wy, all_wy)
    isempty(wy_in_range) && return (total_possible = 0, n_valid = 0, frac = NaN)
    total_possible = maximum(wy_in_range) - start_wy + 1
    n_valid = count(yr -> start_wy <= yr <= end_wy, p.valid_years)
    (total_possible = total_possible, n_valid = n_valid, frac = n_valid / total_possible)
end
WINDOWS = ((1993, 2005), (1993, 2006), (1980, 2025), (2000, 2004))
windows = Dict(m => Dict("$(w[1])-$(w[2])" => window_denominator(frames[m], pp[m], w...) for w in WINDOWS) for m in ("left", "outer"))
check("total_possible, the valid-year count and the fraction are identical between modes for every window",
      all(windows["left"][k] == windows["outer"][k] for k in keys(windows["left"])))
check("a window ending inside the absent WY 2005 (1993–2005) has total_possible 12 in both modes (2005 never enters the frame)",
      windows["left"]["1993-2005"].total_possible == 12 && windows["outer"]["1993-2005"].total_possible == 12)
check("a window ending after it (1993–2006) has total_possible 14 in both modes",
      windows["outer"]["1993-2006"].total_possible == 14)

println("\n(ii) Q-only preprocessor outputs identical between modes")
check("valid_years identical", pp["left"].valid_years == pp["outer"].valid_years)
check("rejected_years identical (WY 1991 and 2020 partial; 2005 absent from both)",
      isequal(pp["left"].rejected_years, pp["outer"].rejected_years) && !(2005 in pp["left"].rejected_years.water_year))
check("WY 1991 and WY 2020 rejected in both modes", all(y in pp["outer"].rejected_years.water_year for y in (1991, 2020)))
check("seasonal_flags identical", isequal(pp["left"].seasonal_flags, pp["outer"].seasonal_flags))
check("diagnostics identical", isequal(pp["left"].diagnostics, pp["outer"].diagnostics))
check("preprocessed Q identical (bitwise, NaN == NaN)",
      nrow(pp["left"].data) == nrow(pp["outer"].data) && all(same.(pp["left"].data.Q, pp["outer"].data.Q)) && pp["left"].data.date == pp["outer"].data.date)
check("the data frames hold the same water years", pp["left"].data.water_year == pp["outer"].data.water_year)

cp = (start_year = CFG_CP_START_WATER_YEAR, end_year = CFG_CP_END_WATER_YEAR,
      min_total_obs = CFG_CP_MIN_TOTAL_OBS, min_segment_obs = CFG_CP_MIN_SEGMENT_OBS)
kw = (trend_completeness = CFG_NA_TREND_MIN_FRACTION, decade_completeness = CFG_NA_DECADE_MIN_FRACTION,
      min_values_for_stats = CFG_MIN_VALUES_FOR_STATS, changepoint = cp)
qonly = Dict(m => calculate_all_signatures(pp[m].data, false; gage_id = GID, seasonal_flags = pp[m].seasonal_flags, kw...)
             for m in ("left", "outer"))
qkeys = sort(collect(keys(qonly["left"])))
check("Q-only orchestrator: the same $(length(qkeys)) keys in both modes", qkeys == sort(collect(keys(qonly["outer"]))))
check("Q-only orchestrator: every value bit-identical between modes",
      all(same(qonly["left"][k], qonly["outer"][k]) for k in qkeys))

println("\n(i) precipitation totals")
ppt_at(m, d) = (df = pp[m].data; i = findfirst(==(d), Date.(df.date)); i === nothing ? NaN : df.PPT[i])
daymet_sum(wy) = begin
    ds = filter(d -> wy_of(d) == wy, clim_dates)
    s = 0.0
    for d in sort(ds); s += pval(d); end
    d31 = Date(wy - 1, 12, 31)
    if isleapyear(d31)   # the Daymet-missing day, interpolated between Dec 30 and Jan 1
        s += (pval(Date(wy - 1, 12, 30)) + pval(Date(wy, 1, 1))) / 2
    end
    s
end
function interp_resid_wy1995()   # what the left frame replaces the five absent wet days with
    r = 0.0
    for k in (45, 135, 243)
        r += pval(wyd(1995, k)) - (pval(wyd(1995, k - 1)) + pval(wyd(1995, k + 1))) / 2
    end
    a, b = pval(wyd(1995, 197)), pval(wyd(1995, 200))
    r += pval(wyd(1995, 198)) - (a + (b - a) / 3)
    r += pval(wyd(1995, 199)) - (a + 2 * (b - a) / 3)
    return r
end
interp_resid_1995 = interp_resid_wy1995()
function run_full(p)   # the production orchestrator call of run_julia_benchmark.jl (new path)
    data = p.data
    cdata = data[in.(data.water_year, Ref(Set(p.valid_climate_years))), :]
    sdata = data[in.(data.water_year, Ref(Set(p.valid_swe_years))), :]
    tdata = data[in.(data.water_year, Ref(Set(p.valid_temp_years))), :]
    c = AnnualCollector()
    r = calculate_all_signatures(data, true; gage_id = GID, seasonal_flags = p.seasonal_flags,
                                 climate_data = cdata, snow_data = sdata, snow_climate_years = p.valid_climate_years,
                                 temp_data = tdata, collector = c, area_normalized = true, kw...)
    return c, r
end
coll = Dict{String, AnnualCollector}(); res = Dict{String, Dict{String, Any}}()
for m in ("left", "outer")
    coll[m], res[m] = run_full(pp[m])
end
annual(m) = DataFrame(signature = String.(coll[m].signature), water_year = Int64.(coll[m].water_year), value = coll[m].value)
ann = Dict(m => annual(m) for m in ("left", "outer"))
prcp(m, wy) = (a = ann[m]; i = findfirst((a.signature .== "prcp_annual") .& (a.water_year .== wy)); i === nothing ? NaN : a.value[i])
check("outer: prcp_annual equals the Daymet water-year sum (Dec 31 rule) on every climate-valid year",
      all(close_to(prcp("outer", wy), daymet_sum(wy)) for wy in pp["outer"].valid_climate_years))
check("left: prcp_annual WY 1995 = Daymet sum − the interpolation residual of the five absent wet days ($(round(interp_resid_1995, digits=4)) mm)",
      close_to(prcp("left", 1995), daymet_sum(1995) - interp_resid_1995) && !close_to(prcp("left", 1995), daymet_sum(1995)))
check("left and outer prcp_annual identical on every climate-valid year without an absent date",
      all(same(prcp("left", wy), prcp("outer", wy)) for wy in pp["left"].valid_climate_years if !(wy in (1995, 2001, 2009, 2013))))

println("\n(iii) Dec 31 neighbour and eligibility flips")
p29, p30, p1 = pval(Date(2000, 12, 29)), pval(Date(2000, 12, 30)), pval(Date(2001, 1, 1))
check("left: 2000-12-31 interpolated across the 2-day gap (Dec 30 absent row + Daymet-missing Dec 31)",
      close_to(ppt_at("left", Date(2000, 12, 31)), p29 + 2 * (p1 - p29) / 3; atol = 1e-12, rtol = 1e-12) &&
      close_to(ppt_at("left", Date(2000, 12, 30)), p29 + (p1 - p29) / 3; atol = 1e-12, rtol = 1e-12))
check("outer: 2000-12-30 restored and 2000-12-31 = mean(Dec 30, Jan 1)",
      ppt_at("outer", Date(2000, 12, 30)) == p30 &&
      close_to(ppt_at("outer", Date(2000, 12, 31)), (p30 + p1) / 2; atol = 1e-12, rtol = 1e-12))
check("WY 2009 (30 absent days + Daymet-missing Dec 31): climate-, temperature- and SWE-valid only in outer mode",
      !(2009 in pp["left"].valid_climate_years) && (2009 in pp["outer"].valid_climate_years) &&
      !(2009 in pp["left"].valid_temp_years) && (2009 in pp["outer"].valid_temp_years) &&
      !(2009 in pp["left"].valid_swe_years) && (2009 in pp["outer"].valid_swe_years))
check("WY 2013 (4-day climate gap around Dec 31): climate-valid only in outer mode; Q-valid in both",
      !(2013 in pp["left"].valid_climate_years) && (2013 in pp["outer"].valid_climate_years) &&
      (2013 in pp["left"].valid_years) && (2013 in pp["outer"].valid_years))
check("WY 2009 and 2013 are Q-valid in both modes", all(y in pp["left"].valid_years for y in (2009, 2013)))
check("every other climate-valid year is valid in both modes",
      setdiff(pp["outer"].valid_climate_years, pp["left"].valid_climate_years) == [2009, 2013] &&
      isempty(setdiff(pp["left"].valid_climate_years, pp["outer"].valid_climate_years)))

# ---------------- reference outputs for the ports ----------------
for m in ("left", "outer")
    Parquet2.writefile(joinpath(out, "climate_frame_fixture_annual_$m.parquet"), ann[m])
    enc = Dict(k => (v isa AbstractFloat && isnan(v) ? "NaN" : v) for (k, v) in res[m])
    open(joinpath(out, "climate_frame_fixture_summary_$m.json"), "w") do f
        JSON.print(f, Dict("mode" => m, "n_values" => length(enc), "values" => enc,
                           "valid_years" => pp[m].valid_years, "valid_climate_years" => pp[m].valid_climate_years,
                           "valid_swe_years" => pp[m].valid_swe_years, "valid_temp_years" => pp[m].valid_temp_years,
                           "n_frame_rows" => nrow(frames[m]), "n_added_rows" => nrow(frames[m]) - nrow(L),
                           "added_dates" => string.(sort(setdiff(Date.(frames[m].date), Date.(L.date)))),
                           "windows" => Dict(k => Dict("total_possible" => v.total_possible, "n_valid" => v.n_valid) for (k, v) in windows[m]),
                           # reference numbers for the Python and rpkg self-tests
                           "checks" => Dict("ppt_2000_12_30" => ppt_at(m, Date(2000, 12, 30)),
                                            "ppt_2000_12_31" => ppt_at(m, Date(2000, 12, 31)),
                                            "prcp_annual_1995" => prcp(m, 1995),
                                            "daymet_sum_1995" => daymet_sum(1995),
                                            "interp_resid_1995" => interp_resid_1995,
                                            "n_qonly_keys" => length(qkeys))), 1)
    end
end
println("\nwrote the fixture inputs and, per mode, $(nrow(ann["left"])) / $(nrow(ann["outer"])) annual values and $(length(res["left"])) summary values to $out")
println("\n$(n_pass[]) passed, $(n_fail[]) failed")
exit(n_fail[] == 0 ? 0 : 1)

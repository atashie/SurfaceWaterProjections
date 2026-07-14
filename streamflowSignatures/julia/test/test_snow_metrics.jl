# Tests for snow metrics (calculate_snow_metrics + SWE plumbing).
# Test map: docs/plans/snow_signatures_plan.md §7 (items 1-14, 16, 17; item 15 lives
# in smoke_test.jl). Synthetic gages are deterministic with hand-derived expectations.
#
# Conventions under test (plan §2): thresholded series SWE* (SWE >= 10 mm else 0),
# spells on SWE* > 0, anchor spell = spell containing the annual max (first-max ties),
# boundary censoring, melt increments m_t attributed to day t, SSM with 60-day
# seasonal spells, swe_max_to_ppt with the runoff-ratio PPT convention.

using Test
using DataFrames
using Dates
using Statistics
using Logging
using StreamflowSignatures

const SNOW_CP = (start_year=1980, end_year=2024, min_total_obs=20, min_segment_obs=10)

# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

"""Build a daily gage frame for `years`; `swe_fn(t)` maps dowy -> SWE (mm)."""
function snow_gage(years; swe_fn, ppt=2.0)
    dfs = DataFrame[]
    for yr in years
        dates = collect(Date(yr - 1, 10, 1):Day(1):Date(yr, 9, 30))
        n = length(dates)
        swe = Float64[swe_fn(t) for t in 1:n]
        ppt_vec = ppt isa Function ? Float64[ppt(t) for t in 1:n] : fill(Float64(ppt), n)
        push!(dfs, DataFrame(gage_id="snowtest", date=dates, Q=fill(1.0, n),
                             water_year=fill(yr, n), month=month.(dates),
                             dowy=collect(1:n), PPT=ppt_vec, SWE=swe))
    end
    return vcat(dfs...)
end

# Triangle: accumulate 2 mm/day dowy 32-152 (peak 242), melt 3 mm/day after.
triangle_swe(t) = 32 <= t <= 152 ? 2.0 * (t - 31) :
                  (t > 152 ? max(0.0, 242.0 - 3.0 * (t - 152)) : 0.0)
# Two peaks with a mid-winter thaw (valley stays >= 10 -> single spell).
twopeak_swe(t) = t < 32 ? 0.0 :
                 t <= 81 ? 2.0 * (t - 31) :
                 t <= 101 ? 100.0 - 4.0 * (t - 81) :
                 t <= 137 ? 20.0 + 5.0 * (t - 101) :
                 max(0.0, 200.0 - 4.0 * (t - 137))
# Three constant blocks: 90 d @ 50 mm, 15 d @ 30 mm, 15 d @ 20 mm.
blocks_swe(t) = (40 <= t <= 129) ? 50.0 : (160 <= t <= 174) ? 30.0 :
                (200 <= t <= 214) ? 20.0 : 0.0
# Five ephemeral 10-day spells at 15 mm.
ephemeral_swe(t) = any(s -> s <= t <= s + 9, (40, 70, 100, 130, 160)) ? 15.0 : 0.0
# Sub-threshold only (peak 8 mm < 10 mm).
trace_swe(t) = 100 <= t <= 149 ? 8.0 : 0.0
# Carryover: SWE = 50 from Oct 1, melts out mid-winter.
carryover_swe(t) = t <= 100 ? 50.0 : max(0.0, 50.0 - 5.0 * (t - 100))
# Persists: accumulates from dowy 201 through Sep 30.
persist_swe(t) = t <= 200 ? 0.0 : 2.0 * (t - 200)
# Snow-covered every day of the year.
fullyear_swe(t) = 50.0

collector_frame(coll) = DataFrame(signature=coll.signature, wy=coll.water_year,
                                  value=coll.value)

"""Single collector value for (signature, water_year); fails the test if not unique."""
function cval(adf::DataFrame, sig::String, wy::Int)
    rows = adf[(adf.signature .== sig) .& (adf.wy .== wy), :value]
    @test length(rows) == 1
    return rows[1]
end

snow_keys_of(d) = [k for k in keys(d)
                   if any(startswith(k, "$(m)_") for m in StreamflowSignatures.SNOW_METRICS)]

# ---------------------------------------------------------------------------
@testset "Snow metrics" begin

@testset "snow_spells helper" begin
    spells = StreamflowSignatures.snow_spells
    @test spells(Bool[]) == UnitRange{Int}[]
    @test spells([false, false]) == UnitRange{Int}[]
    @test spells([true, true, true]) == [1:3]
    @test spells([false, true, true, false, true]) == [2:3, 5:5]
    @test spells([true, false, true]) == [1:1, 3:3]
end

@testset "Triangle gage: exact per-year values (plan §7.1)" begin
    df = snow_gage([1999, 2001, 2004]; swe_fn=triangle_swe)  # 2004 = leap WY (366 d)
    coll = AnnualCollector()
    res = calculate_snow_metrics(df; valid_climate_years=[1999, 2001, 2004],
                                 collector=coll)
    @test length(res) == 14 * 8  # 112 keys without changepoint
    adf = collector_frame(coll)

    expected = Dict(
        "swe_max" => 242.0, "swe_max_dowy" => 152.0, "snow_cover_days" => 194.0,
        "snow_on_dowy" => 36.0, "snow_off_dowy" => 230.0, "melt_season_days" => 78.0,
        "melt_rate" => 242.0 / 78.0, "ssm" => 1.0, "swe_apr1" => 149.0,
        "melt_before_peak" => 0.0, "melt_before_peak_pct" => 0.0,
        "melt_before_peak_to_max_swe" => 0.0, "melt_com_dowy" => 193.0,
        "swe_max_to_ppt" => 242.0 / 730.0
    )
    for (sig, exp) in expected
        @test cval(adf, sig, 2001) ≈ exp atol=1e-12
    end

    # Leap year: April 1 2004 is dowy 184 -> SWE = 242 - 3*32 = 146; PPT total 732
    @test cval(adf, "swe_apr1", 2004) ≈ 146.0 atol=1e-12
    @test cval(adf, "swe_max_to_ppt", 2004) ≈ 242.0 / 732.0 atol=1e-12
    # Leap year has the same t-indexed geometry otherwise
    @test cval(adf, "swe_max", 2004) ≈ 242.0 atol=1e-12
    @test cval(adf, "snow_cover_days", 2004) ≈ 194.0 atol=1e-12
end

@testset "Two-peak gage: mid-winter melt accounting (plan §7.2, §7.16 m_t pin)" begin
    df = snow_gage([2001]; swe_fn=twopeak_swe)
    coll = AnnualCollector()
    calculate_snow_metrics(df; valid_climate_years=[2001], collector=coll)
    adf = collector_frame(coll)

    @test cval(adf, "swe_max", 2001) ≈ 200.0 atol=1e-12
    @test cval(adf, "swe_max_dowy", 2001) ≈ 137.0 atol=1e-12   # LARGER peak anchors
    @test cval(adf, "snow_on_dowy", 2001) ≈ 36.0 atol=1e-12
    @test cval(adf, "snow_off_dowy", 2001) ≈ 185.0 atol=1e-12
    @test cval(adf, "snow_cover_days", 2001) ≈ 149.0 atol=1e-12
    @test cval(adf, "ssm", 2001) ≈ 1.0 atol=1e-12              # one 149-d spell
    @test cval(adf, "melt_season_days", 2001) ≈ 48.0 atol=1e-12
    @test cval(adf, "melt_rate", 2001) ≈ 200.0 / 48.0 atol=1e-12
    # Thaw melt: 20 daily drops of 4 mm on days 82..101, all before the peak
    @test cval(adf, "melt_before_peak", 2001) ≈ 80.0 atol=1e-12
    @test cval(adf, "melt_before_peak_pct", 2001) ≈ 100.0 * 80.0 / 280.0 atol=1e-12
    @test cval(adf, "melt_before_peak_to_max_swe", 2001) ≈ 0.4 atol=1e-12
    # Cumulative melt reaches 140 (half of 280) on day 152
    @test cval(adf, "melt_com_dowy", 2001) ≈ 152.0 atol=1e-12
    @test cval(adf, "swe_apr1", 2001) ≈ 16.0 atol=1e-12        # 200 - 4*(183-137)
end

@testset "Spell arithmetic: SSM mixed / ephemeral (plan §7.3)" begin
    # 90 d + 15 d + 15 d -> SSM = (90 - 30) / 120
    df = snow_gage([2001]; swe_fn=blocks_swe)
    coll = AnnualCollector()
    calculate_snow_metrics(df; collector=coll)
    adf = collector_frame(coll)
    @test cval(adf, "ssm", 2001) ≈ 0.5 atol=1e-12
    @test cval(adf, "snow_cover_days", 2001) ≈ 120.0 atol=1e-12
    # Constant block: max plateau -> FIRST day (tie rule, plan §7.7)
    @test cval(adf, "swe_max_dowy", 2001) ≈ 40.0 atol=1e-12
    @test cval(adf, "snow_on_dowy", 2001) ≈ 40.0 atol=1e-12
    @test cval(adf, "snow_off_dowy", 2001) ≈ 130.0 atol=1e-12
    @test cval(adf, "melt_season_days", 2001) ≈ 90.0 atol=1e-12
    @test cval(adf, "melt_rate", 2001) ≈ 50.0 / 90.0 atol=1e-12
    # Block melt-outs land on days 130 / 175 / 215 (drop-to-zero attribution)
    @test cval(adf, "melt_com_dowy", 2001) ≈ 130.0 atol=1e-12  # 50 of 100 total

    # Five 10-day spells -> SSM = -1
    df2 = snow_gage([2001]; swe_fn=ephemeral_swe)
    coll2 = AnnualCollector()
    calculate_snow_metrics(df2; collector=coll2)
    adf2 = collector_frame(coll2)
    @test cval(adf2, "ssm", 2001) ≈ -1.0 atol=1e-12
    @test cval(adf2, "snow_cover_days", 2001) ≈ 50.0 atol=1e-12
    @test cval(adf2, "swe_max_dowy", 2001) ≈ 40.0 atol=1e-12   # first spell's first day
    @test cval(adf2, "snow_off_dowy", 2001) ≈ 50.0 atol=1e-12  # anchor spell only
    # Spell melt-outs: 15 mm on days 50/80/110/140/170; half of 75 reached on day 110
    @test cval(adf2, "melt_com_dowy", 2001) ≈ 110.0 atol=1e-12
end

@testset "Sub-threshold censoring: peak 8 mm is snow-free (plan §7.4)" begin
    df = snow_gage([1999, 2001, 2003]; swe_fn=trace_swe)
    coll = AnnualCollector()
    res = calculate_snow_metrics(df; valid_climate_years=[1999, 2001, 2003],
                                 collector=coll)
    adf = collector_frame(coll)
    for yr in (1999, 2001, 2003)
        @test cval(adf, "swe_max", yr) == 0.0
        @test cval(adf, "snow_cover_days", yr) == 0.0
        @test cval(adf, "swe_apr1", yr) == 0.0
        @test cval(adf, "swe_max_to_ppt", yr) == 0.0     # 0 / 730, PPT qualifies
        for sig in ("swe_max_dowy", "snow_on_dowy", "snow_off_dowy",
                    "melt_season_days", "melt_rate", "ssm", "melt_before_peak",
                    "melt_before_peak_pct", "melt_before_peak_to_max_swe",
                    "melt_com_dowy")
            @test isnan(cval(adf, sig, yr))
        end
    end
    # Dense zeros -> summary stats compute; timing stays NaN
    @test res["swe_max_mean"] == 0.0
    @test isnan(res["ssm_mean"])
end

@testset "Spell splitting at the threshold boundary (plan §7.5)" begin
    # One 9.9 mm day inside a deep spell splits it (30 d + 29 d, both ephemeral)
    split_a(t) = (100 <= t <= 159) ? (t == 130 ? 9.9 : 50.0) : 0.0
    dfa = snow_gage([2001]; swe_fn=split_a)
    colla = AnnualCollector()
    calculate_snow_metrics(dfa; collector=colla)
    adfa = collector_frame(colla)
    @test cval(adfa, "ssm", 2001) ≈ -1.0 atol=1e-12
    @test cval(adfa, "snow_cover_days", 2001) ≈ 59.0 atol=1e-12
    @test cval(adfa, "snow_off_dowy", 2001) ≈ 130.0 atol=1e-12  # anchor = first spell

    # A 10.0 mm day does NOT split (>= threshold): one 60-d spell -> seasonal
    split_b(t) = (100 <= t <= 159) ? (t == 130 ? 10.0 : 50.0) : 0.0
    dfb = snow_gage([2001]; swe_fn=split_b)
    collb = AnnualCollector()
    calculate_snow_metrics(dfb; collector=collb)
    adfb = collector_frame(collb)
    @test cval(adfb, "ssm", 2001) ≈ 1.0 atol=1e-12
    @test cval(adfb, "snow_cover_days", 2001) ≈ 60.0 atol=1e-12
    @test cval(adfb, "snow_off_dowy", 2001) ≈ 160.0 atol=1e-12
end

@testset "Boundary censoring: carryover / persistence / full-year (plan §7.6, §7.16)" begin
    # Snowpack present on Oct 1 -> snow_on censored; everything else computes
    df = snow_gage([2001]; swe_fn=carryover_swe)
    coll = AnnualCollector()
    calculate_snow_metrics(df; collector=coll)
    adf = collector_frame(coll)
    @test isnan(cval(adf, "snow_on_dowy", 2001))
    @test cval(adf, "swe_max_dowy", 2001) ≈ 1.0 atol=1e-12     # plateau tie -> day 1
    @test cval(adf, "snow_off_dowy", 2001) ≈ 109.0 atol=1e-12  # 50-5*(t-100) < 10 at 109
    @test cval(adf, "melt_season_days", 2001) ≈ 108.0 atol=1e-12
    @test cval(adf, "melt_rate", 2001) ≈ 50.0 / 108.0 atol=1e-12
    @test cval(adf, "melt_com_dowy", 2001) ≈ 105.0 atol=1e-12  # cum 25 of 50 total

    # Snowpack persists through Sep 30 -> snow_off + melt-season metrics censored
    dfp = snow_gage([2001]; swe_fn=persist_swe)
    collp = AnnualCollector()
    calculate_snow_metrics(dfp; collector=collp)
    adfp = collector_frame(collp)
    @test cval(adfp, "snow_on_dowy", 2001) ≈ 205.0 atol=1e-12  # 2*(t-200) >= 10
    @test isnan(cval(adfp, "snow_off_dowy", 2001))
    @test isnan(cval(adfp, "melt_season_days", 2001))
    @test isnan(cval(adfp, "melt_rate", 2001))
    @test cval(adfp, "swe_max_dowy", 2001) ≈ 365.0 atol=1e-12  # peak on last day
    @test cval(adfp, "melt_before_peak", 2001) == 0.0
    @test isnan(cval(adfp, "melt_before_peak_pct", 2001))      # zero total melt
    @test isnan(cval(adfp, "melt_com_dowy", 2001))

    # Anchor spell spans the ENTIRE water year -> both ends censored
    dff = snow_gage([2001]; swe_fn=fullyear_swe)
    collf = AnnualCollector()
    calculate_snow_metrics(dff; collector=collf)
    adff = collector_frame(collf)
    @test isnan(cval(adff, "snow_on_dowy", 2001))
    @test isnan(cval(adff, "snow_off_dowy", 2001))
    @test isnan(cval(adff, "melt_season_days", 2001))
    @test isnan(cval(adff, "melt_rate", 2001))
    @test cval(adff, "snow_cover_days", 2001) ≈ 365.0 atol=1e-12
    @test cval(adff, "ssm", 2001) ≈ 1.0 atol=1e-12
    @test cval(adff, "swe_max", 2001) ≈ 50.0 atol=1e-12
end

@testset "Peak on the spell's last day: 1-day melt season (plan §7.16)" begin
    # Monotone accumulation, then instant melt-out below threshold
    lastday_swe(t) = (100 <= t <= 149) ? 2.0 * (t - 99) : 0.0   # peak 100 at t=149
    df = snow_gage([2001]; swe_fn=lastday_swe)
    coll = AnnualCollector()
    calculate_snow_metrics(df; collector=coll)
    adf = collector_frame(coll)
    @test cval(adf, "swe_max_dowy", 2001) ≈ 149.0 atol=1e-12
    @test cval(adf, "snow_off_dowy", 2001) ≈ 150.0 atol=1e-12
    @test cval(adf, "melt_season_days", 2001) ≈ 1.0 atol=1e-12
    @test cval(adf, "melt_rate", 2001) ≈ 100.0 atol=1e-12       # swe_max / 1 day
end

@testset "swe_max_to_ppt gating (plan §7.9, §7.14)" begin
    # Year not in valid_climate_years -> NaN despite valid SWE
    df = snow_gage([1999, 2001, 2003]; swe_fn=triangle_swe)
    coll = AnnualCollector()
    calculate_snow_metrics(df; valid_climate_years=[1999, 2003], collector=coll)
    adf = collector_frame(coll)
    @test cval(adf, "swe_max_to_ppt", 1999) ≈ 242.0 / 730.0 atol=1e-12
    @test isnan(cval(adf, "swe_max_to_ppt", 2001))
    @test cval(adf, "swe_max_to_ppt", 2003) ≈ 242.0 / 730.0 atol=1e-12
    # Other metrics unaffected by the climate gate
    @test cval(adf, "swe_max", 2001) ≈ 242.0 atol=1e-12

    # Annual PPT at/below the 10 mm floor -> NaN (strict > floor)
    dflow = snow_gage([2001]; swe_fn=triangle_swe, ppt=0.02)    # total 7.3 mm
    coll2 = AnnualCollector()
    calculate_snow_metrics(dflow; valid_climate_years=[2001], collector=coll2)
    @test isnan(cval(collector_frame(coll2), "swe_max_to_ppt", 2001))

    # LEGACY fallback (valid_climate_years=nothing): runoff-ratio convention —
    # sum the year's non-NaN PPT days (5 NaN days -> total 720), floor still applies
    ppt_gappy(t) = (10 <= t <= 14) ? NaN : 2.0
    dfleg = snow_gage([2001]; swe_fn=triangle_swe, ppt=ppt_gappy)
    coll3 = AnnualCollector()
    calculate_snow_metrics(dfleg; valid_climate_years=nothing, collector=coll3)
    @test cval(collector_frame(coll3), "swe_max_to_ppt", 2001) ≈ 242.0 / 720.0 atol=1e-12

    # No PPT column at all -> metric NaN, everything else computes
    dfnop = select(snow_gage([2001]; swe_fn=triangle_swe), Not(:PPT))
    coll4 = AnnualCollector()
    calculate_snow_metrics(dfnop; collector=coll4)
    adf4 = collector_frame(coll4)
    @test isnan(cval(adf4, "swe_max_to_ppt", 2001))
    @test cval(adf4, "swe_max", 2001) ≈ 242.0 atol=1e-12
end

@testset "Preprocessor: valid_swe_years + SWE interpolation (plan §7.10)" begin
    function raw_year(yr; swe_fn, mutate=nothing)
        dates = collect(Date(yr - 1, 10, 1):Day(1):Date(yr, 9, 30))
        n = length(dates)
        d = DataFrame(date=dates, Q=fill(1.0, n), water_year=fill(yr, n),
                      month=month.(dates), dowy=collect(1:n),
                      SWE=Vector{Union{Missing, Float64}}(Float64[swe_fn(t) for t in 1:n]))
        mutate !== nothing && mutate(d)
        return d
    end

    d1 = raw_year(2001; swe_fn=t -> 5.0)                              # clean
    d2 = raw_year(2002; swe_fn=t -> 5.0,
                  mutate = d -> (d.SWE[50:80] .= missing))            # 31 missing > 30
    d3 = raw_year(2003; swe_fn=t -> Float64(t),
                  mutate = d -> (d.SWE[150:152] .= missing))          # 3-day internal gap
    d4 = raw_year(2004; swe_fn=t -> 5.0,
                  mutate = d -> (d.SWE[200] = -1.0))                  # negative SWE
    pp = preprocess_daily_data(vcat(d1, d2, d3, d4))

    @test pp.valid_years == [2001, 2002, 2003, 2004]      # Q is clean everywhere
    @test pp.valid_swe_years == [2001, 2003]              # 2002 too many NAs, 2004 negative
    @test "SWE" in names(pp.data)
    # 3-day gap linearly interpolated between SWE(149)=149 and SWE(153)=153
    y3 = pp.data[pp.data.water_year .== 2003, :]
    @test collect(skipmissing(y3.SWE[150:152])) ≈ [150.0, 151.0, 152.0] atol=1e-12
    # valid_swe_years ⊆ valid_years by construction
    @test issubset(Set(pp.valid_swe_years), Set(pp.valid_years))

    # No SWE column -> field present and empty
    pp2 = preprocess_daily_data(select(d1, Not(:SWE)))
    @test pp2.valid_swe_years == Int[]
    @test !("SWE" in names(pp2.data))
end

@testset "No-fallback gate + schema contract (plan §7.12)" begin
    df = snow_gage([1999, 2001, 2003]; swe_fn=triangle_swe)

    # SWE column in gage_data but snow_data NOT passed -> ZERO snow keys
    sigs = calculate_all_signatures(df, false; gage_id="snowtest")
    @test isempty(snow_keys_of(sigs))

    # snow_data=nothing explicitly -> still zero snow keys
    sigs2 = calculate_all_signatures(df, false; gage_id="snowtest", snow_data=nothing)
    @test isempty(snow_keys_of(sigs2))

    # 0-row snow_data (SWE column present) -> the FULL NaN key set
    empty_snow = df[df.water_year .== -1, :]
    sigs3 = calculate_all_signatures(df, false; gage_id="snowtest", snow_data=empty_snow)
    ks = snow_keys_of(sigs3)
    @test length(ks) == 14 * 8
    @test all(isnan(sigs3[k]) for k in ks)

    # With changepoint config: 16 keys per base (8 stats + 8 Pettitt), all present
    res_cp = calculate_snow_metrics(empty_snow; changepoint=SNOW_CP)
    @test length(res_cp) == 14 * 16
    for m in StreamflowSignatures.SNOW_METRICS
        @test haskey(res_cp, "$(m)_mean")
        @test haskey(res_cp, "$(m)_pettitt_pval")
    end

    # Defensive direct-API misuse: no SWE column -> NaN stat set, no throw
    res_noswe = calculate_snow_metrics(select(df, Not(:SWE)))
    @test length(res_noswe) == 14 * 8
    @test all(isnan, values(res_noswe))

    # ...and the SAME contract with changepoint enabled: full 224-key set
    # (impl-review finding 1: the old empty_stats path dropped the Pettitt keys)
    res_noswe_cp = calculate_snow_metrics(select(df, Not(:SWE)); changepoint=SNOW_CP)
    @test length(res_noswe_cp) == 14 * 16
    @test all(isnan, values(res_noswe_cp))
    for m in StreamflowSignatures.SNOW_METRICS
        @test haskey(res_noswe_cp, "$(m)_pettitt_pval")
    end
end

@testset "Malformed year grid rejected (impl-review finding 2)" begin
    # Duplicate one interior day and drop another: endpoints still read 1..365 and
    # the row count is right, but the sequence is broken -> year must be all-NaN
    df = snow_gage([2001]; swe_fn=triangle_swe)
    dup = df[df.dowy .== 150, :]
    df_bad = sort(vcat(df[df.dowy .!= 151, :], dup), :dowy)
    @test nrow(df_bad) == 365
    @test df_bad.dowy[1] == 1 && df_bad.dowy[end] == 365   # endpoint checks would pass

    coll = AnnualCollector()
    calculate_snow_metrics(df_bad; collector=coll)
    adf = collector_frame(coll)
    for m in StreamflowSignatures.SNOW_METRICS
        @test isnan(cval(adf, m, 2001))
    end

    # Control: the intact frame computes
    coll_ok = AnnualCollector()
    calculate_snow_metrics(df; collector=coll_ok)
    @test cval(collector_frame(coll_ok), "swe_max", 2001) ≈ 242.0 atol=1e-12
end

@testset "SWE-invalid year cannot leak into snow metrics (plan §7.13)" begin
    # Gage A: three clean triangle years + one COMPLETE but SWE-invalid year
    # (deep constant 500 mm with one negative day -> preprocessor rejects it for
    # snow; its values are NOT NaN, so any leak would visibly shift swe_max)
    clean_years = [1999, 2001, 2003]
    dfs = [snow_gage([yr]; swe_fn=triangle_swe) for yr in clean_years]
    bad = snow_gage([2002]; swe_fn=t -> t == 200 ? -5.0 : 500.0)
    df_a = sort(vcat(dfs..., bad), [:water_year, :dowy])
    pp_a = preprocess_daily_data(df_a)
    @test pp_a.valid_swe_years == clean_years
    @test 2002 in pp_a.valid_years                     # Q-valid, so it stays in data

    snow_a = pp_a.data[in.(pp_a.data.water_year, Ref(Set(pp_a.valid_swe_years))), :]
    stats_a = calculate_snow_metrics(snow_a)

    # Gage B: the corrupt year hard-removed before preprocessing
    pp_b = preprocess_daily_data(sort(vcat(dfs...), [:water_year, :dowy]))
    snow_b = pp_b.data[in.(pp_b.data.water_year, Ref(Set(pp_b.valid_swe_years))), :]
    stats_b = calculate_snow_metrics(snow_b)

    @test isequal(stats_a, stats_b)

    # Leak detector: computing on the UNFILTERED frame (the bug the orchestrator
    # gate prevents) must give a different swe_max_mean (2002's 500 mm would count)
    stats_leak = calculate_snow_metrics(pp_a.data)
    @test stats_leak["swe_max_mean"] != stats_a["swe_max_mean"]
end

@testset "Collector + harness invariants (plan §7.11)" begin
    years = collect(1991:2011)   # 21 years, 5 leap WYs
    df = snow_gage(years; swe_fn=triangle_swe)

    # With/without collector: statistics identical
    res_plain = calculate_snow_metrics(df; valid_climate_years=years,
                                       trend_completeness=0.8, decade_completeness=0.8,
                                       changepoint=SNOW_CP)
    coll = AnnualCollector()
    res_coll = calculate_snow_metrics(df; valid_climate_years=years,
                                      trend_completeness=0.8, decade_completeness=0.8,
                                      changepoint=SNOW_CP, collector=coll)
    @test isequal(res_plain, res_coll)
    @test length(res_plain) == 14 * 16

    # Collector: all 14 signatures, dense (one row per year each)
    adf = collector_frame(coll)
    @test Set(adf.signature) == Set(StreamflowSignatures.SNOW_METRICS)
    for m in StreamflowSignatures.SNOW_METRICS
        @test sum(adf.signature .== m) == length(years)
    end
    # NaN placeholders are collected too (identical snowpacks: melt_before is 0, valid)
    @test all(adf[adf.signature .== "melt_before_peak", :value] .== 0.0)

    # Identical years -> exact summary stats; constant series -> zero Theil-Sen slope
    @test res_plain["swe_max_mean"] ≈ 242.0 atol=1e-12
    @test res_plain["swe_max_senn_slp"] == 0.0
    @test res_plain["ssm_mean"] ≈ 1.0 atol=1e-12
    @test res_plain["swe_apr1_median"] ≈ 149.0 atol=1e-12   # 16 non-leap vs 5 leap years

    # Zero-warnings through the full orchestrator with snow wired in
    pp = preprocess_daily_data(df)
    snow_data = pp.data[in.(pp.data.water_year, Ref(Set(pp.valid_swe_years))), :]
    @test_logs min_level=Logging.Warn calculate_all_signatures(
        pp.data, true;
        gage_id="snowtest", seasonal_flags=pp.seasonal_flags,
        trend_completeness=0.8, decade_completeness=0.8,
        climate_data=pp.data[in.(pp.data.water_year, Ref(Set(pp.valid_climate_years))), :],
        snow_data=snow_data, snow_climate_years=pp.valid_climate_years,
        changepoint=SNOW_CP)
end

@testset "QA flag semantics: present-but-missing snow columns count (plan §7.17)" begin
    base = DataFrame(gage_id=["g1"],
                     Qann_mean=[1.0], BFI_Eckhardt_mean=[0.5], D50_day_mean=[150.0],
                     TQmean_mean=[30.0], flashinessRB_mean=[0.2], avg_storage_mean=[10.0])
    snow_missing = DataFrame(swe_max_mean=[missing], ssm_mean=[missing],
                             melt_rate_mean=[missing], swe_apr1_mean=[missing])

    # 4 of 10 numeric columns missing -> 0.4 > 0.30 threshold -> flagged
    flagged = compute_qa_flags(hcat(base, snow_missing))
    @test flagged.flagged_for_high_na[1] == true

    # Same row WITHOUT the snow columns -> 0/6 -> not flagged. Pins that ABSENT
    # columns never enter the denominator (the golden-comparison whitelist for
    # flagged_for_high_na rests on exactly this present-vs-absent distinction).
    unflagged = compute_qa_flags(base)
    @test unflagged.flagged_for_high_na[1] == false
end

end  # @testset "Snow metrics"

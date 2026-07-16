# Tests for the stats floor (min_values_for_stats) — July 2026.
# A metric with fewer non-NaN annual values than the floor emits NaN for ALL 8
# statistics AND its changepoint fields. The annual-values collector still receives
# the full series (collection precedes gating). Recession & elasticity are exempt
# at the orchestration layer (the floor kwarg is never passed to them), mirroring
# the trend_completeness exemption.
# Motivating case: gage 07292500 snow_on_dowy — 4 clustered non-NA years (1982-85)
# passed the own-range trend gate and produced a 4-point Theil-Sen slope.

using Test
using DataFrames
using Dates
using Statistics
using Logging
using StreamflowSignatures

# Permissive changepoint config so cp WOULD compute absent the floor
const FLOOR_CP = (start_year=1900, end_year=2100, min_total_obs=3, min_segment_obs=1)

"""Daily gage frame; `swe_of(yr, t)` maps (water_year, dowy) -> SWE."""
function floor_gage(years; swe_of=nothing, ppt=2.0)
    dfs = DataFrame[]
    for yr in years
        dates = collect(Date(yr - 1, 10, 1):Day(1):Date(yr, 9, 30))
        n = length(dates)
        d = DataFrame(gage_id="floortest", date=dates,
                      Q=[1.0 + 0.3 * sin(2pi * t / 365) for t in 1:n],
                      water_year=fill(yr, n), month=month.(dates),
                      dowy=collect(1:n), PPT=fill(Float64(ppt), n))
        if swe_of !== nothing
            d[!, :SWE] = Float64[swe_of(yr, t) for t in 1:n]
        end
        push!(dfs, d)
    end
    return vcat(dfs...)
end

@testset "Stats floor (min_values_for_stats)" begin

@testset "generate_stats unit behavior" begin
    df19 = DataFrame(water_year=collect(2000:2018), m=collect(1.0:19.0))
    r19 = generate_stats(df19; value_cols=["m"], min_values_for_stats=20,
                         changepoint=FLOOR_CP)
    @test length(r19) == 16                      # 8 stats + 8 cp fields, all present
    @test all(isnan, values(r19))                # incl. mean/median AND cp fields

    df20 = DataFrame(water_year=collect(2000:2019), m=collect(1.0:20.0))
    r20 = generate_stats(df20; value_cols=["m"], min_values_for_stats=20)
    @test r20["m_mean"] ≈ 10.5
    @test r20["m_senn_slp"] ≈ 1.0

    # The floor counts NON-NaN values, not rows: 25 rows, 19 non-NaN -> gated
    v = Vector{Float64}(1.0:25.0); v[1:6] .= NaN
    df_nan = DataFrame(water_year=collect(2000:2024), m=v)
    @test isnan(generate_stats(df_nan; value_cols=["m"],
                               min_values_for_stats=20)["m_mean"])
    v2 = Vector{Float64}(1.0:25.0); v2[1:5] .= NaN   # 20 non-NaN -> computed
    df_ok = DataFrame(water_year=collect(2000:2024), m=v2)
    @test !isnan(generate_stats(df_ok; value_cols=["m"],
                                min_values_for_stats=20)["m_mean"])

    # Legacy behavior with the floor absent: 4 values still compute
    df4 = DataFrame(water_year=collect(2000:2003), m=[1.0, 2.0, 3.0, 4.0])
    @test generate_stats(df4; value_cols=["m"])["m_mean"] ≈ 2.5
end

@testset "Collector receives gated series (collection precedes the floor)" begin
    df19 = DataFrame(water_year=collect(2000:2018), m=collect(1.0:19.0))
    c_floor = AnnualCollector()
    generate_stats(df19; value_cols=["m"], min_values_for_stats=20, collector=c_floor)
    c_plain = AnnualCollector()
    generate_stats(df19; value_cols=["m"], collector=c_plain)
    @test isequal(c_floor.value, c_plain.value)
    @test length(c_floor.value) == 19
end

@testset "Orchestrator: 07292500-shaped snow regression" begin
    years = collect(1994:2015)                       # 22-year record
    snowy = Set([1996, 1997, 1998, 1999])            # 4 clustered snow years
    swe_of = (yr, t) -> (yr in snowy && 100 <= t <= 160) ? 50.0 : 0.0
    df = floor_gage(years; swe_of=swe_of)

    sig_floor = calculate_all_signatures(df, false; gage_id="floortest",
        snow_data=df, min_values_for_stats=20)
    sig_plain = calculate_all_signatures(df, false; gage_id="floortest",
        snow_data=df)

    # Without the floor: the 4 clustered years produce a Sen slope (the leak)
    @test !isnan(sig_plain["snow_on_dowy_senn_slp"])
    @test !isnan(sig_plain["snow_on_dowy_mean"])
    # With the floor: ALL 8 stats NaN for the sparse timing metric
    for suf in ["_mean", "_median", "_senn_slp", "_linear_slp",
                "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval"]
        @test isnan(sig_floor["snow_on_dowy$(suf)"])
        @test isnan(sig_floor["swe_max_dowy$(suf)"])
    end
    # Dense magnitude metric (22 values incl. zeros) is untouched by the floor
    @test !isnan(sig_floor["swe_max_mean"])
    @test !isnan(sig_floor["swe_max_senn_slp"])
    @test sig_floor["swe_max_mean"] ≈ sig_plain["swe_max_mean"]
end

@testset "Recession & elasticity exempt (floor never reaches them)" begin
    # Recession events only in the first 12 of 25 years -> log_a series has < 20
    # values; the exemption must keep its stats computed even with the floor set.
    years = collect(1991:2015)
    event_years = Set(years[1:12])
    function q_of(yr, t)
        if yr in event_years && 100 <= t <= 130
            return 8.0 * 0.9^(t - 100) + 1.0     # monthly-spike decay
        end
        return 1.0 + 0.3 * sin(2pi * t / 365)
    end
    dfs = DataFrame[]
    for yr in years
        dates = collect(Date(yr - 1, 10, 1):Day(1):Date(yr, 9, 30))
        n = length(dates)
        push!(dfs, DataFrame(gage_id="rectest", date=dates,
                             Q=Float64[q_of(yr, t) for t in 1:n],
                             water_year=fill(yr, n), month=month.(dates),
                             dowy=collect(1:n), PPT=fill(2.0, n)))
    end
    dfr = vcat(dfs...)
    s_floor = calculate_all_signatures(dfr, true; gage_id="rectest",
        min_values_for_stats=20)
    s_plain = calculate_all_signatures(dfr, true; gage_id="rectest")
    for k in keys(s_plain)
        if startswith(k, "log_a_") || startswith(k, "b_point") ||
           startswith(k, "b_event") || startswith(k, "concavity") ||
           startswith(k, "elasticity_")
            @test isequal(s_floor[k], s_plain[k])
        end
    end
    # Elasticity rolling: 25 climate years -> 15 rolling values < 20, still computed
    @test !isnan(s_plain["elasticity_rolling_mean"])
    @test !isnan(s_floor["elasticity_rolling_mean"])
end

@testset "Zero warnings with the floor active" begin
    years = collect(1994:2015)
    swe_of = (yr, t) -> (yr == 1996 && 100 <= t <= 160) ? 50.0 : 0.0
    df = floor_gage(years; swe_of=swe_of)
    @test_logs min_level=Logging.Warn calculate_all_signatures(
        df, true; gage_id="floortest", snow_data=df,
        min_values_for_stats=20, changepoint=FLOOR_CP)
end

end  # testset

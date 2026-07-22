# Tests for the snow record-anchored decade gate
# (plan: docs/plans/2026-07-22-snow-record-anchored-decade-gate.md).
#
# Part 1: generate_stats force_skip_trends mechanics (external trend suppression
#         through the same path as trend_completeness).
# Part 2: calculate_snow_metrics record-anchored gate — anchored to the SWE-valid
#         record (swe_max non-NaN rows), threshold LINKED to decade_completeness,
#         10 threshold-dependent metrics gated, magnitude metrics exempt, failure
#         NaNs the 6 trend stats only (mean/median + Pettitt survive).

using Test
using DataFrames
using Dates
using Statistics
using StreamflowSignatures

const RDG_CP = (start_year=1980, end_year=2024, min_total_obs=20, min_segment_obs=10)

"Daily SWE gage: snowy years get a triangle profile (peak 242 mm), others all-zero."
function rdg_gage(years; snowy::Function, ppt=2.0)
    tri(t) = 32 <= t <= 152 ? 2.0 * (t - 31) : (t > 152 ? max(0.0, 242.0 - 3.0 * (t - 152)) : 0.0)
    dfs = DataFrame[]
    for yr in years
        dates = collect(Date(yr - 1, 10, 1):Day(1):Date(yr, 9, 30))
        n = length(dates)
        swe = snowy(yr) ? Float64[tri(t) for t in 1:n] : zeros(Float64, n)
        push!(dfs, DataFrame(gage_id="rdgtest", date=dates, Q=fill(1.0, n),
                             water_year=fill(yr, n), month=month.(dates),
                             dowy=collect(1:n), PPT=fill(Float64(ppt), n), SWE=swe))
    end
    return vcat(dfs...)
end

# The 10 threshold-dependent metrics gated by the record-anchored decade gate
const RDG_GATED = ["swe_max_dowy", "snow_on_dowy", "snow_off_dowy", "melt_season_days",
                   "melt_rate", "ssm", "melt_before_peak", "melt_before_peak_pct",
                   "melt_before_peak_to_max_swe", "melt_com_dowy"]

const RDG_TREND_SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho",
                            "_spearman_pval", "_mk_rho", "_mk_pval"]

@testset "Snow record-anchored decade gate" begin

    @testset "generate_stats force_skip_trends mechanics" begin
        df = DataFrame(water_year=collect(2000:2024),
                       a=collect(1.0:25.0), b=collect(2.0:2.0:50.0))

        # Baseline: kwarg omitted vs explicitly nothing — identical results
        base = generate_stats(df; value_cols=["a", "b"], changepoint=RDG_CP)
        same = generate_stats(df; value_cols=["a", "b"], changepoint=RDG_CP,
                              force_skip_trends=nothing)
        @test keys(base) == keys(same)
        @test all(isequal(base[k], same[k]) for k in keys(base))

        # Skip "a" only: 6 trend stats NaN; mean/median computed; Pettitt still runs
        skipped = generate_stats(df; value_cols=["a", "b"], changepoint=RDG_CP,
                                 force_skip_trends=Set(["a"]))
        for sfx in RDG_TREND_SUFFIXES
            @test isnan(skipped["a$(sfx)"])
            @test !isnan(skipped["b$(sfx)"])   # untouched column unaffected
        end
        @test skipped["a_mean"] ≈ mean(1.0:25.0)
        @test skipped["a_median"] ≈ median(1.0:25.0)
        @test !isnan(skipped["a_pettitt_pval"])            # changepoint independent
        @test all(isequal(skipped["b$(s)"], base["b$(s)"])
                  for s in vcat(RDG_TREND_SUFFIXES, ["_mean", "_median", "_pettitt_pval"]))

        # Forced skip works WITHOUT trend_completeness passed (orthogonal paths)
        solo = generate_stats(df; value_cols=["a"], force_skip_trends=Set(["a"]))
        @test isnan(solo["a_senn_slp"])
        @test !isnan(solo["a_mean"])
    end

    @testset "gate fires when snow vanishes from the last decade" begin
        # 30 SWE-valid years 1981-2010; snowy 1981-2000 (20 values -> passes the
        # stats floor), snow-free 2001-2010. The metric's OWN span is 1981-2000
        # (100% complete, own decades complete) so the pre-existing gate and floor
        # both pass -> without the new gate every trend WOULD be computed
        # (control testset below proves it).
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2000)
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                   min_values_for_stats=20, changepoint=RDG_CP)
        for m in RDG_GATED
            for sfx in RDG_TREND_SUFFIXES
                @test isnan(r["$(m)$(sfx)"])
            end
        end
        @test !isnan(r["snow_on_dowy_mean"])          # mean survives (over snowy years)
        # Onset = first day at/above the 10 mm threshold: tri(t) = 2(t-31), so days
        # 32-35 (2,4,6,8 mm) are sub-threshold and dowy 36 (10 mm) starts the spell.
        @test r["snow_on_dowy_mean"] ≈ 36.0
        @test !isnan(r["snow_on_dowy_pettitt_pval"])  # Pettitt survives (20 obs)
        # Magnitude metrics are NOT gated (dense zero-including series legitimately
        # carry the snow-decline signal)
        for m in ["swe_max", "snow_cover_days", "swe_apr1"]
            @test !isnan(r["$(m)_senn_slp"])
        end
    end

    @testset "control: gate off when decade_completeness not passed" begin
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2000)
        r = calculate_snow_metrics(df; trend_completeness=0.6,
                                   min_values_for_stats=20, changepoint=RDG_CP)
        @test !isnan(r["snow_on_dowy_senn_slp"])  # pre-gate behavior: trend computed
    end

    @testset "gate fires on first-decade absence" begin
        # Snowy only 1991-2010: own span 1991-2010 is complete (pre-existing checks
        # pass) but the record's first decade [1981,1990] has zero computable years.
        df = rdg_gage(1981:2010; snowy = yr -> yr >= 1991)
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                   min_values_for_stats=20)
        @test isnan(r["snow_on_dowy_senn_slp"])
    end

    @testset "denominator is SWE-valid years, not calendar years" begin
        # Record 1981-2010 with 2004-2006 SWE-invalid (absent); snowy every valid
        # year EXCEPT 2010. Run at dc = 0.7 where the three checks separate:
        #   new gate, valid-year denominator: last window [2001,2010] den=7
        #     (2001-03, 2007-10), num=6 (2010 snow-free) -> 6/7 = 0.857 >= 0.7 PASSES
        #   new gate, calendar denominator (wrong design): 6/10 = 0.6 < 0.7 would FIRE
        #   own-span decade check (calendar-denominator, pre-existing): own span
        #     1981-2009, last own decade [2000,2009] has 7 of 10 -> 0.7 >= 0.7 PASSES
        # So a computed trend pins the valid-year denominator.
        yrs = vcat(collect(1981:2003), collect(2007:2010))
        df = rdg_gage(yrs; snowy = yr -> yr != 2010)
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.7,
                                   min_values_for_stats=20)
        @test !isnan(r["snow_on_dowy_senn_slp"])
    end

    @testset "linked threshold: gate obeys the passed decade_completeness" begin
        # Snowy through 2007, snow-free 2008-2010: own span 1981-2007 complete
        # (own checks pass at any dc). New-gate last window [2001,2010]: 7/10 = 0.7.
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2007)
        r60 = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.6,
                                     min_values_for_stats=20)
        r80 = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                     min_values_for_stats=20)
        @test !isnan(r60["snow_on_dowy_senn_slp"])
        @test isnan(r80["snow_on_dowy_senn_slp"])
    end

    @testset "short record (<10-year span): gate skipped" begin
        df = rdg_gage(2001:2008; snowy = yr -> yr <= 2004)  # record span 8
        r = calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8)
        # 4 snowy years >= min_rows=3; own decade check skipped (span < 10);
        # the new gate must NOT fire on a <10-year record either
        @test !isnan(r["snow_on_dowy_senn_slp"])
    end

    @testset "collector unaffected by the gate" begin
        df = rdg_gage(1981:2010; snowy = yr -> yr <= 2000)
        c_on = AnnualCollector()
        c_off = AnnualCollector()
        calculate_snow_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                               min_values_for_stats=20, collector=c_on)
        calculate_snow_metrics(df; trend_completeness=0.6,
                               min_values_for_stats=20, collector=c_off)
        @test c_on.signature == c_off.signature
        @test c_on.water_year == c_off.water_year
        @test all(isequal.(c_on.value, c_off.value))
    end

    @testset "config flag is shipped true" begin
        @test StreamflowSignatures.CFG_SNOW_RECORD_DECADE_GATE === true
    end

end

# Regression test for the seasonal runoff-ratio completeness masking (fixed 2026-08-24).
# The bug (CHANGELOG Known Issues): analyze_Q_PPT_relationships looked up FULL-name
# flags (:winter_complete, …) while the preprocessor emits ABBREVIATED names
# (win_complete/spr_complete/sum_complete/fal_complete — io.jl seasonal_rows; the
# identical, correct lookup lives in flow_volumes.jl). The existence check failed
# silently, so the masking block was dead code and incomplete seasons were never
# NaN'd. The bug was silent precisely because no test asserted the mask FIRES with
# preprocessor-shaped flags — this file does.

using Test
using DataFrames
using StreamflowSignatures

@testset "Seasonal runoff-ratio completeness masking (regression 2026-08-24)" begin
    # Synthetic gage: 10 water years x 12 months x 28 days. Q varies by year so
    # per-year ratios differ (a constant series would hide the mask in the mean).
    wy_v = Int[]; mo_v = Int[]; q_v = Float64[]; ppt_v = Float64[]
    for wy in 2001:2010
        for month in [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            for _ in 1:28
                push!(wy_v, wy); push!(mo_v, month)
                push!(q_v, 1.0 + 0.1 * (wy - 2001)); push!(ppt_v, 2.0)
            end
        end
    end
    df = DataFrame(water_year = wy_v, month = mo_v, Q = q_v, PPT = ppt_v)

    # Preprocessor-shaped flags (abbreviated names, Bool):
    # summer of WY2005 and winter of WY2007 incomplete, everything else complete.
    flags = DataFrame(
        water_year   = collect(2001:2010),
        win_complete = [wy != 2007 for wy in 2001:2010],
        spr_complete = trues(10),
        sum_complete = [wy != 2005 for wy in 2001:2010],
        fal_complete = trues(10),
    )

    c0 = AnnualCollector()
    r0 = analyze_Q_PPT_relationships(df; collector = c0)
    c1 = AnnualCollector()
    r1 = analyze_Q_PPT_relationships(df; seasonal_flags = flags, collector = c1)

    function getval(c, sig, wy)
        idx = findfirst(i -> c.signature[i] == sig && c.water_year[i] == wy,
                        eachindex(c.signature))
        @test idx !== nothing
        idx === nothing ? NaN : c.value[idx]
    end

    # The mask fires exactly where flagged (two different seasons pin the mapping)…
    @test isnan(getval(c1, "summer_runoff_ratio", 2005))
    @test isnan(getval(c1, "winter_runoff_ratio", 2007))
    # …and nowhere else.
    @test !isnan(getval(c1, "summer_runoff_ratio", 2004))
    @test !isnan(getval(c1, "winter_runoff_ratio", 2006))
    @test !isnan(getval(c0, "summer_runoff_ratio", 2005))
    @test !isnan(getval(c0, "winter_runoff_ratio", 2007))

    # Annual ratio and unflagged seasons are untouched by the flags.
    @test getval(c1, "annual_runoff_ratio", 2005) == getval(c0, "annual_runoff_ratio", 2005)
    @test getval(c1, "annual_runoff_ratio", 2007) == getval(c0, "annual_runoff_ratio", 2007)
    @test getval(c1, "spring_runoff_ratio", 2005) == getval(c0, "spring_runoff_ratio", 2005)
    @test getval(c1, "fall_runoff_ratio", 2007) == getval(c0, "fall_runoff_ratio", 2007)

    # Summary statistics reflect the mask (year-varying series ⇒ the mean moves for
    # masked families only).
    @test r1["summer_runoff_ratio_mean"] != r0["summer_runoff_ratio_mean"]
    @test r1["winter_runoff_ratio_mean"] != r0["winter_runoff_ratio_mean"]
    @test r1["fall_runoff_ratio_mean"] == r0["fall_runoff_ratio_mean"]
    @test r1["spring_runoff_ratio_mean"] == r0["spring_runoff_ratio_mean"]
    @test r1["annual_runoff_ratio_mean"] == r0["annual_runoff_ratio_mean"]

    # Guard behavior unchanged: a flags frame missing seasonal columns is ignored
    # (no crash, no masking).
    partial_flags = DataFrame(water_year = collect(2001:2010))
    c2 = AnnualCollector()
    r2 = analyze_Q_PPT_relationships(df; seasonal_flags = partial_flags, collector = c2)
    @test !isnan(getval(c2, "summer_runoff_ratio", 2005))
    @test r2["summer_runoff_ratio_mean"] == r0["summer_runoff_ratio_mean"]
end

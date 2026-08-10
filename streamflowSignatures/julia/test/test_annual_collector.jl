# Tests for the opt-in AnnualCollector (per-year annual values export).
# See docs/plans/annual_values_export_plan.md §4.1.
#
# Key invariants:
#   1. The collector has ZERO effect on computed statistics (frozen contract).
#   2. Collected rows are exactly the series generate_stats() consumed.
#   3. No warnings are emitted with the collector enabled (a botched kwarg
#      thread would be silently swallowed by calculate_all_signatures' try/catch).

using StreamflowSignatures
using Test
using DataFrames
using Dates
using Statistics
using Logging
using Random

# Signatures guaranteed to produce dense annual series for a healthy 25-year
# synthetic gage with climate (caller-pruned/sparse signatures like recession
# metrics are excluded — coverage for those depends on the data).
const EXPECTED_DENSE_SIGNATURES = [
    # Flow volumes (21)
    "Qann", "Qwin", "Qspr", "Qsum", "Qfal",
    "Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40", "Q50",
    "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99", "Q95_Q10",
    # FDC (3)
    "FDCall", "FDC90th", "FDCmid",
    # Baseflow fixed-parameter (2)
    "BFI_Eckhardt", "BFI_LyneHollick",
    # Flashiness
    "flashinessRB",
    # Flow timing (15)
    "D1_day", "D5_day", "D10_day", "D20_day", "D30_day", "D40_day", "D50_day",
    "D60_day", "D70_day", "D80_day", "D90_day", "D95_day", "D99_day",
    "D25_to_D75", "Dmax",
    # Pulses (14)
    "n_high_pulses_year", "n_low_pulses_year", "dur_high_pulses_year", "dur_low_pulses_year",
    "n_high_pulses_all", "n_low_pulses_all", "dur_high_pulses_all", "dur_low_pulses_all",
    "TQmean", "Flow_Reversals_annual", "Flow_Reversals_winter", "Flow_Reversals_spring",
    "Flow_Reversals_summer", "Flow_Reversals_fall",
    # Negative days
    "negative_ann",
    # Recession event count — dense even with ZERO recession events (an annual
    # count of 0 is a valid non-NaN value, so events_df keeps every year)
    "n_recession_events",
    # Runoff ratios (5)
    "annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
    "summer_runoff_ratio", "fall_runoff_ratio",
    # Elasticity series (2)
    "elasticity_rolling", "elasticity_annual",
    # Q-P seasonality (2)
    "qp_slope_sd", "qp_bimodality",
    # Storage
    "avg_storage",
    # Streamflow drought (10): dense — a drought-free year is a valid 0, and the
    # 25-year synthetic gage clears the 10-year threshold floor
    "drought_duration_fixed_p2", "drought_duration_fixed_p5", "drought_duration_fixed_p10",
    "drought_duration_fixed_p20", "drought_duration_fixed_p30",
    "drought_deficit_fixed_p2", "drought_deficit_fixed_p5", "drought_deficit_fixed_p10",
    "drought_deficit_fixed_p20", "drought_deficit_fixed_p30",
]

function collector_synthetic_gage(n_years=25; seed=42)
    Random.seed!(seed)
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years
    dates = [start_date + Day(i-1) for i in 1:n_days]
    doy = [dayofyear(d) for d in dates]
    base_flow = 2.0 .+ 1.5 .* sin.(2π .* (doy .- 100) ./ 365)
    Q = max.(0.01, base_flow .+ rand(n_days) .* 0.5)
    df = DataFrame(gage_id=fill("TESTAC1", n_days), date=dates, Q=Q)
    df = add_water_year_columns(df)
    df.PPT = max.(0.0, df.Q .* 2.5 .+ rand(nrow(df)))
    return df
end

# Deterministic gage with a monthly spike + pure exponential recession
# (Q_{t+1} = 0.95 Q_t). The noisy gage above produces ZERO recession events,
# so this gage is what exercises the recession metrics and the
# recession-parameterized BFI collector paths (alpha scalar = 0.95 exactly).
function collector_recession_gage(n_years=25)
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years
    dates = [start_date + Day(i-1) for i in 1:n_days]
    Q = Vector{Float64}(undef, n_days)
    q = 10.0
    for i in 1:n_days
        if day(dates[i]) == 1
            q = 10.0
        else
            q *= 0.95
        end
        Q[i] = max(q, 0.001)
    end
    df = DataFrame(gage_id=fill("TESTREC1", n_days), date=dates, Q=Q)
    return add_water_year_columns(df)
end

# Assert that for every collected signature whose summary mean is non-NaN,
# the mean/median recomputed from the collected rows reproduce the summary
# exactly, and at least min_rows non-NaN values back it (generate_stats only
# computes a mean from >= 3 values). This is the full completeness contract —
# a collector that dropped or duplicated rows would fail here.
function assert_collector_completeness(c::AnnualCollector, r::Dict{String, Any}; min_rows=3)
    for sig in unique(c.signature)
        mask = c.signature .== sig
        vals = filter(!isnan, c.value[mask])
        smean = get(r, "$(sig)_mean", NaN)
        smedian = get(r, "$(sig)_median", NaN)
        if smean isa Real && !isnan(smean)
            @test length(vals) >= min_rows
            @test isapprox(mean(vals), smean; rtol=1e-9, atol=1e-12)
            @test isapprox(median(vals), smedian; rtol=1e-9, atol=1e-12)
        end
    end
end

@testset "Annual Collector" begin
    df = collector_synthetic_gage()

    @testset "Collector has zero effect on statistics" begin
        r_without = calculate_all_signatures(df, true)
        c = AnnualCollector()
        r_with = calculate_all_signatures(df, true; collector=c)

        # isequal treats NaN == NaN, so this asserts full key AND value identity
        @test isequal(r_without, r_with)
        @test Set(keys(r_without)) == Set(keys(r_with))
        @test length(c.value) > 0
    end

    @testset "No warnings emitted with collector enabled" begin
        # calculate_all_signatures wraps every signature in try/catch that only
        # warns — a broken collector thread would be silently swallowed there.
        c = AnnualCollector()
        @test_logs min_level=Logging.Warn calculate_all_signatures(df, true; collector=c)
    end

    @testset "Signature coverage regression" begin
        c = AnnualCollector()
        r = calculate_all_signatures(df, true; collector=c)
        collected = Set(c.signature)

        # This gage has zero recession events, so the collected set must be
        # EXACTLY the dense list — nothing missing, nothing unexpected.
        @test collected == Set(EXPECTED_DENSE_SIGNATURES)

        # Every collected signature must exist in the summary output
        for sig in collected
            @test haskey(r, "$(sig)_mean")
        end

        # Dense signatures: one row per water year (25); the two elasticity
        # series have their own deterministic lengths (windows / pairs).
        n_years = length(unique(df.water_year))
        for sig in ["Qann", "Q50", "D50_day", "flashinessRB", "BFI_Eckhardt",
                    "annual_runoff_ratio", "qp_slope_sd", "avg_storage", "negative_ann"]
            @test count(==(sig), c.signature) == n_years
        end
        @test count(==("elasticity_rolling"), c.signature) == n_years - 11 + 1
        @test count(==("elasticity_annual"), c.signature) == n_years - 1
    end

    @testset "Collector completeness: recomputed mean/median reproduce summary" begin
        c = AnnualCollector()
        r = calculate_all_signatures(df, true; collector=c)
        assert_collector_completeness(c, r)
    end

    @testset "Recession + parameterized BFI collector paths" begin
        # The main synthetic gage has no recession events; this deterministic
        # exponential-recession gage exercises those collector paths.
        dfr = collector_recession_gage()
        c = AnnualCollector()
        r = calculate_all_signatures(dfr, false; collector=c)

        n_years = length(unique(dfr.water_year))
        for sig in ["BFI_Eckhardt_param", "BFI_LyneHollick_param",
                    "n_recession_events", "log_a_pointcloud", "log_a_events",
                    "b_pointcloud", "b_events", "concavity", "alpha_linear"]
            @test count(==(sig), c.signature) == n_years
        end

        # Q_{t+1} = 0.95 Q_t => the linear-reservoir alpha scalar is exactly 0.95
        @test isapprox(r["recession_alpha_point_cloud_linear_reservoir"], 0.95; atol=1e-9)

        # Full completeness contract on this gage too
        assert_collector_completeness(c, r)

        # Frozen contract + zero warnings on the recession paths as well
        r_plain = calculate_all_signatures(dfr, false)
        @test isequal(r_plain, r)
        c2 = AnnualCollector()
        @test_logs min_level=Logging.Warn calculate_all_signatures(dfr, false; collector=c2)
    end

    @testset "No duplicate (signature, water_year) keys per gage" begin
        c = AnnualCollector()
        calculate_all_signatures(df, true; collector=c)

        seen = Set{Tuple{String, Int32}}()
        has_dup = false
        for i in eachindex(c.signature)
            p = (c.signature[i], c.water_year[i])
            if p in seen
                has_dup = true
                break
            end
            push!(seen, p)
        end
        @test !has_dup
    end

    @testset "Collected Qann rows match independent annual sums" begin
        c = AnnualCollector()
        calculate_flow_vols_by_year(df; collector=c)

        expected = combine(groupby(df, :water_year), :Q => sum => :Qann)
        qann_mask = c.signature .== "Qann"
        @test sum(qann_mask) == nrow(expected)

        for row in eachrow(expected)
            idx = findfirst(qann_mask .& (c.water_year .== Int32(row.water_year)))
            @test idx !== nothing
            @test isapprox(c.value[idx], row.Qann; atol=1e-9)
        end
    end

    @testset "NaN rows collected for seasonal-flag exclusions" begin
        years = sort(unique(df.water_year))
        n = length(years)
        sf = DataFrame(
            water_year = years,
            win_complete = trues(n),
            spr_complete = trues(n),
            sum_complete = trues(n),
            fal_complete = trues(n),
        )
        excluded_year = years[3]
        sf[3, :win_complete] = false

        c = AnnualCollector()
        calculate_flow_vols_by_year(df; seasonal_flags=sf, collector=c)

        idx = findfirst((c.signature .== "Qwin") .& (c.water_year .== Int32(excluded_year)))
        @test idx !== nothing
        @test isnan(c.value[idx])

        # A non-excluded year's Qwin must be non-NaN
        idx_ok = findfirst((c.signature .== "Qwin") .& (c.water_year .== Int32(years[5])))
        @test idx_ok !== nothing
        @test !isnan(c.value[idx_ok])
    end

    @testset "Collection happens before the min_rows early return" begin
        c = AnnualCollector()
        small = DataFrame(water_year=[2000, 2001], m=[1.5, 2.5])
        s = generate_stats(small; value_cols=["m"], collector=c)

        @test isnan(s["m_mean"])          # stats gated (nrow < min_rows)
        @test length(c.value) == 2        # but values were still collected
        @test c.water_year == Int32[2000, 2001]
        @test c.value == [1.5, 2.5]
    end

    @testset "Year-skip guard: NaN and missing years" begin
        # NaN year (via generate_stats early-return path — no sort/convert issues)
        c = AnnualCollector()
        dfnan = DataFrame(water_year=[2000.0, NaN], m=[1.0, 2.0])
        generate_stats(dfnan; value_cols=["m"], collector=c)
        @test c.water_year == Int32[2000]
        @test c.value == [1.0]

        # missing year (internal helper — generate_stats itself would throw on
        # missing years downstream regardless of the collector)
        c2 = AnnualCollector()
        dfm = DataFrame(water_year=Union{Int, Missing}[2000, missing, 2002], m=[1.0, 2.0, 3.0])
        StreamflowSignatures._collect_annual!(c2, dfm, ["m"], "water_year")
        @test c2.water_year == Int32[2000, 2002]
        @test c2.value == [1.0, 3.0]

        # missing VALUE becomes NaN (defensive; year kept)
        c3 = AnnualCollector()
        dfv = DataFrame(water_year=[2000, 2001, 2002], m=Union{Float64, Missing}[1.0, missing, 3.0])
        StreamflowSignatures._collect_annual!(c3, dfv, ["m"], "water_year")
        @test c3.water_year == Int32[2000, 2001, 2002]
        @test c3.value[1] == 1.0
        @test isnan(c3.value[2])
        @test c3.value[3] == 3.0
    end

    @testset "Collector respects explicit year_col" begin
        c = AnnualCollector()
        dfy = DataFrame(yr=[2001, 2002, 2003, 2004], m=[1.0, 2.0, 3.0, 4.0])
        generate_stats(dfy; value_cols=["m"], year_col="yr", collector=c)
        @test c.water_year == Int32[2001, 2002, 2003, 2004]
        @test c.signature == fill("m", 4)
    end

    @testset "Default collector=nothing leaves behavior unchanged" begin
        dfx = DataFrame(water_year=2000:2019, m=rand(20))
        s1 = generate_stats(dfx; value_cols=["m"])
        s2 = generate_stats(dfx; value_cols=["m"], collector=nothing)
        @test isequal(s1, s2)
    end
end

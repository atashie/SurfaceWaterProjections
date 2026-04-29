using Test
using Statistics
using Random

# Load the module
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using StreamflowSignatures

@testset "Changepoint Detection" begin

    @testset "Flat series — should select model 1" begin
        Random.seed!(42)
        years = Float64.(1980:2024)
        values = fill(5.0, 45) .+ randn(45) .* 0.1
        result = detect_changepoint(years, values)
        @test result.best_model == 1.0  # Flat
        @test isnan(result.cp_year)
        @test isnan(result.pre_slope)
        @test isnan(result.post_slope)
    end

    @testset "Strong trend — should select model 2 or 4 (not flat/step)" begin
        Random.seed!(42)
        years = Float64.(1980:2024)
        # Strong linear trend with small noise — should prefer trend or piecewise
        # (piecewise-continuous includes linear as a special case)
        values = collect(1.0:45.0) .+ randn(45) .* 0.5
        result = detect_changepoint(years, values)
        @test result.best_model in (2.0, 4.0)  # Trend or Piecewise (not flat or step)
    end

    @testset "Clear step change — should select model 3" begin
        Random.seed!(42)
        years = Float64.(1980:2024)
        values = vcat(fill(10.0, 22) .+ randn(22) .* 0.3,
                      fill(20.0, 23) .+ randn(23) .* 0.3)
        result = detect_changepoint(years, values)
        @test result.best_model == 3.0  # Step
        @test !isnan(result.cp_year)
        @test 1999.0 <= result.cp_year <= 2003.0  # Near 2001
        @test result.delta_bic > 0  # CP model preferred
        @test result.pre_mean < result.post_mean
        @test abs(result.pre_mean - 10.0) < 1.0
        @test abs(result.post_mean - 20.0) < 1.0
    end

    @testset "Trend reversal — should select model 4" begin
        Random.seed!(42)
        years = Float64.(1980:2024)
        n = length(years)
        mid = 22
        values = Float64[i <= mid ? Float64(i) : Float64(2*mid - i) for i in 1:n]
        values .+= randn(n) .* 0.3
        result = detect_changepoint(years, values)
        @test result.best_model == 4.0  # Piecewise
        @test !isnan(result.cp_year)
        @test result.pre_slope > 0  # Increasing before
        @test result.post_slope < 0  # Decreasing after
    end

    @testset "Too few observations — returns NaN" begin
        years = Float64.(1980:1994)  # 15 years, below min_total_obs=20
        values = randn(15)
        result = detect_changepoint(years, values)
        @test isnan(result.best_model)
        @test isnan(result.cp_year)
        @test isnan(result.delta_bic)
    end

    @testset "NaN handling — filters NaN pairs" begin
        Random.seed!(42)
        years = Float64.(1980:2024)
        values = fill(5.0, 45) .+ randn(45) .* 0.1
        values[5:10] .= NaN  # 6 NaN values, 39 remaining (> 20)
        result = detect_changepoint(years, values)
        @test !isnan(result.best_model)  # Should still work
    end

    @testset "Delta BIC sign convention" begin
        # Step series: delta_bic should be positive
        years = Float64.(1980:2024)
        values_step = vcat(fill(0.0, 22), fill(10.0, 23))
        result_step = detect_changepoint(years, values_step)
        @test result_step.delta_bic > 0  # CP preferred
    end

    @testset "Segment means for step model" begin
        years = Float64.(1980:2024)
        values = vcat(fill(5.0, 22), fill(15.0, 23))
        result = detect_changepoint(years, values)
        @test result.best_model == 3.0
        # τ can land at 2001 or 2002 (step boundary), so allow 1.0 tolerance
        @test abs(result.pre_mean - 5.0) < 1.0
        @test abs(result.post_mean - 15.0) < 1.0
    end

    @testset "Min segment obs respected" begin
        years = Float64.(1980:2024)
        # Changepoint at year 5 — pre-segment too small with min_segment_obs=10
        values = vcat(fill(0.0, 5), fill(10.0, 40))
        result = detect_changepoint(years, values; min_segment_obs=10)
        # CP should NOT be at year 5 (only 5 pre-obs)
        if result.best_model >= 3.0
            pre_count = count(y -> y <= result.cp_year, years)
            @test pre_count >= 10
        end
    end

    @testset "All 7 output fields present" begin
        years = Float64.(1980:2024)
        values = vcat(fill(5.0, 22), fill(15.0, 23))
        result = detect_changepoint(years, values)
        @test haskey(result, :best_model)
        @test haskey(result, :cp_year)
        @test haskey(result, :delta_bic)
        @test haskey(result, :pre_slope)
        @test haskey(result, :post_slope)
        @test haskey(result, :pre_mean)
        @test haskey(result, :post_mean)
    end

    @testset "best_model values are Float64" begin
        years = Float64.(1980:2024)
        values = fill(5.0, 45)
        result = detect_changepoint(years, values)
        @test result.best_model isa Float64
    end

    @testset "Noiseless step — RSS=0 perfect fit" begin
        years = Float64.(1980:2024)
        values = vcat(fill(5.0, 22), fill(15.0, 23))
        result = detect_changepoint(years, values)
        @test result.best_model == 3.0  # Step should win with perfect fit
        @test result.delta_bic > 0
        # Segment means should be exact
        @test result.pre_mean == 5.0 || abs(result.pre_mean - 5.0) < 1.0
        @test result.post_mean == 15.0 || abs(result.post_mean - 15.0) < 1.0
    end

    @testset "Noiseless constant — RSS=0 for flat" begin
        years = Float64.(1980:2024)
        values = fill(7.0, 45)
        result = detect_changepoint(years, values)
        @test result.best_model == 1.0  # Flat wins
        @test isnan(result.cp_year)
    end
end

# ---- Integration tests: generate_stats() with changepoint kwarg ----
using DataFrames

@testset "generate_stats() changepoint integration" begin

    cp_config = (start_year=1980, end_year=2024, min_total_obs=20, min_segment_obs=10)

    @testset "Pettitt columns appended to output" begin
        # Build a simple annual DataFrame
        years = collect(1980:2024)
        vals = vcat(fill(10.0, 22), fill(20.0, 23))
        df = DataFrame(water_year=years, metric=vals)

        result = generate_stats(df; value_cols=[:metric], year_col=:water_year, changepoint=cp_config)

        # Should have 8 standard stats + 8 Pettitt = 16 keys
        @test haskey(result, "metric_mean")
        @test haskey(result, "metric_senn_slp")
        # BIC columns should NOT exist
        @test !haskey(result, "metric_best_model")
        @test !haskey(result, "metric_cp_year")
        @test !haskey(result, "metric_cp_delta_bic")
        # Pettitt (8)
        @test haskey(result, "metric_pettitt_cp_year")
        @test haskey(result, "metric_pettitt_pval")
        @test haskey(result, "metric_pettitt_pre_mean")
        @test haskey(result, "metric_pettitt_post_mean")
        @test haskey(result, "metric_pettitt_delta_mean")
        @test haskey(result, "metric_pettitt_pct_change")
        @test haskey(result, "metric_pettitt_pre_mk_pval")
        @test haskey(result, "metric_pettitt_post_mk_pval")

        # Pettitt should detect a changepoint with low p-value
        @test !isnan(result["metric_pettitt_cp_year"])
        @test result["metric_pettitt_pval"] < 0.05

        # Pettitt differential metrics should be populated
        @test !isnan(result["metric_pettitt_delta_mean"])
        @test result["metric_pettitt_delta_mean"] ≈ 10.0 atol=2.0  # 20-10=10
    end

    @testset "Changepoint independent of trend_completeness" begin
        # Build data with many NaN years to fail trend_completeness
        years = collect(1980:2024)
        vals = Vector{Float64}(undef, 45)
        vals .= NaN
        # Only fill 15 values (below 80% of 45) — trends should be NaN
        # but changepoint should still run if >= 20 non-NA
        vals[1:12] .= 10.0  # 12 pre-CP
        vals[23:34] .= 20.0  # 12 post-CP (total 24 non-NA, but sparse)
        df = DataFrame(water_year=years, metric=vals)

        result = generate_stats(df; value_cols=[:metric], year_col=:water_year,
            trend_completeness=0.8, changepoint=cp_config)

        # Trend stats should be NaN (< 80% non-NA)
        @test isnan(result["metric_senn_slp"])
        @test isnan(result["metric_mk_rho"])

        # Mean/median should still be computed
        @test !isnan(result["metric_mean"])

        # Pettitt should still run (independent of trend completeness)
        @test !isnan(result["metric_pettitt_cp_year"])
    end

    @testset "Changepoint NaN when disabled (nothing)" begin
        years = collect(1980:2024)
        vals = vcat(fill(10.0, 22), fill(20.0, 23))
        df = DataFrame(water_year=years, metric=vals)

        result = generate_stats(df; value_cols=[:metric], year_col=:water_year, changepoint=nothing)

        # Standard stats should exist
        @test haskey(result, "metric_mean")
        # Changepoint columns should NOT exist
        @test !haskey(result, "metric_pettitt_cp_year")
        @test !haskey(result, "metric_pettitt_pval")
    end

    @testset "Changepoint respects WY window" begin
        # Data from 1960-2050, but changepoint window is 1980-2024
        years = collect(1960:2050)
        n = length(years)
        vals = fill(5.0, n)
        # Put a step change at 1970 (outside CP window)
        vals[1:11] .= 50.0
        df = DataFrame(water_year=years, metric=vals)

        result = generate_stats(df; value_cols=[:metric], year_col=:water_year, changepoint=cp_config)

        # Changepoint should not find the 1970 step (outside 1980-2024 window)
        # Within 1980-2024, data is constant at 5.0 — Pettitt p-value should be high
        @test result["metric_pettitt_pval"] > 0.05
    end

    @testset "Changepoint NaN-fill when too few rows" begin
        years = collect(1980:1989)  # Only 10 years
        vals = fill(5.0, 10)
        df = DataFrame(water_year=years, metric=vals)

        result = generate_stats(df; value_cols=[:metric], year_col=:water_year, changepoint=cp_config)

        # With min_rows=3 default, stats should compute but changepoint
        # needs 20 obs — should get NaN
        @test isnan(result["metric_pettitt_cp_year"])
        @test isnan(result["metric_pettitt_pval"])
    end
end

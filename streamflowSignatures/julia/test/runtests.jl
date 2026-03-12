using StreamflowSignatures
using Test
using DataFrames
using Dates
using Statistics

# Generate synthetic data for testing
function generate_synthetic_streamflow(n_years=25)
    # Generate daily dates for n_years
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years

    dates = [start_date + Day(i-1) for i in 1:n_days]

    # Create seasonal flow pattern with noise
    doy = [dayofyear(d) for d in dates]

    # Base flow with seasonal variation (peak in spring)
    base_flow = 2.0 .+ 1.5 .* sin.(2π .* (doy .- 100) ./ 365)
    noise = rand(n_days) .* 0.5
    Q = max.(0.01, base_flow .+ noise)

    df = DataFrame(
        gage_id = fill("TEST001", n_days),
        date = dates,
        Q = Q
    )

    return add_water_year_columns(df)
end

function generate_synthetic_with_climate(n_years=25)
    df = generate_synthetic_streamflow(n_years)

    # Add correlated precipitation
    PPT = max.(0.0, df.Q .* 2.5 .+ rand(nrow(df)))

    df.PPT = PPT
    return df
end


@testset "StreamflowSignatures.jl" begin

    @testset "Stats Module" begin
        @testset "Theil-Sen slope" begin
            x = collect(1.0:20.0)
            y = 2.0 .* x .+ 1.0 .+ randn(20) .* 0.5

            slope, intercept = theil_sen_slope(x, y)

            @test slope > 1.5
            @test slope < 2.5
        end

        @testset "Mann-Kendall test" begin
            y = collect(1.0:30.0) .+ randn(30)

            tau, pval = mann_kendall_test(y)

            @test tau > 0.7
            @test pval < 0.05
        end

        @testset "Spearman handles ties with average rank" begin
            x = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            y = Float64[3, 3, 3, 4, 4, 5, 5, 5, 6, 6]  # many ties

            rho, pval = StreamflowSignatures.spearman_correlation(x, y)

            # With ordinalrank (the bug), rho=1.0 because ties get unique ranks
            # With tiedrank (correct), rho≈0.9692 from Pearson cor of average ranks
            # This matches R's cor.test(method="spearman") which also uses average ranks
            @test abs(rho - 0.9692) < 0.01
            @test rho < 1.0  # Must NOT be 1.0 — that indicates the ordinalrank bug
            @test pval < 0.001
        end

        @testset "generate_stats output format" begin
            df = DataFrame(
                water_year = 2000:2019,
                metric1 = rand(20),
                metric2 = rand(20)
            )

            stats = generate_stats(df; value_cols=["metric1", "metric2"])

            expected_suffixes = [
                "_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                "_mk_rho", "_mk_pval", "_mean", "_median"
            ]

            for metric in ["metric1", "metric2"]
                for suffix in expected_suffixes
                    key = "$(metric)$(suffix)"
                    @test haskey(stats, key)
                    @test stats[key] isa Real
                end
            end
        end
    end

    @testset "Flow Volumes" begin
        df = generate_synthetic_streamflow()
        results = calculate_flow_vols_by_year(df)

        @test haskey(results, "Qann_mean")
        @test haskey(results, "Qann_senn_slp")
        @test haskey(results, "Q50_mean")

        @test results["Qann_mean"] > 0
        @test results["Q50_mean"] > 0
    end

    @testset "Flashiness" begin
        df = generate_synthetic_streamflow()
        results = analyze_flashiness_trends(df)

        @test haskey(results, "flashinessRB_mean")
        @test haskey(results, "flashinessRB_senn_slp")

        if !isnan(results["flashinessRB_mean"])
            @test results["flashinessRB_mean"] >= 0
        end
    end

    @testset "Timing" begin
        df = generate_synthetic_streamflow()
        results = analyze_flow_timing_trends(df)

        @test haskey(results, "D50_day_mean")
        @test haskey(results, "Dmax_mean")

        # Check ordering
        d10 = results["D10_day_mean"]
        d50 = results["D50_day_mean"]
        d90 = results["D90_day_mean"]

        if !any(isnan.([d10, d50, d90]))
            @test d10 < d50 < d90
        end

        @testset "D25_to_D75 computed" begin
            result = analyze_flow_timing_trends(df)
            # D25_to_D75 must NOT be NaN — it should be a positive number
            @test haskey(result, "D25_to_D75_mean")
            @test !isnan(result["D25_to_D75_mean"])
            @test result["D25_to_D75_mean"] > 0
        end
    end

    @testset "FDC" begin
        df = generate_synthetic_streamflow()
        results = analyze_fdc_trends(df)

        @test haskey(results, "FDCall_mean")
        @test haskey(results, "FDC90th_mean")
        @test haskey(results, "FDCmid_mean")
    end

    @testset "Baseflow" begin
        @testset "Eckhardt filter" begin
            Q = max.(0.01, randn(1000) .+ 2.0)
            bf = eckhardt_filter(Q)

            @test all(bf .>= 0)
            @test all(bf .<= Q .+ 1e-10)
        end

        @testset "Lyne-Hollick filter" begin
            Q = max.(0.01, randn(1000) .+ 2.0)
            bf = lyne_hollick_filter(Q)

            @test all(bf .>= 0)
            @test all(bf .<= Q .+ 1e-10)
        end

        @testset "Baseflow indices" begin
            df = generate_synthetic_streamflow()
            results = analyze_baseflow_indices(df)

            @test haskey(results, "BFI_Eckhardt_mean")
            @test haskey(results, "BFI_LyneHollick_mean")

            bfi_e = results["BFI_Eckhardt_mean"]
            bfi_lh = results["BFI_LyneHollick_mean"]

            if !isnan(bfi_e)
                @test 0 <= bfi_e <= 1
            end
            if !isnan(bfi_lh)
                @test 0 <= bfi_lh <= 1
            end
        end
    end

    @testset "Recession" begin
        df = generate_synthetic_streamflow()
        results = analyze_recession_parameters(df)

        @test haskey(results, "log_a_pointcloud_mean")
        @test haskey(results, "b_pointcloud_mean")

        @testset "Seasonality minimum matches R formula" begin
            # Known test case: cos curve with minimum at day 200
            n = 1000
            dowy = rand(1:365, n)
            omega = 2 * pi / 365
            log_a = -0.5 .* cos.(omega .* (dowy .- 200)) .+ randn(n) .* 0.01

            result = StreamflowSignatures.fit_recession_seasonality(
                Float64.(log_a), Int.(dowy)
            )

            # The minimum day should be near 200, NOT near 365-200=165
            @test abs(result["minimum_day"] - 200) < 30  # within 30 days
        end
    end

    @testset "Pulses" begin
        df = generate_synthetic_streamflow()
        results = calculate_pulse_metrics(df)

        @test haskey(results, "n_high_pulses_year_mean")
        @test haskey(results, "TQmean_mean")
        @test haskey(results, "Flow_Reversals_annual_mean")

        tqmean = results["TQmean_mean"]
        if !isnan(tqmean)
            @test 0 <= tqmean <= 100
        end
    end

    @testset "Climate Signatures" begin
        df = generate_synthetic_with_climate()

        @testset "Runoff Ratios" begin
            results = analyze_Q_PPT_relationships(df)

            @test haskey(results, "annual_runoff_ratio_mean")
        end

        @testset "Elasticity" begin
            results = calculate_streamflow_elasticity(df)

            @test haskey(results, "elasticity_static")
        end

        @testset "Q-P Seasonality" begin
            results = calculate_qp_seasonality(df)

            @test haskey(results, "qp_slope_sd_mean")
            @test haskey(results, "qp_bimodality_mean")
        end

        @testset "Average Storage" begin
            results = calculate_average_storage(df)

            @test haskey(results, "avg_storage_mean")
        end
    end

    @testset "Integration" begin
        df = generate_synthetic_streamflow()

        # Collect all non-climate signatures
        all_sigs = Dict{String, Float64}()
        merge!(all_sigs, calculate_flow_vols_by_year(df))
        merge!(all_sigs, analyze_flashiness_trends(df))
        merge!(all_sigs, analyze_flow_timing_trends(df))
        merge!(all_sigs, analyze_fdc_trends(df))
        merge!(all_sigs, analyze_baseflow_indices(df))
        merge!(all_sigs, analyze_recession_parameters(df))
        merge!(all_sigs, calculate_pulse_metrics(df))

        @test length(all_sigs) > 300

        # Count non-NaN values
        non_nan = count(!isnan, values(all_sigs))
        @test non_nan > 200
    end
end

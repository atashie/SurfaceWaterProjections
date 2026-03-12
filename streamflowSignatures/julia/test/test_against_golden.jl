"""
Test Julia signatures against R golden outputs.

This test suite validates that the Julia implementation produces
results within tolerance of the canonical R implementation.
"""

using StreamflowSignatures
using Test
using DataFrames
using CSV

# Configuration
const RELATIVE_TOLERANCE = 0.0001  # 0.01%
const ABSOLUTE_TOLERANCE = 1e-10

# Paths (relative to test directory)
const BASE_DIR = dirname(dirname(@__FILE__))
const TEST_DATA_DIR = joinpath(BASE_DIR, "..", "test-data")
const GOLDEN_OUTPUT_DIR = joinpath(BASE_DIR, "..", "golden-outputs")


function relative_error(actual::Real, expected::Real)
    if isnan(actual) && isnan(expected)
        return 0.0
    end
    if isnan(actual) || isnan(expected)
        return Inf
    end
    if abs(expected) < ABSOLUTE_TOLERANCE
        if abs(actual) < ABSOLUTE_TOLERANCE
            return 0.0
        end
        return Inf
    end
    return abs(actual - expected) / abs(expected)
end


function load_test_data()
    streamflow_path = joinpath(TEST_DATA_DIR, "sample_streamflow.parquet")
    if !isfile(streamflow_path)
        return nothing, nothing
    end

    streamflow = read_parquet(streamflow_path)
    streamflow = add_water_year_columns(streamflow)

    climate_path = joinpath(TEST_DATA_DIR, "sample_climate.parquet")
    climate = nothing
    if isfile(climate_path)
        climate = read_parquet(climate_path)
        climate = add_water_year_columns(climate)
    end

    return streamflow, climate
end


function load_golden_outputs()
    golden_path = joinpath(GOLDEN_OUTPUT_DIR, "expected_signatures.csv")
    if !isfile(golden_path)
        return nothing
    end
    return CSV.read(golden_path, DataFrame)
end


function get_gage_data(df, gage_id)
    return df[df.gage_id .== gage_id, :]
end


function get_gage_golden(golden, gage_id)
    rows = golden[string.(golden.gage_id) .== string(gage_id), :]
    if nrow(rows) == 0
        return nothing
    end
    return Dict(pairs(rows[1, :]))
end


@testset "Golden Output Validation" begin
    streamflow, climate = load_test_data()
    golden = load_golden_outputs()

    if streamflow === nothing
        @warn "Test data not found - skipping golden validation tests"
        @test_skip true
        return
    end

    if golden === nothing
        @warn "Golden outputs not found - skipping golden validation tests"
        @test_skip true
        return
    end

    gages = unique(streamflow.gage_id)

    @testset "Flow Volumes" begin
        for gage_id in gages[1:min(5, length(gages))]
            gage_data = get_gage_data(streamflow, gage_id)
            gage_golden = get_gage_golden(golden, gage_id)

            if gage_golden === nothing
                continue
            end

            results = calculate_flow_vols_by_year(gage_data)

            for metric in ["Qann_mean", "Qann_senn_slp", "Q50_mean"]
                if haskey(gage_golden, Symbol(metric)) && haskey(results, metric)
                    error = relative_error(results[metric], gage_golden[Symbol(metric)])
                    @test error <= RELATIVE_TOLERANCE
                end
            end
        end
    end

    @testset "Flashiness" begin
        for gage_id in gages[1:min(5, length(gages))]
            gage_data = get_gage_data(streamflow, gage_id)
            gage_golden = get_gage_golden(golden, gage_id)

            if gage_golden === nothing
                continue
            end

            results = analyze_flashiness_trends(gage_data)

            for metric in ["flashinessRB_mean", "flashinessRB_senn_slp"]
                if haskey(gage_golden, Symbol(metric)) && haskey(results, metric)
                    error = relative_error(results[metric], gage_golden[Symbol(metric)])
                    @test error <= RELATIVE_TOLERANCE
                end
            end
        end
    end

    @testset "Timing" begin
        for gage_id in gages[1:min(5, length(gages))]
            gage_data = get_gage_data(streamflow, gage_id)
            gage_golden = get_gage_golden(golden, gage_id)

            if gage_golden === nothing
                continue
            end

            results = analyze_flow_timing_trends(gage_data)

            for metric in ["D50_day_mean", "Dmax_mean"]
                if haskey(gage_golden, Symbol(metric)) && haskey(results, metric)
                    error = relative_error(results[metric], gage_golden[Symbol(metric)])
                    @test error <= RELATIVE_TOLERANCE
                end
            end
        end
    end

    @testset "FDC" begin
        for gage_id in gages[1:min(5, length(gages))]
            gage_data = get_gage_data(streamflow, gage_id)
            gage_golden = get_gage_golden(golden, gage_id)

            if gage_golden === nothing
                continue
            end

            results = analyze_fdc_trends(gage_data)

            for metric in ["FDCall_mean", "FDC90th_mean"]
                if haskey(gage_golden, Symbol(metric)) && haskey(results, metric)
                    error = relative_error(results[metric], gage_golden[Symbol(metric)])
                    @test error <= RELATIVE_TOLERANCE
                end
            end
        end
    end

    @testset "Baseflow" begin
        for gage_id in gages[1:min(5, length(gages))]
            gage_data = get_gage_data(streamflow, gage_id)
            gage_golden = get_gage_golden(golden, gage_id)

            if gage_golden === nothing
                continue
            end

            results = analyze_baseflow_indices(gage_data)

            for metric in ["BFI_Eckhardt_mean", "BFI_LyneHollick_mean"]
                if haskey(gage_golden, Symbol(metric)) && haskey(results, metric)
                    error = relative_error(results[metric], gage_golden[Symbol(metric)])
                    @test error <= RELATIVE_TOLERANCE
                end
            end
        end
    end
end

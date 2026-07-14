# Tests for the area_normalized gate on Q-to-PPT signatures (July 2026).
#
# Gages with no drainage area (HYDAT gap — canals, dam outflows, channel splits)
# carry Q in raw m³/s instead of mm/day. Q-to-PPT signatures (runoff ratios,
# elasticity, Q-P seasonality, storage) mix Q and PPT units and are therefore
# skipped when area_normalized=false. Q-only signatures must be unaffected.

using StreamflowSignatures
using Test
using DataFrames
using Dates
using Logging
using Random

# Prefixes of every Q-to-PPT (climate-dependent) output column
const CLIMATE_SIGNATURE_PREFIXES = [
    "annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
    "summer_runoff_ratio", "fall_runoff_ratio", "runoff_ratio_high_count",
    "elasticity_static", "elasticity_rolling", "elasticity_annual",
    "elasticity_years_total", "elasticity_years_low_ppt",
    "qp_slope_sd", "qp_bimodality",
    "avg_storage",
]

# Representative Q-only columns that must survive the gate
const QONLY_SENTINELS = [
    "Qann_mean", "Q50_mean", "BFI_Eckhardt_mean", "flashinessRB_mean",
    "D50_day_mean", "TQmean_mean", "FDCall_mean",
]

function gate_synthetic_gage(n_years=25; seed=7)
    Random.seed!(seed)
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years
    dates = [start_date + Day(i-1) for i in 1:n_days]
    doy = [dayofyear(d) for d in dates]
    base_flow = 2.0 .+ 1.5 .* sin.(2π .* (doy .- 100) ./ 365)
    Q = max.(0.01, base_flow .+ rand(n_days) .* 0.5)
    df = DataFrame(gage_id=fill("TESTGATE1", n_days), date=dates, Q=Q)
    df = add_water_year_columns(df)
    df.PPT = max.(0.0, df.Q .* 2.5 .+ rand(nrow(df)))
    return df
end

is_climate_key(k::String) = any(startswith(k, p) for p in CLIMATE_SIGNATURE_PREFIXES)

@testset "Area-normalized gate on Q-to-PPT signatures" begin
    df = gate_synthetic_gage()

    @testset "Default (area_normalized=true) computes climate signatures" begin
        r = calculate_all_signatures(df, true)
        @test any(is_climate_key(k) for k in keys(r))
        @test haskey(r, "annual_runoff_ratio_mean")
        @test haskey(r, "avg_storage_mean")
        for k in QONLY_SENTINELS
            @test haskey(r, k)
        end
    end

    @testset "area_normalized=false skips ALL climate signatures" begin
        r = calculate_all_signatures(df, true; area_normalized=false)
        climate_keys = [k for k in keys(r) if is_climate_key(k)]
        @test isempty(climate_keys)
        # Q-only signatures unaffected
        for k in QONLY_SENTINELS
            @test haskey(r, k)
        end
    end

    @testset "Gated result identical to gage without climate data" begin
        r_gated = calculate_all_signatures(df, true; area_normalized=false)
        df_noclimate = select(df, Not(:PPT))
        r_noclimate = calculate_all_signatures(df_noclimate, false)
        @test isequal(r_gated, r_noclimate)
    end

    @testset "Gate is inert when area_normalized=true" begin
        r_default = calculate_all_signatures(df, true)
        r_explicit = calculate_all_signatures(df, true; area_normalized=true)
        @test isequal(r_default, r_explicit)
    end

    @testset "Collector receives no climate annual values when gated" begin
        c = AnnualCollector()
        calculate_all_signatures(df, true; area_normalized=false, collector=c)
        @test !any(is_climate_key(s) for s in unique(c.signature))
        # Q-only annual values are still collected
        @test "Qann" in unique(c.signature)
    end

    @testset "Gate skip logs info, not warnings" begin
        @test_logs (:info,) min_level=Logging.Info match_mode=:any calculate_all_signatures(df, true; area_normalized=false)
        @test_logs min_level=Logging.Warn calculate_all_signatures(df, true; area_normalized=false)
    end
end

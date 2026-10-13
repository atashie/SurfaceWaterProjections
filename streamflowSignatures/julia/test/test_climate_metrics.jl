"""
Climate metrics (guidelines §3.9; decisions 2026-10-12) — closed-form unit tests.

Every expected value below is derived by hand from the definitions in
docs/plans/2026-10-12-climate-metrics-plan.md, never from the implementation. Annual
values are read back through the AnnualCollector (collected BEFORE any gating), so a
single hand-built water year is enough to pin a definition.
"""

using StreamflowSignatures
using Test
using DataFrames
using Dates
using Statistics
using Logging

# One complete water year (Oct 1 of wy-1 .. Sep 30 of wy) with optional climate columns.
function climate_year_frame(wy::Int; P=nothing, TN=nothing, TX=nothing, Q=nothing)
    dates = collect(Date(wy - 1, 10, 1):Day(1):Date(wy, 9, 30))
    n = length(dates)
    df = DataFrame(gage_id = fill("TCLIM", n), date = dates,
                   Q = Q === nothing ? fill(1.0, n) : Q)
    df = add_water_year_columns(df)
    P !== nothing && (df.PPT = Float64.(P))
    TN !== nothing && (df.TMIN = Float64.(TN))
    TX !== nothing && (df.TMAX = Float64.(TX))
    return df
end

# Deterministic multi-year frame: P from the month number plus a small dowy-based term,
# temperatures seasonal — no randomness, so expectations are reproducible.
function climate_multi_year(wys; dry_every=7)
    frames = DataFrame[]
    for wy in wys
        base = climate_year_frame(wy)
        n = nrow(base)
        P = [(d % dry_every == 0) ? 0.0 : Float64(base.month[d]) + 0.01 * d for d in 1:n]
        TX = [20.0 + 15.0 * sin(2π * (d - 100) / 365) for d in 1:n]
        TN = TX .- 12.0
        push!(frames, climate_year_frame(wy; P=P, TN=TN, TX=TX))
    end
    return vcat(frames...)
end

function annual_lookup(c::AnnualCollector, sig::String, wy::Int)
    for i in eachindex(c.value)
        if c.signature[i] == sig && c.water_year[i] == wy
            return c.value[i]
        end
    end
    return nothing
end

function collect_precip(df)
    c = AnnualCollector()
    r = calculate_precipitation_metrics(df; collector=c)
    return r, c
end

function collect_temp(df; kwargs...)
    c = AnnualCollector()
    r = calculate_temperature_metrics(df; collector=c, kwargs...)
    return r, c
end

const PRECIP_STAT_KEYS = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                          "_mk_rho", "_mk_pval", "_mean", "_median"]

@testset "Climate metrics (guidelines 3.9)" begin

    @testset "Registry: 20 precipitation + 19 temperature bases, config-derived names" begin
        @test length(PRECIP_METRICS) == 20
        @test length(TEMP_METRICS) == 19
        @test length(unique(vcat(PRECIP_METRICS, TEMP_METRICS))) == 39
        @test "prcp_5d_max" in PRECIP_METRICS           # multi_day_max_window_days = 5
        @test "prcp_50_dowy" in PRECIP_METRICS          # cumulative_fraction = 0.5
        @test "prcp_ndays_abv_10mm" in PRECIP_METRICS   # wet_day_thresholds_mm = [10, 25, 50]
        @test "prcp_ndays_abv_25mm" in PRECIP_METRICS
        @test "prcp_ndays_abv_50mm" in PRECIP_METRICS
        for s in ["win", "spr", "sum", "fal"]             # the Qwin/Qspr/Qsum/Qfal suffixes
            @test "prcp_$(s)" in PRECIP_METRICS
            @test "temp_max_$(s)" in TEMP_METRICS
            @test "temp_min_$(s)" in TEMP_METRICS
            @test "gdd_cum_$(s)" in TEMP_METRICS
        end
        @test !("prcp_sumr" in PRECIP_METRICS) && !("prcp_fall" in PRECIP_METRICS)
    end

    @testset "Single-storm year: timing, concentration and spells by hand" begin
        wy = 2001                                # non-leap: 365 days
        P = zeros(365); P[100] = 100.0           # dowy 100 = Jan 8 2001 (month 1)
        df = climate_year_frame(wy; P=P)
        r, c = collect_precip(df)
        v(sig) = annual_lookup(c, sig, wy)
        @test v("prcp_annual") == 100.0
        @test v("prcp_1d_max") == 100.0
        @test v("prcp_5d_max") == 100.0
        @test v("prcp_max_dowy") == 100.0
        @test v("prcp_50_dowy") == 100.0
        @test v("prcp_cm_dowy") == 100.0
        @test v("prcp_ndays_abv_10mm") == 1.0
        @test v("prcp_ndays_abv_25mm") == 1.0
        @test v("prcp_ndays_abv_50mm") == 1.0
        @test isapprox(v("prcp_annual_sd"), std(P); rtol=1e-12)
        @test isapprox(v("prcp_annual_cv"), std(P) / mean(P); rtol=1e-12)
        # Gini with one wet day = (n - 1)/n; seasonality index with one wet month = 11/6
        @test isapprox(v("prcp_gini"), 364 / 365; rtol=1e-12)
        @test isapprox(v("prcp_seasonality_index"), 11 / 6; rtol=1e-12)
        # Dry runs: dowy 1-99 (99 days) and 101-365 (265 days) -> the longer one
        @test v("prcp_intermit_annual") == 265.0
        @test v("prcp_intermit_annual_dowy_start") == 101.0
        @test v("prcp_intermit_annual_dowy_end") == 365.0
        # Seasonal totals: the storm sits in January (winter)
        @test v("prcp_win") == 100.0
        @test v("prcp_spr") == 0.0 && v("prcp_sum") == 0.0 && v("prcp_fal") == 0.0
        # Key-set contract: 8 statistics per base (no changepoint here)
        for m in PRECIP_METRICS, s in PRECIP_STAT_KEYS
            @test haskey(r, "$(m)$(s)")
        end
        @test length(r) == 20 * 8
    end

    @testset "Uniform year: zero concentration, mid-year timing, no dry spell" begin
        wy = 2001
        P = fill(2.0, 365)
        df = climate_year_frame(wy; P=P)
        _, c = collect_precip(df)
        v(sig) = annual_lookup(c, sig, wy)
        @test v("prcp_annual") == 730.0
        @test isapprox(v("prcp_gini"), 0.0; atol=1e-12)
        # SI from the month lengths of WY 2001 (Oct 2000 .. Sep 2001)
        mlen = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # Jan..Dec of the WY's months
        X = 2.0 .* mlen
        R = sum(X)
        @test isapprox(v("prcp_seasonality_index"), sum(abs.(X .- R / 12)) / R; rtol=1e-12)
        # cum/total >= 0.5 first on day 183 (2*182/730 = 0.4986 < 0.5 <= 2*183/730)
        @test v("prcp_50_dowy") == 183.0
        @test isapprox(v("prcp_cm_dowy"), 183.0; rtol=1e-12)      # Σ d·2 / 730
        @test v("prcp_5d_max") == 10.0
        @test v("prcp_max_dowy") == 1.0                           # ties -> first day
        @test v("prcp_intermit_annual") == 0.0                    # no dry day
        @test isnan(v("prcp_intermit_annual_dowy_start"))
        @test isnan(v("prcp_intermit_annual_dowy_end"))
        @test v("prcp_annual_sd") == 0.0 && v("prcp_annual_cv") == 0.0
        @test v("prcp_ndays_abv_10mm") == 0.0
    end

    @testset "Thresholds are strict/inclusive as decided: dry <= 0.1 mm, wet > t" begin
        wy = 2001
        P = fill(5.0, 365)
        P[10] = 0.1          # dry (<= 0.1)
        P[11] = 0.1000001    # wet
        P[12] = 0.0          # dry
        P[13] = 0.1          # dry -> but separated from day 10 by the wet day 11
        P[50] = 10.0         # NOT above 10
        P[51] = 10.0001      # above 10
        P[52] = 25.0         # not above 25
        P[53] = 50.0         # not above 50
        P[54] = 50.5         # above all three
        df = climate_year_frame(wy; P=P)
        _, c = collect_precip(df)
        v(sig) = annual_lookup(c, sig, wy)
        @test v("prcp_ndays_abv_10mm") == 4.0   # 10.0001, 25, 50, 50.5
        @test v("prcp_ndays_abv_25mm") == 2.0   # 50, 50.5
        @test v("prcp_ndays_abv_50mm") == 1.0   # 50.5
        # Longest dry run is days 12-13 (2 days); day 10 alone is a 1-day run
        @test v("prcp_intermit_annual") == 2.0
        @test v("prcp_intermit_annual_dowy_start") == 12.0
        @test v("prcp_intermit_annual_dowy_end") == 13.0
    end

    @testset "Ties take the first occurrence; the 5-day window stays inside the year" begin
        wy = 2001
        P = zeros(365)
        P[1:5] = [10.0, 0.0, 0.0, 0.0, 0.0]
        P[361:365] = [1.0, 2.0, 3.0, 4.0, 5.0]   # window 361-365 sums to 15
        P[200] = 10.0                           # a second 10 -> max dowy stays 1
        df = climate_year_frame(wy; P=P)
        _, c = collect_precip(df)
        v(sig) = annual_lookup(c, sig, wy)
        @test v("prcp_max_dowy") == 1.0
        @test v("prcp_5d_max") == 15.0
        # Two equally long dry runs (6-199 is longer than 201-360: 194 vs 160 days);
        # make them equal and check the FIRST wins
        P2 = zeros(365); P2[1] = 1.0; P2[101] = 1.0; P2[201] = 1.0   # runs 2-100 and 102-200 (99 each), 202-365 (164)
        P2[301] = 1.0                                                 # splits the tail: 202-300 (99), 302-365 (64)
        df2 = climate_year_frame(wy; P=P2)
        _, c2 = collect_precip(df2)
        @test annual_lookup(c2, "prcp_intermit_annual", wy) == 99.0
        @test annual_lookup(c2, "prcp_intermit_annual_dowy_start", wy) == 2.0
        @test annual_lookup(c2, "prcp_intermit_annual_dowy_end", wy) == 100.0
    end

    @testset "Seasonal totals use the Q seasonal convention (fall = Oct+Nov+Sep)" begin
        wy = 2001
        base = climate_year_frame(wy)
        P = Float64.(base.month)                 # P = calendar month number
        df = climate_year_frame(wy; P=P)
        _, c = collect_precip(df)
        v(sig) = annual_lookup(c, sig, wy)
        @test v("prcp_win") == 12 * 31 + 1 * 31 + 2 * 28       # Dec 2000, Jan, Feb 2001
        @test v("prcp_spr") == 3 * 31 + 4 * 30 + 5 * 31
        @test v("prcp_sum") == 6 * 30 + 7 * 31 + 8 * 31
        @test v("prcp_fal") == 9 * 30 + 10 * 31 + 11 * 30      # Sep 2001 + Oct, Nov 2000
        @test v("prcp_annual") == v("prcp_win") + v("prcp_spr") + v("prcp_sum") + v("prcp_fal")
    end

    @testset "Grid guard: leap year accepted, short year and NaN year left NaN" begin
        df_leap = climate_year_frame(2000; P=fill(1.0, 366))   # Oct 1999 - Sep 2000 has Feb 29
        _, c = collect_precip(df_leap)
        @test annual_lookup(c, "prcp_annual", 2000) == 366.0
        @test annual_lookup(c, "prcp_50_dowy", 2000) == 183.0

        df_short = climate_year_frame(2001; P=fill(1.0, 365))[1:364, :]
        _, c2 = collect_precip(df_short)
        @test isnan(annual_lookup(c2, "prcp_annual", 2001))

        P_nan = fill(1.0, 365); P_nan[50] = NaN
        _, c3 = collect_precip(climate_year_frame(2001; P=P_nan))
        @test isnan(annual_lookup(c3, "prcp_annual", 2001))
        @test isnan(annual_lookup(c3, "prcp_gini", 2001))
    end

    @testset "Temperature: extremes, seasons, frost window, GDD clamp" begin
        wy = 2001
        base = climate_year_frame(wy)
        TX = Float64.(base.month)                 # TMAX = month number
        TN = TX .- 30.0                           # TMIN = month - 30 (always < 0)
        df = climate_year_frame(wy; TN=TN, TX=TX)
        _, c = collect_temp(df)
        v(sig) = annual_lookup(c, sig, wy)
        @test v("temp_max_annual") == 12.0
        @test v("temp_max_annual_dowy") == 62.0   # Dec 1 2000 = dowy 31 + 30 + 1
        @test v("temp_min_annual") == 1.0 - 30.0
        @test v("temp_min_annual_dowy") == 93.0   # Jan 1 2001
        @test v("temp_max_win") == 12.0 && v("temp_max_spr") == 5.0
        @test v("temp_max_sum") == 8.0 && v("temp_max_fal") == 11.0
        @test v("temp_min_win") == 1.0 - 30.0 && v("temp_min_spr") == 3.0 - 30.0
        @test v("temp_min_sum") == 6.0 - 30.0 && v("temp_min_fal") == 9.0 - 30.0
        @test v("frost_first_dowy") == 1.0        # below zero from Oct 1
        @test v("frost_last_dowy") == 273.0       # Jun 30 2001 = dowy 273
        @test v("gdd_cum_annual") == 0.0          # mean temp = month - 15 <= -3 everywhere

        # Explicit frost days: dowy 20, 250 (Jun 7) and 300 (Jul 27, after the June cut-off)
        TN2 = fill(5.0, 365); TN2[20] = -1.0; TN2[250] = -0.5; TN2[300] = -2.0
        TX2 = fill(15.0, 365)
        _, c2 = collect_temp(climate_year_frame(wy; TN=TN2, TX=TX2))
        @test annual_lookup(c2, "frost_first_dowy", wy) == 20.0
        @test annual_lookup(c2, "frost_last_dowy", wy) == 250.0
        # Jun 30 (dowy 273) counts, Jul 1 (dowy 274) does not
        TN3 = fill(5.0, 365); TN3[273] = -1.0; TN3[274] = -1.0
        _, c3 = collect_temp(climate_year_frame(wy; TN=TN3, TX=TX2))
        @test annual_lookup(c3, "frost_first_dowy", wy) == 273.0
        @test annual_lookup(c3, "frost_last_dowy", wy) == 273.0
        # Frost-free year -> both NaN; exactly 0 is NOT frost (strict <)
        TN4 = fill(0.0, 365)
        _, c4 = collect_temp(climate_year_frame(wy; TN=TN4, TX=TX2))
        @test isnan(annual_lookup(c4, "frost_first_dowy", wy))
        @test isnan(annual_lookup(c4, "frost_last_dowy", wy))

        # GDD: (30 + 10)/2 = 20 per day over 365 days; a cold year clamps to 0
        _, c5 = collect_temp(climate_year_frame(wy; TN=fill(10.0, 365), TX=fill(30.0, 365)))
        @test annual_lookup(c5, "gdd_cum_annual", wy) == 20.0 * 365
        @test annual_lookup(c5, "gdd_cum_win", wy) == 20.0 * (31 + 31 + 28)
        @test annual_lookup(c5, "gdd_cum_fal", wy) == 20.0 * (31 + 30 + 30)
        _, c6 = collect_temp(climate_year_frame(wy; TN=fill(-20.0, 365), TX=fill(10.0, 365)))
        @test annual_lookup(c6, "gdd_cum_annual", wy) == 0.0
        # Partial clamp: tmean = -5 on half the days, +5 on the other half
        TNm = [d <= 180 ? -10.0 : 0.0 for d in 1:365]
        TXm = [d <= 180 ? 0.0 : 10.0 for d in 1:365]
        _, c7 = collect_temp(climate_year_frame(wy; TN=TNm, TX=TXm))
        @test annual_lookup(c7, "gdd_cum_annual", wy) == 5.0 * 185

        # A NaN in either series voids the year
        TXn = copy(TX2); TXn[7] = NaN
        _, c8 = collect_temp(climate_year_frame(wy; TN=TN2, TX=TXn))
        @test isnan(annual_lookup(c8, "temp_max_annual", wy))
        @test isnan(annual_lookup(c8, "frost_first_dowy", wy))
    end

    @testset "Record-anchored decade gate on the frost metrics (snow precedent)" begin
        # 20 years; frost every year in the first decade, never in the last
        frames = DataFrame[]
        for (k, wy) in enumerate(2001:2020)
            n = length(Date(wy - 1, 10, 1):Day(1):Date(wy, 9, 30))
            TN = fill(5.0, n)
            if k <= 10
                TN[30 + k] = -1.0     # frost moves a day later each year -> a trend exists
            end
            push!(frames, climate_year_frame(wy; TN=TN, TX=fill(15.0, n)))
        end
        df = vcat(frames...)
        r_gated = calculate_temperature_metrics(df; trend_completeness=0.6, decade_completeness=0.8)
        # Anchor years 2001-2020 (all temperature-valid); last decade 0/10 frost -> gate
        @test isnan(r_gated["frost_first_dowy_senn_slp"])
        @test isnan(r_gated["frost_first_dowy_mk_pval"])
        @test !isnan(r_gated["frost_first_dowy_mean"])        # mean/median still computed
        @test isapprox(r_gated["frost_first_dowy_mean"], mean(31:40); atol=1e-12)
        @test !isnan(r_gated["temp_max_annual_senn_slp"])      # non-gated metrics unaffected
        # Without the decade knob the gate is off: the 10-value series carries its trend
        r_open = calculate_temperature_metrics(df)
        @test !isnan(r_open["frost_first_dowy_senn_slp"])
        @test isapprox(r_open["frost_first_dowy_senn_slp"], 1.0; atol=1e-12)
    end

    @testset "Multi-year statistics: 16 fields with changepoint, stats floor respected" begin
        df = climate_multi_year(1991:2020)
        cp = (start_year=1980, end_year=2024, min_total_obs=20, min_segment_obs=10)
        rp = calculate_precipitation_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                             min_values_for_stats=20, changepoint=cp)
        rt = calculate_temperature_metrics(df; trend_completeness=0.6, decade_completeness=0.8,
                                           min_values_for_stats=20, changepoint=cp)
        @test length(rp) == 20 * 16
        @test length(rt) == 19 * 16
        for m in PRECIP_METRICS
            @test haskey(rp, "$(m)_pettitt_cp_year")
        end
        @test !isnan(rp["prcp_annual_mean"]) && rp["prcp_annual_mean"] > 0
        @test 0.0 <= rp["prcp_gini_mean"] <= 1.0
        @test 0.0 <= rp["prcp_seasonality_index_mean"] <= 11 / 6
        @test !isnan(rt["gdd_cum_annual_mean"]) && rt["gdd_cum_annual_mean"] > 0
        # Identity: the four seasonal GDD sums add up to the annual sum, per year
        c = AnnualCollector()
        calculate_temperature_metrics(df; collector=c)
        for wy in 1991:2020
            s = sum(annual_lookup(c, "gdd_cum_$(k)", wy) for k in ["win", "spr", "sum", "fal"])
            @test isapprox(s, annual_lookup(c, "gdd_cum_annual", wy); rtol=1e-12)
        end
        # Stats floor: 10 years -> every statistic NaN
        r_short = calculate_precipitation_metrics(climate_multi_year(2011:2020); min_values_for_stats=20)
        @test all(isnan, values(r_short))
    end

    @testset "Preprocessor: temperature validity mirrors the SWE policy without a sign check" begin
        frames = DataFrame[]
        for wy in 2001:2004
            n = length(Date(wy - 1, 10, 1):Day(1):Date(wy, 9, 30))
            P = fill(2.0, n)
            TN = fill(-15.0, n)                      # negative temperatures are fine
            TX = fill(5.0, n)
            Q = fill(1.0, n)
            if wy == 2002
                TN[100:101] .= NaN                   # 2-day internal gap -> filled
            elseif wy == 2003
                TX[50:85] .= NaN                     # 36 NaN > 30 -> temperature-invalid
            elseif wy == 2004
                TN[200:203] .= NaN                   # 4-day gap > 3 -> residual NaN -> invalid
            end
            df = climate_year_frame(wy; P=P, TN=TN, TX=TX, Q=Q)
            push!(frames, df)
        end
        df = vcat(frames...)
        pp = preprocess_daily_data(df)
        @test pp.valid_years == [2001, 2002, 2003, 2004]
        @test pp.valid_climate_years == [2001, 2002, 2003, 2004]
        @test pp.valid_temp_years == [2001, 2002]
        @test "TMIN" in names(pp.data) && "TMAX" in names(pp.data)
        # The 2-day gap was linearly filled
        y2002 = pp.data[pp.data.water_year .== 2002, :]
        @test all(!isnan, Float64.(coalesce.(y2002.TMIN, NaN)))
        @test y2002.TMIN[100] == -15.0

        # Without TMIN/TMAX the result is unchanged and valid_temp_years is empty
        pp0 = preprocess_daily_data(select(df, Not([:TMIN, :TMAX])))
        @test pp0.valid_temp_years == Int[]
        @test pp0.valid_years == pp.valid_years
        @test pp0.valid_climate_years == pp.valid_climate_years
    end

    @testset "Orchestrator contract: explicit PPT-qualified and temperature frames only" begin
        df = climate_multi_year(1991:2015)
        pp = preprocess_daily_data(df)
        cdata = pp.data[in.(pp.data.water_year, Ref(Set(pp.valid_climate_years))), :]
        temp_data = pp.data[in.(pp.data.water_year, Ref(Set(pp.valid_temp_years))), :]

        # No explicit climate_data -> NO precipitation keys (no legacy fallback to
        # gage_data: review finding 2); TMIN/TMAX inside gage_data alone -> NO temperature keys
        r_bare = calculate_all_signatures(pp.data, true)
        @test !haskey(r_bare, "prcp_annual_mean")
        @test !haskey(r_bare, "temp_max_annual_mean")

        # Explicit PPT-qualified frame -> precipitation keys
        r_no_temp = calculate_all_signatures(pp.data, true; climate_data=cdata)
        @test haskey(r_no_temp, "prcp_annual_mean")
        @test !haskey(r_no_temp, "temp_max_annual_mean")

        # Explicit temp_data -> temperature keys
        r_temp = calculate_all_signatures(pp.data, true; climate_data=cdata, temp_data=temp_data)
        @test haskey(r_temp, "temp_max_annual_mean")
        @test !isnan(r_temp["temp_max_annual_mean"])
        # Everything else is byte-identical with or without the temperature frame
        shared = [k for k in keys(r_no_temp)]
        @test all(isequal(r_no_temp[k], r_temp[k]) for k in shared)

        # 0-row temperature frame -> the full NaN temperature key set
        r_empty = calculate_all_signatures(pp.data, true; climate_data=cdata, temp_data=temp_data[1:0, :])
        @test haskey(r_empty, "temp_max_annual_mean") && isnan(r_empty["temp_max_annual_mean"])
        @test count(k -> startswith(k, "temp_") || startswith(k, "frost_") || startswith(k, "gdd_"), collect(keys(r_empty))) == 19 * 8

        # No PPT -> no precipitation keys; the un-normalized gate does NOT remove them
        r_q = calculate_all_signatures(select(pp.data, Not([:PPT, :TMIN, :TMAX])), false)
        @test !haskey(r_q, "prcp_annual_mean")
        r_unnorm = calculate_all_signatures(pp.data, true; climate_data=cdata, area_normalized=false)
        @test haskey(r_unnorm, "prcp_annual_mean")
        @test !haskey(r_unnorm, "annual_runoff_ratio_mean")

        # A frame whose precipitation FAILED qualification never reaches the family:
        # negative PPT every day -> valid_climate_years empty -> no precipitation keys
        # through the runner's contract (climate_data only when climate years exist)
        bad = copy(df); bad.PPT .= -1.0
        ppb = preprocess_daily_data(bad)
        @test isempty(ppb.valid_climate_years)
        r_bad = calculate_all_signatures(ppb.data, true)
        @test !haskey(r_bad, "prcp_annual_mean")

        # The per-family try/catch is not tripped: no warnings on the happy path
        @test_logs min_level=Logging.Warn calculate_all_signatures(pp.data, true; climate_data=cdata, temp_data=temp_data)
    end

    @testset "Crossing-tolerance boundary (decision 2026-10-12: 1e-12 of the total)" begin
        # two storms of 1 - d and 1 + d mm (total 2 mm): d = 1e-12 lies inside the allowance
        # (the first storm crosses, day 7); d = 4e-12 lies outside (the second, day 301)
        @test StreamflowSignatures.CLIMATE_CROSSING_REL_TOL == 1e-12
        wy = 2001
        for (d, expected) in ((1e-12, 7.0), (4e-12, 301.0))
            P = zeros(365); P[7] = 1.0 - d; P[301] = 1.0 + d
            _, c = collect_precip(climate_year_frame(wy; P=P))
            @test annual_lookup(c, "prcp_50_dowy", wy) == expected
        end
    end

    @testset "Round-2 regressions: crossing tolerance, non-finite inputs are missing, canonical fraction" begin
        wy = 2001
        # two 7-day blocks of 0.1 mm: sequential prefix 0.7 vs half-total 0.7000000000000001 —
        # the 1e-12 relative tolerance makes the mathematical tie land on day 7
        P = zeros(365); P[1:7] .= 0.1; P[301:307] .= 0.1
        _, c = collect_precip(climate_year_frame(wy; P=P))
        @test annual_lookup(c, "prcp_50_dowy", wy) == 7.0
        # constant 0.1 mm over a leap water year reaches 50 % on day 183
        _, c2 = collect_precip(climate_year_frame(2000; P=fill(0.1, 366)))
        @test annual_lookup(c2, "prcp_50_dowy", 2000) == 183.0
        # a strictly-below-half day is still not a crossing: 49 % on day 1, rest on day 2
        _, c3 = collect_precip(climate_year_frame(wy; P=[49.0, 51.0, zeros(363)...]))
        @test annual_lookup(c3, "prcp_50_dowy", wy) == 2.0
        # ±Inf precipitation is a missing day: interpolated and the year stays PPT-valid
        Pi = fill(2.0, 365); Pi[100] = Inf
        ppi = preprocess_daily_data(climate_year_frame(wy; P=Pi))
        @test ppi.valid_climate_years == [wy]
        @test ppi.data.PPT[100] == 2.0
        Pn = fill(2.0, 365); Pn[100] = -Inf
        @test preprocess_daily_data(climate_year_frame(wy; P=Pn)).valid_climate_years == [wy]
        Pm = fill(2.0, 365); Pm[100:103] .= Inf        # a 4-day run of non-finite values fails
        @test isempty(preprocess_daily_data(climate_year_frame(wy; P=Pm)).valid_climate_years)
        # Q: 31 non-finite days reject the year like 31 NaNs
        Qi = fill(1.0, 365); Qi[1:31] .= Inf
        @test isempty(preprocess_daily_data(climate_year_frame(wy; P=fill(1.0, 365), Q=Qi)).valid_years)
        # the helper never interpolates from an infinite neighbour
        @test isequal(StreamflowSignatures._interpolate_internal_gaps!([0.0, NaN, Inf], 3), [0.0, NaN, Inf])
        # negative policy values are rejected at the API too
        @test_throws ArgumentError preprocess_daily_data(climate_year_frame(wy; P=fill(1.0, 365)); max_gap_temp=-1)
        @test StreamflowSignatures.CFG_CLIMATE_CUM_FRACTION == 0.5
    end

    @testset "Review regressions (2026-10-12): equal blocks, truncated leap year, month column" begin
        wy = 2001
        # Finding 1: two equal 5 mm blocks separated by a dry interval -> the 50 %
        # crossing is the end of the first block, whatever the rounding of 0.1 sums
        P = zeros(365); P[1:50] .= 0.1; P[301:350] .= 0.1
        _, c = collect_precip(climate_year_frame(wy; P=P))
        @test annual_lookup(c, "prcp_50_dowy", wy) == 50.0
        # and the annual total equals the last cumulative value exactly
        @test annual_lookup(c, "prcp_annual", wy) == foldl(+, P)

        # Finding 3: WY 2000 has 366 days; 365 rows with dowy 1..365 is NOT complete
        df_trunc = climate_year_frame(2000; P=fill(2.0, 366), TN=fill(-1.0, 366), TX=fill(5.0, 366))[1:365, :]
        _, c2 = collect_precip(df_trunc)
        @test isnan(annual_lookup(c2, "prcp_annual", 2000))
        _, c3 = collect_temp(df_trunc)
        @test isnan(annual_lookup(c3, "gdd_cum_annual", 2000))

        # Finding 4: seasons come from the date, not from the month column
        df_m = climate_year_frame(wy; P=fill(2.0, 365))
        df_m.month = Vector{Union{Missing, Int}}(df_m.month)
        df_m.month[1] = missing
        _, c4 = collect_precip(df_m)
        parts = sum(annual_lookup(c4, "prcp_$(s)", wy) for s in ["win", "spr", "sum", "fal"])
        @test parts == 730.0
        @test annual_lookup(c4, "prcp_annual", wy) == 730.0
    end
end

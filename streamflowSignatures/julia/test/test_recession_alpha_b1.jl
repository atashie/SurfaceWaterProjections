# Tests for recession alpha under the fixed linear-reservoir assumption (b = 1).
# See docs/plans/recession_alpha_linear_reservoir_plan.md.
#
# Key invariants:
#   1. log_a_pointcloud / log_a_events / seasonality use b = 1:
#      log(a) = median(log(-dQ/dt) - log(Q)) — NOT the free-fit intercept.
#   2. b_pointcloud / b_events / concavity keep their FREE power-law fits.
#   3. For a true linear reservoir the two conventions coincide (continuity).

using StreamflowSignatures
using Test
using DataFrames
using Dates
using Statistics

# Q_{t+1} = alpha * Q_t (linear reservoir): -dQ = (1-alpha)*Q, so
# log_a = log(1 - alpha) exactly, and the free fit recovers b = 1.
function b1_linear_gage(n_years=25; alpha=0.95, peak=10.0)
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years
    dates = [start_date + Day(i-1) for i in 1:n_days]
    Q = Vector{Float64}(undef, n_days)
    q = peak
    for i in 1:n_days
        if day(dates[i]) == 1
            q = peak
        else
            q *= alpha
        end
        Q[i] = max(q, 1e-6)
    end
    df = DataFrame(gage_id=fill("TESTB1L", n_days), date=dates, Q=Q)
    return add_water_year_columns(df)
end

# Q_{t+1} = Q_t - k*Q_t^2 (quadratic reservoir): -dQ = k*Q^2 exactly, so the
# free fit recovers b = 2 while the b=1 convention gives
# log_a = median(log(k*Q_j)) — which VARIES along the recession and differs
# strongly from the free-fit intercept log(k). This is the discriminating case.
function b1_quadratic_gage(n_years=25; k=0.02, peak=10.0)
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years
    dates = [start_date + Day(i-1) for i in 1:n_days]
    Q = Vector{Float64}(undef, n_days)
    q = peak
    for i in 1:n_days
        if day(dates[i]) == 1
            q = peak
        else
            q = q - k * q^2
        end
        Q[i] = max(q, 1e-6)
    end
    df = DataFrame(gage_id=fill("TESTB1Q", n_days), date=dates, Q=Q)
    return add_water_year_columns(df)
end

# Seasonally-varying linear reservoir: alpha = 0.90 Oct-Mar, 0.98 Apr-Sep, so
# per-event log_a alternates between log(0.10) and log(0.02) keyed to the event's
# mid-event day-of-water-year. Used to pin the ALIGNMENT of the seasonality
# inputs (all_log_a vs all_dowy) — with a flat log_a signal, misalignment would
# be invisible to the sinusoid fit.
function b1_seasonal_gage(n_years=25; peak=10.0)
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years
    dates = [start_date + Day(i-1) for i in 1:n_days]
    Q = Vector{Float64}(undef, n_days)
    q = peak
    for i in 1:n_days
        alpha = month(dates[i]) in 4:9 ? 0.98 : 0.90
        if day(dates[i]) == 1
            q = peak
        else
            q *= alpha
        end
        Q[i] = max(q, 1e-9)
    end
    df = DataFrame(gage_id=fill("TESTB1S", n_days), date=dates, Q=Q)
    return add_water_year_columns(df)
end

# Flat baseline with exactly THREE isolated 0.9-exponential recessions across the
# whole record — far below the 25-event gate. Used to pin gate behavior and
# n_recession_events counting.
function b1_sparse_gage(n_years=25)
    start_date = Date(1990, 10, 1)
    n_days = 365 * n_years
    dates = [start_date + Day(i-1) for i in 1:n_days]
    Q = fill(5.0, n_days)
    for y_offset in (2, 10, 20)
        i0 = Dates.value(Date(1990 + y_offset, 11, 15) - start_date) + 1
        q = 10.0
        for j in 0:11
            Q[i0 + j] = q
            q *= 0.9
        end
    end
    df = DataFrame(gage_id=fill("TESTB1X", n_days), date=dates, Q=Q)
    return add_water_year_columns(df)
end

# Independent recompute of the seasonality INPUTS (per-event b=1 log_a paired
# with mid-event dowy), using the same event identification, fit-success gating,
# and pair conventions as the pipeline.
function recompute_event_log_a_dowy(df)
    la = Float64[]
    dw = Int[]
    for yr in unique(df.water_year)
        year_df = sort(df[df.water_year .== yr, :], :dowy)
        Q = Float64.(year_df.Q)
        dowy = Int.(year_df.dowy)
        for (s, e) in StreamflowSignatures.identify_recession_events(Q, dowy)
            Q_event = Q[s:e]
            log_a, b = StreamflowSignatures.fit_recession_power_law(Q_event; remove_first_day=true)
            (isnan(log_a) || isnan(b)) && continue
            push!(la, StreamflowSignatures.event_log_a_b1(Q_event))
            push!(dw, dowy[div(s + e, 2)])
        end
    end
    return (la, dw)
end

# Independent recompute of the per-year b=1 point-cloud log_a, using the same
# event identification and pair conventions as the pipeline.
function recompute_log_a_pc_b1(df)
    expected = Dict{Int, Float64}()
    for yr in unique(df.water_year)
        year_df = sort(df[df.water_year .== yr, :], :dowy)
        Q = Float64.(year_df.Q)
        dowy = Int.(year_df.dowy)
        vals = Float64[]
        for (s, e) in StreamflowSignatures.identify_recession_events(Q, dowy)
            Q_event = Q[s:e]
            # Only events whose free fit succeeds enter the point cloud
            log_a, b = StreamflowSignatures.fit_recession_power_law(Q_event; remove_first_day=true)
            (isnan(log_a) || isnan(b)) && continue
            Q_sub = Q_event[2:end]
            dQ = -diff(Q_sub)
            Q_mid = Q_sub[1:end-1]
            for j in eachindex(dQ)
                if Q_mid[j] > 0 && dQ[j] > 0
                    push!(vals, log(dQ[j]) - log(Q_mid[j]))
                end
            end
        end
        if length(vals) > 10
            expected[Int(yr)] = median(vals)
        end
    end
    return expected
end

@testset "Recession alpha with fixed b = 1" begin

    @testset "event_log_a_b1 helper" begin
        # Pure exponential event, alpha=0.9: every pair gives exactly log(0.1)
        Q_exp = 10.0 .* 0.9 .^ (0:19)
        @test isapprox(StreamflowSignatures.event_log_a_b1(collect(Q_exp)), log(0.1); atol=1e-12)

        # Too short -> NaN; non-recession (increasing) -> NaN
        @test isnan(StreamflowSignatures.event_log_a_b1([5.0, 4.0]))
        @test isnan(StreamflowSignatures.event_log_a_b1([1.0, 2.0, 3.0, 4.0]))
    end

    @testset "Linear reservoir: b=1 convention coincides with the free fit" begin
        df = b1_linear_gage(; alpha=0.95)
        r = analyze_recession_parameters(df)

        @test isapprox(r["log_a_pointcloud_mean"], log(0.05); atol=1e-9)
        @test isapprox(r["log_a_events_mean"], log(0.05); atol=1e-9)
        @test isapprox(r["b_pointcloud_mean"], 1.0; atol=1e-6)
        @test isapprox(r["b_events_mean"], 1.0; atol=1e-6)
        @test abs(r["concavity_mean"]) < 1e-9
        @test isapprox(r["recession_alpha_point_cloud_linear_reservoir"], 0.95; atol=1e-9)
        # log_a = log(1 - alpha) consistency with the discrete alpha
        @test isapprox(r["log_a_pointcloud_mean"],
                       log(1 - r["recession_alpha_point_cloud_linear_reservoir"]); atol=1e-9)
    end

    @testset "Quadratic reservoir: b stays free, alpha is decoupled from b" begin
        k = 0.02
        df = b1_quadratic_gage(; k=k)
        c = AnnualCollector()
        r = analyze_recession_parameters(df; collector=c)

        # FREE fits must still recover b = 2 (-dQ = k*Q^2 exactly by construction)
        @test isapprox(r["b_pointcloud_mean"], 2.0; atol=1e-6)
        @test isapprox(r["b_events_mean"], 2.0; atol=1e-6)

        # b=1 alpha: per-year values must equal the independent recompute exactly
        expected = recompute_log_a_pc_b1(df)
        @test !isempty(expected)
        pc_mask = c.signature .== "log_a_pointcloud"
        n_compared = 0
        for (i, wy) in enumerate(c.water_year)
            pc_mask[i] || continue
            if haskey(expected, Int(wy)) && !isnan(c.value[i])
                @test isapprox(c.value[i], expected[Int(wy)]; atol=1e-12)
                n_compared += 1
            end
        end
        @test n_compared >= 20

        # The OLD behavior (free/fitted-b intercept) would give ~log(k) here;
        # the b=1 value is median(log(k*Q_j)) with Q_j in (Q_min, peak) — a
        # strictly different number. Guard against regression to the old method.
        q_min = minimum(df.Q)
        @test log(k * q_min) < r["log_a_pointcloud_mean"] < log(k * 10.0)
        @test abs(r["log_a_pointcloud_mean"] - log(k)) > 0.3
        @test log(k * q_min) < r["log_a_events_mean"] < log(k * 10.0)
        @test abs(r["log_a_events_mean"] - log(k)) > 0.3

        # Seasonality now runs on the b=1 per-event values — still produced
        @test !isnan(r["log_a_seasonality_amplitude_all"])
        @test !isnan(r["log_a_seasonality_minimum_all"])
    end

    @testset "Seasonality-input alignment (all_log_a paired with all_dowy)" begin
        # Seasonally alternating alpha => strong seasonal log_a signal. The
        # pipeline's seasonality outputs must exactly reproduce a sinusoid fit
        # to independently recomputed (per-event b=1 log_a, mid-event dowy)
        # pairs — a misalignment between the two arrays would change the fit.
        df = b1_seasonal_gage()
        r = analyze_recession_parameters(df)

        la, dw = recompute_event_log_a_dowy(df)
        @test length(la) >= 250   # ~12 events/yr x 25 yr
        expected = StreamflowSignatures.fit_recession_seasonality(la, dw)

        @test isapprox(r["log_a_seasonality_amplitude_all"], expected["amplitude"]; rtol=1e-9)
        @test isapprox(r["log_a_seasonality_minimum_all"], expected["minimum_day"]; rtol=1e-9)

        # The signal must be non-trivial (otherwise alignment isn't exercised):
        # log_a alternates between log(0.10) and log(0.02), ~1.6 apart.
        @test r["log_a_seasonality_amplitude_all"] > 0.4
        # Minimum (most negative log_a = slow-recession season) in Apr-Sep
        # (dowy ~183-365 for an Oct 1 water year start)
        @test 200 < r["log_a_seasonality_minimum_all"] < 350
    end

    @testset "Event counts and min-events gate unchanged" begin
        # Monthly-spike gages: exactly one recession event per month, every year
        df = b1_linear_gage()
        c = AnnualCollector()
        r = analyze_recession_parameters(df; collector=c)
        @test isapprox(r["n_recession_events_mean"], 12.0; atol=1e-12)
        ev_mask = c.signature .== "n_recession_events"
        @test sum(ev_mask) == length(unique(df.water_year))
        @test all(c.value[ev_mask] .== 12.0)

        # Sparse gage: 3 events total < 25 => min-events gate NaNs all fitted
        # recession metrics, while n_recession_events and the alpha scalar
        # (both gate-independent) are still reported.
        dfs = b1_sparse_gage()
        rs = analyze_recession_parameters(dfs)
        @test isnan(rs["log_a_pointcloud_mean"])
        @test isnan(rs["log_a_events_mean"])
        @test isnan(rs["b_pointcloud_mean"])
        @test isnan(rs["b_events_mean"])
        @test isnan(rs["log_a_seasonality_amplitude_all"])
        @test isapprox(rs["n_recession_events_mean"], 3 / 25; atol=1e-12)
        # Alpha scalar: 3 events x 9 pairs = 27 > 10, every ratio exactly 0.9
        @test isapprox(rs["recession_alpha_point_cloud_linear_reservoir"], 0.9; atol=1e-9)
    end
end

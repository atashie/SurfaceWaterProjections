"""
Is the drought family redundant with the existing period-of-record low-pulse metrics?

    julia --project=julia docs/benchmarks/analyze_drought_redundancy.jl \\
        RUN_signatures_annual.parquet RUN_signatures.csv

Answers Codex MAJOR-3 (2026-07-27) properly. The earlier attempt divided
`mean(n_low_pulses_all)` × `mean(dur_low_pulses_all)` into the mean drought duration,
which is invalid: the product of two summary means equals mean annual sub-threshold days
only if annual count and duration are uncorrelated. This script works on the ANNUAL
series instead, where the identity IS exact — `dur_low_pulses_all` is the mean duration
of that year's low pulses, so `n × dur` is exactly that year's sub-threshold day count.

Per gage-year it compares:
  drought_duration_fixed_p10   (7-day smoothed Q, Weibull/type-6 10th percentile, strict <)
  pulse_days = n_low_pulses_all * dur_low_pulses_all   (raw Q, type-7 10th percentile)

and reports Pearson + Spearman correlation and the absolute-difference distribution,
overall and stratified by hydrologic regime (aridity proxy = annual runoff ratio
quartiles; intermittent gages = fixed p10 threshold of exactly 0).
"""

using Parquet2
using CSV
using DataFrames
using Statistics
using StatsBase
using Printf

const LEVELS = (2, 5, 10, 20, 30)
const DUR = "drought_duration_fixed_p10"   # the level with an existing analogue
const NPUL = "n_low_pulses_all"
const DPUL = "dur_low_pulses_all"

function pearson(x, y)
    length(x) < 3 && return NaN
    cor(x, y)
end

function spearman(x, y)
    length(x) < 3 && return NaN
    cor(tiedrank(x), tiedrank(y))
end

"""
Report a comparison. Reports the EXACT-ZERO FRACTION, nonzero count and max |diff|
alongside the quantiles: a median and p90 of 0 prove only that >=90% of pairs agree,
NOT that all do (Codex MAJOR-1, 2026-07-28 — an earlier version of this script printed
only the quantiles and the surrounding docs wrongly claimed "exactly zero").
"""
function report(label, x, y)
    n = length(x)
    if n < 3
        @printf("%-34s n=%-7d (too few)\n", label, n)
        return
    end
    d = x .- y
    ad = abs.(d)
    nz = count(>(0.0), ad)
    @printf("%-34s n=%-7d r=%6.3f rho=%6.3f  med|d|=%5.1f p90|d|=%6.1f MAX|d|=%6.1f  exact-equal=%6.2f%% (%d differ)  mean(dr-pu)=%+6.2f\n",
            label, n, pearson(x, y), spearman(x, y), median(ad), quantile(ad, 0.90),
            maximum(ad), 100 * (n - nz) / n, nz, mean(d))
end

"""
Within-gage (interannual) agreement — the quantity a TREND analysis actually consumes.
Comparing a median absolute difference against the series MEAN is the wrong scale
(Codex MAJOR-2): what matters is whether the two series move together year to year and
how the disagreement compares with the interannual SD.
"""
function within_gage(label, gages, x, y; min_years=10)
    cors = Float64[]
    nrmse = Float64[]
    for g in unique(gages)
        m = gages .== g
        count(m) >= min_years || continue
        xi = x[m]; yi = y[m]
        sx = std(xi)
        sx > 0 || continue                      # constant series: correlation undefined
        std(yi) > 0 || continue
        push!(cors, cor(xi, yi))
        push!(nrmse, sqrt(mean((xi .- yi) .^ 2)) / sx)
    end
    if isempty(cors)
        println("$label: no gage with variable series")
        return
    end
    @printf("%-34s gages=%-6d within-gage r: median %5.3f  p10 %5.3f  |  RMSE/SD: median %5.3f  p90 %5.3f\n",
            label, length(cors), median(cors), quantile(cors, 0.10),
            median(nrmse), quantile(nrmse, 0.90))
end

function main(argv)
    length(argv) == 2 || error("usage: analyze_drought_redundancy.jl ANNUAL.parquet SUMMARY.csv")
    ann_path, sum_path = argv

    println("=" ^ 118)
    println("DROUGHT vs LOW-PULSE REDUNDANCY — annual-series comparison")
    println("  annual : $ann_path")
    println("  summary: $sum_path")
    println("=" ^ 118)

    wanted = Set(vcat([NPUL, DPUL], ["drought_duration_fixed_p$(p)" for p in LEVELS]))
    ann = DataFrame(Parquet2.Dataset(ann_path); copycols=false)
    keep = ann[in.(ann.signature, Ref(wanted)), :]
    wide = unstack(keep, [:gage_id, :water_year], :signature, :value)
    for c in (DUR, NPUL, DPUL)
        c in names(wide) || error("annual parquet lacks signature $c")
    end

    # n × dur is exactly the year's sub-threshold day count; dur is NaN when n == 0
    pulse_days = map(eachrow(wide)) do r
        n = r[NPUL]; d = r[DPUL]
        (ismissing(n) || isnan(n)) && return NaN
        n == 0 && return 0.0
        (ismissing(d) || isnan(d)) && return NaN
        n * d
    end
    wide.pulse_days = pulse_days
    wide.drought_days = [ismissing(v) ? NaN : Float64(v) for v in wide[!, DUR]]
    ok = .!isnan.(wide.pulse_days) .& .!isnan.(wide.drought_days)
    w = wide[ok, :]
    @printf("\ngage-years compared: %d (%d gages)\n\n", nrow(w), length(unique(w.gage_id)))

    println("--- p10 (the only level with an existing analogue) ---")
    report("ALL gage-years", w.drought_days, w.pulse_days)
    within_gage("ALL gages (interannual)", w.gage_id, w.drought_days, w.pulse_days)

    println("\n--- every level vs the SAME p10 pulse pair (do the other levels differ?) ---")
    for p in LEVELS
        col = "drought_duration_fixed_p$(p)"
        col in names(wide) || continue
        v = [ismissing(x) ? NaN : Float64(x) for x in wide[!, col]]
        m = .!isnan.(v) .& .!isnan.(wide.pulse_days)
        report("  p$(p) vs pulse p10", v[m], wide.pulse_days[m])
    end
    println()

    # --- stratify by regime, using the summary CSV -------------------------------
    s = CSV.read(sum_path, DataFrame; types=Dict("gage_id" => String))
    thrcol = "drought_threshold_fixed_p10"
    rrcol = "annual_runoff_ratio_mean"

    if thrcol in names(s)
        interm = Set(s.gage_id[[!ismissing(v) && v isa Number && v == 0 for v in s[!, thrcol]]])
        if !isempty(interm)
            m = in.(w.gage_id, Ref(interm))
            report("  intermittent (thr p10 == 0)", w.drought_days[m], w.pulse_days[m])
            report("  perennial (thr p10 > 0)", w.drought_days[.!m], w.pulse_days[.!m])
        else
            println("  (no gage with a zero p10 threshold in this run)")
        end
    end

    if rrcol in names(s)
        rr = Dict{String, Float64}()
        for r in eachrow(s)
            v = r[rrcol]
            (ismissing(v) || !(v isa Number) || isnan(v)) && continue
            rr[r.gage_id] = Float64(v)
        end
        vals = collect(values(rr))
        if length(vals) >= 40
            qs = quantile(vals, [0.25, 0.5, 0.75])
            @printf("\naridity proxy = %s quartiles at %.3f / %.3f / %.3f (n=%d gages)\n",
                    rrcol, qs[1], qs[2], qs[3], length(vals))
            labels = ("Q1 driest (lowest runoff ratio)", "Q2", "Q3", "Q4 wettest")
            bins = [String[] for _ in 1:4]
            for (g, v) in rr
                b = v <= qs[1] ? 1 : v <= qs[2] ? 2 : v <= qs[3] ? 3 : 4
                push!(bins[b], g)
            end
            for b in 1:4
                m = in.(w.gage_id, Ref(Set(bins[b])))
                report("  " * labels[b], w.drought_days[m], w.pulse_days[m])
            end
        end
    end

    println("\nInterpretation: high correlation is EXPECTED (both are low-flow day counts).")
    println("The question is whether the drought family carries information the pulse pair")
    println("does not — judge that on the absolute differences, on regime-dependence of the")
    println("gap, and on `deficit`, which has no counterpart in the existing output at all.")
    return 0
end

exit(main(ARGS))

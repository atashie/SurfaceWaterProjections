# Tests for streamflow drought signatures (calculate_drought_metrics + helpers).
# Test map: docs/plans/2026-07-27-drought-signatures-plan.md §9.
#
# Conventions under test (plan §2): 7-day CENTERED smoothing applied within maximal
# runs of consecutive dates (never across a gap), fixed whole-record MAGNITUDE
# percentile thresholds via the unbiased Weibull plotting position i/(n+1)
# (Hyndman-Fan definition 6), strict `<` comparison, water-year aggregation, dense
# annual series whose zeros are valid values.
#
# The strongest assertions here are record invariants that need no external numbers:
# summed over water years, duration must equal the number of pooled smoothed values
# below the threshold (which is exactly floor((n+1)p) for a distinct-valued sample),
# and summed deficit must equal the pooled sum of departures. Together they pin the
# threshold definition AND the day -> water-year attribution.

using Test
using DataFrames
using Dates
using Statistics
using Logging
using StreamflowSignatures

const DROUGHT_CP = (start_year=1980, end_year=2024, min_total_obs=20, min_segment_obs=10)
const DROUGHT_LEVELS = [2, 5, 10, 20, 30]

# ---------------------------------------------------------------------------
# Builders / helpers (d-prefixed to avoid clashing with other test files that are
# included into the same top-level scope by runtests.jl)
# ---------------------------------------------------------------------------

"""Daily gage frame over [d0, d1] with `Q = q_fn(date, running_index)`."""
function dgage(d0::Date, d1::Date, q_fn)
    dates = collect(d0:Day(1):d1)
    Q = Float64[q_fn(dates[i], i) for i in eachindex(dates)]
    return add_water_year_columns(DataFrame(gage_id=fill("TESTD", length(dates)),
                                           date=dates, Q=Q))
end

# Seasonal + slow drift => every smoothed value distinct (exact-count invariants)
dseasonal(d, i) = 5.0 + 3.0 * sin(2π * (dayofyear(d) - 100) / 365.0) + 0.0002 * i

"""Water year of a date (Oct 1 - Sep 30), independent of add_water_year_columns."""
dwy(d::Date) = Dates.year(d) + (Dates.month(d) >= 10 ? 1 : 0)

dcoll_frame(coll) = DataFrame(signature=coll.signature, wy=coll.water_year, value=coll.value)

"""Single collected value for (signature, water_year); fails if not unique."""
function dcval(adf::DataFrame, sig::String, wy::Int)
    rows = adf[(adf.signature .== sig) .& (adf.wy .== wy), :value]
    @test length(rows) == 1
    return rows[1]
end

dur_name(p) = "drought_duration_fixed_p$(p)"
def_name(p) = "drought_deficit_fixed_p$(p)"
thr_name(p) = "drought_threshold_fixed_p$(p)"

# ---------------------------------------------------------------------------
@testset "Drought signatures" begin

@testset "config + column-name contract" begin
    # Pins the shipped config and the frozen output column names
    @test StreamflowSignatures.CFG_DROUGHT_ENABLED
    @test StreamflowSignatures.CFG_DROUGHT_PERCENTILES == DROUGHT_LEVELS
    @test StreamflowSignatures.CFG_DROUGHT_THRESHOLD_METHOD == "fixed"
    @test StreamflowSignatures.CFG_DROUGHT_PLOTTING_POSITION == "weibull"
    @test StreamflowSignatures.CFG_DROUGHT_SMOOTH_WINDOW == 7
    @test StreamflowSignatures.CFG_DROUGHT_SMOOTH_ALIGNMENT == "center"
    @test StreamflowSignatures.CFG_DROUGHT_BELOW_RANGE_POLICY == "na"

    @test length(StreamflowSignatures.DROUGHT_METRICS) == 2 * length(DROUGHT_LEVELS)
    for p in DROUGHT_LEVELS
        @test dur_name(p) in StreamflowSignatures.DROUGHT_METRICS
        @test def_name(p) in StreamflowSignatures.DROUGHT_METRICS
        @test thr_name(p) in StreamflowSignatures.DROUGHT_THRESHOLD_SCALARS
    end
end

@testset "weibull_quantile" begin
    x = Float64[1, 2, 3, 4, 5, 6, 7, 8, 9]   # n=9 => plotting positions i/10

    # Hand-derived: h = (n+1)p
    @test weibull_quantile(x, 0.10) ≈ 1.0      # h = 1.0 -> x[1]
    @test weibull_quantile(x, 0.30) ≈ 3.0      # h = 3.0 -> x[3]
    @test weibull_quantile(x, 0.25) ≈ 2.5      # h = 2.5 -> midpoint x[2],x[3]
    @test weibull_quantile(x, 0.90) ≈ 9.0      # h = 9.0 -> h >= n -> max

    # p below the smallest plotting position 1/(n+1) is NOT interpolable
    @test isnan(weibull_quantile(x, 0.05))
    @test isnan(weibull_quantile(x, 0.02))
    @test weibull_quantile(x, 0.05; below_range="clamp") ≈ 1.0

    # Agreement with Julia's built-in definition 6 (alpha = beta = 0), in range
    for p in (0.10, 0.20, 0.25, 0.30, 0.50, 0.75)
        @test weibull_quantile(x, p) ≈ quantile(x, p; alpha=0.0, beta=0.0)
    end

    # NaNs dropped; degenerate inputs
    @test weibull_quantile([NaN, 1.0, 2, 3, 4, 5, 6, 7, 8, 9, NaN], 0.10) ≈ 1.0
    @test isnan(weibull_quantile(Float64[], 0.10))
    @test isnan(weibull_quantile([NaN, NaN], 0.10))
    @test weibull_quantile([5.0], 0.5) ≈ 5.0
end

@testset "smooth_daily_flow" begin
    dates = collect(Date(2000, 1, 1):Day(1):Date(2000, 1, 10))
    Q = Float64.(1:10)

    # Centered, shrinking at run edges (min_valid=4 => day 1 averages days 1-4)
    @test smooth_daily_flow(dates, Q) ≈ [2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 7.5, 8.0, 8.5]

    # Trailing: the first 3 days hold < 4 values -> NaN
    st = smooth_daily_flow(dates, Q; alignment="trailing")
    @test all(isnan, st[1:3])
    @test st[4:10] ≈ [2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0]

    # A window must NEVER blend across a temporal gap: three runs at 10 / 20 / 30,
    # the last shorter than min_valid.
    d2 = vcat(collect(Date(2000, 1, 1):Day(1):Date(2000, 1, 5)),
              collect(Date(2000, 4, 10):Day(1):Date(2000, 4, 14)),
              collect(Date(2000, 7, 1):Day(1):Date(2000, 7, 3)))
    q2 = vcat(fill(10.0, 5), fill(20.0, 5), fill(30.0, 3))
    s2 = smooth_daily_flow(d2, q2)
    @test s2[1:5] ≈ fill(10.0, 5)      # no 20s leaked in
    @test s2[6:10] ≈ fill(20.0, 5)     # no 10s leaked in
    @test all(isnan, s2[11:13])         # run shorter than min_valid

    # A duplicated date also breaks a run (defensive: 99.0 must not contaminate)
    d3 = [Date(2000, 1, 1), Date(2000, 1, 1), Date(2000, 1, 2),
          Date(2000, 1, 3), Date(2000, 1, 4), Date(2000, 1, 5)]
    s3 = smooth_daily_flow(d3, Float64[99, 1, 2, 3, 4, 5])
    @test isnan(s3[1])
    @test s3[2:6] ≈ [2.5, 3.0, 3.0, 3.0, 3.5]

    # NaN days are skipped and do not count toward min_valid: with valid data only on
    # days 1,2,7,8,9,10 the four late days average to 8.5 and everything else is NaN.
    s4 = smooth_daily_flow(dates, Float64[1, 2, NaN, NaN, NaN, NaN, 7, 8, 9, 10])
    @test all(isnan, s4[1:6])
    @test s4[7:10] ≈ fill(8.5, 4)
    # Fewer than min_valid values anywhere in the run -> all NaN
    @test all(isnan, smooth_daily_flow(dates, Float64[1, NaN, NaN, NaN, NaN, NaN, NaN, 8, 9, 10]))

    # Edge cases
    @test isempty(smooth_daily_flow(Date[], Float64[]))
    @test_throws ArgumentError smooth_daily_flow(dates, Float64[1, 2])
end

@testset "record invariants (25-year distinct-valued gage)" begin
    df = dgage(Date(1990, 10, 1), Date(2015, 9, 30), dseasonal)
    n_years = length(unique(df.water_year))
    @test n_years == 25

    coll = AnnualCollector()
    res = calculate_drought_metrics(df; collector=coll)
    adf = dcoll_frame(coll)

    pool = filter(!isnan, smooth_daily_flow(df.date, Float64.(df.Q)))
    n = length(pool)
    @test n == nrow(df)                       # one contiguous run -> every day smoothed
    @test length(unique(pool)) == n           # distinct values (exact-count invariants)

    for p in DROUGHT_LEVELS
        thr = res[thr_name(p)]
        @test !isnan(thr)

        durs = adf[adf.signature .== dur_name(p), :value]
        defs = adf[adf.signature .== def_name(p), :value]
        @test length(durs) == n_years
        @test all(>=(0.0), durs)

        # Every sub-threshold day is attributed to exactly one water year...
        @test sum(durs) == count(<(thr), pool)
        # ...and for a distinct-valued sample that count is exactly floor((n+1)p),
        # which pins the Weibull threshold definition itself.
        @test sum(durs) == floor((n + 1) * p / 100)
        # Deficit is the pooled sum of departures, split across years without loss
        @test sum(defs) ≈ sum(thr - v for v in pool if v < thr) rtol=1e-9
    end

    # Documented self-check (plan §11.6): a 10% threshold sits below ~10% of days,
    # so the mean duration is ~36.5 days/yr. 913/25 = 36.52 exactly here.
    @test res["drought_duration_fixed_p10_mean"] ≈ 36.52 rtol=1e-6

    # Thresholds and metrics increase with the percentile level, every year
    thrs = [res[thr_name(p)] for p in DROUGHT_LEVELS]
    @test issorted(thrs)
    for yr in unique(df.water_year)
        @test issorted([dcval(adf, dur_name(p), yr) for p in DROUGHT_LEVELS])
        @test issorted([dcval(adf, def_name(p), yr) for p in DROUGHT_LEVELS])
    end

    # Dense series -> trend statistics are computed under the shipped gates
    r_gated = calculate_drought_metrics(df; trend_completeness=0.60, decade_completeness=0.80)
    @test !isnan(r_gated["drought_duration_fixed_p10_senn_slp"])
    @test !isnan(r_gated["drought_deficit_fixed_p10_mk_pval"])
end

@testset "water-year attribution (Codex MAJOR-1)" begin
    # The conservation invariant above proves no day is lost or double-counted, but it
    # CANNOT detect a bijective misattribution (shifting every water_year by one keeps
    # the totals identical). These two checks pin group MEMBERSHIP, not just totals.

    # (a) Hand-derived boundary fixture. Baseline Q = 10 everywhere except a 7-day
    # low-flow block 2000-09-28..2000-10-04 straddling the Sep 30 -> Oct 1 boundary.
    # Only 13 of ~7,305 days are affected, so every level's threshold is exactly 10.0
    # and a drought day is simply any smoothed value < 10. The 7-day centered window
    # spreads the block to 2000-09-25..2000-10-07: SIX days in WY 2000, SEVEN in WY 2001.
    lo_start, lo_end = Date(2000, 9, 28), Date(2000, 10, 4)
    dfb = dgage(Date(1990, 10, 1), Date(2010, 9, 30),
                (d, i) -> (lo_start <= d <= lo_end) ? 1.0 : 10.0)
    collb = AnnualCollector()
    resb = calculate_drought_metrics(dfb; collector=collb)
    adfb = dcoll_frame(collb)

    for p in DROUGHT_LEVELS
        @test resb[thr_name(p)] == 10.0
        @test dcval(adfb, dur_name(p), 2000) == 6.0
        @test dcval(adfb, dur_name(p), 2001) == 7.0
        # Deficits are exact sevenths that sum to whole numbers:
        # WY2000 = (9+18+27+36+45+54)/7 = 27, WY2001 = 9 + 189/7 = 36
        @test dcval(adfb, def_name(p), 2000) ≈ 27.0 atol=1e-10
        @test dcval(adfb, def_name(p), 2001) ≈ 36.0 atol=1e-10
        # Every other water year is drought-free
        for yr in unique(dfb.water_year)
            (yr == 2000 || yr == 2001) && continue
            @test dcval(adfb, dur_name(p), yr) == 0.0
            @test dcval(adfb, def_name(p), yr) == 0.0
        end
    end

    # (b) General case: recompute EVERY per-year duration independently from the
    # smoothed series and an independent water-year mapping (`dwy`, not the frame's
    # water_year column). A shift-by-one misattribution fails this immediately.
    df = dgage(Date(1990, 10, 1), Date(2015, 9, 30), dseasonal)
    coll = AnnualCollector()
    res = calculate_drought_metrics(df; collector=coll)
    adf = dcoll_frame(coll)
    qsm = smooth_daily_flow(df.date, Float64.(df.Q))
    wyv = dwy.(df.date)
    for p in DROUGHT_LEVELS
        thr = res[thr_name(p)]
        for yr in unique(df.water_year)
            expected = count(i -> wyv[i] == yr && !isnan(qsm[i]) && qsm[i] < thr,
                             eachindex(qsm))
            @test dcval(adf, dur_name(p), yr) == Float64(expected)
        end
    end
end

@testset "invalid dates are excluded, not smoothed (Codex MINOR-1)" begin
    # A row with a missing/unparseable date has no known position in time: it must not
    # be averaged with its neighbours, and it must not contribute a value of its own.
    dates = Vector{Union{Missing, Date}}(collect(Date(2000, 1, 1):Day(1):Date(2000, 1, 10)))
    dates[6] = missing
    s = smooth_daily_flow(dates, Float64[fill(10.0, 5); 1000.0; fill(10.0, 4)])
    @test isnan(s[6])                       # the invalid row itself
    @test all(==(10.0), s[1:5])             # 1000.0 must not leak backwards
    @test all(==(10.0), s[7:10])            # ...nor forwards
end

@testset "wet years are valid zeros, not NaN" begin
    # First 13 water years low flow, last 11 high flow: the pooled 10% threshold sits
    # inside the dry-year distribution, so wet years must report 0 (a real value).
    df = dgage(Date(1990, 10, 1), Date(2014, 9, 30),
               (d, i) -> (dwy(d) <= 2002 ? 1.0 : 50.0) + 0.001 * dayofyear(d))
    coll = AnnualCollector()
    res = calculate_drought_metrics(df; collector=coll)
    adf = dcoll_frame(coll)

    for p in DROUGHT_LEVELS
        @test dcval(adf, dur_name(p), 2010) == 0.0      # wet year: zero, not NaN
        @test dcval(adf, def_name(p), 2010) == 0.0
        @test !isnan(dcval(adf, dur_name(p), 1995))
    end
    @test dcval(adf, dur_name(10), 1995) > 0.0          # dry year has drought days
    @test !isnan(res["drought_duration_fixed_p10_mean"])
end

@testset "intermittent gage: zero threshold, strict <" begin
    # 300 zero-flow days per water year -> every level's threshold is exactly 0, and a
    # STRICT comparison must report 0 drought days (zero flow is not "below zero").
    df = dgage(Date(1990, 10, 1), Date(2010, 9, 30),
               (d, i) -> (Dates.value(d - Date(dwy(d) - 1, 10, 1)) + 1) <= 300 ? 0.0 : 5.0)
    coll = AnnualCollector()
    res = calculate_drought_metrics(df; collector=coll)
    adf = dcoll_frame(coll)

    for p in DROUGHT_LEVELS
        @test res[thr_name(p)] == 0.0
        @test all(==(0.0), adf[adf.signature .== dur_name(p), :value])
        @test all(==(0.0), adf[adf.signature .== def_name(p), :value])
        @test res[dur_name(p) * "_mean"] == 0.0
    end
end

@testset "threshold data floor (short record)" begin
    df9 = dgage(Date(1990, 10, 1), Date(1999, 9, 30), dseasonal)
    @test length(unique(df9.water_year)) == 9        # below min_years_for_threshold = 10
    res = calculate_drought_metrics(df9)
    for p in DROUGHT_LEVELS
        @test isnan(res[thr_name(p)])
        @test isnan(res[dur_name(p) * "_mean"])
        @test isnan(res[def_name(p) * "_median"])
    end

    df10 = dgage(Date(1990, 10, 1), Date(2000, 9, 30), dseasonal)
    @test length(unique(df10.water_year)) == 10      # exactly at the floor -> computed
    res10 = calculate_drought_metrics(df10)
    @test !isnan(res10[thr_name(10)])
    @test !isnan(res10[dur_name(10) * "_mean"])
end

@testset "stats floor interaction" begin
    df15 = dgage(Date(1995, 10, 1), Date(2010, 9, 30), dseasonal)
    @test length(unique(df15.water_year)) == 15

    floored = calculate_drought_metrics(df15; min_values_for_stats=20, changepoint=DROUGHT_CP)
    for suffix in ("_mean", "_median", "_senn_slp", "_pettitt_pval")
        @test isnan(floored[dur_name(10) * suffix])
        @test isnan(floored[def_name(10) * suffix])
    end
    # ...but the family is not exempt only because of the floor: without it, stats compute
    unfloored = calculate_drought_metrics(df15)
    @test !isnan(unfloored[dur_name(10) * "_mean"])
end

@testset "schema contract + collector invariance" begin
    df = dgage(Date(1990, 10, 1), Date(2015, 9, 30), dseasonal)
    nm = length(StreamflowSignatures.DROUGHT_METRICS)
    nsc = length(StreamflowSignatures.DROUGHT_THRESHOLD_SCALARS)

    base = calculate_drought_metrics(df)
    @test length(base) == nm * 8 + nsc
    with_cp = calculate_drought_metrics(df; changepoint=DROUGHT_CP)
    @test length(with_cp) == nm * 16 + nsc

    full_keys = sort(collect(keys(with_cp)))

    # 0-row input emits the IDENTICAL key set
    empty_df = DataFrame(date=Date[], Q=Float64[], water_year=Int[])
    @test sort(collect(keys(calculate_drought_metrics(empty_df; changepoint=DROUGHT_CP)))) == full_keys

    # Missing required column: warns, and still emits the IDENTICAL key set
    bad = select(df, [:date, :water_year])
    r_bad = with_logger(NullLogger()) do
        calculate_drought_metrics(bad; changepoint=DROUGHT_CP)
    end
    @test sort(collect(keys(r_bad))) == full_keys
    @test all(isnan, values(r_bad))

    # Collector is read-only with respect to the statistics
    coll = AnnualCollector()
    with_coll = calculate_drought_metrics(df; changepoint=DROUGHT_CP, collector=coll)
    @test sort(collect(keys(with_coll))) == full_keys
    for k in full_keys
        a, b = with_cp[k], with_coll[k]
        @test (isnan(a) && isnan(b)) || a == b
    end

    # One row per (metric, water year), no duplicate keys
    n_years = length(unique(df.water_year))
    @test length(coll.signature) == nm * n_years
    @test length(unique(zip(coll.signature, coll.water_year))) == length(coll.signature)
    @test all(!isnan, coll.value)

    # Happy path emits no warnings
    @test_logs min_level=Logging.Warn calculate_drought_metrics(df; changepoint=DROUGHT_CP,
                                                               trend_completeness=0.60,
                                                               decade_completeness=0.80,
                                                               min_values_for_stats=20)
end

@testset "orchestrator wiring" begin
    df = dgage(Date(1990, 10, 1), Date(2015, 9, 30), dseasonal)
    res = calculate_all_signatures(df; gage_id="TESTD", trend_completeness=0.60,
                                   decade_completeness=0.80, min_values_for_stats=20,
                                   changepoint=DROUGHT_CP)
    for p in DROUGHT_LEVELS
        @test haskey(res, dur_name(p) * "_mean")
        @test haskey(res, def_name(p) * "_senn_slp")
        @test haskey(res, dur_name(p) * "_pettitt_pval")
        @test haskey(res, thr_name(p))
    end
    @test !isnan(res["drought_duration_fixed_p10_mean"])

    # The drought family reaches the annual-values export through the shared collector
    coll = AnnualCollector()
    calculate_all_signatures(df; gage_id="TESTD", collector=coll)
    @test "drought_duration_fixed_p10" in unique(coll.signature)
    @test "drought_deficit_fixed_p2" in unique(coll.signature)
end

end # testset Drought signatures

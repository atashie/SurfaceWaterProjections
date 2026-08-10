"""
Additivity gate: prove that a benchmark run ADDED columns without changing any
pre-existing one.

    julia --project=julia docs/benchmarks/check_additivity.jl NEW.csv REFERENCE.csv \\
        [--allow-shift col1,col2] [--expect-added N] [--tol T] [--min-finite-frac F]

Checks (any failure => nonzero exit):
  1. No reference column is missing from the new run.
  2. Gage sets are identical.
  3. Every shared column is value-identical per gage (NaN == NaN, missing == missing);
     columns named by --allow-shift are reported but not failed.
  4. Added columns are POPULATED at or above `--min-finite-frac` (default 0.5), reported
     by statistic suffix. "At least one finite value" is far too weak a test: the
     orchestrator's per-signature try/catch can turn a computation failure into columns
     that are empty for all but a handful of gages.

Motivation: unit tests cannot establish production additivity (Codex MAJOR-4,
2026-07-27; docs/plans/2026-07-27-drought-signatures-plan.md §16). `flagged_for_high_na`
is the one column expected to move whenever the column count changes, since its
denominator counts all signature fields.
"""

using CSV
using DataFrames
using Printf: @printf, @sprintf
using Statistics

function parse_args(argv)
    argv2 = String[]
    allow_shift = String[]
    expect_added = nothing
    max_report = 20
    tol = 0.0
    min_finite_frac = 0.5
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--allow-shift"
            append!(allow_shift, split(argv[i + 1], ","))
            i += 2
        elseif a == "--expect-added"
            expect_added = parse(Int, argv[i + 1]); i += 2
        elseif a == "--tol"
            tol = parse(Float64, argv[i + 1]); i += 2
        elseif a == "--min-finite-frac"
            min_finite_frac = parse(Float64, argv[i + 1]); i += 2
        elseif a == "--max-report"
            max_report = parse(Int, argv[i + 1]); i += 2
        else
            push!(argv2, a); i += 1
        end
    end
    length(argv2) == 2 || error("usage: check_additivity.jl NEW.csv REFERENCE.csv [--allow-shift a,b] [--expect-added N] [--tol T] [--min-finite-frac F] [--max-report K]")
    return argv2[1], argv2[2], allow_shift, expect_added, max_report, tol, min_finite_frac
end

"""Element equality that treats NaN as equal to NaN and missing as equal to missing."""
function same_value(a, b)
    am = ismissing(a); bm = ismissing(b)
    (am || bm) && return am && bm
    if a isa Number && b isa Number
        (isnan(a) || isnan(b)) && return isnan(a) && isnan(b)
        return a == b
    end
    return string(a) == string(b)
end

function main(argv)
    new_path, ref_path, allow_shift, expect_added, max_report, tol, min_finite_frac = parse_args(argv)
    allow = Set(allow_shift)

    println("=" ^ 78)
    println("ADDITIVITY GATE")
    println("  new       : $new_path")
    println("  reference : $ref_path")
    isempty(allow) || println("  allow-shift: $(join(sort(collect(allow)), ", "))")
    if tol > 0
        println("  tolerance  : $tol (relative; differences at or below this are reported")
        println("               as NOISE, not failures — use when the two runs come from")
        println("               different machines or Julia versions)")
    end
    println("=" ^ 78)

    new = CSV.read(new_path, DataFrame; types=Dict("gage_id" => String))
    ref = CSV.read(ref_path, DataFrame; types=Dict("gage_id" => String))
    @printf("new: %d rows x %d cols | reference: %d rows x %d cols\n",
            nrow(new), ncol(new), nrow(ref), ncol(ref))

    failures = String[]

    # --- 1. columns ---------------------------------------------------------
    ncols = Set(names(new)); rcols = Set(names(ref))
    dropped = sort(collect(setdiff(rcols, ncols)))
    added = sort(collect(setdiff(ncols, rcols)))
    if isempty(dropped)
        println("\n[PASS] 1 no reference column dropped")
    else
        println("\n[FAIL] 1 dropped $(length(dropped)) reference column(s): " *
                join(first(dropped, max_report), ", "))
        push!(failures, "dropped columns")
    end
    @printf("       added %d column(s)%s\n", length(added),
            expect_added === nothing ? "" : " (expected $expect_added)")
    if expect_added !== nothing && length(added) != expect_added
        println("[FAIL] 1b added-column count $(length(added)) != expected $expect_added")
        push!(failures, "added-column count")
    end
    isempty(added) || println("       e.g. " * join(first(added, min(6, length(added))), ", "))

    # --- 2. gage sets -------------------------------------------------------
    ng = Set(new.gage_id); rg = Set(ref.gage_id)
    if ng == rg
        println("[PASS] 2 gage sets identical ($(length(ng)) gages)")
    else
        @printf("[FAIL] 2 gage sets differ: new-only %d, ref-only %d\n",
                length(setdiff(ng, rg)), length(setdiff(rg, ng)))
        push!(failures, "gage sets")
    end

    # --- 3. shared columns value-identical ----------------------------------
    shared = sort(collect(intersect(ncols, rcols)))
    # align on gage_id
    nidx = Dict(g => i for (i, g) in enumerate(new.gage_id))
    order = [get(nidx, g, 0) for g in ref.gage_id]
    if any(==(0), order)
        println("[SKIP] 3 cannot align rows (gage sets differ)")
        push!(failures, "row alignment")
    else
        mismatched = Tuple{String, Int, Float64}[]   # col, n_material, max_absdiff
        noise = Tuple{String, Int, Float64}[]        # differ, but within tolerance
        shifted = Tuple{String, Int}[]
        for col in shared
            col == "gage_id" && continue
            a = new[order, col]; b = ref[!, col]
            nm = 0; n_noise = 0; maxd = 0.0; maxd_noise = 0.0
            @inbounds for i in eachindex(b)
                same_value(a[i], b[i]) && continue
                numeric = a[i] isa Number && b[i] isa Number && !ismissing(a[i]) &&
                          !ismissing(b[i]) && !isnan(a[i]) && !isnan(b[i])
                d = numeric ? abs(Float64(a[i]) - Float64(b[i])) : Inf
                # relative tolerance, with an absolute floor so values near 0 are covered
                if tol > 0 && numeric && d <= tol * max(1.0, abs(Float64(b[i])))
                    n_noise += 1
                    maxd_noise = max(maxd_noise, d)
                else
                    nm += 1
                    numeric && (maxd = max(maxd, d))
                end
            end
            if nm > 0
                col in allow ? push!(shifted, (col, nm)) : push!(mismatched, (col, nm, maxd))
            elseif n_noise > 0
                push!(noise, (col, n_noise, maxd_noise))
            end
        end

        if isempty(mismatched)
            println("[PASS] 3 all $(length(shared) - 1) shared columns unchanged" *
                    (isempty(noise) ? " (bitwise)" : " beyond tolerance"))
        else
            println("[FAIL] 3 $(length(mismatched)) shared column(s) changed materially:")
            for (c, n, d) in first(sort(mismatched, by = x -> -x[2]), max_report)
                @printf("         %-52s %6d rows, max|diff| %.6g\n", c, n, d)
            end
            push!(failures, "shared-column values")
        end
        if !isempty(noise)
            worst = maximum(x -> x[3], noise)
            println("[INFO] 3c $(length(noise)) column(s) differ only within tolerance " *
                    "(max|diff| $(@sprintf("%.3g", worst))) — expected when comparing " *
                    "across machines or Julia versions")
        end
        if !isempty(shifted)
            println("[INFO] 3b whitelisted shifts (expected):")
            for (c, n) in shifted
                @printf("         %-52s %6d rows\n", c, n)
            end
        end
    end

    # --- 4. added columns populated -----------------------------------------
    # A "has at least one finite value" test is far too weak: a signature that failed on
    # 6,677 of 6,678 gages via the orchestrator's try/catch would still pass it
    # (Codex MAJOR-3, 2026-07-28). Require a minimum finite fraction, and report the
    # distribution BY STATISTIC SUFFIX so gating (trend/floor) is distinguishable from
    # breakage — trend suffixes legitimately sit lower than _mean/_median.
    if !isempty(added)
        println("[    ] 4 added-column population (min-finite-frac = $min_finite_frac):")
        fracs = Dict(c => count(x -> !ismissing(x) && x isa Number && !isnan(x), new[!, c]) /
                          nrow(new) for c in added)
        by_suffix = Dict{String, Vector{Float64}}()
        for c in added
            m = match(r"_(mean|median|senn_slp|linear_slp|spearman_rho|spearman_pval|mk_rho|mk_pval|pettitt_[a-z_]+)$", c)
            push!(get!(by_suffix, m === nothing ? "(scalar)" : m.captures[1], Float64[]), fracs[c])
        end
        for k in sort(collect(keys(by_suffix)))
            v = by_suffix[k]
            @printf("         %-22s n=%-4d min %.3f  median %.3f  max %.3f\n",
                    k, length(v), minimum(v), median(v), maximum(v))
        end
        allf = collect(values(fracs))
        @printf("         OVERALL                n=%-4d min %.3f  median %.3f  max %.3f\n",
                length(allf), minimum(allf), median(allf), maximum(allf))

        under = sort([c for c in added if fracs[c] < min_finite_frac], by = c -> fracs[c])
        if isempty(under)
            println("[PASS] 4 every added column is populated at or above the threshold")
        else
            println("[FAIL] 4 $(length(under)) added column(s) below min-finite-frac:")
            for c in first(under, max_report)
                @printf("         %-52s %.4f\n", c, fracs[c])
            end
            push!(failures, "under-populated added columns")
        end
    end

    println("\n" * "=" ^ 78)
    if isempty(failures)
        println("ADDITIVITY GATE: PASS")
        return 0
    else
        println("ADDITIVITY GATE: FAIL — " * join(unique(failures), "; "))
        return 1
    end
end

exit(main(ARGS))

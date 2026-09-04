# Recompute the 12 QA/QC flag columns of a delivered summary CSV in place, using the
# canonical `compute_qa_flags()` — needed after any post-hoc masking (e.g.
# apply_stats_floor_mask.py) so flags stay consistent with the masked statistics
# (flagged_for_high_na's NA fraction, range/order flags computed from _mean columns).
#
# Usage: julia refresh_qa_flags.jl <signatures.csv> [--dry-run]
# Prints per-flag flip counts vs the incoming file; --dry-run reports without rewriting.
# NOTE: CSV.write re-serializes every value — for a byte-preserving recompute of
# flagged_for_high_na alone use docs/benchmarks/recompute_high_na_flag.py and treat
# this script (dry-run) as the canonical-Julia cross-check.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "julia"))
using StreamflowSignatures
using CSV
using DataFrames

function main()
    csv_path = ARGS[1]
    dry_run = "--dry-run" in ARGS
    df = CSV.read(csv_path, DataFrame; types=Dict("gage_id" => String))
    flag_cols = [c for c in get_flag_columns() if c in names(df)]
    old_flags = Dict(c => copy(df[!, c]) for c in flag_cols)

    refreshed = compute_qa_flags(df)

    println("QA flag flips after refresh ($(basename(csv_path))):")
    total = 0
    for c in [x for x in get_flag_columns() if x in names(refreshed)]
        if haskey(old_flags, c)
            oldv = coalesce.(old_flags[c], false)
            newv = coalesce.(refreshed[!, c], false)
            n = count(oldv .!= newv)
            total += n
            n > 0 && println("  $(rpad(c, 36)) $(n) gages flipped")
        else
            println("  $(rpad(c, 36)) (new column)")
        end
    end
    println("Total flag flips: $total")

    if dry_run
        println("Dry run — $(csv_path) not rewritten")
    else
        CSV.write(csv_path, refreshed)
        println("Rewrote $(csv_path) ($(nrow(refreshed)) rows, $(ncol(refreshed)) cols)")
    end
end

main()

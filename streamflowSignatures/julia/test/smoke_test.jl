"""
Smoke Test - Runs signature extraction on a small subset of gages.

Usage:
    julia --project=. test/smoke_test.jl
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using StreamflowSignatures
using DataFrames
using Dates

# ── Configuration ──────────────────────────────────────────────────────────────
const STREAMFLOW_PATH = raw"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet"
const CLIMATE_PATH = raw"D:\processedOuts_feb2026\daymet_1980_2023.parquet"

# 10 hardcoded test gages (same as Python smoke test)
const TEST_GAGES = [
    "01011000", "01030500", "01057000", "01073000", "01094400",
    "01118300", "01137500", "01152500", "01169000", "01181000",
]

# Expected minimum signature columns (non-climate)
const MIN_EXPECTED_COLS = 400

# ── Main ───────────────────────────────────────────────────────────────────────
function main()
    println("=" ^ 70)
    println("JULIA SMOKE TEST - Streamflow Signatures")
    println("=" ^ 70)

    t_start = time()

    # Phase 1: Load streamflow data
    println("\nPhase 1: Loading streamflow data...")
    if !isfile(STREAMFLOW_PATH)
        println("ERROR: Streamflow file not found: $STREAMFLOW_PATH")
        println("Edit STREAMFLOW_PATH in this file to point to your data.")
        return 1
    end

    df = read_parquet(STREAMFLOW_PATH)  # auto-normalizes columns
    df.gage_id = string.(df.gage_id)
    test_set = Set(TEST_GAGES)
    df = df[in.(df.gage_id, Ref(test_set)), :]
    df = add_water_year_columns(df)
    println("  Loaded $(nrow(df)) rows for $(length(unique(df.gage_id))) gages")

    # Phase 2: Load and merge climate data (optional)
    has_climate = false
    local climate_grouped
    if isfile(CLIMATE_PATH)
        println("\nPhase 2: Loading climate data...")
        climate = read_parquet(CLIMATE_PATH)  # auto-normalizes columns
        climate.gage_id = string.(climate.gage_id)
        climate = climate[in.(climate.gage_id, Ref(test_set)), :]
        climate = add_water_year_columns(climate)
        climate_grouped = groupby(climate, :gage_id)
        has_climate = true
        println("  Loaded $(nrow(climate)) climate rows")
    else
        println("\nPhase 2: Skipping climate data (file not found)")
    end

    # Phase 3: Process each gage
    println("\nPhase 3: Processing $(length(TEST_GAGES)) gages...")
    results = Dict{String, Any}[]
    n_passed = 0
    n_skipped = 0

    for gage_id in TEST_GAGES
        gage_data = df[df.gage_id .== gage_id, :]
        if nrow(gage_data) == 0
            println("  [$gage_id] Not found in data, skipping")
            n_skipped += 1
            continue
        end

        qual_years, qualifies = filter_qualifying_years(gage_data)
        if !qualifies
            println("  [$gage_id] Only $(length(qual_years)) qualifying years, skipping")
            n_skipped += 1
            continue
        end

        qual_set = Set(qual_years)
        gage_data = gage_data[in.(gage_data.water_year, Ref(qual_set)), :]

        # Merge climate if available
        if has_climate
            try
                gage_climate = DataFrame(climate_grouped[(gage_id,)])[:, [:gage_id, :date, :PPT]]
                gage_data = leftjoin(gage_data, gage_climate, on=[:gage_id, :date])
            catch
                # Climate data not available for this gage
            end
        end

        sigs = calculate_all_signatures(gage_data, has_climate; gage_id=gage_id)
        sigs["gage_id"] = gage_id

        n_cols = length(sigs) - 1  # exclude gage_id
        n_nonnull = sum(1 for (k, v) in sigs if k != "gage_id" && !ismissing(v) && !(v isa Number && isnan(v)))
        println("  [$gage_id] $(length(qual_years)) years, $n_cols signature cols, $n_nonnull non-null")

        push!(results, sigs)
        n_passed += 1
    end

    # Phase 4: Validation
    println("\n" * "=" ^ 70)
    println("VALIDATION")
    println("=" ^ 70)

    all_ok = true

    if n_passed == 0
        println("FAIL: No gages processed")
        return 1
    end

    # Get all unique keys across results
    all_keys = Set{String}()
    for r in results
        union!(all_keys, keys(r))
    end
    n_sig_cols = length(all_keys) - 1  # exclude gage_id

    # Check 1: Minimum column count
    if n_sig_cols >= MIN_EXPECTED_COLS
        println("  [PASS] Column count: $n_sig_cols >= $MIN_EXPECTED_COLS")
    else
        println("  [FAIL] Column count: $n_sig_cols < $MIN_EXPECTED_COLS")
        all_ok = false
    end

    # Check 2: Key signatures present
    key_sigs = ["Qann_mean", "flashinessRB_mean", "BFI_Eckhardt_mean", "D50_day_mean"]
    for sig in key_sigs
        if sig in all_keys
            vals = [r[sig] for r in results if haskey(r, sig) && !ismissing(r[sig]) && !(r[sig] isa Number && isnan(r[sig]))]
            if !isempty(vals)
                println("  [PASS] $sig: $(length(vals))/$n_passed non-NA, range [$(round(minimum(vals), digits=4)), $(round(maximum(vals), digits=4))]")
            else
                println("  [WARN] $sig: all values NA")
            end
        else
            println("  [FAIL] $sig: missing from output")
            all_ok = false
        end
    end

    # Check 3: Range validation
    if "Qann_mean" in all_keys
        vals = [r["Qann_mean"] for r in results if haskey(r, "Qann_mean") && !ismissing(r["Qann_mean"]) && !(r["Qann_mean"] isa Number && isnan(r["Qann_mean"]))]
        if !isempty(vals) && all(v -> 0 <= v <= 100000, vals)
            println("  [PASS] Qann_mean range check (0-100000)")
        elseif !isempty(vals)
            println("  [FAIL] Qann_mean out of range")
            all_ok = false
        end
    end

    if "BFI_Eckhardt_mean" in all_keys
        vals = [r["BFI_Eckhardt_mean"] for r in results if haskey(r, "BFI_Eckhardt_mean") && !ismissing(r["BFI_Eckhardt_mean"]) && !(r["BFI_Eckhardt_mean"] isa Number && isnan(r["BFI_Eckhardt_mean"]))]
        if !isempty(vals) && all(v -> 0 <= v <= 1, vals)
            println("  [PASS] BFI_Eckhardt_mean range check (0-1)")
        elseif !isempty(vals)
            println("  [FAIL] BFI_Eckhardt_mean out of range")
            all_ok = false
        end
    end

    # Check 4: Climate signatures if available
    if has_climate
        climate_sigs = ["annual_runoff_ratio_mean", "elasticity_static", "avg_storage_mean"]
        for sig in climate_sigs
            if sig in all_keys
                vals = [r[sig] for r in results if haskey(r, sig) && !ismissing(r[sig]) && !(r[sig] isa Number && isnan(r[sig]))]
                if !isempty(vals)
                    println("  [PASS] Climate sig $sig present with data")
                else
                    println("  [WARN] Climate sig $sig all-NA")
                end
            else
                println("  [WARN] Climate sig $sig missing")
            end
        end
    end

    elapsed = time() - t_start

    println("\n" * "=" ^ 70)
    println("Gages: $n_passed passed, $n_skipped skipped")
    println("Total time: $(round(elapsed, digits=1))s")

    if all_ok
        println("\nSTATUS: SMOKE TEST PASSED")
        return 0
    else
        println("\nSTATUS: SMOKE TEST FAILED")
        return 1
    end
end

exit(main())

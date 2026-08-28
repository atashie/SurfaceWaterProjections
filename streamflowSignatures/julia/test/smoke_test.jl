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
# Paths default to the Windows data drive but are ENV-overridable, using the SAME
# variable names as docs/benchmarks/run_julia_benchmark.jl, so the smoke test can run
# wherever the feb2026 inputs are mounted (e.g. STREAMFLOW_DATA_PATH=/Volumes/...).
const STREAMFLOW_PATH = get(ENV, "STREAMFLOW_DATA_PATH",
    raw"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet")
const CLIMATE_PATH = get(ENV, "STREAMFLOW_CLIMATE_PATH",
    raw"D:\processedOuts_feb2026\daymet_1980_2023.parquet")

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
        println("Set the STREAMFLOW_DATA_PATH (and STREAMFLOW_CLIMATE_PATH) ENV variables to point to your data.")
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
    has_swe = false
    climate_cols = [:gage_id, :date, :PPT]
    local climate_grouped
    if isfile(CLIMATE_PATH)
        println("\nPhase 2: Loading climate data...")
        climate = read_parquet(CLIMATE_PATH)  # auto-normalizes columns (incl. swe -> SWE)
        climate.gage_id = string.(climate.gage_id)
        climate = climate[in.(climate.gage_id, Ref(test_set)), :]
        climate = add_water_year_columns(climate)
        climate_grouped = groupby(climate, :gage_id)
        has_climate = true
        has_swe = "SWE" in names(climate)
        if has_swe
            push!(climate_cols, :SWE)
        end
        println("  Loaded $(nrow(climate)) climate rows (SWE: $(has_swe ? "yes" : "no"))")
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
                gage_climate = DataFrame(climate_grouped[(gage_id,)])[:, climate_cols]
                gage_data = leftjoin(gage_data, gage_climate, on=[:gage_id, :date])
            catch
                # Climate data not available for this gage
            end
        end

        # Snow input: canonical path — preprocess to obtain valid_swe_years and pass
        # the filtered frame, exercising the same plumbing as the benchmark's
        # non-legacy path (rest of the harness stays legacy-style)
        snow_input = nothing
        if "SWE" in names(gage_data)
            pp_snow = preprocess_daily_data(gage_data)
            snow_input = pp_snow.data[in.(pp_snow.data.water_year, Ref(Set(pp_snow.valid_swe_years))), :]
        end

        sigs = calculate_all_signatures(gage_data, has_climate; gage_id=gage_id, snow_data=snow_input)
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

    # Check 3b: Snow signatures if SWE available (presence + sane ranges)
    if has_swe
        snow_bases = ["swe_max", "swe_max_dowy", "snow_cover_days", "snow_on_dowy",
                      "snow_off_dowy", "melt_season_days", "melt_rate", "ssm",
                      "swe_apr1", "melt_before_peak", "melt_before_peak_pct",
                      "melt_before_peak_to_max_swe", "melt_com_dowy", "swe_max_to_ppt"]
        n_missing_snow = count(b -> !("$(b)_mean" in all_keys), snow_bases)
        if n_missing_snow == 0
            println("  [PASS] Snow signatures: all 14 bases present")
        else
            println("  [FAIL] Snow signatures: $n_missing_snow of 14 bases missing from output")
            all_ok = false
        end
        snow_ranges = [("swe_max_mean", 0.0, 5000.0), ("ssm_mean", -1.0, 1.0),
                       ("snow_cover_days_mean", 0.0, 366.0),
                       ("swe_max_dowy_mean", 1.0, 366.0),
                       ("snow_off_dowy_mean", 1.0, 367.0)]
        for (sig, lo, hi) in snow_ranges
            if sig in all_keys
                vals = [r[sig] for r in results if haskey(r, sig) && !ismissing(r[sig]) && !(r[sig] isa Number && isnan(r[sig]))]
                if isempty(vals)
                    println("  [WARN] $sig: all values NA (no snowy gage in test set?)")
                elseif all(v -> lo <= v <= hi, vals)
                    println("  [PASS] $sig range check ($lo-$hi), $(length(vals)) non-NA")
                else
                    println("  [FAIL] $sig out of range [$lo, $hi]")
                    all_ok = false
                end
            end
        end
    end

    # Check 3c: Drought signatures (presence + ranges + level monotonicity)
    if CFG_DROUGHT_ENABLED
        levels = CFG_DROUGHT_PERCENTILES
        drought_bases = vcat(["drought_duration_fixed_p$(p)" for p in levels],
                             ["drought_deficit_fixed_p$(p)" for p in levels])
        n_missing_drought = count(b -> !("$(b)_mean" in all_keys), drought_bases)
        if n_missing_drought == 0
            println("  [PASS] Drought signatures: all $(length(drought_bases)) bases present")
        else
            println("  [FAIL] Drought signatures: $n_missing_drought of $(length(drought_bases)) bases missing")
            all_ok = false
        end

        finite_vals(sig) = [r[sig] for r in results
                            if haskey(r, sig) && !ismissing(r[sig]) &&
                               !(r[sig] isa Number && isnan(r[sig]))]

        # These gages all carry 45+ water years, so every drought value MUST be finite:
        # an all-NA column is a failure, not a warning (a NaN sweep would otherwise
        # sail through the range and monotonicity checks below).
        n_gages = length(results)
        for p in levels
            for (sig, lo, hi) in (("drought_duration_fixed_p$(p)_mean", 0.0, 366.0),
                                  ("drought_deficit_fixed_p$(p)_mean", 0.0, Inf),
                                  ("drought_threshold_fixed_p$(p)", 0.0, Inf))
                if !(sig in all_keys)
                    println("  [FAIL] $sig missing from output")
                    all_ok = false
                    continue
                end
                vals = finite_vals(sig)
                if length(vals) < n_gages
                    println("  [FAIL] $sig: only $(length(vals))/$n_gages finite " *
                            "(these gages have 45+ years — expected all finite)")
                    all_ok = false
                elseif all(v -> lo <= v <= hi, vals)
                    println("  [PASS] $sig range check, $(length(vals))/$n_gages finite")
                else
                    println("  [FAIL] $sig out of range [$lo, $hi]")
                    all_ok = false
                end
            end
        end

        # Duration, deficit, and the thresholds themselves must be non-decreasing in
        # the percentile level for every gage (structural invariant). An incomplete
        # sequence FAILS — silently skipping it would make the check vacuous.
        mono_ok = true
        for r in results
            for prefix in ("drought_duration_fixed_p", "drought_deficit_fixed_p",
                           "drought_threshold_fixed_p")
                suffix = prefix == "drought_threshold_fixed_p" ? "" : "_mean"
                seq = Float64[]
                for p in levels
                    k = "$(prefix)$(p)$(suffix)"
                    haskey(r, k) && r[k] isa Number && !isnan(r[k]) && push!(seq, Float64(r[k]))
                end
                if length(seq) != length(levels) || !issorted(seq)
                    mono_ok = false
                end
            end
        end
        if mono_ok
            println("  [PASS] Drought metrics non-decreasing across percentile levels")
        else
            println("  [FAIL] Drought metrics not monotone (or incomplete) across levels")
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

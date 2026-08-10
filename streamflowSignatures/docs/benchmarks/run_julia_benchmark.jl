"""
Julia Benchmark Runner for Streamflow Signatures

Runs full signature extraction on production data and captures timing.
Uses per-year quality filtering matching R's process_signatures_from_parquet().
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "julia"))

using StreamflowSignatures
using DataFrames
using CSV
using Parquet2
using Dates
using JSON
using Logging
using SHA

# Configuration
const STREAMFLOW_PATH = get(ENV, "STREAMFLOW_DATA_PATH",
    raw"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet")
const CLIMATE_PATH = get(ENV, "STREAMFLOW_CLIMATE_PATH",
    raw"D:\processedOuts_feb2026\daymet_1980_2023.parquet")
const METADATA_PATH = get(ENV, "STREAMFLOW_METADATA_PATH",
    raw"D:\processedOuts_feb2026\combined_watershed_metadata_09feb2026.csv")
const OUTPUT_DIR = get(ENV, "STREAMFLOW_OUTPUT_DIR", @__DIR__)
const OUTPUT_PREFIX = get(ENV, "STREAMFLOW_OUTPUT_PREFIX", "julia")

function main()
    # Read experiment controls from ENV at runtime (not from module constants)
    # to avoid Julia precompilation caching issues
    local start_water_year = let
        v = get(ENV, "STREAMFLOW_START_WATER_YEAR", "")
        v != "" ? parse(Int, v) : nothing
    end
    local min_qualifying_frac = let
        v = get(ENV, "STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION", "")
        v != "" ? parse(Float64, v) : nothing
    end
    local end_water_year = let
        v = get(ENV, "STREAMFLOW_END_WATER_YEAR", "")
        v != "" ? parse(Int, v) : nothing
    end

    println("=" ^ 70)
    println("JULIA BENCHMARK - Streamflow Signatures")
    println("=" ^ 70)
    println("Start time: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    # Changepoint analysis configuration
    cp_config = if CFG_CHANGEPOINT_ENABLED
        (
            start_year = CFG_CP_START_WATER_YEAR,
            end_year = CFG_CP_END_WATER_YEAR,
            min_total_obs = CFG_CP_MIN_TOTAL_OBS,
            min_segment_obs = CFG_CP_MIN_SEGMENT_OBS
        )
    else
        nothing
    end

    println("Config: MIN_NUM_YEARS=$(CFG_MIN_NUM_YEARS), MIN_FRAC_GOOD_DATA=$(CFG_MIN_FRAC_GOOD_DATA), MIN_Q_VALUE=$(CFG_MIN_Q_VALUE), MIN_DAYS_ABOVE=$(CFG_MIN_DAYS_ABOVE_THRESHOLD)")
    println("Changepoint: enabled=$(CFG_CHANGEPOINT_ENABLED), WY=$(CFG_CP_START_WATER_YEAR)-$(CFG_CP_END_WATER_YEAR), min_obs=$(CFG_CP_MIN_TOTAL_OBS), min_seg=$(CFG_CP_MIN_SEGMENT_OBS)")
    println("Annual values: save=$(CFG_SAVE_ANNUAL_VALUES)")
    println("Stats floor: min_values_for_stats=$(CFG_MIN_VALUES_FOR_STATS) (recession/elasticity exempt)")
    println("Experiment: OUTPUT_PREFIX=$(OUTPUT_PREFIX), START_WY=$(start_water_year), END_WY=$(end_water_year), MIN_QUAL_FRAC=$(min_qualifying_frac)")
    println("Output dir: $(OUTPUT_DIR)")
    println()

    timing = Dict{String, Any}(
        "language" => "Julia",
        "start_time" => Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
        "phases" => Dict{String, Float64}()
    )

    # Phase 1: Load streamflow data
    println("Phase 1: Loading streamflow data...")
    t0 = time()

    if !isfile(STREAMFLOW_PATH)
        println("ERROR: Streamflow file not found: $STREAMFLOW_PATH")
        return 1
    end

    streamflow = read_parquet(STREAMFLOW_PATH)  # auto-normalizes columns

    # Check if water_year columns already exist
    if !("water_year" in names(streamflow))
        streamflow = add_water_year_columns(streamflow)
    end

    t1 = time()
    timing["phases"]["load_streamflow"] = t1 - t0
    println("  Loaded $(nrow(streamflow)) rows in $(round(t1-t0, digits=2))s")
    println("  Unique gages: $(length(unique(streamflow.gage_id)))")

    # Phase 2: Load climate data
    println("\nPhase 2: Loading climate data...")
    t0 = time()

    has_climate = false
    has_swe = false
    climate_cols = [:gage_id, :date, :PPT]
    local climate
    if isfile(CLIMATE_PATH)
        climate = read_parquet(CLIMATE_PATH)  # auto-normalizes columns (incl. swe -> SWE)

        # Memory: keep only the columns this runner uses BEFORE the water-year copy —
        # the full 8-column ~98M-row frame OOMs on 16 GB machines (result-identical:
        # the dropped columns are never read by any signature path in this runner)
        climate_keep = intersect(["gage_id", "date", "PPT", "SWE"], names(climate))
        select!(climate, climate_keep)

        # Add water year columns if needed
        if !("water_year" in names(climate))
            climate = add_water_year_columns(climate)
        end

        has_climate = true
        has_swe = "SWE" in names(climate)
        if has_swe
            push!(climate_cols, :SWE)
        end
        t1 = time()
        timing["phases"]["load_climate"] = t1 - t0
        println("  Loaded $(nrow(climate)) rows in $(round(t1-t0, digits=2))s")
        println("  SWE column: $(has_swe ? "present (snow metrics enabled)" : "absent (snow metrics skipped)")")
    else
        println("  Climate file not found, skipping climate signatures")
        timing["phases"]["load_climate"] = 0.0
    end

    # Phase 2b: Area-normalization lookup (Q-to-PPT gate).
    # Gages with no drainage area (HYDAT gap — canals, dam outflows, channel
    # splits) carry Q in raw m³/s; Q-to-PPT signatures (runoff ratios,
    # elasticity, Q-P seasonality, storage) are skipped for them because Q and
    # PPT units don't match. Keys use the same leading-zero-stripped join id
    # as the Phase 5 metadata merge.
    # Accepts any AbstractString (CSV.jl yields InlineStrings, not String) and
    # always returns a concrete String for stable Dict keys
    normalize_join_id(id::AbstractString) = begin
        jid = all(isdigit, id) ? lstrip(id, '0') : id
        isempty(jid) ? "0" : String(jid)
    end
    # Only an explicit false-like value gates (Bool false, numeric 0, "FALSE"/"0"
    # strings — harmonized with the Python/rpkg runners); missing/absent defaults
    # to normalized (true)
    is_explicit_false(v) = !ismissing(v) && (
        (v isa Real && v == 0) ||
        (v isa AbstractString && uppercase(strip(v)) in ("FALSE", "0"))
    )

    area_norm_lookup = Dict{String, Bool}()
    if isfile(METADATA_PATH)
        meta_header = names(CSV.read(METADATA_PATH, DataFrame; limit=0))
        if "area_normalized" in meta_header && "gage_id" in meta_header
            meta_an = CSV.read(METADATA_PATH, DataFrame; select=["gage_id", "area_normalized"])
            for row in eachrow(meta_an)
                area_norm_lookup[normalize_join_id(string(row.gage_id))] = !is_explicit_false(row.area_normalized)
            end
            n_unnorm = count(!, values(area_norm_lookup))
            println("  Area-normalization lookup: $(length(area_norm_lookup)) gages, $n_unnorm not normalized (Q-to-PPT signatures will be skipped for these)")
        else
            println("  Metadata has no area_normalized column — assuming all gages area-normalized")
        end
    else
        println("  Metadata not found — assuming all gages area-normalized")
    end

    # Phase 3: Per-year quality filtering
    use_legacy = CFG_USE_LEGACY_FILTERING
    println("\nPhase 3: Per-year quality filtering (legacy=$use_legacy)...")
    t0 = time()

    # Pre-group data for O(1) lookups
    println("  Pre-grouping streamflow data...")
    grouped_streamflow = groupby(streamflow, :gage_id)

    grouped_climate = nothing
    if has_climate
        println("  Pre-grouping climate data...")
        grouped_climate = groupby(climate, :gage_id)
    end

    # Per-gage filtering — build list of (original_gage_id, qualifying_years) tuples
    qualifying_entries = Vector{Tuple{Any, Vector{Int}}}()
    total_gages = length(grouped_streamflow)

    # Cache preprocessing results for non-legacy path
    preprocess_cache = Dict{Any, NamedTuple}()

    if use_legacy
        for (key, gage_data_view) in pairs(grouped_streamflow)
            gage_df = DataFrame(gage_data_view)
            qual_years, qualifies = filter_qualifying_years(gage_df)
            if qualifies
                push!(qualifying_entries, (key.gage_id, qual_years))
            end
        end
    else
        for (i_g, (key, gage_data_view)) in enumerate(pairs(grouped_streamflow))
            if i_g % 500 == 0
                println("  Preprocessing gage $i_g/$total_gages...")
            end
            gage_df = DataFrame(gage_data_view)

            # Merge climate data BEFORE preprocessing so PPT (and SWE) are available
            if has_climate && grouped_climate !== nothing
                try
                    gage_climate = DataFrame(grouped_climate[(key.gage_id,)])[:, climate_cols]
                    gage_df = leftjoin(gage_df, gage_climate, on=[:gage_id, :date])
                catch
                    # Climate data not available for this gage
                end
            end

            result = preprocess_daily_data(gage_df)

            # Apply water-year window filter if configured (start and/or end)
            if start_water_year !== nothing || end_water_year !== nothing
                wy_lo = start_water_year === nothing ? typemin(Int) : start_water_year
                wy_hi = end_water_year === nothing ? typemax(Int) : end_water_year
                wy_mask_data = (result.data.water_year .>= wy_lo) .& (result.data.water_year .<= wy_hi)
                wy_mask_sf = (result.seasonal_flags.water_year .>= wy_lo) .& (result.seasonal_flags.water_year .<= wy_hi)
                wy_mask_diag = (result.diagnostics.water_year .>= wy_lo) .& (result.diagnostics.water_year .<= wy_hi)
                result = (
                    data = result.data[wy_mask_data, :],
                    valid_years = filter(yr -> wy_lo <= yr <= wy_hi, result.valid_years),
                    valid_climate_years = filter(yr -> wy_lo <= yr <= wy_hi, result.valid_climate_years),
                    valid_swe_years = filter(yr -> wy_lo <= yr <= wy_hi, result.valid_swe_years),
                    rejected_years = result.rejected_years,
                    seasonal_flags = result.seasonal_flags[wy_mask_sf, :],
                    diagnostics = result.diagnostics[wy_mask_diag, :],
                )
            end

            # Apply min qualifying data fraction filter if configured.
            # Denominator = per-gage possible years within the [start, end] window:
            # window start .. the gage's last in-window year (window-capped).
            if min_qualifying_frac !== nothing
                all_wy = unique(gage_df.water_year)
                start_yr = start_water_year !== nothing ? start_water_year : minimum(all_wy)
                end_yr = end_water_year !== nothing ? end_water_year : typemax(Int)
                wy_in_range = filter(yr -> start_yr <= yr <= end_yr, all_wy)
                if !isempty(wy_in_range)
                    max_wy = maximum(wy_in_range)
                    total_possible = max_wy - start_yr + 1
                    frac = length(result.valid_years) / total_possible
                    if frac < min_qualifying_frac
                        continue  # Skip this gage
                    end
                end
            end

            if length(result.valid_years) >= CFG_MIN_NUM_YEARS
                push!(qualifying_entries, (key.gage_id, result.valid_years))
                preprocess_cache[key.gage_id] = result
            end
        end
    end

    println("  Total gages: $total_gages")
    println("  Qualifying gages: $(length(qualifying_entries))")

    t1 = time()
    timing["phases"]["filter_gages"] = t1 - t0

    # Memory: the non-legacy Phase 4 consumes ONLY preprocess_cache — release the
    # raw frames + groupings (several GB) before per-gage processing. Verified:
    # grouped_streamflow / grouped_climate / streamflow / climate are referenced
    # past this point only in the legacy branch. (Added 2026-07-15 after the WY1980
    # run thrashed at 20.8 GB commit on the 16 GB machine.)
    if !use_legacy
        grouped_streamflow = nothing
        grouped_climate = nothing
        streamflow = DataFrame()
        if has_climate
            climate = DataFrame()
        end
        GC.gc()
    end

    # Phase 4: Process signatures
    n_gages = length(qualifying_entries)
    println("\nPhase 4: Processing $n_gages gages...")
    t0 = time()

    all_results = Vector{Dict{String, Any}}()

    # Annual-values accumulators (long format; one collector per gage, drained
    # here with the gage_id tag). See docs/plans/annual_values_export_plan.md.
    save_annual = CFG_SAVE_ANNUAL_VALUES
    annual_gage_id = String[]
    annual_signature = String[]
    annual_water_year = Int32[]
    annual_value = Float64[]

    for (i, (orig_gage_id, qual_years)) in enumerate(qualifying_entries)
        if i % 100 == 0 || i == 1
            elapsed = time() - t0
            rate = i / elapsed
            eta = (n_gages - i) / rate
            println("  [$i/$n_gages] Processing... ($(round(rate, digits=1)) gages/s, ETA: $(round(eta/60, digits=1)) min)")
        end

        gage_id_str = string(orig_gage_id)
        gage_collector = save_annual ? AnnualCollector() : nothing
        gage_area_normalized = get(area_norm_lookup, normalize_join_id(gage_id_str), true)

        if use_legacy
            # Legacy path: raw data + per-year filter
            gage_data = DataFrame(grouped_streamflow[(orig_gage_id,)])
            if !isempty(qual_years)
                qual_set = Set(qual_years)
                gage_data = gage_data[in.(gage_data.water_year, Ref(qual_set)), :]
            end

            # Merge climate data if available
            if has_climate && grouped_climate !== nothing
                try
                    gage_climate = DataFrame(grouped_climate[(orig_gage_id,)])[:, climate_cols]
                    gage_data = leftjoin(gage_data, gage_climate, on=[:gage_id, :date])
                catch
                end
            end

            # Legacy snow input: explicit and unfiltered (deprecated compat mode —
            # legacy has no valid_swe_years concept, mirroring its climate handling).
            # NOTE: the stats floor (min_values_for_stats), like trend_completeness,
            # is intentionally NOT applied on this deprecated legacy path.
            legacy_snow = "SWE" in names(gage_data) ? gage_data : nothing

            signatures = calculate_all_signatures(gage_data, has_climate; gage_id=gage_id_str, changepoint=cp_config, collector=gage_collector, area_normalized=gage_area_normalized, snow_data=legacy_snow)
        else
            # New path: use preprocessed data with seasonal flags + trend completeness
            pp = preprocess_cache[orig_gage_id]
            gage_data = pp.data
            gage_has_climate = has_climate && length(pp.valid_climate_years) > 0

            # Filter to valid_climate_years for climate signatures (Bug fix:
            # prevent Q-valid/PPT-invalid years from leaking into climate functions)
            climate_data = nothing
            if gage_has_climate
                climate_yr_set = Set(pp.valid_climate_years)
                climate_data = gage_data[in.(gage_data.water_year, Ref(climate_yr_set)), :]
            end

            # Snow data: explicit opt-in frame filtered to SWE-valid years. Possibly
            # 0 rows (gage has SWE but no valid SWE year -> full NaN snow key set);
            # nothing when the gage has no SWE column at all -> snow keys absent.
            # NO fallback to gage_data inside the orchestrator (plan, Codex finding 1).
            snow_data = nothing
            if "SWE" in names(gage_data)
                swe_yr_set = Set(pp.valid_swe_years)
                snow_data = gage_data[in.(gage_data.water_year, Ref(swe_yr_set)), :]
            end

            signatures = calculate_all_signatures(
                gage_data, gage_has_climate;
                gage_id=gage_id_str,
                seasonal_flags=pp.seasonal_flags,
                trend_completeness=CFG_NA_TREND_MIN_FRACTION,
                decade_completeness=CFG_NA_DECADE_MIN_FRACTION,
                min_values_for_stats=CFG_MIN_VALUES_FOR_STATS,
                climate_data=climate_data,
                snow_data=snow_data,
                snow_climate_years=pp.valid_climate_years,
                changepoint=cp_config,
                collector=gage_collector,
                area_normalized=gage_area_normalized
            )
        end

        # Drain per-gage annual values into the global long-format accumulators
        if gage_collector !== nothing
            n_rows = length(gage_collector.value)
            if n_rows > 0
                append!(annual_gage_id, fill(gage_id_str, n_rows))
                append!(annual_signature, gage_collector.signature)
                append!(annual_water_year, gage_collector.water_year)
                append!(annual_value, gage_collector.value)
            end
        end

        signatures["gage_id"] = gage_id_str
        signatures["start_water_year"] = isempty(qual_years) ? NaN : Float64(minimum(qual_years))
        signatures["end_water_year"] = isempty(qual_years) ? NaN : Float64(maximum(qual_years))
        signatures["num_water_years"] = length(qual_years)

        # Add per-gage ice-affected day total from diagnostics
        if !use_legacy && haskey(preprocess_cache, orig_gage_id)
            pp_diag = preprocess_cache[orig_gage_id].diagnostics
            signatures["ice_affected_days_total"] = Float64(sum(pp_diag.na_cause_ice))
        end

        # Memory: this gage's cached preprocess frames are no longer needed —
        # evict so the cache shrinks as Phase 4 progresses
        if !use_legacy
            delete!(preprocess_cache, orig_gage_id)
        end

        push!(all_results, signatures)
    end

    t1 = time()
    timing["phases"]["process_signatures"] = t1 - t0
    println("  Processed $(length(all_results)) gages in $(round(t1-t0, digits=2))s")

    # Zero-gage guard
    if isempty(all_results)
        println("WARNING: Zero gages qualified. Writing empty output.")
        if save_annual
            println("  Annual values parquet skipped (zero gages).")
        end
        CSV.write(joinpath(OUTPUT_DIR, "$(OUTPUT_PREFIX)_signatures.csv"),
                  DataFrame(gage_id=String[]))
        return 0
    end

    # Phase 4b: Write annual values parquet (before metadata merge, so a
    # metadata failure can't lose it). Long format, deterministically sorted.
    if save_annual
        println("\nPhase 4b: Writing annual values parquet...")
        t0 = time()
        annual_df = DataFrame(
            gage_id = annual_gage_id,
            signature = annual_signature,
            water_year = annual_water_year,
            value = annual_value
        )
        sort!(annual_df, [:gage_id, :signature, :water_year])
        annual_path = joinpath(OUTPUT_DIR, "$(OUTPUT_PREFIX)_signatures_annual.parquet")
        Parquet2.writefile(annual_path, annual_df; compression_codec=:zstd)
        t1 = time()
        timing["phases"]["write_annual_values"] = t1 - t0
        timing["n_annual_rows"] = nrow(annual_df)
        println("  Saved $(nrow(annual_df)) rows ($(length(unique(annual_df.signature))) signatures) to $annual_path")
        println("  Size: $(round(filesize(annual_path)/1e6, digits=1)) MB, took $(round(t1-t0, digits=2))s")
    end

    # Phase 5: Merge metadata and compute QA/QC flags
    println("\nPhase 5: Merging metadata and computing QA/QC flags...")
    t0 = time()

    # Convert to DataFrame
    # Get all unique keys
    all_keys = Set{String}()
    for result in all_results
        union!(all_keys, keys(result))
    end

    # Create DataFrame with all columns
    results_dict = Dict{String, Vector{Any}}()
    for key in all_keys
        results_dict[key] = [get(r, key, missing) for r in all_results]
    end

    results_df = DataFrame(results_dict)

    # Load and merge metadata
    if isfile(METADATA_PATH)
        println("  Loading metadata...")
        metadata = CSV.read(METADATA_PATH, DataFrame)

        # Ensure gage_id is String in both DataFrames for merging
        # Metadata stores USGS IDs without leading zeros (e.g. "1011000"),
        # while benchmark output preserves them (e.g. "01011000").
        # Strip leading zeros from NUMERIC-only IDs for matching, then restore originals.
        results_df.gage_id = string.(results_df.gage_id)
        metadata.gage_id = string.(metadata.gage_id)

        # Create normalized join key: strip leading zeros from numeric IDs only
        results_df[!, :_join_id] = [all(isdigit, id) ? lstrip(id, '0') : id for id in results_df.gage_id]
        metadata[!, :_join_id] = [all(isdigit, id) ? lstrip(id, '0') : id for id in metadata.gage_id]
        # Handle edge case: all-zeros ID
        results_df._join_id = [id == "" ? "0" : id for id in results_df._join_id]
        metadata._join_id = [id == "" ? "0" : id for id in metadata._join_id]

        # Rename basin_area_km2 to basin_area to match R output
        if "basin_area_km2" in names(metadata)
            rename!(metadata, :basin_area_km2 => :basin_area)
        end

        # Select basic metadata columns to merge
        metadata_cols = [
            "gage_id", "latitude", "longitude", "basin_area",
            "gage_type", "area_normalized"
        ]

        # Only keep columns that exist in metadata
        available_meta_cols = [c for c in metadata_cols if c in names(metadata)]

        if length(available_meta_cols) > 1  # More than just gage_id
            println("  Merging $(length(available_meta_cols)-1) metadata columns...")
            # Join on normalized key, then drop temp columns
            meta_join_cols = [c == "gage_id" ? "_join_id" : c for c in available_meta_cols]
            results_df = leftjoin(results_df, metadata[:, Symbol.(meta_join_cols)], on=:_join_id)
            n_matched = count(!ismissing, results_df.latitude)
            println("  Matched $n_matched / $(nrow(results_df)) gages with metadata")
            select!(results_df, Not(:_join_id))
        else
            println("  Warning: No metadata columns found to merge")
            select!(results_df, Not(:_join_id))
        end

        # Enrich with human interference metadata from GAGES-II
        println("  Loading GAGES-II interference metadata...")
        gages_ii = load_gages_ii_interference()
        if nrow(gages_ii) > 0 && "STAID" in names(gages_ii)
            gages_ii.STAID = string.(gages_ii.STAID)
            interference_cols = [
                "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
                "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
                "HYDRO_DISTURB_INDX", "CLASS",
            ]
            avail_int_cols = [c for c in interference_cols if c in names(gages_ii)]
            gages_subset = select(gages_ii, vcat(["STAID"], avail_int_cols))
            rename!(gages_subset, :STAID => :gage_id)
            results_df = leftjoin(results_df, gages_subset, on=:gage_id)
            n_matched = sum(.!ismissing.(results_df[!, Symbol(avail_int_cols[1])]))
            println("  Matched $n_matched gages with GAGES-II data")

            # Compute human_interference_class from CLASS
            results_df[!, :human_interference_class] = map(results_df[!, :CLASS]) do cls
                if ismissing(cls)
                    "unknown"
                elseif strip(string(cls)) == "Ref"
                    "reference"
                elseif strip(string(cls)) == "Non-ref"
                    "non-reference"
                else
                    "unknown"
                end
            end

            # Load Canadian HYDAT interference metadata (RHBN, REGULATED)
            println("  Loading Canadian HYDAT interference metadata...")
            canadian = load_canadian_interference()
            if nrow(canadian) > 0
                # Initialize RHBN/REGULATED as Union{Missing, Bool} to allow Bool assignment
                results_df[!, :RHBN] = Vector{Union{Missing, Bool}}(fill(missing, nrow(results_df)))
                results_df[!, :REGULATED] = Vector{Union{Missing, Bool}}(fill(missing, nrow(results_df)))
                # Build lookup from Canadian data
                can_lookup = Dict{String, NamedTuple{(:rhbn, :regulated, :hic), Tuple{Any, Any, String}}}()
                for row in eachrow(canadian)
                    can_lookup[string(row.gage_id)] = (
                        rhbn = hasproperty(row, :RHBN) ? row.RHBN : missing,
                        regulated = hasproperty(row, :REGULATED) ? row.REGULATED : missing,
                        hic = hasproperty(row, :human_interference_class) ? row.human_interference_class : "unknown"
                    )
                end
                n_can_matched = 0
                for i in 1:nrow(results_df)
                    gid = string(results_df.gage_id[i])
                    if haskey(can_lookup, gid)
                        results_df[i, :RHBN] = can_lookup[gid].rhbn
                        results_df[i, :REGULATED] = can_lookup[gid].regulated
                        # Override human_interference_class from Canadian data
                        results_df[i, :human_interference_class] = can_lookup[gid].hic
                        n_can_matched += 1
                    end
                end
                println("  Matched $n_can_matched Canadian gages with HYDAT metadata")
            else
                println("  Warning: No Canadian HYDAT metadata available - RHBN/REGULATED will be missing")
                results_df[!, :RHBN] = fill(missing, nrow(results_df))
                results_df[!, :REGULATED] = fill(missing, nrow(results_df))
            end
            println("  Interference columns added: $(length(avail_int_cols) + 3)")
        else
            println("  Warning: GAGES-II data not available")
        end
    else
        println("  Warning: Metadata file not found: $METADATA_PATH")
    end

    # Compute QA/QC flags
    println("  Computing QA/QC flags...")
    results_df = compute_qa_flags(results_df)

    # Organize columns: gage_id first, then metadata, then signatures, then flags
    metadata_order = [
        "gage_id", "latitude", "longitude", "basin_area",
        "gage_type", "num_water_years", "start_water_year", "end_water_year",
        "area_normalized",
        "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
        "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
        "HYDRO_DISTURB_INDX", "CLASS", "RHBN", "REGULATED",
        "human_interference_class",
    ]
    flag_cols = [c for c in names(results_df) if startswith(c, "flagged_")]
    signature_cols = [c for c in names(results_df)
                      if !(c in metadata_order) && !(c in flag_cols)]

    # Build final column order
    final_cols = Symbol[]
    for c in metadata_order
        if c in names(results_df)
            push!(final_cols, Symbol(c))
        end
    end
    append!(final_cols, Symbol.(sort(signature_cols)))
    append!(final_cols, Symbol.(sort(flag_cols)))

    results_df = results_df[:, final_cols]

    output_path = joinpath(OUTPUT_DIR, "$(OUTPUT_PREFIX)_signatures.csv")
    CSV.write(output_path, results_df)

    t1 = time()
    timing["phases"]["metadata_qaqc_save"] = t1 - t0
    println("  Saved to $output_path")
    println("  Shape: $(size(results_df))")

    # Calculate totals
    timing["end_time"] = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    timing["total_seconds"] = sum(values(timing["phases"]))
    timing["n_gages_processed"] = length(all_results)

    # Count columns: exclude gage_id, metadata cols, and flag cols
    n_meta = length([c for c in metadata_order if c in names(results_df)])
    n_flags = length([c for c in names(results_df) if startswith(c, "flagged_")])
    timing["n_signature_columns"] = ncol(results_df) - n_meta - n_flags
    timing["n_metadata_columns"] = n_meta
    timing["n_qaqc_flags"] = n_flags

    # --- Provenance (added 2026-07-28, Codex MAJOR-4) --------------------------------
    # Counts and timings alone cannot reconstruct a run: record the resolved inputs, the
    # config actually loaded (SHA-256 — it is small), the source revision, the Julia
    # version, and the experiment ENV overrides. Multi-GB inputs record size+mtime rather
    # than a hash by default; set STREAMFLOW_HASH_INPUTS=1 to hash them too (slow).
    timing["provenance"] = let
        hash_big = get(ENV, "STREAMFLOW_HASH_INPUTS", "") == "1"
        function file_info(path)
            isfile(path) || return Dict("path" => path, "exists" => false)
            info = Dict{String, Any}("path" => abspath(path), "exists" => true,
                                     "bytes" => filesize(path),
                                     "mtime" => Dates.format(Dates.unix2datetime(mtime(path)),
                                                             "yyyy-mm-ddTHH:MM:SS"))
            if hash_big || filesize(path) < 50_000_000
                info["sha256"] = bytes2hex(SHA.sha256(read(path)))
            end
            return info
        end
        config_file = get(ENV, "STREAMFLOW_CONFIG",
                          joinpath(@__DIR__, "..", "..", "config", "signatures_config.json"))
        git_rev = try
            strip(read(`git -C $(@__DIR__) rev-parse HEAD`, String))
        catch
            "unavailable"
        end
        git_dirty = try
            !isempty(strip(read(`git -C $(@__DIR__) status --porcelain`, String)))
        catch
            nothing
        end
        Dict{String, Any}(
            "streamflow" => file_info(STREAMFLOW_PATH),
            "climate" => file_info(CLIMATE_PATH),
            "metadata" => file_info(METADATA_PATH),
            "config" => file_info(config_file),
            "git_revision" => git_rev,
            "git_working_tree_dirty" => git_dirty,
            "julia_version" => string(VERSION),
            "hostname" => gethostname(),
            "env_overrides" => Dict(k => get(ENV, k, nothing) for k in
                ["STREAMFLOW_START_WATER_YEAR", "STREAMFLOW_END_WATER_YEAR",
                 "STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION", "STREAMFLOW_CONFIG",
                 "STREAMFLOW_GAGES_II_DIR", "STREAMFLOW_OUTPUT_PREFIX"] if haskey(ENV, k)),
        )
    end

    # Save timing
    timing_path = joinpath(OUTPUT_DIR, "$(OUTPUT_PREFIX)_timing.json")
    open(timing_path, "w") do f
        JSON.print(f, timing, 2)
    end
    println("  Timing saved to $timing_path")

    # Summary
    println("\n" * "=" ^ 70)
    println("BENCHMARK COMPLETE")
    println("=" ^ 70)
    println("Total time: $(round(timing["total_seconds"], digits=2))s ($(round(timing["total_seconds"]/60, digits=2)) min)")
    println("Gages processed: $(timing["n_gages_processed"])")
    println("Total columns: $(ncol(results_df))")
    println("  Signature columns: $(timing["n_signature_columns"])")
    println("  Metadata columns: $(timing["n_metadata_columns"])")
    println("  QA/QC flag columns: $(timing["n_qaqc_flags"])")
    println("Rate: $(round(timing["n_gages_processed"]/timing["total_seconds"], digits=2)) gages/s")

    return 0
end

# Run
exit(main())

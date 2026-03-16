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

# Configuration
const STREAMFLOW_PATH = raw"D:\processedOuts_feb2026\combined_streamflow_data_09feb2026.parquet"
const CLIMATE_PATH = raw"D:\processedOuts_feb2026\daymet_1980_2023.parquet"
const METADATA_PATH = raw"D:\processedOuts_feb2026\combined_watershed_metadata_09feb2026.csv"
const OUTPUT_DIR = @__DIR__



function main()
    println("=" ^ 70)
    println("JULIA BENCHMARK - Streamflow Signatures")
    println("=" ^ 70)
    println("Start time: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println("Config: MIN_NUM_YEARS=$(CFG_MIN_NUM_YEARS), MIN_FRAC_GOOD_DATA=$(CFG_MIN_FRAC_GOOD_DATA), MIN_Q_VALUE=$(CFG_MIN_Q_VALUE), MIN_DAYS_ABOVE=$(CFG_MIN_DAYS_ABOVE_THRESHOLD)")
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
    local climate
    if isfile(CLIMATE_PATH)
        climate = read_parquet(CLIMATE_PATH)  # auto-normalizes columns

        # Add water year columns if needed
        if !("water_year" in names(climate))
            climate = add_water_year_columns(climate)
        end

        has_climate = true
        t1 = time()
        timing["phases"]["load_climate"] = t1 - t0
        println("  Loaded $(nrow(climate)) rows in $(round(t1-t0, digits=2))s")
    else
        println("  Climate file not found, skipping climate signatures")
        timing["phases"]["load_climate"] = 0.0
    end

    # Phase 3: Per-year quality filtering (matching R's three-stage filter)
    println("\nPhase 3: Per-year quality filtering...")
    t0 = time()

    # Pre-group data for O(1) lookups
    println("  Pre-grouping streamflow data...")
    grouped_streamflow = groupby(streamflow, :gage_id)

    grouped_climate = nothing
    if has_climate
        println("  Pre-grouping climate data...")
        grouped_climate = groupby(climate, :gage_id)
    end

    # Per-gage, per-year filtering — build list of (original_gage_id, qualifying_years) tuples
    qualifying_entries = Vector{Tuple{Any, Vector{Int}}}()  # (original gage_id value, qualifying years)
    total_gages = length(grouped_streamflow)

    for (key, gage_data_view) in pairs(grouped_streamflow)
        gage_df = DataFrame(gage_data_view)
        qual_years, qualifies = filter_qualifying_years(gage_df)
        if qualifies
            push!(qualifying_entries, (key.gage_id, qual_years))
        end
    end

    println("  Total gages: $total_gages")
    println("  Qualifying gages: $(length(qualifying_entries))")

    t1 = time()
    timing["phases"]["filter_gages"] = t1 - t0

    # Phase 4: Process signatures
    n_gages = length(qualifying_entries)
    println("\nPhase 4: Processing $n_gages gages...")
    t0 = time()

    all_results = Vector{Dict{String, Any}}()

    for (i, (orig_gage_id, qual_years)) in enumerate(qualifying_entries)
        if i % 100 == 0 || i == 1
            elapsed = time() - t0
            rate = i / elapsed
            eta = (n_gages - i) / rate
            println("  [$i/$n_gages] Processing... ($(round(rate, digits=1)) gages/s, ETA: $(round(eta/60, digits=1)) min)")
        end

        gage_id_str = string(orig_gage_id)

        # O(1) lookup from pre-grouped data using original key
        gage_data = DataFrame(grouped_streamflow[(orig_gage_id,)])

        # Filter to qualifying water years only (matching R's per-year filter)
        if !isempty(qual_years)
            qual_set = Set(qual_years)
            gage_data = gage_data[in.(gage_data.water_year, Ref(qual_set)), :]
        end

        # Merge climate data if available — O(1) lookup
        if has_climate && grouped_climate !== nothing
            try
                gage_climate = DataFrame(grouped_climate[(orig_gage_id,)])[:, [:gage_id, :date, :PPT]]
                gage_data = leftjoin(gage_data, gage_climate, on=[:gage_id, :date])
            catch
                # Climate data not available for this gage
            end
        end

        # Calculate signatures
        signatures = calculate_all_signatures(gage_data, has_climate; gage_id=gage_id_str)
        signatures["gage_id"] = gage_id_str

        # Add computed metadata columns (matching R output)
        signatures["start_water_year"] = isempty(qual_years) ? NaN : Float64(minimum(qual_years))
        signatures["end_water_year"] = isempty(qual_years) ? NaN : Float64(maximum(qual_years))
        signatures["num_water_years"] = length(qual_years)

        push!(all_results, signatures)
    end

    t1 = time()
    timing["phases"]["process_signatures"] = t1 - t0
    println("  Processed $(length(all_results)) gages in $(round(t1-t0, digits=2))s")

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
        results_df.gage_id = string.(results_df.gage_id)
        metadata.gage_id = string.(metadata.gage_id)

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
            results_df = leftjoin(results_df, metadata[:, Symbol.(available_meta_cols)], on=:gage_id)
        else
            println("  Warning: No metadata columns found to merge")
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

            # Add empty RHBN/REGULATED columns (Canadian HYDAT not available from Julia)
            results_df[!, :RHBN] = fill(missing, nrow(results_df))
            results_df[!, :REGULATED] = fill(missing, nrow(results_df))
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

    output_path = joinpath(OUTPUT_DIR, "julia_signatures.csv")
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

    # Save timing
    timing_path = joinpath(OUTPUT_DIR, "julia_timing.json")
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

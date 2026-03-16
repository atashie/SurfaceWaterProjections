"""
Human interference metadata module.

Loads GAGES-II and HYDAT metadata for watershed characterization.
"""

using CSV
using DataFrames

# Configuration loaded from config.jl:
# CFG_INCLUDE_HUMAN_INTERFERENCE, CFG_GAGES_II_DIR, CFG_HYDAT_PATH, CFG_INTERFERENCE_COLUMNS


"""
    load_gages_ii_interference(gages_dir::String=CFG_GAGES_II_DIR) -> DataFrame

Load GAGES-II human interference metadata.

Parameters
----------
gages_dir : String
    Path to GAGES-II metadata directory

Returns
-------
DataFrame
    GAGES-II interference columns for USGS gages
"""
function load_gages_ii_interference(gages_dir::String=CFG_GAGES_II_DIR)
    if isnothing(gages_dir) || gages_dir == ""
        @warn "GAGES-II directory not configured"
        return DataFrame(
            STAID = String[],
            NDAMS_2009 = Float64[],
            MAJ_DDENS_2009 = Float64[],
            STOR_NID_2009 = Float64[],
            IMPNLCD06 = Float64[],
            DEVNLCD06 = Float64[],
            FRESHW_WITHDRAWAL = Float64[],
            HYDRO_DISTURB_INDX = Float64[],
            CLASS = String[]
        )
    end

    if !isdir(gages_dir)
        @warn "GAGES-II directory not found: $gages_dir"
        return DataFrame()
    end

    # CONUS files (matching R's config.R)
    conus_files = Dict(
        "hydromod_dams" => "conterm_hydromod_dams.txt",
        "pop_infrastr" => "conterm_pop_infrastr.txt",
        "hydromod_other" => "conterm_hydromod_other.txt",
        "bas_classif" => "conterm_bas_classif.txt",
        "lc06_basin" => "conterm_lc06_basin.txt",
    )

    # AKHIPR files (Alaska, Hawaii, Puerto Rico)
    akhipr_files = Dict(
        "hydromod_dams" => "AKHIPR_hydromod_dams.txt",
        "pop_infrastr" => "AKHIPR_pop_infrastr.txt",
        "hydromod_other" => "AKHIPR_hydromod_other.txt",
        "bas_classif" => "AKHIPR_bas_classif.txt",
        # Note: AKHIPR does not have lc06_basin file
    )

    combined_data = DataFrame()

    # Load CONUS data
    conus_data = load_gages_ii_region(gages_dir, conus_files)
    if nrow(conus_data) > 0
        combined_data = conus_data
    end

    # Load AKHIPR data
    akhipr_data = load_gages_ii_region(gages_dir, akhipr_files)
    if nrow(akhipr_data) > 0
        if nrow(combined_data) > 0
            combined_data = vcat(combined_data, akhipr_data; cols=:union)
        else
            combined_data = akhipr_data
        end
    end

    # Replace -999 (missing value sentinel) with missing
    for col in names(combined_data)
        if eltype(combined_data[!, col]) <: Union{Missing, Number}
            combined_data[!, col] = replace(combined_data[!, col], -999.0 => missing)
        end
    end

    return combined_data
end


"""
    load_gages_ii_region(gages_dir::String, files::Dict) -> DataFrame

Load GAGES-II data for a specific region (CONUS or AKHIPR).
"""
function load_gages_ii_region(gages_dir::String, files::Dict)
    result = DataFrame()

    # Define columns to extract from each file type
    # Dams file: NDAMS_2009, MAJ_DDENS_2009, STOR_NID_2009
    # Class file: CLASS
    # Landuse file: IMPNLCD06, DEVNLCD06, FRESHW_WITHDRAWAL, HYDRO_DISTURB_INDX

    for (file_type, filename) in files
        filepath = joinpath(gages_dir, filename)
        if !isfile(filepath)
            continue
        end

        try
            df = CSV.read(filepath, DataFrame; missingstring=["-999", "NA", ""],
                          types=Dict("STAID" => String))

            # Ensure STAID is string
            if "STAID" in names(df)
                df.STAID = string.(df.STAID)
            end

            if nrow(result) == 0
                result = df
            else
                # Merge on STAID
                result = leftjoin(result, df; on=:STAID, makeunique=true)
            end
        catch e
            @warn "Failed to load $filename: $e"
        end
    end

    # Select relevant columns
    relevant_cols = ["STAID", "NDAMS_2009", "MAJ_DDENS_2009", "STOR_NID_2009",
                     "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
                     "HYDRO_DISTURB_INDX", "CLASS"]
    available_cols = intersect(relevant_cols, names(result))

    if length(available_cols) > 0
        return select(result, available_cols)
    else
        return DataFrame()
    end
end


"""
    load_canadian_interference(hydat_path=CFG_HYDAT_PATH) -> DataFrame

Load Canadian HYDAT interference metadata.

For Canadian gages, returns RHBN and REGULATED status.
Note: This requires a pre-exported CSV since tidyhydat is an R package.

Parameters
----------
hydat_path : String or Nothing
    Path to HYDAT metadata CSV (pre-exported from tidyhydat)

Returns
-------
DataFrame
    Canadian interference columns
"""
function load_canadian_interference(hydat_path=CFG_HYDAT_PATH)
    if isnothing(hydat_path) || hydat_path == ""
        @warn "HYDAT path not configured - Canadian interference metadata unavailable"
        return DataFrame(
            gage_id = String[],
            RHBN = Union{Bool, Missing}[],
            REGULATED = Union{Bool, Missing}[],
            human_interference_class = String[]
        )
    end

    if !isfile(hydat_path)
        @warn "HYDAT file not found: $hydat_path"
        return DataFrame(
            gage_id = String[],
            RHBN = Union{Bool, Missing}[],
            REGULATED = Union{Bool, Missing}[],
            human_interference_class = String[]
        )
    end

    try
        df = CSV.read(hydat_path, DataFrame)

        # Ensure required columns exist
        if !("STATION_NUMBER" in names(df))
            @warn "HYDAT file missing STATION_NUMBER column"
            return DataFrame()
        end

        # Rename and select columns
        rename!(df, :STATION_NUMBER => :gage_id)

        # Calculate human_interference_class
        df.human_interference_class = map(row -> begin
            rhbn = get(row, :RHBN, missing)
            if !ismissing(rhbn) && rhbn == true
                "reference"
            elseif !ismissing(rhbn) && rhbn == false
                "non-reference"
            else
                "unknown"
            end
        end, eachrow(df))

        return df
    catch e
        @warn "Failed to load HYDAT data: $e"
        return DataFrame()
    end
end


"""
    enrich_signatures_with_metadata(signatures::DataFrame, metadata::DataFrame) -> DataFrame

Add human interference columns to signature output.

Parameters
----------
signatures : DataFrame
    Signature results with gage_id column
metadata : DataFrame
    Watershed metadata with interference columns

Returns
-------
DataFrame
    Signatures with interference columns added
"""
function enrich_signatures_with_metadata(signatures::DataFrame, metadata::DataFrame)
    if !CFG_INCLUDE_HUMAN_INTERFERENCE
        return signatures
    end

    # Ensure both have gage_id column
    if !("gage_id" in names(signatures)) || !("gage_id" in names(metadata))
        @warn "Missing gage_id column - cannot enrich with metadata"
        return signatures
    end

    # Select interference columns from metadata
    interference_cols = intersect(CFG_INTERFERENCE_COLUMNS, names(metadata))
    cols_to_join = vcat(["gage_id"], interference_cols)
    cols_available = intersect(cols_to_join, names(metadata))

    if length(cols_available) <= 1
        @warn "No interference columns found in metadata"
        return signatures
    end

    # Join metadata to signatures
    metadata_subset = select(metadata, cols_available)
    result = leftjoin(signatures, metadata_subset; on=:gage_id)

    return result
end

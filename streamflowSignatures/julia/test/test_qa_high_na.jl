# flagged_for_high_na denominator — regression for the 2026-09-04 finding.
#
# The benchmark runner builds results_df from Dict{String, Vector{Any}}, so every
# signature column has eltype Any. compute_qa_flags used to select denominator
# columns by `eltype(col) <: Union{Missing, Number}`, which silently dropped all of
# them and left only the 16 typed numeric metadata columns joined afterwards — the
# delivered products flagged every Canadian gage and nothing else. The fix selects
# columns by NAME (16 statistic suffixes + the per-gage scalar registry) and tests
# NA at the value level.

@testset "flagged_for_high_na denominator (2026-09-04)" begin
    # Manifest sanity
    @test length(StreamflowSignatures.CFG_QAQC_HIGH_NA_SUFFIXES) == 16
    @test length(StreamflowSignatures.CFG_QAQC_HIGH_NA_SCALARS) == 21
    @test "_mean" in StreamflowSignatures.CFG_QAQC_HIGH_NA_SUFFIXES && "_pettitt_post_mk_pval" in StreamflowSignatures.CFG_QAQC_HIGH_NA_SUFFIXES
    @test "elasticity_static" in StreamflowSignatures.CFG_QAQC_HIGH_NA_SCALARS && "drought_threshold_fixed_p30" in StreamflowSignatures.CFG_QAQC_HIGH_NA_SCALARS

    # 1. Runner-shaped frame: UNTYPED signature columns + typed metadata columns.
    d = Dict{String, Vector{Any}}(
        "gage_id"            => Any["g1", "g2", "g3", "g4"],
        "Qann_mean"          => Any[1.0, 1.0, 1.0, missing],
        "Qann_median"        => Any[1.0, NaN, 1.0, missing],      # NaN counts as NA
        "Q5_pettitt_cp_year" => Any[2000.0, missing, missing, missing],
        "elasticity_static"  => Any[1.2, missing, missing, missing],   # registered scalar
        "swe_max_mean"       => Any[missing, missing, 3.0, missing],   # NA-filled family
    )
    df = DataFrame(d)
    df.latitude        = [40.0, 41.0, 42.0, 43.0]                          # numeric metadata: excluded
    df.NDAMS_2009      = Union{Missing, Float64}[missing, missing, missing, 3.0]
    df.num_water_years = [20, 21, 22, 23]
    df.area_normalized = [true, true, true, false]                         # Bool <: Number: excluded
    df.gage_type       = ["Canada", "Canada", "USGS", "USGS"]

    @test all(c -> eltype(df[!, c]) == Any, ["Qann_mean", "Qann_median", "Q5_pettitt_cp_year"])
    @test Set(high_na_denominator_columns(names(df))) ==
          Set(["Qann_mean", "Qann_median", "Q5_pettitt_cp_year", "elasticity_static", "swe_max_mean"])

    out = compute_qa_flags(df)
    # NA fractions over the 5 signature columns: g1 1/5, g2 4/5, g3 3/5, g4 5/5
    @test out.flagged_for_high_na == [false, true, true, true]
    # g1 would have been flagged by the OLD metadata-only rule (NDAMS_2009 missing =
    # 1 of 4 numeric metadata columns = 25 %... make it decisive: add more missing)
    df_old = copy(df)
    df_old.MAJ_DDENS_2009 = Union{Missing, Float64}[missing, missing, missing, 1.0]
    df_old.STOR_NID_2009  = Union{Missing, Float64}[missing, missing, missing, 1.0]
    out_old = compute_qa_flags(df_old)   # 3 of 7 numeric metadata missing for g1 = 43 % — must NOT flag g1
    @test out_old.flagged_for_high_na == [false, true, true, true]
    # flag columns themselves never enter the denominator on a second pass
    @test compute_qa_flags(out).flagged_for_high_na == out.flagged_for_high_na

    # 2. Typed frame gives the identical answer
    dt = DataFrame(gage_id=["g1", "g2", "g3", "g4"],
                   Qann_mean=Union{Missing, Float64}[1.0, 1.0, 1.0, missing],
                   Qann_median=Union{Missing, Float64}[1.0, NaN, 1.0, missing],
                   Q5_pettitt_cp_year=Union{Missing, Float64}[2000.0, missing, missing, missing],
                   elasticity_static=Union{Missing, Float64}[1.2, missing, missing, missing],
                   swe_max_mean=Union{Missing, Float64}[missing, missing, 3.0, missing],
                   latitude=[40.0, 41.0, 42.0, 43.0],
                   NDAMS_2009=Union{Missing, Float64}[missing, missing, missing, 3.0])
    @test compute_qa_flags(dt).flagged_for_high_na == out.flagged_for_high_na

    # 3. No signature columns at all -> the flag stays false (never metadata-driven)
    meta_only = DataFrame(gage_id=["a", "b"], latitude=[1.0, 2.0],
                          NDAMS_2009=Union{Missing, Float64}[missing, missing],
                          RHBN=Union{Missing, Bool}[missing, true])
    @test isempty(high_na_denominator_columns(names(meta_only)))
    @test compute_qa_flags(meta_only).flagged_for_high_na == [false, false]

    # 4. Full-product-shaped header: 100 bases x 16 + 21 scalars selected, metadata/flags not
    bases = ["base$(i)" for i in 1:100]
    meta = ["gage_id", "latitude", "longitude", "basin_area", "gage_type", "num_water_years",
            "start_water_year", "end_water_year", "area_normalized", "NDAMS_2009", "MAJ_DDENS_2009",
            "STOR_NID_2009", "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL", "HYDRO_DISTURB_INDX",
            "CLASS", "RHBN", "REGULATED", "human_interference_class"]
    header = vcat(meta, [b * s for b in bases for s in StreamflowSignatures.CFG_QAQC_HIGH_NA_SUFFIXES],
                  StreamflowSignatures.CFG_QAQC_HIGH_NA_SCALARS, get_flag_columns())
    @test length(header) == 1653
    @test length(high_na_denominator_columns(header)) == 1621
end

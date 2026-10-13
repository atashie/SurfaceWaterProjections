"""
STANDARD PRODUCT #1, VERSION 2 (12 Oct 2026): WY 1993-2025 @ 60 % qualifying fraction,
on the REPROCESSED Daymet input and with the CLIMATE metrics family (guidelines 3.9,
2026-10-12): 39 new annual signatures (20 precipitation + 19 temperature) x 16 fields =
+624 columns, 1,653 -> 2,277. Supersedes `processedOuts_drought_28jul2026` (same window, config
and qualification; the gage set is streamflow-determined and therefore unchanged).

What differs from the superseded product, by column class (plan:
docs/plans/2026-10-12-climate-metrics-plan.md §6):
  - Q-only columns: identical values (same streamflow input, same machine class);
  - climate-dependent columns: the new input (every gage has a series; calendar 1980-2025,
    so WY 2024-2025 gain climate; last-bit input differences on the 5,965 basins the
    stale file covered);
  - 624 new climate columns; `flagged_for_high_na` regenerated (its denominator grows).

Inputs default to the verified LOCAL copies (the exFAT drive truncates files silently);
override with STREAMFLOW_DATA_PATH / _CLIMATE_PATH / _METADATA_PATH / _OUTPUT_DIR /
_GAGES_II_DIR — `docs/benchmarks/run_batched_julia.py run` sets all of them per batch.

NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (config reads at load).
Verify the precompiled build carries the climate family before trusting the output:
  julia --project=julia -e 'using StreamflowSignatures; println(CFG_CLIMATE_ENABLED)'
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_END_WATER_YEAR"] = "2025"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
get!(ENV, "STREAMFLOW_OUTPUT_PREFIX", "streamflow_1993_2025_60pct_climate_12oct2026")

const _LOCAL = get(ENV, "STREAMFLOW_LOCAL_INPUTS", joinpath(homedir(), "HISSS_data"))
get!(ENV, "STREAMFLOW_DATA_PATH",
     joinpath(_LOCAL, "daymet-replay-01oct2026", "inputs", "combined_streamflow_data_09feb2026.parquet"))
# Reprocessed Daymet (2026-09-30): 8,017 basins, calendar 1980-2025, six variables;
# md5 c059076088cd8abd963c64094468e90f, 5,134,193,404 bytes (EO_data_processing/README_DAYMET.md)
get!(ENV, "STREAMFLOW_CLIMATE_PATH",
     joinpath(_LOCAL, "daymet-processed-29sep2026", "daymet_1980_2025_29sep2026.parquet"))
get!(ENV, "STREAMFLOW_METADATA_PATH",
     joinpath(_LOCAL, "daymet-replay-01oct2026", "inputs", "combined_watershed_metadata_09feb2026.csv"))
get!(ENV, "STREAMFLOW_OUTPUT_DIR",
     joinpath(_LOCAL, "climate-rerun-12oct2026", "processedOuts_1993_2025_climate_12oct2026"))
# Without this the 11 GAGES-II / HYDAT human-interference columns are silently absent
# and the additivity diff against the reference product fails on "dropped columns".
get!(ENV, "STREAMFLOW_GAGES_II_DIR", joinpath(_LOCAL, "daymet-replay-01oct2026", "inputs", "gagesMetadata"))

mkpath(ENV["STREAMFLOW_OUTPUT_DIR"])

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

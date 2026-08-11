"""
PRODUCTION standard run #2, REGENERATED WITH DROUGHT (10 Aug 2026):
WY 1980-2025 ("entire period of record", user decision 2026-07-22) + 60% qualifying
fraction, now carrying the streamflow drought signature family (+165 columns, 1,488
-> 1,653). Supersedes the 22 Jul 2026 product in `processedOuts_1980_2025_22jul2026`,
whose feature stack was otherwise identical.

Companion to the WY 1993-2025 drought product; both standard products carry the
family after this run, so drought values stay comparable within each window (they
are record-dependent and must NOT be compared across the two windows).

Paths default to the thumbdrive as mounted on macOS (/Volumes/Untitled); override
with STREAMFLOW_DATA_PATH / _CLIMATE_PATH / _METADATA_PATH / _OUTPUT_DIR.

NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (config reads at load).
Verify the precompiled build actually carries the drought family before trusting the
output (see docs/DEVELOPMENT.md -> config-variant precompilation gotcha):
  julia --project=julia -e 'using StreamflowSignatures; println(CFG_DROUGHT_ENABLED)'
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1980"
ENV["STREAMFLOW_END_WATER_YEAR"] = "2025"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
get!(ENV, "STREAMFLOW_OUTPUT_PREFIX", "streamflow_1980_2025_60pct_11aug2026")

const _DRIVE = get(ENV, "STREAMFLOW_DRIVE", "/Volumes/Untitled")
get!(ENV, "STREAMFLOW_DATA_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "combined_streamflow_data_09feb2026.parquet"))
# REBUILT climate input (2026-08-11): the canonical `daymet_1980_2023.parquet` on the
# drive is truncated and unreadable (see CHANGELOG -> Known Issues). This file was
# rebuilt from the restored annual CSVs with
# docs/benchmarks/convert_daymet_csvs_to_parquet.py and validated by reproducing the
# WY 1993-2025 product (0 columns added/dropped, identical gage set, max |diff| 3.4e-13
# on 98 of 1,653 columns). That residual is CONSISTENT WITH a last-bit ingestion/
# serialization difference but its cause is NOT established -- the original parquet is
# unreadable, so it cannot be isolated. See that script's docstring.
get!(ENV, "STREAMFLOW_CLIMATE_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "daymet_1980_2023_rebuilt_10aug2026.parquet"))
get!(ENV, "STREAMFLOW_METADATA_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "combined_watershed_metadata_09feb2026.csv"))
get!(ENV, "STREAMFLOW_OUTPUT_DIR",
     joinpath(_DRIVE, "processedOuts_1980_2025_11aug2026"))
# Without this the 11 GAGES-II / HYDAT human-interference columns are silently absent
# (the accessor reads ENV at call time; a precompiled const would ignore it).
get!(ENV, "STREAMFLOW_GAGES_II_DIR", joinpath(_DRIVE, "gagesMetadata"))

mkpath(ENV["STREAMFLOW_OUTPUT_DIR"])

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

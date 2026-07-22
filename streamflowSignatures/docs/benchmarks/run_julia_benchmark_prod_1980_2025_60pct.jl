"""
PRODUCTION standard run #2 (22 Jul 2026): WY 1980-2025 ("entire period of record",
user decision 2026-07-22) + 60% qualifying fraction. Same committed feature stack as
standard run #1 (97e0be0): July features + 60% overall trend gate + snow
record-anchored decade gate. Own experiment folder per the one-folder convention.
Plan + validation gates: docs/plans/production_run_1980_2025_60pct_plan.md.
NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (config reads at load).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1980"
ENV["STREAMFLOW_END_WATER_YEAR"] = "2025"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "streamflow_1980_2025_60pct_22jul2026"
ENV["STREAMFLOW_OUTPUT_DIR"] = raw"D:/processedOuts_1980_2025_22jul2026"

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

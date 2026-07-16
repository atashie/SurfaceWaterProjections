"""
PRODUCTION rerun (14 Jul 2026): WY >= 1993 (no end cap) + 60% qualifying fraction,
with all July 2026 features (annual-values parquet, b=1 recession alpha, snow metrics,
area-normalized gate). EXACTLY the April `startIn1993_60pct` experiment window/filters
(user decision — plan §10), so startIn1993_60pct_signatures.csv serves as the strict
pre-features baseline. Outputs -> D:/processedOuts_14jul2026.
Plan + validation gates: docs/plans/production_rerun_1993_2022_60pct_plan.md.
NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (module loads config at compile time).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "streamflow_1993_60pct_14jul2026"
ENV["STREAMFLOW_OUTPUT_DIR"] = raw"D:/processedOuts_14jul2026"

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

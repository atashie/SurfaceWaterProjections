"""
Experiment: startIn1993 — analysis restricted to water years >= 1993.
All other config identical to baseline.
NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (module loads config at compile time).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "startIn1993"

# Ensure the 60% filter is NOT active
delete!(ENV, "STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION")

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

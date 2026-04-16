"""
Experiment: startIn1993-and-60pct — WY >= 1993 + require 60% of years qualifying.
NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (module loads config at compile time).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "startIn1993_60pct"

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

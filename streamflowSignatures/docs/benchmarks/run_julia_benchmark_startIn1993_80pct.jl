"""
Experiment: startIn1993-and-80pct — WY >= 1993 + require 80% of years qualifying.
NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (module loads config at compile time).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.80"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "startIn1993_80pct"

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

"""
PRODUCTION full-period run (15 Jul 2026): WY >= 1980 (no end cap) + 60% qualifying
fraction, with all July 2026 features. Identical to the validated 1993+ run
(run_julia_benchmark_prod_1993_60pct.jl, Codex results-review GO) except the start
year.

Outputs go to STREAMFLOW_OUTPUT_DIR (set it to the run's own experiment folder per
the one-folder convention; the base runner falls back to docs/benchmarks/ if unset).
Override STREAMFLOW_OUTPUT_PREFIX the same way.

NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (module loads config at compile time).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1980"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
get!(ENV, "STREAMFLOW_OUTPUT_PREFIX", "streamflow_1980_60pct_16jul2026")
# STREAMFLOW_OUTPUT_DIR is honored if set in the environment; no hardcoded default.

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

"""
⚠️ SUPERSEDED (2026-08-10) — this wrapper produced the 1,488-column standard product #1
in `D:/processedOuts_22jul2026`, which no longer carries the delivered feature set. The
current standard product #1 (1,653 columns, incl. the drought family) is built by
`run_julia_benchmark_drought_1993_2025_60pct.jl`. Kept as the historical record of the
22 Jul run; the Windows `D:` paths below also predate the macOS/ENV-override setup.

PRODUCTION standard run #1 (22 Jul 2026): WY 1993-2025 (EXPLICIT end cap — first
production exercise of STREAMFLOW_END_WATER_YEAR) + 60% qualifying fraction. The
first "standard output" window (user decisions 2026-07-21/22; HISSS manuscript
§2.2.2). Feature set = the committed July stack (annual-values parquet, b=1
recession alpha, snow metrics, area-normalized gate, 20-value stats floor) PLUS
the 60% overall trend gate (2026-07-21) and the snow record-anchored decade gate
(2026-07-22). Outputs -> D:/processedOuts_22jul2026.
Plan + validation gates: docs/plans/production_run_1993_2025_60pct_plan.md.
NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (config reads at load).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_END_WATER_YEAR"] = "2025"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "streamflow_1993_2025_60pct_22jul2026"
ENV["STREAMFLOW_OUTPUT_DIR"] = raw"D:/processedOuts_22jul2026"

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

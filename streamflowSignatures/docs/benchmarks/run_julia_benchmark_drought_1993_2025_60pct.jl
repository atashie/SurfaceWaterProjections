"""
DROUGHT VALIDATION run (28 Jul 2026): standard window #1 (WY 1993-2025 @ 60%
qualifying fraction) re-run with the new streamflow drought signature family
(`calculate_drought_metrics`, 10 bases + 5 threshold scalars = +165 columns).

Purpose: the ADDITIVITY GATE. Diffed against the delivered standard product #1
(`processedOuts_22jul2026/streamflow_1993_2025_60pct_22jul2026_signatures.csv`),
every pre-existing column must be unchanged — only `flagged_for_high_na` may shift,
because its denominator counts all signature fields — and the 165 new columns must
be populated. Unit tests cannot establish this: the orchestrator's per-signature
try/catch would turn an unexpected drought failure on some production gage into
silently missing columns (Codex MAJOR-4, 2026-07-27).

Feature set is otherwise IDENTICAL to the 22 Jul standard run #1, so the diff
isolates the drought family. Plan: docs/plans/2026-07-27-drought-signatures-plan.md
(§11 rollout, §15 review targets 7 and 12, §16 review record).

Paths default to the thumbdrive as mounted on macOS (/Volumes/Untitled); override
with STREAMFLOW_DATA_PATH / _CLIMATE_PATH / _METADATA_PATH / _OUTPUT_DIR.

NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (config reads at load).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_END_WATER_YEAR"] = "2025"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "streamflow_1993_2025_60pct_drought_28jul2026"

const _DRIVE = get(ENV, "STREAMFLOW_DRIVE", "/Volumes/Untitled")
get!(ENV, "STREAMFLOW_DATA_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "combined_streamflow_data_09feb2026.parquet"))
get!(ENV, "STREAMFLOW_CLIMATE_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "daymet_1980_2023.parquet"))
get!(ENV, "STREAMFLOW_METADATA_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "combined_watershed_metadata_09feb2026.csv"))
get!(ENV, "STREAMFLOW_OUTPUT_DIR",
     joinpath(_DRIVE, "processedOuts_drought_28jul2026"))
# Without this the 11 GAGES-II / HYDAT human-interference columns are silently absent
# and the additivity diff against the reference product fails on "dropped columns".
get!(ENV, "STREAMFLOW_GAGES_II_DIR", joinpath(_DRIVE, "gagesMetadata"))

mkpath(ENV["STREAMFLOW_OUTPUT_DIR"])

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

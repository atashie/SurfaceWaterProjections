"""
PORT-CAMPAIGN REFERENCE run (24 Aug 2026): standard window #1 (WY 1993-2025 @ 60%
qualifying fraction, drought family on) re-run AFTER the seasonal runoff-ratio
flag-name fix (runoff_ratios.jl win_/spr_/sum_/fal_complete — CHANGELOG Known
Issues, fixed 2026-08-24).

Purpose: the same-machine Julia reference for the Python/rpkg six-feature port
(docs/plans/2026-08-24-port-julia-features-to-python-rpkg-plan.md, Phase 0).
Diffed against delivered standard product #1
(processedOuts_drought_28jul2026/streamflow_1993_2025_60pct_drought_28jul2026_signatures.csv)
the expected delta is ONLY the four seasonal runoff-ratio stat/Pettitt families
(masking now fires) plus possibly `flagged_for_high_na`; the annual parquet's
seasonal runoff-ratio rows shift correspondingly (masking precedes collection —
Codex F4). Annual runoff ratio and `runoff_ratio_high_count` must NOT move.

Differences from the 28 Jul product wrapper: climate defaults to the REBUILT
Daymet parquet (the canonical-name file is truncated), own output folder/prefix,
and STREAMFLOW_HASH_INPUTS=1 for full provenance.

NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (config reads at load).
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_END_WATER_YEAR"] = "2025"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
ENV["STREAMFLOW_OUTPUT_PREFIX"] = "streamflow_1993_2025_60pct_portref_24aug2026"
get!(ENV, "STREAMFLOW_HASH_INPUTS", "1")

const _DRIVE = get(ENV, "STREAMFLOW_DRIVE", "/Volumes/Untitled")
get!(ENV, "STREAMFLOW_DATA_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "combined_streamflow_data_09feb2026.parquet"))
get!(ENV, "STREAMFLOW_CLIMATE_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "daymet_1980_2023_rebuilt_10aug2026.parquet"))
get!(ENV, "STREAMFLOW_METADATA_PATH",
     joinpath(_DRIVE, "processedOuts_feb2026", "combined_watershed_metadata_09feb2026.csv"))
get!(ENV, "STREAMFLOW_OUTPUT_DIR",
     joinpath(_DRIVE, "processedOuts_portref_24aug2026"))
# Without this the 11 GAGES-II / HYDAT human-interference columns are silently absent
# and the diff against the reference product fails on "dropped columns".
get!(ENV, "STREAMFLOW_GAGES_II_DIR", joinpath(_DRIVE, "gagesMetadata"))

mkpath(ENV["STREAMFLOW_OUTPUT_DIR"])

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

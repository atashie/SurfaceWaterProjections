"""
STANDARD PRODUCT #1, VERSION 3 (built 2026-10-27, verified 2026-10-28; public name
"version 2"): WY 1993-2025 @ 60 % qualifying fraction — the version-2 configuration
(reprocessed Daymet input, climate-metrics family, 2,277 columns) rerun on the CORRECTED
streamflow input with the gage set FROZEN at the version-2 cohort (6,678 gages; user
decision 2026-10-23). Supersedes the undeposited `processedOuts_1993_2025_climate_12oct2026`.

What differs from version 2 (plan docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md,
record docs/plans/2026-10-21-qualifier-rerun-results.md):
  - input `combined_streamflow_data_09feb2026_qualifier_patched.parquet`: the 193,542 USGS
    days the February compile masked under its qualifier rule (`A R`, `A <`, `A >`, `A [4]`)
    restored from re-fetched NWIS records; the metadata rows of the 987 patched gages
    updated in place (`combined_watershed_metadata_09feb2026_qualifier_patched.csv`);
  - climate frame `outer` (STREAMFLOW_CLIMATE_FRAME; the default since 2026-10-28): each
    gage's Daymet rows enter on Daymet's own calendar inside the record span, so absent
    streamflow days keep their climate and precipitation totals equal the daily sums;
  - everything else identical: the 5,069 gages untouched by both corrections are
    bit-identical to version 2.

HOW THE PRODUCT WAS BUILT: through `docs/benchmarks/run_batched_julia.py` (8 gage batches
from the FROZEN selection manifest `manifests/selection_1993-2025.txt`, `--climate-frame
outer`, `--qa-schema`, `--stage S2`, `--require-clean-tree`) from the clean worktree of the
tag `qualifier-rerun-baseline-2026-10-26`, with this wrapper as `--wrapper` — the batched
runner overrides the input, climate, metadata and output paths per batch and merges the
batches under the product prefix. Run ALONE on the full input this wrapper selects gages by
qualification and therefore ADDS the 31 gages that qualify only under the corrected record
(listed in the product's `would_qualify_not_added` file): reproduce the product through the
batched runner with the selection manifest, never standalone.

The delivered product folder is never written by this wrapper: it refuses to run when
STREAMFLOW_OUTPUT_DIR is unset and that folder already exists — set STREAMFLOW_OUTPUT_DIR
to run elsewhere.

Inputs default to the verified LOCAL copies (the exFAT drive truncates files silently);
override with STREAMFLOW_DATA_PATH / _CLIMATE_PATH / _METADATA_PATH / _OUTPUT_DIR /
_GAGES_II_DIR — `run_batched_julia.py run` sets them per batch.

NOTE: ENV vars must be set BEFORE `using StreamflowSignatures` (config reads at load).
Verify the precompiled build carries the climate family before trusting the output:
  julia --project=julia -e 'using StreamflowSignatures; println(CFG_CLIMATE_ENABLED)'
"""
ENV["STREAMFLOW_START_WATER_YEAR"] = "1993"
ENV["STREAMFLOW_END_WATER_YEAR"] = "2025"
ENV["STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION"] = "0.60"
# The product's climate frame (also the helpers' default since 2026-10-28); `get!` so the
# batched runner's --climate-frame (e.g. `left` for an S0 baseline) still wins.
get!(ENV, "STREAMFLOW_CLIMATE_FRAME", "outer")
get!(ENV, "STREAMFLOW_HASH_INPUTS", "1")
get!(ENV, "STREAMFLOW_OUTPUT_PREFIX", "streamflow_1993_2025_60pct_climate_v3_26oct2026")

const _LOCAL = get(ENV, "STREAMFLOW_LOCAL_INPUTS", joinpath(homedir(), "HISSS_data"))
const _RERUN = joinpath(_LOCAL, "qualifier-rerun-26oct2026")
# Corrected streamflow input (splice of 2026-10-26; sha256 2d44c1af…, input/provenance.json)
get!(ENV, "STREAMFLOW_DATA_PATH",
     joinpath(_RERUN, "input", "combined_streamflow_data_09feb2026_qualifier_patched.parquet"))
# Reprocessed Daymet (2026-09-30): 8,017 basins, calendar 1980-2025, six variables;
# md5 c059076088cd8abd963c64094468e90f, 5,134,193,404 bytes (EO_data_processing/README_DAYMET.md)
get!(ENV, "STREAMFLOW_CLIMATE_PATH",
     joinpath(_LOCAL, "daymet-processed-29sep2026", "daymet_1980_2025_29sep2026.parquet"))
# Metadata with the 987 patched gages' num_years / num_days / end_year updated (sha256 ca335b2d…)
get!(ENV, "STREAMFLOW_METADATA_PATH",
     joinpath(_RERUN, "input", "combined_watershed_metadata_09feb2026_qualifier_patched.csv"))
# Without this the 11 GAGES-II / HYDAT human-interference columns are silently absent
# and the additivity diff against the reference product fails on "dropped columns".
get!(ENV, "STREAMFLOW_GAGES_II_DIR", joinpath(_LOCAL, "daymet-replay-01oct2026", "inputs", "gagesMetadata"))

# Delivered products are never rewritten (CLAUDE.md): refuse the product folder itself.
const _PRODUCT_DIR = joinpath(_RERUN, "processedOuts_1993_2025_climate_v3_26oct2026")
if !haskey(ENV, "STREAMFLOW_OUTPUT_DIR") && isdir(_PRODUCT_DIR) && !isempty(readdir(_PRODUCT_DIR))
    error("refusing to run: $_PRODUCT_DIR holds the delivered product — set STREAMFLOW_OUTPUT_DIR to run elsewhere")
end
get!(ENV, "STREAMFLOW_OUTPUT_DIR", _PRODUCT_DIR)
mkpath(ENV["STREAMFLOW_OUTPUT_DIR"])

include(joinpath(@__DIR__, "run_julia_benchmark.jl"))

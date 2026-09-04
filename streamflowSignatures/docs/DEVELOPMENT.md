# Development Guide

This guide covers development workflows, architecture, and common tasks for the Streamflow Signatures project.

## Architecture

### Data Flow

```
Data Sources                    Processing                     Output
────────────                    ──────────                     ──────
USGS (dataRetrieval)  ──┐
  (streamflow only)     ├──> Parquet Storage ──┐
                        │                       │
Canadian HYDAT  ────────┘                       ├──> Signature ──> CSV Summary
  (streamflow only)                             │    Extraction     (1653 columns)
                                                │        │
Daymet (climate data)  ──> Parquet Storage ─────┘        └──> Annual Values Parquet
  (PPT, temp, SWE)         (joined at runtime)                (gage_id, signature,
                                                               water_year, value)

Alternative Pipeline:
Caravan NetCDF ──────────────> Direct Processing ──> Caravan Output
  (bundled Q + climate)        (annualized CSVs)
```

### Design Principles

1. **Plain-English Guardrails** — Domain experts define signature methodology in `SIGNATURE_GUIDELINES.md` (auto-synced from a shared Google Doc). Code implements those definitions. This separates hydrological expertise from implementation. Note: Notion was evaluated as an alternative hosting platform (March 2026) but rejected — Notion's JS-rendered published pages are unreliable for automated fetching, and Google Docs' static HTML + permanent URLs are better suited to the auto-sync workflow.

2. **Julia Canonical, Others Follow** — All changes are made in Julia first, validated via benchmark (~27 min), then propagated to Python and rpkg. Golden outputs from Julia serve as the reference for cross-language validation.

3. **Strict Output Schema** — The CSV output format (column names, ordering) is a contract. Downstream tools depend on exact column names. Every signature produces exactly 8 statistics via `generate_stats()`.

4. **Per-Year Quality Filtering** — Each water year is independently evaluated by `preprocess_daily_data()` against data completeness thresholds (>30 raw NAs, >3-day gaps, residual boundary NAs). Negative Q rejection is config-driven (`reject_negative_flow: false` by default). Constant-SD is a QA flag only. Years that fail are excluded; gages need 20+ qualifying years.

5. **Centralized NA Handling** — Missing data is handled once per gage via `preprocess_daily_data()` before any signatures are computed. This replaces ad-hoc per-signature NA handling (fillna(0), forward-fill, etc.) with a standardized pipeline: daily grid normalization → interpolation (<=3 day internal gaps) → year rejection → residual check. Configuration is in `config/signatures_config.json` under `na_handling`.

### NA Handling Architecture

```
Raw gage data (may have NAs, gaps, duplicates)
    │
    ▼
preprocess_daily_data()           ◄── Called ONCE per gage, BEFORE all signatures
    │
    ├── 1. Daily grid normalization (one row/day, sorted, unique)
    ├── 2. Raw diagnostics (NA count, max run, seasonal completeness, negative check)
    ├── 3. Year rejection (>30 raw NAs, >3-day gaps, negative Q if config-enabled)
    ├── 4. Interpolation (internal gaps <=3 days, linear, no extrapolation)
    ├── 5. Residual check (boundary NAs → reject year)
    ├── 6. PPT handling (same rules, tracked separately)
    ├── 7. SWE handling (same rules, tracked separately — July 2026)
    │
    └── Returns: cleaned data, valid_years, valid_climate_years,
                 valid_swe_years, seasonal_flags, diagnostics, rejected_years
    │
    ▼
Signature functions receive clean data
    ├── Flow volumes, timing, baseflow, recession, etc.
    ├── Climate signatures use valid_climate_years subset
    ├── Snow signatures use an EXPLICIT snow_data frame filtered to valid_swe_years
    │   (opt-in — never derived implicitly from an SWE column in the gage frame)
    ├── Drought signatures smooth Q within contiguous date runs (never across a
    │   rejected-year gap) before thresholding — July 2026; all 3 languages
    └── Seasonal signatures respect seasonal_flags (incomplete → NA)
```

**Key design decisions:**
- `config/signatures_config.json` → `na_handling` section is the single source of truth
- `use_legacy_filtering: false` — new preprocessing is the default
- Negative Q rejection is conditional on `reject_negative_flow` (default: false); `negative_ann` signature counts Q<0 days instead
- Constant-SD detection is a QA flag, not a year rejection criterion
- `ice_affected_days_total` per-gage output aggregates ice-related NA days from diagnostics
- Seasonal completeness is computed from RAW observations (pre-interpolation)
- `generate_stats()` has optional `trend_completeness` / `decade_completeness` params

### Annual Values Export (July 2026; all three languages since August 2026)

Every signature's per-year annual values — previously discarded after
`generate_stats()` collapsed them into the 8 statistics — are exported as one
long-format parquet alongside the summary CSV:
`{output_dir}/{prefix}_signatures_annual.parquet` — i.e. the run's own experiment folder
(the one-folder convention, CLAUDE.md Critical Constraint #5), NOT `docs/benchmarks/`.

| Column | Type | Notes |
|---|---|---|
| `gage_id` | String | Zero-padded original ID (same as summary CSV) |
| `signature` | String | Signature base name (matches summary column prefixes) |
| `water_year` | Int32 | Oct 1 – Sep 30. `elasticity_rolling` = END year of the 11-yr window; `elasticity_annual` = later year of the consecutive pair |
| `value` | Float64 | Annual value; **NaN and absent-row are equivalent** ("not computable that year") |

**Mechanism**: an opt-in `AnnualCollector` (defined in `julia/src/stats.jl`) is
threaded through `calculate_all_signatures()` and all 14 signature functions into
`generate_stats()`, which appends the incoming series — exactly as passed, BEFORE
any min_rows or trend-completeness gating — as long-format triplets. With
`collector=nothing` (the default) behavior is byte-identical to before; the summary
CSV contract is untouched.

**Semantics caveat**: the export records the exact series `generate_stats()`
consumed. Dense signatures (flow volumes, FDC, timing, pulses, runoff ratios,
recession base metrics) carry NaN placeholder rows; caller-pruned signatures
(`qp_bimodality`, `n_recession_events`) and sparse-built frames (flashiness,
storage, baseflow, elasticity) omit non-computable years entirely.

**Config**: `config/signatures_config.json` → `annual_values.save` (repo default
`true`; absent section → disabled). Julia constant: `CFG_SAVE_ANNUAL_VALUES`.

**Validation**: `julia/test/test_annual_collector.jl` (collector has zero effect on
stats, zero warnings, signature coverage, NaN/missing-year guards, no duplicate
keys); `docs/benchmarks/validate_annual_values.py` recomputes mean/median from the
parquet and cross-checks the summary CSV after each benchmark run.

**Ports**: implemented in Python and rpkg by the August 2026 port campaign.
The collector kwarg still defaults to nothing/NULL, so passing no collector is
byte-identical to the pre-collector behaviour in every language. Python's
collector is a small class over list accumulators; rpkg's is an environment
with pre-grown chunk lists (`c()` append would be O(n^2) in R). Cross-language
validation of the exported parquet: identical 100-signature sets, 0 NA-pattern
mismatches, 18,898,405 of 18,898,406 rows shared with the Julia reference
(docs/CROSS_LANGUAGE_STATUS.md).

Design + Codex review record: `docs/plans/annual_values_export_plan.md`.

## File Structure

```
streamflowSignatures/
├── README.md                    # User entry point
├── CHANGELOG.md                 # Bug fixes, version history
├── CLAUDE.md                    # Claude Code instructions
├── .gitignore                   # Git ignore patterns
│
│ ## R Workflow (unchanged at root for current users)
├── config.R                     # Centralized configuration parameters
├── run_full_processing.R        # PRIMARY - Full signature extraction with climate data
├── run_ingest_usgs_hydat.R      # Raw USGS/HYDAT data ingestion to parquet
├── run_caravan_processing.R     # Caravan data processing
├── run_restricted_processing.R  # Restricted processing
│
├── R/                           # Legacy R (deprecated for signatures — ingestion utilities still active)
│   ├── helperFunctions.R        # DEPRECATED — Legacy shim (ingestion utilities only)
│   ├── load_config.R            # Config loader
│   ├── run_conversion.R         # Daymet ZIP to Parquet conversion
│   ├── run_enrich_metadata.R    # Human interference metadata enrichment
│   ├── run_regenerate_metadata.R # Regenerate combined metadata
│   ├── precompute_cross_signature_analysis.R  # Offline computation
│   ├── export_hydat_metadata.R  # One-time: export HYDAT RHBN/REGULATED to CSV
│   └── tests/                   # R test suite
│       ├── smoke_test.R         # Quick validation on subset (10 gages)
│       ├── smoke_test_reorganization.R
│       ├── qa_qc_signatures.R   # Output validation and QA/QC checks
│       ├── visualize_qa_qc.R    # QA/QC visualization plots
│       ├── test_climate_functions.R # Climate function tests
│       ├── test_climate_signatures.R
│       ├── test_na_handling.R       # Legacy NA-handling suite (6 known pre-existing failures)
│       ├── generate_golden_outputs.R
│       └── verify_no_regression.R  # Golden output regression check
│
├── rpkg/                        # R port (production-ready, mirrors Julia/Python structure)
│   ├── DESCRIPTION
│   ├── NAMESPACE
│   ├── README.md
│   ├── R/                       # 20 modules (+ zzz.R)
│   └── tests/                   # testthat tests + real-data smoke test
│
├── python/                      # Python package (production-ready)
│   ├── README.md
│   ├── pyproject.toml
│   ├── streamflow_signatures/   # 20 modules
│   └── tests/                   # Python tests
│
├── julia/                       # Julia canonical implementation (production-ready)
│   ├── README.md
│   ├── Project.toml
│   ├── src/                     # 20 modules + StreamflowSignatures.jl
│   └── test/                    # Julia tests
│
├── config/                      # Cross-language configuration
│   └── signatures_config.json
│
├── golden-outputs/              # Julia canonical (April 2026) + R historical (Feb 2026)
│   ├── README.md
│   ├── streamflow_signatures_full_10feb2026.csv
│   └── combined_watershed_metadata_09feb2026.csv
│   # Julia is the canonical golden output (April 2026). R golden outputs (Feb 2026)
│   # are historical — they pre-date trend_completeness, Section 3 signatures,
│   # and recession fix.
│
├── docs/                        # Extended documentation
│   ├── DEVELOPMENT.md           # This file
│   ├── SIGNATURES.md            # Detailed signature documentation
│   ├── SIGNATURE_GUIDELINES.md  # Collaborative guidelines (auto-synced)
│   ├── WORKFLOW_REVIEW.md       # Workflow review
│   ├── CROSS_LANGUAGE_STATUS.md # Cross-language alignment detail
│   ├── CODE_REVIEW.md          # Cross-language code review findings
│   ├── benchmarks/              # Benchmark runners and results
│   │   ├── run_python_benchmark.py
│   │   ├── run_julia_benchmark.jl
│   │   ├── run_r_benchmark.R
│   │   ├── compare_three_way.py
│   │   ├── compare_julia_vs_golden_r.py  # 6-tier Julia vs Golden R comparison
│   │   ├── build_julia_vs_golden_r_dashboard.py  # Interactive HTML dashboard
│   │   ├── build_section3_dashboard.py   # Section 3 pre/post dashboard
│   │   ├── compare_experiment_vs_julia.py  # Parameterized experiment comparison
│   │   ├── build_experiment_vs_julia_dashboard.py  # Experiment dashboard builder
│   │   ├── build_new_vs_golden_julia_dashboard.py  # New benchmark vs golden Julia validation
│   │   ├── run_julia_benchmark_startIn1993.jl  # Experiment: WY >= 1993
│   │   ├── run_julia_benchmark_startIn1993_60pct.jl  # Experiment: WY >= 1993 + 60% filter
│   │   ├── run_julia_benchmark_startIn1993_80pct.jl  # Experiment: WY >= 1993 + 80% filter
│   │   ├── comparison_report.md
│   │   ├── julia_vs_golden_r_summary.md  # Generated comparison report
│   │   └── diagnostics/        # Archived diagnostic scripts
│   └── plans/                   # Planning notes
│
├── claude-skill/                # Claude AI skill
│   └── streamflow-signatures.md
│
├── streamflowAndClimateVisualizationApp/  # Shiny dashboard
│   ├── app.R                    # Main Shiny application
│   └── helperFunctions.R        # App-specific utilities
│
├── metadata/                    # Basin and gage metadata (43 files)
│   └── canadian_hydat_interference.csv  # RHBN + REGULATED for 8,012 Canadian stations
│
├── data_out/                    # Processed outputs (gitignored)
├── test_output/                 # Test outputs (gitignored)
└── logs/                        # Processing logs (gitignored)
```

## Common Tasks

### Publish to the public mirror (CZ-Sync/HISSS)

The public code repository for the HISSS manuscript is
**https://github.com/CZ-Sync/HISSS** — a snapshot mirror of this project's
git-tracked content, kept current by:

```bash
./publish_to_hisss.sh
```

The script copies exactly the tracked files minus a fixed exclusion list
(internal/Claude tooling, the unpublished manuscript snapshot, `docs/plans/`,
the defunct Shiny app, and the two >50 MB golden-output CSVs — see the script
header), commits them on HISSS `main` with the source revision in the message,
and pushes. **Run it after every merge to main.** Never commit to HISSS
directly — the next publish overwrites it.

### Run Signature Extraction

**Julia (canonical)**:

```julia
using StreamflowSignatures
# See docs/benchmarks/run_julia_benchmark.jl for full pipeline
```

**R (via rpkg)**:

```r
library(streamflowsignatures)
# See docs/benchmarks/run_rpkg_benchmark.R for full pipeline
```

### Run Full Processing Pipeline

The fastest way to run a complete signature extraction:

```bash
# Julia (canonical, ~27 min for 7,313 gages):
julia docs/benchmarks/run_julia_benchmark.jl

# R legacy (deprecated — use Julia or rpkg instead):
Rscript run_full_processing.R
```

**Prerequisites:** Data paths are configured in `config/signatures_config.json` and `config.R`.

### Run a Config-Variant Experiment (⚠️ precompilation gotcha)

`STREAMFLOW_CONFIG` points the Julia package at an alternative
`signatures_config.json`, but **on its own it does nothing to an already-precompiled
package**. Every `CFG_*` constant in `julia/src/config.jl` is evaluated when the module is
PRECOMPILED, and Julia does not invalidate the precompile cache when an environment
variable changes — so the run silently uses the *previous* config. Verified 2026-07-28: a
benchmark launched with a drought-disabled config produced an output byte-identical in
size to the drought-enabled run, all 165 drought columns still present.

Force a recompile so the new config is read. **`touch` is NOT enough** — Julia validates
the cache by file *content*, so an mtime-only change is ignored (verified 2026-07-28: a
`touch` + `STREAMFLOW_CONFIG` run still produced all the drought columns). Delete the
package's compiled cache instead:

```bash
rm -rf ~/.julia/compiled/v1.12/StreamflowSignatures       # adjust the version dir
STREAMFLOW_CONFIG=/path/to/variant_config.json \
  julia --project=julia docs/benchmarks/run_julia_benchmark.jl
```

(An actual source edit also invalidates it, and `julia --compiled-modules=no` works but
makes every run much slower.) **Verify the variant took effect before trusting the
result** — the cheap probe is:

```bash
STREAMFLOW_CONFIG=/path/to/variant_config.json \
  julia --project=julia -e 'using StreamflowSignatures; println(CFG_DROUGHT_ENABLED)'
```

Note the cache then holds the VARIANT build: purge it again afterwards, or subsequent
"normal" runs silently keep using the variant config.

**Audit of past config-variant results (2026-07-28)**: no previously recorded result is
invalidated by this discovery. The July 2026 60 % trend-gate change has empirical output
evidence that the gate took effect (gained 45 / lost 0 trend gages), and the snow gates
involved source-file changes, which do invalidate the cache. But **any future or
historical experiment that relied only on switching `STREAMFLOW_CONFIG`** — without a
cache purge, a source-content change, a config probe, or an observed expected delta — is
untrustworthy and should be re-run. The durable fix would be to read config values at
runtime rather than into module constants. This is the same
hazard that led the benchmark runner to read `STREAMFLOW_START_WATER_YEAR`,
`STREAMFLOW_END_WATER_YEAR`, and `STREAMFLOW_MIN_QUALIFYING_DATA_FRACTION` at *runtime*
inside `main()` instead of from module constants.

### Verify input integrity before a production run

The inputs live on an **exFAT** thumbdrive and can be silently truncated: on 2026-08-10
`daymet_1980_2023.parquet` was found at 1,261,436,928 bytes instead of 4,125,630,653 —
**with an unchanged mtime**, so only the byte count betrayed it (Julia then failed at
Phase 2 with `invalid parquet: final bytes … expect "PAR1"`). Check sizes against the
`provenance` block of the most recent timing JSON, and check that each parquet still ends
in its magic bytes, before committing to a long run:

```bash
for f in combined_streamflow_data_09feb2026.parquet daymet_1980_2023.parquet; do
  p=/Volumes/Untitled/processedOuts_feb2026/$f
  echo "$f  $(stat -f%z "$p")  footer=$(tail -c 4 "$p")"
done
```

This is what the provenance block (July 2026) exists for — record it on every run.

**If the climate parquet is the casualty**, rebuild it from the 44 annual Daymet CSVs
(no R needed):

```bash
python docs/benchmarks/convert_daymet_csvs_to_parquet.py \
  --csv-dir /Volumes/Untitled/daymet_1980_2023 \
  --out    /Volumes/Untitled/processedOuts_feb2026/daymet_1980_2023_rebuilt_10aug2026.parquet
```

then prove the rebuilt input before building anything new on it, by re-running a
delivered product's config against it and diffing:

```bash
julia --project=julia docs/benchmarks/check_additivity.jl NEW.csv DELIVERED.csv --expect-added 0
```

The 2026-08-11 rebuild reproduced the WY 1993–2025 product with 0 columns added/dropped,
an identical gage set, and ≤ 3.4e-13 on 98 of 1,653 columns. Note the source CSVs carry
**no day column** and Daymet uses a **365-day calendar** (leap years drop Dec 31, not
Feb 29) — see the converter's docstring and CHANGELOG → August 2026.

### Validate Output Quality

After processing, run QA/QC validation:

```r
source("config.R")
source("R/helperFunctions.R")
source("R/tests/qa_qc_signatures.R")
```

Or run the visualization script for diagnostic plots:

```r
source("R/tests/visualize_qa_qc.R")
# Outputs to data_out/qa_plots/
```

QA/QC checks include:
- Range validation (e.g., BFI in [0,1])
- Baseflow consistency (BFI_Eckhardt < BFI_LyneHollick)
- Elasticity constraints
- Correlation checks between related metrics

### Process Caravan Data

```r
source("R/helperFunctions.R")

process_caravan_to_annual(
  caravan_directory = "path/to/caravan/netcdf",
  data_project = "camels",
  min_num_years_data = 30,
  start_date_filter = as.Date("1979-09-01"),
  end_date_filter = as.Date("2025-06-01"),
  output_dir = "annualized_caravan_data"
)
```

## Watershed HydroATLAS Metadata

A standalone, per-gage metadata product describing the **hydro-geophysical character of
each gage's entire upstream watershed** (climate, hydrology, terrain, land cover,
soils/geology, anthropogenic), drawn from HydroATLAS / BasinATLAS v10. It sits alongside
the signature outputs and joins by `gage_id`. This is metadata/ingestion work (R) — **not**
a cross-language signature, so it is not ported to Julia/Python.

### Run

```bash
Rscript run_hydroatlas_metadata.R                 # full run (success gages) -> data_out/
Rscript run_hydroatlas_metadata.R --subset 40     # dry-run on 40 spread gages -> test_output/
Rscript run_hydroatlas_metadata.R --status all    # all gages with valid coordinates
Rscript run_hydroatlas_metadata.R --metadata <path/to/combined_watershed_metadata.csv>
```

Outputs (in `data_out/`): `watershed_hydroatlas_metadata_{date}.csv` + `.parquet` (~211
columns, one row per gage), `watershed_hydroatlas_metadata_dictionary.csv`
(column → theme / unit / aggregation / description), and a reusable
`hydroatlas_member_attrs_{date}.parquet` attribute cache. `combined_watershed_metadata.csv`
is left untouched.

### How it works

1. **gage → outlet basin** — one vectorized `st_join` of gage points against
   `basinAt_NorAm_polys.gpkg` (HydroBASINS lev12, 167k NorAm basins, WGS84). s2 is disabled
   (GEOS planar) to tolerate minor invalid loops in the source polygons.
2. **outlet → watershed** — member basin set (outlet + all upstream) from the cached
   `upstream_hydrobasins.rds` (keyed by outlet HYBAS_ID); outlets not yet cached are filled
   via a `NEXT_DOWN` BFS over the topology.
3. **aggregate (hybrid)** — rules built programmatically by `classify_hydroatlas_attributes()`:
   - `_u` / `_p` attributes (91): **passthrough** from the outlet basin. The delineated
     upstream set equals the extent HydroATLAS already accumulates into the outlet basin's
     `_u` fields, so this is the rigorous watershed value (no re-aggregation).
   - `_s`-only attributes (105): **SUB_AREA-weighted** — area-weighted mean for continuous,
     percentage, and monthly-climate fields (HydroATLAS `-9999` NoData masked; percentages
     always mean, never max); spatial **min/max** for elevation only; area-weighted
     **majority** for categoricals (glc/pnv via argmax of the outlet's upstream `_u` fractions;
     wet and the other 8 via area-weighted mode of the source class).
   - `_s` attributes that have a `_u` twin (85) are dropped in favor of the `_u` value.
4. **assemble + write** — one row per gage, keyed by zero-padded `gage_id`. Join to the
   signatures with `join_hydroatlas_to_signatures()` (leading-zero-safe canonical matching).

### Memory / size

Aggregation runs per **unique outlet** (~7,100) over O(1) keyed member slices; the
gage × basin × attribute long table is never materialized. The attribute table is subset to
only the member basins actually used (~96k) and cached as parquet. The full output is ~13 MB
CSV — all artifacts stay far under 1 GB.

Config/source: `HYDROATLAS_GPKG` in `config.R`; module `R/aggregate_hydroatlas_metadata.R`;
entry `run_hydroatlas_metadata.R`. Validation (`validate_hydroatlas_metadata()`) cross-checks
member area vs. outlet `UP_AREA` and area-weighted `_s` vs. HydroATLAS `_u`.

## Static HTML Explorer

A self-contained, double-clickable map for exploring the gages and their watersheds:
`data_out/streamflow_explorer.html` (Leaflet + Canvas). 8,014 gage points are colored by a
selected variable; a dataset switcher toggles between the **HydroATLAS watershed metadata**
(~207 vars) and the **streamflow signatures** (~205 vars: `_mean`, `_senn_slp`, `_mk_pval` per
base + scalars). Clicking a gage draws its **watershed boundary** (border only).

### Build

```bash
Rscript build_explorer.R     # dissolve+simplify member basins per outlet -> data_out/watershed_borders.geojson (~12 MB)
Rscript assemble_explorer.R  # join points + inject into explorer_template.html -> data_out/streamflow_explorer.html (~31 MB)
```

- `build_explorer.R` reuses `build_upstream_members()` (cache + NEXT_DOWN BFS) so every gage
  outlet (incl. those not in the cached RDS) gets a border. It **unions the RAW lev12 basins**
  (identical shared edges → clean dissolve, no sliver holes/spikes — do NOT pre-simplify each
  basin first), strips interior holes + tiny sliver parts, simplifies the clean outline to ~1 km
  (`preserveTopology=TRUE`) with an adaptive **~250-vertex cap** on large basins, then runs a final
  `st_make_valid` + hole-strip on the rounded GeoJSON. Flags: `--tol`, `--maxv`, `--n` (subset).
- `assemble_explorer.R` rounds values, builds a column-oriented points object, and injects
  `__POINTS__` / `__VARS__` / `__BORDERS__` into `explorer_template.html`.

### Size note

Full-resolution borders for all ~7,100 *nested* watersheds are **~111 MB** (un-inline-able); the
clean-dissolve + simplified + vertex-capped build is **~9.5 MB**, giving a ~28 MB self-contained HTML. The page
needs internet to load Leaflet + the CARTO basemap from CDN (points/borders still render offline
once loaded). For a lighter file, trim the variable list in `assemble_explorer.R` or lower `--maxv`.

## Adding a New Signature

### Step-by-Step Process (Julia-First Workflow)

1. **Create calculation function** in the appropriate `julia/src/` module file
   - Function should accept daily data and return annual values as a DataFrame
   - Use the existing module structure (e.g., `flow_volumes.jl`, `pulses.jl`, `recession.jl`)

2. **Apply the 8-statistic rule** using `generate_stats()`:
   ```julia
   # Your function should return annual values, then call:
   stats = generate_stats(annual_df, [:metric_name], :water_year)
   # This produces 8 columns: metric_senn_slp, metric_linear_slp,
   # metric_spearman_rho, metric_spearman_pval, metric_mk_rho,
   # metric_mk_pval, metric_mean, metric_median
   ```

3. **Add call in `julia/src/signatures.jl`** in the `calculate_all_signatures()` function

4. **Register the signature** in `config/signatures_config.json` and `config.R` (for R tests):
   - Add base name to `EXPECTED_SIGNATURE_BASES` in `config.R`
   - Per-gage scalars (not 8-stat) need their own `EXPECTED_*` constant in `config.R`,
     wired into `validate_output_schema()`
   - `EXPECTED_DENSE_SIGNATURES` in `julia/test/test_annual_collector.jl` — this asserts
     **set equality** of the collected annual series, so a new dense signature fails the
     suite until it is listed there
   - The signature-count gate in `docs/benchmarks/validate_production_run.py`
     (`ann.signature.nunique() == N`)
   - The **expected total summary-column count** in the docs — count the 8 stat + 8
     Pettitt fields per base AND any non-8-stat scalars

5. **Run the unit suite, then the Julia benchmark** (~27 min) to validate:
   ```bash
   julia --project=julia julia/test/runtests.jl
   julia docs/benchmarks/run_julia_benchmark.jl
   ```

6. **Prove additivity against the previous canonical run** — every pre-existing column
   unchanged (in principle `flagged_for_high_na` may shift because its denominator is
   meant to count all fields — but see CHANGELOG → Known Issues, 2026-09-04: the Julia
   runner feeds it `Vector{Any}` columns, so no signature field currently enters the
   denominator and the flag cannot move), and the new columns fully populated. The orchestrator catches per-signature
   exceptions, so a failure on some production gage shows up as silently missing columns
   that a green unit suite will not reveal. Smoke checks should require new values to be
   FINITE, not merely present as keys.

6. **Port to Python and rpkg**:
   - Python: Add function in the appropriate `python/streamflow_signatures/` module, call from `signatures.py`
   - rpkg: Add function in the appropriate `rpkg/R/` module, call from `signatures.R`, export in `NAMESPACE`

7. **Run cross-language comparison** to verify alignment:
   ```bash
   python docs/benchmarks/compare_three_way.py
   python docs/benchmarks/compare_rpkg.py
   ```

### Output Column Naming Convention

All signatures follow the pattern: `{metric}_{stat}`

Examples:
- `Qann_senn_slp` - Theil-Sen slope for annual mean flow
- `BFI_Eckhardt_mk_pval` - Mann-Kendall p-value for Eckhardt baseflow index

## Testing

### Quick Validation (Smoke Test)

```r
source("R/tests/smoke_test.R")
# Runs on 10 gages, validates output schema
```

### Climate Function Tests

```r
source("R/tests/test_climate_functions.R")
# Tests climate signatures with synthetic data
```

### Unit Tests

```r
# Run all tests in R/tests/ directory
testthat::test_dir("R/tests/")
```

## Cross-Language Benchmarks

Python, rpkg, and R canonical implementations are validated against Julia golden outputs using the R² of the identity line (y = x). This measures whether implementations produce identical values (not just correlated). Spearman rank correlation is reported as a secondary diagnostic.

### Current Status (August 2026 — FULL PARITY)

The August 2026 port campaign closed the six-feature gap: Python and rpkg now
produce the same **1,653-column** product as canonical Julia (Pettitt fields,
20-value stats floor, annual-values collector, b=1 recession alpha, snow,
drought), and **both are validated at full scale** — Python 2026-08-26, rpkg
2026-08-27 — each passing strict schema equality with no waivers, the
swallowed-family-failure gate, and the cross-language annual-parquet gate over
the WY 1993–2025 window:

| | columns | gages | identity-R² tiers | mean R² | annual rows shared |
|---|---|---|---|---|---|
| Python | 1,653 | 6,678 | 1,615 Perfect / 5 Good / **0 below 0.99** | 0.999988 | 18,898,405 / 18,898,406 |
| rpkg | 1,653 | 6,678 | 1,601 Perfect / 10 Good / 9 Poor / **0 below 0.95** | 0.999843 | 18,898,405 / 18,898,406 |

Both share 0 NA-pattern mismatches. rpkg's 19 non-Perfect columns are entirely
the pre-existing irreducible classes (FDC90th near-zero-tail OLS plus its
downstream Pettitt tie flips, and 3 recession Spearman p-values). Full detail,
including the characterised residuals, the defects the gates caught that unit
tests could not, and the convention decisions taken (Mann-Kendall
p-value method; sparse-family row emission; canonical metric names; R's
`mean()` vs Julia's sequential summation in the drought smoother), is in
[`CROSS_LANGUAGE_STATUS.md`](CROSS_LANGUAGE_STATUS.md) and CHANGELOG →
August 2026. Benchmark timings are no longer comparable to the April table
below (the ports now compute ~2.6x the columns).

**Acceptance gating changed too**: `compare_*_vs_golden_julia.py` intersects
columns and exits 0 on missing schemas, so they are DIAGNOSTICS. The gate is
`docs/benchmarks/check_schema_equality.py` (strict column + gage set equality
with explicit, logged waivers).

### Historical status (April 28, 2026 — Pettitt changepoint only)

*(April 2026 snapshot — superseded, retained for provenance.)* Julia produced 1,264 columns (656 base/metadata + 608 changepoint = 76 sigs × 8 Pettitt fields) across 7,313 gages. Python and rpkg were at 624 (changepoint port pending). R canonical still has the old recession algorithm (46 poor columns, sync pending).

| Metric | rpkg | Julia | Python | R (canonical) |
|--------|------|-------|--------|---------------|
| Total Time | 244.6 min* | ~15 min | 150.4 min | 874 min** |
| Gages Processed | 7,313 | 7,313 | 7,313 | 5,707 |
| Signature Columns | 624 | 1,264 | 624 | 551 |
| Processing Rate | 0.5/s* | ~8/s | 0.81/s | 0.11/s** |

*rpkg April 27 re-run had I/O contention (concurrent with Python benchmark). Previous solo run: ~114 min.
**R canonical March 16-17 re-run, also with I/O contention. Previous solo R runs: ~1-2 hours.

**Pettitt changepoint signal summary**: 13.4% overall significance rate (2.7x null expectation). Strongest signal in Flashiness (19.4%), Flow Percentiles (18.4%), Baseflow (15.1%); weakest in Flow Timing (3.7%, below null). Effective independence ~17/76 signatures (77% redundancy). After BH-FDR correction, ~3.5% of evaluations survive. See `docs/SIGNATURES.md` → Changepoint Detection → Signal Robustness for full analysis.

#### Synced Implementations (rpkg, Python, Julia — 624 columns, April 27, 2026)

| Pair | Perfect (>=0.999) | Good (0.99-0.999) | Poor (<0.99) | Min R² |
|------|-------------------|-------------------|-------------|--------|
| Python vs Julia | **615** | 5 | 3 | 0.984 |
| rpkg vs Julia | **614** | 4 | 5 | 0.969 |

Python-Julia: 620 of 623 columns (99.5%) have R² >= 0.99. All 25 new recession-parameterized BFI columns are Perfect.

rpkg-Julia: 618 of 623 columns (99.2%) have R² >= 0.99. 24 of 25 new recession-parameterized BFI columns are Perfect (1 poor: `alpha_linear_spearman_pval` — Spearman p-value library difference).

#### R Canonical vs Synced Implementations (April 14, 2026)

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 |
|------|---------|-----------|--------|-------------|
| R vs Python | 0.6018 | 0.9999 | -49.72 | 49 |
| R vs Julia | 0.6022 | 0.9999 | -49.72 | 49 |
| rpkg vs R | 0.6791 | 0.9999 | -50.04 | 46 |

All 46-49 poor columns are recession metrics — R canonical still has the old look-ahead algorithm. Non-recession categories show near-perfect agreement (0 poor columns).

#### Known Remaining Divergences (7 columns max, April 2026)

All poor columns across the 3 synced implementations are irreducible library-level differences:
- 2 recession pointcloud p-values: Spearman p-value calculation differences (exact permutation vs t-approximation for small n)
- 2 FDC90th p-values: 28-gage NA mismatch from floating-point precision in near-zero regression
- 3 recession seasonality minimum: Sinusoidal fit sensitivity (Python-Julia only, R²=0.984-0.989) — **RESOLVED August 2026**: the true cause was an off-by-one in Python's mid-event day-of-water-year (upper-middle index vs Julia's floor midpoint), not fit sensitivity; after the fix the seasonality outputs agree to ~1e-15

#### Julia Post-Section 3 vs Golden R (Feb 2026)

Julia's April 2026 output includes Guidelines Section 3 changes (new signatures, recession algorithm fix, trend_completeness) that the Feb 2026 Golden R reference predates. A separate comparison pipeline (`compare_julia_vs_golden_r.py`) uses 6-tier R² classification. See `docs/benchmarks/julia_vs_golden_r_summary.md` for the full report. Key divergence drivers:

| Root Cause | Cols Affected | Type |
|------------|--------------|------|
| trend_completeness (80% gate) | 220 | Temporal gap — Golden R predates feature |
| Recession algorithm rewrite | 46 | Intentional — R canonical sync pending |
| Elasticity operator bug (fixed) | 9 | Bug fix applied |
| dur_low_pulses_all NAs | 6 | Under investigation |

#### Alignment Progress (historical — Spearman rho through Round 6, then identity R²)

| Pair | Round 0 | Round 6 | Post-tau-b (Apr) | Post-Section 3 (Apr 15) |
|------|---------|---------|-----------------|------------------------|
| R vs Python | 323 | **4** | 4 | 49* |
| R vs Julia | 321 | **4** | 4 | 49* |
| Python vs Julia | 73 | **0** | 0 | **3** |
| rpkg vs Python | — | — | 6 | **7** |
| rpkg vs Julia | — | — | 8 | **4** |

*R canonical recession divergence (46 cols) is intentional — old algorithm, sync pending. Non-recession: 3 poor cols (FDC90th p-values + dur_low_pulses).

### Running Benchmarks

```bash
# All scripts use __file__/@__DIR__ for paths, so they work from any directory.
# Run from project root for consistency:

# R canonical benchmark (~1-5 hours)
Rscript docs/benchmarks/run_r_benchmark.R

# rpkg benchmark (~2 hours)
Rscript docs/benchmarks/run_rpkg_benchmark.R

# Python benchmark (~70-130 min)
python docs/benchmarks/run_python_benchmark.py

# Julia benchmark (~10 min)
julia docs/benchmarks/run_julia_benchmark.jl

# Primary comparisons: Python/rpkg vs Julia golden (canonical)
python docs/benchmarks/compare_python_vs_golden_julia.py
python docs/benchmarks/compare_rpkg_vs_golden_julia.py

# Julia golden dashboards
python docs/benchmarks/build_julia_golden_dashboard.py python
python docs/benchmarks/build_julia_golden_dashboard.py rpkg

# Legacy three-way comparison (historical — uses R as baseline)
python docs/benchmarks/compare_three_way.py
python docs/benchmarks/compare_rpkg.py

# Sensitivity experiments (~10-20 min each)
julia docs/benchmarks/run_julia_benchmark_startIn1993.jl
julia docs/benchmarks/run_julia_benchmark_startIn1993_60pct.jl
julia docs/benchmarks/run_julia_benchmark_startIn1993_80pct.jl

# Experiment comparisons + dashboards
python docs/benchmarks/compare_experiment_vs_julia.py startIn1993
python docs/benchmarks/compare_experiment_vs_julia.py startIn1993_60pct
python docs/benchmarks/compare_experiment_vs_julia.py startIn1993_80pct
python docs/benchmarks/build_experiment_vs_julia_dashboard.py startIn1993
python docs/benchmarks/build_experiment_vs_julia_dashboard.py startIn1993_60pct
python docs/benchmarks/build_experiment_vs_julia_dashboard.py startIn1993_80pct
```

### Benchmark Files

| File | Description |
|------|-------------|
| `docs/benchmarks/run_r_benchmark.R` | R canonical full signature extraction |
| `docs/benchmarks/run_rpkg_benchmark.R` | rpkg package full signature extraction |
| `docs/benchmarks/run_python_benchmark.py` | Python full signature extraction |
| `docs/benchmarks/run_julia_benchmark.jl` | **CANONICAL** — Julia full signature extraction (reference output) |
| `docs/benchmarks/compare_python_vs_golden_julia.py` | **PRIMARY** — Python vs Julia golden, 6-tier R² classification |
| `docs/benchmarks/compare_rpkg_vs_golden_julia.py` | **PRIMARY** — rpkg vs Julia golden, 6-tier R² classification |
| `docs/benchmarks/build_julia_golden_dashboard.py` | Interactive HTML dashboard: Python or rpkg vs Julia golden |
| `docs/benchmarks/validate_annual_values.py` | Annual values parquet vs summary CSV consistency check (run after each benchmark). **Self-referential** — necessary, not sufficient; pair it with the cross-language gate below |
| `docs/benchmarks/check_schema_equality.py` | **GATE** — strict column-set + gage-set equality vs a reference, with waivers that must be named on the command line. The `compare_*` scripts intersect columns and exit 0 on a missing family, so they are diagnostics, not gates |
| `docs/benchmarks/check_annual_parquet_equality.py` | **GATE** — cross-language annual parquet vs the Julia reference: signature set, key set, duplicates and NA patterns exactly; values to a stated tolerance. This is the *sufficient* annual check `validate_annual_values.py` cannot be |
| `docs/benchmarks/check_signature_failures.py` | **GATE** — scans a run log for the per-gage signature-family exceptions all three orchestrators swallow into warnings (they resurface as ordinary NA in the rectangular CSV, invisible to every column-based check). Refuses to certify an R log ending in the deferred-warning truncation banner |
| `docs/benchmarks/run_rpkg_acceptance_gates.sh` | Runs the four gates above against an rpkg benchmark in one command |
| `docs/benchmarks/compare_three_way.py` | Legacy — Three-way comparison (historical, uses R as baseline) |
| `docs/benchmarks/compare_rpkg.py` | Legacy — rpkg vs all other implementations |
| `docs/benchmarks/compare_julia_vs_golden_r.py` | Historical — Julia vs Golden R (Feb 2026) — 6-tier R² classification |
| `docs/benchmarks/build_comparison_dashboard.py` | Historical — Interactive HTML dashboard (maps + scatterplot) |
| `docs/benchmarks/build_julia_vs_golden_r_dashboard.py` | Historical — Julia vs Golden R with dual maps + scatterplot |
| `docs/benchmarks/build_new_vs_golden_julia_dashboard.py` | New benchmark vs golden Julia validation dashboard |
| `docs/benchmarks/build_section3_dashboard.py` | Section 3 pre/post comparison dashboard |
| `docs/benchmarks/compare_experiment_vs_julia.py` | Parameterized experiment vs Julia baseline comparison |
| `docs/benchmarks/build_experiment_vs_julia_dashboard.py` | Parameterized experiment dashboard (HTML) |
| `docs/benchmarks/run_julia_benchmark_startIn1993.jl` | Experiment: WY >= 1993 wrapper |
| `docs/benchmarks/run_julia_benchmark_startIn1993_60pct.jl` | Experiment: WY >= 1993 + 60% qualifying fraction wrapper |
| `docs/benchmarks/run_julia_benchmark_startIn1993_80pct.jl` | Experiment: WY >= 1993 + 80% qualifying fraction wrapper |
| `docs/benchmarks/run_julia_benchmark_drought_1993_2025_60pct.jl` | **STANDARD PRODUCT #1** wrapper — WY 1993-2025 @ 60% with the drought family (promoted 2026-08-10) |
| `docs/benchmarks/run_julia_benchmark_prod_1980_2025_60pct_drought.jl` | **STANDARD PRODUCT #2** wrapper — WY 1980-2025 @ 60% with the drought family (run 2026-08-11; defaults to the rebuilt climate parquet) |
| `docs/benchmarks/check_additivity.jl` | Proves a run ADDED columns without changing pre-existing ones (column/gage set, per-gage value identity, population gate) |
| `docs/benchmarks/analyze_drought_redundancy.jl` | Measures drought-duration overlap with the pulse metrics on the annual series |
| `docs/benchmarks/convert_daymet_csvs_to_parquet.py` | Rebuilds `daymet_1980_2023.parquet` from the 44 annual Daymet CSVs (Python equivalent of the R `convert_daymet_zip_to_parquet()`) |
| `docs/benchmarks/comparison_report.md` | Generated comparison report |
| `docs/benchmarks/julia_vs_golden_r_summary.md` | Generated Julia vs Golden R detailed report |

For implementation details, alignment history, and known divergences, see [`CROSS_LANGUAGE_STATUS.md`](CROSS_LANGUAGE_STATUS.md).

## Data Sources

### USGS (via dataRetrieval)
- Parameter code: `00060` (Discharge, cubic feet per second)
- Quality codes accepted: `A`, `A e`, `P`, `P e`
- Conversion: cfs -> mm/day using drainage area

### Canadian HYDAT (via tidyhydat)
- Parameter: Flow (m3/s)
- Includes both regulated and unregulated stations (REGULATED flag tracked in metadata)
- Conversion: m3/s -> mm/day using drainage area (`DRAINAGE_AREA_GROSS` from the HYDAT STATIONS table)
- Interference metadata (RHBN, REGULATED) exported to `metadata/canadian_hydat_interference.csv` for cross-language use

#### Missing drainage areas (`area_normalized = FALSE`)

HYDAT publishes **no drainage area at all** (neither `DRAINAGE_AREA_GROSS` nor
`DRAINAGE_AREA_EFFECT`) for 1,601 Canadian stations in the Feb 2026 metadata — of
which **73 are `processing_status == success` gages**, and **37 survive the 20-year
filter into the Julia canonical signature output**. This was verified directly against
`Hydat.sqlite3` (July 2026): the ingestion code is correct; the source data simply
does not exist.

Station names explain why. Of the 73:
- **~40 irrigation/diversion canals, ditches, and drains** (Alberta/Saskatchewan
  irrigation districts, Welland Canal diversion) — a canal has no drainage basin, so
  watershed area is genuinely undefined
- **~15 dam/powerhouse outflows** (Revelstoke, Arrow Reservoir, Duncan Dam, Kemano,
  Jenpeg, Ghost Tailrace, Skins Lake Spillway) — a real upstream area exists but WSC
  does not publish it
- **Several huge-river channel splits / lake outlets** (St. Lawrence at LaSalle,
  Mackenzie at Strong Point + East Channel at Inuvik, both Nelson River channels, the
  three Lake of the Woods outlets) — each channel carries only part of the flow, so
  the full upstream area would mis-normalize
- **~8 apparently natural streams** (e.g., Ewart Creek 08NL076, Dash Creek 08MD035,
  Bridge Creek 08LA027, Cedar Creek 08MH166, Blackstone River 10ED007)

62 of the 73 are `REGULATED = TRUE`.

**DECISION (user, July 2026)**: keep these gages in the product with **raw m³/s
flow — no area backfill** (see feasibility assessment below for why backfill was
rejected). Instead, all **Q-to-PPT signatures are structurally gated**:
`calculate_all_signatures()` in Julia, Python, and rpkg takes an
`area_normalized = true` argument; when false, runoff ratios (+
`runoff_ratio_high_count`), elasticity (+ diagnostics), Q-P seasonality, and
`avg_storage` are skipped entirely — Q (m³/s) and PPT (mm) units don't match, so
Q/P, dQ/dP, and cumsum(P − Q) are meaningless. Q-only signatures (Q##/D##
percentiles, recession, BFI, FDC, timing, pulses, flashiness) still run. The three
benchmark runners read `area_normalized` from the metadata CSV (leading-zero-safe
join; only an explicit FALSE gates) and pass it per gage. Tests:
`julia/test/test_area_normalized_gate.jl`, `python/tests/test_area_normalized_gate.py`,
`rpkg/tests/testthat/test-area_normalized_gate.R`.

**Downstream consequence — mixed units in the signature CSV**: unit-carrying Q-only
signatures (Qann/seasonal totals, flow percentiles, Q95_Q10, log_a) remain in raw
m³/s for these 37 rows and are incomparable with the mm/day-based gages. Q-to-PPT
signatures are NA by design (as of July 2026 — before the gate, 16 of the 37 carried
mixed-unit runoff ratios/elasticity/storage values; those become NA at the next
benchmark). `flagged_for_qann_range` catches only 27 of the 37 (the other 10 have
m³/s totals that happen to land inside [0, 2000]); the `area_normalized` column is
the only reliable discriminator. **Downstream users should filter on
`area_normalized == TRUE` before any cross-gage comparison of unit-carrying
signatures.** Tracked as a Known Issue in [CHANGELOG.md](../CHANGELOG.md).

**HydroBasins backfill feasibility (assessed July 2026)**: `UP_AREA` from the lev12
outlet basin (`basinAt_NorAm_polys.gpkg`, lookup via the metadata's
`Downstream_HB_ID`, available for 41/73) was validated against the 1,383 Canadian
success gages that have BOTH a published HYDAT area and a snapped outlet basin.
Accuracy is strongly size-dependent: median |rel. error| is 1.3% for basins >10,000
km², 2.4% for 2,000–10,000, 6.4% for 500–2,000, 22% for 100–500, and **549% for <100
km²** (lev12 polygon granularity ≈ small-basin size → unusable). Ground-truthing the
missing gages against neighboring published stations confirms the pattern: main-stem
dam/reservoir outflows are near-exact (Revelstoke 26,250 vs 26,700 published; Arrow
36,471 vs 36,500), while channel splits are structurally wrong regardless of snapping
(Nelson East Channel → 927,013 km², West Channel → 2,957 km²; each carries only part
of the flow — WSC leaves these blank deliberately) and canals/diversions get the area
of whatever basin they happen to cross (physically meaningless). HydroATLAS adds
nothing over HydroBasins here (same lev12 framework, same `UP_AREA`). **Verdict:
backfill is defensible ONLY for main-stem dam-outflow stations (name-reviewed,
~10–15 gages); small natural streams should instead take areas from the official ECCC
MDA_ADP polygons already delineated in the EO geometry layer
(`watershed_polygons_26jun2026.parquet`, `geom_area_km2`); canals, diversions, and
channel splits should remain un-normalized (or be excluded from the signature
product).** Outcome: the user opted to skip backfill entirely and gate the Q-to-PPT
signatures instead (see DECISION above).

### Caravan
- NetCDF format with daily streamflow + climate variables
- Includes: PPT, SWE, temperature
- Trade-off: Shorter records (ends ~2018-2020) but has climate data

## Parquet Data Files

### Active Parquet Files (Use These)

| File | Location | Created | Description |
|------|----------|---------|-------------|
| `combined_streamflow_data_09feb2026.parquet` | `D:/processedOuts_feb2026/` | Feb 2026 | Current streamflow parquet with bug fixes. 111,624,189 rows / 8,014 gages, **858 MiB** (899,935,349 B) |
| `combined_watershed_metadata_09feb2026.csv` | `D:/processedOuts_feb2026/` | Feb 2026 | Corresponding metadata, 2.1 MiB |
| **`daymet_1980_2023_rebuilt_10aug2026.parquet`** | `D:/processedOuts_feb2026/` | Aug 2026 | **USE THIS** for climate (PPT, temp, SWE). 97,757,220 rows / 6,087 sites / WY-covering CY 1980–2023, **3.76 GiB** (4,040,997,608 B). Rebuilt from the annual CSVs after the original was truncated — see below |
| `daymet_1980_2023/` (44 annual CSVs) | `D:/` (drive root) | Jan 2026 | **Source of record** for Daymet: `site_id, month, year, prcp, tmin, tmax, swe, vp, srad`; no day column. 10,825,714,689 B = 10.08 GiB total, ~234.6 MiB/yr. Rebuild the parquet with `docs/benchmarks/convert_daymet_csvs_to_parquet.py` |

> ⚠️ **`daymet_1980_2023.parquet` (the canonical name) is TRUNCATED and unreadable** —
> 1,261,436,928 bytes with no `PAR1` footer, against 4,125,630,653 recorded in the 28 Jul
> run's provenance block (mtime unchanged, so only the byte count reveals it). It is left
> in place; **set `STREAMFLOW_CLIMATE_PATH` to the `_rebuilt_10aug2026` file** (the
> WY 1980–2025 wrapper already defaults to it). Standard product #1 was built from the
> ORIGINAL, #2 from the rebuild; a controlled replay bounds the difference at ≤ 3.4e-13.
> See CHANGELOG → August 2026.

> ⚠️ **Provenance limitation on both delivered products.** Each run's `timing.json` records
> `git_working_tree_dirty = true` (#1 at `b7e8988`, #2 at `0487bbd`) and no patch or
> source-tree hash was retained, so neither product can be tied byte-for-byte to a
> reproducible source tree. The multi-GB inputs carry size/mtime but no hash unless
> `STREAMFLOW_HASH_INPUTS=1`. **For future standard products: run from a clean tree, or
> retain the diff, and set `STREAMFLOW_HASH_INPUTS=1`.**

> **To be unambiguous about product #1 and the truncated file — #1 is NOT built on
> truncated data.** It ran 2026-07-28 against the intact 4,125,630,653-byte parquet (its
> provenance records that size), and a truncated parquet cannot be read at all: the footer
> is at the end of the file, so the reader fails to open it outright — which is exactly how
> the corruption surfaced on 2026-08-10. The corruption occurred **after** the run, with the
> mtime unchanged (consistent with exFAT damage, not a rewrite). Confirmed independently:
> replaying #1's config against the rebuilt, verified-complete climate input reproduced it
> with 0 columns added/dropped, an identical gage set, and ≤ 3.4e-13 — impossible had #1
> been computed from a partial climate record. The consequence is confined to
> **reproducibility**: #1 cannot be regenerated byte-for-byte (nor its 28 Jul additivity
> PASS recreated), because those exact input bytes no longer exist. The delivered product
> itself is sound.

**Boundary polygons** (source layers, on the same drive):

| Layer | Location | Size |
|---|---|---|
| HydroBASINS lev12 polygons (167,665 N. American basins) | `geospatial_derivedData/basinAt_NorAm_polys.gpkg` | 462 MiB (484,356,096 B) |
| HydroBASINS lev12 centroids (*not* boundaries) | `geospatial_derivedData/basinAt_NorAm_centroids.gpkg` | 131 MiB (136,949,760 B) |
| Official US basin boundaries (GAGES-II) | `official_watershed_polygons/gagesII_bnd/` | **120 MiB extracted**, plus an 18 MiB ZIP of the same data (138 MiB on disk) |
| Official Canadian basin boundaries (ECCC/WSC MDA_ADP, 11 gpkg) | `official_watershed_polygons/ca/` | 4.75 MiB apparent (5.0 MiB allocated) |
| BasinATLAS v10 geodatabase (281 attributes; HydroATLAS metadata source, not boundaries) | `BasinATLAS_Data_v10/` | 6.1 GiB |
| ⚠️ `unified_watersheds_simplified.gpkg` — **corrupt CSV, not a polygon layer** | `geospatial_derivedData/` | 2.0 MiB (do not use) |

Boundary layers proper (excluding centroids, BasinATLAS attributes, and the redundant ZIP)
total **≈ 587 MiB**. Folder figures are `du` (exFAT-allocated) and run slightly above the
summed file bytes — e.g. `gagesII_bnd` is 124.7 MiB apparent vs 138 MiB allocated.

The DERIVED per-gage geometry product (`watershed_polygons_26jun2026.{gpkg,parquet}`,
~42/36 MB, 7,964 basins) was delivered to S3, not this drive — but **S3 access was lost
on 2026-08-24; retrieve it from the project Google Drive backup instead** (CHANGELOG →
August 2026) — see `EO_data_processing/README.md`.

### Deprecated Parquet Files (DO NOT USE)

| File | Location | Issue |
|------|----------|-------|
| `combined_streamflow_data.parquet` | `D:/combined_streamflow_output/` | **CORRUPTED** - Contains 99999 multiplier bug for Canadian gages without basin area |
| `combined_watershed_metadata.csv` | `D:/combined_streamflow_output/` | Outdated metadata |
| `streamflowSignature_summaryData_OCT2025.csv` | `D:/combined_streamflow_output/` | Generated from corrupted parquet |

### The 99999 Bug

The October 2025 parquet was created with buggy code that applied `conversion = 99999` for Canadian gages without basin area, instead of keeping values in raw units. This resulted in Q values ~100,000x too high.

**Example (gage 08ND025):**
- Raw HYDAT: 785.6 m³/s
- Corrupted parquet: 78,557,297 "mm/day" (785.6 × 99999)
- Fixed parquet: 785.6 (raw m³/s, flagged as `area_normalized = FALSE`)

See docs/CHANGELOG_ARCHIVE.md entry for "H5 follow-up" for full details of the fix.

## Human Interference Metadata

Watershed metadata is automatically enriched with human interference indicators when `concatenate_with_metadata()` is called during data processing.

### Data Sources

**USGS Gages (GAGES-II):**
- Location: `D:/gagesMetadata/` (configured via `GAGES_II_DIR` in config.R)
- Files: `conterm_hydromod_dams.txt`, `conterm_bas_classif.txt`, etc.
- Columns extracted: NDAMS_2009, MAJ_DDENS_2009, STOR_NID_2009, IMPNLCD06, DEVNLCD06, FRESHW_WITHDRAWAL, HYDRO_DISTURB_INDX, CLASS

**Canadian Gages (HYDAT via tidyhydat):**
- RHBN: Reference Hydrometric Basin Network designation
- REGULATED: Station regulation status from `hy_stn_regulation()`

### Unified Classification

The `human_interference_class` column provides a unified classification:
- **reference**: USGS gages with CLASS="Ref" or Canadian gages with RHBN=TRUE
- **non-reference**: USGS gages with CLASS="Non-ref" or Canadian gages with RHBN=FALSE
- **unknown**: Gages without classification data

### Manual Enrichment

To re-enrich existing metadata (one-time use):

```r
source("config.R")
source("R/helperFunctions.R")
source("R/run_enrich_metadata.R")
```

## Dependencies

```r
# Core data handling
library(data.table)
library(arrow)        # Parquet I/O
library(lubridate)    # Date handling

# Data retrieval
library(dataRetrieval) # USGS NWIS API
library(tidyhydat)     # Canadian HYDAT database

# Statistics
library(zyp)          # Theil-Sen slope estimation
library(Kendall)      # Mann-Kendall trend test
library(mblm)         # Alternative Theil-Sen

# Spatial (for basin delineation)
library(sf)
library(terra)

# Visualization app
library(shiny)
library(leaflet)
library(plotly)
library(aws.s3)       # S3 data storage
```

# HISSS — Hydrologic Information Signatures and Summary Statistics

A framework for extracting hydrological signatures from daily streamflow data, with implementations in **Julia (canonical), Python, and R** that produce validated, near-identical results. It processes data from USGS NWIS, Canadian HYDAT, and Caravan sources to calculate 100+ signatures characterizing flow regimes, trends, changepoints, and variability — including snow and streamflow-drought metrics — for thousands of gages across the United States and Canada.

This repository contains the code behind the **HISSS dataset** (Kaiser et al., manuscript in preparation for *Scientific Data*). The delivered data products live on HydroShare: [HISSS collection](https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/).

## Project Goals

**Primary:**
1. **Data Processing** — Ingest raw streamflow data (USGS, HYDAT, Caravan), clean and filter by configurable quality thresholds, collate metadata, and produce standardized parquet/CSV outputs.
2. **Signature Extraction** — Calculate 100+ hydrological signatures under strict guardrails. Methodology is defined by domain experts in a [plain-English guidelines document](https://docs.google.com/document/d/e/2PACX-1vQnt7OCPm19vnWF4yynXL9JTzTvq9CrGoEaDv7yFSngLoFsypiWsx6fZLKWwaO5YQ/pub) (synced to `docs/SIGNATURE_GUIDELINES.md`) — hydrologists define what to calculate, and code implements those definitions.

**Secondary:**
3. **Visualization** — Interactive exploration of signatures, trends, and cross-signature relationships across thousands of gages (self-contained HTML explorers; dashboard rebuild planned).
4. **Cross-Language Implementations** — Python and R ports validated at full scale against the Julia canonical implementation (August 2026: identical 1,653-column output for the same 6,678 gages; see [Cross-Language Alignment](#cross-language-alignment)).

## Features

- **Multi-source data ingestion**: USGS NWIS, Canadian HYDAT, Caravan/HYSETS datasets
- **Comprehensive signature suite**: flow volumes and percentiles, flow timing, flow duration curve slopes, baseflow indices (fixed and recession-parameterized), recession parameters, pulse metrics, flashiness, runoff ratios, elasticity, Q-P seasonality, storage, snow metrics (Daymet SWE), and streamflow drought
- **Robust trend analysis**: Theil-Sen and linear slopes, Mann-Kendall and Spearman tests for every signature's annual series
- **Changepoint detection**: non-parametric Pettitt test with pre/post segment diagnostics for every signature
- **Annual values export**: the per-year series behind every statistic, as a long-format parquet
- **Centralized quality control**: per-water-year qualification via a single preprocessor (`preprocess_daily_data()`), configurable via `config/signatures_config.json`

## Quick Start

Each implementation has its own README with full API details, input format, and climate-dependent signatures: [`julia/README.md`](julia/README.md) · [`python/README.md`](python/README.md) · [`rpkg/README.md`](rpkg/README.md).

Input expectations: a table with columns `gage_id`, `Date` (or `date`), and `Q` (mm/day; see [Units](#critical-conventions)); optional climate columns (`PPT`, `SWE`) enable the climate/snow families.

### Julia (canonical)

```julia
# julia> using Pkg; Pkg.activate("julia"); Pkg.instantiate()
using StreamflowSignatures, DataFrames

df = read_parquet("path/to/streamflow.parquet")
df = add_water_year_columns(df)          # auto-detects "Date" or "date"
gage_data = df[df.gage_id .== "01011000", :]

results = calculate_all_signatures(gage_data, false)   # has_climate = false
```

### Python

```python
# pip install -e python/
import pandas as pd
from streamflow_signatures import add_water_year_columns, calculate_all_signatures

df = pd.read_parquet("path/to/streamflow.parquet")
df = add_water_year_columns(df)
gage_data = df[df["gage_id"] == "01011000"]

results = calculate_all_signatures(gage_data, has_climate=False)
```

### R (rpkg)

```r
# R CMD INSTALL rpkg
library(streamflowsignatures)

df <- arrow::read_parquet("path/to/streamflow.parquet")
df <- add_water_year_columns(df)
gage_data <- df[df$gage_id == "01011000", ]

results <- calculate_all_signatures(gage_data, has_climate = FALSE)
```

**Producing the full product**: the bare `calculate_all_signatures` call returns the 8 statistics per signature. The full 1,653-column product additionally requires the `changepoint` configuration (Pettitt fields), an explicit `snow_data` frame (snow metrics), climate data (`has_climate = true` + PPT), and the annual-values collector. The benchmark runners wire all of this up and are the reference entry points:

```bash
julia --project=julia docs/benchmarks/run_julia_benchmark.jl   # canonical full pipeline
python docs/benchmarks/run_python_benchmark.py
Rscript docs/benchmarks/run_rpkg_benchmark.R
```

**Year qualification** is handled centrally by `preprocess_daily_data()` (the production path; `use_legacy_filtering: false`): internal gaps ≤ 3 days are interpolated, years with > 30 raw NAs or > 3-day gaps are rejected, and gages need 20+ qualifying water years. The per-function `filter_qualifying_years()` shown in some examples is the legacy path only.

## Signature Categories

| Category | Metrics | Climate needed | Description |
|----------|---------|----------------|-------------|
| **Flow Volumes** | Qann, Qwin, Qspr, Qsum, Qfal, Q1–Q99 (15 percentiles), Q95_Q10 | No | Annual/seasonal totals and percentiles (21 metrics) |
| **FDC** | FDCall, FDC90th, FDCmid | No | Flow duration curve slopes |
| **Baseflow** | BFI_Eckhardt, BFI_LyneHollick (+ recession-parameterized variants) | No | Groundwater contribution indices |
| **Recession** | log_a, b, concavity, alpha_linear, n_recession_events, seasonality | No | Recession parameters (alpha under the b=1 linear-reservoir convention) |
| **Pulse Metrics** | n/dur high & low pulses (per-year and period-of-record), TQmean, flow reversals | No | High/low flow event characteristics |
| **Flashiness** | flashinessRB | No | Richards-Baker flashiness index |
| **Flow Timing** | D1–D99_day (13 percentiles), D25_to_D75, Dmax | No | Cumulative flow timing (15 metrics) |
| **Negative Flow Days** | negative_ann | No | Days with Q < 0 per year |
| **Runoff Ratios** | annual + 4 seasonal | Yes | Q/P ratios |
| **Elasticity** | elasticity_rolling, elasticity_annual, elasticity_static | Yes | Streamflow sensitivity to precipitation |
| **Q-P Seasonality** | qp_slope_sd, qp_bimodality | Yes | Precipitation-runoff relationship |
| **Average Storage** | avg_storage | Yes | Mean catchment storage (computed, but omitted from major analyses — no ET term) |
| **Snow Metrics** | swe_max, snow on/off timing, melt rate/timing, SSM, … (14 metrics) | SWE (Daymet) | Snowpack magnitude, timing, and regime |
| **Streamflow Drought** | duration + deficit at 5 fixed severity levels (+ 5 threshold scalars) | No | Days below and departure from fixed percentile thresholds (USDM-analog ladder) |

Every time-series signature produces **16 statistics**: 8 trend/summary statistics — `_senn_slp` (Theil-Sen), `_linear_slp`, `_spearman_rho`, `_spearman_pval`, `_mk_rho` (Mann-Kendall tau), `_mk_pval`, `_mean`, `_median` — plus 8 Pettitt changepoint fields (`_pettitt_cp_year`, `_pettitt_pval`, pre/post means, delta, percent change, pre/post MK p-values). A 20-value stats floor and 60%/80% trend-completeness gates govern when statistics are emitted; see `docs/SIGNATURES.md`.

## Cross-Language Alignment

Julia is canonical; all changes start there. Python and rpkg were validated at full scale against the Julia reference output (WY 1993–2025 standard configuration, 6,678 gages × 1,653 columns) in August 2026:

| | Columns | Shared signature columns with R² ≥ 0.999 ("Perfect") | Below R² 0.99 | Mean R² | Validated |
|---|---|---|---|---|---|
| **Python** | 1,653 | 1,615 / 1,620 (99.7%) | 0 | 0.999988 | 2026-08-26 |
| **rpkg** | 1,653 | 1,601 / 1,620 (98.8%) | 9 | 0.999843 | 2026-08-27 |

Both ports pass strict schema equality (identical column and gage sets), a swallowed-failure log gate, and cross-language annual-parquet equality (18,898,405 of 18,898,406 rows shared, 0 NA-pattern mismatches). The few remaining sub-0.99 columns are characterized, irreducible library-level differences (near-zero-tail OLS in FDC90th and downstream Pettitt tie flips; Spearman p-value methods at small n). All implementations share configuration via `config/signatures_config.json`.

See [`docs/CROSS_LANGUAGE_STATUS.md`](docs/CROSS_LANGUAGE_STATUS.md) for the full methodology, gate definitions, and alignment history.

## Output Format

A full run writes **two** files to its output folder:

### 1. Summary CSV — one row per gage (1,653 columns)

| Column Type | Examples | Description |
|-------------|----------|-------------|
| Metadata | gage_id, latitude, longitude, basin_area, area_normalized | Gage identification |
| Record Info | num_water_years, start_water_year, end_water_year | Data coverage |
| Human Interference | NDAMS_2009, HYDRO_DISTURB_INDX, CLASS, RHBN, REGULATED, human_interference_class | Watershed disturbance indicators |
| Signatures | Qann_senn_slp, …, Qann_pettitt_cp_year, … | 16 statistics per signature base + 21 per-gage scalars |
| QA/QC flags | flagged_for_qann_range, flagged_for_high_na, … | 12 automated quality flags |

> **Un-normalized gages**: 37 Canadian gages have no published drainage area and carry flow in raw m³/s (`area_normalized = FALSE`); their Q-to-PPT signatures are NA by design. **Filter on `area_normalized == TRUE` before any cross-gage comparison of unit-carrying signatures.**

### 2. Annual values parquet — `{prefix}_signatures_annual.parquet`

The per-year value of every signature, **before** aggregation into trends — long format (all three languages):

| Column | Type | Notes |
|---|---|---|
| `gage_id` | String | zero-padded; joins the summary CSV |
| `signature` | String | matches the summary's column prefixes |
| `water_year` | Int32 | Oct 1 – Sep 30 |
| `value` | Float64 | NaN and absent-row both mean "not computable that year" |

Covers all 100 time-series signatures (per-gage scalars have no annual series and correctly do not appear). Controlled by `annual_values.save` in `config/signatures_config.json` (on by default) and cross-validated against the summary CSV by `docs/benchmarks/validate_annual_values.py`.

> ⚠️ **Do not re-aggregate record-dependent signatures over a different window.** `drought_*`, the `*_all` pulse metrics, elasticity, and the recession-parameterized BFIs are computed against thresholds or record means from the ORIGINAL run window. Re-run the pipeline for a new window instead. Within-year-computable signatures (flow volumes/percentiles, timing, BFI, FDC, flashiness, runoff ratios, snow) re-aggregate safely.

## Standard Data Products

Two delivered products (both 1,653 columns, 60% qualifying fraction), published via the [HISSS HydroShare collection](https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/):

| Window | Gages | Annual parquet rows |
|---|---|---|
| WY 1993–2025 | 6,678 | 18,898,406 |
| WY 1980–2025 | 6,250 | 24,366,487 |

Neither is a subset of the other, and record-dependent signatures must never be compared across them (thresholds come from each run's own window).

## Human Interference Metadata

Watershed metadata is automatically enriched with human interference indicators:

**USGS gages (GAGES-II)**: NDAMS_2009, MAJ_DDENS_2009, STOR_NID_2009, IMPNLCD06, DEVNLCD06, FRESHW_WITHDRAWAL, HYDRO_DISTURB_INDX, CLASS (Ref/Non-ref).
**Canadian gages (HYDAT)**: RHBN (Reference Hydrometric Basin Network), REGULATED.
**Unified**: `human_interference_class` ∈ {reference, non-reference, unknown}.

## File Structure

```
├── julia/                         # Julia implementation (CANONICAL): 20 modules + package entry
├── python/                        # Python implementation (validated port)
├── rpkg/                          # R package (validated port)
├── config/signatures_config.json  # Shared cross-language config (single source of truth)
├── docs/                          # Extended documentation
│   ├── DEVELOPMENT.md             #   architecture, workflows, common tasks
│   ├── SIGNATURES.md              #   detailed signature documentation (14 categories)
│   ├── SIGNATURE_GUIDELINES.md    #   domain-expert guidelines (auto-synced)
│   ├── CROSS_LANGUAGE_STATUS.md   #   cross-language validation detail
│   └── benchmarks/                #   benchmark runners, comparison tools, validation gates
├── golden-outputs/                # Reference outputs for cross-language validation
├── metadata/                      # Basin/gage metadata (GAGES-II tables, HYDAT interference)
├── EO_data_processing/            # Earth Observation companions (MODIS LAI/LULC, Annual NLCD)
├── config.R                       # R-side configuration (ingestion scripts)
├── run_ingest_usgs_hydat.R        # Raw USGS/HYDAT data ingestion (R)
├── run_caravan_processing.R       # Caravan processing (R)
├── run_hydroatlas_metadata.R      # Per-gage HydroATLAS watershed metadata (R)
└── R/                             # Legacy R (deprecated for signatures; ingestion utilities active)
```

## Data Sources

| Source | Access | Use |
|---|---|---|
| USGS NWIS | `dataRetrieval::readNWISdv()`, parameter 00060 | Daily discharge, US gages |
| Canadian HYDAT | `tidyhydat::hy_daily()` (local HYDAT database) | Daily discharge, Canadian gages |
| Daymet V4 | ORNL DAAC; basin-aggregated | Precipitation, temperature, SWE (1 km) |
| Caravan | NetCDF bundles | Alternative pipeline: streamflow + climate |
| GAGES-II / WSC | USGS ScienceBase; ECCC | Watershed boundaries, disturbance metadata |
| HydroATLAS v10 | BasinATLAS | Static watershed attributes (~211 columns) |

See `docs/DATA_SOURCES.md` for the complete inventory.

## Configuration

Cross-language configuration lives in **`config/signatures_config.json`** — the single source of truth for NA handling, trend-completeness gates (60% overall, 80% first/last decade), the 20-value stats floor, changepoint window, snow and drought parameters, and the annual-values export. All three implementations read it at load.

Key conventions:

| Convention | Value |
|---|---|
| Water year | October 1 – September 30 |
| Flow units | mm/day (converted from cfs / m³/s via published drainage area) |
| Minimum record | 20+ qualifying water years per gage |
| Year rejection | > 30 raw NAs or > 3-day gap in a water year |
| Interpolation | internal gaps ≤ 3 days, linear |

(The legacy R path's `min_frac_good_data = 0.95` filter applies only when `use_legacy_filtering: true`.)

## Development

- `docs/DEVELOPMENT.md` — architecture, workflows, and common tasks
- `docs/SIGNATURES.md` — detailed signature documentation
- `CHANGELOG.md` — change history and known issues
- `docs/CROSS_LANGUAGE_STATUS.md` — cross-language alignment detail
- Per-language quickstarts: `julia/README.md`, `python/README.md`, `rpkg/README.md`

### Adding New Signatures

1. Implement in the appropriate `julia/src/` module (canonical), returning annual values, and produce statistics via `generate_stats()`
2. Register in `julia/src/signatures.jl` orchestration
3. Register in the validation surfaces: `EXPECTED_SIGNATURE_BASES` (config.R), `EXPECTED_DENSE_SIGNATURES` (`julia/test/test_annual_collector.jl`), the signature-count gate in `docs/benchmarks/validate_production_run.py`, and the documented column totals
4. Run the Julia unit suite, then the full benchmark; **prove additivity** against the previous canonical run (`docs/benchmarks/check_additivity.jl`)
5. Port to Python (`python/streamflow_signatures/`) and rpkg (`rpkg/R/`), then run the cross-language gates

See `docs/DEVELOPMENT.md` → "Adding a New Signature" for the full checklist.

## References (methodology)

- Eckhardt (2005); Lyne & Hollick (1979); Collischonn & Fan (2013) — baseflow filters
- Sen (1968); Mann (1945) & Kendall — trend estimation; Pettitt (1979) — changepoint test
- Baker et al. (2004) — Richards-Baker flashiness index
- Sawicz et al. (2011) — elasticity; Yilmaz et al. (2008) — FDC
- Hatchett (2021); Petersky & Harpold (2018) — snow seasonality
- Adelsperger et al. (in review) — streamflow drought severity ladder

Full citations per signature: `docs/SIGNATURES.md`.

## License

Mozilla Public License 2.0 — see [LICENSE](LICENSE).

## Contact

For questions or collaboration, please open an issue in the repository.

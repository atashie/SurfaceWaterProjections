# Streamflow Signatures Assistant

A skill for helping users understand, use, and interpret the streamflow signatures libraries.

## Overview

This skill assists scientists working with hydrological signature extraction from daily streamflow time series. It covers the Julia (canonical), Python, and R implementations.

## What Are Streamflow Signatures?

Streamflow signatures are quantitative metrics that characterize watershed hydrology. Each signature captures a different aspect of streamflow behavior:

- **Flow Volumes**: How much water flows annually and seasonally
- **Baseflow**: Groundwater contribution to streamflow
- **Recession**: How quickly flow decreases after precipitation
- **Flashiness**: How rapidly flow rises and falls
- **Timing**: When during the year flow occurs
- **Elasticity**: How streamflow responds to precipitation changes

## The 8-Statistic Rule

Every signature metric produces **8 statistics** via `generate_stats()`:

| Suffix | Statistic | Interpretation |
|--------|-----------|----------------|
| `_senn_slp` | Theil-Sen slope | Robust trend (units/year) |
| `_linear_slp` | Linear slope | Parametric trend (units/year) |
| `_spearman_rho` | Spearman rho | Correlation with time (-1 to 1) |
| `_spearman_pval` | Spearman p-value | Significance of correlation |
| `_mk_rho` | Mann-Kendall tau | Monotonic trend (-1 to 1) |
| `_mk_pval` | Mann-Kendall p-value | Significance of trend |
| `_mean` | Arithmetic mean | Central tendency |
| `_median` | Median | Robust central tendency |

### Recession Alpha Convention (July 2026)

`log_a_pointcloud`, `log_a_events`, and the `log_a_seasonality_*` scalars assume a
**linear reservoir (b = 1)**: `log(a) = median(log(-dQ/dt) - log(Q))`. They are NOT
free-fit intercepts — free-fit alpha is convolved with b, so alpha is decoupled by
fixing b = 1 everywhere. `b_pointcloud`/`b_events`/`concavity` remain free fits.
Interpretation: per-year `log_a_pointcloud` approximates `log(1 - alpha_linear)` —
approximate, not exact, because the two draw on slightly different event pools
(`alpha_linear` includes events whose free power-law fit failed).

### Per-Year Annual Values (July 2026)

The annual series behind these statistics are exported alongside the summary CSV
as a long-format parquet (`{prefix}_signatures_annual.parquet`): one row per
`(gage_id, signature, water_year)` with the annual `value`. Use it for per-year
QA, custom trend windows, or time-series plots without re-running extraction.
Interpretation notes:
- NaN value and absent row both mean "year not computable" — treat them the same.
- `elasticity_rolling` rows are keyed to the END year of each 11-year window;
  `elasticity_annual` to the later year of each consecutive-year pair.
- Recomputing `mean`/`median` from this file reproduces the summary `_mean`/`_median`
  columns exactly (validated by `docs/benchmarks/validate_annual_values.py`).
- Julia only (controlled by `annual_values.save` in `config/signatures_config.json`).

### Snow Metrics (July 2026)

Fourteen per-water-year metrics from Daymet daily SWE (`calculate_snow_metrics`,
Julia only): `swe_max`, `swe_max_dowy`, `snow_cover_days`, `snow_on_dowy`,
`snow_off_dowy`, `melt_season_days`, `melt_rate`, `ssm`, `swe_apr1`,
`melt_before_peak(+_pct, _to_max_swe)`, `melt_com_dowy`, `swe_max_to_ppt`.
Interpretation keys:
- ALL metrics use SWE ≥ 10 mm as the snow-day definition — including magnitudes
  (a year peaking at 8 mm reads as snow-free: swe_max = 0, timing NaN).
- Snow-on/off bracket the spell containing the annual peak; NaN means censored
  (snowpack predates Oct 1 / persists past Sep 30), NOT missing data.
- `ssm` ∈ [−1, +1]: +1 fully seasonal (≥ 60-day spells), −1 fully ephemeral
  (Hatchett 2021, Hydrology 8(1):32).
- Magnitude metrics emit valid ZEROS at snow-free gages; timing/melt metrics emit
  NaN there. Daymet covers Canadian gages too (July 2026 correction — ~1,100
  Canadian gages carry snow values); only gages without Daymet rows are all-NA.
- Trend/stat NaNs with values visible in the annual parquet usually mean the
  20-value stats floor (July 2026): metrics with <20 non-NA annual values report
  NaN for ALL 8 statistics (recession/elasticity exempt).
- Daymet SWE is modeled, not observed — prefer timing/trend signals over absolute
  magnitudes, especially in mountain terrain.

### Streamflow Drought Metrics (July 2026)

Ten per-water-year metrics (`calculate_drought_metrics`, Julia only) after
Adelsperger et al. (in review): `drought_duration_fixed_p{2,5,10,20,30}` (days below
threshold) and `drought_deficit_fixed_p{2,5,10,20,30}` (summed departures below
threshold, mm), plus five `drought_threshold_fixed_p{n}` per-gage scalars holding the
threshold values (mm/day). Interpretation keys:
- The five levels are **magnitude (non-exceedance) percentiles** of the 7-day-smoothed
  flow distribution, mirroring the U.S. Drought Monitor ladder: p30≈D0, p20≈D1,
  p10≈D2, p5≈D3, p2≈D4. The 10 % flow is the flow *exceeded* 90 % of the time —
  same convention as the `Q{n}` columns, opposite to how drought papers label
  exceedance thresholds.
- They are a **severity ladder, not independent measurements**: duration and deficit
  are non-decreasing in the level by construction and strongly correlated across
  levels. Pick one or two levels for headline analysis rather than treating all five
  as evidence.
- **Duration vs deficit**: rising deficit with flat duration = droughts deepening but
  not lengthening; the reverse = longer but shallower. Deficit is the only
  magnitude-weighted low-flow measure in the output.
- **`drought_duration_fixed_p10` is largely redundant with the existing low-pulse
  columns** — measured over 200,834 gage-years: r = 0.979 against
  `n_low_pulses_all × dur_low_pulses_all`, and within a gage the annual series track at
  median r = 0.994 (disagreement ≈ 12 % of the interannual SD), so don't present the two
  as independent evidence. The other four levels ARE distinct (r = 0.71–0.90 against the
  same pulse pair). The family's new content is `drought_deficit_*` (no counterpart
  anywhere) and the four non-p10 severity levels — the pulse metrics exist only at the
  10th percentile. **p10 is kept deliberately** (user decision, July 2026: don't prune a
  metric merely because a neighbouring family overlaps it, and keep the severity ladder
  complete) — so treat the overlap as a documentation caveat, not a defect.
- **Zeros are real values**: a year with no sub-threshold day reports 0, and that
  dense zero-including series is what carries the drought trend. Do not read 0 as
  missing.
- **Intermittent gages**: where most days are zero flow the low-level thresholds are
  exactly 0, and the strict `<` comparison then gives 0 duration/deficit in every
  year. Check `drought_threshold_fixed_p{n}` before interpreting those zeros.
- **Record-dependent** (thresholds come from the run's own window) — see the standard
  products section below; never compare drought values across the two windows.
- `drought_deficit_*` is unit-carrying (mm only where `area_normalized = TRUE`); the
  durations are scale-invariant and valid everywhere.
- **Two documented deviations from the source paper**: only the FIXED (whole-record)
  thresholds are implemented — the paper's variable day-of-year thresholds draw on one
  value per calendar day per year, so at 20–46 years the low levels carry very large
  sampling uncertainty and 2 % falls below the smallest plotting position `1/(n+1)`
  altogether (needs ≥ 49 years), which this project declines to extrapolate through —
  and aggregation is by WATER year, not the paper's climate year (Apr–Mar), so a
  drought peaking in September and running into October is split across two annual
  values.
- Weak sanity anchor: a level-`p` threshold sits below ~`p` % of all days by
  construction, so `drought_duration_fixed_p10_mean` averages ≈36.5 days/year. Treat a
  large departure as a prompt to investigate, not proof of a bug — and note the check is
  near-circular (the threshold is a percentile of the very series it is counted
  against), so passing it does NOT confirm the smoothing or plotting position.
  Intermittent gages legitimately fall short of it.

### The Two Standard Output Products (July 2026; #1 re-promoted Aug 2026)

Two production windows are the project's standard products (both @ 60% qualifying
fraction, 20-year floor, one self-contained folder each):
- **#1 WY 1993–2025** (6,678 gages, **1,653 columns**) — recent-era coverage.
  Promoted 2026-08-10 from the 28 Jul drought validation run
  (`processedOuts_drought_28jul2026`), **superseding the 22 Jul folder**
  (`processedOuts_22jul2026`, 1,488 columns). Same window, config, and gage set;
  the only difference is the +165-column streamflow drought family.
- **#2 WY 1980–2025** "entire period of record" (6,250 gages, **1,653 columns**) —
  long-established gages. Regenerated 2026-08-11 in
  `processedOuts_1980_2025_11aug2026`, **superseding the 22 Jul folder**
  (`processedOuts_1980_2025_22jul2026`, 1,488 columns). Both standards now carry
  the drought family.

Provenance caveat: #1 was built from the original Daymet parquet, #2 from a
parquet rebuilt out of the annual Daymet CSVs after the original was found
truncated. A controlled replay showed the rebuild shifts results by ≤ 3.4e-13
(CHANGELOG → August 2026), so this does not affect interpretation — but the two
products' `timing.json` provenance blocks record different climate inputs.

Interpretation keys:
- **Neither is a subset of the other**: 5,771 gages shared; 907 only in #1 (late
  starters — the 60% fraction is anchored at the WINDOW start, so e.g. 25 valid
  years / 46 possible = 0.54 fails the 1980 window while 25/33 = 0.76 passes the
  1993 window); 479 only in #2 (1980s records that ended before accumulating 20
  years inside 1993+). The longer window has FEWER gages because a longer window
  demands a longer record (60% × 46 ≈ 28 years vs the 20-year floor binding in
  #1).
- **Annual values are mutually consistent** (verified 2026-07-23): for shared
  gages, #2's annual parquet contains every (gage, signature, year) row of #1's,
  and overlap-year values are bit-identical for all within-year-computable
  signatures (flow volumes/percentiles, timing, BFI, FDC, flashiness, runoff
  ratios, snow).
- **Record-dependent signatures legitimately differ between the windows even in
  shared years**: `*_all` pulse metrics (period-of-record percentile
  thresholds), elasticity (record-mean Q̄/P̄ normalization), the
  recession-parameterized BFIs (whole-record alpha), and the drought metrics
  (whole-record percentile thresholds). Never mix these across window products.

## Interpreting Results

### Trend Statistics

- **Positive slope**: Increasing trend over time
- **Negative slope**: Decreasing trend over time
- **p-value < 0.05**: Statistically significant trend
- **Theil-Sen vs Linear**: Theil-Sen is more robust to outliers

### Common Interpretations

| Signature | Positive Trend | Negative Trend |
|-----------|----------------|----------------|
| Qann | Increasing annual flow | Decreasing annual flow |
| BFI | More groundwater contribution | Less groundwater contribution |
| flashinessRB | Becoming more flashy | Becoming less flashy |
| D50_day | Flow center shifting later | Flow center shifting earlier |

## Data Requirements

### Input Format

Daily streamflow data with columns:
- `gage_id`: Unique gage identifier
- `date`: Observation date
- `Q`: Daily discharge in **mm/day**
- `water_year`: Water year (Oct 1 - Sep 30)
- `month`: Calendar month (1-12)
- `dowy`: Day of water year (1-366)

### Quality Thresholds

All year qualification is handled centrally by `preprocess_daily_data()` (single source of truth):
- **Minimum years**: 20+ water years per gage
- **Year rejection**: >30 raw NAs, >3-day gaps, negative Q, boundary NAs
- **Interpolation**: Internal gaps ≤3 days filled by linear interpolation
- **Climate years**: Separate `valid_climate_years` set (same rules applied to PPT)
- **Trend completeness**: ≥60% non-NA annual values overall (lowered from 80% in July 2026) and ≥80% in first/last decade for trend stats
- **Snow record-anchored decade gate** (July 2026): the 10 threshold-dependent snow metrics (timing/melt/regime) also require ≥ decade_min_fraction of the SWE-valid years in the record's first AND last decade to be snowy/computable, else their 6 trend stats are NaN (mean/median + Pettitt survive). Linked to the same decade knob as the streamflow gate; anchored to the SWE record, not the metric's own span
- **Seasonal completeness**: Seasons with <80% raw observations flagged → seasonal metrics set to NA
- **No per-signature min_days**: Preprocessor is the only year filter (no additional gates in signature functions)

## Troubleshooting

### Common Issues

**Q: Why are all my statistics NA?**
A: Check that you have enough data. Each metric requires at least 3 non-NA annual values for statistics.

**Q: Why do Theil-Sen and linear slopes differ significantly?**
A: Theil-Sen is robust to outliers. Large differences indicate outliers or non-linear trends in your data.

**Q: How do I interpret very small p-values (e.g., 1e-10)?**
A: These indicate highly significant trends. However, statistical significance doesn't imply practical importance.

**Q: My timing metrics (D50_day) seem wrong. What's happening?**
A: Check that `dowy` (day of water year) is calculated correctly. Day 1 = October 1, not January 1.

**Q: Why does a Canadian gage have an absurdly large Qann (tens of thousands to millions)?**
A: Its `area_normalized` column is FALSE — HYDAT publishes no drainage area for that station (mostly canals, dam outflows, and channel splits like the St. Lawrence or Mackenzie), so flow was left in raw m³/s instead of mm/day. 37 such gages exist in the current output, by design (no backfill). Filter on `area_normalized == TRUE` before comparing unit-carrying signatures (volumes, percentiles, log_a) across gages; dimensionless Q-only signatures (BFI, FDC slopes, timing, TQmean, recession b) are still valid. All Q-to-PPT signatures (runoff ratios, elasticity, Q-P seasonality, avg_storage) are NA **by design** for these gages — a unit gate (`area_normalized` argument to `calculate_all_signatures()`, July 2026) skips them because Q (m³/s) and PPT (mm) units don't match. Note: 10 of the 37 have m³/s totals small enough to pass `flagged_for_qann_range` unflagged — the flag alone is not sufficient.

### QA/QC Flags

The output may include flags for:
- `flagged_for_qann_range`: Annual flow outside expected range (0-2000 mm) — does NOT catch all un-normalized (raw m³/s) Canadian gages; check `area_normalized` too
- `flagged_for_bfi_range`: Baseflow index outside [0, 1]
- `area_normalized`: FALSE = flow in raw m³/s (no drainage area in HYDAT) — exclude from cross-gage comparison of unit-carrying signatures
- `processing_status`: "success" or error description

## Language-Specific Notes

### Julia (Canonical)

```julia
using StreamflowSignatures
df = CSV.read("streamflow_data.csv", DataFrame)
df = add_water_year_columns(df)
signatures = calculate_flow_vols_by_year(df)
```

### Python

```python
from streamflow_signatures import calculate_flow_vols_by_year, add_water_year_columns
df = add_water_year_columns(df, date_col="date")
signatures = calculate_flow_vols_by_year(df)
```

### R

```r
source("config.R")
source("R/helperFunctions.R")
result <- process_signatures_from_parquet(parquet_file, metadata_file, output_file)
```

## References

- **SIGNATURES.md**: Detailed signature documentation
- **SIGNATURE_GUIDELINES.md**: Collaborative guidelines from hydrology team
- **METHODOLOGY.md**: Mathematical specifications (planned)

## Cross-Language Alignment

When debugging divergences between Julia, Python, and R implementations, use the **cross-language-alignment** skill. It provides a systematic 4-phase workflow, a 9-pattern divergence taxonomy, and the Iron Rule (one fix at a time, benchmark after each, revert if worse).

## Getting Help

If this skill doesn't answer your question, check:
1. The SIGNATURES.md file for detailed metric documentation
2. The CHANGELOG.md for known issues and recent fixes
3. The DEVELOPMENT.md for architecture and workflow details
4. The **cross-language-alignment** skill for debugging cross-language divergences

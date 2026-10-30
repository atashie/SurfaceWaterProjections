# Streamflow Signatures Reference

Detailed documentation of all hydrological signatures calculated by this project.

> **Comprehensive Documentation**: See `SIGNATURE_GUIDELINES.md` for mathematical formulations, parameter choices, and scientific rationale. This file is auto-synced from the [collaborative Google Doc](https://docs.google.com/document/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub).

## Overview

### Signature categories

Every output belongs to one of **nine signature categories** (the eight adopted
2026-09-10 from the co-authors' categorization plus Climate, 2026-10-12; canonical mapping in
[`signature_categories.csv`](signature_categories.csv)). The numbered sections of this
file are **function families** — one per computing function — and each category
gathers one or more of them:

| Category | Function families (sections below) | Annual signatures | Per-gage scalars |
|----------|-------------------------------------|-------------------|------------------|
| Flow Volume | 1 Flow Volumes; Negative Flow Days (`calculate_negative_days`); TQmean from 5 Pulse Metrics | 21 + 1 + 1 | 4 (`season_excluded_years_*`) |
| Flow Duration | 2 Flow Duration Curve | 3 | — |
| Storage | 3 Baseflow (both variants); 4 Recession; 8 Q-PPT Relationships (runoff ratios); 11 Average Storage | 4 + 7 + 5 + 1 | 1 + 6 + 1 |
| Flashiness | 5 Pulse Metrics (pulses + reversals); 6 Flashiness; 9 Streamflow Elasticity | 13 + 1 + 2 | 3 (elasticity) |
| Drought | 13 Streamflow Drought | 10 | 5 (thresholds) |
| Flow Timing | 7 Flow Timing | 15 | — |
| Precipitation Streamflow | 10 Q-P Seasonality | 2 | — |
| Snow | 12 Snow Metrics | 14 | — |
| Climate | 14 Climate Metrics (precipitation + temperature families; October 2026) | 20 + 19 | — |

Totals: 139 annual signatures (100 through August 2026 + 39 climate, 2026-10-12) and
20 scalars; the 21st scalar, `ice_affected_days_total`, is a preprocessing diagnostic
outside the categories. Nine categories, 17 computing function families.
*Category* is the term used by the manuscript, the guidelines document, and the
HydroShare data dictionary; *function family* names the computing function.

Each signature produces **8 statistics** via `generate_stats()`:

> **Per-year annual values** (July 2026): the annual series behind these statistics
> are also exported as a long-format parquet (`{prefix}_signatures_annual.parquet`:
> `gage_id, signature, water_year, value`) written by the Julia benchmark alongside
> the summary CSV. NaN and absent-row are equivalent ("year not computable") — some
> signatures emit NaN placeholders, others omit the row. See
> `docs/DEVELOPMENT.md` → Annual Values Export for schema and semantics.

| Suffix | Statistic | Method |
|--------|-----------|--------|
| `_senn_slp` | Theil-Sen slope | `zyp::zyp.sen` |
| `_linear_slp` | Linear regression slope | `lm()` |
| `_spearman_rho` | Spearman's rho | `cor.test` |
| `_spearman_pval` | Spearman p-value | `cor.test` |
| `_mk_rho` | Mann-Kendall tau | `Kendall::MannKendall` |
| `_mk_pval` | Mann-Kendall p-value | `Kendall::MannKendall` |
| `_mean` | Arithmetic mean | `mean()` |
| `_median` | Median | `median()` |

---

## 1. Flow Volumes

**Function**: `calculate_flow_vols_by_year`

### Metrics

| Metric | Description |
|--------|-------------|
| **Qann** | Annual total streamflow |
| **Qwin** | Winter total (Dec-Feb) |
| **Qspr** | Spring total (Mar-May) |
| **Qsum** | Summer total (Jun-Aug) |
| **Qfal** | Fall total (Sep-Nov) |
| **Q1-Q99** | Flow percentiles (Q1, Q5, Q10, Q20, Q25, Q30, Q40, Q50, Q60, Q70, Q75, Q80, Q90, Q95, Q99) |
| **Q95_Q10** | Difference of high to low flow percentiles (Q95 - Q10) |

### Notes
- **Q95_Q10 naming**: The standard column name is `Q95_Q10` (underscore), used by Julia, Python, and rpkg.

### Units
- Qann, Qwin, Qspr, Qsum, Qfal: mm (total over period, summed from daily mm/day values)
- Q percentiles: mm/day (daily flow at each percentile)

### Data Quality
Year qualification is handled centrally by `preprocess_daily_data()` before any signature functions run.

### Diagnostics

| Metric | Description |
|--------|-------------|
| **season_excluded_years_winter** | Count of years where winter failed 80% completeness threshold (per-gage scalar) |
| **season_excluded_years_spring** | Count of years where spring failed 80% completeness threshold (per-gage scalar) |
| **season_excluded_years_summer** | Count of years where summer failed 80% completeness threshold (per-gage scalar) |
| **season_excluded_years_fall** | Count of years where fall failed 80% completeness threshold (per-gage scalar) |

---

## 2. Flow Duration Curve (FDC)

**Function**: `analyze_fdc_trends`

Constructs an empirical flow duration curve (FDC) for each water year and fits the slope of the curve (in log space) over three exceedance-probability ranges. FDC slope characterizes flow variability: steeper (more negative) slopes indicate flashier, more variable regimes, while flatter slopes indicate damped, baseflow-dominated regimes. Follows Yilmaz et al. (2008) and Sawicz et al. (2011), though the segment definitions and fitting method differ (see Method).

### Metrics

| Metric | Description | Exceedance range |
|--------|-------------|------------------|
| **FDCall** | Slope of the annual FDC across the full range | 0.0–1.0 |
| **FDC90th** | Slope over the low-flow tail (driest ~10% of days) | ≥ 0.90 |
| **FDCmid** | Slope over the mid-segment | 0.20–0.80 |

### Method

1. Within each water year, drop NaN and negative daily flows (zeros are retained); require **≥10 valid daily values**, otherwise all three slopes for that year are NA.
2. Sort the valid daily flows in descending order and assign exceedance probabilities via the Weibull plotting position, `p = i / (n + 1)` for `i = 1…n`.
3. Log-transform flows as `log10(Q + 1e-10)`; the small constant keeps zero-flow days finite.
4. For each segment, subset to its exceedance range and fit the slope as the coefficient of an ordinary least-squares regression of `log10(Q)` on exceedance probability. A segment with **<3 points** yields NA for that year.

This differs from the two-point percentile-difference formula of Sawicz et al. (2011): regressing over a full segment is less sensitive to the choice of endpoint percentiles.

### Units

`log10(mm/day)` per unit exceedance probability. Slopes are typically negative (flow decreases as exceedance probability increases); a trend toward more negative values indicates increasing flow variability.

### Data Quality

Year qualification is handled centrally by `preprocess_daily_data()` before this function runs. FDC applies no additional per-year completeness filter beyond the ≥10-valid-value and ≥3-points-per-segment guards inside the slope fit.

### Reference
Yilmaz, K.K., Gupta, H.V., & Wagener, T. (2008). A process-based diagnostic approach to model evaluation: application to the NWS distributed hydrologic model. *Water Resources Research*, 44(9), W09417.

---

## 3. Baseflow

**Function**: `analyze_baseflow_indices`

### Metrics

| Metric | Description | Parameters |
|--------|-------------|------------|
| **BFI_Eckhardt** | Baseflow index using Eckhardt recursive digital filter | BFImax=0.8, a=0.98 |
| **BFI_LyneHollick** | Baseflow index using Lyne-Hollick filter | alpha=0.925, 2 passes |
| **BFI_Eckhardt_param** | Baseflow index using Eckhardt filter with recession-derived alpha | a=recession_alpha, BFImax=0.8 |
| **BFI_LyneHollick_param** | Baseflow index using Lyne-Hollick filter with recession-derived alpha (heuristic) | alpha=recession_alpha, 2 passes |

### Expected Relationship
BFI_Eckhardt < BFI_LyneHollick (validated in QA/QC)

### Recession-Parameterized BFI

**Function**: `analyze_baseflow_indices_with_parameters`

The parameterized BFI signatures use a recession-derived discrete filter constant instead of fixed defaults. The filter constant `alpha` is estimated as the median of `Q_{i+1}/Q_i` ratios across all recession events in the entire record (stored as `recession_alpha_point_cloud_linear_reservoir`).

- **Eckhardt**: `alpha` directly replaces the fixed `a=0.98`. BFImax remains fixed at 0.8.
- **Lyne-Hollick**: `alpha` replaces the fixed `alpha=0.925`. This is a heuristic — the L-H filter parameter has no physical derivation from recession analysis (Nathan & McMahon, 1990).
- Gages with insufficient recession data (<=10 alpha pairs in the whole record) produce NA for all parameterized BFI values.

### References
- Eckhardt, K. (2005). How to construct recursive digital filters for baseflow separation.
- Lyne, V., & Hollick, M. (1979). Stochastic time-variable rainfall-runoff modelling.
- Collischonn, W., & Fan, F.M. (2013). Defining parameters for Eckhardt's digital baseflow filter.

---

## 4. Recession

**Function**: `analyze_recession_parameters`

Analyzes recession curve behavior using the dQ/dt = -a*Q^b relationship (event
identification, free-fit exponent, concavity and log(a) seasonality follow Tashie et
al. 2020).

**Parameter conventions (July 2026)**: the exponent **b** (and concavity) are estimated
with FREE power-law fits, but all **alpha** outputs assume a **linear reservoir (b fixed
at 1)** across all locations and periods: `log(a) = median(log(-dQ/dt) - log(Q))` — no
regression involved. Rationale: log(a) is the intercept of a regression whose slope is
b, so free-fit alpha estimates are convolved with b (uncertainty and temporal variation
in b leak into alpha); fixing b = 1 decouples them. Under the forward-difference
discretization each b=1 point equals `log(1 - Q_{i+1}/Q_i)`, so per-year
`log_a_pointcloud` approximates `log(1 - alpha_linear)`. The relation is approximate,
not a same-sample identity: `alpha_linear` pools pairs from ALL identified events,
while the point cloud only pools events whose free power-law fit succeeded, and the
median only commutes exactly with the monotone transform at odd pair counts.

### Metrics

| Metric | Description |
|--------|-------------|
| **log_a_pointcloud** | Recession rate parameter, b=1 fixed (point cloud method) |
| **log_a_events** | Recession rate parameter, b=1 fixed (event-based method) |
| **b_pointcloud** | Recession exponent, free fit (point cloud method) |
| **b_events** | Recession exponent, free fit (event-based method) |
| **concavity** | Difference in b between first and second halves of recession |
| **n_recession_events** | Count of recession events per water year (independent of min_events gate) |
| **alpha_linear** | Discrete recession constant under linear reservoir (b=1); per-year median of Q_{i+1}/Q_i ratios from point cloud |
| **log_a_seasonality_amplitude_all** | Seasonal amplitude in recession rate (all data) |
| **log_a_seasonality_amplitude_first_half** | Seasonal amplitude (first half of record) |
| **log_a_seasonality_amplitude_last_half** | Seasonal amplitude (last half of record) |
| **log_a_seasonality_minimum_all** | Seasonal minimum day for recession rate (all data) |
| **log_a_seasonality_minimum_first_half** | Seasonal minimum day (first half of record) |
| **log_a_seasonality_minimum_last_half** | Seasonal minimum day (last half of record) |

### Notes
- Seasonality metrics are single values (exceptions to 8-statistic rule); the sinusoid
  is fit to the per-event **b=1** log_a values (July 2026 — previously per-event
  free-fit intercepts), so seasonal alpha signals are no longer confounded by seasonal
  variation in b
- Documented in `config.R` as `EXPECTED_RECESSION_SEASONALITY`
- Requires minimum 25 recession events (`RECESSION_MIN_EVENTS`) for all metrics except `n_recession_events`
- `n_recession_events` is computed independently of the min_events gate (useful for gages with few events)
- `alpha_linear` per-year values require >10 valid alpha pairs per year; same 25-event gate as other recession metrics for trend statistics
- `recession_alpha_point_cloud_linear_reservoir` is a per-gage scalar (whole-record median of Q_{i+1}/Q_i ratios across all recession events). Computed independently of the 25-event gate (only requires >10 alpha pairs). Used to parameterize `BFI_Eckhardt_param` and `BFI_LyneHollick_param`.
- Recession events are identified as contiguous windows where both Q and |dQ/dt| are monotonically decreasing, with a minimum length of 5 days

### Reference
Tashie, A., Pavelsky, T., & Band, L.E. (2020). An empirical reevaluation of streamflow recession analysis at the continental scale. *Water Resources Research*, 56(1), e2019WR025448.

---

## 5. Pulse Metrics

**Function**: `calculate_pulse_metrics`

### Metrics

Two sets of pulse metrics are computed: `*_year` variants use per-year percentile thresholds, while `*_all` variants use period-of-record percentiles (as specified in the guidelines). Both are retained for scientific flexibility.

| Metric | Description | Threshold |
|--------|-------------|-----------|
| **n_high_pulses_year** | Count of high pulse events (per-year percentiles) | > year-specific 90th percentile |
| **n_low_pulses_year** | Count of low pulse events (per-year percentiles) | < year-specific 10th percentile |
| **dur_high_pulses_year** | Mean duration of high pulses (per-year percentiles) | days |
| **dur_low_pulses_year** | Mean duration of low pulses (per-year percentiles) | days |
| **n_high_pulses_all** | Count of high pulse events (period-of-record percentiles) | > period-of-record 90th percentile |
| **n_low_pulses_all** | Count of low pulse events (period-of-record percentiles) | < period-of-record 10th percentile |
| **dur_high_pulses_all** | Mean duration of high pulses (period-of-record percentiles) | days |
| **dur_low_pulses_all** | Mean duration of low pulses (period-of-record percentiles) | days |
| **TQmean** | Percentage of days with flow above annual mean | % |
| **Flow_Reversals_annual** | Annual count of flow direction changes | count |
| **Flow_Reversals_winter** | Winter (Dec-Feb) flow direction changes | count |
| **Flow_Reversals_spring** | Spring (Mar-May) flow direction changes | count |
| **Flow_Reversals_summer** | Summer (Jun-Aug) flow direction changes | count |
| **Flow_Reversals_fall** | Fall (Sep-Nov) flow direction changes | count |

### Notes
- The `*_all` variants match the guidelines' requirement to "calculate annual metrics using period-of-record percentiles"
- The `*_year` variants provide complementary per-year sensitivity analysis
- Seasonal flow reversals complement the annual count with seasonal decomposition

---

## 6. Flashiness

**Function**: `analyze_flashiness_trends`

### Metrics

| Metric | Description |
|--------|-------------|
| **flashinessRB** | Richards-Baker flashiness index |

### Formula
```
R-B Index = sum(|Q_i - Q_{i-1}|) / sum(Q_i)
```

Sum of absolute day-to-day changes divided by total flow.

### Reference
Baker, D.B., et al. (2004). A new flashiness index: characteristics and applications to midwestern rivers and streams.

---

## 7. Flow Timing

**Function**: `analyze_flow_timing_trends`

### Metrics

| Metric | Description |
|--------|-------------|
| **D1_day** | Day of water year when cumulative flow reaches 1% |
| **D5_day** | Day of water year when cumulative flow reaches 5% |
| **D10_day** | Day of water year when cumulative flow reaches 10% |
| **D20_day** | Day of water year when cumulative flow reaches 20% |
| **D30_day** | Day of water year when cumulative flow reaches 30% |
| **D40_day** | Day of water year when cumulative flow reaches 40% |
| **D50_day** | Day of water year when cumulative flow reaches 50% (center of mass) |
| **D60_day** | Day of water year when cumulative flow reaches 60% |
| **D70_day** | Day of water year when cumulative flow reaches 70% |
| **D80_day** | Day of water year when cumulative flow reaches 80% |
| **D90_day** | Day of water year when cumulative flow reaches 90% |
| **D95_day** | Day of water year when cumulative flow reaches 95% |
| **D99_day** | Day of water year when cumulative flow reaches 99% |
| **D25_to_D75** | Duration between 25% and 75% cumulative flow (days) |
| **Dmax** | Day of maximum flow |

### Notes
- Day 1 = October 1 (start of water year)
- Day 365/366 = September 30 (end of water year)
- D1 may be near-constant (day 1-2) for flashy streams; D99 may be ~364-365 — trend statistics may be unreliable for constant series
- Percentiles are config-driven via `config/signatures_config.json` → `timing.d_percentiles`

---

## 8. Q-PPT Relationships

**Function**: `analyze_Q_PPT_relationships`

**Requires Climate Data** (PPT column via Daymet integration or Caravan)

### Metrics

| Metric | Code Name | Description |
|--------|-----------|-------------|
| **annual_runoff_ratio** | `annual_runoff_ratio` | Annual runoff ratio (Q/P) |
| **winter_runoff_ratio** | `winter_runoff_ratio` | Winter runoff ratio |
| **spring_runoff_ratio** | `spring_runoff_ratio` | Spring runoff ratio |
| **summer_runoff_ratio** | `summer_runoff_ratio` | Summer runoff ratio |
| **fall_runoff_ratio** | `fall_runoff_ratio` | Fall runoff ratio |

### Diagnostics

| Metric | Description |
|--------|-------------|
| **runoff_ratio_high_count** | Count of years where annual runoff ratio > 2.0 (per-gage scalar) |

### PPT Thresholds
- Annual: minimum 10mm PPT required (below returns NA)
- Seasonal: minimum 1mm PPT required (below returns NA)

### Data Quality
- Year qualification is handled centrally by `preprocess_daily_data()` (PPT-valid years).
- Climate frame: since the version-3 products (built 2026-10-27; `STREAMFLOW_CLIMATE_FRAME=outer`,
  the runners' default since 2026-10-28) PPT enters the screening on Daymet's own calendar, so
  a date absent from the streamflow record keeps its precipitation. Every earlier product
  joined PPT onto the streamflow ROWS and interpolated such days (storm days → biased low);
  the ratios inherit it there. See §14 → Data quality and caveats.

---

## 9. Streamflow Elasticity

**Function**: `calculate_streamflow_elasticity`

**Requires Climate Data**

### Metrics

| Metric | Description |
|--------|-------------|
| **elasticity_static** | Overall catchment elasticity (single value) |
| **elasticity_rolling** | Rolling window (11-year) elasticity trend (Sawicz departure-from-mean) |
| **elasticity_annual** | Year-over-year elasticity trend (consecutive-year differences) |

### Formulas

**Departure-from-mean (elasticity_static, elasticity_rolling)**:
```
E_i = ((Q_i - Q_mean) / (P_i - P_mean)) / (Q_mean / P_mean)
```

**Year-over-year (elasticity_annual)**:
```
E_t = ((Q_t - Q_{t-1}) / (P_t - P_{t-1})) / (Q_mean / P_mean)
```

### Diagnostics

| Metric | Description |
|--------|-------------|
| **elasticity_years_total** | Total years available for this gage (per-gage scalar) |
| **elasticity_years_low_ppt** | Years excluded due to annual PPT < 10mm (per-gage scalar) |

### Interpretation
- E ~ 1.0: Proportional response to precipitation changes
- E > 1.0: Amplified response (Q changes more than P)
- E < 1.0: Dampened response

### Notes
- `elasticity_static` is a single value (exception to 8-statistic rule)
- `elasticity_rolling` uses departure-from-mean within 11-year windows (Sawicz et al. 2011)
- `elasticity_annual` uses consecutive-year differences for year-to-year sensitivity (NOT Sawicz method)
- Diagnostics are per-gage scalars documented in `config.R` as `EXPECTED_ELASTICITY_DIAGNOSTICS`

### Reference
Sawicz, K., et al. (2011). Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA.

---

## 10. Q-P Seasonality

**Function**: `calculate_qp_seasonality`

**Requires Climate Data**

### Metrics

| Metric | Description |
|--------|-------------|
| **qp_slope_sd** | Standard deviation of monthly cumulative Q-P slopes |
| **qp_bimodality** | Bimodality coefficient of Q-P relationship |

### Interpretation
- **qp_slope_sd**: Higher values indicate stronger seasonality
- **qp_bimodality**: Values > 0.555 suggest seasonal/bimodal patterns

### Method
Calculated from 30-day rolling slopes of cumulative Q vs cumulative P

### Reference
Wrede, S., et al. (2015). Towards a common classification framework for hydrological models.

---

## 11. Average Storage

**Function**: `calculate_average_storage`

**Requires Climate Data**

### Metrics

| Metric | Description |
|--------|-------------|
| **avg_storage** | Mean annual catchment storage (mm) |

### Method
1. Calculate cumulative storage: S = cumsum(P - Q)
2. Interpolate annual storage at mean discharge for each water year
3. Average across years

### Known Limitation

> **Warning**: This calculation ignores evapotranspiration (ET), using only P - Q for the water balance. This simplification may overestimate storage in watersheds with significant ET losses.

Future improvements should incorporate:
- Hargreaves-Samani ET estimation from temperature data
- External ET products (e.g., MODIS ET)

### Reference
Peters, N.E., & Aulenbach, B.T. (2011). Water storage at the Panola Mountain Research Watershed, Georgia, USA.

---

## 12. Snow Metrics

**Function**: `calculate_snow_metrics`

**Requires SWE Data** (Daymet `swe` via the climate parquet; July 2026; all three languages since the August 2026 port campaign)

Fourteen per-water-year snow metrics computed from daily snow water equivalent.
ALL metrics operate on a thresholded series `SWE* = (SWE >= 10 mm) ? SWE : 0` —
days below the threshold are treated as snow-free for durations AND magnitudes
(domain decision, July 2026). Configuration: `config/signatures_config.json` → `snow`
(threshold, 60-day seasonal-spell rule, melt center-of-mass fraction, PPT floor).

### Metrics

| Metric | Description | Units | No-snow year |
|--------|-------------|-------|--------------|
| **swe_max** | Maximum daily SWE | mm | 0 (valid) |
| **swe_max_dowy** | Day of water year of the peak (ties → first day) | day | NA |
| **snow_cover_days** | Days with SWE ≥ threshold | days | 0 (valid) |
| **snow_on_dowy** | First day of the anchor spell (the continuous spell containing the peak) | day | NA |
| **snow_off_dowy** | First snow-free day after the anchor spell | day | NA |
| **melt_season_days** | snow_off_dowy − swe_max_dowy | days | NA |
| **melt_rate** | swe_max ÷ melt_season_days (net ablation rate between peak and snow-off) | mm/day | NA |
| **ssm** | Snow seasonality metric: (seasonal − ephemeral days) ÷ total snow days, where seasonal spells last ≥ 60 continuous days | −1…+1 | NA |
| **swe_apr1** | SWE on calendar April 1 (leap-safe) | mm | 0 (valid) |
| **melt_before_peak** | Total melt (daily SWE decreases) before the peak day | mm | NA |
| **melt_before_peak_pct** | melt_before_peak ÷ total water-year melt × 100 | % | NA |
| **melt_before_peak_to_max_swe** | melt_before_peak ÷ swe_max | – | NA |
| **melt_com_dowy** | Day cumulative melt reaches 50% of the water-year total | day | NA |
| **swe_max_to_ppt** | swe_max ÷ water-year total PPT (PPT-qualified years only; > 10 mm floor) | – | 0 (valid) |

### Method

1. Daily SWE for each SWE-valid water year (preprocessor `valid_swe_years` — same
   per-year NA rules as PPT: >30 raw NAs reject, internal gaps ≤ 3 days linearly
   interpolated, negative SWE rejects).
2. Threshold to SWE*; snow days are SWE* > 0; spells are maximal consecutive runs
   of snow days.
3. The **anchor spell** is the spell containing the annual max (first-max tie
   rule). Snow-on/off are **censored (NA)** when the anchor spell touches Oct 1
   (carryover snowpack predates the water year) or Sep 30 (snowpack persists past
   it) — a spell spanning the whole year censors both ends.
4. Melt increments `m_t = max(0, SWE*_{t−1} − SWE*_t)` are attributed to day `t`,
   within-year only (no attribution across the Sep 30 → Oct 1 boundary). The
   final melt-out step (crossing below the threshold) attributes the full
   remaining SWE* to the crossing day.
5. SSM pools ALL spells in the year (Hatchett 2021, Section 2.2; the 60-day
   seasonal threshold follows Petersky & Harpold 2018): +1 fully seasonal,
   −1 fully ephemeral.
6. **Record-anchored decade gate** (July 2026): the 10 timing/melt/regime metrics
   (all except `swe_max`, `snow_cover_days`, `swe_apr1`, `swe_max_to_ppt`)
   additionally require ≥ `decade_min_fraction` (the SAME
   `na_handling.trend_completeness` knob the streamflow decade gate uses — linked
   by design; 0.80 as shipped) of the SWE-valid years in BOTH the first and last
   decade of the gage's SWE record to be computable (snowy). Failing metrics have
   their 6 trend statistics set to NaN; mean/median and Pettitt fields are still
   computed. Anchored to the SWE record — NOT the metric's own non-NA span — so a
   gage whose snow disappears (or appears) mid-record cannot carry a trend
   conditioned only on its snow-present years (the clustered-years hole that the
   own-span gate and the 20-value stats floor leave open). Denominators count
   SWE-valid years in each window, so missing Daymet years don't count against
   snow presence; records spanning < 10 years skip the gate. Magnitude metrics
   are exempt: their dense zero-including series carry the snow-decline signal
   legitimately. Config: `snow.record_anchored_decade_gate` (absent → disabled).

### Data Quality & Caveats

- **Source**: Daymet V4 `swe` (kg/m² ≡ mm), calendar 1980–2023. **CORRECTION
  (July 2026 rerun finding): the Daymet parquet covers Canadian gages too** —
  e.g. the WY1980+ run has snow values for 5,622 gages (4,533 US + 1,089
  Canadian); earlier "US-only Daymet" claims were wrong. Gages without Daymet
  rows (any nationality) get all-NA snow metrics. Daymet SWE is a model product
  (mass-balance bookkeeping from Daymet's own precipitation/temperature) with
  documented biases, especially in mountain terrain.
- The spatial support of the per-gage Daymet series (basin-average vs
  point/centroid extraction) is not recorded in-repo — treat absolute magnitudes
  with caution; timing and trend signals are more robust.
- Daily ΔSWE understates melt when snowfall and melt co-occur within a day.
- The 10 mm operational threshold intentionally diverges from Hatchett's literal
  SWE > 0 rule (config-adjustable for sensitivity runs). Magnitudes are censored
  at the threshold: a year peaking at 8 mm is operationally snow-free
  (swe_max = 0), and sub-threshold April 1 values report 0.
- Zero-vs-NA policy: magnitude metrics emit valid zeros (dense series at
  snow-free gages, like D1_day's near-constant behavior); timing/melt/regime
  metrics emit NA — the trend-completeness gate (applied; NOT exempt)
  suppresses trends at marginal-snow gages, and the record-anchored decade gate
  (Method step 6) additionally suppresses trends where snowy years are absent
  from the first or last decade of the SWE record.
- WY 1980 and WY 2024 are partial in the Daymet parquet → automatically
  SWE-invalid (same boundary behavior as the climate signatures).
- **Daymet uses a 365-day calendar** (confirmed against the source CSVs,
  August 2026): leap years keep Feb 29 and instead **drop Dec 31**. Every
  Daymet year is exactly 365 days, so the SWE (and PPT) series carries a
  one-day hole on Dec 31 of each leap year, which the preprocessor absorbs as
  an internal gap ≤ 3 days (linearly interpolated). Melt increments are
  therefore interpolated across that boundary in leap years rather than
  observed. Longstanding behavior — not introduced by the 2026-08-11 parquet
  rebuild, which reproduces it faithfully.
- Days absent from the streamflow record keep their SWE (and PPT) since the version-3
  products (`STREAMFLOW_CLIMATE_FRAME=outer`); in every earlier product they were absent from
  SWE as well and interpolated although Daymet has them (§14 → Data quality and caveats).
- Snow metrics only run on an explicitly provided, SWE-valid-year-filtered frame
  (`snow_data`); an SWE column in the main gage frame is never used implicitly.

### References

- Hatchett, B.J. (2021). Seasonal and Ephemeral Snowpacks of the Conterminous United States. *Hydrology*, 8(1), 32.
- Petersky, R., & Harpold, A. (2018). Now you see it, now you don't: a case study of ephemeral snowpacks and soil moisture response in the Great Basin, USA. *Hydrology and Earth System Sciences*, 22, 4891–4906.

---

## 13. Streamflow Drought

**Function**: `calculate_drought_metrics` (July 2026; all three languages since the August 2026 port campaign)

Per-water-year drought duration and deficit against **fixed** percentile thresholds,
following Adelsperger et al. (in review). Ten metrics: two measures × five severity
levels. Configuration: `config/signatures_config.json` → `drought`.

### Metrics

| Metric | Description | Units |
|--------|-------------|-------|
| **drought_duration_fixed_p{2,5,10,20,30}** | Days in the water year with 7-day-smoothed Q below the threshold | days |
| **drought_deficit_fixed_p{2,5,10,20,30}** | Sum of flow departures below the threshold, `Σ (Q_thr − Q_smooth)` | mm |

### Diagnostics

| Metric | Description |
|--------|-------------|
| **drought_threshold_fixed_p{2,5,10,20,30}** | The threshold value itself (mm/day, per-gage scalar) — makes the thresholds auditable and identifies gages whose low-level threshold is exactly 0 |

### Severity levels

Thresholds are **magnitude (non-exceedance) percentiles**, matching the paper and this
document's `Q{n}` columns: the 10 % flow is the flow *exceeded* 90 % of the time. The
five levels mirror the operational U.S. Drought Monitor classes:

| Level | Percentile | USDM analogue |
|---|---|---|
| `p30` | 30 % | D0 — abnormally dry |
| `p20` | 20 % | D1 — moderate |
| `p10` | 10 % | D2 — severe |
| `p5` | 5 % | D3 — extreme |
| `p2` | 2 % | D4 — exceptional |

### Method

1. **7-day centered smoothing** of daily Q, applied to the *continuous* date-indexed
   series so the Oct 1 boundary creates no artificial discontinuity. The window is
   applied within maximal runs of consecutive dates and never averages across a
   temporal gap (rejected years, duplicate or missing dates all break a run); it
   shrinks at run edges, and a day whose in-run window holds fewer than
   `smoothing_min_valid_days` (4) values is NaN. Because the preprocessor supplies
   complete daily grids for qualifying years, only runs shorter than the minimum
   produce NaN in practice.
2. **Thresholds** are percentiles of the smoothed values pooled over the whole record
   (the paper's *fixed* method), computed with the **unbiased Weibull plotting
   position** `p_i = i/(n+1)` — Hyndman & Fan definition 6, i.e.
   `quantile(x, p; alpha=0, beta=0)`. This deliberately differs from the type-7
   (linear-interpolation) default used elsewhere in this project, and the two diverge
   most in exactly the low tail these thresholds occupy. A probability below the
   smallest plotting position `1/(n+1)` is not interpolable and yields NaN rather than
   silently clamping to the record minimum (unreachable for the fixed method, whose
   pool holds ~7,000–17,000 values; retained as a defensive invariant).
3. **Per water year and level**: count days with `Q_smooth < Q_thr` (strict), and sum
   `Q_thr − Q_smooth` over those days. Gages with fewer than
   `min_years_for_threshold` (10) qualifying years produce NaN for all drought metrics.

### Interpretation

- Rising `deficit` with flat `duration` = droughts deepening but not lengthening;
  the reverse = longer but shallower droughts.
- The five levels are a severity ladder, not independent measurements — duration and
  deficit are non-decreasing in the level by construction and strongly correlated
  across levels. Prefer one or two levels for headline analysis.
- **`drought_duration_fixed_p10` is largely redundant with the existing pulse metrics.**
  Measured over the full WY 1993–2025 benchmark (200,834 gage-years): it correlates with
  `n_low_pulses_all × dur_low_pulses_all` at r = 0.979 (ρ = 0.971), and — the test that
  matters for trend work — *within* a gage the two annual series track each other at
  median r = 0.994 (p10 of the distribution 0.971), disagreeing by a median of 11.7 % of
  the gage's own interannual SD. The other four levels are measurably distinct
  (r = 0.712 / 0.902 / 0.846 / 0.731 for p2 / p5 / p20 / p30 against the same pulse pair),
  so only p10 collapses onto the existing signal. Redundancy is an aggregate statement:
  32.5 % of gage-years agree exactly and the maximum disagreement is 318 days. On
  intermittent gages the two measures agree exactly in 99.87 % of gage-years (13 of 9,652
  differ, max 8 days).
  The genuinely new content of this family is (a) **`drought_deficit_*`**, the only
  magnitude-weighted low-flow measure in the output — it separates a long shallow drought
  from a short severe one — and (b) the **four non-p10 levels**, since the pulse metrics
  use a single 10th-percentile threshold. Prefer those for new analysis, and do not
  present `drought_duration_fixed_p10` and the pulse pair as independent evidence for the
  same conclusion.
  **`drought_duration_fixed_p10` is retained deliberately** (user decision, July 2026):
  cross-family redundancy is not grounds for removal — the same quantity computed by two
  documented methods has independent value, and the drought family stays internally
  complete across all five USDM rungs rather than carrying a hole at D2.
- **Weak internal consistency check** (not a correctness proof): because a level-`p`
  threshold is a percentile *of the same series it is then counted against*,
  `drought_duration_fixed_p10_mean` averages ≈36.5 days/year (`p × 365.25`) across a
  large gage sample. This is largely circular — it would also hold if the threshold and
  the comparison both used the wrong series (raw instead of smoothed), and switching
  plotting position moves it by only a few observations. It catches gross mismatches
  (e.g. a raw-flow threshold compared against smoothed flow) and nothing subtler. It can
  also legitimately fall short of `p × 365.25` on tied/intermittent records, where the
  strict `<` excludes a mass of equal values at the threshold.

### Data Quality & Caveats

- **Record-dependent**: thresholds are computed from the run's own window (matching the
  paper's per-analysis-period thresholds), so values are comparable *within* a product
  but **must not be compared across** the WY 1993–2025 and WY 1980–2025 standard
  products — same caveat as `*_all` pulses, elasticity, and the parameterized BFI.
- **Water year vs the paper's climate year**: the source aggregates over climate years
  (Apr 1 – Mar 31) so the summer–autumn low-flow season is not split. This
  implementation uses the water year (Oct 1 – Sep 30) for consistency with every other
  signature and with the annual-values parquet key, so a drought peaking in September
  and persisting into October is **split across two annual values**. Day-level
  indicators — and therefore record totals — are identical under either grouping; the
  annual series and its trend statistics are not.
- **The paper's variable (day-of-year) threshold method is not implemented.** Its
  per-day sample is one value per year, so at 20–46 years of record the low levels are
  estimated with very large sampling uncertainty, and the 2 % level falls below the
  smallest Weibull plotting position `1/(n+1)` entirely (2 % would need ≥49 years, 5 %
  ≥19). Standard type-6 behaviour there is to clamp to the sample minimum; this project
  instead refuses to extrapolate (`below_plotting_range_policy: "na"`), which is a
  deliberate estimability policy layered on top of type 6, not part of type 6 itself.
  So the variable method is *too uncertain under that policy*, not mathematically
  impossible.
- **Units**: `drought_deficit_*` is mm only where Q is area-normalized; for the
  `area_normalized = FALSE` gages it carries m³/s·day and is not cross-gage comparable.
  The duration metrics are scale-invariant (the threshold is a percentile of the same
  series) and valid everywhere.
- **Zeros are values, not gaps**: a year with no sub-threshold days reports 0, and the
  dense zero-including series legitimately carries the drought signal. No snow-style
  record-anchored trend gate is applied.
- **Intermittent gages**: where >~70 % of days are zero flow the low-level thresholds
  are exactly 0, and the strict `<` comparison then reports 0 duration and 0 deficit
  for every year. Use `drought_threshold_fixed_p{n}` to identify these gages rather
  than reading the zeros as "no drought".
- Trend completeness and the 20-value stats floor both apply (this family is **not**
  exempt).

### References

- Adelsperger, S., et al. (in review). A novel severity-based approach for assessing streamflow drought characteristics and drivers.
- Laaha, G., et al. (2017). Unbiased plotting positions for low-flow frequency analysis.
- U.S. Drought Monitor percentile classes: https://droughtmonitor.unl.edu

---

## 14. Climate Metrics

**Functions**: `calculate_precipitation_metrics` (20 bases), `calculate_temperature_metrics`
(19 bases) — `julia/src/climate.jl` (October 2026, Julia first; Python and rpkg ports
cross-checked bit-for-bit on the annual values the same day)

**Requires Daymet `prcp`, `tmin`, `tmax`** via the climate parquet. Category **Climate**
(the ninth; guidelines §3.9, synced 2026-10-12). Pure-climate series: no Q enters them,
so they are computed for the un-normalized gages too, and no Q-derived seasonal mask
applies. Configuration: `config/signatures_config.json` → `climate` (dry-day threshold,
wet-day thresholds, multi-day window, cumulative fraction, frost temperature and cut-off
month, GDD base, decade gate) and `na_handling.temperature_na_policy`.

Conventions: seasons are the calendar months DJF / MAM / JJA / SON within the water year
(the `Qwin/Qspr/Qsum/Qfal` convention — `fal` = Oct + Nov of the water-year start plus
Sep of its end; suffixes `win`, `spr`, `sum`, `fal`); dowy 1 = Oct 1; ties take the first
occurrence; a year is used only on a complete grid (the water year's exact day count, no
NaN after preprocessing). Every sum is a plain sequential double accumulation in all
three languages, so a cumulative sum's last element equals the total it is compared
against and the annual values are bit-identical across Julia, Python and R.

### Precipitation metrics (PPT-qualified years — the `climate_data` frame)

| Metric | Description | Units | Year NaN when |
|--------|-------------|-------|---------------|
| **prcp_annual** | Water-year total | mm | — |
| **prcp_win / prcp_spr / prcp_sum / prcp_fal** | Seasonal totals | mm | — |
| **prcp_1d_max** | Maximum daily precipitation | mm | — |
| **prcp_5d_max** | Maximum 5-consecutive-day total (window inside the water year) | mm | — |
| **prcp_max_dowy** | Day of the 1-day maximum (first on ties; an all-zero year gives day 1 — none in the input) | day | — |
| **prcp_50_dowy** | First day whose cumulative precipitation reaches 50 % of the total | day | total ≤ 0 |
| **prcp_cm_dowy** | Centre of mass Σ(dowy·P) ÷ ΣP | day | total ≤ 0 |
| **prcp_ndays_abv_10mm / _25mm / _50mm** | Days with P strictly above 10 / 25 / 50 mm | days | — |
| **prcp_annual_sd** | Sample sd (n − 1) of the daily values within the year | mm | — |
| **prcp_annual_cv** | sd ÷ mean of the daily values | – | mean ≤ 0 |
| **prcp_intermit_annual** | Longest run of dry days (P ≤ 0.1 mm); 0 when every day is wet | days | — |
| **prcp_intermit_annual_dowy_start / _dowy_end** | First / last day of the first longest dry run | day | no dry day |
| **prcp_gini** | Gini index of the daily values over all days, zeros included: G = Σᵢ(2i − n − 1)·y₍ᵢ₎ ÷ (n·Σy), y ascending (0 uniform … (n−1)/n one wet day) | – | total ≤ 0 |
| **prcp_seasonality_index** | Walsh & Lawler (1981): SI = (1/R)·Σₘ₌₁¹² \|Xₘ − R/12\| over the water year's calendar-month totals Xₘ (0 … 1.83) | – | R ≤ 0 |

### Temperature metrics (temperature-qualified years — the `temp_data` frame)

| Metric | Description | Units | Year NaN when |
|--------|-------------|-------|---------------|
| **temp_max_annual / temp_max_annual_dowy** | Maximum of daily TMAX and the first day it occurs | °C, day | — |
| **temp_min_annual / temp_min_annual_dowy** | Minimum of daily TMIN and the first day it occurs | °C, day | — |
| **temp_max_win / _spr / _sum / _fal** | Seasonal maximum of daily TMAX | °C | — |
| **temp_min_win / _spr / _sum / _fal** | Seasonal minimum of daily TMIN | °C | — |
| **frost_first_dowy** | First day with basin-mean TMIN < 0 °C | day | no frost day |
| **frost_last_dowy** | Last day with TMIN < 0 °C among Oct 1 – Jun 30 | day | no such day |
| **gdd_cum_annual** | Growing degree days Σ max(0, (TMAX + TMIN)/2 − 0 °C) | °C·day | — |
| **gdd_cum_win / _spr / _sum / _fal** | The same over the season's months | °C·day | — |

### Method and gates

1. The runner keeps six climate columns (`prcp`, `swe`, `tmin`, `tmax` + keys). The
   preprocessor qualifies the temperature pair per water year with the PPT/SWE rules
   (≤ 30 raw NaN, internal gaps ≤ 3 days linearly filled, no residual NaN; no sign
   check; both series must pass; non-finite values count as missing) and returns
   `valid_temp_years`. As for PPT and SWE, temperature is evaluated only for
   streamflow-valid years (the snow-family precedent, user decision 2026-10-12).
2. Precipitation metrics run ONLY on the explicit PPT-qualified `climate_data` frame;
   temperature metrics ONLY on the explicit `temp_data` frame (never on columns inside
   the gage frame). Both are config-gated: an absent `climate` section emits no climate
   columns (the additivity switch).
3. Standard gates: 60 % / 80 % trend completeness, the 20-value stats floor, the 16
   statistics via `generate_stats()`. The two frost dates, NaN in frost-free years, carry
   the snow family's record-anchored decade gate (anchored to the temperature-valid
   years): their six trend statistics are suppressed unless ≥ 80 % of the anchor years
   in both the first and the last decade have a frost.
4. Column names derive from the config (`prcp_5d_max`, `prcp_50_dowy`,
   `prcp_ndays_abv_<t>mm`), the `D{p}_day` pattern; the config validation rejects a
   window outside 2..366 and a cumulative fraction that is not a whole percentage.

### Data quality and caveats

- **The interpolated Dec 31.** Daymet's 365-day calendar omits Dec 31 of every leap
  CALENDAR year, so the water year that follows one (WY 1981, 1985, …, 2025: 12 of the 45
  climate years of the WY 1980–2025 product, 9 of 33 for WY 1993–2025) has one day that was
  never observed: the preprocessor fills it linearly from Dec 30 and Jan 1, as it does for
  PPT and SWE in every family. That day enters every precipitation and temperature metric
  of that year (totals, `prcp_win`, the 5-day window, wet-day counts, dry runs, the 50 %
  and centre-of-mass dates, sd, Gini, SI, `temp_*_win`, `gdd_cum_win/annual`). Size (40
  batch-1 gages, 480 such gage-years, Sonnet review 2026-10-12): mean 0.25 % of the annual
  total, median 0.1 %, max 2.1 %; above 10 mm in 5 % of those gage-years, so
  `prcp_ndays_abv_10mm` can count one unobserved day.
- **Days absent from the streamflow record (found 2026-10-16; FIXED in the version-3
  products, built 2026-10-27).** Until then the benchmark runners joined each gage's Daymet
  series onto its streamflow ROWS, so a date with no streamflow row at all (not a row with
  NA) reached the preprocessor without PPT either and was linearly interpolated like a
  climate gap (≤ 3 days; a longer run rejects the year for Q anyway). Gages drop out
  preferentially in storms, so that year's totals, wet-day counts and 5-day maximum ran low.
  Measured on the WY 1993–2025 version-2 product: 919 of 200,832 `prcp_annual` gage-years
  (0.46 %), median 6 mm, max −606 mm; 8 gages' `prcp_annual_mean` off by more than 1 %
  (02146285 −24.5 %). TMIN/TMAX and SWE behaved the same. Since the version-3 products the
  runners build each gage's frame through `docs/benchmarks/climate_frame.{jl,py,R}` under
  `STREAMFLOW_CLIMATE_FRAME=outer` (the default since 2026-10-28; `left` reproduces the
  earlier join): Daymet's own calendar inside the record span and the gage's water years, so
  annual totals equal the daily sums. Every product before version 3 carries the defect —
  CHANGELOG → Known Issues; `docs/plans/2026-10-16-precipitation-review.md`.
- **Non-finite inputs** (NaN, ±Inf) are missing observations in preprocessing for every
  series, and the climate kernels require finite values (a year with a residual
  non-finite value is NaN in every metric).
- The climate input starts on 1980-01-01, so WY 1980 has no climate metrics in the
  WY 1980–2025 product.
- Values are basin means: a dry day is a basin-mean day ≤ 0.1 mm, and frost is a
  basin-mean minimum below 0 °C, not "any part of the basin".
- **The 50 % crossing carries a 1e-12 relative tolerance** (`cumulative ≥ 0.5 × total −
  1e-12 × total`): two mathematically equal rainfall blocks are tied at 50 %, but their
  sequential sums differ in the last bits and the strict test moved the date by hundreds
  of days in such synthetic cases. 1e-12 of an annual total is physically nothing, and the
  closest real margin in the 360,585 production water years is 1.3e-9, so no production
  value depends on it. `D50_day` keeps its strict comparison (its own precedent).
- Review record: `docs/plans/2026-10-12-climate-reviews/` (blind Codex and Sonnet
  reviews per language; the Julia round-1 findings and their fixes are in CHANGELOG
  [October 2026]).

### References

- Walsh, R.P.D. & Lawler, D.M. (1981). Rainfall seasonality: description, spatial patterns and change through time. *Weather*, 36(7), 201–208 (via Sumner, G. (1988). *Precipitation: Process and Analysis*).
- Sahbeni, G. et al. (2023). Rainfall Gini index as a drought indicator. *Atmospheric Science Letters* (the co-authors' citation for `prcp_gini`).

---

## Summary Table

| Category | Function family | Function | Requires Climate | Notes |
|----------|-----------------|----------|------------------|-------|
| Flow Volume | Flow Volumes | `calculate_flow_vols_by_year` | No | 21 metrics (5 totals + 15 percentiles + Q95_Q10) + 4 season exclusion diagnostics |
| Flow Volume | Negative Flow Days | `calculate_negative_days` | No | 1 metric (negative_ann) |
| Flow Volume | Pulse Metrics (TQmean only) | `calculate_pulse_metrics` | No | TQmean (the rest of the function's output is Flashiness) |
| Flow Duration | FDC | `analyze_fdc_trends` | No | 3 metrics (FDCall, FDC90th, FDCmid) |
| Storage | Baseflow (Fixed) | `analyze_baseflow_indices` | No | 2 metrics |
| Storage | Baseflow (Recession-Parameterized) | `analyze_baseflow_indices_with_parameters` | No | 2 metrics + 1 scalar |
| Storage | Recession | `analyze_recession_parameters` | No | 7 metrics + 6 seasonality (the per-gage recession alpha scalar is counted under Baseflow (Recession-Parameterized)) |
| Storage | Q-PPT Relationships | `analyze_Q_PPT_relationships` | Yes | 5 metrics + 1 diagnostic |
| Storage | Average Storage | `calculate_average_storage` | Yes | 1 metric (computed; omitted from major analyses — no ET term) |
| Flashiness | Pulse Metrics | `calculate_pulse_metrics` | No | 13 metrics (8 pulse + 5 flow reversals); TQmean is counted under Flow Volume |
| Flashiness | Flashiness | `analyze_flashiness_trends` | No | 1 metric |
| Flashiness | Elasticity | `calculate_streamflow_elasticity` | Yes | 3 metrics (rolling + annual + static) + 2 diagnostics |
| Drought | Streamflow Drought | `calculate_drought_metrics` | No | 10 metrics (2 measures × 5 severity levels) + 5 threshold scalars (July 2026; all 3 languages since Aug 2026) |
| Flow Timing | Flow Timing | `analyze_flow_timing_trends` | No | 15 metrics |
| Precipitation Streamflow | Q-P Seasonality | `calculate_qp_seasonality` | Yes | 2 metrics |
| Snow | Snow Metrics | `calculate_snow_metrics` | SWE (Daymet) | 14 metrics (July 2026; all 3 languages since Aug 2026) |
| Climate | Precipitation | `calculate_precipitation_metrics` | PPT (Daymet) | 20 metrics (October 2026; all 3 languages) |
| Climate | Temperature | `calculate_temperature_metrics` | TMIN + TMAX (Daymet) | 19 metrics (October 2026; all 3 languages) |

---

## Automated QA/QC Flags (12 columns)

Twelve boolean `flagged_for_*` columns are computed on the assembled output table by
`compute_qa_flags()` in all three languages, from the thresholds in
`config/signatures_config.json` → `qa_qc`. They never remove a gage, and a missing
input never triggers a flag. Definitions are in the guidelines document (Part 5); the
one that needs care is:

- **`flagged_for_high_na`** — more than 30 % of the gage's **signature columns** are
  NA/NaN. The denominator is every signature-derived column present in the table —
  everything except the gage metadata and the flag columns (1,621 of the 1,653 columns of
  the 2026 products; 2,245 of 2,277 with the climate family: the 16 statistics of each
  signature plus the single-value-per-gage outputs) —
  selected by name from the config manifest `qa_qc.high_na_denominator`. Families
  NA-filled for a gage (no Daymet, no SWE) DO count, which is why the flag mostly marks
  missing climate coverage rather than bad streamflow.
- **Known issue in the delivered products (cataloged 2026-09-04).** Until 2026-09-04 each
  language counted something different, and the Julia-built standard products counted
  the 16 numeric metadata columns only — so in both delivered products this one column
  is TRUE for every Canadian gage and for USGS gages lacking GAGES-II attributes (1,224
  and 1,243 gages) and carries no information about signature completeness. Decision:
  the products are not rewritten; the column will be regenerated the next time any
  portion of the data is rerun (see CHANGELOG → Known Issues and → September 2026).

## Changepoint Detection

The non-parametric Pettitt test is applied to every signature's annual time series via `generate_stats()`.

### Pettitt Test

**Function**: `pettitt_test`

Non-parametric rank-based test (Pettitt 1979) that detects a single changepoint in the mean:

- **Test statistic**: K_T = max|U_t| where U_t = U_{t-1} + V_t and V_t = Σ_{j=1}^{n} sgn(x_t - x_j)
- **P-value**: Asymptotic approximation p ≈ 2·exp(-6K²/(T³+T²))
- **Always identifies a location**: Pettitt always returns a cp_year; use p-value for significance
- **Complexity**: O(n²) per signature — fast for annual series (n ≤ 45)

### Differential Metrics

**Function**: `segment_differential_metrics`

Applied at the Pettitt changepoint location:

- **delta_mean**: post_mean − pre_mean
- **pct_change**: delta_mean / |pre_mean| × 100 (NaN if |pre_mean| < 1e-10)
- **pre_mk_pval**: Mann-Kendall p-value for the pre-changepoint segment
- **post_mk_pval**: Mann-Kendall p-value for the post-changepoint segment

### Output Columns (per signature)

Each signature gains 8 changepoint columns:

| Suffix | Description |
|--------|-------------|
| `_pettitt_cp_year` | Most likely changepoint year |
| `_pettitt_pval` | Asymptotic p-value (< 0.05 = significant) |
| `_pettitt_pre_mean` | Mean before Pettitt changepoint |
| `_pettitt_post_mean` | Mean after Pettitt changepoint |
| `_pettitt_delta_mean` | Post minus pre mean at Pettitt changepoint |
| `_pettitt_pct_change` | Percent change at Pettitt changepoint |
| `_pettitt_pre_mk_pval` | Mann-Kendall p-value for pre-Pettitt segment |
| `_pettitt_post_mk_pval` | Mann-Kendall p-value for post-Pettitt segment |

### Configuration

From `config/signatures_config.json` → `changepoint` section:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `enabled` | `true` | Enable/disable changepoint analysis |
| `start_water_year` | 1980 | Start of analysis window |
| `end_water_year` | 2024 | End of analysis window |
| `min_total_obs` | 20 | Minimum non-NA observations required |
| `min_segment_obs` | 10 | Minimum non-NA observations per segment |

### Design Decisions

- **Scope**: Applied to ALL signatures producing annual time series — redundancy across ~76 base signatures (100 from July 2026: + 14 snow metrics + 10 drought metrics) serves as a pseudo-robustness test
- **Independence**: Changepoint analysis runs independently of the trend completeness gate
- **Non-parametric**: No distributional assumptions; robust to outliers
- **Per-signature, not per-gage**: Each signature gets its own changepoint analysis — different signatures may identify changepoints at different years

### Interpretation

- p < 0.05: Significant changepoint detected
- p ≥ 0.05: No significant changepoint (cp_year still reported but not reliable)

### Signal Robustness & Known Limitations

**Overall significance rate**: ~13.4% of evaluations have p < 0.05 (vs 5% expected under null). The excess is concentrated in physically interpretable function families (the April 2026 analysis grouping, finer than the eight signature categories):

| Function family | Sig. Rate | Notes |
|----------|-----------|-------|
| Flow Timing | 3.7% | Below null — timing is stationary; useful calibration anchor |
| Q-P Seasonality | 7.7% | Near null — weak signal |
| Storage | 6.5% | Near null |
| Flow Volumes | 12.4% | Modest excess |
| Baseflow | 15.1% | Real signal likely |
| Flow Percentiles | 18.4% | Strong signal |
| Flashiness | 19.4% | Strong signal |
| Pulses | 15.9% | Moderate signal |
| Elasticity | 25.8% | Inflated — elasticity_rolling at ~49% due to short series from 11-year rolling window |

**Redundancy across signatures**: The effective number of independent tests is ~17 out of 76 signatures (~77% redundancy). Flow percentiles (Q5–Q95) are highly correlated with Qann: P(Q30 sig | Qann sig) ≈ 75%. Users should focus on category-level summaries rather than treating each signature independently.

**Multiple testing**: With 76 tests per gage at alpha=0.05, the expected false positive count is 3.8 per gage and the family-wise error rate is ~98%. No correction is applied in the output. After Benjamini-Hochberg FDR at 0.05, approximately 3.5% of evaluations survive — a core of robust detections. Users performing per-gage analysis should apply their own multiple testing correction.

**CP year center bias**: The 10-year minimum segment constraint limits cp_year to [1989, 2014], with maximum test power near the center of the window. The observed concentration around 2000–2003 (~50% of significant detections in 1998–2006) is partly a geometric artifact and partly a real hydroclimatic signal (North American drought/pluvial transitions).

**Changepoint clustering within gages**: Many signatures detect the same cp_year for a given gage (median 7 significant signatures per gage, ~40% sharing the same modal year). This is expected — a real shift in Qann propagates to flow percentiles, baseflow, flashiness — but means the 608 columns overstate the independent information content.

**elasticity_rolling caveat**: The 11-year rolling window produces n−10 annual values from n input years, so gages with 25 qualifying years yield only 15 elasticity values. The Pettitt test has inflated false positive rates on these short series. The near-uniform CP year distribution for this signature (vs concentrated for other signatures) is consistent with noise.

**Limitations of the Pettitt test**:
- Detects shifts in central tendency only — misses gradual trends, variance changes, and seasonal pattern shifts
- The asymptotic p-value is conservative for small n (significance rate is 4.3% for n=20–25, below the 5% null expectation)
- Does not distinguish step changes from trend breaks (the BIC 4-model comparison, retained in `changepoint.jl` but not called, was designed for this)

**Differential metrics**: The pct_change metric discriminates well — median |pct_change| is ~40% for significant detections vs ~15% for non-significant (2.65x ratio). The pre/post MK p-values are largely uninformative: ~92% of segments show no significant within-segment trend, confirming the step-change model is appropriate but adding little actionable information.

### References

- Pettitt, A.N. (1979). A non-parametric approach to the change-point problem. Applied Statistics, 28(2), 126-135.

---

## Implementation Design Notes

The following design decisions affect how signatures are calculated and which gages produce valid results.

### Recession Analysis Requirements

The recession analysis (`analyze_recession_parameters`) requires a **minimum of 25 recession events** across the entire record. This is a conservative threshold that ensures robust statistical estimation but may cause many gages in arid or flashy climates to produce NA values for recession metrics.

**Impact**: Gages with sparse recession events (e.g., intermittent streams, heavily regulated rivers) will have all recession metrics set to NA.

### Centralized NA Handling (April 2026)

Missing data is now handled centrally by `preprocess_daily_data()`, called once per gage before any signature computation. This replaces the previous per-signature ad-hoc handling:

1. **Daily grid normalization**: One row per calendar day per water year, sorted, duplicates removed
2. **Interpolation**: Internal gaps of <=3 consecutive days are linearly interpolated. Leading/trailing (boundary) NAs are NOT interpolated.
3. **Year rejection**: Years with >30 raw NAs, >3-day gaps, or residual boundary NAs are excluded. Negative Q rejection is conditional on `reject_negative_flow` config (default: false)
4. **Seasonal completeness**: Computed from raw (pre-interpolation) observation counts. Seasons below 80% completeness are flagged; affected seasonal metrics are set to NA
5. **Climate NA policy**: PPT/SWE/temp have their own NA limits, tracked separately. A year can be valid for Q-only signatures but invalid for Q-PPT signatures.
6. **Constant-SD flag**: Detects months with a single unique non-zero Q value (sensor flatline). This is a QA flag, not a rejection criterion.
7. **Trend completeness**: `generate_stats()` optionally requires >=60% non-NA annual metrics overall (lowered from 80% in July 2026 per updated guidelines) and >=80% in first/last decade before computing trend statistics. Short series (<10 years) skip the decade check.

**Configuration**: `config/signatures_config.json` → `na_handling` section. `use_legacy_filtering: true` preserves the old 95%-non-NA-days rule.

**Previous behavior** (preserved when `use_legacy_filtering: true`): The flow timing function replaced missing daily flow values with zero before computing cumulative sums. This was a conservative approach that prevented NA propagation through cumsum. The new preprocessing eliminates this need by ensuring valid years have no NAs.

### Trend Completeness Exemptions

The trend completeness gate (requiring >=60% non-NA annual values overall — lowered from 80% in July 2026 — and >=80% in first/last decades) is **not applied** to two signature categories:

- **Recession**: Event-based and inherently sparse. Many years may have no recession events, so annual value counts naturally fall well below 80%. Applying the gate would eliminate most gages.
- **Elasticity**: Rolling-window method produces N-10 values from N input years, inherently falling below the 80% threshold relative to the full record length.

These exemptions are implemented consistently in Python, Julia, and rpkg orchestration functions.

### Seasonality Model Periodicity

The recession seasonality analysis fits a sinusoidal model with a **hardcoded 365-day period**. This does not account for leap years (366 days), which may introduce minor phase shifts in multi-decadal records.

**Impact**: Minimal for most analyses; consider using day-of-year instead of calendar day for improved precision in future versions.

### Baseflow Index Relationship

The expected relationship **BFI_Eckhardt < BFI_LyneHollick** is validated in QA/QC but not enforced during calculation. Some gages may violate this relationship due to:
- Unusual flow regimes
- Data quality issues
- Edge cases in filter initialization

**Recommendation**: Flag gages where this relationship is violated for manual review.

### Recession-Informed BFI (Implemented April 2026)

The guidelines describe `analyze_baseflow_indices_with_parameters()` — calculating BFI using recession parameters as inputs to the Eckhardt and Lyne-Hollick filters. Implemented in all three languages (Julia canonical April 24, Python and rpkg ported April 25) using the discrete recession constant `alpha = Q_{i+1}/Q_i` under the linear reservoir assumption (b=1). BFImax remains fixed at 0.8; only the recession constant `a` is parameterized. The Lyne-Hollick parameterization is heuristic. See Section 3 (Baseflow) for details.

**Known limitation — narrow BFI_Eckhardt_param range**: The parameterized Eckhardt BFI has low discriminating power across gages (range [0.47, 0.80], std=0.036) compared to the fixed-parameter version (range [0.14, 0.80], std=0.119). This is a mathematical property of the Eckhardt filter, not a bug. The `a` parameter controls filter dynamics (how fast baseflow tracks total flow), while `BFImax` controls the steady-state BFI level. Since 92% of gages have recession alpha < 0.95, the filter converges rapidly to the BFImax=0.8 ceiling. With the fixed `a=0.98`, the filter's high inertia causes BFI to deviate significantly from BFImax during storm events, producing wider variation. The parameterized L-H filter does not have this issue (range [0.25, 0.97]) because the alpha parameter plays a structurally different role.

### BFImax Estimation via Backward Filter (Not Yet Implemented)

To improve BFI_Eckhardt_param discriminating power, BFImax should also be parameterized per gage. Collischonn & Fan (2013) proposed a backward filter method:

1. **Estimate recession constant `a`** — already available as `recession_alpha_point_cloud_linear_reservoir`
2. **Run backward filter** — exploits the linear reservoir equation in reverse (`b_{k-1} = b_k / a`) to reconstruct the maximum possible baseflow at each timestep:
   ```
   Initialize: b_N = Q_N
   For i = N-1 down to 1:
       b_i = b_{i+1} / a
       if b_i > Q_i: b_i = Q_i   (baseflow cannot exceed total flow)
   ```
3. **Compute BFImax** — `BFImax = sum(b) / sum(Q)`, or max of annual BFI values with a 0.9 safeguard

This would replace the fixed BFImax=0.8 with a per-gage value (expected range ~0.3–0.9 depending on geology/hydrogeology), giving the Eckhardt filter real per-gage variation.

**Known limitations of the backward filter** (Zhang et al. 2021, HESS):
- Overestimates BFImax in humid catchments (backward baseflow approaches total flow)
- Single recession constant assumption (no seasonal/non-stationary variation)
- Produces an upper-bound estimate, not necessarily actual baseflow
- Reference implementation (`xiejx5/baseflow` Python package) uses ad-hoc cap at BFImax >= 0.9

**References**: Collischonn & Fan (2013), Hydrological Processes 27(18):2614-2622; Zhang et al. (2021), HESS 25:1747.

### Un-Normalized Canadian Gages (`area_normalized = FALSE`)

37 Canadian gages in the signature output have **no drainage area in HYDAT** (neither
gross nor effective — verified against `Hydat.sqlite3`, July 2026), so their flow is
in **raw m³/s instead of mm/day** (`area_normalized = FALSE`). Most are not natural
watersheds: irrigation/diversion canals, dam/powerhouse outflows, and huge-river
channel splits (St. Lawrence, Mackenzie, Nelson) where drainage area is undefined or
unpublished. Decision (July 2026): these gages are retained with raw flow — no area
backfill.

**Q-to-PPT gate (July 2026)**: because Q (m³/s) and PPT (mm) units don't match for
these gages, ALL Q-to-PPT signatures — runoff ratios (+ `runoff_ratio_high_count`),
elasticity (+ diagnostics), Q-P seasonality, and `avg_storage` — are **skipped and
emit NA by design**. `calculate_all_signatures()` takes an `area_normalized`
argument (all three languages); the benchmark runners set it from the metadata.
Elasticity and `qp_bimodality` are technically scale-invariant in Q but are gated
with the rest of the family (mixed-unit inputs are conceptually invalid, and the
internal PPT thresholds and P − Q terms are not scale-invariant).

**Impact on remaining signatures**: unit-carrying Q-only signatures (Qann/seasonal
totals, flow percentiles, Q95_Q10, log_a) stay in raw m³/s — incomparable with
mm/day gages (Qann_mean reaches 3.18M for the St. Lawrence). Dimensionless Q-only
signatures (BFI, FDC slopes, timing, TQmean, recession b, pulse counts, flashiness)
remain valid.

**Flag gap**: `flagged_for_qann_range` catches only 27 of the 37; 10 small
canals/creeks land inside [0, 2000] unflagged. **Filter on `area_normalized == TRUE`
before cross-gage comparison of unit-carrying signatures.** Details:
`docs/DEVELOPMENT.md` → Canadian HYDAT → Missing drainage areas; tracked in
`CHANGELOG.md` → Known Issues.

### Elasticity 30% Diagnostic (Pending)

The guidelines request a counterfactual diagnostic: "what the number of excluded years would be if the value was <30% data missing." This is **pending clarification with domain experts**. Current diagnostics (`elasticity_years_total`, `elasticity_years_low_ppt`) track actual exclusions but not this hypothetical threshold.

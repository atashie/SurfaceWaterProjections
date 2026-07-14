# Streamflow Signatures Reference

Detailed documentation of all hydrological signatures calculated by this project.

> **Comprehensive Documentation**: See `SIGNATURE_GUIDELINES.md` for mathematical formulations, parameter choices, and scientific rationale. This file is auto-synced from the [collaborative Google Doc](https://docs.google.com/document/d/e/2PACX-1vQnt7OCPm19vnWF4yynXL9JTzTvq9CrGoEaDv7yFSngLoFsypiWsx6fZLKWwaO5YQ/pub).

## Overview

Each signature produces **8 statistics** via `generate_stats()`:

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

Analyzes recession curve behavior using dQ/dt = a*Q^b relationship.

### Metrics

| Metric | Description |
|--------|-------------|
| **log_a_pointcloud** | Recession rate parameter (point cloud method) |
| **log_a_events** | Recession rate parameter (event-based method) |
| **b_pointcloud** | Recession exponent (point cloud method) |
| **b_events** | Recession exponent (event-based method) |
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
- Seasonality metrics are single values (exceptions to 8-statistic rule)
- Documented in `config.R` as `EXPECTED_RECESSION_SEASONALITY`
- Requires minimum 25 recession events (`RECESSION_MIN_EVENTS`) for all metrics except `n_recession_events`
- `n_recession_events` is computed independently of the min_events gate (useful for gages with few events)
- `alpha_linear` per-year values require >10 valid alpha pairs per year; same 25-event gate as other recession metrics for trend statistics
- `recession_alpha_point_cloud_linear_reservoir` is a per-gage scalar (whole-record median of Q_{i+1}/Q_i ratios across all recession events). Computed independently of the 25-event gate (only requires >10 alpha pairs). Used to parameterize `BFI_Eckhardt_param` and `BFI_LyneHollick_param`.
- Recession events are identified as contiguous windows where both Q and |dQ/dt| are monotonically decreasing, with a minimum length of 5 days

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

## Summary Table

| Category | Function | Requires Climate | Notes |
|----------|----------|------------------|-------|
| Flow Volumes | `calculate_flow_vols_by_year` | No | 22 metrics (5 totals + 16 percentiles + Q95-Q10) + 4 season exclusion diagnostics |
| FDC | `analyze_fdc_trends` | No | 3 metrics (FDCall, FDC90th, FDCmid) |
| Baseflow (Fixed) | `analyze_baseflow_indices` | No | 2 metrics |
| Baseflow (Recession-Parameterized) | `analyze_baseflow_indices_with_parameters` | No | 2 metrics + 1 scalar |
| Recession | `analyze_recession_parameters` | No | 8 metrics + 6 seasonality |
| Pulse Metrics | `calculate_pulse_metrics` | No | 15 metrics (8 pulse + TQmean + 5 flow reversals + negative_ann) |
| Flashiness | `analyze_flashiness_trends` | No | 1 metric |
| Flow Timing | `analyze_flow_timing_trends` | No | 15 metrics |
| Q-PPT Relationships | `analyze_Q_PPT_relationships` | Yes | 5 metrics + 1 diagnostic |
| Elasticity | `calculate_streamflow_elasticity` | Yes | 3 metrics (rolling + annual + static) + 2 diagnostics |
| Q-P Seasonality | `calculate_qp_seasonality` | Yes | 2 metrics |
| Average Storage | `calculate_average_storage` | Yes | 1 metric |
| Negative Flow Days | `calculate_negative_days` | No | 1 metric (negative_ann) |

---

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

- **Scope**: Applied to ALL signatures producing annual time series — redundancy across ~76 base signatures serves as a pseudo-robustness test
- **Independence**: Changepoint analysis runs independently of the 80% trend completeness gate
- **Non-parametric**: No distributional assumptions; robust to outliers
- **Per-signature, not per-gage**: Each signature gets its own changepoint analysis — different signatures may identify changepoints at different years

### Interpretation

- p < 0.05: Significant changepoint detected
- p ≥ 0.05: No significant changepoint (cp_year still reported but not reliable)

### Signal Robustness & Known Limitations

**Overall significance rate**: ~13.4% of evaluations have p < 0.05 (vs 5% expected under null). The excess is concentrated in physically interpretable categories:

| Category | Sig. Rate | Notes |
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
7. **Trend completeness**: `generate_stats()` optionally requires >=80% non-NA annual metrics and >=80% in first/last decade before computing trend statistics. Short series (<10 years) skip the decade check.

**Configuration**: `config/signatures_config.json` → `na_handling` section. `use_legacy_filtering: true` preserves the old 95%-non-NA-days rule.

**Previous behavior** (preserved when `use_legacy_filtering: true`): The flow timing function replaced missing daily flow values with zero before computing cumulative sums. This was a conservative approach that prevented NA propagation through cumsum. The new preprocessing eliminates this need by ensuring valid years have no NAs.

### Trend Completeness Exemptions

The 80% trend completeness gate (requiring >=80% non-NA annual values and >=80% in first/last decades) is **not applied** to two signature categories:

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

### Elasticity 30% Diagnostic (Pending)

The guidelines request a counterfactual diagnostic: "what the number of excluded years would be if the value was <30% data missing." This is **pending clarification with domain experts**. Current diagnostics (`elasticity_years_total`, `elasticity_years_low_ppt`) track actual exclusions but not this hypothetical threshold.

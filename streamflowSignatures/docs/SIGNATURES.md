# Streamflow Signatures Reference

Detailed documentation of all hydrological signatures calculated by this project.

> **Comprehensive Documentation**: See `SIGNATURE_GUIDELINES.md` for mathematical formulations, parameter choices, and scientific rationale. This file is auto-synced from the [collaborative Google Doc](https://docs.google.com/document/u/1/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub).

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
| **Q95-Q10** | Ratio of high to low flow percentiles |

### Units
- Qann, Qwin, Qspr, Qsum, Qfal: mm (total over period, summed from daily mm/day values)
- Q percentiles: mm/day (daily flow at each percentile)

### Data Quality
Years with fewer than 250 non-NA days are excluded.

---

## 2. Baseflow

**Function**: `analyze_baseflow_indices`

### Metrics

| Metric | Description | Parameters |
|--------|-------------|------------|
| **BFI_Eckhardt** | Baseflow index using Eckhardt recursive digital filter | BFImax=0.8, a=0.98 |
| **BFI_LyneHollick** | Baseflow index using Lyne-Hollick filter | alpha=0.925, 2 passes |

### Expected Relationship
BFI_Eckhardt < BFI_LyneHollick (validated in QA/QC)

### References
- Eckhardt, K. (2005). How to construct recursive digital filters for baseflow separation.
- Lyne, V., & Hollick, M. (1979). Stochastic time-variable rainfall-runoff modelling.

---

## 3. Recession

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
| **log_a_seasonality_amplitude_all** | Seasonal amplitude in recession rate (all data) |
| **log_a_seasonality_amplitude_first_half** | Seasonal amplitude (first half of record) |
| **log_a_seasonality_amplitude_last_half** | Seasonal amplitude (last half of record) |
| **log_a_seasonality_minimum_all** | Seasonal minimum day for recession rate (all data) |
| **log_a_seasonality_minimum_first_half** | Seasonal minimum day (first half of record) |
| **log_a_seasonality_minimum_last_half** | Seasonal minimum day (last half of record) |

### Notes
- Seasonality metrics are single values (exceptions to 8-statistic rule)
- Documented in `config.R` as `EXPECTED_RECESSION_SEASONALITY`
- Requires minimum 25 recession events (`RECESSION_MIN_EVENTS`)

---

## 4. Pulse Metrics

**Function**: `calculate_pulse_metrics`

### Metrics

| Metric | Description | Threshold |
|--------|-------------|-----------|
| **n_high_pulses_year** | Count of high pulse events | > 90th percentile |
| **n_low_pulses_year** | Count of low pulse events | < 10th percentile |
| **dur_high_pulses_year** | Mean duration of high pulses | days |
| **dur_low_pulses_year** | Mean duration of low pulses | days |
| **TQmean** | Percentage of days with flow above annual mean | % |
| **Flow_Reversals** | Direction changes in flow | annual and seasonal |

---

## 5. Flashiness

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

## 6. Flow Timing

**Function**: `analyze_flow_timing_trends`

### Metrics

| Metric | Description |
|--------|-------------|
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
| **D25_to_D75** | Duration between 25% and 75% cumulative flow (days) |
| **Dmax** | Day of maximum flow |

### Notes
- Day 1 = October 1 (start of water year)
- Day 365/366 = September 30 (end of water year)
- NAs in daily flow are treated as zero for cumulative sum calculations

---

## 7. Q-PPT Relationships

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

### PPT Thresholds
- Annual: minimum 10mm PPT required (below returns NA)
- Seasonal: minimum 1mm PPT required (below returns NA)

---

## 8. Streamflow Elasticity

**Function**: `calculate_streamflow_elasticity`

**Requires Climate Data**

### Metrics

| Metric | Description |
|--------|-------------|
| **elasticity_static** | Overall catchment elasticity (single value) |
| **elasticity** | Rolling window (11-year) elasticity trend |

### Formula
```
E = (dQ/dP) / (Q_mean/P_mean)
```

### Interpretation
- E ~ 1.0: Proportional response to precipitation changes
- E > 1.0: Amplified response (Q changes more than P)
- E < 1.0: Dampened response

### Notes
- `elasticity_static` is a single value (exception to 8-statistic rule)
- Documented in `config.R` as `EXPECTED_ELASTICITY_STATIC`

### Reference
Sawicz, K., et al. (2011). Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA.

---

## 9. Q-P Seasonality

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

## 10. Average Storage

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
| Flow Volumes | `calculate_flow_vols_by_year` | No | 22 metrics (5 totals + 16 percentiles + Q95-Q10) |
| FDC | `analyze_fdc_trends_from_streamflow` | No | 3 metrics (FDCall, FDC90th, FDCmid) |
| Baseflow | `analyze_baseflow_indices` | No | 2 metrics |
| Recession | `analyze_recession_parameters` | No | 5 metrics + 6 seasonality |
| Pulse Metrics | `calculate_pulse_metrics` | No | 15 metrics |
| Flashiness | `analyze_flashiness_trends` | No | 1 metric |
| Flow Timing | `analyze_flow_timing_trends` | No | 13 metrics |
| Q-PPT Relationships | `analyze_Q_PPT_relationships` | Yes | 5 metrics |
| Elasticity | `calculate_streamflow_elasticity` | Yes | 1 metric + 1 static |
| Q-P Seasonality | `calculate_qp_seasonality` | Yes | 2 metrics |
| Average Storage | `calculate_average_storage` | Yes | 1 metric |

---

## Implementation Design Notes

The following design decisions affect how signatures are calculated and which gages produce valid results.

### Recession Analysis Requirements

The recession analysis (`analyze_recession_parameters`) requires a **minimum of 25 recession events** across the entire record. This is a conservative threshold that ensures robust statistical estimation but may cause many gages in arid or flashy climates to produce NA values for recession metrics.

**Impact**: Gages with sparse recession events (e.g., intermittent streams, heavily regulated rivers) will have all recession metrics set to NA.

### Flow Timing NA Handling

The flow timing function (`analyze_flow_timing_trends`) replaces missing daily flow values with **zero before computing cumulative sums**. This is a conservative approach:
- Prevents a single NA from corrupting the entire year's D-day calculations
- May slightly delay calculated timing metrics (D50, D90, etc.) in years with data gaps

**Alternative approaches**: Interpolation would provide smoother estimates but could introduce bias in gages with systematic data gaps.

### Seasonality Model Periodicity

The recession seasonality analysis fits a sinusoidal model with a **hardcoded 365-day period**. This does not account for leap years (366 days), which may introduce minor phase shifts in multi-decadal records.

**Impact**: Minimal for most analyses; consider using day-of-year instead of calendar day for improved precision in future versions.

### Baseflow Index Relationship

The expected relationship **BFI_Eckhardt < BFI_LyneHollick** is validated in QA/QC but not enforced during calculation. Some gages may violate this relationship due to:
- Unusual flow regimes
- Data quality issues
- Edge cases in filter initialization

**Recommendation**: Flag gages where this relationship is violated for manual review.

# Streamflow Signatures Reference

Detailed documentation of all hydrological signatures calculated by this project.

> **Comprehensive Documentation**: See "Summary Documentation for Streamflow Signatures.docx" for mathematical formulations, parameter choices, and scientific rationale.

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
| **Qann** | Annual mean streamflow |
| **Qwin** | Winter mean (Dec-Feb) |
| **Qspr** | Spring mean (Mar-May) |
| **Qsum** | Summer mean (Jun-Aug) |
| **Qfal** | Fall mean (Sep-Nov) |
| **Q05-Q95** | Flow percentiles (Q05, Q10, Q20, Q25, Q30, Q40, Q50, Q60, Q70, Q75, Q80, Q90, Q95) |

### Units
mm/day (after conversion from raw API units)

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
| **log_a_seasonality_amplitude** | Seasonal amplitude in recession rate |
| **log_a_seasonality_minimum** | Seasonal minimum day for recession rate |

### Notes
- Seasonality metrics are single values (exceptions to 8-statistic rule)
- Documented in `config.R` as `EXPECTED_RECESSION_SEASONALITY`

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
| **R-B Index** | Richards-Baker flashiness index |

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
| **D05_day** | Day of water year when cumulative flow reaches 5% |
| **D10_day** | Day of water year when cumulative flow reaches 10% |
| **D25_day** | Day of water year when cumulative flow reaches 25% |
| **D50_day** | Day of water year when cumulative flow reaches 50% (center of mass) |
| **D75_day** | Day of water year when cumulative flow reaches 75% |
| **D90_day** | Day of water year when cumulative flow reaches 90% |
| **D95_day** | Day of water year when cumulative flow reaches 95% |
| **D25_to_D75** | Duration between 25% and 75% cumulative flow (days) |
| **Dmax** | Day of maximum flow |

### Notes
- Day 1 = October 1 (start of water year)
- Day 365/366 = September 30 (end of water year)

---

## 7. Q-PPT Relationships

**Function**: `analyze_Q_PPT_relationships`

**Requires Climate Data** (PPT column via Daymet integration or Caravan)

### Metrics

| Metric | Description |
|--------|-------------|
| **RR_ann** | Annual runoff ratio (Q/P) |
| **RR_win** | Winter runoff ratio |
| **RR_spr** | Spring runoff ratio |
| **RR_sum** | Summer runoff ratio |
| **RR_fal** | Fall runoff ratio |

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
| Flow Volumes | `calculate_flow_vols_by_year` | No | 18 metrics |
| Baseflow | `analyze_baseflow_indices` | No | 2 metrics |
| Recession | `analyze_recession_parameters` | No | 5 metrics + 2 seasonality |
| Pulse Metrics | `calculate_pulse_metrics` | No | 6 metrics |
| Flashiness | `analyze_flashiness_trends` | No | 1 metric |
| Flow Timing | `analyze_flow_timing_trends` | No | 9 metrics |
| Q-PPT Relationships | `analyze_Q_PPT_relationships` | Yes | 5 metrics |
| Elasticity | `calculate_streamflow_elasticity` | Yes | 1 metric + 1 static |
| Q-P Seasonality | `calculate_qp_seasonality` | Yes | 2 metrics |
| Average Storage | `calculate_average_storage` | Yes | 1 metric |

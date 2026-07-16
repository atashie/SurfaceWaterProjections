# Experiment 'prod14jul_1993_60pct' vs Julia Baseline: Comparison Report

**Generated**: 2026-07-15 08:40:45

## Experiment Description

Water years >= 1993 AND at least 60% of possible years (WY1993 to gage max) must have qualifying data.

## Input Files

| Dataset | File | Gages | Columns |
|---------|------|-------|---------|
| Julia Baseline | `startIn1993_60pct_signatures.csv` | 6,579 | 631 |
| Experiment (prod14jul_1993_60pct) | `streamflow_1993_60pct_14jul2026_signatures.csv` | 6,579 | 1488 |

**Common gages**: 6,579
**Dropped gages** (in baseline, not in experiment): 0
**Added gages** (in experiment, not in baseline): 0

### Years Per Gage (Common Gages)

```
Years per gage (common gages):
    Baseline:   mean=30.2, median=32.0
    Experiment: mean=30.2, median=32.0
    Mean diff:  0.0 years
```

## High-Level Alignment Summary

### Distribution Statistics

| Metric | Value |
|--------|-------|
| Columns compared | 594 |
| Mean R2 (identity) | -3.105844 |
| Median R2 | 1.000000 |
| SD of R2 | 50.858593 |
| Min R2 | -827.5216 |

### Agreement Tiers

| Tier | Threshold | Count | % |
|------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 572 | 96.3% |
| Good | 0.99 <= R2 < 0.999 | 0 | 0.0% |
| Poor | 0.95 <= R2 < 0.99 | 0 | 0.0% |
| Low | 0.90 <= R2 < 0.95 | 0 | 0.0% |
| Very Low | 0.50 <= R2 < 0.90 | 0 | 0.0% |
| Extremely Low | R2 < 0.50 | 22 | 3.7% |
| **Total** | | **594** | **100%** |

## Agreement by Signature Category

| Category | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| Baseflow | 16 | 16 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Elasticity | 19 | 19 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| FDC | 24 | 24 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Flashiness | 8 | 8 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Flow Percentiles | 128 | 128 | 0 | 0 | 0 | 0 | 0 | 0.999999 | 1.0000 |
| Flow Timing | 120 | 120 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Flow Volumes | 40 | 40 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Negative Flow | 8 | 8 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Pulse Metrics | 112 | 112 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Q-P Seasonality | 16 | 16 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Recession | 54 | 32 | 0 | 0 | 0 | 0 | 22 | -44.164286 | -827.5216 |
| Runoff Ratios | 41 | 41 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |
| Storage | 8 | 8 | 0 | 0 | 0 | 0 | 0 | 1.000000 | 1.0000 |

## Agreement by Statistic Type

| Stat Type | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|-----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| mean | 73 | 71 | 0 | 0 | 0 | 0 | 2 | 0.288257 | -48.4408 |
| median | 73 | 71 | 0 | 0 | 0 | 0 | 2 | 0.630124 | -23.9176 |
| senn_slp | 73 | 71 | 0 | 0 | 0 | 0 | 2 | 0.222606 | -47.8463 |
| linear_slp | 73 | 71 | 0 | 0 | 0 | 0 | 2 | -1.514944 | -176.2219 |
| spearman_rho | 73 | 71 | 0 | 0 | 0 | 0 | 2 | 0.948722 | -1.4011 |
| spearman_pval | 73 | 71 | 0 | 0 | 0 | 0 | 2 | 0.953363 | -0.7273 |
| mk_rho | 73 | 71 | 0 | 0 | 0 | 0 | 2 | 0.948484 | -1.3905 |
| mk_pval | 73 | 71 | 0 | 0 | 0 | 0 | 2 | 0.955476 | -0.7280 |
| scalar | 10 | 4 | 0 | 0 | 0 | 0 | 6 | -209.541399 | -827.5216 |

## Columns with R2 < 0.99 (22 columns)

| Column | Category | Stat | R2 | Spearman | MAD | Max Diff | Baseline NAs | Experiment NAs | NA Mismatch |
|--------|----------|------|----|----------|-----|----------|--------------|----------------|-------------|
| `log_a_seasonality_amplitude_last_half` | Recession | scalar | -827.5216 | 0.4636 | 3.3559 | 141.3139 | 1218 | 1218 | 0 |
| `log_a_seasonality_amplitude_all` | Recession | scalar | -711.3948 | 0.4522 | 2.9573 | 90.8460 | 1216 | 1216 | 0 |
| `log_a_seasonality_amplitude_first_half` | Recession | scalar | -558.3776 | 0.4457 | 2.8078 | 89.5238 | 1218 | 1218 | 0 |
| `log_a_events_linear_slp` | Recession | linear_slp | -176.2219 | -0.2165 | 0.0594 | 3.3224 | 1216 | 1216 | 0 |
| `log_a_events_mean` | Recession | mean | -48.4408 | 0.2588 | 1.8804 | 47.1275 | 1216 | 1216 | 0 |
| `log_a_events_senn_slp` | Recession | senn_slp | -47.8463 | -0.1484 | 0.0368 | 1.4917 | 1216 | 1216 | 0 |
| `log_a_events_median` | Recession | median | -23.9176 | 0.3808 | 1.3858 | 28.0682 | 1216 | 1216 | 0 |
| `log_a_pointcloud_senn_slp` | Recession | senn_slp | -6.9034 | 0.3186 | 0.0166 | 2.0509 | 1342 | 1342 | 0 |
| `log_a_pointcloud_linear_slp` | Recession | linear_slp | -5.3690 | 0.2802 | 0.0171 | 0.7078 | 1342 | 1342 | 0 |
| `log_a_pointcloud_mean` | Recession | mean | -1.5164 | 0.7248 | 0.4381 | 11.1600 | 1342 | 1342 | 0 |
| `log_a_events_spearman_rho` | Recession | spearman_rho | -1.4011 | -0.0870 | 0.3066 | 1.4138 | 1216 | 1216 | 0 |
| `log_a_events_mk_rho` | Recession | mk_rho | -1.3905 | -0.0888 | 0.2135 | 1.0493 | 1216 | 1216 | 0 |
| `log_a_pointcloud_median` | Recession | median | -1.0834 | 0.7622 | 0.3930 | 10.5794 | 1342 | 1342 | 0 |
| `log_a_events_mk_pval` | Recession | mk_pval | -0.7280 | 0.1386 | 0.3276 | 0.9998 | 1216 | 1216 | 0 |
| `log_a_events_spearman_pval` | Recession | spearman_pval | -0.7273 | 0.1389 | 0.3231 | 0.9975 | 1216 | 1216 | 0 |
| `log_a_seasonality_minimum_last_half` | Recession | scalar | -0.7213 | -0.1870 | 113.9424 | 342.7690 | 1218 | 1218 | 0 |
| `log_a_seasonality_minimum_all` | Recession | scalar | -0.7211 | -0.2072 | 111.0826 | 334.3436 | 1216 | 1216 | 0 |
| `log_a_seasonality_minimum_first_half` | Recession | scalar | -0.6775 | -0.1715 | 105.8926 | 344.5595 | 1218 | 1218 | 0 |
| `log_a_pointcloud_spearman_pval` | Recession | spearman_pval | -0.6772 | 0.1831 | 0.3181 | 0.9950 | 1342 | 1342 | 0 |
| `log_a_pointcloud_mk_pval` | Recession | mk_pval | -0.5221 | 0.2509 | 0.3162 | 0.9956 | 1342 | 1342 | 0 |
| `log_a_pointcloud_mk_rho` | Recession | mk_rho | -0.3701 | 0.3583 | 0.1995 | 2.0000 | 1342 | 1342 | 0 |
| `log_a_pointcloud_spearman_rho` | Recession | spearman_rho | -0.3421 | 0.3696 | 0.2701 | 2.0000 | 1342 | 1342 | 0 |

## NA Mismatch Analysis

Columns where the number of NAs differs by >50 gages.

No columns with >50 NA mismatches.

## Summary

| Agreement Tier | Threshold | Count | % |
|----------------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 572 | 96.3% |
| Good | 0.99 <= R2 < 0.999 | 0 | 0.0% |
| Poor | 0.95 <= R2 < 0.99 | 0 | 0.0% |
| Low | 0.90 <= R2 < 0.95 | 0 | 0.0% |
| Very Low | 0.50 <= R2 < 0.90 | 0 | 0.0% |
| Extremely Low | R2 < 0.50 | 22 | 3.7% |
| **Total compared** | | **594** | **100%** |

Gages dropped by experiment filter: **0**

---
*Generated by `docs/benchmarks/compare_experiment_vs_julia.py prod14jul_1993_60pct`*

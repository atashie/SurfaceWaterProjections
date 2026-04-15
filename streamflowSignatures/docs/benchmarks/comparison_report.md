# Three-Way Cross-Language Comparison Report

**Generated**: 2026-04-14 15:18:00

## Input Files

| File | Timestamp | Gages | Columns |
|------|-----------|-------|---------|
| R (`r_signatures.csv`) | 2026-04-11 17:33:03 | 5687 | 572 |
| Python (`python_signatures.csv`) | 2026-04-14 14:38:22 | 7313 | 627 |
| Julia (`julia_signatures.csv`) | 2026-04-14 14:28:44 | 7313 | 627 |

**Common gages (all 3)**: 5,687
**Common signature columns**: 543

## Performance Comparison

| Metric | Python | Julia | Ratio (Py/Jl) |
|--------|--------|-------|---------------|
| Total time | 10531s (175.5 min) | 1650s (27.5 min) | 6.4x |
| Gages processed | 7,313 | 7,313 | - |
| Processing rate | 0.69/s | 4.43/s | 6.4x faster |
| R benchmark | ~1-2 hours (estimated) | - | - |

### Phase Breakdown

| Phase | Python | Julia |
|-------|--------|-------|
| filter_gages | 2015s | 614s |
| load_climate | 65s | 84s |
| load_streamflow | 12s | 28s |
| metadata_qaqc_save | 5s | 20s |
| process_signatures | 8435s | 903s |

## Overall Identity R² Summary

R² of the identity line (y = x): measures whether implementations produce identical values, not just correlated values.

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 | Cols < 0.95 | Cols < 0.90 |
|------|---------|-----------|--------|-------------|-------------|-------------|
| R vs Python | 0.601828 | 0.999890 | -49.7162 | 49 | 46 | 46 |
| R vs Julia | 0.602189 | 0.999887 | -49.7162 | 49 | 46 | 46 |
| Python vs Julia | 0.999867 | 1.000000 | 0.9863 | 3 | 0 | 0 |

### Spearman Rank Correlation (secondary diagnostic)

| Pair | Mean rho | Min rho | Cols < 0.99 |
|------|----------|---------|-------------|
| R vs Python | 0.944434 | 0.0421 | 46 |
| R vs Julia | 0.944432 | 0.0423 | 46 |
| Python vs Julia | 0.999976 | 0.9975 | 0 |

## Agreement by Signature Category

| Category | Total Cols | Perfect (>=0.999) | Good (>=0.99) | Poor (<0.99) | Min R² |
|----------|-----------|-------------------|---------------|-------------|--------|
| Baseflow | 16 | 16 | 0 | 0 | 0.9995 |
| Elasticity | 1 | 0 | 1 | 0 | 0.9987 |
| FDC | 24 | 20 | 2 | 2 | 0.9792 |
| Flashiness | 8 | 8 | 0 | 0 | 0.9997 |
| Flow Percentiles | 128 | 128 | 0 | 0 | 0.9991 |
| Flow Timing | 104 | 98 | 6 | 0 | 0.9979 |
| Flow Volumes | 40 | 40 | 0 | 0 | 0.9992 |
| Pulse Metrics | 112 | 104 | 7 | 1 | 0.9882 |
| Q-P Seasonality | 16 | 16 | 0 | 0 | 0.9991 |
| Recession | 46 | 0 | 0 | 46 | -49.7162 |
| Runoff Ratios | 40 | 40 | 0 | 0 | 0.9990 |
| Storage | 8 | 6 | 2 | 0 | 0.9969 |

## Columns with min Identity R² < 0.99 (49 columns)

| Column | Category | R-Py R² | R-Jl R² | Py-Jl R² | R NA | Py NA | Jl NA |
|--------|----------|---------|---------|----------|------|-------|-------|
| `b_pointcloud_senn_slp` | Recession | -49.7162 | -49.7162 | 1.0000 | 4718 | 1087 | 1087 |
| `b_pointcloud_linear_slp` | Recession | -41.9700 | -41.9700 | 1.0000 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_senn_slp` | Recession | -24.7931 | -24.7931 | 1.0000 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_linear_slp` | Recession | -14.3452 | -14.3452 | 1.0000 | 4718 | 1087 | 1087 |
| `b_pointcloud_mk_rho` | Recession | -6.3540 | -6.3540 | 1.0000 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_mk_rho` | Recession | -5.6652 | -5.6652 | 1.0000 | 4718 | 1087 | 1087 |
| `b_pointcloud_spearman_rho` | Recession | -4.0661 | -4.0661 | 1.0000 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_spearman_rho` | Recession | -3.5203 | -3.5203 | 1.0000 | 4718 | 1087 | 1087 |
| `b_pointcloud_mk_pval` | Recession | -1.6328 | -1.5892 | 0.9997 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_mk_pval` | Recession | -1.5950 | -1.5502 | 0.9997 | 4718 | 1087 | 1087 |
| `concavity_mk_rho` | Recession | -1.3424 | -1.3424 | 1.0000 | 4340 | 919 | 919 |
| `log_a_pointcloud_spearman_pval` | Recession | -1.3022 | -1.3022 | 1.0000 | 4718 | 1087 | 1087 |
| `b_pointcloud_spearman_pval` | Recession | -1.2769 | -1.2769 | 1.0000 | 4718 | 1087 | 1087 |
| `concavity_senn_slp` | Recession | -1.2660 | -1.2660 | 1.0000 | 4340 | 919 | 919 |
| `concavity_mean` | Recession | -1.2651 | -1.2651 | 1.0000 | 4340 | 919 | 919 |
| `concavity_spearman_rho` | Recession | -1.2255 | -1.2255 | 1.0000 | 4340 | 919 | 919 |
| `concavity_mk_pval` | Recession | -0.9702 | -0.9463 | 0.9997 | 4340 | 919 | 919 |
| `concavity_median` | Recession | -0.9683 | -0.9683 | 1.0000 | 4340 | 919 | 919 |
| `log_a_seasonality_minimum_first_half` | Recession | -0.9531 | -0.9441 | 0.9863 | 4341 | 920 | 920 |
| `b_events_linear_slp` | Recession | -0.9023 | -0.9023 | 1.0000 | 4340 | 919 | 919 |
| `concavity_spearman_pval` | Recession | -0.8702 | -0.8702 | 1.0000 | 4340 | 919 | 919 |
| `b_events_mk_pval` | Recession | -0.8092 | -0.7899 | 0.9997 | 4340 | 919 | 919 |
| `log_a_seasonality_minimum_last_half` | Recession | -0.7838 | -0.7756 | 0.9873 | 4340 | 923 | 923 |
| `b_events_spearman_pval` | Recession | -0.7836 | -0.7836 | 1.0000 | 4340 | 919 | 919 |
| `log_a_events_mk_pval` | Recession | -0.7640 | -0.7475 | 0.9997 | 4340 | 919 | 919 |
| `log_a_events_spearman_pval` | Recession | -0.7211 | -0.7211 | 1.0000 | 4340 | 919 | 919 |
| `b_events_senn_slp` | Recession | -0.6871 | -0.6871 | 1.0000 | 4340 | 919 | 919 |
| `log_a_seasonality_minimum_all` | Recession | -0.6697 | -0.6601 | 0.9866 | 4340 | 919 | 919 |
| `concavity_linear_slp` | Recession | -0.6259 | -0.6259 | 1.0000 | 4340 | 919 | 919 |
| `b_events_mk_rho` | Recession | -0.4479 | -0.4479 | 1.0000 | 4340 | 919 | 919 |
| `b_events_spearman_rho` | Recession | -0.4291 | -0.4291 | 1.0000 | 4340 | 919 | 919 |
| `log_a_events_mk_rho` | Recession | -0.3771 | -0.3771 | 1.0000 | 4340 | 919 | 919 |
| `log_a_events_spearman_rho` | Recession | -0.3189 | -0.3189 | 1.0000 | 4340 | 919 | 919 |
| `b_events_mean` | Recession | -0.2075 | -0.2075 | 1.0000 | 4340 | 919 | 919 |
| `b_events_median` | Recession | -0.1939 | -0.1939 | 1.0000 | 4340 | 919 | 919 |
| `log_a_seasonality_amplitude_last_half` | Recession | -0.1169 | -0.1169 | 1.0000 | 4340 | 923 | 923 |
| `log_a_seasonality_amplitude_all` | Recession | -0.1035 | -0.1034 | 1.0000 | 4340 | 919 | 919 |
| `log_a_events_senn_slp` | Recession | -0.0906 | -0.0906 | 1.0000 | 4340 | 919 | 919 |
| `log_a_seasonality_amplitude_first_half` | Recession | -0.0702 | -0.0701 | 1.0000 | 4341 | 920 | 920 |
| `log_a_events_linear_slp` | Recession | -0.0328 | -0.0328 | 1.0000 | 4340 | 919 | 919 |
| `b_pointcloud_median` | Recession | 0.4502 | 0.4502 | 1.0000 | 4718 | 1087 | 1087 |
| `b_pointcloud_mean` | Recession | 0.4599 | 0.4599 | 1.0000 | 4718 | 1087 | 1087 |
| `log_a_events_mean` | Recession | 0.7580 | 0.7580 | 1.0000 | 4340 | 919 | 919 |
| `log_a_events_median` | Recession | 0.8240 | 0.8240 | 1.0000 | 4340 | 919 | 919 |
| `log_a_pointcloud_median` | Recession | 0.8445 | 0.8445 | 1.0000 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_mean` | Recession | 0.8720 | 0.8720 | 1.0000 | 4718 | 1087 | 1087 |
| `FDC90th_mk_pval` | FDC | 0.9792 | 0.9812 | 0.9979 | 569 | 592 | 592 |
| `FDC90th_spearman_pval` | FDC | 0.9808 | 0.9819 | 0.9985 | 569 | 592 | 592 |
| `n_low_pulses_all_mk_rho` | Pulse Metrics | 0.9883 | 0.9882 | 0.9999 | 569 | 836 | 836 |

## NA Mismatch Analysis

Columns with >100 NA mismatches in any pair: **48**

| Column | NA R-Py | NA R-Jl | NA Py-Jl | R NA | Py NA | Jl NA |
|--------|---------|---------|----------|------|-------|-------|
| `b_pointcloud_spearman_rho` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_pointcloud_spearman_pval` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_pointcloud_senn_slp` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_pointcloud_mk_rho` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_pointcloud_mk_pval` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_pointcloud_median` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_pointcloud_mean` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_pointcloud_linear_slp` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_linear_slp` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_mean` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_median` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_mk_pval` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_mk_rho` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_senn_slp` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_spearman_pval` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `log_a_pointcloud_spearman_rho` | 3631 | 3631 | 0 | 4718 | 1087 | 1087 |
| `b_events_linear_slp` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `b_events_mean` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `b_events_median` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `b_events_mk_pval` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `b_events_mk_rho` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `b_events_senn_slp` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `b_events_spearman_pval` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `b_events_spearman_rho` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_linear_slp` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_mean` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_median` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_mk_pval` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_mk_rho` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_senn_slp` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_spearman_pval` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `concavity_spearman_rho` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_spearman_rho` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_spearman_pval` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_senn_slp` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_mk_rho` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_mk_pval` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_median` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_mean` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_events_linear_slp` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_seasonality_amplitude_all` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_seasonality_amplitude_first_half` | 3421 | 3421 | 0 | 4341 | 920 | 920 |
| `log_a_seasonality_minimum_all` | 3421 | 3421 | 0 | 4340 | 919 | 919 |
| `log_a_seasonality_minimum_first_half` | 3421 | 3421 | 0 | 4341 | 920 | 920 |
| `log_a_seasonality_amplitude_last_half` | 3417 | 3417 | 0 | 4340 | 923 | 923 |
| `log_a_seasonality_minimum_last_half` | 3417 | 3417 | 0 | 4340 | 923 | 923 |
| `n_low_pulses_all_mk_pval` | 277 | 277 | 0 | 569 | 836 | 836 |
| `n_low_pulses_all_mk_rho` | 277 | 277 | 0 | 569 | 836 | 836 |

## Column Coverage

| Scope | Count |
|-------|-------|
| All 3 languages | 543 |
| R + Python only | 0 |
| R + Julia only | 0 |
| Python + Julia only | 52 |
| Only R | 8 |
| Only Python | 0 |
| Only Julia | 0 |

### Columns only in R (8)

- `elasticity_linear_slp`
- `elasticity_mean`
- `elasticity_median`
- `elasticity_mk_pval`
- `elasticity_mk_rho`
- `elasticity_senn_slp`
- `elasticity_spearman_pval`
- `elasticity_spearman_rho`

## Deep Dive: Worst 10 Metrics

### `b_pointcloud_senn_slp`

**Identity R²**: R-Py=-49.7162, R-Jl=-49.7162, Py-Jl=1.0000
**Spearman rho**: R-Py=0.1338, R-Jl=0.1338, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.1721 | -0.6156 | -0.6156 |
| Median | 0.0003 | 0.0007 | 0.0007 |
| Max | 0.2992 | 0.7484 | 0.7484 |
| NAs | 4718 | 1087 | 1087 |

### `b_pointcloud_linear_slp`

**Identity R²**: R-Py=-41.9700, R-Jl=-41.9700, Py-Jl=1.0000
**Spearman rho**: R-Py=0.1645, R-Jl=0.1645, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.1634 | -0.6156 | -0.6156 |
| Median | 0.0002 | 0.0005 | 0.0005 |
| Max | 0.2806 | 0.4334 | 0.4334 |
| NAs | 4718 | 1087 | 1087 |

### `log_a_pointcloud_senn_slp`

**Identity R²**: R-Py=-24.7931, R-Jl=-24.7931, Py-Jl=1.0000
**Spearman rho**: R-Py=0.2208, R-Jl=0.2208, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.2877 | -0.8285 | -0.8285 |
| Median | 0.0001 | -0.0006 | -0.0006 |
| Max | 0.2818 | 1.9481 | 1.9481 |
| NAs | 4718 | 1087 | 1087 |

### `log_a_pointcloud_linear_slp`

**Identity R²**: R-Py=-14.3452, R-Jl=-14.3452, Py-Jl=1.0000
**Spearman rho**: R-Py=0.2541, R-Jl=0.2541, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.2877 | -0.6290 | -0.6290 |
| Median | 0.0004 | -0.0005 | -0.0005 |
| Max | 0.3363 | 0.6902 | 0.6902 |
| NAs | 4718 | 1087 | 1087 |

### `b_pointcloud_mk_rho`

**Identity R²**: R-Py=-6.3540, R-Jl=-6.3540, Py-Jl=1.0000
**Spearman rho**: R-Py=0.1088, R-Jl=0.1088, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -1.0000 | -1.0000 | -1.0000 |
| Median | 0.0000 | 0.0219 | 0.0219 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 1087 | 1087 |

### `log_a_pointcloud_mk_rho`

**Identity R²**: R-Py=-5.6652, R-Jl=-5.6652, Py-Jl=1.0000
**Spearman rho**: R-Py=0.2392, R-Jl=0.2392, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -1.0000 | -1.0000 | -1.0000 |
| Median | 0.0000 | -0.0172 | -0.0172 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 1087 | 1087 |

### `b_pointcloud_spearman_rho`

**Identity R²**: R-Py=-4.0661, R-Jl=-4.0661, Py-Jl=1.0000
**Spearman rho**: R-Py=0.1146, R-Jl=0.1146, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -1.0000 | -1.0000 | -1.0000 |
| Median | 0.0000 | 0.0344 | 0.0344 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 1087 | 1087 |

### `log_a_pointcloud_spearman_rho`

**Identity R²**: R-Py=-3.5203, R-Jl=-3.5203, Py-Jl=1.0000
**Spearman rho**: R-Py=0.2361, R-Jl=0.2361, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -1.0000 | -1.0000 | -1.0000 |
| Median | 0.0118 | -0.0269 | -0.0269 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 1087 | 1087 |

### `b_pointcloud_mk_pval`

**Identity R²**: R-Py=-1.6328, R-Jl=-1.5892, Py-Jl=0.9997
**Spearman rho**: R-Py=0.0739, R-Jl=0.0736, Py-Jl=0.9998

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0008 | 0.0000 | 0.0000 |
| Median | 0.7105 | 0.4320 | 0.4298 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 1087 | 1087 |

### `log_a_pointcloud_mk_pval`

**Identity R²**: R-Py=-1.5950, R-Jl=-1.5502, Py-Jl=0.9997
**Spearman rho**: R-Py=0.0834, R-Jl=0.0835, Py-Jl=0.9998

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0007 | 0.0000 | 0.0000 |
| Median | 0.7071 | 0.4167 | 0.4173 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 1087 | 1087 |

## Summary

| Agreement Level | Count | % |
|-----------------|-------|---|
| Perfect (R² >= 0.999) | 476 | 87.7% |
| Good (0.99 <= R² < 0.999) | 18 | 3.3% |
| Poor (R² < 0.99) | 49 | 9.0% |
| **Total compared** | **543** | **100%** |

---
*Generated by `docs/benchmarks/compare_three_way.py`*

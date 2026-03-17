# Three-Way Cross-Language Comparison Report

**Generated**: 2026-03-17 08:45:26

## Input Files

| File | Timestamp | Gages | Columns |
|------|-----------|-------|---------|
| R (`r_signatures.csv`) | 2026-03-17 07:47:44 | 5707 | 572 |
| Python (`python_signatures.csv`) | 2026-03-16 18:29:56 | 7369 | 583 |
| Julia (`julia_signatures.csv`) | 2026-03-16 17:30:01 | 7369 | 583 |

**Common gages (all 3)**: 5,707
**Common signature columns**: 551

## Performance Comparison

| Metric | Python | Julia | Ratio (Py/Jl) |
|--------|--------|-------|---------------|
| Total time | 4732s (78.9 min) | 1173s (19.6 min) | 4.0x |
| Gages processed | 7,369 | 7,369 | - |
| Processing rate | 1.56/s | 6.28/s | 4.0x faster |
| R benchmark | ~1-2 hours (estimated) | - | - |

### Phase Breakdown

| Phase | Python | Julia |
|-------|--------|-------|
| filter_gages | 218s | 184s |
| load_climate | 216s | 167s |
| load_streamflow | 45s | 17s |
| metadata_qaqc_save | 6s | 23s |
| process_signatures | 4248s | 782s |

## Overall Identity R² Summary

R² of the identity line (y = x): measures whether implementations produce identical values, not just correlated values.

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 | Cols < 0.95 | Cols < 0.90 |
|------|---------|-----------|--------|-------------|-------------|-------------|
| R vs Python | 0.996777 | 0.999995 | 0.6745 | 17 | 5 | 5 |
| R vs Julia | 0.996497 | 0.999994 | 0.6745 | 19 | 5 | 5 |
| Python vs Julia | 0.999666 | 1.000000 | 0.9771 | 3 | 0 | 0 |

### Spearman Rank Correlation (secondary diagnostic)

| Pair | Mean rho | Min rho | Cols < 0.99 |
|------|----------|---------|-------------|
| R vs Python | 0.998843 | 0.8498 | 4 |
| R vs Julia | 0.998777 | 0.8355 | 4 |
| Python vs Julia | 0.999943 | 0.9976 | 0 |

## Agreement by Signature Category

| Category | Total Cols | Perfect (>=0.999) | Good (>=0.99) | Poor (<0.99) | Min R² |
|----------|-----------|-------------------|---------------|-------------|--------|
| Baseflow | 16 | 8 | 6 | 2 | 0.9778 |
| Elasticity | 9 | 2 | 5 | 2 | 0.9841 |
| FDC | 24 | 14 | 3 | 7 | 0.6745 |
| Flashiness | 8 | 5 | 2 | 1 | 0.9856 |
| Flow Percentiles | 128 | 120 | 8 | 0 | 0.9918 |
| Flow Timing | 104 | 104 | 0 | 0 | 0.9994 |
| Flow Volumes | 40 | 40 | 0 | 0 | 0.9997 |
| Pulse Metrics | 112 | 84 | 26 | 2 | 0.9768 |
| Q-P Seasonality | 16 | 11 | 5 | 0 | 0.9974 |
| Recession | 46 | 35 | 7 | 4 | 0.6846 |
| Runoff Ratios | 40 | 40 | 0 | 0 | 0.9996 |
| Storage | 8 | 4 | 2 | 2 | 0.9851 |

## Columns with min Identity R² < 0.99 (20 columns)

| Column | Category | R-Py R² | R-Jl R² | Py-Jl R² | R NA | Py NA | Jl NA |
|--------|----------|---------|---------|----------|------|-------|-------|
| `FDC90th_senn_slp` | FDC | 0.6745 | 0.6745 | 0.9998 | 1 | 0 | 0 |
| `log_a_pointcloud_mk_pval` | Recession | 0.7114 | 0.6846 | 0.9983 | 4732 | 4732 | 4732 |
| `b_pointcloud_mk_pval` | Recession | 0.7244 | 0.6996 | 0.9983 | 4732 | 4732 | 4732 |
| `b_pointcloud_spearman_pval` | Recession | 0.7587 | 0.7587 | 1.0000 | 4732 | 4732 | 4732 |
| `log_a_pointcloud_spearman_pval` | Recession | 0.7709 | 0.7709 | 1.0000 | 4732 | 4732 | 4732 |
| `FDC90th_linear_slp` | FDC | 0.9597 | 0.9597 | 1.0000 | 1 | 0 | 0 |
| `FDC90th_mk_pval` | FDC | 0.9762 | 0.9782 | 0.9968 | 1 | 36 | 36 |
| `n_low_pulses_year_mk_rho` | Pulse Metrics | 0.9998 | 0.9768 | 0.9771 | 0 | 34 | 34 |
| `BFI_LyneHollick_mk_pval` | Baseflow | 0.9778 | 0.9783 | 0.9997 | 0 | 0 | 0 |
| `FDC90th_spearman_pval` | FDC | 0.9781 | 0.9792 | 0.9975 | 1 | 36 | 36 |
| `BFI_LyneHollick_spearman_pval` | Baseflow | 0.9784 | 0.9784 | 1.0000 | 0 | 0 | 0 |
| `n_low_pulses_all_mk_rho` | Pulse Metrics | 0.9999 | 0.9793 | 0.9795 | 0 | 309 | 309 |
| `FDC90th_median` | FDC | 0.9832 | 0.9832 | 0.9999 | 1 | 0 | 0 |
| `elasticity_mk_pval` | Elasticity | 0.9841 | 0.9841 | 1.0000 | 4 | 4 | 4 |
| `elasticity_spearman_pval` | Elasticity | 0.9843 | 0.9843 | 1.0000 | 4 | 4 | 4 |
| `avg_storage_mk_pval` | Storage | 0.9851 | 0.9858 | 0.9996 | 4 | 4 | 4 |
| `FDC90th_mean` | FDC | 0.9852 | 0.9852 | 1.0000 | 1 | 0 | 0 |
| `flashinessRB_linear_slp` | Flashiness | 0.9998 | 0.9856 | 0.9859 | 0 | 0 | 0 |
| `avg_storage_spearman_pval` | Storage | 0.9861 | 0.9861 | 1.0000 | 4 | 4 | 4 |
| `FDC90th_mk_rho` | FDC | 0.9900 | 0.9935 | 0.9912 | 1 | 36 | 36 |

## NA Mismatch Analysis

Columns with >100 NA mismatches in any pair: **2**

| Column | NA R-Py | NA R-Jl | NA Py-Jl | R NA | Py NA | Jl NA |
|--------|---------|---------|----------|------|-------|-------|
| `n_low_pulses_all_mk_pval` | 309 | 309 | 0 | 0 | 309 | 309 |
| `n_low_pulses_all_mk_rho` | 309 | 309 | 0 | 0 | 309 | 309 |

## Column Coverage

| Scope | Count |
|-------|-------|
| All 3 languages | 551 |
| R + Python only | 0 |
| R + Julia only | 0 |
| Python + Julia only | 0 |
| Only R | 0 |
| Only Python | 0 |
| Only Julia | 0 |

## Deep Dive: Worst 10 Metrics

### `FDC90th_senn_slp`

**Identity R²**: R-Py=0.6745, R-Jl=0.6745, Py-Jl=0.9998
**Spearman rho**: R-Py=0.9990, R-Jl=0.9990, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -2.3065 | -1.4014 | -1.4014 |
| Median | 0.0000 | 0.0000 | 0.0000 |
| Max | 1.0044 | 1.0044 | 1.0044 |
| NAs | 1 | 0 | 0 |

### `log_a_pointcloud_mk_pval`

**Identity R²**: R-Py=0.7114, R-Jl=0.6846, Py-Jl=0.9983
**Spearman rho**: R-Py=0.8498, R-Jl=0.8355, Py-Jl=0.9990

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0007 | 0.0005 | 0.0007 |
| Median | 0.7071 | 0.5484 | 0.5376 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `b_pointcloud_mk_pval`

**Identity R²**: R-Py=0.7244, R-Jl=0.6996, Py-Jl=0.9983
**Spearman rho**: R-Py=0.8555, R-Jl=0.8428, Py-Jl=0.9990

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0008 | 0.0002 | 0.0008 |
| Median | 0.7105 | 0.6259 | 0.6204 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `b_pointcloud_spearman_pval`

**Identity R²**: R-Py=0.7587, R-Jl=0.7587, Py-Jl=1.0000
**Spearman rho**: R-Py=0.9040, R-Jl=0.9043, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0001 | 0.0000 | 0.0000 |
| Median | 0.5639 | 0.5457 | 0.5457 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `log_a_pointcloud_spearman_pval`

**Identity R²**: R-Py=0.7709, R-Jl=0.7709, Py-Jl=1.0000
**Spearman rho**: R-Py=0.9071, R-Jl=0.9072, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0001 | 0.0000 | 0.0000 |
| Median | 0.5032 | 0.4986 | 0.4986 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `FDC90th_linear_slp`

**Identity R²**: R-Py=0.9597, R-Jl=0.9597, Py-Jl=1.0000
**Spearman rho**: R-Py=0.9963, R-Jl=0.9963, Py-Jl=0.9997

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -2.5128 | -2.5128 | -2.5128 |
| Median | 0.0000 | 0.0000 | 0.0000 |
| Max | 2.5873 | 2.5873 | 2.5873 |
| NAs | 1 | 0 | 0 |

### `FDC90th_mk_pval`

**Identity R²**: R-Py=0.9762, R-Jl=0.9782, Py-Jl=0.9968
**Spearman rho**: R-Py=0.9907, R-Jl=0.9912, Py-Jl=0.9987

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0000 | 0.0000 | 0.0000 |
| Median | 0.3065 | 0.2926 | 0.2962 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 1 | 36 | 36 |

### `n_low_pulses_year_mk_rho`

**Identity R²**: R-Py=0.9998, R-Jl=0.9768, Py-Jl=0.9771
**Spearman rho**: R-Py=0.9998, R-Jl=0.9974, Py-Jl=0.9976

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.6562 | -0.6562 | -0.6367 |
| Median | 0.0032 | 0.0010 | 0.0010 |
| Max | 1.0000 | 0.7340 | 0.5877 |
| NAs | 0 | 34 | 34 |

### `BFI_LyneHollick_mk_pval`

**Identity R²**: R-Py=0.9778, R-Jl=0.9783, Py-Jl=0.9997
**Spearman rho**: R-Py=0.9929, R-Jl=0.9929, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0000 | 0.0000 | 0.0000 |
| Median | 0.2804 | 0.2792 | 0.2804 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 0 | 0 | 0 |

### `FDC90th_spearman_pval`

**Identity R²**: R-Py=0.9781, R-Jl=0.9792, Py-Jl=0.9975
**Spearman rho**: R-Py=0.9911, R-Jl=0.9916, Py-Jl=0.9989

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0000 | 0.0000 | 0.0000 |
| Median | 0.2998 | 0.2898 | 0.2892 |
| Max | 1.0000 | 0.9990 | 0.9990 |
| NAs | 1 | 36 | 36 |

## Summary

| Agreement Level | Count | % |
|-----------------|-------|---|
| Perfect (R² >= 0.999) | 467 | 84.8% |
| Good (0.99 <= R² < 0.999) | 64 | 11.6% |
| Poor (R² < 0.99) | 20 | 3.6% |
| **Total compared** | **551** | **100%** |

---
*Generated by `docs/benchmarks/compare_three_way.py`*

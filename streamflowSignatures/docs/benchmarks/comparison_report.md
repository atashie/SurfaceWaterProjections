# Three-Way Cross-Language Comparison Report

**Generated**: 2026-04-11 17:33:23

## Input Files

| File | Timestamp | Gages | Columns |
|------|-----------|-------|---------|
| R (`r_signatures.csv`) | 2026-04-11 17:33:03 | 5687 | 572 |
| Python (`python_signatures.csv`) | 2026-04-10 07:59:28 | 7297 | 583 |
| Julia (`julia_signatures.csv`) | 2026-04-07 14:57:37 | 7297 | 583 |

**Common gages (all 3)**: 5,687
**Common signature columns**: 551

## Performance Comparison

| Metric | Python | Julia | Ratio (Py/Jl) |
|--------|--------|-------|---------------|
| Total time | 48494s (808.2 min) | 1308s (21.8 min) | 37.1x |
| Gages processed | 7,297 | 7,297 | - |
| Processing rate | 0.15/s | 5.58/s | 37.1x faster |
| R benchmark | ~1-2 hours (estimated) | - | - |

### Phase Breakdown

| Phase | Python | Julia |
|-------|--------|-------|
| filter_gages | 1550s | 374s |
| load_climate | 82s | 24s |
| load_streamflow | 19s | 18s |
| metadata_qaqc_save | 4s | 24s |
| process_signatures | 46840s | 868s |

## Overall Identity R² Summary

R² of the identity line (y = x): measures whether implementations produce identical values, not just correlated values.

| Pair | Mean R² | Median R² | Min R² | Cols < 0.99 | Cols < 0.95 | Cols < 0.90 |
|------|---------|-----------|--------|-------------|-------------|-------------|
| R vs Python | 0.997992 | 1.000000 | 0.7152 | 6 | 4 | 4 |
| R vs Julia | 0.997742 | 1.000000 | 0.6886 | 8 | 4 | 4 |
| Python vs Julia | 0.999715 | 1.000000 | 0.9797 | 2 | 0 | 0 |

### Spearman Rank Correlation (secondary diagnostic)

| Pair | Mean rho | Min rho | Cols < 0.99 |
|------|----------|---------|-------------|
| R vs Python | 0.999077 | 0.8529 | 4 |
| R vs Julia | 0.999011 | 0.8374 | 4 |
| Python vs Julia | 0.999948 | 0.9975 | 0 |

## Agreement by Signature Category

| Category | Total Cols | Perfect (>=0.999) | Good (>=0.99) | Poor (<0.99) | Min R² |
|----------|-----------|-------------------|---------------|-------------|--------|
| Baseflow | 16 | 16 | 0 | 0 | 0.9997 |
| Elasticity | 9 | 9 | 0 | 0 | 1.0000 |
| FDC | 24 | 20 | 2 | 2 | 0.9796 |
| Flashiness | 8 | 8 | 0 | 0 | 0.9998 |
| Flow Percentiles | 128 | 120 | 8 | 0 | 0.9916 |
| Flow Timing | 104 | 104 | 0 | 0 | 0.9994 |
| Flow Volumes | 40 | 40 | 0 | 0 | 0.9996 |
| Pulse Metrics | 112 | 101 | 9 | 2 | 0.9797 |
| Q-P Seasonality | 16 | 16 | 0 | 0 | 0.9996 |
| Recession | 46 | 42 | 0 | 4 | 0.6886 |
| Runoff Ratios | 40 | 40 | 0 | 0 | 0.9994 |
| Storage | 8 | 6 | 2 | 0 | 0.9977 |

## Columns with min Identity R² < 0.99 (8 columns)

| Column | Category | R-Py R² | R-Jl R² | Py-Jl R² | R NA | Py NA | Jl NA |
|--------|----------|---------|---------|----------|------|-------|-------|
| `log_a_pointcloud_mk_pval` | Recession | 0.7152 | 0.6886 | 0.9983 | 4718 | 4718 | 4718 |
| `b_pointcloud_mk_pval` | Recession | 0.7246 | 0.6994 | 0.9983 | 4718 | 4718 | 4718 |
| `b_pointcloud_spearman_pval` | Recession | 0.7621 | 0.7621 | 1.0000 | 4718 | 4718 | 4718 |
| `log_a_pointcloud_spearman_pval` | Recession | 0.7738 | 0.7738 | 1.0000 | 4718 | 4718 | 4718 |
| `FDC90th_mk_pval` | FDC | 0.9796 | 0.9816 | 0.9979 | 569 | 597 | 597 |
| `n_low_pulses_year_mk_rho` | Pulse Metrics | 1.0000 | 0.9797 | 0.9797 | 569 | 595 | 595 |
| `n_low_pulses_all_mk_rho` | Pulse Metrics | 1.0000 | 0.9802 | 0.9802 | 569 | 844 | 844 |
| `FDC90th_spearman_pval` | FDC | 0.9812 | 0.9823 | 0.9985 | 569 | 597 | 597 |

## NA Mismatch Analysis

Columns with >100 NA mismatches in any pair: **2**

| Column | NA R-Py | NA R-Jl | NA Py-Jl | R NA | Py NA | Jl NA |
|--------|---------|---------|----------|------|-------|-------|
| `n_low_pulses_all_mk_pval` | 275 | 275 | 0 | 569 | 844 | 844 |
| `n_low_pulses_all_mk_rho` | 275 | 275 | 0 | 569 | 844 | 844 |

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

## Deep Dive: Worst 8 Metrics

### `log_a_pointcloud_mk_pval`

**Identity R²**: R-Py=0.7152, R-Jl=0.6886, Py-Jl=0.9983
**Spearman rho**: R-Py=0.8529, R-Jl=0.8374, Py-Jl=0.9988

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0007 | 0.0005 | 0.0007 |
| Median | 0.7071 | 0.5484 | 0.5362 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 4718 | 4718 |

### `b_pointcloud_mk_pval`

**Identity R²**: R-Py=0.7246, R-Jl=0.6994, Py-Jl=0.9983
**Spearman rho**: R-Py=0.8554, R-Jl=0.8434, Py-Jl=0.9992

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0008 | 0.0002 | 0.0008 |
| Median | 0.7105 | 0.6122 | 0.6022 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 4718 | 4718 |

### `b_pointcloud_spearman_pval`

**Identity R²**: R-Py=0.7621, R-Jl=0.7621, Py-Jl=1.0000
**Spearman rho**: R-Py=0.9073, R-Jl=0.9077, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0001 | 0.0000 | 0.0000 |
| Median | 0.5639 | 0.5441 | 0.5441 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 4718 | 4718 |

### `log_a_pointcloud_spearman_pval`

**Identity R²**: R-Py=0.7738, R-Jl=0.7738, Py-Jl=1.0000
**Spearman rho**: R-Py=0.9095, R-Jl=0.9096, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0001 | 0.0000 | 0.0000 |
| Median | 0.4991 | 0.4986 | 0.4986 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4718 | 4718 | 4718 |

### `FDC90th_mk_pval`

**Identity R²**: R-Py=0.9796, R-Jl=0.9816, Py-Jl=0.9979
**Spearman rho**: R-Py=0.9924, R-Jl=0.9929, Py-Jl=0.9990

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0000 | 0.0000 | 0.0000 |
| Median | 0.3078 | 0.2929 | 0.2976 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 569 | 597 | 597 |

### `n_low_pulses_year_mk_rho`

**Identity R²**: R-Py=1.0000, R-Jl=0.9797, Py-Jl=0.9797
**Spearman rho**: R-Py=1.0000, R-Jl=0.9977, Py-Jl=0.9977

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.6614 | -0.6614 | -0.6367 |
| Median | -0.0011 | -0.0021 | -0.0019 |
| Max | 1.0000 | 0.7183 | 0.5605 |
| NAs | 569 | 595 | 595 |

### `n_low_pulses_all_mk_rho`

**Identity R²**: R-Py=1.0000, R-Jl=0.9802, Py-Jl=0.9802
**Spearman rho**: R-Py=1.0000, R-Jl=0.9992, Py-Jl=0.9992

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.7581 | -0.7581 | -0.6831 |
| Median | -0.0207 | -0.0343 | -0.0319 |
| Max | 1.0000 | 0.7391 | 0.6302 |
| NAs | 569 | 844 | 844 |

### `FDC90th_spearman_pval`

**Identity R²**: R-Py=0.9812, R-Jl=0.9823, Py-Jl=0.9985
**Spearman rho**: R-Py=0.9927, R-Jl=0.9931, Py-Jl=0.9992

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0000 | 0.0000 | 0.0000 |
| Median | 0.3000 | 0.2901 | 0.2892 |
| Max | 1.0000 | 0.9990 | 0.9990 |
| NAs | 569 | 597 | 597 |

## Summary

| Agreement Level | Count | % |
|-----------------|-------|---|
| Perfect (R² >= 0.999) | 522 | 94.7% |
| Good (0.99 <= R² < 0.999) | 21 | 3.8% |
| Poor (R² < 0.99) | 8 | 1.5% |
| **Total compared** | **551** | **100%** |

---
*Generated by `docs/benchmarks/compare_three_way.py`*

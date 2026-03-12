# Three-Way Cross-Language Comparison Report

**Generated**: 2026-03-12 14:08:06

## Input Files

| File | Timestamp | Gages | Columns |
|------|-----------|-------|---------|
| R (`r_signatures.csv`) | 2026-03-10 08:05:24 | 5707 | 572 |
| Python (`python_signatures.csv`) | 2026-03-12 13:22:10 | 7369 | 571 |
| Julia (`julia_signatures.csv`) | 2026-03-12 11:22:07 | 7369 | 571 |

**Common gages (all 3)**: 5,707
**Common signature columns**: 551

## Performance Comparison

| Metric | Python | Julia | Ratio (Py/Jl) |
|--------|--------|-------|---------------|
| Total time | 6321s (105.3 min) | 855s (14.2 min) | 7.4x |
| Gages processed | 7,369 | 7,369 | - |
| Processing rate | 1.17/s | 8.62/s | 7.4x faster |
| R benchmark | ~1-2 hours (estimated) | - | - |

### Phase Breakdown

| Phase | Python | Julia |
|-------|--------|-------|
| filter_gages | 81s | 90s |
| load_climate | 46s | 20s |
| load_streamflow | 12s | 14s |
| metadata_qaqc_save | 3s | 12s |
| process_signatures | 6179s | 719s |

## Overall Spearman Correlation Summary

| Pair | Mean rho | Median rho | Min rho | Cols < 0.99 | Cols < 0.95 | Cols < 0.90 |
|------|----------|------------|---------|-------------|-------------|-------------|
| R vs Python | 0.998558 | 0.999998 | 0.8498 | 7 | 5 | 2 |
| R vs Julia | 0.998782 | 0.999997 | 0.8355 | 4 | 4 | 2 |
| Python vs Julia | 0.999631 | 1.000000 | 0.9137 | 3 | 1 | 0 |

## Agreement by Signature Category

| Category | Total Cols | Perfect (>=0.999) | Good (>=0.99) | Poor (<0.99) | Min rho |
|----------|-----------|-------------------|---------------|-------------|--------|
| Baseflow | 16 | 3 | 10 | 3 | 0.9136 |
| Elasticity | 9 | 3 | 6 | 0 | 0.9956 |
| FDC | 24 | 18 | 6 | 0 | 0.9912 |
| Flashiness | 8 | 8 | 0 | 0 | 0.9992 |
| Flow Percentiles | 128 | 123 | 5 | 0 | 0.9976 |
| Flow Timing | 104 | 104 | 0 | 0 | 1.0000 |
| Flow Volumes | 40 | 40 | 0 | 0 | 0.9999 |
| Pulse Metrics | 112 | 105 | 7 | 0 | 0.9974 |
| Q-P Seasonality | 16 | 14 | 2 | 0 | 0.9989 |
| Recession | 46 | 38 | 4 | 4 | 0.8355 |
| Runoff Ratios | 40 | 40 | 0 | 0 | 0.9998 |
| Storage | 8 | 2 | 6 | 0 | 0.9929 |

## Columns with min Spearman < 0.99 (7 columns)

| Column | Category | R-Py | R-Jl | Py-Jl | R NA | Py NA | Jl NA |
|--------|----------|------|------|-------|------|-------|-------|
| `log_a_pointcloud_mk_pval` | Recession | 0.8498 | 0.8355 | 0.9990 | 4732 | 4732 | 4732 |
| `b_pointcloud_mk_pval` | Recession | 0.8555 | 0.8428 | 0.9990 | 4732 | 4732 | 4732 |
| `b_pointcloud_spearman_pval` | Recession | 0.9040 | 0.9043 | 1.0000 | 4732 | 4732 | 4732 |
| `log_a_pointcloud_spearman_pval` | Recession | 0.9071 | 0.9072 | 1.0000 | 4732 | 4732 | 4732 |
| `BFI_Eckhardt_linear_slp` | Baseflow | 0.9136 | 1.0000 | 0.9137 | 0 | 0 | 0 |
| `BFI_Eckhardt_mk_pval` | Baseflow | 0.9781 | 0.9996 | 0.9782 | 0 | 0 | 0 |
| `BFI_Eckhardt_spearman_pval` | Baseflow | 0.9785 | 0.9996 | 0.9785 | 0 | 0 | 0 |

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

## Deep Dive: Worst 7 Metrics

### `log_a_pointcloud_mk_pval`

**Spearman**: R-Py=0.8498, R-Jl=0.8355, Py-Jl=0.9990

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0007 | 0.0005 | 0.0007 |
| Median | 0.7071 | 0.5484 | 0.5376 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `b_pointcloud_mk_pval`

**Spearman**: R-Py=0.8555, R-Jl=0.8428, Py-Jl=0.9990

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0008 | 0.0002 | 0.0008 |
| Median | 0.7105 | 0.6259 | 0.6204 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `b_pointcloud_spearman_pval`

**Spearman**: R-Py=0.9040, R-Jl=0.9043, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0001 | 0.0000 | 0.0000 |
| Median | 0.5639 | 0.5457 | 0.5457 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `log_a_pointcloud_spearman_pval`

**Spearman**: R-Py=0.9071, R-Jl=0.9072, Py-Jl=1.0000

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0001 | 0.0000 | 0.0000 |
| Median | 0.5032 | 0.4986 | 0.4986 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 4732 | 4732 | 4732 |

### `BFI_Eckhardt_linear_slp`

**Spearman**: R-Py=0.9136, R-Jl=1.0000, Py-Jl=0.9137

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | -0.0125 | -0.0125 | -0.0125 |
| Median | -0.0000 | -0.0000 | -0.0000 |
| Max | 0.0135 | 0.0135 | 0.0135 |
| NAs | 0 | 0 | 0 |

### `BFI_Eckhardt_mk_pval`

**Spearman**: R-Py=0.9781, R-Jl=0.9996, Py-Jl=0.9782

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0000 | 0.0000 | 0.0000 |
| Median | 0.3044 | 0.2933 | 0.3065 |
| Max | 1.0000 | 1.0000 | 1.0000 |
| NAs | 0 | 0 | 0 |

### `BFI_Eckhardt_spearman_pval`

**Spearman**: R-Py=0.9785, R-Jl=0.9996, Py-Jl=0.9785

| Stat | R | Python | Julia |
|------|---|--------|-------|
| Min | 0.0000 | 0.0000 | 0.0000 |
| Median | 0.2988 | 0.2970 | 0.2999 |
| Max | 1.0000 | 0.9997 | 0.9997 |
| NAs | 0 | 0 | 0 |

## Summary

| Agreement Level | Count | % |
|-----------------|-------|---|
| Perfect (rho >= 0.999) | 498 | 90.4% |
| Good (0.99 <= rho < 0.999) | 46 | 8.3% |
| Poor (rho < 0.99) | 7 | 1.3% |
| **Total compared** | **551** | **100%** |

---
*Generated by `benchmarks/compare_three_way.py`*

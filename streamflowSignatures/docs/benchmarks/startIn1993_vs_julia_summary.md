# Experiment 'startIn1993' vs Julia Baseline: Comparison Report

**Generated**: 2026-04-16 10:48:54

## Experiment Description

Water years >= 1993 only (restricts analysis period from ~1980).

## Input Files

| Dataset | File | Gages | Columns |
|---------|------|-------|---------|
| Julia Baseline | `julia_signatures.csv` | 7,313 | 627 |
| Experiment (startIn1993) | `startIn1993_signatures.csv` | 6,678 | 631 |

**Common gages**: 6,678
**Dropped gages** (in baseline, not in experiment): 635
**Added gages** (in experiment, not in baseline): 0

### Gage Differences

```

  Dropped (baseline only): 635 gages
    USGS: 485
    Canada: 150
    non-reference: 503
    reference: 132
    Examples: 01021200, 01060000, 01064000, 01064118, 01064300
```

### Years Per Gage (Common Gages)

```
Years per gage (common gages):
    Baseline:   mean=39.4, median=44.0
    Experiment: mean=30.1, median=32.0
    Mean diff:  -9.3 years
```

## High-Level Alignment Summary

### Distribution Statistics

| Metric | Value |
|--------|-------|
| Columns compared | 594 |
| Mean R2 (identity) | 0.433786 |
| Median R2 | 0.614506 |
| SD of R2 | 0.882145 |
| Min R2 | -17.3977 |

### Agreement Tiers

| Tier | Threshold | Count | % |
|------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 35 | 5.9% |
| Good | 0.99 <= R2 < 0.999 | 41 | 6.9% |
| Poor | 0.95 <= R2 < 0.99 | 58 | 9.8% |
| Low | 0.90 <= R2 < 0.95 | 9 | 1.5% |
| Very Low | 0.50 <= R2 < 0.90 | 228 | 38.4% |
| Extremely Low | R2 < 0.50 | 223 | 37.5% |
| **Total** | | **594** | **100%** |

## Agreement by Signature Category

| Category | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| Baseflow | 16 | 0 | 4 | 0 | 0 | 8 | 4 | 0.568517 | -0.0416 |
| Elasticity | 19 | 1 | 0 | 0 | 1 | 5 | 12 | -0.609689 | -17.3977 |
| FDC | 24 | 0 | 0 | 4 | 1 | 13 | 6 | 0.570051 | -0.0806 |
| Flashiness | 8 | 0 | 2 | 0 | 0 | 4 | 2 | 0.584195 | -0.0207 |
| Flow Percentiles | 128 | 21 | 11 | 0 | 0 | 32 | 64 | 0.324141 | -1.2369 |
| Flow Timing | 120 | 0 | 9 | 20 | 1 | 50 | 40 | 0.521255 | -0.2098 |
| Flow Volumes | 40 | 7 | 3 | 0 | 0 | 10 | 20 | 0.336564 | -1.5983 |
| Negative Flow | 8 | 1 | 1 | 0 | 0 | 4 | 2 | 0.624749 | 0.0345 |
| Pulse Metrics | 112 | 0 | 5 | 17 | 4 | 56 | 30 | 0.538062 | -1.4812 |
| Q-P Seasonality | 16 | 0 | 0 | 4 | 0 | 5 | 7 | 0.427636 | -0.2852 |
| Recession | 54 | 0 | 1 | 11 | 2 | 28 | 12 | 0.612385 | -0.0149 |
| Runoff Ratios | 41 | 5 | 4 | 1 | 0 | 12 | 19 | 0.381746 | -2.2003 |
| Storage | 8 | 0 | 1 | 1 | 0 | 1 | 5 | 0.434434 | -0.2649 |

## Agreement by Statistic Type

| Stat Type | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|-----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| mean | 73 | 20 | 25 | 24 | 3 | 1 | 0 | 0.985116 | 0.7739 |
| median | 73 | 14 | 16 | 34 | 4 | 5 | 0 | 0.976154 | 0.8269 |
| senn_slp | 73 | 0 | 0 | 0 | 0 | 45 | 28 | 0.326986 | -1.4812 |
| linear_slp | 73 | 0 | 0 | 0 | 0 | 46 | 27 | 0.295428 | -2.2003 |
| spearman_rho | 73 | 0 | 0 | 0 | 0 | 62 | 11 | 0.582436 | 0.2022 |
| spearman_pval | 73 | 0 | 0 | 0 | 0 | 0 | 73 | -0.049668 | -0.7360 |
| mk_rho | 73 | 0 | 0 | 0 | 0 | 63 | 10 | 0.595499 | 0.2354 |
| mk_pval | 73 | 0 | 0 | 0 | 0 | 0 | 73 | -0.044401 | -0.7466 |
| scalar | 10 | 1 | 0 | 0 | 2 | 6 | 1 | -1.006200 | -17.3977 |

## Columns with R2 < 0.99 (518 columns)

| Column | Category | Stat | R2 | Spearman | MAD | Max Diff | Baseline NAs | Experiment NAs | NA Mismatch |
|--------|----------|------|----|----------|-----|----------|--------------|----------------|-------------|
| `elasticity_years_total` | Elasticity | scalar | -17.3977 | 0.7523 | 10.3153 | 12.0000 | 1165 | 1165 | 0 |
| `spring_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | -2.2003 | 0.7479 | 0.0047 | 4.0210 | 1612 | 1566 | 420 |
| `Qfal_linear_slp` | Flow Volumes | linear_slp | -1.5983 | 0.6947 | 1.4328 | 3175.9690 | 645 | 544 | 513 |
| `dur_low_pulses_all_senn_slp` | Pulse Metrics | senn_slp | -1.4812 | 0.7648 | 0.0698 | 15.4390 | 4364 | 4141 | 755 |
| `Q10_senn_slp` | Flow Percentiles | senn_slp | -1.2369 | 0.7803 | 0.0116 | 32.7018 | 645 | 544 | 513 |
| `Q5_senn_slp` | Flow Percentiles | senn_slp | -1.2301 | 0.7788 | 0.0114 | 28.2927 | 645 | 544 | 513 |
| `Q1_senn_slp` | Flow Percentiles | senn_slp | -1.2166 | 0.7812 | 0.0103 | 27.0480 | 645 | 544 | 513 |
| `Q10_linear_slp` | Flow Percentiles | linear_slp | -1.1988 | 0.7940 | 0.0102 | 25.4106 | 645 | 544 | 513 |
| `Q20_linear_slp` | Flow Percentiles | linear_slp | -1.1582 | 0.8069 | 0.0106 | 27.6171 | 645 | 544 | 513 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | linear_slp | -1.0606 | 0.7613 | 0.0913 | 16.7031 | 4364 | 4141 | 755 |
| `Q25_senn_slp` | Flow Percentiles | senn_slp | -1.0546 | 0.7936 | 0.0139 | 43.0240 | 645 | 544 | 513 |
| `Q20_senn_slp` | Flow Percentiles | senn_slp | -1.0091 | 0.7975 | 0.0134 | 42.1578 | 645 | 544 | 513 |
| `Q30_senn_slp` | Flow Percentiles | senn_slp | -0.8716 | 0.7771 | 0.0150 | 45.2000 | 645 | 544 | 513 |
| `Q25_linear_slp` | Flow Percentiles | linear_slp | -0.8691 | 0.8102 | 0.0110 | 29.3314 | 645 | 544 | 513 |
| `Q30_linear_slp` | Flow Percentiles | linear_slp | -0.7940 | 0.8057 | 0.0118 | 31.5025 | 645 | 544 | 513 |
| `Q5_linear_slp` | Flow Percentiles | linear_slp | -0.7795 | 0.7948 | 0.0099 | 24.7992 | 645 | 544 | 513 |
| `elasticity_rolling_mk_pval` | Elasticity | mk_pval | -0.7466 | 0.1866 | 0.2782 | 1.0000 | 1165 | 1165 | 0 |
| `elasticity_rolling_spearman_pval` | Elasticity | spearman_pval | -0.7360 | 0.1961 | 0.2648 | 0.9967 | 1165 | 1165 | 0 |
| `Q1_linear_slp` | Flow Percentiles | linear_slp | -0.7139 | 0.7872 | 0.0086 | 19.2563 | 645 | 544 | 513 |
| `Qfal_senn_slp` | Flow Volumes | senn_slp | -0.6873 | 0.6813 | 1.5104 | 4084.4444 | 645 | 544 | 513 |
| `Q50_linear_slp` | Flow Percentiles | linear_slp | -0.5057 | 0.8181 | 0.0135 | 35.8286 | 645 | 544 | 513 |
| `Q40_linear_slp` | Flow Percentiles | linear_slp | -0.5009 | 0.8107 | 0.0131 | 34.0502 | 645 | 544 | 513 |
| `Q40_senn_slp` | Flow Percentiles | senn_slp | -0.4442 | 0.7855 | 0.0163 | 48.1071 | 645 | 544 | 513 |
| `spring_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | -0.4340 | 0.7318 | 0.0035 | 2.9330 | 1612 | 1566 | 420 |
| `elasticity_annual_spearman_pval` | Elasticity | spearman_pval | -0.3733 | 0.3169 | 0.2658 | 0.9896 | 1165 | 1165 | 0 |
| `elasticity_annual_mk_pval` | Elasticity | mk_pval | -0.3513 | 0.3256 | 0.2666 | 0.9650 | 1165 | 1165 | 0 |
| `summer_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.3258 | 0.3995 | 0.2599 | 0.9958 | 1786 | 1762 | 420 |
| `Qann_linear_slp` | Flow Volumes | linear_slp | -0.3230 | 0.7958 | 5.1932 | 11771.1424 | 645 | 544 | 513 |
| `summer_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.3088 | 0.4038 | 0.2627 | 0.9988 | 1786 | 1762 | 420 |
| `qp_bimodality_spearman_pval` | Q-P Seasonality | spearman_pval | -0.2852 | 0.3788 | 0.2538 | 0.9854 | 1615 | 1576 | 423 |
| `Q50_senn_slp` | Flow Percentiles | senn_slp | -0.2773 | 0.7956 | 0.0153 | 46.9455 | 645 | 544 | 513 |
| `Q60_linear_slp` | Flow Percentiles | linear_slp | -0.2757 | 0.8046 | 0.0139 | 33.0408 | 645 | 544 | 513 |
| `qp_bimodality_mk_pval` | Q-P Seasonality | mk_pval | -0.2651 | 0.3867 | 0.2547 | 0.9893 | 1615 | 1576 | 423 |
| `avg_storage_spearman_pval` | Storage | spearman_pval | -0.2649 | 0.3883 | 0.2544 | 0.9481 | 1627 | 1583 | 422 |
| `qp_slope_sd_spearman_pval` | Q-P Seasonality | spearman_pval | -0.2554 | 0.4278 | 0.2491 | 0.9614 | 1608 | 1568 | 422 |
| `qp_slope_sd_mk_pval` | Q-P Seasonality | mk_pval | -0.2513 | 0.4291 | 0.2521 | 0.9946 | 1608 | 1568 | 422 |
| `avg_storage_mk_pval` | Storage | mk_pval | -0.2484 | 0.3914 | 0.2558 | 0.9798 | 1627 | 1583 | 422 |
| `spring_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.2339 | 0.4517 | 0.2482 | 0.9842 | 1612 | 1566 | 420 |
| `Qwin_senn_slp` | Flow Volumes | senn_slp | -0.2284 | 0.7764 | 1.4426 | 3283.7500 | 645 | 544 | 513 |
| `spring_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.2213 | 0.4543 | 0.2499 | 0.9931 | 1612 | 1566 | 420 |
| `Qann_senn_slp` | Flow Volumes | senn_slp | -0.2101 | 0.7668 | 5.8436 | 14990.1222 | 645 | 544 | 513 |
| `D99_day_mk_pval` | Flow Timing | mk_pval | -0.2098 | 0.4653 | 0.2286 | 0.9976 | 653 | 554 | 509 |
| `fall_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.2048 | 0.4663 | 0.2571 | 0.9994 | 1611 | 1566 | 421 |
| `D99_day_spearman_pval` | Flow Timing | spearman_pval | -0.2047 | 0.4669 | 0.2243 | 0.9828 | 653 | 554 | 509 |
| `fall_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.1932 | 0.4704 | 0.2536 | 0.9962 | 1611 | 1566 | 421 |
| `D20_day_spearman_pval` | Flow Timing | spearman_pval | -0.1773 | 0.4306 | 0.2131 | 0.9800 | 653 | 554 | 509 |
| `winter_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.1769 | 0.5071 | 0.2468 | 0.9945 | 1609 | 1564 | 417 |
| `D20_day_mk_pval` | Flow Timing | mk_pval | -0.1709 | 0.4365 | 0.2144 | 0.9628 | 653 | 554 | 509 |
| `winter_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.1465 | 0.5184 | 0.2467 | 0.9968 | 1609 | 1564 | 417 |
| `Q95_Q10_spearman_pval` | Flow Percentiles | spearman_pval | -0.1395 | 0.4808 | 0.2201 | 0.9854 | 645 | 544 | 513 |
| `Q95_spearman_pval` | Flow Percentiles | spearman_pval | -0.1295 | 0.4874 | 0.2187 | 0.9928 | 645 | 544 | 513 |
| `D10_day_spearman_pval` | Flow Timing | spearman_pval | -0.1214 | 0.4523 | 0.2115 | 0.9928 | 653 | 554 | 509 |
| `Q90_spearman_pval` | Flow Percentiles | spearman_pval | -0.1183 | 0.4906 | 0.2192 | 0.9856 | 645 | 544 | 513 |
| `Q90_mk_pval` | Flow Percentiles | mk_pval | -0.1183 | 0.4858 | 0.2227 | 0.9939 | 645 | 544 | 513 |
| `Dmax_spearman_pval` | Flow Timing | spearman_pval | -0.1166 | 0.4330 | 0.2137 | 0.9788 | 653 | 554 | 509 |
| `Dmax_mk_pval` | Flow Timing | mk_pval | -0.1147 | 0.4333 | 0.2157 | 0.9727 | 653 | 554 | 509 |
| `Q95_Q10_mk_pval` | Flow Percentiles | mk_pval | -0.1138 | 0.4839 | 0.2222 | 0.9835 | 645 | 544 | 513 |
| `Q95_mk_pval` | Flow Percentiles | mk_pval | -0.1125 | 0.4897 | 0.2208 | 0.9835 | 645 | 544 | 513 |
| `Q99_spearman_pval` | Flow Percentiles | spearman_pval | -0.1115 | 0.4864 | 0.2177 | 0.9889 | 645 | 544 | 513 |
| `D10_day_mk_pval` | Flow Timing | mk_pval | -0.1011 | 0.4636 | 0.2127 | 0.9576 | 653 | 554 | 509 |
| `D5_day_spearman_pval` | Flow Timing | spearman_pval | -0.1003 | 0.4612 | 0.2079 | 0.9631 | 653 | 554 | 509 |
| `D5_day_mk_pval` | Flow Timing | mk_pval | -0.1000 | 0.4604 | 0.2106 | 0.9771 | 653 | 554 | 509 |
| `D30_day_mk_pval` | Flow Timing | mk_pval | -0.0966 | 0.4693 | 0.2052 | 0.9521 | 653 | 554 | 509 |
| `annual_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.0960 | 0.5295 | 0.2394 | 0.9970 | 1607 | 1563 | 418 |
| `D50_day_spearman_pval` | Flow Timing | spearman_pval | -0.0931 | 0.4737 | 0.2060 | 0.9357 | 653 | 554 | 509 |
| `Q99_mk_pval` | Flow Percentiles | mk_pval | -0.0924 | 0.4917 | 0.2190 | 0.9902 | 645 | 544 | 513 |
| `D50_day_mk_pval` | Flow Timing | mk_pval | -0.0873 | 0.4767 | 0.2067 | 0.9652 | 653 | 554 | 509 |
| `annual_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.0869 | 0.5344 | 0.2346 | 0.9957 | 1607 | 1563 | 418 |
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.0827 | 0.5045 | 0.2116 | 0.9796 | 4364 | 4141 | 755 |
| `dur_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.0809 | 0.4828 | 0.2167 | 0.9789 | 1487 | 1352 | 673 |
| `FDC90th_spearman_pval` | FDC | spearman_pval | -0.0806 | 0.5250 | 0.2184 | 0.9944 | 690 | 608 | 514 |
| `Qwin_linear_slp` | Flow Volumes | linear_slp | -0.0793 | 0.8129 | 1.3998 | 2572.4661 | 645 | 544 | 513 |
| `D60_day_mk_pval` | Flow Timing | mk_pval | -0.0767 | 0.4844 | 0.2079 | 0.9975 | 653 | 554 | 509 |
| `Qsum_mk_pval` | Flow Volumes | mk_pval | -0.0765 | 0.5220 | 0.2125 | 0.9965 | 646 | 545 | 513 |
| `Qsum_spearman_pval` | Flow Volumes | spearman_pval | -0.0761 | 0.5243 | 0.2103 | 0.9994 | 645 | 544 | 513 |
| `D60_day_spearman_pval` | Flow Timing | spearman_pval | -0.0759 | 0.4841 | 0.2058 | 0.9651 | 653 | 554 | 509 |
| `FDC90th_mk_pval` | FDC | mk_pval | -0.0753 | 0.5245 | 0.2208 | 0.9972 | 690 | 608 | 514 |
| `dur_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.0752 | 0.5067 | 0.2192 | 0.9880 | 1186 | 1058 | 550 |
| `dur_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0740 | 0.4861 | 0.2183 | 0.9968 | 1487 | 1352 | 673 |
| `D30_day_spearman_pval` | Flow Timing | spearman_pval | -0.0731 | 0.4780 | 0.2018 | 0.9597 | 653 | 554 | 509 |
| `D40_day_spearman_pval` | Flow Timing | spearman_pval | -0.0716 | 0.4856 | 0.2034 | 0.9458 | 653 | 554 | 509 |
| `n_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.0706 | 0.5215 | 0.2162 | 0.9899 | 688 | 606 | 514 |
| `D40_day_mk_pval` | Flow Timing | mk_pval | -0.0700 | 0.4826 | 0.2064 | 0.9694 | 653 | 554 | 509 |
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0700 | 0.5114 | 0.2140 | 0.9792 | 4364 | 4141 | 755 |
| `n_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.0684 | 0.5215 | 0.2194 | 0.9920 | 688 | 606 | 514 |
| `n_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0640 | 0.5118 | 0.2190 | 0.9998 | 645 | 544 | 513 |
| `dur_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.0638 | 0.5124 | 0.2219 | 0.9975 | 1186 | 1058 | 550 |
| `D25_to_D75_spearman_pval` | Flow Timing | spearman_pval | -0.0564 | 0.4964 | 0.2100 | 0.9799 | 653 | 554 | 509 |
| `n_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.0551 | 0.5173 | 0.2154 | 0.9700 | 645 | 544 | 513 |
| `D70_day_spearman_pval` | Flow Timing | spearman_pval | -0.0490 | 0.4952 | 0.2063 | 0.9484 | 653 | 554 | 509 |
| `Q95_Q10_linear_slp` | Flow Percentiles | linear_slp | -0.0468 | 0.7536 | 0.0233 | 20.6690 | 645 | 544 | 513 |
| `D25_to_D75_mk_pval` | Flow Timing | mk_pval | -0.0437 | 0.5031 | 0.2114 | 0.9914 | 653 | 554 | 509 |
| `D70_day_mk_pval` | Flow Timing | mk_pval | -0.0422 | 0.4979 | 0.2065 | 0.9679 | 653 | 554 | 509 |
| `BFI_Eckhardt_spearman_pval` | Baseflow | spearman_pval | -0.0416 | 0.5475 | 0.2139 | 0.9912 | 653 | 554 | 509 |
| `BFI_Eckhardt_mk_pval` | Baseflow | mk_pval | -0.0415 | 0.5451 | 0.2169 | 0.9980 | 653 | 554 | 509 |
| `Qann_mk_pval` | Flow Volumes | mk_pval | -0.0408 | 0.5184 | 0.2166 | 0.9986 | 645 | 544 | 513 |
| `D1_day_mk_pval` | Flow Timing | mk_pval | -0.0394 | 0.4986 | 0.2038 | 0.9579 | 653 | 554 | 509 |
| `Qann_spearman_pval` | Flow Volumes | spearman_pval | -0.0390 | 0.5219 | 0.2148 | 0.9974 | 645 | 544 | 513 |
| `D1_day_spearman_pval` | Flow Timing | spearman_pval | -0.0294 | 0.5028 | 0.2000 | 0.9455 | 653 | 554 | 509 |
| `Qspr_spearman_pval` | Flow Volumes | spearman_pval | -0.0257 | 0.5400 | 0.2105 | 0.9620 | 645 | 544 | 513 |
| `D95_day_spearman_pval` | Flow Timing | spearman_pval | -0.0251 | 0.5386 | 0.2049 | 0.9599 | 653 | 554 | 509 |
| `Q80_mk_pval` | Flow Percentiles | mk_pval | -0.0240 | 0.5377 | 0.2157 | 0.9978 | 645 | 544 | 513 |
| `Q80_spearman_pval` | Flow Percentiles | spearman_pval | -0.0220 | 0.5431 | 0.2122 | 0.9956 | 645 | 544 | 513 |
| `flashinessRB_spearman_pval` | Flashiness | spearman_pval | -0.0207 | 0.5685 | 0.2089 | 0.9810 | 653 | 554 | 509 |
| `flashinessRB_mk_pval` | Flashiness | mk_pval | -0.0198 | 0.5655 | 0.2106 | 0.9973 | 653 | 554 | 509 |
| `Qsum_linear_slp` | Flow Volumes | linear_slp | -0.0193 | 0.7765 | 1.8591 | 2872.1548 | 645 | 544 | 513 |
| `n_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0176 | 0.5763 | 0.2082 | 0.9973 | 979 | 889 | 512 |
| `D95_day_mk_pval` | Flow Timing | mk_pval | -0.0163 | 0.5458 | 0.2063 | 0.9937 | 653 | 554 | 509 |
| `Q95_linear_slp` | Flow Percentiles | linear_slp | -0.0155 | 0.7515 | 0.0298 | 32.3488 | 645 | 544 | 513 |
| `b_events_spearman_pval` | Recession | spearman_pval | -0.0149 | 0.5467 | 0.2095 | 0.9795 | 1068 | 1232 | 164 |
| `Qspr_mk_pval` | Flow Volumes | mk_pval | -0.0149 | 0.5422 | 0.2114 | 0.9901 | 645 | 544 | 513 |
| `Q5_mk_pval` | Flow Percentiles | mk_pval | -0.0070 | 0.6023 | 0.2060 | 0.9997 | 698 | 620 | 510 |
| `Q75_mk_pval` | Flow Percentiles | mk_pval | -0.0057 | 0.5531 | 0.2118 | 0.9933 | 645 | 544 | 513 |
| `Q60_senn_slp` | Flow Percentiles | senn_slp | -0.0053 | 0.7847 | 0.0152 | 40.3571 | 645 | 544 | 513 |
| `b_events_mk_pval` | Recession | mk_pval | -0.0049 | 0.5502 | 0.2116 | 0.9963 | 1068 | 1232 | 164 |
| `FDCall_spearman_pval` | FDC | spearman_pval | -0.0044 | 0.5585 | 0.2099 | 0.9920 | 645 | 544 | 513 |
| `Qfal_mk_pval` | Flow Volumes | mk_pval | -0.0038 | 0.5520 | 0.2174 | 0.9902 | 645 | 544 | 513 |
| `Q1_mk_pval` | Flow Percentiles | mk_pval | -0.0019 | 0.6033 | 0.2038 | 0.9894 | 715 | 641 | 510 |
| `Q1_spearman_pval` | Flow Percentiles | spearman_pval | -0.0009 | 0.6060 | 0.2010 | 0.9960 | 715 | 641 | 510 |
| `Q5_spearman_pval` | Flow Percentiles | spearman_pval | 0.0001 | 0.6087 | 0.2020 | 0.9938 | 698 | 620 | 510 |
| `FDCall_mk_pval` | FDC | mk_pval | 0.0028 | 0.5591 | 0.2121 | 0.9972 | 645 | 544 | 513 |
| `n_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | 0.0035 | 0.5855 | 0.2029 | 0.9980 | 979 | 889 | 512 |
| `log_a_events_spearman_pval` | Recession | spearman_pval | 0.0039 | 0.5429 | 0.2075 | 0.9910 | 1068 | 1232 | 164 |
| `BFI_LyneHollick_mk_pval` | Baseflow | mk_pval | 0.0041 | 0.5686 | 0.2131 | 0.9867 | 653 | 554 | 509 |
| `Q75_spearman_pval` | Flow Percentiles | spearman_pval | 0.0041 | 0.5589 | 0.2084 | 0.9887 | 645 | 544 | 513 |
| `Qfal_spearman_pval` | Flow Volumes | spearman_pval | 0.0067 | 0.5563 | 0.2145 | 0.9777 | 645 | 544 | 513 |
| `BFI_LyneHollick_spearman_pval` | Baseflow | spearman_pval | 0.0069 | 0.5715 | 0.2109 | 0.9921 | 653 | 554 | 509 |
| `n_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.0135 | 0.5388 | 0.2077 | 0.9686 | 645 | 544 | 513 |
| `FDCmid_mk_pval` | FDC | mk_pval | 0.0163 | 0.5540 | 0.2082 | 0.9834 | 645 | 544 | 513 |
| `log_a_events_mk_pval` | Recession | mk_pval | 0.0169 | 0.5478 | 0.2092 | 0.9867 | 1068 | 1232 | 164 |
| `dur_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.0173 | 0.5398 | 0.2069 | 0.9939 | 655 | 557 | 508 |
| `FDCmid_spearman_pval` | FDC | spearman_pval | 0.0179 | 0.5539 | 0.2056 | 0.9838 | 645 | 544 | 513 |
| `dur_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.0223 | 0.5417 | 0.2092 | 0.9820 | 655 | 557 | 508 |
| `Q10_mk_pval` | Flow Percentiles | mk_pval | 0.0234 | 0.6052 | 0.2052 | 0.9999 | 687 | 605 | 514 |
| `n_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.0254 | 0.5427 | 0.2100 | 0.9822 | 645 | 544 | 513 |
| `Flow_Reversals_summer_mk_pval` | Pulse Metrics | mk_pval | 0.0276 | 0.6141 | 0.1999 | 1.0000 | 645 | 544 | 513 |
| `Q10_spearman_pval` | Flow Percentiles | spearman_pval | 0.0293 | 0.6118 | 0.2022 | 0.9866 | 687 | 605 | 514 |
| `Flow_Reversals_summer_spearman_pval` | Pulse Metrics | spearman_pval | 0.0341 | 0.6153 | 0.1970 | 0.9921 | 645 | 544 | 513 |
| `Q70_spearman_pval` | Flow Percentiles | spearman_pval | 0.0345 | 0.5730 | 0.2048 | 0.9884 | 646 | 545 | 513 |
| `negative_ann_spearman_pval` | Negative Flow | spearman_pval | 0.0345 | 0.6108 | 0.1504 | 0.8662 | 6634 | 6629 | 5 |
| `TQmean_spearman_pval` | Pulse Metrics | spearman_pval | 0.0363 | 0.5536 | 0.2036 | 0.9880 | 645 | 544 | 513 |
| `Q70_mk_pval` | Flow Percentiles | mk_pval | 0.0366 | 0.5708 | 0.2075 | 0.9991 | 646 | 545 | 513 |
| `concavity_spearman_pval` | Recession | spearman_pval | 0.0373 | 0.5607 | 0.2040 | 0.9879 | 1068 | 1232 | 164 |
| `D80_day_spearman_pval` | Flow Timing | spearman_pval | 0.0382 | 0.5459 | 0.1960 | 0.9692 | 653 | 554 | 509 |
| `Flow_Reversals_winter_mk_pval` | Pulse Metrics | mk_pval | 0.0390 | 0.6028 | 0.1976 | 0.9996 | 647 | 549 | 516 |
| `Q20_mk_pval` | Flow Percentiles | mk_pval | 0.0452 | 0.6083 | 0.2056 | 0.9972 | 671 | 578 | 513 |
| `TQmean_mk_pval` | Pulse Metrics | mk_pval | 0.0462 | 0.5573 | 0.2064 | 0.9816 | 645 | 544 | 513 |
| `b_pointcloud_spearman_pval` | Recession | spearman_pval | 0.0480 | 0.5506 | 0.1965 | 0.9824 | 1249 | 1358 | 109 |
| `concavity_mk_pval` | Recession | mk_pval | 0.0493 | 0.5666 | 0.2051 | 0.9902 | 1068 | 1232 | 164 |
| `D80_day_mk_pval` | Flow Timing | mk_pval | 0.0512 | 0.5513 | 0.1972 | 0.9805 | 653 | 554 | 509 |
| `Flow_Reversals_winter_spearman_pval` | Pulse Metrics | spearman_pval | 0.0516 | 0.6082 | 0.1953 | 0.9939 | 647 | 549 | 516 |
| `Qwin_spearman_pval` | Flow Volumes | spearman_pval | 0.0540 | 0.5814 | 0.2061 | 0.9911 | 645 | 544 | 513 |
| `Q20_spearman_pval` | Flow Percentiles | spearman_pval | 0.0595 | 0.6164 | 0.2023 | 0.9919 | 671 | 578 | 513 |
| `Q60_spearman_pval` | Flow Percentiles | spearman_pval | 0.0623 | 0.5921 | 0.2017 | 0.9981 | 646 | 545 | 513 |
| `Q99_linear_slp` | Flow Percentiles | linear_slp | 0.0629 | 0.7744 | 0.0468 | 22.7524 | 645 | 544 | 513 |
| `Qwin_mk_pval` | Flow Volumes | mk_pval | 0.0699 | 0.5863 | 0.2069 | 0.9874 | 645 | 544 | 513 |
| `Q25_mk_pval` | Flow Percentiles | mk_pval | 0.0703 | 0.6144 | 0.2018 | 0.9997 | 669 | 576 | 515 |
| `Q60_mk_pval` | Flow Percentiles | mk_pval | 0.0704 | 0.5919 | 0.2038 | 0.9905 | 646 | 545 | 513 |
| `D90_day_spearman_pval` | Flow Timing | spearman_pval | 0.0731 | 0.5759 | 0.1952 | 0.9696 | 653 | 554 | 509 |
| `Q25_spearman_pval` | Flow Percentiles | spearman_pval | 0.0746 | 0.6190 | 0.2001 | 0.9925 | 669 | 576 | 515 |
| `Flow_Reversals_spring_mk_pval` | Pulse Metrics | mk_pval | 0.0748 | 0.6068 | 0.2013 | 0.9991 | 646 | 546 | 512 |
| `negative_ann_mk_pval` | Negative Flow | mk_pval | 0.0779 | 0.6021 | 0.1547 | 0.8502 | 6634 | 6629 | 5 |
| `D90_day_mk_pval` | Flow Timing | mk_pval | 0.0805 | 0.5778 | 0.1953 | 0.9891 | 653 | 554 | 509 |
| `Flow_Reversals_spring_spearman_pval` | Pulse Metrics | spearman_pval | 0.0842 | 0.6113 | 0.1994 | 0.9675 | 646 | 546 | 512 |
| `log_a_pointcloud_spearman_pval` | Recession | spearman_pval | 0.0882 | 0.5780 | 0.1931 | 0.9742 | 1249 | 1358 | 109 |
| `Q30_mk_pval` | Flow Percentiles | mk_pval | 0.0968 | 0.6245 | 0.1982 | 0.9859 | 663 | 566 | 513 |
| `Q30_spearman_pval` | Flow Percentiles | spearman_pval | 0.0968 | 0.6274 | 0.1952 | 0.9866 | 663 | 566 | 513 |
| `b_pointcloud_mk_pval` | Recession | mk_pval | 0.1070 | 0.5777 | 0.1982 | 0.9813 | 1249 | 1358 | 109 |
| `Flow_Reversals_fall_mk_pval` | Pulse Metrics | mk_pval | 0.1100 | 0.6670 | 0.1803 | 0.9997 | 645 | 544 | 513 |
| `Flow_Reversals_fall_spearman_pval` | Pulse Metrics | spearman_pval | 0.1115 | 0.6675 | 0.1777 | 0.9882 | 645 | 544 | 513 |
| `Q50_spearman_pval` | Flow Percentiles | spearman_pval | 0.1262 | 0.6259 | 0.1912 | 0.9917 | 651 | 551 | 512 |
| `Q90_linear_slp` | Flow Percentiles | linear_slp | 0.1308 | 0.7494 | 0.0247 | 40.2369 | 645 | 544 | 513 |
| `log_a_pointcloud_mk_pval` | Recession | mk_pval | 0.1392 | 0.6016 | 0.1973 | 0.9847 | 1249 | 1358 | 109 |
| `Flow_Reversals_annual_mk_pval` | Pulse Metrics | mk_pval | 0.1404 | 0.6909 | 0.1678 | 0.9986 | 645 | 544 | 513 |
| `Q50_mk_pval` | Flow Percentiles | mk_pval | 0.1405 | 0.6277 | 0.1928 | 0.9984 | 651 | 551 | 512 |
| `n_recession_events_mk_pval` | Recession | mk_pval | 0.1418 | 0.6730 | 0.1841 | 0.9826 | 55 | 65 | 10 |
| `Q40_spearman_pval` | Flow Percentiles | spearman_pval | 0.1459 | 0.6411 | 0.1908 | 0.9831 | 657 | 558 | 513 |
| `n_recession_events_spearman_pval` | Recession | spearman_pval | 0.1473 | 0.6745 | 0.1808 | 0.9955 | 55 | 65 | 10 |
| `Flow_Reversals_annual_spearman_pval` | Pulse Metrics | spearman_pval | 0.1475 | 0.6945 | 0.1648 | 0.9956 | 645 | 544 | 513 |
| `Qsum_senn_slp` | Flow Volumes | senn_slp | 0.1489 | 0.7709 | 1.7799 | 2665.1399 | 645 | 544 | 513 |
| `Q40_mk_pval` | Flow Percentiles | mk_pval | 0.1525 | 0.6406 | 0.1922 | 0.9996 | 657 | 558 | 513 |
| `Q70_linear_slp` | Flow Percentiles | linear_slp | 0.1545 | 0.8002 | 0.0138 | 32.4460 | 645 | 544 | 513 |
| `Q99_senn_slp` | Flow Percentiles | senn_slp | 0.1861 | 0.7535 | 0.0427 | 25.0000 | 645 | 544 | 513 |
| `Q90_senn_slp` | Flow Percentiles | senn_slp | 0.1869 | 0.7403 | 0.0254 | 51.3101 | 645 | 544 | 513 |
| `Qspr_linear_slp` | Flow Volumes | linear_slp | 0.1933 | 0.7723 | 1.5415 | 3150.5525 | 645 | 544 | 513 |
| `Q95_senn_slp` | Flow Percentiles | senn_slp | 0.1975 | 0.7452 | 0.0310 | 41.8500 | 645 | 544 | 513 |
| `elasticity_rolling_spearman_rho` | Elasticity | spearman_rho | 0.2022 | 0.5520 | 0.3306 | 1.4199 | 1165 | 1165 | 0 |
| `Qspr_senn_slp` | Flow Volumes | senn_slp | 0.2046 | 0.7649 | 1.7509 | 4849.7345 | 645 | 544 | 513 |
| `Q80_linear_slp` | Flow Percentiles | linear_slp | 0.2074 | 0.7867 | 0.0174 | 37.0060 | 645 | 544 | 513 |
| `elasticity_rolling_mk_rho` | Elasticity | mk_rho | 0.2354 | 0.5548 | 0.2420 | 1.1086 | 1165 | 1165 | 0 |
| `annual_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.2389 | 0.7683 | 0.0020 | 2.0050 | 1607 | 1563 | 418 |
| `Q75_linear_slp` | Flow Percentiles | linear_slp | 0.2408 | 0.7990 | 0.0149 | 35.2612 | 645 | 544 | 513 |
| `Q95_Q10_senn_slp` | Flow Percentiles | senn_slp | 0.2606 | 0.7498 | 0.0224 | 21.3900 | 645 | 544 | 513 |
| `Q80_senn_slp` | Flow Percentiles | senn_slp | 0.2802 | 0.7684 | 0.0179 | 41.8299 | 645 | 544 | 513 |
| `elasticity_rolling_senn_slp` | Elasticity | senn_slp | 0.2951 | 0.5226 | 0.0318 | 0.4640 | 1165 | 1165 | 0 |
| `elasticity_rolling_linear_slp` | Elasticity | linear_slp | 0.3319 | 0.5426 | 0.0332 | 0.4354 | 1165 | 1165 | 0 |
| `fall_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.3405 | 0.7388 | 0.0041 | 3.7436 | 1611 | 1566 | 421 |
| `annual_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.3613 | 0.7825 | 0.0018 | 1.2977 | 1607 | 1563 | 418 |
| `elasticity_annual_spearman_rho` | Elasticity | spearman_rho | 0.3622 | 0.6273 | 0.1225 | 0.6638 | 1165 | 1165 | 0 |
| `elasticity_annual_mk_rho` | Elasticity | mk_rho | 0.3725 | 0.6290 | 0.0846 | 0.4223 | 1165 | 1165 | 0 |
| `Q75_senn_slp` | Flow Percentiles | senn_slp | 0.3765 | 0.7801 | 0.0164 | 41.5000 | 645 | 544 | 513 |
| `Q70_senn_slp` | Flow Percentiles | senn_slp | 0.3780 | 0.7764 | 0.0159 | 42.7495 | 645 | 544 | 513 |
| `elasticity_annual_senn_slp` | Elasticity | senn_slp | 0.3982 | 0.6031 | 0.0333 | 0.4016 | 1165 | 1165 | 0 |
| `summer_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.4012 | 0.6891 | 0.1199 | 0.6362 | 1786 | 1762 | 420 |
| `summer_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.4162 | 0.6908 | 0.0821 | 0.4238 | 1786 | 1762 | 420 |
| `avg_storage_senn_slp` | Storage | senn_slp | 0.4182 | 0.6704 | 2.6698 | 1403.6658 | 1627 | 1583 | 422 |
| `qp_bimodality_spearman_rho` | Q-P Seasonality | spearman_rho | 0.4253 | 0.6856 | 0.1105 | 0.5937 | 1615 | 1576 | 423 |
| `fall_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.4381 | 0.6903 | 0.0026 | 2.4167 | 1611 | 1566 | 421 |
| `qp_bimodality_mk_rho` | Q-P Seasonality | mk_rho | 0.4419 | 0.6902 | 0.0760 | 0.3817 | 1615 | 1576 | 423 |
| `avg_storage_spearman_rho` | Storage | spearman_rho | 0.4575 | 0.6768 | 0.1153 | 0.7562 | 1627 | 1583 | 422 |
| `D20_day_spearman_rho` | Flow Timing | spearman_rho | 0.4621 | 0.7020 | 0.0847 | 0.5317 | 653 | 554 | 509 |
| `D99_day_spearman_rho` | Flow Timing | spearman_rho | 0.4622 | 0.7308 | 0.0961 | 0.8971 | 653 | 554 | 509 |
| `avg_storage_mk_rho` | Storage | mk_rho | 0.4736 | 0.6795 | 0.0797 | 0.4969 | 1627 | 1583 | 422 |
| `D99_day_mk_rho` | Flow Timing | mk_rho | 0.4749 | 0.7354 | 0.0682 | 0.5850 | 653 | 554 | 509 |
| `D20_day_mk_rho` | Flow Timing | mk_rho | 0.4824 | 0.7101 | 0.0576 | 0.3671 | 653 | 554 | 509 |
| `Dmax_spearman_rho` | Flow Timing | spearman_rho | 0.4827 | 0.6911 | 0.0875 | 0.5226 | 653 | 554 | 509 |
| `D50_day_spearman_rho` | Flow Timing | spearman_rho | 0.4902 | 0.7156 | 0.0838 | 0.5642 | 653 | 554 | 509 |
| `D40_day_spearman_rho` | Flow Timing | spearman_rho | 0.4916 | 0.7191 | 0.0818 | 0.5760 | 653 | 554 | 509 |
| `Dmax_mk_rho` | Flow Timing | mk_rho | 0.4926 | 0.6932 | 0.0601 | 0.3662 | 653 | 554 | 509 |
| `D40_day_mk_rho` | Flow Timing | mk_rho | 0.4942 | 0.7155 | 0.0567 | 0.4243 | 653 | 554 | 509 |
| `qp_slope_sd_linear_slp` | Q-P Seasonality | linear_slp | 0.4949 | 0.5711 | 0.1193 | 11.6260 | 1608 | 1568 | 422 |
| `fall_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.4955 | 0.7316 | 0.1233 | 0.8972 | 1611 | 1566 | 421 |
| `D50_day_mk_rho` | Flow Timing | mk_rho | 0.4959 | 0.7140 | 0.0576 | 0.3867 | 653 | 554 | 509 |
| `dur_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.5005 | 0.7312 | 0.0150 | 0.8182 | 655 | 557 | 508 |
| `D30_day_spearman_rho` | Flow Timing | spearman_rho | 0.5023 | 0.7296 | 0.0786 | 0.4901 | 653 | 554 | 509 |
| `qp_bimodality_senn_slp` | Q-P Seasonality | senn_slp | 0.5040 | 0.6972 | 0.0014 | 0.0092 | 1615 | 1576 | 423 |
| `fall_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5044 | 0.7286 | 0.0851 | 0.5975 | 1611 | 1566 | 421 |
| `D30_day_mk_rho` | Flow Timing | mk_rho | 0.5062 | 0.7298 | 0.0544 | 0.3708 | 653 | 554 | 509 |
| `D10_day_spearman_rho` | Flow Timing | spearman_rho | 0.5101 | 0.7232 | 0.0833 | 0.6612 | 653 | 554 | 509 |
| `qp_bimodality_linear_slp` | Q-P Seasonality | linear_slp | 0.5106 | 0.7030 | 0.0014 | 0.0105 | 1615 | 1576 | 423 |
| `qp_slope_sd_spearman_rho` | Q-P Seasonality | spearman_rho | 0.5150 | 0.7210 | 0.1241 | 0.5930 | 1608 | 1568 | 422 |
| `D60_day_spearman_rho` | Flow Timing | spearman_rho | 0.5188 | 0.7333 | 0.0833 | 0.5444 | 653 | 554 | 509 |
| `D5_day_spearman_rho` | Flow Timing | spearman_rho | 0.5214 | 0.7242 | 0.0832 | 0.6196 | 653 | 554 | 509 |
| `D10_day_mk_rho` | Flow Timing | mk_rho | 0.5261 | 0.7289 | 0.0572 | 0.4028 | 653 | 554 | 509 |
| `qp_slope_sd_mk_rho` | Q-P Seasonality | mk_rho | 0.5264 | 0.7245 | 0.0855 | 0.3987 | 1608 | 1568 | 422 |
| `D60_day_mk_rho` | Flow Timing | mk_rho | 0.5289 | 0.7348 | 0.0574 | 0.3839 | 653 | 554 | 509 |
| `D5_day_mk_rho` | Flow Timing | mk_rho | 0.5330 | 0.7258 | 0.0578 | 0.4139 | 653 | 554 | 509 |
| `winter_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5449 | 0.7519 | 0.1255 | 0.8049 | 1609 | 1564 | 417 |
| `dur_low_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.5455 | 0.7408 | 0.1099 | 0.7304 | 1186 | 1058 | 550 |
| `dur_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.5460 | 0.7115 | 0.0400 | 1.7500 | 1487 | 1352 | 673 |
| `Q99_spearman_rho` | Flow Percentiles | spearman_rho | 0.5542 | 0.7673 | 0.0981 | 0.6567 | 645 | 544 | 513 |
| `D1_day_spearman_rho` | Flow Timing | spearman_rho | 0.5564 | 0.7418 | 0.0823 | 0.6098 | 653 | 554 | 509 |
| `D95_day_spearman_rho` | Flow Timing | spearman_rho | 0.5580 | 0.7687 | 0.0871 | 0.7388 | 653 | 554 | 509 |
| `D70_day_spearman_rho` | Flow Timing | spearman_rho | 0.5586 | 0.7555 | 0.0824 | 0.5860 | 653 | 554 | 509 |
| `winter_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5587 | 0.7554 | 0.0863 | 0.5450 | 1609 | 1564 | 417 |
| `dur_low_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5589 | 0.7496 | 0.1114 | 0.8850 | 4364 | 4141 | 755 |
| `dur_low_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.5612 | 0.7451 | 0.0778 | 0.5009 | 1186 | 1058 | 550 |
| `spring_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5615 | 0.7669 | 0.1198 | 0.6438 | 1612 | 1566 | 420 |
| `log_a_seasonality_minimum_first_half` | Recession | scalar | 0.5637 | 0.8049 | 19.2114 | 361.1328 | 1069 | 1234 | 165 |
| `Q95_Q10_spearman_rho` | Flow Percentiles | spearman_rho | 0.5658 | 0.7726 | 0.0986 | 0.6229 | 645 | 544 | 513 |
| `D25_to_D75_spearman_rho` | Flow Timing | spearman_rho | 0.5667 | 0.7667 | 0.0885 | 0.6135 | 653 | 554 | 509 |
| `Q90_spearman_rho` | Flow Percentiles | spearman_rho | 0.5671 | 0.7742 | 0.1013 | 0.7223 | 645 | 544 | 513 |
| `Q95_spearman_rho` | Flow Percentiles | spearman_rho | 0.5674 | 0.7727 | 0.0999 | 0.6783 | 645 | 544 | 513 |
| `Q99_mk_rho` | Flow Percentiles | mk_rho | 0.5693 | 0.7697 | 0.0676 | 0.4176 | 645 | 544 | 513 |
| `D1_day_mk_rho` | Flow Timing | mk_rho | 0.5704 | 0.7452 | 0.0585 | 0.4372 | 653 | 554 | 509 |
| `D70_day_mk_rho` | Flow Timing | mk_rho | 0.5708 | 0.7584 | 0.0567 | 0.3716 | 653 | 554 | 509 |
| `D95_day_mk_rho` | Flow Timing | mk_rho | 0.5723 | 0.7743 | 0.0601 | 0.4861 | 653 | 554 | 509 |
| `D25_to_D75_mk_rho` | Flow Timing | mk_rho | 0.5725 | 0.7668 | 0.0612 | 0.4248 | 653 | 554 | 509 |
| `dur_low_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5726 | 0.7511 | 0.0776 | 0.6852 | 4364 | 4141 | 755 |
| `n_low_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.5728 | 0.7557 | 0.1139 | 0.7950 | 688 | 606 | 514 |
| `dur_high_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5756 | 0.7579 | 0.0989 | 0.5810 | 1487 | 1352 | 673 |
| `negative_ann_linear_slp` | Negative Flow | linear_slp | 0.5763 | 0.7728 | 0.0003 | 0.3527 | 645 | 544 | 513 |
| `qp_slope_sd_senn_slp` | Q-P Seasonality | senn_slp | 0.5764 | 0.6398 | 0.0138 | 1.7796 | 1608 | 1568 | 422 |
| `spring_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5783 | 0.7729 | 0.0822 | 0.4432 | 1612 | 1566 | 420 |
| `dur_high_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5811 | 0.7600 | 0.0690 | 0.4106 | 1487 | 1352 | 673 |
| `D60_day_senn_slp` | Flow Timing | senn_slp | 0.5825 | 0.7079 | 0.2340 | 3.2725 | 653 | 554 | 509 |
| `Qann_spearman_rho` | Flow Volumes | spearman_rho | 0.5827 | 0.7924 | 0.1002 | 0.6248 | 645 | 544 | 513 |
| `n_low_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.5829 | 0.7575 | 0.0834 | 0.6395 | 688 | 606 | 514 |
| `FDC90th_spearman_rho` | FDC | spearman_rho | 0.5853 | 0.7667 | 0.1179 | 0.8541 | 690 | 608 | 514 |
| `Q95_Q10_mk_rho` | Flow Percentiles | mk_rho | 0.5866 | 0.7788 | 0.0681 | 0.4035 | 645 | 544 | 513 |
| `n_high_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5874 | 0.7834 | 0.0725 | 0.5558 | 645 | 544 | 513 |
| `n_high_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5892 | 0.7863 | 0.1002 | 0.7458 | 645 | 544 | 513 |
| `Q95_mk_rho` | Flow Percentiles | mk_rho | 0.5901 | 0.7784 | 0.0691 | 0.4400 | 645 | 544 | 513 |
| `dur_high_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.5902 | 0.7186 | 0.0349 | 0.6946 | 655 | 557 | 508 |
| `D5_day_senn_slp` | Flow Timing | senn_slp | 0.5909 | 0.7078 | 0.1906 | 4.3639 | 653 | 554 | 509 |
| `Qfal_spearman_rho` | Flow Volumes | spearman_rho | 0.5910 | 0.7522 | 0.1080 | 0.8724 | 645 | 544 | 513 |
| `elasticity_annual_linear_slp` | Elasticity | linear_slp | 0.5920 | 0.6110 | 0.3267 | 26.9528 | 1165 | 1165 | 0 |
| `Q90_mk_rho` | Flow Percentiles | mk_rho | 0.5923 | 0.7807 | 0.0703 | 0.5087 | 645 | 544 | 513 |
| `concavity_linear_slp` | Recession | linear_slp | 0.5923 | 0.7442 | 0.0279 | 0.9138 | 1068 | 1232 | 164 |
| `annual_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5943 | 0.7867 | 0.1254 | 0.8608 | 1607 | 1563 | 418 |
| `FDC90th_mk_rho` | FDC | mk_rho | 0.5966 | 0.7689 | 0.0821 | 0.5246 | 690 | 608 | 514 |
| `D50_day_senn_slp` | Flow Timing | senn_slp | 0.5977 | 0.7005 | 0.2361 | 3.0440 | 653 | 554 | 509 |
| `Qann_mk_rho` | Flow Volumes | mk_rho | 0.5995 | 0.7953 | 0.0694 | 0.4039 | 645 | 544 | 513 |
| `D70_day_senn_slp` | Flow Timing | senn_slp | 0.6011 | 0.7403 | 0.2316 | 4.0111 | 653 | 554 | 509 |
| `D80_day_spearman_rho` | Flow Timing | spearman_rho | 0.6016 | 0.7776 | 0.0805 | 0.6596 | 653 | 554 | 509 |
| `Qsum_spearman_rho` | Flow Volumes | spearman_rho | 0.6056 | 0.7847 | 0.0983 | 0.7338 | 645 | 544 | 513 |
| `BFI_Eckhardt_spearman_rho` | Baseflow | spearman_rho | 0.6066 | 0.7927 | 0.1072 | 0.9792 | 653 | 554 | 509 |
| `Q80_spearman_rho` | Flow Percentiles | spearman_rho | 0.6089 | 0.8022 | 0.1006 | 0.7460 | 645 | 544 | 513 |
| `annual_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.6090 | 0.7876 | 0.0878 | 0.5443 | 1607 | 1563 | 418 |
| `concavity_spearman_rho` | Recession | spearman_rho | 0.6101 | 0.7743 | 0.1030 | 0.6158 | 1068 | 1232 | 164 |
| `b_events_linear_slp` | Recession | linear_slp | 0.6106 | 0.7343 | 0.0181 | 0.9240 | 1068 | 1232 | 164 |
| `Qfal_mk_rho` | Flow Volumes | mk_rho | 0.6112 | 0.7529 | 0.0749 | 0.5377 | 645 | 544 | 513 |
| `b_events_spearman_rho` | Recession | spearman_rho | 0.6116 | 0.7784 | 0.1089 | 0.6756 | 1068 | 1232 | 164 |
| `log_a_events_spearman_rho` | Recession | spearman_rho | 0.6123 | 0.7793 | 0.1034 | 0.8533 | 1068 | 1232 | 164 |
| `D40_day_senn_slp` | Flow Timing | senn_slp | 0.6127 | 0.7098 | 0.2541 | 3.2889 | 653 | 554 | 509 |
| `flashinessRB_spearman_rho` | Flashiness | spearman_rho | 0.6139 | 0.8019 | 0.1148 | 0.8695 | 653 | 554 | 509 |
| `n_recession_events_mk_rho` | Recession | mk_rho | 0.6142 | 0.8061 | 0.0793 | 0.5869 | 55 | 65 | 10 |
| `BFI_LyneHollick_spearman_rho` | Baseflow | spearman_rho | 0.6148 | 0.8004 | 0.1080 | 0.8944 | 653 | 554 | 509 |
| `n_recession_events_spearman_rho` | Recession | spearman_rho | 0.6152 | 0.8093 | 0.1050 | 0.8510 | 55 | 65 | 10 |
| `D90_day_spearman_rho` | Flow Timing | spearman_rho | 0.6168 | 0.7851 | 0.0818 | 0.6579 | 653 | 554 | 509 |
| `D80_day_mk_rho` | Flow Timing | mk_rho | 0.6179 | 0.7832 | 0.0552 | 0.4318 | 653 | 554 | 509 |
| `Qsum_mk_rho` | Flow Volumes | mk_rho | 0.6185 | 0.7846 | 0.0677 | 0.4679 | 646 | 545 | 513 |
| `D60_day_linear_slp` | Flow Timing | linear_slp | 0.6186 | 0.7306 | 0.2320 | 3.3036 | 653 | 554 | 509 |
| `D10_day_senn_slp` | Flow Timing | senn_slp | 0.6198 | 0.7214 | 0.2453 | 2.8933 | 653 | 554 | 509 |
| `FDCmid_spearman_rho` | FDC | spearman_rho | 0.6199 | 0.7804 | 0.1000 | 0.7671 | 645 | 544 | 513 |
| `BFI_Eckhardt_mk_rho` | Baseflow | mk_rho | 0.6199 | 0.7955 | 0.0741 | 0.6716 | 653 | 554 | 509 |
| `n_low_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.6203 | 0.7908 | 0.1259 | 1.1851 | 979 | 889 | 512 |
| `D50_day_linear_slp` | Flow Timing | linear_slp | 0.6215 | 0.7087 | 0.2461 | 3.3354 | 653 | 554 | 509 |
| `n_low_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.6217 | 0.7888 | 0.0933 | 0.8463 | 979 | 889 | 512 |
| `concavity_mk_rho` | Recession | mk_rho | 0.6223 | 0.7779 | 0.0713 | 0.4163 | 1068 | 1232 | 164 |
| `flashinessRB_mk_rho` | Flashiness | mk_rho | 0.6233 | 0.8045 | 0.0798 | 0.6253 | 653 | 554 | 509 |
| `FDCall_spearman_rho` | FDC | spearman_rho | 0.6244 | 0.7849 | 0.1083 | 0.8746 | 645 | 544 | 513 |
| `b_events_mk_rho` | Recession | mk_rho | 0.6248 | 0.7822 | 0.0754 | 0.4705 | 1068 | 1232 | 164 |
| `log_a_events_mk_rho` | Recession | mk_rho | 0.6249 | 0.7822 | 0.0715 | 0.5748 | 1068 | 1232 | 164 |
| `Dmax_senn_slp` | Flow Timing | senn_slp | 0.6254 | 0.7074 | 0.4938 | 8.3080 | 653 | 554 | 509 |
| `BFI_LyneHollick_mk_rho` | Baseflow | mk_rho | 0.6261 | 0.8031 | 0.0749 | 0.6150 | 653 | 554 | 509 |
| `Q80_mk_rho` | Flow Percentiles | mk_rho | 0.6268 | 0.8038 | 0.0701 | 0.5182 | 645 | 544 | 513 |
| `Flow_Reversals_summer_spearman_rho` | Pulse Metrics | spearman_rho | 0.6283 | 0.8083 | 0.1164 | 0.8152 | 645 | 544 | 513 |
| `Q75_spearman_rho` | Flow Percentiles | spearman_rho | 0.6285 | 0.8140 | 0.1017 | 0.7774 | 645 | 544 | 513 |
| `b_pointcloud_spearman_rho` | Recession | spearman_rho | 0.6291 | 0.7853 | 0.1149 | 1.4455 | 1249 | 1358 | 109 |
| `dur_high_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.6311 | 0.7939 | 0.0935 | 0.6076 | 655 | 557 | 508 |
| `D80_day_senn_slp` | Flow Timing | senn_slp | 0.6318 | 0.7709 | 0.2313 | 3.4375 | 653 | 554 | 509 |
| `D90_day_mk_rho` | Flow Timing | mk_rho | 0.6318 | 0.7912 | 0.0562 | 0.4465 | 653 | 554 | 509 |
| `Flow_Reversals_fall_spearman_rho` | Pulse Metrics | spearman_rho | 0.6325 | 0.8129 | 0.1190 | 0.9492 | 645 | 544 | 513 |
| `FDCmid_mk_rho` | FDC | mk_rho | 0.6327 | 0.7818 | 0.0692 | 0.5080 | 645 | 544 | 513 |
| `Dmax_linear_slp` | Flow Timing | linear_slp | 0.6328 | 0.7179 | 0.5985 | 9.4977 | 653 | 554 | 509 |
| `D5_day_linear_slp` | Flow Timing | linear_slp | 0.6329 | 0.7169 | 0.2288 | 3.9052 | 653 | 554 | 509 |
| `D70_day_linear_slp` | Flow Timing | linear_slp | 0.6336 | 0.7630 | 0.2303 | 3.2429 | 653 | 554 | 509 |
| `concavity_senn_slp` | Recession | senn_slp | 0.6351 | 0.7440 | 0.0180 | 0.5672 | 1068 | 1232 | 164 |
| `n_high_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.6354 | 0.7961 | 0.0941 | 0.6214 | 645 | 544 | 513 |
| `Q70_spearman_rho` | Flow Percentiles | spearman_rho | 0.6362 | 0.8167 | 0.1030 | 0.8040 | 646 | 545 | 513 |
| `D40_day_linear_slp` | Flow Timing | linear_slp | 0.6369 | 0.7061 | 0.2651 | 3.0155 | 653 | 554 | 509 |
| `FDCall_mk_rho` | FDC | mk_rho | 0.6371 | 0.7865 | 0.0752 | 0.5447 | 645 | 544 | 513 |
| `TQmean_spearman_rho` | Pulse Metrics | spearman_rho | 0.6372 | 0.7927 | 0.0955 | 0.6937 | 645 | 544 | 513 |
| `dur_high_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.6377 | 0.7955 | 0.0670 | 0.4136 | 655 | 557 | 508 |
| `winter_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.6379 | 0.7377 | 0.0035 | 3.5605 | 1609 | 1564 | 417 |
| `dur_high_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.6392 | 0.7246 | 0.0546 | 1.8284 | 1487 | 1352 | 673 |
| `n_high_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.6399 | 0.7973 | 0.0686 | 0.4684 | 645 | 544 | 513 |
| `Flow_Reversals_summer_mk_rho` | Pulse Metrics | mk_rho | 0.6407 | 0.8120 | 0.0827 | 0.5916 | 645 | 544 | 513 |
| `D20_day_linear_slp` | Flow Timing | linear_slp | 0.6413 | 0.7031 | 0.3010 | 3.2940 | 653 | 554 | 509 |
| `Qwin_spearman_rho` | Flow Volumes | spearman_rho | 0.6424 | 0.8100 | 0.0990 | 0.9064 | 645 | 544 | 513 |
| `D10_day_linear_slp` | Flow Timing | linear_slp | 0.6428 | 0.7250 | 0.2757 | 3.6144 | 653 | 554 | 509 |
| `D1_day_senn_slp` | Flow Timing | senn_slp | 0.6438 | 0.6801 | 0.0847 | 2.9195 | 653 | 554 | 509 |
| `log_a_pointcloud_spearman_rho` | Recession | spearman_rho | 0.6439 | 0.7953 | 0.1155 | 1.4818 | 1249 | 1358 | 109 |
| `Q75_mk_rho` | Flow Percentiles | mk_rho | 0.6447 | 0.8150 | 0.0709 | 0.5212 | 645 | 544 | 513 |
| `Qspr_spearman_rho` | Flow Volumes | spearman_rho | 0.6448 | 0.8009 | 0.0943 | 0.6402 | 645 | 544 | 513 |
| `Flow_Reversals_fall_mk_rho` | Pulse Metrics | mk_rho | 0.6453 | 0.8161 | 0.0851 | 0.6473 | 645 | 544 | 513 |
| `D20_day_senn_slp` | Flow Timing | senn_slp | 0.6465 | 0.7131 | 0.2785 | 3.6327 | 653 | 554 | 509 |
| `Flow_Reversals_spring_spearman_rho` | Pulse Metrics | spearman_rho | 0.6467 | 0.8048 | 0.1095 | 0.7447 | 646 | 546 | 512 |
| `D30_day_senn_slp` | Flow Timing | senn_slp | 0.6468 | 0.7327 | 0.2650 | 3.1180 | 653 | 554 | 509 |
| `TQmean_mk_rho` | Pulse Metrics | mk_rho | 0.6472 | 0.7949 | 0.0660 | 0.4600 | 645 | 544 | 513 |
| `D80_day_linear_slp` | Flow Timing | linear_slp | 0.6475 | 0.7728 | 0.2336 | 3.4114 | 653 | 554 | 509 |
| `D30_day_linear_slp` | Flow Timing | linear_slp | 0.6492 | 0.7135 | 0.2832 | 3.0639 | 653 | 554 | 509 |
| `FDC90th_senn_slp` | FDC | senn_slp | 0.6505 | 0.7456 | 0.0190 | 1.9703 | 645 | 544 | 513 |
| `log_a_events_senn_slp` | Recession | senn_slp | 0.6505 | 0.7539 | 0.0176 | 0.7576 | 1068 | 1232 | 164 |
| `n_low_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.6508 | 0.7219 | 0.0232 | 0.6667 | 645 | 544 | 513 |
| `Q70_mk_rho` | Flow Percentiles | mk_rho | 0.6512 | 0.8171 | 0.0717 | 0.5478 | 646 | 545 | 513 |
| `b_events_senn_slp` | Recession | senn_slp | 0.6512 | 0.7532 | 0.0110 | 0.3449 | 1068 | 1232 | 164 |
| `Flow_Reversals_winter_spearman_rho` | Pulse Metrics | spearman_rho | 0.6517 | 0.8147 | 0.1220 | 0.8464 | 647 | 549 | 516 |
| `Flow_Reversals_spring_mk_rho` | Pulse Metrics | mk_rho | 0.6524 | 0.8059 | 0.0782 | 0.4926 | 646 | 546 | 512 |
| `Flow_Reversals_winter_senn_slp` | Pulse Metrics | senn_slp | 0.6534 | 0.8000 | 0.0812 | 1.2102 | 645 | 544 | 513 |
| `b_pointcloud_mk_rho` | Recession | mk_rho | 0.6545 | 0.7910 | 0.0816 | 1.2364 | 1249 | 1358 | 109 |
| `avg_storage_linear_slp` | Storage | linear_slp | 0.6546 | 0.6904 | 2.6342 | 1573.3937 | 1627 | 1583 | 422 |
| `Q60_spearman_rho` | Flow Percentiles | spearman_rho | 0.6564 | 0.8170 | 0.1040 | 0.7858 | 646 | 545 | 513 |
| `Flow_Reversals_winter_mk_rho` | Pulse Metrics | mk_rho | 0.6570 | 0.8150 | 0.0876 | 0.5784 | 647 | 549 | 516 |
| `D1_day_linear_slp` | Flow Timing | linear_slp | 0.6573 | 0.7049 | 0.1269 | 4.4434 | 653 | 554 | 509 |
| `Qwin_mk_rho` | Flow Volumes | mk_rho | 0.6607 | 0.8151 | 0.0683 | 0.5950 | 645 | 544 | 513 |
| `D95_day_senn_slp` | Flow Timing | senn_slp | 0.6616 | 0.7547 | 0.1959 | 2.7325 | 653 | 554 | 509 |
| `n_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.6625 | 0.7410 | 0.0187 | 0.4706 | 645 | 544 | 513 |
| `Qspr_mk_rho` | Flow Volumes | mk_rho | 0.6630 | 0.8059 | 0.0649 | 0.4148 | 645 | 544 | 513 |
| `Q1_spearman_rho` | Flow Percentiles | spearman_rho | 0.6631 | 0.8112 | 0.1276 | 1.2065 | 715 | 641 | 510 |
| `log_a_pointcloud_mk_rho` | Recession | mk_rho | 0.6633 | 0.8002 | 0.0826 | 1.2727 | 1249 | 1358 | 109 |
| `D25_to_D75_linear_slp` | Flow Timing | linear_slp | 0.6647 | 0.7594 | 0.3161 | 4.1164 | 653 | 554 | 509 |
| `D90_day_linear_slp` | Flow Timing | linear_slp | 0.6652 | 0.7717 | 0.2285 | 2.9680 | 653 | 554 | 509 |
| `D25_to_D75_senn_slp` | Flow Timing | senn_slp | 0.6665 | 0.7628 | 0.3046 | 4.5714 | 653 | 554 | 509 |
| `D90_day_senn_slp` | Flow Timing | senn_slp | 0.6674 | 0.7766 | 0.2229 | 2.7449 | 653 | 554 | 509 |
| `Q60_mk_rho` | Flow Percentiles | mk_rho | 0.6714 | 0.8169 | 0.0726 | 0.5332 | 646 | 545 | 513 |
| `Flow_Reversals_winter_linear_slp` | Pulse Metrics | linear_slp | 0.6718 | 0.8126 | 0.0830 | 0.9813 | 645 | 544 | 513 |
| `BFI_LyneHollick_linear_slp` | Baseflow | linear_slp | 0.6752 | 0.7959 | 0.0007 | 0.0206 | 653 | 554 | 509 |
| `BFI_LyneHollick_senn_slp` | Baseflow | senn_slp | 0.6762 | 0.7974 | 0.0006 | 0.0158 | 653 | 554 | 509 |
| `D95_day_linear_slp` | Flow Timing | linear_slp | 0.6768 | 0.7614 | 0.2042 | 2.7095 | 653 | 554 | 509 |
| `Flow_Reversals_annual_spearman_rho` | Pulse Metrics | spearman_rho | 0.6769 | 0.8309 | 0.1286 | 0.9938 | 645 | 544 | 513 |
| `log_a_seasonality_amplitude_first_half` | Recession | scalar | 0.6772 | 0.8605 | 0.9288 | 75.8573 | 1069 | 1234 | 165 |
| `n_low_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.6781 | 0.7674 | 0.0293 | 0.6891 | 645 | 544 | 513 |
| `Q5_spearman_rho` | Flow Percentiles | spearman_rho | 0.6783 | 0.8145 | 0.1226 | 1.2552 | 698 | 620 | 510 |
| `Q1_mk_rho` | Flow Percentiles | mk_rho | 0.6803 | 0.8123 | 0.0897 | 0.7992 | 715 | 641 | 510 |
| `BFI_Eckhardt_linear_slp` | Baseflow | linear_slp | 0.6819 | 0.7932 | 0.0005 | 0.0148 | 653 | 554 | 509 |
| `log_a_events_linear_slp` | Recession | linear_slp | 0.6825 | 0.7346 | 0.0298 | 1.8091 | 1068 | 1232 | 164 |
| `Q10_spearman_rho` | Flow Percentiles | spearman_rho | 0.6842 | 0.8159 | 0.1188 | 1.0311 | 687 | 605 | 514 |
| `Q50_spearman_rho` | Flow Percentiles | spearman_rho | 0.6846 | 0.8259 | 0.1023 | 0.8183 | 651 | 551 | 512 |
| `BFI_Eckhardt_senn_slp` | Baseflow | senn_slp | 0.6857 | 0.7902 | 0.0005 | 0.0110 | 653 | 554 | 509 |
| `Flow_Reversals_annual_mk_rho` | Pulse Metrics | mk_rho | 0.6863 | 0.8333 | 0.0916 | 0.6879 | 645 | 544 | 513 |
| `Q30_spearman_rho` | Flow Percentiles | spearman_rho | 0.6894 | 0.8191 | 0.1087 | 0.8754 | 663 | 566 | 513 |
| `FDC90th_linear_slp` | FDC | linear_slp | 0.6901 | 0.6915 | 0.0630 | 3.1092 | 645 | 544 | 513 |
| `Q25_spearman_rho` | Flow Percentiles | spearman_rho | 0.6918 | 0.8189 | 0.1106 | 0.8751 | 669 | 576 | 515 |
| `Q5_mk_rho` | Flow Percentiles | mk_rho | 0.6924 | 0.8130 | 0.0864 | 0.8233 | 698 | 620 | 510 |
| `Q20_spearman_rho` | Flow Percentiles | spearman_rho | 0.6933 | 0.8200 | 0.1128 | 0.8192 | 671 | 578 | 513 |
| `Q40_spearman_rho` | Flow Percentiles | spearman_rho | 0.6940 | 0.8249 | 0.1043 | 0.8798 | 657 | 558 | 513 |
| `winter_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.6950 | 0.7194 | 0.0034 | 2.1723 | 1609 | 1564 | 417 |
| `Q10_mk_rho` | Flow Percentiles | mk_rho | 0.6961 | 0.8127 | 0.0837 | 0.6869 | 687 | 605 | 514 |
| `Flow_Reversals_fall_senn_slp` | Pulse Metrics | senn_slp | 0.6981 | 0.8159 | 0.0706 | 1.0217 | 645 | 544 | 513 |
| `TQmean_senn_slp` | Pulse Metrics | senn_slp | 0.6987 | 0.7818 | 0.0578 | 0.8088 | 645 | 544 | 513 |
| `D99_day_senn_slp` | Flow Timing | senn_slp | 0.6998 | 0.6913 | 0.1000 | 4.0679 | 653 | 554 | 509 |
| `Q50_mk_rho` | Flow Percentiles | mk_rho | 0.7004 | 0.8270 | 0.0716 | 0.5226 | 651 | 551 | 512 |
| `Q30_mk_rho` | Flow Percentiles | mk_rho | 0.7037 | 0.8178 | 0.0764 | 0.5776 | 663 | 566 | 513 |
| `n_high_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.7056 | 0.8175 | 0.0232 | 0.4075 | 645 | 544 | 513 |
| `Q20_mk_rho` | Flow Percentiles | mk_rho | 0.7060 | 0.8186 | 0.0796 | 0.5712 | 671 | 578 | 513 |
| `Q25_mk_rho` | Flow Percentiles | mk_rho | 0.7060 | 0.8185 | 0.0778 | 0.5877 | 669 | 576 | 515 |
| `Q40_mk_rho` | Flow Percentiles | mk_rho | 0.7099 | 0.8250 | 0.0732 | 0.5644 | 657 | 558 | 513 |
| `n_low_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.7120 | 0.7803 | 0.0462 | 1.8651 | 645 | 544 | 513 |
| `TQmean_linear_slp` | Pulse Metrics | linear_slp | 0.7126 | 0.7937 | 0.0578 | 0.7953 | 645 | 544 | 513 |
| `D99_day_linear_slp` | Flow Timing | linear_slp | 0.7131 | 0.7127 | 0.1216 | 4.0688 | 653 | 554 | 509 |
| `summer_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.7137 | 0.6487 | 0.0055 | 4.9540 | 1786 | 1762 | 420 |
| `Flow_Reversals_fall_linear_slp` | Pulse Metrics | linear_slp | 0.7140 | 0.8212 | 0.0705 | 1.0094 | 645 | 544 | 513 |
| `n_low_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.7157 | 0.7575 | 0.0324 | 1.4148 | 645 | 544 | 513 |
| `Flow_Reversals_summer_senn_slp` | Pulse Metrics | senn_slp | 0.7179 | 0.8132 | 0.0645 | 1.3000 | 645 | 544 | 513 |
| `log_a_pointcloud_linear_slp` | Recession | linear_slp | 0.7232 | 0.7439 | 0.0080 | 0.4060 | 1249 | 1358 | 109 |
| `Flow_Reversals_summer_linear_slp` | Pulse Metrics | linear_slp | 0.7252 | 0.8186 | 0.0662 | 1.2650 | 645 | 544 | 513 |
| `Flow_Reversals_annual_senn_slp` | Pulse Metrics | senn_slp | 0.7308 | 0.8390 | 0.2043 | 3.4220 | 645 | 544 | 513 |
| `Flow_Reversals_spring_senn_slp` | Pulse Metrics | senn_slp | 0.7323 | 0.8098 | 0.0564 | 0.8920 | 645 | 544 | 513 |
| `flashinessRB_linear_slp` | Flashiness | linear_slp | 0.7352 | 0.7849 | 0.0007 | 0.0318 | 653 | 554 | 509 |
| `Flow_Reversals_annual_linear_slp` | Pulse Metrics | linear_slp | 0.7386 | 0.8360 | 0.2082 | 3.3095 | 645 | 544 | 513 |
| `Flow_Reversals_spring_linear_slp` | Pulse Metrics | linear_slp | 0.7398 | 0.8160 | 0.0586 | 0.8214 | 645 | 544 | 513 |
| `negative_ann_spearman_rho` | Negative Flow | spearman_rho | 0.7399 | 0.9094 | 0.0610 | 0.3326 | 6634 | 6629 | 5 |
| `n_recession_events_senn_slp` | Recession | senn_slp | 0.7447 | 0.8144 | 0.0157 | 0.2229 | 0 | 0 | 0 |
| `flashinessRB_senn_slp` | Flashiness | senn_slp | 0.7485 | 0.7933 | 0.0007 | 0.0174 | 653 | 554 | 509 |
| `n_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.7500 | 0.7726 | 0.0271 | 0.5516 | 645 | 544 | 513 |
| `negative_ann_mk_rho` | Negative Flow | mk_rho | 0.7515 | 0.9148 | 0.0476 | 0.2671 | 6634 | 6629 | 5 |
| `FDCall_senn_slp` | FDC | senn_slp | 0.7522 | 0.7581 | 0.0067 | 0.3203 | 645 | 544 | 513 |
| `b_pointcloud_linear_slp` | Recession | linear_slp | 0.7568 | 0.7571 | 0.0056 | 0.2029 | 1249 | 1358 | 109 |
| `FDCmid_senn_slp` | FDC | senn_slp | 0.7600 | 0.7533 | 0.0060 | 0.4185 | 645 | 544 | 513 |
| `FDCmid_linear_slp` | FDC | linear_slp | 0.7621 | 0.7491 | 0.0092 | 0.5769 | 645 | 544 | 513 |
| `n_recession_events_linear_slp` | Recession | linear_slp | 0.7631 | 0.8529 | 0.0184 | 0.2544 | 0 | 0 | 0 |
| `n_high_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.7648 | 0.8277 | 0.0311 | 0.6639 | 645 | 544 | 513 |
| `summer_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.7694 | 0.5734 | 0.0088 | 4.8511 | 1786 | 1762 | 420 |
| `elasticity_annual_mean` | Elasticity | mean | 0.7739 | 0.8236 | 1.9697 | 199.3512 | 1165 | 1165 | 0 |
| `FDCall_linear_slp` | FDC | linear_slp | 0.7786 | 0.7504 | 0.0081 | 0.2994 | 645 | 544 | 513 |
| `log_a_seasonality_minimum_last_half` | Recession | scalar | 0.7790 | 0.9241 | 9.9258 | 360.9091 | 1072 | 1234 | 162 |
| `b_pointcloud_senn_slp` | Recession | senn_slp | 0.7836 | 0.7940 | 0.0050 | 0.4522 | 1249 | 1358 | 109 |
| `elasticity_static` | Elasticity | scalar | 0.8107 | 0.9001 | 0.2109 | 2.9901 | 1165 | 1165 | 0 |
| `dur_low_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.8114 | 0.7434 | 0.0534 | 1.0872 | 1186 | 1058 | 550 |
| `log_a_pointcloud_senn_slp` | Recession | senn_slp | 0.8198 | 0.7894 | 0.0070 | 0.6651 | 1249 | 1358 | 109 |
| `negative_ann_senn_slp` | Negative Flow | senn_slp | 0.8202 | 0.8660 | 0.0000 | 0.1786 | 645 | 544 | 513 |
| `log_a_seasonality_minimum_all` | Recession | scalar | 0.8254 | 0.9451 | 7.1119 | 363.3182 | 1068 | 1232 | 164 |
| `dur_low_pulses_all_median` | Pulse Metrics | median | 0.8269 | 0.9327 | 1.4599 | 116.3139 | 383 | 375 | 40 |
| `elasticity_annual_median` | Elasticity | median | 0.8485 | 0.9194 | 0.1624 | 3.1040 | 1165 | 1165 | 0 |
| `runoff_ratio_high_count` | Runoff Ratios | scalar | 0.8562 | 0.9662 | 0.0363 | 12.0000 | 1165 | 1165 | 0 |
| `FDC90th_median` | FDC | median | 0.8787 | 0.9786 | 0.1963 | 40.4516 | 0 | 0 | 0 |
| `elasticity_rolling_median` | Elasticity | median | 0.8817 | 0.9363 | 0.1590 | 2.4302 | 1165 | 1165 | 0 |
| `n_low_pulses_all_median` | Pulse Metrics | median | 0.8862 | 0.9362 | 0.3803 | 11.5000 | 0 | 0 | 0 |
| `dur_low_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.8882 | 0.7200 | 0.0292 | 0.7195 | 1186 | 1058 | 550 |
| `log_a_seasonality_amplitude_last_half` | Recession | scalar | 0.9035 | 0.9616 | 0.5728 | 55.2083 | 1072 | 1234 | 162 |
| `dur_low_pulses_all_mean` | Pulse Metrics | mean | 0.9107 | 0.9642 | 1.4675 | 85.8383 | 383 | 375 | 40 |
| `elasticity_rolling_mean` | Elasticity | mean | 0.9164 | 0.9528 | 0.1545 | 1.9362 | 1165 | 1165 | 0 |
| `log_a_seasonality_amplitude_all` | Recession | scalar | 0.9200 | 0.9792 | 0.3869 | 39.6246 | 1068 | 1232 | 164 |
| `FDCmid_median` | FDC | median | 0.9384 | 0.9908 | 0.0776 | 10.1806 | 0 | 0 | 0 |
| `dur_low_pulses_year_median` | Pulse Metrics | median | 0.9406 | 0.9794 | 0.3638 | 18.5000 | 85 | 103 | 18 |
| `D1_day_median` | Flow Timing | median | 0.9427 | 0.9877 | 1.0095 | 133.5000 | 0 | 0 | 0 |
| `n_low_pulses_all_mean` | Pulse Metrics | mean | 0.9474 | 0.9707 | 0.2757 | 10.7591 | 0 | 0 | 0 |
| `dur_high_pulses_all_median` | Pulse Metrics | median | 0.9483 | 0.9916 | 0.6224 | 21.6250 | 0 | 0 | 0 |
| `concavity_median` | Recession | median | 0.9551 | 0.9850 | 0.1397 | 5.0101 | 1068 | 1232 | 164 |
| `qp_slope_sd_mean` | Q-P Seasonality | mean | 0.9582 | 0.9670 | 0.6407 | 66.8511 | 1165 | 1165 | 0 |
| `D5_day_median` | Flow Timing | median | 0.9601 | 0.9864 | 2.0479 | 116.0000 | 0 | 0 | 0 |
| `n_low_pulses_year_median` | Pulse Metrics | median | 0.9603 | 0.9749 | 0.2540 | 14.5000 | 0 | 0 | 0 |
| `dur_high_pulses_all_mean` | Pulse Metrics | mean | 0.9619 | 0.9961 | 0.6314 | 35.3472 | 0 | 0 | 0 |
| `FDC90th_mean` | FDC | mean | 0.9640 | 0.9786 | 0.4511 | 14.5611 | 0 | 0 | 0 |
| `concavity_mean` | Recession | mean | 0.9650 | 0.9902 | 0.1849 | 5.9918 | 1068 | 1232 | 164 |
| `b_pointcloud_median` | Recession | median | 0.9653 | 0.9820 | 0.0308 | 2.0134 | 1249 | 1358 | 109 |
| `Dmax_median` | Flow Timing | median | 0.9656 | 0.9805 | 4.3053 | 163.5000 | 0 | 0 | 0 |
| `D99_day_median` | Flow Timing | median | 0.9676 | 0.9890 | 1.3091 | 101.5000 | 0 | 0 | 0 |
| `dur_high_pulses_year_median` | Pulse Metrics | median | 0.9686 | 0.9960 | 0.2135 | 18.0000 | 0 | 0 | 0 |
| `FDCall_median` | FDC | median | 0.9688 | 0.9952 | 0.0767 | 6.6265 | 0 | 0 | 0 |
| `D10_day_median` | Flow Timing | median | 0.9700 | 0.9854 | 2.7347 | 93.0000 | 0 | 0 | 0 |
| `b_events_median` | Recession | median | 0.9714 | 0.9924 | 0.0888 | 3.9277 | 1068 | 1232 | 164 |
| `Flow_Reversals_fall_median` | Pulse Metrics | median | 0.9737 | 0.9859 | 0.7398 | 18.0000 | 0 | 0 | 0 |
| `b_events_mean` | Recession | mean | 0.9754 | 0.9941 | 0.1189 | 4.7767 | 1068 | 1232 | 164 |
| `b_pointcloud_mean` | Recession | mean | 0.9761 | 0.9865 | 0.0319 | 1.1528 | 1249 | 1358 | 109 |
| `D25_to_D75_median` | Flow Timing | median | 0.9762 | 0.9869 | 2.8645 | 88.0000 | 0 | 0 | 0 |
| `Flow_Reversals_winter_median` | Pulse Metrics | median | 0.9771 | 0.9866 | 0.7050 | 16.5000 | 0 | 0 | 0 |
| `log_a_events_mean` | Recession | mean | 0.9772 | 0.9902 | 0.1934 | 9.2042 | 1068 | 1232 | 164 |
| `log_a_events_median` | Recession | median | 0.9772 | 0.9944 | 0.1287 | 9.9365 | 1068 | 1232 | 164 |
| `qp_bimodality_median` | Q-P Seasonality | median | 0.9781 | 0.9888 | 0.0118 | 0.1988 | 1165 | 1165 | 0 |
| `D1_day_mean` | Flow Timing | mean | 0.9782 | 0.9941 | 0.9714 | 45.5784 | 0 | 0 | 0 |
| `FDCmid_mean` | FDC | mean | 0.9790 | 0.9958 | 0.0749 | 4.7324 | 0 | 0 | 0 |
| `Dmax_mean` | Flow Timing | mean | 0.9790 | 0.9860 | 3.9167 | 35.0664 | 0 | 0 | 0 |
| `n_recession_events_median` | Recession | median | 0.9791 | 0.9915 | 0.2260 | 2.5000 | 0 | 0 | 0 |
| `dur_low_pulses_year_mean` | Pulse Metrics | mean | 0.9794 | 0.9861 | 0.3916 | 8.4316 | 85 | 103 | 18 |
| `Flow_Reversals_summer_median` | Pulse Metrics | median | 0.9797 | 0.9896 | 0.6744 | 15.5000 | 0 | 0 | 0 |
| `D20_day_median` | Flow Timing | median | 0.9800 | 0.9824 | 2.9477 | 79.0000 | 0 | 0 | 0 |
| `n_low_pulses_year_mean` | Pulse Metrics | mean | 0.9803 | 0.9887 | 0.2263 | 4.0621 | 0 | 0 | 0 |
| `log_a_pointcloud_median` | Recession | median | 0.9805 | 0.9944 | 0.0448 | 2.4972 | 1249 | 1358 | 109 |
| `Flow_Reversals_fall_mean` | Pulse Metrics | mean | 0.9811 | 0.9906 | 0.6988 | 9.4646 | 0 | 0 | 0 |
| `qp_slope_sd_median` | Q-P Seasonality | median | 0.9812 | 0.9893 | 0.1341 | 24.5498 | 1165 | 1165 | 0 |
| `Flow_Reversals_spring_median` | Pulse Metrics | median | 0.9814 | 0.9895 | 0.5566 | 16.0000 | 0 | 0 | 0 |
| `D30_day_median` | Flow Timing | median | 0.9814 | 0.9852 | 2.6963 | 145.0000 | 0 | 0 | 0 |
| `D40_day_median` | Flow Timing | median | 0.9816 | 0.9868 | 2.5557 | 97.0000 | 0 | 0 | 0 |
| `D95_day_median` | Flow Timing | median | 0.9840 | 0.9916 | 2.0646 | 52.0000 | 0 | 0 | 0 |
| `TQmean_median` | Pulse Metrics | median | 0.9840 | 0.9915 | 0.5075 | 24.8245 | 0 | 0 | 0 |
| `D50_day_median` | Flow Timing | median | 0.9846 | 0.9907 | 2.1951 | 109.5000 | 0 | 0 | 0 |
| `D5_day_mean` | Flow Timing | mean | 0.9847 | 0.9946 | 1.6237 | 38.9594 | 0 | 0 | 0 |
| `Flow_Reversals_annual_median` | Pulse Metrics | median | 0.9847 | 0.9914 | 2.0960 | 57.0000 | 0 | 0 | 0 |
| `D90_day_median` | Flow Timing | median | 0.9852 | 0.9915 | 2.1610 | 66.0000 | 0 | 0 | 0 |
| `Flow_Reversals_winter_mean` | Pulse Metrics | mean | 0.9854 | 0.9915 | 0.6166 | 9.5007 | 0 | 0 | 0 |
| `D80_day_median` | Flow Timing | median | 0.9857 | 0.9910 | 2.0862 | 101.0000 | 0 | 0 | 0 |
| `Flow_Reversals_summer_mean` | Pulse Metrics | mean | 0.9859 | 0.9933 | 0.6218 | 8.2964 | 0 | 0 | 0 |
| `FDCall_mean` | FDC | mean | 0.9861 | 0.9962 | 0.0703 | 2.9295 | 0 | 0 | 0 |
| `D60_day_median` | Flow Timing | median | 0.9863 | 0.9921 | 2.0469 | 86.5000 | 0 | 0 | 0 |
| `summer_runoff_ratio_median` | Runoff Ratios | median | 0.9866 | 0.9945 | 0.0577 | 65.5529 | 1165 | 1165 | 0 |
| `D70_day_median` | Flow Timing | median | 0.9866 | 0.9920 | 1.9977 | 94.0000 | 0 | 0 | 0 |
| `qp_bimodality_mean` | Q-P Seasonality | mean | 0.9872 | 0.9933 | 0.0092 | 0.0752 | 1165 | 1165 | 0 |
| `Flow_Reversals_spring_mean` | Pulse Metrics | mean | 0.9876 | 0.9933 | 0.5113 | 8.0909 | 0 | 0 | 0 |
| `D25_to_D75_mean` | Flow Timing | mean | 0.9876 | 0.9927 | 2.2960 | 34.0380 | 0 | 0 | 0 |
| `Flow_Reversals_annual_mean` | Pulse Metrics | mean | 0.9876 | 0.9934 | 2.0162 | 35.7075 | 0 | 0 | 0 |
| `D10_day_mean` | Flow Timing | mean | 0.9881 | 0.9938 | 1.9950 | 27.5154 | 0 | 0 | 0 |
| `log_a_pointcloud_mean` | Recession | mean | 0.9885 | 0.9954 | 0.0451 | 1.6429 | 1249 | 1358 | 109 |
| `D99_day_mean` | Flow Timing | mean | 0.9887 | 0.9948 | 1.0509 | 30.2521 | 0 | 0 | 0 |
| `n_high_pulses_all_median` | Pulse Metrics | median | 0.9891 | 0.9918 | 0.2854 | 7.5000 | 0 | 0 | 0 |
| `avg_storage_median` | Storage | median | 0.9893 | 0.9827 | 22.0267 | 18500.4942 | 1165 | 1165 | 0 |

## NA Mismatch Analysis

Columns where the number of NAs differs by >50 gages.

| Column | Category | Baseline NAs | Experiment NAs | Mismatch | R2 |
|--------|----------|--------------|----------------|----------|----|
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | 4364 | 4141 | 755 | -0.0700 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | 4364 | 4141 | 755 | -1.0606 |
| `dur_low_pulses_all_senn_slp` | Pulse Metrics | 4364 | 4141 | 755 | -1.4812 |
| `dur_low_pulses_all_mk_rho` | Pulse Metrics | 4364 | 4141 | 755 | 0.5726 |
| `dur_low_pulses_all_spearman_rho` | Pulse Metrics | 4364 | 4141 | 755 | 0.5589 |
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | 4364 | 4141 | 755 | -0.0827 |
| `dur_high_pulses_all_spearman_rho` | Pulse Metrics | 1487 | 1352 | 673 | 0.5756 |
| `dur_high_pulses_all_spearman_pval` | Pulse Metrics | 1487 | 1352 | 673 | -0.0809 |
| `dur_high_pulses_all_senn_slp` | Pulse Metrics | 1487 | 1352 | 673 | 0.5460 |
| `dur_high_pulses_all_mk_rho` | Pulse Metrics | 1487 | 1352 | 673 | 0.5811 |
| `dur_high_pulses_all_mk_pval` | Pulse Metrics | 1487 | 1352 | 673 | -0.0740 |
| `dur_high_pulses_all_linear_slp` | Pulse Metrics | 1487 | 1352 | 673 | 0.6392 |
| `dur_low_pulses_year_spearman_pval` | Pulse Metrics | 1186 | 1058 | 550 | -0.0752 |
| `dur_low_pulses_year_spearman_rho` | Pulse Metrics | 1186 | 1058 | 550 | 0.5455 |
| `dur_low_pulses_year_linear_slp` | Pulse Metrics | 1186 | 1058 | 550 | 0.8114 |
| `dur_low_pulses_year_mk_pval` | Pulse Metrics | 1186 | 1058 | 550 | -0.0638 |
| `dur_low_pulses_year_senn_slp` | Pulse Metrics | 1186 | 1058 | 550 | 0.8882 |
| `dur_low_pulses_year_mk_rho` | Pulse Metrics | 1186 | 1058 | 550 | 0.5612 |
| `Flow_Reversals_winter_spearman_rho` | Pulse Metrics | 647 | 549 | 516 | 0.6517 |
| `Flow_Reversals_winter_mk_pval` | Pulse Metrics | 647 | 549 | 516 | 0.0390 |
| `Flow_Reversals_winter_spearman_pval` | Pulse Metrics | 647 | 549 | 516 | 0.0516 |
| `Flow_Reversals_winter_mk_rho` | Pulse Metrics | 647 | 549 | 516 | 0.6570 |
| `Q25_mk_rho` | Flow Percentiles | 669 | 576 | 515 | 0.7060 |
| `Q25_mk_pval` | Flow Percentiles | 669 | 576 | 515 | 0.0703 |
| `Q25_spearman_rho` | Flow Percentiles | 669 | 576 | 515 | 0.6918 |
| `Q25_spearman_pval` | Flow Percentiles | 669 | 576 | 515 | 0.0746 |
| `n_low_pulses_year_spearman_rho` | Pulse Metrics | 688 | 606 | 514 | 0.5728 |
| `n_low_pulses_year_spearman_pval` | Pulse Metrics | 688 | 606 | 514 | -0.0706 |
| `Q10_mk_pval` | Flow Percentiles | 687 | 605 | 514 | 0.0234 |
| `Q10_mk_rho` | Flow Percentiles | 687 | 605 | 514 | 0.6961 |

## Summary

| Agreement Tier | Threshold | Count | % |
|----------------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 35 | 5.9% |
| Good | 0.99 <= R2 < 0.999 | 41 | 6.9% |
| Poor | 0.95 <= R2 < 0.99 | 58 | 9.8% |
| Low | 0.90 <= R2 < 0.95 | 9 | 1.5% |
| Very Low | 0.50 <= R2 < 0.90 | 228 | 38.4% |
| Extremely Low | R2 < 0.50 | 223 | 37.5% |
| **Total compared** | | **594** | **100%** |

Gages dropped by experiment filter: **635**

---
*Generated by `docs/benchmarks/compare_experiment_vs_julia.py startIn1993`*

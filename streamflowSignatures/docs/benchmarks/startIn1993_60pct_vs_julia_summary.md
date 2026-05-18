# Experiment 'startIn1993_60pct' vs Julia Baseline: Comparison Report

**Generated**: 2026-04-24 14:06:01

## Experiment Description

Water years >= 1993 AND at least 60% of possible years (WY1993 to gage max) must have qualifying data.

## Input Files

| Dataset | File | Gages | Columns |
|---------|------|-------|---------|
| Julia Baseline | `julia_signatures.csv` | 7,313 | 656 |
| Experiment (startIn1993_60pct) | `startIn1993_60pct_signatures.csv` | 6,579 | 656 |

**Common gages**: 6,579
**Dropped gages** (in baseline, not in experiment): 734
**Added gages** (in experiment, not in baseline): 0

### Gage Differences

```

  Dropped (baseline only): 734 gages
    USGS: 584
    Canada: 150
    non-reference: 582
    reference: 152
    Examples: 01021200, 01060000, 01064000, 01064118, 01064300
```

### Years Per Gage (Common Gages)

```
Years per gage (common gages):
    Baseline:   mean=39.7, median=44.0
    Experiment: mean=30.2, median=32.0
    Mean diff:  -9.4 years
```

## High-Level Alignment Summary

### Distribution Statistics

| Metric | Value |
|--------|-------|
| Columns compared | 619 |
| Mean R2 (identity) | 0.433049 |
| Median R2 | 0.608309 |
| SD of R2 | 0.913739 |
| Min R2 | -18.6911 |

### Agreement Tiers

| Tier | Threshold | Count | % |
|------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 34 | 5.5% |
| Good | 0.99 <= R2 < 0.999 | 42 | 6.8% |
| Poor | 0.95 <= R2 < 0.99 | 65 | 10.5% |
| Low | 0.90 <= R2 < 0.95 | 9 | 1.5% |
| Very Low | 0.50 <= R2 < 0.90 | 236 | 38.1% |
| Extremely Low | R2 < 0.50 | 233 | 37.6% |
| **Total** | | **619** | **100%** |

## Agreement by Signature Category

| Category | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| Baseflow | 32 | 0 | 4 | 4 | 0 | 16 | 8 | 0.567489 | -0.0530 |
| Elasticity | 19 | 1 | 0 | 0 | 1 | 5 | 12 | -0.676324 | -18.6911 |
| FDC | 24 | 0 | 0 | 4 | 1 | 13 | 6 | 0.564987 | -0.0938 |
| Flashiness | 8 | 0 | 2 | 0 | 0 | 4 | 2 | 0.578405 | -0.0318 |
| Flow Percentiles | 128 | 21 | 11 | 0 | 0 | 32 | 64 | 0.319715 | -1.2369 |
| Flow Timing | 120 | 0 | 9 | 20 | 1 | 47 | 43 | 0.513275 | -0.2241 |
| Flow Volumes | 40 | 7 | 3 | 0 | 0 | 10 | 20 | 0.331756 | -1.5985 |
| Negative Flow | 8 | 0 | 2 | 0 | 0 | 4 | 2 | 0.612747 | -0.0104 |
| Pulse Metrics | 112 | 0 | 5 | 17 | 4 | 55 | 31 | 0.526517 | -1.8187 |
| Q-P Seasonality | 16 | 0 | 0 | 4 | 0 | 5 | 7 | 0.427829 | -0.2851 |
| Recession | 63 | 0 | 1 | 14 | 2 | 32 | 14 | 0.617799 | -0.0237 |
| Runoff Ratios | 41 | 5 | 4 | 1 | 0 | 12 | 19 | 0.381886 | -2.2003 |
| Storage | 8 | 0 | 1 | 1 | 0 | 1 | 5 | 0.434787 | -0.2640 |

## Agreement by Statistic Type

| Stat Type | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|-----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| mean | 76 | 19 | 26 | 27 | 3 | 1 | 0 | 0.985121 | 0.7739 |
| median | 76 | 14 | 16 | 37 | 4 | 5 | 0 | 0.976275 | 0.8237 |
| senn_slp | 76 | 0 | 0 | 0 | 0 | 47 | 29 | 0.333459 | -1.8187 |
| linear_slp | 76 | 0 | 0 | 0 | 0 | 49 | 27 | 0.304128 | -2.2003 |
| spearman_rho | 76 | 0 | 0 | 0 | 0 | 63 | 13 | 0.578512 | 0.2041 |
| spearman_pval | 76 | 0 | 0 | 0 | 0 | 0 | 76 | -0.055726 | -0.7371 |
| mk_rho | 76 | 0 | 0 | 0 | 0 | 65 | 11 | 0.591720 | 0.2375 |
| mk_pval | 76 | 0 | 0 | 0 | 0 | 0 | 76 | -0.049778 | -0.7465 |
| scalar | 11 | 1 | 0 | 1 | 2 | 6 | 1 | -0.944057 | -18.6911 |

## Columns with R2 < 0.99 (543 columns)

| Column | Category | Stat | R2 | Spearman | MAD | Max Diff | Baseline NAs | Experiment NAs | NA Mismatch |
|--------|----------|------|----|----------|-----|----------|--------------|----------------|-------------|
| `elasticity_years_total` | Elasticity | scalar | -18.6911 | 0.7494 | 10.3147 | 12.0000 | 1085 | 1085 | 0 |
| `spring_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | -2.2003 | 0.7478 | 0.0047 | 4.0210 | 1514 | 1469 | 419 |
| `dur_low_pulses_all_senn_slp` | Pulse Metrics | senn_slp | -1.8187 | 0.7588 | 0.0710 | 15.4390 | 4299 | 4080 | 751 |
| `Qfal_linear_slp` | Flow Volumes | linear_slp | -1.5985 | 0.6894 | 1.4479 | 3175.9690 | 607 | 514 | 505 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | linear_slp | -1.3108 | 0.7555 | 0.0928 | 16.7031 | 4299 | 4080 | 751 |
| `Q10_senn_slp` | Flow Percentiles | senn_slp | -1.2369 | 0.7765 | 0.0118 | 32.7018 | 607 | 514 | 505 |
| `Q5_senn_slp` | Flow Percentiles | senn_slp | -1.2301 | 0.7749 | 0.0115 | 28.2927 | 607 | 514 | 505 |
| `Q1_senn_slp` | Flow Percentiles | senn_slp | -1.2166 | 0.7777 | 0.0104 | 27.0480 | 607 | 514 | 505 |
| `Q10_linear_slp` | Flow Percentiles | linear_slp | -1.1989 | 0.7903 | 0.0103 | 25.4106 | 607 | 514 | 505 |
| `Q20_linear_slp` | Flow Percentiles | linear_slp | -1.1582 | 0.8035 | 0.0107 | 27.6171 | 607 | 514 | 505 |
| `Q25_senn_slp` | Flow Percentiles | senn_slp | -1.0546 | 0.7901 | 0.0140 | 43.0240 | 607 | 514 | 505 |
| `Q20_senn_slp` | Flow Percentiles | senn_slp | -1.0091 | 0.7941 | 0.0135 | 42.1578 | 607 | 514 | 505 |
| `Q30_senn_slp` | Flow Percentiles | senn_slp | -0.8716 | 0.7731 | 0.0152 | 45.2000 | 607 | 514 | 505 |
| `Q25_linear_slp` | Flow Percentiles | linear_slp | -0.8692 | 0.8067 | 0.0111 | 29.3314 | 607 | 514 | 505 |
| `Q30_linear_slp` | Flow Percentiles | linear_slp | -0.7940 | 0.8022 | 0.0119 | 31.5025 | 607 | 514 | 505 |
| `Q5_linear_slp` | Flow Percentiles | linear_slp | -0.7796 | 0.7911 | 0.0100 | 24.7992 | 607 | 514 | 505 |
| `elasticity_rolling_mk_pval` | Elasticity | mk_pval | -0.7465 | 0.1869 | 0.2779 | 1.0000 | 1085 | 1085 | 0 |
| `elasticity_rolling_spearman_pval` | Elasticity | spearman_pval | -0.7371 | 0.1965 | 0.2647 | 0.9967 | 1085 | 1085 | 0 |
| `Q1_linear_slp` | Flow Percentiles | linear_slp | -0.7139 | 0.7838 | 0.0087 | 19.2563 | 607 | 514 | 505 |
| `Qfal_senn_slp` | Flow Volumes | senn_slp | -0.6873 | 0.6757 | 1.5263 | 4084.4444 | 607 | 514 | 505 |
| `Q50_linear_slp` | Flow Percentiles | linear_slp | -0.5057 | 0.8149 | 0.0136 | 35.8286 | 607 | 514 | 505 |
| `Q40_linear_slp` | Flow Percentiles | linear_slp | -0.5009 | 0.8073 | 0.0132 | 34.0502 | 607 | 514 | 505 |
| `Q40_senn_slp` | Flow Percentiles | senn_slp | -0.4442 | 0.7815 | 0.0165 | 48.1071 | 607 | 514 | 505 |
| `spring_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | -0.4340 | 0.7319 | 0.0035 | 2.9330 | 1514 | 1469 | 419 |
| `elasticity_annual_spearman_pval` | Elasticity | spearman_pval | -0.3734 | 0.3174 | 0.2658 | 0.9896 | 1085 | 1085 | 0 |
| `elasticity_annual_mk_pval` | Elasticity | mk_pval | -0.3510 | 0.3266 | 0.2664 | 0.9650 | 1085 | 1085 | 0 |
| `summer_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.3263 | 0.3992 | 0.2599 | 0.9958 | 1688 | 1665 | 419 |
| `Qann_linear_slp` | Flow Volumes | linear_slp | -0.3231 | 0.7918 | 5.2480 | 11771.1424 | 607 | 514 | 505 |
| `summer_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.3092 | 0.4035 | 0.2627 | 0.9988 | 1688 | 1665 | 419 |
| `qp_bimodality_spearman_pval` | Q-P Seasonality | spearman_pval | -0.2851 | 0.3790 | 0.2537 | 0.9854 | 1517 | 1479 | 422 |
| `Q50_senn_slp` | Flow Percentiles | senn_slp | -0.2773 | 0.7917 | 0.0155 | 46.9455 | 607 | 514 | 505 |
| `Q60_linear_slp` | Flow Percentiles | linear_slp | -0.2757 | 0.8012 | 0.0141 | 33.0408 | 607 | 514 | 505 |
| `qp_bimodality_mk_pval` | Q-P Seasonality | mk_pval | -0.2648 | 0.3869 | 0.2546 | 0.9893 | 1517 | 1479 | 422 |
| `avg_storage_spearman_pval` | Storage | spearman_pval | -0.2640 | 0.3887 | 0.2543 | 0.9481 | 1529 | 1486 | 421 |
| `qp_slope_sd_spearman_pval` | Q-P Seasonality | spearman_pval | -0.2554 | 0.4277 | 0.2491 | 0.9614 | 1510 | 1471 | 421 |
| `qp_slope_sd_mk_pval` | Q-P Seasonality | mk_pval | -0.2512 | 0.4291 | 0.2521 | 0.9946 | 1510 | 1471 | 421 |
| `avg_storage_mk_pval` | Storage | mk_pval | -0.2471 | 0.3920 | 0.2557 | 0.9798 | 1529 | 1486 | 421 |
| `spring_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.2335 | 0.4519 | 0.2482 | 0.9842 | 1514 | 1469 | 419 |
| `Qwin_senn_slp` | Flow Volumes | senn_slp | -0.2285 | 0.7719 | 1.4578 | 3283.7500 | 607 | 514 | 505 |
| `D99_day_mk_pval` | Flow Timing | mk_pval | -0.2241 | 0.4594 | 0.2310 | 0.9976 | 615 | 524 | 501 |
| `spring_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.2208 | 0.4545 | 0.2498 | 0.9931 | 1514 | 1469 | 419 |
| `D99_day_spearman_pval` | Flow Timing | spearman_pval | -0.2169 | 0.4621 | 0.2266 | 0.9828 | 615 | 524 | 501 |
| `Qann_senn_slp` | Flow Volumes | senn_slp | -0.2102 | 0.7619 | 5.9053 | 14990.1222 | 607 | 514 | 505 |
| `fall_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.2041 | 0.4665 | 0.2570 | 0.9994 | 1513 | 1469 | 420 |
| `fall_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.1928 | 0.4704 | 0.2535 | 0.9962 | 1513 | 1469 | 420 |
| `D20_day_spearman_pval` | Flow Timing | spearman_pval | -0.1917 | 0.4240 | 0.2154 | 0.9800 | 615 | 524 | 501 |
| `D20_day_mk_pval` | Flow Timing | mk_pval | -0.1840 | 0.4308 | 0.2167 | 0.9628 | 615 | 524 | 501 |
| `winter_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.1757 | 0.5077 | 0.2467 | 0.9945 | 1511 | 1467 | 416 |
| `Q95_Q10_spearman_pval` | Flow Percentiles | spearman_pval | -0.1510 | 0.4756 | 0.2224 | 0.9854 | 607 | 514 | 505 |
| `winter_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.1454 | 0.5189 | 0.2465 | 0.9968 | 1511 | 1467 | 416 |
| `Q95_spearman_pval` | Flow Percentiles | spearman_pval | -0.1423 | 0.4815 | 0.2209 | 0.9928 | 607 | 514 | 505 |
| `D10_day_spearman_pval` | Flow Timing | spearman_pval | -0.1365 | 0.4452 | 0.2137 | 0.9928 | 615 | 524 | 501 |
| `Q90_spearman_pval` | Flow Percentiles | spearman_pval | -0.1310 | 0.4847 | 0.2214 | 0.9856 | 607 | 514 | 505 |
| `Dmax_spearman_pval` | Flow Timing | spearman_pval | -0.1296 | 0.4263 | 0.2160 | 0.9788 | 615 | 524 | 501 |
| `Q90_mk_pval` | Flow Percentiles | mk_pval | -0.1292 | 0.4808 | 0.2250 | 0.9939 | 607 | 514 | 505 |
| `Dmax_mk_pval` | Flow Timing | mk_pval | -0.1266 | 0.4271 | 0.2180 | 0.9727 | 615 | 524 | 501 |
| `Q95_mk_pval` | Flow Percentiles | mk_pval | -0.1259 | 0.4835 | 0.2231 | 0.9835 | 607 | 514 | 505 |
| `Q95_Q10_mk_pval` | Flow Percentiles | mk_pval | -0.1257 | 0.4784 | 0.2245 | 0.9835 | 607 | 514 | 505 |
| `Q99_spearman_pval` | Flow Percentiles | spearman_pval | -0.1250 | 0.4802 | 0.2200 | 0.9889 | 607 | 514 | 505 |
| `D10_day_mk_pval` | Flow Timing | mk_pval | -0.1159 | 0.4566 | 0.2150 | 0.9576 | 615 | 524 | 501 |
| `D5_day_spearman_pval` | Flow Timing | spearman_pval | -0.1122 | 0.4555 | 0.2101 | 0.9631 | 615 | 524 | 501 |
| `D5_day_mk_pval` | Flow Timing | mk_pval | -0.1108 | 0.4552 | 0.2128 | 0.9771 | 615 | 524 | 501 |
| `D30_day_mk_pval` | Flow Timing | mk_pval | -0.1079 | 0.4642 | 0.2073 | 0.9521 | 615 | 524 | 501 |
| `Q99_mk_pval` | Flow Percentiles | mk_pval | -0.1069 | 0.4850 | 0.2213 | 0.9902 | 607 | 514 | 505 |
| `D50_day_spearman_pval` | Flow Timing | spearman_pval | -0.1037 | 0.4687 | 0.2081 | 0.9357 | 615 | 524 | 501 |
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.1004 | 0.4969 | 0.2152 | 0.9796 | 4299 | 4080 | 751 |
| `D50_day_mk_pval` | Flow Timing | mk_pval | -0.0973 | 0.4720 | 0.2089 | 0.9652 | 615 | 524 | 501 |
| `annual_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.0951 | 0.5299 | 0.2392 | 0.9970 | 1509 | 1466 | 417 |
| `FDC90th_spearman_pval` | FDC | spearman_pval | -0.0938 | 0.5198 | 0.2207 | 0.9944 | 652 | 578 | 506 |
| `dur_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.0930 | 0.4772 | 0.2191 | 0.9789 | 1443 | 1316 | 665 |
| `Qsum_spearman_pval` | Flow Volumes | spearman_pval | -0.0917 | 0.5171 | 0.2125 | 0.9994 | 607 | 514 | 505 |
| `Qsum_mk_pval` | Flow Volumes | mk_pval | -0.0916 | 0.5151 | 0.2147 | 0.9965 | 608 | 515 | 505 |
| `FDC90th_mk_pval` | FDC | mk_pval | -0.0881 | 0.5194 | 0.2230 | 0.9972 | 652 | 578 | 506 |
| `dur_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.0878 | 0.5010 | 0.2216 | 0.9880 | 1146 | 1024 | 542 |
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0876 | 0.5040 | 0.2176 | 0.9792 | 4299 | 4080 | 751 |
| `D60_day_mk_pval` | Flow Timing | mk_pval | -0.0872 | 0.4794 | 0.2100 | 0.9975 | 615 | 524 | 501 |
| `dur_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0867 | 0.4803 | 0.2208 | 0.9968 | 1443 | 1316 | 665 |
| `annual_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.0864 | 0.5347 | 0.2345 | 0.9957 | 1509 | 1466 | 417 |
| `D60_day_spearman_pval` | Flow Timing | spearman_pval | -0.0859 | 0.4794 | 0.2079 | 0.9651 | 615 | 524 | 501 |
| `D30_day_spearman_pval` | Flow Timing | spearman_pval | -0.0853 | 0.4723 | 0.2040 | 0.9597 | 615 | 524 | 501 |
| `n_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.0841 | 0.5159 | 0.2185 | 0.9899 | 650 | 576 | 506 |
| `n_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.0818 | 0.5160 | 0.2217 | 0.9920 | 650 | 576 | 506 |
| `D40_day_spearman_pval` | Flow Timing | spearman_pval | -0.0818 | 0.4808 | 0.2055 | 0.9458 | 615 | 524 | 501 |
| `D40_day_mk_pval` | Flow Timing | mk_pval | -0.0810 | 0.4772 | 0.2085 | 0.9694 | 615 | 524 | 501 |
| `Qwin_linear_slp` | Flow Volumes | linear_slp | -0.0794 | 0.8095 | 1.4146 | 2572.4661 | 607 | 514 | 505 |
| `dur_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.0773 | 0.5064 | 0.2244 | 0.9975 | 1146 | 1024 | 542 |
| `n_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0754 | 0.5067 | 0.2213 | 0.9998 | 607 | 514 | 505 |
| `D25_to_D75_spearman_pval` | Flow Timing | spearman_pval | -0.0667 | 0.4917 | 0.2122 | 0.9799 | 615 | 524 | 501 |
| `n_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.0663 | 0.5121 | 0.2176 | 0.9700 | 607 | 514 | 505 |
| `D70_day_spearman_pval` | Flow Timing | spearman_pval | -0.0595 | 0.4903 | 0.2085 | 0.9484 | 615 | 524 | 501 |
| `D70_day_mk_pval` | Flow Timing | mk_pval | -0.0535 | 0.4925 | 0.2087 | 0.9679 | 615 | 524 | 501 |
| `D25_to_D75_mk_pval` | Flow Timing | mk_pval | -0.0532 | 0.4989 | 0.2136 | 0.9914 | 615 | 524 | 501 |
| `BFI_Eckhardt_spearman_pval` | Baseflow | spearman_pval | -0.0530 | 0.5432 | 0.2160 | 0.9912 | 615 | 524 | 501 |
| `BFI_Eckhardt_mk_pval` | Baseflow | mk_pval | -0.0527 | 0.5409 | 0.2191 | 0.9980 | 615 | 524 | 501 |
| `Qann_mk_pval` | Flow Volumes | mk_pval | -0.0502 | 0.5140 | 0.2188 | 0.9986 | 607 | 514 | 505 |
| `Qann_spearman_pval` | Flow Volumes | spearman_pval | -0.0491 | 0.5170 | 0.2170 | 0.9974 | 607 | 514 | 505 |
| `D1_day_mk_pval` | Flow Timing | mk_pval | -0.0490 | 0.4943 | 0.2060 | 0.9579 | 615 | 524 | 501 |
| `Q95_Q10_linear_slp` | Flow Percentiles | linear_slp | -0.0473 | 0.7487 | 0.0235 | 20.6690 | 607 | 514 | 505 |
| `D1_day_spearman_pval` | Flow Timing | spearman_pval | -0.0397 | 0.4980 | 0.2021 | 0.9455 | 615 | 524 | 501 |
| `Qspr_spearman_pval` | Flow Volumes | spearman_pval | -0.0376 | 0.5348 | 0.2127 | 0.9620 | 607 | 514 | 505 |
| `D95_day_spearman_pval` | Flow Timing | spearman_pval | -0.0357 | 0.5344 | 0.2071 | 0.9599 | 615 | 524 | 501 |
| `Q80_mk_pval` | Flow Percentiles | mk_pval | -0.0326 | 0.5339 | 0.2178 | 0.9978 | 607 | 514 | 505 |
| `Q80_spearman_pval` | Flow Percentiles | spearman_pval | -0.0318 | 0.5388 | 0.2144 | 0.9956 | 607 | 514 | 505 |
| `flashinessRB_mk_pval` | Flashiness | mk_pval | -0.0318 | 0.5610 | 0.2127 | 0.9973 | 615 | 524 | 501 |
| `flashinessRB_spearman_pval` | Flashiness | spearman_pval | -0.0316 | 0.5644 | 0.2110 | 0.9810 | 615 | 524 | 501 |
| `n_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.0307 | 0.5711 | 0.2104 | 0.9973 | 940 | 857 | 505 |
| `BFI_LyneHollick_param_mk_pval` | Baseflow | mk_pval | -0.0276 | 0.5589 | 0.2147 | 0.9825 | 725 | 667 | 520 |
| `D95_day_mk_pval` | Flow Timing | mk_pval | -0.0273 | 0.5414 | 0.2084 | 0.9937 | 615 | 524 | 501 |
| `Qspr_mk_pval` | Flow Volumes | mk_pval | -0.0264 | 0.5371 | 0.2135 | 0.9901 | 607 | 514 | 505 |
| `BFI_LyneHollick_param_spearman_pval` | Baseflow | spearman_pval | -0.0242 | 0.5631 | 0.2120 | 0.9862 | 725 | 667 | 520 |
| `b_events_spearman_pval` | Recession | spearman_pval | -0.0237 | 0.5431 | 0.2115 | 0.9795 | 1052 | 1216 | 164 |
| `Qsum_linear_slp` | Flow Volumes | linear_slp | -0.0193 | 0.7721 | 1.8787 | 2872.1548 | 607 | 514 | 505 |
| `Q5_mk_pval` | Flow Percentiles | mk_pval | -0.0181 | 0.5981 | 0.2080 | 0.9997 | 660 | 589 | 503 |
| `Q95_linear_slp` | Flow Percentiles | linear_slp | -0.0157 | 0.7466 | 0.0301 | 32.3488 | 607 | 514 | 505 |
| `FDCall_spearman_pval` | FDC | spearman_pval | -0.0156 | 0.5540 | 0.2121 | 0.9920 | 607 | 514 | 505 |
| `Q75_mk_pval` | Flow Percentiles | mk_pval | -0.0154 | 0.5489 | 0.2139 | 0.9933 | 607 | 514 | 505 |
| `b_events_mk_pval` | Recession | mk_pval | -0.0151 | 0.5462 | 0.2136 | 0.9963 | 1052 | 1216 | 164 |
| `Q1_mk_pval` | Flow Percentiles | mk_pval | -0.0139 | 0.5992 | 0.2059 | 0.9894 | 677 | 610 | 503 |
| `Qfal_mk_pval` | Flow Volumes | mk_pval | -0.0137 | 0.5478 | 0.2196 | 0.9902 | 607 | 514 | 505 |
| `Q1_spearman_pval` | Flow Percentiles | spearman_pval | -0.0125 | 0.6019 | 0.2031 | 0.9960 | 677 | 610 | 503 |
| `Q5_spearman_pval` | Flow Percentiles | spearman_pval | -0.0120 | 0.6042 | 0.2041 | 0.9938 | 660 | 589 | 503 |
| `negative_ann_spearman_pval` | Negative Flow | spearman_pval | -0.0104 | 0.5989 | 0.1539 | 0.8662 | 6536 | 6531 | 5 |
| `n_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.0094 | 0.5802 | 0.2050 | 0.9980 | 940 | 857 | 505 |
| `FDCall_mk_pval` | FDC | mk_pval | -0.0087 | 0.5544 | 0.2143 | 0.9972 | 607 | 514 | 505 |
| `log_a_events_spearman_pval` | Recession | spearman_pval | -0.0067 | 0.5382 | 0.2096 | 0.9910 | 1052 | 1216 | 164 |
| `Q75_spearman_pval` | Flow Percentiles | spearman_pval | -0.0059 | 0.5545 | 0.2105 | 0.9887 | 607 | 514 | 505 |
| `Q60_senn_slp` | Flow Percentiles | senn_slp | -0.0053 | 0.7808 | 0.0154 | 40.3571 | 607 | 514 | 505 |
| `BFI_LyneHollick_mk_pval` | Baseflow | mk_pval | -0.0053 | 0.5652 | 0.2153 | 0.9867 | 615 | 524 | 501 |
| `Qfal_spearman_pval` | Flow Volumes | spearman_pval | -0.0029 | 0.5521 | 0.2166 | 0.9777 | 607 | 514 | 505 |
| `BFI_LyneHollick_spearman_pval` | Baseflow | spearman_pval | -0.0028 | 0.5679 | 0.2130 | 0.9921 | 615 | 524 | 501 |
| `BFI_Eckhardt_param_spearman_pval` | Baseflow | spearman_pval | -0.0021 | 0.5533 | 0.2108 | 0.9828 | 725 | 667 | 520 |
| `BFI_Eckhardt_param_mk_pval` | Baseflow | mk_pval | -0.0016 | 0.5523 | 0.2131 | 0.9900 | 725 | 667 | 520 |
| `n_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.0034 | 0.5344 | 0.2099 | 0.9686 | 607 | 514 | 505 |
| `log_a_events_mk_pval` | Recession | mk_pval | 0.0052 | 0.5426 | 0.2114 | 0.9867 | 1052 | 1216 | 164 |
| `FDCmid_mk_pval` | FDC | mk_pval | 0.0065 | 0.5496 | 0.2104 | 0.9834 | 607 | 514 | 505 |
| `FDCmid_spearman_pval` | FDC | spearman_pval | 0.0077 | 0.5494 | 0.2077 | 0.9838 | 607 | 514 | 505 |
| `dur_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.0083 | 0.5360 | 0.2091 | 0.9939 | 617 | 527 | 500 |
| `dur_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.0123 | 0.5374 | 0.2114 | 0.9820 | 617 | 527 | 500 |
| `Q10_mk_pval` | Flow Percentiles | mk_pval | 0.0132 | 0.6011 | 0.2072 | 0.9999 | 649 | 575 | 506 |
| `n_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.0145 | 0.5379 | 0.2122 | 0.9822 | 607 | 514 | 505 |
| `Flow_Reversals_summer_mk_pval` | Pulse Metrics | mk_pval | 0.0180 | 0.6115 | 0.2020 | 1.0000 | 607 | 514 | 505 |
| `Q10_spearman_pval` | Flow Percentiles | spearman_pval | 0.0183 | 0.6074 | 0.2043 | 0.9866 | 649 | 575 | 506 |
| `Flow_Reversals_summer_spearman_pval` | Pulse Metrics | spearman_pval | 0.0242 | 0.6126 | 0.1990 | 0.9921 | 607 | 514 | 505 |
| `Q70_spearman_pval` | Flow Percentiles | spearman_pval | 0.0251 | 0.5688 | 0.2069 | 0.9884 | 608 | 515 | 505 |
| `TQmean_spearman_pval` | Pulse Metrics | spearman_pval | 0.0262 | 0.5490 | 0.2056 | 0.9880 | 607 | 514 | 505 |
| `D80_day_spearman_pval` | Flow Timing | spearman_pval | 0.0275 | 0.5412 | 0.1980 | 0.9692 | 615 | 524 | 501 |
| `Q70_mk_pval` | Flow Percentiles | mk_pval | 0.0276 | 0.5668 | 0.2096 | 0.9991 | 608 | 515 | 505 |
| `Flow_Reversals_winter_mk_pval` | Pulse Metrics | mk_pval | 0.0287 | 0.5998 | 0.1997 | 0.9996 | 609 | 519 | 508 |
| `concavity_spearman_pval` | Recession | spearman_pval | 0.0294 | 0.5574 | 0.2060 | 0.9879 | 1052 | 1216 | 164 |
| `Q20_mk_pval` | Flow Percentiles | mk_pval | 0.0336 | 0.6034 | 0.2076 | 0.9972 | 633 | 548 | 505 |
| `TQmean_mk_pval` | Pulse Metrics | mk_pval | 0.0368 | 0.5531 | 0.2085 | 0.9816 | 607 | 514 | 505 |
| `b_pointcloud_spearman_pval` | Recession | spearman_pval | 0.0392 | 0.5468 | 0.1984 | 0.9824 | 1233 | 1342 | 109 |
| `D80_day_mk_pval` | Flow Timing | mk_pval | 0.0399 | 0.5462 | 0.1993 | 0.9805 | 615 | 524 | 501 |
| `concavity_mk_pval` | Recession | mk_pval | 0.0411 | 0.5632 | 0.2071 | 0.9902 | 1052 | 1216 | 164 |
| `Flow_Reversals_winter_spearman_pval` | Pulse Metrics | spearman_pval | 0.0413 | 0.6052 | 0.1973 | 0.9939 | 609 | 519 | 508 |
| `Qwin_spearman_pval` | Flow Volumes | spearman_pval | 0.0462 | 0.5782 | 0.2081 | 0.9911 | 607 | 514 | 505 |
| `negative_ann_mk_pval` | Negative Flow | mk_pval | 0.0472 | 0.5946 | 0.1583 | 0.8502 | 6536 | 6531 | 5 |
| `Q20_spearman_pval` | Flow Percentiles | spearman_pval | 0.0478 | 0.6114 | 0.2043 | 0.9919 | 633 | 548 | 505 |
| `Q60_spearman_pval` | Flow Percentiles | spearman_pval | 0.0526 | 0.5878 | 0.2037 | 0.9981 | 608 | 515 | 505 |
| `Q25_mk_pval` | Flow Percentiles | mk_pval | 0.0586 | 0.6094 | 0.2038 | 0.9997 | 631 | 546 | 507 |
| `Q60_mk_pval` | Flow Percentiles | mk_pval | 0.0614 | 0.5880 | 0.2059 | 0.9905 | 608 | 515 | 505 |
| `Q99_linear_slp` | Flow Percentiles | linear_slp | 0.0616 | 0.7701 | 0.0473 | 22.7524 | 607 | 514 | 505 |
| `D90_day_spearman_pval` | Flow Timing | spearman_pval | 0.0624 | 0.5714 | 0.1972 | 0.9696 | 615 | 524 | 501 |
| `Qwin_mk_pval` | Flow Volumes | mk_pval | 0.0625 | 0.5833 | 0.2090 | 0.9874 | 607 | 514 | 505 |
| `Q25_spearman_pval` | Flow Percentiles | spearman_pval | 0.0633 | 0.6142 | 0.2021 | 0.9925 | 631 | 546 | 507 |
| `Flow_Reversals_spring_mk_pval` | Pulse Metrics | mk_pval | 0.0650 | 0.6032 | 0.2035 | 0.9991 | 608 | 516 | 504 |
| `D90_day_mk_pval` | Flow Timing | mk_pval | 0.0696 | 0.5732 | 0.1973 | 0.9891 | 615 | 524 | 501 |
| `Flow_Reversals_spring_spearman_pval` | Pulse Metrics | spearman_pval | 0.0738 | 0.6075 | 0.2015 | 0.9675 | 608 | 516 | 504 |
| `log_a_pointcloud_spearman_pval` | Recession | spearman_pval | 0.0765 | 0.5728 | 0.1954 | 0.9742 | 1233 | 1342 | 109 |
| `Q30_spearman_pval` | Flow Percentiles | spearman_pval | 0.0864 | 0.6230 | 0.1971 | 0.9866 | 625 | 536 | 505 |
| `Q30_mk_pval` | Flow Percentiles | mk_pval | 0.0864 | 0.6200 | 0.2002 | 0.9859 | 625 | 536 | 505 |
| `b_pointcloud_mk_pval` | Recession | mk_pval | 0.0970 | 0.5733 | 0.2002 | 0.9813 | 1233 | 1342 | 109 |
| `Flow_Reversals_fall_mk_pval` | Pulse Metrics | mk_pval | 0.0983 | 0.6644 | 0.1821 | 0.9997 | 607 | 514 | 505 |
| `Flow_Reversals_fall_spearman_pval` | Pulse Metrics | spearman_pval | 0.0992 | 0.6648 | 0.1795 | 0.9882 | 607 | 514 | 505 |
| `Q50_spearman_pval` | Flow Percentiles | spearman_pval | 0.1185 | 0.6225 | 0.1931 | 0.9917 | 613 | 521 | 504 |
| `Flow_Reversals_annual_mk_pval` | Pulse Metrics | mk_pval | 0.1263 | 0.6878 | 0.1696 | 0.9986 | 607 | 514 | 505 |
| `log_a_pointcloud_mk_pval` | Recession | mk_pval | 0.1280 | 0.5967 | 0.1995 | 0.9847 | 1233 | 1342 | 109 |
| `Q90_linear_slp` | Flow Percentiles | linear_slp | 0.1308 | 0.7448 | 0.0250 | 40.2369 | 607 | 514 | 505 |
| `Q50_mk_pval` | Flow Percentiles | mk_pval | 0.1331 | 0.6244 | 0.1947 | 0.9984 | 613 | 521 | 504 |
| `n_recession_events_mk_pval` | Recession | mk_pval | 0.1334 | 0.6709 | 0.1857 | 0.9826 | 54 | 64 | 10 |
| `Flow_Reversals_annual_spearman_pval` | Pulse Metrics | spearman_pval | 0.1336 | 0.6913 | 0.1665 | 0.9956 | 607 | 514 | 505 |
| `Q40_spearman_pval` | Flow Percentiles | spearman_pval | 0.1375 | 0.6374 | 0.1926 | 0.9831 | 619 | 528 | 505 |
| `n_recession_events_spearman_pval` | Recession | spearman_pval | 0.1391 | 0.6724 | 0.1824 | 0.9955 | 54 | 64 | 10 |
| `Q40_mk_pval` | Flow Percentiles | mk_pval | 0.1442 | 0.6370 | 0.1941 | 0.9996 | 619 | 528 | 505 |
| `Qsum_senn_slp` | Flow Volumes | senn_slp | 0.1489 | 0.7662 | 1.7988 | 2665.1399 | 607 | 514 | 505 |
| `alpha_linear_spearman_pval` | Recession | spearman_pval | 0.1502 | 0.6443 | 0.1840 | 0.9832 | 1234 | 1342 | 108 |
| `Q70_linear_slp` | Flow Percentiles | linear_slp | 0.1545 | 0.7967 | 0.0140 | 32.4460 | 607 | 514 | 505 |
| `Q99_senn_slp` | Flow Percentiles | senn_slp | 0.1851 | 0.7491 | 0.0432 | 25.0000 | 607 | 514 | 505 |
| `Q90_senn_slp` | Flow Percentiles | senn_slp | 0.1869 | 0.7352 | 0.0256 | 51.3101 | 607 | 514 | 505 |
| `Qspr_linear_slp` | Flow Volumes | linear_slp | 0.1933 | 0.7686 | 1.5577 | 3150.5525 | 607 | 514 | 505 |
| `Q95_senn_slp` | Flow Percentiles | senn_slp | 0.1974 | 0.7405 | 0.0314 | 41.8500 | 607 | 514 | 505 |
| `alpha_linear_mk_pval` | Recession | mk_pval | 0.2018 | 0.6620 | 0.1881 | 0.9825 | 1234 | 1342 | 108 |
| `elasticity_rolling_spearman_rho` | Elasticity | spearman_rho | 0.2041 | 0.5527 | 0.3302 | 1.4199 | 1085 | 1085 | 0 |
| `Qspr_senn_slp` | Flow Volumes | senn_slp | 0.2045 | 0.7606 | 1.7694 | 4849.7345 | 607 | 514 | 505 |
| `Q80_linear_slp` | Flow Percentiles | linear_slp | 0.2074 | 0.7830 | 0.0176 | 37.0060 | 607 | 514 | 505 |
| `elasticity_rolling_mk_rho` | Elasticity | mk_rho | 0.2375 | 0.5557 | 0.2418 | 1.1086 | 1085 | 1085 | 0 |
| `annual_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.2389 | 0.7684 | 0.0020 | 2.0050 | 1509 | 1466 | 417 |
| `Q75_linear_slp` | Flow Percentiles | linear_slp | 0.2408 | 0.7955 | 0.0151 | 35.2612 | 607 | 514 | 505 |
| `Q95_Q10_senn_slp` | Flow Percentiles | senn_slp | 0.2603 | 0.7453 | 0.0226 | 21.3900 | 607 | 514 | 505 |
| `Q80_senn_slp` | Flow Percentiles | senn_slp | 0.2802 | 0.7640 | 0.0181 | 41.8299 | 607 | 514 | 505 |
| `elasticity_rolling_senn_slp` | Elasticity | senn_slp | 0.2966 | 0.5236 | 0.0317 | 0.4640 | 1085 | 1085 | 0 |
| `elasticity_rolling_linear_slp` | Elasticity | linear_slp | 0.3345 | 0.5440 | 0.0331 | 0.4354 | 1085 | 1085 | 0 |
| `fall_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.3405 | 0.7387 | 0.0041 | 3.7436 | 1513 | 1469 | 420 |
| `annual_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.3613 | 0.7825 | 0.0018 | 1.2977 | 1509 | 1466 | 417 |
| `elasticity_annual_spearman_rho` | Elasticity | spearman_rho | 0.3643 | 0.6289 | 0.1223 | 0.6638 | 1085 | 1085 | 0 |
| `elasticity_annual_mk_rho` | Elasticity | mk_rho | 0.3746 | 0.6304 | 0.0845 | 0.4223 | 1085 | 1085 | 0 |
| `Q75_senn_slp` | Flow Percentiles | senn_slp | 0.3765 | 0.7760 | 0.0166 | 41.5000 | 607 | 514 | 505 |
| `Q70_senn_slp` | Flow Percentiles | senn_slp | 0.3780 | 0.7723 | 0.0161 | 42.7495 | 607 | 514 | 505 |
| `summer_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.4008 | 0.6889 | 0.1199 | 0.6362 | 1688 | 1665 | 419 |
| `elasticity_annual_senn_slp` | Elasticity | senn_slp | 0.4038 | 0.6045 | 0.0332 | 0.3522 | 1085 | 1085 | 0 |
| `summer_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.4159 | 0.6906 | 0.0821 | 0.4238 | 1688 | 1665 | 419 |
| `avg_storage_senn_slp` | Storage | senn_slp | 0.4182 | 0.6705 | 2.6678 | 1403.6658 | 1529 | 1486 | 421 |
| `qp_bimodality_spearman_rho` | Q-P Seasonality | spearman_rho | 0.4254 | 0.6858 | 0.1104 | 0.5937 | 1517 | 1479 | 422 |
| `fall_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.4381 | 0.6903 | 0.0026 | 2.4167 | 1513 | 1469 | 420 |
| `qp_bimodality_mk_rho` | Q-P Seasonality | mk_rho | 0.4420 | 0.6904 | 0.0760 | 0.3817 | 1517 | 1479 | 422 |
| `D20_day_spearman_rho` | Flow Timing | spearman_rho | 0.4496 | 0.6968 | 0.0856 | 0.5317 | 615 | 524 | 501 |
| `D99_day_spearman_rho` | Flow Timing | spearman_rho | 0.4567 | 0.7284 | 0.0971 | 0.8971 | 615 | 524 | 501 |
| `avg_storage_spearman_rho` | Storage | spearman_rho | 0.4578 | 0.6770 | 0.1153 | 0.7562 | 1529 | 1486 | 421 |
| `D99_day_mk_rho` | Flow Timing | mk_rho | 0.4697 | 0.7331 | 0.0689 | 0.5850 | 615 | 524 | 501 |
| `D20_day_mk_rho` | Flow Timing | mk_rho | 0.4700 | 0.7049 | 0.0582 | 0.3671 | 615 | 524 | 501 |
| `Dmax_spearman_rho` | Flow Timing | spearman_rho | 0.4721 | 0.6859 | 0.0884 | 0.5226 | 615 | 524 | 501 |
| `avg_storage_mk_rho` | Storage | mk_rho | 0.4740 | 0.6796 | 0.0796 | 0.4969 | 1529 | 1486 | 421 |
| `D50_day_spearman_rho` | Flow Timing | spearman_rho | 0.4817 | 0.7115 | 0.0847 | 0.5642 | 615 | 524 | 501 |
| `Dmax_mk_rho` | Flow Timing | mk_rho | 0.4825 | 0.6878 | 0.0607 | 0.3662 | 615 | 524 | 501 |
| `D40_day_spearman_rho` | Flow Timing | spearman_rho | 0.4846 | 0.7154 | 0.0827 | 0.5760 | 615 | 524 | 501 |
| `dur_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.4855 | 0.7273 | 0.0152 | 0.8182 | 617 | 527 | 500 |
| `D40_day_mk_rho` | Flow Timing | mk_rho | 0.4871 | 0.7117 | 0.0573 | 0.4243 | 615 | 524 | 501 |
| `D50_day_mk_rho` | Flow Timing | mk_rho | 0.4874 | 0.7099 | 0.0582 | 0.3867 | 615 | 524 | 501 |
| `D30_day_spearman_rho` | Flow Timing | spearman_rho | 0.4929 | 0.7255 | 0.0795 | 0.4901 | 615 | 524 | 501 |
| `qp_slope_sd_linear_slp` | Q-P Seasonality | linear_slp | 0.4949 | 0.5710 | 0.1193 | 11.6260 | 1510 | 1471 | 421 |
| `fall_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.4955 | 0.7315 | 0.1233 | 0.8972 | 1513 | 1469 | 420 |
| `D30_day_mk_rho` | Flow Timing | mk_rho | 0.4970 | 0.7257 | 0.0549 | 0.3708 | 615 | 524 | 501 |
| `D10_day_spearman_rho` | Flow Timing | spearman_rho | 0.4988 | 0.7185 | 0.0841 | 0.6612 | 615 | 524 | 501 |
| `qp_bimodality_senn_slp` | Q-P Seasonality | senn_slp | 0.5043 | 0.6975 | 0.0014 | 0.0092 | 1517 | 1479 | 422 |
| `fall_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5045 | 0.7286 | 0.0851 | 0.5975 | 1513 | 1469 | 420 |
| `qp_bimodality_linear_slp` | Q-P Seasonality | linear_slp | 0.5107 | 0.7031 | 0.0014 | 0.0105 | 1517 | 1479 | 422 |
| `D60_day_spearman_rho` | Flow Timing | spearman_rho | 0.5108 | 0.7294 | 0.0841 | 0.5444 | 615 | 524 | 501 |
| `D5_day_spearman_rho` | Flow Timing | spearman_rho | 0.5131 | 0.7199 | 0.0841 | 0.6196 | 615 | 524 | 501 |
| `D10_day_mk_rho` | Flow Timing | mk_rho | 0.5147 | 0.7242 | 0.0578 | 0.4028 | 615 | 524 | 501 |
| `qp_slope_sd_spearman_rho` | Q-P Seasonality | spearman_rho | 0.5151 | 0.7209 | 0.1240 | 0.5930 | 1510 | 1471 | 421 |
| `D60_day_mk_rho` | Flow Timing | mk_rho | 0.5210 | 0.7309 | 0.0580 | 0.3839 | 615 | 524 | 501 |
| `D5_day_mk_rho` | Flow Timing | mk_rho | 0.5246 | 0.7213 | 0.0584 | 0.4139 | 615 | 524 | 501 |
| `qp_slope_sd_mk_rho` | Q-P Seasonality | mk_rho | 0.5267 | 0.7244 | 0.0855 | 0.3987 | 1510 | 1471 | 421 |
| `dur_low_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.5383 | 0.7372 | 0.1110 | 0.7304 | 1146 | 1024 | 542 |
| `dur_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.5415 | 0.7071 | 0.0404 | 1.7500 | 1443 | 1316 | 665 |
| `Q99_spearman_rho` | Flow Percentiles | spearman_rho | 0.5439 | 0.7635 | 0.0992 | 0.6567 | 607 | 514 | 505 |
| `winter_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5454 | 0.7522 | 0.1254 | 0.8049 | 1511 | 1467 | 416 |
| `dur_low_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5471 | 0.7443 | 0.1132 | 0.8850 | 4299 | 4080 | 751 |
| `D1_day_spearman_rho` | Flow Timing | spearman_rho | 0.5497 | 0.7380 | 0.0831 | 0.6098 | 615 | 524 | 501 |
| `D95_day_spearman_rho` | Flow Timing | spearman_rho | 0.5518 | 0.7659 | 0.0880 | 0.7388 | 615 | 524 | 501 |
| `D70_day_spearman_rho` | Flow Timing | spearman_rho | 0.5518 | 0.7520 | 0.0833 | 0.5860 | 615 | 524 | 501 |
| `dur_low_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.5543 | 0.7415 | 0.0787 | 0.5009 | 1146 | 1024 | 542 |
| `Q90_spearman_rho` | Flow Percentiles | spearman_rho | 0.5572 | 0.7704 | 0.1024 | 0.7223 | 607 | 514 | 505 |
| `Q95_Q10_spearman_rho` | Flow Percentiles | spearman_rho | 0.5575 | 0.7692 | 0.0996 | 0.6229 | 607 | 514 | 505 |
| `Q95_spearman_rho` | Flow Percentiles | spearman_rho | 0.5579 | 0.7693 | 0.1009 | 0.6783 | 607 | 514 | 505 |
| `Q99_mk_rho` | Flow Percentiles | mk_rho | 0.5589 | 0.7660 | 0.0684 | 0.4176 | 607 | 514 | 505 |
| `winter_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5591 | 0.7556 | 0.0863 | 0.5450 | 1511 | 1467 | 416 |
| `D25_to_D75_spearman_rho` | Flow Timing | spearman_rho | 0.5593 | 0.7633 | 0.0894 | 0.6135 | 615 | 524 | 501 |
| `dur_low_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5613 | 0.7459 | 0.0789 | 0.6852 | 4299 | 4080 | 751 |
| `log_a_seasonality_minimum_first_half` | Recession | scalar | 0.5616 | 0.8041 | 19.3656 | 361.1328 | 1053 | 1218 | 165 |
| `spring_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5616 | 0.7670 | 0.1198 | 0.6438 | 1514 | 1469 | 419 |
| `D1_day_mk_rho` | Flow Timing | mk_rho | 0.5637 | 0.7413 | 0.0591 | 0.4372 | 615 | 524 | 501 |
| `D70_day_mk_rho` | Flow Timing | mk_rho | 0.5639 | 0.7548 | 0.0573 | 0.3716 | 615 | 524 | 501 |
| `negative_ann_linear_slp` | Negative Flow | linear_slp | 0.5644 | 0.7675 | 0.0003 | 0.3527 | 607 | 514 | 505 |
| `D25_to_D75_mk_rho` | Flow Timing | mk_rho | 0.5649 | 0.7634 | 0.0618 | 0.4248 | 615 | 524 | 501 |
| `D95_day_mk_rho` | Flow Timing | mk_rho | 0.5664 | 0.7714 | 0.0607 | 0.4861 | 615 | 524 | 501 |
| `n_low_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.5666 | 0.7526 | 0.1151 | 0.7950 | 650 | 576 | 506 |
| `D60_day_senn_slp` | Flow Timing | senn_slp | 0.5670 | 0.7017 | 0.2363 | 3.2725 | 615 | 524 | 501 |
| `dur_high_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5686 | 0.7542 | 0.1000 | 0.5810 | 1443 | 1316 | 665 |
| `Qann_spearman_rho` | Flow Volumes | spearman_rho | 0.5721 | 0.7885 | 0.1013 | 0.6248 | 607 | 514 | 505 |
| `dur_high_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5743 | 0.7563 | 0.0697 | 0.4106 | 1443 | 1316 | 665 |
| `qp_slope_sd_senn_slp` | Q-P Seasonality | senn_slp | 0.5764 | 0.6398 | 0.0138 | 1.7796 | 1510 | 1471 | 421 |
| `n_low_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.5769 | 0.7544 | 0.0843 | 0.6395 | 650 | 576 | 506 |
| `Q95_Q10_mk_rho` | Flow Percentiles | mk_rho | 0.5783 | 0.7755 | 0.0688 | 0.4035 | 607 | 514 | 505 |
| `spring_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5783 | 0.7729 | 0.0822 | 0.4432 | 1514 | 1469 | 419 |
| `n_high_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5787 | 0.7801 | 0.0732 | 0.5558 | 607 | 514 | 505 |
| `concavity_linear_slp` | Recession | linear_slp | 0.5789 | 0.7405 | 0.0283 | 0.9138 | 1052 | 1216 | 164 |
| `FDC90th_spearman_rho` | FDC | spearman_rho | 0.5795 | 0.7638 | 0.1191 | 0.8541 | 652 | 578 | 506 |
| `n_high_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5804 | 0.7830 | 0.1013 | 0.7458 | 607 | 514 | 505 |
| `Q95_mk_rho` | Flow Percentiles | mk_rho | 0.5808 | 0.7751 | 0.0698 | 0.4400 | 607 | 514 | 505 |
| `D5_day_senn_slp` | Flow Timing | senn_slp | 0.5818 | 0.7022 | 0.1926 | 4.3639 | 615 | 524 | 501 |
| `BFI_Eckhardt_param_spearman_rho` | Baseflow | spearman_rho | 0.5827 | 0.7838 | 0.1063 | 0.9776 | 725 | 667 | 520 |
| `Qfal_spearman_rho` | Flow Volumes | spearman_rho | 0.5837 | 0.7485 | 0.1091 | 0.8724 | 607 | 514 | 505 |
| `Q90_mk_rho` | Flow Percentiles | mk_rho | 0.5838 | 0.7771 | 0.0710 | 0.5087 | 607 | 514 | 505 |
| `D50_day_senn_slp` | Flow Timing | senn_slp | 0.5840 | 0.6940 | 0.2385 | 3.0440 | 615 | 524 | 501 |
| `dur_high_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.5856 | 0.7161 | 0.0352 | 0.6946 | 617 | 527 | 500 |
| `D70_day_senn_slp` | Flow Timing | senn_slp | 0.5880 | 0.7348 | 0.2339 | 4.0111 | 615 | 524 | 501 |
| `Qann_mk_rho` | Flow Volumes | mk_rho | 0.5899 | 0.7916 | 0.0701 | 0.4039 | 607 | 514 | 505 |
| `FDC90th_mk_rho` | FDC | mk_rho | 0.5905 | 0.7660 | 0.0830 | 0.5246 | 652 | 578 | 506 |
| `BFI_Eckhardt_param_mk_rho` | Baseflow | mk_rho | 0.5916 | 0.7854 | 0.0735 | 0.6833 | 725 | 667 | 520 |
| `Qsum_spearman_rho` | Flow Volumes | spearman_rho | 0.5938 | 0.7807 | 0.0994 | 0.7338 | 607 | 514 | 505 |
| `annual_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5946 | 0.7868 | 0.1254 | 0.8608 | 1509 | 1466 | 417 |
| `b_events_linear_slp` | Recession | linear_slp | 0.5953 | 0.7305 | 0.0183 | 0.9240 | 1052 | 1216 | 164 |
| `D80_day_spearman_rho` | Flow Timing | spearman_rho | 0.5959 | 0.7749 | 0.0813 | 0.6596 | 615 | 524 | 501 |
| `D40_day_senn_slp` | Flow Timing | senn_slp | 0.5981 | 0.7040 | 0.2565 | 3.2889 | 615 | 524 | 501 |
| `elasticity_annual_linear_slp` | Elasticity | linear_slp | 0.5991 | 0.6125 | 0.3242 | 26.9528 | 1085 | 1085 | 0 |
| `Q80_spearman_rho` | Flow Percentiles | spearman_rho | 0.6007 | 0.7990 | 0.1016 | 0.7460 | 607 | 514 | 505 |
| `BFI_Eckhardt_spearman_rho` | Baseflow | spearman_rho | 0.6010 | 0.7901 | 0.1083 | 0.9792 | 615 | 524 | 501 |
| `Qfal_mk_rho` | Flow Volumes | mk_rho | 0.6043 | 0.7494 | 0.0756 | 0.5377 | 607 | 514 | 505 |
| `BFI_LyneHollick_param_spearman_rho` | Baseflow | spearman_rho | 0.6056 | 0.7975 | 0.1073 | 0.8648 | 725 | 667 | 520 |
| `concavity_spearman_rho` | Recession | spearman_rho | 0.6056 | 0.7718 | 0.1040 | 0.6158 | 1052 | 1216 | 164 |
| `D60_day_linear_slp` | Flow Timing | linear_slp | 0.6058 | 0.7249 | 0.2344 | 3.3036 | 615 | 524 | 501 |
| `b_events_spearman_rho` | Recession | spearman_rho | 0.6061 | 0.7759 | 0.1099 | 0.6756 | 1052 | 1216 | 164 |
| `n_recession_events_mk_rho` | Recession | mk_rho | 0.6065 | 0.8032 | 0.0800 | 0.5869 | 54 | 64 | 10 |
| `Qsum_mk_rho` | Flow Volumes | mk_rho | 0.6066 | 0.7806 | 0.0684 | 0.4679 | 608 | 515 | 505 |
| `Dmax_senn_slp` | Flow Timing | senn_slp | 0.6066 | 0.7009 | 0.4989 | 8.3080 | 615 | 524 | 501 |
| `D10_day_senn_slp` | Flow Timing | senn_slp | 0.6070 | 0.7155 | 0.2478 | 2.8933 | 615 | 524 | 501 |
| `flashinessRB_spearman_rho` | Flashiness | spearman_rho | 0.6081 | 0.7994 | 0.1160 | 0.8695 | 615 | 524 | 501 |
| `n_recession_events_spearman_rho` | Recession | spearman_rho | 0.6081 | 0.8065 | 0.1058 | 0.8510 | 54 | 64 | 10 |
| `D50_day_linear_slp` | Flow Timing | linear_slp | 0.6083 | 0.7024 | 0.2487 | 3.3354 | 615 | 524 | 501 |
| `BFI_LyneHollick_spearman_rho` | Baseflow | spearman_rho | 0.6089 | 0.7977 | 0.1091 | 0.8944 | 615 | 524 | 501 |
| `log_a_events_spearman_rho` | Recession | spearman_rho | 0.6089 | 0.7780 | 0.1042 | 0.8533 | 1052 | 1216 | 164 |
| `annual_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.6093 | 0.7878 | 0.0878 | 0.5443 | 1509 | 1466 | 417 |
| `D90_day_spearman_rho` | Flow Timing | spearman_rho | 0.6114 | 0.7825 | 0.0826 | 0.6579 | 615 | 524 | 501 |
| `D80_day_mk_rho` | Flow Timing | mk_rho | 0.6120 | 0.7804 | 0.0558 | 0.4318 | 615 | 524 | 501 |
| `Dmax_linear_slp` | Flow Timing | linear_slp | 0.6126 | 0.7114 | 0.6048 | 9.4977 | 615 | 524 | 501 |
| `n_low_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.6130 | 0.7876 | 0.1273 | 1.1851 | 940 | 857 | 505 |
| `FDCmid_spearman_rho` | FDC | spearman_rho | 0.6134 | 0.7771 | 0.1011 | 0.7671 | 607 | 514 | 505 |
| `n_low_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.6139 | 0.7856 | 0.0943 | 0.8463 | 940 | 857 | 505 |
| `BFI_Eckhardt_mk_rho` | Baseflow | mk_rho | 0.6145 | 0.7930 | 0.0749 | 0.6716 | 615 | 524 | 501 |
| `BFI_LyneHollick_param_mk_rho` | Baseflow | mk_rho | 0.6170 | 0.8003 | 0.0743 | 0.5597 | 725 | 667 | 520 |
| `flashinessRB_mk_rho` | Flashiness | mk_rho | 0.6174 | 0.8021 | 0.0806 | 0.6253 | 615 | 524 | 501 |
| `D80_day_senn_slp` | Flow Timing | senn_slp | 0.6175 | 0.7666 | 0.2336 | 3.4375 | 615 | 524 | 501 |
| `concavity_mk_rho` | Recession | mk_rho | 0.6176 | 0.7754 | 0.0720 | 0.4163 | 1052 | 1216 | 164 |
| `b_events_mk_rho` | Recession | mk_rho | 0.6189 | 0.7797 | 0.0761 | 0.4705 | 1052 | 1216 | 164 |
| `FDCall_spearman_rho` | FDC | spearman_rho | 0.6190 | 0.7822 | 0.1095 | 0.8746 | 607 | 514 | 505 |
| `Q80_mk_rho` | Flow Percentiles | mk_rho | 0.6195 | 0.8007 | 0.0708 | 0.5182 | 607 | 514 | 505 |
| `BFI_LyneHollick_mk_rho` | Baseflow | mk_rho | 0.6199 | 0.8003 | 0.0757 | 0.6150 | 615 | 524 | 501 |
| `D40_day_linear_slp` | Flow Timing | linear_slp | 0.6201 | 0.7000 | 0.2678 | 3.0155 | 615 | 524 | 501 |
| `log_a_events_mk_rho` | Recession | mk_rho | 0.6207 | 0.7808 | 0.0721 | 0.5748 | 1052 | 1216 | 164 |
| `Q75_spearman_rho` | Flow Percentiles | spearman_rho | 0.6210 | 0.8110 | 0.1027 | 0.7774 | 607 | 514 | 505 |
| `D70_day_linear_slp` | Flow Timing | linear_slp | 0.6223 | 0.7584 | 0.2327 | 3.2429 | 615 | 524 | 501 |
| `b_pointcloud_spearman_rho` | Recession | spearman_rho | 0.6232 | 0.7825 | 0.1159 | 1.4455 | 1233 | 1342 | 109 |
| `Flow_Reversals_summer_spearman_rho` | Pulse Metrics | spearman_rho | 0.6240 | 0.8065 | 0.1176 | 0.8152 | 607 | 514 | 505 |
| `concavity_senn_slp` | Recession | senn_slp | 0.6258 | 0.7396 | 0.0182 | 0.5672 | 1052 | 1216 | 164 |
| `dur_high_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.6259 | 0.7914 | 0.0945 | 0.6076 | 617 | 527 | 500 |
| `FDCmid_mk_rho` | FDC | mk_rho | 0.6261 | 0.7784 | 0.0699 | 0.5080 | 607 | 514 | 505 |
| `D90_day_mk_rho` | Flow Timing | mk_rho | 0.6263 | 0.7887 | 0.0568 | 0.4465 | 615 | 524 | 501 |
| `D20_day_linear_slp` | Flow Timing | linear_slp | 0.6268 | 0.6971 | 0.3039 | 3.2940 | 615 | 524 | 501 |
| `D5_day_linear_slp` | Flow Timing | linear_slp | 0.6274 | 0.7120 | 0.2312 | 3.9052 | 615 | 524 | 501 |
| `D20_day_senn_slp` | Flow Timing | senn_slp | 0.6279 | 0.7067 | 0.2813 | 3.6327 | 615 | 524 | 501 |
| `Flow_Reversals_fall_spearman_rho` | Pulse Metrics | spearman_rho | 0.6287 | 0.8110 | 0.1202 | 0.9492 | 607 | 514 | 505 |
| `Q70_spearman_rho` | Flow Percentiles | spearman_rho | 0.6288 | 0.8137 | 0.1041 | 0.8040 | 608 | 515 | 505 |
| `n_high_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.6304 | 0.7937 | 0.0950 | 0.6214 | 607 | 514 | 505 |
| `D1_day_senn_slp` | Flow Timing | senn_slp | 0.6308 | 0.6751 | 0.0856 | 2.9195 | 615 | 524 | 501 |
| `FDCall_mk_rho` | FDC | mk_rho | 0.6317 | 0.7838 | 0.0760 | 0.5447 | 607 | 514 | 505 |
| `TQmean_spearman_rho` | Pulse Metrics | spearman_rho | 0.6317 | 0.7898 | 0.0965 | 0.6937 | 607 | 514 | 505 |
| `dur_high_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.6325 | 0.7930 | 0.0677 | 0.4136 | 617 | 527 | 500 |
| `D30_day_senn_slp` | Flow Timing | senn_slp | 0.6336 | 0.7271 | 0.2675 | 3.1180 | 615 | 524 | 501 |
| `b_events_senn_slp` | Recession | senn_slp | 0.6339 | 0.7493 | 0.0111 | 0.3449 | 1052 | 1216 | 164 |
| `dur_high_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.6339 | 0.7206 | 0.0552 | 1.8284 | 1443 | 1316 | 665 |
| `D30_day_linear_slp` | Flow Timing | linear_slp | 0.6346 | 0.7079 | 0.2860 | 3.0639 | 615 | 524 | 501 |
| `n_high_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.6348 | 0.7949 | 0.0693 | 0.4684 | 607 | 514 | 505 |
| `D10_day_linear_slp` | Flow Timing | linear_slp | 0.6349 | 0.7200 | 0.2783 | 3.6144 | 615 | 524 | 501 |
| `Qwin_spearman_rho` | Flow Volumes | spearman_rho | 0.6364 | 0.8074 | 0.1000 | 0.9064 | 607 | 514 | 505 |
| `Flow_Reversals_summer_mk_rho` | Pulse Metrics | mk_rho | 0.6365 | 0.8103 | 0.0835 | 0.5916 | 607 | 514 | 505 |
| `D80_day_linear_slp` | Flow Timing | linear_slp | 0.6375 | 0.7685 | 0.2360 | 3.4114 | 615 | 524 | 501 |
| `Qspr_spearman_rho` | Flow Volumes | spearman_rho | 0.6377 | 0.7980 | 0.0953 | 0.6402 | 607 | 514 | 505 |
| `log_a_pointcloud_spearman_rho` | Recession | spearman_rho | 0.6378 | 0.7928 | 0.1167 | 1.4818 | 1233 | 1342 | 109 |
| `winter_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.6379 | 0.7378 | 0.0035 | 3.5605 | 1511 | 1467 | 416 |
| `Q75_mk_rho` | Flow Percentiles | mk_rho | 0.6381 | 0.8122 | 0.0716 | 0.5212 | 607 | 514 | 505 |
| `n_low_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.6391 | 0.7171 | 0.0234 | 0.6667 | 607 | 514 | 505 |
| `Flow_Reversals_spring_spearman_rho` | Pulse Metrics | spearman_rho | 0.6411 | 0.8019 | 0.1107 | 0.7447 | 608 | 516 | 504 |
| `TQmean_mk_rho` | Pulse Metrics | mk_rho | 0.6415 | 0.7920 | 0.0666 | 0.4600 | 607 | 514 | 505 |
| `Flow_Reversals_fall_mk_rho` | Pulse Metrics | mk_rho | 0.6415 | 0.8141 | 0.0860 | 0.6473 | 607 | 514 | 505 |
| `log_a_events_senn_slp` | Recession | senn_slp | 0.6437 | 0.7505 | 0.0178 | 0.7576 | 1052 | 1216 | 164 |
| `Q70_mk_rho` | Flow Percentiles | mk_rho | 0.6445 | 0.8142 | 0.0724 | 0.5478 | 608 | 515 | 505 |
| `FDC90th_senn_slp` | FDC | senn_slp | 0.6447 | 0.7408 | 0.0192 | 1.9703 | 607 | 514 | 505 |
| `D95_day_senn_slp` | Flow Timing | senn_slp | 0.6469 | 0.7507 | 0.1980 | 2.7325 | 615 | 524 | 501 |
| `Flow_Reversals_spring_mk_rho` | Pulse Metrics | mk_rho | 0.6470 | 0.8030 | 0.0790 | 0.4926 | 608 | 516 | 504 |
| `Flow_Reversals_winter_senn_slp` | Pulse Metrics | senn_slp | 0.6478 | 0.7972 | 0.0821 | 1.2102 | 607 | 514 | 505 |
| `b_pointcloud_mk_rho` | Recession | mk_rho | 0.6486 | 0.7883 | 0.0822 | 1.2364 | 1233 | 1342 | 109 |
| `Flow_Reversals_winter_spearman_rho` | Pulse Metrics | spearman_rho | 0.6492 | 0.8131 | 0.1232 | 0.8464 | 609 | 519 | 508 |
| `Q60_spearman_rho` | Flow Percentiles | spearman_rho | 0.6499 | 0.8140 | 0.1051 | 0.7858 | 608 | 515 | 505 |
| `n_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.6502 | 0.7367 | 0.0189 | 0.4706 | 607 | 514 | 505 |
| `D25_to_D75_linear_slp` | Flow Timing | linear_slp | 0.6517 | 0.7555 | 0.3192 | 4.1164 | 615 | 524 | 501 |
| `D1_day_linear_slp` | Flow Timing | linear_slp | 0.6519 | 0.7011 | 0.1282 | 4.4434 | 615 | 524 | 501 |
| `D90_day_senn_slp` | Flow Timing | senn_slp | 0.6525 | 0.7725 | 0.2252 | 2.7449 | 615 | 524 | 501 |
| `D90_day_linear_slp` | Flow Timing | linear_slp | 0.6539 | 0.7679 | 0.2308 | 2.9680 | 615 | 524 | 501 |
| `D25_to_D75_senn_slp` | Flow Timing | senn_slp | 0.6543 | 0.7586 | 0.3075 | 4.5714 | 615 | 524 | 501 |
| `Flow_Reversals_winter_mk_rho` | Pulse Metrics | mk_rho | 0.6544 | 0.8134 | 0.0885 | 0.5784 | 609 | 519 | 508 |
| `avg_storage_linear_slp` | Storage | linear_slp | 0.6546 | 0.6907 | 2.6335 | 1573.3937 | 1529 | 1486 | 421 |
| `Qwin_mk_rho` | Flow Volumes | mk_rho | 0.6553 | 0.8125 | 0.0690 | 0.5950 | 607 | 514 | 505 |
| `Qspr_mk_rho` | Flow Volumes | mk_rho | 0.6560 | 0.8031 | 0.0656 | 0.4148 | 607 | 514 | 505 |
| `log_a_pointcloud_mk_rho` | Recession | mk_rho | 0.6574 | 0.7981 | 0.0835 | 1.2727 | 1233 | 1342 | 109 |
| `Q1_spearman_rho` | Flow Percentiles | spearman_rho | 0.6581 | 0.8086 | 0.1290 | 1.2065 | 677 | 610 | 503 |
| `D95_day_linear_slp` | Flow Timing | linear_slp | 0.6652 | 0.7579 | 0.2063 | 2.7095 | 615 | 524 | 501 |
| `Q60_mk_rho` | Flow Percentiles | mk_rho | 0.6655 | 0.8140 | 0.0733 | 0.5332 | 608 | 515 | 505 |
| `BFI_LyneHollick_senn_slp` | Baseflow | senn_slp | 0.6663 | 0.7935 | 0.0006 | 0.0158 | 615 | 524 | 501 |
| `Flow_Reversals_winter_linear_slp` | Pulse Metrics | linear_slp | 0.6666 | 0.8099 | 0.0838 | 0.9813 | 607 | 514 | 505 |
| `BFI_LyneHollick_linear_slp` | Baseflow | linear_slp | 0.6672 | 0.7925 | 0.0007 | 0.0206 | 615 | 524 | 501 |
| `n_low_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.6680 | 0.7632 | 0.0296 | 0.6891 | 607 | 514 | 505 |
| `Q5_spearman_rho` | Flow Percentiles | spearman_rho | 0.6728 | 0.8116 | 0.1239 | 1.2552 | 660 | 589 | 503 |
| `Flow_Reversals_annual_spearman_rho` | Pulse Metrics | spearman_rho | 0.6733 | 0.8290 | 0.1300 | 0.9938 | 607 | 514 | 505 |
| `BFI_Eckhardt_linear_slp` | Baseflow | linear_slp | 0.6752 | 0.7900 | 0.0005 | 0.0148 | 615 | 524 | 501 |
| `Q1_mk_rho` | Flow Percentiles | mk_rho | 0.6756 | 0.8097 | 0.0906 | 0.7992 | 677 | 610 | 503 |
| `log_a_seasonality_amplitude_first_half` | Recession | scalar | 0.6766 | 0.8597 | 0.9419 | 75.8573 | 1053 | 1218 | 165 |
| `BFI_Eckhardt_senn_slp` | Baseflow | senn_slp | 0.6777 | 0.7871 | 0.0005 | 0.0110 | 615 | 524 | 501 |
| `Q10_spearman_rho` | Flow Percentiles | spearman_rho | 0.6786 | 0.8130 | 0.1200 | 1.0311 | 649 | 575 | 506 |
| `Q50_spearman_rho` | Flow Percentiles | spearman_rho | 0.6786 | 0.8230 | 0.1033 | 0.8183 | 613 | 521 | 504 |
| `log_a_events_linear_slp` | Recession | linear_slp | 0.6791 | 0.7308 | 0.0301 | 1.8091 | 1052 | 1216 | 164 |
| `alpha_linear_spearman_rho` | Recession | spearman_rho | 0.6824 | 0.8180 | 0.1156 | 1.6364 | 1234 | 1342 | 108 |
| `Flow_Reversals_annual_mk_rho` | Pulse Metrics | mk_rho | 0.6830 | 0.8315 | 0.0926 | 0.6879 | 607 | 514 | 505 |
| `Q30_spearman_rho` | Flow Percentiles | spearman_rho | 0.6836 | 0.8161 | 0.1098 | 0.8754 | 625 | 536 | 505 |
| `Q25_spearman_rho` | Flow Percentiles | spearman_rho | 0.6858 | 0.8159 | 0.1118 | 0.8751 | 631 | 546 | 507 |
| `FDC90th_linear_slp` | FDC | linear_slp | 0.6869 | 0.6877 | 0.0636 | 3.1092 | 607 | 514 | 505 |
| `Q5_mk_rho` | Flow Percentiles | mk_rho | 0.6871 | 0.8101 | 0.0873 | 0.8233 | 660 | 589 | 503 |
| `Q20_spearman_rho` | Flow Percentiles | spearman_rho | 0.6874 | 0.8169 | 0.1139 | 0.8192 | 633 | 548 | 505 |
| `TQmean_senn_slp` | Pulse Metrics | senn_slp | 0.6878 | 0.7777 | 0.0584 | 0.8088 | 607 | 514 | 505 |
| `Q40_spearman_rho` | Flow Percentiles | spearman_rho | 0.6884 | 0.8220 | 0.1054 | 0.8798 | 619 | 528 | 505 |
| `Q10_mk_rho` | Flow Percentiles | mk_rho | 0.6907 | 0.8097 | 0.0846 | 0.6869 | 649 | 575 | 506 |
| `Flow_Reversals_fall_senn_slp` | Pulse Metrics | senn_slp | 0.6929 | 0.8132 | 0.0713 | 1.0217 | 607 | 514 | 505 |
| `D99_day_senn_slp` | Flow Timing | senn_slp | 0.6948 | 0.6883 | 0.1011 | 4.0679 | 615 | 524 | 501 |
| `Q50_mk_rho` | Flow Percentiles | mk_rho | 0.6949 | 0.8242 | 0.0723 | 0.5226 | 613 | 521 | 504 |
| `winter_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.6950 | 0.7196 | 0.0034 | 2.1723 | 1511 | 1467 | 416 |
| `n_high_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.6977 | 0.8145 | 0.0235 | 0.4075 | 607 | 514 | 505 |
| `Q30_mk_rho` | Flow Percentiles | mk_rho | 0.6980 | 0.8148 | 0.0772 | 0.5776 | 625 | 536 | 505 |
| `n_low_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.7001 | 0.7760 | 0.0466 | 1.8651 | 607 | 514 | 505 |
| `Q20_mk_rho` | Flow Percentiles | mk_rho | 0.7002 | 0.8156 | 0.0804 | 0.5712 | 633 | 548 | 505 |
| `Q25_mk_rho` | Flow Percentiles | mk_rho | 0.7003 | 0.8156 | 0.0786 | 0.5877 | 631 | 546 | 507 |
| `TQmean_linear_slp` | Pulse Metrics | linear_slp | 0.7025 | 0.7899 | 0.0584 | 0.7953 | 607 | 514 | 505 |
| `n_low_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.7036 | 0.7525 | 0.0328 | 1.4148 | 607 | 514 | 505 |
| `Q40_mk_rho` | Flow Percentiles | mk_rho | 0.7046 | 0.8221 | 0.0739 | 0.5644 | 619 | 528 | 505 |
| `alpha_linear_mk_rho` | Recession | mk_rho | 0.7051 | 0.8199 | 0.0833 | 1.4182 | 1234 | 1342 | 108 |
| `BFI_LyneHollick_param_senn_slp` | Baseflow | senn_slp | 0.7064 | 0.7969 | 0.0006 | 0.0108 | 725 | 667 | 520 |
| `BFI_LyneHollick_param_linear_slp` | Baseflow | linear_slp | 0.7078 | 0.7964 | 0.0006 | 0.0158 | 725 | 667 | 520 |
| `D99_day_linear_slp` | Flow Timing | linear_slp | 0.7079 | 0.7095 | 0.1228 | 4.0688 | 615 | 524 | 501 |
| `Flow_Reversals_fall_linear_slp` | Pulse Metrics | linear_slp | 0.7095 | 0.8185 | 0.0713 | 1.0094 | 607 | 514 | 505 |
| `Flow_Reversals_summer_senn_slp` | Pulse Metrics | senn_slp | 0.7127 | 0.8104 | 0.0651 | 1.3000 | 607 | 514 | 505 |
| `summer_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.7137 | 0.6486 | 0.0056 | 4.9540 | 1688 | 1665 | 419 |
| `alpha_linear_linear_slp` | Recession | linear_slp | 0.7155 | 0.7832 | 0.0005 | 0.0372 | 1234 | 1342 | 108 |
| `log_a_pointcloud_linear_slp` | Recession | linear_slp | 0.7155 | 0.7403 | 0.0080 | 0.4060 | 1233 | 1342 | 109 |
| `Flow_Reversals_summer_linear_slp` | Pulse Metrics | linear_slp | 0.7205 | 0.8159 | 0.0669 | 1.2650 | 607 | 514 | 505 |
| `Flow_Reversals_spring_senn_slp` | Pulse Metrics | senn_slp | 0.7248 | 0.8062 | 0.0570 | 0.8920 | 607 | 514 | 505 |
| `BFI_Eckhardt_param_linear_slp` | Baseflow | linear_slp | 0.7248 | 0.7690 | 0.0002 | 0.0033 | 725 | 667 | 520 |
| `Flow_Reversals_annual_senn_slp` | Pulse Metrics | senn_slp | 0.7261 | 0.8368 | 0.2065 | 3.4220 | 607 | 514 | 505 |
| `flashinessRB_linear_slp` | Flashiness | linear_slp | 0.7295 | 0.7811 | 0.0007 | 0.0318 | 615 | 524 | 501 |
| `Flow_Reversals_spring_linear_slp` | Pulse Metrics | linear_slp | 0.7320 | 0.8123 | 0.0592 | 0.8214 | 607 | 514 | 505 |
| `Flow_Reversals_annual_linear_slp` | Pulse Metrics | linear_slp | 0.7345 | 0.8337 | 0.2103 | 3.3095 | 607 | 514 | 505 |
| `n_recession_events_senn_slp` | Recession | senn_slp | 0.7382 | 0.8122 | 0.0158 | 0.2229 | 0 | 0 | 0 |
| `n_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.7392 | 0.7674 | 0.0274 | 0.5516 | 607 | 514 | 505 |
| `negative_ann_spearman_rho` | Negative Flow | spearman_rho | 0.7393 | 0.9104 | 0.0624 | 0.3326 | 6536 | 6531 | 5 |
| `BFI_Eckhardt_param_senn_slp` | Baseflow | senn_slp | 0.7417 | 0.7776 | 0.0002 | 0.0032 | 725 | 667 | 520 |
| `flashinessRB_senn_slp` | Flashiness | senn_slp | 0.7425 | 0.7896 | 0.0007 | 0.0174 | 615 | 524 | 501 |
| `b_pointcloud_linear_slp` | Recession | linear_slp | 0.7451 | 0.7530 | 0.0057 | 0.2029 | 1233 | 1342 | 109 |
| `alpha_linear_senn_slp` | Recession | senn_slp | 0.7501 | 0.8090 | 0.0005 | 0.0350 | 1234 | 1342 | 108 |
| `FDCall_senn_slp` | FDC | senn_slp | 0.7505 | 0.7545 | 0.0067 | 0.3203 | 607 | 514 | 505 |
| `negative_ann_mk_rho` | Negative Flow | mk_rho | 0.7506 | 0.9146 | 0.0487 | 0.2671 | 6536 | 6531 | 5 |
| `n_high_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.7548 | 0.8240 | 0.0315 | 0.6639 | 607 | 514 | 505 |
| `n_recession_events_linear_slp` | Recession | linear_slp | 0.7565 | 0.8507 | 0.0185 | 0.2544 | 0 | 0 | 0 |
| `FDCmid_senn_slp` | FDC | senn_slp | 0.7588 | 0.7491 | 0.0061 | 0.4185 | 607 | 514 | 505 |
| `FDCmid_linear_slp` | FDC | linear_slp | 0.7603 | 0.7452 | 0.0093 | 0.5769 | 607 | 514 | 505 |
| `summer_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.7694 | 0.5732 | 0.0088 | 4.8511 | 1688 | 1665 | 419 |
| `elasticity_annual_mean` | Elasticity | mean | 0.7739 | 0.8239 | 1.9649 | 199.3512 | 1085 | 1085 | 0 |
| `b_pointcloud_senn_slp` | Recession | senn_slp | 0.7745 | 0.7902 | 0.0051 | 0.4522 | 1233 | 1342 | 109 |
| `log_a_seasonality_minimum_last_half` | Recession | scalar | 0.7767 | 0.9232 | 10.0177 | 360.9091 | 1056 | 1218 | 162 |
| `FDCall_linear_slp` | FDC | linear_slp | 0.7770 | 0.7466 | 0.0082 | 0.2994 | 607 | 514 | 505 |
| `dur_low_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.8104 | 0.7394 | 0.0540 | 1.0872 | 1146 | 1024 | 542 |
| `elasticity_static` | Elasticity | scalar | 0.8119 | 0.9004 | 0.2104 | 2.9901 | 1085 | 1085 | 0 |
| `log_a_pointcloud_senn_slp` | Recession | senn_slp | 0.8167 | 0.7868 | 0.0070 | 0.6651 | 1233 | 1342 | 109 |
| `negative_ann_senn_slp` | Negative Flow | senn_slp | 0.8202 | 0.8660 | 0.0000 | 0.1786 | 607 | 514 | 505 |
| `dur_low_pulses_all_median` | Pulse Metrics | median | 0.8237 | 0.9319 | 1.4787 | 116.3139 | 376 | 368 | 40 |
| `log_a_seasonality_minimum_all` | Recession | scalar | 0.8262 | 0.9453 | 7.1526 | 363.3182 | 1052 | 1216 | 164 |
| `elasticity_annual_median` | Elasticity | median | 0.8495 | 0.9203 | 0.1619 | 3.1040 | 1085 | 1085 | 0 |
| `runoff_ratio_high_count` | Runoff Ratios | scalar | 0.8562 | 0.9662 | 0.0364 | 12.0000 | 1085 | 1085 | 0 |
| `FDC90th_median` | FDC | median | 0.8778 | 0.9787 | 0.1978 | 40.4516 | 0 | 0 | 0 |
| `elasticity_rolling_median` | Elasticity | median | 0.8822 | 0.9368 | 0.1586 | 2.4302 | 1085 | 1085 | 0 |
| `n_low_pulses_all_median` | Pulse Metrics | median | 0.8835 | 0.9353 | 0.3842 | 11.5000 | 0 | 0 | 0 |
| `dur_low_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.8878 | 0.7160 | 0.0295 | 0.7195 | 1146 | 1024 | 542 |
| `log_a_seasonality_amplitude_last_half` | Recession | scalar | 0.9031 | 0.9613 | 0.5808 | 55.2083 | 1056 | 1218 | 162 |
| `dur_low_pulses_all_mean` | Pulse Metrics | mean | 0.9095 | 0.9637 | 1.4866 | 85.8383 | 376 | 368 | 40 |
| `elasticity_rolling_mean` | Elasticity | mean | 0.9170 | 0.9534 | 0.1542 | 1.9362 | 1085 | 1085 | 0 |
| `log_a_seasonality_amplitude_all` | Recession | scalar | 0.9197 | 0.9793 | 0.3920 | 39.6246 | 1052 | 1216 | 164 |
| `FDCmid_median` | FDC | median | 0.9363 | 0.9907 | 0.0778 | 10.1806 | 0 | 0 | 0 |
| `dur_low_pulses_year_median` | Pulse Metrics | median | 0.9401 | 0.9791 | 0.3681 | 18.5000 | 83 | 100 | 17 |
| `D1_day_median` | Flow Timing | median | 0.9421 | 0.9876 | 1.0125 | 133.5000 | 0 | 0 | 0 |
| `n_low_pulses_all_mean` | Pulse Metrics | mean | 0.9461 | 0.9702 | 0.2787 | 10.7591 | 0 | 0 | 0 |
| `dur_high_pulses_all_median` | Pulse Metrics | median | 0.9480 | 0.9915 | 0.6304 | 21.6250 | 0 | 0 | 0 |
| `concavity_median` | Recession | median | 0.9547 | 0.9848 | 0.1412 | 5.0101 | 1052 | 1216 | 164 |
| `qp_slope_sd_mean` | Q-P Seasonality | mean | 0.9594 | 0.9673 | 0.6321 | 66.8511 | 1085 | 1085 | 0 |
| `n_low_pulses_year_median` | Pulse Metrics | median | 0.9594 | 0.9745 | 0.2562 | 14.5000 | 0 | 0 | 0 |
| `D5_day_median` | Flow Timing | median | 0.9606 | 0.9865 | 2.0575 | 116.0000 | 0 | 0 | 0 |
| `dur_high_pulses_all_mean` | Pulse Metrics | mean | 0.9619 | 0.9961 | 0.6375 | 35.3472 | 0 | 0 | 0 |
| `FDC90th_mean` | FDC | mean | 0.9640 | 0.9790 | 0.4532 | 14.5611 | 0 | 0 | 0 |
| `b_pointcloud_median` | Recession | median | 0.9644 | 0.9817 | 0.0312 | 2.0134 | 1233 | 1342 | 109 |
| `concavity_mean` | Recession | mean | 0.9646 | 0.9901 | 0.1871 | 5.9918 | 1052 | 1216 | 164 |
| `Dmax_median` | Flow Timing | median | 0.9654 | 0.9804 | 4.3271 | 163.5000 | 0 | 0 | 0 |
| `dur_high_pulses_year_median` | Pulse Metrics | median | 0.9684 | 0.9959 | 0.2160 | 18.0000 | 0 | 0 | 0 |
| `FDCall_median` | FDC | median | 0.9698 | 0.9952 | 0.0763 | 6.6265 | 0 | 0 | 0 |
| `D10_day_median` | Flow Timing | median | 0.9701 | 0.9853 | 2.7554 | 93.0000 | 0 | 0 | 0 |
| `D99_day_median` | Flow Timing | median | 0.9701 | 0.9891 | 1.3036 | 101.5000 | 0 | 0 | 0 |
| `b_events_median` | Recession | median | 0.9710 | 0.9924 | 0.0900 | 3.9277 | 1052 | 1216 | 164 |
| `Flow_Reversals_fall_median` | Pulse Metrics | median | 0.9735 | 0.9859 | 0.7451 | 18.0000 | 0 | 0 | 0 |
| `recession_alpha_point_cloud_linear_reservoir` | Recession | scalar | 0.9743 | 0.9907 | 0.0043 | 0.2971 | 129 | 167 | 38 |
| `b_events_mean` | Recession | mean | 0.9752 | 0.9940 | 0.1205 | 4.7767 | 1052 | 1216 | 164 |
| `b_pointcloud_mean` | Recession | mean | 0.9756 | 0.9863 | 0.0323 | 1.1528 | 1233 | 1342 | 109 |
| `D25_to_D75_median` | Flow Timing | median | 0.9761 | 0.9869 | 2.8880 | 88.0000 | 0 | 0 | 0 |
| `Flow_Reversals_winter_median` | Pulse Metrics | median | 0.9769 | 0.9864 | 0.7108 | 16.5000 | 0 | 0 | 0 |
| `log_a_events_mean` | Recession | mean | 0.9770 | 0.9901 | 0.1959 | 9.2042 | 1052 | 1216 | 164 |
| `log_a_events_median` | Recession | median | 0.9771 | 0.9944 | 0.1301 | 9.9365 | 1052 | 1216 | 164 |
| `qp_bimodality_median` | Q-P Seasonality | median | 0.9781 | 0.9888 | 0.0118 | 0.1988 | 1085 | 1085 | 0 |
| `D1_day_mean` | Flow Timing | mean | 0.9785 | 0.9940 | 0.9735 | 45.5784 | 0 | 0 | 0 |
| `n_recession_events_median` | Recession | median | 0.9788 | 0.9913 | 0.2284 | 2.5000 | 0 | 0 | 0 |
| `Dmax_mean` | Flow Timing | mean | 0.9790 | 0.9860 | 3.9483 | 35.0664 | 0 | 0 | 0 |
| `dur_low_pulses_year_mean` | Pulse Metrics | mean | 0.9793 | 0.9859 | 0.3963 | 8.4316 | 83 | 100 | 17 |
| `FDCmid_mean` | FDC | mean | 0.9794 | 0.9957 | 0.0745 | 4.7324 | 0 | 0 | 0 |
| `Flow_Reversals_summer_median` | Pulse Metrics | median | 0.9796 | 0.9896 | 0.6792 | 15.5000 | 0 | 0 | 0 |
| `D20_day_median` | Flow Timing | median | 0.9799 | 0.9823 | 2.9781 | 79.0000 | 0 | 0 | 0 |
| `n_low_pulses_year_mean` | Pulse Metrics | mean | 0.9801 | 0.9886 | 0.2279 | 4.0621 | 0 | 0 | 0 |
| `log_a_pointcloud_median` | Recession | median | 0.9804 | 0.9943 | 0.0452 | 2.4972 | 1233 | 1342 | 109 |
| `Flow_Reversals_fall_mean` | Pulse Metrics | mean | 0.9809 | 0.9906 | 0.7045 | 9.4646 | 0 | 0 | 0 |
| `D30_day_median` | Flow Timing | median | 0.9813 | 0.9852 | 2.7254 | 145.0000 | 0 | 0 | 0 |
| `Flow_Reversals_spring_median` | Pulse Metrics | median | 0.9813 | 0.9894 | 0.5613 | 16.0000 | 0 | 0 | 0 |
| `qp_slope_sd_median` | Q-P Seasonality | median | 0.9813 | 0.9893 | 0.1336 | 24.5498 | 1085 | 1085 | 0 |
| `D40_day_median` | Flow Timing | median | 0.9814 | 0.9867 | 2.5803 | 97.0000 | 0 | 0 | 0 |
| `BFI_Eckhardt_param_median` | Baseflow | median | 0.9821 | 0.9898 | 0.0021 | 0.1044 | 129 | 167 | 38 |
| `BFI_Eckhardt_param_mean` | Baseflow | mean | 0.9835 | 0.9908 | 0.0021 | 0.0916 | 129 | 167 | 38 |
| `D95_day_median` | Flow Timing | median | 0.9838 | 0.9916 | 2.0863 | 52.0000 | 0 | 0 | 0 |
| `TQmean_median` | Pulse Metrics | median | 0.9839 | 0.9915 | 0.5114 | 24.8245 | 0 | 0 | 0 |
| `D50_day_median` | Flow Timing | median | 0.9845 | 0.9907 | 2.2107 | 109.5000 | 0 | 0 | 0 |
| `Flow_Reversals_annual_median` | Pulse Metrics | median | 0.9847 | 0.9913 | 2.1113 | 57.0000 | 0 | 0 | 0 |
| `D5_day_mean` | Flow Timing | mean | 0.9848 | 0.9945 | 1.6315 | 38.9594 | 0 | 0 | 0 |
| `D90_day_median` | Flow Timing | median | 0.9850 | 0.9914 | 2.1804 | 66.0000 | 0 | 0 | 0 |
| `BFI_LyneHollick_param_median` | Baseflow | median | 0.9852 | 0.9908 | 0.0071 | 0.2367 | 129 | 167 | 38 |
| `Flow_Reversals_winter_mean` | Pulse Metrics | mean | 0.9854 | 0.9915 | 0.6219 | 9.5007 | 0 | 0 | 0 |
| `alpha_linear_median` | Recession | median | 0.9854 | 0.9937 | 0.0037 | 0.1199 | 1234 | 1342 | 108 |
| `D80_day_median` | Flow Timing | median | 0.9857 | 0.9910 | 2.1036 | 101.0000 | 0 | 0 | 0 |
| `Flow_Reversals_summer_mean` | Pulse Metrics | mean | 0.9858 | 0.9933 | 0.6268 | 8.2964 | 0 | 0 | 0 |
| `FDCall_mean` | FDC | mean | 0.9862 | 0.9961 | 0.0704 | 2.9295 | 0 | 0 | 0 |
| `D60_day_median` | Flow Timing | median | 0.9862 | 0.9920 | 2.0635 | 86.5000 | 0 | 0 | 0 |
| `D70_day_median` | Flow Timing | median | 0.9865 | 0.9920 | 2.0144 | 94.0000 | 0 | 0 | 0 |
| `summer_runoff_ratio_median` | Runoff Ratios | median | 0.9866 | 0.9945 | 0.0579 | 65.5529 | 1085 | 1085 | 0 |
| `qp_bimodality_mean` | Q-P Seasonality | mean | 0.9873 | 0.9933 | 0.0091 | 0.0752 | 1085 | 1085 | 0 |
| `Flow_Reversals_spring_mean` | Pulse Metrics | mean | 0.9875 | 0.9933 | 0.5160 | 8.0909 | 0 | 0 | 0 |
| `Flow_Reversals_annual_mean` | Pulse Metrics | mean | 0.9876 | 0.9934 | 2.0323 | 35.7075 | 0 | 0 | 0 |
| `D25_to_D75_mean` | Flow Timing | mean | 0.9876 | 0.9927 | 2.3134 | 34.0380 | 0 | 0 | 0 |
| `BFI_LyneHollick_param_mean` | Baseflow | mean | 0.9878 | 0.9922 | 0.0065 | 0.2179 | 129 | 167 | 38 |
| `D10_day_mean` | Flow Timing | mean | 0.9882 | 0.9938 | 2.0097 | 27.5154 | 0 | 0 | 0 |
| `log_a_pointcloud_mean` | Recession | mean | 0.9884 | 0.9954 | 0.0456 | 1.6429 | 1233 | 1342 | 109 |
| `D99_day_mean` | Flow Timing | mean | 0.9887 | 0.9947 | 1.0550 | 30.2521 | 0 | 0 | 0 |
| `alpha_linear_mean` | Recession | mean | 0.9890 | 0.9947 | 0.0036 | 0.0854 | 1234 | 1342 | 108 |
| `n_high_pulses_all_median` | Pulse Metrics | median | 0.9890 | 0.9917 | 0.2881 | 7.5000 | 0 | 0 | 0 |
| `avg_storage_median` | Storage | median | 0.9893 | 0.9828 | 22.0474 | 18500.4942 | 1085 | 1085 | 0 |

## NA Mismatch Analysis

Columns where the number of NAs differs by >50 gages.

| Column | Category | Baseline NAs | Experiment NAs | Mismatch | R2 |
|--------|----------|--------------|----------------|----------|----|
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | 4299 | 4080 | 751 | -0.1004 |
| `dur_low_pulses_all_spearman_rho` | Pulse Metrics | 4299 | 4080 | 751 | 0.5471 |
| `dur_low_pulses_all_mk_rho` | Pulse Metrics | 4299 | 4080 | 751 | 0.5613 |
| `dur_low_pulses_all_senn_slp` | Pulse Metrics | 4299 | 4080 | 751 | -1.8187 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | 4299 | 4080 | 751 | -1.3108 |
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | 4299 | 4080 | 751 | -0.0876 |
| `dur_high_pulses_all_linear_slp` | Pulse Metrics | 1443 | 1316 | 665 | 0.6339 |
| `dur_high_pulses_all_mk_pval` | Pulse Metrics | 1443 | 1316 | 665 | -0.0867 |
| `dur_high_pulses_all_senn_slp` | Pulse Metrics | 1443 | 1316 | 665 | 0.5415 |
| `dur_high_pulses_all_mk_rho` | Pulse Metrics | 1443 | 1316 | 665 | 0.5743 |
| `dur_high_pulses_all_spearman_pval` | Pulse Metrics | 1443 | 1316 | 665 | -0.0930 |
| `dur_high_pulses_all_spearman_rho` | Pulse Metrics | 1443 | 1316 | 665 | 0.5686 |
| `dur_low_pulses_year_spearman_rho` | Pulse Metrics | 1146 | 1024 | 542 | 0.5383 |
| `dur_low_pulses_year_spearman_pval` | Pulse Metrics | 1146 | 1024 | 542 | -0.0878 |
| `dur_low_pulses_year_senn_slp` | Pulse Metrics | 1146 | 1024 | 542 | 0.8878 |
| `dur_low_pulses_year_mk_rho` | Pulse Metrics | 1146 | 1024 | 542 | 0.5543 |
| `dur_low_pulses_year_mk_pval` | Pulse Metrics | 1146 | 1024 | 542 | -0.0773 |
| `dur_low_pulses_year_linear_slp` | Pulse Metrics | 1146 | 1024 | 542 | 0.8104 |
| `BFI_LyneHollick_param_mk_rho` | Baseflow | 725 | 667 | 520 | 0.6170 |
| `BFI_LyneHollick_param_senn_slp` | Baseflow | 725 | 667 | 520 | 0.7064 |
| `BFI_LyneHollick_param_spearman_pval` | Baseflow | 725 | 667 | 520 | -0.0242 |
| `BFI_LyneHollick_param_spearman_rho` | Baseflow | 725 | 667 | 520 | 0.6056 |
| `BFI_LyneHollick_param_mk_pval` | Baseflow | 725 | 667 | 520 | -0.0276 |
| `BFI_LyneHollick_param_linear_slp` | Baseflow | 725 | 667 | 520 | 0.7078 |
| `BFI_Eckhardt_param_spearman_rho` | Baseflow | 725 | 667 | 520 | 0.5827 |
| `BFI_Eckhardt_param_spearman_pval` | Baseflow | 725 | 667 | 520 | -0.0021 |
| `BFI_Eckhardt_param_senn_slp` | Baseflow | 725 | 667 | 520 | 0.7417 |
| `BFI_Eckhardt_param_mk_rho` | Baseflow | 725 | 667 | 520 | 0.5916 |
| `BFI_Eckhardt_param_mk_pval` | Baseflow | 725 | 667 | 520 | -0.0016 |
| `BFI_Eckhardt_param_linear_slp` | Baseflow | 725 | 667 | 520 | 0.7248 |

## Summary

| Agreement Tier | Threshold | Count | % |
|----------------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 34 | 5.5% |
| Good | 0.99 <= R2 < 0.999 | 42 | 6.8% |
| Poor | 0.95 <= R2 < 0.99 | 65 | 10.5% |
| Low | 0.90 <= R2 < 0.95 | 9 | 1.5% |
| Very Low | 0.50 <= R2 < 0.90 | 236 | 38.1% |
| Extremely Low | R2 < 0.50 | 233 | 37.6% |
| **Total compared** | | **619** | **100%** |

Gages dropped by experiment filter: **734**

---
*Generated by `docs/benchmarks/compare_experiment_vs_julia.py startIn1993_60pct`*

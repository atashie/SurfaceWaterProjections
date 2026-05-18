# Experiment 'startIn1993_80pct' vs Julia Baseline: Comparison Report

**Generated**: 2026-04-24 14:04:30

## Experiment Description

Water years >= 1993 AND at least 80% of possible years (WY1993 to gage max) must have qualifying data.

## Input Files

| Dataset | File | Gages | Columns |
|---------|------|-------|---------|
| Julia Baseline | `julia_signatures.csv` | 7,313 | 656 |
| Experiment (startIn1993_80pct) | `startIn1993_80pct_signatures.csv` | 5,431 | 656 |

**Common gages**: 5,431
**Dropped gages** (in baseline, not in experiment): 1,882
**Added gages** (in experiment, not in baseline): 0

### Gage Differences

```

  Dropped (baseline only): 1882 gages
    USGS: 1553
    Canada: 329
    non-reference: 1500
    reference: 382
    Examples: 01017550, 01017960, 01021200, 01021480, 01027200
```

### Years Per Gage (Common Gages)

```
Years per gage (common gages):
    Baseline:   mean=42.3, median=45.0
    Experiment: mean=31.6, median=33.0
    Mean diff:  -10.7 years
```

## High-Level Alignment Summary

### Distribution Statistics

| Metric | Value |
|--------|-------|
| Columns compared | 619 |
| Mean R2 (identity) | 0.340258 |
| Median R2 | 0.543115 |
| SD of R2 | 1.894589 |
| Min R2 | -45.0300 |

### Agreement Tiers

| Tier | Threshold | Count | % |
|------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 36 | 5.8% |
| Good | 0.99 <= R2 < 0.999 | 37 | 6.0% |
| Poor | 0.95 <= R2 < 0.99 | 67 | 10.8% |
| Low | 0.90 <= R2 < 0.95 | 9 | 1.5% |
| Very Low | 0.50 <= R2 < 0.90 | 199 | 32.1% |
| Extremely Low | R2 < 0.50 | 271 | 43.8% |
| **Total** | | **619** | **100%** |

## Agreement by Signature Category

| Category | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| Baseflow | 32 | 0 | 4 | 4 | 0 | 16 | 8 | 0.502880 | -0.1866 |
| Elasticity | 19 | 1 | 0 | 0 | 1 | 5 | 12 | -2.059938 | -45.0300 |
| FDC | 24 | 0 | 0 | 4 | 1 | 12 | 7 | 0.483990 | -0.2254 |
| Flashiness | 8 | 0 | 2 | 0 | 0 | 4 | 2 | 0.518245 | -0.1508 |
| Flow Percentiles | 128 | 23 | 9 | 0 | 0 | 24 | 72 | 0.271554 | -1.2371 |
| Flow Timing | 120 | 0 | 9 | 20 | 1 | 28 | 62 | 0.438452 | -0.3644 |
| Flow Volumes | 40 | 7 | 3 | 0 | 0 | 9 | 21 | 0.282973 | -1.6012 |
| Negative Flow | 8 | 0 | 1 | 1 | 0 | 1 | 5 | 0.208237 | -0.7043 |
| Pulse Metrics | 112 | 0 | 4 | 18 | 3 | 50 | 37 | 0.498576 | -0.3036 |
| Q-P Seasonality | 16 | 0 | 0 | 4 | 0 | 5 | 7 | 0.428559 | -0.2806 |
| Recession | 63 | 0 | 0 | 14 | 3 | 32 | 14 | 0.567051 | -0.1412 |
| Runoff Ratios | 41 | 5 | 4 | 1 | 0 | 12 | 19 | 0.381676 | -2.2016 |
| Storage | 8 | 0 | 1 | 1 | 0 | 1 | 5 | 0.434561 | -0.2652 |

## Agreement by Statistic Type

| Stat Type | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|-----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| mean | 76 | 21 | 22 | 29 | 2 | 2 | 0 | 0.984157 | 0.7805 |
| median | 76 | 14 | 15 | 37 | 5 | 5 | 0 | 0.974502 | 0.8112 |
| senn_slp | 76 | 0 | 0 | 0 | 0 | 39 | 37 | 0.306612 | -1.2371 |
| linear_slp | 76 | 0 | 0 | 0 | 0 | 45 | 31 | 0.274201 | -2.2016 |
| spearman_rho | 76 | 0 | 0 | 0 | 0 | 48 | 28 | 0.520534 | 0.2075 |
| spearman_pval | 76 | 0 | 0 | 0 | 0 | 0 | 76 | -0.173348 | -0.7493 |
| mk_rho | 76 | 0 | 0 | 0 | 0 | 54 | 22 | 0.534734 | 0.2385 |
| mk_pval | 76 | 0 | 0 | 0 | 0 | 0 | 76 | -0.166610 | -0.7608 |
| scalar | 11 | 1 | 0 | 1 | 2 | 6 | 1 | -3.340323 | -45.0300 |

## Columns with R2 < 0.99 (546 columns)

| Column | Category | Stat | R2 | Spearman | MAD | Max Diff | Baseline NAs | Experiment NAs | NA Mismatch |
|--------|----------|------|----|----------|-----|----------|--------------|----------------|-------------|
| `elasticity_years_total` | Elasticity | scalar | -45.0300 | 0.7011 | 10.3781 | 12.0000 | 290 | 290 | 0 |
| `spring_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | -2.2016 | 0.7471 | 0.0047 | 4.0210 | 509 | 458 | 243 |
| `Qfal_linear_slp` | Flow Volumes | linear_slp | -1.6012 | 0.6220 | 1.6375 | 3175.9690 | 255 | 192 | 259 |
| `Q10_senn_slp` | Flow Percentiles | senn_slp | -1.2371 | 0.7368 | 0.0133 | 32.7018 | 255 | 192 | 259 |
| `Q5_senn_slp` | Flow Percentiles | senn_slp | -1.2303 | 0.7350 | 0.0131 | 28.2927 | 255 | 192 | 259 |
| `Q1_senn_slp` | Flow Percentiles | senn_slp | -1.2167 | 0.7389 | 0.0117 | 27.0480 | 255 | 192 | 259 |
| `Q10_linear_slp` | Flow Percentiles | linear_slp | -1.1993 | 0.7537 | 0.0117 | 25.4106 | 255 | 192 | 259 |
| `Q20_linear_slp` | Flow Percentiles | linear_slp | -1.1587 | 0.7689 | 0.0121 | 27.6171 | 255 | 192 | 259 |
| `Q25_senn_slp` | Flow Percentiles | senn_slp | -1.0547 | 0.7512 | 0.0159 | 43.0240 | 255 | 192 | 259 |
| `Q20_senn_slp` | Flow Percentiles | senn_slp | -1.0093 | 0.7562 | 0.0153 | 42.1578 | 255 | 192 | 259 |
| `Q30_senn_slp` | Flow Percentiles | senn_slp | -0.8718 | 0.7279 | 0.0172 | 45.2000 | 255 | 192 | 259 |
| `Q25_linear_slp` | Flow Percentiles | linear_slp | -0.8695 | 0.7722 | 0.0126 | 29.3314 | 255 | 192 | 259 |
| `Q30_linear_slp` | Flow Percentiles | linear_slp | -0.7943 | 0.7661 | 0.0135 | 31.5025 | 255 | 192 | 259 |
| `Q5_linear_slp` | Flow Percentiles | linear_slp | -0.7798 | 0.7540 | 0.0113 | 24.7992 | 255 | 192 | 259 |
| `elasticity_rolling_mk_pval` | Elasticity | mk_pval | -0.7608 | 0.1829 | 0.2780 | 1.0000 | 290 | 290 | 0 |
| `elasticity_rolling_spearman_pval` | Elasticity | spearman_pval | -0.7493 | 0.1923 | 0.2652 | 0.9967 | 290 | 290 | 0 |
| `Q1_linear_slp` | Flow Percentiles | linear_slp | -0.7142 | 0.7466 | 0.0099 | 19.2563 | 255 | 192 | 259 |
| `negative_ann_spearman_pval` | Negative Flow | spearman_pval | -0.7043 | 0.3291 | 0.2256 | 0.8662 | 5402 | 5398 | 4 |
| `Qfal_senn_slp` | Flow Volumes | senn_slp | -0.6879 | 0.6055 | 1.7288 | 4084.4444 | 255 | 192 | 259 |
| `negative_ann_mk_pval` | Negative Flow | mk_pval | -0.5650 | 0.3453 | 0.2315 | 0.8502 | 5402 | 5398 | 4 |
| `Q50_linear_slp` | Flow Percentiles | linear_slp | -0.5060 | 0.7822 | 0.0154 | 35.8286 | 255 | 192 | 259 |
| `Q40_linear_slp` | Flow Percentiles | linear_slp | -0.5012 | 0.7729 | 0.0150 | 34.0502 | 255 | 192 | 259 |
| `Q40_senn_slp` | Flow Percentiles | senn_slp | -0.4444 | 0.7377 | 0.0187 | 48.1071 | 255 | 192 | 259 |
| `spring_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | -0.4342 | 0.7312 | 0.0035 | 2.9330 | 509 | 458 | 243 |
| `elasticity_annual_spearman_pval` | Elasticity | spearman_pval | -0.3706 | 0.3200 | 0.2645 | 0.9896 | 290 | 290 | 0 |
| `D99_day_mk_pval` | Flow Timing | mk_pval | -0.3644 | 0.4084 | 0.2594 | 0.9976 | 262 | 200 | 256 |
| `D99_day_spearman_pval` | Flow Timing | spearman_pval | -0.3520 | 0.4130 | 0.2546 | 0.9828 | 262 | 200 | 256 |
| `elasticity_annual_mk_pval` | Elasticity | mk_pval | -0.3468 | 0.3291 | 0.2651 | 0.9625 | 290 | 290 | 0 |
| `D20_day_spearman_pval` | Flow Timing | spearman_pval | -0.3432 | 0.3560 | 0.2423 | 0.9800 | 262 | 200 | 256 |
| `D20_day_mk_pval` | Flow Timing | mk_pval | -0.3341 | 0.3641 | 0.2438 | 0.9628 | 262 | 200 | 256 |
| `summer_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.3259 | 0.3985 | 0.2599 | 0.9958 | 682 | 654 | 244 |
| `Qann_linear_slp` | Flow Volumes | linear_slp | -0.3236 | 0.7488 | 5.9398 | 11771.1424 | 255 | 192 | 259 |
| `Q95_Q10_spearman_pval` | Flow Percentiles | spearman_pval | -0.3109 | 0.4062 | 0.2502 | 0.9854 | 255 | 192 | 259 |
| `summer_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.3109 | 0.4017 | 0.2630 | 0.9988 | 682 | 654 | 244 |
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.3036 | 0.4112 | 0.2511 | 0.9796 | 3452 | 3289 | 675 |
| `Q95_spearman_pval` | Flow Percentiles | spearman_pval | -0.2952 | 0.4154 | 0.2485 | 0.9928 | 255 | 192 | 259 |
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.2930 | 0.4184 | 0.2537 | 0.9792 | 3452 | 3289 | 675 |
| `Q95_Q10_mk_pval` | Flow Percentiles | mk_pval | -0.2808 | 0.4086 | 0.2526 | 0.9835 | 255 | 192 | 259 |
| `qp_bimodality_spearman_pval` | Q-P Seasonality | spearman_pval | -0.2806 | 0.3806 | 0.2537 | 0.9854 | 512 | 467 | 245 |
| `Q90_spearman_pval` | Flow Percentiles | spearman_pval | -0.2786 | 0.4207 | 0.2488 | 0.9856 | 255 | 192 | 259 |
| `D10_day_spearman_pval` | Flow Timing | spearman_pval | -0.2779 | 0.3810 | 0.2405 | 0.9928 | 262 | 200 | 256 |
| `Q50_senn_slp` | Flow Percentiles | senn_slp | -0.2774 | 0.7488 | 0.0175 | 46.9455 | 255 | 192 | 259 |
| `Q95_mk_pval` | Flow Percentiles | mk_pval | -0.2763 | 0.4170 | 0.2509 | 0.9835 | 255 | 192 | 259 |
| `Q60_linear_slp` | Flow Percentiles | linear_slp | -0.2760 | 0.7662 | 0.0160 | 33.0408 | 255 | 192 | 259 |
| `Q99_spearman_pval` | Flow Percentiles | spearman_pval | -0.2740 | 0.4155 | 0.2475 | 0.9889 | 255 | 192 | 259 |
| `Q90_mk_pval` | Flow Percentiles | mk_pval | -0.2738 | 0.4167 | 0.2527 | 0.9939 | 255 | 192 | 259 |
| `Dmax_spearman_pval` | Flow Timing | spearman_pval | -0.2695 | 0.3569 | 0.2432 | 0.9788 | 262 | 200 | 256 |
| `Dmax_mk_pval` | Flow Timing | mk_pval | -0.2693 | 0.3555 | 0.2455 | 0.9727 | 262 | 200 | 256 |
| `avg_storage_spearman_pval` | Storage | spearman_pval | -0.2652 | 0.3884 | 0.2546 | 0.9481 | 524 | 475 | 243 |
| `qp_bimodality_mk_pval` | Q-P Seasonality | mk_pval | -0.2591 | 0.3889 | 0.2544 | 0.9893 | 512 | 467 | 245 |
| `qp_slope_sd_spearman_pval` | Q-P Seasonality | spearman_pval | -0.2591 | 0.4265 | 0.2498 | 0.9614 | 505 | 460 | 245 |
| `D10_day_mk_pval` | Flow Timing | mk_pval | -0.2588 | 0.3920 | 0.2418 | 0.9576 | 262 | 200 | 256 |
| `Q99_mk_pval` | Flow Percentiles | mk_pval | -0.2562 | 0.4191 | 0.2489 | 0.9902 | 255 | 192 | 259 |
| `qp_slope_sd_mk_pval` | Q-P Seasonality | mk_pval | -0.2560 | 0.4276 | 0.2526 | 0.9946 | 505 | 460 | 245 |
| `D5_day_spearman_pval` | Flow Timing | spearman_pval | -0.2526 | 0.3904 | 0.2365 | 0.9631 | 262 | 200 | 256 |
| `D5_day_mk_pval` | Flow Timing | mk_pval | -0.2494 | 0.3902 | 0.2395 | 0.9771 | 262 | 200 | 256 |
| `avg_storage_mk_pval` | Storage | mk_pval | -0.2487 | 0.3914 | 0.2559 | 0.9798 | 524 | 475 | 243 |
| `dur_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.2468 | 0.4060 | 0.2499 | 0.9789 | 1010 | 927 | 451 |
| `dur_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.2449 | 0.4074 | 0.2518 | 0.9968 | 1010 | 927 | 451 |
| `D30_day_mk_pval` | Flow Timing | mk_pval | -0.2356 | 0.4054 | 0.2332 | 0.9521 | 262 | 200 | 256 |
| `spring_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.2340 | 0.4516 | 0.2483 | 0.9842 | 509 | 458 | 243 |
| `Qsum_spearman_pval` | Flow Volumes | spearman_pval | -0.2328 | 0.4623 | 0.2385 | 0.9994 | 255 | 192 | 259 |
| `Qsum_mk_pval` | Flow Volumes | mk_pval | -0.2325 | 0.4598 | 0.2409 | 0.9965 | 256 | 193 | 259 |
| `Qwin_senn_slp` | Flow Volumes | senn_slp | -0.2291 | 0.7244 | 1.6494 | 3283.7500 | 255 | 192 | 259 |
| `D50_day_spearman_pval` | Flow Timing | spearman_pval | -0.2267 | 0.4139 | 0.2340 | 0.9357 | 262 | 200 | 256 |
| `FDC90th_spearman_pval` | FDC | spearman_pval | -0.2254 | 0.4670 | 0.2481 | 0.9944 | 293 | 244 | 265 |
| `D60_day_mk_pval` | Flow Timing | mk_pval | -0.2237 | 0.4177 | 0.2360 | 0.9975 | 262 | 200 | 256 |
| `D60_day_spearman_pval` | Flow Timing | spearman_pval | -0.2233 | 0.4176 | 0.2338 | 0.9651 | 262 | 200 | 256 |
| `spring_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.2216 | 0.4544 | 0.2498 | 0.9931 | 509 | 458 | 243 |
| `FDC90th_mk_pval` | FDC | mk_pval | -0.2204 | 0.4663 | 0.2506 | 0.9972 | 293 | 244 | 265 |
| `D50_day_mk_pval` | Flow Timing | mk_pval | -0.2192 | 0.4171 | 0.2348 | 0.9652 | 262 | 200 | 256 |
| `dur_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.2166 | 0.4474 | 0.2493 | 0.9880 | 709 | 617 | 326 |
| `n_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.2135 | 0.4469 | 0.2487 | 0.9998 | 255 | 192 | 259 |
| `n_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.2127 | 0.4631 | 0.2492 | 0.9920 | 291 | 242 | 265 |
| `n_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.2118 | 0.4642 | 0.2456 | 0.9899 | 291 | 242 | 265 |
| `Qann_senn_slp` | Flow Volumes | senn_slp | -0.2104 | 0.7096 | 6.6875 | 14990.1222 | 255 | 192 | 259 |
| `dur_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.2093 | 0.4516 | 0.2524 | 0.9975 | 709 | 617 | 326 |
| `D30_day_spearman_pval` | Flow Timing | spearman_pval | -0.2087 | 0.4158 | 0.2292 | 0.9597 | 262 | 200 | 256 |
| `fall_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.2046 | 0.4669 | 0.2569 | 0.9994 | 508 | 457 | 245 |
| `D40_day_spearman_pval` | Flow Timing | spearman_pval | -0.1996 | 0.4286 | 0.2310 | 0.9458 | 262 | 200 | 256 |
| `n_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.1991 | 0.4549 | 0.2445 | 0.9700 | 255 | 192 | 259 |
| `D40_day_mk_pval` | Flow Timing | mk_pval | -0.1987 | 0.4235 | 0.2341 | 0.9694 | 262 | 200 | 256 |
| `D25_to_D75_spearman_pval` | Flow Timing | spearman_pval | -0.1934 | 0.4364 | 0.2382 | 0.9799 | 262 | 200 | 256 |
| `fall_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.1919 | 0.4716 | 0.2533 | 0.9962 | 508 | 457 | 245 |
| `BFI_Eckhardt_mk_pval` | Baseflow | mk_pval | -0.1866 | 0.4888 | 0.2470 | 0.9980 | 262 | 200 | 256 |
| `D70_day_spearman_pval` | Flow Timing | spearman_pval | -0.1864 | 0.4343 | 0.2343 | 0.9484 | 262 | 200 | 256 |
| `BFI_Eckhardt_spearman_pval` | Baseflow | spearman_pval | -0.1857 | 0.4920 | 0.2434 | 0.9912 | 262 | 200 | 256 |
| `Qspr_spearman_pval` | Flow Volumes | spearman_pval | -0.1836 | 0.4732 | 0.2392 | 0.9620 | 255 | 192 | 259 |
| `D70_day_mk_pval` | Flow Timing | mk_pval | -0.1812 | 0.4361 | 0.2345 | 0.9679 | 262 | 200 | 256 |
| `winter_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.1809 | 0.5069 | 0.2472 | 0.9945 | 505 | 456 | 241 |
| `D25_to_D75_mk_pval` | Flow Timing | mk_pval | -0.1791 | 0.4437 | 0.2398 | 0.9914 | 262 | 200 | 256 |
| `Qann_mk_pval` | Flow Volumes | mk_pval | -0.1762 | 0.4590 | 0.2457 | 0.9986 | 255 | 192 | 259 |
| `Qann_spearman_pval` | Flow Volumes | spearman_pval | -0.1753 | 0.4625 | 0.2436 | 0.9974 | 255 | 192 | 259 |
| `D1_day_mk_pval` | Flow Timing | mk_pval | -0.1737 | 0.4391 | 0.2313 | 0.9579 | 262 | 200 | 256 |
| `Qspr_mk_pval` | Flow Volumes | mk_pval | -0.1691 | 0.4769 | 0.2398 | 0.9901 | 255 | 192 | 259 |
| `n_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | -0.1687 | 0.5211 | 0.2371 | 0.9973 | 548 | 484 | 280 |
| `D1_day_spearman_pval` | Flow Timing | spearman_pval | -0.1644 | 0.4432 | 0.2270 | 0.9455 | 262 | 200 | 256 |
| `Q80_spearman_pval` | Flow Percentiles | spearman_pval | -0.1628 | 0.4847 | 0.2411 | 0.9956 | 255 | 192 | 259 |
| `Q80_mk_pval` | Flow Percentiles | mk_pval | -0.1613 | 0.4798 | 0.2452 | 0.9978 | 255 | 192 | 259 |
| `BFI_LyneHollick_param_mk_pval` | Baseflow | mk_pval | -0.1577 | 0.5125 | 0.2420 | 0.9825 | 361 | 335 | 278 |
| `D95_day_spearman_pval` | Flow Timing | spearman_pval | -0.1574 | 0.4896 | 0.2332 | 0.9599 | 262 | 200 | 256 |
| `BFI_LyneHollick_param_spearman_pval` | Baseflow | spearman_pval | -0.1551 | 0.5165 | 0.2390 | 0.9862 | 361 | 335 | 278 |
| `D95_day_mk_pval` | Flow Timing | mk_pval | -0.1514 | 0.4962 | 0.2347 | 0.9937 | 262 | 200 | 256 |
| `flashinessRB_mk_pval` | Flashiness | mk_pval | -0.1508 | 0.5195 | 0.2393 | 0.9973 | 262 | 200 | 256 |
| `winter_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.1493 | 0.5189 | 0.2469 | 0.9968 | 505 | 456 | 241 |
| `FDCall_spearman_pval` | FDC | spearman_pval | -0.1490 | 0.5015 | 0.2383 | 0.9920 | 255 | 192 | 259 |
| `flashinessRB_spearman_pval` | Flashiness | spearman_pval | -0.1478 | 0.5242 | 0.2373 | 0.9810 | 262 | 200 | 256 |
| `Q5_mk_pval` | Flow Percentiles | mk_pval | -0.1419 | 0.5560 | 0.2340 | 0.9997 | 301 | 253 | 264 |
| `Q1_mk_pval` | Flow Percentiles | mk_pval | -0.1416 | 0.5557 | 0.2313 | 0.9894 | 316 | 270 | 266 |
| `n_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | -0.1414 | 0.5318 | 0.2311 | 0.9980 | 548 | 484 | 280 |
| `b_events_spearman_pval` | Recession | spearman_pval | -0.1412 | 0.4920 | 0.2372 | 0.9795 | 870 | 1023 | 153 |
| `Q75_mk_pval` | Flow Percentiles | mk_pval | -0.1385 | 0.4979 | 0.2404 | 0.9933 | 255 | 192 | 259 |
| `negative_ann_linear_slp` | Negative Flow | linear_slp | -0.1380 | 0.7239 | 0.0003 | 0.3527 | 255 | 192 | 259 |
| `Q1_spearman_pval` | Flow Percentiles | spearman_pval | -0.1360 | 0.5599 | 0.2281 | 0.9960 | 316 | 270 | 266 |
| `FDCall_mk_pval` | FDC | mk_pval | -0.1356 | 0.5041 | 0.2409 | 0.9972 | 255 | 192 | 259 |
| `b_events_mk_pval` | Recession | mk_pval | -0.1334 | 0.4948 | 0.2400 | 0.9963 | 870 | 1023 | 153 |
| `BFI_LyneHollick_mk_pval` | Baseflow | mk_pval | -0.1331 | 0.5180 | 0.2423 | 0.9867 | 262 | 200 | 256 |
| `Qfal_mk_pval` | Flow Volumes | mk_pval | -0.1319 | 0.5006 | 0.2470 | 0.9902 | 255 | 192 | 259 |
| `Q5_spearman_pval` | Flow Percentiles | spearman_pval | -0.1311 | 0.5640 | 0.2295 | 0.9938 | 301 | 253 | 264 |
| `BFI_LyneHollick_spearman_pval` | Baseflow | spearman_pval | -0.1293 | 0.5212 | 0.2398 | 0.9921 | 262 | 200 | 256 |
| `Q75_spearman_pval` | Flow Percentiles | spearman_pval | -0.1266 | 0.5052 | 0.2366 | 0.9887 | 255 | 192 | 259 |
| `log_a_events_spearman_pval` | Recession | spearman_pval | -0.1248 | 0.4878 | 0.2352 | 0.9910 | 870 | 1023 | 153 |
| `FDCmid_mk_pval` | FDC | mk_pval | -0.1210 | 0.4970 | 0.2363 | 0.9834 | 255 | 192 | 259 |
| `n_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.1189 | 0.4818 | 0.2361 | 0.9686 | 255 | 192 | 259 |
| `BFI_Eckhardt_param_spearman_pval` | Baseflow | spearman_pval | -0.1175 | 0.5082 | 0.2376 | 0.9828 | 361 | 335 | 278 |
| `BFI_Eckhardt_param_mk_pval` | Baseflow | mk_pval | -0.1175 | 0.5068 | 0.2402 | 0.9900 | 361 | 335 | 278 |
| `Qfal_spearman_pval` | Flow Volumes | spearman_pval | -0.1156 | 0.5078 | 0.2435 | 0.9777 | 255 | 192 | 259 |
| `FDCmid_spearman_pval` | FDC | spearman_pval | -0.1152 | 0.4983 | 0.2333 | 0.9838 | 255 | 192 | 259 |
| `dur_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | -0.1132 | 0.4829 | 0.2353 | 0.9939 | 264 | 203 | 255 |
| `log_a_events_mk_pval` | Recession | mk_pval | -0.1103 | 0.4924 | 0.2368 | 0.9867 | 870 | 1023 | 153 |
| `n_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.1087 | 0.4844 | 0.2384 | 0.9822 | 255 | 192 | 259 |
| `dur_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | -0.1083 | 0.4849 | 0.2377 | 0.9820 | 264 | 203 | 255 |
| `Q10_mk_pval` | Flow Percentiles | mk_pval | -0.1075 | 0.5576 | 0.2331 | 0.9999 | 290 | 241 | 265 |
| `concavity_spearman_pval` | Recession | spearman_pval | -0.1054 | 0.5014 | 0.2322 | 0.9879 | 870 | 1023 | 153 |
| `Q10_spearman_pval` | Flow Percentiles | spearman_pval | -0.1021 | 0.5649 | 0.2297 | 0.9866 | 290 | 241 | 265 |
| `b_pointcloud_spearman_pval` | Recession | spearman_pval | -0.0994 | 0.4861 | 0.2251 | 0.9824 | 1031 | 1133 | 102 |
| `annual_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | -0.0951 | 0.5305 | 0.2390 | 0.9970 | 504 | 455 | 241 |
| `Q70_spearman_pval` | Flow Percentiles | spearman_pval | -0.0943 | 0.5206 | 0.2323 | 0.9884 | 256 | 193 | 259 |
| `Flow_Reversals_summer_mk_pval` | Pulse Metrics | mk_pval | -0.0936 | 0.5799 | 0.2267 | 1.0000 | 255 | 192 | 259 |
| `Flow_Reversals_winter_mk_pval` | Pulse Metrics | mk_pval | -0.0936 | 0.5589 | 0.2247 | 0.9996 | 257 | 197 | 262 |
| `concavity_mk_pval` | Recession | mk_pval | -0.0914 | 0.5089 | 0.2330 | 0.9902 | 870 | 1023 | 153 |
| `Q70_mk_pval` | Flow Percentiles | mk_pval | -0.0909 | 0.5184 | 0.2353 | 0.9991 | 256 | 193 | 259 |
| `D80_day_spearman_pval` | Flow Timing | spearman_pval | -0.0879 | 0.4951 | 0.2225 | 0.9692 | 262 | 200 | 256 |
| `Q20_mk_pval` | Flow Percentiles | mk_pval | -0.0867 | 0.5594 | 0.2336 | 0.9972 | 276 | 218 | 262 |
| `TQmean_spearman_pval` | Pulse Metrics | spearman_pval | -0.0858 | 0.5005 | 0.2310 | 0.9880 | 255 | 192 | 259 |
| `annual_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | -0.0854 | 0.5355 | 0.2340 | 0.9957 | 504 | 455 | 241 |
| `Flow_Reversals_summer_spearman_pval` | Pulse Metrics | spearman_pval | -0.0852 | 0.5812 | 0.2233 | 0.9921 | 255 | 192 | 259 |
| `Flow_Reversals_winter_spearman_pval` | Pulse Metrics | spearman_pval | -0.0819 | 0.5637 | 0.2221 | 0.9939 | 257 | 197 | 262 |
| `Qwin_linear_slp` | Flow Volumes | linear_slp | -0.0801 | 0.7709 | 1.6002 | 2572.4661 | 255 | 192 | 259 |
| `TQmean_mk_pval` | Pulse Metrics | mk_pval | -0.0785 | 0.5030 | 0.2343 | 0.9816 | 255 | 192 | 259 |
| `D80_day_mk_pval` | Flow Timing | mk_pval | -0.0745 | 0.5003 | 0.2239 | 0.9805 | 262 | 200 | 256 |
| `Q20_spearman_pval` | Flow Percentiles | spearman_pval | -0.0734 | 0.5670 | 0.2300 | 0.9919 | 276 | 218 | 262 |
| `Qwin_spearman_pval` | Flow Volumes | spearman_pval | -0.0663 | 0.5347 | 0.2341 | 0.9911 | 255 | 192 | 259 |
| `Q60_spearman_pval` | Flow Percentiles | spearman_pval | -0.0625 | 0.5424 | 0.2293 | 0.9981 | 256 | 193 | 259 |
| `Q25_mk_pval` | Flow Percentiles | mk_pval | -0.0619 | 0.5646 | 0.2292 | 0.9997 | 274 | 217 | 263 |
| `Q25_spearman_pval` | Flow Percentiles | spearman_pval | -0.0594 | 0.5686 | 0.2273 | 0.9925 | 274 | 217 | 263 |
| `Q95_Q10_linear_slp` | Flow Percentiles | linear_slp | -0.0577 | 0.6964 | 0.0265 | 20.6690 | 255 | 192 | 259 |
| `Q60_mk_pval` | Flow Percentiles | mk_pval | -0.0546 | 0.5416 | 0.2317 | 0.9905 | 256 | 193 | 259 |
| `log_a_pointcloud_spearman_pval` | Recession | spearman_pval | -0.0483 | 0.5183 | 0.2205 | 0.9506 | 1031 | 1133 | 102 |
| `D90_day_spearman_pval` | Flow Timing | spearman_pval | -0.0466 | 0.5296 | 0.2217 | 0.9696 | 262 | 200 | 256 |
| `Qwin_mk_pval` | Flow Volumes | mk_pval | -0.0461 | 0.5411 | 0.2351 | 0.9874 | 255 | 192 | 259 |
| `Flow_Reversals_spring_mk_pval` | Pulse Metrics | mk_pval | -0.0400 | 0.5671 | 0.2284 | 0.9991 | 256 | 194 | 258 |
| `b_pointcloud_mk_pval` | Recession | mk_pval | -0.0368 | 0.5143 | 0.2269 | 0.9813 | 1031 | 1133 | 102 |
| `D90_day_mk_pval` | Flow Timing | mk_pval | -0.0368 | 0.5320 | 0.2217 | 0.9891 | 262 | 200 | 256 |
| `Flow_Reversals_spring_spearman_pval` | Pulse Metrics | spearman_pval | -0.0294 | 0.5718 | 0.2262 | 0.9675 | 256 | 194 | 258 |
| `Flow_Reversals_fall_mk_pval` | Pulse Metrics | mk_pval | -0.0282 | 0.6361 | 0.2049 | 0.9997 | 255 | 192 | 259 |
| `Q30_spearman_pval` | Flow Percentiles | spearman_pval | -0.0279 | 0.5811 | 0.2215 | 0.9866 | 270 | 209 | 261 |
| `Flow_Reversals_fall_spearman_pval` | Pulse Metrics | spearman_pval | -0.0248 | 0.6369 | 0.2019 | 0.9882 | 255 | 192 | 259 |
| `Q30_mk_pval` | Flow Percentiles | mk_pval | -0.0241 | 0.5793 | 0.2249 | 0.9859 | 270 | 209 | 261 |
| `Qsum_linear_slp` | Flow Volumes | linear_slp | -0.0199 | 0.7348 | 2.1249 | 2872.1548 | 255 | 192 | 259 |
| `Q95_linear_slp` | Flow Percentiles | linear_slp | -0.0194 | 0.6935 | 0.0340 | 32.3488 | 255 | 192 | 259 |
| `Q60_senn_slp` | Flow Percentiles | senn_slp | -0.0054 | 0.7377 | 0.0174 | 40.3571 | 255 | 192 | 259 |
| `Flow_Reversals_annual_mk_pval` | Pulse Metrics | mk_pval | 0.0063 | 0.6604 | 0.1898 | 0.9986 | 255 | 192 | 259 |
| `log_a_pointcloud_mk_pval` | Recession | mk_pval | 0.0085 | 0.5440 | 0.2249 | 0.9847 | 1031 | 1133 | 102 |
| `Flow_Reversals_annual_spearman_pval` | Pulse Metrics | spearman_pval | 0.0174 | 0.6643 | 0.1865 | 0.9956 | 255 | 192 | 259 |
| `Q50_spearman_pval` | Flow Percentiles | spearman_pval | 0.0193 | 0.5843 | 0.2169 | 0.9917 | 261 | 198 | 259 |
| `Q99_linear_slp` | Flow Percentiles | linear_slp | 0.0330 | 0.7201 | 0.0532 | 22.7524 | 255 | 192 | 259 |
| `Q50_mk_pval` | Flow Percentiles | mk_pval | 0.0331 | 0.5850 | 0.2187 | 0.9984 | 261 | 198 | 259 |
| `Q40_spearman_pval` | Flow Percentiles | spearman_pval | 0.0369 | 0.5997 | 0.2164 | 0.9831 | 264 | 202 | 260 |
| `alpha_linear_spearman_pval` | Recession | spearman_pval | 0.0448 | 0.6057 | 0.2071 | 0.9832 | 1032 | 1133 | 101 |
| `Q40_mk_pval` | Flow Percentiles | mk_pval | 0.0450 | 0.5991 | 0.2182 | 0.9996 | 264 | 202 | 260 |
| `n_recession_events_mk_pval` | Recession | mk_pval | 0.0472 | 0.6540 | 0.2041 | 0.9826 | 47 | 56 | 9 |
| `n_recession_events_spearman_pval` | Recession | spearman_pval | 0.0548 | 0.6558 | 0.2006 | 0.9955 | 47 | 56 | 9 |
| `negative_ann_senn_slp` | Negative Flow | senn_slp | 0.0863 | 0.8164 | 0.0000 | 0.1786 | 255 | 192 | 259 |
| `alpha_linear_mk_pval` | Recession | mk_pval | 0.1050 | 0.6265 | 0.2116 | 0.9825 | 1032 | 1133 | 101 |
| `Q90_linear_slp` | Flow Percentiles | linear_slp | 0.1301 | 0.6926 | 0.0282 | 40.2369 | 255 | 192 | 259 |
| `Qsum_senn_slp` | Flow Volumes | senn_slp | 0.1486 | 0.7215 | 2.0353 | 2665.1399 | 255 | 192 | 259 |
| `Q70_linear_slp` | Flow Percentiles | linear_slp | 0.1543 | 0.7597 | 0.0158 | 32.4460 | 255 | 192 | 259 |
| `Q99_senn_slp` | Flow Percentiles | senn_slp | 0.1674 | 0.6936 | 0.0487 | 25.0000 | 255 | 192 | 259 |
| `Q90_senn_slp` | Flow Percentiles | senn_slp | 0.1864 | 0.6777 | 0.0290 | 51.3101 | 255 | 192 | 259 |
| `Qspr_linear_slp` | Flow Volumes | linear_slp | 0.1930 | 0.7243 | 1.7597 | 3150.5525 | 255 | 192 | 259 |
| `Q95_senn_slp` | Flow Percentiles | senn_slp | 0.1959 | 0.6842 | 0.0354 | 41.8500 | 255 | 192 | 259 |
| `Qspr_senn_slp` | Flow Volumes | senn_slp | 0.2044 | 0.7108 | 2.0016 | 4849.7345 | 255 | 192 | 259 |
| `Q80_linear_slp` | Flow Percentiles | linear_slp | 0.2071 | 0.7421 | 0.0199 | 37.0060 | 255 | 192 | 259 |
| `elasticity_rolling_spearman_rho` | Elasticity | spearman_rho | 0.2075 | 0.5553 | 0.3272 | 1.4199 | 290 | 290 | 0 |
| `elasticity_rolling_mk_rho` | Elasticity | mk_rho | 0.2385 | 0.5583 | 0.2390 | 1.1086 | 290 | 290 | 0 |
| `annual_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.2388 | 0.7688 | 0.0021 | 2.0050 | 504 | 455 | 241 |
| `Q75_linear_slp` | Flow Percentiles | linear_slp | 0.2405 | 0.7578 | 0.0171 | 35.2612 | 255 | 192 | 259 |
| `Q95_Q10_senn_slp` | Flow Percentiles | senn_slp | 0.2554 | 0.6907 | 0.0255 | 21.3900 | 255 | 192 | 259 |
| `Q80_senn_slp` | Flow Percentiles | senn_slp | 0.2801 | 0.7165 | 0.0205 | 41.8299 | 255 | 192 | 259 |
| `elasticity_rolling_senn_slp` | Elasticity | senn_slp | 0.3064 | 0.5315 | 0.0304 | 0.3517 | 290 | 290 | 0 |
| `elasticity_rolling_linear_slp` | Elasticity | linear_slp | 0.3383 | 0.5482 | 0.0318 | 0.3057 | 290 | 290 | 0 |
| `fall_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.3405 | 0.7410 | 0.0041 | 3.7436 | 508 | 457 | 245 |
| `D20_day_spearman_rho` | Flow Timing | spearman_rho | 0.3482 | 0.6476 | 0.0963 | 0.5317 | 262 | 200 | 256 |
| `annual_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.3613 | 0.7838 | 0.0018 | 1.2977 | 504 | 455 | 241 |
| `D20_day_mk_rho` | Flow Timing | mk_rho | 0.3703 | 0.6565 | 0.0654 | 0.3671 | 262 | 200 | 256 |
| `elasticity_annual_spearman_rho` | Elasticity | spearman_rho | 0.3711 | 0.6327 | 0.1205 | 0.6638 | 290 | 290 | 0 |
| `Q75_senn_slp` | Flow Percentiles | senn_slp | 0.3764 | 0.7313 | 0.0188 | 41.5000 | 255 | 192 | 259 |
| `Q70_senn_slp` | Flow Percentiles | senn_slp | 0.3779 | 0.7270 | 0.0182 | 42.7495 | 255 | 192 | 259 |
| `elasticity_annual_mk_rho` | Elasticity | mk_rho | 0.3799 | 0.6333 | 0.0832 | 0.4223 | 290 | 290 | 0 |
| `D99_day_spearman_rho` | Flow Timing | spearman_rho | 0.3829 | 0.6942 | 0.1090 | 0.8971 | 262 | 200 | 256 |
| `Dmax_spearman_rho` | Flow Timing | spearman_rho | 0.3868 | 0.6375 | 0.0992 | 0.5226 | 262 | 200 | 256 |
| `dur_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.3881 | 0.6866 | 0.0171 | 0.8182 | 264 | 203 | 255 |
| `Dmax_mk_rho` | Flow Timing | mk_rho | 0.3963 | 0.6388 | 0.0681 | 0.3662 | 262 | 200 | 256 |
| `D99_day_mk_rho` | Flow Timing | mk_rho | 0.3972 | 0.6989 | 0.0773 | 0.5850 | 262 | 200 | 256 |
| `summer_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.4017 | 0.6896 | 0.1197 | 0.6362 | 682 | 654 | 244 |
| `elasticity_annual_senn_slp` | Elasticity | senn_slp | 0.4017 | 0.6071 | 0.0322 | 0.3489 | 290 | 290 | 0 |
| `D50_day_spearman_rho` | Flow Timing | spearman_rho | 0.4088 | 0.6716 | 0.0951 | 0.5642 | 262 | 200 | 256 |
| `D10_day_spearman_rho` | Flow Timing | spearman_rho | 0.4090 | 0.6741 | 0.0946 | 0.6612 | 262 | 200 | 256 |
| `D30_day_spearman_rho` | Flow Timing | spearman_rho | 0.4098 | 0.6853 | 0.0891 | 0.4901 | 262 | 200 | 256 |
| `D40_day_spearman_rho` | Flow Timing | spearman_rho | 0.4114 | 0.6761 | 0.0928 | 0.5760 | 262 | 200 | 256 |
| `D40_day_mk_rho` | Flow Timing | mk_rho | 0.4140 | 0.6713 | 0.0643 | 0.4243 | 262 | 200 | 256 |
| `D50_day_mk_rho` | Flow Timing | mk_rho | 0.4142 | 0.6691 | 0.0654 | 0.3867 | 262 | 200 | 256 |
| `D30_day_mk_rho` | Flow Timing | mk_rho | 0.4144 | 0.6847 | 0.0617 | 0.3708 | 262 | 200 | 256 |
| `summer_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.4168 | 0.6912 | 0.0820 | 0.4238 | 682 | 654 | 244 |
| `avg_storage_senn_slp` | Storage | senn_slp | 0.4182 | 0.6703 | 2.6723 | 1403.6658 | 524 | 475 | 243 |
| `qp_bimodality_spearman_rho` | Q-P Seasonality | spearman_rho | 0.4266 | 0.6861 | 0.1104 | 0.5937 | 512 | 467 | 245 |
| `D10_day_mk_rho` | Flow Timing | mk_rho | 0.4278 | 0.6813 | 0.0650 | 0.4028 | 262 | 200 | 256 |
| `dur_low_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.4334 | 0.6939 | 0.1321 | 0.8850 | 3452 | 3289 | 675 |
| `D5_day_spearman_rho` | Flow Timing | spearman_rho | 0.4337 | 0.6786 | 0.0945 | 0.6196 | 262 | 200 | 256 |
| `FDC90th_senn_slp` | FDC | senn_slp | 0.4372 | 0.6876 | 0.0216 | 1.9703 | 255 | 192 | 259 |
| `fall_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.4381 | 0.6917 | 0.0026 | 2.4167 | 508 | 457 | 245 |
| `D60_day_spearman_rho` | Flow Timing | spearman_rho | 0.4423 | 0.6933 | 0.0945 | 0.5444 | 262 | 200 | 256 |
| `dur_low_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.4426 | 0.6945 | 0.0920 | 0.6852 | 3452 | 3289 | 675 |
| `qp_bimodality_mk_rho` | Q-P Seasonality | mk_rho | 0.4431 | 0.6906 | 0.0760 | 0.3817 | 512 | 467 | 245 |
| `D60_day_senn_slp` | Flow Timing | senn_slp | 0.4453 | 0.6417 | 0.2640 | 2.7841 | 262 | 200 | 256 |
| `D5_day_mk_rho` | Flow Timing | mk_rho | 0.4475 | 0.6800 | 0.0656 | 0.4139 | 262 | 200 | 256 |
| `D60_day_mk_rho` | Flow Timing | mk_rho | 0.4531 | 0.6942 | 0.0652 | 0.3839 | 262 | 200 | 256 |
| `Dmax_senn_slp` | Flow Timing | senn_slp | 0.4551 | 0.6352 | 0.5585 | 8.3080 | 262 | 200 | 256 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.4555 | 0.6933 | 0.0959 | 2.2569 | 3452 | 3289 | 675 |
| `avg_storage_spearman_rho` | Storage | spearman_rho | 0.4583 | 0.6772 | 0.1153 | 0.7562 | 524 | 475 | 243 |
| `Q99_spearman_rho` | Flow Percentiles | spearman_rho | 0.4603 | 0.7290 | 0.1115 | 0.6567 | 255 | 192 | 259 |
| `Q90_spearman_rho` | Flow Percentiles | spearman_rho | 0.4623 | 0.7357 | 0.1148 | 0.7223 | 255 | 192 | 259 |
| `D50_day_senn_slp` | Flow Timing | senn_slp | 0.4661 | 0.6287 | 0.2670 | 3.0440 | 262 | 200 | 256 |
| `Q95_Q10_spearman_rho` | Flow Percentiles | spearman_rho | 0.4666 | 0.7356 | 0.1119 | 0.6229 | 255 | 192 | 259 |
| `Q95_spearman_rho` | Flow Percentiles | spearman_rho | 0.4686 | 0.7356 | 0.1132 | 0.6783 | 255 | 192 | 259 |
| `dur_low_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.4692 | 0.6992 | 0.1249 | 0.7304 | 709 | 617 | 326 |
| `avg_storage_mk_rho` | Storage | mk_rho | 0.4745 | 0.6800 | 0.0796 | 0.4969 | 524 | 475 | 243 |
| `Q99_mk_rho` | Flow Percentiles | mk_rho | 0.4756 | 0.7314 | 0.0768 | 0.4176 | 255 | 192 | 259 |
| `Dmax_linear_slp` | Flow Timing | linear_slp | 0.4768 | 0.6464 | 0.6771 | 9.4977 | 262 | 200 | 256 |
| `D70_day_senn_slp` | Flow Timing | senn_slp | 0.4777 | 0.6830 | 0.2608 | 4.0111 | 262 | 200 | 256 |
| `D40_day_senn_slp` | Flow Timing | senn_slp | 0.4849 | 0.6430 | 0.2866 | 2.7679 | 262 | 200 | 256 |
| `dur_low_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.4853 | 0.7036 | 0.0885 | 0.5009 | 709 | 617 | 326 |
| `dur_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.4857 | 0.6617 | 0.0458 | 1.7500 | 1010 | 927 | 451 |
| `D1_day_spearman_rho` | Flow Timing | spearman_rho | 0.4868 | 0.7036 | 0.0933 | 0.6098 | 262 | 200 | 256 |
| `Qann_spearman_rho` | Flow Volumes | spearman_rho | 0.4871 | 0.7583 | 0.1136 | 0.6248 | 255 | 192 | 259 |
| `D5_day_senn_slp` | Flow Timing | senn_slp | 0.4874 | 0.6448 | 0.2154 | 4.3639 | 262 | 200 | 256 |
| `dur_high_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.4876 | 0.7105 | 0.1139 | 0.5810 | 1010 | 927 | 451 |
| `Q95_Q10_mk_rho` | Flow Percentiles | mk_rho | 0.4903 | 0.7423 | 0.0772 | 0.4035 | 255 | 192 | 259 |
| `D60_day_linear_slp` | Flow Timing | linear_slp | 0.4905 | 0.6709 | 0.2615 | 2.9759 | 262 | 200 | 256 |
| `D50_day_linear_slp` | Flow Timing | linear_slp | 0.4908 | 0.6380 | 0.2779 | 3.3354 | 262 | 200 | 256 |
| `dur_high_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.4936 | 0.7132 | 0.0795 | 0.4106 | 1010 | 927 | 451 |
| `Q90_mk_rho` | Flow Percentiles | mk_rho | 0.4941 | 0.7433 | 0.0796 | 0.5087 | 255 | 192 | 259 |
| `Q95_mk_rho` | Flow Percentiles | mk_rho | 0.4947 | 0.7417 | 0.0783 | 0.4400 | 255 | 192 | 259 |
| `D25_to_D75_spearman_rho` | Flow Timing | spearman_rho | 0.4954 | 0.7329 | 0.1003 | 0.6135 | 262 | 200 | 256 |
| `negative_ann_spearman_rho` | Negative Flow | spearman_rho | 0.4955 | 0.8448 | 0.0884 | 0.3326 | 5402 | 5398 | 4 |
| `D70_day_spearman_rho` | Flow Timing | spearman_rho | 0.4958 | 0.7210 | 0.0934 | 0.5860 | 262 | 200 | 256 |
| `D95_day_spearman_rho` | Flow Timing | spearman_rho | 0.4968 | 0.7378 | 0.0988 | 0.7388 | 262 | 200 | 256 |
| `fall_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.4974 | 0.7328 | 0.1228 | 0.8972 | 508 | 457 | 245 |
| `D25_to_D75_mk_rho` | Flow Timing | mk_rho | 0.4993 | 0.7319 | 0.0693 | 0.4248 | 262 | 200 | 256 |
| `qp_slope_sd_linear_slp` | Q-P Seasonality | linear_slp | 0.4996 | 0.5713 | 0.1190 | 11.6260 | 505 | 460 | 245 |
| `n_high_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5006 | 0.7473 | 0.0822 | 0.5558 | 255 | 192 | 259 |
| `D1_day_mk_rho` | Flow Timing | mk_rho | 0.5019 | 0.7060 | 0.0663 | 0.4372 | 262 | 200 | 256 |
| `D40_day_linear_slp` | Flow Timing | linear_slp | 0.5023 | 0.6386 | 0.2992 | 3.0155 | 262 | 200 | 256 |
| `n_low_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.5033 | 0.7161 | 0.1294 | 0.7950 | 291 | 242 | 265 |
| `qp_bimodality_senn_slp` | Q-P Seasonality | senn_slp | 0.5043 | 0.6979 | 0.0014 | 0.0092 | 512 | 467 | 245 |
| `n_high_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5045 | 0.7512 | 0.1137 | 0.7458 | 255 | 192 | 259 |
| `fall_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5060 | 0.7296 | 0.0848 | 0.5975 | 508 | 457 | 245 |
| `Qann_mk_rho` | Flow Volumes | mk_rho | 0.5073 | 0.7613 | 0.0787 | 0.4039 | 255 | 192 | 259 |
| `dur_low_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.5080 | 0.6995 | 0.0707 | 2.0207 | 3452 | 3289 | 675 |
| `D70_day_mk_rho` | Flow Timing | mk_rho | 0.5106 | 0.7243 | 0.0642 | 0.3716 | 262 | 200 | 256 |
| `FDC90th_spearman_rho` | FDC | spearman_rho | 0.5107 | 0.7271 | 0.1337 | 0.8541 | 293 | 244 | 265 |
| `log_a_seasonality_minimum_first_half` | Recession | scalar | 0.5117 | 0.7806 | 21.7206 | 361.1328 | 870 | 1025 | 155 |
| `D95_day_mk_rho` | Flow Timing | mk_rho | 0.5119 | 0.7435 | 0.0682 | 0.4861 | 262 | 200 | 256 |
| `qp_bimodality_linear_slp` | Q-P Seasonality | linear_slp | 0.5120 | 0.7053 | 0.0014 | 0.0105 | 512 | 467 | 245 |
| `D20_day_senn_slp` | Flow Timing | senn_slp | 0.5133 | 0.6420 | 0.3147 | 3.6327 | 262 | 200 | 256 |
| `D10_day_senn_slp` | Flow Timing | senn_slp | 0.5137 | 0.6569 | 0.2771 | 2.8933 | 262 | 200 | 256 |
| `n_low_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.5144 | 0.7181 | 0.0947 | 0.6395 | 291 | 242 | 265 |
| `qp_slope_sd_spearman_rho` | Q-P Seasonality | spearman_rho | 0.5148 | 0.7202 | 0.1241 | 0.5930 | 505 | 460 | 245 |
| `negative_ann_mk_rho` | Negative Flow | mk_rho | 0.5150 | 0.8695 | 0.0688 | 0.2671 | 5402 | 5398 | 4 |
| `b_events_linear_slp` | Recession | linear_slp | 0.5156 | 0.6940 | 0.0205 | 0.9240 | 870 | 1023 | 153 |
| `BFI_Eckhardt_param_spearman_rho` | Baseflow | spearman_rho | 0.5178 | 0.7518 | 0.1197 | 0.9776 | 361 | 335 | 278 |
| `Qfal_spearman_rho` | Flow Volumes | spearman_rho | 0.5199 | 0.7095 | 0.1223 | 0.8724 | 255 | 192 | 259 |
| `Q80_spearman_rho` | Flow Percentiles | spearman_rho | 0.5205 | 0.7692 | 0.1140 | 0.7460 | 255 | 192 | 259 |
| `dur_high_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.5208 | 0.6843 | 0.0395 | 0.6946 | 264 | 203 | 255 |
| `D70_day_linear_slp` | Flow Timing | linear_slp | 0.5222 | 0.7133 | 0.2594 | 2.6946 | 262 | 200 | 256 |
| `D20_day_linear_slp` | Flow Timing | linear_slp | 0.5222 | 0.6357 | 0.3401 | 3.2940 | 262 | 200 | 256 |
| `D30_day_linear_slp` | Flow Timing | linear_slp | 0.5224 | 0.6498 | 0.3199 | 3.0639 | 262 | 200 | 256 |
| `concavity_linear_slp` | Recession | linear_slp | 0.5228 | 0.7081 | 0.0314 | 0.9138 | 870 | 1023 | 153 |
| `FDC90th_mk_rho` | FDC | mk_rho | 0.5238 | 0.7296 | 0.0931 | 0.5246 | 293 | 244 | 265 |
| `D30_day_senn_slp` | Flow Timing | senn_slp | 0.5243 | 0.6695 | 0.2992 | 2.9658 | 262 | 200 | 256 |
| `BFI_Eckhardt_param_mk_rho` | Baseflow | mk_rho | 0.5267 | 0.7535 | 0.0828 | 0.6833 | 361 | 335 | 278 |
| `qp_slope_sd_mk_rho` | Q-P Seasonality | mk_rho | 0.5267 | 0.7238 | 0.0855 | 0.3987 | 505 | 460 | 245 |
| `Qsum_spearman_rho` | Flow Volumes | spearman_rho | 0.5274 | 0.7507 | 0.1113 | 0.7338 | 255 | 192 | 259 |
| `D80_day_senn_slp` | Flow Timing | senn_slp | 0.5332 | 0.7222 | 0.2604 | 2.7619 | 262 | 200 | 256 |
| `FDC90th_linear_slp` | FDC | linear_slp | 0.5349 | 0.6285 | 0.0714 | 3.1092 | 255 | 192 | 259 |
| `BFI_Eckhardt_spearman_rho` | Baseflow | spearman_rho | 0.5354 | 0.7583 | 0.1219 | 0.9792 | 262 | 200 | 256 |
| `Qsum_mk_rho` | Flow Volumes | mk_rho | 0.5408 | 0.7496 | 0.0766 | 0.4679 | 256 | 193 | 259 |
| `Qfal_mk_rho` | Flow Volumes | mk_rho | 0.5410 | 0.7085 | 0.0848 | 0.5377 | 255 | 192 | 259 |
| `D25_to_D75_senn_slp` | Flow Timing | senn_slp | 0.5431 | 0.7125 | 0.3435 | 4.5714 | 262 | 200 | 256 |
| `D95_day_senn_slp` | Flow Timing | senn_slp | 0.5432 | 0.7044 | 0.2217 | 2.7325 | 262 | 200 | 256 |
| `b_events_spearman_rho` | Recession | spearman_rho | 0.5433 | 0.7424 | 0.1223 | 0.6756 | 870 | 1023 | 153 |
| `concavity_spearman_rho` | Recession | spearman_rho | 0.5436 | 0.7391 | 0.1163 | 0.6158 | 870 | 1023 | 153 |
| `Q80_mk_rho` | Flow Percentiles | mk_rho | 0.5438 | 0.7711 | 0.0795 | 0.5182 | 255 | 192 | 259 |
| `winter_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5449 | 0.7523 | 0.1253 | 0.8049 | 505 | 456 | 241 |
| `BFI_LyneHollick_param_spearman_rho` | Baseflow | spearman_rho | 0.5461 | 0.7686 | 0.1208 | 0.8648 | 361 | 335 | 278 |
| `BFI_LyneHollick_spearman_rho` | Baseflow | spearman_rho | 0.5491 | 0.7690 | 0.1227 | 0.8944 | 262 | 200 | 256 |
| `dur_high_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.5495 | 0.6794 | 0.0625 | 1.8284 | 1010 | 927 | 451 |
| `FDCmid_spearman_rho` | FDC | spearman_rho | 0.5495 | 0.7457 | 0.1135 | 0.7671 | 255 | 192 | 259 |
| `BFI_Eckhardt_mk_rho` | Baseflow | mk_rho | 0.5496 | 0.7609 | 0.0843 | 0.6716 | 262 | 200 | 256 |
| `Q75_spearman_rho` | Flow Percentiles | spearman_rho | 0.5497 | 0.7836 | 0.1153 | 0.7774 | 255 | 192 | 259 |
| `D25_to_D75_linear_slp` | Flow Timing | linear_slp | 0.5514 | 0.7122 | 0.3566 | 4.1164 | 262 | 200 | 256 |
| `D80_day_spearman_rho` | Flow Timing | spearman_rho | 0.5521 | 0.7491 | 0.0912 | 0.6596 | 262 | 200 | 256 |
| `flashinessRB_spearman_rho` | Flashiness | spearman_rho | 0.5541 | 0.7725 | 0.1304 | 0.8695 | 262 | 200 | 256 |
| `D90_day_senn_slp` | Flow Timing | senn_slp | 0.5541 | 0.7279 | 0.2515 | 2.7449 | 262 | 200 | 256 |
| `n_low_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.5546 | 0.7557 | 0.1058 | 0.8463 | 548 | 484 | 280 |
| `n_low_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.5553 | 0.7584 | 0.1429 | 1.1851 | 548 | 484 | 280 |
| `D10_day_linear_slp` | Flow Timing | linear_slp | 0.5563 | 0.6624 | 0.3119 | 3.6144 | 262 | 200 | 256 |
| `D80_day_linear_slp` | Flow Timing | linear_slp | 0.5566 | 0.7274 | 0.2625 | 2.6198 | 262 | 200 | 256 |
| `b_events_mk_rho` | Recession | mk_rho | 0.5568 | 0.7467 | 0.0847 | 0.4705 | 870 | 1023 | 153 |
| `concavity_mk_rho` | Recession | mk_rho | 0.5569 | 0.7436 | 0.0803 | 0.4163 | 870 | 1023 | 153 |
| `FDCall_spearman_rho` | FDC | spearman_rho | 0.5579 | 0.7540 | 0.1228 | 0.8746 | 255 | 192 | 259 |
| `D5_day_linear_slp` | Flow Timing | linear_slp | 0.5581 | 0.6562 | 0.2588 | 3.9052 | 262 | 200 | 256 |
| `BFI_LyneHollick_param_mk_rho` | Baseflow | mk_rho | 0.5581 | 0.7712 | 0.0836 | 0.5597 | 361 | 335 | 278 |
| `log_a_events_spearman_rho` | Recession | spearman_rho | 0.5586 | 0.7492 | 0.1155 | 0.8533 | 870 | 1023 | 153 |
| `winter_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5587 | 0.7557 | 0.0861 | 0.5450 | 505 | 456 | 241 |
| `n_low_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.5603 | 0.6603 | 0.0262 | 0.6667 | 255 | 192 | 259 |
| `spring_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5605 | 0.7667 | 0.1198 | 0.6438 | 509 | 458 | 243 |
| `Q70_spearman_rho` | Flow Percentiles | spearman_rho | 0.5605 | 0.7878 | 0.1168 | 0.8040 | 256 | 193 | 259 |
| `n_recession_events_mk_rho` | Recession | mk_rho | 0.5608 | 0.7848 | 0.0876 | 0.5869 | 47 | 56 | 9 |
| `BFI_LyneHollick_mk_rho` | Baseflow | mk_rho | 0.5608 | 0.7715 | 0.0851 | 0.6150 | 262 | 200 | 256 |
| `FDCmid_mk_rho` | FDC | mk_rho | 0.5621 | 0.7469 | 0.0785 | 0.5080 | 255 | 192 | 259 |
| `n_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.5633 | 0.6854 | 0.0212 | 0.4706 | 255 | 192 | 259 |
| `n_recession_events_spearman_rho` | Recession | spearman_rho | 0.5633 | 0.7890 | 0.1160 | 0.8510 | 47 | 56 | 9 |
| `flashinessRB_mk_rho` | Flashiness | mk_rho | 0.5645 | 0.7751 | 0.0906 | 0.6253 | 262 | 200 | 256 |
| `b_pointcloud_spearman_rho` | Recession | spearman_rho | 0.5654 | 0.7493 | 0.1305 | 1.4455 | 1031 | 1133 | 102 |
| `D90_day_spearman_rho` | Flow Timing | spearman_rho | 0.5661 | 0.7559 | 0.0927 | 0.6579 | 262 | 200 | 256 |
| `D90_day_linear_slp` | Flow Timing | linear_slp | 0.5675 | 0.7273 | 0.2573 | 2.9680 | 262 | 200 | 256 |
| `Qspr_spearman_rho` | Flow Volumes | spearman_rho | 0.5684 | 0.7682 | 0.1068 | 0.6402 | 255 | 192 | 259 |
| `Q75_mk_rho` | Flow Percentiles | mk_rho | 0.5701 | 0.7845 | 0.0804 | 0.5212 | 255 | 192 | 259 |
| `FDCall_mk_rho` | FDC | mk_rho | 0.5705 | 0.7554 | 0.0853 | 0.5447 | 255 | 192 | 259 |
| `log_a_events_mk_rho` | Recession | mk_rho | 0.5708 | 0.7518 | 0.0799 | 0.5748 | 870 | 1023 | 153 |
| `D80_day_mk_rho` | Flow Timing | mk_rho | 0.5711 | 0.7552 | 0.0625 | 0.4318 | 262 | 200 | 256 |
| `dur_high_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.5712 | 0.7617 | 0.1063 | 0.6076 | 264 | 203 | 255 |
| `D95_day_linear_slp` | Flow Timing | linear_slp | 0.5722 | 0.7161 | 0.2309 | 2.7095 | 262 | 200 | 256 |
| `Flow_Reversals_summer_spearman_rho` | Pulse Metrics | spearman_rho | 0.5730 | 0.7819 | 0.1321 | 0.8152 | 255 | 192 | 259 |
| `TQmean_spearman_rho` | Pulse Metrics | spearman_rho | 0.5742 | 0.7578 | 0.1083 | 0.6937 | 255 | 192 | 259 |
| `qp_slope_sd_senn_slp` | Q-P Seasonality | senn_slp | 0.5766 | 0.6381 | 0.0139 | 1.7796 | 505 | 460 | 245 |
| `spring_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.5771 | 0.7723 | 0.0822 | 0.4432 | 509 | 458 | 243 |
| `n_high_pulses_year_spearman_rho` | Pulse Metrics | spearman_rho | 0.5782 | 0.7649 | 0.1069 | 0.6214 | 255 | 192 | 259 |
| `dur_high_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.5782 | 0.7635 | 0.0762 | 0.4136 | 264 | 203 | 255 |
| `Q70_mk_rho` | Flow Percentiles | mk_rho | 0.5791 | 0.7881 | 0.0813 | 0.5478 | 256 | 193 | 259 |
| `D1_day_senn_slp` | Flow Timing | senn_slp | 0.5799 | 0.6178 | 0.0958 | 2.9195 | 262 | 200 | 256 |
| `concavity_senn_slp` | Recession | senn_slp | 0.5809 | 0.7040 | 0.0203 | 0.5672 | 870 | 1023 | 153 |
| `Qwin_spearman_rho` | Flow Volumes | spearman_rho | 0.5809 | 0.7822 | 0.1123 | 0.9064 | 255 | 192 | 259 |
| `D90_day_mk_rho` | Flow Timing | mk_rho | 0.5822 | 0.7619 | 0.0637 | 0.4465 | 262 | 200 | 256 |
| `n_high_pulses_year_mk_rho` | Pulse Metrics | mk_rho | 0.5826 | 0.7663 | 0.0780 | 0.4684 | 255 | 192 | 259 |
| `Flow_Reversals_fall_spearman_rho` | Pulse Metrics | spearman_rho | 0.5826 | 0.7890 | 0.1349 | 0.9492 | 255 | 192 | 259 |
| `Flow_Reversals_spring_spearman_rho` | Pulse Metrics | spearman_rho | 0.5837 | 0.7721 | 0.1240 | 0.7447 | 256 | 194 | 258 |
| `b_events_senn_slp` | Recession | senn_slp | 0.5840 | 0.7103 | 0.0124 | 0.3449 | 870 | 1023 | 153 |
| `TQmean_mk_rho` | Pulse Metrics | mk_rho | 0.5840 | 0.7600 | 0.0748 | 0.4600 | 255 | 192 | 259 |
| `Flow_Reversals_winter_senn_slp` | Pulse Metrics | senn_slp | 0.5846 | 0.7635 | 0.0920 | 1.2102 | 255 | 192 | 259 |
| `Flow_Reversals_summer_mk_rho` | Pulse Metrics | mk_rho | 0.5858 | 0.7853 | 0.0939 | 0.5916 | 255 | 192 | 259 |
| `Q60_spearman_rho` | Flow Percentiles | spearman_rho | 0.5879 | 0.7881 | 0.1179 | 0.7858 | 256 | 193 | 259 |
| `Flow_Reversals_spring_mk_rho` | Pulse Metrics | mk_rho | 0.5898 | 0.7728 | 0.0885 | 0.4926 | 256 | 194 | 258 |
| `Qspr_mk_rho` | Flow Volumes | mk_rho | 0.5898 | 0.7744 | 0.0735 | 0.4148 | 255 | 192 | 259 |
| `log_a_pointcloud_spearman_rho` | Recession | spearman_rho | 0.5933 | 0.7658 | 0.1307 | 1.4818 | 1031 | 1133 | 102 |
| `n_low_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.5937 | 0.7140 | 0.0333 | 0.6891 | 255 | 192 | 259 |
| `annual_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.5943 | 0.7873 | 0.1253 | 0.8608 | 504 | 455 | 241 |
| `Flow_Reversals_fall_mk_rho` | Pulse Metrics | mk_rho | 0.5965 | 0.7920 | 0.0965 | 0.6473 | 255 | 192 | 259 |
| `D1_day_linear_slp` | Flow Timing | linear_slp | 0.5968 | 0.6505 | 0.1436 | 4.4434 | 262 | 200 | 256 |
| `b_pointcloud_mk_rho` | Recession | mk_rho | 0.5971 | 0.7558 | 0.0927 | 1.2364 | 1031 | 1133 | 102 |
| `BFI_LyneHollick_senn_slp` | Baseflow | senn_slp | 0.5999 | 0.7526 | 0.0007 | 0.0158 | 262 | 200 | 256 |
| `BFI_LyneHollick_linear_slp` | Baseflow | linear_slp | 0.6027 | 0.7530 | 0.0007 | 0.0206 | 262 | 200 | 256 |
| `Qwin_mk_rho` | Flow Volumes | mk_rho | 0.6032 | 0.7876 | 0.0775 | 0.5950 | 255 | 192 | 259 |
| `Flow_Reversals_winter_spearman_rho` | Pulse Metrics | spearman_rho | 0.6036 | 0.7888 | 0.1388 | 0.8464 | 257 | 197 | 262 |
| `log_a_events_senn_slp` | Recession | senn_slp | 0.6037 | 0.7120 | 0.0201 | 0.7576 | 870 | 1023 | 153 |
| `Q60_mk_rho` | Flow Percentiles | mk_rho | 0.6052 | 0.7877 | 0.0823 | 0.5332 | 256 | 193 | 259 |
| `Q1_spearman_rho` | Flow Percentiles | spearman_rho | 0.6060 | 0.7821 | 0.1447 | 1.2065 | 316 | 270 | 266 |
| `n_low_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.6061 | 0.7016 | 0.0365 | 1.4148 | 255 | 192 | 259 |
| `Flow_Reversals_winter_linear_slp` | Pulse Metrics | linear_slp | 0.6085 | 0.7793 | 0.0940 | 0.9813 | 255 | 192 | 259 |
| `Flow_Reversals_winter_mk_rho` | Pulse Metrics | mk_rho | 0.6091 | 0.7888 | 0.0996 | 0.5784 | 257 | 197 | 262 |
| `annual_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.6091 | 0.7883 | 0.0876 | 0.5443 | 504 | 455 | 241 |
| `BFI_Eckhardt_senn_slp` | Baseflow | senn_slp | 0.6092 | 0.7429 | 0.0005 | 0.0110 | 262 | 200 | 256 |
| `BFI_Eckhardt_linear_slp` | Baseflow | linear_slp | 0.6100 | 0.7475 | 0.0005 | 0.0148 | 262 | 200 | 256 |
| `n_low_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.6127 | 0.7296 | 0.0521 | 1.8651 | 255 | 192 | 259 |
| `n_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.6141 | 0.7118 | 0.0307 | 0.5516 | 255 | 192 | 259 |
| `log_a_events_linear_slp` | Recession | linear_slp | 0.6162 | 0.7030 | 0.0332 | 1.8091 | 870 | 1023 | 153 |
| `log_a_pointcloud_mk_rho` | Recession | mk_rho | 0.6176 | 0.7710 | 0.0936 | 1.2727 | 1031 | 1133 | 102 |
| `TQmean_senn_slp` | Pulse Metrics | senn_slp | 0.6212 | 0.7365 | 0.0655 | 0.8088 | 255 | 192 | 259 |
| `D99_day_senn_slp` | Flow Timing | senn_slp | 0.6228 | 0.6298 | 0.1133 | 4.0679 | 262 | 200 | 256 |
| `Q5_spearman_rho` | Flow Percentiles | spearman_rho | 0.6230 | 0.7861 | 0.1389 | 1.2552 | 301 | 253 | 264 |
| `elasticity_annual_linear_slp` | Elasticity | linear_slp | 0.6232 | 0.6234 | 0.3129 | 26.9528 | 290 | 290 | 0 |
| `Q1_mk_rho` | Flow Percentiles | mk_rho | 0.6239 | 0.7825 | 0.1017 | 0.7992 | 316 | 270 | 266 |
| `n_high_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.6253 | 0.7782 | 0.0264 | 0.4075 | 255 | 192 | 259 |
| `Q50_spearman_rho` | Flow Percentiles | spearman_rho | 0.6257 | 0.7982 | 0.1159 | 0.8183 | 261 | 198 | 259 |
| `Q10_spearman_rho` | Flow Percentiles | spearman_rho | 0.6286 | 0.7884 | 0.1346 | 1.0311 | 290 | 241 | 265 |
| `Flow_Reversals_annual_spearman_rho` | Pulse Metrics | spearman_rho | 0.6292 | 0.8075 | 0.1459 | 0.9938 | 255 | 192 | 259 |
| `Q30_spearman_rho` | Flow Percentiles | spearman_rho | 0.6321 | 0.7904 | 0.1232 | 0.8754 | 270 | 209 | 261 |
| `D99_day_linear_slp` | Flow Timing | linear_slp | 0.6345 | 0.6582 | 0.1377 | 4.0688 | 262 | 200 | 256 |
| `Q25_spearman_rho` | Flow Percentiles | spearman_rho | 0.6350 | 0.7905 | 0.1254 | 0.8751 | 274 | 217 | 263 |
| `Flow_Reversals_fall_senn_slp` | Pulse Metrics | senn_slp | 0.6368 | 0.7833 | 0.0797 | 1.0217 | 255 | 192 | 259 |
| `Q20_spearman_rho` | Flow Percentiles | spearman_rho | 0.6372 | 0.7916 | 0.1278 | 0.8192 | 276 | 218 | 262 |
| `Q5_mk_rho` | Flow Percentiles | mk_rho | 0.6375 | 0.7836 | 0.0980 | 0.8233 | 301 | 253 | 264 |
| `winter_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.6379 | 0.7380 | 0.0036 | 3.5605 | 505 | 456 | 241 |
| `Q40_spearman_rho` | Flow Percentiles | spearman_rho | 0.6381 | 0.7972 | 0.1181 | 0.8798 | 264 | 202 | 260 |
| `TQmean_linear_slp` | Pulse Metrics | linear_slp | 0.6391 | 0.7525 | 0.0654 | 0.7953 | 255 | 192 | 259 |
| `alpha_linear_spearman_rho` | Recession | spearman_rho | 0.6393 | 0.7938 | 0.1292 | 1.6364 | 1032 | 1133 | 101 |
| `Flow_Reversals_annual_mk_rho` | Pulse Metrics | mk_rho | 0.6393 | 0.8098 | 0.1039 | 0.6879 | 255 | 192 | 259 |
| `Q10_mk_rho` | Flow Percentiles | mk_rho | 0.6403 | 0.7836 | 0.0949 | 0.6869 | 290 | 241 | 265 |
| `FDCall_senn_slp` | FDC | senn_slp | 0.6408 | 0.7125 | 0.0075 | 0.3203 | 255 | 192 | 259 |
| `BFI_LyneHollick_param_linear_slp` | Baseflow | linear_slp | 0.6417 | 0.7562 | 0.0007 | 0.0158 | 361 | 335 | 278 |
| `BFI_LyneHollick_param_senn_slp` | Baseflow | senn_slp | 0.6418 | 0.7558 | 0.0007 | 0.0108 | 361 | 335 | 278 |
| `Q50_mk_rho` | Flow Percentiles | mk_rho | 0.6433 | 0.7985 | 0.0811 | 0.5226 | 261 | 198 | 259 |
| `Q30_mk_rho` | Flow Percentiles | mk_rho | 0.6479 | 0.7878 | 0.0866 | 0.5776 | 270 | 209 | 261 |
| `Q20_mk_rho` | Flow Percentiles | mk_rho | 0.6502 | 0.7893 | 0.0902 | 0.5712 | 276 | 218 | 262 |
| `Q25_mk_rho` | Flow Percentiles | mk_rho | 0.6507 | 0.7895 | 0.0881 | 0.5877 | 274 | 217 | 263 |
| `BFI_Eckhardt_param_linear_slp` | Baseflow | linear_slp | 0.6517 | 0.7215 | 0.0002 | 0.0033 | 361 | 335 | 278 |
| `avg_storage_linear_slp` | Storage | linear_slp | 0.6547 | 0.6916 | 2.6359 | 1573.3937 | 524 | 475 | 243 |
| `Q40_mk_rho` | Flow Percentiles | mk_rho | 0.6555 | 0.7960 | 0.0828 | 0.5644 | 264 | 202 | 260 |
| `flashinessRB_linear_slp` | Flashiness | linear_slp | 0.6563 | 0.7377 | 0.0008 | 0.0318 | 262 | 200 | 256 |
| `Flow_Reversals_fall_linear_slp` | Pulse Metrics | linear_slp | 0.6566 | 0.7888 | 0.0798 | 1.0094 | 255 | 192 | 259 |
| `Flow_Reversals_summer_senn_slp` | Pulse Metrics | senn_slp | 0.6590 | 0.7768 | 0.0730 | 1.3000 | 255 | 192 | 259 |
| `log_a_pointcloud_linear_slp` | Recession | linear_slp | 0.6596 | 0.7053 | 0.0091 | 0.4060 | 1031 | 1133 | 102 |
| `FDCmid_senn_slp` | FDC | senn_slp | 0.6600 | 0.7058 | 0.0068 | 0.4185 | 255 | 192 | 259 |
| `Flow_Reversals_spring_senn_slp` | Pulse Metrics | senn_slp | 0.6622 | 0.7677 | 0.0637 | 0.8083 | 255 | 192 | 259 |
| `n_high_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.6636 | 0.7870 | 0.0352 | 0.6639 | 255 | 192 | 259 |
| `FDCmid_linear_slp` | FDC | linear_slp | 0.6651 | 0.7028 | 0.0104 | 0.5769 | 255 | 192 | 259 |
| `Flow_Reversals_summer_linear_slp` | Pulse Metrics | linear_slp | 0.6652 | 0.7834 | 0.0750 | 1.2650 | 255 | 192 | 259 |
| `alpha_linear_mk_rho` | Recession | mk_rho | 0.6661 | 0.7958 | 0.0932 | 1.4182 | 1032 | 1133 | 101 |
| `BFI_Eckhardt_param_senn_slp` | Baseflow | senn_slp | 0.6708 | 0.7293 | 0.0002 | 0.0032 | 361 | 335 | 278 |
| `Flow_Reversals_spring_linear_slp` | Pulse Metrics | linear_slp | 0.6710 | 0.7757 | 0.0662 | 0.7084 | 255 | 192 | 259 |
| `flashinessRB_senn_slp` | Flashiness | senn_slp | 0.6782 | 0.7452 | 0.0008 | 0.0174 | 262 | 200 | 256 |
| `Flow_Reversals_annual_senn_slp` | Pulse Metrics | senn_slp | 0.6823 | 0.8110 | 0.2313 | 3.4220 | 255 | 192 | 259 |
| `n_recession_events_senn_slp` | Recession | senn_slp | 0.6830 | 0.7862 | 0.0171 | 0.2229 | 0 | 0 | 0 |
| `alpha_linear_linear_slp` | Recession | linear_slp | 0.6884 | 0.7468 | 0.0006 | 0.0121 | 1032 | 1133 | 101 |
| `Flow_Reversals_annual_linear_slp` | Pulse Metrics | linear_slp | 0.6899 | 0.8075 | 0.2355 | 3.3095 | 255 | 192 | 259 |
| `b_pointcloud_senn_slp` | Recession | senn_slp | 0.6914 | 0.7511 | 0.0058 | 0.4522 | 1031 | 1133 | 102 |
| `FDCall_linear_slp` | FDC | linear_slp | 0.6922 | 0.7044 | 0.0092 | 0.2994 | 255 | 192 | 259 |
| `log_a_pointcloud_senn_slp` | Recession | senn_slp | 0.6936 | 0.7515 | 0.0080 | 0.6651 | 1031 | 1133 | 102 |
| `b_pointcloud_linear_slp` | Recession | linear_slp | 0.6937 | 0.7055 | 0.0064 | 0.2029 | 1031 | 1133 | 102 |
| `winter_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.6950 | 0.7195 | 0.0034 | 2.1723 | 505 | 456 | 241 |
| `n_recession_events_linear_slp` | Recession | linear_slp | 0.7107 | 0.8304 | 0.0200 | 0.2544 | 0 | 0 | 0 |
| `summer_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.7137 | 0.6493 | 0.0056 | 4.9540 | 682 | 654 | 244 |
| `log_a_seasonality_amplitude_first_half` | Recession | scalar | 0.7193 | 0.8458 | 1.0360 | 75.8573 | 870 | 1025 | 155 |
| `alpha_linear_senn_slp` | Recession | senn_slp | 0.7242 | 0.7768 | 0.0005 | 0.0133 | 1032 | 1133 | 101 |
| `summer_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.7699 | 0.5723 | 0.0088 | 4.8511 | 682 | 654 | 244 |
| `log_a_seasonality_minimum_last_half` | Recession | scalar | 0.7753 | 0.9219 | 10.8305 | 360.9091 | 874 | 1025 | 151 |
| `elasticity_annual_mean` | Elasticity | mean | 0.7805 | 0.8232 | 1.9303 | 199.3512 | 290 | 290 | 0 |
| `dur_low_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.7934 | 0.6976 | 0.0607 | 1.0872 | 709 | 617 | 326 |
| `log_a_seasonality_minimum_all` | Recession | scalar | 0.8105 | 0.9403 | 7.9284 | 363.3182 | 870 | 1023 | 153 |
| `dur_low_pulses_all_median` | Pulse Metrics | median | 0.8112 | 0.9253 | 1.6185 | 116.3139 | 312 | 303 | 37 |
| `elasticity_static` | Elasticity | scalar | 0.8118 | 0.9001 | 0.2081 | 2.9901 | 290 | 290 | 0 |
| `elasticity_annual_median` | Elasticity | median | 0.8556 | 0.9222 | 0.1581 | 3.1040 | 290 | 290 | 0 |
| `runoff_ratio_high_count` | Runoff Ratios | scalar | 0.8575 | 0.9638 | 0.0372 | 12.0000 | 290 | 290 | 0 |
| `FDC90th_median` | FDC | median | 0.8612 | 0.9751 | 0.2230 | 40.4516 | 0 | 0 | 0 |
| `n_low_pulses_all_median` | Pulse Metrics | median | 0.8691 | 0.9275 | 0.4292 | 11.5000 | 0 | 0 | 0 |
| `dur_low_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.8824 | 0.6673 | 0.0333 | 0.7195 | 709 | 617 | 326 |
| `elasticity_rolling_median` | Elasticity | median | 0.8862 | 0.9393 | 0.1553 | 2.4302 | 290 | 290 | 0 |
| `dur_low_pulses_all_mean` | Pulse Metrics | mean | 0.8993 | 0.9596 | 1.6634 | 85.8383 | 312 | 303 | 37 |
| `log_a_seasonality_amplitude_last_half` | Recession | scalar | 0.9082 | 0.9565 | 0.6323 | 49.1892 | 874 | 1025 | 151 |
| `elasticity_rolling_mean` | Elasticity | mean | 0.9182 | 0.9543 | 0.1520 | 1.9362 | 290 | 290 | 0 |
| `log_a_seasonality_amplitude_all` | Recession | scalar | 0.9223 | 0.9775 | 0.4314 | 39.6246 | 870 | 1023 | 153 |
| `dur_low_pulses_year_median` | Pulse Metrics | median | 0.9293 | 0.9767 | 0.4084 | 18.5000 | 69 | 84 | 15 |
| `FDCmid_median` | FDC | median | 0.9296 | 0.9895 | 0.0855 | 10.1806 | 0 | 0 | 0 |
| `n_low_pulses_all_mean` | Pulse Metrics | mean | 0.9394 | 0.9670 | 0.3096 | 10.7591 | 0 | 0 | 0 |
| `dur_high_pulses_all_median` | Pulse Metrics | median | 0.9450 | 0.9905 | 0.6939 | 21.6250 | 0 | 0 | 0 |
| `D1_day_median` | Flow Timing | median | 0.9487 | 0.9866 | 1.0934 | 133.5000 | 0 | 0 | 0 |
| `concavity_median` | Recession | median | 0.9492 | 0.9826 | 0.1603 | 5.0101 | 870 | 1023 | 153 |
| `n_low_pulses_year_median` | Pulse Metrics | median | 0.9536 | 0.9709 | 0.2868 | 14.5000 | 0 | 0 | 0 |
| `FDC90th_mean` | FDC | mean | 0.9580 | 0.9776 | 0.4987 | 14.5611 | 0 | 0 | 0 |
| `b_pointcloud_median` | Recession | median | 0.9584 | 0.9785 | 0.0358 | 2.0134 | 1031 | 1133 | 102 |
| `qp_slope_sd_mean` | Q-P Seasonality | mean | 0.9598 | 0.9679 | 0.6421 | 66.8511 | 290 | 290 | 0 |
| `D5_day_median` | Flow Timing | median | 0.9602 | 0.9851 | 2.2926 | 104.0000 | 0 | 0 | 0 |
| `Dmax_median` | Flow Timing | median | 0.9605 | 0.9780 | 4.8063 | 163.5000 | 0 | 0 | 0 |
| `dur_high_pulses_all_mean` | Pulse Metrics | mean | 0.9606 | 0.9958 | 0.7011 | 31.9194 | 0 | 0 | 0 |
| `concavity_mean` | Recession | mean | 0.9621 | 0.9892 | 0.2114 | 5.5824 | 870 | 1023 | 153 |
| `dur_high_pulses_year_median` | Pulse Metrics | median | 0.9670 | 0.9951 | 0.2381 | 18.0000 | 0 | 0 | 0 |
| `FDCall_median` | FDC | median | 0.9678 | 0.9948 | 0.0845 | 6.6265 | 0 | 0 | 0 |
| `b_events_median` | Recession | median | 0.9691 | 0.9915 | 0.1018 | 3.9277 | 870 | 1023 | 153 |
| `D10_day_median` | Flow Timing | median | 0.9691 | 0.9843 | 3.0748 | 93.0000 | 0 | 0 | 0 |
| `recession_alpha_point_cloud_linear_reservoir` | Recession | scalar | 0.9698 | 0.9899 | 0.0048 | 0.2971 | 105 | 140 | 35 |
| `Flow_Reversals_fall_median` | Pulse Metrics | median | 0.9709 | 0.9842 | 0.8310 | 14.0000 | 0 | 0 | 0 |
| `b_pointcloud_mean` | Recession | mean | 0.9716 | 0.9840 | 0.0370 | 1.1528 | 1031 | 1133 | 102 |
| `D99_day_median` | Flow Timing | median | 0.9728 | 0.9879 | 1.4115 | 84.0000 | 0 | 0 | 0 |
| `D25_to_D75_median` | Flow Timing | median | 0.9737 | 0.9856 | 3.2067 | 88.0000 | 0 | 0 | 0 |
| `b_events_mean` | Recession | mean | 0.9742 | 0.9936 | 0.1351 | 4.7767 | 870 | 1023 | 153 |
| `Flow_Reversals_winter_median` | Pulse Metrics | median | 0.9749 | 0.9852 | 0.7872 | 15.0000 | 0 | 0 | 0 |
| `log_a_events_median` | Recession | median | 0.9754 | 0.9939 | 0.1479 | 9.9365 | 870 | 1023 | 153 |
| `n_recession_events_median` | Recession | median | 0.9756 | 0.9903 | 0.2542 | 2.5000 | 0 | 0 | 0 |
| `dur_low_pulses_year_mean` | Pulse Metrics | mean | 0.9762 | 0.9845 | 0.4400 | 8.4316 | 69 | 84 | 15 |
| `FDCmid_mean` | FDC | mean | 0.9766 | 0.9951 | 0.0831 | 4.7324 | 0 | 0 | 0 |
| `log_a_events_mean` | Recession | mean | 0.9767 | 0.9899 | 0.2181 | 9.2042 | 870 | 1023 | 153 |
| `n_low_pulses_year_mean` | Pulse Metrics | mean | 0.9772 | 0.9873 | 0.2549 | 4.0621 | 0 | 0 | 0 |
| `BFI_Eckhardt_param_median` | Baseflow | median | 0.9773 | 0.9882 | 0.0023 | 0.1044 | 105 | 140 | 35 |
| `Flow_Reversals_summer_median` | Pulse Metrics | median | 0.9773 | 0.9884 | 0.7538 | 15.5000 | 0 | 0 | 0 |
| `Dmax_mean` | Flow Timing | mean | 0.9774 | 0.9850 | 4.3448 | 35.0664 | 0 | 0 | 0 |
| `log_a_pointcloud_median` | Recession | median | 0.9779 | 0.9938 | 0.0514 | 2.4972 | 1031 | 1133 | 102 |
| `D1_day_mean` | Flow Timing | mean | 0.9785 | 0.9931 | 1.0745 | 45.5784 | 0 | 0 | 0 |
| `Flow_Reversals_fall_mean` | Pulse Metrics | mean | 0.9786 | 0.9896 | 0.7854 | 9.0455 | 0 | 0 | 0 |
| `qp_bimodality_median` | Q-P Seasonality | median | 0.9788 | 0.9890 | 0.0117 | 0.1988 | 290 | 290 | 0 |
| `D20_day_median` | Flow Timing | median | 0.9789 | 0.9822 | 3.2837 | 79.0000 | 0 | 0 | 0 |
| `Flow_Reversals_spring_median` | Pulse Metrics | median | 0.9791 | 0.9882 | 0.6220 | 16.0000 | 0 | 0 | 0 |
| `BFI_Eckhardt_param_mean` | Baseflow | mean | 0.9791 | 0.9895 | 0.0023 | 0.0916 | 105 | 140 | 35 |
| `D30_day_median` | Flow Timing | median | 0.9796 | 0.9846 | 3.0123 | 145.0000 | 0 | 0 | 0 |
| `D40_day_median` | Flow Timing | median | 0.9798 | 0.9861 | 2.8425 | 97.0000 | 0 | 0 | 0 |
| `TQmean_median` | Pulse Metrics | median | 0.9815 | 0.9901 | 0.5720 | 24.8245 | 0 | 0 | 0 |
| `qp_slope_sd_median` | Q-P Seasonality | median | 0.9815 | 0.9891 | 0.1363 | 24.5498 | 290 | 290 | 0 |
| `BFI_LyneHollick_param_median` | Baseflow | median | 0.9825 | 0.9895 | 0.0079 | 0.2367 | 105 | 140 | 35 |
| `D95_day_median` | Flow Timing | median | 0.9828 | 0.9908 | 2.2822 | 45.0000 | 0 | 0 | 0 |
| `D50_day_median` | Flow Timing | median | 0.9830 | 0.9903 | 2.4420 | 109.5000 | 0 | 0 | 0 |
| `D90_day_median` | Flow Timing | median | 0.9832 | 0.9902 | 2.4178 | 66.0000 | 0 | 0 | 0 |
| `negative_ann_median` | Negative Flow | median | 0.9832 | 0.9259 | 0.0013 | 5.0000 | 0 | 0 | 0 |
| `Flow_Reversals_annual_median` | Pulse Metrics | median | 0.9833 | 0.9906 | 2.3407 | 41.0000 | 0 | 0 | 0 |
| `D80_day_median` | Flow Timing | median | 0.9836 | 0.9897 | 2.3532 | 101.0000 | 0 | 0 | 0 |
| `Flow_Reversals_winter_mean` | Pulse Metrics | mean | 0.9836 | 0.9905 | 0.6903 | 9.5007 | 0 | 0 | 0 |
| `alpha_linear_median` | Recession | median | 0.9838 | 0.9931 | 0.0041 | 0.1199 | 1032 | 1133 | 101 |
| `Flow_Reversals_summer_mean` | Pulse Metrics | mean | 0.9840 | 0.9926 | 0.6952 | 8.2964 | 0 | 0 | 0 |
| `D5_day_mean` | Flow Timing | mean | 0.9844 | 0.9940 | 1.8146 | 38.2747 | 0 | 0 | 0 |
| `FDCall_mean` | FDC | mean | 0.9844 | 0.9956 | 0.0780 | 2.9295 | 0 | 0 | 0 |
| `D60_day_median` | Flow Timing | median | 0.9849 | 0.9917 | 2.2719 | 86.5000 | 0 | 0 | 0 |
| `BFI_LyneHollick_param_mean` | Baseflow | mean | 0.9855 | 0.9912 | 0.0071 | 0.2179 | 105 | 140 | 35 |
| `D70_day_median` | Flow Timing | median | 0.9857 | 0.9915 | 2.2222 | 94.0000 | 0 | 0 | 0 |
| `Flow_Reversals_spring_mean` | Pulse Metrics | mean | 0.9860 | 0.9927 | 0.5710 | 8.0909 | 0 | 0 | 0 |
| `Flow_Reversals_annual_mean` | Pulse Metrics | mean | 0.9860 | 0.9927 | 2.2566 | 35.7075 | 0 | 0 | 0 |
| `n_high_pulses_all_median` | Pulse Metrics | median | 0.9862 | 0.9901 | 0.3162 | 7.5000 | 0 | 0 | 0 |
| `D25_to_D75_mean` | Flow Timing | mean | 0.9865 | 0.9920 | 2.5489 | 34.0380 | 0 | 0 | 0 |
| `summer_runoff_ratio_median` | Runoff Ratios | median | 0.9866 | 0.9944 | 0.0601 | 65.5529 | 290 | 290 | 0 |
| `log_a_pointcloud_mean` | Recession | mean | 0.9868 | 0.9948 | 0.0517 | 1.6429 | 1031 | 1133 | 102 |
| `alpha_linear_mean` | Recession | mean | 0.9871 | 0.9940 | 0.0040 | 0.0854 | 1032 | 1133 | 101 |
| `D10_day_mean` | Flow Timing | mean | 0.9876 | 0.9934 | 2.2371 | 25.7745 | 0 | 0 | 0 |
| `qp_bimodality_mean` | Q-P Seasonality | mean | 0.9878 | 0.9936 | 0.0090 | 0.0752 | 290 | 290 | 0 |
| `D99_day_mean` | Flow Timing | mean | 0.9880 | 0.9940 | 1.1652 | 22.0198 | 0 | 0 | 0 |
| `avg_storage_median` | Storage | median | 0.9893 | 0.9828 | 22.3325 | 18500.4942 | 290 | 290 | 0 |
| `TQmean_mean` | Pulse Metrics | mean | 0.9895 | 0.9944 | 0.4800 | 7.6642 | 0 | 0 | 0 |
| `n_recession_events_mean` | Recession | mean | 0.9895 | 0.9974 | 0.2216 | 2.1667 | 0 | 0 | 0 |

## NA Mismatch Analysis

Columns where the number of NAs differs by >50 gages.

| Column | Category | Baseline NAs | Experiment NAs | Mismatch | R2 |
|--------|----------|--------------|----------------|----------|----|
| `dur_low_pulses_all_spearman_rho` | Pulse Metrics | 3452 | 3289 | 675 | 0.4334 |
| `dur_low_pulses_all_senn_slp` | Pulse Metrics | 3452 | 3289 | 675 | 0.5080 |
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | 3452 | 3289 | 675 | -0.2930 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | 3452 | 3289 | 675 | 0.4555 |
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | 3452 | 3289 | 675 | -0.3036 |
| `dur_low_pulses_all_mk_rho` | Pulse Metrics | 3452 | 3289 | 675 | 0.4426 |
| `dur_high_pulses_all_mk_pval` | Pulse Metrics | 1010 | 927 | 451 | -0.2449 |
| `dur_high_pulses_all_mk_rho` | Pulse Metrics | 1010 | 927 | 451 | 0.4936 |
| `dur_high_pulses_all_spearman_rho` | Pulse Metrics | 1010 | 927 | 451 | 0.4876 |
| `dur_high_pulses_all_spearman_pval` | Pulse Metrics | 1010 | 927 | 451 | -0.2468 |
| `dur_high_pulses_all_senn_slp` | Pulse Metrics | 1010 | 927 | 451 | 0.4857 |
| `dur_high_pulses_all_linear_slp` | Pulse Metrics | 1010 | 927 | 451 | 0.5495 |
| `dur_low_pulses_year_mk_pval` | Pulse Metrics | 709 | 617 | 326 | -0.2093 |
| `dur_low_pulses_year_linear_slp` | Pulse Metrics | 709 | 617 | 326 | 0.7934 |
| `dur_low_pulses_year_spearman_pval` | Pulse Metrics | 709 | 617 | 326 | -0.2166 |
| `dur_low_pulses_year_senn_slp` | Pulse Metrics | 709 | 617 | 326 | 0.8824 |
| `dur_low_pulses_year_mk_rho` | Pulse Metrics | 709 | 617 | 326 | 0.4853 |
| `dur_low_pulses_year_spearman_rho` | Pulse Metrics | 709 | 617 | 326 | 0.4692 |
| `n_low_pulses_all_spearman_pval` | Pulse Metrics | 548 | 484 | 280 | -0.1414 |
| `n_low_pulses_all_spearman_rho` | Pulse Metrics | 548 | 484 | 280 | 0.5553 |
| `n_low_pulses_all_mk_rho` | Pulse Metrics | 548 | 484 | 280 | 0.5546 |
| `n_low_pulses_all_mk_pval` | Pulse Metrics | 548 | 484 | 280 | -0.1687 |
| `BFI_LyneHollick_param_mk_pval` | Baseflow | 361 | 335 | 278 | -0.1577 |
| `BFI_LyneHollick_param_senn_slp` | Baseflow | 361 | 335 | 278 | 0.6418 |
| `BFI_LyneHollick_param_mk_rho` | Baseflow | 361 | 335 | 278 | 0.5581 |
| `BFI_LyneHollick_param_spearman_rho` | Baseflow | 361 | 335 | 278 | 0.5461 |
| `BFI_LyneHollick_param_spearman_pval` | Baseflow | 361 | 335 | 278 | -0.1551 |
| `BFI_LyneHollick_param_linear_slp` | Baseflow | 361 | 335 | 278 | 0.6417 |
| `BFI_Eckhardt_param_spearman_rho` | Baseflow | 361 | 335 | 278 | 0.5178 |
| `BFI_Eckhardt_param_spearman_pval` | Baseflow | 361 | 335 | 278 | -0.1175 |

## Summary

| Agreement Tier | Threshold | Count | % |
|----------------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 36 | 5.8% |
| Good | 0.99 <= R2 < 0.999 | 37 | 6.0% |
| Poor | 0.95 <= R2 < 0.99 | 67 | 10.8% |
| Low | 0.90 <= R2 < 0.95 | 9 | 1.5% |
| Very Low | 0.50 <= R2 < 0.90 | 199 | 32.1% |
| Extremely Low | R2 < 0.50 | 271 | 43.8% |
| **Total compared** | | **619** | **100%** |

Gages dropped by experiment filter: **1882**

---
*Generated by `docs/benchmarks/compare_experiment_vs_julia.py startIn1993_80pct`*

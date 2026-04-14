# Julia (Post-Section 3) vs Golden R Reference: Comparison Report

**Generated**: 2026-04-14 08:24:53

## Input Files

| Dataset | File | Gages | Columns | Date |
|---------|------|-------|---------|------|
| Golden R | `golden-outputs/streamflow_signatures_full_10feb2026.csv` | 5,707 | 561 | 2026-03-12 |
| Julia (post-Section 3) | `docs/benchmarks/julia_signatures.csv` | 7,313 | 627 | 2026-04-13 |

**Common gages**: 5,697
**Common signature columns**: 551

## Context

The Golden R output (Feb 2026) predates several code changes:
1. **Section 3 Guidelines changes** (Apr 2026, Julia only): D1/D99 timing, recession algorithm fix, n_recession_events, elasticity rename/annual, runoff ratio flag, elasticity diagnostics
2. **Constant-SD flag-only** (Apr 2026): Years with constant monthly SD are flagged but no longer rejected -- adds ~64 extra qualifying gages
3. **Negative-Q conditional rejection** (Apr 2026): Q<0 rejection now config-driven (default: off)
4. **R canonical alignment fixes** (Apr 2026): Removed leftover min_Q filter, FDC negative Q guard, seasonal_flags passthrough
5. **Pre-existing cross-language divergences**: Recession OLS library diffs, FDC90th floating-point precision, Mann-Kendall implementation diffs

Differences reflect ALL of the above, not just Section 3 changes.

## High-Level Alignment Summary

### Distribution Statistics

| Metric | Value |
|--------|-------|
| Columns compared | 551 |
| Mean R2 (identity) | 0.602827 |
| Median R2 | 0.989446 |
| SD of R2 | 2.909206 |
| Min R2 | -44.2829 |

### Agreement Tiers

| Tier | Threshold | Count | % |
|------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 118 | 21.4% |
| Good | 0.99 <= R2 < 0.999 | 145 | 26.3% |
| Poor | 0.95 <= R2 < 0.99 | 208 | 37.7% |
| Low | 0.90 <= R2 < 0.95 | 15 | 2.7% |
| Very Low | 0.50 <= R2 < 0.90 | 22 | 4.0% |
| Extremely Low | R2 < 0.50 | 43 | 7.8% |
| **Total** | | **551** | **100%** |

## Agreement by Signature Category

| Category | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Median R2 | SD R2 | Min R2 |
|----------|------|---------|------|------|-----|----------|---------------|---------|-----------|-------|--------|
| Baseflow | 16 | 1 | 7 | 8 | 0 | 0 | 0 | 0.983964 | 0.984054 | 0.0127 | 0.9637 |
| Elasticity | 9 | 0 | 0 | 4 | 3 | 2 | 0 | 0.916606 | 0.947508 | 0.0750 | 0.7862 |
| FDC | 24 | 0 | 3 | 15 | 4 | 2 | 0 | 0.955863 | 0.973746 | 0.0700 | 0.6472 |
| Flashiness | 8 | 0 | 4 | 3 | 1 | 0 | 0 | 0.977528 | 0.985211 | 0.0252 | 0.9252 |
| Flow Percentiles | 128 | 64 | 24 | 40 | 0 | 0 | 0 | 0.993419 | 0.998950 | 0.0097 | 0.9632 |
| Flow Timing | 104 | 6 | 32 | 66 | 0 | 0 | 0 | 0.983255 | 0.987936 | 0.0121 | 0.9543 |
| Flow Volumes | 40 | 20 | 10 | 10 | 0 | 0 | 0 | 0.991136 | 0.997008 | 0.0116 | 0.9678 |
| Pulse Metrics | 112 | 10 | 55 | 47 | 0 | 0 | 0 | 0.988682 | 0.991347 | 0.0098 | 0.9636 |
| Q-P Seasonality | 16 | 3 | 7 | 6 | 0 | 0 | 0 | 0.987560 | 0.991307 | 0.0129 | 0.9644 |
| Recession | 46 | 0 | 0 | 0 | 0 | 4 | 42 | -3.501952 | -0.783162 | 9.1982 | -44.2829 |
| Runoff Ratios | 40 | 10 | 1 | 7 | 7 | 14 | 1 | 0.874469 | 0.939413 | 0.1775 | 0.1037 |
| Storage | 8 | 4 | 2 | 2 | 0 | 0 | 0 | 0.990020 | 0.995819 | 0.0138 | 0.9673 |

## Agreement by Statistic Type

| Stat Type | Cols | Perfect | Good | Poor | Low | Very Low | Extremely Low | Mean R2 | Min R2 |
|-----------|------|---------|------|------|-----|----------|---------------|---------|--------|
| mean | 68 | 40 | 20 | 3 | 0 | 2 | 3 | 0.934069 | -1.2615 |
| median | 68 | 32 | 26 | 4 | 1 | 2 | 3 | 0.937554 | -0.9618 |
| senn_slp | 68 | 23 | 9 | 24 | 3 | 4 | 5 | -0.148947 | -44.2829 |
| linear_slp | 68 | 23 | 10 | 23 | 4 | 2 | 6 | 0.100353 | -37.8971 |
| spearman_rho | 68 | 0 | 40 | 20 | 3 | 0 | 5 | 0.774624 | -4.0713 |
| spearman_pval | 68 | 0 | 0 | 57 | 0 | 6 | 5 | 0.808096 | -1.2977 |
| mk_rho | 68 | 0 | 40 | 19 | 4 | 0 | 5 | 0.705435 | -6.3793 |
| mk_pval | 68 | 0 | 0 | 57 | 0 | 6 | 5 | 0.798277 | -1.5738 |
| scalar | 7 | 0 | 0 | 1 | 0 | 0 | 6 | -0.240776 | -0.9329 |

## Columns with R2 < 0.99 (288 columns)

| Column | Category | Stat | R2 | Spearman | MAD | Max Diff | R NAs | Jl NAs | NA Mismatch |
|--------|----------|------|----|----------|-----|----------|-------|--------|-------------|
| `b_pointcloud_senn_slp` | Recession | senn_slp | -44.2829 | 0.1436 | 0.0135 | 0.2666 | 4722 | 1095 | 3627 |
| `b_pointcloud_linear_slp` | Recession | linear_slp | -37.8971 | 0.1689 | 0.0128 | 0.2779 | 4722 | 1095 | 3627 |
| `log_a_pointcloud_senn_slp` | Recession | senn_slp | -25.0780 | 0.2186 | 0.0176 | 0.2846 | 4722 | 1095 | 3627 |
| `log_a_pointcloud_linear_slp` | Recession | linear_slp | -14.6680 | 0.2553 | 0.0161 | 0.3391 | 4722 | 1095 | 3627 |
| `b_pointcloud_mk_rho` | Recession | mk_rho | -6.3793 | 0.1122 | 0.3220 | 1.2756 | 4722 | 1095 | 3627 |
| `log_a_pointcloud_mk_rho` | Recession | mk_rho | -5.7573 | 0.2389 | 0.3210 | 1.2492 | 4722 | 1095 | 3627 |
| `b_pointcloud_spearman_rho` | Recession | spearman_rho | -4.0713 | 0.1181 | 0.4073 | 1.4163 | 4722 | 1095 | 3627 |
| `log_a_pointcloud_spearman_rho` | Recession | spearman_rho | -3.5642 | 0.2357 | 0.4045 | 1.3808 | 4722 | 1095 | 3627 |
| `b_pointcloud_mk_pval` | Recession | mk_pval | -1.5738 | 0.0805 | 0.4102 | 1.0000 | 4722 | 1095 | 3627 |
| `log_a_pointcloud_mk_pval` | Recession | mk_pval | -1.5428 | 0.0841 | 0.4112 | 1.0000 | 4722 | 1095 | 3627 |
| `concavity_mk_rho` | Recession | mk_rho | -1.3240 | 0.1328 | 0.1746 | 0.9142 | 4346 | 927 | 3419 |
| `log_a_pointcloud_spearman_pval` | Recession | spearman_pval | -1.2977 | 0.0432 | 0.3868 | 0.9999 | 4722 | 1095 | 3627 |
| `b_pointcloud_spearman_pval` | Recession | spearman_pval | -1.2825 | 0.0760 | 0.3767 | 1.0000 | 4722 | 1095 | 3627 |
| `concavity_senn_slp` | Recession | senn_slp | -1.2692 | 0.1429 | 0.0149 | 0.2691 | 4346 | 927 | 3419 |
| `concavity_mean` | Recession | mean | -1.2615 | 0.5608 | 1.0593 | 4.6119 | 4346 | 927 | 3419 |
| `concavity_spearman_rho` | Recession | spearman_rho | -1.2111 | 0.1404 | 0.2480 | 1.2658 | 4346 | 927 | 3419 |
| `concavity_mk_pval` | Recession | mk_pval | -0.9884 | 0.0296 | 0.3570 | 0.9999 | 4346 | 927 | 3419 |
| `concavity_median` | Recession | median | -0.9618 | 0.6254 | 0.8425 | 4.1101 | 4346 | 927 | 3419 |
| `log_a_seasonality_minimum_first_half` | Recession | scalar | -0.9329 | 0.2305 | 47.7435 | 338.5619 | 4347 | 928 | 3419 |
| `concavity_spearman_pval` | Recession | spearman_pval | -0.9095 | 0.0342 | 0.3494 | 0.9993 | 4346 | 927 | 3419 |
| `b_events_linear_slp` | Recession | linear_slp | -0.8547 | 0.2221 | 0.0086 | 0.1941 | 4346 | 927 | 3419 |
| `log_a_seasonality_minimum_last_half` | Recession | scalar | -0.8070 | 0.3166 | 47.6658 | 339.2441 | 4346 | 931 | 3415 |
| `b_events_mk_pval` | Recession | mk_pval | -0.7863 | 0.1025 | 0.3337 | 0.9968 | 4346 | 927 | 3419 |
| `b_events_spearman_pval` | Recession | spearman_pval | -0.7800 | 0.0986 | 0.3311 | 0.9965 | 4346 | 927 | 3419 |
| `log_a_events_mk_pval` | Recession | mk_pval | -0.7518 | 0.1417 | 0.3227 | 0.9946 | 4346 | 927 | 3419 |
| `log_a_events_spearman_pval` | Recession | spearman_pval | -0.7299 | 0.1466 | 0.3208 | 0.9847 | 4346 | 927 | 3419 |
| `b_events_senn_slp` | Recession | senn_slp | -0.6577 | 0.2367 | 0.0065 | 0.0669 | 4346 | 927 | 3419 |
| `log_a_seasonality_minimum_all` | Recession | scalar | -0.6437 | 0.3288 | 42.2591 | 357.8427 | 4346 | 927 | 3419 |
| `concavity_linear_slp` | Recession | linear_slp | -0.6347 | 0.1217 | 0.0188 | 0.4581 | 4346 | 927 | 3419 |
| `b_events_mk_rho` | Recession | mk_rho | -0.4420 | 0.2692 | 0.1449 | 0.6521 | 4346 | 927 | 3419 |
| `b_events_spearman_rho` | Recession | spearman_rho | -0.4247 | 0.2682 | 0.2086 | 0.8524 | 4346 | 927 | 3419 |
| `log_a_events_mk_rho` | Recession | mk_rho | -0.3931 | 0.4293 | 0.1294 | 0.7446 | 4346 | 927 | 3419 |
| `log_a_events_spearman_rho` | Recession | spearman_rho | -0.3343 | 0.4393 | 0.1846 | 1.1002 | 4346 | 927 | 3419 |
| `b_events_mean` | Recession | mean | -0.2042 | 0.8112 | 0.5022 | 2.7303 | 4346 | 927 | 3419 |
| `b_events_median` | Recession | median | -0.1932 | 0.8207 | 0.4527 | 2.8010 | 4346 | 927 | 3419 |
| `log_a_seasonality_amplitude_last_half` | Recession | scalar | -0.1169 | 0.3979 | 1.3226 | 26.7376 | 4346 | 931 | 3415 |
| `log_a_seasonality_amplitude_all` | Recession | scalar | -0.1045 | 0.4724 | 1.2080 | 37.5587 | 4346 | 927 | 3419 |
| `log_a_events_senn_slp` | Recession | senn_slp | -0.1009 | 0.4058 | 0.0092 | 0.1022 | 4346 | 927 | 3419 |
| `log_a_seasonality_amplitude_first_half` | Recession | scalar | -0.0586 | 0.4026 | 1.1123 | 75.9254 | 4347 | 928 | 3419 |
| `log_a_events_linear_slp` | Recession | linear_slp | -0.0372 | 0.3483 | 0.0127 | 0.2486 | 4346 | 927 | 3419 |
| `fall_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.1037 | 0.9421 | 0.0013 | 1.1852 | 4 | 556 | 552 |
| `b_pointcloud_median` | Recession | median | 0.4536 | 0.7959 | 0.1328 | 0.9733 | 4722 | 1095 | 3627 |
| `b_pointcloud_mean` | Recession | mean | 0.4701 | 0.8210 | 0.1275 | 1.2068 | 4722 | 1095 | 3627 |
| `annual_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | 0.5970 | 0.8696 | 0.1383 | 0.8954 | 4 | 552 | 548 |
| `annual_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | 0.6038 | 0.8709 | 0.1386 | 0.9161 | 4 | 552 | 548 |
| `winter_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | 0.6328 | 0.8812 | 0.1316 | 0.7913 | 4 | 554 | 550 |
| `winter_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | 0.6331 | 0.8817 | 0.1326 | 0.7781 | 4 | 554 | 550 |
| `FDC90th_senn_slp` | FDC | senn_slp | 0.6472 | 0.9962 | 0.0013 | 2.3134 | 1 | 570 | 569 |
| `fall_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | 0.6865 | 0.8849 | 0.1212 | 0.6858 | 4 | 556 | 552 |
| `fall_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | 0.6884 | 0.8849 | 0.1213 | 0.6748 | 4 | 556 | 552 |
| `spring_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.7385 | 0.9809 | 0.0007 | 0.8684 | 4 | 557 | 553 |
| `log_a_events_mean` | Recession | mean | 0.7568 | 0.9391 | 0.4294 | 6.8595 | 4346 | 927 | 3419 |
| `fall_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.7843 | 0.9482 | 0.0006 | 0.4844 | 4 | 556 | 552 |
| `elasticity_rolling_mk_pval` | Elasticity | mk_pval | 0.7862 | 0.9387 | 0.0612 | 0.9985 | 4 | 4 | 0 |
| `elasticity_rolling_spearman_pval` | Elasticity | spearman_pval | 0.7991 | 0.9419 | 0.0560 | 0.9781 | 4 | 4 | 0 |
| `log_a_events_median` | Recession | median | 0.8235 | 0.9524 | 0.3481 | 4.6172 | 4346 | 927 | 3419 |
| `summer_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | 0.8382 | 0.9359 | 0.0831 | 0.7392 | 4 | 733 | 729 |
| `summer_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | 0.8392 | 0.9372 | 0.0821 | 0.7455 | 4 | 733 | 729 |
| `log_a_pointcloud_median` | Recession | median | 0.8429 | 0.9521 | 0.1983 | 1.7515 | 4722 | 1095 | 3627 |
| `spring_runoff_ratio_spearman_pval` | Runoff Ratios | spearman_pval | 0.8544 | 0.9436 | 0.0784 | 0.8414 | 5 | 557 | 552 |
| `spring_runoff_ratio_mk_pval` | Runoff Ratios | mk_pval | 0.8570 | 0.9444 | 0.0780 | 0.8400 | 4 | 557 | 553 |
| `log_a_pointcloud_mean` | Recession | mean | 0.8718 | 0.9574 | 0.1843 | 1.5427 | 4722 | 1095 | 3627 |
| `spring_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.8781 | 0.9822 | 0.0007 | 0.7146 | 4 | 557 | 553 |
| `annual_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.8789 | 0.9291 | 0.0007 | 0.2762 | 4 | 552 | 548 |
| `FDCmid_senn_slp` | FDC | senn_slp | 0.8918 | 0.9824 | 0.0009 | 0.2669 | 0 | 570 | 570 |
| `elasticity_rolling_linear_slp` | Elasticity | linear_slp | 0.9126 | 0.9691 | 0.0035 | 0.3033 | 4 | 4 | 0 |
| `elasticity_rolling_senn_slp` | Elasticity | senn_slp | 0.9133 | 0.9703 | 0.0036 | 0.2682 | 4 | 4 | 0 |
| `flashinessRB_linear_slp` | Flashiness | linear_slp | 0.9252 | 0.9892 | 0.0001 | 0.0202 | 0 | 575 | 575 |
| `annual_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.9331 | 0.9509 | 0.0005 | 0.2438 | 4 | 552 | 548 |
| `annual_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.9355 | 0.9674 | 0.0563 | 0.4103 | 4 | 552 | 548 |
| `annual_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.9368 | 0.9674 | 0.0392 | 0.2928 | 4 | 552 | 548 |
| `FDCall_senn_slp` | FDC | senn_slp | 0.9371 | 0.9822 | 0.0009 | 0.1483 | 0 | 570 | 570 |
| `winter_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.9375 | 0.9612 | 0.0530 | 0.2702 | 4 | 554 | 550 |
| `winter_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.9376 | 0.9609 | 0.0369 | 0.2226 | 4 | 554 | 550 |
| `FDCmid_linear_slp` | FDC | linear_slp | 0.9407 | 0.9852 | 0.0013 | 0.2549 | 0 | 570 | 570 |
| `fall_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.9412 | 0.9651 | 0.0324 | 0.2193 | 4 | 556 | 552 |
| `fall_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.9414 | 0.9655 | 0.0469 | 0.3283 | 4 | 556 | 552 |
| `FDCall_linear_slp` | FDC | linear_slp | 0.9450 | 0.9797 | 0.0011 | 0.1625 | 0 | 570 | 570 |
| `FDCmid_median` | FDC | median | 0.9451 | 0.9928 | 0.0391 | 10.6113 | 0 | 0 | 0 |
| `elasticity_rolling_mk_rho` | Elasticity | mk_rho | 0.9475 | 0.9742 | 0.0338 | 0.8889 | 4 | 4 | 0 |
| `elasticity_rolling_spearman_rho` | Elasticity | spearman_rho | 0.9507 | 0.9754 | 0.0433 | 1.2228 | 4 | 4 | 0 |
| `D5_day_linear_slp` | Flow Timing | linear_slp | 0.9543 | 0.9882 | 0.0168 | 2.5214 | 0 | 575 | 575 |
| `FDC90th_linear_slp` | FDC | linear_slp | 0.9583 | 0.9893 | 0.0035 | 2.3362 | 1 | 570 | 569 |
| `flashinessRB_senn_slp` | Flashiness | senn_slp | 0.9594 | 0.9925 | 0.0001 | 0.0111 | 0 | 575 | 575 |
| `D5_day_senn_slp` | Flow Timing | senn_slp | 0.9602 | 0.9902 | 0.0125 | 1.7679 | 0 | 575 | 575 |
| `FDC90th_mk_pval` | FDC | mk_pval | 0.9614 | 0.9852 | 0.0167 | 0.9181 | 1 | 598 | 597 |
| `FDC90th_spearman_pval` | FDC | spearman_pval | 0.9629 | 0.9857 | 0.0164 | 0.9602 | 1 | 598 | 597 |
| `Q95_Q10_mk_pval` | Flow Percentiles | mk_pval | 0.9632 | 0.9826 | 0.0144 | 0.9663 | 0 | 570 | 570 |
| `D10_day_linear_slp` | Flow Timing | linear_slp | 0.9634 | 0.9901 | 0.0181 | 2.3807 | 0 | 575 | 575 |
| `TQmean_mk_pval` | Pulse Metrics | mk_pval | 0.9636 | 0.9839 | 0.0147 | 0.9491 | 0 | 570 | 570 |
| `BFI_Eckhardt_linear_slp` | Baseflow | linear_slp | 0.9637 | 0.9899 | 0.0000 | 0.0055 | 0 | 575 | 575 |
| `Q95_Q10_spearman_pval` | Flow Percentiles | spearman_pval | 0.9643 | 0.9831 | 0.0144 | 0.9360 | 0 | 570 | 570 |
| `qp_slope_sd_spearman_pval` | Q-P Seasonality | spearman_pval | 0.9644 | 0.9853 | 0.0178 | 0.8369 | 4 | 553 | 549 |
| `Q95_mk_pval` | Flow Percentiles | mk_pval | 0.9646 | 0.9833 | 0.0146 | 0.9663 | 0 | 570 | 570 |
| `Q95_spearman_pval` | Flow Percentiles | spearman_pval | 0.9650 | 0.9835 | 0.0140 | 0.9360 | 0 | 570 | 570 |
| `qp_slope_sd_mk_pval` | Q-P Seasonality | mk_pval | 0.9650 | 0.9854 | 0.0181 | 0.7923 | 4 | 553 | 549 |
| `D50_day_spearman_pval` | Flow Timing | spearman_pval | 0.9653 | 0.9821 | 0.0150 | 0.7501 | 0 | 575 | 575 |
| `n_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | 0.9657 | 0.9850 | 0.0150 | 0.7409 | 0 | 570 | 570 |
| `n_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | 0.9658 | 0.9852 | 0.0144 | 0.7853 | 1 | 570 | 569 |
| `Q99_mk_pval` | Flow Percentiles | mk_pval | 0.9658 | 0.9841 | 0.0139 | 0.8705 | 0 | 570 | 570 |
| `TQmean_spearman_pval` | Pulse Metrics | spearman_pval | 0.9659 | 0.9848 | 0.0141 | 0.9049 | 0 | 570 | 570 |
| `Dmax_mk_pval` | Flow Timing | mk_pval | 0.9663 | 0.9822 | 0.0147 | 0.6319 | 0 | 575 | 575 |
| `Dmax_spearman_pval` | Flow Timing | spearman_pval | 0.9667 | 0.9824 | 0.0145 | 0.6263 | 0 | 575 | 575 |
| `D50_day_mk_pval` | Flow Timing | mk_pval | 0.9671 | 0.9829 | 0.0150 | 0.8000 | 0 | 575 | 575 |
| `summer_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.9671 | 0.9458 | 0.0019 | 0.7078 | 4 | 733 | 729 |
| `D30_day_spearman_pval` | Flow Timing | spearman_pval | 0.9672 | 0.9829 | 0.0144 | 0.7220 | 0 | 575 | 575 |
| `avg_storage_mk_pval` | Storage | mk_pval | 0.9673 | 0.9839 | 0.0245 | 0.8099 | 4 | 572 | 568 |
| `dur_high_pulses_all_mk_pval` | Pulse Metrics | mk_pval | 0.9674 | 0.9852 | 0.0152 | 0.7829 | 0 | 1316 | 1316 |
| `dur_high_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | 0.9674 | 0.9850 | 0.0149 | 0.7610 | 0 | 1316 | 1316 |
| `Q99_spearman_pval` | Flow Percentiles | spearman_pval | 0.9677 | 0.9849 | 0.0139 | 0.8223 | 0 | 570 | 570 |
| `Qann_mk_pval` | Flow Volumes | mk_pval | 0.9678 | 0.9854 | 0.0133 | 0.9706 | 0 | 570 | 570 |
| `FDCmid_spearman_pval` | FDC | spearman_pval | 0.9678 | 0.9858 | 0.0141 | 0.8355 | 0 | 570 | 570 |
| `D30_day_mk_pval` | Flow Timing | mk_pval | 0.9678 | 0.9832 | 0.0145 | 0.6996 | 0 | 575 | 575 |
| `Qann_spearman_pval` | Flow Volumes | spearman_pval | 0.9679 | 0.9854 | 0.0137 | 0.9642 | 0 | 570 | 570 |
| `n_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.9679 | 0.9864 | 0.0144 | 0.8168 | 0 | 570 | 570 |
| `FDCmid_mk_pval` | FDC | mk_pval | 0.9683 | 0.9860 | 0.0138 | 0.7943 | 0 | 570 | 570 |
| `D60_day_mk_pval` | Flow Timing | mk_pval | 0.9684 | 0.9834 | 0.0151 | 0.6332 | 0 | 575 | 575 |
| `n_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.9684 | 0.9866 | 0.0140 | 0.8559 | 0 | 570 | 570 |
| `D40_day_mk_pval` | Flow Timing | mk_pval | 0.9685 | 0.9834 | 0.0147 | 0.6993 | 0 | 575 | 575 |
| `D20_day_linear_slp` | Flow Timing | linear_slp | 0.9687 | 0.9925 | 0.0181 | 2.3428 | 0 | 575 | 575 |
| `D25_to_D75_spearman_pval` | Flow Timing | spearman_pval | 0.9689 | 0.9853 | 0.0139 | 0.8385 | 0 | 575 | 575 |
| `summer_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.9689 | 0.9822 | 0.0199 | 0.2838 | 4 | 733 | 729 |
| `D60_day_spearman_pval` | Flow Timing | spearman_pval | 0.9691 | 0.9839 | 0.0147 | 0.6279 | 0 | 575 | 575 |
| `D10_day_mk_pval` | Flow Timing | mk_pval | 0.9691 | 0.9841 | 0.0142 | 0.7381 | 0 | 575 | 575 |
| `qp_bimodality_mk_pval` | Q-P Seasonality | mk_pval | 0.9691 | 0.9859 | 0.0234 | 0.5887 | 4 | 560 | 556 |
| `summer_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.9694 | 0.9824 | 0.0285 | 0.3973 | 4 | 733 | 729 |
| `D25_to_D75_mk_pval` | Flow Timing | mk_pval | 0.9694 | 0.9855 | 0.0140 | 0.8283 | 0 | 575 | 575 |
| `avg_storage_spearman_pval` | Storage | spearman_pval | 0.9696 | 0.9849 | 0.0225 | 0.7350 | 4 | 572 | 568 |
| `D10_day_senn_slp` | Flow Timing | senn_slp | 0.9697 | 0.9914 | 0.0148 | 2.2583 | 0 | 575 | 575 |
| `D40_day_spearman_pval` | Flow Timing | spearman_pval | 0.9698 | 0.9844 | 0.0142 | 0.6172 | 0 | 575 | 575 |
| `D10_day_spearman_pval` | Flow Timing | spearman_pval | 0.9698 | 0.9844 | 0.0138 | 0.7855 | 0 | 575 | 575 |
| `D95_day_spearman_pval` | Flow Timing | spearman_pval | 0.9698 | 0.9869 | 0.0129 | 0.9222 | 0 | 575 | 575 |
| `BFI_LyneHollick_linear_slp` | Baseflow | linear_slp | 0.9700 | 0.9911 | 0.0000 | 0.0065 | 0 | 575 | 575 |
| `D20_day_spearman_pval` | Flow Timing | spearman_pval | 0.9702 | 0.9843 | 0.0138 | 0.8159 | 0 | 575 | 575 |
| `D5_day_mk_pval` | Flow Timing | mk_pval | 0.9704 | 0.9850 | 0.0141 | 0.6959 | 0 | 575 | 575 |
| `qp_bimodality_spearman_pval` | Q-P Seasonality | spearman_pval | 0.9705 | 0.9865 | 0.0222 | 0.6502 | 4 | 560 | 556 |
| `BFI_Eckhardt_senn_slp` | Baseflow | senn_slp | 0.9707 | 0.9915 | 0.0000 | 0.0052 | 0 | 575 | 575 |
| `D95_day_mk_pval` | Flow Timing | mk_pval | 0.9707 | 0.9874 | 0.0132 | 0.9252 | 0 | 575 | 575 |
| `D20_day_mk_pval` | Flow Timing | mk_pval | 0.9707 | 0.9846 | 0.0140 | 0.8213 | 0 | 575 | 575 |
| `D95_day_linear_slp` | Flow Timing | linear_slp | 0.9711 | 0.9903 | 0.0124 | 1.5185 | 0 | 575 | 575 |
| `D70_day_mk_pval` | Flow Timing | mk_pval | 0.9712 | 0.9855 | 0.0145 | 0.4835 | 0 | 575 | 575 |
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | 0.9713 | 0.9886 | 0.0121 | 0.6275 | 313 | 3782 | 3469 |
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | 0.9713 | 0.9886 | 0.0121 | 0.5688 | 313 | 3782 | 3469 |
| `Q90_spearman_pval` | Flow Percentiles | spearman_pval | 0.9714 | 0.9867 | 0.0132 | 0.9743 | 0 | 570 | 570 |
| `Q90_mk_pval` | Flow Percentiles | mk_pval | 0.9716 | 0.9866 | 0.0137 | 0.9453 | 0 | 570 | 570 |
| `Qsum_spearman_pval` | Flow Volumes | spearman_pval | 0.9716 | 0.9876 | 0.0131 | 0.9769 | 0 | 570 | 570 |
| `D60_day_linear_slp` | Flow Timing | linear_slp | 0.9717 | 0.9897 | 0.0142 | 1.5701 | 0 | 575 | 575 |
| `D80_day_mk_pval` | Flow Timing | mk_pval | 0.9717 | 0.9862 | 0.0140 | 0.6755 | 0 | 575 | 575 |
| `D90_day_linear_slp` | Flow Timing | linear_slp | 0.9719 | 0.9905 | 0.0134 | 1.8918 | 0 | 575 | 575 |
| `D5_day_spearman_pval` | Flow Timing | spearman_pval | 0.9719 | 0.9857 | 0.0136 | 0.7002 | 0 | 575 | 575 |
| `Qsum_mk_pval` | Flow Volumes | mk_pval | 0.9721 | 0.9878 | 0.0126 | 0.9712 | 0 | 570 | 570 |
| `D80_day_spearman_pval` | Flow Timing | spearman_pval | 0.9721 | 0.9862 | 0.0138 | 0.5754 | 0 | 575 | 575 |
| `Qwin_mk_pval` | Flow Volumes | mk_pval | 0.9722 | 0.9888 | 0.0125 | 0.9255 | 0 | 570 | 570 |
| `spring_runoff_ratio_spearman_rho` | Runoff Ratios | spearman_rho | 0.9723 | 0.9840 | 0.0273 | 0.3097 | 5 | 557 | 552 |
| `D70_day_spearman_pval` | Flow Timing | spearman_pval | 0.9723 | 0.9861 | 0.0140 | 0.5040 | 0 | 575 | 575 |
| `Qwin_spearman_pval` | Flow Volumes | spearman_pval | 0.9725 | 0.9890 | 0.0128 | 0.8764 | 0 | 570 | 570 |
| `D70_day_linear_slp` | Flow Timing | linear_slp | 0.9727 | 0.9902 | 0.0141 | 1.2471 | 0 | 575 | 575 |
| `dur_low_pulses_all_median` | Pulse Metrics | median | 0.9729 | 0.9974 | 0.1548 | 95.9167 | 313 | 315 | 6 |
| `spring_runoff_ratio_mk_rho` | Runoff Ratios | mk_rho | 0.9731 | 0.9845 | 0.0186 | 0.2124 | 4 | 557 | 553 |
| `D80_day_linear_slp` | Flow Timing | linear_slp | 0.9734 | 0.9900 | 0.0138 | 1.7457 | 0 | 575 | 575 |
| `Qfal_spearman_pval` | Flow Volumes | spearman_pval | 0.9735 | 0.9895 | 0.0129 | 0.7329 | 0 | 570 | 570 |
| `BFI_LyneHollick_senn_slp` | Baseflow | senn_slp | 0.9736 | 0.9913 | 0.0000 | 0.0067 | 0 | 575 | 575 |
| `D90_day_spearman_pval` | Flow Timing | spearman_pval | 0.9736 | 0.9879 | 0.0132 | 0.5806 | 0 | 575 | 575 |
| `FDCall_mk_pval` | FDC | mk_pval | 0.9737 | 0.9889 | 0.0122 | 0.9020 | 0 | 570 | 570 |
| `D90_day_mk_pval` | Flow Timing | mk_pval | 0.9738 | 0.9880 | 0.0133 | 0.7103 | 0 | 575 | 575 |
| `FDCall_spearman_pval` | FDC | spearman_pval | 0.9738 | 0.9889 | 0.0125 | 0.9013 | 0 | 570 | 570 |
| `BFI_LyneHollick_mk_pval` | Baseflow | mk_pval | 0.9742 | 0.9906 | 0.0118 | 0.9355 | 0 | 575 | 575 |
| `Qfal_mk_pval` | Flow Volumes | mk_pval | 0.9743 | 0.9897 | 0.0125 | 0.7144 | 0 | 570 | 570 |
| `D30_day_linear_slp` | Flow Timing | linear_slp | 0.9743 | 0.9917 | 0.0166 | 2.0716 | 0 | 575 | 575 |
| `dur_high_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.9744 | 0.9891 | 0.0133 | 0.9161 | 0 | 577 | 577 |
| `BFI_LyneHollick_spearman_pval` | Baseflow | spearman_pval | 0.9744 | 0.9907 | 0.0120 | 0.8620 | 0 | 575 | 575 |
| `D50_day_linear_slp` | Flow Timing | linear_slp | 0.9745 | 0.9906 | 0.0145 | 1.5541 | 0 | 575 | 575 |
| `TQmean_linear_slp` | Pulse Metrics | linear_slp | 0.9746 | 0.9915 | 0.0034 | 0.4481 | 0 | 570 | 570 |
| `Qspr_mk_pval` | Flow Volumes | mk_pval | 0.9748 | 0.9884 | 0.0124 | 0.9109 | 0 | 570 | 570 |
| `BFI_Eckhardt_mk_pval` | Baseflow | mk_pval | 0.9750 | 0.9902 | 0.0119 | 0.9209 | 0 | 575 | 575 |
| `TQmean_senn_slp` | Pulse Metrics | senn_slp | 0.9752 | 0.9915 | 0.0034 | 0.4777 | 0 | 570 | 570 |
| `dur_high_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.9752 | 0.9892 | 0.0128 | 0.9497 | 0 | 577 | 577 |
| `Q80_mk_pval` | Flow Percentiles | mk_pval | 0.9753 | 0.9891 | 0.0126 | 0.9321 | 0 | 570 | 570 |
| `BFI_Eckhardt_spearman_pval` | Baseflow | spearman_pval | 0.9754 | 0.9904 | 0.0122 | 0.7804 | 0 | 575 | 575 |
| `D60_day_senn_slp` | Flow Timing | senn_slp | 0.9755 | 0.9899 | 0.0128 | 1.1173 | 0 | 575 | 575 |
| `Qspr_spearman_pval` | Flow Volumes | spearman_pval | 0.9755 | 0.9887 | 0.0126 | 0.7826 | 1 | 570 | 569 |
| `flashinessRB_mk_pval` | Flashiness | mk_pval | 0.9758 | 0.9920 | 0.0115 | 0.8621 | 0 | 575 | 575 |
| `flashinessRB_spearman_pval` | Flashiness | spearman_pval | 0.9760 | 0.9921 | 0.0117 | 0.9145 | 0 | 575 | 575 |
| `winter_runoff_ratio_linear_slp` | Runoff Ratios | linear_slp | 0.9762 | 0.9292 | 0.0008 | 0.5219 | 4 | 554 | 550 |
| `Q80_spearman_pval` | Flow Percentiles | spearman_pval | 0.9764 | 0.9896 | 0.0120 | 0.8663 | 0 | 570 | 570 |
| `Q75_spearman_pval` | Flow Percentiles | spearman_pval | 0.9766 | 0.9901 | 0.0117 | 0.8085 | 0 | 570 | 570 |
| `Flow_Reversals_winter_mk_pval` | Pulse Metrics | mk_pval | 0.9766 | 0.9932 | 0.0107 | 0.9023 | 0 | 571 | 571 |
| `Q75_mk_pval` | Flow Percentiles | mk_pval | 0.9766 | 0.9901 | 0.0118 | 0.9290 | 0 | 570 | 570 |
| `Flow_Reversals_winter_spearman_pval` | Pulse Metrics | spearman_pval | 0.9770 | 0.9931 | 0.0103 | 0.8974 | 1 | 571 | 570 |
| `D40_day_linear_slp` | Flow Timing | linear_slp | 0.9772 | 0.9901 | 0.0152 | 1.4636 | 0 | 575 | 575 |
| `D20_day_senn_slp` | Flow Timing | senn_slp | 0.9772 | 0.9932 | 0.0156 | 1.8694 | 0 | 575 | 575 |
| `Q70_spearman_pval` | Flow Percentiles | spearman_pval | 0.9773 | 0.9908 | 0.0116 | 0.7641 | 0 | 570 | 570 |
| `D50_day_senn_slp` | Flow Timing | senn_slp | 0.9774 | 0.9902 | 0.0132 | 1.2967 | 0 | 575 | 575 |
| `elasticity_rolling_median` | Elasticity | median | 0.9775 | 0.9931 | 0.0328 | 2.7079 | 4 | 4 | 0 |
| `Flow_Reversals_summer_spearman_pval` | Pulse Metrics | spearman_pval | 0.9777 | 0.9936 | 0.0102 | 0.8910 | 0 | 570 | 570 |
| `n_low_pulses_all_mk_pval` | Pulse Metrics | mk_pval | 0.9779 | 0.9936 | 0.0104 | 0.7869 | 0 | 842 | 842 |
| `elasticity_static` | Elasticity | scalar | 0.9781 | 0.9912 | 0.0373 | 2.3373 | 4 | 4 | 0 |
| `Q70_mk_pval` | Flow Percentiles | mk_pval | 0.9781 | 0.9911 | 0.0119 | 0.6991 | 0 | 570 | 570 |
| `Flow_Reversals_summer_mk_pval` | Pulse Metrics | mk_pval | 0.9782 | 0.9938 | 0.0104 | 0.7807 | 0 | 570 | 570 |
| `D70_day_senn_slp` | Flow Timing | senn_slp | 0.9782 | 0.9893 | 0.0125 | 0.9675 | 0 | 575 | 575 |
| `D80_day_senn_slp` | Flow Timing | senn_slp | 0.9783 | 0.9900 | 0.0124 | 1.3727 | 0 | 575 | 575 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.9784 | 0.9956 | 0.0035 | 1.2720 | 313 | 3782 | 3469 |
| `dur_low_pulses_all_mean` | Pulse Metrics | mean | 0.9789 | 0.9982 | 0.1351 | 94.1875 | 313 | 315 | 6 |
| `Flow_Reversals_annual_mk_pval` | Pulse Metrics | mk_pval | 0.9791 | 0.9955 | 0.0084 | 0.9208 | 0 | 570 | 570 |
| `D90_day_senn_slp` | Flow Timing | senn_slp | 0.9794 | 0.9920 | 0.0122 | 1.5398 | 0 | 575 | 575 |
| `D30_day_senn_slp` | Flow Timing | senn_slp | 0.9799 | 0.9926 | 0.0148 | 1.7215 | 0 | 575 | 575 |
| `D40_day_senn_slp` | Flow Timing | senn_slp | 0.9801 | 0.9915 | 0.0138 | 1.3662 | 0 | 575 | 575 |
| `n_low_pulses_all_spearman_pval` | Pulse Metrics | spearman_pval | 0.9802 | 0.9941 | 0.0099 | 0.7098 | 309 | 842 | 535 |
| `Q60_spearman_pval` | Flow Percentiles | spearman_pval | 0.9802 | 0.9923 | 0.0109 | 0.9129 | 0 | 570 | 570 |
| `Flow_Reversals_annual_spearman_pval` | Pulse Metrics | spearman_pval | 0.9803 | 0.9956 | 0.0080 | 0.8877 | 0 | 570 | 570 |
| `n_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.9807 | 0.9869 | 0.0009 | 0.0909 | 0 | 570 | 570 |
| `Q60_mk_pval` | Flow Percentiles | mk_pval | 0.9807 | 0.9923 | 0.0111 | 0.8011 | 0 | 570 | 570 |
| `FDC90th_median` | FDC | median | 0.9808 | 0.9991 | 0.0279 | 25.8866 | 1 | 0 | 1 |
| `dur_high_pulses_year_senn_slp` | Pulse Metrics | senn_slp | 0.9808 | 0.9878 | 0.0007 | 0.1420 | 0 | 577 | 577 |
| `D95_day_senn_slp` | Flow Timing | senn_slp | 0.9808 | 0.9918 | 0.0109 | 1.1524 | 0 | 575 | 575 |
| `Flow_Reversals_spring_spearman_pval` | Pulse Metrics | spearman_pval | 0.9812 | 0.9931 | 0.0104 | 0.5856 | 2 | 571 | 569 |
| `Flow_Reversals_spring_mk_pval` | Pulse Metrics | mk_pval | 0.9812 | 0.9931 | 0.0105 | 0.6119 | 0 | 571 | 571 |
| `winter_runoff_ratio_senn_slp` | Runoff Ratios | senn_slp | 0.9818 | 0.9393 | 0.0008 | 0.5011 | 4 | 554 | 550 |
| `Flow_Reversals_fall_spearman_pval` | Pulse Metrics | spearman_pval | 0.9820 | 0.9953 | 0.0085 | 0.7107 | 0 | 570 | 570 |
| `Flow_Reversals_fall_mk_pval` | Pulse Metrics | mk_pval | 0.9821 | 0.9953 | 0.0086 | 0.7732 | 0 | 570 | 570 |
| `Q50_spearman_pval` | Flow Percentiles | spearman_pval | 0.9822 | 0.9935 | 0.0104 | 0.5840 | 0 | 570 | 570 |
| `Q50_mk_pval` | Flow Percentiles | mk_pval | 0.9823 | 0.9935 | 0.0108 | 0.5336 | 0 | 570 | 570 |
| `Q40_mk_pval` | Flow Percentiles | mk_pval | 0.9831 | 0.9940 | 0.0103 | 0.6791 | 0 | 572 | 572 |
| `Q40_spearman_pval` | Flow Percentiles | spearman_pval | 0.9834 | 0.9941 | 0.0098 | 0.7098 | 4 | 572 | 568 |
| `D25_to_D75_senn_slp` | Flow Timing | senn_slp | 0.9840 | 0.9927 | 0.0153 | 1.6559 | 0 | 575 | 575 |
| `FDC90th_mean` | FDC | mean | 0.9842 | 0.9996 | 0.0642 | 36.7766 | 1 | 0 | 1 |
| `elasticity_rolling_mean` | Elasticity | mean | 0.9843 | 0.9945 | 0.0302 | 2.3568 | 4 | 4 | 0 |
| `dur_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.9843 | 0.9902 | 0.0016 | 0.4605 | 0 | 1316 | 1316 |
| `D25_to_D75_linear_slp` | Flow Timing | linear_slp | 0.9845 | 0.9929 | 0.0152 | 1.6645 | 0 | 575 | 575 |
| `n_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.9852 | 0.9944 | 0.0098 | 0.4536 | 0 | 596 | 596 |
| `Flow_Reversals_spring_senn_slp` | Pulse Metrics | senn_slp | 0.9853 | 0.9942 | 0.0029 | 0.4367 | 0 | 570 | 570 |
| `Q30_mk_pval` | Flow Percentiles | mk_pval | 0.9853 | 0.9948 | 0.0098 | 0.5586 | 0 | 577 | 577 |
| `n_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.9858 | 0.9945 | 0.0096 | 0.4285 | 34 | 596 | 562 |
| `Q30_spearman_pval` | Flow Percentiles | spearman_pval | 0.9860 | 0.9951 | 0.0092 | 0.6068 | 9 | 577 | 568 |
| `Dmax_linear_slp` | Flow Timing | linear_slp | 0.9861 | 0.9922 | 0.0244 | 1.4346 | 0 | 575 | 575 |
| `FDC90th_mk_rho` | FDC | mk_rho | 0.9861 | 0.9934 | 0.0053 | 0.5394 | 1 | 598 | 597 |
| `Q95_Q10_mk_rho` | Flow Percentiles | mk_rho | 0.9866 | 0.9940 | 0.0039 | 0.3263 | 0 | 570 | 570 |
| `Dmax_senn_slp` | Flow Timing | senn_slp | 0.9868 | 0.9916 | 0.0207 | 1.1586 | 0 | 575 | 575 |
| `Flow_Reversals_spring_linear_slp` | Pulse Metrics | linear_slp | 0.9868 | 0.9965 | 0.0029 | 0.3777 | 0 | 570 | 570 |
| `FDCmid_spearman_rho` | FDC | spearman_rho | 0.9869 | 0.9930 | 0.0059 | 0.5007 | 0 | 570 | 570 |
| `Q20_mk_pval` | Flow Percentiles | mk_pval | 0.9869 | 0.9956 | 0.0088 | 0.6534 | 0 | 582 | 582 |
| `Q95_Q10_spearman_rho` | Flow Percentiles | spearman_rho | 0.9869 | 0.9940 | 0.0054 | 0.4701 | 0 | 570 | 570 |
| `qp_slope_sd_mk_rho` | Q-P Seasonality | mk_rho | 0.9870 | 0.9938 | 0.0051 | 0.3483 | 4 | 553 | 549 |
| `Q25_mk_pval` | Flow Percentiles | mk_pval | 0.9872 | 0.9954 | 0.0090 | 0.6298 | 0 | 580 | 580 |
| `TQmean_mk_rho` | Pulse Metrics | mk_rho | 0.9872 | 0.9940 | 0.0040 | 0.3383 | 0 | 570 | 570 |
| `Q95_mk_rho` | Flow Percentiles | mk_rho | 0.9872 | 0.9942 | 0.0039 | 0.3263 | 0 | 570 | 570 |
| `Q10_mk_pval` | Flow Percentiles | mk_pval | 0.9874 | 0.9962 | 0.0084 | 0.4289 | 0 | 595 | 595 |
| `qp_slope_sd_spearman_rho` | Q-P Seasonality | spearman_rho | 0.9874 | 0.9938 | 0.0070 | 0.4973 | 4 | 553 | 549 |
| `FDCmid_mk_rho` | FDC | mk_rho | 0.9874 | 0.9932 | 0.0041 | 0.3390 | 0 | 570 | 570 |
| `Q25_spearman_pval` | Flow Percentiles | spearman_pval | 0.9874 | 0.9955 | 0.0085 | 0.6729 | 12 | 580 | 568 |
| `TQmean_spearman_rho` | Pulse Metrics | spearman_rho | 0.9875 | 0.9940 | 0.0057 | 0.5096 | 0 | 570 | 570 |
| `Q95_spearman_rho` | Flow Percentiles | spearman_rho | 0.9876 | 0.9942 | 0.0053 | 0.4701 | 0 | 570 | 570 |
| `Q20_spearman_pval` | Flow Percentiles | spearman_pval | 0.9877 | 0.9959 | 0.0083 | 0.5962 | 16 | 582 | 566 |
| `n_high_pulses_all_mk_rho` | Pulse Metrics | mk_rho | 0.9877 | 0.9936 | 0.0045 | 0.2718 | 0 | 570 | 570 |
| `Q1_mk_pval` | Flow Percentiles | mk_pval | 0.9878 | 0.9968 | 0.0078 | 0.5305 | 0 | 619 | 619 |
| `dur_low_pulses_year_mk_pval` | Pulse Metrics | mk_pval | 0.9880 | 0.9954 | 0.0083 | 0.3807 | 57 | 1016 | 959 |
| `dur_low_pulses_year_spearman_pval` | Pulse Metrics | spearman_pval | 0.9881 | 0.9954 | 0.0082 | 0.3647 | 57 | 1016 | 959 |
| `Q99_mk_rho` | Flow Percentiles | mk_rho | 0.9882 | 0.9942 | 0.0037 | 0.3186 | 0 | 570 | 570 |
| `Q99_spearman_rho` | Flow Percentiles | spearman_rho | 0.9882 | 0.9943 | 0.0052 | 0.4645 | 0 | 570 | 570 |
| `n_high_pulses_all_spearman_rho` | Pulse Metrics | spearman_rho | 0.9882 | 0.9938 | 0.0060 | 0.3953 | 1 | 570 | 569 |
| `n_high_pulses_all_senn_slp` | Pulse Metrics | senn_slp | 0.9884 | 0.9889 | 0.0014 | 0.0882 | 0 | 570 | 570 |
| `FDC90th_spearman_rho` | FDC | spearman_rho | 0.9885 | 0.9940 | 0.0069 | 0.7302 | 1 | 598 | 597 |
| `Q1_spearman_pval` | Flow Percentiles | spearman_pval | 0.9885 | 0.9970 | 0.0076 | 0.4284 | 61 | 619 | 558 |
| `Q5_mk_pval` | Flow Percentiles | mk_pval | 0.9886 | 0.9967 | 0.0080 | 0.4986 | 0 | 605 | 605 |
| `Q10_spearman_pval` | Flow Percentiles | spearman_pval | 0.9886 | 0.9965 | 0.0077 | 0.4212 | 32 | 595 | 563 |
| `n_high_pulses_year_linear_slp` | Pulse Metrics | linear_slp | 0.9889 | 0.9959 | 0.0011 | 0.0984 | 0 | 570 | 570 |
| `FDCall_spearman_rho` | FDC | spearman_rho | 0.9890 | 0.9944 | 0.0059 | 0.4483 | 0 | 570 | 570 |
| `D20_day_mk_rho` | Flow Timing | mk_rho | 0.9891 | 0.9947 | 0.0032 | 0.2256 | 0 | 575 | 575 |
| `Dmax_mk_rho` | Flow Timing | mk_rho | 0.9891 | 0.9935 | 0.0031 | 0.1188 | 0 | 575 | 575 |
| `D10_day_mk_rho` | Flow Timing | mk_rho | 0.9891 | 0.9947 | 0.0032 | 0.1809 | 0 | 575 | 575 |
| `D30_day_mk_rho` | Flow Timing | mk_rho | 0.9892 | 0.9946 | 0.0033 | 0.1909 | 0 | 575 | 575 |
| `Dmax_spearman_rho` | Flow Timing | spearman_rho | 0.9892 | 0.9935 | 0.0044 | 0.1795 | 0 | 575 | 575 |
| `D30_day_spearman_rho` | Flow Timing | spearman_rho | 0.9894 | 0.9947 | 0.0047 | 0.2717 | 0 | 575 | 575 |
| `D50_day_mk_rho` | Flow Timing | mk_rho | 0.9894 | 0.9939 | 0.0033 | 0.1438 | 0 | 575 | 575 |
| `D20_day_spearman_rho` | Flow Timing | spearman_rho | 0.9894 | 0.9948 | 0.0045 | 0.3258 | 0 | 575 | 575 |
| `n_low_pulses_all_median` | Pulse Metrics | median | 0.9894 | 0.9961 | 0.0309 | 14.0000 | 0 | 0 | 0 |
| `FDCall_mk_rho` | FDC | mk_rho | 0.9895 | 0.9946 | 0.0041 | 0.3048 | 0 | 570 | 570 |
| `dur_high_pulses_all_linear_slp` | Pulse Metrics | linear_slp | 0.9895 | 0.9900 | 0.0021 | 0.3742 | 0 | 1316 | 1316 |
| `D10_day_spearman_rho` | Flow Timing | spearman_rho | 0.9895 | 0.9947 | 0.0046 | 0.2627 | 0 | 575 | 575 |
| `Q5_spearman_pval` | Flow Percentiles | spearman_pval | 0.9897 | 0.9969 | 0.0074 | 0.4906 | 46 | 605 | 559 |
| `Q90_mk_rho` | Flow Percentiles | mk_rho | 0.9897 | 0.9949 | 0.0039 | 0.3147 | 0 | 570 | 570 |
| `D40_day_mk_rho` | Flow Timing | mk_rho | 0.9897 | 0.9945 | 0.0032 | 0.1522 | 0 | 575 | 575 |
| `D50_day_spearman_rho` | Flow Timing | spearman_rho | 0.9897 | 0.9943 | 0.0047 | 0.2124 | 0 | 575 | 575 |
| `D25_to_D75_spearman_rho` | Flow Timing | spearman_rho | 0.9898 | 0.9942 | 0.0048 | 0.3014 | 0 | 575 | 575 |
| `D60_day_mk_rho` | Flow Timing | mk_rho | 0.9898 | 0.9938 | 0.0033 | 0.1688 | 0 | 575 | 575 |
| `Flow_Reversals_annual_senn_slp` | Pulse Metrics | senn_slp | 0.9899 | 0.9966 | 0.0091 | 1.2922 | 0 | 570 | 570 |
| `Q90_spearman_rho` | Flow Percentiles | spearman_rho | 0.9899 | 0.9950 | 0.0053 | 0.4824 | 0 | 570 | 570 |
| `D60_day_spearman_rho` | Flow Timing | spearman_rho | 0.9899 | 0.9940 | 0.0047 | 0.2390 | 0 | 575 | 575 |

## Substantially Changed Signatures

Columns where mean or median values shifted meaningfully (outside typical cross-language noise).

| Column | Category | R2 | R Mean | Julia Mean | Rel Change | Max Diff |
|--------|----------|----|--------|------------|------------|----------|
| `concavity_mean` | Recession | -1.2615 | 0.4868 | 1.5435 | 217.1% | 4.6119 |
| `concavity_median` | Recession | -0.9618 | 0.4870 | 1.3279 | 172.7% | 4.1101 |
| `b_events_mean` | Recession | -0.2042 | 1.7290 | 2.2270 | 28.8% | 2.7303 |
| `b_events_median` | Recession | -0.1932 | 1.6750 | 2.1266 | 27.0% | 2.8010 |
| `b_pointcloud_median` | Recession | 0.4536 | 1.4256 | 1.3732 | 3.7% | 0.9733 |
| `b_pointcloud_mean` | Recession | 0.4701 | 1.4423 | 1.3963 | 3.2% | 1.2068 |
| `log_a_events_mean` | Recession | 0.7568 | -1.7553 | -1.6517 | 5.9% | 6.8595 |
| `log_a_events_median` | Recession | 0.8235 | -1.8051 | -1.7756 | 1.6% | 4.6172 |
| `log_a_pointcloud_median` | Recession | 0.8429 | -1.8503 | -1.9463 | 5.2% | 1.7515 |
| `log_a_pointcloud_mean` | Recession | 0.8718 | -1.8558 | -1.9393 | 4.5% | 1.5427 |
| `FDCmid_median` | FDC | 0.9451 | -1.4874 | -1.4521 | 2.4% | 10.6113 |
| `dur_low_pulses_all_median` | Pulse Metrics | 0.9729 | 10.7110 | 10.7212 | 0.1% | 95.9167 |
| `elasticity_rolling_median` | Elasticity | 0.9775 | 1.6466 | 1.6495 | 0.2% | 2.7079 |
| `dur_low_pulses_all_mean` | Pulse Metrics | 0.9789 | 14.0086 | 14.0123 | 0.0% | 94.1875 |
| `FDC90th_median` | FDC | 0.9808 | -2.2589 | -2.2447 | 0.6% | 25.8866 |
| `FDC90th_mean` | FDC | 0.9842 | -4.3050 | -4.2654 | 0.9% | 36.7766 |
| `elasticity_rolling_mean` | Elasticity | 0.9843 | 1.6524 | 1.6563 | 0.2% | 2.3568 |
| `n_low_pulses_all_median` | Pulse Metrics | 0.9894 | 2.9253 | 2.9281 | 0.1% | 14.0000 |

## NA Mismatch Analysis

Columns where the number of NAs differs between R and Julia by >50 gages.

| Column | Category | R NAs | Julia NAs | Mismatch | R2 |
|--------|----------|-------|-----------|----------|----|
| `b_pointcloud_senn_slp` | Recession | 4722 | 1095 | 3627 | -44.2829 |
| `b_pointcloud_mk_rho` | Recession | 4722 | 1095 | 3627 | -6.3793 |
| `b_pointcloud_mk_pval` | Recession | 4722 | 1095 | 3627 | -1.5738 |
| `b_pointcloud_median` | Recession | 4722 | 1095 | 3627 | 0.4536 |
| `b_pointcloud_mean` | Recession | 4722 | 1095 | 3627 | 0.4701 |
| `b_pointcloud_linear_slp` | Recession | 4722 | 1095 | 3627 | -37.8971 |
| `b_pointcloud_spearman_pval` | Recession | 4722 | 1095 | 3627 | -1.2825 |
| `b_pointcloud_spearman_rho` | Recession | 4722 | 1095 | 3627 | -4.0713 |
| `log_a_pointcloud_linear_slp` | Recession | 4722 | 1095 | 3627 | -14.6680 |
| `log_a_pointcloud_mean` | Recession | 4722 | 1095 | 3627 | 0.8718 |
| `log_a_pointcloud_median` | Recession | 4722 | 1095 | 3627 | 0.8429 |
| `log_a_pointcloud_mk_pval` | Recession | 4722 | 1095 | 3627 | -1.5428 |
| `log_a_pointcloud_mk_rho` | Recession | 4722 | 1095 | 3627 | -5.7573 |
| `log_a_pointcloud_senn_slp` | Recession | 4722 | 1095 | 3627 | -25.0780 |
| `log_a_pointcloud_spearman_pval` | Recession | 4722 | 1095 | 3627 | -1.2977 |
| `log_a_pointcloud_spearman_rho` | Recession | 4722 | 1095 | 3627 | -3.5642 |
| `dur_low_pulses_all_senn_slp` | Pulse Metrics | 313 | 3782 | 3469 | 0.9921 |
| `dur_low_pulses_all_mk_rho` | Pulse Metrics | 313 | 3782 | 3469 | 0.9921 |
| `dur_low_pulses_all_spearman_rho` | Pulse Metrics | 313 | 3782 | 3469 | 0.9920 |
| `dur_low_pulses_all_spearman_pval` | Pulse Metrics | 313 | 3782 | 3469 | 0.9713 |
| `dur_low_pulses_all_linear_slp` | Pulse Metrics | 313 | 3782 | 3469 | 0.9784 |
| `dur_low_pulses_all_mk_pval` | Pulse Metrics | 313 | 3782 | 3469 | 0.9713 |
| `log_a_seasonality_amplitude_all` | Recession | 4346 | 927 | 3419 | -0.1045 |
| `log_a_seasonality_minimum_first_half` | Recession | 4347 | 928 | 3419 | -0.9329 |
| `log_a_events_mk_pval` | Recession | 4346 | 927 | 3419 | -0.7518 |
| `log_a_events_median` | Recession | 4346 | 927 | 3419 | 0.8235 |
| `log_a_events_senn_slp` | Recession | 4346 | 927 | 3419 | -0.1009 |
| `log_a_events_mk_rho` | Recession | 4346 | 927 | 3419 | -0.3931 |
| `log_a_events_mean` | Recession | 4346 | 927 | 3419 | 0.7568 |
| `log_a_events_linear_slp` | Recession | 4346 | 927 | 3419 | -0.0372 |

*... and 388 more*

## New Julia-Only Columns (not in Golden R)

**43 columns** present in Julia but not in Golden R reference:

### Elasticity (10 columns)
- `elasticity_annual_linear_slp`
- `elasticity_annual_mean`
- `elasticity_annual_median`
- `elasticity_annual_mk_pval`
- `elasticity_annual_mk_rho`
- `elasticity_annual_senn_slp`
- `elasticity_annual_spearman_pval`
- `elasticity_annual_spearman_rho`
- `elasticity_years_low_ppt`
- `elasticity_years_total`

### Flow Timing (16 columns)
- `D1_day_linear_slp`
- `D1_day_mean`
- `D1_day_median`
- `D1_day_mk_pval`
- `D1_day_mk_rho`
- `D1_day_senn_slp`
- `D1_day_spearman_pval`
- `D1_day_spearman_rho`
- `D99_day_linear_slp`
- `D99_day_mean`
- `D99_day_median`
- `D99_day_mk_pval`
- `D99_day_mk_rho`
- `D99_day_senn_slp`
- `D99_day_spearman_pval`
- `D99_day_spearman_rho`

### Negative Flow (8 columns)
- `negative_ann_linear_slp`
- `negative_ann_mean`
- `negative_ann_median`
- `negative_ann_mk_pval`
- `negative_ann_mk_rho`
- `negative_ann_senn_slp`
- `negative_ann_spearman_pval`
- `negative_ann_spearman_rho`

### Recession (8 columns)
- `n_recession_events_linear_slp`
- `n_recession_events_mean`
- `n_recession_events_median`
- `n_recession_events_mk_pval`
- `n_recession_events_mk_rho`
- `n_recession_events_senn_slp`
- `n_recession_events_spearman_pval`
- `n_recession_events_spearman_rho`

### Runoff Ratios (1 columns)
- `runoff_ratio_high_count`


## Deep Dive: Worst 10 Columns

### `b_pointcloud_senn_slp` (Recession)

- **Identity R2**: -44.2829
- **Spearman rho**: 0.1436

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | 0.001943 | 0.000945 |
| Median | 0.000290 | 0.000677 |
| SD | 0.026432 | 0.028493 |
| Min | -0.172141 | -0.615552 |
| Max | 0.266543 | 0.748443 |
| NAs | 4722 | 1095 |

### `b_pointcloud_linear_slp` (Recession)

- **Identity R2**: -37.8971
- **Spearman rho**: 0.1689

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | 0.002438 | 0.000912 |
| Median | 0.000267 | 0.000514 |
| SD | 0.026142 | 0.025425 |
| Min | -0.163360 | -0.615552 |
| Max | 0.279412 | 0.433374 |
| NAs | 4722 | 1095 |

### `log_a_pointcloud_senn_slp` (Recession)

- **Identity R2**: -25.0780
- **Spearman rho**: 0.2186

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | -0.001032 | -0.000993 |
| Median | 0.000177 | -0.000633 |
| SD | 0.033524 | 0.051782 |
| Min | -0.287713 | -0.828539 |
| Max | 0.281845 | 1.948086 |
| NAs | 4722 | 1095 |

### `log_a_pointcloud_linear_slp` (Recession)

- **Identity R2**: -14.6680
- **Spearman rho**: 0.2553

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | -0.000949 | -0.000877 |
| Median | 0.000417 | -0.000505 |
| SD | 0.032212 | 0.037808 |
| Min | -0.287713 | -0.629002 |
| Max | 0.336293 | 0.690223 |
| NAs | 4722 | 1095 |

### `b_pointcloud_mk_rho` (Recession)

- **Identity R2**: -6.3793
- **Spearman rho**: 0.1122

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | 0.018581 | 0.027821 |
| Median | 0.000000 | 0.021739 |
| SD | 0.409393 | 0.218781 |
| Min | -1.000000 | -1.000000 |
| Max | 1.000000 | 1.000000 |
| NAs | 4722 | 1095 |

### `log_a_pointcloud_mk_rho` (Recession)

- **Identity R2**: -5.7573
- **Spearman rho**: 0.2389

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | -0.003835 | -0.024845 |
| Median | 0.000000 | -0.017483 |
| SD | 0.425033 | 0.228329 |
| Min | -1.000000 | -1.000000 |
| Max | 1.000000 | 1.000000 |
| NAs | 4722 | 1095 |

### `b_pointcloud_spearman_rho` (Recession)

- **Identity R2**: -4.0713
- **Spearman rho**: 0.1181

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | 0.022468 | 0.037882 |
| Median | 0.000000 | 0.034274 |
| SD | 0.481912 | 0.289699 |
| Min | -1.000000 | -1.000000 |
| Max | 1.000000 | 1.000000 |
| NAs | 4722 | 1095 |

### `log_a_pointcloud_spearman_rho` (Recession)

- **Identity R2**: -3.5642
- **Spearman rho**: 0.2357

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | -0.000562 | -0.034045 |
| Median | 0.006061 | -0.027009 |
| SD | 0.500537 | 0.299756 |
| Min | -1.000000 | -1.000000 |
| Max | 1.000000 | 1.000000 |
| NAs | 4722 | 1095 |

### `b_pointcloud_mk_pval` (Recession)

- **Identity R2**: -1.5738
- **Spearman rho**: 0.0805

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | 0.628393 | 0.459527 |
| Median | 0.710523 | 0.429795 |
| SD | 0.328492 | 0.323459 |
| Min | 0.000779 | 0.000000 |
| Max | 1.000000 | 1.000000 |
| NAs | 4722 | 1095 |

### `log_a_pointcloud_mk_pval` (Recession)

- **Identity R2**: -1.5428
- **Spearman rho**: 0.0841

| Stat | Golden R | Julia |
|------|----------|-------|
| Count | 975 | 4602 |
| Mean | 0.615529 | 0.449166 |
| Median | 0.707114 | 0.417304 |
| SD | 0.337692 | 0.327243 |
| Min | 0.000699 | 0.000000 |
| Max | 1.000000 | 1.000000 |
| NAs | 4722 | 1095 |

## Summary

| Agreement Tier | Threshold | Count | % |
|----------------|-----------|-------|---|
| Perfect | R2 >= 0.999 | 118 | 21.4% |
| Good | 0.99 <= R2 < 0.999 | 145 | 26.3% |
| Poor | 0.95 <= R2 < 0.99 | 208 | 37.7% |
| Low | 0.90 <= R2 < 0.95 | 15 | 2.7% |
| Very Low | 0.50 <= R2 < 0.90 | 22 | 4.0% |
| Extremely Low | R2 < 0.50 | 43 | 7.8% |
| **Total compared** | | **551** | **100%** |

New Julia-only columns (Section 3 additions): **43**

---
*Generated by `docs/benchmarks/compare_julia_vs_golden_r.py`*

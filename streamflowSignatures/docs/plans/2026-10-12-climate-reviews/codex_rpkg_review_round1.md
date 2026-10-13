**NO-GO.** The port can change annual frost dates, fails the annual bit-identity contract after interpolation, and produces summary statistics outside its documented tolerances. Repository files were not modified; reproductions are under [/tmp/climate_adversarial](/tmp/climate_adversarial).

1. **Severity: BLOCKER — interpolation changes annual results.**  
   **File:line:** [rpkg/R/io.R:456](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/io.R:456), also precipitation interpolation at line 370.

   **Julia/spec:** Julia evaluates `left + (right-left)*offset/span` in that order. Frost requires strictly negative TMIN; annual values must match exactly.

   **Port:** `stats::approx()` uses different floating-point arithmetic.

   **Minimal reproduction:** Complete WY 2001, Q=1, TMAX=15, TMIN=1 except days 99–102: `[-30, NA, NA, 15]`. Both preprocessors accept the year. Julia fills days 100–101 with `[-15, 0]`, giving `frost_last_dowy=100`. R fills them with `[-15, -1.6653345369377348e-15]`, giving **101**.

   The same problem affects precipitation: zeros except days 99–102=`[0, NA, NA, 2.1]` give annual totals **4.2000000000000002** in Julia versus **4.1999999999999993** in R.

   **Fix:** Implement Julia’s interpolation expression and operation order explicitly. Do not repair this by adding a frost tolerance. Add both reproductions as regressions.

2. **Severity: BLOCKER — Spearman p-values use a different statistical method.**  
   **File:line:** [rpkg/R/stats.R:181](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/stats.R:181), reached from `climate.R:232` and `:311`.

   **Julia/spec:** Julia derives the two-sided p-value from the Student-t approximation to Spearman’s rho. The climate summaries must reproduce it within tolerance.

   **Port:** `cor.test(..., method="spearman")` defaults to an exact/AS89 calculation when applicable.

   **Minimal reproduction:** Three complete years with constant daily TMAX values `[1,3,2]`. Annual maxima and rho=0.5 match. Expected `temp_max_annual_spearman_pval` is **0.6666666666666667**; R returns **1**.

   This also affects production-length records: my 30-year fixture produced `prcp_annual_spearman_pval` **0.0190424339** in Julia versus **0.0198555719** in R, despite identical annual values.

   **Fix:** Use `exact=FALSE`, or implement Julia’s calculation explicitly. This is an inherited statistics defect, but it invalidates the climate port.

3. **Severity: BLOCKER — Mann–Kendall loses distinctions and exceeds its relaxed tolerance.**  
   **File:line:** [rpkg/R/stats.R:23](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/stats.R:23).

   **Julia/spec:** Comparisons and tie counts operate on the original doubles. The documented R allowance is `atol=1e-12`, `mk_rtol=1e-6`.

   **Port:** The single-precision `Kendall::MannKendall` backend can collapse distinct doubles into ties and round small p-values to zero. The double-precision constant-series guard does not prevent this.

   **Minimal reproductions:** Twenty complete years with daily TMAX constant within each year:

   - Annual maxima `10 + (1:20)*1e-8`: Julia returns tau=1, p=**8.67976357454836e-10**. R returns tau=1, **p=1**, with `IFAULT=12`.
   - Annual maxima `1:20`: Julia returns the same p-value; R returns **p=0**, exceeding even the relaxed tolerance.

   **Fix:** Calculate S, tie corrections, variance, tau and the normal-tail probability in double precision. Apply the correction to Pettitt segment MK statistics too.

4. **Severity: MAJOR — the complete-grid guard accepts fractional day numbers.**  
   **File:line:** [rpkg/R/climate.R:84](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/climate.R:84).

   **Julia/spec:** Every sorted `dowy[i]` must equal `i` exactly.

   **Port:** `as.integer(dowy)` truncates before comparison.

   **Minimal reproduction:** WY 2001, 365 rows, `dowy=(1:365)+0.5`, PPT=1, TMIN=−1, TMAX=5. Julia rejects the grid and emits missing annual metrics. R accepts it: `prcp_annual=365`, `prcp_50_dowy=183.5`, `frost_first_dowy=1.5`, `frost_last_dowy=272.5`.

   **Fix:** Check finiteness and compare the original numeric values directly with `seq_len(n)`.

5. **Severity: MAJOR — explicit temperature-policy overrides are ignored.**  
   **File:line:** [rpkg/R/io.R:436](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/io.R:436).

   **Julia/spec:** Temperature qualification uses the supplied raw-missing and interpolation-gap limits; Julia exposes both as preprocessing arguments.

   **Port:** The temperature block reads only `pkg_env`, ignoring `config$temperature_na_policy`. The default config construction also omits this section.

   **Minimal reproduction:** Complete WY 2001 with two internal missing TMIN days. Call:

   ```r
   preprocess_daily_data(d, config = list(
     temperature_na_policy = list(
       max_raw_na_per_year_temp = 0,
       max_interpolation_gap_temp = 0
     )
   ))
   ```

   Expected, and Julia’s result with corresponding arguments: no valid temperature years. Actual R: **2001 is accepted and interpolated**.

   **Fix:** Resolve these settings from the supplied config, with package defaults used only when absent.

6. **Severity: MAJOR — config coercion silently changes settings.**  
   **File:line:** [rpkg/R/config.R:214](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/config.R:214), lines 217, 219 and 162.

   **Julia/spec versus port, independently tested using temporary JSON configs:**

   | Input | Julia | R |
   |---|---|---|
   | `multi_day_max_window_days=5.9` | Rejects | Accepts as 5 |
   | `frost_last_end_month=6.9` | Rejects | Accepts as 6 |
   | `max_interpolation_gap_temp=3.5` | Rejects | Accepts |
   | `record_anchored_decade_gate=1` | Enables gate | Silently disables gate |

   The existing range checks run after lossy conversion. Neither implementation rejects a negative temperature gap setting.

   **Fix:** Validate scalar types, finiteness, integrality and ranges before conversion. Define boolean handling consistently; reject invalid temperature-policy values in both implementations.

7. **Severity: MAJOR — threshold names diverge and can collide, overwriting metrics.**  
   **File:line:** [rpkg/R/climate.R:47](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/climate.R:47), assignment at line 204.

   **Julia/spec:** Names identify distinct configured thresholds.

   **Port:** R’s default numeric formatting loses significant digits; integer conversion is limited to 32 bits.

   **Minimal reproductions:**

   | Threshold | Julia label | R label |
   |---|---|---|
   | `1e-7` | `1p0e-7` | `1e-07` |
   | `2147483648` | `2147483648` | `NA` |

   More seriously, valid ascending thresholds `[1.234567890123451, 1.234567890123452]` both become `1p23456789012345` in R. With one precipitation day equal to the larger threshold, expected counts are **1 and 0**. R overwrites the shared column and emits **0 twice**.

   **Fix:** Match Julia’s numeric-label serialization and validate uniqueness of generated names before calculating metrics.

8. **Severity: MAJOR — configured season definitions are ignored.**  
   **File:line:** [rpkg/R/climate.R:25](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/climate.R:25), used at lines 175 and 287.

   **Julia/spec:** Season membership comes from `na_handling.seasonal_completeness.season_definitions`.

   **Port:** Climate calculations always use the hard-coded month lists.

   **Minimal reproduction:** Configure winter as `[1]`, then supply WY 2001 with PPT=1 daily. Julia uses the configured mapping and returns `prcp_win=31`. R accepts the config but returns **90**. Temperature seasonal extrema and GDD have the same wiring problem.

   **Fix:** Load and use the configured mappings. Preserve DJF/MAM/JJA/SON as defaults.

9. **Severity: MAJOR — the runner silently drops temperature-only input.**  
   **File:line:** [docs/benchmarks/run_rpkg_benchmark.R:272](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_rpkg_benchmark.R:272).

   **Julia/spec:** Temperature qualification and the explicit `temp_data` route do not require precipitation.

   **Port:** The runner merges all climate columns only when PPT exists. TMIN/TMAX can be successfully read and normalized, then discarded.

   **Minimal reproduction:** I ran the actual benchmark with one qualifying gage, complete streamflow, and a climate parquet containing `site_id, Date, tmin, tmax`, without precipitation. Expected: temperature output. Actual: **zero temperature columns**, with `processed=1`, `errors=0`. Adding `prcp` to the same file produced **304 temperature summary columns**.

   **Fix:** Merge available climate variables independently; derive precipitation availability separately from temperature availability. The library’s explicit temperature-only call already works.

10. **Severity: NOTE — both implementations admit infinite precipitation into qualified years.**  
    **File:line:** [rpkg/R/io.R:343](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/io.R:343), [rpkg/R/climate.R:169](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/rpkg/R/climate.R:169); corresponding Julia guards have the same weakness.

    **Julia/spec:** This is a shared numerical robustness problem, not a port discrepancy. A fraction-of-total timing is not meaningful when the total is infinite.

    **Minimal reproduction:** Complete WY 2001, Q=1, PPT=1 except day 100=`Inf`. Both qualify the precipitation year and return `prcp_annual=Inf`, **`prcp_50_dowy=100`**, and undefined centre of mass/Gini.

    **Fix:** Treat non-finite precipitation as missing during preprocessing and require finite inputs in the annual kernels. Coordinate this change with Julia.

**What I verified and found correct**

- **Complete finite annual calculations:** On an independently constructed, varied 30-year fixture, all **1,158 finite annual values were byte-identical**, with 12 matching missing values. Nevertheless, **22 summary fields failed** the documented tolerances.
- **Single storm:** WY 2001, 100 mm on day 100 gives maximum, median-volume timing and centre of mass at day 100; longest dry spell=265 days, days 101–365; sample SD=`100/√365`≈5.234239225902137; CV=`√365`; Gini=`364/365`; SI=`11/6`.
- **Zero precipitation:** WY 2000 gives longest dry spell=366, start=1, end=366, maximum timing=1; fractional timing, centre of mass, CV, Gini and SI are missing.
- **Calendar boundaries:** Truncated leap years and duplicated interior days are rejected. All-frost WY 2000 ends its permitted frost search at day **274**; TMIN=−1/TMAX=5 gives annual GDD=732 and winter GDD=182.
- **Daymet omission:** Removing December 31, 2000, from a linear daily series reconstructs day 92 correctly in both preprocessors; all 39 annual outputs match.
- **Seasons and ties under shipped config:** Fall includes October+November+September; bogus `month` values and shuffled rows do not change results. First maxima and first longest dry spells win ties. Wet-day thresholds are strict.
- **Sequential accumulation:** The explicit summation kernels preserve left-to-right double accumulation. The two separated 50-day blocks of 0.1 mm reach half the annual total on day **50**.
- **Temperature qualification:** Three-day internal gaps and 30 isolated missing days pass; four-day gaps, 31 missing days and boundary gaps fail. Internal infinities are treated as missing. Both series must pass, write-back is atomic, and streamflow-rejected years cannot qualify.
- **Record gate:** Seven qualifying frost observations out of ten suppress trends; eight pass. Frost-free record endings suppress frost trends while preserving the descriptive statistics.
- **Orchestration and additivity:** Explicit-frame requirements, empty temperature frames and independence from `area_normalized` work. Adding, removing or invalidating temperature columns left compared existing-family outputs unchanged.
- **Runner window wiring with PPT present:** The actual runner retained temperature years 1995–2015, excluding deliberately invalid WY 2000 and every out-of-window year.

The supplied test file passed **81 assertions**; those assertions do not establish the required mirror contract.

**Verdict: NO-GO for use as a validated mirror of Julia.**
# Summary Documentation for Streamflow Signatures (in helperFunctions.R)

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/d/e/2PACX-1vQnt7OCPm19vnWF4yynXL9JTzTvq9CrGoEaDv7yFSngLoFsypiWsx6fZLKWwaO5YQ/pub)
> **Last synced**: 2026-08-28 (content unchanged since 2026-07-21)

github repo [here](https://github.com) or [here](https://github.com)

Summary data audit sheet [here](https://docs.google.com)

## NA Handling (for the whole analysis):

- Do not replace NAs with zeros; streamflow that is zero is okay
- NAs can be for various reasons, we need to use the USGS codes that describe what creates those NAs (e.g. Ice affected, etc.) for each site that has ice effected days we need a total count of the number of days that are ice affected
- If there is up to 3 days of continuous missing streamflow data, interpolate between the points; tell us if there is a data driven answer to the number of days that would not create problems for trend analysis
- From USGS data release: remove years with (i) more than 3 consecutive days of NAs, (ii) gaps of more than 30 days in a year, (iii) negative values, (iv) constant monthly standard deviations during periods of non-zero streamflow indicative of data errors. Items i, iii, and iv are set in the config to flag (not remove) by default.
- Flag years with constant monthly standard deviations during periods of non-zero streamflow indicative highly controlled discharges.
- Negative_ann: Calculate the number of days with negative values per gage per year
- Flow_reversals_ann: Count number of flow reversals annually: changes in flow direction exceeding 2% of current flow.
- For a trend to be calculated, at least 60% of the entire annual streamflow metric time series must be complete, at least 80% of the first decade must be complete, and at least 80% of the last decade of all trend periods
- For any seasonal or annual metrics the period must have at least 80% of the data.
- Seasonal periods are defined as:
  - Winter: December-February
  - Spring: March-May
  - Summer: June-August
  - Fall: September-November
- All annual metrics should be calculated on the water year October 1 - September 30

## 1. Glossary of Streamflow Signatures

### Flow Volume Metrics

- Qann: Annual total streamflow.
- Qwin: Winter total streamflow.
- Qspr: Spring total streamflow.
- Qsum: Summer total streamflow.
- Qfal: Fall total streamflow.
- Qxx: Flow percentiles (e.g., Q10, Q50, Q90).

### Baseflow Metrics

- BFI_Eckhardt: Baseflow index calculated using the Eckhardt filter.
- BFI_LyneHollick: Baseflow index calculated using the Lyne-Hollick filter.

### Flow Duration Curve (FDC) Metrics

- FDCall: Slope of the annual flow duration curve across the full range of exceedance probabilities.
- FDC90th: Slope of the annual flow duration curve over the low-flow tail (exceedance probability >= 0.90).
- FDCmid: Slope of the annual flow duration curve over the midsegment (exceedance probability 0.20-0.80).

### Recession Metrics

- log_a_pointcloud: Median log(a) from point cloud analysis of all recession events.
- log_a_events: Median log(a) from individual recession events.
- b_pointcloud: Median recession exponent b from point cloud analysis.
- b_events: Median b from individual events.
- concavity: Difference in b between the first and second halves of recession events.
- log_a_seasonality_amplitude_all: Seasonal amplitude of log(a) for all years.
- log_a_seasonality_minimum_all: water year date of minimum log(a) for all years.

### Pulse Metrics

- n_high_pulses_year: Number of high-flow pulses per year.
- n_low_pulses_year: Number of low-flow pulses per year.
- dur_high_pulses_year: Mean duration of high-flow pulses per year.
- dur_low_pulses_year: Mean duration of low-flow pulses per year.
- TQmean: Percentage of days with flow above the annual mean.
- Flow_Reversals: Number of flow reversals (annual or seasonal).

### Flashiness Metrics

- R-B Index: Richards-Baker flashiness index.

### Flow Timing Metrics

- Dxx_day: water year day when cumulative flow reaches a given percentile (e.g., D10, D50, D90).
- D25_to_D75: Days between 25% and 75% cumulative flow.
- Dmax: water year day of maximum flow.

### Storage Metric

- avg_storage: Mean annual catchment storage (mm) derived from water balance.

### Elasticity Metrics:

- elasticity_static: Overall catchment elasticity - median of annual elasticity values. Measures how sensitively streamflow responds to precipitation changes. Values ~1.0 indicate proportional response; >1 indicates amplified response.
- elasticity_rolling: Rolling window (11-year) elasticity values, used to detect trends in catchment sensitivity over time.

### Runoff Ratio Metrics:

- annual_runoff_ratio: Annual water year streamflow divided by annual precipitation (Q/P).
- winter_runoff_ratio: Winter (Dec-Feb) streamflow divided by winter precipitation.
- spring_runoff_ratio: Spring (Mar-May) streamflow divided by spring precipitation.
- summer_runoff_ratio: Summer (Jun-Aug) streamflow divided by summer precipitation.
- fall_runoff_ratio: Fall (Sep-Nov) streamflow divided by fall precipitation.

## 2. Glossary of Statistical Metrics

- slp: Theil-Sen slope, a robust linear trend estimator.
- rho: Spearman's rank correlation coefficient, a non-parametric measure of monotonic trends.
- pval: P-value associated with rho.
- mean: Arithmetic mean of the signature across all years.
- median: Median value of the signature across all years.

## 3. Summary of Functions

### process_gages_rawData()

**Purpose:** Processes raw gage data (USGS/Canadian) to calculate streamflow metrics for individual gages and save results to a file. Data are available up to the present day, but no climate data are directly available.

**Requirements and Decisions:**
- A folder of USGS and HYDAT (Canadian) metadata must be stored locally in order to process watersheds.
- A minimum number of valid years (min_num_years) is required for a gage to be included.
- Data is filtered by a minimum flow threshold (min_Q_value_and_days). A watershed with fewer than the minimum days of > min_Q_value is not processed, but if the threshold is surpassed then non-zero flow is still passed to subsequent functions for processing.
- Missing or invalid years are excluded from analysis.

### process_caravan_gages()

**Purpose:** Processes Caravan NetCDF timeseries data for watersheds, calculates streamflow metrics, and saves results. Coincident climate data are available, but the record of streamflow is truncated to as early as the early 2000s (as late as 2018 for HYSETS).

**Requirements and Decisions:**
- Minimum years of valid data (min_num_years) are required for processing.
- Filters years based on minimum flow thresholds and valid days.
- Handles redundancy between datasets (e.g., CAMELS and HYSETS).

### generate_streamflow_dt() and generate_streamflow_dt_caravan()

**Purpose:** Converts raw gage data (USGS/Canadian) or Caravan NetCDF data into a standardized streamflow data.table.

**Requirements and Decisions:**
- Streamflow is provided in appropriate units (or converted, e.g., m³/s to mm/day).
- Filters data by specified date ranges and required years of data.
- Handles missing or invalid values (e.g., flagged data).
- Adds derived columns (e.g., year, month, day of year).

### calculate_flow_vols_by_year()

**Purpose:** Calculates annual and seasonal flow volumes, as well as flow percentiles (e.g., Q1, Q10, Q50, Q90).

**Requirements and Decisions:**
- Requires year, Q, month, and doy columns in input data.
- frozen days are okay (USGS data flag and Hydat data flag requirement); 0 flow days are okay
- Any season with less than 80% of the data should not be included, count the number of years where the season does not meet that threshold

### analyze_baseflow_indices()

**Purpose:** Calculates baseflow indices (BFI) using Eckhardt and Lyne-Hollick digital filters.

**Requirements and Decisions:**
- T-1 is part of the moving window that the Lyne-Hollick uses for baseflow separation. Data continuity could be an issue with the way it is calculated
- linear interpolation of streamflow for data gaps less than 3days;
- Years with data gaps >3 days are rejected by the centralized preprocessor, so baseflow filters never encounter gaps in valid years. No filter restart or gap counting is needed.
- Default parameters for Eckhardt filter: BFImax = 0.8, a = 0.98.
- Once this has been run, make sure each day of the baseflow time series is equal to or less than the daily streamflow
- BFI must be between 0 and 1

### analyze_baseflow_incidies_with_parameters()

**Purpose:** Calculate baseflow indices using the recession parameters below as inputs to (BFI) using Eckhardt and Lyne-Hollick digital filters.

**Requirements and Decisions:** same as the baseflow index decisions above

### analyze_recession_parameters()

**Purpose:** Analyzes recession events to calculate parameters like log(a), b, and concavity.

**Requirements and Decisions:**
- Recession events are defined as periods of at least 5 consecutive days with monotonic decreases in flow (Q) and its derivative (|dQ/dt|).
- Recession events are identified as the longest contiguous windows where both Q and |dQ/dt| are monotonically decreasing, with a minimum length of 5 days. The first day of each event is removed during power law fitting (standard practice — storm peak influence).
- Split recession events temporally (first vs. second half) for concavity analysis.
- Fits sinusoidal models to analyze seasonal variation in log(a).
- Count the number of recession events in each year
- b is determined using the line of best fit to the log-log plot of −dq̂/dt versus q̂.
- a is determined by setting b = 1 (linear aquifer) and using the line of best fit to the log-log plot of −dq̂/dt versus q̂.
- To control the quality of fitted parameters, calculate recession fits and create a flag for any R2 < 0.8

### analyze_flashiness()

**Purpose:** Calculates the Richards-Baker (RB) flashiness index by year.

- Uses day to day changes

### calculate_pulse_metrics()

**Purpose:** Calculates metrics for high/low flow pulses (e.g., frequency, duration) and flow reversals.

**Requirements and Decisions:**
- This metric requires data from across the period of record, confirm that the code can actually calculate the 90th/10th percentiles from the full period of record before the annual metrics.
- High/low flow pulses are defined using the 90th and 10th percentiles of flow (Q).
- Calculate annual metrics using period-of-record percentiles

### analyze_flow_timing()

**Purpose:** Analyzes the timing of flow events (e.g., day of maximum flow, cumulative flow percentiles).

**Requirements and Decisions:**
- Calculate day of water year for cumulative flow percentiles, and the day with maximum flow
- Dxx_day: day of water year when cumulative flow reaches a given percentile (e.g., D10, D50, D90). D1,5,10,20,30,40,50,60,70,80,90,95,99
- D25_to_D75: Days between 25% and 75% cumulative flow.
- Dmax: water year day of maximum flow.

### analyze_runoff_ratios()

**Purpose:** Analyzes runoff ratios (streamflow divided by precipitation) and their trends.

**Runoff Ratio Definition:**
- Runoff_ratio: calculated as annual or seasonal periods with total streamflow divided by total precipitation (Q/P).

**Requirements and Decisions:**
- Precipitation (PPT) data must be available alongside streamflow (Q).
- Same QAQC procedures for Q should apply to NAs (missing data) in precipitation; zero is okay in precipitation data
- Calculates separate runoff ratios for annual and seasonal periods (winter, spring, summer, fall).
- Cases where annual PPT < 10mm or seasonal PPT < 1mm the ratio is set to NA
- Runoff_ratio_high: flag runoff ratios that are greater than two

### calculate_streamflow_elasticity() (Kendra and Deni)

**Purpose:** Calculates streamflow elasticity, which measures how sensitive streamflow is to changes in precipitation. Based on Sawicz et al. (2011).

**Definitions:**
- elasticity_static: Overall catchment elasticity - median of annual elasticity values. Measures how sensitively streamflow responds to precipitation changes. Values ~1.0 indicate proportional response; >1 indicates amplified response.
- elasticity: Rolling window (11-year) elasticity values, used to detect trends in catchment sensitivity over time.
- elastacity_annual: measures year-to-year sensitivity rather than deviation-from-mean sensitivity

**Requirements and Decisions:**
- Requires daily streamflow (Q) and precipitation (PPT) data. Calculating differences between years, so need minimum of two years of "good data"
- Minimum of 15 years of valid data required for elasticity calculation.
- Flag years with annual precipitation <10mm and exclude from analysis.
- Calculates both static elasticity (median of all annual values) and rolling window elasticity (11-year window by default).
- Annual elasticity formula: E = (dQ/dP) / (Q_mean/P_mean), where dQ = Q_year - Q_mean and dP = P_year - P_mean. This is like an anomaly direction, while the approach from the paper is accounting for anomaly (lets create a second calculation to do that) — achieved as of 4/16
- Distance from average v.s. Accounting for antecedent conditions
- Values ~1.0 indicate proportional response; values >1 indicate amplified streamflow response to precipitation changes.
- From Sawicz 2011 slightly different interpretation from above based on dq= t1-t0 :  A value of 1 indicates that a 1 % precipitation change leads to a 1 % change in streamflow. A value greater or less than 1 would, respectively, define the catchment as being elastic, i.e., sensitive to change of precipitation, or inelastic, i.e., insensitive to a change of precipitation.
- to go forward with the additional calculation of using the previous years precip/Q then you would need to pass the minimum data requirements for both years (consecutive years)
- Year qualification handled by centralized preprocessor (config).
- Count the number of years each site doesn't have sufficient data within a year, and what the number of excluded years would be if the value was <30% data missing — added as documentation not filter
- Output: Returns 'elasticity_static' (single value) plus the 8 standard metric statistics (mean, median, linear trend, theil sen, p value); 'elasticity_rolling' (one value per 11 year windows); 'elasticity_annual' plus the 8 standard stats

### calculate_qp_seasonality()

**Purpose:** Quantifies the seasonality in the relationship between cumulative streamflow (Q) and cumulative precipitation (P). Based on Wrede et al. (2015), however rolling window approach described below may differ from original publication.

**Q-P Seasonality Metric Definitions**
- qp_slope_sd: Standard deviation of monthly cumulative Q-P slopes. Higher values indicate stronger seasonal variation in the streamflow-precipitation relationship.
- qp_bimodality: Bimodality coefficient of the Q-P slope distribution. Values >0.555 suggest bimodal or strongly seasonal patterns.

**Requirements and Decisions:**
- Requires daily streamflow (Q), precipitation (PPT), month, and day-of-water-year (dowy) columns.
- Year qualification handled by centralized preprocessor (config).
- Calculates 30-day rolling slope (backward-looking) of cumulative Q vs cumulative P for each water year.
- Aggregates rolling slopes to 12 monthly mean values per year.
- Two metrics are calculated per year:
  - qp_slope_sd: Standard deviation of the 12 monthly slopes (higher values = more seasonal variation).
  - qp_bimodality: Bimodality coefficient = (skewness² + 1) / kurtosis. Values >0.555 suggest bimodal/seasonal patterns across all water years.

### calculate_average_storage() — OMITTING VARIABLE FROM MAJOR ANALYSES 4/23/26

**NOTES From 4/16:** 1) we see 3 options here, go with metric as is (low interpretability/maybe not that useful), 2) create water balance w/ ET data from independent source (GLEAM; 0.1 degree res.), 3) edit metric below but need to define the nitty gritty of dormant season, considering snow.
- "Annual Water Balance" approach would use actual ET, daily precip (liquid + solid), and runoff
- Caveat: Negative values would have to be considered and dealt with

**Purpose:** Estimates mean annual catchment storage using a simplified water balance approach. Based on Peters & Aulenbach (2011). Note: current implementation does not account for AET

**Storage Metric Definition:** avg_storage: Mean annual catchment storage (mm) derived from water balance.

**Requirements and Decisions:**
- Requires daily streamflow (Q), precipitation (PPT), and evapotranspiration (PET) and day-of-water-year (dowy) columns.
- Year qualification handled by centralized preprocessor (config).
- Need to define dormant season
- Need to have parameter for initial storage?
- Annual storage can't be greater than annual precipitation
- Daily water balance: dS = P - Q (evapotranspiration is not estimated).
- Cumulative storage: S = cumsum(dS) for each water year.
- [account for aet losses??!!??]
- Units: millimeters (mm).

**Notes from Erin:** storage defined by the range of change in cumulative catchment storage that begins during dormant season streamflow recession
- Change in storage is simple WB: P-R-E
- E is calculated as PET
- Calculated as daily time step
- Storage for storage-discharge relation derived using WB equation added to an initial estimate of storage
- Storage was only estimated for baseflow/recession days (defined as days having no precip on the day/for X days)
- Although storage-discharge relation was developed using only recession days, Annual storage from WB was tracked continuously during dormant season, when ET was at minimum
- Had to make assumption about initial storage for each dormant season (assumption was that baseflow is increased single-valued function of storage described by exponential function between storage and Q)
- If we want storage (not change in storage) need to add P to an intiial estimate of storage and substract PET and runoff → what is this initial value of storage for all catchments?
- Questions: 1) how would snow be factored into this?
- 2) is the cumulative storage estimated only during dormant season/baseflow?
- 3) how do we estimate the initial storage value to get to a storage estimate (rather than just a delta storage estimate)

## Part 2: utility functions

### integrate_daymet_with_streamflow()

**Purpose:** Joins Daymet climate data (precipitation, temperature) to streamflow data based on date matching.

**Requirements and Decisions:**
- Requires a pre-processed Daymet parquet file (`data_out/daymet_1980_2023.parquet`).
- Matches streamflow gage IDs to Daymet site IDs.
- Joins climate variables (PPT, tmin, tmax, swe, vp, srad) to streamflow data.table on Date.
- Reports coverage percentage; warns if <95% of dates have matching climate data.
- Renames `prcp` column to `PPT` for compatibility with existing Q-PPT analysis functions.

### convert_daymet_zip_to_parquet()

**Purpose:** One-time conversion of Daymet ZIP archive (containing 44 annual CSV files) to a single optimized parquet file.

**Requirements and Decisions:**
- Input: ZIP file containing daymet_1980.csv through daymet_2023.csv.
- Reconstructs the Date column from year, month, and row order (original CSVs lack a day column).
- Validates that each site has 365 or 366 days per year.
- Uses snappy compression for efficient storage (~70% reduction vs uncompressed).
- Output columns: site_id, Date, prcp, tmin, tmax, swe, vp, srad.
- Output:
- Single parquet file (~3.9 GB) containing ~98 million rows for 6,087 sites.

## Part 3: Definitions

### Statistical Metrics

| Suffix | Metric | Method | Purpose |
|--------|--------|--------|---------|
| _senn_slp | Theil-Sen slope | zyp::zyp.sen | Robust non-parametric trend detection |
| _linear_slp | Linear regression slope | lm() | Parametric trend for comparison |
| _spearman_rho | Spearman's rho | cor.test | Rank correlation with time |
| _spearman_pval | Spearman p-value | cor.test | Significance of Spearman correlation |
| _mk_rho | Mann-Kendall tau | Kendall::MannKendall | Non-parametric trend strength |
| _mk_pval | Mann-Kendall p-value | Kendall::MannKendall | Significance of Mann-Kendall test |
| _mean | Arithmetic mean | mean() | Central tendency across all years |
| _median | Median | median() | Robust central tendency |

## Part 3: References

- Peters, N. E., & Aulenbach, B. T. (2011). Water storage at the Panola Mountain Research Watershed, Georgia, USA. Hydrological Processes, 25(25), 3878-3889.
- Sawicz, K., Wagener, T., Sivapalan, M., Troch, P. A., & Carrillo, G. (2011). Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA. Hydrology and Earth System Sciences, 15(9), 2895-2911.
- Wrede, S., Fenicia, F., Martinez-Carreras, N., Juilleret, J., Hissler, C., Krein, A., ... & Pfister, L. (2015). Towards more systematic perceptual model development: a case study using 3 Luxembourgish catchments. Hydrological Processes, 29(12), 2731-2750.

## Part 4: Automated Data Quality Flags

| Flag Column | Condition |
|---|---|
| flagged_for_qann_range | Qann_mean outside [0, 2000] |
| flagged_for_bfi_eckhardt_range | BFI outside [0, 1] |
| flagged_for_bfi_lynehollick_range | BFI outside [0, 1] |
| flagged_for_flashiness_range | Flashiness outside [0, 2] |
| flagged_for_tqmean_range | TQmean outside [0, 100] |
| flagged_for_d50_range | D50 outside [1, 366] |
| flagged_for_elasticity_range | Elasticity outside [0.1, 5] |
| flagged_for_runoff_ratio_range | Runoff ratio outside [0.01, 1.5] |
| flagged_for_seasonal_sum | Sum not within 20% of Qann |
| flagged_for_percentile_order | Q5<Q25<Q50<Q75<Q95 violated |
| flagged_for_timing_order | D5<D50<D95 violated |
| flagged_for_high_na | >30% signatures are NA |

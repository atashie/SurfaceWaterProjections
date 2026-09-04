# Summary Documentation for Streamflow Signatures

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub)
> **Last synced**: 2026-09-04 (second sync of the day, after the Parts 4–5 accuracy review;
> three text edits since the morning sync are mirrored: the
> `convert_daymet_zip_to_parquet()` block was removed from Part 4, the Daymet join
> line now lists `prcp` rather than `PPT`, and the Part 5 preamble now reads "flags never
> remove a gage from the downstream analysis (data are flagged only)". The legacy text
> is now preceded by a "START: ORIGINAL DOCUMENT" marker and a banner of "!" lines,
> not reproduced below.) The Google Doc now carries the RESTRUCTURED document
> (pasted from `docs/plans/2026-08-31-guidelines-restructure-draft.md` with edits)
> between "START: NEW DOCUMENT" / "END: NEW DOCUMENT" markers, followed by the
> previous (legacy) document text, which is unchanged since the 2026-08-31 sync.
> Both parts are mirrored below in that order. Differences between the pasted
> new document and the restructure draft (48 regions) are summarized in
> CHANGELOG → Guidelines Document TODOs (2026-09-04); the notable ones: the
> Storage module (draft §3.11) was NOT carried over (modules renumbered
> 3.11 Negative Flow Days … 3.13 Streamflow Drought — 13 modules, not 14), the
> Pettitt evaluation-window paragraph (WY 1980–2024, ≥20 obs / ≥10 per segment)
> was dropped, the header repo/audit-sheet links are still generic "here" links,
> and several snow definitions were reworded.

---

## NEW DOCUMENT (between the doc's START/END markers)

# Summary Documentation for Streamflow Signatures

GitHub repo here· Summary data audit sheet here

How this document is organized. Part 1 holds the global rules that apply to every signature (time conventions, NA handling, completeness requirements). Part 2 defines the 16 statistics computed for every signature. Part 3 documents each signature family as a self-contained module with a fixed structure — Metrics (name: definition, units), Function, Method and parameters, and Requirements and decisions — so a reader can take in one family without cross-referencing. Part 4 covers data ingestion and utility functions, Part 5 the automated quality flags, and Part 6 the references.

## Part 1: Global Rules (apply to every signature)

### 1.1 Time conventions

All annual metrics are calculated on the water year, October 1 – September 30.

Seasonal periods are defined as: Winter: December–February; Spring: March–May; Summer: June–August; Fall: September–November.

### 1.2 NA handling and per-year screening

Do not replace NAs with zeros; streamflow that is zero is okay.

NAs can be for various reasons, we need to use the USGS codes that describe what creates those NAs (e.g. Ice affected, etc.); for each site that has ice affected days we need a total count of the number of days that are ice affected.

If there is up to 3 days of continuous missing streamflow data, interpolate between the points; tell us if there is a data driven answer to the number of days that would not create problems for trend analysis.

Following the USGS data release QA criteria, four conditions are screened per water year: (i) more than 3 consecutive days of missing data, (ii) more than 30 total missing days in the year, (iii) negative flow values, and (iv) constant monthly standard deviations during periods of non-zero streamflow (indicative of data errors or highly controlled discharge). Items i and ii REMOVE the year: any gap longer than 3 consecutive days, or more than 30 missing days in total, excludes that water year from all calculations. Internal gaps of up to 3 days are filled by linear interpolation for years that pass, but interpolated days still count toward the 30-missing-day limit (the count is taken on the raw data, before interpolation). Items iii and iv FLAG but do not remove: negative values are retained by default (config option reject_negative_flow, default false) and are counted in Negative_ann (see 3.12); a constant monthly standard deviation sets a QA flag only and never rejects a year.

Flag years with constant monthly standard deviations during periods of non-zero streamflow, indicative of highly controlled discharges.

### 1.3 Completeness requirements

For any seasonal or annual metrics the period must have at least 80% of the data.

For a trend to be calculated, at least 80% of the first decade must be complete, AND at least 80% of the last decade of all trend periods, AND at least 60% of the entire annual streamflow metric time series must be complete.

Statistics floor: a signature must have at least 20 non-NA annual values before ANY of its 16 statistics are computed (below 20, all statistics are NA). This harmonizes with the Pettitt test's own 20-observation minimum.

Exemptions: the recession and elasticity families are exempt from the trend-completeness requirements and the 20-value floor — both are inherently sparse (recession is event-based; the rolling elasticity window consumes 10 years of record by construction).

Gage inclusion (standard products): a gage is included only if it has at least 20 qualifying water years AND at least 60% of the possible water years in the analysis window.

## Part 2: Statistical Metrics (16 per signature)

The following 16 statistics are calculated for every streamflow signature metric (except where otherwise noted):

senn_slp: Theil-Sen slope — robust non-parametric trend estimate (median of pairwise slopes), in signature units per year.

linear_slp: Ordinary least-squares linear regression slope, in signature units per year; parametric companion to the Theil-Sen slope.

spearman_rho: Spearman's rank correlation between the annual values and water year — non-parametric measure of monotonic trend strength (−1 to +1).

spearman_pval: P-value of the Spearman correlation.

mk_rho: Mann-Kendall tau (tie-corrected tau-b) — non-parametric measure of monotonic trend strength (−1 to +1).

mk_pval: P-value of the Mann-Kendall test (continuity-corrected normal approximation with tie-corrected variance).

mean: Arithmetic mean of the signature's annual values across all qualifying water years.

median: Median of the signature's annual values across all qualifying water years.

pettitt_cp_year: Most likely changepoint year. Always reported when the test runs; use pettitt_pval to judge significance.

pettitt_pval: Asymptotic p-value of the Pettitt test (< 0.05 indicates a significant changepoint).

pettitt_pre_mean: Mean of the annual values before the changepoint year.

pettitt_post_mean: Mean of the annual values after the changepoint year.

pettitt_delta_mean: pettitt_post_mean − pettitt_pre_mean (signature units).

pettitt_pct_change: pettitt_delta_mean as a percentage of |pettitt_pre_mean|.

pettitt_pre_mk_pval: Mann-Kendall p-value for the pre-changepoint segment (tests for a residual trend within the segment).

pettitt_post_mk_pval: Mann-Kendall p-value for the post-changepoint segment.

Exceptions to the 16-statistic rule. A small number of outputs are single-valued per gage and carry no statistics: elasticity_static, the recession seasonality values, the recession-derived filter constant, the drought thresholds, and the diagnostic counts (each is noted in its module below).

## Part 3: Signature Modules

Each module follows the same structure: Metrics (name: definition, units) → Function → Method and parameters→ Requirements and decisions.

### 3.1 Flow Volumes and Percentiles

Function: calculate_flow_vols_by_year()

Purpose: Calculates annual and seasonal flow volumes, as well as flow percentiles.

Metrics:

Qann: Annual total streamflow, in mm (sum of daily mm/day values).

Qwin, Qspr, Qsum, Qfal: Winter, spring, summer, and fall total streamflow, in mm.

Qxx: Flow percentiles (Q1, Q5, Q10, Q20, Q25, Q30, Q40, Q50, Q60, Q70, Q75, Q80, Q90, Q95, Q99), in mm/day.

Q95_Q10: Difference between the high- and low-flow percentiles (Q95 − Q10), in mm/day.

season_excluded_years_winter / _spring / _summer / _fall: Count of years where that season failed the 80% completeness threshold (per-gage diagnostic scalars).

Requirements and decisions:

Requires year, Q, month, and day of year columns in input data.

Frozen days are okay (USGS data flag and HYDAT data flag requirement); 0 flow days are okay.

Any season with less than 80% of the data should not be included; count the number of years where the season does not meet that threshold.

### 3.2 Flow Duration Curve (FDC)

Function: analyze_fdc_trends()

Purpose: Fits the slope of each water year's empirical flow duration curve, in log space, over three exceedance-probability ranges. Steeper (more negative) slopes indicate flashier, more variable regimes.

Metrics:

FDCall: Slope of the annual flow duration curve across the full range of exceedance probabilities.

FDC90th: Slope of the annual flow duration curve over the low-flow tail (exceedance probability ≥ 0.90).

FDCmid: Slope of the annual flow duration curve over the midsegment (exceedance probability 0.20–0.80).

Method and parameters:

Within each water year, NaN and negative daily flows are dropped (zeros retained); at least 10 valid daily values are required, otherwise all three slopes are NA for that year.

Valid flows are sorted descending and assigned exceedance probabilities via the Weibull plotting position p = i / (n + 1).

Flows are log-transformed as log10(Q + 1e-10) so zero-flow days remain finite.

Each segment's slope is the coefficient of an ordinary least-squares regression of log10(Q) on exceedance probability; a segment with fewer than 3 points yields NA. Regressing over a full segment (rather than a two-point percentile difference) is less sensitive to endpoint choice.

Units: log10(mm/day) per unit exceedance probability; slopes are typically negative.

### 3.3 Baseflow

Functions: analyze_baseflow_indices(), analyze_baseflow_indices_with_parameters()

Purpose: Calculates baseflow indices (BFI) using the Eckhardt and Lyne-Hollick digital filters, each in a default and a recession-parameterized variant.

Metrics:

BFI_Eckhardt: Baseflow index from the Eckhardt filter with default parameters. BFImax = 0.8, a = 0.98.

BFI_LyneHollick: Baseflow index from the Lyne-Hollick filter with default parameter a = 0.925.

BFI_Eckhardt_param: Eckhardt filter with the site-specific recession-derived constant in place of a = 0.98.

BFI_LyneHollick_param: Lyne-Hollick filter with the recession-derived alpha in place of a = 0.925. This parameterization is heuristic — the L-H parameter has no formal derivation from recession analysis.

recession_alpha_point_cloud: Median of day-to-day recession ratios Q(i+1)/Q(i) across all recession events, under the linear-reservoir assumption (b = 1) for each gage. This constant parameterizes the two _param filters.

Method and parameters:

Filters are applied per water year, each year's time series (sorted by day of water year) is filtered independently, so the filter re-initializes every Oct 1 rather than running continuously over the record.

BFI = Σ(daily baseflow) / Σ(daily Q) within each water year, constrained to [0, 1]; a year with total Q ≤ 0 is skipped.

Eckhardt filter: BFImax = 0.8 (maximum baseflow index), a = 0.98 (recession constant). Single forward pass; initialized at baseflow₁ = BFImax · Q₁, and constrained to [0, Q].

Lyne-Hollick filter: a= 0.925, passes = 2 (each pass is a forward sweep followed by a backward sweep, so 2 passes = 4 sweeps total). Quickflow constrained to [0, Q] in both directions; baseflow = Q − quickflow, must be ≥0.

Parameterized versions: gages with insufficient recession data (10 or fewer alpha pairs in the whole record), or a derived constant outside (0, 1), produce NA for all parameterized BFI values.

Requirements and decisions:

Linear interpolation of streamflow for data gaps less than 3 days; years with data gaps > 3 days are rejected by the centralized preprocessor, so baseflow filters never encounter gaps in valid years. No filter restart or gap counting is needed.

Each day of the baseflow time series is equal to or less than the daily streamflow.

BFI must be between 0 and 1.

### 3.4 Recession

Function: analyze_recession_parameters()

Purpose: Analyzes recession events to calculate parameters like log(a), b, and concavity, from the relationship −dQ/dt = a·Q^b.

Metrics:

log_a_pointcloud: Median log(a) from point cloud analysis of all recession events (b fixed at 1).

log_a_events: Median log(a) from individual recession events (b fixed at 1).

b_pointcloud: Median recession exponent b from point cloud analysis (free fit).

b_events: Median b from individual events (free fit).

concavity: Difference in b between the first and second halves of recession events.

n_recession_events: Count of recession events per water year (computed independently of the 25-event minimum below).

log_a_seasonality_amplitude_all: Seasonal amplitude of log(a) for all years (single value; also computed for the first and last halves of the record).

log_a_seasonality_minimum_all: Water year date of minimum log(a) for all years (single value; also computed for the first and last halves of the record).

Method and parameters:

Recession events are identified as the longest contiguous windows where both Q and |dQ/dt| are monotonically decreasing, with a minimum length of 5 days. The first day of each event is removed during power law fitting (standard practice — storm peak influence).

Parameter b is determined using the line of best fit to the log-log plot of −dQ/dt versus Q (free power-law fit).

Parameter a is determined with the exponent constrained to b = 1 (linear reservoir assumption, applied at all locations and periods): log(a) = median of [log(−dQ/dt) − log(Q)] — a median, not a regression intercept. Fixing b = 1 decouples trends and seasonality in a from changes in b.

Split recession events temporally (first vs. second half) for concavity analysis.

Fits sinusoidal models to analyze seasonal variation in log(a) (fit to the per-event b = 1 values).

All recession metrics except n_recession_events require at least 25 recession events across the record; gages below the minimum report NA.

The recession and parameterized-BFI families are exempt from the trend-completeness requirements and the 20-value statistics floor (Part 1.3).

Requirements and decisions:

Count the number of recession events in each year.

To control the quality of fitted parameters, calculate recession fits and create a flag for any R² < 0.8.

### 3.5 Pulses and Flow Reversals

Function: calculate_pulse_metrics()

Purpose: Calculates metrics for high/low flow pulses (frequency, duration) and flow reversals.

Metrics:

n_high_pulses_year / n_high_pulses_all: Number of high-flow pulse events per year. A high-pulse day is a day with flow above the 90th percentile of flow; _year and _all calculate percentiles from Q for the specified year or the period of record, respectively. A pulse event is a maximal run of consecutive pulse days.

n_low_pulses_year / n_low_pulses_all: Number of low-flow pulse events per year. A low-pulse day is a day with flow below the 10th percentile of flow; _year and _all as above. Zero-flow days are included in all calculations; on strongly intermittent streams the 10th-percentile threshold can be exactly 0, in which case zero-flow days do not count as low-pulse days because the comparison is strictly below the threshold.

dur_high_pulses_year / dur_high_pulses_all: Mean duration in days of high-flow pulses per year.

dur_low_pulses_year / dur_low_pulses_all: Mean duration in days of low-flow pulses per year.

TQmean: Percentage of days with flow above the annual mean.

Flow_Reversals_annual / _winter / _spring / _summer / _fall: Number of flow reversals (annual and seasonal). Only includes changes in flow direction exceeding 2% of current flow: a reversal day requires a sign change between the incoming and outgoing day-to-day flow differences, with the outgoing change larger than 2% of that day's flow.

Requirements and decisions:

This metric requires data from across the period of record; confirm that the code can actually calculate the 90th/10th percentiles from the full period of record before the annual metrics.

High/low flow pulses are defined using the 90th and 10th percentiles of flow (Q).

Calculate annual metrics using period-of-record percentiles (_all); the _year variants are retained as a complementary per-year sensitivity.

### 3.6 Flashiness

Function: analyze_flashiness_trends()

Purpose: Calculates the Richards-Baker (R-B) flashiness index by year, from day-to-day changes.

Metrics:

flashinessRB (R-B Index): Sum of absolute differences in consecutive daily flows divided by the sum of total flow: R-B = Σ|Q(i) − Q(i−1)| / ΣQ(i) (Baker et al. 2004). Dimensionless; higher values indicate flashier response.

### 3.7 Flow Timing

Function: analyze_flow_timing_trends()

Purpose: Analyzes the timing of flow events (day of maximum flow, cumulative flow percentiles). Day 1 = October 1.

Metrics:

Dxx_day: Day of water year when cumulative flow reaches a given percentile. Labeled as D1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99.

D25_to_D75: Count of days between 25% and 75% cumulative flow.

Dmax: Water year day of maximum flow.

### 3.8 Runoff Ratios

Function: analyze_Q_PPT_relationships() · Requires precipitation data

Purpose: Analyzes runoff ratios (streamflow divided by precipitation) and their trends.

Metrics:

annual_runoff_ratio: Annual water year total streamflow divided by annual precipitation (Q/P).

winter_runoff_ratio: Winter (Dec–Feb) total streamflow divided by winter precipitation.

spring_runoff_ratio: Spring (Mar–May) total streamflow divided by spring precipitation.

summer_runoff_ratio: Summer (Jun–Aug) total streamflow divided by summer precipitation.

fall_runoff_ratio: Fall (Sep–Nov) total streamflow divided by fall precipitation.

runoff_ratio_high_count: Count of years where the annual runoff ratio is greater than 2 (per-gage diagnostic scalar).

Requirements and decisions:

Precipitation (PPT) data must be available alongside streamflow (Q).

Same QAQC procedures for Q should apply to NAs (missing data) in precipitation; zero is okay in precipitation data.

Calculates separate runoff ratios for annual and seasonal periods (winter, spring, summer, fall).

Cases where annual PPT < 10 mm or seasonal PPT < 1 mm: the ratio is set to NA.

### 3.9 Streamflow Elasticity

Function: calculate_streamflow_elasticity() · Requires precipitation data

Purpose: Calculates streamflow elasticity, which measures how sensitive streamflow is to changes in precipitation. Based on Sawicz et al. (2011). Values around 1.0 indicate proportional response; values greater than 1 define the catchment as elastic (sensitive to precipitation change), less than 1 as inelastic.

Metrics:

elasticity_static (departure-from-mean, Sawicz et al. 2011): Overall catchment elasticity — median of annual elasticity values (single per-gage value, no statistics). For each qualifying water year i:

E_i = [ (Q_i − Q_mean) / (P_i − P_mean) ] / (Q_mean / P_mean)

where Q_i and P_i are that water year's total streamflow and precipitation (mm), and Q_mean and P_mean are the means over all qualifying years. Years where |P_i − P_mean| ≤ 0.1 mm are skipped. Then elasticity_static = median(E_i).

elasticity_rolling: Rolling-window elasticity, used to detect trends in catchment sensitivity over time. Same formula, but computed within each window of 11 consecutive qualifying-year observations, using the window's own means (Q_mean_w, P_mean_w):

E_j = [ (Q_j − Q_mean_w) / (P_j − P_mean_w) ] / (Q_mean_w / P_mean_w)

elasticity_rolling(t) = median of the E_j within the window ending in year t.

elasticity_annual (year-over-year differences): Measures year-to-year sensitivity (accounting for antecedent conditions) rather than deviation-from-mean sensitivity:

E_t = [ (Q_t − Q_t−1) / (P_t − P_t−1) ] / (Q_mean / P_mean)

where the differences are between adjacent qualifying years, Q_mean / P_mean are the whole-record means, and the value is assigned to the later year t. Pairs with |P_t − P_t−1| ≤ 0.1 mm are skipped.

elasticity_years_total / elasticity_years_low_ppt: Diagnostic scalars — the number of years available for the gage, and the number excluded for annual PPT below 10 mm.

Requirements and decisions:

Requires daily streamflow (Q) and precipitation (PPT) data; calculating differences between years, so a minimum of two years of "good data" is needed for any pair, and a minimum of 15 years of valid data for the elasticity calculation overall.

Years with annual precipitation < 10 mm are flagged and excluded from analysis.

The year-over-year calculation requires both years of each pair to pass the minimum data requirements (consecutive qualifying years).

Year qualification is handled by the centralized preprocessor (config).

Count the number of years each site doesn't have sufficient data within a year, and what the number of excluded years would be if the threshold were < 30% data missing — added as documentation, not a filter.

Output: elasticity_static is a single value; elasticity_rolling and elasticity_annual each carry the 16 standard statistics. The elasticity family is exempt from the trend-completeness requirements and the 20-value statistics floor (Part 1.3).

### 3.10 Q-P Seasonality

Function: calculate_qp_seasonality() · Requires precipitation data

Purpose: Quantifies the seasonality in the relationship between cumulative streamflow (Q) and cumulative precipitation (P). Based on Wrede et al. (2015); however, the rolling window approach described below may differ from the original publication.

Metrics:

qp_slope_sd: Standard deviation of the 12 monthly mean cumulative Q-P slopes. Higher values indicate stronger seasonal variation in the streamflow-precipitation relationship.

qp_bimodality: Bimodality coefficient of the Q-P slope distribution = (skewness² + 1) / kurtosis. Values > 0.555 suggest bimodal or strongly seasonal patterns.

Method and parameters:

Requires daily streamflow (Q), precipitation (PPT), month, and day-of-water-year columns; year qualification is handled by the centralized preprocessor.

Calculates a 30-day rolling slope (backward-looking) of cumulative Q vs cumulative P for each water year, then aggregates the rolling slopes to 12 monthly mean values per year.

### 3.11 Negative Flow Days

Function: calculate_negative_days()

Metrics:

Negative_ann: The number of days with negative flow values per gage per year. Negative-flow days are retained in the data by default (reject_negative_flow = false, Part 1.2) — this metric is the accounting for them.

### 3.12 Snow

Function: calculate_snow_metrics() · Requires SWE data (Daymet)

Purpose: Fourteen metrics are computed from daily snow water equivalent (SWE; Daymet V4) for water years 1980–2023. All metrics operate on a thresholded series SWE* = SWE if SWE ≥ 10 mm, else SWE*=0. Days with SWE below 10 mm are treated as snow-free for both durations and magnitudes. A "snow day" is a day with SWE* > 0. A "spell" is a maximal run of consecutive snow days and the "anchor spell" is the spell containing the annual SWE maximum. If the annual SWE maximum falls across more than one day, then the anchor spell will correspond to the period with the first day with the annual SWE maximum. Gages without Daymet coverage carry NA for all snow metrics.

Metrics:

swe_max: Maximum daily SWE in the water year (mm). A snow-free year reports a valid 0.

swe_max_dowy: Day of water year of the SWE peak. NA in snow-free years.

snow_cover_days: Number of days with SWE ≥ 10 mm. A snow-free year reports a valid 0.

snow_on_dowy: First day of the anchor spell.

snow_off_dowy: First snow-free day after the anchor spell.

melt_season_days: Number of days during the melt season (or between peak SWE and no snow. Calculated as: snow_off_dowy − swe_max_dowy.

melt_rate: rate of snow melt during the melt season: swe_max divided by melt_season_days (mm/day; net ablation rate between peak and snow-off).

ssm: Snow seasonality metric = (seasonal snow days − ephemeral snow days) / total snow days, where a seasonal spell lasts at least 60 continuous days (Hatchett 2021; Petersky & Harpold 2018). Ranges −1 (fully ephemeral) to +1 (fully seasonal).

swe_apr1: SWE on calendar April 1.

melt_before_peak: Total melt (sum of daily SWE decreases) occurring before the peak day (mm).

melt_before_peak_pct: melt_before_peak as a percentage of total water-year melt.

melt_before_peak_to_max_swe: Ratio of total melt to maximum SWE: melt_before_peak divided by swe_max.

melt_com_dowy: Day of water year when cumulative melt reaches 50% of the water-year total.

swe_max_to_ppt: swe_max divided by total water-year precipitation (years with precipitation ≤ 10 mm return NA).

Requirements and decisions:

snow_on_dowy and snow_off_dowy are censored (NA) when the anchor spell touches Oct 1 or Sep 30, since the snowpack predates or outlasts the water year.

Timing and melt metrics are NA in snow-free years; magnitude metrics (swe_max, snow_cover_days, swe_apr1, swe_max_to_ppt) report valid zeros.

In addition to the standard trend-completeness requirements, the ten timing/melt/regime metrics require at least 80% of SWE-valid years in both the first and last decade of the gage's SWE record to be snowy before trend statistics are reported (prevents trends conditioned only on snow-present years at gages whose snow is disappearing or appearing).

### 3.13 Streamflow Drought

Function: calculate_drought_metrics()

Purpose: Drought duration and deficit against fixed percentile thresholds calculated per-water-year, following Adelsperger et al. (in review). Daily flow is first smoothed with a 7-day centered moving average, applied within continuous runs of consecutive dates only (the window never averages across a gap left by a rejected year). Thresholds are magnitude percentiles of the smoothed flow pooled over the gage's whole record, computed with the unbiased Weibull plotting position p = i/(n+1) (Laaha et al. 2017); a gage needs at least 10 qualifying years to receive thresholds. Drought severity levels mirror the U.S. Drought Monitor classes: p30 (D0, abnormally dry), p20 (D1), p10 (D2), p5 (D3), p2 (D4, exceptional).

Metrics:

drought_duration_fixed_p{2,5,10,20,30}: Number of days in the water year with smoothed flow strictly below the threshold (days).

drought_deficit_fixed_p{2,5,10,20,30}: Sum of (threshold − smoothed flow) over those days — the magnitude-weighted counterpart to duration, separating long shallow droughts from short severe ones. Units: mm for area-normalized gages (raw flow units otherwise; see notes).

drought_threshold_fixed_p{2,5,10,20,30}: The threshold values themselves (mm/day for area-normalized gages; one per gage, no statistics), reported so thresholds are auditable.

Requirements and decisions:

A year with no sub-threshold days reports a valid 0, not NA.

Thresholds come from each analysis window's own record, so drought values are comparable within a product but must not be compared across the WY 1993–2025 and WY 1980–2025 products.

The five levels form a severity ladder — duration and deficit are non-decreasing from p2 up to p30 by construction — not five independent measurements.

On strongly intermittent streams the low-level thresholds can be exactly 0; the strict below-threshold comparison then reports 0 duration and 0 deficit every year. Check drought_threshold_fixed_p{n} before reading zeros as "no drought".

drought_duration_fixed_p10 is highly correlated with the low-pulse metrics (both use a 10th-percentile threshold) and should not be presented alongside them as independent evidence.

Units: the duration metrics are scale-invariant and valid for all gages; deficits and thresholds carry flow units, so for the small set of un-normalized gages (area_normalized = false) they are in raw m³/s-based units and are not comparable across gages.

## Part 4: Data Ingestion and Utility Functions

process_gages_rawData()

Purpose: Processes raw gage data (USGS/Canadian) to calculate streamflow metrics for individual gages and save results to a file. Data are available up to the present day, but no climate data are directly available.

Requirements and decisions:

A folder of USGS and HYDAT (Canadian) metadata must be stored locally in order to process watersheds.

A minimum number of valid years (min_num_years) is required for a gage to be included.

Legacy note: the minimum-flow filter (min_Q_value_and_days) applied only to the legacy processing path and was removed from the current pipeline in April 2026 — the centralized preprocessor is now the single source of year qualification.

Missing or invalid years are excluded from analysis.

process_caravan_gages()

Purpose: Processes Caravan NetCDF timeseries data for watersheds, calculates streamflow metrics, and saves results. Coincident climate data are generally available, but the record of streamflow is truncated to as early as the early 2000s (as late as 2018 for HYSETS).

Requirements and decisions:

Minimum years of valid data (min_num_years) are required for processing.

Filters years based on minimum flow thresholds and valid days.

Handles redundancy between datasets (e.g., CAMELS and HYSETS).

generate_streamflow_dt() and generate_streamflow_dt_caravan()

Purpose: Converts raw gage data (USGS/Canadian) or Caravan NetCDF data into a standardized streamflow data.table.

Requirements and decisions:

Streamflow is provided in appropriate units (or converted, e.g., m³/s to mm/day).

Filters data by specified date ranges and required years of data.

Handles missing or invalid values (e.g., flagged data).

Adds derived columns (e.g., year, month, day of year).

integrate_daymet_with_streamflow()

Purpose: Joins Daymet climate data (precipitation, temperature) to streamflow data based on date matching.

Requirements and decisions:

Requires a pre-processed Daymet parquet file (data_out/daymet_1980_2023.parquet).

Matches streamflow gage IDs to Daymet site IDs.

Joins climate variables (prcp, tmin, tmax, swe, vp, srad) to streamflow data.table on Date.

Reports coverage percentage; warns if < 95% of dates have matching climate data.

Renames prcp column to PPT for compatibility with existing Q-PPT analysis functions.

## Part 5: Automated Data Quality Flags

All checks run on per-gage summary values (the _mean statistic, except where noted). Thresholds are config-driven (qa_qc section of signatures_config.json) and identical across the Julia, Python, and R implementations. A value of true marks a potential quality issue for review; flags never remove a gage from the downstream analysis (data are flagged only), and a missing input value never triggers a flag.

flagged_for_qann_range: Qann_mean outside [0, 2000] mm

flagged_for_bfi_eckhardt_range: BFI_Eckhardt_mean outside [0, 1]. Defensive check — per-year BFI values are clamped to [0, 1] during calculation, so this flag cannot fire under current code.

flagged_for_bfi_lynehollick_range: BFI_LyneHollick_mean outside [0, 1]. Same defensive note as above.

flagged_for_flashiness_range: flashinessRB_mean outside [0, 2]

flagged_for_tqmean_range: TQmean_mean outside [0, 100]

flagged_for_d50_range: D50_day_mean outside [1, 366]

flagged_for_elasticity_range: elasticity_static outside [0.1, 5] (the per-gage scalar — the only flag not based on a _mean column)

flagged_for_runoff_ratio_range: annual_runoff_ratio_mean outside [0.01, 1.5]

flagged_for_seasonal_sum: the sum of the four seasonal totals (Qwin_mean + Qspr_mean + Qsum_mean + Qfal_mean) deviates from Qann_mean by more than 20%

flagged_for_percentile_order: non-decreasing order Q5 ≤ Q25 ≤ Q50 ≤ Q75 ≤ Q95 violated (on _mean values; ties are allowed and do not flag)

flagged_for_timing_order: non-decreasing order D5 ≤ D50 ≤ D95 violated (on _day_mean values; ties are allowed and do not flag)

flagged_for_high_na: more than 30% of a gage's numeric output columns are NA — often indicates missing climate/SWE coverage rather than bad streamflow data

## Part 6: References

Adelsperger, S., et al. (in review). A novel severity-based approach for assessing streamflow drought characteristics and drivers.

Baker, D.B., Richards, R.P., Loftus, T.T., & Kramer, J.W. (2004). A new flashiness index: characteristics and applications to midwestern rivers and streams. Journal of the American Water Resources Association, 40(2), 503–522.

Collischonn, W., & Fan, F.M. (2013). Defining parameters for Eckhardt's digital baseflow filter. Hydrological Processes, 27(18), 2614–2622.

Eckhardt, K. (2005). How to construct recursive digital filters for baseflow separation. Hydrological Processes, 19(2), 507–515.

Hatchett, B.J. (2021). Seasonal and Ephemeral Snowpacks of the Conterminous United States. Hydrology, 8(1), 32.

Laaha, G., et al. (2017). Unbiased plotting positions for low-flow frequency analysis.

Lyne, V., & Hollick, M. (1979). Stochastic time-variable rainfall-runoff modelling. Institute of Engineers Australia National Conference, 89–93.

Peters, N.E., & Aulenbach, B.T. (2011). Water storage at the Panola Mountain Research Watershed, Georgia, USA. Hydrological Processes, 25(25), 3878–3889.

Petersky, R., & Harpold, A. (2018). Now you see it, now you don't: a case study of ephemeral snowpacks and soil moisture response in the Great Basin, USA. Hydrology and Earth System Sciences, 22, 4891–4906.

Pettitt, A.N. (1979). A non-parametric approach to the change-point problem. Applied Statistics, 28(2), 126–135.

Sawicz, K., Wagener, T., Sivapalan, M., Troch, P. A., & Carrillo, G. (2011). Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA. Hydrology and Earth System Sciences, 15(9), 2895–2911.

Wrede, S., Fenicia, F., Martinez-Carreras, N., Juilleret, J., Hissler, C., Krein, A., ... & Pfister, L. (2015). Towards more systematic perceptual model development: a case study using 3 Luxembourgish catchments. Hydrological Processes, 29(12), 2731–2750.

Yilmaz, K.K., Gupta, H.V., & Wagener, T. (2008). A process-based diagnostic approach to model evaluation: application to the NWS distributed hydrologic model. Water Resources Research, 44(9), W09417.

---

## LEGACY DOCUMENT (still trailing the new one in the Google Doc; unchanged since 2026-08-31)

# Summary Documentation for Streamflow Signatures

github repo here or here

Summary data audit sheet here

NA Handling (for the whole analysis):

Do not replace NAs with zeros; streamflow that is zero is okay

NAs can be for various reasons, we need to use the USGS codes that describe what creates those NAs (e.g. Ice affected, etc.) for each site that has ice effected days we need a total count of the number of days that are ice affected

If there is up to 3 days of continuous missing streamflow data, interpolate between the points; tell us if there is a data driven answer to the number of days that would not create problems for trend analysis

Following the USGS data release QA criteria, four conditions are screened per water year: (i) more than 3 consecutive days of missing data, (ii) more than 30 total missing days in the year, (iii) negative flow values, and (iv) constant monthly standard deviations during periods of non-zero streamflow (indicative of data errors or highly controlled discharge). Items i and ii REMOVE the year: any gap longer than 3 consecutive days, or more than 30 missing days in total, excludes that water year from all calculations. Internal gaps of up to 3 days are filled by linear interpolation for years that pass, but interpolated days still count toward the 30-missing-day limit (the count is taken on the raw data, before interpolation). Items iii and iv FLAG but do not remove: negative values are retained by default (config option reject_negative_flow, default false) and are counted in Negative_ann; a constant monthly standard deviation sets a QA flag only and never rejects a year.

Flag years with constant monthly standard deviations during periods of non-zero streamflow indicative highly controlled discharges.

Negative_ann: Calculate the number of days with negative values per gage per year

For any seasonal or annual metrics the period must have at least 80% of the data.

Seasonal periods are defined as:

Winter: December-February

Spring: March-May

Summer: June-August

Fall: September-November

All annual metrics should be calculated on the water year October 1 - September 30

For a trend to be calculated, at least 80% of the first decade must be complete, AND at least 80% of the last decade of all trend periods AND at least 60% of the entire annual streamflow metric time series must be complete

## 1. Glossary of Streamflow Signatures

Flow Volume Metrics

Qann: Annual total streamflow in mm.

Qwin: Winter total streamflow in mm.

Qspr: Spring total streamflow in mm.

Qsum: Summer total streamflow in mm.

Qfal: Fall total streamflow in mm.

Qxx: Flow percentiles (e.g., Q10, Q50, Q90) in mm.

Baseflow Metrics

Filters are applied per water year, independently — each year's daily series (sorted by day of water year) is filtered on its own, so the filter re-initializes every Oct 1 rather than running continuously over the record.

BFI itself = Σ(daily baseflow) / Σ(daily Q) within each water year, clamped to [0, 1]; a year with total Q ≤ 0 is skipped.

BFI_Eckhardt: Baseflow index calculated using the Eckhardt filter:

BFImax = 0.8 (maximum baseflow index)

a = 0.98 (filter/recession constant)

Single forward pass; initialized at baseflow₁ = BFImax · Q₁, and constrained to 0 ≤ baseflow ≤ Q at every step

BFI_LyneHollick: Baseflow index calculated using the Lyne-Hollick filter.

alpha = 0.925 (filter parameter)

passes = 2 (each pass is a forward sweep followed by a backward sweep, so 2 passes = 4 sweeps total)

Quickflow constrained to [0, Q] in both sweep directions; baseflow = Q − quickflow, floored at 0

Flow Duration Curve (FDC) Metrics

FDCall: Slope of the annual flow duration curve across the full range of exceedance probabilities.

FDC90th: Slope of the annual flow duration curve over the low-flow tail (exceedance probability ≥ 0.90).

FDCmid: Slope of the annual flow duration curve over the midsegment (exceedance probability 0.20-0.80).

________________

Recession Metrics

log_a_pointcloud: Median log(a) from point cloud analysis of all recession events.

log_a_events: Median log(a) from individual recession events.

b_pointcloud: Median recession exponent b from point cloud analysis.

b_events: Median b from individual events.

concavity: Difference in b between the first and second halves of recession events.

log_a_seasonality_amplitude_all: Seasonal amplitude of log(a) for all years.

log_a_seasonality_minimum_all: water year date of minimum log(a) for all years.

Pulse Metrics

n_high_pulses_year: Number of high-flow pulses per year:

A high-pulse day is defined as a day with flow above the 90th percentile of flow; *_year and *_all calculate percentiles based on Q for the specified year or the period of record, respectively

n_low_pulses_year: Number of low-flow pulses per year:

A low-pulse day is defined as a day with flow below the 10th percentile of flow; *_year and *_all calculate percentiles based on Q for the specified year or the period of record, respectively. Note: for intermittent streams, 0 flow days are excluded

dur_high_pulses_year: Mean duration in days of high-flow pulses per year.

dur_low_pulses_year: Mean duration in days of low-flow pulses per year.

TQmean: Percentage of days with flow above the annual mean.

Flow_Reversals: Number of flow reversals (annual or seasonal). Only included changes in flow direction exceeding 2% of current flow.

Flashiness Metrics

R-B Index: Richards-Baker flashiness index calculated as the sum of absolute differences in consecutive flows divided by the sum of total flow.

Flow Timing Metrics

Dxx_day: water year day when cumulative flow reaches a given percentile (e.g., D10, D50, D90).

D25_to_D75: Count of days between 25% and 75% cumulative flow.

Dmax: water year day of maximum flow.

Elasticity Metrics:

elasticity_static (departure-from-mean, Sawicz et al. 2011): Overall catchment elasticity - median of annual elasticity values. Measures how sensitively streamflow responds to precipitation changes. Values around 1.0 indicate proportional response; >1 indicates amplified response.

For each qualifying water year i:

E_i = [ (Q_i − Q_mean) / (P_i − P_mean) ] / (Q_mean / P_mean)

where Q_i and P_i are that water year's total streamflow and precipitation (mm), and Q_mean and P_mean are the means over all qualifying years. Years where |P_i − P_mean| ≤ 0.1 mm are skipped. Then:

elasticity_static = median(E_i)

elasticity_rolling: Rolling window (11-year) of elasticity values, used to detect trends in catchment sensitivity over time.

Same formula, but computed within each 11-year window using the window's own means (Q_mean_w, P_mean_w):

E_j = [ (Q_j − Q_mean_w) / (P_j − P_mean_w) ] / (Q_mean_w / P_mean_w)

elasticity_rolling(t) = median of the E_j within the window ending in year t

elasticity_annual: (year-over-year differences)

E_t = [ (Q_t − Q_t−1) / (P_t − P_t−1) ] / (Q_mean / P_mean)

where Q_t − Q_t−1 and P_t − P_t−1 are differences between adjacent qualifying years, Q_mean / P_mean are the whole-record means, and the value is assigned to the later year t. Pairs with |P_t − P_t−1| ≤ 0.1 mm are skipped.

Runoff Ratio Metrics:

annual_runoff_ratio: Annual water year streamflow divided by annual precipitation (Q/P).

winter_runoff_ratio: Winter (Dec-Feb) total streamflow divided by winter precipitation.

spring_runoff_ratio: Spring (Mar-May) total streamflow divided by spring precipitation.

summer_runoff_ratio: Summer (Jun-Aug) total streamflow divided by summer precipitation.

fall_runoff_ratio: Fall (Sep-Nov) total streamflow divided by fall precipitation.

Storage Metric:

avg_storage: Mean annual catchment storage (mm), from the simplified water balance S = cumulative sum of (P − Q) within each water year; annual storage is interpolated at mean discharge and averaged across years (Peters & Aulenbach 2011). Computed only for gages with valid precipitation data and area-normalized flow (NA otherwise), and included in the data product, but OMITTED FROM MAJOR ANALYSES: the water balance carries no evapotranspiration term, so values overestimate storage wherever ET losses are significant. Retained in the output for completeness; interpret with caution.

Snow Metrics

Fourteen per-water-year metrics computed from daily snow water equivalent (SWE; Daymet V4, a model product, calendar years 1980–2023). All metrics operate on a thresholded series SWE* = SWE if SWE ≥ 10 mm, else 0; days below 10 mm are treated as snow-free for both durations and magnitudes. A "snow day" is a day with SWE* > 0; a "spell" is a maximal run of consecutive snow days; the "anchor spell" is the spell containing the annual SWE maximum (ties resolved to the first day). Gages without Daymet coverage carry NA for all snow metrics.

swe_max: Maximum daily SWE in the water year (mm). A snow-free year reports a valid 0.

swe_max_dowy: Day of water year of the SWE peak. NA in snow-free years.

snow_cover_days: Number of days with SWE ≥ 10 mm. A snow-free year reports a valid 0.

snow_on_dowy: First day of the anchor spell.

snow_off_dowy: First snow-free day after the anchor spell.

melt_season_days: snow_off_dowy − swe_max_dowy.

melt_rate: swe_max divided by melt_season_days (mm/day; net ablation rate between peak and snow-off).

ssm: Snow seasonality metric = (seasonal snow days − ephemeral snow days) / total snow days, where a seasonal spell lasts at least 60 continuous days (Hatchett 2021; Petersky & Harpold 2018). Ranges −1 (fully ephemeral) to +1 (fully seasonal).

swe_apr1: SWE on calendar April 1 (leap-year safe).

melt_before_peak: Total melt (sum of daily SWE decreases) occurring before the peak day (mm).

melt_before_peak_pct: melt_before_peak as a percentage of total water-year melt.

melt_before_peak_to_max_swe: melt_before_peak divided by swe_max.

melt_com_dowy: Day of water year when cumulative melt reaches 50% of the water-year total.

swe_max_to_ppt: swe_max divided by total water-year precipitation (years with PPT ≤ 10 mm return NA).

Notes: snow_on_dowy and snow_off_dowy are censored (NA) when the anchor spell touches Oct 1 or Sep 30, since the snowpack predates or outlasts the water year. Timing and melt metrics are NA in snow-free years; magnitude metrics (swe_max, snow_cover_days, swe_apr1, swe_max_to_ppt) report valid zeros. In addition to the standard trend-completeness requirements, the ten timing/melt/regime metrics require at least 80% of SWE-valid years in both the first and last decade of the gage's SWE record to be snowy before trend statistics are reported (prevents trends conditioned only on snow-present years at gages whose snow is disappearing or appearing).

Streamflow Drought Metrics

Per-water-year drought duration and deficit against fixed percentile thresholds. Daily flow is first smoothed with a 7-day centered moving average, applied within continuous runs of consecutive dates only (the window never averages across a gap left by a rejected year). Thresholds are magnitude percentiles of the smoothed flow pooled over the gage's whole record, computed with the unbiased Weibull plotting position p = i/(n+1) (Laaha et al. 2017). Five severity levels mirror the U.S. Drought Monitor classes: p30 (D0, abnormally dry), p20 (D1), p10 (D2), p5 (D3), p2 (D4, exceptional).

drought_duration_fixed_p{2,5,10,20,30}: Number of days in the water year with smoothed flow strictly below the threshold (days).

drought_deficit_fixed_p{2,5,10,20,30}: Sum of (threshold − smoothed flow) over those days (mm) — the magnitude-weighted counterpart to duration, separating long shallow droughts from short severe ones.

drought_threshold_fixed_p{2,5,10,20,30}: The threshold values themselves (mm/day; one per gage), reported so thresholds are auditable.

Notes: a year with no sub-threshold days reports a valid 0, not NA. Thresholds come from each analysis window's own record, so drought values are comparable within a product but must not be compared across the WY 1993–2025 and WY 1980–2025 products.ß

## 2. Glossary of Statistical Metrics

The following 16 statistics are calculated for every streamflow signature metric (except where otherwise noted):

senn_slp: Theil-Sen slope — robust non-parametric trend estimate (median of pairwise slopes), in signature units per year.

linear_slp: Ordinary least-squares linear regression slope, in signature units per year; parametric companion to the Theil-Sen slope.

spearman_rho: Spearman's rank correlation between the annual values and water year — non-parametric measure of monotonic trend strength (−1 to +1).

spearman_pval: P-value of the Spearman correlation.

mk_rho: Mann-Kendall tau (tie-corrected tau-b) — non-parametric measure of monotonic trend strength (−1 to +1).

mk_pval: P-value of the Mann-Kendall test (continuity-corrected normal approximation with tie-corrected variance).

mean: Arithmetic mean of the signature's annual values across all qualifying water years.

median: Median of the signature's annual values across all qualifying water years.

pettitt_cp_year: Most likely changepoint year. Always reported when the test runs; use pettitt_pval to judge significance.

pettitt_pval: Asymptotic p-value of the Pettitt test (< 0.05 indicates a significant changepoint).

pettitt_pre_mean: Mean of the annual values before the changepoint year.

pettitt_post_mean: Mean of the annual values after the changepoint year.

pettitt_delta_mean: pettitt_post_mean − pettitt_pre_mean (signature units).

pettitt_pct_change: pettitt_delta_mean as a percentage of |pettitt_pre_mean|.

pettitt_pre_mk_pval: Mann-Kendall p-value for the pre-changepoint segment (tests for a residual trend within the segment).

pettitt_post_mk_pval: Mann-Kendall p-value for the post-changepoint segment.

## 3. Summary of Functions

process_gages_rawData()

Purpose: Processes raw gage data (USGS/Canadian) to calculate streamflow metrics for individual gages and save results to a file. Data are available up to the present day, but no climate data are directly available.

Requirements and Decisions:

A folder of USGS and HYDAT (Canadian) metadata must be stored locally in order to process watersheds.

A minimum number of valid years (min_num_years) is required for a gage to be included.

Data is filtered by a minimum flow threshold (min_Q_value_and_days). A watershed with fewer than the minimum days of > min_Q_value is not processed, but if the threshold is surpassed then non-zero flow is still passed to subsequent functions for processing.

Missing or invalid years are excluded from analysis.

process_caravan_gages()

Purpose: Processes Caravan NetCDF timeseries data for watersheds, calculates streamflow metrics, and saves results. Coincident climate data are available, but the record of streamflow is truncated to as early as the early 2000s (as late as 2018 for HYSETS).

Requirements and Decisions:

Minimum years of valid data (min_num_years) are required for processing.

Filters years based on minimum flow thresholds and valid days.

Handles redundancy between datasets (e.g., CAMELS and HYSETS).

generate_streamflow_dt() and generate_streamflow_dt_caravan()

Purpose: Converts raw gage data (USGS/Canadian) or Caravan NetCDF data into a standardized streamflow data.table.

Requirements and Decisions:

Streamflow is provided in appropriate units (or converted, e.g., m³/s to mm/day).

Filters data by specified date ranges and required years of data.

Handles missing or invalid values (e.g., flagged data).

Adds derived columns (e.g., year, month, day of year).

calculate_flow_vols_by_year()

Purpose: Calculates annual and seasonal flow volumes, as well as flow percentiles (e.g., Q1, Q10, Q50, Q90).

Requirements and Decisions:

Requires year, Q, month, and day of year columns in input data.

Frozen days are okay (USGS data flag and Hydat data flag requirement); 0 flow days are okay.

Any season with less than 80% of the data should not be included, count the number of years where the season does not meet that threshold.

analyze_baseflow_indices()

Purpose: Calculates baseflow indices (BFI) using Eckhardt and Lyne-Hollick digital filters.

Requirements and Decisions:

T-1 is part of the moving window that the Lyne-Hollick uses for baseflow separation. Data continuity could be an issue with the way it is calculated.

Linear interpolation of streamflow for data gaps less than 3days;

Years with data gaps >3 days are rejected by the centralized preprocessor, so baseflow filters never encounter gaps in valid years. No filter restart or gap counting is needed.

Default parameters for Eckhardt filter: BFImax = 0.8, a = 0.98.

Once this has been run, make sure each day of the baseflow time series is equal to or less than the daily streamflow.

BFI must be between 0 and 1.

analyze_baseflow_incidies_with_parameters()

Purpose: Calculate baseflow indices using the recession parameters below as inputs to (BFI) using Eckhardt and Lyne-Hollick digital filters.

Requirements and Decisions: Same as the baseflow index decisions above.

analyze_recession_parameters()

Purpose: Analyzes recession events to calculate parameters like log(a), b, and concavity.

Requirements and Decisions:

Recession events are defined as periods of at least 5 consecutive days with monotonic decreases in flow (Q) and its derivative (|dQ/dt|).

Recession events are identified as the longest contiguous windows where both Q and |dQ/dt| are monotonically decreasing, with a minimum length of 5 days. The first day of each event is removed during power law fitting (standard practice — storm peak influence).

Split recession events temporally (first vs. second half) for concavity analysis.

Fits sinusoidal models to analyze seasonal variation in log(a).

Count the number of recession events in each year.

Parameter b is determined using the line of best ﬁt to the log-log plot of −dq∕dt versus q; parameter similarly determined, but constrained with b = 1.

To control the quality of ﬁtted parameters, calculate recession ﬁts and create a flag for any R2 < 0.8.

analyze_flashiness()

Purpose: Calculates the Richards-Baker (RB) flashiness index by year.

Uses day to day changes.

calculate_pulse_metrics()

Purpose: Calculates metrics for high/low flow pulses (e.g., frequency, duration) and flow reversals.

Requirements and Decisions:

This metric requires data from across the period of record, confirm that the code can actually calculate the 90th/10th percentiles from the full period of record before the annual metrics.

High/low flow pulses are defined using the 90th and 10th percentiles of flow (Q).

Calculate annual metrics using period-of-record percentiles.

analyze_flow_timing()

Purpose: Analyzes the timing of flow events (e.g., day of maximum flow, cumulative flow percentiles).

Requirements and Decisions:

Calculate day of water year for cumulative flow percentiles, and the day with maximum flow.

Dxx_day: Day of water year when cumulative flow reaches a given percentile (e.g., D10, D50, D90). Labeled as D1,5,10,20,30,40,50,60,70,80,90,95,99.

D25_to_D75: Days between 25% and 75% cumulative flow.

Dmax: Water year day of maximum flow.

analyze_runoff_ratios()

Purpose: Analyzes runoff ratios (streamflow divided by precipitation) and their trends.

Runoff Ratio Definition:

Runoff_ratio: calculated as annual or seasonal periods with total streamflow divided by total precipitation (Q/P).

Requirements and Decisions:

Precipitation (PPT) data must be available alongside streamflow (Q).

Same QAQC procedures for Q should apply to NAs (missing data) in precipitation; zero is okay in precipitation data.

Calculates separate runoff ratios for annual and seasonal periods (winter, spring, summer, fall).

Cases where annual PPT < 10mm or seasonal PPT < 1mm the ratio is set to NA.

Runoff_ratio_high: flag runoff ratios that are greater than two.

calculate_streamflow_elasticity()

Purpose: Calculates streamflow elasticity, which measures how sensitive streamflow is to changes in precipitation. Based on Sawicz et al. (2011).

Definitions:

elasticity_static: Overall catchment elasticity - median of annual elasticity values. Measures how sensitively streamflow responds to precipitation changes. Values ~1.0 indicate proportional response; >1 indicates amplified response.

elasticity_rolling: Rolling window (11-year) elasticity values, used to detect trends in catchment sensitivity over time.

elasticity_annual: measures year-to-year sensitivity rather than deviation-from-mean sensitivity

Requirements and Decisions:

Requires daily streamflow (Q) and precipitation (PPT) data. Calculating differences between years, so need minimum of two years of “good data.”

Minimum of 15 years of valid data required for elasticity calculation.

Flag years with annual precipitation <10mm and exclude from analysis.

Calculates both static elasticity (median of all annual values) and rolling window elasticity (11-year window by default).

Annual elasticity formula: E = (dQ/dP) / (Q_mean/P_mean), where dQ = Q_year - Q_mean and dP = P_year - P_mean. This is like an anomaly direction, while the approach from the paper is accounting for anomaly (lets create a second calculation to do that) – achieved as of 4/16

Distance from average v.s. Accounting for antecedent conditions

Values ~1.0 indicate proportional response; values >1 indicate amplified streamflow response to precipitation changes.

From Sawicz 2011 slightly different interpretation from above based on dq= t1-t0 : A value of 1 indicates that a 1 % precipitation change leads to a 1 % change in streamflow. A value greater or less than 1 would, respectively, define the catchment as being elastic, i.e., sensitive to change of precipitation, or inelastic, i.e., insensitive to a change of precipitation.

To proceed with the additional calculation of using the previous years precip/Q then you would need to pass the minimum data requirements for both years (consecutive years)

Year qualification handled by centralized preprocessor (config).

Count the number of years each site doesn’t have sufficient data within a year, and what the number of excluded years would be if the value was <30% data missing – added as documentation not filter

Output: Returns ‘elasticity_static` (single value) plus the 8 standard metric statistics (mean, median, linear trend, theil sen, p value); ‘elasticity_rolling’ (one value per 11 year windows); ‘elasticity_annual’ plus the 8 standard stats

calculate_qp_seasonality()

Purpose: Quantifies the seasonality in the relationship between cumulative streamflow (Q) and cumulative precipitation (P). Based on Wrede et al. (2015), however the rolling window approach described below may differ from the original publication.

Q-P Seasonality Metric Definitions

qp_slope_sd: Standard deviation of monthly cumulative Q-P slopes. Higher values indicate stronger seasonal variation in the streamflow-precipitation relationship.

qp_bimodality: Bimodality coefficient of the Q-P slope distribution. Values >0.555 suggest bimodal or strongly seasonal patterns.

Requirements and Decisions:

Requires daily streamflow (Q), precipitation (PPT), month, and day-of-water-year (dowy) columns.

Year qualification handled by centralized preprocessor (config).

Calculates 30-day rolling slope (backward-looking) of cumulative Q vs cumulative P for each water year.

Aggregates rolling slopes to 12 monthly mean values per year.

Two metrics are calculated per year:

qp_slope_sd: Standard deviation of the 12 monthly slopes (higher values = more seasonal variation).

qp_bimodality: Bimodality coefficient = (skewness² + 1) / kurtosis. Values >0.555 suggest bimodal/seasonal patterns across all water years.

## Part 2: utility functions

integrate_daymet_with_streamflow()

Purpose: Joins Daymet climate data (precipitation, temperature) to streamflow data based on date matching.

Requirements and Decisions:

Requires a pre-processed Daymet parquet file (`data_out/daymet_1980_2023.parquet`).

Matches streamflow gage IDs to Daymet site IDs.

Joins climate variables (PPT, tmin, tmax, swe, vp, srad) to streamflow data.table on Date.

Reports coverage percentage; warns if <95% of dates have matching climate data.

Renames `prcp` column to `PPT` for compatibility with existing Q-PPT analysis functions.

convert_daymet_zip_to_parquet()

Purpose: One-time conversion of Daymet ZIP archive (containing 44 annual CSV files) to a single optimized parquet file.

Requirements and Decisions:

Input: ZIP file containing daymet_1980.csv through daymet_2023.csv.

Reconstructs the Date column from year, month, and row order (original CSVs lack a day column).

Validates that each site has 365 or 366 days per year.

Uses snappy compression for efficient storage (~70% reduction vs uncompressed).

Output columns: site_id, Date, prcp, tmin, tmax, swe, vp, srad.

Output:

Single parquet file (~3.9 GB) containing ~98 million rows for 6,087 sites.

## Part 3: References

Sawicz, K., Wagener, T., Sivapalan, M., Troch, P. A., & Carrillo, G. (2011). Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA. Hydrology and Earth System Sciences, 15(9), 2895-2911.

Wrede, S., Fenicia, F., Martinez-Carreras, N., Juilleret, J., Hissler, C., Krein, A., ... & Pfister, L. (2015). Towards more systematic perceptual model development: a case study using 3 Luxembourgish catchments. Hydrological Processes, 29(12), 2731-2750.

Hatchett, B.J. (2021). Seasonal and Ephemeral Snowpacks of the Conterminous United States. Hydrology, 8(1), 32.

Petersky, R., & Harpold, A. (2018). Now you see it, now you don't: a case study of ephemeral snowpacks and soil moisture response in the Great Basin, USA. Hydrology and Earth System Sciences, 22, 4891–4906.

Adelsperger, S., et al. (in review). A novel severity-based approach for assessing streamflow drought characteristics and drivers.

Laaha, G., et al. (2017). Unbiased plotting positions for low-flow frequency analysis.

Peters, N.E., & Aulenbach, B.T. (2011). Water storage at the Panola Mountain Research Watershed, Georgia, USA. Hydrological Processes, 25(25), 3878–3889.

## Part 4: Automated Data Quality Flags

All checks run on per-gage summary values (the _mean statistic, except where noted). Thresholds are config-driven (qa_qc section of signatures_config.json) and identical across the Julia, Python, and R implementations. A value of true marks a potential quality issue for review; flags never remove a gage, and a missing input value never triggers a flag.

flagged_for_qann_range: Qann_mean outside [0, 2000] mm

flagged_for_bfi_eckhardt_range: BFI_Eckhardt_mean outside [0, 1]. Defensive check — per-year BFI values are clamped to [0, 1] during calculation, so this flag cannot fire under current code.

flagged_for_bfi_lynehollick_range: BFI_LyneHollick_mean outside [0, 1]. Same defensive note as above.

flagged_for_flashiness_range: flashinessRB_mean outside [0, 2]

flagged_for_tqmean_range: TQmean_mean outside [0, 100]

flagged_for_d50_range: D50_day_mean outside [1, 366]

flagged_for_elasticity_range: elasticity_static outside [0.1, 5] (the per-gage scalar — the only flag not based on a _mean column)

flagged_for_runoff_ratio_range: annual_runoff_ratio_mean outside [0.01, 1.5]

flagged_for_seasonal_sum: the sum of the four seasonal totals (Qwin_mean + Qspr_mean + Qsum_mean + Qfal_mean) deviates from Qann_mean by more than 20%

flagged_for_percentile_order: non-decreasing order Q5 ≤ Q25 ≤ Q50 ≤ Q75 ≤ Q95 violated (on _mean values; ties are allowed and do not flag)

flagged_for_timing_order: non-decreasing order D5 ≤ D50 ≤ D95 violated (on _day_mean values; ties are allowed and do not flag)

flagged_for_high_na: more than 30% of a gage's numeric output columns are NA — often indicates missing climate/SWE coverage rather than bad streamflow data

var k,aa=typeof Object.create=="function"?Object.create:function(a){function b(){}b.prototype=a;return new b},ba=typeof Object.defineProperties=="function"?Object.defineProperty:function(a,b,c){if(a==Array.prototype||a==Object.prototype)return a;a[b]=c.value;return a};

function ca(a){a=["object"==typeof globalThis&&globalThis,a,"object"==typeof window&&window,"object"==typeof self&&self,"object"==typeof global&&global];for(var b=0;b>>0)+"_",e=0;return b});l("Symbol.dispose",function(a){return a?a:Symbol("f")});l("Object.is",function(a){return a?a:function(b,c){return b===c?b!==0||1/b===1/c:b!==b&&c!==c}});

l("Array.prototype.includes",function(a){return a?a:function(b,c){var d=this;d instanceof String&&(d=String(d));var e=d.length;c=c||0;for(c1114111||e!==Math.floor(e))throw new RangeError("g`"+e);e>>10&1023|55296),c+=String.fromCharCode(e&1023|56320))}return c}});

l("Array.prototype.find",function(a){return a?a:function(b,c){a:{var d=this;d instanceof String&&(d=String(d));for(var e=d.length,f=0;f>>0).toString(16))};r.prototype.j=["java.lang.Object",0];function w(){}m(w,r);function x(a,b){a.h=b;qa(b,a)}function y(a){ra(a.h)&&(Error.captureStackTrace?Error.captureStackTrace(z(a.h,ra,sa)):z(a.h,ra,sa).stack=Error().stack)}w.prototype.toString=function(){var a=u(v(this.constructor)),b=this.i;return b==null?a:t(a)+": "+t(b)};function ta(a){if(a!=null){var b=a.B;if(b)return b}a instanceof TypeError?b=ua():(b=new va,y(b),x(b,Error(b)));b.i=a==null?"null":a.toString();x(b,a);return b}w.prototype.j=["java.lang.Throwable",0];function wa(){}m(wa,w);wa.prototype.j=["java.lang.Exception",0];function A(){}m(A,wa);A.prototype.j=["java.lang.RuntimeException",0];function xa(){}m(xa,A);xa.prototype.j=["java.lang.IndexOutOfBoundsException",0];function oa(a,b){return Object.is(a,b)||a==null&&b==null};function ya(){}m(ya,A);ya.prototype.j=["java.lang.ClassCastException",0];function va(){}m(va,A);va.prototype.j=["java.lang.JsException",0];function za(){}m(za,va);function ua(){var a=new za;y(a);x(a,new TypeError(a));return a}za.prototype.j=["java.lang.NullPointerException",0];function Aa(){}m(Aa,xa);function Ba(a){var b=new Aa;b.i=a;y(b);x(b,Error(b));return b}Aa.prototype.j=["java.lang.StringIndexOutOfBoundsException",0];function Ca(){}m(Ca,r);Ca.prototype.j=["java.lang.Number",0];function Da(){}m(Da,Ca);Da.prototype.j=["java.lang.Double",0];function Ea(){}m(Ea,r);Ea.prototype.j=["java.lang.Boolean",0];function z(a,b,c){if(a!=null&&!b(a))throw a=t(u(Fa(a)))+" cannot be cast to "+t(u(v(c))),b=new ya,b.i=a,y(b),x(b,Error(b)),b.h;return a};function Ga(a,b,c){if(Object.prototype.hasOwnProperty.call(a.prototype,b))return a.prototype[b];c=c();return a.prototype[b]=c};var pa=0;function Fa(a){switch(Ha(typeof a)){case "number":return v(Da);case "boolean":return v(Ea);case "string":return v(B);case "function":return v(Ia)}if(a instanceof r)a=v(a.constructor);else if(Array.isArray(a))a=(a=a.H)?v(a.K,a.J):v(r,1);else if(a!=null)a=v(Ja);else throw new TypeError("h");return a};function Ia(){}Ia.prototype.j=["",1];function Ja(){}m(Ja,r);Ja.prototype.j=["",0];function Ka(a,b){C(a)?(La(b,a.length),a=a.charCodeAt(b)):a=a.h(b);return a};function Ha(a){if(a==null)throw ua().h;return a}function La(a,b){if(a=b)throw Ba("i`"+a+"`"+b).h;};function sa(){}function ra(a){return a instanceof Error}sa.prototype.j=["Error",0];function qa(a,b){if(a instanceof Object)try{a.B=b,Object.defineProperties(a,{cause:{get:function(){return b.o&&b.o.h}}})}catch(c){}};function B(){}m(B,r);function t(a){return a==null?"null":a.toString()}function Ma(a,b){var c=b,d=a.length,e;b=Ka(a,(e=c,c=c+1|0,e));var f;if(e=b>=55296&&b=56320&&a0&&d=a.length)break a;f=Ma(a,d+e|0);if(f==61||f==38||f==35)break a}d=d+(e+1)|0}d=-1}if(dc)e=c;d=d+(b.length+1)|0;b=decodeURIComponent;c=Math.min(a.length,d);d=Math.min(a.length,Math.max(d,e));e=a.length;if(ce||d/g,hb=/"/g,ib=/'/g,jb=/\x00/g,kb=/[\x00&<>"']/;/*

Copyright Google LLC

SPDX-License-Identifier: Apache-2.0

*/

var lb=globalThis.trustedTypes,mb;function nb(){var a=null;if(!lb)return a;try{var b=function(c){return c};a=lb.createPolicy("goog#html",{createHTML:b,createScript:b,createScriptURL:b})}catch(c){}return a};function D(a){this.h=a}D.prototype.toString=function(){return this.h};var ob=new D("about:invalid#zClosurez");function pb(a){if(a instanceof D)return a.h;throw Error("l");};function qb(a){this.C=a}function E(a){return new qb(function(b){return b.substr(0,a.length+1).toLowerCase()===a+":"})}var rb=[E("data"),E("http"),E("https"),E("mailto"),E("ftp"),new qb(function(a){return/^[^:]*([/?#]|$)/.test(a)})],sb=/^\s*(?!javascript:)(?:[\w+.-]+:|[^:/?#]*(?:[/?#]|$))/i;function tb(a){a instanceof D?a=pb(a):a=sb.test(a)?a:void 0;return a};function ub(a){this.h=a}ub.prototype.toString=function(){return this.h+""};function vb(a,b,c,d){b=tb(b);return b!==void 0?a.open(b,c,d):null};var wb=Array.prototype.indexOf?function(a,b){return Array.prototype.indexOf.call(a,b,void 0)}:function(a,b){if(typeof a==="string")return typeof b!=="string"||b.length!=1?-1:a.indexOf(b,0);for(var c=0;c=0;e--){var f=a[e]|0;d&&f==b||(c[e]=f,d=!1)}this.h=c}var ab={};function xb(a){return-128=c;d++)b[d]=a/c|0,c*=4294967296;return new F(b,0)}var H=xb(0),J=xb(1),yb=xb(16777216);

function K(a){if(M(a))return-K(I(a));for(var b=0,c=1,d=0;d=0?e:4294967296+e)*c;c*=4294967296}return b}k=F.prototype;k.toString=function(a){a=a||10;if(a0?c.h[0]:c.i)>>>0).toString(a);c=e;if(O(c))return f+d;for(;f.length>>16)+(N(this,e)>>>16)+(N(a,e)>>>16);d=g>>>16;f&=65535;g&=65535;c[e]=g>>16,g=N(this,d)&65535,h=N(a,e)>>>16,p=N(a,e)&65535;c[2*d+2*e]+=g*p;Ab(c,2*d+2*e);c[2*d+2*e+1]+=f*p;Ab(c,2*d+2*e+1);c[2*d+2*e+1]+=g*h;Ab(c,2*d+2*e+1);

c[2*d+2*e+2]+=f*h;Ab(c,2*d+2*e+2)}for(a=0;a>>16,a[b]&=65535,b++}function Q(a,b){this.h=a;this.i=b}

function zb(a,b){if(O(b))throw Error("o");if(O(a))return new Q(H,H);if(M(a))return b=zb(I(a),b),new Q(I(b.h),I(b.i));if(M(b))return b=zb(a,I(b)),new Q(I(b.h),b.i);if(a.h.length>30){if(M(a)||M(b))throw Error("n");for(var c=J,d=b;d.compare(a)=0;){c=Math.max(1,Math.floor(K(a)/K(b)));d=

Math.ceil(Math.log(c)/Math.LN2);d=d0;)c-=d,f=G(c),g=f.multiply(b);O(f)&&(f=J);e=e.add(f);a=P(a,g)}return new Q(e,a)}k.and=function(a){for(var b=Math.max(this.h.length,a.h.length),c=[],d=0;d>5;a%=32;for(var c=this.h.length+b+(a>0?1:0),d=[],e=0;e0?N(this,e-b)>>32-a:N(this,e-b);return new F(d,this.i)};function R(a,b){var c=b>>5;b%=32;for(var d=a.h.length-c,e=[],f=0;f0?N(a,f+c)>>>b|N(a,f+c+1)>24&1),Gb=!!(q[0]>>19&1),Hb=!!(q[0]>>26&1),Ib=!!(q[0]&128);var Jb=Eb?Hb:ka(610401301,!1),Kb=Eb?Gb||!Ib:ka(1331761403,!0);var Lb=null,S,Mb=n.navigator;S=Mb?Mb.userAgentData||null:null;function Nb(){var a;if(Lb==null)a:{if(a=n.navigator)if(a=a.userAgent)break a;a=""}else a=Lb;return a}function T(a){return Nb().indexOf(a)!=-1};function Ob(){if(!Jb)return!1;var a=S;return!!a&&a.brands.length>0}function Pb(a){var b={};a.forEach(function(c){b[c[0]]=c[1]});return function(c){return b[c.find(function(d){return d in b})]||""}}function Qb(){for(var a=Nb(),b=RegExp("([A-Z][\\w ]+)/([^\\s]+)\\s*(?:\\((.*?)\\))?","g"),c=[],d;d=b.exec(a);)c.push([d[1],d[2],d[3]||void 0]);Pb(c);return T("Firefox")||T("FxiOS")?(a=c[2])&&a[1]||"":""}function Db(){Ob()?S.brands.find(function(a){return a.brand==="Firefox"}):Qb()};var Cb=Nb().toLowerCase().indexOf("webkit")!=-1&&!T("Edge"),Rb;if(Jb){var Sb=S;Rb=!!Sb&&!!Sb.platform}else Rb=!1;var Tb=Rb?S.platform==="macOS":T("Macintosh");function Ub(a){a&&typeof a.dispose=="function"&&a.dispose()};function U(){this.i=this.i;this.h=this.h}U.prototype.i=!1;U.prototype.dispose=function(){this.i||(this.i=!0,this.o())};U.prototype[Symbol.dispose]=function(){this.dispose()};U.prototype.o=function(){if(this.h)for(;this.h.length;)this.h.shift()()};function Vb(a,b){this.type=a;this.target=b;this.l=!1}Vb.prototype.h=function(){this.l=!0};var Wb=function(){if(!n.addEventListener||!Object.defineProperty)return!1;var a=!1,b=Object.defineProperty({},"passive",{get:function(){a=!0}});try{var c=function(){};n.addEventListener("test",c,b);n.removeEventListener("test",c,b)}catch(d){}return a}();function V(a,b){Vb.call(this,a?a.type:"");this.relatedTarget=this.target=null;this.button=this.screenY=this.screenX=this.clientY=this.clientX=0;this.key="";this.metaKey=this.shiftKey=this.altKey=this.ctrlKey=!1;this.state=null;this.pointerId=0;this.pointerType="";this.i=null;a&&this.init(a,b)}na(V,Vb);

V.prototype.init=function(a){var b=this.type=a.type,c=a.changedTouches&&a.changedTouches.length?a.changedTouches[0]:null;this.target=a.target||a.srcElement;var d=a.relatedTarget;d||(b=="mouseover"?d=a.fromElement:b=="mouseout"&&(d=a.toElement));this.relatedTarget=d;c?(this.clientX=c.clientX!==void 0?c.clientX:c.pageX,this.clientY=c.clientY!==void 0?c.clientY:c.pageY,this.screenX=c.screenX||0,this.screenY=c.screenY||0):(this.clientX=a.clientX!==void 0?a.clientX:a.pageX,this.clientY=a.clientY!==void 0?

a.clientY:a.pageY,this.screenX=a.screenX||0,this.screenY=a.screenY||0);this.button=a.button;this.key=a.key||"";this.ctrlKey=a.ctrlKey;this.altKey=a.altKey;this.shiftKey=a.shiftKey;this.metaKey=a.metaKey;this.pointerId=a.pointerId||0;this.pointerType=a.pointerType;this.state=a.state;this.i=a;a.defaultPrevented&&V.u.h.call(this)};V.prototype.h=function(){V.u.h.call(this);var a=this.i;a.preventDefault?a.preventDefault():a.returnValue=!1};var Xb="closure_listenable_"+(Math.random()*1E6|0);var Yb=0;function Zb(a,b,c,d,e){this.listener=a;this.proxy=null;this.src=b;this.type=c;this.capture=!!d;this.handler=e;this.key=++Yb;this.h=this.v=!1}function $b(a){a.h=!0;a.listener=null;a.proxy=null;a.src=null;a.handler=null};function ac(a){this.src=a;this.h={};this.i=0}ac.prototype.add=function(a,b,c,d,e){var f=a.toString();a=this.h[f];a||(a=this.h[f]=[],this.i++);var g;a:{for(g=0;g-1?(b=a[g],c||(b.v=!1)):(b=new Zb(b,this.src,f,!!d,e),b.v=c,a.push(b));return b};var bc="closure_lm_"+(Math.random()*1E6|0),cc={},dc=0;function ec(a,b,c,d,e){if(d&&d.once)return fc(a,b,c,d,e);if(Array.isArray(b)){for(var f=0;f=0)&&Array.prototype.splice.call(e,f,1);g&&($b(a),c.h[d].length==0&&(delete c.h[d],c.i--))}c.i==0&&(c.src=null,b[bc]=null)}else $b(a)}}}

function lc(a){return a in cc?cc[a]:cc[a]="on"+a}function mc(a,b){if(a.h)a=!0;else{b=new V(b,this);var c=a.listener,d=a.handler||a.src;a.v&&nc(a);a=c.call(d,b)}return a}function jc(a){a=a[bc];return a instanceof ac?a:null}var oc="__closure_events_fn_"+(Math.random()*1E9>>>0);function hc(a){if(typeof a==="function")return a;a[oc]||(a[oc]=function(b){return a.handleEvent(b)});return a[oc]};function X(a){U.call(this);this.m=a;this.l={}}na(X,U);var pc=[];function qc(a){db(a.l,function(b,c){this.l.hasOwnProperty(c)&&nc(b)},a);a.l={}}X.prototype.o=function(){X.u.o.call(this);qc(this)};X.prototype.handleEvent=function(){throw Error("r");};function rc(a,b,c,d){b=b===void 0?!1:b;c=c===void 0?!1:c;d=d===void 0?"editors":d;U.call(this);this.l=a||document.body;this.D=!!b;this.F=!!c;this.G=d;this.m=new X(this);a=ma(Ub,this.m);this.i?a():(this.h||(this.h=[]),this.h.push(a));a=this.m;b=this.l;c=this.A;d="click";Array.isArray(d)||(d&&(pc[0]=d.toString()),d=pc);for(var e=0;e0&&!e&&!Fb)L.length===2&&(L=["noreferrer"]),vb(d,b,f,L[0]);else if(e=p.join(","),(T("iPhone")&&!T("iPod")&&!T("iPad")||T("iPad")||T("iPod"))&&d.navigator&&d.navigator.standalone&&f&&f!="_self")e="A",h=document,e=String(e),h.contentType==="application/xhtml+xml"&&(e=e.toLowerCase()),h=e=h.createElement(e),b=tb(b),b!==void 0&&(h.href=b),e.target=f,g&&(e.rel="noreferrer"),((c=c.attributionsrc)||c==="")&&e.setAttribute("attributionsrc",c),c=document.createEvent("MouseEvent"),

c.initMouseEvent("click",!0,!0,d,1),e.dispatchEvent(c);else if(g){if(d=vb(d,"",f,e),c=pb(b),d&&(d.opener=null,c===""&&(c="javascript:''"),kb.test(c)&&(c.indexOf("&")!=-1&&(c=c.replace(eb,"&")),c.indexOf("")!=-1&&(c=c.replace(gb,">")),c.indexOf('"')!=-1&&(c=c.replace(hb,""")),c.indexOf("'")!=-1&&(c=c.replace(ib,"'")),c.indexOf("\x00")!=-1&&(c=c.replace(jb,"�"))),c='',mb===void 0&&(mb=nb()),c=(b=mb)?b.createHTML(c):c,c=new ub(c),(d=d.document)&&d.write)){b=d;e=b.write;if(c instanceof ub)c=c.h;else throw Error("l");e.call(b,c);d.close()}}else(d=vb(d,b,f,e))&&c.noopener&&(d.opener=null),d&&c.noreferrer&&(d.opener=null);a.h();break}}b=b.parentNode}};function sc(a,b,c,d){new rc(a,b===void 0?!1:b,c===void 0?!1:c,d)}

for(var tc=["DOCS_installLinkReferrerSanitizer"],Y=n,Z;tc.length&&(Z=tc.shift());)tc.length||sc===void 0?Y[Z]&&Y[Z]!==Object.prototype[Z]?Y=Y[Z]:Y=Y[Z]={}:Y[Z]=sc;

// Google Inc.

//# sourceMappingURL=linkreferrer_binary_linkreferrer_binary_chunk.sourcemap

DOCS_installLinkReferrerSanitizer( undefined, true, true, 'docs');_docs_flag_initialData={"ilcm":{"eui":"ADFN-cu604P0d6Dt3moDCv2PtlhO3FmXC79a6CqzCmTb-hJOtIsaTPKBt2-v9q2AnFntlqAAk5jr","je":1,"sstu":1788547843613942,"si":"CPfJh_fL1ZYDFSqqIAEdOq4fnw","gsc":0,"ei":[5703839,5703890,5704571,5704621,5704883,5704905,5705891,5706270,5707047,5707075,5707204,5708365,5708574,5709476,5710587,5710692,5711550,5712635,5712909,5713195,5713554,5714310,5715711,5715791,5716179,5716196,5717567,5723473,5724199,5724215,5724880,5724896,5725855,5727701,5729076,5729092,5734616,5734632,5735236,5735252,5737784,5737800,5742636,5742652,5743771,5743787,5746708,5746724,5753665,5753681,5754923,5754939,5758658,5758674,5760151,5760167,5760454,5760470,5763814,5763830,5780496,5780512,5781974,5781990,5782329,5782345,5787796,5792429,5799519,5799527,5799571,5799579,5799644,5799652,5799822,5799830,5799875,5799883,5799906,5799914,5799927,5799935,13702623,48966234,48966242,49398661,49398669,49472103,49472111,49491737,49491745,49498973,49498981,49623441,49623449,49629214,49629222,49644055,49644063,49661340,49661348,49769417,49769425,49822841,49822849,49833562,49833570,49842875,49842883,49894860,49894868,49904439,49904447,49924626,49924634,49926153,49926161,49943099,49943107,49971993,49972001,49974154,49974162,49979538,49979546,50223794,50223802,50266102,50266110,50273468,50273476,50293588,50293596,50297136,50297144,50389210,50389218,50438805,50439240,50439248,50513094,50538784,50538792,50550111,50550119,50561463,50561471,50562741,50562752,50578543,50578551,50587122,50587130,50596469,50602181,50602189,70971116,70971124,70973202,70973210,71061469,71079958,71079966,71085281,71085289,71120928,71120936,71145354,71145365,71200393,71200401,71289801,71289809,71291305,71291313,71316507,71316515,71330194,71330202,71331358,71331366,71376156,71376164,71387387,71387398,71387809,71387817,71406857,71406865,71451142,71451150,71465967,71505700,71505708,71520490,71520498,71528637,71528645,71530223,71530231,71544794,71544802,71546405,71558484,71558492,71608300,71608308,71628337,71628345,71628367,71628375,71628457,71628465,71638483,71638491,71657940,71657948,71659933,71659941,71675235,71679620,71679628,71681950,71689900,71689908,71710060,71710068,71710693,71710709,71825483,71825491,71832993,71833032,71854950,71854958,71897947,71897955,71899323,71899334,71960420,71960428,94326699,94327571,94327579,94353348,94353356,94386894,94394912,94394920,94409697,94434317,94434325,94460799,94460815,94461479,94629737,94629745,94640319,94640327,94661722,94661730,94692338,94692346,94733477,94733485,94813624,94813635,94871310,94871318,94904209,94919072,94919080,94940334,95087166,95087174,95093652,95093668,95099833,95099841,95104414,95104425,95111885,95111893,95112733,95112741,95113216,95113224,95131173,95131181,95260593,95299606,95299614,95309227,95309235,95314722,95314730,99251823,99251831,99257767,99257775,99291459,99338860,99338868,99342838,99342846,99343448,99343456,99368732,99368740,99393819,99400262,99400270,99402231,99402239,99440933,99440941,99457586,99457594,99457827,99457835,99778386,99778391,100640000,101448020,101448275,101448280,101478046,101478054,101489651,101489656,101489906,101489911,101492831,101492839,101540517,101540525,101561672,101561680,101562346,101562354,101575587,101575592,101631231,101631239,101632054,101632070,101672701,101687057,101687065,101718603,101739377,101755356,101755361,101755547,101755552,101781895,101781911,101788143,101788151,101793871,101823263,101823271,101847134,101858214,101858219,101860627,101860635,101867803,101867811,101868111,101868133,101887634,101887642,101889011,101896291,101896379,101917205,101917213,101919498,101919506,101920169,101920174,101922599,101922607,101922980,101922988,101933591,101933599,101962010,101962022,101962589,102013387,102013395,102030502,102030510,102059460,102059468,102070596,102070604,102070890,102070898,102083186,102083194,102146707,102146715,102161225,102161233,102161587,102161595,102178352,102178360,102198442,102198450,102201965,102205353,102205361,102205628,102205636,102208252,102208260,102236046,102236054,102244469,102244477,102280748,102280756,102287344,102287352,102292136,102292144,102322635,102322651,102343360,102343368,102388302,102388307,102399901,102399909,102402779,102402787,102448099,102448115,102449888,102449896,102469780,102469788,102517165,102517170,102548661,102548669,102554627,102554632,102554869,102554874,102587753,102587761,102600150,102600155,102618668,102618676,102628048,102641327,102641332,102649643,102649648,102649928,102649933,102672672,102672677,102673089,102673094,102673355,102685397,102685402,102685916,102685921,102691121,102691126,102691386,102691391,102718422,102721112,102721117,102727389,102727394,102752415,102761409,102761414,102762335,102762343,102783389,102783397,102784888,102784896,102787523,102807730,102807738,102811903,102811911,102858946,102858954,102863247,102863252,102864324,102864332,102867795,102867803,102883798,102883806,102902676,102902684,102903580,102903588,102909802,102909807,102925899,102925907,102926483,102926491,102932527,102934198,102944137,102944142,102944346,102944351,102944425,102944430,102944555,102944560,102944693,102944698,102944811,102944816,102956645,102956653,102972729,102972737,102986857,102986865,102988337,102988342,102988441,102988446,102988747,102988752,102988916,102988921,103011599,103011607,103015173,103018775,103018781,103020595,103020603,103039144,103039152,103058195,103058211,103163004,103163012,103175930,103175938,103209325,103209333,103285926,103285934,103286093,103286101,103286424,103286432,103288266,103288274,103289142,103289147,103289221,103289226,103299198,103299206,103326447,103326455,103339776,103339784,103343810,103343818,104530118,104574039,104574047,104575397,104575402,104575543,104575548,104575639,104575644,104575920,104575925,104615555,104615563,104627936,104627942,104654326,104654334,104699524,104699532,104717880,104717896,104746802,104746807,104761731,104761739,104790721,104794503,104798707,104798723,104907771,104907779,104959355,104963462,104963470,104976557,104976565,104983089,104983097,105032492,105082057,105082065,105083768,105083776,105084989,105084997,105087628,105087636,105111908,105111916,105125337,105125345,105158864,105180813,105180821,105197504,105217109,105217117,105250945,105250953,105255698,105255706,105283747,105283755,105292976,105292984,105345084,105345089,105346125,105346133,105360298,105360306,105371963,105374004,105374012,105393831,105393839,105439291,105439299,115492195,115492211,115511311,115511316,115601169,115601177,115626338,115626346,115626399,115626407,115665132,115665137,115665205,115665210,115666017,115666025,115669418,115669426,115680658,115738038,115738046,115748334,115748340,115755264,115755272,115780021,115780029,115780478,115780486,115797565,115797570,115822485,115822493,115849188,115849193,115854828,115854836,115902264,115902269,115909781,115909789,115916250,115916255,115917358,115917366,115928899,115928907,115986507,115986513,116036386,116036402,116062770,116062778,116068758,116068766,116076925,116076931,116091118,116150705,116150713,116176889,116176897,116214049,116214057,116221354,116221362,116229714,116229722,116250785,116250793,116276281,116276289,116307748,116339544,116339549,116361234,116372457,116372465,116407665,116416014,116416022,116418247,116418255,116425644,116426540,116448407,116449517,116459613,116470832,116479751,116480001,116495787,116496398,116502867,116503963,116508642,116508647,116514183,116514189,116524154,116538556,116538564,116542683,116550348,116567535,116567540,116593745,116609391,116609399,116613094,116613710,116685399,116685789,116697553,116697559,116699409,116699417,116704082,116704120,116714309,116714312,116731449,116731455,116731463,116731469,116731477,116731483,116731491,116731497,116731505,116731511,116731519,116731525,116757124,116757439,116757444,116771270,116771276,116778957,116778962,116780649,116780655,116780663,116780669,116780677,116780683,116780691,116780697,116780705,116780711,116780719,116780725,116783786,116783791,116840720,116840725,116845807,116845812,116847573,116847578,116874389,116874397,116874913,116877099,116882876,116882883,116913250,116918027,116919494,116923675,116924992,116980527,116980535,116982432,116982440,116982914,116982922,116988224,116988232,117010853,117010861,117026763,117026878,117026880,117037168,117037176,117037622,117037626,117091085,117124349,117124356,117136346,117143742,117143747,117145792,117145800,117148946,117148952,117166647,117224099,117224107,117236066,117249330,117249338,117267913,117269057,117269065,117284159,117284164,117286731,117286739,117290310,117290318,117297869,117297874,117309234,117309242,117310170,117310178,117357715,117357721,117367088,117367094,117382832,117387079,117387085,117411961,117411966,117457918,117457926,117474409,117474414,117503097,117503105,117504343,117504351,117528468,117528471,117530134,117547521,117547529,117556935,117557280,117586569,117586574,117605990,117605995,117617889,117623323,117623329,117632590,117642488,117642592,117642597,117642980,117642988,117643748,117643756,117683635,117706017,117711597,117723591,117723597,117728473,117728481,117734392,117736031,117740479,117740485,117752075,117752081,117755703,117755716,117755820,117756394,117756755,117756763,117775368,117782123,117852663,117852668,117855973,117855979,117879413,117879418,117895920,117895928,117930842,117930847,117952935,117970097,117970102,117983432,117983437,117995685,117995689,117997204,117999163,118010935,118026977,118029473,118029479,118043880,118043885,118044535,118048893,118048900,118057575,118057581,118086120,118086126,118094462,118094467,118098329,118098337,118109056,118110654,118110662,118119539,118119545,118127424,118127429,118132041,118132049,118144499,118144504,118144972,118144980,118155108,118155115,118155211,118155217,118155412,118171214,118172011,118179259,118181389,118181395,118182951,118183989,118199432,118199438,118210158,118210163,118213580,118213585,118219914,118222565,118227602,118227608,118237155,118237163,118240038,118240044,118240393,118246247,118254618,118255415,118255420,118256780,118256787,118269183,118269270,118269275,118275356,118281124,118281132,118292262,118292268,118292276,118292282,118311122,118311127,118325328,118342369,118346423,118347706,118349323,118349326,118349356,118364657,118367501,118367507,118367791,118368037,118368045,118368696,118368704,118373667,118419982,118419988,118427971,118427977,118449075,118462069,118482248,118482256,118489188,118494715,118497529,118497535,118509494,118511786,118515459,118515465,118517998,118518004,118519229,118521972,118521977,118528660,118528668,118531377,118532140,118532142,118541625,118541633,118544270,118544277,118555977,118555985,118566888,118566900,118567234,118567242,118575738,118575744,118576741,118579779,118579784,118596018,118596023,118599481,118599489,118600863,118623871,118623879,118627810,118627821,118654314,118656213,118656218,118656880,118656885,118669236,118670990,118670996,118683648,118683653,118717389,118728532,118728537,118735214,118735220,118736895,118753224,118753232,118780368,118781276,118781281,118786501,118786509,118792272,118795502,118795508,118796216,118803374,118803380,118807052,118809999,118810531,118824720,118824725,118825246,118830481,118830620,118830622,118830624,118830626,118830628,118830630,118830632,118830634,118830636,118830640,118835027,118855606,118870266,118871022,118884772,118896670,118896675,118910283,118910289,118925905,118925913,118931884,118938363,118938369,118952709,118952717,118957772,118957777,118960868,118960876,118965149,118965154,118970617,119001820,119001825,119012212,119012218,119012461,119012477,119022212,119022220,119036798,119036806,119038738,119038746,119058017,119058022,119061297,119062070,119062078,119070701,119071842,119071848,119076592,119076597,119077195,119077203,119078243,119087275,119092108,119096521,119096529,119106920,119126674,119127375,119140584,119140592,119142923,119159957,119164990,119164992,119164994,119168749,119170344,119170349,119178779,119187339,119189027,119189032,119190842,119190847,119194829,119194834,119202558,119202566,119206371,119206377,119207504,119207509,119209341,119209377,119219121,119219126,119220909,119220911,119220913,119220915,119220917,119220919,119221438,119235564,119235569,119237704,119237709,119239792,119249911,119255437,119261740,119261748,119262908,119265733,119275421,119287095,119287101,119290941,119290946,119293882,119293887,119311194,119311199,119317348,119317356,119335774,119335779,119346731,119350203,119350208,119370737,119372604,119372610,119374110,119374116,119382576,119400586,119400592,119412173,119412189,119412195,119419025,119419030,119419128,119419133,119422183,119422199,119431313,119446111,119446117,119453372,119454198,119454203,119454286,119454291,119455720,119455728,119455933,119455939,119464287,119464292,119465389,119465397,119471147,119474814,119474822,119486429,119497026,119497031,119502682,119502687,119520108,119522295,119522311,119525128,119525133,119528282,119528290,119530180,119530188,119538769,119538776,119546473,119546479,119546481,119546486,119552616,119566068,119566076,119566747,119566753,119576796,119576803,119577660,119577668,119588321,119588326,119590765,119590807,119590821,119590826,119590935,119592860,119592868,119595187,119595195,119599143,119605170,119605178,119608194,119608200,119612419,119612423,119614697,119614703,119615541,119615546,119621403,119621408,119621918,119621936,119636046,119643427,119643433,119666107,119666115,119680681,119680688,119681665,119681673,119684288,119684294,119684496,119685915,119685921,119702649,119702657,119703104,119703112,119707152,119707157,119707161,119727967,119727972,119741495,119741503,119745075,119761002,119761007,119765383,119765391,119768442,119770485,119770493,119773246,119773713,119776352,119794118,119794312,119794317,119804024,119804032,119804531,119805088,119805096,119807183,119809651,119809656,119820261,119826320,119826325,119836044,119840989,119842248,119842254,119854316,119880872,119880877,119885476,119914240,119914248,119914253,119914638,119914643,119919341,119919349,119930504,119930509,119932762,119932767,119941748,119941750,119944463,119976748,119976753,119989280,119999033,120006922,120006930,120008223,120008231,120024054,120025498,120025503,120049011,120049019,120050620,120050628,120051513,120051518,120073901,120073907,120074338,120074344,120082244,120082249,120098782,120098787,120107754,120118521,120118527,120120566,120120575,120120580,120129974,120129980,120129988,120129994,120130002,120130008,120140748,120140753,120140755,120140959,120140964,120142704,120142711,120142990,120142997,120143451,120143457,120161430,120161435,120163559,120163564,120167720,120167725,120169620,120169628,120175687,120175695,120179202,120179208,120182621,120182626,120182928,120182944,120184008,120187257,120187262,120190042,120190047,120191501,120191506,120195953,120221882,120221887,120222179,120226251,120226258,120226267,120226273,120232064,120239562,120239570,120250425,120250430,120256334,120256340,120257682,120257690,120268393,120272513,120275434,120275439,120286666,120286671,120292177,120300727,120308624,120343165,120381605,120381613,120388576,120389577,120389585,120395666,120408009,120408014,120419550,120419556,120426040,120426046,120426054,120426060,120426068,120426074,120426082,120426088,120426090,120426095,120426110,120426116,120426124,120426130,120426138,120426144,120436312,120436318,120437412,120439867,120439872,120446673,120456509,120456514,120459269,120465957,120465965,120483765,120483771,120493262,120493270,120495316,120495323,120540098,120554201,120554207,120563795,120573339,120576479,120579829,120579837,120596839,120596845,120603282,120603287,120614405,120614412,120652083,120652090,120661035,120661040,120666070,120666077,120666798,120672603,120674270,120674277,120674439,120674441,120674443,120674444,120674447,120674449,120678228,120705172,120711524,120711531,120712193,120712199,120713457,120715203,120715208,120716393,120730641,120730646,120730864,120730872,120731452,120731459,120732417,120732422,120735838,120735843,120772217,120777459,120777464,120794198,120794203,120797333,120797339,120822134,120827764,120827769,120844533,120844538,120844732,120844737,120845137,120845142,120864326,120887325,120887330,120894578,120931070,120931079,120933033,120951758,120957398,120957404],"crc":0,"cvi":[]},"docs-cclt":56}; _docs_flag_cek=''; if (window['DOCS_timing']) {DOCS_timing['ifdld']=new Date().getTime();}DOCS_timing['sjl'] = performance.now();this._pubi=this._pubi||{};(function(_){var window=this;

_._F_toggles_initialize=function(a){(typeof globalThis!=="undefined"?globalThis:typeof self!=="undefined"?self:this)._F_toggles__pubi=a||[]};(0,_._F_toggles_initialize)([]);

/*

Copyright The Closure Library Authors.

SPDX-License-Identifier: Apache-2.0

*/

var aa=function(){return function(a){return a}};var r=function(){return function(){}};var ba=function(a){return function(b){this[a]=b}};var ca=function(a){return function(){return a}};

var t,ea=function(a){a=["object"==typeof globalThis&&globalThis,a,"object"==typeof window&&window,"object"==typeof self&&self,"object"==typeof global&&global];for(var b=0;b2){var d=Array.prototype.slice.call(arguments,2);return function(){var e=Array.prototype.slice.call(arguments);Array.prototype.unshift.apply(e,d);return a.apply(b,e)}}return function(){return a.apply(b,arguments)}},Ma=function(a,b,c){Ma=Function.prototype.bind&&Function.prototype.bind.toString().indexOf("native code")!=

-1?Ka:La;return Ma.apply(null,arguments)},Na=function(a,b){var c=Array.prototype.slice.call(arguments,1);return function(){var d=c.slice();d.push.apply(d,arguments);return a.apply(this,d)}},Oa=function(a,b){z[a]=b},Pa=aa(),Qa=function(a,b){function c(){}c.prototype=b.prototype;a.T=b.prototype;a.prototype=new c;a.prototype.constructor=a;a.jb=function(d,e,f){for(var g=Array(arguments.length-2),h=2;h0:!1},Ya=function(){return Xa()?!1:C("Opera")},Za=function(){return Xa()?

Wa("Microsoft Edge"):C("Edg/")},$a=function(){return C("Firefox")||C("FxiOS")},ab=function(){return Xa()?Wa("Chromium"):(C("Chrome")||C("CriOS"))&&!(Xa()?0:C("Edge"))||C("Silk")},bb=function(){return Ua?!!Va&&!!Va.platform:!1},cb=function(){return C("iPhone")&&!C("iPod")&&!C("iPad")},db=function(){return cb()||C("iPad")||C("iPod")},eb=function(){return bb()?Va.platform==="macOS":C("Macintosh")},fb=function(){return bb()?Va.platform==="Chrome OS":C("CrOS")},gb=function(a){gb[" "](a);return a},hb=function(a,

b){a.__closure__error__context__984382||(a.__closure__error__context__984382={});a.__closure__error__context__984382.severity=b},ib=function(a){a=Error(a);hb(a,"warning");return a},kb=function(a,b){if(a!=null){var c;var d=(c=jb)!=null?c:jb={};c=d[a]||0;c>=b||(d[a]=c+1,a=Error(),hb(a,"incident"),Sa(a))}},lb=function(){return typeof BigInt==="function"},mb=function(a,b,c){return typeof Symbol==="function"&&typeof Symbol()==="symbol"?(c===void 0?0:c)&&Symbol.for&&a?Symbol.for(a):a!=null?Symbol(a):Symbol():

b},rb=function(a,b){ob||D in a||pb(a,qb);a[D]|=b},F=function(a,b){ob||D in a||pb(a,qb);a[D]=b},tb=function(a,b){return b===void 0?a.h!==sb&&!!(2&(a.m[D]|0)):!!(2&b)&&a.h!==sb},ub=function(a){a.mb=!0;return a},zb=function(a){var b=a;if(vb(b)){if(!/^\s*(?:-?[1-9]\d*|0)?\s*$/.test(b))throw Error(String(b));}else if(wb(b)&&!Number.isSafeInteger(b))throw Error(String(b));return xb?BigInt(a):a=yb(a)?a?"1":"0":vb(a)?a.trim()||"0":String(a)},Ab=function(a,b){if(a.length>b.length)return!1;if(a.lengthe)return!1;if(d>>0;H=b;Bb=(a-b)/4294967296>>>0},Eb=function(a){if(a>>0;Bb=b>>>0}else Cb(a)},Gb=function(a,b){b>>>=0;a>>>=0;if(b>>24|b>16&65535,a=(a&16777215)+c*6777216+b*6710656,c+=b*8147497,b*=2,a>=1E7&&(c+=a/1E7>>>0,a%=1E7),

c>=1E7&&(b+=c/1E7>>>0,c%=1E7),c=b+Fb(c)+Fb(a));return c},Fb=function(a){a=String(a);return"0000000".slice(a.length)+a},Hb=function(){var a=H,b=Bb;b&2147483648?lb()?a=""+(BigInt(b|0)>>0)):(b=x(Db(a,b)),a=b.next().value,b=b.next().value,a="-"+Gb(a,b)):a=Gb(a,b);return a},Db=function(a,b){b=~b;a?a=~a+1:b+=1;return[a,b]},Ib=function(a){return Array.prototype.slice.call(a)},Jb=function(a){return a.displayName||a.name||"unknown type name"},Mb=function(a){switch(typeof a){case "bigint":return!0;

case "number":return Kb(a);case "string":return Lb.test(a);default:return!1}},Nb=function(a){if(typeof a!=="number")throw ib("int32");if(!Kb(a))throw ib("int32");return a|0},Ob=function(a){if(a==null)return a;if(typeof a==="string"&&a)a=+a;else if(typeof a!=="number")return;return Kb(a)?a|0:void 0},Vb=function(a){var b=void 0;b!=null||(b=1024);if(!Mb(a))throw ib("int64");var c=typeof a;switch(b){case 512:switch(c){case "string":return Pb(a);case "bigint":return String(Qb(64,a));default:return Rb(a)}case 1024:switch(c){case "string":return Sb(a);

case "bigint":return zb(Qb(64,a));default:return Tb(a)}case 0:switch(c){case "string":return Pb(a);case "bigint":return zb(Qb(64,a));default:return Ub(a)}default:throw Error("Unknown format requested type for int64");}},Wb=function(a){var b=a.length;if(a[0]==="-"?b>>0,Bb=Number(a>>BigInt(32)&BigInt(4294967295));else{b=+(a[0]===

"-");Bb=H=0;for(var c=a.length,d=b,e=(c-b)%6+b;e=4294967296&&(Bb+=Math.trunc(H/4294967296),Bb>>>=0,H>>>=0);b&&(b=x(Db(H,Bb)),a=b.next().value,b=b.next().value,H=a,Bb=b)}return Hb()},Ub=function(a){Mb(a);a=Xb(a);if(!Yb(a)){Eb(a);var b=H,c=Bb;if(a=c&2147483648)b=~b+1>>>0,c=~c>>>0,b==0&&(c=c+1>>>0);var d=c*4294967296+(b>>>0);b=Number.isSafeInteger(d)?d:Gb(b,c);a=typeof b==="number"?a?-b:b:a?"-"+b:b}return a},Rb=function(a){Mb(a);a=Xb(a);Yb(a)?a=

String(a):(Eb(a),a=Hb());return a},Pb=function(a){Mb(a);var b=Xb(Number(a));if(Yb(b))return String(b);b=a.indexOf(".");b!==-1&&(a=a.substring(0,b));return Wb(a)},Sb=function(a){var b=Xb(Number(a));if(Yb(b))return zb(b);b=a.indexOf(".");b!==-1&&(a=a.substring(0,b));return lb()?zb(Qb(64,BigInt(a))):zb(Wb(a))},Tb=function(a){return Yb(a)?zb(Ub(a)):zb(Rb(a))},Zb=function(a){var b=typeof a;if(a==null)return a;if(b==="bigint")return zb(Qb(64,a));if(Mb(a))return b==="string"?Sb(a):Tb(a)},ec=function(a,b){if(!(a instanceof

b))throw Error("x`"+Jb(b)+"`"+(a&&Jb(a.constructor)));return a},ic=function(a,b,c){if(a!=null&&a[fc]===hc)return a;if(Array.isArray(a)){var d=a[D]|0;c=d|c&32|c&2;c!==d&&F(a,c);return new b(a)}},jc=aa(),lc=function(a,b){b=g){var n=e-m,u=void 0;((u=b)!=null?u:b={})[n]=q}else f[e]=q;if(p)for(var G in p)a=p[G],a!=null&&(a=c(a,d))!=null&&(h=+G,e=void 0,l&&!Number.isNaN(h)&&(e=h+m)=1024)throw Error("D");for(var k in h)f=+k,f1024)throw Error("E");e=e&-16760833|(k&1023)=f){var g=a[f];if(g!=null&&typeof g==="object"&&g.constructor===Object){c=g[b];var h=!0}else if(e===f)c=g;else return}else c=a[e];if(d&&c!=null){d=d(c);if(d==null)return d;if(!Object.is(d,c))return h?g[b]=d:a[e]=d,d}return c}},Hc=function(a,b,c){Dc(a);var d=a.m;Gc(d,d[D]|0,b,c);return a},Gc=function(a,b,c,d,e){var f=c+(e?0:-1),g=a.length-1;if(g>=1+(e?0:-1)&&f>=g){var h=a[g];if(h!=null&&typeof h==="object"&&h.constructor===

Object)return h[c]=d,b}if(f>14&1023||536870912;c>=g?d!=null&&(f={},a[g+(e?0:-1)]=(f[c]=d,f)):a[f]=d}return b},K=function(a,b,c){a=a.m;return Ic(a,a[D]|0,b,c)!==void 0},Lc=function(a,b,c,d,e,f,g,h){var k=b;f===1||(f!==4?0:2&b||!(16&b)&&32&d)?Jc(b)||(b|=!a.length||g&&!(4096&b)||32&d&&!(4096&b||16&b)?2:256,b!==k&&F(a,b),Object.freeze(a)):(f===2&&Jc(b)&&(a=Ib(a),k=0,b=Kc(b,d),d=Gc(c,d,e,a)),Jc(b)||(h||(b|=16),b!==k&&F(a,b)));2&b||!(4096&

b||16&b)||Ec(c,d);return a},Nc=function(a,b){a=Fc(a,b);return Array.isArray(a)?a:Mc},Oc=function(a,b){2&b&&(a|=2);return a|1},Jc=function(a){return!!(2&a)&&!!(4&a)||!!(256&a)},Ic=function(a,b,c,d,e){var f=!1;d=Fc(a,d,e,function(g){var h=ic(g,c,b);f=h!==g&&h!=null;return h});if(d!=null)return f&&!tb(d)&&Ec(a,b),d},L=function(a,b,c,d){var e=a.m,f=e[D]|0;b=Ic(e,f,b,c,d);if(b==null)return b;f=e[D]|0;if(!tb(a,f)){var g=Ac(b);g!==b&&(Cc(a)&&(e=a.m,f=e[D]|0),b=g,f=Gc(e,f,c,b,d),Ec(e,f))}return b},Pc=function(a,

b,c,d,e,f,g,h){var k=tb(a,c);f=k?1:f;g=!!g||f===3;k=h&&!k;(f===2||k)&&Cc(a)&&(b=a.m,c=b[D]|0);a=Nc(b,e);var l=a===Mc?7:a[D]|0,m=Oc(l,c);if(h=!(4&m)){var p=a,q=c,n=!!(2&m);n&&(q|=2);for(var u=!n,G=!0,E=0,da=0;E>>0)},Fd=function(a){return a.B==0&&a.A==0},Gd=function(a){var b=~a.B+1|0;return V(b,~a.A+!b|0)},Kd=function(a){return a>0?a>=0x7fffffffffffffff?Hd:new Dd(a,a/4294967296):a>>0).toString(16);b=W("0".repeat(Math.max(0,8-b.length|0)))+W(b);a=(a(2147483647)>>>0).toString(16);return W(a)+W(b)},Xd=function(){this.h=new Wd},Yd=function(a){this.m=I(a)},Wd=function(a){this.m=I(a)},Zd=function(a){this.m=I(a)},be=function(a){a==null||$d(a);return a==null?null:ae(a)},ae=function(a){$d(a);

pc(a);return pc(a)?Number(a):String(a)},ce=function(a){this.m=I(a)},de=function(a){this.m=I(a)},ee=function(a){this.m=I(a)},he=function(){if(!fe){var a=new ge(null);fe=function(){return a}}var b;return b=fe,b()},ie=r(),ge=function(a){this.h=new ie;this.i=null;if(a!=null)for(var b in a){var c=b,d=a[b];if(this.i)throw Nd("P").h;var e=this.h.get();e[c]=d!=null?d:null}},je=function(a,b){a=a.h.get();return b in a},X=function(a,b){a=a.get(b);return typeof a=="string"?a=="true"||a=="1":!!a},me=function(a,

b){ke(a,b);if(!je(a,b)||a.get(b)==null)return NaN;try{var c=W(a.get(b));le||(le=RegExp("^\\s*[+-]?(NaN|Infinity|((\\d+\\.?\\d*)|(\\.\\d+))([eE][+-]?\\d+)?[dDfF]?)\\s*$"));if(!le.test(c)){var d=new Rd;d.i="L`"+W(c);xd(d);wd(d,Error(d));throw d.h;}return parseFloat(c)}catch(f){var e=Ad(f);if(e instanceof Rd)return NaN;throw e.h;}},ne=function(a,b){ke(a,b);if(!je(a,b))return"";a=a.get(b);if(a==null)return"";var c;if(b="number"===typeof a&&(c=a,!0))b=Kd(c).equals(Kd(c));var d;b?d=""+Kd(c):d=W(a);return d},

ke=function(a,b){if(a.i){try{var c=a.h.get()[b]}catch(h){var d=Ad(h);if(d instanceof Cd)c="injection-failed";else throw d.h;}try{var e=a.i;if(e==null)throw yd().h;var f=e.get()[b]}catch(h){var g=Ad(h);if(g instanceof Cd)f="injection-failed";else throw g.h;}a=c;!(b=Qd(a,f))&&(b=a!=null)&&(b=a.equals?a.equals(f):Object.is(a,f));if(!b)throw Nd("Q").h;}},oe=function(){function a(){e[0]=1732584193;e[1]=4023233417;e[2]=2562383102;e[3]=271733878;e[4]=3285377520;m=l=0}function b(p){for(var q=g,n=0;n>>31)&4294967295;p=e[0];var u=e[1],G=e[2],E=e[3],da=e[4];for(n=0;n>>27)&4294967295)+A+da+O+q[n]&4294967295;da=E;E=G;G=(u>>2)&4294967295;u=p;p=A}e[0]=e[0]+p&4294967295;e[1]=e[1]+u&4294967295;e[2]=e[2]+G&4294967295;e[3]=e[3]+E&4294967295;

e[4]=e[4]+da&4294967295}function c(p,q){if(typeof p==="string"){p=unescape(encodeURIComponent(p));for(var n=[],u=0,G=p.length;u=56;n--)f[n]=q&255,q>>>=8;b(f);for(n=q=0;n=0;u-=8)p[q++]=e[n]>>u&255;return p}for(var e=

[],f=[],g=[],h=[128],k=1;k=0){var f=a[c].substring(0,d);e=a[c].substring(d+1)}else f=a[c];b(f,e?decodeURIComponent(e.replace(/\+/g," ")):"")}}},Ce=function(a,b,c,d,e,f,g,h,k){this.v=a;this.h=b;this.C=c;this.D=d;this.i=e;this.j=f;this.l=g;

this.o=h;this.u=k},Ee=function(){var a=he(),b=a.get("ilcm");if(b==null)return null;var c=b.eui,d=b.je,e=b.sstu;if(De)var f=De;else{f=he();var g=f.get("ilcm");g==null?f=null:De=f=X(f,"icso")?Vd():g.si}g=b.ei;a=a.get("buildLabel");return new Ce(c,d,e,f,g,a,b.crc||0,b.cvi||[],b.gsc||null)},Fe=function(a){this.m=I(a)},Ge=function(a){this.m=I(a)},He=function(a){this.m=I(a)},Ie=function(a){this.m=I(a)},Je=function(a,b,c){Y.call(this);this.u=a;this.i=typeof c==="number"?c:null;this.o=(a=Ee())?a.h:0;var d;

this.j=(d=a==null?void 0:a.i)!=null?d:[];this.h=null;this.l=b},Ke=ba("i"),Le=function(a){this.m=I(a)},Me=function(a){this.m=I(a)},Ne=function(a){this.m=I(a,4)},Oe=function(a){this.m=I(a,37)},Pe=function(a,b){return P(a,8,b)},Qe=ba("i"),Re=function(a){Ra.call(this);this.h=a},Te=function(a){var b=b===void 0?3E4:b;this.o=a;this.h=this.l=this.i=0;this.j=b;for(a=Se;a1)));g=g.next)e||(f=g);e&&(c.h==0&&d==1?Hf(c,b):(f?(d=f,d.next==c.l&&(c.l=d),d.next=

d.next.next):If(c),Jf(c,e,3,b)))}a.j=null}else wf(a,3,b)},Ff=function(a,b){a.i||a.h!=2&&a.h!=3||Kf(a);a.l?a.l.next=b:a.i=b;a.l=b},Mf=function(a,b,c,d){var e=Af(null,null,null);e.h=new xf(function(f,g){e.o=b?function(h){try{var k=b.call(d,h);f(k)}catch(l){g(l)}}:f;e.i=c?function(h){try{var k=c.call(d,h);k===void 0&&h instanceof Lf?g(h):f(k)}catch(l){g(l)}}:g});e.h.j=a;Ff(a,e);return e.h},wf=function(a,b,c){a.h==0&&(a===c&&(b=3,c=new TypeError("U")),a.h=1,Bf(c,a.Ra,a.Sa,a)||(a.u=c,a.h=b,a.j=null,Kf(a),

b!=3||c instanceof Lf||Nf(a,c)))},Bf=function(a,b,c,d){if(a instanceof xf)return Ff(a,Af(b||vf,c||null,d)),!0;if(a)try{var e=!!a.$goog_Thenable}catch(g){e=!1}else e=!1;if(e)return a.then(b,c,d),!0;if(Ga(a))try{var f=a.then;if(typeof f==="function")return Of(a,f,b,c,d),!0}catch(g){return c.call(d,g),!0}return!1},Of=function(a,b,c,d,e){function f(k){h||(h=!0,d.call(e,k))}function g(k){h||(h=!0,c.call(e,k))}var h=!1;try{b.call(a,g,f)}catch(k){f(k)}},Kf=function(a){a.v||(a.v=!0,df(a.Ha,a))},If=function(a){var b=

null;a.i&&(b=a.i,a.i=b.next,b.next=null);a.i||(a.l=null);return b},Jf=function(a,b,c,d){if(c==3&&b.i&&!b.l)for(;a&&a.o;a=a.j)a.o=!1;if(b.h)b.h.j=null,Pf(b,c,d);else try{b.l?b.o.call(b.j):Pf(b,c,d)}catch(e){Qf.call(null,e)}We(zf,b)},Pf=function(a,b,c){b==2?a.o.call(a.j,c):a.i&&a.i.call(a.j,c)},Nf=function(a,b){a.o=!0;df(function(){a.o&&Qf.call(null,b)})},Lf=function(a){Ra.call(this,a)},Rf=function(a,b){this.type=a;this.h=this.target=b;this.defaultPrevented=!1},Sf=function(a,b){Rf.call(this,a?a.type:

"");this.relatedTarget=this.h=this.target=null;this.button=this.screenY=this.screenX=this.clientY=this.clientX=0;this.key="";this.metaKey=this.shiftKey=this.altKey=this.ctrlKey=!1;this.state=null;this.pointerId=0;this.pointerType="";this.j=null;a&&this.init(a,b)},Uf=function(a,b,c,d,e){this.listener=a;this.proxy=null;this.src=b;this.type=c;this.capture=!!d;this.handler=e;this.key=++Tf;this.aa=this.ja=!1},Vf=function(a){a.aa=!0;a.listener=null;a.proxy=null;a.src=null;a.handler=null},Wf=function(a,

b,c){for(var d in a)b.call(c,a[d],d,a)},Yf=function(a,b){for(var c,d,e=1;e=0)&&Array.prototype.splice.call(d,e,1);f&&(Vf(b),a.h[c].length==0&&(delete a.h[c],a.i--))}},bg=function(a,b,c,d){for(var e=0;e-1&&(Vf(f[c]),Array.prototype.splice.call(f,c,1),f.length==0&&(delete a.h[b],a.i--)))):a&&(a=hg(a))&&(b=a.h[b.toString()],a=-1,b&&(a=bg(b,c,d,e)),(c=a>-1?b[a]:null)&&pg(c))},pg=function(a){if(typeof a!=="number"&&a&&!a.aa){var b=a.src;if(b&&b[fg])ag(b.h,a);else{var c=a.type,d=a.proxy;b.removeEventListener?b.removeEventListener(c,d,a.capture):b.detachEvent?

b.detachEvent(lg(c),d):b.addListener&&b.removeListener&&b.removeListener(d);mg--;(c=hg(b))?(ag(c,a),c.i==0&&(c.src=null,b[ig]=null)):Vf(a)}}},lg=function(a){return a in qg?qg[a]:qg[a]="on"+a},ng=function(a,b){if(a.aa)a=!0;else{b=new Sf(b,this);var c=a.listener,d=a.handler||a.src;a.ja&&pg(a);a=c.call(d,b)}return a},hg=function(a){a=a[ig];return a instanceof Zf?a:null},eg=function(a){if(typeof a==="function")return a;a[rg]||(a[rg]=function(b){return a.handleEvent(b)});return a[rg]},sg=function(){Y.call(this);

this.h=new Zf(this);this.D=this;this.l=null},ug=function(a,b){var c=a.l;if(c){var d=[];for(var e=1;c;c=c.l)d.push(c),++e}a=a.D;c=b.type||b;typeof b==="string"?b=new Rf(b,a):b instanceof Rf?b.target=b.target||a:(e=b,b=new Rf(c,a),Yf(b,e));e=!0;var f;if(d)for(f=d.length-1;f>=0;f--){var g=b.h=d[f];e=tg(g,c,!0,b)&&e}g=b.h=a;e=tg(g,c,!0,b)&&e;e=tg(g,c,!1,b)&&e;if(d)for(f=0;f2147483647?-1:z.setTimeout(a,b||0)},yg=function(a){var b=

null;return(new xf(function(c,d){b=xg(function(){c(void 0)},a);b==-1&&d(Error("Y"))})).na(function(c){z.clearTimeout(b);throw c;})},zg=function(a,b){this.h=a;this.i=b;a=this.h;b=me(this.i,"docs-clibs");a.U=b;this.h.X=2E4},Cg=function(a,b,c){var d=Pe(new Oe,bd(b));Ag(a.h,d);return(new xf(function(e,f){var g=Date.now();c.i++;c.l=g;Bg(a,e,f)})).na(function(e){if(typeof e==="number"&&(500>4&15).toString(16)+(a&15).toString(16)},Mg=function(a,b){this.i=this.h=null;this.j=a||null;this.l=!!b},Rg=function(a){a.h||(a.h=new Map,a.i=0,a.j&&Be(a.j,function(b,c){a.add(decodeURIComponent(b.replace(/\+/g," ")),c)}))},Tg=function(a,b){Rg(a);b=Sg(a,b);a.h.has(b)&&(a.j=null,a.i=a.i-a.h.get(b).length,a.h.delete(b))},Ug=function(a,b){Rg(a);b=Sg(a,b);return a.h.has(b)},Jg=function(a){var b=new Mg;

b.j=a.j;a.h&&(b.h=new Map(a.h),b.i=a.i);return b},Sg=function(a,b){b=String(b);a.l&&(b=b.toLowerCase());return b},Ng=function(a,b){b&&!a.l&&(Rg(a),a.j=null,a.h.forEach(function(c,d){var e=d.toLowerCase();if(d!=e&&(Tg(this,d),Tg(this,e),c.length>0)){this.j=null;d=this.h;var f=d.set;e=Sg(this,e);var g=c.length;if(g>0){for(var h=Array(g),k=0;k0?d:void 0);d=ad(d,4,f>0?f:void 0);d=ad(d,5,g>0?

g:void 0);d=Bc(d);M(h,jh,10,d)}a=kd(a.h);h=Date.now().toString();a=Hc(a,4,h==null?h:Vb(h));b=Sc(a,Oe,3,b.slice());e&&(a=new Le,e=ad(a,13,e),a=new Me,e=M(a,Le,2,e),a=new Ne,e=M(a,Me,1,e),e=T(e,2,9),M(b,Ne,18,e));c&&N(b,14,c);return b},Ch=function(a){this.i=this.h=this.j=a},ld=function(a){this.m=I(a,8)},md=function(a){this.m=I(a)},Gh=function(a){Y.call(this);var b=this;this.h=[];this.ia="";this.ca=!1;this.va=this.O=-1;this.G=null;this.u=this.l=0;this.L=null;this.M=this.N=0;this.Ca=1;this.X=0;this.W=

a.W;this.V=a.V||r();this.j=new wh(a.W,a.J);this.F=a.F||null;this.S=a.S||null;this.U=1E3;this.D=a.Ta||null;this.R=a.R||null;this.Y=a.Y||!1;this.withCredentials=!a.xa;this.J=a.J||!1;this.ha=!this.J&&!!window&&!!window.navigator&&window.navigator.sendBeacon!==void 0;this.ga=typeof URLSearchParams!=="undefined"&&!!(new URL(Dh())).searchParams&&!!(new URL(Dh())).searchParams.set;var c=th(new sh);vh(this.j,c);this.o=new Ch(1E4);a=Eh(this,a.wa);this.i=new eh(this.o.h,a);this.da=new eh(6E5,a);this.Y||this.da.start();

if(!this.J){document.addEventListener("visibilitychange",function(){if(document.visibilityState==="hidden"){Fh(b);var f;(f=b.L)==null||f.flush()}});var d,e;(d=window)==null||(e=d.addEventListener)==null||e.call(d,"pagehide",function(){Fh(b);var f;(f=b.L)==null||f.flush()})}},Eh=function(a,b){function c(){a.flush()}return a.ga?b?function(){b().then(c)}:c:r()},Hh=function(a){a.D||(a.D=Dh());try{return(new URL(a.D)).toString()}catch(b){return(new URL(a.D,window.location.origin)).toString()}},Ag=function(a,

b){if(b instanceof Oe)a.log(b);else try{var c=Pe(new Oe,bd(b));a.log(c)}catch(d){Ih(a,4,1)}},Ih=function(a,b,c){var d;(d=a.L)==null||d.lb(b,c)},Jh=function(a,b,c){c=c===void 0?null:c;var d=d===void 0?a.withCredentials:d;var e={},f=new URL(Hh(a));c&&(e.Authorization=c);a.R&&(e["X-Goog-AuthUser"]=a.R,f.searchParams.set("authuser",a.R));return{url:f.toString(),body:b,Fa:1,ua:e,Oa:"POST",withCredentials:d,X:a.X}},Fh=function(a){a.j.l=!0;a.ca&&(a.j.i=3,Kh(a));a.flush();a.j.l=!1},Kh=function(a){Lh(a,function(b,

c){b=new URL(b);b.searchParams.set("format","json");var d=!1;try{d=window.navigator.sendBeacon(b.toString(),bd(c))}catch(e){}d||(a.ha=!1);return d})},Lh=function(a,b){if(a.h.length!==0){var c=new URL(Hh(a));c.searchParams.delete("format");var d=a.V();d&&c.searchParams.set("auth",d);c.searchParams.set("authuser",a.R||"0");for(d=0;d=0;--c)a.push(typeof b[c],b[c]);return a.join("\v")},mj=function(a){sg.call(this);a||(a=lj||(lj=new ih));this.i=a;if(this.j=this.Ja())this.o=

dg(this.i.h,this.j,Ma(this.La,this))},nj=function(a,b){Y.call(this);this.i=a;this.h=new mj(b);Ae(this,this.h);this.j=new Yi(this);Ae(this,this.j);this.h.Z()&&$i(this.j,this.h,"visibilitychange",this.l)},oj=function(a,b,c){c=c===void 0?!1:c;Y.call(this);this.h=a;this.i=b;Ae(this,this.i);this.j=c},pj=function(a,b,c){Y.call(this);this.u=c!=null?a.bind(c):a;this.o=b;this.i=null;this.l=!1;this.h=null},qj=function(a){a.h=xg(function(){a.h=null;a.l&&(a.l=!1,qj(a))},a.o);var b=a.i;a.i=null;a.u.apply(null,

b)},rj=function(a,b,c,d,e){Y.call(this);this.h=a;this.D=b;this.i=new pj(this.l,3E3,this);this.j=new Set;this.o=d;this.u=e||6E4},sj=function(a){Y.call(this);this.h=a;Ae(this,this.h)},tj=function(a){this.m=I(a)},uj=function(a){this.m=I(a)},vj=function(a){this.m=I(a)},wj=function(a){this.m=I(a)},xj=function(a){this.m=I(a)},yj=function(a){this.m=I(a)},zj=function(a){this.m=I(a)},Aj=function(a){this.m=I(a)},Bj=function(a){this.m=I(a)},Cj=function(a){this.m=I(a)},Ej=function(){var a=void 0,b=this;a=a===

void 0?!0:a;this.o=new Dj;this.j=this.i=0;(this.l=a)&&document.addEventListener("mousemove",function(c){b.i=c.clientX;b.j=c.clientY})},Dj=function(){this.i=this.h=null},Gj=function(a){if(a.h!==null)return a.h;if(a.i===null&&(a.i=Fj(),a.i===null))return null;a.h=a.i.chrome.runtime.connect("ciiamoeghmklpofjbocenebdfbgjapaa");a.h.onDisconnect.addListener(function(){a.h=null});return a.h},Fj=function(){for(var a=0;a=0;f--){var g=L(c[f],$g,5);if(g&&L(g,ee,1)){g=L(g,ee,1);if(Uc(g,12)!=null&&l===void 0){var h=void 0,k=void 0;k=k===void 0?!1:k;var l=(h=Uc(g,12))!=null?h:k}g=L(g,ce,20);g!==void 0&&e===void 0&&(e=new Fe,h=Uc(g,2,Si),h!==void 0&&$c(e,2,h),g=Uc(g,1,Si),g!==void 0&&$c(e,1,g));if(l!==void 0&&e!==void 0)break}}d=d.l?kd(d.l):null;if(l!==void 0||e!==void 0)d||(d=new Ie),l!==void 0&&$c(d,6,l),e!==void 0&&M(d,Fe,13,e);(l=d)&&M(b,

Ie,3,l);a=kd(a.i.o);M(b,Nj,4,a);Sc(b,Mj,1,c);return b},Vj=function(a,b){var c=dj(Rj(a,b),0,a.v++);var d=a.u;var e=Object.keys(d.h);if(e.length==0)d=null;else{for(var f=[],g=0;g0;){var a=this.h.pop();if(a in this.j)return a}return null};ra.prototype.getNext=ra.prototype.i;v("Symbol",function(a){function b(f){if(this instanceof b)throw new TypeError("h");return new c(d+(f||"")+"_"+e++,f)}function c(f,g){this.h=f;ha(this,"description",{configurable:!0,writable:!0,value:g})}if(a)return a;c.prototype.toString=function(){return this.h};var d="jscomp_symbol_"+(Math.random()*1E9>>>0)+"_",e=0;return b});

v("Symbol.iterator",function(a){if(a)return a;a=Symbol("i");ha(Array.prototype,a,{configurable:!0,writable:!0,value:function(){return Aa(ka(this))}});return a});

v("Promise",function(a){function b(g){this.h=0;this.j=void 0;this.i=[];this.u=!1;var h=this.l();try{g(h.resolve,h.reject)}catch(k){h.reject(k)}}function c(){this.h=null}function d(g){return g instanceof b?g:new b(function(h){h(g)})}if(a)return a;c.prototype.i=function(g){if(this.h==null){this.h=[];var h=this;this.j(function(){h.o()})}this.h.push(g)};var e=fa.setTimeout;c.prototype.j=function(g){e(g,0)};c.prototype.o=function(){for(;this.h&&this.h.length;){var g=this.h;this.h=[];for(var h=0;h=f}});v("Symbol.dispose",function(a){return a?a:Symbol("o")});

v("WeakMap",function(a){function b(k){this.h=(h+=Math.random()+1).toString();if(k){k=x(k);for(var l;!(l=k.next()).done;)l=l.value,this.set(l[0],l[1])}}function c(){}function d(k){var l=typeof k;return l==="object"&&k!==null||l==="function"}function e(k){if(!na(k,g)){var l=new c;ha(k,g,{value:l})}}function f(k){var l=Object[k];l&&(Object[k]=function(m){if(m instanceof c)return m;Object.isExtensible(m)&&e(m);return l(m)})}if(function(){if(!a||!Object.seal)return!1;try{var k=Object.seal({}),l=Object.seal({}),

m=new a([[k,2],[l,3]]);if(m.get(k)!=2||m.get(l)!=3)return!1;m.delete(k);m.set(l,4);return!m.has(k)&&m.get(l)==4}catch(p){return!1}}())return a;var g="$jscomp_hidden_"+Math.random();f("freeze");f("preventExtensions");f("seal");var h=0;b.prototype.set=function(k,l){if(!d(k))throw Error("p");e(k);if(!na(k,g))throw Error("q`"+k);k[g][this.h]=l;return this};b.prototype.get=function(k){return d(k)&&na(k,g)?k[g][this.h]:void 0};b.prototype.has=function(k){return d(k)&&na(k,g)&&na(k[g],this.h)};b.prototype.delete=

function(k){return d(k)&&na(k,g)&&na(k[g],this.h)?delete k[g][this.h]:!1};return b});

v("Map",function(a){function b(){var h={};return h.K=h.next=h.head=h}function c(h,k){var l=h[1];return Aa(function(){if(l){for(;l.head!=h[1];)l=l.K;for(;l.next!=l.head;)return l=l.next,{done:!1,value:k(l)};l=null}return{done:!0,value:void 0}})}function d(h,k){var l=k&&typeof k;l=="object"||l=="function"?f.has(k)?l=f.get(k):(l=""+ ++g,f.set(k,l)):l="p_"+k;var m=h[0][l];if(m&&na(h[0],l))for(h=0;h1342177279)throw new RangeError("r");b|=0;for(var d="";b;)if(b&1&&(d+=c),b>>>=1)c+=c;return d}});var z=this||self,xk=z._F_toggles__pubi||[],Ha="closure_uid_"+(Math.random()*1E9>>>0),Ia=0;Qa(Ra,Error);Ra.prototype.name="CustomError";var lj;var ve=String.prototype.trim?function(a){return a.trim()}:function(a){return/^[\s\xa0]*([\s\S]*?)[\s\xa0]*$/.exec(a)[1]};var yk=!!(xk[0]>>24&1),zk=!!(xk[0]>>19&1),Ak=!!(xk[0]>>26&1),Bk=!!(xk[0]&4096);var Ua=yk?Ak:Fa(610401301,!1),tc=yk?zk||!Bk:Fa(748402147,!0);var Va,Ck=z.navigator;Va=Ck?Ck.userAgentData||null:null;var $f=Array.prototype.indexOf?function(a,b){return Array.prototype.indexOf.call(a,b,void 0)}:function(a,b){if(typeof a==="string")return typeof b!=="string"||b.length!=1?-1:a.indexOf(b,0);for(var c=0;cparseFloat(ll)){Zk=String(Gl);break a}}Zk=ll}var Il=Zk;var lk=$a(),Jl=cb()||C("iPod"),Kl=C("iPad"),Ll=C("Android")&&!(ab()||$a()||Ya()||C("Silk")),Rh=ab(),Ml=C("Safari")&&!(ab()||(Xa()?0:C("Coast"))||Ya()||(Xa()?0:C("Edge"))||Za()||(Xa()?Wa("Opera"):C("OPR"))||$a()||C("Silk")||C("Android"))&&!db();var jb=void 0;var ob=typeof Symbol==="function"&&typeof Symbol()==="symbol",Nl=mb("jas",void 0,!0),mc=mb(void 0,Symbol()),Ol=mb(void 0,"0ub"),kc=mb(void 0,"0ubs"),vc=mb(void 0,"0actk"),fc=mb("m_m","ob",!0),Pl=mb();Math.max.apply(Math,ma(Object.values({bb:1,ab:2,Za:4,gb:8,ib:16,eb:32,Ua:64,Xa:128,Va:256,hb:512,Wa:1024,Ya:2048,fb:4096,cb:8192})));var qb={Ma:{value:0,configurable:!0,writable:!0,enumerable:!1}},pb=Object.defineProperties,D=ob?Nl:"Ma",Mc,Ql=[];F(Ql,7);Mc=Object.freeze(Ql);var hc={},sb={},Qc=Object.freeze({}),nd={};var wb=ub(function(a){return typeof a==="number"}),vb=ub(function(a){return typeof a==="string"}),yb=ub(function(a){return typeof a==="boolean"}),Rl=ub(function(a){return typeof a==="bigint"});var xb=typeof z.BigInt==="function"&&typeof z.BigInt(0)==="bigint";var $d=ub(function(a){return xb?Rl(a):vb(a)&&/^(?:-?[1-9]\d*|0)$/.test(a)}),pc=ub(function(a){return xb?a>=Sl&&a>>0).toString(16))};w(ud,td);ud.prototype.toString=function(){var a=Ud(Td(this.constructor)),b=this.i;return b==null?a:W(a)+": "+W(b)};w(Bd,ud);w(Cd,Bd);t=Dd.prototype;t.isSafeInteger=function(){var a=this.A>>21;return a==0||a==-1&&!(this.B==0&&this.A==-2097152)};

t.toString=function(a){a=a||10;if(a>2);var c=Math.pow(a,b),d=V(c,c/4294967296);c=this.div(d);var e=Math,f=e.abs;d=c.multiply(d);d=this.add(Gd(d));e=f.call(e,Ed(d));f=a==10?""+e:e.toString(a);f.length>>0>a.B>>>0?1:-1:this.A>a.A?1:-1};t.add=function(a){var b=this.A>>>16,c=this.A&65535,d=this.B>>>16,e=a.A>>>16,f=a.A&65535,g=a.B>>>16;a=(this.B&65535)+(a.B&65535);g=(a>>>16)+(d+g);d=g>>>16;d+=c+f;return V((g&65535)>>16)+(b+e)&65535)>>16,c=this.A&65535,d=this.B>>>16,e=this.B&65535,f=a.A>>>16,g=a.A&65535,h=a.B>>>16;a=a.B&65535;var k=e*a;var l=(k>>>16)+d*a;var m=l>>>16;l=(l&65535)+e*h;m+=l>>>16;m+=c*a;var p=m>>>16;m=(m&65535)+d*h;p+=m>>>16;m=(m&65535)+e*g;p=p+(m>>>16)+(b*a+c*h+d*g+e*f)&65535;return V((l&65535)>>1|b>1);b=b.div(a).shiftLeft(1);if(b.equals(Jd))return a.A=0;){var d=Math.max(1,Math.floor(Ed(c)/Ed(a))),

e=Math.ceil(Math.log(d)/Math.LN2);e=e0;)d-=e,f=Kd(d),g=f.multiply(a);Fd(f)&&(f=Yl);b=b.add(f);c=c.add(Gd(g))}return b};t.and=function(a){return V(this.B&a.B,this.A&a.A)};t.or=function(a){return V(this.B|a.B,this.A|a.A)};t.xor=function(a){return V(this.B^a.B,this.A^a.A)};t.shiftLeft=function(a){a&=63;if(a==0)return this;var b=this.B;return a>>32-a):V(0,b=0;b--){var c=a[b];this.get(c);this.set(c,"",{Na:0,path:void 0,domain:void 0})}};Y.prototype.C=!1;Y.prototype.dispose=function(){this.C||(this.C=!0,this.H())};Y.prototype[Symbol.dispose]=function(){this.dispose()};Y.prototype.H=function(){if(this.v)for(;this.v.length;)this.v.shift()()};var Kg=RegExp("^(?:([^:/?#.]+):)?(?://(?:([^\\\\/?#]*)@)?([^\\\\/?#]*?)(?::([0-9]+))?(?=[\\\\/?#]|$))?([^?#]+)?(?:\\?([^#]*))?(?:#([\\s\\S]*))?$");var De=null;w(Fe,U);w(Ge,U);w(He,U);w(Ie,U);w(Je,Y);

Je.prototype.get=function(){if(this.h)return this.h;var a=new Ie;a=P(a,1,"en");a=P(a,2,Ta());typeof this.i==="number"&&ad(a,11,this.i);var b=new He;b=$c(b,2,this.u);var c=X(this.l,"icso");b=$c(b,1,c);M(a,He,5,b);T(a,9,this.o);b=new Ge;c=this.j;Dc(b);var d=b.m;var e=d[D]|0,f=tb(b,e)?1:2;f===2&&Cc(b)&&(d=b.m,e=d[D]|0);var g=Nc(d,1),h=g===Mc?7:g[D]|0,k=Oc(h,e);var l=4&k?!1:!0;if(l){4&k&&(g=Ib(g),h=0,k=Kc(k,e),e=Gc(d,e,1,g));for(var m=0,p=0;m0){this.i--;var a=this.h;this.h=a.next;a.next=null}else a=this.j();return a};Xe.prototype.add=function(a,b){var c=uf.get();c.set(a,b);this.i?this.i.next=c:this.h=c;this.i=c};var uf=new Ve(function(){return new $e},function(a){return a.reset()});$e.prototype.set=function(a,b){this.i=a;this.h=b;this.next=null};$e.prototype.reset=function(){this.next=this.h=this.i=null};var af,cf=!1,Ye=new Xe;yf.prototype.reset=function(){this.j=this.i=this.o=this.h=null;this.l=!1};var zf=new Ve(function(){return new yf},function(a){a.reset()});xf.prototype.then=function(a,b,c){return Mf(this,Ef(typeof a==="function"?a:null),Ef(typeof b==="function"?b:null),c)};xf.prototype.$goog_Thenable=!0;t=xf.prototype;t.na=function(a,b){return Mf(this,null,Ef(a),b)};t.catch=xf.prototype.na;t.cancel=function(a){if(this.h==0){var b=new Lf(a);df(function(){Hf(this,b)},this)}};t.Ra=function(a){this.h=0;wf(this,2,a)};

t.Sa=function(a){this.h=0;wf(this,3,a)};t.Ha=function(){for(var a;a=If(this);)Jf(this,a,this.h,this.u);this.v=!1};var Qf=Sa;Qa(Lf,Ra);Lf.prototype.name="cancel";Rf.prototype.i=function(){this.defaultPrevented=!0};var jk=!!(z.navigator&&z.navigator.maxTouchPoints||"ontouchstart"in z||z.document&&document.documentElement&&"ontouchstart"in document.documentElement||z.navigator&&z.navigator.msMaxTouchPoints),kk="PointerEvent"in z,kg=function(){if(!z.addEventListener||!Object.defineProperty)return!1;var a=!1,b=Object.defineProperty({},"passive",{get:function(){a=!0}});try{var c=r();z.addEventListener("test",c,b);z.removeEventListener("test",c,b)}catch(d){}return a}();Qa(Sf,Rf);

Sf.prototype.init=function(a,b){var c=this.type=a.type,d=a.changedTouches&&a.changedTouches.length?a.changedTouches[0]:null;this.target=a.target||a.srcElement;this.h=b;b=a.relatedTarget;b||(c=="mouseover"?b=a.fromElement:c=="mouseout"&&(b=a.toElement));this.relatedTarget=b;d?(this.clientX=d.clientX!==void 0?d.clientX:d.pageX,this.clientY=d.clientY!==void 0?d.clientY:d.pageY,this.screenX=d.screenX||0,this.screenY=d.screenY||0):(this.clientX=a.clientX!==void 0?a.clientX:a.pageX,this.clientY=a.clientY!==

void 0?a.clientY:a.pageY,this.screenX=a.screenX||0,this.screenY=a.screenY||0);this.button=a.button;this.key=a.key||"";this.ctrlKey=a.ctrlKey;this.altKey=a.altKey;this.shiftKey=a.shiftKey;this.metaKey=a.metaKey;this.pointerId=a.pointerId||0;this.pointerType=a.pointerType;this.state=a.state;this.j=a;a.defaultPrevented&&Sf.T.i.call(this)};Sf.prototype.i=function(){Sf.T.i.call(this);var a=this.j;a.preventDefault?a.preventDefault():a.returnValue=!1};var fg="closure_listenable_"+(Math.random()*1E6|0);var Tf=0;var Xf="constructor hasOwnProperty isPrototypeOf propertyIsEnumerable toLocaleString toString valueOf".split(" ");Zf.prototype.add=function(a,b,c,d,e){var f=a.toString();a=this.h[f];a||(a=this.h[f]=[],this.i++);var g=bg(a,b,d,e);g>-1?(b=a[g],c||(b.ja=!1)):(b=new Uf(b,this.src,f,!!d,e),b.ja=c,a.push(b));return b};var ig="closure_lm_"+(Math.random()*1E6|0),qg={},mg=0,rg="__closure_events_fn_"+(Math.random()*1E9>>>0);Qa(sg,Y);sg.prototype[fg]=!0;sg.prototype.addEventListener=function(a,b,c,d){dg(this,a,b,c,d)};sg.prototype.removeEventListener=function(a,b,c,d){og(this,a,b,c,d)};sg.prototype.H=function(){sg.T.H.call(this);if(this.h){var a=this.h,b=0,c;for(c in a.h){for(var d=a.h[c],e=0;e0&&a1||f.length==1&&

f[0]!="")&&f.pop(),d&&g==e.length&&f.push("")):(f.push(h),d=!0)}d=f.join("/")}else d=e}c?b.i=d:c=a.u.toString()!=="";c?Ig(b,Jg(a.u)):c=!!a.v;c&&(b.v=a.v);return b};var $l=/[#\/\?@]/g,bm=/[#\?:]/g,am=/[#\?]/g,Pg=/[#\?@]/g,cm=/#/g;t=Mg.prototype;t.add=function(a,b){Rg(this);this.j=null;a=Sg(this,a);var c=this.h.get(a);c||this.h.set(a,c=[]);c.push(b);this.i=this.i+1;return this};t.clear=function(){this.h=this.j=null;this.i=0};

t.forEach=function(a,b){Rg(this);this.h.forEach(function(c,d){c.forEach(function(e){a.call(b,e,d,this)},this)},this)};t.ta=function(a){Rg(this);var b=[];if(typeof a==="string")Ug(this,a)&&(b=b.concat(this.h.get(Sg(this,a))));else{a=Array.from(this.h.values());for(var c=0;c0?String(a[0]):b};t.toString=function(){if(this.j)return this.j;if(!this.h)return"";for(var a=[],b=Array.from(this.h.keys()),c=0;c0&&(this.h.splice(0,b),this.l+=b,Ih(this,3,b));this.h.push(a);this.Y||this.i.i||this.i.start()}};

Gh.prototype.flush=function(a,b){var c=this;if(this.h.length===0)a&&a();else{var d=Date.now();if(this.va>d&&this.O0&&(c.O=Date.now(),c.va=c.O+q);q=em.h;u=Pa(mc);var G;ob&&u&&((G=n.m[u])==null?void 0:G[q])!=null&&kb(Ol,3);a:{G=em.h;var E=E===void 0?!1:E;if(Pa(Pl)&&Pa(mc)&&void 0===Pl){q=n.m;u=q[mc];if(!u)break a;if(u=u.qb)try{u(q,G,Wl);break a}catch(A){Sa(A)}}E&&(E=n.m,(q=Pa(mc))&&q in E&&(E=E[q])&&delete E[G])}E=em.i?em.j(n,em.i,em.h,em.l):em.j(n,em.h,null,em.l);if(n=E===null?void 0:E){E=-1;

E=E===void 0?0:E;var da;n=(da=Ob(J(n,1)))!=null?da:E;n!==-1&&(c.o=new Ch(n0?setTimeout(function(){e.abort()},a.X):void 0,u.O(2,3),g=Object.assign({},{method:a.Oa,headers:Object.assign({},a.ua)},a.body&&{body:a.body},a.withCredentials&&{credentials:"include"},{signal:a.X&&e?e.signal:null}),u.o(fetch(a.url,g),5);case 5:h=u.u;if(h.status!==200){(k=c)==null||k(h.status);u.P(3);break}if((l=b)==null){u.P(7);

break}return u.o(h.text(),8);case 8:l(u.u);case 7:case 3:u.M();clearTimeout(f);u.N(0);break;case 2:m=u.L();switch((p=m)==null?void 0:p.name){case "AbortError":(q=c)==null||q(408);break;default:(n=c)==null||n(400)}u.P(3)}})))};Oh.prototype.sa=ca(4);w(Ph,Y);Ph.prototype.xa=function(){this.o=!0;return this};(function(){if(lk)return Qh(/Firefox\/([0-9.]+)/);if(Ek||Fk||Dk)return Il;if(Rh){if(db()||eb()){var a=Qh(/CriOS\/([0-9.]+)/);if(a)return a}return Qh(/Chrome\/([0-9.]+)/)}if(Ml&&!db())return Qh(/Version\/([0-9.]+)/);if(Jl||Kl){if(a=/Version\/(\S+).*Mobile\/(\S+)/.exec(Ta()))return a[1]+"."+a[2]}else if(Ll)return(a=Qh(/Android\s+([0-9.]+)/))?a:Qh(/Version\/([0-9.]+)/);return""})();Xh.prototype.toString=function(){var a=Yh(this);if(a===null)throw Error("ea`K1cgmc");return a};w(ki,U);w(li,U);w(mi,U);w(zi,U);var ik=function(){if(dk){var a=/Windows NT ([0-9.]+)/;return(a=a.exec(Ta()))?a[1]:"0"}return ek?(a=/1[0|1][_.][0-9_.]+/,(a=a.exec(Ta()))?a[0].replace(/_/g,"."):"10"):gk?(a=/Android\s+([^\);]+)(\)|;)/,(a=a.exec(Ta()))?a[1]:""):Mk||Xk||Yk?(a=/(?:iPhone|CPU)\s+OS\s+(\S+)/,(a=a.exec(Ta()))?a[1].replace(/_/g,"."):""):""}();w(Pi,U);w(Qi,U);w(Ri,U);w(Xi,Ui);var ck=new Oi("high_frequency_builder");Qa(Yi,Y);var Zi=[];Yi.prototype.H=function(){Yi.T.H.call(this);aj(this)};Yi.prototype.handleEvent=function(){throw Error("ga");};w(cj,Ui);var Qj=new Oi("system_builder");w(hj,Rf);var jj=new WeakMap;Qa(mj,sg);t=mj.prototype;t.Ja=kj(function(){var a=this.Z(),b=this.la()!="hidden";if(a){var c;b?c=(((ej?"Webkit":fj?"Moz":null)||"")+"visibilitychange").toLowerCase():c="visibilitychange";a=c}else a=null;return a});t.la=kj(function(){return gj("hidden",this.i.h)});t.Ka=kj(function(){return gj("visibilityState",this.i.h)});t.Z=function(){return!!this.la()};t.La=function(){var a=this.Z()?this.i.h[this.Ka()]:null;a=new hj(!!this.i.h[this.la()],a);ug(this,a)};t.H=function(){pg(this.o);mj.T.H.call(this)};w(nj,Y);nj.prototype.l=function(){if(this.i.pa()){var a=this.h;a=!!a.i.h[a.la()];a=this.i.fa(a?102001:102E3,0);this.i.ma(a)}};w(oj,Y);t=oj.prototype;t.ma=function(a){var b=this.h;N(a.h,6,b.j);b.l=!0;a=Wi(a);b.h.add(a);b=this.i;b.h.h.h.length>=3&&b.i.j();return new Ni};t.fa=function(a,b){a=dj(Rj(this.h,a),b,this.h.v++);b==1&&(b=this.h,Ti(L(a.h,Ri,8)),b.u.add(a));return a};t.oa=function(){return this.h.i};t.qa=function(){var a=this.h,b=Vj(a,716);Uj(a,b);b=Wi(b);a.h.add(b);a.D=!0;a.C=!0;a=this.i;bj(a.u,a.i.j,a.i);bj(36E5,a.G,a);this.i.i.j();this.j&&new nj(this)};t.ra=function(){this.i.l();return Df(Array.from(this.i.j)).then()};

t.pa=function(){var a=this.h;return a.D&&a.C&&!0};w(pj,Y);pj.prototype.j=function(a){this.i=arguments;this.h?this.l=!0:qj(this)};pj.prototype.H=function(){Y.prototype.H.call(this);this.h&&(z.clearTimeout(this.h),this.h=null,this.l=!1,this.i=null)};w(rj,Y);rj.prototype.l=function(){var a=this;if(this.h.h.h.length!=0&&(!this.o||this.h.l)){var b=Tj(this.h),c=this.D.h(b);c&&(Gf(c,function(){return void a.j.delete(c)}),this.j.add(c))}};rj.prototype.G=function(){var a=this.h,b=Vj(a,1153);b=Wi(b);a.h.add(b);this.i.j()};w(sj,Y);t=sj.prototype;t.ma=function(a){a=this.h.ma(a);this.ra();return a};t.fa=function(a,b){return this.h.fa(a,b)};t.oa=function(){return this.h.oa()};t.qa=function(){return this.h.qa()};t.ra=function(){return this.h.ra()};t.pa=function(){return this.h.pa()};w(tj,U);w(uj,U);w(vj,U);w(wj,U);w(xj,U);w(yj,U);w(zj,U);w(Aj,U);w(Bj,U);w(Cj,U);var Hj=z.window?[z.window,z.window.opener,z.window.parent]:[];Ej.prototype.h=function(a){var b=Gj(this.o);if(b){var c=L(a,tj,2)||new tj,d=L(a,Cj,5)||new Cj,e=L(a,Ie,3)||new Ie;a=Rc(a,Mj,1);a=x(a);var f=a.next(),g;try{for(;!f.done;f=a.next())b.postMessage({detail:{impression:bd(f.value),session_info:bd(c),session_invariants:bd(d),client_info:bd(e)}})}finally{f&&!f.done&&(g=a.return)&&g.call(a)}this.l&&document.dispatchEvent(new CustomEvent("ripple",{detail:{clientX:this.i,clientY:this.j}}))}};Ij.prototype.za=function(){return new Xi};Kj.prototype.add=function(a){this.h.push(a)};Lj.prototype.add=function(a){Ti(L(a.h,Ri,8));var b=ae(Yc(a.h,12));this.h[b]=a};w(Mj,U);w(Nj,U);w(Oj,U);Wj.prototype.za=function(){return new cj};nk.prototype.start=function(){var a=z.DOCS_drawing_load,b=z.DOCS_drawing_decode,c;for(c in a)rk(this,c);for(var d in b)qk(this,d);(a=z.DOCS_timing)&&(a=a.ejl-a.sjl)&&ok(this,29035,a);ok(this,29031,performance.now()-this.h)};z.DOCS_initPublishImpressionTracker=function(a,b,c,d){var e=document.querySelectorAll("img[id^='ed.']").length;a=new mk(a,b,e,c,d);var f=new nk(a,e);Oa("DOCS_notifyDrawingLoad",function(g){return rk(f,g)});Oa("DOCS_notifyDrawingDecode",function(g){return qk(f,g)});f.start();return f};

})(this._pubi);

// Google Inc.

//# sourceMappingURL=publish_binary_core.sourcemap

DOCS_timing['ejl'] = performance.now(); DOCS_initPublishImpressionTracker( 142.0 , 1.0 , 1.0 ); DOCS_timing['epr'] = performance.now();

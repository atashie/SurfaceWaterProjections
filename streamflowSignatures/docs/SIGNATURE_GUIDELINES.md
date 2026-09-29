# Summary Documentation for Streamflow Signatures

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/d/e/2PACX-1vSVjtqLKk1r9TczxLEBhlnzfBWbm1TQVfvqERm-jEwLISZTEWx73ofV4Ng9H0JaXA/pub)
> **Last synced**: 2026-09-29 (277 extracted paragraphs vs 270 on 2026-09-10). Structure
> unchanged (Parts 1–6, eight Part-3 modules); 20 wording regions changed. Substantive:
> Part 1.2 now cites "Hecht et al. 2024" for the four QA conditions (not in Part 6) and
> drops the "tell us if there is a data-driven answer" request; Part 3.3 rewrites the
> b = 1 recession paragraph (same rule), cites Tashie et al. 2020a/b (both added to
> Part 6), names `recession_alpha_point_cloud` in the L-H param line, and changes the
> parameterized-BFI validity interval from (0, 1) to **[0, 1] — the code rejects
> alpha ≤ 0 or ≥ 1, so the open interval was right**; Part 3.4 adds "Threshold of 2%
> selected to remove error / noise"; Part 3.8 rewords four snow definitions (no rule
> change); **Part 4 rewritten** — entry-point intro, `Function:` headings, the compilation
> gate restated ("at least 20 calendar years each have at least 30 days with daily
> discharge above 0.0001 mm … applies only at compilation"), "All daily records in the
> retrieval window are kept …", the USGS qualifier mask / HYDAT-never-masked sentence
> (three 2026-09-04 relay items now landed), the 95 % coverage-warning bullet dropped;
> **Part 5 rewritten** as a `compute_qa_flags()` block — high_na redefined as ">30 % of a
> gage's signature columns (the 16 statistics of every signature plus the 21 per-gage
> scalars)", "confirmatory" BFI notes, and a new Requirements block ("recomputed from an
> output CSV without rerunning any signature"; "aligned 2026-09-04"). NOT landed: the
> released-product high_na known-issue sentence; "(see 3.12)"; avg_storage;
> ice_affected_days_total; the shipped `_linear_reservoir` column name. See
> `docs/reconciliation/guidelines_todos.md` (2026-09-29).

---

# Summary Documentation for Streamflow Signatures

GitHub repo here· Summary data audit sheet here

How this document is organized. Part 1 holds the global rules that apply to every signature (time conventions, NA handling, completeness requirements). Part 2 defines the 16 statistics computed for every signature. Part 3 documents each signature family as a self-contained module with a fixed structure — Metrics (name: definition, units), Function, Method and parameters, and Requirements and decisions — so a reader can take in one family without cross-referencing. Part 4 covers data ingestion and utility functions, Part 5 the automated quality flags, and Part 6 the references.

## Part 1: Global Rules (apply to every signature)

### 1.1 Time conventions

All annual metrics are calculated on the water year, October 1 – September 30.

Seasonal periods are defined as: Winter: December–February; Spring: March–May; Summer: June–August; Fall: September–November.

### 1.2 NA handling and per-year screening

Do not replace NAs with zeros; streamflow that is zero is okay.

NAs can be for various reasons, we need to use the U.S. Geological Survey (USGS) codes that describe what creates those NAs (e.g. Ice affected, etc.); for each site that has ice affected days we need a total count of the number of days that are ice affected.

If there are up to 3 days of continuous missing streamflow data, interpolate between the points.

Following Hecht et al. 2024 QA criteria, four conditions are screened per water year: (i) more than 3 consecutive days of missing data, (ii) more than 30 total missing days in the year, (iii) negative flow values, and (iv) constant monthly standard deviations during periods of non-zero streamflow (indicative of data errors or highly controlled discharge). Items i and ii REMOVE the year: any gap longer than 3 consecutive days, or more than 30 missing days in total, excludes that water year from all calculations. Internal gaps of up to 3 days are filled by linear interpolation for years that pass, but interpolated days still count toward the 30-missing-day limit (the count is taken on the raw data, before interpolation). Items iii and iv FLAG but do not remove: negative values are retained by default (config option reject_negative_flow, default false) and are counted in Negative_ann (see 3.12); a constant monthly standard deviation sets a QA flag only and never rejects a year.

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

### 3.1 Flow Volume

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

Function: calculate_negative_days()

Metrics:

Negative_ann: The number of days with negative flow values per gage per year. Negative-flow days are retained in the data by default (reject_negative_flow = false, Part 1.2); this metric accounts for them.

### 3.2 Flow Duration

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

### 3.3 Storage

Functions: analyze_baseflow_indices(), analyze_baseflow_indices_with_parameters()

Purpose: Calculates baseflow indices (BFI) using the Eckhardt and Lyne-Hollick digital filters, each in a default and a recession-parameterized variant.

Metrics:

BFI_Eckhardt: Baseflow index from the Eckhardt filter with default parameters. BFImax = 0.8, a = 0.98.

BFI_LyneHollick: Baseflow index from the Lyne-Hollick filter with default parameter a = 0.925.

BFI_Eckhardt_param: Eckhardt filter with the site-specific recession-derived constant in place of a = 0.98.

BFI_LyneHollick_param: Lyne-Hollick filter with the recession-derived constant (recession_alpha_point_cloud) in place of a = 0.925. This parameterization is heuristic — the L-H parameter has no formal derivation from recession analysis.

recession_alpha_point_cloud: Median of day-to-day recession ratios Q(i+1)/Q(i) across all recession events, under the linear-reservoir assumption (b = 1) for each gage. This constant parameterizes the two _param filters.

Method and parameters:

Filters are applied per water year, each year's time series (sorted by day of water year) is filtered independently, so the filter re-initializes every Oct 1 rather than running continuously over the record.

BFI = Σ(daily baseflow) / Σ(daily Q) within each water year, constrained to [0, 1]; a year with total Q ≤ 0 is skipped.

Eckhardt filter: BFImax = 0.8 (maximum baseflow index), a = 0.98 (recession constant). Single forward pass; initialized at baseflow₁ = BFImax · Q₁, and constrained to [0, Q].

Lyne-Hollick filter: a= 0.925, passes = 2 (each pass is a forward sweep followed by a backward sweep, so 2 passes = 4 sweeps total). Quickflow constrained to [0, Q] in both directions; baseflow = Q − quickflow, must be ≥0.

Parameterized versions: gages with insufficient recession data (10 or fewer alpha pairs in the whole record), or a derived constant outside [0, 1], produce NA for all parameterized BFI values.

Requirements and decisions:

Linear interpolation of streamflow for data gaps less than 3 days; years with data gaps > 3 days are rejected by the centralized preprocessor, so baseflow filters never encounter gaps in valid years. No filter restart or gap counting is needed.

Each day of the baseflow time series is equal to or less than the daily streamflow.

BFI must be between 0 and 1.

Function: analyze_recession_parameters()

Purpose: Analyzes recession events to calculate parameters like log(a), b, concavity, and recession seasonality from the relationship −dQ/dt = a·Q^b following Tashie et al. (2020a).

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

Recession events are identified as the longest contiguous windows where both Q and |dQ/dt| are monotonically decreasing, with a minimum length of 5 days (Tashie et al. 2020b). The first day of each event is removed during power law fitting (standard practice — storm peak influence).

Parameter b is determined using the line of best fit to the log-log plot of −dQ/dt versus Q (free power-law fit).

Parameter a is determined with the exponent fixed at b = 1, the linear reservoir assumption. The same fixed exponent is used at every gage and in every period. With b = 1, each day of recession gives one estimate of log(a), equal to log(-dQ/dt) - log(Q). The reported log(a) is the median of these daily estimates. The linear reservoir assumption is applied because in a free fit, log(a) is the intercept and b is the slope of the same line, and any change in b therefore also shifts (a). Fixing b = 1 removes this coupling, and therefore trends and seasonal cycles in log(a) reflect changes in log(a) alone, not changes in b.

Splits recession events temporally (first vs. second half) for concavity analysis.

Fits sinusoidal models to analyze seasonal variation in log(a) (fit to the event-scale b = 1 values).

All recession metrics except n_recession_events require at least 25 recession events across the record; gages below the minimum report NA.

The recession and parameterized-BFI families are exempt from the trend-completeness requirements and the 20-value statistics floor (Part 1.3).

Requirements and decisions:

Count the number of recession events in each year.

To control the quality of fitted parameters, calculate recession fits and create a flag for any R² < 0.8.

Function: analyze_Q_PPT_relationships(). Requires precipitation data

Purpose: Analyzes runoff ratios (streamflow divided by precipitation).

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

### 3.4 Flashiness

Function: calculate_pulse_metrics()

Purpose: Calculates metrics for high/low flow pulses (frequency, duration) and flow reversals.

Metrics:

n_high_pulses_year / n_high_pulses_all: Number of high-flow pulse events per year. A high-pulse day is a day with flow above the 90th percentile of flow; _year and _all calculate percentiles from Q for the specified year or the period of record, respectively. A pulse event is a maximal run of consecutive pulse days.

n_low_pulses_year / n_low_pulses_all: Number of low-flow pulse events per year. A low-pulse day is a day with flow below the 10th percentile of flow; _year and _all as above. Zero-flow days are included in all calculations; on strongly intermittent streams the 10th-percentile threshold can be exactly 0, in which case zero-flow days do not count as low-pulse days because the comparison is strictly below the threshold.

dur_high_pulses_year / dur_high_pulses_all: Mean duration in days of high-flow pulses per year.

dur_low_pulses_year / dur_low_pulses_all: Mean duration in days of low-flow pulses per year.

TQmean: Percentage of days with flow above the annual mean.

Flow_Reversals_annual / _winter / _spring / _summer / _fall: Number of flow reversals (annual and seasonal). Only includes changes in flow direction exceeding 2% of current flow: a reversal day requires a sign change between the incoming and outgoing day-to-day flow differences, with the outgoing change larger than 2% of that day's flow. Threshold of 2% selected to remove error / noise.

Requirements and decisions:

This metric requires data from across the period of record; confirm that the code can actually calculate the 90th/10th percentiles from the full period of record before the annual metrics.

High/low flow pulses are defined using the 90th and 10th percentiles of flow (Q).

Calculate annual metrics using period-of-record percentiles (_all); the _year variants are retained as a complementary per-year sensitivity.

Function: analyze_flashiness_trends()

Purpose: Calculates the Richards-Baker (R-B) flashiness index by water year from daily streamflow data.

Metrics:

flashinessRB (R-B Index): Sum of absolute differences in consecutive daily streamflow divided by the sum of total flow: R-B = Σ|Q(i) − Q(i−1)| / ΣQ(i) (Baker et al. 2004). Dimensionless; higher values indicate flashier streamflow response. I.e., higher day-to-day change in streamflow relative to total streamflow.

Function: calculate_streamflow_elasticity(). Requires precipitation data

Purpose: Calculates streamflow elasticity, which measures how sensitive streamflow is to changes in precipitation. From Sawicz et al. (2011). Values around 1.0 indicate proportional response; values greater than 1 define the catchment as elastic (sensitive to precipitation change), less than 1 as inelastic.

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

elasticity_years_total / elasticity_years_low_ppt: Diagnostic scalars: the number of years available for the gage, and the number excluded for annual PPT below 10 mm.

Requirements and decisions:

Requires daily streamflow (Q) and precipitation (PPT) data; calculating differences between years, so a minimum of two years of "good data" is needed for any pair, and a minimum of 15 years of valid data for the elasticity calculation overall.

Years with annual precipitation < 10 mm are flagged and excluded from analysis.

The year-over-year calculation requires both years of each pair to pass the minimum data requirements (consecutive qualifying years).

Year qualification is handled by the centralized preprocessor (config).

Count the number of years each site doesn't have sufficient data within a year, and what the number of excluded years would be if the threshold were < 30% data missing — added as documentation, not a filter.

Output: elasticity_static is a single value; elasticity_rolling and elasticity_annual each carry the 16 standard statistics (see Part 2). The elasticity family is exempt from the trend-completeness requirements and the 20-value statistics floor (Part 1.3).

### 3.5 Drought

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

### 3.6 Flow Timing

Function: analyze_flow_timing_trends()

Purpose: Analyzes the timing of flow events (day of maximum flow, cumulative flow percentiles). Day 1 = start of water year (Oct 1).

Metrics:

Dxx_day: Day of water year when cumulative flow reaches xx percentile. Labeled, for the 20th percentile flow for example, as D20_day. Calculated for percentiles = {1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99}.

D25_to_D75: Count of days between 25% and 75% cumulative flow.

Dmax: Water year day of maximum daily flow.

### 3.7 Precipitation Streamflow

Function: calculate_qp_seasonality(). Requires precipitation data

Purpose: Quantifies the seasonality in the relationship between cumulative streamflow (Q) and cumulative precipitation (P). Based on Wrede et al. (2015).

Metrics:

qp_slope_sd: Standard deviation of the 12 monthly mean cumulative Q-P slopes. Higher values indicate stronger seasonal variation in the streamflow-precipitation cumulative relationship.

qp_bimodality: Bimodality coefficient of the Q-P slope distribution = (skewness² + 1) / kurtosis. Values > 0.555 suggest bimodal or strongly seasonal patterns.

Method and parameters:

Requires daily streamflow (Q), precipitation (PPT), month, and day-of-water-year columns; year qualification is handled by the centralized preprocessor.

Calculates a 30-day rolling slope (backward-looking) of cumulative Q vs cumulative P for each water year, then aggregates the rolling slopes to 12 monthly mean values per year.

### 3.8 Snow

Function: calculate_snow_metrics(). Requires SWE data (Daymet)

Purpose: Fourteen metrics are computed from daily snow water equivalent (SWE; Daymet V4) for water years 1980–2023. All metrics operate on a thresholded series SWE* = SWE if SWE ≥ 10 mm, else SWE* = 0. Days with SWE below 10 mm are treated as snow-free for both durations and magnitudes. A "snow day" is a day with SWE* > 0. A "spell" is a maximal run of consecutive snow days and the "anchor spell" is the spell containing the annual SWE maximum. If the annual SWE maximum falls across more than one day, then the anchor spell will correspond to the period with the first day with the annual SWE maximum. Gages without Daymet coverage carry NA for all snow metrics.

Metrics:

swe_max: Maximum daily SWE in the water year (mm). A snow-free year reports a valid 0.

swe_max_dowy: Day of water year when maximum SWE (i.e., swe_max) occurs . NA in snow-free years.

snow_cover_days: Number of days with SWE ≥ 10 mm. A snow-free year reports a valid 0.

snow_on_dowy: First day of the anchor spell.

snow_off_dowy: First snow-free day after the anchor spell.

melt_season_days: Number of days during the melt season (i.e., between the day of max SWE and the first subsequent day of no snow). Calculated as: snow_off_dowy − swe_max_dowy.

melt_rate: rate of snowmelt during the melt season: swe_max divided by melt_season_days (mm/day; net ablation rate between peak and snow-off).

ssm: Snow seasonality metric = (seasonal snow days − ephemeral snow days) / total snow days, where a seasonal spell lasts at least 60 continuous days (Hatchett 2021; Petersky & Harpold 2018). Ranges −1 (fully ephemeral) to +1 (fully seasonal).

swe_apr1: SWE on April 1.

melt_before_peak: Total melt (sum of daily SWE decreases) occurring before the peak day (mm).

melt_before_peak_pct: melt_before_peak as a percentage of total water-year melt.

melt_before_peak_to_max_swe: Ratio of total melt to maximum SWE: melt_before_peak divided by swe_max.

melt_com_dowy: Day of water year when cumulative melt reaches 50% of the water-year total (melt center of mass).

swe_max_to_ppt: swe_max divided by total water-year precipitation (years with precipitation ≤ 10 mm return NA).

Requirements and decisions:

snow_on_dowy and snow_off_dowy are censored (NA) when the anchor spell touches Oct 1 or Sep 30, since the snowpack predates or outlasts the water year.

Timing and melt metrics are NA in snow-free years; magnitude metrics (swe_max, snow_cover_days, swe_apr1, swe_max_to_ppt) report valid zeros.

In addition to the standard trend-completeness requirements, the ten timing/melt/regime metrics require at least 80% of SWE-valid years in both the first and last decade of the gage's SWE record to be snowy before trend statistics are reported (prevents trends conditioned only on snow-present years at gages whose snow is disappearing or appearing).

## Part 4: Data Ingestion and Utility Functions

### The ingestion utility functions used to compile the harmonized input tables. They run once, upstream of the signature library.

### Function: process_gages_rawToRaw()

Purpose: Processes raw gage data (United States / Canadian) to calculate streamflow metrics for individual gages and save results to a file. Streamflow data are available up to the present day, but climate data must be processed independently.

Requirements and decisions:

A folder of USGS (United States) and HYDAT (Canadian) metadata must be stored locally in order to process watersheds.

A gage is compiled only if at least 20 calendar years each have at least 30 days with daily discharge above 0.0001 mm. This rule applies only at compilation.

All daily records in the retrieval window are kept once a gage qualifies; year rejection happens later in the preprocessor (Part 1.2).

### Function: process_caravan_gages()

Purpose: Processes Caravan NetCDF timeseries data for watersheds, calculates streamflow metrics, and saves results. Coincident climate data are generally available in Caravan, but Caravan’s record of streamflow is truncated to as early as the early 2000s (as late as 2018 for HYSETS).

Requirements and decisions:

Minimum years of valid data (min_num_years) are required for processing.

Filters years based on minimum flow thresholds and valid days (0.0001 mm / d for 30 days).

Handles redundancy between datasets (e.g., CAMELS and HYSETS) via common identifiers.

### Function: generate_streamflow_dt() and generate_streamflow_dt_caravan()

Purpose: Converts raw streamflow gage data (from both USGS/Canadian) or Caravan NetCDF data into a standardized streamflow data.table.

Requirements and decisions:

Streamflow is provided in appropriate units (or converted, e.g., m³/s to mm/day).

Filters data by specified date ranges and required years of data.

USGS values are kept only for qualifiers A, A e, P, P e (others set to NA); HYDAT values are never masked, their symbol is carried in the flag column.

Adds derived columns (e.g., year, month, day of year).

### Function: integrate_daymet_with_streamflow()

Purpose: Joins Daymet climate data (precipitation, temperature) to streamflow data based on date matching.

Requirements and decisions:

Requires a pre-processed Daymet parquet file (data_out/daymet_1980_2023.parquet).

Matches streamflow gage IDs to Daymet site IDs.

Joins climate variables (prcp, tmin, tmax, swe, vp, srad) to streamflow data.table on Date.

Renames prcp column to PPT for compatibility with existing Q-PPT analysis functions.

## Part 5: Automated Data Quality Flags

Function: compute_qa_flags()

Purpose: Adds boolean flagged_for_* columns to the per-gage summary table so that implausible or internally inconsistent results can be reviewed. A value of true marks a potential quality issue; flags never remove a gage from the downstream analysis (data are flagged only).

flagged_for_qann_range: Qann_mean outside [0, 2000] mm.

flagged_for_bfi_eckhardt_range: BFI_Eckhardt_mean outside [0, 1]. This is a confirmatory check; per-year BFI values are constrained to [0, 1] during calculation, so this flag should not appear in current code configuration.

flagged_for_bfi_lynehollick_range: BFI_LyneHollick_mean outside [0, 1]. Same confirmatory note as eckhardt bfi above.

flagged_for_flashiness_range: flashinessRB_mean outside [0, 2].

flagged_for_tqmean_range: TQmean_mean outside [0, 100].

flagged_for_d50_range: D50_day_mean outside [1, 366].

flagged_for_elasticity_range: elasticity_static outside [0.1, 5] (the per-gage scalar — the only flag not based on a _mean column).

flagged_for_runoff_ratio_range: annual_runoff_ratio_mean outside [0.01, 1.5].

flagged_for_seasonal_sum: the sum of the four seasonal totals (Qwin_mean + Qspr_mean + Qsum_mean + Qfal_mean) deviates from Qann_mean by more than 20%.

flagged_for_percentile_order: non-decreasing order Q5 ≤ Q25 ≤ Q50 ≤ Q75 ≤ Q95 violated (on _mean values; ties are allowed and do not flag).

flagged_for_timing_order: non-decreasing order D5 ≤ D50 ≤ D95 violated (on _day_mean values; ties are allowed and do not flag).

flagged_for_high_na: more than 30% of a gage's signature columns (the 16 statistics of every signature plus the 21 per-gage scalars) are NA — often indicates missing climate/SWE coverage rather than bad streamflow data.

Requirements and decisions:

A missing input value never triggers a flag.

Because the flags depend only on the summary table, they can be recomputed from an output CSV without rerunning any signature.

The same function, thresholds, and high-NA denominator are implemented in the Julia, Python, and R libraries (aligned 2026-09-04).

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

Tashie, A., Pavelsky, T., & Band, L. E. (2020)a. An empirical reevaluation of streamflow recession analysis at the continental scale. Water Resources Research, 56(1), e2019WR025448.

Tashie, A., Pavelsky, T., & Emanuel, R. E. (2020)b. Spatial and temporal patterns in baseflow recession in the continental United States. Water Resources Research, 56(3), e2019WR025425.

Wrede, S., Fenicia, F., Martinez-Carreras, N., Juilleret, J., Hissler, C., Krein, A., ... & Pfister, L. (2015). Towards more systematic perceptual model development: a case study using 3 Luxembourgish catchments. Hydrological Processes, 29(12), 2731–2750.

Yilmaz, K.K., Gupta, H.V., & Wagener, T. (2008). A process-based diagnostic approach to model evaluation: application to the NWS distributed hydrologic model. Water Resources Research, 44(9), W09417.

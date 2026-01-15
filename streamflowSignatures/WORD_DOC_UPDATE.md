# Content to Add to "Summary Documentation for Streamflow Signatures.docx"

**Instructions**: Copy the following sections into the Word document in the appropriate locations.

---

## Section 1: Summary of Functions (Add after `analyze_Q_PPT_relationships()`)

---

### calculate_streamflow_elasticity()

**Purpose**: Calculates streamflow elasticity, which measures how sensitive streamflow is to changes in precipitation. Based on Sawicz et al. (2011).

**Requirements and Decisions**:
- Requires daily streamflow (Q) and precipitation (PPT) data.
- Minimum of 15 years of valid data required for elasticity calculation.
- Years with annual precipitation <10mm are excluded.
- Calculates both static elasticity (median of all annual values) and rolling window elasticity (11-year window by default).
- Annual elasticity formula: E = (dQ/dP) / (Q_mean/P_mean), where dQ = Q_year - Q_mean and dP = P_year - P_mean.
- Values ~1.0 indicate proportional response; values >1 indicate amplified streamflow response to precipitation changes.
- Handles NA values by skipping years with >10% missing data.

**Output**: Returns 6 columns: `elasticity_static` (single value), plus 5 standard statistics for rolling elasticity (`elasticity_slp`, `elasticity_rho`, `elasticity_pval`, `elasticity_mean`, `elasticity_median`).

---

### calculate_qp_seasonality()

**Purpose**: Quantifies the seasonality in the relationship between cumulative streamflow (Q) and cumulative precipitation (P). Based on Wrede et al. (2015).

**Requirements and Decisions**:
- Requires daily streamflow (Q), precipitation (PPT), month, and day-of-water-year (dowy) columns.
- Minimum of 10 years of valid data required.
- Years with >10% missing Q or PPT data are excluded.
- Years with <300 valid days are excluded.
- Calculates 30-day rolling slope of cumulative Q vs cumulative P for each water year.
- Aggregates slopes to 12 monthly mean values per year.
- Two metrics are calculated per year:
  - **qp_slope_sd**: Standard deviation of the 12 monthly slopes (higher values = more seasonal variation).
  - **qp_bimodality**: Bimodality coefficient = (skewness² + 1) / kurtosis. Values >0.555 suggest bimodal/seasonal patterns.

**Output**: Returns 10 columns: 5 standard statistics for each metric (`qp_slope_sd_slp`, `qp_slope_sd_rho`, `qp_slope_sd_pval`, `qp_slope_sd_mean`, `qp_slope_sd_median`, `qp_bimodality_slp`, `qp_bimodality_rho`, `qp_bimodality_pval`, `qp_bimodality_mean`, `qp_bimodality_median`).

---

### calculate_average_storage()

**Purpose**: Estimates mean annual catchment storage using a simplified water balance approach. Based on Peters & Aulenbach (2011).

**Requirements and Decisions**:
- Requires daily streamflow (Q), precipitation (PPT), and day-of-water-year (dowy) columns.
- Minimum of 10 years of valid data required.
- Years with >10% missing Q or PPT data are excluded.
- Years with <300 valid days are excluded.
- Daily water balance: dS = P - Q (evapotranspiration is not estimated).
- Cumulative storage: S = cumsum(dS) for each water year.
- Annual average storage is interpolated from the storage-discharge relationship at mean annual discharge.
- Units: millimeters (mm).

**Output**: Returns 5 columns: `avg_storage_slp`, `avg_storage_rho`, `avg_storage_pval`, `avg_storage_mean`, `avg_storage_median`.

---

### integrate_daymet_with_streamflow()

**Purpose**: Joins Daymet climate data (precipitation, temperature) to streamflow data based on date matching.

**Requirements and Decisions**:
- Requires a pre-processed Daymet parquet file (`data_out/daymet_1980_2023.parquet`).
- Matches streamflow gage IDs to Daymet site IDs.
- Joins climate variables (PPT, tmin, tmax, swe, vp, srad) to streamflow data.table on Date.
- Reports coverage percentage; warns if <95% of dates have matching climate data.
- Renames `prcp` column to `PPT` for compatibility with existing Q-PPT analysis functions.

**Output**: Returns the input streamflow data.table with added climate columns (PPT, tmin, tmax, swe, vp, srad).

---

### convert_daymet_zip_to_parquet()

**Purpose**: One-time conversion of Daymet ZIP archive (containing 44 annual CSV files) to a single optimized parquet file.

**Requirements and Decisions**:
- Input: ZIP file containing daymet_1980.csv through daymet_2023.csv.
- Reconstructs the Date column from year, month, and row order (original CSVs lack a day column).
- Validates that each site has 365 or 366 days per year.
- Uses snappy compression for efficient storage (~70% reduction vs uncompressed).
- Output columns: site_id, Date, prcp, tmin, tmax, swe, vp, srad.

**Output**: Single parquet file (~3.9 GB) containing ~98 million rows for 6,087 sites.

---

## Section 2: Glossary of Streamflow Signatures (Add at end of glossary)

---

### Climate-Dependent Signatures (Require Precipitation Data)

**Runoff Ratio Metrics**
- annual_runoff_ratio: Annual streamflow divided by annual precipitation (Q/P).
- winter_runoff_ratio: Winter (Dec-Feb) streamflow divided by winter precipitation.
- spring_runoff_ratio: Spring (Mar-May) streamflow divided by spring precipitation.
- summer_runoff_ratio: Summer (Jun-Aug) streamflow divided by summer precipitation.
- fall_runoff_ratio: Fall (Sep-Nov) streamflow divided by fall precipitation.

**Streamflow Elasticity Metrics** (Sawicz et al. 2011)
- elasticity_static: Overall catchment elasticity - median of annual elasticity values. Measures how sensitively streamflow responds to precipitation changes. Values ~1.0 indicate proportional response; >1 indicates amplified response.
- elasticity: Rolling window (11-year) elasticity values, used to detect trends in catchment sensitivity over time.

**Q-P Seasonality Metrics** (Wrede et al. 2015)
- qp_slope_sd: Standard deviation of monthly cumulative Q-P slopes. Higher values indicate stronger seasonal variation in the streamflow-precipitation relationship.
- qp_bimodality: Bimodality coefficient of the Q-P slope distribution. Values >0.555 suggest bimodal or strongly seasonal patterns.

**Storage Metrics** (Peters & Aulenbach 2011)
- avg_storage: Mean annual catchment storage (mm) derived from water balance. Represents the average amount of water stored in the catchment at mean discharge conditions.

---

## Section 3: Glossary of Statistical Metrics (No changes needed)

The existing glossary already documents _slp, _rho, _pval, _mean, and _median.

---

## References

- Peters, N. E., & Aulenbach, B. T. (2011). Water storage at the Panola Mountain Research Watershed, Georgia, USA. Hydrological Processes, 25(25), 3878-3889.
- Sawicz, K., Wagener, T., Sivapalan, M., Troch, P. A., & Carrillo, G. (2011). Catchment classification: empirical analysis of hydrologic similarity based on catchment function in the eastern USA. Hydrology and Earth System Sciences, 15(9), 2895-2911.
- Wrede, S., Fenicia, F., Martinez-Carreras, N., Juilleret, J., Hissler, C., Krein, A., ... & Pfister, L. (2015). Towards more systematic perceptual model development: a case study using 3 Luxembourgish catchments. Hydrological Processes, 29(12), 2731-2750.

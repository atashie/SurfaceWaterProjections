# HISSS Manuscript Draft (auto-synced snapshot)

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/d/e/2PACX-1vS7j4FRp7SEwlXoBUVA8NA7cj_I0XzyS0u58r3bl8SOz4BfpZPrdPJge4RMcFocnX8Gnllkc1M-CTJ3/pub)
> **Last synced**: 2026-10-30 (227 paragraphs vs 106 on 2026-10-21; 15 regions). MAJOR
> co-author revision: Table 1 (data sources), Table 3 (land-cover datasets) and Table 4
> (resources) are now in the text as one-cell-per-paragraph tables; §2.1.1 regains the funnel
> paragraph (ephemeral filter, < 20 years, "85,000 km²" → 7,964) and the compile paragraph;
> §2.2.1 says 139 signatures and NINE families while the §2 preamble still says 121 / 100 + 21;
> §5 Technical Validation (three-level audit) is new and Usage Notes are §6.1–6.3; the
> funding paragraph is filled; references renumbered 1–43. Every number and method statement
> is under review: `docs/plans/2026-10-30-manuscript-accuracy-review-plan.md` (claim ledger
> §1); sync entry in `docs/reconciliation/manuscript_log.md` (2026-10-30).
>
> This is a read-only snapshot of the collaborative manuscript draft used for
> change detection and reconciliation review (the `/sync-docs` skill). Do NOT hand-edit the body below the header — it is overwritten at
> each sync. Discrepancies between the manuscript's methods and the code/docs
> are logged in `docs/reconciliation/manuscript_log.md`; manuscript-side corrections must be made in the Google Doc by the
> co-authors (@ Arik convention per the draft's own notes).

---




Journal: Scientific Data - Special Collection on Water Storage due Nov 9, 2026

Kaiser, K.E., A. Tashie, L. Lowman, K. Jennings, G. Gorski, D. Murray, Dymond, S.F.

FILL out Authorship contributions including Hydroshare user name; Update author list on hydroshare collection

Deleted Material

Active Signatures Doc

Abstract (modified from AGU)

The hydroclimate across North America is changing, with widespread transitions from snow-dominated to increasingly rain-dominated precipitation regimes. While such shifts in snowpack and winter hydroclimate are well documented, how these changes propagate streamflow behavior is neither uniform nor well characterized across climate, physiography, and critical zone structure, partly because of difficulties with standardizing and harmonizing disparate datasets. To meet this need, we present HISSS: Hydrologic Information Signatures and Summary Statistics - a new integrated dataset that includes annual snow, climate, and streamflow data for 8,014 North American watersheds, with trends and physiographic attributes for 6,041-7,964 locations depending on the subset. These products are intended for critical-zone scientists and beyond, who study watershed processes, snow hydrology, and vegetation dynamics, among others, and provide an extensive dataset for exploring changing hydroclimate and streamflow responses, and how these relationships play out across various physiographies, climate zones, and critical zone structures.

1 Background & Summary

The critical zone is the region of the Earth’s surface and subsurface where water is stored and cycled in the environment and is used by humans1. Water in the critical zone is stored and exchanged in and between soil, vegetation, snow, surface reservoirs, and groundwater. At the Earth’s surface, streamflow is an important metric of water storage and availability that indicates the severity of stress for human and natural systems2,3. The composition of the critical zone determines water availability for streamflow dynamics under normal conditions, disturbances (e.g., drought, fire, land conversion), and increasing atmospheric temperatures and CO2 4,5. Open questions in Earth system science include how system characteristics (i.e., soil texture and depth, vegetation type and rooting depth, hydroclimate, etc.) mediate water availability and surface and subsurface routing in the critical zone. Understanding the critical zone’s control on streamflow responses to climate is challenging because subsurface and aboveground critical zone components vary spatially and temporally1,2. Further, datasets describing critical zone features and climate are not readily available at spatial and temporal scales relevant for assessing impacts to streamflow. The Hydrologic Information Signatures and Summary Statistics (HISSS) dataset presented in this manuscript overcomes these issues by harmonizing climate, vegetation, soil, and human impact data at the watershed scale for >6,000 streamgages across the United States and Canada. HISSS provides a comprehensive streamflow signature dataset calculated from streamflow and gridded precipitation records spanning the years 1980 through 2025.

Many studies have considered how the critical zone mediates water exchange between surface and sub-surface storage components6. For instance, soil depth and texture influence the partitioning of rain and snow in the critical zone7, as well as runoff generation 8,9. The role of vegetation in mediating runoff responses has been explored previously for different drought scenarios10, with regard to the loss of global forests11 and land conversion in tropical forests12. Snow depth, snow water equivalent, and snowmelt timing influence runoff and groundwater recharge in the critical zone3. However, while prior studies explore relationships at individual watersheds or specific ecosystems, data are not available at continental scales that allow for cross-comparisons and syntheses across climate zones and ecoregions.

One challenge in harmonizing climate and vegetation datasets with streamflow is mismatches in temporal and spatial scales13. Streamflow represents point-scale measurements that integrate water fluxes across a watershed, while climate and vegetation data are typically provided as gridded products representing average values across a defined pixel at a point in time. Several large-sample datasets have addressed this mismatch by aggregating gridded and geospatial data to the watershed scale, with products focused on hydroclimate and hydrological signatures14,15, biogeochemistry13, stream chemistry16, as well as broad multi-source archives that now cover much of North America17,18. Among the products that report hydrological signatures, data are given only as single summary values rather than as a time series. Thus, these products do not provide per-signature trend or changepoint statistics, and the harmonized multi-region archives typically derive their climate forcing from comparatively coarse reanalyses (e.g., ERA5-Land at ~11 km or coarser). The recent CAMELS-SPAT19 dataset spans the United States and Canada and includes high spatial resolution climate data, time-varying MODIS leaf area index (LAI), and land use and land cover (LULC). However, it is focused on a subset of available watersheds (n < 1500) to support hydrological model development and therefore precomputes only a narrow subset of hydroclimatological signatures as period-of-record statistics.

The HISSS dataset and software library presented here provide holistic hydroclimate, vegetation, snow, and geophysical information at the scale of individual watersheds across the continental United States and Canada. For over 6,000 watersheds in the study, we calculate roughly 100 hydrological signatures and report the underlying annual values along with a consistent suite of computed metrics: central tendency statistics (mean, median), linear and robust (Theil-Sen) trend slopes, Spearman and Mann-Kendall trend significance, Pettitt changepoint detection, and other supporting diagnostics. These are paired with annually resolved Moderate Resolution Imaging Spectroradiometer (MODIS) and National Land Cover Dataset (NLCD) LULC and monthly resolved MODIS LAI, watershed climate features derived from ~1 km Daymet, and HydroATLAS physiographic attributes. Because every signature is computed by an open, cross-language library (Julia, Python, and R), each metric is reproducible and extensible. These products are intended for critical-zone scientists and beyond, who study watershed processes, snow hydrology, and vegetation dynamics, among others.

2 Methods

To build the HISSS dataset, we accessed numerous time-series and geospatial datasets across the United States (US) and Canada (Table 1) using an open-source, reproducible processing workflow (Figure 1). We assembled daily mean streamflow records for 16,994 gages across the United States (9,154; USGS National Water Information System, NWIS) 20and Canada (7,840; Hydrometric Historical Data, HYDAT)21, of which 8,014 (6,160 USGS; 1,854 HYDAT) yielded usable daily records for watersheds smaller than 85,000 km2, totalling more than 111 million observations. We defined the watershed boundary for each gage using official agency basin polygons and the HydroBASINS dataset (Sect. 2.1.1). These polygons defined the basin boundaries over which we aggregated gridded climate, land-cover, and basin-attribute data (Sects. 2.1.3–2.1.4). For each basin and its suite of hydroclimatic and geospatial data, we computed 121 signatures across nine signature families (drought, flow volume, flashiness, precipitation-streamflow, flow duration, flow timing, snow, storage, climate): 100 annually resolved signatures, each summarized with a common set of 16 statistics, plus 21 per-gage diagnostics. Trends on these signatures and statistics were calculated for 1980-2025 (6,250 sites) and for 1993-2025 (6,678 sites) (Sect. 2.1.2). The following subsections detail the input data, filtering procedures, metric calculations, and trend analyses. The code for the full workflow is available in an open library in Julia, Python, and R (https://github.com/CZ-Sync/HISSS).

Figure 1. Workflow used to generate the Hydrologic Information Signatures and Summary Statistics (HISSS) dataset.

Table 1. Summary of data sources used to generate the Hydrologic Information Signatures and Summary Statistics (HISSS) and the number of sites, temporal resolution, and temporal coverage included. <TABLE OF ALL INPUT DATA SOURCES>

Dataset type

Spatial Unit

Temporal Resolution

Data Sources

Number of Spatial Units in Dataset

Temporal Coverage

Streamflow Timeseries

Agency gage ID

Daily

NWIS20

6,160

1980-2025

HYDAT21

1,854

Cathment Polygons

HydroBASINS22; USGS23 & Canada24

8,014

Static Basin Characteristics

Catchment Polygons

-

HydroAtlas25

8,014

-

Caculated Derivatives of Streamflow, Climate & Landcover

Annual Signatures & Summary Statistics

Gage IDs assigned to polygons

Annual & Multi-decadal

Daymet21,26, MODIS27,28, NLCD29,30, NWIS, HYDAT

6,041 - 7,964*

variable

Trends**

6,250

1980-2025

6,678

1993-2025

*Number of sites with annual signatures and summary statistics is variable across datasets. Landcover datasets include MODIS (7,964; 2021 -2024) and NLCD (6,119, 1985-2024); Climatic data from DAYMET (6,041)

**trends require at least X% of data in the first and last decade, the shorter timeseries is not necessarily a subset of the longer, as the data sufficiency of an individual site will be different in each time period

2.1 Input data

2.1.1 Streamflow data

We retrieved daily mean observed streamflow from 01 October 1979 to 30 September 2025 for all candidate gages (n = 16,994) using the USGS’s dataRetrieval R package31 for US gages and the Water Survey of Canada’s tidyhydat R package32 for Canadian gages. To ensure we retained only the most useful observations for time series analysis, we applied standardized quality assurance criteria at three stages during the calculation of annualized data: 1) retrieval, 2) compilation, and 3) per-year qualification.

From the candidate gages, we applied a filter to remove ephemeral streams (perennial defined as >30 days of streamflow data over 0.0001 mm/day)23, and locations with less than 20 years of streamflow data. We also subset the gages to those located in watersheds with polygons smaller than 85,000 km2 (≤ largest HUC8 size). This resulted in 7,964 gages across the United States and Canada.

These data were compiled into a single dataset of 111.6 million observations keyed to the respective agency’s gage ID. We retained original agency quality flags, including ice-affected and estimated values. We stored missing days as NA rather than zero-filled, treated zero flow as a valid observation, and area-normalized discharge (mm d-1) using agency-published drainage areas for all but 73 gages, for which no drainage area is published (see Usage Notes). We aggregated each gage’s record to a continuous daily time series and applied the following data quality screening by water year (WY; October 01 - September 30): gaps of up to three consecutive days were filled by linear interpolation, any WY with more than 30 missing values were set to NA, and if any data gap exceeded three consecutive days those days were set to NA. Seasonal periods are defined as: Winter: January-March, Spring: April-June, Summer: July - September, Fall: October-December. Annual metrics were calculated on the water year October 1 - September 30.

2.1.2 Watershed Boundaries

To define the spatial extents of the contributing area to the gaged locations, we combined official agency basin polygons from the US and Canada. For the US, we used watershed polygons from the USGS GAGES-II dataset33 and merged the individual shapefiles for reference and non-reference gages across all US regions into a single US-wide shapefile (covering 100% of US gages). For Canada, we used the Water Survey of Canada’s Hydrometric Network Basin Polygons dataset24 and merged them into a single Canada-wide layer (covering 98.3% of Canadian gages). Gages lacking an official polygon were delineated from the HydroBASINS level-12 drainage network34 (detailed methods below) by aggregating all basins upstream of the gage’s outlet. We leveraged HydroATLAS25 to define static basin attributes within the data product domain at the gages with suitable streamflow records. The source of the watershed geometry was recorded for each watershed. We excluded basins that exceeded 85,000 km2, transformed the remaining watershed boundaries into a common projection, and then merged them into one full-domain layer. The final watershed boundaries were used to aggregate climate and landcover data only where we had sufficient streamflow data coverage for trends analyses (see Section 2.1.3).

2.1.3 Climate and snow data

Meteorological data came from Daymet26, a daily long-term 1 km2 gridded continuous product covering North America. Daily Daymet estimates of minimum and maximum air temperature, precipitation, vapor pressure, shortwave radiation, snow water equivalent, and day length covering 1980–2023 were downloaded from the ORNL Distributed Active Archive Center (DAAC)26. The watershed size threshold corresponds to the largest HUC8 basin in the Watershed Boundary Dataset35, which was used to ensure that hydrologic response to weather was meaningfully well represented with basin-averaged climate inputs. We performed this aggregation using the gdptools package36, which performs area-weighted zonal aggregation of gridded data to polygon features across all variables and time steps. The result is a single area-weighted daily value for each variable for each basin.

2.1.4 Basin attributes and land cover

The HISSS dataset includes both static physiographic characteristics and landcover characterization through time (Table 3). We used static basin attributes from HydroATLAS (BasinATLAS version 10)25, a global compilation of hydro-environmental descriptors from HydroBASINS level-12 subbasins34. HydroBASINS is a global, hierarchically nested set of subbasin polygons derived from a 15 arcsecond drainage network. Level 12 is its finest tier (mean area of 110 km2), and each unit carries a downstream neighbor identifier so that a gage’s full contributing area can be assembled by aggregating the upstream network. For every gage, we identified the intersecting level-12 subbasin and the complete set of upstream subbasins, then aggregated the HydroATLAS attributes over that contributing area. Continuous and percentage attributes were aggregated as area-weighted means, elevation as the spatial minimum and maximum, and categorical attributes as the area-weighted majority. The resulting table is keyed to gage ID so that it joins directly to streamflow signatures. For headwater basins smaller than the intersecting level-12 HydroATLAS polygons, the static basin attributes include downstream components of the larger watershed, in addition to the headwater basin of interest. In total, we derived approximately 210 attributes per watershed describing hydrology, physiography, climate, LULC, soils and geology, and anthropogenic influence.

Table 3. Summary of landcover datasets included in HISSS {tab in the tables xls}

Data Source

Characteristics

Details

Timespan

Temporal Resolution

Spatial Coverage

Spatial Resolution

Hydro

ATLAS22

210 attributes per watershed

6 categories: hydrology, physiography, climate, LULC, soils and geology, and anthropogenic influence.

static

Global

Basin Polygon

MODIS MCD12Q128

Land Cover Type

5 Land Cover Type schemes: 1-International Geosphere-Biosphere Programme (IGBP), 2-University of Maryland (UMD), 3- Leaf Area Index (LAI) classification, 4-BIOME-Biogeochemical Cycles (BGC) classification, and 5-Plant Functional Types (PFT) describing land cover, land use, and surface hydrology

Benefits: Provides land cover classifications for CONUS and Canada using 5 different classification schemes that are relevant for different applications.

Limitations: Data record begins in 2021.

2021-2024

Annual

Global

500m

MODIS MCD15A3H27

Leaf Area Index (LAI)

Leaf area index (LAI), a measure of canopy density and indicator of plant phenology

Benefits: Provides vegetation state information for CONUS and Canada.

Limitations: Data record begins in 2022.

2022-2024

Monthly

Global

500m

USGS National Land Cover Database (NLCD)29,30

Land Cover Type

16 land cover classes together with basin-mean fractional impervious surface, a direct measure of urbanization intensity that the MODIS products do not offer

Benefits: Extends the land cover record 16 years further back and resolves land cover heterogeneity within small watersheds that 500 m MODIS pixels may miss.

Limitations: Spatial coverage limited to CONUS.

1985-2025

Annual

CONUS

30m

To obtain information on how land surface characteristics change through time, the HISSS dataset also includes land cover and vegetation characteristics from MODIS and NLCD (Table 3). This includes all land cover classification schemes from MODIS MCD12Q1 Land Cover Type Yearly product and LAI from the MODIS MCD15A3H, both v06127. The MODIS LAI product represents the seasonal phenology of vegetation on the land surface. It is estimated from a look-up table derived from a 3-D radiative transfer model relating surface reflectance to observed LAI values37. It also relies on the MODIS Land Cover Type 3 (LAI Classification) to determine the most likely LAI value for a pixel given its vegetation type. MODIS LAI has valid values between 0 and 10 m2 m-2, which describe the one-sided projection of leaf area per unit ground area.

To complement the MODIS land cover record, HISSS also includes land cover from the USGS Annual National Land Cover Database (Annual NLCD, Collection 1), a Landsat-derived product mapped at 30 m resolution annually from 1985 through 202529,30. For each gaged watershed we provide annual percent coverage of the 16 land cover classes together with basin-mean fractional impervious surface, a direct measure of urbanization intensity that the MODIS products do not offer. Relative to MCD12Q1, Annual NLCD extends the land cover record 16 years further back and resolves land cover heterogeneity within small watersheds that 500 m MODIS pixels may miss. Its principal limitation is its spatial coverage: Annual NLCD spans CONUS only, so it is provided for only those 6,119 watersheds while Alaskan and Canadian watersheds retain MODIS-only land cover. The two products are intended to be complementary, with MCD12Q1 supplying consistent coverage across all HISSS watersheds, while Annual NLCD supplies a longer, higher-resolution record for CONUS watersheds.

2.2 Metrics for analysis

2.2.1 Computing Streamflow and Hydroclimate Signatures

We calculated a comprehensive suite of 139 streamflow and hydroclimate signatures for each gage and its associated watershed, including signature-specific diagnostics and data-quality flags. The selection of metrics is consistent with analysis in catchment hydrology38 and snow hydrology39,40 and span the nine signature families: flow volume, flow duration, storage, flashiness, drought, flow timing, precipitation-to-streamflow, snow, and climate. Details for each metric, including the definition, requirements, and relevant citations, are in the data dictionary (hisss_data_dictionary.csv) included with Resources 1 and 2.

We calculated signatures for gages only if they retained at least 20 qualifying WYs and spanned at least 60% of the possible WYs in the analysis window (additional per-signature completeness requirements determine whether trend statistics are reported for a given signature, as described in Section 2.2.2). Water year and seasonal metrics required having at least 80% of raw observations passing quality assurance, if this was not met, the affected metrics were set to NA for that WY or season. We calculated signatures (Sect. 2.2) over two standard analysis windows (WYs 1993-2025 and 1980-2025), yielding 6,678 gages for WYs 1993–2025 and 6,250 for WYs 1980–2025. The shorter time range (1993-2025) was included as it is the date range with the most gages with the same data record length, allowing for direct comparison of signatures across a standardized date range.

2.2.2 Computing Trends and Summary Statistics

For the trends analysis, we removed years with more than 3 consecutive days of NAs or more than 30 days of NAs throughout the water year41. We flagged (but did not remove) water years that contained negative streamflow values or any calendar month in which at least 15 days of non-zero streamflow held at a single constant value, patterns that are indicative of data errors or highly regulated discharge. For a given signature, trend statistics were computed only if its annual series had at least 20 valid values, at least 60% of the entire series had valid values, and at least 80% of both the first and last decades of the trend had valid values; the event-based recession and elasticity metrics are exempt from these completeness requirements.

3 Data Record

The HISSS dataset is published in HydroShare, the repository for water data operated by the Consortium of Universities for the Advancement of Hydrological Science, Inc. (CUAHSI; https://www.hydroshare.org). The dataset is organized as a public collection (DOI: http://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768) grouping five resources, each independently citable with its own DOI (Table 4). All resources are released under a CC-BY 4.0 license, and each contains a README describing its files and a machine readable data dictionary defining every column. Tabular data are provided as UTF-8 CSV where file sizes permit and as Apache Parquet for large tables; Parquet files are readable with standard open-source libraries in R (arrow), Python (pyarrow/pandas), and Julia. The processing code is archived separately (see Code Availability).

Table 4. Summary of resources within the HISSS dataset. <hisss_data_resources_table>

Resource

DOI

Number of Gages

Contents

Formats

1. Streamflow signatures, WY 1993–2025

Generated upon acceptance

6,678

Signature summary table (6,678 gages), per-signature annual values, interactive explorer, run provenance and validation reports

CSV, Parquet, HTML, JSON

2. Streamflow signatures, WY 1980–2025

Generated upon acceptance

6,250

As Resource 1, for the full-record window

CSV, Parquet, HTML, JSON

3. Harmonized hydroclimate data

Generated upon acceptance

8,014

Daily streamflow, gage metadata (metadata included for all 16,994 candidate gages)

Parquet, CSV

6,087

Daily basin-averaged Daymet climate

4. Watershed geometry and basin attributes

Generated upon acceptance

8,014

HydroATLAS-derived watershed attributes

GeoPackage, Parquet, CSV

7,964

Watershed boundary polygons

5. Vegetation and land cover

Generated upon acceptance

7,964

Monthly MODIS LAI, annual MODIS land cover, interactive explorers

Parquet, CSV, HTML

6,119

Annual NLCD land cover and impervious surface (CONUS), interactive explorer

Resources 1 and 2 (Hydroclimate streamflow signatures). These parallel resources implement the two analysis windows described in Section 2.1.2 (water years 1993-2025 and 1980-2025); gages qualify independently under each window’s completeness criteria, so neither product is a subset of the other. Each resource contains: (i) a signature summary table (hiss_signatures_wy{window}.csv) with one row per gage and 1,653 columns comprising gage metadata, the 100 signatures x 16 statistics each, 21 signatures and related outputs that do not carry statistics, and 12 automated quality-assurance flags; (ii) a long-format annual-values table (Parquet) holding the yearly value behind every signature-statistic pair, with columns gage_id, signature, water_year, and value (18,898,406 and 24,366,487 rows for the two windows, respectively); (iii) a self-contained interactive HTML explorer that maps every signature-statistic combination and plots per-gage annual series with fitted checksums, configuration, and software versions) and validation reports.

Resource 3 (harmonized input data). Daily discharge for all 8,014 processed gages (111,624,189 records; 1980 through 2025) compiled from the USGS National Water Information System via dataRetrieval31 and the Water Survey of Canada HYDAT database via tidyhydat21,32, expressed in mm d-1 (see area_normalized below). The accompanying gage metadata table (16,994 candidate gages) records station coordinates, published drainage areas, processing status, and human-interference indicators drawn from GAGES-II attributes23,33 and Reference Hydrometric Basin Network (RHBN) and regulation status via HYDAT. Daily Basin-averaged Daymet Version 4 R126 for the variables precipitation, temperature, snow water equivalent (SWE), vapor pressure, and shortwave radiation is provided for 6,087 basins (97,757,220 records, calendar years 1980-2023) as produced by the area-weighted aggregation described in Section 2.1.3.

Resource 4 (watershed geometry and basin attributes). Watershed boundary polygons for 7,964 gaged watersheds (GeoPackage and GeoParquet, WGS84), assembled from GAGES-II basin boundaries23,33, the Water Survey of Canada hydrometric basin polygons24, and HydroBASINS-derived delineations where no official polygon exists34; the provenance of each polygon is recorded in watershed_geom_source. A companion table provides 211 watershed-scale attributes (climate, hydrology, terrain, land cover, soils, and geology, and anthropogenic influence) aggregated over each gage’s full upstream area from BasinAtlas v1022,25, with a data dictionary mapping each column to its theme, unit, and aggregation method.

Resource 5 (vegetation and land cover). Three per-watershed time-series tables keyed to watershed polygons (Resource 4): (i) monthly MODIS LAI (MCD15A3H v061; Myneni et al. 2015) for 2002-2024 (2,150,280 rows; 7,964 watersheds x 270 months) with basin-mean monthly LAI, spatial-distribution statistics (quantiles, standard deviation, extrema), and per-month quality fractions; (ii) annual MODIS land cover (MCD12Q1 v061)28 for 2001-2024 (191,136 rows with 102 per-class percent-coverage columns across all eight classification schemes, as described in Section 2.1.4); and (iii) annual NLCD land cover and fractional impervious surface (Annual NLCD Collection 1)30 for 1985-2025 (250,879 rows covering 6,119 CONUS watersheds, with 16 per-class percent-coverage columns and basin-mean imperviousness. Each table is accompanied by its data dictionary, granule-level provenance manifests for the MODIS products, and a self-contained interactive HTML explorer.

Conventions and key fields. All tables join on gage_id, the agency station identifier (USGS site numbers; HYDAT station numbers); read it as a character string and join with leading zeroes removed; the geometry and land-cover tables also carry this form as canon_id. Files derived from the geometry layer also carry canon_id, a zero-stripped variant (for USGS sites) retained for joins to legacy metadata. water_year denotes the period 1 October - 30 September, labeled by the ending calendar year. Signature columns in Resources 1-2 follow the pattern {signature}_{statistic}, where the sixteen statistic suffixes are: _senn_slp (Theil-Sen slope), _linear_slp (ordinary least-squares slope), _spearman_rho / _spearman_pval (Spearman rank correlation with time and its p-value), _mk_rho / _mk_pval (Mann-Kendall tau and p-value), _mean, _median, and eight Pettitt changepoint fields beginning with _pettitt including _cp_year (most likely changepoint year), _pval, _pre_mean / _post_mean / _delta_mean / _pct_change (segment means and their differences), and _pre_mk_pval / _post_mk_pval (within-segment trend tests). A small number of signatures are single-valued per gage s(e.g., elasticity_static, the recession-derived filter constant, and diagnostic counts) and appear without the suffix set; these are enumerated in the data dictionary. In the annual-values tables, a missing row and a stored NaN are equivalent, both meaning “not computable for that gage-year.” The Boolean area-normalized flag (see unnormalized gages section 5.1.1) identifies whether a gage’s discharge is expressed in mm d-1 (true) or retained in native m3 s-1 because no published drainage area exists (false). Columns prefixed flagged_for_ are automated range- and consistency-check indicators defined in the dictionary. Land cover columns follow {scheme}_c{code}_pct (MODIS) and nlcd_c{code}_pct} (NLCD), with pct indicating the percent of the watershed area that is in the class code.

4 Data Overview

Sites included in the HISSS dataset span a broad geographic range across the United States and Canada, covering 12 Level I North American Ecoregions42. A majority (62%) of the sites are located in the Eastern Temperate Forests or the Northwestern Forested Mountains, while there are fewer than 10 sites within each of the Southern Semiarid Highlands, Tundra, and the Hudson Plain ecoregions (Figure 2). HISSS sites comprise small basins encompassing headwater streams up to large river basins with a median watershed size of 698 km2, 25% of sites having drainage areas greater than 2,700 km2, and 25% of sites having drainage areas less than 180 km2. Streamflow signature data coverage was lowest for the family of storage metrics which include baseflow indices which have more specific data requirements (Figure 3A). The majority of sites have greater than 80% streamflow signature coverage across all families of signatures, with the lowest coverage at sites within the ecoregions taiga and Husdon plain (Figure 3B).

Figure 2. Distribution of sites with delineated catchment areas and paired streamflow and meteorological data for the 1993-2025 period; sites are color-coded by their Level I North American Ecoregion42 with the number of sites falling within each ecoregion shown in parentheses in the map legend. Inset shows distribution of basin size in km2.

Figure 3. Boxplots showing distributions of coverage among sites. For both panels, coverage is calculated as the number of metric-years/total possible metric-years for each site, and darker boxes show distributions for the time period 1980-2025, while lighter boxes show distributions for 1993-2025. a) Coverage by signature family, where the number below each signature family is the number of signatures within that family. b) Coverage by ecoregions, beneath each ecoregion label, the number of sites within that ecoregion is shown for 1980-2025 (first number) and 1993-2025 (second number), which add up to 6,250 and 6,678, respectively.

5 Technical Validation

To generate the HISSS data resource, we aggregated existing, freely-available datasets that have previously undergone thorough validation independently (Table 1). That is, the aggregation of the climate and streamflow data into basin polygons and calculation of summary statistics and signatures were all completed using previously validated data. All code developed to conduct the data aggregation and calculate summary statistics and signatures underwent manual code review by the authors. For thoroughness, we implemented a three-level technical audit to validate a subset of streamflow and climate signatures, as it was not feasible to validate all 139 signatures for every watershed. This audit involved 1) a visual inspection of the spatial distribution of signatures across the US and Canada by ecoregion and climate, 2) a spot check of individual gauges by the author team to ensure that streamflow and climate signatures had expected trends and magnitudes, and 3) qualitative evaluations of computed streamflow and climate signatures for individual watersheds against values published the literature.

6 Usage Notes

6.1 Unnormalized gages

73 of the 8,014 processed gages have no agency-published drainage area (32 and 28 of them qualify for the WY 1993–2025 and WY 1980–2025 products, respectively). Most are irrigation or diversion canals, dam and powerhouse outflows, or channel splits of large rivers, where a contributing area is undefined or unpublished. These gages are retained with discharge in native m3 s-1 and flagged area_normalized = FALSE. Their unit-carrying signatures (flow volumes, flow percentiles, recession log(a), drought deficits) are not comparable with the mm d-1 values at other gages, and their precipitation-dependent signatures (runoff ratios, elasticity, Q-P seasonality, storage) are structurally NA. Users should filter on area_normalized == TRUE before any cross-gage comparison of unit-carrying signatures.

6.2 How to join HISSS with other datasets

Every HISSS table is keyed by the agency station identifier named “gage_id” (USGS site number, e.g. 0103500; HYDAT station number, e.g. 01AD002). This station identifier or gage_id may be read as a character string and joined on the leading-zero-stripped form for each dataset (most US site numbers begin with a zero that numeric parsing silently ignores (see resource README for details)). Any dataset keyed to the same agency identifiers then joins directly. Datasets that are not gage-keyed join spatially through the watershed polygons of Resource 4: gridded products such as the Daymet-derived VPD of Corak et al. (2025)43, which shares Daymet's 1 km grid and 365-day calendar, can be aggregated with the area-weighted workflow used for Daymet (Sect. 2.1.3), and point records such as MacroSheds sites13 can be located within the nested polygons or matched to the nearest gage by coordinates. HydroATLAS-based products join on the HydroBASINS level-12 outlet identifier (Downstream_HB_ID) in the Resource 4 attribute table. When combining data, aggregate partner records to the water year (October–September, labeled by the ending year); re-normalize fluxes to a common drainage area, since HISSS reports both the agency-published area used for normalization (basin_area) and the polygon area (geom_area_km2); and pair record-dependent signatures (Sect. 5.1.3) with the product whose window matches the partner record.

6.3 Record-dependent signatures

Hydroclimate signatures whose definition uses thresholds or means from the full analysis window are valid within the product for that time window and should not be compared against the WY 1993–2025 and WY 1980–2025 products, nor re-derived from the annual values over a different window. The specific signatures include the period-of-record (*_all) pulse metrics, elasticity, the recession-parameterized baseflow indices, and all drought metrics.

Data availability

The HISSS dataset is openly published in HydroShare at http://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768 licensed under a Creative Commons Attribution 4.0 International (CC BY 4.0) licenceCC-BY 4.0 license.

Code availability

The code for the full workflow is available in an open library in Julia, Python, and R, https://github.com/CZ-Sync/HISSS. Details for each metric, including the definition, requirements, and relevant citations, are in the Streamflow Signatures Reference on GitHub (https://github.com/CZ-Sync/HISSS/blob/main/docs/SIGNATURES.md). And data resources and their associated metadata and guidance for usage are in https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/ [TEMPORARY: WILL CHANGE WHEN WE PUBLISH].

Acknowledgements

AI-assisted coding tools (Claude Code 0.145.0, Anthropic) were employed to generate code used in data ingestion and processing, hydrological signature extraction, cross-language validation and benchmarking, and interactive visualization. All code was reviewed, tested, and validated by the authors to ensure correctness and reproducibility. Generative AI tools (Claude Code 0.145.0, Anthropic) were used to support data analysis and visualization. These tools were applied under the supervision of the authors, and all outputs were reviewed and validated against established scientific methods to ensure reproducibility and transparency.

The authors sincerely thank Maggie Phillips and Steven Schmitz for assisting with the validation of the trends and signatures. Additionally, the authors would like to thank Mike Wieczorek of the USGS for providing early access to Daymet data for initial prototyping of the workflow.

Funding

This work was conducted as a part of the “Critical zone as a mediator of hydroclimate-ecosystem asynchrony” Working Group supported by the John Wesley Powell Center for Analysis and Synthesis, funded by the U.S. Geological Survey [and other funding sources if appropriate].

Support for DSM was provided by the United States National Science Foundation (NSF) 626 Earth Sciences Postdoctoral Fellowship (9500317412). Support for KEK was provided through Idaho State Board of Education appropriations Award Number CW8836.

Disclaimers

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

References

1. Lin, H. Earth’s Critical Zone and hydropedology: concepts, characteristics, and advances. Hydrol. Earth Syst. Sci. 14, 25–45 (2010).

2. Brooks, P. D. et al. Hydrological partitioning in the critical zone: Recent advances and opportunities for developing transferable understanding of water cycle dynamics. Water Resour. Res. 51, 6973–6987 (2015).

3. Wlostowski, A. N. et al. Signatures of Hydrologic Function Across the Critical Zone Observatory Network. Water Resour. Res. 57, e2019WR026635 (2021).

4. Keller, C. K. Carbon Exports from Terrestrial Ecosystems: A Critical-Zone Framework. Ecosystems 22, 1691–1705 (2019).

5. Condon, L. E. et al. Where Is the Bottom of a Watershed? Water Resour. Res. 56, e2019WR026010 (2020).

6. Tague, C. et al. James Buttle Review: Dynamic Water Storage Shapes Critical Zone Function in Snow-Dominated Mountain Watersheds. Hydrol. Process. 39, e70325 (2025).

7. Hammond, J. C., Harpold, A. A., Weiss, S. & Kampf, S. K. Partitioning snowmelt and rainfall in the critical zone: effects of climate type and soil properties. Hydrol. Earth Syst. Sci. 23, 3553–3570 (2019).

8. Zimmer, M. A. & Gannon, J. P. Run-off processes from mountains to foothills: The role of soil stratigraphy and structure in influencing run-off characteristics across high to low relief landscapes. Hydrol. Process. 32, 1546–1560 (2018).

9. Chen, H., Niu, Q., McNamara, J. P. & Flores, A. N. Influence of Subsurface Critical Zone Structure on Hydrological Partitioning in Mountainous Headwater Catchments. Geophys. Res. Lett. 51, e2023GL106964 (2024).

10. Corak, N. K., Otkin, J. A., Ford, T. W. & Lowman, L. E. L. Unraveling phenological and stomatal responses to flash drought and implications for water and carbon budgets. Hydrol. Earth Syst. Sci. 28, 1827–1851 (2024).

11. Zhang, M. et al. A global review on hydrological responses to forest change across multiple spatial scales: Importance of scale, climate, forest type and hydrological regime. J. Hydrol. 546, 44–59 (2017).

12. Muñoz-Villers, L. E. & McDonnell, J. J. Land use change effects on runoff generation in a humid tropical montane cloud forest region. Hydrol. Earth Syst. Sci. 17, 3543–3560 (2013).

13. Vlah, M. J. et al. MacroSheds: A synthesis of long-term biogeochemical, hydroclimatic, and geospatial data from small watershed ecosystem studies. Limnol. Oceanogr. Lett. 8, 419–452 (2023).

14. Newman, A. J. et al. Development of a large-sample watershed-scale hydrometeorological data set for the contiguous USA: data set characteristics and assessment of regional variability in hydrologic model performance. Hydrol. Earth Syst. Sci. 19, 209–223 (2015).

15. Addor, N., Newman, A. J., Mizukami, N. & Clark, M. P. The CAMELS data set: catchment attributes and meteorology for large-sample studies. Hydrol. Earth Syst. Sci. 21, 5293–5313 (2017).

16. Sterle, G. et al. CAMELS-Chem: augmenting CAMELS (Catchment Attributes and Meteorology for Large-sample Studies) with atmospheric and stream water chemistry data. Hydrol. Earth Syst. Sci. 28, 611–630 (2024).

17. Arsenault, R. et al. A comprehensive, multisource database for hydrometeorological modeling of 14,425 North American watersheds. Sci. Data 7, 243 (2020).

18. Kratzert, F. et al. Caravan - A global community dataset for large-sample hydrology. Sci. Data 10, 61 (2023).

19. Knoben, W. J. M. et al. Catchment Attributes and MEteorology for Large-Sample SPATially distributed analysis (CAMELS-SPAT): streamflow observations, forcing data and geospatial data for hydrologic studies across North America. Hydrol. Earth Syst. Sci. 29, 5791–5833 (2025).

20. U.S. Geological Survey. U.S. Geological Survey National Water Information System database. https://doi.org/https://doi.org/10.5066/F7P55KJN (2024).

21. Environment and Climate Change Canada. HYDAT: National Water Data Archive (hydrometric database). Water Survey of Canada (2025).

22. Lehner, B., Messager, M. L., Korver, M. C. & Linke, S. Global hydro-environmental lake characteristics at high spatial resolution. Sci. Data 9, 351 (2022).

23. Falcone, J. A. GAGES-II: Geospatial Attributes of Gages for Evaluating Streamflow. U.S. Geological Survey https://doi.org/https://www.sciencebase.gov/vocab/category/item/identifier (2011).

24. Environment and Climate Change Canada. National hydrometric network basin polygons. (2016).

25. Linke, S. et al. Global hydro-environmental sub-basin and river reach characteristics at high spatial resolution. Sci. Data 6, 283 (2019).

26. Thornton, M., Shrestha, R., Thornton, P. E. & Kao, S.-C. Daymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 4 R1 | NASA Earthdata. ORNL Distributed Active Archive Center https://doi.org/https://doi.org/10.3334/ORNLDAAC/2129 (2025).

27. Myneni, R. B., Knyazikhin, Yuri & Park, T. MODIS/Terra Leaf Area Index/FPAR 8-Day L4 Global 500m SIN Grid V061 | NASA Earthdata. NASA Land Processes Distributed Active Archive Center https://doi.org/https://doi.org/10.5067/MODIS/MOD15A2H.061 (2025).

28. Friedl, M. A. & Sulla-Menashe, D. MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid V061 | NASA Earthdata. NASA Land Processes Distributed Active Archive Center https://doi.org/https://doi.org/10.5067/MODIS/MCD12Q1.061 (2022).

29. Anderson, J. R., Hardy, E. E., Roach, J. T. & Witmer, R. E. A Land Use and Land Cover Classification System for Use with Remote Sensor Data. doi:10.3133/pp964.

30. U.S. Geological Survey. Annual NLCD Collection 1 Science Products: U.S. Geological Survey data release. Earth Resources Observation and Science (EROS) Center https://doi.org/https://doi.org/10.5066/P94UXNTS.

31. DeCicco, L. et al. dataRetrieval: R packages for discovering and retrieving water data available from U.S. federal hydrologic web services. U.S. Geological Survey https://doi.org/https://doi.org/10.5066/P9X4L3GE (2026).

32. Albers, S. J. tidyhydat: Extract and Tidy Canadian Hydrometric Data. J. Open Source Softw. 2, 511 (2017).

33. Falcone, J. A. U.S. Geological Survey GAGES-II time series data from consistent sources of land use, water use, agriculture, timber activities, dam removals, and other historical anthropogenic influences. U.S. Geological Survey https://doi.org/10.5066/F7HQ3XS4 (2017).

34. Lehner, B. & Grill, G. Global river hydrography and network routing: baseline data and new approaches to study the world’s large river systems. Hydrol. Process. 27, 2171–2186 (2013).

35. Jones, K. A., Niknami, L. S., Buto, S. G. & Decker, D. Federal Standards and Procedures for the National Watershed Boundary Dataset (WBD). 54 https://pubs.usgs.gov/tm/11/a3/ (2022).

36. U.S. Geological Survey. gdptools. U.S. Geological Survey Water Mission Area.

37. Knyazikhin, Y., Martonchik, J. V., Myneni, R. B., Diner, D. J. & Running, S. W. Synergistic algorithm for estimating vegetation canopy leaf area index and fraction of absorbed photosynthetically active radiation from MODIS and MISR data. J. Geophys. Res. Atmos. 103, 32257–32275 (1998).

38. McMillan, H. K. A review of hydrologic signatures and their applications. WIREs Water 8, e1499 (2021).

39. Petersky, R. & Harpold, A. Now you see it, now you don’t: a case study of ephemeral snowpacks and soil moisture response in the Great Basin, USA. Hydrol. Earth Syst. Sci. 22, 4891–4906 (2018).

40. Hatchett, B. J. Seasonal and Ephemeral Snowpacks of the Conterminous United States. Hydrology 8, 32 (2021).

41. Hecht, J. S., Johnson, Z. C. & Dunne, K. A. Long-term monotonic trends in annual and monthly streamflow metrics at streamgages in the United States (ver. 2.0, October 2024). U.S. Geological Survey https://doi.org/10.5066/P9VBR38I (2024).

42. Omernik, J. M. & Griffith, G. E. Ecoregions of the Conterminous United States: Evolution of a Hierarchical Spatial Framework. Environ. Manage. 54, 1249–1266 (2014).

43. Corak, N. K., Thornton, P. E. & Lowman, L. E. L. A high resolution, gridded product for vapor pressure deficit using Daymet. Sci. Data 12, 256 (2025).

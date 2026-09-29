# HISSS Manuscript Draft (auto-synced snapshot)

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/d/e/2PACX-1vS7j4FRp7SEwlXoBUVA8NA7cj_I0XzyS0u58r3bl8SOz4BfpZPrdPJge4RMcFocnX8Gnllkc1M-CTJ3/pub)
> **Last synced**: 2026-09-29 (110 paragraphs vs 105 on 2026-09-10). **Major co-author
> revision**: front-matter to-do lists deleted; abstract added; title line now
> "Hydroclimate Information, Signatures and Summary Statistics" (body still says
> "Hydrologic …" four times); S.F. Dymond added to the author list; §1 tightened (old
> para 5 removed); §2 preamble rewritten (eight signature families; trends windows);
> **§2.1.1 Streamflow and §2.1.2 Boundaries swapped** (cross-references not updated);
> new pre-2.1.1 paragraph with counts; §2.1.3 Daymet rewritten (gap explained as
> "processed before the streamflow requirements were finalized"; "6,087 or 6,041");
> §2.1.4 HydroATLAS expanded; §2.2.1 uses the eight families and absorbs the seasonal-flag
> and gage-inclusion sentences; §2.2.2 states the constant-flow rule explicitly; §3
> citations resolved (no "tbd" left); **§4 drafted** (ecoregions, size quartiles, Figures
> 2–3 captions); Data/Code availability + Funding sections added; references reformatted
> and extended (27 → 39). §5, Acknowledgements ("Claude Code 0.145.0") and Disclaimers
> unchanged. Methods claims that now DISAGREE with the code are logged in
> `docs/reconciliation/manuscript_log.md` (2026-09-29).
>
> This is a read-only snapshot of the collaborative manuscript draft used for
> change detection and reconciliation review (the `/sync-docs` skill). Do NOT hand-edit the body below the header — it is overwritten at
> each sync. Discrepancies between the manuscript's methods and the code/docs
> are logged in `docs/reconciliation/manuscript_log.md`; manuscript-side corrections must be made in the Google Doc by the
> co-authors (@ Arik convention per the draft's own notes).

---

Journal: Scientific Data - Special Collection on Water Storage due Nov 9, 2026

Draft Title: Harmonized dataset for analyzing flows in the critical zone

HISSS: Hydroclimate Information, Signatures and Summary Statistics

Kaiser, K.E., A. Tashie, L. Lowman, K. Jennings, G. Gorski, D. Murray, Dymond, S.F.

FILL out Authorship contributions including Hydroshare user name; Update author list on hydroshare collection

Deleted Material

Active Signatures Doc

Abstract (modified from AGU)

The hydroclimate across North America is changing, with widespread transitions from snow-dominated to increasingly rain-dominated precipitation regimes. While such shifts in snowpack and winter hydroclimate are well documented, how these changes propagate streamflow behavior is neither uniform nor well characterized across climate, physiography, and critical zone structure, partly because of difficulties with standardizing and harmonizing disparate datasets. To meet this need, we present HISSS: Hydrologic Information Signatures and Summary Statistics - a new integrated dataset that includes annual snow, climate, and streamflow data for 8,014 North American watersheds, with trends and physiographic attributes for 6,041-7,964 locations depending on the subset. These products are intended for critical-zone scientists and beyond, who study watershed processes, snow hydrology, and vegetation dynamics, among others, and provide an extensive dataset for exploring changing hydroclimate and streamflow responses, and how these relationships play out across various physiographies, climate zones, and critical zone structures.

1 Background & Summary

The critical zone is the region of the Earth’s surface and subsurface where water is stored and cycled in the environment and is used by humans (Lin 2010). Water in the critical zone is stored and exchanged in and between soil, vegetation, snow, surface reservoirs, and groundwater. At the Earth’s surface, streamflow is an important metric of water storage and availability that indicates the severity of stress for human and natural systems (Brooks et al. 2015; Wlostowski et al. 2021). The composition of the critical zone determines water availability for streamflow dynamics under normal conditions, disturbances (e.g., drought, fire, land conversion), and increasing atmospheric temperatures and CO2 (Keller et al. 2019; Condon et al. 2020). Open questions in Earth system science include how system characteristics (i.e., soil texture and depth, vegetation type and rooting depth, hydroclimate, etc.) mediate water availability and surface and subsurface routing in the critical zone. Understanding the critical zone’s control on streamflow responses to climate is challenging because subsurface and aboveground critical zone components vary spatially and temporally (Lin 2010; Brooks et al. 2015). Further, datasets describing critical zone features and climate are not readily available at spatial and temporal scales relevant for assessing impacts to streamflow. The Hydrologic Information Signatures and Summary Statistics (HISSS) dataset presented in this manuscript overcomes these issues by harmonizing climate, vegetation, soil, and human impact data at the watershed scale for >6,000 streamgages across the United States and Canada. HISSS provides a comprehensive streamflow signature dataset calculated from streamflow and gridded precipitation records spanning the years 1980 through 2025.

Many studies have considered how the critical zone mediates water exchange between surface and sub-surface storage components (Tague et al. 2025). For instance, soil depth and texture influence the partitioning of rain and snow in the critical zone (Hammond et al. 2019), as well as runoff generation (Zimmer and Gannon 2018; Chen et al. 2024). The role of vegetation in mediating runoff responses has been explored previously for different drought scenarios (Corak et al. 2024), with regard to the loss of global forests (Zhang et al. 2017) and land conversion in tropical forests (Muñoz-Villers and McDonnell 2013). Snow depth, snow water equivalent, and snowmelt timing influence runoff and groundwater recharge in the critical zone (Wlostowski et al. 2021). However, while prior studies explore relationships at individual watersheds or specific ecosystems, data are not available at continental scales that allow for cross-comparisons and syntheses across climate zones and ecoregions.

One challenge in harmonizing climate and vegetation datasets with streamflow is mismatches in temporal and spatial scales (Vlah et al. 2023). Streamflow represents point-scale measurements that integrate water fluxes across a watershed, while climate and vegetation data are typically provided as gridded products representing average values across a defined pixel at a point in time. Several large-sample datasets have addressed this mismatch by aggregating gridded and geospatial data to the watershed scale, with products focused on hydroclimate and hydrological signatures (Newman et al. 2015; Addor et al. 2017), biogeochemistry (Vlah et al. 2023), stream chemistry (Sterle et al. 2024), as well as broad multi-source archives that now cover much of North America (HYSETS, Arsenault et al. 2020; Kratzert et al. 2023). Among the products that report hydrological signatures, data are given only as single period-of-record values rather than as a time series. Thus, these products do not provide per-signature trend or changepoint statistics, and the harmonized multi-region archives typically derive their climate forcing from comparatively coarse reanalyses (e.g., ERA5-Land at ~11 km or coarser). The recent CAMELS-SPAT (Knoben et al. 2025) dataset spans the United States and Canada and includes high spatial resolution climate data, time-varying MODIS leaf area index (LAI), and land use and land cover (LULC). However, it is focused on a subset of available watersheds (n < 1500) to support hydrological model development and therefore precomputes only a narrow subset of hydroclimatological signatures as period-of-record statistics.

The HISSS dataset and software library presented here provide holistic hydroclimate, vegetation, snow, and geophysical information at the scale of individual watersheds across the continental United States and Canada. For over 6,000 watersheds in the study, we calculate roughly 100 hydrological signatures and report the underlying annual values along with a consistent suite of computed metrics: central tendency statistics (mean, median), linear and robust (Theil-Sen) trend slopes, Spearman and Mann-Kendall trend significance, Pettitt changepoint detection, and other supporting diagnostics. These are paired with annually resolved Moderate Resolution Imaging Spectroradiometer (MODIS) and National Land Cover Dataset (NLCD) LULC and monthly resolved MODIS LAI, watershed climate features derived from ~1 km Daymet, and HydroATLAS physiographic attributes. Because every signature is computed by an open, cross-language library (Julia, Python, and R), each metric is reproducible and extensible. These products are intended for critical-zone scientists and beyond, who study watershed processes, snow hydrology, and vegetation dynamics, among others.

2 Methods

To build the HISSS dataset, we accessed numerous time-series and geospatial datasets across the United States (US) and Canada using an open-source, reproducible processing workflow (Figure 1). We assembled daily mean streamflow records for 16,994 gages across the United States (9,154; USGS Water Data for the Nation) and Canada (7,840; HYDAT), of which 8,014 (6,160 USGS; 1,854 HYDAT) yielded usable daily records for watersheds smaller than 85,000 km2, totalling more than 111 million observations. We defined the watershed boundary for each gage using official agency basin polygons and the HydroBASINS dataset (Sect. 2.1.1). These polygons defined the basin boundaries over which we aggregated gridded climate, land-cover, and basin-attribute data (Sects. 2.1.3–2.1.4). For each basin and its suite of hydroclimatic and geospatial data, we computed 121 signatures across eight signature families (drought, flow volume, flashiness, precipitation-streamflow, flow duration, flow timing, snow, storage): 100 annually resolved signatures, each summarized with a common set of 16 statistics, plus 21 per-gage diagnostics. Trends on these signatures and statistics were calculated for 1980-2025 (6,250 sites) and for 1993-2025 (6,678 sites) (Sect. 2.1.2). The following subsections detail the input data, filtering procedures, metric calculations, and trend analyses. The code for the full workflow is available in an open library in Julia, Python, and R (https://github.com/CZ-Sync/HISSS).

Figure 1. Workflow used to generate the Hydrologic Information Signatures and Summary Statistics (HISSS) dataset.

2.1 Input data

Table 1. Summary of data sources used to generate the Hydrologic Information Signatures and Summary Statistics (HISSS) and the number of sites, temporal resolution, and temporal coverage included.

<TABLE OF ALL INPUT DATA SOURCES>

From the candidate 16,994 gages with streamflow data, we filtered to those that had at least 20 years of valid streamflow data with >30 days over .0001 mm/day, yielding 8,014 gages across the United States and Canada. Of these gages, 7,964 watersheds had polygons and were smaller than 85,000 km (≤ largest HUC8 size). The total number of gages with streamflow, Daymet and MODIS LULC is 5,965.

2.1.1 Streamflow data

We retrieved daily mean observed streamflow from 01 October 1979 to 30 September 2025 for all candidate gages (n = 16,994) using the USGS’s dataRetrieval R package (DeCicco et al., 2026) for US gages and the Water Survey of Canada’s tidyhydat R package (Albers et al., 2026) for Canadian gages. To ensure we retained only the most useful observations for time series analysis, we applied standardized quality assurance criteria at three stages during the calculation of annualized data: 1) retrieval, 2) compilation, and 3) per-year qualification. Gages that returned fewer than 20 total years of valid daily observations were excluded, leaving 8,014 gages (6,160 US; 1,854 Canadian) with usable daily records. These data were compiled into a single dataset of 111.6 million observations keyed to the respective agency’s gage ID. We retained original agency quality flags, including ice-affected and estimated values. We stored missing days as NA rather than zero-filled, treated zero flow as a valid observation, and area-normalized discharge (mm d-1) using agency-published drainage areas for all but 73 gages, for which no drainage area is published (see Usage Notes). We aggregated each gage’s record to a continuous daily time series and applied the following data quality screening by water year (WY; October 01 - September 30): gaps of up to three consecutive days were filled by linear interpolation, any WY with more than 30 missing values were set to NA, and if any data gap exceeded three consecutive days those days were set to NA.

2.1.2 Watershed Boundaries

To define the spatial extents of the contributing area to the gaged locations, we combined official agency basin polygons from the US and Canada. For the US, we used watershed polygons from the USGS GAGES-II dataset (Falcone 2011) and merged the individual shapefiles for reference and non-reference gages across all US regions into a single US-wide shapefile (covering 100% of US gages). For Canada, we used the Water Survey of Canada’s Hydrometric Network Basin Polygons dataset (Environment and Climate Change Canada, 2016) and merged them into a single Canada-wide layer (covering 98.3% of Canadian gages). Gages lacking an official polygon were delineated from the HydroBASINS level-12 drainage network (Lehner and Grill 2013; detailed methods below) by aggregating all basins upstream of the gage’s outlet. We leveraged HydroATLAS (Linke et al. 2019) to define static basin attributes within the data product domain at the 8,014 gages with suitable streamflow records. The source of the watershed geometry was recorded for each watershed. We excluded basins that exceeded 85,000 km2, transformed the remaining watershed boundaries into a common projection, and then merged them into one full-domain layer. We used these boundaries to aggregate remotely sensed LULC and LAI data, and then aggregated gridded climate data only where we had sufficient streamflow data coverage for trends analyses (see Section 2.1.3).

2.1.3 Climate and snow data

Meteorological data came from Daymet (Thornton et al. 2022), a daily long-term 1 km x 1 km gridded continuous product covering North America. Daily Daymet estimates of minimum and maximum air temperature, precipitation, vapor pressure, shortwave radiation, snow water equivalent, and day length covering 1980–2023 were downloaded from the ORNL Distributed Active Archive Center (DAAC, Thornton et al. 2022). Daymet aggregation was performed for a subset of the gages with usable streamflow records, these were processed before the streamflow time series requirements were finalized, which created a discrepancy between sites with polygons (7,964) and those with basin-averaged Daymet series (6,087 or 6,041); climate-, snow-, and precipitation-dependent signatures are reported as NA for the remaining basins. The watershed size threshold corresponds to the largest HUC8 basin in the Watershed Boundary Dataset (Jones et al. 2022), which was used to ensure that hydrologic response to weather was meaningfully well represented with basin-averaged climate inputs. We performed this aggregation using the gdptools package (U.S. Geological Survey 2026), which performs area-weighted zonal aggregation of gridded data to polygon features across all variables and time steps. The result is a single area-weighted daily value for each variable for each basin.

2.1.4 Basin attributes and land cover

The HISSS dataset includes both static physiographic characteristics and landcover characterization through time (Table 3). We used static basin attributes from HydroATLAS (BasinATLAS version 10; Linke et al. 2019), a global compilation of hydro-environmental descriptors from HydroBASINS level-12 subbasins (Lehner and Grill 2013). HydroBASINS is a global, hierarchically nested set of subbasin polygons derived from a 15 arcsecond drainage network. Level 12 is its finest tier (mean area of 110 km2), and each unit carries a downstream neighbor identifier so that a gage’s full contributing area can be assembled by aggregating the upstream network. For every gage, we identified the intersecting level-12 subbasin and the complete set of upstream subbasins, then aggregated the HydroATLAS attributes over that contributing area. Continuous and percentage attributes were aggregated as area-weighted means, elevation as the spatial minimum and maximum, and categorical attributes as the area-weighted majority. The resulting table is keyed to gage ID so that it joins directly to streamflow signatures. For headwater basins smaller than the intersecting level-12 HydroATLAS polygons (average ~ 110 km2), the static basin attributes include downstream components of the larger watershed, in addition to the headwater basin of interest. In total, we derived approximately 210 attributes per watershed spanning six thematic categories: hydrology, physiography, climate, LULC, soils and geology, and anthropogenic influence.

Table 3. All Landcover datasets included in HISSS {tab in the tables xls}

To obtain information on how land surface characteristics change through time, the HISSS dataset also includes all eight land cover classification schemes from MODIS MCD12Q1 Land Cover Type Yearly product and LAI from the MODIS MCD15A3H, both v061 (Myneni et al. 2021). Both have a native spatial resolution of 500 m and global coverage, with MCD12Q1 produced as annual composites and MCD15A3H as 4-day composites. MCD12Q1 classifies each pixel under eight schemes, all of which HISSS retains: five Land Cover Type schemes based on the International Geosphere-Biosphere Programme (IGBP), University of Maryland (UMD), LAI/fPAR biome, Biome/Biogeochemical Cycles (BGC), and Plant Functional Types (PFT), and 3 FAO Land Cover Classification System property layers describing land cover, land use, and surface hydrology, respectively. Together, these yield 102 per-class percent-coverage values for each gage-specific watershed and year.

The MODIS MCD15A3H LAI product represents the seasonal phenology of vegetation on the land surface. It is estimated from a look-up table derived from a 3-D radiative transfer model relating surface reflectance to observed LAI values (Knyazikhin et al. 1998). It also relies on the MODIS Land Cover Type 3 (LAI Classification) to determine the most likely LAI value for a pixel given its vegetation type. MODIS LAI has valid values between 0 and 10 m2 m-2, which describe the one-sided projection of leaf area per unit ground area.

To complement the MODIS land cover record, HISSS also includes land cover from the USGS Annual National Land Cover Database (Annual NLCD, Collection 1), a Landsat-derived product mapped at 30 m resolution annually from 1985 through 2025. For each gaged watershed we provide annual percent coverage of the 16 land cover classes together with basin-mean fractional impervious surface, a direct measure of urbanization intensity that the MODIS products do not offer. Relative to MCD12Q1, Annual NLCD extends the land cover record 16 years further back and resolves land cover heterogeneity within small watersheds that 500 m MODIS pixels may miss. Its principal limitation is its spatial coverage: Annual NLCD spans CONUS only, so it is provided for only those 6,119 watersheds while Alaskan and Canadian watersheds retain MODIS-only land cover. The two products are intended to be complementary, with MCD12Q1 supplying consistent coverage across all HISSS watersheds, while Annual NLCD supplies a longer, higher-resolution record for CONUS watersheds.

2.2 Metrics for analysis

2.2.1 Computing Streamflow and Hydroclimate Signatures

We calculated a comprehensive suite of 121 streamflow and hydroclimate signatures for each gage and its associated watershed, comprising 100 annually resolved metrics with statistics and 21 stand-alone diagnostic metrics. The selection of metrics is consistent with analysis in catchment hydrology (McMillan 2021) and snow hydrology (Petersky and Harpold 2018; Hatchett 2021) and span the eight signature families: flow volume, flow duration, storage, flashiness, drought, flow timing, precipitation-to-streamflow, and snow. Details for each metric, including the definition, requirements, and relevant citations, are in the data dictionary (hisss_data_dictionary.csv) included with Resources 1 and 2. Seasonal periods are defined as: Winter: December-February, Spring: March-May, Summer: June-August, Fall: September-November. Annual metrics were calculated on the water year October 1 - September 30.

[MOVED from above]

We further flagged any WY-season (fall: September - November, winter: December - February, spring: March - May, summer: June - August) with fewer than 80% of raw observations passing quality assurance, setting the affected seasonal metrics to NA for that WY.

Finally, we calculated signatures for gages only if they retained at least 20 qualifying WYs and spanned at least 60% of the possible WYs in the analysis window (additional per-signature completeness requirements determine whether trend statistics are reported for a given signature, as described in Section 2.2.2). We calculated signatures (Sect. 2.2) over two standard analysis windows (WYs 1993-2025 and 1980-2025), yielding 6,678 gages for WYs 1993–2025 and 6,250 for WYs 1980–2025. The shorter time range (1993-2025) was included as it is the date range with the most gages with the same data record length, allowing for direct comparison of trends across a standardized date range.

2.2.2 Computing Trends and Summary Statistics

For the trends analysis, we removed years with more than 3 consecutive days of NAs or more than 30 days of NAs throughout the water year. We flagged (but did not remove) water years that contained negative streamflow values or any calendar month in which at least 15 days of non-zero streamflow held at a single constant value, patterns that are indicative of data errors or highly regulated discharge. A gage was included only if it had at least 20 qualifying water years across its full period of record and at least 60% of the water years in the period of interest. For a given signature, trend statistics were then computed only if its annual series had at least 20 valid values, at least 60% of the entire series had valid values, and at least 80% of both the first and last decades of the trend had valid values; the event-based recession and elasticity metrics are exempt from these completeness requirements. Any seasonal metric additionally required at least 80% of the period’s observations had valid values. We also computed 12 automated data-quality flags (range and consistency checks) on the resulting signatures.

3 Data Record

The HISSS dataset is published in HydroShare, the repository for water data operated by the Consortium of Universities for the Advancement of Hydrological Science, Inc. (CUAHSI; https://www.hydroshare.org). The dataset is organized as a public collection (DOI: http://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768) grouping five resources, each independently citable with its own DOI (Table X). All resources are released under a CC-BY 4.0 license, and each contains a README describing its files and a machine readable data dictionary defining every column. Tabular data are provided as UTF-8 CSV where file sizes permit and as Apache Parquet for large tables; Parquet files are readable with standard open-source libraries in R (arrow), Python (pyarrow/pandas), and Julia. The processing code is archived separately (see Code Availability).

Table [2]: hisss_data_resources_table

Resources 1 and 2 (streamflow signatures). These parallel resources implement the two analysis windows described in Section 2.1.2 (water years 1993-2025 and 1980-2025); gages qualify independently under each window’s completeness criteria, so neither product is a subset of the other. Each resource contains: (i) a signature summary table (hiss_signatures_wy{window}.csv) with one row per gage and 1,653 columns comprising gage metadata, the 100 signatures x 16 statistics each, 21 signatures and related outputs that do not carry statistics, and 12 automated quality-assurance flags; (ii) a long-format annual-values table (Parquet) holding the yearly value behind every signature-statistic pair, with columns gage_id, signature, water_year, and value (18,898,406 and 24,366,487 rows for the two windows, respectively); (iii) a self-contained interactive HTML explorer that maps every signature-statistic combination and plots per-gage annual series with fitted checksums, configuration, and software versions) and validation reports.

Resource 3 (harmonized input data). Daily discharge for all 8,014 processed gages (111,624,189 records; 1980 through 2025) compiled from the USGS National Water Information System via dataRetrieval (DeCicco et al. 2026) and the Water Survey of Canada HYDAT database via tidyhydat (Albers 2017; ECCC 2025), expressed in mm d-1 (see area_normalized below). The accompanying gage metadata table (16,994 candidate gages) records station coordinates, published drainage areas, processing status, and human-interference indicators drawn from GAGES-II attributes (Falcone 2011; 2017) and Reference Hydrometric Basin Network (RHBN) and regulation status via HYDAT. Daily Basin-averaged Daymet Version 4 R1 (Thornton et al. 2022) for the variables precipitation, temperature, snow water equivalent (SWE), vapor pressure, and shortwave radiation is provided for 6,087 basins (97,757,220 records, calendar years 1980-2023) as produced by the area-weighted aggregation described in Section 2.1.3.

Resource 4 (watershed geometry and basin attributes). Watershed boundary polygons for 7,964 gaged watersheds (GeoPackage and GeoParquet, WGS84), assembled from GAGES-II basin boundaries (Falcone 2011), the Water Survey of Canada hydrometric basin polygons (Environment and Climate Change Canada 2016), and HydroBASINS-derived delineations where no official polygon exists (Lehner and Grill 2013); the provenance of each polygon is recorded in watershed_geom_source. A companion table provides 211 watershed-scale attributes (climate, hydrology, terrain, land cover, soils, and geology, and anthropogenic influence) aggregated over each gage’s full upstream area from BasinAtlas v10 (Linke et al. 2019; Lehner et al. 2022), with a data dictionary mapping each column to its theme, unit, and aggregation method.

Resource 5 (vegetation and land cover). Three per-watershed time-series tables keyed to watershed polygons (Resource 4): (i) monthly MODIS LAI (MCD15A3H v061; Myneni et al. 2015) for 2002-2024 (2,150,280 rows; 7,964 watersheds x 270 months) with basin-mean monthly LAI, spatial-distribution statistics (quantiles, standard deviation, extrema), and per-month quality fractions; (ii) annual MODIS land cover (MCD12Q1 v061; Friedl and Sulla-Menashe 2015) for 2001-2024 (191,136 rows with 102 per-class percent-coverage columns across all eight classification schemes, as described in Section 2.1.4); and (iii) annual NLCD land cover and fractional impervious surface (Annual NLCD Collection 1) for 1985-2025 (250,879 rows covering 6,119 CONUS watersheds, with 16 per-class percent-coverage columns and basin-mean imperviousness. Each table is accompanied by its data dictionary, granule-level provenance manifests for the MODIS products, and a self-contained interactive HTML explorer.

Conventions and key fields. All tables join on gage_id, the agency station identifier (USGS site numbers; HYDAT station numbers); read it as a character string and join with leading zeroes removed; the geometry and land-cover tables also carry this form as canon_id. Files derived from the geometry layer also carry canon_id, a zero-stripped variant (for USGS sites) retained for joins to legacy metadata. water_year denotes the period 1 October - 30 September, labeled by the ending calendar year. Signature columns in Resources 1-2 follow the pattern {signature}_{statistic}, where the sixteen statistic suffixes are: _senn_slp (Theil-Sen slope), _linear_slp (ordinary least-squares slope), _spearman_rho / _spearman_pval (Spearman rank correlation with time and its p-value), _mk_rho / _mk_pval (Mann-Kendall tau and p-value), _mean, _median, and eight Pettitt changepoint fields beginning with _pettitt including _cp_year (most likely changepoint year), _pval, _pre_mean / _post_mean / _delta_mean / _pct_change (segment means and their differences), and _pre_mk_pval / _post_mk_pval (within-segment trend tests). A small number of signatures are single-valued per gage s(e.g., elasticity_static, the recession-derived filter constant, and diagnostic counts) and appear without the suffix set; these are enumerated in the data dictionary. In the annual-values tables, a missing row and a stored NaN are equivalent, both meaning “not computable for that gage-year.” The Boolean area-normalized flag (see unnormalized gages section 5.1.1) identifies whether a gage’s discharge is expressed in mm d-1 (true) or retained in native m3 s-1 because no published drainage area exists (false). Columns prefixed flagged_for_ are automated range- and consistency-check indicators defined in the dictionary. Land cover columns follow {scheme}_c{code}_pct (MODIS) and nlcd_c{code}_pct} (NLCD), with pct indicating the percent of the watershed area that is in the class code.

4 Data Overview

Sites included in the HISSS dataset span a broad geographic range across the United States and Canada, covering 12 Level I North American Ecoregions (Omernik & Griffith 2014). A majority (62%) of the sites are located in the Eastern Temperate Forests or the Northwestern Forested Mountains, while there are fewer than 10 sites within each of the Southern Semiarid Highlands, Tundra, and the Hudson Plain ecoregions. HISSS sites comprise small basins encompassing headwater streams up to large river basins with a median watershed size of 698 km2, 25% of sites having drainage areas greater than 2,700 km2, and 25% of sites having drainage areas less than 190 km2. Streamflow signature data coverage was lowest for the family of storage metrics (Figure 3A). The majority of sites have greater than 80% streamflow signature coverage across all families of signatures, with the lowest coverage at sites within the ecoregions taiga and Husdon plain (Figure 3B).

Figure 2. Distribution of sites with delineated catchment areas and paired streamflow and meteorological data for the 1993-2025 period; sites are color-coded by their Level I North American Ecoregion (Omernik and Griffith 2014) with the number of sites falling within each ecoregion shown in parentheses in the map legend.

Figure 3. Boxplots showing distribution of coverage among sites. For both panels, coverage is calculated as the number of metric-years/total possible metric-years for each site, and darker boxes show distributions for the time period 1980-2025, while lighter boxes show distributions for 1993-2025. A) Coverage by signature family, where the number below each signature family is the number of signatures within that family. B) Coverage by ecoregions, beneath each ecoregion label, the number of sites within that ecoregion is shown for 1980-2025 (first number) and 1993-2025 (second number), which add up to 6,250 and 6,678, respectively.

5 Usage Notes

5.1.1 Unnormalized gages

73 of the 8,014 processed gages have no agency-published drainage area (32 and 28 of them qualify for the WY 1993–2025 and WY 1980–2025 products, respectively). Most are irrigation or diversion canals, dam and powerhouse outflows, or channel splits of large rivers, where a contributing area is undefined or unpublished. These gages are retained with discharge in native m3 s-1 and flagged area_normalized = FALSE. Their unit-carrying signatures (flow volumes, flow percentiles, recession log(a), drought deficits) are not comparable with the mm d-1 values at other gages, and their precipitation-dependent signatures (runoff ratios, elasticity, Q-P seasonality, storage) are structurally NA. Users should filter on area_normalized == TRUE before any cross-gage comparison of unit-carrying signatures.

5.1.2 How to join HISSS with other datasets

Every HISSS table is keyed by the agency station identifier named “gage_id” (USGS site number, e.g. 0103500; HYDAT station number, e.g. 01AD002). This station identifier or gage_id may be read as a character string and joined on the leading-zero-stripped form on both sides (most US site numbers begin with a zero that numeric parsing silently ignores (see resource READMEs for details). Any dataset keyed to the same agency identifiers then joins directly. Datasets that are not gage-keyed join spatially through the watershed polygons of Resource 4: gridded products such as the Daymet-derived VPD of Corak et al. (2025), which shares Daymet's 1 km grid and 365-day calendar, can be aggregated with the area-weighted workflow used for Daymet (Sect. 2.1.3), and point records such as MacroSheds sites (Vlah et al. 2023) can be located within the nested polygons or matched to the nearest gage by coordinates. HydroATLAS-based products join on the HydroBASINS level-12 outlet identifier (Downstream_HB_ID) in the Resource 4 attribute table. When combining data, aggregate partner records to the water year (October–September, labeled by the ending year); re-normalize fluxes to a common drainage area, since HISSS reports both the agency-published area used for normalization (basin_area) and the polygon area (geom_area_km2); and pair record-dependent signatures (Sect. 5.1.3) with the product whose window matches the partner record.

5.1.3 Record-dependent signatures

Hydroclimate signatures whose definition uses thresholds or means from the full analysis window are valid within the product for that time window and should not be compared against the WY 1993–2025 and WY 1980–2025 products, nor re-derived from the annual values over a different window. The specific signatures include the period-of-record (*_all) pulse metrics, elasticity, the recession-parameterized baseflow indices, and all drought metrics.

Data availability

The HISSS dataset is openly published in HydroShare at http://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768 licensed under a Creative Commons Attribution 4.0 International (CC BY 4.0) licenceCC-BY 4.0 license.

Code availability

The code for the full workflow is available in an open library in Julia, Python, and R, https://github.com/CZ-Sync/HISSS. Details for each metric, including the definition, requirements, and relevant citations, are in the Streamflow Signatures Reference on GitHub (https://github.com/CZ-Sync/HISSS/blob/main/docs/SIGNATURES.md). And data resources and their associated metadata and guidance for usage are in https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/ [TEMPORARY: WILL CHANGE WHEN WE PUBLISH].

Acknowledgements

AI-assisted coding tools (Claude Code 0.145.0, Anthropic) were employed to generate code used in data ingestion and processing, hydrological signature extraction, cross-language validation and benchmarking, and interactive visualization. All code was reviewed, tested, and validated by the authors to ensure correctness and reproducibility. Generative AI tools (Claude Code 0.145.0, Anthropic) were used to support data analysis and visualization. These tools were applied under the supervision of the authors, and all outputs were reviewed and validated against established scientific methods to ensure reproducibility and transparency.

Funding

Powell Center, others?

Disclaimers

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

References

Addor, N., Newman, A. J., Mizukami, N., & Clark, M. P. The CAMELS Data Set: Catchment Attributes and Meteorology for Large-Sample Studies. Hydrol. and Earth Syst. Sci. 21, 10, 5293–313 (2017). https://doi.org/10.5194/hess-21-5293-2017.

Albers, S. tidyhydat: Extract and Tidy Canadian Hydrometric Data. J. Open Source Softw. , 2, 20, 1-4 (2017). https://doi.org/10.21105/joss.00511.

Arsenault, R., Brissette, F., Martel, J. L.,Troin, M., Lévesque, G., Davidson-Chaput, J., Castañeda Gonzalez, M., Ameli, A., & Poulin, A. A comprehensive, multisource database for hydrometeorological modeling of 14,425 North American watersheds. Sci. Data 7, 243 (2020). https://doi.org/10.1038/s41597-020-00583-2.

Brooks, P. D., et al. Hydrological Partitioning in the Critical Zone: Recent Advances and Opportunities for Developing Transferable Understanding of Water Cycle Dynamics. Water Resour. Res. 51, 9, 6973–87 (2015). https://doi.org/10.1002/2015WR017039.

Chen, H., Niu, Q., McNamara, J. P., & Alejandro N. Flores, A. N. Influence of Subsurface Critical Zone Structure on Hydrological Partitioning in Mountainous Headwater Catchments. Geophys. Res. Lett. 51, 6, e2023GL106964 (2024). https://doi.org/10.1029/2023GL106964.

Condon, L. E., et al. Where Is the Bottom of a Watershed? Water Resour. Res. 56, 3, e2019WR026010 (2010). https://doi.org/10.1029/2019WR026010.

Corak, N. K., Otkin, J. A., Ford, T. E., & Lowman, L. E. Unraveling Phenological and Stomatal Responses to Flash Drought and Implications for Water and Carbon Budgets. Hydrol. and Earth Syst. Sci., 28, 8, 1827–51 (2024). https://doi.org/10.5194/hess-28-1827-2024.

Corak, N. K., Thornton, P. E. & Lowman, L. E. A High Resolution, Gridded Product for Vapor Pressure Deficit Using Daymet.”Sci. Data 12, 1, 256 (2025). https://doi.org/10.1038/s41597-025-04544-5.

DeCicco L., Hirsch R, Lorenz D, Read J, Walker J, Platt L, Watkins D, Blodgett D, Johnson M, Krall A, Stanish L, Zemmels J, Hinman E, & Mahoney M. dataRetrieval: R packages for discovering and retrieving water data available from U.S. federal hydrologic web services. (2026) U.S. Geological Survey. https://doi.org/10.5066/P9X4L3GE.

Environment and Climate Change Canada (ECCC). National hydrometric network basin polygons. Water Survey of Canada (2016). https://open.canada.ca/data/en/dataset/0c121878-ac23-46f5-95df-eb9960753375

Environment and Climate Change Canada (ECCC). HYDAT: National Water Data Archive (hydrometric database), release 2025-10-14 [SQLite database]. Water Survey of Canada (2025). https://wateroffice.ec.gc.ca/ (retrieved 7 February 2026 via tidyhydat).

Falcone, J. GAGES-II: Geospatial Attributes of Gages for Evaluating Streamflow: U.S. Geological Survey data release (2011) https://doi.org/10.5066/P96CPHOT.

Falcone, J. A. US Geological Survey GAGES-II time series data from consistent sources of land use, water use, agriculture, timber activities, dam removals, and other historical anthropogenic influences. US Geological Survey (USGS) Data Release, p.740 (2017). https://doi.org/10.5066/F7HQ3XS4

Friedl, M. & Sulla-Menashe, D. Mcd12q1 modis. Terra+ aqua land cover type yearly l3 global 500m SIN grid 6 (2015).

Hammond, J. C., Harpold, A. A., Weiss, S., & Kampf, S. K. Partitioning Snowmelt and Rainfall in the Critical Zone: Effects of Climate Type and Soil Properties. Hydrol. and Earth Syst. Sci. 23, 9, 3553–70 (2019). https://doi.org/10.5194/hess-23-3553-2019.

Hatchett, B. J. Seasonal and Ephemeral Snowpacks of the Conterminous United States. Hydrol. 8, 1, 32. https://doi.org/10.3390/hydrology8010032.

Jones, K. A., Niknami, L. S., Buto, S. G., & Decker, D. Federal standards and procedures for the national Watershed Boundary Dataset (WBD) (5 ed.): U.S. Geological Survey Techniques and Methods 11-A3, 54 p. (2022). https://pubs.usgs.gov/tm/11/a3/.

Keller, C. K. Carbon Exports from Terrestrial Ecosystems: A Critical-Zone Framework. Ecosystems 22, 8,1691–705 (2019). https://doi.org/10.1007/s10021-019-00375-9.

Knyazikhin, Y., Martonchik, J. V., Myneni, R. B., Diner, D. J.,& Running, S. W. Synergistic Algorithm for Estimating Vegetation Canopy Leaf Area Index and Fraction of Absorbed Photosynthetically Active Radiation from MODIS and MISR Data. J. Geophys. Res.: Atmospheres 103, D24, 32257–75 (1998). https://doi.org/10.1029/98JD02462.

Knoben, W. J. M., Thébault, C., Keshavarz, K., Torres-Rojas, L., Chaney, N. W., Pietroniro, A., and Clark, M. P.: Catchment Attributes and MEteorology for Large-Sample SPATially distributed analysis (CAMELS-SPAT): streamflow observations, forcing data and geospatial data for hydrologic studies across North America, Hydrol. Earth Syst. Sci., 29, 5791–5833, (2025). https://doi.org/10.5194/hess-29-5791-2025.

Kratzert, F., Nearing, G., Addor, N., Erickson, T., Gauch, M., Gilon, O., Gudmundsson, L., Hassidim, A., Klotz, D., Nevo, S., Shaley, G., & Matias, Y. Caravan - A global community dataset for large-sample hydrology. Sci. Data 10, 61 (2023). https://doi.org/10.1038/s41597-023-01975-w.

Lehner, B., & Grill, G. Global River Hydrography and Network Routing: Baseline Data and New Approaches to Study the World’s Large River Systems. Hydrol. Process. 27, 15 2171–86 (2013). https://doi.org/10.1002/hyp.9740.

Lehner, B., Messager, M. L., Korver, M. C., &Linke, S. Global Hydro-Environmental Lake Characteristics at High Spatial Resolution.” Sci. Data 9, 1, 351 (2022). https://doi.org/10.1038/s41597-022-01425-z.

Lin, H. Earth’s Critical Zone and Hydropedology: Concepts, Characteristics, and Advances. Hydrol. Earth Syst. Sci. 14, 1, 25–45 (2010). https://doi.org/10.5194/hess-14-25-2010.

Linke, S. et al. Global Hydro-Environmental Sub-Basin and River Reach Characteristics at High Spatial Resolution. Sci. Data 6, 1, 283 (2019). https://doi.org/10.1038/s41597-019-0300-6.

McMillan, H. K. A review of hydrologic signatures and their applications. WIREs Water 8, 1, (2021). https://doi.org/10.1002/wat2.1499.

Muñoz-Villers, L. E. & McDonnell, J. J. Land Use Change Effects on Runoff Generation in a Humid Tropical Montane Cloud Forest Region. Hydrol. Earth Syst. Sci. 17, 9, 3543–60 (2013). https://doi.org/10.5194/hess-17-3543-2013.

Myneni, R., Knyazikhin, Y., & Park, T. MODIS/Terra Leaf Area Index/FPAR 8-Day L4 Global 500m SIN Grid V061 [Dataset]. NASA Land Processes Distributed Active Archive Center. (2021). https://doi.org/10.5067/MODIS/MOD15A2H.061 Date Accessed: 2026-09-25

Newman, A. J., et al. Development of a Large-Sample Watershed-Scale Hydrometeorological Data Set for the Contiguous USA: Data Set Characteristics and Assessment of Regional Variability in Hydrologic Model Performance. Hydrol. Earth Syst. Sci.19, 1, 209–23 (2015). https://doi.org/10.5194/hess-19-209-2015.

Omernik, J. M., & Griffith, G. E. Ecoregions of the conterminous United States: evolution of a hierarchical spatial framework. J. Environ. Mgmt. 54, 6, 1249-1266 (2014). https://doi.org/10.1007/s00267-014-0364-1.

Petersky, R., & Harpold, A. (2018). Now you see it, now you don't: a case study of ephemeral snowpacks and soil moisture response in the Great Basin, USA. Hydrology and Earth System Sciences, 22, 4891–4906.

Sterle, G., et al. CAMELS-Chem: Augmenting CAMELS (Catchment Attributes and Meteorology for Large-Sample Studies) with Atmospheric and Stream Water Chemistry Data. Hydrol. Earth Syst. Sci. 28, 3, 611–30 (2024). https://doi.org/10.5194/hess-28-611-2024.

Tague, C., et al. James Buttle Review: Dynamic Water Storage Shapes Critical Zone Function in Snow-Dominated Mountain Watersheds. Hydrol. Process. 39, 11, e70325 (2025). https://doi.org/10.1002/hyp.70325.,

Thornton, M. M., Shrestha, R., Wei, Y., Thornton, P. E., Kao, S-C. Daymet: Daily Surface Weather Data on a 1-Km Grid for North America, Version 4 R1. ORNL Distributed Active Archive Center (2022). htpps://doi.org/10.3334/ORNLDAAC/2129. Date Accessed: 2026-07-16

U.S. Geological Survey. gdptools (Version 0.3.11). U.S. Geological Survey Water Mission Area. [Software]. https://gdptools.readthedocs.io/en/latest/index.html

Vlah, M. J., et al. MacroSheds: A Synthesis of Long-Term Biogeochemical, Hydroclimatic, and Geospatial Data from Small Watershed Ecosystem Studies. Limnology and Oceanography Letters 8, 3, 419–52 (2023). https://doi.org/10.1002/lol2.10325.

Wlostowski, A. N., et al. “Signatures of Hydrologic Function Across the Critical Zone Observatory Network.” Water Resourc. Res. 57, 3, e2019WR026635 (2021). https://doi.org/10.1029/2019WR026635.

Zhang, M., et al. A Global Review on Hydrological Responses to Forest Change across Multiple Spatial Scales: Importance of Scale, Climate, Forest Type and Hydrological Regime. J. of Hydrol. 546, 44–59 (2017). https://doi.org/10.1016/j.jhydrol.2016.12.040.

Zimmer, M. A., & Gannon, J. P. Run-off Processes from Mountains to Foothills: The Role of Soil Stratigraphy and Structure in Influencing Run-off Characteristics across High to Low Relief Landscapes. Hydrol. Process. 32, 11, (2018). https://doi.org/10.1002/hyp.11488.

# HISSS Manuscript Draft (auto-synced snapshot)

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/d/e/2PACX-1vS7j4FRp7SEwlXoBUVA8NA7cj_I0XzyS0u58r3bl8SOz4BfpZPrdPJge4RMcFocnX8Gnllkc1M-CTJ3/pub)
> **Last synced**: 2026-08-24
>
> This is a read-only snapshot of the collaborative manuscript draft used for
> change detection and reconciliation review (see CLAUDE.md → Session-Start
> Workflow). Do NOT hand-edit the body below the header — it is overwritten at
> each sync. Discrepancies between the manuscript's methods and the code/docs
> are logged in CHANGELOG.md → `[Unreleased]` → `### Manuscript Reconciliation
> Log`; manuscript-side corrections must be made in the Google Doc by the
> co-authors (@ Arik convention per the draft's own notes).

---

Journal: Scientific Data - Special Collection on Water Storage due Nov 9, 2026

Draft Title: Harmonized dataset for a nalyzing flows in the critical zone

HISSS: Hydrologic Information Signatures and Summary Statistics

Alt: Snow to Flow: Dataset of harmonized metrics of watershed hydroclimate and stream flow

CZ-SYNC Critical Zone - SYNChronizing hydrologic signatures

CZ-SynCH: Critical Zone - Synchrony of Climate and Hydrology

HydroCZAR: Hydrologic Critical Zone Attribute Repository

Workflow Diagram & figure ideas

Shooting for August 14th paper being done (internally), submit by Nov. 9th (after some data analysis from other projects)

Ways to help:

- Fine line writing edits (especially where comments are located)
- Review to ask questions about the process or the datasets (summary meta doc)
- ** if there is a question for Arik @ him so he can look it up and fill it in, or if we need to update something then @ him and ask him to modify
- Add/edit Citations
- Support on making the workflow diagram nice, or example figures
- Sill diagram example

To Do:

- (PIs) Email Cuashi – ask if it is possible to host the dashboard that Arik is making [done]
- [done] (Galen) - asking around with USGS re: reviewing the paper (Roy Sando? Ask for mid-september)
- (Anyone) Add acknowledgement of AI coding and our review process (check the MS policy)
- Our novel approach - re group review of code, outputs and development of rules
- ReadME — Arik's Claude can handle this once dataset is finalized
- (Arik) move code into the CZ github – (purpose was more for people to check the way things are calculated) potential for it to get turned into a package?
- Check in on dataset, too, after reading through manuscript
- (Kendra) copy over authorship table ; Arik is first author on data release, kek as 1st on data paper
- (everyone) Limitations of data (please add to usage notes section of MS)
- Nic is compiling all the data files into one nice zip file - they have code to streamline this
- (everyone) Text in sections 3, 4, and 5

## 1 Background & Summary

The critical zone defines the region of the Earth's surface and subsurface where water is used by humans and the environment (Lin 2010). Water in the critical zone is stored in soil and vegetation and can be actively exchanged with the atmosphere above and groundwater below. The makeup of the critical zone determines how streamflow responds under normal conditions, disturbances (i.e. drought, fire, land conversion), and increasing atmospheric temperatures and CO2 (Condon et al. 2020, Keller et al. 2019). Open questions in Earth system science involve understanding how water availability in the critical zone is mediated by its characteristics (i.e. soil texture and depth, rooting depth, vegetation type, hydroclimate, etc.). Streamflow is an important metric of water availability that indicates the severity of stress for human and natural systems (Brooks et al. 2015, Wlostowski et al. 2021). However, understanding how streamflow responses to weather conditions are mediated by the critical is challenging because of the spatial and temporal heterogeneity of subsurface and aboveground critical zone components (Lin 2010, Brooks et al. 2015). Further, datasets describing critical zone features and climate are not readily available at spatial and temporal scales relevant for streamflow. The Hydrologic Information Signatures and Summary Statistics (HISSS) dataset presented in this manuscript overcomes this issue by harmonizing climate, vegetation, soil, and human impact data at the watershed scale for watersheds across the US and Canada.

While precipitation provides the primary water input to the critical zone, how that precipitation is stored drives where and when water is available (Tague et al. 2025). Water in the critical zone may be stored in soil, plants, and snow. Many studies have considered how the critical zone mediates the exchanges of water from these water stores. Soil depth and texture influences the partitioning of rain and snow in the critical zone (Hammon et al, 2019) and runoff generation (Zimmer and Gannon 2018, Chen et al. 2024). The role of vegetation in mediating runoff responses has been explored previously for different drought scenarios (Corak et al. 2024), with regards to the loss of global forests (Zhang et al. 2017), and for land conversion in tropical forests (Muñoz-Villers and McDonnell 2013). Snow depth, snow water equivalent, and timing of melt influence runoff and groundwater recharge in the critical zone (Wlostowski et al. 2021) However, while these studies explore relationships at individual watersheds or specific ecosystems, data are not available at continental scales that allow for cross comparisons and synthesis across climate and ecoregion.

One of the challenges of harmonizing climate and vegetation datasets with streamflow are mismatches in temporal and spatial scale (Vlah et al. 2023). Streamflow represents pointscale measurements that integrate water fluxes across a watershed, whereas climate and vegetation data are typically provided as gridded products representing average values across a defined pixel. Several large-sample datasets have addressed this mismatch by aggregating gridded and geospatial data to the watershed scale, with products focused on hydroclimate and hydrological signatures (Newman et al. 2015; Addor et al. 2017), biogeochemistry (Vlah et al. 2023), stream chemistry (Sterle et al. 2024), as well as broad multi-source archives that now cover much of North America (HYSETS, Arsenault et al. 2020; Caravan, Kratzert et al. 2023). Among the products that do report hydrological signatures, however, these are given only as single period-of-record values rather than as a time series, these products do not provide per-signature trend or changepoint statistics, and the harmonized multi-region archives typically derive their climate forcing from comparatively coarse reanalyses in most cases (ERA5-Land at ~11 km or coarser). While the recent CAMELS-SPAT (Knoben et al. 2025) spans the US and Canada and does include high spatial resolution climate data and time varying MODIS leaf area index (LAI) and land use and land cover (LULC), it is focused on a subset of available watersheds (n < 1500) with the goal of supporting hydrological model development and therefore precomputes only a narrow subset of hydroclimatological signatures as period-of-record statistics.

The dataset and software library presented here represent an attempt to provide holistic hydroclimate, vegetation, snow, and geophysical information at the scale of individual watersheds across the continental United States and Canada. For each of the ~7000 watersheds in the study, we calculate roughly 100 hydrological signatures and distribute the underlying annual values together with a consistent suite of computed metrics: central tendency statistics (mean, median), linear and robust (Theil-Sen) trend slopes, Spearman and Mann Kendall trend significance, Pettitt changepoint detection, and other supporting diagnostics. These are paired with annually resolved MODIS and NLCD LULC and monthly resolved MODIS LAI, watershed climate features derived from ~1 km Daymet, and HydroATLAS physiographic attributes. Because every signature is computed by an open, cross-language library (Julia, Python, and R), each metric is reproducible and extensible. These products are indented for critical-zone scientists studying catchment processes, snow hydrology, and vegetation dynamics, among others, and are structured to join readily with existing datasets of stream chemistry (e.g., CAMELS-Chem, biogeochemistry (e.g., MacroSheds), and vapor pressure deficit (e.g., Corack et al. 2025).

In this manuscript, we describe the HISSS data and code repository to explore how the critical zone mediates hydrologic responses at watershed scale. For over 6000 watersheds across the US and Canada, we curated data and summary statistics for key variables including streamflow, snow water equivalent (SWE), meteorological condition, vegetation state, subsurface characteristics, and human interventions. Watershed boundaries and streamflow for the US and Canada were obtained from the US Geological Survey (USGS) and the Water Survey of Canada, respectively (Falcone, 2011, Environment and Climate Change Canada, 2016). Basin characteristics including subsurface composition and human interventions are from HydroAtlas (Linke et al., 2013; Lehner et al., 2022). We provide both raw and processed data for each watershed. Daily meteorological conditions and SWE at 1 km were obtained from Daymet (Thornton et al. 2022). Land surface characteristics were derived from the MODIS MCD12Q1 annual land cover classification product and the 4-day MCD15A3H LAI product at 500 m (Friedl and Sulla-Menashe 2015, Myneni et al. 2015). The raw data are at the temporal resolution of the original data product. For the processed data, the hydrologic and meteorological data are available at a daily timestep, leaf area index is available at a monthly timestep, and basin characteristics are static fields. Annualized trends are provided for all non-static variables. In the HISSS repository, all data and code are freely available.

## 2 Methods

We assembled daily streamflow records for 16,994 gages across the United States (9,154; USGS National Water Information System) and Canada (7,840; HYDAT), of which 8,014 (6,160 usgs; 1,854 HYDAT) yielded usable daily records totalling more than 111 million observations. We delineated each gage's contributing watershed, area-normalized discharge, and applied an automated quality-control procedure to every record. From each standardized series we computed 106 hydrological signatures across 13 categories, and summarized every signature's annual time series with a common set of descriptive statistics, trend estimators, and changepoint tests. We then paired each watershed with harmonized physiographic attributes, time-varying satellite vegetation and land-cover records, and gridded meteorological forcing. The following subsections detail the specific methods and filtering procedures applied, and the code for the full workflow is available in an open library in Julia, Python, and R \<LINK TO REPO\>.

\<SCHEMATIC SHOWING PROCESSING WORKFLOW FROM KENDRA\>

### 2.1 Input data

\<TABLE OF ALL INPUT DATA SOURCES\>

#### 2.1.1 Watershed Boundaries

To define the spatial extents of gaged watersheds in this work, we combined authoritative basin polygons from the United States and Canada. For the former, we used watershed polygons from the USGS GAGES-II dataset (Falcone, 2011) and merged the individual shapefiles for reference and non-reference gages across all US regions into a single US-wide shapefile (covering 100% of US gages). For the latter, we used the Water Survey of Canada's Hydrometric Network Basin Polygons dataset (Environment and Climate Change Canada, 2016) and merged them into one Canada-wide layer (covering 98.3% of Canadian gages). Gages lacking an official polygon were delineated from the HydroBASINS level-12 drainage network (Lehner and Grill, 2013; detailed methods below) by aggregating all basins upstream of the gage's outlet. The geometry source is recorded for every watershed. We excluded 54 basins that exceeded 100,000 km^2, then transformed all remaining 7,964 watershed boundaries into a common projection and merged them into one full-domain shapefile. We used these boundaries to aggregate gridded climate data and remotely sensed LULC and LAI data. We also utilized HydroATLAS (Linke et al., 2019) to define static basin attributes within our study domain.

#### 2.1.2 Streamflow data

We applied quality assurance at three stages in the development of the annualized data: retrieval, compilation, and per-year qualification. We retrieved daily observed streamflow from 01 October 1979 to 30 September 2025 for all candidate gages (n = 16,994) using the USGS's dataRetrieval R package (DeCicco et al., 2026) for US gages and the Water Survey of Canada's tidyhydat R package (Albers et al., 2026) for Canadian gages. Gages that returned fewer than 20 years containing valid daily observations (8,980) were excluded. The remaining 8,014 gages (6,160 US; 1,854 Canadian) returned usable daily records that we compiled into a single dataset of 111.6 million observations keyed to the respective agency's gage ID. We retained original agency quality flags, including ice-affected and estimated values. We stored missing days as NA rather than zero-filled, treated zero flow as a valid observation, and area-normalized discharge (mm d-1) for all but 73 gages, as noted above (Sect. 2.1.1). Before computing signatures, we standardized each gage's record to a continuous daily grid and screened by WY: internal gaps of up to three consecutive days were filled by linear interpolation, and we set to NA any WY with more than 30 missing values or any single data gap exceeding three consecutive days. We further flagged any WY-season with less than 80% of raw observations passing quality assurance, setting the affected seasonal metrics to NA for that WY. Finally, we calculated signatures for gages only if they retained at least 20 qualifying WYs, spanned at least 60% of the possible WYs in the analysis window, and included 80% of the possible WYs in the first and final decade of the window. We calculated signatures (Sect. 2.2) over two standard analysis windows (WYs 1993-2025 and 1980-2025) yielding 6,678 gages for WYs 1993–2025 and 6,250 for WYs 1980–2025.

#### 2.1.3 Climate and snow data

Meteorological data were supplied by Daymet (Thornton et al. 2022), a research product developed by Oak Ridge National Laboratory, which provides 1 km x 1 km gridded, long-term, continuous estimates of daily meteorological conditions across North America. Daily Daymet estimates of minimum and maximum temperature, precipitation, vapor pressure, shortwave radiation, snow water equivalent, and day length covering 1980–2023 were downloaded. The gridded Daymet data were spatially aggregated to the 6,041 basins using the gdptools package², which performs area-weighted zonal aggregation of gridded data to polygon features. Specifically, the agg_gen function was used to compute the fractional overlap between each 1 km Daymet grid cell and each basin polygon, with grid cell values weighted proportionally by the fraction of the cell area falling within each basin boundary. This intersection and weighting was computed once and applied consistently across all variables and time steps. The result is a single area-weighted daily value for each variable for each basin, effectively summarizing the spatial distribution of each meteorological variable across the contributing area of each basin into a single representative daily estimate.

#### 2.1.4 Basin attributes and land cover

To characterize the physiographic characteristics of each watershed, we join basin attributes from HydroATLAS (BasinATLAS version 10; Linke et al. 2019), a global compilation of hydro-environmental descriptors built on the HydroBASINS level-12 subbasin framework (Lehner and Grill 2013). For every gage, we identify the intersecting level-12 subbasin and the complete set of upstream subbasins, then aggregate the HydroATLAS attributes over that contributing area. This yields approximately 210 attributes per watershed spanning six thematic categories: hydrology, physiography, climate, LULC, soils and geology, and anthropogenic influence. Continuous and percentage attributes are aggregated as area-weighted means, elevation as its spatial minimum and maximum, and categorical attributes as the area-weighted majority class. The resulting table is keyed to gage ID so that it joins directly to streamflow signatures.

The land surface characteristics in the HISSS dataset include land cover classifications from MODIS MCD12Q1 Land Cover Type Yearly product and LAI from the MODIS MCD15A3H version 6.1. Both have a native spatial resolution of 500 m and global coverage, with MCD12Q1 aggregated annually and MCD15A3H as a 4-day composite. MCD12Q1 classifies each pixel under eight schemes all of which HISSS retains: five Land Cover Type schemes based on the International Geosphere-Biosphere Programme (IGBP), University of Maryland (UMD), LAI/fPAR biome, Biome/Biogeochemical Cycles (BGC), and Plant Functional Types (PFT), and 3 FAO Land Cover Classification System property layers describing land cover, land use, and surface hydrology, respectively. Together these yield 102 per-class coverage fractions for each watershed and year.

The LAI product represents the seasonal phenology of vegetation on the land surface. It is estimated from a look-up table derived from a 3-D radiative transfer model relating surface reflectance to observed LAI values (Knyazikhin et al. 1998). It also relies on the MODIS Land Cover Type 3 to determine the most likely LAI value for a pixel given its vegetation type. MODIS LAI has valid values between 0 and 10 m2m-2 which describe the one-sided projection of leaf area per unit ground area.

### 2.2 Metrics for analysis

#### 2.2.1 Computing Summary Statistics

We calculated a comprehensive suite of 106 streamflow and hydroclimate signatures for each gage and its associated watershed, comprising 90 annually resolved metrics and 16 per-gage scalar diagnostics. The selection of metrics are consistent with analysis in catchment hydrology (McMillan, 2021) and snow hydrology (??Hatchett, 2021; Petersky and Harpold, 2018??) and spans 13 categories: flow volumes and percentiles, flow timing, flow duration curve slopes, baseflow, recession behavior, high-and low-flow pulses, flashiness, runoff ratios, streamflow elasticity, precipitation-streamflow (P-Q) seasonality, catchment storage, negative-flow days, and snow metrics. Catchment storage (avg_storage) is computed but omitted from major analyses because its simplified water balance does not account for evapotranspiration. Details for each metric, including the definition, requirements and relevant citations are in the metadata file (name).

We calculated baseflow indices (BFI) using both the Lyne-Hollick and the Eckhardt filters, each in two variants: a default parameterization (Eckhardt: BFImax = 0.8, a = 0.98; Lyne-Hollick: a = 0.925) and a recession-parameterized variant in which the filter constant is replaced by a site-specific recession constant, estimated as the median of event-based recession ratios under the assumption of a linear reservoir. In the parameterized Eckhardt filter only the recession constant is replaced (BFImax remains 0.8); the Lyne-Hollick parameterization is heuristic, as that filter's parameter has no formal derivation from recession analysis.

For the Eckhardt filter, we used both default parameters (BFImax = 0.8, a=0.98) and a parametrized version that includes recession parameters that were calculated for each site. We then calculated statistics for each metric across both the full period of record and a subset from 1993-2025 wateryears which represent the highest (proportion of sites) including slope, Spearman's Rho, p-value, mean and median.

#### 2.2.2 Computing Long-Term Signatures and Trends

For the trends analysis, we removed years with more than 3 consecutive days of NAs or more than 30 days of NAs throughout the year. We flagged (but did not remove) years with negative values or constant monthly standard deviations during periods of non-zero streamflow indicative of data errors. We also computed 12 automated data-quality flags (range and consistency checks) on the resulting signatures. A gage was included only if it had at least 20 qualifying water years and at least 60% of the water years in the period of interest. For a given signature, trend statistics were then computed only if its annual series had at least 20 valid values, at least 60% of the entire series had valid values, and at least 80% of both the first and last decades of the trend had valid values; the event-based recession and elasticity metrics are exempt from these completeness requirements. Any seasonal metric additionally required at least 80% of the period's observations had valid values. Seasonal periods are defined as: Winter: December-February, Spring: March-May, Summer: June-August, Fall: September-November, Annual metrics were calculated on the water year October 1 - September 30.

## 3 Data Record

Point to metadata and DOI and README

## 4 Data Overview

Figures summarizing data (map of streamflow, data by ecoregions, CDF of years and sites of data, Whittaker plot)

Rough draft of figure summarizing sites and signature categories, needs to be cleaned up

## 5 Usage Notes

- How to join with other datasets (CAMELS-Chem, MacroShed, Daymet VPD)

## References

- Addor, Nans, Andrew J. Newman, Naoki Mizukami, and Martyn P. Clark. "The CAMELS Data Set: Catchment Attributes and Meteorology for Large-Sample Studies." Hydrology and Earth System Sciences 21, no. 10 (2017): 5293–313. https://doi.org/10.5194/hess-21-5293-2017.
- Brooks, Paul D., Jon Chorover, Ying Fan, et al. "Hydrological Partitioning in the Critical Zone: Recent Advances and Opportunities for Developing Transferable Understanding of Water Cycle Dynamics." Water Resources Research 51, no. 9 (2015): 6973–87. https://doi.org/10.1002/2015WR017039.
- Chen, Hang, Qifei Niu, James P. McNamara, and Alejandro N. Flores. "Influence of Subsurface Critical Zone Structure on Hydrological Partitioning in Mountainous Headwater Catchments." Geophysical Research Letters 51, no. 6 (2024): e2023GL106964. https://doi.org/10.1029/2023GL106964.
- Condon, Laura E., Katherine H. Markovich, Christa A. Kelleher, Jeffrey J. McDonnell, Grant Ferguson, and Jennifer C. McIntosh. "Where Is the Bottom of a Watershed?" Water Resources Research 56, no. 3 (2020): e2019WR026010. https://doi.org/10.1029/2019WR026010.
- Corak, Nicholas K., Jason A. Otkin, Trent W. Ford, and Lauren E. L. Lowman. "Unraveling Phenological and Stomatal Responses to Flash Drought and Implications for Water and Carbon Budgets." Hydrology and Earth System Sciences 28, no. 8 (2024): 1827–51. https://doi.org/10.5194/hess-28-1827-2024.
- Corak, Nicholas K., Peter E. Thornton, and Lauren E. L. Lowman. "A High Resolution, Gridded Product for Vapor Pressure Deficit Using Daymet." Scientific Data 12, no. 1 (2025): 256. https://doi.org/10.1038/s41597-025-04544-5.
- Falcone, J., 2011, GAGES-II: Geospatial Attributes of Gages for Evaluating Streamflow: U.S. Geological Survey data release, https://doi.org/10.5066/P96CPHOT.
- Falcone, J.A., 2017. US Geological Survey GAGES-II time series data from consistent sources of land use, water use, agriculture, timber activities, dam removals, and other historical anthropogenic influences. US Geological Survey (USGS) Data Release, p.740. https://doi.org/10.5066/F7HQ3XS4
- Friedl, M., Sulla-Menashe, D., 2015. Mcd12q1 modis. Terra+ aqua land cover type yearly l3 global 500m SIN grid 6.
- Hammond, John C., Adrian A. Harpold, Sydney Weiss, and Stephanie K. Kampf. "Partitioning Snowmelt and Rainfall in the Critical Zone: Effects of Climate Type and Soil Properties." Hydrology and Earth System Sciences 23, no. 9 (2019): 3553–70. https://doi.org/10.5194/hess-23-3553-2019.
- Keller, C. Kent. "Carbon Exports from Terrestrial Ecosystems: A Critical-Zone Framework." Ecosystems 22, no. 8 (2019): 1691–705. https://doi.org/10.1007/s10021-019-00375-9.
- Knyazikhin, Y., J. V. Martonchik, R. B. Myneni, D. J. Diner, and S. W. Running. "Synergistic Algorithm for Estimating Vegetation Canopy Leaf Area Index and Fraction of Absorbed Photosynthetically Active Radiation from MODIS and MISR Data." Journal of Geophysical Research: Atmospheres 103, no. D24 (1998): 32257–75. https://doi.org/10.1029/98JD02462.
- Lehner, Bernhard, and Günther Grill. "Global River Hydrography and Network Routing: Baseline Data and New Approaches to Study the World's Large River Systems." Hydrological Processes 27, no. 15 (2013): 2171–86. https://doi.org/10.1002/hyp.9740.
- Lehner, Bernhard, Mathis L. Messager, Maartje C. Korver, and Simon Linke. "Global Hydro-Environmental Lake Characteristics at High Spatial Resolution." Scientific Data 9, no. 1 (2022): 351. https://doi.org/10.1038/s41597-022-01425-z.
- Lin, H. "Earth's Critical Zone and Hydropedology: Concepts, Characteristics, and Advances." Hydrology and Earth System Sciences 14, no. 1 (2010): 25–45. https://doi.org/10.5194/hess-14-25-2010.
- Linke, Simon, Bernhard Lehner, Camille Ouellet Dallaire, et al. "Global Hydro-Environmental Sub-Basin and River Reach Characteristics at High Spatial Resolution." Scientific Data 6, no. 1 (2019): 283. https://doi.org/10.1038/s41597-019-0300-6.
- Muñoz-Villers, L. E., and J. J. McDonnell. "Land Use Change Effects on Runoff Generation in a Humid Tropical Montane Cloud Forest Region." Hydrology and Earth System Sciences 17, no. 9 (2013): 3543–60. https://doi.org/10.5194/hess-17-3543-2013.
- Myneni, R., Knyazikhin, Y., Park, T., 2015. Mod15a2h modis/terra leaf area index/fpar 8-day l4 global 500m sin grid v006. NASA EOSDIS Land Processes DAAC.
- Newman, A. J., M. P. Clark, K. Sampson, et al. "Development of a Large-Sample Watershed-Scale Hydrometeorological Data Set for the Contiguous USA: Data Set Characteristics and Assessment of Regional Variability in Hydrologic Model Performance." Hydrology and Earth System Sciences 19, no. 1 (2015): 209–23. https://doi.org/10.5194/hess-19-209-2015.
- Sterle, Gary, Julia Perdrial, Dustin W. Kincaid, et al. "CAMELS-Chem: Augmenting CAMELS (Catchment Attributes and Meteorology for Large-Sample Studies) with Atmospheric and Stream Water Chemistry Data." Hydrology and Earth System Sciences 28, no. 3 (2024): 611–30. https://doi.org/10.5194/hess-28-611-2024.
- Tague, Christina, Holly R. Barnard, Adrian A. Harpold, et al. "James Buttle Review: Dynamic Water Storage Shapes Critical Zone Function in Snow-Dominated Mountain Watersheds." Hydrological Processes 39, no. 11 (2025): e70325. https://doi.org/10.1002/hyp.70325.
- Thornton, M.M., R. Shrestha, Y. Wei, P.E. Thornton, and S-C. Kao. "Daymet: Daily Surface Weather Data on a 1-Km Grid for North America, Version 4 R1." ORNL Distributed Active Archive Center, January 1, 2022. doi:10.3334/ORNLDAAC/2129. Date Accessed: 2026-07-16
- Vlah, Michael J., Spencer Rhea, Emily S. Bernhardt, et al. "MacroSheds: A Synthesis of Long-Term Biogeochemical, Hydroclimatic, and Geospatial Data from Small Watershed Ecosystem Studies." Limnology and Oceanography Letters 8, no. 3 (2023): 419–52. https://doi.org/10.1002/lol2.10325.
- Wlostowski, Adam N., Noah Molotch, Suzanne P. Anderson, et al. "Signatures of Hydrologic Function Across the Critical Zone Observatory Network." Water Resources Research 57, no. 3 (2021): e2019WR026635. https://doi.org/10.1029/2019WR026635.
- Zhang, Mingfang, Ning Liu, Richard Harper, et al. "A Global Review on Hydrological Responses to Forest Change across Multiple Spatial Scales: Importance of Scale, Climate, Forest Type and Hydrological Regime." Journal of Hydrology 546 (March 2017): 44–59. https://doi.org/10.1016/j.jhydrol.2016.12.040.
- Zimmer, Margaret A., and John P. Gannon. "Run-off Processes from Mountains to Foothills: The Role of Soil Stratigraphy and Structure in Influencing Run-off Characteristics across High to Low Relief Landscapes." Hydrological Processes 32, no. 11 (2018): 1546–60. https://doi.org/10.1002/hyp.11488.

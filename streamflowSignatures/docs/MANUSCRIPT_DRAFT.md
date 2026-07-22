# HISSS Manuscript Draft (auto-synced snapshot)

> **Auto-synced from**: [Google Doc](https://docs.google.com/document/d/e/2PACX-1vS7j4FRp7SEwlXoBUVA8NA7cj_I0XzyS0u58r3bl8SOz4BfpZPrdPJge4RMcFocnX8Gnllkc1M-CTJ3/pub)
> **Last synced**: 2026-07-21
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

Draft Title: Harmonized dataset for analyzing flows in the critical zone

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

- (PIs) Email Cuashi – ask if it is possible to host the dashboard that Arik is making
- (Galen) - asking around with USGS re: reviewing the paper (Roy Sando? Ask for mid-september)
- (Anyone) Add acknowledgement of AI coding and our review process (check the MS policy)
- Our novel approach - re group review of code, outputs and development of rules
- ReadME?
- (Arik) move code into the CZ github – (purpose was more for people to check the way things are calculated) potential for it to get turned into a package?
- (Kendra) copy over authorship table ; Arik is first author on data release, kek as 1st on data paper
- (everyone) Limitations of data (please add to usage notes section of MS)
- Nic is compiling all the data files into one nice zip file - they have code to streamline this

## 1 Background & Summary

The critical zone defines the region of the Earth's surface and subsurface where water is used by humans and the environment (Lin 2010). Water in the critical zone is stored in soil and vegetation and can be actively exchanged with the atmosphere above and groundwater below. The makeup of the critical zone determines how streamflow responds under normal conditions, disturbances (i.e. drought, fire, land conversion), and increasing atmospheric temperatures and CO2 (Condon et al. 2020, Keller et al. 2019). Open questions in Earth system science involve understanding how water availability in the critical zone is mediated by its characteristics (i.e. soil texture and depth, rooting depth, vegetation type, hydroclimate, etc.). Streamflow is an important metric of water availability that indicates the severity of stress for human and natural systems (Brooks et al. 2015, Wlostowski et al. 2021). However, understanding how streamflow responses to weather conditions are mediated by the critical is challenging because of the spatial and temporal heterogeneity of subsurface and aboveground critical zone components (Lin 2010, Brooks et al. 2015). Further, datasets describing critical zone features and climate are not readily available at spatial and temporal scales relevant for streamflow. The Hydrologic Information Signatures and Summary Statistics (HISSS) dataset presented in this manuscript overcomes this issue by harmonizing climate, vegetation, soil, and human impact data at the watershed scale for watersheds across the US and Canada.

While precipitation provides the primary water input to the critical zone, how that precipitation is stored drives where and when water is available (Tague et al. 2025). Water in the critical zone may be stored in soil, plants, and snow. Many studies have considered how the critical zone mediates the exchanges of water from these water stores. Soil depth and texture influences the partitioning of rain and snow in the critical zone (Hammon et al, 2019) and runoff generation (Zimmer and Gannon 2018, Chen et al. 2024). The role of vegetation in mediating runoff responses has been explored previously for different drought scenarios (Corak et al. 2024), with regards to the loss of global forests (Zhang et al. 2017), and for land conversion in tropical forests (Muñoz-Villers and McDonnell 2013). Snow depth, snow water equivalent, and timing of melt influence runoff and groundwater recharge in the critical zone (Wlostowski et al. 2021) However, while these studies explore relationships at individual watersheds or specific ecosystems, data are not available at continental scales that allow for cross comparisons and synthesis across climate and ecoregion.

One of the challenges of harmonizing climate and vegetation datasets with streamflow are mismatches in spatial scale (Vlah et al. 2023). While streamflow represents pointscale measurements that integrate water transfers across a watershed, climate and vegetation data are typically provided as gridded products representing average values across a defined pixel. Other datasets have been generated to overcome this issue to provide harmonized data at the watershed scale hydroclimate (Newman et al. 2015, Addor et al. 2017), for anthropogenic influences (Falcone 2017), and biogeochemical processes (Vlah et a. 2023), and stream chemistry (Sterle et al. 2024). However, while prior datasets provide harmonized data at watershed scale, they do not include summary statistics and tend to not include land cover and leaf area index (LAI). Additionally, hydroclimate data tends to be derived from coarse reanalysis products. The dataset and code presented in this manuscript are comprehensive of hydroclimate, vegetation, and subsurface at the HUC-8 scale, and provide computed metrics and trend analyses relevant for critical zone scientists interested in catchment processes, snow hydrology, and vegetation dynamics, among others. The dataset and library are flexible and may easily be joined with existing datasets of stream chemistry (i.e., CAMELS-Chem) and biogeochemistry (i.e., MacroSheds) and vapor pressure deficit (Corak et al. 2025) to extend its utility to questions regarding water quality.

In this manuscript, we describe the HISSS data and code repository to explore how the critical zone mediates hydrologic responses at watershed scale. For over 6000 watersheds across the US and Canada, we curated data and summary statistics for key variables including streamflow, snow water equivalent (SWE), meteorological condition, vegetation state, subsurface characteristics, and human interventions. Watershed boundaries and streamflow for the US and Canada were obtained from the US Geological Survey (USGS) and the Water Survey of Canada, respectively (Falcone, 2011, Environment and Climate Change Canada, 2016). Basin characteristics including subsurface composition and human interventions are from HydroAtlas (Linke et al., 2013; Lehner et al., 2022). We provide both raw and processed data for each watershed. Daily meteorological conditions and SWE at 1 km were obtained from Daymet (Thornton et al. 2022). Land surface characteristics were derived from the MODIS MCD12Q1 annual land cover classification product and the 4-day MCD15A3H LAI product at 500 m (Friedl and Sulla-Menashe 2015, Myneni et al. 2015). The raw data are at the temporal resolution of the original data product. For the processed data, the hydrologic and meteorological data are available at a daily timestep, leaf area index is available at a monthly timestep, and basin characteristics are static fields. Annualized trends are provided for all non-static variables. In the HISSS repository, all data and code are freely available.

## 2 Methods

\<SUMMARY PARAGRAPH OF METHODS\>

\<SCHEMATIC SHOWING PROCESSING WORKFLOW FROM KENDRA\>

### 2.1 Input data

\<TABLE OF ALL INPUT DATA SOURCES\>

#### 2.1.1 Watershed Boundaries

To define the spatial extents of gaged watersheds in this work, we used several data sources from the US and Canada. For the former, we accessed watershed polygons from the USGS GAGES II (Falcone, 2011) dataset and merged the individual shapefiles for reference and non-reference gages across all US regions into a single US-wide shapefile. For the latter, we downloaded per-region gaged watershed boundary shapefiles from Water Survey of Canada's Hydrometric Network Basin Polygons dataset (Environment and Climate Change Canada, 2016) and merged them into one Canada-wide shapefile. We transformed the US and Canada shapefiles into a common projection and merged them into one full-domain shapefile that we used to aggregate the gridded climate data into watershed-scale summaries. We also utilized the HydroBASINS (Lehner and Grill, 2013) and HydroATLAS (Linke et al., 2013; Lehner et al., 2022) products from HydroSHEDS to define static basin attributes within our study domain.

#### 2.1.2 Streamflow data

We downloaded daily observed streamflow from 1980 to present for each gaged watershed in the merged US-Canada shapefile using the USGS's dataRetrieval R package (DeCicco et al., 2026) for US gages and the Water Survey of Canada's tidyhydat R package (Albers et al., 2026) for Canadian gages. We compiled the per-gage records into a single dataset with all measured observations for every US and Canadian gage, coded by the respective agency's gage ID. As part of the aggregation process, we converted all observations to a common unit (mm d-1) where we normalized volumetric measurements by watershed area as defined in the watershed boundary shapefile (Sect. 2.1.1). We computed streamflow and derived metrics (Sect. 2.2) for the full period of record from each gage, even when they vary in record length.

Streamflow QC: record length and missing data considerations

We created a series of requirements for handling missing streamflow values, (i) NA's are retained, and (ii) we interpolate streamflow if up to 3 days of continuous streamflow data is missing.

#### 2.1.3 Climate and snow data

Meteorological data were supplied by Daymet (Thornton et al. 2022), a research product developed by Oak Ridge National Laboratory, which provides 1 km x 1 km gridded, long-term, continuous estimates of daily meteorological conditions across North America. Daily Daymet estimates of minimum and maximum temperature, precipitation, vapor pressure, shortwave radiation, snow water equivalent, and day length covering 1980–2023 were downloaded. The gridded Daymet data were spatially aggregated to the 6,041 basins using the gdptools package², which performs area-weighted zonal aggregation of gridded data to polygon features. Specifically, the agg_gen function was used to compute the fractional overlap between each 1 km Daymet grid cell and each basin polygon, with grid cell values weighted proportionally by the fraction of the cell area falling within each basin boundary. This intersection and weighting was computed once and applied consistently across all variables and time steps. The result is a single area-weighted daily value for each variable for each basin, effectively summarizing the spatial distribution of each meteorological variable across the contributing area of each basin into a single representative daily estimate.

#### 2.1.4 Basin attributes and land cover

The land surface characteristics in the HISSS dataset include land cover classifications from MODIS MCD12Q1 Land Cover Type Yearly product and LAI from the MODIS MCD15A3H version 6.1. The native spatial resolution of both datasets is 500 m with global spatial coverage. The MODIS land cover data provide five distinct classification schemes for the land surface based on: 1) Annual International Geosphere-Biosphere Programme (IGBP), 2) Annual University of Maryland (UMD), 3) Annual LAI classification, 4) Annual BIOME-Biogeochemical Cycles (BGC), and 5) Annual Plant Functional Types. Each classification system categorizes the land surface using its own distinct classes. The LAI product represents the seasonal phenology of vegetation on the land surface. It is estimated by using a look-up table that simulates a 3-D radiative transfer model relating surface reflectance to observed LAI values (Knyazikhin et al. 1998). It also relies on the MODIS Land Cover Type 3 to determine the most likely LAI value for a pixel given its vegetation type. MODIS LAI has valid values between 0 and 10 m2m-2 which describe the one-sided projection of leaf area over 1 sq meter of ground.

### 2.2 Metrics for analysis

#### 2.2.1 Basin, Gage, and Climate Data Filtering

#### 2.2.2 Computing Annual Metrics

We calculated a comprehensive suite of streamflow and climatic metrics for each gage and its associated watershed (n=xx). The selection of metrics are consistent with analysis in catchment hydrology (McMillan), snow hydrology (cites) and include flow volume, timing, flow duration curve, baseflow, recession, baseflow, pulse, flashiness, elasticity metrics. Details for each metric, including the definition, requirements and relevant citations are in the metadata file (name). We calculated baseflow indices (BFI) using both the Lyne-Hollick and the Eckhardt filter, with two sets of parameters for the Eckhardt filter. For the Eckhardt filter, we used both default parameters (BFImax = 0.8, a=0.98) and a parametrized version that includes recession parameters that were calculated for each site. We then calculated statistics for each metric across both the full period of record and a subset from 1993-2025 wateryears which represent the highest (proportion of sites) including slope, Spearman's Rho, p-value, mean and median.

#### 2.2.3 Computing Long-Term Signatures and Trends

For the trends analysis, we removed years with (i) more than 3 consecutive days of NAs, (ii) gaps of more than 30 days in a year, (iii) negative values, (iv) constant monthly standard deviations during periods of non-zero streamflow indicative of data errors. Items i, iii, and iv are set in the config to flag (not remove) by default (cite USGS trend data release). Requirements For a trend to be calculated included (i) at least 60% of the entire annual streamflow metric time series was complete, (ii) at least 80% of the first and last decade of the trend period must be complete. For any seasonal or annual metrics the period must have at least 80% of the data. Seasonal periods are defined as: Winter: December-February, Spring: March-May, Summer: June-August, Fall: September-November, Annual metrics were calculated on the water year October 1 - September 30.

## 3 Data Record

Point to metadata and DOI and README

## 4 Data Overview

Figures summarizing data (map of streamflow, CDF of years and sites of data, Whittaker plot)

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

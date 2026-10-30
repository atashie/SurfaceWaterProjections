# Consolidated relay note for the co-authors — manuscript + guidelines doc (2026-10-30)

STATUS: reviewed and complete 2026-10-30 (blind reviews A, B, C and Codex applied,
`docs/plans/2026-10-30-manuscript-reviews/`) but NOT SENT — DECISION (user, 2026-10-30): the
manuscript and the guidelines doc are updated together with the user, section by section, in the
terminal; the deposited readmes and dictionary (section G) are left as deposited. Every number below was measured on the HydroShare files as deposited
(`docs/plans/2026-10-30-manuscript-claim-ledger.md`, ids in brackets refer to its rows).
Format: one line per item, "replace X with Y" (user preference 2026-09-04); the three tables
are given as corrected rows (user decision 2026-10-30). This note folds in the open items of
`docs/plans/2026-10-12-coauthor-relay.md` and the 2026-10-28 entries of
`docs/reconciliation/manuscript_log.md`; it supersedes both.

## A. Manuscript — numbers

1. §1: "roughly 100 hydrological signatures" → "139 hydrological and climate signatures" [B3].
2. §2 preamble: "121 signatures across nine signature families (…): 100 annually resolved signatures, each summarized with a common set of 16 statistics, plus 21 per-gage diagnostics" → "160 signatures across nine signature categories (…): 139 annually resolved signatures, each summarized with 16 statistics, plus 21 per-gage diagnostics (2,277 columns)" [C5]; §2.2.1 already says 139 — keep both consistent.
3. §2 preamble: "of which 8,014 (6,160 USGS; 1,854 HYDAT) yielded usable daily records for watersheds smaller than 85,000 km2" → "of which 8,014 (…) yielded usable daily records (at least 20 calendar years with at least 30 days of flow above 0.0001 mm/day)"; 58 of the 8,014 drain more than 85,000 km² [C3]; "totalling more than 111 million observations" stands (111,624,189 daily records, of which 111,229,979 carry a discharge value) [C4].
4. Abstract: "6,041-7,964 locations depending on the subset" — 6,041 matches no deposited file; the deposit supports "6,250–8,014 locations" (trend products 6,250 / 6,678 → every compiled gage 8,014; polygons 7,964, 8,017 at full resolution) [A2]. Authors' call.
5. §2.1.1: "This resulted in 7,964 gages across the United States and Canada" → "This resulted in 8,014 gages. The simplified boundary layer (Sect. 2.1.2) holds 7,964 polygons of 100,000 km² or less, covering 7,960 of them plus four stations that were not compiled"; and "polygons smaller than 85,000 km2 (≤ largest HUC8 size)" → "100,000 km²" (6 polygons lie between 85,000 and 100,000 km²) [D4].
6. §2.1.1: "perennial defined as >30 days" → "at least 30 days" (the rule is ≥ 30 days per calendar year, ≥ 20 such years); please confirm that reference 23 (GAGES-II) is the intended source of this rule [D3, Q4].
7. §2.1.1 "from 01 October 1979 to 30 September 2025" and §3 Resource 3 "1980 through 2025" → "water years 1980–2025 (1 October 1979 – 1 October 2025; the last day belongs to WY 2026 and is excluded downstream)" [D1].
8. §2.1.2: "covering 98.3% of Canadian gages" — 1,823 of 1,853 with a boundary = 98.4 % (98.3 % counts the one mis-delineated station); either is fine [E1].
9. §2.1.3: "covering 1980–2023" → "covering 1980–2025" [F2]; delete "and day length" (six variables are deposited: precipitation, minimum and maximum temperature, SWE, vapour pressure, shortwave radiation) [F1].
10. §2.1.4: "mean area of 110 km2" — not verifiable from the deposit (it describes the global level-12 tier); cite its source [G1]; "approximately 210 attributes" (also Table 3 "210", §3 Resource 4 "211 watershed-scale attributes") → "198 attributes (211 columns including identifiers and diagnostics)" [G3].
11. §3 Resources 1 and 2: "1,653 columns comprising gage metadata, the 100 signatures x 16 statistics each, 21 signatures and related outputs that do not carry statistics, and 12 automated quality-assurance flags" → "2,277 columns comprising 20 gage-metadata columns, the 139 signatures × 16 statistics, 21 per-gage scalars and 12 automated quality-assurance flags" [J6]; "(18,898,406 and 24,366,487 rows for the two windows, respectively)" → "(27,766,578 and 34,813,881 rows)" [J7]; "hiss_signatures_wy{window}.csv" → "hisss_signatures_wy{window}.csv" [J9].
12. §3 Resource 3: "is provided for 6,087 basins (97,757,220 records, calendar years 1980-2023)" → "is provided for 8,017 basins (134,605,430 records, calendar years 1980–2025, on Daymet's 365-day calendar)"; "the variables precipitation, temperature, snow water equivalent (SWE), vapor pressure, and shortwave radiation" → "precipitation, minimum and maximum temperature, snow water equivalent, vapour pressure and shortwave radiation" [J11].
13. §6.1: counts (73 of 8,014; 32 and 28) are confirmed [M1].

## B. Manuscript — methods

14. §2.1.1: "gaps of up to three consecutive days were filled by linear interpolation, any WY with more than 30 missing values were set to NA, and if any data gap exceeded three consecutive days those days were set to NA" → "gaps of up to three consecutive days inside the water year were filled by linear interpolation; any water year with more than 30 missing days (counted before interpolation), with a gap longer than three consecutive days, or with a gap at its first or last day was excluded" [D7, D8, D9].
15. §2.1.1: "Seasonal periods are defined as: Winter: January-March, Spring: April-June, Summer: July - September, Fall: October-December" → "Winter: December–February, Spring: March–May, Summer: June–August, Fall: September–November (fall = October–November of the water year's start plus September of its end)" [D10].
16. §2.1.1: after "including ice-affected and estimated values" add: "A USGS daily value is retained whenever its agency status is approved (A) or provisional (P) and a finite value other than the no-data sentinel (−999999) is served; remark codes (R revised, e estimated, < and > censored at a reporting limit, [4], Ice) are carried as flags and never mask a value. An earlier compile had masked approved values that carried a remark code; the 193,542 affected days on 987 gages were re-retrieved from NWIS on 26 October 2026 and restored before the signatures were computed (the ledger is in Resource 3)" [D5].
17. §2.1.2: "We excluded basins that exceeded 85,000 km2, transformed …" → "The simplified boundary layer excludes the 53 basins larger than 100,000 km²; the full-resolution layer (Resource 4) retains them. We transformed …" [E3].
18. §2.1.2: "The final watershed boundaries were used to aggregate climate and landcover data only where we had sufficient streamflow data coverage for trends analyses (see Section 2.1.3)" → "Climate was aggregated over all 8,017 full-resolution polygons and land cover over all 7,964 simplified polygons (Sects. 2.1.3–2.1.4)" [E4].
19. §2.1.3: replace "The watershed size threshold corresponds to the largest HUC8 basin … basin-averaged climate inputs. We performed this aggregation using the gdptools package36, which performs area-weighted zonal aggregation of gridded data to polygon features across all variables and time steps. The result is a single area-weighted daily value for each variable for each basin." with "Each variable was aggregated to each basin as the mean of the 1 km cells weighted by the fraction of each cell inside the full-resolution basin polygon times the cell's true area, renormalized over valid cells (the HISSS Daymet toolchain in the code repository), giving one daily value per variable for 8,017 basins — every compiled gage with a boundary plus the 53 basins larger than 100,000 km² — for calendar years 1980–2025 on Daymet's 365-day calendar (31 December of leap years is filled by interpolation downstream). On the 5,965 basins of the original 6,087-basin gdptools36 product used for validation, the two methods agree to within 0.02 mm/day of precipitation and 0.003 °C of temperature." [F3, F4, F5].
20. §2.1.4: "Continuous and percentage attributes were aggregated as area-weighted means, elevation as the spatial minimum and maximum, and categorical attributes as the area-weighted majority" → "Upstream-accumulated and pour-point attributes (91 columns) are taken from the gage's outlet subbasin, which HydroATLAS has already accumulated over the same upstream area; subbasin-scale attributes are re-aggregated over the member subbasins as area-weighted means (92 columns), elevation minimum and maximum, and dominant classes (11 columns)" [G2].
21. §2.2.1: "Water year and seasonal metrics required having at least 80% of raw observations passing quality assurance, if this was not met, the affected metrics were set to NA for that WY or season" → "The eight seasonal metrics (the four seasonal flow totals and the four seasonal runoff ratios) required at least 80 % of the season's raw observations; otherwise the seasonal value was set to NA. Water years themselves were screened by the 30-day and 3-day rules of Sect. 2.1.1" [H5].
22. §2.2.1: "and spanned at least 60% of the possible WYs in the analysis window" → "and whose qualifying water years were at least 60 % of the water years from the window start to the gage's last water year with data within the window" [H4]; add "The gage sets were fixed before the remark-code days of Sect. 2.1.1 were restored; 31 (WY 1993–2025) and 29 (WY 1980–2025) further gages qualify only with those days and are listed in Resource 3 rather than added" [H4].
23. §2.2.2: "We flagged (but did not remove) water years that contained negative streamflow values or any calendar month in which at least 15 days of non-zero streamflow held at a single constant value, patterns that are indicative of data errors or highly regulated discharge" → "Water years containing negative values were retained, and the number of negative-flow days per water year is reported as a signature (negative_ann); calendar months in which at least 15 days of non-zero streamflow held a single constant value — a pattern indicative of data errors or highly regulated discharge — are recorded in the processing diagnostics" (the constant-value flag is not a column of the published tables) [I2].
24. §2.2.2: add "Pettitt changepoint fields are computed over WY 1980–2024 and require at least 20 valid values with at least 10 on each side of the candidate changepoint; the 20-year minimum is counted within the analysis window" [I4].
24a. §2.2.2: "trend statistics were computed only if its annual series had at least 20 valid values, at least 60% of the entire series had valid values, and at least 80% of both the first and last decades of the trend had valid values" → "All 16 statistics of a signature require at least 20 valid annual values (below that every statistic is NA); the six trend statistics additionally require valid values in at least 60 % of the years between the signature's first and last valid year and in at least 80 % of both its first and last decade"; keep the recession / elasticity exemption [I5].
24b. §2.1.3 (end) add: "Each gage's Daymet series was joined on Daymet's own calendar within the gage's record span, so a day missing from the streamflow record keeps its climate values; precipitation, snow water equivalent and temperature were screened per water year with the same 30-day and 3-day rules as streamflow" [I6].
24c. §2.1.4: "MCD12Q1 supplying consistent coverage across all HISSS watersheds" → "across the 7,964 watershed polygons (the 53 basins larger than 100,000 km² carry no land-cover record, so 44 and 45 gages of the two signature products lack MODIS and NLCD coverage)" [G7].
24d. §2.2.1: "the date range with the most gages with the same data record length" cannot be verified as an optimum (3,267 of the 6,678 gages carry all 33 years) → "a shorter common window that retains more qualifying gages" [H6].
25. §6.1: "(runoff ratios, elasticity, Q-P seasonality, storage) are structurally NA" → "(runoff ratios, elasticity, Q–P seasonality, average storage) are structurally NA; the baseflow and recession signatures are computed in native units" [M2].

## C. Tables — corrected rows

26. Table 1 (replace the body; citations as in §G below):

| Dataset type | Spatial unit | Temporal resolution | Data sources | Number of spatial units | Temporal coverage |
|---|---|---|---|---|---|
| Streamflow time series | Agency gage ID | Daily | NWIS20 | 6,160 | 1 Oct 1979 – 1 Oct 2025 |
| | | | HYDAT21 | 1,854 | 1 Oct 1979 – 31 Dec 2024 |
| Catchment polygons | Catchment polygon | – | GAGES-II23,33; WSC24; HydroBASINS34 | 7,964 (8,017 at full resolution) | – |
| Static basin characteristics | Catchment polygon | – | HydroATLAS25 | 8,014 | – |
| Basin-averaged climate | Catchment polygon | Daily | Daymet26 | 8,017 | 1980–2025 |
| Land cover and vegetation | Catchment polygon | Annual / monthly | MODIS27,28, NLCD29,30 | 7,964 (MODIS); 6,119 (NLCD) | 2001–2024 (MCD12Q1); Jul 2002 – Dec 2024 (MCD15A3H); 1985–2025 (NLCD) |
| Annual signatures and summary statistics | Gage ID | Annual | NWIS, HYDAT, Daymet | 6,678 / 6,250 | WY 1993–2025 / WY 1980–2025 |
| Trends* | Gage ID | Multi-decadal | as above | 6,678 / 6,250 | WY 1993–2025 / WY 1980–2025 |

Footnote * → "Trend statistics require at least 20 valid annual values, 60 % of the series and 80 % of both the first and last decade; the two gage sets qualify independently and neither is a subset of the other." Delete the current "6,041 – 7,964*" row and both footnotes [T2, T4, T5]; fix "Cathment", "Caculated".

27. Table 3 (replace three cells and two timespans):

| Data source | Characteristics | Details | Timespan | Resolution | Coverage | Spatial resolution |
|---|---|---|---|---|---|---|
| HydroATLAS25 | 198 attributes per watershed (211 columns) | 6 categories (unchanged) | static | – | Global | basin polygon |
| MODIS MCD12Q128 | Land cover type | 8 classification schemes: the 5 Land Cover Type layers (IGBP, UMD, LAI, BGC, PFT) and the 3 FAO-LCCS property layers, 102 class columns. Limitations: the provider flags reduced confidence from 2021 onward | 2001–2024 | Annual | Global | 500 m |
| MODIS MCD15A3H27 | Leaf area index | (unchanged). Limitations: record begins July 2002 | Jul 2002 – Dec 2024 | Monthly (from 4-day composites) | Global | 500 m |
| NLCD29,30 | (unchanged) | (unchanged) | 1985–2025 | Annual | CONUS | 30 m |

[G3, G4, G5, G6]. Also §2.1.4 "(Table 3)" while no Table 2 exists → renumber [G8].

28. Table 4 (replace the body):

| Resource | DOI | Number of gages | Contents | Formats |
|---|---|---|---|---|
| 1. Streamflow signatures, WY 1993–2025 | upon acceptance | 6,678 | Signature summary table (2,277 columns), annual values (27,766,578 rows), data dictionary, signature-category table, interactive explorer with its time-series archive | CSV, Parquet, HTML, ZIP, MD |
| 2. Streamflow signatures, WY 1980–2025 | upon acceptance | 6,250 | As Resource 1 (34,813,881 annual rows) | CSV, Parquet, HTML, ZIP, MD |
| 3. Harmonized hydroclimate data | upon acceptance | 8,014 | Daily streamflow with agency quality codes; gage metadata for all 16,994 candidates; HYDAT interference table (8,012 stations); ledger of the 193,542 restored USGS days; the two products' gage lists; the 31 + 29 would-qualify gages | Parquet, CSV, TXT, MD |
| | | 8,017 | Daily basin-averaged Daymet climate, 1980–2025; table of 101 flagged basins | |
| 4. Watershed geometry and basin attributes | upon acceptance | 8,014 | HydroATLAS-derived attributes (198 attributes, 211 columns) with dictionary | GeoPackage, GeoParquet, CSV, Parquet, MD |
| | | 7,964 / 8,017 | Watershed boundary polygons, simplified / full resolution, with QA table | |
| 5. Vegetation and land cover | upon acceptance | 7,964 | Monthly MODIS LAI, annual MODIS land cover, dictionary, interactive explorers | Parquet, CSV, HTML, MD |
| | | 6,119 | Annual NLCD land cover and impervious surface (CONUS), list of the 45 excluded Alaska gages, explorer | |

Delete "run provenance and validation reports" and "JSON" (not deposited) [J2, J3, J4].

## D. Manuscript — §3 text, cross-references, usage notes

29. §3 Resources 1 and 2: "(iii) a self-contained interactive HTML explorer that maps every signature-statistic combination and plots per-gage annual series with fitted checksums, configuration, and software versions) and validation reports" → "(iii) an interactive HTML map explorer of every signature-statistic combination with a sidecar archive of the per-gage annual series; and (iv) the data dictionary and a signature-category table" [J8]; "described in Section 2.1.2" → "Section 2.2.1" [J9].
30. §3 Resource 3: "The accompanying gage metadata table (16,994 candidate gages) records station coordinates, published drainage areas, processing status, and human-interference indicators drawn from GAGES-II attributes23,33 and Reference Hydrometric Basin Network (RHBN) and regulation status via HYDAT" → "The accompanying gage metadata table (16,994 candidate gages) records station coordinates, published drainage areas, retrieval extent and processing status. RHBN membership and regulation status for 8,012 Canadian stations form a separate HYDAT-derived table, and the GAGES-II human-interference indicators23,33 are carried as metadata columns of Resources 1 and 2" [J10]; add "Resource 3 also holds the ledger of the restored USGS days, the 101 basins with low-confidence outlines or few Daymet cells, the two products' gage lists and the gages that qualify only under the restored record" [J12].
31. §3 Resource 4: add "A second GeoPackage holds the same outlines at full resolution together with the 53 basins larger than 100,000 km² (8,017 polygons); the Daymet forcing of Resource 3 was aggregated over it" [J13]; "from BasinAtlas v1022,25" → "v1025".
32. §3 Resource 5: "(MCD15A3H v061; Myneni et al. 2015)" → "(MCD15A3H v061)27" [J14]; "granule-level provenance manifests for the MODIS products" → "the list of the 45 Alaska gages outside the NLCD footprint" (manifests are not deposited) [J17].
33. §3 Conventions: "(see unnormalized gages section 5.1.1)" → "(Sect. 6.1)"; delete the duplicate sentence "Files derived from the geometry layer also carry canon_id, a zero-stripped variant …"; "nlcd_c{code}_pct}" → "nlcd_c{code}_pct"; "All tables join on gage_id" → "All tables join on gage_id (named site_id in the Daymet table)" [J18].
34. §2 preamble cross-references: "(Sect. 2.1.1)" after the basin polygons → "(Sect. 2.1.2)"; "(Sect. 2.1.2)" after the two windows → "(Sect. 2.2.1)" [C6].
35. §6.2: "e.g. 0103500" → "e.g. 01030500"; "may be read as a character string" → "must be read as a character string"; "(Sect. 5.1.3)" → "(Sect. 6.3)" [N1, N2, N4].
36. §6.3: "should not be compared against the WY 1993–2025 and WY 1980–2025 products" → "between" [N4]; extend the list: "… all drought metrics, and the per-gage scalars derived from the window's record (elasticity_static, recession_alpha_point_cloud_linear_reservoir, the log_a_seasonality scalars, the drought thresholds)" [O1].
37. Data availability: "licenceCC-BY 4.0 license" → "license" [P1].
38. Acknowledgements: "Claude Code 0.145.0" (twice) → "Claude Code (Anthropic), releases 1.x–2.x, 2025–2026" [P2].
39. §4 (information only, no recommendation): "median watershed size of 698 km2, 25% … greater than 2,700 km2, … less than 180 km2" are the WY 1980–2025 product's quartiles (184 / 698 / 2,718 km²); the WY 1993–2025 set shown in Figure 2 gives 166 / 645 / 2,500 km² [K2].

## E. Citations

40. Table 1 "HydroBASINS22" → 34 (Lehner & Grill 2013); "Daymet21,26" → 26; Table 3 "Hydro ATLAS22" → 25 (Linke 2019); §3 Resource 4 "22,25" → 25 — reference 22 is the 2022 lake paper [Q1, Q3, G6].
41. Reference 27: the DOI is MOD15A2H (Terra, 8-day); the product used is MCD15A3H v061 (Terra + Aqua, 4-day): 10.5067/MODIS/MCD15A3H.061 [Q2].
42. Reference 23 (GAGES-II) is cited for the perennial-stream rule in §2.1.1 — confirm or cite the project's own rule [Q4].
43. References 29, 30, 36 carry no year; 36 (gdptools) no version [Q5].

## F. Guidelines doc ("Summary Documentation for Streamflow Signatures")

44. Part 1.1: "We filtered for gages with watersheds smaller than 85,000 km2 and at least 20 total years of valid daily flow observations, resulting in 8,014 gages" → "We kept gages with at least 20 calendar years each having at least 30 days of flow above 0.0001 mm/day, resulting in 8,014 gages" (no size filter) [R1]; "Albers et al., 2026" → "Albers, 2017".
45. Part 1.3: "at least 60% of the possible water years in the analysis window" → "at least 60 % of the water years from the window start to the gage's last water year within the window" [R2, H4].
46. §3.6 Snow: "for water years 1980–2023" → "for water years 1981–2025 (the Daymet input begins 1 January 1980, so WY 1980 has no climate-valid year)" [R5].
47. Part 4, generate_streamflow_dt(): "USGS values are kept only for qualifiers A, A e, P, P e (others set to NA)" → "A USGS value is kept whenever its qualification status is A or P and the value is finite; the remark code (e, R, <, >, [4], Ice, …) is carried in the flag column and never masks a value (censored values enter as the reported bound)" [R7 — DUE since the version-3 products].
48. Part 4, integrate_daymet_with_streamflow(): "data_out/daymet_1980_2023.parquet" → "the 1980–2025 Daymet table (daymet_1980_2025_29sep2026.parquet, 8,017 basins)" [R11].
49. Part 6 "Resource Collections": replace the working notes ("Insert the table …", "Here is the original …") with the Table 4 rows above, or delete [R9].
50. §3.1 "Negative_annual" → "negative_ann" (the shipped column base, lower case) — corrects item 1 of the 2026-10-12 relay, which wrote "Negative_ann" [R12].
51. Carried from the 2026-10-12 relay, all still open: naming-convention revert (item 1: Q_annual → Qann, Qsumr → Qsum, Qfall → Qfal, runoff_ratio_annual / _Xseason → annual/winter/spring/summer/fall_runoff_ratio, Flow_Reversals_xSeason → _winter/…); §3.7 names (precip_intermit_* → prcp_intermit_*, GDD_cum_* → gdd_cum_*, _xSeason → _win/_spr/_sum/_fal); dry day ≤ 0.1 mm and "strictly above" on the wet-day lines; prcp_5d_max window inside the water year and first-occurrence ties; prcp_annual_sd / _cv as within-year daily spread; Gini and seasonality-index formulas typed (pasted images today); frost "(NA when none)" and the Oct 1 – Jun 30 window; §3.7 header (functions, inputs, Requirements block); duplicate "1.1" numbering; "daicly"; parameterized-BFI interval "[0, 1]" → "(0, 1)"; Hecht et al. 2024 in the references; "recession_alpha_point_cloud" → "recession_alpha_point_cloud_linear_reservoir"; Part 1.2 "(see 3.12)" → "(see 3.1)"; prcp_50_dowy tie rule; the three convention lines (Gini range, all-zero year, Daymet leap-day fill).
52. Part 3 numbering is now 3.3 Flow Timing … 3.7 Climate Metrics … 3.9 Storage (reorder of 2026-10-30) — fine; nothing to change.
52a. Part 1.3: "For any seasonal or annual metrics the period must have at least 80% of the data" → "The four seasonal flow totals and the four seasonal runoff ratios require at least 80 % of the season's raw observations; water years are screened by the 30-day and 3-day rules of Part 1.2" [R2].
52b. §3.9: "The recession and parameterized-BFI families are exempt from the trend-completeness requirements and the 20-value statistics floor" → "The recession family is exempt …; the parameterized BFI signatures carry both gates" (the code applies them; BFI_Eckhardt_param has 6,507 means and 6,111 trend slopes in the WY 1993–2025 product) [R4].
52c. §3.9: "All recession metrics except n_recession_events require at least 25 recession events" → add "the whole-record constant recession_alpha_point_cloud_linear_reservoir needs more than 10 day-pairs and a value strictly between 0 and 1, and is computed independently of the 25-event gate" [R4].
52d. Part 2: "pettitt_pre_mean: Mean of the annual values before the changepoint year" → "through and including the changepoint year (the reported year is the last year of the pre-segment)"; "pettitt_post_mean: … after" → "strictly after it" [R3].

## G. Deposited data dictionary (for the HydroShare deposit, not the co-authors) — reviewer C, 0 HIGH
53. Seven MEDIUM omissions (`docs/plans/2026-10-30-manuscript-reviews/review_C_definitions.md`): fall-season origin on the four Q fall rows; "NA by design for un-normalized gages" on the 14 Q-to-PPT rows; the 25-event gate on the 12 recession rows; the parameterized-BFI NA rule; "NA in snow-free years" on 7 snow rows; avg_storage's interpolation step; drought "strictly below" / "0, not NA". Eight LOW wording items. Applying them means re-rendering and re-staging the dictionary of Resources 1 and 2 (user decision).
53a. Codex (MAJOR, beyond the `_mean` rows C reviewed): every `_pettitt_pct_change` definition ("percentage of the pre-changepoint mean") must say "of the ABSOLUTE pre-changepoint mean; NA when |pre-mean| ≤ 1e-10" — the sign matters for the 6,578 negative minimum-temperature pre-means in the WY 1993–2025 product (01011000: pre-mean −34.13 °C, change +4.03, reported +11.8 %); and the `_pettitt_pre_mean` / `_post_mean` rows should state that the pre-segment includes the changepoint year.
54. Readmes of Resources 1 and 2 (sibling session): "≥ 20 valid water years that are ≥ 60 % of the window's years" states the wrong denominator (325 gages of the 1980–2025 product have fewer than 60 % of 46 years) [H4]; "Precipitation-dependent signatures (runoff ratios, elasticity, Q–P seasonality, storage)" → "average storage" [M2].

# Review: Hydroclimate Information, Signatures and Summary Statistics - Vegetation and Land Cover

**Verdict:** FINDINGS (0 BLOCKER / 4 MAJOR / 4 MINOR / 0 NOTE)

Locations below refer to deposited files unless explicitly marked **repo**. README = `hisss_readme_eo.md`; dictionary = `hisss_eo_data_dictionary.csv`.

## Findings (ranked)

| # | Severity | File (location) | Claim or expectation | Evidence (what you computed/quoted) | Recommended fix |
|---|---|---|---|---|---|
| 1 | MAJOR | Dictionary, lines 100–136; `hisss_modis_lulc_explorer.html`, line 126, `BANDS.names` | The dictionary defines every column, including the meaning of classification codes. | **37 of 102 class columns lack meaningful class definitions.** Their labels merely repeat codes: `LCCS1_1`, `LCCS2_20`, `LCCS3_51`, etc. The explorer repeats these placeholders. All 37 columns contain nonzero observations: **2,887,885 nonzero cells** altogether. Example: gage `05DA010`, 2023, `lc_prop1_c1_pct = 59.777963004009074`; the deposit does not explain what that approximately 60% represents. | Supply the authoritative v061 class names and definitions for all 37 codes in the dictionary and explorer legends. |
| 2 | MAJOR | Dictionary, line 23, `n_composites`; `hisss_modis_lai_monthly.parquet` | Composite count is described as the rounded mean over “the basin’s tiles.” | **Repo** `EO_data_processing/eo_processing/lai_pipeline.py:173` instead assigns one scalar calculated over **all successfully processed tiles in that run** to every basin. The backfill does likewise at `lai_backfill_lpdaac.py:160`. Consistent with this, **242 of 270 months have one count across all 7,964 gages**; the other 28 have two counts. For May 2024, 7,872 gages receive 10 and 92 receive 8. This is not the documented basin-specific count; it also feeds `partial_month` and `good_coverage_frac`. | Correct the definition to disclose the run/cohort-wide proxy and its limitations, or recompute counts using each basin’s intersecting tiles. |
| 3 | MAJOR | Dictionary, line 26; README, lines 82–84; `hisss_modis_lai_explorer.html`, lines 91–97 | `good_coverage_frac` represents valid observations divided by all expected basin pixel-observations. | All **2,150,280 values** reproduce the implementation’s formula, rounded to four decimals, using a **month-of-year modal composite count**. That denominator ignores year-specific calendar differences. Example: `01011000`, March 2004: `n_composites=8`, `mean_obs_valid_frac=0.1926677390604162`, deposited coverage **0.1713**, because the March mode is 9. Enumerating the nominal four-day dates with the repository’s overlap rule gives **8** expected composites that month. Substituting calendar-specific counts changes **299,707 rows by more than 0.0001**, with maximum absolute difference **0.1111561**. These differences are not rounding errors. The run-wide numerator proxy in finding 2 further prevents interpreting this as an exact basin observation fraction. | Describe the deposited field explicitly as an empirical coverage proxy, including the modal denominator, rounding, and missing-value treatment; alternatively recompute the literal fraction from basin-specific observations and expected dates. |
| 4 | MAJOR | README, lines 140–149; `hisss_nlcd_explorer.html`, line 83, `POINTS.v.low_confidence`, `POINTS.lowconf`, and `VARS` | All explorers use the deposited boundary layer’s `low_confidence` flag. | The NLCD explorer instead embeds the **table’s union QA flag** and defines it as “geometry low-confidence OR <100 eff px OR partial coverage.” It marks **three gages**: `01154500`, `01170500`, `01184000`. Their deposited boundary flags are all **false**; their NLCD coverage fractions are respectively **0.9792286612, 0.9852934874, 0.9880295673**. Thus the explorer has three flagged points where the boundary-based flag would have zero among its 6,119 gages. LAI and LULC explorers correctly use the boundary flag. | Rebuild the NLCD explorer using the deposited boundary flag, consistent with the owner’s decision. Retain `partial_coverage` as its separate variable. |
| 5 | MINOR | Dictionary, header and lines 17–24; LAI Parquet schema and QA columns | Every column has a documented type and NA convention. | The dictionary has only `table,column,units,definition`, with types inconsistently placed in `units`. For example, LAI `low_confidence` is labelled `bool` but physically stored as **int64, values 0/1**. Six QA fractions—`valid_pixel_frac`, `mean_obs_valid_frac`, `frac_cloud`, `frac_snow_ice`, `frac_aerosol`, `frac_backup_algorithm`—each contain **5,361 NAs across 65 gages**, all at zero `n_basin_pixels`. Their definitions give ranges but omit this zero-denominator NA rule. | Add explicit data types and document zero-denominator NAs, distinguishing them from `good_coverage_frac`, which is zero in those rows. |
| 6 | MINOR | `hisss_modis_lai_explorer.html`, lines 96–97 | Public-facing wording should avoid suggesting an earlier version. | The explanation says `partial_month` “**is the older version of this idea**.” This is conceptual version-history wording rather than an explicit release claim, but it is avoidable under the stated publication rule. | Say: “`partial_month` is a complementary binary flag for fewer than six composites.” |
| 7 | MINOR | README, line 3 — **only occurrence found across this resource’s text files and explorers** | HISSS should use the renamed expansion. | README still says “**Hydrologic Information Signatures and Summary Statistics**.” The supplied live title and abstract use “Hydroclimate Information, Signatures and Summary Statistics.” | Update this one occurrence to the current expansion. |
| 8 | MINOR | **Repo** `EO_data_processing/README.md`, lines 64 and 330; deposited README, line 21 | Repository reference material and deposit should agree about LULC width. | Repository documentation says approximately **200 class columns**. The deposited table has exactly **102 class columns**, within **109 total columns**; the deposited README is correct. | Update the repository’s current-product description to 102 class columns, or clearly label the approximately 200-column text as an internal planning estimate. |

## Checks performed

All checks were offline and read-only, using the specified Python interpreter. Parquet tables received one complete batched scan, followed by targeted column projections. Observed peak Python memory was approximately **1.48 GB**.

### File, schema, and panel completeness

All **nine listed resource files exist**. Every substantive file is described in the README; referenced signature, metadata, and geometry files exist in the named companion resources. No ZIP is present or promised.

| Table | Rows | Gages | Physical columns | Period | Missing or duplicate gage–period rows |
|---|---:|---:|---:|---|---:|
| MODIS LAI | 2,150,280 | 7,964 | 26 | July 2002–December 2024; 270 months | 0 |
| MODIS LULC | 191,136 | 7,964 | 109 | 2001–2024; 24 years | 0 |
| NLCD | 250,879 | 6,119 | 29 | 1985–2025; 41 years | 0 |
| NLCD exclusions | 45 | 45 | 5 | Not applicable | 0 duplicate IDs |

These counts agree with the README and supplied live abstract.

The dictionary contains **173 rows**: **169 actual-column definitions plus four explicitly identified caveats**. There are **zero missing column names, nonexistent non-caveat columns, duplicate definitions, or empty dictionary cells**. Semantic completeness issues are findings 1 and 5.

### Values, sums, flags, and caveats

- **MODIS LULC:** eight schemes contain respectively **17, 16, 11, 9, 12, 16, 11, 10** class columns. Every scheme sums to 100 in every row; maximum absolute error **1.14 × 10⁻¹³ percentage points**. No missing values or percentages outside 0–100.
- **NLCD:** all 16-class compositions sum to 100; maximum error **1.44 × 10⁻¹² percentage points**. No missing values. `impervious_mean_pct` ranges **0–68.8237638841**. Every `nlcd_collection` value is **C1V2**.
- **NLCD formulas:** zero failures for the low-pixel, partial-coverage, and union flags. `valid_area_km2` agrees with `n_nlcd_pixels × 0.0009` within **0.00005 km²**, consistent with four-decimal rounding.
- **NLCD coverage:** range **0.9792286612–1**; no within-gage changes across years. Exactly **three** partially covered gages, yielding **123 flagged annual rows**. Zero geometry-low-confidence or low-pixel-support rows.
- **LAI:** all declared fraction ranges hold; zero quantile-order violations, means outside min/max, or inconsistent missingness among the nine LAI statistics. `valid_pixel_frac = n_eff_pixels / n_basin_pixels` exactly wherever defined.
- Each LAI statistic has **62,799 missing values**. Exactly **17 gages** are all-NA across 270 months, as documented. Their MODIS IGBP urban percentage is **100 in every annual row**, supporting the urban explanation.
- `partial_month` exactly equals `n_composites < 6`: **31,580 rows**, comprising April 2022 (**7,964**) and October 2022, October 2024, December 2024 (**7,872 each**).
- LAI flags **57 gages**; LULC flags **62**, including **six** low-pixel-support gages. LULC’s flag equals LAI’s geometry flag OR low-pixel support in every row.
- Boundary flags total **56**. LAI differs only at documented gage `04HA001`; its stored area discrepancy is **24.3081761%**.
- The dictionary’s NLCD artifact statistics are supported: shrub/grass annual-change correlation **−0.7143602**; maximum absolute annual changes **77.4967** and **89.8240 percentage points**. Developed percentage versus impervious percentage correlation is **0.9631177**.

### IDs and companion-resource joins

Every EO `canon_id` equals the documented zero-stripped key. There are **zero canonical-key collisions**.

| EO table → companion | Matched EO gages | Unmatched EO gages | Companion gages without this EO table |
|---|---:|---:|---:|
| LAI or LULC → boundary QA | 7,964 | 0 | 0 |
| NLCD → boundary QA | 6,119 | 0 | 1,845 |
| LAI or LULC → gage metadata | 7,964 | 0 | 9,030 |
| NLCD → gage metadata | 6,119 | 0 | 10,875 |
| LAI or LULC → WY1993–2025 signatures | 6,634 | 1,330 | 44 |
| NLCD → WY1993–2025 signatures | 5,454 | 665 | 1,224 |
| LAI or LULC → WY1980–2025 signatures | 6,205 | 1,759 | 45 |
| NLCD → WY1980–2025 signatures | 5,007 | 1,112 | 1,243 |

- The metadata table contains **16,994 unique IDs**, with zero padded USGS IDs.
- Confirmed **44 unpadded long EO IDs** against the companion attribute table: **14 originally nine-digit and 30 originally ten-digit IDs**. Exact-string signature joins lose **35 / 9 matches**, exactly as documented. Exact-string joins to boundary QA lose **35**; canonical joins lose none.
- All **44 / 45** signature gages lacking MODIS records belong to the **53 large basins** identified in the unsimplified layer.
- The EO universe contains **6,164 US gages**. Its **45 US gages absent from NLCD exactly equal the exclusion list**. Exclusion coordinates match both boundary QA and metadata exactly; their ranges are **55.176906–70.281385° latitude**, **−163.677373–−131.195583° longitude**.
- Simplified and unsimplified GeoPackages each have layer `watersheds`, CRS **EPSG:4326**, with **7,964 / 8,017 unique gages**. The latter adds exactly **53 polygons above 100,000 km²**. `05KH009` is absent from both and present in the companion attribute table.
- A **20-feature** simplified-geometry sample had zero invalid geometries. Recomputed areas in **EPSG:6933** matched stored areas with maximum relative difference **3.94 × 10⁻¹⁵**.

### Explorers

| Explorer | Embedded unique gages | Map variables | Time axis |
|---|---:|---:|---|
| LAI | 7,964 | 17 | 270 months |
| MODIS LULC | 7,964 | 27 | 24 years; eight schemes |
| NLCD | 6,119 | 32 | 41 years |

- All embedded gage sets match their tables exactly; all checked array dimensions and variable definitions match their payloads.
- Coordinates agree with boundary QA within the displayed **five-decimal rounding**.
- Sampled **12,960 LAI**, **14,688 MODIS class**, and **3,280 NLCD class** time-series values match their tables within the explorers’ **0.1-unit encoding**. Sampled LAI coverage strips also match within integer-percentage rounding.
- All **7,947 nonempty LAI record means** agree within **0.0005**. All **250,879 NLCD impervious trace values** agree within **0.05 percentage points**; all 6,119 displayed impervious means agree within **0.005**.
- Data are embedded. External dependencies are Leaflet **1.9.4 JS/CSS from unpkg**, **Esri World Light Gray Canvas**, and **OpenStreetMap tiles**, consistent with the README’s connectivity caveat. No CARTO references were found.
- LAI/LULC explorer boundary flags match exactly. The NLCD exception is finding 4.

## Not checked (and why)

- Original raster-to-basin aggregation, source-granule identities, and historical processing/retrieval dates: source rasters and complete original manifests are not deposited; the README acknowledges the missing MODIS manifest.
- Authoritative external class-code meanings, product release histories, and provider methodological claims: network use was prohibited. Finding 1 concerns demonstrably absent definitions, without inventing replacement meanings.
- Actual browser rendering, external tile availability, and automatic basemap fallback: inspected embedded data and code only.
- The README’s **17,166-station collection-wide union**: validated this resource’s joins and the 16,994-row metadata table, but did not scan every companion daily-data identifier.
- Geometry validity beyond the 20-feature sample, and exhaustive recomputation of every derived explorer summary variable.

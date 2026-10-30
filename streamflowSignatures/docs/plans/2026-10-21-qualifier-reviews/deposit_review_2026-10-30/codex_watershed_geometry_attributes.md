# Review: Hydroclimate Information, Signatures and Summary Statistics - Watershed Geometry and Basin Attributes

**Verdict:** FINDINGS (2 BLOCKER / 7 MAJOR / 4 MINOR / 0 NOTE)

## Findings (ranked)

Locations below use these abbreviations:

- **README:** `hisss_readme_geometry_attributes.md`
- **Dictionary:** `hisss_hydroatlas_attributes_dictionary.csv`
- **Attributes:** `hisss_watershed_hydroatlas_attributes.csv` and `.parquet`; CSV line numbers include the header.
- **Abstract:** supplied `live_metadata/watershed_geometry_attributes.txt`, line 4.
- **Catalog:** locally available [BasinATLAS_Catalog_v10.pdf](/Users/arik/Downloads/HydroATLAS/BasinATLAS_Catalog_v10.pdf), using PDF page numbers. No network was used.

| # | Severity | File (location) | Claim or expectation | Evidence (what you computed/quoted) | Recommended fix |
|---|---|---|---|---|---|
| 1 | **BLOCKER** | Dictionary lines 21, 24, 98, 101, 103–104; corresponding Attributes columns | Six columns omit essential scaling factors. | Catalog pp. 5, 8, 49, 52, 54 specifies: `lka_pc_use` and `dor_pc_pva` are **percent ×10**; `pop_ct_usu` is **thousands of people**; `nli_ix_uav` is **index ×100**; both `hft_ix_u93/u09` are **index ×10**. All six deposited columns exactly match the local source values for all **8,014 rows**, confirming they remain scaled. Examples: `04087204` lake value **757 = 75.7%**; `05NB036` regulation **10000 = 1000%**; `02OA016` population **40620.859375 thousand**; `01100600` lights **6300 = 63**; `01208950` footprint **451 = 45.1** in both years. | Correct the dictionary’s stored units and conversion factors. Preserve the numeric values unless deliberately converting the entire columns. |
| 2 | **BLOCKER** | Attributes, especially CSV lines 4880, 4896, 4831; README lines 8–9; Abstract | Attributes describe each gage’s full upstream contributing area. Some instead describe vastly larger networks than the deposited watershed. | For `07001910`, reported area is **3 km²**, deposited boundary area **3.0136158 km²**, but attribute aggregation covers **1,828,065.8 km² / 14,079 subbasins**. For `07010208`, boundary **6.7360 km²**, aggregation **1,829,180.3 km²**. **36 gages** with reported area ≤100 km² receive aggregations exceeding 100,000 km². Overall, **2,020** attribute rows have `watershed_area_rel_diff > 0.5`. The README acknowledges possible mis-snapped outlets, but the actual extreme assignments remain populated as watershed attributes. | Correct the affected outlet assignments or explicitly mark their attributes unusable pending correction. Quantify this limitation prominently and qualify the full-watershed claim. This does not require changing the frozen gage set. |
| 3 | **MAJOR** | README lines 103–105; Attributes `watershed_area_sum_SUB_AREA`, `outlet_UP_AREA` | Outlet attributes are accumulated over the “same upstream extent” as the re-aggregated attributes. | Relative difference between these two area fields exceeds **1% for 118 gages**, **5% for 32**, and **25% for 5**. `05HE001`: **996.5 versus 494.2 km²**, a **101.639%** difference. `11530500`: **40,894.9 versus 31,130.4 km²**, a **31.366%** difference. Independently traversing local `NEXT_DOWN` topology reproduced the deposited member counts and sums for both examples. | Remove the unconditional same-extent assertion. Document and flag these discrepancies; reconcile the spatial supports before presenting the two attribute families as directly comparable. |
| 4 | **MAJOR** | Dictionary line 12, `outlet_ORDER` | “Strahler stream order at the outlet basin.” | Every deposited value matches source HydroBASINS `ORDER_`. Source topology contains **55,428 downstream order decreases and zero increases**, incompatible with Strahler ordering. Concrete deposited pair: `09510200`, outlet `7120689160`, order **4**, drains directly to `09511300`, outlet `7120690850`, order **3**. There are also **116 order-zero rows**, all coastal. | Describe the actual HydroBASINS ordering convention and its zero code. Do not label this field Strahler order. |
| 5 | **MAJOR** | Dictionary lines 105 and 196, `gdp_ud_usu` / `gdp_ud_sav` | Both are described simply as “Gross domestic product (USD).” | Catalog p. 56 distinguishes **GDP totals** (`su`) from **GDP per capita** (`av`), using **2011 international dollars at PPP**, for **2015**. Deposited `gdp_ud_sav` ranges **26,724.684–71,226**; the independently reproduced value for `01648000` is **71,226**, an area-weighted per-capita measure, not watershed GDP total. | Distinguish total GDP from per-capita GDP; specify PPP, currency basis, reference year and spatial weighting. |
| 6 | **MAJOR** | Dictionary header and coverage; boundary GeoPackages/GeoParquet; `hisss_watershed_boundaries_qa.csv` | Every table column needs a usable definition, type and unit. | All **211 attribute columns** have dictionary rows, but **seven boundary/QA fields** do not: `canon_id`, `watershed_geom_source`, `geom_area_km2`, `area_rel_diff`, `area_flag`, `geom_simplified`, `low_confidence`; `geometry` is also absent. The dictionary has **no type field**, and **10 numeric diagnostic columns** have blank unit cells. The README never states the operational `area_flag` threshold: computed data show exactly **28 rows with relative difference >0.5** in each layer. The two `_cl_smj_frac` rows incorrectly retain unit code `cl` although their values are percentages. | Add a boundary/QA dictionary and complete types, units, formulas and flag thresholds. Distinguish boundary confidence from EO pixel-support flags. Change the fraction columns’ units to percent. |
| 7 | **MAJOR** | README lines 24, 117–119; Dictionary categorical rows | The BasinATLAS catalog supplies the class-code legends. | The catalog itself directs users to **`HydroATLAS_v10_Legends.xlsx`**—for example pp. 15, 16, 24, 26 and 28. Neither that workbook nor equivalent legends are deposited. **11 categorical columns** and **48 class-percentage/group columns** lack deposited code mappings. Actual values include `08ME023` with `tbi_cl_smj=99`, `tec_cl_smj=0`; neither meaning is supplied. `outlet_ENDO` lists 0/1/2 without defining their meanings. | Deposit the relevant versioned legends, or provide precise links and complete mappings. Explain observed exceptional codes, wetland groupings, and the month/suffix convention for the **84 monthly columns**. |
| 8 | **MAJOR** | Dictionary line 202, `wet_cl_smj`; boundary/attribute missing-value conventions | Missing values must distinguish absence of a feature from unavailable information. | `wet_cl_smj` is null for **2,038/8,014 gages (25.43%)**. Every one has `wet_pc_ug1=wet_pc_ug2=0`. Catalog p. 28 explains that source −9999 denotes **no wetland areas**; repository aggregation masks this to NA. This meaning is absent from the deposited definition. Separately, boundary area diagnostics are null for **68 simplified / 73 unsimplified rows**, yet their `area_flag` is false because no comparison can be made. | Document the wetland-absence NA convention and missing drainage-area behavior explicitly. Explain that a false area flag does not establish area agreement when its diagnostic is NA. |
| 9 | **MAJOR** | Abstract | The live description should identify the deposited boundary products and their coverage. | The abstract describes **7,964 boundaries and 29 fallbacks**, but omits `hisss_watershed_boundaries_unsimplified.gpkg`: **8,017 features, 30 fallbacks**, including **53 additional large basins**. This omitted layer supplies polygons for the **44 / 45 signature gages** absent from the simplified layer. Its categorical aggregation description also omits the documented `glc`/`pnv` upstream-fraction method. | Update the abstract to distinguish both boundary layers, their coverage and Daymet role; describe the categorical aggregation exceptions and normalized join key. |
| 10 | **MINOR** | Dictionary line 2 | USGS identifiers are “8-digit zero-padded.” | Of **6,160 USGS attribute IDs**, **73** are longer: **21 nine-digit, 48 ten-digit, 4 fifteen-digit**. README lines 152–154 already acknowledge this dictionary error. | Correct the dictionary itself to specify string IDs with variable length and the zero-stripping join rule. |
| 11 | **MINOR** | README lines 23, 59–60 | CSV decimals agree within `5e-10`; median area difference for the remaining simplified basins is 0.057%. | Exact decimal-versus-binary comparison found **one** violation: Attributes CSV line **4881**, `07001985`, `watershed_area_rel_diff`, error **5.2427579463e-10**. Separately, the **7,823 simplified basins** have median area difference **0.0589876%**, rounding to **0.059%**. The quoted **0.057%** is the median across **all 7,964 shared basins**, including the 141 retained at full resolution. | Relax the export tolerance slightly and identify the median’s denominator, or report 0.059% for the simplified subset. |
| 12 | **MINOR** | README line 3 | Dataset expansion should match the renamed resource. | The sole occurrence in this resource’s text files is **“Hydrologic Information Signatures and Summary Statistics”** at README line 3. The supplied live title uses **“Hydroclimate Information, Signatures and Summary Statistics.”** | Update that expansion. |
| 13 | **MINOR** | Repository context: `EO_data_processing/geometry/README.md` lines 4–5, 37, 48–49; `EO_data_processing/README_DAYMET.md` lines 33–35; `R/aggregate_hydroatlas_metadata.R` lines 346, 350 | Repository reference material should agree with the deposit. | Geometry documentation calls the four extra USGS gages “signature-bearing”; all four are absent from both current products. It promises zero-padded IDs despite the **nine documented exceptions**, and identifies `04HA001` as low-confidence although deposited `area_rel_diff=0.24308176`, `low_confidence=False`. Daymet documentation says Resource 4 omits the 53 large basins now deposited in the unsimplified layer. The dictionary generator still contains erosion/HDI unit placeholders, unlike the deposited definitions. | Synchronize these reference documents and the dictionary generator with the deposited products. |

## Checks performed

**Inventory and schemas**

- All **eight deposited files** exist and are described in the README.
- Attribute CSV and Parquet: **8,014 rows × 211 columns**, identical column order, identifiers, strings and missing-value locations. Exact numeric export exception is finding 11.
- Dictionary: **211 unique rows**, no missing or surplus attribute-column names.
- Aggregation definitions: **91 passthrough**, **92 weighted means**, **9 categorical modes**, **2 categorical argmax**, **2 dominant fractions**, **2 elevation extrema**, **13 diagnostics**.
- Both GeoPackages contain layer **`watersheds`**, CRS **EPSG:4326**, with 12 nongeometry data fields. GeoParquet also declares EPSG:4326.

| Boundary check | Simplified | Unsimplified |
|---|---:|---:|
| Features / unique gage IDs / unique normalized IDs | 7,964 | 8,017 |
| GAGES-II | 6,164 | 6,164 |
| WSC/ECCC | 1,771 | 1,823 |
| HydroBASINS fallback | 29 | 30 |
| `geom_simplified=True` | 7,823 | 0 |
| `area_flag=True` | 28 | 28 |
| `low_confidence=True` | 56 | 57 |
| Missing published area / area-relative difference | 68 | 73 |

- Canadian official-polygon shares are **98.3889% / 98.3810%**, consistent with rounded 98%.
- The 53 added polygons are all Canadian, with geometry areas **100,761.3187–1,696,281.5035 km²**. Their sole low-confidence feature is `02OA024`.
- `05KH009` is absent from both boundary layers.
- All `canon_id` values equal the documented normalized key. Flag formulas reproduce exactly; no unexpected boolean or source codes.
- Simplified area agreement: **7,726/7,896 = 97.8470%** within 10%; median relative difference **0.00100766**.
- QA CSV matches the simplified layer apart from floating-point text rounding: maximum area difference **7.28e-12 km²**.

**Geometry**

- Checked **all 7,964 simplified geometries**: zero invalid, empty or null.
- Simplified GeoPackage and GeoParquet geometries match exactly as WKB after promoting single polygons to multipolygons; all nongeometry attributes match.
- Recomputed every simplified area in **EPSG:6933**: maximum relative discrepancy **2.03e-12**.
- Checked **230 unsimplified geometries**, including all 141 retained-resolution basins, fallback basins and selected large basins: zero invalid, empty or null; recomputed EPSG:6933 areas match stored values.
- All **141 retained-resolution basins** have matching coordinate-array shapes; maximum coordinate difference **1.27898e-13 degrees**, within the README bound.
- Simplified duplicate outlines: **8 groups / 10 additional gages**, as documented. Their sampled unsimplified counterparts form **7 duplicate groups**.
- Maximum shared-layer area difference: **1.990269%**, below 2%.

**Cross-resource joins**

Applied zero stripping to numeric identifiers on **both sides**, retaining Canadian identifiers. Each cell below is **companion-only / resource-table-only**.

| Companion | Simplified boundaries | Unsimplified boundaries | Attributes |
|---|---:|---:|---:|
| Gage metadata: 16,994 IDs | 9,030 / 0 | 8,977 / 0 | 8,980 / 0 |
| Daily streamflow: 8,014 IDs | 54 / 4 | 1 / 4 | 0 / 0 |
| Daymet: 8,017 IDs | 53 / 0 | 0 / 0 | 4 / 1 |
| WY1993–2025 signatures: 6,678 IDs | 44 / 1,330 | 0 / 1,339 | 0 / 1,336 |
| WY1980–2025 signatures: 6,250 IDs | 45 / 1,759 | 0 / 1,767 | 0 / 1,764 |
| MODIS LAI and LULC: 7,964 IDs each | 0 / 0 | 0 / 53 | 4 / 54 |
| NLCD: 6,119 IDs | 0 / 1,845 | 0 / 1,898 | 4 / 1,899 |

- **No unexplained join failures.** The four boundary-only IDs are `01591000`, `01591400`, `01591610`, `01591700`; the single attribute-only ID relative to unsimplified boundaries is `05KH009`.
- Attributes exactly match the metadata’s **8,014 successful gages**.
- Verified all **nine** documented boundary/attribute spelling differences join after normalization.
- The overall ID union reaches the documented **17,166**: metadata’s 16,994 plus **172 additional HYDAT interference-table stations**.
- NLCD’s missing US boundary IDs exactly equal the **45-gage Alaska exclusion file**.
- All **101 Daymet flags** join to unsimplified boundaries: **57 geometry-confidence flags**, **45 low-pixel-support flags**, one overlapping. Areas, sources and geometry-confidence values agree.
- EO confidence differences are expected: LULC differs from boundary confidence for **6 gages**, NLCD for **3**; NLCD’s separate `geom_low_confidence` agrees throughout.

**Attributes and numerical claims**

- Independently reproduced member counts, areas, three weighted attributes, elevation extrema and wetland modes for **10 gages** using local HydroBASINS topology. Differences were floating-point rounding only.
- Both dominant-class columns reproduce the argmax of their upstream fractions for **all 8,014 rows**; both dominant fractions equal the corresponding maximum.
- No numeric infinity or residual **−9999** values.
- Discharge, elevation and temperature minimum/mean/maximum ordering checks passed.
- HDI range **793–955**, consistent with its documented ×1000 scale.
- Integer class-percentage sums: land cover **97–103** (**3,120 rows not exactly 100**); potential vegetation **98–102** (**826 not exactly 100**). Soil fractions sum **93–101**. No exact-100 claim is made in the deposited text.
- Recomputed the 2023 simplified-versus-full-resolution Daymet comparison over **7,964 basins × 365 days**: precipitation p01/p99 **−0.09888% / +0.10032%**, maximum **+1.21358%**; SWE p01/p99 **−2.52711% / +2.07786%** among 6,311 basins with positive reference SWE. These support the README’s approximate claims.
- GeoPackage timestamps are consistent with the stated August simplified build and October 30 unsimplified export.
- No prohibited earlier-public-release wording found. Internal build provenance was not treated as a release-history error.

## Not checked (and why)

- Full validity, duplicate search and source-WKB equality for all **8,017 unsimplified geometries**: sampled 230, following the requested geometry-reading constraint.
- Historical byte identity with the lost June layer, or actual agency download dates: cannot establish independently from the deposited files.
- Every attribute aggregation for every gage: checked selected aggregations independently and all values of the six scaling-error columns against the local source.
- Full Daymet daily values and original raster aggregation: projected identifiers from the large deposited table; used the local 2023 comparison products for the quoted sensitivity claims.
- Signature statistics, explorers and ZIP contents: none are deposited in this resource; companion signature tables were used for ID coverage only.
- Remote HydroShare state and external legend availability: network prohibited; reviewed the supplied live-metadata snapshot.

No files were created or modified. All Python used the specified interpreter; DuckDB used a 1 GB memory limit and four threads.

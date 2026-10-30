# Review: Hydroclimate Information, Signatures and Summary Statistics - Harmonized Input Data

**Verdict:** FINDINGS (1 BLOCKER / 3 MAJOR / 2 MINOR / 0 NOTE)

Locations below refer to this resource unless a companion resource is named. `README` means `hisss_readme_inputs.md`; `dictionary` means `hisss_input_data_dictionary.csv`.

## Findings (ranked)

| # | Severity | File (location) | Claim or expectation | Evidence (what you computed/quoted) | Recommended fix |
|---|---|---|---|---|---|
| 1 | BLOCKER | Live metadata, `harmonized_input_data.txt`, line 4 | Daymet covers “6,087 gaged basins” over “1980–2023”; dictionary covers “four tables.” | Full scan found **8,017 sites × 46 years × 365 dates = 134,605,430 rows**, covering **1980–2025**. The abstract omits **1,930 basins and two years** from its coverage description. There are **seven data tables**, plus two selection lists. README lines 21–27 correctly describe these deposits. | Update the live abstract’s coverage and inventory to match the deposited files. |
| 2 | MAJOR | README lines 101–106, 125–128; dictionary line 55; companion signature READMEs lines 117–120 | Qualification requires ≥20 valid water years covering ≥60% of “the window’s years.” The dictionary instead says possible years are “bounded by [the gage’s] record,” without an exact formula. | Independent replay reproduced both products’ year counts, but **325 of 6,250** selected gages fail 60% of the full 46-year window. Example: `15051010` has **26 valid years**, **26/46 = 56.52%**, but qualifies using **26/38 = 68.42%**. Its last in-window year with rows is 2017. The would-qualify table also contradicts the literal rule: line 55 gives `10308800`, **23/38 = 60.53%**, whereas **23/46 = 50%**. The actual denominator is **last in-window water year with rows − window start + 1**, not the full window or the span between first and last observations. The deposit also omits the numeric missing-day rejection threshold. | State the exact denominator, distinguishing raw-row years from valid years. Define valid-year screening locally: ≤30 missing days on the complete water-year calendar, no missing run >3 days, internal interpolation only, and no residual endpoint missing values. Preserve the frozen cohorts. |
| 3 | MAJOR | Dictionary lines 65–69; `hisss_daymet_basin_flags.csv`, line 90; README lines 113–115 | `no_area_reference` means “No published drainage area to compare the polygon against.” Other flags refer to an unspecified area tolerance and “few” cells. | `08HD026` has missing `metadata_area_km2` but **`no_area_reference=False`**. Across the unsimplified layer, **73** basins lack published area: **27 HydroBASINS** and **46 official ECCC** outlines. Only the 27 fallback basins receive this flag. Of the 46 official outlines, 45 are absent from the flags table; `08HD026` appears for low pixel support. Repository code explicitly restricts this flag to fallback polygons. Computed flags exactly follow **relative area difference >0.5** and **`coverage_sum <4`**; these thresholds are absent from the deposited definitions. | Define `no_area_reference` as the conjunction of HydroBASINS fallback and missing published area. Explain the official-outline exception and publish both numeric thresholds. |
| 4 | MAJOR | Dictionary line 10; `hisss_streamflow_daily.parquet`, `flag` | Agency quality codes need an interpretable legend. | Seven observed tokens are listed without their meanings: **`Dis`, `Eqp`, `Fld`, `Ssn`, `Mnt`, `***`, `[4]`**. They occur on **6,868 rows**, including **419 finite observations carrying `A [4]`** and 6,449 missing-value rows. Merely retaining or enumerating these codes does not define them. | Add authoritative, agency-specific definitions for every observed token, retaining the documented distinction between USGS and Canadian `A`. |
| 5 | MINOR | Dictionary header and lines 51–56 | Every table column should have its own name, type, unit, and definition. | The seven data tables contain **87 columns**; only **69** have exact `(table, column)` dictionary matches. The other **18**, all in `hisss_would_qualify_not_added.csv`, are mentioned inside four combined column-name entries. The two selection filenames likewise share one table-name entry. The dictionary has no explicit datatype field; grouped units mix counts, year lists, fractions, and booleans. These fields are mentioned, but not individually specified. | Expand to one row per actual table/column and add explicit types, units, and year-list/empty-value conventions. |
| 6 | MINOR | README line 3 | HISSS should use the renamed expansion. | The resource’s sole occurrence of the old expansion is **“Hydrologic Information Signatures and Summary Statistics (HISSS)”** at this line. The supplied live title and abstract already use the new name. | Replace with **“Hydroclimate Information, Signatures and Summary Statistics.”** |

## Checks performed

- **Inventory:** All **11 deposited files** exist. The README describes every accompanying file; no named input file is missing. Both signature READMEs’ pointers to their selection list and `hisss_would_qualify_not_added.csv` resolve correctly.

- **Streamflow — complete scan of 111,624,189 rows:**
  - **8,014 unique gages**; zero duplicate gage/date pairs and zero out-of-order dates within gages.
  - USGS dates: **1979-10-01–2025-10-01**; Canadian dates: **1979-10-01–2024-12-31**.
  - **394,210 missing Q values**, zero infinities, **21,518 negative values on 82 gages**, and **3,496,220 zeros**.
  - **7,028,100 absent dates within individual record spans**, affecting **3,553 gages**; absent rows are explicitly allowed by the README.
  - **5,381 rows** occur on 2025-10-01, the single included date of WY 2026.
  - All five derived calendar columns matched independently calculated dates: **zero discrepancies**.
  - **15,257 provisional-remark rows** have missing discharge. The **12,103 USGS rows with null flags** also all have missing discharge.

- **Gage metadata:**
  - **16,994 rows**, with no duplicate original or zero-stripped keys.
  - Status counts exactly match the dictionary: **8,014 success; 8,787 no_data; 80 insufficient_years; 113 processing**.
  - Recomputed `num_years`, `start_year`, `end_year`, and `num_days` for all compiled gages: **zero discrepancies**.
  - All USGS metadata IDs are zero-stripped.
  - **1,601** candidates lack area and have `area_normalized=False`; **73** are compiled. The products contain **32 / 28**, with **37 distinct gages** across their union.
  - Product coordinates, drainage areas, agency labels, and normalization flags match metadata: **zero discrepancies**.

- **Qualifier correction ledger:**
  - **193,542 unique gage/date rows on 987 gages**.
  - Every ledger date exists in streamflow; **zero stored-Q or retained-old-flag mismatches**.
  - Every `old_Q` is missing; every status is `restored`; every served unit is `ft3/s`.
  - Numeric and string served values agree; `Q_new = value_today_cfs × conversion` agrees exactly.
  - Drainage areas match metadata exactly; conversion-factor discrepancies are at most **2.22 × 10⁻¹⁶**.
  - Retrieval timestamps span **2026-10-26 18:29:06–18:49:21 UTC**.
  - Current status counts are **193,446 approved / 96 provisional**, matching the README.

- **Qualification and frozen selections:**
  - Independent full-calendar screening reproduced all selected gages’ valid-year counts and first/last valid years: **zero discrepancies**.
  - Qualifying totals: **6,709 / 6,279**; subtracting frozen products leaves exactly **31 / 29** listed exclusions.
  - Reconstructed baseline years by masking ledger dates. All **60 rows × 19 columns** of the would-qualify table matched.
  - Selection lists contain **6,678 / 6,250 unique IDs** and exactly reproduce product row order.
  - No selected Canadian gage has a qualifying year after **WY 2024**.

- **Daymet — complete scan of 134,605,430 rows:**
  - **8,017 sites**, each with exactly **16,790 dates** over 1980–2025.
  - **Zero missing expected dates, unexpected dates, duplicate site/date pairs, null values, or infinities.**
  - Leap years retain February 29 and omit December 31.
  - Zero instances of `tmin > tmax`; precipitation, SWE, vapor pressure, and radiation are nonnegative.
  - Observed ranges: precipitation **0–436.642 mm/day**; minimum temperature **−58.415–30.443 °C**; maximum temperature **−53.939–47.472 °C**; SWE **0–3,684.502 mm**; vapor pressure **1.890–4,091.874 Pa**; radiation **0–901.724 W/m²**.

- **Daymet flags:**
  - **101 unique basins**: **57 low-confidence**, **30 fallback**, **28 area mismatch**, **27 no-area-reference**, **45 low pixel support**; **one** overlaps the low-confidence and low-support groups.
  - Product membership: **46 / 33**, **53 distinct** across products; all membership and normalization fields match the products.
  - Flagged geometry sources, confidence flags, and areas agree with unsimplified geometry; maximum area discrepancy is **4.55 × 10⁻¹³ km²**.
  - Low-support basins intersect **2–12 cells**; only **five** intersect fewer than four. The operative threshold is fractional coverage, as noted in finding 3.

- **HYDAT interference:** **8,012 unique stations**; RHBN counts **316 true / 7,696 false**; regulation counts **3,116 true / 4,525 false / 371 missing**. All **7,840 Canadian metadata candidates** join; the interference table contains **172 additional stations**.

**Cross-resource joins**, after stripping leading zeros:

| Comparison | Unmatched on first side | Unmatched on second side |
|---|---:|---:|
| Compiled streamflow ↔ successful metadata | 0 | 0 |
| Compiled streamflow ↔ Daymet | 1: `05KH009` | 4: `01591000`, `01591400`, `01591610`, `01591700` |
| 1993–2025 product ↔ compiled streamflow | 0 | 1,336 |
| 1980–2025 product ↔ compiled streamflow | 0 | 1,764 |
| 1993–2025 product ↔ Daymet | 0 | 1,339 |
| 1980–2025 product ↔ Daymet | 0 | 1,767 |
| Daymet ↔ unsimplified boundaries | 0 | 0 |
| Compiled streamflow ↔ simplified boundaries | 54 | 4 |
| Compiled streamflow ↔ MODIS LAI, separately LULC | 54 | 4 |
| Compiled streamflow ↔ NLCD | 1,899 | 4 |

- The unsimplified GeoPackage contains **8,017 unique normalized IDs**, layer **`watersheds`**, **EPSG:4326**, including **53 polygons above 100,000 km²**. Nine literal ID spellings differ from Daymet; normalization resolves all nine.
- Naive streamflow-to-metadata ID matching gives exactly **3,118 matches**; padding metadata USGS IDs to eight digits leaves exactly **66 unmatched**, confirming the documented warning.
- No wording implying an earlier public release was found in this resource’s text files.

## Not checked (and why)

- Original agency observations, HYDAT release identity, source raster values, and the claimed comparison tolerances against the 5,965-basin historical aggregation: these require source archives beyond the deposited tables; network access was prohibited.
- The distinct **263-gage** count for provisional remarks was not separately tallied; the **15,257 affected rows** were verified.
- Signature formulas, annual-to-summary statistics, explorers, ZIP contents, and geometry validity/area calculations: these belong to companion-resource reviews. This resource contains none of those file types; companion checks here covered IDs, metadata, selection, and qualification.
- No files were changed. Each large input Parquet received one complete scan; the streamflow audit’s measured peak memory was approximately **564 MB**.

# Review A: numeric claims vs deposited data

Verdict: about 75 numeric claims checked; 40 MATCH, 29 MISMATCH, 6 NOT MEASURABLE. The mismatches fall into four groups. (1) The Resource 1/2 description is stale (100 signatures, 1,653 columns, 121 signatures, old annual row counts). (2) The Daymet description (6,087 / 6,041 basins, 97,757,220 records, 1980-2023, 7 variables) does not match the deposited file. (3) The MODIS date ranges in Table 3 and the Table 1 footnote (2021/2022) and the "5 schemes" statement disagree with the data and with Section 3. (4) "Watersheds smaller than 85,000 km2" and the 8,014 vs 7,964 polygon/gage counts are inconsistent with the files. Several internal two-value conflicts are listed at the end.

Abbreviations: R1/R2 = signature CSV/parquet for WY1993-2025 / WY1980-2025; R3 = harmonized inputs; R4 = gpkg + attributes; R5 = EO parquets. Gage ids compared zero-stripped.

| Location | Claim | Measured | File / method | Verdict |
|---|---|---|---|---|
| §2 para 1 | 16,994 candidate gages | 16,994 | R3 gage_metadata rows | MATCH |
| §2 / Table 1 | 9,154 USGS, 7,840 HYDAT candidates | 9,154 / 7,840 | metadata gage_type counts | MATCH |
| §2 / Table 1 | 8,014 gages with usable records (6,160 USGS; 1,854 HYDAT) | 8,014 (6,160 / 1,854); streamflow parquet has the same 8,014 ids | metadata processing_status=success; unique gage_id in streamflow parquet | MATCH |
| §2 para 1 | the 8,014 are "watersheds smaller than 85,000 km2" | 58 of the 8,014 have basin_area > 85,000 km2 (101,000-1,679,100); 59 by polygon area. 49 of them are inside BOTH signature products (R1 and R2: basin_area > 85,000 for 49 gages, all Canadian) | metadata basin_area; R1/R2 CSV basin_area; unsimplified gpkg geom_area_km2 | MISMATCH |
| §2.1.1 | filter "resulted in 7,964 gages across US and Canada" (vs 8,014 in §2, §3, Table 1) | 8,014 gages have flow data; 7,964 is the simplified-polygon / MODIS LAI site count. Of those 7,964, 7,960 are in the 8,014 and 4 are USGS ids not in it; 54 of the 8,014 (all Canada, > 85,000 km2 except 5) have no simplified polygon | R3 metadata vs R4 hisss_watershed_boundaries.gpkg | MISMATCH (internal conflict 8,014 vs 7,964; 7,964 is not a gage count) |
| §2 para 1 | >111 million observations | 111,624,189 rows | streamflow parquet metadata | MATCH |
| §2.1.1 / §3 R3 | 111.6 M / 111,624,189 records | 111,624,189 | pyarrow metadata | MATCH |
| §2.1.1 | daily Q 1 Oct 1979 to 30 Sep 2025 | min date 1979-10-01, max 2025-10-01; 5,381 rows after 2025-09-30 (1 Oct 2025) | streamflow parquet Date | MISMATCH (minor: runs one day past 30 Sep 2025) |
| §3 R3 | streamflow "1980 through 2025" | dates 1979-10-01 to 2025-10-01; 492,956 rows before 1 Jan 1980 | same | MISMATCH (minor; WY1980 starts Oct 1979) |
| §2.1.1 / §6.1 | area-normalized for all but 73 gages without published area | 73 (all Canadian) with area_normalized=False among 8,014; none USGS | metadata | MATCH |
| §6.1 | 73 of 8,014; 32 and 28 qualify for WY1993 / WY1980 products | 73; R1 32 (area_normalized False); R2 28 | R1/R2 CSV area_normalized | MATCH |
| §2.1.1 | "at least 20 years" | min num_water_years 20 (both products) | R1/R2 CSV | MATCH |
| §2.2.1 | gages need >= 20 WYs AND >= 60 % of window years | WY1993-2025: min fraction 0.606, OK. WY1980-2025 (46 yr window, 60 % = 27.6): min num_water_years = 20 (43.5 %); 325 gages have fewer than 28 valid WYs, 282 have a start-to-end span under 60 % | R2 CSV num_water_years, start/end_water_year | MISMATCH (WY1980-2025 only) |
| §2 / §2.2.1 | 6,678 gages (WY1993-2025), 6,250 (WY1980-2025) | 6,678 rows / 6,678 ids in annual parquet; 6,250 / 6,250 | R1, R2 CSV + parquet; selection txt files 6,678 / 6,250 | MATCH |
| Table 1 / §2 | Trends 1980-2025 = 6,250 sites; 1993-2025 = 6,678 | gages carrying summaries 6,250 / 6,678; gages with any non-NA trend slope in any family: 5,885 / 6,266 (climate family alone) | R1/R2 CSV | MATCH (as gage counts) |
| Fig. 3 caption | ecoregion counts sum to 6,250 and 6,678 | consistent with the two products | R1/R2 | MATCH |
| §2 para 1 | 121 signatures across nine families = 100 annual x 16 stats + 21 diagnostics | 139 annual signatures x 16 stats + 21 scalars = 160 signature names | categories csv (kind: annual 139, scalar 21); CSV has 16 stat columns for each of 139 | MISMATCH |
| §2.2.1 | 139 signatures | 139 annual (+ 21 scalars) | categories csv | MATCH (conflicts with 121 / 100 in §2 and §3: internal inconsistency) |
| §2 / §2.2.1 | nine signature families | 9 distinct categories (Climate 39, Flow Volume 27, Storage 25, Flashiness 19, Flow Timing 15, Drought 15, Snow 14, Flow Duration 3, Precip-Streamflow 2) + 1 row "(none: preprocessing diagnostic)"; 17 distinct function_family | hisss_signature_categories.csv | MATCH |
| §3 R1/R2 | 100 signatures x 16 statistics | 139 x 16 | R1/R2 CSV columns (16 columns for each of 139 annual signatures; 8 + 8 Pettitt; 1,112 Pettitt columns) | MISMATCH |
| §3 R1/R2 | 21 signatures without statistics | 21 scalar signatures present | categories csv + CSV | MATCH |
| §3 R1/R2 | 12 QA flags | 12 flagged_for_* columns | CSV | MATCH |
| §3 R1/R2 | 1,653 columns | 2,277 (20 metadata + 139x16 + 21 + 12 flags + ... = 2,277) | R1/R2 CSV shape | MISMATCH |
| §3 R1/R2 | annual table 18,898,406 and 24,366,487 rows | 27,766,578 (R1) and 34,813,881 (R2) | parquet metadata num_rows | MISMATCH |
| §3 R1/R2 | annual table columns gage_id, signature, water_year, value | exactly these | parquet schema | MATCH |
| §3 conventions | sixteen statistic suffixes | 16 per signature (8 + 8 Pettitt) | CSV | MATCH |
| §3 conventions | "elasticity_static" is single-valued | listed as scalar? (21 scalars present; name check not done) | categories | NOT MEASURABLE (not checked by name) |
| §3 R3 | gage metadata for all 16,994 candidates | 16,994 rows | metadata | MATCH |
| Table 4 R3 | 8,014 gages with daily streamflow | 8,014 | parquet ids | MATCH |
| Table 4 / §3 R3 | Daymet for 6,087 basins | 8,017 sites | hisss_daymet_basin_daily.parquet site_id unique over all batches | MISMATCH |
| §3 R3 | Daymet 97,757,220 records | 134,605,430 rows | parquet metadata | MISMATCH |
| §3 R3 / §2.1.3 | Daymet calendar years 1980-2023 | 1980-01-01 to 2025-12-31 | Date min/max | MISMATCH |
| §2.1.3 | seven daily variables incl. day length | 6 value columns: prcp, tmin, tmax, swe, vp, srad (no dayl) | parquet schema | MISMATCH |
| §3 R3 | variables precipitation, temperature, SWE, vapor pressure, shortwave | prcp, tmin, tmax, swe, vp, srad | schema | MATCH |
| Table 1 footnote | Daymet (6,041) sites; trends/summaries 6,041-7,964 | Daymet 8,017 sites. Gages with climate summaries: 6,678 (R1), 6,223 (R2); with a climate trend: 6,266 / 5,885. No file gives 6,041 | R3 parquet; R1/R2 CSV climate-family columns | MISMATCH (also conflicts with 6,087 in Table 4) |
| Table 1 | Catchment polygons 8,014 | simplified layer 7,964; unsimplified layer 8,017 (8,013 of the 8,014 processed gages + 4 USGS ids not processed) | R4 gpkg tables | MISMATCH |
| Table 1 / Table 4 | HydroATLAS 8,014 | attribute table 8,014 rows, 8,014 unique gage_id | R4 attributes csv/parquet | MATCH |
| Table 4 R4 | 7,964 boundary polygons (gpkg + GeoParquet) | 7,964 in both; WGS84 (EPSG 4326). Unsimplified gpkg (also in R4) has 8,017 | gpkg `watersheds` count; parquet num_rows; gpkg_geometry_columns | MATCH (Table does not mention the 8,017-polygon layer) |
| §2.1.2 | US official polygons cover 100 % of US gages | 6,160 of 6,160 processed USGS gages have a gagesii polygon (6,164 polygons in total) | unsimplified gpkg watershed_geom_source | MATCH |
| §2.1.2 | Canada 98.3 % covered by official polygons | 1,823 wsc_eccc / 1,854 = 98.33 % (30 hydrobasins, 1 gage without polygon). Simplified layer: 1,771 / 1,800 = 98.4 % | unsimplified gpkg | MATCH |
| §3 R4 | 211 watershed-scale attributes | attribute table has 211 columns, but 7 are identifiers / metadata (gage_id, latitude, longitude, basin_area, gage_type, Downstream_HB_ID, n_member_basins); 204 are attributes | attributes csv header | MISMATCH |
| §2.1.4 / Table 3 | ~210 attributes per watershed | 211 columns / ~204 attributes (see above) | same | MATCH (approximate) |
| Table 3 | 6 attribute categories | cannot be derived from column names without the dictionary | - | NOT MEASURABLE |
| §2.1.4 | HydroBASINS level 12 mean area 110 km2 | no file contains it | - | NOT MEASURABLE |
| §3 R5 | LAI 2002-2024, 2,150,280 rows, 7,964 watersheds x 270 months | 2,150,280 rows; 7,964 gages x 270 months each; 2002-07 to 2024-12 | LAI parquet | MATCH |
| Table 3 | MODIS LAI timespan 2022-2024; "record begins in 2022" | data span 2002-2024 | LAI parquet year | MISMATCH |
| §3 R5 | MODIS land cover 2001-2024, 191,136 rows | 191,136 rows = 7,964 x 24; years 2001-2024 | LULC parquet | MATCH |
| §3 R5 | 102 per-class percent columns, eight classification schemes | 102 `_pct` columns across 8 schemes (lc_type1-5, lc_prop1-3) | LULC parquet schema | MATCH |
| Table 3 | MODIS land cover timespan 2021-2024; "record begins in 2021" | 2001-2024 | LULC parquet | MISMATCH |
| Table 1 footnote | MODIS (7,964; 2021-2024) | 7,964 sites OK; years 2001-2024 (LULC) / 2002-2024 (LAI) | R5 parquets | MISMATCH (years) |
| Table 3 | 5 Land Cover Type schemes | 8 schemes in the deposited table (type1-5 + prop1-3); conflicts with "eight" in §3 | LULC parquet | MISMATCH (internal: 5 vs eight) |
| §3 R5 / §2.1.4 / Table 3 | NLCD 1985-2025; 250,879 rows; 6,119 CONUS watersheds | 250,879 = 6,119 x 41; years 1985-2025 | NLCD parquet | MATCH |
| Table 1 footnote | NLCD (6,119, 1985-2024) | 1985-2025 | NLCD parquet | MISMATCH (years) |
| §2.1.4 / §3 R5 | 16 NLCD classes + impervious | 16 `nlcd_c*_pct` + impervious_mean_pct | NLCD parquet | MATCH |
| §2.1.4 | NLCD extends record 16 years further back than MCD12Q1 | 2001 - 1985 = 16 | parquets | MATCH |
| §2.1.4 | NLCD only for 6,119 CONUS; Alaska/Canada MODIS-only | 6,119 = all 6,119 CONUS-US polygons with NLCD; 45 US (Alaska) excluded (6,164 US polygons) | NLCD parquet; excluded csv | MATCH |
| Table 3 | LAI valid range 0-10 | lai values within physical range (checked min/max columns) | LAI parquet | MATCH (data fall inside; not a count) |
| §4 | median watershed 698 km2, Q75 2,700, Q25 180 | R2 (WY1980-2025): median 697.55, Q75 2,718, Q25 183.5. R1 (WY1993-2025, Fig. 2 period): 644.65 / 2,500 / 165.75. All 8,014: 613 / 2,403 / 159 | basin_area in R1/R2 CSV; metadata | MATCH for WY1980-2025 only (rounded); MISMATCH if the 1993-2025 sites are meant |
| §4 | 12 Level I ecoregions; 62 % in two ecoregions; < 10 sites in three | no Level I ecoregion field in any file (HydroATLAS fec/cls are different classes) | - | NOT MEASURABLE |
| §4 | coverage statements (80 % of sites, lowest for storage) | not tested | - | NOT MEASURABLE |
| §2.2.2 | trend stats need >= 20 valid annual values (recession/elasticity exempt) | Non-exempt signatures: no trend statistic with < 20 annual values. 8,224 (R1) / 5,292 (R2) trend cells with < 20 values, all in alpha_linear, b_*, concavity, elasticity_*, log_a_* | annual parquet counts vs mk_pval | MATCH |
| §2.1.1 | gaps <= 3 days filled, > 30 NA years set to NA, ephemeral filter | streamflow parquet is the pre-screening series (394,210 NaN Q) | - | NOT MEASURABLE |
| §6.2 | example id 0103500 | USGS ids in the deposit are 8 digits zero-padded in the signature tables (e.g. 01035000), 7 digits zero-stripped in R3; 0103500 does not exist | R3 metadata / R1 CSV | MISMATCH (example id) |
| §6.2 | HYDAT example 01AD002 | present | metadata | MATCH |
| §6.2 | columns basin_area, geom_area_km2, Downstream_HB_ID exist | basin_area (R1/R2/R3/R4), geom_area_km2 (gpkg), Downstream_HB_ID (attributes) all present | headers | MATCH |
| §3 conventions | file name hiss_signatures_wy{window}.csv | deposited name is hisss_signatures_wy{window}.csv (extra s) | file listing | MISMATCH (typo) |
| §3 | five resources | R1-R5 as described | - | MATCH |

## Internal inconsistencies (same quantity, two values)
- Signatures: 100 x 16 + 21 = 121 (§2, §3) vs 139 (§2.2.1, §5) vs measured 139 annual + 21 scalar.
- Gage / polygon count: 8,014 (§2, Table 1, Table 4 R3/R4) vs 7,964 gages (§2.1.1) vs "6,041-7,964" (Table 1) vs measured 8,014 gages / 7,964 simplified polygons / 8,017 unsimplified polygons.
- Daymet basins: 6,041 (Table 1 footnote) vs 6,087 (Table 4, §3) vs measured 8,017; period 1980-2023 (§2.1.3, §3) vs 1980-2025 measured; the intro says 1980 through 2025.
- MODIS years: 2021-2024 / 2022-2024 (Table 3, Table 1 footnote) vs 2001-2024 / 2002-2024 (§3).
- NLCD end year: 1985-2024 (Table 1 footnote) vs 1985-2025 (Table 3, §2.1.4, §3).
- MODIS land-cover schemes: 5 (Table 3) vs eight (§3); data has 8.
- Attribute count: 210 (§2.1.4, Table 3) vs 211 (§3); measured 211 columns, 204 attributes.
- Median watershed 698 km2 is the WY1980-2025 value; Fig. 2 describes the 1993-2025 sites (measured 645).

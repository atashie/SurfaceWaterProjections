# Joining HISSS with other datasets — verified facts, audit findings, draft §5 text (2026-09-04)

Purpose: ground the manuscript's §5 Usage Notes item "How to join with other datasets
(CAMELS-Chem etc.)" in verified facts about (a) the partner datasets' identifiers and
(b) the actual join behaviour of the five staged HISSS HydroShare resources
(`~/Downloads/Signatures/`). Everything below was checked on 2026-09-04 against the
partner datasets' own files/docs (three research passes) and against the staged
deposit files (scripts run in-session; numbers reproducible from the staged files).

## 1. HISSS join keys — what the staged resources actually contain

| Table | Key column(s) | Form | Distinct ids |
|---|---|---|---|
| R1/R2 `hisss_signatures_wy*.csv`, `*_annual.parquet` | `gage_id` | zero-padded string | 6,678 / 6,250 |
| R3 `hisss_streamflow_daily.parquet` | `gage_id` | zero-padded (authoritative set) | 8,014 |
| R3 `hisss_daymet_basin_daily.parquet` | `site_id` | zero-padded string | 6,087 (122 never-compiled candidates) |
| R3 `hisss_gage_metadata.csv` | `gage_id` | **zero-STRIPPED** for US ids (7,046 seven-char, 72 nine-char, 9 fifteen-char) | 16,994 |
| R4 `hisss_watershed_boundaries.{gpkg,parquet}` + `_qa.csv` | `gage_id`, `canon_id` | padded + stripped | 7,964 |
| R4 `hisss_watershed_hydroatlas_attributes.*` | `gage_id`, `Downstream_HB_ID` | padded; HYBAS_ID present for **every** product gage | 8,014 |
| R5 LAI / LULC / NLCD parquets | `gage_id`, `canon_id` | padded + stripped | 7,964 / 7,964 / 6,119 |

US `gage_id` lengths in the products: 8 digits 5,432 (#1) / 5,009 (#2); **9–10 digits 59 / 31**;
15 digits 4 / 3. Canadian ids: 7-char `ddAAddd`, 1,183 / 1,207. `area_normalized = FALSE`:
32 / 28. Product metadata columns available for joins/filters: `latitude`, `longitude`
(WGS84), `basin_area` (agency-published km²), `gage_type` (USGS / Canada), `CLASS`
(GAGES-II Ref / Non-ref: 1,121 / 1,041 Ref), `RHBN` (245 / 247 true),
`human_interference_class`.

### Defect 1 — mis-padded `gage_id` in R5 (44 ids) and R4 boundaries (9 ids)
Root cause: `EO_data_processing/geometry/build_v3.py:25` and
`~/Downloads/geometry_rebuild/rebuild_watershed_polygons.py:137` —
`padded(c) = amap.get(c, c.zfill(8) if c.isdigit() else c)`. When a canonical
(stripped) id was not in the signature files used to build `amap`, the fallback
`zfill(8)` re-pads to 8 digits, which is wrong for every 9–10-digit USGS site number
that starts with 0 (stripped form has 8–9 digits, so `zfill(8)` is a no-op):
`011055566` → stored as `11055566`. `canon_id` is correct in all rows (0 internal
inconsistencies), so only the padded key is wrong.

Measured impact (direct `gage_id` join vs canonical join):

| Product | Table | direct | canonical | gages silently lost |
|---|---|---|---|---|
| WY 1993–2025 | LAI / LULC | 6,599 | 6,634 | **35** |
| WY 1993–2025 | NLCD | 5,419 | 5,454 | **35** |
| WY 1993–2025 | boundaries | 6,634 | 6,634 | 0 |
| WY 1980–2025 | LAI / LULC | 6,196 | 6,205 | **9** |
| WY 1980–2025 | NLCD | 4,998 | 5,007 | **9** |
| WY 1980–2025 | boundaries | 6,203 | 6,205 | **2** (`054310157`, `0214253830`) |

The R5 README's coverage sentence ("6,599 of 6,678 … 5,419 and 4,998") and its
"no identifier re-padding is needed" claim are therefore wrong by exactly these counts;
the remaining product gages without EO/polygon coverage are the 44 / 45 basins larger
than 100,000 km² that were excluded from delineation.

> **DECISION (user, 2026-09-04, given in the census session): the staged data files are NOT
> rewritten — they stay byte-identical to their validated builds. The join rule (read ids
> as text; strip leading zeros on both sides; never re-pad) is documented instead, in the
> root README (mirrored to CZ-Sync/HISSS) and the Resource 3/4/5 READMEs + dictionaries.
> The recommendation below is retained for the record; the §4 draft text must be revised
> accordingly (it assumes zero-padded ids everywhere).**

Fix (recommended, unambiguous — superseded by the decision above): rewrite `gage_id` in the three R5 parquets, the R4
boundary gpkg/parquet/QA csv (and the ids embedded in the three R5 explorers) from the
map canonical-id → padded id built on the authoritative R3 streamflow `gage_id` set
(every mis-padded id canonicalizes to exactly one streamflow id; the 4 inclusive-universe
extras `01591000/01591400/01591610/01591700` are already correctly padded). Then re-run
the census above (expect 0 lost), update the R5 README numbers, and fix the `zfill(8)`
fallback in both geometry scripts so a future rebuild cannot reintroduce it.

### Defect 2 — R3 README join recipe is wrong for 9–10-digit ids
`hisss_readme_inputs.md` recommends `s.zfill(8) if s.isdigit() and len(s) <= 8`
to re-pad the zero-stripped metadata ids. That recovers 6,621 / 6,678 (#1) and
6,221 / 6,250 (#2) product gages — **57 and 29 gages are lost**, all 9–10-digit ids with a
leading zero. Stripping leading zeros on BOTH sides (the repo's own
`canonical_gage_id` rule, `R/aggregate_hydroatlas_metadata.R:60`) matches 100 %.
Fix: replace the recipe (Python and R blocks) with strip-both-sides.

## 2. Partner datasets — verified identifier conventions

| Dataset | Key | Form / pitfall | Coverage | Source |
|---|---|---|---|---|
| CAMELS-Chem (Sterle et al. 2024) | `gauge_id` | **integer, zero-stripped** (`1013500`); re-pad `zfill(8)` (all CAMELS ids are 8-digit, so unambiguous) | 516 US catchments ⊂ CAMELS 671 ⊂ HCDN-2009 | hs.841f5e85085c423f889ac809c1bed4ac |
| CAMELS (Addor et al. 2017) | `gauge_id` | zero-padded 8-digit in the `camels_*.txt` files; basin shapefile `hru_id` is integer | 671 US | Zenodo 15529996 |
| HYSETS (Arsenault et al. 2020) | `Official_ID` (`Watershed_ID` is an internal index) | agency id; read as string | 14,425 US/CA/MX | OSF rpc3w |
| CAMELS-SPAT (Knoben et al. 2025) | `Station_id` (+ `Country`); files `USA_01013500`, `CAN_01AD002` | agency id, `{gauge_id:08}` | 1,426 US+CA | 10.20383/103.01306 |
| Caravan (Kratzert et al. 2023) | `gauge_id` = `{source}_{id}` (`camels_01013500`, `hysets_01AD003`, `hysets_0137449480`) | strip the prefix; HYSETS entries use the agency id | US/CA/MX via HYSETS | Zenodo 10968468 |
| CANOPEX (Arsenault et al. 2016) | HYDAT station number | column names UNVERIFIED (site unreachable) | 698 CA | — |
| GAGES-II (Falcone 2011/2017) | `STAID` (attributes), `GAGE_ID` (shapefiles, C15) | zero-padded text | 9,322 US | 10.5066/P96CPHOT |
| Water Quality Portal / NWIS | `MonitoringLocationIdentifier` = `USGS-01013500`; dataRetrieval `site_no` character | prefix | US | waterqualitydata.us |
| USGS NLDI → NHDPlus V2 `comid` (StreamCat key) | `/linked-data/nwissite/USGS-{site}` | live: USGS-01013500 → comid 724696; **AK/HI/Canada return null** | CONUS only | api.water.usgs.gov/nldi |
| HydroBASINS/HydroATLAS | `HYBAS_ID` (10-digit int) | HISSS R4 attributes carry it as `Downstream_HB_ID` for all 8,014 gages | global | HydroSHEDS |
| RHBN | HYDAT station number | 1,027 stations | CA | ECCC |
| MacroSheds (Vlah et al. 2023) | `network`/`domain`/`site_code` | **not gage-keyed**; one `usgs` site (Black Earth Creek = NWIS 05406457, not among HISSS candidates); no NWIS crosswalk column; join spatially | 383 sites / 224 stream gauges, US (+Sweden, Antarctica); no Canada | figshare per domain; R pkg 2.0.3 |
| Daymet VPD (Corak et al. 2025) | gridded only | HydroShare hs.de74b0a457c74deca09f9a41afa03c8f, 154 GB; 23 annual CONUS+ NetCDF mosaics **2001–2023** (IGBP and KG variants); `VPD` ushort in Pa (×1e-3 = kPa); Daymet V4 1 km LCC grid and 365-day calendar, but **no CRS/x/y variables and swapped dim names** — reconstruct before gdptools/exactextract | CONUS+ | HydroShare |
| NADP NTN / SNOTEL | `siteID` (`NC41`) / station triplet | station-keyed → spatial join | US | — |

## 3. Overlap measured against the staged products

- **CAMELS-Chem**: 507 / 516 catchments in WY 1993–2025; **515 / 516** in WY 1980–2025;
  all overlapping gages are GAGES-II `CLASS = Ref`. The 9 absent from #1 mostly ended
  their record ~2009–2011 (fail the 60 % window rule); 1 is `processing` status.
  CAMELS: 657 / 671 and 669 / 671. Chemistry = 76,284 dated grab samples 1980–2018
  (wide, 24 constituents, `measure_unit_code` varies); NADP deposition per gauge-year
  1985–2018 (calendar vs water year unstated); paired `q_daily` in cfs; no boundaries
  shipped (CAMELS uses Geospatial Fabric + GAGES-II areas: `area_gages2` matches
  HISSS `basin_area` source for US gages).
- **MacroSheds** (224 stream gauges with coordinates): 146 lie inside ≥ 1 WY 1993–2025
  product watershed (144 for WY 1980–2025), but nested — median 3 containing polygons,
  smallest containing polygon median **799 km²** — so the relation is almost always
  *nesting*, not equivalence. Co-located (≤ 500 m from a product gage): 11 / 9
  (Baltimore LTER ×5, HJ Andrews ×2–3, Sleepers River ×2, Loch Vale, Santa Barbara MC06).
  Discharge L/s daily; `ms_status`/`ms_interp` QA flags.
- **Polygons**: 6,634 / 6,678 and 6,203 (6,205 after Defect 1) / 6,250 product gages have a
  R4 polygon; the rest are the > 100,000 km² exclusions (agency polygons must be
  fetched directly for those). Daymet aggregation used the 6,087 basins < 85,000 km².

## 4. Draft §5.1.3 text (first draft, delivered 2026-09-04; assumes Defect 1 is fixed)

> 5.1.3 Joining HISSS with other datasets. All HISSS tables share one key, gage_id: the
> agency station identifier stored as a zero-padded character string (USGS site number,
> e.g. 01013500; HYDAT station number, e.g. 01AD002). Read it as text — numeric parsing
> silently drops the leading zero carried by most US sites. Any dataset keyed to the same
> agency identifiers joins directly once its identifiers are normalized to this form:
> CAMELS and CAMELS-Chem (Addor et al. 2017; Sterle et al. 2024) store the USGS site
> number as an integer that must be re-padded to eight digits; HYSETS (Official_ID),
> CAMELS-SPAT (USA_/CAN_ prefix) and Caravan (camels_/hysets_ prefix) carry the same
> identifier under other names; the Water Quality Portal uses USGS-{site}; and the USGS
> Network-Linked Data Index maps CONUS sites to the NHDPlus COMID used by StreamCat and
> other reach-based products. Of the 516 CAMELS-Chem catchments, 507 are in the
> WY 1993–2025 product and 515 in WY 1980–2025. Datasets that are not gage-keyed join
> spatially through the watershed polygons of Resource 4: gridded products such as the
> Daymet-derived VPD of Corak et al. (2025), which shares Daymet's 1 km grid and 365-day
> calendar, can be aggregated with the same area-weighted workflow used for Daymet
> (Sect. 2.1.3), and point records such as MacroSheds sites (Vlah et al. 2023) can be
> located within the (nested) polygons or matched to the nearest gage by coordinates.
> HydroATLAS-based products join on the HydroBASINS level-12 outlet identifier
> (Downstream_HB_ID) carried in the Resource 4 attribute table. When combining data,
> aggregate partner records to the water year (October–September, labeled by the ending
> year); re-normalize fluxes to a common drainage area, since HISSS reports both the
> agency-published area used for normalization (basin_area) and the polygon area
> (geom_area_km2); and pair record-dependent signatures (Sect. 5.1.2) with the product
> whose window matches the partner record.

References the paragraph introduces (check the list): Arsenault et al. 2020 (HYSETS,
already cited in §1); Knoben et al. 2025 (CAMELS-SPAT, cited in §1); Kratzert et al.
2023 (Caravan, cited in §1); Hill et al. 2016 (StreamCat) if named; USGS NLDI
(api.water.usgs.gov/nldi) as a URL; Sterle et al. 2025 HydroShare data DOI
10.4211/hs.841f5e85085c423f889ac809c1bed4ac if the data (not only the paper) is cited.

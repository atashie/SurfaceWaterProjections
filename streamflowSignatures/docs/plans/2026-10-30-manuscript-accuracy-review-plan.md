# Manuscript + guidelines accuracy review — plan (2026-10-30)

**Request (user, 2026-10-30):** assess the HISSS manuscript and the accompanying signatures
document for (a) methodological accuracy against what the code actually did and (b) every
"number of X" claim against the finalized, deposited datasets. A sibling session is checking
the HydroShare resources themselves; this review uses those files as the ground truth and
does not re-review them.

**Status:** plan drafted; execution waits for the numbered decisions in §6.

## 0. Scope, ground truth, assumptions

| Item | Choice |
|---|---|
| Manuscript | the live Google Doc, synced 2026-10-30 (`docs/MANUSCRIPT_DRAFT.md`, 227 paragraphs; a major revision since 2026-10-21 — Tables 1/3/4 now in the text, §5 Technical Validation new, references 1–43) |
| "Accompanying signatures document" | ASSUMED to be the co-authors' guidelines doc "Summary Documentation for Streamflow Signatures" (`docs/SIGNATURE_GUIDELINES.md`, synced 2026-10-30 — a pure reorder since 2026-10-21: 3.3 Flow Timing … 3.7 Climate Metrics … 3.9 Storage). The HydroShare readmes / dictionaries are the sibling session's scope (decision 1) |
| Ground truth for NUMBERS | the files AS UPLOADED to HydroShare: the staged folders `~/Downloads/hisss_R{1,2,3}_upload_2026-10-29/`, `hisss_R4_upload_final_2026-10-30/`, `hisss_vegetation_and_land_cover_upload_2026-10-30/` (+ the four unchanged R4/R5 tables of the August staging) — MD5-matched to the live resources 2026-10-29/30; the R1 summary CSV is byte-identical to `~/HISSS_data/qualifier-rerun-26oct2026/processedOuts_1993_2025_climate_v3_26oct2026/…_signatures.csv` (md5 5211f8bc…). **Every number is RECOMPUTED from these files by script; nothing is copied from a readme.** |
| Ground truth for METHODS | `julia/src/` + `config/signatures_config.json` (canonical), the guidelines doc, `docs/SIGNATURES.md`, `docs/DEVELOPMENT.md`, `EO_data_processing/README*.md`, `docs/QUALIFIER_RERUN_RESULTS.md` |
| Prior findings reused, not re-derived | `docs/plans/2026-09-04-filtering-stage-census.md` (funnel counts; its gage-inclusion logic reproduced both products exactly), `docs/reconciliation/manuscript_log.md` (2026-09-29 → 2026-10-28 items), `docs/plans/2026-10-12-coauthor-relay.md` (open relay list) |
| Out of scope | prose quality, citation formatting beyond wrong-paper / wrong-product mappings, figure aesthetics |

## 1. Claim ledger — every "number of X" and every method statement

Columns: **ID** · where · claim as written (2026-10-30 text) · what the deposited files / code say ·
how it is verified · **status**. Status key: ✔ consistent (still re-measured in Phase 1);
✘ wrong or stale (replacement known); ≈ imprecise wording; ? not yet measured; ⚠ needs a
user decision. "readme" = the deposited resource readme of 2026-10-29/30; "census" = the
2026-09-04 census; "config" = `signatures_config.json`.

### 1.1 Abstract and §1

| ID | Where | Claim | Deposited / code | Verify by | Status |
|---|---|---|---|---|---|
| A1 | Abstract | "annual snow, climate, and streamflow data for 8,014 watersheds" | 8,014 = compiled DAILY streamflow gages; annual signatures exist for 6,678 / 6,250; Daymet covers 8,013 of them (8,017 basins) | R3 metadata + R1/R2 CSV gage counts | ≈ wording: "daily streamflow and climate for 8,014 (8,013) … annual signatures for 6,678 / 6,250" |
| A2 | Abstract | "trends and physiographic attributes for 6,041–7,964 locations depending on the subset" | 6,041 traced to no file (2026-09-18); trends 6,250 / 6,678; attributes 8,014; polygons 7,964 | counts in Phase 1 | ✘ → "6,250–8,014" (decision 6) |
| B1 | §1 | ">6,000 streamgages", "over 6,000 watersheds" | 6,678 / 6,250 | R1/R2 CSV | ✔ |
| B2 | §1 | "spanning the years 1980 through 2025" | WY 1980–2025 product; climate 1980–2025 | R2 CSV, R3 Daymet | ✔ |
| B3 | §1 | "roughly 100 hydrological signatures" | 139 annual (+ 21 scalars = 160) | R1 dictionary / categories file | ✘ → "139" |
| B4 | §1 | "~1 km Daymet", "annually resolved MODIS and NLCD LULC, monthly MODIS LAI", "Julia, Python, and R" | as stated | — | ✔ |

### 1.2 §2 preamble and Table 1

| ID | Where | Claim | Deposited / code | Verify by | Status |
|---|---|---|---|---|---|
| C1 | §2 | 16,994 candidates (9,154 USGS; 7,840 HYDAT) | census: 16,994 (9,154 / 7,840) | R3 metadata rows × gage_type | ✔ |
| C2 | §2 | 8,014 usable (6,160 USGS; 1,854 HYDAT) | census: 8,014 (6,160 / 1,854) | R3 streamflow parquet distinct gage_id × metadata | ✔ |
| C3 | §2 | "usable daily records for watersheds smaller than 85,000 km²" | the 8,014 are NOT size-filtered (58 exceed 85,000 km² by published area; the boundary layer's cut was 100,000 km², 54 dropped) | R3 metadata basin_area; R4 QA csv | ✘ carried (2026-09-29 item 3) |
| C4 | §2, §2.1.1, §3 R3 | "more than 111 million" / "111.6 million" / "111,624,189 observations" | readme: 111,624,189 rows after the qualifier patch (restored days were NA ROWS, so the count is unchanged — confirm) | pyarrow row count of the deposited streamflow parquet | ✔? re-measure |
| C5 | §2 | "121 signatures across nine signature families … 100 annually resolved … plus 21 per-gage diagnostics" | 139 annual + 21 scalars = 160; nine categories; 2,277 columns | categories file (139 annual / 21 scalar), dictionary (2,277) | ✘ → "160 signatures (139 annually resolved + 21 per-gage diagnostics)" — must match §2.2.1's 139 |
| C6 | §2 | windows 1980–2025 (6,250) / 1993–2025 (6,678) "(Sect. 2.1.2)" | counts ✔; the windows are defined in §2.2.1 | — | ✔ numbers; ✘ cross-ref → 2.2.1 (also "(Sect. 2.1.1)" for boundaries → 2.1.2) |
| T1 | Table 1 | Streamflow NWIS 6,160 / HYDAT 1,854, "1980-2025" | ✔; Canadian records end 2024-12-31 (no WY 2025) — footnote | R3 parquet max date by gage_type | ✔ (+ footnote) |
| T2 | Table 1 | "Cathment Polygons … 8,014" | simplified layer 7,964; unsimplified 8,017 (universe 8,018 incl. 4 geometry-only gages) | R4 gpkg feature counts (sqlite3) | ✘ → 7,964 (+ 8,017 full-resolution) — decision 4 on which to cite; typo "Cathment" |
| T3 | Table 1 | HydroATLAS 8,014 | 8,014 × 211 | R4 attribute CSV | ✔ |
| T4 | Table 1 | signatures row "6,041 – 7,964*"; footnote "MODIS (7,964; 2021-2024) and NLCD (6,119, 1985-2024); Climatic data from DAYMET (6,041)" | signatures 6,678 / 6,250; MODIS LULC 2001–2024, LAI Jul 2002–Dec 2024; NLCD 1985–2025; Daymet 8,017 basins (8,013 compiled gages) | R5 parquets year ranges; R3 Daymet distinct site_id | ✘ every cell; "Caculated" typo; "Daymet21,26" cites HYDAT (21) → 26 |
| T5 | Table 1 | Trends 6,250 / 6,678; footnote "at least X% of data in the first and last decade" | 80 % first/last decade + 60 % overall + ≥ 20 values (config) | config | ✔ counts; ✘ "X%" → 80 % |

### 1.3 §2.1 Input data

| ID | Where | Claim | Deposited / code | Verify by | Status |
|---|---|---|---|---|---|
| D1 | §2.1.1 | retrieval 1979-10-01 to 2025-09-30, n = 16,994, dataRetrieval / tidyhydat | requests ran to 2025-10-01 (one WY 2026 day present, rejected) | R3 parquet date range | ✔ (note) |
| D2 | §2.1.1 | three QA stages: retrieval, compilation, per-year | matches `R/` ingestion + `preprocess_daily_data` | docs | ✔ |
| D3 | §2.1.1 | "perennial defined as >30 days of streamflow data over 0.0001 mm/day" + citation 23 (GAGES-II) | code: a calendar year counts when ≥ 30 days have Q > 0.0001; ≥ 20 such years at compilation | `R/` compile rule; census §10 | ≈ "at least 30 days"; the GAGES-II citation for this rule needs the authors' confirmation |
| D4 | §2.1.1 | "polygons smaller than 85,000 km² … resulted in 7,964 gages" | 8,014 from the streamflow rule; 7,964 = boundary layer after the 100,000 km² cut; 6 polygons 87,150–97,182 km² remain | R4 QA csv areas | ✘ carried (2026-10-21 item 1): "…yielding 8,014 gages; 7,964 of these have a watershed polygon ≤ 100,000 km²" |
| D5 | §2.1.1 | "retained original agency quality flags, including ice-affected and estimated values" | ✔, but the QUALIFIER RULE sentence is missing (relay 2026-10-28 item 1): A/P status + finite value kept; remark codes never mask; February compile had masked `A R`/`A <`/`A >`/`A [4]` — 193,542 days on 987 gages re-fetched 2026-10-26 and restored | R3 ledger row count / gages | ✘ add sentence (counts re-measured from the ledger) |
| D6 | §2.1.1 | missing = NA, zero valid, "all but 73 gages" area-normalized | readme: 73; negatives kept as served (21,518 rows, 82 gages) | R3 metadata area_normalized; parquet Q < 0 | ✔ |
| D7 | §2.1.1 | "gaps of up to three consecutive days filled by linear interpolation" | config `interpolation.max_gap_days` 3, linear, internal only | config | ✔ |
| D8 | §2.1.1 | "any WY with more than 30 missing values set to NA" | year rejected (config `max_raw_na_per_year` 30) | config | ✔ |
| D9 | §2.1.1 | "if any data gap exceeded three consecutive days those days were set to NA" | the WHOLE water year is rejected (`reject_residual_na` true; `io.jl` l. 411–414) | code | ✘ carried → "the water year was excluded" |
| D10 | §2.1.1 | seasons Winter Jan–Mar, Spring Apr–Jun, Summer Jul–Sep, Fall Oct–Dec | config: winter [12,1,2], spring [3,4,5], summer [6,7,8], fall [9,10,11] | config | ✘ carried → DJF / MAM / JJA / SON |
| E1 | §2.1.2 | GAGES-II covers 100 % of US gages; WSC polygons 98.3 % of Canadian | census: 6,164 / 6,164; 1,823 / 1,854 | R4 gpkg `watershed_geom_source` counts | ✔ |
| E2 | §2.1.2 | HydroBASINS level-12 fallback | 29 (simplified) / 30 (unsimplified) | R4 gpkg | ✔ (counts to add if wanted) |
| E3 | §2.1.2 | "We excluded basins that exceeded 85,000 km²" | simplified layer excludes 53 basins > 100,000 km² + 05KH009; the UNSIMPLIFIED layer (deposited 2026-10-30) KEEPS the 53 | R4 gpkg counts / areas | ✘ → 100,000 km², and say both layers exist |
| E4 | §2.1.2 | "boundaries were used to aggregate climate and landcover data only where we had sufficient streamflow data coverage for trends analyses" | Daymet aggregated for ALL 8,017 basins; EO products for all 7,964 polygons — not restricted to trend gages | R3 Daymet, R5 tables | ✘ → delete the restriction |
| F1 | §2.1.3 | Daymet "1 km²", variables tmin, tmax, prcp, vp, srad, swe, **day length** | deposited table: prcp, tmin, tmax, swe, vp, srad (six; no dayl) | R3 parquet schema | ✘ drop "day length" |
| F2 | §2.1.3 | "covering 1980–2023" | 1980–2025 | R3 parquet year range | ✘ carried → 1980–2025 |
| F3 | §2.1.3 | "watershed size threshold … largest HUC8 … 85,000 km²" (for Daymet) | no size cap in the deposited Daymet (8,017 basins incl. the 53 > 100,000 km²) | R3 Daymet flags / R4 | ✘ sentence now moot |
| F4 | §2.1.3 | "gdptools … area-weighted zonal aggregation" | the deposited table is the project's own fractional cell-coverage (exactextract) weighting over the full-resolution polygons; gdptools = the co-authors' original 5,965-basin run, which the new build reproduces (daily diffs ≤ 0.016 mm prcp) | `EO_data_processing/README_DAYMET.md` | ✘ carried → rewrite; keep gdptools as the validation reference |
| F5 | §2.1.3 | (missing) | coverage: 8,017 basins = 8,013 of 8,014 compiled gages + 4; 365-day calendar; no missing days; 101 flagged basins | R3 flags csv | ≈ add one sentence |
| G1 | §2.1.4 | BasinATLAS v10; level-12 "mean area of 110 km²" | readme: "roughly 100 km²" | R4 attribute table / HydroBASINS docs (SUB_AREA) | ? measure the level-12 mean over the domain |
| G2 | §2.1.4 | aggregation: continuous/percentage as area-weighted means, elevation min/max, categorical area-weighted majority | readme: HYBRID — 91 `_u`/`_p` columns pass through from the outlet subbasin; `_s` columns re-aggregated (area-weighted means; elevation min/max; argmax of outlet `_u` fractions for glc/pnv; SUB_AREA-weighted mode for nine others) | R4 dictionary `aggregation` column counts | ≈ relay: add the pass-through rule |
| G3 | §2.1.4 / Table 3 / §3 R4 | "approximately 210 attributes" / "210" / "211 watershed-scale attributes" | 211 columns = 198 attributes + 13 identifiers/diagnostics | R4 CSV header | ≈ → "198 attributes (211 columns)" |
| G4 | Table 3 | MCD12Q1 "5 Land Cover Type schemes" | 8 schemes (5 types + 3 FAO-LCCS property layers), 102 columns — as §3 R5 says | R5 LULC parquet column prefixes | ✘ internal contradiction → 8 |
| G5 | Table 3 | MCD12Q1 "2021-2024 / Data record begins in 2021"; LAI "2022-2024 / begins in 2022" | LULC 2001–2024 (24 years); LAI Jul 2002–Dec 2024 (270 months) | R5 parquets | ✘ both |
| G6 | Table 3 | "Hydro ATLAS22" | ref 22 = Lehner 2022 (lakes); HydroATLAS = Linke 2019 (25); HydroBASINS = Lehner & Grill 2013 (34) | reference list | ✘ citation |
| G7 | Table 3 / §2.1.4 | NLCD 16 classes, 30 m, 1985–2025, CONUS, 6,119; "16 years further back"; Alaska/Canada MODIS-only | ✔ (2001 − 1985 = 16; 45 Alaska gages excluded) | R5 NLCD parquet, excluded csv | ✔ |
| G8 | §2.1.4 | "(Table 3)" while no Table 2 exists | tables are 1, 3, 4 | — | ✘ numbering |

### 1.4 §2.2 Metrics

| ID | Where | Claim | Deposited / code | Verify by | Status |
|---|---|---|---|---|---|
| H1 | §2.2.1 | "139 streamflow and hydroclimate signatures … including signature-specific diagnostics and data-quality flags" | 139 annual; + 21 scalars; + 12 flags | categories / dictionary | ✔ (but contradicts C5) |
| H2 | §2.2.1 | nine families: flow volume, flow duration, storage, flashiness, drought, flow timing, precipitation-to-streamflow, snow, climate | = `docs/signature_categories.csv` (23/3/17/16/10/15/2/14/39 annual) | categories file | ✔ |
| H3 | §2.2.1 | `hisss_data_dictionary.csv` with Resources 1 and 2 | ✔ | R1/R2 folders | ✔ |
| H4 | §2.2.1 | ≥ 20 qualifying WYs and ≥ 60 % of the window → 6,678 / 6,250 | census reproduced both exactly | `qualification_census.jl` (re-run on the corrected input only if cheap — the gage sets were FROZEN, so the rule no longer reproduces them: 31 / 29 would-qualify gages) | ✔ counts; ⚠ add the frozen-cohort sentence (relay 2026-10-28 item 2) |
| H5 | §2.2.1 | "Water year and seasonal metrics required at least 80 % of raw observations passing QA" | SEASONAL only (config `seasonal_completeness` 0.8 on raw observations; seasonal flow reversals not masked); the water year uses the ≤ 30-missing-day rule | config, `io.jl` l. 336–339 | ≈ → "seasonal metrics" |
| H6 | §2.2.1 | 1993–2025 "the date range with the most gages with the same data record length" | census: defensible (2,909 of 6,678 complete 33-year records) | re-measure on v3 annual parquet | ✔ |
| I1 | §2.2.2 | years removed: > 3 consecutive NAs or > 30 NAs | config (+ residual boundary NA, 838 of 33,732 at the census) | config | ✔ |
| I2 | §2.2.2 | flagged (not removed): negative flow; ≥ 15 non-zero days at one constant value in a month | config `constant_sd_flag` 15 / 1; `reject_negative_flow` false; `Negative_ann` signature | WHERE these surface in the products (a column? the annual parquet?) — ? | ? verify surfacing |
| I3 | §2.2.2 | trend gate ≥ 20 values, ≥ 60 % series, ≥ 80 % first/last decade; recession + elasticity exempt | config `trend_completeness` 0.6 / 0.8; `stats_floor`; exemptions | config; non-NA counts `_mean` vs `_senn_slp` in the v3 CSVs | ✔ |
| I4 | §2.2.2 | (missing) Pettitt: fixed WY 1980–2024 window, ≥ 20 values with ≥ 10 each side; 20-year floor counted IN-WINDOW | readme / config `changepoint` | config | ≈ add (carried in-window item) |

### 1.5 §3 Data Record

| ID | Where | Claim | Deposited / code | Verify by | Status |
|---|---|---|---|---|---|
| J1 | §3 | five resources, CC-BY 4.0, README + dictionary each, CSV/Parquet | sibling session | — | ✔ (not re-checked here) |
| J2 | Table 4 R1/R2 | contents "run provenance and validation reports", formats "JSON" | DECISION 2026-10-29: seven files each (CSV, Parquet, dictionary, categories, readme, HTML, ZIP); no provenance / validation files deposited | upload notes | ✘ → drop; add categories file; formats CSV, Parquet, HTML, ZIP, MD |
| J3 | Table 4 R3 | 8,014 streamflow; 16,994 metadata; Daymet "6,087" | 8,017 basins | R3 Daymet | ✘ 6,087 → 8,017 |
| J4 | Table 4 R4 | 8,014 attributes; 7,964 polygons | + unsimplified 8,017-polygon layer (DECISION 2026-10-30) | R4 | ≈ add |
| J5 | Table 4 R5 | 7,964 MODIS; 6,119 NLCD | ✔ | R5 | ✔ |
| J6 | §3 R1–R2 | "1,653 columns … 100 signatures × 16 … 21 … 12 flags" | 2,277 = 20 metadata + 139 × 16 + 21 + 12 | R1/R2 CSV header | ✘ |
| J7 | §3 R1–R2 | annual rows "18,898,406 and 24,366,487" | 27,766,578 / 34,813,881 | R1/R2 parquet metadata | ✘ |
| J8 | §3 R1–R2 | "(iii) … explorer … with fitted checksums, configuration, and software versions) and validation reports" | garbled; explorer = HTML + sidecar ZIP; nothing else deposited | upload notes | ✘ rewrite (iii); add (iv) dictionary + categories |
| J9 | §3 R1–R2 | "hiss_signatures_wy{window}.csv"; "Section 2.1.2" | `hisss_…`; §2.2.1 | — | ✘ typo + cross-ref |
| J10 | §3 R3 | "metadata table (16,994) records … human-interference indicators from GAGES-II and RHBN / regulation via HYDAT" | metadata CSV columns: id, lat/lon, basin_area, gage_type, num_years/days, processing_status, area_normalized, Downstream_HB_ID — NO GAGES-II columns (those 20 metadata columns are in the R1/R2 summary tables); RHBN / REGULATED are the separate `hisss_hydat_interference.csv` (8,012 stations) | R3 CSV headers | ✘ wording |
| J11 | §3 R3 | Daymet "6,087 basins (97,757,220 records, calendar years 1980-2023)", variables "precipitation, temperature, SWE, vapor pressure, shortwave radiation" | 8,017 basins, 134,605,430 rows, 1980–2025, six variables | R3 Daymet parquet | ✘ |
| J12 | §3 R3 | (missing) five companion files | ledger (193,542 rows), Daymet flags (101), two gage lists (6,678 / 6,250), would-qualify list (31 / 29) — DECISION 2026-10-29 | R3 folder | ≈ add |
| J13 | §3 R4 | 7,964 polygons GeoPackage + GeoParquet WGS84; `watershed_geom_source`; "211 attributes"; "BasinAtlas v10 22,25" | ✔ / ✔ / 198 (+13) / ref 22 wrong; + unsimplified layer, QA csv, dictionary | R4 | ≈ |
| J14 | §3 R5 | LAI 2002–2024, 2,150,280 = 7,964 × 270; "Myneni et al. 2015" | ✔ arithmetic; author-year cite in a numbered list; ref 27 is the MOD15A2H DOI, product is MCD15A3H | R5 LAI parquet | ✔ counts; ✘ citation (carried) |
| J15 | §3 R5 | LULC 2001–2024, 191,136 rows, 102 columns, "all eight classification schemes" | ✔ (7,964 × 24) | R5 LULC parquet | ✔ (contradicts Table 3's five) |
| J16 | §3 R5 | NLCD 1985–2025, 250,879 rows, 6,119, 16 classes | ✔ (6,119 × 41) | R5 NLCD parquet | ✔ |
| J17 | §3 R5 | "granule-level provenance manifests for the MODIS products" | NOT deposited (readme: manifest not retained; provenance in the code repo) | R5 upload note | ✘ drop; add the Alaska exclusion list |
| J18 | §3 Conventions | join on zero-stripped ids; canon_id in geometry/land-cover tables; 16 suffix names; single-valued examples; NaN ≡ missing row; `flagged_for_` | all ✔ against the readmes; "section 5.1.1" → 6.1; duplicate canon_id sentence; stray "}" in `nlcd_c{code}_pct}` | dictionary | ✔ / ✘ cross-ref |

### 1.6 §4–§6, availability, acknowledgements, references

| ID | Where | Claim | Deposited / code | Verify by | Status |
|---|---|---|---|---|---|
| K1 | §4 / Fig. 2 | 12 Level I ecoregions; 62 % in Eastern Temperate Forests + Northwestern Forested Mountains; < 10 sites in Southern Semiarid Highlands, Tundra, Hudson Plain; per-ecoregion counts | never verified locally (needs the CEC/EPA Level I layer) | spatial join of the 6,678 / 6,250 gage coordinates to the Level I polygons (decision 2) | ? |
| K2 | §4 | median 698 km², 25 % > 2,700, 25 % < 180 km² | = product #2's `basin_area` quartiles (184 / 698 / 2,718; product #1: 166 / 645 / 2,500) while Figure 2 shows the 1993–2025 set | re-measure on both v3 CSVs | ⚠ name the product (decision 5) |
| K3 | §4 / Fig. 3 | coverage lowest for storage; most sites > 80 %; lowest in Taiga / Hudson Plain; Fig. 3a per-family signature counts; Fig. 3b counts sum to 6,250 / 6,678 | coverage CHANGES with the v3 rerun (619 / 502 gages regained 1,177 / 950 years); counts per category are now 23/3/17/16/10/15/2/14/39 | per-gage, per-category metric-years / possible metric-years from the v3 annual parquets (decision 3) | ? |
| L1 | §5 | "all 139 signatures"; three-level audit | ✔; the cross-language full-scale validation (Python + R ports), the change census and the blind reviews are not mentioned | — | ≈ optional addition |
| M1 | §6.1 | 73 of 8,014 un-normalized; 32 / 28 in the products | readme 73; 32 / 28 (the canonical full-record run has 37 — products hold 32 / 28) | `area_normalized` counts in R3 metadata and both v3 CSVs | ✔ re-measure |
| M2 | §6.1 | precipitation-dependent signatures "(runoff ratios, elasticity, Q-P seasonality, storage) are structurally NA" | readme repeats it; does the STORAGE family really go NA for `area_normalized = FALSE` (BFI needs no PPT)? | NA pattern of `BFI_*`, `avg_storage*`, recession columns on the 32 gages in the v3 CSV | ? verify |
| M3 | §6.1 | unit-carrying list: flow volumes, percentiles, recession log(a), drought deficits | ✔ (readme) | — | ✔ |
| N1 | §6.2 | "e.g. 0103500" | 7 digits; USGS ids ≥ 8 (e.g. 01030500) | — | ✘ carried |
| N2 | §6.2 | "may be read as a character string" | must | — | ✘ carried |
| N3 | §6.2 | `Downstream_HB_ID` "in the Resource 4 attribute table"; `basin_area`; `geom_area_km2` | Downstream_HB_ID is in the R3 metadata CSV — is it also in the R4 attribute table? basin_area (R1/R2), geom_area_km2 (R4) ✔ | R4 CSV header | ? |
| N4 | §6.2 / §6.3 | "(Sect. 5.1.3)"; "compared against the WY 1993–2025 and WY 1980–2025 products" | 6.3; "between" | — | ✘ carried |
| O1 | §6.3 | record-dependent list: `*_all` pulses, elasticity, parameterized BFI, all drought metrics | = readme | — | ✔ |
| P1 | Data / code availability | HydroShare URL; GitHub mirror URL; `docs/SIGNATURES.md` path; "licenceCC-BY 4.0 license" | mirror snapshot pending the next publish; typo | `gh`/curl after publish | ✔ / ✘ typo |
| P2 | Acknowledgements | "Claude Code 0.145.0" | the project ran on many releases (local CLI today: 2.1.314) | — | ⚠ decision 7 |
| Q1 | References | 22 (Lehner 2022, lakes) cited for HydroBASINS (Table 1) and BasinAtlas (§3 R4) | 34 (Lehner & Grill 2013) / 25 (Linke 2019) | — | ✘ |
| Q2 | References | 27 Myneni = MOD15A2H DOI; product MCD15A3H v061 | DOI 10.5067/MODIS/MCD15A3H.061 | — | ✘ carried |
| Q3 | References | Table 1 "Daymet21,26" | 26 only | — | ✘ |
| Q4 | References | 23 (GAGES-II) cited for the ephemeral-stream definition (§2.1.1) | authors to confirm the source of the > 30-day rule | — | ? ask |
| Q5 | References | 29 Anderson 1976 (classification system) for NLCD; 36 gdptools undated; 29/30/36 without years | fine as secondary; formatting | — | ≈ |

### 1.7 Guidelines doc (signatures document)

| ID | Where | Claim | Code / deposit | Verify by | Status |
|---|---|---|---|---|---|
| R1 | Part 1.1 (sources) | 16,994 / 9,154 / 7,840; "watersheds smaller than 85,000 km² and at least 20 total years" → 8,014 (6,160 / 1,854); "> 111 million"; "Albers et al., 2026" | same conflation as C3 / D4; Albers 2017 | as C1–C4 | ≈ / ✘ |
| R2 | Part 1.2–1.3 | 30 NA days; 3-day gap; 60 % / 80 % decades / 20 values; seasonal 80 %; season months | config | script: extract every number in Parts 1–3 and print the config value beside it | ✔ expected; re-verify mechanically |
| R3 | Part 2 | 16 statistics, names | dictionary suffixes | dictionary | ✔ |
| R4 | §3.1–3.9 | every threshold (0.0001 mm/day; 90th/10th pulse percentiles; 2 % reversal; 11-year window; 15 years; 0.1 / 10 / 1 mm; SWE 10 mm; 60-day spell; 7-day MA; p2–p30; 0.555 bimodality; R² 0.8 …) | config sections | same script | ? mechanical pass |
| R5 | §3.6 Snow | "fourteen metrics … for water years 1980–2023" | 14 ✔; years 1980–2025 | — | ✘ → 1980–2025 |
| R6 | §3.7 Climate | 39 metrics (20 prcp + 19 temp); names and conventions | implemented; open relay items 2–9, 14–15 of 2026-10-12 | — | carried |
| R7 | Part 4 | "USGS values are kept only for qualifiers A, A e, P, P e (others set to NA)" | the HIGH defect; rule now A/P + finite value, remark codes carried | `R/usgs_qualifiers.R` | ✘ DUE (relay) |
| R8 | Part 5 | 12 automated flags | dictionary's 12 `flagged_for_*` names | dictionary | ? names match |
| R9 | Part 6 | "Resource Collections" working notes | — | — | ✘ carried (relay 12) |
| R10 | numbering | duplicate "1.1"; the new 3.3–3.9 order vs every repo note that says "§3.9 Climate" | repo docs: STATUS, guidelines_todos, CHANGELOG, SIGNATURES.md? | grep "§3.9\|3.9 Climate" | repo-side chore after the review |

## 2. Verification tooling (Phase 1)

One script, `docs/benchmarks/manuscript_claim_census.py` (new; `docs/benchmarks/` holds tools —
rule `benchmarks.md`), run with the repo `.venv` (pandas + pyarrow; sqlite3 for the
GeoPackages — no geopandas is installed). It takes the five staged folders as arguments and
writes `docs/plans/2026-10-30-manuscript-claim-ledger.csv` (ID, claim, measured value, file,
method) so every number in §1 has a reproducible origin. Measurements:

| Files | Measured |
|---|---|
| R1 / R2 summary CSVs | gage count; column count and composition (20 + 139 × 16 + 21 + 12); US/CA split; `area_normalized` FALSE; `flagged_for_*` counts; `basin_area` quartiles; `_mean` vs `_senn_slp` vs `_pettitt_cp_year` non-NA counts per base; NA pattern of the storage family on the un-normalized gages (M2); gages with complete 33- / 46-year records (H6) |
| R1 / R2 annual parquets | row count; distinct signatures (139); water-year range; per-gage × per-category coverage = metric-years / (signatures × valid years) for Figure 3 (K3, if decision 3) |
| R1 dictionary + categories | 2,277 rows; category counts; the 21 scalars; the 12 flags (R8) |
| R3 streamflow parquet | rows (C4), distinct gages, date range by agency (T1), negative rows, NA rows |
| R3 metadata CSV | 16,994 rows; `processing_status` × `gage_type`; `area_normalized` FALSE (73 compiled); column list (J10, N3); ≥ 85,000 km² count among compiled (C3) |
| R3 Daymet parquet (5.1 GB) | row count from file metadata; schema (F1); distinct `site_id` and year range from a column-only scan |
| R3 ledger / flags / lists | 193,542 rows & 987 gages; 101; 6,678 / 6,250; 31 / 29 |
| R4 GeoPackages (sqlite3) | feature counts 7,964 / 8,017; `watershed_geom_source` counts (E1–E2); areas > 85,000 / > 100,000 km² (D4, E3); attribute CSV 8,014 × 211, `aggregation` counts in the dictionary (G2), `Downstream_HB_ID` presence (N3); level-12 mean SUB_AREA if available (G1) |
| R5 parquets | rows / distinct gages / year or month ranges (G5, J14–J16); 102 LULC columns by scheme prefix (G4); 16 NLCD classes; 45 excluded |
| Ecoregions (decision 2) | download the CEC North America Level I layer; point-in-polygon of gage coordinates (pure-Python ray casting or `R sf` — `sf` availability to be checked) → counts per ecoregion for both products, the 62 % and the < 10 claims |

Runtime: minutes for everything but the Daymet column scan (~5 min) and the two annual
parquets (~2 min each). Memory: read columns, not tables; stagger with the sibling session
if it is loading the same parquets.

## 3. Methodology review (Phase 2)

For each method statement (D2, D3, D5, D7–D10, E4, F4, G2, H4–H5, I1–I4, M2 and the
guidelines' Parts 1–3) the check is: manuscript sentence → guidelines sentence → the
config knob or the code line → the deposited product behaves that way (one measured
witness where cheap, e.g. a gage-year with a 4-day gap absent from the annual parquet for
D9; a seasonal `Qwin` NA with 79 % raw coverage for H5). New since the last pass and not
yet verified: G2 (HydroATLAS hybrid aggregation vs the manuscript's three-rule
description), I2 (where the negative / constant-value flags surface in the products),
M2 (storage NA for un-normalized gages), F1 (six Daymet variables, no day length), E4
(aggregation not restricted to trend gages). The guidelines doc gets the mechanical
number-vs-config pass of R2 / R4 (script output appended to the ledger).

## 4. Review and outputs (Phases 3–4)

1. `docs/plans/2026-10-30-manuscript-claim-ledger.md` — §1 with the measured column filled,
   every ✔ / ✘ backed by the CSV; the `?` rows resolved.
2. Consolidated co-author relay note `docs/reconciliation/2026-10-30-coauthor-relay.md`:
   one line per item, "replace X with Y" (user preference 2026-09-04), folding the open
   2026-10-12 list and the 2026-10-28 items; Tables 1 / 3 / 4 as corrected rows (decision 4).
   Grouped: manuscript numbers · manuscript methods · guidelines doc · citations.
3. Blind reviews before the relay (decision 8, the close-out workflow): reviewer A re-derives
   20 randomly drawn numbers from the deposited files without seeing the ledger; reviewer B
   reads the Methods sections against `julia/src/` + config without seeing §1; Codex on the
   ledger as a whole. Findings resolved, ledger re-issued.
4. Logs: `manuscript_log.md` and `guidelines_todos.md` entries (started 2026-10-30), STATUS
   one-liner, CHANGELOG only if code or repo docs change (R10 renumbering chore). Commit →
   push → `./publish_to_hisss.sh` with the pending 2026-10-30 batch.

## 5. Order and effort

| Phase | Work | Est. |
|---|---|---|
| 0 | sync, snapshots, this plan (DONE 2026-10-30) | — |
| 1 | census script + run; ecoregion join if decided | 2–3 h |
| 2 | methods checks with witnesses; guidelines number pass | 1–2 h |
| 3 | blind reviews (A numbers, B methods, Codex) and resolution | 1–2 h |
| 4 | relay note, logs, STATUS, commit / push / publish | 1 h |

## 6. Decisions wanted (recommended default first)

1. **Scope of "signatures document"**: (a) the guidelines Google Doc only [recommended — the
   HydroShare readmes / dictionaries are the sibling session's]; (b) also the deposited
   `hisss_data_dictionary.csv` definitions (139 + 21 rows against `docs/SIGNATURES.md`).
2. **Ecoregion claims (K1)**: (a) verify — download the CEC Level I layer and join the gage
   coordinates [recommended; ~30 min; the only §4 numbers never checked]; (b) leave to the
   figure's author and mark "not verified".
3. **Figure 3 coverage (K3)**: (a) recompute per-category coverage from the v3 annual parquets
   and hand the authors the numbers behind a regenerated figure [recommended — the figure
   must be redrawn anyway after the rerun]; (b) only flag that it must be regenerated.
4. **Relay format**: (a) one-line replace items + the three tables as corrected ROWS
   [recommended — cell-by-cell replace items for a 20-cell table are unreadable]; (b) one-line
   items only.
5. **§4 basin-size quartiles (K2)**: (a) report both products' quartiles and let the authors
   pick [recommended]; (b) replace with the 1993–2025 set's (Figure 2's window).
6. **Abstract range (A2)**: (a) "6,250–8,014 locations" (trend products → compiled gages)
   [recommended]; (b) "6,250–7,964" (trend products → polygons); (c) authors' call, flag only.
7. **"Claude Code 0.145.0" (P2)**: (a) suggest "Claude Code (Anthropic), releases 1.x–2.x,
   2025–2026" [recommended]; (b) leave as the authors wrote it.
8. **Blind reviews**: (a) two Sonnet reviewers + Codex as in §4.3 [recommended, the standing
   close-out workflow]; (b) Codex only; (c) none.

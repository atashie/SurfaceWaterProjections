# Manuscript + guidelines claim ledger — measured (2026-10-30)

The filled version of §1 of `docs/plans/2026-10-30-manuscript-accuracy-review-plan.md`.
Every number in the **Measured** column comes from `docs/benchmarks/manuscript_claim_census.py`
run on the HydroShare files AS DEPOSITED (the staged upload folders, MD5-matched to the live
resources 2026-10-29/30); the raw measurements, 205 rows keyed by the ids in brackets, are in
`docs/plans/2026-10-30-manuscript-claim-ledger.csv`. Nothing is copied from a readme.
Decisions applied (user, 2026-10-30): the guidelines doc AND the deposited dictionary are in
scope (1a+b); ecoregions (K1), Figure 3 coverage (K3) and the §4 quartile wording (K2) are OUT
of scope; relay = one-line items + corrected table rows (4a); "Claude Code 0.145.0" gets a
suggested wording (7a); blind reviews A / B / C + Codex (8a). Item 6 (abstract range) is
measured and flagged, not decided.

Status key: ✔ consistent · ✘ wrong or stale · ≈ imprecise · ⚠ open · — out of scope.

## 1. Abstract and §1

| ID | Where | Claim | Measured | Status |
|---|---|---|---|---|
| A1 | Abstract | "annual snow, climate, and streamflow data for 8,014 watersheds" | 8,014 = compiled DAILY streamflow gages [r3.compiled]; Daymet covers 8,013 of them [r3.daymet_vs_compiled]; ANNUAL signatures exist for 6,678 / 6,250 [wy*.gages] | ≈ "daily streamflow and climate for 8,014 (8,013) watersheds, annual signatures and trends for 6,678 / 6,250" |
| A2 | Abstract, Table 1 | "trends and physiographic attributes for 6,041–7,964 locations" | nothing in the deposit yields 6,041 (also found 2026-09-18). Candidates: trends 6,250 / 6,678; attributes 8,014 [r4.attr_rows_cols]; polygons 7,964 / 8,017 [r4.*.features]; product gages with Daymet AND MODIS 6,634 / 6,205, either product 7,112 [r5.sig_daymet_modis] | ✘ flagged for the authors (item 6 not decided); the defensible range is 6,250–8,014 |
| B1 | §1 | "> 6,000 streamgages"; "over 6,000 watersheds" | 6,678 / 6,250 | ✔ |
| B2 | §1 | "1980 through 2025" | WY 1980–2025 product; Daymet 1980-01-01 – 2025-12-31 [r3.daymet_dates] | ✔ |
| B3 | §1 | "roughly 100 hydrological signatures" | 139 annual signatures (+ 21 per-gage scalars) [wy*.bases, .scalars] | ✘ → 139 |
| B4 | §1 | ~1 km Daymet; annual MODIS + NLCD LULC; monthly LAI; Julia / Python / R | as stated | ✔ |

## 2. §2 preamble and Table 1

| ID | Where | Claim | Measured | Status |
|---|---|---|---|---|
| C1 | §2 | 16,994 candidates (9,154 USGS; 7,840 HYDAT) | 16,994 (9,154 / 7,840) [r3.meta_by_type] | ✔ |
| C2 | §2 | 8,014 usable (6,160 USGS; 1,854 HYDAT) | 8,014 (6,160 / 1,854) in the metadata AND the streamflow table [r3.compiled, r3.q_gages] | ✔ |
| C3 | §2 | "usable daily records for watersheds smaller than 85,000 km²" | 58 compiled gages have a published area > 85,000 km² (52 > 100,000; 73 have none) [r3.compiled_gt85k] — the 8,014 are not size-filtered | ✘ carried |
| C4 | §2, §2.1.1, §3 R3 | "> 111 million" / "111.6 million" / "111,624,189 observations" | 111,624,189 rows (metadata and scan agree) [r3.q_rows_sum], of which 394,210 carry no discharge value → 111,229,979 observations [r3.q_nonmissing]; the qualifier patch filled NA rows, so the row count did not move | ✔ "more than 111 million"; ≈ call 111.6 M "daily records" (Codex) |
| C5 | §2 | "121 signatures across nine families … 100 annually resolved … plus 21 per-gage diagnostics" | 139 annual × 16 = 2,224 columns + 21 scalars + 12 flags + 20 metadata = 2,277 [wy*.composition]; nine categories [categories_by_cat] | ✘ → "160 signatures (139 annually resolved, each with 16 statistics, plus 21 per-gage diagnostics) across nine categories" — must agree with §2.2.1's 139 |
| C6 | §2 | windows "(Sect. 2.1.2)"; boundaries "(Sect. 2.1.1)" | windows are defined in §2.2.1, boundaries in §2.1.2 | ✘ cross-references |
| T1 | Table 1 | NWIS 6,160 / HYDAT 1,854; "1980-2025" | 6,160 / 1,854; streamflow 1979-10-01 – 2025-10-01; Canadian records end 2024-12-31 (no WY 2025 at any Canadian gage) [r3.q_dates_by_agency] | ✔ + footnote the Canadian end date |
| T2 | Table 1 | "Cathment Polygons … 8,014" | simplified layer 7,964; unsimplified 8,017 [r4.simplified.features, r4.unsimplified.features] | ✘ → 7,964 (8,017 at full resolution); typo |
| T3 | Table 1 | HydroATLAS 8,014 | 8,014 × 211 [r4.attr_rows_cols] | ✔ |
| T4 | Table 1 | "6,041 – 7,964*"; footnote "MODIS (7,964; 2021-2024) and NLCD (6,119, 1985-2024); DAYMET (6,041)" | signatures 6,678 / 6,250; MODIS LULC 7,964 gages, 2001–2024 [r5.lulc_rows]; LAI 7,964, Jul 2002 – Dec 2024 [r5.lai_months]; NLCD 6,119, 1985–2025 [r5.nlcd_rows]; Daymet 8,017 basins [r3.daymet_sites] | ✘ every cell; "Caculated" typo; "Daymet21,26" cites HYDAT (21) |
| T5 | Table 1 | Trends 6,250 / 6,678; "at least X% of data in the first and last decade" | 6,250 / 6,678; the gate is 80 % of each of the first and last decades, 60 % overall, ≥ 20 values (config `trend_completeness`, `stats_floor`) | ✔ counts; ✘ "X%" → 80 % |

## 3. §2.1 Input data

| ID | Where | Claim | Measured | Status |
|---|---|---|---|---|
| D1 | §2.1.1 | retrieval 1979-10-01 to 2025-09-30; n = 16,994; §3 R3 "1980 through 2025" | table spans 1979-10-01 (start of WY 1980) – 2025-10-01; 5,381 US gages carry the single 2025-10-01 day (rejected downstream) [r3.q_dates*] | ✔ (minor: "water years 1980–2025, 1 Oct 1979 – 1 Oct 2025") |
| D2 | §2.1.1 | three QA stages | retrieval → compilation (`R/`) → per-year screening (`preprocess_daily_data`) | ✔ |
| D3 | §2.1.1 | "perennial defined as > 30 days of streamflow data over 0.0001 mm/day" + ref 23 | compile rule: a CALENDAR year counts when ≥ 30 days have Q > 0.0001; ≥ 20 such years (config `filtering`; guidelines Part 4; metadata `num_years` 20–47 [r3.meta_years]) | ≈ "at least 30 days"; ref 23 (GAGES-II) as the source of this rule needs the authors' confirmation |
| D4 | §2.1.1 | "polygons smaller than 85,000 km² … resulted in 7,964 gages" | 8,014 gages from the streamflow rule; 7,964 = simplified boundary layer after the 100,000 km² cut; 6 polygons of 87,150–97,182 km² remain in it [r4.simplified.areas] | ✘ carried → "yielding 8,014 gages; the simplified boundary layer holds 7,964 polygons ≤ 100,000 km², covering 7,960 of them plus four stations that were never compiled (01591000 / 400 / 610 / 700)" [r4.attr_overlap_layers] (Codex) |
| D5 | §2.1.1 | "retained original agency quality flags, including ice-affected and estimated values" | true, but the qualifier RULE is unstated: A/P status + finite value kept, remark codes carried; the February compile had masked them — ledger 193,542 days on 987 gages (`A R` 181,461, `A <` 11,444, `A [4]` 419, `A >` 218; 96 provisional today), re-fetched 2026-10-26 [r3.ledger_*] | ≈ true but incomplete → add one sentence (relay 2026-10-28 item 1); the rule also drops the −999999 no-data sentinel (`R/usgs_qualifiers.R` l. 39–47, Codex) |
| D6 | §2.1.1 | missing = NA, zero valid, "all but 73 gages" | 73 compiled gages `area_normalized = FALSE`, all Canadian [r3.compiled_unnormalized]; NA stored as 394,210 null rows; negatives kept: 21,518 rows on 82 gages [r3.q_na_rows] | ✔ |
| D7 | §2.1.1 | gaps ≤ 3 days linearly interpolated | config `interpolation` 3 / linear / internal | ✔ |
| D8 | §2.1.1 | WY with > 30 missing values "set to NA" | the water year is rejected (`io.jl` l. 363) | ✔ (wording: excluded) |
| D9 | §2.1.1 | gap > 3 days: "those days were set to NA" | the WHOLE water year is rejected (`io.jl` l. 367 `gap_too_long`, l. 414 `residual_na`) | ✘ carried → "the water year was excluded" |
| D10 | §2.1.1 | seasons Jan–Mar / Apr–Jun / Jul–Sep / Oct–Dec | config: winter 12,1,2 · spring 3,4,5 · summer 6,7,8 · fall 9,10,11 | ✘ carried → DJF / MAM / JJA / SON |
| E1 | §2.1.2 | GAGES-II 100 % of US gages; WSC polygons 98.3 % of Canadian | gagesii 6,164 = all US; wsc_eccc 1,823 of 1,853 Canadian = 98.4 % (unsimplified), 1,771 of 1,800 = 98.4 % (simplified) [r4.*.sources] | ✔ (98.3 % was 1,823 / 1,854 incl. 05KH009) |
| E2 | §2.1.2 | HydroBASINS level-12 fallback | 29 (simplified) / 30 (unsimplified) | ✔ |
| E3 | §2.1.2 | "We excluded basins that exceeded 85,000 km²" | simplified layer: 0 polygons > 100,000 km², max 97,182; unsimplified: 53 > 100,000 (max 1,696,282) and 59 > 85,000 [r4.*.areas] | ✘ → "the simplified layer excludes the 53 basins > 100,000 km²; the full-resolution layer keeps them" |
| E4 | §2.1.2 | "used to aggregate climate and landcover data only where we had sufficient streamflow data coverage for trends" | Daymet aggregated for 8,017 basins [r3.daymet_sites]; LAI / LULC for 7,964 [r5.lai_gages]; NLCD 6,119 — far beyond the 6,678 / 6,250 trend gages | ✘ delete the restriction |
| F1 | §2.1.3 | Daymet variables incl. "day length" | deposited: prcp, tmin, tmax, swe, vp, srad — six, no dayl [r3.daymet_schema] | ✘ drop "day length" |
| F2 | §2.1.3 | "covering 1980–2023" | 1980-01-01 – 2025-12-31, 16,790 rows per basin, no NaN [r3.daymet_dates, r3.daymet_nan_prcp] | ✘ carried → 1980–2025 |
| F3 | §2.1.3 | size threshold "85,000 km² … largest HUC8" for Daymet | no size cap: the 53 basins > 100,000 km² are in the table [r4.unsimplified.areas] | ✘ sentence moot |
| F4 | §2.1.3 | "gdptools … area-weighted zonal aggregation" | deposited table = the project's fractional cell-coverage (exactextract) weighting over the full-resolution polygons (`EO_data_processing/README_DAYMET.md`); gdptools was the co-authors' original 5,965-basin run, reproduced | ✘ carried → rewrite; keep gdptools as the validation reference |
| F5 | §2.1.3 | (missing) coverage | 8,017 basins = 8,013 of 8,014 compiled gages (05KH009 without) + 4 never-compiled [r3.daymet_vs_compiled]; 101 flagged basins (57 low-confidence outlines, 45 few cells) [r3.flags_rows] | ≈ add |
| G1 | §2.1.4 | level-12 "mean area of 110 km²" | HISSS domain: Σ SUB_AREA / Σ member subbasins = 131 km²; median per-gage mean 139 km² [r4.level12_mean_area] | ⚠ not verifiable from the deposit: 110 km² describes the global level-12 tier (Codex); the 131 km² here is weighted by gage membership of nested basins, not a domain mean — ask for the source, do not replace |
| G2 | §2.1.4 | aggregation: area-weighted means / elevation min–max / categorical majority | dictionary: 91 columns PASS THROUGH from the outlet subbasin (HydroATLAS `_u` / `_p`), 92 area-weighted means, 9 SUB_AREA-weighted modes, 2 argmax of outlet class fractions, 2 areal fractions, 1 min, 1 max, 13 keys [r4.attr_dictionary] | ≈ add the pass-through rule |
| G3 | §2.1.4, Table 3, §3 R4 | "approximately 210" / "210" / "211 attributes" | 211 columns = 198 attributes + 13 identifiers / diagnostics | ≈ → "198 attributes (211 columns)" |
| G4 | Table 3 | MCD12Q1 "5 Land Cover Type schemes" | 8 schemes: lc_type1–5 + lc_prop1–3, 102 class columns [r5.lulc_schemes] — as §3 R5 says | ✘ → 8 |
| G5 | Table 3 | MCD12Q1 "2021-2024 / begins in 2021"; LAI "2022-2024 / begins in 2022" | LULC 2001–2024 (24 years); LAI Jul 2002 – Dec 2024 (270 months) | ✘ both |
| G6 | Table 3 | "Hydro ATLAS22" | ref 22 = Lehner 2022 (lakes); HydroATLAS = Linke 2019 (25) | ✘ |
| G7 | Table 3 / §2.1.4 | NLCD 16 classes, 30 m, 1985–2025, CONUS, 6,119; "16 years further back"; Alaska / Canada MODIS-only; "MCD12Q1 supplying consistent coverage across all HISSS watersheds" | 6,119 gages × 41 years = 250,879; 16 class columns; C1V2; 45 Alaska gages excluded; 6,119 + 45 = 6,164 US [r5.nlcd_*, r5.excluded]; MODIS covers the 7,964 polygons only — 44 / 45 signature-product gages (the large basins) have no land-cover record [r5.lulc_product_overlap] | ✔ NLCD; ✘ "all HISSS watersheds" (Codex) |
| G8 | §2.1.4 | "(Table 3)" | tables are numbered 1, 3, 4 | ✘ numbering |

## 4. §2.2 Metrics

| ID | Where | Claim | Measured | Status |
|---|---|---|---|---|
| H1 | §2.2.1 | "139 streamflow and hydroclimate signatures … diagnostics and data-quality flags" | 139 annual [wy*.bases]; 21 scalars; 12 flags | ✔ (C5 must match it) |
| H2 | §2.2.1 | nine families listed | = `docs/signature_categories.csv` (deposited copy identical on signature / kind / category [categories_vs_canonical]): Climate 39, Flow Volume 23, Storage 17, Flashiness 16, Flow Timing 15, Snow 14, Drought 10, Flow Duration 3, Precipitation Streamflow 2 — note the co-author sheet puts the runoff ratios and avg_storage under Storage and elasticity under Flashiness | ✔ |
| H3 | §2.2.1 | `hisss_data_dictionary.csv` with Resources 1 and 2 | present, 2,277 rows in column order [wy*.dictionary_matches_csv] | ✔ |
| H4 | §2.2.1 | "retained at least 20 qualifying WYs and spanned at least 60 % of the possible WYs in the analysis window" | runner (`run_julia_benchmark.jl` l. 259–269): possible years = window START to the gage's LAST water year WITH DATA in the window (the raw record, l. 261–269 — not the exported `end_water_year`, which is the last QUALIFYING year), not the window length. Products: min 20 valid years; minimum fraction under that rule 0.606 / 0.600 (Codex: 10287655 has 27 qualifying years and a last raw year 2024 → 27 / 45); by window length 325 gages of product #2 fall below 60 % (min 0.435, e.g. 02083800: 20 years 1980–2001) [wy*.num_water_years, min_fraction] | ≈ state the denominator: "at least 60 % of the water years from the window start to the gage's last year of record in the window" (the deposited readmes say "of the window's years" — WRONG for 325 gages of the 1980–2025 product; tell the sibling session) + the frozen-cohort sentence (31 / 29 would-qualify gages listed, not added [r3.would_qualify]) |
| H5 | §2.2.1 | "Water year and seasonal metrics required at least 80 % of raw observations passing QA" | SEASONAL only (config `seasonal_completeness` 0.8 raw; `io.jl` l. 336–339; applied to Qwin/Qspr/Qsum/Qfal and the four seasonal runoff ratios, not to the seasonal reversals). The water year uses the ≤ 30-missing-day / ≤ 3-day-gap rule. In BOTH products the seasonal rule masked NO valid water year (every Qwin … Qfal row finite: 202,011 / 255,145 = the Qann rows [wy*.annual_seasonal_rows]); the `season_excluded_years_*` diagnostics (Σ 3,532 / 3,300 / 2,766 / 3,938 in #1) count seasons of ALL water years including rejected ones (`io.jl` l. 442 pushes before `continue`) | ≈ → name the eight metrics (four seasonal totals, four seasonal runoff ratios); the diagnostic's meaning belongs in the dictionary (Codex) |
| H6 | §2.2.1 | 1993–2025 = "the date range with the most gages with the same data record length" | 3,267 of 6,678 carry the complete 33-year record (4,730 ≥ 30 years); 2,629 of 6,250 the complete 46-year record [wy*.complete_records] — two windows cannot show an optimum, and only half the gages share the full length | ⚠ not verifiable as stated (Codex); suggest "a shorter common window that retains more qualifying gages" |
| I1 | §2.2.2 | years removed: > 3 consecutive NAs or > 30 NAs | config + `io.jl` l. 363 / 367 (+ residual boundary NA l. 414) | ✔ |
| I2 | §2.2.2 | "flagged (but did not remove)" negative flow / constant monthly value | both are preprocessor diagnostics (`io.jl` l. 342–353, `diag_rows`); negative days surface as the signature `negative_ann` (16 columns); the constant-value flag is NOT exported in the products [wy*.negative_constant_cols] | ≈ say the count of negative-flow days is a signature and the constant-value check is internal (or export it) |
| I3 | §2.2.2 | trend gate ≥ 20 values, ≥ 60 % series, ≥ 80 % first/last decade; recession + elasticity exempt | 130 of 139 bases carry fewer trend than mean values; the 9 exempt = 7 recession + 2 elasticity bases [wy*.trend_gate_bases]; Qann: 6,678 means / 6,266 trends / 6,607 Pettitt (#1), 6,250 / 5,907 / 6,250 (#2) | ✔ |
| I5 | §2.2.2 | "trend statistics were computed only if its annual series had at least 20 valid values, at least 60% of the entire series … 80% of both the first and last decades" | the 20-value floor blanks ALL 16 statistics, not only the trends (`stats.jl` l. 395–401; config `stats_floor`; guidelines Part 1.3); the 60 % / 80 % gates are relative to the signature's own first-to-last valid year (`stats.jl` l. 415–443), not the window; recession + elasticity exempt from both | ≈ (found by reviewer B) |
| I6 | §2.1.3 / §2.2.1 | (missing) climate join and climate QA | each gage's Daymet series enters on Daymet's own calendar inside the record span (`climate_frame.jl`, default `outer` since 2026-10-28), so a day absent from the streamflow record keeps its climate values; PPT / SWE / temperature years screened with the same ≤ 30-day / ≤ 3-day rule (config `*_na_policy`) | ≈ add (relay 2026-10-28 item 3; found again by reviewer B) |
| I4 | §2.2.2 | (missing) Pettitt window; in-window 20-year floor | config `changepoint` 1980–2024, ≥ 20 obs, ≥ 10 per segment → 71 gages of product #1 have no Pettitt fields; the 20 valid years are counted inside the window (H4) | ≈ add |

## 5. §3 Data Record

| ID | Where | Claim | Measured | Status |
|---|---|---|---|---|
| J1 | §3 | five resources, CC-BY 4.0, README + dictionary each | sibling session's scope | ✔ (not re-checked) |
| J2 | Table 4 R1 / R2 | "run provenance and validation reports"; formats "JSON" | seven files each: CSV, Parquet, dictionary, categories, readme, explorer HTML + ZIP (decision 2026-10-29) | ✘ → drop; formats CSV, Parquet, HTML, ZIP, MD |
| J3 | Table 4 R3 | Daymet "6,087" | 8,017 basins | ✘ |
| J4 | Table 4 R4 | 8,014 attributes; 7,964 polygons | + the full-resolution layer, 8,017 polygons (decision 2026-10-30) | ≈ add |
| J5 | Table 4 R5 | 7,964 MODIS; 6,119 NLCD | 7,964 / 6,119 | ✔ |
| J6 | §3 R1–R2 | "1,653 columns … 100 signatures × 16 … 21 … 12 flags" | 2,277 = 20 + 139 × 16 + 21 + 12 | ✘ |
| J7 | §3 R1–R2 | "18,898,406 and 24,366,487 rows" | 27,766,578 / 34,813,881 [wy*.annual_rows] | ✘ |
| J8 | §3 R1–R2 | "(iii) … explorer … with fitted checksums, configuration, and software versions) and validation reports" | garbled; the explorer is an HTML + a sidecar ZIP; no provenance or validation files are deposited | ✘ rewrite; add the dictionary + categories file as (iv) |
| J9 | §3 R1–R2 | "hiss_signatures_wy{window}.csv"; "Section 2.1.2" | `hisss_…`; §2.2.1 | ✘ |
| J10 | §3 R3 | metadata table "records … human-interference indicators from GAGES-II … RHBN and regulation status via HYDAT" | metadata columns: gage_id, latitude, longitude, basin_area, gage_type, num_years, start_year, end_year, num_days, processing_status, error_message, area_normalized, source_dir, Downstream_HB_ID [r3.meta_cols] — the 11 GAGES-II / HYDAT indicators are the R1 / R2 summary tables' metadata; RHBN / REGULATED are `hisss_hydat_interference.csv` (8,012 stations; RHBN 316; REGULATED 3,116 / 4,525 / 371 empty) [r3.interference] | ✘ wording |
| J11 | §3 R3 | Daymet "6,087 basins (97,757,220 records, 1980-2023)"; five variables | 8,017 basins; 134,605,430 rows; 1980–2025; six variables (prcp, tmin, tmax, swe, vp, srad) | ✘ |
| J12 | §3 R3 | (missing) companion files | ledger 193,542 rows / 987 gages; Daymet flags 101; gage lists 6,678 / 6,250 (identical to the product gage sets [r3.list_*]); would-qualify 31 / 29 | ≈ add |
| J13 | §3 R4 | 7,964 polygons GeoPackage + GeoParquet; `watershed_geom_source`; "211 attributes"; "BasinAtlas v10 22,25" | 7,964 / 13 columns [r4.geoparquet]; source column present; 198 + 13; ref 22 wrong; + unsimplified layer, QA csv, dictionary | ≈ |
| J14 | §3 R5 | LAI 2002–2024, 2,150,280 rows = 7,964 × 270; "Myneni et al. 2015" | 2,150,280; 7,964; 270 months Jul 2002 – Dec 2024 [r5.lai_*] | ✔ counts; ✘ author-year cite (ref 27 is the MOD15A2H DOI; product MCD15A3H) |
| J15 | §3 R5 | LULC 2001–2024, 191,136 rows, 102 columns, eight schemes | 191,136 = 7,964 × 24; 102; 8 | ✔ (contradicts Table 3) |
| J16 | §3 R5 | NLCD 1985–2025, 250,879 rows, 6,119, 16 classes | 250,879 = 6,119 × 41; 16; C1V2 | ✔ |
| J17 | §3 R5 | "granule-level provenance manifests for the MODIS products" | not deposited (manifest not retained) | ✘ drop; mention the Alaska exclusion list |
| J18 | §3 Conventions | join rule, canon_id, 16 suffixes, scalar examples, NaN ≡ missing row, flagged_for_; "All tables join on gage_id" | ✔ against the files except the key name: the Daymet table's key is `site_id` [r3.daymet_schema] (Codex); "section 5.1.1" → 6.1; duplicate canon_id sentence; stray "}" | ✔ / ✘ |

## 6. §4–§6, availability, acknowledgements, references

| ID | Where | Claim | Measured | Status |
|---|---|---|---|---|
| K1 | §4 / Fig. 2 | 12 Level I ecoregions; 62 %; < 10 sites in three | out of scope (decision 2) | — not verified |
| K2 | §4 | median 698 km², 25 % > 2,700, 25 % < 180 km² | measured anyway: product #2 184 / 698 / 2,718 (n = 6,222 with area); product #1 166 / 645 / 2,500 (n = 6,646) [wy*.basin_area_quartiles] — the text matches the 1980–2025 product while Figure 2 shows the 1993–2025 set | — (decision 5: no recommendation; the measurement is here) |
| K3 | §4 / Fig. 3 | coverage statements; per-family counts | out of scope (decision 3); the per-category counts of H2 are what Figure 3a needs | — |
| L1 | §5 | "all 139 signatures"; three-level audit | ✔ (the cross-language full-scale validation and the change census are unmentioned) | ✔ |
| M1 | §6.1 | 73 of 8,014; 32 / 28 in the products | 73; 32 / 28, all Canadian [wy*.area_normalized_false]; `flagged_for_qann_range` catches 22 / 25 of them | ✔ |
| M2 | §6.1 | "precipitation-dependent signatures (runoff ratios, elasticity, Q-P seasonality, storage) are structurally NA" | on the 32 / 28 gages: runoff ratios (5) and avg_storage 0 non-NA; elasticity_rolling / _annual 0; qp_* 0; BUT BFI_Eckhardt, BFI_LyneHollick, n_recession_events present for all 32 / 28, parameterized BFI 16 / 10, recession 9–10 / 2 (sparse, not gated) [wy*.unnorm_na.storage_detail] | ≈ → "average storage", not "storage" (the baseflow and recession signatures are computed) |
| M3 | §6.1 | unit-carrying list | ✔ | ✔ |
| N1 | §6.2 | "e.g. 0103500" | 7 digits; USGS ids are 8–15 digits | ✘ carried |
| N2 | §6.2 | "may be read as a character string" | must | ✘ carried |
| N3 | §6.2 | `Downstream_HB_ID` in the R4 attribute table; `basin_area`; `geom_area_km2` | Downstream_HB_ID is column 6 of the attribute table (and in the R3 metadata) [r4.attr_has_downstream]; basin_area in R1 / R2; geom_area_km2 in R4 | ✔ |
| N4 | §6.2 / §6.3 | "(Sect. 5.1.3)"; "compared against the WY … and WY … products" | 6.3; "between" | ✘ carried |
| O1 | §6.3 | record-dependent list | matches the readme for the annual signatures; it omits the per-gage scalars that are also computed from the window's record: `elasticity_static`, `recession_alpha_point_cloud_linear_reservoir`, the six `log_a_seasonality_*` scalars and the five `drought_threshold_fixed_*` (reviewer B) | ≈ |
| P1 | Availability | HydroShare URL; mirror URL; `docs/SIGNATURES.md` path; "licenceCC-BY 4.0 license" | mirror snapshot pending the next publish; typo | ✔ / ✘ |
| P2 | Acknowledgements | "Claude Code 0.145.0" | the project ran on releases from 1.x to 2.x (today's CLI 2.1.314) | ⚠ decision 7a: suggest "Claude Code (Anthropic), releases 1.x–2.x, 2025–2026" |
| Q1 | References | 22 (Lehner 2022, lakes) for HydroBASINS / BasinAtlas | 34 Lehner & Grill 2013 (HydroBASINS); 25 Linke 2019 (HydroATLAS) | ✘ |
| Q2 | References | 27 Myneni = MOD15A2H DOI | MCD15A3H v061, DOI 10.5067/MODIS/MCD15A3H.061 | ✘ carried |
| Q3 | References | Table 1 "Daymet21,26" | 26 | ✘ |
| Q4 | References | 23 (GAGES-II) for the ephemeral-stream rule | authors to confirm | ⚠ |
| Q5 | References | 29 Anderson 1976 for NLCD; 29 / 30 / 36 undated | fine as a secondary cite; formatting | ≈ |

## 7. Guidelines doc (synced 2026-10-30; reorder only) and the deposited dictionary

| ID | Where | Claim | Measured | Status |
|---|---|---|---|---|
| R1 | Part 1.1 | 16,994 / 9,154 / 7,840; "watersheds smaller than 85,000 km² and at least 20 total years … 8,014 (6,160 / 1,854) … > 111 million"; "Albers et al., 2026" | counts ✔ (C1, C2, C4); the 85,000 km² clause is the C3 conflation; tidyhydat is Albers 2017 | ≈ / ✘ |
| R2 | Part 1.2–1.3 | 3-day gap; 30 missing days; negative / constant-value screens; seasonal 80 %; trend 80 % / 80 % / 60 %; 20-value floor; recession + elasticity exempt; gage inclusion 20 years + 60 % | every number = config (`na_handling`, `stats_floor`); "60 % of the possible water years in the analysis window" has the same denominator ambiguity as H4 | ✔ (+ H4 wording) |
| R3 | Part 2 | 16 statistics, names; "pettitt_pre_mean: … before the changepoint year" | names = the dictionary suffixes ✔; the pre-segment INCLUDES the changepoint year (`changepoint.jl` l. 305–311, 333–334: year ≤ cp_year vs > cp_year; the reported year is the last pre-segment year) — Codex | ✔ / ≈ |
| R4 | Part 3 | every parameter: Q percentiles (15), D percentiles (13), FDC (10 days, Weibull, 1e-10, 3 points), pulses 90 / 10 %, reversal 2 %, elasticity 11-year window / 15 years / 10 mm / 0.1 mm, Q-P 30-day slope, runoff ratios 10 / 1 mm and > 2, snow 10 mm / 60 days / 50 % / 10 mm, climate 10 / 25 / 50 mm, 5 days, 50 %, 0 °C, Oct–Jun, GDD base 0, drought 7-day centered / Weibull / p2–p30, BFI 0.8 / 0.98 / 0.925 / 2 passes, recession 5 days / 25 events; §3.9 "recession and parameterized-BFI families are exempt from the trend-completeness requirements and the 20-value floor" | all parameters = `config/signatures_config.json` or the code (`fdc.jl` l. 29–51, `elasticity.jl` l. 130, `runoff_ratios.jl` l. 143, `baseflow.jl` l. 283); BUT the parameterized BFI signatures carry both gates (`signatures.jl` l. 124–128; `BFI_Eckhardt_param` 6,507 means vs 6,111 slopes in #1, 6,093 vs 5,762 in #2); the alpha interval is strictly (0, 1) (carried); the whole-record alpha scalar needs > 10 day-pairs and is computed independently of the 25-event gate (`recession.jl` l. 508–525) — Codex | ✔ parameters; ✘ the BFI-param exemption |
| R5 | §3.6 Snow | "for water years 1980–2023" | Daymet begins 1980-01-01, so WY 1980 has no climate-valid year: `swe_max` and `prcp_annual` start in WY 1981 (250,296 rows vs 255,145 for Qann in #2; Codex); the products run to WY 2025 | ✘ → "1981–2025" |
| R6 | §3.7 Climate | 39 metrics (20 prcp + 19 temp) | 39 [categories_by_cat]; the open naming / convention items of the 2026-10-12 relay stand | ✔ (carried items) |
| R7 | Part 4 | "USGS values are kept only for qualifiers A, A e, P, P e (others set to NA)" | the HIGH defect of 2026-10-21; rule now A/P status + finite value, remark codes carried (`R/usgs_qualifiers.R`; ledger D5) | ✘ DUE |
| R8 | Part 5 | 12 flags; high-NA denominator "16 statistics of every signature plus the 21 scalars" | 12 names = the dictionary's [wy*.flags]; 2,224 + 21 = 2,245 | ✔ |
| R9 | Part 6 | "Insert the table …", "Here is the original …" | working notes still on the page | ✘ carried |
| R10 | numbering | duplicate "1.1"; the new 3.3–3.9 order | repo notes citing "§3.9 Climate" now mean §3.7 — chore after this review | ⚠ repo-side |
| R11 | Part 4 | "Requires a pre-processed Daymet parquet file (data_out/daymet_1980_2023.parquet)" | the deposited input is the 1980–2025 table | ✘ → 1980–2025 |
| R12 | Part 3.1 / relay 2026-10-12 | "Negative_annual" → shipped name | the shipped column base is `negative_ann` (lower case) [wy*.negative_constant_cols; categories] — the 2026-10-12 relay item 1 wrote "Negative_ann" | ✘ correct the relay |
| S1 | dictionary (decision 1b) | 160 definitions vs SIGNATURES.md + guidelines | blind reviewer C (Sonnet, 2026-10-30, `docs/plans/2026-10-30-manuscript-reviews/review_C_definitions.md`): 0 HIGH; 7 MEDIUM rows = missing caveats (fall-season origin; "NA by design for un-normalized gages" on 14 Q-to-PPT rows; 25-event gate on 12 recession rows; parameterized-BFI NA rule; "NA in snow-free years" on 7 snow rows; avg_storage's interpolation step; drought "strictly below" / "0 not NA"); 8 LOW. Codex adds a MAJOR beyond the `_mean` rows: every `_pettitt_pct_change` definition says "percentage of the pre-changepoint mean" while the code uses the ABSOLUTE pre-mean (`changepoint.jl` l. 345–348; NA when |pre| ≤ 1e-10) — sign-relevant for the 6,578 negative tmin pre-means in #1 (01011000: pre −34.13, delta +4.03, reported +11.8 %); and the pre / post segment rule of R3 | ≈ dictionary-side list for the user (changing it means re-staging R1 / R2) |
| X1 | §1, §2.1.4 | literature and external-product facts (ERA5-Land ~11 km; CAMELS-SPAT n < 1,500; HydroBASINS 15 arc-seconds; HydroATLAS six categories; MODIS 500 m; LAI algorithm; LAI 0–10) | not data claims; not assessed | — |

## 8. Blind reviews (2026-10-30, `docs/plans/2026-10-30-manuscript-reviews/`)
- **A — numbers (Sonnet, blind to this ledger):** ~75 claims, ~40 match, 29 mismatch, 6 not measurable. Every mismatch is a ✘ row above; no contradiction with the census. Extra: 492,956 streamflow rows precede 1980-01-01 (WY 1980 starts 1979-10-01 — folded into D1); A counts 204 attributes (7 id columns) where the deposited dictionary counts 198 (+13 keys / diagnostics) — the ledger keeps the dictionary's own split.
- **B — methods (Sonnet, blind):** 7 WRONG (= D9, D10, D4 / E3, E4, F1 / F2, F4, H5, I2 and the cross-references), 10 IMPRECISE; independently confirms the H4 denominator; three new rows added from it: I5, I6, O1.
- **C — dictionary definitions (Sonnet):** 0 HIGH, 7 MEDIUM omissions, 8 LOW (row S1).
- **Codex — adversarial read of this ledger (gpt-6-astra, read-only):** NO-GO as first written; 11 MAJOR + 8 MINOR, every one checked and applied above (H4 raw-record denominator and 0.600 minimum; D4 7,960 + 4; R4 parameterized-BFI exemption wrong in the guidelines; S1 `_pettitt_pct_change` absolute denominator; H6 and G1 downgraded to ⚠; R5 annual values start WY 1981; H5 eight named metrics; G7 "all HISSS watersheds"; J18 `site_id`; D7 / D9 internal-only interpolation and boundary gaps; D5 sentinel and ≈; C4 records vs observations; R3 segment rule). Census script tightened on its request and re-run (205 measurements): finite-value filter, (gage, signature, year) keys unique in both annual parquets, no duplicate dates per gage, all six Daymet variables free of null / NaN / inf on a complete 46 × 365 grid, `negative_ann` selector, the suffix measure dropped.

## 9. Not verified / left open
- K1 ecoregions, K3 Figure 3 coverage (out of scope by decision).
- Q4 (source of the > 30-day rule) and P2 wording: the authors'.
- J1: the live-resource file lists are the sibling session's; this ledger measured the staged copies.
- The dictionary's 15 non-`_mean` rows per signature were not read by reviewer C.

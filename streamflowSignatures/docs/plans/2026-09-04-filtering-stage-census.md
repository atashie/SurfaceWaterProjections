# Filtering-stage watershed census (2026-09-04)

Counts of watersheds (gages) at every filtering stage of the HISSS workflow schematic
(`dataset_workflow_schematic.{md,png}`), measured directly from the delivered inputs and
products on 2026-09-04, and cross-checked against the manuscript text published the same
day. Every number below was recomputed — none is copied from earlier docs.

**Method.** (1) Candidate → usable: `combined_watershed_metadata_09feb2026.csv` (16,994
rows) crossed by `processing_status × gage_type`, and the streamflow parquet's own gage
set. (2) Year rejection + gage inclusion: `docs/benchmarks/qualification_census.jl` runs
the canonical `preprocess_daily_data` on **all 8,014 gages** and re-implements the
runner's window / 60 %-fraction / 20-year gates for both standard windows; its
include/exclude decisions reproduce the delivered products **exactly** (6,678 and 6,250
gages, 0 mismatches either way). Summarized with
`docs/benchmarks/summarize_qualification_census.py`. (3) Coverage branches: exact and
canonical (zero-stripped) id joins between the staged HydroShare resources
(`~/Downloads/Signatures/resource{3,4,5}_*`), the Daymet parquet's distinct `site_id`s,
and each product's gage set. (4) Trend gate: non-NA counts of `_mean` vs `_senn_slp` vs
`_pettitt_cp_year` per signature base in each product CSV.

## 1. Streamflow path (the schematic's vertical chain)

| Stage (schematic node) | Count | US / Canada | Source / note |
|---|---|---|---|
| Candidate gages retrieved (USGS NWIS + HYDAT) | **16,994** | 9,154 / 7,840 | metadata CSV rows |
| ↳ returned no valid daily data (`no_data`) | 8,787 | 2,812 / 5,975 | "No valid streamflow data retrieved" |
| ↳ fewer than 20 years of valid data (`insufficient_years`) | 80 | 69 / 11 | 0–19 years (8 with 0, 6 with 19) |
| ↳ retrieval never completed (`processing`) | 113 | 113 / 0 | no year counts, no error message; 110 CONUS + 3 AK |
| **Usable daily records (Standardized daily streamflow)** | **8,014** | 6,160 / 1,854 | = streamflow parquet gage set; 111,624,189 rows; dates 1979-10-01 → 2025-10-01 |
| ↳ of which NOT area-normalized (no published drainage area) | 73 | 0 / 73 | `area_normalized = FALSE`; 1,601 across all 16,994 candidates |
| Water years with data, WY 1980–2025, all 8,014 gages | 317,182 gage-years | | + 5,381 one-day WY 2026 stubs (2025-10-01 only), all rejected |
| ↳ **rejected by the preprocessor** | **33,732** (10.6 %) | | 31,091 >30 raw NA days · 1,803 NA run >3 days · 838 residual boundary NA (see §5.3) |
| ↳ kept (valid water years) | 283,450 | | |
| Gages with ≥20 valid years over the full record (no window) | 7,313 | | = the April 2026 full-record canonical gage count |
| **Include gage? — WY 1993–2025 window** | **6,678** | 5,495 / 1,183 | product #1 (`processedOuts_drought_28jul2026`) |
| ↳ excluded | 1,336 | 665 / 671 | 924 fail both gates · 412 fail only the 20-year floor · **0** fail only the 60 % rule · 374 have zero valid years in-window |
| ↳ valid years per included gage | 20–33, median 32 | | 158 gages at exactly 20; 200,834 gage-years in the annual parquet |
| ↳ rejected years within the included gages | 6,668 of 207,502 (3.2 %) | | 5,268 NA-count · 1,059 gap · 341 boundary; 3,336 gages lose no year at all |
| **Include gage? — WY 1980–2025 window** | **6,250** | 5,043 / 1,207 | product #2 (`processedOuts_1980_2025_11aug2026`) |
| ↳ excluded | 1,764 | 1,117 / 647 | **1,063 fail only the 60 % rule (they have ≥20 valid years)** · 674 fail both · 27 fail only the 20-year floor · 341 zero valid years |
| ↳ valid years per included gage | 20–46, median 44 | | 38 gages at exactly 20; 254,195 gage-years |
| ↳ rejected years within the included gages | 7,689 of 261,884 (2.9 %) | | 6,095 NA-count · 981 gap · 613 boundary; 2,866 gages lose no year |
| **Area-normalized gage?** — un-normalized in product | 32 (#1) · 28 (#2) | all Canada | Q-to-PPT signatures skipped; 37 in the full-record canonical |
| Trend gate (per signature) — dense signatures, e.g. `Qann` | 6,183 of 6,678 (#1) · 5,851 of 6,250 (#2) carry trend statistics | | 495 / 399 gages fail the 60 %-series / 80 %-decade / 20-value gate; Pettitt fields present for 6,590 / 6,250 (WY 1980–2024 window ⇒ 88 gages in #1 lack them) |
| Trend gate — sparsest signatures | `dur_low_pulses_all` 2,614 of 4,727 (#1) · 2,449 of 5,293 (#2); snow timing/melt 2,763–3,373 (#1) · 2,741–3,544 (#2) | | recession + elasticity are exempt (trend count == mean count); modal count 6,183 / 5,851 on 45 of 100 bases |

Signature counts per product: 121 signatures (100 annual series + 21 per-gage scalars)
across 14 families → 1,653 columns (20 metadata + 100 × 16 + 21 + 12 QA flags); annual
parquets 18,898,406 / 24,366,487 rows, 100 signatures each. QA flags fired (#1 / #2):
`high_na` 1,224 / 1,243, `qann_range` 159 / 151, `elasticity_range` 91 / 81,
`runoff_ratio_range` 88 / 86, all others 0.

## 2. Coverage branches (the schematic's right-hand column)

| Layer | Universe | In product #1 (6,678) | In product #2 (6,250) | Note |
|---|---|---|---|---|
| Watershed boundaries (Resource 4) | **7,964** = 8,018 − 54 basins > 100,000 km² | 6,634 (44 of the 54 big basins qualify) | 6,205 (45 qualify) | sources gagesii 6,164 / wsc_eccc 1,771 / hydrobasins 29; 4 geometry-only gages (01591000/400/610/700, "processing" status, never compiled); 56 `low_confidence` |
| HydroATLAS attributes (Resource 4) | **8,014** × 211 columns | 6,678 | 6,250 | 211 = 198 aggregated attributes + 13 key/diagnostic columns |
| Daymet basin series (Resource 3) | **6,087** sites; **5,965** of the 8,014 | **5,517** (4,489 US / 1,028 CAN) | **5,638** (4,549 US / 1,089 CAN) | 122 of the 6,087 belong to never-compiled candidates; CY 1980-01-01 → 2023-12-31, 97,757,220 rows |
| ↳ snow metrics actually computed (`swe_max_mean` non-NA) | | 5,485 (4,457 / 1,028) | 5,635 (4,546 / 1,089) | 32 / 3 Daymet gages have no SWE-valid year in window |
| ↳ runoff ratios computed (`annual_runoff_ratio_mean` non-NA) | | 5,471 | 5,619 | Daymet ∧ area-normalized ∧ PPT-valid years |
| MODIS LAI / LULC (Resource 5) | 7,964 (2,150,280 monthly rows; 191,136 annual rows) | **6,634** | **6,205** | = the boundary layer; exact-string join gives only 6,599 / 6,196 (see §5.1) |
| Annual NLCD (Resource 5) | **6,119** CONUS = 6,164 US − 45 Alaska | **5,454** (41 US gages uncovered) | **5,007** (36 uncovered) | exact-string join gives only 5,419 / 4,998 (§5.1) |

## 3. Manuscript cross-check (published text as of 2026-09-04)

| Manuscript claim | Measured | Verdict |
|---|---|---|
| 16,994 candidates (9,154 USGS; 7,840 HYDAT) | 16,994 (9,154 / 7,840) | ✔ |
| 8,980 excluded "fewer than 20 total years"; 8,014 usable (6,160 US; 1,854 CAN) | 8,980 = 8,787 no data + 80 <20 yr + **113 never-completed**; 8,014 (6,160 / 1,854) | ✔ numbers; wording glosses the 113 (§5.2) |
| "more than 111 million observations" / 111.6 M / 111,624,189 | 111,624,189 | ✔ |
| retrieval "01 October 1979 to 30 September 2025"; Resource 3 "1980 through 2025" | parquet ends **2025-10-01** (5,381 gages carry that single WY 2026 day; all rejected) | ✔ in effect; note §5.4 |
| 73 gages not area-normalized; 32 / 28 in the products (§5.1.1) | 73; 32 / 28 | ✔ |
| WY rejection: >30 missing values or a gap >3 days | + a third rule: residual boundary NA (838 of 33,732) | ≈ (§5.3) |
| ≥20 qualifying WYs and ≥60 % of window → 6,678 / 6,250 | 6,678 / 6,250, reproduced exactly | ✔ |
| "1993-2025 … date range with the most gages with the same record length" (new sentence) | 6,678 > 6,250; 2,909 of the 6,678 carry the complete 33-year record (4,660 have ≥30 years) vs 2,359 of the 6,250 with the complete 46-year record | ✔ (defensible) |
| 54 basins > 100,000 km² excluded; 7,964 boundaries | 54 (49 with published area > 100k, 5 with no published area but delineated > 100k); 7,964 | ✔ |
| GAGES-II covers 100 % of US gages; ECCC 98.3 % of Canadian | 6,164 / 6,164; 1,823 / 1,854 = 98.3 % | ✔ |
| Daymet: 6,087 candidate basins; 5,965 of 8,014; 5,517 / 5,638 in products | 6,087; 5,965; 5,517 / 5,638 | ✔ |
| HydroATLAS "approximately 210 attributes" (§2.1.4) / "211 watershed-scale attributes" (§3) | 211 columns = 198 attributes + 13 key/diagnostic | ≈ (§5.5) |
| MODIS 102 class columns; LAI 2,150,280 rows = 7,964 × 270; LULC 191,136 | 102; 2,150,280 (270 months, 2002–2024); 191,136 (2001–2024) | ✔ |
| Annual NLCD 6,119 CONUS watersheds; 250,879 rows; 16 classes; 1985–2025 | 6,119; 250,879; 16; 1985–2025 | ✔ |
| 121 signatures / 100 annual / 21 scalars / 14 categories / 16 statistics / 1,653 columns / 12 flags | all match | ✔ |
| annual-values rows 18,898,406 / 24,366,487 | match | ✔ |
| "All tables join on gage_id, the zero-padded agency station identifier" | FALSE for 9 gages in Resource 4 and 44 in Resource 5 (§5.1) | ✘ until the staged files are re-padded |

## 4. Schematic

Nothing in `dataset_workflow_schematic.{md,png}` is contradicted. If counts are wanted on
the figure, the defensible per-node numbers are: 16,994 → 8,014 (usable) → 7,964
boundaries; year rejection 33,732 of 317,182 gage-years (10.6 %); include gage 6,678 /
6,250; un-normalized 32 / 28; trend gate ~92.6 % / 93.6 % of gages on dense signatures.

## 5. Findings

### 5.1 DEPOSIT DEFECT — `gage_id` is not zero-padded for long USGS ids in Resources 4 and 5
Nine gages in `hisss_watershed_boundaries.{gpkg,parquet}` + `_qa.csv` and **44** in each of
`hisss_modis_lai_monthly.parquet`, `hisss_modis_lulc_annual.parquet`, and
`hisss_nlcd_annual.parquet` carry the zero-stripped form (`212433550`, `21650905`,
`72632962`, `11055566`, …) where the streamflow, signature and HydroATLAS tables carry the
agency form (`0212433550`, `021650905`, `072632962`, `011055566`, …). Cause: the EO
pipelines' `zfill(8)` rule cannot restore leading zeros on ids that are 8–10 digits long
after stripping (9- and 10-digit USGS ids with a leading zero); the geometry rebuild
fixed most of them but missed 9. Consequences: (a) the manuscript's and the Resource 4/5
READMEs' "join directly on gage_id, no re-padding needed" claim is false for those
gages; (b) the R5 README's coverage counts (6,599 / 6,196 MODIS; 5,419 / 4,998 NLCD)
were measured with exact-string joins and undercount — canonical values are **6,634 /
6,205 and 5,454 / 5,007**. `canon_id` joins correctly everywhere. HydroATLAS (8,014) and
Resource 3 are clean. **DECISION (user, 2026-09-04): the staged files stay as delivered; the join rule is
documented instead** (root README → "Reading and Joining the Data"; the Resource 3/4/5
READMEs and dictionaries; the pad-to-8 recipe that Resource 3's README carried also fails
for 66 of the 8,014 compiled gages and was replaced by strip-on-both-sides). Superseded
recommendation kept for the record: rewrite `gage_id` in the five
affected staged files from `canon_id` using the streamflow parquet's ids as the
authority (`hisss_nlcd_excluded_gages.csv` is already correct), then re-run the
coverage counts and update the R5 README + the CHANGELOG 2026-08-25 R5 entry.** The
affected-id list is reproducible from the census script's inputs (44 MODIS/NLCD ids ⊃ the
9 geometry ids).

### 5.2 The "8,980 excluded" bucket contains 113 gages that were never processed
`processing_status == "processing"` (110 CONUS + 3 Alaska USGS gages, no year counts, no
error message) — retrieval never completed rather than "fewer than 20 years". Four of
them (01591000, 01591400, 01591610, 01591700) do have official GAGES-II polygons and sit
in the geometry/EO layers (the "8,018 universe"). Manuscript §2.1.2 should either say
"returned fewer than 20 years of valid daily observations or could not be retrieved
(8,980)" or footnote the 113.

### 5.3 A third year-rejection rule exists and is not in the manuscript
Interpolation is per water year and never extrapolates, so a ≤3-day NA run touching Oct 1
or Sep 30 survives interpolation and rejects the year (`residual_na`): **838** of the
33,732 rejected gage-years (2.5 %). Documented in DEVELOPMENT.md (step 5) and correct in
code; §2.1.2 / §2.2.2 list only the >30-day and >3-day rules. One clause fixes it: "…or
any missing day at the start or end of the water year that could not be interpolated".

### 5.4 Record end date
The streamflow parquet's last date is 2025-10-01 (one WY 2026 day for 5,381 gages, all
rejected). The manuscript's "to 30 September 2025" describes the analysed record
correctly; Resource 3's "1980 through 2025" is fine. No change needed unless the deposit
parquet is trimmed to WY 2025 (then 111,618,808 rows — not recommended, keep the
delivered bytes).

### 5.5 HydroATLAS attribute count
The staged table has 211 columns: 198 HydroATLAS-derived attributes (92 area-weighted
means, 91 `_u`/`_p` passthroughs, 11 categorical majorities + 2 fractions, elevation
min/max) plus 13 key/diagnostic columns (`gage_id`, lat/lon, `basin_area`, `gage_type`,
`Downstream_HB_ID`, `watershed_area_rel_diff`, …). §2.1.4's "approximately 210" is fine;
§3's "211 watershed-scale attributes" should read "211 columns (198 attributes plus keys
and diagnostics)".

### 5.6 Still open from earlier passes (verified present in the 2026-09-04 text)
Reference list lacks entries for McMillan 2021, Hatchett 2021, Petersky & Harpold 2018,
Adelsperger (in review), Laaha 2017, Peters & Aulenbach 2011, HYSETS (Arsenault 2020),
Caravan (Kratzert 2023), CAMELS-SPAT (Knoben 2025), DeCicco, Albers, gdptools, Annual
NLCD; "Linke et al., 2013" is still cited (nonexistent; 2019 is the HydroATLAS paper);
Myneni / Friedl "[tbd]"; "Claude Code 0.145.0" ×2 (that is the codex-cli version);
"fitted checksums, configuration, and software versions)" garble in §3; typos "wtih",
"diagnostics.We" (missing space), "s(e.g.,", "MOIDS"; §3 CC-BY licence and the
collection "DOI" label to confirm; Table [X] / schematic / input-data table placeholders.

The stage table is also provided as a CSV for sharing: `docs/filtering_stage_census_2026-09-04.csv`
(kept outside `docs/plans/` so it is included in the public CZ-Sync/HISSS mirror).

## Reproduce
```bash
STREAMFLOW_DATA_PATH=/Volumes/Untitled/processedOuts_feb2026/combined_streamflow_data_09feb2026.parquet \
  julia --project=julia docs/benchmarks/qualification_census.jl /path/census.csv        # ~3 min
python docs/benchmarks/summarize_qualification_census.py /path/census.csv \
  /Volumes/Untitled/processedOuts_drought_28jul2026/streamflow_1993_2025_60pct_drought_28jul2026_signatures.csv \
  /Volumes/Untitled/processedOuts_1980_2025_11aug2026/streamflow_1980_2025_60pct_11aug2026_signatures.csv \
  /Volumes/Untitled/processedOuts_feb2026/combined_watershed_metadata_09feb2026.csv
```


## 7. Addendum 2026-09-15 — coverage intersections for the manuscript's §2.1 placeholder

Measured on the staged deposit (`~/Downloads/Signatures`), joining every table on the
zero-stripped `gage_id` (script: session scratch `intersections.py` / `daymet_why.py`;
Daymet ids = distinct `site_id` of `hisss_daymet_basin_daily.parquet`, 6,087; MODIS ids =
distinct `gage_id` of `hisss_modis_lulc_annual.parquet`, 7,964; NLCD 6,119).

| Set | WY 1993–2025 | WY 1980–2025 |
|---|---|---|
| signature gages | 6,678 | 6,250 |
| … with Daymet | 5,517 | 5,638 |
| … with MODIS LULC/LAI | 6,634 | 6,205 |
| … with Daymet AND MODIS | **5,517** | **5,638** |
| … with Daymet AND MODIS — in EITHER product / in BOTH (added 2026-09-18) | 5,665 / 5,490 | (union / intersection of the two gage sets: 7,157 / 5,771) |
| … with Daymet AND MODIS — all 8,014 usable gages (added 2026-09-18) | 5,965 | — |
| … with MODIS but no Daymet | 1,117 | 567 |
| … with NLCD | 5,454 | 5,007 |
| … with Daymet, MODIS and NLCD | 4,464 | 4,524 |
| … without Daymet | 1,161 | 612 |
| ↳ polygon area > 85,000 km² | 6 | 6 |
| ↳ agency `basin_area` > 85,000 km² | 46 | 47 |

**2026-09-18 check of the manuscript's "6,041" (the §2.1 XXXX fill-in)**: matches nothing
above. Ceiling is 5,969 = Daymet sites with a polygon (6,087 − 118 without one); 122 of the
6,087 Daymet sites are not usable-streamflow gages (109 `processing`, 13 `no_data`); every
site has all 16,060 rows (5 carry some NaN prcp/swe). Script: session scratch `check_6041.py`.

Every Daymet site is also a MODIS gage (8,014-level: MODIS 7,960, Daymet 5,965, both
5,965), so "signatures + climate + MODIS" is simply the Daymet count. **Why 2,049 of the
8,014 usable gages have no Daymet series is NOT basin size**: only 55 of them exceed
85,000 km² by agency area. The 6,087-basin Daymet request list predates the Feb 2026
compilation — all 5,707 gages of `golden-outputs/streamflow_signatures_full_10feb2026.csv`
have Daymet (Daymet has 380 more, 122 of them never-compiled candidates); Daymet gages are
the long-record ones (median 46 record years, 75 % starting 1979) and the no-Daymet gages
are the shorter/later records the 20-year rule admits (median 27 years, median start 1992,
end year < 2023 for 543 of them). The exact selection rule of the January 2026 request is
not recorded in this repo.

## 8. Addendum 2026-09-15 — measured counts of every staged deposit file (for the manuscript's resources table)

Review of the co-authors' resources sheet (`1Hx8hnr02c…`, tab 1694936652, column C "Data
Rows"). Counted directly on `~/Downloads/Signatures` with pyarrow/pandas; ids joined
zero-stripped.

| Resource / file | Gages (distinct) | Rows |
|---|---|---|
| R1 `hisss_signatures_wy1993_2025.csv` | 6,678 | 6,678 (1,653 columns) |
| R1 `hisss_signatures_annual_wy1993_2025.parquet` | 6,678 | 18,898,406 |
| R2 `hisss_signatures_wy1980_2025.csv` | 6,250 | 6,250 |
| R2 `hisss_signatures_annual_wy1980_2025.parquet` | 6,250 | 24,366,487 |
| R1/R2 `hisss_data_dictionary.csv` / `hisss_signature_categories.csv` | — | 1,653 / 100 |
| R3 `hisss_streamflow_daily.parquet` | 8,014 | 111,624,189 |
| R3 `hisss_gage_metadata.csv` | 16,994 candidates (8,014 success, 113 processing, 80 insufficient_years, 8,787 no_data) | 16,994 |
| R3 `hisss_daymet_basin_daily.parquet` | 6,087 sites (5,965 of the 8,014) | 97,757,220 |
| R3 `hisss_hydat_interference.csv` | 8,012 stations | 8,012 |
| R4 `hisss_watershed_boundaries.{gpkg,parquet}` + `_qa.csv` | 7,964 (7,958 ≤ 85,000 km² by polygon area; max 97,182 km²; the exclusion rule was 100,000 km², 54 dropped) | 7,964 |
| R4 `hisss_watershed_hydroatlas_attributes.{csv,parquet}` | 8,014 (every success gage; 54 have no polygon) | 8,014 × 211 columns |
| R5 `hisss_modis_lai_monthly.parquet` | 7,964 | 2,150,280 |
| R5 `hisss_modis_lulc_annual.parquet` | 7,964 | 191,136 |
| R5 `hisss_nlcd_annual.parquet` (+ 45 excluded Alaska gages) | 6,119 | 250,879 |

Sheet values that match nothing measured: **6,041** (R3), **8,044** ("total streamflow
timeseries" → 8,014), **6,014** ("sites with daymet and modis" → 5,969 in the boundary
layer, 5,965 of the 8,014 usable gages, 5,517 / 5,638 within the products). R4's **8,014**
is the attribute table, not the boundary layer (7,964).

## 9. Addendum 2026-09-28 — why 1,995 polygon gages have no Daymet (the rule is now recovered)

Co-author question: is the gap between the 7,964 polygons and the Daymet sites entirely due to
processing, with the Daymet sites under stricter data requirements? **Mostly yes.** The Daymet
request list was built from the ORIGINAL compilation, preserved in git as
`combined_streamflow_output/combined_watershed_metadata_og.csv` (added in `27c38d9`,
2026-02-12; 16,994 rows: 6,046 success, 134 processing, 81 insufficient_years, 10,733 no_data).
That run used **`min_num_years = 30`** (raw-data script and `config.R`; lowered to 20 in
`27c38d9`, which produced the final 8,014). Measured on the staged deposit, ids zero-stripped:

- Daymet list = 5,956 og-success + 131 og-processing (retrieval unfinished) = 6,087.
- **Polygon gages without Daymet = 1,995:**
  - **1,952 failed the 30-year compilation**: 1,910 og `no_data` ("No valid streamflow data
    retrieved"), all with < 365 × 30 = 10,950 daily records (max 10,943), so they failed the
    retrieval pre-check `nrow > 365 × min_num_years`; 42 og `insufficient_years` (20–29
    qualifying calendar years). All passed the final 20-year compilation.
  - **6 > 85,000 km²** (Canadian; polygons 87,150–97,182 km²) — the Daymet size rule.
  - **22 USGS gages with 9–10-digit site numbers** — none of the 66 such usable gages has
    Daymet; the same leading-zero id-format loss as §5.1 (the request list almost certainly
    mis-keyed them).
  - **15 Canadian HydroBASINS-fallback polygons** — none of the 29 fallback polygons has
    Daymet; the Daymet run used an earlier boundary set without them.
- Usable gages WITH Daymet (5,965): 5,952 passed the 30-year compilation, 13 were still in
  retrieval in the og run. The other direction: 118 Daymet sites have no polygon (105
  `processing`, 13 `no_data` in the final metadata).
- "Stricter" means record LENGTH at compilation only (≥ 30 vs ≥ 20 calendar years with ≥ 30
  days of Q > 0.0001); per-year QA, the 60 % window rule and the trend gates are identical
  for all gages in the products. Separation in the final metadata: 5,960 of 5,965 Daymet
  gages have ≥ 30 qualifying years vs 564 of 2,049 no-Daymet gages.
- **Correction to §7**: "all 5,707 golden gages have Daymet" is circular, not evidence of
  vintage — `process_signatures_from_parquet()` iterates over the Daymet `site_id`s when a
  Daymet path is given, so the R golden run could only contain Daymet gages.
- The unexplained "6,041" (§7) is within 5 of the og run's 6,046 successes; probably a count
  from that run, not verified.


## 10. Addendum 2026-10-21 — what the compilation's flow rule removed (user question)

Question: how many candidate gages did "at least 20 years of valid streamflow data with
>30 days over .0001 mm/day" remove, and how many of those because of the flow threshold itself?

**Rule as coded** (`R/helperFunctions.R` → `process_gages_rawToRaw` + `generate_streamflow_dt`;
`config.R` `MIN_Q_VALUE_AND_DAYS = c(0.0001, 30)`, `MIN_NUM_YEARS = 20`): (i) retrieval
pre-check — more than 365 × 20 = 7,300 daily rows inside 1979-10-01 … 2025-10-01, else
`no_data` (no flow value is examined); (ii) USGS values whose qualifier string is not exactly
`A`, `A e`, `P` or `P e` are set to NA (HYDAT symbols are never masked); (iii) a CALENDAR year
counts when at least 30 non-NA values exceed 0.0001 (mm/day after area normalization; raw
units for gages without a drainage area); (iv) fewer than 20 counted years →
`insufficient_years`. Once a gage qualifies every day in the window is kept, and no flow
threshold is applied downstream (`preprocess_daily_data` rejects water years on NA count and
gap length only). The `config.R` comment "values below are set to NA" does not describe this
path — values are only counted.

**Method.** The 80 `insufficient_years` gages (69 USGS, 11 HYDAT) were re-fetched 2026-10-21
(NWIS legacy daily-values service, the one `readNWISdv` used in February; ECCC GeoMet
`hydrometric-daily-mean`, served from HYDAT) and the gate replayed under four variants of rule
(iii) with `docs/benchmarks/min_flow_rule_replay.py`; the 8,014 compiled gages were scanned
with `docs/benchmarks/min_flow_rule_footprint.py` on the local copy of
`combined_streamflow_data_09feb2026.parquet` (899,935,349 bytes, `PAR1` footer, 111,624,189
rows). Replay check: 76 of the 80 recorded year counts reproduce exactly, 3 differ by one year
(records extended since February), and 02441390 now has 37 qualifying years (16 in February —
its record was revised; it would pass today). Per-gage table:
`docs/filtering_min_flow_rule_2026-10-21.csv`.

| Excluded candidates (8,980) | Gages | US / Canada | What removed them |
|---|---|---|---|
| record-length pre-check (`no_data`) | 8,787 | 2,812 / 5,975 | ≤ 7,300 daily rows in the window — the flow threshold never ran |
| retrieval never completed (`processing`) | 113 | 113 / 0 | — |
| per-year rule (`insufficient_years`) | 80 | 69 / 11 | split below |
| ↳ **zero flow** | **67** | 57 / 10 | ≥ 20 years with 30 valid days but < 20 years with 30 days of ANY positive flow — dry channels, diversions, canals (Pecos River 08398500–08405150, Trent Canal lock 02EC016, 09513860, 07232250/900 …); here the threshold acts as "> 0" |
| ↳ **the 0.0001 mm/day level itself** | **5** | 5 / 0 | ≥ 20 years with 30 positive days, < 20 above the threshold: 06863000 (21 → 19 years), 06871800 (23 → 7), 07155590 (22 → 8), 07305500 (20 → 19), 09386300 (20 → 18); 0.0001 mm/day is 0.01–0.62 cfs at their drainage areas |
| ↳ qualifier-code mask | 4 | 4 / 0 | 02201230, 07069000, 07074420, 07188653: 23–27 qualifying years without the mask, 15–19 with it — their revised periods (`A R`) were NA |
| ↳ record length | 3 | 2 / 1 | 02479310, 04084445, 05NB039: fewer than 20 years with 30 valid days under any value rule |
| ↳ now passes | 1 | 1 / 0 | 02441390 (revised since February) |

**Answer.** The flow threshold proper removed **72** gages (67 zero-flow + 5 at the level); the
flow-and-length rule as a whole removed 8,867 (8,980 less the 113 never retrieved), 8,787 of
them on record length alone.

**Footprint among the 8,014 compiled gages** (calendar years 1979–2025, 322,232 gage-years):
2,269 gage-years (0.7 %) in 489 gages have ≥ 30 valid days but < 30 days above 0.0001 — 2,044
of them (455 gages) are zero-flow years, 225 (100 gages) sit in the band 0 < Q ≤ 0.0001 mm/day.
These years did not count towards the 20, but every one of them is in the compiled data and
in the products (admission is unaffected by construction: 7,525 gages have no such year; one
gage has 23).

### 10.1 FINDING (HIGH, every product) — the qualifier-code mask drops approved, revised values
`generate_streamflow_dt` keeps a USGS value only when its qualifier string is exactly `A`,
`A e`, `P` or `P e`. NWIS's `R` means "Records for these data have been revised" — an
APPROVED value — yet every `A R` day is NA in the compiled record: **181,461 days on 755
gages** (1,571 gage-water-years; 934 with > 30 masked days; 277 whole years). Replaying the
preprocessor's two rules (> 30 NA days, NA run > 3 days) on those gage-years with and without
the mask: 1,215 are rejected with it, 1 without — **1,214 water years in 613 gages are lost
solely to the mask** (worst: 07362000 16 years, 03250322 15, 07056000 12, 02313700 12,
01585200 11). 553 of the 613 gages are in product #1 v2 and 447 in #2 v2; 71 / 193 of the 755
`A R` gages are outside the products (how many because of the mask was not separated).
Revisions cluster in the most recent years (2013–2025 for most of the worst gages), so record
ends, trend last-decade completeness and the drought/pulse `*_all` thresholds are hit hardest.
Smaller masks: `A <` 11,444 days / 139 gages and `A >` 218 (values at a reporting limit),
`A [4]` 419 / 101, and 15,257 days of provisional values with remark codes (`P Dis` 3,873,
`P Rat` 3,763, `P Ice` 3,693, `P ***` 1,372, `P Bkw` 889, `P Eqp` 416, `P Fld`, `P [4]`, `P <`,
`P ZFL` 193 — zero flow is a valid value —, `P Ssn`, `P Mnt`, `P >`), most of which are
defensible to drop. The manuscript's "we retained original agency quality flags, including
ice-affected and estimated values" is true of the `flag` column but not of the values behind
`A R`. **Fix (NOT applied — user decision): accept `R` (and decide `<`, `>`, `ZFL`) in
`keep_codes`, re-ingest the affected USGS gages, rerun.** CORRECTION 2026-10-23 (Codex review of
the rerun plan): the 1,214 / 613 budget above counted NAs among existing rows only; on the full
October–September calendar (missing days count as NA, as the preprocessor does) the loss is
**1,197 water years in 601 gages for `A R` alone and 1,370 in 676 gages with `A <`, `A >`,
`A [4]` included** (1,177 in product #1 v2 / 950 in #2 v2; 899 / 757 product gages hold a masked
approved day inside their window); the union of the four approved codes is 193,542 days on 987
gages; 11129800 is a fifth compile-floor casualty (20 qualifying calendar years without the mask).
Plan and decisions (2026-10-26, gage sets frozen):
`docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md`. Tables: `masked_days_by_gage_wy_flag.csv`
and `AR_gage_wy_budget.csv`, written by the footprint tool. Logged: CHANGELOG → Known Issues,
`docs/STATUS.md`.

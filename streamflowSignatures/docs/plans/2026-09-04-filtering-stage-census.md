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

# HISSS manuscript — reconciliation log

Dated entries, newest first, written by `/sync-docs` after each sync of the manuscript
Google Doc (`docs/MANUSCRIPT_DRAFT.md` is the read-only snapshot). Each entry groups the
findings by direction of fix: manuscript-side (relay to the co-authors) vs code/docs-side
(implementation TODO). Moved here VERBATIM from `CHANGELOG.md → [Unreleased] →
Manuscript Reconciliation Log` on 2026-09-29.

---

Session-start reconciliation of the HISSS manuscript draft (Scientific Data,
submission target Nov 9 2026) against code + repo docs. Snapshot:
`docs/MANUSCRIPT_DRAFT.md`; workflow: `/sync-docs` §3.

**Earlier passes are in changelog-old.md → Manuscript Reconciliation Log** — the
2026-07-21 baseline with the eight queued manuscript edits, the 2026-08-24 second pass and
§2.1.4 / §3 drafting support, the 2026-08-28 third pass with its numbered open-item list
1–12, and the 2026-09-01 fourth pass with the paste-ready correction blocks 1a–7
(`docs/plans/2026-09-01-manuscript-correction-blocks.md`) and the §2 workflow schematic
(`docs/plans/dataset_workflow_schematic.md`). The queued-edit and item numbers cited below
refer to those lists. One thread carried forward from 2026-09-01: the committed schematic
fails GitHub's rich rendering because of a GitHub-side mermaid bundle crash (verified
against GitHub's own mermaid README, not this file) — re-check after GitHub ships a fixed
bundle; until then the 3x PNG beside the doc is the review copy.

**2026-10-21 — co-author edits synced (106 vs 109 paragraphs; 14 regions); reconciliation
pass.** Snapshot overwritten. Thirteen of the fourteen regions are the switch from author–year
to numbered citations (1–36), the removal of the three title lines, the §5.1.x → §5.1–5.3
renumbering and a Figure 2 inset sentence. Checked against code, config and the census
(`docs/plans/2026-09-04-filtering-stage-census.md`, new §10). Findings by direction of fix:
- **Manuscript wrong or imprecise, code/docs right (relay):**
  1. **§2 funnel sentence** now reads "…at least 20 years of valid streamflow data with >30 days
     over .0001 mm/day and were located in watersheds with polygons smaller than 85,000 km
     (≤ largest HUC8 size), yielding 7,964 gages" — two stages fused: the streamflow rule yields
     8,014 (§2.1.1 and §3 Resource 3 still say so), the boundary layer 7,964; the exclusion
     rule applied was 100,000 km² (54 basins dropped; 6 of the 7,964 polygons are 87,150–97,182
     km²); "km" → "km²"; double full stop. In code the rule is "at least 30 days" (`>= 30`) with
     Q strictly `> 0.0001`, counted per CALENDAR year at compilation (water years only from the
     per-year screening onwards). Measured 2026-10-21 (census §10): of the 8,980 excluded
     candidates 8,787 failed the record-length pre-check before any flow value was examined,
     113 were never retrieved and 80 failed the per-year rule — 67 of those lack 20 years with
     30 days of ANY positive flow, 5 fail only the 0.0001 mm/day level, 4 were lost to the
     qualifier-code mask ("A R" revised values set to NA), 3 to record length, 1 now passes.
     Suggested: "…at least 20 calendar years with at least 30 days of flow above 0.0001 mm/day,
     yielding 8,014 gages; 7,964 of these have a watershed polygon ≤ 100,000 km²".
  2. §2.2.1 still "121 streamflow and hydroclimate signatures … eight signature families" (the
     "100 annually resolved … 21 stand-alone" clause was cut) while the §2 preamble says nine —
     2026-10-12 item 2 unchanged; once the version-2 products are deposited the numbers are
     139 annual + 21 scalars = 160 signatures, nine categories, 2,277 columns.
  3. §2.1.4: Annual NLCD now carries "(CITE)" — add the Annual NLCD Collection 1 reference; the
     MODIS 500 m / composite-interval sentence was cut (no loss). "Falcone 2011; 2017" collapsed
     to reference 22 alone — check the 2017 GAGES-II update is still cited for the
     human-interference indicators (§3 Resource 3).
  4. Carried, all still present: §2.1.1 seasons Jan–Mar / Apr–Jun / Jul–Sep / Oct–Dec (code:
     DJF / MAM / JJA / SON); §2.1.3 "1980–2023", gdptools; gap-rule wording; in-window 20-year
     floor; the 80 % seasonal rule and the 12-flag sentences; §4 quartiles are product #2's;
     reference entries and typos (2026-10-12 item 7).
- **Verified consistent:** the numbered citations map one-to-one onto the previous author–year
  citations where checked (1 Lin, 20 DeCicco, 21 Albers, 22 Falcone, 23 ECCC 2016, 24 Lehner &
  Grill, 25 Linke, 26 Thornton, 29 Myneni, 31 McMillan, 32–33 Petersky & Harpold / Hatchett,
  34 ECCC 2025, 35 Lehner 2022, 36 Friedl).
- **Code/docs-side:** the qualifier-code mask finding (USGS "A R" = approved-revised days set to
  NA by `generate_streamflow_dt`; 4 candidates lost, footprint in the compiled records in census
  §10) is a legacy-ingestion item logged in CHANGELOG → Known Issues, not a manuscript item.

**2026-10-12 — co-author edits synced (109 vs 110 paragraphs; nine regions); reconciliation
pass.** Snapshot overwritten. The edits track the guidelines' new §3.9 Climate Metrics (same
day): the §2 preamble now lists NINE signature families. Every changed methods claim was
checked against the code, the config and the delivered products. Findings by direction of fix:
- **Manuscript wrong, code/docs right (relay):**
  1. **§2.1.1 now defines the seasons as "Winter: January–March, Spring: April–June, Summer:
     July–September, Fall: October–December"** — the water-year quarters. The code
     (`config/signatures_config.json` → `na_handling.seasonal_completeness.season_definitions`,
     `CFG_NA_SEASON_MONTHS`, `SEASONS` in `flow_volumes.jl` / `runoff_ratios.jl` / `pulses.jl`),
     the guidelines Part 1.1 and the manuscript's own previous §2.2.1 sentence all use
     Winter Dec–Feb, Spring Mar–May, Summer Jun–Aug, Fall Sep–Nov. The correct sentence was
     deleted from §2.2.1 in the same edit, so the manuscript now states only the wrong
     definition. Restore DJF / MAM / JJA / SON.
  2. **Family count and signature count disagree within the text**: §2 preamble "121
     signatures across nine signature families (…, snow, storage, climate)" vs §2.2.1 "eight
     signature families" and "121 … 100 annually resolved metrics … 21 stand-alone". No climate
     family exists in the code or the products yet (guidelines §3.9 is a request, logged
     2026-10-12 in `guidelines_todos.md`); if it is implemented before submission the counts
     become ~160 annual signatures (+39 climate, +5 snowfall; final numbers at implementation)
     and both passages, Table 1 and the data dictionary change together. Until then: eight.
  3. **§2.1.3 Daymet**: the sentence explaining the 6,087 / 6,041 gap was removed (good — the
     reprocessed input covers every gage) but the paragraph still says "covering 1980–2023",
     "gdptools", and keeps no statement of coverage. If the products are rerun on
     `daymet_1980_2025_29sep2026.parquet` (planned, 2026-10-12) the paragraph must say
     1980–2025, 8,017 basins (all gages with a boundary + 53 > 100,000 km²), exactextract
     true-area weights over the full-resolution polygons, masked mean over valid cells
     (`EO_data_processing/README_DAYMET.md`); if not, "6,087 basins (5,965 usable gages)"
     must come back.
  4. **§2.2.1 and §2.2.2 dropped two rules the code still applies**: (a) "any WY-season with
     fewer than 80 % of raw observations passing QA → seasonal metrics NA" — implemented
     (`seasonal_completeness.min_fraction` 0.80 on raw observations; applied to Qwin/Qspr/
     Qsum/Qfal and the four seasonal runoff ratios; the seasonal flow reversals are NOT
     masked, a code nuance worth one clause); (b) "We also computed 12 automated data-quality
     flags" — implemented (`compute_qa_flags`, 12 `flagged_for_*` columns in both products).
     Either restore both sentences or state deliberately that they were cut.
  5. **§4 size quartiles describe product #2 only**: "median 698 km², 25 % > 2,700 km², 25 % <
     180 km²" (190 → 180 in this edit) = product #2's `basin_area` quartiles 184 / 698 / 2,718
     (6,222 of 6,250 gages carry an area). Product #1 gives 166 / 645 / 2,500 (6,646 of 6,678).
     Say "for the WY 1980–2025 product" (or report both). "(Figure 2)" added — fine.
  6. §5.1.2 "joined on the leading-zero-stripped form for each dataset … (see resource README
     for details))" — double closing parenthesis; otherwise consistent with
     `docs/plans/2026-09-04-dataset-join-guidance.md`.
  7. Still open from 2026-09-29: the gap-rule wording ("those days were set to NA" — the
     year is rejected); the 20-year floor is counted in-window; 85,000 vs 100,000 km²;
     cross-references after the §2.1.1/§2.1.2 swap; title/abstract naming; citations
     (Albers 2026 vs 2017, Myneni DOI, Hatchett year, Condon 2010 → 2020, "htpps://");
     "Claude Code 0.145.0"; §5.1.2 "0103500".
- **Verified consistent:** the §2.1.4 edit (dropped "(average ~ 110 km²)") changes nothing
  of substance; the funding sentence (NSF EAR Postdoctoral Fellowship 9500317412, DSM) needs
  no check.
- **Code/docs-side: none from this sync.** The climate family itself is tracked in
  `guidelines_todos.md` (2026-10-12) and `docs/STATUS.md` → Pending user decisions.

**2026-09-29 (later) — data finding from the Daymet single-year test (relay; no doc change
yet).** Five of the 6,087 stale Daymet basins are NaN on every day of 1980–2023 (02234500,
02236000, 02236125, 02244040, 01372058), so §2.1.3's "6,087 basins" holds 6,082 with data and
the "5,965" streamflow + Daymet + MODIS gages hold 5,961 with a usable climate series. Moot if
the Daymet reprocess lands before submission (every count in §2.1.3, §3 Resource 3 and the
abstract changes then — action plan §0); otherwise add the caveat.
- [ ] relay only if the reprocess misses submission — D1–D4 decided 2026-09-29 (full-resolution
  polygons, all six variables, + the 53 basins > 100,000 km² if RAM allows); if it lands, §2.1.3
  also changes method (own exactextract area weighting), basin count and years (1980–2025)

**2026-09-29 — MAJOR co-author revision synced (110 vs 105 paragraphs); reconciliation
pass.** Abstract, §4 Data Overview, Data/Code availability and Funding sections added;
front-matter to-do lists removed; §2.1.1 (streamflow) and §2.1.2 (boundaries) swapped;
references reformatted and extended to 39 (Albers 2017, ECCC HYDAT release 2025-10-14 and
the 2026-09-10 HYDAT relay all landed; Arsenault, Knoben, Kratzert, McMillan, Omernik,
Petersky, gdptools 0.3.11 added). Snapshot overwritten. Findings by direction of fix:
- **Manuscript wrong, code/docs right (relay):**
  1. §2.1.1 gap rule now reads "if any data gap exceeded three consecutive days those
     days were set to NA" — the preprocessor REJECTS the whole water year (guidelines
     Part 1.2 items i–ii; `preprocess_daily_data`); the 2026-09-04 wording was correct.
  2. §2.2.2 "at least 20 qualifying water years across its full period of record" — the
     20-year floor is counted WITHIN the analysis window (`run_julia_benchmark.jl`
     applies `length(valid_years) >= 20` after the window filter; census 2026-09-04: 412
     gages with ≥ 20 full-record years are excluded from product #1 by the in-window
     floor). Drop "across its full period of record".
  3. §2.1.2 "We excluded basins that exceeded 85,000 km2" and the new pre-§2.1.1
     paragraph "7,964 watersheds had polygons and were smaller than 85,000 km" — the
     boundary layer (Resource 4) excluded the 54 basins > 100,000 km²; 85,000 km² was the
     Daymet aggregation cap only (58 success gages exceed 85,000 km² per the metadata, so
     the 7,964 include basins between the two thresholds). Also "km" → km². "The total
     number of gages with streamflow, Daymet and MODIS LULC is 5,965" — 5,965 is
     streamflow ∧ Daymet; the ∧ MODIS count needs checking against Resource 5.
  4. "6,041" is back (abstract "6,041–7,964"; §2.1.3 "6,087 or 6,041") — settled
     2026-09-04 as 6,087 (distinct `site_id`s in the parquet).
  5. §2 preamble "8,014 … usable daily records for watersheds smaller than 85,000 km2"
     conflates the gage count with the size cut; its cross-references pre-date the swap:
     "(Sect. 2.1.1)" for boundaries → 2.1.2, "(Sect. 2.1.2)" for the trend windows → 2.2.1;
     §3 "described in Section 2.1.2" → 2.2.1. The agency-drainage-area normalization
     sentence dropped from the preamble survives in §2.1.1 and §5.1.1 (fine).
  6. Title line "Hydroclimate Information, Signatures and Summary Statistics" vs
     "Hydrologic Information Signatures and Summary Statistics" in the abstract, §1,
     Figure 1 and Table 1 captions.
  7. Citations: §2.1.1 still "Albers et al., 2026" (list has Albers 2017); Myneni — the new
     reference is the MOD15A2H (Terra 8-day) DOI while the product used is MCD15A3H v061
     (DOI 10.5067/MODIS/MCD15A3H.061), and §3 cites "Myneni et al. 2015" vs §2.1.4 "2021";
     Hatchett has no year; "Condon … (2010)" → 2020; "htpps://"; "licenceCC-BY 4.0
     license"; "Chen … & Alejandro N. Flores, A. N."; several trailing-punctuation slips.
  8. Still open from earlier passes: "Claude Code 0.145.0" (item 7), §5.1.2 seven-digit
     "0103500", "may be read" → "must", unclosed parenthesis, §5.1.3 "against" →
     "between".
- **Verified consistent:** §2.2.2 constant-flow rule ("any calendar month in which at
  least 15 days of non-zero streamflow held at a single constant value") = config
  `constant_sd_flag` (`min_nonzero_days_per_month` 15, `max_unique_values` 1); the
  eight-family list in §2 and §2.2.1 = `docs/signature_categories.csv` (the manuscript
  says "signature families" where the repo says "categories" — vocabulary only, but the
  repo reserves "function family" for the 15 computing functions; user to decide whether
  to align); the Daymet-gap explanation in §2.1.3 ("processed before the streamflow time
  series requirements were finalized") matches the 2026-09-29 finding that basin size
  explains almost none of the 2,049-gage hole (see `docs/plans/2026-09-29-daymet-
  reprocessing-options.md`).
- **Not verifiable locally:** §4's 12 Level I ecoregions / 62 % / median 698 km² /
  quartiles 190 and 2,700 km² for the product sets (the 8,014-gage set gives 159 / 613 /
  2,403 km², so a larger-basin product subset is plausible); check against the product
  CSVs on the drive.
- **Code/docs-side: none.** Note for later: §2.1.3 and Resource 3 will need rewriting if
  the Daymet input is reprocessed (plan above).

**2026-09-28 — Daymet selection rule RECOVERED (co-author question: is the 7,964-polygon vs
Daymet gap all processing?).** Mostly yes. The 6,087-site request list came from the
ORIGINAL compilation (`combined_watershed_metadata_og.csv`, in git at `27c38d9`), run with
`min_num_years = 30` before the final 20-year recompile: list = 5,956 og-success + 131
og-processing. Of the 1,995 polygon gages without Daymet, **1,952 failed that 30-year
compilation** (1,910 at the retrieval pre-check `nrow > 365 × 30`, 42 on qualifying years);
**6** exceed 85,000 km²; **22** are 9–10-digit USGS ids (none of the 66 such usable gages has
Daymet — the leading-zero id loss); **15** are HydroBASINS-fallback polygons (none of the 29
has Daymet). "Stricter" = record length at compilation only (≥ 30 vs ≥ 20 years); per-year
QA and trend gates are identical. Relay: the manuscript §2.1 placeholder's "too large to
aggregate" explanation → "gages added when the compilation's minimum record was lowered from
30 to 20 years". Correction: the census §7 "all 5,707 golden gages have Daymet" evidence is
circular (the R path iterates over Daymet sites). Record: census §9.

**2026-09-18 — §2.1 "signatures + Daymet + MODIS/LULC = 6,041" checked (user request):
NOT supported by any file in the deposit (relay).** Measured on the staged deposit, ids
joined zero-stripped: signature gages with Daymet AND MODIS LULC = **5,517** (WY
1993–2025), **5,638** (WY 1980–2025), **5,665** in either product, 5,490 in both; 5,965 if
"signatures" is read as the 8,014 usable-streamflow gages. 6,041 exceeds the hard ceiling
of 5,969 (Daymet sites that have a watershed polygon, hence MODIS): 122 of the 6,087 Daymet
sites never reached the usable set (109 `processing`, 13 `no_data`), and all 6,087 are
complete (16,060 rows each), so no "complete-Daymet" subset yields it either. It is the
same untraceable figure that appeared in §2.1.3 in August ("6,041 basins via gdptools") and
in the R3 row of the resources sheet on 2026-09-15 — presumably a count from the
co-authors' upstream gdptools run, not from any delivered file — and it contradicts the
manuscript's own §2.1.3 (5,965 / 5,517 / 5,638). Suggested fill: 5,517 (product #1, the
window the Abstract's trends refer to), or 5,665 across both products. Census record:
`docs/plans/2026-09-04-filtering-stage-census.md` §7.

**2026-09-18 (later) — §2.1.1 boundary sentence checked (user request).** "54 excluded … 7,964
remaining" is right but the base is **8,018**, not 8,014: the boundary/EO layers carry four
`processing`-status USGS gages (01591000/400/610/700 — GAGES-II polygons, retrieval never
completed, no streamflow or signatures). Verified on the staged deposit: the 54 absent
`success` gages are all Canadian, 49 with published area > 100,000 km² and 5 delineated
> 100,000 km² with no published area; **no gage was dropped for lacking a boundary** (31
Canadian fallbacks → 29 survive). Relay: "shapefile" → "layer" (GeoPackage/GeoParquet,
WGS84); drop the claim that climate was aggregated with these boundaries "only where we had
sufficient streamflow data coverage for trends analyses" — Daymet used the earlier 6,087-basin
list (≤ 85,000 km²; 118 sites have no polygon in the layer, 122 never compiled; selection rule
unrecorded). Replacement sentence delivered in-session.

**2026-09-15 (later) — resources-table sheet reviewed (user request; the manuscript's
Table [x]).** Column C "Data Rows" checked file-by-file against the staged deposit (full
count table: `docs/plans/2026-09-04-filtering-stage-census.md` §8). Correct: R1 6,678, R2
6,250, R5 7,964 (MODIS; NLCD is 6,119). Wrong (relay): R3 "6,041" matches no file — it is
8,014 gages / 111,624,189 daily rows of streamflow, 16,994 metadata rows, 6,087 Daymet
sites / 97,757,220 rows, 8,012 HYDAT rows; R4 "8,014" is the HydroATLAS attribute table,
the boundary layer has 7,964; the notes "8,044" → 8,014, "polygons < 85,000 km²: 8,014" →
7,964 polygons (exclusion at 100,000 km²; 7,958 of them ≤ 85,000 km²), "with daymet and
modis: 6,014" → 5,969 (5,965 of the 8,014; 5,517 / 5,638 within the products). The
column mixes gage counts with row counts; suggest a "Gages" column plus per-file rows.
**Re-checked after the user's edit the same day: every value now matches the deposit**
(column renamed "Number of Gages"; R3 split 8,014 / 6,087, R4 8,014 / 7,964, R5 7,964 /
6,119; notes 8,014 / 7,964 / 5,965). Only nuance left: the R3 metadata table lists all
16,994 candidates, and the 8,012-row HYDAT interference table is not named in Contents.

**2026-09-15 — co-author revision synced; reconciliation pass (user request).** Changes:
new co-author (Dymond, S.F.); a new Abstract; §1 rewritten; §2 preamble and §2.2.1 now
carry the EIGHT categories (category-edit catalogue items 1, 2 and 8 RESOLVED; item 5, the
§4 figure scheme, still unknown); a new placeholder paragraph under §2.1; §3 join sentence
revised; a Figure 2 caption. Snapshot overwritten. Every changed number was checked against
the code, the config and the staged deposit files (`~/Downloads/Signatures`, joins on
zero-stripped ids; script record in `docs/plans/2026-09-04-filtering-stage-census.md` §7):
- **WRONG (relay) — gage counts in the Abstract and §1.** ">8,000 North American
  watersheds" / ">8,000 stream gages" / "For over 8000 watersheds … we calculate roughly
  100 hydrological streamflow signatures": 8,014 is the number of gages with USABLE DAILY
  STREAMFLOW (Resource 3); boundaries and EO tables cover 7,964 (below 8,000); signatures
  and trends exist for 6,678 (WY 1993–2025) and 6,250 (WY 1980–2025); Daymet for 5,965 of
  the 8,014. "Trends (1993 to 2025)" describes product #1 only.
- **WRONG (relay) — the new §2.1 placeholder paragraph, three numbers.** (a) ">30 days
  over .001mm/day": the compilation gate is ≥ 30 days per calendar year with Q >
  **0.0001** mm/d in ≥ 20 calendar years (`config.R` `MIN_Q_VALUE_AND_DAYS = c(0.0001,
  30)`, `MIN_NUM_YEARS = 20`; strict `>` in `R/helperFunctions.R`). (b) "smaller than
  8,500 km (about average ~ HUC8)": the manuscript's own §2.1.3 says **85,000 km²** (the
  LARGEST HUC8); the geometry exclusion is 100,000 km²; "km" → "km²". (c) The XXXX fill-ins,
  measured on the staged deposit: gages with signatures AND Daymet AND MODIS LULC =
  **5,517** (WY 1993–2025) / **5,638** (WY 1980–2025) — every Daymet gage also has MODIS,
  so this equals §2.1.3's Daymet counts; signatures only 6,678 / 6,250; MODIS but no
  Daymet 1,117 / 567; NLCD 5,454 / 5,007; all four 4,464 / 4,524. At the 8,014 level:
  MODIS 7,960, Daymet 5,965, both 5,965.
- **WRONG (relay) — "additional gages … too large to reasonably aggregate climatic
  metrics".** Of the 1,161 product-#1 gages without Daymet only **6** have polygon area
  > 85,000 km² (46 by agency drainage area; median 230 km²). The cause is that the
  6,087-basin Daymet request list predates the Feb 2026 compilation: all 5,707 gages of the
  Feb 2026 R golden set have Daymet; Daymet gages are the long-record ones (median 46 record
  years, 75 % starting 1979) and the no-Daymet gages the shorter/later records admitted by
  the 20-year rule (median 27 years, median start 1992). Suggested wording: Daymet was
  aggregated for the 6,087 candidate basins (≤ 85,000 km²) of the earlier compilation; the
  ~2,000 gages added by the final 20-year compilation carry no climate series.
- Minor (relay): §2 preamble category order differs from §2.2.1; "Precipitation-Streamflow"
  (hyphen) vs the sheet's "Precipitation Streamflow"; "(FIG REG TO WORKFLOW SCHEMATIC)"
  placeholder — the schematic is `docs/plans/dataset_workflow_schematic.png`; §3
  Conventions now states canon_id twice (the new "join with leading zeroes removed" clause
  makes the older canon_id sentence redundant); Figure 2 "ecoregion" — HISSS carries no
  ecoregion column except the HydroATLAS majority classes `tec_cl_smj` / `fec_cl_smj`
  (Resource 4); name the ecoregion layer in the caption.
- RESOLVED by the co-authors' edits: ">6,000stream", "7, 964", "wtih", "Linke et al.,
  2013" (paragraph deleted), the unclosed "(e.g., CAMELS-Chem," parenthesis (sentence
  deleted). Still open from earlier passes: §5.1.2 "0103500" (seven digits) and "may be
  read" → "must"; §5.1.3 "compared against" → "between"; "s(e.g.,"; "[tbd]" citations;
  "fitted checksums …" garble; "Claude Code 0.145.0" (item 7); references (item 5);
  CC-BY/DOI labelling (item 11); §2.1.2 end date for Canadian gages (HYDAT, 2026-09-10).
- Code/docs-side: none — no methods claim changed. The eight-category wording matches
  `docs/signature_categories.csv`.

**2026-09-10 (evening) — HYDAT citation check (user request).** The manuscript's
"tbd: ECC hydat citation" and the proposed reference "ECCC (2026). HYDAT: National
Hydrometric Database [2026-07-17]" name the WRONG release: 2026-07-17 is the release
current today (the only one WSC serves), five months after the 2026-02-07 retrieval. The
release actually used is **2025-10-14** (inferred: tidyhydat 1.0.0's CRAN build of
2026-02-03 still queried "HYDAT released on 2025-10-14", the next release was 2026-04-17,
and the compiled Canadian records end 2024-12-31 with no 2025 data, which the 2026-04-17
release already carries). Definitive confirmation = `tidyhydat::hy_version()` on the
Windows ingestion machine. Relay: (1) cite the 2025-10-14 release with a 2025 year and the
February 2026 retrieval date; (2) `Albers et al., 2026` → Albers, S. (2017), JOSS 2(20),
511, doi:10.21105/joss.00511 (single author; add the package version if desired — 0.7.2 or
1.0.0, whichever was installed on 2026-02-07); (3) §2.1.2 "retrieved … to 30 September
2025 for all candidate gages" holds for USGS only — Canadian records end 31 December
2024, so **no Canadian gage has a qualifying WY 2025 in either product** (end_water_year
≤ 2024; 813 / 772 Canadian gages end in WY 2023). Recorded in docs/DATA_SOURCES.md (HYDAT
row); HydroShare R3 README/dictionary need the same statement (Planned → HydroShare
documentation updates, H3).

**2026-09-10 (pm) — re-synced: manuscript unchanged.** **Catalogue of the eight manuscript
locations that must change to the eight categories (user decision, later the same day):
`docs/plans/2026-09-10-manuscript-category-edits.md` §A.** Category crosswalk against the
restructured guidelines doc and the staged HydroShare deposit: the manuscript's 14-category
wording (§2 preamble, §2.2.1, §4 figure) is now the outlier — see Guidelines Document
TODOs (2026-09-10 pm) and `docs/plans/2026-09-10-signature-category-review.md`.

**2026-09-10 — sync; no methods change.** Since the 2026-09-04 second sync: §3
Resources 1–2 now reads "21 signatures and related outputs that do not carry statistics"
(was "per-gage scalar diagnostics"; consistent with §2.2.1's "21 stand-alone diagnostic
metrics"); **§5.1.2 "How to join HISSS with other datasets" landed** (the
`docs/plans/2026-09-04-dataset-join-guidance.md` §4 paragraph with co-author edits) and the
record-dependent note moved to §5.1.3; the Acknowledgements double period is fixed. Snapshot
overwritten. Relay (manuscript-side): §5.1.2's example id "0103500" has seven digits (USGS
site numbers have eight or more — e.g. 01013500); "may be read as a character string" should
be "must"; an unclosed parenthesis "(most US site numbers begin with a zero … (see resource
READMEs for details)."; §5.1.3 "should not be compared *against* the WY 1993–2025 and
WY 1980–2025 products" → "*between*"; new §2.1.1 typo "7, 964". Still open from earlier
passes: avg_storage sentence (item 6), "Claude Code 0.145.0" (item 7), references +
"Linke et al. 2013" (item 5), CC-BY/DOI labelling (item 11), §2.1.2 end date,
"(e.g., CAMELS-Chem,", ">6,000stream". **Category consistency**: §2.2.1's 14 categories
match the code and README exactly; the guidelines doc has 13 modules and the shipped
dictionary a 9-class grouping — see Guidelines Document TODOs (2026-09-10) for the
proposed two-axis resolution; if the §4 "sites and signature categories" figure uses the
colleague's 8 groups, §2.2.1 or §4 needs one sentence saying the 14 categories are
aggregated into 8 analysis groups for presentation.

**2026-09-04 — filtering-stage watershed census + fifth reconciliation pass (user
request: count the watersheds at every stage of the workflow schematic and confirm the
manuscript matches).** Full record: `docs/plans/2026-09-04-filtering-stage-census.md`;
new tools `docs/benchmarks/qualification_census.jl` (runs the canonical preprocessor on
ALL 8,014 gages and re-implements the runner's inclusion gates — reproduces both products'
gage sets EXACTLY, 0 mismatches) and `summarize_qualification_census.py`. The manuscript
was revised by the co-authors since 2026-09-01 (snapshot overwritten): correction blocks
1a/1b/1c/2a/2b/3/4/5/6 are APPLIED (121 signatures / 14 categories incl. drought,
agency-published drainage areas, Daymet coverage 5,965 / 5,517 / 5,638, duplicate BFI
paragraph deleted, avg_storage listed, new Usage Notes 5.1.1 + 5.1.2, HISSS repo link,
"eight" MODIS schemes); block 7 (Acknowledgements) is NOT — "Claude Code 0.145.0" still
appears twice. **Every checkable count in the published text is CORRECT** (16,994 =
9,154 + 7,840; 8,980 excluded; 8,014 = 6,160 + 1,854; 111,624,189; 73 / 32 / 28;
6,678 / 6,250; 54 / 7,964; 100 % / 98.3 %; 6,087 / 5,965 / 5,517 / 5,638; 6,119;
1,653; 18,898,406 / 24,366,487; 2,150,280; 191,136; 250,879). Measured stage counts not
in the text, for the figure/Data Overview: 33,732 of 317,182 gage-years (10.6 %)
rejected by the preprocessor over WY 1980–2025 (31,091 >30 NA days, 1,803 gap >3 d, 838
boundary NA); within the products only 3.2 % / 2.9 % of gage-years are rejected;
exclusions 1,336 (924 both gates, 412 only <20 yr, 0 only <60 %) and 1,764 (1,063 only
<60 % with ≥20 valid years, 674 both, 27 only <20 yr); trend statistics survive for
6,183 / 5,851 gages on the dense signatures. Findings, by direction of fix:
- **DEPOSIT DEFECT (staging, fix before upload)**: `gage_id` is NOT zero-padded for 9
  gages in Resource 4's boundary layer (+ its QA CSV) and 44 gages in each Resource 5
  table (`zfill(8)` cannot restore the leading zero of 9–10-digit USGS ids), so the
  manuscript's and the READMEs' "join directly on gage_id" claim fails for those gages
  and the R5 README's coverage counts (6,599 / 6,196 MODIS, 5,419 / 4,998 NLCD — also
  quoted in this file's 2026-08-25 R5 entry) undercount; canonical values are 6,634 /
  6,205 and 5,454 / 5,007. `canon_id`, HydroATLAS and Resource 3 are clean. User decision same day: keep the
  files, document the strip-on-both-sides join rule in both repos' README and the
  resource READMEs (done — see Known Issues).
- Manuscript (relay): §2.1.2's "8,980 … fewer than 20 total years" bucket includes 113
  USGS gages whose retrieval never completed (`processing` status; 4 of them carry
  GAGES-II polygons and sit in the EO layers) — add "or could not be retrieved"; the
  year-rejection rule list omits the third rule (≤3-day NA run touching Oct 1 / Sep 30
  is not interpolated ⇒ `residual_na`, 838 gage-years); §3's "211 watershed-scale
  attributes" = 211 columns (198 attributes + 13 keys/diagnostics); new typos "wtih",
  "diagnostics.We", plus the surviving "s(e.g.,", "MOIDS", "fitted checksums …" garble,
  "Linke et al., 2013", "[tbd]" citations and the missing references (McMillan,
  Hatchett, Petersky & Harpold, Adelsperger, Laaha, Peters & Aulenbach, HYSETS, Caravan,
  CAMELS-SPAT, DeCicco, Albers, gdptools, Annual NLCD).
- Code/docs-side: none. The schematic is not contradicted by any measured count.

**2026-09-04 — co-author revision synced; dataset-join guidance drafted (§5 stub).**
The published draft changed since 2026-09-01: §2 preamble and §2.1.2 now state the
agency-published drainage areas as the normalization source (queued #1 core error
RESOLVED; the §5 caveat also landed as 5.1.1); counts corrected to **121 signatures /
14 categories incl. drought** and the "(name)" placeholder resolved to
`hisss_data_dictionary.csv` (queued #8 + open items 1a/1b RESOLVED; the drought
*methods* paragraph 1c is still absent); §2.1.3 Daymet coverage numbers added
(6,087 candidates / 5,965 of 8,014 / 5,517 and 5,638 — matches the repo); the 47 %
sentence replaced by the explicit "8,980 of 16,994 excluded" (item 3 RESOLVED); §2.2.1
/ §2.2.2 rewritten and the duplicate BFI-statistics paragraph deleted (item 4
RESOLVED); §2.1.4 "all five" → "all eight" (the 2026-09-01 finding RESOLVED); §3
mechanical typos fixed (`wy{window}`, WGS84, teh/ot, "[tbd])"); `<LINK TO REPO>` filled
with https://github.com/CZ-Sync/HISSS; **§5 Usage Notes drafted** — 5.1.1 unnormalized
gages and 5.1.2 record-dependent signatures (both verbatim from block 2b), with the
dataset-join item still a one-line stub. Still open: avg_storage sentence (item 6),
Acknowledgements "Claude Code 0.145.0" (item 7), missing references + "Linke et al.
2013" (item 5), CC-BY/DOI labelling (item 11), §2.1.2 end date, "(e.g., CAMELS-Chem,"
unclosed paren, ">6,000stream".

*Dataset-join guidance (user request, same day)*: CAMELS-Chem, MacroSheds, and the
Daymet-VPD product were reviewed against their own documentation and files, the
identifier conventions of HYSETS / CAMELS / CAMELS-SPAT / Caravan / CANOPEX / GAGES-II
/ WQP / NLDI→COMID / HydroBASINS were verified, and the actual join behaviour of the
five staged HydroShare resources was audited. Measured overlaps: CAMELS-Chem 507 / 516
(WY 1993–2025) and 515 / 516 (WY 1980–2025), all GAGES-II Ref; MacroSheds stream
gauges 146 / 224 inside ≥ 1 product watershed but almost always *nested* (median
smallest containing polygon 799 km²), only 11 co-located within 500 m. The audit found
two deposit defects (Known Issues, HIGH: mis-padded `gage_id` in R5/R4; R3 README
recipe). Verified facts, numbers, and the first-draft §5.1.3 paragraph:
`docs/plans/2026-09-04-dataset-join-guidance.md`.

**2026-09-04 (second sync, afternoon) — user edits; manuscript otherwise unchanged.**
Two typo fixes only (a missing space after "diagnostics." in the §2 preamble; "MOIDS" →
"MODIS" in §2.1.4); no methods claim changed, nothing to reconcile. The §5 "How to join
with other datasets" item is still the one-line stub — the paste-ready §5.1.3 text
(revised to the strip-leading-zeros-on-both-sides join rule after the user's decision
to leave the staged deposit files as delivered) is in
`docs/plans/2026-09-04-dataset-join-guidance.md` §4. The snapshot rebuild also removed
~170 lines of Google's page JavaScript that the morning rebuild had appended below the
body (the extractor did not skip a `<script>` block inside the contents div; the
paragraph-level diff was unaffected because the JS lines had no counterpart on either
side).
Later the same day the user asked for a replacement for §3's "Conventions and key
fields" join sentence (their interim edit — "All tables join on a zero-stripped gage_id,
the zero-padded agency station identifier …" — contradicts itself); two variants were
delivered and recorded in `docs/plans/2026-09-04-dataset-join-guidance.md` §4b. A
six-minute re-poll of the published copy showed no change, so the snapshot stands and
the next sync will confirm what landed.

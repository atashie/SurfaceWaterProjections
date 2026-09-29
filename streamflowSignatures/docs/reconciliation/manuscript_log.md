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

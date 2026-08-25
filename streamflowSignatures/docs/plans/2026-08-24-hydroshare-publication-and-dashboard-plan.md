# HydroShare Publication + CUAHSI-Hosted Dashboard Plan

**Date**: 2026-08-24 · **Status**: ACTIVE · **Owner**: Arik (+ HISSS working group)
**Drives**: (1) publishing the HISSS data on HydroShare, (2) a CUAHSI-hosted dashboard,
(A) the S3→local data staging pass, (B) finishing the Scientific Data manuscript.
**Hard deadline**: Scientific Data submission **Nov 9, 2026**.

Background: the Aug 24 repository-policy review (HydroShare vs Zenodo, verified live;
report artifact `claude.ai/code/artifact/d0a75c48-94cd-4a88-a745-06cac89cce09`) concluded
HydroShare primary for data + Zenodo for the code DOI. CUAHSI has since replied to the
PIs' hosting inquiry (email received before 2026-08-24), with three offers and one caution:

1. **CUAHSI shinyApps server** — can host **beyond the life of the grant**, but the app
   **must be Shiny**.
2. **Duration-limited app hosting** on CUAHSI's Google Cloud billing — needs design +
   estimated hosting costs up front.
3. CUAHSI could **help create a new view of the dataset** as a HydroShare showcase.
- Caution: data preservation on HydroShare is long-term; **apps are not** — build in
  **sunsetting expectations**. Also: the CZHub award got a second NCE → operations to
  **end of August 2027** (little funding in the final year); CUAHSI can cover publication
  costs + minor hosting fees.

## ⚠ Update — 2026-08-24 (later same day)

Two conditions changed after this plan was written; the edits are folded into the
sections below (marked), this block is the summary:

- **S3 access is LOST** (both project buckets). Most delivered data were backed up
  first to the project Google Drive:
  https://drive.google.com/drive/u/0/folders/1DVuq4nC5j_Y01sBaDj9cwjbv7S8sndjj.
  Workstream A is now a **Drive-based** staging pass and begins with an inventory
  verification ("most" ≠ all). The S3 option for the dashboard backend no longer
  exists — D3 is settled by events. See CHANGELOG → August 2026.
- **Dashboard: scope set, work DEFERRED.** The original Shiny app is judged **very low
  performance** — plan a performant rebuild, not a port. The final app must showcase
  (i) streamflow signatures and metrics, (ii) annualized values, (iii) raw time series
  values, and (iv) per-watershed MODIS LAI/LULC + Annual NLCD annual/monthly values
  and metrics (aggregated tables only — **no raw LAI/LULC**). Workstream 2 is
  **decoupled from the Nov 9 track** and will be revisited after data publication; the
  CUAHSI reply (D1/D2) should still go out to hold the hosting option.

## Decisions (proposed — confirm with group, then reply to CUAHSI)

- [ ] **D1 — Hosting option: choose #1 (CUAHSI shinyApps server), and offer #3 as a
      collaboration on top of it.** Rationale: we already have a working Shiny app
      (`streamflowAndClimateVisualizationApp/app.R`, 1,738 lines, leaflet + plotly), and
      option 1 is the only offer that outlives the grant. Decline #2 (GCP) unless the
      Shiny server can't take the app — it has the shortest lifetime and the most
      up-front cost engineering.
- [ ] **D2 — Sunset terms to propose**: app runs through the promotion/discovery phase;
      reviewed annually with CUAHSI; minimum lifetime through **2 years
      post-paper-publication**; at sunset, the HydroShare resource landing pages gain a
      note pointing to the archived self-contained HTML explorers (which are deposited
      IN the resources and work forever as downloads). This makes the app disposable
      by design — the durable interactive artifacts are the archived explorers.
- [x] **D3 — Dashboard data backend: HydroShare public-resource URLs, not S3.**
      *(Settled by events 2026-08-24: S3 access was lost, so a non-S3 backend is now
      mandatory, not preferential.)* Once
      resources are public, their files have direct unauthenticated download URLs — the
      app fetches (or ships pre-bundled) compact extracts, and the `aws.s3` dependency +
      credentials disappear. This is also exactly CUAHSI's option-3 story ("how
      HydroShare can be used in similar cases"). The 4.6 GB daily parquets do NOT move
      to the app backend (see Workstream 2).
- [ ] **D4 — Publish timing**: keep resources **private through peer review** (share via
      HydroShare private links; formal deposition required from review round 2), publish
      (immutable + DOI activation) **at acceptance**. The DOI strings
      `10.4211/hs.<resource-id>` are known from resource creation, so the manuscript's
      Data Records section can cite them before they resolve. Trade-off vs
      publish-before-submission: reviewers-requested data changes would force new
      DOIs on published resources; private-until-acceptance avoids that.
- [ ] **D5 — Raw inputs go in the deposit** (streamflow parquet + metadata + rebuilt
      per-basin Daymet parquet, ~4.7 GB): makes the descriptor self-contained and
      reproducible; still well inside the free 20 GB quota. The 44 raw Daymet CSVs
      (10 GB, drive root) are the co-authors' upstream product — DECIDE: include as a
      single ZIP (~2–4 GB compressed) or document the converter + ORNL provenance and
      omit. Default: omit, revisit if reviewers push.
- [ ] **D6 — Resource architecture** (each gets its own DOI; grouped in a published
      Collection): see table in Workstream 1. Confirm the 5-resource split with the
      group — it maps 1:1 onto the manuscript's Data Records table, and coordinates
      with **Nic's "one nice zip" task** in the draft to-do list (that task becomes
      "per-resource packaging" — sync with Nic before they build a monolith).

---

## Workstream 1 — HydroShare data publication

### Resource structure (proposed)

| # | Resource (own DOI) | Contents | Size | Source location today |
|---|---|---|---|---|
| R1 | HISSS signatures — WY 1993–2025 standard product | signatures CSV (159 MB), annual parquet (86 MB), explorer HTML + 100-file `_annual/` sidecar, comparison dashboard + CSV + summary, RUN_NOTES, timing/provenance JSON, validation reports | 434 MB | `/Volumes/Untitled/processedOuts_drought_28jul2026/` |
| R2 | HISSS signatures — WY 1980–2025 standard product | same layout as R1 | 581 MB | `/Volumes/Untitled/processedOuts_1980_2025_11aug2026/` |
| R3 | Input data | `combined_streamflow_data_09feb2026.parquet` (858 MB), `combined_watershed_metadata_09feb2026.csv`, `daymet_1980_2023_rebuilt_10aug2026.parquet` (3.8 GB), `metadata/canadian_hydat_interference.csv` | ~4.7 GB | thumbdrive `processedOuts_feb2026/` + repo; streamflow parquet + metadata also on S3 |
| R4 | Watershed geometry + basin attributes | `watershed_polygons_26jun2026.{gpkg,parquet}` + QA CSV, HydroATLAS metadata CSV/parquet + dictionary, GAGES-II/HYDAT human-interference columns doc | ~120 MB (verify) | geometry on S3; **HydroATLAS product location TBD** (built in R, `data_out/` on the machine that ran it — not on this Mac) |
| R5 | EO companion products | MODIS LAI monthly parquet, MODIS LULC annual parquet, NLCD annual parquet + reconstructed exclusions CSV, merged dictionary + README, LAI/LULC/NLCD explorers | ~395 MB (measured) | **STAGED + verified + user-signed-off 2026-08-25** (`resource5_eo_products/`); NA-basins CSV + granule manifest absent from backup (documented); NLCD human QA complete |
| — | Collection: "HISSS" | links R1–R5 + paper + code DOI | — | title + abstract + keywords drafted 2026-08-25 (`Signatures/hisss_abstract_collection.txt`); create on HydroShare at deposit |

Total ≈ **6–8 GB** → fits the free 20 GB quota, every file under the 25 GB/file limit; no
quota request needed. Code is NOT a HydroShare resource: archived via **Zenodo GitHub
integration** at acceptance (software DOI), after the "move code into the CZ github" step.

### Checklist

**Account + metadata prerequisites**
- [ ] All co-authors create/complete HydroShare profiles (publishing is blocked on
      incomplete owner profiles). Collect HydroShare usernames alongside Kendra's
      authorship table.
- [ ] Decide author ORDER per resource now (authorship is FROZEN at publication;
      draft note says Arik first on the data release, Kendra first on the paper).
- [ ] Write per-resource abstracts (≥ 150 chars — publish gate) — include per-product
      window, gage count, column count, and the **record-dependent-signatures
      cross-product caveat** in R1/R2 abstracts; ≥ 3 keywords incl. geographic
      ("United States", "Canada", "North America"); funding agency entries.
- [ ] Per-resource README.md + data dictionary (column name → definition/units).
      R1/R2 largely exists (SIGNATURES.md, RUN_NOTES.md, claude-skill); needs assembly
      into a standalone README per resource. R3 README must document the Daymet
      365-day-calendar quirk, the rebuilt-parquet provenance, and `area_normalized`.
- [ ] License: **CC-BY 4.0** on every resource (journal prohibits -NC/-SA).
- [ ] Sensitive-content pass: confirm no credentials/keys/emails in RUN_NOTES,
      timing JSONs (provenance blocks carry hostnames/paths — acceptable? review),
      validation logs.

**Resource creation + upload**
- [x] **COLLECTION CREATED (2026-08-25)**: "Hydrologic Information Signatures and
      Summary Statistics" —
      https://www.hydroshare.org/resource/f702201faa5d46069a5ee83ffa4c9768/
      → collection DOI at publication:
      **https://doi.org/10.4211/hs.f702201faa5d46069a5ee83ffa4c9768**
      (verified live 2026-08-25; recorded in memory, CHANGELOG, and all five
      staged resource READMEs).
- [ ] Create the 5 member resources PRIVATE; record the 5 resource ids →
      derive DOI strings for the manuscript; add each to the collection.
- [x] **Stage all files locally — DONE for ALL FIVE resources (2026-08-25)**:
      `~/Downloads/Signatures/resource{1..5}_*/` with deposit names, per-resource
      READMEs, dictionaries, abstracts, all Codex-reviewed (R1–R2, R4 this session;
      R3 + R5 sister session, R5 with the user's NLCD QA sign-off; R4's boundary
      layer is the 2026-08-25 rebuild). Upload path: web UI fine at these sizes;
      `hsclient` for scripted/resumable upload of R3's 3.8 GB file.
- [ ] Verify upload integrity: HydroShare checksums / re-download spot-check of the
      two parquets; confirm folder structure survived (R1/R2 `_annual/` sidecars).
- [ ] Set Geographic (bounding box: CONUS+Canada) and temporal coverage per resource.
- [ ] Cross-link: collection ↔ resources; "related resources" → paper (when DOI known)
      and Zenodo code DOI.
- [ ] Generate a **private link** per resource for the journal reviewer package.

**Publication (at acceptance — not before, per D4)**
- [ ] Freeze check: no pending fixes; RUN_NOTES/READMEs final.
- [ ] Flip public → Publish each resource (expect CUAHSI light staff review), then
      publish the collection.
- [ ] Confirm DOIs resolve; update manuscript proofs if resource ids changed.
- [ ] Tag GitHub release; archive via Zenodo GitHub integration → software DOI
      (add `CITATION.cff` + LICENSE to repo first; note: no DOI pre-reservation via
      the GitHub route).
- [ ] CHANGELOG entry + claude-skill update (publication URLs, DOIs, citation text).

---

## Workstream 2 — CUAHSI-hosted dashboard (Shiny)

**Constraint discovered in review**: the current `app.R` pre-opens Arrow datasets
directly from OUR S3 bucket (`climate-ai-data-science-shiny-app-data`): the 858 MB daily
streamflow parquet, `daymet_1980_2023.parquet`, metadata CSV, a signature summary CSV,
and a precomputed cross-signature CSV. As-is it (a) needs our AWS credentials on
CUAHSI's server, (b) depends on our bucket outliving the app, and (c) leans on multi-GB
daily data a shared Shiny host shouldn't scan.

**Superseded 2026-08-24 (see Update block)**: S3 access is gone, so the current app is
non-functional as deployed and the backend re-point is mandatory. The group also judged
the original app **very low performance** — plan a performant REBUILD, not a slimming
port — and fixed the final app's required scope: signatures + metrics, annualized
values, raw time series values, and the per-watershed MODIS LAI/LULC + NLCD
annual/monthly values and metrics (no raw LAI/LULC rasters). **This workstream is
DEFERRED** — decoupled from the Nov 9 submission; only the CUAHSI reply proceeds now.

- [ ] **Reply to CUAHSI** (PIs; draft below) — choose option 1 + offer option 3,
      propose sunset terms (D2), and ask their shinyApps server specs: RAM per app,
      app bundle size limit, R version, allowed outbound network, deploy mechanism
      (rsconnect?), and whether they want the HydroShare-backed design as a showcase.
- [x] **Define the hosted app's scope** — SET by the group (2026-08-24), superseding
      the earlier drop-the-daily-view recommendation: the app must showcase
      (i) streamflow signatures and metrics, (ii) annualized values, (iii) **raw time
      series values** — required, so dropping the daily view is off the table;
      precomputed per-gage compressed series fetched on click is the leading design —
      and (iv) per-watershed MODIS LAI/LULC + Annual NLCD annual/monthly values and
      metrics (aggregated tables only, no raw LAI/LULC). S3 is no longer an option.
- [ ] **Re-point data access to HydroShare public URLs** (D3): app startup fetches (or
      ships bundled) the R1/R2 signature CSVs + compact per-gage annual extracts;
      no `aws.s3`, no creds. Cache in-app; measure cold-start time.
- [ ] Precompute the app-ready extracts (columns actually visualized, rounded, zstd
      parquet or fst) — target **< 300 MB resident, < 1 GB RAM** under profiling.
- [ ] Local test: `shiny::runApp()` against HydroShare URLs (use a public test resource
      or private-link URLs during development); memory/startup profile.
- [ ] Add an "About / data" page: cites the DOIs, links the resources, states the
      sunset policy, and credits CUAHSI hosting + funding.
- [ ] Deploy with CUAHSI; smoke-test; add the app URL to the HydroShare resource
      landing pages ("app connector" or README link) and to the paper's Usage Notes
      if timing allows.
- [ ] Post-sunset artifact check: confirm the deposited self-contained explorers open
      cleanly from a fresh download (they are the app's permanent fallback).

---

## Workstream A — data staging inventory (what to pull from where)

**Rewritten 2026-08-24: S3 access is LOST.** Both project buckets are unreachable; most
delivered data were backed up first to the project **Google Drive**:
https://drive.google.com/drive/u/0/folders/1DVuq4nC5j_Y01sBaDj9cwjbv7S8sndjj.
"Most" is doing real work in that sentence — the staging pass now BEGINS with an
inventory of the Drive folder against the manifest below.

This Mac currently has: the two product folders + inputs on the thumbdrive; repo
`data_out/` is nearly empty.

- [ ] **Inventory the Google Drive backup** against the manifest below: file names +
      sizes vs the documented delivery sizes; record what is present, what is missing.
      Also check whether the two delivered product folders and the thumbdrive inputs
      were themselves backed up (if not, uploading them to Drive doubles as the
      off-exFAT backup this plan already required).
- [ ] **⚠ FIRST: look for `daymet_1980_2023.parquet` (the ORIGINAL) in the Drive
      backup** — the Shiny app read it from S3, so a copy existed there and may have
      been backed up. If present: compare size to the recorded **4,125,630,653 bytes**
      (28 Jul provenance) and verify the `PAR1` footer. Intact → **resolves the
      CHANGELOG "restore the original" item** (product #1 reproducible, the ≤ 3.4e-13
      rebuild residual attributable, the truncated thumbdrive file deletable; decide
      whether R3 ships original, rebuild, or both with provenance notes). Absent or
      truncated → record it; the rebuild remains the deposit copy, and the original is
      likely unrecoverable (S3 was the last known location).
- [ ] Retrieve from the Drive backup (verify sizes on arrival):
  - [ ] `watershed_polygons_26jun2026.gpkg` + `.parquet` + `_qa.csv` (R4) — while here,
        tally `watershed_geom_source` from the QA CSV to settle the doc conflict
        (README_NLCD's "~19 US HB-fallback" vs README §11's "all 29 HB fallbacks
        Canadian"; manuscript §2.1.1's "100% of US gages" claim rides on it)
  - [x] `watershed_modis_lai_monthly_30jun2026.parquet` + `_dictionary_*.csv` (R5) —
        retrieved 2026-08-25 (`~/Downloads/LAI-LULCC/`), verified, staged. The
        `_na_basins_*.csv` was NOT in the backup — moot: the 17 urban basins are
        derivable from the parquet (documented in `hisss_readme_eo.md`)
  - [x] `watershed_modis_lulc_annual_29jun2026.parquet` + `.csv` + `_dictionary.csv`
        (R5) — retrieved, verified, staged. The `_granule_manifest_29jun2026.parquet`
        was NOT in the backup — documented in the README as repo-level provenance
  - [x] LAI + LULC explorers (R5) — retrieved (as `lai_explorer-30jun2026.html` /
        `lulc_explorer-30jun2026.html`), truncation-checked, staged
  - [ ] `combined_streamflow_data_09feb2026.parquet` + metadata CSV — **no pull
        needed** (intact local copies on the thumbdrive); just confirm Drive holds
        copies as the off-site backup.
- [x] **NLCD — backup CONFIRMED (user, 2026-08-24)**: the finalized product is in the
      Google Drive backup — no re-run needed. Remaining NLCD work joins the normal
      staging flow: retrieve `watershed_nlcd_annual*.{parquet,csv}` + dictionary +
      exclusions CSV + `nlcd_explorer.html` from Drive, and complete the still-owed
      **human QA via the explorer** before the R5 deposit (the old S3 publish step is
      moot; publication target is HydroShare).
      **DONE 2026-08-25**: parquet + dictionary + explorer retrieved, verified
      (finalized build: 250,879 rows / 6,119 gages / C1V2), and staged; the exclusions
      CSV was not in the backup and was exactly reconstructed
      (`hisss_nlcd_excluded_gages.csv`, 45 gages, all verified Alaska). **Human QA via
      the explorer COMPLETE — user signed off 2026-08-25.** R5 staging
      (`resource5_eo_products/`: 3 parquets + exclusions CSV + merged dictionary +
      README + 3 explorers, ~395 MB) is deposit-ready; see CHANGELOG → August 2026.
- [ ] **Locate the HydroATLAS product** `watershed_hydroatlas_metadata_{date}.{csv,parquet}`
      + `_dictionary.csv` (R4): not on this Mac; check the Windows machine's `data_out/`
      and the Drive backup; if lost, re-run `run_hydroatlas_metadata.R` (R-only, needs
      the Windows box or an R env + the gpkg on the thumbdrive).
- [ ] Cross-signature-analysis CSV (formerly read from S3 by the app) — check the Drive
      backup only if the future app keeps that panel; otherwise skip.
- [ ] Assemble a staging directory mirroring the R1–R5 layout; record SHA-256 of every
      file (the thumbdrive is exFAT and has silently truncated a file before — hash
      BEFORE upload, verify after).
- [ ] Back up the staging directory off the thumbdrive — Google Drive is now the
      natural off-site home; the deposit copy must not live only on exFAT.

---

## Workstream B — finish the Scientific Data manuscript

Target journal: Scientific Data (Special Collection on Water Storage, **due Nov 9,
2026**). Draft's internal "done" date (Aug 14) has passed — propose a revised internal
deadline of **Sept 30** for a full draft, co-author review through October.

**Structural gaps (from the live draft)**
- [ ] §2 opening: write the methods summary paragraph.
- [ ] §2 schematic: workflow diagram (**Kendra**) — offer to generate a draft figure
      from the repo's data-flow docs as a starting point.
- [ ] §2.1: input-data table — generate directly from `docs/DATA_SOURCES.md` (11
      sources + access modes already tabulated there).
- [ ] §2.2.1 "Basin, Gage, and Climate Data Filtering" — currently EMPTY. Write from
      the preprocessor docs: 20+ qualifying water years, >30-NA / >3-day-gap year
      rejection, ≤3-day interpolation, seasonal 80 % completeness, 60 % qualifying
      fraction anchored to each window, climate/SWE-valid year tracking.
- [ ] §2.2.2: fill `n=xx` → **6,678 gages (WY 1993–2025) and 6,250 (WY 1980–2025)**;
      name the metadata file(s) once R1–R5 file names are frozen.
- [ ] §3 Data Record: write it around the 5 HydroShare resources — DOI strings
      (predictable pre-publication), per-resource file manifests, formats, dictionary
      pointers, and the code's Zenodo DOI.
- [ ] §4 Data Overview figures: gage map (colored by record length or a headline
      signature), CDF of record years × sites, Whittaker plot (MAT/MAP from HydroATLAS
      or Daymet per basin) — all generable from the staged products; draft in a
      notebook, hand to co-authors for styling.
- [ ] §5 Usage Notes — the limitations everyone was asked to add; seed with the
      documented ones: filter `area_normalized == TRUE` before cross-gage comparison
      of unit-carrying signatures (37 raw-m³/s Canadian gages; flag gap); **never
      compare record-dependent signatures across the two products** (drought
      thresholds, `*_all` pulses, elasticity, parameterized BFI); drought
      `duration_fixed_p10` ↔ pulse-metric redundancy (not independent evidence);
      Daymet 365-day calendar (leap-year Dec 31 hole) + model-product SWE caveats;
      Pettitt multiple-testing guidance (~3.5 % survive BH-FDR); avg_storage omitted
      from major analyses; trend-gate semantics (60 % overall / 80 % decades);
      snow record-anchored gate; how to join companion products (`gage_id`/`canon_id`
      leading-zero rule).

**Queued reconciliation edits (CHANGELOG → Manuscript Reconciliation Log, apply in the
Google Doc — flag each to co-authors via @Arik convention)**
- [ ] 1. §2.1.2 — normalization uses agency-published drainage areas (not shapefile
      areas); add the 37-gage `area_normalized = FALSE` caveat.
- [ ] 2. §2.2.2 — fix metric-family list (de-dupe baseflow; add runoff ratios, Q-P
      seasonality, snow; note storage intentionally omitted).
- [ ] 3. §2.2.2 — complete the statistics list (8 stats + 8 Pettitt fields).
- [ ] 4. §2.2.2 — describe both standard windows explicitly (WY 1993–2025;
      "entire period" = WY 1980–2025; both @ 60 %).
- [ ] 5. §2.1.3 — reconcile 6,041 vs 6,087 Daymet basin count with co-authors.
- [ ] 6. §1 — "HUC-8 scale" → gaged-watershed scale.
- [ ] 7. §2.2.3 — items i/iii/iv flag-vs-remove wording (blocked on the domain-expert
      decision for item i; chase in guidelines doc).
- [ ] 8. §2.2.2 — add the drought family + its two documented deviations from
      Adelsperger et al. (fixed thresholds only; water-year aggregation).

**Draft to-do list items (owners from the doc)**
- [x] Email CUAHSI re: hosting — DONE, response received (this plan responds to it).
- [x] Galen: USGS reviewer ask (Roy Sando?) — marked done in the draft (2026-08-24 sync).
- [ ] Add AI-coding acknowledgement + description of the group's review process
      (check Springer Nature's AI policy for required wording).
- [ ] Kendra: authorship table into the draft.
- [ ] Arik: move code into the CZ GitHub org (then point the Zenodo integration at
      that repo for the release archive).
- [ ] Nic: re-scope the "one nice zip" compilation to the R1–R5 per-resource
      packaging (coordinate before duplicate work happens).
- [ ] Everyone: limitations → §5 (seed list above).

**Submission mechanics**
- [ ] Reviewer data-access package: HydroShare private links per resource (anonymous,
      never expire) listed in the cover letter / data-availability form.
- [ ] Data Records DOIs in the manuscript = the predictable `10.4211/hs.<id>` strings
      (note to editors that they activate at publication).
- [ ] Pre-submission check against Scientific Data requirements: CC-BY, DOI,
      confidential review route, direct download — all satisfied by the plan.

---

## Draft reply to CUAHSI (talking points for the PIs)

- Thank you + we'd like **option 1 (shinyApps server)**: our viewer is already a Shiny
  app; we'll slim it to run against HydroShare-published resources (no external
  credentials, modest memory), which also delivers **option 3** — the app doubles as a
  worked example of building a viewer on published HydroShare resources.
- Sunset expectations (per your caution): host through the promotion phase; annual
  review; minimum ~2 years post-paper-publication; at sunset, resource landing pages
  point to the archived self-contained HTML explorers deposited with the data.
- Questions: shinyApps server specs (RAM/app, bundle size cap, R version, outbound
  network policy, deploy workflow)? Any constraints we should design to now?
- Data side: ~6–8 GB across ~5 resources + a collection under one account — fine under
  standard quota, but flagging since you offered to help with publication costs; also,
  two clarifications: do published resources' storage get released from quota today,
  and does a new version of a published resource mint its own DOI?

## Sequencing (back-planned from Nov 9)

| When | Milestone |
|---|---|
| Week of Aug 25 | Group sign-off on D1–D6; PIs reply to CUAHSI; Google Drive backup inventory + daymet-original check; start Workstream A retrievals; NLCD backup status verified |
| Early Sept | Staging directory complete + hashed; HydroShare resources created (private); resource ids → DOI strings recorded |
| Mid Sept | Per-resource READMEs/dictionaries done; §2.2.1, §3, input table, Usage Notes drafted; reconciliation edits 1–8 applied in the Google Doc |
| **Sept 30** | Full internal manuscript draft (revised internal deadline) |
| October | Co-author review; §4 figures finalized; reviewer data-access package prepared *(app work deferred — removed from this track)* |
| Early Nov | Freeze data package; reviewer private links generated; **submit by Nov 9** |
| At acceptance | Publish resources + collection (DOIs live); Zenodo code release; CHANGELOG + skill updates *(app revisited separately, per the Update block)* |

## Risks / open questions

- **NLCD data-loss exposure — RESOLVED (user, 2026-08-24)**: the finalized NLCD
  product is confirmed in the Drive backup. Its human-QA pass is still owed before
  the R5 deposit. **QA pass COMPLETE — user signed off 2026-08-25; fully closed.**
- **S3 access is gone and the Drive backup is unverified** — until the Workstream A
  inventory is done, the thumbdrive (exFAT, with one prior silent truncation) may be
  the only home of several artifacts; treat the inventory as urgent.
- **HydroATLAS product location unknown** on this Mac — if it's only on the Windows
  box, retrieving it is a manual step; if lost, the re-run needs R.
- **HydroShare publication is immutable** (content/title/authorship) with no correction
  window — everything rides on the freeze check before the publish click.
- **exFAT thumbdrive** has silently truncated a file before — hash at staging time,
  keep a second copy off-drive.
- The two unverified HydroShare policy points (version-DOI behavior; quota release on
  publish) are folded into the CUAHSI reply rather than blocking anything.

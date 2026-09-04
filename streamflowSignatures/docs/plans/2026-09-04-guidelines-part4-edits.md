# Guidelines doc — Part 4 targeted edits (2026-09-04)

Six one-line replacements for the Part 4 findings of the 2026-09-04 accuracy review
(CHANGELOG → Guidelines Document TODOs; Codex-confirmed 8/8). Verified against
`R/helperFunctions.R`, `run_ingest_usgs_hydat.R`, `config.R`, the Julia runner, and the
production metadata.

1. **Entry point** — rename the module `process_gages_rawToRaw()`; Purpose: "Compiles raw
   daily discharge (USGS, HYDAT) into the daily parquet and gage metadata table; computes
   no signatures. process_gages_rawData() is the superseded end-to-end R path."
2. **Minimum-flow gate** — replace the legacy note with: "A gage is compiled only if at
   least 20 calendar years each have at least 30 days with discharge above 0.0001. This
   rule now applies only at compilation; its per-year use inside the signature path was
   removed in April 2026." (Manuscript §2.1.2 glosses the same rule.)
3. **Years excluded** — replace "Missing or invalid years are excluded from analysis."
   with: "All daily records in the retrieval window are kept once a gage qualifies; year
   rejection happens later in the preprocessor (Part 1.2)."
4. **Qualifier mask** — replace "Handles missing or invalid values (e.g., flagged data)."
   with: "USGS values are kept only for qualifiers A, A e, P, P e (others set to NA); HYDAT
   values are never masked, their symbol is carried in the flag column."
5. **Caravan** — replace "Handles redundancy between datasets (e.g., CAMELS and HYSETS)."
   with: "Skips already-processed watersheds on resume; gages present in both CAMELS and
   HYSETS are processed twice. Caravan is not part of the HISSS products."
6. **Daymet join** — add to the Purpose: "Legacy helper, not used by the signature runs,
   which join precipitation and SWE from the Daymet parquet directly and judge climate
   completeness in the preprocessor."

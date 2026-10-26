# Codex final scope-discipline review (2026-10-26, gpt-6-astra, high reasoning, read-only) — "nothing changed but the two facets, and every output change follows from them"

Brief: `scratchpad/codex/review_brief_final_scope.md` of the session (the user's question of 2026-10-26 evening). Verdict: the input patch and the production code change are within scope (independently re-verified in full); the unconditional claim is NO-GO until the products run; GO WITH CHANGES to proceed with S2 and the full-product validation; three reporting corrections (plan §11). Verbatim report:

155,603
**The input patch and tested production behavior stay within the intended scope. The unconditional claim that every output change is physically meaningful is not established.** I found two reporting errors and numerical differences that should not be described as hydrological effects. Read-only review; no files changed.

**Q1 — Input change: PASS on cell-level scope; qualified PASS on interpretation.**

I independently scanned both complete parquets: **111,624,189 rows, 909 row groups, identical schema and row order**. Exactly **193,542 Q cells** changed, all from NA, all matching the ledger’s keys and values. Every non-Q cell is unchanged. Originally nonmissing Q values retain their exact floating-point bits, including 4,470 negative zeros.

I reopened the NWIS cache and independently derived every restored value using the served discharge, February basin area and existing conversion:

`Q = discharge_cfs × 86400 / (basin_area_km² × 3280.84³) × 10⁶`

All 193,542 values match exactly. Site, parameter, daily statistic, units and the unqualified method-block selection checked out. This establishes faithful use of the **October retrieval**, not reconstruction of the values hidden by February’s mask.

Metadata comparisons and independent compile-rule recalculation also pass:

- Only `num_years` (88 rows), `num_days` (88) and `end_year` (2) actually change.
- `start_year` and `error_message` do not change. All 987 targeted error messages were already blank.
- No other cells, rows, ordering or identifiers change.

Other inputs:

- Daymet MD5 matches **`c059076088cd8abd963c64094468e90f`**.
- All **16 reused streamflow/climate batch files** match the hashes recorded by the v2 runs.
- `config/` is unchanged; `config.R` changes comments only.
- Current loaded GAGES-II/HYDAT attributes match both complete v2 products, null-aware. HYDAT also matches HEAD. I could not establish historical byte identity of every unused file in the external GAGES-II directory.
- Both selection manifests exactly match their v2 gage lists: **6,678 / 6,250**.

I independently reproduced all **322,563 baseline/final year statuses**, the Julia patched census, and all **16,028 product-window decisions**, with zero mismatches.

**Q2 — Code change: PASS for the successful production numerical path, subject to S2 verification.**

`git diff HEAD` shows no changes under `julia/src`, `python/streamflow_signatures` or `rpkg/R`.

The three climate helpers preserve the former left join. Outer mode appends climate rows only within the gage’s streamflow span and original water-year set; added Q and flags are missing. Existing streamflow rows retain their joined climate. The unchanged preprocessor still supplies the full calendar and applies the existing eligibility/interpolation rules.

There are additional behavioral changes, but they are explicit safeguards or separate authorized work:

- Date/duplicate checks and the narrowed Julia catch alter failure handling.
- Schema padding can alter subset QA denominators. It is **conditionally**, not inherently, a full-product no-op.
- Batch selection, provenance and validation behavior change.
- The R ingestion rule accepts a broader set of future A/P records under decision 1d; these runners consume the patched parquet and do not execute that ingestion rule.

For batch 1, schema padding added **zero columns**. I independently hashed the v2, left-mode and pinned-schema CSVs **and annual parquets**: each corresponding trio is byte-identical. Full-product no-op behavior remains to be checked across S2.

**Q3 — Output evidence: PASS for the tested climate attribution; INCOMPLETE for universal attribution.**

I independently verified:

- **824-gage batch 1:** byte-identical to v2 at both output layers.
- **40-gage smoke test:** its streamflow and climate subsets exactly match their parent inputs; left annual results exactly match v2. Outer mode changes **599 columns / 4,079 summary cells**, comprising **437 climate-only and 162 coupled columns**, with **no Q-only changes**. Exactly **1,810 annual cells** change.
- Replaying the actual smoke inputs reproduces all raw, cleaned and eligible daily layers for **571,981 dates per mode**, without mismatches.
- **11,800 annual/seasonal precipitation and GDD totals per mode** agree with sums calculated from those daily series.
- **Python comparison:** identical annual keys, with one beyond-tolerance difference: `02160390`, WY2001, `D40_day`, **167 versus 168**. The summary evidence also has the documented propagated/tie differences and `ice_affected_days_total` waiver—Julia **0**, Python **23** at `01011000`.
- **Bounds sensitivity:** the two inputs differ only in Q on the **11,662 censored-bound days**. Qualification is **156 versus 155**; the additional gage is `12202310`, excluded from the frozen product. The reported shared-gage year-count difference is wrong; see finding 2.
- **Five-gage splice:** independently confirmed **1,391 Q-only replacements** throughout the full parquet. This is input-patching evidence, not a signature-output validation.

The climate changes have a sound mechanism: Daymet values replace interpolated climate, then totals, timing, thresholds and record-normalized metrics follow. Changes in distant years through record-wide calculations are expected. The sensitivity census, however, explicitly lacks attribution/recomputation gates; its “unexplained” entries cannot be treated as individually validated.

**Q4 — Recorded exceptions: within the agreed patch scope, with necessary qualifications.**

| Exception | Assessment |
|---|---|
| 97 recoded days | Within scope selected by **February flags**. Today, 96 are `P` and one is plain `A`. Calling every restored value currently approved is false. |
| 554 negative values | Faithfully preserved agency values; consistent with the unchanged configuration. Their sign alone is not a scope violation or proof of invalidity. |
| 31,676 revised February values | Independently reproduced; deliberately unchanged in the patch. |
| 1,991 February-NA days now carrying kept codes | Independently reproduced; unchanged. |
| 8,129 February provisional-remark days accepted today | Independently reproduced; unchanged, consistent with decision 1e. |
| Clearing metadata errors | Authorized by decision 7 and an actual no-op: all 987 were already blank. |

The bounds remain **bounds used as values by policy**, not exact measured discharges. Likewise, restoring Daymet precipitation does not recover an unobserved runoff peak: Q interpolation on absent streamflow dates remains the existing policy.

The findings requiring correction or limiting the claim are:

1. **MAJOR — “Every restored value is approved” is false.**  
   Evidence: [correction ledger](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/input/correction_ledger.csv) and independently reopened cache: **96 `P` values at `12344000`, 2024-05-15–2024-08-18**, retain February `A R` flags. The audit nevertheless calls restored values approved in [splice_nwis_qualifiers.py:899](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/splice_nwis_qualifiers.py:899).  
   **Fix:** describe selection by February approved flags and restoration from October A/P records, explicitly naming the 96 provisional values. Preserve the agreed patch and historical flags.

2. **MINOR — Sensitivity report overcounts changed year counts: 87 should be 79.**  
   [bounds_sensitivity.py:114](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/bounds_sensitivity.py:114) compares all 163 rows. It counts seven `NaN != NaN` cases and the one newly qualifying gage as “shared.” Direct comparison of the 155 shared output rows gives **79**.  
   **Fix:** restrict the count to gages qualifying on both sides; correct `bounds_sensitivity.md` and the plan’s progress log.

3. **MODERATE — Exact numerical changes are overstated as physical changes.**  
   The sensitivity census includes **141 changed annual cells within `rtol=1e-9, atol=1e-12`**. For example, `03289200` drought deficit changes from `31.51028714964012` to `31.510287149640106`. Such differences carry no meaningful hydrological interpretation. Cross-language threshold/rank residuals also require explicit qualification.  
   **Fix:** separate exact differences, numerical residuals and substantive hydrological effects in reporting. This does not warrant an unrelated library change.

4. **RELEASE-BLOCKING EVIDENCE GAP — Universal output attribution remains untested.**  
   S2 does not exist yet; the smoke census is explicitly **INCOMPLETE**, and the sensitivity run omits several attribution gates.  
   **Fix:** complete the full-product checks before asserting universal scope or physical validity.

After the products run, establish full S0–v2 equality for both windows; exact frozen emitted cohorts; unchanged immutable metadata; zero schema padding in every full-product batch; equality for all unaffected gages; oracle-consistent year gains; and complete annual, summary, flag and interaction attribution with explicitly bounded port residuals. Source agreement alone also cannot certify the physical truth of every agency observation or censored bound.

**Final verdict: NO-GO for the unconditional claim as written today.** The verified input scope and tested runner behavior support **GO WITH CHANGES to proceed with S2**, followed by full-product validation and corrected reporting.

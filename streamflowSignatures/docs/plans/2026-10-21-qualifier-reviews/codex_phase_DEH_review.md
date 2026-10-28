# Codex review of Phases D / E / H (2026-10-28, gpt-6-astra, read-only sandbox) — the executed rerun, the post-tag verification tools, the staged deposit

Brief: `review_brief_phase_DEH.md` (same folder). Verdict GO WITH CHANGES, no BLOCKER; dispositions in plan §11 (row "G, round 1 (Codex)"). Verbatim report:

The corrected products are well supported, but the attribution report and staged deposit need changes before release. I found no demonstrated corruption of the product values.

This was read-only: no files edited and no production jobs launched. Checks used the specified Python environment, DuckDB, pandas and PyArrow.

Path shorthand below: **R** = `/Users/arik/HISSS_data/qualifier-rerun-26oct2026`; **D** = `/Users/arik/Downloads/Signatures_v3_28oct2026`; **R1/R2/R3** = its three resource folders.

**Ranked findings**

| Rank | Finding | Required disposition |
|---|---|---|
| **MAJOR** | The claimed full interaction cohort was not executed in S1; the attribution tool silently excludes missing gages. | Correct the coverage gate and regenerate or explicitly narrow the attribution claims. |
| **MAJOR** | Attribution reports conflate absent keys with present-but-NaN values. | Separate key presence from value availability and regenerate counts. |
| **MAJOR** | Staged dictionaries misdescribe the repaired high-NA flag, the unusable ice diagnostic, and identifier joining. | Correct all deposited dictionary copies before deposit. |
| **MAJOR** | Resource 3’s checksum manifest does not match its reconstruction recipe. | Regenerate and verify the manifest after final edits. |
| **MINOR** | Several staged counts, explanations and repository status statements are stale or inaccurate. | Correct the specific items below before committing the release documentation. |
| **NOTE** | The post-tag census and summary exceptions are defensible for these outputs, but PASS means agreement under documented allowances—not unconditional equality or complete causal proof. | Preserve that qualification in public descriptions. |

No **BLOCKER** establishing a need to rebuild the scientific products was found.

**A. Claims and independent numerical checks**

The main numerical claims are supported. My independent checks produced:

| Check | WY 1993–2025 | WY 1980–2025 | Method |
|---|---:|---:|---|
| Product annual rows | 27,766,578 | 34,813,881 | Read parquet metadata. |
| Restored Q-valid years | **1,177** | **950** | Extract old/new `Qann` key sets per gage; compare both sets against `acceptance_oracle.csv`. |
| Lost Q-valid years | **0** | **0** | Same set comparison. |
| Gages disagreeing with either oracle year set | **0 / 6,678** | **0 / 6,250** | Exact year-set comparisons. |
| Untouched gages | **5,069** | **4,894** | Read expected-identical manifests. |
| Untouched annual cells compared | **21,041,343** | **27,034,893** | Chunked full outer joins on `(gage_id, signature, water_year)`, checking missing keys, values, null state and signed zero. |
| Untouched annual differences | **0** | **0** | Same comparison. |
| Untouched summary rows differing | **0** | **0** | Compare complete CSV record bytes against v2; headers also identical. |
| E7 finite values differing by more than **absolute** `1e-9` | **117** | **188** | Independent keyed join of Julia and Python annual parquets. |
| Summary gate-state cells reproduced | **8 / 8** | **7 / 7** | Recompute the relevant statistics from each side’s own annual series. |

These reproduce the central claims in [plan §11](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md:871) and the census evidence at [R, first-window census:36](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/census_1993_2025/change_census_v3_vs_v2_1993-2025.md:36) and [second-window census:35](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/census_1980_2025/change_census_v3_vs_v2_1980-2025.md:35). Annual equality here is keyed value equality, including zero signs; it is not a claim that the parquet files themselves have identical bytes.

The E7 residual breakdown was:

| Residual family | 1993–2025 | 1980–2025 |
|---|---:|---:|
| `D*_day`, including `D25_to_D75` | 65 | 93 |
| `TQmean` | 3 | 4 |
| `dur_low_pulses_all` | 0 | 26 |
| `n_low_pulses_all` | 0 | 11 |
| `elasticity_annual` | 27 | 26 |
| `qp_slope_sd` | 22 | 28 |
| **Total** | **117** | **188** |

Maximum absolute differences were 4.0 / 5.6 overall, approximately `8.25e-8` / `8.39e-8` for elasticity and `6.46e-8` for `qp_slope_sd`. The first window actually has residuals in **17 signatures**; 19 signatures are allowed by name. The second has residuals in all 19. This agrees with the recorded [annual gate](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/e7_1980_2025/e7_annual.txt:17). I verified the counts and family membership; I did not independently reconstruct every threshold crossing from daily data.

Additional confirmations:

- Both S0 summary CSVs have the same SHA-256 as their v2 counterparts.
- The worktree currently reports clean at `54f0483b52fbb1bc743df26f5ee2753a43575be1`; all 16 S2 batch fingerprints record that revision and `git_working_tree_dirty=False`.
- All four E7 reports’ reference/candidate hashes match the current files.
- The annual-change ledgers contain **387,055 / 360,065** records and zero `UNEXPLAINED` labels. This independently checks the ledger totals, not the correctness of every explanation.
- Daily-ledger cause counts sum to **2,702,639 / 2,217,736**. The second includes **18 secondary-interpolation records**, omitted from the three-category breakdown in plan line 874.
- E1/E2/E4 PASS statements have supporting logs; the 16 batch logs contain no `warning`, `error` or `exception` matches. E6 logs support successful basemap checks, although explorers loaded **60/60** tiles and dashboards **24/24**, rather than universally 24/24. See [E6 log:5](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/phase_e6.log:5).

The unsupported or inaccurate claims are principally the interaction coverage, dictionary/staging statements and smaller reporting errors below.

**B. Do the patches weaken gates?**

**(i) Record-dependent key emission/withdrawal — broader explanation, but backed by a meaningful independent gate.**

The new branches accept an old-year key change based on affected-gage membership and dependency classification, without demonstrating the particular emission threshold. See [assessor:288](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:288). Consequently, the explanation alone could conceal an erroneous key omission in a record-dependent family.

For this execution, E7 requires identical Julia/Python key sets with **zero permitted key differences**, and its current hashes and affected-gage coverage check out. That substantially limits the risk, although two implementations sharing a conceptual mistake would still pass.

I independently explained the sole lost key in both windows:

- Gage **05527500**, `elasticity_annual`, WY **2018**.
- Summing eligible daily precipitation gives the 2018-minus-2017 difference as **23.2406334779** in S0, versus **0.0768047762** in S2.
- The implementation requires `abs(dP) > 0.1`; withdrawal is therefore expected. See [elasticity.jl:198](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/elasticity.jl:198).

The 116 first-window old-year exceptions comprise 92 parameterized-BFI gains, 5 annual-elasticity gains, 18 rolling-elasticity gains and that one loss. Thus “0 unexplained” should mean **classified and independently checked**, not that the census itself supplies a mechanistic derivation for every cell.

**(ii) `--allow-gate-state-column` — the executed exceptions pass the stated contract.**

It recomputes Python `generate_stats` twice: once from the **Julia annual parquet**, once from the **Python annual parquet**. Each reconstructed result must match its own summary cell. It also requires membership in a tie-sensitive pair and a matching column regex. See [recomputation:67](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/check_summary_equality.py:67) and [acceptance:217](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/check_summary_equality.py:217).

I reproduced all **8 / 7** exceptions, with zero failures. All observed state exceptions are `FDC90th` fields. I found no wrongly accepted cell among them.

This establishes propagation consistency, not agreement between the two annual series or correctness of their physical calculation. Those depend on the separate annual gate. Also, ordinary finite-value column waivers remain broader: [lines 238–240](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/check_summary_equality.py:238) bypass their numeric mismatch gate.

**(iii) Qann fallback — it does not replace the oracle comparison.**

The fallback changes how non-replayed gages’ valid-year sets are obtained. Gate 5f still compares their baseline and final sets against the oracle and checks summary year metadata. See [assessor:373](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:373).

My all-gage comparison found zero year-set disagreements. This is acceptable for these untouched gages; it is not equivalent to having replayed their preprocessing state.

**(iv) Scoped `total_possible` freeze — a narrower check, adequately supported here.**

[Lines 394–402](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/assess_rerun_changes.py:394) permit absent replay evidence only for class U. This is weaker than demanding replay tables for every gage, and should be described as such.

Independently, baseline and final oracle denominators agree for **all 6,678 / 6,250 selected gages**, and untouched outputs agree. I found no resulting denominator defect.

The climate-invalid-NA patch also matches the stated distinction between raw/clean and eligible values: it keeps an invalid variable’s missing day NaN while still checking finite raw values and eligible masking. See [daily ledger checker:303](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/preprocess_state_ledger.py:303).

**C. Stage attribution**

The executed identity tests are meaningful. I reran the annual comparisons:

| Identity | 1993–2025 shared rows | 1980–2025 shared rows | Differences |
|---|---:|---:|---:|
| S1 = S0, 200 controls | 833,690 | 1,101,466 | 0 |
| S1′ = S2, 197 / 199 join-only gages | 813,718 | 1,136,764 | 0 |

Excluding 3 / 1 sampled join-only gages because they belong to the interaction cohort is legitimate. The reported summary identities are consistent with the corresponding artifacts: [first window:5](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/attribution_1993_2025/attribution_1993-2025.md:5), [second window:5](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/attribution_1980_2025/attribution_1980-2025.md:5).

**MAJOR: the full interaction-execution claim is false.** Comparing cohort JSON membership with both manifests and actual S1 CSVs shows:

- First window: **7 of 178 absent from S1** — `0214269560`, `02294491`, `02484000`, `03380500`, `06187915`, `07285500`, `14301500`.
- Second window: **5 of 149 absent** — `02294491`, `02484000`, `03380500`, `07285500`, `14301500`.

The cohort builder constructs S1 from `AQ + control`, while adding eligibility and December-31 cases only to S1′: [builder:75](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/build_stage_cohorts.py:75). The attribution tool checks each manifest’s execution but never requires the interaction union to be contained in S1; its summary calculation silently intersects available gages: [attribution:96](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/stage_attribution.py:96), [line 128](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/stage_attribution.py:128). Hence the actual four-stage summary cohorts are **171 / 144**.

**MAJOR: “keys absent” is incorrectly calculated from `.notna()`.** [Lines 115–118](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/stage_attribution.py:115) treat stored NaNs as absent keys. Independent full joins with explicit presence markers give:

| Measure | 1993–2025 | 1980–2025 |
|---|---:|---:|
| Reported keys absent from some stage | 88,004 | 86,782 |
| Actual keys absent from some stage | **63,188** | **57,991** |
| Keys present in all four, NaN in all four | **24,747** | **28,734** |

Use separate key-presence and finite-value masks. Either execute the missing S1 cases or identify their S1=S0 substitution explicitly as inference; do not describe them as executed.

The observed finite-cell interaction residual is nevertheless plausible. Recomputing `(S2−S1)−(S1′−S0)` reproduces **2,650 / 2,779** cells outside tolerance, maxima **21.2322 / 20.1614**, including **2,329 / 2,533** annual-elasticity and **292 / 227** rolling-elasticity cells.

The plan’s explanation that elasticity “is a regression across years” is inaccurate. It uses ratios of differences/anomalies, record or window means, medians and precipitation thresholds: [elasticity.jl:115](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/elasticity.jl:115). Those operations are nonlinear, so non-additivity is expected. The residual itself is not evidence of an additional product defect.

**D. Staged documentation, files and reconstruction**

**MAJOR — dictionary semantics must be fixed.**

- [R1 dictionary:2272](/Users/arik/Downloads/Signatures_v3_28oct2026/resource1_signatures_wy1993_2025/hisss_data_dictionary.csv:2272), also copied to R2 and the collection root, still says `flagged_for_high_na` is broken and must not be used. I recomputed `NA fraction > 0.30` across the **2,245 signature-derived columns**: **zero flag disagreements**, with **471 / 345 TRUE**. The README describes the repaired definition correctly.
- [Dictionary:2254](/Users/arik/Downloads/Signatures_v3_28oct2026/resource1_signatures_wy1993_2025/hisss_data_dictionary.csv:2254) describes `ice_affected_days_total` as an actual diagnostic count. Both staged products contain **only zero**. The known structural cause is documented in [Julia io.jl:419](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/io.jl:419). Users must be told zero does not establish absence of ice-affected days.
- [R3 dictionary:19](/Users/arik/Downloads/Signatures_v3_28oct2026/resource3_input_data/hisss_input_data_dictionary.csv:19) recommends padding numeric IDs to eight digits. Independently joining identifiers yields **66 unmatched compiled gages with that advice**, versus **zero** with normalization by stripping leading zeros. The README correctly prohibits padding.
- [R3 dictionary:74](/Users/arik/Downloads/Signatures_v3_28oct2026/resource3_input_data/hisss_input_data_dictionary.csv:74) describes metadata changes as `gage_id,column,old,new`; the actual file is wide, with `old_num_years`, `new_num_years`, etc. Grouped and omitted definitions also do not support the claim that every new column is individually defined.

**MAJOR — stale reconstruction checksum.**

[R3 staging manifest:26](/Users/arik/Downloads/Signatures_v3_28oct2026/resource3_input_data/hisss_staging_manifest_inputs.csv:26) records:

- **2,683 bytes**, SHA-256 `00c6529aa100db03203be592895345194816d99ce61c5c2bb01bfe9bf851282c`.

The current recipe is:

- **2,858 bytes**, SHA-256 `cf50508eff761d3bc28efd59160f8c9d393ca25481fc8de5be0e5eb9947c8c60`.

All other checked top-level manifest entries matched: **20 files each in R1/R2, and 24 other files in R3**. I excluded the 5.1-GB Daymet file and explorer directories from checksum recomputation. Directory entries are counts/sizes rather than member hashes.

**MINOR — concrete text corrections.**

| Location | Correction |
|---|---|
| Both run notes, [line 69](/Users/arik/Downloads/Signatures_v3_28oct2026/resource1_signatures_wy1993_2025/hisss_run_notes_wy1993_2025.md:69) | Replace `? values` with verified **117 / 188**. The generator requests JSON fields the annual checker does not emit. |
| Both run notes, [line 19](/Users/arik/Downloads/Signatures_v3_28oct2026/resource1_signatures_wy1993_2025/hisss_run_notes_wy1993_2025.md:19) | “Restored days in valid years” occurs on **891 / 749 gages**, not 899 / 757. The latter are restored-day gages anywhere inside the window. |
| [R3 README:84](/Users/arik/Downloads/Signatures_v3_28oct2026/resource3_input_data/hisss_readme_inputs.md:84) | Corrected-table missing-Q count is **394,210**, not approximately 588,000. The latter describes the February table. |
| [R3 README:115](/Users/arik/Downloads/Signatures_v3_28oct2026/resource3_input_data/hisss_readme_inputs.md:115) | The one compiled gage without forcing, **05KH009**, is in **neither product**, contrary to the sentence saying it carries signatures there. |
| [R3 abstract:1](/Users/arik/Downloads/Signatures_v3_28oct2026/resource3_input_data/hisss_abstract_resource3.txt:1) | “Every compiled basin” overstates Daymet coverage: it covers **8,013 of 8,014**. |
| [Plan:876](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/plans/2026-10-21-qualifier-mask-fix-rerun-plan.md:876), [cross-language status:55](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/CROSS_LANGUAGE_STATUS.md:55) | The **103 finite-value column-waiver cells** and **8 / 7 gate-state cells** are separate counts; “of which” is incorrect. |
| [STATUS:21](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/STATUS.md:21) | Distinguish the earlier 29 / 26 forecast from the measured **31 / 29** would-qualify lists; update completed E6/H work currently listed as NEXT. |

The public name **“version 2”** is clearly reconciled with internal build v3 in run-note lines 1 and 8. That naming is locally consistent; actual publication history was not externally verified.

The two streamflow tables are clearly distinguished. Both have **111,624,189 rows, 909 row groups and the same Arrow schema**. The ledger has **193,542 unique gage/date keys on 987 gages**. Its current statuses are **193,446 A and 96 P**; all 96 P records are gage `12344000`, 2024-05-15 through 2024-08-18. The README correctly says the stored flag is the February flag, but the deposited audit’s “agency’s approved values” wording should explicitly acknowledge those 96 currently provisional values.

The would-qualify lists are correct: **31 / 29**, all outside the frozen cohorts and all failing baseline qualification while passing final qualification.

Other checked input counts agree: **16,994** metadata rows; **134,605,430** Daymet rows over **8,017** basins; **101** basin flags, including 57 low-confidence outlines; **8,012** HYDAT rows; control comparison **13,755,058 / 13,786,734** equal, leaving **31,676** revisions.

**Reconstruction contract:** the current recipe matches the repository source, uses round-trip float parsing, applies ledger values only over missing baseline Q, and compares other columns. See [recipe:18](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/reconstruct_qualifier_patched_streamflow.py:18). My independent checks of row groups 0, 200, 500 and 908 compared **406,171 rows**, applied **1,353 corrections**, and found no Q, other-column or signed-zero differences. The existing [full-run log:1](/Users/arik/HISSS_data/qualifier-rerun-26oct2026/recipe_check.log:1) reports all 111,624,189 rows passing; I did not rerun that full reconstruction.

As a general validator, the recipe could be hardened: duplicate ledger keys collapse into a dictionary, schema equality is not explicitly asserted, and Q equality is numeric rather than bitwise. Those limitations did not produce a failure in the staged-file checks.

**E. Before commit/deposit versus later**

Before committing release claims or depositing:

1. Fix interaction coverage accounting and distinguish absent keys from NaNs; regenerate attribution reports and corresponding plan/CHANGELOG statements.
2. Correct all dictionary copies, especially high-NA, ice-count and identifier semantics.
3. Correct the specific staged/documentation inaccuracies above.
4. Regenerate checksum manifests **after** the final recipe and documentation edits, then verify them.

Can wait without rebuilding these products: generalized recipe hardening, narrower cell-specific waiver inventories, and future fixes to the already documented scientific/port residuals.

**Unverified:** a full repeat of source-cache restoration, every daily attribution mechanism, every scientific dictionary definition, historical Daymet equivalence to the co-authors’ aggregation, the complete Daymet checksum, every explorer sidecar/browser interaction, and external deposit/publication state. Recorded PASS artifacts support several of these, but this review does not independently certify them.

**Verdict: GO WITH CHANGES.**
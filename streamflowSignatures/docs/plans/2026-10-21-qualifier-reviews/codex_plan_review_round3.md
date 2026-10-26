# Codex plan review, round 3 — plan version 4 (the climate-frame fix folded in; 2026-10-26; gpt-6-astra, reasoning high, read-only)

Brief: §10 (fix design, footprint forecast, three-stage attribution, change census, decisions 12–17, §10.8 coupled Q–P metrics) and its interactions with §0–§9. Verdict and findings as returned:

**NO-GO as written.** The span-only join changes the preprocessing year set and Q-only diagnostics; several proposed acceptance gates reject correct results. Before implementation, preserve the original water-year support, expand the join ledger beyond inserted dates, replace the universal delta/direction gates, and specify consistent attribution trees and dependency rules.

1. **BLOCKER — §10.2: span restriction does not preserve the year loop.**

   [`io.jl:256`](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/io.jl:256) iterates the years present in the input. An outer join inserts entirely absent interior years. Those years fail Q qualification, but their seasonal flags survive; [`signatures.jl:62`](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/signatures.jl:62) counts their exclusions.

   Independent census analysis found **3,562 / 5,116 absent interior gage-years**, across **759 / 894 product gages**. Of those gages, **560 / 677 currently belong to forecast class U**. This is consequential: `01029500` has no streamflow rows in WY1993, but Daymet supplies 364 rows.

   A Julia fixture containing WY2000 and WY2002 confirmed: the span-restricted outer join leaves Q-valid years unchanged but adds rejected WY2001 and changes winter exclusions **0→1**. Thus the forecast and “Q-only columns identical” gate fail.

   **Change:** restrict added climate rows to both the original date span **and the original set of water years**. Test valid years, rejected years, seasonal flags, diagnostics, and qualification denominators.

   The unrestricted-denominator concern is real, but overstated: the census counterfactual gives **0 / 325 product drops**, respectively. Span restriction preserves denominators for these actual products; it is not a general proof. In the same fixture, a window ending in WY2001 changes `total_possible` **1→2** despite unchanged global endpoints.

2. **HIGH — §10.2: the two Julia joins have different contracts.**

   [`run_julia_benchmark.jl:320`](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark.jl:320) places the second join inside the **legacy branch**, after filtering qualifying years. Production instead consumes the preprocessing cache. A span-only helper at the legacy join reintroduces rejected interior years.

   For production, the column claim is substantially correct: `io.jl:297` selects date, Q and climate variables; the grid constructs month/dowy. I found no production-family dependency on input `year` or `doy`. Input `flag` is dropped before the existing ice diagnostic—do not silently repair that unrelated behavior.

   Legacy functions do read input month/dowy, including runoff ratios, seasonal pulses, timing, baseflow and storage. Calling `add_water_year_columns` after joining supplies them; assigning only `water_year` does not.

   **Change:** explicitly distinguish production and legacy paths. For legacy, retain the already-selected year set and regenerate all temporal columns. Narrow the Julia join’s broad `catch` at lines 218–223 so helper failures cannot masquerade as unavailable climate. A runner helper remains my choice for this release; a separate climate argument is a reasonable later library contract.

3. **HIGH — §§10.3/10.5: the ledger and annual-delta gate are mathematically incomplete.**

   Changes extend beyond added rows. Independent replay found **three existing-date PPT changes in each window**, all on Daymet’s missing December 31. For example, `07285500`, 2016-12-31 changes **0.3328632184→2.7973501605 mm** because an interpolation neighbour is restored. TMIN/TMAX also change there.

   Moreover, [`climate.jl:304`](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/climate.jl:304) distinguishes totals from nonlinear metrics. Gini, seasonality index, timing/centre-of-mass, dry spells, maxima—including the five-day maximum—SD/CV and threshold counts cannot equal the sum of raw daily precipitation differences. Temperature and snow bases require their own formulas.

   **Change:** ledger **all changed preprocessed daily values**, including secondary interpolation changes, with old/new eligibility and missingness states. For shared valid years, annual/seasonal precipitation totals must equal the appropriately scoped sum of daily deltas. Other bases must equal recomputation from the reconstructed daily series. Use absolute-plus-relative tolerances; relative tolerance alone is unsuitable near zero. Treat newly eligible years as gains, not numerical differences from an invented old interpolate.

4. **HIGH — §10.8: retain decision 17’s default, but remove its direction gate and accuracy guarantee.**

   In the production path, coupled calculations receive Q-valid, fully interpolated years; PPT-dependent families receive PPT-valid years as well ([runner:348](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark.jl:348)). Thus the concern about an unfilled partial-year Q numerator is addressed **for production**, not for unrestricted library/legacy calls.

   The family-specific qualifications still matter:

   - Runoff ratios use paired finite days, precipitation floors and raw-Q seasonal completeness flags (`runoff_ratios.jl:83–136`). These reduce eligibility; production supplies complete paired series.
   - Elasticity filters low-P years. Its “annual” calculation compares adjacent **retained observations**, potentially across missing calendar years, and divides by record-wide `Q_mean/P_mean` (`elasticity.jl:205`). Even distant annual values can change.
   - Q–P seasonality uses rolling-window precipitation thresholds and minimum monthly coverage, without the runoff-ratio seasonal mask.
   - Storage requires sufficient years and at least ten unique Q values. Zero/intermittent-flow cases can produce unchanged zeros, missing elasticity, or no storage value.
   - `swe_max_to_ppt` requires SWE/PPT eligibility; its SWE numerator also changes. Snow’s record-anchored gate affects trend eligibility when SWE-valid years change (`snow.jl:264–321`).

   The October artefact contains **501 precipitation decreases and 418 increases** under correction. `02168504`, WY2016, decreases by **167.941288 mm**. Missing rows do not establish a universal storm-exclusion direction.

   Even positive annual ΔP does not imply increased storage. A Julia experiment using the actual storage function, with ΔP=−10 on day 10 and +20 on day 350, produced **annual ΔP=+10 but Δstorage=−10**. Storage interpolates cumulative balances at mean discharge; it is not the year-end balance ([storage.jl:85](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/storage.jl:85)).

   **Change:** gate against formulas, not directions. Ratios decrease only with fixed nonnegative numerator and increased corresponding denominator; SWE ratios and seasonality have no general monotonicity. Storage is monotone under pointwise nonnegative P changes with fixed Q, not merely positive annual ΔP.

   Keep full-frame P. Interpolating Q can nevertheless miss most runoff in a short flashy event; “≤30 days” does not bound volume error or prove the residual is small. Artificially interpolated P can accidentally offset Q error. That is not a justification for decision 17’s alternative.

5. **HIGH — §10.4: tree selection and the “interaction check” are underspecified.**

   D-S1 explicitly uses the **Phase-A commit**, before A′, while selection manifests, fingerprinting, validators and other harness changes arrive in C. D0 merely says “the committed tree”; the implementation order puts A′ before D0. Consequently, neither the intended old join nor otherwise identical execution code is assured.

   **Change:** use one clean, identified tree containing the necessary tooling, with an explicit fingerprinted old/new join mode for all stages; alternatively, identify clean trees differing only in the join implementation. Record matching configuration, cohort, metadata and batch identities.

   S0→S1→S2 gives **ordered attribution**. It does not uniquely separate biological/statistical interactions. Checking only non-spliced gages cannot measure interaction. On those gages, test endpoint equality **S0=S1 and S1′=S2** directly.

   For overlapping causes, report `(S2−S1)−(S1′−S0)` on a deliberately selected batch containing overlaps, eligibility gains and interpolation-neighbour cases. Do not require zero. Handle gained/lost cells separately and bound floating-point residuals using endpoint scales.

   A cheaper rigorous S1 is possible: rerun complete records only for gages whose Q changed, reuse verified S0 results elsewhere, and record that provenance. Arbitrarily restricting S1 to batch 1 is not equivalent.

6. **HIGH — §10.5/E8 and §5: “0 unexplained changes” needs an output-level dependency contract.**

   A family-level gage-year label is insufficient. Record dependence includes static/annual elasticity’s normalization, rolling windows, family eligibility thresholds, summary statistics and trend gates. Adding years can enable the stats floor or disable trend completeness by extending the record span. Snow additionally changes its anchor-year denominator.

   [`stats.jl:395`](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/stats.jl:395) applies the statistics floor after annual collection. [`qa_qc.jl:183`](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/qa_qc.jl:183) computes `flagged_for_high_na` from the final signature fields. Neither is adequately represented by daily ledger keys alone.

   **Change:** define dependency closure separately for annual values, summaries, scalar diagnostics, eligibility and QA flags. Persist per-stage Q/PPT/SWE/temperature year sets and gate reasons. “Covered by the map” should mean **permitted to change**, not numerically validated; affected values still need replay/parity checks.

   Resolve the remaining contradictory instructions: §0 and E5 still compare splice-only class U and exact climate equality against v2. Apply E5(iv) only to **shared-year climate-only annual values across S0→S1**; exclude coupled `swe_max_to_ppt`. E8’s “climate-only columns identical” cannot include summaries when years are added. Generate final U/A membership from realized effects, not fixed forecast counts.

7. **MEDIUM — §10.3: totals reproduce, but several labels and reconciliations are wrong.**

   Using product gage sets and `365 + isleap(water_year) − n_rows` for census-valid years, I reproduced every forecast row and these totals:

   | Quantity | 1993–2025 | 1980–2025 |
   |---|---:|---:|
   | Q-valid years | 200,834 | 254,195 |
   | Years with absent dates | 1,174 | 938 |
   | Gages / absent dates | 857 / 3,186 | 722 / 2,198 |
   | Splice overlap / union | 147 / 1,609 | 122 / 1,357 |

   The October difference is **eligibility**, not counting bookkeeping: `0214269560/WY1997` has 30 absent streamflow dates plus Daymet’s missing December 31, exceeding the climate NA limit; `02484000/WY2017` has a four-day climate gap. Replay restores PPT, SWE and temperature eligibility for both. Product #2 gains one year in each climate family.

   Further corrections:

   - The long-window decade counts are **14/47/59/419/399 for 1980s–2020s**; the table omits the 2000s and mislabels earlier counts.
   - Three absent dates per product are themselves Daymet-missing December 31 dates. Baseline inserted climate-row counts are therefore **3,183 / 2,195**, not the absent-date totals.
   - Forecast splice-restored years contain another **26 years/59 absent dates** and **15 years/26 dates**, respectively. Include these in realized attribution.
   - No Canadian gage occurs in the forecast.

   “Reproduce the Daymet sum” must include the documented calendar interpolation and specify the year set. Historical gage means such as 892→1,182 are not fixed acceptance targets after splice-restored years enter the mean.

8. **MEDIUM — census implementation, paper-facing evidence and release integration.**

   The census is feasible below 3 GB using keyed gage partitions and column-wise accumulation; it does not require materializing annual-table outer merges. However, the named ledgers alone cannot explain eligibility/gate transitions without the additional records in finding 6.

   Specify relative-change handling at zero, identity-R² versus correlation², constant-series handling and cross-unit aggregation for “largest aggregate change.” Inventory p-values from the schema: these products contain **695 `_pval` columns across five suffix types**, not eight fields.

   A data-paper assessment should distinguish same-year value corrections from changed year composition, report regional/seasonal distributions, and quantify interpolated-Q exposure alongside rainfall on those days. Product-wide Q-only means may change from the splice; §10.5’s expectation otherwise is unsupported.

   Extend E7 with long-window-only cases and the new validity/interpolation fixtures. Add the join ledger, stage manifests and gate-state census explicitly to Phase H’s deposit inventory. Document climate support as **Q-valid water years**, not an unrestricted full climate record.

   Finally, §0 still says seven person-days/20 hours; §8 says ten/~28. Update the headline, enumerate download and S1′ time, and include rerun contingency before retaining the November 9 claim. Decisions 12–17 are also still described elsewhere as 12–16.

**Verified independently**

Read-only runner/library inspection; forecast reaggregation and exact row comparison; October artefact reconciliation; 65,536-row batch scans and climate-validity/interpolation replay for forecast years; interior-year census; Julia year-loop, denominator and storage counterexamples. No internet or file changes.

**Taken on trust**

The supplied census’s full-domain Q qualification labels, source-retrieval assumptions, prior accepted splice findings, and runtime estimates. No production stages or full cross-language reruns were executed.

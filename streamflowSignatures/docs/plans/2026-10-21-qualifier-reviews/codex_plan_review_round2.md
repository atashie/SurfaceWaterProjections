# Codex plan review, round 2 — delta review of plan version 2 (2026-10-23; gpt-6-astra, reasoning high, read-only)

Brief: RESOLVED / PARTLY / NOT RESOLVED per round-1 finding, new findings on version 2 (freeze mechanism, decisions 1–2, year-level oracle, §5 criteria, Phase C, effort). Verdict and findings as returned:

**GO WITH CHANGES.** The frozen cohort and patched-snapshot approach are defensible. The acceptance contract, reproducibility mechanism, and Phase C scope need the changes below before production.

The round-1 findings stand as follows; “resolved” here means resolved in the proposed design, not implemented.

1. **PARTLY RESOLVED — Oracle.** All **1,656 CSV rows** match the supplied replay and masked-day file; an independent parquet replay matched **259 sampled gage-years**. The forecast still assumes targeted values become finite, and the CSV lacks the final valid-year fields demanded by B7.
2. **PARTLY RESOLVED — Provisional values.** The false premise is withdrawn. However, fetching only 987 approved-code gages cannot support the promised audit of all provisional-remark days.
3. **RESOLVED — Candidates and metadata keys.** The freeze removes the admission question; `11129800` is acknowledged, and metadata updates now require unique normalized keys and in-place changes.
4. **RESOLVED — Footprint.** Reaggregation confirms **193,542 approved masked days across 987 gages** and the per-code totals.
5. **RESOLVED — Climate attribution.** The join explanation is corrected. Code inspection and a Julia experiment support exact shared-year climate annual equality, including a shared year with restored Q.
6. **PARTLY RESOLVED — Qualifier policy.** Storage and analytical decisions are explicit, but the proposed “sensitivity table” measures exposure rather than numerical sensitivity; future unknown remarks remain blanket-accepted.
7. **PARTLY RESOLVED — Vintage attribution.** The patched-snapshot label, ledger, raw responses, and alternative ingestion strategy address provenance. The proposed control cannot establish historical values on originally masked dates.
8. **PARTLY RESOLVED — Input gates.** Exact cached-record verification replaces plausibility as the correctness gate. B2 still invokes unavailable historical raw values, and B8 cannot be established from a flag/NA tally.
9. **PARTLY RESOLVED — A/U comparison.** The partition is exhaustive and disjoint when based on actual changed Q keys. Blanket metadata equality conflicts with restored years, and several affected-result checks remain diagnostics.
10. **PARTLY RESOLVED — Harness/provenance.** Phase C identifies substantial work, but its resume fingerprint, baseline control, qualification auditor, and bounded-memory validation remain insufficiently specified.
11. **RESOLVED — Equality semantics.** Keyed equality and NaN/null/signed-zero semantics are now explicit.
12. **PARTLY RESOLVED — Staging.** A qualifier-release mode is specified. The deposit inventory still needs the frozen membership manifests, executable selection recipe, and explicit immutable-base inclusion.
13. **PARTLY RESOLVED — Documentation.** Historical preservation is addressed, but current rule documentation remains omitted—for example, `docs/SIGNATURE_GUIDELINES.md:540`.
14. **PARTLY RESOLVED — Effort.** Tooling and review are budgeted, but the table totals **6.2–6.7 hands-on person-days**, inconsistent with the headline without explicit staffing and overlap.

New findings, most severe first:

1. **BLOCKER — E5(ii) and Phase C’s “metadata equality” reject the intended result.**

   **Evidence:** The [runner](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/run_julia_benchmark.jl:404) derives `num_water_years`, `start_water_year`, and `end_water_year` from qualifying years and explicitly classifies them as metadata. The four-code replay increases year counts for **619 / 502 existing product gages**, so blanket equality necessarily fails.

   **Change:** Enumerate immutable site attributes separately. Gate year count and endpoints against the final oracle. Define affected membership using **actual Q changes**, excluding audit-only or unchanged ledger records. R²/max-difference reports are diagnostics: require a documented review disposition and an independent numerical check of affected results, rather than treating report generation as acceptance.

2. **MAJOR — The freeze needs a deposited selection contract and a qualification-aware auditor.**

   **Evidence:** Filtering is computationally sound, but the full corrected R3 input plus existing wrappers would admit additional gages. “Filtered input, recorded in provenance” does not specify how a reader reconstructs that input. Moreover, [audit_qualification.jl](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/docs/benchmarks/audit_qualification.jl:69) equates qualification with inclusion: excluded candidates either fail against the full input or are unavailable against the filtered input.

   **Change:** Deposit canonical, unique per-product membership CSVs, their hashes, the originating v2 reference hashes, and executable commands. Prefer a first-class `--gage-list` selection in the splitter/runner against the full patched input. Make the auditor distinguish **qualified**, **selected**, and **emitted** sets; require emitted = selected and every selected gage qualified. Describe the products as a frozen cohort, not all qualifying gages.

   Filtering changes sorted-rank batch assignment for **5,927 / 5,610 retained gages**. It does not inherently change their calculations or normalized metadata joins. `flagged_for_high_na` depends on each batch’s available columns, so retain exact schema/denominator checks; all existing batches have identical **2,277-column** headers.

3. **MAJOR — The oracle is a verified counterfactual, not yet an exact source-backed gate.**

   **Evidence:** The replay matches all **39,743 supplied baseline statuses**. I reproduced §1’s **1,197/601** and **1,370/676** year/gage totals; product gains **1,016/550**, **810/444**, **1,177/619**, **950/502**; and **500/403** already-valid touched years. All §6 lists match `fullcal_newly.csv`: **29/26**, or **31/29** with four codes.

   Example oracle rows also reconcile: `01096500` has 22 `A <` days restoring WY2012; `01109070` has 72 `A R` days restoring WY2019–2020; `06400875` restores nine years with four codes versus one with `A R`.

   These results assume restoration availability. The CSV contains neither `n_valid_new` nor complete baseline/final year lists, although B7 demands them.

   **Change:** Preserve this forecast separately. Generate the acceptance oracle independently from **base + eligible cached source records**, then compare the spliced parquet against it. Include baseline/final year sets and counts, source-unavailable targets, selected policy, and configuration hashes. Recompute would-qualify lists from realized restorations instead of requiring an unconditional forecast.

   I found no additional Q-admission discrepancy in the sampled replay under current settings. It does not model cleaned values, seasonal completeness flags, climate qualification, or diagnostics; these need separate checks.

4. **MAJOR — Fetch coverage, B2, B8, and the control do not share a computable contract.**

   **Evidence:** The masked-day file contains **198 provisional-only gages with 12,277 provisional-remark days** outside the 987-gage fetch. Historical raw values were not retained, so B2 cannot identify values “identical in raw units” independently. A `flag × is.na(Q)` tally cannot establish current source availability.

   **Change:** Separate the **audit-fetch universe** from the **splice universe**; fetch all historically masked gages needed by B3 even when 1e excludes their changes. Verify conversion against cached records and metadata, reporting differences in normalized units without asserting unavailable historical raw identity. Compute B8 through a daily source/ledger/output join.

   Define the control as three datasets: immutable base, October data under the old rule, and October data under the new rule. The latter pair isolates rule effects **at October vintage**. It cannot isolate historical vintage effects on February-masked dates; state that limit.

5. **MAJOR — Decisions 1b–1e still blur policy and evidence.**

   **Evidence:** Counts, low-flow exposure, and a list of dependent signatures do not measure the effect of substituting censoring bounds. Also, February availability is unknown for approved masked days as well as provisional ones; that uncertainty does not establish that one restoration is a bug fix and the other necessarily a refresh.

   **Change:** Describe 1e as an explicit scope exception. Report numerical qualification/signature sensitivity for bound substitution versus the named alternative. Scope 1d explicitly to ingestion versus this historical patch, and require review of newly encountered remark categories before analytical inclusion. Keeping unknown finite values in raw storage is a separate decision.

6. **MAJOR — Phase C still permits stale results and overstates D0.**

   **Evidence:** A wrapper’s git blob can remain unchanged while included Julia code or configuration changes. Original-input batch 1 demonstrates reproducibility only for its sampled gages; class U cannot establish the missing baseline source for affected gages elsewhere.

   **Change:** Bind resume to the complete source revision, runtime/project environment, configuration, metadata, selection manifests, and actual batch-input hashes; verify output identity too. Production should reject dirty or unknown provenance. Expand the original-input control to cover affected gages and both windows, or explicitly document its limited coverage. An untracked-file list does not preserve file contents.

7. **MAJOR — Validation and staging adaptations remain incomplete.**

   **Evidence:** Annual validators still materialize entire tables; the existing products contain **27,602,161 / 34,681,331 annual rows**, and cross-language comparison loads both sides. Parameterizing only the rpkg reference and column count leaves its candidate directory hard-coded. Phase H’s inventory omits the immutable base promised by decision 2.

   **Change:** Budget and implement partitioned validation/comparison, explicit production floor arguments, and rpkg candidate/reference/log paths with named residual allowances. Stage and verify every required artifact: immutable base, patched input, raw responses, ledger, controls, final oracle, selection manifests, and executable reconstruction recipe.

8. **MINOR — Effort and documentation need another reconciliation.**

   **Evidence:** The table requires **6.2–6.7 person-days** before additional reruns; B/C overlap requires separate capacity. `docs/SIGNATURE_GUIDELINES.md:540` still states the old mask, while the current CHANGELOG issue retains the superseded 1,214-year estimate.

   **Change:** Distinguish person-days from elapsed days, name parallel owners, and budget the expanded audit/control work. Add omitted current documentation to the checklist while preserving dated history.

E5(iv) itself should remain a gate. The runner constructs precipitation, temperature, and SWE frames from their respective valid-year sets; adding years does not change shared-year climate inputs. A canonical Julia experiment restoring Q in one shared year and admitting another produced **53 shared annual values, zero differences**, and 53 new values. Summary completeness gates operate separately. Use an explicit family registry—including snow metrics without `swe_` prefixes—and avoid claiming every inequality uniquely diagnoses a moved join.

**Verified independently:** oracle reconciliation; §1 totals and §6 lists; 259 sampled parquet/calendar replays; baseline-status agreement; provisional audit omissions; batch reassignment and headers; annual row counts; relevant runner, auditor, validator, staging, and provenance code; climate experiment; effort arithmetic.

**Taken on trust:** exhaustive cached-response claims from round 1, February retrieval provenance, full-copy byte identity, future service availability, prospective restoration availability, and production/review runtimes. No internet or file changes were used.

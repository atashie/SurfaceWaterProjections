[Full round-3 report](/tmp/codex_julia_review_round3.md). No repository files were changed.

**Round-2 reproductions**

| Finding | Result | Actual round-3 values |
|---|---|---|
| 1. Halfway crossing | **PASS** | Seven-day rainfall blocks: **day 7**. Original 50-day blocks: **day 50**. Uniform 0.1 mm in leap WY 2000: **day 183**. Totals respectively **1.4000000000000001**, **9.99999999999998**, **36.60000000000025**. |
| 2. PPT non-finite preprocessing | **PASS** | Isolated `Inf`, `−Inf`, and `NaN` each interpolate to **2 mm**, retain WY 2001, and produce annual precipitation **730 mm**. Four-day runs and boundary infinities give **`valid_climate_years=[]`**. |
| 3. Boolean configuration loopholes | **PASS** | Boolean dry threshold, frost-end month, and temperature allowances are rejected. Numeric **0/1** flags are rejected. Threshold list **`[true,25,50]`** is rejected. |
| 4. Fraction canonicalisation | **PASS** | `0.500000000001` becomes **0.5**. Endpoint storms return **`prcp_50_dowy=1`**. Original `0.504` remains rejected. |
| 5. Finite interpolation neighbours | **PASS** | `[0,NaN,Inf]` remains **`[0,NaN,Inf]`**. Original temperature preprocessing fixture yields **`[0,1.6666666666666667,3.3333333333333335,5]`**, a qualifying year, and GDD **725.5**. |

1. **Severity: MINOR — canonicalisation can admit a fraction outside its validated range.**

   **Location:** [config.jl:210](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/config.jl:210), canonicalisation at [config.jl:222](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/config.jl:222).

   **Contract:** Computation uses the canonical whole percentage. Configuration validation also requires the fraction strictly inside `(0,1)`, implying supported percentages **1–99**.

   **Behaviour:** The range check applies before canonicalisation. Near-endpoint inputs pass and subsequently become **0.0 or 1.0**, without another range check.

   **Minimal reproduction:** Complete WY 2001; P[365]=2 mm, zero elsewhere. Change only `cumulative_fraction`.

   | Input | Expected | Actual |
   |---|---|---|
   | `1e-13` | Reject canonical 0% | Accepted as **0.0**; **`prcp_0_dowy=1`**, before any precipitation |
   | `1-1e-13` | Reject canonical 100% | Accepted as **1.0**; **`prcp_100_dowy=365`** |

   **Proposed fix:** Compute the canonical integer percentage once, require `1 <= percentage <= 99`, then derive both fraction and name from it. Test both endpoints.

   This affects custom configurations; the shipped **0.5** setting is unaffected.

**What I verified and found correct**

- **Remaining round-1 fixes:** Negative-P years cannot enter through an implicit precipitation fallback; explicitly filtered empty input gives mean **NaN**. Truncated leap WY 2000 produces **20/20 precipitation and 19/19 temperature NaNs**. Missing month metadata is reconstructed: annual and seasonal-sum precipitation both **730**, as are GDD totals. Windows **1/367** are rejected; the shipped precipitation registry has **20 distinct names**.

- **Tolerance boundary:** With total **2 mm**, place `1−δ` on day 7 and `1+δ` on day 301. For δ=`1e-12`, crossing is **7**; for δ=`4e-12`, it is **301**. These correctly straddle the specified **2e-12** allowance.

- **Missing-observation contract:** **280 cases passed** across Q, PPT, SWE, TMIN and TMAX, using NaN, ±Inf, `missing`, strings, complex values with zero/nonzero imaginary parts, and `nothing`. Isolated and three-day internal gaps and 30 separated gaps pass; four-day gaps, boundary gaps and 31 separated gaps fail. Corresponding internal day-100 fills are **Q=100, PPT=2, SWE=1, Tmin=−1, Tmax=5**. Q diagnostics count negative infinity as missing without setting the negative-flow flag.

- **Policy validation:** All **eight negative count/gap overrides** throw `ArgumentError`. Season definitions containing fractional months, **0**, **13**, or Boolean values reject. January-only winter produces precipitation **62** and GDD **62**. Temperature interpolation allowance **0** rejects one internal missing day. PPT and temperature qualification remain independent; both require streamflow-valid years.

- **Seasons and GDD:** WY 2001 with P=month produces winter/spring/summer/fall totals **459, 368, 645, 910 mm**. Tmax=month and Tmin=month−4 produces GDD **310, 184, 461, 728**, annual **1683**. Fall correctly includes October–November of the starting calendar year and September of the ending year.

- **Calendar boundaries:** WYs **1900/2001/2100** have 365 days; WY **2000** has 366. June 30 frost gives first=last **273/274**; July 1 gives first **274/275**, last **NaN**. October 1 gives **1**; September 30 gives first **365/366**, last **NaN**. Missing December 31, 2000 is restored at dowy **92**, Tmin **−2**, Tmax **6**, within a qualifying **365-row WY 2001**.

- **Independent numerical reference:** All **39 metrics across 24 synthetic years—936 annual values—passed** rational-arithmetic comparisons. Temperature values matched exactly; maximum absolute precipitation-metric difference was **1.4551915228366852e-11**. A single 100 mm storm gives Gini **0.9972602739726028**, sample SD **5.234239225902137**, CV **19.1049731745428**, SI **1.8333333333333337**. Equal monthly totals give SI **0**.

- **Zeros, ties and thresholds:** Zero-rain years give dry length **365/366**, start **1**, end **365/366**. Equal longest two-day runs select **1–2**. Endpoint storms do not wrap: five-day maximum **10**, not 20. Exact-threshold/next-float probes give wet counts **5,3,1** above 10/25/50 mm. Signed-zero extrema select day **1**. Shuffling preserves results; duplicate dowy invalidates the year.

- **Gates and statistics:** First-decade **7/10** frost coverage gives slope **NaN**; **8/10** gives **1.0**. Means remain **36.22222222222222 / 35.92857142857143**, Pettitt year **2016**. Removing three temperature-invalid anchor years gives **7/7**, slope **1.0**. Last-decade probes also give **NaN/1.0** at **7/10 versus 8/10**, for both frost metrics. With the 20-value floor, 19 years give temperature-maximum mean **NaN**; 20 give **5.0**. Annual collection remains intact.

- **Existing families and integration:** Temperatures present versus absent preserve **976 shared summary keys**, adding **152** temperature keys. Climate enabled versus disabled preserves **704 existing keys**, adding **312** climate fields without changepoints. Common cleaned data match with both temperature columns, neither, or only one. Parquet aliases normalize correctly. Explicit qualified inputs with `area_normalized=false` give precipitation and GDD means **730**, consistent with the brief. Changepoints expand the combined climate families to **624 fields**.

- **Implementer tests:** **366/366 passed**. They do not detect the new endpoint-validation defect.

**Verdict: GO for a production rerun with the shipped configuration.** All five round-2 findings pass their original reproductions. The remaining minor finding affects custom fractions near 0% or 100%.
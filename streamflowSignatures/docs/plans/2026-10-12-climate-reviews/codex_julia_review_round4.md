**No open Julia findings.** Reviewed Julia only; no repository files changed.

1. **MINOR — round-3 canonical-percentage defect: PASS, resolved.**  
   Location: [config.jl:223](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/config.jl:223).  
   Contract: the canonical percentage must be **1–99**. Original reproduction: complete WY2001, P[365]=2 mm, zero elsewhere; fractions `1e-13` and `1−1e-13`. Expected: reject both. **Actual: both throw `LoadError` stating that the canonical percentage must be between 1 and 99.** Neither reaches calculation. Valid endpoints `0.01` and `0.99` remain accepted. The proposed post-canonicalisation range check is implemented; no further fix required.

2. **NOTE — round-3 direct complex-observation discrepancy: PASS, resolved.**  
   Location: [climate.jl:129](/Users/arik/github/SurfaceWaterProjections/streamflowSignatures/julia/src/climate.jl:129).  
   Contract: Real observations, including Bool, are accepted; other types become NaN. Original complete-WY2001 inputs TMIN=`1+0im`, TMAX=`10+0im`, and the corresponding inputs with imaginary part 2, now produce **19/19 annual temperature NaNs**, without exceptions; `valid_temp_years=[]`. Complex precipitation produces **20/20 NaNs**. Boolean controls give precipitation **365 mm** and GDD **182.5**. The proposed climate-local conversion is implemented; no further fix required.

**What I verified and found correct**

All five earlier reproductions retained in the round-3 report still pass:

| Reproduction | Result and actual values |
|---|---|
| Halfway crossing | **PASS:** seven-day blocks → **7**; 50-day blocks → **50**; uniform leap-year precipitation → **183**. |
| Non-finite PPT | **PASS:** isolated NaN/±Inf fills to **2 mm**, annual **730 mm**; four-day and boundary gaps reject the year. |
| Boolean configuration loopholes | **PASS:** Boolean numeric settings, `[true,25,50]` thresholds, and numeric 0/1 flags reject. |
| Near-whole fraction | **PASS:** `0.500000000001` → **0.5**, endpoint-storm crossing **1**; `0.504` rejects. |
| Interpolation neighbours | **PASS:** `[0,NaN,Inf]` remains unchanged; original temperature fixture gives `[0,1.6666666666666667,3.3333333333333335,5]`, GDD **725.5**. |

Additional regression checks:

- **936 independent annual comparisons passed** across 24 synthetic years. Maximum absolute precipitation error: **1.4551915228366852e−11**; temperature error: **0**. All **280 missing-data policy cases** passed.
- Seasonal precipitation: **459,368,645,910**; seasonal GDD: **310,184,461,728**, total **1683**. Single-storm Gini: **0.9972602739726028**; SI: **1.8333333333333337**.
- June 30 frost: **273/274**; July 1 last frost: **NaN**. Missing Daymet December 31 restores dowy **92**, Tmin **−2**, Tmax **6**. Dry-year, ties, thresholds, duplicate-day and truncated-leap checks passed.
- Both decade gates give slope **NaN at 7/10**, **1.0 at 8/10**. The 19/20-year statistics floor remains correct.
- Temperature presence preserves **976 shared fields**; disabling climate preserves **704 existing fields**. With `area_normalized=false`, precipitation and GDD means remain **730**. Combined climate schema: **624 fields**.
- Implementer suite: **366/366 passed**, supplementary to the independent checks.

**Verdict: GO for the Julia production rerun.**
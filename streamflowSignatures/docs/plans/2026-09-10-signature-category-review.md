# Signature category review — 2026-09-10

User request: the manuscript lists **14** signature categories, the guidelines doc **13**
modules, and a colleague proposed an **8-group** categorization
(Google Sheet `1sFVNL9bpbGXznBBfp41IjT0RqTbpRvewq1i8J0Qj5tA`, two columns
`category,signature`, 100 rows = the 100 annual signature bases, no scalars). Compare
every description of the categories and find a consistent solution.

Per-output matrix (121 rows = 100 annual bases + 21 per-gage scalars, one column per
source): `2026-09-10-signature-category-matrix.csv`. Rows are ordered by the **producing
function** (`function`, `module`, `function_order`, `signature_order` columns — the order
`calculate_all_signatures()` calls the 15 family functions in `julia/src/signatures.jl`,
signatures in each function's emission order). Two outputs come from outside the family
functions: `season_excluded_years_*` are computed in the orchestrator from the
preprocessor's seasonal flags (`signatures.jl:57`), and `ice_affected_days_total` is
written by the benchmark runner from the preprocessor's `na_cause_ice` diagnostic
(`run_julia_benchmark.jl:393`). Sources read: manuscript §2 preamble
+ §2.2.1 (live 2026-09-10), guidelines Part 3 (live 2026-09-10), README.md "Signature
Categories", docs/SIGNATURES.md summary table, `docs/plans/dataset_workflow_schematic.md`,
`julia/src/signatures.jl` family functions, `docs/benchmarks/build_signature_explorer.py`
`category_of`, `compare_experiment_vs_julia.py` `categorize_metric`,
`build_experiment_vs_julia_dashboard.py` `SIGNATURE_GROUPS` (shared by the other five
dashboard builders), the staged HydroShare `hisss_signature_categories.csv` and the
`category` column of `hisss_data_dictionary.csv` (md5-identical copies in Resources 1 and 2).

## Findings

1. **The 14-lists all agree**: manuscript, README, SIGNATURES.md summary table, schematic,
   and the code's family functions (baseflow's two functions counted once) name the same
   14: flow volumes & percentiles, flow duration curve, baseflow, recession, pulses &
   reversals, flashiness, flow timing, runoff ratios, elasticity, Q-P seasonality,
   storage, negative flow days, snow, streamflow drought.
2. **Guidelines Part 3 = the same list minus Storage** (13 modules). `avg_storage` survives
   only in the legacy tail. Already queued 2026-09-04 ("Storage module dropped"). Side
   effect: Part 1.2's "(see 3.12)" pointer for Negative_ann is off by one (3.12 is Snow) —
   the cross-reference was written against the draft's numbering with Storage at 3.11.
3. **The colleague's 8 groups are a different axis**, not a re-cut of the 14:

   | 14-category (annual bases) | → 8-group |
   |---|---|
   | Baseflow (4), Recession (7), Runoff ratios (5), Storage (1) | Storage (17) |
   | Pulses & reversals (13 of 14), Elasticity (2), Flashiness (1) | Flashiness (16) |
   | Flow volumes & percentiles (21), TQmean (1), Negative flow days (1) | Flow Volume (23) |
   | Flow timing (15) | Flow Timing (15) |
   | Snow (14) | Snow (14) |
   | Streamflow drought (10) | Drought (10) |
   | Flow duration curve (3) | Flow Duration (3) |
   | Q-P seasonality (2) | Precipitation Streamflow (2) |

   The only 14-category split across two groups is Pulses & reversals (TQmean → Flow Volume).
4. **The sheet is a revision of the staged 9-class exploratory grouping** (identical
   100-base set; 7 rows differ): `flashinessRB` Storage → Flashiness (the staged placement
   is an error), and the six "Snow Timing" bases (`swe_max_dowy`, `snow_on_dowy`,
   `snow_off_dowy`, `melt_season_days`, `ssm`, `melt_com_dowy`) fold into "Snow".
5. **The inconsistency that matters**: manuscript §2.2.1 says "14 categories … details
   in hisss_data_dictionary.csv", but the dictionary's only `category` column carries the
   9-class grouping. The Resource 1/2 READMEs call the CSV a "coarse nine-class exploratory
   grouping" and the paper's categories "a separate taxonomy" — but nothing shipped
   carries the 14.
6. Repo tooling drift (minor, all internal): dashboard `SIGNATURE_GROUPS` puts
   `negative_ann` in "Pulses" and cannot select the 21 scalars; the explorer and compare
   categorizers split "Flow percentiles" from "Flow volumes" (16 categories) and drop
   `ice_affected_days_total`, `season_excluded_years_*` (and, in the explorer,
   `recession_alpha_point_cloud_linear_reservoir`) to "Other"; the claude-skill overview
   lists six categories; SIGNATURES.md has no numbered Negative Flow Days section (summary
   table only).

## Proposed solution (pending user decision)

Keep both axes under distinct names and use each consistently:

- **`category`** (14, methodological) — manuscript §2.2.1, guidelines Part 3 modules,
  code families, README, SIGNATURES.md, and a `category` column in the dictionary.
- **`analysis_group`** (8, the colleague's sheet) — figures, exploratory dashboards,
  `hisss_signature_categories.csv`, and an `analysis_group` column in the dictionary.

Steps: (1) restore a Storage module in guidelines Part 3 (→ 14; renumber 3.11–3.13 to
3.12–3.14, which repairs the "(see 3.12)" pointer); (2) adopt the sheet over the staged
9-class file and extend it to the 21 scalars by inheriting the family's group (drought
thresholds → Drought; recession seasonality, recession alpha, runoff_ratio_high_count →
Storage; elasticity scalars → Flashiness; season_excluded_years_* → Flow Volume;
ice_affected_days_total → none, it is a preprocessing diagnostic); (3) regenerate the
dictionary with both columns and reword the Resource 1/2 README rows; (4) one bridging
sentence in the manuscript where the §4 figure uses the 8 groups; (5) align the repo
tooling listed in finding 6.

Judgment calls to raise with the colleague (not errors): elasticity under "Flashiness"
while runoff ratios are "Storage" and Q-P seasonality is "Precipitation Streamflow" —
three precipitation–streamflow families in three groups; TQmean (a duration measure) and
negative_ann (a data-quality count) under "Flow Volume".

## Afternoon crosswalk (2026-09-10 pm) — after the user restructured the guidelines doc to 8 modules

The guidelines doc now has eight Part 3 modules matching the colleague's sheet (3.1 Flow
Volume, 3.2 Flow Duration, 3.3 Storage, 3.4 Flashiness, 3.5 Drought, 3.6 Flow Timing,
3.7 Precipitation Streamflow, 3.8 Snow); the legacy tail is gone. Manuscript unchanged
since the morning sync. Live HydroShare pages are private (403) and Chrome was not
connected, so the deposit side is the STAGED upload folder (`~/Downloads/Signatures/`).
Per-output crosswalk: `2026-09-10-signature-category-crosswalk-pm.csv` (121 rows:
manuscript 14-category, guidelines module, sheet group, dictionary class, flags).

Inconsistencies found, by product:

**Manuscript (14-language throughout)** — §2 preamble "121 signatures across 14
streamflow and hydroclimate categories"; §2.2.1 "spans 14 categories: flow volumes and
percentiles, flow timing, flow duration curve slopes, baseflow, recession behavior,
high-and low-flow pulses, flashiness, runoff ratios, streamflow elasticity, P-Q
seasonality, catchment storage, negative-flow days, snow metrics, and streamflow
drought" + "details … in hisss_data_dictionary.csv"; §4 figure "sites and signature
categories" (scheme unknown). None of these matches the doc's eight.

**HydroShare deposit (staged R1/R2)** —
- `hisss_data_dictionary.csv` `category` column and `hisss_signature_categories.csv`
  (identical copies at the root, R1 and R2): the 9-class grouping. Versus the doc/sheet,
  7 bases differ (flashinessRB → Storage; swe_max_dowy, snow_on_dowy, snow_off_dowy,
  melt_season_days, ssm, melt_com_dowy → "Snow Timing"). The 21 scalars are
  `category = "Per-gage scalar"` in the dictionary and absent from the categories CSV;
  the doc places 20 of them in modules (3.1 ×4, 3.3 ×8, 3.4 ×3, 3.5 ×5).
- README lede (R1 and R2, identical): "spanning flow volumes and percentiles, timing,
  flow-duration curves, baseflow, recession, pulses, flashiness, runoff ratios,
  elasticity, precipitation–streamflow seasonality, storage, negative-flow days, snow,
  and streamflow drought" — the 14 list.
- README file table: "Coarse nine-class exploratory grouping of the 100 signature bases
  (the data paper's finer methodological categories are a separate taxonomy)" — nine is
  now wrong, and the "finer taxonomy" no longer exists in the doc.
- Shipped explorer HTML (`hisss_signature_explorer_*.html`): category picker uses the
  repo's 16-class `category_of` scheme (Flow percentiles split from Flow volumes, Pulses,
  Negative flow, Other …). Shipped validation summary MD and dashboard HTML: the
  comparator's 16 / 15-group schemes ("Flow Percentiles", "Pulse Metrics", "Other").
- Collection abstract: "roughly 100 hydrological signatures" — no category count (fine).

**Guidelines doc, internal** —
- `avg_storage` / `calculate_average_storage()` absent, while the sheet and dictionary
  place avg_storage in Storage and the manuscript lists catchment storage.
- TQmean: doc 3.4 Flashiness (with calculate_pulse_metrics) vs sheet + dictionary Flow
  Volume. The only annual base where doc and sheet disagree.
- Part 1.2 "counted in Negative_ann (see 3.12)" → 3.1.
- `ice_affected_days_total`: shipped scalar with no module (Part 1.2 asks for the count).
- Names vs shipped columns: "recession_alpha_point_cloud" (column is
  `recession_alpha_point_cloud_linear_reservoir`); "Negative_ann" (column `negative_ann`).
- "family" is used for the sub-blocks (recession, elasticity, parameterized-BFI) inside
  modules; intro says "each signature family as a self-contained module" — define once.

**Repo docs on the public mirror (linked from the manuscript)** — README.md "Signature
Categories" table (14), SIGNATURES.md (13 numbered sections + 14-row summary table),
claude-skill overview (6), CLAUDE.md "(14 categories)"; dashboards 15/16 groups.

Resolution requires one vocabulary decision: if "category" = the eight, then the
manuscript (§2 preamble, §2.2.1, §4 figure), the dictionary/categories CSV (regenerate
with the 8 groups + scalars), the README lede/file-table rows, and the repo docs all
change; the 15 functions can still be named as the computing units. TQmean and
avg_storage need a ruling before the CSV is regenerated.

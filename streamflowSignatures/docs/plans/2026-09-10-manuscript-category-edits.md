# Signature-category alignment — manuscript and guidelines-doc edits (2026-09-10)

**Reference**: the co-authors' 8-category sheet (`1sFVNL9bp…`), now the canonical
grouping — repo copy `docs/signature_categories.csv` (121 rows, scalars assigned by their
producing function; `ice_affected_days_total` stays outside as a preprocessing
diagnostic). **Vocabulary**: *category* = one of the eight; *function family* = the
computing function (15 functions of `calculate_all_signatures()`). Repo docs were aligned
in this session (README.md, docs/SIGNATURES.md, CLAUDE.md, claude-skill, the workflow
schematic, CROSS_LANGUAGE_STATUS note, explorer builder). The two Google Docs cannot be
edited from here — catalogue below. **HydroShare files are NOT changed this session** —
see the pending list at the end.

## A. Manuscript — every location that must change

| # | Location | Current text | Change to |
|---|----------|--------------|-----------|
| 1 | §2 preamble | "we computed 121 signatures across 14 streamflow and hydroclimate categories: 100 annually resolved signatures … plus 21 per-gage diagnostics" | "… across eight signature categories (Flow Volume, Flow Duration, Storage, Flashiness, Drought, Flow Timing, Precipitation Streamflow, Snow): 100 annually resolved signatures … plus 21 per-gage diagnostics" |
| 2 | §2.2.1, sentence 2 | "spans 14 categories: flow volumes and percentiles, flow timing, flow duration curve slopes, baseflow, recession behavior, high-and low-flow pulses, flashiness, runoff ratios, streamflow elasticity, precipitation-streamflow (P-Q) seasonality, catchment storage, negative-flow days, snow metrics, and streamflow drought" | "spans eight categories: flow volume (annual and seasonal totals, percentiles, negative-flow days), flow duration, storage (baseflow, recession, runoff ratios, catchment storage), flashiness (pulses, reversals, the Richards-Baker index, elasticity), drought, flow timing, precipitation–streamflow seasonality, and snow" |
| 3 | §2.2.1, "Details for each metric … are in the data dictionary" | (unchanged) | keep, but only once the dictionary's `category` column carries the eight (HydroShare pending item H1) |
| 4 | §3 Resources 1–2, item (i) | "the 100 annually resolved signatures x 16 statistics each, 21 signatures and related outputs that do not carry statistics" | keep; optionally add "each labeled with its signature category in the data dictionary" |
| 5 | §4 Data Overview, figure "summarizing sites and signature categories" | scheme unknown | use the eight categories and their names verbatim; the workflow schematic (`docs/plans/dataset_workflow_schematic.md`) now reads "121 signatures in 8 categories" |
| 6 | §5.1.1 | "precipitation-dependent signatures (runoff ratios, elasticity, Q-P seasonality, storage)" | keep (metric names, not categories); if "storage" is read as the category, write "catchment storage (avg_storage)" |
| 7 | §5.1.3 | "the recession-parameterized baseflow indices" | keep (metric names) |
| 8 | §1, "roughly 100 hydrological streamflow signatures" | — | optional: "roughly 100 … signatures in eight categories" |

No other sentence in the current draft names a category count or list.

## B. Guidelines doc — residual edits so Part 3 matches the sheet exactly

| # | Location | Current | Change to |
|---|----------|---------|-----------|
| 1 | 3.4 Flashiness, `calculate_pulse_metrics()` block | `TQmean: Percentage of days with flow above the annual mean.` listed under Flashiness | move the TQmean line to 3.1 Flow Volume (add a note "computed by calculate_pulse_metrics(); reported under Flow Volume") — the sheet and the dictionary put TQmean in Flow Volume |
| 2 | 3.3 Storage | no `calculate_average_storage()` block | add the block: `avg_storage`: mean annual catchment storage from S = cumsum(P − Q), interpolated at mean discharge; requires PPT + area-normalized flow; computed but omitted from major analyses (no ET term). Sheet, dictionary and manuscript all carry it |
| 3 | Part 1.2 | "counted in Negative_ann (see 3.12)" | "(see 3.1)" |
| 4 | 3.1 Flow Volume (or Part 1.2) | `ice_affected_days_total` not documented | add one line: per-gage diagnostic scalar, count of ice-flagged NA days from the preprocessor (structurally 0 in the delivered products — known issue) |
| 5 | 3.3 Storage, baseflow block | `recession_alpha_point_cloud` | `recession_alpha_point_cloud_linear_reservoir` (the shipped column name) |
| 6 | 3.1 and Part 1.2 | `Negative_ann` | `negative_ann` (shipped column name) |
| 7 | Intro + Part 1.3 | "each signature family as a self-contained module"; "recession and elasticity families" | say once that a *module* = one of the eight categories and a *family* = one computing function inside it |
| 8 | Part 5 (already queued 2026-09-04) | `flagged_for_high_na` "numeric output columns"; "identical across the Julia, Python, and R implementations" | the three queued high-NA edits |

## C. HydroShare documents — pending, NOT changed this session

Terminology and content updates needed in the staged deposit (and on HydroShare once
uploaded), to be applied in a later session:

- **H1 — category terminology.** `hisss_data_dictionary.csv` `category` column and
  `hisss_signature_categories.csv` (root, R1, R2 — identical copies) carry the old 9-class
  grouping ("Snow Timing" separate; flashinessRB under Storage; scalars as "Per-gage
  scalar"): regenerate both from `docs/signature_categories.csv` (121 rows incl. the
  scalars). README R1/R2: the lede's 14-family list → the eight categories; the file-table
  row "Coarse nine-class exploratory grouping … finer methodological categories are a
  separate taxonomy" → "Signature category of each of the 121 outputs (eight categories;
  the same grouping the data paper uses)". The shipped explorer HTML (category picker)
  and the validation summary/dashboard embed the repo's function-family groupings — rebuild
  the explorer with the updated builder; relabel the validation tables "by function
  family" or leave them as diagnostics with a note.
- **H2 — the other pending update already on record**: the `flagged_for_high_na`
  known-issue text (README R1/R2 conventions + the dictionary row) must be rewritten when
  the column is regenerated at the next data rerun (CHANGELOG → Planned / Known Issues).
  [The user noted a further aspect of the HydroShare docs also needs updating — add it
  here when specified.]
- **H3 — HYDAT provenance (added 2026-09-10 evening).** R3 README ("Provenance and
  conventions") and the input dictionary should state the HYDAT release used —
  **2025-10-14**, the quarterly release current at the 2026-02-07 retrieval (inferred; see
  docs/DATA_SOURCES.md) — and that Canadian records end 2024-12-31, so no Canadian gage has
  a qualifying WY 2025 in either product. Optional one-liner in the R1/R2 READMEs.

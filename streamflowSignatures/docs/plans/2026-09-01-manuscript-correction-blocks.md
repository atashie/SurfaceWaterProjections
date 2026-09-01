# Manuscript correction blocks — items 1–7 of the 2026-09-01 review

Paste-ready replacement text for the seven substantive technical inaccuracies in the
HISSS manuscript draft (Google Doc), drafted 2026-09-01 at the user's request after
the fourth reconciliation pass (CHANGELOG → Manuscript Reconciliation Log →
2026-09-01). Every number below is verified against the repo record (SIGNATURES.md,
DEVELOPMENT.md, CHANGELOG August 2026). The "Current text" quotes are from the
published doc as of 2026-09-01 (curly quotes/dashes normalized); locate by the
opening words.

Also delivered as a paste-ready artifact page (same content).

---

## Block 1a — §2 preamble: signature count

**REPLACE** (mid-paragraph, "We assembled daily mean streamflow records…"):

> From each standardized series, we computed 106 hydrological signatures across 13
> categories and summarized each signature's annual time series using a common set
> of descriptive statistics, trend estimators, and changepoint tests.

**WITH:**

> From each standardized series, we computed 121 hydrological signatures across 14
> categories — 100 annually resolved signatures, each summarized with a common set
> of 16 statistics (descriptive statistics, trend estimators, and Pettitt
> changepoint tests), plus 21 per-gage scalar diagnostics.

*Why: the 106/90/16/13 counts describe the pre-drought product. Delivered products:
100 annually resolved bases + 21 scalars = 121, across 14 categories — matching the
manuscript's own (correct) §3.*

## Block 1b — §2.2.1 first paragraph: counts, categories, dictionary placeholder

**REPLACE** (entire paragraph):

> We calculated a comprehensive suite of 106 streamflow and hydroclimate signatures
> for each gage and its associated watershed, comprising 90 annually resolved
> metrics with statistics and 16 stand-alone diagnostic metrics. The selection of
> metrics is consistent with analysis in catchment hydrology (McMillan, 2021) and
> snow hydrology (??Hatchett, 2021; Petersky and Harpold, 2018??) and spans 13
> categories: flow volumes and percentiles, flow timing, flow duration curve
> slopes, baseflow, recession behavior, high-and low-flow pulses, flashiness,
> runoff ratios, streamflow elasticity, precipitation-streamflow (P-Q) seasonality,
> negative-flow days, and snow metrics. Details for each metric, including the
> definition, requirements, and relevant citations, are in the metadata file (name).

**WITH:**

> We calculated a comprehensive suite of 121 streamflow and hydroclimate signatures
> for each gage and its associated watershed, comprising 100 annually resolved
> metrics with statistics and 21 stand-alone per-gage diagnostic metrics. The
> selection of metrics is consistent with analyses in catchment hydrology
> (McMillan, 2021), snow hydrology (Hatchett, 2021; Petersky and Harpold, 2018),
> and streamflow drought characterization (Adelsperger et al., in review), and
> spans 14 categories: flow volumes and percentiles, flow timing, flow duration
> curve slopes, baseflow, recession behavior, high- and low-flow pulses,
> flashiness, runoff ratios, streamflow elasticity, precipitation-streamflow (P-Q)
> seasonality, catchment storage, negative-flow days, snow metrics, and streamflow
> drought. Details for each metric, including the definition, requirements, and
> relevant citations, are in the machine-readable data dictionary
> (hisss_data_dictionary.csv) included with Resources 1 and 2.

*Why: fixes the counts, adds the two missing categories (storage, drought), removes
the "??" placeholder markers (both snow citations are correct as written), and
resolves the "(name)" placeholder with the actual dictionary filename.*

## Block 1c — §2.2.1: NEW streamflow-drought methods paragraph

**ADD** as a new paragraph in §2.2.1 (after the BFI paragraph; before Block 6's
storage paragraph):

> We characterized streamflow drought following the severity-based approach of
> Adelsperger et al. (in review): per-water-year drought duration (days below
> threshold) and drought deficit (the summed departure below the threshold, in mm)
> at five fixed severity levels — the 2, 5, 10, 20, and 30% magnitude percentiles
> of flow, mirroring the U.S. Drought Monitor classes D4 (exceptional drought)
> through D0 (abnormally dry). Daily flow was first smoothed with a 7-day centered
> moving average applied within continuous runs of consecutive dates (never across
> a gap left by a rejected year), and thresholds were computed from the smoothed
> flow pooled over each analysis window's own record, using the unbiased Weibull
> plotting position (Laaha et al. 2017). The five threshold values are reported per
> gage as auditable scalar diagnostics. We deviate from the source method in two
> documented ways: only the fixed (whole-record) thresholds are implemented — the
> variable day-of-year thresholds are too uncertain at the low severity levels for
> records of 20–46 years, with the 2% level falling below the smallest plotting
> position entirely — and aggregation follows the water year (October–September)
> used throughout this dataset rather than the source's April–March climate year,
> so a drought spanning 1 October is split across two annual values. Because
> thresholds derive from each window's own record, drought values are comparable
> within an analysis window but must not be compared across the two windows.

## Block 2a — §2 preamble + §2.1.2: area-normalization source

**REPLACE** (§2 preamble, third sentence of "We assembled daily mean streamflow
records…"):

> The contributing watershed area was extracted for each gage and was used to
> area-normalize discharge.

**WITH:**

> Discharge was area-normalized using the drainage area published by each agency
> (GAGES-II attributes for USGS gages; the HYDAT station metadata for Canadian
> gages), rather than areas derived from the boundary polygons; the polygons were
> used to aggregate gridded climate, land-cover, and basin-attribute data
> (Sects. 2.1.3–2.1.4).

**REPLACE** (§2.1.2, mid-paragraph):

> …and area-normalized discharge (mm d-1) for all but 73 gages, as noted above
> (Sect. 2.1.1).

**WITH:**

> …and area-normalized discharge (mm d-1) using the agency-published drainage
> areas for all but 73 gages, for which no drainage area is published (see Usage
> Notes).

## Block 2b — §5 Usage Notes: un-normalized gages + record-dependent signatures

**ADD** to §5 (currently a stub):

> Un-normalized gages. 73 of the 8,014 processed gages have no agency-published
> drainage area (32 and 28 of them qualify for the WY 1993–2025 and WY 1980–2025
> products, respectively) — most are irrigation or diversion canals, dam and
> powerhouse outflows, or channel splits of large rivers, where a contributing
> area is undefined or unpublished. These gages are retained with discharge in
> native m3 s-1 and flagged area_normalized = FALSE. Their unit-carrying
> signatures (flow volumes, flow percentiles, recession log(a), drought deficits)
> are not comparable with the mm d-1 values at other gages, and their
> precipitation-dependent signatures (runoff ratios, elasticity, Q-P seasonality,
> storage) are structurally NA. Users should filter on area_normalized == TRUE
> before any cross-gage comparison of unit-carrying signatures; the automated
> range flags do not catch all of these gages.
>
> Record-dependent signatures. Signatures whose definition uses thresholds or
> means from the full analysis window — the period-of-record (\*_all) pulse
> metrics, elasticity, the recession-parameterized baseflow indices, and all
> drought metrics — are valid within a product but must not be compared between
> the WY 1993–2025 and WY 1980–2025 products, nor re-derived from the annual
> values over a different window.

## Block 3 — §2.1.2: 47% arithmetic

**REPLACE:**

> Gages that returned fewer than 20 total years containing valid daily
> observations (8,980) were excluded. This resulted in a 47% reduction in valid
> gages, with the remaining 8,014 gages (6,160 US; 1,854 Canadian) returning
> usable daily records.

**WITH:**

> Gages that returned fewer than 20 total years containing valid daily
> observations (8,980 of 16,994; a 53% reduction) were excluded, leaving 8,014
> gages (6,160 US; 1,854 Canadian) with usable daily records.

## Block 4 — §2.2.1: DELETE the duplicate BFI/statistics paragraph

**DELETE** (entire paragraph — everything it says is stated correctly elsewhere:
BFI variants in the paragraph directly above it, the two windows in §2.1.2, the 16
statistics in §3):

> For the Eckhardt filter, we used both default parameters (BFImax = 0.8, a =
> 0.98) and a parametrized version that includes recession parameters calculated
> for each site. We then calculated statistics for each metric across both the
> full period of record and a subset from 1993-2025 water years, which represent
> the highest proportion of sites, including slope, Spearman's Rho, p-value, mean,
> and median.

**OPTIONAL replacement** if §2.2.1 should still enumerate the statistics:

> For every annually resolved signature we computed the same 16 statistics over
> each analysis window: mean and median; Theil-Sen and ordinary least-squares
> trend slopes; Spearman and Mann-Kendall correlations with their p-values; and
> eight Pettitt changepoint fields (changepoint year and p-value; pre- and
> post-changepoint means, their difference, and percent change; and within-segment
> Mann-Kendall p-values).

## Block 5 — §2.1.3: Daymet coverage

**REPLACE:**

> Daymet aggregation was performed for a subset of 6,087 basins that both met both
> the aforementioned data quality standards and were less than 85,000 km2.

**WITH:**

> Daymet aggregation was performed for a subset of 6,087 candidate basins smaller
> than 85,000 km2. Of the 8,014 gages with usable streamflow records, 5,965 have
> basin-averaged Daymet series (5,517 and 5,638 of the gages in the WY 1993–2025
> and WY 1980–2025 signature products, respectively); climate-, snow-, and
> precipitation-dependent signatures are reported as NA for the remainder.

*Why: 122 of the 6,087 aggregated basins belong to candidate gages that never
yielded usable streamflow records, so "met the data quality standards" overstates;
Daymet coverage is a strict subset of each gage set. Also fixes the "both met
both" typo. The following HUC8-threshold sentence stands unchanged.*

## Block 6 — §2.2.1: reinstate avg_storage

**ADD** as a new paragraph in §2.2.1 (e.g., after Block 1c's drought paragraph):

> We additionally computed a mean annual catchment storage signature (avg_storage)
> from the simplified water balance S = Σ(P − Q), with annual storage interpolated
> at mean discharge and averaged across years (Peters and Aulenbach 2011). Because
> this water balance carries no evapotranspiration term, it overestimates storage
> wherever ET losses are significant; avg_storage is therefore retained in the
> data product for completeness but omitted from the major analyses. It is
> computed only for gages with valid precipitation data and area-normalized flow.

## Block 7 — Acknowledgements: AI-tool attribution

**REPLACE** (both duplicated passages — the "AI-assisted coding tools…" sentence
pair AND the "Generative AI tools…" sentence pair):

> AI-assisted coding tools (Claude Code 0.145.0, Anthropic) were employed to
> generate code used in data ingestion and processing, hydrological signature
> extraction, cross-language validation and benchmarking, and interactive
> visualization.. All code was reviewed, tested, and validated by the authors to
> ensure correctness and reproducibility. Generative AI tools (Claude Code
> 0.145.0, Anthropic) were used to support data analysis and visualization. These
> tools were applied under the supervision of the authors, and all outputs were
> reviewed and validated against established scientific methods to ensure
> reproducibility and transparency.

**WITH** (single merged paragraph):

> AI-assisted coding tools (Claude Code, Anthropic) were employed under the
> supervision of the authors to generate code used in data ingestion and
> processing, hydrological signature extraction, cross-language validation and
> benchmarking, data analysis, and interactive visualization. Independent
> adversarial reviews of the code and outputs were additionally performed with a
> second AI tool (Codex CLI 0.145.0, OpenAI). All code and outputs were reviewed,
> tested, and validated by the authors against established scientific methods to
> ensure correctness, reproducibility, and transparency.

*Why: 0.145.0 is the codex-cli version from this repo's adversarial-review records,
not a Claude Code version, and the two AI sentences were near-duplicates. If the
journal requires a version for Claude Code, the current one on the project machine
is 2.1.252 — but the work spanned many releases, so naming none (or "v1.x–2.x") is
also defensible.*

---

## Reference-list entries these blocks introduce

The blocks above cite four works not yet in the reference list (part of open item
8 from the 2026-08-28 pass):

- Adelsperger, S., et al. (in review). A novel severity-based approach for
  assessing streamflow drought characteristics and drivers. *(update when the
  paper lands)*
- Hatchett, B.J. (2021). Seasonal and Ephemeral Snowpacks of the Conterminous
  United States. Hydrology, 8(1), 32. https://doi.org/10.3390/hydrology8010032
- Laaha, G., et al. (2017). Unbiased plotting positions for low-flow frequency
  analysis. *(journal details to be completed — cited in the guidelines doc in
  this short form)*
- Peters, N.E., & Aulenbach, B.T. (2011). Water storage at the Panola Mountain
  Research Watershed, Georgia, USA. Hydrological Processes, 25(25), 3878–3889.
- Petersky, R., & Harpold, A. (2018). Now you see it, now you don't: a case study
  of ephemeral snowpacks and soil moisture response in the Great Basin, USA.
  Hydrology and Earth System Sciences, 22, 4891–4906.
  https://doi.org/10.5194/hess-22-4891-2018

(McMillan 2021 is retained by Block 1b and is also missing from the list — it is
already tracked under open item 8 with the other missing references.)

# Adversarial review of one HydroShare resource of the HISSS dataset — ACCURACY and COMPLETENESS

You are an independent, adversarial reviewer. Your job is to find every way in which the deposited files of ONE
resource are inaccurate (a claim in a document that the data contradicts; a value, unit, count, id format,
period or column that is wrong or inconsistent) or incomplete (a file, column, gage, year, caveat or definition
that is promised or needed but absent). Assume mistakes exist and go looking for them with the data, not by
reading alone. Every finding must be backed by evidence you computed or quoted (file, line/column, numbers).

## Hard rules
- READ-ONLY. Do not create, modify, move or delete any file anywhere. Do not run anything that writes outside
  /tmp. Do not use the network.
- Python: use exactly `/Users/arik/HISSS_data/venvs/daymet/bin/python` (pandas 3, pyarrow 25, geopandas 1.2,
  pyogrio, duckdb 1.5). If you use DuckDB, run `SET memory_limit='1GB'; SET threads=4;` first. `sqlite3` and
  `unzip -l` are available. Keep peak memory under ~3 GB: read Parquet by row group / column projection, never
  load the 5 GB Daymet table or a whole explorer HTML naively into pandas; one full pass over a large file at most.
- Work in the review tree (symlinks to the verified upload copies — these ARE the files now live on HydroShare,
  byte-verified today): /private/tmp/claude-501/-Users-arik-github-SurfaceWaterProjections-streamflowSignatures/be13f632-088d-489b-9c6c-c121142b0d4a/scratchpad/codex_reviews/tree/
  Your working directory is the project repository (code, docs, config): use it as reference material
  (docs/SIGNATURES.md, docs/signature_categories.csv, config/signatures_config.json, docs/DATA_SOURCES.md,
  docs/DEVELOPMENT.md, docs/QUALIFIER_RERUN_RESULTS.md, EO_data_processing/README*.md, the staging tools
  docs/benchmarks/stage_hydroshare_*.py and the readme template docs/benchmarks/hydroshare_readme_signatures.tmpl.md).
  The deposited files are the subject of the review; the repo is context. Where the repo and a deposited file
  disagree, report it — do not assume either side is right.
- Budget: aim to finish within ~45 minutes of wall time. Prioritise checks by the damage a wrong public claim
  would do (counts, units, ids, join rules, missing files) over cosmetics.

## The dataset (context)
HISSS = five HydroShare resources (private, about to be published as the FIRST AND ONLY public version):
  1. Signatures Product 1993-2025 (84edce76f22e4aaebd3f50e508bacfbe)
  2. Signatures Product 1980-2025 (e86b8c60fa454c18bce642b035b39ece)
  3. Harmonized Input Data (e238a594bfcd4719ab4bcab232ac48dc)
  4. Watershed Geometry and Basin Attributes (1283981f7f6a451a91f745902256443e)
  5. Vegetation and Land Cover (27da62083ba440ebbda8b4a805aa2b1d)
  collection f702201faa5d46069a5ee83ffa4c9768.
Streamflow signatures: 139 annual signatures (100 hydrological + 39 climate) per gage, each summarised by 16
statistics (8 statistics + 8 Pettitt fields) — the "wide" summary CSV — plus a long annual-values Parquet.
Water year = Oct 1–Sep 30; flow in mm/day except Canadian gages with no published drainage area
(`area_normalized = FALSE`, raw m³/s, Q-to-PPT signatures NA by design). Gage qualification: ≥ 20 water years and
≥ 60 % of the window. Signature categories: the canonical eight (Flow Volume, Flow Duration, Storage,
Flashiness, Drought, Flow Timing, Precipitation Streamflow, Snow) plus Climate; "function family" is the sub-level.

## Decisions already taken by the data owner (do NOT report these as errors; DO report files that contradict them)
- First and only public version: no user-facing text may say "previous build", "supersedes", "corrected (this
  version)", "regenerated", "earlier release", "version 2/3" or otherwise imply an earlier public release.
  Internal provenance (run dates, pipeline names) is allowed. FLAG any such wording you find.
- Gage sets are frozen: gages that would newly qualify are LISTED (Harmonized Input Data:
  hisss_would_qualify_not_added.csv) but not added.
- USGS gage ids: the input-data gage-metadata table stores EVERY USGS id without leading zeros; the signature
  tables and the geometry/EO tables keep ids as delivered (some 8-digit ids with their leading zero, 44 (+9)
  long ids un-padded). The documented join rule is "strip leading zeros on both sides". Test that this rule
  actually joins the tables; report any id that still fails to join.
- Climate signatures are computed on the OUTER streamflow–climate date frame (annual totals = sums of all daily
  values, not only of days with streamflow).
- Daymet: basin-averaged over the UNSIMPLIFIED polygons (true-area weights), 8,017 sites, 1980–2025, with a flags
  table for questionable basins.
- EO products: MODIS LAI/LULC cover 7,964 watersheds; NLCD excludes 45 Alaska gages (listed in a companion
  file); the boundary layer's `low_confidence` (HydroBASINS fallback OR area outlier) is a different flag from
  the LAI/LULC tables' own `low_confidence` columns (which also fold in pixel support); explorers show the
  deposited boundary layer's flag.
- Known and accepted: `flagged_for_high_na` was regenerated; `ice_affected_days_total` is 0 for every gage
  (documented); 37 Canadian gages carry raw m³/s; the 1,653-column / 100-signature / "14 categories" wording
  belongs to an OLDER internal build and is WRONG for these files (2,277 columns, 139 signatures, 8 + Climate).
- TODAY (2026-10-30 ~11:45 EDT) the owner renamed all five resources and the collection to
  "Hydroclimate Information, Signatures and Summary Statistics - …". The files still expand HISSS as
  "Hydrologic Information Signatures and Summary Statistics". List every file/line that carries the old
  expansion as ONE finding (severity MINOR) and move on.

## What to check — ACCURACY
1. Every quantitative or categorical claim in the readme, the data dictionary, the categories file and the LIVE
   ABSTRACT (given below; it is the text on the resource page) against the data files: gage counts, row counts,
   column counts, periods/years, months, units, id formats, number of flags/categories/statistics, percentages,
   named exceptions (e.g. counts of Canadian raw-m³/s gages, HydroBASINS fallbacks, excluded gages, large basins).
2. Column-level: every column of every data table is defined in the dictionary (name, type, unit, description)
   and every dictionary row corresponds to an existing column; stated value ranges / allowed values / NA rules
   hold in the data (test them); stated sums (e.g. percent columns summing to 100) hold.
3. Signature tables: the 16-statistic column pattern holds for every annual signature; single-valued
   exceptions match docs/SIGNATURES.md and config; the categories file maps every signature exactly once and
   its categories/function families agree with the dictionary and with docs/signature_categories.csv; the
   annual Parquet's signatures, gages and water years are consistent with the summary CSV (same gage set, the
   summary statistics reproducible for a sample of gage×signature pairs, e.g. `_mean`/`_median` from the annual values).
4. Join rules and cross-resource consistency: the gage ids of this resource against the other resources' key
   files (paths below) under the documented join rule; counts of unmatched ids each way.
5. Explorer HTML files (if any): self-contained? which tile/basemap services do they reference (no keyless CARTO
   tiles; expect Esri/OSM)? does the embedded data agree with the tables (gage count, variable list, a sampled
   value)? any "previous version" wording? The sidecar ZIP (if any): list its contents and check the readme's
   description of it.
6. GeoPackage/GeoParquet (if any): layer names, feature counts, CRS, geometry validity on a sample, id
   uniqueness, attribute columns vs dictionary, area fields vs geometry area on a sample (document the CRS used).
7. Dates/times: every "as of" date, release note (e.g. HYDAT release), period of record and
   retrieval-extent claim is internally consistent across the resource's files.

## What to check — COMPLETENESS
1. Every file the readme/abstract mentions exists in THIS resource or is correctly attributed to a named
   companion resource (companion resource file lists are given below). Every file in the resource is described.
2. Every gage in the resource's gage set appears in every table that claims to cover "all gages"; missing years /
   months / gages are either documented or reported by you (quantify).
3. Required caveats are present and correct (units exceptions, flags, id padding, low-confidence semantics,
   Alaska exclusions, HYDAT release note, qualification rule, climate frame, NA conventions).
4. Dictionary completeness: no undefined columns, no undefined code values (e.g. flag codes, source codes,
   class codes, scheme names), units given wherever a column is numeric.

## Output — print this report as your final answer (Markdown, ≤ 300 lines)
# Review: <resource title>
**Verdict:** CLEAR | FINDINGS (n BLOCKER / n MAJOR / n MINOR / n NOTE)
## Findings (ranked)
| # | Severity | File (location) | Claim or expectation | Evidence (what you computed/quoted) | Recommended fix |
  Severity: BLOCKER = a public user would be misled or could not use the data (wrong counts/units/ids/missing
  file); MAJOR = material inaccuracy or omission; MINOR = wording/consistency; NOTE = observation, no change needed.
## Checks performed (with the numbers you obtained — so the owner can audit you)
## Not checked (and why)
Do not pad the report; do not restate the brief; no speculation without evidence.

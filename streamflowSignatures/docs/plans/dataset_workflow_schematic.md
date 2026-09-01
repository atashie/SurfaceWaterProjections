# HISSS Dataset Workflow Schematic

Workflow figure for the HISSS manuscript (§2 Methods) — from raw sources through
QA screening, signature extraction, and watershed aggregation to the five
HydroShare resources' core tables. GitHub renders the diagram below natively.

- **Corrected 2026-09-01** (from the draft schematic): 121 signatures
  (100 annual series + 21 scalars) across **14** families (was "~100 across 13");
  trend gate stated as ≥20 values / **≥60% of series** / ≥80% of first and last
  decade (was "≥80% of series and each decade"); **Annual NLCD** added as a raw
  source and to the land-cover output table; Daymet source noted as pre-computed
  area-weighted basin averages; "gauge" → "gage" (USGS convention); the open
  library linked as "reproducible via" rather than a downstream consumer.
- **GitHub rendering note (2026-09-01)**: "Unable to render rich display" on this
  file was NOT caused by this source. Verified directly: the diagram renders
  cleanly under stock mermaid 10.9.3, 11.12.0, AND 11.17.2 — GitHub's exact
  pinned version, extracted from its live viewscreen bundle — in strict security
  mode. GitHub's production bundle
  (`viewscreen.githubusercontent.com/static/assets/mermaidMarkdown-7fad2455….js`)
  currently crashes during its own bootstrap (TypeError in its octicons module
  init, before mermaid ever runs), reproduced on GitHub's own public mermaid
  README where ALL diagrams fail with the identical symptom — a platform-wide
  outage; re-check after GitHub ships a fixed bundle (hard refresh, the old
  bundle may be cached). The source is nonetheless kept directive-free and
  edge-label-free as the maximally compatible form (an earlier GitHub-pinned
  mermaid had a dagre bug with `%%{init}%%` + labeled edges,
  mermaid-js/mermaid#6452 / #6022); branch outcomes live in node text instead
  ("Rejected:", "Kept:", "Not normalized:").
- **PNG render**: `dataset_workflow_schematic.png` (same folder, 3× resolution,
  white ground) — regenerate with `render_mermaid.py` (same folder; needs a venv
  with `playwright` + `playwright install chromium`, plus `mermaid.min.js`
  downloaded beside the script — cdnjs mermaid 10.9.3 UMD):
  `python render_mermaid.py <in.mmd> <out.png> 3`. The theme/spacing/curve config
  the init directive used to carry is injected by the renderer via
  `mermaid.initialize()`, so the PNG keeps the styled look. The script also works
  around two headless-shell quirks: mermaid's font CSS never reaches the
  foreignObject labels (falls back to Times), and the SVG collapses to container
  width unless pinned to its viewBox size.
- This folder (`docs/plans/`) is excluded from the public CZ-Sync/HISSS mirror,
  so the figure stays private until the manuscript ships.

```mermaid
flowchart TB

classDef raw fill:#dbeafe,stroke:#2563eb,stroke-width:1.5px,color:#14274e;
classDef proc fill:#cffafe,stroke:#0891b2,stroke-width:1.5px,color:#0b3b45;
classDef filt fill:#fef3c7,stroke:#d97706,stroke-width:1.5px,color:#5c2c06;
classDef inter fill:#ede9fe,stroke:#7c3aed,stroke-width:1.5px,color:#3b1d6e;
classDef out fill:#bbf7d0,stroke:#15803d,stroke-width:2.5px,color:#0c3a1e;
classDef flag fill:#fee2e2,stroke:#dc2626,stroke-width:1.5px,color:#6a1414;
classDef drop fill:#eef2f7,stroke:#7c8aa0,stroke-width:1px,color:#425064,stroke-dasharray:4 3;

subgraph SRC["RAW DATA SOURCES"]
  direction LR
  USGS[("USGS NWIS<br/>US daily discharge")]:::raw
  HYDAT[("HYDAT<br/>Canada daily discharge")]:::raw
  BND[("Basin boundaries<br/>GAGES-II · ECCC · HydroBASINS")]:::raw
  DAY[("Daymet ~1 km<br/>precip · temp · SWE<br/>area-weighted basin averages")]:::raw
  ATL[("HydroATLAS v10<br/>physiographic attributes")]:::raw
  LAI[("MODIS MCD15A3H<br/>LAI, 4-day")]:::raw
  LC[("MODIS MCD12Q1<br/>LULC, annual")]:::raw
  NLCD[("Annual NLCD<br/>land cover · impervious<br/>30 m, CONUS")]:::raw
end

BND --> DEL["Delineate contributing<br/>watershed per gage"]:::proc
DEL --> GEOM[/"Per-gage<br/>watershed geometry"/]:::inter

USGS --> CONV["Convert to mm/day ·<br/>assign water year Oct–Sep"]:::proc
HYDAT --> CONV
GEOM --> CONV
CONV --> DDP[/"Standardized<br/>daily streamflow"/]:::inter
DDP --> QC["Preprocess per gage:<br/>daily grid · interpolate gaps ≤3 d"]:::proc

QC -.-> FL{{"Flag only:<br/>negative Q ·<br/>constant monthly SD"}}:::flag
QC --> YR{"Reject year?<br/>over 3 consecutive NA days<br/>or over 30 NA days total"}:::filt
YR --> XY["Rejected:<br/>excluded years"]:::drop
YR --> WYS[/"Kept: clean<br/>water-year series"/]:::inter

WYS --> GG{"Include gage?<br/>≥20 qualifying years ·<br/>≥60% of window qualifying"}:::filt
GG --> XG["Excluded gages"]:::drop
GG --> SIG["Compute 121 signatures<br/>(100 annual series + 21 scalars)<br/>across 14 families"]:::proc

DAY --> SIG
SIG --> AN{"Area-normalized<br/>gage?"}:::filt
AN --> SKIP["Not normalized:<br/>skip Q-to-PPT signatures"]:::drop
AN --> ST["Summarize each annual series:<br/>mean · median · Theil-Sen ·<br/>linear · Spearman · Mann-Kendall ·<br/>Pettitt changepoint"]:::proc
SKIP --> ST

ST -.-> QAF{{"12 QA/QC range and<br/>consistency flags"}}:::flag
ST --> TG{"Trend gate:<br/>≥20 annual values · ≥60% of series ·<br/>≥80% of first and last decade"}:::filt

GEOM --> ZON["Zonal / area-weighted<br/>aggregation to watershed"]:::proc
ATL --> ZON
LAI --> ZON
LC --> ZON
NLCD --> ZON

TG --> O1[["Signature summary table<br/>+ annual-values parquet"]]:::out
ZON --> O2[["Watershed attribute table"]]:::out
ZON --> O3[["MODIS LAI monthly ·<br/>MODIS + NLCD land cover annual"]]:::out
O1 -.-> LIB[["Reproducible via open library<br/>Julia · Python · R"]]:::out
```

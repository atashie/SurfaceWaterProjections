# HISSS Dataset Workflow Schematic

Workflow figure for the HISSS manuscript (§2 Methods) — from raw sources through
QA screening, signature extraction, and watershed aggregation to the five
HydroShare resources' core tables. GitHub renders the diagram below natively.

- **Corrected 2026-09-01** (from the draft schematic): 121 signatures
  (100 annual series + 21 scalars) across **14** families (was "~100 across 13");
  trend gate stated as ≥20 values / **≥60% of series** / ≥80% of first and last
  decade (was "≥80% of series and each decade"); **Annual NLCD** added as a raw
  source and to the land-cover output table; Daymet edge labeled as pre-computed
  area-weighted basin averages; "gauge" → "gage" (USGS convention); the open
  library is linked as "reproducible via" rather than a downstream consumer.
- **PNG render**: `dataset_workflow_schematic.png` (same folder, 3× resolution,
  white ground) — regenerate with the Playwright/Chromium mermaid renderer
  (scratchpad tool `render_mermaid.py`; any mermaid-cli works too).
- This folder (`docs/plans/`) is excluded from the public CZ-Sync/HISSS mirror,
  so the figure stays private until the manuscript ships.

```mermaid
%%{init: {"theme":"base","themeVariables":{"fontFamily":"ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, sans-serif","fontSize":"13px","lineColor":"#7f8ea4"},"flowchart":{"htmlLabels":true,"nodeSpacing":40,"rankSpacing":52,"curve":"basis"}}}%%
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
  DAY[("Daymet ~1 km<br/>precip · temp · SWE")]:::raw
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
YR -->|reject| XY["Excluded years"]:::drop
YR -->|keep| WYS[/"Clean<br/>water-year series"/]:::inter

WYS --> GG{"Include gage?<br/>≥20 qualifying years ·<br/>≥60% of window qualifying"}:::filt
GG -->|no| XG["Excluded gages"]:::drop
GG -->|yes| SIG["Compute 121 signatures<br/>(100 annual series + 21 scalars)<br/>across 14 families"]:::proc

DAY -->|"area-weighted<br/>basin averages"| SIG
SIG --> AN{"Area-normalized<br/>gage?"}:::filt
AN -->|"no"| SKIP["Skip Q-to-PPT<br/>signatures"]:::drop
AN -->|"yes"| ST["Summarize each annual series:<br/>mean · median · Theil-Sen ·<br/>linear · Spearman · Mann-Kendall ·<br/>Pettitt changepoint"]:::proc
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
O1 -. reproducible via .-> LIB[["Open library<br/>Julia · Python · R"]]:::out
```

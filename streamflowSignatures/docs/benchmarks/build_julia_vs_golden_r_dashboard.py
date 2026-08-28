"""
Build an interactive HTML dashboard comparing Julia (post-Section 3)
output against the Golden R reference.

Usage:
    python docs/benchmarks/build_julia_vs_golden_r_dashboard.py

Output:
    docs/benchmarks/julia_vs_golden_r_dashboard.html
"""
import pandas as pd
import numpy as np
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

GOLDEN_R_PATH = os.path.join(PROJECT_ROOT, "golden-outputs",
                              "streamflow_signatures_full_10feb2026.csv")
JULIA_PATH = os.path.join(SCRIPT_DIR, "julia_signatures.csv")
OUTPUT_HTML = os.path.join(SCRIPT_DIR, "julia_vs_golden_r_dashboard.html")

SIGNATURE_GROUPS = {
    "Flow Volumes": ["Qann", "Qwin", "Qspr", "Qsum", "Qfal"],
    "Flow Percentiles": [
        "Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40", "Q50",
        "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99", "Q95_Q10",
    ],
    "Flow Timing": [
        "D1_day", "D5_day", "D10_day", "D20_day", "D30_day", "D40_day",
        "D50_day", "D60_day", "D70_day", "D80_day", "D90_day", "D95_day",
        "D99_day", "D25_to_D75", "Dmax",
    ],
    "Baseflow": ["BFI_Eckhardt", "BFI_LyneHollick"],
    "Flashiness": ["flashinessRB"],
    "Pulses": ["TQmean", "n_low_pulses_all", "dur_low_pulses_all",
               "n_high_pulses_all", "dur_high_pulses_all",
               "Flow_Reversals_annual", "Flow_Reversals_winter",
               "Flow_Reversals_spring", "Flow_Reversals_summer",
               "Flow_Reversals_fall"],
    "FDC": ["FDCall", "FDC90th", "FDCmid"],
    "Recession": ["log_a_pointcloud", "log_a_events", "b_pointcloud",
                  "b_events", "concavity", "n_recession_events"],
    "Runoff Ratios": ["annual_runoff_ratio", "winter_runoff_ratio",
                      "spring_runoff_ratio", "summer_runoff_ratio",
                      "fall_runoff_ratio"],
    "Elasticity": ["elasticity_rolling", "elasticity_annual"],
    "Q-P Seasonality": ["qp_slope_sd", "qp_bimodality"],
    "Storage": ["avg_storage"],
    "Negative Flow": ["negative_ann"],
}

STATS = ["_mean", "_median", "_senn_slp", "_linear_slp",
         "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval"]

SINGLE_VALUE_SIGS = [
    "elasticity_static",
    "runoff_ratio_high_count",
    "elasticity_years_total",
    "elasticity_years_low_ppt",
    "log_a_seasonality_minimum_all",
    "log_a_seasonality_minimum_first_half",
    "log_a_seasonality_minimum_last_half",
    "log_a_seasonality_amplitude_all",
    "log_a_seasonality_amplitude_first_half",
    "log_a_seasonality_amplitude_last_half",
]

STAT_SUFFIXES = ["_senn_slp", "_linear_slp", "_spearman_rho", "_spearman_pval",
                 "_mk_rho", "_mk_pval", "_mean", "_median"]


def normalize_col(col):
    """Normalize R column names to Julia convention."""
    col = col.replace("-", "_")
    col = col.replace(".", "_")
    col = col.replace("FDC_all", "FDCall")
    col = col.replace("FDC_90th", "FDC90th")
    col = col.replace("FDC_mid", "FDCmid")
    return col


# R's "elasticity_*" was renamed to "elasticity_rolling_*" in Julia
RENAMED_R_TO_JL = {}
for suf in STAT_SUFFIXES:
    RENAMED_R_TO_JL[f"elasticity{suf}"] = f"elasticity_rolling{suf}"


def build_data():
    print("Loading Golden R output...")
    r_df = pd.read_csv(GOLDEN_R_PATH, low_memory=False)
    r_df["gage_id"] = r_df["gage_id"].astype(str).str.strip()
    print(f"  {r_df.shape[0]} gages x {r_df.shape[1]} cols")

    print("Loading Julia post-Section 3 output...")
    jl_df = pd.read_csv(JULIA_PATH, low_memory=False)
    jl_df["gage_id"] = jl_df["gage_id"].astype(str).str.strip()
    print(f"  {jl_df.shape[0]} gages x {jl_df.shape[1]} cols")

    common_gages = sorted(set(r_df["gage_id"]) & set(jl_df["gage_id"]))
    print(f"  Common gages: {len(common_gages)}")

    r = r_df[r_df["gage_id"].isin(common_gages)].set_index("gage_id").loc[common_gages]
    j = jl_df[jl_df["gage_id"].isin(common_gages)].set_index("gage_id").loc[common_gages]

    # Normalize R column names and apply renames
    r_rename_map = {}
    for c in r.columns:
        nc = normalize_col(c)
        if nc in RENAMED_R_TO_JL:
            r_rename_map[c] = RENAMED_R_TO_JL[nc]
        elif nc != c:
            r_rename_map[c] = nc
    r = r.rename(columns=r_rename_map)

    # Get coordinates from Julia output (already has metadata merged)
    # If Julia is missing coords, backfill from golden-R or metadata
    for df_src, df_tgt in [(r, j), (j, r)]:
        if "latitude" in df_tgt.columns:
            n_missing = df_tgt["latitude"].isna().sum()
            if n_missing > 0 and "latitude" in df_src.columns:
                mask = df_tgt["latitude"].isna() & df_src["latitude"].notna()
                df_tgt.loc[mask, "latitude"] = df_src.loc[mask, "latitude"]
                df_tgt.loc[mask, "longitude"] = df_src.loc[mask, "longitude"]

    # Use Julia coords as primary (post-metadata-fix)
    lat_col = j["latitude"] if "latitude" in j.columns else r.get("latitude")
    lon_col = j["longitude"] if "longitude" in j.columns else r.get("longitude")

    # Build target column list
    target_cols = []
    col_to_display = {}  # normalized col -> display name for dropdown

    for sigs in SIGNATURE_GROUPS.values():
        for sig in sigs:
            for stat in STATS:
                col = sig + stat
                if col in r.columns and col in j.columns:
                    target_cols.append(col)
                elif col in j.columns and col not in r.columns:
                    target_cols.append(col)  # Julia-only

    for col in SINGLE_VALUE_SIGS:
        jl_col = col
        if col in r.columns and jl_col in j.columns:
            target_cols.append(col)
        elif jl_col in j.columns:
            target_cols.append(jl_col)

    # Deduplicate
    target_cols = list(dict.fromkeys(target_cols))

    n_common = sum(1 for c in target_cols if c in r.columns and c in j.columns)
    n_jl_only = sum(1 for c in target_cols if c in j.columns and c not in r.columns)
    print(f"  Target columns: {len(target_cols)} ({n_common} common, {n_jl_only} Julia-only)")

    # Build JSON data
    data = {
        "g": list(common_gages),
        "lt": [round(float(v), 4) if pd.notna(v) else None
               for v in lat_col.values],
        "ln": [round(float(v), 4) if pd.notna(v) else None
               for v in lon_col.values],
        "r": {},
        "j": {},
    }

    for col in target_cols:
        if col in r.columns:
            rv = pd.to_numeric(r[col], errors="coerce").values
            data["r"][col] = [round(float(v), 6) if np.isfinite(v) else None for v in rv]
        else:
            data["r"][col] = [None] * len(common_gages)

        if col in j.columns:
            jv = pd.to_numeric(j[col], errors="coerce").values
            data["j"][col] = [round(float(v), 6) if np.isfinite(v) else None for v in jv]
        else:
            data["j"][col] = [None] * len(common_gages)

    return data, target_cols


def build_html(data, target_cols):
    data_json = json.dumps(data, separators=(",", ":"))
    print(f"  JSON payload size: {len(data_json) / 1024 / 1024:.1f} MB")

    # Build signature options HTML
    optgroups = []
    for group_name, sigs in SIGNATURE_GROUPS.items():
        opts = []
        for sig in sigs:
            if any(sig + s in target_cols for s in STATS):
                opts.append(f'<option value="{sig}">{sig}</option>')
        if opts:
            optgroups.append(
                f'<optgroup label="{group_name}">{"".join(opts)}</optgroup>'
            )

    sv_opts = []
    for sig in SINGLE_VALUE_SIGS:
        if sig in target_cols:
            sv_opts.append(f'<option value="{sig}">{sig}</option>')
    if sv_opts:
        optgroups.append(
            f'<optgroup label="Single-Value (Scalars)">{"".join(sv_opts)}</optgroup>'
        )

    sig_options_html = "\n".join(optgroups)

    html = HTML_TEMPLATE.replace("__DATA_JSON__", data_json)
    html = html.replace("__SIG_OPTIONS__", sig_options_html)

    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"  Written to: {OUTPUT_HTML}")
    print(f"  File size: {os.path.getsize(OUTPUT_HTML) / 1024 / 1024:.1f} MB")


HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Julia (Post-Section 3) vs Golden R Reference</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #1a1a2e; color: #e0e0e0; }
.toolbar { display: flex; align-items: center; gap: 16px; padding: 12px 20px; background: #16213e; border-bottom: 1px solid #0f3460; flex-wrap: wrap; }
.toolbar label { font-size: 13px; color: #a0a0c0; }
.toolbar select { padding: 6px 10px; border-radius: 4px; border: 1px solid #0f3460; background: #1a1a2e; color: #e0e0e0; font-size: 14px; cursor: pointer; }
.toolbar select:disabled { opacity: 0.4; cursor: default; }
#summary { font-size: 14px; padding: 8px 20px; background: #16213e; border-bottom: 1px solid #0f3460; line-height: 1.6; }
#summary .r2-val { font-weight: bold; font-size: 16px; }
.tier-perfect { color: #4ade80; }
.tier-good { color: #86efac; }
.tier-poor { color: #facc15; }
.tier-low { color: #fb923c; }
.tier-verylow { color: #f87171; }
.tier-extremelylow { color: #ef4444; font-weight: bold; }
#tier-bar { font-size: 12px; padding: 4px 20px; background: #0f3460; display: flex; gap: 14px; flex-wrap: wrap; }
#tier-bar span { white-space: nowrap; }
#selection-info { font-size: 13px; padding: 4px 20px; background: #0f3460; color: #facc15; display: none; }
.maps-row { display: flex; height: 40vh; min-height: 280px; }
.map-col { flex: 1; position: relative; }
.map-col .map-label { position: absolute; top: 10px; left: 50%; transform: translateX(-50%); z-index: 1000; background: rgba(22,33,62,0.9); color: #e0e0e0; padding: 4px 12px; border-radius: 4px; font-size: 13px; font-weight: 600; pointer-events: none; }
.map-container { width: 100%; height: 100%; }
.legend { position: absolute; bottom: 20px; right: 20px; z-index: 1000; background: rgba(22,33,62,0.92); padding: 8px 12px; border-radius: 6px; font-size: 11px; }
.legend .grad { width: 20px; height: 120px; margin: 4px 0; }
.legend .labels { display: flex; flex-direction: column; justify-content: space-between; height: 120px; margin-left: 4px; }
.legend-row { display: flex; align-items: stretch; }
#scatter-container { height: 50vh; min-height: 380px; padding: 0 10px 10px 10px; background: #1a1a2e; }
</style>
</head>
<body>

<div class="toolbar">
  <div>
    <label>Signature:</label>
    <select id="sig-select">
      __SIG_OPTIONS__
    </select>
  </div>
  <div>
    <label>Statistic:</label>
    <select id="stat-select">
      <option value="_mean">Mean</option>
      <option value="_median">Median</option>
      <option value="_senn_slp">Theil-Sen Slope</option>
      <option value="_linear_slp">Linear Slope</option>
      <option value="_spearman_rho">Spearman Rho</option>
      <option value="_spearman_pval">Spearman P-value</option>
      <option value="_mk_rho">Mann-Kendall Tau</option>
      <option value="_mk_pval">Mann-Kendall P-value</option>
    </select>
  </div>
</div>
<div id="summary"></div>
<div id="tier-bar"></div>
<div id="selection-info"></div>

<div class="maps-row">
  <div class="map-col">
    <div class="map-label">Golden R (Feb 2026)</div>
    <div id="map-r" class="map-container"></div>
  </div>
  <div class="map-col">
    <div class="map-label">Julia Post-Section 3 (Apr 2026)</div>
    <div id="map-julia" class="map-container"></div>
  </div>
</div>

<div id="scatter-container"></div>

<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<script>
const DATA = __DATA_JSON__;

const SINGLE_VALUE = new Set([
  "elasticity_static",
  "runoff_ratio_high_count",
  "elasticity_years_total",
  "elasticity_years_low_ppt",
  "log_a_seasonality_minimum_all", "log_a_seasonality_minimum_first_half",
  "log_a_seasonality_minimum_last_half", "log_a_seasonality_amplitude_all",
  "log_a_seasonality_amplitude_first_half", "log_a_seasonality_amplitude_last_half"
]);

const VIRIDIS = [
  [68,1,84],[72,35,116],[64,67,135],[33,144,140],
  [53,183,121],[143,215,68],[210,226,27],[253,231,37]
];
const PVAL_PALETTE = [
  [255,255,204],[255,237,160],[254,217,118],[254,178,76],
  [253,141,60],[252,78,42],[227,26,28],[177,0,38]
];

function interpColor(palette, t) {
  t = Math.max(0, Math.min(1, t));
  let idx = t * (palette.length - 1);
  let lo = Math.floor(idx), hi = Math.min(lo + 1, palette.length - 1);
  let f = idx - lo;
  let r = Math.round(palette[lo][0] + (palette[hi][0] - palette[lo][0]) * f);
  let g = Math.round(palette[lo][1] + (palette[hi][1] - palette[lo][1]) * f);
  let b = Math.round(palette[lo][2] + (palette[hi][2] - palette[lo][2]) * f);
  return `rgb(${r},${g},${b})`;
}

function r2Identity(rArr, jArr) {
  let pairs = [];
  for (let i = 0; i < rArr.length; i++) {
    if (rArr[i] !== null && jArr[i] !== null) pairs.push([rArr[i], jArr[i]]);
  }
  if (pairs.length < 10) return null;
  let meanY = pairs.reduce((s, p) => s + p[1], 0) / pairs.length;
  let ssRes = pairs.reduce((s, p) => s + (p[1] - p[0]) ** 2, 0);
  let ssTot = pairs.reduce((s, p) => s + (p[1] - meanY) ** 2, 0);
  if (ssTot === 0) return ssRes === 0 ? 1.0 : null;
  return 1.0 - ssRes / ssTot;
}

function classifyR2(r2) {
  if (r2 === null) return { tier: "N/A", cls: "" };
  if (r2 >= 0.999) return { tier: "Perfect", cls: "tier-perfect" };
  if (r2 >= 0.99) return { tier: "Good", cls: "tier-good" };
  if (r2 >= 0.95) return { tier: "Poor", cls: "tier-poor" };
  if (r2 >= 0.9) return { tier: "Low", cls: "tier-low" };
  if (r2 >= 0.5) return { tier: "Very Low", cls: "tier-verylow" };
  return { tier: "Extremely Low", cls: "tier-extremelylow" };
}

function percentile(arr, p) {
  let s = arr.slice().sort((a, b) => a - b);
  let idx = (p / 100) * (s.length - 1);
  let lo = Math.floor(idx), hi = Math.ceil(idx);
  return s[lo] + (s[hi] - s[lo]) * (idx - lo);
}

// Compute global tier counts across all columns for the tier bar
function computeGlobalTiers() {
  let tiers = { "Perfect": 0, "Good": 0, "Poor": 0, "Low": 0, "Very Low": 0, "Extremely Low": 0, "N/A": 0 };
  let allCols = Object.keys(DATA.r);
  for (let col of allCols) {
    let rArr = DATA.r[col], jArr = DATA.j[col];
    if (!rArr || !jArr) continue;
    // Skip Julia-only columns (all R values are null)
    if (rArr.every(v => v === null)) continue;
    let r2 = r2Identity(rArr, jArr);
    let t = classifyR2(r2).tier;
    tiers[t] = (tiers[t] || 0) + 1;
  }
  return tiers;
}

// Initialize maps
const mapR = L.map("map-r", { preferCanvas: true }).setView([41.3, -95.6], 4);
const mapJ = L.map("map-julia", { preferCanvas: true }).setView([41.3, -95.6], 4);

L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapR);
L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapJ);

let syncing = false;
function syncMap(source, target) {
  source.on("moveend zoomend", function () {
    if (syncing) return;
    syncing = true;
    target.setView(source.getCenter(), source.getZoom(), { animate: false });
    syncing = false;
  });
}
syncMap(mapR, mapJ);
syncMap(mapJ, mapR);

const markersR = [];
const markersJ = [];
const N = DATA.g.length;

for (let i = 0; i < N; i++) {
  if (DATA.lt[i] === null || DATA.ln[i] === null) continue;
  let ll = [DATA.lt[i], DATA.ln[i]];
  let mR = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapR);
  let mJ = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapJ);
  mR._idx = i;
  mJ._idx = i;
  markersR.push(mR);
  markersJ.push(mJ);
}

// Legend
const legendDiv = document.createElement("div");
legendDiv.className = "legend";
legendDiv.innerHTML = `
  <div style="font-weight:600;margin-bottom:4px;">Value</div>
  <div class="legend-row">
    <canvas id="legend-canvas" width="20" height="120"></canvas>
    <div class="labels">
      <span id="leg-max"></span>
      <span id="leg-mid"></span>
      <span id="leg-min"></span>
    </div>
  </div>`;
const LegendControl = L.Control.extend({
  onAdd: function() { return legendDiv; }
});
new LegendControl({ position: "bottomright" }).addTo(mapJ);

function drawLegend(vMin, vMax, isPval) {
  let canvas = document.getElementById("legend-canvas");
  let ctx = canvas.getContext("2d");
  let pal = isPval ? PVAL_PALETTE : VIRIDIS;
  for (let y = 0; y < 120; y++) {
    let t = 1 - y / 119;
    ctx.fillStyle = interpColor(pal, t);
    ctx.fillRect(0, y, 20, 1);
  }
  let fmt = (v) => Math.abs(v) >= 1000 ? v.toExponential(1) : v.toPrecision(4);
  document.getElementById("leg-max").textContent = fmt(vMax);
  document.getElementById("leg-mid").textContent = fmt((vMin + vMax) / 2);
  document.getElementById("leg-min").textContent = fmt(vMin);
}

const scatterDiv = document.getElementById("scatter-container");
const selInfo = document.getElementById("selection-info");
const tierBar = document.getElementById("tier-bar");

let highlightR = null;
let highlightJ = null;

const markerByIdxR = {};
const markerByIdxJ = {};
for (let m of markersR) markerByIdxR[m._idx] = m;
for (let m of markersJ) markerByIdxJ[m._idx] = m;

let scatterToData = [];

function clearHighlight() {
  if (highlightR) { mapR.removeLayer(highlightR); highlightR = null; }
  if (highlightJ) { mapJ.removeLayer(highlightJ); highlightJ = null; }
  selInfo.style.display = "none";
}

function highlightGage(dataIdx) {
  clearHighlight();
  let lat = DATA.lt[dataIdx], lon = DATA.ln[dataIdx];
  if (lat === null || lon === null) return;

  let ll = [lat, lon];
  highlightR = L.circleMarker(ll, {
    radius: 10, weight: 3, color: "#facc15", fillColor: "#facc15",
    fillOpacity: 0.3, pane: "markerPane"
  }).addTo(mapR);
  highlightJ = L.circleMarker(ll, {
    radius: 10, weight: 3, color: "#facc15", fillColor: "#facc15",
    fillOpacity: 0.3, pane: "markerPane"
  }).addTo(mapJ);

  syncing = true;
  mapR.setView(ll, 8, { animate: true });
  mapJ.setView(ll, 8, { animate: true });
  setTimeout(() => { syncing = false; }, 300);

  let mR = markerByIdxR[dataIdx], mJ = markerByIdxJ[dataIdx];
  if (mR) mR.openPopup();
  if (mJ) mJ.openPopup();

  let colKey = getCurrentColKey();
  let rV = DATA.r[colKey] ? DATA.r[colKey][dataIdx] : null;
  let jV = DATA.j[colKey] ? DATA.j[colKey][dataIdx] : null;
  let diff = (rV !== null && jV !== null) ? (jV - rV).toPrecision(4) : "N/A";
  selInfo.innerHTML =
    `Selected: <strong>${DATA.g[dataIdx]}</strong> &nbsp;|&nbsp; ` +
    `Golden R = ${rV !== null ? rV : "NA"} &nbsp;|&nbsp; ` +
    `Julia = ${jV !== null ? jV : "NA"} &nbsp;|&nbsp; ` +
    `Diff = ${diff} &nbsp; ` +
    `<span style="cursor:pointer;text-decoration:underline;" onclick="clearSelection()">[clear]</span>`;
  selInfo.style.display = "block";
}

function clearSelection() {
  clearHighlight();
  let colors = scatterToData.map(() => "rgba(99,179,237,0.5)");
  let sizes = scatterToData.map(() => 3);
  Plotly.restyle(scatterDiv, { "marker.color": [colors], "marker.size": [sizes] }, [0]);
}

function getCurrentColKey() {
  let sig = document.getElementById("sig-select").value;
  let stat = document.getElementById("stat-select").value;
  return SINGLE_VALUE.has(sig) ? sig : sig + stat;
}

// Show global tier bar
function showGlobalTiers() {
  let tiers = computeGlobalTiers();
  let total = Object.values(tiers).reduce((a, b) => a + b, 0) - (tiers["N/A"] || 0);
  let parts = [
    ["Perfect", "tier-perfect", ">= 0.999"],
    ["Good", "tier-good", ">= 0.99"],
    ["Poor", "tier-poor", ">= 0.95"],
    ["Low", "tier-low", ">= 0.90"],
    ["Very Low", "tier-verylow", ">= 0.50"],
    ["Extremely Low", "tier-extremelylow", "< 0.50"],
  ];
  let html = '<span style="color:#a0a0c0;">All columns:</span> ';
  for (let [name, cls, desc] of parts) {
    let c = tiers[name] || 0;
    let pct = total > 0 ? (100 * c / total).toFixed(1) : "0";
    html += `<span class="${cls}">${name}: ${c} (${pct}%)</span> `;
  }
  tierBar.innerHTML = html;
}

function update() {
  clearHighlight();
  let sig = document.getElementById("sig-select").value;
  let stat = document.getElementById("stat-select").value;
  let statSel = document.getElementById("stat-select");
  let isSingle = SINGLE_VALUE.has(sig);
  statSel.disabled = isSingle;
  let colKey = isSingle ? sig : sig + stat;

  let rVals = DATA.r[colKey];
  let jVals = DATA.j[colKey];
  if (!rVals || !jVals) {
    document.getElementById("summary").innerHTML = `<em>Column "${colKey}" not found.</em>`;
    return;
  }

  let nR_NA = rVals.filter(v => v === null).length;
  let nJ_NA = jVals.filter(v => v === null).length;
  let r2 = r2Identity(rVals, jVals);
  let nPaired = 0;
  for (let i = 0; i < N; i++) if (rVals[i] !== null && jVals[i] !== null) nPaired++;

  let r2Str = r2 !== null ? r2.toFixed(4) : "N/A";
  let { tier, cls } = classifyR2(r2);
  let jlOnlyNote = nR_NA === N ? " <em>(NEW column -- no Golden R data)</em>" : "";
  document.getElementById("summary").innerHTML =
    `<strong>${colKey}</strong> &mdash; ` +
    `R&sup2; = <span class="r2-val ${cls}">${r2Str}</span> ` +
    `<span class="${cls}">[${tier}]</span> &nbsp;|&nbsp; ` +
    `N paired = ${nPaired} &nbsp;|&nbsp; ` +
    `Golden R NA = ${nR_NA}, Julia NA = ${nJ_NA}` + jlOnlyNote;

  // Color scale
  let allVals = [];
  for (let i = 0; i < N; i++) {
    if (rVals[i] !== null) allVals.push(rVals[i]);
    if (jVals[i] !== null) allVals.push(jVals[i]);
  }
  let isPval = colKey.endsWith("_pval");
  let vMin, vMax;
  if (isPval) {
    vMin = 0; vMax = 1;
  } else if (allVals.length > 10) {
    vMin = percentile(allVals, 2);
    vMax = percentile(allVals, 98);
  } else {
    vMin = Math.min(...allVals);
    vMax = Math.max(...allVals);
  }
  if (vMax === vMin) vMax = vMin + 1;

  let pal = isPval ? PVAL_PALETTE : VIRIDIS;
  function colorFor(v) {
    if (v === null) return "#444";
    let t = (v - vMin) / (vMax - vMin);
    return interpColor(pal, t);
  }

  for (let m of markersR) {
    let v = rVals[m._idx];
    m.setStyle({ fillColor: colorFor(v) });
    m.unbindPopup();
    m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>Golden R: ${v !== null ? v : "NA"}`);
  }
  for (let m of markersJ) {
    let v = jVals[m._idx];
    m.setStyle({ fillColor: colorFor(v) });
    m.unbindPopup();
    m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>Julia: ${v !== null ? v : "NA"}`);
  }

  drawLegend(vMin, vMax, isPval);

  // Scatterplot
  let xVals = [], yVals = [], hoverText = [];
  scatterToData = [];
  for (let i = 0; i < N; i++) {
    if (rVals[i] !== null && jVals[i] !== null) {
      xVals.push(rVals[i]);
      yVals.push(jVals[i]);
      hoverText.push(DATA.g[i]);
      scatterToData.push(i);
    }
  }

  let pLo, pHi;
  if (isPval) {
    pLo = 0; pHi = 1;
  } else if (xVals.length > 20) {
    let combined = xVals.concat(yVals);
    pLo = percentile(combined, 1);
    pHi = percentile(combined, 99);
    let pad = (pHi - pLo) * 0.05;
    pLo -= pad; pHi += pad;
  } else {
    pLo = Math.min(...xVals, ...yVals);
    pHi = Math.max(...xVals, ...yVals);
  }

  Plotly.react(scatterDiv, [{
    x: xVals, y: yVals, text: hoverText,
    mode: "markers",
    type: "scattergl",
    marker: { size: 3, color: "rgba(99,179,237,0.5)" },
    hovertemplate: "<b>%{text}</b><br>Golden R: %{x:.4f}<br>Julia: %{y:.4f}<extra></extra>"
  }], {
    shapes: [{
      type: "line", x0: pLo, y0: pLo, x1: pHi, y1: pHi,
      line: { color: "rgba(250,200,50,0.5)", width: 1.5, dash: "dash" }
    }],
    xaxis: { title: "Golden R (Feb 2026)", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: "Julia Post-Section 3 (Apr 2026)", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `${colKey}  (R\u00b2 = ${r2Str},  N = ${nPaired},  Tier: ${tier})`, font: { color: "#e0e0e0", size: 14 } },
    plot_bgcolor: "#1a1a2e",
    paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 50, l: 60, r: 20 },
    font: { color: "#a0a0c0" }
  }, { responsive: true });
}

document.getElementById("sig-select").addEventListener("change", update);
document.getElementById("stat-select").addEventListener("change", update);

showGlobalTiers();
update();

scatterDiv.on("plotly_click", function(eventData) {
  if (!eventData || !eventData.points || eventData.points.length === 0) return;
  let ptIdx = eventData.points[0].pointIndex;
  let dataIdx = scatterToData[ptIdx];
  if (dataIdx === undefined) return;

  let colors = scatterToData.map(() => "rgba(99,179,237,0.5)");
  let sizes = scatterToData.map(() => 3);
  colors[ptIdx] = "rgba(250,200,50,1)";
  sizes[ptIdx] = 8;
  Plotly.restyle(scatterDiv, { "marker.color": [colors], "marker.size": [sizes] }, [0]);

  highlightGage(dataIdx);
});
</script>
</body>
</html>
'''

if __name__ == "__main__":
    data, target_cols = build_data()
    build_html(data, target_cols)
    print("Done.")

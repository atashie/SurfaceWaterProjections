"""
Build an interactive HTML dashboard comparing Python or rpkg output
against the Golden Julia (canonical) reference.

Usage:
    python docs/benchmarks/build_julia_golden_dashboard.py python
    python docs/benchmarks/build_julia_golden_dashboard.py rpkg

Output:
    docs/benchmarks/python_vs_golden_julia_dashboard.html
    docs/benchmarks/rpkg_vs_golden_julia_dashboard.html
"""
import pandas as pd
import numpy as np
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

GOLDEN_JULIA_PATH = os.path.join(PROJECT_ROOT, "golden-outputs",
                                  "streamflow_signatures_julia_apr2026.csv")

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
    "Baseflow": ["BFI_Eckhardt", "BFI_LyneHollick",
                 "BFI_Eckhardt_param", "BFI_LyneHollick_param"],
    "Flashiness": ["flashinessRB"],
    "Pulses": ["TQmean", "n_low_pulses_all", "dur_low_pulses_all",
               "n_high_pulses_all", "dur_high_pulses_all",
               "n_low_pulses_year", "dur_low_pulses_year",
               "n_high_pulses_year", "dur_high_pulses_year",
               "Flow_Reversals_annual", "Flow_Reversals_winter",
               "Flow_Reversals_spring", "Flow_Reversals_summer",
               "Flow_Reversals_fall", "negative_ann"],
    "FDC": ["FDCall", "FDC90th", "FDCmid"],
    "Recession": ["log_a_pointcloud", "log_a_events", "b_pointcloud",
                  "b_events", "concavity", "n_recession_events", "alpha_linear"],
    "Runoff Ratios": ["annual_runoff_ratio", "winter_runoff_ratio",
                      "spring_runoff_ratio", "summer_runoff_ratio",
                      "fall_runoff_ratio"],
    "Elasticity": ["elasticity_rolling", "elasticity_annual"],
    "Q-P Seasonality": ["qp_slope_sd", "qp_bimodality"],
    "Storage": ["avg_storage"],
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
    "recession_alpha_point_cloud_linear_reservoir",
    "season_excluded_years_winter",
    "season_excluded_years_spring",
    "season_excluded_years_summer",
    "season_excluded_years_fall",
]


def normalize_col(col):
    col = col.replace("-", "_")
    col = col.replace(".", "_")
    col = col.replace("FDC_all", "FDCall")
    col = col.replace("FDC_90th", "FDC90th")
    col = col.replace("FDC_mid", "FDCmid")
    return col


def build_data(impl_name, impl_path):
    print("Loading Golden Julia output...")
    jl_df = pd.read_csv(GOLDEN_JULIA_PATH, low_memory=False)
    jl_df["gage_id"] = jl_df["gage_id"].astype(str).str.strip()
    print(f"  {jl_df.shape[0]} gages x {jl_df.shape[1]} cols")

    print(f"Loading {impl_name} benchmark output...")
    impl_df = pd.read_csv(impl_path, low_memory=False)
    impl_df["gage_id"] = impl_df["gage_id"].astype(str).str.strip()
    print(f"  {impl_df.shape[0]} gages x {impl_df.shape[1]} cols")

    common_gages = sorted(set(jl_df["gage_id"]) & set(impl_df["gage_id"]))
    print(f"  Common gages: {len(common_gages)}")

    jl = jl_df[jl_df["gage_id"].isin(common_gages)].set_index("gage_id").loc[common_gages]
    impl = impl_df[impl_df["gage_id"].isin(common_gages)].set_index("gage_id").loc[common_gages]

    # Normalize impl column names
    impl_rename = {}
    for c in impl.columns:
        nc = normalize_col(c)
        if nc != c:
            impl_rename[c] = nc
    impl = impl.rename(columns=impl_rename)

    # Coordinates from Julia (canonical has metadata)
    lat_col = jl["latitude"] if "latitude" in jl.columns else impl.get("latitude")
    lon_col = jl["longitude"] if "longitude" in jl.columns else impl.get("longitude")

    # Build target column list
    target_cols = []
    for sigs in SIGNATURE_GROUPS.values():
        for sig in sigs:
            for stat in STATS:
                col = sig + stat
                if col in jl.columns and col in impl.columns:
                    target_cols.append(col)

    for col in SINGLE_VALUE_SIGS:
        if col in jl.columns and col in impl.columns:
            target_cols.append(col)

    target_cols = list(dict.fromkeys(target_cols))

    n_common = sum(1 for c in target_cols if c in jl.columns and c in impl.columns)
    n_jl_only = sum(1 for c in target_cols if c in jl.columns and c not in impl.columns)
    print(f"  Target columns: {len(target_cols)} ({n_common} common, {n_jl_only} Julia-only)")

    data = {
        "g": list(common_gages),
        "lt": [round(float(v), 4) if pd.notna(v) else None for v in lat_col.values],
        "ln": [round(float(v), 4) if pd.notna(v) else None for v in lon_col.values],
        "jl": {},
        "impl": {},
    }

    for col in target_cols:
        if col in jl.columns:
            jv = pd.to_numeric(jl[col], errors="coerce").values
            data["jl"][col] = [round(float(v), 6) if np.isfinite(v) else None for v in jv]
        else:
            data["jl"][col] = [None] * len(common_gages)

        if col in impl.columns:
            iv = pd.to_numeric(impl[col], errors="coerce").values
            data["impl"][col] = [round(float(v), 6) if np.isfinite(v) else None for v in iv]
        else:
            data["impl"][col] = [None] * len(common_gages)

    return data, target_cols


def build_html(data, target_cols, impl_name, output_path):
    data_json = json.dumps(data, separators=(",", ":"))
    print(f"  JSON payload size: {len(data_json) / 1024 / 1024:.1f} MB")

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
    html = html.replace("__IMPL_NAME__", impl_name)
    html = html.replace("__IMPL_UPPER__", impl_name.upper())

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"  Written to: {output_path}")
    print(f"  File size: {os.path.getsize(output_path) / 1024 / 1024:.1f} MB")


HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>__IMPL_UPPER__ vs Golden Julia (Canonical)</title>
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
    <div class="map-label">Golden Julia (Canonical)</div>
    <div id="map-jl" class="map-container"></div>
  </div>
  <div class="map-col">
    <div class="map-label">__IMPL_NAME__</div>
    <div id="map-impl" class="map-container"></div>
  </div>
</div>

<div id="scatter-container"></div>

<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<script>
const DATA = __DATA_JSON__;
const IMPL_NAME = "__IMPL_NAME__";

// Keep in sync with SINGLE_VALUE_SIGS in this file's Python section.
const SINGLE_VALUE = new Set([
  "elasticity_static", "runoff_ratio_high_count",
  "elasticity_years_total", "elasticity_years_low_ppt",
  "log_a_seasonality_minimum_all", "log_a_seasonality_minimum_first_half",
  "log_a_seasonality_minimum_last_half", "log_a_seasonality_amplitude_all",
  "log_a_seasonality_amplitude_first_half", "log_a_seasonality_amplitude_last_half",
  "recession_alpha_point_cloud_linear_reservoir",
  "season_excluded_years_winter", "season_excluded_years_spring",
  "season_excluded_years_summer", "season_excluded_years_fall"
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

function r2Identity(a, b) {
  let pairs = [];
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== null && b[i] !== null) pairs.push([a[i], b[i]]);
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

function computeGlobalTiers() {
  let tiers = { "Perfect": 0, "Good": 0, "Poor": 0, "Low": 0, "Very Low": 0, "Extremely Low": 0, "N/A": 0 };
  let allCols = Object.keys(DATA.jl);
  for (let col of allCols) {
    let jlArr = DATA.jl[col], implArr = DATA.impl[col];
    if (!jlArr || !implArr) continue;
    if (implArr.every(v => v === null)) continue;
    let r2 = r2Identity(jlArr, implArr);
    let t = classifyR2(r2).tier;
    tiers[t] = (tiers[t] || 0) + 1;
  }
  return tiers;
}

const mapJl = L.map("map-jl", { preferCanvas: true }).setView([41.3, -95.6], 4);
const mapImpl = L.map("map-impl", { preferCanvas: true }).setView([41.3, -95.6], 4);

L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapJl);
L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapImpl);

let syncing = false;
function syncMap(source, target) {
  source.on("moveend zoomend", function () {
    if (syncing) return;
    syncing = true;
    target.setView(source.getCenter(), source.getZoom(), { animate: false });
    syncing = false;
  });
}
syncMap(mapJl, mapImpl);
syncMap(mapImpl, mapJl);

const markersJl = [];
const markersImpl = [];
const N = DATA.g.length;

for (let i = 0; i < N; i++) {
  if (DATA.lt[i] === null || DATA.ln[i] === null) continue;
  let ll = [DATA.lt[i], DATA.ln[i]];
  let mJl = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapJl);
  let mImpl = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapImpl);
  mJl._idx = i;
  mImpl._idx = i;
  markersJl.push(mJl);
  markersImpl.push(mImpl);
}

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
const LegendControl = L.Control.extend({ onAdd: function() { return legendDiv; } });
new LegendControl({ position: "bottomright" }).addTo(mapImpl);

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

let highlightJl = null, highlightImpl = null;
const markerByIdxJl = {}, markerByIdxImpl = {};
for (let m of markersJl) markerByIdxJl[m._idx] = m;
for (let m of markersImpl) markerByIdxImpl[m._idx] = m;

let scatterToData = [];

function clearHighlight() {
  if (highlightJl) { mapJl.removeLayer(highlightJl); highlightJl = null; }
  if (highlightImpl) { mapImpl.removeLayer(highlightImpl); highlightImpl = null; }
  selInfo.style.display = "none";
}

function highlightGage(dataIdx) {
  clearHighlight();
  let lat = DATA.lt[dataIdx], lon = DATA.ln[dataIdx];
  if (lat === null || lon === null) return;
  let ll = [lat, lon];
  highlightJl = L.circleMarker(ll, { radius: 10, weight: 3, color: "#facc15", fillColor: "#facc15", fillOpacity: 0.3, pane: "markerPane" }).addTo(mapJl);
  highlightImpl = L.circleMarker(ll, { radius: 10, weight: 3, color: "#facc15", fillColor: "#facc15", fillOpacity: 0.3, pane: "markerPane" }).addTo(mapImpl);
  syncing = true;
  mapJl.setView(ll, 8, { animate: true });
  mapImpl.setView(ll, 8, { animate: true });
  setTimeout(() => { syncing = false; }, 300);

  let colKey = getCurrentColKey();
  let jlV = DATA.jl[colKey] ? DATA.jl[colKey][dataIdx] : null;
  let implV = DATA.impl[colKey] ? DATA.impl[colKey][dataIdx] : null;
  let diff = (jlV !== null && implV !== null) ? (implV - jlV).toPrecision(4) : "N/A";
  selInfo.innerHTML =
    `Selected: <strong>${DATA.g[dataIdx]}</strong> &nbsp;|&nbsp; ` +
    `Golden Julia = ${jlV !== null ? jlV : "NA"} &nbsp;|&nbsp; ` +
    `${IMPL_NAME} = ${implV !== null ? implV : "NA"} &nbsp;|&nbsp; ` +
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

function showGlobalTiers() {
  let tiers = computeGlobalTiers();
  let total = Object.values(tiers).reduce((a, b) => a + b, 0) - (tiers["N/A"] || 0);
  let parts = [
    ["Perfect", "tier-perfect"], ["Good", "tier-good"], ["Poor", "tier-poor"],
    ["Low", "tier-low"], ["Very Low", "tier-verylow"], ["Extremely Low", "tier-extremelylow"],
  ];
  let html = '<span style="color:#a0a0c0;">All columns:</span> ';
  for (let [name, cls] of parts) {
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

  let jlVals = DATA.jl[colKey];
  let implVals = DATA.impl[colKey];
  if (!jlVals || !implVals) {
    document.getElementById("summary").innerHTML = `<em>Column "${colKey}" not found.</em>`;
    return;
  }

  let nJl_NA = jlVals.filter(v => v === null).length;
  let nImpl_NA = implVals.filter(v => v === null).length;
  let r2 = r2Identity(jlVals, implVals);
  let nPaired = 0;
  for (let i = 0; i < N; i++) if (jlVals[i] !== null && implVals[i] !== null) nPaired++;

  let r2Str = r2 !== null ? r2.toFixed(4) : "N/A";
  let { tier, cls } = classifyR2(r2);
  document.getElementById("summary").innerHTML =
    `<strong>${colKey}</strong> &mdash; ` +
    `R&sup2; = <span class="r2-val ${cls}">${r2Str}</span> ` +
    `<span class="${cls}">[${tier}]</span> &nbsp;|&nbsp; ` +
    `N paired = ${nPaired} &nbsp;|&nbsp; ` +
    `Golden Julia NA = ${nJl_NA}, ${IMPL_NAME} NA = ${nImpl_NA}`;

  let allVals = [];
  for (let i = 0; i < N; i++) {
    if (jlVals[i] !== null) allVals.push(jlVals[i]);
    if (implVals[i] !== null) allVals.push(implVals[i]);
  }
  let isPval = colKey.endsWith("_pval");
  let vMin, vMax;
  if (isPval) { vMin = 0; vMax = 1; }
  else if (allVals.length > 10) { vMin = percentile(allVals, 2); vMax = percentile(allVals, 98); }
  else { vMin = Math.min(...allVals); vMax = Math.max(...allVals); }
  if (vMax === vMin) vMax = vMin + 1;

  let pal = isPval ? PVAL_PALETTE : VIRIDIS;
  function colorFor(v) { return v === null ? "#444" : interpColor(pal, (v - vMin) / (vMax - vMin)); }

  for (let m of markersJl) {
    let v = jlVals[m._idx];
    m.setStyle({ fillColor: colorFor(v) });
    m.unbindPopup();
    m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>Golden Julia: ${v !== null ? v : "NA"}`);
  }
  for (let m of markersImpl) {
    let v = implVals[m._idx];
    m.setStyle({ fillColor: colorFor(v) });
    m.unbindPopup();
    m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>${IMPL_NAME}: ${v !== null ? v : "NA"}`);
  }

  drawLegend(vMin, vMax, isPval);

  let xVals = [], yVals = [], hoverText = [];
  scatterToData = [];
  for (let i = 0; i < N; i++) {
    if (jlVals[i] !== null && implVals[i] !== null) {
      xVals.push(jlVals[i]);
      yVals.push(implVals[i]);
      hoverText.push(DATA.g[i]);
      scatterToData.push(i);
    }
  }

  if (xVals.length === 0) {
    Plotly.react(scatterDiv, [], {
      title: { text: `${colKey}  (no paired data)`, font: { color: "#e0e0e0", size: 14 } },
      plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
      margin: { t: 40, b: 50, l: 60, r: 20 }, font: { color: "#a0a0c0" }
    }, { responsive: true });
    return;
  }

  let pLo, pHi;
  if (isPval) { pLo = 0; pHi = 1; }
  else if (xVals.length > 20) {
    let combined = xVals.concat(yVals);
    pLo = percentile(combined, 1);
    pHi = percentile(combined, 99);
    let pad = (pHi - pLo) * 0.05;
    pLo -= pad; pHi += pad;
  } else { pLo = Math.min(...xVals, ...yVals); pHi = Math.max(...xVals, ...yVals); }

  Plotly.react(scatterDiv, [{
    x: xVals, y: yVals, text: hoverText,
    mode: "markers", type: "scattergl",
    marker: { size: 3, color: "rgba(99,179,237,0.5)" },
    hovertemplate: "<b>%{text}</b><br>Golden Julia: %{x:.4f}<br>" + IMPL_NAME + ": %{y:.4f}<extra></extra>"
  }], {
    shapes: [{ type: "line", x0: pLo, y0: pLo, x1: pHi, y1: pHi,
      line: { color: "rgba(250,200,50,0.5)", width: 1.5, dash: "dash" } }],
    xaxis: { title: "Golden Julia (Canonical)", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: IMPL_NAME, range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `${colKey}  (R\u00b2 = ${r2Str},  N = ${nPaired},  Tier: ${tier})`, font: { color: "#e0e0e0", size: 14 } },
    plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
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
    if len(sys.argv) < 2 or sys.argv[1] not in ("python", "rpkg"):
        print("Usage: python build_julia_golden_dashboard.py <python|rpkg>")
        sys.exit(1)

    impl = sys.argv[1]
    if impl == "python":
        impl_path = os.path.join(SCRIPT_DIR, "python_signatures.csv")
        output_path = os.path.join(SCRIPT_DIR, "python_vs_golden_julia_dashboard.html")
        impl_name = "Python"
    else:
        impl_path = os.path.join(SCRIPT_DIR, "rpkg_signatures.csv")
        output_path = os.path.join(SCRIPT_DIR, "rpkg_vs_golden_julia_dashboard.html")
        impl_name = "rpkg"

    if not os.path.exists(GOLDEN_JULIA_PATH):
        print(f"ERROR: Golden Julia file not found: {GOLDEN_JULIA_PATH}")
        sys.exit(1)
    if not os.path.exists(impl_path):
        print(f"ERROR: {impl_name} benchmark output not found: {impl_path}")
        sys.exit(1)

    data, target_cols = build_data(impl_name, impl_path)
    build_html(data, target_cols, impl_name, output_path)
    print("Done.")

"""
Build an interactive HTML dashboard comparing Julia benchmark output
against Golden R reference outputs.

Usage:
    python docs/benchmarks/build_comparison_dashboard.py

Output:
    docs/benchmarks/signature_comparison.html
"""
import pandas as pd
import numpy as np
import json
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

GOLDEN_PATH = os.path.join(PROJECT_ROOT, "golden-outputs", "streamflow_signatures_full_10feb2026.csv")
JULIA_PATH = os.path.join(SCRIPT_DIR, "julia_signatures.csv")
OUTPUT_PATH = os.path.join(SCRIPT_DIR, "signature_comparison.html")

# Signatures to include, grouped by category
SIGNATURE_GROUPS = {
    "Flow Volumes": ["Qann", "Qwin", "Qspr", "Qsum", "Qfal"],
    "Flow Timing": [
        "D5_day", "D10_day", "D20_day", "D30_day", "D40_day",
        "D50_day", "D60_day", "D70_day", "D80_day", "D90_day",
        "D95_day", "D25_to_D75", "Dmax",
    ],
    "Baseflow": ["BFI_Eckhardt", "BFI_LyneHollick"],
    "Flashiness": ["flashinessRB"],
    "Pulses": ["TQmean", "n_low_pulses_all", "dur_low_pulses_all",
               "n_high_pulses_all", "dur_high_pulses_all"],
    "FDC": ["FDCall", "FDC90th", "FDCmid"],
    "Runoff Ratios": ["annual_runoff_ratio", "winter_runoff_ratio",
                      "spring_runoff_ratio", "summer_runoff_ratio",
                      "fall_runoff_ratio"],
    "Elasticity": ["elasticity"],
    "Q-P Seasonality": ["qp_slope_sd", "qp_bimodality"],
    "Storage": ["avg_storage"],
}

STATS = ["_mean", "_median", "_senn_slp", "_mk_pval"]

SINGLE_VALUE_SIGS = [
    "elasticity_static",
    "log_a_seasonality_minimum_all",
    "log_a_seasonality_minimum_first_half",
    "log_a_seasonality_minimum_last_half",
    "log_a_seasonality_amplitude_all",
    "log_a_seasonality_amplitude_first_half",
    "log_a_seasonality_amplitude_last_half",
]

def build_data():
    print("Loading golden R output...")
    golden = pd.read_csv(GOLDEN_PATH, low_memory=False)
    print(f"  {golden.shape[0]} gages x {golden.shape[1]} cols")

    print("Loading Julia output...")
    julia = pd.read_csv(JULIA_PATH, low_memory=False)
    print(f"  {julia.shape[0]} gages x {julia.shape[1]} cols")

    golden["gage_id"] = golden["gage_id"].astype(str).str.strip()
    julia["gage_id"] = julia["gage_id"].astype(str).str.strip()

    common_gages = sorted(set(golden["gage_id"]) & set(julia["gage_id"]))
    print(f"  Common gages: {len(common_gages)}")

    g = golden.set_index("gage_id").loc[common_gages]
    j = julia.set_index("gage_id").loc[common_gages]

    # Build target column list
    target_cols = []
    for sigs in SIGNATURE_GROUPS.values():
        for sig in sigs:
            for stat in STATS:
                col = sig + stat
                if col in g.columns and col in j.columns:
                    target_cols.append(col)

    for col in SINGLE_VALUE_SIGS:
        if col in g.columns and col in j.columns:
            target_cols.append(col)

    print(f"  Target columns: {len(target_cols)}")

    # Build columnar data
    data = {
        "g": list(common_gages),
        "lt": [round(v, 4) if pd.notna(v) else None
               for v in g["latitude"].astype(float).values],
        "ln": [round(v, 4) if pd.notna(v) else None
               for v in g["longitude"].astype(float).values],
        "r": {},
        "j": {},
    }

    for col in target_cols:
        gv = pd.to_numeric(g[col], errors="coerce").values
        jv = pd.to_numeric(j[col], errors="coerce").values
        data["r"][col] = [round(float(v), 4) if np.isfinite(v) else None for v in gv]
        data["j"][col] = [round(float(v), 4) if np.isfinite(v) else None for v in jv]

    return data, target_cols


def build_html(data, target_cols):
    data_json = json.dumps(data, separators=(",", ":"))
    print(f"  JSON payload size: {len(data_json) / 1024 / 1024:.1f} MB")

    # Build signature options HTML
    optgroups = []
    for group_name, sigs in SIGNATURE_GROUPS.items():
        opts = []
        for sig in sigs:
            # Check at least one stat column exists
            if any(sig + s in target_cols for s in STATS):
                opts.append(f'<option value="{sig}">{sig}</option>')
        if opts:
            optgroups.append(
                f'<optgroup label="{group_name}">{"".join(opts)}</optgroup>'
            )

    # Single-value group
    sv_opts = []
    for sig in SINGLE_VALUE_SIGS:
        if sig in target_cols:
            sv_opts.append(f'<option value="{sig}">{sig}</option>')
    if sv_opts:
        optgroups.append(
            f'<optgroup label="Single-Value (Recession/Elasticity)">{"".join(sv_opts)}</optgroup>'
        )

    sig_options_html = "\n".join(optgroups)

    html = HTML_TEMPLATE.replace("__DATA_JSON__", data_json)
    html = html.replace("__SIG_OPTIONS__", sig_options_html)

    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"  Written to: {OUTPUT_PATH}")
    print(f"  File size: {os.path.getsize(OUTPUT_PATH) / 1024 / 1024:.1f} MB")


HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Julia vs Golden R - Signature Comparison</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #1a1a2e; color: #e0e0e0; }
.toolbar { display: flex; align-items: center; gap: 16px; padding: 12px 20px; background: #16213e; border-bottom: 1px solid #0f3460; flex-wrap: wrap; }
.toolbar label { font-size: 13px; color: #a0a0c0; }
.toolbar select { padding: 6px 10px; border-radius: 4px; border: 1px solid #0f3460; background: #1a1a2e; color: #e0e0e0; font-size: 14px; cursor: pointer; }
.toolbar select:disabled { opacity: 0.4; cursor: default; }
#summary { font-size: 14px; padding: 8px 20px; background: #16213e; border-bottom: 1px solid #0f3460; }
#summary .r2-val { font-weight: bold; font-size: 16px; }
#summary .r2-perfect { color: #4ade80; }
#summary .r2-good { color: #facc15; }
#summary .r2-poor { color: #f87171; }
#selection-info { font-size: 13px; padding: 4px 20px; background: #0f3460; color: #facc15; display: none; }
.maps-row { display: flex; height: 45vh; min-height: 300px; }
.map-col { flex: 1; position: relative; }
.map-col .map-label { position: absolute; top: 10px; left: 50%; transform: translateX(-50%); z-index: 1000; background: rgba(22,33,62,0.9); color: #e0e0e0; padding: 4px 12px; border-radius: 4px; font-size: 13px; font-weight: 600; pointer-events: none; }
.map-container { width: 100%; height: 100%; }
.legend { position: absolute; bottom: 20px; right: 20px; z-index: 1000; background: rgba(22,33,62,0.92); padding: 8px 12px; border-radius: 6px; font-size: 11px; }
.legend .grad { width: 20px; height: 120px; margin: 4px 0; }
.legend .labels { display: flex; flex-direction: column; justify-content: space-between; height: 120px; margin-left: 4px; }
.legend-row { display: flex; align-items: stretch; }
#scatter-container { height: 45vh; min-height: 350px; padding: 0 10px 10px 10px; background: #1a1a2e; }
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
      <option value="_mk_pval">Mann-Kendall P-value</option>
    </select>
  </div>
</div>
<div id="summary"></div>
<div id="selection-info"></div>

<div class="maps-row">
  <div class="map-col">
    <div class="map-label">Golden R Reference (Feb 2026)</div>
    <div id="map-r" class="map-container"></div>
  </div>
  <div class="map-col">
    <div class="map-label">Julia Benchmark</div>
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
  "log_a_seasonality_minimum_all", "log_a_seasonality_minimum_first_half",
  "log_a_seasonality_minimum_last_half", "log_a_seasonality_amplitude_all",
  "log_a_seasonality_amplitude_first_half", "log_a_seasonality_amplitude_last_half"
]);

// Viridis-ish palette (8 stops)
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

function percentile(arr, p) {
  let s = arr.slice().sort((a, b) => a - b);
  let idx = (p / 100) * (s.length - 1);
  let lo = Math.floor(idx), hi = Math.ceil(idx);
  return s[lo] + (s[hi] - s[lo]) * (idx - lo);
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

// Sync maps manually
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

// Create circle markers
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

// Legend (on Julia map)
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

// Plotly scatter
const scatterDiv = document.getElementById("scatter-container");
const selInfo = document.getElementById("selection-info");

// Highlight ring markers for selected gage
let highlightR = null;
let highlightJ = null;

// Build idx->marker lookup for fast highlight
const markerByIdxR = {};
const markerByIdxJ = {};
for (let m of markersR) markerByIdxR[m._idx] = m;
for (let m of markersJ) markerByIdxJ[m._idx] = m;

// Scatter point index -> DATA index mapping (rebuilt each update)
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

  // Zoom both maps
  syncing = true;
  mapR.setView(ll, 8, { animate: true });
  mapJ.setView(ll, 8, { animate: true });
  setTimeout(() => { syncing = false; }, 300);

  // Open popups
  let mR = markerByIdxR[dataIdx], mJ = markerByIdxJ[dataIdx];
  if (mR) mR.openPopup();
  if (mJ) mJ.openPopup();

  // Show selection info bar
  let sig = document.getElementById("sig-select").value;
  let stat = document.getElementById("stat-select").value;
  let isSingle = SINGLE_VALUE.has(sig);
  let colKey = isSingle ? sig : sig + stat;
  let rV = DATA.r[colKey] ? DATA.r[colKey][dataIdx] : null;
  let jV = DATA.j[colKey] ? DATA.j[colKey][dataIdx] : null;
  let diff = (rV !== null && jV !== null) ? (jV - rV).toPrecision(4) : "N/A";
  selInfo.innerHTML =
    `Selected: <strong>${DATA.g[dataIdx]}</strong> &nbsp;|&nbsp; ` +
    `Original = ${rV !== null ? rV : "NA"} &nbsp;|&nbsp; ` +
    `New = ${jV !== null ? jV : "NA"} &nbsp;|&nbsp; ` +
    `Diff = ${diff} &nbsp; ` +
    `<span style="cursor:pointer;text-decoration:underline;" onclick="clearSelection()">[clear]</span>`;
  selInfo.style.display = "block";
}

function clearSelection() {
  clearHighlight();
  // Reset scatter highlight
  let colKey = getCurrentColKey();
  let rVals = DATA.r[colKey], jVals = DATA.j[colKey];
  if (!rVals) return;
  let colors = scatterToData.map(() => "rgba(99,179,237,0.5)");
  let sizes = scatterToData.map(() => 3);
  Plotly.restyle(scatterDiv, { "marker.color": [colors], "marker.size": [sizes] }, [0]);
}

function getCurrentColKey() {
  let sig = document.getElementById("sig-select").value;
  let stat = document.getElementById("stat-select").value;
  return SINGLE_VALUE.has(sig) ? sig : sig + stat;
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
    document.getElementById("summary").innerHTML = `<em>Column "${colKey}" not found in data.</em>`;
    return;
  }

  // Stats
  let nR_NA = rVals.filter(v => v === null).length;
  let nJ_NA = jVals.filter(v => v === null).length;
  let r2 = r2Identity(rVals, jVals);
  let nPaired = 0;
  for (let i = 0; i < N; i++) if (rVals[i] !== null && jVals[i] !== null) nPaired++;

  let r2Str = r2 !== null ? r2.toFixed(4) : "N/A";
  let r2Class = r2 === null ? "" : r2 >= 0.999 ? "r2-perfect" : r2 >= 0.99 ? "r2-good" : "r2-poor";
  document.getElementById("summary").innerHTML =
    `<strong>${colKey}</strong> &mdash; ` +
    `R&sup2; = <span class="r2-val ${r2Class}">${r2Str}</span> &nbsp;|&nbsp; ` +
    `N paired = ${nPaired} &nbsp;|&nbsp; ` +
    `Original NA = ${nR_NA}, New NA = ${nJ_NA}`;

  // Compute color scale
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

  // Update markers
  for (let m of markersR) {
    let v = rVals[m._idx];
    m.setStyle({ fillColor: colorFor(v) });
    m.unbindPopup();
    m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>Original: ${v !== null ? v : "NA"}`);
  }
  for (let m of markersJ) {
    let v = jVals[m._idx];
    m.setStyle({ fillColor: colorFor(v) });
    m.unbindPopup();
    m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>New: ${v !== null ? v : "NA"}`);
  }

  drawLegend(vMin, vMax, isPval);

  // Scatterplot - build index mapping
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
    hovertemplate: "<b>%{text}</b><br>Original: %{x:.4f}<br>New: %{y:.4f}<extra></extra>"
  }], {
    shapes: [{
      type: "line", x0: pLo, y0: pLo, x1: pHi, y1: pHi,
      line: { color: "rgba(250,200,50,0.5)", width: 1.5, dash: "dash" }
    }],
    xaxis: { title: "Original Filter (R)", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: "New Filter (Julia)", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `${colKey}  (R\u00b2 = ${r2Str},  N = ${nPaired})`, font: { color: "#e0e0e0", size: 14 } },
    plot_bgcolor: "#1a1a2e",
    paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 50, l: 60, r: 20 },
    font: { color: "#a0a0c0" }
  }, { responsive: true });
}

document.getElementById("sig-select").addEventListener("change", update);
document.getElementById("stat-select").addEventListener("change", update);

// Initial render (must happen before registering plotly_click)
update();

// Click on scatter point -> zoom maps to that gage
// Registered AFTER first Plotly.react() so the .on() method exists
scatterDiv.on("plotly_click", function(eventData) {
  if (!eventData || !eventData.points || eventData.points.length === 0) return;
  let ptIdx = eventData.points[0].pointIndex;
  let dataIdx = scatterToData[ptIdx];
  if (dataIdx === undefined) return;

  // Highlight the clicked point on the scatter
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

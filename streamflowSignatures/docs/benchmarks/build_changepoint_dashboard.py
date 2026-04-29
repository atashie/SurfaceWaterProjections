"""
Build an interactive HTML dashboard for Pettitt changepoint detection results.

Visualizes Pettitt test results across all signatures:
- Geographic maps: significance (p<0.05), changepoint year, pre/post means, etc.
- CP year histogram (all vs significant gages)
- Pre vs post mean scatterplot colored by significance

Usage:
    python build_changepoint_dashboard.py
"""

import pandas as pd
import numpy as np
import json
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

# ---------------------------------------------------------------------------
# Signature groups (same as other dashboards)
# ---------------------------------------------------------------------------
SIGNATURE_GROUPS = {
    "Flow Volumes": ["Qann", "Qwin", "Qspr", "Qsum", "Qfal"],
    "Flow Percentiles": ["Q1", "Q5", "Q10", "Q20", "Q25", "Q30", "Q40", "Q50",
                         "Q60", "Q70", "Q75", "Q80", "Q90", "Q95", "Q99", "Q95_Q10"],
    "Flow Timing": ["D1_day", "D5_day", "D10_day", "D20_day", "D30_day", "D40_day",
                    "D50_day", "D60_day", "D70_day", "D80_day", "D90_day", "D95_day",
                    "D99_day", "D25_to_D75", "Dmax"],
    "Baseflow": ["BFI_Eckhardt", "BFI_LyneHollick", "BFI_Eckhardt_param", "BFI_LyneHollick_param"],
    "Flashiness": ["flashinessRB"],
    "Pulses": ["TQmean", "n_low_pulses_all", "dur_low_pulses_all", "n_high_pulses_all",
               "dur_high_pulses_all", "n_low_pulses_year", "dur_low_pulses_year",
               "n_high_pulses_year", "dur_high_pulses_year",
               "Flow_Reversals_annual", "Flow_Reversals_winter", "Flow_Reversals_spring",
               "Flow_Reversals_summer", "Flow_Reversals_fall", "negative_ann"],
    "FDC": ["FDCall", "FDC90th", "FDCmid"],
    "Recession": ["log_a_pointcloud", "log_a_events", "b_pointcloud", "b_events",
                  "concavity", "n_recession_events", "alpha_linear"],
    "Runoff Ratios": ["annual_runoff_ratio", "winter_runoff_ratio", "spring_runoff_ratio",
                      "summer_runoff_ratio", "fall_runoff_ratio"],
    "Elasticity": ["elasticity_rolling", "elasticity_annual"],
    "Q-P Seasonality": ["qp_slope_sd", "qp_bimodality"],
    "Storage": ["avg_storage"],
}

CP_SUFFIXES = [
    # Pettitt test (8)
    "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean", "_pettitt_post_mean",
    "_pettitt_delta_mean", "_pettitt_pct_change", "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval",
]


def safe_val(v):
    """Convert value for JSON: NaN/Inf -> None, else round to 6 digits."""
    if v is None or (isinstance(v, float) and (np.isnan(v) or np.isinf(v))):
        return None
    return round(float(v), 6)


def main():
    # Load benchmark output
    csv_path = os.path.join(SCRIPT_DIR, "julia_signatures.csv")
    print(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path, low_memory=False)
    df["gage_id"] = df["gage_id"].astype(str).str.strip()
    print(f"  Loaded {len(df)} gages, {len(df.columns)} columns")

    # Extract lat/lon
    lat_col = "dec_lat_va" if "dec_lat_va" in df.columns else "latitude"
    lon_col = "dec_long_va" if "dec_long_va" in df.columns else "longitude"
    lats = [safe_val(v) for v in df[lat_col].values] if lat_col in df.columns else [None] * len(df)
    lons = [safe_val(v) for v in df[lon_col].values] if lon_col in df.columns else [None] * len(df)

    # Find all signatures that have Pettitt columns
    sigs_with_cp = []
    for group_name, sigs in SIGNATURE_GROUPS.items():
        for sig in sigs:
            pval_col = f"{sig}_pettitt_pval"
            if pval_col in df.columns:
                sigs_with_cp.append((group_name, sig))

    print(f"  Found {len(sigs_with_cp)} signatures with Pettitt columns")

    # Build columnar data for all changepoint columns
    cp_data = {}
    for _, sig in sigs_with_cp:
        for suf in CP_SUFFIXES:
            col = sig + suf
            if col in df.columns:
                cp_data[col] = [safe_val(v) for v in df[col].values]

    # Build signature options HTML
    sig_options = []
    current_group = None
    for group_name, sig in sigs_with_cp:
        if group_name != current_group:
            if current_group is not None:
                sig_options.append("</optgroup>")
            sig_options.append(f'<optgroup label="{group_name}">')
            current_group = group_name
        sig_options.append(f'<option value="{sig}">{sig}</option>')
    if current_group is not None:
        sig_options.append("</optgroup>")
    sig_options_html = "\n".join(sig_options)

    # Compute global Pettitt summary stats
    pettitt_sig_count = 0
    pettitt_sig_total = 0
    for _, sig in sigs_with_cp:
        pettitt_pval_col = f"{sig}_pettitt_pval"
        if pettitt_pval_col in df.columns:
            pvals = df[pettitt_pval_col].dropna()
            pettitt_sig_total += len(pvals)
            pettitt_sig_count += (pvals < 0.05).sum()

    summary_json = json.dumps({
        "total_sigs": len(sigs_with_cp),
        "total_gages": len(df),
        "pettitt_sig_pct": round(pettitt_sig_count / pettitt_sig_total * 100, 1) if pettitt_sig_total > 0 else 0,
        "pettitt_sig_count": int(pettitt_sig_count),
        "pettitt_sig_total": int(pettitt_sig_total),
    })

    # Build data payload
    data = {
        "g": df["gage_id"].tolist(),
        "lt": lats,
        "ln": lons,
        "cp": cp_data,
    }
    data_json = json.dumps(data, separators=(",", ":"))
    print(f"  JSON payload size: {len(data_json) / 1024 / 1024:.1f} MB")

    # Build HTML
    html = HTML_TEMPLATE
    html = html.replace("__DATA_JSON__", data_json)
    html = html.replace("__SIG_OPTIONS__", sig_options_html)
    html = html.replace("__SUMMARY_JSON__", summary_json)

    output_path = os.path.join(SCRIPT_DIR, "changepoint_dashboard.html")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"  Dashboard saved to {output_path}")
    print(f"  File size: {os.path.getsize(output_path) / 1024 / 1024:.1f} MB")


# ---------------------------------------------------------------------------
# HTML Template
# ---------------------------------------------------------------------------
HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Pettitt Changepoint Dashboard</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
* { box-sizing: border-box; margin: 0; padding: 0; }
body { background: #1a1a2e; color: #e0e0e0; font-family: 'Segoe UI', system-ui, sans-serif; }
.toolbar {
  background: #16213e; padding: 8px 16px; display: flex; align-items: center; gap: 16px;
  border-bottom: 1px solid #0f3460; flex-wrap: wrap;
}
.toolbar label { font-size: 13px; color: #a0a0c0; }
.toolbar select {
  background: #0f3460; color: #e0e0e0; border: 1px solid #1a3a6e; border-radius: 4px;
  padding: 4px 8px; font-size: 13px; min-width: 160px;
}
.toolbar .title { font-size: 16px; font-weight: 700; color: #63b3ed; margin-right: 20px; }
#summary {
  background: #16213e; padding: 6px 16px; display: flex; align-items: center; gap: 16px;
  border-bottom: 1px solid #0f3460; font-size: 12px; flex-wrap: wrap;
}
#summary .stat-label { color: #a0a0c0; }
#summary .stat-value { color: #63b3ed; font-weight: 600; }
.pettitt-badge { background: #0e4429; color: #56d364; padding: 2px 8px; border-radius: 3px; font-weight: 600; font-size: 11px; }

#global-summary {
  background: #0f1a2e; padding: 8px 16px; border-bottom: 1px solid #0f3460;
  display: flex; gap: 20px; align-items: center; font-size: 12px; flex-wrap: wrap;
}
#global-summary .gs-title { color: #a0a0c0; font-weight: 600; }
#global-summary .gs-bar { display: flex; height: 16px; border-radius: 3px; overflow: hidden; width: 260px; }
#global-summary .gs-seg { display: flex; align-items: center; justify-content: center; font-size: 10px; font-weight: 600; }

.maps-row { display: flex; height: 38vh; min-height: 260px; }
.map-panel { flex: 1; position: relative; }
.map-panel .map-label {
  position: absolute; top: 6px; left: 50%; transform: translateX(-50%);
  z-index: 999; background: rgba(22,33,62,0.85); padding: 3px 10px; border-radius: 4px;
  font-size: 12px; font-weight: 600; color: #63b3ed; pointer-events: none;
  white-space: nowrap;
}
.map-div { width: 100%; height: 100%; }
.legend {
  background: rgba(22,33,62,0.9); padding: 6px 8px; border-radius: 4px;
  font-size: 11px; color: #e0e0e0;
}
.legend .legend-row { display: flex; gap: 4px; align-items: stretch; }
.legend .labels { display: flex; flex-direction: column; justify-content: space-between; font-size: 10px; }

.model-legend {
  background: rgba(22,33,62,0.9); padding: 8px; border-radius: 4px;
  font-size: 11px; color: #e0e0e0;
}
.model-legend .ml-title { font-weight: 600; margin-bottom: 4px; }
.model-legend .ml-item { display: flex; align-items: center; gap: 6px; margin: 2px 0; }
.model-legend .ml-dot { width: 12px; height: 12px; border-radius: 50%; }

.charts-row { display: flex; height: 47vh; min-height: 360px; }
.chart-panel { flex: 1; min-width: 0; }

#selection-info {
  position: fixed; bottom: 0; left: 0; right: 0; z-index: 9999;
  background: rgba(22,33,62,0.95); padding: 6px 16px; color: #facc15;
  font-size: 12px; display: none; border-top: 1px solid #0f3460;
  line-height: 1.6;
}
</style>
</head>
<body>

<div class="toolbar">
  <span class="title">Pettitt Changepoint Dashboard</span>
  <label>Signature:
    <select id="sig-select">__SIG_OPTIONS__</select>
  </label>
  <label>Map right:
    <select id="map-metric">
      <option value="pettitt_cp_year">CP Year</option>
      <option value="pettitt_pval">p-value</option>
      <option value="pettitt_pre_mean">Pre-CP Mean</option>
      <option value="pettitt_post_mean">Post-CP Mean</option>
      <option value="pettitt_delta_mean">Delta Mean</option>
      <option value="pettitt_pct_change">% Change</option>
      <option value="pettitt_pre_mk_pval">Pre MK p-val</option>
      <option value="pettitt_post_mk_pval">Post MK p-val</option>
    </select>
  </label>
</div>

<div id="summary"></div>
<div id="global-summary"></div>

<div class="maps-row">
  <div class="map-panel">
    <div class="map-label" id="map-label-left">Pettitt Significance</div>
    <div class="map-div" id="map-model"></div>
  </div>
  <div class="map-panel">
    <div class="map-label" id="map-label-right">Changepoint Year</div>
    <div class="map-div" id="map-metric-div"></div>
  </div>
</div>

<div class="charts-row">
  <div class="chart-panel" id="chart-left"></div>
  <div class="chart-panel" id="chart-right"></div>
</div>

<div id="selection-info"></div>

<script>
const DATA = __DATA_JSON__;
const GLOBAL = __SUMMARY_JSON__;
const N = DATA.g.length;

const SIG_COLOR = "#56d364";  // Significant (p<0.05)
const NONSIG_COLOR = "#555";  // Non-significant
const NA_COLOR = "#333";

const VIRIDIS = [
  [68,1,84],[72,35,116],[64,67,135],[52,95,141],
  [33,144,140],[53,183,121],[143,215,68],[253,231,37]
];
const DIVERGING = [
  [33,102,172],[67,147,195],[146,197,222],[209,229,240],
  [253,219,199],[244,165,130],[214,96,77],[178,24,43]
];
const YEAR_PALETTE = [
  [253,231,37],[143,215,68],[53,183,121],[33,144,140],
  [52,95,141],[64,67,135],[72,35,116],[68,1,84]
];
const PVAL_PALETTE = [
  [178,24,43],[214,96,77],[244,165,130],[253,219,199],
  [209,229,240],[146,197,222],[67,147,195],[33,102,172]
];

function interpColor(pal, t) {
  t = Math.max(0, Math.min(1, t));
  let idx = t * (pal.length - 1);
  let lo = Math.floor(idx), hi = Math.min(lo + 1, pal.length - 1);
  let f = idx - lo;
  let r = Math.round(pal[lo][0] + (pal[hi][0] - pal[lo][0]) * f);
  let g = Math.round(pal[lo][1] + (pal[hi][1] - pal[lo][1]) * f);
  let b = Math.round(pal[lo][2] + (pal[hi][2] - pal[lo][2]) * f);
  return `rgb(${r},${g},${b})`;
}

function percentile(arr, p) {
  let s = arr.slice().sort((a, b) => a - b);
  let idx = (p / 100) * (s.length - 1);
  let lo = Math.floor(idx), hi = Math.ceil(idx);
  return lo === hi ? s[lo] : s[lo] + (s[hi] - s[lo]) * (idx - lo);
}

function fmt(v, digits) {
  if (v === null || v === undefined) return "N/A";
  if (Math.abs(v) >= 1000) return v.toFixed(1);
  if (Math.abs(v) >= 1) return v.toFixed(digits || 3);
  return v.toFixed(digits || 4);
}

// ---- Maps ----
const mapModel = L.map("map-model", { preferCanvas: true }).setView([41.3, -95.6], 4);
const mapMetric = L.map("map-metric-div", { preferCanvas: true }).setView([41.3, -95.6], 4);
L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapModel);
L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapMetric);

let syncing = false;
function syncMap(src, tgt) {
  src.on("moveend zoomend", function () {
    if (syncing) return;
    syncing = true;
    tgt.setView(src.getCenter(), src.getZoom(), { animate: false });
    syncing = false;
  });
}
syncMap(mapModel, mapMetric);
syncMap(mapMetric, mapModel);

let markersL = [], markersR = [];
for (let i = 0; i < N; i++) {
  if (DATA.lt[i] === null || DATA.ln[i] === null) continue;
  let ll = [DATA.lt[i], DATA.ln[i]];
  let mL = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapModel);
  let mR = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapMetric);
  mL._idx = i; mR._idx = i;
  mL.on("click", () => highlightGage(i));
  mR.on("click", () => highlightGage(i));
  markersL.push(mL);
  markersR.push(mR);
}

const mlDiv = document.createElement("div");
mlDiv.className = "model-legend";
mlDiv.innerHTML = `
  <div class="ml-title">Pettitt Test</div>
  <div class="ml-item"><div class="ml-dot" style="background:${SIG_COLOR}"></div> p &lt; 0.05 (significant)</div>
  <div class="ml-item"><div class="ml-dot" style="background:${NONSIG_COLOR}"></div> p &ge; 0.05</div>
  <div class="ml-item"><div class="ml-dot" style="background:${NA_COLOR}"></div> N/A</div>
`;
const MLControl = L.Control.extend({ onAdd: function() { return mlDiv; } });
new MLControl({ position: "bottomright" }).addTo(mapModel);

const legDiv = document.createElement("div");
legDiv.className = "legend";
legDiv.innerHTML = `
  <div style="font-weight:600;margin-bottom:4px;" id="leg-title">Value</div>
  <div class="legend-row">
    <canvas id="legend-canvas" width="20" height="120"></canvas>
    <div class="labels">
      <span id="leg-max"></span>
      <span id="leg-mid"></span>
      <span id="leg-min"></span>
    </div>
  </div>`;
const LegControl = L.Control.extend({ onAdd: function() { return legDiv; } });
new LegControl({ position: "bottomright" }).addTo(mapMetric);

function drawLegend(vMin, vMax, pal, title) {
  document.getElementById("leg-title").textContent = title || "Value";
  let canvas = document.getElementById("legend-canvas");
  let ctx = canvas.getContext("2d");
  for (let y = 0; y < 120; y++) {
    let t = 1 - y / 119;
    ctx.fillStyle = interpColor(pal, t);
    ctx.fillRect(0, y, 20, 1);
  }
  document.getElementById("leg-max").textContent = fmt(vMax, 1);
  document.getElementById("leg-mid").textContent = fmt((vMin + vMax) / 2, 1);
  document.getElementById("leg-min").textContent = fmt(vMin, 1);
}

// ---- Highlight ----
let highlightRings = [];
function highlightGage(idx) {
  highlightRings.forEach(r => { r.remove(); });
  highlightRings = [];
  if (DATA.lt[idx] === null || DATA.ln[idx] === null) return;
  let ll = [DATA.lt[idx], DATA.ln[idx]];
  let r1 = L.circleMarker(ll, { radius: 10, weight: 3, color: "#facc15", fill: false }).addTo(mapModel);
  let r2 = L.circleMarker(ll, { radius: 10, weight: 3, color: "#facc15", fill: false }).addTo(mapMetric);
  highlightRings.push(r1, r2);

  let sig = document.getElementById("sig-select").value;
  let pY = (DATA.cp[sig+"_pettitt_cp_year"]||[])[idx];
  let pp = (DATA.cp[sig+"_pettitt_pval"]||[])[idx];
  let preM = (DATA.cp[sig+"_pettitt_pre_mean"]||[])[idx];
  let postM = (DATA.cp[sig+"_pettitt_post_mean"]||[])[idx];
  let pdm = (DATA.cp[sig+"_pettitt_delta_mean"]||[])[idx];
  let ppct = (DATA.cp[sig+"_pettitt_pct_change"]||[])[idx];
  let info = `<b>${DATA.g[idx]}</b> &mdash; CP year: ${fmt(pY,0)}, p=${fmt(pp)}, pre=${fmt(preM)}, post=${fmt(postM)}, &Delta;mean=${fmt(pdm)}, %chg=${fmt(ppct,1)}`;

  let sel = document.getElementById("selection-info");
  sel.innerHTML = info;
  sel.style.display = "block";
}

// ---- Global summary bar ----
function renderGlobalSummary() {
  let el = document.getElementById("global-summary");
  let sigPct = GLOBAL.pettitt_sig_pct;
  let nonPct = (100 - sigPct).toFixed(1);
  el.innerHTML = `
    <span class="gs-title">${GLOBAL.total_sigs} sigs &times; ${GLOBAL.total_gages} gages (${GLOBAL.pettitt_sig_total.toLocaleString()} evaluations)</span>
    <div class="gs-bar">
      <div class="gs-seg" style="width:${sigPct}%;background:#0e4429;color:#56d364">p&lt;0.05 ${sigPct}%</div>
      <div class="gs-seg" style="width:${nonPct}%;background:#333;color:#888">n.s. ${nonPct}%</div>
    </div>
    <span style="color:#56d364;font-weight:600">Significant: ${GLOBAL.pettitt_sig_count.toLocaleString()} / ${GLOBAL.pettitt_sig_total.toLocaleString()} (${sigPct}%)</span>
  `;
}
renderGlobalSummary();

// ---- Main update ----
function update() {
  let sig = document.getElementById("sig-select").value;
  let metricKey = document.getElementById("map-metric").value;
  let pvalCol = sig + "_pettitt_pval";
  let metricCol = sig + "_" + metricKey;

  document.getElementById("map-label-left").textContent = sig + " — Pettitt Significance";
  let metricLabel = document.getElementById("map-metric").selectedOptions[0].text;
  document.getElementById("map-label-right").textContent = sig + " — " + metricLabel;

  let pvalVals = DATA.cp[pvalCol] || [];
  let metVals = DATA.cp[metricCol] || [];

  // Left map: green = significant (p<0.05), grey = not significant
  for (let m of markersL) {
    let v = pvalVals[m._idx];
    m.setStyle({ fillColor: v === null ? NA_COLOR : (v < 0.05 ? SIG_COLOR : NONSIG_COLOR) });
  }

  let validMetVals = metVals.filter(v => v !== null && v !== undefined);
  let vMin, vMax, pal;
  let isDivergent = metricKey.includes("delta") || metricKey.includes("slope") || metricKey.includes("pct_change");
  let isYear = metricKey.includes("cp_year");
  let isPval = metricKey.includes("pval");

  if (isYear) {
    pal = YEAR_PALETTE;
    vMin = validMetVals.length > 0 ? Math.min(...validMetVals) : 1990;
    vMax = validMetVals.length > 0 ? Math.max(...validMetVals) : 2014;
  } else if (isPval) {
    pal = PVAL_PALETTE;
    vMin = 0; vMax = 1;
  } else if (isDivergent && validMetVals.length > 10) {
    pal = DIVERGING;
    let absMax = percentile(validMetVals.map(Math.abs), 98);
    vMin = -absMax; vMax = absMax;
  } else if (validMetVals.length > 10) {
    pal = VIRIDIS;
    vMin = percentile(validMetVals, 2);
    vMax = percentile(validMetVals, 98);
  } else {
    pal = VIRIDIS;
    vMin = validMetVals.length > 0 ? Math.min(...validMetVals) : 0;
    vMax = validMetVals.length > 0 ? Math.max(...validMetVals) : 1;
  }
  if (vMax === vMin) vMax = vMin + 1;

  for (let m of markersR) {
    let v = metVals[m._idx];
    m.setStyle({ fillColor: v === null ? NA_COLOR : interpColor(pal, (v - vMin) / (vMax - vMin)) });
  }
  drawLegend(vMin, vMax, pal, metricLabel);

  // -- Summary bar --
  let pettittPvals = DATA.cp[sig + "_pettitt_pval"] || [];
  let pSigCount = 0, pTotal = 0;
  for (let v of pettittPvals) { if (v !== null) { pTotal++; if (v < 0.05) pSigCount++; } }

  let cpYears = DATA.cp[sig + "_pettitt_cp_year"] || [];
  let sigCpYears = [];
  for (let i = 0; i < N; i++) {
    if (cpYears[i] !== null && pettittPvals[i] !== null && pettittPvals[i] < 0.05) {
      sigCpYears.push(cpYears[i]);
    }
  }
  let medianCp = sigCpYears.length > 0 ? percentile(sigCpYears, 50) : null;

  document.getElementById("summary").innerHTML = `
    <span><span class="stat-label">Sig: </span><span class="stat-value">${sig}</span></span>
    <span class="pettitt-badge">p&lt;0.05: ${pSigCount}/${pTotal} (${pTotal>0?(pSigCount/pTotal*100).toFixed(1):"0"}%)</span>
    <span style="color:#a0a0c0;font-size:11px">Median CP year (sig): ${medianCp !== null ? medianCp.toFixed(0) : "N/A"}</span>
  `;

  // -- Left chart: CP year histogram (significant vs all) --
  let allCpYears = cpYears.filter(v => v !== null);

  let histAll = {
    x: allCpYears, type: "histogram", opacity: 0.35,
    marker: { color: "rgba(160,160,192,0.5)" },
    nbinsx: 25, name: "All gages"
  };
  let histSig = {
    x: sigCpYears, type: "histogram", opacity: 0.7,
    marker: { color: "rgba(86,211,100,0.7)" },
    nbinsx: 25, name: "p < 0.05"
  };

  Plotly.react("chart-left", [histAll, histSig], {
    barmode: "overlay",
    xaxis: { title: "Changepoint Year", gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: "Count", gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `${sig}: CP Year Distribution (N=${allCpYears.length})`, font: { color: "#e0e0e0", size: 13 } },
    plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 50, l: 50, r: 20 },
    font: { color: "#a0a0c0" },
    legend: { x: 0.65, y: 0.95, bgcolor: "rgba(22,33,62,0.8)", font: { size: 10, color: "#e0e0e0" } },
    showlegend: true
  }, { responsive: true });

  // -- Right chart: Pre vs Post Mean scatter (colored by significance) --
  let sx = [], sy = [], sColors = [], sText = [];
  let preVals = DATA.cp[sig + "_pettitt_pre_mean"] || [];
  let postVals = DATA.cp[sig + "_pettitt_post_mean"] || [];
  for (let i = 0; i < N; i++) {
    if (preVals[i] !== null && postVals[i] !== null) {
      sx.push(preVals[i]); sy.push(postVals[i]);
      let pp = pettittPvals[i];
      sColors.push(pp !== null && pp < 0.05 ? "#56d364" : "#555");
      sText.push(`${DATA.g[i]}<br>Pre: ${fmt(preVals[i])}<br>Post: ${fmt(postVals[i])}<br>p=${fmt(pp)}`);
    }
  }

  let allXY = sx.concat(sy).filter(v => isFinite(v));
  let pLo, pHi;
  if (allXY.length > 10) {
    pLo = percentile(allXY, 1); pHi = percentile(allXY, 99);
    let pad = (pHi - pLo) * 0.05;
    pLo -= pad; pHi += pad;
  } else {
    pLo = allXY.length > 0 ? Math.min(...allXY) - 1 : -1;
    pHi = allXY.length > 0 ? Math.max(...allXY) + 1 : 1;
  }

  Plotly.react("chart-right", [{
    x: sx, y: sy, text: sText,
    mode: "markers", type: "scattergl",
    marker: { size: 4, color: sColors, opacity: 0.7 },
    hovertemplate: "%{text}<extra></extra>"
  }], {
    shapes: [{
      type: "line", x0: pLo, y0: pLo, x1: pHi, y1: pHi,
      line: { color: "rgba(250,200,50,0.4)", width: 1.5, dash: "dash" }
    }],
    xaxis: { title: "Pre-CP Mean", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: "Post-CP Mean", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `${sig}: Pre vs Post Mean (N=${sx.length})`, font: { color: "#e0e0e0", size: 13 } },
    plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 50, l: 60, r: 20 },
    font: { color: "#a0a0c0" },
    showlegend: false
  }, { responsive: true });

  // Click scatter to highlight gage on maps
  document.getElementById("chart-right").on("plotly_click", function(ed) {
    if (!ed || !ed.points || !ed.points.length) return;
    let px = ed.points[0].x, py = ed.points[0].y;
    for (let j = 0; j < N; j++) {
      if (preVals[j] !== null && postVals[j] !== null &&
          Math.abs(preVals[j] - px) < 1e-4 && Math.abs(postVals[j] - py) < 1e-4) {
        highlightGage(j);
        if (DATA.lt[j] !== null) mapModel.setView([DATA.lt[j], DATA.ln[j]], 8, { animate: true });
        break;
      }
    }
  });
}

document.getElementById("sig-select").addEventListener("change", update);
document.getElementById("map-metric").addEventListener("change", update);

update();
</script>
</body>
</html>
'''

if __name__ == "__main__":
    main()

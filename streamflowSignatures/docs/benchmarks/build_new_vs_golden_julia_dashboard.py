"""
Build an interactive HTML dashboard comparing a new Julia benchmark against
the Golden Julia (canonical) reference.

Shows:
  - R² scatterplot for all common columns (regression check)
  - Histogram + geographic map for new-only columns (distribution check)
  - Global tier summary across all common columns

Usage:
    python docs/benchmarks/build_new_vs_golden_julia_dashboard.py

Output:
    docs/benchmarks/new_vs_golden_julia_dashboard.html
"""
import pandas as pd
import numpy as np
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

GOLDEN_PATH = os.path.join(PROJECT_ROOT, "golden-outputs",
                           "streamflow_signatures_julia_apr2026.csv")
NEW_PATH = os.path.join(SCRIPT_DIR, "julia_signatures.csv")
OUTPUT_PATH = os.path.join(SCRIPT_DIR, "new_vs_golden_julia_dashboard.html")

# Metadata / non-signature columns to exclude from comparison
META_COLS = {
    "gage_id", "latitude", "longitude", "drainage_area_km2",
    "area_normalized", "country", "state_province", "gage_name",
    "min_year", "max_year", "n_years", "n_qualifying_years",
    "human_interference_class", "CLASS", "NDAMS_2009", "MAJ_DDENS_2009",
    "STOR_NID_2009", "IMPNLCD06", "DEVNLCD06", "FRESHW_WITHDRAWAL",
    "HYDRO_DISTURB_INDX", "RHBN", "REGULATED",
    "ice_affected_days_total",
    # QA flags
    "flagged_for_qann_range", "flagged_for_bfi_eckhardt_range",
    "flagged_for_bfi_lynehollick_range", "flagged_for_flashiness_range",
    "flagged_for_tqmean_range", "flagged_for_d50_range",
    "flagged_for_elasticity_range", "flagged_for_runoff_ratio_range",
    "flagged_for_seasonal_sum", "flagged_for_percentile_order",
    "flagged_for_timing_order", "flagged_for_high_na",
}

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
    "Baseflow (Fixed)": ["BFI_Eckhardt", "BFI_LyneHollick"],
    "Baseflow (Recession-Parameterized)": ["BFI_Eckhardt_param", "BFI_LyneHollick_param"],
    "Flashiness": ["flashinessRB"],
    "Pulses": [
        "TQmean", "n_low_pulses_all", "dur_low_pulses_all",
        "n_high_pulses_all", "dur_high_pulses_all",
        "n_low_pulses_year", "dur_low_pulses_year",
        "n_high_pulses_year", "dur_high_pulses_year",
        "Flow_Reversals_annual", "Flow_Reversals_winter",
        "Flow_Reversals_spring", "Flow_Reversals_summer",
        "Flow_Reversals_fall", "negative_ann",
    ],
    "FDC": ["FDCall", "FDC90th", "FDCmid"],
    "Recession": [
        "log_a_pointcloud", "log_a_events", "b_pointcloud",
        "b_events", "concavity", "n_recession_events", "alpha_linear",
    ],
    "Runoff Ratios": [
        "annual_runoff_ratio", "winter_runoff_ratio",
        "spring_runoff_ratio", "summer_runoff_ratio",
        "fall_runoff_ratio",
    ],
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


def r2_identity(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10:
        return np.nan
    xm, ym = x[mask], y[mask]
    ss_res = np.sum((ym - xm) ** 2)
    ss_tot = np.sum((ym - np.mean(ym)) ** 2)
    if ss_tot == 0:
        return 1.0 if ss_res == 0 else np.nan
    return 1.0 - ss_res / ss_tot


def build_data():
    print("Loading Golden Julia output...")
    golden_df = pd.read_csv(GOLDEN_PATH, low_memory=False)
    golden_df["gage_id"] = golden_df["gage_id"].astype(str).str.strip()
    print(f"  {golden_df.shape[0]} gages x {golden_df.shape[1]} cols")

    print("Loading New Julia benchmark output...")
    new_df = pd.read_csv(NEW_PATH, low_memory=False)
    new_df["gage_id"] = new_df["gage_id"].astype(str).str.strip()
    print(f"  {new_df.shape[0]} gages x {new_df.shape[1]} cols")

    # Common gages
    common_gages = sorted(set(golden_df["gage_id"]) & set(new_df["gage_id"]))
    print(f"  Common gages: {len(common_gages)}")

    golden = golden_df[golden_df["gage_id"].isin(common_gages)].set_index("gage_id").loc[common_gages]
    new = new_df[new_df["gage_id"].isin(common_gages)].set_index("gage_id").loc[common_gages]

    # Identify signature columns
    golden_sig = set(golden.columns) - META_COLS
    new_sig = set(new.columns) - META_COLS
    common_cols = sorted(golden_sig & new_sig)
    new_only_cols = sorted(new_sig - golden_sig)

    print(f"  Common signature columns: {len(common_cols)}")
    print(f"  New-only signature columns: {len(new_only_cols)}")

    # Build column lists organized by group
    target_common = []
    target_new_only = []

    for sigs in SIGNATURE_GROUPS.values():
        for sig in sigs:
            for stat in STATS:
                col = sig + stat
                if col in common_cols:
                    target_common.append(col)
                elif col in new_only_cols:
                    target_new_only.append(col)

    for col in SINGLE_VALUE_SIGS:
        if col in common_cols:
            target_common.append(col)
        elif col in new_only_cols:
            target_new_only.append(col)

    # Deduplicate preserving order
    target_common = list(dict.fromkeys(target_common))
    target_new_only = list(dict.fromkeys(target_new_only))

    # Pick up any remaining columns not in groups
    for col in common_cols:
        if col not in target_common:
            target_common.append(col)
    for col in new_only_cols:
        if col not in target_new_only:
            target_new_only.append(col)

    print(f"  Target common: {len(target_common)}, Target new-only: {len(target_new_only)}")

    # Coordinates
    lat_col = golden["latitude"] if "latitude" in golden.columns else new.get("latitude")
    lon_col = golden["longitude"] if "longitude" in golden.columns else new.get("longitude")

    data = {
        "g": list(common_gages),
        "lt": [round(float(v), 4) if pd.notna(v) else None for v in lat_col.values],
        "ln": [round(float(v), 4) if pd.notna(v) else None for v in lon_col.values],
        "golden": {},
        "new": {},
        "common_cols": target_common,
        "new_only_cols": target_new_only,
    }

    all_cols = target_common + target_new_only
    for col in all_cols:
        if col in golden.columns:
            gv = pd.to_numeric(golden[col], errors="coerce").values
            data["golden"][col] = [round(float(v), 6) if np.isfinite(v) else None for v in gv]

        if col in new.columns:
            nv = pd.to_numeric(new[col], errors="coerce").values
            data["new"][col] = [round(float(v), 6) if np.isfinite(v) else None for v in nv]

    # Compute tier summary for common cols
    tier_counts = {"Perfect": 0, "Good": 0, "Poor": 0, "Low": 0,
                   "Very Low": 0, "Extremely Low": 0, "N/A": 0}
    for col in target_common:
        if col not in data["golden"] or col not in data["new"]:
            tier_counts["N/A"] += 1
            continue
        gv = np.array([v if v is not None else np.nan for v in data["golden"][col]])
        nv = np.array([v if v is not None else np.nan for v in data["new"][col]])
        r2 = r2_identity(gv, nv)
        if np.isnan(r2):
            tier_counts["N/A"] += 1
        elif r2 >= 0.999:
            tier_counts["Perfect"] += 1
        elif r2 >= 0.99:
            tier_counts["Good"] += 1
        elif r2 >= 0.95:
            tier_counts["Poor"] += 1
        elif r2 >= 0.9:
            tier_counts["Low"] += 1
        elif r2 >= 0.5:
            tier_counts["Very Low"] += 1
        else:
            tier_counts["Extremely Low"] += 1

    data["tier_counts"] = tier_counts

    print(f"\n  Tier summary (common columns):")
    for tier, count in tier_counts.items():
        if count > 0:
            print(f"    {tier}: {count}")

    return data


def build_html(data):
    data_json = json.dumps(data, separators=(",", ":"))
    print(f"  JSON payload size: {len(data_json) / 1024 / 1024:.1f} MB")

    # Build optgroups for common columns
    common_set = set(data["common_cols"])
    new_only_set = set(data["new_only_cols"])

    optgroups_common = []
    optgroups_new = []

    for group_name, sigs in SIGNATURE_GROUPS.items():
        opts_common = []
        opts_new = []
        for sig in sigs:
            has_common = any((sig + s) in common_set for s in STATS)
            has_new = any((sig + s) in new_only_set for s in STATS)
            if has_common:
                opts_common.append(f'<option value="{sig}">{sig}</option>')
            if has_new:
                opts_new.append(f'<option value="{sig}">{sig}</option>')
        if opts_common:
            optgroups_common.append(
                f'<optgroup label="{group_name}">{"".join(opts_common)}</optgroup>')
        if opts_new:
            optgroups_new.append(
                f'<optgroup label="{group_name} (NEW)">{"".join(opts_new)}</optgroup>')

    # Single-value sigs
    sv_common = []
    sv_new = []
    for sig in SINGLE_VALUE_SIGS:
        if sig in common_set:
            sv_common.append(f'<option value="{sig}">{sig}</option>')
        elif sig in new_only_set:
            sv_new.append(f'<option value="{sig}">{sig}</option>')
    if sv_common:
        optgroups_common.append(
            f'<optgroup label="Single-Value (Scalars)">{"".join(sv_common)}</optgroup>')
    if sv_new:
        optgroups_new.append(
            f'<optgroup label="Single-Value (Scalars) (NEW)">{"".join(sv_new)}</optgroup>')

    all_optgroups = optgroups_common + optgroups_new
    sig_options_html = "\n".join(all_optgroups)

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
<title>New Julia Benchmark vs Golden Julia (Validation)</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #1a1a2e; color: #e0e0e0; }
.toolbar { display: flex; align-items: center; gap: 16px; padding: 12px 20px; background: #16213e; border-bottom: 1px solid #0f3460; flex-wrap: wrap; }
.toolbar label { font-size: 13px; color: #a0a0c0; }
.toolbar select { padding: 6px 10px; border-radius: 4px; border: 1px solid #0f3460; background: #1a1a2e; color: #e0e0e0; font-size: 14px; cursor: pointer; }
.toolbar select:disabled { opacity: 0.4; cursor: default; }
#mode-badge { padding: 4px 10px; border-radius: 12px; font-size: 12px; font-weight: 700; }
.mode-comparison { background: #0f3460; color: #63b3ed; }
.mode-distribution { background: #1e3a29; color: #4ade80; }
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
.maps-row { display: flex; height: 35vh; min-height: 260px; }
.map-col { flex: 1; position: relative; }
.map-col .map-label { position: absolute; top: 10px; left: 50%; transform: translateX(-50%); z-index: 1000; background: rgba(22,33,62,0.9); color: #e0e0e0; padding: 4px 12px; border-radius: 4px; font-size: 13px; font-weight: 600; pointer-events: none; }
.map-container { width: 100%; height: 100%; }
.legend { position: absolute; bottom: 20px; right: 20px; z-index: 1000; background: rgba(22,33,62,0.92); padding: 8px 12px; border-radius: 6px; font-size: 11px; }
.legend .grad { width: 20px; height: 120px; margin: 4px 0; }
.legend .labels { display: flex; flex-direction: column; justify-content: space-between; height: 120px; margin-left: 4px; }
.legend-row { display: flex; align-items: stretch; }
#charts-row { display: flex; height: 45vh; min-height: 340px; }
#scatter-container { flex: 1; padding: 0 5px 10px 10px; background: #1a1a2e; }
#hist-container { flex: 1; padding: 0 10px 10px 5px; background: #1a1a2e; }
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
  <span id="mode-badge"></span>
</div>
<div id="summary"></div>
<div id="tier-bar"></div>
<div id="selection-info"></div>

<div class="maps-row">
  <div class="map-col">
    <div class="map-label" id="map-left-label">Golden Julia (Canonical)</div>
    <div id="map-left" class="map-container"></div>
  </div>
  <div class="map-col">
    <div class="map-label" id="map-right-label">New Benchmark</div>
    <div id="map-right" class="map-container"></div>
  </div>
</div>

<div id="charts-row">
  <div id="scatter-container"></div>
  <div id="hist-container"></div>
</div>

<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<script>
const DATA = __DATA_JSON__;

const COMMON_COLS = new Set(DATA.common_cols);
const NEW_ONLY_COLS = new Set(DATA.new_only_cols);

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
const GREEN_PALETTE = [
  [247,252,245],[229,245,224],[199,233,192],[161,217,155],
  [116,196,118],[65,171,93],[35,139,69],[0,109,44]
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

// Maps
const mapLeft = L.map("map-left", { preferCanvas: true }).setView([41.3, -95.6], 4);
const mapRight = L.map("map-right", { preferCanvas: true }).setView([41.3, -95.6], 4);

L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapLeft);
L.tileLayer("https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png", {
  attribution: '&copy; OSM &copy; CARTO', maxZoom: 18
}).addTo(mapRight);

let syncing = false;
function syncMap(source, target) {
  source.on("moveend zoomend", function () {
    if (syncing) return;
    syncing = true;
    target.setView(source.getCenter(), source.getZoom(), { animate: false });
    syncing = false;
  });
}
syncMap(mapLeft, mapRight);
syncMap(mapRight, mapLeft);

const markersLeft = [];
const markersRight = [];
const N = DATA.g.length;

for (let i = 0; i < N; i++) {
  if (DATA.lt[i] === null || DATA.ln[i] === null) continue;
  let ll = [DATA.lt[i], DATA.ln[i]];
  let mL = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapLeft);
  let mR = L.circleMarker(ll, { radius: 3, weight: 0.5, color: "#555", fillOpacity: 0.8 }).addTo(mapRight);
  mL._idx = i;
  mR._idx = i;
  markersLeft.push(mL);
  markersRight.push(mR);
}

// Legend on right map
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
new LegendControl({ position: "bottomright" }).addTo(mapRight);

function drawLegend(vMin, vMax, isPval, isNew) {
  let canvas = document.getElementById("legend-canvas");
  let ctx = canvas.getContext("2d");
  let pal = isPval ? PVAL_PALETTE : (isNew ? GREEN_PALETTE : VIRIDIS);
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
const histDiv = document.getElementById("hist-container");
const selInfo = document.getElementById("selection-info");
const tierBar = document.getElementById("tier-bar");
const modeBadge = document.getElementById("mode-badge");

let highlightL = null, highlightR = null;

let scatterToData = [];

function clearHighlight() {
  if (highlightL) { mapLeft.removeLayer(highlightL); highlightL = null; }
  if (highlightR) { mapRight.removeLayer(highlightR); highlightR = null; }
  selInfo.style.display = "none";
}

function getCurrentColKey() {
  let sig = document.getElementById("sig-select").value;
  let stat = document.getElementById("stat-select").value;
  return SINGLE_VALUE.has(sig) ? sig : sig + stat;
}

function isNewOnly(colKey) {
  return NEW_ONLY_COLS.has(colKey);
}

function highlightGage(dataIdx) {
  clearHighlight();
  let lat = DATA.lt[dataIdx], lon = DATA.ln[dataIdx];
  if (lat === null || lon === null) return;
  let ll = [lat, lon];
  highlightL = L.circleMarker(ll, { radius: 10, weight: 3, color: "#facc15", fillColor: "#facc15", fillOpacity: 0.3, pane: "markerPane" }).addTo(mapLeft);
  highlightR = L.circleMarker(ll, { radius: 10, weight: 3, color: "#facc15", fillColor: "#facc15", fillOpacity: 0.3, pane: "markerPane" }).addTo(mapRight);
  syncing = true;
  mapLeft.setView(ll, 8, { animate: true });
  mapRight.setView(ll, 8, { animate: true });
  setTimeout(() => { syncing = false; }, 300);

  let colKey = getCurrentColKey();
  let newOnly = isNewOnly(colKey);
  let newV = DATA.new[colKey] ? DATA.new[colKey][dataIdx] : null;

  if (newOnly) {
    selInfo.innerHTML =
      `Selected: <strong>${DATA.g[dataIdx]}</strong> &nbsp;|&nbsp; ` +
      `Value = ${newV !== null ? newV : "NA"} &nbsp; ` +
      `<span style="cursor:pointer;text-decoration:underline;" onclick="clearSelection()">[clear]</span>`;
  } else {
    let goldenV = DATA.golden[colKey] ? DATA.golden[colKey][dataIdx] : null;
    let diff = (goldenV !== null && newV !== null) ? (newV - goldenV).toPrecision(4) : "N/A";
    selInfo.innerHTML =
      `Selected: <strong>${DATA.g[dataIdx]}</strong> &nbsp;|&nbsp; ` +
      `Golden = ${goldenV !== null ? goldenV : "NA"} &nbsp;|&nbsp; ` +
      `New = ${newV !== null ? newV : "NA"} &nbsp;|&nbsp; ` +
      `Diff = ${diff} &nbsp; ` +
      `<span style="cursor:pointer;text-decoration:underline;" onclick="clearSelection()">[clear]</span>`;
  }
  selInfo.style.display = "block";
}

function clearSelection() {
  clearHighlight();
  let colors = scatterToData.map(() => "rgba(99,179,237,0.5)");
  let sizes = scatterToData.map(() => 3);
  if (scatterToData.length > 0) {
    Plotly.restyle(scatterDiv, { "marker.color": [colors], "marker.size": [sizes] }, [0]);
  }
}

function showGlobalTiers() {
  let tiers = DATA.tier_counts;
  let total = Object.values(tiers).reduce((a, b) => a + b, 0) - (tiers["N/A"] || 0);
  let parts = [
    ["Perfect", "tier-perfect"], ["Good", "tier-good"], ["Poor", "tier-poor"],
    ["Low", "tier-low"], ["Very Low", "tier-verylow"], ["Extremely Low", "tier-extremelylow"],
  ];
  let html = `<span style="color:#a0a0c0;">Common columns (${total}):</span> `;
  for (let [name, cls] of parts) {
    let c = tiers[name] || 0;
    if (c === 0) continue;
    let pct = total > 0 ? (100 * c / total).toFixed(1) : "0";
    html += `<span class="${cls}">${name}: ${c} (${pct}%)</span> `;
  }
  html += `<span style="color:#4ade80;margin-left:20px;">New-only columns: ${DATA.new_only_cols.length}</span>`;
  tierBar.innerHTML = html;
}

function updateMaps(colKey, newOnly) {
  let newVals = DATA.new[colKey];
  let goldenVals = newOnly ? null : DATA.golden[colKey];
  if (!newVals) return;

  let isPval = colKey.endsWith("_pval");

  // Compute value range from new values (and golden if comparison)
  let allVals = [];
  for (let i = 0; i < N; i++) {
    if (newVals[i] !== null) allVals.push(newVals[i]);
    if (!newOnly && goldenVals && goldenVals[i] !== null) allVals.push(goldenVals[i]);
  }
  let vMin, vMax;
  if (isPval) { vMin = 0; vMax = 1; }
  else if (allVals.length > 10) { vMin = percentile(allVals, 2); vMax = percentile(allVals, 98); }
  else if (allVals.length > 0) { vMin = Math.min(...allVals); vMax = Math.max(...allVals); }
  else { vMin = 0; vMax = 1; }
  if (vMax === vMin) vMax = vMin + 1;

  let pal = isPval ? PVAL_PALETTE : (newOnly ? GREEN_PALETTE : VIRIDIS);
  function colorFor(v) { return v === null ? "#444" : interpColor(pal, (v - vMin) / (vMax - vMin)); }

  // Update map labels
  document.getElementById("map-left-label").textContent = newOnly ? "New Benchmark (Value)" : "Golden Julia (Canonical)";
  document.getElementById("map-right-label").textContent = newOnly ? "NA Coverage (grey = NA)" : "New Benchmark";

  if (newOnly) {
    // Left map: value coloring, Right map: NA coverage
    for (let m of markersLeft) {
      let v = newVals[m._idx];
      m.setStyle({ fillColor: colorFor(v) });
      m.unbindPopup();
      m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>Value: ${v !== null ? v : "NA"}`);
    }
    for (let m of markersRight) {
      let v = newVals[m._idx];
      let color = v !== null ? "#4ade80" : "#555";
      m.setStyle({ fillColor: color });
      m.unbindPopup();
      m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>${v !== null ? "Has value" : "NA"}`);
    }
  } else {
    // Both maps: value coloring
    for (let m of markersLeft) {
      let v = goldenVals ? goldenVals[m._idx] : null;
      m.setStyle({ fillColor: colorFor(v) });
      m.unbindPopup();
      m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>Golden: ${v !== null ? v : "NA"}`);
    }
    for (let m of markersRight) {
      let v = newVals[m._idx];
      m.setStyle({ fillColor: colorFor(v) });
      m.unbindPopup();
      m.bindPopup(`<b>${DATA.g[m._idx]}</b><br>New: ${v !== null ? v : "NA"}`);
    }
  }

  drawLegend(vMin, vMax, isPval, newOnly);
}

function updateComparison(colKey) {
  let goldenVals = DATA.golden[colKey];
  let newVals = DATA.new[colKey];
  if (!goldenVals || !newVals) {
    document.getElementById("summary").innerHTML = `<em>Column "${colKey}" not found in both outputs.</em>`;
    return;
  }

  let nGolden_NA = goldenVals.filter(v => v === null).length;
  let nNew_NA = newVals.filter(v => v === null).length;
  let r2 = r2Identity(goldenVals, newVals);
  let nPaired = 0;
  for (let i = 0; i < N; i++) if (goldenVals[i] !== null && newVals[i] !== null) nPaired++;

  let r2Str = r2 !== null ? r2.toFixed(6) : "N/A";
  let { tier, cls } = classifyR2(r2);
  document.getElementById("summary").innerHTML =
    `<strong>${colKey}</strong> &mdash; ` +
    `R&sup2; = <span class="r2-val ${cls}">${r2Str}</span> ` +
    `<span class="${cls}">[${tier}]</span> &nbsp;|&nbsp; ` +
    `N paired = ${nPaired} &nbsp;|&nbsp; ` +
    `Golden NA = ${nGolden_NA}, New NA = ${nNew_NA}`;

  // Scatterplot
  let xVals = [], yVals = [], hoverText = [];
  scatterToData = [];
  for (let i = 0; i < N; i++) {
    if (goldenVals[i] !== null && newVals[i] !== null) {
      xVals.push(goldenVals[i]);
      yVals.push(newVals[i]);
      hoverText.push(DATA.g[i]);
      scatterToData.push(i);
    }
  }

  let isPval = colKey.endsWith("_pval");
  let pLo, pHi;
  if (isPval) { pLo = 0; pHi = 1; }
  else if (xVals.length > 20) {
    let combined = xVals.concat(yVals);
    pLo = percentile(combined, 1);
    pHi = percentile(combined, 99);
    let pad = (pHi - pLo) * 0.05;
    pLo -= pad; pHi += pad;
  } else if (xVals.length > 0) {
    pLo = Math.min(...xVals, ...yVals);
    pHi = Math.max(...xVals, ...yVals);
  } else { pLo = 0; pHi = 1; }

  Plotly.react(scatterDiv, [{
    x: xVals, y: yVals, text: hoverText,
    mode: "markers", type: "scattergl",
    marker: { size: 3, color: "rgba(99,179,237,0.5)" },
    hovertemplate: "<b>%{text}</b><br>Golden: %{x:.4f}<br>New: %{y:.4f}<extra></extra>"
  }], {
    shapes: [{ type: "line", x0: pLo, y0: pLo, x1: pHi, y1: pHi,
      line: { color: "rgba(250,200,50,0.5)", width: 1.5, dash: "dash" } }],
    xaxis: { title: "Golden Julia", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: "New Benchmark", range: [pLo, pHi], gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `${colKey}  (R\u00b2 = ${r2Str},  N = ${nPaired})`, font: { color: "#e0e0e0", size: 14 } },
    plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 50, l: 60, r: 20 }, font: { color: "#a0a0c0" }
  }, { responsive: true });

  // Histogram of differences
  let diffs = [];
  for (let i = 0; i < xVals.length; i++) diffs.push(yVals[i] - xVals[i]);

  let maxAbsDiff = diffs.length > 0 ? Math.max(...diffs.map(Math.abs)) : 0;
  let diffLabel = maxAbsDiff < 1e-10 ? "Max |diff| < 1e-10 (identical)" :
                  `Max |diff| = ${maxAbsDiff.toExponential(2)}`;

  Plotly.react(histDiv, [{
    x: diffs, type: "histogram",
    marker: { color: "rgba(99,179,237,0.6)", line: { color: "rgba(99,179,237,0.9)", width: 1 } },
    nbinsx: Math.min(50, Math.max(10, Math.ceil(diffs.length / 20))),
  }], {
    xaxis: { title: "Difference (New - Golden)", gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: "Count", gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `Differences: ${diffLabel}`, font: { color: "#e0e0e0", size: 14 } },
    plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 50, l: 60, r: 20 }, font: { color: "#a0a0c0" }
  }, { responsive: true });
}

function updateDistribution(colKey) {
  let newVals = DATA.new[colKey];
  if (!newVals) {
    document.getElementById("summary").innerHTML = `<em>Column "${colKey}" not found.</em>`;
    return;
  }

  let nonNA = newVals.filter(v => v !== null);
  let nNA = newVals.filter(v => v === null).length;

  let stats = {};
  if (nonNA.length > 0) {
    let sorted = nonNA.slice().sort((a, b) => a - b);
    stats.min = sorted[0];
    stats.max = sorted[sorted.length - 1];
    stats.median = percentile(nonNA, 50);
    stats.q25 = percentile(nonNA, 25);
    stats.q75 = percentile(nonNA, 75);
    stats.mean = nonNA.reduce((a, b) => a + b, 0) / nonNA.length;
  }

  document.getElementById("summary").innerHTML =
    `<strong style="color:#4ade80;">NEW</strong> <strong>${colKey}</strong> &mdash; ` +
    `Non-NA: <span style="color:#4ade80;font-weight:bold;">${nonNA.length}</span> / ${N} &nbsp;|&nbsp; ` +
    `NA: ${nNA} &nbsp;|&nbsp; ` +
    (nonNA.length > 0 ?
      `Range: [${stats.min.toPrecision(4)}, ${stats.max.toPrecision(4)}] &nbsp;|&nbsp; ` +
      `Median: ${stats.median.toPrecision(4)} &nbsp;|&nbsp; ` +
      `IQR: [${stats.q25.toPrecision(4)}, ${stats.q75.toPrecision(4)}]`
      : 'All NA');

  // Left chart: histogram
  scatterToData = [];
  Plotly.react(scatterDiv, [{
    x: nonNA, type: "histogram",
    marker: { color: "rgba(74,222,128,0.6)", line: { color: "rgba(74,222,128,0.9)", width: 1 } },
    nbinsx: Math.min(80, Math.max(15, Math.ceil(nonNA.length / 15))),
  }], {
    xaxis: { title: colKey, gridcolor: "#2a2a4a", color: "#a0a0c0" },
    yaxis: { title: "Count", gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `Distribution (N = ${nonNA.length})`, font: { color: "#e0e0e0", size: 14 } },
    plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 50, l: 60, r: 20 }, font: { color: "#a0a0c0" }
  }, { responsive: true });

  // Right chart: box + violin
  Plotly.react(histDiv, [{
    y: nonNA, type: "violin", box: { visible: true },
    meanline: { visible: true },
    fillcolor: "rgba(74,222,128,0.3)",
    line: { color: "rgba(74,222,128,0.8)" },
    marker: { color: "rgba(74,222,128,0.4)", size: 2 },
    points: nonNA.length < 500 ? "all" : "outliers",
    jitter: 0.3,
  }], {
    yaxis: { title: colKey, gridcolor: "#2a2a4a", color: "#a0a0c0" },
    title: { text: `Violin + Box Plot`, font: { color: "#e0e0e0", size: 14 } },
    plot_bgcolor: "#1a1a2e", paper_bgcolor: "#1a1a2e",
    margin: { t: 40, b: 30, l: 70, r: 20 }, font: { color: "#a0a0c0" },
    showlegend: false,
  }, { responsive: true });
}

function update() {
  clearHighlight();
  let sig = document.getElementById("sig-select").value;
  let stat = document.getElementById("stat-select").value;
  let statSel = document.getElementById("stat-select");
  let isSingle = SINGLE_VALUE.has(sig);
  statSel.disabled = isSingle;
  let colKey = isSingle ? sig : sig + stat;

  let newOnly = isNewOnly(colKey);

  if (newOnly) {
    modeBadge.textContent = "NEW SIGNATURE";
    modeBadge.className = "mode-distribution";
    updateDistribution(colKey);
  } else {
    modeBadge.textContent = "COMPARISON";
    modeBadge.className = "mode-comparison";
    updateComparison(colKey);
  }

  updateMaps(colKey, newOnly);
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
    if not os.path.exists(GOLDEN_PATH):
        print(f"ERROR: Golden Julia file not found: {GOLDEN_PATH}")
        sys.exit(1)
    if not os.path.exists(NEW_PATH):
        print(f"ERROR: New benchmark output not found: {NEW_PATH}")
        sys.exit(1)

    data = build_data()
    build_html(data)
    print("Done.")

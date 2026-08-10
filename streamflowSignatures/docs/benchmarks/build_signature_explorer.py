#!/usr/bin/env python
"""Build a self-contained signature map explorer (Leaflet + Canvas HTML).

Maps every signature base x statistic combination from a benchmark summary CSV:
gage points colored by the selected variable, robust quantile color scaling,
per-gage click panel with the full 8-stat row, gage search. Simpler sibling of
data_out/streamflow_explorer.html, focused purely on signatures + stats.

Usage:
    python build_signature_explorer.py --csv <signatures.csv> --out <explorer.html>
        [--annual <signatures_annual.parquet>] [--title "..."] [--limit N]

Needs internet at view time for the Leaflet/CARTO CDN (same as the other explorers);
points and summary data are fully inlined.

With --annual, per-signature sidecar files (<out-stem>_annual/<signature>.js, ~90 files)
are written next to the HTML and lazily loaded on gage click to draw the annual
time-series plot (works from file:// via <script> injection — fetch() would be blocked).
The sidecar folder must travel with the HTML.
"""
import argparse
import json
import math
import os
from datetime import date

import numpy as np
import pandas as pd

STAT_SUFFIXES = ["_mean", "_median", "_senn_slp", "_linear_slp",
                 "_spearman_rho", "_spearman_pval", "_mk_rho", "_mk_pval",
                 "_pettitt_cp_year", "_pettitt_pval", "_pettitt_pre_mean",
                 "_pettitt_post_mean", "_pettitt_delta_mean", "_pettitt_pct_change",
                 "_pettitt_pre_mk_pval", "_pettitt_post_mk_pval"]
STAT_LABELS = {"_mean": "Mean", "_median": "Median", "_senn_slp": "Theil-Sen slope",
               "_linear_slp": "Linear slope", "_spearman_rho": "Spearman rho",
               "_spearman_pval": "Spearman p-value", "_mk_rho": "Mann-Kendall tau",
               "_mk_pval": "Mann-Kendall p-value",
               "_pettitt_cp_year": "Pettitt changepoint year",
               "_pettitt_pval": "Pettitt p-value",
               "_pettitt_pre_mean": "Pettitt pre-CP mean",
               "_pettitt_post_mean": "Pettitt post-CP mean",
               "_pettitt_delta_mean": "Pettitt delta mean (post-pre)",
               "_pettitt_pct_change": "Pettitt % change",
               "_pettitt_pre_mk_pval": "Pettitt pre-CP MK p-value",
               "_pettitt_post_mk_pval": "Pettitt post-CP MK p-value"}
# Documented per-gage scalar signatures (no 8-stat suffixes)
SCALARS = ["elasticity_static", "recession_alpha_point_cloud_linear_reservoir",
           "log_a_seasonality_amplitude_all", "log_a_seasonality_amplitude_first_half",
           "log_a_seasonality_amplitude_last_half", "log_a_seasonality_minimum_all",
           "log_a_seasonality_minimum_first_half", "log_a_seasonality_minimum_last_half",
           "runoff_ratio_high_count", "elasticity_years_total", "elasticity_years_low_ppt",
           "ice_affected_days_total", "season_excluded_years_winter",
           "season_excluded_years_spring", "season_excluded_years_summer",
           "season_excluded_years_fall",
           # drought family thresholds (per-gage scalars, July 2026)
           "drought_threshold_fixed_p2", "drought_threshold_fixed_p5",
           "drought_threshold_fixed_p10", "drought_threshold_fixed_p20",
           "drought_threshold_fixed_p30"]


def category_of(base):
    b = base.lower()
    # drought first: `dur_*`/`d*` rules below would otherwise not match, but the
    # explicit rule keeps the family together regardless of future renames
    if b.startswith("drought_"):
        return "Streamflow drought"
    if b.startswith(("swe_", "snow_", "melt_")) or b == "ssm":
        return "Snow"
    if b.startswith("q") and (b[1:].replace("_q10", "").replace("95_q10", "").isdigit()
                              or b in ("q95_q10",)):
        return "Flow percentiles"
    if b in ("qann", "qwin", "qspr", "qsum", "qfal"):
        return "Flow volumes"
    if b.startswith("fdc"):
        return "Flow duration curve"
    if b.startswith("bfi"):
        return "Baseflow"
    if b.startswith(("log_a", "b_point", "b_event", "concavity", "alpha_linear",
                     "n_recession")):
        return "Recession"
    if b.startswith(("n_high_pulses", "n_low_pulses", "dur_high", "dur_low", "tqmean",
                     "flow_reversals")):
        return "Pulses"
    if b.startswith("flashiness"):
        return "Flashiness"
    if b.startswith("d") and ("_day" in b or b in ("d25_to_d75", "dmax")):
        return "Flow timing"
    if "runoff_ratio" in b:
        return "Runoff ratios"
    if b.startswith("elasticity"):
        return "Elasticity"
    if b.startswith("qp_"):
        return "Q-P seasonality"
    if b.startswith("avg_storage"):
        return "Storage"
    if b.startswith("negative"):
        return "Negative flow"
    return "Other"


def sig4(x):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return None
    if x == 0:
        return 0
    return float(f"{x:.4g}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--annual", default=None,
                    help="signatures_annual.parquet -> per-signature sidecar files + time-series plot")
    ap.add_argument("--title", default="Streamflow Signature Explorer")
    ap.add_argument("--limit", type=int, default=None, help="debug: first N gages")
    args = ap.parse_args()

    df = pd.read_csv(args.csv, dtype={"gage_id": str}, low_memory=False)
    if args.limit:
        df = df.head(args.limit)
    df = df.dropna(subset=["latitude", "longitude"])

    bases = sorted({c[:-5] for c in df.columns
                    if c.endswith("_mean") and "_pettitt_" not in c})
    scalars = [s for s in SCALARS if s in df.columns]

    cols = {}
    for b in bases:
        for s in STAT_SUFFIXES:
            c = b + s
            if c in df.columns:
                vals = pd.to_numeric(df[c], errors="coerce")
                cols[c] = [sig4(v) for v in vals.to_numpy(dtype=float)]
    for s in scalars:
        vals = pd.to_numeric(df[s], errors="coerce")
        cols[s] = [sig4(v) for v in vals.to_numpy(dtype=float)]
    # Pettitt changepoint year + p-value per base (plot markers; not in the pickers)
    for b in bases:
        for suf in ("_pettitt_cp_year", "_pettitt_pval"):
            c = b + suf
            if c in df.columns:
                cols[c] = [sig4(v) for v in
                           pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)]

    # Annual-values sidecar files (lazy-loaded time series)
    annual_meta = {}
    if args.annual:
        ann = pd.read_parquet(args.annual,
                              columns=["gage_id", "signature", "water_year", "value"])
        ann["gage_id"] = ann["gage_id"].astype(str)
        keep_ids = set(df["gage_id"].astype(str))
        ann = ann[ann["gage_id"].isin(keep_ids)]
        y0, y1 = int(ann["water_year"].min()), int(ann["water_year"].max())
        span = y1 - y0 + 1
        out_dir_name = os.path.splitext(os.path.basename(args.out))[0] + "_annual"
        out_dir = os.path.join(os.path.dirname(os.path.abspath(args.out)), out_dir_name)
        os.makedirs(out_dir, exist_ok=True)
        n_files = 0
        for sig, sub in ann.groupby("signature", sort=True):
            gmap = {}
            for gid, g in sub.groupby("gage_id", sort=False):
                arr = [None] * span
                for yy, vv in zip(g["water_year"].to_numpy(),
                                  g["value"].to_numpy(dtype=float)):
                    arr[int(yy) - y0] = sig4(vv)
                gmap[gid] = arr
            blob = json.dumps(gmap, separators=(",", ":"), allow_nan=False)
            with open(os.path.join(out_dir, f"{sig}.js"), "w", encoding="utf-8") as f:
                f.write(f"ANNUAL_CB({json.dumps(sig)},{y0},{blob});")
            n_files += 1
        annual_meta = {"dir": out_dir_name, "y0": y0, "y1": y1}
        print(f"Wrote {n_files} annual sidecar files ({y0}-{y1}) to {out_dir}")

    cats = {}
    for b in bases:
        cats.setdefault(category_of(b), []).append(b)

    payload = {
        "meta": {
            "title": args.title,
            "source": os.path.basename(args.csv),
            "generated": str(date.today()),
            "n_gages": int(len(df)),
            "stats": [[s, STAT_LABELS[s]] for s in STAT_SUFFIXES],
            "categories": cats,
            "scalars": scalars,
            "annual": annual_meta,
        },
        "id": df["gage_id"].astype(str).tolist(),
        "lat": [sig4(v) for v in df["latitude"].to_numpy(dtype=float)],
        "lon": [sig4(v) for v in df["longitude"].to_numpy(dtype=float)],
        "nyr": pd.to_numeric(df["num_water_years"], errors="coerce")
                 .fillna(0).astype(int).tolist(),
        "wy0": pd.to_numeric(df["start_water_year"], errors="coerce")
                 .fillna(0).astype(int).tolist(),
        "wy1": pd.to_numeric(df["end_water_year"], errors="coerce")
                 .fillna(0).astype(int).tolist(),
        "anorm": [str(v).upper() != "FALSE" for v in
                  df.get("area_normalized", pd.Series(True, index=df.index))],
        "cols": cols,
    }
    data_json = json.dumps(payload, separators=(",", ":"), allow_nan=False)

    html = TEMPLATE.replace("__TITLE__", args.title).replace("__DATA__", data_json)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Wrote {args.out}: {os.path.getsize(args.out)/1e6:.1f} MB, "
          f"{len(df)} gages, {len(bases)} bases x {len(STAT_SUFFIXES)} stats "
          f"+ {len(scalars)} scalars = {len(cols)} mapped variables")


TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css">
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  html,body{margin:0;height:100%;font-family:system-ui,Segoe UI,Arial,sans-serif}
  #map{position:absolute;top:0;bottom:0;left:0;right:0}
  .panel{position:absolute;z-index:1000;background:rgba(255,255,255,.96);
         border-radius:8px;box-shadow:0 1px 8px rgba(0,0,0,.25);padding:10px 12px}
  #controls{top:10px;left:10px;width:295px}
  #controls h1{font-size:15px;margin:0 0 6px}
  #controls .sub{font-size:11px;color:#666;margin-bottom:8px}
  #controls label{font-size:11px;color:#444;display:block;margin-top:6px}
  #controls select,#controls input[type=text]{width:100%;font-size:12px;padding:3px}
  #controls .row{display:flex;gap:6px;align-items:center;margin-top:6px}
  #controls .row label{margin:0}
  #legend{bottom:18px;left:10px;width:430px;font-size:11px}
  #legend .bar{height:12px;border-radius:3px;margin:4px 0}
  #legend .lab{display:flex;justify-content:space-between;color:#333}
  #hist{display:block;width:100%;height:110px;margin-top:6px;cursor:crosshair}
  #histinfo{font-size:10px;color:#666;min-height:13px}
  #tsplot{bottom:18px;left:450px;width:640px;display:none}
  #tsplot h3{font-size:12px;margin:0 0 2px}
  #tsplot .sub{font-size:10px;color:#666;margin-bottom:4px}
  #ts{display:block;width:100%;height:170px;cursor:crosshair}
  #tsinfo{font-size:10px;color:#666;min-height:13px}
  #tsplot .close{float:right;cursor:pointer;color:#888;font-size:14px}
  #info{top:10px;right:10px;width:305px;max-height:82vh;overflow-y:auto;display:none}
  #info h2{font-size:13px;margin:0 0 4px}
  #info table{width:100%;border-collapse:collapse;font-size:11px}
  #info td{padding:2px 4px;border-bottom:1px solid #eee}
  #info td:last-child{text-align:right;font-variant-numeric:tabular-nums}
  #info .close{float:right;cursor:pointer;color:#888;font-size:14px}
  #stats-line{font-size:11px;color:#333;margin-top:6px}
  .hl{stroke:#ff2d00 !important}
</style>
</head>
<body>
<div id="map"></div>
<div class="panel" id="controls">
  <h1>__TITLE__</h1>
  <div class="sub" id="subtitle"></div>
  <label>Signature</label>
  <select id="base"></select>
  <label>Statistic</label>
  <select id="stat"></select>
  <div class="row">
    <input type="checkbox" id="clip" checked>
    <label for="clip">Robust color scale (p5&ndash;p95)</label>
  </div>
  <div class="row">
    <input type="checkbox" id="hideNA">
    <label for="hideNA">Hide NA gages</label>
  </div>
  <label>Find gage</label>
  <input type="text" id="search" placeholder="gage_id, Enter to zoom">
  <div id="stats-line"></div>
</div>
<div class="panel" id="legend">
  <div id="legvar" style="font-weight:600"></div>
  <canvas id="hist" width="406" height="110"></canvas>
  <div id="histinfo"></div>
  <div class="bar" id="legbar"></div>
  <div class="lab"><span id="legmin"></span><span id="legmid"></span><span id="legmax"></span></div>
  <div style="color:#777;margin-top:3px">grey = NA &nbsp;&middot;&nbsp; dashed ring = area_normalized FALSE (raw m&sup3;/s) &nbsp;&middot;&nbsp; grey end bars = values beyond the p5&ndash;p95 clip &nbsp;&middot;&nbsp; dashed line = median</div>
</div>
<div class="panel" id="tsplot">
  <span class="close" onclick="document.getElementById('tsplot').style.display='none'">&#10005;</span>
  <h3 id="ts-title"></h3>
  <div class="sub" id="ts-sub"></div>
  <canvas id="ts" width="616" height="170"></canvas>
  <div id="tsinfo"></div>
</div>
<div class="panel" id="info">
  <span class="close" onclick="document.getElementById('info').style.display='none'">&#10005;</span>
  <h2 id="info-title"></h2>
  <div id="info-meta" style="font-size:11px;color:#555;margin-bottom:6px"></div>
  <table id="info-table"></table>
</div>
<script>
const D = __DATA__;
const VIRIDIS = ["#440154","#46085c","#471063","#481769","#481d6f","#482475","#472a7a","#46307e","#453781","#433d84","#414287","#3f4889","#3d4e8a","#3a538b","#38598c","#355e8d","#33638d","#31688e","#2e6d8e","#2c718e","#2a768e","#297b8e","#27808e","#25848e","#23898e","#218e8d","#20928c","#1f978b","#1e9c89","#1fa188","#21a585","#24aa83","#28ae80","#2eb37c","#35b779","#3dbc74","#46c06f","#50c46a","#5ac864","#65cb5e","#70cf57","#7cd250","#89d548","#95d840","#a2da37","#b0dd2f","#bddf26","#cae11f","#d8e219","#e5e419","#f1e51d","#fde725"];
function colorFor(t){const i=Math.max(0,Math.min(VIRIDIS.length-1,Math.floor(t*(VIRIDIS.length-1))));return VIRIDIS[i];}
function quantile(sorted,q){if(!sorted.length)return NaN;const p=(sorted.length-1)*q,lo=Math.floor(p),hi=Math.ceil(p);return sorted[lo]+(sorted[hi]-sorted[lo])*(p-lo);}
function fmt(v){if(v===null||v===undefined||Number.isNaN(v))return "NA";const a=Math.abs(v);if(a!==0&&(a<0.001||a>=100000))return v.toExponential(2);return (Math.round(v*10000)/10000).toString();}

const map=L.map('map',{preferCanvas:true}).setView([48,-98],4);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
  {attribution:'&copy; OpenStreetMap &copy; CARTO',maxZoom:12}).addTo(map);
const renderer=L.canvas({padding:0.3});

document.getElementById('subtitle').textContent =
  D.meta.source+" — "+D.meta.n_gages+" gages — generated "+D.meta.generated;

const baseSel=document.getElementById('base'),statSel=document.getElementById('stat');
for(const [cat,bases] of Object.entries(D.meta.categories)){
  const og=document.createElement('optgroup');og.label=cat;
  for(const b of bases){const o=document.createElement('option');o.value=b;o.textContent=b;og.appendChild(o);}
  baseSel.appendChild(og);
}
if(D.meta.scalars.length){
  const og=document.createElement('optgroup');og.label="Per-gage scalars";
  for(const s of D.meta.scalars){const o=document.createElement('option');o.value="SCALAR:"+s;o.textContent=s;og.appendChild(o);}
  baseSel.appendChild(og);
}
for(const [suf,lab] of D.meta.stats){const o=document.createElement('option');o.value=suf;o.textContent=lab;statSel.appendChild(o);}

const markers=[];let curIdx=-1;
for(let i=0;i<D.id.length;i++){
  const m=L.circleMarker([D.lat[i],D.lon[i]],{renderer,radius:4,weight:D.anorm[i]?0.5:1.5,
    color:D.anorm[i]?"#333":"#d62728",dashArray:D.anorm[i]?null:"2,2",fillOpacity:0.85,fillColor:"#bbb"});
  m.on('click',()=>showInfo(i));
  m.addTo(map);markers.push(m);
}

function currentCol(){
  const b=baseSel.value;
  if(b.startsWith("SCALAR:"))return b.slice(7);
  return b+statSel.value;
}
function recolor(){
  const col=currentCol();const vals=D.cols[col]||[];
  const nonNull=vals.filter(v=>v!==null).sort((a,b)=>a-b);
  const clip=document.getElementById('clip').checked;
  const hideNA=document.getElementById('hideNA').checked;
  let lo=nonNull.length?nonNull[0]:0,hi=nonNull.length?nonNull[nonNull.length-1]:1;
  if(clip&&nonNull.length>20){lo=quantile(nonNull,0.05);hi=quantile(nonNull,0.95);}
  if(hi===lo)hi=lo+1e-12;
  for(let i=0;i<markers.length;i++){
    const v=vals[i];
    if(v===null){markers[i].setStyle({fillColor:"#c8c8c8",fillOpacity:hideNA?0:0.35,opacity:hideNA?0:1});}
    else{const t=Math.max(0,Math.min(1,(v-lo)/(hi-lo)));
      markers[i].setStyle({fillColor:colorFor(t),fillOpacity:0.85,opacity:1});}
  }
  document.getElementById('legvar').textContent=col;
  document.getElementById('legbar').style.background=
    "linear-gradient(to right,"+[0,.25,.5,.75,1].map(t=>colorFor(t)).join(",")+")";
  document.getElementById('legmin').textContent=fmt(lo)+(clip?" (p5)":"");
  document.getElementById('legmid').textContent=fmt((lo+hi)/2);
  document.getElementById('legmax').textContent=fmt(hi)+(clip?" (p95)":"");
  document.getElementById('stats-line').textContent=
    nonNull.length+" of "+vals.length+" gages non-NA — median "+fmt(quantile(nonNull,0.5));
  statSel.disabled=baseSel.value.startsWith("SCALAR:");
  drawHist(nonNull,lo,hi);
  if(curIdx>=0)showInfo(curIdx);
}

// ---- histogram of the selected variable (bars viridis-colored by bin position,
// grey end bars = values outside the clip domain, dashed line = median) ----
let histBins=null;
function drawHist(nonNull,lo,hi){
  const cv=document.getElementById('hist'),ctx=cv.getContext('2d');
  const W=cv.width,H=cv.height,padB=14,padT=12;
  ctx.clearRect(0,0,W,H);
  histBins=null;
  if(!nonNull.length){ctx.fillStyle="#999";ctx.font="11px sans-serif";
    ctx.fillText("no non-NA values",8,H/2);return;}
  const NB=40,edgeW=9,plotX=edgeW+2,plotW=W-2*(edgeW+2);
  const counts=new Array(NB).fill(0);let under=0,over=0;
  for(const v of nonNull){
    if(v<lo){under++;continue;}
    if(v>hi){over++;continue;}
    let k=Math.floor((v-lo)/(hi-lo)*NB);if(k>=NB)k=NB-1;if(k<0)k=0;
    counts[k]++;
  }
  const maxC=Math.max(1,...counts,under,over);
  const y=(c)=>H-padB-(H-padB-padT)*(c/maxC);
  // underflow / overflow edge bars (grey)
  ctx.fillStyle="#9e9e9e";
  if(under>0)ctx.fillRect(0,y(under),edgeW,H-padB-y(under));
  if(over>0)ctx.fillRect(W-edgeW,y(over),edgeW,H-padB-y(over));
  // main bins, viridis-colored at bin centers
  const bw=plotW/NB;
  for(let k=0;k<NB;k++){
    ctx.fillStyle=colorFor((k+0.5)/NB);
    ctx.fillRect(plotX+k*bw,y(counts[k]),Math.max(1,bw-1),H-padB-y(counts[k]));
  }
  // median marker (within clip domain)
  const med=quantile(nonNull,0.5);
  if(med>=lo&&med<=hi){
    const mx=plotX+(med-lo)/(hi-lo)*plotW;
    ctx.strokeStyle="#222";ctx.setLineDash([3,3]);ctx.beginPath();
    ctx.moveTo(mx,padT);ctx.lineTo(mx,H-padB);ctx.stroke();ctx.setLineDash([]);
  }
  // baseline + count label
  ctx.strokeStyle="#ccc";ctx.beginPath();ctx.moveTo(0,H-padB+0.5);ctx.lineTo(W,H-padB+0.5);ctx.stroke();
  ctx.fillStyle="#555";ctx.font="10px sans-serif";
  ctx.fillText("max bin: "+maxC+" gages",8,10);
  histBins={lo,hi,NB,counts,under,over,edgeW,plotX,plotW};
  document.getElementById('histinfo').textContent="hover bars for bin ranges";
}
document.getElementById('hist').addEventListener('mousemove',e=>{
  if(!histBins)return;
  const cv=document.getElementById('hist');
  const x=(e.offsetX/cv.clientWidth)*cv.width;
  const {lo,hi,NB,counts,under,over,edgeW,plotX,plotW}=histBins;
  const info=document.getElementById('histinfo');
  if(x<edgeW+1){info.textContent="< "+fmt(lo)+" (below clip): "+under+" gages";return;}
  if(x>cv.width-edgeW-1){info.textContent="> "+fmt(hi)+" (above clip): "+over+" gages";return;}
  let k=Math.floor((x-plotX)/plotW*NB);k=Math.max(0,Math.min(NB-1,k));
  const a=lo+(hi-lo)*k/NB,b=lo+(hi-lo)*(k+1)/NB;
  info.textContent="["+fmt(a)+", "+fmt(b)+"): "+counts[k]+" gages";
});
document.getElementById('hist').addEventListener('mouseleave',()=>{
  if(histBins)document.getElementById('histinfo').textContent="hover bars for bin ranges";
});
function showInfo(i){
  curIdx=i;
  const b=baseSel.value;const panel=document.getElementById('info');
  document.getElementById('info-title').textContent="Gage "+D.id[i];
  document.getElementById('info-meta').innerHTML=
    D.nyr[i]+" water years ("+D.wy0[i]+"–"+D.wy1[i]+")"+
    (D.anorm[i]?"":" — <b style='color:#d62728'>area_normalized = FALSE (raw m&sup3;/s)</b>");
  const t=document.getElementById('info-table');t.innerHTML="";
  const addRow=(k,v)=>{const tr=t.insertRow();tr.insertCell().textContent=k;tr.insertCell().textContent=fmt(v);};
  if(b.startsWith("SCALAR:")){
    const s=b.slice(7);addRow(s,(D.cols[s]||[])[i]);
  }else{
    for(const [suf,lab] of D.meta.stats){
      const col=b+suf;if(D.cols[col])addRow(b+suf,(D.cols[col])[i]);
    }
  }
  const hr=t.insertRow();const c=hr.insertCell();c.colSpan=2;c.style.paddingTop="6px";
  c.style.color="#777";c.textContent="Per-gage scalars:";
  for(const s of D.meta.scalars){addRow(s,(D.cols[s]||[])[i]);}
  panel.style.display="block";
  markers.forEach(m=>m.setStyle({weight:D.anorm[markers.indexOf(m)]?0.5:1.5}));
  requestTS();
}

// ---- annual time series (lazy sidecar loading; plan: sidecar folder next to HTML) ----
const annualCache={};let pendingSig=null;
window.ANNUAL_CB=function(sig,y0,data){annualCache[sig]={y0,data};
  if(pendingSig===sig){pendingSig=null;drawTS();}};
function requestTS(){
  const A=D.meta.annual;const panel=document.getElementById('tsplot');
  if(!A||!A.dir||curIdx<0)return;
  const b=baseSel.value;
  if(b.startsWith("SCALAR:")){
    panel.style.display="block";
    document.getElementById('ts-title').textContent=b.slice(7)+" — gage "+D.id[curIdx];
    document.getElementById('ts-sub').textContent="per-gage scalar — no annual series";
    const ctx=document.getElementById('ts').getContext('2d');
    ctx.clearRect(0,0,616,170);tsData=null;return;
  }
  if(annualCache[b]){drawTS();return;}
  panel.style.display="block";
  document.getElementById('ts-title').textContent=b+" — gage "+D.id[curIdx];
  document.getElementById('ts-sub').textContent="loading annual values…";
  pendingSig=b;
  const s=document.createElement('script');
  s.src=A.dir+"/"+b+".js";
  s.onerror=()=>{document.getElementById('ts-sub').textContent=
    "annual sidecar not found ("+A.dir+"/"+b+".js — keep the folder next to this HTML)";};
  document.head.appendChild(s);
}
let tsData=null;
function drawTS(){
  const A=D.meta.annual;if(!A||curIdx<0)return;
  const b=baseSel.value;if(b.startsWith("SCALAR:"))return;
  const entry=annualCache[b];if(!entry)return;
  const series=entry.data[D.id[curIdx]];
  const panel=document.getElementById('tsplot');panel.style.display="block";
  const cv=document.getElementById('ts'),ctx=cv.getContext('2d');
  const W=cv.width,H=cv.height,mL=52,mR=10,mT=8,mB=22;
  ctx.clearRect(0,0,W,H);tsData=null;
  document.getElementById('ts-title').textContent=b+" — gage "+D.id[curIdx];
  const slope=(D.cols[b+"_senn_slp"]||[])[curIdx];
  const mkp=(D.cols[b+"_mk_pval"]||[])[curIdx];
  const cpy=(D.cols[b+"_pettitt_cp_year"]||[])[curIdx];
  const cpp=(D.cols[b+"_pettitt_pval"]||[])[curIdx];
  document.getElementById('ts-sub').textContent=
    "Theil-Sen slope "+fmt(slope)+"/yr (MK p "+fmt(mkp)+")"+
    (cpy!==null&&cpy!==undefined?" — Pettitt CP "+Math.round(cpy)+" (p "+fmt(cpp)+")":"");
  if(!series||series.every(v=>v===null)){
    ctx.fillStyle="#999";ctx.font="11px sans-serif";
    ctx.fillText("no annual values for this gage",mL,H/2);return;}
  const y0=entry.y0,n=series.length;
  const pts=[];for(let k=0;k<n;k++)if(series[k]!==null)pts.push([y0+k,series[k]]);
  let vmin=Math.min(...pts.map(p=>p[1])),vmax=Math.max(...pts.map(p=>p[1]));
  if(vmin===vmax){vmin-=1;vmax+=1;}
  const pad=(vmax-vmin)*0.07;vmin-=pad;vmax+=pad;
  const x=(yr)=>mL+(yr-(y0-0.5))/((y0+n-0.5)-(y0-0.5))*(W-mL-mR);
  const y=(v)=>mT+(1-(v-vmin)/(vmax-vmin))*(H-mT-mB);
  // axes
  ctx.strokeStyle="#ddd";ctx.fillStyle="#666";ctx.font="9px sans-serif";
  for(let t=0;t<4;t++){const v=vmin+(vmax-vmin)*t/3;
    ctx.beginPath();ctx.moveTo(mL,y(v));ctx.lineTo(W-mR,y(v));ctx.stroke();
    ctx.fillText(fmt(v),4,y(v)+3);}
  const step=n>25?5:(n>12?2:1);
  for(let yr=Math.ceil(y0/step)*step;yr<=y0+n-1;yr+=step){
    ctx.fillText(yr,x(yr)-11,H-6);}
  // Pettitt changepoint marker
  if(cpy!==null&&cpy!==undefined){
    const sig05=cpp!==null&&cpp!==undefined&&cpp<0.05;
    ctx.strokeStyle=sig05?"#e07b00":"#aaa";ctx.setLineDash([4,3]);
    ctx.beginPath();ctx.moveTo(x(cpy),mT);ctx.lineTo(x(cpy),H-mB);ctx.stroke();
    ctx.setLineDash([]);ctx.fillStyle=sig05?"#e07b00":"#999";
    ctx.fillText("CP "+Math.round(cpy),Math.min(x(cpy)+3,W-58),mT+8);
  }
  // Theil-Sen trend line (intercept = median(y - slope*x) over the plotted points)
  if(slope!==null&&slope!==undefined&&pts.length>=2){
    const resid=pts.map(p=>p[1]-slope*p[0]).sort((a,b)=>a-b);
    const b0=quantile(resid,0.5);
    const xa=pts[0][0],xb=pts[pts.length-1][0];
    ctx.strokeStyle="#d62728";ctx.lineWidth=1.5;
    ctx.beginPath();ctx.moveTo(x(xa),y(b0+slope*xa));ctx.lineTo(x(xb),y(b0+slope*xb));
    ctx.stroke();ctx.lineWidth=1;
  }
  // series: line segments between consecutive years, points everywhere
  ctx.strokeStyle="#31688e";ctx.fillStyle="#31688e";
  ctx.beginPath();let started=false;
  for(let k=0;k<n;k++){
    if(series[k]===null){started=false;continue;}
    const px=x(y0+k),py=y(series[k]);
    if(started)ctx.lineTo(px,py);else{ctx.moveTo(px,py);started=true;}
  }
  ctx.stroke();
  for(const [yr,v] of pts){ctx.beginPath();ctx.arc(x(yr),y(v),2.6,0,2*Math.PI);ctx.fill();}
  tsData={y0,n,series,x,vmin,vmax};
  document.getElementById('tsinfo').textContent=
    pts.length+" annual values ("+pts[0][0]+"–"+pts[pts.length-1][0]+"); gaps = NA years";
}
document.getElementById('ts').addEventListener('mousemove',e=>{
  if(!tsData)return;
  const cv=document.getElementById('ts');
  const px=(e.offsetX/cv.clientWidth)*cv.width;
  const yr=Math.round((y0=>{const {n}=tsData;
    return (px-52)/(cv.width-52-10)*((y0+n-0.5)-(y0-0.5))+(y0-0.5);})(tsData.y0));
  const k=yr-tsData.y0;
  if(k>=0&&k<tsData.n){
    const v=tsData.series[k];
    document.getElementById('tsinfo').textContent="WY "+yr+": "+(v===null?"NA":fmt(v));
  }
});
document.getElementById('search').addEventListener('keydown',e=>{
  if(e.key!=="Enter")return;
  const q=e.target.value.trim();
  const i=D.id.indexOf(q)>=0?D.id.indexOf(q):D.id.findIndex(g=>g.includes(q));
  if(i>=0){map.setView([D.lat[i],D.lon[i]],9);markers[i].setStyle({weight:3,color:"#ff2d00"});showInfo(i);}
});
baseSel.addEventListener('change',recolor);
statSel.addEventListener('change',recolor);
document.getElementById('clip').addEventListener('change',recolor);
document.getElementById('hideNA').addEventListener('change',recolor);
baseSel.value="Qann";recolor();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()

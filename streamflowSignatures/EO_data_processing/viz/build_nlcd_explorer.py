"""Self-contained HTML QA/QC explorer for the per-watershed Annual NLCD product (CONUS 1985-2025).

Parallels build_lulc_explorer.py but for NLCD's SINGLE 30 m land-cover scheme (16 classes) + a
continuous fractional-impervious layer. Map: ~6,119 in-CONUS gage points colored by a selectable
variable; click a gage -> stacked-area land-cover composition 1985-2025 (official NLCD palette) with
an overlaid impervious-mean trend line; the annual-update seam (2024-2025) is shaded. Per-year hover.

Design reconciled with the Codex review (EO_data_processing/README_NLCD.md):
  - deltas are SMOOTHED endpoint differences mean(1985-87) vs mean(2023-25), LABELLED as an endpoint
    difference (not a trend) with the annual-update seam disclosed;
  - developed% (classes 21-24) and impervious% are shown as related-but-NON-ADDITIVE (separate);
  - a shrub<->grass / largest-one-year-change "volatility" variable surfaces the documented annual
    model-noise artifact;
  - partial-coverage / low-support basins get a visible badge; the 45 excluded Alaska gages are
    disclosed (count + listed in the glossary), never shown as zero-valued points;
  - normalized Shannon diversity (H/log 16); compact 0.1%-integer series encoding.

Usage: /home/sagemaker-user/.conda/envs/geo/bin/python EO_data_processing/viz/build_nlcd_explorer.py
Output: data_out/eo_nlcd/nlcd_out/nlcd_explorer.html
"""
import os, json
from collections import Counter
import numpy as np, pandas as pd, geopandas as gpd

ROOT = "/home/sagemaker-user/streamflowSignatures/data_out/eo_nlcd"
NLCD = f"{ROOT}/nlcd_out/watershed_nlcd_annual_24jul2026.parquet"
LEG = "/home/sagemaker-user/streamflowSignatures/EO_data_processing/eo_processing/nlcd_legends.csv"
GEOM = f"{ROOT}/watershed_polygons_26jun2026.parquet"
OUT = f"{ROOT}/nlcd_out/nlcd_explorer.html"
SEAM_FROM = 2024  # C1V2 annual-update-path years (methodological seam)

CLASS_LABEL = {11:"Open Water",12:"Perennial Ice/Snow",21:"Developed, Open Space",
 22:"Developed, Low",23:"Developed, Medium",24:"Developed, High",31:"Barren Land",
 41:"Deciduous Forest",42:"Evergreen Forest",43:"Mixed Forest",52:"Shrub/Scrub",
 71:"Grassland/Herb.",81:"Pasture/Hay",82:"Cultivated Crops",90:"Woody Wetlands",
 95:"Emergent Herb. Wetlands"}
# aggregate groups (lump-robust) — code -> group
GROUP = {11:"water",12:"snow/ice",21:"developed",22:"developed",23:"developed",24:"developed",
 31:"barren",41:"forest",42:"forest",43:"forest",52:"shrubland",71:"grassland",
 81:"agriculture",82:"agriculture",90:"wetland",95:"wetland"}


def main():
    print("loading finalized NLCD panel + legend + geometry ...", flush=True)
    leg = pd.read_csv(LEG)
    codes = [int(c) for c in leg.code]
    cols = list(leg.output_column)
    colors = list(leg.hex_color)
    names = [CLASS_LABEL[c] for c in codes]
    df = pd.read_parquet(NLCD); df["gage_id"] = df.gage_id.astype(str)
    geo = gpd.read_parquet(GEOM); geo["gage_id"] = geo.gage_id.astype(str)
    geo = geo[["gage_id", "canon_id", "latitude", "longitude"]].drop_duplicates("gage_id")

    years = sorted(int(y) for y in df.year.unique()); nyear = len(years)
    df = df.merge(geo, on="gage_id", how="left", suffixes=("", "_g"))
    gages = sorted(df.gage_id.unique()); ngage = len(gages)
    print(f"  {len(df):,} rows | {nyear} years {years[0]}..{years[-1]} | {ngage} CONUS gages", flush=True)

    full = pd.MultiIndex.from_product([gages, years], names=["gage_id", "year"])
    d2 = df.set_index(["gage_id", "year"]).reindex(full)
    def mat(col): return d2[col].to_numpy(dtype=float).reshape(ngage, nyear)
    meta1 = df.drop_duplicates("gage_id").set_index("gage_id")

    C = {c: np.nan_to_num(mat(col)) for c, col in zip(codes, cols)}   # code -> (ngage,nyear) %
    imp = np.nan_to_num(mat("impervious_mean_pct"))
    f3, l3 = slice(0, 3), slice(nyear - 3, nyear)
    endpt = lambda M: np.nanmean(M[:, l3], axis=1) - np.nanmean(M[:, f3], axis=1)
    tmean = lambda M: np.nanmean(M, axis=1)

    # aggregate group matrices
    groups = {}
    for g in ["water", "snow/ice", "developed", "barren", "forest", "shrubland", "grassland", "agriculture", "wetland"]:
        gcodes = [c for c in codes if GROUP[c] == g]
        groups[g] = np.sum([C[c] for c in gcodes], axis=0)

    comp_mean = np.stack([tmean(C[c]) for c in codes], axis=1)          # (ngage,16) temporal-mean %
    dominant_class = [names[j] for j in np.argmax(comp_mean, axis=1)]
    gnames = list(groups)
    gmean = np.stack([tmean(groups[g]) for g in gnames], axis=1)
    dominant_group = [gnames[j] for j in np.argmax(gmean, axis=1)]
    p = np.clip(comp_mean / 100.0, 1e-12, 1); shannon = -(p * np.log(p)).sum(axis=1) / np.log(len(codes))  # normalized 0..1

    # volatility / artifact indicator: mean over years of total-variation 0.5*sum_c|dC_c| (pp/yr)
    stack = np.stack([C[c] for c in codes], axis=0)                     # (16,ngage,nyear)
    tv = 0.5 * np.abs(np.diff(stack, axis=2)).sum(axis=0)              # (ngage,nyear-1) pp changed per year
    volatility = tv.mean(axis=1); max_oneyr = tv.max(axis=1)
    shrub_grass_swap = np.minimum(np.abs(np.diff(C[52], axis=1)), np.abs(np.diff(C[71], axis=1))).max(axis=1)

    dev = groups["developed"]                                          # 21-24 combined
    def r(a, nd=2):
        return [None if (v is None or (isinstance(v, float) and not np.isfinite(v))) else round(float(v), nd) for v in a]

    minc = df.groupby("gage_id").valid_coverage_frac.min().reindex(gages).to_numpy()
    medpx = df.groupby("gage_id").n_nlcd_pixels.median().reindex(gages).to_numpy()

    DNOTE = " Endpoint difference = mean(2023-25) - mean(1985-87); NOT a fitted trend. The 2024-25 endpoint is in NLCD's rule-based annual-update path (a methodological seam) and annual NLCD carries model noise, so read change cautiously."
    cont = [
     ("class_diversity", r(shannon, 3), "Normalized Shannon diversity (H/log16, 0-1) of the temporal-mean land-cover composition; higher = more mixed."),
     ("impervious_mean", r(tmean(imp), 2), "Temporal-mean basin fractional impervious surface % (1985-2025). Continuous sub-pixel; NOT a land-cover class - do NOT add to class %s."),
     ("delta_impervious", r(endpt(imp), 2), "Endpoint difference in impervious % (+urbanization)." + DNOTE),
     ("delta_developed", r(endpt(dev), 1), "Endpoint difference in % developed (NLCD 21-24 combined)." + DNOTE),
     ("delta_forest", r(endpt(groups["forest"]), 1), "Endpoint difference in % forest (41-43)." + DNOTE),
     ("delta_agriculture", r(endpt(groups["agriculture"]), 1), "Endpoint difference in % agriculture (81 pasture + 82 crop)." + DNOTE),
     ("annual_volatility", r(volatility, 2), "Mean year-to-year land-cover turnover (pp/yr, half total-variation). High values often flag annual model noise (esp. shrub<->grass), not real change."),
     ("max_oneyr_change", r(max_oneyr, 1), "Largest single-year land-cover turnover (pp). Screens for artifacts / abrupt disturbance; QA-inspect the time series."),
     ("shrub_grass_swap", r(shrub_grass_swap, 1), "Max single-year min(|d shrub|,|d grass|) - a shrub<->grassland swap signature (documented annual-NLCD artifact in the arid West)."),
     ("median_n_pixels", r(medpx, 0), "Median effective (coverage-weighted) 30 m pixel count (basin-size proxy)."),
     ("min_coverage_frac", r(minc, 3), "Minimum-over-years valid_coverage_frac (1 = fully inside the CONUS footprint; <1 = edge basin)."),
    ]
    cat = [
     ("dominant_class", dominant_class, "Largest single NLCD class (temporal mean).", names),
     ("dominant_group", dominant_group, "Largest aggregated land-cover group (lump-robust).", gnames),
     ("watershed_geom_source", meta1.watershed_geom_source.reindex(gages).fillna("NA").tolist(),
      "Basin polygon source: gagesii / wsc_eccc / hydrobasins.", None),
     ("low_confidence", ["yes" if bool(v) else "no" for v in meta1.low_confidence.reindex(gages).fillna(False)],
      "Union QA flag: geometry low-confidence OR <100 eff px OR partial coverage.", ["no", "yes"]),
     ("partial_coverage", ["yes" if bool(v) else "no" for v in meta1.partial_coverage.reindex(gages).fillna(False)],
      "Basin partly outside the CONUS raster footprint (valid_coverage_frac<0.99).", ["no", "yes"]),
    ]

    VARS = []; V = {}
    def add_cont(key, vals, desc, lo_zero=False):
        fin = [v for v in vals if v is not None]
        if not fin: dom = [0, 1]
        elif lo_zero:
            hi = float(np.nanpercentile(fin, 98)); dom = [0, round(hi if hi > 0 else max(fin) or 1, 2)]
        else:
            dom = [round(float(np.nanpercentile(fin, 2)), 2), round(float(np.nanpercentile(fin, 98)), 2)]
        if dom[0] == dom[1]: dom[1] = dom[0] + 1
        VARS.append({"key": key, "label": key, "group": "Summary", "type": "continuous", "domain": dom, "desc": desc}); V[key] = vals
    def add_cat(key, vals, desc, cats):
        cs = cats if cats else sorted(set(x for x in vals if x is not None))
        VARS.append({"key": key, "label": key, "group": "Summary", "type": "categorical", "cats": cs, "desc": desc}); V[key] = vals
    for key, vals, desc in cont: add_cont(key, vals, desc, lo_zero=key in ("impervious_mean", "annual_volatility", "max_oneyr_change", "shrub_grass_swap"))
    for key, vals, desc, cats in cat: add_cat(key, vals, desc, cats)
    # each NLCD class as a temporal-mean % map variable
    for j, c in enumerate(codes):
        VARS.append({"key": f"class:{names[j]}", "label": names[j], "group": "NLCD classes (temporal-mean %)",
                     "type": "continuous", "domain": [0, max(1.0, round(float(np.nanpercentile(comp_mean[:, j], 98)), 1))],
                     "desc": f"NLCD class {c} \"{CLASS_LABEL[c]}\" - temporal-mean % of the watershed (1985-2025)."})
        V[f"class:{names[j]}"] = r(comp_mean[:, j], 1)
    print(f"  {len(VARS)} map variables ({len(cont)+len(cat)} summary + {len(codes)} classes)", flush=True)

    POINTS = {"n": ngage, "gage_id": gages,
              "lat": r(meta1.latitude.reindex(gages).tolist(), 5), "lon": r(meta1.longitude.reindex(gages).tolist(), 5),
              "partial": [1 if bool(v) else 0 for v in meta1.partial_coverage.reindex(gages).fillna(False)],
              "lowconf": [1 if bool(v) else 0 for v in meta1.low_confidence.reindex(gages).fillna(False)],
              "v": V}

    # compact per-gage series: class composition (16 x nyear, 0.1% ints) + impervious (nyear, 0.1% ints)
    print("  encoding composition + impervious series (0.1%) ...", flush=True)
    Ai = np.clip(np.round(stack * 10), 0, 1000).astype(int)            # (16,ngage,nyear)
    Ii = np.clip(np.round(imp * 10), 0, 1000).astype(int)              # (ngage,nyear)
    ser = [";".join(",".join(map(str, Ai[:, i, y].tolist())) for y in range(nyear)) for i in range(ngage)]
    iser = [",".join(map(str, Ii[i].tolist())) for i in range(ngage)]
    TS = {"years": years, "names": names, "codes": codes, "colors": colors,
          "series": ser, "imperv": iser, "seam_from": SEAM_FROM}
    EXCLUDED = {"n": 45, "note": "45 Alaska gages are outside the CONUS Annual NLCD footprint and are excluded "
                "(see nlcd_out_of_footprint CSV). They are NOT shown here and are NOT zero-valued."}

    html = (TEMPLATE.replace("__POINTS__", json.dumps(POINTS)).replace("__VARS__", json.dumps(VARS))
            .replace("__TS__", json.dumps(TS)).replace("__EXCLUDED__", json.dumps(EXCLUDED)))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f: f.write(html)
    print(f"\nwrote {OUT}  ({os.path.getsize(OUT)/1e6:.1f} MB)")
    print("  dominant-class mix:", dict(Counter(dominant_class).most_common(6)))
    print("  dominant-group mix:", dict(Counter(dominant_group).most_common(6)))


TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"/><meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Annual NLCD — Watershed Explorer</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
 html,body{margin:0;font-family:system-ui,Arial,sans-serif;color:#222}
 #app{display:flex;height:100vh}
 #side{width:350px;min-width:350px;padding:14px 16px;box-sizing:border-box;overflow-y:auto;border-right:1px solid #ddd;background:#fafafa}
 #main{flex:1;display:flex;flex-direction:column;min-width:0;position:relative}
 #map{flex:1}
 #ts{height:360px;border-top:1px solid #ccc;background:#fff;display:none;position:relative}
 #ts.show{display:block}
 h1{font-size:16px;margin:0 0 4px}.sub{color:#666;font-size:12px;margin:0 0 12px}
 label{font-size:12px;font-weight:600;color:#333;display:block;margin:10px 0 3px}
 select,input{width:100%;box-sizing:border-box;padding:5px;font-size:13px}
 #legend{margin-top:14px;font-size:12px}
 #legend .bar{height:12px;border-radius:2px;margin:4px 0;background:linear-gradient(to right,#440154,#3b528b,#21918c,#5ec962,#fde725)}
 .legrow{display:flex;align-items:center;gap:6px;margin:2px 0}
 .sw{width:12px;height:12px;border-radius:2px;display:inline-block;border:1px solid #999}
 #count{font-size:11px;color:#888;margin-top:12px}
 #vdesc{font-size:11px;color:#555;margin-top:6px;line-height:1.4;font-style:italic}
 #tech{font-size:11px;color:#555;line-height:1.5;margin-top:16px;border-top:1px solid #ddd;padding-top:10px}
 #tech h3{font-size:12px;margin:0 0 6px}#tech b{color:#333}.warn{color:#9a3412}
 #tsclose{position:absolute;top:6px;right:10px;cursor:pointer;font-size:18px;color:#888;border:none;background:none;z-index:5}
 #tstitle{position:absolute;top:8px;left:14px;font-size:13px;font-weight:600;z-index:5}
 #tssub{position:absolute;top:26px;left:14px;right:40px;font-size:11px;color:#666;z-index:5}
 #tslegend{position:absolute;top:44px;left:14px;right:14px;font-size:10px;z-index:5;display:flex;flex-wrap:wrap;gap:3px 10px;max-height:32px;overflow:hidden}
 #tslegend .lg{display:inline-flex;align-items:center;gap:3px}#tslegend .sw{width:10px;height:10px}
 .hbox{position:fixed;pointer-events:none;background:#fff;border:1px solid #aaa;border-radius:4px;padding:6px 9px;font-size:11px;line-height:1.5;box-shadow:0 2px 8px rgba(0,0,0,.18);z-index:1000;display:none;max-width:300px}
 .hbox b{color:#111}.hbox .k{color:#777}.hbox .hi{background:#fff7cc}
 #maphint{position:absolute;top:10px;left:50%;transform:translateX(-50%);z-index:1001;background:#2b6cb0;color:#fff;font-size:12px;padding:6px 14px;border-radius:16px;box-shadow:0 2px 6px rgba(0,0,0,.25)}
 #defsbox{position:absolute;top:10px;right:10px;z-index:1001;text-align:right}
 #defsbtn{background:#fff;border:1px solid #aaa;border-radius:5px;padding:6px 10px;font-size:12px;font-weight:600;cursor:pointer;box-shadow:0 1px 4px rgba(0,0,0,.18)}
 #defspanel{display:none;margin-top:6px;width:500px;max-width:46vw;max-height:76vh;overflow-y:auto;background:#fff;border:1px solid #aaa;border-radius:6px;box-shadow:0 4px 14px rgba(0,0,0,.22);padding:12px 14px;text-align:left}
 #defspanel.show{display:block}#defspanel h2{font-size:14px;margin:0 0 4px}
 #defspanel .dnote{color:#555;font-size:11px;margin:0 0 10px;line-height:1.45;background:#fff7ed;border:1px solid #fed7aa;border-radius:5px;padding:8px 10px}
 .defrow{display:flex;gap:10px;font-size:12px;padding:3px 0;border-bottom:1px solid #f0f0f0}
 .defrow code{color:#2b6cb0;font-weight:600;min-width:150px;display:inline-block}.defrow span{color:#444}
 .defcat{font-size:12px;font-weight:700;margin:12px 0 4px;color:#333}
 .clrow{display:inline-flex;align-items:center;gap:4px;margin:2px 10px 2px 0;font-size:11px}
</style></head><body>
<div id="app">
 <div id="side">
  <h1>Annual NLCD — Watershed Explorer</h1>
  <p class="sub">~6,119 CONUS watersheds, annual land cover + impervious 1985–2025 (Collection 1 C1V2, 30 m).
     Points colored by a summary feature or any NLCD class; click a gage for its composition over time.</p>
  <label>Find a gage</label>
  <input id="gagePick" list="gageList" placeholder="type a gage id..."/><datalist id="gageList"></datalist>
  <label>Filter variables</label>
  <input id="filter" type="text" placeholder="forest, impervious, delta, dominant, volatility..."/>
  <label>Variable (color)</label><select id="varSel" size="1"></select>
  <div id="vdesc"></div>
  <label style="margin-top:8px"><input type="checkbox" id="hiCov" style="width:auto"> only fully-covered basins (hide partial)</label>
  <div id="legend"></div><div id="count"></div>
  <div id="tech">
   <h3>How this NLCD product is made</h3>
   <p><b>Source.</b> USGS/MRLC <b>Annual NLCD Collection 1 (C1V2)</b>, 30 m, CONUS, 1985–2025. Land cover
      (16 classes) + Fractional Impervious Surface, aggregated to each gage's upstream watershed with
      coverage-weighted <code>exactextract</code> (basins reprojected into the native WGS84-Albers grid;
      the categorical raster is never resampled). Class %s are of the basin's <b>valid</b> pixels (fill=250
      excluded), so the 16 classes sum to ~100.</p>
   <p class="warn"><b>Read cautiously.</b> (1) Annual NLCD is model-derived — raw year-to-year change can be
      <b>model noise</b>, not real change (documented shrub↔grassland swaps in the arid West; developed/barren
      over-water near the Great Lakes). Use the <code>annual_volatility</code> / <code>shrub_grass_swap</code>
      maps to spot it. (2) Deltas are <b>endpoint differences</b> (mean 2023–25 vs 1985–87), not trends, and
      2024–25 is a rule-based <b>annual-update seam</b> (shaded in the time series). (3) <b>Developed %</b>
      (classes 21–24) and <b>impervious %</b> are related but <b>non-additive</b> — don't sum them.
      (4) 45 Alaska gages are outside CONUS and excluded (see glossary), not shown as zeros.</p>
  </div>
 </div>
 <div id="main">
  <div id="map"></div>
  <div id="maphint">👆 Click any gage (or use “Find a gage”) to see its land-cover composition + impervious over time</div>
  <div id="defsbox"><button id="defsbtn">▸ Variables &amp; classes</button><div id="defspanel"></div></div>
  <div id="ts"><div id="tstitle"></div><div id="tssub"></div><div id="tslegend"></div>
   <button id="tsclose" title="close">&times;</button><svg id="tssvg" width="100%" height="100%"></svg></div>
 </div>
</div>
<div class="hbox" id="hover"></div><div class="hbox" id="tshover"></div>
<script>
const POINTS=__POINTS__, VARS=__VARS__, TS=__TS__, EXCLUDED=__EXCLUDED__;
</script><script>
const cvs=L.canvas({padding:0.5});
const map=L.map('map',{preferCanvas:true}).setView([39.5,-98],4);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
 {attribution:'&copy; OpenStreetMap &copy; CARTO',subdomains:'abcd',maxZoom:19}).addTo(map);
const VIR=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
const lerp=(a,b,t)=>a+(b-a)*t;
function viridis(t){t=Math.max(0,Math.min(1,t));const s=t*(VIR.length-1),i=Math.floor(s),f=s-i,a=VIR[i],b=VIR[Math.min(i+1,VIR.length-1)];
 return `rgb(${Math.round(lerp(a[0],b[0],f))},${Math.round(lerp(a[1],b[1],f))},${Math.round(lerp(a[2],b[2],f))})`;}
const CATPAL=['#e6194b','#3cb44b','#ffe119','#4363d8','#f58231','#911eb4','#42d4f4','#f032e6','#bfef45','#fabed4','#469990','#dcbeff','#9a6324','#800000','#aaffc3','#808000','#ffd8b1','#000075','#a9a9a9','#222'];
const NAMECOLOR={}; TS.names.forEach((n,i)=>NAMECOLOR[n]=TS.colors[i]);
function colorFor(meta,val){
 if(val===null||val===undefined||val===''||(typeof val==='number'&&!isFinite(val)))return null;
 if(meta.type==='categorical'){ if(meta.key==='dominant_class'&&NAMECOLOR[val])return NAMECOLOR[val];
   const i=meta.cats.indexOf(val);return CATPAL[(i<0?CATPAL.length-1:i)%CATPAL.length];}
 const lo=meta.domain[0],hi=meta.domain[1];return viridis((val-lo)/((hi-lo)||1));}
const fmtv=v=>(v===null||v===undefined||(typeof v==='number'&&!isFinite(v)))?'NA':v;
const VARMAP={};VARS.forEach(v=>VARMAP[v.key]=v);
const SUMMARY=VARS.filter(v=>v.group==='Summary');
const hoverEl=document.getElementById('hover');
function pointInfoHTML(i){
 let h=`<b>${POINTS.gage_id[i]}</b> <span class="k">(${(+POINTS.lat[i]).toFixed(3)}, ${(+POINTS.lon[i]).toFixed(3)})</span>`;
 if(POINTS.partial[i])h+=` <span style="color:#b45309">◐ partial coverage</span>`;
 if(curVar&&curVar.group!=='Summary')h+=`<br><span class="hi"><span class="k">${curVar.key}:</span> <b>${fmtv(POINTS.v[curVar.key][i])}</b></span>`;
 ['dominant_class','dominant_group','impervious_mean','class_diversity','annual_volatility'].forEach(k=>{
   const sel=curVar&&k===curVar.key;h+=`<br><span class="${sel?'hi':''}"><span class="k">${k}:</span> ${fmtv(POINTS.v[k][i])}</span>`;});
 return h;}
function showHover(i,e){hoverEl.innerHTML=pointInfoHTML(i);hoverEl.style.display='block';
 const x=e.originalEvent.clientX,y=e.originalEvent.clientY;
 hoverEl.style.left=Math.min(x+14,window.innerWidth-hoverEl.offsetWidth-8)+'px';
 hoverEl.style.top=Math.min(Math.max(8,y+14),window.innerHeight-hoverEl.offsetHeight-8)+'px';}
const markers=[];
for(let i=0;i<POINTS.n;i++){
 const m=L.circleMarker([POINTS.lat[i],POINTS.lon[i]],{renderer:cvs,radius:3,weight:0.4,color:'#444',fillOpacity:0.85,fillColor:'#888'});
 m._i=i;m.on('click',()=>selectGage(i));m.on('mouseover',e=>showHover(i,e));m.on('mousemove',e=>showHover(i,e));
 m.on('mouseout',()=>hoverEl.style.display='none');m.addTo(map);markers.push(m);}
let curVar=null,selI=null,hiCovOnly=false;
function visible(i){return !(hiCovOnly&&POINTS.partial[i]);}
function recolor(key){
 const meta=VARMAP[key];if(!meta)return;curVar=meta;const vals=POINTS.v[key];let n=0;
 for(let i=0;i<markers.length;i++){
   if(!visible(i)){markers[i].setStyle({opacity:0,fillOpacity:0});continue;}
   const c=colorFor(meta,vals?vals[i]:null);if(c)n++;
   markers[i].setStyle({opacity:1,color:POINTS.partial[i]?'#b45309':'#444',weight:POINTS.partial[i]?1.2:0.4,
     fillColor:c||'#ccc',fillOpacity:c?0.85:0.12});}
 document.getElementById('count').textContent=`${POINTS.n} gages | ${n} with a value | ${POINTS.partial.reduce((a,b)=>a+b,0)} partial-coverage (orange ring) | ${EXCLUDED.n} Alaska excluded`;
 document.getElementById('vdesc').textContent=meta.desc||'';updateLegend(meta);
 if(selI!==null&&visible(selI))markers[selI].setStyle({weight:2,color:'#d00'});}
function updateLegend(meta){const el=document.getElementById('legend');
 if(meta.type==='categorical'){let h=`<b>${meta.key}</b><div style="margin-top:4px">`;
   meta.cats.forEach((c,k)=>{const sw=(meta.key==='dominant_class'&&NAMECOLOR[c])?NAMECOLOR[c]:CATPAL[k%CATPAL.length];
     h+=`<div class="legrow"><span class="sw" style="background:${sw}"></span>${c}</div>`;});el.innerHTML=h+'</div>';}
 else el.innerHTML=`<b>${meta.key}</b><div class="bar"></div><div style="display:flex;justify-content:space-between"><span>${meta.domain[0]}</span><span>${meta.domain[1]}</span></div>`;}
const NS='http://www.w3.org/2000/svg';let curComp=null,curImp=null;
function parse(i){const yrs=TS.series[i].split(';'),ncls=TS.names.length,comp=Array.from({length:ncls},()=>new Array(yrs.length).fill(0));
 yrs.forEach((yr,y)=>{const v=yr.split(',');for(let c=0;c<ncls;c++)comp[c][y]=(+v[c]||0)/10;});
 const imp=TS.imperv[i].split(',').map(x=>(+x||0)/10);return{comp,imp};}
const tshoverEl=document.getElementById('tshover');
function selectGage(i){if(selI!==null&&selI!==i)recolorOne(selI);selI=i;
 if(visible(i))markers[i].setStyle({weight:2,color:'#d00'});
 document.getElementById('gagePick').value=POINTS.gage_id[i];
 const mh=document.getElementById('maphint');if(mh)mh.style.display='none';drawTS(i);}
function recolorOne(i){if(!visible(i))return;const c=curVar?colorFor(curVar,POINTS.v[curVar.key][i]):null;
 markers[i].setStyle({weight:POINTS.partial[i]?1.2:0.4,color:POINTS.partial[i]?'#b45309':'#444',fillColor:c||'#ccc',fillOpacity:c?0.85:0.12});}
function drawTS(i){const {comp,imp}=parse(i);curComp=comp;curImp=imp;const Y0=TS.years,n=Y0.length,ncls=TS.names.length;
 document.getElementById('ts').classList.add('show');
 const svg=document.getElementById('tssvg');while(svg.firstChild)svg.removeChild(svg.firstChild);
 const W=svg.clientWidth||svg.parentNode.clientWidth,H=360,mL=42,mR=44,mT=88,mB=28,pw=W-mL-mR,ph=H-mT-mB;
 const X=k=>mL+(n<=1?0:pw*k/(n-1)),Y=v=>mT+ph*(1-v/100);
 function add(t,at){const e=document.createElementNS(NS,t);for(const k in at)e.setAttribute(k,at[k]);svg.appendChild(e);return e;}
 // annual-update seam shading
 const seamK=Y0.indexOf(TS.seam_from);
 if(seamK>0){add('rect',{x:X(seamK-0.5<0?0:seamK)-0,y:mT,width:Math.max(2,X(n-1)-X(seamK)),height:ph,fill:'#000',opacity:0.05});
   const st=add('text',{x:X(seamK),y:mT-2,'font-size':9,fill:'#999'});st.textContent='update seam →';}
 add('line',{x1:mL,y1:mT+ph,x2:mL+pw,y2:mT+ph,stroke:'#999'});add('line',{x1:mL,y1:mT,x2:mL,y2:mT+ph,stroke:'#999'});
 for(let g=0;g<=4;g++){const yv=25*g,yy=Y(yv);add('line',{x1:mL,y1:yy,x2:mL+pw,y2:yy,stroke:'#eee'});
   const t=add('text',{x:mL-5,y:yy+3,'text-anchor':'end','font-size':10,fill:'#666'});t.textContent=yv+'%';}
 for(let k=0;k<n;k++)if(k%5===0||k===n-1){const xx=X(k);add('line',{x1:xx,y1:mT+ph,x2:xx,y2:mT+ph+4,stroke:'#999'});
   const t=add('text',{x:xx,y:mT+ph+18,'text-anchor':'middle','font-size':9,fill:'#666'});t.textContent=Y0[k];}
 const present=[];for(let c=0;c<ncls;c++)if(comp[c].some(v=>v>0.5))present.push(c);
 const lower=new Array(n).fill(0);
 present.forEach(c=>{const upper=lower.map((lo,k)=>lo+comp[c][k]);let pts='';
   for(let k=0;k<n;k++)pts+=`${X(k)},${Y(upper[k])} `;for(let k=n-1;k>=0;k--)pts+=`${X(k)},${Y(lower[k])} `;
   add('polygon',{points:pts.trim(),fill:TS.colors[c],stroke:'#fff','stroke-width':0.3,opacity:0.92});
   for(let k=0;k<n;k++)lower[k]=upper[k];});
 // impervious mean overlaid as a line (same 0-100 %)
 let ip='';for(let k=0;k<n;k++)ip+=`${X(k)},${Y(imp[k])} `;
 add('polyline',{points:ip.trim(),fill:'none',stroke:'#111','stroke-width':1.8,'stroke-dasharray':'4 2',opacity:0.9});
 const il=add('text',{x:mL+pw,y:Y(imp[n-1])-4,'text-anchor':'end','font-size':9,fill:'#111'});il.textContent='impervious %';
 const lg=document.getElementById('tslegend');const meanp=c=>comp[c].reduce((a,b)=>a+b,0)/n;
 const ord=present.slice().sort((a,b)=>meanp(b)-meanp(a));
 lg.innerHTML=ord.map(c=>`<span class="lg"><span class="sw" style="background:${TS.colors[c]}"></span>${TS.names[c]} ${meanp(c).toFixed(0)}%</span>`).join('')
   +`<span class="lg"><span class="sw" style="background:#111"></span>impervious (dashed)</span>`;
 document.getElementById('tstitle').textContent=`Gage ${POINTS.gage_id[i]} — land-cover composition + impervious`;
 const q=k=>fmtv(POINTS.v[k]?POINTS.v[k][i]:null);
 document.getElementById('tssub').innerHTML=`stacked % by NLCD class · impervious dashed · 1985–2025 &nbsp;|&nbsp; dominant=${q('dominant_class')} · group=${q('dominant_group')} · vol=${q('annual_volatility')} pp/yr`
   +(POINTS.partial[i]?` &nbsp;|&nbsp; <span style="color:#b45309">◐ partial coverage</span>`:'')+` &nbsp;|&nbsp; <i>2024–25 = update seam</i>`;
 const guide=add('line',{x1:mL,y1:mT,x2:mL,y2:mT+ph,stroke:'#222','stroke-width':1,opacity:0});
 const ov=add('rect',{x:mL,y:mT,width:pw,height:ph,fill:'transparent','pointer-events':'all'});
 ov.addEventListener('mousemove',ev=>{const rect=svg.getBoundingClientRect();let k=Math.round((ev.clientX-rect.left-mL)/(pw||1)*(n-1));
   k=Math.max(0,Math.min(n-1,k));guide.setAttribute('x1',X(k));guide.setAttribute('x2',X(k));guide.setAttribute('opacity',0.6);
   const rows=present.map(c=>[c,curComp[c][k]]).filter(r=>r[1]>0.05).sort((a,b)=>b[1]-a[1]);
   let h=`<b>${Y0[k]}</b>${Y0[k]>=TS.seam_from?' <span style="color:#b45309">(update seam)</span>':''}`;
   rows.forEach(r=>{h+=`<br><span class="sw" style="display:inline-block;background:${TS.colors[r[0]]}"></span> <span class="k">${TS.names[r[0]]}:</span> ${r[1].toFixed(1)}%`;});
   h+=`<br><span class="sw" style="display:inline-block;background:#111"></span> <span class="k">impervious:</span> ${curImp[k].toFixed(1)}%`;
   tshoverEl.innerHTML=h;tshoverEl.style.display='block';
   tshoverEl.style.left=Math.min(ev.clientX+14,window.innerWidth-tshoverEl.offsetWidth-8)+'px';
   tshoverEl.style.top=Math.max(8,ev.clientY-tshoverEl.offsetHeight-10)+'px';});
 ov.addEventListener('mouseleave',()=>{guide.setAttribute('opacity',0);tshoverEl.style.display='none';});}
document.getElementById('tsclose').onclick=()=>{document.getElementById('ts').classList.remove('show');tshoverEl.style.display='none';if(selI!==null){recolorOne(selI);selI=null;}};
window.addEventListener('resize',()=>{if(selI!==null&&document.getElementById('ts').classList.contains('show'))drawTS(selI);});
function buildDefs(){const el=document.getElementById('defspanel');
 let h='<h2>Variables &amp; classes</h2><div class="dnote"><b>Read me first.</b> '
  +'Deltas are <b>endpoint differences</b> mean(2023–25)−mean(1985–87), not fitted trends; 2024–25 is a rule-based '
  +'<b>annual-update seam</b>. Annual NLCD carries model noise — <code>annual_volatility</code>, <code>max_oneyr_change</code> '
  +'and <code>shrub_grass_swap</code> help spot the documented shrub↔grassland artifact. <b>Developed %</b> (classes 21–24) '
  +'and <b>impervious %</b> are related but <b>non-additive</b>. Diversity is normalized Shannon (H/log 16, 0–1). '
  +'<b>'+EXCLUDED.n+' Alaska gages</b> are outside CONUS and excluded from this product (not shown, not zero) — '+EXCLUDED.note+'</div>';
 h+='<div class="defcat">Summary variables</div>';
 SUMMARY.forEach(v=>{const t=v.type==='categorical'?' (categorical)':' (continuous; color 2nd–98th pct)';
   h+=`<div class="defrow"><code>${v.key}</code><span>${v.desc}${t}</span></div>`;});
 h+='<div class="defcat">NLCD land-cover classes (each also a temporal-mean % map variable)</div><div style="padding:2px 0 8px">';
 TS.names.forEach((nm,c)=>{h+=`<span class="clrow"><span class="sw" style="background:${TS.colors[c]}"></span>${nm} (${TS.codes[c]})</span>`;});
 h+='</div>';el.innerHTML=h;}
const varSel=document.getElementById('varSel'),filterEl=document.getElementById('filter');
function populateVars(){const f=filterEl.value.toLowerCase();varSel.innerHTML='';const groups=[],gmap={};let first=null,count=0;
 VARS.forEach(v=>{if(!(v.key.toLowerCase().includes(f)||(v.label||'').toLowerCase().includes(f)||(v.desc||'').toLowerCase().includes(f)))return;
   if(!(v.group in gmap)){gmap[v.group]=[];groups.push(v.group);}gmap[v.group].push(v);});
 groups.forEach(g=>{const og=document.createElement('optgroup');og.label=g;
   gmap[g].forEach(v=>{const o=document.createElement('option');o.value=v.key;o.textContent=(v.group==='Summary'?v.key:'  '+v.label);
     og.appendChild(o);if(first===null)first=v.key;count++;});varSel.appendChild(og);});
 varSel.size=Math.max(2,Math.min(22,count));if(first!==null){varSel.value=first;recolor(first);}}
varSel.onchange=()=>recolor(varSel.value);filterEl.oninput=populateVars;
document.getElementById('hiCov').onchange=function(){hiCovOnly=this.checked;if(curVar)recolor(curVar.key);};
const dl=document.getElementById('gageList');POINTS.gage_id.forEach(g=>{const o=document.createElement('option');o.value=g;dl.appendChild(o);});
const gidx={};POINTS.gage_id.forEach((g,i)=>gidx[g]=i);
document.getElementById('gagePick').onchange=function(){const i=gidx[this.value.trim()];if(i!==undefined){map.setView([POINTS.lat[i],POINTS.lon[i]],Math.max(map.getZoom(),6));selectGage(i);}};
const defsbtn=document.getElementById('defsbtn'),defspanel=document.getElementById('defspanel');
defsbtn.onclick=()=>{const o=defspanel.classList.toggle('show');defsbtn.textContent=(o?'▾':'▸')+' Variables & classes';};
buildDefs();populateVars();
</script></body></html>
"""

if __name__ == "__main__":
    main()

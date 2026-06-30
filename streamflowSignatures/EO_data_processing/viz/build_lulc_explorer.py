"""Build a self-contained HTML QA/QC explorer for the per-watershed annual MODIS LULC product.

Parallels build_lai_explorer.py. Map: 7,964 gage points colored by a selectable variable. Two
families of map variables are offered: (1) IGBP-derived SUMMARY features (dominant class/group,
diversity, first->last change, QA) and (2) every individual IGBP (LC_Type1) class as its own
record-mean % (the five forest types, woody savanna, etc. each directly mappable). Summary vars
that merely restate an IGBP class (pct_water/urban/etc.) are intentionally omitted to avoid
overlap. Hover a point for the summary values; click a gage (or pick from the dropdown) -> a
STACKED-AREA chart of land-cover composition (%) over 2001-2024 for a selectable band; the chart
can still switch across all 8 MODIS classification schemes. Hover a year for the full class
breakdown. An expandable glossary documents variables, the per-band class palettes, and the MODIS
interpretation caveats; a sidebar panel documents the processing.

Design mirrors the LAI explorer (Leaflet + Canvas, viridis/categorical scales). Per-gage per-band
annual composition is encoded as compact LAI-style integer strings (% x10, 0.1% resolution) parsed
lazily on click, so the whole product inlines into one double-clickable file.

Usage: /home/sagemaker-user/.conda/envs/geo/bin/python EO_data_processing/viz/build_lulc_explorer.py
Output: data_out/eo_lai/lulc_explorer.html
"""
import os, json
from collections import Counter
import numpy as np, pandas as pd, geopandas as gpd

ROOT="/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures"
LULC=f"{ROOT}/data_out/eo_lai/lulc_out/watershed_modis_lulc_annual_29jun2026.parquet"
LEG =f"{ROOT}/EO_data_processing/eo_processing/lulc_legends.csv"
GEOM=f"{ROOT}/data_out/eo_lai/geometry_7964.gpkg"
OUT =f"{ROOT}/data_out/eo_lai/lulc_explorer.html"
SUMMARY_GROUP="Summary (IGBP-derived)"

BAND_LABEL={
 "LC_Type1":"IGBP global vegetation (17)","LC_Type2":"UMD (16)",
 "LC_Type3":"MODIS LAI/fPAR biomes (11)","LC_Type4":"BGC / NPP (9)",
 "LC_Type5":"Plant Functional Types (12)","LC_Prop1":"FAO-LCCS land cover (16)",
 "LC_Prop2":"FAO-LCCS land use (11)","LC_Prop3":"FAO-LCCS surface hydrology (10)"}
BAND_ORDER=["LC_Type1","LC_Type2","LC_Type3","LC_Type4","LC_Type5","LC_Prop1","LC_Prop2","LC_Prop3"]

# class-name -> color: official MODIS MCD12Q1 IGBP palette, matched across LC_Type* bands by name
NAME_COLOR={
 "Water":"#1c0dff","EvergreenNeedleleaf":"#05450a","EvergreenBroadleaf":"#086a10",
 "DeciduousNeedleleaf":"#54a708","DeciduousBroadleaf":"#78d203","MixedForest":"#009900",
 "ClosedShrub":"#c6b044","OpenShrub":"#dcd159","Shrub":"#cdb33b","WoodySavanna":"#dade48",
 "Savanna":"#fbff13","Grassland":"#b6ff05","Grass":"#b6ff05","Grass_Cereal":"#d2e878",
 "PermanentWetland":"#27ff87","Cropland":"#c24f44","CerealCrop":"#c98a4b","BroadleafCrop":"#e8a87c",
 "Urban":"#a5a5a5","CropNatMosaic":"#ff6d4c","SnowIce":"#69fff8","Barren":"#f9ffa4",
 "NonVegetated":"#cdc5b4","EvergreenNeedleleafTree":"#05450a","EvergreenBroadleafTree":"#086a10",
 "DeciduousNeedleleafTree":"#54a708","DeciduousBroadleafTree":"#78d203",
 "EvergreenNeedleleafVeg":"#05450a","EvergreenBroadleafVeg":"#086a10",
 "DeciduousNeedleleafVeg":"#54a708","DeciduousBroadleafVeg":"#78d203",
 "AnnualBroadleafVeg":"#e3a3d6","AnnualGrassVeg":"#c6ff8a"}
# deterministic fallback palette for LCCS / unmatched codes
FALLBACK=['#e6194b','#3cb44b','#ffe119','#4363d8','#f58231','#911eb4','#42d4f4','#f032e6',
 '#bfef45','#fabed4','#469990','#dcbeff','#9a6324','#800000','#aaffc3','#808000','#ffd8b1',
 '#000075','#a9a9a9','#666666']

def band_colors(names):
    out=[]; fi=0
    for nm in names:
        if nm in NAME_COLOR: out.append(NAME_COLOR[nm])
        else: out.append(FALLBACK[fi%len(FALLBACK)]); fi+=1
    return out

# IGBP (LC_Type1) class code -> aggregate group, for the derived map summary variables
IGBP_GROUP={1:"forest",2:"forest",3:"forest",4:"forest",5:"forest",6:"shrub",7:"shrub",
 8:"savanna",9:"savanna",10:"grass",11:"wetland",12:"cropland",13:"urban",14:"cropland",
 15:"snow_ice",16:"barren",17:"water"}

def main():
    print("loading LULC panel + legends + geometry ...",flush=True)
    leg=pd.read_csv(LEG)
    df=pd.read_parquet(LULC); df["gage_id"]=df.gage_id.astype(str)
    geo=gpd.read_file(GEOM,ignore_geometry=True); geo["gage_id"]=geo.gage_id.astype(str)
    geo=geo[["gage_id","canon_id","latitude","longitude","watershed_geom_source","low_confidence"]].drop_duplicates("gage_id")

    years=sorted(int(y) for y in df.year.unique()); nyear=len(years)
    gages=geo.gage_id.tolist(); ngage=len(gages)
    print(f"  {len(df):,} rows | {nyear} years {years[0]}..{years[-1]} | {ngage} gages",flush=True)

    # reshape every column to (ngage, nyear) on a complete grid
    full=pd.MultiIndex.from_product([gages,years],names=["gage_id","year"])
    df2=df.set_index(["gage_id","year"]).reindex(full)
    def mat(col): return df2[col].to_numpy().reshape(ngage,nyear)

    rec =lambda M:np.nanmean(M,axis=1)                                   # record-mean over years
    f3=slice(0,3); l3=slice(nyear-3,nyear)                               # robust 3-yr endpoints
    endpt=lambda M:np.nanmean(M[:,l3],axis=1)-np.nanmean(M[:,f3],axis=1) # mean(last3)-mean(first3)

    # ---- IGBP group matrices for the derived summary variables ----
    print("  computing IGBP-derived summary variables ...",flush=True)
    t1=leg[leg.band=="LC_Type1"]
    grp={}
    for g in ["forest","shrub","savanna","grass","wetland","cropland","urban","snow_ice","barren","water"]:
        cols=[r.output_column for r in t1.itertuples() if IGBP_GROUP[r.code]==g]
        grp[g]=np.nansum([np.nan_to_num(mat(c)) for c in cols],axis=0)   # (ngage,nyear) %
    d_forest=endpt(grp["forest"]); d_urban=endpt(grp["urban"]); d_crop=endpt(grp["cropland"])
    npix=mat("n_modis_pixels")[:,0]

    # record-mean IGBP composition -> dominant class, Shannon diversity
    igbp_names=[r.class_name for r in t1.itertuples()]
    igbp_mats=[np.nan_to_num(mat(r.output_column)) for r in t1.itertuples()]  # list (ngage,nyear)
    comp_mean=np.stack([m.mean(axis=1) for m in igbp_mats],axis=1)            # (ngage,17) %
    dominant_class=[igbp_names[j] for j in np.argmax(comp_mean,axis=1)]
    p=np.clip(comp_mean/100.0,1e-12,1); shannon=(-(p*np.log(p)).sum(axis=1))  # nats
    # robust dominant-change: argmax of mean(first3) composition vs mean(last3) composition
    cf=np.stack([m[:,f3].mean(axis=1) for m in igbp_mats],axis=1)
    cl=np.stack([m[:,l3].mean(axis=1) for m in igbp_mats],axis=1)
    dom_change=["yes" if igbp_names[np.argmax(cf[i])]!=igbp_names[np.argmax(cl[i])] else "no" for i in range(ngage)]
    # dominant GROUP (issue 2 fix): argmax over aggregated, lump-robust groups
    GROUPS={"forest":grp["forest"],"savanna":grp["savanna"],"shrubland":grp["shrub"],
            "grassland":grp["grass"],"cropland":grp["cropland"],"wetland":grp["wetland"],
            "water":grp["water"],"urban":grp["urban"],"barren/snow":grp["barren"]+grp["snow_ice"]}
    gnames=list(GROUPS); gmeans=np.stack([rec(GROUPS[g]) for g in gnames],axis=1)
    dominant_group=[gnames[j] for j in np.argmax(gmeans,axis=1)]

    def r(arr,nd=2):
        return [None if (v is None or (isinstance(v,float) and not np.isfinite(v))) else round(float(v),nd) for v in arr]

    DELTA_NOTE=" Endpoints are mean(2001-03) vs mean(2022-24); the recent endpoint is in MODIS's 2021+ reduced-confidence era, so treat changes cautiously."
    # Summary vars are deliberately limited to NON-class-level quantities (diversity, temporal change,
    # categorical dominance, QA). Static class-% roll-ups (pct_forest/cropland/water/urban/...) are
    # omitted because each overlaps the now-directly-mappable individual IGBP classes.
    cont=[  # (key, values, desc)
     ("class_diversity",r(shannon,2),"Shannon diversity (nats) of the RECORD-MEAN IGBP composition; higher = more mixed land cover"),
     ("delta_forest",r(d_forest,1),"Change in % forest (IGBP 1-5 combined; +gain/-loss). Temporal signal, not a static class level."+DELTA_NOTE),
     ("delta_urban",r(d_urban,1),"Change in % urban (+urbanization). Temporal signal, not a static class level."+DELTA_NOTE),
     ("delta_cropland",r(d_crop,1),"Change in % cropland (IGBP cropland + cropland/natural mosaic). Temporal signal, not a static class level."+DELTA_NOTE),
     ("n_modis_pixels",r(npix,0),"Effective coverage-weighted count of 500 m MODIS pixels in the basin (size proxy)"),
    ]
    cat=[  # (key, values, desc, cats-or-None)
     ("dominant_class",dominant_class,"Single largest IGBP class (record-mean). CAVEAT: IGBP splits trees across 5 forest classes but keeps Woody Savanna & Grassland as single buckets, so a tree-dominated basin can read 'Woody Savanna'. Use dominant_group for a lump-robust view",igbp_names),
     ("dominant_group",dominant_group,"Largest AGGREGATED land-cover group (forest/savanna/shrubland/grassland/cropland/wetland/water/urban/barren-snow) - robust to IGBP's class splitting",gnames),
     ("dominant_change",dom_change,"Did the dominant IGBP class differ between mean(2001-03) and mean(2022-24)? (turnover; recent years are 2021+ reduced-confidence)",["no","yes"]),
     ("watershed_geom_source",geo.watershed_geom_source.fillna("NA").tolist(),"Basin polygon source: gagesii / wsc_eccc / hydrobasins",None),
     ("low_confidence",[("yes" if int(v)==1 else "no") for v in geo.low_confidence.fillna(0)],"Geometry/coverage low-confidence flag (HydroBasins fallback or area outlier)",None),
    ]
    VARS=[]; V={}
    def add_cont(key,vals,desc,group,label,lo_zero=False):
        finite=[v for v in vals if v is not None]
        if not finite: dom=[0,1]
        elif lo_zero:
            hi=float(np.nanpercentile(finite,98));  hi=hi if hi>0 else (max(finite) if max(finite)>0 else 1)
            dom=[0,round(hi,2)]
        else:
            dom=[round(float(np.nanpercentile(finite,2)),2),round(float(np.nanpercentile(finite,98)),2)]
        if dom[0]==dom[1]: dom[1]=dom[0]+1
        VARS.append({"key":key,"label":label,"group":group,"type":"continuous","domain":dom,"desc":desc}); V[key]=vals
    def add_cat(key,vals,desc,group,label,cats):
        cs=cats if cats else sorted(set(x for x in vals if x is not None))
        VARS.append({"key":key,"label":label,"group":group,"type":"categorical","cats":cs,"desc":desc}); V[key]=vals

    for key,vals,desc in cont: add_cont(key,vals,desc,SUMMARY_GROUP,key,lo_zero=False)
    for key,vals,desc,cats in cat: add_cat(key,vals,desc,SUMMARY_GROUP,key,cats)

    # ---- every individual IGBP (LC_Type1) class as its own record-mean % map variable ----
    # (only IGBP; the other 7 bands remain available in the time-series band switcher, not as map vars)
    print("  adding individual IGBP class variables (record-mean %) ...",flush=True)
    band="LC_Type1"; sub=leg[leg.band==band]; lbl=BAND_LABEL[band]
    for row in sub.itertuples():
        arr=r(rec(np.nan_to_num(mat(row.output_column))),1)
        key=f"{band}:{row.class_name}"
        add_cont(key,arr,f"{lbl} - record-mean % of watershed in class '{row.class_name}' (col {row.output_column})",
                 lbl,row.class_name,lo_zero=True)
    print(f"  total map variables: {len(VARS)} ({len(cont)+len(cat)} summary + {len(VARS)-len(cont)-len(cat)} IGBP classes)",flush=True)

    POINTS={"n":ngage,"lat":r(geo.latitude.tolist(),5),"lon":r(geo.longitude.tolist(),5),
            "gage_id":gages,"v":V}

    # ---- per-band composition strings (% x10 int; ';'-years, ','-classes) + band meta ----
    print("  encoding band composition time series (0.1% resolution) ...",flush=True)
    BANDS=[]; SER={}
    for band in BAND_ORDER:
        sub=leg[leg.band==band]
        names=[row.class_name for row in sub.itertuples()]
        codes=[int(row.code) for row in sub.itertuples()]
        cols=[row.output_column for row in sub.itertuples()]
        colors=band_colors(names)
        A=np.stack([np.nan_to_num(mat(c)) for c in cols])              # (ncls,ngage,nyear)
        Ai=np.clip(np.round(A*10),0,1000).astype(int)                  # x10 -> 0.1% resolution
        ser=[";".join(",".join(map(str,Ai[:,i,y].tolist())) for y in range(nyear)) for i in range(ngage)]
        BANDS.append({"key":band,"label":BAND_LABEL[band],"names":names,"codes":codes,"colors":colors})
        SER[band]=ser
    TSB={"years":years,"series":SER}

    html=TEMPLATE.replace("__POINTS__",json.dumps(POINTS)) \
                 .replace("__VARS__",json.dumps(VARS)) \
                 .replace("__BANDS__",json.dumps(BANDS)) \
                 .replace("__TSB__",json.dumps(TSB)) \
                 .replace("__SUMMARY_GROUP__",json.dumps(SUMMARY_GROUP))
    os.makedirs(os.path.dirname(OUT),exist_ok=True)
    with open(OUT,"w") as f: f.write(html)
    mb=os.path.getsize(OUT)/1e6
    print(f"\nwrote {OUT}  ({mb:.1f} MB)")
    print(f"  {ngage} gages | {len(VARS)} map variables | {len(BANDS)} bands x {nyear} years")
    print("  dominant-class mix:",dict(Counter(dominant_class).most_common(6)))
    print("  dominant-group mix:",dict(Counter(dominant_group).most_common(6)))

TEMPLATE=r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>MODIS LULC — Watershed Explorer</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  html, body { margin:0; font-family: system-ui, Arial, sans-serif; color:#222; }
  #app { display:flex; height:100vh; }
  #side { width:340px; min-width:340px; padding:14px 16px; box-sizing:border-box;
          overflow-y:auto; border-right:1px solid #ddd; background:#fafafa; }
  #main { flex:1; display:flex; flex-direction:column; min-width:0; position:relative; }
  #map { flex:1; }
  #ts { height:344px; border-top:1px solid #ccc; background:#fff; display:none; position:relative; }
  #ts.show { display:block; }
  h1 { font-size:16px; margin:0 0 4px; }
  .sub { color:#666; font-size:12px; margin:0 0 12px; }
  label { font-size:12px; font-weight:600; color:#333; display:block; margin:10px 0 3px; }
  select, input { width:100%; box-sizing:border-box; padding:5px; font-size:13px; }
  #legend { margin-top:14px; font-size:12px; }
  #legend .bar { height:12px; border-radius:2px; margin:4px 0;
                 background:linear-gradient(to right,#440154,#3b528b,#21918c,#5ec962,#fde725); }
  .legrow { display:flex; align-items:center; gap:6px; margin:2px 0; }
  .sw { width:12px; height:12px; border-radius:2px; display:inline-block; border:1px solid #999; }
  #count { font-size:11px; color:#888; margin-top:12px; }
  #vdesc { font-size:11px; color:#555; margin-top:6px; line-height:1.4; font-style:italic; }
  #tech { font-size:11px; color:#555; line-height:1.5; margin-top:16px; border-top:1px solid #ddd; padding-top:10px; }
  #tech h3 { font-size:12px; margin:0 0 6px; color:#333; }
  #tech b { color:#333; }
  #tech .warn { color:#9a3412; }
  #tsclose { position:absolute; top:6px; right:10px; cursor:pointer; font-size:18px; color:#888; border:none; background:none; z-index:5; }
  #tstitle { position:absolute; top:8px; left:14px; font-size:13px; font-weight:600; z-index:5; }
  #tssub { position:absolute; top:26px; left:14px; right:40px; font-size:11px; color:#666; z-index:5; }
  #tslegend { position:absolute; top:44px; left:14px; right:14px; font-size:10px; z-index:5;
              display:flex; flex-wrap:wrap; gap:3px 10px; max-height:34px; overflow:hidden; }
  #tslegend .lg { display:inline-flex; align-items:center; gap:3px; }
  #tslegend .sw { width:10px; height:10px; }
  .hbox { position:fixed; pointer-events:none; background:#fff; border:1px solid #aaa; border-radius:4px;
          padding:6px 9px; font-size:11px; line-height:1.5; box-shadow:0 2px 8px rgba(0,0,0,.18);
          z-index:1000; display:none; max-width:300px; }
  .hbox b { color:#111; } .hbox .k { color:#777; } .hbox .hi { background:#fff7cc; }
  #maphint { position:absolute; top:10px; left:50%; transform:translateX(-50%); z-index:1001;
             background:#2b6cb0; color:#fff; font-size:12px; padding:6px 14px; border-radius:16px;
             box-shadow:0 2px 6px rgba(0,0,0,.25); }
  #defsbox { position:absolute; top:10px; right:10px; z-index:1001; text-align:right; }
  #defsbtn { background:#fff; border:1px solid #aaa; border-radius:5px; padding:6px 10px; font-size:12px;
             font-weight:600; cursor:pointer; box-shadow:0 1px 4px rgba(0,0,0,.18); }
  #defspanel { display:none; margin-top:6px; width:480px; max-width:46vw; max-height:76vh; overflow-y:auto;
               background:#fff; border:1px solid #aaa; border-radius:6px; box-shadow:0 4px 14px rgba(0,0,0,.22);
               padding:12px 14px; text-align:left; }
  #defspanel.show { display:block; }
  #defspanel h2 { font-size:14px; margin:0 0 4px; }
  #defspanel .dnote { color:#555; font-size:11px; margin:0 0 10px; line-height:1.45;
                      background:#fff7ed; border:1px solid #fed7aa; border-radius:5px; padding:8px 10px; }
  .defrow { display:flex; gap:10px; font-size:12px; padding:3px 0; border-bottom:1px solid #f0f0f0; }
  .defrow code { color:#2b6cb0; font-weight:600; min-width:150px; display:inline-block; }
  .defrow span { color:#444; }
  .defcat { font-size:12px; font-weight:700; margin:12px 0 4px; color:#333; }
  .clrow { display:inline-flex; align-items:center; gap:4px; margin:2px 10px 2px 0; font-size:11px; }
</style>
</head>
<body>
<div id="app">
  <div id="side">
    <h1>MODIS LULC — Watershed Explorer</h1>
    <p class="sub">7,964 watersheds, annual land cover 2001–2024. Points colored by an IGBP summary
       feature <i>or</i> any individual IGBP class; hover a point for its summary values, click (or
       use the picker) for its full composition over time (switchable across all 8 MODIS schemes).</p>

    <label>Find a gage</label>
    <input id="gagePick" list="gageList" placeholder="type a gage id..." />
    <datalist id="gageList"></datalist>

    <label>Filter variables</label>
    <input id="filter" type="text" placeholder="type to filter, e.g. forest, savanna, urban, delta, dominant..." />
    <label>Variable (color)</label>
    <select id="varSel" size="1"></select>
    <div id="vdesc"></div>

    <label>Time-series classification band</label>
    <select id="bandSel"></select>

    <div id="legend"></div>
    <div id="count"></div>

    <div id="tech">
      <h3>How this LULC product is made</h3>
      <p><b>Source.</b> MODIS <b>MCD12Q1</b> v061 (500 m, annual), 2001–2024, aggregated to each
         gage's upstream watershed.</p>
      <p><b>Eight classification schemes (bands).</b> Every pixel is labeled under 8 schemes:
         <b>LC_Type1</b> IGBP (17), <b>LC_Type2</b> UMD (16), <b>LC_Type3</b> LAI/fPAR biomes (11),
         <b>LC_Type4</b> BGC (9), <b>LC_Type5</b> Plant Functional Types (12), and <b>LC_Prop1/2/3</b>
         FAO-LCCS land cover / land use / surface hydrology. Map variables are the IGBP summary
         features plus every individual IGBP (LC_Type1) class; the time-series chart can switch
         between all 8 schemes.</p>
      <p><b>Per-class %.</b> For each basin and year, the coverage-weighted fraction of the basin's
         500 m pixels in each class is computed with <code>exactextract</code> (partial pixels weighted
         by overlap). Fill (255) is excluded from the denominator, so each band's class percentages
         sum to ~100 over valid pixels. Polygons are reprojected into the native MODIS sinusoidal grid;
         the categorical raster is never resampled.</p>
      <p class="warn"><b>Interpretation caveats.</b> (1) MODIS IGBP labels much of the SE US (and other
         open/managed-canopy forest) as <b>Woody Savanna</b>, not Evergreen Needleleaf — flip the band
         to PFT to see the same area as broadleaf tree. (2) NASA flags <b>reduced land-cover confidence
         for 2021+</b>, so recent change should be read cautiously.</p>
    </div>
  </div>
  <div id="main">
    <div id="map"></div>
    <div id="maphint">👆 Click any gage (or use “Find a gage” at left) to see its land-cover composition over time</div>
    <div id="defsbox"><button id="defsbtn">▸ Variables &amp; classes</button><div id="defspanel"></div></div>
    <div id="ts">
      <div id="tstitle"></div><div id="tssub"></div><div id="tslegend"></div>
      <button id="tsclose" title="close">&times;</button>
      <svg id="tssvg" width="100%" height="100%"></svg>
    </div>
  </div>
</div>
<div class="hbox" id="hover"></div>
<div class="hbox" id="tshover"></div>

<script>
const POINTS = __POINTS__;
const VARS   = __VARS__;
const BANDS  = __BANDS__;
const TSB    = __TSB__;
const SUMMARY_GROUP = __SUMMARY_GROUP__;
</script>
<script>
// ---------- basemap ----------
const cvs = L.canvas({ padding: 0.5 });
const map = L.map('map', { preferCanvas: true }).setView([48, -100], 4);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
  { attribution:'&copy; OpenStreetMap &copy; CARTO', subdomains:'abcd', maxZoom:19 }).addTo(map);

// ---------- color scales (map points) ----------
const VIR=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
const lerp=(a,b,t)=>a+(b-a)*t;
function viridis(t){ t=Math.max(0,Math.min(1,t)); const s=t*(VIR.length-1);
  const i=Math.floor(s),f=s-i,a=VIR[i],b=VIR[Math.min(i+1,VIR.length-1)];
  return `rgb(${Math.round(lerp(a[0],b[0],f))},${Math.round(lerp(a[1],b[1],f))},${Math.round(lerp(a[2],b[2],f))})`; }
const CATPAL=['#e6194b','#3cb44b','#ffe119','#4363d8','#f58231','#911eb4','#42d4f4','#f032e6',
  '#bfef45','#fabed4','#469990','#dcbeff','#9a6324','#800000','#aaffc3','#808000','#ffd8b1',
  '#000075','#a9a9a9','#222222'];
function colorFor(meta,val){
  if (val===null||val===undefined||val===''||(typeof val==='number'&&!isFinite(val))) return null;
  if (meta.type==='categorical'){ const idx=meta.cats.indexOf(val); return CATPAL[(idx<0?CATPAL.length-1:idx)%CATPAL.length]; }
  const lo=meta.domain[0],hi=meta.domain[1]; return viridis((val-lo)/((hi-lo)||1));
}
const fmtv=v=>(v===null||v===undefined||(typeof v==='number'&&!isFinite(v)))?'NA':v;
const BANDMAP={}; BANDS.forEach(b=>BANDMAP[b.key]=b);
const VARMAP={}; VARS.forEach(v=>VARMAP[v.key]=v);
const SUMMARY_VARS=VARS.filter(v=>v.group===SUMMARY_GROUP);

// ---------- markers + hover (summary vars + the active variable) ----------
const hoverEl=document.getElementById('hover');
function pointInfoHTML(i){
  let h=`<b>${POINTS.gage_id[i]}</b> <span class="k">(${POINTS.lat[i].toFixed(3)}, ${POINTS.lon[i].toFixed(3)})</span>`;
  if (curVar && curVar.group!==SUMMARY_GROUP)
    h+=`<br><span class="hi"><span class="k">${curVar.key}:</span> <b>${fmtv(POINTS.v[curVar.key][i])}</b></span>`;
  SUMMARY_VARS.forEach(v=>{ const sel=curVar&&v.key===curVar.key;
    h+=`<br><span class="${sel?'hi':''}"><span class="k">${v.key}:</span> ${fmtv(POINTS.v[v.key][i])}</span>`; });
  return h;
}
function showHover(i,e){ hoverEl.innerHTML=pointInfoHTML(i); hoverEl.style.display='block';
  const x=e.originalEvent.clientX, y=e.originalEvent.clientY;
  hoverEl.style.left=Math.min(x+14, window.innerWidth-hoverEl.offsetWidth-8)+'px';
  hoverEl.style.top=Math.min(Math.max(8,y+14), window.innerHeight-hoverEl.offsetHeight-8)+'px'; }
const markers=[];
for (let i=0;i<POINTS.n;i++){
  const m=L.circleMarker([POINTS.lat[i],POINTS.lon[i]],
    { renderer:cvs, radius:3, weight:0.4, color:'#444', fillOpacity:0.85, fillColor:'#888' });
  m._i=i; m.on('click',()=>selectGage(i));
  m.on('mouseover',e=>showHover(i,e)); m.on('mousemove',e=>showHover(i,e));
  m.on('mouseout',()=>{ hoverEl.style.display='none'; });
  m.addTo(map); markers.push(m);
}
let curVar=null, selI=null;
function recolor(key){
  const meta=VARMAP[key]; if(!meta) return; curVar=meta;
  const vals=POINTS.v[key]; let n=0;
  for (let i=0;i<markers.length;i++){ const c=colorFor(meta,vals?vals[i]:null);
    if(c) n++; markers[i].setStyle({fillColor:c||'#ccc',fillOpacity:c?0.85:0.12}); }
  document.getElementById('count').textContent=`${POINTS.n} gages | ${n} with a value`;
  document.getElementById('vdesc').textContent=meta.desc||'';
  updateLegend(meta);
  if (selI!==null) markers[selI].setStyle({weight:2,color:'#d00'});
}
function updateLegend(meta){
  const el=document.getElementById('legend');
  if (meta.type==='categorical'){ let h=`<b>${meta.key}</b><div style="margin-top:4px">`;
    meta.cats.forEach((c,k)=>{ h+=`<div class="legrow"><span class="sw" style="background:${CATPAL[k%CATPAL.length]}"></span>${c}</div>`; });
    el.innerHTML=h+'</div>';
  } else { el.innerHTML=`<b>${meta.key}</b><div class="bar"></div>`
    +`<div style="display:flex;justify-content:space-between"><span>${meta.domain[0]}</span><span>${meta.domain[1]}</span></div>`; }
}

// ---------- time series (stacked-area composition) ----------
const NS='http://www.w3.org/2000/svg';
let curBand="LC_Type1", curComp=null, tsGeom=null;
function parseComp(i,band){                       // -> comp[classIdx][yearIdx] in % (x10 ints)
  const b=BANDMAP[band]; const ncls=b.names.length;
  const yrs=TSB.series[band][i].split(';');
  const comp=Array.from({length:ncls},()=>new Array(yrs.length).fill(0));
  yrs.forEach((yr,y)=>{ const v=yr.split(','); for(let c=0;c<ncls;c++) comp[c][y]=(+v[c]||0)/10; });
  return comp;
}
const tshoverEl=document.getElementById('tshover');
function selectGage(i){
  if (selI!==null && selI!==i) recolorOne(selI);
  selI=i; markers[i].setStyle({weight:2,color:'#d00'});
  document.getElementById('gagePick').value=POINTS.gage_id[i];
  const mh=document.getElementById('maphint'); if(mh) mh.style.display='none';
  drawTS(i);
}
function recolorOne(i){ const vals=curVar?POINTS.v[curVar.key]:null; const c=curVar?colorFor(curVar,vals?vals[i]:null):null;
  markers[i].setStyle({weight:0.4,color:'#444',fillColor:c||'#ccc',fillOpacity:c?0.85:0.12}); }

function drawTS(i){
  const band=BANDMAP[curBand]; const comp=parseComp(i,curBand); curComp=comp;
  const Y0=TSB.years, n=Y0.length, ncls=band.names.length;
  document.getElementById('ts').classList.add('show');
  const svg=document.getElementById('tssvg'); while(svg.firstChild) svg.removeChild(svg.firstChild);
  const W=svg.clientWidth||svg.parentNode.clientWidth, H=344;
  const mL=42,mR=12,mT=86,mB=28, pw=W-mL-mR, ph=H-mT-mB;
  const X=k=> mL + (n<=1?0:pw*k/(n-1));
  const Y=v=> mT + ph*(1-v/100);
  tsGeom={mL,mT,pw,ph,n,X};
  function add(tag,at){ const e=document.createElementNS(NS,tag); for(const k in at) e.setAttribute(k,at[k]); svg.appendChild(e); return e; }
  add('line',{x1:mL,y1:mT+ph,x2:mL+pw,y2:mT+ph,stroke:'#999'});
  add('line',{x1:mL,y1:mT,x2:mL,y2:mT+ph,stroke:'#999'});
  for (let g=0; g<=4; g++){ const yv=25*g, yy=Y(yv);
    add('line',{x1:mL,y1:yy,x2:mL+pw,y2:yy,stroke:'#eee'});
    const t=add('text',{x:mL-5,y:yy+3,'text-anchor':'end','font-size':10,fill:'#666'}); t.textContent=yv+'%'; }
  for (let k=0;k<n;k++){ if (k%4===0||k===n-1){ const xx=X(k);
    add('line',{x1:xx,y1:mT+ph,x2:xx,y2:mT+ph+4,stroke:'#999'});
    const t=add('text',{x:xx,y:mT+ph+18,'text-anchor':'middle','font-size':9,fill:'#666'}); t.textContent=Y0[k]; } }
  const present=[]; for(let c=0;c<ncls;c++){ if (comp[c].some(v=>v>0.5)) present.push(c); }
  const lower=new Array(n).fill(0);
  present.forEach(c=>{
    const upper=lower.map((lo,k)=>lo+comp[c][k]);
    let pts=''; for(let k=0;k<n;k++) pts+=`${X(k)},${Y(upper[k])} `;
    for(let k=n-1;k>=0;k--) pts+=`${X(k)},${Y(lower[k])} `;
    add('polygon',{points:pts.trim(),fill:band.colors[c],stroke:'#fff','stroke-width':0.3,opacity:0.92});
    for(let k=0;k<n;k++) lower[k]=upper[k];
  });
  const lg=document.getElementById('tslegend');
  const meanp=c=>comp[c].reduce((a,b)=>a+b,0)/n;
  const ord=present.slice().sort((a,b)=>meanp(b)-meanp(a));
  lg.innerHTML=ord.map(c=>`<span class="lg"><span class="sw" style="background:${band.colors[c]}"></span>${band.names[c]} ${meanp(c).toFixed(0)}%</span>`).join('');
  document.getElementById('tstitle').textContent=`Gage ${POINTS.gage_id[i]} — land-cover composition`;
  const q=k=>fmtv(POINTS.v[k]?POINTS.v[k][i]:null);
  document.getElementById('tssub').innerHTML=
    `${band.label} · stacked % by class · 2001–2024 &nbsp;|&nbsp; dominant=${q('dominant_class')} · group=${q('dominant_group')} · n_pixels=${q('n_modis_pixels')} `
    +`&nbsp;|&nbsp; <i>2021+ reduced confidence</i>`;
  const guide=add('line',{x1:mL,y1:mT,x2:mL,y2:mT+ph,stroke:'#222','stroke-width':1,opacity:0,'pointer-events':'none'});
  const ov=add('rect',{x:mL,y:mT,width:pw,height:ph,fill:'transparent','pointer-events':'all'});
  ov.addEventListener('mousemove',ev=>{
    const rect=svg.getBoundingClientRect(); const px=ev.clientX-rect.left;
    let k=Math.round((px-mL)/(pw||1)*(n-1)); k=Math.max(0,Math.min(n-1,k));
    guide.setAttribute('x1',X(k)); guide.setAttribute('x2',X(k)); guide.setAttribute('opacity',0.6);
    const rows=present.map(c=>[c,curComp[c][k]]).filter(rr=>rr[1]>0.05).sort((a,b)=>b[1]-a[1]);
    let h=`<b>${Y0[k]}</b>`;
    rows.forEach(rw=>{ h+=`<br><span class="sw" style="display:inline-block;background:${band.colors[rw[0]]}"></span> `
      +`<span class="k">${band.names[rw[0]]}:</span> ${rw[1].toFixed(1)}%`; });
    tshoverEl.innerHTML=h; tshoverEl.style.display='block';
    tshoverEl.style.left=Math.min(ev.clientX+14, window.innerWidth-tshoverEl.offsetWidth-8)+'px';
    tshoverEl.style.top=Math.max(8, ev.clientY-tshoverEl.offsetHeight-10)+'px';
  });
  ov.addEventListener('mouseleave',()=>{ guide.setAttribute('opacity',0); tshoverEl.style.display='none'; });
}
document.getElementById('tsclose').onclick=()=>{ document.getElementById('ts').classList.remove('show');
  tshoverEl.style.display='none'; if(selI!==null){ recolorOne(selI); selI=null; } };
window.addEventListener('resize',()=>{ if(selI!==null && document.getElementById('ts').classList.contains('show')) drawTS(selI); });

// ---------- glossary (expandable, top-right) ----------
function buildDefs(){
  const el=document.getElementById('defspanel');
  let h='<h2>Variables &amp; classes</h2>'
   +'<div class="dnote"><b>Read me first.</b> '
   +'Map variables are the IGBP summary features <i>plus</i> every individual IGBP (LC_Type1) class. '
   +'Summary vars that merely restate an IGBP class (e.g. % water / urban / wetland / cropland) are '
   +'omitted — map those directly from the IGBP class list below. '
   +'<b>dominant_class</b> is the single largest of 17 IGBP classes; because trees are split across 5 forest '
   +'classes while Woody Savanna &amp; Grassland are single buckets, a tree-covered basin can read "Woody Savanna" — '
   +'use <b>dominant_group</b> for a lump-robust view, or map the individual forest + Woody Savanna classes. '
   +'MODIS classes much of the SE US pine region as <b>Woody Savanna</b> (IGBP/UMD) or broadleaf tree (PFT), '
   +'not Evergreen Needleleaf — a known product behavior, visible by switching the time-series band. '
   +'<b>delta_*</b> compare mean(2001–03) vs mean(2022–24); the recent endpoint is in MODIS\'s <b>2021+ '
   +'reduced-confidence</b> era. Individual IGBP class % map values are full-precision record means; the '
   +'time-series chart is encoded at 0.1% resolution.</div>';
  h+='<div class="defcat">Summary variables (IGBP-derived; shown in point hover)</div>';
  SUMMARY_VARS.forEach(v=>{ const t=v.type==='categorical'?' (categorical)':' (continuous; color clipped to 2nd–98th pct)';
    h+=`<div class="defrow"><code>${v.key}</code><span>${v.desc}${t}</span></div>`; });
  const ig=BANDMAP['LC_Type1'];
  h+='<div class="defcat">IGBP classes (LC_Type1) — each is also a map variable (record-mean %)</div>';
  h+='<div style="padding:2px 0 8px">';
  ig.names.forEach((nm,c)=>{ h+=`<span class="clrow"><span class="sw" style="background:${ig.colors[c]}"></span>${nm}</span>`; });
  h+='</div>';
  h+='<div class="defcat">Other classification bands (time series only — switch via the band selector)</div>';
  BANDS.filter(b=>b.key!=='LC_Type1').forEach(b=>{ h+=`<div class="defrow"><code>${b.key}</code><span>${b.label}</span></div><div style="padding:2px 0 8px">`;
    b.names.forEach((nm,c)=>{ h+=`<span class="clrow"><span class="sw" style="background:${b.colors[c]}"></span>${nm}</span>`; });
    h+='</div>'; });
  el.innerHTML=h;
}

// ---------- controls (optgroup-grouped, filterable) ----------
const varSel=document.getElementById('varSel'), filterEl=document.getElementById('filter');
function populateVars(){ const f=filterEl.value.toLowerCase();
  varSel.innerHTML=''; const groups=[], gmap={}; let first=null, count=0;
  VARS.forEach(v=>{ if(!(v.key.toLowerCase().includes(f)||(v.label||'').toLowerCase().includes(f)
      ||(v.desc||'').toLowerCase().includes(f)||(v.group||'').toLowerCase().includes(f))) return;
    if(!(v.group in gmap)){ gmap[v.group]=[]; groups.push(v.group); } gmap[v.group].push(v); });
  groups.forEach(g=>{ const og=document.createElement('optgroup'); og.label=g;
    gmap[g].forEach(v=>{ const o=document.createElement('option'); o.value=v.key;
      o.textContent=(v.group===SUMMARY_GROUP?v.key:'  '+v.label)+(v.type==='categorical'?'  (class)':'');
      og.appendChild(o); if(first===null) first=v.key; count++; });
    varSel.appendChild(og); });
  varSel.size=Math.max(2,Math.min(20,count));
  if(first!==null){ varSel.value=first; recolor(first); }
}
varSel.onchange=()=>recolor(varSel.value); filterEl.oninput=populateVars;
const bandSel=document.getElementById('bandSel');
BANDS.forEach(b=>{ const o=document.createElement('option'); o.value=b.key; o.textContent=b.label; bandSel.appendChild(o); });
bandSel.value=curBand;
bandSel.onchange=()=>{ curBand=bandSel.value; if(selI!==null) drawTS(selI); };
const dl=document.getElementById('gageList');
POINTS.gage_id.forEach(g=>{ const o=document.createElement('option'); o.value=g; dl.appendChild(o); });
const gidx={}; POINTS.gage_id.forEach((g,i)=>gidx[g]=i);
document.getElementById('gagePick').onchange=function(){ const i=gidx[this.value.trim()];
  if(i!==undefined){ map.setView([POINTS.lat[i],POINTS.lon[i]],Math.max(map.getZoom(),6)); selectGage(i); } };
const defsbtn=document.getElementById('defsbtn'), defspanel=document.getElementById('defspanel');
defsbtn.onclick=()=>{ const open=defspanel.classList.toggle('show');
  defsbtn.textContent=(open?'▾':'▸')+' Variables & classes'; };
buildDefs(); populateVars();
</script>
</body>
</html>
"""

if __name__=="__main__": main()

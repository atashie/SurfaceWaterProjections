"""Build a self-contained HTML QA/QC explorer for the per-watershed monthly MODIS LAI product.

Map: 7,964 gage points colored by a selectable QA/summary variable (mean LAI, seasonal
amplitude, % NA months, trend, peak/trough month, data-quality fractions, ...). Hover a point
for all its variable values; click a gage (or pick from the searchable dropdown) -> a time-series
panel draws the full 270-month record (2002-07..2024-12): mean line + min/max + quantile bands,
with NA gaps broken and a per-month hover readout. A glossary footer defines every variable and a
sidebar panel documents the processing.

Design mirrors the project's streamflow_explorer.html (Leaflet + Canvas, viridis/categorical
scales). The per-gage monthly series are encoded as compact LAI*10 integer strings and parsed
lazily on click, so the whole product inlines into one double-clickable file.

Usage: /home/sagemaker-user/.conda/envs/geo/bin/python EO_data_processing/viz/build_lai_explorer.py
Output: data_out/eo_lai/lai_explorer.html
"""
import os, json
import numpy as np, pandas as pd, geopandas as gpd

ROOT="/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures"
LAI=f"{ROOT}/data_out/eo_lai/out_qa/watershed_modis_lai_monthly.parquet"
GEOM=f"{ROOT}/data_out/eo_lai/geometry_7964.gpkg"
OUT=f"{ROOT}/data_out/eo_lai/lai_explorer.html"
MONTHS_ABBR=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
# time-series stats to inline (order matters; the JS reads this order). 1-dp (LAI*10 int) encoding.
TS_STATS=["lai_min","lai_q05","lai_q25","lai_q50","lai_mean","lai_q75","lai_q95","lai_max"]

def main():
    print("loading LAI panel + geometry ...",flush=True)
    df=pd.read_parquet(LAI); df["gage_id"]=df.gage_id.astype(str)
    geo=gpd.read_file(GEOM,ignore_geometry=True)
    geo["gage_id"]=geo.gage_id.astype(str)
    geo=geo[["gage_id","canon_id","latitude","longitude","watershed_geom_source","low_confidence"]].drop_duplicates("gage_id")

    # ---- month grid (shared across all gages) ----
    df["ym"]=df.year*12+(df.month-1)
    ym0,ym1=int(df.ym.min()),int(df.ym.max())
    grid=list(range(ym0,ym1+1)); nmon=len(grid)
    labels=[f"{m//12}-{m%12+1:02d}" for m in grid]
    print(f"  {len(df):,} rows | {nmon} months {labels[0]}..{labels[-1]} | {df.gage_id.nunique()} gages",flush=True)

    gages=geo.gage_id.tolist()                      # map/order driven by geometry (has coords for all 7,964)
    ngage=len(gages)

    # ---- pivot each stat to gage x month matrix on the full grid ----
    def matrix(col):
        p=df.pivot_table(index="gage_id",columns="ym",values=col,aggfunc="mean")
        p=p.reindex(index=gages,columns=grid)
        return p.to_numpy()                          # (ngage, nmon), NaN where missing
    print("  pivoting stats ...",flush=True)
    M={c:matrix(c) for c in set(TS_STATS+["lai_mean"])}
    mean_m=M["lai_mean"]                             # (ngage,nmon)
    qual={c:matrix(c) for c in ["mean_obs_valid_frac","frac_cloud","frac_snow_ice","n_eff_pixels","partial_month","n_composites"]}

    # ---- % GOOD COVERAGE per basin-month ----
    # = (valid pixel-observations) / (basin pixels x ALL expected composite dates in the month)
    # = mean_obs_valid_frac  [areal x among-present-dates]  x  (present composites / expected composites) [temporal]
    # expected composites = the normal full count for that calendar month (mode; robust to the few gap months)
    nexp={m:int(df.loc[df.month==m,"n_composites"].mode().iloc[0]) for m in range(1,13)}
    months_of=np.array([(mm%12)+1 for mm in grid])
    nexp_row=np.array([nexp[mo] for mo in months_of],dtype=float)
    with np.errstate(all="ignore"):
        temporal=np.clip(qual["n_composites"]/nexp_row[None,:],0,1)
        good_cov=np.where(np.isnan(qual["mean_obs_valid_frac"]),0.0,qual["mean_obs_valid_frac"])*temporal  # 0..1

    # ---- per-gage SUMMARY variables (map colors) ----
    print("  computing summary variables ...",flush=True)
    yrs=np.array([m//12 for m in grid]); mon=np.array([m%12 for m in grid])  # mon 0-11
    with np.errstate(all="ignore"):
        lai_mean_overall=np.nanmean(mean_m,axis=1)
        lai_range=np.nanmax(M["lai_max"],axis=1)-np.nanmin(M["lai_min"],axis=1)
        pct_na_months=np.isnan(mean_m).mean(axis=1)*100.0
        mean_valid_frac=np.nanmean(qual["mean_obs_valid_frac"],axis=1)
        mean_frac_cloud=np.nanmean(qual["frac_cloud"],axis=1)
        mean_frac_snow=np.nanmean(qual["frac_snow_ice"],axis=1)
        mean_n_eff=np.nanmean(qual["n_eff_pixels"],axis=1)
        mean_good_coverage=good_cov.mean(axis=1)*100.0
        pct_low_coverage=(good_cov<0.5).mean(axis=1)*100.0

    uyrs=np.unique(yrs); annual=np.full((ngage,len(uyrs)),np.nan)
    for j,y in enumerate(uyrs):
        cols=np.where(yrs==y)[0]; annual[:,j]=np.nanmean(mean_m[:,cols],axis=1)
    def sen_slope(a):                                # robust slope of annual means (LAI/yr); NA-tolerant
        ok=~np.isnan(a)
        if ok.sum()<5: return np.nan
        x=uyrs[ok].astype(float); y=a[ok]; sl=[]
        for p in range(len(x)):
            dx=x[p+1:]-x[p]; dy=y[p+1:]-y[p]; sl.extend((dy/dx).tolist())
        return float(np.median(sl)) if sl else np.nan
    lai_trend_slope=np.array([sen_slope(annual[i]) for i in range(ngage)])
    def first_last_delta(a):
        ok=np.where(~np.isnan(a))[0]
        return float(a[ok[-1]]-a[ok[0]]) if len(ok)>=2 else np.nan
    lai_delta_last_first=np.array([first_last_delta(annual[i]) for i in range(ngage)])
    def yr_of(a,fn):
        ok=~np.isnan(a)
        if ok.sum()<3: return np.nan
        return float(uyrs[ok][fn(a[ok])])
    year_max_lai=np.array([yr_of(annual[i],np.argmax) for i in range(ngage)])
    year_min_lai=np.array([yr_of(annual[i],np.argmin) for i in range(ngage)])

    clim=np.full((ngage,12),np.nan)
    for k in range(12):
        cols=np.where(mon==k)[0]; clim[:,k]=np.nanmean(mean_m[:,cols],axis=1)
    def clim_month(c,fn):
        if np.isnan(c).all(): return None
        return MONTHS_ABBR[int(fn(np.where(np.isnan(c),(np.inf if fn is np.argmin else -np.inf),c)))]
    peak_month=[clim_month(clim[i],np.argmax) for i in range(ngage)]
    trough_month=[clim_month(clim[i],np.argmin) for i in range(ngage)]

    # ---- assemble POINTS + VARS ----
    def r(arr,nd=3):
        return [None if (v is None or (isinstance(v,float) and not np.isfinite(v))) else round(float(v),nd) for v in arr]
    cont={
        "lai_mean_overall":(r(lai_mean_overall),"Mean monthly LAI averaged over the whole record (m2/m2)"),
        "lai_seasonal_range":(r(lai_range),"Record max(lai_max) minus record min(lai_min) — overall LAI amplitude"),
        "pct_na_months":(r(pct_na_months,1),"% of the 270 months with no valid LAI (winter snow/darkness, urban, gaps)"),
        "mean_good_coverage":(r(mean_good_coverage,1),"Mean % good coverage = valid pixel-obs / (all basin pixels x all expected composite dates); low in the snowy north (legit winter no-data)"),
        "pct_months_low_coverage":(r(pct_low_coverage,1),"% of months with <50% good coverage (poorly observed months; includes legitimate winter snow)"),
        "mean_valid_frac":(r(mean_valid_frac),"Mean fraction of land observations that were clear/usable (0-1)"),
        "mean_frac_cloud":(r(mean_frac_cloud),"Mean fraction of observations flagged cloud/cirrus/shadow (0-1)"),
        "mean_frac_snow_ice":(r(mean_frac_snow),"Mean fraction of observations flagged snow/ice (0-1)"),
        "mean_n_eff_pixels":(r(mean_n_eff,1),"Mean count of usable-LAI 500 m pixels in the basin (size/coverage proxy)"),
        "lai_trend_slope":(r(lai_trend_slope,4),"Theil-Sen slope of annual mean LAI (LAI per year; + greening, - browning)"),
        "lai_delta_last_first":(r(lai_delta_last_first),"Final-year mean LAI minus first-year-of-record mean LAI (net change over the record)"),
        "year_max_lai":(r(year_max_lai,0),"Calendar year with the highest annual mean LAI"),
        "year_min_lai":(r(year_min_lai,0),"Calendar year with the lowest annual mean LAI"),
    }
    cat={
        "peak_month":(peak_month,"Climatological month of peak LAI (greenest month, averaged over years)",MONTHS_ABBR),
        "trough_month":(trough_month,"Climatological month of minimum LAI (least-green month)",MONTHS_ABBR),
        "watershed_geom_source":(geo.watershed_geom_source.fillna("NA").tolist(),"Basin polygon source: gagesii / wsc_eccc / hydrobasins",None),
        "low_confidence":([("yes" if int(v)==1 else "no") for v in geo.low_confidence.fillna(0)],"Geometry/coverage low-confidence flag (HydroBasins fallback or area outlier)",None),
    }
    VARS=[]; V={}
    for k,(vals,desc) in cont.items():
        finite=[v for v in vals if v is not None]
        dom=[round(float(np.nanpercentile(finite,2)),3),round(float(np.nanpercentile(finite,98)),3)] if finite else [0,1]
        if dom[0]==dom[1]: dom[1]=dom[0]+1
        VARS.append({"key":k,"type":"continuous","domain":dom,"desc":desc}); V[k]=vals
    for k,(vals,desc,cats) in cat.items():
        cs=cats if cats else sorted(set(x for x in vals if x is not None))
        VARS.append({"key":k,"type":"categorical","cats":cs,"desc":desc}); V[k]=vals

    POINTS={"n":ngage,"lat":r(geo.latitude.tolist(),5),"lon":r(geo.longitude.tolist(),5),
            "gage_id":gages,"v":V}

    # ---- per-gage time-series strings (LAI*10 int; ';'-joined stats, ','-joined months, ''=NA) ----
    print("  encoding time series ...",flush=True)
    def enc_row(vec):
        return ",".join("" if np.isnan(x) else str(int(round(x*10))) for x in vec)
    stat_mats=[M[s] for s in TS_STATS]
    series=[";".join(enc_row(sm[i]) for sm in stat_mats) for i in range(ngage)]
    cov_pct=(np.where(np.isnan(good_cov),0,good_cov)*100).round().astype(int)   # per-month good-coverage %
    cov=[",".join(map(str,cov_pct[i].tolist())) for i in range(ngage)]
    TS={"t":labels,"stats":TS_STATS,"s":series,"cov":cov}

    html=TEMPLATE.replace("__POINTS__",json.dumps(POINTS)) \
                 .replace("__VARS__",json.dumps(VARS)) \
                 .replace("__TS__",json.dumps(TS))
    os.makedirs(os.path.dirname(OUT),exist_ok=True)
    with open(OUT,"w") as f: f.write(html)
    mb=os.path.getsize(OUT)/1e6
    print(f"\nwrote {OUT}  ({mb:.1f} MB)")
    print(f"  {ngage} gages | {len(VARS)} variables | {nmon}-month series x {len(TS_STATS)} stats")
    print(f"  all-NA gages (urban etc.): {int((pct_na_months>=99.9).sum())} | median %NA months: {np.nanmedian(pct_na_months):.1f}")

TEMPLATE=r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>MODIS LAI — Watershed QA/QC Explorer</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  html, body { margin:0; font-family: system-ui, Arial, sans-serif; color:#222; }
  #app { display:flex; height:100vh; }
  #side { width:340px; min-width:340px; padding:14px 16px; box-sizing:border-box;
          overflow-y:auto; border-right:1px solid #ddd; background:#fafafa; }
  #main { flex:1; display:flex; flex-direction:column; min-width:0; position:relative; }
  #map { flex:1; }
  #ts { height:300px; border-top:1px solid #ccc; background:#fff; display:none; position:relative; }
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
  #tsclose { position:absolute; top:6px; right:10px; cursor:pointer; font-size:18px; color:#888; border:none; background:none; z-index:5; }
  #tstitle { position:absolute; top:8px; left:14px; font-size:13px; font-weight:600; z-index:5; }
  #tssub { position:absolute; top:26px; left:14px; right:40px; font-size:11px; color:#666; z-index:5; }
  .hbox { position:fixed; pointer-events:none; background:#fff; border:1px solid #aaa; border-radius:4px;
          padding:6px 9px; font-size:11px; line-height:1.5; box-shadow:0 2px 8px rgba(0,0,0,.18);
          z-index:1000; display:none; max-width:260px; }
  .hbox b { color:#111; } .hbox .k { color:#777; }
  #maphint { position:absolute; top:10px; left:50%; transform:translateX(-50%); z-index:1001;
             background:#2b6cb0; color:#fff; font-size:12px; padding:6px 14px; border-radius:16px;
             box-shadow:0 2px 6px rgba(0,0,0,.25); }
  #defsbox { position:absolute; top:10px; right:10px; z-index:1001; text-align:right; }
  #defsbtn { background:#fff; border:1px solid #aaa; border-radius:5px; padding:6px 10px; font-size:12px;
             font-weight:600; cursor:pointer; box-shadow:0 1px 4px rgba(0,0,0,.18); }
  #defspanel { display:none; margin-top:6px; width:440px; max-width:46vw; max-height:72vh; overflow-y:auto;
               background:#fff; border:1px solid #aaa; border-radius:6px; box-shadow:0 4px 14px rgba(0,0,0,.22);
               padding:12px 14px; text-align:left; }
  #defspanel.show { display:block; }
  #defspanel h2 { font-size:14px; margin:0 0 4px; }
  #defspanel .dnote { color:#666; font-size:11px; margin:0 0 10px; line-height:1.4; }
  .defrow { display:flex; gap:10px; font-size:12px; padding:3px 0; border-bottom:1px solid #f0f0f0; }
  .defrow code { color:#2b6cb0; font-weight:600; min-width:150px; display:inline-block; }
  .defrow span { color:#444; }
  .defcat { font-size:12px; font-weight:700; margin:12px 0 4px; color:#333; }
</style>
</head>
<body>
<div id="app">
  <div id="side">
    <h1>MODIS LAI — Watershed QA/QC</h1>
    <p class="sub">7,964 watersheds, monthly LAI 2002–2024. Points colored by the selected
       variable; hover a point for all its values, click (or use the picker) for its time series.</p>

    <label>Find a gage</label>
    <input id="gagePick" list="gageList" placeholder="type a gage id..." />
    <datalist id="gageList"></datalist>

    <label>Filter variables</label>
    <input id="filter" type="text" placeholder="type to filter, e.g. na, month, trend..." />
    <label>Variable (color)</label>
    <select id="varSel" size="1"></select>
    <div id="vdesc"></div>

    <div id="legend"></div>
    <div id="count"></div>

    <div id="tech">
      <h3>How this LAI product is made</h3>
      <p><b>Source.</b> MODIS <b>MCD15A3H</b> (500 m, 4-day composites), 2002–2024, aggregated to each
         gage's upstream watershed.</p>
      <p><b>Identifying problematic data.</b> Every 4-day composite is screened pixel-by-pixel on the
         <code>FparLai_QC</code> band. A pixel-observation is kept only when <b>MODLAND_QC</b> = good
         (bit&nbsp;0&nbsp;=&nbsp;0) <b>and</b> <b>SCF_QC</b> (bits&nbsp;5–7) marks the <b>main radiative-transfer
         algorithm</b> (values 0–1), <b>and</b> it is a real retrieval — raw DN&nbsp;&lt;&nbsp;101 (codes
         &ge;101 are fills, incl. 250 = urban, 255 = fill), masked <b>before</b> the &times;0.1 scaling.
         Cloud / cirrus / shadow / snow-ice / aerosol flags from <code>FparExtra_QC</code> are recorded as
         diagnostic fractions (<code>frac_cloud</code>, <code>frac_snow_ice</code>, …) rather than a separate
         drop step — cloud-contaminated retrievals already fail the MODLAND / main-algorithm screen.</p>
      <p><b>Cleaning &amp; aggregation.</b> Surviving composites are combined into a day-weighted
         per-pixel <b>monthly mean</b>, then summarized across the basin's pixels with area-
         (coverage-) weighting to give the spatial mean, std and quantiles; min/max are unweighted extrema.</p>
      <p><b>Inclusion threshold &amp; coverage.</b> A month is reported when at least one QA-good composite
         covers the basin. The continuous <b>good coverage</b> metric (map variables + the red&rarr;green strip
         under each time series) makes completeness explicit: valid pixel-observations &divide; (<i>all</i> basin
         pixels &times; <i>all</i> expected composite dates), so a month reads ~5% when mostly cloudy/missing and
         ~95% when fully observed and clear. (<code>partial_month</code>, a binary &lt;6-composite flag retained
         in the product parquet, is the older version of this idea.)</p>
      <p><b>Why some months are NA but a location still has data.</b> At high latitudes, winter
         snow/darkness yields no valid LAI retrieval, so those months are <b>legitimately NA</b> (not an
         error) while the warm-season record is complete — hence northern basins show many NA months yet
         are still delivered. <b>17 urban basins</b> are NA for the entire record because MODIS does not
         retrieve LAI over built-up land.</p>
    </div>
  </div>
  <div id="main">
    <div id="map"></div>
    <div id="maphint">👆 Click any gage (or use “Find a gage” at left) to see its full LAI time series</div>
    <div id="defsbox"><button id="defsbtn">▸ Variable definitions</button><div id="defspanel"></div></div>
    <div id="ts">
      <div id="tstitle"></div><div id="tssub"></div>
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
const TS     = __TS__;
</script>
<script>
// ---------- basemap ----------
const cvs = L.canvas({ padding: 0.5 });
const map = L.map('map', { preferCanvas: true }).setView([48, -100], 4);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
  { attribution:'&copy; OpenStreetMap &copy; CARTO', subdomains:'abcd', maxZoom:19 }).addTo(map);

// ---------- color scales ----------
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

// ---------- markers + hover ----------
const hoverEl=document.getElementById('hover');
function pointInfoHTML(i){
  let h=`<b>${POINTS.gage_id[i]}</b> <span class="k">(${POINTS.lat[i].toFixed(3)}, ${POINTS.lon[i].toFixed(3)})</span>`;
  VARS.forEach(v=>{ h+=`<br><span class="k">${v.key}:</span> ${fmtv(POINTS.v[v.key][i])}`; });
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
  const meta=VARS.find(v=>v.key===key); if(!meta) return; curVar=meta;
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

// ---------- time series ----------
const NS='http://www.w3.org/2000/svg';
function parseSeries(i){                          // lazy decode: ';'-stats, ','-months, ''=NA, ints are LAI*10
  const blocks=TS.s[i].split(';'); const out={};
  TS.stats.forEach((st,k)=>{ out[st]=blocks[k].split(',').map(x=> x===''?null:(+x)/10 ); });
  out._cov = TS.cov ? TS.cov[i].split(',').map(x=>+x) : null;   // per-month % good coverage (0-100)
  return out;
}
function covColor(p){ const f=Math.max(0,Math.min(100,p))/100; return `hsl(${Math.round(120*f)},65%,45%)`; } // red->green
let curParsed=null, tsGeom=null;                  // for the per-month hover readout
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
  const d=parseSeries(i); curParsed=d; const T=TS.t, n=T.length;
  document.getElementById('ts').classList.add('show');
  const svg=document.getElementById('tssvg'); while(svg.firstChild) svg.removeChild(svg.firstChild);
  const W=svg.clientWidth||svg.parentNode.clientWidth, H=300;
  const mL=42,mR=12,mT=44,mB=46, pw=W-mL-mR, ph=H-mT-mB;
  let ymax=0; (d.lai_max||[]).forEach(v=>{ if(v!==null&&v>ymax) ymax=v; });
  if (ymax<=0) ymax=1; ymax=Math.ceil(ymax*1.05*10)/10;
  const X=k=> mL + (n<=1?0:pw*k/(n-1));
  const Y=v=> mT + ph*(1-v/ymax);
  tsGeom={mL,mT,pw,ph,n,X};
  function add(tag,at){ const e=document.createElementNS(NS,tag); for(const k in at) e.setAttribute(k,at[k]); svg.appendChild(e); return e; }
  add('line',{x1:mL,y1:mT+ph,x2:mL+pw,y2:mT+ph,stroke:'#999'});
  add('line',{x1:mL,y1:mT,x2:mL,y2:mT+ph,stroke:'#999'});
  for (let g=0; g<=4; g++){ const yv=ymax*g/4, yy=Y(yv);
    add('line',{x1:mL,y1:yy,x2:mL+pw,y2:yy,stroke:'#eee'});
    const t=add('text',{x:mL-5,y:yy+3,'text-anchor':'end','font-size':10,fill:'#666'}); t.textContent=yv.toFixed(1); }
  for (let k=0;k<n;k++){ if (T[k].endsWith('-01')){ const xx=X(k);
    add('line',{x1:xx,y1:mT+ph,x2:xx,y2:mT+ph+4,stroke:'#999'});
    const t=add('text',{x:xx,y:mT+ph+32,'text-anchor':'middle','font-size':9,fill:'#666'}); t.textContent=T[k].slice(0,4); } }
  function band(lo,hi,fill,op){ let k=0;
    while(k<n){ if(lo[k]===null||hi[k]===null){k++;continue;} let j=k; const up=[],dn=[];
      while(j<n && lo[j]!==null && hi[j]!==null){ up.push(`${X(j)},${Y(hi[j])}`); dn.push(`${X(j)},${Y(lo[j])}`); j++; }
      if (up.length>1) add('polygon',{points:up.concat(dn.reverse()).join(' '),fill:fill,opacity:op,stroke:'none'}); k=j+1; } }
  function line(arr,stroke,w){ let k=0;
    while(k<n){ if(arr[k]===null){k++;continue;} let j=k; const pts=[];
      while(j<n && arr[j]!==null){ pts.push(`${X(j)},${Y(arr[j])}`); j++; }
      if (pts.length>1) add('polyline',{points:pts.join(' '),fill:'none',stroke:stroke,'stroke-width':w});
      else if (pts.length===1){ const p=pts[0].split(','); add('circle',{cx:p[0],cy:p[1],r:1.3,fill:stroke}); }
      k=j+1; } }
  band(d.lai_min,d.lai_max,'#2b6cb0',0.12);
  band(d.lai_q05,d.lai_q95,'#2b6cb0',0.18);
  band(d.lai_q25,d.lai_q75,'#2b6cb0',0.30);
  line(d.lai_q50,'#1a4971',1);
  line(d.lai_mean,'#000',1.6);
  // coverage strip: per-month % good coverage (red=low -> green=high)
  if (d._cov){ const cy=mT+ph+6, ch=9, cw=pw/(n-1||1);
    for(let k=0;k<n;k++) add('rect',{x:X(k)-cw/2,y:cy,width:Math.max(cw,1.2),height:ch,fill:covColor(d._cov[k]),stroke:'none'});
    const t=add('text',{x:mL-5,y:cy+ch,'text-anchor':'end','font-size':8,fill:'#888'}); t.textContent='cov'; }
  // title + summary (incl. per-gage quality stats)
  const q=k=>fmtv(POINTS.v[k]?POINTS.v[k][i]:null);
  document.getElementById('tstitle').textContent=`Gage ${POINTS.gage_id[i]} — monthly LAI`;
  const nNA=(d.lai_mean||[]).filter(v=>v===null).length;
  document.getElementById('tssub').innerHTML=
    `bands: min–max / 5–95 / 25–75 pct · bold = mean · thin = median · strip = % good coverage (red→green) &nbsp;|&nbsp; ${nNA}/${n} months NA`
    + (curVar?` &nbsp;|&nbsp; ${curVar.key} = ${q(curVar.key)}`:'')
    + ` &nbsp;|&nbsp; valid_frac=${q('mean_valid_frac')} · cloud=${q('mean_frac_cloud')} · snow_ice=${q('mean_frac_snow_ice')} · n_eff_px=${q('mean_n_eff_pixels')}`;
  // top-right band swatches
  const lx=W-mR-168, ly=mT-30;
  [['min–max',0.12],['5–95',0.18],['25–75',0.30]].forEach((s,rr)=>{
    add('rect',{x:lx+rr*56,y:ly,width:10,height:10,fill:'#2b6cb0',opacity:s[1]});
    const t=add('text',{x:lx+rr*56+13,y:ly+9,'font-size':9,fill:'#666'}); t.textContent=s[0]; });
  // hover guideline + overlay
  const guide=add('line',{x1:mL,y1:mT,x2:mL,y2:mT+ph,stroke:'#d00','stroke-width':1,opacity:0,'pointer-events':'none'});
  const ov=add('rect',{x:mL,y:mT,width:pw,height:ph,fill:'transparent','pointer-events':'all'});
  ov.addEventListener('mousemove',ev=>{
    const rect=svg.getBoundingClientRect(); const px=ev.clientX-rect.left;
    let k=Math.round((px-mL)/(pw||1)*(n-1)); k=Math.max(0,Math.min(n-1,k));
    guide.setAttribute('x1',X(k)); guide.setAttribute('x2',X(k)); guide.setAttribute('opacity',0.6);
    const rows=[['mean','lai_mean'],['max','lai_max'],['q95','lai_q95'],['q75','lai_q75'],
                ['q50','lai_q50'],['q25','lai_q25'],['q05','lai_q05'],['min','lai_min']];
    let h=`<b>${T[k]}</b>`;
    rows.forEach(rw=>{ const val=curParsed[rw[1]][k]; h+=`<br><span class="k">${rw[0]}:</span> ${val===null?'NA':val.toFixed(2)}`; });
    if (curParsed._cov) h+=`<br><span class="k">good coverage:</span> <b>${curParsed._cov[k]}%</b>`;
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
  let h='<h2>Variable definitions</h2><p class="dnote">Map variables are one value per watershed. '
       +'Time-series stats are the within-watershed spatial distribution of that month\'s mean LAI.</p>'
       +'<div class="defgrid">';
  h+='<div class="defcat">Map variables (point color &amp; hover)</div>';
  VARS.forEach(v=>{ const t=v.type==='categorical'?' (categorical)':' (continuous; color clipped to 2nd–98th pct)';
    h+=`<div class="defrow"><code>${v.key}</code><span>${v.desc}${t}</span></div>`; });
  h+='<div class="defcat">Data coverage (time-series strip &amp; hover)</div>';
  h+=`<div class="defrow"><code>good coverage</code><span>Valid pixel-observations &divide; (all basin pixels &times; all expected composite dates) for the month. `
   +`100% = fully observed &amp; clear; ~0% = no usable data (cloud, snow/darkness, or missing composites). `
   +`Generalizes the binary partial-month flag into a continuous % — see the red&rarr;green strip under each chart.</span></div>`;
  const ts=[['lai_mean','Basin coverage-weighted mean of the per-pixel monthly-mean LAI'],
    ['lai_min','Spatial minimum across basin pixels (unweighted extremum — a sliver can drive it)'],
    ['lai_q05','Coverage-weighted spatial 5th percentile'],['lai_q25','Coverage-weighted spatial 25th percentile'],
    ['lai_q50','Coverage-weighted spatial median'],['lai_q75','Coverage-weighted spatial 75th percentile'],
    ['lai_q95','Coverage-weighted spatial 95th percentile'],
    ['lai_max','Spatial maximum across basin pixels (unweighted extremum)']];
  h+='<div class="defcat">Time-series stats (per month, on click)</div>';
  ts.forEach(r=>{ h+=`<div class="defrow"><code>${r[0]}</code><span>${r[1]} (m2/m2)</span></div>`; });
  el.innerHTML=h+'</div>';
}

// ---------- controls ----------
const varSel=document.getElementById('varSel'), filterEl=document.getElementById('filter');
function populateVars(){ const f=filterEl.value.toLowerCase();
  const list=VARS.filter(v=>v.key.toLowerCase().includes(f)||(v.desc||'').toLowerCase().includes(f));
  varSel.innerHTML=''; list.forEach(v=>{ const o=document.createElement('option');
    o.value=v.key; o.textContent=v.key+(v.type==='categorical'?'  (class)':''); varSel.appendChild(o); });
  varSel.size=Math.max(2,Math.min(18,list.length)); if(list.length) recolor(list[0].key); }
varSel.onchange=()=>recolor(varSel.value); filterEl.oninput=populateVars;
const dl=document.getElementById('gageList');
POINTS.gage_id.forEach(g=>{ const o=document.createElement('option'); o.value=g; dl.appendChild(o); });
const gmap={}; POINTS.gage_id.forEach((g,i)=>gmap[g]=i);
document.getElementById('gagePick').onchange=function(){ const i=gmap[this.value.trim()];
  if(i!==undefined){ map.setView([POINTS.lat[i],POINTS.lon[i]],Math.max(map.getZoom(),6)); selectGage(i); } };
const defsbtn=document.getElementById('defsbtn'), defspanel=document.getElementById('defspanel');
defsbtn.onclick=()=>{ const open=defspanel.classList.toggle('show');
  defsbtn.textContent=(open?'▾':'▸')+' Variable definitions'; };
buildDefs(); populateVars();
</script>
</body>
</html>
"""

if __name__=="__main__": main()

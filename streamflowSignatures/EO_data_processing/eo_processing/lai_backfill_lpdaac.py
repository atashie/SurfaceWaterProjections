"""LP DAAC backfill for the MODIS LAI gaps the ClimateAI catalog can't serve.

Two gap sets (see lai_always_na_basins.csv + Codex review):
  * --target farN     : 92 far-N basins (6 tiles absent from the catalog), full 2002-2024 record
  * --target m2024_11 : Nov-2024 (catalog has zero composites that month), all basins

Reads MCD15A3H HDF4 granules from NASA LP DAAC (CMR search + curl Bearer download + pyhdf,
since the catalog COG path / earthaccess-S3 path is unavailable here and pip-GDAL lacks HDF4),
then applies the SAME QC/mask/day-weighted-monthly-mean + coverage-weighted exactextract +
QA-fraction logic as lai_pipeline.py, so backfilled rows match the main table schema exactly.

Usage:
  python lai_backfill_lpdaac.py --target farN --workers 4
  python lai_backfill_lpdaac.py --target m2024_11 --workers 4
"""
import os, re, sys, glob, time, json, argparse, subprocess, datetime as dt
os.environ.setdefault("PROJ_LIB","/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA","/home/sagemaker-user/.conda/envs/geo/share/proj")
import numpy as np, pandas as pd, geopandas as gpd, rasterio, pyproj
from rasterio.transform import Affine
from concurrent.futures import ProcessPoolExecutor, as_completed

SINU="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
T=1111950.5196669998; X0=-20015109.354; YTOP=10007554.6770; PX=T/2400.0
LAI_OPS=["count","mean","stdev","min","max",
     "quantile(q=0.05)","quantile(q=0.25)","quantile(q=0.5)","quantile(q=0.75)","quantile(q=0.95)"]
LAI_RENAME={"count":"n_eff_pixels","mean":"lai_mean","stdev":"lai_std","min":"lai_min","max":"lai_max",
        "quantile_5":"lai_q05","quantile_25":"lai_q25","quantile_50":"lai_q50","quantile_75":"lai_q75","quantile_95":"lai_q95"}
QA_BANDS=["nland","nvalid","ncloud","nsnow","naero","nbackup"]
ROOT="/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures/data_out/eo_lai"
TOK=open(os.path.expanduser("~/.edl_token")).read().strip()
KEEP=["gage_id","canon_id","year","month","lai_min","lai_q05","lai_q25","lai_q50","lai_q75",
      "lai_q95","lai_max","lai_mean","lai_std","n_eff_pixels","n_basin_pixels","valid_pixel_frac",
      "mean_obs_valid_frac","frac_cloud","frac_snow_ice","frac_aerosol","frac_backup_algorithm",
      "n_composites","low_confidence"]

def tile_center_lonlat(tile):
    h=int(tile[1:3]); v=int(tile[4:6])
    cx=X0+(h+0.5)*T; cy=YTOP-(v+0.5)*T
    return pyproj.Transformer.from_crs(SINU,"EPSG:4326",always_xy=True).transform(cx,cy)

def tile_transform(tile):
    h=int(tile[1:3]); v=int(tile[4:6])
    return Affine(PX,0,X0+h*T,0,-PX,YTOP-v*T)

def build_granule_cache(tiles,t0,t1,cache_path):
    if os.path.exists(cache_path): return pd.read_parquet(cache_path)
    import earthaccess
    rows=[]
    for tid in tiles:
        lon,lat=tile_center_lonlat(tid)
        r=earthaccess.search_data(short_name="MCD15A3H",version="061",temporal=(t0,t1),
                                  bounding_box=(lon-0.01,lat-0.01,lon+0.01,lat+0.01))
        n=0
        for g in r:
            for u in g.data_links():
                if not (u.endswith(".hdf") and f".{tid}." in u): continue
                m=re.search(r"\.A(\d{4})(\d{3})\.",u)
                d=(dt.date(int(m.group(1)),1,1)+dt.timedelta(days=int(m.group(2))-1)).isoformat()
                rows.append((tid,d,u)); n+=1
        print(f"  {tid}: {n} granules",flush=True)
    c=pd.DataFrame(rows,columns=["tile","date","uri"]).drop_duplicates()
    c.to_parquet(cache_path); return c

def comp_weight(cdate,mstart,mend):
    c0=dt.date.fromisoformat(cdate); c1=c0+dt.timedelta(days=4)
    return max(0,(min(c1,mend)-max(c0,mstart)).days)

def _dl(uri,dst,n=7):
    """Download a granule, lock-aware. Earthdata locks the account on request bursts
    ('invalid_account_status ... currently locked'); on a non-HDF response we back off
    (30s,60s,...,up to 10min) instead of hammering. Gentle throttle on success."""
    for i in range(n):
        rc=subprocess.run(["curl","-sL","-H",f"Authorization: Bearer {TOK}",uri,"-o",dst,"--max-time","300"]).returncode
        if rc==0 and os.path.exists(dst) and os.path.getsize(dst)>1000:
            with open(dst,"rb") as f: magic=f.read(4)
            if magic==b"\x0e\x03\x13\x01":
                time.sleep(0.3); return True
            head=""
            try: head=open(dst,"rb").read(400).decode("utf-8","ignore").lower()
            except Exception: pass
            if os.path.exists(dst): os.remove(dst)
            if "locked" in head or "invalid_account" in head or "rate" in head or "error" in head:
                time.sleep(min(600,30*(2**i)))   # lock/auth: exponential backoff
                continue
        if os.path.exists(dst): os.remove(dst)
        time.sleep(min(120,5*(i+1)))             # transient/network backoff
    return False

def _read_hdf(path):
    from pyhdf.SD import SD,SDC
    sd=SD(path,SDC.READ)
    lai=sd.select("Lai_500m")[:]; qc=sd.select("FparLai_QC")[:]; xqc=sd.select("FparExtra_QC")[:]
    sd.end(); return lai,qc,xqc

def _w(path,arr,tr,dtype,nodata):
    prof=dict(driver="GTiff",height=arr.shape[0],width=arr.shape[1],count=1,dtype=dtype,
              crs=SINU,transform=tr,nodata=nodata,compress="deflate")
    with rasterio.open(path,"w",**prof) as ds: ds.write(arr.astype(dtype),1)

def monthly_mean_tile(tid,tcache,year,month,tmpdir):
    mstart=dt.date(year,month,1); mend=dt.date(year+(month==12),(month%12)+1,1)
    g=tcache.copy(); g["w"]=g["date"].map(lambda d: comp_weight(d,mstart,mend)); g=g[g.w>0]
    if g.empty: return False,0
    tr=tile_transform(tid)
    acc=wsum=None; cnt={b:None for b in QA_BANDS}; ncomp=0
    for _,row in g.iterrows():
        hdf=os.path.join(tmpdir,os.path.basename(row.uri))
        if not _dl(row.uri,hdf): continue
        try: lai,qc,xqc=_read_hdf(hdf)
        except Exception:
            if os.path.exists(hdf): os.remove(hdf)
            continue
        os.remove(hdf)
        v=lai.astype(np.float32); w=row.w
        scf=(qc>>5)&7
        keep=((qc&1)==0)&(scf<=1)&(v<101); backup=(scf>=2)&(scf<=3)
        land=((xqc&3)==0)|(v<101)
        cloud=(((xqc>>4)&1)|((xqc>>5)&1)|((xqc>>6)&1)).astype(bool)
        snow=((xqc>>2)&1).astype(bool); aero=((xqc>>3)&1).astype(bool)
        fl={"nland":land,"nvalid":keep,"ncloud":cloud&land,"nsnow":snow&land,"naero":aero&land,"nbackup":backup&land}
        if acc is None:
            acc=np.zeros(v.shape,np.float32); wsum=np.zeros(v.shape,np.float32)
            cnt={b:np.zeros(v.shape,np.uint16) for b in QA_BANDS}
        for b in QA_BANDS: cnt[b]+=fl[b].astype(np.uint16)
        val=np.where(keep,v*0.1,np.nan); m=np.isfinite(val); acc[m]+=w*val[m]; wsum[m]+=w; ncomp+=1
    if acc is None: return False,0
    mm=np.where(wsum>0,acc/np.maximum(wsum,1e-9),np.nan).astype(np.float32)
    _w(os.path.join(tmpdir,f"{tid}_mean.tif"),mm,tr,"float32",float("nan"))
    for b in QA_BANDS: _w(os.path.join(tmpdir,f"{tid}_{b}.tif"),cnt[b],tr,"uint16",0)
    return True,ncomp

def process_month(args):
    year,month,cache_path,basins_path,tiles,outdir=args
    outp=os.path.join(outdir,f"lai_{year}_{month:02d}.parquet")
    if os.path.exists(outp): return outp,"cached"
    try:
        from exactextract import exact_extract
        from osgeo import gdal
        cache=pd.read_parquet(cache_path); basins=gpd.read_parquet(basins_path)
        tmpdir=os.path.join(outdir,f"_tmp_{year}_{month:02d}"); os.makedirs(tmpdir,exist_ok=True)
        ncomp={}
        for tid in tiles:
            tc=cache[cache.tile==tid]
            if tc.empty: continue
            ok,nc=monthly_mean_tile(tid,tc,year,month,tmpdir)
            if ok: ncomp[tid]=nc
        if not ncomp:
            pd.DataFrame(columns=["gage_id"]).to_parquet(outp); return outp,"empty"
        def vrt(suf):
            tifs=sorted(glob.glob(os.path.join(tmpdir,f"*_{suf}.tif")))
            vp=os.path.join(tmpdir,f"m_{suf}.vrt"); gdal.BuildVRT(vp,tifs); return vp
        res=exact_extract(vrt("mean"),basins,LAI_OPS,include_cols=["gage_id","canon_id","low_confidence"],output="pandas").rename(columns=LAI_RENAME)
        def qa(suf,count=False):
            return exact_extract(vrt(suf),basins,(["count","sum"] if count else ["sum"]),include_cols=["gage_id"],output="pandas")
        res=res.merge(qa("nland",count=True).rename(columns={"count":"n_basin_pixels","sum":"sum_land"}),on="gage_id",how="left")
        for suf,col in [("nvalid","sum_valid"),("ncloud","sum_cloud"),("nsnow","sum_snow"),("naero","sum_aero"),("nbackup","sum_backup")]:
            res=res.merge(qa(suf).rename(columns={"sum":col}),on="gage_id",how="left")
        res["year"]=year; res["month"]=month
        res["n_composites"]=int(round(np.mean(list(ncomp.values())))) if ncomp else 0
        dp=res["n_basin_pixels"].replace(0,np.nan); do=res["sum_land"].replace(0,np.nan)
        res["valid_pixel_frac"]=res["n_eff_pixels"]/dp
        res["mean_obs_valid_frac"]=res["sum_valid"]/do
        res["frac_cloud"]=res["sum_cloud"]/do; res["frac_snow_ice"]=res["sum_snow"]/do
        res["frac_aerosol"]=res["sum_aero"]/do; res["frac_backup_algorithm"]=res["sum_backup"]/do
        res[KEEP].to_parquet(outp)
        for f in glob.glob(os.path.join(tmpdir,"*")): os.remove(f)
        os.rmdir(tmpdir)
        return outp,f"{int(res['lai_mean'].notna().sum())} basins"
    except Exception as e:
        return outp,f"ERROR {year}-{month:02d}: {type(e).__name__} {str(e)[:90]}"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--target",required=True,choices=["farN","m2024_11"])
    ap.add_argument("--workers",type=int,default=2)  # keep LOW: Earthdata locks on request bursts
    a=ap.parse_args()
    geo=gpd.read_file(f"{ROOT}/geometry_7964.gpkg")
    if a.target=="farN":
        outdir=f"{ROOT}/backfill_farN"
        tiles=["h08v03","h09v03","h10v02","h11v02","h12v01","h12v02"]  # tiles absent from the ClimateAI catalog
        # the 92 far-N basins (catalog had no coverage) — committed manifest beside this module, reproducible
        man=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"lai_backfill_farN_manifest.csv"),dtype=str)
        gids=set(man.gage_id)
        ym=[(y,m) for y in range(2002,2025) for m in range(1,13) if not(y==2002 and m<7)]
        t0,t1="2002-07-01","2024-12-31"
    else:
        outdir=f"{ROOT}/backfill_2024_11"; tiles=None  # all tiles touched by all basins
        gids=set(geo.gage_id); ym=[(2024,11)]; t0,t1="2024-11-01","2024-11-30"
    os.makedirs(outdir,exist_ok=True)
    g=geo[geo.gage_id.isin(gids)].to_crs(SINU)[["gage_id","canon_id","low_confidence","geometry"]]
    if tiles is None:
        b=g.geometry.bounds; ts=set()
        for r in b.itertuples():
            for h in range(int((r.minx-X0)//T),int((r.maxx-X0)//T)+1):
                for v in range(int((YTOP-r.maxy)//T),int((YTOP-r.miny)//T)+1): ts.add(f"h{h:02d}v{v:02d}")
        tiles=sorted(ts)
    basins_path=os.path.join(outdir,"_basins_sinu.parquet"); g.to_parquet(basins_path)
    print(f"target {a.target} | basins {len(g)} | tiles {len(tiles)} {tiles} | months {len(ym)} | workers {a.workers}",flush=True)
    cache_path=os.path.join(outdir,"_granule_cache.parquet")
    t=time.time(); build_granule_cache(tiles,t0,t1,cache_path); print(f"granule cache: {time.time()-t:.0f}s",flush=True)
    tasks=[(y,m,cache_path,basins_path,tiles,outdir) for y,m in ym]
    t=time.time()
    if a.workers>1:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for f in as_completed([ex.submit(process_month,x) for x in tasks]):
                p,msg=f.result(); print(f"  {os.path.basename(p)}: {msg}",flush=True)
    else:
        for x in tasks:
            p,msg=process_month(x); print(f"  {os.path.basename(p)}: {msg}",flush=True)
    parts=[pd.read_parquet(p) for p in sorted(glob.glob(os.path.join(outdir,"lai_*.parquet"))) if pd.read_parquet(p).shape[1]>1]
    if parts:
        out=pd.concat(parts,ignore_index=True).sort_values(["gage_id","year","month"])
        out.to_parquet(os.path.join(outdir,"backfill_lai_monthly.parquet"))
        print(f"FINAL {a.target}: {len(out)} rows, {out.gage_id.nunique()} gages, {out[['year','month']].drop_duplicates().shape[0]} months",flush=True)

if __name__=="__main__": main()

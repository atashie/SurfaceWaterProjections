"""Tile-batched per-watershed monthly MODIS LAI (MCD15A3H) pipeline.

Reads MCD15A3H B02 (Lai_500m) + B03 (FparLai_QC) + B04 (FparExtra_QC) COGs from the
ClimateAI Geospatial Asset Catalog (GEOLATELY, us-east-2), QC+fill-masks, builds a
day-weighted per-pixel monthly-mean LAI field per MODIS sinusoidal tile plus per-pixel
QA counters (land/valid/cloud/snow/aerosol/backup observation counts), mosaics the tiles
for each month into seamless VRTs, and runs coverage-weighted exactextract over every
watershed (so basins spanning tile boundaries are handled correctly).
Output: one row per (gage_id, year, month) with LAI stats + data-quality columns.

Usage:
  python lai_pipeline.py --months 2020-07                 # one month (validation)
  python lai_pipeline.py --years 2002-2024 --workers 3    # full record
"""
import os, re, sys, glob, time, argparse, datetime as dt
os.environ.setdefault("PROJ_LIB","/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA","/home/sagemaker-user/.conda/envs/geo/share/proj")
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from concurrent.futures import ProcessPoolExecutor, as_completed

SINU="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
T=1111950.5196669998; X0=-20015109.354; YTOP=10007554.6770
LAI_OPS=["count","mean","stdev","min","max",
     "quantile(q=0.05)","quantile(q=0.25)","quantile(q=0.5)","quantile(q=0.75)","quantile(q=0.95)"]
LAI_RENAME={"count":"n_eff_pixels","mean":"lai_mean","stdev":"lai_std","min":"lai_min","max":"lai_max",
        "quantile_5":"lai_q05","quantile_25":"lai_q25","quantile_50":"lai_q50","quantile_75":"lai_q75","quantile_95":"lai_q95"}
# per-pixel QA counter bands written per tile (uint8; #composites/month <= ~9)
QA_BANDS=["nland","nvalid","ncloud","nsnow","naero","nbackup"]

def creds():
    import boto3
    fc=boto3.Session().get_credentials().get_frozen_credentials()
    os.environ.update(AWS_ACCESS_KEY_ID=fc.access_key, AWS_SECRET_ACCESS_KEY=fc.secret_key)
    if fc.token: os.environ["AWS_SESSION_TOKEN"]=fc.token

def basin_tiles(gs):
    b=gs.geometry.bounds; tiles=set()
    for r in b.itertuples():
        h0=int((r.minx-X0)//T); h1=int((r.maxx-X0)//T)
        v0=int((YTOP-r.maxy)//T); v1=int((YTOP-r.miny)//T)
        for h in range(h0,h1+1):
            for v in range(v0,v1+1): tiles.add((h,v))
    return sorted(tiles)

def build_uri_cache(tiles, start, end, cache_path):
    """Query the catalog once per tile (full range); cache {tile,granule,date,band,uri} to parquet."""
    if os.path.exists(cache_path):
        return pd.read_parquet(cache_path)
    creds()
    from climateai.geolately.db_utils import Geospatial_DB
    from shapely.geometry import box
    db=Geospatial_DB(); rows=[]
    for (h,v) in tiles:
        xmin=X0+h*T; ymax=YTOP-v*T
        poly=gpd.GeoDataFrame(geometry=[box(xmin,ymax-T,xmin+T,ymax)],crs=SINU).to_crs(4326)
        df=db.query_asset_catalog(product="MCD15A3H",start_date=start,end_date=end,query_gdf=poly)
        tid=f"h{h:02d}v{v:02d}"
        if df is None or len(df)==0:
            print(f"  skip {tid}: no catalog data", flush=True); continue
        df=df[df["URI"].str.contains(tid)]
        if len(df)==0:
            print(f"  skip {tid}: no {tid} URIs in result", flush=True); continue
        for u in df["URI"]:
            key=re.sub(r"_B0\d\.tif$","",u); band=u[-7:-4]
            rows.append((tid, key, str(df.loc[df.URI==u,"START_DATE"].iloc[0])[:10], band, u))
        print(f"  cached {tid}: {df['URI'].str.endswith('_B02.tif').sum()} B02", flush=True)
    cache=pd.DataFrame(rows, columns=["tile","granule","date","band","uri"])
    cache.to_parquet(cache_path); return cache

def comp_weight(cdate, mstart, mend):
    c0=dt.date.fromisoformat(cdate); c1=c0+dt.timedelta(days=4)
    return max(0,(min(c1,mend)-max(c0,mstart)).days)

def _read_retry(uri, n=4):
    for _ in range(n):
        try:
            with rasterio.open(uri) as s: return s.read(1), s.profile
        except Exception: continue
    return None, None

def _write(path, arr, prof, dtype, nodata):
    p=prof.copy(); p.update(count=1,dtype=dtype,nodata=nodata,driver="GTiff",compress="deflate")
    with rasterio.open(path,"w",**p) as ds: ds.write(arr.astype(dtype),1)

def monthly_mean_tile(tile_cache, year, month, tmpdir):
    """Day-weighted per-pixel monthly-mean LAI + per-pixel QA counters for one tile.
    Writes {tile}_mean.tif (float32) and {tile}_<band>.tif (uint8) for each QA band.
    Returns (ok, ncomp)."""
    mstart=dt.date(year,month,1); mend=dt.date(year+(month==12),(month%12)+1,1)
    g=tile_cache.copy(); g["w"]=g["date"].map(lambda d: comp_weight(d,mstart,mend))
    g=g[g["w"]>0]
    if g.empty: return False,0
    acc=wsum=prof=None; ncomp=0
    cnt={b:None for b in QA_BANDS}
    for gran,grp in g.groupby("granule"):
        b02=grp.loc[grp.band=="B02","uri"]; b03=grp.loc[grp.band=="B03","uri"]; b04=grp.loc[grp.band=="B04","uri"]
        if b02.empty or b03.empty: continue
        w=grp["w"].iloc[0]
        v,pr=_read_retry(b02.iloc[0])
        if v is None: continue
        qc,_=_read_retry(b03.iloc[0])
        if qc is None: continue
        v=v.astype(np.float32)
        qce=None
        if not b04.empty: qce,_=_read_retry(b04.iloc[0])
        if prof is None: prof=pr
        scf=(qc>>5)&7
        keep=((qc&1)==0)&(scf<=1)&(v<101)              # MODLAND good & main RT algo & not fill
        backup=(scf>=2)&(scf<=3)                         # empirical/backup retrieval
        if qce is not None:
            land=((qce&3)==0)|(v<101)                   # land per LandSea OR any non-fill Lai (valid⊆land)
            cloud=(((qce>>4)&1)|((qce>>5)&1)|((qce>>6)&1)).astype(bool)  # cirrus|cloud|shadow
            snow=((qce>>2)&1).astype(bool); aero=((qce>>3)&1).astype(bool)
        else:
            land=(v<101); cloud=snow=aero=np.zeros(v.shape,bool)
        flags={"nland":land,"nvalid":keep,"ncloud":cloud&land,"nsnow":snow&land,
               "naero":aero&land,"nbackup":backup&land}
        if acc is None:
            acc=np.zeros(v.shape,np.float32); wsum=np.zeros(v.shape,np.float32)
            cnt={b:np.zeros(v.shape,np.uint16) for b in QA_BANDS}
        for b in QA_BANDS: cnt[b]+=flags[b].astype(np.uint16)
        val=np.where(keep, v*0.1, np.nan)
        m=np.isfinite(val); acc[m]+=w*val[m]; wsum[m]+=w; ncomp+=1
    if acc is None: return False,0
    mm=np.where(wsum>0, acc/np.maximum(wsum,1e-9), np.nan).astype(np.float32)
    tid=tile_cache['tile'].iloc[0]
    _write(os.path.join(tmpdir,f"{tid}_mean.tif"), mm, prof, "float32", float("nan"))
    for b in QA_BANDS:
        _write(os.path.join(tmpdir,f"{tid}_{b}.tif"), cnt[b], prof, "uint16", 0)
    return True,ncomp

def process_month(args):
    year,month,cache_path,basins_path,outdir = args
    outp=os.path.join(outdir,f"lai_{year}_{month:02d}.parquet")
    if os.path.exists(outp): return outp,"cached"
    try:
        import boto3
        from rasterio.session import AWSSession
        from exactextract import exact_extract
        from osgeo import gdal
        # SageMaker role rotates ~30 min. creds() materialized STATIC env creds before the
        # workers forked; inherited static creds expire mid-run -> S3 reads auth-fail -> empty
        # months. Drop them so boto3 resolves the auto-refreshing role provider; a fresh
        # Session+AWSSession per month then always uses current creds.
        for _k in ("AWS_ACCESS_KEY_ID","AWS_SECRET_ACCESS_KEY","AWS_SESSION_TOKEN"): os.environ.pop(_k,None)
        cache=pd.read_parquet(cache_path)
        basins=gpd.read_parquet(basins_path)        # SINU, has gage_id,canon_id,low_confidence
        tmpdir=os.path.join(outdir,f"_tmp_{year}_{month:02d}"); os.makedirs(tmpdir,exist_ok=True)
        ncomp_by_tile={}
        with rasterio.Env(AWSSession(boto3.Session())):
            for tid,tc in cache.groupby("tile"):
                ok,nc=monthly_mean_tile(tc,year,month,tmpdir)
                if ok: ncomp_by_tile[tid]=nc
        if not ncomp_by_tile:
            pd.DataFrame(columns=["gage_id"]).to_parquet(outp); return outp,"empty"
        def vrt(suffix):
            tifs=sorted(glob.glob(os.path.join(tmpdir,f"*_{suffix}.tif")))
            vp=os.path.join(tmpdir,f"m_{suffix}.vrt"); gdal.BuildVRT(vp,tifs); return vp
        # LAI distribution stats (count here = n_eff_pixels = pixels with usable monthly LAI)
        res=exact_extract(vrt("mean"), basins, LAI_OPS,
                          include_cols=["gage_id","canon_id","low_confidence"], output="pandas").rename(columns=LAI_RENAME)
        # QA counter bands -> coverage-weighted sums (and count of land pixels = denominator)
        def qa(suffix, count=False):
            d=exact_extract(vrt(suffix), basins, (["count","sum"] if count else ["sum"]),
                            include_cols=["gage_id"], output="pandas")
            return d
        nl=qa("nland",count=True).rename(columns={"count":"n_basin_pixels","sum":"sum_land"})
        res=res.merge(nl,on="gage_id",how="left")
        for suf,col in [("nvalid","sum_valid"),("ncloud","sum_cloud"),("nsnow","sum_snow"),
                        ("naero","sum_aero"),("nbackup","sum_backup")]:
            res=res.merge(qa(suf).rename(columns={"sum":col}),on="gage_id",how="left")
        res["year"]=year; res["month"]=month
        res["n_composites"]=int(round(np.mean(list(ncomp_by_tile.values())))) if ncomp_by_tile else 0
        den_pix=res["n_basin_pixels"].replace(0,np.nan)
        den_obs=res["sum_land"].replace(0,np.nan)
        res["valid_pixel_frac"]=res["n_eff_pixels"]/den_pix
        res["mean_obs_valid_frac"]=res["sum_valid"]/den_obs
        res["frac_cloud"]=res["sum_cloud"]/den_obs
        res["frac_snow_ice"]=res["sum_snow"]/den_obs
        res["frac_aerosol"]=res["sum_aero"]/den_obs
        res["frac_backup_algorithm"]=res["sum_backup"]/den_obs
        keep=["gage_id","canon_id","year","month",
              "lai_min","lai_q05","lai_q25","lai_q50","lai_q75","lai_q95","lai_max","lai_mean","lai_std",
              "n_eff_pixels","n_basin_pixels","valid_pixel_frac","mean_obs_valid_frac",
              "frac_cloud","frac_snow_ice","frac_aerosol","frac_backup_algorithm",
              "n_composites","low_confidence"]
        res[keep].to_parquet(outp)
        for f in glob.glob(os.path.join(tmpdir,"*")): os.remove(f)
        os.rmdir(tmpdir)
        return outp,f"{int(res['lai_mean'].notna().sum())} basins, valid_frac~{res['valid_pixel_frac'].median():.2f}"
    except Exception as e:
        return outp, f"ERROR {year}-{month:02d}: {type(e).__name__} {str(e)[:90]}"

def main():
    ROOT="/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures/data_out/eo_lai"
    ap=argparse.ArgumentParser()
    ap.add_argument("--years"); ap.add_argument("--months"); ap.add_argument("--workers",type=int,default=1)
    ap.add_argument("--geometry",default=f"{ROOT}/geometry_7964.gpkg")
    ap.add_argument("--outdir",default=f"{ROOT}/out_qa")
    a=ap.parse_args()
    os.makedirs(a.outdir,exist_ok=True)
    ym=[]
    if a.months:
        for tok in a.months.split(","):
            y,m=tok.split("-"); ym.append((int(y),int(m)))
    else:
        y0,y1=(int(x) for x in a.years.split("-")) if "-" in a.years else (int(a.years),int(a.years))
        for y in range(y0,y1+1):
            for m in range(1,13):
                if y==2002 and m<7: continue          # MCD15A3H combined starts ~2002-07
                ym.append((y,m))
    g=gpd.read_file(a.geometry).to_crs(SINU)[["gage_id","canon_id","low_confidence","geometry"]]
    basins_path=os.path.join(a.outdir,"_basins_sinu.parquet"); g.to_parquet(basins_path)
    tiles=basin_tiles(g)
    print(f"basins {len(g)} | tiles {len(tiles)} | months {len(ym)} | workers {a.workers}",flush=True)
    start=f"{min(y for y,_ in ym)}-01-01"; end=f"{max(y for y,_ in ym)}-12-31"
    cache_path=os.path.join(a.outdir,"_uri_cache.parquet")
    t0=time.time(); build_uri_cache(tiles,start,end,cache_path); print(f"uri cache: {time.time()-t0:.1f}s",flush=True)
    tasks=[(y,m,cache_path,basins_path,a.outdir) for y,m in ym]
    t0=time.time()
    if a.workers>1:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for f in as_completed([ex.submit(process_month,t) for t in tasks]):
                p,msg=f.result(); print(f"  {os.path.basename(p)}: {msg}",flush=True)
    else:
        for t in tasks:
            p,msg=process_month(t); print(f"  {os.path.basename(p)}: {msg}",flush=True)
    print(f"processing: {(time.time()-t0)/60:.1f} min",flush=True)
    parts=[pd.read_parquet(p) for p in sorted(glob.glob(os.path.join(a.outdir,"lai_*.parquet"))) if pd.read_parquet(p).shape[1]>1]
    out=pd.concat(parts,ignore_index=True).sort_values(["gage_id","year","month"])
    out.to_parquet(os.path.join(a.outdir,"watershed_modis_lai_monthly.parquet"))
    out.to_csv(os.path.join(a.outdir,"watershed_modis_lai_monthly.csv"),index=False)
    print(f"FINAL: {len(out)} rows, {out['gage_id'].nunique()} gages, {out[['year','month']].drop_duplicates().shape[0]} months",flush=True)

if __name__=="__main__": main()

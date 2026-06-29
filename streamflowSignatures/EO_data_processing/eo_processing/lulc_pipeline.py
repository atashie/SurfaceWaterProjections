"""Per-watershed annual MODIS LULC (MCD12Q1 v061) % class coverage.

For every QA-passing watershed (geometry_7964.gpkg), per year 2001-2024: download the MCD12Q1
HDF4 granule for each intersecting MODIS tile from NASA LP DAAC (CMR search + curl Bearer +
pyhdf, same path as the LAI backfill — pip-GDAL has no HDF4 driver), write each LC band to a
sinusoidal GeoTIFF, VRT-mosaic the tiles per band (lossless for aligned categorical MODIS), and
run coverage-weighted exactextract `unique`/`frac` to get per-class area fractions. Fractions are
reindexed to the committed class manifest (lulc_legends.csv) -> fixed ~102-col wide schema.
Output: one row per (gage_id, year); columns {band}_c{code}_pct (% over valid non-fill pixels).

Usage:
  python lulc_pipeline.py --years 2020        # one year (validation)
  python lulc_pipeline.py --years 2001-2024 --workers 2
"""
import os, re, glob, time, argparse, subprocess, datetime as dt
os.environ.setdefault("PROJ_LIB","/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA","/home/sagemaker-user/.conda/envs/geo/share/proj")
import numpy as np, pandas as pd, geopandas as gpd, rasterio, pyproj
from rasterio.transform import Affine
from concurrent.futures import ProcessPoolExecutor, as_completed

SINU="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
T=1111950.5196669998; X0=-20015109.354; YTOP=10007554.6770; PX=T/2400.0
ROOT="/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures/data_out/eo_lai"
HERE=os.path.dirname(os.path.abspath(__file__))
TOK=open(os.path.expanduser("~/.edl_token")).read().strip()
LEG=pd.read_csv(os.path.join(HERE,"lulc_legends.csv"))
BANDS=list(LEG.band.unique())
BAND_CODES={b:LEG[LEG.band==b].code.tolist() for b in BANDS}
COL={ (r.band,r.code):r.output_column for r in LEG.itertuples() }

def tile_center_lonlat(tile):
    h=int(tile[1:3]); v=int(tile[4:6])
    return pyproj.Transformer.from_crs(SINU,"EPSG:4326",always_xy=True).transform(X0+(h+0.5)*T, YTOP-(v+0.5)*T)
def tile_transform(tile):
    h=int(tile[1:3]); v=int(tile[4:6]); return Affine(PX,0,X0+h*T,0,-PX,YTOP-v*T)
def basin_tiles(gs):
    b=gs.geometry.bounds; ts=set()
    for r in b.itertuples():
        for h in range(int((r.minx-X0)//T),int((r.maxx-X0)//T)+1):
            for v in range(int((YTOP-r.maxy)//T),int((YTOP-r.miny)//T)+1): ts.add(f"h{h:02d}v{v:02d}")
    return sorted(ts)

def build_granule_cache(tiles,cache_path):
    if os.path.exists(cache_path): return pd.read_parquet(cache_path)
    import earthaccess, time
    rows=[]
    for tid in tiles:
        lon,lat=tile_center_lonlat(tid)
        r=None
        for attempt in range(5):                 # Earthdata locks on rapid search bursts -> retry w/ backoff
            try:
                r=earthaccess.search_data(short_name="MCD12Q1",version="061",temporal=("2001-01-01","2024-12-31"),
                                          bounding_box=(lon-0.01,lat-0.01,lon+0.01,lat+0.01))
                if r is not None and len(r)>0: break
            except Exception: pass
            time.sleep(min(300,20*(2**attempt)))
        if not r: print(f"  {tid}: 0 granules (search failed/empty)",flush=True); continue
        n=0
        for g in r:
            for u in g.data_links():
                if not (u.endswith(".hdf") and f".{tid}." in u): continue
                yr=int(re.search(r"\.A(\d{4})\d{3}\.",u).group(1)); rows.append((tid,yr,u)); n+=1
        print(f"  {tid}: {n} granules",flush=True)
    c=pd.DataFrame(rows,columns=["tile","year","uri"]).drop_duplicates()
    c.to_parquet(cache_path); return c

def _dl(uri,dst,n=7):
    import time
    for i in range(n):
        rc=subprocess.run(["curl","-sL","-H",f"Authorization: Bearer {TOK}",uri,"-o",dst,"--max-time","300"]).returncode
        if rc==0 and os.path.exists(dst) and os.path.getsize(dst)>1000:
            with open(dst,"rb") as f: magic=f.read(4)
            if magic==b"\x0e\x03\x13\x01": time.sleep(0.3); return True
            head=""
            try: head=open(dst,"rb").read(400).decode("utf-8","ignore").lower()
            except Exception: pass
            if os.path.exists(dst): os.remove(dst)
            if any(s in head for s in ("locked","invalid_account","rate","error")):
                time.sleep(min(600,30*(2**i))); continue
        if os.path.exists(dst): os.remove(dst)
        time.sleep(min(120,5*(i+1)))
    return False

def _write_band(path,arr,tr):
    prof=dict(driver="GTiff",height=arr.shape[0],width=arr.shape[1],count=1,dtype="uint8",
              crs=SINU,transform=tr,nodata=255,compress="deflate")
    with rasterio.open(path,"w",**prof) as ds: ds.write(arr.astype("uint8"),1)

def process_year(args):
    year,cache_path,basins_path,outdir=args
    outp=os.path.join(outdir,f"lulc_{year}.parquet")
    if os.path.exists(outp): return outp,"cached"
    try:
        from pyhdf.SD import SD,SDC
        from exactextract import exact_extract
        from osgeo import gdal
        cache=pd.read_parquet(cache_path); basins=gpd.read_parquet(basins_path)
        tmpdir=os.path.join(outdir,f"_tmp_{year}"); os.makedirs(tmpdir,exist_ok=True)
        ntile=0
        for tid,sub in cache[cache.year==year].groupby("tile"):
            hdf=os.path.join(tmpdir,f"{tid}.hdf")
            if not _dl(sub.uri.iloc[0],hdf): continue
            ok=False
            try:
                sd=SD(hdf,SDC.READ); tr=tile_transform(tid)
                for b in BANDS: _write_band(os.path.join(tmpdir,f"{tid}_{b}.tif"), sd.select(b)[:], tr)
                sd.end()
                ok=all(os.path.exists(os.path.join(tmpdir,f"{tid}_{b}.tif")) for b in BANDS)  # all 8 bands present
            except Exception as e:
                print(f"    {year} {tid}: decode/write failed: {type(e).__name__} {str(e)[:60]}",flush=True)
                for b in BANDS:                          # drop any partial band TIFFs so vrt() can't mosaic a truncated tile
                    p=os.path.join(tmpdir,f"{tid}_{b}.tif")
                    if os.path.exists(p): os.remove(p)
            finally:
                if os.path.exists(hdf): os.remove(hdf)
            if ok: ntile+=1                              # count a tile ONLY if all 8 bands decoded+written (else re-run fills)
        req=int(cache[cache.year==year].tile.nunique())
        if ntile==0:
            return outp,f"INCOMPLETE 0/{req} tiles (all downloads failed) -> skipped (re-run to fill)"
        if ntile<req:
            # lock/download gaps -> do NOT write (skip-if-exists would lock in a partial year); leave for re-run
            for f in glob.glob(os.path.join(tmpdir,"*")):
                try: os.remove(f)
                except OSError: pass
            try: os.rmdir(tmpdir)
            except OSError: pass
            return outp,f"INCOMPLETE {ntile}/{req} tiles (locks) -> skipped (re-run to fill)"
        def vrt(b):
            tifs=sorted(glob.glob(os.path.join(tmpdir,f"*_{b}.tif")))
            assert len(tifs)==ntile, f"{b} year {year}: {len(tifs)} tile TIFFs != {ntile} complete tiles"
            vp=os.path.join(tmpdir,f"m_{b}.vrt")
            # explicit nearest + nodata so the categorical mosaic is provably lossless (no resampling/averaging)
            gdal.BuildVRT(vp,tifs,options=gdal.BuildVRTOptions(resampleAlg="nearest",srcNodata=255,VRTNodata=255))
            return vp
        out=basins[["gage_id","canon_id"]].copy(); out["year"]=year
        for bi,b in enumerate(BANDS):
            ops=["unique","frac"]+(["count"] if bi==0 else [])
            res=exact_extract(vrt(b),basins,ops,include_cols=["gage_id"],output="pandas")
            colmap={c:COL[(b,c)] for c in BAND_CODES[b]}
            recs=[]; unknown=set()
            for gid,u,fr in zip(res["gage_id"],res["unique"],res["frac"]):
                for code,f in zip(u,fr):
                    code=int(code)
                    if code==255: continue                # fill/NoData
                    if code in colmap: recs.append((gid,colmap[code],float(f)*100.0))
                    else: unknown.add(code)               # real class absent from manifest -> would bias band total low
            if unknown:                                   # hard-fail (caught by outer except -> year skipped, logged) rather than silently drop
                raise ValueError(f"{b} year {year}: class codes not in lulc_legends.csv manifest: {sorted(unknown)}")
            if recs:
                long=pd.DataFrame(recs,columns=["gage_id","col","pct"])
                wide=long.pivot_table(index="gage_id",columns="col",values="pct",fill_value=0.0).reset_index()
                out=out.merge(wide,on="gage_id",how="left")
            if bi==0:
                out=out.merge(res[["gage_id","count"]].rename(columns={"count":"n_modis_pixels"}),on="gage_id",how="left")
        allcols=[COL[(bb,c)] for bb in BANDS for c in BAND_CODES[bb]]
        for col in allcols:
            if col not in out.columns: out[col]=0.0
        out[allcols]=out[allcols].fillna(0.0)
        out.to_parquet(outp)
        for f in glob.glob(os.path.join(tmpdir,"*")): os.remove(f)
        os.rmdir(tmpdir)
        return outp,f"{ntile} tiles, {len(out)} basins"
    except Exception as e:
        return outp,f"ERROR {year}: {type(e).__name__} {str(e)[:90]}"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--years",required=True); ap.add_argument("--workers",type=int,default=2)
    ap.add_argument("--geometry",default=f"{ROOT}/geometry_7964.gpkg")
    ap.add_argument("--outdir",default=f"{ROOT}/lulc_out")
    a=ap.parse_args(); os.makedirs(a.outdir,exist_ok=True)
    y0,y1=(int(x) for x in a.years.split("-")) if "-" in a.years else (int(a.years),int(a.years))
    years=list(range(y0,y1+1))
    g=gpd.read_file(a.geometry).to_crs(SINU)[["gage_id","canon_id","geometry"]]
    basins_path=os.path.join(a.outdir,"_basins_sinu.parquet"); g.to_parquet(basins_path)
    tiles=basin_tiles(g)
    print(f"LULC | basins {len(g)} | tiles {len(tiles)} | years {years[0]}-{years[-1]} | workers {a.workers} | class cols {len(LEG)}",flush=True)
    cache_path=os.path.join(a.outdir,"_granule_cache.parquet")
    t=time.time(); build_granule_cache(tiles,cache_path); print(f"granule cache: {time.time()-t:.0f}s",flush=True)
    tasks=[(y,cache_path,basins_path,a.outdir) for y in years]
    if a.workers>1:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for f in as_completed([ex.submit(process_year,x) for x in tasks]):
                p,msg=f.result(); print(f"  {os.path.basename(p)}: {msg}",flush=True)
    else:
        for x in tasks:
            p,msg=process_year(x); print(f"  {os.path.basename(p)}: {msg}",flush=True)
    parts=[pd.read_parquet(p) for p in sorted(glob.glob(os.path.join(a.outdir,"lulc_*.parquet"))) if pd.read_parquet(p).shape[1]>1]
    if parts:
        out=pd.concat(parts,ignore_index=True).sort_values(["gage_id","year"])
        out.to_parquet(os.path.join(a.outdir,"watershed_modis_lulc_annual.parquet"))
        print(f"FINAL: {len(out)} rows, {out.gage_id.nunique()} gages, {out.year.nunique()} years",flush=True)

if __name__=="__main__": main()

import os, glob, re
os.environ.setdefault("PROJ_LIB", "/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA", "/home/sagemaker-user/.conda/envs/geo/share/proj")
import pandas as pd, geopandas as gpd, pyogrio
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__))
def canon(s):
    s=str(s).strip()
    return (s.lstrip("0") or "0") if s.isdigit() else s

# --- passing gages with lat/lon/basin_area ---
m=pd.read_csv(os.path.join(HERE,"combined_watershed_metadata_09feb2026.csv"), dtype=str, low_memory=False)
m=m[m["processing_status"]=="success"].copy()
m["basin_area"]=pd.to_numeric(m["basin_area"], errors="coerce")
ids=m["gage_id"].astype(str)
is_can=ids.str.match(r"^\d{2}[A-Za-z]")
us_meta=m[~is_can].copy()
us_canon={canon(g):g for g in us_meta["gage_id"]}
print(f"US passing gages: {len(us_meta)}")

# --- read all GAGES-II boundary shapefiles ---
shps=sorted(glob.glob(os.path.join(HERE,"official","gagesII_bnd","boundaries-shapefiles-by-aggeco","*.shp")))
parts=[]
for f in shps:
    g=pyogrio.read_dataframe(f)
    g["aggeco"]=os.path.basename(f).replace("bas_","").replace(".shp","")
    parts.append(g[["GAGE_ID","geometry","aggeco"]])
g2=gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=parts[0].crs)
print(f"GAGES-II polygons: {len(g2)} | CRS: {g2.crs.to_string()[:40]}")
# area (km2) in native Albers (equal-area) BEFORE reproject
g2["geom_area_km2"]=g2.geometry.area/1e6
g2["gid_canon"]=g2["GAGE_ID"].map(canon)

# subset to US passing
g2u=g2[g2["gid_canon"].isin(us_canon)].copy()
# dedup (a gage could appear once per shapefile; keep largest)
g2u=g2u.sort_values("geom_area_km2", ascending=False).drop_duplicates("gid_canon", keep="first")
g2u["gage_id"]=g2u["gid_canon"].map(us_canon)
print(f"matched US gages: {g2u['gage_id'].nunique()} / {len(us_meta)}")

# reproject to EPSG:4326
g2u_4326=g2u.to_crs(4326)
g2u_4326["watershed_geom_source"]="gagesii"

# --- area QA vs reported basin_area ---
qa=g2u.merge(us_meta[["gage_id","basin_area"]], on="gage_id", how="left")
qa["area_rel_diff"]=(qa["geom_area_km2"]-qa["basin_area"]).abs()/qa["basin_area"]
good=qa["area_rel_diff"].replace([np.inf,-np.inf],np.nan).dropna()
print(f"\nArea QA (geom vs reported basin_area), n={len(good)}:")
print(f"  median rel.diff: {good.median():.3f} | <10%: {(good<0.10).mean()*100:.1f}% | <25%: {(good<0.25).mean()*100:.1f}% | >50%: {(good>0.50).sum()} gages")

# geometry validity
inv=(~g2u_4326.geometry.is_valid).sum()
print(f"invalid geometries: {inv} (will st_make_valid if >0)")

# save
out=os.path.join(HERE,"official","us_gagesII_4326.gpkg")
cols=["gage_id","watershed_geom_source","geom_area_km2","aggeco","geometry"]
g2u_4326[cols].to_file(out, driver="GPKG", layer="watersheds")
print(f"\nwrote {out}: {len(g2u_4326)} US watersheds (EPSG:4326)")
# list any US passing gage NOT matched
miss=set(us_meta["gage_id"]) - set(g2u_4326["gage_id"])
print(f"US passing NOT matched by GAGES-II: {len(miss)}", list(miss)[:10])

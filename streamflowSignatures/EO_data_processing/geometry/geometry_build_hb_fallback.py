import os
os.environ.setdefault("PROJ_LIB", "/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA", "/home/sagemaker-user/.conda/envs/geo/share/proj")
import pandas as pd, geopandas as gpd, numpy as np
from collections import defaultdict
from shapely import union_all

HERE=os.path.dirname(os.path.abspath(__file__))
BASIN=os.path.join(HERE,"basinAt_NorAm_polys.gpkg")

resid=pd.read_csv(os.path.join(HERE,"residual_for_hb.csv"), dtype=str)
resid["latitude"]=pd.to_numeric(resid["latitude"],errors="coerce")
resid["longitude"]=pd.to_numeric(resid["longitude"],errors="coerce")
print(f"residual gages: {len(resid)}")

print("loading basinAt (167k lev12 basins)...")
b=gpd.read_file(BASIN)
print("  basinAt CRS:", b.crs.to_string(), "| features:", len(b))
if b.crs is None or b.crs.to_epsg()!=4326:
    b=b.to_crs(4326)
b["HYBAS_ID"]=b["HYBAS_ID"].astype("int64")
b["NEXT_DOWN"]=b["NEXT_DOWN"].astype("int64")
hyb=b["HYBAS_ID"].to_numpy(); nxt=b["NEXT_DOWN"].to_numpy()
row_of={h:i for i,h in enumerate(hyb)}
child_of=defaultdict(list)
for h,n in zip(hyb,nxt):
    if n>0: child_of[n].append(h)
print("  adjacency built")

def bfs(outlet):
    if outlet not in row_of: return None
    seen=set(); st=[outlet]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        st.extend(child_of.get(x,()))
    return seen

# derive outlet for the 15 without Downstream_HB_ID via point-in-basin
need_pt=resid[resid["Downstream_HB_ID"].fillna("").str.match(r"^\d+$")==False].copy()
outlet_map={}  # gage_id -> outlet HYBAS_ID
for _,r in resid.iterrows():
    hb=str(r["Downstream_HB_ID"])
    if hb.isdigit(): outlet_map[r["gage_id"]]=int(hb)
if len(need_pt):
    pts=gpd.GeoDataFrame(need_pt, geometry=gpd.points_from_xy(need_pt.longitude, need_pt.latitude), crs=4326)
    j=gpd.sjoin(pts, b[["HYBAS_ID","geometry"]], predicate="within", how="left")
    for _,r in j.iterrows():
        if pd.notna(r.get("HYBAS_ID")): outlet_map[r["gage_id"]]=int(r["HYBAS_ID"])
print(f"outlets resolved: {len(outlet_map)}/{len(resid)} (pt-in-basin for {len(need_pt)})")

rows=[]; failed=[]
for _,r in resid.iterrows():
    g=r["gage_id"]; o=outlet_map.get(g)
    if o is None: failed.append((g,"no outlet")); continue
    mem=bfs(o)
    if not mem: failed.append((g,"outlet not in basinAt")); continue
    idx=[row_of[h] for h in mem if h in row_of]
    geom=union_all(b.geometry.values[idx])
    rows.append({"gage_id":g,"watershed_geom_source":"hydrobasins",
                 "geom_area_km2":np.nan,"n_members":len(idx),"geometry":geom})

hb=gpd.GeoDataFrame(rows, crs=4326)
# equal-area for area
hb["geom_area_km2"]=hb.to_crs("ESRI:102008").geometry.area/1e6
print(f"\nHB fallback delineated: {len(hb)} | failed: {len(failed)} {failed}")
print("invalid geoms:", int((~hb.geometry.is_valid).sum()))
hb.to_file(os.path.join(HERE,"official","hb_fallback_4326.gpkg"), driver="GPKG", layer="watersheds")
print("wrote hb_fallback_4326.gpkg")

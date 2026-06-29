import os, re
os.environ.setdefault("PROJ_LIB", "/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA", "/home/sagemaker-user/.conda/envs/geo/share/proj")
import pandas as pd, geopandas as gpd, pyogrio, numpy as np

HERE=os.path.dirname(os.path.abspath(__file__))
LAYER="DrainageBasin_BassinDeDrainage"

m=pd.read_csv(os.path.join(HERE,"combined_watershed_metadata_09feb2026.csv"), dtype=str, low_memory=False)
m=m[m["processing_status"]=="success"].copy()
m["basin_area"]=pd.to_numeric(m["basin_area"], errors="coerce")
ids=m["gage_id"].astype(str)
ca_meta=m[ids.str.match(r"^\d{2}[A-Za-z]")].copy()
print(f"CAN passing gages: {len(ca_meta)}")
ca_set=set(ca_meta["gage_id"])

parts=[]
for nn in [f"{i:02d}" for i in range(1,12)]:
    f=os.path.join(HERE,"official","ca",f"MDA_ADP_{nn}.gpkg")
    want=sorted(g for g in ca_set if g[:2]==nn)
    if not want: continue
    where="StationNum IN (" + ",".join("'%s'"%s for s in want) + ")"
    try:
        g=pyogrio.read_dataframe(f, layer=LAYER, where=where)
    except Exception as e:
        print(f"  {nn}: where failed ({e}); reading full + filter")
        g=pyogrio.read_dataframe(f, layer=LAYER)
        g=g[g["StationNum"].astype(str).isin(want)]
    g["mda_file"]=nn
    parts.append(g)
    print(f"  MDA_{nn}: wanted {len(want)}, got {len(g)}")

ca=gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=parts[0].crs)
print(f"\nCanada polygons read: {len(ca)} | src CRS: {ca.crs.to_string()[:50]}")
ca["geom_area_km2"]=ca.geometry.area/1e6   # src CRS is Canada Albers (equal-area)
ca["gage_id"]=ca["StationNum"].astype(str).str.strip()
ca=ca.sort_values("geom_area_km2", ascending=False).drop_duplicates("gage_id", keep="first")

ca_4326=ca.to_crs(4326)
ca_4326["watershed_geom_source"]="wsc_eccc"
print(f"matched CAN gages: {ca_4326['gage_id'].nunique()} / {len(ca_meta)}")

# area QA
qa=ca.merge(ca_meta[["gage_id","basin_area"]], on="gage_id", how="left")
qa["rel"]=(qa["geom_area_km2"]-qa["basin_area"]).abs()/qa["basin_area"]
good=qa["rel"].replace([np.inf,-np.inf],np.nan).dropna()
print(f"Area QA vs basin_area, n={len(good)}: median {good.median():.3f} | <10% {(good<0.10).mean()*100:.1f}% | <25% {(good<0.25).mean()*100:.1f}% | >50% {(good>0.5).sum()}")
inv=(~ca_4326.geometry.is_valid).sum(); print("invalid geometries:", inv)

out=os.path.join(HERE,"official","ca_mda_4326.gpkg")
ca_4326[["gage_id","watershed_geom_source","geom_area_km2","mda_file","geometry"]].to_file(out, driver="GPKG", layer="watersheds")
print(f"wrote {out}: {len(ca_4326)} CAN watersheds (EPSG:4326)")
miss=ca_set - set(ca_4326["gage_id"])
print(f"CAN passing NOT matched: {len(miss)} ->", sorted(miss)[:35])

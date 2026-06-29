import os, glob
os.environ.setdefault("PROJ_LIB", "/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA", "/home/sagemaker-user/.conda/envs/geo/share/proj")
import pandas as pd, geopandas as gpd, numpy as np
from shapely import make_valid

HERE=os.path.dirname(os.path.abspath(__file__))
EQ="EPSG:6933"   # global cylindrical equal-area -> correct areas worldwide
TOL=0.002        # ~200 m
MAXREL=0.02      # if simplify changes area >2%, keep full-res (protects small/narrow basins)

def canon(x):
    s=str(x).strip(); return (str(int(s)) if s.isdigit() else s.upper())

# authoritative zero-padded ids from signature files
auth=set()
for p in ["/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures/golden-outputs/streamflow_signatures_full_10feb2026.csv",
          os.path.join(HERE,"official","sig_jan2026.csv")]:
    auth|=set(pd.read_csv(p,usecols=[0],dtype=str).iloc[:,0].astype(str).str.strip())
amap={}
for a in auth:
    c=canon(a); amap.setdefault(c,set()).add(a)
assert all(len(v)==1 for v in amap.values()), "canon collision in authoritative signature ids!"
amap={k:next(iter(v)) for k,v in amap.items()}
def padded(c): return amap.get(c, c.zfill(8) if c.isdigit() else c)

# full-res per-source layers
def load(fn,src):
    g=gpd.read_file(os.path.join(HERE,"official",fn))[["gage_id","geometry"]].copy()
    g["watershed_geom_source"]=src; return g
us=load("us_gagesII_4326.gpkg","gagesii")
ca=load("ca_mda_4326.gpkg","wsc_eccc")
hb=gpd.read_file(os.path.join(HERE,"official","hb_fallback_4326.gpkg"))[["gage_id","geometry"]].copy()
hb["watershed_geom_source"]="hydrobasins"
# the 4 inclusive-universe gages from GAGES-II (full-res)
NEW={"1591000","1591400","1591610","1591700"}; parts=[]
for f in glob.glob(os.path.join(HERE,"official","gagesII_bnd","boundaries-shapefiles-by-aggeco","*.shp")):
    gg=gpd.read_file(f); gg=gg[gg["GAGE_ID"].map(canon).isin(NEW)]
    if len(gg): parts.append(gg[["GAGE_ID","geometry"]].rename(columns={"GAGE_ID":"gage_id"}))
new=gpd.GeoDataFrame(pd.concat(parts,ignore_index=True),crs=parts[0].crs).to_crs(4326)
new["watershed_geom_source"]="gagesii"

final=gpd.GeoDataFrame(pd.concat([us,new,ca,hb],ignore_index=True),crs="EPSG:4326")
final["canon_id"]=final["gage_id"].map(canon)
final["gage_id"]=final["canon_id"].map(padded)
assert final["canon_id"].is_unique, "duplicate canon_id across sources!"
assert final["gage_id"].is_unique, "duplicate padded gage_id!"

# metadata attrs by canon
m=pd.read_csv(os.path.join(HERE,"combined_watershed_metadata_09feb2026.csv"),dtype=str,low_memory=False)
for c in ("basin_area","latitude","longitude"): m[c]=pd.to_numeric(m[c],errors="coerce")
m["canon"]=m["gage_id"].map(canon)
final=final.merge(m[["canon","basin_area","latitude","longitude","gage_type"]],left_on="canon_id",right_on="canon",how="left").drop(columns=["canon"])

# adaptive simplify: keep full-res where simplify would distort area > MAXREL
full=final.geometry
simp=full.simplify(TOL,preserve_topology=True)
simp=gpd.GeoSeries([make_valid(s) if (s is not None and not s.is_empty and not s.is_valid) else s for s in simp],crs=4326)
af=full.to_crs(EQ).area/1e6; asi=simp.to_crs(EQ).area/1e6
rel=((asi-af).abs()/af).replace([np.inf,-np.inf],np.nan)
use_full=(rel>MAXREL)|asi.isna()|simp.is_empty
final["geometry"]=gpd.GeoSeries([f if uf else s for uf,f,s in zip(use_full,full,simp)],crs=4326)
final["geom_simplified"]=~use_full.values
# area FROM DELIVERED geometry
final["geom_area_km2"]=final.geometry.to_crs(EQ).area/1e6
final["area_rel_diff"]=((final["geom_area_km2"]-final["basin_area"]).abs()/final["basin_area"]).replace([np.inf,-np.inf],np.nan)
final["area_flag"]=final["area_rel_diff"]>0.5
final["low_confidence"]=(final["watershed_geom_source"]=="hydrobasins")|final["area_flag"].fillna(False)

print("=== QA v3 ===")
print("rows:",len(final),"| unique gage_id:",final["gage_id"].nunique(),"| by source:",final["watershed_geom_source"].value_counts().to_dict())
print("dups:",int(final["gage_id"].duplicated().sum()),"| invalid:",int((~final.geometry.is_valid).sum()),"| empty:",int(final.geometry.is_empty.sum()))
print("kept full-res (simplify>2% area):",int((~final['geom_simplified']).sum()),"| simplified:",int(final['geom_simplified'].sum()))
g=final["area_rel_diff"].dropna()
print(f"area(delivered) vs basin_area n={len(g)}: median {g.median():.4f} | <10% {(g<0.1).mean()*100:.1f}% | <25% {(g<0.25).mean()*100:.1f}% | >50% {int((g>0.5).sum())}")
print("low_confidence:",int(final['low_confidence'].sum()),"(hydrobasins",int((final['watershed_geom_source']=='hydrobasins').sum()),"+ area_flag)")
print("US gage_id min len:",final[~final.canon_id.str.match(r'^\d{2}[A-Za-z]')]['gage_id'].str.len().min())

out=os.path.join(HERE,"watershed_polygons_26jun2026.gpkg")
if os.path.exists(out): os.remove(out)
cols=["gage_id","canon_id","watershed_geom_source","geom_area_km2","basin_area","area_rel_diff","area_flag","geom_simplified","low_confidence","latitude","longitude","gage_type","geometry"]
final[cols].to_file(out,driver="GPKG",layer="watersheds")
final[cols].to_parquet(os.path.join(HERE,"watershed_polygons_26jun2026.parquet"))
final[cols].drop(columns=["geometry"]).to_csv(os.path.join(HERE,"watershed_polygons_26jun2026_qa.csv"),index=False)
print(f"wrote gpkg {os.path.getsize(out)/1e6:.1f} MB + parquet + qa.csv")

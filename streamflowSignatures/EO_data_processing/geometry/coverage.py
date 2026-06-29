import os, re
os.environ.setdefault("PROJ_LIB", "/home/sagemaker-user/.conda/envs/geo/share/proj")
import pandas as pd, pyogrio

HERE=os.path.dirname(__file__)
def canon(s):
    s=str(s).strip()
    return s.lstrip("0") or "0" if s.isdigit() else s

# passing gages
m=pd.read_csv(os.path.join(HERE,"metadata_slim.csv"), dtype=str)
m=m[m["processing_status"]=="success"].copy()
ids=m["gage_id"].astype(str)
is_can=ids.str.match(r"^\d{2}[A-Za-z]")
us=set(ids[~is_can]); ca=set(ids[is_can])
hb_numeric=set(ids[m["Downstream_HB_ID"].fillna("").str.match(r"^\d+$")])
print(f"passing gages: {len(ids)}  | US-like: {len(us)}  | CAN-like: {len(ca)}")

# GAGES-II ids
g=pyogrio.read_dataframe(os.path.join(HERE,"official","gagesII_point.zip"), read_geometry=False)
idcol=[c for c in g.columns if c.upper() in ("STAID","GAGE_ID","GAGEID")][0]
gages2=set(g[idcol].astype(str))
print(f"\nGAGES-II stations: {len(gages2)} (key col '{idcol}')")

# Canada MDA_ADP included stations
ca_tab=pd.read_csv(os.path.join(HERE,"official","ca_included_stations.txt"), sep="\t", dtype=str)
ca_inc=set(ca_tab.iloc[:,0].astype(str).str.strip())
print(f"Canada MDA_ADP included stations: {len(ca_inc)}")

def cov(passing, source, label):
    pc={canon(x):x for x in passing}; sc={canon(x) for x in source}
    hit=[pc[k] for k in pc if k in sc]
    miss=[pc[k] for k in pc if k not in sc]
    print(f"  {label}: {len(hit)}/{len(passing)} covered ({100*len(hit)/max(1,len(passing)):.1f}%)  | miss {len(miss)}")
    return set(hit), set(miss)

print("\n=== US coverage ===")
us_g2_hit, us_g2_miss = cov(us, gages2, "GAGES-II")
# of the GAGES-II misses, how many would NLDI plausibly cover (CONUS) vs HB fallback
us_miss_hb = us_g2_miss & hb_numeric
print(f"    of {len(us_g2_miss)} GAGES-II misses -> NLDI-fill candidates (also have HB fallback: {len(us_miss_hb)})")

print("\n=== Canada coverage ===")
ca_hit, ca_miss = cov(ca, ca_inc, "ECCC MDA_ADP")
ca_miss_hb = ca_miss & hb_numeric
print(f"    of {len(ca_miss)} MDA_ADP misses -> have HB fallback: {len(ca_miss_hb)}, no geometry at all: {len(ca_miss - hb_numeric)}")

print("\n=== summary (official primary, HB fallback) ===")
official_hit = us_g2_hit | ca_hit
print(f"covered by official (GAGES-II + MDA_ADP) directly: {len(official_hit)}/{len(ids)} ({100*len(official_hit)/len(ids):.1f}%)")
print(f"US needing NLDI fill: {len(us_g2_miss)}")
remaining = (us_g2_miss | ca_miss)
no_any = remaining - hb_numeric
print(f"still uncovered after official: {len(remaining)} (US {len(us_g2_miss)}, CAN {len(ca_miss)})")
print(f"  of those, NO HB fallback either (need NLDI/StreamStats or drop): {len(no_any)}")
print("  sample no-geometry-at-all:", list(no_any)[:15])

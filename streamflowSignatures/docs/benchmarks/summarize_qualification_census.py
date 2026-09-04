"""Summarize a qualification_census.jl run against the two standard products.

Usage: python docs/benchmarks/summarize_qualification_census.py <census.csv> <product1_signatures.csv> <product2_signatures.csv> <combined_watershed_metadata.csv>
Cross-checks the census include/exclude decisions against each product's gage set (must be
0/0 mismatches), then reports exclusion reasons and rejected-year counts by reason.
Written 2026-09-04 (docs/plans/2026-09-04-filtering-stage-census.md).
"""
import pandas as pd, numpy as np, json
import sys; CENSUS=sys.argv[1]; P1=sys.argv[2]; P2=sys.argv[3]; META=sys.argv[4]
c=pd.read_csv(CENSUS,dtype={"gage_id":str})
y=pd.read_csv(CENSUS.replace(".csv","_years.csv"),dtype={"gage_id":str})
meta=pd.read_csv(META,dtype=str,low_memory=False)
canon=lambda x: x.lstrip("0") if x.isdigit() else x.upper()
gt=dict(zip(meta.gage_id.map(canon),meta.gage_type)); c["gage_type"]=c.gage_id.map(canon).map(gt)
p1=set(pd.read_csv(P1,usecols=["gage_id"],dtype=str).gage_id)
p2=set(pd.read_csv(P2,usecols=["gage_id"],dtype=str).gage_id)
print("gages in census:",c.gage_id.nunique())
print("\n=== WHOLE-RECORD gage-year census (WY1980-2026 rows) ===")
print(y.status.value_counts().to_string()); print("gage-years total:",len(y))
y26=y[y.water_year==2026]; print("WY2026 rows:",len(y26), y26.status.value_counts().to_dict())
yy=y[y.water_year<=2025]; print("WY1980-2025 gage-years:",len(yy)," valid:",(yy.status=='valid').sum()," rejected:",(yy.status!='valid').sum())
print(yy.status.value_counts().to_string())
for w,prod,lo in [("1993-2025",p1,1993),("1980-2025",p2,1980)]:
    d=c[c.window==w].copy(); inc=set(d[d.included].gage_id)
    print(f"\n=== WINDOW {w} ===")
    print(f"census included: {len(inc)}  product: {len(prod)}  census-only: {len(inc-prod)}  product-only: {len(prod-inc)}")
    if inc^prod: print("   mismatch ids:", sorted(inc^prod)[:20])
    exc=d[~d.included]
    print("excluded:",len(exc))
    r=[]
    for _,row in exc.iterrows():
        if row.n_wy_with_rows==0: r.append("no rows in window")
        elif row.fails_frac and row.fails_min_years: r.append("both: <60% of window AND <20 valid years")
        elif row.fails_frac: r.append("<60% of window only (has ≥20 valid years)")
        elif row.fails_min_years: r.append("<20 valid years only (≥60% of window)")
        else: r.append("?")
    exc=exc.assign(reason=r); print(exc.reason.value_counts().to_string())
    print("excluded by gage_type:", exc.gage_type.value_counts().to_dict(), " included by type:", d[d.included].gage_type.value_counts().to_dict())
    # of excluded, how many had zero valid years in window
    print("excluded with 0 valid years in window:",(exc.n_valid==0).sum(), " excluded with ≥20 valid years:",(exc.n_valid>=20).sum())
    # year-level accounting for INCLUDED gages
    di=d[d.included]
    print(f"included gages: WYs with rows in window = {int(di.n_wy_with_rows.sum())}, valid = {int(di.n_valid.sum())}, rejected = {int(di.n_rejected.sum())} "
          f"(no_data {int(di.rej_no_data.sum())}, too_many_na {int(di.rej_too_many_na.sum())}, gap {int(di.rej_gap.sum())}, negative {int(di.rej_negative.sum())}, residual {int(di.rej_residual.sum())}, unlisted {int(di.rej_unlisted.sum())})")
    print(f"ALL 8,014 gages: WYs with rows in window = {int(d.n_wy_with_rows.sum())}, valid = {int(d.n_valid.sum())}, rejected = {int(d.n_rejected.sum())} "
          f"(no_data {int(d.rej_no_data.sum())}, too_many_na {int(d.rej_too_many_na.sum())}, gap {int(d.rej_gap.sum())}, negative {int(d.rej_negative.sum())}, residual {int(d.rej_residual.sum())}, unlisted {int(d.rej_unlisted.sum())})")
    print(f"valid gage-years lost with excluded gages: {int(exc.n_valid.sum())}")
    print("included n_valid: min",di.n_valid.min()," max",di.n_valid.max()," median",di.n_valid.median(), " gages with 0 rejected years:",int((di.n_rejected==0).sum()))
    print("total_possible for included: min",di.total_possible.min()," max",di.total_possible.max(), " frac min", round(di.frac.min(),4))
d=c[c.window=="full"]; print("\n=== FULL RECORD (no window) ===")
print("valid years per gage: min",d.n_valid.min()," max",d.n_valid.max()," gages with ≥20 valid:",(d.n_valid>=20).sum(), " total valid gage-years:",int(d.n_valid.sum()), " rejected:",int(d.n_rejected.sum()))
c.to_csv(CENSUS.replace(".csv","_typed.csv"),index=False)

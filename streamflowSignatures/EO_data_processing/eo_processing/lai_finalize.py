"""Finalize the monthly MODIS LAI product: merge the catalog-derived table with the LP DAAC
backfills, add the partial_month flag, regenerate the always-NA sidecar, and write the data
dictionary. Idempotent; run after lai_pipeline.py + both lai_backfill_lpdaac.py targets.

Inputs (under data_out/eo_lai/):
  out_qa/watershed_modis_lai_monthly.parquet   (catalog-only, from lai_pipeline.py)   [or _prebackfill.parquet]
  backfill_farN/backfill_lai_monthly.parquet   (92 far-N basins, all months)
  backfill_2024_11/lai_2024_11.parquet         (all basins, Nov-2024)
Outputs (out_qa/): watershed_modis_lai_monthly.{parquet,csv} + _{date}.parquet,
  lai_always_na_basins.csv, watershed_modis_lai_dictionary.csv

Usage: python lai_finalize.py [--date 29jun2026]
"""
import os, glob, argparse, shutil
import pandas as pd, geopandas as gpd
ROOT="/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures/data_out/eo_lai"
OUT=f"{ROOT}/out_qa"

DICT=[("gage_id","string","zero-padded gage id; joins to streamflow_signatures"),
("canon_id","string","zero-stripped canonical id; joins to combined_watershed_metadata"),
("year","int","calendar year"),("month","int","calendar month 1-12"),
("lai_mean","m2/m2","watershed coverage-weighted mean of the per-pixel monthly-mean LAI"),
("lai_std","m2/m2","spatial std of monthly-mean LAI across watershed pixels (heterogeneity)"),
("lai_min","m2/m2","spatial min (UNWEIGHTED extremum; prefer lai_q05 for robust low)"),
("lai_q05","m2/m2","coverage-weighted spatial 5th pctile"),("lai_q25","m2/m2","coverage-weighted spatial 25th pctile"),
("lai_q50","m2/m2","coverage-weighted spatial median"),("lai_q75","m2/m2","coverage-weighted spatial 75th pctile"),
("lai_q95","m2/m2","coverage-weighted spatial 95th pctile"),("lai_max","m2/m2","spatial max (UNWEIGHTED extremum; prefer lai_q95)"),
("n_eff_pixels","count","coverage-weighted # pixels with usable monthly LAI"),
("n_basin_pixels","count","coverage-weighted # land pixels in basin (valid_pixel_frac denominator)"),
("valid_pixel_frac","0-1","n_eff_pixels / n_basin_pixels (spatial % of basin with usable LAI; saturates ~1)"),
("mean_obs_valid_frac","0-1","fraction of land observations (pixel x composite) that were clear/usable"),
("frac_cloud","0-1","fraction of land obs flagged cloud/cirrus/shadow (FparExtra_QC)"),
("frac_snow_ice","0-1","fraction of land obs flagged snow/ice"),
("frac_aerosol","0-1","fraction of land obs flagged high aerosol"),
("frac_backup_algorithm","0-1","fraction of land obs from empirical/backup retrieval (SCF_QC 2-3) vs main RT"),
("n_composites","count","# QA-good 4-day MCD15A3H composites contributing to the month (month-level)"),
("low_confidence","bool","geometry low-confidence (HydroBasins fallback / area outlier)"),
("partial_month","bool","True if n_composites<6 (thin month: fewer obs, less robust; incl. late-2024 record edge)"),
("good_coverage_frac","0-1","valid pixel-observations / (all basin pixels x ALL expected composite dates in the month) "
 "= mean_obs_valid_frac x min(1, n_composites/expected_for_calendar_month). 0 = no usable data "
 "(urban/deep winter/cloud/missing), 1 = fully observed & clear. Continuous generalization of partial_month.")]

def load_concat(d):
    f=f"{d}/backfill_lai_monthly.parquet"
    if os.path.exists(f): return pd.read_parquet(f)
    ps=[p for p in glob.glob(f"{d}/lai_*.parquet") if pd.read_parquet(p).shape[1]>1]
    return pd.concat([pd.read_parquet(p) for p in ps],ignore_index=True) if ps else pd.DataFrame()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--date",default="29jun2026"); a=ap.parse_args()
    main_f=f"{OUT}/watershed_modis_lai_monthly_prebackfill.parquet"
    if not os.path.exists(main_f): main_f=f"{OUT}/watershed_modis_lai_monthly.parquet"
    m=pd.read_parquet(main_f); m["gage_id"]=m.gage_id.astype(str)
    bf=pd.concat([load_concat(f"{ROOT}/backfill_farN"),load_concat(f"{ROOT}/backfill_2024_11")],ignore_index=True)
    bf["gage_id"]=bf.gage_id.astype(str)
    bf=bf.drop_duplicates(["gage_id","year","month"],keep="first")
    k=lambda d:d.gage_id+"_"+d.year.astype(str)+"_"+d.month.astype(str)
    merged=pd.concat([m[~k(m).isin(set(k(bf)))],bf],ignore_index=True).sort_values(["gage_id","year","month"])
    merged["partial_month"]=merged.n_composites<6
    # % good coverage = valid pixel-obs / (all basin pixels x ALL expected composite dates in the month)
    #   = mean_obs_valid_frac (areal x among-present-dates) x (present composites / expected composites)
    # expected = normal full count for that calendar month (mode; robust to the few gap months)
    nexp={mm:int(merged.loc[merged.month==mm,"n_composites"].mode().iloc[0]) for mm in range(1,13)}
    temporal=(merged.n_composites/merged.month.map(nexp)).clip(0,1)
    merged["good_coverage_frac"]=(merged.mean_obs_valid_frac.fillna(0)*temporal).round(4)
    assert merged.duplicated(["gage_id","year","month"]).sum()==0, "duplicate keys after merge"
    # write product
    if os.path.exists(f"{OUT}/watershed_modis_lai_monthly.parquet") and not os.path.exists(f"{OUT}/watershed_modis_lai_monthly_prebackfill.parquet"):
        shutil.copy(f"{OUT}/watershed_modis_lai_monthly.parquet", f"{OUT}/watershed_modis_lai_monthly_prebackfill.parquet")
    merged.to_parquet(f"{OUT}/watershed_modis_lai_monthly.parquet")
    merged.to_parquet(f"{OUT}/watershed_modis_lai_monthly_{a.date}.parquet")
    merged.to_csv(f"{OUT}/watershed_modis_lai_monthly.csv",index=False)
    # regenerate always-NA sidecar (basins with 0 valid LAI months = the urban basins)
    valid=merged.groupby("gage_id").lai_mean.apply(lambda s:int(s.notna().sum()))
    na_ids=set(valid[valid==0].index)
    geo=gpd.read_file(f"{ROOT}/geometry_7964.gpkg")[["gage_id","canon_id","latitude","longitude","watershed_geom_source"]]
    geo["gage_id"]=geo.gage_id.astype(str)
    sc=geo[geo.gage_id.isin(na_ids)].copy()
    sc["na_reason"]="urban_builtup_no_modis_lai"; sc["backfillable_from_lpdaac"]=False
    sc=sc.sort_values("gage_id"); sc.to_csv(f"{OUT}/lai_always_na_basins.csv",index=False)
    sc.to_csv(f"{OUT}/watershed_modis_lai_na_basins_{a.date}.csv",index=False)
    dd=pd.DataFrame(DICT,columns=["column","unit","description"])
    dd.to_csv(f"{OUT}/watershed_modis_lai_dictionary.csv",index=False)
    dd.to_csv(f"{OUT}/watershed_modis_lai_dictionary_{a.date}.csv",index=False)
    print(f"finalized: {len(merged)} rows, {merged.gage_id.nunique()} gages, "
          f"{merged[['year','month']].drop_duplicates().shape[0]} months, partial={int(merged.partial_month.sum())}, always-NA={len(na_ids)}")

if __name__=="__main__": main()

"""Finalize the annual MODIS LULC product: take the assembled per-year table
(watershed_modis_lulc_annual.parquet from lulc_pipeline.py), attach geometry provenance
+ a low_confidence flag (geometry low-confidence OR sub-5-pixel MODIS support), write the
date-stamped parquet + CSV, and generate the data dictionary. Idempotent.

Inputs (under data_out/eo_lai/lulc_out/):
  watershed_modis_lulc_annual.parquet   (one row per gage_id x year, from lulc_pipeline.py)
Geometry: data_out/eo_lai/geometry_7964.gpkg  (watershed_geom_source, low_confidence)
Legends:  EO_data_processing/eo_processing/lulc_legends.csv

Outputs (lulc_out/): watershed_modis_lulc_annual_{date}.{parquet,csv},
  watershed_modis_lulc_dictionary.csv

Usage: python lulc_finalize.py [--date 29jun2026]
"""
import os, argparse
import pandas as pd, geopandas as gpd
ROOT="/home/sagemaker-user/SurfaceWaterProjections/streamflowSignatures/data_out/eo_lai"
OUT=f"{ROOT}/lulc_out"
HERE=os.path.dirname(os.path.abspath(__file__))
LOWPX=5.0  # < this many effective MODIS pixels -> low-confidence support

BAND_DESC={
 "LC_Type1":"IGBP global vegetation classification (17 classes)",
 "LC_Type2":"University of Maryland (UMD) classification (16 classes)",
 "LC_Type3":"MODIS-derived LAI/fPAR classification (11 classes)",
 "LC_Type4":"BIOME-BGC classification (9 classes)",
 "LC_Type5":"Plant Functional Type (PFT) classification (12 classes)",
 "LC_Prop1":"FAO-LCCS1 land cover layer (sparse codes)",
 "LC_Prop2":"FAO-LCCS2 land use layer (sparse codes)",
 "LC_Prop3":"FAO-LCCS3 surface hydrology layer (sparse codes)"}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--date",default="29jun2026"); a=ap.parse_args()
    df=pd.read_parquet(f"{OUT}/watershed_modis_lulc_annual.parquet"); df["gage_id"]=df.gage_id.astype(str)
    df=df.drop(columns=[c for c in ["watershed_geom_source","low_confidence","low_pixel_support"] if c in df.columns])  # idempotent: strip prior finalize cols
    leg=pd.read_csv(os.path.join(HERE,"lulc_legends.csv"))

    # geometry provenance + low_confidence (read attribute table only; avoid PROJ db)
    g=gpd.read_file(f"{ROOT}/geometry_7964.gpkg",ignore_geometry=True)[["gage_id","watershed_geom_source","low_confidence"]]
    g["gage_id"]=g.gage_id.astype(str)
    df=df.merge(g,on="gage_id",how="left",validate="many_to_one")   # guard against geometry dup-key row blow-up
    assert df.watershed_geom_source.notna().all(), "null watershed_geom_source after geometry join (gage_id mismatch)"
    df["geom_low_confidence"]=df.low_confidence.fillna(False).astype(bool)
    df["low_pixel_support"]=df.n_modis_pixels<LOWPX
    df["low_confidence"]=df.geom_low_confidence | df.low_pixel_support
    df=df.drop(columns=["geom_low_confidence"])

    # column order: keys -> QA -> class columns (band, code order)
    classcols=[c for b in leg.band.unique() for c in leg[leg.band==b].output_column if c in df.columns]
    keys=["gage_id","canon_id","year"]
    qa=["watershed_geom_source","n_modis_pixels","low_pixel_support","low_confidence"]
    df=df[keys+qa+classcols].sort_values(["gage_id","year"]).reset_index(drop=True)

    df.to_parquet(f"{OUT}/watershed_modis_lulc_annual_{a.date}.parquet")
    df.to_parquet(f"{OUT}/watershed_modis_lulc_annual.parquet")
    df.to_csv(f"{OUT}/watershed_modis_lulc_annual_{a.date}.csv",index=False)

    # data dictionary
    rows=[("gage_id","string","zero-padded gage id; joins to streamflow_signatures"),
          ("canon_id","string","zero-stripped canonical id; joins to combined_watershed_metadata"),
          ("year","int","calendar year (MCD12Q1 annual product, Jan-1 stamped)"),
          ("watershed_geom_source","string","gagesii | wsc_eccc | hydrobasins (basin polygon provenance)"),
          ("n_modis_pixels","count","coverage-weighted # of valid (non-fill) 500m MODIS pixels in the basin"),
          ("low_pixel_support","bool",f"True if n_modis_pixels < {int(LOWPX)} (sub-pixel/tiny basin; class %s less robust)"),
          ("low_confidence","bool","True if geometry low-confidence (HydroBasins fallback/area outlier) OR low_pixel_support")]
    for r in leg.itertuples():
        nm=getattr(r,"class_name","") or ""
        rows.append((r.output_column,"% (0-100)",
                     f"{BAND_DESC.get(r.band,r.band)} - class {int(r.code)}"+(f' "{nm}"' if nm else "")+
                     "; coverage-weighted % of valid pixels"))
    pd.DataFrame(rows,columns=["column","unit","description"]).to_csv(f"{OUT}/watershed_modis_lulc_dictionary.csv",index=False)

    print(f"finalized: {len(df)} rows, {df.gage_id.nunique()} gages, {df.year.nunique()} years, "
          f"{len(classcols)} class cols | low_confidence basin-years={int(df.low_confidence.sum())} "
          f"(geom {int(g.low_confidence.fillna(False).sum())} basins, low_px {int(df.low_pixel_support.sum())} basin-years)")
    print(f"  -> watershed_modis_lulc_annual_{a.date}.{{parquet,csv}} + _dictionary.csv")

if __name__=="__main__": main()

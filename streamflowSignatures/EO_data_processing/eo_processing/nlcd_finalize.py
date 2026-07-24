"""Finalize the annual NLCD product (Annual NLCD Collection 1 C1V2, CONUS 1985-2025).

Takes the assembled per-year table (watershed_nlcd_annual.parquet from nlcd_pipeline.py), which
covers all 6,164 US gages incl. 45 out-of-footprint Alaska gages (all-fill -> zeros), and:
  - DROPS the out-of-footprint gages (identified conservatively by max(n_nlcd_pixels)==0 over all
    years, so their published rows never masquerade as real all-zero land cover) to a traceable
    exclusion table; the CONUS product keeps ~6,119 gages x 41 yr.
  - attaches geometry provenance (watershed_geom_source) + SEPARATE, documented QA flags
    (geom_low_confidence / low_pixel_support / partial_coverage) and their union low_confidence.
  - adds provenance (nlcd_collection, valid_area_km2) + a machine-readable metadata JSON.
  - writes date-stamped parquet + CSV + a data dictionary + an all-gage QA companion.
  - S3 upload is OPT-IN (--upload); default OFF so the product is reviewed via the explorer first.

Design + review reconciliation: EO_data_processing/README_NLCD.md. Mirrors lulc_finalize.py.

Usage:
  python nlcd_finalize.py --date 24jul2026            # local finalize only
  python nlcd_finalize.py --date 24jul2026 --upload   # + upload to S3 (after QA)
"""
import os, json, argparse, subprocess
import numpy as np, pandas as pd, geopandas as gpd

ROOT = "/home/sagemaker-user/streamflowSignatures/data_out/eo_nlcd"
OUT = f"{ROOT}/nlcd_out"
HERE = os.path.dirname(os.path.abspath(__file__))
GEOM = f"{ROOT}/watershed_polygons_26jun2026.parquet"
S3_DEST = "s3://climate-ai-data-science-shiny-app-data/streamflow"
COLLECTION = "C1V2"
YEARS = (1985, 2025)
LOWPX = 100.0        # < this many effective 30 m pixels (~0.09 km2) -> low pixel support (QA guard)
PARTIAL = 0.99       # valid_coverage_frac < this -> basin partly outside the raster footprint
# annual-update-path years for C1V2 (rule-based ICC update, not the full baseline) — methodological seam
ANNUAL_UPDATE_YEARS = "2024-2025"

CLASS_DESC = {  # NLCD CONUS legend (modified Anderson Level II)
 11:"Open Water",12:"Perennial Ice/Snow",21:"Developed, Open Space",22:"Developed, Low Intensity",
 23:"Developed, Medium Intensity",24:"Developed, High Intensity",31:"Barren Land",
 41:"Deciduous Forest",42:"Evergreen Forest",43:"Mixed Forest",52:"Shrub/Scrub",
 71:"Grassland/Herbaceous",81:"Pasture/Hay",82:"Cultivated Crops",90:"Woody Wetlands",
 95:"Emergent Herbaceous Wetlands"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--date", default="24jul2026")
    ap.add_argument("--upload", action="store_true", help="upload finalized artifacts to S3 (do AFTER QA)")
    a = ap.parse_args()
    leg = pd.read_csv(os.path.join(HERE, "nlcd_legends.csv"))
    classcols = [c for c in leg.output_column]  # committed order

    df = pd.read_parquet(f"{OUT}/watershed_nlcd_annual.parquet")
    df["gage_id"] = df.gage_id.astype(str)
    # idempotent: strip any prior finalize-added columns
    df = df.drop(columns=[c for c in ["watershed_geom_source", "nlcd_collection", "valid_area_km2",
                                      "geom_low_confidence", "low_pixel_support", "partial_coverage",
                                      "low_confidence"] if c in df.columns])

    # ---- geometry provenance + coordinates (attribute table only) ----
    g = gpd.read_parquet(GEOM)[["gage_id", "canon_id", "watershed_geom_source", "low_confidence",
                                "latitude", "longitude"]].copy()
    g["gage_id"] = g.gage_id.astype(str)
    assert g.gage_id.is_unique, "geometry gage_id not unique -> join would blow up rows"
    # canon_id integrity: annual table's canon_id is authoritative; confirm it agrees with geometry
    cid = g.set_index("gage_id").canon_id.astype(str)
    chk = df[["gage_id", "canon_id"]].drop_duplicates()
    mism = chk[chk.gage_id.map(cid).notna() & (chk.canon_id.astype(str) != chk.gage_id.map(cid))]
    assert len(mism) == 0, f"canon_id disagrees between annual table and geometry for {len(mism)} gages"

    # ---- outside-footprint (conservative: no valid pixels in ANY year) -> exclude + report ----
    maxpix = df.groupby("gage_id").n_nlcd_pixels.max()
    oof_ids = set(maxpix[maxpix == 0].index)
    excl = g[g.gage_id.isin(oof_ids)].copy()
    # classify jurisdiction from coordinates (QA confirmed all Alaska; be explicit, don't assume PR)
    def juris(r):
        return "Alaska" if r.latitude >= 51 else ("Hawaii" if (r.latitude < 24 and r.longitude < -140)
                                                  else ("Puerto Rico" if r.longitude > -68 else "non-CONUS"))
    excl["jurisdiction"] = excl.apply(juris, axis=1)
    excl["max_n_nlcd_pixels"] = excl.gage_id.map(maxpix).fillna(0)
    excl["reason"] = excl.jurisdiction + " — outside CONUS Annual NLCD footprint (all pixels fill=250)"
    excl[["gage_id", "canon_id", "latitude", "longitude", "watershed_geom_source",
          "jurisdiction", "max_n_nlcd_pixels", "reason"]].to_csv(f"{OUT}/nlcd_out_of_footprint_{a.date}.csv", index=False)
    print(f"excluded (out-of-footprint): {len(excl)} gages "
          f"[{dict(excl.jurisdiction.value_counts())}] -> nlcd_out_of_footprint_{a.date}.csv", flush=True)

    df = df[~df.gage_id.isin(oof_ids)].copy()  # CONUS product only

    # ---- join provenance/flags (validate one-to-one on the per-year rows) ----
    df = df.merge(g[["gage_id", "watershed_geom_source", "low_confidence"]],
                  on="gage_id", how="left", validate="many_to_one", indicator=True)
    assert (df._merge == "both").all(), "unmatched gage_id after geometry join"
    assert df.watershed_geom_source.notna().all(), "null watershed_geom_source after join"
    df = df.drop(columns="_merge")

    # ---- separate, documented QA flags + union ----
    df["nlcd_collection"] = COLLECTION
    df["valid_area_km2"] = (df.n_nlcd_pixels * 900.0 / 1e6).round(4)  # 30 m equal-area pixels
    df["geom_low_confidence"] = df.low_confidence.fillna(False).astype(bool)
    df["low_pixel_support"] = df.n_nlcd_pixels < LOWPX
    df["partial_coverage"] = df.valid_coverage_frac < PARTIAL
    df["low_confidence"] = df.geom_low_confidence | df.low_pixel_support | df.partial_coverage

    # ---- column order: keys -> provenance/QA -> class % -> impervious ----
    order = (["gage_id", "canon_id", "year", "nlcd_collection", "watershed_geom_source",
              "n_nlcd_pixels", "valid_area_km2", "valid_coverage_frac",
              "geom_low_confidence", "low_pixel_support", "partial_coverage", "low_confidence"]
             + classcols + ["impervious_mean_pct"])
    df = df[order].sort_values(["gage_id", "year"]).reset_index(drop=True)

    # completeness assertions
    ng, ny = df.gage_id.nunique(), df.year.nunique()
    assert len(df) == ng * ny, f"row count {len(df)} != {ng} gages x {ny} yr"
    assert not df.duplicated(["gage_id", "year"]).any()

    # write the finalized CONUS product date-stamped; do NOT clobber the raw pipeline intermediate
    # (watershed_nlcd_annual.parquet keeps all 6,164 gages incl. Alaska, so finalize stays re-runnable)
    df.to_parquet(f"{OUT}/watershed_nlcd_annual_{a.date}.parquet")
    df.to_csv(f"{OUT}/watershed_nlcd_annual_{a.date}.csv", index=False)

    # ---- all-gage QA companion (traceability: covered + excluded) ----
    qa = g[["gage_id", "canon_id", "latitude", "longitude", "watershed_geom_source"]].copy()
    qa["in_nlcd_footprint"] = ~qa.gage_id.isin(oof_ids)
    qa = qa.merge(maxpix.rename("max_n_nlcd_pixels"), on="gage_id", how="left")
    qa.to_csv(f"{OUT}/watershed_nlcd_gage_qa_{a.date}.csv", index=False)

    # ---- data dictionary ----
    rows = [
     ("gage_id","string","zero-padded gage id; joins to streamflow_signatures"),
     ("canon_id","string","zero-stripped canonical id; joins to combined_watershed_metadata"),
     ("year","int",f"calendar year ({YEARS[0]}-{YEARS[1]}), Annual NLCD land-cover mapping year"),
     ("nlcd_collection","string",f"source product version = Annual NLCD Collection 1 {COLLECTION}"),
     ("watershed_geom_source","string","basin polygon source: gagesii | wsc_eccc | hydrobasins"),
     ("n_nlcd_pixels","count","COVERAGE-WEIGHTED effective count of valid (non-fill) 30 m pixels in the basin; may be fractional"),
     ("valid_area_km2","km2","n_nlcd_pixels x 900 m2 (valid classified area; equal-area grid)"),
     ("valid_coverage_frac","0-1","valid pixels / total intersecting pixels (fill=250 excluded from the denominator); ~1 for CONUS basins"),
     ("geom_low_confidence","bool","basin polygon flagged low-confidence upstream (HydroBasins fallback / area outlier)"),
     ("low_pixel_support","bool",f"n_nlcd_pixels < {int(LOWPX)} (~0.09 km2; class %s less robust)"),
     ("partial_coverage","bool",f"valid_coverage_frac < {PARTIAL} (basin partly outside the CONUS raster footprint)"),
     ("low_confidence","bool","documented UNION of geom_low_confidence | low_pixel_support | partial_coverage"),
     ("impervious_mean_pct","% (0-100)","basin coverage-weighted MEAN of the Fractional Impervious Surface layer over valid [0,100] pixels; a continuous sub-pixel measure, NOT a land-cover class (do not add to class %s)"),
    ]
    for r in leg.itertuples():
        rows.append((r.output_column, "% (0-100)",
                     f"NLCD land cover class {int(r.code)} \"{CLASS_DESC.get(int(r.code), r.class_name)}\"; "
                     f"coverage-weighted % of the basin's VALID classified pixels (the 16 class %s sum to ~100)"))
    dd = pd.DataFrame(rows, columns=["column", "unit", "description"])
    caveats = pd.DataFrame([
     ("CAVEAT_footprint","note","CONUS only; 45 Alaska gages excluded to nlcd_out_of_footprint (see companion CSV). Non-CONUS gages simply have no NLCD rows."),
     ("CAVEAT_annual_noise","note","Annual NLCD uses a deep-learning + change-detection model; raw year-to-year deltas can reflect MODEL NOISE, not real change. Documented shrub/scrub<->grassland confusion (desert SW) and developed/barren-over-water artifacts near the Great Lakes; verified here as strong shrub<->grass year-to-year anti-correlation with 77-90 pp single-year swings in some basins. Prefer multi-year-smoothed change."),
     (f"CAVEAT_update_seam","note",f"C1V2 years {ANNUAL_UPDATE_YEARS} come from the rule-based 'annual update' path (not the full baseline model) — a methodological seam; treat the record tail cautiously."),
     ("CAVEAT_developed_vs_impervious","note","developed classes (21-24) and impervious_mean_pct are related but NON-ADDITIVE (measured corr ~0.96); impervious is continuous sub-pixel, the developed classes are categorical."),
    ], columns=["column","unit","description"])
    pd.concat([dd, caveats], ignore_index=True).to_csv(f"{OUT}/watershed_nlcd_dictionary.csv", index=False)

    # ---- machine-readable provenance metadata ----
    meta = {
     "product": "Per-watershed Annual NLCD (land cover % + fractional impervious surface)",
     "source_product": "Annual NLCD Collection 1", "collection_version": COLLECTION,
     "coverage_years": list(YEARS), "region": "CONUS (US gages; Alaska/Hawaii excluded)",
     "native_resolution_m": 30, "native_crs": "Albers Equal Area on WGS84 datum (Landsat ARD CONUS Albers; no standard EPSG code)",
     "nodata": 250, "class_pct_denominator": "valid (non-fill) classified pixels; 16 classes sum to ~100",
     "impervious": "coverage-weighted mean of Fractional Impervious Surface over valid [0,100] pixels",
     "extraction": "coverage-weighted exactextract over basins reprojected to the granule CRS",
     "geometry_source": "watershed_polygons_26jun2026", "annual_update_seam_years": ANNUAL_UPDATE_YEARS,
     "n_gages_conus": int(ng), "n_gages_excluded_alaska": len(oof_ids), "years": int(ny),
     "rows": int(len(df)), "finalize_date": a.date,
     "staged_source_mosaics": f"{S3_DEST}/temp_lulc_conus/ (TEMPORARY — delete after QA)",
    }
    with open(f"{OUT}/watershed_nlcd_metadata_{a.date}.json", "w") as f:
        json.dump(meta, f, indent=2)

    print(f"finalized: {len(df):,} rows, {ng} CONUS gages, {ny} years, {len(classcols)} class cols", flush=True)
    print(f"  flags: geom_low_conf={int(df.geom_low_confidence.sum())} low_px={int(df.low_pixel_support.sum())} "
          f"partial_cov={int(df.partial_coverage.sum())} low_confidence(union)={int(df.low_confidence.sum())} basin-years", flush=True)
    print(f"  -> watershed_nlcd_annual_{a.date}.{{parquet,csv}} + _dictionary.csv + _metadata.json "
          f"+ out_of_footprint + gage_qa CSVs", flush=True)

    if a.upload:
        arts = [f"watershed_nlcd_annual_{a.date}.parquet", f"watershed_nlcd_annual_{a.date}.csv",
                "watershed_nlcd_dictionary.csv", f"watershed_nlcd_metadata_{a.date}.json",
                f"nlcd_out_of_footprint_{a.date}.csv", f"watershed_nlcd_gage_qa_{a.date}.csv"]
        for fn in arts:
            subprocess.run(["aws", "s3", "cp", f"{OUT}/{fn}", f"{S3_DEST}/{fn}"], check=True)
        print(f"uploaded {len(arts)} artifacts to {S3_DEST}/", flush=True)
    else:
        print("  (S3 upload skipped — run with --upload after QA/QC via the explorer)", flush=True)


if __name__ == "__main__":
    main()

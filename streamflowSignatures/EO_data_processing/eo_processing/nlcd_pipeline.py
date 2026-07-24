"""Per-watershed annual NLCD (Annual NLCD Collection 1, C1V2, CONUS 1985-2025) land-cover %
+ basin-mean fractional impervious surface, for every US (gage_type=='USGS') watershed.

Mirrors lulc_pipeline.py (MODIS), adapted for NLCD (design + pilot: EO_data_processing/README_NLCD.md):
  - Source = 30 m COG mosaics staged in-region at s3://.../streamflow/temp_lulc_conus/
    (Annual_NLCD_{LndCov,FctImp}_{year}_CU_C1V2.tif). Per year: download both mosaics locally
    (~13 s, in-region), extract, delete (peak disk ~5 GB/worker).
  - CRS = the granule's WGS84-Albers WKT read from src.crs (NOT EPSG:5070); basins reprojected to it.
  - Land cover: coverage-weighted exact_extract unique/frac/count -> per-class % over valid (non-250)
    pixels, reindexed to the committed 16-class manifest (nlcd_legends.csv). Unknown code -> hard-fail.
  - Impervious: coverage-weighted mean over [0,100]-valid pixels. C1 guard: NLCD FIS has a documented
    bug (untruncated regression -> values >100 / uint8 underflow that are NOT the 250 fill). Pilot
    showed V2 clean (0 out-of-range across all 6,164 basins), but if any basin's max>100 we recompute
    those basins on a [0,100]-clipped window (memory-safe; dormant when clean).
  - One row per (gage_id, year); per-year parquet checkpoint (ATOMIC write), skip-if-exists (restartable).
    AK/HI/out-of-CONUS USGS basins extract ~all fill -> valid_coverage_frac~0; dropped/reported in
    nlcd_finalize.py, not here.

Integrity gates (added after adversarial review, 22 Jul 2026):
  - Per-year output validated before write (id-set == basins, unique, year, no null counts, class %s
    sum to 100 over covered basins / 0 over uncovered, impervious in [0,100], coverage-frac sane).
  - Checkpoints written atomically (tmp + os.replace) so an interrupted write is never seen as cached.
  - Final table built ONLY from the requested years' validated checkpoints; if any requested year is
    missing/invalid the final table is NOT (re)written and the process exits nonzero.
  - Downloads validated against S3 ContentLength (truncation-proof); both rasters schema-checked.

Usage:
  python nlcd_pipeline.py --years 2023               # one year (validation)
  python nlcd_pipeline.py --years 1985-2025 --workers 4
"""
import os, sys, glob, time, shutil, argparse
os.environ.setdefault("PROJ_LIB", "/home/sagemaker-user/.conda/envs/geo/share/proj")
os.environ.setdefault("PROJ_DATA", "/home/sagemaker-user/.conda/envs/geo/share/proj")
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed

BUCKET = "climate-ai-data-science-shiny-app-data"
S3_DIR = "streamflow/temp_lulc_conus"
KEY = {"LndCov": "Annual_NLCD_LndCov_{y}_CU_C1V2.tif", "FctImp": "Annual_NLCD_FctImp_{y}_CU_C1V2.tif"}
ROOT = "/home/sagemaker-user/streamflowSignatures/data_out/eo_nlcd"
HERE = os.path.dirname(os.path.abspath(__file__))
FILL = 250
PIX_M2 = 900.0  # 30 m x 30 m
AWS_ENV = ("AWS_ACCESS_KEY_ID", "AWS_SECRET_ACCESS_KEY", "AWS_SESSION_TOKEN")

# ---- committed class manifest (validated at import) ----
LEG = pd.read_csv(os.path.join(HERE, "nlcd_legends.csv"))
if LEG.code.duplicated().any() or LEG.output_column.duplicated().any():
    raise ValueError("nlcd_legends.csv has duplicate code / output_column")
if LEG[["code", "output_column"]].isna().any().any() or (LEG.code == FILL).any():
    raise ValueError("nlcd_legends.csv has null fields or the fill code 250")
CODES = [int(c) for c in LEG.code]
COL = {int(r.code): r.output_column for r in LEG.itertuples()}
VALID = set(CODES)
CLASSCOLS = [COL[c] for c in CODES]
ORDERED = ["gage_id", "canon_id", "year"] + CLASSCOLS + ["n_nlcd_pixels", "impervious_mean_pct", "valid_coverage_frac"]


def _drop_static_creds():
    """Use the auto-refreshing container role, not static env creds that expire ~30 min mid-run."""
    for k in AWS_ENV:
        os.environ.pop(k, None)


def _download(year, dstdir):
    """Download the year's two mosaics from the in-region temp S3 (our bucket, not requester-pays).
    Validates each file against the S3 ContentLength (truncation-proof); raises after retries."""
    import boto3
    s3 = boto3.client("s3", region_name="us-east-2")
    out = {}
    for layer, tmpl in KEY.items():
        fn = tmpl.format(y=year)
        key = f"{S3_DIR}/{fn}"
        dst = os.path.join(dstdir, fn)
        size = int(s3.head_object(Bucket=BUCKET, Key=key)["ContentLength"])
        ok, last_err = False, ""
        for attempt in range(4):
            try:
                s3.download_file(BUCKET, key, dst)
                if os.path.exists(dst) and os.path.getsize(dst) == size:
                    ok = True
                    break
                last_err = f"size {os.path.getsize(dst) if os.path.exists(dst) else 0}/{size}"
            except Exception as e:
                last_err = f"{type(e).__name__}: {str(e)[:80]}"
            time.sleep(5 * (attempt + 1))
        if not ok:
            got = os.path.getsize(dst) if os.path.exists(dst) else 0
            raise IOError(f"{year} {layer}: download failed ({got}/{size} bytes; last: {last_err})")
        out[layer] = dst
    return out


def _check_raster(path, want_crs=None):
    """Schema-check a mosaic; return its CRS. Explicit exceptions (not asserts -> survive python -O)."""
    with rasterio.open(path) as ds:
        if ds.count != 1:
            raise ValueError(f"{os.path.basename(path)}: {ds.count} bands != 1")
        if str(ds.dtypes[0]) != "uint8":
            raise ValueError(f"{os.path.basename(path)}: dtype {ds.dtypes[0]} != uint8")
        if ds.nodata != FILL:
            raise ValueError(f"{os.path.basename(path)}: nodata {ds.nodata} != {FILL}")
        if want_crs is not None and ds.crs != want_crs:
            raise ValueError(f"{os.path.basename(path)}: CRS differs from the land-cover granule")
        return ds.crs


def _clipped_mean(fi_path, basins_flagged):
    """C1 fallback (dormant when the FIS raster is clean): coverage-weighted mean over [0,100]-clipped
    pixels, per flagged basin, on a small windowed MemoryFile (memory-safe; never reads the 16 GB raster
    whole). Window is rounded + intersected with the raster so read array and transform always agree."""
    from rasterio.io import MemoryFile
    from exactextract import exact_extract
    out = {}
    with rasterio.open(fi_path) as fi:
        full = Window(0, 0, fi.width, fi.height)
        for row in basins_flagged.itertuples():
            geom = row.geometry
            win = fi.window(*geom.bounds).round_offsets().round_lengths().intersection(full)
            if win.width <= 0 or win.height <= 0:
                out[row.gage_id] = np.nan
                continue
            arr = fi.read(1, window=win)
            tr = fi.window_transform(win)
            clip = np.where((arr > 100) & (arr != FILL), FILL, arr).astype("uint8")  # OOB/underflow -> fill
            prof = dict(driver="GTiff", height=clip.shape[0], width=clip.shape[1], count=1,
                        dtype="uint8", crs=fi.crs, transform=tr, nodata=FILL)
            g1 = gpd.GeoDataFrame({"gage_id": [row.gage_id]}, geometry=[geom], crs=basins_flagged.crs)
            with MemoryFile() as mf:
                with mf.open(**prof) as ds:
                    ds.write(clip, 1)
                with mf.open() as ds:
                    r = exact_extract(ds, g1, ["mean"], include_cols=["gage_id"], output="pandas")
            out[row.gage_id] = float(r["mean"].iloc[0]) if len(r) and pd.notna(r["mean"].iloc[0]) else np.nan
    return out


def _validate_year(out, basins, year):
    """Hard-fail (raise) if the year's table is not publish-clean. Returns the max class-sum deviation."""
    ids = set(basins.gage_id)
    if not out.gage_id.is_unique:
        raise ValueError("duplicate gage_id in year output")
    if set(out.gage_id) != ids or len(out) != len(basins):
        raise ValueError(f"gage_id set mismatch ({len(out)} rows vs {len(basins)} basins)")
    if not (out.year == year).all():
        raise ValueError("year column mismatch")
    if out.n_nlcd_pixels.isna().any():
        raise ValueError("null n_nlcd_pixels (land-cover extraction failure)")
    covered = out.n_nlcd_pixels.to_numpy() > 0
    s = out[CLASSCOLS].sum(axis=1).to_numpy()
    dev = 0.0
    if covered.any():
        dev = float(np.abs(s[covered] - 100.0).max())
        if dev > 0.5:
            raise ValueError(f"class %% sum != 100 on covered basins (max dev {dev:.4f})")
    if (~covered).any() and (s[~covered] > 0.5).any():
        raise ValueError("nonzero class %% on a zero-coverage basin")
    # wholesale impervious-extraction failure guard (tolerate rare edge NaN where LC covers but FIS is fill)
    n_cov = int(covered.sum())
    if n_cov:
        imp_null_cov = int(pd.isna(out.impervious_mean_pct.to_numpy()[covered]).sum())
        if imp_null_cov > max(5, 0.01 * n_cov):
            raise ValueError(f"impervious null on {imp_null_cov}/{n_cov} covered basins (extraction failure)")
    imp = out.impervious_mean_pct.to_numpy()
    fin = np.isfinite(imp)
    if fin.any() and ((imp[fin] < 0).any() or (imp[fin] > 100).any()):
        raise ValueError("impervious_mean_pct outside [0,100] after the C1 guard")
    vf = out.valid_coverage_frac.to_numpy()
    if not np.isfinite(vf).all() or (vf < 0).any() or (vf > 1).any():
        raise ValueError("valid_coverage_frac not in [0,1]")
    return dev


def process_year(args):
    year, basins_path, outdir, scratch = args
    outp = os.path.join(outdir, f"nlcd_{year}.parquet")
    if os.path.exists(outp):
        return year, outp, "cached"
    _drop_static_creds()
    tmpdir = os.path.join(scratch, f"_tmp_{year}")                             # bulky mosaics -> big scratch fs
    tmp_out = os.path.join(outdir, f".nlcd_{year}.parquet.tmp.{os.getpid()}")  # atomic write on the outdir fs
    shutil.rmtree(tmpdir, ignore_errors=True)  # clear any stale temp from a hard-killed prior run
    os.makedirs(tmpdir, exist_ok=True)
    try:
        from exactextract import exact_extract
        basins = gpd.read_parquet(basins_path)
        basins["gage_id"] = basins["gage_id"].astype(str)
        if not basins.gage_id.is_unique:
            raise ValueError("duplicate gage_id in basins geometry")

        paths = _download(year, tmpdir)  # raises on failure/truncation -> year skipped, not written
        lc_path, fi_path = paths["LndCov"], paths["FctImp"]
        lc_crs = _check_raster(lc_path)
        _check_raster(fi_path, want_crs=lc_crs)
        if basins.crs != lc_crs:
            basins = basins.to_crs(lc_crs)  # defensive: extraction requires basins in the raster CRS

        out = basins[["gage_id", "canon_id"]].copy()
        out["year"] = year

        # ---- LAND COVER: per-class %, sum-to-100 over valid, unknown-code hard-fail ----
        lc = exact_extract(lc_path, basins, ["unique", "frac", "count"],
                           include_cols=["gage_id"], output="pandas")
        lc["gage_id"] = lc["gage_id"].astype(str)
        recs, unknown = [], set()
        for gid, u, fr in zip(lc["gage_id"], lc["unique"], lc["frac"]):
            for code, f in zip(u, fr):
                code = int(code)
                if code == FILL:
                    continue
                if code in VALID:
                    recs.append((gid, COL[code], float(f) * 100.0))
                else:
                    unknown.add(code)
        if unknown:
            raise ValueError(f"{year}: NLCD codes not in nlcd_legends.csv manifest: {sorted(unknown)}")
        if recs:
            long = pd.DataFrame(recs, columns=["gage_id", "col", "pct"])
            wide = long.pivot_table(index="gage_id", columns="col", values="pct", fill_value=0.0).reset_index()
            out = out.merge(wide, on="gage_id", how="left", validate="one_to_one")
        for c in CLASSCOLS:
            if c not in out.columns:
                out[c] = 0.0
        out[CLASSCOLS] = out[CLASSCOLS].fillna(0.0)
        out = out.merge(lc[["gage_id", "count"]].rename(columns={"count": "n_nlcd_pixels"}),
                        on="gage_id", how="left", validate="one_to_one")

        # ---- IMPERVIOUS: coverage-weighted mean over [0,100]; C1 guard (recompute flagged basins clipped) ----
        fi = exact_extract(fi_path, basins, ["mean", "max"], include_cols=["gage_id"], output="pandas")
        fi["gage_id"] = fi["gage_id"].astype(str)
        fi = fi.rename(columns={"mean": "impervious_mean_pct", "max": "_imp_max"})
        n_oob = int((fi["_imp_max"].astype(float) > 100).sum())
        if n_oob:  # documented FIS bug present this year -> recompute flagged basins on a [0,100]-clipped window
            flagged = fi.loc[fi["_imp_max"].astype(float) > 100, "gage_id"].tolist()
            fixed = _clipped_mean(fi_path, basins[basins.gage_id.isin(flagged)])
            if set(fixed) != set(flagged):
                raise ValueError(f"{year}: C1 recompute covered {len(fixed)}/{len(flagged)} flagged basins")
            mask = fi.gage_id.isin(flagged)
            fi.loc[mask, "impervious_mean_pct"] = fi.loc[mask, "gage_id"].map(fixed)
        out = out.merge(fi[["gage_id", "impervious_mean_pct"]], on="gage_id", how="left", validate="one_to_one")

        # ---- QA: valid_coverage_frac = valid pixels / total intersecting pixels (area/900 m^2, equal-area) ----
        area = basins.set_index("gage_id").geometry.area
        am = out["gage_id"].map(area).to_numpy()
        if not np.isfinite(am).all() or (am <= 0).any():
            raise ValueError("non-finite or non-positive basin area")
        raw = out["n_nlcd_pixels"].to_numpy() * PIX_M2 / am
        if np.nanmax(raw) > 1.1:  # >1 beyond FP tolerance => count/area inconsistency, not a clip candidate
            raise ValueError(f"valid_coverage_frac raw max {np.nanmax(raw):.3f} > 1.1 (count/area inconsistency)")
        out["valid_coverage_frac"] = np.clip(raw, 0.0, 1.0)

        out = out.reindex(columns=ORDERED)
        dev = _validate_year(out, basins, year)

        out.to_parquet(tmp_out)      # tmp_out is on the outdir fs -> os.replace is atomic (same device)
        os.replace(tmp_out, outp)
        n_cov = int((out.n_nlcd_pixels > 0).sum())
        return year, outp, f"OK {len(out)} basins ({n_cov} in-CONUS), sum100 dev<={dev:.1e}, imp_oob={n_oob}"
    except Exception as e:
        return year, outp, f"ERROR {year}: {type(e).__name__} {str(e)[:120]}"
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
        if os.path.exists(tmp_out):
            try:
                os.remove(tmp_out)
            except OSError:
                pass


def _load_checkpoint(cp, year, expected_ids):
    """Load + validate a per-year checkpoint for the final aggregation; raise if not publish-clean."""
    df = pd.read_parquet(cp)
    df["gage_id"] = df["gage_id"].astype(str)
    if list(df.columns) != ORDERED:
        raise ValueError(f"{os.path.basename(cp)}: schema/column-order mismatch")
    if not df.gage_id.is_unique or set(df.gage_id) != expected_ids:
        raise ValueError(f"{os.path.basename(cp)}: gage_id set mismatch")
    if not (df.year == year).all():
        raise ValueError(f"{os.path.basename(cp)}: unexpected year values")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--geometry", default=f"{ROOT}/watershed_polygons_26jun2026.parquet")
    ap.add_argument("--outdir", default=f"{ROOT}/nlcd_out")
    ap.add_argument("--scratch", default="/tmp/nlcd_scratch",
                    help="big-filesystem dir for transient per-year mosaic downloads; must NOT be the small "
                         "home volume (~1.4 GB LndCov + ~1.0 GB FctImp per concurrent worker)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(a.scratch, exist_ok=True)
    y0, y1 = (int(x) for x in a.years.split("-")) if "-" in a.years else (int(a.years), int(a.years))
    years = list(range(y0, y1 + 1))

    # preflight: scratch must hold ~2.6 GB per concurrent worker; outdir only needs room for tiny checkpoints
    need_gb = a.workers * 2.6
    free_scratch = shutil.disk_usage(a.scratch).free / 1e9
    if free_scratch < need_gb:
        raise SystemExit(f"scratch {a.scratch} has {free_scratch:.1f} GB free, needs ~{need_gb:.1f} GB for "
                         f"{a.workers} workers -> point --scratch at a bigger fs or lower --workers")
    free_out = shutil.disk_usage(a.outdir).free / 1e9
    if free_out < 0.5:
        raise SystemExit(f"outdir {a.outdir} has only {free_out:.2f} GB free (need room for checkpoints)")
    print(f"disk: scratch {a.scratch} {free_scratch:.1f} GB free | outdir {a.outdir} {free_out:.1f} GB free", flush=True)

    _drop_static_creds()  # parent uses the role too (consistency with workers)
    # US basins only (gage_type=='USGS'); reproject once to the granule WGS84-Albers WKT
    g = gpd.read_parquet(a.geometry)
    g["gage_id"] = g.gage_id.astype(str)
    g = g[g.gage_type == "USGS"][["gage_id", "canon_id", "geometry"]].copy()
    if not g.gage_id.is_unique:
        raise ValueError("duplicate gage_id among US basins")
    import boto3
    from rasterio.session import AWSSession
    with rasterio.Env(AWSSession(boto3.Session()), AWS_REGION="us-east-2"):
        with rasterio.open(f"/vsis3/{BUCKET}/{S3_DIR}/{KEY['LndCov'].format(y=years[0])}") as ds:
            rcrs = ds.crs
    g = g.to_crs(rcrs)
    basins_path = os.path.join(a.outdir, "_basins_aea.parquet")
    g.to_parquet(basins_path)
    expected_ids = set(g.gage_id)
    print(f"NLCD | US basins {len(g)} | years {years[0]}-{years[-1]} | workers {a.workers} | "
          f"class cols {len(CLASSCOLS)} | CRS={rcrs.to_wkt()[:40]}...", flush=True)

    tasks = [(y, basins_path, a.outdir, a.scratch) for y in years]
    status = {}
    if a.workers > 1:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for f in as_completed([ex.submit(process_year, x) for x in tasks]):
                yr, p, msg = f.result()
                status[yr] = msg
                print(f"  nlcd_{yr}.parquet: {msg}", flush=True)
    else:
        for x in tasks:
            yr, p, msg = process_year(x)
            status[yr] = msg
            print(f"  nlcd_{yr}.parquet: {msg}", flush=True)

    # ---- publish ONLY if every requested year has a valid checkpoint (scoped to requested years) ----
    parts, bad = [], []
    for y in years:
        cp = os.path.join(a.outdir, f"nlcd_{y}.parquet")
        if not os.path.exists(cp):
            bad.append((y, "no checkpoint"))
            continue
        try:
            parts.append(_load_checkpoint(cp, y, expected_ids))
        except Exception as e:
            bad.append((y, f"{type(e).__name__}: {e}"))
    if bad:
        print(f"\nINCOMPLETE: {len(bad)}/{len(years)} years missing/invalid -> final table NOT written:", flush=True)
        for y, why in bad:
            print(f"  {y}: {why}  (status: {status.get(y, 'n/a')})", flush=True)
        print("Re-run the SAME command to fill failed years (valid checkpoints are skipped).", flush=True)
        sys.exit(1)

    out = pd.concat(parts, ignore_index=True).sort_values(["gage_id", "year"])
    if out.duplicated(["gage_id", "year"]).any():
        raise ValueError("duplicate (gage_id, year) after concatenation")
    if len(out) != len(years) * len(expected_ids):
        raise ValueError(f"row count {len(out)} != {len(years)} years x {len(expected_ids)} basins")
    final = os.path.join(a.outdir, "watershed_nlcd_annual.parquet")
    tmp_final = final + f".tmp.{os.getpid()}"
    out.to_parquet(tmp_final)
    os.replace(tmp_final, final)
    print(f"FINAL: {len(out)} rows, {out.gage_id.nunique()} gages, {out.year.nunique()} years -> {final}", flush=True)


if __name__ == "__main__":
    main()

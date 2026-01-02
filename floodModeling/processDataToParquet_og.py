import io
import re
import json
import boto3
import numpy as np
import pandas as pd
from typing import Tuple, Dict, List
from osgeo import gdal
import pyarrow as pa
import pyarrow.parquet as pq

# ---------- S3 helpers ----------

def parse_s3_uri(uri: str) -> Tuple[str, str]:
    assert uri.startswith("s3://")
    p = uri[5:]
    return p.split("/", 1)[0], p.split("/", 1)[1]

def to_vsis3(uri: str) -> str:
    b, k = parse_s3_uri(uri)
    return f"/vsis3/{b}/{k}"

def s3_exists(s3_uri: str) -> bool:
    s3 = boto3.client('s3')
    b, k = parse_s3_uri(s3_uri)
    try:
        s3.head_object(Bucket=b, Key=k)
        return True
    except Exception:
        return False

def s3_read_csv(s3_uri: str, dtype=None) -> pd.DataFrame:
    s3 = boto3.client('s3')
    b, k = parse_s3_uri(s3_uri)
    obj = s3.get_object(Bucket=b, Key=k)
    return pd.read_csv(io.BytesIO(obj['Body'].read()), dtype=dtype)

def s3_read_text(s3_uri: str) -> str:
    s3 = boto3.client('s3')
    b, k = parse_s3_uri(s3_uri)
    obj = s3.get_object(Bucket=b, Key=k)
    return obj['Body'].read().decode('utf-8')

# ---------- Generic helpers ----------

def fmt4(x: float) -> str:
    return f"{float(x):.4f}"

def sanitize_name(name: str) -> str:
    name = (name or "").strip()
    name = re.sub(r"[^\w\-]+", "_", name)
    name = re.sub(r"_+", "_", name)
    return name

def parse_band_date(desc: str) -> pd.Timestamp:
    # Expect 'DSWx_S1_YYYY-MM-DD' and tolerate minor variations
    m = re.search(r"\d{4}-\d{2}-\d{2}", desc or "")
    if not m:
        raise ValueError(f"Could not parse date from band description: '{desc}'")
    return pd.to_datetime(m.group(0)).normalize()

# Rolling-window specs
STREAMFLOW_WINDOWS = [2, 4, 8, 16, 32]   # mean over current day + previous (w-1) days
CLIMATE_WINDOWS    = [3, 9]              # for t2m and tp (not ssrd)

def add_rolling_ts_features(df: pd.DataFrame) -> pd.DataFrame:
    """
    Given a timeseries DataFrame indexed by Date with columns:
      - streamflow_pred, t2m, ssrd, tp
    add rolling mean columns:
      - streamflow_pred_{w} for w in STREAMFLOW_WINDOWS
      - t2m_{w}, tp_{w} for w in CLIMATE_WINDOWS
    Uses min_periods = window size (NaN until enough history exists).
    """
    if df.empty:
        return df
    df = df.sort_index()
    # Streamflow rolling means
    if "streamflow_pred" in df.columns:
        for w in STREAMFLOW_WINDOWS:
            df[f"streamflow_pred_{w}"] = df["streamflow_pred"].rolling(window=w, min_periods=w).mean()
    # Climate rolling means (t2m, tp)
    for col in ("t2m", "tp"):
        if col in df.columns:
            for w in CLIMATE_WINDOWS:
                df[f"{col}_{w}"] = df[col].rolling(window=w, min_periods=w).mean()
    return df


# Soil features list (ordering used when writing multiband rasters)
SOIL_FEATURES = [
    'CLAY_0_TO_5CM', 'CLAY_5_TO_15CM', 'CLAY_15_TO_30CM', 
    'CLAY_30_TO_60CM', 'CLAY_60_TO_100CM', 'CLAY_100_TO_200CM',
    'SAND_0_TO_5CM', 'SAND_5_TO_15CM', 'SAND_15_TO_30CM',
    'SAND_30_TO_60CM', 'SAND_60_TO_100CM', 'SAND_100_TO_200CM',
    'SILT_0_TO_5CM', 'SILT_5_TO_15CM', 'SILT_15_TO_30CM',
    'SILT_30_TO_60CM', 'SILT_60_TO_100CM', 'SILT_100_TO_200CM'
]

# Expected band-name sequence used when writting multiband TIFs (if we must synthesize)
EXPECTED_STATIC_NAME_ORDER = (
    ['elevation'] +
    # Lexicographic order of keys like elevation_anomaly_{window}x{window} when sorted during writing:
    ['elevation_anomaly_17x17', 'elevation_anomaly_37x37', 'elevation_anomaly_3x3',
     'elevation_anomaly_5x5', 'elevation_anomaly_9x9'] +
    ['HAND', 'IO_LULC'] + SOIL_FEATURES
)

# ---------- Band template discovery ----------

def try_read_sidecar_band_names(training_path: str, lat4: str, lon4: str) -> List[str]:
    # Optional sidecar if you adopt the writer patch to emit JSON
    sidecar = f"{training_path}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}_bands.json"
    if not s3_exists(sidecar):
        return []
    try:
        data = json.loads(s3_read_text(sidecar))
        names = data.get("band_names", [])
        if names and all(isinstance(n, str) and n for n in names):
            return [sanitize_name(n) for n in names]
    except Exception:
        pass
    return []

def read_band_descriptions(static_tif_s3: str) -> List[str]:
    ds = gdal.Open(to_vsis3(static_tif_s3), gdal.GA_ReadOnly)
    if ds is None or ds.RasterCount == 0:
        return []
    names = []
    for i in range(1, ds.RasterCount+1):
        d = ds.GetRasterBand(i).GetDescription()
        names.append(sanitize_name(d) if d else "")
    ds = None
    return names

def choose_canonical_template(training_path: str, watersheds: List[Tuple[str,str]]) -> List[str]:
    # 1) Prefer a sidecar JSON if found for any watershed
    for lat4, lon4 in watersheds:
        names = try_read_sidecar_band_names(training_path, lat4, lon4)
        if names:
            return names

    # 2) Otherwise pick the watershed whose static TIF has the most non-empty band descriptions
    best = None
    best_nonempty = -1
    best_total = 0
    for lat4, lon4 in watersheds:
        static_tif = f"{training_path}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"
        if not s3_exists(static_tif):
            continue
        names = read_band_descriptions(static_tif)
        if not names:
            continue
        nonempty = sum(1 for n in names if n)
        if nonempty > best_nonempty or (nonempty == best_nonempty and len(names) > best_total):
            best = names
            best_nonempty = nonempty
            best_total = len(names)

    if best:
        # Fill blanks deterministically by index
        templ = []
        for i, nm in enumerate(best, start=1):
            templ.append(nm if nm else f"static_band_{i:02d}")
        return templ

    # 3) Last resort: synthesize by band count of first available static
    for lat4, lon4 in watersheds:
        static_tif = f"{training_path}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"
        if not s3_exists(static_tif):
            continue
        ds = gdal.Open(to_vsis3(static_tif), gdal.GA_ReadOnly)
        if ds and ds.RasterCount > 0:
            C = ds.RasterCount
            ds = None
            # Take the first C names from expected order; if fewer, append generic names
            base = list(EXPECTED_STATIC_NAME_ORDER[:C])
            if len(base) < C:
                base += [f"static_band_{i:02d}" for i in range(len(base)+1, C+1)]
            return base

    raise RuntimeError("Could not derive canonical static band names (no static rasters found).")

# ---------- Static loader ----------

def load_static_as_map(static_tif_s3: str, template_names: List[str]) -> Tuple[Dict[str, np.ndarray], Tuple[int,int]]:
    """
    Read all bands and map them to template_names:
    - If the band's description matches a template name, use that name.
    - Otherwise, map by band index to template_names[i-1] (or synthesize).
    - Convert declared NoData, -9999, and -32768 to NaN.
    Returns (column_map, (H, W))
    """
    ds = gdal.Open(to_vsis3(static_tif_s3), gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f"Cannot open static raster: {static_tif_s3}")
    H, W, C = ds.RasterYSize, ds.RasterXSize, ds.RasterCount

    col_map: Dict[str, np.ndarray] = {}
    for i in range(1, C+1):
        band = ds.GetRasterBand(i)
        arr = band.ReadAsArray().astype(np.float32)
        nd = band.GetNoDataValue()
        # Map common sentinels to NaN (-9999 padding; -32768 SoilGrids)
        mask = np.zeros(arr.shape, dtype=bool)
        if nd is not None:
            mask |= (arr == float(nd))
        mask |= (arr == -9999.0)
        mask |= (arr == -32768.0)
        arr = np.where(mask, np.nan, arr)

        # Choose column name
        desc = sanitize_name(band.GetDescription() or "")
        if desc and desc in template_names:
            col = desc
        else:
            col = template_names[i-1] if i-1 < len(template_names) else f"static_band_{i:02d}"

        # In case of duplicate names (shouldn't happen), last one wins
        col_map[col] = arr.ravel()

    ds = None

    # Ensure all template columns exist
    n_pix = H * W
    for col in template_names:
        if col not in col_map:
            col_map[col] = np.full(n_pix, np.nan, dtype=np.float32)

    return col_map, (H, W)

# ---------- Time series loader ----------

def load_ts_indexed(csv_s3: str) -> pd.DataFrame:
    if s3_exists(csv_s3):
        df = s3_read_csv(csv_s3)
        if 'Date' not in df.columns:
            raise ValueError(f"{csv_s3} missing Date column")
        keep = ['Date', 'streamflow_pred', 't2m', 'ssrd', 'tp']
        df = df[[c for c in keep if c in df.columns]].copy()
        df['Date'] = pd.to_datetime(df['Date']).dt.normalize()
        df = df.set_index('Date').sort_index()
        # Add derived rolling features
        df = add_rolling_ts_features(df)
        return df
    # If missing, return empty so we emit NaNs on join as requested
    idx = pd.DatetimeIndex([], name='Date')
    cols = ['streamflow_pred','t2m','ssrd','tp']
    # Also include the rolling columns to keep schema consistent
    cols += [f"streamflow_pred_{w}" for w in STREAMFLOW_WINDOWS]
    cols += [f"t2m_{w}" for w in CLIMATE_WINDOWS]
    cols += [f"tp_{w}" for w in CLIMATE_WINDOWS]
    return pd.DataFrame(index=idx, columns=cols, dtype='float32')

# ---------- Main builder ----------

def build_flood_training_parquet_v3(
    base_lat_4: str,
    base_lon_4: str,
    bucket: str = "climate-ai-data-science-datasets",
    base_prefix: str = "arrakis-data/floodOutputs",
    parquet_relpath: str = "postProcessedTable/flood_training.parquet",
    include_pixel_coords: bool = True,
    gdal_quiet: bool = True
):
    """
    Build a single Parquet with one row per valid DSWx pixel per date across all watersheds
    defined in list_of_training_locations.csv under the base training folder.
    """
    # GDAL config
    gdal.UseExceptions()
    if gdal_quiet:
        gdal.PushErrorHandler("CPLQuietErrorHandler")
    gdal.SetConfigOption('AWS_NO_SIGN_REQUEST', 'NO')
    gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
    gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS', 'tif,tiff,geojson,gpkg')

    # Paths and inputs
    training_folder = f"trainingFolderFor_Lat_{base_lat_4}_Lon_{base_lon_4}"
    training_path = f"s3://{bucket}/{base_prefix}/{training_folder}"

    # Load watershed list (lat_str/lon_str are the canonical keys for file names) [[2]]
    locations_csv = f"{training_path}/list_of_training_locations.csv"
    if not s3_exists(locations_csv):
        raise FileNotFoundError(f"Missing list_of_training_locations.csv at {locations_csv} [[2]]")
    loc_df = s3_read_csv(locations_csv, dtype={'lat_str': str, 'lon_str': str})
    if not {'lat_str','lon_str'}.issubset(loc_df.columns):
        raise ValueError("list_of_training_locations.csv missing columns lat_str/lon_str [[2]]")
    watersheds: List[Tuple[str,str]] = []
    for _, row in loc_df[['lat_str','lon_str']].dropna().iterrows():
        watersheds.append((fmt4(float(row['lat_str'])), fmt4(float(row['lon_str']))))

    # Discover canonical static band template (names in column order)
    template_names = choose_canonical_template(training_path, watersheds)
    print(f"Static template columns ({len(template_names)}): {template_names[:10]}{' ...' if len(template_names)>10 else ''}")

    # Time series columns (constant in space per watershed/date)
    # Base 8 features
    ts_base = [
        "subwatershed_streamflow_pred", "subwatershed_t2m", "subwatershed_ssrd", "subwatershed_tp",
        "river_basin_streamflow_pred",  "river_basin_t2m",  "river_basin_ssrd",  "river_basin_tp"
    ]
    # Derived rolling features
    ts_derived = []
    # streamflow rolling
    for w in STREAMFLOW_WINDOWS:
        ts_derived += [
            f"subwatershed_streamflow_pred_{w}",
            f"river_basin_streamflow_pred_{w}",
        ]
    # climate rolling (t2m, tp)
    for w in CLIMATE_WINDOWS:
        ts_derived += [
            f"subwatershed_t2m_{w}",
            f"river_basin_t2m_{w}",
            f"subwatershed_tp_{w}",
            f"river_basin_tp_{w}",
        ]

    ts_cols = ts_base + ts_derived

    # Column order
    base_cols = ["Target", "Date", "Watershed_Loc"]
    if include_pixel_coords:
        base_cols += ["PixelRow", "PixelCol"]
    column_order = base_cols + template_names + ts_cols

    # Build Arrow schema
    fields = [
        pa.field("Target", pa.uint8()),
        pa.field("Date", pa.timestamp('ns')),
        pa.field("Watershed_Loc", pa.string())
    ]
    if include_pixel_coords:
        fields += [pa.field("PixelRow", pa.int32()), pa.field("PixelCol", pa.int32())]
    fields += [pa.field(n, pa.float32()) for n in template_names]
    fields += [pa.field(c, pa.float32()) for c in ts_cols]
    schema = pa.schema(fields)
    
    # Open Parquet writer (S3FileSystem expects 'bucket/key', not 's3://…')
    dest_uri = f"{training_path}/{parquet_relpath}"
    s3fs = pa.fs.S3FileSystem()
    bkt, key = parse_s3_uri(dest_uri)
    sink = s3fs.open_output_stream(f"{bkt}/{key}")
    writer = pq.ParquetWriter(sink, schema, compression="snappy", write_statistics=True)

    total_rows = 0
    watersheds_processed = 0

    try:
        for lat4, lon4 in watersheds:
            ws_loc = f"Lat_{lat4}_Lon_{lon4}"
            static_tif = f"{training_path}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"  # [[1]]
            flood_tif  = f"{training_path}/floodMaps/dswx_s1_timeseries_subwatershed_Lon{lon4}_Lat{lat4}.tif"  # [[1]]

            sub_csv = f"{training_path}/streamflowPredictions/subwatershed_predictions_Lat{lat4}_Lon{lon4}.csv"
            riv_csv = f"{training_path}/streamflowPredictions/river_basin_predictions_Lat{lat4}_Lon{lon4}.csv"

            if not s3_exists(static_tif):
                print(f"Skip {ws_loc}: static raster not found")
                continue
            if not s3_exists(flood_tif):
                print(f"Skip {ws_loc}: flood raster not found")
                continue

            # Load static features
            static_map, (Hs, Ws) = load_static_as_map(static_tif, template_names)

            # Load time series (constant in space)
            sub_df = load_ts_indexed(sub_csv)
            riv_df = load_ts_indexed(riv_csv)

            # Load flood DSWx TIFF and validate shape, nodata, dates [[1]]
            ds_flood = gdal.Open(to_vsis3(flood_tif), gdal.GA_ReadOnly)
            if ds_flood is None or ds_flood.RasterCount == 0:
                print(f"Skip {ws_loc}: cannot open DSWx TIFF or empty [[1]]")
                continue

            Hf, Wf = ds_flood.RasterYSize, ds_flood.RasterXSize
            if (Hf != Hs) or (Wf != Ws):
                print(f"Skip {ws_loc}: static and flood size mismatch ({Hs}x{Ws} vs {Hf}x{Wf})")
                ds_flood = None
                continue

            flood_nodata = 255  # set during DSWx writing [[1]]
            bands = ds_flood.RasterCount

            watersheds_processed += 1
            print(f"Processing {ws_loc}: {bands} DSWx dates, {Hs}x{Ws} pixels")

            for b in range(1, bands+1):
                band = ds_flood.GetRasterBand(b)
                desc = band.GetDescription()  # 'DSWx_S1_YYYY-MM-DD' [[1]]
                try:
                    date_ts = parse_band_date(desc)
                except Exception:
                    # Skip bands we can't date
                    continue

                flood_arr = band.ReadAsArray()
                valid_mask = (flood_arr != flood_nodata)
                if not valid_mask.any():
                    continue

                target_vec = (flood_arr != 0) & (flood_arr != flood_nodata)
                idx = np.where(valid_mask.ravel())[0]
                n_rows = idx.size
                if n_rows == 0:
                    continue

                # Core columns
                data_cols = {
                    "Target": target_vec.ravel()[idx].astype(np.uint8),
                    "Date": np.repeat(date_ts, n_rows),
                    "Watershed_Loc": np.repeat(ws_loc, n_rows),
                }
                if include_pixel_coords:
                    rows = idx // Ws
                    cols = idx %  Ws
                    data_cols["PixelRow"] = rows.astype(np.int32)
                    data_cols["PixelCol"] = cols.astype(np.int32)

                # Static feature columns (subset by idx)
                for name in template_names:
                    data_cols[name] = static_map[name][idx]

                # Time series values (fill NaN if missing date)
                sv = sub_df.loc[date_ts] if date_ts in sub_df.index else pd.Series(dtype='float32')
                rv = riv_df.loc[date_ts] if date_ts in riv_df.index else pd.Series(dtype='float32')
                def g(s, col):
                    try:
                        return float(s[col])
                    except Exception:
                        return np.nan

                # Base 8
                data_cols["subwatershed_streamflow_pred"] = np.full(n_rows, g(sv, 'streamflow_pred'), dtype=np.float32)
                data_cols["subwatershed_t2m"]           = np.full(n_rows, g(sv, 't2m'),            dtype=np.float32)
                data_cols["subwatershed_ssrd"]          = np.full(n_rows, g(sv, 'ssrd'),           dtype=np.float32)
                data_cols["subwatershed_tp"]            = np.full(n_rows, g(sv, 'tp'),             dtype=np.float32)
                data_cols["river_basin_streamflow_pred"] = np.full(n_rows, g(rv, 'streamflow_pred'), dtype=np.float32)
                data_cols["river_basin_t2m"]            = np.full(n_rows, g(rv, 't2m'),            dtype=np.float32)
                data_cols["river_basin_ssrd"]           = np.full(n_rows, g(rv, 'ssrd'),           dtype=np.float32)
                data_cols["river_basin_tp"]             = np.full(n_rows, g(rv, 'tp'),             dtype=np.float32)

                # Derived rolling time series features
                # streamflow rolling: windows 2,4,8,16,32
                for w in STREAMFLOW_WINDOWS:
                    data_cols[f"subwatershed_streamflow_pred_{w}"] = np.full(
                        n_rows, g(sv, f"streamflow_pred_{w}"), dtype=np.float32
                    )
                    data_cols[f"river_basin_streamflow_pred_{w}"] = np.full(
                        n_rows, g(rv, f"streamflow_pred_{w}"), dtype=np.float32
                    )

                # climate rolling (t2m, tp): windows 3, 9
                for w in CLIMATE_WINDOWS:
                    data_cols[f"subwatershed_t2m_{w}"] = np.full(
                        n_rows, g(sv, f"t2m_{w}"), dtype=np.float32
                    )
                    data_cols[f"river_basin_t2m_{w}"] = np.full(
                        n_rows, g(rv, f"t2m_{w}"), dtype=np.float32
                    )
                    data_cols[f"subwatershed_tp_{w}"] = np.full(
                        n_rows, g(sv, f"tp_{w}"), dtype=np.float32
                    )
                    data_cols[f"river_basin_tp_{w}"] = np.full(
                        n_rows, g(rv, f"tp_{w}"), dtype=np.float32
                    )

                # Build Arrow table in fixed column order and write
                arrays = []
                for col in column_order:
                    v = data_cols[col]
                    if col == "Date":
                        arr = pa.array(pd.to_datetime(v))
                    elif col == "Target":
                        arr = pa.array(v, type=pa.uint8())
                    elif col in ("Watershed_Loc",):
                        arr = pa.array(v, type=pa.string())
                    elif col in ("PixelRow","PixelCol"):
                        arr = pa.array(v, type=pa.int32())
                    else:
                        arr = pa.array(v, type=pa.float32())
                    arrays.append(arr)

                table = pa.Table.from_arrays(arrays, names=column_order)
                writer.write_table(table)
                total_rows += n_rows

            ds_flood = None
            print(f"  -> completed {ws_loc}")

        print(f"\nDONE. Total rows written: {total_rows:,} across {watersheds_processed} watersheds.")
        print(f"Parquet: {dest_uri}")

    finally:
        writer.close()
        sink.close()
        if gdal_quiet:
            gdal.PopErrorHandler()

    return dest_uri
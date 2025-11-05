"""
U-Net flood prediction (pixel-wise) with strict LOSO validation, memory-safe pipeline, and speed optimizations.

Highlights
- PROJ/GDAL configured once; warnings silenced to avoid memory growth.
- Streaming feature statistics (no large lists in memory).
- Truly lazy dataset: store static stack once per watershed; per-date dynamic vectors and labels; assemble tiles on-the-fly.
- Masked BCE + Dice loss with optional class weighting; AMP for stability/speed.
- Single-process DataLoaders (num_workers=0) to avoid GDAL/PROJ duplication.
- Overlap-aware tiling at inference; no full [C,H,W] allocation per date.
- Speed optimizations: cuDNN autotune, TF32, channels_last memory format, optional torch.compile.
- Outputs:
    assessUNetModel/unet_predictions_Lat_{lat}_Lon_{lon}.tif
    assessUNetModel/unet_confusion_matrix_{date}.tif
    assessUNetModel/unet_confusion_matrix_{date}.png
"""

import os
import io
import re
import gc
import boto3
import tempfile
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass

import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler

from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, average_precision_score
)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from osgeo import gdal
gdal.UseExceptions()


# ------------------------ GDAL/PROJ configuration ------------------------

def _configure_gdal_proj():
    """
    Ensure GDAL/PROJ can find projection databases and silence repeated warnings.
    """
    # Try pyproj data dir
    try:
        from pyproj import datadir as _pydatadir
        proj_dir = _pydatadir.get_data_dir()
        if proj_dir and os.path.isdir(proj_dir):
            os.environ["PROJ_LIB"] = proj_dir
            gdal.SetConfigOption("PROJ_LIB", proj_dir)
    except Exception:
        pass

    # Fallback common locations
    for cand in ["/opt/conda/share/proj", "/usr/share/proj", "/usr/local/share/proj"]:
        if os.path.isdir(cand):
            os.environ.setdefault("PROJ_LIB", cand)
            gdal.SetConfigOption("PROJ_LIB", os.environ["PROJ_LIB"])
            break

    # GDAL data
    try:
        import osgeo
        gdal_data_dir = os.path.join(os.path.dirname(osgeo.__file__), "data")
        if os.path.isdir(gdal_data_dir):
            os.environ["GDAL_DATA"] = gdal_data_dir
            gdal.SetConfigOption("GDAL_DATA", gdal_data_dir)
    except Exception:
        pass

    # Network-enabled grids if needed
    gdal.SetConfigOption("PROJ_NETWORK", "ON")
    # Limit internal cache (MB) and silence debug
    gdal.SetConfigOption("GDAL_CACHEMAX", "64")
    gdal.SetConfigOption("CPL_DEBUG", "OFF")
    gdal.SetConfigOption("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", "tif,tiff,csv,gpkg,geojson")

    # Silence repeated CPL warnings
    try:
        gdal.PushErrorHandler("CPLQuietErrorHandler")
    except Exception:
        pass


try:
    _configure_gdal_proj()
except Exception:
    pass


# ------------------------ Speed setup ------------------------

def speed_setup():
    """
    Enable fast kernels and lower-precision matmul where supported.
    """
    try:
        torch.backends.cudnn.benchmark = True
    except Exception:
        pass
    try:
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        # PyTorch 2.0+
        torch.set_float32_matmul_precision("high")
    except Exception:
        pass


# ------------------------ S3 helpers ------------------------

def _parse_s3_uri(uri: str) -> Tuple[str, str]:
    assert uri.startswith("s3://"), f"Not an s3 uri: {uri}"
    p = uri[5:]
    return p.split("/", 1)[0], p.split("/", 1)[1]


def _to_vsis3(uri: str) -> str:
    b, k = _parse_s3_uri(uri)
    return f"/vsis3/{b}/{k}"


def _s3_exists(s3_uri: str) -> bool:
    s3 = boto3.client("s3")
    b, k = _parse_s3_uri(s3_uri)
    try:
        s3.head_object(Bucket=b, Key=k)
        return True
    except Exception:
        return False


def _s3_read_csv(s3_uri: str, dtype=None) -> pd.DataFrame:
    s3 = boto3.client("s3")
    b, k = _parse_s3_uri(s3_uri)
    obj = s3.get_object(Bucket=b, Key=k)
    return pd.read_csv(io.BytesIO(obj["Body"].read()), dtype=dtype)


def _s3_upload(local_path: str, s3_uri: str):
    s3 = boto3.client("s3")
    b, k = _parse_s3_uri(s3_uri)
    s3.upload_file(local_path, b, k)


def _fmt4(x: float) -> str:
    return f"{float(x):.4f}"


def _parse_latlon_id(ws_id: str) -> Tuple[str, str]:
    """
    Accepts 'Lat_XX.XXXX_Lon_YY.YYYY' or 'LatXX.XXXX_LonYY.YYYY'.
    Returns (lat4, lon4) as strings with 4 decimals.
    """
    m = re.search(r"Lat[_]?([-\d\.]+)_?Lon[_]?([-\d\.]+)", ws_id)
    if not m:
        raise ValueError(f"Cannot parse watershed id: {ws_id}")
    return _fmt4(float(m.group(1))), _fmt4(float(m.group(2)))


# ------------------------ Config ------------------------

@dataclass
class UNetConfig:
    # Architecture
    in_channels: Optional[int] = None
    base_channels: int = 32
    depth: int = 4

    # Training
    epochs: int = 20
    batch_size: int = 4
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4

    # Early stopping
    patience: int = 5
    min_delta: float = 0.001

    # Data
    tile_size: int = 256
    stride_train: int = 128
    stride_test: int = 256
    min_valid_fraction: float = 0.2

    # Class balancing
    pos_weight: float = 3.0

    # System
    num_workers: int = 0    # single-process loaders to avoid GDAL/PROJ duplication
    device: Optional[str] = None
    seed: int = 42

    # Speed options
    channels_last: bool = True      # use NHWC memory format for faster convs on CUDA
    compile_model: bool = True      # torch.compile for graph capture (PyTorch 2+)

    def __post_init__(self):
        if self.device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"


# ------------------------ U-Net ------------------------

class ConvBlock(nn.Module):
    def __init__(self, in_ch: int, out_ch: int, residual: bool = False):
        super().__init__()
        self.residual = residual
        self.conv1 = nn.Conv2d(in_ch, out_ch, 3, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_ch)
        self.conv2 = nn.Conv2d(out_ch, out_ch, 3, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_ch)
        self.relu = nn.ReLU(inplace=True)
        self.match = nn.Conv2d(in_ch, out_ch, 1) if (residual and in_ch != out_ch) else None

    def forward(self, x):
        identity = x
        out = self.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        if self.residual:
            if self.match is not None:
                identity = self.match(identity)
            out = out + identity
        return self.relu(out)


class OptimizedUNet(nn.Module):
    def __init__(self, in_channels: int, base_ch: int = 32, depth: int = 4):
        super().__init__()
        self.depth = depth
        self.enc = nn.ModuleList()
        self.pool = nn.MaxPool2d(2)
        ch = base_ch
        self.enc.append(ConvBlock(in_channels, ch))
        for _ in range(depth - 1):
            self.enc.append(ConvBlock(ch, ch * 2))
            ch *= 2
        self.bottleneck = ConvBlock(ch, ch * 2)
        ch *= 2
        self.up = nn.ModuleList()
        self.dec = nn.ModuleList()
        for _ in range(depth):
            self.up.append(nn.ConvTranspose2d(ch, ch // 2, 2, stride=2))
            self.dec.append(ConvBlock(ch, ch // 2))
            ch //= 2
        self.outc = nn.Conv2d(ch, 1, 1)

    def forward(self, x):
        skips = []
        for i, block in enumerate(self.enc):
            x = block(x)
            skips.append(x)
            if i < self.depth - 1:
                x = self.pool(x)
        x = self.bottleneck(x)
        for i in range(self.depth):
            x = self.up[i](x)
            skip = skips[-(i + 1)]
            if x.shape[-2:] != skip.shape[-2:]:
                dy = skip.size(2) - x.size(2)
                dx = skip.size(3) - x.size(3)
                x = F.pad(x, [dx // 2, dx - dx // 2, dy // 2, dy - dy // 2])
            x = torch.cat([x, skip], dim=1)
            x = self.dec[i](x)
        return self.outc(x)


# ------------------------ IO helpers for rasters/CSV ------------------------

def load_watershed_list(training_folder: str) -> List[str]:
    csv_path = f"{training_folder}/list_of_training_locations.csv"
    if not _s3_exists(csv_path):
        raise FileNotFoundError(f"Missing list_of_training_locations.csv at {csv_path}")
    df = _s3_read_csv(csv_path, dtype={"lat_str": str, "lon_str": str})
    if not {"lat_str", "lon_str"}.issubset(df.columns):
        raise ValueError("list_of_training_locations.csv missing lat_str/lon_str columns")
    ws_ids = []
    for _, row in df[["lat_str", "lon_str"]].dropna().iterrows():
        lat4 = _fmt4(float(row["lat_str"])); lon4 = _fmt4(float(row["lon_str"]))
        ws_ids.append(f"Lat_{lat4}_Lon_{lon4}")
    return sorted(list(set(ws_ids)))


def _read_static_stack(static_tif_s3: str) -> Tuple[np.ndarray, List[str], Tuple[int, int], Any, Any]:
    ds = gdal.Open(_to_vsis3(static_tif_s3), gdal.GA_ReadOnly)
    if ds is None or ds.RasterCount == 0:
        raise RuntimeError(f"Cannot open static raster: {static_tif_s3}")
    H, W, C = ds.RasterYSize, ds.RasterXSize, ds.RasterCount
    stack = np.zeros((C, H, W), dtype=np.float32)
    names = []
    for i in range(1, C + 1):
        band = ds.GetRasterBand(i)
        arr = band.ReadAsArray().astype(np.float32)
        nd = band.GetNoDataValue()
        if nd is not None:
            arr[arr == float(nd)] = np.nan
        arr[arr == -9999.0] = np.nan
        arr[arr == -32768.0] = np.nan
        stack[i - 1] = arr
        names.append(band.GetDescription() or f"static_band_{i:02d}")
    gt, prj = ds.GetGeoTransform(), ds.GetProjection()
    ds = None
    return stack, names, (H, W), gt, prj


def _parse_band_date(desc: str) -> Optional[pd.Timestamp]:
    m = re.search(r"\d{4}-\d{2}-\d{2}", desc or "")
    return pd.to_datetime(m.group(0)).normalize() if m else None


def _read_flood_timeseries(flood_tif_s3: str) -> Tuple[Dict[pd.Timestamp, np.ndarray], List[pd.Timestamp], Tuple[int, int], Any, Any]:
    ds = gdal.Open(_to_vsis3(flood_tif_s3), gdal.GA_ReadOnly)
    if ds is None or ds.RasterCount == 0:
        raise RuntimeError(f"Cannot open flood TIFF: {flood_tif_s3}")
    H, W, B = ds.RasterYSize, ds.RasterXSize, ds.RasterCount
    labels_by_date: Dict[pd.Timestamp, np.ndarray] = {}
    dates: List[pd.Timestamp] = []
    for b in range(1, B + 1):
        band = ds.GetRasterBand(b)
        d = _parse_band_date(band.GetDescription() or "")
        if d is None:
            continue
        labels_by_date[d] = band.ReadAsArray().astype(np.uint8)
        dates.append(d)
    dates = sorted(list(set(dates)))
    gt, prj = ds.GetGeoTransform(), ds.GetProjection()
    ds = None
    return labels_by_date, dates, (H, W), gt, prj


def _read_timeseries(training_folder: str, lat4: str, lon4: str) -> Dict[pd.Timestamp, Dict[str, float]]:
    """
    Build date -> {feature_name: value} from subwatershed_* and river_basin_* CSVs.
    Includes rolling columns if present.
    """
    sub_csv = f"{training_folder}/streamflowPredictions/subwatershed_predictions_Lat{lat4}_Lon{lon4}.csv"
    riv_csv = f"{training_folder}/streamflowPredictions/river_basin_predictions_Lat{lat4}_Lon{lon4}.csv"
    ts_by_date: Dict[pd.Timestamp, Dict[str, float]] = {}

    def add(df: Optional[pd.DataFrame], prefix: str):
        if df is None or df.empty or "Date" not in df.columns:
            return
        df = df.copy()
        df["Date"] = pd.to_datetime(df["Date"]).dt.normalize()
        for _, row in df.iterrows():
            d = row["Date"]
            rec = ts_by_date.setdefault(d, {})
            for col, val in row.items():
                if col == "Date":
                    continue
                if isinstance(val, (int, float, np.integer, np.floating)) and pd.notna(val):
                    rec[f"{prefix}_{col}"] = float(val)

    try:
        add(_s3_read_csv(sub_csv), "subwatershed")
    except Exception:
        pass
    try:
        add(_s3_read_csv(riv_csv), "river_basin")
    except Exception:
        pass
    return ts_by_date


def _plan_dynamic_channels(ts_by_date: Dict[pd.Timestamp, Dict[str, float]]) -> List[str]:
    base = [
        "subwatershed_streamflow_pred", "subwatershed_t2m", "subwatershed_ssrd", "subwatershed_tp",
        "river_basin_streamflow_pred", "river_basin_t2m", "river_basin_ssrd", "river_basin_tp"
    ]
    extra = set()
    for _, d in ts_by_date.items():
        for k in d.keys():
            if k not in base:
                extra.add(k)
        if extra:
            break
    return base + sorted(list(extra))


# ------------------------ Streaming feature statistics ------------------------

def compute_feature_statistics(training_folder: str, watersheds: List[str]) -> Dict:
    """
    Streaming (memory-efficient) computation of means/stds for static and dynamic features.
    - Statics: per-band sums/sumsq/count over finite pixels
    - Dynamics: per-feature sums/sumsq/count across all dates
    """
    static_names_ref = None
    C_static = None
    s_stat = None
    ss_stat = None
    cnt_stat = None

    dynamic_channels: Optional[List[str]] = None
    s_dyn = None
    ss_dyn = None
    cnt_dyn = None

    for ws_id in watersheds:
        lat4, lon4 = _parse_latlon_id(ws_id)
        static_tif = f"{training_folder}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"
        flood_tif = f"{training_folder}/floodMaps/dswx_s1_timeseries_subwatershed_Lon{lon4}_Lat{lat4}.tif"
        if not _s3_exists(static_tif) or not _s3_exists(flood_tif):
            continue

        # Statics
        arr, names, _, _, _ = _read_static_stack(static_tif)
        if static_names_ref is None:
            static_names_ref = names
            C_static = len(names)
            s_stat = np.zeros(C_static, dtype=np.float64)
            ss_stat = np.zeros(C_static, dtype=np.float64)
            cnt_stat = np.zeros(C_static, dtype=np.float64)
        A = arr.reshape(C_static, -1)
        mask = np.isfinite(A)
        s_stat += np.where(mask, A, 0.0).sum(axis=1)
        ss_stat += np.where(mask, A * A, 0.0).sum(axis=1)
        cnt_stat += mask.sum(axis=1)

        # Dynamics
        ts_by_date = _read_timeseries(training_folder, lat4, lon4)
        if dynamic_channels is None and ts_by_date:
            dynamic_channels = _plan_dynamic_channels(ts_by_date)
            C_dyn = len(dynamic_channels)
            s_dyn = np.zeros(C_dyn, dtype=np.float64)
            ss_dyn = np.zeros(C_dyn, dtype=np.float64)
            cnt_dyn = np.zeros(C_dyn, dtype=np.float64)
        if dynamic_channels:
            for _, rec in ts_by_date.items():
                vals = np.array([rec.get(ch, np.nan) for ch in dynamic_channels], dtype=np.float64)
                m = ~np.isnan(vals)
                s_dyn += np.where(m, vals, 0.0)
                ss_dyn += np.where(m, vals * vals, 0.0)
                cnt_dyn += m.astype(np.float64)

        del arr, A, ts_by_date
        gc.collect()

    if static_names_ref is None:
        raise RuntimeError("No static rasters found to compute statistics.")
    if dynamic_channels is None:
        dynamic_channels = [
            "subwatershed_streamflow_pred", "subwatershed_t2m", "subwatershed_ssrd", "subwatershed_tp",
            "river_basin_streamflow_pred", "river_basin_t2m", "river_basin_ssrd", "river_basin_tp"
        ]
        C_dyn = len(dynamic_channels)
        s_dyn = np.zeros(C_dyn, dtype=np.float64)
        ss_dyn = np.zeros(C_dyn, dtype=np.float64)
        cnt_dyn = np.ones(C_dyn, dtype=np.float64)  # avoid zero div

    # Finalize statics
    sm_stat = s_stat / np.maximum(cnt_stat, 1.0)
    var_stat = (ss_stat / np.maximum(cnt_stat, 1.0)) - sm_stat ** 2
    var_stat[var_stat < 0] = 0
    sd_stat = np.sqrt(var_stat); sd_stat[sd_stat == 0] = 1.0

    # Finalize dynamics
    sm_dyn = s_dyn / np.maximum(cnt_dyn, 1.0)
    var_dyn = (ss_dyn / np.maximum(cnt_dyn, 1.0)) - sm_dyn ** 2
    var_dyn[var_dyn < 0] = 0
    sd_dyn = np.sqrt(var_dyn); sd_dyn[sd_dyn == 0] = 1.0

    return {
        "static_names": static_names_ref,
        "static_means": sm_stat.astype(np.float32),
        "static_stds": sd_stat.astype(np.float32),
        "dynamic_channels": dynamic_channels,
        "dynamic_means": sm_dyn.astype(np.float32),
        "dynamic_stds": sd_dyn.astype(np.float32),
        "n_static": len(static_names_ref),
        "n_temporal": len(dynamic_channels),
    }


# ------------------------ Truly lazy dataset ------------------------

class StaticDynamicTileDataset(Dataset):
    """
    Truly lazy dataset:
      - Store normalized static stack once per watershed.
      - Store dynamic vectors per date (length C_dyn).
      - Store Y/M per date.
      - Build flat index (ws_idx, date_idx, rr, cc, pos_frac).
      - Assemble tile features on-the-fly for __getitem__.
    """
    def __init__(self, training_folder: str, watersheds: List[str],
                 feature_stats: Dict, tile_size: int, stride: int,
                 min_valid_fraction: float, pos_tile_threshold: float = 0.02):
        super().__init__()
        self.training_folder = training_folder
        self.tile = int(tile_size)
        self.stride = int(stride)
        self.min_valid = float(min_valid_fraction)
        self.pos_thr = float(pos_tile_threshold)

        self.static_means = feature_stats["static_means"]
        self.static_stds = feature_stats["static_stds"]
        self.dyn_channels = feature_stats["dynamic_channels"]
        self.dynamic_means = feature_stats["dynamic_means"]
        self.dynamic_stds = feature_stats["dynamic_stds"]

        self.ws_statics: List[np.ndarray] = []     # [C_static,H,W]
        self.ws_dates: List[List[pd.Timestamp]] = []
        self.ws_dyn: List[List[np.ndarray]] = []   # list per ws: list of [C_dyn]
        self.ws_labels: List[List[Dict[str, np.ndarray]]] = []  # list per ws: list of dicts Y/M

        self.index: List[Tuple[int, int, slice, slice, float]] = []

        for ws_id in watersheds:
            lat4, lon4 = _parse_latlon_id(ws_id)
            static_tif = f"{training_folder}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"
            flood_tif = f"{training_folder}/floodMaps/dswx_s1_timeseries_subwatershed_Lon{lon4}_Lat{lat4}.tif"
            if not _s3_exists(static_tif) or not _s3_exists(flood_tif):
                continue

            static_stack, _, (H, W), _, _ = _read_static_stack(static_tif)
            # normalize statics
            stat = static_stack.astype(np.float32, copy=True)
            for c in range(stat.shape[0]):
                m, s = self.static_means[c], self.static_stds[c]
                band = stat[c]
                band = np.where(np.isfinite(band), band, m)
                stat[c] = (band - m) / s

            labels_by_date, dates, (Hf, Wf), _, _ = _read_flood_timeseries(flood_tif)
            if (H, W) != (Hf, Wf):
                del static_stack, stat, labels_by_date
                gc.collect()
                continue

            ts_by_date = _read_timeseries(training_folder, lat4, lon4)
            if not ts_by_date:
                del static_stack, stat, labels_by_date
                gc.collect()
                continue

            dyn_list: List[np.ndarray] = []
            lab_list: List[Dict[str, np.ndarray]] = []
            date_list: List[pd.Timestamp] = []

            for d in dates:
                feats = ts_by_date.get(d)
                if not feats:
                    continue
                vals = np.array([feats.get(ch, np.nan) for ch in self.dyn_channels], dtype=np.float32)
                vals = (np.where(np.isnan(vals), self.dynamic_means, vals) - self.dynamic_means) / self.dynamic_stds
                dyn_list.append(vals.astype(np.float32))

                Y = labels_by_date[d].astype(np.uint8)
                M = (Y != 255).astype(np.uint8)
                lab_list.append({"Y": Y, "M": M})
                date_list.append(d)

            if not dyn_list:
                del static_stack, stat, labels_by_date, ts_by_date
                gc.collect()
                continue

            ws_idx = len(self.ws_statics)
            self.ws_statics.append(stat)
            self.ws_dates.append(date_list)
            self.ws_dyn.append(dyn_list)
            self.ws_labels.append(lab_list)

            # build tile index per date
            for di, lab in enumerate(lab_list):
                M = lab["M"]; Y = lab["Y"]
                for y0 in range(0, max(H - self.tile + 1, 1), self.stride):
                    for x0 in range(0, max(W - self.tile + 1, 1), self.stride):
                        rr = slice(y0, min(y0 + self.tile, H))
                        cc = slice(x0, min(x0 + self.tile, W))
                        m = M[rr, cc]
                        if m.size == 0:
                            continue
                        valid_frac = m.sum() / m.size
                        if valid_frac < self.min_valid:
                            continue
                        pos_frac = ((Y[rr, cc] == 1) & (m == 1)).sum() / (m.sum() + 1e-9)
                        self.index.append((ws_idx, di, rr, cc, float(pos_frac)))

            del static_stack, labels_by_date, ts_by_date
            gc.collect()

    def __len__(self):
        return len(self.index)

    def __getitem__(self, i: int):
        ws_idx, di, rr, cc, _ = self.index[i]
        stat = self.ws_statics[ws_idx]
        dyn_vec = self.ws_dyn[ws_idx][di]
        Y = self.ws_labels[ws_idx][di]["Y"]
        M = self.ws_labels[ws_idx][di]["M"]

        x_static = stat[:, rr, cc]  # [C_static,h,w]
        h = rr.stop - rr.start; w = cc.stop - cc.start
        x_dyn = np.broadcast_to(dyn_vec[:, None, None], (dyn_vec.shape[0], h, w)).astype(np.float32)
        x = np.concatenate([x_static, x_dyn], axis=0)
        y = (Y[rr, cc] > 0).astype(np.float32)
        m = M[rr, cc].astype(np.float32)
        return torch.from_numpy(x), torch.from_numpy(y), torch.from_numpy(m)

    def make_weighted_sampler(self, pos_floor: float = 0.02) -> Optional[WeightedRandomSampler]:
        if len(self.index) == 0:
            return None
        weights = np.array([pf for *_, pf in self.index], dtype=np.float32)
        weights = weights + float(pos_floor)
        weights = weights / (weights.sum() + 1e-9)
        return WeightedRandomSampler(weights, num_samples=len(weights), replacement=True)


# ------------------------ Losses & Trainer ------------------------

def masked_dice_loss(logits: torch.Tensor, targets: torch.Tensor, mask: torch.Tensor, eps: float = 1e-6) -> torch.Tensor:
    # logits [B,H,W], targets [B,H,W], mask [B,H,W]
    probs = torch.sigmoid(logits)
    probs = probs * mask
    targets = targets * mask
    inter = (probs * targets).sum(dim=(1, 2))
    denom = probs.sum(dim=(1, 2)) + targets.sum(dim=(1, 2)) + eps
    dice = 2.0 * inter / denom
    return 1.0 - dice.mean()


class Trainer:
    def __init__(self, model: nn.Module, config: UNetConfig):
        self.model = model.to(config.device)
        self.config = config
        self.device = config.device
        self.optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, mode="min", patience=3, factor=0.5)
        self.use_amp = (self.device == "cuda")
        self.scaler = torch.amp.GradScaler('cuda', enabled=self.use_amp)
        self.best_loss = float("inf")
        self.patience_counter = 0

    def train_epoch(self, dataloader: DataLoader) -> Dict:
        self.model.train()
        total = 0.0; n = 0
        for xb, yb, mb in dataloader:
            # Optionally convert to channels_last for faster convs
            if self.config.channels_last and xb.device.type == "cpu":
                xb = xb.contiguous(memory_format=torch.channels_last)
            xb = xb.to(self.device, non_blocking=True)
            yb = yb.to(self.device, non_blocking=True)
            mb = mb.to(self.device, non_blocking=True)

            self.optimizer.zero_grad(set_to_none=True)
            with torch.amp.autocast('cuda', enabled=self.use_amp):
                logits = self.model(xb).squeeze(1)
                pos_w = torch.tensor(self.config.pos_weight, device=self.device) if self.config.pos_weight else None
                bce = F.binary_cross_entropy_with_logits(logits, yb, reduction="none", pos_weight=pos_w)
                bce = (bce * mb).sum() / (mb.sum() + 1e-8)
                dice = masked_dice_loss(logits, yb, mb)
                loss = 0.5 * bce + 0.5 * dice
            self.scaler.scale(loss).backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
            self.scaler.step(self.optimizer)
            self.scaler.update()
            total += float(loss.item()); n += 1

            del xb, yb, mb, logits, bce, dice, loss

        # one empty_cache at epoch boundary is enough
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return {"loss": total / max(n, 1)}

    def validate(self, dataloader: DataLoader) -> Dict:
        self.model.eval()
        total_loss = 0.0; n_batches = 0
        nbins = 256
        pos_hist = np.zeros(nbins, dtype=np.int64)
        neg_hist = np.zeros(nbins, dtype=np.int64)
        total_pos = 0
        with torch.no_grad():
            for xb, yb, mb in dataloader:
                if self.config.channels_last and xb.device.type == "cpu":
                    xb = xb.contiguous(memory_format=torch.channels_last)
                xb = xb.to(self.device, non_blocking=True)
                yb = yb.to(self.device, non_blocking=True)
                mb = mb.to(self.device, non_blocking=True)
                logits = self.model(xb).squeeze(1)
                per_pixel = F.binary_cross_entropy_with_logits(logits, yb, reduction="none")
                loss = (per_pixel * mb).sum() / (mb.sum() + 1e-8)
                total_loss += float(loss.item()); n_batches += 1

                probs = torch.sigmoid(logits)
                valid = (mb > 0)
                if valid.any():
                    y_true = yb[valid].float().cpu().numpy().astype(np.uint8)
                    y_prob = probs[valid].float().cpu().numpy().astype(np.float32)
                    idx = np.clip((y_prob * (nbins - 1)).astype(np.int32), 0, nbins - 1)
                    pos_idx = idx[y_true == 1]
                    if pos_idx.size:
                        pos_hist += np.bincount(pos_idx, minlength=nbins)
                        total_pos += int(pos_idx.size)
                    neg_idx = idx[y_true == 0]
                    if neg_idx.size:
                        neg_hist += np.bincount(neg_idx, minlength=nbins)

                del xb, yb, mb, logits, per_pixel, loss, probs

        if total_pos > 0 and (pos_hist.sum() + neg_hist.sum()) > 0:
            pos_cum = np.cumsum(pos_hist[::-1])[::-1]
            neg_cum = np.cumsum(neg_hist[::-1])[::-1]
            TP = pos_cum; FP = neg_cum; FN = total_pos - TP
            precision = TP / np.maximum(TP + FP, 1)
            recall = TP / np.maximum(TP + FN, 1)
            f1 = 2 * (precision * recall) / np.maximum(precision + recall, 1e-10)
            best_k = int(np.nanargmax(f1))
            thr = best_k / (nbins - 1)
            f1_val = float(f1[best_k])
        else:
            thr = 0.5; f1_val = 0.0

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        return {"loss": total_loss / max(n_batches, 1), "f1": f1_val, "threshold": float(thr)}


# ------------------------ Training pipeline (LOSO) ------------------------

def train_unet_loso(training_folder: str, config: Optional[UNetConfig] = None) -> Dict:
    config = config or UNetConfig()
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    _configure_gdal_proj()
    speed_setup()

    print("=" * 70)
    print("U-NET FLOOD MODEL - OPTIMIZED LOSO TRAINING")
    print("=" * 70)

    watersheds = load_watershed_list(training_folder)
    if len(watersheds) < 3:
        raise ValueError(f"Need at least 3 watersheds for LOSO, found {len(watersheds)}")

    test_watershed = watersheds[0]
    cv_watersheds = watersheds[1:]
    print(f"Test watershed: {test_watershed}")
    print(f"CV watersheds: {cv_watersheds}")

    print("\nComputing normalization statistics...")
    feature_stats = compute_feature_statistics(training_folder, cv_watersheds)
    n_channels = feature_stats["n_static"] + feature_stats["n_temporal"]
    config.in_channels = n_channels

    print(f"\nRunning LOSO CV on {len(cv_watersheds)} folds...")
    cv_thresholds: List[float] = []

    for fold, val_watershed in enumerate(cv_watersheds):
        print(f"\nFold {fold+1}/{len(cv_watersheds)}: Validating on {val_watershed}")
        train_watersheds = [w for w in cv_watersheds if w != val_watershed]

        train_dataset = StaticDynamicTileDataset(
            training_folder=training_folder,
            watersheds=train_watersheds,
            feature_stats=feature_stats,
            tile_size=config.tile_size,
            stride=config.stride_train,
            min_valid_fraction=config.min_valid_fraction
        )
        val_dataset = StaticDynamicTileDataset(
            training_folder=training_folder,
            watersheds=[val_watershed],
            feature_stats=feature_stats,
            tile_size=config.tile_size,
            stride=config.stride_test,
            min_valid_fraction=config.min_valid_fraction
        )

        if len(train_dataset) == 0 or len(val_dataset) == 0:
            print("  Warning: empty train/val tiles in this fold; using default threshold 0.5")
            cv_thresholds.append(0.5)
            del train_dataset, val_dataset
            gc.collect()
            if torch.cuda.is_available(): torch.cuda.empty_cache()
            continue

        sampler = train_dataset.make_weighted_sampler(pos_floor=0.02)  # oversample watery tiles
        pin = (config.device == "cuda")
        train_loader = DataLoader(
            train_dataset,
            batch_size=config.batch_size,
            shuffle=(sampler is None),
            sampler=sampler,
            num_workers=0,
            pin_memory=pin,
            persistent_workers=False
        )
        val_loader = DataLoader(
            val_dataset,
            batch_size=config.batch_size,
            shuffle=False,
            num_workers=0,
            pin_memory=pin,
            persistent_workers=False
        )

        model = OptimizedUNet(n_channels, config.base_channels, config.depth)

        # Speed: channels_last + optional compile
        if config.channels_last and config.device == "cuda":
            model = model.to(memory_format=torch.channels_last)
        if config.compile_model:
            try:
                model = torch.compile(model, mode="max-autotune")
            except Exception:
                pass

        model = model.to(config.device)
        trainer = Trainer(model, config)

        best_state = None
        best_val_loss = float("inf")
        best_val_threshold = 0.5
        patience_ctr = 0

        for epoch in range(config.epochs):
            tr = trainer.train_epoch(train_loader)
            va = trainer.validate(val_loader)
            trainer.scheduler.step(va["loss"])

            print(f"  Epoch {epoch+1}: Train loss={tr['loss']:.4f}, Val loss={va['loss']:.4f}, Val F1={va['f1']:.4f}")

            if va["loss"] < best_val_loss - config.min_delta:
                best_val_loss = va["loss"]
                best_val_threshold = va["threshold"]
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                patience_ctr = 0
            else:
                patience_ctr += 1
                if patience_ctr >= config.patience:
                    print(f"  Early stopping at epoch {epoch+1}")
                    break

            gc.collect()
            if torch.cuda.is_available(): torch.cuda.empty_cache()

        if best_state is not None:
            model.load_state_dict(best_state)
        cv_thresholds.append(best_val_threshold)

        del model, trainer, train_loader, val_loader, train_dataset, val_dataset
        gc.collect()
        if torch.cuda.is_available(): torch.cuda.empty_cache()

    print("\n" + "=" * 70)
    print("FINAL MODEL TRAINING")
    print("=" * 70)
    final_threshold = float(np.mean(cv_thresholds)) if cv_thresholds else 0.5
    print(f"Threshold from CV: {final_threshold:.3f}")

    train_dataset = StaticDynamicTileDataset(
        training_folder=training_folder,
        watersheds=cv_watersheds,
        feature_stats=feature_stats,
        tile_size=config.tile_size,
        stride=config.stride_train,
        min_valid_fraction=config.min_valid_fraction
    )
    if len(train_dataset) == 0:
        raise RuntimeError("No training tiles available for final training.")

    pin = (config.device == "cuda")
    sampler = train_dataset.make_weighted_sampler(pos_floor=0.02)
    train_loader = DataLoader(
        train_dataset,
        batch_size=config.batch_size,
        shuffle=(sampler is None),
        sampler=sampler,
        num_workers=0,
        pin_memory=pin,
        persistent_workers=False
    )

    final_model = OptimizedUNet(n_channels, config.base_channels, config.depth)
    if config.channels_last and config.device == "cuda":
        final_model = final_model.to(memory_format=torch.channels_last)
    if config.compile_model:
        try:
            final_model = torch.compile(final_model, mode="max-autotune")
        except Exception:
            pass
    final_model = final_model.to(config.device)

    final_trainer = Trainer(final_model, config)

    for epoch in range(config.epochs):
        try:
            tr = final_trainer.train_epoch(train_loader)
            print(f"Epoch {epoch+1}: Loss={tr['loss']:.4f}")
        except RuntimeError as e:
            if "CUDA out of memory" in str(e) and config.batch_size > 1:
                config.batch_size = max(1, config.batch_size // 2)
                print(f"  OOM encountered. Reducing batch_size to {config.batch_size} and rebuilding DataLoader...")
                train_loader = DataLoader(
                    train_dataset,
                    batch_size=config.batch_size,
                    shuffle=(sampler is None),
                    sampler=sampler,
                    num_workers=0,
                    pin_memory=pin,
                    persistent_workers=False
                )
                torch.cuda.empty_cache()
                gc.collect()
                continue
            else:
                raise
        finally:
            gc.collect()
            if torch.cuda.is_available(): torch.cuda.empty_cache()

    print(f"\nEvaluating on test watershed: {test_watershed}")
    test_metrics = evaluate_watershed(final_model, training_folder, test_watershed, feature_stats, final_threshold, config)

    print("\nTest Metrics:")
    for k, v in test_metrics.items():
        if isinstance(v, float):
            print(f"  {k}: {v:.4f}")
        else:
            print(f"  {k}: {v}")

    return {
        "model": final_model,
        "threshold": final_threshold,
        "feature_stats": feature_stats,
        "test_watershed": test_watershed,
        "test_metrics": test_metrics,
        "config": config
    }


# ------------------------ Evaluation & outputs ------------------------

def evaluate_watershed(model: nn.Module, training_folder: str,
                       watershed: str, feature_stats: Dict,
                       threshold: float, config: UNetConfig) -> Dict:
    """
    Full-frame inference per date with tilewise assembly (no full [C,H,W] per date).
    Writes:
      assessUNetModel/unet_predictions_Lat_{lat}_Lon_{lon}.tif
      assessUNetModel/unet_confusion_matrix_{date}.tif
      assessUNetModel/unet_confusion_matrix_{date}.png
    Returns aggregated metrics across dates.
    """
    _configure_gdal_proj()
    lat4, lon4 = _parse_latlon_id(watershed)
    static_means = feature_stats["static_means"]
    static_stds = feature_stats["static_stds"]
    dyn_channels = feature_stats["dynamic_channels"]
    dynamic_means = feature_stats["dynamic_means"]
    dynamic_stds = feature_stats["dynamic_stds"]

    static_tif = f"{training_folder}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"
    flood_tif = f"{training_folder}/floodMaps/dswx_s1_timeseries_subwatershed_Lon{lon4}_Lat{lat4}.tif"
    assess_folder = f"{training_folder}/assessUNetModel"

    static_stack, _, (H, W), gt, prj = _read_static_stack(static_tif)
    labels_by_date, dates, (Hf, Wf), _, _ = _read_flood_timeseries(flood_tif)
    assert (H, W) == (Hf, Wf)
    ts_by_date = _read_timeseries(training_folder, lat4, lon4)

    # Normalize statics once
    stat = static_stack.astype(np.float32, copy=True)
    for c in range(stat.shape[0]):
        m, s = static_means[c], static_stds[c]
        band = stat[c]
        band = np.where(np.isfinite(band), band, m)
        stat[c] = (band - m) / s

    ys_all: List[np.ndarray] = []
    ps_all: List[np.ndarray] = []

    with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp_pred:
        pred_tif_local = tmp_pred.name
    drv = gdal.GetDriverByName("GTiff")
    ds_out = drv.Create(pred_tif_local, W, H, len(dates), gdal.GDT_Byte,
                        options=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"])
    try: ds_out.SetGeoTransform(gt)
    except Exception: pass
    try: ds_out.SetProjection(prj)
    except Exception: pass

    for i, d in enumerate(dates, start=1):
        vals = np.array([ts_by_date.get(d, {}).get(ch, np.nan) for ch in dyn_channels], dtype=np.float32)
        vals = (np.where(np.isnan(vals), dynamic_means, vals) - dynamic_means) / dynamic_stds

        # Predict with tilewise assembly
        prob = np.zeros((H, W), dtype=np.float32)
        hits = np.zeros((H, W), dtype=np.float32)
        with torch.no_grad():
            for y0 in range(0, max(H - config.tile_size + 1, 1), config.stride_test):
                for x0 in range(0, max(W - config.tile_size + 1, 1), config.stride_test):
                    rr = slice(y0, min(y0 + config.tile_size, H))
                    cc = slice(x0, min(x0 + config.tile_size, W))
                    h = rr.stop - rr.start; w = cc.stop - cc.start
                    x_stat = stat[:, rr, cc]
                    x_dyn = np.broadcast_to(vals[:, None, None], (len(vals), h, w)).astype(np.float32)
                    x = np.concatenate([x_stat, x_dyn], axis=0)
                    xt = torch.from_numpy(x).unsqueeze(0).to(config.device)
                    logits = model(xt)
                    p = torch.sigmoid(logits).squeeze().detach().cpu().numpy().astype(np.float32)
                    prob[rr, cc] += p
                    hits[rr, cc] += 1.0
        hits[hits == 0] = 1.0
        prob = prob / hits

        lab = labels_by_date[d]
        valid = (lab != 255)
        pred_bin = np.full((H, W), 255, dtype=np.uint8)
        pred_bin[valid] = (prob[valid] >= threshold).astype(np.uint8)

        # write prediction band
        b = ds_out.GetRasterBand(i)
        b.SetNoDataValue(255)
        b.SetDescription(f"UNet_Pred_{d.strftime('%Y-%m-%d')}")
        b.WriteArray(pred_bin)
        b.FlushCache()

        # Confusion map
        cm = np.full((H, W), 255, dtype=np.uint8)
        cm_valid = np.where(
            (lab[valid] == 0),
            np.where(pred_bin[valid] == 0, 0, 3),
            np.where(pred_bin[valid] == 1, 1, 2)
        ).astype(np.uint8)
        cm[valid] = cm_valid

        # Save CM tif
        with tempfile.NamedTemporaryFile(suffix=".tif", delete=False) as tmp_cm:
            cm_local = tmp_cm.name
        ds_cm = drv.Create(cm_local, W, H, 1, gdal.GDT_Byte, options=["COMPRESS=DEFLATE", "TILED=YES"])
        try: ds_cm.SetGeoTransform(gt)
        except Exception: pass
        try: ds_cm.SetProjection(prj)
        except Exception: pass
        rb = ds_cm.GetRasterBand(1)
        rb.SetNoDataValue(255)
        rb.SetDescription(f"UNet_CM_{d.strftime('%Y-%m-%d')}_0TN_1TP_2FN_3FP")
        rb.WriteArray(cm); rb.FlushCache(); ds_cm = None
        cm_tif_s3 = f"{assess_folder}/unet_confusion_matrix_{d.strftime('%Y-%m-%d')}.tif"
        _s3_upload(cm_local, cm_tif_s3); os.remove(cm_local)

        # Save CM png
        fig, ax = plt.subplots(figsize=(12, 10))
        disp = cm.astype(float); disp[cm == 255] = np.nan
        cmap = plt.matplotlib.colors.ListedColormap(["green", "blue", "orange", "red"])
        ax.imshow(disp, cmap=cmap, vmin=0, vmax=3)
        legend = [
            mpatches.Patch(color="green", label="TN"),
            mpatches.Patch(color="blue", label="TP"),
            mpatches.Patch(color="orange", label="FN"),
            mpatches.Patch(color="red", label="FP"),
        ]
        ax.legend(handles=legend, loc="upper right", fontsize=10)
        y_true = (lab[valid] > 0).astype(np.uint8)
        y_prob = prob[valid]
        y_pred = (y_prob >= threshold).astype(np.uint8)
        try: auc = roc_auc_score(y_true, y_prob)
        except Exception: auc = np.nan
        text = (
            f"Date: {d.strftime('%Y-%m-%d')}\n"
            f"F1: {f1_score(y_true,y_pred,zero_division=0):.3f}\n"
            f"Acc: {accuracy_score(y_true,y_pred):.3f}\n"
            f"Prec: {precision_score(y_true,y_pred,zero_division=0):.3f}\n"
            f"Rec: {recall_score(y_true,y_pred,zero_division=0):.3f}\n"
            f"AUC: {auc:.3f}" if np.isfinite(auc) else f"AUC: N/A"
        )
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.8)
        ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11, va="top", bbox=props)
        ax.axis("off"); plt.tight_layout()
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_png:
            png_local = tmp_png.name
        plt.savefig(png_local, dpi=150, bbox_inches="tight"); plt.close()
        cm_png_s3 = f"{assess_folder}/unet_confusion_matrix_{d.strftime('%Y-%m-%d')}.png"
        _s3_upload(png_local, cm_png_s3); os.remove(png_local)

        ys_all.append(y_true.ravel())
        ps_all.append(y_prob.ravel())

        gc.collect()
        if torch.cuda.is_available(): torch.cuda.empty_cache()

    ds_out = None
    preds_tif_s3 = f"{assess_folder}/unet_predictions_Lat_{lat4}_Lon_{lon4}.tif"
    _s3_upload(pred_tif_local, preds_tif_s3); os.remove(pred_tif_local)
    print(f"  Wrote {preds_tif_s3}")

    if ys_all and ps_all:
        y_true = np.concatenate(ys_all)
        y_prob = np.concatenate(ps_all)
        y_pred = (y_prob >= threshold).astype(np.uint8)
        try: auc = roc_auc_score(y_true, y_prob)
        except Exception: auc = np.nan
        try: ap = average_precision_score(y_true, y_prob)
        except Exception: ap = np.nan
        metrics = {
            "threshold": float(threshold),
            "accuracy": accuracy_score(y_true, y_pred),
            "precision": precision_score(y_true, y_pred, zero_division=0),
            "recall": recall_score(y_true, y_pred, zero_division=0),
            "f1": f1_score(y_true, y_pred, zero_division=0),
            "roc_auc": float(auc) if np.isfinite(auc) else np.nan,
            "pr_auc": float(ap) if np.isfinite(ap) else np.nan,
        }
    else:
        metrics = {}

    return metrics

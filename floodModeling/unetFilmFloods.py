# unetFilmFloods.py
"""
U-Net with FiLM conditioning for flood prediction.

Enhanced with:
- Optuna hyperparameter tuning with PR-AUC objective
- SageMaker deployment capability
- Rasterio-based I/O (no system GDAL dependency)
- Authenticated S3 access via AWSSession
- Spatial k-fold CV with per-fold normalization
- FiLM conditioning for temporal features
- Full visualization and per-date metrics (GLM-parity)
- Memory-constrained search space for GPU instances
- Intermediate checkpointing
- Reproducible DataLoader worker seeding
"""

import os
import sys
import io
import re
import gc
import json
import pickle
import argparse
import warnings
import tempfile
import random
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, asdict

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler

import optuna
from optuna.trial import TrialState

from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, average_precision_score, cohen_kappa_score,
    matthews_corrcoef, precision_recall_curve, confusion_matrix
)

# Rasterio for GDAL-free raster I/O
import rasterio
from rasterio.env import Env
from rasterio.session import AWSSession

# Visualization
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import boto3


def str2bool(v: str) -> bool:
    return str(v).lower() in ("1", "true", "t", "yes", "y", "on")


# ======================== Speed setup ========================

def speed_setup():
    """Enable fast kernels and lower-precision matmul where supported."""
    try:
        torch.backends.cudnn.benchmark = True
        torch.backends.cuda.matmul.allow_tf32 = True
        torch.backends.cudnn.allow_tf32 = True
        torch.set_float32_matmul_precision("high")
    except Exception:
        pass

speed_setup()


# ======================== Configuration ========================

@dataclass
class UNetConfig:
    # Architecture
    in_channels: Optional[int] = None
    n_temporal: Optional[int] = None
    base_channels: int = 32
    depth: int = 4
    use_film: bool = True

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
    stride_test: int = 128
    min_valid_fraction: float = 0.2
    
    # Augmentation
    use_augmentation: bool = True
    augment_flip: bool = True
    augment_rotate: bool = True

    # Class balancing
    use_weighted_sampler: bool = True
    sampler_pos_floor: float = 0.02
    use_bce_pos_weight: bool = False
    bce_pos_weight: float = 3.0

    # CV
    n_folds: int = 5
    train_frac: float = 0.6
    val_frac: float = 0.2
    force_test_ws: Optional[str] = None

    # System
    num_workers: int = 4
    device: Optional[str] = None
    seed: int = 42
    use_amp: bool = True
    channels_last: bool = True
    
    # Checkpointing
    checkpoint_dir: Optional[str] = None

    def __post_init__(self):
        if self.device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"


# ======================== S3 helpers ========================

def _parse_s3_uri(uri: str) -> Tuple[str, str]:
    assert uri.startswith("s3://"), f"Not an s3 uri: {uri}"
    p = uri[5:]
    parts = p.split("/", 1)
    return parts[0], parts[1] if len(parts) > 1 else ""


def _fmt4(x: float) -> str:
    return f"{float(x):.4f}"


def _parse_latlon_id(ws_id: str) -> Tuple[str, str]:
    """Parse watershed ID to get lat/lon strings."""
    m = re.search(r"Lat[_]?([-\d\.]+)_?Lon[_]?([-\d\.]+)", ws_id)
    if not m:
        raise ValueError(f"Cannot parse watershed id: {ws_id}")
    return _fmt4(float(m.group(1))), _fmt4(float(m.group(2)))


def _s3_exists(s3_uri: str) -> bool:
    """Check if S3 object or prefix exists (handles both files and folders)."""
    s3 = boto3.client("s3")
    b, k = _parse_s3_uri(s3_uri)
    try:
        # Try head_object for exact key match
        s3.head_object(Bucket=b, Key=k)
        return True
    except Exception:
        # Fallback to prefix listing (for folders or partial keys)
        try:
            resp = s3.list_objects_v2(Bucket=b, Prefix=k, MaxKeys=1)
            return 'Contents' in resp
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


# ======================== Rasterio-based I/O with S3 auth ========================

def load_watershed_list(training_folder: str) -> List[str]:
    csv_path = f"{training_folder}/list_of_training_locations.csv"
    if not _s3_exists(csv_path):
        raise FileNotFoundError(f"Missing list_of_training_locations.csv at {csv_path}")
    df = _s3_read_csv(csv_path, dtype={"lat_str": str, "lon_str": str})
    ws_ids = []
    for _, row in df[["lat_str", "lon_str"]].dropna().iterrows():
        lat4 = _fmt4(float(row["lat_str"]))
        lon4 = _fmt4(float(row["lon_str"]))
        ws_ids.append(f"Lat_{lat4}_Lon_{lon4}")
    return sorted(list(set(ws_ids)))


def _read_static_stack(static_tif_s3: str) -> Tuple[np.ndarray, List[str], Tuple[int, int], Any, Any]:
    """Read static multiband raster using rasterio with S3 authentication."""
    session = AWSSession(boto3.Session())
    with Env(session=session, GDAL_DISABLE_READDIR_ON_OPEN='YES'):
        with rasterio.open(static_tif_s3) as src:
            H, W = src.height, src.width
            C = src.count
            stack = np.zeros((C, H, W), dtype=np.float32)
            names = []
            
            for i in range(1, C + 1):
                arr = src.read(i, out_dtype=np.float32)
                nd = src.nodatavals[i - 1]
                if nd is not None:
                    arr[arr == float(nd)] = np.nan
                # Map common sentinel values to NaN
                arr[arr == -9999.0] = np.nan
                arr[arr == -32768.0] = np.nan
                stack[i - 1] = arr
                names.append(src.descriptions[i - 1] or f"static_band_{i:02d}")
            
            gt = src.transform
            prj = src.crs.to_wkt() if src.crs else None
            
    return stack, names, (H, W), gt, prj


def _align_static_stack_to_ref(arr: np.ndarray, names: List[str], ref_names: List[str]) -> Tuple[np.ndarray, List[str]]:
    """Align static stack to reference names."""
    C, H, W = arr.shape
    name_to_idx = {n: i for i, n in enumerate(names)}
    aligned = np.full((len(ref_names), H, W), np.nan, dtype=np.float32)
    for i, ref in enumerate(ref_names):
        if ref in name_to_idx:
            aligned[i] = arr[name_to_idx[ref]]
    return aligned, ref_names


def _read_flood_timeseries(flood_tif_s3: str) -> Tuple[Dict[pd.Timestamp, np.ndarray], List[pd.Timestamp], Tuple[int, int], Any, Any]:
    """Read flood label time series using rasterio with S3 authentication."""
    session = AWSSession(boto3.Session())
    with Env(session=session, GDAL_DISABLE_READDIR_ON_OPEN='YES'):
        with rasterio.open(flood_tif_s3) as src:
            H, W = src.height, src.width
            B = src.count
            labels_by_date = {}
            dates = []
            
            for b in range(1, B + 1):
                desc = src.descriptions[b - 1] or ""
                m = re.search(r"\d{4}-\d{2}-\d{2}", desc)
                if m:
                    d = pd.to_datetime(m.group(0)).normalize()
                    lab = src.read(b, out_dtype=np.uint8)
                    
                    # Honor per-band nodata (usually 255 for flood products)
                    nd = src.nodatavals[b - 1]
                    if nd is not None:
                        lab[lab == int(nd)] = 255
                    
                    labels_by_date[d] = lab
                    dates.append(d)
            
            dates = sorted(list(set(dates)))
            gt = src.transform
            prj = src.crs.to_wkt() if src.crs else None
            
    return labels_by_date, dates, (H, W), gt, prj


def _read_timeseries(training_folder: str, lat4: str, lon4: str) -> Dict[pd.Timestamp, Dict[str, float]]:
    """Read temporal features from CSVs."""
    sub_csv = f"{training_folder}/streamflowPredictions/subwatershed_predictions_Lat{lat4}_Lon{lon4}.csv"
    riv_csv = f"{training_folder}/streamflowPredictions/river_basin_predictions_Lat{lat4}_Lon{lon4}.csv"
    ts_by_date = {}

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


def _discover_all_dynamic_channels(training_folder: str, watersheds: List[str]) -> List[str]:
    """Discover union of all dynamic channels across watersheds."""
    all_keys = set()
    for ws_id in watersheds:
        lat4, lon4 = _parse_latlon_id(ws_id)
        ts_by_date = _read_timeseries(training_folder, lat4, lon4)
        for _, rec in ts_by_date.items():
            all_keys.update(rec.keys())
    
    base = [
        "subwatershed_streamflow_pred", "subwatershed_t2m", "subwatershed_ssrd", "subwatershed_tp",
        "river_basin_streamflow_pred", "river_basin_t2m", "river_basin_ssrd", "river_basin_tp"
    ]
    
    extras = sorted(list(all_keys - set(base)))
    return base + extras


# ======================== Feature statistics ========================

def compute_feature_statistics(training_folder: str, watersheds: List[str]) -> Dict:
    """Compute normalization statistics using only training watersheds."""
    static_names_ref: Optional[List[str]] = None
    C_static: Optional[int] = None
    s_stat = ss_stat = cnt_stat = None

    dynamic_channels = _discover_all_dynamic_channels(training_folder, watersheds)
    C_dyn = len(dynamic_channels)
    s_dyn = np.zeros(C_dyn, dtype=np.float64)
    ss_dyn = np.zeros(C_dyn, dtype=np.float64)
    cnt_dyn = np.zeros(C_dyn, dtype=np.float64)

    for ws_id in watersheds:
        lat4, lon4 = _parse_latlon_id(ws_id)
        static_tif = f"{training_folder}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"

        if not _s3_exists(static_tif):
            continue

        arr, names, _, _, _ = _read_static_stack(static_tif)

        if static_names_ref is None:
            static_names_ref = names[:]
            C_static = len(static_names_ref)
            s_stat = np.zeros(C_static, dtype=np.float64)
            ss_stat = np.zeros(C_static, dtype=np.float64)
            cnt_stat = np.zeros(C_static, dtype=np.float64)

        arr_aligned, _ = _align_static_stack_to_ref(arr, names, static_names_ref)
        A = arr_aligned.reshape(len(static_names_ref), -1)
        mask = np.isfinite(A)
        s_stat += np.where(mask, A, 0.0).sum(axis=1)
        ss_stat += np.where(mask, A * A, 0.0).sum(axis=1)
        cnt_stat += mask.sum(axis=1)

        ts_by_date = _read_timeseries(training_folder, lat4, lon4)
        for _, rec in ts_by_date.items():
            vals = np.array([rec.get(ch, np.nan) for ch in dynamic_channels], dtype=np.float64)
            m = ~np.isnan(vals)
            s_dyn += np.where(m, vals, 0.0)
            ss_dyn += np.where(m, vals * vals, 0.0)
            cnt_dyn += m.astype(np.float64)

        del arr, arr_aligned, A, ts_by_date
        gc.collect()

    if static_names_ref is None:
        raise RuntimeError("No valid static rasters found")

    sm_stat = s_stat / np.maximum(cnt_stat, 1.0)
    var_stat = (ss_stat / np.maximum(cnt_stat, 1.0)) - sm_stat ** 2
    var_stat[var_stat < 0] = 0
    sd_stat = np.sqrt(var_stat)
    sd_stat[sd_stat == 0] = 1.0

    sm_dyn = s_dyn / np.maximum(cnt_dyn, 1.0)
    var_dyn = (ss_dyn / np.maximum(cnt_dyn, 1.0)) - sm_dyn ** 2
    var_dyn[var_dyn < 0] = 0
    sd_dyn = np.sqrt(var_dyn)
    sd_dyn[sd_dyn == 0] = 1.0
    
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


# ======================== FiLM Conditioning ========================

class FiLMBlock(nn.Module):
    """Feature-wise Linear Modulation for conditioning on temporal features."""
    def __init__(self, n_temporal: int, n_channels: int):
        super().__init__()
        hidden = max(64, n_channels // 2)
        self.mlp = nn.Sequential(
            nn.Linear(n_temporal, hidden),
            nn.LayerNorm(hidden),
            nn.ReLU(inplace=True),
            nn.Linear(hidden, n_channels * 2)
        )
        
    def forward(self, features: torch.Tensor, temporal: torch.Tensor) -> torch.Tensor:
        params = self.mlp(temporal)
        gamma, beta = params.chunk(2, dim=1)
        gamma = gamma.unsqueeze(-1).unsqueeze(-1)
        beta = beta.unsqueeze(-1).unsqueeze(-1)
        return features * (1 + gamma) + beta


# ======================== U-Net with FiLM ========================

class ConvBlock(nn.Module):
    def __init__(self, in_ch: int, out_ch: int, n_temporal: int = 0, use_film: bool = False):
        super().__init__()
        self.conv1 = nn.Conv2d(in_ch, out_ch, 3, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_ch)
        self.conv2 = nn.Conv2d(out_ch, out_ch, 3, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_ch)
        self.relu = nn.ReLU(inplace=True)
        
        self.film = FiLMBlock(n_temporal, out_ch) if use_film and n_temporal > 0 else None

    def forward(self, x, temporal=None):
        out = self.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        
        if self.film is not None and temporal is not None:
            out = self.film(out, temporal)
        
        return self.relu(out)


class UNetWithFiLM(nn.Module):
    def __init__(self, n_static: int, n_temporal: int, base_ch: int = 32, 
                 depth: int = 4, use_film: bool = True):
        super().__init__()
        self.depth = depth
        self.n_temporal = n_temporal
        self.use_film = use_film
        
        # Encoder
        self.enc = nn.ModuleList()
        self.pool = nn.MaxPool2d(2)
        ch = base_ch
        
        self.enc.append(ConvBlock(n_static, ch, n_temporal if use_film else 0, use_film))
        
        for i in range(depth - 1):
            next_ch = ch * 2
            self.enc.append(ConvBlock(ch, next_ch, n_temporal if use_film else 0, use_film))
            ch = next_ch
        
        # Bottleneck
        self.bottleneck = ConvBlock(ch, ch * 2, n_temporal if use_film else 0, use_film)
        ch *= 2
        
        # Decoder
        self.up = nn.ModuleList()
        self.dec = nn.ModuleList()
        for i in range(depth):
            self.up.append(nn.ConvTranspose2d(ch, ch // 2, 2, stride=2))
            self.dec.append(ConvBlock(ch, ch // 2, n_temporal if use_film else 0, use_film))
            ch //= 2
        
        self.outc = nn.Conv2d(ch, 1, 1)

    def forward(self, x_static, x_temporal=None):
        x = x_static
        skips = []
        
        for i, block in enumerate(self.enc):
            x = block(x, x_temporal)
            skips.append(x)
            if i < self.depth - 1:
                x = self.pool(x)
        
        x = self.bottleneck(x, x_temporal)
        
        for i in range(self.depth):
            x = self.up[i](x)
            skip = skips[-(i + 1)]
            
            if x.shape[-2:] != skip.shape[-2:]:
                dy = skip.size(2) - x.size(2)
                dx = skip.size(3) - x.size(3)
                x = F.pad(x, [dx // 2, dx - dx // 2, dy // 2, dy - dy // 2])
            
            x = torch.cat([x, skip], dim=1)
            x = self.dec[i](x, x_temporal)
        
        return self.outc(x)


# ======================== Dataset ========================

def _seed_worker(worker_id):
    """Seed worker for reproducibility."""
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)


class SpatioTemporalTileDataset(Dataset):
    """Dataset that separates static and temporal features for FiLM conditioning."""
    def __init__(self, training_folder: str, watersheds: List[str],
                 feature_stats: Dict, config: UNetConfig, mode: str = "train"):
        super().__init__()
        self.config = config
        self.mode = mode
        self.stride = config.stride_train if mode == "train" else config.stride_test
        
        self.static_means = feature_stats["static_means"]
        self.static_stds = feature_stats["static_stds"]
        self.dyn_channels = feature_stats["dynamic_channels"]
        self.dynamic_means = feature_stats["dynamic_means"]
        self.dynamic_stds = feature_stats["dynamic_stds"]
        
        self.ws_statics = []
        self.ws_dates = []
        self.ws_dyn = []
        self.ws_labels = []
        self.index = []
        
        for ws_id in watersheds:
            lat4, lon4 = _parse_latlon_id(ws_id)
            static_tif = f"{training_folder}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"
            flood_tif = f"{training_folder}/floodMaps/dswx_s1_timeseries_subwatershed_Lon{lon4}_Lat{lat4}.tif"
            
            if not _s3_exists(static_tif) or not _s3_exists(flood_tif):
                continue
            
            static_stack, names, (H, W), _, _ = _read_static_stack(static_tif)
            ref_names = feature_stats["static_names"]
            static_stack, _ = _align_static_stack_to_ref(static_stack, names, ref_names)
            
            stat = static_stack.astype(np.float32, copy=True)
            for c in range(stat.shape[0]):
                m, s = self.static_means[c], self.static_stds[c]
                band = stat[c]
                band = np.where(np.isfinite(band), band, m)
                stat[c] = (band - m) / s
            
            labels_by_date, dates, (Hf, Wf), _, _ = _read_flood_timeseries(flood_tif)
            if (H, W) != (Hf, Wf):
                continue
            
            ts_by_date = _read_timeseries(training_folder, lat4, lon4)
            if not ts_by_date:
                continue
            
            dyn_list = []
            lab_list = []
            date_list = []
            
            for d in dates:
                feats = ts_by_date.get(d)
                if not feats:
                    continue
                
                vals = np.array([feats.get(ch, np.nan) for ch in self.dyn_channels], dtype=np.float32)
                vals = (np.where(np.isnan(vals), self.dynamic_means, vals) - self.dynamic_means) / self.dynamic_stds
                dyn_list.append(vals)
                
                Y = labels_by_date[d].astype(np.uint8)
                M = (Y != 255).astype(np.uint8)
                lab_list.append({"Y": Y, "M": M})
                date_list.append(d)
            
            if not dyn_list:
                continue
            
            ws_idx = len(self.ws_statics)
            self.ws_statics.append(stat)
            self.ws_dates.append(date_list)
            self.ws_dyn.append(dyn_list)
            self.ws_labels.append(lab_list)
            
            for di, lab in enumerate(lab_list):
                M = lab["M"]
                Y = lab["Y"]
                for y0 in range(0, max(H - config.tile_size + 1, 1), self.stride):
                    for x0 in range(0, max(W - config.tile_size + 1, 1), self.stride):
                        # Store start/end as integers, not slice objects
                        y_start, y_end = y0, min(y0 + config.tile_size, H)
                        x_start, x_end = x0, min(x0 + config.tile_size, W)
                        
                        # Create temporary slices for validation only
                        rr_temp = slice(y_start, y_end)
                        cc_temp = slice(x_start, x_end)
                        m = M[rr_temp, cc_temp]
                        
                        if m.size == 0:
                            continue
                        valid_frac = m.sum() / (m.size + 1e-9)
                        if valid_frac < config.min_valid_fraction:
                            continue
                        pos_frac = (((Y[rr_temp, cc_temp] > 0) & (m == 1)).sum()) / (m.sum() + 1e-9)
                        
                        # Store integers, NOT slice objects
                        self.index.append((ws_idx, di, y_start, y_end, x_start, x_end, float(pos_frac)))
            
            gc.collect()
    
    # CRITICAL FIX: These methods must be at class level, not inside __init__
    def __len__(self):
        return len(self.index)

    @staticmethod
    def _pad_tile_to_size(x_static: np.ndarray, y: np.ndarray, m: np.ndarray, tile_size: int):
        """Pad tile to fixed size with shape validation."""
        # Validate shapes
        if x_static.ndim != 3:
            raise ValueError(f"x_static must be 3D, got shape {x_static.shape}")
        if y.ndim != 2:
            raise ValueError(f"y must be 2D, got shape {y.shape}")
        if m.ndim != 2:
            raise ValueError(f"m must be 2D, got shape {m.shape}")
        
        C, h, w = x_static.shape
        th, tw = tile_size, tile_size

        pad_h = max(0, th - h)
        pad_w = max(0, tw - w)

        if pad_h > 0 or pad_w > 0:
            x_static = np.pad(
                x_static,
                pad_width=((0, 0), (0, pad_h), (0, pad_w)),
                mode="constant",
                constant_values=0.0
            )
            y = np.pad(y, pad_width=((0, pad_h), (0, pad_w)), mode="constant", constant_values=0.0)
            m = np.pad(m, pad_width=((0, pad_h), (0, pad_w)), mode="constant", constant_values=0.0)

        # Crop if oversized
        if x_static.shape[1] > th or x_static.shape[2] > tw:
            x_static = x_static[:, :th, :tw]
        if y.shape[0] > th or y.shape[1] > tw:
            y = y[:th, :tw]
        if m.shape[0] > th or m.shape[1] > tw:
            m = m[:th, :tw]

        return x_static.astype(np.float32, copy=False), y.astype(np.float32, copy=False), m.astype(np.float32, copy=False)        
    def __getitem__(self, i: int):
        ws_idx, di, y_start, y_end, x_start, x_end, _ = self.index[i]
        
        # Reconstruct slices from integers
        rr = slice(y_start, y_end)
        cc = slice(x_start, x_end)
        
        try:
            stat = self.ws_statics[ws_idx]
            x_static = stat[:, rr, cc].copy()
            
            x_temporal = self.ws_dyn[ws_idx][di].copy()
            
            Y = self.ws_labels[ws_idx][di]["Y"]
            M = self.ws_labels[ws_idx][di]["M"]
            y = (Y[rr, cc] > 0).astype(np.float32)
            m = M[rr, cc].astype(np.float32)
            
            # Validate shapes before padding
            if x_static.ndim != 3:
                raise ValueError(f"x_static has {x_static.ndim} dims, expected 3. Shape: {x_static.shape}")
            if y.ndim != 2:
                raise ValueError(f"y has {y.ndim} dims, expected 2. Shape: {y.shape}")
            if m.ndim != 2:
                raise ValueError(f"m has {m.ndim} dims, expected 2. Shape: {m.shape}")
            
            if (y.shape[0] != self.config.tile_size) or (y.shape[1] != self.config.tile_size):
                x_static, y, m = self._pad_tile_to_size(x_static, y, m, self.config.tile_size)
            
            if self.mode == "train" and self.config.use_augmentation:
                x_static, y, m = self._augment(x_static, y, m)
            
            return (torch.from_numpy(x_static), 
                    torch.from_numpy(x_temporal),
                    torch.from_numpy(y),
                    torch.from_numpy(m))
        
        except Exception as e:
            # Provide detailed error message for debugging
            raise ValueError(
                f"Error in __getitem__ at index {i}: ws_idx={ws_idx}, di={di}, "
                f"y_slice=({y_start},{y_end}), x_slice=({x_start},{x_end}). "
                f"Original error: {type(e).__name__}: {str(e)}"
            )        
    def _augment(self, x_static, y, m):
        """Apply spatial augmentation with shape validation."""
        # Validate input shapes
        assert x_static.ndim == 3, f"x_static must be 3D, got {x_static.shape}"
        assert y.ndim == 2, f"y must be 2D, got {y.shape}"
        assert m.ndim == 2, f"m must be 2D, got {m.shape}"
        
        # Flip (pick axis for 3D static, adjust for 2D labels/mask)
        if self.config.augment_flip and np.random.rand() > 0.5:
            if np.random.rand() > 0.5:
                # Vertical flip
                x_static = np.ascontiguousarray(np.flip(x_static, axis=1))
                y = np.ascontiguousarray(np.flip(y, axis=0))
                m = np.ascontiguousarray(np.flip(m, axis=0))
            else:
                # Horizontal flip
                x_static = np.ascontiguousarray(np.flip(x_static, axis=2))
                y = np.ascontiguousarray(np.flip(y, axis=1))
                m = np.ascontiguousarray(np.flip(m, axis=1))
        
        # Rotation (90 degree only)
        if self.config.augment_rotate and np.random.rand() > 0.5:
            k = np.random.choice([1, 3])  # 90 or 270 degrees
            x_static = np.ascontiguousarray(np.rot90(x_static, k, axes=(1, 2)))
            y = np.ascontiguousarray(np.rot90(y, k, axes=(0, 1)))
            m = np.ascontiguousarray(np.rot90(m, k, axes=(0, 1)))
        
        return x_static, y, m
    
    def make_weighted_sampler(self) -> Optional[WeightedRandomSampler]:
        """Create weighted sampler for class imbalance."""
        if not self.config.use_weighted_sampler or len(self.index) == 0:
            return None
        
        weights = np.array([pf for *_, pf in self.index], dtype=np.float32)
        weights = weights + self.config.sampler_pos_floor
        return WeightedRandomSampler(weights, num_samples=len(weights), replacement=True)


# ======================== Training ========================

def masked_dice_loss(logits, targets, mask, eps=1e-6):
    probs = torch.sigmoid(logits)
    probs = probs * mask
    targets = targets * mask
    inter = (probs * targets).sum(dim=(1, 2))
    denom = probs.sum(dim=(1, 2)) + targets.sum(dim=(1, 2)) + eps
    dice = 2.0 * inter / denom
    return 1.0 - dice.mean()


class Trainer:
    def __init__(self, model, config):
        self.model = model.to(config.device)
        self.config = config
        self.device = config.device

        self.optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=config.learning_rate,
            weight_decay=config.weight_decay
        )

        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer, mode="min", patience=3, factor=0.5
        )

        self.use_amp = config.use_amp and (self.device == "cuda")
        self.scaler = torch.cuda.amp.GradScaler() if self.use_amp else None

        self.best_loss = float("inf")
        self.patience_counter = 0
    
    def train_epoch(self, dataloader):
        self.model.train()
        total_loss = 0.0
        n_batches = 0

        for x_static, x_temporal, y, mask in dataloader:
            if self.config.channels_last and self.device == "cuda":
                x_static = x_static.contiguous(memory_format=torch.channels_last).to(self.device, non_blocking=True)
            else:
                x_static = x_static.to(self.device, non_blocking=True)

            x_temporal = x_temporal.to(self.device, non_blocking=True)
            y = y.to(self.device, non_blocking=True)
            mask = mask.to(self.device, non_blocking=True)

            self.optimizer.zero_grad(set_to_none=True)

            if self.use_amp:
                with torch.cuda.amp.autocast():
                    logits = self.model(x_static, x_temporal).squeeze(1)
                    
                    if self.config.use_bce_pos_weight:
                        pos_weight = torch.tensor(self.config.bce_pos_weight, device=self.device)
                        bce = F.binary_cross_entropy_with_logits(
                            logits, y, reduction="none", pos_weight=pos_weight
                        )
                    else:
                        bce = F.binary_cross_entropy_with_logits(logits, y, reduction="none")

                    bce = (bce * mask).sum() / (mask.sum() + 1e-8)
                    dice = masked_dice_loss(logits, y, mask)
                    loss = 0.5 * bce + 0.5 * dice

                self.scaler.scale(loss).backward()
                self.scaler.unscale_(self.optimizer)
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
                self.scaler.step(self.optimizer)
                self.scaler.update()
            else:
                logits = self.model(x_static, x_temporal).squeeze(1)
                
                if self.config.use_bce_pos_weight:
                    pos_weight = torch.tensor(self.config.bce_pos_weight, device=self.device)
                    bce = F.binary_cross_entropy_with_logits(
                        logits, y, reduction="none", pos_weight=pos_weight
                    )
                else:
                    bce = F.binary_cross_entropy_with_logits(logits, y, reduction="none")

                bce = (bce * mask).sum() / (mask.sum() + 1e-8)
                dice = masked_dice_loss(logits, y, mask)
                loss = 0.5 * bce + 0.5 * dice
                
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
                self.optimizer.step()

            total_loss += float(loss.item())
            n_batches += 1

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        return {"loss": total_loss / max(n_batches, 1)}
    
    def validate(self, dataloader):
        self.model.eval()
        total_loss = 0.0
        all_probs = []
        all_targets = []
        n_batches = 0

        with torch.no_grad():
            for x_static, x_temporal, y, mask in dataloader:
                if self.config.channels_last and self.device == "cuda":
                    x_static = x_static.contiguous(memory_format=torch.channels_last).to(self.device, non_blocking=True)
                else:
                    x_static = x_static.to(self.device, non_blocking=True)

                x_temporal = x_temporal.to(self.device, non_blocking=True)
                y = y.to(self.device, non_blocking=True)
                mask = mask.to(self.device, non_blocking=True)

                if self.use_amp:
                    with torch.cuda.amp.autocast():
                        logits = self.model(x_static, x_temporal).squeeze(1)
                        per_pixel = F.binary_cross_entropy_with_logits(logits, y, reduction="none")
                        loss = (per_pixel * mask).sum() / (mask.sum() + 1e-8)
                else:
                    logits = self.model(x_static, x_temporal).squeeze(1)
                    per_pixel = F.binary_cross_entropy_with_logits(logits, y, reduction="none")
                    loss = (per_pixel * mask).sum() / (mask.sum() + 1e-8)

                total_loss += float(loss.item())
                n_batches += 1

                probs = torch.sigmoid(logits)
                valid = mask > 0
                if valid.any():
                    all_probs.extend(probs[valid].cpu().numpy())
                    all_targets.extend(y[valid].cpu().numpy())

        if all_probs and all_targets:
            all_probs = np.array(all_probs)
            all_targets = np.array(all_targets)

            # Use PR-AUC (consistent with Optuna objective)
            try:
                pr_auc = average_precision_score(all_targets, all_probs)
            except:
                pr_auc = 0.0
            
            # Also find best threshold
            best_f1 = 0
            best_thresh = 0.5
            for t in np.linspace(0.2, 0.8, 13):
                preds = (all_probs >= t).astype(int)
                f1 = f1_score(all_targets, preds, zero_division=0)
                if f1 > best_f1:
                    best_f1 = f1
                    best_thresh = t
        else:
            best_thresh = 0.5
            best_f1 = 0.0
            pr_auc = 0.0

        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        return {
            "loss": total_loss / max(n_batches, 1),
            "f1": best_f1,
            "pr_auc": pr_auc,
            "threshold": best_thresh
        }
    
    def save_checkpoint(self, epoch: int, path: str):
        """Save training checkpoint."""
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'scheduler_state_dict': self.scheduler.state_dict(),
            'best_loss': self.best_loss,
            'patience_counter': self.patience_counter
        }
        if self.scaler:
            checkpoint['scaler_state_dict'] = self.scaler.state_dict()
        torch.save(checkpoint, path)
    
    def load_checkpoint(self, path: str):
        """Load training checkpoint."""
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        self.scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
        self.best_loss = checkpoint['best_loss']
        self.patience_counter = checkpoint['patience_counter']
        if self.scaler and 'scaler_state_dict' in checkpoint:
            self.scaler.load_state_dict(checkpoint['scaler_state_dict'])
        return checkpoint['epoch']


# ======================== K-fold CV ========================

def spatial_kfold_split(watersheds: List[str], config: UNetConfig):
    """Create k-fold spatial splits."""
    rng = np.random.RandomState(config.seed)
    
    if config.force_test_ws and config.force_test_ws in watersheds:
        test_ws = [config.force_test_ws]
        remaining = [w for w in watersheds if w != config.force_test_ws]
    else:
        n_test = max(1, int(len(watersheds) * 0.2))
        test_ws = rng.choice(watersheds, n_test, replace=False).tolist()
        remaining = [w for w in watersheds if w not in test_ws]
    
    folds = []
    for fold_idx in range(config.n_folds):
        rng_fold = np.random.RandomState(config.seed + fold_idx)
        shuffled = remaining[:]
        rng_fold.shuffle(shuffled)
        
        n_train = max(1, int(config.train_frac * len(shuffled)))
        n_val = max(1, int(config.val_frac * len(shuffled)))
        
        train_ws = shuffled[:n_train]
        val_ws = shuffled[n_train:n_train + n_val]
        
        folds.append({'train': train_ws, 'val': val_ws, 'test': test_ws})
    
    return folds


# ======================== Main Training ========================

def train_unet_flood_model(training_folder: str, config: Optional[UNetConfig] = None,
                           save_path: Optional[str] = None) -> Dict:
    """Train U-Net+FiLM flood model with k-fold CV."""
    
    config = config or UNetConfig()
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    random.seed(config.seed)
    
    print("=" * 70)
    print("U-NET + FiLM FLOOD MODEL - K-FOLD SPATIAL CV")
    print("=" * 70)
    
    watersheds = load_watershed_list(training_folder)
    if len(watersheds) < 3:
        raise ValueError(f"Need at least 3 watersheds, found {len(watersheds)}")
    
    folds = spatial_kfold_split(watersheds, config)
    print(f"Test watersheds: {folds[0]['test']}")
    
    cv_thresholds = []
    
    # Create worker seed generator
    g = torch.Generator()
    g.manual_seed(config.seed)
    
    for fold_idx, fold in enumerate(folds):
        print(f"\nFold {fold_idx + 1}/{len(folds)}")
        print(f"  Train: {len(fold['train'])} watersheds")
        print(f"  Val: {len(fold['val'])} watersheds")
        
        feature_stats = compute_feature_statistics(training_folder, fold['train'])
        
        train_dataset = SpatioTemporalTileDataset(
            training_folder, fold['train'], feature_stats, config, mode="train"
        )
        val_dataset = SpatioTemporalTileDataset(
            training_folder, fold['val'], feature_stats, config, mode="val"
        )
        
        if len(train_dataset) == 0 or len(val_dataset) == 0:
            print("  Warning: empty dataset")
            cv_thresholds.append(0.5)
            continue
        
        num_train_workers = min(8, os.cpu_count() or 6)  # More workers

        sampler = train_dataset.make_weighted_sampler()
        train_loader = DataLoader(
            train_dataset,
            batch_size=config.batch_size,
            shuffle=(sampler is None),
            sampler=sampler,
            num_workers=num_train_workers,  # <-- Increased from config.num_workers
            pin_memory=(config.device == "cuda"),
            persistent_workers=(num_train_workers > 0),
            prefetch_factor=4,  # <-- ADDED: Prefetch more batches
            worker_init_fn=_seed_worker,
            generator=g
        )

        num_val_workers = min(6, os.cpu_count() or 4)
        val_loader = DataLoader(
            val_dataset,
            batch_size=config.batch_size * 2,  # Larger val batches
            shuffle=False,
            num_workers=min(6, os.cpu_count() or 4),  # <-- More workers
            pin_memory=(config.device == "cuda"),
            persistent_workers=(num_val_workers > 0),
            prefetch_factor=4  # <-- ADDED
        )
        
        model = UNetWithFiLM(
            n_static=feature_stats["n_static"],
            n_temporal=feature_stats["n_temporal"],
            base_ch=config.base_channels,
            depth=config.depth,
            use_film=config.use_film
        )
        
        if config.channels_last and config.device == "cuda":
            model = model.to(memory_format=torch.channels_last)
        
        trainer = Trainer(model, config)
        best_threshold = 0.5
        best_val_loss = float("inf")
        
        for epoch in range(min(config.epochs, 15)):
            train_metrics = trainer.train_epoch(train_loader)
            val_metrics = trainer.validate(val_loader)
            trainer.scheduler.step(val_metrics["loss"])
            
            print(f"  Epoch {epoch+1}: Train loss={train_metrics['loss']:.4f}, "
                  f"Val loss={val_metrics['loss']:.4f}, Val F1={val_metrics['f1']:.4f}, "
                  f"Val PR-AUC={val_metrics['pr_auc']:.4f}")
            
            # Log GPU memory
            if config.device == "cuda" and torch.cuda.is_available():
                allocated = torch.cuda.memory_allocated(0) / 1e9
                reserved = torch.cuda.memory_reserved(0) / 1e9
                if (epoch + 1) % 5 == 0:
                    print(f"    GPU mem: {allocated:.2f}GB alloc, {reserved:.2f}GB reserved")
            
            if val_metrics["loss"] < best_val_loss - config.min_delta:
                best_val_loss = val_metrics["loss"]
                best_threshold = val_metrics["threshold"]
                trainer.patience_counter = 0
            else:
                trainer.patience_counter += 1
                if trainer.patience_counter >= config.patience:
                    print(f"  Early stopping at epoch {epoch+1}")
                    break
        
        cv_thresholds.append(best_threshold)
        
        del model, trainer
        gc.collect()
        torch.cuda.empty_cache()
    
    # Final model
    print("\n" + "=" * 70)
    print("FINAL MODEL TRAINING")
    print("=" * 70)
    
    final_threshold = float(np.mean(cv_thresholds)) if cv_thresholds else 0.5
    print(f"Average threshold from CV: {final_threshold:.3f}")
    
    test_ws = folds[0]['test']
    train_ws = [w for w in watersheds if w not in test_ws]
    
    final_stats = compute_feature_statistics(training_folder, train_ws)
    
    final_dataset = SpatioTemporalTileDataset(
        training_folder, train_ws, final_stats, config, mode="train"
    )
    
    sampler = final_dataset.make_weighted_sampler()
    final_loader = DataLoader(
        final_dataset,
        batch_size=config.batch_size,
        shuffle=(sampler is None),
        sampler=sampler,
        num_workers=config.num_workers,
        pin_memory=(config.device == "cuda"),
        persistent_workers=(config.num_workers > 0),
        worker_init_fn=_seed_worker,
        generator=g
    )
    
    final_model = UNetWithFiLM(
        n_static=final_stats["n_static"],
        n_temporal=final_stats["n_temporal"],
        base_ch=config.base_channels,
        depth=config.depth,
        use_film=config.use_film
    )
    
    if config.channels_last and config.device == "cuda":
        final_model = final_model.to(memory_format=torch.channels_last)
    
    final_trainer = Trainer(final_model, config)
    
    for epoch in range(config.epochs):
        train_metrics = final_trainer.train_epoch(final_loader)
        print(f"Epoch {epoch+1}: Loss={train_metrics['loss']:.4f}")
        
        # Save checkpoint every 5 epochs
        if config.checkpoint_dir and (epoch + 1) % 5 == 0:
            os.makedirs(config.checkpoint_dir, exist_ok=True)
            ckpt_path = os.path.join(config.checkpoint_dir, f"checkpoint_epoch_{epoch+1}.pth")
            final_trainer.save_checkpoint(epoch, ckpt_path)
            print(f"  Saved checkpoint: {ckpt_path}")
    
    # Package results
    results = {
        'model': final_model,
        'threshold': final_threshold,
        'feature_stats': final_stats,
        'config': config,
        'test_watersheds': test_ws,
        'training_folder': training_folder
    }
    
    if save_path:
        torch.save(results, save_path)
        print(f"Model saved to {save_path}")
    
    return results


# ======================== Optuna Objective ========================

def objective(trial, training_folder, base_config):
    """Optuna objective for U-Net hyperparameter tuning (PR-AUC objective)."""
    
    # Memory check
    try:
        import psutil
        if psutil.virtual_memory().percent > 85:
            gc.collect()
            torch.cuda.empty_cache()
            if psutil.virtual_memory().percent > 90:
                raise optuna.TrialPruned("Memory usage too high")
    except Exception:
        pass
    
    # Explicit batch_size suggestion first, then derive tile_size
    bs = trial.suggest_categorical('batch_size', [4, 8, 16])
    ts = 256 if bs <= 8 else 128  # Use 256x256 for batch sizes up to 8
    
    # Constrained search space for 24GB GPU (ml.g5.2xlarge)
    config = UNetConfig(
        # Architecture - constrained to avoid OOM
        base_channels=trial.suggest_categorical('base_channels', [16, 32]),
        depth=trial.suggest_int('depth', 3, 4),
        use_film=trial.suggest_categorical('use_film', [True, False]),
        
        # Training
        learning_rate=trial.suggest_float('learning_rate', 1e-4, 1e-2, log=True),
        weight_decay=trial.suggest_float('weight_decay', 1e-5, 1e-3, log=True),
        batch_size=bs,
        
        # Data
        tile_size=ts,
        stride_train=trial.suggest_int('stride_train', 64, 128, step=32),
        use_augmentation=trial.suggest_categorical('use_augmentation', [True, False]),
        
        # Class balancing
        use_weighted_sampler=trial.suggest_categorical('use_weighted_sampler', [True, False]),
        use_bce_pos_weight=trial.suggest_categorical('use_bce_pos_weight', [True, False]),
        bce_pos_weight=trial.suggest_float('bce_pos_weight', 2.0, 10.0),
        
        # Fixed
        epochs=base_config.epochs,
        patience=base_config.patience,
        n_folds=2,
        num_workers=0,  # Disable workers during tuning to save memory
        device=base_config.device,
        seed=base_config.seed + trial.number
    )
    
    watersheds = load_watershed_list(training_folder)
    
    # Use subset for tuning
    if len(watersheds) > 10:
        np.random.seed(config.seed)
        watersheds = np.random.choice(watersheds, 10, replace=False).tolist()
    
    # Simple train/val split
    np.random.seed(config.seed)
    np.random.shuffle(watersheds)
    n_train = int(0.7 * len(watersheds))
    train_ws = watersheds[:n_train]
    val_ws = watersheds[n_train:]
    
    # Compute stats
    feature_stats = compute_feature_statistics(training_folder, train_ws)
    
    # Create datasets
    train_dataset = SpatioTemporalTileDataset(
        training_folder, train_ws, feature_stats, config, mode="train"
    )
    val_dataset = SpatioTemporalTileDataset(
        training_folder, val_ws, feature_stats, config, mode="val"
    )
    
    if len(train_dataset) == 0 or len(val_dataset) == 0:
        return 0.0
    
    # Create dataloaders
    sampler = train_dataset.make_weighted_sampler()
    train_loader = DataLoader(
        train_dataset,
        batch_size=config.batch_size,
        shuffle=(sampler is None),
        sampler=sampler,
        num_workers=2,  # <-- Changed from 0 to 2
        prefetch_factor=2  # <-- ADDED
    )
    val_loader = DataLoader(
        val_dataset, 
        batch_size=config.batch_size * 2,  # <-- Larger val batches
        shuffle=False, 
        num_workers=2,  # <-- Changed from 0 to 2
        prefetch_factor=2  # <-- ADDED
    )
    
    # Train model
    model = UNetWithFiLM(
        n_static=feature_stats["n_static"],
        n_temporal=feature_stats["n_temporal"],
        base_ch=config.base_channels,
        depth=config.depth,
        use_film=config.use_film
    )

    # Apply channels_last optimization during tuning too
    if config.channels_last and config.device == "cuda":
        model = model.to(memory_format=torch.channels_last)

    trainer = Trainer(model, config)    
    # Quick training
    for epoch in range(config.epochs):
        trainer.train_epoch(train_loader)
        val_metrics = trainer.validate(val_loader)
        
        # Report PR-AUC for pruning (consistent with objective)
        trial.report(val_metrics['pr_auc'], epoch)
        
        if trial.should_prune():
            raise optuna.TrialPruned()
    
    # Return final PR-AUC
    pr_auc = val_metrics['pr_auc']
    
    del model, trainer
    torch.cuda.empty_cache()
    
    return pr_auc


# ======================== Visualization (GLM-parity) ========================

class UNetFloodVisualizer:
    """Generate visualizations matching GLM/MLP format."""
    
    def __init__(self, model_results: Dict, training_folder: str):
        self.model_results = model_results
        self.model = model_results['model']
        self.config = model_results['config']
        self.model = self.model.to(self.config.device)
        self.threshold = model_results['threshold']
        self.feature_stats = model_results['feature_stats']
        self.training_folder = training_folder
        
        # Parse watershed info
        test_ws = model_results['test_watersheds'][0]
        lat4, lon4 = _parse_latlon_id(test_ws)
        self.test_lat = lat4
        self.test_lon = lon4
        self.test_watershed = test_ws
        
        self.s3_client = boto3.client('s3')
    
    def generate_all_outputs(self, output_folder: Optional[str] = None):
        """Generate all visualization outputs."""
        
        if output_folder is None:
            output_folder = f"{self.training_folder}/assessUNetModel"
        
        print("\n" + "="*60)
        print("GENERATING U-NET MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        print(f"Test watershed: {self.test_watershed}")
        print(f"Output folder: {output_folder}")
        
        # Generate predictions
        predictions, confusions, per_date_metrics = self._generate_predictions()
        
        # Load reference TIFF
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        
        session = AWSSession(boto3.Session())
        with Env(session=session, GDAL_DISABLE_READDIR_ON_OPEN='YES'):
            with rasterio.open(flood_tiff_path) as ref_src:
                dates = sorted(predictions.keys())
                
                # Create outputs
                pred_tiff = f"{output_folder}/unet_predictions_{self.test_watershed}.tif"
                cm_tiff = f"{output_folder}/unet_confusion_matrix_{self.test_watershed}.tif"
                
                self._write_multiband_tiff(ref_src, predictions, dates, pred_tiff, "UNet_Prediction")
                self._write_multiband_tiff(ref_src, confusions, dates, cm_tiff, "UNet_CM")
        
        # Create PNGs and CSVs
        png_paths = self._create_pngs(confusions, dates, per_date_metrics, output_folder)
        csv_paths = self._write_csvs(output_folder, per_date_metrics)
        
        print("\n✅ U-Net visualization complete!")
        
        return {
            'prediction_tiff': pred_tiff,
            'confusion_matrix_tiff': cm_tiff,
            'png_visualizations': png_paths,
            'csv_outputs': csv_paths
        }
    
    def _generate_predictions(self):
        """Generate predictions for test watershed."""
        lat4, lon4 = self.test_lat, self.test_lon
        static_tif = f"{self.training_folder}/rasters/multiband_subwatershed_Lon{lon4}_Lat{lat4}.tif"
        flood_tif = f"{self.training_folder}/floodMaps/dswx_s1_timeseries_subwatershed_Lon{lon4}_Lat{lat4}.tif"
        
        # Load data
        static_stack, names, (H, W), _, _ = _read_static_stack(static_tif)
        ref_names = self.feature_stats["static_names"]
        static_stack, _ = _align_static_stack_to_ref(static_stack, names, ref_names)
        
        # Normalize static
        stat = static_stack.astype(np.float32, copy=True)
        for c in range(stat.shape[0]):
            m, s = self.feature_stats["static_means"][c], self.feature_stats["static_stds"][c]
            band = stat[c]
            band = np.where(np.isfinite(band), band, m)
            stat[c] = (band - m) / s
        
        labels_by_date, dates, _, _, _ = _read_flood_timeseries(flood_tif)
        ts_by_date = _read_timeseries(self.training_folder, lat4, lon4)
        
        predictions = {}
        confusions = {}
        per_date_metrics = {}
        
        self.model.eval()
        
        print(f"  Processing {len(dates)} dates...")
        for date_idx, date_ts in enumerate(dates):
            if (date_idx + 1) % 10 == 0:
                print(f"    Date {date_idx + 1}/{len(dates)}")
            
            feats = ts_by_date.get(date_ts)
            if not feats:
                continue
            
            # Prepare temporal features
            vals = np.array([feats.get(ch, np.nan) for ch in self.feature_stats["dynamic_channels"]], dtype=np.float32)
            vals = (np.where(np.isnan(vals), self.feature_stats["dynamic_means"], vals) - self.feature_stats["dynamic_means"]) / self.feature_stats["dynamic_stds"]
            
            # Predict on tiles
            Y = labels_by_date[date_ts].astype(np.uint8)
            M = (Y != 255).astype(np.uint8)
            
            pred_full = np.zeros((H, W), dtype=np.float32)
            count_full = np.zeros((H, W), dtype=np.float32)
            
            tile_size = self.config.tile_size
            stride = self.config.stride_test
            
            with torch.no_grad():
                for y0 in range(0, max(H - tile_size + 1, 1), stride):
                    for x0 in range(0, max(W - tile_size + 1, 1), stride):
                        rr = slice(y0, min(y0 + tile_size, H))
                        cc = slice(x0, min(x0 + tile_size, W))
                        
                        tile_static = stat[:, rr, cc]
                        tile_h, tile_w = tile_static.shape[1], tile_static.shape[2]
                        
                        # Pad if needed
                        if tile_h < tile_size or tile_w < tile_size:
                            pad_h = tile_size - tile_h
                            pad_w = tile_size - tile_w
                            tile_static = np.pad(tile_static, ((0,0), (0, pad_h), (0, pad_w)))
                        
                        x_static = torch.from_numpy(tile_static).unsqueeze(0).float().to(self.config.device)
                        x_temporal = torch.from_numpy(vals).unsqueeze(0).float().to(self.config.device)
                        
                        logits = self.model(x_static, x_temporal).squeeze()
                        probs = torch.sigmoid(logits).cpu().numpy()
                        
                        # Crop to actual size
                        probs = probs[:tile_h, :tile_w]
                        
                        pred_full[rr, cc] += probs
                        count_full[rr, cc] += 1.0
            
            # Average overlapping predictions
            pred_full = pred_full / np.maximum(count_full, 1.0)
            y_pred = (pred_full >= self.threshold).astype(np.uint8)
            
            # Create images
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            pred_img[M == 1] = y_pred[M == 1]
            
            # Use ISO date string as key for consistency
            iso_date = date_ts.date().isoformat()
            predictions[iso_date] = pred_img
            
            # Confusion matrix
            cm_img = np.full((H, W), 255, dtype=np.uint8)
            y_true = Y[M == 1]
            y_pred_valid = y_pred[M == 1]
            cm_vals = np.where(
                y_true == 0,
                np.where(y_pred_valid == 0, 0, 3),
                np.where(y_pred_valid == 1, 1, 2)
            )
            cm_img[M == 1] = cm_vals
            confusions[iso_date] = cm_img
            
            # Metrics - matching GLM format
            probs_valid = pred_full[M == 1]
            
            try:
                roc_auc = roc_auc_score(y_true, probs_valid)
            except:
                roc_auc = np.nan
            
            try:
                pr_auc = average_precision_score(y_true, probs_valid)
            except:
                pr_auc = np.nan
            
            if len(np.unique(y_true)) > 1 and len(np.unique(y_pred_valid)) > 1:
                mcc = matthews_corrcoef(y_true, y_pred_valid)
            else:
                mcc = 0.0
            
            per_date_metrics[iso_date] = {
                'date': iso_date,
                'F1': f1_score(y_true, y_pred_valid, zero_division=0),
                'Accuracy': accuracy_score(y_true, y_pred_valid),
                'Precision': precision_score(y_true, y_pred_valid, zero_division=0),
                'Recall': recall_score(y_true, y_pred_valid, zero_division=0),
                'ROC_AUC': float(roc_auc) if np.isfinite(roc_auc) else np.nan,
                'PR_AUC': float(pr_auc) if np.isfinite(pr_auc) else np.nan,
                'Kappa': cohen_kappa_score(y_true, y_pred_valid),
                'MCC': float(mcc),
                'actual_pos_frac': float(y_true.mean()),
                'pred_pos_frac': float(y_pred_valid.mean())
            }
        
        return predictions, confusions, per_date_metrics
    
    def _write_multiband_tiff(self, ref_src, arrays, dates, output_path, prefix):
        """Write multiband TIFF using rasterio."""
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
            temp_path = tmp.name
        
        profile = ref_src.profile.copy()
        profile.update({
            'count': len(dates),
            'dtype': 'uint8',
            'compress': 'deflate',
            'tiled': True,
            'nodata': 255
        })
        
        with rasterio.open(temp_path, 'w', **profile) as dst:
            for i, date_str in enumerate(dates, 1):
                dst.write(arrays.get(date_str, np.full((ref_src.height, ref_src.width), 255, dtype=np.uint8)), i)
                dst.set_band_description(i, f"{prefix}_{date_str}")
        
        _s3_upload(temp_path, output_path)
        os.remove(temp_path)
        print(f"  Created: {output_path}")
    
    def _create_pngs(self, confusions, dates, metrics, output_folder):
        """Create PNG visualizations."""
        png_paths = []
        
        for date_str in sorted(dates)[:50]:  # Limit to 50 PNGs
            cm = confusions.get(date_str)
            if cm is None or np.count_nonzero(cm != 255) < 100:
                continue
            
            m = metrics.get(date_str, {})
            
            fig, ax = plt.subplots(figsize=(12, 10))
            cm_disp = cm.astype(float)
            cm_disp[cm == 255] = np.nan
            
            cmap = matplotlib.colors.ListedColormap(['green', 'blue', 'orange', 'red'])
            ax.imshow(cm_disp, cmap=cmap, vmin=0, vmax=3)
            
            legend_elements = [
                mpatches.Patch(color='green', label='True Negative'),
                mpatches.Patch(color='blue', label='True Positive'),
                mpatches.Patch(color='orange', label='False Negative'),
                mpatches.Patch(color='red', label='False Positive')
            ]
            ax.legend(handles=legend_elements, loc='upper right')
            ax.set_title(f'U-Net - {date_str}', fontweight='bold')
            
            text = (
                f"Date: {date_str}\n"
                f"F1: {m.get('F1', 0):.3f}\n"
                f"Accuracy: {m.get('Accuracy', 0):.3f}\n"
                f"Precision: {m.get('Precision', 0):.3f}\n"
                f"Recall: {m.get('Recall', 0):.3f}\n"
                f"ROC AUC: {m.get('ROC_AUC', np.nan):.3f}\n"
                f"PR AUC: {m.get('PR_AUC', np.nan):.3f}"
            )
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            ax.text(0.02, 0.98, text, transform=ax.transAxes, va='top', bbox=props)
            ax.axis('off')
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()
            
            s3_path = f"{output_folder}/unet_confusion_matrix_{date_str}.png"
            _s3_upload(temp_png, s3_path)
            os.remove(temp_png)
            png_paths.append(s3_path)
        
        print(f"  Created {len(png_paths)} PNGs")
        return png_paths
    
    def _write_csvs(self, output_folder, per_date_metrics):
        """Write CSV outputs matching GLM format."""
        paths = {}
        
        if per_date_metrics:
            df = pd.DataFrame.from_dict(per_date_metrics, orient='index')
            column_order = ['date', 'F1', 'Accuracy', 'Precision', 'Recall', 'ROC_AUC', 
                           'PR_AUC', 'Kappa', 'MCC', 'actual_pos_frac', 'pred_pos_frac']
            df = df[column_order]
            
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
            df.to_csv(tmp_csv, index=False)
            s3_path = f"{output_folder}/unet_per_date_metrics_{self.test_watershed}.csv"
            _s3_upload(tmp_csv, s3_path)
            os.remove(tmp_csv)
            paths['per_date_metrics'] = s3_path
            print(f"  Created: {s3_path}")
        
        return paths


# ======================== Main Entry Point ========================

def main():
    parser = argparse.ArgumentParser()
    
    # SageMaker arguments
    parser.add_argument('--model-dir', type=str, default=os.environ.get('SM_MODEL_DIR', '/opt/ml/model'))
    parser.add_argument('--output-dir', type=str, default=os.environ.get('SM_OUTPUT_DATA_DIR', '/opt/ml/output'))
    
    # Training arguments
    parser.add_argument('--training-folder', type=str, required=True)
    parser.add_argument('--n-trials', type=int, default=30)
    parser.add_argument('--epochs', type=int, default=15)
    parser.add_argument('--patience', type=int, default=5)
    parser.add_argument('--use-gpu', type=str2bool, default=True)
    parser.add_argument('--num-workers', type=int, default=6)
    
    args = parser.parse_args()
    
    # Configure device
    device = "cuda" if args.use_gpu and torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    
    if torch.cuda.is_available():
        print(f"GPU: {torch.cuda.get_device_name(0)}")
        props = torch.cuda.get_device_properties(0)
        print(f"Memory: {props.total_memory / 1e9:.2f} GB")
        print(f"Allocated: {torch.cuda.memory_allocated(0) / 1e9:.2f} GB")
    
    # Base configuration
    base_config = UNetConfig(
        epochs=args.epochs,
        patience=args.patience,
        device=device,
        num_workers=args.num_workers,
        checkpoint_dir=os.environ.get('SM_CHECKPOINT_DIR', '/opt/ml/checkpoints')
    )
    
    # Create Optuna study
    study = optuna.create_study(
        direction='maximize',
        sampler=optuna.samplers.TPESampler(seed=42),
        pruner=optuna.pruners.MedianPruner(n_startup_trials=5)
    )
    
    # Run optimization
    print(f"Starting Optuna optimization with {args.n_trials} trials...")
    study.optimize(
        lambda trial: objective(trial, args.training_folder, base_config),
        n_trials=args.n_trials,
        show_progress_bar=True
    )
    
    # Results
    print("\n" + "="*70)
    print("HYPERPARAMETER OPTIMIZATION COMPLETE")
    print("="*70)
    
    best_trial = study.best_trial
    print(f"Best PR-AUC: {best_trial.value:.4f}")
    print("\nBest parameters:")
    for key, value in best_trial.params.items():
        print(f"  {key}: {value}")
    
    # Save results
    results = {
        'best_params': best_trial.params,
        'best_value': best_trial.value,
        'n_trials': len(study.trials),
        'study_stats': {
            'n_pruned': len([t for t in study.trials if t.state == TrialState.PRUNED]),
            'n_complete': len([t for t in study.trials if t.state == TrialState.COMPLETE]),
            'n_failed': len([t for t in study.trials if t.state == TrialState.FAIL])
        }
    }
    
    with open(os.path.join(args.model_dir, 'optuna_results.json'), 'w') as f:
        json.dump(results, f, indent=2)
    
    import joblib
    joblib.dump(study, os.path.join(args.model_dir, 'optuna_study.pkl'))
    
    # Train final model
    print("\n" + "="*70)
    print("TRAINING FINAL MODEL WITH BEST PARAMETERS")
    print("="*70)
    
    best_params = best_trial.params
    
    # Derive tile_size from batch_size
    bs = best_params['batch_size']
    ts = 256 if bs <= 8 else 128  # Use 256x256 for batch sizes up to 8
    
    best_config = UNetConfig(
        base_channels=best_params['base_channels'],
        depth=best_params['depth'],
        use_film=best_params.get('use_film', True),  # <-- FIXED: Use .get() with default
        learning_rate=best_params['learning_rate'],
        weight_decay=best_params['weight_decay'],
        batch_size=bs,
        tile_size=ts,
        stride_train=best_params['stride_train'],
        use_augmentation=best_params['use_augmentation'],
        use_weighted_sampler=best_params['use_weighted_sampler'],
        use_bce_pos_weight=best_params['use_bce_pos_weight'],
        bce_pos_weight=best_params['bce_pos_weight'],
        epochs=base_config.epochs * 2,
        patience=base_config.patience,
        device=device,
        num_workers=args.num_workers,
        n_folds=5,
        checkpoint_dir=base_config.checkpoint_dir
    )
    
    with open(os.path.join(args.model_dir, 'best_config.json'), 'w') as f:
        json.dump(asdict(best_config), f, indent=2)
    
    final_model_results = train_unet_flood_model(
        args.training_folder,
        best_config,
        save_path=os.path.join(args.model_dir, 'unet_final_model.pth')
    )
    
    print("\nTraining complete!")


if __name__ == "__main__" and os.environ.get("RUN_TRAINING", "0") == "1":
    main()

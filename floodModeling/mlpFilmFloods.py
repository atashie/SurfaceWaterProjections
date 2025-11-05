"""
mlpFilmFloods.py - ORIGINAL WITH MINIMAL ADDITIONS

Pixel-wise MLP with FiLM conditioning for flood prediction.

Added features:
- Optional Optuna hyperparameter tuning (n=20 trials)
- Diagnostic plots matching glmFloods.py
"""

import os
import io
import re
import gc
import time
import pickle
import warnings
import tempfile
from typing import List, Tuple, Dict, Optional, Any
from dataclasses import dataclass

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

# NEW: Optuna for hyperparameter tuning
import optuna
from optuna.samplers import TPESampler

# Parquet I/O
import pyarrow as pa
import pyarrow.dataset as ds

# Visualization
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from osgeo import gdal
gdal.UseExceptions()
import boto3

from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, average_precision_score, cohen_kappa_score,
    matthews_corrcoef, precision_recall_curve, confusion_matrix
)


# ======================== Data reading ===============================

def arrow_batch_to_numpy(batch, columns: List[str]) -> np.ndarray:
    """Efficiently convert Arrow batch to NumPy array without pandas."""
    if batch.num_rows == 0:
        return np.empty((0, len(columns)), dtype=np.float32)
    
    arrays = []
    for col in columns:
        arr = batch.column(col).to_numpy(zero_copy_only=False)
        if arr.dtype != np.float32:
            arr = arr.astype(np.float32, copy=False)
        arrays.append(arr)
    
    return np.column_stack(arrays) if arrays else np.empty((batch.num_rows, 0), dtype=np.float32)


# ======================== Configuration ========================

@dataclass
class MLPConfig:
    # Architecture
    hidden_dims: List[int] = None  # Will default to [256, 128, 64]
    dropout: float = 0.3
    use_film: bool = True
    use_interaction_features: bool = True
    n_interaction_features: int = 50

    # Normalization for hidden layers: "layer" (robust) or "batch"
    norm_type: str = "layer"

    # Training
    epochs: int = 30
    batch_size: int = 512
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4

    # Early stopping
    patience: int = 5
    min_delta: float = 0.01

    # Class balancing
    use_weighted_loss: bool = True
    pos_weight: float = 10.0
    focal_loss: bool = True
    focal_alpha: float = 0.25
    focal_gamma: float = 2.0

    # Data
    sample_fraction: float = 0.1
    max_pixels_per_watershed: int = 200000
    max_val_pixels_per_watershed: int = 1000000  # NEW: Cap validation for speed

    # NEW: Hyperparameter tuning
    tune_hyperparameters: bool = False
    n_tuning_iterations: int = 20

    # Fold strategy
    n_folds: int = 5
    train_frac: float = 0.6
    val_frac: float = 0.2
    force_test_ws: Optional[str] = None

    # System
    device: Optional[str] = None
    seed: int = 42
    use_amp: bool = True
    num_workers: int = 0

    # Memory
    batch_size_inference: int = 50000

    # Visualization
    max_png_dates: Optional[int] = None

    def __post_init__(self):
        if self.device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"
        if self.hidden_dims is None:
            self.hidden_dims = [256, 128, 64]

# ======================== GDAL Configuration ========================

def configure_gdal():
    """Configure GDAL/PROJ once."""
    try:
        from pyproj import datadir as _pydatadir
        proj_dir = _pydatadir.get_data_dir()
        if proj_dir and os.path.isdir(proj_dir):
            os.environ["PROJ_LIB"] = proj_dir
            gdal.SetConfigOption("PROJ_LIB", proj_dir)
    except:
        pass
    
    for cand in ["/opt/conda/share/proj", "/usr/share/proj"]:
        if os.path.isdir(cand):
            os.environ.setdefault("PROJ_LIB", cand)
            gdal.SetConfigOption("PROJ_LIB", os.environ["PROJ_LIB"])
            break
    
    gdal.SetConfigOption("PROJ_NETWORK", "ON")
    gdal.SetConfigOption("GDAL_CACHEMAX", "64")
    gdal.SetConfigOption("CPL_DEBUG", "OFF")

configure_gdal()


# ======================== PyArrow Helpers ========================

def get_s3fs():
    try:
        import pyarrow.fs as pafs
        return pafs.S3FileSystem()
    except:
        from pyarrow import filesystem as pafs
        return pafs.S3FileSystem()


def open_parquet_dataset(parquet_uri: str):
    if parquet_uri.startswith("s3://"):
        bkt, key = parquet_uri[5:].split("/", 1)
        s3fs = get_s3fs()
        return ds.dataset(f"{bkt}/{key}", format="parquet", filesystem=s3fs)
    else:
        return ds.dataset(parquet_uri, format="parquet")


def make_scanner(dataset, columns=None, filter=None, batch_size=None):
    kwargs = {}
    if columns is not None:
        kwargs["columns"] = columns
    if filter is not None:
        kwargs["filter"] = filter
    if batch_size is not None:
        try:
            kwargs["batch_size"] = batch_size
        except:
            pass
    return dataset.scanner(**kwargs)


def scanner_batches(scanner):
    if hasattr(scanner, "to_batches"):
        return scanner.to_batches()
    elif hasattr(scanner, "scan_batches"):
        return scanner.scan_batches()
    else:
        tbl = scanner.to_table()
        return tbl.to_batches(max_chunksize=65536)


# ======================== FiLM Module ========================

class FiLMBlock(nn.Module):
    """Feature-wise Linear Modulation for temporal conditioning."""
    def __init__(self, n_temporal: int, n_hidden: int):
        super().__init__()
        hidden_film = max(64, n_hidden)
        self.mlp = nn.Sequential(
            nn.Linear(n_temporal, hidden_film),
            nn.LayerNorm(hidden_film),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_film, n_hidden * 2)
        )
    
    def forward(self, features: torch.Tensor, temporal: torch.Tensor) -> torch.Tensor:
        params = self.mlp(temporal)
        gamma, beta = params.chunk(2, dim=1)
        return features * (1 + torch.tanh(gamma)) + beta

# ======================== MLP + FiLM Model ========================

class MLPWithFiLM(nn.Module):
    """Pixel-wise MLP with FiLM conditioning."""
    def __init__(self, n_static: int, n_temporal: int, config: MLPConfig):
        super().__init__()
        self.n_static = n_static
        self.n_temporal = n_temporal
        self.config = config
        
        if config.use_interaction_features:
            self.n_interactions = min(config.n_interaction_features, n_static * 5)
            self.interaction_weight = nn.Parameter(
                torch.randn(n_static, self.n_interactions) * 0.01
            )
            gate_hidden = max(64, n_temporal)
            self.temporal_gate = nn.Sequential(
                nn.Linear(n_temporal, gate_hidden),
                nn.ReLU(inplace=True),
                nn.Linear(gate_hidden, self.n_interactions),
                nn.Sigmoid()
            )
            input_dim = n_static + self.n_interactions
        else:
            self.n_interactions = 0
            self.temporal_gate = None
            input_dim = n_static
        
        layers = []
        prev_dim = input_dim
        self.hidden_dims = config.hidden_dims[:]
        
        film_layer_indices = []
        if config.use_film and len(self.hidden_dims) >= 1:
            first_idx = 0
            penultimate_idx = max(0, len(self.hidden_dims) - 2)
            film_layer_indices = [first_idx]
            if penultimate_idx != first_idx:
                film_layer_indices.append(penultimate_idx)
        self._film_layer_indices = film_layer_indices
        
        film_added = 0
        
        for i, hidden_dim in enumerate(self.hidden_dims):
            layers.append(nn.Linear(prev_dim, hidden_dim))
            if self.config.norm_type.lower() == "batch":
                layers.append(nn.BatchNorm1d(hidden_dim))
            else:
                layers.append(nn.LayerNorm(hidden_dim))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.Dropout(config.dropout))
            prev_dim = hidden_dim
        
            if config.use_film and (i in film_layer_indices) and (i < len(self.hidden_dims) - 1):
                self.add_module(f'film_{film_added}', FiLMBlock(n_temporal, hidden_dim))
                film_added += 1
        
        layers.append(nn.Linear(prev_dim, 1))
        self.mlp = nn.ModuleList(layers)
        
        self.film_indices = []
        layer_count = 0
        film_counter = 0
        for i in range(len(self.hidden_dims)):
            layer_count += 4
            if config.use_film and (i in film_layer_indices) and (i < len(self.hidden_dims) - 1):
                self.film_indices.append(layer_count - 1)
                film_counter += 1
                
    def forward(self, static_features: torch.Tensor, temporal_features: torch.Tensor) -> torch.Tensor:
        if self.config.use_interaction_features:
            static_proj = torch.matmul(static_features, self.interaction_weight)
            gates = self.temporal_gate(temporal_features)
            interactions = static_proj * gates
            x = torch.cat([static_features, interactions], dim=1)
        else:
            x = static_features

        film_ptr = 0
        for idx, module in enumerate(self.mlp):
            x = module(x)
            if (self.config.use_film
                and film_ptr < len(self.film_indices)
                and idx == self.film_indices[film_ptr]):
                film_module = getattr(self, f'film_{film_ptr}')
                x = film_module(x, temporal_features)
                film_ptr += 1

        return x

# ======================== Focal Loss ========================

class FocalLoss(nn.Module):
    """Focal loss for handling extreme class imbalance."""
    def __init__(self, alpha: float = 0.25, gamma: float = 2.0):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma
    
    def forward(self, logits: torch.Tensor, targets: torch.Tensor) -> torch.Tensor:
        bce = F.binary_cross_entropy_with_logits(logits, targets, reduction='none')
        pt = torch.exp(-bce)
        focal = self.alpha * (1 - pt) ** self.gamma * bce
        return focal.mean()


# ======================== Data Loading ========================

def get_unique_watersheds(dataset) -> List[str]:
    """Extract unique watershed IDs from dataset."""
    uniques = set()
    scanner = make_scanner(dataset, columns=["Watershed_Loc"])
    for batch in scanner_batches(scanner):
        if batch.num_columns == 0 or batch.num_rows == 0:
            continue
        s = pa.Table.from_batches([batch]).column(0).to_pandas()
        uniques.update(s.dropna().unique().tolist())
    return sorted(uniques)


def get_feature_columns(dataset) -> Tuple[List[str], List[str]]:
    """Get static and temporal feature columns."""
    exclude = {"Target", "Date", "Watershed_Loc", "PixelRow", "PixelCol"}
    
    static_cols = []
    temporal_cols = []
    
    for field in dataset.schema:
        if field.name in exclude:
            continue
        t = field.type
        if pa.types.is_floating(t) or pa.types.is_integer(t):
            if field.name.startswith("subwatershed_") or field.name.startswith("river_basin_"):
                temporal_cols.append(field.name)
            else:
                static_cols.append(field.name)
    
    return static_cols, temporal_cols


def compute_normalization_stats(dataset, feature_cols: List[str], 
                               watersheds: List[str], batch_size: int = 250000):
    """Compute mean and std for normalization."""
    k = len(feature_cols)
    s = np.zeros(k, dtype=np.float64)
    ss = np.zeros(k, dtype=np.float64)
    cnt = np.zeros(k, dtype=np.float64)
    
    filt = ds.field("Watershed_Loc").isin(pa.array(watersheds))
    scanner = make_scanner(dataset, columns=feature_cols, filter=filt, batch_size=batch_size)
    
    for batch in scanner_batches(scanner):
        if batch.num_rows == 0:
            continue
        X = arrow_batch_to_numpy(batch, feature_cols).astype(np.float64, copy=False)
        
        mask = ~np.isnan(X)
        s += np.nansum(X, axis=0)
        ss += np.nansum(X**2, axis=0)
        cnt += mask.sum(axis=0)
    
    means = s / np.maximum(cnt, 1.0)
    var = (ss / np.maximum(cnt, 1.0)) - means**2
    var[var < 0] = 0
    stds = np.sqrt(var)
    stds[stds == 0] = 1.0
    
    return means.astype(np.float32), stds.astype(np.float32)
    
    
class MLPFloodDataset(Dataset):
    """Optimized dataset with watershed-level sampling."""
    def __init__(self, parquet_dataset, watersheds: List[str],
                 static_cols: List[str], temporal_cols: List[str],
                 means: np.ndarray, stds: np.ndarray,
                 config: MLPConfig, mode: str = "train"):
        
        self.static_cols = static_cols
        self.temporal_cols = temporal_cols
        self.all_cols = static_cols + temporal_cols
        self.means = means
        self.stds = stds
        self.config = config
        self.mode = mode
        
        self.static_data = []
        self.temporal_data = []
        self.targets = []
        
        print(f"  Loading {mode} data from {len(watersheds)} watersheds...")
        
        for ws in watersheds:
            filt = ds.field("Watershed_Loc") == ws
            cols = self.all_cols + ["Target"]
            scanner = make_scanner(parquet_dataset, columns=cols, filter=filt, batch_size=500000)
            
            ws_static_list = []
            ws_temporal_list = []
            ws_target_list = []
            
            for batch in scanner_batches(scanner):
                if batch.num_rows == 0:
                    continue
                
                X_static = arrow_batch_to_numpy(batch, static_cols) if static_cols else np.empty((batch.num_rows, 0), dtype=np.float32)
                X_temporal = arrow_batch_to_numpy(batch, temporal_cols) if temporal_cols else np.empty((batch.num_rows, 0), dtype=np.float32)
                y = batch.column("Target").to_numpy(zero_copy_only=False).astype(np.uint8)
                
                ws_static_list.append(X_static)
                ws_temporal_list.append(X_temporal)
                ws_target_list.append(y)
            
            if not ws_static_list:
                continue
                
            ws_static = np.vstack(ws_static_list)
            ws_temporal = np.vstack(ws_temporal_list)
            ws_target = np.hstack(ws_target_list)
            
            # Training: positive-aware sampling (unchanged)
            if mode == "train":
                n_samples = len(ws_target)
                
                pos_idx = np.where(ws_target == 1)[0]
                neg_idx = np.where(ws_target == 0)[0]
                
                if config.max_pixels_per_watershed and n_samples > config.max_pixels_per_watershed:
                    n_pos = len(pos_idx)
                    if n_pos < config.max_pixels_per_watershed:
                        n_neg_keep = config.max_pixels_per_watershed - n_pos
                        if n_neg_keep > 0 and len(neg_idx) > 0:
                            neg_idx_keep = np.random.choice(neg_idx, min(n_neg_keep, len(neg_idx)), replace=False)
                            keep_idx = np.concatenate([pos_idx, neg_idx_keep])
                        else:
                            keep_idx = pos_idx
                    else:
                        keep_idx = np.random.choice(n_samples, config.max_pixels_per_watershed, replace=False)
                    
                    ws_static = ws_static[keep_idx]
                    ws_temporal = ws_temporal[keep_idx]
                    ws_target = ws_target[keep_idx]
                    
                elif config.sample_fraction < 1.0:
                    n_keep = int(n_samples * config.sample_fraction)
                    if n_keep > 0:
                        n_pos = len(pos_idx)
                        if n_pos < n_keep:
                            n_neg_keep = n_keep - n_pos
                            if n_neg_keep > 0 and len(neg_idx) > 0:
                                neg_idx_keep = np.random.choice(neg_idx, min(n_neg_keep, len(neg_idx)), replace=False)
                                keep_idx = np.concatenate([pos_idx, neg_idx_keep])
                            else:
                                keep_idx = pos_idx
                        else:
                            keep_idx = np.random.choice(n_samples, n_keep, replace=False)
                        
                        ws_static = ws_static[keep_idx]
                        ws_temporal = ws_temporal[keep_idx]
                        ws_target = ws_target[keep_idx]
            
            # NEW: Validation - stratified downsampling for speed
            elif mode == "val":
                n_samples = len(ws_target)
                if config.max_val_pixels_per_watershed and n_samples > config.max_val_pixels_per_watershed:
                    # Check if we have both classes (required for stratification)
                    if len(np.unique(ws_target)) > 1:
                        try:
                            from sklearn.model_selection import train_test_split
                            # Stratified sample to preserve exact class balance
                            indices = np.arange(n_samples)
                            _, keep_idx = train_test_split(
                                indices,
                                test_size=config.max_val_pixels_per_watershed,
                                stratify=ws_target,
                                random_state=config.seed
                            )
                            ws_static = ws_static[keep_idx]
                            ws_temporal = ws_temporal[keep_idx]
                            ws_target = ws_target[keep_idx]
                            print(f"    {ws}: stratified val sample {n_samples:,} → {len(keep_idx):,} "
                                  f"(preserved {ws_target.mean():.4f} positive rate)")
                        except ValueError:
                            # Fallback to random if stratification fails (very rare edge case)
                            keep_idx = np.random.choice(n_samples, config.max_val_pixels_per_watershed, replace=False)
                            ws_static = ws_static[keep_idx]
                            ws_temporal = ws_temporal[keep_idx]
                            ws_target = ws_target[keep_idx]
                    else:
                        # Single class - just random sample
                        keep_idx = np.random.choice(n_samples, config.max_val_pixels_per_watershed, replace=False)
                        ws_static = ws_static[keep_idx]
                        ws_temporal = ws_temporal[keep_idx]
                        ws_target = ws_target[keep_idx]
            
            # Test mode: no sampling (full data for accurate metrics)
            # (implicit - no elif needed)
            
            self.static_data.append(ws_static)
            self.temporal_data.append(ws_temporal)
            self.targets.append(ws_target)
                   
            self.static_data.append(ws_static)
            self.temporal_data.append(ws_temporal)
            self.targets.append(ws_target)
        
        if self.static_data:
            self.static_data = np.vstack(self.static_data)
            self.temporal_data = np.vstack(self.temporal_data)
            self.targets = np.hstack(self.targets)
            
            static_idx = len(static_cols)
            self.static_data = np.where(np.isnan(self.static_data), 
                                       means[:static_idx], self.static_data)
            self.static_data = (self.static_data - means[:static_idx]) / stds[:static_idx]
            
            self.temporal_data = np.where(np.isnan(self.temporal_data),
                                         means[static_idx:], self.temporal_data)
            self.temporal_data = (self.temporal_data - means[static_idx:]) / stds[static_idx:]
        else:
            self.static_data = np.empty((0, len(static_cols)), dtype=np.float32)
            self.temporal_data = np.empty((0, len(temporal_cols)), dtype=np.float32)
            self.targets = np.empty((0,), dtype=np.uint8)
        
        print(f"    Loaded {len(self.targets):,} pixels ({self.targets.mean():.4f} positive rate)")
    
    def __len__(self):
        return len(self.targets)
    
    def __getitem__(self, idx):
        return (torch.from_numpy(self.static_data[idx].copy()).float(),
                torch.from_numpy(self.temporal_data[idx].copy()).float(),
                torch.tensor(self.targets[idx], dtype=torch.float32))

def compute_metrics(y_true: np.ndarray, y_proba: np.ndarray, threshold: float) -> Dict[str, float]:
    """Compute comprehensive metrics."""
    y_pred = (y_proba >= threshold).astype(np.uint8)

    prevalence = float(np.mean(y_true))
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    kappa = cohen_kappa_score(y_true, y_pred)
    mcc = matthews_corrcoef(y_true, y_pred) if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1 else 0.0

    try:
        roc_auc = roc_auc_score(y_true, y_proba)
    except Exception:
        roc_auc = np.nan
    try:
        pr_auc = average_precision_score(y_true, y_proba)
    except Exception:
        pr_auc = np.nan

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    tnr = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    return {
        "threshold": float(threshold),
        "prevalence": prevalence,
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "kappa": kappa,
        "mcc": mcc,
        "roc_auc": roc_auc,
        "pr_auc": pr_auc,
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "tpr": tpr,
        "tnr": tnr,
    }


def evaluate_dataset(model: nn.Module,
                     dataset: Dataset,
                     device: str,
                     batch_size: int,
                     threshold: float) -> Dict[str, float]:
    """Run inference and compute metrics."""
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0, pin_memory=(device == "cuda"))
    model.eval()

    all_probs = []
    all_targets = []

    with torch.no_grad():
        for static, temporal, targets in loader:
            static = static.to(device, non_blocking=True)
            temporal = temporal.to(device, non_blocking=True)
            logits = model(static, temporal).squeeze()
            probs = torch.sigmoid(logits).float().cpu().numpy()
            all_probs.append(probs)
            all_targets.append(targets.float().cpu().numpy())

    if not all_probs:
        return compute_metrics(np.array([], dtype=np.uint8), np.array([], dtype=np.float32), threshold)

    y_proba = np.concatenate(all_probs)
    y_true = np.concatenate(all_targets).astype(np.uint8)
    return compute_metrics(y_true, y_proba, threshold)


# ======================== Training ========================

class Trainer:
    """Handles model training."""
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
        
        if config.focal_loss:
            self.criterion = FocalLoss(config.focal_alpha, config.focal_gamma)
        else:
            self.criterion = nn.BCEWithLogitsLoss(
                pos_weight=torch.tensor(config.pos_weight) if config.use_weighted_loss else None
            )
        
        self.scaler = torch.amp.GradScaler('cuda') if config.use_amp and config.device == 'cuda' else None
    
    def train_epoch(self, dataloader):
        self.model.train()
        total_loss = 0
        n_batches = 0
        
        for static, temporal, targets in dataloader:
            static = static.to(self.device, non_blocking=True)
            temporal = temporal.to(self.device, non_blocking=True)
            targets = targets.to(self.device, non_blocking=True)
            
            self.optimizer.zero_grad(set_to_none=True)
            
            if self.config.use_amp and self.scaler:
                with torch.amp.autocast('cuda'):
                    logits = self.model(static, temporal).squeeze()
                    loss = self.criterion(logits, targets)
                
                self.scaler.scale(loss).backward()
                self.scaler.unscale_(self.optimizer)
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
                self.scaler.step(self.optimizer)
                self.scaler.update()
            else:
                logits = self.model(static, temporal).squeeze()
                loss = self.criterion(logits, targets)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
                self.optimizer.step()
            
            total_loss += loss.item()
            n_batches += 1
        
        return total_loss / max(n_batches, 1)
    
    def validate(self, dataloader):
        self.model.eval()
        all_probs = []
        all_targets = []
        total_loss = 0.0
        n_batches = 0

        with torch.no_grad():
            for static, temporal, targets in dataloader:
                static = static.to(self.device, non_blocking=True)
                temporal = temporal.to(self.device, non_blocking=True)
                targets = targets.to(self.device, non_blocking=True)

                if self.scaler is not None:
                    with torch.amp.autocast('cuda'):
                        logits = self.model(static, temporal).squeeze()
                        loss = self.criterion(logits, targets)
                else:
                    logits = self.model(static, temporal).squeeze()
                    loss = self.criterion(logits, targets)

                total_loss += float(loss.item())
                n_batches += 1

                probs = torch.sigmoid(logits).float().cpu().numpy()
                all_probs.append(probs)
                all_targets.append(targets.float().cpu().numpy())

        if n_batches == 0:
            return {"loss": 0.0, "f1": 0.0, "threshold": 0.5}

        all_probs = np.concatenate(all_probs) if len(all_probs) else np.array([])
        all_targets = np.concatenate(all_targets) if len(all_targets) else np.array([])

        best_thresh = 0.5
        best_f1 = 0.0
        if all_targets.size > 0 and np.unique(all_targets).size > 1:
            precisions, recalls, thresholds = precision_recall_curve(all_targets, all_probs)
            if thresholds.size > 0:
                f1 = 2 * precisions[:-1] * recalls[:-1] / np.clip(precisions[:-1] + recalls[:-1], 1e-8, None)
                best_idx = int(np.nanargmax(f1))
                best_f1 = float(f1[best_idx]) if np.isfinite(f1[best_idx]) else 0.0
                best_thresh = float(thresholds[best_idx])

        return {
            "loss": total_loss / max(n_batches, 1),
            "f1": best_f1,
            "threshold": best_thresh
        }


# ======================== NEW: Hyperparameter Tuning with Optuna ========================

def tune_hyperparameters_optuna(dataset, fold, static_cols, temporal_cols, 
                                all_cols, config: MLPConfig) -> Dict[str, Any]:
    """
    Tune hyperparameters using Optuna (n=20 trials).
    Returns best hyperparameters.
    """
    print("  Hyperparameter tuning with Optuna...")
    
    # Compute normalization on train fold
    means, stds = compute_normalization_stats(dataset, all_cols, fold['train'])
    
    # Load train/val for tuning
    train_dataset = MLPFloodDataset(
        dataset, fold['train'], static_cols, temporal_cols,
        means, stds, config, mode='train'
    )
    val_dataset = MLPFloodDataset(
        dataset, fold['val'], static_cols, temporal_cols,
        means, stds, config, mode='val'
    )
    
    if len(train_dataset) == 0 or len(val_dataset) == 0:
        print("    Warning: insufficient data for tuning")
        return {}
    
    def objective(trial):
        # Suggest hyperparameters
        trial_config = MLPConfig(
            hidden_dims=trial.suggest_categorical('hidden_dims', [
                [256, 128, 64],
                [256, 128],
                [128, 64],
                [256, 256, 128, 64]
            ]),
            dropout=trial.suggest_float('dropout', 0.2, 0.5),
            learning_rate=trial.suggest_float('learning_rate', 5e-4, 2e-3, log=True),
            weight_decay=trial.suggest_float('weight_decay', 1e-5, 5e-4, log=True),
            use_interaction_features=trial.suggest_categorical('use_interaction_features', [True, False]),
            n_interaction_features=trial.suggest_int('n_interaction_features', 32, 80),
            # Keep other params from base config
            use_film=config.use_film,
            focal_loss=config.focal_loss,
            focal_alpha=config.focal_alpha,
            focal_gamma=config.focal_gamma,
            batch_size=config.batch_size,
            epochs=10,  # Short training for tuning
            patience=3,
            device=config.device,
            seed=config.seed,
            use_amp=config.use_amp
        )
        
        # Create model and trainer
        model = MLPWithFiLM(len(static_cols), len(temporal_cols), trial_config)
        trainer = Trainer(model, trial_config)
        
        # Create loaders
        train_loader = DataLoader(
            train_dataset, batch_size=trial_config.batch_size,
            shuffle=True, num_workers=0, pin_memory=(config.device == "cuda"),
            drop_last=True
        )
        val_loader = DataLoader(
            val_dataset, batch_size=trial_config.batch_size * 2,
            shuffle=False, num_workers=0, pin_memory=(config.device == "cuda")
        )
        
        # Train briefly
        best_val_f1 = 0.0
        for epoch in range(trial_config.epochs):
            trainer.train_epoch(train_loader)
            val_metrics = trainer.validate(val_loader)
            if val_metrics['f1'] > best_val_f1:
                best_val_f1 = val_metrics['f1']
        
        # Cleanup
        del model, trainer
        gc.collect()
        torch.cuda.empty_cache()
        
        return best_val_f1
    
    # Run Optuna optimization
    study = optuna.create_study(
        direction='maximize',
        sampler=TPESampler(seed=config.seed)
    )
    study.optimize(objective, n_trials=config.n_tuning_iterations, show_progress_bar=False)
    
    print(f"    Best F1: {study.best_value:.4f}")
    print(f"    Best params: {study.best_params}")
    
    return study.best_params


# ======================== NEW: Diagnostic Plots ========================

def create_diagnostic_plots(model_results: Dict, output_folder: str):
    """Create diagnostic plots matching glmFloods.py style."""
    print("\nCreating diagnostic plots...")
    
    cv_results = model_results.get('cv_results')
    
    # Helper to save and upload
    def save_and_upload(fig, s3_path):
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            temp_path = tmp.name
        fig.savefig(temp_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        if s3_path.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_path, bucket, key)
        else:
            import shutil
            os.makedirs(os.path.dirname(s3_path), exist_ok=True)
            shutil.move(temp_path, s3_path)
        os.remove(temp_path) if os.path.exists(temp_path) else None
        print(f"  Created: {s3_path}")
    
    # 1. CV Performance
    if cv_results is not None and not cv_results.empty:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        metrics_to_plot = ['mcc', 'f1', 'precision', 'recall']
        titles = ['MCC', 'F1 Score', 'Precision', 'Recall']
        
        for idx, (metric, title) in enumerate(zip(metrics_to_plot, titles)):
            ax = axes[idx // 2, idx % 2]
            train_vals = cv_results[cv_results['set'] == 'train'][metric].values if metric in cv_results.columns else []
            val_vals = cv_results[cv_results['set'] == 'val'][metric].values if metric in cv_results.columns else []
            
            x = np.arange(max(len(train_vals), len(val_vals)))
            width = 0.35
            
            tv = train_vals if len(train_vals) == len(x) else np.zeros_like(x, dtype=float)
            vv = val_vals if len(val_vals) == len(x) else np.zeros_like(x, dtype=float)
            
            ax.bar(x - width/2, tv, width, label='Train', color='steelblue', alpha=0.8)
            ax.bar(x + width/2, vv, width, label='Val', color='coral', alpha=0.8)
            ax.set_xlabel('Fold')
            ax.set_ylabel(metric.upper())
            ax.set_title(title)
            ax.set_xticks(x)
            ax.set_xticklabels([f'{i+1}' for i in range(len(x))])
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        save_and_upload(fig, f"{output_folder}/mlp_cv_performance.png")
    
    # 2. Threshold Distribution
    if cv_results is not None and 'threshold' in cv_results.columns:
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(cv_results['threshold'], bins=30, color='steelblue', alpha=0.8)
        ax.set_title('Distribution of CV Thresholds (F1-optimal)')
        ax.set_xlabel('Threshold')
        ax.set_ylabel('Count')
        ax.grid(alpha=0.3)
        plt.tight_layout()
        save_and_upload(fig, f"{output_folder}/mlp_threshold_distribution.png")
    
    print("✅ Diagnostic plots complete!")


# ======================== Main Training Pipeline ========================

def spatial_kfold_split(watersheds: List[str], config: MLPConfig):
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
        
        n_train = int(config.train_frac * len(shuffled))
        n_val = int(config.val_frac * len(shuffled))
        
        train_ws = shuffled[:n_train]
        val_ws = shuffled[n_train:n_train + n_val]
        
        folds.append({
            'train': train_ws,
            'val': val_ws,
            'test': test_ws
        })
    
    return folds


def train_mlp_flood_model(parquet_uri: str, config: Optional[MLPConfig] = None,
                         save_path: Optional[str] = None) -> Dict:
    """Train MLP+FiLM flood model with k-fold CV."""
    
    config = config or MLPConfig()
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)
    
    print("=" * 70)
    print("MLP + FiLM FLOOD MODEL - K-FOLD SPATIAL CV")
    print("=" * 70)
    print(f"Hyperparameter tuning: {'ENABLED' if config.tune_hyperparameters else 'DISABLED'}")
    
    dataset = open_parquet_dataset(parquet_uri)
    watersheds = get_unique_watersheds(dataset)
    static_cols, temporal_cols = get_feature_columns(dataset)
    all_cols = static_cols + temporal_cols
    
    print(f"Watersheds: {len(watersheds)}")
    print(f"Static features: {len(static_cols)}")
    print(f"Temporal features: {len(temporal_cols)}")
    
    folds = spatial_kfold_split(watersheds, config)
    print(f"Test watersheds: {folds[0]['test']}")
    
    # NEW: Optional hyperparameter tuning
    best_params = {}
    if config.tune_hyperparameters:
        best_params = tune_hyperparameters_optuna(
            dataset, folds[0], static_cols, temporal_cols, all_cols, config
        )
        # Update config with best params
        if best_params:
            for key, value in best_params.items():
                if hasattr(config, key):
                    setattr(config, key, value)
    
    cv_thresholds = []
    fold_rows = [] 
    cv_metrics = []
    
    for fold_idx, fold in enumerate(folds):
        print(f"\nFold {fold_idx + 1}/{len(folds)}")
        print(f"  Train: {len(fold['train'])} watersheds")
        print(f"  Val: {len(fold['val'])} watersheds")
        
        means, stds = compute_normalization_stats(dataset, all_cols, fold['train'])
        
        train_dataset = MLPFloodDataset(
            dataset, fold['train'], static_cols, temporal_cols,
            means, stds, config, mode='train'
        )
        
        val_dataset = MLPFloodDataset(
            dataset, fold['val'], static_cols, temporal_cols,
            means, stds, config, mode='val'
        )
        
        if len(train_dataset) == 0 or len(val_dataset) == 0:
            print("  Warning: empty dataset")
            cv_thresholds.append(0.5)
            continue
        
        train_loader = DataLoader(
            train_dataset,
            batch_size=config.batch_size,
            shuffle=True,
            num_workers=2 if len(train_dataset) > 10000 else 0,
            pin_memory=(config.device == "cuda"),
            persistent_workers=True if len(train_dataset) > 10000 else False,
            prefetch_factor=2 if len(train_dataset) > 10000 else None,
            drop_last=True
        )

        val_loader = DataLoader(
            val_dataset,
            batch_size=config.batch_size * 2,
            shuffle=False,
            num_workers=2 if len(val_dataset) > 10000 else 0,
            pin_memory=(config.device == "cuda"),
            persistent_workers=True if len(val_dataset) > 10000 else False,
            prefetch_factor=2 if len(val_dataset) > 10000 else None,
            drop_last=False
        )
        
        model = MLPWithFiLM(len(static_cols), len(temporal_cols), config)
        trainer = Trainer(model, config)
                
        best_val_loss = float('inf')
        patience_counter = 0
        best_threshold = 0.5
        best_model_state = None  # NEW: save best weights
        
        for epoch in range(config.epochs):
            train_loss = trainer.train_epoch(train_loader)
            val_metrics = trainer.validate(val_loader)
            
            trainer.scheduler.step(val_metrics['loss'])
            
            print(f"  Epoch {epoch+1}: Train loss={train_loss:.4f}, "
                  f"Val loss={val_metrics['loss']:.4f}, Val F1={val_metrics['f1']:.4f}")
            
            # NEW: Use RELATIVE improvement threshold
            improvement_threshold = best_val_loss * (1 - config.min_delta)
            
            if val_metrics['loss'] < improvement_threshold:
                best_val_loss = val_metrics['loss']
                best_threshold = val_metrics['threshold']
                best_model_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}  # NEW
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= config.patience:
                    print(f"  Early stopping at epoch {epoch+1}")
                    break
        
        # NEW: Restore best weights
        if best_model_state is not None:
            model.load_state_dict(best_model_state)
        
        cv_metrics.append(val_metrics)
        train_metrics_row = evaluate_dataset(trainer.model, train_dataset, config.device, config.batch_size, best_threshold)
        val_metrics_row = evaluate_dataset(trainer.model, val_dataset, config.device, config.batch_size, best_threshold)
        
        train_metrics_row = {"fold": fold_idx + 1, "set": "train", **train_metrics_row}
        val_metrics_row = {"fold": fold_idx + 1, "set": "val", **val_metrics_row}
        fold_rows.extend([train_metrics_row, val_metrics_row])
        
        cv_thresholds.append(best_threshold)
        
        del model, trainer
        gc.collect()
        torch.cuda.empty_cache()

    folds_df = pd.DataFrame(fold_rows) if fold_rows else pd.DataFrame()
    summary_metrics = {}
    
    if not folds_df.empty:
        for s in ["train", "val"]:
            sub = folds_df[folds_df["set"] == s]
            if not sub.empty:
                summary_metrics[s] = sub.drop(columns=["fold", "set"]).mean(numeric_only=True).to_dict()
    
    print("\n" + "=" * 70)
    print("FINAL MODEL TRAINING")
    print("=" * 70)
    
    final_threshold = np.mean(cv_thresholds)
    print(f"Average threshold from CV: {final_threshold:.3f}")
    
    test_ws = folds[0]['test']
    train_ws = [w for w in watersheds if w not in test_ws]
    
    means, stds = compute_normalization_stats(dataset, all_cols, train_ws)
    
    final_train_dataset = MLPFloodDataset(
        dataset, train_ws, static_cols, temporal_cols,
        means, stds, config, mode='train'
    )
    
    final_loader = DataLoader(
        final_train_dataset,
        batch_size=config.batch_size,
        shuffle=True,
        num_workers=config.num_workers,
        pin_memory=(config.device == "cuda"),
        drop_last=True
    )
    
    final_model = MLPWithFiLM(len(static_cols), len(temporal_cols), config)
    final_trainer = Trainer(final_model, config)
    
    for epoch in range(config.epochs):
        train_loss = final_trainer.train_epoch(final_loader)
        print(f"Epoch {epoch+1}: Loss={train_loss:.4f}")
    
    test_dataset = MLPFloodDataset(
        dataset, test_ws, static_cols, temporal_cols,
        means, stds, config, mode='val'
    )
    test_metrics_row = evaluate_dataset(final_model, test_dataset, config.device, config.batch_size, final_threshold)
    summary_metrics["test"] = test_metrics_row
    
    results = {
        'model': final_model,
        'threshold': final_threshold,
        'static_cols': static_cols,
        'temporal_cols': temporal_cols,
        'means': means,
        'stds': stds,
        'config': config,
        'test_watersheds': test_ws,
        'cv_metrics': cv_metrics,
        'metrics_summary': summary_metrics,
        'cv_results': folds_df,  # NEW: for diagnostic plots
        'best_hyperparameters': best_params  # NEW
    }
    
    if save_path:
        torch.save(results, save_path)
        print(f"Model saved to {save_path}")
    
    return results


# ======================== Visualization ========================

class MLPFloodVisualizer:
    """Generate visualizations matching XGBoost format."""
    
    def __init__(self, model_results: Dict, parquet_uri: str, training_folder: str):
        self.model_results = model_results
        self.model = model_results['model']
        self.threshold = model_results['threshold']
        self.static_cols = model_results['static_cols']
        self.temporal_cols = model_results['temporal_cols']
        self.means = model_results['means']
        self.stds = model_results['stds']
        self.config = model_results['config']
        
        self.parquet_uri = parquet_uri
        self.training_folder = training_folder
        
        m = re.search(r"Lat_([-\d\.]+)_Lon_([-\d\.]+)", training_folder)
        if m:
            self.test_lat = m.group(1)
            self.test_lon = m.group(2)
            self.test_watershed = f"Lat_{self.test_lat}_Lon_{self.test_lon}"
        
        self.s3_client = boto3.client('s3')
        configure_gdal()
    
    def generate_all_outputs(self, output_folder: Optional[str] = None):
        """Generate all visualization outputs."""
        
        if output_folder is None:
            output_folder = f"{self.training_folder}/assessMLPModel"
        
        print("\n" + "="*60)
        print("GENERATING MLP MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        print(f"Test watershed: {self.test_watershed}")
        print(f"Output folder: {output_folder}")
        
        # Create diagnostic plots
        create_diagnostic_plots(self.model_results, output_folder)
        
        predictions, confusions, per_date_metrics = self._generate_predictions()
        
        # Load reference TIFF for georeferencing
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
        
        dates = sorted(predictions.keys())
        
        # Create output paths
        pred_tiff = f"{output_folder}/mlp_predictions_{self.test_watershed}.tif"
        cm_tiff = f"{output_folder}/mlp_confusion_matrix_{self.test_watershed}.tif"  # ADD THIS LINE
        
        self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "MLP_Prediction")
        self._write_multiband_tiff(ref_ds, confusions, dates, cm_tiff, "MLP_CM")
        
        # Create PNGs
        png_paths = self._create_pngs(confusions, dates, per_date_metrics, output_folder)
        
        # Create CSVs
        csv_paths = self._write_csvs(output_folder, per_date_metrics)
        
        print("\n✅ MLP visualization complete!")
        
        return {
            'prediction_tiff': pred_tiff,
            'confusion_matrix_tiff': cm_tiff,
            'png_visualizations': png_paths,
            'csv_outputs': csv_paths
        }


    def _generate_predictions(self):
        """Generate predictions for test watershed."""
        
        dataset = open_parquet_dataset(self.parquet_uri)
        
        filt = ds.field("Watershed_Loc") == self.test_watershed
        scanner = make_scanner(dataset, filter=filt)
        
        pixel_data = []
        for batch in scanner_batches(scanner):
            if batch.num_rows > 0:
                df = pa.Table.from_batches([batch]).to_pandas()
                pixel_data.append(df)
        
        if not pixel_data:
            return {}, {}, {}
        
        df_all = pd.concat(pixel_data, ignore_index=True)
        
        dates = sorted(df_all['Date'].unique())
        
        H = df_all['PixelRow'].max() + 1
        W = df_all['PixelCol'].max() + 1
        
        predictions = {}
        confusions = {}
        per_date_metrics = {}
        
        self.model.eval()
        
        print(f"  Processing {len(dates)} dates...")
        for date_idx, date_str in enumerate(dates):
            if (date_idx + 1) % 10 == 0:
                print(f"    Date {date_idx + 1}/{len(dates)}")
            
            date_data = df_all[df_all['Date'] == date_str]
            
            X_static = date_data[self.static_cols].to_numpy(dtype=np.float32)
            X_temporal = date_data[self.temporal_cols].to_numpy(dtype=np.float32)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8)
            rows = date_data['PixelRow'].to_numpy()
            cols = date_data['PixelCol'].to_numpy()
            
            static_idx = len(self.static_cols)
            X_static = np.where(np.isnan(X_static), self.means[:static_idx], X_static)
            X_static = (X_static - self.means[:static_idx]) / self.stds[:static_idx]
            
            X_temporal = np.where(np.isnan(X_temporal), self.means[static_idx:], X_temporal)
            X_temporal = (X_temporal - self.means[static_idx:]) / self.stds[static_idx:]
            
            all_probs = []
            
            with torch.no_grad():
                for i in range(0, len(X_static), self.config.batch_size_inference):
                    batch_static = torch.from_numpy(X_static[i:i+self.config.batch_size_inference])
                    batch_temporal = torch.from_numpy(X_temporal[i:i+self.config.batch_size_inference])
                    
                    batch_static = batch_static.to(self.config.device)
                    batch_temporal = batch_temporal.to(self.config.device)
                    
                    logits = self.model(batch_static, batch_temporal).squeeze()
                    probs = torch.sigmoid(logits).cpu().numpy()
                    all_probs.append(probs)
            
            all_probs = np.hstack(all_probs)
            y_pred = (all_probs >= self.threshold).astype(np.uint8)
            
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            pred_img[rows, cols] = y_pred
            predictions[date_str] = pred_img
            
            cm_img = np.full((H, W), 255, dtype=np.uint8)
            cm_vals = np.where(
                y_true == 0,
                np.where(y_pred == 0, 0, 3),
                np.where(y_pred == 1, 1, 2)
            )
            cm_img[rows, cols] = cm_vals
            confusions[date_str] = cm_img
            
            try:
                roc_auc = roc_auc_score(y_true, all_probs)
            except:
                roc_auc = np.nan
            
            try:
                pr_auc = average_precision_score(y_true, all_probs)
            except:
                pr_auc = np.nan
            
            if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1:
                mcc = matthews_corrcoef(y_true, y_pred)
            else:
                mcc = 0.0
            
            per_date_metrics[date_str] = {
                'date': date_str,
                'F1': f1_score(y_true, y_pred, zero_division=0),
                'Accuracy': accuracy_score(y_true, y_pred),
                'Precision': precision_score(y_true, y_pred, zero_division=0),
                'Recall': recall_score(y_true, y_pred, zero_division=0),
                'ROC_AUC': float(roc_auc) if np.isfinite(roc_auc) else np.nan,
                'PR_AUC': float(pr_auc) if np.isfinite(pr_auc) else np.nan,
                'Kappa': cohen_kappa_score(y_true, y_pred),
                'MCC': float(mcc),
                'actual_pos_frac': float(y_true.mean()),
                'pred_pos_frac': float(y_pred.mean())
            }
        
        return predictions, confusions, per_date_metrics

    def _write_multiband_tiff(self, ref_ds, arrays, dates, output_path, prefix):
        """Write multiband TIFF."""
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
            temp_path = tmp.name
        
        driver = gdal.GetDriverByName('GTiff')
        out = driver.Create(
            temp_path, ref_ds.RasterXSize, ref_ds.RasterYSize, len(dates),
            gdal.GDT_Byte, options=["COMPRESS=DEFLATE", "TILED=YES"]
        )
        
        out.SetGeoTransform(ref_ds.GetGeoTransform())
        out.SetProjection(ref_ds.GetProjection())
        
        for i, date_str in enumerate(dates, 1):
            band = out.GetRasterBand(i)
            band.SetNoDataValue(255)
            band.SetDescription(f"{prefix}_{date_str}")
            band.WriteArray(arrays.get(date_str, np.full((ref_ds.RasterYSize, ref_ds.RasterXSize), 255)))
        
        out = None
        
        # Upload to S3
        bucket, key = output_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(temp_path, bucket, key)
        os.remove(temp_path)
        print(f"  Created: {output_path}")
    
    def _create_pngs(self, confusions, dates, metrics, output_folder):
        """Create PNG visualizations (one per date), labeled with metrics."""
        png_paths = []
    
        for date_str in sorted(dates):
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
            ax.set_title(f'MLP+FiLM Confusion Matrix - {date_str}')
    
            text = (
                f"Date: {date_str}\n"
                f"Model: MLP+FiLM\n"
                f"F1: {m.get('F1', 0):.3f}\n"
                f"Precision: {m.get('Precision', 0):.3f}\n"
                f"Recall: {m.get('Recall', 0):.3f}\n"
                f"ROC AUC: {m.get('ROC_AUC', 0):.3f}"
            )
            ax.text(0.02, 0.98, text, transform=ax.transAxes, va='top',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
            ax.axis('off')
            plt.tight_layout()
    
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()
    
            s3_path = f"{output_folder}/mlp_confusion_matrix_{date_str}.png"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(temp_png, bucket, key)
            os.remove(temp_png)
            png_paths.append(s3_path)
    
        print(f"  Created {len(png_paths)} PNG visualizations")
        return png_paths

    def _write_csvs(self, output_folder, per_date_metrics):
        """Write CSV outputs matching GLM format."""
        paths = {}
        
        # Per-date metrics - matching GLM column order exactly
        if per_date_metrics:
            df = pd.DataFrame.from_dict(per_date_metrics, orient='index')
            # Ensure consistent column order matching GLM
            column_order = ['date', 'F1', 'Accuracy', 'Precision', 'Recall', 'ROC_AUC', 
                           'PR_AUC', 'Kappa', 'MCC', 'actual_pos_frac', 'pred_pos_frac']
            df = df[column_order]
            
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
            df.to_csv(tmp_csv, index=False)
            s3_path = f"{output_folder}/mlp_per_date_metrics_{self.test_watershed}.csv"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(tmp_csv, bucket, key)
            os.remove(tmp_csv)
            paths['per_date_metrics'] = s3_path
            print(f"  Created: {s3_path}")
        
        # Summary metrics (train/val/test), aligned to XGBoost-like schema
        metrics_summary = self.model_results.get("metrics_summary", {})
        if metrics_summary:
            rows = []
            for split in ["train", "val", "test"]:
                m = metrics_summary.get(split)
                if not m:
                    continue
                rows.append({
                    "set": split,
                    "threshold": m.get("threshold", self.threshold),
                    "prevalence": m.get("prevalence", np.nan),
                    "accuracy": m.get("accuracy", np.nan),
                    "precision": m.get("precision", np.nan),
                    "recall": m.get("recall", np.nan),
                    "f1": m.get("f1", np.nan),
                    "kappa": m.get("kappa", np.nan),
                    "mcc": m.get("mcc", np.nan),
                    "roc_auc": m.get("roc_auc", np.nan),
                    "pr_auc": m.get("pr_auc", np.nan),
                    "tp": m.get("tp", np.nan),
                    "fp": m.get("fp", np.nan),
                    "tn": m.get("tn", np.nan),
                    "fn": m.get("fn", np.nan),
                    "tpr": m.get("tpr", np.nan),
                    "tnr": m.get("tnr", np.nan),
                })
            if rows:
                df_sum = pd.DataFrame(rows)
                with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                    tmp_csv = tmp.name
                df_sum.to_csv(tmp_csv, index=False)
                s3_path = f"{output_folder}/mlp_metrics_summary.csv"
                bucket, key = s3_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(tmp_csv, bucket, key)
                os.remove(tmp_csv)
                paths['metrics_summary'] = s3_path
        
        # Feature importance placeholder (unchanged)
        importance_df = pd.DataFrame({
            'feature': ['static_features', 'temporal_features', 'interaction_features'],
            'importance': [0.5, 0.3, 0.2]
        })
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
            tmp_csv = tmp.name
        importance_df.to_csv(tmp_csv, index=False)
        s3_path = f"{output_folder}/mlp_feature_importances.csv"
        bucket, key = s3_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(tmp_csv, bucket, key)
        os.remove(tmp_csv)
        paths['feature_importances'] = s3_path
        
        return paths

"""
mlpFilmFloods.py

Pixel-wise MLP with FiLM conditioning for flood prediction.
Combines XGBoost's data loading from parquet with FiLM temporal modulation.
Treats each pixel independently without spatial convolutions.
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
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler


# Parquet I/O
import pyarrow as pa
import pyarrow.dataset as ds

# Visualization
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
    n_interaction_features: int = 50  # Cap on interaction features

    # Normalization for hidden layers: "layer" (robust) or "batch"
    norm_type: str = "layer"

    # Training
    epochs: int = 30
    batch_size: int = 512
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4

    # Early stopping
    patience: int = 5
    min_delta: float = 0.001

    # Class balancing
    use_weighted_loss: bool = True
    pos_weight: float = 10.0
    focal_loss: bool = True
    focal_alpha: float = 0.25
    focal_gamma: float = 2.0

    # Data
    sample_fraction: float = 0.1  # Sample fraction per watershed
    max_pixels_per_watershed: int = 100000

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
        hidden_film = max(64, n_hidden)  # a bit more capacity than n_hidden//2
        self.mlp = nn.Sequential(
            nn.Linear(n_temporal, hidden_film),
            nn.LayerNorm(hidden_film),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_film, n_hidden * 2)  # gamma and beta
        )
    
    def forward(self, features: torch.Tensor, temporal: torch.Tensor) -> torch.Tensor:
        """
        Args:
            features: [B, H] hidden features
            temporal: [B, T] temporal features
        Returns:
            Modulated features [B, H]
        """
        params = self.mlp(temporal)
        gamma, beta = params.chunk(2, dim=1)
        # gentle gating to avoid overwhelming the static embedding
        return features * (1 + torch.tanh(gamma)) + beta

# ======================== MLP + FiLM Model ========================

class MLPWithFiLM(nn.Module):
    """
    Pixel-wise MLP with FiLM conditioning.
    Treats each pixel independently with explicit feature interactions.
    """
    def __init__(self, n_static: int, n_temporal: int, config: MLPConfig):
        super().__init__()
        self.n_static = n_static
        self.n_temporal = n_temporal
        self.config = config
        
        # Interaction feature generator
        if config.use_interaction_features:
            self.n_interactions = min(config.n_interaction_features, n_static * 5)
            self.interaction_weight = nn.Parameter(
                torch.randn(n_static, self.n_interactions) * 0.01
            )
            # NEW: temporal MLP producing per-interaction gates (0..1)
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
        
        # Build MLP hidden layers
        layers = []
        prev_dim = input_dim
        self.hidden_dims = config.hidden_dims[:]
        
        # Choose strategic FiLM layers: first and penultimate (if present and distinct)
        film_layer_indices = []
        if config.use_film and len(self.hidden_dims) >= 1:
            first_idx = 0
            penultimate_idx = max(0, len(self.hidden_dims) - 2)
            film_layer_indices = [first_idx]
            if penultimate_idx != first_idx:
                film_layer_indices.append(penultimate_idx)
        self._film_layer_indices = film_layer_indices  # e.g., [0] or [0, penultimate]
        
        # We will add FiLM modules under names film_0, film_1 in this order
        film_added = 0
        
        for i, hidden_dim in enumerate(self.hidden_dims):
            layers.append(nn.Linear(prev_dim, hidden_dim))
            if self.config.norm_type.lower() == "batch":
                layers.append(nn.BatchNorm1d(hidden_dim))
            else:
                # Default to LayerNorm to avoid batch-size=1 failures in training
                layers.append(nn.LayerNorm(hidden_dim))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.Dropout(config.dropout))
            prev_dim = hidden_dim
        
            # Add FiLM only on strategic layers (and not after the last hidden)
            if config.use_film and (i in film_layer_indices) and (i < len(self.hidden_dims) - 1):
                self.add_module(f'film_{film_added}', FiLMBlock(n_temporal, hidden_dim))
                film_added += 1
        
        # Output layer
        layers.append(nn.Linear(prev_dim, 1))
        self.mlp = nn.ModuleList(layers)
        
        # Compute positions in self.mlp where FiLM should be applied
        # We apply FiLM right after finishing a hidden block (after Dropout), before the next Linear.
        # Each hidden block contributes 4 modules: Linear, Norm, ReLU, Dropout
        self.film_indices = []
        layer_count = 0
        film_counter = 0
        for i in range(len(self.hidden_dims)):
            layer_count += 4  # consumed modules for this hidden block
            if config.use_film and (i in film_layer_indices) and (i < len(self.hidden_dims) - 1):
                # apply FiLM after the Dropout (i.e., before the next Linear)
                self.film_indices.append(layer_count - 1)
                film_counter += 1
                
    def forward(self, static_features: torch.Tensor, temporal_features: torch.Tensor) -> torch.Tensor:
        """
        Args:
            static_features: [B, n_static] spatially varying features
            temporal_features: [B, n_temporal] temporally varying features
        Returns:
            logits: [B, 1]
        """
        # Generate interaction features with temporal gates (if enabled)
        if self.config.use_interaction_features:
            static_proj = torch.matmul(static_features, self.interaction_weight)  # [B, n_interactions]
            gates = self.temporal_gate(temporal_features)                         # [B, n_interactions] in [0,1]
            interactions = static_proj * gates
            x = torch.cat([static_features, interactions], dim=1)
        else:
            x = static_features

        film_ptr = 0
        # Iterate through the MLP layers, injecting FiLM after the chosen hidden blocks.
        # Each hidden block in self.mlp is [Linear, BatchNorm1d, ReLU, Dropout].
        # self.film_indices stores the index of Dropout for blocks where FiLM should be applied.
        for idx, module in enumerate(self.mlp):
            x = module(x)
            if (self.config.use_film
                and film_ptr < len(self.film_indices)
                and idx == self.film_indices[film_ptr]):
                film_module = getattr(self, f'film_{film_ptr}')
                x = film_module(x, temporal_features)
                film_ptr += 1

        return x  # final Linear in self.mlp returns logits [B, 1]

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
    """Compute mean and std for normalization - using helper."""
    k = len(feature_cols)
    s = np.zeros(k, dtype=np.float64)
    ss = np.zeros(k, dtype=np.float64)
    cnt = np.zeros(k, dtype=np.float64)
    
    filt = ds.field("Watershed_Loc").isin(pa.array(watersheds))
    scanner = make_scanner(dataset, columns=feature_cols, filter=filt, batch_size=batch_size)
    
    for batch in scanner_batches(scanner):
        if batch.num_rows == 0:
            continue
        # USE THE HELPER FUNCTION
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
    """Optimized dataset with watershed-level sampling and direct Arrow conversion."""
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
        
        # Load and sample data
        self.static_data = []
        self.temporal_data = []
        self.targets = []
        
        print(f"  Loading {mode} data from {len(watersheds)} watersheds...")
        
        for ws in watersheds:
            filt = ds.field("Watershed_Loc") == ws
            cols = self.all_cols + ["Target"]
            scanner = make_scanner(parquet_dataset, columns=cols, filter=filt, batch_size=500000)
            
            # Collect ALL data for this watershed first
            ws_static_list = []
            ws_temporal_list = []
            ws_target_list = []
            
            for batch in scanner_batches(scanner):
                if batch.num_rows == 0:
                    continue
                
                # USE THE HELPER
                X_static = arrow_batch_to_numpy(batch, static_cols) if static_cols else np.empty((batch.num_rows, 0), dtype=np.float32)
                X_temporal = arrow_batch_to_numpy(batch, temporal_cols) if temporal_cols else np.empty((batch.num_rows, 0), dtype=np.float32)
                y = batch.column("Target").to_numpy(zero_copy_only=False).astype(np.uint8)
                
                ws_static_list.append(X_static)
                ws_temporal_list.append(X_temporal)
                ws_target_list.append(y)
            
            if not ws_static_list:
                continue
                
            # Concatenate all batches for this watershed
            ws_static = np.vstack(ws_static_list)
            ws_temporal = np.vstack(ws_temporal_list)
            ws_target = np.hstack(ws_target_list)
            
            # NOW apply sampling at WATERSHED LEVEL (not per batch!)
            if mode == "train":
                n_samples = len(ws_target)
                
                # Positive-aware sampling like TabNet
                pos_idx = np.where(ws_target == 1)[0]
                neg_idx = np.where(ws_target == 0)[0]
                
                if config.max_pixels_per_watershed and n_samples > config.max_pixels_per_watershed:
                    # Keep all positives if possible, sample negatives
                    n_pos = len(pos_idx)
                    if n_pos < config.max_pixels_per_watershed:
                        # Keep all positives, sample negatives
                        n_neg_keep = config.max_pixels_per_watershed - n_pos
                        if n_neg_keep > 0 and len(neg_idx) > 0:
                            neg_idx_keep = np.random.choice(neg_idx, min(n_neg_keep, len(neg_idx)), replace=False)
                            keep_idx = np.concatenate([pos_idx, neg_idx_keep])
                        else:
                            keep_idx = pos_idx
                    else:
                        # Too many positives, random sample
                        keep_idx = np.random.choice(n_samples, config.max_pixels_per_watershed, replace=False)
                    
                    ws_static = ws_static[keep_idx]
                    ws_temporal = ws_temporal[keep_idx]
                    ws_target = ws_target[keep_idx]
                    
                elif config.sample_fraction < 1.0:
                    n_keep = int(n_samples * config.sample_fraction)
                    if n_keep > 0:
                        # Positive-aware sampling
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
            
            self.static_data.append(ws_static)
            self.temporal_data.append(ws_temporal)
            self.targets.append(ws_target)
        
        # Concatenate all watersheds
        if self.static_data:
            self.static_data = np.vstack(self.static_data)
            self.temporal_data = np.vstack(self.temporal_data)
            self.targets = np.hstack(self.targets)
            
            # Normalize
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
    """Compute a rich set of metrics at a given threshold, similar to the XGBoost summary."""
    y_pred = (y_proba >= threshold).astype(np.uint8)

    # Basic metrics
    prevalence = float(np.mean(y_true))
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    kappa = cohen_kappa_score(y_true, y_pred)
    mcc = matthews_corrcoef(y_true, y_pred) if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1 else 0.0

    # AUCs
    try:
        roc_auc = roc_auc_score(y_true, y_proba)
    except Exception:
        roc_auc = np.nan
    try:
        pr_auc = average_precision_score(y_true, y_proba)
    except Exception:
        pr_auc = np.nan

    # Confusion components
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
    """Run inference over a dataset and compute metrics at the given threshold."""
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


        # Full PR curve thresholding for F1-optimal selection
        best_thresh = 0.5
        best_f1 = 0.0
        if all_targets.size > 0 and np.unique(all_targets).size > 1:
            precisions, recalls, thresholds = precision_recall_curve(all_targets, all_probs)
            # thresholds has length N-1 compared to precisions/recalls length N
            if thresholds.size > 0:
                f1 = 2 * precisions[:-1] * recalls[:-1] / np.clip(precisions[:-1] + recalls[:-1], 1e-8, None)
                best_idx = int(np.nanargmax(f1))
                best_f1 = float(f1[best_idx]) if np.isfinite(f1[best_idx]) else 0.0
                best_thresh = float(thresholds[best_idx])
            else:
                # degenerate case; fall back to 0.5
                best_f1 = 0.0
                best_thresh = 0.5
        else:
            # single-class validation set; keep default threshold
            best_f1 = 0.0
            best_thresh = 0.5

        return {
            "loss": total_loss / max(n_batches, 1),
            "f1": best_f1,
            "threshold": best_thresh
        }

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
    
    # Load dataset
    dataset = open_parquet_dataset(parquet_uri)
    watersheds = get_unique_watersheds(dataset)
    static_cols, temporal_cols = get_feature_columns(dataset)
    all_cols = static_cols + temporal_cols
    
    print(f"Watersheds: {len(watersheds)}")
    print(f"Static features: {len(static_cols)}")
    print(f"Temporal features: {len(temporal_cols)}")
    
    # Create folds
    folds = spatial_kfold_split(watersheds, config)
    print(f"Test watersheds: {folds[0]['test']}")
    
    # Run CV
    cv_thresholds = []
    fold_rows = [] 
    cv_metrics = []
    
    for fold_idx, fold in enumerate(folds):
        print(f"\nFold {fold_idx + 1}/{len(folds)}")
        print(f"  Train: {len(fold['train'])} watersheds")
        print(f"  Val: {len(fold['val'])} watersheds")
        
        # Compute normalization on train only
        means, stds = compute_normalization_stats(dataset, all_cols, fold['train'])
        
        # Create datasets
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
        
        # CV train loader
        train_loader = DataLoader(
            train_dataset,
            batch_size=config.batch_size,
            shuffle=True,
            num_workers=2 if len(train_dataset) > 10000 else 0,  # Use workers for large datasets
            pin_memory=(config.device == "cuda"),
            persistent_workers=True if len(train_dataset) > 10000 else False,
            prefetch_factor=2 if len(train_dataset) > 10000 else None,
            drop_last=True
        )

        # For validation loader  
        val_loader = DataLoader(
            val_dataset,
            batch_size=config.batch_size * 2,  # Larger batch for validation (no gradients)
            shuffle=False,
            num_workers=2 if len(val_dataset) > 10000 else 0,
            pin_memory=(config.device == "cuda"),
            persistent_workers=True if len(val_dataset) > 10000 else False,
            prefetch_factor=2 if len(val_dataset) > 10000 else None,
            drop_last=False
        )
        
        # Create model
        model = MLPWithFiLM(len(static_cols), len(temporal_cols), config)
        trainer = Trainer(model, config)
        
        # Train
        best_val_loss = float('inf')
        patience_counter = 0
        best_threshold = 0.5
        
        for epoch in range(config.epochs):
            train_loss = trainer.train_epoch(train_loader)
            val_metrics = trainer.validate(val_loader)
            
            trainer.scheduler.step(val_metrics['loss'])
            
            print(f"  Epoch {epoch+1}: Train loss={train_loss:.4f}, "
                  f"Val loss={val_metrics['loss']:.4f}, Val F1={val_metrics['f1']:.4f}")
            
            if val_metrics['loss'] < best_val_loss - config.min_delta:
                best_val_loss = val_metrics['loss']
                best_threshold = val_metrics['threshold']
                patience_counter = 0
            else:
                patience_counter += 1
                if patience_counter >= config.patience:
                    print(f"  Early stopping at epoch {epoch+1}")
                    break
        
        cv_metrics.append(val_metrics)
        train_metrics_row = evaluate_dataset(trainer.model, train_dataset, config.device, config.batch_size, best_threshold)
        val_metrics_row = evaluate_dataset(trainer.model, val_dataset, config.device, config.batch_size, best_threshold)
        
        # Tag and store fold metrics
        train_metrics_row = {"fold": fold_idx + 1, "set": "train", **train_metrics_row}
        val_metrics_row = {"fold": fold_idx + 1, "set": "val", **val_metrics_row}
        fold_rows.extend([train_metrics_row, val_metrics_row])
        
        # keep the original cv_thresholds append
        cv_thresholds.append(best_threshold)

        
        # Cleanup
        del model, trainer
        gc.collect()
        torch.cuda.empty_cache()

    folds_df = pd.DataFrame(fold_rows) if fold_rows else pd.DataFrame()
    summary_metrics = {}
    
    if not folds_df.empty:
        for s in ["train", "val"]:
            sub = folds_df[folds_df["set"] == s]
            if not sub.empty:
                # average numeric columns
                summary_metrics[s] = sub.drop(columns=["fold", "set"]).mean(numeric_only=True).to_dict()

    
    # Train final model
    print("\n" + "=" * 70)
    print("FINAL MODEL TRAINING")
    print("=" * 70)
    
    final_threshold = np.mean(cv_thresholds)
    print(f"Average threshold from CV: {final_threshold:.3f}")
    
    test_ws = folds[0]['test']
    train_ws = [w for w in watersheds if w not in test_ws]
    
    # Final normalization
    means, stds = compute_normalization_stats(dataset, all_cols, train_ws)
    
    # Final training dataset and loader
    final_train_dataset = MLPFloodDataset(
        dataset, train_ws, static_cols, temporal_cols,
        means, stds, config, mode='train'
    )
    # Final training loader
    final_loader = DataLoader(
        final_train_dataset,
        batch_size=config.batch_size,
        shuffle=True,
        num_workers=config.num_workers,
        pin_memory=(config.device == "cuda"),
        drop_last=True  # avoid batch size 1 in final training as well
    )
    
    # Final model and trainer
    final_model = MLPWithFiLM(len(static_cols), len(temporal_cols), config)
    final_trainer = Trainer(final_model, config)
    
    for epoch in range(config.epochs):
        train_loss = final_trainer.train_epoch(final_loader)
        print(f"Epoch {epoch+1}: Loss={train_loss:.4f}")
    
    # Build TEST dataset (normalized with final TRAIN stats) and evaluate
    test_dataset = MLPFloodDataset(
        dataset, test_ws, static_cols, temporal_cols,
        means, stds, config, mode='val'  # mode='val' to avoid any training-time sampling
    )
    test_metrics_row = evaluate_dataset(final_model, test_dataset, config.device, config.batch_size, final_threshold)
    summary_metrics["test"] = test_metrics_row
    
    # Package results
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
        'metrics_summary': summary_metrics  # NEW
    }
    
    if save_path:
        torch.save(results, save_path)
        print(f"Model saved to {save_path}")
    
    return results


# ======================== Visualization ========================

class MLPFloodVisualizer:
    """Generate visualizations matching XGBoost format."""
    
    def __init__(self, model_results: Dict, parquet_uri: str, training_folder: str):
        self.model_results = model_results  # NEW: store full results
        self.model = model_results['model']
        self.threshold = model_results['threshold']
        self.static_cols = model_results['static_cols']
        self.temporal_cols = model_results['temporal_cols']
        self.means = model_results['means']
        self.stds = model_results['stds']
        self.config = model_results['config']
        
        self.parquet_uri = parquet_uri
        self.training_folder = training_folder
        
        # Parse watershed info
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
        
        # Generate predictions
        predictions, confusions, per_date_metrics = self._generate_predictions()
        
        # Load reference TIFF for georeferencing
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
        
        dates = sorted(predictions.keys())
        
        # Create outputs
        pred_tiff = f"{output_folder}/mlp_predictions_{self.test_watershed}.tif"
        cm_tiff = f"{output_folder}/mlp_confusion_matrix_{self.test_watershed}.tif"
        
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
        """Generate predictions for test watershed - matching GLM metrics."""
        
        # Load parquet data for test watershed
        dataset = open_parquet_dataset(self.parquet_uri)
        
        filt = ds.field("Watershed_Loc") == self.test_watershed
        scanner = make_scanner(dataset, filter=filt)
        
        # Load all data
        pixel_data = []
        for batch in scanner_batches(scanner):
            if batch.num_rows > 0:
                df = pa.Table.from_batches([batch]).to_pandas()
                pixel_data.append(df)
        
        if not pixel_data:
            return {}, {}, {}
        
        df_all = pd.concat(pixel_data, ignore_index=True)
        
        # Get unique dates
        dates = sorted(df_all['Date'].unique())
        
        # Get image dimensions
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
            
            # Extract features
            X_static = date_data[self.static_cols].to_numpy(dtype=np.float32)
            X_temporal = date_data[self.temporal_cols].to_numpy(dtype=np.float32)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8)
            rows = date_data['PixelRow'].to_numpy()
            cols = date_data['PixelCol'].to_numpy()
            
            # Normalize
            static_idx = len(self.static_cols)
            X_static = np.where(np.isnan(X_static), self.means[:static_idx], X_static)
            X_static = (X_static - self.means[:static_idx]) / self.stds[:static_idx]
            
            X_temporal = np.where(np.isnan(X_temporal), self.means[static_idx:], X_temporal)
            X_temporal = (X_temporal - self.means[static_idx:]) / self.stds[static_idx:]
            
            # Predict in batches
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
            
            # Create images
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            pred_img[rows, cols] = y_pred
            predictions[date_str] = pred_img
            
            # Confusion matrix
            cm_img = np.full((H, W), 255, dtype=np.uint8)
            cm_vals = np.where(
                y_true == 0,
                np.where(y_pred == 0, 0, 3),  # TN=0, FP=3
                np.where(y_pred == 1, 1, 2)   # TP=1, FN=2
            )
            cm_img[rows, cols] = cm_vals
            confusions[date_str] = cm_img
            
            # Calculate metrics - MATCHING GLM EXACTLY
            
            # ROC AUC with error handling
            try:
                roc_auc = roc_auc_score(y_true, all_probs)
            except:
                roc_auc = np.nan
            
            # PR AUC with error handling (ADDED)
            try:
                pr_auc = average_precision_score(y_true, all_probs)
            except:
                pr_auc = np.nan
            
            # MCC with error handling (ADDED)
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
                'PR_AUC': float(pr_auc) if np.isfinite(pr_auc) else np.nan,  # ADDED
                'Kappa': cohen_kappa_score(y_true, y_pred),  # ADDED
                'MCC': float(mcc),  # ADDED
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
        """Create PNG visualizations (one per date), labeled with ROC AUC."""
        png_paths = []
    
        for date_str in sorted(dates):  # no limit; produce for all dates in chronological order
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

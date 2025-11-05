"""
tabnetFloods.py

TabNet architecture for flood prediction using pytorch-tabnet.
Rewritten with:
- Positive-aware watershed-level sampling (keep all positives; sample negatives to cap/fraction)
- Per-sample weights passed directly to TabNet fit()
- Strict scheduler param mapping by type (cosine vs step)
- Full PR-curve thresholding, rich metrics, and all-date PNGs
- TRAIN-only normalization per fold and final training (no leakage)
"""

import os
import re
import gc
import time
import warnings
import tempfile
from typing import List, Tuple, Dict, Optional, Any
from dataclasses import dataclass

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import torch

# TabNet imports
from pytorch_tabnet.tab_model import TabNetClassifier  # [[3]] [[5]] [[10]]
from pytorch_tabnet.augmentations import ClassificationSMOTE  # [[5]] [[3]]

from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, average_precision_score, cohen_kappa_score,
    matthews_corrcoef, precision_recall_curve, confusion_matrix
)

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
class TabNetConfig:
    # TabNet architecture parameters
    n_d: int = 64
    n_a: int = 64
    n_steps: int = 5
    gamma: float = 1.5
    n_independent: int = 2
    n_shared: int = 2
    epsilon: float = 1e-8

    # Virtual batch size for Ghost BatchNorm
    virtual_batch_size: int = 256
    momentum: float = 0.98

    # Attention and sparsity
    mask_type: str = "sparsemax"  # or "entmax"
    lambda_sparse: float = 1e-3

    # Training parameters
    epochs: int = 50
    batch_size: int = 2048
    learning_rate: float = 2e-2
    weight_decay: float = 1e-5

    # Scheduler (strict param mapping by type)
    scheduler_type: str = "cosine"  # "cosine" or "step"
    scheduler_step_size: int = 10   # for StepLR
    scheduler_gamma_lr: float = 0.9 # for StepLR
    scheduler_t0: int = 50          # for CosineAnnealingWarmRestarts
    scheduler_t_mult: int = 1       # for CosineAnnealingWarmRestarts

    # Early stopping
    patience: int = 10
    min_delta: float = 0.001

    # Class balancing
    use_class_weights: bool = True  # enable per-sample weights
    # (we compute per-sample weights from class frequencies; no SMOTE by default)

    # Data augmentation
    use_augmentation: bool = False
    aug_p: float = 0.2

    # Watershed-level sampling (uniform policy)
    sample_fraction: float = 0.2
    max_pixels_per_watershed: int = 200000
    enforce_at_watershed_level: bool = True
    positive_priority: bool = True  # keep all positives first, then sample negatives
    neg_floor_frac: float = 0.1  # NEW: ensure at least this fraction of negatives in samples

    # Fold strategy
    n_folds: int = 5
    train_frac: float = 0.6
    val_frac: float = 0.2
    force_test_ws: Optional[str] = None

    # System
    device: Optional[str] = None
    seed: int = 42
    num_workers: int = 4
    clip_value: float = 2.0  # gradient clipping

    # Metrics/logging
    compute_pr_auc_during_training: bool = True

    # Inference
    batch_size_inference: int = 50000

    # Visualization
    max_png_dates: Optional[int] = None  # None = process all dates

    def __post_init__(self):
        if self.device is None:
            self.device = "cuda" if torch.cuda.is_available() else "cpu"

    def get_scheduler_params(self) -> dict:
        """Return scheduler params based on scheduler type."""
        if self.scheduler_type == "cosine":
            return {"T_0": self.scheduler_t0, "T_mult": self.scheduler_t_mult}
        elif self.scheduler_type == "step":
            return {"step_size": self.scheduler_step_size, "gamma": self.scheduler_gamma_lr}
        else:
            return {}


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


def get_feature_columns(dataset) -> List[str]:
    """Get all numeric feature columns (static + temporal)."""
    exclude = {"Target", "Date", "Watershed_Loc", "PixelRow", "PixelCol"}
    feature_cols = []
    for field in dataset.schema:
        if field.name in exclude:
            continue
        t = field.type
        if pa.types.is_floating(t) or pa.types.is_integer(t):
            feature_cols.append(field.name)
    return feature_cols


def compute_normalization_stats(dataset, feature_cols: List[str],
                                watersheds: List[str], batch_size: int = 250000):
    """Optimized normalization computation with direct Arrow conversion."""
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

def _positive_aware_sample(ws_X: np.ndarray,
                           ws_y: np.ndarray,
                           sample_fraction: float,
                           cap: int,
                           neg_floor_frac: float = 0.1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Watershed-level sampling that keeps all positives (or a large portion),
    but ALWAYS includes some negatives to avoid all-positive batches.

    - desired_size = min(cap, floor(n * sample_fraction)) if either is set; else n.
    - If positives exceed desired_size, downsample positives to fit alongside a small
      negative floor (neg_floor_frac of desired_size, if available).
    - Otherwise, include all positives and fill the remaining slots with negatives.

    Args:
        ws_X: [N, F] features
        ws_y: [N] binary labels (0/1)
        sample_fraction: fraction of watershed rows to keep (<=1.0)
        cap: absolute cap on rows to keep (>0 to enable)
        neg_floor_frac: ensure at least this fraction of negatives in the sample (if available)

    Returns:
        (X_sel, y_sel) sampled arrays
    """
    n = len(ws_y)
    if n == 0:
        return ws_X, ws_y

    # Determine desired sample size
    desired = n
    if sample_fraction < 1.0:
        desired = int(max(1, np.floor(n * sample_fraction)))
    if cap > 0:
        desired = min(desired, cap)
    desired = max(1, desired)

    pos_idx = np.where(ws_y == 1)[0]
    neg_idx = np.where(ws_y == 0)[0]

    n_pos = len(pos_idx)
    n_neg = len(neg_idx)

    if n_neg == 0:
        # No negatives available; downsample positives if needed
        if n_pos > desired:
            sel_pos = np.random.choice(pos_idx, desired, replace=False)
            return ws_X[sel_pos], ws_y[sel_pos]
        else:
            return ws_X[pos_idx], ws_y[pos_idx]

    # There are some negatives; enforce a negative floor
    neg_floor = max(1, int(round(desired * neg_floor_frac)))
    neg_floor = min(neg_floor, n_neg)

    if n_pos >= desired:
        # Too many positives: reserve space for a small set of negatives
        pos_budget = max(1, desired - neg_floor)
        pos_budget = min(pos_budget, n_pos)
        sel_pos = np.random.choice(pos_idx, pos_budget, replace=False)
        sel_neg = np.random.choice(neg_idx, desired - pos_budget, replace=False)
        sel_idx = np.concatenate([sel_pos, sel_neg])
    else:
        # Keep all positives, fill the rest with negatives
        remaining = desired - n_pos
        neg_take = min(remaining, n_neg)
        # But ensure at least a small floor of negatives if possible
        neg_take = max(neg_take, min(neg_floor, n_neg))
        neg_take = min(neg_take, desired - min(desired - neg_floor, n_pos))  # limit to budget
        neg_take = min(neg_take, n_neg)
        sel_neg = np.random.choice(neg_idx, neg_take, replace=False)
        sel_idx = np.concatenate([pos_idx, sel_neg])

        # If still short due to rounding, top off with more negatives if available
        if len(sel_idx) < desired and n_neg > neg_take:
            top_off = min(desired - len(sel_idx), n_neg - neg_take)
            extra_neg = np.random.choice(np.setdiff1d(neg_idx, sel_neg, assume_unique=False), top_off, replace=False)
            sel_idx = np.concatenate([sel_idx, extra_neg])

        # If we somehow exceeded desired (shouldn't), trim randomly
        if len(sel_idx) > desired:
            sel_idx = np.random.choice(sel_idx, desired, replace=False)

    np.random.shuffle(sel_idx)
    return ws_X[sel_idx], ws_y[sel_idx]

def _build_tabnet_classifier(config: TabNetConfig, seed: int, verbose: int = 0):
    """
    Instantiate TabNetClassifier with best-effort compatibility across pytorch-tabnet versions.
    Some versions don't accept virtual_batch_size/momentum/clip_value in __init__.
    We try a rich signature first, then fall back to a minimal one.
    """
    scheduler_params = config.get_scheduler_params()
    scheduler_fn = (torch.optim.lr_scheduler.CosineAnnealingWarmRestarts
                    if config.scheduler_type == "cosine"
                    else torch.optim.lr_scheduler.StepLR)

    base_kwargs = dict(
        n_d=config.n_d,
        n_a=config.n_a,
        n_steps=config.n_steps,
        gamma=config.gamma,
        n_independent=config.n_independent,
        n_shared=config.n_shared,
        epsilon=config.epsilon,
        mask_type=config.mask_type,
        lambda_sparse=config.lambda_sparse,
        optimizer_fn=torch.optim.AdamW,
        optimizer_params={"lr": config.learning_rate, "weight_decay": config.weight_decay},
        scheduler_fn=scheduler_fn,
        scheduler_params=scheduler_params,
        seed=seed,
        verbose=verbose,
        device_name=config.device
    )

    # Try with virtual_batch_size/momentum/clip_value if supported
    try:
        clf = TabNetClassifier(
            virtual_batch_size=config.virtual_batch_size,
            momentum=config.momentum,
            clip_value=config.clip_value,
            **base_kwargs
        )
        return clf, True  # has_vbnbn
    except TypeError:
        pass

    # Fallback: try without virtual_batch_size and momentum/clip_value in constructor
    try:
        clf = TabNetClassifier(**base_kwargs)
        return clf, False  # no vbnbn in ctor
    except TypeError as e:
        # As a last resort: remove lambda_sparse/mask_type if also unsupported
        minimal_kwargs = base_kwargs.copy()
        minimal_kwargs.pop("mask_type", None)
        minimal_kwargs.pop("lambda_sparse", None)
        clf = TabNetClassifier(**minimal_kwargs)
        return clf, False

def load_data_for_watersheds(parquet_dataset, watersheds: List[str],
                             feature_cols: List[str], means: np.ndarray,
                             stds: np.ndarray, config: TabNetConfig,
                             mode: str = "train") -> Tuple[np.ndarray, np.ndarray]:
    """Optimized data loading with direct Arrow conversion."""
    X_parts = []
    y_parts = []

    print(f"  Loading {mode} data from {len(watersheds)} watersheds...")

    for ws in watersheds:
        filt = ds.field("Watershed_Loc") == ws
        cols = feature_cols + ["Target"]
        scanner = make_scanner(parquet_dataset, columns=cols, filter=filt, batch_size=500000)

        ws_X_batches = []
        ws_y_batches = []

        for batch in scanner_batches(scanner):
            if batch.num_rows == 0:
                continue
            
            X = arrow_batch_to_numpy(batch, feature_cols)
            y = batch.column("Target").to_numpy(zero_copy_only=False).astype(np.uint8)
            
            ws_X_batches.append(X)
            ws_y_batches.append(y)

        if ws_X_batches:
            ws_X = np.vstack(ws_X_batches)
            ws_y = np.hstack(ws_y_batches)

            # Apply watershed-level sampling only for TRAIN
            if mode == "train" and config.enforce_at_watershed_level:
                ws_X, ws_y = _positive_aware_sample(
                    ws_X, ws_y,
                    sample_fraction=config.sample_fraction,
                    cap=config.max_pixels_per_watershed,
                    neg_floor_frac=config.neg_floor_frac
                )

            X_parts.append(ws_X)
            y_parts.append(ws_y)

    if X_parts:
        X = np.vstack(X_parts)
        y = np.hstack(y_parts)

        # Normalize
        X = np.where(np.isnan(X), means, X)
        X = (X - means) / stds

        print(f"    Loaded {len(y):,} pixels ({y.mean():.4f} positive rate)")
        return X, y
    else:
        return np.empty((0, len(feature_cols)), dtype=np.float32), np.empty((0,), dtype=np.uint8)

# ======================== Metrics ========================

def compute_metrics(y_true: np.ndarray, y_proba: np.ndarray, threshold: float) -> Dict[str, float]:
    """Compute comprehensive metrics."""
    if len(y_true) == 0:
        return {
            "threshold": float(threshold), "prevalence": np.nan,
            "accuracy": np.nan, "precision": np.nan, "recall": np.nan, "f1": np.nan,
            "kappa": np.nan, "mcc": np.nan, "roc_auc": np.nan, "pr_auc": np.nan,
            "tp": np.nan, "fp": np.nan, "tn": np.nan, "fn": np.nan, "tpr": np.nan, "tnr": np.nan
        }

    y_pred = (y_proba >= threshold).astype(np.uint8)
    prevalence = float(np.mean(y_true))
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    kappa = cohen_kappa_score(y_true, y_pred)

    # MCC only when both classes present in y_true and y_pred
    mcc = 0.0
    if (len(np.unique(y_true)) > 1) and (len(np.unique(y_pred)) > 1):
        mcc = matthews_corrcoef(y_true, y_pred)

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


def find_optimal_threshold(y_true: np.ndarray, y_proba: np.ndarray) -> float:
    """Find F1-optimal threshold using PR curve (robust to class imbalance)."""
    if len(np.unique(y_true)) < 2:
        return 0.5
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_proba)
    if thresholds.size == 0:
        return 0.5
    f1_scores = 2 * precisions[:-1] * recalls[:-1] / np.clip(precisions[:-1] + recalls[:-1], 1e-8, None)
    best_idx = int(np.nanargmax(f1_scores))
    return float(thresholds[best_idx]) if best_idx < thresholds.size else 0.5


# ======================== CV Splits ========================

def spatial_kfold_split(watersheds: List[str], config: TabNetConfig):
    """Create k-fold splits with a fixed TEST set and repeated random Train/Val."""
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

        folds.append({"train": train_ws, "val": val_ws, "test": test_ws})

    return folds


# ======================== Training ========================

def _compute_per_sample_weights(y_train: np.ndarray) -> np.ndarray:
    """
    Compute per-sample weights for binary classification:
    balance classes so each contributes equally to the loss.
    """
    n_pos = max(1, int((y_train == 1).sum()))
    n_neg = max(1, int((y_train == 0).sum()))
    total = n_pos + n_neg
    w0 = total / (2.0 * n_neg)
    w1 = total / (2.0 * n_pos)

    weights = np.ones(len(y_train), dtype=np.float32)
    weights[y_train == 0] = w0
    weights[y_train == 1] = w1
    return weights


def train_tabnet_fold(X_train, y_train, X_val, y_val, config: TabNetConfig, fold_idx: int):
    """Train TabNet for one fold with per-sample weights and robust constructor."""

    # Per-sample weights (balanced by class frequency)
    if config.use_class_weights:
        n_pos = max(1, int((y_train == 1).sum()))
        n_neg = max(1, int((y_train == 0).sum()))
        total = n_pos + n_neg
        w0 = total / (2.0 * n_neg)
        w1 = total / (2.0 * n_pos)
        sample_weights = np.where(y_train == 1, w1, w0).astype(np.float32)
    else:
        sample_weights = None

    clf, has_vbnbn = _build_tabnet_classifier(config, seed=config.seed + fold_idx, verbose=0)

    aug = ClassificationSMOTE(p=config.aug_p, seed=config.seed + fold_idx) if config.use_augmentation else None

    clf.fit(
        X_train=X_train,
        y_train=y_train,
        eval_set=[(X_val, y_val)],
        eval_metric=['auc'],
        max_epochs=config.epochs,
        patience=config.patience,
        batch_size=config.batch_size,
        # Pass virtual batch size during fit (works across versions)
        virtual_batch_size=config.virtual_batch_size if not has_vbnbn else config.virtual_batch_size,
        num_workers=config.num_workers,
        weights=sample_weights,  # per-sample weights (or None)
        drop_last=False,
        augmentations=aug
    )

    if config.compute_pr_auc_during_training:
        y_val_proba = clf.predict_proba(X_val)[:, 1]
        try:
            pr_auc_val = average_precision_score(y_val, y_val_proba)
            print(f"    Validation PR-AUC: {pr_auc_val:.4f}")
        except Exception:
            pass

    return clf

def train_tabnet_flood_model(parquet_uri: str, config: Optional[TabNetConfig] = None,
                             save_path: Optional[str] = None) -> Dict:
    """Train TabNet flood model with k-fold CV and produce a trained final model."""
    config = config or TabNetConfig()
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)

    print("=" * 70)
    print("TABNET FLOOD MODEL - K-FOLD SPATIAL CV")
    print("=" * 70)

    # Load dataset and metadata
    dataset = open_parquet_dataset(parquet_uri)
    watersheds = get_unique_watersheds(dataset)
    feature_cols = get_feature_columns(dataset)

    print(f"Watersheds: {len(watersheds)}")
    print(f"Features: {len(feature_cols)}")

    folds = spatial_kfold_split(watersheds, config)
    print(f"Test watersheds: {folds[0]['test']}")

    cv_thresholds = []
    fold_rows = []

    # Cross-validation loop
    for fold_idx, fold in enumerate(folds):
        print(f"\nFold {fold_idx + 1}/{len(folds)}")
        print(f"  Train: {len(fold['train'])} watersheds")
        print(f"  Val: {len(fold['val'])} watersheds")

        # TRAIN-only normalization
        means, stds = compute_normalization_stats(dataset, feature_cols, fold['train'])

        # Load data (watershed-level, positive-aware sampling)
        X_train, y_train = load_data_for_watersheds(dataset, fold['train'], feature_cols, means, stds, config, mode='train')
        X_val, y_val = load_data_for_watersheds(dataset, fold['val'], feature_cols, means, stds, config, mode='val')

        if len(X_train) == 0 or len(X_val) == 0:
            print("  Warning: empty dataset")
            cv_thresholds.append(0.5)
            continue

        # Train
        print(f"  Training TabNet (pos rate: train={y_train.mean():.4f}, val={y_val.mean():.4f})")
        clf = train_tabnet_fold(X_train, y_train, X_val, y_val, config, fold_idx)

        # Predictions and best threshold (F1-optimal)
        y_train_proba = clf.predict_proba(X_train)[:, 1]
        y_val_proba = clf.predict_proba(X_val)[:, 1]
        best_threshold = find_optimal_threshold(y_val, y_val_proba)
        cv_thresholds.append(best_threshold)

        # Metrics
        train_metrics = compute_metrics(y_train, y_train_proba, best_threshold)
        val_metrics = compute_metrics(y_val, y_val_proba, best_threshold)
        print(f"  Val F1={val_metrics['f1']:.4f}, PR-AUC={val_metrics['pr_auc']:.4f}, Threshold={best_threshold:.3f}")

        # Store metrics
        train_metrics['fold'] = fold_idx + 1
        train_metrics['set'] = 'train'
        val_metrics['fold'] = fold_idx + 1
        val_metrics['set'] = 'val'
        fold_rows.extend([train_metrics, val_metrics])

        # Cleanup
        del clf, X_train, y_train, X_val, y_val
        gc.collect()
        if config.device == "cuda":
            torch.cuda.empty_cache()

    # Aggregate CV metrics
    folds_df = pd.DataFrame(fold_rows) if fold_rows else pd.DataFrame()
    summary_metrics = {}
    if not folds_df.empty:
        for s in ["train", "val"]:
            sub = folds_df[folds_df["set"] == s]
            if not sub.empty:
                summary_metrics[s] = sub.drop(columns=["fold", "set"]).mean(numeric_only=True).to_dict()

    # Final training
    print("\n" + "=" * 70)
    print("FINAL MODEL TRAINING")
    print("=" * 70)

    final_threshold = float(np.mean(cv_thresholds)) if cv_thresholds else 0.5
    print(f"Average threshold from CV: {final_threshold:.3f}")

    test_ws = folds[0]['test']
    train_ws = [w for w in watersheds if w not in test_ws]

    # Final normalization (all non-TEST watersheds)
    means, stds = compute_normalization_stats(dataset, feature_cols, train_ws)

    # Final training data
    X_train_final, y_train_final = load_data_for_watersheds(dataset, train_ws, feature_cols, means, stds, config, mode='train')
    print(f"Final training on {len(X_train_final):,} pixels")

    # Final per-sample weights (if enabled)
    final_sample_weights = _compute_per_sample_weights(y_train_final) if config.use_class_weights else None

    scheduler_params = config.get_scheduler_params()
    scheduler_fn = (torch.optim.lr_scheduler.CosineAnnealingWarmRestarts
                    if config.scheduler_type == "cosine"
                    else torch.optim.lr_scheduler.StepLR)

    # Build final classifier (version-tolerant)
    final_clf, has_vbnbn = _build_tabnet_classifier(config, seed=config.seed, verbose=1)

    # Per-sample weights for final training (optional)
    if config.use_class_weights:
        n_pos = max(1, int((y_train_final == 1).sum()))
        n_neg = max(1, int((y_train_final == 0).sum()))
        total = n_pos + n_neg
        w0 = total / (2.0 * n_neg)
        w1 = total / (2.0 * n_pos)
        final_sample_weights = np.where(y_train_final == 1, w1, w0).astype(np.float32)
    else:
        final_sample_weights = None

    final_clf.fit(
        X_train=X_train_final,
        y_train=y_train_final,
        max_epochs=config.epochs,
        batch_size=config.batch_size,
        virtual_batch_size=config.virtual_batch_size if not has_vbnbn else config.virtual_batch_size,
        num_workers=config.num_workers,
        weights=final_sample_weights,
        drop_last=False
    )

    # Test evaluation
    if test_ws:
        X_test, y_test = load_data_for_watersheds(dataset, test_ws, feature_cols, means, stds, config, mode='test')
        if len(X_test) > 0:
            y_test_proba = final_clf.predict_proba(X_test)[:, 1]
            test_metrics = compute_metrics(y_test, y_test_proba, final_threshold)
            summary_metrics['test'] = test_metrics
            print(f"Test F1={test_metrics['f1']:.4f}, PR-AUC={test_metrics['pr_auc']:.4f}")

    # Feature importance (TabNet mask importances) [[5]] [[10]]
    feature_importance_df = pd.DataFrame({
        'feature': feature_cols,
        'importance': final_clf.feature_importances_
    }).sort_values('importance', ascending=False)

    # Package results
    results = {
        'model': final_clf,
        'threshold': final_threshold,
        'feature_cols': feature_cols,
        'means': means,
        'stds': stds,
        'config': config,
        'test_watersheds': test_ws,
        'metrics_summary': summary_metrics,
        'feature_importance': feature_importance_df
    }

    if save_path:
        torch.save(results, save_path)
        print(f"Model saved to {save_path}")

    return results


# ======================== Visualization ========================

class TabNetFloodVisualizer:
    """Generate visualizations matching XGBoost/MLP format."""

    def __init__(self, model_results: Dict, parquet_uri: str, training_folder: str):
        self.model_results = model_results
        self.model = model_results['model']
        self.threshold = model_results['threshold']
        self.feature_cols = model_results['feature_cols']
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
        else:
            self.test_lat = None
            self.test_lon = None
            self.test_watershed = None

        self.s3_client = boto3.client('s3')
        configure_gdal()

    def generate_all_outputs(self, output_folder: Optional[str] = None):
        """Generate all visualization outputs."""
        if self.test_watershed is None:
            raise ValueError("Could not parse test watershed from training_folder path.")

        if output_folder is None:
            output_folder = f"{self.training_folder}/assessTabNetModel"

        print("\n" + "="*60)
        print("GENERATING TABNET MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        print(f"Test watershed: {self.test_watershed}")
        print(f"Output folder: {output_folder}")

        # Generate predictions
        predictions, confusions, per_date_metrics = self._generate_predictions()

        # Reference TIFF for georeferencing
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)

        dates = sorted(predictions.keys())

        # Write TIFFs
        pred_tiff = f"{output_folder}/tabnet_predictions_{self.test_watershed}.tif"
        cm_tiff = f"{output_folder}/tabnet_confusion_matrix_{self.test_watershed}.tif"
        self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "TabNet_Prediction")
        self._write_multiband_tiff(ref_ds, confusions, dates, cm_tiff, "TabNet_CM")

        # PNGs and CSVs
        png_paths = self._create_pngs(confusions, dates, per_date_metrics, output_folder)
        csv_paths = self._write_csvs(output_folder, per_date_metrics)

        print("\nTabNet visualization complete!")
        return {
            'prediction_tiff': pred_tiff,
            'confusion_matrix_tiff': cm_tiff,
            'png_visualizations': png_paths,
            'csv_outputs': csv_paths
        }

    def _generate_predictions(self):
        """Generate predictions for the test watershed."""
        dataset = open_parquet_dataset(self.parquet_uri)
        filt = ds.field("Watershed_Loc") == self.test_watershed
        scanner = make_scanner(dataset, filter=filt)

        # Load all rows for this watershed
        pixel_data = []
        for batch in scanner_batches(scanner):
            if batch.num_rows > 0:
                df = pa.Table.from_batches([batch]).to_pandas()
                pixel_data.append(df)
        if not pixel_data:
            return {}, {}, {}

        df_all = pd.concat(pixel_data, ignore_index=True)
        dates = sorted(df_all['Date'].unique())

        H = int(df_all['PixelRow'].max()) + 1
        W = int(df_all['PixelCol'].max()) + 1

        predictions = {}
        confusions = {}
        per_date_metrics = {}

        print(f"  Processing {len(dates)} dates...")

        for i, date_str in enumerate(dates):
            if (i + 1) % 10 == 0:
                print(f"    Date {i + 1}/{len(dates)}")

            date_data = df_all[df_all['Date'] == date_str]

            X = date_data[self.feature_cols].to_numpy(dtype=np.float32)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8)
            rows = date_data['PixelRow'].to_numpy()
            cols = date_data['PixelCol'].to_numpy()

            # Normalize
            X = np.where(np.isnan(X), self.means, X)
            X = (X - self.means) / self.stds

            # Predict
            y_proba = self.model.predict_proba(X)[:, 1]
            y_pred = (y_proba >= self.threshold).astype(np.uint8)

            # Prediction image
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            pred_img[rows, cols] = y_pred
            predictions[date_str] = pred_img

            # Confusion image (0 TN, 1 TP, 2 FN, 3 FP, 255 nodata)
            cm_img = np.full((H, W), 255, dtype=np.uint8)
            cm_vals = np.where(
                y_true == 0,
                np.where(y_pred == 0, 0, 3),
                np.where(y_pred == 1, 1, 2)
            )
            cm_img[rows, cols] = cm_vals
            confusions[date_str] = cm_img

            # Per-date metrics
            try:
                auc = roc_auc_score(y_true, y_proba)
            except Exception:
                auc = np.nan
            try:
                pr_auc = average_precision_score(y_true, y_proba)
            except Exception:
                pr_auc = np.nan

            per_date_metrics[date_str] = {
                'date': date_str,
                'F1': f1_score(y_true, y_pred, zero_division=0),
                'Accuracy': accuracy_score(y_true, y_pred),
                'Precision': precision_score(y_true, y_pred, zero_division=0),
                'Recall': recall_score(y_true, y_pred, zero_division=0),
                'ROC_AUC': auc,
                'PR_AUC': pr_auc,
                'Kappa': cohen_kappa_score(y_true, y_pred),
                'MCC': matthews_corrcoef(y_true, y_pred) if (len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1) else 0.0,
                'actual_pos_frac': y_true.mean(),
                'pred_pos_frac': y_pred.mean()
            }

        return predictions, confusions, per_date_metrics

    def _write_multiband_tiff(self, ref_ds, arrays, dates, output_path, prefix):
        """Write multiband GeoTIFF, one band per date."""
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
            temp_path = tmp.name

        driver = gdal.GetDriverByName('GTiff')
        out = driver.Create(
            temp_path, ref_ds.RasterXSize, ref_ds.RasterYSize, len(dates),
            gdal.GDT_Byte, options=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"]
        )

        out.SetGeoTransform(ref_ds.GetGeoTransform())
        out.SetProjection(ref_ds.GetProjection())

        for i, date_str in enumerate(dates, 1):
            band = out.GetRasterBand(i)
            band.SetNoDataValue(255)
            band.SetDescription(f"{prefix}_{date_str}")
            arr = arrays.get(date_str, np.full((ref_ds.RasterYSize, ref_ds.RasterXSize), 255, dtype=np.uint8))
            band.WriteArray(arr)

        out = None

        # Upload to S3
        bucket, key = output_path.replace('s3://', '').split('/', 1)
        boto3.client('s3').upload_file(temp_path, bucket, key)
        os.remove(temp_path)
        print(f"  Created: {output_path}")

    def _create_pngs(self, confusions, dates, metrics, output_folder):
        """Create PNG visualizations."""
        png_paths = []

        dates_to_process = sorted(dates)
        if self.config.max_png_dates is not None:
            dates_to_process = dates_to_process[:self.config.max_png_dates]

        for date_str in dates_to_process:
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
                mpatches.Patch(color='green', label='True Negative (Land correct)'),
                mpatches.Patch(color='blue', label='True Positive (Water correct)'),
                mpatches.Patch(color='orange', label='False Negative (Missed water)'),
                mpatches.Patch(color='red', label='False Positive (False water)')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
            ax.set_title(f'TabNet Confusion Matrix - {date_str}', fontsize=14, fontweight='bold')

            auc_val = m.get('ROC_AUC', np.nan)
            pr_auc_val = m.get('PR_AUC', np.nan)
            auc_str = f"{auc_val:.3f}" if np.isfinite(auc_val) else "N/A"
            pr_auc_str = f"{pr_auc_val:.3f}" if np.isfinite(pr_auc_val) else "N/A"

            text = (
                f"Date: {date_str}\n"
                f"Model: TabNet\n"
                f"F1 Score: {m.get('F1', 0):.3f}\n"
                f"Accuracy: {m.get('Accuracy', 0):.3f}\n"
                f"Precision: {m.get('Precision', 0):.3f}\n"
                f"Recall: {m.get('Recall', 0):.3f}\n"
                f"ROC AUC: {auc_str}\n"
                f"PR AUC: {pr_auc_str}\n"
                f"Kappa: {m.get('Kappa', 0):.3f}\n"
                f"MCC: {m.get('MCC', 0):.3f}"
            )
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11, va='top', bbox=props)

            ax.axis('off')
            plt.tight_layout()

            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()

            s3_path = f"{output_folder}/tabnet_confusion_matrix_{date_str}.png"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(temp_png, bucket, key)
            os.remove(temp_png)
            png_paths.append(s3_path)

        print(f"  Created {len(png_paths)} PNG visualizations")
        return png_paths

    def _write_csvs(self, output_folder, per_date_metrics):
        """Write CSV outputs."""
        paths = {}

        # Per-date metrics
        if per_date_metrics:
            df = pd.DataFrame.from_dict(per_date_metrics, orient='index')
            df = df.sort_values("date") if "date" in df.columns else df
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
            df.to_csv(tmp_csv, index=False)
            s3_path = f"{output_folder}/tabnet_per_date_metrics_{self.test_watershed}.csv"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(tmp_csv, bucket, key)
            os.remove(tmp_csv)
            paths['per_date_metrics'] = s3_path

        # Summary metrics
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
                s3_path = f"{output_folder}/tabnet_metrics_summary.csv"
                bucket, key = s3_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(tmp_csv, bucket, key)
                os.remove(tmp_csv)
                paths['metrics_summary'] = s3_path

        # Feature importance
        feature_importance = self.model_results.get('feature_importance')
        if feature_importance is not None:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
            feature_importance.to_csv(tmp_csv, index=False)
            s3_path = f"{output_folder}/tabnet_feature_importances.csv"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(tmp_csv, bucket, key)
            os.remove(tmp_csv)
            paths['feature_importances'] = s3_path

        return paths


# ======================== Main Entry Point ========================

if __name__ == "__main__":
    # Example configuration
    config = TabNetConfig(
        n_d=64, n_a=64, n_steps=5, gamma=1.5,
        lambda_sparse=1e-3,
        epochs=50,
        batch_size=2048,
        learning_rate=2e-2,
        patience=10,
        use_class_weights=True,
        use_augmentation=False,  # SMOTE disabled by default
        sample_fraction=0.2,
        enforce_at_watershed_level=True,
        positive_priority=True,
        compute_pr_auc_during_training=True,
        max_png_dates=None
    )

    parquet_uri = "s3://your-bucket/flood-data.parquet"
    results = train_tabnet_flood_model(parquet_uri, config)

    training_folder = "s3://your-bucket/trainingFolderFor_Lat_XX_Lon_YY"
    viz = TabNetFloodVisualizer(results, parquet_uri, training_folder)
    outputs = viz.generate_all_outputs()

    print("\nTabNet training and visualization complete!")
    print("Generated outputs:")
    for key, path in outputs.items():
        if isinstance(path, list):
            print(f"  {key}: {len(path)} files")
        else:
            print(f"  {key}: {path}")



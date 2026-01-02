"""
xgboostFloods.py - FULLY OPTIMIZED VERSION

XGBoost flood prediction with all GLM improvements applied.
Uses proper sampling, hyperparameter tuning, and comprehensive diagnostics.

Key features:
- Adaptive target-based sampling for cross-watershed balance
- Optional hyperparameter tuning (RandomizedSearchCV)
- Sampling-only class balancing (no scale_pos_weight)
- Streaming data loading for memory efficiency
- MCC-based threshold optimization
- Comprehensive diagnostic plots
- Feature importance analysis
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
from collections import defaultdict

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import xgboost as xgb
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from scipy.stats import randint, uniform, loguniform

from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, average_precision_score, cohen_kappa_score,
    matthews_corrcoef, precision_recall_curve, confusion_matrix,
    make_scorer
)

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


# ======================== Configuration ========================

@dataclass
class XGBoostConfig:
    """Configuration for XGBoost flood model."""
    
    # XGBoost parameters (can be tuned)
    n_estimators: int = 320 #300
    max_depth: int = 8 #6
    learning_rate: float = 0.1
    subsample: float = 0.8
    colsample_bytree: float = 0.8
    min_child_weight: int = 3
    gamma: float = 0.1
    reg_alpha: float = 0.05 #0.01
    reg_lambda: float = 1.0
    
    # Data sampling - ADAPTIVE TARGET-BASED
    use_adaptive_sampling: bool = True
    min_target_rate: float = 0.05
    max_target_rate: float = 0.50
    max_pixels_per_watershed: int = 200000
    min_samples_per_class: int = 20
    
    # Hyperparameter tuning
    tune_hyperparameters: bool = False
    n_tuning_iterations: int = 20
    tuning_cv_folds: int = 3
    
    # Fold strategy
    n_folds: int = 5
    force_test_ws: Optional[str] = None
    
    # System
    seed: int = 42
    n_jobs: int = -1
    
    # Streaming
    batch_size_loading: int = 250000
    
    # XGBoost specific
    tree_method: str = "hist"
    objective: str = "binary:logistic"
    eval_metric: str = "logloss"
    early_stopping_rounds: int = 20
    
    # Visualization
    max_png_dates: Optional[int] = None


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


def arrow_batch_to_numpy(batch, columns: List[str]) -> np.ndarray:
    """Efficiently convert Arrow batch to NumPy array."""
    if batch.num_rows == 0:
        return np.empty((0, len(columns)), dtype=np.float32)
    
    arrays = []
    for i, col_name in enumerate(columns):
        arr = batch.column(i).to_numpy(zero_copy_only=False)
        if arr.dtype != np.float32:
            arr = arr.astype(np.float32, copy=False)
        arrays.append(arr)
    
    return np.column_stack(arrays) if arrays else np.empty((batch.num_rows, 0), dtype=np.float32)


# ======================== Data Utilities ========================

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
    """Get all feature columns (both static and temporal)."""
    exclude = {"Target", "Date", "Watershed_Loc", "PixelRow", "PixelCol"}
    
    feature_cols = []
    for field in dataset.schema:
        if field.name in exclude:
            continue
        t = field.type
        if pa.types.is_floating(t) or pa.types.is_integer(t):
            feature_cols.append(field.name)
    
    return feature_cols


def compute_watershed_prevalence(parquet_dataset, watersheds: List[str], 
                                batch_size: int = 250000) -> Dict[str, float]:
    """Compute natural positive prevalence for each watershed."""
    print(f"  Computing natural prevalence for {len(watersheds)} watersheds...")
    
    prevalence_map = {}
    
    for ws in watersheds:
        filt = ds.field("Watershed_Loc") == ws
        scanner = make_scanner(parquet_dataset, columns=["Target"], 
                             filter=filt, batch_size=batch_size)
        
        pos_count = 0
        total_count = 0
        
        for batch in scanner_batches(scanner):
            if batch.num_rows == 0:
                continue
            
            y_batch = batch.column(0).to_numpy(zero_copy_only=False).astype(np.uint8)
            pos_count += np.sum(y_batch == 1)
            total_count += len(y_batch)
        
        prevalence = pos_count / total_count if total_count > 0 else 0.0
        prevalence_map[ws] = prevalence
        print(f"    {ws}: {prevalence:.4f} ({pos_count:,}/{total_count:,})")
    
    return prevalence_map


def compute_adaptive_target_rate_sqrt(natural_prevalence: float, 
                                     min_rate: float = 0.05,
                                     max_rate: float = 0.50) -> float:
    """Compute adaptive target rate using square root compression."""
    sqrt_prev = np.sqrt(natural_prevalence)
    sqrt_min = np.sqrt(0.01)
    sqrt_max = np.sqrt(0.50)
    
    normalized = (sqrt_prev - sqrt_min) / (sqrt_max - sqrt_min)
    normalized = np.clip(normalized, 0, 1)
    
    target = min_rate + normalized * (max_rate - min_rate)
    
    return float(target)


# ======================== Sampling Strategy ========================

def balanced_watershed_sample(X: np.ndarray, y: np.ndarray,
                              target_pos_rate: float = 0.25,
                              max_samples: int = 200000,
                              min_samples_per_class: int = 20) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """Sample to achieve TARGET positive rate per watershed."""
    n = len(y)
    if n == 0:
        return None, None
    
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    
    n_pos = len(pos_idx)
    n_neg = len(neg_idx)
    
    if n_pos < min_samples_per_class or n_neg < min_samples_per_class:
        return None, None
    
    n_neg_target = int(n_pos * (1 - target_pos_rate) / target_pos_rate)
    n_neg_sample = min(n_neg_target, n_neg)
    
    total = n_pos + n_neg_sample
    
    if total > max_samples:
        scale = max_samples / total
        n_pos_sample = max(min_samples_per_class, int(n_pos * scale))
        n_neg_sample = max_samples - n_pos_sample
        
        if n_pos_sample < min_samples_per_class or n_neg_sample < min_samples_per_class:
            return None, None
        
        sel_pos = np.random.choice(pos_idx, n_pos_sample, replace=False)
        sel_neg = np.random.choice(neg_idx, n_neg_sample, replace=False)
    else:
        sel_pos = pos_idx
        sel_neg = np.random.choice(neg_idx, n_neg_sample, replace=False)
    
    sel_idx = np.concatenate([sel_pos, sel_neg])
    np.random.shuffle(sel_idx)
    
    return X[sel_idx], y[sel_idx]


# ======================== Data Loading - OPTIMIZED ========================

def compute_normalization_stats(dataset, feature_cols: List[str], 
                               watersheds: List[str], batch_size: int = 250000):
    """Compute mean and std for normalization using streaming."""
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


def load_data_streaming(parquet_dataset, watersheds: List[str],
                       feature_cols: List[str], means: np.ndarray, 
                       stds: np.ndarray, config: XGBoostConfig,
                       mode: str = "train",
                       prevalence_map: Optional[Dict[str, float]] = None) -> Tuple[np.ndarray, np.ndarray]:
    """OPTIMIZED: Single-pass data loading with adaptive streaming sampling."""
    
    print(f"  Loading {mode} data from {len(watersheds)} watersheds...")
    
    filt = ds.field("Watershed_Loc").isin(pa.array(watersheds))
    cols = feature_cols + ["Target", "Watershed_Loc"]
    scanner = make_scanner(parquet_dataset, columns=cols, filter=filt, 
                          batch_size=config.batch_size_loading)
    
    ws_data = defaultdict(lambda: {'X': [], 'y': []})
    
    for batch in scanner_batches(scanner):
        if batch.num_rows == 0:
            continue
        
        X_batch = arrow_batch_to_numpy(batch, feature_cols)
        
        target_idx = len(feature_cols)
        ws_idx = len(feature_cols) + 1
        y_batch = batch.column(target_idx).to_numpy(zero_copy_only=False).astype(np.uint8)
        ws_batch = batch.column(ws_idx).to_pandas().values
        
        for ws in np.unique(ws_batch):
            mask = ws_batch == ws
            ws_data[ws]['X'].append(X_batch[mask])
            ws_data[ws]['y'].append(y_batch[mask])
    
    X_parts = []
    y_parts = []
    skipped = []
    
    for ws in watersheds:
        if ws not in ws_data:
            skipped.append(ws)
            continue
        
        X_ws = np.vstack(ws_data[ws]['X'])
        y_ws = np.hstack(ws_data[ws]['y'])
        
        if mode == "train":
            if config.use_adaptive_sampling and prevalence_map is not None:
                natural_prev = prevalence_map.get(ws, 0.25)
                target_rate = compute_adaptive_target_rate_sqrt(
                    natural_prev,
                    min_rate=config.min_target_rate,
                    max_rate=config.max_target_rate
                )
                boost = target_rate / natural_prev if natural_prev > 0 else 1.0
                print(f"    {ws}: natural={natural_prev:.4f} → target={target_rate:.4f} (boost={boost:.2f}x)")
            else:
                target_rate = (config.min_target_rate + config.max_target_rate) / 2
            
            X_sampled, y_sampled = balanced_watershed_sample(
                X_ws, y_ws,
                target_pos_rate=target_rate,
                max_samples=config.max_pixels_per_watershed,
                min_samples_per_class=config.min_samples_per_class
            )
            
            if X_sampled is None:
                skipped.append(ws)
                print(f"    Skipped {ws}: insufficient samples")
                continue
            
            X_ws = X_sampled
            y_ws = y_sampled
        
        X_parts.append(X_ws)
        y_parts.append(y_ws)
    
    if skipped:
        print(f"    Skipped {len(skipped)} watersheds with insufficient data")
    
    if not X_parts:
        return np.empty((0, len(feature_cols)), dtype=np.float32), np.empty((0,), dtype=np.uint8)
    
    X = np.vstack(X_parts)
    y = np.hstack(y_parts)
    
    # Normalize (XGBoost benefits from normalization for regularization)
    X = np.where(np.isnan(X), means, X)
    X = (X - means) / stds
    
    print(f"    Loaded {len(y):,} pixels ({y.mean():.4f} positive rate)")
    
    return X, y


# ======================== Metrics - MCC OPTIMIZED ========================

def compute_metrics(y_true: np.ndarray, y_proba: np.ndarray, threshold: float) -> Dict[str, float]:
    """Compute comprehensive metrics."""
    y_pred = (y_proba >= threshold).astype(np.uint8)
    
    prevalence = float(np.mean(y_true))
    accuracy = accuracy_score(y_true, y_pred)
    precision = precision_score(y_true, y_pred, zero_division=0)
    recall = recall_score(y_true, y_pred, zero_division=0)
    f1 = f1_score(y_true, y_pred, zero_division=0)
    kappa = cohen_kappa_score(y_true, y_pred)
    
    if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1:
        mcc = matthews_corrcoef(y_true, y_pred)
    else:
        mcc = 0.0
    
    try:
        roc_auc = roc_auc_score(y_true, y_proba)
    except:
        roc_auc = np.nan
    try:
        pr_auc = average_precision_score(y_true, y_proba)
    except:
        pr_auc = np.nan
    
    if len(np.unique(y_true)) > 1:
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    else:
        if y_true[0] == 0:
            tn, fp, fn, tp = len(y_true), 0, 0, 0
        else:
            tn, fp, fn, tp = 0, 0, 0, len(y_true)
    
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


def find_optimal_threshold_mcc(y_true: np.ndarray, y_proba: np.ndarray) -> float:
    """Find MCC-optimal threshold."""
    if len(np.unique(y_true)) < 2:
        return 0.5
    
    thresholds = np.linspace(0.01, 0.99, 100)
    best_mcc = -1
    best_thresh = 0.5
    
    for thresh in thresholds:
        y_pred = (y_proba >= thresh).astype(int)
        if len(np.unique(y_pred)) < 2:
            continue
        try:
            mcc = matthews_corrcoef(y_true, y_pred)
            if mcc > best_mcc:
                best_mcc = mcc
                best_thresh = thresh
        except:
            continue
    
    return float(best_thresh)


# ======================== Hyperparameter Tuning ========================

def get_hyperparameter_grid():
    """Define hyperparameter distributions for RandomizedSearch."""
    return {
        'n_estimators': randint(100, 500),
        'max_depth': randint(3, 10),
        'learning_rate': loguniform(0.01, 0.3),
        'subsample': uniform(0.6, 0.4),  # 0.6 to 1.0
        'colsample_bytree': uniform(0.6, 0.4),  # 0.6 to 1.0
        'min_child_weight': randint(1, 7),
        'gamma': loguniform(0.01, 1.0),
        'reg_alpha': loguniform(0.001, 1.0),
        'reg_lambda': loguniform(0.1, 10.0)
    }


def tune_hyperparameters(X_train: np.ndarray, y_train: np.ndarray, 
                        config: XGBoostConfig) -> Dict[str, Any]:
    """Perform randomized hyperparameter search."""
    print("  Performing hyperparameter tuning...")
    print(f"    Trials: {config.n_tuning_iterations}")
    print(f"    Inner CV folds: {config.tuning_cv_folds}")
    
    base_model = xgb.XGBClassifier(
        objective=config.objective,
        tree_method=config.tree_method,
        random_state=config.seed,
        n_jobs=1,
        eval_metric=config.eval_metric
    )
    
    param_distributions = get_hyperparameter_grid()
    mcc_scorer = make_scorer(matthews_corrcoef)
    
    search = RandomizedSearchCV(
        estimator=base_model,
        param_distributions=param_distributions,
        n_iter=config.n_tuning_iterations,
        cv=StratifiedKFold(n_splits=config.tuning_cv_folds, shuffle=True, random_state=config.seed),
        scoring=mcc_scorer,
        n_jobs=config.n_jobs,
        random_state=config.seed,
        verbose=1
    )
    
    search.fit(X_train, y_train)
    
    print(f"    Best MCC: {search.best_score_:.4f}")
    print(f"    Best params: {search.best_params_}")
    
    return search.best_params_


# ======================== Training Pipeline ========================

def spatial_kfold_split(watersheds: List[str], config: XGBoostConfig):
    """Create k-fold spatial splits for leave-one-watershed-out CV."""
    
    rng = np.random.RandomState(config.seed)
    
    if config.force_test_ws is not None and config.force_test_ws in watersheds:
        test_ws = [config.force_test_ws]
    else:
        if config.force_test_ws is not None:
            print(f"Warning: force_test_ws '{config.force_test_ws}' not found in watersheds")
        test_ws = [rng.choice(watersheds)]
        print(f"Randomly selected test watershed: {test_ws[0]}")
    
    remaining = [w for w in watersheds if w not in test_ws]
    
    shuffled = remaining[:]
    rng.shuffle(shuffled)
    
    folds = []
    for fold_idx in range(config.n_folds):
        val_ws = [shuffled[fold_idx]]
        train_ws = [w for w in shuffled if w != val_ws[0]]
        
        folds.append({
            'train': train_ws,
            'val': val_ws,
            'test': test_ws
        })
    
    return folds


def train_xgboost_flood_model(parquet_uri: str, training_folder: str, 
                              config: Optional[XGBoostConfig] = None) -> Dict:
    """Train XGBoost model with k-fold CV and adaptive sampling.
    
    Args:
        parquet_uri: Path to parquet file
        training_folder: S3 folder path for saving outputs (e.g., s3://bucket/path/trainingFolderFor_Lat_X_Lon_Y)
        config: Optional configuration object
    """
    
    config = config or XGBoostConfig()
    np.random.seed(config.seed)
    
    # Construct output paths
    model_folder = f"{training_folder}/modelXGBoost"
    save_path = f"{model_folder}/xgb_model.pkl"
    
    print("=" * 70)
    print("XGBOOST FLOOD MODEL - OPTIMIZED")
    print("=" * 70)
    
    if config.use_adaptive_sampling:
        print(f"Adaptive sampling: ENABLED (sqrt compression)")
        print(f"Target range: {config.min_target_rate:.1%} - {config.max_target_rate:.1%}")
    else:
        print(f"Fixed sampling: {(config.min_target_rate + config.max_target_rate)/2:.1%}")
    
    print(f"Min samples per class: {config.min_samples_per_class}")
    print(f"Hyperparameter tuning: {'ENABLED' if config.tune_hyperparameters else 'DISABLED'}")
    print(f"Optimization metric: MCC (Matthews Correlation Coefficient)")
    
    dataset = open_parquet_dataset(parquet_uri)
    watersheds = get_unique_watersheds(dataset)
    feature_cols = get_feature_columns(dataset)
    
    print(f"Watersheds: {len(watersheds)}")
    print(f"Features: {len(feature_cols)}")
    
    folds = spatial_kfold_split(watersheds, config)
    print(f"Test watersheds: {folds[0]['test']}")
    
    prevalence_map = None
    if config.use_adaptive_sampling:
        train_val_ws = [w for w in watersheds if w not in folds[0]['test']]
        prevalence_map = compute_watershed_prevalence(dataset, train_val_ws)
    
    cv_thresholds = []
    fold_rows = []
    best_params = None
    
    for fold_idx, fold in enumerate(folds):
        print(f"\n{'='*60}")
        print(f"Fold {fold_idx + 1}/{len(folds)}")
        print(f"{'='*60}")
        print(f"  Train: {len(fold['train'])} watersheds")
        print(f"  Val: {len(fold['val'])} watersheds")
        
        means, stds = compute_normalization_stats(dataset, feature_cols, fold['train'])
        
        X_train, y_train = load_data_streaming(
            dataset, fold['train'], feature_cols, means, stds, config, 
            mode='train', prevalence_map=prevalence_map
        )
        
        X_val, y_val = load_data_streaming(
            dataset, fold['val'], feature_cols, means, stds, config, 
            mode='val', prevalence_map=None
        )
        
        if len(X_train) == 0 or len(X_val) == 0:
            print("  Warning: empty dataset, skipping fold")
            cv_thresholds.append(0.5)
            continue
        
        if len(np.unique(y_train)) < 2:
            print(f"  Warning: single-class training data, skipping fold")
            cv_thresholds.append(0.5)
            continue
        
        if config.tune_hyperparameters and fold_idx == 0:
            best_params = tune_hyperparameters(X_train, y_train, config)
        
        if best_params:
            xgb_params = best_params.copy()
            xgb_params['objective'] = config.objective
            xgb_params['tree_method'] = config.tree_method
            xgb_params['eval_metric'] = config.eval_metric
            xgb_params['random_state'] = config.seed + fold_idx
            xgb_params['n_jobs'] = config.n_jobs
        else:
            xgb_params = {
                'n_estimators': config.n_estimators,
                'max_depth': config.max_depth,
                'learning_rate': config.learning_rate,
                'subsample': config.subsample,
                'colsample_bytree': config.colsample_bytree,
                'min_child_weight': config.min_child_weight,
                'gamma': config.gamma,
                'reg_alpha': config.reg_alpha,
                'reg_lambda': config.reg_lambda,
                'objective': config.objective,
                'tree_method': config.tree_method,
                'eval_metric': config.eval_metric,
                'random_state': config.seed + fold_idx,
                'n_jobs': config.n_jobs
            }
        
        # NO scale_pos_weight - using sampling only
        
        print(f"  Training XGBoost...")
        print(f"  Positive rates: train={y_train.mean():.4f}, val={y_val.mean():.4f}")
        
        xgb_model = xgb.XGBClassifier(**xgb_params)
        xgb_model.fit(
            X_train, y_train,
            eval_set=[(X_val, y_val)],
            verbose=False
        )
        
        y_train_proba = xgb_model.predict_proba(X_train)[:, 1]
        y_val_proba = xgb_model.predict_proba(X_val)[:, 1]
        
        best_threshold = find_optimal_threshold_mcc(y_val, y_val_proba)
        cv_thresholds.append(best_threshold)
        
        train_metrics = compute_metrics(y_train, y_train_proba, best_threshold)
        val_metrics = compute_metrics(y_val, y_val_proba, best_threshold)
        
        print(f"  Val MCC={val_metrics['mcc']:.4f}, F1={val_metrics['f1']:.4f}, "
              f"PR-AUC={val_metrics['pr_auc']:.4f}, Threshold={best_threshold:.3f}")
        
        train_metrics['fold'] = fold_idx + 1
        train_metrics['set'] = 'train'
        val_metrics['fold'] = fold_idx + 1
        val_metrics['set'] = 'val'
        fold_rows.extend([train_metrics, val_metrics])
        
        del xgb_model, X_train, y_train, X_val, y_val
        gc.collect()
    
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
    
    final_threshold = np.mean(cv_thresholds) if cv_thresholds else 0.5
    print(f"Average threshold from CV: {final_threshold:.3f}")
    
    test_ws = folds[0]['test']
    train_ws = [w for w in watersheds if w not in test_ws]
    
    means, stds = compute_normalization_stats(dataset, feature_cols, train_ws)
    
    X_train_final, y_train_final = load_data_streaming(
        dataset, train_ws, feature_cols, means, stds, config, 
        mode='train', prevalence_map=prevalence_map
    )
    
    print(f"Final training on {len(X_train_final):,} pixels")
    
    if len(np.unique(y_train_final)) < 2:
        print("WARNING: Single-class final training data")
    
    if best_params:
        final_xgb_params = best_params.copy()
        final_xgb_params['objective'] = config.objective
        final_xgb_params['tree_method'] = config.tree_method
        final_xgb_params['eval_metric'] = config.eval_metric
        final_xgb_params['random_state'] = config.seed
        final_xgb_params['n_jobs'] = config.n_jobs
    else:
        final_xgb_params = {
            'n_estimators': config.n_estimators,
            'max_depth': config.max_depth,
            'learning_rate': config.learning_rate,
            'subsample': config.subsample,
            'colsample_bytree': config.colsample_bytree,
            'min_child_weight': config.min_child_weight,
            'gamma': config.gamma,
            'reg_alpha': config.reg_alpha,
            'reg_lambda': config.reg_lambda,
            'objective': config.objective,
            'tree_method': config.tree_method,
            'eval_metric': config.eval_metric,
            'random_state': config.seed,
            'n_jobs': config.n_jobs
        }
    
    final_xgb = xgb.XGBClassifier(**final_xgb_params)
    final_xgb.fit(X_train_final, y_train_final)
    
    if test_ws:
        X_test, y_test = load_data_streaming(
            dataset, test_ws, feature_cols, means, stds, config, mode='test'
        )
        if len(X_test) > 0:
            y_test_proba = final_xgb.predict_proba(X_test)[:, 1]
            test_metrics = compute_metrics(y_test, y_test_proba, final_threshold)
            summary_metrics['test'] = test_metrics
            print(f"Test MCC={test_metrics['mcc']:.4f}, F1={test_metrics['f1']:.4f}, "
                  f"PR-AUC={test_metrics['pr_auc']:.4f}")
    
    feature_importance_df = pd.DataFrame({
        'feature': feature_cols,
        'importance': final_xgb.feature_importances_
    }).sort_values('importance', ascending=False)
    
    results = {
        'model': final_xgb,
        'threshold': final_threshold,
        'feature_cols': feature_cols,
        'means': means,
        'stds': stds,
        'config': config,
        'test_watersheds': test_ws,
        'metrics_summary': summary_metrics,
        'feature_importance': feature_importance_df,
        'best_hyperparameters': best_params,
        'cv_results': folds_df
    }
    
    # Save model and results
    with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as tmp:
        temp_path = tmp.name
    
    with open(temp_path, 'wb') as f:
        pickle.dump(results, f)
    
    if save_path.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = save_path.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_path, bucket, key)
        print(f"\nModel saved to {save_path}")
    else:
        import shutil
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        shutil.move(temp_path, save_path)
        print(f"\nModel saved to {save_path}")
    
    os.remove(temp_path) if os.path.exists(temp_path) else None
    
    # Create diagnostic plots in same folder
    create_diagnostic_plots(results, model_folder)
    
    return results



# ======================== Diagnostic Plots ========================


def create_diagnostic_plots(model_results: Dict, output_folder: str):
    """Create comprehensive diagnostic plots."""
    
    print("\nCreating diagnostic plots...")
    
    feature_importance = model_results['feature_importance']
    cv_results = model_results.get('cv_results')
    
    def save_and_upload_plot(fig, s3_path):
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            temp_path = tmp.name
        
        fig.savefig(temp_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        
        if s3_path.startswith('s3://'):
            s3_client = boto3.client('s3')
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            s3_client.upload_file(temp_path, bucket, key)
            print(f"  Created: {s3_path}")
        else:
            import shutil
            os.makedirs(os.path.dirname(s3_path), exist_ok=True)
            shutil.move(temp_path, s3_path)
            print(f"  Created: {s3_path}")
        
        os.remove(temp_path) if os.path.exists(temp_path) else None
    
    # 1. Feature Importance (Top 20)
    fig, ax = plt.subplots(figsize=(12, 10))
    top_features = feature_importance.head(20)
    
    ax.barh(range(len(top_features)), top_features['importance'], color='steelblue')
    ax.set_yticks(range(len(top_features)))
    ax.set_yticklabels(top_features['feature'])
    ax.set_xlabel('Feature Importance (Gain)', fontsize=12)
    ax.set_title('XGBoost Feature Importance (Top 20)', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    importance_path = f"{output_folder}/xgboost_feature_importance.png"
    save_and_upload_plot(fig, importance_path)
    
    # 2. Importance Distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].hist(feature_importance['importance'], bins=50, color='steelblue', alpha=0.7)
    axes[0].set_xlabel('Importance Value', fontsize=12)
    axes[0].set_ylabel('Count', fontsize=12)
    axes[0].set_title('Distribution of Feature Importances', fontsize=13, fontweight='bold')
    axes[0].grid(alpha=0.3)
    
    sorted_imp = feature_importance['importance'].sort_values()
    axes[1].plot(range(len(sorted_imp)), sorted_imp, color='steelblue', linewidth=2)
    axes[1].set_xlabel('Feature Rank', fontsize=12)
    axes[1].set_ylabel('Importance Value', fontsize=12)
    axes[1].set_title('Sorted Feature Importances', fontsize=13, fontweight='bold')
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    dist_path = f"{output_folder}/xgboost_importance_distribution.png"
    save_and_upload_plot(fig, dist_path)
    
    # 3. CV Performance
    if cv_results is not None and not cv_results.empty:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        metrics_to_plot = ['mcc', 'f1', 'precision', 'recall']
        titles = ['MCC (Primary Metric)', 'F1 Score', 'Precision', 'Recall']
        
        for idx, (metric, title) in enumerate(zip(metrics_to_plot, titles)):
            ax = axes[idx // 2, idx % 2]
            
            train_vals = cv_results[cv_results['set'] == 'train'][metric].values
            val_vals = cv_results[cv_results['set'] == 'val'][metric].values
            
            x = np.arange(len(train_vals))
            width = 0.35
            
            ax.bar(x - width/2, train_vals, width, label='Train', color='steelblue', alpha=0.8)
            ax.bar(x + width/2, val_vals, width, label='Val', color='coral', alpha=0.8)
            
            ax.set_xlabel('Fold', fontsize=11)
            ax.set_ylabel(metric.upper(), fontsize=11)
            ax.set_title(title, fontsize=12, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f'{i+1}' for i in range(len(train_vals))])
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        cv_path = f"{output_folder}/xgboost_cv_performance.png"
        save_and_upload_plot(fig, cv_path)
    
    print("✅ Diagnostic plots complete!")


# ======================== Visualization Class ========================

class XGBoostFloodVisualizer:
    """Generate visualizations matching GLM/RF format."""
    
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
        
        m = re.search(r"Lat_([-\d\.]+)_Lon_([-\d\.]+)", training_folder)
        self.test_lat = m.group(1)
        self.test_lon = m.group(2)
        self.test_watershed = f"Lat_{self.test_lat}_Lon_{self.test_lon}"
        
        self.s3_client = boto3.client('s3')
        configure_gdal()

    def _validate_and_align_features(self, df: pd.DataFrame) -> np.ndarray:
        """Extract features in correct order, adding missing columns as NaN.
        
        Ensures prediction data matches training feature order exactly.
        """
        # Check for missing features
        missing_features = [f for f in self.feature_cols if f not in df.columns]
        
        if missing_features:
            print(f"  Warning: {len(missing_features)} features missing from parquet, filling with NaN")
            print(f"    Missing: {missing_features[:5]}{'...' if len(missing_features) > 5 else ''}")
            for feat in missing_features:
                df[feat] = np.nan
        
        # Extract in exact training order
        X = df[self.feature_cols].to_numpy(dtype=np.float32, copy=False)
        
        return X

    def _create_target_watershed_diagnostic_png(self, predictions, output_folder):
        """Create diagnostic PNG showing ONLY target watershed with boundary overlay."""
        
        # Get flood raster path for target watershed
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        
        # Open to get dimensions
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
        if ref_ds is None:
            print("  ⚠️  Cannot open reference raster for diagnostic")
            return
        
        H, W = ref_ds.RasterYSize, ref_ds.RasterXSize
        geotransform = ref_ds.GetGeoTransform()
        
        # Load boundary for target watershed only
        gpkg_path = (
            f"{self.training_folder}/watershedBoundaries/"
            f"processedWatershedAndBasins_Lon{self.test_lon}_Lat{self.test_lat}.gpkg"
        )
        
        # Generate prediction for target watershed only
        dataset = open_parquet_dataset(self.parquet_uri)
        filt = ds.field("Watershed_Loc") == self.test_watershed
        scanner = make_scanner(dataset, filter=filt)
        
        pixel_data = []
        for batch in scanner_batches(scanner):
            if batch.num_rows > 0:
                df = pa.Table.from_batches([batch]).to_pandas()
                pixel_data.append(df)
        
        if not pixel_data:
            print("  ⚠️  No data for target watershed diagnostic")
            return
        
        df_all = pd.concat(pixel_data, ignore_index=True)
        df_all['Date_Normalized'] = df_all['Date'].apply(self._normalize_date_format)
        
        # Use first available date
        dates = sorted(df_all['Date_Normalized'].dropna().unique())
        if not dates:
            return
        
        date_str = dates[0]
        date_data = df_all[df_all['Date_Normalized'] == date_str]
        
        # Generate prediction for target watershed
        X = self._validate_and_align_features(date_data)
        rows = date_data['PixelRow'].to_numpy(dtype=np.int32)
        cols = date_data['PixelCol'].to_numpy(dtype=np.int32)
        
        valid_mask = ~np.all(np.isnan(X), axis=1)
        pred_img = np.full((H, W), 255, dtype=np.uint8)
        
        if valid_mask.sum() > 0:
            X_valid = X[valid_mask]
            rows_valid = rows[valid_mask]
            cols_valid = cols[valid_mask]
            
            X_valid = np.where(np.isnan(X_valid), self.means, X_valid)
            X_valid = (X_valid - self.means) / self.stds
            
            y_proba_valid = self.model.predict_proba(X_valid)[:, 1]
            y_pred_valid = (y_proba_valid >= self.threshold).astype(np.uint8)
            pred_img[rows_valid, cols_valid] = y_pred_valid
        
        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))
        
        pred_disp = pred_img.astype(float)
        pred_disp[pred_img == 255] = np.nan
        
        cmap = matplotlib.colors.ListedColormap(['tan', 'blue'])
        ax.imshow(pred_disp, cmap=cmap, vmin=0, vmax=1, 
                 extent=[0, W, H, 0], origin='upper')
        
        # Load and overlay boundary
        try:
            import geopandas as gpd
            from matplotlib.patches import Polygon
            
            vsi_gpkg = gpkg_path.replace('s3://', '/vsis3/')
            gdf = gpd.read_file(vsi_gpkg)
            
            def geo_to_pixel(x_geo, y_geo, gt):
                px = (x_geo - gt[0]) / gt[1]
                py = (y_geo - gt[3]) / gt[5]
                return px, py
            
            for geom in gdf.geometry:
                if geom is None or geom.is_empty:
                    continue
                
                polys = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
                
                for poly in polys:
                    coords = list(poly.exterior.coords)
                    pixel_coords = [geo_to_pixel(x, y, geotransform) for x, y in coords]
                    
                    polygon = Polygon(pixel_coords, fill=False, edgecolor='red',
                                    linewidth=2.5, linestyle='-', zorder=10)
                    ax.add_patch(polygon)
            
            print(f"  Added {len(gdf)} boundary polygons")
        except Exception as e:
            print(f"  Could not load boundaries: {e}")
        
        ax.set_xlim(0, W)
        ax.set_ylim(H, 0)
        ax.set_title(f'Target Watershed Only: {self.test_watershed}\nDate: {date_str}', 
                    fontsize=12, fontweight='bold')
        ax.axis('off')
        plt.tight_layout()
        
        import tempfile
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            temp_png = tmp.name
        plt.savefig(temp_png, dpi=150, bbox_inches='tight')
        plt.close()
        
        s3_path = f"{output_folder}/diagnostic_target_watershed.png"
        bucket, key = s3_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(temp_png, bucket, key)
        os.remove(temp_png)
        
        ref_ds = None
        print(f"  Created target watershed diagnostic: {s3_path}")

    
    def generate_all_outputs(self, output_folder: Optional[str] = None,
                            select_dates: Optional[List[str]] = None,
                            generate_forecasts_for_all_subwatersheds: bool = False):
        """Generate all visualization outputs.
        
        Args:
            output_folder: Output folder path (defaults to training_folder/assessXGBoostModel)
            select_dates: List of date strings (YYYYMMDD format) to generate predictions for.
                         If provided, generates prediction maps (not confusion matrices).
            generate_forecasts_for_all_subwatersheds: If True, generates predictions for all
                                                      watersheds in parquet file. If False,
                                                      only generates for test watershed.
        """
        
        # Import parquet generation function
        from processDataToParquet import build_flood_training_parquet_v3
        
        if output_folder is None:
            output_folder = f"{self.training_folder}/assessXGBoostModel"
        
        print("\n" + "="*60)
        print("GENERATING XGBOOST MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        
        if generate_forecasts_for_all_subwatersheds:
            print(f"Mode: All watersheds")
        else:
            print(f"Test watershed: {self.test_watershed}")
        
        if select_dates:
            print(f"Selected dates: {select_dates}")
        
        print(f"Output folder: {output_folder}")
        
        # Parse training_folder to extract parameters for parquet generation
        # Format: s3://bucket/prefix/trainingFolderFor_Lat_{lat}_Lon_{lon}
        training_folder_match = re.match(
            r's3://([^/]+)/(.+)/trainingFolderFor_Lat_([-\d\.]+)_Lon_([-\d\.]+)',
            self.training_folder
        )
        
        if not training_folder_match:
            raise ValueError(f"Cannot parse training_folder format: {self.training_folder}")
        
        bucket = training_folder_match.group(1)
        base_prefix = training_folder_match.group(2)
        base_lat_4 = training_folder_match.group(3)
        base_lon_4 = training_folder_match.group(4)
        
        # Generate parquet file
        is_prediction_mode = (select_dates is not None)
        
        if is_prediction_mode:
            print("\nGenerating prediction parquet for selected dates...")
            try:
                parquet_uri = build_flood_training_parquet_v3(
                    base_lat_4=base_lat_4,
                    base_lon_4=base_lon_4,
                    bucket=bucket,
                    base_prefix=base_prefix,
                    select_dates=select_dates
                )
                print("✅ Prediction parquet generation complete")
            except Exception as e:
                raise RuntimeError(f"Failed to generate prediction parquet: {str(e)}")
        else:
            print("\nGenerating training parquet from DSWx data...")
            try:
                parquet_uri = build_flood_training_parquet_v3(
                    base_lat_4=base_lat_4,
                    base_lon_4=base_lon_4,
                    bucket=bucket,
                    base_prefix=base_prefix
                )
                print("✅ Training parquet generation complete")
            except Exception as e:
                raise RuntimeError(f"Failed to generate training parquet: {str(e)}")
        
        # Update parquet_uri to use the newly generated file
        self.parquet_uri = parquet_uri
        
        # Only create diagnostic plots if using default behavior
        if select_dates is None and not generate_forecasts_for_all_subwatersheds:
            create_diagnostic_plots(self.model_results, output_folder)
        
        # Determine which watersheds to process
        if generate_forecasts_for_all_subwatersheds:
            dataset = open_parquet_dataset(self.parquet_uri)
            watersheds_to_process = get_unique_watersheds(dataset)
            watershed_suffix = "all_watersheds"
        else:
            watersheds_to_process = [self.test_watershed]
            watershed_suffix = self.test_watershed
        
        # Generate predictions or confusion matrices based on mode
        if generate_forecasts_for_all_subwatersheds:
            predictions, confusions, per_date_metrics = self._generate_predictions_multi_watershed(
                watersheds_to_process, select_dates, is_prediction_mode
            )
        else:
            predictions, confusions, per_date_metrics = self._generate_predictions(
                select_dates, is_prediction_mode
            )
        
        # Check if predictions were generated
        if not predictions:
            print("\n⚠️  No predictions generated - no data found for specified parameters")
            return {
                'prediction_tiff': None,
                'confusion_matrix_tiff': None if not is_prediction_mode else None,
                'png_visualizations': [],
                'csv_outputs': {}
            }
        
        # Get reference dataset for geotransform
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
        
        dates = sorted(predictions.keys())
        
        # Create appropriate output files
        if is_prediction_mode:
            # Prediction mode - only save predictions
            pred_tiff = f"{output_folder}/xgboost_predictions_{watershed_suffix}.tif"
            self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "XGBoost_Prediction")
            
            png_paths = self._create_prediction_pngs(predictions, dates, per_date_metrics, output_folder)
            csv_paths = self._write_csvs(output_folder, per_date_metrics)

            self._create_target_watershed_diagnostic_png(predictions, output_folder)
            
            # Existing PNG generation
            png_paths = self._create_prediction_pngs(predictions, dates, per_date_metrics, output_folder)
            csv_paths = self._write_csvs(output_folder, per_date_metrics)
            print("\n✅ XGBoost prediction visualization complete!")
            
            return {
                'prediction_tiff': pred_tiff,
                'png_visualizations': png_paths,
                'csv_outputs': csv_paths
            }
        else:
            # Confusion matrix mode (default behavior)
            pred_tiff = f"{output_folder}/xgboost_predictions_{watershed_suffix}.tif"
            cm_tiff = f"{output_folder}/xgboost_confusion_matrix_{watershed_suffix}.tif"
            
            self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "XGBoost_Prediction")
            self._write_multiband_tiff(ref_ds, confusions, dates, cm_tiff, "XGBoost_CM")
            
            png_paths = self._create_pngs(confusions, dates, per_date_metrics, output_folder)
            csv_paths = self._write_csvs(output_folder, per_date_metrics)
            
            print("\n✅ XGBoost visualization complete!")
            
            return {
                'prediction_tiff': pred_tiff,
                'confusion_matrix_tiff': cm_tiff,
                'png_visualizations': png_paths,
                'csv_outputs': csv_paths
            }
    
    def _generate_predictions(self, select_dates: Optional[List[str]] = None, 
                             is_prediction_mode: bool = False):
        """Generate predictions for test watershed with proper spatial alignment.
        
        Args:
            select_dates: Optional list of dates (YYYYMMDD) to filter to
            is_prediction_mode: If True, only generate predictions (no confusion matrices)
        """
        
        # Load data from parquet
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
        
        # Normalize date column to YYYYMMDD format
        df_all['Date_Normalized'] = df_all['Date'].apply(self._normalize_date_format)
        
        # Filter to selected dates if provided
        if select_dates is not None:
            select_dates_normalized = [self._normalize_date_format(d) for d in select_dates]
            available_dates = sorted(df_all['Date_Normalized'].dropna().unique())
            print(f"  Available dates: {len(available_dates)} total")
            if available_dates:
                print(f"  Date range: {available_dates[0]} to {available_dates[-1]}")
            print(f"  Requested dates: {select_dates_normalized}")
            
            df_all = df_all[df_all['Date_Normalized'].isin(select_dates_normalized)]
            
            if df_all.empty:
                print(f"  ⚠️  No data found for selected dates in watershed {self.test_watershed}")
                return {}, {}, {}
            
            print(f"  Found {len(df_all)} pixels across {df_all['Date_Normalized'].nunique()} matching dates")
        
        dates = sorted(df_all['Date_Normalized'].dropna().unique())
        
        # CRITICAL FIX: Always use reference raster dimensions, not data-inferred dimensions
        # This ensures consistent array sizes regardless of which pixels have data
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
        
        if ref_ds is not None:
            H, W = ref_ds.RasterYSize, ref_ds.RasterXSize
            print(f"  Using reference raster dimensions: {H}x{W}")
            ref_ds = None
        else:
            # Fallback to data-inferred dimensions (should not happen)
            H = int(df_all['PixelRow'].max()) + 1
            W = int(df_all['PixelCol'].max()) + 1
            print(f"  Warning: Using data-inferred dimensions: {H}x{W}")
        
        predictions = {}
        confusions = {} if not is_prediction_mode else None
        per_date_metrics = {}
        
        print(f"  Processing {len(dates)} dates...")
        
        for date_idx, date_str in enumerate(dates):
            if (date_idx + 1) % 10 == 0:
                print(f"    Date {date_idx + 1}/{len(dates)}")
            
            date_data = df_all[df_all['Date_Normalized'] == date_str]
            
            # CRITICAL FIX: Validate feature columns exist and are in correct order
            missing_features = [f for f in self.feature_cols if f not in date_data.columns]
            if missing_features:
                print(f"    Warning: {len(missing_features)} features missing for date {date_str}, filling with NaN")
                for feat in missing_features:
                    date_data[feat] = np.nan
            
            # Extract features in EXACT training order
            X = date_data[self.feature_cols].to_numpy(dtype=np.float32, copy=False)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8, copy=False) if not is_prediction_mode else None
            rows = date_data['PixelRow'].to_numpy(dtype=np.int32, copy=False)
            cols = date_data['PixelCol'].to_numpy(dtype=np.int32, copy=False)
            
            # Identify valid pixels: at least one non-NaN feature AND within bounds
            valid_mask = ~np.all(np.isnan(X), axis=1)
            in_bounds = (rows >= 0) & (rows < H) & (cols >= 0) & (cols < W)
            valid_mask = valid_mask & in_bounds
            
            # Initialize prediction array with NoData (full reference dimensions)
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            
            if valid_mask.sum() > 0:
                # Extract valid pixels only
                X_valid = X[valid_mask]
                rows_valid = rows[valid_mask]
                cols_valid = cols[valid_mask]
                
                # Fill remaining NaNs with training means and normalize (SAME AS TRAINING)
                X_valid = np.where(np.isnan(X_valid), self.means, X_valid)
                X_valid = (X_valid - self.means) / self.stds
                
                # Make predictions
                y_proba_valid = self.model.predict_proba(X_valid)[:, 1]
                y_pred_valid = (y_proba_valid >= self.threshold).astype(np.uint8)
                
                # Place predictions in output array
                pred_img[rows_valid, cols_valid] = y_pred_valid
            
            predictions[date_str] = pred_img
            
            # Create confusion matrices if not in prediction mode
            if not is_prediction_mode and y_true is not None:
                cm_img = np.full((H, W), 255, dtype=np.uint8)
                
                if valid_mask.sum() > 0:
                    y_true_valid = y_true[valid_mask]
                    y_pred_valid = pred_img[rows_valid, cols_valid]
                    
                    # Confusion matrix values: 0=TN, 1=TP, 2=FN, 3=FP
                    cm_vals = np.where(
                        y_true_valid == 0,
                        np.where(y_pred_valid == 0, 0, 3),
                        np.where(y_pred_valid == 1, 1, 2)
                    )
                    cm_img[rows_valid, cols_valid] = cm_vals
                
                confusions[date_str] = cm_img
                
                # Calculate metrics on valid pixels only
                if valid_mask.sum() > 0:
                    X_valid = X[valid_mask]
                    y_true_valid = y_true[valid_mask]
                    X_valid = np.where(np.isnan(X_valid), self.means, X_valid)
                    X_valid = (X_valid - self.means) / self.stds
                    y_proba_valid = self.model.predict_proba(X_valid)[:, 1]
                    y_pred_valid = (y_proba_valid >= self.threshold).astype(np.uint8)
                    
                    try:
                        auc = roc_auc_score(y_true_valid, y_proba_valid)
                    except:
                        auc = np.nan
                    
                    try:
                        pr_auc = average_precision_score(y_true_valid, y_proba_valid)
                    except:
                        pr_auc = np.nan
                    
                    per_date_metrics[date_str] = {
                        'date': date_str,
                        'MCC': matthews_corrcoef(y_true_valid, y_pred_valid) if len(np.unique(y_true_valid)) > 1 and len(np.unique(y_pred_valid)) > 1 else 0.0,
                        'F1': f1_score(y_true_valid, y_pred_valid, zero_division=0),
                        'Accuracy': accuracy_score(y_true_valid, y_pred_valid),
                        'Precision': precision_score(y_true_valid, y_pred_valid, zero_division=0),
                        'Recall': recall_score(y_true_valid, y_pred_valid, zero_division=0),
                        'ROC_AUC': auc,
                        'PR_AUC': pr_auc,
                        'Kappa': cohen_kappa_score(y_true_valid, y_pred_valid),
                        'actual_pos_frac': y_true_valid.mean(),
                        'pred_pos_frac': y_pred_valid.mean(),
                        'valid_pixel_count': int(valid_mask.sum())
                    }
            else:
                # Prediction mode - only store prediction fraction
                if valid_mask.sum() > 0:
                    pred_fraction = pred_img[pred_img != 255].mean() if (pred_img != 255).any() else 0.0
                    per_date_metrics[date_str] = {
                        'date': date_str,
                        'pred_pos_frac': pred_fraction,
                        'valid_pixel_count': int(valid_mask.sum())
                    }
        
        return predictions, confusions if confusions is not None else {}, per_date_metrics


    def _generate_predictions_multi_watershed(self, watersheds: List[str],
                                             select_dates: Optional[List[str]] = None,
                                             is_prediction_mode: bool = False):
        """Generate predictions for multiple watersheds.
        
        Args:
            watersheds: List of watershed IDs to process
            select_dates: Optional list of dates to filter to
            is_prediction_mode: If True, only generate predictions (no confusion matrices)
        """
        
        dataset = open_parquet_dataset(self.parquet_uri)
        
        # Load all data for specified watersheds
        filt = ds.field("Watershed_Loc").isin(pa.array(watersheds))
        scanner = make_scanner(dataset, filter=filt)
        
        pixel_data = []
        for batch in scanner_batches(scanner):
            if batch.num_rows > 0:
                df = pa.Table.from_batches([batch]).to_pandas()
                pixel_data.append(df)
        
        if not pixel_data:
            return {}, {}, {}
        
        df_all = pd.concat(pixel_data, ignore_index=True)
        
        # Normalize date column to YYYYMMDD format
        df_all['Date_Normalized'] = df_all['Date'].apply(self._normalize_date_format)
        
        # Filter to selected dates if provided
        if select_dates is not None:
            select_dates_normalized = [self._normalize_date_format(d) for d in select_dates]
            available_dates = sorted(df_all['Date_Normalized'].unique())
            print(f"  Available dates across all watersheds: {len(available_dates)} total")
            print(f"  Date range: {available_dates[0]} to {available_dates[-1]}")
            print(f"  Requested dates: {select_dates_normalized}")
            
            found_dates = set(select_dates_normalized) & set(available_dates)
            missing_dates = set(select_dates_normalized) - set(available_dates)
            
            if missing_dates:
                print(f"  ⚠️  Missing dates: {sorted(missing_dates)}")
            if found_dates:
                print(f"  ✓ Found dates: {sorted(found_dates)}")
            
            df_all = df_all[df_all['Date_Normalized'].isin(select_dates_normalized)]
            
            if df_all.empty:
                print(f"  ⚠️  No data found for any selected dates across {len(watersheds)} watersheds")
                return {}, {}, {}
            
            print(f"  Found {len(df_all)} pixels across {df_all['Date_Normalized'].nunique()} matching dates")
        
        dates = sorted(df_all['Date_Normalized'].unique())
        H = df_all['PixelRow'].max() + 1
        W = df_all['PixelCol'].max() + 1
        
        predictions = {}
        confusions = {} if not is_prediction_mode else None
        per_date_metrics = {}
        
        print(f"  Processing {len(dates)} dates across {len(watersheds)} watersheds...")
        
        for date_idx, date_str in enumerate(dates):
            if (date_idx + 1) % 10 == 0:
                print(f"    Date {date_idx + 1}/{len(dates)}")
            
            date_data = df_all[df_all['Date_Normalized'] == date_str]
            
            X = self._validate_and_align_features(date_data)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8) if not is_prediction_mode else None
            rows = date_data['PixelRow'].to_numpy()
            cols = date_data['PixelCol'].to_numpy()
            
            # Identify pixels with valid data (at least some non-NaN features)
            valid_mask = ~np.all(np.isnan(X), axis=1)
            
            # Initialize prediction array with NoData
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            
            if valid_mask.sum() > 0:
                # Only process pixels with valid data
                X_valid = X[valid_mask]
                rows_valid = rows[valid_mask]
                cols_valid = cols[valid_mask]
                
                # Fill remaining NaNs with means and normalize
                X_valid = np.where(np.isnan(X_valid), self.means, X_valid)
                X_valid = (X_valid - self.means) / self.stds
                
                # Make predictions only on valid pixels
                y_proba_valid = self.model.predict_proba(X_valid)[:, 1]
                y_pred_valid = (y_proba_valid >= self.threshold).astype(np.uint8)
                
                # Place predictions in output array
                pred_img[rows_valid, cols_valid] = y_pred_valid
            
            predictions[date_str] = pred_img
            
            # Only create confusion matrices if not in prediction mode
            if not is_prediction_mode and y_true is not None:
                cm_img = np.full((H, W), 255, dtype=np.uint8)
                
                if valid_mask.sum() > 0:
                    y_true_valid = y_true[valid_mask]
                    y_pred_valid = pred_img[rows_valid, cols_valid]
                    
                    cm_vals = np.where(
                        y_true_valid == 0,
                        np.where(y_pred_valid == 0, 0, 3),
                        np.where(y_pred_valid == 1, 1, 2)
                    )
                    cm_img[rows_valid, cols_valid] = cm_vals
                
                confusions[date_str] = cm_img
                
                # Calculate metrics (aggregated across watersheds, only on valid pixels)
                if valid_mask.sum() > 0:
                    pred_fraction = pred_img[pred_img != 255].mean() if (pred_img != 255).any() else 0.0
                    per_date_metrics[date_str] = {
                        'date': date_str,
                        'pred_pos_frac': pred_fraction,
                        'valid_pixel_count': int(valid_mask.sum())
                    }
                    
                    if len(np.unique(y_true[valid_mask])) > 1:
                        y_true_valid = y_true[valid_mask]
                        per_date_metrics[date_str].update({
                            'actual_pos_frac': y_true_valid.mean(),
                            'Accuracy': accuracy_score(y_true_valid, pred_img[rows_valid, cols_valid])
                        })
            else:
                # Prediction mode - only store prediction fraction for valid pixels
                if valid_mask.sum() > 0:
                    pred_fraction = pred_img[pred_img != 255].mean() if (pred_img != 255).any() else 0.0
                    per_date_metrics[date_str] = {
                        'date': date_str,
                        'pred_pos_frac': pred_fraction,
                        'valid_pixel_count': int(valid_mask.sum())
                    }
        
        return predictions, confusions if confusions is not None else {}, per_date_metrics

    def _create_prediction_pngs(self, predictions, dates, metrics, output_folder):
        """Create PNG visualizations with ALL watershed boundaries overlaid."""
        png_paths = []
        
        # Get all watershed locations from the data
        dataset = open_parquet_dataset(self.parquet_uri)
        all_watersheds = get_unique_watersheds(dataset)
        
        # Load boundaries for ALL watersheds
        all_pixel_polygons = []
        for ws_loc in all_watersheds:
            # Parse lat/lon from watershed location
            match = re.match(r"Lat_([-\d\.]+)_Lon_([-\d\.]+)", ws_loc)
            if not match:
                continue
            
            lat4, lon4 = match.groups()
            
            # Get reference raster for this watershed
            flood_tiff_path = (
                f"{self.training_folder}/floodMaps/"
                f"dswx_s1_timeseries_subwatershed_Lon{lon4}_Lat{lat4}.tif"
            )
            ws_ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
            
            if ws_ref_ds is None:
                continue
            
            geotransform = ws_ref_ds.GetGeoTransform()
            ws_ref_ds = None
            
            # Load boundary for this watershed
            gpkg_path = (
                f"{self.training_folder}/watershedBoundaries/"
                f"processedWatershedAndBasins_Lon{lon4}_Lat{lat4}.gpkg"
            )
            
            try:
                import geopandas as gpd
                vsi_gpkg = gpkg_path.replace('s3://', '/vsis3/')
                gdf = gpd.read_file(vsi_gpkg)
                
                def geo_to_pixel(x_geo, y_geo, gt):
                    px = (x_geo - gt[0]) / gt[1]
                    py = (y_geo - gt[3]) / gt[5]
                    return px, py
                
                for geom in gdf.geometry:
                    if geom is None or geom.is_empty:
                        continue
                    
                    polys = [geom] if geom.geom_type == 'Polygon' else list(geom.geoms)
                    
                    for poly in polys:
                        coords = list(poly.exterior.coords)
                        pixel_coords = [geo_to_pixel(x, y, geotransform) for x, y in coords]
                        all_pixel_polygons.append(pixel_coords)
            except Exception as e:
                print(f"  Could not load boundaries for {ws_loc}: {e}")
        
        print(f"  Loaded {len(all_pixel_polygons)} boundary polygons across all watersheds")
        
        # Create PNGs with all boundaries
        for date_str in dates:
            pred = predictions.get(date_str)
            if pred is None or np.count_nonzero(pred != 255) < 100:
                continue
            
            m = metrics.get(date_str, {})
            
            fig, ax = plt.subplots(figsize=(12, 10))
            
            pred_disp = pred.astype(float)
            pred_disp[pred == 255] = np.nan
            
            cmap = matplotlib.colors.ListedColormap(['tan', 'blue'])
            ax.imshow(pred_disp, cmap=cmap, vmin=0, vmax=1)
            
            # Overlay ALL watershed boundaries
            if all_pixel_polygons:
                from matplotlib.patches import Polygon
                for poly_coords in all_pixel_polygons:
                    polygon = Polygon(poly_coords, fill=False, edgecolor='red',
                                    linewidth=2.0, linestyle='-', zorder=10)
                    ax.add_patch(polygon)
            
            legend_elements = [
                mpatches.Patch(color='tan', label='Predicted Land'),
                mpatches.Patch(color='blue', label='Predicted Water')
            ]
            
            if all_pixel_polygons:
                from matplotlib.lines import Line2D
                legend_elements.append(
                    Line2D([0], [0], color='red', linewidth=2.0, 
                          label='Watershed Boundaries')
                )
            
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
            ax.set_title(f'XGBoost Flood Prediction - {date_str}\nAll Watersheds', 
                        fontsize=14, fontweight='bold')
            
            text = (
                f"Date: {date_str}\n"
                f"Model: XGBoost\n"
                f"Watersheds: {len(all_watersheds)}\n"
                f"Predicted Flood Fraction: {m.get('pred_pos_frac', 0):.3f}"
            )
            props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
            ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11, va='top', bbox=props)
            
            ax.axis('off')
            plt.tight_layout()
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()
            
            s3_path = f"{output_folder}/xgboost_prediction_{date_str}.png"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(temp_png, bucket, key)
            os.remove(temp_png)
            png_paths.append(s3_path)
        
        print(f"  Created {len(png_paths)} prediction PNG visualizations with all watershed boundaries")
        return png_paths
    
    
    def _write_multiband_tiff(self, ref_ds, arrays, dates, output_path, prefix):
        """Write multiband TIFF using prediction array dimensions."""
        
        # Get dimensions from first prediction array
        first_array = arrays[dates[0]]
        height, width = first_array.shape
        
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
            temp_path = tmp.name
        
        driver = gdal.GetDriverByName('GTiff')
        out = driver.Create(
            temp_path, width, height, len(dates),  # Use prediction dimensions
            gdal.GDT_Byte, options=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"]
        )
        
        # Copy geotransform and projection from reference
        out.SetGeoTransform(ref_ds.GetGeoTransform())
        out.SetProjection(ref_ds.GetProjection())
        
        for i, date_str in enumerate(dates, 1):
            band = out.GetRasterBand(i)
            band.SetNoDataValue(255)
            band.SetDescription(f"{prefix}_{date_str}")
            arr = arrays.get(date_str, np.full((height, width), 255, dtype=np.uint8))
            band.WriteArray(arr)
        
        out = None
        
        bucket, key = output_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(temp_path, bucket, key)
        os.remove(temp_path)
        print(f"  Created: {output_path}")

    
    def _create_pngs(self, confusions, dates, metrics, output_folder):
        """Create PNG visualizations."""
        png_paths = []
        
        dates_to_process = sorted(dates)
        if hasattr(self.config, 'max_png_dates') and self.config.max_png_dates is not None:
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
            ax.set_title(f'XGBoost Confusion Matrix - {date_str}', fontsize=14, fontweight='bold')
            
            auc_val = m.get('ROC_AUC', np.nan)
            pr_auc_val = m.get('PR_AUC', np.nan)
            auc_str = f"{auc_val:.3f}" if np.isfinite(auc_val) else "N/A"
            pr_auc_str = f"{pr_auc_val:.3f}" if np.isfinite(pr_auc_val) else "N/A"
            
            text = (
                f"Date: {date_str}\n"
                f"Model: XGBoost\n"
                f"MCC: {m.get('MCC', 0):.3f}\n"
                f"F1 Score: {m.get('F1', 0):.3f}\n"
                f"Accuracy: {m.get('Accuracy', 0):.3f}\n"
                f"Precision: {m.get('Precision', 0):.3f}\n"
                f"Recall: {m.get('Recall', 0):.3f}\n"
                f"ROC AUC: {auc_str}\n"
                f"PR AUC: {pr_auc_str}\n"
                f"Kappa: {m.get('Kappa', 0):.3f}"
            )
            props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
            ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11, va='top', bbox=props)
            
            ax.axis('off')
            plt.tight_layout()
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()
            
            s3_path = f"{output_folder}/xgboost_confusion_matrix_{date_str}.png"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(temp_png, bucket, key)
            os.remove(temp_png)
            png_paths.append(s3_path)
        
        print(f"  Created {len(png_paths)} PNG visualizations")
        return png_paths
    
    def _write_csvs(self, output_folder, per_date_metrics):
        """Write CSV outputs."""
        paths = {}
        
        if per_date_metrics:
            df = pd.DataFrame.from_dict(per_date_metrics, orient='index')
            df = df.sort_values("date") if "date" in df.columns else df
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
            df.to_csv(tmp_csv, index=False)
            s3_path = f"{output_folder}/xgboost_per_date_metrics_{self.test_watershed}.csv"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(tmp_csv, bucket, key)
            os.remove(tmp_csv)
            paths['per_date_metrics'] = s3_path
        
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
                    "pr_auc": m.get("pr_auc", np.nan)
                })
            
            if rows:
                df_sum = pd.DataFrame(rows)
                with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                    tmp_csv = tmp.name
                df_sum.to_csv(tmp_csv, index=False)
                s3_path = f"{output_folder}/xgboost_metrics_summary.csv"
                bucket, key = s3_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(tmp_csv, bucket, key)
                os.remove(tmp_csv)
                paths['metrics_summary'] = s3_path
        
        feature_importance = self.model_results.get('feature_importance')
        if feature_importance is not None:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
            feature_importance.to_csv(tmp_csv, index=False)
            s3_path = f"{output_folder}/xgboost_feature_importances.csv"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(tmp_csv, bucket, key)
            os.remove(tmp_csv)
            paths['feature_importances'] = s3_path
        
        return paths
    def _normalize_date_format(self, date_value) -> str:
        """Normalize various date formats to YYYYMMDD string."""
        if pd.isna(date_value):
            return None
        
        # If it's already a string, try to normalize it
        if isinstance(date_value, str):
            # Remove dashes/slashes
            clean = date_value.replace('-', '').replace('/', '')
            if len(clean) == 8 and clean.isdigit():
                return clean
            return date_value
        
        # If it's an integer
        if isinstance(date_value, (int, np.integer)):
            date_str = str(date_value)
            if len(date_str) == 8:
                return date_str
        
        # If it's a datetime
        if isinstance(date_value, (pd.Timestamp, np.datetime64)):
            return pd.Timestamp(date_value).strftime('%Y%m%d')
        
        # Try converting to string
        return str(date_value)


def load_model_from_s3(model_path: str) -> Dict:
    """Load trained XGBoost model from S3 pickle file.
    
    Args:
        model_path: S3 path to xgb_model.pkl (e.g., s3://bucket/path/modelXGBoost/xgb_model.pkl)
    
    Returns:
        Dictionary containing model, threshold, normalization params, etc.
    """
    import tempfile
    import pickle
    import boto3
    
    if model_path.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = model_path.replace('s3://', '').split('/', 1)
        
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as tmp:
            temp_path = tmp.name
        
        s3_client.download_file(bucket, key, temp_path)
        
        with open(temp_path, 'rb') as f:
            results = pickle.load(f)
        
        os.remove(temp_path)
        print(f"✅ Model loaded from {model_path}")
    else:
        with open(model_path, 'rb') as f:
            results = pickle.load(f)
        print(f"✅ Model loaded from {model_path}")
    
    return results

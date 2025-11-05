"""
randomForestFloods.py - FULLY OPTIMIZED VERSION

Random Forest flood prediction with all GLM improvements applied.
Uses proper sampling, hyperparameter tuning, and comprehensive diagnostics.

Key features:
- Adaptive target-based sampling for cross-watershed balance
- Optional hyperparameter tuning (RandomizedSearchCV)
- Sampling-only class balancing (no class_weight)
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

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from scipy.stats import randint, uniform

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
class RandomForestConfig:
    """Configuration for Random Forest flood model."""
    
    # Random Forest parameters (can be tuned)
    n_estimators: int = 460 # 300
    max_depth: Optional[int] = 30 #20
    min_samples_split: int = 3 #5
    min_samples_leaf: int = 1 # 2
    max_features: str = 0.7 #"sqrt"
    bootstrap: bool = True
    
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
    """
    Compute adaptive target rate using square root compression.
    
    Parameters:
    -----------
    natural_prevalence : float
        Natural positive rate in basin (0-1)
    min_rate, max_rate : float
        Bounds for target rate
    
    Returns:
    --------
    float : Target positive rate for sampling
    """
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
    """
    Sample to achieve TARGET positive rate per watershed.
    
    Strategy:
    1. Keep all positives (minority class)
    2. Sample negatives to achieve target positive rate
    3. Enforce minimum samples per class
    4. Return None if insufficient data
    """
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
                       stds: np.ndarray, config: RandomForestConfig,
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
    
    # Normalize (Random Forests don't require this, but we do it for consistency)
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
        'max_depth': [10, 15, 20, 25, 30, None],
        'min_samples_split': randint(2, 11),
        'min_samples_leaf': randint(1, 5),
        'max_features': ['sqrt', 'log2', 0.5, 0.7]
    }


def tune_hyperparameters(X_train: np.ndarray, y_train: np.ndarray, 
                        config: RandomForestConfig) -> Dict[str, Any]:
    """Perform randomized hyperparameter search."""
    print("  Performing hyperparameter tuning...")
    print(f"    Trials: {config.n_tuning_iterations}")
    print(f"    Inner CV folds: {config.tuning_cv_folds}")
    
    base_model = RandomForestClassifier(
        random_state=config.seed,
        n_jobs=1,
        bootstrap=config.bootstrap
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

def spatial_kfold_split(watersheds: List[str], config: RandomForestConfig):
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


def train_random_forest_flood_model(parquet_uri: str, config: Optional[RandomForestConfig] = None,
                                   save_path: Optional[str] = None) -> Dict:
    """Train Random Forest model with k-fold CV and adaptive sampling."""
    
    config = config or RandomForestConfig()
    np.random.seed(config.seed)
    
    print("=" * 70)
    print("RANDOM FOREST FLOOD MODEL - OPTIMIZED")
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
            rf_params = best_params.copy()
            rf_params['random_state'] = config.seed + fold_idx
            rf_params['n_jobs'] = config.n_jobs
            rf_params['bootstrap'] = config.bootstrap
        else:
            rf_params = {
                'n_estimators': config.n_estimators,
                'max_depth': config.max_depth,
                'min_samples_split': config.min_samples_split,
                'min_samples_leaf': config.min_samples_leaf,
                'max_features': config.max_features,
                'bootstrap': config.bootstrap,
                'random_state': config.seed + fold_idx,
                'n_jobs': config.n_jobs
            }
        
        # NO class_weight - using sampling only
        
        print(f"  Training Random Forest...")
        print(f"  Positive rates: train={y_train.mean():.4f}, val={y_val.mean():.4f}")
        
        rf = RandomForestClassifier(**rf_params)
        rf.fit(X_train, y_train)
        
        y_train_proba = rf.predict_proba(X_train)[:, 1]
        y_val_proba = rf.predict_proba(X_val)[:, 1]
        
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
        
        del rf, X_train, y_train, X_val, y_val
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
        final_rf_params = best_params.copy()
        final_rf_params['random_state'] = config.seed
        final_rf_params['n_jobs'] = config.n_jobs
        final_rf_params['bootstrap'] = config.bootstrap
    else:
        final_rf_params = {
            'n_estimators': config.n_estimators,
            'max_depth': config.max_depth,
            'min_samples_split': config.min_samples_split,
            'min_samples_leaf': config.min_samples_leaf,
            'max_features': config.max_features,
            'bootstrap': config.bootstrap,
            'random_state': config.seed,
            'n_jobs': config.n_jobs
        }
    
    final_rf = RandomForestClassifier(**final_rf_params)
    final_rf.fit(X_train_final, y_train_final)
    
    if test_ws:
        X_test, y_test = load_data_streaming(
            dataset, test_ws, feature_cols, means, stds, config, mode='test'
        )
        if len(X_test) > 0:
            y_test_proba = final_rf.predict_proba(X_test)[:, 1]
            test_metrics = compute_metrics(y_test, y_test_proba, final_threshold)
            summary_metrics['test'] = test_metrics
            print(f"Test MCC={test_metrics['mcc']:.4f}, F1={test_metrics['f1']:.4f}, "
                  f"PR-AUC={test_metrics['pr_auc']:.4f}")
    
    feature_importance_df = pd.DataFrame({
        'feature': feature_cols,
        'importance': final_rf.feature_importances_
    }).sort_values('importance', ascending=False)
    
    results = {
        'model': final_rf,
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
    
    if save_path:
        with open(save_path, 'wb') as f:
            pickle.dump(results, f)
        print(f"\nModel saved to {save_path}")
    
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
    ax.set_xlabel('Feature Importance', fontsize=12)
    ax.set_title('Random Forest Feature Importance (Top 20)', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    importance_path = f"{output_folder}/randomforest_feature_importance.png"
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
    dist_path = f"{output_folder}/randomforest_importance_distribution.png"
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
        cv_path = f"{output_folder}/randomforest_cv_performance.png"
        save_and_upload_plot(fig, cv_path)
    
    print("✅ Diagnostic plots complete!")


# ======================== Visualization Class ========================

class RandomForestFloodVisualizer:
    """Generate visualizations matching GLM format."""
    
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
        if m:
            self.test_lat = m.group(1)
            self.test_lon = m.group(2)
            self.test_watershed = f"Lat_{self.test_lat}_Lon_{self.test_lon}"
        
        self.s3_client = boto3.client('s3')
        configure_gdal()
    
    def generate_all_outputs(self, output_folder: Optional[str] = None):
        """Generate all visualization outputs."""
        
        if output_folder is None:
            output_folder = f"{self.training_folder}/assessRandomForestModel"
        
        print("\n" + "="*60)
        print("GENERATING RANDOM FOREST MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        print(f"Test watershed: {self.test_watershed}")
        print(f"Output folder: {output_folder}")
        
        create_diagnostic_plots(self.model_results, output_folder)
        
        predictions, confusions, per_date_metrics = self._generate_predictions()
        
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
        
        dates = sorted(predictions.keys())
        
        pred_tiff = f"{output_folder}/randomforest_predictions_{self.test_watershed}.tif"
        cm_tiff = f"{output_folder}/randomforest_confusion_matrix_{self.test_watershed}.tif"
        
        self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "RandomForest_Prediction")
        self._write_multiband_tiff(ref_ds, confusions, dates, cm_tiff, "RandomForest_CM")
        
        png_paths = self._create_pngs(confusions, dates, per_date_metrics, output_folder)
        csv_paths = self._write_csvs(output_folder, per_date_metrics)
        
        print("\n✅ Random Forest visualization complete!")
        
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
        
        print(f"  Processing {len(dates)} dates...")
        
        for date_idx, date_str in enumerate(dates):
            if (date_idx + 1) % 10 == 0:
                print(f"    Date {date_idx + 1}/{len(dates)}")
            
            date_data = df_all[df_all['Date'] == date_str]
            
            X = date_data[self.feature_cols].to_numpy(dtype=np.float32)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8)
            rows = date_data['PixelRow'].to_numpy()
            cols = date_data['PixelCol'].to_numpy()
            
            X = np.where(np.isnan(X), self.means, X)
            X = (X - self.means) / self.stds
            
            y_proba = self.model.predict_proba(X)[:, 1]
            y_pred = (y_proba >= self.threshold).astype(np.uint8)
            
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
                auc = roc_auc_score(y_true, y_proba)
            except:
                auc = np.nan
            
            try:
                pr_auc = average_precision_score(y_true, y_proba)
            except:
                pr_auc = np.nan
            
            per_date_metrics[date_str] = {
                'date': date_str,
                'MCC': matthews_corrcoef(y_true, y_pred) if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1 else 0.0,
                'F1': f1_score(y_true, y_pred, zero_division=0),
                'Accuracy': accuracy_score(y_true, y_pred),
                'Precision': precision_score(y_true, y_pred, zero_division=0),
                'Recall': recall_score(y_true, y_pred, zero_division=0),
                'ROC_AUC': auc,
                'PR_AUC': pr_auc,
                'Kappa': cohen_kappa_score(y_true, y_pred),
                'actual_pos_frac': y_true.mean(),
                'pred_pos_frac': y_pred.mean()
            }
        
        return predictions, confusions, per_date_metrics
    
    def _write_multiband_tiff(self, ref_ds, arrays, dates, output_path, prefix):
        """Write multiband TIFF."""
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
            band.WriteArray(arrays.get(date_str, np.full((ref_ds.RasterYSize, ref_ds.RasterXSize), 255, dtype=np.uint8)))
        
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
            ax.set_title(f'Random Forest Confusion Matrix - {date_str}', fontsize=14, fontweight='bold')
            
            auc_val = m.get('ROC_AUC', np.nan)
            pr_auc_val = m.get('PR_AUC', np.nan)
            auc_str = f"{auc_val:.3f}" if np.isfinite(auc_val) else "N/A"
            pr_auc_str = f"{pr_auc_val:.3f}" if np.isfinite(pr_auc_val) else "N/A"
            
            text = (
                f"Date: {date_str}\n"
                f"Model: Random Forest\n"
                f"MCC: {m.get('MCC', 0):.3f}\n"
                f"F1 Score: {m.get('F1', 0):.3f}\n"
                f"Accuracy: {m.get('Accuracy', 0):.3f}\n"
                f"Precision: {m.get('Precision', 0):.3f}\n"
                f"Recall: {m.get('Recall', 0):.3f}\n"
                f"ROC AUC: {auc_str}\n"
                f"PR AUC: {pr_auc_str}\n"
                f"Kappa: {m.get('Kappa', 0):.3f}"
            )
            props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.8)
            ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11, va='top', bbox=props)
            
            ax.axis('off')
            plt.tight_layout()
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()
            
            s3_path = f"{output_folder}/randomforest_confusion_matrix_{date_str}.png"
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
            s3_path = f"{output_folder}/randomforest_per_date_metrics_{self.test_watershed}.csv"
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
                s3_path = f"{output_folder}/randomforest_metrics_summary.csv"
                bucket, key = s3_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(tmp_csv, bucket, key)
                os.remove(tmp_csv)
                paths['metrics_summary'] = s3_path
        
        feature_importance = self.model_results.get('feature_importance')
        if feature_importance is not None:
            with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
            feature_importance.to_csv(tmp_csv, index=False)
            s3_path = f"{output_folder}/randomforest_feature_importances.csv"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(tmp_csv, bucket, key)
            os.remove(tmp_csv)
            paths['feature_importances'] = s3_path
        
        return paths

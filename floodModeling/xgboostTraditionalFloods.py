"""
xgboostTraditionalFloods.py - TRADITIONAL FEATURES ONLY

XGBoost flood prediction using ONLY features available to traditional models
like Google Flood Hub: HAND, elevation, and river_basin_streamflow_pred.

This serves as a baseline to test H1: whether hydrologically-informed features
improve performance compared to traditional discharge-centric approaches.

Key differences from xgboostFloods.py:
- Uses only 3 features: HAND, elevation, river_basin_streamflow_pred
- Saves outputs to assessXGBoostTraditionalModel/
- Labels outputs as "XGBoost-Traditional" for clarity
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
#from statsmodels.stats.outliers_influence import variance_inflation_factor
import statsmodels.api as sm

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

from xgboostFloods import XGBoostFloodVisualizer as BaseVisualizer

class XGBoostFloodVisualizer(BaseVisualizer):
    def generate_all_outputs(self, output_folder: Optional[str] = None, **kwargs):
        if output_folder is None:
            output_folder = f"{self.training_folder}/assessXGBoostTraditionalModel"
        
        print(f"Generating outputs to: {output_folder}")
        return super().generate_all_outputs(output_folder=output_folder, **kwargs)


# ======================== Configuration ========================

@dataclass
class XGBoostConfig:
    """Configuration for XGBoost flood model."""
    
    # XGBoost parameters (can be tuned)
    n_estimators: int = 320
    max_depth: int = 8
    learning_rate: float = 0.1
    subsample: float = 0.8
    colsample_bytree: float = 0.8
    min_child_weight: int = 3
    gamma: float = 0.1
    reg_alpha: float = 0.05
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


def get_traditional_feature_columns(dataset) -> List[str]:
    """Get ONLY traditional features used by models like Google Flood Hub.
    
    Traditional models rely on:
    - HAND (Height Above Nearest Drainage)
    - elevation (base DEM elevation)
    - river_basin_streamflow_pred (basin-scale discharge)
    
    This deliberately excludes:
    - Multi-scale terrain features
    - Soil properties
    - Land cover
    - Local streamflow
    - Weather variables
    
    Returns:
        List containing exactly 3 features
    """
    # Define the exact traditional features we want
    target_features = ['HAND', 'elevation', 'river_basin_streamflow_pred']
    
    # Get all available columns
    all_columns = [field.name for field in dataset.schema 
                   if field.name not in {"Target", "Date", "Watershed_Loc", "PixelRow", "PixelCol"}]
    
    # Find matching features (case-sensitive exact match)
    traditional_features = []
    for target in target_features:
        if target in all_columns:
            traditional_features.append(target)
        else:
            # Try case-insensitive fallback
            for col in all_columns:
                if col.lower() == target.lower():
                    traditional_features.append(col)
                    print(f"  Warning: Using '{col}' for '{target}' (case mismatch)")
                    break
    
    if len(traditional_features) != 3:
        # Provide helpful error message
        available_candidates = [col for col in all_columns if 
                               'hand' in col.lower() or 
                               'elevation' in col.lower() or 
                               'streamflow' in col.lower()]
        raise ValueError(
            f"Expected 3 traditional features, found {len(traditional_features)}.\n"
            f"Required: {target_features}\n"
            f"Found: {traditional_features}\n"
            f"Available candidates in dataset: {available_candidates}"
        )
    
    print(f"  ✅ Traditional features selected: {traditional_features}")
    print(f"     This represents the baseline approach used by Google Flood Hub")
    
    return traditional_features


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
    
    # Normalize
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
        'subsample': uniform(0.6, 0.4),
        'colsample_bytree': uniform(0.6, 0.4),
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
    """Train XGBoost-Traditional model with k-fold CV and adaptive sampling.
    
    Args:
        parquet_uri: Path to parquet file
        training_folder: S3 folder path for saving outputs
        config: Optional configuration object
    """
    
    config = config or XGBoostConfig()
    np.random.seed(config.seed)
    
    # Construct output paths for TRADITIONAL model
    model_folder = f"{training_folder}/modelXGBoostTraditional"
    save_path = f"{model_folder}/xgb_traditional_model.pkl"
    
    print("=" * 70)
    print("XGBOOST-TRADITIONAL FLOOD MODEL")
    print("Using ONLY traditional features (HAND, elevation, basin streamflow)")
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
    
    # ✅ KEY CHANGE: Use traditional features only
    feature_cols = get_traditional_feature_columns(dataset)
    
    print(f"Watersheds: {len(watersheds)}")
    print(f"Features: {len(feature_cols)} (TRADITIONAL BASELINE)")
    
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
            xgb_params['early_stopping_rounds'] = config.early_stopping_rounds
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
                'early_stopping_rounds': config.early_stopping_rounds,
                'random_state': config.seed + fold_idx,
                'n_jobs': config.n_jobs
            }
        
        print(f"  Training XGBoost-Traditional...")
        print(f"  Positive rates: train={y_train.mean():.4f}, val={y_val.mean():.4f}")
        
        xgb_model = xgb.XGBClassifier(**xgb_params)
        xgb_model.fit(
            X_train, y_train,
            eval_set=[(X_train, y_train), (X_val, y_val)],
            verbose=False
        )
        
        evals_results = xgb_model.evals_result()
        
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
        
        plot_training_curves(evals_results, fold_idx + 1, model_folder)
        plot_roc_curves(y_train, y_train_proba, y_val, y_val_proba, fold_idx + 1, model_folder)
        
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
        final_xgb_params['early_stopping_rounds'] = config.early_stopping_rounds
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
            'early_stopping_rounds': config.early_stopping_rounds,
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
    
    # Save model
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
    
    # Feature Importance (should only be 3 features)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.barh(range(len(feature_importance)), feature_importance['importance'], color='steelblue')
    ax.set_yticks(range(len(feature_importance)))
    ax.set_yticklabels(feature_importance['feature'])
    ax.set_xlabel('Feature Importance (Gain)', fontsize=12)
    ax.set_title('XGBoost-Traditional Feature Importance\n(HAND, Elevation, Basin Streamflow Only)', 
                fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    importance_path = f"{output_folder}/xgboost_traditional_feature_importance.png"
    save_and_upload_plot(fig, importance_path)
    
    # CV Performance
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
            ax.set_title(f'{title}\n(Traditional Features Only)', fontsize=12, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f'{i+1}' for i in range(len(train_vals))])
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        cv_path = f"{output_folder}/xgboost_traditional_cv_performance.png"
        save_and_upload_plot(fig, cv_path)
    
    print("✅ Diagnostic plots complete!")


def plot_training_curves(evals_results: dict, fold_idx: int, output_folder: str):
    """Plot training and validation loss curves."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    train_loss = evals_results['validation_0']['logloss']
    val_loss = evals_results['validation_1']['logloss']
    epochs = range(len(train_loss))
    
    ax.plot(epochs, train_loss, label='Training Loss', linewidth=2, color='steelblue')
    ax.plot(epochs, val_loss, label='Validation Loss', linewidth=2, color='coral')
    
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Log Loss', fontsize=12)
    ax.set_title(f'Training History (Traditional Features) - Fold {fold_idx}', 
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_path = tmp.name
    plt.savefig(temp_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    s3_path = f"{output_folder}/training_curve_fold_{fold_idx}.png"
    if s3_path.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = s3_path.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_path, bucket, key)
    else:
        import shutil
        os.makedirs(os.path.dirname(s3_path), exist_ok=True)
        shutil.move(temp_path, s3_path)
    
    os.remove(temp_path) if os.path.exists(temp_path) else None


def plot_roc_curves(y_train, y_train_proba, y_val, y_val_proba, fold_idx: int, output_folder: str):
    """Plot ROC curves for training and validation sets."""
    from sklearn.metrics import roc_curve, auc
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    if len(np.unique(y_train)) > 1:
        fpr_train, tpr_train, _ = roc_curve(y_train, y_train_proba)
        roc_auc_train = auc(fpr_train, tpr_train)
        ax.plot(fpr_train, tpr_train, label=f'Train (AUC = {roc_auc_train:.3f})', 
                linewidth=2, color='steelblue')
    
    if len(np.unique(y_val)) > 1:
        fpr_val, tpr_val, _ = roc_curve(y_val, y_val_proba)
        roc_auc_val = auc(fpr_val, tpr_val)
        ax.plot(fpr_val, tpr_val, label=f'Validation (AUC = {roc_auc_val:.3f})', 
                linewidth=2, color='coral')
    
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, label='Random (AUC = 0.500)')
    
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title(f'ROC Curves (Traditional Features) - Fold {fold_idx}', 
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='lower right')
    ax.grid(alpha=0.3)
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    
    plt.tight_layout()
    
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        temp_path = tmp.name
    plt.savefig(temp_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    
    s3_path = f"{output_folder}/roc_curve_fold_{fold_idx}.png"
    if s3_path.startswith('s3://'):
        s3_client = boto3.client('s3')
        bucket, key = s3_path.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_path, bucket, key)
    else:
        import shutil
        os.makedirs(os.path.dirname(s3_path), exist_ok=True)
        shutil.move(temp_path, s3_path)
    
    os.remove(temp_path) if os.path.exists(temp_path) else None


# ======================== Main Entry Point ========================

def load_model_from_s3(model_path: str) -> Dict:
    """Load trained XGBoost-Traditional model from S3 pickle file."""
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
        print(f"✅ XGBoost-Traditional model loaded from {model_path}")
    else:
        with open(model_path, 'rb') as f:
            results = pickle.load(f)
        print(f"✅ XGBoost-Traditional model loaded from {model_path}")
    
    return results

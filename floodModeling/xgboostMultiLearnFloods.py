"""
xgboostMultiLearnFloods.py - MULTITASK LEARNING VERSION

XGBoost flood prediction with multitask learning for multiple related objectives.
Uses proper sampling, hyperparameter tuning, and comprehensive diagnostics.

Key features:
- Multitask learning with 5 objectives:
  * Primary: Flood/not flood (binary classification) - 50% weight
  * Auxiliary 1: river_basin_streamflow_pred (regression) - 12.5% weight
  * Auxiliary 2: elevation (regression) - 12.5% weight
  * Auxiliary 3: SILT_5_TO_15CM (regression) - 12.5% weight
  * Auxiliary 4: IO_LULC (classification) - 12.5% weight
- Adaptive target-based sampling for cross-watershed balance
- Optional hyperparameter tuning (RandomizedSearchCV)
- Sampling-only class balancing (no scale_pos_weight)
- Streaming data loading for memory efficiency
- MCC-based threshold optimization (flood task only)
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
    make_scorer, mean_squared_error, mean_absolute_error, r2_score
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
class XGBoostMultiLearnConfig:
    """Configuration for XGBoost MultiLearn flood model."""
    
    # XGBoostMultiLearn parameters (can be tuned)
    n_estimators: int = 320
    max_depth: int = 8
    learning_rate: float = 0.1
    subsample: float = 0.8
    colsample_bytree: float = 0.8
    min_child_weight: int = 3
    gamma: float = 0.1
    reg_alpha: float = 0.05
    reg_lambda: float = 1.0
    
    # Multitask learning configuration
    task_weights: Dict[str, float] = None  # Will be set in __post_init__
    auxiliary_targets: List[str] = None  # Will be set in __post_init__
    
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
    
    # XGBoostMultiLearn specific
    tree_method: str = "hist"
    objective: str = "binary:logistic"  # For primary flood task
    eval_metric: str = "logloss"
    early_stopping_rounds: int = 20
    
    # Visualization
    max_png_dates: Optional[int] = None
    
    def __post_init__(self):
        """Initialize task weights and auxiliary targets."""
        if self.task_weights is None:
            self.task_weights = {
                'flood': 0.50,  # Primary task
                'river_basin_streamflow_pred': 0.25,
                'elevation': 0.125,
#                'SILT_5_TO_15CM': 0.125,
                'IO_LULC': 0.125
            }
        
        if self.auxiliary_targets is None:
            self.auxiliary_targets = [
                'river_basin_streamflow_pred',
                'elevation',
#                'SILT_5_TO_15CM',
                'IO_LULC'
            ]


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


def get_feature_columns(dataset, exclude_auxiliary_targets: List[str]) -> List[str]:
    """Get all feature columns excluding auxiliary targets and standard exclusions."""
    exclude = {"Target", "Date", "Watershed_Loc", "PixelRow", "PixelCol"}
    exclude.update(exclude_auxiliary_targets)
    
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
                              min_samples_per_class: int = 20) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """Sample to achieve TARGET positive rate per watershed.
    
    Returns: X_sampled, y_sampled, sampled_indices
    """
    n = len(y)
    if n == 0:
        return None, None, None
    
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    
    n_pos = len(pos_idx)
    n_neg = len(neg_idx)
    
    if n_pos < min_samples_per_class or n_neg < min_samples_per_class:
        return None, None, None
    
    n_neg_target = int(n_pos * (1 - target_pos_rate) / target_pos_rate)
    n_neg_sample = min(n_neg_target, n_neg)
    
    total = n_pos + n_neg_sample
    
    if total > max_samples:
        scale = max_samples / total
        n_pos_sample = max(min_samples_per_class, int(n_pos * scale))
        n_neg_sample = max_samples - n_pos_sample
        
        if n_pos_sample < min_samples_per_class or n_neg_sample < min_samples_per_class:
            return None, None, None
        
        sel_pos = np.random.choice(pos_idx, n_pos_sample, replace=False)
        sel_neg = np.random.choice(neg_idx, n_neg_sample, replace=False)
    else:
        sel_pos = pos_idx
        sel_neg = np.random.choice(neg_idx, n_neg_sample, replace=False)
    
    sel_idx = np.concatenate([sel_pos, sel_neg])
    np.random.shuffle(sel_idx)
    
    # ✅ Return indices directly - NO SLOW LOOP
    return X[sel_idx], y[sel_idx], sel_idx


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


def load_data_streaming_multitask(parquet_dataset, watersheds: List[str],
                                  feature_cols: List[str], auxiliary_targets: List[str],
                                  means: np.ndarray, stds: np.ndarray, 
                                  config: XGBoostMultiLearnConfig,
                                  mode: str = "train",
                                  prevalence_map: Optional[Dict[str, float]] = None) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """OPTIMIZED: Single-pass multitask data loading with adaptive streaming sampling."""
    
    print(f"  Loading {mode} data from {len(watersheds)} watersheds...")
    
    filt = ds.field("Watershed_Loc").isin(pa.array(watersheds))
    cols = feature_cols + ["Target"] + auxiliary_targets + ["Watershed_Loc"]
    scanner = make_scanner(parquet_dataset, columns=cols, filter=filt, 
                          batch_size=config.batch_size_loading)
    
    ws_data = defaultdict(lambda: {'X': [], 'y_flood': [], **{f'y_{t}': [] for t in auxiliary_targets}})
    
    for batch in scanner_batches(scanner):
        if batch.num_rows == 0:
            continue
        
        X_batch = arrow_batch_to_numpy(batch, feature_cols)
        
        # Get target indices
        target_idx = len(feature_cols)
        y_flood_batch = batch.column(target_idx).to_numpy(zero_copy_only=False).astype(np.uint8)
        
        # Get auxiliary targets
        aux_targets_batch = {}
        for i, aux_name in enumerate(auxiliary_targets):
            aux_idx = target_idx + 1 + i
            aux_targets_batch[aux_name] = batch.column(aux_idx).to_numpy(zero_copy_only=False).astype(np.float32)
        
        # Get watershed
        ws_idx = len(feature_cols) + 1 + len(auxiliary_targets)
        ws_batch = batch.column(ws_idx).to_pandas().values
        
        for ws in np.unique(ws_batch):
            mask = ws_batch == ws
            ws_data[ws]['X'].append(X_batch[mask])
            ws_data[ws]['y_flood'].append(y_flood_batch[mask])
            for aux_name in auxiliary_targets:
                ws_data[ws][f'y_{aux_name}'].append(aux_targets_batch[aux_name][mask])
    
    X_parts = []
    y_flood_parts = []
    y_aux_parts = {aux_name: [] for aux_name in auxiliary_targets}
    skipped = []
    
    for ws in watersheds:
        if ws not in ws_data:
            skipped.append(ws)
            continue
        
        X_ws = np.vstack(ws_data[ws]['X'])
        y_flood_ws = np.hstack(ws_data[ws]['y_flood'])
        y_aux_ws = {aux_name: np.hstack(ws_data[ws][f'y_{aux_name}']) for aux_name in auxiliary_targets}
        
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
            
            # ✅ Get indices directly - FAST!
            X_sampled, y_flood_sampled, sampled_idx = balanced_watershed_sample(
                X_ws, y_flood_ws,
                target_pos_rate=target_rate,
                max_samples=config.max_pixels_per_watershed,
                min_samples_per_class=config.min_samples_per_class
            )
            
            if X_sampled is None:
                skipped.append(ws)
                print(f"    Skipped {ws}: insufficient samples")
                continue
            
            # ✅ Apply same indices to auxiliary targets - INSTANT!
            y_aux_ws_sampled = {aux_name: y_aux_ws[aux_name][sampled_idx] for aux_name in auxiliary_targets}
            
            X_ws = X_sampled
            y_flood_ws = y_flood_sampled
            y_aux_ws = y_aux_ws_sampled
        
        X_parts.append(X_ws)
        y_flood_parts.append(y_flood_ws)
        for aux_name in auxiliary_targets:
            y_aux_parts[aux_name].append(y_aux_ws[aux_name])
    
    if skipped:
        print(f"    Skipped {len(skipped)} watersheds with insufficient data")
    
    if not X_parts:
        return np.empty((0, len(feature_cols)), dtype=np.float32), {'flood': np.empty((0,), dtype=np.uint8)}
    
    X = np.vstack(X_parts)
    y_flood = np.hstack(y_flood_parts)
    y_aux = {aux_name: np.hstack(y_aux_parts[aux_name]) for aux_name in auxiliary_targets}
    
    # Normalize (XGBoostMultiLearn benefits from normalization for regularization)
    X = np.where(np.isnan(X), means, X)
    X = (X - means) / stds
    
    # Combine all targets into dictionary
    y_all = {'flood': y_flood}
    y_all.update(y_aux)
    
    print(f"    Loaded {len(y_flood):,} pixels ({y_flood.mean():.4f} flood positive rate)")
    
    return X, y_all


# ======================== Metrics - MULTITASK ========================

def compute_metrics_classification(y_true: np.ndarray, y_proba: np.ndarray, threshold: float) -> Dict[str, float]:
    """Compute comprehensive metrics for classification tasks."""
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


def compute_metrics_regression(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Compute metrics for regression tasks."""
    # Remove NaNs
    mask = ~(np.isnan(y_true) | np.isnan(y_pred))
    y_true_clean = y_true[mask]
    y_pred_clean = y_pred[mask]
    
    if len(y_true_clean) == 0:
        return {
            "mse": np.nan,
            "rmse": np.nan,
            "mae": np.nan,
            "r2": np.nan
        }
    
    mse = mean_squared_error(y_true_clean, y_pred_clean)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true_clean, y_pred_clean)
    
    try:
        r2 = r2_score(y_true_clean, y_pred_clean)
    except:
        r2 = np.nan
    
    return {
        "mse": mse,
        "rmse": rmse,
        "mae": mae,
        "r2": r2
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
                        config: XGBoostMultiLearnConfig) -> Dict[str, Any]:
    """Perform randomized hyperparameter search on flood task."""
    print("  Performing hyperparameter tuning on primary flood task...")
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


# ======================== Training Pipeline - MULTITASK ========================

def spatial_kfold_split(watersheds: List[str], config: XGBoostMultiLearnConfig):
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


def train_xgboost_multilearn_flood_model(parquet_uri: str, training_folder: str, 
                                        config: Optional[XGBoostMultiLearnConfig] = None) -> Dict:
    """Train XGBoostMultiLearn model with CASCADING architecture.
    
    Architecture:
    1. Train auxiliary models on base features (53 features)
    2. Generate auxiliary predictions
    3. Augment features with auxiliary predictions (53 + 4 = 57 features)
    4. Train flood model on augmented features
    
    This implements TRUE multitask learning where auxiliary tasks directly
    contribute to flood prediction through feature augmentation.
    
    Args:
        parquet_uri: Path to parquet file
        training_folder: S3 folder path for saving outputs
        config: Optional configuration object
    """
    
    config = config or XGBoostMultiLearnConfig()
    np.random.seed(config.seed)
    
    # Construct output paths
    model_folder = f"{training_folder}/modelXGBoostMultiLearn"
    save_path = f"{model_folder}/xgboostmultilearn_model.pkl"
    
    print("=" * 70)
    print("XGBOOST MULTILEARN FLOOD MODEL - CASCADING ARCHITECTURE")
    print("=" * 70)
    
    print(f"\n🔗 Multitask learning: CASCADING (auxiliary → flood)")
    print(f"📊 Architecture:")
    print(f"   Step 1: Train auxiliary models on 53 base features")
    print(f"   Step 2: Generate auxiliary predictions")
    print(f"   Step 3: Train flood model on 53 + 4 = 57 features")
    print(f"\n🎯 Primary task: Flood prediction (weight: {config.task_weights['flood']:.1%})")
    print(f"🔧 Auxiliary tasks (used as feature generators):")
    for aux_name in config.auxiliary_targets:
        print(f"   - {aux_name} (weight: {config.task_weights[aux_name]:.1%})")
    
    if config.use_adaptive_sampling:
        print(f"\n📈 Adaptive sampling: ENABLED (sqrt compression)")
        print(f"   Target range: {config.min_target_rate:.1%} - {config.max_target_rate:.1%}")
    else:
        print(f"\n📈 Fixed sampling: {(config.min_target_rate + config.max_target_rate)/2:.1%}")
    
    print(f"📏 Min samples per class: {config.min_samples_per_class}")
    print(f"🔍 Hyperparameter tuning: {'ENABLED' if config.tune_hyperparameters else 'DISABLED'}")
    print(f"📐 Primary optimization metric: MCC (Matthews Correlation Coefficient)")
    
    dataset = open_parquet_dataset(parquet_uri)
    watersheds = get_unique_watersheds(dataset)
    feature_cols = get_feature_columns(dataset, config.auxiliary_targets)
    
    print(f"\n📂 Data: {len(watersheds)} watersheds")
    print(f"🔢 Base features: {len(feature_cols)} (excluding {len(config.auxiliary_targets)} converted to auxiliary targets)")
    print(f"🎯 Auxiliary targets: {config.auxiliary_targets}")
    
    folds = spatial_kfold_split(watersheds, config)
    print(f"🧪 Test watersheds: {folds[0]['test']}")
    
    prevalence_map = None
    if config.use_adaptive_sampling:
        train_val_ws = [w for w in watersheds if w not in folds[0]['test']]
        prevalence_map = compute_watershed_prevalence(dataset, train_val_ws)
    
    cv_thresholds = []
    fold_rows_flood = []
    fold_rows_aux = {aux_name: [] for aux_name in config.auxiliary_targets}
    best_params = None
    
    # Store models for each task
    task_models = {
        'flood': [],
        **{aux_name: [] for aux_name in config.auxiliary_targets}
    }
    
    for fold_idx, fold in enumerate(folds):
        print(f"\n{'='*70}")
        print(f"FOLD {fold_idx + 1}/{len(folds)}")
        print(f"{'='*70}")
        print(f"  📍 Train: {len(fold['train'])} watersheds")
        print(f"  📍 Val: {len(fold['val'])} watersheds")
        
        # Compute normalization stats
        means, stds = compute_normalization_stats(dataset, feature_cols, fold['train'])
        
        # Load data with base features only
        X_train, y_train_all = load_data_streaming_multitask(
            dataset, fold['train'], feature_cols, config.auxiliary_targets,
            means, stds, config, mode='train', prevalence_map=prevalence_map
        )
        
        X_val, y_val_all = load_data_streaming_multitask(
            dataset, fold['val'], feature_cols, config.auxiliary_targets,
            means, stds, config, mode='val', prevalence_map=None
        )
        
        if len(X_train) == 0 or len(X_val) == 0:
            print("  ⚠️  Warning: empty dataset, skipping fold")
            cv_thresholds.append(0.5)
            continue
        
        if len(np.unique(y_train_all['flood'])) < 2:
            print(f"  ⚠️  Warning: single-class training data, skipping fold")
            cv_thresholds.append(0.5)
            continue
        
        # Hyperparameter tuning on primary flood task (optional)
        if config.tune_hyperparameters and fold_idx == 0:
            best_params = tune_hyperparameters(X_train, y_train_all['flood'], config)
        
        # Set base parameters
        if best_params:
            base_xgb_params = best_params.copy()
            base_xgb_params['random_state'] = config.seed + fold_idx
            base_xgb_params['n_jobs'] = config.n_jobs
        else:
            base_xgb_params = {
                'n_estimators': config.n_estimators,
                'max_depth': config.max_depth,
                'learning_rate': config.learning_rate,
                'subsample': config.subsample,
                'colsample_bytree': config.colsample_bytree,
                'min_child_weight': config.min_child_weight,
                'gamma': config.gamma,
                'reg_alpha': config.reg_alpha,
                'reg_lambda': config.reg_lambda,
                'random_state': config.seed + fold_idx,
                'n_jobs': config.n_jobs
            }
        
        print(f"\n  {'─'*60}")
        print(f"  STEP 1: Training Auxiliary Models on {X_train.shape[1]} Base Features")
        print(f"  {'─'*60}")
        
        # Dictionary to store auxiliary models for this fold
        fold_aux_models = {}
        
        # Train auxiliary models (regression and classification)
        for aux_name in config.auxiliary_targets:
            is_classification = (aux_name == 'IO_LULC')
            
            if is_classification:
                print(f"  🔧 Training {aux_name} model (classification)...")
                aux_params = base_xgb_params.copy()
                
                # Clean training data
                train_mask = ~np.isnan(y_train_all[aux_name])
                val_mask = ~np.isnan(y_val_all[aux_name])
                
                if train_mask.sum() < 10 or val_mask.sum() < 10:
                    print(f"    ⚠️  Skipping {aux_name}: insufficient valid data")
                    fold_aux_models[aux_name] = None
                    task_models[aux_name].append(None)
                    continue
                
                # Get all unique classes from BOTH train and val
                train_classes = np.unique(y_train_all[aux_name][train_mask])
                val_classes = np.unique(y_val_all[aux_name][val_mask])
                all_classes = np.unique(np.concatenate([train_classes, val_classes]))
                n_classes = len(all_classes)
                
                # Create label mapping to consecutive integers
                class_to_idx = {cls: idx for idx, cls in enumerate(all_classes)}
                
                # Map classes to consecutive integers
                y_train_mapped = np.array([class_to_idx[int(c)] for c in y_train_all[aux_name][train_mask]])
                y_val_mapped = np.array([class_to_idx[int(c)] for c in y_val_all[aux_name][val_mask]])
                
                if n_classes <= 2:
                    aux_params['objective'] = 'binary:logistic'
                    aux_params['eval_metric'] = 'logloss'
                    aux_model = xgb.XGBClassifier(**aux_params)
                else:
                    aux_params['objective'] = 'multi:softprob'
                    aux_params['eval_metric'] = 'mlogloss'
                    aux_params['num_class'] = int(n_classes)
                    aux_model = xgb.XGBClassifier(**aux_params)
                
                # Train with mapped classes
                aux_model.fit(
                    X_train[train_mask], y_train_mapped,
                    eval_set=[(X_val[val_mask], y_val_mapped)],
                    verbose=False
                )
                
                # Store class mapping with model
                aux_model.class_mapping_ = class_to_idx
                fold_aux_models[aux_name] = aux_model
                task_models[aux_name].append(aux_model)
                
                # Evaluate
                train_acc = accuracy_score(y_train_mapped, aux_model.predict(X_train[train_mask]))
                val_acc = accuracy_score(y_val_mapped, aux_model.predict(X_val[val_mask]))
                
                print(f"    ✅ {aux_name} - Train Acc={train_acc:.4f}, Val Acc={val_acc:.4f} ({n_classes} classes)")
                
                fold_rows_aux[aux_name].append({
                    'fold': fold_idx + 1,
                    'set': 'train',
                    'accuracy': train_acc,
                    'n_classes': n_classes
                })
                fold_rows_aux[aux_name].append({
                    'fold': fold_idx + 1,
                    'set': 'val',
                    'accuracy': val_acc,
                    'n_classes': n_classes
                })
                
            else:  # Regression
                print(f"  🔧 Training {aux_name} model (regression)...")
                aux_params = base_xgb_params.copy()
                aux_params['objective'] = 'reg:squarederror'
                aux_params['tree_method'] = config.tree_method
                aux_params['eval_metric'] = 'rmse'
                
                aux_model = xgb.XGBRegressor(**aux_params)
                
                # Clean training data (remove NaNs)
                train_mask = ~np.isnan(y_train_all[aux_name])
                val_mask = ~np.isnan(y_val_all[aux_name])
                
                if train_mask.sum() < 10 or val_mask.sum() < 10:
                    print(f"    ⚠️  Skipping {aux_name}: insufficient valid data")
                    fold_aux_models[aux_name] = None
                    task_models[aux_name].append(None)
                    continue
                
                aux_model.fit(
                    X_train[train_mask], y_train_all[aux_name][train_mask],
                    eval_set=[(X_val[val_mask], y_val_all[aux_name][val_mask])],
                    verbose=False
                )
                fold_aux_models[aux_name] = aux_model
                task_models[aux_name].append(aux_model)
                
                # Evaluate
                y_train_aux_pred = aux_model.predict(X_train[train_mask])
                y_val_aux_pred = aux_model.predict(X_val[val_mask])
                
                train_aux_metrics = compute_metrics_regression(y_train_all[aux_name][train_mask], y_train_aux_pred)
                val_aux_metrics = compute_metrics_regression(y_val_all[aux_name][val_mask], y_val_aux_pred)
                
                print(f"    ✅ {aux_name} - Train RMSE={train_aux_metrics['rmse']:.4f}, Val RMSE={val_aux_metrics['rmse']:.4f}, Val R²={val_aux_metrics['r2']:.4f}")
                
                train_aux_metrics['fold'] = fold_idx + 1
                train_aux_metrics['set'] = 'train'
                val_aux_metrics['fold'] = fold_idx + 1
                val_aux_metrics['set'] = 'val'
                fold_rows_aux[aux_name].extend([train_aux_metrics, val_aux_metrics])
        
        print(f"\n  {'─'*60}")
        print(f"  STEP 2: Generating Auxiliary Predictions")
        print(f"  {'─'*60}")
        
        # Create augmented feature matrices
        X_train_aug = X_train.copy()
        X_val_aug = X_val.copy()
        auxiliary_feature_names = []
        
        # Generate predictions from auxiliary models and append as features
        for aux_name in config.auxiliary_targets:
            aux_model = fold_aux_models.get(aux_name)
            
            if aux_model is None:
                # If aux model failed, use zeros as placeholder
                train_aux_feature = np.zeros(len(X_train), dtype=np.float32)
                val_aux_feature = np.zeros(len(X_val), dtype=np.float32)
                print(f"    ⚠️  {aux_name}: model unavailable, using zeros")
            else:
                is_classification = (aux_name == 'IO_LULC')
                
                # Initialize prediction arrays
                train_aux_feature = np.zeros(len(X_train), dtype=np.float32)
                val_aux_feature = np.zeros(len(X_val), dtype=np.float32)
                
                # Generate predictions where data is valid
                train_mask = ~np.isnan(y_train_all[aux_name])
                val_mask = ~np.isnan(y_val_all[aux_name])
                
                if is_classification:
                    # For classification, use predicted class as feature
                    if train_mask.sum() > 0:
                        train_aux_feature[train_mask] = aux_model.predict(X_train[train_mask]).astype(np.float32)
                    if val_mask.sum() > 0:
                        val_aux_feature[val_mask] = aux_model.predict(X_val[val_mask]).astype(np.float32)
                else:
                    # For regression, use predicted value directly
                    if train_mask.sum() > 0:
                        train_aux_feature[train_mask] = aux_model.predict(X_train[train_mask])
                    if val_mask.sum() > 0:
                        val_aux_feature[val_mask] = aux_model.predict(X_val[val_mask])
                
                print(f"    ✅ {aux_name}: generated predictions (train mean={train_aux_feature.mean():.4f}, val mean={val_aux_feature.mean():.4f})")
            
            # Append as new feature column
            X_train_aug = np.column_stack([X_train_aug, train_aux_feature])
            X_val_aug = np.column_stack([X_val_aug, val_aux_feature])
            auxiliary_feature_names.append(f"{aux_name}_pred")
        
        n_base_features = X_train.shape[1]
        n_aug_features = X_train_aug.shape[1]
        print(f"\n  ✅ Feature augmentation complete: {n_base_features} → {n_aug_features} features")
        print(f"     Added features: {auxiliary_feature_names}")
        
        print(f"\n  {'─'*60}")
        print(f"  STEP 3: Training Flood Model on Augmented Features")
        print(f"  {'─'*60}")
        
        # Train flood model on augmented features
        print(f"  🎯 Training flood prediction model on {n_aug_features} features...")
        flood_params = base_xgb_params.copy()
        flood_params['objective'] = config.objective
        flood_params['tree_method'] = config.tree_method
        flood_params['eval_metric'] = config.eval_metric
        
        flood_model = xgb.XGBClassifier(**flood_params)
        flood_model.fit(
            X_train_aug,  # ✅ Use augmented features
            y_train_all['flood'],
            eval_set=[(X_val_aug, y_val_all['flood'])],  # ✅ Use augmented validation
            verbose=False
        )
        task_models['flood'].append(flood_model)
        
        # Evaluate flood task on augmented features
        y_train_flood_proba = flood_model.predict_proba(X_train_aug)[:, 1]
        y_val_flood_proba = flood_model.predict_proba(X_val_aug)[:, 1]
        
        flood_threshold = find_optimal_threshold_mcc(y_val_all['flood'], y_val_flood_proba)
        cv_thresholds.append(flood_threshold)
        
        train_flood_metrics = compute_metrics_classification(y_train_all['flood'], y_train_flood_proba, flood_threshold)
        val_flood_metrics = compute_metrics_classification(y_val_all['flood'], y_val_flood_proba, flood_threshold)
        
        print(f"  ✅ Flood - Val MCC={val_flood_metrics['mcc']:.4f}, F1={val_flood_metrics['f1']:.4f}, "
              f"PR-AUC={val_flood_metrics['pr_auc']:.4f}, Threshold={flood_threshold:.3f}")
        
        train_flood_metrics['fold'] = fold_idx + 1
        train_flood_metrics['set'] = 'train'
        val_flood_metrics['fold'] = fold_idx + 1
        val_flood_metrics['set'] = 'val'
        fold_rows_flood.extend([train_flood_metrics, val_flood_metrics])
        
        # Clean up memory
        del X_train, X_val, X_train_aug, X_val_aug, y_train_all, y_val_all
        gc.collect()
    
    # Create DataFrames from fold results
    folds_df_flood = pd.DataFrame(fold_rows_flood) if fold_rows_flood else pd.DataFrame()
    folds_df_aux = {aux_name: pd.DataFrame(fold_rows_aux[aux_name]) if fold_rows_aux[aux_name] else pd.DataFrame() 
                    for aux_name in config.auxiliary_targets}
    
    summary_metrics = {'flood': {}}
    
    if not folds_df_flood.empty:
        for s in ["train", "val"]:
            sub = folds_df_flood[folds_df_flood["set"] == s]
            if not sub.empty:
                summary_metrics['flood'][s] = sub.drop(columns=["fold", "set"]).mean(numeric_only=True).to_dict()
    
    # Summarize auxiliary metrics
    for aux_name in config.auxiliary_targets:
        summary_metrics[aux_name] = {}
        if not folds_df_aux[aux_name].empty:
            for s in ["train", "val"]:
                sub = folds_df_aux[aux_name][folds_df_aux[aux_name]["set"] == s]
                if not sub.empty:
                    summary_metrics[aux_name][s] = sub.drop(columns=["fold", "set"]).mean(numeric_only=True).to_dict()
    
    print("\n" + "=" * 70)
    print("FINAL MODEL TRAINING")
    print("=" * 70)
    
    final_threshold = np.mean(cv_thresholds) if cv_thresholds else 0.5
    print(f"✅ Average flood threshold from CV: {final_threshold:.3f}")
    
    test_ws = folds[0]['test']
    train_ws = [w for w in watersheds if w not in test_ws]
    
    # Load final training data
    means, stds = compute_normalization_stats(dataset, feature_cols, train_ws)
    
    X_train_final, y_train_final_all = load_data_streaming_multitask(
        dataset, train_ws, feature_cols, config.auxiliary_targets,
        means, stds, config, mode='train', prevalence_map=prevalence_map
    )
    
    print(f"\n📊 Final training on {len(X_train_final):,} pixels")
    
    # Set final parameters
    if best_params:
        final_base_params = best_params.copy()
        final_base_params['random_state'] = config.seed
        final_base_params['n_jobs'] = config.n_jobs
    else:
        final_base_params = {
            'n_estimators': config.n_estimators,
            'max_depth': config.max_depth,
            'learning_rate': config.learning_rate,
            'subsample': config.subsample,
            'colsample_bytree': config.colsample_bytree,
            'min_child_weight': config.min_child_weight,
            'gamma': config.gamma,
            'reg_alpha': config.reg_alpha,
            'reg_lambda': config.reg_lambda,
            'random_state': config.seed,
            'n_jobs': config.n_jobs
        }
    
    print(f"\n{'─'*60}")
    print(f"STEP 1: Training Final Auxiliary Models")
    print(f"{'─'*60}")
    
    # Train final auxiliary models
    final_aux_models = {}
    for aux_name in config.auxiliary_targets:
        is_classification = (aux_name == 'IO_LULC')
        
        if is_classification:
            print(f"🔧 Training final {aux_name} model (classification)...")
            aux_params = final_base_params.copy()
            
            train_mask = ~np.isnan(y_train_final_all[aux_name])
            unique_classes = np.unique(y_train_final_all[aux_name][train_mask])
            n_classes = len(unique_classes)
            
            # Create label mapping
            class_to_idx = {cls: idx for idx, cls in enumerate(unique_classes)}
            y_train_mapped = np.array([class_to_idx[int(c)] for c in y_train_final_all[aux_name][train_mask]])
            
            if n_classes <= 2:
                aux_params['objective'] = 'binary:logistic'
                aux_params['eval_metric'] = 'logloss'
                final_aux_model = xgb.XGBClassifier(**aux_params)
            else:
                aux_params['objective'] = 'multi:softprob'
                aux_params['eval_metric'] = 'mlogloss'
                aux_params['num_class'] = int(n_classes)
                final_aux_model = xgb.XGBClassifier(**aux_params)
            
            final_aux_model.fit(X_train_final[train_mask], y_train_mapped)
            
            # Store class mapping
            final_aux_model.class_mapping_ = class_to_idx
            final_aux_models[aux_name] = final_aux_model
            print(f"  ✅ {aux_name} trained ({n_classes} classes)")
            
        else:
            print(f"🔧 Training final {aux_name} model (regression)...")
            aux_params = final_base_params.copy()
            aux_params['objective'] = 'reg:squarederror'
            aux_params['tree_method'] = config.tree_method
            aux_params['eval_metric'] = 'rmse'
            
            final_aux_model = xgb.XGBRegressor(**aux_params)
            
            train_mask = ~np.isnan(y_train_final_all[aux_name])
            final_aux_model.fit(
                X_train_final[train_mask],
                y_train_final_all[aux_name][train_mask]
            )
            final_aux_models[aux_name] = final_aux_model
            print(f"  ✅ {aux_name} trained")
    
    print(f"\n{'─'*60}")
    print(f"STEP 2: Generating Final Auxiliary Predictions")
    print(f"{'─'*60}")
    
    # Generate auxiliary predictions for final training
    X_train_final_aug = X_train_final.copy()
    
    for aux_name in config.auxiliary_targets:
        aux_model = final_aux_models[aux_name]
        is_classification = (aux_name == 'IO_LULC')
        
        train_aux_feature = np.zeros(len(X_train_final), dtype=np.float32)
        train_mask = ~np.isnan(y_train_final_all[aux_name])
        
        if is_classification:
            if train_mask.sum() > 0:
                train_aux_feature[train_mask] = aux_model.predict(X_train_final[train_mask]).astype(np.float32)
        else:
            if train_mask.sum() > 0:
                train_aux_feature[train_mask] = aux_model.predict(X_train_final[train_mask])
        
        X_train_final_aug = np.column_stack([X_train_final_aug, train_aux_feature])
        print(f"  ✅ {aux_name}: prediction added (mean={train_aux_feature.mean():.4f})")
    
    print(f"\n{'─'*60}")
    print(f"STEP 3: Training Final Flood Model on Augmented Features")
    print(f"{'─'*60}")
    
    # Train final flood model on augmented features
    print(f"🎯 Training final flood model on {X_train_final_aug.shape[1]} features...")
    final_flood_params = final_base_params.copy()
    final_flood_params['objective'] = config.objective
    final_flood_params['tree_method'] = config.tree_method
    final_flood_params['eval_metric'] = config.eval_metric
    
    final_flood_model = xgb.XGBClassifier(**final_flood_params)
    final_flood_model.fit(X_train_final_aug, y_train_final_all['flood'])
    print(f"  ✅ Final flood model trained")
    
    # Evaluate on test set if available
    if test_ws:
        print(f"\n{'─'*60}")
        print(f"EVALUATING ON TEST SET")
        print(f"{'─'*60}")
        
        X_test, y_test_all = load_data_streaming_multitask(
            dataset, test_ws, feature_cols, config.auxiliary_targets,
            means, stds, config, mode='test'
        )
        
        if len(X_test) > 0:
            # Generate auxiliary predictions for test set
            X_test_aug = X_test.copy()
            
            for aux_name in config.auxiliary_targets:
                aux_model = final_aux_models[aux_name]
                is_classification = (aux_name == 'IO_LULC')
                
                test_aux_feature = np.zeros(len(X_test), dtype=np.float32)
                test_mask = ~np.isnan(y_test_all[aux_name])
                
                if test_mask.sum() > 0:
                    if is_classification:
                        test_aux_feature[test_mask] = aux_model.predict(X_test[test_mask]).astype(np.float32)
                    else:
                        test_aux_feature[test_mask] = aux_model.predict(X_test[test_mask])
                
                X_test_aug = np.column_stack([X_test_aug, test_aux_feature])
            
            # Evaluate flood task
            y_test_flood_proba = final_flood_model.predict_proba(X_test_aug)[:, 1]
            test_flood_metrics = compute_metrics_classification(y_test_all['flood'], y_test_flood_proba, final_threshold)
            summary_metrics['flood']['test'] = test_flood_metrics
            print(f"🎯 Test Flood - MCC={test_flood_metrics['mcc']:.4f}, F1={test_flood_metrics['f1']:.4f}, "
                  f"PR-AUC={test_flood_metrics['pr_auc']:.4f}")
            
            # Evaluate auxiliary tasks
            for aux_name in config.auxiliary_targets:
                is_classification = (aux_name == 'IO_LULC')
                test_mask = ~np.isnan(y_test_all[aux_name])
                
                if test_mask.sum() < 10:
                    continue
                
                if is_classification:
                    y_test_aux_pred = final_aux_models[aux_name].predict(X_test[test_mask])
                    test_acc = accuracy_score(
                        y_test_all[aux_name][test_mask].astype(int),
                        y_test_aux_pred
                    )
                    summary_metrics[aux_name]['test'] = {'accuracy': test_acc}
                    print(f"🔧 Test {aux_name} - Accuracy={test_acc:.4f}")
                else:
                    y_test_aux_pred = final_aux_models[aux_name].predict(X_test[test_mask])
                    test_aux_metrics = compute_metrics_regression(y_test_all[aux_name][test_mask], y_test_aux_pred)
                    summary_metrics[aux_name]['test'] = test_aux_metrics
                    print(f"🔧 Test {aux_name} - RMSE={test_aux_metrics['rmse']:.4f}, R²={test_aux_metrics['r2']:.4f}")
    
    # Create augmented feature list for interpretability
    feature_cols_augmented = feature_cols + [f"{aux}_pred" for aux in config.auxiliary_targets]
    
    # Feature importance from flood model (includes auxiliary predictions)
    feature_importance_df = pd.DataFrame({
        'feature': feature_cols_augmented,
        'importance': final_flood_model.feature_importances_
    }).sort_values('importance', ascending=False)
    
    print(f"\n{'─'*60}")
    print(f"TOP 10 MOST IMPORTANT FEATURES (including auxiliary predictions)")
    print(f"{'─'*60}")
    for idx, row in feature_importance_df.head(10).iterrows():
        is_aux = '_pred' in row['feature']
        marker = '🔧' if is_aux else '📊'
        print(f"  {marker} {row['feature']}: {row['importance']:.4f}")
    
    results = {
        'model_flood': final_flood_model,
        'models_auxiliary': final_aux_models,
        'threshold': final_threshold,
        'feature_cols': feature_cols,  # Base features (53)
        'feature_cols_augmented': feature_cols_augmented,  # Base + auxiliary predictions (57)
        'auxiliary_targets': config.auxiliary_targets,
        'means': means,
        'stds': stds,
        'config': config,
        'test_watersheds': test_ws,
        'metrics_summary': summary_metrics,
        'feature_importance': feature_importance_df,
        'best_hyperparameters': best_params,
        'cv_results_flood': folds_df_flood,
        'cv_results_auxiliary': folds_df_aux
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
        print(f"\n💾 Model saved to {save_path}")
    else:
        import shutil
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        shutil.move(temp_path, save_path)
        print(f"\n💾 Model saved to {save_path}")
    
    os.remove(temp_path) if os.path.exists(temp_path) else None
    
    # Create diagnostic plots in same folder
    create_diagnostic_plots_multilearn(results, model_folder)
    
    print("\n" + "=" * 70)
    print("✅ CASCADING MULTITASK LEARNING COMPLETE")
    print("=" * 70)
    print(f"📊 Architecture: Auxiliary models → Feature augmentation → Flood model")
    print(f"🔢 Feature space: {len(feature_cols)} base + {len(config.auxiliary_targets)} auxiliary = {len(feature_cols_augmented)} total")
    print(f"🎯 Primary task MCC: {summary_metrics['flood'].get('val', {}).get('mcc', 'N/A'):.4f} (validation)")
    
    return results


# ======================== Diagnostic Plots - MULTITASK ========================

def create_diagnostic_plots_multilearn(model_results: Dict, output_folder: str):
    """Create comprehensive diagnostic plots for multitask model."""
    
    print("\nCreating diagnostic plots...")
    
    feature_importance = model_results['feature_importance']
    cv_results_flood = model_results.get('cv_results_flood')
    cv_results_aux = model_results.get('cv_results_auxiliary', {})
    
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
    
    # 1. Feature Importance (Top 20) - from flood model
    fig, ax = plt.subplots(figsize=(12, 10))
    top_features = feature_importance.head(20)
    
    ax.barh(range(len(top_features)), top_features['importance'], color='steelblue')
    ax.set_yticks(range(len(top_features)))
    ax.set_yticklabels(top_features['feature'])
    ax.set_xlabel('Feature Importance (Gain)', fontsize=12)
    ax.set_title('XGBoostMultiLearn Feature Importance - Flood Task (Top 20)', fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)
    
    plt.tight_layout()
    importance_path = f"{output_folder}/xgboostmultilearn_feature_importance.png"
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
    dist_path = f"{output_folder}/xgboostmultilearn_importance_distribution.png"
    save_and_upload_plot(fig, dist_path)
    
    # 3. CV Performance - Flood Task
    if cv_results_flood is not None and not cv_results_flood.empty:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        metrics_to_plot = ['mcc', 'f1', 'precision', 'recall']
        titles = ['MCC (Primary Metric)', 'F1 Score', 'Precision', 'Recall']
        
        for idx, (metric, title) in enumerate(zip(metrics_to_plot, titles)):
            ax = axes[idx // 2, idx % 2]
            
            train_vals = cv_results_flood[cv_results_flood['set'] == 'train'][metric].values
            val_vals = cv_results_flood[cv_results_flood['set'] == 'val'][metric].values
            
            x = np.arange(len(train_vals))
            width = 0.35
            
            ax.bar(x - width/2, train_vals, width, label='Train', color='steelblue', alpha=0.8)
            ax.bar(x + width/2, val_vals, width, label='Val', color='coral', alpha=0.8)
            
            ax.set_xlabel('Fold', fontsize=11)
            ax.set_ylabel(metric.upper(), fontsize=11)
            ax.set_title(f'Flood Task - {title}', fontsize=12, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f'{i+1}' for i in range(len(train_vals))])
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        cv_path = f"{output_folder}/xgboostmultilearn_cv_performance_flood.png"
        save_and_upload_plot(fig, cv_path)
    
    # 4. CV Performance - Auxiliary Tasks
    for aux_name, aux_df in cv_results_aux.items():
        if aux_df is None or aux_df.empty:
            continue
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Determine metric to plot
        if 'rmse' in aux_df.columns:
            metric = 'rmse'
            ylabel = 'RMSE'
            title_suffix = 'RMSE'
        elif 'accuracy' in aux_df.columns:
            metric = 'accuracy'
            ylabel = 'Accuracy'
            title_suffix = 'Accuracy'
        else:
            continue
        
        train_vals = aux_df[aux_df['set'] == 'train'][metric].values
        val_vals = aux_df[aux_df['set'] == 'val'][metric].values
        
        x = np.arange(len(train_vals))
        width = 0.35
        
        ax.bar(x - width/2, train_vals, width, label='Train', color='seagreen', alpha=0.8)
        ax.bar(x + width/2, val_vals, width, label='Val', color='orange', alpha=0.8)
        
        ax.set_xlabel('Fold', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(f'{aux_name} Task - {title_suffix}', fontsize=12, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels([f'{i+1}' for i in range(len(train_vals))])
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        aux_path = f"{output_folder}/xgboostmultilearn_cv_performance_{aux_name}.png"
        save_and_upload_plot(fig, aux_path)
    
    print("✅ Diagnostic plots complete!")


# ======================== Visualization Class - MULTITASK ========================

class XGBoostMultiLearnFloodVisualizer:
    """Generate visualizations for XGBoostMultiLearn model."""
    
    def __init__(self, model_results: Dict, parquet_uri: str, training_folder: str):
        self.model_results = model_results
        self.model_flood = model_results['model_flood']
        self.models_auxiliary = model_results['models_auxiliary']
        self.threshold = model_results['threshold']
        self.feature_cols = model_results['feature_cols']
        self.auxiliary_targets = model_results['auxiliary_targets']
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
    
    def generate_all_outputs(self, output_folder: Optional[str] = None,
                            select_dates: Optional[List[str]] = None,
                            generate_forecasts_for_all_subwatersheds: bool = False):
        """Generate all visualization outputs - same interface as single-task version."""
        
        from processDataToParquet import build_flood_training_parquet_v3
        
        if output_folder is None:
            output_folder = f"{self.training_folder}/assessXGBoostMultiLearnModel"
        
        print("\n" + "="*60)
        print("GENERATING XGBOOSTMULTILEARN MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        
        if generate_forecasts_for_all_subwatersheds:
            print(f"Mode: All watersheds")
        else:
            print(f"Test watershed: {self.test_watershed}")
        
        if select_dates:
            print(f"Selected dates: {select_dates}")
        
        print(f"Output folder: {output_folder}")
        
        # ✅ Determine mode based on select_dates
        is_prediction_mode = (select_dates is not None)
        
        if is_prediction_mode:
            # ✅ PREDICTION MODE: Generate new parquet for specific dates
            from processDataToParquet import build_flood_training_parquet_v3
            
            # Parse training_folder
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
            
            print(f"\n📝 Generating prediction parquet for selected dates: {select_dates}")
            try:
                new_parquet_uri = build_flood_training_parquet_v3(
                    base_lat_4=base_lat_4,
                    base_lon_4=base_lon_4,
                    bucket=bucket,
                    base_prefix=base_prefix,
                    select_dates=select_dates
                )
                self.parquet_uri = new_parquet_uri  # ✅ UPDATE to use new parquet
                print(f"✅ Prediction parquet generated: {self.parquet_uri}")
            except Exception as e:
                raise RuntimeError(f"Failed to generate prediction parquet: {str(e)}")
        else:
            # ✅ TRAINING/ASSESSMENT MODE: Use pre-existing parquet from training
            # The parquet already contains ALL watersheds (including test watershed)
            # Test watershed was only excluded from TRAINING, not from the parquet file
            print(f"\n📂 Using pre-existing training parquet (contains all watersheds): {self.parquet_uri}")
            
            # ✅ Verify parquet exists and contains data
            try:
                dataset = open_parquet_dataset(self.parquet_uri)
                all_watersheds = get_unique_watersheds(dataset)
                print(f"   Parquet contains {len(all_watersheds)} watersheds: {all_watersheds}")
                
                if self.test_watershed not in all_watersheds:
                    print(f"\n⚠️  WARNING: Test watershed '{self.test_watershed}' not found in parquet!")
                    print(f"   Available watersheds: {all_watersheds}")
                    print(f"   This suggests the parquet file may be incomplete or corrupted.")
            except Exception as e:
                print(f"\n⚠️  WARNING: Could not verify parquet contents: {str(e)}")
        
        # Only flood predictions are visualized (same as single-task)
        # Auxiliary task outputs are used during training but not visualized
        
        if select_dates is None and not generate_forecasts_for_all_subwatersheds:
            create_diagnostic_plots_multilearn(self.model_results, output_folder)
        
        # Determine watersheds to process
        if generate_forecasts_for_all_subwatersheds:
            dataset = open_parquet_dataset(self.parquet_uri)
            watersheds_to_process = get_unique_watersheds(dataset)
            watershed_suffix = "all_watersheds"
        else:
            watersheds_to_process = [self.test_watershed]
            watershed_suffix = self.test_watershed
        
        # Generate predictions (only flood task)
        if generate_forecasts_for_all_subwatersheds:
            predictions, confusions, per_date_metrics = self._generate_predictions_multi_watershed(
                watersheds_to_process, select_dates, is_prediction_mode
            )
        else:
            predictions, confusions, per_date_metrics = self._generate_predictions(
                select_dates, is_prediction_mode
            )
        
        if not predictions:
            print("\n⚠️  No predictions generated - no data found for specified parameters")
            print(f"   Parquet URI: {self.parquet_uri}")
            print(f"   Test watershed: {self.test_watershed}")
            print(f"   Try regenerating the parquet file or check watershed naming")
            return {
                'prediction_tiff': None,
                'confusion_matrix_tiff': None if not is_prediction_mode else None,
                'png_visualizations': [],
                'csv_outputs': {}
            }
       
        # Get reference dataset
        flood_tiff_path = (
            f"{self.training_folder}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        ref_ds = gdal.Open(flood_tiff_path.replace('s3://', '/vsis3/'), gdal.GA_ReadOnly)
        
        dates = sorted(predictions.keys())
        
        # Create output files
        if is_prediction_mode:
            pred_tiff = f"{output_folder}/xgboostmultilearn_predictions_{watershed_suffix}.tif"
            self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "XGBoostMultiLearn_Prediction")
            
            png_paths = self._create_prediction_pngs(predictions, dates, per_date_metrics, output_folder)
            csv_paths = self._write_csvs(output_folder, per_date_metrics)
            
            print("\n✅ XGBoostMultiLearn prediction visualization complete!")
            
            return {
                'prediction_tiff': pred_tiff,
                'png_visualizations': png_paths,
                'csv_outputs': csv_paths
            }
        else:
            pred_tiff = f"{output_folder}/xgboostmultilearn_predictions_{watershed_suffix}.tif"
            cm_tiff = f"{output_folder}/xgboostmultilearn_confusion_matrix_{watershed_suffix}.tif"
            
            self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "XGBoostMultiLearn_Prediction")
            self._write_multiband_tiff(ref_ds, confusions, dates, cm_tiff, "XGBoostMultiLearn_CM")
            
            png_paths = self._create_pngs(confusions, dates, per_date_metrics, output_folder)
            csv_paths = self._write_csvs(output_folder, per_date_metrics)
            
            print("\n✅ XGBoostMultiLearn visualization complete!")
            
            return {
                'prediction_tiff': pred_tiff,
                'confusion_matrix_tiff': cm_tiff,
                'png_visualizations': png_paths,
                'csv_outputs': csv_paths
            }
    
    # Rest of methods are identical to single-task version but use model_flood
    # I'll include the key methods with minimal changes
    
    def _generate_predictions(self, select_dates, is_prediction_mode):
        """Generate predictions using CASCADING architecture (auxiliary → flood)."""
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
        df_all['Date_Normalized'] = df_all['Date'].apply(self._normalize_date_format)
        
        if select_dates is not None:
            select_dates_normalized = [self._normalize_date_format(d) for d in select_dates]
            df_all = df_all[df_all['Date_Normalized'].isin(select_dates_normalized)]
            
            if df_all.empty:
                return {}, {}, {}
        
        dates = sorted(df_all['Date_Normalized'].unique())
        H = df_all['PixelRow'].max() + 1
        W = df_all['PixelCol'].max() + 1
        
        predictions = {}
        confusions = {} if not is_prediction_mode else None
        per_date_metrics = {}
        
        for date_str in dates:
            date_data = df_all[df_all['Date_Normalized'] == date_str]
            
            X = date_data[self.feature_cols].to_numpy(dtype=np.float32)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8) if not is_prediction_mode else None
            rows = date_data['PixelRow'].to_numpy()
            cols = date_data['PixelCol'].to_numpy()
            
            valid_mask = ~np.all(np.isnan(X), axis=1)
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            
            if valid_mask.sum() > 0:
                X_valid = X[valid_mask]
                X_valid = np.where(np.isnan(X_valid), self.means, self.means)
                X_valid = (X_valid - self.means) / self.stds
                
                # ✅ CASCADING: Generate auxiliary predictions first
                X_valid_aug = X_valid.copy()
                for aux_name in self.auxiliary_targets:
                    aux_model = self.models_auxiliary.get(aux_name)
                    
                    if aux_model is None:
                        # Use zeros if model unavailable
                        aux_pred = np.zeros(len(X_valid), dtype=np.float32)
                    else:
                        is_classification = (aux_name == 'IO_LULC')
                        if is_classification:
                            aux_pred = aux_model.predict(X_valid).astype(np.float32)
                        else:
                            aux_pred = aux_model.predict(X_valid)
                    
                    # Append auxiliary prediction as feature
                    X_valid_aug = np.column_stack([X_valid_aug, aux_pred])
                
                # ✅ Use augmented features for flood prediction
                y_proba_valid = self.model_flood.predict_proba(X_valid_aug)[:, 1]
                y_pred_valid = (y_proba_valid >= self.threshold).astype(np.uint8)
                
                pred_img[rows[valid_mask], cols[valid_mask]] = y_pred_valid
            
            predictions[date_str] = pred_img
            
            if not is_prediction_mode and y_true is not None:
                cm_img = np.full((H, W), 255, dtype=np.uint8)
                
                if valid_mask.sum() > 0:
                    y_true_valid = y_true[valid_mask]
                    y_pred_valid = pred_img[rows[valid_mask], cols[valid_mask]]
                    
                    cm_vals = np.where(
                        y_true_valid == 0,
                        np.where(y_pred_valid == 0, 0, 3),
                        np.where(y_pred_valid == 1, 1, 2)
                    )
                    cm_img[rows[valid_mask], cols[valid_mask]] = cm_vals
                
                confusions[date_str] = cm_img
                
                if valid_mask.sum() > 0:
                    X_valid = X[valid_mask]
                    y_true_valid = y_true[valid_mask]
                    X_valid = np.where(np.isnan(X_valid), self.means, X_valid)
                    X_valid = (X_valid - self.means) / self.stds
                    
                    # Generate auxiliary predictions for metrics
                    X_valid_aug = X_valid.copy()
                    for aux_name in self.auxiliary_targets:
                        aux_model = self.models_auxiliary.get(aux_name)
                        if aux_model is None:
                            aux_pred = np.zeros(len(X_valid), dtype=np.float32)
                        else:
                            is_classification = (aux_name == 'IO_LULC')
                            if is_classification:
                                aux_pred = aux_model.predict(X_valid).astype(np.float32)
                            else:
                                aux_pred = aux_model.predict(X_valid)
                        X_valid_aug = np.column_stack([X_valid_aug, aux_pred])
                    
                    y_proba_valid = self.model_flood.predict_proba(X_valid_aug)[:, 1]
                    y_pred_valid = (y_proba_valid >= self.threshold).astype(np.uint8)
                    
                    metrics = compute_metrics_classification(y_true_valid, y_proba_valid, self.threshold)
                    metrics['date'] = date_str
                    metrics['valid_pixel_count'] = int(valid_mask.sum())
                    per_date_metrics[date_str] = metrics
            else:
                if valid_mask.sum() > 0:
                    pred_fraction = pred_img[pred_img != 255].mean() if (pred_img != 255).any() else 0.0
                    per_date_metrics[date_str] = {
                        'date': date_str,
                        'pred_pos_frac': pred_fraction,
                        'valid_pixel_count': int(valid_mask.sum())
                    }
        
        return predictions, confusions if confusions is not None else {}, per_date_metrics
        
    def _create_pngs(self, confusions, dates, metrics, output_folder):
        """Create confusion matrix PNG visualizations matching mlpFilmFloods.py style."""
        png_paths = []
        
        # Limit number of PNGs if configured
        if self.config.max_png_dates and len(dates) > self.config.max_png_dates:
            dates_to_plot = sorted(dates)[:self.config.max_png_dates]
            print(f"  Limiting to {self.config.max_png_dates} PNGs (of {len(dates)} dates)")
        else:
            dates_to_plot = dates
        
        for date_str in dates_to_plot:
            cm = confusions.get(date_str)
            if cm is None:
                continue
            
            # Skip dates with almost no valid pixels
            if np.count_nonzero(cm != 255) < 100:
                continue
            
            m = metrics.get(date_str, {})
            
            fig, ax = plt.subplots(figsize=(12, 10))
            
            # ✅ Use ListedColormap like mlpFilmFloods.py
            cm_disp = cm.astype(float)
            cm_disp[cm == 255] = np.nan  # ✅ NaN for nodata (renders as white)
            
            # ✅ Color scheme: 0=TN (green), 1=TP (blue), 2=FN (orange), 3=FP (red)
            cmap = matplotlib.colors.ListedColormap(['green', 'blue', 'orange', 'red'])
            ax.imshow(cm_disp, cmap=cmap, vmin=0, vmax=3, interpolation='nearest')
            
            # ✅ Simple legend matching mlpFilmFloods.py
            legend_elements = [
                mpatches.Patch(color='green', label='True Negative'),
                mpatches.Patch(color='blue', label='True Positive'),
                mpatches.Patch(color='orange', label='False Negative'),
                mpatches.Patch(color='red', label='False Positive')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
            
            # ✅ Title matching mlpFilmFloods.py format
            ax.set_title(f'XGBoostMultiLearn Confusion Matrix - {date_str}', fontsize=14, fontweight='bold')
            
            # ✅ Text annotation matching mlpFilmFloods.py format and position
            text = (
                f"Date: {date_str}\n"
                f"Model: XGBoostMultiLearn\n"
                f"F1: {m.get('f1', 0):.3f}\n"
                f"Precision: {m.get('precision', 0):.3f}\n"
                f"Recall: {m.get('recall', 0):.3f}\n"
                f"ROC AUC: {m.get('roc_auc', 0):.3f}"
            )
            ax.text(0.02, 0.98, text, transform=ax.transAxes, va='top',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
                    fontsize=11)
            
            ax.axis('off')
            plt.tight_layout()
            
            # Save PNG
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_path = tmp.name
            
            fig.savefig(temp_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            
            # Upload to S3
            png_path = f"{output_folder}/xgboostmultilearn_cm_{date_str}.png"
            
            if png_path.startswith('s3://'):
                bucket, key = png_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(temp_path, bucket, key)
                os.remove(temp_path)
            else:
                import shutil
                os.makedirs(os.path.dirname(png_path), exist_ok=True)
                shutil.move(temp_path, png_path)
            
            png_paths.append(png_path)
        
        print(f"  ✅ Created {len(png_paths)} confusion matrix PNGs")
        return png_paths
    
    
    def _write_csvs(self, output_folder, per_date_metrics):
        """Write CSV outputs matching mlpFilmFloods.py format."""
        paths = {}
        
        # ✅ Per-date metrics - matching mlpFilmFloods.py column order exactly
        if per_date_metrics:
            df = pd.DataFrame.from_dict(per_date_metrics, orient='index')
            
            # ✅ Ensure consistent column order matching mlpFilmFloods.py
            column_order = ['date', 'f1', 'accuracy', 'precision', 'recall', 'roc_auc', 
                           'pr_auc', 'kappa', 'mcc']
            
            # Rename columns to match mlpFilmFloods.py format (capitalize first letter)
            df = df.rename(columns={
                'f1': 'F1',
                'accuracy': 'Accuracy',
                'precision': 'Precision',
                'recall': 'Recall',
                'roc_auc': 'ROC_AUC',
                'pr_auc': 'PR_AUC',
                'kappa': 'Kappa',
                'mcc': 'MCC'
            })
            
            # Add actual_pos_frac and pred_pos_frac if not present
            if 'actual_pos_frac' not in df.columns:
                df['actual_pos_frac'] = df['prevalence'] if 'prevalence' in df.columns else np.nan
            if 'pred_pos_frac' not in df.columns:
                # Calculate from confusion matrix
                df['pred_pos_frac'] = (df['tp'] + df['fp']) / (df['tp'] + df['fp'] + df['tn'] + df['fn'])
            
            # Final column order matching mlpFilmFloods.py
            final_columns = ['date', 'F1', 'Accuracy', 'Precision', 'Recall', 'ROC_AUC', 
                            'PR_AUC', 'Kappa', 'MCC', 'actual_pos_frac', 'pred_pos_frac']
            df = df[[col for col in final_columns if col in df.columns]]
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
                df.to_csv(tmp_csv, index=False)
            
            s3_path = f"{output_folder}/xgboostmultilearn_per_date_metrics_{self.test_watershed}.csv"
            
            if s3_path.startswith('s3://'):
                bucket, key = s3_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(tmp_csv, bucket, key)
                os.remove(tmp_csv)
            else:
                import shutil
                os.makedirs(os.path.dirname(s3_path), exist_ok=True)
                shutil.move(tmp_csv, s3_path)
            
            paths['per_date_metrics'] = s3_path
            print(f"  ✅ Created: {s3_path}")
        
        # ✅ Summary metrics (train/val/test), aligned to mlpFilmFloods.py schema
        metrics_summary = self.model_results.get("metrics_summary", {})
        if metrics_summary:
            rows = []
            for split in ["train", "val", "test"]:
                m = metrics_summary.get('flood', {}).get(split)
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
                    "tp": m.get("tp", 0),
                    "fp": m.get("fp", 0),
                    "tn": m.get("tn", 0),
                    "fn": m.get("fn", 0),
                    "tpr": m.get("tpr", np.nan),
                    "tnr": m.get("tnr", np.nan),
                })
            if rows:
                df_sum = pd.DataFrame(rows)
                
                with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
                    tmp_csv = tmp.name
                    df_sum.to_csv(tmp_csv, index=False)
                
                s3_path = f"{output_folder}/xgboostmultilearn_metrics_summary.csv"
                
                if s3_path.startswith('s3://'):
                    bucket, key = s3_path.replace('s3://', '').split('/', 1)
                    self.s3_client.upload_file(tmp_csv, bucket, key)
                    os.remove(tmp_csv)
                else:
                    import shutil
                    os.makedirs(os.path.dirname(s3_path), exist_ok=True)
                    shutil.move(tmp_csv, s3_path)
                
                paths['metrics_summary'] = s3_path
                print(f"  ✅ Created: {s3_path}")
        
        # ✅ Feature importance (real XGBoost importances, not placeholder)
        feature_importance = self.model_results.get('feature_importance')
        if feature_importance is not None and not feature_importance.empty:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
                tmp_csv = tmp.name
                feature_importance.to_csv(tmp_csv, index=False)
            
            s3_path = f"{output_folder}/xgboostmultilearn_feature_importances.csv"
            
            if s3_path.startswith('s3://'):
                bucket, key = s3_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(tmp_csv, bucket, key)
                os.remove(tmp_csv)
            else:
                import shutil
                os.makedirs(os.path.dirname(s3_path), exist_ok=True)
                shutil.move(tmp_csv, s3_path)
            
            paths['feature_importances'] = s3_path
            print(f"  ✅ Created: {s3_path}")
        
        return paths
    
    
    def _write_multiband_tiff(self, ref_ds, arrays, dates, output_path, prefix):
        """Write multiband TIFF matching mlpFilmFloods.py format."""
        if not arrays or not dates:
            return
        
        # Get dimensions from reference
        H = ref_ds.RasterYSize
        W = ref_ds.RasterXSize
        
        # Create output dataset
        driver = gdal.GetDriverByName('GTiff')
        
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
            temp_path = tmp.name
        
        # ✅ Use DEFLATE compression like mlpFilmFloods.py
        out_ds = driver.Create(
            temp_path,
            W, H,
            len(dates),
            gdal.GDT_Byte,
            options=['COMPRESS=DEFLATE', 'TILED=YES']  # ✅ Changed from LZW to DEFLATE
        )
        
        # Copy georeferencing
        out_ds.SetGeoTransform(ref_ds.GetGeoTransform())
        out_ds.SetProjection(ref_ds.GetProjection())
        
        # Write bands
        for band_idx, date_str in enumerate(dates, start=1):
            band = out_ds.GetRasterBand(band_idx)
            band.WriteArray(arrays[date_str])
            band.SetDescription(f"{prefix}_{date_str}")
            band.SetNoDataValue(255)
            band.FlushCache()
        
        out_ds.FlushCache()
        out_ds = None
        
        # Upload to S3 if needed
        if output_path.startswith('s3://'):
            bucket, key = output_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(temp_path, bucket, key)
            os.remove(temp_path)
        else:
            import shutil
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            shutil.move(temp_path, output_path)
        
        print(f"  ✅ Created: {output_path}")

    
    def _create_prediction_pngs(self, predictions, dates, metrics, output_folder):
        """Create prediction PNG visualizations (for prediction mode)."""
        png_paths = []
        
        # Limit number of PNGs if configured
        if self.config.max_png_dates and len(dates) > self.config.max_png_dates:
            dates_to_plot = sorted(dates)[:self.config.max_png_dates]
            print(f"  Limiting to {self.config.max_png_dates} PNGs (of {len(dates)} dates)")
        else:
            dates_to_plot = dates
        
        # Color map for predictions
        # 0=no flood (blue), 1=flood (red), 255=nodata (white)
        colors = {
            0: (0.2, 0.2, 0.8),  # No flood - blue
            1: (0.9, 0.1, 0.1),  # Flood - red
            255: (1.0, 1.0, 1.0)  # nodata - white
        }
        
        for date_str in dates_to_plot:
            pred_img = predictions[date_str]
            metrics_date = metrics.get(date_str, {})
            
            fig, ax = plt.subplots(figsize=(10, 10))
            
            # Create RGB image
            H, W = pred_img.shape
            rgb = np.ones((H, W, 3))
            for val, color in colors.items():
                mask = (pred_img == val)
                rgb[mask] = color
            
            ax.imshow(rgb, interpolation='nearest')
            ax.set_title(f'Flood Prediction - {date_str}', fontsize=14, fontweight='bold')
            ax.axis('off')
            
            # Add legend
            n_flood = np.sum(pred_img == 1)
            n_no_flood = np.sum(pred_img == 0)
            legend_elements = [
                mpatches.Patch(color=colors[1], label=f"Predicted Flood: {n_flood:,}"),
                mpatches.Patch(color=colors[0], label=f"Predicted No-Flood: {n_no_flood:,}"),
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
            
            # Add summary text
            pred_frac = metrics_date.get('pred_pos_frac', 0)
            n_valid = metrics_date.get('valid_pixel_count', 0)
            
            textstr = f'Flood Fraction: {pred_frac:.3f}\nValid Pixels: {n_valid:,}'
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=11,
                    verticalalignment='top', bbox=props)
            
            plt.tight_layout()
            
            # Save PNG
            png_path = f"{output_folder}/xgboostmultilearn_pred_{date_str}.png"
            
            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_path = tmp.name
            
            fig.savefig(temp_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            
            # Upload to S3 if needed
            if png_path.startswith('s3://'):
                bucket, key = png_path.replace('s3://', '').split('/', 1)
                self.s3_client.upload_file(temp_path, bucket, key)
                os.remove(temp_path)
            else:
                import shutil
                os.makedirs(os.path.dirname(png_path), exist_ok=True)
                shutil.move(temp_path, png_path)
            
            png_paths.append(png_path)
        
        print(f"  ✅ Created {len(png_paths)} prediction PNGs")
        return png_paths
    
    def create_diagnostic_plots_multilearn(model_results: Dict, output_folder: str):
        """Create comprehensive diagnostic plots matching mlpFilmFloods.py style."""
        
        print("\nCreating diagnostic plots...")
        
        feature_importance = model_results['feature_importance']
        cv_results_flood = model_results.get('cv_results_flood')
        cv_results_aux = model_results.get('cv_results_auxiliary', {})
        
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
        ax.set_title('XGBoostMultiLearn Feature Importance (Top 20)', fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        save_and_upload_plot(fig, f"{output_folder}/xgboostmultilearn_feature_importance.png")
        
        # 2. ✅ CV Performance - matching mlpFilmFloods.py 2×2 grid
        if cv_results_flood is not None and not cv_results_flood.empty:
            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            
            metrics_to_plot = ['mcc', 'f1', 'precision', 'recall']
            titles = ['MCC', 'F1 Score', 'Precision', 'Recall']
            
            for idx, (metric, title) in enumerate(zip(metrics_to_plot, titles)):
                ax = axes[idx // 2, idx % 2]
                
                train_vals = cv_results_flood[cv_results_flood['set'] == 'train'][metric].values
                val_vals = cv_results_flood[cv_results_flood['set'] == 'val'][metric].values
                
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
            save_and_upload_plot(fig, f"{output_folder}/xgboostmultilearn_cv_performance.png")
        
        # 3. ✅ NEW: Threshold Distribution (matching mlpFilmFloods.py)
        if cv_results_flood is not None and 'threshold' in cv_results_flood.columns:
            fig, ax = plt.subplots(figsize=(10, 5))
            thresholds = cv_results_flood[cv_results_flood['set'] == 'val']['threshold'].values
            ax.hist(thresholds, bins=30, color='steelblue', alpha=0.8, edgecolor='black')
            ax.set_title('Distribution of CV Thresholds (MCC-optimal)', fontsize=14, fontweight='bold')
            ax.set_xlabel('Threshold', fontsize=12)
            ax.set_ylabel('Count', fontsize=12)
            ax.grid(alpha=0.3)
            ax.axvline(np.mean(thresholds), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(thresholds):.3f}')
            ax.legend()
            plt.tight_layout()
            save_and_upload_plot(fig, f"{output_folder}/xgboostmultilearn_threshold_distribution.png")
        
        # 4. Auxiliary task performance (keep as-is, unique to multitask)
        for aux_name, aux_df in cv_results_aux.items():
            if aux_df is None or aux_df.empty:
                continue
            
            fig, ax = plt.subplots(figsize=(10, 6))
            
            if 'rmse' in aux_df.columns:
                metric = 'rmse'
                ylabel = 'RMSE'
                title_suffix = 'RMSE'
            elif 'accuracy' in aux_df.columns:
                metric = 'accuracy'
                ylabel = 'Accuracy'
                title_suffix = 'Accuracy'
            else:
                continue
            
            train_vals = aux_df[aux_df['set'] == 'train'][metric].values
            val_vals = aux_df[aux_df['set'] == 'val'][metric].values
            
            x = np.arange(len(train_vals))
            width = 0.35
            
            ax.bar(x - width/2, train_vals, width, label='Train', color='seagreen', alpha=0.8)
            ax.bar(x + width/2, val_vals, width, label='Val', color='orange', alpha=0.8)
            
            ax.set_xlabel('Fold', fontsize=11)
            ax.set_ylabel(ylabel, fontsize=11)
            ax.set_title(f'{aux_name} Task - {title_suffix}', fontsize=12, fontweight='bold')
            ax.set_xticks(x)
            ax.set_xticklabels([f'{i+1}' for i in range(len(train_vals))])
            ax.legend()
            ax.grid(axis='y', alpha=0.3)
            
            plt.tight_layout()
            save_and_upload_plot(fig, f"{output_folder}/xgboostmultilearn_cv_performance_{aux_name}.png")
        
        print("✅ Diagnostic plots complete!")
    
    
    
    def _normalize_date_format(self, date_value) -> str:
        """Normalize date formats to YYYYMMDD."""
        if isinstance(date_value, str):
            # Remove any non-digit characters
            digits = ''.join(c for c in date_value if c.isdigit())
            if len(digits) == 8:
                return digits
            return date_value
        elif isinstance(date_value, (int, np.integer)):
            return str(date_value)
        else:
            return str(date_value)
    
    
    def _generate_predictions_multi_watershed(self, watersheds, select_dates, is_prediction_mode):
        """Generate predictions for multiple watersheds (not yet implemented)."""
        print("  ⚠️  Multi-watershed prediction not yet implemented")
        return {}, {}, {}





def load_model_from_s3(model_path: str) -> Dict:
    """Load trained XGBoostMultiLearn model from S3 pickle file."""
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

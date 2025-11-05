"""
Random Forest flood prediction with LOSO cross-validation (strict).

Simplifications:
- Always requires sufficient watersheds for true LOSO CV:
  - At least 3 total watersheds: 1 held-out test + >=2 non-test for LOSO folds
  - If not satisfied, the routine raises a clear error and exits (no internal split fallback)
- Features = all numeric columns except Target, Date, Watershed_Loc, PixelRow, PixelCol
- Normalization = mean/std computed over all non-test watersheds
- Handles class imbalance with optional negative downsampling + balanced class weights
- Reports F1 (primary), accuracy, precision, recall, ROC AUC, PR AUC, kappa, MCC, and confusion matrix details
- Minimal, robust PyArrow + S3 helpers
"""

import time
import pickle
from typing import List, Tuple, Dict, Optional
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, average_precision_score, confusion_matrix,
    cohen_kappa_score, matthews_corrcoef, precision_recall_curve
)
from sklearn.utils.class_weight import compute_class_weight

import pyarrow as pa
import pyarrow.dataset as ds

import warnings
warnings.filterwarnings("ignore")

import os
import io
import re
import tempfile
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import pandas as pd
import pickle

# Geospatial
from osgeo import gdal
gdal.UseExceptions()

# Image generation
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# AWS
import boto3

# ML
from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score, roc_auc_score
)


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

# Call it once at module load
configure_gdal()


# ----------------------------- Config -----------------------------

class TrainingConfig:
    def __init__(
        self,
        n_estimators: int = 300,
        max_depth: Optional[int] = 20,
        min_samples_split: int = 5,
        min_samples_leaf: int = 2,
        max_features: str = "sqrt",
        criterion: str = "gini",
        min_impurity_decrease: float = 0.0,
        max_leaf_nodes: Optional[int] = None,
        random_state: int = 42,
        neg_pos_ratio: Optional[float] = 20.0,
        batch_size: int = 250_000,
        n_jobs: int = -1,
        use_precalculated_hyperparameters: Optional[str] = None
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.max_features = max_features
        self.criterion = criterion
        self.min_impurity_decrease = min_impurity_decrease
        self.max_leaf_nodes = max_leaf_nodes
        self.random_state = random_state
        self.neg_pos_ratio = neg_pos_ratio
        self.batch_size = batch_size
        self.n_jobs = n_jobs
        self.use_precalculated_hyperparameters = use_precalculated_hyperparameters
        
        # Load and override hyperparameters if S3 URI provided
        if self.use_precalculated_hyperparameters:
            self._load_hyperparameters_from_s3()
    
    def _load_hyperparameters_from_s3(self):
        """Load hyperparameters from S3 JSON and override defaults."""
        try:
            import boto3
            import json
            
            uri = self.use_precalculated_hyperparameters
            if not uri.startswith('s3://'):
                raise ValueError(f"Expected s3:// URI, got: {uri}")
            
            bucket, key = uri[5:].split('/', 1)
            s3 = boto3.client('s3')
            
            obj = s3.get_object(Bucket=bucket, Key=key)
            hyperparams = json.load(obj['Body'])
            
            print("\n" + "="*70)
            print("LOADING PRECALCULATED HYPERPARAMETERS")
            print("="*70)
            print(f"Source: {uri}\n")
            
            # Map of JSON keys to attribute names
            param_mapping = {
                'n_estimators': 'n_estimators',
                'max_depth': 'max_depth',
                'min_samples_split': 'min_samples_split',
                'min_samples_leaf': 'min_samples_leaf',
                'max_features': 'max_features',
                'criterion': 'criterion',
                'min_impurity_decrease': 'min_impurity_decrease',
                'max_leaf_nodes': 'max_leaf_nodes'
            }
            
            overridden = []
            for json_key, attr_name in param_mapping.items():
                if json_key in hyperparams:
                    old_val = getattr(self, attr_name)
                    new_val = hyperparams[json_key]
                    setattr(self, attr_name, new_val)
                    overridden.append((attr_name, old_val, new_val))
            
            if overridden:
                print("Overridden hyperparameters:")
                for name, old, new in overridden:
                    print(f"  {name}: {old} → {new}")
                print()
            
        except Exception as e:
            print(f"\n⚠️ Warning: Failed to load hyperparameters from S3: {e}")
            print("Proceeding with default hyperparameters.\n")


# --------------------- PyArrow + S3 compatibility ---------------------

def _get_s3fs():
    # Try modern API
    try:
        import pyarrow.fs as pafs
        return pafs.S3FileSystem()
    except Exception:
        pass
    # Fallback older API
    try:
        from pyarrow import filesystem as pafs
        return pafs.S3FileSystem()
    except Exception:
        raise RuntimeError("Could not create a PyArrow S3FileSystem. Please install/upgrade pyarrow or use a local Parquet path.")


def _open_dataset(parquet_uri: str):
    if parquet_uri.startswith("s3://"):
        bkt, key = parquet_uri[5:].split("/", 1)
        s3fs = _get_s3fs()
        return ds.dataset(f"{bkt}/{key}", format="parquet", filesystem=s3fs)
    else:
        return ds.dataset(parquet_uri, format="parquet")


def _make_scanner(dataset, columns=None, filter=None, batch_size: Optional[int] = None):
    # Prefer dataset.scanner(...) when available
    try:
        kwargs = {}
        if columns is not None:
            kwargs["columns"] = columns
        if filter is not None:
            kwargs["filter"] = filter
        try:
            if batch_size is not None:
                kwargs["batch_size"] = batch_size
            return dataset.scanner(**kwargs)
        except TypeError:
            kwargs.pop("batch_size", None)
            return dataset.scanner(**kwargs)
    except AttributeError:
        # Fallback to loading fully (last resort)
        table = dataset.to_table(columns=columns, filter=filter)
        class _TableScanner:
            def to_batches(self_inner):
                if len(table) == 0:
                    yield pa.RecordBatch.from_arrays([], [])
                else:
                    for b in table.to_batches(max_chunksize=batch_size or 65536):
                        yield b
        return _TableScanner()


def _scanner_batches(scanner):
    if hasattr(scanner, "to_batches"):
        return scanner.to_batches()
    if hasattr(scanner, "scan_batches"):
        return scanner.scan_batches()
    if hasattr(scanner, "to_table"):
        tbl = scanner.to_table()
        return tbl.to_batches(max_chunksize=65536)
    raise RuntimeError("Unable to iterate scanner batches for this pyarrow version.")


# --------------------------- Data utilities ---------------------------

def get_unique_watersheds(dataset) -> List[str]:
    uniques = set()
    scanner = _make_scanner(dataset, columns=["Watershed_Loc"])
    for batch in _scanner_batches(scanner):
        if batch.num_columns == 0 or batch.num_rows == 0:
            continue
        s = pa.Table.from_batches([batch]).column(0).to_pandas()
        uniques.update(s.dropna().unique().tolist())
    return sorted(uniques)


def get_feature_columns(dataset, extra_exclude: Optional[List[str]] = None) -> List[str]:
    """
    Return numeric feature columns, excluding Target, Date, Watershed_Loc, PixelRow, PixelCol, and extra_exclude.
    """
    exclude = set(extra_exclude or [])
    exclude |= {"Target", "Date", "Watershed_Loc", "PixelRow", "PixelCol"}

    feats = []
    for field in dataset.schema:
        if field.name in exclude:
            continue
        t = field.type
        if pa.types.is_floating(t) or pa.types.is_integer(t):
            feats.append(field.name)
    return feats
def compute_metrics(y_true: np.ndarray, y_proba: np.ndarray, threshold: float) -> Dict[str, float]:
    """Compute comprehensive metrics including confusion matrix components."""
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
    except:
        roc_auc = np.nan
    try:
        pr_auc = average_precision_score(y_true, y_proba)
    except:
        pr_auc = np.nan
    
    # Get confusion matrix components
    if len(np.unique(y_true)) > 1:
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    else:
        # Handle single-class case
        if y_true[0] == 0:
            tn = len(y_true)
            fp = fn = tp = 0
        else:
            tp = len(y_true)
            tn = fp = fn = 0
    
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


def compute_non_test_mean_std(
    dataset, feature_cols: List[str], non_test_ws: List[str], batch_size: int
) -> Tuple[np.ndarray, np.ndarray]:
    k = len(feature_cols)
    s = np.zeros(k, dtype=np.float64)
    ss = np.zeros(k, dtype=np.float64)
    cnt = np.zeros(k, dtype=np.float64)

    filt = ds.field("Watershed_Loc").isin(pa.array(non_test_ws))
    scanner = _make_scanner(dataset, columns=feature_cols, filter=filt, batch_size=batch_size)

    for batch in _scanner_batches(scanner):
        if batch.num_rows == 0:
            continue
        X = pa.Table.from_batches([batch]).to_pandas().to_numpy(dtype=np.float64)
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


def load_Xy_for_watersheds(
    dataset,
    watersheds: List[str],
    feature_cols: List[str],
    means: np.ndarray,
    stds: np.ndarray,
    batch_size: int,
    neg_pos_ratio: Optional[float] = None,
    random_state: int = 42
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load features/labels for watersheds.
    - Impute NaN with means and standardize (x - mean) / std
    - Optionally downsample negatives (keep <= neg_pos_ratio * positives) per watershed
    """
    rng = np.random.RandomState(random_state)
    X_parts, y_parts = [], []

    for ws in watersheds:
        filt = ds.field("Watershed_Loc") == ws
        cols = feature_cols + ["Target"]
        scanner = _make_scanner(dataset, columns=cols, filter=filt, batch_size=batch_size)

        X_ws, y_ws = [], []
        for batch in _scanner_batches(scanner):
            if batch.num_rows == 0:
                continue
            pdf = pa.Table.from_batches([batch]).to_pandas()
            yb = pdf["Target"].to_numpy(dtype=np.uint8, copy=False)
            Xb = pdf[feature_cols].to_numpy(dtype=np.float32, copy=False)

            # Impute NaN then standardize
            Xb = np.where(np.isnan(Xb), means, Xb)
            Xb = (Xb - means) / stds

            X_ws.append(Xb)
            y_ws.append(yb)

        if X_ws:
            Xw = np.concatenate(X_ws, axis=0)
            yw = np.concatenate(y_ws, axis=0)

            if neg_pos_ratio is not None:
                pos_idx = np.where(yw == 1)[0]
                neg_idx = np.where(yw == 0)[0]
                if len(pos_idx) > 0 and len(neg_idx) > 0:
                    keep_neg = min(len(neg_idx), int(len(pos_idx) * neg_pos_ratio))
                    if keep_neg < len(neg_idx):
                        # Systematic selection for diversity without relying on PixelRow/Col
                        step = max(len(neg_idx) // keep_neg, 1)
                        neg_idx = neg_idx[::step][:keep_neg]
                    keep_idx = np.concatenate([pos_idx, neg_idx])
                    rng.shuffle(keep_idx)
                    Xw = Xw[keep_idx]
                    yw = yw[keep_idx]

            X_parts.append(Xw)
            y_parts.append(yw)

    if not X_parts:
        return np.empty((0, len(feature_cols)), dtype=np.float32), np.empty((0,), dtype=np.uint8)

    X = np.concatenate(X_parts, axis=0).astype(np.float32, copy=False)
    y = np.concatenate(y_parts, axis=0).astype(np.uint8, copy=False)
    return X, y


# --------------------------- Metrics utilities ---------------------------

def find_optimal_threshold(y_true: np.ndarray, y_proba: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return 0.5
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_proba)
    f1 = 2 * (precisions * recalls) / (precisions + recalls + 1e-10)
    if f1.size == 0 or np.all(np.isnan(f1)):
        return 0.5
    best_idx = int(np.nanargmax(f1))
    return float(thresholds[best_idx]) if best_idx < len(thresholds) else 0.5



# ---------------------------- Main training ----------------------------

def train_random_forest_loso(
    parquet_uri: str,
    config: TrainingConfig = None,
    save_path: Optional[str] = None
) -> Dict:
    """
    End-to-end training (strict LOSO CV only) with proper per-fold normalization.
    """
    t0 = time.time()
    cfg = config or TrainingConfig()

    dataset = _open_dataset(parquet_uri)
    watersheds = get_unique_watersheds(dataset)
    if len(watersheds) < 3:
        raise RuntimeError(
            f"Insufficient watersheds: need >= 3 for LOSO (1 test + >=2 non-test). Found {len(watersheds)}."
        )

    test_ws = watersheds[0]
    non_test_ws = watersheds[1:]

    print("=" * 72)
    print("RANDOM FOREST FLOOD MODEL — LOSO (STRICT)")
    print("=" * 72)
    print(f"Watersheds: {watersheds}")
    print(f"Test watershed: {test_ws}")
    print(f"Non-test for CV: {non_test_ws}")

    feature_cols = get_feature_columns(dataset)
    print(f"\nUsing {len(feature_cols)} features")
    print("First 10:", feature_cols[:10])

    # LOSO CV across non-test watersheds
    print("\nRunning LOSO CV on non-test watersheds...")
    cv_thresholds = []
    cv_rows = []
    
    # Track all metrics for summary
    train_metrics_list = []
    val_metrics_list = []

    for i, val_ws in enumerate(non_test_ws):
        train_ws = [w for w in non_test_ws if w != val_ws]
        print(f"  Fold {i+1}/{len(non_test_ws)}: val={val_ws} | train={len(train_ws)} watersheds")

        # FIX: Compute normalization PER FOLD on train watersheds only (no leakage!)
        fold_means, fold_stds = compute_non_test_mean_std(
            dataset, feature_cols, train_ws, batch_size=cfg.batch_size
        )

        # Train set
        X_train, y_train = load_Xy_for_watersheds(
            dataset, train_ws, feature_cols, fold_means, fold_stds,
            batch_size=cfg.batch_size, neg_pos_ratio=cfg.neg_pos_ratio,
            random_state=cfg.random_state + i
        )
        if len(y_train) == 0:
            raise RuntimeError(f"Empty training set for fold with val={val_ws}.")

        cw = compute_class_weight('balanced', classes=np.array([0, 1]), y=y_train)
        rf = RandomForestClassifier(
            n_estimators=cfg.n_estimators,
            max_depth=cfg.max_depth,
            min_samples_split=cfg.min_samples_split,
            min_samples_leaf=cfg.min_samples_leaf,
            max_features=cfg.max_features,
            criterion=cfg.criterion,
            min_impurity_decrease=cfg.min_impurity_decrease,
            max_leaf_nodes=cfg.max_leaf_nodes,
            class_weight={0: cw[0], 1: cw[1]},
            bootstrap=True,
            oob_score=True,
            n_jobs=cfg.n_jobs,
            random_state=cfg.random_state + i
        )
        rf.fit(X_train, y_train)

        # Get train predictions for metrics
        proba_train = rf.predict_proba(X_train)[:, 1]

        # Validation set (no downsampling) - use same normalization as training
        X_val, y_val = load_Xy_for_watersheds(
            dataset, [val_ws], feature_cols, fold_means, fold_stds,
            batch_size=cfg.batch_size, neg_pos_ratio=None
        )
        if len(y_val) == 0:
            raise RuntimeError(f"Empty validation set for fold with val={val_ws}.")

        proba_val = rf.predict_proba(X_val)[:, 1]
        t_star = find_optimal_threshold(y_val, proba_val)
        cv_thresholds.append(t_star)
        
        # Use compute_metrics instead of evaluate_metrics
        train_metrics = compute_metrics(y_train, proba_train, t_star)
        val_metrics = compute_metrics(y_val, proba_val, t_star)
        
        train_metrics["fold"] = i + 1
        train_metrics["set"] = "train"
        val_metrics["fold"] = i + 1
        val_metrics["set"] = "val"
        
        train_metrics_list.append(train_metrics)
        val_metrics_list.append(val_metrics)
        
        cv_rows.append(val_metrics)

    cv_df = pd.DataFrame(cv_rows)
    all_metrics_df = pd.DataFrame(train_metrics_list + val_metrics_list)
    
    print("\nCV summary (validation):")
    print(cv_df[["f1", "precision", "recall", "roc_auc", "pr_auc", "threshold"]].round(4))
    
    # Compute mean metrics for summary
    metrics_summary = {}
    for set_name in ["train", "val"]:
        set_df = all_metrics_df[all_metrics_df["set"] == set_name]
        if not set_df.empty:
            metrics_summary[set_name] = set_df.drop(columns=["fold", "set"]).mean().to_dict()

    final_threshold = float(np.mean(cv_thresholds))
    print(f"\nChosen threshold from LOSO CV: {final_threshold:.3f}")

    # Train final model on all non-test watersheds
    print("\nTraining final model on all non-test watersheds...")
    
    # Final normalization on all non-test watersheds
    final_means, final_stds = compute_non_test_mean_std(
        dataset, feature_cols, non_test_ws, batch_size=cfg.batch_size
    )
    
    X_train_all, y_train_all = load_Xy_for_watersheds(
        dataset, non_test_ws, feature_cols, final_means, final_stds,
        batch_size=cfg.batch_size, neg_pos_ratio=cfg.neg_pos_ratio,
        random_state=cfg.random_state
    )
    cw = compute_class_weight('balanced', classes=np.array([0, 1]), y=y_train_all)
    rf_final = RandomForestClassifier(
        n_estimators=cfg.n_estimators,
        max_depth=cfg.max_depth,
        min_samples_split=cfg.min_samples_split,
        min_samples_leaf=cfg.min_samples_leaf,
        max_features=cfg.max_features,
        criterion=cfg.criterion,
        min_impurity_decrease=cfg.min_impurity_decrease,
        max_leaf_nodes=cfg.max_leaf_nodes,
        class_weight={0: cw[0], 1: cw[1]},
        bootstrap=True,
        oob_score=True,
        n_jobs=cfg.n_jobs,
        random_state=cfg.random_state
    )

    rf_final.fit(X_train_all, y_train_all)
    print(f"OOB score: {getattr(rf_final, 'oob_score_', np.nan):.4f}")

    # Evaluate on held-out test watershed
    print("\nEvaluating on held-out test watershed...")
    X_test, y_test = load_Xy_for_watersheds(
        dataset, [test_ws], feature_cols, final_means, final_stds,
        batch_size=cfg.batch_size, neg_pos_ratio=None
    )
    
    if len(y_test) > 0:
        proba_test = rf_final.predict_proba(X_test)[:, 1]
        test_metrics = compute_metrics(y_test, proba_test, final_threshold)
        metrics_summary["test"] = test_metrics
    else:
        test_metrics = {}
        
    print("\nTest metrics:")
    for m in ["accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc", "kappa", "mcc", "tpr", "tnr"]:
        print(f"  {m:10s}: {test_metrics.get(m, np.nan):.4f}")

    # Feature importances
    importances = rf_final.feature_importances_
    importance_df = pd.DataFrame({"feature": feature_cols, "importance": importances}) \
                      .sort_values("importance", ascending=False)
    print("\nTop 10 features:")
    print(importance_df.head(10).to_string(index=False))

    # Optional save bundle
    if save_path:
        with open(save_path, "wb") as f:
            pickle.dump({
                "model": rf_final,
                "feature_cols": feature_cols,
                "means": final_means,  # Use final normalization stats
                "stds": final_stds,
                "threshold": final_threshold,
                "config": cfg,
                "metrics_summary": metrics_summary,
                "test_watersheds": [test_ws],
                "feature_importance": importance_df
            }, f)
        print(f"\nSaved model bundle to: {save_path}")

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed/60:.1f} minutes")

    return {
        "model": rf_final,
        "feature_cols": feature_cols,
        "means": final_means,  # Use final normalization stats
        "stds": final_stds,
        "threshold": final_threshold,
        "cv_results": cv_df,
        "test_metrics": test_metrics,
        "metrics_summary": metrics_summary,
        "feature_importance": importance_df,
        "watersheds": watersheds,
        "test_watersheds": [test_ws],
        "config": cfg
    }




"""
Visualization utilities for Random Forest flood prediction model outputs.
Creates:
  1) predictions TIFF (0=land, 1=water, 255=nodata),
  2) confusion TIFF (0=TN, 1=TP, 2=FN, 3=FP, 255=nodata),
  3) per-date PNGs with Date, F1, Accuracy, Precision, Recall, ROC AUC.

Notes:
- Accept latitude/longitude as strings; we normalize to 4 decimals.
- Does NOT require 'test_watershed' in the model bundle. We parse Lat/Lon from the base_s3_path.
- Computes AUC from probabilities (not hard labels).
- Preloads static/time-series once; vectorized confusion construction.
"""


# -------------- small helpers --------------

def _normalize_coord_str(value: str) -> str:
    """
    Accept a string or number-like value and return a 4-decimal string.
    Example: "-99.3920"
    """
    s = str(value).strip()
    try:
        return f"{float(s):.4f}"
    except Exception:
        # If not a number, just return as-is (best effort)
        return s

def _extract_lat_lon_from_training_path(base_s3_path: str) -> Tuple[str, str]:
    """
    Parse 'trainingFolderFor_Lat_{lat}_Lon_{lon}' from the given base path.
    Returns (lat_str, lon_str).
    """
    # Expect .../trainingFolderFor_Lat_30.0750_Lon_-99.3920
    m = re.search(r"trainingFolderFor_Lat_([-\d\.]+)_Lon_([-\d\.]+)$", base_s3_path)
    if not m:
        # try without end anchor, in case trailing slash is omitted
        m = re.search(r"trainingFolderFor_Lat_([-\d\.]+)_Lon_([-\d\.]+)", base_s3_path)
    if not m:
        raise ValueError(f"Could not parse Lat/Lon from base_s3_path: {base_s3_path}")
    return m.group(1), m.group(2)

def _configure_gdal_proj():
    """
    Ensure GDAL/PROJ can find its data directories. Safe to call multiple times.
    """
    # Try pyproj packaged data dir (most reliable under conda)
    try:
        from pyproj import datadir as _pydatadir
        proj_dir = _pydatadir.get_data_dir()
        if proj_dir and os.path.isdir(proj_dir):
            os.environ["PROJ_LIB"] = proj_dir
            gdal.SetConfigOption("PROJ_LIB", proj_dir)
    except Exception:
        pass

    # Fallback candidates (common paths)
    if "PROJ_LIB" not in os.environ or not os.path.isdir(os.environ["PROJ_LIB"]):
        for cand in ["/opt/conda/share/proj", "/usr/share/proj"]:
            if os.path.isdir(cand):
                os.environ["PROJ_LIB"] = cand
                gdal.SetConfigOption("PROJ_LIB", cand)
                break

    # GDAL data dir (for EPSG, coordinate system definitions)
    try:
        gdal_mod_dir = os.path.dirname(gdal.__file__)
        gdal_data_dir = os.path.join(gdal_mod_dir, "data")
        if os.path.isdir(gdal_data_dir):
            os.environ["GDAL_DATA"] = gdal_data_dir
            gdal.SetConfigOption("GDAL_DATA", gdal_data_dir)
    except Exception:
        pass

    # Optional: allow network access for PROJ if local grids/db unavailable
    gdal.SetConfigOption("PROJ_NETWORK", "ON")

class FloodModelVisualizer:
    """Generate visual outputs for flood prediction model results."""

    def __init__(self, model_results: Dict, parquet_uri: str, base_s3_path: str):
        """
        Parameters:
        - model_results: dict containing model, feature_cols, means, stds, threshold (from training)
        - parquet_uri: S3 URI to parquet (not used here, kept for backward compatibility)
        - base_s3_path: root training folder, e.g.,
          s3://.../trainingFolderFor_Lat_{lat}_Lon_{lon}
        """
        self.model = model_results["model"]
        self.feature_cols = model_results["feature_cols"]
        self.means = model_results["means"]
        self.stds = model_results["stds"]
        self.threshold = model_results["threshold"]

        self.parquet_uri = parquet_uri
        self.base_s3_path = base_s3_path.rstrip('/')

        # Derive watershed coordinates from base_s3_path; do not rely on model bundle
        self.test_lat, self.test_lon = _extract_lat_lon_from_training_path(self.base_s3_path)
        self.test_watershed = f"Lat_{self.test_lat}_Lon_{self.test_lon}"

        # S3 client
        self.s3_client = boto3.client('s3')

        # GDAL S3 config
        gdal.SetConfigOption('AWS_NO_SIGN_REQUEST', 'NO')
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS', 'tif,tiff')

        # Caches
        self._static_map: Optional[Dict[str, np.ndarray]] = None   # name -> flattened array
        self._ts_by_date: Optional[Dict[str, Dict[str, float]]] = None  # date_str -> dict(feature_name->value)

    @staticmethod
    def _sanitize(name: str) -> str:
        """Sanitize band descriptions to match training feature names."""
        name = (name or "").strip()
        name = re.sub(r"[^\w\-]+", "_", name)
        name = re.sub(r"_+", "_", name)
        return name

    def generate_all_outputs(self, output_s3_folder: Optional[str] = None) -> Dict[str, Any]:
        """
        Generate predictions TIFF, confusion-matrix TIFF, and per-band PNGs.
        Returns s3 paths dict.
        """
        if output_s3_folder is None:
            output_s3_folder = f"{self.base_s3_path}/assessModel"

        print("\n" + "="*60)
        print("GENERATING MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        print(f"Test watershed: {self.test_watershed}")
        print(f"Output folder: {output_s3_folder}")

        # Flood time-series TIFF (multi-band, one band per date, nodata=255)
        flood_tiff_path = (
            f"{self.base_s3_path}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        print(f"\nLoading flood maps from: {flood_tiff_path}")
        flood_ds, flood_data, dates = self._load_flood_tiff(flood_tiff_path)

        # Preload static features and time series for this watershed
        print("Preloading static features from rasters/ and time series CSVs...")
        self._preload_static_map()
        self._preload_timeseries()

        # Generate predictions/confusion + per-band metrics
        print("\nGenerating predictions and confusion matrices...")
        predictions, confusions, per_date_metrics = self._generate_predictions(flood_data, dates)

        # Write TIFFs
        print("\nCreating output TIFFs...")
        pred_tiff_path = f"{output_s3_folder}/predictions_{self.test_watershed}.tif"
        cm_tiff_path = f"{output_s3_folder}/confusion_matrix_{self.test_watershed}.tif"
        self._write_multiband_tiff(flood_ds, predictions, dates, pred_tiff_path, prefix="Prediction")
        self._write_multiband_tiff(flood_ds, confusions, dates, cm_tiff_path, prefix="CM_0TN_1TP_2FN_3FP")

        # PNG overlays
        print("\nCreating PNG visualizations...")
        png_paths = self._create_png_visualizations(confusions, dates, flood_data, predictions, per_date_metrics, output_s3_folder)

        print("\n✅ Visualization complete!")
        return {
            "prediction_tiff": pred_tiff_path,
            "confusion_matrix_tiff": cm_tiff_path,
            "png_visualizations": png_paths
        }

    def _load_flood_tiff(self, s3_path: str) -> Tuple[gdal.Dataset, Dict[str, np.ndarray], List[str]]:
        vsi = s3_path.replace('s3://', '/vsis3/')
        ds = gdal.Open(vsi, gdal.GA_ReadOnly)
        if ds is None:
            raise ValueError(f"Could not open flood TIFF: {s3_path}")

        n_bands = ds.RasterCount
        n_rows = ds.RasterYSize
        n_cols = ds.RasterXSize
        print(f"  Dimensions: {n_cols}x{n_rows}, {n_bands} bands")

        flood_data: Dict[str, np.ndarray] = {}
        dates: List[str] = []
        for b in range(1, n_bands + 1):
            band = ds.GetRasterBand(b)
            desc = band.GetDescription() or ""
            # Expect "DSWx_S1_YYYY-MM-DD"
            m = re.search(r"\d{4}-\d{2}-\d{2}", desc)
            date_str = m.group(0) if m else f"band_{b}"
            dates.append(date_str)
            flood_data[date_str] = band.ReadAsArray()

        return ds, flood_data, dates

    def _preload_static_map(self):
        """Load all static bands into a dict name->flattened array with NoData mapped to NaN."""
        if self._static_map is not None:
            return
        static_tiff_path = (
            f"{self.base_s3_path}/rasters/"
            f"multiband_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        vsi = static_tiff_path.replace('s3://', '/vsis3/')
        ds = gdal.Open(vsi, gdal.GA_ReadOnly)
        static_map = {}
        if ds is None or ds.RasterCount == 0:
            print(f"  Warning: Could not load static raster: {static_tiff_path}")
            self._static_map = {}
            return

        for b in range(1, ds.RasterCount + 1):
            band = ds.GetRasterBand(b)
            nm = self._sanitize(band.GetDescription() or f"static_band_{b:02d}")
            arr = band.ReadAsArray().astype(np.float32)
            nd = band.GetNoDataValue()
            mask = np.zeros(arr.shape, dtype=bool)
            if nd is not None:
                mask |= (arr == float(nd))
            mask |= (arr == -9999.0)     # alignment sentinel
            mask |= (arr == -32768.0)    # SoilGrids NoData
            arr = np.where(mask, np.nan, arr)
            static_map[nm] = arr.ravel()
        ds = None
        self._static_map = static_map

    def _preload_timeseries(self):
        """Load both time series CSVs once into a dict date->feature_name->value."""
        if self._ts_by_date is not None:
            return
        ts_by_date: Dict[str, Dict[str, float]] = {}

        def add_rows(df: pd.DataFrame, prefix: str):
            if 'Date' not in df.columns:
                return
            df = df.copy()
            df['Date'] = pd.to_datetime(df['Date']).dt.strftime('%Y-%m-%d')
            for _, row in df.iterrows():
                date = row['Date']
                d = ts_by_date.setdefault(date, {})
                for col in ('streamflow_pred', 't2m', 'ssrd', 'tp'):
                    if col in row:
                        d[f"{prefix}_{col}"] = float(row[col]) if pd.notna(row[col]) else np.nan

        # subwatershed
        sub_csv = (
            f"{self.base_s3_path}/streamflowPredictions/"
            f"subwatershed_predictions_Lat{self.test_lat}_Lon{self.test_lon}.csv"
        )
        # river basin
        riv_csv = (
            f"{self.base_s3_path}/streamflowPredictions/"
            f"river_basin_predictions_Lat{self.test_lat}_Lon{self.test_lon}.csv"
        )

        for path, prefix in [(sub_csv, "subwatershed"), (riv_csv, "river_basin")]:
            try:
                df = self._read_csv_from_s3(path)
                add_rows(df, prefix)
            except Exception:
                pass

        self._ts_by_date = ts_by_date

    def _read_csv_from_s3(self, s3_path: str) -> pd.DataFrame:
        bucket, key = s3_path.replace('s3://', '').split('/', 1)
        obj = self.s3_client.get_object(Bucket=bucket, Key=key)
        return pd.read_csv(io.BytesIO(obj['Body'].read()))

    def _generate_predictions(
        self,
        flood_data: Dict[str, np.ndarray],
        dates: List[str]
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Dict[str, Dict[str, float]]]:
        """Generate predictions/confusions for all dates; return per-date metrics using probabilities."""
        predictions: Dict[str, np.ndarray] = {}
        confusions: Dict[str, np.ndarray] = {}
        per_date_metrics: Dict[str, Dict[str, float]] = {}

        first = next(iter(flood_data.values()))
        n_rows, n_cols = first.shape
        n_pix = n_rows * n_cols
        print(f"  Processing {len(dates)} dates, {n_pix:,} pixels per date")

        # Identify which feature_cols are time-series vs static
        ts_cols = [c for c in self.feature_cols if c.startswith("subwatershed_") or c.startswith("river_basin_")]
        static_cols = [c for c in self.feature_cols if c not in ts_cols]

        for idx, date_str in enumerate(dates):
            if (idx + 1) % 10 == 0:
                print(f"  Date {idx+1}/{len(dates)}: {date_str}")

            arr = flood_data[date_str]
            valid_mask = (arr != 255)  # flood NoData
            valid_idx = np.where(valid_mask.ravel())[0]
            if valid_idx.size == 0:
                # No valid pixels for this date
                predictions[date_str] = np.full((n_rows, n_cols), 255, dtype=np.uint8)
                confusions[date_str] = np.full((n_rows, n_cols), 255, dtype=np.uint8)
                per_date_metrics[date_str] = {"F1": np.nan, "Accuracy": np.nan, "Precision": np.nan, "Recall": np.nan, "ROC_AUC": np.nan}
                continue

            # Build X for valid pixels efficiently:
            # - static: pull from self._static_map[name][valid_idx]
            # - timeseries: broadcast scalar from self._ts_by_date[date_str]
            X = np.empty((valid_idx.size, len(self.feature_cols)), dtype=np.float32)

            # Fill static columns
            for j, name in enumerate(self.feature_cols):
                if name in static_cols:
                    if self._static_map and (name in self._static_map):
                        X[:, j] = self._static_map[name][valid_idx]
                    else:
                        X[:, j] = np.nan

            # Fill time-series columns (constant per date across pixels)
            ts_vals = self._ts_by_date.get(date_str, {}) if self._ts_by_date else {}
            for j, name in enumerate(self.feature_cols):
                if name in ts_cols:
                    val = ts_vals.get(name, np.nan)
                    X[:, j] = val

            # Impute -> standardize with training stats
            X = np.where(np.isnan(X), self.means, X)
            X = (X - self.means) / self.stds

            # Predict
            y_proba = self.model.predict_proba(X)[:, 1]
            y_pred = (y_proba >= self.threshold).astype(np.uint8)

            # Ground truth (0 land, 1 water)
            y_true = (arr.ravel()[valid_idx] > 0).astype(np.uint8)

            # Predictions array
            pred_full = np.full(n_pix, 255, dtype=np.uint8)
            pred_full[valid_idx] = y_pred
            predictions[date_str] = pred_full.reshape(n_rows, n_cols)

            # Confusion array (vectorized): 0=TN, 1=TP, 2=FN, 3=FP
            cm_full = np.full(n_pix, 255, dtype=np.uint8)
            cm_valid = np.where(
                y_true == 0,
                np.where(y_pred == 0, 0, 3),  # TN or FP
                np.where(y_pred == 1, 1, 2)   # TP or FN
            ).astype(np.uint8)
            cm_full[valid_idx] = cm_valid
            confusions[date_str] = cm_full.reshape(n_rows, n_cols)

            # Per-date metrics (using probabilities for ROC AUC)
            try:
                auc = roc_auc_score(y_true, y_proba)
            except Exception:
                auc = np.nan
            metrics = {
                "F1": f1_score(y_true, y_pred, zero_division=0),
                "Accuracy": accuracy_score(y_true, y_pred),
                "Precision": precision_score(y_true, y_pred, zero_division=0),
                "Recall": recall_score(y_true, y_pred, zero_division=0),
                "ROC_AUC": float(auc)
            }
            per_date_metrics[date_str] = metrics

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
            # FIX: Explicitly set dtype=np.uint8
            band.WriteArray(arrays.get(date_str, np.full((ref_ds.RasterYSize, ref_ds.RasterXSize), 255, dtype=np.uint8)))
        
        out = None
        
        # Upload to S3
        bucket, key = output_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(temp_path, bucket, key)
        os.remove(temp_path)
        print(f"  Created: {output_path}")

    def _create_png_visualizations(
        self,
        confusion_matrices: Dict[str, np.ndarray],
        dates: List[str],
        flood_data: Dict[str, np.ndarray],
        predictions: Dict[str, np.ndarray],
        per_date_metrics: Dict[str, Dict[str, float]],
        output_s3_folder: str
    ) -> List[str]:
        png_paths = []
        for i, date_str in enumerate(dates):
            cm = confusion_matrices.get(date_str)
            if cm is None:
                continue

            # Require enough valid pixels
            if np.count_nonzero(cm != 255) < 100:
                continue

            metrics = per_date_metrics.get(date_str, {})
            auc_val = metrics.get("ROC_AUC", np.nan)
            auc_str = f"{auc_val:.3f}" if np.isfinite(auc_val) else "N/A"

            fig, ax = plt.subplots(figsize=(12, 10))

            cm_disp = cm.astype(float)
            cm_disp[cm == 255] = np.nan
            cmap = plt.matplotlib.colors.ListedColormap(['green', 'blue', 'orange', 'red'])
            ax.imshow(cm_disp, cmap=cmap, vmin=0, vmax=3)

            legend_elements = [
                mpatches.Patch(color='green', label='True Negative (Land correct)'),
                mpatches.Patch(color='blue', label='True Positive (Water correct)'),
                mpatches.Patch(color='orange', label='False Negative (Missed water)'),
                mpatches.Patch(color='red', label='False Positive (False water)')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10)
            ax.set_title(f'Confusion Matrix - {date_str}', fontsize=14, fontweight='bold')

            text = (
                f"Date: {date_str}\n"
                f"F1 Score: {metrics.get('F1', np.nan):.3f}\n"
                f"Accuracy: {metrics.get('Accuracy', np.nan):.3f}\n"
                f"Precision: {metrics.get('Precision', np.nan):.3f}\n"
                f"Recall: {metrics.get('Recall', np.nan):.3f}\n"
                f"ROC AUC: {auc_str}"
            )
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
            ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11, va='top', bbox=props)
            ax.axis('off')
            plt.tight_layout()

            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()

            s3_path = f"{output_s3_folder}/confusion_matrix_{date_str}.png"
            self._upload_to_s3(temp_png, s3_path)
            os.remove(temp_png)
            png_paths.append(s3_path)

            if (i + 1) % 10 == 0:
                print(f"  Created {i + 1} PNGs")

        print(f"  Created {len(png_paths)} PNG visualizations")
        return png_paths

    def _upload_to_s3(self, local_path: str, s3_path: str):
        bucket, key = s3_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(local_path, bucket, key)


class RandomForestFloodVisualizer:
    """Generate visualizations matching other model formats."""
    
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
        pred_tiff = f"{output_folder}/randomforest_predictions_{self.test_watershed}.tif"
        cm_tiff = f"{output_folder}/randomforest_confusion_matrix_{self.test_watershed}.tif"
        
        self._write_multiband_tiff(ref_ds, predictions, dates, pred_tiff, "RandomForest_Prediction")
        self._write_multiband_tiff(ref_ds, confusions, dates, cm_tiff, "RandomForest_CM")
        
        # Create PNGs
        png_paths = self._create_pngs(confusions, dates, per_date_metrics, output_folder)
        
        # Create CSVs
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
        
        print(f"  Processing {len(dates)} dates...")
        
        for date_idx, date_str in enumerate(dates):
            if (date_idx + 1) % 10 == 0:
                print(f"    Date {date_idx + 1}/{len(dates)}")
            
            date_data = df_all[df_all['Date'] == date_str]
            
            # Extract features
            X = date_data[self.feature_cols].to_numpy(dtype=np.float32)
            y_true = date_data['Target'].to_numpy(dtype=np.uint8)
            rows = date_data['PixelRow'].to_numpy()
            cols = date_data['PixelCol'].to_numpy()
            
            # Normalize
            X = np.where(np.isnan(X), self.means, X)
            X = (X - self.means) / self.stds
            
            # Predict with Random Forest
            y_proba = self.model.predict_proba(X)[:, 1]
            y_pred = (y_proba >= self.threshold).astype(np.uint8)
            
            # Create images
            pred_img = np.full((H, W), 255, dtype=np.uint8)
            pred_img[rows, cols] = y_pred
            predictions[date_str] = pred_img
            
            # Confusion matrix (0=TN, 1=TP, 2=FN, 3=FP, 255=nodata)
            cm_img = np.full((H, W), 255, dtype=np.uint8)
            cm_vals = np.where(
                y_true == 0,
                np.where(y_pred == 0, 0, 3),  # TN=0, FP=3
                np.where(y_pred == 1, 1, 2)   # TP=1, FN=2
            )
            cm_img[rows, cols] = cm_vals
            confusions[date_str] = cm_img
            
            # Metrics
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
                'F1': f1_score(y_true, y_pred, zero_division=0),
                'Accuracy': accuracy_score(y_true, y_pred),
                'Precision': precision_score(y_true, y_pred, zero_division=0),
                'Recall': recall_score(y_true, y_pred, zero_division=0),
                'ROC_AUC': auc,
                'PR_AUC': pr_auc,
                'Kappa': cohen_kappa_score(y_true, y_pred),
                'MCC': matthews_corrcoef(y_true, y_pred) if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1 else 0.0,
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
        
        # Upload to S3
        bucket, key = output_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(temp_path, bucket, key)
        os.remove(temp_path)
        print(f"  Created: {output_path}")
    
    def _create_pngs(self, confusions, dates, metrics, output_folder):
        """Create PNG visualizations."""
        png_paths = []
        
        # Process all dates or up to max_png_dates
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
            
            s3_path = f"{output_folder}/randomforest_confusion_matrix_{date_str}.png"
            bucket, key = s3_path.replace('s3://', '').split('/', 1)
            self.s3_client.upload_file(temp_png, bucket, key)
            os.remove(temp_png)
            png_paths.append(s3_path)
        
        print(f"  Created {len(png_paths)} PNG visualizations")
        return png_paths
    
    def _write_csvs(self, output_folder, per_date_metrics):
        """Write CSV outputs with GLM-style parity (no confusion matrix components in summary)."""
        paths = {}
        
        # Per-date metrics
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
            print(f"  Created: {s3_path}")
        
        # Summary metrics - GLM parity (no tp/fp/tn/fn/tpr/tnr)
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
                    # Removed tp, fp, tn, fn, tpr, tnr for GLM parity
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
                print(f"  Created: {s3_path}")
        
        # Feature importance
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
            print(f"  Created: {s3_path}")
        
        return paths


def visualize_randomforest_outputs(model_results: Dict, parquet_uri: str, 
                                  training_folder: str) -> Dict[str, Any]:
    """
    Generate all visualization outputs for Random Forest model.
    """
    # Ensure we have the RandomForestFloodVisualizer class
    viz = RandomForestFloodVisualizer(model_results, parquet_uri, training_folder)
    return viz.generate_all_outputs()
    
    
    
    
    
# temporary: public wrappers
# ======================== Public API Wrappers ========================

def open_parquet_dataset(parquet_uri: str):
    """Public wrapper for _open_dataset."""
    return _open_dataset(parquet_uri)

def make_scanner(dataset, columns=None, filter=None, batch_size=None):
    """Public wrapper for _make_scanner."""
    return _make_scanner(dataset, columns, filter, batch_size)

def scanner_batches(scanner):
    """Public wrapper for _scanner_batches."""
    return _scanner_batches(scanner)

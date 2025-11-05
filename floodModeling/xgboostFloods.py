"""
xgboostFloods.py

Memory-efficient XGBoost flood prediction with:
- Spatial K-Folds CV (random 60/20/20 per fold) across watersheds
- Streaming training via PyArrow batches with bounded memory
- Adaptive scale_pos_weight tuning per fold based on actual class distribution
- Imputation + standardization computed on TRAIN-only per fold
- Outputs and visualization compatible with existing pipeline
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

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

import xgboost as xgb
from xgboost.core import DataIter
from xgboost.callback import EarlyStopping as XGBEarlyStopping

from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, average_precision_score, confusion_matrix,
    cohen_kappa_score, matthews_corrcoef, precision_recall_curve
)


# Parquet I/O
import pyarrow as pa
import pyarrow.dataset as ds

# Geospatial viz deps
from osgeo import gdal
gdal.UseExceptions()
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import boto3


# ============================== CONFIG ==============================

class XGBoostConfig:
    def __init__(
        self,
        # XGBoost model params
        n_estimators: int = 300,
        max_depth: int = 6,
        learning_rate: float = 0.1,
        subsample: float = 0.8,
        colsample_bytree: float = 0.8,
        min_child_weight: int = 3,
        gamma: float = 0.1,
        reg_alpha: float = 0.01,
        reg_lambda: float = 1.0,
        objective: str = "binary:logistic",
        eval_metric: str = "logloss",
        tree_method: str = "hist",
        max_bin: int = 256,
        grow_policy: str = "depthwise",

        # Adaptive weight tuning
        scale_pos_weight_method: str = "sqrt",  # 'sqrt', 'log', 'capped', 'linear', or 'none'
        scale_pos_weight_cap: float = 50.0,     # Maximum allowed scale_pos_weight
        
        # General training
        random_state: int = 42,
        batch_size: int = 250_000,
        n_jobs: int = -1,
        early_stopping_rounds: Optional[int] = 20,
        verbosity: int = 0,

        # Memory + streaming controls
        use_external_memory: bool = True,
        max_train_rows: Optional[int] = 3_000_000,
        max_val_rows: Optional[int] = 2_000_000,
        max_test_rows: Optional[int] = 5_000_000,
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.learning_rate = learning_rate
        self.subsample = subsample
        self.colsample_bytree = colsample_bytree
        self.min_child_weight = min_child_weight
        self.gamma = gamma
        self.reg_alpha = reg_alpha
        self.reg_lambda = reg_lambda
        self.objective = objective
        self.eval_metric = eval_metric
        self.tree_method = tree_method
        self.max_bin = max_bin
        self.grow_policy = grow_policy

        self.scale_pos_weight_method = scale_pos_weight_method
        self.scale_pos_weight_cap = scale_pos_weight_cap

        self.random_state = random_state
        self.batch_size = batch_size
        self.n_jobs = n_jobs
        self.early_stopping_rounds = early_stopping_rounds
        self.verbosity = verbosity

        self.use_external_memory = use_external_memory
        self.max_train_rows = max_train_rows
        self.max_val_rows = max_val_rows
        self.max_test_rows = max_test_rows


# ======================= PICKLE-SAFE BOOSTER WRAPPER =======================

class BoosterPredictor:
    """
    Pickle-safe wrapper around xgboost.Booster with:
      - predict_proba(X) API (to match scikit estimators)
      - feature name validation during prediction
      - custom __getstate__/__setstate__ for robust pickling
    """
    def __init__(self, booster: xgb.Booster, best_iteration: Optional[int], feature_names: List[str], importance_map: Optional[Dict[str, float]] = None):
        self.booster = booster
        self.best_iteration = best_iteration
        self.feature_names = feature_names or []
        if importance_map is None:
            importance_map = booster.get_score(importance_type="gain")
        self.feature_importances_ = np.array([importance_map.get(f, 0.0) for f in self.feature_names], dtype=np.float32)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        dm = xgb.DMatrix(X, feature_names=self.feature_names)
        it_end = (self.best_iteration or 0) + 1 if self.best_iteration is not None else 0
        if it_end > 0:
            proba = self.booster.predict(dm, iteration_range=(0, it_end))
        else:
            proba = self.booster.predict(dm)
        return np.vstack([1.0 - proba, proba]).T

    def __getstate__(self):
        return {
            "booster_raw": self.booster.save_raw(),
            "best_iteration": self.best_iteration,
            "feature_names": self.feature_names,
            "feature_importances_": self.feature_importances_,
        }

    def __setstate__(self, state):
        booster = xgb.Booster()
        booster.load_model(bytearray(state["booster_raw"]))
        self.booster = booster
        self.best_iteration = state["best_iteration"]
        self.feature_names = state["feature_names"]
        self.feature_importances_ = state["feature_importances_"]


# ======================= PYARROW + S3 HELPERS =======================

def _get_s3fs():
    try:
        import pyarrow.fs as pafs
        return pafs.S3FileSystem()
    except Exception:
        pass
    try:
        from pyarrow import filesystem as pafs
        return pafs.S3FileSystem()
    except Exception:
        raise RuntimeError("Could not create a PyArrow S3FileSystem. Install/upgrade pyarrow or use local Parquet path.")


def _open_dataset(parquet_uri: str):
    if parquet_uri.startswith("s3://"):
        bkt, key = parquet_uri[5:].split("/", 1)
        s3fs = _get_s3fs()
        return ds.dataset(f"{bkt}/{key}", format="parquet", filesystem=s3fs)
    else:
        return ds.dataset(parquet_uri, format="parquet")


def _make_scanner(dataset, columns=None, filter=None, batch_size: Optional[int] = None):
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


# ============================ DATA UTILS ============================

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
    random_state: int = 42
) -> Tuple[np.ndarray, np.ndarray]:
    """Load and normalize data without any negative downsampling."""
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

            Xb = np.where(np.isnan(Xb), means, Xb)
            Xb = (Xb - means) / stds

            X_ws.append(Xb)
            y_ws.append(yb)

        if X_ws:
            Xw = np.concatenate(X_ws, axis=0)
            yw = np.concatenate(y_ws, axis=0)
            X_parts.append(Xw)
            y_parts.append(yw)

        del X_ws, y_ws
        gc.collect()

    if not X_parts:
        return np.empty((0, len(feature_cols)), dtype=np.float32), np.empty((0,), dtype=np.uint8)

    X = np.concatenate(X_parts, axis=0).astype(np.float32, copy=False)
    y = np.concatenate(y_parts, axis=0).astype(np.uint8, copy=False)
    X_parts.clear(); y_parts.clear(); gc.collect()
    return X, y


# ===================== STREAMING/EXTERNAL MEMORY =====================

def iter_Xy_batches(
    dataset,
    watersheds: List[str],
    feature_cols: List[str],
    means: np.ndarray,
    stds: np.ndarray,
    batch_size: int,
    random_state: int = 42
):
    """Stream data batches without any negative downsampling."""
    accum_X, accum_y = [], []
    rows = 0

    for ws in watersheds:
        filt = ds.field("Watershed_Loc") == ws
        cols = feature_cols + ["Target"]
        scanner = _make_scanner(dataset, columns=cols, filter=filt, batch_size=batch_size)

        for batch in _scanner_batches(scanner):
            if batch.num_rows == 0:
                continue
            pdf = pa.Table.from_batches([batch]).to_pandas()
            yb = pdf["Target"].to_numpy(dtype=np.uint8, copy=False)
            Xb = pdf[feature_cols].to_numpy(dtype=np.float32, copy=False)

            Xb = np.where(np.isnan(Xb), means, Xb)
            Xb = (Xb - means) / stds

            accum_X.append(Xb)
            accum_y.append(yb)
            rows += Xb.shape[0]

            if rows >= batch_size:
                X_batch = np.concatenate(accum_X, axis=0).astype(np.float32, copy=False)
                y_batch = np.concatenate(accum_y, axis=0).astype(np.uint8, copy=False)
                yield X_batch, y_batch
                accum_X.clear(); accum_y.clear(); rows = 0
                gc.collect()

    if rows > 0:
        X_batch = np.concatenate(accum_X, axis=0).astype(np.float32, copy=False)
        y_batch = np.concatenate(accum_y, axis=0).astype(np.uint8, copy=False)
        yield X_batch, y_batch
    accum_X.clear(); accum_y.clear()
    gc.collect()


class ArrowDataIter(DataIter):
    def __init__(self, batch_yield_fn, feature_names: List[str]):
        super().__init__()
        self._yield_fn = batch_yield_fn
        self._feature_names = feature_names
        self._gen = None

    def reset(self):
        self._gen = self._yield_fn()

    def next(self, input_data):
        if self._gen is None:
            self.reset()
        try:
            Xb, yb = next(self._gen)
            input_data(data=Xb, label=yb, feature_names=self._feature_names)
            return 0
        except StopIteration:
            return 1


def build_quantile_dmatrix_from_iter(
    dataset, watersheds, feature_cols, means, stds,
    batch_size, random_state, max_bin
) -> xgb.QuantileDMatrix:
    def _fn():
        return iter_Xy_batches(
            dataset=dataset,
            watersheds=watersheds,
            feature_cols=feature_cols,
            means=means, stds=stds,
            batch_size=batch_size,
            random_state=random_state
        )
    it = ArrowDataIter(_fn, feature_cols)
    qdm = xgb.QuantileDMatrix(it, feature_names=feature_cols, max_bin=max_bin)
    return qdm


def build_dmatrix_from_iter(
    dataset, watersheds, feature_cols, means, stds,
    batch_size, random_state, row_cap: Optional[int]
) -> xgb.DMatrix:
    X_parts, y_parts = [], []
    rows = 0
    for Xb, yb in iter_Xy_batches(dataset, watersheds, feature_cols, means, stds, batch_size, random_state=random_state):
        X_parts.append(Xb)
        y_parts.append(yb)
        rows += Xb.shape[0]
        if (row_cap is not None) and (rows >= row_cap):
            break
    if X_parts:
        X = np.concatenate(X_parts, axis=0)
        y = np.concatenate(y_parts, axis=0)
        dm = xgb.DMatrix(X, label=y, feature_names=feature_cols)
        del X_parts, y_parts, X, y; gc.collect()
        return dm
    return xgb.DMatrix(
        np.empty((0, len(feature_cols)), dtype=np.float32),
        label=np.empty((0,), dtype=np.uint8),
        feature_names=feature_cols
    )


# ============================= ADAPTIVE WEIGHTING =============================

def compute_adaptive_scale_pos_weight(
    y_train: np.ndarray,
    method: str = 'sqrt',
    cap: float = 50.0
) -> float:
    """
    Compute adaptive scale_pos_weight based on actual class distribution.
    
    Methods:
    - 'linear': natural ratio (can be extreme)
    - 'sqrt': square root dampening
    - 'log': logarithmic dampening
    - 'capped': linear with cap
    - 'none': no weighting (1.0)
    """
    n_pos = (y_train == 1).sum()
    n_neg = (y_train == 0).sum()
    
    if n_pos == 0:
        return 1.0
    
    natural_ratio = n_neg / n_pos
    
    if method == 'sqrt':
        # Square root dampening - recommended for extreme imbalance
        scale_pos_weight = np.sqrt(natural_ratio)
        
    elif method == 'log':
        # Logarithmic dampening
        scale_pos_weight = 1 + np.log(natural_ratio)
        
    elif method == 'capped' or method == 'linear':
        # Linear with optional cap
        scale_pos_weight = natural_ratio
        
    else:  # 'none' or unknown
        scale_pos_weight = 1.0
    
    # Apply cap
    scale_pos_weight = min(scale_pos_weight, cap)
    
    return float(scale_pos_weight)


# ============================= METRICS =============================

def find_optimal_threshold(y_true: np.ndarray, y_proba: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return 0.5
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_proba)
    f1 = 2 * (precisions * recalls) / (precisions + recalls + 1e-10)
    if f1.size == 0 or np.all(np.isnan(f1)):
        return 0.5
    best_idx = int(np.nanargmax(f1))
    return float(thresholds[best_idx]) if best_idx < len(thresholds) else 0.5


def evaluate_metrics(y_true: np.ndarray, y_proba: np.ndarray, threshold: float) -> Dict:
    y_pred = (y_proba >= threshold).astype(np.uint8)
    out = {
        "threshold": float(threshold),
        "prevalence": float(y_true.mean()),
        "accuracy": accuracy_score(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "kappa": cohen_kappa_score(y_true, y_pred),
        "mcc": matthews_corrcoef(y_true, y_pred),
    }
    try:
        out["roc_auc"] = roc_auc_score(y_true, y_proba)
    except Exception:
        out["roc_auc"] = np.nan
    try:
        out["pr_auc"] = average_precision_score(y_true, y_proba)
    except Exception:
        out["pr_auc"] = np.nan

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    out.update({
        "tp": int(tp), "fp": int(fp), "tn": int(tn), "fn": int(fn),
        "tpr": tp / (tp + fn) if (tp + fn) > 0 else 0.0,
        "tnr": tn / (tn + fp) if (tn + fp) > 0 else 0.0,
    })
    return out


# ============================ SPLIT LOGIC ============================

def _spatial_split_kfolds(
    watersheds: List[str],
    n_folds: int,
    random_state: int,
    train_frac: float = 0.6,
    val_frac: float = 0.2,
    force_test_ws: Optional[str] = None
) -> List[Tuple[List[str], List[str], List[str]]]:
    """
    Create n_folds random spatial splits (60/20/20) across watersheds.
    If force_test_ws is provided and present, it is included in the TEST set for every fold.
    """
    splits = []
    base_ws = watersheds[:]
    if force_test_ws and force_test_ws in base_ws:
        base_ws = [w for w in base_ws if w != force_test_ws]

    for f in range(n_folds):
        rng = np.random.RandomState(random_state + f)
        ws = base_ws[:]
        rng.shuffle(ws)
        n = len(ws)
        n_train = max(1, int(round(train_frac * n)))
        n_val = max(1, int(round(val_frac * n)))
        if n_train + n_val >= n:
            n_val = max(1, n - n_train - 1)
        train_ws = ws[:n_train]
        val_ws = ws[n_train:n_train+n_val]
        test_ws = ws[n_train+n_val:]
        if force_test_ws:
            test_ws = [force_test_ws] + test_ws
        if len(test_ws) == 0:
            test_ws = [val_ws.pop()] if len(val_ws) > 1 else [train_ws.pop()]
        splits.append((train_ws, val_ws, test_ws))
    return splits


# ============================== TRAINERS ==============================

def _train_one_split(
    dataset,
    feature_cols: List[str],
    cfg: XGBoostConfig,
    train_ws: List[str],
    val_ws: List[str],
    test_ws: List[str],
    fold_seed: int
) -> Tuple[Any, float, Dict, Dict, Dict, pd.DataFrame]:
    """
    Train a model for a single spatial split using memory-bounded streaming.
    Now computes adaptive scale_pos_weight based on actual class distribution.
    Returns: (model_obj, threshold, train_metrics, val_metrics, test_metrics, importance_df)
    """
    means, stds = compute_non_test_mean_std(dataset, feature_cols, train_ws, batch_size=cfg.batch_size)

    if cfg.use_external_memory:
        try:
            dtrain = build_quantile_dmatrix_from_iter(
                dataset, train_ws, feature_cols, means, stds,
                batch_size=cfg.batch_size,
                random_state=fold_seed, max_bin=cfg.max_bin
            )
            qdm_ok = True
        except ValueError:
            dtrain = build_dmatrix_from_iter(
                dataset, train_ws, feature_cols, means, stds,
                batch_size=cfg.batch_size, random_state=fold_seed,
                row_cap=cfg.max_train_rows
            )
            qdm_ok = False

        dval = build_dmatrix_from_iter(
            dataset, val_ws, feature_cols, means, stds,
            batch_size=cfg.batch_size, random_state=fold_seed, row_cap=cfg.max_val_rows
        )
        dtest = build_dmatrix_from_iter(
            dataset, test_ws, feature_cols, means, stds,
            batch_size=cfg.batch_size, random_state=fold_seed, row_cap=cfg.max_test_rows
        )

        ytr = dtrain.get_label()
        
        # Compute adaptive scale_pos_weight
        spw = compute_adaptive_scale_pos_weight(
            y_train=ytr,
            method=cfg.scale_pos_weight_method,
            cap=cfg.scale_pos_weight_cap
        )
        
        # Print class distribution info
        n_pos = int((ytr == 1).sum())
        n_neg = int((ytr == 0).sum())
        natural_ratio = n_neg / max(n_pos, 1)
        print(f"    Natural ratio: {natural_ratio:.1f}:1, Adaptive weight: {spw:.2f}")

        params = {
            "objective": cfg.objective,
            "eval_metric": cfg.eval_metric,
            "tree_method": cfg.tree_method or "hist",
            "max_depth": cfg.max_depth,
            "learning_rate": cfg.learning_rate,
            "subsample": cfg.subsample,
            "colsample_bytree": cfg.colsample_bytree,
            "min_child_weight": cfg.min_child_weight,
            "gamma": cfg.gamma,
            "reg_alpha": cfg.reg_alpha,
            "reg_lambda": cfg.reg_lambda,
            "scale_pos_weight": spw,
            "verbosity": cfg.verbosity,
            "grow_policy": cfg.grow_policy
        }
        if qdm_ok:
            params["max_bin"] = cfg.max_bin

        evals = [(dval, "validation")] if dval.num_row() > 0 else []
        booster = xgb.train(
            params=params,
            dtrain=dtrain,
            num_boost_round=cfg.n_estimators,
            evals=evals,
            early_stopping_rounds=cfg.early_stopping_rounds
        )

        if dval.num_row() > 0:
            proba_val = booster.predict(dval, iteration_range=(0, (getattr(booster, "best_iteration", 0) or 0) + 1))
            y_val = dval.get_label().astype(np.uint8, copy=False)
            threshold = find_optimal_threshold(y_val, proba_val)
        else:
            threshold = 0.5

        proba_train = booster.predict(dtrain, iteration_range=(0, (getattr(booster, "best_iteration", 0) or 0) + 1))
        y_train_arr = dtrain.get_label().astype(np.uint8, copy=False)
        train_metrics = evaluate_metrics(y_train_arr, proba_train, threshold); del y_train_arr, proba_train; gc.collect()

        if dval.num_row() > 0:
            y_val = dval.get_label().astype(np.uint8, copy=False)
            val_metrics = evaluate_metrics(y_val, proba_val, threshold); del y_val, proba_val; gc.collect()
        else:
            val_metrics = {}

        if dtest.num_row() > 0:
            proba_test = booster.predict(dtest, iteration_range=(0, (getattr(booster, "best_iteration", 0) or 0) + 1))
            y_test_arr = dtest.get_label().astype(np.uint8, copy=False)
            test_metrics = evaluate_metrics(y_test_arr, proba_test, threshold); del y_test_arr, proba_test; gc.collect()
        else:
            test_metrics = {}

        score = booster.get_score(importance_type="gain")
        importance_df = pd.DataFrame({
            "feature": feature_cols,
            "importance": [score.get(f, 0.0) for f in feature_cols]
        }).sort_values("importance", ascending=False)

        model_obj = BoosterPredictor(
            booster=booster,
            best_iteration=getattr(booster, "best_iteration", None),
            feature_names=feature_cols,
            importance_map=score
        )

        del dtrain, dval, dtest; gc.collect()
        return model_obj, threshold, train_metrics, val_metrics, test_metrics, importance_df

    # In-memory fallback
    X_train, y_train = load_Xy_for_watersheds(
        dataset, train_ws, feature_cols, means, stds,
        batch_size=cfg.batch_size, random_state=fold_seed
    )
    X_val, y_val = load_Xy_for_watersheds(
        dataset, val_ws, feature_cols, means, stds,
        batch_size=cfg.batch_size
    )
    X_test, y_test = load_Xy_for_watersheds(
        dataset, test_ws, feature_cols, means, stds,
        batch_size=cfg.batch_size
    )

    # Compute adaptive scale_pos_weight
    spw = compute_adaptive_scale_pos_weight(
        y_train=y_train,
        method=cfg.scale_pos_weight_method,
        cap=cfg.scale_pos_weight_cap
    )
    
    n_pos = int((y_train == 1).sum())
    n_neg = int((y_train == 0).sum())
    natural_ratio = n_neg / max(n_pos, 1)
    print(f"    Natural ratio: {natural_ratio:.1f}:1, Adaptive weight: {spw:.2f}")

    xgb_model = xgb.XGBClassifier(
        n_estimators=cfg.n_estimators,
        max_depth=cfg.max_depth,
        learning_rate=cfg.learning_rate,
        subsample=cfg.subsample,
        colsample_bytree=cfg.colsample_bytree,
        min_child_weight=cfg.min_child_weight,
        gamma=cfg.gamma,
        reg_alpha=cfg.reg_alpha,
        reg_lambda=cfg.reg_lambda,
        objective=cfg.objective,
        eval_metric=cfg.eval_metric,
        tree_method=(cfg.tree_method or "hist"),
        max_bin=cfg.max_bin,
        grow_policy=cfg.grow_policy,
        scale_pos_weight=spw,
        n_jobs=cfg.n_jobs,
        random_state=fold_seed,
        verbosity=cfg.verbosity
    )
    try:
        xgb_model.fit(
            X_train, y_train,
            eval_set=[(X_val, y_val)],
            early_stopping_rounds=cfg.early_stopping_rounds
        )
    except TypeError:
        callbacks = [XGBEarlyStopping(rounds=cfg.early_stopping_rounds, save_best=True, maximize=False)]
        xgb_model.fit(
            X_train, y_train, eval_set=[(X_val, y_val)], callbacks=callbacks
        )

    proba_val = xgb_model.predict_proba(X_val)[:, 1] if len(X_val) else np.array([])
    threshold = find_optimal_threshold(y_val, proba_val) if len(X_val) else 0.5

    proba_train = xgb_model.predict_proba(X_train)[:, 1] if len(X_train) else np.array([])
    train_metrics = evaluate_metrics(y_train, proba_train, threshold) if len(X_train) else {}
    val_metrics = evaluate_metrics(y_val, proba_val, threshold) if len(X_val) else {}
    proba_test = xgb_model.predict_proba(X_test)[:, 1] if len(X_test) else np.array([])
    test_metrics = evaluate_metrics(y_test, proba_test, threshold) if len(X_test) else {}

    importance_df = pd.DataFrame({
        "feature": feature_cols,
        "importance": xgb_model.feature_importances_
    }).sort_values("importance", ascending=False)

    del X_train, y_train, X_val, y_val, X_test, y_test; gc.collect()

    return xgb_model, threshold, train_metrics, val_metrics, test_metrics, importance_df


def train_xgboost_spatial_kfold(
    parquet_uri: str,
    config: XGBoostConfig = None,
    n_folds: int = 5,
    force_test_ws: Optional[str] = None,
    save_path: Optional[str] = None
) -> Dict:
    """
    K-folds spatial CV over watersheds with random 60/20/20 splits per fold.
    Memory-bounded training per fold with adaptive scale_pos_weight tuning.

    Returns:
      - model: final model trained on all watersheds except force_test_ws (if provided)
      - threshold: mean of fold validation-optimal thresholds
      - metrics: dict with mean train/val/test metrics across folds
      - fold_details: per-fold metrics as DataFrame
      - feature_importances: from the final model
    """
    t0 = time.time()
    cfg = config or XGBoostConfig()

    dataset = _open_dataset(parquet_uri)
    watersheds = get_unique_watersheds(dataset)
    if len(watersheds) < 3:
        raise RuntimeError("Insufficient watersheds: need >=3 for k-fold spatial CV.")

    print("=" * 72)
    print(f"XGBOOST FLOOD MODEL – Spatial K-Folds ({n_folds} folds, 60/20/20 per fold)")
    print("=" * 72)
    print(f"All watersheds ({len(watersheds)}): {watersheds}")
    print(f"Adaptive weight method: {cfg.scale_pos_weight_method}")
    if force_test_ws:
        print(f"Force this watershed into TEST for each fold: {force_test_ws}")

    feature_cols = get_feature_columns(dataset)
    print(f"\nUsing {len(feature_cols)} features")
    print("First 10:", feature_cols[:10])

    # build splits
    splits = _spatial_split_kfolds(
        watersheds=watersheds,
        n_folds=n_folds,
        random_state=cfg.random_state,
        train_frac=0.6,
        val_frac=0.2,
        force_test_ws=force_test_ws
    )

    fold_rows = []
    fold_thresholds = []

    for f, (train_ws, val_ws, test_ws) in enumerate(splits, start=1):
        print(f"\nFold {f}/{n_folds}")
        print(f"  Train ({len(train_ws)}): {train_ws}")
        print(f"  Val   ({len(val_ws)}): {val_ws}")
        print(f"  Test  ({len(test_ws)}): {test_ws}")

        model_obj, thr, tr_m, va_m, te_m, imp_df = _train_one_split(
            dataset=dataset,
            feature_cols=feature_cols,
            cfg=cfg,
            train_ws=train_ws,
            val_ws=val_ws,
            test_ws=test_ws,
            fold_seed=cfg.random_state + f
        )
        fold_thresholds.append(thr)

        def _row(label, m):
            row = {"fold": f, "set": label, "threshold": thr}
            row.update({k: float(v) if isinstance(v, (int, float, np.floating)) else v for k, v in m.items()})
            return row

        if tr_m: fold_rows.append(_row("train", tr_m))
        if va_m: fold_rows.append(_row("val", va_m))
        if te_m: fold_rows.append(_row("test", te_m))

        gc.collect()

    folds_df = pd.DataFrame(fold_rows)
    print("\nPer-fold metrics (head):")
    if not folds_df.empty:
        print(folds_df.head())

    # aggregate means
    metrics_mean = {}
    for s in ["train", "val", "test"]:
        sub = folds_df[folds_df["set"] == s]
        if not sub.empty:
            metrics_mean[s] = sub.drop(columns=["fold", "set"]).mean(numeric_only=True).to_dict()
        else:
            metrics_mean[s] = {}

    final_threshold = float(np.mean(fold_thresholds)) if fold_thresholds else 0.5
    print(f"\nMean threshold across folds: {final_threshold:.3f}")

    # Final model training on all watersheds except the forced test ws (if provided)
    if force_test_ws and force_test_ws in watersheds:
        final_train_ws = [w for w in watersheds if w != force_test_ws]
        print(f"Final training on {len(final_train_ws)} watersheds (excluding forced test: {force_test_ws})")
    else:
        final_train_ws = watersheds[:]
        print(f"Final training on all {len(final_train_ws)} watersheds")

    means_final, stds_final = compute_non_test_mean_std(dataset, feature_cols, final_train_ws, batch_size=cfg.batch_size)

    if cfg.use_external_memory:
        try:
            dtrain_all = build_quantile_dmatrix_from_iter(
                dataset, final_train_ws, feature_cols, means_final, stds_final,
                batch_size=cfg.batch_size, random_state=cfg.random_state, max_bin=cfg.max_bin
            )
            qdm_ok = True
        except ValueError:
            dtrain_all = build_dmatrix_from_iter(
                dataset, final_train_ws, feature_cols, means_final, stds_final,
                batch_size=cfg.batch_size, random_state=cfg.random_state,
                row_cap=cfg.max_train_rows
            )
            qdm_ok = False

        ytr = dtrain_all.get_label()
        
        # Compute adaptive scale_pos_weight for final model
        spw = compute_adaptive_scale_pos_weight(
            y_train=ytr,
            method=cfg.scale_pos_weight_method,
            cap=cfg.scale_pos_weight_cap
        )
        
        n_pos = int((ytr == 1).sum())
        n_neg = int((ytr == 0).sum())
        natural_ratio = n_neg / max(n_pos, 1)
        print(f"Final model - Natural ratio: {natural_ratio:.1f}:1, Adaptive weight: {spw:.2f}")
        
        params = {
            "objective": cfg.objective,
            "eval_metric": cfg.eval_metric,
            "tree_method": cfg.tree_method or "hist",
            "max_depth": cfg.max_depth,
            "learning_rate": cfg.learning_rate,
            "subsample": cfg.subsample,
            "colsample_bytree": cfg.colsample_bytree,
            "min_child_weight": cfg.min_child_weight,
            "gamma": cfg.gamma,
            "reg_alpha": cfg.reg_alpha,
            "reg_lambda": cfg.reg_lambda,
            "scale_pos_weight": spw,
            "verbosity": cfg.verbosity,
            "grow_policy": cfg.grow_policy
        }
        if qdm_ok:
            params["max_bin"] = cfg.max_bin

        booster = xgb.train(params=params, dtrain=dtrain_all, num_boost_round=cfg.n_estimators)
        score = booster.get_score(importance_type="gain")
        importance_df = pd.DataFrame({
            "feature": feature_cols,
            "importance": [score.get(f, 0.0) for f in feature_cols]
        }).sort_values("importance", ascending=False)

        final_model = BoosterPredictor(
            booster=booster,
            best_iteration=None,
            feature_names=feature_cols,
            importance_map=score
        )
        del dtrain_all, ytr; gc.collect()
    else:
        X_train_all, y_train_all = load_Xy_for_watersheds(
            dataset, final_train_ws, feature_cols, means_final, stds_final, 
            batch_size=cfg.batch_size, random_state=cfg.random_state
        )
        
        spw = compute_adaptive_scale_pos_weight(
            y_train=y_train_all,
            method=cfg.scale_pos_weight_method,
            cap=cfg.scale_pos_weight_cap
        )
        
        n_pos = int((y_train_all == 1).sum())
        n_neg = int((y_train_all == 0).sum())
        natural_ratio = n_neg / max(n_pos, 1)
        print(f"Final model - Natural ratio: {natural_ratio:.1f}:1, Adaptive weight: {spw:.2f}")
        
        xgb_final = xgb.XGBClassifier(
            n_estimators=cfg.n_estimators,
            max_depth=cfg.max_depth,
            learning_rate=cfg.learning_rate,
            subsample=cfg.subsample,
            colsample_bytree=cfg.colsample_bytree,
            min_child_weight=cfg.min_child_weight,
            gamma=cfg.gamma,
            reg_alpha=cfg.reg_alpha,
            reg_lambda=cfg.reg_lambda,
            objective=cfg.objective,
            eval_metric=cfg.eval_metric,
            tree_method=(cfg.tree_method or "hist"),
            max_bin=cfg.max_bin,
            grow_policy=cfg.grow_policy,
            scale_pos_weight=spw,
            n_jobs=cfg.n_jobs,
            random_state=cfg.random_state,
            verbosity=cfg.verbosity
        )
        xgb_final.fit(X_train_all, y_train_all)
        importance_df = pd.DataFrame({
            "feature": feature_cols,
            "importance": xgb_final.feature_importances_
        }).sort_values("importance", ascending=False)
        final_model = xgb_final
        del X_train_all, y_train_all; gc.collect()

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed/60:.1f} minutes")

    results = {
        "model": final_model,
        "feature_cols": feature_cols,
        "means": means_final,
        "stds": stds_final,
        "threshold": final_threshold,
        "feature_importances": importance_df,
        "watersheds": watersheds,
        "split": {"kfolds": n_folds, "force_test_ws": force_test_ws},
        "metrics": metrics_mean,
        "fold_details": folds_df,
        "config": cfg
    }
    if save_path:
        with open(save_path, "wb") as f:
            pickle.dump(results, f)
        print(f"Saved k-fold model bundle to: {save_path}")
    return results


# ============================= VISUALIZATION =============================

def _normalize_coord_str(value: str) -> str:
    s = str(value).strip()
    try:
        return f"{float(s):.4f}"
    except Exception:
        return s


def _extract_lat_lon_from_training_path(base_s3_path: str) -> Tuple[str, str]:
    m = re.search(r"trainingFolderFor_Lat_([-\d\.]+)_Lon_([-\d\.]+)", base_s3_path)
    if not m:
        raise ValueError(f"Could not parse Lat/Lon from base_s3_path: {base_s3_path}")
    return m.group(1), m.group(2)


def _configure_gdal_proj():
    try:
        from pyproj import datadir as _pydatadir
        proj_dir = _pydatadir.get_data_dir()
        if proj_dir and os.path.isdir(proj_dir):
            os.environ["PROJ_LIB"] = proj_dir
            gdal.SetConfigOption("PROJ_LIB", proj_dir)
    except Exception:
        pass

    if "PROJ_LIB" not in os.environ or not os.path.isdir(os.environ.get("PROJ_LIB", "")):
        for cand in ["/opt/conda/share/proj", "/usr/share/proj"]:
            if os.path.isdir(cand):
                os.environ["PROJ_LIB"] = cand
                gdal.SetConfigOption("PROJ_LIB", cand)
                break

    try:
        gdal_mod_dir = os.path.dirname(gdal.__file__)
        gdal_data_dir = os.path.join(gdal_mod_dir, "data")
        if os.path.isdir(gdal_data_dir):
            os.environ["GDAL_DATA"] = gdal_data_dir
            gdal.SetConfigOption("GDAL_DATA", gdal_data_dir)
    except Exception:
        pass

    gdal.SetConfigOption("PROJ_NETWORK", "ON")


class XGBoostFloodVisualizer:
    """Generate visual outputs for XGBoost flood prediction model results."""

    def __init__(self, model_results: Dict, parquet_uri: str, base_s3_path: str):
        self.model_results = model_results
        self.model = model_results["model"]
        self.feature_cols = model_results["feature_cols"]
        self.means = model_results["means"]
        self.stds = model_results["stds"]
        self.threshold = model_results["threshold"]

        self.parquet_uri = parquet_uri
        self.base_s3_path = base_s3_path.rstrip('/')

        self.test_lat, self.test_lon = _extract_lat_lon_from_training_path(self.base_s3_path)
        self.test_watershed = f"Lat_{self.test_lat}_Lon_{self.test_lon}"

        self.s3_client = boto3.client('s3')

        _configure_gdal_proj()
        gdal.SetConfigOption('AWS_NO_SIGN_REQUEST', 'NO')
        gdal.SetConfigOption('GDAL_DISABLE_READDIR_ON_OPEN', 'YES')
        gdal.SetConfigOption('CPL_VSIL_CURL_ALLOWED_EXTENSIONS', 'tif,tiff')

        self._static_map: Optional[Dict[str, np.ndarray]] = None
        self._ts_by_date: Optional[Dict[str, Dict[str, float]]] = None

    @staticmethod
    def _sanitize(name: str) -> str:
        name = (name or "").strip()
        name = re.sub(r"[^\w\-]+", "_", name)
        name = re.sub(r"_+", "_", name)
        return name

    def generate_all_outputs(self, output_s3_folder: Optional[str] = None) -> Dict[str, Any]:
        if output_s3_folder is None:
            output_s3_folder = f"{self.base_s3_path}/assessXGBoostModel"

        print("\n" + "="*60)
        print("GENERATING XGBOOST MODEL VISUALIZATION OUTPUTS")
        print("="*60)
        print(f"Test watershed: {self.test_watershed}")
        print(f"Output folder: {output_s3_folder}")

        flood_tiff_path = (
            f"{self.base_s3_path}/floodMaps/"
            f"dswx_s1_timeseries_subwatershed_Lon{self.test_lon}_Lat{self.test_lat}.tif"
        )
        print(f"\nLoading flood maps from: {flood_tiff_path}")
        flood_ds, flood_data, dates = self._load_flood_tiff(flood_tiff_path)

        print("Preloading static features from rasters/ and time series CSVs...")
        self._preload_static_map()
        self._preload_timeseries()

        print("\nGenerating XGBoost predictions and confusion matrices...")
        predictions, confusions, per_date_metrics = self._generate_predictions(flood_data, dates)

        print("\nCreating output TIFFs...")
        pred_tiff_path = f"{output_s3_folder}/xgboost_predictions_{self.test_watershed}.tif"
        cm_tiff_path = f"{output_s3_folder}/xgboost_confusion_matrix_{self.test_watershed}.tif"
        self._write_multiband_tiff(flood_ds, predictions, dates, pred_tiff_path, prefix="XGB_Prediction")
        self._write_multiband_tiff(flood_ds, confusions, dates, cm_tiff_path, prefix="XGB_CM_0TN_1TP_2FN_3FP")

        print("\nCreating PNG visualizations...")
        png_paths = self._create_png_visualizations(confusions, dates, flood_data, predictions, per_date_metrics, output_s3_folder)

        print("\nCreating CSV summaries (metrics and importances)...")
        metrics_summary = self.model_results.get("metrics", {})
        feature_importances = self.model_results.get("feature_importances", pd.DataFrame(columns=["feature","importance"]))
        csv_paths = self._write_csv_results(output_s3_folder, metrics_summary, feature_importances, per_date_metrics)

        print("\n✅ XGBoost visualization complete!")
        return {
            "prediction_tiff": pred_tiff_path,
            "confusion_matrix_tiff": cm_tiff_path,
            "png_visualizations": png_paths,
            "csv_outputs": csv_paths
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
            m = re.search(r"\d{4}-\d{2}-\d{2}", desc)
            date_str = m.group(0) if m else f"band_{b}"
            dates.append(date_str)
            flood_data[date_str] = band.ReadAsArray()

        return ds, flood_data, dates

    def _preload_static_map(self):
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
            mask |= (arr == -9999.0)
            mask |= (arr == -32768.0)
            arr = np.where(mask, np.nan, arr)
            static_map[nm] = arr.ravel()
        ds = None
        self._static_map = static_map

    def _preload_timeseries(self):
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
                    if col in df.columns:
                        d[f"{prefix}_{col}"] = float(row[col]) if pd.notna(row[col]) else np.nan
                for col in df.columns:
                    if col != 'Date' and col not in ('streamflow_pred', 't2m', 'ssrd', 'tp'):
                        d[f"{prefix}_{col}"] = float(row[col]) if pd.notna(row[col]) else np.nan

        sub_csv = (
            f"{self.base_s3_path}/streamflowPredictions/"
            f"subwatershed_predictions_Lat{self.test_lat}_Lon{self.test_lon}.csv"
        )
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
        """
        Generate per-date predictions and confusion matrices for the target watershed.
        Updated to match GLM metrics exactly.
        """
        predictions: Dict[str, np.ndarray] = {}
        confusions: Dict[str, np.ndarray] = {}
        per_date_metrics: Dict[str, Dict[str, float]] = {}

        first = next(iter(flood_data.values()))
        n_rows, n_cols = first.shape
        n_pix = n_rows * n_cols
        print(f"  Processing {len(dates)} dates with XGBoost, {n_pix:,} pixels per date")

        ts_cols = [c for c in self.feature_cols if c.startswith("subwatershed_") or c.startswith("river_basin_")]
        static_cols = [c for c in self.feature_cols if c not in ts_cols]

        for idx, date_str in enumerate(dates):
            if (idx + 1) % 10 == 0:
                print(f"  Date {idx+1}/{len(dates)}: {date_str}")

            arr = flood_data[date_str]
            valid_mask = (arr != 255)
            valid_idx = np.where(valid_mask.ravel())[0]
            if valid_idx.size == 0:
                predictions[date_str] = np.full((n_rows, n_cols), 255, dtype=np.uint8)
                confusions[date_str] = np.full((n_rows, n_cols), 255, dtype=np.uint8)
                per_date_metrics[date_str] = {
                    "date": date_str,  # ADDED
                    "F1": np.nan, 
                    "Accuracy": np.nan, 
                    "Precision": np.nan, 
                    "Recall": np.nan, 
                    "ROC_AUC": np.nan,
                    "PR_AUC": np.nan,  # ADDED
                    "Kappa": np.nan,  # ADDED
                    "MCC": np.nan,  # ADDED
                    "actual_pos_frac": np.nan, 
                    "pred_pos_frac": np.nan
                }
                continue

            # Assemble feature matrix for valid pixels
            X = np.empty((valid_idx.size, len(self.feature_cols)), dtype=np.float32)

            # Static features (per-pixel)
            for j, name in enumerate(self.feature_cols):
                if name in static_cols:
                    if (self._static_map is not None) and (name in self._static_map):
                        X[:, j] = self._static_map[name][valid_idx]
                    else:
                        X[:, j] = np.nan

            # Dynamic (per-date) features (broadcast across valid pixels)
            ts_vals = self._ts_by_date.get(date_str, {}) if self._ts_by_date else {}
            for j, name in enumerate(self.feature_cols):
                if name in ts_cols:
                    X[:, j] = ts_vals.get(name, np.nan)

            # Impute + standardize
            X = np.where(np.isnan(X), self.means, X)
            X = (X - self.means) / self.stds

            # Predict probabilities
            if hasattr(self.model, "predict_proba"):
                y_proba = self.model.predict_proba(X)
                if y_proba.ndim == 2:
                    y_proba = y_proba[:, 1]
            else:
                dm = xgb.DMatrix(X, feature_names=self.feature_cols)
                y_proba = self.model.predict(dm)

            # Threshold
            y_pred = (y_proba >= self.threshold).astype(np.uint8)
            y_true = (arr.ravel()[valid_idx] > 0).astype(np.uint8)

            # Full-frame predictions (with nodata=255)
            pred_full = np.full(n_pix, 255, dtype=np.uint8)
            pred_full[valid_idx] = y_pred
            predictions[date_str] = pred_full.reshape(n_rows, n_cols)

            # Confusion map (0 TN, 1 TP, 2 FN, 3 FP, 255 nodata)
            cm_full = np.full(n_pix, 255, dtype=np.uint8)
            cm_valid = np.where(
                y_true == 0,
                np.where(y_pred == 0, 0, 3),
                np.where(y_pred == 1, 1, 2)
            ).astype(np.uint8)
            cm_full[valid_idx] = cm_valid
            confusions[date_str] = cm_full.reshape(n_rows, n_cols)

            # Calculate metrics - MATCHING GLM EXACTLY
            
            # ROC AUC with error handling
            try:
                roc_auc = roc_auc_score(y_true, y_proba)
            except:
                roc_auc = np.nan

            # PR AUC with error handling (ADDED)
            try:
                pr_auc = average_precision_score(y_true, y_proba)
            except:
                pr_auc = np.nan

            # MCC with error handling (ADDED)
            if len(np.unique(y_true)) > 1 and len(np.unique(y_pred)) > 1:
                mcc = matthews_corrcoef(y_true, y_pred)
            else:
                mcc = 0.0

            total_valid = float(y_true.size)
            actual_pos_frac = float(y_true.mean()) if total_valid > 0 else np.nan
            pred_pos_frac = float(y_pred.mean()) if total_valid > 0 else np.nan

            per_date_metrics[date_str] = {
                "date": date_str,  # ADDED
                "F1": f1_score(y_true, y_pred, zero_division=0),
                "Accuracy": accuracy_score(y_true, y_pred),
                "Precision": precision_score(y_true, y_pred, zero_division=0),
                "Recall": recall_score(y_true, y_pred, zero_division=0),
                "ROC_AUC": float(roc_auc) if np.isfinite(roc_auc) else np.nan,
                "PR_AUC": float(pr_auc) if np.isfinite(pr_auc) else np.nan,  # ADDED
                "Kappa": cohen_kappa_score(y_true, y_pred),  # ADDED
                "MCC": float(mcc),  # ADDED
                "actual_pos_frac": actual_pos_frac,
                "pred_pos_frac": pred_pos_frac
            }

        return predictions, confusions, per_date_metrics
        
    def _write_multiband_tiff(
        self,
        reference_ds: gdal.Dataset,
        band_arrays: Dict[str, np.ndarray],
        dates: List[str],
        output_s3_path: str,
        prefix: str
    ):
        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmp:
            temp_path = tmp.name

        n_bands = len(dates)
        n_rows = reference_ds.RasterYSize
        n_cols = reference_ds.RasterXSize

        driver = gdal.GetDriverByName('GTiff')
        out = driver.Create(
            temp_path, n_cols, n_rows, n_bands, gdal.GDT_Byte,
            options=["COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=IF_SAFER"]
        )
        try:
            out.SetGeoTransform(reference_ds.GetGeoTransform())
        except Exception:
            pass
        try:
            out.SetProjection(reference_ds.GetProjection())
        except Exception:
            pass

        for i, date_str in enumerate(dates, start=1):
            band = out.GetRasterBand(i)
            band.SetNoDataValue(255)
            band.SetDescription(f"{prefix}_{date_str}")
            band.WriteArray(band_arrays.get(date_str, np.full((n_rows, n_cols), 255, dtype=np.uint8)))
            band.FlushCache()

        out.FlushCache()
        out = None

        self._upload_to_s3(temp_path, output_s3_path)
        os.remove(temp_path)
        print(f"  Created TIFF: {output_s3_path}")

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
            ax.set_title(f'XGBoost Confusion Matrix - {date_str}', fontsize=14, fontweight='bold')

            text = (
                f"Date: {date_str}\n"
                f"Model: XGBoost\n"
                f"F1 Score: {metrics.get('F1', np.nan):.3f}\n"
                f"Accuracy: {metrics.get('Accuracy', np.nan):.3f}\n"
                f"Precision: {metrics.get('Precision', np.nan):.3f}\n"
                f"Recall: {metrics.get('Recall', np.nan):.3f}\n"
                f"ROC AUC: {auc_str}"
            )
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            ax.text(0.02, 0.98, text, transform=ax.transAxes, fontsize=11, va='top', bbox=props)
            ax.axis('off')
            plt.tight_layout()

            with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
                temp_png = tmp.name
            plt.savefig(temp_png, dpi=150, bbox_inches='tight')
            plt.close()

            s3_path = f"{output_s3_folder}/xgboost_confusion_matrix_{date_str}.png"
            self._upload_to_s3(temp_png, s3_path)
            os.remove(temp_png)
            png_paths.append(s3_path)

            if (i + 1) % 10 == 0:
                print(f"  Created {i + 1} PNGs")

        print(f"  Created {len(png_paths)} XGBoost PNG visualizations")
        return png_paths

    def _upload_to_s3(self, local_path: str, s3_path: str):
        bucket, key = s3_path.replace('s3://', '').split('/', 1)
        self.s3_client.upload_file(local_path, bucket, key)

    def _write_csv_results(
        self,
        output_s3_folder: str,
        metrics_summary: Dict[str, Dict[str, float]],
        feature_importances: pd.DataFrame,
        per_date_metrics: Dict[str, Dict[str, float]]
    ) -> Dict[str, str]:
        paths = {}

        metrics_df = pd.DataFrame.from_dict(metrics_summary, orient="index").reset_index().rename(columns={"index": "set"})
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
            tmp_csv = tmp.name
        metrics_df.to_csv(tmp_csv, index=False)
        s3_path = f"{output_s3_folder}/xgb_metrics_summary.csv"
        self._upload_to_s3(tmp_csv, s3_path)
        os.remove(tmp_csv)
        paths["metrics_summary_csv"] = s3_path

        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
            tmp_csv = tmp.name
        feature_importances.to_csv(tmp_csv, index=False)
        s3_path = f"{output_s3_folder}/xgb_feature_importances.csv"
        self._upload_to_s3(tmp_csv, s3_path)
        os.remove(tmp_csv)
        paths["feature_importances_csv"] = s3_path

        if per_date_metrics:
            per_date_df = pd.DataFrame.from_dict(per_date_metrics, orient="index")
            # Ensure consistent column order matching GLM
            column_order = ['date', 'F1', 'Accuracy', 'Precision', 'Recall', 'ROC_AUC', 
                           'PR_AUC', 'Kappa', 'MCC', 'actual_pos_frac', 'pred_pos_frac']
            per_date_df = per_date_df[column_order]
            
            with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
                tmp_csv = tmp.name
            per_date_df.to_csv(tmp_csv, index=False)
            s3_path = f"{output_s3_folder}/xgboost_per_date_metrics_{self.test_watershed}.csv"
            self._upload_to_s3(tmp_csv, s3_path)
            os.remove(tmp_csv)
            paths["per_date_metrics_csv"] = s3_path

        return paths


def visualize_xgboost_outputs(
    model_results: Dict,
    parquet_uri: str,
    latitude: str,
    longitude: str
) -> Dict[str, Any]:
    lat_str = _normalize_coord_str(latitude)
    lon_str = _normalize_coord_str(longitude)

    base_s3_path = (
        f"s3://climate-ai-data-science-datasets/arrakis-data/floodOutputs/"
        f"trainingFolderFor_Lat_{lat_str}_Lon_{lon_str}"
    )
    output_folder = f"{base_s3_path}/assessXGBoostModel"

    _configure_gdal_proj()
    viz = XGBoostFloodVisualizer(model_results=model_results, parquet_uri=parquet_uri, base_s3_path=base_s3_path)
    return viz.generate_all_outputs(output_folder)

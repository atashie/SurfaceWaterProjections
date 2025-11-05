# decision_tree_pipeline.py

import numpy as np
import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import (
    f1_score, accuracy_score, precision_score, recall_score,
    roc_auc_score, confusion_matrix, classification_report
)
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
import json
from datetime import datetime

logger = logging.getLogger(__name__)


class DecisionTreeFloodModel:
    """
    Decision Tree model for flood prediction with leave-one-watershed-out CV.
    """
    
    def __init__(self, 
                 max_depth: Optional[int] = None,
                 min_samples_split: int = 100,
                 min_samples_leaf: int = 50,
                 random_state: int = 42):
        """
        Initialize Decision Tree model.
        
        Parameters:
        -----------
        max_depth : int or None
            Maximum depth of the tree
        min_samples_split : int
            Minimum samples required to split internal node
        min_samples_leaf : int
            Minimum samples required at leaf node
        random_state : int
            Random seed for reproducibility
        """
        self.model = DecisionTreeClassifier(
            max_depth=max_depth,
            min_samples_split=min_samples_split,
            min_samples_leaf=min_samples_leaf,
            random_state=random_state,
            class_weight='balanced'  # Handle class imbalance
        )
        
        self.scaler = StandardScaler()
        self.imputer = SimpleImputer(strategy='median')
        self.feature_columns = None
        self.metrics_history = []
        
    def prepare_features(self, df: pd.DataFrame, 
                        fit_preprocessors: bool = False) -> np.ndarray:
        """
        Prepare features for model training/inference.
        
        Parameters:
        -----------
        df : pd.DataFrame
            Feature dataframe
        fit_preprocessors : bool
            Whether to fit preprocessors (True for training)
            
        Returns:
        --------
        Processed feature array
        """
        # Select feature columns (exclude metadata and target)
        exclude_cols = ['pixel_idx', 'timestep_idx', 'date', 
                       'watershed_id', 'target']
        
        if self.feature_columns is None:
            self.feature_columns = [col for col in df.columns 
                                   if col not in exclude_cols]
            logger.info(f"Selected {len(self.feature_columns)} features")
        
        X = df[self.feature_columns].values
        
        # Impute missing values
        if fit_preprocessors:
            X = self.imputer.fit_transform(X)
        else:
            X = self.imputer.transform(X)
        
        # Scale features
        if fit_preprocessors:
            X = self.scaler.fit_transform(X)
        else:
            X = self.scaler.transform(X)
        
        return X
    
    def train(self, train_df: pd.DataFrame):
        """Train the model on training data."""
        
        logger.info(f"Training on {len(train_df)} samples")
        
        # Prepare features and target
        X_train = self.prepare_features(train_df, fit_preprocessors=True)
        y_train = train_df['target'].values
        
        # Remove samples with invalid targets
        valid_mask = ~np.isnan(y_train)
        X_train = X_train[valid_mask]
        y_train = y_train[valid_mask]
        
        # Train model
        self.model.fit(X_train, y_train)
        
        # Store feature importance
        self.feature_importance = pd.DataFrame({
            'feature': self.feature_columns,
            'importance': self.model.feature_importances_
        }).sort_values('importance', ascending=False)
        
        logger.info("Model training complete")
    
    def predict(self, test_df: pd.DataFrame) -> Dict:
        """
        Make predictions and compute metrics.
        
        Returns:
        --------
        Dictionary with predictions and metrics
        """
        # Prepare features
        X_test = self.prepare_features(test_df, fit_preprocessors=False)
        y_true = test_df['target'].values
        
        # Remove invalid samples
        valid_mask = ~np.isnan(y_true)
        X_test = X_test[valid_mask]
        y_true = y_true[valid_mask]
        
        # Make predictions
        y_pred = self.model.predict(X_test)
        y_prob = self.model.predict_proba(X_test)[:, 1]
        
        # Compute metrics
        metrics = {
            'f1_score': f1_score(y_true, y_pred, average='binary'),
            'accuracy': accuracy_score(y_true, y_pred),
            'precision': precision_score(y_true, y_pred, average='binary'),
            'recall': recall_score(y_true, y_pred, average='binary'),
            'roc_auc': roc_auc_score(y_true, y_prob) if len(np.unique(y_true)) > 1 else 0.0,
            'n_samples': len(y_true),
            'n_positive': np.sum(y_true == 1),
            'n_negative': np.sum(y_true == 0)
        }
        
        # Add confusion matrix
        cm = confusion_matrix(y_true, y_pred)
        metrics['confusion_matrix'] = cm.tolist()
        
        return {
            'y_true': y_true,
            'y_pred': y_pred,
            'y_prob': y_prob,
            'metrics': metrics,
            'test_indices': test_df.index[valid_mask].tolist()
        }
    
    def save_model(self, output_dir: Path):
        """Save model and preprocessors."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save model
        joblib.dump(self.model, output_dir / 'decision_tree_model.pkl')
        
        # Save preprocessors
        joblib.dump(self.scaler, output_dir / 'scaler.pkl')
        joblib.dump(self.imputer, output_dir / 'imputer.pkl')
        
        # Save feature columns and importance
        with open(output_dir / 'feature_columns.json', 'w') as f:
            json.dump(self.feature_columns, f)
        
        self.feature_importance.to_csv(
            output_dir / 'feature_importance.csv', index=False
        )
        
        logger.info(f"Model saved to {output_dir}")
    
    def load_model(self, model_dir: Path):
        """Load saved model and preprocessors."""
        model_dir = Path(model_dir)
        
        self.model = joblib.load(model_dir / 'decision_tree_model.pkl')
        self.scaler = joblib.load(model_dir / 'scaler.pkl')
        self.imputer = joblib.load(model_dir / 'imputer.pkl')
        
        with open(model_dir / 'feature_columns.json', 'r') as f:
            self.feature_columns = json.load(f)
        
        logger.info(f"Model loaded from {model_dir}")


class LeaveOneWatershedOutCV:
    """
    Implements leave-one-watershed-out cross-validation.
    """
    
    def __init__(self, watershed_data: Dict[str, pd.DataFrame]):
        """
        Initialize CV splitter.
        
        Parameters:
        -----------
        watershed_data : Dict[str, pd.DataFrame]
            Dictionary mapping watershed IDs to feature DataFrames
        """
        self.watershed_data = watershed_data
        self.watershed_ids = list(watershed_data.keys())
        
        # Identify initial watershed (test set)
        self.test_watershed = None
        for wid, df in watershed_data.items():
            if 'is_initial' in df.columns and df['is_initial'].any():
                self.test_watershed = wid
                break
        
        if self.test_watershed is None:
            # Use first watershed as test if not specified
            self.test_watershed = self.watershed_ids[0]
        
        self.train_watersheds = [w for w in self.watershed_ids 
                                if w != self.test_watershed]
        
        logger.info(f"Test watershed: {self.test_watershed}")
        logger.info(f"Training watersheds: {self.train_watersheds}")
    
    def get_cv_splits(self) -> List[Dict]:
        """
        Generate CV splits.
        
        Returns:
        --------
        List of split dictionaries
        """
        splits = []
        
        # Create leave-one-out splits for training watersheds
        for i, val_watershed in enumerate(self.train_watersheds):
            train_watersheds = [w for w in self.train_watersheds 
                              if w != val_watershed]
            
            split = {
                'fold': i,
                'train': train_watersheds,
                'val': val_watershed,
                'test': self.test_watershed
            }
            splits.append(split)
        
        return splits
    
    def run_cv(self, model_class=DecisionTreeFloodModel, 
              model_params: Dict = None) -> Dict:
        """
        Run complete cross-validation.
        
        Parameters:
        -----------
        model_class : class
            Model class to instantiate
        model_params : Dict
            Parameters for model initialization
            
        Returns:
        --------
        CV results dictionary
        """
        if model_params is None:
            model_params = {}
        
        splits = self.get_cv_splits()
        cv_results = {
            'fold_results': [],
            'test_results': None,
            'summary': {}
        }
        
        # Run each fold
        for split in splits:
            logger.info(f"\n=== Fold {split['fold']} ===")
            logger.info(f"Train: {split['train']}")
            logger.info(f"Val: {split['val']}")
            
            # Prepare data
            train_dfs = [self.watershed_data[w] for w in split['train']]
            train_df = pd.concat(train_dfs, ignore_index=True)
            
            val_df = self.watershed_data[split['val']]
            
            # Train model
            model = model_class(**model_params)
            model.train(train_df)
            
            # Validate
            val_results = model.predict(val_df)
            val_results['fold'] = split['fold']
            val_results['val_watershed'] = split['val']
            
            cv_results['fold_results'].append(val_results)
            
            logger.info(f"Validation F1: {val_results['metrics']['f1_score']:.4f}")
        
        # Train final model on all training data
        logger.info("\n=== Final Model Training ===")
        all_train_dfs = [self.watershed_data[w] for w in self.train_watersheds]
        all_train_df = pd.concat(all_train_dfs, ignore_index=True)
        
        final_model = model_class(**model_params)
        final_model.train(all_train_df)
        
        # Test on held-out watershed
        test_df = self.watershed_data[self.test_watershed]
        test_results = final_model.predict(test_df)
        test_results['test_watershed'] = self.test_watershed
        
        cv_results['test_results'] = test_results
        cv_results['final_model'] = final_model
        
        # Compute summary statistics
        cv_results['summary'] = self._compute_summary(cv_results)
        
        logger.info(f"\nTest F1: {test_results['metrics']['f1_score']:.4f}")
        
        return cv_results
    
    def _compute_summary(self, cv_results: Dict) -> Dict:
        """Compute summary statistics across folds."""
        
        # Validation metrics across folds
        val_metrics = {}
        for metric in ['f1_score', 'accuracy', 'precision', 'recall', 'roc_auc']:
            values = [fold['metrics'][metric] 
                     for fold in cv_results['fold_results']]
            val_metrics[f'val_{metric}_mean'] = np.mean(values)
            val_metrics[f'val_{metric}_std'] = np.std(values)
        
        # Test metrics
        test_metrics = cv_results['test_results']['metrics']
        
        summary = {
            'n_folds': len(cv_results['fold_results']),
            'validation': val_metrics,
            'test': test_metrics,
            'timestamp': datetime.now().isoformat()
        }
        
        return summary

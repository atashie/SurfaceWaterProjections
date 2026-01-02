###############################################################
###############################################################
###############################################################
### Out of Basin Transfer learning
###############################################################
###############################################################
###############################################################




"""
visualize_out_of_basin_transfer.py

Generate boxplots comparing within-basin performance to out-of-basin transfer performance.
Creates 7 figures (one per metric), each with 2x4 subplots (one per target location).
Each subplot shows 8 boxplots: 1 within-basin + 7 out-of-basin sources.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from typing import Dict, List, Tuple
import boto3
import tempfile


# ======================== Configuration ========================

basin_locations = [
    {"latitude": "10.7875",  "longitude": "0.7479",  "name": "EMSR470"},
    {"latitude": "48.1747",  "longitude": "23.3978", "name": "EMSR444"},
    {"latitude": "53.6583",  "longitude": "-1.1398", "name": "EMSR407"},
    {"latitude": "-16.5646", "longitude": "46.7857", "name": "EMSR424"},
    {"latitude": "-1.8038",  "longitude": "33.3857", "name": "EMSR438"},
    {"latitude": "40.9437",  "longitude": "24.7669", "name": "EMSR292"},
    {"latitude": "27.2981",  "longitude": "68.1875", "name": "EMSR629"},
    {"latitude": "41.9275",  "longitude": "-1.3792", "name": "EMSR279"},
]

METRICS = ['MCC', 'F1', 'Accuracy', 'Precision', 'Recall', 'ROC_AUC', 'PR_AUC']


# ======================== Data Loading Functions ========================

def normalize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names to handle both uppercase and lowercase."""
    # Create mapping from lowercase to preferred case
    column_mapping = {}
    for col in df.columns:
        col_lower = col.lower()
        if col_lower == 'mcc':
            column_mapping[col] = 'MCC'
        elif col_lower == 'f1':
            column_mapping[col] = 'F1'
        elif col_lower == 'accuracy':
            column_mapping[col] = 'Accuracy'
        elif col_lower == 'precision':
            column_mapping[col] = 'Precision'
        elif col_lower == 'recall':
            column_mapping[col] = 'Recall'
        elif col_lower == 'roc_auc':
            column_mapping[col] = 'ROC_AUC'
        elif col_lower == 'pr_auc':
            column_mapping[col] = 'PR_AUC'
        elif col_lower == 'date':
            column_mapping[col] = 'date'
    
    return df.rename(columns=column_mapping)


def load_s3_csv(s3_path: str) -> pd.DataFrame:
    """Load CSV from S3 and normalize column names."""
    s3_client = boto3.client('s3')
    bucket, key = s3_path.replace('s3://', '').split('/', 1)
    
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
        temp_path = tmp.name
    
    try:
        s3_client.download_file(bucket, key, temp_path)
        df = pd.read_csv(temp_path)
        df = normalize_column_names(df)
        os.remove(temp_path)
        return df
    except Exception as e:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        print(f"    ✗ Failed to load {s3_path}: {e}")
        return None


def load_within_basin_results(location: Dict, base_prefix: str = "arrakis-data/floodOutputs") -> pd.DataFrame:
    """Load within-basin test results."""
    lat = location['latitude']
    lon = location['longitude']
    
    csv_path = (
        f"s3://climate-ai-data-science-datasets/{base_prefix}/"
        f"trainingFolderFor_Lat_{lat}_Lon_{lon}/"
        f"assessXGBoostModel/"
        f"xgboost_per_date_metrics_Lat_{lat}_Lon_{lon}.csv"
    )
    
    return load_s3_csv(csv_path)


def load_out_of_basin_results(target_location: Dict, 
                               all_locations: List[Dict],
                               base_prefix: str = "arrakis-data/floodOutputs") -> Dict[str, pd.DataFrame]:
    """Load out-of-basin transfer test results."""
    target_lat = target_location['latitude']
    target_lon = target_location['longitude']
    
    source_locations = [loc for loc in all_locations if loc != target_location]
    
    results = {}
    
    for source_loc in source_locations:
        source_lat = source_loc['latitude']
        source_lon = source_loc['longitude']
        
        csv_path = (
            f"s3://climate-ai-data-science-datasets/{base_prefix}/"
            f"trainingFolderFor_Lat_{target_lat}_Lon_{target_lon}/"
            f"outOfBasinsTests/"
            f"xgboost_Lat_{source_lat}_Lon_{source_lon}_per_day_metrics_Lat_{target_lat}_Lon_{target_lon}.csv"
        )
        
        df = load_s3_csv(csv_path)
        
        if df is not None:
            results[source_loc['name']] = df
    
    return results


# ======================== Data Collection ========================

def collect_all_data(basin_locations: List[Dict], 
                     base_prefix: str = "arrakis-data/floodOutputs") -> Dict:
    """Collect all within-basin and out-of-basin data for all locations.
    
    Returns:
        Dict with structure:
        {
            'EMSR470': {
                'within': DataFrame,
                'out_of_basin': {'EMSR444': DataFrame, 'EMSR407': DataFrame, ...}
            },
            ...
        }
    """
    
    print("=" * 70)
    print("COLLECTING ALL TRANSFER LEARNING DATA")
    print("=" * 70)
    
    all_data = {}
    
    for target_loc in basin_locations:
        print(f"\nTarget: {target_loc['name']}")
        
        # Load within-basin
        print(f"  Loading within-basin results...")
        within_df = load_within_basin_results(target_loc, base_prefix)
        
        # Load out-of-basin
        print(f"  Loading out-of-basin results...")
        out_of_basin_results = load_out_of_basin_results(target_loc, basin_locations, base_prefix)
        
        all_data[target_loc['name']] = {
            'within': within_df,
            'out_of_basin': out_of_basin_results
        }
    
    return all_data


# ======================== Statistical Testing ========================

def perform_statistical_test(within_values: np.ndarray, 
                             out_values: np.ndarray,
                             alpha: float = 0.05) -> Tuple[bool, float]:
    """Perform Mann-Whitney U test."""
    if len(within_values) < 2 or len(out_values) < 2:
        return False, 1.0
    
    within_clean = within_values[~np.isnan(within_values)]
    out_clean = out_values[~np.isnan(out_values)]
    
    if len(within_clean) < 2 or len(out_clean) < 2:
        return False, 1.0
    
    try:
        statistic, p_value = stats.mannwhitneyu(within_clean, out_clean, alternative='two-sided')
        return p_value < alpha, p_value
    except:
        return False, 1.0


# ======================== Visualization ========================

def create_aggregate_transfer_figure(all_data: Dict,
                                     metric: str,
                                     basin_locations: List[Dict],
                                     output_path: str):
    """Create one figure with 2x4 subplots (one subplot per target location).
    
    Each subplot shows 8 boxplots: 1 within-basin + 7 out-of-basin sources.
    
    Args:
        all_data: Dict with all collected data
        metric: Which metric to plot
        basin_locations: List of all basin locations
        output_path: Where to save the figure
    """
    
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    axes = axes.flatten()
    
    for idx, target_loc in enumerate(basin_locations):
        ax = axes[idx]
        target_name = target_loc['name']
        
        # Get data for this target location
        target_data = all_data.get(target_name, {})
        within_df = target_data.get('within')
        out_of_basin_dfs = target_data.get('out_of_basin', {})
        
        # Prepare data for this subplot
        boxplot_data = []
        boxplot_labels = []
        edge_colors = []
        edge_widths = []
        
        # 1. Within-basin data
        if within_df is not None and metric in within_df.columns:
            within_values = within_df[metric].dropna().values
            boxplot_data.append(within_values)
            boxplot_labels.append(target_name)
            edge_colors.append('black')
            edge_widths.append(2)
        else:
            within_values = np.array([])
            boxplot_data.append([])
            boxplot_labels.append(target_name)
            edge_colors.append('black')
            edge_widths.append(2)
        
        # 2. Out-of-basin data (7 source models)
        # Get source names in order (excluding target)
        source_names = [loc['name'] for loc in basin_locations if loc['name'] != target_name]
        
        for source_name in source_names:
            out_df = out_of_basin_dfs.get(source_name)
            
            if out_df is not None and metric in out_df.columns:
                out_values = out_df[metric].dropna().values
                boxplot_data.append(out_values)
                
                # Statistical test
                if len(within_values) > 0:
                    is_sig, p_val = perform_statistical_test(within_values, out_values)
                    edge_colors.append('red' if is_sig else 'black')
                    edge_widths.append(3 if is_sig else 2)
                else:
                    edge_colors.append('black')
                    edge_widths.append(2)
            else:
                boxplot_data.append([])
                edge_colors.append('black')
                edge_widths.append(2)
            
            boxplot_labels.append(source_name)
        
        # Create boxplots
        positions = list(range(len(boxplot_data)))
        
        bp = ax.boxplot(boxplot_data, 
                        positions=positions,
                        widths=0.6,
                        patch_artist=True,
                        showfliers=True)
        
        # Style each box
        for patch_idx, (patch, edge_color, edge_width) in enumerate(zip(bp['boxes'], edge_colors, edge_widths)):
            if patch_idx == 0:
                # Within-basin: light blue
                patch.set_facecolor('lightblue')
            else:
                # Out-of-basin: light coral
                patch.set_facecolor('lightcoral')
            
            patch.set_edgecolor(edge_color)
            patch.set_linewidth(edge_width)
        
        # Style medians, whiskers, caps
        for element in bp['medians']:
            element.set_color('darkred')
            element.set_linewidth(2)
        
        for whisker, edge_color, edge_width in zip(bp['whiskers'], 
                                                    [c for c in edge_colors for _ in range(2)],
                                                    [w for w in edge_widths for _ in range(2)]):
            whisker.set_color(edge_color)
            whisker.set_linewidth(edge_width)
        
        for cap, edge_color, edge_width in zip(bp['caps'],
                                               [c for c in edge_colors for _ in range(2)],
                                               [w for w in edge_widths for _ in range(2)]):
            cap.set_color(edge_color)
            cap.set_linewidth(edge_width)
        
        # Labels and formatting
        ax.set_xticks(positions)
        ax.set_xticklabels(boxplot_labels, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel(metric, fontsize=11, fontweight='bold')
        ax.set_title(target_name, fontsize=12, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        
        # Y-axis limits
        if metric == 'MCC':
            ax.set_ylim(-0.1, 1.0)
        else:
            ax.set_ylim(0, 1.0)
    
    # Overall title
    fig.suptitle(
        f'Transfer Learning Analysis: {metric}\n'
        f'Each subplot = 1 target basin | X-axis = source basin (where model was trained) | '
        f'Red borders = statistically significant difference from within-basin (p<0.05)',
        fontsize=16, fontweight='bold', y=0.98
    )
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save figure
    if output_path.startswith('s3://'):
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
            temp_path = tmp.name
        plt.savefig(temp_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        s3_client = boto3.client('s3')
        bucket, key = output_path.replace('s3://', '').split('/', 1)
        s3_client.upload_file(temp_path, bucket, key)
        os.remove(temp_path)
        print(f"  ✓ Saved {metric} to S3")
    else:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved {metric} locally")


# ======================== Main Execution ========================

def generate_transfer_learning_analysis(
    basin_locations: List[Dict],
    base_prefix: str = "arrakis-data/floodOutputs",
    metrics: List[str] = METRICS,
    output_folder: str = None
):
    """Generate transfer learning analysis figures.
    
    Creates 7 PNG files (one per metric), each with 2x4 subplots.
    
    Args:
        basin_locations: List of basin location dicts
        base_prefix: S3 base prefix
        metrics: List of metrics to analyze
        output_folder: Optional output folder (defaults to transfer_learning_analysis/)
    """
    
    # Collect all data
    all_data = collect_all_data(basin_locations, base_prefix)
    
    # Determine output folder
    if output_folder is None:
        output_folder = (
            f"s3://climate-ai-data-science-datasets/{base_prefix}/"
            f"transfer_learning_analysis"
        )
    
    print("\n" + "=" * 70)
    print(f"GENERATING {len(metrics)} AGGREGATE FIGURES")
    print(f"Output: {output_folder}")
    print("=" * 70)
    
    # Create one figure per metric
    for metric in metrics:
        print(f"\nGenerating {metric} figure...")
        output_path = f"{output_folder}/transfer_learning_{metric}.png"
        
        create_aggregate_transfer_figure(
            all_data,
            metric,
            basin_locations,
            output_path
        )
    
    print("\n" + "=" * 70)
    print("✅ TRANSFER LEARNING ANALYSIS COMPLETE")
    print(f"Generated {len(metrics)} figures")
    print("=" * 70)


# ======================== Example Usage ========================

if __name__ == "__main__":
    
    # Run analysis for all metrics
    generate_transfer_learning_analysis(basin_locations)
    
    # Or specify custom output location
    # generate_transfer_learning_analysis(
    #     basin_locations,
    #     output_folder="s3://my-bucket/custom/path"
    # )
    
    # Or analyze specific metrics only
    # generate_transfer_learning_analysis(
    #     basin_locations,
    #     metrics=['MCC', 'F1']
    # )
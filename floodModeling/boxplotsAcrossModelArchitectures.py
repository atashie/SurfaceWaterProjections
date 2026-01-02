###############################################################
###############################################################
###############################################################
### Boxplots across model architectures - Enhanced Text
###############################################################
###############################################################
###############################################################



import os
import tempfile
from typing import List, Dict, Optional

import s3fs
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import numpy as np
from matplotlib.patches import PathPatch

# --------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------
BASE_S3_PREFIX = "climate-ai-data-science-datasets/arrakis-data/floodOutputs"

# Locations (EMSR ID, type, and lat/lon strings as used in S3 folder names)
locations = [
    {"emsr_id": "EMSR470", "flood_type": "Flash",    "lat": "10.7875",  "lon": "0.7479"},
    {"emsr_id": "EMSR444", "flood_type": "Riverine", "lat": "48.1747",  "lon": "23.3978"},
    {"emsr_id": "EMSR407", "flood_type": "Flash",    "lat": "53.6583",  "lon": "-1.1398"},
    {"emsr_id": "EMSR424", "flood_type": "Riverine", "lat": "-16.5646", "lon": "46.7857"},
    {"emsr_id": "EMSR438", "flood_type": "Coastal",  "lat": "-1.8038",  "lon": "33.3857"},
    {"emsr_id": "EMSR292", "flood_type": "Flash",    "lat": "40.9437",  "lon": "24.7669"},
    {"emsr_id": "EMSR629", "flood_type": "Riverine", "lat": "27.2981",  "lon": "68.1875"},
    {"emsr_id": "EMSR279", "flood_type": "Riverine", "lat": "41.9275",  "lon": "-1.3792"},
]

# Map logical model names to their subfolder names
MODEL_SUBFOLDERS = {
    "GLM":        "assessGLMModel",
    "RandomForest": "assessRandomForestModel",
    "XGBoost":    "assessXGBoostModel",
    "MLP-FiLM":   "assessMLPModel",
    "U-Net-FiLM": "assessUNetModel",
    "TabNet":     "assessTabNetModel",
}

# The order in which models should appear on x-axis
MODEL_ORDER = ["GLM", "RandomForest", "XGBoost", "MLP-FiLM", "U-Net-FiLM", "TabNet"]

# Metrics to summarize
METRICS = ["MCC", "F1", "Accuracy", "Precision", "Recall", "ROC_AUC", "PR_AUC"]

# Output folder
OUTPUT_DIR = "modelComparisons"

# Significance threshold for Mann–Whitney U test
ALPHA = 0.05

fs = s3fs.S3FileSystem(anon=False)


# --------------------------------------------------------------------
# Utilities: S3 reading and metrics loading
# --------------------------------------------------------------------
def _list_per_date_csvs_in_model_folder(
    base_prefix: str, lat: str, lon: str, model_subfolder: str
) -> List[str]:
    """
    List CSV files in a given model subfolder that contain 'per_date_' in the file name.
    Returns S3 keys (bucket/key without 's3://').
    """
    folder_name = f"trainingFolderFor_Lat_{lat}_Lon_{lon}"
    s3_prefix = f"{base_prefix}/{folder_name}/{model_subfolder}"
    pattern = f"{s3_prefix}/*.csv"
    all_csvs = fs.glob(pattern)
    per_date_csvs = [key for key in all_csvs if "per_date_" in os.path.basename(key)]
    return per_date_csvs


def _read_csv_from_s3(s3_key: str) -> pd.DataFrame:
    """
    Download CSV from S3 and read into a DataFrame.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        local_path = os.path.join(tmpdir, os.path.basename(s3_key))
        fs.download(s3_key, local_path)
        df = pd.read_csv(local_path)
    return df


def load_all_per_date_metrics(
    locations: List[Dict],
    base_s3_prefix: str = BASE_S3_PREFIX,
    model_subfolders: Dict[str, str] = MODEL_SUBFOLDERS,
) -> pd.DataFrame:
    """
    Read all per-date test metrics from all subwatersheds for all available models.

    Returns a single DataFrame with columns:
      - emsr_id, flood_type, lat, lon, model
      - plus all metric columns present in the per-date CSV files.
    """
    all_records = []

    for loc in locations:
        emsr_id = loc["emsr_id"]
        flood_type = loc["flood_type"]
        lat = loc["lat"]
        lon = loc["lon"]

        for model_name, subfolder in model_subfolders.items():
            csv_keys = _list_per_date_csvs_in_model_folder(
                base_prefix=base_s3_prefix,
                lat=lat,
                lon=lon,
                model_subfolder=subfolder,
            )
            if len(csv_keys) == 0:
                print(f"[WARN] No per_date_ CSVs for {emsr_id}, model {model_name}")
                continue

            dfs = []
            for key in csv_keys:
                df = _read_csv_from_s3(key)
                dfs.append(df)

            model_df = pd.concat(dfs, ignore_index=True)
            model_df["emsr_id"] = emsr_id
            model_df["flood_type"] = flood_type
            model_df["lat"] = float(lat)
            model_df["lon"] = float(lon)
            model_df["model"] = model_name

            all_records.append(model_df)

    if not all_records:
        raise RuntimeError("No per_date_ metrics found for any location/model.")

    all_metrics_df = pd.concat(all_records, ignore_index=True)
    return all_metrics_df


# --------------------------------------------------------------------
# Plot helpers
# --------------------------------------------------------------------

def _get_box_patches(ax):
    """
    Return the list of box patches corresponding to the categories on the x-axis.
    This is robust to differences in where Seaborn/Matplotlib store the boxes
    (ax.artists vs ax.patches).
    """
    # Prefer ax.artists if it contains PathPatch boxes
    boxes = [a for a in ax.artists if isinstance(a, PathPatch)]
    if not boxes:
        # Fall back to ax.patches and take PathPatch instances
        boxes = [p for p in ax.patches if isinstance(p, PathPatch)]

    # There may be more patches than categories (whiskers, medians, etc.),
    # so keep only the first n categories.
    n_cats = len(ax.get_xticklabels())
    return boxes[:n_cats]


def _highlight_boxes_vs_glm(ax, df_plot: pd.DataFrame, metric: str, alpha: float = ALPHA):
    """
    For each model (except GLM baseline), run Mann–Whitney U tests against GLM:
    - If significantly BETTER (alternative='greater', p < alpha): thick GOLD outline, blue fill
    - If significantly WORSE (alternative='less', p < alpha): thick RED outline, blue fill
    
    Boxes are mapped to models using the x-tick labels.
    """
    # Keep only rows for models we care about and valid metric values
    df_plot = df_plot[df_plot["model"].isin(MODEL_ORDER)].dropna(subset=[metric])

    # GLM (baseline) distribution
    glm_vals = df_plot.loc[df_plot["model"] == "GLM", metric].dropna()
    if glm_vals.empty:
        return

    # Get actual category labels and box patches from the axis
    xticklabels = [t.get_text() for t in ax.get_xticklabels()]
    box_patches = _get_box_patches(ax)

    if len(box_patches) != len(xticklabels):
        return

    # Map model label -> box patch
    model_to_box = {label: box for label, box in zip(xticklabels, box_patches)}

    # For each model, compute p-values vs GLM and apply highlighting
    for model in MODEL_ORDER:
        if model == "GLM":
            continue
        if model not in df_plot["model"].unique():
            continue
        if model not in model_to_box:
            continue

        model_vals = df_plot.loc[df_plot["model"] == model, metric].dropna()
        if model_vals.empty:
            continue

        box = model_to_box[model]
        
        # Test if significantly BETTER (greater)
        try:
            _, p_greater = mannwhitneyu(model_vals, glm_vals, alternative="greater")
        except ValueError:
            p_greater = 1.0
        
        # Test if significantly WORSE (less)
        try:
            _, p_less = mannwhitneyu(model_vals, glm_vals, alternative="less")
        except ValueError:
            p_less = 1.0

        # Apply styling based on results
        if p_greater < alpha:
            # Significantly better: thick gold outline, keep blue fill
            box.set_edgecolor("#FFD700")  # gold
            box.set_linewidth(3.0)
        elif p_less < alpha:
            # Significantly worse: thick red outline, keep blue fill
            box.set_edgecolor("#DC143C")  # crimson red
            box.set_linewidth(3.0)


def plot_per_location_per_metric(
    metrics_df: pd.DataFrame,
    output_dir: str = OUTPUT_DIR,
    metrics: List[str] = METRICS,
):
    """
    Step 2:
    For each location and each specified metric, generate a PNG with a boxplot of
    metric values per model. y in [0, 1], x: 6 models.
    Gold outline = significantly better than GLM
    Red outline = significantly worse than GLM
    """
    os.makedirs(output_dir, exist_ok=True)

    for loc in locations:
        emsr_id = loc["emsr_id"]
        df_loc = metrics_df[metrics_df["emsr_id"] == emsr_id].copy()

        if df_loc.empty:
            print(f"[WARN] No metrics for location {emsr_id}")
            continue

        for metric in metrics:
            if metric not in df_loc.columns:
                print(f"[WARN] Metric {metric} not found for {emsr_id}, skipping.")
                continue

            df_plot = df_loc[["model", metric]].dropna()
            if df_plot.empty:
                print(f"[WARN] No valid data for {emsr_id}, metric {metric}")
                continue

            # Restrict to models in MODEL_ORDER
            df_plot = df_plot[df_plot["model"].isin(MODEL_ORDER)]
            plt.figure(figsize=(8, 4))
            ax = sns.boxplot(
                data=df_plot,
                x="model",
                y=metric,
                order=MODEL_ORDER,
                showfliers=False,
            )
            ax.set_ylim(0.0, 1.0)
            # ✅ CHANGE 3: Remove "Model" label
            ax.set_xlabel("", fontsize=17)  # ✅ CHANGE 1: 70% bigger (10 * 1.7)
            # ✅ CHANGE 2: Remove metric name from y-axis
            ax.set_ylabel("", fontsize=20)  # ✅ CHANGE 1: 70% bigger (12 * 1.7)
            ax.set_title(f"{emsr_id} – per-date {metric}", fontsize=20)  # ✅ CHANGE 1: 70% bigger (12 * 1.7)
            ax.tick_params(axis='x', rotation=45, labelsize=15)  # ✅ CHANGE 1: 70% bigger (9 * 1.7)
            ax.tick_params(axis='y', labelsize=15)  # ✅ CHANGE 1: 70% bigger

            _highlight_boxes_vs_glm(ax, df_plot, metric, alpha=ALPHA)

            plt.tight_layout()
            out_path = os.path.join(output_dir, f"{emsr_id}_{metric}_per_date_by_model.png")
            plt.savefig(out_path, dpi=300, bbox_inches="tight")
            plt.close()


def plot_per_metric_all_locations_panels(
    metrics_df: pd.DataFrame,
    output_dir: str = OUTPUT_DIR,
    metrics: List[str] = METRICS,
):
    """
    Step 3:
    For each metric, generate a PNG figure with multiple subplots, one per location.
    Each subplot: y in [0, 1], x: 6 models, boxplots for that metric in that location.
    Now also applies highlighting (gold=better, red=worse) vs GLM in each panel.
    """
    os.makedirs(output_dir, exist_ok=True)

    n_locs = len(locations)
    ncols = 4
    nrows = (n_locs + ncols - 1) // ncols  # ceil

    for metric in metrics:
        plt.figure(figsize=(5 * ncols, 3 * nrows))

        for idx, loc in enumerate(locations):
            emsr_id = loc["emsr_id"]
            df_loc = metrics_df[metrics_df["emsr_id"] == emsr_id].copy()
            if metric not in df_loc.columns:
                ax = plt.subplot(nrows, ncols, idx + 1)
                ax.set_title(f"{emsr_id} (no {metric})", fontsize=20)  # ✅ CHANGE 1
                ax.set_xticks([])
                ax.set_yticks([])
                continue

            # Keep only rows with valid metric and models in MODEL_ORDER
            df_plot = df_loc[["model", metric]].dropna()
            df_plot = df_plot[df_plot["model"].isin(MODEL_ORDER)]

            ax = plt.subplot(nrows, ncols, idx + 1)

            if df_plot.empty:
                ax.set_title(f"{emsr_id} (no data)", fontsize=20)  # ✅ CHANGE 1
                ax.set_xticks([])
                ax.set_yticks([])
                continue

            sns.boxplot(
                data=df_plot,
                x="model",
                y=metric,
                order=MODEL_ORDER,
                showfliers=False,
                ax=ax,
            )
            ax.set_ylim(0.0, 1.0)
            # ✅ CHANGE 3: Remove "Model" label
            ax.set_xlabel("", fontsize=17)  # ✅ CHANGE 1
            # ✅ CHANGE 2: Remove metric name from y-axis (only show on leftmost)
            ax.set_ylabel("" if idx % ncols != 0 else "", fontsize=20)  # ✅ CHANGE 1 & 2
            ax.set_title(emsr_id, fontsize=20)  # ✅ CHANGE 1
            ax.tick_params(axis="x", rotation=45, labelsize=15)  # ✅ CHANGE 1
            ax.tick_params(axis="y", labelsize=15)  # ✅ CHANGE 1

            # Apply highlighting vs GLM within this panel
            _highlight_boxes_vs_glm(ax, df_plot, metric, alpha=ALPHA)

        plt.suptitle(f"Per-date {metric} by model for each location", 
                    y=0.98, fontsize=27)  # ✅ CHANGE 1: 70% bigger (16 * 1.7)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        out_path = os.path.join(
            output_dir,
            f"all_locations_{metric}_per_date_by_model.png",
        )
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()


def plot_aggregated_all_locations_and_save_stats(
    metrics_df: pd.DataFrame,
    output_dir: str = OUTPUT_DIR,
    metrics: List[str] = METRICS,
    alpha: float = ALPHA,
):
    """
    Aggregated plots + summary CSV.

    For each metric:
      - Aggregate per-date metrics across all locations.
      - Produce a boxplot: y in [0, 1], x: 6 models.
      - Run Mann–Whitney U tests vs GLM across all locations/dates.
      - Highlight significantly better (gold outline) and worse (red outline) models.

    Additionally:
      - Build a summary table with percentiles and p-values.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Percentile levels to compute
    quantiles = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # Row-oriented stats collector
    rows_dict: Dict[str, Dict[str, float]] = {}

    for metric in metrics:
        if metric not in metrics_df.columns:
            print(f"[WARN] Metric {metric} not found in DataFrame columns; skipping.")
            continue

        df_plot = metrics_df[["model", metric]].dropna()
        df_plot = df_plot[df_plot["model"].isin(MODEL_ORDER)]
        if df_plot.empty:
            print(f"[WARN] No data for metric {metric} across all locations")
            continue

        # 1) Make aggregated boxplot (with highlighting)
        plt.figure(figsize=(8, 4))
        ax = sns.boxplot(
            data=df_plot,
            x="model",
            y=metric,
            order=MODEL_ORDER,
            showfliers=False,
        )
        ax.set_ylim(0.0, 1.0)
        # ✅ CHANGE 3: Remove "Model" label
        ax.set_xlabel("", fontsize=17)  # ✅ CHANGE 1
        # ✅ CHANGE 2: Remove metric name from y-axis
        ax.set_ylabel("", fontsize=20)  # ✅ CHANGE 1 & 2
        ax.set_title(f"All locations – per-date {metric} by model", 
                    fontsize=20)  # ✅ CHANGE 1
        ax.tick_params(axis='x', rotation=45, labelsize=15)  # ✅ CHANGE 1
        ax.tick_params(axis='y', labelsize=15)  # ✅ CHANGE 1

        _highlight_boxes_vs_glm(ax, df_plot, metric, alpha=alpha)

        plt.tight_layout()
        out_path = os.path.join(
            output_dir,
            f"all_locations_aggregated_{metric}_per_date_by_model.png",
        )
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()

        # 2) Percentiles and p-values vs GLM
        glm_vals = df_plot.loc[df_plot["model"] == "GLM", metric].dropna()

        for model in MODEL_ORDER:
            model_vals = df_plot.loc[df_plot["model"] == model, metric].dropna()
            if model_vals.empty:
                for q in quantiles:
                    row_label = f"{model}_p{int(q * 100)}"
                    if row_label not in rows_dict:
                        rows_dict[row_label] = {}
                    rows_dict[row_label][metric] = np.nan
                    pcol = f"{metric}_pvalue"
                    rows_dict[row_label][pcol] = np.nan
                continue

            # Percentiles
            for q in quantiles:
                row_label = f"{model}_p{int(q * 100)}"
                if row_label not in rows_dict:
                    rows_dict[row_label] = {}
                q_val = np.quantile(model_vals, q)
                rows_dict[row_label][metric] = float(q_val)

            # Mann–Whitney U p-value vs GLM
            pcol = f"{metric}_pvalue"
            if model == "GLM" or glm_vals.empty:
                pval = np.nan
            else:
                try:
                    _, pval = mannwhitneyu(model_vals, glm_vals, alternative="greater")
                except ValueError:
                    pval = np.nan

            for q in quantiles:
                row_label = f"{model}_p{int(q * 100)}"
                if row_label not in rows_dict:
                    rows_dict[row_label] = {}
                rows_dict[row_label][pcol] = float(pval) if pval is not np.nan else np.nan

    # 3) Build and save summary CSV
    summary_df = pd.DataFrame.from_dict(rows_dict, orient="index")
    summary_df.index.name = "model_percentile"

    ordered_cols = []
    for metric in metrics:
        if metric in summary_df.columns:
            ordered_cols.append(metric)
    for metric in metrics:
        pcol = f"{metric}_pvalue"
        if pcol in summary_df.columns:
            ordered_cols.append(pcol)

    summary_df = summary_df.reindex(columns=ordered_cols)

    csv_path = os.path.join(
        output_dir,
        "aggregated_metrics_percentiles_and_mannwhitney_pvalues.csv",
    )
    summary_df.to_csv(csv_path, index=True)
    print(f"[INFO] Saved aggregated summary CSV to: {csv_path}")


def plot_all_metrics_summary(
    metrics_df: pd.DataFrame,
    output_dir: str = OUTPUT_DIR,
    metrics: List[str] = METRICS,
    alpha: float = ALPHA,
):
    """
    NEW: Create a single summary figure with all 7 metrics in a 2x4 grid.
    Bottom-right subplot contains legend.
    
    Each subplot shows aggregated boxplots (across all locations) for one metric,
    with gold/red outlining for significant differences vs GLM baseline.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        
        if metric not in metrics_df.columns:
            ax.set_title(f"{metric} (no data)", fontsize=20)  # ✅ CHANGE 1
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        
        df_plot = metrics_df[["model", metric]].dropna()
        df_plot = df_plot[df_plot["model"].isin(MODEL_ORDER)]
        
        if df_plot.empty:
            ax.set_title(f"{metric} (no data)", fontsize=20)  # ✅ CHANGE 1
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        
        # Create boxplot
        sns.boxplot(
            data=df_plot,
            x="model",
            y=metric,
            order=MODEL_ORDER,
            showfliers=False,
            ax=ax,
        )
        
        ax.set_ylim(0.0, 1.0)
        # ✅ CHANGE 3: Remove "Model" label
        ax.set_xlabel("", fontsize=17)  # ✅ CHANGE 1
        # ✅ CHANGE 2: Remove metric name from y-axis
        ax.set_ylabel("", fontsize=20)  # ✅ CHANGE 1 & 2
        ax.set_title(f"{metric}", fontsize=20, fontweight='bold')  # ✅ CHANGE 1
        ax.tick_params(axis='x', rotation=45, labelsize=15)  # ✅ CHANGE 1
        ax.tick_params(axis='y', labelsize=15)  # ✅ CHANGE 1
        ax.grid(True, alpha=0.3, axis='y')
        
        # Apply highlighting
        _highlight_boxes_vs_glm(ax, df_plot, metric, alpha=alpha)
    
    # Hide the 8th subplot (bottom right) and add legend
    axes[7].axis('off')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='lightblue', edgecolor='black', linewidth=1, label='No significant difference'),
        Patch(facecolor='lightblue', edgecolor='#FFD700', linewidth=3, label='Significantly better (p<0.05)'),
        Patch(facecolor='lightblue', edgecolor='#DC143C', linewidth=3, label='Significantly worse (p<0.05)'),
    ]
    axes[7].legend(handles=legend_elements, loc='center', 
                   fontsize=20,  # ✅ CHANGE 1: 70% bigger (12 * 1.7)
                   frameon=True, 
                   title='vs. GLM Baseline', 
                   title_fontsize=20)  # ✅ CHANGE 1
    
    # Overall title
    fig.suptitle('Model Architecture Comparison\nAggregated Performance Across All Locations',
                fontsize=27, fontweight='bold', y=0.98)  # ✅ CHANGE 1: 70% bigger (16 * 1.7)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    out_path = os.path.join(output_dir, "all_metrics_summary.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[INFO] Saved all-metrics summary figure to: {out_path}")


# --------------------------------------------------------------------
# Main driver
# --------------------------------------------------------------------
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1) Load all per-date metrics
    all_metrics = load_all_per_date_metrics(locations)

    # 2) Per-location, per-metric boxplots (with significance vs GLM)
    plot_per_location_per_metric(
        all_metrics,
        output_dir=OUTPUT_DIR,
        metrics=METRICS,
    )

    # 3) Per-metric, all-locations multi-panel figures
    plot_per_metric_all_locations_panels(
        all_metrics,
        output_dir=OUTPUT_DIR,
        metrics=METRICS,
    )

    # 4) Aggregated across locations, per-metric boxplots + CSV
    plot_aggregated_all_locations_and_save_stats(
        all_metrics,
        output_dir=OUTPUT_DIR,
        metrics=METRICS,
        alpha=ALPHA,
    )

    # 5) NEW: All-metrics summary figure (2x4 grid)
    plot_all_metrics_summary(
        all_metrics,
        output_dir=OUTPUT_DIR,
        metrics=METRICS,
        alpha=ALPHA,
    )

# train_flood_model.py

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List
import pandas as pd
import boto3

from flood_data_container import FloodDataContainer, WatershedMetadata
from feature_processor import TabularFeatureProcessor
from decision_tree_pipeline import DecisionTreeFloodModel, LeaveOneWatershedOutCV
from visualization import PerformanceVisualizer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_watershed_metadata(training_folder: str) -> List[WatershedMetadata]:
    """Load watershed metadata from training folder."""
    
    # Load list of training locations
    locations_file = Path(training_folder) / 'list_of_training_locations.csv'
    
    if locations_file.exists():
        df = pd.read_csv(locations_file)
    else:
        # Try S3
        s3 = boto3.client('s3')
        bucket = training_folder.split('/')[2]
        key = '/'.join(training_folder.split('/')[3:]) + '/list_of_training_locations.csv'
        
        obj = s3.get_object(Bucket=bucket, Key=key)
        df = pd.read_csv(obj['Body'])
    
    # Create metadata objects
    metadata_list = []
    for _, row in df.iterrows():
        metadata = WatershedMetadata(
            hybas_id=row['hybas_id'],
            latitude=row['latitude'],
            longitude=row['longitude'],
            lat_str=row['lat_str'],
            lon_str=row['lon_str'],
            is_initial=row.get('is_initial', False),
            area_km2=row.get('area_km2', 0)
        )
        metadata_list.append(metadata)
    
    return metadata_list


def prepare_watershed_data(watershed_metadata: WatershedMetadata,
                          training_folder: str,
                          processor: TabularFeatureProcessor) -> pd.DataFrame:
    """Prepare data for a single watershed."""
    
    logger.info(f"Processing watershed {watershed_metadata.hybas_id}")
    
    # Create container
    container = FloodDataContainer(watershed_metadata, training_folder)
    
    # Load data
    container.load_static_raster()
    container.load_dynamic_features()
    container.load_flood_targets()
    
    # Process features
    df = processor.process_watershed(container)
    
    # Add metadata
    df['is_initial'] = watershed_metadata.is_initial
    
    return df


def main(args):
    """Main training pipeline."""
    
    logger.info("Starting flood model training pipeline")
    logger.info(f"Training folder: {args.training_folder}")
    
    # Load watershed metadata
    watershed_metadata_list = load_watershed_metadata(args.training_folder)
    logger.info(f"Found {len(watershed_metadata_list)} watersheds")
    
    # Initialize feature processor
    processor = TabularFeatureProcessor(
        lag_days=args.lag_days,
        window_sizes=[3, 5, 7],
        include_temporal_features=True,
        include_spatial_context=True
    )
    
    # Process all watersheds
    watershed_data = {}
    for metadata in watershed_metadata_list:
        try:
            df = prepare_watershed_data(metadata, args.training_folder, processor)
            if not df.empty:
                watershed_data[str(metadata.hybas_id)] = df
                logger.info(f"Processed {len(df)} samples for watershed {metadata.hybas_id}")
        except Exception as e:
            logger.error(f"Failed to process watershed {metadata.hybas_id}: {e}")
    
    if not watershed_data:
        logger.error("No watersheds successfully processed")
        return
    
    logger.info(f"Successfully processed {len(watershed_data)} watersheds")
    
    # Run cross-validation
    cv = LeaveOneWatershedOutCV(watershed_data)
    
    model_params = {
        'max_depth': args.max_depth,
        'min_samples_split': args.min_samples_split,
        'min_samples_leaf': args.min_samples_leaf,
        'random_state': args.random_state
    }
    
    logger.info("Running cross-validation...")
    cv_results = cv.run_cv(
        model_class=DecisionTreeFloodModel,
        model_params=model_params
    )
    
    # Save model and results
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save model
    cv_results['final_model'].save_model(output_dir / 'model')
    
    # Save CV results
    with open(output_dir / 'cv_results.json', 'w') as f:
        # Remove non-serializable objects
        results_to_save = {
            'summary': cv_results['summary'],
            'test_metrics': cv_results['test_results']['metrics']
        }
        json.dump(results_to_save, f, indent=2, default=str)
    
    # Create visualizations
    logger.info("Creating visualizations...")
    visualizer = PerformanceVisualizer(output_dir / 'visualizations')
    visualizer.plot_cv_results(cv_results)
    visualizer.create_interactive_metrics_dashboard(cv_results)
    visualizer.save_metrics_report(cv_results)
    
    # Print summary
    print("\n" + "=" * 60)
    print("TRAINING COMPLETE")
    print("=" * 60)
    print(f"\nValidation F1 Score: {cv_results['summary']['validation']['val_f1_score_mean']:.4f} "
          f"(± {cv_results['summary']['validation']['val_f1_score_std']:.4f})")
    print(f"Test F1 Score: {cv_results['test_results']['metrics']['f1_score']:.4f}")
    print(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train flood prediction model')
    
    parser.add_argument(
        '--training_folder',
        type=str,
        required=True,
        help='Path to training folder (S3 or local)'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default='./flood_model_output',
        help='Directory to save outputs'
    )
    
    parser.add_argument(
        '--lag_days',
        type=int,
        default=3,
        help='Number of lag days for dynamic features'
    )
    
    parser.add_argument(
        '--max_depth',
        type=int,
        default=None,
        help='Maximum depth of decision tree'
    )
    
    parser.add_argument(
        '--min_samples_split',
        type=int,
        default=100,
        help='Minimum samples to split node'
    )
    
    parser.add_argument(
        '--min_samples_leaf',
        type=int,
        default=50,
        help='Minimum samples in leaf'
    )
    
    parser.add_argument(
        '--random_state',
        type=int,
        default=42,
        help='Random seed'
    )
    
    args = parser.parse_args()
    main(args)

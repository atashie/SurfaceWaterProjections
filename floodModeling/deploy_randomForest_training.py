# deploy_randomForest_training.py
"""
Deploy Random Forest hyperparameter tuning to SageMaker CPU instance.
Uses standard instances only (no spot instances).
Random Forest is CPU-based, so we use compute-optimized instances.
"""

import sagemaker
from sagemaker.sklearn import SKLearn
from sagemaker import get_execution_role
import os
import shutil


def prepare_training_code():
    """Prepare training code directory with necessary files."""
    if os.path.exists('training_code'):
        shutil.rmtree('training_code')
    os.makedirs('training_code', exist_ok=True)
    
    # Copy the Random Forest module
    if os.path.exists('randomForestFloods.py'):
        shutil.copy('randomForestFloods.py', 'training_code/')
        print("Copied randomForestFloods.py to training_code/")
    else:
        raise FileNotFoundError("randomForestFloods.py not found")
    
    # Create minimal requirements (no GDAL for training)
    minimal_requirements = [
        "numpy>=1.23.0",
        "pandas>=1.5.0",
        "optuna>=3.5.0",
        "pyarrow>=10.0.0",
        "boto3>=1.26.0",
        "matplotlib>=3.6.0",
        "psutil>=5.8.0",
        "joblib>=1.0.0"
    ]
    
    with open('training_code/requirements.txt', 'w') as f:
        f.write("\n".join(minimal_requirements) + "\n")
    print("Created training_code/requirements.txt (no GDAL)")
    
    # Create entry point wrapper - EXACTLY like MLP
    entry_content = """
import os
os.environ['RUN_TRAINING'] = '1'
import sys
from randomForestFloods import main

if __name__ == '__main__':
    main()
"""
    
    with open('training_code/train_rf_optuna.py', 'w') as f:
        f.write(entry_content)
    print("Created training_code/train_rf_optuna.py")
    
    print("Training code prepared successfully")


def deploy_training_job(parquet_uri, experiment_name='rf-flood-optuna'):
    """Deploy Random Forest training with Optuna to SageMaker."""
    
    # Initialize SageMaker session
    sagemaker_session = sagemaker.Session()
    role = get_execution_role()
    bucket = sagemaker_session.default_bucket()
    
    # Prepare training code
    prepare_training_code()
    
    # Create SKLearn estimator for CPU training
    estimator = SKLearn(
        entry_point='train_rf_optuna.py',
        source_dir='training_code',
        role=role,
        instance_type='ml.m5.12xlarge',  # CPU-optimized: 48 vCPUs, 192GB RAM
        instance_count=1,
        framework_version='1.2-1',
        py_version='py3',
        hyperparameters={
            'parquet-uri': parquet_uri,
            'n-trials': 50
        },
        max_run=43200,  # 12 hours max
        use_spot_instances=False,
        volume_size=100,
        output_path=f's3://{bucket}/randomforest-outputs/{experiment_name}',
        code_location=f's3://{bucket}/randomforest-code/{experiment_name}',
        environment={'PYTHONUNBUFFERED': '1'}
    )
    
    # Start training job
    print("Starting SageMaker training job...")
    estimator.fit(wait=False)
    
    # Get the training job name
    training_job_name = estimator.latest_training_job.name
    
    print("\n" + "="*70)
    print("RANDOM FOREST TRAINING JOB STARTED SUCCESSFULLY")
    print("="*70)
    print(f"Job name: {training_job_name}")
    print(f"Instance type: ml.m5.12xlarge (48 vCPUs, 192GB RAM)")
    print(f"Monitor: https://console.aws.amazon.com/sagemaker/home?region={sagemaker_session.boto_region_name}#/jobs/{training_job_name}")
    print(f"Outputs: s3://{bucket}/randomforest-outputs/{experiment_name}/{training_job_name}/output/")
    
    return estimator, training_job_name


def attach_to_training_job(training_job_name):
    """Attach to an existing training job."""
    estimator = SKLearn.attach(training_job_name)
    return estimator


if __name__ == "__main__":
    # Configuration
    PARQUET_URI = "s3://climate-ai-data-science-datasets/your-flood-data.parquet"  # Update this
    EXPERIMENT_NAME = "rf-flood-optuna-001"
    
    # Deploy the training job
    estimator, job_name = deploy_training_job(PARQUET_URI, EXPERIMENT_NAME)
    
    print("\nTo attach later:")
    print(f"from deploy_randomForest_training import attach_to_training_job")
    print(f"estimator = attach_to_training_job('{job_name}')")
    print("estimator.wait()")

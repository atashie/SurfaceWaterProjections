"""
deploy_mlpFilm_training.py
Deploy MLP+FiLM hyperparameter tuning to SageMaker GPU instance.
Uses standard instances only (no spot instances).
"""

import sagemaker
from sagemaker.pytorch import PyTorch
from sagemaker import get_execution_role
import os
import shutil


def prepare_training_code():
    """Prepare training code directory with necessary files."""
    if os.path.exists('training_code'):
        shutil.rmtree('training_code')
    os.makedirs('training_code', exist_ok=True)
    
    # Copy the MLP+FiLM module
    if os.path.exists('mlpFilmFloods.py'):
        shutil.copy('mlpFilmFloods.py', 'training_code/')
        print("Copied mlpFilmFloods.py to training_code/")
    else:
        raise FileNotFoundError("mlpFilmFloods.py not found")
    
    # Create minimal requirements (no GDAL for training, like TabNet)
    minimal_requirements = [
        "torch>=2.0.0",
        "numpy>=1.23.0",
        "pandas>=1.5.0",
        "scikit-learn>=1.2.0",
        "optuna>=3.5.0",
        "pyarrow>=10.0.0",
        "boto3>=1.26.0",
        "matplotlib>=3.6.0",
        "psutil>=5.8.0",
        "joblib>=1.0.0"
    ]
    
    # Filter out geospatial dependencies if copying from root requirements
    root_req_path = 'requirements.txt'
    filtered = []
    if os.path.exists(root_req_path):
        with open(root_req_path, 'r') as f:
            for line in f:
                ln = line.strip()
                if not ln or ln.startswith('#'):
                    continue
                lower = ln.lower()
                # Skip packages that pull native GDAL toolchain
                if any(bad in lower for bad in ["gdal", "rasterio", "fiona", "geopandas", "osgeo"]):
                    print(f"Skipping geospatial dependency: {ln}")
                    continue
                filtered.append(ln)
    
    # Ensure required deps are present
    present = {pkg.split('=')[0].split('<')[0].split('>')[0] for pkg in filtered}
    for pkg in minimal_requirements:
        base = pkg.split('=')[0].split('<')[0].split('>')[0]
        if base not in present:
            filtered.append(pkg)
    
    # Write training-only requirements
    with open('training_code/requirements.txt', 'w') as f:
        f.write("\n".join(filtered) + "\n")
    print("Created training_code/requirements.txt (no GDAL)")
    
    # Create entry point wrapper
    entry_content = """
import os
os.environ['RUN_TRAINING'] = '1'
import sys
from mlpFilmFloods import main

if __name__ == '__main__':
    main()
"""
    
    with open('training_code/train_mlp_optuna.py', 'w') as f:
        f.write(entry_content)
    print("Created training_code/train_mlp_optuna.py")
    
    print("Training code prepared successfully")


def deploy_training_job(parquet_uri, experiment_name='mlp-film-optuna'):
    """Deploy MLP+FiLM training with Optuna to SageMaker."""
    
    # Initialize SageMaker session
    sagemaker_session = sagemaker.Session()
    role = get_execution_role()
    bucket = sagemaker_session.default_bucket()
    
    # Prepare training code
    prepare_training_code()
    
    # Create PyTorch estimator for GPU training
    estimator = PyTorch(
        entry_point='train_mlp_optuna.py',
        source_dir='training_code',
        role=role,
        instance_type='ml.p3.2xlarge',  # Same GPU as TabNet
        instance_count=1,
        framework_version='2.1.0',
        py_version='py310',
        hyperparameters={
            'parquet-uri': parquet_uri,
            'n-trials': 50,  # Number of Optuna trials
            'epochs': 20,    # Epochs per trial (reduced for tuning)
            'patience': 5,   # Early stopping
            'use-gpu': True
        },
        max_run=86400,  # 24 hours max
        use_spot_instances=False,  # Standard instances only
        volume_size=100,  # 100GB storage
        output_path=f's3://{bucket}/mlp-film-outputs/{experiment_name}',
        checkpoint_s3_uri=f's3://{bucket}/mlp-film-checkpoints/{experiment_name}',
        checkpoint_local_path='/opt/ml/checkpoints',
        environment={'PYTHONUNBUFFERED': '1'}
    )
    
    # Start training job
    print("Starting SageMaker training job...")
    estimator.fit(wait=False)
    
    # Get the training job name
    training_job_name = estimator.latest_training_job.name
    
    print("\n" + "="*70)
    print("MLP+FiLM TRAINING JOB STARTED SUCCESSFULLY")
    print("="*70)
    print(f"Job name: {training_job_name}")
    print(f"Instance type: ml.p3.2xlarge")
    print(f"Monitor: https://console.aws.amazon.com/sagemaker/home?region={sagemaker_session.boto_region_name}#/jobs/{training_job_name}")
    print(f"Outputs: s3://{bucket}/mlp-film-outputs/{experiment_name}/{training_job_name}/output/")
    print(f"Checkpoints: s3://{bucket}/mlp-film-checkpoints/{experiment_name}/")
    
    return estimator, training_job_name


def attach_to_training_job(training_job_name):
    """Attach to an existing training job."""
    estimator = PyTorch.attach(training_job_name)
    return estimator


if __name__ == "__main__":
    # Configuration
    PARQUET_URI = "s3://climate-ai-data-science-datasets/your-flood-data.parquet"  # Update this
    EXPERIMENT_NAME = "mlp-film-flood-optuna-001"
    
    # Deploy the training job
    estimator, job_name = deploy_training_job(PARQUET_URI, EXPERIMENT_NAME)
    
    print("\nTo attach later:")
    print(f"from deploy_mlpFilm_training import attach_to_training_job")
    print(f"estimator = attach_to_training_job('{job_name}')")
    print("estimator.wait()")
    
    print("\nTo download results:")
    print("import boto3")
    print("s3 = boto3.client('s3')")
    print(f"bucket = '{estimator.output_path.split('/')[2]}'")
    print(f"s3.download_file(bucket, 'mlp-film-outputs/{EXPERIMENT_NAME}/{job_name}/output/model.tar.gz', 'model.tar.gz')")

"""
deploy_tabnet_training.py
Deploy TabNet hyperparameter tuning to SageMaker GPU instance.
Uses standard instances only (no spot instances).
"""

import sagemaker
from sagemaker.pytorch import PyTorch
from sagemaker import get_execution_role
import os
import shutil

def prepare_training_code():
    """Prepare training code directory with all necessary files."""
    if os.path.exists('training_code'):
        shutil.rmtree('training_code')
    os.makedirs('training_code', exist_ok=True)

    # Always copy the training module
    if os.path.exists('tabnetFloods.py'):
        shutil.copy('tabnetFloods.py', 'training_code/')
        print("Copied tabnetFloods.py to training_code/")
    else:
        print("Warning: tabnetFloods.py not found")

    # Create a minimal, training-only requirements.txt to avoid GDAL build failures
    # Add only the packages needed for tuning/training (no GDAL/rasterio/fiona/geopandas)
    minimal_requirements = [
        "pytorch-tabnet>=4.1.0",
        "optuna>=3.5.0",
        "pyarrow>=10.0.0",
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "scikit-learn>=1.2.0",
        "boto3>=1.26.0",
        # matplotlib is used only if you run visualization inside training; keep lightweight
        "matplotlib>=3.6.0"
    ]

    # If a root requirements.txt exists, filter out heavy geospatial deps (GDAL/rasterio/fiona/geopandas)
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
                    print(f"Skipping geospatial dependency in training environment: {ln}")
                    continue
                filtered.append(ln)

    # Ensure required training deps are present
    present = {pkg.split('=')[0].split('<')[0].split('>')[0] for pkg in filtered}
    for pkg in minimal_requirements:
        base = pkg.split('=')[0].split('<')[0].split('>')[0]
        if base not in present:
            filtered.append(pkg)

    # Write the training-only requirements.txt
    train_req_path = os.path.join('training_code', 'requirements.txt')
    with open(train_req_path, 'w') as f:
        f.write("\n".join(filtered) + "\n")
    print("Created training_code/requirements.txt (training-only, no GDAL)")

    # Create the entry point wrapper for Optuna tuning
    entry_path = os.path.join('training_code', 'train_tabnet_optuna.py')
    with open(entry_path, 'w') as f:
        f.write(
            "import sys\n"
            "from tabnetFloods import main\n"
            "if __name__ == '__main__':\n"
            "    main()\n"
        )
    print("Created training_code/train_tabnet_optuna.py")

    print("Training code prepared successfully")
    
    
def deploy_training_job(parquet_uri, experiment_name='tabnet-optuna'):
    """Deploy TabNet training with Optuna to SageMaker using standard instances."""
    
    # Initialize SageMaker session
    sagemaker_session = sagemaker.Session()
    role = get_execution_role()
    bucket = sagemaker_session.default_bucket()
    
    # Prepare training code
    prepare_training_code()
    
    # Create PyTorch estimator for GPU training
    estimator = PyTorch(
        entry_point='train_tabnet_optuna.py',   # matches the file we just created
        source_dir='training_code',
        role=role,
        instance_type='ml.g5.4xlarge',
        instance_count=1,
        framework_version='2.1.0',
        py_version='py310',
        hyperparameters={
            'parquet-uri': parquet_uri,
            'n-trials': 50,
            'epochs': 30,
            'patience': 5,
            'use-gpu': True
        },
        max_run=86400,
        use_spot_instances=False,
        volume_size=100,
        output_path=f's3://{bucket}/tabnet-outputs/{experiment_name}',
        checkpoint_s3_uri=f's3://{bucket}/tabnet-checkpoints/{experiment_name}',
        checkpoint_local_path='/opt/ml/checkpoints',
        environment={'PYTHONUNBUFFERED': '1'}
    )    
    
    # Start training job
    print("Starting SageMaker training job...")
    estimator.fit(wait=False)
    
    # Get the training job name
    training_job_name = estimator.latest_training_job.name
    
    print("\n" + "="*70)
    print("TRAINING JOB STARTED SUCCESSFULLY")
    print("="*70)
    print(f"Job name: {training_job_name}")
    print(f"Instance type: ml.g5.4xlarge (standard instance)")
    print(f"Monitor at: https://console.aws.amazon.com/sagemaker/home?region={sagemaker_session.boto_region_name}#/jobs/{training_job_name}")
    print(f"Outputs will be saved to: s3://{bucket}/tabnet-outputs/{experiment_name}/{training_job_name}/output/")
    print(f"Checkpoints will be saved to: s3://{bucket}/tabnet-checkpoints/{experiment_name}/")
    
    return estimator, training_job_name


def attach_to_training_job(training_job_name):
    """Attach to an existing training job to monitor or download results."""
    estimator = PyTorch.attach(training_job_name)
    return estimator


if __name__ == "__main__":
    # Configuration
    PARQUET_URI = "s3://climate-ai-data-science-datasets/your-flood-data.parquet"  # Update this
    EXPERIMENT_NAME = "tabnet-optuna-flood-001"
    
    # Deploy the training job
    estimator, job_name = deploy_training_job(PARQUET_URI, EXPERIMENT_NAME)
    
    print("\nTo attach to this job later and wait for completion:")
    print(f"from deploy_tabnet_training import attach_to_training_job")
    print(f"estimator = attach_to_training_job('{job_name}')")
    print("estimator.wait()")
    
    print("\nTo download results after completion:")
    print("import boto3")
    print("s3 = boto3.client('s3')")
    print(f"bucket = '{estimator.output_path.split('/')[2]}'")
    print(f"s3.download_file(bucket, 'tabnet-outputs/{EXPERIMENT_NAME}/{job_name}/output/model.tar.gz', 'model.tar.gz')")

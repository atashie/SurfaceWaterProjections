# deploy_unetFilm_training.py
"""
Deploy U-Net+FiLM hyperparameter tuning to SageMaker GPU instance.
Uses ml.g5.2xlarge (24GB A10G) for better cost/performance than P3.
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
    
    # Copy the U-Net+FiLM module
    if os.path.exists('unetFilmFloods.py'):
        shutil.copy('unetFilmFloods.py', 'training_code/')
        print("Copied unetFilmFloods.py to training_code/")
    else:
        raise FileNotFoundError("unetFilmFloods.py not found")
    
    # Training-only requirements (rasterio instead of GDAL)
    minimal_requirements = [
        "numpy>=1.23.0",
        "pandas>=1.5.0",
        "scikit-learn>=1.2.0",
        "optuna>=3.5.0",
        "boto3>=1.26.0",
        "matplotlib>=3.6.0",
        "psutil>=5.8.0",
        "joblib>=1.0.0",
        "rasterio[s3]>=1.3.0"  # Includes S3 extras
    ]
    
    with open('training_code/requirements.txt', 'w') as f:
        f.write("\n".join(minimal_requirements) + "\n")
    print("Created training_code/requirements.txt (rasterio[s3], no system GDAL)")
    
    # Create entry point wrapper
    entry_content = """
import os
os.environ['RUN_TRAINING'] = '1'
import sys
from unetFilmFloods import main

if __name__ == '__main__':
    main()
"""
    
    with open('training_code/train_unet_optuna.py', 'w') as f:
        f.write(entry_content)
    print("Created training_code/train_unet_optuna.py")
    
    print("Training code prepared successfully")


def deploy_training_job(training_folder, experiment_name='unet-film-optuna'):
    """Deploy U-Net+FiLM training with Optuna to SageMaker."""
    
    # Initialize SageMaker session
    sagemaker_session = sagemaker.Session()
    role = get_execution_role()
    bucket = sagemaker_session.default_bucket()
    
    # Prepare training code
    prepare_training_code()
    
    # Create PyTorch estimator for GPU training
    # ml.g5.2xlarge: 1 A10G (24GB), 8 vCPUs, 32GB RAM - better than P3 for this workload
    estimator = PyTorch(
        entry_point='train_unet_optuna.py',
        source_dir='training_code',
        role=role,
        instance_type='ml.g5.2xlarge',  # A10G GPU with 24GB VRAM
        instance_count=1,
        framework_version='2.1.0',
        py_version='py310',
        hyperparameters={
            'training-folder': training_folder,
            'n-trials': 30,      # Reduced from 50 to fit 24hr budget
            'epochs': 15,        # Epochs per trial (reduced for tuning)
            'patience': 5,       # Early stopping
            'use-gpu': True,
            'num-workers': 8     # Good for 8 vCPU instance
        },
        max_run=86400,  # 24 hours max
        use_spot_instances=False,  # Standard instances only
        volume_size=100,  # 100GB storage
        output_path=f's3://{bucket}/unet-film-outputs/{experiment_name}',
        checkpoint_s3_uri=f's3://{bucket}/unet-film-checkpoints/{experiment_name}',
        checkpoint_local_path='/opt/ml/checkpoints',
        environment={'PYTHONUNBUFFERED': '1'}
    )
    
    # Start training job
    print("Starting SageMaker training job...")
    estimator.fit(wait=False)
    
    # Get the training job name
    training_job_name = estimator.latest_training_job.name
    
    print("\n" + "="*70)
    print("U-NET+FiLM TRAINING JOB STARTED SUCCESSFULLY")
    print("="*70)
    print(f"Job name: {training_job_name}")
    print(f"Instance type: ml.g5.2xlarge (1 A10G GPU, 24GB VRAM, 8 vCPUs, 32GB RAM)")
    print(f"Monitor: https://console.aws.amazon.com/sagemaker/home?region={sagemaker_session.boto_region_name}#/jobs/{training_job_name}")
    print(f"Outputs: s3://{bucket}/unet-film-outputs/{experiment_name}/{training_job_name}/output/")
    print(f"Checkpoints: s3://{bucket}/unet-film-checkpoints/{experiment_name}/")
    
    return estimator, training_job_name


def attach_to_training_job(training_job_name):
    """Attach to an existing training job."""
    estimator = PyTorch.attach(training_job_name)
    return estimator


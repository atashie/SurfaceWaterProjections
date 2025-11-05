"""
deploy_glm_training.py
Deploy GLM hyperparameter tuning to SageMaker CPU instance.
Uses standard instances only (no spot instances).
GLM is CPU-based and lightweight, so we use smaller instances.
"""

import os
import shutil
import textwrap
import sagemaker
from sagemaker.sklearn import SKLearn
from sagemaker import get_execution_role

def prepare_training_code():
    """Prepare training code directory with necessary files."""
    if os.path.exists('training_code'):
        shutil.rmtree('training_code')
    os.makedirs('training_code', exist_ok=True)

    # Copy GLM module
    if not os.path.exists('glmFloods.py'):
        raise FileNotFoundError("glmFloods.py not found")
    shutil.copy('glmFloods.py', 'training_code/')
    print("Copied glmFloods.py to training_code/")

    # Minimal training-only requirements (no GDAL, no sklearn)
    requirements = textwrap.dedent("""\
        numpy>=1.23.0
        pandas>=1.5.0
        optuna>=3.5.0
        pyarrow>=10.0.0
        boto3>=1.26.0
        matplotlib>=3.6.0
        psutil>=5.8.0
        joblib>=1.0.0
    """)
    with open('training_code/requirements.txt', 'w') as f:
        f.write(requirements.strip() + "\n")
    print("Created training_code/requirements.txt (training-only; no GDAL, no sklearn)")

    # Robust entry wrapper: parse known args, ignore unknowns, rebuild sys.argv
    entry_content = textwrap.dedent("""\
        import os
        import sys
        import argparse
        from glmFloods import main as glm_main

        def _parse_known():
            p = argparse.ArgumentParser(add_help=False)
            p.add_argument('--parquet-uri', type=str, required=True)
            p.add_argument('--n-trials', type=int, default=50)
            # SageMaker standard flags we accept (for dir wiring)
            p.add_argument('--model-dir', type=str, default=os.environ.get('SM_MODEL_DIR', '/opt/ml/model'))
            p.add_argument('--output-data-dir', type=str, default=os.environ.get('SM_OUTPUT_DATA_DIR', '/opt/ml/output'))
            p.add_argument('--output-dir', type=str, default=os.environ.get('SM_OUTPUT_DATA_DIR', '/opt/ml/output'))
            args, unknown = p.parse_known_args()
            if unknown:
                print(f"[Info] Wrapper ignoring unknown container args: {unknown}", flush=True)
            return args

        if __name__ == '__main__':
            args = _parse_known()
            out_dir = args.output_data_dir if getattr(args, 'output_data_dir', None) else args.output_dir
            model_dir = args.model_dir

            # Rebuild argv strictly for glmFloods.main()
            sys.argv = [
                'glmFloods.py',
                '--parquet-uri', args.parquet_uri,
                '--n-trials', str(args.n_trials),
                '--model-dir', model_dir,
                '--output-dir', out_dir,
            ]
            glm_main()
    """)
    with open('training_code/train_glm_optuna.py', 'w') as f:
        f.write(entry_content)
    print("Created training_code/train_glm_optuna.py")

    print("Training code prepared successfully")


def deploy_training_job(parquet_uri, experiment_name='glm-flood-optuna-001'):
    """Deploy GLM training with Optuna to SageMaker."""
    sagemaker_session = sagemaker.Session()
    role = get_execution_role()
    bucket = sagemaker_session.default_bucket()

    prepare_training_code()

    estimator = SKLearn(
        entry_point='train_glm_optuna.py',
        source_dir='training_code',
        role=role,
        instance_type='ml.m5.4xlarge',
        instance_count=1,
        framework_version='1.2-1',   # SKLearn container version
        py_version='py3',
        hyperparameters={
            'parquet-uri': parquet_uri,
            'n-trials': 50
        },
        max_run=6 * 3600,
        use_spot_instances=False,
        volume_size=50,
        output_path=f's3://{bucket}/glm-outputs/{experiment_name}',
        code_location=f's3://{bucket}/glm-code/{experiment_name}',
        environment={'PYTHONUNBUFFERED': '1'}
    )

    print("Starting SageMaker training job...")
    estimator.fit(wait=False)

    job_name = estimator.latest_training_job.name
    region = sagemaker_session.boto_region_name

    print("\n" + "="*70)
    print("GLM TRAINING JOB STARTED SUCCESSFULLY")
    print("="*70)
    print(f"Job name: {job_name}")
    print(f"Instance type: ml.m5.4xlarge (16 vCPUs, 64 GB RAM)")
    print(f"Monitor: https://console.aws.amazon.com/sagemaker/home?region={region}#/jobs/{job_name}")
    print(f"Outputs: s3://{bucket}/glm-outputs/{experiment_name}/{job_name}/output/")

    return estimator, job_name



def attach_to_training_job(training_job_name: str):
    """Attach to an existing GLM training job to monitor/wait."""
    return SKLearn.attach(training_job_name)


if __name__ == "__main__":
    # Configuration
    PARQUET_URI = "s3://climate-ai-data-science-datasets/your-flood-data.parquet"  # Update this
    EXPERIMENT_NAME = "glm-flood-optuna-001"
    
    # Deploy the training job
    estimator, job_name = deploy_training_job(PARQUET_URI, EXPERIMENT_NAME)
    
    print("\nTo attach later:")
    print(f"from deploy_glm_training import attach_to_training_job")
    print(f"estimator = attach_to_training_job('{job_name}')")
    print("estimator.wait()")
    
    print("\nTo download results:")
    print("import boto3")
    print("s3 = boto3.client('s3')")
    print(f"bucket = '{estimator.output_path.split('/')[2]}'")
    print(f"s3.download_file(bucket, 'glm-outputs/{EXPERIMENT_NAME}/{job_name}/output/model.tar.gz', 'model.tar.gz')")

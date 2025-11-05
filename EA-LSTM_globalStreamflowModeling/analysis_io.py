# analysis_io.py
import os
import shutil
from typing import List, Optional

def build_local_analysis_bundle(model_dir: str, output_data_dir: str, bundle_dir: str) -> None:
    """
    Merge analysis outputs from:
      - {output_data_dir}/analysis/
      - {model_dir}/output_files/analysis/
    into a single local folder (bundle_dir) for inspection.
    """
    os.makedirs(bundle_dir, exist_ok=True)
    sources = [
        os.path.join(output_data_dir or '', 'analysis'),
        os.path.join(model_dir or '', 'output_files', 'analysis'),
    ]
    for src in sources:
        if not src or not os.path.exists(src):
            continue
        for entry in os.listdir(src):
            s = os.path.join(src, entry)
            d = os.path.join(bundle_dir, entry)
            if os.path.isdir(s):
                if os.path.exists(d):
                    shutil.rmtree(d)
                shutil.copytree(s, d)
            else:
                os.makedirs(os.path.dirname(d), exist_ok=True)
                shutil.copy2(s, d)
    print(f"Bundled analysis outputs into: {bundle_dir}")


def download_analysis_bundle_from_s3(
    bucket: str,
    s3_prefixes: List[str],
    local_dir: str,
    aws_profile: Optional[str] = None
) -> None:
    """
    Download all objects under the provided S3 prefixes into local_dir (merging).
    s3_prefixes examples:
      - 's3://my-bucket/path/to/model_dir/output_files/analysis/'
      - 's3://my-bucket/path/to/output_data_dir/analysis/'
    """
    import boto3

    # Support optional profiles for local use
    if aws_profile:
        import botocore
        session = boto3.Session(profile_name=aws_profile)
        s3 = session.client('s3', config=botocore.config.Config(signature_version='s3v4'))
    else:
        s3 = boto3.client('s3')

    os.makedirs(local_dir, exist_ok=True)

    def _split_s3_uri(uri: str):
        # Accept either 's3://bucket/key' or 'bucket/key'
        if uri.startswith('s3://'):
            _, _, rest = uri.partition('s3://')
            b, _, k = rest.partition('/')
            return b, k
        else:
            b, _, k = uri.partition('/')
            return b, k

    for prefix in s3_prefixes:
        bkt, key_prefix = _split_s3_uri(prefix)
        paginator = s3.get_paginator('list_objects_v2')
        pages = paginator.paginate(Bucket=bkt, Prefix=key_prefix)
        found = False
        for page in pages:
            for obj in page.get('Contents', []):
                found = True
                key = obj['Key']
                rel = os.path.relpath(key, key_prefix) if key.startswith(key_prefix) else os.path.basename(key)
                local_path = os.path.join(local_dir, rel)
                os.makedirs(os.path.dirname(local_path), exist_ok=True)
                s3.download_file(bkt, key, local_path)
        if not found:
            print(f"Warning: No objects under s3://{bkt}/{key_prefix}")

    print(f"Downloaded analysis bundle from S3 into: {local_dir}")

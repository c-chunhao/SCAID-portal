"""Explicit local inputs and relative registered asset resolution."""
from pathlib import Path
import json
import os

def required_path(name):
    value = os.environ.get(name, '').strip()
    if not value:
        raise ValueError(f'Missing required local-input variable: {name}')
    return Path(value).expanduser()

def registered_datasets():
    registry = required_path('SCAID_ASSET_REGISTRY')
    datasets = json.loads(registry.read_text())['datasets']
    base = os.environ.get('SCAID_DATA_H5_ROOT', '').strip()
    for dataset in datasets:
        for method, value in dataset['h5_paths'].items():
            path = Path(value)
            if not path.is_absolute():
                if not base:
                    raise ValueError('Relative asset registry requires SCAID_DATA_H5_ROOT')
                if '..' in path.parts:
                    raise ValueError('Registered relative asset path must not traverse parents')
                path = Path(base).expanduser() / path
            dataset['h5_paths'][method] = str(path)
    return datasets

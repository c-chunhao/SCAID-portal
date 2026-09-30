"""Create a synthetic, nonclinical fixture for testing a clean SCAID install."""
import argparse
import json
from pathlib import Path

import h5py
from PIL import Image

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output_directory')
args = parser.parse_args()
root = Path(args.output_directory).resolve()
figure = root / 'figures/DEMO/DEMO1/PBMC/GeneUmap/Umap_DEMO.png'
figure.parent.mkdir(parents=True, exist_ok=True)
Image.new('RGB', (64, 48), (240, 235, 230)).save(figure)
h5 = root / 'h5/RNA/DEMO_DEMO1_PBMC.h5'
h5.parent.mkdir(parents=True, exist_ok=True)
with h5py.File(h5, 'w') as file:
    file.attrs['fixture_type'] = 'synthetic nonclinical demonstration; no patient data'
    file.create_dataset('illustration_values', data=[0, 1])
manifest = {
    'schema_version': 1,
    'release_id': 'synthetic-demo-v1',
    'source_url': 'https://example.org/synthetic-demo',
    'datasets': [{
        'catalog': {'disease_name': 'Synthetic demonstration (not a clinical dataset)', 'abbreviation': 'DEMO', 'dataset_source': 'Demo', 'dataset_id': 'DEMO1', 'tissue': 'PBMC', 'sample_size': 0, 'cell_count': 0, 'batch': '20260930'},
        'figures': [{'path': 'DEMO/DEMO1/PBMC/GeneUmap/Umap_DEMO.png'}],
        'h5_files': [{'path': 'RNA/DEMO_DEMO1_PBMC.h5', 'method': 'RNA'}],
    }],
}
(root / 'public-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
print(root / 'public-manifest.json')

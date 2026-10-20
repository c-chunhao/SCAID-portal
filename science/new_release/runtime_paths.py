"""Portable local inputs/outputs; no private data, accounts or credentials."""
import json
import os
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parent
REPOSITORY_ROOT = CODE_ROOT.parents[1]
def configured(name, default):
    return Path(os.environ.get(name, str(default))).expanduser().resolve()
WORK = configured('SCAID_SCIENCE_WORK_ROOT', CODE_ROOT / 'work')
INPUT_ROOT = configured('SCAID_SCIENCE_INPUT_ROOT', CODE_ROOT / 'inputs')
PLAN_PATH = configured('SCAID_RELEASE_PLAN', CODE_ROOT / 'examples/release_plan.example.json')
SEURAT_ROOT = configured('SCAID_SEURAT_INPUT_ROOT', INPUT_ROOT / '00.data')
ANNOTATION_ROOT = configured('SCAID_ANNOTATION_INPUT_ROOT', INPUT_ROOT / '01.ann')
RAW_ROOT = configured('SCAID_GSE279086_RAW_ROOT', INPUT_ROOT / 'GSE279086')
EXPORT_ROOT = configured('SCAID_EXPORT_ROOT', WORK / 'exports')
FIGURE_ROOT = configured('SCAID_FIGURE_ROOT', WORK / 'figures')
INDEX_ROOT = configured('SCAID_GENE_INDEX_ROOT', WORK / 'gene_cache/indexes')
GENE_MANIFEST = configured('SCAID_GENE_CACHE_MANIFEST', WORK / 'gene_cache/lazy_gene_plots.json')
BACKEND_ROOT = configured('SCAID_BACKEND_SOURCE', REPOSITORY_ROOT / 'backend')
def load_plan():
    value = json.loads(PLAN_PATH.read_text(encoding='utf-8'))
    if value.get('schema_version') != 1 or not isinstance(value.get('objects'), list):
        raise ValueError('Invalid release plan; provide a reviewed schema_version=1 object list.')
    for entry in value['objects']:
        parts = Path(entry['object_id']).parts
        if len(parts) != 3 or any(x in ('.', '..') for x in parts) or Path(entry['object_id']).is_absolute():
            raise ValueError('Each object_id must be condition/accession/tissue without path traversal.')
    return value

# R subprocesses receive the same reviewed configuration as Python callers.
for name, value in {'SCAID_RELEASE_PLAN': PLAN_PATH, 'SCAID_SCIENCE_WORK_ROOT': WORK,
                    'SCAID_SEURAT_INPUT_ROOT': SEURAT_ROOT, 'SCAID_EXPORT_ROOT': EXPORT_ROOT}.items():
    os.environ.setdefault(name, str(value))

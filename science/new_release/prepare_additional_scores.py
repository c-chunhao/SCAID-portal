"""Validate and index the three previously withheld pathway image families."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from collections import Counter
import hashlib
import json
from PIL import Image

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
SOURCE = ANNOTATION_ROOT
PLAN = PLAN_PATH
METHODS = ('AUCell', 'UCell', 'singscore')


def inspect_file(item):
    object_id, method, path = item
    before = path.stat()
    if path.is_symlink() or path.resolve() != path or before.st_size < 100:
        raise ValueError(f'Unsafe or empty image: {path}')
    with Image.open(path) as image:
        image.load()
        width, height = image.size
        if image.format != 'PNG' or min(width, height) < 100:
            raise ValueError(f'Invalid image: {path}')
    digest = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'Source changed during validation: {path}')
    return {
        'object_id': object_id, 'module': 'kegg_score', 'method': method,
        'name': path.name, 'relative_path': f'{object_id}/KEGGScore/{method}',
        'physical_path': str(path), 'width': width, 'height': height,
        'bytes': before.st_size, 'mtime_ns': before.st_mtime_ns, 'sha256': digest,
        'page_number': 1,
    }


def main():
    objects = [obj for obj in json.loads(PLAN.read_text())['objects']
               if obj['tier'] == 'main_figures_only']
    assert len(objects) == 3
    jobs, scopes = [], []
    for obj in objects:
        key = obj['object_id']
        families = {}
        for method in METHODS:
            root = SOURCE / key / 'KEGGScore' / method
            files = sorted(root.glob('*.png'))
            if not files:
                raise ValueError(f'Missing score family {root}')
            families[method] = {path.name for path in files}
            jobs.extend((key, method, path) for path in files)
        if not all(value == families[METHODS[0]] for value in families.values()):
            raise ValueError(f'Method pathway names differ for {key}')
        scopes.append({'object_id': key, 'catalog': obj['catalog'],
                       'pathway_count': len(families[METHODS[0]]),
                       'approved_root': str(SOURCE / key / 'KEGGScore')})
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(inspect_file, jobs))
    out = WORK / 'validation/additional_scores'
    out.mkdir(parents=True, exist_ok=True)
    index = out / 'images.jsonl'
    index.write_text(''.join(json.dumps(row, ensure_ascii=False) + '\n' for row in rows))
    report = {'passed': True, 'objects': scopes, 'image_count': len(rows),
              'total_bytes': sum(row['bytes'] for row in rows),
              'folder_counts': dict(Counter(row['relative_path'] for row in rows)),
              'index_sha256': hashlib.file_digest(index.open('rb'), 'sha256').hexdigest(),
              'validation': 'Every PNG fully decoded; SHA256 and stable source stat recorded; identical pathway names across all three methods.'}
    (out / 'receipt.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(report, ensure_ascii=False))


if __name__ == '__main__':
    main()

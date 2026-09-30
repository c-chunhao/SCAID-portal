"""Conservative mappings between catalog labels and existing asset names.

Only spelling/format aliases live here. Unknown accessions must be corrected in
catalog metadata or explicitly reviewed; they never fall back to another dataset.
"""
import itertools
from functools import lru_cache
import re
from pathlib import Path

from django.conf import settings
from django.db.models import Q
from django.http import Http404


TISSUE_ALIASES = {
    'epidermal': ('Skin_Epidermal',),
    'dermal': ('Skin_Dermal',),
    'skin biopsy': ('SkinBiopsy',),
    '皮肤': ('Skin',),
    '肺': ('Lung',),
    '外周血': ('PBMC',),
    '水泡': ('Blister',),
    '眼前房液（房水）': ('aqueous_humour',),
}
DISEASE_FOLDERS = {
    'AIH': ('Autoimmune_hepatitis',),
    'BP': ('bullous_pemphigoid',),
    'SSc': ('systemic_sclerosis',),
    'PN': ('prurigo_nodularis',),
    'LS': ('Localized_scleroderma',),
}
H5_DISEASE_NAMES = {
    'systemic_sclerosis': ('SSs',),
    'type_1_diabetes_mellitus': ('T1DM',),
    'IBD_Colitis_inflamed': ('IBD_CI', 'IBD'),
    'IBD_Crohns_Disease': ('IBD_CD', 'IBD'),
    'IBD_Ulcerative_Colitis': ('IBD_UC', 'IBD'),
}


def component(value):
    return (value or '').strip().replace('/', '_').replace('\\', '_').replace(' ', '_')


def tissue_names(value):
    return tuple(dict.fromkeys((component(value), *TISSUE_ALIASES.get((value or '').strip().lower(), ()))))


def dataset_names(dataset_id, source=None):
    names = [component(dataset_id)]
    if (source or '').upper() == 'E-MTAB' and not names[0].upper().startswith('E-MTAB-'):
        names.append('E-MTAB-' + names[0])
    return tuple(names)


def disease_names(record):
    names = [component(record.abbreviation)]
    names.extend(DISEASE_FOLDERS.get(record.abbreviation, ()))
    # These are stored catalog labels, not inferred disease associations.
    names.extend((component(record.disease_full_name), component(record.disease_name)))
    return tuple(dict.fromkeys(name for name in names if name))


def path_sequence_condition(field, *groups):
    """Match complete adjacent directory components, treating input literally.

    Explicit separators prevent AS matching another disease, PBMC matching
    PBMC_SFMC, and GSE1 matching GSE10. re.escape prevents regex injection.
    """
    if not groups or any(not group or not all(group) for group in groups):
        return Q(pk__in=[])
    separator = r'[/\\]'
    patterns = ['(?:' + '|'.join(re.escape(value) for value in group) + ')' for group in groups]
    return Q(**{field + '__iregex': r'(?:^|[/\\])' + separator.join(patterns) + separator})


def record_image_condition(record):
    if not record.dataset_id or not record.tissue:
        return Q(pk__in=[])
    return path_sequence_condition(
        'image', disease_names(record), dataset_names(record.dataset_id, record.dataset_source), tissue_names(record.tissue)
    )


def h5_names(record):
    if not record.dataset_id or not record.tissue:
        return ()
    diseases = [component(record.abbreviation)] if record.abbreviation else []
    for name in disease_names(record):
        diseases.extend(H5_DISEASE_NAMES.get(name, ()))
    if not diseases:
        return ()
    datasets = list(dataset_names(record.dataset_id, record.dataset_source))
    datasets.extend(name.replace('-', '.') for name in datasets if name.upper().startswith('E-MTAB-'))
    stems = {
        '_'.join(parts)
        for parts in itertools.product(set(diseases), set(datasets), tissue_names(record.tissue))
    }
    # The SV files additionally include their actual scoring method in the name.
    suffixes = ('', '_AUCell', '_RNA', '_UCell', '_singscore')
    return tuple(stem + suffix + '.h5' for stem in stems for suffix in suffixes)


def resolved_h5_path(record):
    """Allow only existing H5 files whose resolved path is inside a configured root."""
    try:
        path = Path(record.path)
        if not path.is_absolute():
            raise Http404('File unavailable.')
        path = path.resolve(strict=True)
        roots = [Path(root).resolve(strict=True) for root in getattr(settings, 'H5_DOWNLOAD_ROOTS', ())]
        if path.suffix.lower() not in {'.h5', '.hdf5', '.h5ad'}:
            raise Http404('File unavailable.')
        if not path.is_file() or not any(path.is_relative_to(root) for root in roots):
            raise Http404('File unavailable.')
        return path
    except (OSError, RuntimeError, ValueError, TypeError):
        raise Http404('File unavailable.') from None


# The project-supplied manifest retains the original labels and source-row
# provenance. Figure components and complete H5 paths provide explicit asset
# relationships, including distinct IBD subtypes with the same accession.
@lru_cache(maxsize=1)
def load_manifest():
    import json
    manifest_path = Path(__file__).with_name('data') / 'dataset_assets.json'
    entries = json.loads(manifest_path.read_text(encoding='utf-8'))['datasets']
    figure_root = Path(settings.SCAID_FIGURE_ROOT)
    h5_roots = getattr(settings, 'H5_DOWNLOAD_ROOTS', ())
    h5_root = Path(h5_roots[0]) if h5_roots else Path('/srv/scaid/data/Data.H5')
    for entry in entries:
        if not Path(entry['figure_path']).is_absolute():
            entry['figure_path'] = str(figure_root / entry['figure_path'])
        entry['figure_path_aliases'] = [str(figure_root / path) if not Path(path).is_absolute() else path for path in entry.get('figure_path_aliases', [])]
        entry['h5_paths'] = {method: str(h5_root / path) if path and not Path(path).is_absolute() else path for method, path in entry['h5_paths'].items()}
    return entries


def label_key(value):
    return re.sub(r'[\s_\-/]+', '', value or '').casefold()


def manifest_entry(record):
    record_diseases = {label_key(value) for value in disease_names(record)}
    record_datasets = {label_key(value) for value in dataset_names(record.dataset_id, record.dataset_source)}
    record_tissues = {label_key(value) for value in tissue_names(record.tissue)}
    matches = []
    for entry in load_manifest():
        aliases = entry.get('catalog_aliases', [])
        if any(all(getattr(record, field, None) == value for field, value in alias.items()) for alias in aliases if alias):
            matches.append(entry)
            continue
        figure = Path(entry['figure_path']).parts[-3:]
        diseases = {label_key(entry['disease_name']), label_key(entry['abbreviation']), label_key(entry.get('canonical_abbreviation')), label_key(figure[0])}
        datasets = {label_key(entry.get('asset_dataset_id', entry['dataset_id'])), label_key(figure[1])}
        if ':' in entry['dataset_id'] and not entry['dataset_id'].startswith('http'):
            datasets.add(label_key(entry['dataset_id'].split(':', 1)[1]))
        tissues = {label_key(entry['tissue']), label_key(figure[2])}
        if record_diseases & diseases and record_datasets & datasets and record_tissues & tissues:
            matches.append(entry)
    return matches[0] if len(matches) == 1 else None


def manifest_folders(record):
    """'disease/dataset/tissue' directory prefixes of a manifest-linked record.

    Ingestion stores this prefix (plus the figure family) in relative_path, so
    exact or prefix matches on the indexed column replace regex scans of the
    1.1M-row image table. Returns [] when the record has no manifest entry.
    """
    entry = manifest_entry(record)
    if entry is None:
        return []
    return ['/'.join(Path(path).parts[-3:]) for path in [entry['figure_path'], *entry.get('figure_path_aliases', [])]]


def folder_condition(folders, *suffixes, prefix=False):
    values = [f'{base}/{suffix}' if suffix else base for base in folders for suffix in (suffixes or ('',))]
    if prefix:
        condition = Q(pk__in=[])
        for value in values:
            condition |= Q(relative_path__startswith=value.rstrip('/') + '/')
        return condition
    return Q(relative_path__in=values)


def manifest_image_condition(record):
    entry = manifest_entry(record)
    if entry is None:
        return record_image_condition(record)
    condition = Q(pk__in=[])
    for path in [entry['figure_path'], *entry.get('figure_path_aliases', [])]:
        parts = Path(path).parts[-3:]
        condition |= path_sequence_condition('image', *((part,) for part in parts))
    return condition


def manifest_h5_condition(record):
    entry = manifest_entry(record)
    if entry is not None:
        return Q(path__in=[value for value in entry['h5_paths'].values() if value])
    condition = Q(pk__in=[])
    for name in h5_names(record):
        condition |= Q(name__iexact=name)
    return condition


@lru_cache(maxsize=1)
def figure_manifest_folders():
    by_folder = {}
    for entry in load_manifest():
        for path in [entry['figure_path'], *entry.get('figure_path_aliases', [])]:
            folder = '/'.join(label_key(value) for value in Path(path).parts[-3:])
            if folder in by_folder and by_folder[folder] != entry:
                by_folder[folder] = None
            else:
                by_folder[folder] = entry
    return by_folder


def figure_manifest_entry(record):
    parts = str(record.relative_path or '').replace('\\', '/').strip('/').split('/')
    if len(parts) < 3:
        return None
    return figure_manifest_folders().get('/'.join(label_key(value) for value in parts[:3]))


def effective_abbreviation(record):
    """Public condition code; legacy storage labels remain untouched."""
    entry = manifest_entry(record)
    if entry:
        return entry.get('canonical_abbreviation') or entry['abbreviation']
    return record.abbreviation or record.disease_full_name or record.disease_name or 'Unclassified'


def disease_filter_labels(record, entry=None):
    entry = entry or manifest_entry(record)
    if not entry:
        return {label_key(effective_abbreviation(record))}
    canonical = entry.get('canonical_abbreviation') or entry['abbreviation']
    labels = {label_key(canonical)}
    labels.update(label_key(value) for value in entry.get('search_abbreviations', []))
    # A reclassified object's historical directory code must not masquerade as
    # the current diagnosis: AS means AS, while E-MTAB-8207 is explicitly PsA.
    if not entry.get('canonical_abbreviation'):
        labels.update((label_key(record.abbreviation), label_key(entry['abbreviation']), label_key(Path(entry['figure_path']).parts[-3])))
    return labels


def matching_catalog_records(records, abbreviation=None, dataset_id=None, tissue=None):
    result = []
    for record in records:
        entry = manifest_entry(record)
        disease_labels = disease_filter_labels(record, entry)
        dataset_labels = {label_key(value) for value in dataset_names(record.dataset_id, record.dataset_source)}
        tissue_labels = {label_key(value) for value in tissue_names(record.tissue)}
        if entry:
            figure = Path(entry['figure_path']).parts[-3:]
            dataset_labels.update((label_key(entry.get('asset_dataset_id', entry['dataset_id'])), label_key(figure[1])))
            tissue_labels.update((label_key(entry['tissue']), label_key(figure[2])))
        if abbreviation and label_key(abbreviation) not in disease_labels:
            continue
        if dataset_id and label_key(dataset_id) not in dataset_labels:
            continue
        if tissue and not ({label_key(value) for value in tissue_names(tissue)} & tissue_labels):
            continue
        result.append(record)
    return result


def canonical_query_ids(records, query):
    """An exact public condition code uses the same scope as its filter."""
    records = list(records)
    code = label_key(query)
    known = set().union(*(disease_filter_labels(record) for record in records)) if records else set()
    if code not in known:
        return None
    return [record.pk for record in records if code in disease_filter_labels(record)]


def manifest_search_ids(records, query):
    """Search public catalog labels and actual basenames, never private paths."""
    query = query.strip().casefold()
    matches = []
    for record in records:
        entry = manifest_entry(record)
        if entry is None:
            continue
        labels = [entry.get('display_name') or entry['disease_name'], effective_abbreviation(record), *entry.get('search_abbreviations', [])]
        labels.extend(Path(path).name for path in entry['h5_paths'].values() if path)
        if any(query in value.casefold() for value in labels):
            matches.append(record.pk)
    return matches

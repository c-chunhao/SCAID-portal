"""Register reviewed public catalogue metadata and existing assets, atomically.

Default mode is dry-run. This command accepts aggregate catalogue rows only,
never patient-level tables, counts matrices embedded in JSON, or private IDs.
"""
import hashlib
import json
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from PIL import Image

from tisch_api.models import CellData, H5File, PDFImage, GeneKeggEnum

CATALOG_FIELDS = {'disease_name', 'disease_full_name', 'abbreviation', 'dataset_source', 'dataset_id', 'tissue', 'sample_size', 'cell_count', 'batch'}


def within_root(root, value):
    relative = Path(value)
    if relative.is_absolute() or '..' in relative.parts:
        raise CommandError('Asset paths must be relative to the supplied root, without traversal.')
    try:
        path = (root / relative).resolve(strict=True)
    except OSError as exc:
        raise CommandError('A referenced asset is unavailable.') from exc
    if not path.is_file() or not path.is_relative_to(root):
        raise CommandError('An asset resolves outside the supplied root.')
    return path


def reject_unknown(row, allowed, description):
    if not isinstance(row, dict) or set(row) - allowed:
        raise CommandError(f'{description} has unsupported fields; patient-level metadata is not accepted.')


class Command(BaseCommand):
    help = 'Validate a reviewed public asset manifest; dry-run unless --apply. No patient tables are imported.'

    def add_arguments(self, parser):
        parser.add_argument('manifest')
        parser.add_argument('--figure-root', required=True)
        parser.add_argument('--h5-root', required=True)
        parser.add_argument('--apply', action='store_true', help='Apply the complete validated plan in one database transaction.')
        parser.add_argument('--report', help='Optional JSON audit report; contains aggregate identifiers and manifest digest only.')

    def handle(self, *args, **options):
        try:
            raw = Path(options['manifest']).read_bytes()
            manifest = json.loads(raw)
        except (OSError, ValueError) as exc:
            raise CommandError('Manifest must be a readable JSON file.') from exc
        reject_unknown(manifest, {'schema_version', 'release_id', 'source_url', 'datasets'}, 'Manifest')
        if manifest.get('schema_version') != 1 or not isinstance(manifest.get('datasets'), list) or not manifest['datasets']:
            raise CommandError('Expected schema_version 1 and a nonempty datasets array.')
        if not isinstance(manifest.get('release_id'), str) or not manifest['release_id'].strip():
            raise CommandError('A release_id is required for provenance.')
        if not str(manifest.get('source_url', '')).startswith('https://'):
            raise CommandError('A public HTTPS source_url is required.')
        roots = {}
        for kind in ('figure', 'h5'):
            try:
                roots[kind] = Path(options[kind + '_root']).resolve(strict=True)
            except OSError as exc:
                raise CommandError('Both asset roots must exist.') from exc
            if not roots[kind].is_dir():
                raise CommandError('Asset roots must be directories.')
        plan = []
        identities = set()
        for dataset in manifest['datasets']:
            reject_unknown(dataset, {'catalog', 'figures', 'h5_files'}, 'Dataset')
            catalog = dataset.get('catalog')
            reject_unknown(catalog, CATALOG_FIELDS, 'Catalogue')
            for field in ('abbreviation', 'disease_name', 'dataset_source', 'dataset_id', 'tissue'):
                if not isinstance(catalog.get(field), str) or not catalog[field].strip():
                    raise CommandError(f'Catalogue {field} is required.')
            for field in ('sample_size', 'cell_count'):
                value = catalog.get(field)
                if value is not None and (not isinstance(value, int) or isinstance(value, bool) or not 0 <= value <= 2147483647):
                    raise CommandError(f'{field} must be a nonnegative aggregate count.')
            for field, value in catalog.items():
                model_field = CellData._meta.get_field(field)
                if value is not None and model_field.max_length and (not isinstance(value, str) or len(value) > model_field.max_length):
                    raise CommandError(f'Catalogue {field} exceeds its schema limit.')
            key = {field: catalog[field] for field in ('abbreviation', 'dataset_id', 'tissue')}
            self.add_plan(plan, identities, CellData, key, catalog)
            for key in ('figures', 'h5_files'):
                if not isinstance(dataset.get(key, []), list):
                    raise CommandError(f'{key} must be an array.')
            for figure in dataset.get('figures', []):
                reject_unknown(figure, {'path', 'name', 'page_number', 'is_pdf'}, 'Figure')
                path = within_root(roots['figure'], figure.get('path', ''))
                if path.suffix.lower() not in {'.png', '.jpg', '.jpeg'}:
                    raise CommandError('Only prepared PNG/JPEG figures can be registered.')
                try:
                    with Image.open(path) as image:
                        width, height = image.size
                        image.verify()
                except (OSError, ValueError) as exc:
                    raise CommandError('A figure is not a readable valid image.') from exc
                if not isinstance(figure.get('is_pdf', False), bool):
                    raise CommandError('Figure is_pdf must be a boolean.')
                page = figure.get('page_number', 1)
                if not isinstance(page, int) or isinstance(page, bool) or page < 1:
                    raise CommandError('Figure page_number must be a positive integer.')
                values = {'name': figure.get('name') or path.name, 'path': str(path), 'relative_path': path.parent.relative_to(roots['figure']).as_posix(), 'image': str(path), 'page_number': page, 'image_width': str(width), 'image_height': str(height), 'is_pdf': str(bool(figure.get('is_pdf', False)))}
                if len(values['name']) > 255 or len(str(path)) > 512:
                    raise CommandError('Figure name or path exceeds schema limits.')
                self.add_plan(plan, identities, PDFImage, {'image': str(path), 'page_number': page}, values)
                parts = path.parent.relative_to(roots['figure']).parts
                if 'GeneUmap' in parts and path.stem.startswith('Umap_'):
                    self.add_plan(plan, identities, GeneKeggEnum, {'tp': 'gene', 'name': path.stem[5:]}, {}, allow_repeat=True)
                elif 'KEGGScore' in parts:
                    self.add_plan(plan, identities, GeneKeggEnum, {'tp': 'kegg', 'name': path.stem}, {}, allow_repeat=True)
            methods = set()
            for asset in dataset.get('h5_files', []):
                reject_unknown(asset, {'path', 'method'}, 'H5 asset')
                method = asset.get('method')
                if method not in {'RNA', 'AUCell', 'UCell', 'singscore'} or method in methods:
                    raise CommandError('H5 methods must be distinct RNA/AUCell/UCell/singscore values.')
                methods.add(method)
                path = within_root(roots['h5'], asset.get('path', ''))
                with path.open('rb') as stream:
                    signature = stream.read(8)
                if path.suffix.lower() not in {'.h5', '.hdf5', '.h5ad'} or signature != b'\x89HDF\r\n\x1a\n':
                    raise CommandError('An H5 asset has an invalid extension or HDF5 signature.')
                if len(str(path)) > 1024 or len(path.name) > 255:
                    raise CommandError('H5 path exceeds schema limits.')
                self.add_plan(plan, identities, H5File, {'path': str(path)}, {'name': path.name, 'path': str(path), 'relative_path': method, 'file_size': path.stat().st_size})
        report = {'schema_version': 1, 'release_id': manifest['release_id'], 'source_url': manifest['source_url'], 'manifest_sha256': hashlib.sha256(raw).hexdigest(), 'mode': 'apply' if options['apply'] else 'dry-run', 'patient_tables_imported': False, 'changes': []}
        # Lock existing rows, reject ambiguous duplicates, validate every row and
        # perform the complete import in a single transaction. Repeated imports
        # update by natural identity rather than appending duplicate rows.
        with transaction.atomic():
            for model, key, values in plan:
                matches = list(model.objects.select_for_update().filter(**key)[:2])
                if len(matches) > 1:
                    raise CommandError(f'Ambiguous existing {model.__name__} identity; resolve duplicates first.')
                obj = matches[0] if matches else model(**key)
                changed = {field: value for field, value in values.items() if str(getattr(obj, field)) != str(value)}
                action = 'update' if matches and changed else 'unchanged' if matches else 'create'
                for field, value in values.items():
                    setattr(obj, field, value)
                try:
                    obj.full_clean(validate_unique=False, validate_constraints=False)
                except Exception as exc:
                    raise CommandError(f'{model.__name__} row does not match the public schema: {exc}') from exc
                if options['apply'] and action != 'unchanged':
                    obj.save()
                public_key = key if model in (CellData, GeneKeggEnum) else {'name': values['name'], 'page_number': values.get('page_number'), 'method': values.get('relative_path') if model is H5File else None}
                report['changes'].append({'model': model.__name__, 'identity': public_key, 'action': action})
            if not options['apply']:
                transaction.set_rollback(True)
        report['counts'] = {action: sum(change['action'] == action for change in report['changes']) for action in ('create', 'update', 'unchanged')}
        serialized = json.dumps(report, ensure_ascii=False, indent=2)
        self.stdout.write(serialized)
        if options['report']:
            Path(options['report']).write_text(serialized + '\n', encoding='utf-8')

    @staticmethod
    def add_plan(plan, identities, model, key, values, allow_repeat=False):
        identity = (model.__name__, tuple(sorted(key.items())))
        if identity in identities:
            if allow_repeat:
                return
            raise CommandError('Manifest repeats an asset or catalogue identity.')
        identities.add(identity)
        plan.append((model, key, values))

"""Read-only, exact-object access to validated cell-type marker tables.

Set SCAID_MARKERS_ROOT to the directory containing manifest.json and tables/.
No request parameter is interpreted as a filesystem path. Original source
values and columns are retained in TSV downloads.
"""
import csv
import hashlib
import io
import json
import math
import os
import re
from functools import lru_cache
from pathlib import Path

from django.conf import settings
from django.http import FileResponse, Http404, HttpResponse
from django.shortcuts import get_object_or_404
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.exceptions import APIException, ValidationError
from rest_framework.permissions import AllowAny
from rest_framework.response import Response

from .models import CellData

FIELDS = ('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'gene',
          'group1', 'group2', 'test_group_number', 'test_group')
NUMERIC = frozenset(('p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'test_group_number'))
COMPARISON = 'cell_type_vs_remaining_cells'


class MarkerAssetsUnavailable(APIException):
    status_code = 503
    default_detail = 'Marker files could not be read. Please retry later.'


def _root():
    return Path(getattr(settings, 'SCAID_MARKERS_ROOT',
                        os.environ.get('SCAID_MARKERS_ROOT',
                                       Path(getattr(settings, 'BASE_DIR', Path(__file__).resolve().parents[1])) / 'marker_assets'))).resolve()


@lru_cache(maxsize=4)
def _manifest(path, mtime_ns, size):
    value = json.loads(Path(path).read_text(encoding='utf-8'))
    if (value.get('schema_version') != 1 or value.get('comparison') != COMPARISON
            or tuple(value.get('columns', ())) != FIELDS or not isinstance(value.get('objects'), dict)):
        raise ValueError('Invalid marker manifest')
    return value['objects']


@lru_cache(maxsize=8)
def _rows(path, mtime_ns, size, checksum):
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != checksum:
        raise ValueError('Marker checksum mismatch')
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8'), newline=''), delimiter='\t')
    if tuple(reader.fieldnames or ()) != FIELDS:
        raise ValueError('Marker column mismatch')
    rows = []
    for row in reader:
        if None in row or any(row.get(key) is None or row[key] == '' for key in FIELDS):
            raise ValueError('Invalid marker row')
        for key in NUMERIC:
            number = float(row[key])
            if not math.isfinite(number) or (key in {'p_val', 'p_val_adj', 'pct.1', 'pct.2'} and not 0 <= number <= 1):
                raise ValueError('Invalid marker number')
        if row['group2'] != 'others':
            raise ValueError('Unexpected marker contrast')
        rows.append(row)
    return tuple(rows)


def _integer(params, name, default=None, maximum=100000):
    raw = params.get(name)
    if raw is None:
        if default is not None:
            return default
        raise ValidationError({name: 'This parameter is required.'})
    if not re.fullmatch(r'[0-9]{1,19}', str(raw)) or not 1 <= int(raw) <= maximum:
        raise ValidationError({name: f'Enter an integer from 1 to {maximum}.'})
    return int(raw)


def _text(params, name, limit):
    value = params.get(name, '')
    if len(value) > limit or any(ord(char) < 32 for char in value):
        raise ValidationError({name: f'Enter at most {limit} characters without control characters.'})
    return value.strip()


def _query(params):
    allowed = {'dataset_record_id', 'gene', 'cell_type', 'ordering', 'page', 'page_size', 'format'}
    unknown = set(params) - allowed
    if unknown:
        raise ValidationError({key: 'Unknown parameter.' for key in sorted(unknown)})
    record_id = _integer(params, 'dataset_record_id', maximum=9223372036854775807)
    ordering = params.get('ordering', 'p_val_adj')
    if ordering.lstrip('-') not in FIELDS or ordering.startswith('--'):
        raise ValidationError({'ordering': 'Choose a marker column, optionally prefixed with - for descending order.'})
    return {'record_id': record_id, 'gene': _text(params, 'gene', 128),
            'cell_type': _text(params, 'cell_type', 200), 'ordering': ordering,
            'page': _integer(params, 'page', default=1),
            'page_size': _integer(params, 'page_size', default=25, maximum=100)}


def _assets(record):
    try:
        root = _root()
        path = root / 'manifest.json'
        stat = path.stat()
        entry = _manifest(str(path), stat.st_mtime_ns, stat.st_size).get(str(record.pk))
        if entry is None:
            return None, None, ()
        # IDs are meaningful only with the original accession/condition/tissue.
        # A remapped database must never expose another object's results.
        if int(entry['dataset_record_id']) != record.pk or any(
                entry[key] != getattr(record, key) for key in ('abbreviation', 'dataset_id', 'tissue')):
            raise Http404('No marker table is registered for this object.')
        relative = Path(entry['relative_path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Invalid marker path')
        source = (root / relative).resolve()
        if not source.is_relative_to(root) or not source.is_file():
            raise ValueError('Invalid marker file')
        stat = source.stat()
        if stat.st_size != entry['bytes']:
            raise ValueError('Marker size mismatch')
        rows = _rows(str(source), stat.st_mtime_ns, stat.st_size, entry['sha256'])
        if len(rows) != entry['row_count']:
            raise ValueError('Marker row count mismatch')
        return entry, source, rows
    except Http404:
        raise
    except (OSError, ValueError, KeyError, TypeError, UnicodeError, csv.Error):
        raise MarkerAssetsUnavailable() from None


def _filtered(rows, query):
    gene, group = query['gene'].casefold(), query['cell_type']
    selected = [row for row in rows if (not gene or gene in row['gene'].casefold())
                and (not group or row['group1'] == group)]
    column = query['ordering'].lstrip('-')
    key = (lambda row: float(row[column])) if column in NUMERIC else (lambda row: row[column].casefold())
    return sorted(selected, key=key, reverse=query['ordering'].startswith('-'))


class MarkerViewSet(viewsets.ViewSet):
    """Precomputed cell-type-versus-rest markers, scoped to one data object."""
    http_method_names = ['get', 'head', 'options']
    permission_classes = [AllowAny]

    def list(self, request):
        query = _query(request.query_params)
        record = get_object_or_404(CellData, pk=query['record_id'])
        entry, _, rows = _assets(record)
        selected = _filtered(rows, query)
        start = (query['page'] - 1) * query['page_size']
        if query['page'] > 1 and start >= len(selected):
            raise Http404('This marker results page does not exist.')
        output = []
        for row in selected[start:start + query['page_size']]:
            output.append({key: (int(row[key]) if key == 'test_group_number' else float(row[key]))
                           if key in NUMERIC else row[key] for key in FIELDS})
        return Response({'dataset_record_id': record.pk, 'available': entry is not None,
                         'comparison': COMPARISON, 'columns': FIELDS,
                         'cell_types': entry['cell_types'] if entry else [],
                         'total_rows': len(rows), 'count': len(selected),
                         'page': query['page'], 'page_size': query['page_size'],
                         'ordering': query['ordering'], 'results': output,
                         'download_url': f'/api/markers/download/?dataset_record_id={record.pk}' if entry else None,
                         'note': ('Each cell type is compared with the remaining cells in this object. '
                                  'These are cell-level comparisons; source p values of zero are retained as reported.')
                         if entry else 'A searchable marker table has not been published for this object.'})

    @action(detail=False, methods=['get'], url_path='download')
    def download(self, request):
        query = _query(request.query_params)
        record = get_object_or_404(CellData, pk=query['record_id'])
        entry, source, rows = _assets(record)
        if entry is None:
            raise Http404('A marker table has not been published for this object.')
        filename = re.sub(r'[^A-Za-z0-9._-]', '_', entry['filename'])
        if not query['gene'] and not query['cell_type'] and 'ordering' not in request.query_params:
            response = FileResponse(source.open('rb'), as_attachment=True, filename=filename,
                                    content_type='text/tab-separated-values; charset=utf-8')
            response['X-Content-SHA256'] = entry['sha256']
        else:
            text = io.StringIO(newline='')
            writer = csv.DictWriter(text, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(_filtered(rows, query))
            response = HttpResponse(text.getvalue(), content_type='text/tab-separated-values; charset=utf-8')
            response['Content-Disposition'] = f'attachment; filename="filtered_{filename}"'
        response['X-Content-Type-Options'] = 'nosniff'
        return response

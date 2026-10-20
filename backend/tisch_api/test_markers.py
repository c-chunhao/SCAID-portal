"""Marker browser scope, numeric ordering, downloads and input safety."""
import csv
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from django.test import TestCase, override_settings
from django.urls import include, path
from rest_framework.routers import DefaultRouter
from rest_framework.test import APIClient

from .marker_views import COMPARISON, FIELDS, MarkerViewSet, _manifest, _rows
from .models import CellData

router = DefaultRouter()
router.register('markers', MarkerViewSet, basename='markers')
urlpatterns = [path('api/', include(router.urls))]


@override_settings(ROOT_URLCONF='tisch_api.test_markers')
class MarkerBrowserTests(TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.override = override_settings(SCAID_MARKERS_ROOT=self.root)
        self.override.enable()
        self.addCleanup(self.override.disable)
        self.client = APIClient()
        self.first = CellData.objects.create(abbreviation='A', disease_name='A', dataset_id='GSE123', tissue='PBMC')
        self.second = CellData.objects.create(abbreviation='B', disease_name='B', dataset_id='GSE123', tissue='PBMC')
        self.old = CellData.objects.create(abbreviation='Old', disease_name='Old', dataset_id='GSE001', tissue='PBMC')
        self.rows = [
            ['1e-251', '2.5', '.9', '.1', '1e-250', 'NKG7', 'NK cells', 'others', '2', 'NK cells;T cells'],
            ['0', '-1.75', '.05', '.8', '0', 'CD3D', 'NK cells', 'others', '1', 'NK cells'],
            ['1e-6', '10', '.7', '.2', '2e-5', 'NKG7', 'T cells', 'others', '2', 'NK cells;T cells'],
            ['.01', '1.25', '.6', '.4', '.05', 'C2orf15', 'T cells', 'others', '1', 'T cells'],
        ]
        text = io.StringIO(newline='')
        writer = csv.writer(text, delimiter='\t', lineterminator='\n')
        writer.writerow(FIELDS)
        writer.writerows(self.rows)
        self.raw = text.getvalue().encode()
        self.manifest = {'schema_version': 1, 'comparison': COMPARISON, 'columns': FIELDS, 'objects': {}}
        for record, rows in ((self.first, self.raw), (self.second, self.raw.replace(b'NKG7', b'SECOND'))):
            target = self.root / f'{record.pk}.tsv'
            target.write_bytes(rows)
            self.manifest['objects'][str(record.pk)] = {
                'dataset_record_id': record.pk, 'abbreviation': record.abbreviation,
                'dataset_id': record.dataset_id, 'tissue': record.tissue,
                'relative_path': target.name, 'row_count': 4, 'bytes': len(rows),
                'sha256': hashlib.sha256(rows).hexdigest(),
                'cell_types': ['NK cells', 'T cells'], 'filename': f'{record.abbreviation}_markers.tsv'}
        self.write_manifest()

    def write_manifest(self):
        _manifest.cache_clear()
        _rows.cache_clear()
        (self.root / 'manifest.json').write_text(json.dumps(self.manifest))

    def get(self, **params):
        return self.client.get('/api/markers/', {'dataset_record_id': self.first.pk, **params})

    def test_numeric_sort_preserves_source_zero_negative_and_extreme_values(self):
        response = self.get(ordering='-avg_log2FC')
        self.assertEqual(response.status_code, 200)
        self.assertEqual([row['avg_log2FC'] for row in response.json()['results']], [10, 2.5, 1.25, -1.75])
        response = self.get(ordering='p_val_adj')
        self.assertEqual([row['p_val_adj'] for row in response.json()['results']], [0, 1e-250, 2e-5, .05])
        self.assertEqual(response.json()['comparison'], COMPARISON)

    def test_gene_and_cell_type_filters_are_combined_within_one_exact_object(self):
        data = self.get(gene='nkg', cell_type='NK cells').json()
        self.assertEqual(data['count'], 1)
        self.assertEqual(data['results'][0]['avg_log2FC'], 2.5)
        second = self.client.get('/api/markers/', {'dataset_record_id': self.second.pk, 'gene': 'NKG7'}).json()
        self.assertEqual(second['count'], 0)
        self.assertTrue(second['available'])
        self.assertEqual(self.get(gene='C2orf15').json()['results'][0]['gene'], 'C2orf15')

    def test_pagination_is_stable_and_does_not_repeat_previous_object(self):
        first = self.get(ordering='-avg_log2FC', page_size=2, page=1).json()
        second = self.get(ordering='-avg_log2FC', page_size=2, page=2).json()
        self.assertEqual(first['count'], 4)
        self.assertEqual([row['avg_log2FC'] for row in first['results'] + second['results']], [10, 2.5, 1.25, -1.75])
        self.assertEqual(self.get(page_size=2, page=3).status_code, 404)

    def test_legacy_object_returns_explicit_empty_state(self):
        response = self.client.get('/api/markers/', {'dataset_record_id': self.old.pk})
        self.assertEqual(response.status_code, 200)
        self.assertFalse(response.json()['available'])
        self.assertEqual(response.json()['results'], [])
        self.assertIsNone(response.json()['download_url'])
        self.assertEqual(self.client.get('/api/markers/download/', {'dataset_record_id': self.old.pk}).status_code, 404)

    def test_invalid_inputs_and_write_methods_are_rejected(self):
        self.assertEqual(self.client.get('/api/markers/').status_code, 400)
        for params in ({'dataset_record_id': '../1'}, {'dataset_record_id': '-1'},
                       {'page_size': 101}, {'page_size': 0}, {'page': 'NaN'},
                       {'gene': 'x' * 129}, {'cell_type': 'a\n'},
                       {'ordering': '__class__'}, {'ordering': '--gene'},
                       {'relative_path': '../../etc/passwd'}):
            with self.subTest(params=params):
                self.assertEqual(self.get(**params).status_code, 400)
        self.assertEqual(self.get(dataset_record_id=999999).status_code, 404)
        self.assertEqual(self.client.post('/api/markers/', {}).status_code, 405)
        self.assertEqual(self.client.delete('/api/markers/download/').status_code, 405)

    def test_database_id_metadata_mismatch_never_returns_other_object_results(self):
        self.manifest['objects'][str(self.first.pk)]['abbreviation'] = 'B'
        self.write_manifest()
        self.assertEqual(self.get().status_code, 404)

    def test_table_path_traversal_and_changed_content_fail_closed(self):
        entry = self.manifest['objects'][str(self.first.pk)]
        entry['relative_path'] = '../outside.tsv'
        self.write_manifest()
        self.assertEqual(self.get().status_code, 503)
        entry['relative_path'] = f'{self.first.pk}.tsv'
        self.write_manifest()
        (self.root / entry['relative_path']).write_bytes(self.raw.replace(b'NKG7', b'EVIL'))
        self.assertEqual(self.get().status_code, 503)

    def test_original_download_is_byte_exact_and_filtered_download_keeps_precision(self):
        response = self.client.get('/api/markers/download/', {'dataset_record_id': self.first.pk})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(b''.join(response.streaming_content), self.raw)
        self.assertEqual(response['X-Content-SHA256'], hashlib.sha256(self.raw).hexdigest())
        response = self.client.get('/api/markers/download/', {'dataset_record_id': self.first.pk, 'gene': 'NKG7', 'ordering': '-avg_log2FC'})
        data = list(csv.DictReader(io.StringIO(response.content.decode()), delimiter='\t'))
        self.assertEqual(len(data), 2)
        self.assertEqual([row['avg_log2FC'] for row in data], ['10', '2.5'])
        self.assertEqual(data[1]['p_val_adj'], '1e-250')
        self.assertEqual(list(data[0]), list(FIELDS))

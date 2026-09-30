"""Regression for repeated full-catalog queries in disease filter options."""
from django.test import TestCase
from rest_framework.test import APIClient

from .models import CellData


class CatalogTreeQueryTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        CellData.objects.create(abbreviation='AS', disease_name='Ankylosing Spondylitis',
                                dataset_source='GEO', dataset_id='GSE1', tissue='PBMC',
                                batch='20260101', sample_size=3, cell_count=100, remarks='first')
        CellData.objects.create(abbreviation='AS', disease_name='Ankylosing Spondylitis',
                                dataset_source='GEO', dataset_id='GSE1', tissue='SFMC',
                                batch='20260102', sample_size=4, cell_count=200, remarks='second')
        CellData.objects.create(abbreviation='AS', disease_name='Ankylosing Spondylitis',
                                dataset_source='GEO', dataset_id='GSE2', tissue='PBMC',
                                batch='20260103', sample_size=5, cell_count=300)
        CellData.objects.create(abbreviation='RA', disease_name='Rheumatoid Arthritis',
                                dataset_source='GEO', dataset_id='GSE3', tissue='PBMC',
                                batch='20260104', sample_size=6, cell_count=400)
        CellData.objects.create(disease_full_name='Unmapped condition A',
                                dataset_source='GEO', dataset_id='GSE4', tissue='Skin',
                                sample_size=2, cell_count=50)
        CellData.objects.create(disease_full_name='Unmapped condition B',
                                dataset_source='GEO', dataset_id='GSE4', tissue='Skin',
                                sample_size=2, cell_count=60)

    def test_filter_tree_uses_one_query_and_preserves_disease_dataset_tissue_batches(self):
        with self.assertNumQueries(1):
            response = APIClient().get('/api/cell-data/', {'format': 'json'})
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertEqual(data['code'], 200)
        self.assertEqual(data['message'], 'success')
        self.assertEqual([r['abbreviation'] for r in data['data']],
                         ['AS', 'RA', 'Unmapped condition A', 'Unmapped condition B'])
        self.assertEqual(data['data'][0], {
            'abbreviation': 'AS', 'disease_label': 'Ankylosing Spondylitis', 'disease_name': 'Ankylosing Spondylitis',
            'disease_full_name': None, 'dataset_source': 'GEO',
            'children': [
                {'dataset_id': 'GSE1', 'dataset_label': 'GSE1', 'children': [
                    {'tissue': 'PBMC', 'tissue_label': 'PBMC', 'children': [
                        {'batch': '20260101', 'sample_size': 3, 'cell_count': 100, 'remarks': 'first'}]},
                    {'tissue': 'SFMC', 'tissue_label': 'SFMC', 'children': [
                        {'batch': '20260102', 'sample_size': 4, 'cell_count': 200, 'remarks': 'second'}]},
                ]},
                {'dataset_id': 'GSE2', 'dataset_label': 'GSE2', 'children': [
                    {'tissue': 'PBMC', 'tissue_label': 'PBMC', 'children': [
                        {'batch': '20260103', 'sample_size': 5, 'cell_count': 300, 'remarks': None}]},
                ]},
            ],
        })
        self.assertEqual(data['data'][1]['children'], [
            {'dataset_id': 'GSE3', 'dataset_label': 'GSE3', 'children': [
                {'tissue': 'PBMC', 'tissue_label': 'PBMC', 'children': [
                    {'batch': '20260104', 'sample_size': 6, 'cell_count': 400, 'remarks': None}]},
            ]},
        ])
        # A shared accession and null abbreviation must not merge distinct conditions.
        for row, expected_cells in zip(data['data'][2:], [50, 60]):
            self.assertEqual(row['children'], [
                {'dataset_id': 'GSE4', 'dataset_label': 'GSE4', 'children': [
                    {'tissue': 'Skin', 'tissue_label': 'Skin', 'children': [
                        {'batch': None, 'sample_size': 2, 'cell_count': expected_cells, 'remarks': None}]},
                ]},
            ])


class DisplayLabelTests(TestCase):
    """English display labels for records stored with Chinese tissue names, legacy
    disease names or truncated accessions; stored values stay unchanged for filters."""

    def test_reconciled_disease_and_tissue_labels(self):
        psa = CellData.objects.create(abbreviation='AS', disease_name='强直性脊柱炎', dataset_source='E-MTAB', dataset_id='8207', tissue='PBMC')
        bd = CellData.objects.create(abbreviation='SV', disease_name='血管炎', dataset_id='GSE198616', tissue='外周血')
        client = APIClient()
        rows = {row['id']: row for row in client.get('/api/cell-data-all/', {'format': 'json'}).json()}
        self.assertEqual(rows[psa.pk]['disease_label'], 'Psoriatic Arthritis')
        self.assertEqual(rows[psa.pk]['dataset_label'], 'E-MTAB-8207')
        self.assertEqual(rows[psa.pk]['dataset_id'], '8207')
        self.assertEqual(rows[bd.pk]['disease_label'], "Behçet's Disease")
        self.assertEqual((rows[bd.pk]['tissue'], rows[bd.pk]['tissue_label']), ('外周血', 'PBMC'))
        tree = {node['abbreviation']: node for node in client.get('/api/cell-data/', {'format': 'json'}).json()['data']}
        self.assertEqual(tree['BD']['disease_label'], "Behçet's Disease")
        CellData.objects.create(abbreviation='AS', disease_name='强直性脊柱炎', dataset_id='GSE216883', tissue='PBMC/SFMC')
        tree = {node['abbreviation']: node for node in client.get('/api/cell-data/', {'format': 'json'}).json()['data']}
        self.assertEqual(tree['AS']['disease_label'], 'Ankylosing Spondylitis')
        self.assertEqual(tree['PsA']['disease_label'], 'Psoriatic Arthritis')
        self.assertEqual(tree['BD']['children'][0]['children'][0]['tissue_label'], 'PBMC')
        self.assertEqual(tree['PsA']['children'][0]['dataset_label'], 'E-MTAB-8207')
        # Corrected public codes stay separate; the documented SV alias still works.
        self.assertEqual([r['id'] for r in client.get('/api/cell-data-all/', {'name': 'psoriatic'}).json()], [psa.pk])
        self.assertEqual([r['id'] for r in client.get('/api/cell-data-all/', {'name': 'Behçet'}).json()], [bd.pk])
        self.assertEqual([r['id'] for r in client.get('/api/cell-data-all/', {'abbreviation': 'SV'}).json()], [bd.pk])

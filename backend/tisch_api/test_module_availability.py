"""Publication gates must not offer unfinished cohorts in exploration filters."""
from unittest.mock import patch

from django.test import TestCase
from rest_framework.test import APIClient

from .models import CellData, PDFImage


class ModuleAvailabilityTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        for condition in ('Legacy', 'Ready', 'MainOnly'):
            CellData.objects.create(abbreviation=condition, disease_name=condition,
                                    dataset_source='GEO', dataset_id='GSE-test-' + condition,
                                    tissue='PBMC')

    @staticmethod
    def entry(record):
        if record.abbreviation == 'Legacy':
            return None
        return {'module_availability': {
            'gene_umap': record.abbreviation == 'Ready',
            'kegg_score': record.abbreviation == 'Ready',
            'cellchat': False, 'downloads': False,
        }, 'availability_notes': ['Prepared cohort release.']}

    def test_tree_excludes_unpublished_modules_and_retains_legacy(self):
        with patch('tisch_api.views.manifest_entry', side_effect=self.entry):
            for module in ('gene_umap', 'kegg_score'):
                with self.assertNumQueries(1):
                    response = APIClient().get('/api/cell-data/', {'module': module})
                self.assertEqual(response.status_code, 200)
                self.assertEqual([r['abbreviation'] for r in response.json()['data']], ['Legacy', 'Ready'])
            response = APIClient().get('/api/cell-data/', {'module': 'cellchat'})
            self.assertEqual([r['abbreviation'] for r in response.json()['data']], ['Legacy'])

    def test_unfiltered_tree_keeps_main_plot_only_cohorts(self):
        response = APIClient().get('/api/cell-data/')
        self.assertEqual(len(response.json()['data']), 3)

    def test_unknown_module_does_not_return_full_catalogue(self):
        self.assertEqual(APIClient().get('/api/cell-data/', {'module': 'typo'}).status_code, 400)

    def test_serialized_flags_and_notes(self):
        from .serializers import CellDataSerializer
        record = CellData.objects.get(abbreviation='MainOnly')
        with patch('tisch_api.serializers.manifest_entry', side_effect=self.entry):
            serializer = CellDataSerializer()
            self.assertFalse(serializer.get_module_availability(record)['gene_umap'])
            self.assertEqual(serializer.get_availability_notes(record), ['Prepared cohort release.'])


class CaseDistinctFeatureTests(TestCase):
    def test_source_case_variants_remain_separately_selectable(self):
        record = CellData.objects.create(abbreviation='DM', disease_name='Dermatomyositis',
                                         dataset_source='GEO', dataset_id='GSE303609', tissue='Skin')
        for gene in ('C2orf15', 'C2ORF15'):
            filename = 'Umap_' + gene + '.png'
            folder = 'Dermatomyositis/GSE303609/Skin/GeneUmap'
            PDFImage.objects.create(name=filename, path='/fig/' + folder + '/' + filename,
                                    image='/fig/' + folder + '/' + filename,
                                    relative_path=folder, page_number=1,
                                    image_width='1650', image_height='1500', is_pdf='False')
        with patch('tisch_api.views.case_sensitive_gene_symbols', return_value={'c2orf15'}):
            for gene in ('C2orf15', 'C2ORF15'):
                response = APIClient().get('/api/pdf-images/', {
                    'dataset_record_id': record.pk, 'GeneUmap': '1', 'pdf_name': gene,
                })
                self.assertEqual(response.status_code, 200)
                self.assertEqual(response.json()['count'], 1)
                self.assertEqual(response.json()['results']['data'][0]['name'], 'Umap_' + gene + '.png')
            # An ambiguous lower-case spelling must not silently combine the
            # two source features into one expression result.
            response = APIClient().get('/api/pdf-images/', {
                'dataset_record_id': record.pk, 'GeneUmap': '1', 'pdf_name': 'c2orf15',
            })
            self.assertEqual(response.json()['count'], 0)

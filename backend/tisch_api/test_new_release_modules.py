"""New release families must remain distinct and scoped to the selected object."""
from unittest.mock import patch
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from django.test import TestCase
from django.test import override_settings
from rest_framework.test import APIClient
from .models import CellData, PDFImage, H5File
from .serializers import H5FileSerializer


class NewModuleScopeTests(TestCase):
    def test_cached_matrix_gene_does_not_bypass_index_validation(self):
        from .thumbnails import resolved_figure_path
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            folder = 'Test/GSE1/PBMC/GeneUmap'
            path = root/folder/'Umap_CD3D.png'
            path.parent.mkdir(parents=True)
            path.write_bytes(b'existing cached bytes')
            record = SimpleNamespace(relative_path=folder,image=SimpleNamespace(name=str(path)))
            with override_settings(SCAID_FIGURE_ROOT=root), patch('tisch_api.asset_paths.figure_manifest_entry',return_value={'gene_render_mode':'rna_cache'}), patch('tisch_api.lazy_gene_plots.ensure_registered_gene_figure',return_value=None) as validate:
                self.assertIsNone(resolved_figure_path(record))
                validate.assert_called_once_with(record)

    def test_matrix_gene_http_cache_revalidates_changed_png(self):
        with TemporaryDirectory() as temporary:
            from PIL import Image
            path = Path(temporary)/'gene.png'
            Image.new('RGB',(20,20),'white').save(path)
            record = PDFImage.objects.create(name='Umap_CD3D.png',path=str(path),image=str(path),
                        relative_path='Test/GSE1/PBMC/GeneUmap',page_number=1,
                        image_width='20',image_height='20',is_pdf='False')
            with patch('tisch_api.views.resolved_figure_path',return_value=path), patch('tisch_api.views.figure_manifest_entry',return_value={'gene_render_mode':'rna_cache'}):
                first = APIClient().get(f'/api/pdf-images/{record.pk}/original/')
                self.assertEqual(first.status_code,200)
                self.assertEqual(first['Cache-Control'],'public, no-cache')
                etag = first['ETag']
                first.close()
                cached = APIClient().get(f'/api/pdf-images/{record.pk}/original/',HTTP_IF_NONE_MATCH=etag)
                self.assertEqual(cached.status_code,304)
                Image.new('RGB',(21,21),'black').save(path)
                changed = APIClient().get(f'/api/pdf-images/{record.pk}/original/',HTTP_IF_NONE_MATCH=etag)
                self.assertEqual(changed.status_code,200)
                self.assertNotEqual(changed['ETag'],etag)
                changed.close()

    def test_cellchat_overview_is_not_dataset_overview_or_another_object(self):
        rows = []
        for i in (1, 2):
            record = CellData.objects.create(abbreviation='Test'+str(i),
                       disease_name='Test'+str(i), dataset_source='GEO',
                       dataset_id='GSE-test-'+str(i), tissue='PBMC')
            folder = f'Test{i}/GSE-test-{i}/PBMC'
            for family in ('Dataset/Overview', 'CellChat/Overview'):
                PDFImage.objects.create(name='Count.png', relative_path=folder+'/'+family,
                    image='/fig/'+folder+'/'+family+'/Count.png', path='/fig/'+folder+'/'+family+'/Count.png',
                    page_number=1, image_width='1000', image_height='1000', is_pdf='False')
            rows.append((record, folder))
        record, folder = rows[0]
        with patch('tisch_api.views.manifest_folders', return_value=[folder]):
            response = APIClient().get('/api/pdf-images/', {
                'dataset_record_id': record.pk, 'pdf_name': 'CellChatOverview'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()['count'], 1)
        self.assertEqual(response.json()['results']['data'][0]['relative_path'], folder+'/CellChat/Overview')

    def test_marker_tree_respects_explicit_publication_flag(self):
        for name in ('ReadyMarkers', 'NoMarkers', 'Legacy'):
            CellData.objects.create(abbreviation=name, disease_name=name, dataset_source='GEO', dataset_id=name, tissue='PBMC')
        with patch('tisch_api.views.manifest_entry', side_effect=lambda row:
                   None if row.abbreviation == 'Legacy' else {'module_availability': {'marker_table': row.abbreviation == 'ReadyMarkers'}}):
            response = APIClient().get('/api/cell-data/', {'module': 'marker_table'})
        self.assertEqual(response.status_code, 200)
        self.assertEqual([row['abbreviation'] for row in response.json()['data']], ['ReadyMarkers'])

    def test_downloads_distinguish_h5ad_and_legacy_h5(self):
        for suffix, expected in (('.h5ad', 'H5AD'), ('.h5', 'H5')):
            record = H5File(name='Test_RNA'+suffix, path='/approved/Test_RNA'+suffix,
                            relative_path='new_aid_20261020/RNA', file_size=1024)
            serializer = H5FileSerializer()
            self.assertEqual(serializer.get_method(record), 'RNA')
            self.assertEqual(serializer.get_file_format(record), expected)

from pathlib import Path
from tempfile import TemporaryDirectory

from django.test import TestCase, override_settings
from rest_framework.test import APIClient

from .asset_paths import manifest_entry
from .models import CellData, GeneKeggEnum, H5File, PDFImage


class PublicApiTests(TestCase):
    def setUp(self):
        self.client = APIClient()
        self.record = CellData.objects.create(abbreviation='AS', disease_name='Ankylosing Spondylitis', dataset_id='GSE1', tissue='PBMC/SFMC', sample_size=4, cell_count=123)

    def image(self, path, name=None):
        # Ingestion stores the directory relative to the figure root
        # ('disease/dataset/tissue/family[/sub]'), which the API matches exactly.
        parts = list(Path(path).parent.parts)
        while parts and parts[0] in ('/', 'data', 'Fig', 'home', 'KWQ', '20250706', 'pdf_images'):
            parts.pop(0)
        return PDFImage.objects.create(name=name or Path(path).name, path=path, relative_path='/'.join(parts), page_number=1, image=path)

    def test_cellchat_tab_is_scoped_by_indexed_folder_for_manifest_records(self):
        record = CellData.objects.create(abbreviation='SS', dataset_id='GSE157278', tissue='PBMC')
        wanted = self.image('/data/Fig/SS/GSE157278/PBMC/CellChat/LRcircle/CCL_LRcircle.png')
        self.image('/data/Fig/SS/GSE157278/PBMC/CellChat/LRcontribution/CCL_LRcontribution.png')
        self.image('/data/Fig/SS/GSE157278/SFMC/CellChat/LRcircle/CCL_LRcircle.png')
        self.image('/data/Fig/RA/GSE159117/PBMC/CellChat/LRcircle/CCL_LRcircle.png')
        with self.assertNumQueries(3):  # record lookup, count, page
            rows = self.images(dataset_record_id=record.pk, pdf_name='LRcircle')
        self.assertEqual([r['id'] for r in rows], [wanted.pk])
        # Legacy alias directory (pSS) is covered through the manifest aliases.
        alias = self.image('pdf_images/pSS/GSE157278/PBMC/CellChat/Single/Single.png')
        self.assertEqual([r['id'] for r in self.images(dataset_record_id=record.pk, pdf_name='Single')], [alias.pk])
        kegg = self.image('/data/Fig/SS/GSE157278/PBMC/KEGGScore/UCell/hsa00010.png')
        self.image('/data/Fig/SS/GSE157278/PBMC/KEGGScore/AUCell/hsa00010.png')
        self.assertEqual([r['id'] for r in self.images(dataset_record_id=record.pk, pdf_name='hsa00010', KEGGScore=1, kegg_param='UCell')], [kegg.pk])

    def images(self, **params):
        response = self.client.get('/api/pdf-images/', params)
        self.assertEqual(response.status_code, 200)
        return response.json()['results']['data']

    def test_all_public_model_routes_reject_writes(self):
        routes = ('cell-data', 'cell-data-all', 'pdf-images', 'h5-file', 'gene-kegg-enum')
        for route in routes:
            for verb, suffix in [('post', ''), ('put', '1/'), ('patch', '1/'), ('delete', '1/')]:
                with self.subTest(route=route, verb=verb):
                    response = getattr(self.client, verb)(f'/api/{route}/{suffix}', {'name': 'overwritten'}, format='json')
                    self.assertEqual(response.status_code, 405)
        self.assertEqual(CellData.objects.get(pk=self.record.pk).dataset_id, 'GSE1')
        self.assertEqual(CellData.objects.count(), 1)

    def test_name_search_matches_catalog_and_rejects_unknown(self):
        for value in ('GSE1', 'Ankylosing', 'PBMC'):
            self.assertEqual(len(self.client.get('/api/cell-data-all/', {'name': value}).json()), 1)
        self.assertEqual(self.client.get('/api/cell-data-all/', {'name': 'nonexistent-identifier'}).json(), [])

    def test_invalid_numeric_values_are_client_errors(self):
        for route, key, value in [
            ('cell-data-all', 'sample_size', 'oops'), ('cell-data-all', 'cell_count', '-1'),
            ('cell-data-all', 'cell_count', '99999999999999999999999999'),
            ('h5-file', 'min_size', 'oops'), ('h5-file', 'max_size', '1.2'),
            ('h5-file', 'dataset_record_id', 'null'), ('pdf-images', 'dataset_record_id', '0'),
            ('pdf-images', 'page_size', 'bad'), ('pdf-images', 'page', '0'),
            ('gene-kegg-enum', 'page_size', '0'),
        ]:
            with self.subTest(route=route, key=key):
                self.assertEqual(self.client.get(f'/api/{route}/', {key: value}).status_code, 400)

    def test_image_match_is_exact_and_handles_slash_tissue(self):
        match = self.image('pdf_images/AS/GSE1/PBMC_SFMC/Dataset/Overview/Umap_Sample_page_1.jpg', 'Umap_Sample.pdf')
        for path in [
            'pdf_images/AS/GSE10/PBMC_SFMC/Dataset/Overview/Umap_Sample_page_1.jpg',
            'pdf_images/AS/GSE1/PBMC/Dataset/Overview/Umap_Sample_page_1.jpg',
            'pdf_images/OTHER_AS/GSE1/PBMC_SFMC/Dataset/Overview/Umap_Sample_page_1.jpg',
            'pdf_images/AS/GSE1/PBMC_SFMC/Dataset/Marker/Main_Marker_page_1.jpg',
        ]:
            self.image(path)
        rows = self.images(dataset_record_id=self.record.pk, pdf_name='Overview', batch='obsolete-batch')
        self.assertEqual([r['id'] for r in rows], [match.pk])
        self.assertEqual(len(self.images(abbreviation='AS', dataset_id='GSE1', tissue='PBMC/SFMC', pdf_name='Overview')), 1)

    def test_metacharacters_are_literal_and_cannot_broaden_scope(self):
        self.image('pdf_images/AS/GSE1/PBMC_SFMC/GeneUmap/Umap_CD3D.png')
        for key in ('abbreviation', 'dataset_id', 'tissue'):
            self.assertEqual(self.images(**{key: '.*'}), [])

    def test_gene_name_and_kegg_alias_filter_exact_result_type(self):
        gene = self.image('pdf_images/AS/GSE1/PBMC/GeneUmap/Umap_CD3D.png')
        self.image('pdf_images/AS/GSE1/PBMC/GeneUmap/Umap_OTHER_CD3D.png')
        kegg = self.image('pdf_images/AS/GSE1/PBMC/KEGGScore/AUCell/hsa00010.png')
        self.image('pdf_images/AS/GSE1/PBMC/Other/AUCell/hsa00010.png')
        self.assertEqual([r['id'] for r in self.images(name='CD3D.png', GeneUmap=1)], [gene.pk])
        self.assertEqual([r['id'] for r in self.images(name='hsa00010', KEGGUmap=1, kegg_param='AUCell')], [kegg.pk])

    def test_filename_filter_does_not_match_other_basename(self):
        match = self.image('pdf_images/AS/GSE1/PBMC/CellChat/Count_page_1.jpg', 'Count.pdf')
        self.image('pdf_images/AS/GSE1/PBMC/CellChat/OtherCount_page_1.jpg', 'OtherCount.pdf')
        self.assertEqual([r['id'] for r in self.images(pdf_name='Count')], [match.pk])

    def test_single_cellchat_menu_loads_the_category_not_a_single_named_file(self):
        match = self.image('pdf_images/AS/GSE1/PBMC_SFMC/CellChat/Single/CXCL12_page_1.jpg', 'CXCL12.pdf')
        self.image('pdf_images/AS/GSE1/PBMC_SFMC/CellChat/Other/CXCL12_page_1.jpg', 'CXCL12.pdf')
        self.assertEqual([r['id'] for r in self.images(dataset_record_id=self.record.pk, pdf_name='Single')], [match.pk])

    def test_enumeration_type_is_case_insensitive(self):
        GeneKeggEnum.objects.create(tp='gene', name='CD3D')
        GeneKeggEnum.objects.create(tp='kegg', name='hsa00010')
        for value in ('Gene', 'GENE', ' gene '):
            result = self.client.get('/api/gene-kegg-enum/', {'tp': value}).json()
            self.assertEqual([r['name'] for r in result['results']['data']], ['CD3D'])

    def test_pdf_pagination_preserves_all_results(self):
        for i in range(12):
            self.image(f'pdf_images/AS/GSE{i}/PBMC/GeneUmap/Umap_CD3D.png')
        first = self.client.get('/api/pdf-images/', {'name': 'CD3D', 'GeneUmap': 1}).json()
        self.assertEqual(first['count'], 12)
        self.assertEqual(len(first['results']['data']), 10)
        second = self.client.get(first['links']['next']).json()
        self.assertEqual(len(second['results']['data']), 2)

    def test_h5_listing_requires_explicit_dataset_or_filename(self):
        H5File.objects.create(name='AS_GSE1_PBMC_SFMC.h5', path='/missing/file.h5', relative_path='RNA', file_size=1)
        for params in ({}, {'name': ''}, {'name': 'null'}, {'name': 'undefined'}, {'path': '/'}):
            self.assertEqual(self.client.get('/api/h5-file/', params).json(), [])

    def test_h5_listing_only_contains_current_dataset_and_tissue(self):
        wanted = H5File.objects.create(name='AS_GSE1_PBMC_SFMC.h5', path='/data/wanted.h5', relative_path='RNA', file_size=1)
        for name in ('AS_GSE10_PBMC_SFMC.h5', 'AS_GSE1_PBMC.h5', 'OTHER_GSE1_PBMC_SFMC.h5'):
            H5File.objects.create(name=name, path='/data/' + name, relative_path='RNA', file_size=1)
        result = self.client.get('/api/h5-file/', {'dataset_record_id': self.record.pk}).json()
        self.assertEqual([r['id'] for r in result], [wanted.pk])
        self.assertEqual(result[0]['download_url'], f'/api/h5-file/{wanted.pk}/download/')
        self.assertNotIn('path', result[0])
        self.assertEqual(self.client.get('/api/h5-file/', {'dataset_record_id': 999999}).status_code, 404)

    def test_h5_invalid_dates_and_inverted_ranges_are_400(self):
        for params in ({'start_date': 'yesterday'}, {'end_date': '9999-12-31'}, {'min_size': 20, 'max_size': 1}, {'start_date': '2026-09-09', 'end_date': '2026-09-01'}):
            self.assertEqual(self.client.get('/api/h5-file/', params).status_code, 400)

    def test_manifest_uses_explicit_disease_for_shared_accession(self):
        ci = CellData(abbreviation=None, disease_full_name='IBD_Colitis_inflamed', dataset_id='HRA000072', tissue='colon')
        cd = CellData(abbreviation=None, disease_full_name='IBD_Crohns_Disease', dataset_id='HRA000072', tissue='colon')
        self.assertIn('IBD_CI', manifest_entry(ci)['h5_paths']['RNA'])
        self.assertIn('IBD_CD', manifest_entry(cd)['h5_paths']['RNA'])
        uncertain = CellData(abbreviation='AIH', dataset_id='GSE216064', tissue='外周血')
        self.assertEqual(manifest_entry(uncertain)['source_row'], 2)
        self.assertIn('ncbi.nlm.nih.gov', manifest_entry(uncertain)['mapping_source'])
        self.assertIsNone(manifest_entry(CellData(abbreviation='AIH', dataset_id='GSE-unknown', tissue='PBMC')))

    def test_null_abbreviations_keep_disease_trees_separate(self):
        CellData.objects.create(disease_full_name='IBD_Colitis_inflamed', dataset_id='CI1', tissue='colon')
        CellData.objects.create(disease_full_name='IBD_Crohns_Disease', dataset_id='CD1', tissue='colon')
        trees = self.client.get('/api/cell-data/').json()['data']
        groups = {r['disease_full_name']: {c['dataset_id'] for c in r['children']} for r in trees}
        self.assertEqual(groups['IBD_Colitis_inflamed'], {'CI1'})
        self.assertEqual(groups['IBD_Crohns_Disease'], {'CD1'})

    def test_pdf_name_search_uses_exact_ingested_gene_filename(self):
        wanted = self.image('pdf_images/AS/GSE1/PBMC/GeneUmap/Umap_CD3D.png')
        self.image('pdf_images/AS/GSE1/PBMC/GeneUmap/Umap_OTHER_CD3D.png')
        self.assertEqual([r['id'] for r in self.images(pdf_name='CD3D', GeneUmap=1)], [wanted.pk])
        kegg = self.image('pdf_images/AS/GSE1/PBMC/KEGGScore/AUCell/hsa00010.png')
        self.assertEqual([r['id'] for r in self.images(pdf_name='hsa00010', KEGGUmap=1, kegg_param='AUCell')], [kegg.pk])

    def test_home_canonical_abbreviations_select_legacy_catalog_rows(self):
        fixtures = [
            ('pSS', dict(abbreviation='SS', dataset_id='GSE157278', tissue='PBMC')),
            ('AD', dict(abbreviation='Atopic Dermatitis', disease_full_name='Atopic Dermatitis', dataset_id='GSE222840', tissue='Skin Biopsy')),
            ('CD', dict(disease_full_name='IBD_Crohns_Disease', dataset_id='HRA000072', tissue='colon')),
            ('CI', dict(disease_full_name='IBD_Colitis_inflamed', dataset_id='HRA000072', tissue='colon')),
            ('UC', dict(disease_full_name='IBD_Ulcerative_Colitis', dataset_id='HRA000072', tissue='colon')),
            ('T1DM', dict(disease_name='type_1_diabetes_mellitus', dataset_id='GSE221297', tissue='PBMC')),
            ('BD(Uveitis)', dict(abbreviation='AIR', disease_full_name='Behcets_disease', dataset_id='bAF7fe', tissue='眼前房液（房水）')),
            ('VKHD(Uveitis)', dict(abbreviation='AIR', disease_full_name='Vogt_Koyanag_Harada_disease', dataset_id='bAF7fe', tissue='眼前房液（房水）')),
        ]
        for abbreviation, data in fixtures:
            wanted = CellData.objects.create(**data)
            for key in ('abbreviation', 'name'):
                result = self.client.get('/api/cell-data-all/', {key: abbreviation}).json()
                self.assertEqual([row['id'] for row in result], [wanted.pk])

    def test_null_catalog_abbreviation_displays_manifest_short_name(self):
        record = CellData.objects.create(disease_name='type_1_diabetes_mellitus', dataset_id='GSE221297', tissue='PBMC')
        row = self.client.get('/api/cell-data-all/', {'dataset_id': 'GSE221297'}).json()[0]
        self.assertEqual(row['abbreviation'], 'T1DM')
        self.assertEqual(row['asset_status'], 'linked')
        record.refresh_from_db()
        self.assertIsNone(record.abbreviation)

    def test_gene_canonical_scope_uses_catalog_asset_relationship(self):
        CellData.objects.create(disease_full_name='IBD_Crohns_Disease', dataset_id='HRA000072', tissue='colon')
        wanted = self.image('pdf_images/IBD_Crohns_Disease/HRA000072/colon/GeneUmap/Umap_CD3D.png')
        self.image('pdf_images/IBD_Ulcerative_Colitis/HRA000072/colon/GeneUmap/Umap_CD3D.png')
        self.assertEqual([row['id'] for row in self.images(abbreviation='CD', dataset_id='HRA000072', GeneUmap=1, pdf_name='CD3D')], [wanted.pk])

    def test_name_search_includes_actual_h5_basename(self):
        record = CellData.objects.create(abbreviation='AS', dataset_id='GSE216883', tissue='PBMC/SFMC')
        result = self.client.get('/api/cell-data-all/', {'name': 'AS_GSE216883_PBMC_SFMC.h5'}).json()
        self.assertEqual([row['id'] for row in result], [record.pk])
        self.assertEqual(self.client.get('/api/cell-data-all/', {'name': '/home/data/KWQ/'}).json(), [])

    def test_ss_manifest_supports_original_and_rendered_disease_directories(self):
        record = CellData.objects.create(abbreviation='SS', dataset_id='GSE157278', tissue='PBMC')
        overview = self.image('pdf_images/pSS/GSE157278/PBMC/Dataset/Overview/Umap_Sample_page_1.jpg', 'Umap_Sample.pdf')
        gene = self.image('/data/Fig/SS/GSE157278/PBMC/GeneUmap/Umap_CD3D.png')
        self.image('pdf_images/pSS/GSE157278/SFMC/Dataset/Overview/Umap_Sample_page_1.jpg', 'Umap_Sample.pdf')
        self.assertEqual([row['id'] for row in self.images(dataset_record_id=record.pk, pdf_name='Overview')], [overview.pk])
        self.assertEqual([row['id'] for row in self.images(abbreviation='pSS', dataset_id='GSE157278', GeneUmap=1, pdf_name='CD3D')], [gene.pk])


class DownloadSafetyTests(TestCase):
    def setUp(self):
        self.client = APIClient()
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'h5'
        self.root.mkdir()
        self.settings_override = override_settings(H5_DOWNLOAD_ROOTS=[str(self.root)], H5_DOWNLOAD_ACCEL_PREFIX='')
        self.settings_override.enable()
        self.addCleanup(self.settings_override.disable)

    def file_record(self, path):
        return H5File.objects.create(name=path.name, path=str(path), relative_path='RNA', file_size=10)

    def test_get_returns_file_and_head_reports_length(self):
        path = self.root / 'sample.h5'
        path.write_bytes(b'\x89HDF\r\n\x1a\nreal-data')
        record = self.file_record(path)
        url = f'/api/h5-file/{record.pk}/download/'
        response = self.client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(b''.join(response.streaming_content), path.read_bytes())
        self.assertEqual(response['Content-Type'], 'application/x-hdf5')
        self.assertIn('attachment;', response['Content-Disposition'])
        response.close()
        head = self.client.head(url)
        self.assertEqual(head.status_code, 200)
        self.assertEqual(int(head['Content-Length']), path.stat().st_size)
        self.assertEqual(b''.join(head.streaming_content), b'')
        head.close()

    def test_rejects_outside_root_prefix_traversal_and_symlink(self):
        sibling = self.root.with_name('h5-private')
        sibling.mkdir()
        secret = sibling / 'secret.h5'
        secret.write_bytes(b'private')
        link = self.root / 'symlink.h5'
        link.symlink_to(secret)
        for path in (secret, self.root / '../h5-private/secret.h5', link):
            record = self.file_record(path)
            self.assertEqual(self.client.get(f'/api/h5-file/{record.pk}/download/').status_code, 404)

    def test_missing_non_h5_and_relative_paths_are_not_served(self):
        text = self.root / 'credentials.txt'
        text.write_text('private')
        for path in (self.root / 'missing.h5', text, Path('relative.h5')):
            record = self.file_record(path)
            self.assertEqual(self.client.get(f'/api/h5-file/{record.pk}/download/').status_code, 404)

    def test_missing_root_configuration_denies_download(self):
        path = self.root / 'sample.h5'
        path.write_bytes(b'bytes')
        record = self.file_record(path)
        with override_settings(H5_DOWNLOAD_ROOTS=[]):
            self.assertEqual(self.client.get(f'/api/h5-file/{record.pk}/download/').status_code, 404)

    def test_nginx_accel_uses_encoded_validated_relative_path(self):
        path = self.root / 'RNA' / 'sample space.h5'
        path.parent.mkdir()
        path.write_bytes(b'HDFbytes')
        record = self.file_record(path)
        with override_settings(H5_DOWNLOAD_ACCEL_PREFIX='/_protected_h5/'):
            response = self.client.get(f'/api/h5-file/{record.pk}/download/')
            self.assertEqual(response.status_code, 200)
            self.assertEqual(response['X-Accel-Redirect'], '/_protected_h5/RNA/sample%20space.h5')
            self.assertEqual(int(response['Content-Length']), 8)

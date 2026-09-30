"""Regressions for canonical conditions, public image URLs and atomic imports."""
import io
import json
import multiprocessing
from pathlib import Path
from tempfile import TemporaryDirectory

from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase, override_settings
from PIL import Image
from rest_framework.test import APIClient

from .models import CellData, PDFImage, H5File, GeneKeggEnum
from .thumbnails import make_thumbnail


def concurrent_preview_request(pk, barrier, queue):
    try:
        barrier.wait(timeout=20)
        response = APIClient().get(f'/api/pdf-images/{pk}/thumb/')
        body = b''.join(response.streaming_content)
        response.close()
        with Image.open(io.BytesIO(body)) as image:
            image.verify()
        queue.put({'status': response.status_code, 'bytes': len(body)})
    except Exception as exc:
        queue.put({'error': type(exc).__name__, 'detail': str(exc)})


class CanonicalConditionTests(TestCase):
    def setUp(self):
        self.client = APIClient()
        self.psa = [CellData.objects.create(abbreviation='AS', dataset_id='E-MTAB-8207', tissue=tissue) for tissue in ('PBMC', 'SFMC')]
        self.as_record = CellData.objects.create(abbreviation='AS', dataset_id='GSE216883', tissue='PBMC/SFMC')
        self.bd = CellData.objects.create(abbreviation='SV', dataset_id='GSE198616', tissue='PBMC')
        self.uveitis = CellData.objects.create(abbreviation='AIR', disease_full_name='Behcets_disease', dataset_id='bAF7fe', tissue='眼前房液（房水）')

    def test_psa_and_true_as_are_distinct_in_filter_search_and_tree(self):
        for key in ('abbreviation', 'name'):
            self.assertEqual({row['id'] for row in self.client.get('/api/cell-data-all/', {key: 'PsA'}).json()}, {row.pk for row in self.psa})
            self.assertEqual([row['id'] for row in self.client.get('/api/cell-data-all/', {key: 'AS'}).json()], [self.as_record.pk])
        tree = {row['abbreviation']: row for row in self.client.get('/api/cell-data/').json()['data']}
        self.assertEqual([row['dataset_id'] for row in tree['AS']['children']], ['GSE216883'])
        self.assertEqual([row['dataset_id'] for row in tree['PsA']['children']], ['E-MTAB-8207'])
        self.psa[0].refresh_from_db()
        self.assertEqual(self.psa[0].abbreviation, 'AS')

    def test_bd_excludes_uveitis_and_keeps_documented_legacy_alias(self):
        for value in ('BD', 'SV'):
            row = self.client.get('/api/cell-data-all/', {'abbreviation': value}).json()[0]
            self.assertEqual(row['id'], self.bd.pk)
            self.assertEqual(row['abbreviation'], 'BD')
            self.assertEqual(row['disease_name'], "Behçet's Disease")
        rows = self.client.get('/api/cell-data-all/', {'abbreviation': 'BD(Uveitis)'}).json()
        self.assertEqual([row['id'] for row in rows], [self.uveitis.pk])
        self.assertEqual(rows[0]['abbreviation'], 'BD(Uveitis)')

    def test_canonical_gene_scope_uses_exact_original_storage_folder(self):
        for code, dataset, tissue in [('AS', 'E-MTAB-8207', 'PBMC'), ('AS', 'GSE216883', 'PBMC_SFMC'), ('SV', 'GSE198616', 'PBMC')]:
            path = f'/missing/{code}/{dataset}/{tissue}/GeneUmap/Umap_CD3D.png'
            PDFImage.objects.create(name='Umap_CD3D.png', path=path, image=path, relative_path=f'{code}/{dataset}/{tissue}/GeneUmap', page_number=1)
        for code, expected in [('PsA', 'AS/E-MTAB-8207/PBMC/GeneUmap'), ('AS', 'AS/GSE216883/PBMC_SFMC/GeneUmap'), ('BD', 'SV/GSE198616/PBMC/GeneUmap')]:
            data = self.client.get('/api/pdf-images/', {'abbreviation': code, 'name': 'CD3D', 'GeneUmap': 1}).json()['results']['data']
            self.assertEqual([row['relative_path'] for row in data], [expected])
            self.assertEqual(data[0]['condition_abbreviation'], code)
            self.assertNotIn('path', data[0])
            self.assertEqual(data[0]['image'], f"/api/pdf-images/{data[0]['id']}/original/")
        mismatched = self.client.get('/api/pdf-images/', {'abbreviation': 'AS', 'dataset_id': 'E-MTAB-8207', 'name': 'CD3D', 'GeneUmap': 1}).json()
        self.assertEqual(mismatched['count'], 0)


class ConcurrentPreviewAndOriginalTests(TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'figures'
        self.thumbs = Path(self.temp.name) / 'thumbs'
        self.root.mkdir()
        self.source = self.root / 'AS/GSE1/PBMC/GeneUmap/Umap_CD3D.png'
        self.source.parent.mkdir(parents=True)
        Image.new('RGBA', (2000, 1800), (20, 90, 150, 255)).save(self.source)
        self.record = PDFImage.objects.create(name=self.source.name, image=str(self.source), path=str(self.source), relative_path='AS/GSE1/PBMC/GeneUmap', page_number=1)
        settings = override_settings(SCAID_FIGURE_ROOT=self.root, THUMBNAIL_ROOT=self.thumbs, THUMBNAIL_ACCEL_PREFIX='')
        settings.enable()
        self.addCleanup(settings.disable)

    def test_simultaneous_first_requests_from_four_processes_return_valid_jpegs(self):
        context = multiprocessing.get_context('fork')
        for round_number in range(3):
            target = self.thumbs / f'{self.record.pk % 1000:03d}' / f'{self.record.pk}.jpg'
            target.unlink(missing_ok=True)
            barrier = context.Barrier(4)
            queue = context.Queue()
            processes = [context.Process(target=concurrent_preview_request, args=(self.record.pk, barrier, queue)) for _ in range(4)]
            for process in processes:
                process.start()
            results = [queue.get(timeout=30) for _ in processes]
            for process in processes:
                process.join(timeout=30)
                self.assertFalse(process.is_alive())
                self.assertEqual(process.exitcode, 0)
            self.assertEqual([result.get('status') for result in results], [200] * 4, results)
            self.assertFalse(list(self.thumbs.rglob('*.tmp')))
            with Image.open(target) as image:
                image.verify()
            queue.close()

    def test_original_streams_safe_image_and_rejects_traversal_symlinks(self):
        client = APIClient()
        url = f'/api/pdf-images/{self.record.pk}/original/'
        response = client.get(url)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(b''.join(response.streaming_content), self.source.read_bytes())
        response.close()
        self.assertEqual(client.head(url).status_code, 200)
        outside = Path(self.temp.name) / 'outside.png'
        Image.new('RGB', (10, 10)).save(outside)
        link = self.root / 'outside.png'
        link.symlink_to(outside)
        self.record.image = str(link)
        self.record.save()
        self.assertEqual(client.get(url).status_code, 404)
        self.assertEqual(client.post(url).status_code, 405)

    def test_temporary_file_is_removed_after_failed_generation(self):
        from unittest.mock import patch
        target = self.thumbs / 'test.jpg'
        with patch('tisch_api.thumbnails.os.replace', side_effect=OSError('failed publication')):
            with self.assertRaises(OSError):
                make_thumbnail(self.source, target)
        self.assertFalse(list(self.thumbs.rglob('*.tmp')))

    def test_legacy_rendered_page_uses_only_explicit_extra_root(self):
        extra = Path(self.temp.name) / 'legacy/pdf_images'
        source = extra / 'AS/GSE1/PBMC/Dataset/Overview/Umap_Sample_page_1.jpg'
        source.parent.mkdir(parents=True)
        Image.new('RGB', (500, 300)).save(source)
        self.record.image = 'pdf_images/AS/GSE1/PBMC/Dataset/Overview/Umap_Sample_page_1.jpg'
        self.record.save()
        url = f'/api/pdf-images/{self.record.pk}/original/'
        self.assertEqual(APIClient().get(url).status_code, 404)
        with override_settings(SCAID_EXTRA_FIGURE_ROOTS=[str(extra)]):
            response = APIClient().get(url)
            self.assertEqual(response.status_code, 200)
            self.assertEqual(b''.join(response.streaming_content), source.read_bytes())
            response.close()


class PublicManifestImportTests(TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.figure = self.root / 'DEMO/DEMO1/PBMC/GeneUmap/Umap_CD3D.png'
        self.figure.parent.mkdir(parents=True)
        Image.new('RGB', (16, 12), (30, 40, 50)).save(self.figure)
        self.h5 = self.root / 'RNA/DEMO_DEMO1_PBMC.h5'
        self.h5.parent.mkdir()
        self.h5.write_bytes(b'\x89HDF\r\n\x1a\nfixture-header-only')
        self.manifest = {'schema_version': 1, 'release_id': 'synthetic-demo', 'source_url': 'https://example.org/synthetic-demo', 'datasets': [{'catalog': {'abbreviation': 'DEMO', 'disease_name': 'Synthetic demonstration', 'dataset_source': 'Demo', 'dataset_id': 'DEMO1', 'tissue': 'PBMC', 'sample_size': 0, 'cell_count': 0}, 'figures': [{'path': str(self.figure.relative_to(self.root))}], 'h5_files': [{'path': str(self.h5.relative_to(self.root)), 'method': 'RNA'}]}]}
        self.file = self.root / 'manifest.json'

    def command(self, apply=False):
        self.file.write_text(json.dumps(self.manifest))
        output = io.StringIO()
        call_command('import_public_manifest', str(self.file), figure_root=str(self.root), h5_root=str(self.root), apply=apply, stdout=output)
        return json.loads(output.getvalue())

    def test_default_dry_run_then_apply_and_repeat_are_idempotent(self):
        report = self.command()
        self.assertEqual(report['counts']['create'], 4)
        self.assertEqual(CellData.objects.count(), 0)
        self.assertFalse(report['patient_tables_imported'])
        self.assertNotIn(str(self.root), json.dumps(report))
        self.command(apply=True)
        self.assertEqual(CellData.objects.count(), 1)
        self.assertEqual(PDFImage.objects.count(), 1)
        self.assertEqual(H5File.objects.count(), 1)
        self.assertEqual(GeneKeggEnum.objects.count(), 1)
        report = self.command(apply=True)
        self.assertEqual(report['counts'], {'create': 0, 'update': 0, 'unchanged': 4})
        self.manifest['datasets'][0]['catalog']['cell_count'] = 1
        self.assertEqual(self.command(apply=True)['counts']['update'], 1)
        self.assertEqual(CellData.objects.get().cell_count, 1)

    def test_private_fields_and_repeated_catalogues_are_rejected_before_write(self):
        self.manifest['datasets'][0]['catalog']['patient_id'] = 'unsupported'
        with self.assertRaises(CommandError):
            self.command(apply=True)
        self.assertEqual(CellData.objects.count(), 0)
        del self.manifest['datasets'][0]['catalog']['patient_id']
        self.manifest['datasets'].append(self.manifest['datasets'][0])
        with self.assertRaises(CommandError):
            self.command(apply=True)
        self.assertEqual(CellData.objects.count(), 0)

    def test_invalid_file_and_symlink_outside_root_abort_import(self):
        self.h5.write_bytes(b'not-HDF5')
        with self.assertRaises(CommandError):
            self.command(apply=True)
        self.assertEqual(CellData.objects.count(), 0)
        self.h5.unlink()
        other = TemporaryDirectory()
        self.addCleanup(other.cleanup)
        foreign = Path(other.name) / 'foreign.h5'
        foreign.write_bytes(b'\x89HDF\r\n\x1a\n')
        self.h5.symlink_to(foreign)
        with self.assertRaises(CommandError):
            self.command(apply=True)
        self.assertEqual(CellData.objects.count(), 0)

    def test_duplicate_existing_asset_rolls_back_entire_transaction(self):
        for _ in range(2):
            H5File.objects.create(name=self.h5.name, path=str(self.h5), relative_path='RNA', file_size=25)
        with self.assertRaises(CommandError):
            self.command(apply=True)
        self.assertEqual(CellData.objects.count(), 0)
        self.assertEqual(PDFImage.objects.count(), 0)
        self.assertEqual(GeneKeggEnum.objects.count(), 0)

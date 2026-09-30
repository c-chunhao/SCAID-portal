"""Previews are generated once, cached, and never read files outside the figure roots."""
import shutil
from pathlib import Path

from django.conf import settings
from django.test import TestCase, override_settings
from PIL import Image
from rest_framework.test import APIClient

from .models import PDFImage

FIG = Path('/tmp/scaid-isolated-test-figures')
THUMBS = Path('/tmp/scaid-isolated-test-thumbs')


class ThumbnailTests(TestCase):
    def setUp(self):
        for path in (FIG, THUMBS):
            shutil.rmtree(path, ignore_errors=True)
            path.mkdir(parents=True)
        self.addCleanup(shutil.rmtree, FIG, True)
        self.addCleanup(shutil.rmtree, THUMBS, True)
        big = FIG / 'AS/GSE1/PBMC/GeneUmap/Umap_CD3D.png'
        big.parent.mkdir(parents=True)
        Image.new('RGBA', (3300, 3000), (200, 30, 30, 255)).save(big)
        self.png = PDFImage.objects.create(name='Umap_CD3D.png', path=str(big), relative_path='AS/GSE1/PBMC/GeneUmap', page_number=1, image=str(big))
        small = FIG / 'AS/GSE1/PBMC/Dataset/Overview/Umap_Sample_page_1.jpg'
        small.parent.mkdir(parents=True)
        Image.new('RGB', (561, 360)).save(small)
        self.jpg = PDFImage.objects.create(name='Umap_Sample.pdf', path=str(small), relative_path='AS/GSE1/PBMC/Dataset/Overview', page_number=1, image=str(small))
        outside = Path('/tmp/scaid-isolated-outside.png')
        Image.new('RGB', (10, 10)).save(outside)
        self.addCleanup(outside.unlink)
        self.outside = PDFImage.objects.create(name='x.png', path=str(outside), relative_path='x', page_number=1, image=str(outside))
        self.client = APIClient()

    def test_png_preview_is_generated_once_and_bounded(self):
        response = self.client.get(f'/api/pdf-images/{self.png.pk}/thumb/')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response['Content-Type'], 'image/jpeg')
        self.assertIn('max-age', response['Cache-Control'])
        body = b''.join(response.streaming_content)
        cached = THUMBS / f'{self.png.pk % 1000:03d}' / f'{self.png.pk}.jpg'
        self.assertTrue(cached.exists())
        with Image.open(cached) as im:
            self.assertLessEqual(max(im.size), settings.THUMBNAIL_MAX_SIDE)
        self.assertLess(len(body), 200_000)
        stamp = cached.stat().st_mtime_ns
        self.assertEqual(self.client.get(f'/api/pdf-images/{self.png.pk}/thumb/').status_code, 200)
        self.assertEqual(cached.stat().st_mtime_ns, stamp)

    def test_small_jpeg_and_foreign_paths_have_no_preview(self):
        self.assertEqual(self.client.get(f'/api/pdf-images/{self.jpg.pk}/thumb/').status_code, 404)
        self.assertEqual(self.client.get(f'/api/pdf-images/{self.outside.pk}/thumb/').status_code, 404)
        self.assertEqual(self.client.get('/api/pdf-images/999999/thumb/').status_code, 404)
        self.assertEqual(self.client.post(f'/api/pdf-images/{self.png.pk}/thumb/').status_code, 405)

    def test_list_exposes_thumb_url_only_for_png(self):
        rows = {row['id']: row for row in self.client.get('/api/pdf-images/', {'name': 'CD3D', 'GeneUmap': 1}).json()['results']['data']}
        self.assertEqual(rows[self.png.pk]['thumb_url'], f'/api/pdf-images/{self.png.pk}/thumb/')
        overview = self.client.get('/api/pdf-images/', {'dataset_record_id': ''}).status_code
        self.assertIn(overview, (200, 400))

    @override_settings(THUMBNAIL_ACCEL_PREFIX='/_thumbs/')
    def test_accel_redirect_when_configured(self):
        response = self.client.get(f'/api/pdf-images/{self.png.pk}/thumb/')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response['X-Accel-Redirect'], f'/_thumbs/{self.png.pk % 1000:03d}/{self.png.pk}.jpg')

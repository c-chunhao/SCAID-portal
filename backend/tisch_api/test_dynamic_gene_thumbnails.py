"""Dynamic gene previews must revalidate across source and browser caches."""
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from django.test import TestCase, override_settings
from PIL import Image
from rest_framework.test import APIClient

from .models import PDFImage


class DynamicGeneThumbnailTests(TestCase):
    def test_preview_revalidates_changed_source_and_conditional_requests(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'Umap_CD3D.png'
            Image.new('RGB', (80, 60), 'white').save(source)
            record = PDFImage.objects.create(
                name=source.name, path=str(source), image=str(source),
                relative_path='Test/GSE1/PBMC/GeneUmap', page_number=1,
                image_width='80', image_height='60', is_pdf='False',
            )
            client = APIClient()
            endpoint = f'/api/pdf-images/{record.pk}/thumb/'
            with override_settings(THUMBNAIL_ROOT=root / 'thumbs',
                                   THUMBNAIL_ACCEL_PREFIX='/internal/thumbs/'), \
                    patch('tisch_api.thumbnails.resolved_figure_path', return_value=source), \
                    patch('tisch_api.views.figure_manifest_entry',
                          return_value={'gene_render_mode': 'rna_cache'}):
                first = client.get(endpoint)
                self.assertEqual(first.status_code, 200)
                self.assertEqual(first['Cache-Control'], 'public, no-cache')
                self.assertNotIn('X-Accel-Redirect', first)
                etag = first['ETag']
                original_preview = b''.join(first.streaming_content)
                first.close()

                cached = client.get(endpoint, HTTP_IF_NONE_MATCH=etag)
                self.assertEqual(cached.status_code, 304)
                self.assertEqual(cached['ETag'], etag)
                self.assertEqual(cached['Cache-Control'], 'public, no-cache')
                self.assertEqual(cached.content, b'')

                Image.new('RGB', (81, 61), 'black').save(source)
                changed = client.get(endpoint, HTTP_IF_NONE_MATCH=etag)
                self.assertEqual(changed.status_code, 200)
                self.assertNotEqual(changed['ETag'], etag)
                self.assertNotEqual(b''.join(changed.streaming_content), original_preview)
                changed_etag = changed['ETag']
                changed.close()

                current = client.get(endpoint, HTTP_IF_NONE_MATCH=changed_etag)
                self.assertEqual(current.status_code, 304)
                self.assertEqual(current['ETag'], changed_etag)

from django.test import TestCase
from rest_framework.test import APIClient

from .models import PDFImage


class WithheldCellChatTests(TestCase):
    def test_withheld_tissue_is_not_listed_or_served_but_other_figures_remain(self):
        client = APIClient()
        def image(path):
            return PDFImage.objects.create(name='Count.pdf', path=path,
                relative_path=path.rsplit('/', 1)[0], page_number=1, image=path)
        withheld = image('AS/E-MTAB-8207/SFMC/CellChat/Count_page_1.jpg')
        nested = image('AS/E-MTAB-8207/SFMC/CellChat/LRcircle/Count_page_1.jpg')
        other = image('AS/E-MTAB-8207/PBMC/CellChat/Count_page_1.jpg')
        gene = image('AS/E-MTAB-8207/SFMC/GeneUmap/Count_page_1.jpg')
        listed = client.get('/api/pdf-images/', {'pdf_name': 'Count', 'page_size': 100}).json()
        self.assertEqual({row['id'] for row in listed['results']['data']}, {other.pk, gene.pk})
        for pk in (withheld.pk, nested.pk):
            for suffix in ('', 'original/', 'thumb/'):
                self.assertEqual(client.get(f'/api/pdf-images/{pk}/{suffix}').status_code, 404)
        self.assertEqual(client.get(f'/api/pdf-images/{other.pk}/').status_code, 200)

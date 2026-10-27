from django.test import TestCase
from rest_framework.test import APIClient

from .models import CellData, PDFImage


class WithheldCellChatTests(TestCase):
    def test_only_reviewed_sfmc_release_is_listed_and_old_images_remain_hidden(self):
        client = APIClient()
        record = CellData.objects.create(abbreviation='AS', dataset_id='8207', dataset_source='E-MTAB', tissue='SFMC')
        def image(path, folder, name='Count.pdf'):
            return PDFImage.objects.create(name=name, path=path, relative_path=folder, page_number=1, image=path)
        root = 'AS/E-MTAB-8207/SFMC/CellChat'
        old = image(root+'/Count_page_1.jpg',root)
        nested = image(root+'/LRcircle/Count_page_1.jpg',root+'/LRcircle')
        correct = image('pdf_images/r20261027sfmc/'+'a'*64+'.png',root+'/Overview')
        wrong_name = image('pdf_images/r20261027sfmc/'+'b'*64+'.png',root+'/Overview','OtherCount.pdf')
        other = image('AS/E-MTAB-8207/PBMC/CellChat/Count_page_1.jpg','AS/E-MTAB-8207/PBMC/CellChat')
        gene = image('AS/E-MTAB-8207/SFMC/GeneUmap/Umap_CD3D.png','AS/E-MTAB-8207/SFMC/GeneUmap','Umap_CD3D.png')
        listed = client.get('/api/pdf-images/', {'pdf_name':'Count', 'page_size':100}).json()
        self.assertEqual({row['id'] for row in listed['results']['data']}, {correct.pk,other.pk})
        scoped = client.get('/api/pdf-images/', {'dataset_record_id':record.pk,'pdf_name':'Count'}).json()
        self.assertEqual([row['id'] for row in scoped['results']['data']], [correct.pk])
        for pk in (old.pk,nested.pk):
            for suffix in ('','original/','thumb/'):
                self.assertEqual(client.get(f'/api/pdf-images/{pk}/{suffix}').status_code,404)
        self.assertEqual(client.get(f'/api/pdf-images/{correct.pk}/').status_code,200)
        self.assertEqual(client.get(f'/api/pdf-images/{other.pk}/').status_code,200)

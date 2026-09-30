from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from tisch_api.catalog import collect_catalog


class CatalogTests(TestCase):
    def test_catalog_combines_datasets_and_score_methods_without_duplicate_terms(self):
        with TemporaryDirectory() as temp:
            root = Path(temp)
            for path in [
                'AS/GSE1/PBMC/GeneUmap/Umap_CD3D.png',
                'AS/GSE2/PBMC/GeneUmap/Umap_CD3D.png',
                'AS/GSE2/PBMC/GeneUmap/Umap_HLA-DRA.png',
                'AS/GSE1/PBMC/KEGGScore/AUCell/hsa00010.png',
                'AS/GSE1/PBMC/KEGGScore/UCell/hsa00010.png',
                'AS/GSE1/PBMC/Dataset/Overview/Umap_Sample.png',
                'AS/GSE1/PBMC/GeneUmap/readme.txt',
            ]:
                file = root / path
                file.parent.mkdir(parents=True, exist_ok=True)
                file.touch()
            self.assertEqual(collect_catalog(root), {('gene', 'CD3D'), ('gene', 'HLA-DRA'), ('kegg', 'hsa00010')})

    def test_missing_figure_root_does_not_silently_create_empty_catalog(self):
        with TemporaryDirectory() as temp:
            with self.assertRaises(ValueError):
                collect_catalog(Path(temp) / 'missing')

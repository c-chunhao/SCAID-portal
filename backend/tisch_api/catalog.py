"""Build the search vocabulary from published figure names, without opening images."""
from pathlib import Path


def collect_catalog(figure_root):
    root = Path(figure_root)
    if not root.is_dir():
        raise ValueError(f"Figure directory does not exist: {root}")
    terms = set()
    for disease in root.iterdir():
        if not disease.is_dir():
            continue
        for dataset in disease.iterdir():
            if not dataset.is_dir():
                continue
            for tissue in dataset.iterdir():
                if not tissue.is_dir():
                    continue
                genes = tissue / 'GeneUmap'
                if genes.is_dir():
                    for image in genes.iterdir():
                        if image.suffix.lower() in {'.png', '.jpg', '.jpeg'} and image.stem.startswith('Umap_'):
                            name = image.stem[5:]
                            if name:
                                terms.add(('gene', name))
                scores = tissue / 'KEGGScore'
                if scores.is_dir():
                    for method in scores.iterdir():
                        if method.is_dir():
                            for image in method.iterdir():
                                if image.suffix.lower() in {'.png', '.jpg', '.jpeg'}:
                                    terms.add(('kegg', image.stem))
    return terms

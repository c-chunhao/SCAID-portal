#!/usr/bin/env python3
"""Declare the audited GSE279086 processed-only exception before publication."""
import json
from pathlib import Path
import h5py
import numpy as np

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
P=WORK/'matrices/objects/T1D_GSE279086_Kidney/export_report.json'
RAW_STATUS='unavailable_exact_source_axis_mapping_incomplete'

def main():
    r=json.loads(P.read_text())
    if r['object_id']!='T1D/GSE279086/Kidney' or r['status']!='exported_awaiting_raw_count_recovery':
        raise ValueError('Processed-only declaration requires the explicit, original fractional-matrix export')
    path=Path(r['matrices']['RNA']['temporary_path'])
    with h5py.File(path,'r+') as h:
        if 'source_processed_counts' not in h['layers'] or 'counts' in h['layers']:
            raise ValueError('Incomplete or mislabeled count layer')
        text={'counts_representation':'No raw-count layer is provided. Original fractional processed RNA values are preserved as layers/source_processed_counts.',
              'nCount_RNA_semantics':'Original nCount_RNA and nFeature_RNA describe the preserved fractional source_processed_counts matrix. These are not declared raw UMI counts.',
              'raw_count_status':RAW_STATUS,
              'expression_status':'Original RNA normalized data layer retained as X; original fractional processed expression retained separately.',
              'raw_count_source_scope':'Raw and processed barcode axes could not be matched for all retained source cells; no missing raw counts were filled with zero or inferred from other samples.'}
        g=h['uns/scaid']
        for name,value in text.items():
            if name in g:del g[name]
            d=g.create_dataset(name,data=np.asarray([value],dtype=h5py.string_dtype('utf-8')))
            d.attrs['encoding-type']='string-array';d.attrs['encoding-version']='0.2.0'
        d=g.create_dataset('counts_available',data=np.asarray([False],dtype=bool))
        d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
    item=r['matrices']['RNA']
    item.update(counts_available=False,raw_count_status=RAW_STATUS,
                expression_status='source_normalized_data_and_fractional_processed_expression',
                content_description='Processed RNA expression; a complete integer raw-count layer is unavailable.',
                source_processed_counts_semantics='Original noninteger Seurat RNA counts layer retained exactly; not represented as raw UMI counts.')
    r.update(counts_available=False,raw_count_status=RAW_STATUS,
             expression_status=item['expression_status'],
             source_axis_audit_directory=str(WORK/'matrices/source_recovery_GSE279086'),
             status='exported_awaiting_python_validation')
    temp=P.with_suffix('.json.tmp');temp.write_text(json.dumps(r,ensure_ascii=False,indent=2));temp.replace(P)
    print('GSE279086_PROCESSED_ONLY_DECLARED',path)

if __name__=='__main__':main()

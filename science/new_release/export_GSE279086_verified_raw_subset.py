#!/usr/bin/env python3
"""Export only the 24 fully mapped GSMs, without inventing counts for the other two."""
import json,os,time,gzip
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import anndata as ad
from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from recover_GSE279086_raw_counts import sha,matrix,lines,unique_r,one

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
AUDIT=WORK/'matrices/source_recovery_GSE279086/axis_audit.json'
OUTPUT=EXPORT_ROOT/'RNA_counts_subset/T1D_GSE279086_Kidney_integer_raw_counts_24GSM_23979cells.h5ad'

def main():
    audit=json.loads(AUDIT.read_text())
    if audit['fully_mappable_GSM_count']!=24:raise ValueError('The independent 24-GSM audit is required')
    source_rna=Path(audit['rna_axis_path'])
    if not source_rna.exists():source_rna=Path(str(source_rna).removesuffix('.tmp'))
    with h5py.File(source_rna,'r') as h:obs=read_elem(h['obs']);target_genes=read_elem(h['var']).index
    blocks=[];positions=[];records=[];template=None
    for sample in audit['samples']:
        if not sample['safe_for_full_count_recovery']:continue
        gsm=sample['GSM'];print('RAW_SUBSET',gsm,flush=True)
        mapping=Path(sample['npz_mapping_path'])
        if sha(mapping)!=sample['npz_mapping_sha256']:raise ValueError('Independent axis mapping changed')
        with np.load(mapping) as z:
            p=z['target_cell_indices'];ci=z['raw_barcode_indices'];gi=z['raw_feature_indices']
        directory=Path(audit['source_root'])/gsm
        f=one(directory,'*_features.tsv.gz');bc=one(directory,'*_barcodes.tsv.gz');mtx=one(directory,'*_matrix.mtx.gz')
        features=[x.split('\t') for x in lines(f)];barcodes=lines(bc)
        if template is None:template=features
        elif features!=template:raise ValueError('Raw gene templates differ')
        unique_genes=unique_r([x[1] for x in features])
        if list(np.asarray(unique_genes)[gi])!=list(target_genes):raise ValueError('Independent gene axis mapping differs')
        if any(obs.index[k].rsplit('_',1)[-1]!=barcodes[c] for k,c in zip(p,ci)):raise ValueError('Raw barcode axis differs')
        if not obs.iloc[p]['GSM'].astype(str).eq(gsm).all():raise ValueError('GSM metadata differs from source mapping')
        raw,header=matrix(mtx)
        if 'coordinate integer general' not in header or raw.shape!=(len(features),len(barcodes)):
            raise ValueError('Invalid raw source declaration or shape')
        if not np.isfinite(raw.data).all() or (raw.data<0).any() or not np.equal(raw.data,np.rint(raw.data)).all():
            raise ValueError('Raw count values are not nonnegative integers')
        block=raw[:,ci].T.tocsr();block.sort_indices();block.eliminate_zeros()
        if block.data.max(initial=0)>np.iinfo(np.int32).max:raise ValueError('int32 overflow')
        block=block.astype(np.int32)
        del raw
        blocks.append(block);positions.extend(p.tolist())
        records.append({'GSM':gsm,'cells':len(p),'nnz':block.nnz,
                        'axis_mapping':str(mapping),'axis_mapping_sha256':sha(mapping),
                        'source_files':[{ 'path':str(x),'bytes':x.stat().st_size,'sha256':sha(x)} for x in (mtx,f,bc)]})
    if len(records)!=24 or len(positions)!=23979 or len(set(positions))!=23979:
        raise ValueError('Verified raw subset scope differs from 24 GSM / 23,979 retained cells')
    x=sparse.vstack(blocks,format='csr');del blocks
    ix=np.argsort(np.asarray(positions));x=x[ix,:];selected=np.asarray(positions)[ix]
    meta=obs.iloc[selected].copy()
    meta['source_raw_barcode']=[s.rsplit('_',1)[-1] for s in meta.index]
    meta['nCount_source_raw']=np.asarray(x.sum(1)).ravel().astype(float)
    meta['nFeature_source_raw']=np.diff(x.indptr).astype(float)
    mt=np.asarray([f[1].upper().startswith('MT-') for f in template])
    meta['percent_mt_source_raw']=np.asarray(x[:,mt].sum(1)).ravel()*100/np.maximum(meta['nCount_source_raw'],1)
    var=pd.DataFrame({'gene_id_source_raw':[f[0] for f in template],
                      'gene_symbol_source_raw':[f[1] for f in template],
                      'feature_type_source_raw':[f[2] if len(f)>2 else 'Gene Expression' for f in template]},
                     index=unique_r([f[1] for f in template]))
    a=ad.AnnData(X=x,obs=meta,var=var)
    scope={'source':'GSE279086 integer raw matrices supplied with the study',
           'content':'Raw counts for the 24 GSMs whose retained cell identifiers have complete, unique source-axis mappings.',
           'cells':23979,'genes':36601,'GSMs':24,'original_complete_processed_object_cells':24663,
           'omitted_GSMs':['GSM8561122','GSM8561123'],'omitted_processed_cells':684,
           'original_metadata':'Original source metadata retained; nCount_RNA and nFeature_RNA still describe the processed matrix. *_source_raw columns describe this raw-count matrix.',
           'raw_genes':'All 36,601 source gene rows retained; duplicate symbols are distinguished by the source-reader make.unique rule and explicit gene IDs.',
           'sample_labels':'Original sample labels are preserved and are not authenticated donor IDs.',
           'relationship_to_main_download':'A supplemental subset, not a complete replacement for the 24,663-cell processed RNA file.',
           'excluded_bad_sample_coincidences':'No cells from the two blocked GSMs were included, including three coincidentally matching barcodes.'}
    a.uns['source_scope_json']=json.dumps(scope,ensure_ascii=False)
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    temp=Path(str(OUTPUT)+'.tmp')
    a.write_h5ad(temp,compression='gzip',compression_opts=1)
    # A hard link avoids storing the same integer matrix twice on disk.
    with h5py.File(temp,'r+') as h:h['layers']['counts']=h['X']
    b=ad.read_h5ad(temp)
    if not a.obs_names.equals(b.obs_names) or not a.var_names.equals(b.var_names):raise ValueError('Subset axis readback differs')
    for matrix_out in [b.X,b.layers['counts']]:
        for name in ['data','indices','indptr']:
            if not np.array_equal(getattr(x,name),getattr(matrix_out,name)):raise ValueError('Subset full sparse readback differs')
    if not np.array_equal(np.asarray(b.X.sum(1)).ravel(),b.obs['nCount_source_raw']):raise ValueError('Raw subset QC counts differ')
    if OUTPUT.exists():raise FileExistsError('Refusing to overwrite existing verified raw subset')
    os.replace(temp,OUTPUT)
    manifest={'status':'validated_supplementary_subset','path':str(OUTPUT),'bytes':OUTPUT.stat().st_size,'sha256':sha(OUTPUT),
        'object_id':'T1D/GSE279086/Kidney','method':'RNA_counts_subset','counts_available':True,
        'raw_count_status':'verified_raw_subset_24_GSM_23979_cells','complete_object_counts':False,
        'shape':list(a.shape),'format':'AnnData H5AD','main_124_file_manifest_member':False,
        'validation':{'independent_axis_audit':True,'all_raw_values_nonnegative_integer':True,
                      'AnnData_reopen':True,'full_sparse_components_equal':True,'exact_cell_and_gene_axes':True,
                      'raw_QC_matches':True},'scope':scope,'source_samples':records,
        'completed':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    dest=WORK/'matrices/supplementary_raw_subset_manifest.json'
    dest.write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps({k:manifest[k] for k in ['status','path','bytes','sha256','shape']}),flush=True)

if __name__=='__main__':main()

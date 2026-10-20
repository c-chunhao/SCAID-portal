#!/usr/bin/env python3
"""Independent AnnData validation and atomic promotion of a single object."""
import argparse, gc, hashlib, json, os, time
from importlib.metadata import version
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while b := f.read(8 * 1024 * 1024): h.update(b)
    return h.hexdigest()

def vector(x):
    return x if isinstance(x, list) else [x]

def check_probes(matrix, probes):
    rows = np.asarray(vector(probes['cell_index_zero']),dtype=int)
    cols = np.asarray(vector(probes['feature_index_zero']),dtype=int)
    expected = np.asarray(vector(probes['values']),dtype=float)
    got = np.zeros(len(rows), dtype=float)
    for row in np.unique(rows):
        mask = rows == row
        z = matrix[int(row):int(row)+1, :][:, cols[mask]]
        got[mask] = z.toarray().ravel() if sparse.issparse(z) else np.asarray(z).ravel()
    if not np.allclose(got,expected,rtol=1e-13,atol=1e-14,equal_nan=False):
        raise ValueError(f'Independent AnnData probes differ; maximum abs error {np.max(np.abs(got-expected))}')
    return {'source_value_probes':len(rows),'passed':True,'maximum_absolute_error':float(np.max(np.abs(got-expected)))}

def validate(report_path):
    report_path = Path(report_path)
    report = json.loads(report_path.read_text())
    if report['status'] == 'validated':
        # Resumption never trusts a status flag alone.
        for item in report['matrices'].values():
            if not Path(item['path']).is_file() or sha256(item['path']) != item['sha256']:
                raise ValueError('Previously validated file missing or changed: '+item['path'])
        return report
    if report['status'] != 'exported_awaiting_python_validation':
        raise ValueError('Object export is incomplete: '+report['status'])
    started = time.time()
    meta_path = report_path.parent / 'cell_metadata.tsv.gz'
    original = pd.read_csv(meta_path, sep='\t', index_col=0, low_memory=False,
                           keep_default_na=False, na_values=['__SCAID_SOURCE_MISSING__'])
    axes = None
    for method, item in report['matrices'].items():
        if item.get('validated'):
            if not Path(item['path']).is_file() or sha256(item['path']) != item['sha256']:
                raise ValueError('Previously promoted file missing or changed')
            continue
        temp = Path(item['temporary_path'])
        if not temp.is_file(): raise ValueError('Missing staging file: '+str(temp))
        print('VALIDATE',report['object_id'],method,flush=True)
        a = ad.read_h5ad(temp,backed='r')
        if tuple(a.shape) != tuple(item['shape']): raise ValueError('AnnData dimensions differ')
        if not a.obs_names.is_unique or not a.var_names.is_unique: raise ValueError('Duplicate axis identifiers')
        if not a.obs_names.equals(original.index): raise ValueError('Source cell identifiers/order differ')
        if axes is not None and not a.obs_names.equals(axes): raise ValueError('Cell axes differ across methods')
        axes = a.obs_names.copy()
        expected_additional=item.get('additional_derived_metadata',[])
        if list(a.obs.columns) != list(original.columns)+expected_additional: raise ValueError('Metadata columns/order differ')
        comparisons = {}
        for key in original:
            source = original[key]
            actual = a.obs[key]
            if isinstance(actual.dtype,pd.CategoricalDtype) or pd.api.types.is_object_dtype(actual):
                s = source.astype('string').fillna('__MISSING__')
                t = actual.astype('string').fillna('__MISSING__')
                if not np.array_equal(s.to_numpy(),t.to_numpy()): raise ValueError('Metadata values differ: '+key)
            else:
                if not np.allclose(pd.to_numeric(source).to_numpy(dtype=float),actual.to_numpy(dtype=float),
                                   equal_nan=True,rtol=1e-12,atol=1e-12):
                    raise ValueError('Numeric metadata differs: '+key)
            comparisons[key] = True
        checks = {'AnnData_read_h5ad': True,'shape':list(a.shape),'unique_cell_ids':True,
                  'unique_feature_ids':True,'full_metadata_comparison':comparisons,
                  'X':check_probes(a.X,item['matrix_checks']['X']['source_value_probes'])}
        with h5py.File(temp,'r') as h:
            x = h['X']
            if isinstance(x,h5py.Group):
                if tuple(x.attrs['shape']) != tuple(a.shape): raise ValueError('CSR orientation mismatch')
                if x['indptr'].shape != (a.n_obs+1,): raise ValueError('CSR row pointer length mismatch')
                nnz = len(x['data'])
                if int(x['indptr'][-1]) != nnz or len(x['indices']) != nnz: raise ValueError('CSR component length mismatch')
                checks['X']['nnz'] = nnz
            if set(a.obsm) != set(h['obsm']): raise ValueError('Embedding keys differ')
            for key in a.obsm:
                if a.obsm[key].shape[0] != a.n_obs or not np.isfinite(a.obsm[key]).all():
                    raise ValueError('Invalid embedding: '+key)
                if not np.array_equal(np.asarray(a.obsm[key]),h['obsm'][key][:]):
                    raise ValueError('Embedding readback differs: '+key)
            checks['embeddings'] = {k:list(a.obsm[k].shape) for k in a.obsm}
        processed_only=method=='RNA' and item.get('counts_available') is False
        if processed_only:
            if report['object_id']!='T1D/GSE279086/Kidney' or 'counts' in a.layers or 'source_processed_counts' not in a.layers:
                raise ValueError('Invalid explicitly processed-only RNA policy')
            if item.get('raw_count_status')!='unavailable_exact_source_axis_mapping_incomplete':
                raise ValueError('Missing explicit raw count scope')
            processed=a.layers['source_processed_counts']
            checks['source_processed_counts']=check_probes(processed,item['matrix_checks']['source_processed_counts']['source_value_probes'])
            if not np.allclose(np.asarray(processed.sum(1)).ravel(),a.obs['nCount_RNA'],rtol=1e-12,atol=1e-12):
                raise ValueError('Original QC is inconsistent with processed counts')
            checks['source_processed_counts']['original_nCount_RNA_matches']=True
            checks['counts_available']=False
            checks['raw_count_status']=item['raw_count_status']
        if method == 'RNA' and not processed_only:
            c = a.layers['counts']
            if tuple(c.shape) != tuple(a.shape): raise ValueError('RNA counts shape differs')
            checks['counts'] = check_probes(c,item['matrix_checks']['counts']['source_value_probes'])
            values = c.data if sparse.issparse(c) else np.asarray(c).ravel()
            for pos in range(0,len(values),2_000_000):
                chunk = values[pos:pos+2_000_000]
                if not np.isfinite(chunk).all() or (chunk<0).any() or not np.equal(chunk,np.rint(chunk)).all():
                    raise ValueError('Invalid nonnegative integer counts')
            checks['counts']['all_values_nonnegative_integer'] = True
            checks['counts']['nnz'] = int(c.nnz) if sparse.issparse(c) else int(np.count_nonzero(c))
            if sparse.issparse(c):
                library = np.asarray(c.sum(axis=1)).ravel()
                detected = np.asarray((c>0).sum(axis=1)).ravel()
                checks['counts']['zero_count_cells'] = int((library==0).sum())
                if 'nCount_RNA' in a.obs:
                    checks['counts']['nCount_RNA_matches'] = bool(np.allclose(library,a.obs['nCount_RNA'],rtol=0,atol=1e-8))
                if 'nFeature_RNA' in a.obs:
                    checks['counts']['nFeature_RNA_matches'] = bool(np.array_equal(detected,a.obs['nFeature_RNA']))
                if expected_additional:
                    if not np.array_equal(library,a.obs['nCount_source_raw']):raise ValueError('Recovered raw count QC mismatch')
                    if not np.array_equal(detected,a.obs['nFeature_source_raw']):raise ValueError('Recovered raw feature QC mismatch')
                    checks['counts']['nCount_source_raw_matches']=True
                    checks['counts']['nFeature_source_raw_matches']=True
            if 'source_processed_counts' in a.layers:
                processed=a.layers['source_processed_counts']
                checks['source_processed_counts']=check_probes(processed,item['matrix_checks']['source_processed_counts']['source_value_probes'])
                if not np.allclose(np.asarray(processed.sum(1)).ravel(),a.obs['nCount_RNA'],rtol=1e-12,atol=1e-12):
                    raise ValueError('Original QC is inconsistent with preserved processed counts')
                checks['source_processed_counts']['original_nCount_RNA_matches']=True
        a.file.close()
        del a
        gc.collect()
        digest = sha256(temp)
        # Files are promoted only after their independent reader and source probes pass.
        final = Path(item['path'])
        if final.exists(): raise FileExistsError('Refusing to overwrite: '+str(final))
        os.replace(temp,final)
        item.update(validated=True,validation_status='passed',bytes=final.stat().st_size,
                    sha256=digest,independent_validation=checks,format='AnnData H5AD')
        pending = report_path.with_suffix('.json.tmp')
        pending.write_text(json.dumps(report,ensure_ascii=False,indent=2))
        pending.replace(report_path)
        print('PASSED',method,item['bytes'],digest,flush=True)
    report['source_sha256'] = sha256(report['source_rds'])
    report['metadata_tsv_gz'] = str(meta_path)
    report['metadata_tsv_gz_sha256'] = sha256(meta_path)
    report['status'] = 'validated'
    report['validation_seconds'] = round(time.time()-started,3)
    report['validation_versions'] = {'anndata':version('anndata'),'h5py':h5py.__version__,
                                    'numpy':np.__version__,'pandas':pd.__version__}
    pending = report_path.with_suffix('.json.tmp')
    pending.write_text(json.dumps(report,ensure_ascii=False,indent=2))
    pending.replace(report_path)
    return report

if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('report_path');args=p.parse_args()
    result=validate(args.report_path)
    print(json.dumps({'object_id':result['object_id'],'status':result['status'],
                      'files':len(result['matrices'])},ensure_ascii=False),flush=True)

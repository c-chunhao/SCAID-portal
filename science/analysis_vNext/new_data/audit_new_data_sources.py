"""Full sparse numeric/schema audit; eligibility is distinct from integer values.

No matrix is densified. MatrixMarket data are streamed in bounded blocks.
Clinical identities/captures are never inferred from barcode suffixes or titles.
"""
from pathlib import Path
import argparse, collections, concurrent.futures, csv, gzip, hashlib, json, os, subprocess
import h5py
import numpy as np
import pandas as pd

NATIVE_MTX=None

def checksum(path):
    with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def valid_vector(v):
    return bool(np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.floor(v)).all())

def sparse_check(m,shape,csc=False):
    assert len(shape)==2 and min(shape)>0
    data,indices,ptr=m['data'],m['indices'],m['indptr']
    assert len(data)==len(indices)
    axis=1 if csc else 0;bound=int(shape[1-axis])
    assert len(ptr)==int(shape[axis])+1
    last=0
    for start in range(0,len(ptr),1000000):
        p=ptr[start:start+1000000]
        assert p[0]>=last and np.all(p[1:]>=p[:-1]);last=int(p[-1])
    assert ptr[0]==0 and last==len(data)
    numeric=True
    for start in range(0,len(data),1000000):
        x=data[start:start+1000000];ix=indices[start:start+1000000]
        assert np.all((ix>=0)&(ix<bound));numeric &= valid_vector(x)
    return numeric,len(data)

def strings(dataset):
    return dataset.asstr()[:].astype(str)

def h5_audit(path):
    with h5py.File(path) as f:
        if 'matrix' in f and 'shape' in f['matrix']:
            m=f['matrix'];shape=m['shape'][:].astype(int);integer,nnz=sparse_check(m,shape,csc=True)
            assert len(m['barcodes'])==shape[1] and len(m['features/id'])==shape[0]
            assert len(set(strings(m['barcodes'])))==shape[1]
            types=strings(m['features/feature_type']) if 'features/feature_type' in m else np.array([])
            return dict(format='10x_HDF5',cells=int(shape[1]),features=int(shape[0]),nnz=nnz,
                        full_numeric='integer_nonnegative' if integer else 'noninteger_or_invalid',
                        GEX_features=int((types=='Gene Expression').sum()),schema_pass=True,
                        provenance='source_10x_candidate; requires channel/GSM/clinical mapping')
        loc='layers/counts' if 'layers/counts' in f else 'raw/X' if 'raw/X' in f else 'X'
        if loc not in f:return dict(format='HDF5_other',schema_pass=False,full_numeric='no_expression_schema')
        m=f[loc];shape=m.attrs.get('shape',m.shape if isinstance(m,h5py.Dataset) else [])
        if isinstance(m,h5py.Group):
            enc=m.attrs.get('encoding-type','');enc=enc.decode() if isinstance(enc,bytes) else str(enc)
            assert enc in ['csr_matrix','csc_matrix'],f'Unknown sparse encoding {enc}'
            integer,nnz=sparse_check(m,shape,enc=='csc_matrix')
        else:
            assert len(shape)==2;integer=True;nnz=0
            block=max(1,1000000//int(shape[1]))
            for start in range(0,shape[0],block):
                x=m[start:start+block];integer &= valid_vector(x);nnz+=int(np.count_nonzero(x))
        for axis,name in enumerate(['obs','raw/var' if loc=='raw/X' else 'var']):
            g=f[name];key=g.attrs.get('_index','_index');key=key.decode() if isinstance(key,bytes) else key
            ids=strings(g[key]);assert len(ids)==shape[axis] and len(set(ids))==len(ids)
        return dict(format='AnnData_h5ad',cells=int(shape[0]),features=int(shape[1]),nnz=nnz,
                    expression_slot=loc,full_numeric='integer_nonnegative' if integer else 'noninteger_or_invalid',
                    schema_pass=True,provenance='integer representation does not establish original raw-count provenance')

def companion(path,kind):
    stem=path.name.split('matrix.mtx')[0]
    choices=[p for p in path.parent.iterdir() if p.is_file() and p.name.startswith(stem)
             and kind in p.name.lower() and p.name.endswith(('.tsv','.tsv.gz','.txt','.txt.gz'))]
    if len(choices)!=1:return None
    p=choices[0];op=gzip.open if p.name.endswith('.gz') else open
    with op(p,'rt') as f:rows=[line.rstrip('\n').split('\t') for line in f]
    return p,rows

def mtx_audit(path):
    if NATIVE_MTX:
        result=subprocess.run([str(NATIVE_MTX),str(path)],capture_output=True,text=True,check=False)
        if result.returncode:raise ValueError(result.stderr.strip()[:300])
        full=json.loads(result.stdout);assert full['schema_pass']
        nr,nc,nnz=full['features'],full['cells'],full['nnz']
        integer=full['full_numeric']=='integer_nonnegative'
    else:
        op=gzip.open if path.name.endswith('.gz') else open
        n=0;integer=True
        with op(path,'rt') as f:
            header=f.readline().strip();assert header.lower().startswith('%%matrixmarket matrix coordinate')
            assert ' general' in header.lower() and not any(x in header.lower() for x in ['complex','pattern'])
            line=f.readline()
            while line.startswith('%'):line=f.readline()
            nr,nc,nnz=map(int,line.split());assert min(nr,nc)>0 and nnz>=0
            while True:
                lines=f.readlines(4*1024*1024)
                if not lines:break
                x=np.fromstring(''.join(lines),sep=' ')
                assert len(x)%3==0,'Malformed MatrixMarket entry';x=x.reshape(-1,3)
                assert np.equal(x[:,:2],np.floor(x[:,:2])).all()
                assert ((x[:,0]>=1)&(x[:,0]<=nr)&(x[:,1]>=1)&(x[:,1]<=nc)).all()
                integer &= valid_vector(x[:,2]);n+=len(x)
        assert n==nnz,f'Declared {nnz}, observed {n}'
    bar=companion(path,'barcode');features=companion(path,'feature') or companion(path,'gene')
    if bar:assert len(bar[1])==nc and len({r[0] for r in bar[1]})==nc
    if features:assert len(features[1])==nr
    modality='unresolved_features';gex=''
    if features and all(len(r)>=3 for r in features[1]):
        kinds=collections.Counter(r[2] for r in features[1]);gex=kinds.get('Gene Expression',0)
        modality=json.dumps(dict(kinds))
    return dict(format='MatrixMarket',cells=nc,features=nr,nnz=nnz,full_numeric='integer_nonnegative' if integer else 'noninteger_or_invalid',
                schema_pass=True,barcode_companion_verified=bool(bar),feature_companion_verified=bool(features),
                GEX_features=gex,feature_modalities=modality,
                provenance='full numeric/schema validation; source counts and clinical/capture mapping still required')

def inspect(item):
    item=dict(item);p=Path(item['path']);name=p.name.lower()
    try:
        if name.endswith(('.rds','.rds.gz','.rdata','.rdata.gz','.rda','.rda.gz')):
            item.update(format='R_serialized_object',full_numeric='counts_slot_not_verified',schema_pass=False)
        elif name.endswith(('.h5','.hdf5','.h5ad')):item.update(h5_audit(p))
        elif name.endswith(('.mtx','.mtx.gz')):item.update(mtx_audit(p))
        else:return item
        item['sha256']=checksum(p)
    except (OSError,ValueError,AssertionError,KeyError) as e:
        item.update(schema_pass=False,full_numeric='audit_failed',error=f'{type(e).__name__}: {str(e)[:300]}')
    return item

def run(args):
    global NATIVE_MTX
    NATIVE_MTX=args.native_mtx
    if NATIVE_MTX:assert NATIVE_MTX.is_file() and os.access(NATIVE_MTX,os.X_OK)
    args.out.mkdir(parents=True,exist_ok=True)
    old=pd.read_csv(args.inventory,sep='\t').fillna('')
    # Inventory discovery covers gzipped RData previously incorrectly ancillary.
    rows=old.to_dict('records')
    jobs=[r for r in rows if str(r['path']).lower().endswith(('.mtx','.mtx.gz','.h5','.hdf5','.h5ad','.rds','.rds.gz','.rdata','.rdata.gz','.rda','.rda.gz'))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results=[]
        for i,r in enumerate(pool.map(inspect,jobs)):
            results.append(r)
            if i%10==0:print('Audited',i+1,'/',len(jobs),r['accession'],flush=True)
    replacements={r['path']:r for r in results};rows=[replacements.get(r['path'],r) for r in rows]
    frame=pd.DataFrame(rows)
    # Keep discovery-era sampled evidence separate from the complete audit.
    frame['legacy_numeric_validation']=frame['numeric_validation']
    inspected=frame['full_numeric'].notna()
    frame.loc[inspected,'numeric_validation']=frame.loc[inspected,'full_numeric'].map({
        'integer_nonnegative':'full_expression_vector_integer_nonnegative',
        'noninteger_or_invalid':'full_expression_vector_noninteger_or_invalid',
        'counts_slot_not_verified':'serialized_counts_slot_not_verified',
        'audit_failed':'audit_failed'})
    frame.to_csv(args.out/'full_file_representation_audit.tsv',sep='\t',index=False)
    registry=pd.read_csv(args.registry,sep='\t');registry=registry[registry.classification.eq('geo_public_candidate')]
    published=json.loads(args.assets.read_text())['datasets'];ids={r['dataset_id'].split(':')[-1] for r in published}
    datasets=[]
    for _,r in registry.iterrows():
        acc=r.accession;files=[f for f in results if f['accession']==acc]
        valid=[f for f in files if f.get('schema_pass') and f.get('full_numeric')=='integer_nonnegative']
        datasets.append({'accession':acc,'dataset_name':r.dataset_name,'local_expression_or_R_files':len(files),
                         'full_integer_schema_pass_files':len(valid),'failed_or_unverified_files':len(files)-len(valid),
                         'already_released_accession':acc in ids,
                         'novelty_gate':'sample/barcode/source duplicate audit before extension' if acc in ids else 'new accession; participant/capture mapping still required',
                         'release_eligible':False,'clinical_capture_gate':'not automatically verified by numeric audit',
                         'readiness':'APS_processed_QC_annotation_staging' if acc=='GSE262240' else 'full_sparse_candidates_require_source_clinical_capture_validation' if valid else 'no_verified_full_integer_expression',
                         'date_gate':'GEO release, publication, and download dates must be recorded separately'})
    pd.DataFrame(datasets).to_csv(args.out/'59_candidate_release_gates.tsv',sep='\t',index=False)
    report={'candidates':len(datasets),'expression_or_R_files':len(jobs),
            'full_integer_schema_pass_files':sum(bool(r.get('schema_pass')) and r.get('full_numeric')=='integer_nonnegative' for r in results),
            'accession_overlap_with_release':[r['accession'] for r in datasets if r['already_released_accession']],
            'numeric_or_schema_errors':[{'path':r['path'],'error':r['error']} for r in results if r.get('error')],
            'automatically_released':0,'bounded_stream_block_bytes':4194304,
            'counts_provenance_is_not_inferred_from_integral_values':True,
            'MatrixMarket_engine':'full native gzip stream' if NATIVE_MTX else 'full bounded Python/numpy stream',
            'native_scanner_sha256':checksum(NATIVE_MTX) if NATIVE_MTX else None,
            'RData_gz_recognition_fixed':True,'clinical_and_technical_capture_labels_not_inferred_from_suffix':True}
    (args.out/'new_data_audit_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inventory',type=Path,required=True);p.add_argument('--registry',type=Path,required=True)
    p.add_argument('--assets',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=4);p.add_argument('--native-mtx',type=Path)
    run(p.parse_args())

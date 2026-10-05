#!/usr/bin/env python3
"""Regression checks for observed missing-singleton count export failures."""
from pathlib import Path
import tempfile,json
import numpy as np,h5py
from scipy import sparse
from sparse_export_guard import write_raw_blocks,validate_raw_blocks

results=[]
for n in [0,1,4999,5000,5001,10000,10001]:
    x=sparse.csr_matrix((np.arange(1,n+1,dtype='int64'),(np.arange(n),np.arange(n)%17)),shape=(n,17))
    with tempfile.NamedTemporaryFile(suffix='.h5')as tmp:
        with h5py.File(tmp.name,'w')as h:write_raw_blocks(h.create_group('rawdata'),x)
        with h5py.File(tmp.name)as h:
            parts=[]
            for z in h['rawdata'].values():parts.append(sparse.csr_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=tuple(z.attrs['shape'])))
            y=sparse.vstack(parts,format='csr')if parts else sparse.csr_matrix((0,17))
            assert y.shape==x.shape and (x-y).nnz==0
            results.append({'cells':n,'blocks':len(parts),'all_values_exact':True,'tail_preserved':True})
        if n==5001:
            with h5py.File(tmp.name,'r+')as broken:
                del broken['rawdata/rawdata_00001.npz']
                try:validate_raw_blocks(broken['rawdata'],n,17)
                except AssertionError as e:results.append({'dropped_singleton_tail_detected':True,'error':str(e)})
                else:raise AssertionError('Missing final cell not detected')
Path(__file__).with_name('outputs').mkdir(parents=True,exist_ok=True)
report_dir=Path(__file__).with_name('outputs');report_dir.mkdir(parents=True,exist_ok=True)
report_dir.joinpath('sparse_chunk_boundary_validation.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))

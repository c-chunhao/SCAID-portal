"""Boundary-safe rawdata block writer/reader for the project's released H5 format."""
import numpy as np
import h5py
from scipy import sparse

def write_raw_blocks(group, counts_cells_by_genes, block_size=5000):
    """Write sparse cells×genes, retaining every singleton tail and zero row."""
    assert sparse.issparse(counts_cells_by_genes)
    n,g=counts_cells_by_genes.shape
    assert block_size>0
    for i,start in enumerate(range(0,n,block_size)):
        stop=min(start+block_size,n)
        x=counts_cells_by_genes[start:stop,:].tocsr()
        x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
        z=group.create_group(f'rawdata_{i:05d}.npz')
        z.attrs.update({'encoding-type':'csr_matrix','encoding-version':'0.1.0','shape':np.asarray(x.shape,dtype=np.int64)})
        z.create_dataset('data',data=x.data,compression='gzip')
        z.create_dataset('indices',data=x.indices.astype(np.int32),compression='gzip')
        z.create_dataset('indptr',data=x.indptr.astype(np.int64),compression='gzip')
    validate_raw_blocks(group,n,g)

def validate_raw_blocks(group,n_cells,n_genes):
    offset=0
    for i,key in enumerate(sorted(group)):
        assert key==f'rawdata_{i:05d}.npz','Non-contiguous rawdata block names'
        z=group[key]
        assert isinstance(z,h5py.Group),'RNA rawdata cannot be a pathway/dense matrix'
        shape=tuple(int(a)for a in z.attrs['shape'])
        ptr=z['indptr'][:];ix=z['indices'][:]
        assert len(shape)==2 and shape[1]==n_genes and shape[0]>0
        assert len(ptr)==shape[0]+1 and ptr[0]==0 and ptr[-1]==len(ix)==len(z['data'])
        assert np.all(np.diff(ptr)>=0) and np.all((ix>=0)&(ix<n_genes))
        offset+=shape[0]
    assert offset==n_cells,f'Rawdata has {offset} rows for {n_cells} declared cells'
    return offset

def assert_expression_contract(h5):
    """Reject malformed RNA assets before loading cells, including missing tails."""
    required=['names_obs','names_var','obs','assay/RNA/layers/rawdata']
    for name in required:assert name in h5,f'RNA asset missing {name}'
    return validate_raw_blocks(h5['assay/RNA/layers/rawdata'],len(h5['names_obs']),len(h5['names_var']))

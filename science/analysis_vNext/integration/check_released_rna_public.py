"""Self-contained checks for released RNA axis/normalization/corruption guards."""
from pathlib import Path
import tempfile,json
import h5py,numpy as np
from scipy import sparse
from released_rna import iter_raw_counts,iter_full_log1p,read_dense_hvg,inspect_schema,validate_raw_csr_arrays

def run():
 checks={}
 with tempfile.TemporaryDirectory() as td:
  h5=Path(td)/'toy.h5';raw=np.array([[1,2,7],[2,6,2]],dtype=float);x=sparse.csr_matrix(raw)
  with h5py.File(h5,'w') as h:
   for key,ids in [('names_obs',['c1','c2']),('names_var',['g1','g2','g3']),('var/rawvar/_index',['g1','g2','g3']),('var/var/_index',['g3','g1'])]:h.create_dataset(key,data=np.array(ids,dtype=object),dtype=h5py.string_dtype())
   g=h.create_group('assay/RNA/layers/rawdata/rawdata_00000.npz');g.attrs['shape']=x.shape
   for k,v in [('data',x.data),('indices',x.indices),('indptr',x.indptr)]:g.create_dataset(k,data=v)
   z=h.create_dataset('assay/RNA/layers/data/data_00000.npz',data=np.array([[30,10],[31,11]],dtype=float));z.attrs['shape']=[2,2]
  ids,genes,c=next(iter_raw_counts(h5));assert np.array_equal(c.toarray(),raw)
  ids,genes,log=next(iter_full_log1p(h5));assert np.allclose(np.expm1(log.toarray()).sum(axis=1),10000)
  ids,genes,sub=next(iter_full_log1p(h5,features=['g1']));assert np.allclose(np.expm1(sub.toarray()).ravel(),[1000,2000])
  checks.update(raw_and_dense_feature_axes_distinguished=True,full_library_denominator_preserved_after_subset=True)
  try:read_dense_hvg(h5,[0],['g1']);raise AssertionError('Dense must fail closed')
  except ValueError:checks['unverified_dense_axis_rejected']=True
  z=read_dense_hvg(h5,[1,0],['g1','g3'],verified_source_feature_names=['g3','g1']);assert np.array_equal(z,[[11,31],[10,30]])
  checks['verified_axis_and_permuted_cells_read_correctly']=True
  for key,args in {
   'large_fractional_counts_rejected':(np.array([100000.1]),np.array([0]),np.array([0,1]),(1,1)),
   'bad_pointer_rejected':(np.array([1]),np.array([0]),np.array([1,1]),(1,1)),
   'bad_index_rejected':(np.array([1]),np.array([1]),np.array([0,1]),(1,1)),
   'duplicate_CSR_coordinate_rejected':(np.array([1,2]),np.array([0,0]),np.array([0,2]),(1,1)),
  }.items():
   try:validate_raw_csr_arrays(*args);raise AssertionError(key)
   except ValueError:checks[key]=True
  with h5py.File(h5,'r+') as h:h['var/rawvar/_index'][0]='wrong'
  try:next(iter_raw_counts(h5));raise AssertionError('Raw axis mismatch accepted')
  except ValueError:checks['raw_axis_mismatch_rejected']=True
 print(json.dumps(checks,indent=2));return checks
if __name__=='__main__':run()

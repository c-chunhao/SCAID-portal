"""Explicit, bounded readers for the released SCAID multi-feature HDF5 schema.

Raw RNA and dense HVG matrices have DIFFERENT feature axes. Each reader checks
the matching feature index. A normalized annotation matrix is computed from
full raw counts; the historically named dense ``data`` layer is never silently
treated as full log-normalized RNA or as raw counts.
"""
from pathlib import Path
import argparse, json
import h5py, numpy as np
from scipy import sparse

SCHEMA_VERSION='SCAID-RNA-sidecar-1.0'

def _strings(node): return np.asarray(node.asstr()[:],dtype=str)

def validate_raw_csr_arrays(data,indices,indptr,shape):
 """Strict numerical/structural count checks before constructing a CSR."""
 data=np.asarray(data);indices=np.asarray(indices);indptr=np.asarray(indptr)
 if len(shape)!=2 or any(int(v)!=v or v<0 for v in shape):raise ValueError('Invalid raw CSR shape')
 rows,cols=map(int,shape)
 if data.ndim!=1 or indices.ndim!=1 or indptr.ndim!=1:raise ValueError('CSR components must be 1D')
 if indices.dtype.kind not in 'iu' or indptr.dtype.kind not in 'iu':raise ValueError('CSR indices and pointers must have integer dtype')
 if len(indptr)!=rows+1 or len(data)!=len(indices):raise ValueError('CSR component lengths disagree')
 if indptr[0]!=0 or indptr[-1]!=len(data) or np.any(indptr[1:]<indptr[:-1]):raise ValueError('CSR row pointers must start at zero, end at nnz, and be monotone')
 if np.any(indices<0) or np.any(indices>=cols):raise ValueError('CSR feature index out of bounds')
 if data.dtype.kind not in 'iuf' or not np.isfinite(data).all() or np.any(data<0) or not np.equal(data,np.floor(data)).all():raise ValueError('Raw RNA counts must be finite, nonnegative, exactly integral values')
 x=sparse.csr_matrix((data,indices,indptr),shape=(rows,cols))
 x.sort_indices()
 if not x.has_canonical_format:raise ValueError('Repeated feature coordinates in a raw CSR row')
 return x

def inspect_schema(path):
 with h5py.File(path,'r') as h:
  ids=_strings(h['names_obs']);full=_strings(h['names_var'])
  if len(np.unique(ids))!=len(ids):raise ValueError('Duplicate cell IDs')
  if len(np.unique(full))!=len(full):raise ValueError('Duplicate full raw feature IDs')
  if 'var/rawvar/_index' in h and not np.array_equal(full,_strings(h['var/rawvar/_index'])):raise ValueError('names_var and rawvar feature index disagree')
  raw=h['assay/RNA/layers/rawdata'];blocks=[];n=0
  for k in sorted(raw):
   g=raw[k];shape=tuple(int(x) for x in g.attrs['shape'])
   if len(shape)!=2 or shape[1]!=len(full):raise ValueError(f'{k}: raw CSR features do not match full raw axis')
   if g['indptr'].shape!=(shape[0]+1,):raise ValueError(f'{k}: CSR row pointers malformed')
   if g['data'].shape!=g['indices'].shape:raise ValueError(f'{k}: CSR data/index lengths differ')
   blocks.append({'node':g.name,'shape_cells_by_features':list(shape),'cell_offset':n,'feature_index':'names_var','semantics':'unintegrated raw RNA counts','encoding':'csr_matrix'})
   n+=shape[0]
  if n!=len(ids):raise ValueError('Raw block rows do not cover the cell axis exactly')
  result={'schema_version':SCHEMA_VERSION,'source_file':str(Path(path).resolve()),'cell_index':'names_obs','cells':len(ids),'raw_features':len(full),'raw_feature_index':'names_var','raw_feature_metadata_index':'var/rawvar/_index' if 'var/rawvar/_index' in h else None,'raw_blocks':blocks,'normalization_policy':'sum all raw genes per cell, multiply by 10000, log1p, THEN optional feature subset','hvg_layers':[],'legacy_dense_feature_order_independently_verified':False,'legacy_dense_default_read_policy':'fail closed until corresponding source scale.data row names and export association are verified; raw counts remain readable'}
  if 'assay/RNA/layers/data' in h:
   layer=h['assay/RNA/layers/data']
   if 'var/var/_index' not in h:raise ValueError('Dense/HVG layer has no explicit var/var/_index; cannot infer order from full names_var')
   hvgs=_strings(h['var/var/_index'])
   if len(np.unique(hvgs))!=len(hvgs):raise ValueError('Duplicate HVG feature IDs')
   # Some historical files renamed raw RNA symbols without changing the
   # stored HVG names. Keep both explicit axes; do not guess alias matches or
   # prevent valid raw-count reads because an auxiliary layer needs repair.
   unmatched_hvg=hvgs[~np.isin(hvgs,full)].tolist()
   result['hvg_feature_names_absent_from_raw_axis']=unmatched_hvg
   result['cross_layer_gene_alias_reconciliation_required']=bool(unmatched_hvg)
   for k in sorted(layer):
    z=layer[k]
    if not isinstance(z,h5py.Dataset):raise ValueError('Only explicitly mapped dense HVG blocks are accepted by this reader')
    shape=z.shape
    if len(shape)!=2 or shape[1]!=len(hvgs):raise ValueError(f'{k}: actual dense feature axis does not match HVG index')
    saved=[int(v) for v in np.asarray(z.attrs.get('shape',shape),dtype=int)]
    example=z[:min(12,shape[0]),:min(20,shape[1])]
    result['hvg_layers'].append({'node':z.name,'actual_shape_cells_by_features':list(shape),'legacy_shape_attribute':saved,'legacy_attribute_reversed':saved==list(shape[::-1]) and saved!=list(shape),'feature_index':'var/var/_index','features':len(hvgs),'semantics':'historical dense HVG data; inspect provenance before assuming scaling/normalization','sample_min':float(example.min()),'sample_max':float(example.max()),'has_negative_sample_values':bool((example<0).any()),'valid_for_raw_count_models':False,'valid_for_full_RNA_CellTypist_input':False})
  return result

def iter_raw_counts(path):
 """Yield (cell_ids, full_gene_ids, counts CSR) with all axis checks."""
 inspect_schema(path)
 with h5py.File(path,'r') as h:
  ids=_strings(h['names_obs']);genes=_strings(h['names_var']);offset=0
  for k in sorted(h['assay/RNA/layers/rawdata']):
   g=h['assay/RNA/layers/rawdata'][k];shape=tuple(int(v) for v in g.attrs['shape'])
   x=validate_raw_csr_arrays(g['data'][:],g['indices'][:],g['indptr'][:],shape)
   yield ids[offset:offset+shape[0]],genes,x
   offset+=shape[0]

def iter_full_log1p(path,features=None,target_sum=10000):
 """Normalize the full library, then restrict genes in requested order."""
 if not np.isfinite(target_sum) or target_sum<=0:raise ValueError('Normalization target_sum must be finite and positive')
 if features is not None and len(set(features))!=len(features):raise ValueError('Requested features must be unique')
 for ids,genes,counts in iter_raw_counts(path):
  total=np.asarray(counts.sum(axis=1)).ravel()
  if (total<=0).any():raise ValueError('Zero raw library; quarantine before annotation')
  x=counts.astype(np.float32).multiply((target_sum/total)[:,None]).tocsr();x.data=np.log1p(x.data)
  if features is not None:
   lookup={g:i for i,g in enumerate(genes)}
   missing=[g for g in features if g not in lookup]
   if missing:raise KeyError(f'Requested full RNA features absent: {missing[:10]}')
   idx=[lookup[g] for g in features];x=x[:,idx];genes=np.asarray(features)
  yield ids,genes,x

def read_dense_hvg(path,cell_rows,features,verified_source_feature_names=None):
 """Read using independently verified source scale.data row names.

 A legacy axis-length match or correlation fingerprint alone does not authorize
 this read. The caller must provide the matching export source's scale.data
 row names, not merely repeat the H5 var index. New writers should save that
 exact axis with explicit source provenance.
 """
 if verified_source_feature_names is None:
  raise ValueError('Legacy dense/HVG feature order is not independently verified. Use raw RNA counts, or provide the corresponding export-source scale.data row names with verified export association.')
 schema=inspect_schema(path)
 if len(schema['hvg_layers'])!=1:raise ValueError('Explicit block mapping required for multiple dense HVG blocks')
 if len(cell_rows)>10000:raise ValueError('Bounded dense reader accepts at most 10000 selected rows')
 with h5py.File(path,'r') as h:
  declared_hvgs=_strings(h['var/var/_index'])
  hvgs=np.asarray(verified_source_feature_names,dtype=str)
  if len(hvgs)!=len(declared_hvgs) or len(np.unique(hvgs))!=len(hvgs):raise ValueError('Verified source feature axis length/uniqueness mismatch')
  lookup={g:i for i,g in enumerate(hvgs)}
  missing=[g for g in features if g not in lookup]
  if missing:raise KeyError(f'Genes absent from dense HVG axis (they may exist in raw RNA): {missing[:10]}')
  original_rows=np.asarray(cell_rows)
  if original_rows.ndim!=1 or original_rows.dtype.kind not in 'iuf' or not np.isfinite(original_rows).all() or not np.equal(original_rows,np.floor(original_rows)).all():raise ValueError('Requested cell rows must be finite integral indices')
  rows=original_rows.astype(int)
  if (rows<0).any() or (rows>=schema['cells']).any():raise IndexError('Requested cell rows are outside the released cell axis')
  if len(np.unique(rows))!=len(rows):raise ValueError('Duplicate requested cell rows')
  sort=np.argsort(rows);inverse=np.argsort(sort)
  z=h[schema['hvg_layers'][0]['node']][rows[sort],:][inverse,:]
  return z[:,[lookup[g] for g in features]]

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('h5');p.add_argument('--out',required=True);a=p.parse_args()
 Path(a.out).write_text(json.dumps(inspect_schema(a.h5),indent=2,ensure_ascii=False)+'\n')

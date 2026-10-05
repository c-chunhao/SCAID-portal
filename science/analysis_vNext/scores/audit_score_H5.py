"""Stream every archived score matrix; distinguish scores from RNA counts.
Do not mutate original assets or impute missing cells/pathways as zero.
"""
from pathlib import Path
import sys, os
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime_paths import required_path, registered_datasets
import json,hashlib
import h5py
import numpy as np
import pandas as pd
R=required_path('SCAID_SCORE_AUDIT_DIR');T=R/'tables'
datasets=registered_datasets();rows=[]
prior=T/'Score_H5_144_full_matrix_QA.tsv'
if prior.exists():rows=pd.read_csv(prior,sep='\t').to_dict('records')
completed={z['path'] for z in rows}
def strings(d):return d.asstr()[:]
def shastr(a):return hashlib.sha256('\0'.join(a).encode()).hexdigest()
for i,item in enumerate(datasets):
 rnap=Path(item['h5_paths']['RNA']);obj=rnap.stem
 with h5py.File(rnap,'r') as h:
  obs=strings(h['names_obs']);genes=strings(h['names_var'])
  assert len(obs)==len(set(obs)) and len(genes)==len(set(genes))
 ref=None
 for method in ['AUCell','UCell','singscore']:
  f=Path(item['h5_paths'][method]);print(i+1,method,obj,flush=True)
  if str(f) in completed:
   with h5py.File(f,'r') as h:
    pp=strings(h['matrix/pathway'])
    if ref is None:ref=pp.copy()
   continue
  with h5py.File(f,'r') as h:
   assert set(h.keys())=={'matrix'}
   g=h['matrix'];cell=strings(g['cell']);pathways=strings(g['pathway']);ind=g['x'];ptr=np.asarray(g['y']);values=g['value']
   ncell=len(cell);npath=len(pathways);nnz=len(values);pointers_ok=(len(ptr)==ncell+1 and ptr[0]==0 and ptr[-1]==nnz and np.all(np.diff(ptr)>=0))
   assert len(ind)==nnz and pointers_ok
   finite=True;inbounds=True;vmin=0.;vmax=0.;nonfinite=0;neg=0;high=0;integerlike=True
   for start in range(0,nnz,2_000_000):
    v=np.asarray(values[start:start+2_000_000]);ix=np.asarray(ind[start:start+2_000_000])
    z=np.isfinite(v);finite &= z.all();nonfinite+=int((~z).sum());inbounds &= bool(((ix>=0)&(ix<npath)).all())
    if z.any():vmin=min(vmin,float(v[z].min()));vmax=max(vmax,float(v[z].max()));integerlike &= bool(np.allclose(v[z],np.rint(v[z])))
    neg+=int((v<0).sum());high+=int((v>1+1e-10).sum())
   # Duplicate indices in one cell would cause scipy CSR summation and change
   # the displayed score; test per-row canonical uniqueness/sort in chunks.
   duplicates=0;unsorted=0
   for start in range(0,ncell,4096):
    end=min(ncell,start+4096);x=np.asarray(ind[ptr[start]:ptr[end]])
    if len(x)<2:continue
    diff=np.diff(x);boundary=ptr[start+1:end]-ptr[start]-1
    boundary=boundary[(boundary>=0)&(boundary<len(diff))];diff[boundary]=1
    bad_positions=np.flatnonzero(diff<=0)
    if len(bad_positions):
     bad_rows=np.unique(np.searchsorted(ptr[start:end+1]-ptr[start],bad_positions,side='right')-1)+start
     for j in bad_rows:
      ix=x[ptr[j]-ptr[start]:ptr[j+1]-ptr[start]]
      duplicates+=int(len(ix)!=len(np.unique(ix)));unsorted+=int(len(ix)>1 and np.any(np.diff(ix)<0))
   aligned=np.array_equal(cell,obs);setmatch=len(cell)==len(obs) and set(cell)==set(obs)
   feature_equal=ref is None or np.array_equal(pathways,ref)
   if ref is None:ref=pathways.copy()
   schema_pass=finite and inbounds and pointers_ok and duplicates==0 and len(cell)==len(set(cell)) and len(pathways)==len(set(pathways))
   expected_range=(-.5,.5) if method=='singscore' else (0.,1.)
   range_ok=vmin>=expected_range[0]-1e-8 and vmax<=expected_range[1]+1e-8
   rows.append({'object':obj,'method':method,'path':str(f),'n_cells':ncell,'n_pathways':npath,'nnz':nnz,'shape_semantics':'cell x pathway CSR; value is a rank/enrichment score, not RNA counts','RNA_n_obs':len(obs),'RNA_obs_order_equal':aligned,'RNA_obs_set_equal':setmatch,'pathway_order_equal_across_3_methods':feature_equal,'pointers_valid':pointers_ok,'index_bounds_valid':inbounds,'n_rows_duplicate_indices':duplicates,'n_rows_unsorted_indices':unsorted,'n_nonfinite_stored':nonfinite,'stored_and_implicit_zero_min':vmin,'stored_and_implicit_zero_max':vmax,'expected_range':str(expected_range),'range_gate_pass':range_ok,'schema_gate_pass':schema_pass,'n_explicit_negative':neg,'n_values_gt1':high,'all_values_integer_like':integerlike,'obs_sha256':shastr(cell),'pathway_order_sha256':shastr(pathways),'RNA_gene_universe_sha256':shastr(genes),'gene_set_provenance_gate_pass':False,'gene_set_gate_reason':'exact frozen KEGG membership/hash and score generation provenance not embedded in score H5; frozen source list recovery required'})
  pd.DataFrame(rows).to_csv(T/'Score_H5_144_full_matrix_QA.tsv',sep='\t',index=False)
x=pd.DataFrame(rows);assert len(x)==144
# Resolve method-specific score semantics before declaring any range defect.
# In particular irGSEA uses uncentred singscore and minimum-rank ties.
import runpy
runpy.run_path(str(R/'scripts/resolve_score_method_semantics.py'),run_name='__main__')
print('144 score matrix audit completed.',flush=True)

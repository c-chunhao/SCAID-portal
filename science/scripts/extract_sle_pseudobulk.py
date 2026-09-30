import os
"""Read the archived sparse source MTX once; write donor/cell-type aggregates only.
Select one sequencing batch per verified source individual before target labels/genes
are tested. Largest total source cell count in cohorts 2/3/4; lexical batch tie break.
"""
from pathlib import Path
import gzip,json,time,hashlib
import numpy as np,pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
from threadpoolctl import threadpool_limits
B=Path(__file__).resolve().parents[1];T=B/'tables';S=Path(os.environ['SCAID_SLE_SOURCE_DIR'])
T.mkdir(exist_ok=True)
IFN=['ISG15','IFI6','IFI44L','MX1','OAS1','IFIT1','IFIT3','RSAD2','ISG20','XAF1','LY6E','IFI44']
print('metadata',flush=True)
m=pd.read_csv(S/'GSE174188_metadata.csv');m=m.rename(columns={m.columns[0]:'cell_id'})
m['condition']=np.where(m.SLE_status.eq('Healthy'),'HD','SLE');m['cohort']=m.Processing_Cohort.astype(int)
assert m.groupby('ind_cov').condition.nunique().max()==1
invariant={c:int((m.groupby('ind_cov')[c].nunique(dropna=True)>1).sum()) for c in ['Sex','pop_cov','Age','Status']}
unit=m.groupby(['ind_cov','batch_cov','cohort','condition'],observed=True).agg(n_cells=('cell_id','size'),Age=('Age','median'),Sex=('Sex','first'),pop_cov=('pop_cov','first'),Status=('Status',lambda s:';'.join(sorted(s.dropna().unique())))).reset_index()
unit.to_csv(T/'GSE174188_all_donor_batch_inventory.tsv',sep='\t',index=False)
eligible=unit[unit.cohort.isin([2,3,4])].copy()
chosen=eligible.sort_values(['ind_cov','n_cells','batch_cov'],ascending=[True,False,True]).drop_duplicates('ind_cov').copy()
assert chosen.ind_cov.is_unique
chosen['unit_id']=chosen.ind_cov+'|'+chosen.batch_cov
chosen.to_csv(T/'GSE174188_frozen_one_batch_per_donor.tsv',sep='\t',index=False)
bar=pd.read_csv(S/'10x/barcodes.tsv.gz',header=None,sep='\t')[0].astype(str).to_numpy();assert len(bar)==len(m)
assert np.array_equal(bar,m.cell_id.astype(str).to_numpy()),'Do not assume matrix/meta row order'
genes=pd.read_csv(S/'10x/features.tsv.gz',header=None,sep='\t')[0].astype(str)
print('Reading sparse source MTX (no full-cell copy is written)',flush=True)
with threadpool_limits(limits=8):
 with gzip.open(S/'10x/matrix.mtx.gz','rb') as f:x=mmread(f)
print('loaded',x.shape,x.nnz,flush=True)
assert x.shape==(len(genes),len(m))
assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.rint(x.data))
x=x.tocsr();print('csr',flush=True)
# Match one source biological individual + one selected sequencing batch; never
# combine batches/visits into independent pseudoreplicates.
key=m.ind_cov.astype(str)+'|'+m.batch_cov.astype(str)
choice_keys=chosen.unit_id.tolist();ids=pd.Index(choice_keys).get_indexer(key)
watch_idx=[i for i,g in enumerate(genes) if g in IFN]
lib=np.asarray(x.sum(axis=0)).ravel();watch=x[watch_idx,:].toarray()
ifn=np.log1p(watch/np.maximum(lib,1)*1e4).mean(axis=0)
scores=m[['ind_cov','batch_cov','cohort','condition','cg_cov']].copy();scores['IFN12_log1p10k_mean']=ifn;scores['chosen_batch']=ids>=0
score=scores[scores.chosen_batch].groupby(['ind_cov','cohort','condition','cg_cov'],observed=True).agg(IFN12_log1p10k_mean=('IFN12_log1p10k_mean','mean'),n_cells=('IFN12_log1p10k_mean','size')).reset_index()
score.to_csv(T/'GSE174188_target_cell_IFN_by_unique_donor.tsv',sep='\t',index=False)
for ct in ['cM','T4','B']:
 for mode in ['primary_unique_batch','cohort4_sensitivity']:
  if mode=='primary_unique_batch':
   pick=(ids>=0)&m.cg_cov.eq(ct).to_numpy();pick_ids=ids[pick];md=chosen.set_index('unit_id').copy()
  else:
   pick=m.cohort.eq(4).to_numpy()&m.cg_cov.eq(ct).to_numpy()
   # Within this single cohort, choose one largest source batch for each donor.
   ch4=eligible[eligible.cohort.eq(4)].sort_values(['ind_cov','n_cells','batch_cov'],ascending=[True,False,True]).drop_duplicates('ind_cov').copy();ch4['unit_id']=ch4.ind_cov+'|'+ch4.batch_cov
   ix4=pd.Index(ch4.unit_id).get_indexer(key);pick=pick&(ix4>=0);pick_ids=ix4[pick];md=ch4.set_index('unit_id').copy()
  selected=np.flatnonzero(pick);one=csr_matrix((np.ones(len(selected),dtype=np.int64),(selected,pick_ids)),shape=(len(m),len(md)))
  pb=(x@one).toarray().T
  nc=np.bincount(pick_ids,minlength=len(md));md['target_cells']=nc
  keep=nc>=20;md=md.iloc[np.flatnonzero(keep)].copy();pb=pb[keep]
  assert md.ind_cov.is_unique
  df=pd.DataFrame(pb,index=md.index,columns=genes)
  df=df.T.groupby(level=0,sort=False).sum().T.astype(np.int64)
  prefix=f'GSE174188_{ct}_{mode}'
  df.to_csv(T/(prefix+'_counts.tsv'),sep='\t');md.to_csv(T/(prefix+'_metadata.tsv'),sep='\t')
  print(prefix,len(md),md.condition.value_counts().to_dict(),int(md.target_cells.sum()),flush=True)
summary={'source_n_cells':len(m),'source_n_individuals':m.ind_cov.nunique(),'eligible_unique_individuals':len(chosen),'selected_n_HD':int(chosen.condition.eq('HD').sum()),'selected_n_SLE':int(chosen.condition.eq('SLE').sum()),'invariant_checks_multiple_values':invariant,'one_batch_selection_rule':'Largest total source cells among source Processing_Cohort 2/3/4 per ind_cov; lexical batch_cov tie-break; source cell counts, not expression or target abundance, determine selection.','source_metadata_SHA256':hashlib.sha256((S/'GSE174188_metadata.csv').read_bytes()).hexdigest(),'source_counts_note':'Locally remapped archival processed source matrix, not primary sequencing reads or unfiltered droplets. Barcode order exactly checked against source metadata.','source_counts_dimensions':list(x.shape),'source_counts_nnz':int(x.nnz),'min_target_cells':20,'IFN12_genes':IFN,'excludes_cohort1':'HD-only; no disease contrast possible'}
(B/'pseudobulk_extraction_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)

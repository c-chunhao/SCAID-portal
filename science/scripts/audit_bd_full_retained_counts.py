import os
"""Read-only, bounded, sparse full retained-BD count audit against four 10x sources.
Replays archived CreateSeurat min.features=200 before min.cells=3 inclusion, harmonizes date-corrupted
feature symbols through original stable Ensembl IDs. Does not repeat QC or doublets.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
B=Path(__file__).resolve().parents[1];T=B/'tables'
SOURCE=Path(os.environ['SCAID_BD_SOURCE_DIR'])
PUB=Path(os.environ['SCAID_BD_RELEASE_H5'])
api=json.loads((T/'BD_date_gene_ensembl_lookup.json').read_text())
def sha(p):
 z=hashlib.sha256()
 with open(p,'rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),''): 
   if not chunk:break
   z.update(chunk)
 return z.hexdigest()
def csr_sha(a,barcodes,genes):
 a=a.copy();a.sum_duplicates();a.eliminate_zeros();a.sort_indices();z=hashlib.sha256()
 z.update(('\n'.join(barcodes)+'\n'+'\n'.join(genes)).encode())
 for v in [np.array(a.shape),a.indptr,a.indices,a.data]:z.update(np.asarray(v,dtype='<i8').tobytes())
 return z.hexdigest()
def hcol(h,k):
 n=h['obs'][k]
 if isinstance(n,h5py.Group):return n['categories'].asstr()[:][n['codes'][:]]
 return n[:]
start=time.time()
with h5py.File(PUB) as h:
 ids=h['names_obs'].asstr()[:];genes=h['names_var'].asstr()[:];sample=hcol(h,'Sample');ncount=hcol(h,'nCount_RNA');nfeature=hcol(h,'nFeature_RNA')
 parts=[]
 for key in sorted(h['assay/RNA/layers/rawdata']):
  bl=h['assay/RNA/layers/rawdata'][key];ptr=bl['indptr'][:];parts.append(sparse.csr_matrix((bl['data'][:],bl['indices'][:],ptr),shape=(len(ptr)-1,len(genes))))
 p=sparse.vstack(parts,format='csr');p.sum_duplicates();p.sort_indices();assert p.shape==(20375,16973)
assert np.unique(ids).size==len(ids) and np.unique(genes).size==len(genes)
assert np.all(np.isfinite(p.data)) and np.array_equal(p.data,np.rint(p.data))
gene_ix={g:i for i,g in enumerate(genes)};results=[];blocks=[];index=[];featuremaps=[]
for label in sorted(set(sample)):
 d=SOURCE/label;f=pd.read_csv(d/'features.tsv.gz',sep='\t',header=None,names=['ensembl_id','source_symbol','feature_type']);seen={};uniq=[]
 for x in f.source_symbol:
  n=seen.get(x,0);seen[x]=n+1;uniq.append(x if not n else x+'.'+str(n))
 f['source_unique_symbol']=uniq;f['canonical_symbol']=[(api.get(e)or{}).get('display_name',g).replace('_','-') for e,g in zip(f.ensembl_id,uniq)]
 assert f.canonical_symbol.is_unique, 'Cannot collapse source genes'
 featuremaps.append(f.assign(sample=label))
 bars=pd.read_csv(d/'barcodes.tsv.gz',sep='\t',header=None)[0].to_numpy();assert np.unique(bars).size==len(bars);bi={x:i for i,x in enumerate(bars)}
 with threadpool_limits(limits=1):a=mmread(d/'matrix.mtx.gz').tocsr()
 a.sum_duplicates();a.eliminate_zeros();assert a.shape==(len(f),len(bars));assert np.all(a.data>=0) and np.array_equal(a.data,np.rint(a.data))
 # The historical script used min.cells=3 independently before specimen merge.
 precell_keep=np.diff(a.tocsc().indptr)>=200;keep=np.asarray((a[:,precell_keep]>0).sum(axis=1)).ravel()>=3;missing=(~f.canonical_symbol.isin(genes)) & keep
 if missing.any():print(f.loc[missing].to_string(index=False),flush=True)
 assert not missing.any(), 'Source min.cells-kept gene lacks a defensible released-name match'
 ix=np.where(sample==label)[0];raw=[x.rsplit('_',1)[0]for x in ids[ix]];assert all(x in bi for x in raw)
 cells=np.array([bi[x]for x in raw]);kept=np.where(keep)[0];sub=a[kept,:][:,cells].T.tocoo();mapped=np.array([gene_ix[g]for g in f.canonical_symbol.iloc[kept]])
 aligned=sparse.csr_matrix((sub.data,(sub.row,mapped[sub.col])),shape=(len(ix),len(genes)));aligned.sum_duplicates();aligned.sort_indices()
 published=p[ix,:];diff=published-aligned;diff.eliminate_zeros()
 if diff.nnz:
  dd=diff.tocoo();diag=[]
  for cr,gc,dv in zip(dd.row,dd.col,dd.data):
   sf=np.where(f.canonical_symbol.to_numpy()==genes[gc])[0];diag.append({'sample':label,'cell':ids[ix[cr]],'gene':genes[gc],'published':int(published[cr,gc]),'source_aligned':int(aligned[cr,gc]),'source_gene_row':int(sf[0]) if len(sf)==1 else None,'source_gene_support':int(np.diff(a.indptr)[sf[0]]) if len(sf)==1 else None})
  print(json.dumps(diag[:30],indent=2),flush=True);pd.DataFrame(diag).to_csv(T/'BD_full_retained_numeric_difference_diagnostic.tsv',sep='\t',index=False)
 assert diff.nnz==0, f'{label}: {diff.nnz} real numeric differences'
 assert np.array_equal(np.asarray(aligned.sum(axis=1)).ravel(),ncount[ix]);assert np.array_equal(np.diff(aligned.indptr),nfeature[ix])
 srch=csr_sha(aligned,ids[ix],genes);pubh=csr_sha(published,ids[ix],genes);assert srch==pubh
 results.append({'sample':label,'source_cells':len(bars),'source_features':len(f),'source_min3_features':int(keep.sum()),'source_min200_cells':int(precell_keep.sum()),'retained_cells_exact_barcode_matches':len(ix),'retained_gene_universe':len(genes),'retained_sparse_nonzero_values':int(aligned.nnz),'remaining_numeric_differences':int(diff.nnz),'canonical_sparse_SHA256':srch,'published_canonical_sparse_SHA256':pubh,'nCount_exact_all_cells':True,'nFeature_exact_all_cells':True,'gene_name_collisions':0,'source_matrix_file_SHA256':sha(d/'matrix.mtx.gz'),'source_features_file_SHA256':sha(d/'features.tsv.gz'),'source_barcodes_file_SHA256':sha(d/'barcodes.tsv.gz')})
 blocks.append(aligned);index.extend(ix.tolist());print(label,len(ix),aligned.nnz,'EXACT',flush=True)
 assert time.time()-start<590, 'Bounded audit exceeded time budget'
r=sparse.vstack(blocks,format='csr')[np.argsort(index),:];assert np.array_equal(np.sort(index),np.arange(len(ids)))
pdigest=csr_sha(p,ids,genes);rdigest=csr_sha(r,ids,genes);assert pdigest==rdigest
pd.concat(featuremaps).to_csv(T/'BD_full_source_feature_stable_ID_map.tsv',sep='\t',index=False)
pd.DataFrame(results).to_csv(T/'BD_full_retained_count_audit_by_specimen.tsv',sep='\t',index=False)
out={'source_object':'SV_GSE198616_PBMC_RNA','scope':'All 20,375 retained case-object barcodes, all 16,973 released RNA genes and all retained sparse count values across four PBMC_BD source specimens. Healthy-source cells are not included in this released case object.','passed':True,'source_original_numeric_counts_reconciled':True,'source_gene_inclusion':'Original archived CreateSeuratObject min.features=200 cell inclusion followed by min.cells=3 per source specimen; stable original Ensembl IDs correct date-corrupted feature symbols; the archived Seurat CheckFeaturesNames underscore-to-hyphen normalization is replayed (two source feature symbols in PBMC_BD_4). Source sparse counts at genes excluded by this archived gene filter are outside the released RNA universe.','retained_cells':len(ids),'retained_genes':len(genes),'retained_sparse_nonzero_values':int(r.nnz),'remaining_numeric_differences':0,'nCount_and_nFeature_exact_all_retained_cells':True,'canonical_sparse_SHA256':rdigest,'published_canonical_sparse_SHA256':pdigest,'published_H5_file_SHA256':sha(PUB),'per_specimen':results,'time_seconds':round(time.time()-start,2),'count_values_changed':0,'source_or_release_objects_mutated':False,'limits':'Establishes count-value provenance for this released BD case object, not full raw/unfiltered sequencing authenticity, independent donor diagnosis beyond source labels, original doublet/ambient execution, fresh uniform QC of 48 objects, or healthy/control objects.'}
(B/'BD_full_retained_source_count_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

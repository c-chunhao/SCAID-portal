import os
"""Read-only full retained-pSS sparse counts/source provenance; no fresh QC.
Uses archived min.features=200 then min.cells=3 on all original source groups.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np,pandas as pd,h5py
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
B=Path(__file__).resolve().parents[1];T=B/'tables';SOURCE=Path(os.environ['SCAID_PSS_SOURCE_DIR']);PUB=Path(os.environ['SCAID_PSS_RELEASE_H5'])
def sha(p):
 z=hashlib.sha256()
 with open(p,'rb')as f:
  for v in iter(lambda:f.read(1024*1024),b''):z.update(v)
 return z.hexdigest()
def csr_sha(a,barcodes,genes):
 a=a.copy();a.sum_duplicates();a.eliminate_zeros();a.sort_indices();z=hashlib.sha256();z.update(('\n'.join(barcodes)+'\n'+'\n'.join(genes)).encode())
 for v in [np.array(a.shape),a.indptr,a.indices,a.data]:z.update(np.asarray(v,dtype='<i8').tobytes())
 return z.hexdigest()
start=time.time()
with h5py.File(PUB)as h:
 ids=h['names_obs'].asstr()[:];genes=h['names_var'].asstr()[:];bls=h['assay/RNA/layers/rawdata'];parts=[]
 for key in sorted(bls):
  bl=bls[key];ptr=bl['indptr'][:];parts.append(sparse.csr_matrix((bl['data'][:],bl['indices'][:],ptr),shape=(len(ptr)-1,len(genes))))
 p=sparse.vstack(parts,format='csr');p.sum_duplicates();p.eliminate_zeros();p.sort_indices();assert p.shape==(29567,21176)
 ncount=h['obs/nCount_RNA'][:];nfeature=h['obs/nFeature_RNA'][:];n=h['obs/Sample'];sample=n['categories'].asstr()[:][n['codes'][:]]
f=pd.read_csv(SOURCE/'features.tsv.gz',sep='\t',header=None);raw=f.iloc[:,1].astype(str).to_numpy();seen={};uniq=[]
for gene in raw:
 n=seen.get(gene,0);seen[gene]=n+1;uniq.append((gene if not n else gene+'.'+str(n)).replace('_','-'))
assert np.unique(uniq).size==len(uniq)
bars=pd.read_csv(SOURCE/'barcodes.tsv.gz',sep='\t',header=None)[0].to_numpy();assert np.unique(bars).size==len(bars);bi={x:i for i,x in enumerate(bars)};assert all(x in bi for x in ids);ix=np.array([bi[x]for x in ids]);gi={x:i for i,x in enumerate(genes)}
with threadpool_limits(limits=1):a=mmread(SOURCE/'matrix.mtx.gz').tocsr()
a.sum_duplicates();a.eliminate_zeros();assert a.shape==(len(f),len(bars));assert (a.data>=0).all()and np.array_equal(a.data,np.rint(a.data));pre=np.diff(a.tocsc().indptr)>=200;keep=np.asarray((a[:,pre]>0).sum(axis=1)).ravel()>=3;assert all(g in gi for g in np.asarray(uniq)[keep]), 'Unresolved gene names';kept=np.where(keep)[0]
sub=a[kept,:][:,ix].T.tocoo();mapped=np.array([gi[g]for g in np.asarray(uniq)[keep]]);r=sparse.csr_matrix((sub.data,(sub.row,mapped[sub.col])),shape=p.shape);r.sum_duplicates();r.eliminate_zeros();r.sort_indices();diff=p-r;diff.eliminate_zeros();assert diff.nnz==0,f'{diff.nnz} numeric differences';assert np.array_equal(np.asarray(r.sum(axis=1)).ravel(),ncount);assert np.array_equal(np.diff(r.indptr),nfeature)
rdigest=csr_sha(r,ids,genes);pdigest=csr_sha(p,ids,genes);assert rdigest==pdigest
rows=[]
for s in sorted(set(sample)):
 j=np.where(sample==s)[0];rows.append({'sample':s,'retained_cells':len(j),'retained_sparse_nonzero_values':r[j,:].nnz,'numeric_differences':0,'canonical_sparse_SHA256':csr_sha(r[j,:],ids[j],genes)})
pd.DataFrame(rows).to_csv(T/'pSS_full_retained_count_audit_by_specimen.tsv',sep='\t',index=False)
out={'source_object':'SS_GSE157278_PBMC','passed':True,'source_original_numeric_counts_reconciled':True,'retained_cells':len(ids),'retained_genes':len(genes),'retained_sparse_nonzero_values':int(r.nnz),'remaining_numeric_differences':0,'nCount_and_nFeature_exact_all_retained_cells':True,'canonical_sparse_SHA256':rdigest,'published_canonical_sparse_SHA256':pdigest,'published_H5_file_SHA256':sha(PUB),'source_cells_all_groups':len(bars),'source_features_all_groups':len(f),'min_features200_import_cells':int(pre.sum()),'min_cells3_genes':int(keep.sum()),'source_original_files_SHA256':{fn:sha(SOURCE/fn) for fn in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz','GSE157278_cell_batch.tsv.gz']},'per_specimen':rows,'gene_filter':'Archived CreateSeuratObject applies min.features=200 before min.cells=3 on all source groups. Released disease cells matched by exact source barcodes; no expression harmonization or fresh doublet calling.','source_or_release_objects_mutated':False,'count_values_changed':0,'time_seconds':round(time.time()-start,2),'limits':'Count provenance for all retained pSS case-object cells, not healthy arms, raw read authenticity, all48 uniform upstreamQC or independent clinical/biological validation. Disease-versus-healthy design remains confounded with source heat-shock library preparation and cannot establish shared cross-disease programmes.'}
(B/'pSS_full_retained_source_count_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

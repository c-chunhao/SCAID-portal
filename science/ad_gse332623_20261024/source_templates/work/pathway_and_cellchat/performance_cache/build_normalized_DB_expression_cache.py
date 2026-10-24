"""Internal source-preserving CellChat DB subset, copied from the verified gene CSC."""
import csv,gzip,hashlib,json,os,resource,time
from datetime import datetime,timezone
from pathlib import Path
import h5py,numpy as np
W=Path('@WORK@')
A=Path('@PROJECT@/00.corrected_analysis/GSE332623')
D=W/'pathway_and_cellchat/performance_cache';OUT=D/'normalized_DB_expression_cache.h5';TMP=OUT.with_suffix('.h5.tmp');RECEIPT=D/'normalized_DB_expression_cache_receipt.json'
RNA_SHA='b9aa963123cb8f010d3e092e03507b70d02170afd93d32a6cdea2efe2a2c09ec';INDEX_SHA='a1bf380298c897810fa2a8fdac835251f437f42a10a6890c16999bae7ffdf0fa'
ALIAS=W/'analysis/validation/hgnc_alias_review/cellchat_symbol_aliases_verified.tsv';VIEW=A/'pathway_and_cellchat/cellchat/CellChat_RNA_symbol_matching_view.tsv';INCLUSION=A/'pathway_and_cellchat/cellchat/cell_inclusion.tsv.gz'
START=time.monotonic()
def progress(**kw):print(json.dumps(dict(elapsed_seconds=round(time.monotonic()-START,2),**kw),ensure_ascii=False),flush=True)
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def axis_sha(x):return hashlib.sha256(json.dumps(list(x),ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
def save_report(d):
 p=RECEIPT.with_suffix('.json.tmp');p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n');os.replace(p,RECEIPT)
def wait_memory():
 while True:
  mem={r.split(':')[0]:int(r.split()[1])*1024 for r in Path('/proc/meminfo').read_text().splitlines() if r.startswith('MemAvailable:')}
  available=mem['MemAvailable']
  if available>=320*1024**3:return
  progress(status='paused',reason='system MemAvailable below320GiB',available_bytes=available);time.sleep(10)
def table(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def strings(x):
 if isinstance(x,h5py.Group):
  cats=x['categories'].asstr()[:];codes=x['codes'][:];assert np.all(codes>=0);return cats[codes]
 return x.asstr()[:]
def putstr(g,n,v):g.create_dataset(n,data=np.asarray(v,dtype=object),dtype=h5py.string_dtype('utf-8'),compression='gzip',compression_opts=1)
def main():
 assert not OUT.exists(), 'Do not overwrite a validated cache';resume_validation=TMP.exists()
 wait_memory();entry=json.loads((W/'website/prepared/gene_entry.json').read_text());index=Path(entry['index_path']);rna=Path(entry['source_h5ad'])
 assert entry['index_sha256']==INDEX_SHA and entry['source_sha256']==RNA_SHA
 assert index.stat().st_size==entry['index_bytes']
 assert sha(index)==INDEX_SHA
 assert sha(ALIAS)=='51873d0026a7cb6657298baf266f7de7c4f9b1b7ab132e4b2193c932703b18cd'
 views=table(VIEW);aliases=table(ALIAS);amap={r['source_symbol']:r['canonical_symbol'] for r in aliases};assert len(amap)==9 and len(views)==1434
 sources=[r['source_symbol'] for r in views];canonical=[r['canonical_symbol'] for r in views]
 assert len(set(sources))==len(sources)==len(set(canonical))
 assert sum(s!=c for s,c in zip(sources,canonical))==9
 assert all(c==s or amap.get(s)==c for s,c in zip(sources,canonical))
 progress(status='source_hashes_validated')
 with h5py.File(index,'r') as src,h5py.File(rna,'r') as raw:
  assert src.attrs['source_sha256']==RNA_SHA and bool(src.attrs['full_source_component_readback'])
  genes=list(src['genes'].asstr()[:]);cells=src['cells'].asstr()[:];n=len(cells);assert n==581725
  mapping={s:i for i,s in enumerate(genes)};cols=np.asarray([mapping[s] for s in sources],dtype=np.int64);assert np.all(np.diff(cols)>0)
  assert axis_sha(cells)==entry['cell_axis_sha256']==src.attrs['cell_axis_sha256']
  assert np.array_equal(cells,strings(raw['obs'][raw['obs'].attrs['_index']]))
  assert genes==list(strings(raw['var'][raw['var'].attrs['_index']]))
  meta={k:strings(raw['obs'][k]) for k in ['MajorCellType','source_library_id','source_pool','source_barcode','source_rna_gsm','source_adt_gsm','Group','source_condition','donor_id']}
  assert np.array_equal(meta['source_library_id'],strings(raw['obs/library_id'])) and np.array_equal(meta['source_pool'],strings(raw['obs/pool_id']))
  assert all(c==l+'_'+b for c,l,b in zip(cells,meta['source_library_id'],meta['source_barcode']))
  assert len(set(meta['source_library_id']))==27 and len(set(meta['source_pool']))==7
  assert all(set(meta[k])=={'Unresolved'} for k in ('Group','source_condition','donor_id'))
  assert not np.asarray(raw['obs/clinical_comparison_eligible'][:],dtype=bool).any()
  include=np.zeros(n,dtype=np.uint8);reasons=[]
  with gzip.open(INCLUSION,'rt') as f:
   for i,r in enumerate(csv.DictReader(f,delimiter='\t')):
    assert i<n and r['cell_id']==cells[i] and r['MajorCellType']==meta['MajorCellType'][i]
    assert r['included'] in ('TRUE','FALSE');include[i]=(r['included']=='TRUE');reasons.append(r['reason'])
  assert i+1==n and int(include.sum())==528739
  ptr=src['expression_csc/indptr'][:];sizes=ptr[cols+1]-ptr[cols];outptr=np.r_[np.int64(0),np.cumsum(sizes,dtype=np.int64)];nnz=int(outptr[-1])
  assert tuple(src['expression_csc'].attrs['shape'])==(n,len(genes))
  if not resume_validation:
   with h5py.File(TMP,'w') as dst:
    dst.attrs.update(schema='scaid-internal-normalized-rna-db-csc-v1',internal_only=True,public_RNA_axis_modified=False,source_RNA_sha256=RNA_SHA,source_gene_index_sha256=INDEX_SHA,cellchat_9_alias_sha256=sha(ALIAS),matching_view_sha256=sha(VIEW),cell_inclusion_sha256=sha(INCLUSION),expression_values='Original RNA X log1p(total10000), copied unchanged; not raw counts and not decontaminated values',cell_count=n,gene_count=len(cols),purpose='Internal CellChat DB resource-view cache; no public RNA gene renaming and no model run')
    g=dst.create_group('expression_csc');g.attrs.update(encoding_type='csc_matrix',shape=np.asarray([n,len(cols)],dtype=np.int64))
    data=g.create_dataset('data',shape=(nnz,),dtype='float64',chunks=(min(nnz,131072),),compression='gzip',compression_opts=1,shuffle=True)
    inds=g.create_dataset('indices',shape=(nnz,),dtype='int32',chunks=(min(nnz,131072),),compression='gzip',compression_opts=1,shuffle=True)
    g.create_dataset('indptr',data=outptr,compression='gzip',compression_opts=1)
    gg=dst.create_group('gene');putstr(gg,'source_symbols',sources);putstr(gg,'canonical_symbols',canonical);gg.create_dataset('source_gene_index',data=cols)
    obs=dst.create_group('obs');putstr(obs,'cell_id',cells)
    for key,val in meta.items():putstr(obs,key,val)
    putstr(obs,'cellchat_inclusion_reason',reasons);obs.create_dataset('include_cellchat',data=include,compression='gzip',compression_opts=1);obs.create_dataset('source_cell_index',data=np.arange(n,dtype=np.int64),compression='gzip',compression_opts=1);obs.create_dataset('clinical_comparison_eligible',data=np.zeros(n,dtype=np.uint8),compression='gzip',compression_opts=1)
    for j,col in enumerate(cols):
     if j%16==0:wait_memory()
     a,b=int(ptr[col]),int(ptr[col+1]);v=src['expression_csc/data'][a:b];ix=src['expression_csc/indices'][a:b]
     assert np.all(np.isfinite(v)) and np.all(v>0) and np.all(np.diff(ix)>0) and (not len(ix) or (ix[0]>=0 and ix[-1]<n))
     first,last=int(outptr[j]),int(outptr[j+1]);data[first:last]=v;inds[first:last]=ix
     if j%100==0:progress(status='copying',genes_complete=j+1,genes_total=len(cols),nnz_copied=last,nnz_total=nnz)
  # Independently read back every extracted component against its complete source column.
  with h5py.File(TMP,'r') as dst:
   assert dst.attrs['source_RNA_sha256']==RNA_SHA and dst.attrs['source_gene_index_sha256']==INDEX_SHA
   assert dst.attrs['matching_view_sha256']==sha(VIEW) and dst.attrs['cellchat_9_alias_sha256']==sha(ALIAS) and dst.attrs['cell_inclusion_sha256']==sha(INCLUSION)
   assert np.array_equal(strings(dst['gene/source_symbols']),sources) and np.array_equal(strings(dst['gene/canonical_symbols']),canonical)
   assert np.array_equal(dst['gene/source_gene_index'][:],cols) and tuple(dst['expression_csc'].attrs['shape'])==(n,len(cols))
   outdatahash=hashlib.sha256();outindexhash=hashlib.sha256()
   assert np.array_equal(dst['expression_csc/indptr'][:],outptr)
   for j,col in enumerate(cols):
    if j%16==0:wait_memory()
    a,b=int(ptr[col]),int(ptr[col+1]);p,q=int(outptr[j]),int(outptr[j+1]);v=dst['expression_csc/data'][p:q];ix=dst['expression_csc/indices'][p:q]
    assert np.array_equal(v,src['expression_csc/data'][a:b]) and np.array_equal(ix,src['expression_csc/indices'][a:b])
    outdatahash.update(v.tobytes());outindexhash.update(ix.tobytes())
   assert np.array_equal(strings(dst['obs/cell_id']),cells) and np.array_equal(dst['obs/include_cellchat'][:],include)
   for k in meta:assert np.array_equal(strings(dst['obs/'+k]),meta[k])
   # Direct original RNA CSR sample, covering every library and all9 renamed genes.
   sample_cells=[]
   for lib in dict.fromkeys(meta['source_library_id']):sample_cells.append(int(np.flatnonzero(meta['source_library_id']==lib)[0]))
   sample_cells=sorted(set(sample_cells+list(np.linspace(0,n-1,48,dtype=int))))
   sample_cols=sorted(set(list(np.linspace(0,len(cols)-1,48,dtype=int))+[j for j,(s,c) in enumerate(zip(sources,canonical)) if s!=c]))
   original_ptr=raw['X/indptr'][:];wanted=cols[sample_cols];expected=np.zeros((len(sample_cells),len(sample_cols)),dtype=np.float64);actual=expected.copy()
   for ir,cell in enumerate(sample_cells):
    a,b=map(int,(original_ptr[cell],original_ptr[cell+1]));ix=raw['X/indices'][a:b];v=raw['X/data'][a:b];pos=np.searchsorted(ix,wanted);valid=pos<len(ix);safe=np.minimum(pos,max(0,len(ix)-1));valid&=(ix[safe]==wanted)
    expected[ir,valid]=v[pos[valid]]
   wanted_cells=np.asarray(sample_cells)
   for jc,col in enumerate(sample_cols):
    a,b=map(int,(outptr[col],outptr[col+1]));ix=dst['expression_csc/indices'][a:b];v=dst['expression_csc/data'][a:b]
    if len(ix)==0:continue
    pos=np.searchsorted(ix,wanted_cells);valid=pos<len(ix);safe=np.minimum(pos,max(0,len(ix)-1));valid&=(ix[safe]==wanted_cells)
    actual[valid,jc]=v[pos[valid]]
   assert np.array_equal(expected,actual)
  os.replace(TMP,OUT);OUT.chmod(0o444)
  report={'status':'validated_internal_cache_not_a_model_run','created_utc':datetime.now(timezone.utc).isoformat(),'schema':'scaid-internal-normalized-rna-db-csc-v1','path':str(OUT),'bytes':OUT.stat().st_size,'sha256':sha(OUT),'shape':[n,len(cols)],'nnz':nnz,'zero_nnz_genes_retained':int((sizes==0).sum()),'resumed_existing_complete_component_cache':resume_validation,'modeled_cells_mask_sum':int(include.sum()),'source_RNA_path':str(rna),'source_RNA_sha256':RNA_SHA,'source_gene_index_path':str(index),'source_gene_index_sha256':INDEX_SHA,'source_gene_index_actual_full_file_SHA_verified':True,'matching_view_path':str(VIEW),'matching_view_sha256':sha(VIEW),'cellchat_9_alias_path':str(ALIAS),'cellchat_9_alias_sha256':sha(ALIAS),'cell_inclusion_path':str(INCLUSION),'cell_inclusion_sha256':sha(INCLUSION),'cell_axis_JSON_sha256':axis_sha(cells),'source_gene_axis_JSON_sha256':axis_sha(sources),'canonical_gene_axis_JSON_sha256':axis_sha(canonical),'data_component_sha256':outdatahash.hexdigest(),'indices_component_sha256':outindexhash.hexdigest(),'indptr_component_sha256':hashlib.sha256(outptr.tobytes()).hexdigest(),'source_gene_order_strictly_increasing':True,'every_extracted_column_sorted_unique_finite_and_exact_source_component_readback':True,'entire_cell_library_pool_barcode_and_inclusion_axes_verified':True,'original_RNA_direct_sample':{'cells':len(sample_cells),'genes':len(sample_cols),'pairs':int(expected.size),'all27_libraries':True,'all9_restored_aliases':True,'exact_equality':True},'resources':{'memory_hard_limit_GiB':8,'CPU_quota':1,'pause_system_available_floor_GiB':320,'peak_process_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},'elapsed_seconds':time.monotonic()-START,'source_RNA_or_public_indices_modified':False,'model_started':False,'interpretation':'Internal normalized DB matching view only. Original public RNA feature names/axis are preserved. Candidate optimized CellChat may consume only after full100bootstrap equivalence approval.'}
  save_report(report);progress(status=report['status'],shape=report['shape'],nnz=nnz,receipt=str(RECEIPT))
if __name__=='__main__':
 try:main()
 except BaseException as e:
  save_report({'status':'failed','created_utc':datetime.now(timezone.utc).isoformat(),'error_type':type(e).__name__,'error':str(e),'model_started':False})
  raise

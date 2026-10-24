"""Bounded-memory, source-preserving CSR-to-CSC gene index conversion."""
from pathlib import Path
from urllib.parse import quote
import collections,hashlib,json,mmap,os,shutil,tempfile,time
import h5py
import numpy as np

def digest(path):
 with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def axis_sha(values):return hashlib.sha256(json.dumps(list(values),ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
def gene_filename(gene):
    filename='Umap_'+quote(gene,safe='-_.~')+'.png'
    return filename if len(filename.encode('utf-8'))<=255 else 'Umap_feature_'+hashlib.sha256(gene.encode('utf-8')).hexdigest()+'.png'
def safe_name(g):return bool(g) and g==g.strip() and g.isprintable() and len('Umap_'+g+'.png')<=255 and len(gene_filename(g).encode())<=255


def array_chunks(dataset,limit):
 for start in range(0,len(dataset),limit):yield start,np.asarray(dataset[start:start+limit])


def release_pages(array):
 array.flush()
 if hasattr(array._mmap,'madvise'):array._mmap.madvise(mmap.MADV_DONTNEED)


def row_batches(indptr,max_nnz):
 n=len(indptr)-1;start=0
 while start<n:
  end=int(np.searchsorted(indptr,indptr[start]+max_nnz,side='right'))-1
  end=max(start+1,min(end,n));yield start,end;start=end


def build_index(source,output,object_id,source_sha256,*,max_nnz_chunk=2000000):
 started=time.monotonic();source=Path(source);output=Path(output)
 if output.exists():raise ValueError('Index output exists; refuse overwrite')
 if digest(source)!=source_sha256:raise ValueError('RNA H5AD differs from validated checksum')
 output.parent.mkdir(parents=True,exist_ok=True)
 with h5py.File(source,'r') as f:
  group=f['X']
  if group.attrs.get('encoding-type')!='csr_matrix':raise ValueError('RNA X must be CSR cells-by-genes')
  ncell,ngene=map(int,group.attrs['shape']);ptr=np.asarray(group['indptr'],dtype=np.int64);nnz=len(group['data'])
  if len(ptr)!=ncell+1 or ptr[0]!=0 or ptr[-1]!=nnz or np.any(np.diff(ptr)<0) or len(group['indices'])!=nnz:raise ValueError('Invalid source CSR structure')
  cells=list(f['obs'][f['obs'].attrs['_index']].asstr()[:]);genes=list(f['var'][f['var'].attrs['_index']].asstr()[:]);umap=np.asarray(f['obsm/X_umap'],dtype=np.float64)
  if len(cells)!=ncell or len(set(cells))!=ncell or len(genes)!=ngene or len(set(genes))!=ngene or not all(safe_name(g) for g in genes):raise ValueError('Invalid, ambiguous or unsafe source axes')
  if umap.shape!=(ncell,2) or not np.isfinite(umap).all():raise ValueError('Original UMAP is not complete and finite')
  obj=f['uns/scaid/object_id'].asstr()[:]
  if len(obj)!=1 or obj[0]!=object_id:raise ValueError('RNA object provenance does not match the intended cache')
  if nnz==0:raise ValueError('RNA object has no expression entries for a gene index')
  nnz_by_gene=np.zeros(ngene,dtype=np.int64)
  for _,columns in array_chunks(group['indices'],max_nnz_chunk):
   if len(columns) and (columns.min()<0 or columns.max()>=ngene):raise ValueError('Source feature index outside axis')
   nnz_by_gene+=np.bincount(columns,minlength=ngene)
  csc_ptr=np.empty(ngene+1,dtype=np.int64);csc_ptr[0]=0;np.cumsum(nnz_by_gene,out=csc_ptr[1:]);cursor=csc_ptr[:-1].copy()
  need=nnz*12+max_nnz_chunk*128
  if shutil.disk_usage(output.parent).free<need*1.5:raise ValueError('Insufficient temporary disk for bounded-memory conversion')
  # Prefer private tmpfs scratch when its total footprint is safely below
  # currently available RAM. Final output remains on the persistent data disk.
  scratch_parent=output.parent
  shm=Path('/dev/shm')
  try:
   available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
   if nnz*12<=16*1024**3 and shm.is_dir() and shutil.disk_usage(shm).free>nnz*12*2 and available>max(nnz*12*4+2*1024**3,320*1024**3+nnz*12+512*1024**2):
    scratch_parent=shm
  except (OSError,ValueError,StopIteration):pass
  with tempfile.TemporaryDirectory(prefix='scaid-gene-csc-',dir=scratch_parent) as temporary:
   temporary=Path(temporary);data=np.memmap(temporary/'data.float64',mode='w+',dtype=np.float64,shape=(nnz,));rows=np.memmap(temporary/'rows.int32',mode='w+',dtype=np.int32,shape=(nnz,))
   written=0;batches=0;source_value_hash=hashlib.sha256();source_column_hash=hashlib.sha256()
   for first,last in row_batches(ptr,max_nnz_chunk):
    if scratch_parent==Path('/dev/shm'):
     remaining=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
     if remaining<320*1024**3+max_nnz_chunk*12:raise RuntimeError('Memory reserve below 320 GiB: stop and clean tmpfs scratch before further expansion')
    a,b=map(int,(ptr[first],ptr[last]));values=np.asarray(group['data'][a:b],dtype=np.float64);columns=np.asarray(group['indices'][a:b],dtype=np.int64)
    if not np.isfinite(values).all() or (len(values) and values.min()<0):raise ValueError('RNA data contains invalid normalized expression')
    # Existing source CSR rows must contain unique sorted gene indices.
    if len(columns)>1:
     increasing=np.diff(columns)>0;boundaries=ptr[first+1:last]-a;boundaries=boundaries[(boundaries>0)&(boundaries<len(columns))];increasing[boundaries-1]=True
     if not increasing.all():raise ValueError('Source CSR row has duplicate or unsorted gene indices')
    source_value_hash.update(values.tobytes());source_column_hash.update(columns.tobytes())
    cell_rows=np.repeat(np.arange(first,last,dtype=np.int32),np.diff(ptr[first:last+1]))
    order=np.argsort(columns,kind='stable');counts=np.bincount(columns,minlength=ngene);starts=np.cumsum(counts)-counts
    destinations=np.repeat(cursor,counts)+np.arange(len(order),dtype=np.int64)-np.repeat(starts,counts)
    expected_values=values[order];expected_rows=cell_rows[order]
    data[destinations]=expected_values;rows[destinations]=expected_rows
    if not np.array_equal(data[destinations],expected_values) or not np.array_equal(rows[destinations],expected_rows):raise ValueError('Full source-component write/readback failed')
    cursor+=counts;written+=len(values);batches+=1
    release_pages(data);release_pages(rows)
    if batches%10==0:print(json.dumps({'object_id':object_id,'phase':'csc_index_fill','source_entries_written':written,'total_nnz':nnz,'elapsed_seconds':round(time.monotonic()-started,1)}),flush=True)
   if written!=nnz or not np.array_equal(cursor,csc_ptr[1:]):raise ValueError('CSC did not receive every source entry exactly once')
   temp=output.with_suffix('.h5.tmp')
   try:
    with h5py.File(temp,'w') as out:
     out.attrs.update(schema_version=1,object_id=object_id,source_sha256=source_sha256,source_h5ad=str(source),source_layer='RNA data',cell_count=ncell,feature_count=ngene,full_source_component_readback=True)
     g=out.create_group('expression_csc');g.attrs.update(shape=np.array([ncell,ngene],dtype=np.int64),encoding_type='csc_matrix')
     # Empty CSC structures are allowed but real release matrices have nnz>0.
     for name,source_array,dtype in [('data',data,np.float64),('indices',rows,np.int32)]:
      ds=g.create_dataset(name,shape=(nnz,),dtype=dtype,chunks=(min(131072,max(nnz,1)),),compression='gzip',compression_opts=1,shuffle=True)
      for first in range(0,nnz,131072):
       ds[first:first+131072]=source_array[first:first+131072]
       if first%(131072*16)==0 and hasattr(source_array._mmap,'madvise'):source_array._mmap.madvise(mmap.MADV_DONTNEED)
      if hasattr(source_array._mmap,'madvise'):source_array._mmap.madvise(mmap.MADV_DONTNEED)
     g.create_dataset('indptr',data=csc_ptr,compression='gzip',compression_opts=1,shuffle=True)
     out.create_dataset('genes',data=np.asarray(genes,dtype=object),dtype=h5py.string_dtype('utf-8'))
     out.create_dataset('cells',data=np.asarray(cells,dtype=object),dtype=h5py.string_dtype('utf-8'))
     out.create_dataset('umap',data=umap,compression='gzip',compression_opts=1,shuffle=True)
     out.attrs.update(gene_axis_sha256=axis_sha(genes),cell_axis_sha256=axis_sha(cells),umap_sha256=hashlib.sha256(umap.tobytes()).hexdigest())
    # Read every final compressed component against the fully checked memmaps.
    with h5py.File(temp,'r') as out:
     for first in range(0,nnz,131072):
      if not np.array_equal(out['expression_csc/data'][first:first+131072],data[first:first+131072]) or not np.array_equal(out['expression_csc/indices'][first:first+131072],rows[first:first+131072]):raise ValueError('Final compressed CSC differs from source-checked components')
      if first%(131072*16)==0:
       for array in (data,rows):
        if hasattr(array._mmap,'madvise'):array._mmap.madvise(mmap.MADV_DONTNEED)
     if not np.array_equal(out['expression_csc/indptr'][:],csc_ptr) or list(out['genes'].asstr()[:])!=genes or list(out['cells'].asstr()[:])!=cells or not np.array_equal(out['umap'][:],umap):raise ValueError('Final CSC axes or UMAP differ from source')
    os.replace(temp,output);output.chmod(0o444)
   finally:
    if temp.exists():temp.unlink()
   del data,rows
 case_groups=collections.defaultdict(list)
 for gene in genes:case_groups[gene.casefold()].append(gene)
 case_collisions=[v for v in case_groups.values() if len(v)>1]
 stat=output.stat()
 return {'object_id':object_id,'relative_path':object_id+'/GeneUmap','index_path':str(output),'index_bytes':stat.st_size,'index_sha256':digest(output),'index_mtime_ns':stat.st_mtime_ns,
  'source_h5ad':str(source),'source_sha256':source_sha256,'cell_count':ncell,'feature_count':ngene,'nnz':nnz,'genes':genes,'filenames':[gene_filename(g) for g in genes],'gene_axis_sha256':axis_sha(genes),'cell_axis_sha256':axis_sha(cells),
  'umap_sha256':hashlib.sha256(umap.tobytes()).hexdigest(),'source_csr_value_sha256':source_value_hash.hexdigest(),'source_csr_column_sha256':source_column_hash.hexdigest(),
  'source_case_collisions':case_collisions,'filename_encoded_features':[{'gene':g,'filename':gene_filename(g)} for g in genes if gene_filename(g)!='Umap_'+g+'.png'],'all_cells_and_original_coordinates_retained':True,'full_source_component_readback':True,'final_compressed_component_readback':True,'batches':batches,'max_nnz_chunk':max_nnz_chunk,
  'scratch_storage':'tmpfs' if scratch_parent==Path('/dev/shm') else 'persistent_disk','scratch_peak_bytes':nnz*12,'elapsed_seconds':round(time.monotonic()-started,3),'figure_width':1650,'figure_height':1500,'expression_description':'Normalized RNA data-layer values preserved without recomputation'}


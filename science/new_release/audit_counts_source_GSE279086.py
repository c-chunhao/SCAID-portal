"""Independent raw/processed axis audit. Reads no Seurat RDS or matrix body."""
import collections,csv,gzip,hashlib,io,json,os,re,subprocess
from pathlib import Path
import h5py
import numpy as np
from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
OUT=WORK/'matrices/source_recovery_GSE279086'
OUT.mkdir(parents=True,exist_ok=True)
SRC=RAW_ROOT
META=WORK/'matrices/objects/T1D_GSE279086_Kidney/cell_metadata.tsv.gz'
H5=EXPORT_ROOT/'RNA/T1D_GSE279086_Kidney_RNA.h5ad.tmp'
with h5py.File(H5,'r') as h:
 genes=list(h['var/_index'].asstr()[:]); h5cells=list(h['obs/_index'].asstr()[:])
with gzip.open(META,'rt',newline='') as h: meta=list(csv.DictReader(h,delimiter='\t'))
assert h5cells==[r['cell_id'] for r in meta]
assert len(genes)==len(set(genes)) and len(h5cells)==len(set(h5cells))
by_gsm=collections.defaultdict(list)
for i,row in enumerate(meta): by_gsm[row['GSM']].append((i,row))

def lines(path):
 with gzip.open(path,'rt') as h: return [r.rstrip('\r\n') for r in h]
def sha(x): return hashlib.sha256(x).hexdigest()
def header(path):
 with gzip.open(path,'rt') as h:
  magic=h.readline().strip(); comments=[]
  for line in h:
   if line.startswith('%'):comments.append(line.strip());continue
   dimensions=[int(x) for x in line.split()];break
 return {'header':magic,'comments':comments,'dimensions':dimensions,'path':str(path),'compressed_bytes':path.stat().st_size}
def save_tsv(path,fields,rows):
 with gzip.open(path,'wt',newline='') as h:
  w=csv.writer(h,delimiter='\t',lineterminator='\n');w.writerow(fields);w.writerows(rows)
 return sha(path.read_bytes())
report={'status':'in_progress','retained_cells':len(meta),'retained_genes':len(genes),'retained_GSM_count':len(by_gsm),'matrix_bodies_read':False,'terminal_RDS_read':False,'source_root':str(SRC),'rna_axis_path':str(H5),'metadata_path':str(META),'samples':[]}
unique_cache={}
for gsm,group in sorted(by_gsm.items()):
 p=SRC/gsm
 rawf=next(p.glob('*_features.tsv.gz'));prof=next(p.glob('*_features_processed.tsv.gz'))
 rawb=next(p.glob('*_barcodes.tsv.gz'));prob=next(p.glob('*_barcodes_processed.tsv.gz'))
 rawm=next(p.glob('*_matrix.mtx.gz'));prom=next(p.glob('*_matrix_processed.mtx.gz'))
 rawft=[x.split('\t') for x in lines(rawf)]
 assert all(len(x)==3 for x in rawft)
 symbols=[x[1] for x in rawft]
 feature_sha=sha(('\n'.join('\t'.join(x) for x in rawft)+'\n').encode())
 if feature_sha not in unique_cache:
  # R's original make.unique, rather than a guessed Python duplicate algorithm.
  result=subprocess.run([os.environ.get('SCAID_RSCRIPT','Rscript'),'--vanilla','-e','x <- readLines(file("stdin")); cat(make.unique(x), sep="\\n")'],input='\n'.join(symbols)+'\n',text=True,check=True,capture_output=True)
  unique=result.stdout.splitlines()
  assert len(unique)==len(symbols) and len(set(unique))==len(unique)
  unique_cache[feature_sha]=unique
 else: unique=unique_cache[feature_sha]
 raw_to_idx={g:i for i,g in enumerate(unique)}
 processed=lines(prof)
 assert len(processed)==len(set(processed))
 missing_processed=sorted(set(processed)-set(unique))
 missing_genes=sorted(set(genes)-set(unique))
 assert not missing_processed and not missing_genes,(gsm,missing_processed,missing_genes)
 rawbc=lines(rawb);probc=lines(prob)
 assert len(rawbc)==len(set(rawbc)) and len(probc)==len(set(probc))
 raw_bc_idx={g:i for i,g in enumerate(rawbc)}
 pro_bc_idx={g:i for i,g in enumerate(probc)}
 pro_suffix=[x.rsplit('_',1)[-1] for x in probc]
 assert all(re.fullmatch('[ACGT]+-[0-9]+',x) for x in pro_suffix)
 assert len(pro_suffix)==len(set(pro_suffix))
 processed_missing_raw=sorted(set(pro_suffix)-set(rawbc))
 missing_retained=[]
 cell_indices=[]; retained_raw_indices=[]; processed_indices=[];cellrows=[]
 for i,row in group:
  cell=row['cell_id'];assert cell.startswith(gsm+'_')
  procell=cell[len(gsm)+1:]
  rawcell=cell.rsplit('_',1)[-1]
  assert procell in pro_bc_idx,(gsm,cell)
  rawindex=raw_bc_idx.get(rawcell)
  if rawindex is None: missing_retained.append(cell)
  cell_indices.append(i);retained_raw_indices.append(rawindex);processed_indices.append(pro_bc_idx[procell])
  cellrows.append((i,cell,gsm,rawcell,rawindex if rawindex is not None else '',procell,pro_bc_idx[procell]))
 mapped=[i for i in retained_raw_indices if i is not None]
 assert len(set(mapped))==len(mapped)
 rhead,phead=header(rawm),header(prom)
 assert rhead['header']=='%%MatrixMarket matrix coordinate integer general'
 assert phead['header']=='%%MatrixMarket matrix coordinate real general'
 assert rhead['dimensions'][:2]==[len(rawft),len(rawbc)]
 assert phead['dimensions'][:2]==[len(processed),len(probc)]
 gene_indices=[raw_to_idx[g] for g in genes]
 assert len(set(gene_indices))==len(genes)
 gene_rows=[(i,g,idx,rawft[idx][0],rawft[idx][1],unique[idx]) for i,(g,idx) in enumerate(zip(genes,gene_indices))]
 gene_file=OUT/f'{gsm}_gene_mapping.tsv.gz'
 cell_file=OUT/f'{gsm}_cell_mapping.tsv.gz'
 gh=save_tsv(gene_file,['target_feature_index_0based','target_feature','raw_feature_index_0based','raw_ensembl_id','raw_symbol','raw_symbol_R_make_unique'],gene_rows)
 ch=save_tsv(cell_file,['target_cell_index_0based','target_cell_id','GSM','raw_barcode','raw_barcode_index_0based','processed_barcode','processed_barcode_index_0based'],cellrows)
 npfile=OUT/f'{gsm}_axis_mapping.npz'
 if not missing_retained:
  np.savez_compressed(npfile,target_feature_indices=np.arange(len(genes)),raw_feature_indices=np.array(gene_indices),target_cell_indices=np.array(cell_indices),raw_barcode_indices=np.array(retained_raw_indices),processed_barcode_indices=np.array(processed_indices))
 report['samples'].append({'GSM':gsm,'retained_cells':len(group),'retained_genes':len(genes),'raw_feature_rows':len(rawft),'raw_distinct_symbols':len(set(symbols)),'raw_features_sha256_uncompressed':feature_sha,'processed_feature_rows':len(processed),'raw_barcode_rows':len(rawbc),'processed_barcode_rows':len(probc),'processed_barcodes_raw_coverage':not processed_missing_raw,'processed_barcodes_missing_raw':len(processed_missing_raw),'retained_gene_raw_coverage':True,'retained_cell_raw_coverage':not missing_retained,'retained_cells_missing_raw':len(missing_retained),'missing_retained_cells':missing_retained,'safe_for_full_count_recovery':not missing_retained,'retained_cell_processed_coverage':True,'raw_matrix':rhead,'processed_matrix':phead,'gene_mapping_path':str(gene_file),'gene_mapping_sha256':gh,'cell_mapping_path':str(cell_file),'cell_mapping_sha256':ch,'npz_mapping_path':str(npfile) if not missing_retained else None,'npz_mapping_sha256':sha(npfile.read_bytes()) if not missing_retained else None})
 print(gsm,len(group),'cells and',len(genes),'features; missing raw cell axes:',len(missing_retained),flush=True)
report.update(status='passed' if all(r['safe_for_full_count_recovery'] for r in report['samples']) else 'partial_blocked_raw_barcode_mismatch',all_26_samples_unique_raw_axis_mapping=all(r['safe_for_full_count_recovery'] for r in report['samples']),fully_mappable_GSM_count=sum(r['safe_for_full_count_recovery'] for r in report['samples']),missing_retained_cell_count=sum(r['retained_cells_missing_raw'] for r in report['samples']),r_make_unique_used=True,unique_raw_feature_templates=len(unique_cache),raw_matrix_values_verified=False,note='Integer matrix headers and exact axes verified independently. Matrix bodies, integer values and final recovered count data require the recovery/export task verification; no rounding or source modification was performed.')
report['axis_gene_sha256']=sha(('\n'.join(genes)+'\n').encode())
report['axis_cell_sha256']=sha(('\n'.join(h5cells)+'\n').encode())
(OUT/'axis_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
with (OUT/'axis_audit.tsv').open('w',newline='') as h:
 writer=csv.writer(h,delimiter='\t');writer.writerow(['GSM','retained_cells','retained_genes','raw_features','raw_barcodes','processed_barcodes','all_axes_unique_and_covered'])
 for r in report['samples']:writer.writerow([r['GSM'],r['retained_cells'],r['retained_genes'],r['raw_feature_rows'],r['raw_barcode_rows'],r['processed_barcode_rows'],r['safe_for_full_count_recovery']])
print(json.dumps({k:v for k,v in report.items() if k not in {'samples'}},ensure_ascii=False),flush=True)

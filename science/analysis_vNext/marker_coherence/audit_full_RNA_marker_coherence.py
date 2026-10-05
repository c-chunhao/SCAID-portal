"""Read-only source/specimen evidence for original integrated fine labels.

Only a predeclared 40-gene marker matrix is made dense, in <=2048-cell chunks.
Full sparse libraries provide normalization. Markers are coherence evidence,
not annotation accuracy, donor disease tests or permission to relabel cells.
"""
from pathlib import Path
import sys,json,time,hashlib,argparse
import numpy as np,pandas as pd,h5py
from marker_input_paths import resolve_h5_path,prepare_specimen_source,load_object_crosswalk,exact_left_join
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--input-tables',type=Path,required=True,help='Released and original integrated metadata tables directory.')
ap.add_argument('--registry',type=Path,required=True,help='dataset_assets.json with effective raw RNA H5 paths.')
ap.add_argument('--data-h5-root',type=Path,help='Explicit filesystem root for relative registry h5_paths; absolute paths remain explicit.')
ap.add_argument('--reader',type=Path,required=True,help='Canonical integration/released_rna.py.')
ap.add_argument('--crosswalk',type=Path,required=True,help='Explicit legacy object to effective H5 object identity crosswalk.')
ap.add_argument('--output',type=Path,required=True,help='New output directory, separate from immutable inputs.')
args=ap.parse_args();T=args.input_tables.resolve();OUT=args.output.resolve()
if OUT==T or T in OUT.parents:raise ValueError('Output must be separate from immutable input-tables')
OUT.mkdir(parents=True,exist_ok=True)
reader=args.reader.resolve()
if reader.name!='released_rna.py' or not reader.exists():raise ValueError('Supply canonical released_rna.py')
sys.path.insert(0,str(reader.parent))
from released_rna import inspect_schema,validate_raw_csr_arrays
MARKERS=['TRAC','CD3D','CD4','CD8A','CD8B','CCR7','IL7R','FOXP3','IL2RA','RORC','TBX21','CXCR5','PDCD1','LAG3','MKI67','ISG15','SOX4','GZMK','NR4A2','NKG7','GNLY','FCGR3A','KLRD1','KLRC2','KLRB1','ZBTB16','TRDC','TRDV2','MS4A1','CD79A','TCL1A','BCL6','AICDA','JCHAIN','MZB1','SDC1','IGHA1','IGHG1','TRAV1-2','GZMH']
assert len(MARKERS)==40 and len(set(MARKERS))==40
config={'markers':MARKERS,'declared_before_raw_count_scan':True,'canonical_reader':str(reader),'chunk_cells':2048,'mean_expression':'mean log1p(raw gene count / complete raw RNA library size *10000)','detection':'fraction of included cells with original gene count >0','minimum_review_source_or_specimen_cells':20,'missing_feature_policy':'unavailable gene metrics are NA, never zero','scope':'Original integrated fine labels. Evidence is conditional on source, tissue, retained cells and label; no marker-as-accuracy claim, no pooled disease p-values or automatic relabeling.'}
(OUT/'predeclared_marker_panel.json').write_text(json.dumps(config,indent=2))
crosswalk=load_object_crosswalk(args.crosswalk)
m=pd.read_csv(T/'all_release_original_metadata.tsv.gz',sep='\t',usecols=['OldCells','Sample','Dataset','Data_object'],keep_default_na=False)
parts=[]
for lineage in ['T_NK','B_Plasma']:
 d=pd.read_csv(T/f'{lineage}_original_integration_metadata.tsv.gz',sep='\t',usecols=['OldCells','Sample','Dataset','celltype','celltype1','celltype2'],keep_default_na=False)
 d['Lineage']=lineage
 d=exact_left_join(d,m,on=['OldCells','Sample','Dataset'],validate='one_to_one',context=f'{lineage} integrated-to-release');d.Data_object=d.Data_object.replace(crosswalk);parts.append(d)
targets=pd.concat(parts,ignore_index=True);del m,parts,d
assert targets[['Data_object','OldCells']].duplicated().sum()==0
source=prepare_specimen_source(pd.read_csv(T/'composition_specimen_denominators_and_eligibility.tsv',sep='\t',keep_default_na=False),crosswalk)
assert not source[['Data_object','Sample']].duplicated().any()
registry=args.registry;items=json.loads(registry.read_text())['datasets'];coverage=[];qa=[];all_source=[];all_specimen=[]
assert set(targets.Data_object).issubset({resolve_h5_path(item['h5_paths']['RNA'],args.data_h5_root).stem for item in items})
# Prioritize rare/specialist labels so a bounded run still supplies their
# provenance even if later objects exceed a resource/time budget.
priority=targets[targets.celltype2.str.contains('NKT|ILC|MAIT|Bgc|Plasmablast',regex=True)].groupby('Data_object').size().to_dict()
items=sorted(items,key=lambda item:-priority.get(Path(item['h5_paths']['RNA']).stem,0));start_time=time.time()
for num,item in enumerate(items,1):
 path=resolve_h5_path(item['h5_paths']['RNA'],args.data_h5_root);obj=path.stem;sel=targets[targets.Data_object.eq(obj)].copy();print(f'{num}/{len(items)} {obj}: {len(sel)} target cells',flush=True)
 schema=inspect_schema(path)
 if sel.empty:
  with h5py.File(path,'r')as h:genes=np.asarray(h['names_var'].asstr()[:])
  ix={gene:i for i,gene in enumerate(genes)}
  for gene in MARKERS:coverage.append({'Data_object':obj,'Marker':gene,'available_full_RNA':gene in ix,'full_RNA_features':len(genes),'raw_gene_axis_index':ix.get(gene,np.nan)})
  qa.append({'Data_object':obj,'path':str(path),'raw_cells':schema['cells'],'full_RNA_genes':len(genes),'target_cells':0,'audited_target_cells':0,'raw_rows_loaded':0,'fine_labels':0,'specimens':0,'available_markers':sum(g in ix for g in MARKERS),'missing_markers':';'.join(g for g in MARKERS if g not in ix),'status':'schema/gene axes inspected; no integrated members'});continue
 with h5py.File(path,'r')as h:
  ids=np.asarray(h['names_obs'].asstr()[:]);genes=np.asarray(h['names_var'].asstr()[:]);ix={g:i for i,g in enumerate(genes)}
  present=[g for g in MARKERS if g in ix];colidx=[ix[g]for g in present];outcols=[MARKERS.index(g)for g in present]
  for gene in MARKERS:coverage.append({'Data_object':obj,'Marker':gene,'available_full_RNA':gene in ix,'full_RNA_features':len(genes),'raw_gene_axis_index':ix.get(gene,np.nan)})
  lookup=sel.set_index('OldCells');aligned=lookup.reindex(ids);keep=aligned.Lineage.notna().to_numpy();assert int(keep.sum())==len(sel)
  grouping=aligned.loc[keep,['Lineage','celltype','celltype1','celltype2','Sample']].reset_index(drop=True)
  codes,uniques=pd.factorize(pd.MultiIndex.from_frame(grouping));G=len(uniques);groupnames=pd.DataFrame(uniques.tolist(),columns=grouping.columns)
  rowcode=np.full(len(ids),-1,dtype=np.int32);rowcode[keep]=codes
  N=np.bincount(codes,minlength=G);rawsum=np.zeros((G,40),np.float64);logsum=rawsum.copy();det=rawsum.copy();seen=0;offset=0;raw_rows_loaded=0
  for block in sorted(h['assay/RNA/layers/rawdata']):
   g=h['assay/RNA/layers/rawdata'][block];nrows,nfeatures=map(int,g.attrs['shape']);pointer=np.asarray(g['indptr']);assert nfeatures==len(genes)
   for rs in range(0,nrows,2048):
    re=min(rs+2048,nrows);code=rowcode[offset+rs:offset+re];take=code>=0
    if not take.any():continue
    ptr=pointer[rs:re+1];a,b=int(ptr[0]),int(ptr[-1]);data=np.asarray(g['data'][a:b]);index=np.asarray(g['indices'][a:b])
    x=validate_raw_csr_arrays(data,index,ptr-a,(re-rs,nfeatures));lib=np.asarray(x.sum(axis=1,dtype=np.float64)).ravel()[take]
    if np.any(lib<=0):raise ValueError(f'{obj}: target cell has zero complete RNA library')
    # Only requested marker coordinates are dense; all genes contribute to lib.
    small=x[take][:,colidx].toarray().astype(np.float64);norm=np.log1p(small/lib[:,None]*10000)
    for j,col in enumerate(outcols):
     np.add.at(rawsum[:,col],code[take],small[:,j]);np.add.at(logsum[:,col],code[take],norm[:,j]);np.add.at(det[:,col],code[take],small[:,j]>0)
    seen+=int(take.sum());raw_rows_loaded+=re-rs
   offset+=nrows
  assert seen==len(sel) and offset==len(ids)
  for col,gene in enumerate(MARKERS):
   if gene not in ix:rawsum[:,col]=np.nan;logsum[:,col]=np.nan;det[:,col]=np.nan
  groupnames['Cells']=N;groupnames['Data_object']=obj;groupnames=exact_left_join(groupnames,source,on=['Data_object','Sample'],validate='many_to_one',context=f'{obj} specimen source');assert groupnames.Condition.notna().all()
  # Long per-specimen evidence; source donor mapping remains explicitly separate.
  long=groupnames.loc[groupnames.index.repeat(40)].reset_index(drop=True);long['Marker']=MARKERS*G;long['Mean_raw_count']=(rawsum/N[:,None]).ravel();long['Mean_log1p_10k']=(logsum/N[:,None]).ravel();long['Detection_fraction']=(det/N[:,None]).ravel();long['Review_min20_cells']=long.Cells.ge(20)
  all_specimen.append(long);long.to_csv(OUT/f'{obj}_specimen_marker_evidence.tsv.gz',sep='\t',index=False,compression='gzip')
  # Aggregate sums/counts to source x fine label, without equal-weighting tiny
  # specimen groups or fabricating donor independence.
  key=['Data_object','Dataset','Condition','Tissue','Lineage','celltype','celltype1','celltype2','Marker']
  long['Raw_sum']=long.Mean_raw_count*long.Cells;long['Log_sum']=long.Mean_log1p_10k*long.Cells;long['Detected_cells']=long.Detection_fraction*long.Cells
  z=long.groupby(key,dropna=False,observed=True).agg(Cells=('Cells','sum'),Specimens=('Sample','nunique'),Raw_sum=('Raw_sum',lambda v:v.sum(min_count=1)),Log_sum=('Log_sum',lambda v:v.sum(min_count=1)),Detected_cells=('Detected_cells',lambda v:v.sum(min_count=1))).reset_index()
  z['Mean_raw_count']=z.Raw_sum/z.Cells;z['Mean_log1p_10k']=z.Log_sum/z.Cells;z['Detection_fraction']=z.Detected_cells/z.Cells;z['Review_min20_cells']=z.Cells.ge(20);all_source.append(z)
  z.to_csv(OUT/f'{obj}_source_marker_evidence.tsv',sep='\t',index=False)
  qa.append({'Data_object':obj,'path':str(path),'raw_cells':len(ids),'full_RNA_genes':len(genes),'target_cells':len(sel),'audited_target_cells':seen,'raw_rows_loaded':raw_rows_loaded,'fine_labels':sel.celltype2.nunique(),'specimens':sel.Sample.nunique(),'available_markers':len(present),'missing_markers':';'.join(g for g in MARKERS if g not in ix),'status':'all target cells exact joined and audited'})
 pd.DataFrame(qa).to_csv(OUT/'object_processing_QA.tsv',sep='\t',index=False)
pd.DataFrame(qa).to_csv(OUT/'object_processing_QA.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(OUT/'full_RNA_marker_gene_coverage.tsv',sep='\t',index=False)
pd.concat(all_source,ignore_index=True).to_csv(OUT/'all_source_fine_label_marker_evidence.tsv.gz',sep='\t',index=False,compression='gzip')
pd.concat(all_specimen,ignore_index=True).to_csv(OUT/'all_specimen_fine_label_marker_evidence.tsv.gz',sep='\t',index=False,compression='gzip')
report={'passed':True,'objects_schema_gene_axes_read':len(qa),'objects_counts_streamed':sum(x['target_cells']>0 for x in qa),'full_RNA_marker_coverage_rows':len(coverage),'expected_target_cells':len(targets),'audited_target_cells':sum(x['audited_target_cells']for x in qa),'fine_labels':targets.celltype2.nunique(),'T_NK_cells':int(targets.Lineage.eq('T_NK').sum()),'B_Plasma_cells':int(targets.Lineage.eq('B_Plasma').sum()),'elapsed_seconds':time.time()-start_time,'chunk_cells':2048,'marker_genes':40,'scope':config['scope'],'no_labels_or_raw_assets_modified':True,'source_registry_SHA256':hashlib.sha256(registry.read_bytes()).hexdigest()}
assert report['objects_schema_gene_axes_read']==len(items) and report['audited_target_cells']==len(targets)
(OUT/'marker_coherence_QA.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

"""Measure actual cell-level marker co-detection in source-conditional targets.

Includes original GC/plasmablast/NKT/ILC/MAIT labels and source fine groups
selected by the prior review. Only 30 marker columns are dense in small chunks.
No cells, raw files or production annotation labels are modified.
"""
from pathlib import Path
import argparse,sys,json,time
import h5py,numpy as np,pandas as pd
from marker_input_paths import resolve_h5_path,prepare_specimen_source,load_object_crosswalk,exact_left_join
ap=argparse.ArgumentParser(description=__doc__)
for arg in ['input-tables','registry','reader','crosswalk','marker-audit-dir','output']:ap.add_argument('--'+arg,type=Path,required=True)
ap.add_argument('--data-h5-root',type=Path,help='Explicit filesystem root for relative registry h5_paths.')
a=ap.parse_args();T=a.input_tables.resolve();O=a.output.resolve();O.mkdir(parents=True,exist_ok=True)
if O==T or T in O.parents:raise ValueError('Output must be separate from immutable inputs')
sys.path.insert(0,str(a.reader.resolve().parent));from released_rna import inspect_schema,validate_raw_csr_arrays
GENES=['CD3D','CD3E','CD3G','CD247','TRAC','TRDC','TRGC1','TRGC2','TRDV2','NKG7','GNLY','KLRD1','FCGR3A','MS4A1','CD79A','JCHAIN','IGJ','MZB1','SDC1','BCL6','AICDA','MKI67','TRAV1-2','KLRB1','ZBTB16','IGHA1','IGHG1','CD4','CD8A','CD8B']
rules={'genes':GENES,'predeclared_before_second_raw_scan':True,'scope':'All original GC/plasmablast/NKT/ILC/MAIT labels plus source fine groups selected by prior review; not the whole 63-label second pass.','T_complex':'CD3D/CD3E/CD3G >=2 and all3 detection separately; CD247 is measured but not sufficient for T identity.','gamma_joint':'T-complex >=2 AND TRDC AND (TRGC1 or TRGC2), plus stricter all3 CD3 version; RNA review only, not confirmed gamma-delta identity.','plasma_joint':'at least2 of JCHAIN (explicit IGJ fallback), MZB1, SDC1 in same cell; all feature axes must be available.','marker_alias':'JCHAIN uses IGJ only when canonical JCHAIN absent; resolution is recorded explicitly, never summed with JCHAIN.','missing_feature':'Unavailable component makes its joint program NA, not zero.','thresholds':'No calibrated annotation accuracy or disease tests; actual detection only.','chunk_cells':2048}
(O/'predeclared_actual_joint_panel.json').write_text(json.dumps(rules,indent=2))
cross=load_object_crosswalk(a.crosswalk)
m=pd.read_csv(T/'all_release_original_metadata.tsv.gz',sep='\t',usecols=['OldCells','Sample','Dataset','Data_object'],keep_default_na=False)
parts=[]
for lin in ['T_NK','B_Plasma']:
 d=pd.read_csv(T/f'{lin}_original_integration_metadata.tsv.gz',sep='\t',usecols=['OldCells','Sample','Dataset','celltype2'],keep_default_na=False);d['Lineage']=lin;d=exact_left_join(d,m,on=['OldCells','Sample','Dataset'],validate='one_to_one',context=f'{lin} integrated-to-release');d.Data_object=d.Data_object.replace(cross);parts.append(d)
d=pd.concat(parts,ignore_index=True);del m,parts
prior=pd.read_csv(a.marker_audit_dir/'source_review_candidates.tsv',sep='\t');selected=set(zip(prior.Data_object,prior.celltype2))
special=d.celltype2.str.contains('NKT_|^ILC$|^MAIT$|^Bgc_|^Plasmablast$',regex=True)
select=special|np.array([(obj,label)in selected for obj,label in zip(d.Data_object,d.celltype2)])
targets=d[select].copy();del d;assert not targets[['Data_object','OldCells']].duplicated().any();expected=len(targets)
source=prepare_specimen_source(pd.read_csv(T/'composition_specimen_denominators_and_eligibility.tsv',sep='\t',keep_default_na=False),cross)
items=json.loads(a.registry.read_text())['datasets'];qa=[];coverage=[];outsource=[];outspec=[];st=time.time()
for item in items:
 path=resolve_h5_path(item['h5_paths']['RNA'],a.data_h5_root);obj=path.stem;sel=targets[targets.Data_object.eq(obj)].copy()
 if sel.empty:continue
 schema=inspect_schema(path);print(obj,len(sel),flush=True)
 with h5py.File(path,'r')as h:
  ids=np.asarray(h['names_obs'].asstr()[:]);genes=np.asarray(h['names_var'].asstr()[:]);ix={g:i for i,g in enumerate(genes)};present=[g for g in GENES if g in ix];pidx=[ix[g]for g in present];pi={g:i for i,g in enumerate(present)}
  jchain='JCHAIN'if'JCHAIN'in ix else'IGJ'if'IGJ'in ix else None
  for g in GENES:coverage.append({'Data_object':obj,'Gene':g,'full_RNA_available':g in ix,'axis_index':ix.get(g,np.nan),'effective_JCHAIN_feature':jchain or'unavailable'})
  aligned=sel.set_index('OldCells').reindex(ids);keep=aligned.Lineage.notna().to_numpy();assert keep.sum()==len(sel);group=aligned.loc[keep,['Lineage','celltype2','Sample']].reset_index(drop=True);codes,uni=pd.factorize(pd.MultiIndex.from_frame(group));names=pd.DataFrame(uni.tolist(),columns=group.columns);G=len(uni);N=np.bincount(codes,minlength=G);rowcodes=np.full(len(ids),-1,np.int32);rowcodes[keep]=codes
  specs={'CD3D_TRDC_actual_joint':['CD3D','TRDC'],'CD3D_TRDV2_actual_joint':['CD3D','TRDV2'],'T_CD3_DEG_any':['CD3D','CD3E','CD3G'],'T_CD3_DEG_2of3':['CD3D','CD3E','CD3G'],'T_CD3_DEG_all3':['CD3D','CD3E','CD3G'],'TRD_C_gamma_any':['TRDC','TRGC1','TRGC2'],'T_NK_actual_joint':['CD3D','CD3E','CD3G','NKG7','GNLY','KLRD1'],'gamma_program_joint':['CD3D','CD3E','CD3G','TRDC','TRGC1','TRGC2'],'gamma_program_strict_joint':['CD3D','CD3E','CD3G','TRDC','TRGC1','TRGC2'],'plasma_2of3_joint':[jchain,'MZB1','SDC1'],'MZB1_SDC1_actual_joint':['MZB1','SDC1'],'B_plasma_actual_joint':['MS4A1','CD79A',jchain,'MZB1','SDC1'],'B_GC_actual_joint':['MS4A1','CD79A','BCL6','AICDA'],'B_GC_cycle_actual_joint':['MS4A1','CD79A','BCL6','AICDA','MKI67'],'MAIT_TCR_RNA_actual_joint':['CD3D','CD3E','CD3G','TRAV1-2','KLRB1']}
  available={k:all(g in pi for g in gg)for k,gg in specs.items()};sums={k:np.zeros(G,np.int64)for k in specs};gene_det=np.zeros((G,len(GENES)),np.int64);raw=np.zeros((G,len(GENES)),np.float64);logs=raw.copy();seen=0;off=0
  for block in sorted(h['assay/RNA/layers/rawdata']):
   b=h['assay/RNA/layers/rawdata'][block];nr,nf=map(int,b.attrs['shape']);ptr=np.asarray(b['indptr'])
   for rs in range(0,nr,2048):
    re=min(nr,rs+2048);code=rowcodes[off+rs:off+re];take=code>=0
    if not take.any():continue
    p=ptr[rs:re+1];lo,hi=int(p[0]),int(p[-1]);x=validate_raw_csr_arrays(np.asarray(b['data'][lo:hi]),np.asarray(b['indices'][lo:hi]),p-lo,(re-rs,nf));lib=np.asarray(x.sum(axis=1,dtype=np.float64)).ravel()[take];assert(lib>0).all();z=x[take][:,pidx].toarray().astype(np.float64);yes=z>0;log=np.log1p(z/lib[:,None]*10000);cs=code[take]
    for g,j in pi.items():
     k=GENES.index(g);np.add.at(gene_det[:,k],cs,yes[:,j]);np.add.at(raw[:,k],cs,z[:,j]);np.add.at(logs[:,k],cs,log[:,j])
    def has(g):return yes[:,pi[g]]if g in pi else np.zeros(len(cs),bool)
    t=has('CD3D').astype(int)+has('CD3E')+has('CD3G');nk=has('NKG7')|has('GNLY')|has('KLRD1');trd=has('TRDC')&(has('TRGC1')|has('TRGC2'));plasma=has(jchain).astype(int)+has('MZB1')+has('SDC1');bc=has('MS4A1')|has('CD79A');gc=has('BCL6')|has('AICDA')
    v={'CD3D_TRDC_actual_joint':has('CD3D')&has('TRDC'),'CD3D_TRDV2_actual_joint':has('CD3D')&has('TRDV2'),'T_CD3_DEG_any':t>=1,'T_CD3_DEG_2of3':t>=2,'T_CD3_DEG_all3':t==3,'TRD_C_gamma_any':has('TRDC')|has('TRGC1')|has('TRGC2'),'T_NK_actual_joint':(t>=2)&nk,'gamma_program_joint':(t>=2)&trd,'gamma_program_strict_joint':(t==3)&trd,'plasma_2of3_joint':plasma>=2,'MZB1_SDC1_actual_joint':has('MZB1')&has('SDC1'),'B_plasma_actual_joint':bc&(plasma>=2),'B_GC_actual_joint':bc&gc,'B_GC_cycle_actual_joint':bc&gc&has('MKI67'),'MAIT_TCR_RNA_actual_joint':(t>=2)&has('TRAV1-2')&has('KLRB1')}
    for k in specs:
     if available[k]:np.add.at(sums[k],cs,v[k])
    seen+=int(take.sum())
   off+=nr
  assert seen==len(sel) and off==len(ids);names['Cells']=N;names['Data_object']=obj;names=exact_left_join(names,source,on=['Data_object','Sample'],validate='many_to_one',context=f'{obj} specimen source');names['effective_JCHAIN_feature']=jchain or'unavailable'
  metrics={}
  for k in specs:metrics[k+'_Detected_cells']=sums[k]if available[k]else np.full(G,np.nan);metrics[k+'_Fraction']=sums[k]/N if available[k]else np.full(G,np.nan);metrics[k+'_Available']=np.full(G,available[k])
  for g in GENES:
   k=GENES.index(g);metrics[g+'_Detection_fraction']=gene_det[:,k]/N if g in pi else np.full(G,np.nan);metrics[g+'_Mean_log1p10k']=logs[:,k]/N if g in pi else np.full(G,np.nan)
  names=pd.concat([names,pd.DataFrame(metrics)],axis=1);names['Review_min20_cells']=names.Cells.ge(20);outspec.append(names);names.to_csv(O/f'{obj}_actual_joint_specimen.tsv.gz',sep='\t',index=False,compression='gzip')
  # Weighted exact totals, not means of specimen fractions.
  q=names.groupby(['Data_object','Dataset','Condition','Tissue','Lineage','celltype2','effective_JCHAIN_feature'],dropna=False).agg(Cells=('Cells','sum'),Specimens=('Sample','nunique'),**{k+'_Detected_cells':(k+'_Detected_cells',lambda v:v.sum(min_count=1))for k in specs}).reset_index()
  for k in specs:q[k+'_Fraction']=q[k+'_Detected_cells']/q.Cells;q[k+'_Available']=available[k]
  for g in GENES:
   k=GENES.index(g);tmp=names[['Lineage','celltype2','Cells']].copy();tmp['D']=gene_det[:,k]if g in pi else np.nan;tmp['L']=logs[:,k]if g in pi else np.nan;tot=tmp.groupby(['Lineage','celltype2'],dropna=False).agg(D=('D',lambda v:v.sum(min_count=1)),L=('L',lambda v:v.sum(min_count=1))).reset_index();q=q.merge(tot,on=['Lineage','celltype2'],validate='one_to_one');q[g+'_Detection_fraction']=q.D/q.Cells;q[g+'_Mean_log1p10k']=q.L/q.Cells;q=q.drop(columns=['D','L'])
  q['Review_min20_cells']=q.Cells.ge(20);outsource.append(q);qa.append({'Data_object':obj,'selected_cells':seen,'selected_fine_groups':len(q),'effective_JCHAIN_feature':jchain or'unavailable'})
S=pd.concat(outsource,ignore_index=True);P=pd.concat(outspec,ignore_index=True);S.to_csv(O/'actual_joint_source_evidence.tsv',sep='\t',index=False);P.to_csv(O/'actual_joint_specimen_evidence.tsv.gz',sep='\t',index=False,compression='gzip');pd.DataFrame(coverage).to_csv(O/'actual_joint_gene_coverage.tsv',sep='\t',index=False);pd.DataFrame(qa).to_csv(O/'actual_joint_object_QA.tsv',sep='\t',index=False)
assert S.Cells.sum()==expected and P.Cells.sum()==expected
report={'passed':True,'selected_cells':expected,'streamed_cells':int(S.Cells.sum()),'selected_source_fine_groups':len(S),'specimen_fine_groups':len(P),'objects_counts_streamed':len(qa),'original_plasmablast_cells':int(S[S.celltype2.eq('Plasmablast')].Cells.sum()),'marker_genes':len(GENES),'elapsed_seconds':time.time()-st,'actual_joint_not_marginal_bound':True,'no_labels_or_raw_files_modified':True,'scope':rules['scope']};(O/'actual_joint_QA.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

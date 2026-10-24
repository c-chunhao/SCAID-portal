#!/usr/bin/env python3
"""Recompute GSE332623 RNA analysis from the corrected, QC-retained counts only."""
from pathlib import Path
import os,sys,json,hashlib,time,argparse,gc,subprocess,importlib.metadata as im
import numpy as np,pandas as pd,scipy.sparse as sp,h5py,anndata as ad
ROOT=Path('@PROJECT@/00.corrected_analysis/GSE332623')
WORK=Path('@WORK@/analysis')
OBJECT='Atopic_Dermatitis/GSE332623/PBMC';SEED=2024
PROTECTED=['CD3D','CD3E','TRAC','CD4','CD8A','CD8B','CCR7','IL7R','LTB','TCF7','NKG7','GNLY','KLRD1','FCGR3A','GZMB','PRF1','MS4A1','CD79A','CD74','CD37','CD19','CD27','MZB1','JCHAIN','XBP1','CD38','IGKC','IGLC2','CD14','LYZ','LST1','S100A8','S100A9','FCN1','MS4A7','CTSS','FCER1A','CLEC10A','CD1C','CLEC9A','XCR1','TCF4','IL3RA','LILRA4','SERPINF1','CST3','FCER1G','LILRB1','PPBP','PF4','GP9','HBB','HBA1','HBA2','TRDC','TRGC1','TRGC2','TRAV1-2','KLRB1','SLC4A10','MKI67','TOP2A','STMN1','FCGR3B','CSF3R','AHSP','VWF','PECAM1','EPCAM','KRT8','KRT18','COL1A1','DCN']
STATE_PANELS={'IFN':['ISG15','IFI6','IFIT1','IFIT3','MX1','OAS1','RSAD2','IRF7','STAT1','ISG20'],'Cytotoxic':['NKG7','GNLY','PRF1','GZMB','GZMH','CCL5'],'Immediate_early':['FOS','JUN','JUNB','EGR1','IER2'],'Cell_cycle':['MKI67','TOP2A','STMN1','TYMS','PCNA']}
LINEAGE_PROTECTED=PROTECTED.copy()
PROTECTED+=list(dict.fromkeys(g for genes in STATE_PANELS.values()for g in genes))
PANELS={'T cells':['CD3D','CD3E','TRAC'],'NK cells':['NKG7','GNLY','KLRD1'],'B cells':['CD79A','MS4A1','CD37','CD19'],'Plasma cells':['MZB1','JCHAIN','XBP1'],'Monocytes':['LYZ','LST1','FCN1','S100A8','S100A9','MS4A7'],'Dendritic cells':['FCER1A','CD1C','CLEC10A','CLEC9A'],'pDC':['LILRA4','IL3RA','TCF4'],'Platelets':['PF4','PPBP','GP9'],'Erythroid':['HBB','HBA1','HBA2','AHSP'],'Neutrophils':['FCGR3B','CSF3R','S100A8'],'Mast cells':['KIT','TPSAB1','TPSB2'],'Endothelial cells':['PECAM1','VWF','KDR'],'Epithelial cells':['EPCAM','KRT8','KRT18'],'Fibroblasts':['COL1A1','DCN','COL3A1']}
ADT_PANELS={'T cells':['Hu.CD3_UCHT1'],'NK cells':['Hu.CD56','Hu.CD16'],'B cells':['Hu.CD19','Hu.CD20_2H7'],'Plasma cells':['Hu.CD38_HIT2'],'Monocytes':['Hu.CD14_M5E2','Hu.CD16'],'Dendritic cells':['Hu.CD1c','Hu.CD11c'],'pDC':['Hu.CD123','Hu.CD303'],'Platelets':['Hu.CD41','Hu.CD42b']}

assert PANELS['T cells']==['CD3D','CD3E','TRAC'] and ADT_PANELS['T cells']==['Hu.CD3_UCHT1']

FINE_PANELS={
 'Naive B cells':['TCL1A','IGHD','IGHM','FCER2'], 'Memory B cells':['CD27','TNFRSF13B','CD44'], 'Age-associated B cells':['ITGAX','TBX21','FCRL5','FCRL3'], 'Transitional B cells':['TCL1A','MME','CD24','CD38'],
 'Plasma cells':['MZB1','JCHAIN','XBP1','SDC1'], 'Plasmablasts':['MZB1','JCHAIN','XBP1','SDC1'],
 'Classical monocytes':['CD14','FCN1','S100A8','S100A9','VCAN'], 'Non-classical monocytes':['FCGR3A','MS4A7','LILRB1','IFITM3'],
 'CD16+ NK cells':['FCGR3A','NKG7','GNLY','PRF1'], 'CD16- NK cells':['XCL1','XCL2','NCAM1','KLRC1'],
 'Tcm/Naive helper T cells':['CCR7','IL7R','TCF7','LEF1'], 'Tcm/Naive cytotoxic T cells':['CCR7','TCF7','LEF1','CD8A','CD8B'],
 'Tem/Effector helper T cells':['IL7R','GZMK','CCL5','IL32'], 'Tem/Effector helper T cells PD1+':['IL7R','GZMK','PDCD1'],
 'Tem/Temra cytotoxic T cells':['FGFBP2','GZMH','GZMB','PRF1'], 'Memory CD4+ cytotoxic T cells':['CD4','GZMK','CCL5','NKG7'],
 'MAIT cells':['SLC4A10','TRAV1-2','KLRB1','ZBTB16'], 'gamma-delta T cells':['TRDC','TRGC1','TRGC2'], 'CRTAM+ gamma-delta T cells':['TRDC','TRGC1','TRGC2','CRTAM'],
 'Regulatory T cells':['FOXP3','IL2RA','CTLA4','TIGIT'],
 'Type 1 helper T cells':['TBX21','CXCR3','IFNG'], 'Type 17 helper T cells':['RORC','CCR6','KLRB1','IL23R'], 'Follicular helper T cells':['CXCR5','PDCD1','ICOS','MAF'],
 'DC1':['CLEC9A','XCR1','CADM1','BATF3'], 'DC2':['CD1C','FCER1A','CLEC10A','CD1E'], 'DC3':['CD1C','CD14','CD163'], 'pDC':['LILRA4','IL3RA','TCF4','GZMB'],
 'DC precursor':['AXL','SIGLEC6','CD5'], 'Migratory DCs':['CCR7','FSCN1','LAMP3'], 'Transitional DC':['AXL','SIGLEC6','IRF8'],
 'Megakaryocytes/platelets':['PF4','PPBP','GP9','ITGA2B'], 'Neutrophils':['FCGR3B','CSF3R','S100A8','S100A9'], 'Mast cells':['TPSAB1','TPSB2','KIT','MS4A2']}
for _fine in ['Cycling B cells','Cycling T cells','Cycling NK cells','Cycling monocytes','Cycling DCs','Cycling gamma-delta T cells']:FINE_PANELS[_fine]=['MKI67','TOP2A','STMN1']

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb')as f:
  for z in iter(lambda:f.read(16*1024**2),b''):h.update(z)
 return h.hexdigest()
def axis_sha(names):return hashlib.sha256(('\n'.join(map(str,names))+'\n').encode()).hexdigest()
def write_json(path,value):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,default=lambda x:x.item()if isinstance(x,np.generic)else str(x)));os.replace(tmp,p)
def prune_none(value):
 # Older public AnnData readers do not understand the newer null encoding.
 if isinstance(value,dict):return {k:prune_none(v)for k,v in value.items()if v is not None}
 if isinstance(value,(list,tuple)):return [prune_none(v)for v in value if v is not None]
 return value

def validate_legacy_h5ad(rna,adt,out):
 oldpy='@WEB_PYTHON@'
 code="""import anndata as ad,numpy as np,json,sys,hashlib
r=ad.read_h5ad(sys.argv[1],backed='r');d=ad.read_h5ad(sys.argv[2],backed='r')
assert r.obs_names.equals(d.obs_names) and r.var_names.is_unique and d.var_names.is_unique
assert r.obsm['X_umap'].shape==(r.n_obs,2) and np.isfinite(r.obsm['X_umap']).all()
assert np.array_equal(r.obsm['X_umap'],d.obsm['X_umap'])
assert r.layers['counts'].dtype.kind in 'iu' and d.layers['counts'].dtype.kind in 'iu'
assert set(r.obs.Group.astype(str))=={'Unresolved'} and set(r.obs.source_condition.astype(str))=={'Unresolved'}
assert np.asarray(r.uns['scaid']['object_id']).ravel().tolist()==['Atopic_Dermatitis/GSE332623/PBMC']
report={'status':'validated','anndata_version':ad.__version__,'reader':sys.executable,'RNA_path':sys.argv[1],'ADT_path':sys.argv[2],'RNA_shape':list(r.shape),'ADT_shape':list(d.shape),'axes_equal':True,'UMAP_equal_finite':True,'counts_integer_dtype':True,'clinical_status':'Unresolved','log_transform':'Natural log1p; base=None key omitted for cross-version compatibility'}
open(sys.argv[3],'w').write(json.dumps(report,indent=2));r.file.close();d.file.close();print(json.dumps(report))
"""
 path=out/'provenance/legacy_AnnData_readback.json';run=subprocess.run([oldpy,'-c',code,str(rna),str(adt),str(path)],text=True,capture_output=True)
 (out/'provenance/legacy_AnnData_readback.log').write_text(run.stdout+'\n'+run.stderr)
 if run.returncode:raise RuntimeError('Legacy AnnData readback failed: '+run.stderr[-4000:])
 return json.loads(path.read_text())

def setup(out):
 for x in ['checkpoints','plots','annotation_evidence','labels','celltypist_chunks','provenance']: (out/x).mkdir(parents=True,exist_ok=True)
def strings(x):
 a=np.asarray(x);return np.asarray([z.decode()if isinstance(z,(bytes,np.bytes_))else str(z)for z in a.ravel()],dtype=object)
def read_bridge(path):
 with h5py.File(path,'r')as f:
  mats={};features={}
  for assay in ['RNA','ADT']:
   g=f[assay];shape=tuple(np.asarray(g['shape']).ravel().astype(int));v=np.asarray(g['data']);i=np.asarray(g['indices'],dtype=np.int32);p=np.asarray(g['indptr'],dtype=np.int64)
   assert len(shape)==2 and len(p)==shape[0]+1 and p[0]==0 and p[-1]==len(v)==len(i)
   assert np.all(np.isfinite(v))and np.all(v>=0)and np.all(v==np.rint(v))and (not len(v)or v.max()<np.iinfo(np.int32).max)
   x=sp.csr_matrix((v.astype(np.int32),i,p),shape=shape);x.sum_duplicates();x.sort_indices();mats[assay]=x
   features[assay]=pd.DataFrame({k:strings(g[k])for k in ['feature_ids','feature_names','source_symbols']})
   assert len(features[assay])==shape[1]
  obs=pd.DataFrame({k:strings(v)if v.dtype.kind in 'OSU'else np.asarray(v).ravel()for k,v in f['obs'].items()})
 assert len(obs)==mats['RNA'].shape[0]==mats['ADT'].shape[0]
 assert {'cell_id','library_id','capture_id','pool_id','source_gsm'}.issubset(obs.columns)
 assert not obs.cell_id.duplicated().any();obs.index=pd.Index(obs.pop('cell_id').astype(str),name='cell_id')
 assert np.all(np.asarray(mats['RNA'].sum(axis=1)).ravel()>0)
 for c in ['Group','source_condition']:
  if c in obs:assert set(obs[c].astype(str))=={'Unresolved'}
 return mats,features,obs

def coarse(label):
 s=str(label)
 if s in ['B cells','B-cell lineage']or ' B cells'in s or s.startswith(('Naive B','Memory B','Transitional B','Follicular B','Age-associated B')):return 'B cells'
 if s in ['Plasma cells','Plasmablasts']:return 'Plasma cells'
 if 'NK'in s and 'NKT'not in s:return 'NK cells'
 if s=='T cells'or ' T cells'in s or s.startswith(('CD8','Tcm/','Tem/','Trm ','Treg(','T(agonist)','Type 1 helper','Type 17 helper','Regulatory T','Follicular helper','MAIT','NKT','gamma-delta')):return 'T cells'
 if s in ['Monocytes','Classical monocytes','Non-classical monocytes','Cycling monocytes','Monocyte precursor','Mono-mac','MNP','Macrophages']:return 'Monocytes'
 if s.startswith('pDC'):return 'pDC'
 if s in ['DC','DC1','DC2','DC3','Cycling DCs','Migratory DCs','Transitional DC','DC precursor']:return 'Dendritic cells'
 if s in ['Megakaryocytes/platelets','Megakaryocyte precursor','Early MK']:return 'Platelets'
 if any(z in s.lower()for z in ['erythroid','erythrocyte']):return 'Erythroid'
 if s in ['Neutrophils','Granulocytes']:return 'Neutrophils'
 if s in ['Mast cells','Endothelial cells','Epithelial cells','Fibroblasts']:return s
 return 'Unassigned'

def assemble(out,bridges):
 setup(out);xs=[];aa=[];obslist=[];base=None;abase=None;inputs=[]
 for path in bridges:
  print('ASSEMBLE',path,flush=True);m,f,o=read_bridge(path)
  if base is None:base=f['RNA'];abase=f['ADT']
  else:assert base.equals(f['RNA'])and abase.equals(f['ADT'])
  xs.append(m['RNA']);aa.append(m['ADT']);obslist.append(o);inputs.append({'bridge':str(path),'sha256':sha(path),'n_cells':len(o)})
 counts=sp.vstack(xs,format='csr');adt=sp.vstack(aa,format='csr');obs=pd.concat(obslist);del xs,aa,obslist;gc.collect()
 assert not obs.index.duplicated().any()
 for alias,source in [('source_library_id','library_id'),('source_pool','pool_id'),('source_rna_gsm','source_gsm')]:
  if alias in obs:assert np.array_equal(obs[alias].astype(str),obs[source].astype(str))
  else:obs[alias]=obs[source].copy()
 capture=pd.read_csv(ROOT/'capture_manifest.tsv',sep='\t').set_index('library_id');obs['source_adt_gsm']=obs.library_id.map(capture.adt_gsm).to_numpy();assert not pd.isna(obs.source_adt_gsm).any()
 symbols=base.source_symbols.astype(str);assert symbols.str.len().min()>0
 unique=pd.Index(pd.unique(symbols),name='gene_symbol');mapping=unique.get_indexer(symbols)
 before=np.asarray(counts.sum(axis=1),dtype=np.int64).ravel()
 transform=sp.csr_matrix((np.ones(len(mapping),dtype=np.int32),(np.arange(len(mapping)),mapping)),shape=(len(mapping),len(unique)))
 counts=(counts@transform).tocsr();counts.sum_duplicates();counts.sort_indices();after=np.asarray(counts.sum(axis=1),dtype=np.int64).ravel();assert np.array_equal(before,after)
 base.assign(symbol_index=mapping).to_csv(out/'gene_feature_to_symbol.tsv',sep='\t',index=False)
 var=pd.DataFrame(index=unique);var['source_symbol']=unique
 ids=base.groupby('source_symbols',sort=False).feature_ids.agg(lambda v:';'.join(v.astype(str))).reindex(unique)
 var['source_feature_ids']=ids.to_numpy();var['n_source_features']=np.bincount(mapping,minlength=len(unique));var['feature_type']='Gene Expression'
 for c in ['Group','source_condition','Disease_detail','Disease_subtype']:obs[c]='Unresolved'
 obs['scaid_object_id']=OBJECT;obs['scaid_source']='GSE332623';obs['scaid_tissue']='PBMC';obs['scaid_cohort']='mixed AD/HC study; cell-level clinical status unresolved'
 obs['n_genes_symbol']=np.diff(counts.indptr);obs['n_UMI_RNA_symbol']=after
 for c in obs:
  if obs[c].dtype.kind in 'OUS':obs[c]=obs[c].astype(str).astype('category')
 x=counts.astype(np.float32);sums=np.asarray(x.sum(axis=1)).ravel();x.data*=np.repeat((1e4/sums).astype(np.float32),np.diff(x.indptr));np.log1p(x.data,out=x.data);x.sort_indices()
 a=ad.AnnData(x,obs=obs,var=var);a.layers['counts']=counts
 a.uns['log1p']={'base':None};a.uns['scaid']={'object_id':np.asarray([OBJECT],dtype=str),'cohort':'mixed AD and matched healthy controls; unresolved per-cell donor/clinical mapping','observations':'QC-retained library-barcode records','normalization':'RNA exact-symbol sum; total10000; natural log1p','counts_origin':'Published MEX RNA integer counts, QC-retained; exact same-symbol features aggregated','annotation_status':'not yet annotated'}
 a.write_h5ad(out/'checkpoints/normalized_RNA.h5ad',compression='gzip',compression_opts=1)
 assert abase.source_symbols.is_unique
 avar=pd.DataFrame({'source_feature_id':abase.feature_ids.to_numpy(),'source_symbol':abase.source_symbols.to_numpy(),'feature_type':'Antibody Capture'},index=pd.Index(abase.source_symbols.astype(str),name='antibody'))
 ax=adt.astype(np.float32);sparse_log=ax.copy();np.log1p(sparse_log.data,out=sparse_log.data);denom=np.exp(np.asarray(sparse_log.sum(axis=1)).ravel()/ax.shape[1]);ax.data/=np.repeat(denom.astype(np.float32),np.diff(ax.indptr));np.log1p(ax.data,out=ax.data)
 b=ad.AnnData(ax,obs=obs.copy(),var=avar);b.layers['counts']=adt;b.uns['normalization']='Per-cell CLR: log1p(count / exp(mean(log1p(count)))) across all140ADT features including zeros';b.uns['scaid']={'object_id':np.asarray([OBJECT],dtype=str),'cohort':'mixed AD/HC; clinical status unresolved'};b.write_h5ad(out/'ADT_normalized.h5ad',compression='gzip',compression_opts=1)
 report={'status':'validated','n_cells':a.n_obs,'n_RNA_source_features':len(base),'n_RNA_unique_symbols':a.n_vars,'n_ADT_features':b.n_vars,'cell_axis_sha256':axis_sha(a.obs_names),'RNA_feature_axis_sha256':axis_sha(a.var_names),'RNA_cell_UMI_preserved':True,'sources':inputs,'normalization':a.uns['scaid']['normalization'],'ADT_axis_equal':a.obs_names.equals(b.obs_names)}
 write_json(out/'checkpoints/assembly.json',report);print('ASSEMBLY_READY',json.dumps({k:v for k,v in report.items()if k!='sources'}),flush=True)

def classify(out):
 import celltypist
 from celltypist.models import Model
 a=ad.read_h5ad(out/'checkpoints/normalized_RNA.h5ad',backed='r');result=[]
 for name in ['Immune_All_High','Immune_All_Low']:
  model_path=Path('@CELLTYPIST_MODELS@')/(name+'.pkl');model=Model.load(str(model_path));parts=[]
  overlap=sum(a.var_names.isin(model.features));assert overlap>=2000
  for start in range(0,a.n_obs,10000):
   end=min(a.n_obs,start+10000);path=out/'celltypist_chunks'/f'{name}_{start}_{end}.tsv.gz'
   if path.exists():p=pd.read_csv(path,sep='\t',index_col=0)
   else:
    chunk=ad.AnnData(a.X[start:end].tocsr().astype(np.float32),obs=pd.DataFrame(index=a.obs_names[start:end]),var=a.var.copy())
    pred=celltypist.annotate(chunk,model=model,majority_voting=False,mode='best match');prob=pred.probability_matrix.to_numpy();top=np.sort(prob,axis=1)[:,-2:]
    p=pd.DataFrame({'label':pred.predicted_labels.predicted_labels.astype(str),'probability':top[:,-1],'margin':top[:,-1]-top[:,-2]},index=chunk.obs_names);p.to_csv(path,sep='\t',compression='gzip');del chunk,pred,prob;gc.collect()
   assert list(p.index)==list(a.obs_names[start:end]);parts.append(p);print('CLASSIFY',name,end,'/',a.n_obs,flush=True)
  p=pd.concat(parts);p.columns=[('high_'if name.endswith('High')else'low_')+x for x in p.columns];result.append(p)
  write_json(out/'checkpoints'/f'{name}_model.json',{'path':str(model_path),'sha256':sha(model_path),'n_model_genes':len(model.features),'matched_genes':overlap,'classes':model.cell_types.tolist(),'method':'Per-cell logistic prediction in10000-cell chunks; no inherited labels; cluster majority added separately'})
 r=map_reference_predictions(pd.concat(result,axis=1))
 r.to_csv(out/'labels/reference_predictions.tsv.gz',sep='\t',compression='gzip');write_json(out/'checkpoints/reference_classification.json',reference_mapping_receipt(r));a.file.close()

def map_reference_predictions(r):
 r=r.copy();r['high_coarse']=r.high_label.map(coarse);r['low_coarse']=r.low_label.map(coarse)
 exact=r.high_coarse.eq(r.low_coarse)&r.high_coarse.ne('Unassigned')
 ilc_nk=r.high_label.eq('ILC')&r.low_coarse.eq('NK cells')
 compatible=exact|ilc_nk
 r['reference_agreement']=exact;r['reference_compatible']=compatible
 r['reference_rule']=np.select([exact,ilc_nk,r.high_label.eq('Cycling cells')],['exact_coarse_agreement','high_ILC_low_NK_compatible','high_Cycling_state_only'],default='incompatible_references')
 r['high_cycling_state']=r.high_label.eq('Cycling cells');r['reference_coarse']=np.where(compatible,r.low_coarse,'Unassigned')
 r['reference_confident']=compatible&(r.low_probability>=.5)&(r.high_probability>=.5)
 return r

def reference_mapping_receipt(r):
 return {'status':'validated','n_cells':len(r),'cell_axis_sha256':axis_sha(r.index),'fraction_reference_agreement':float(r.reference_agreement.mean()),'fraction_reference_compatible':float(r.reference_compatible.mean()),'fraction_reference_confident':float(r.reference_confident.mean()),'mapping_policy_version':2,'reference_rule_counts':r.reference_rule.value_counts().to_dict(),'mapping_rule':'Agreement means exact coarse equality. High ILC is compatible only with Low NK; explicitly distinguished from exact agreement. High Cycling cells denotes a state and is not used as independent lineage agreement.'}

def representation_diagnostics(pca,harmony,obs,pred,state_scores,out):
 from sklearn.neighbors import NearestNeighbors
 rng=np.random.default_rng(SEED);eligible=(pred.reference_compatible.astype(bool)&(pred.high_probability>=.5)).to_numpy();selected=[]
 library=obs.library_id.astype(str).to_numpy()
 for lib in pd.unique(library):
  idx=np.flatnonzero((library==lib)&eligible);selected.extend(rng.choice(idx,min(600,len(idx)),replace=False))
 selected=np.asarray(selected,dtype=int)
 if len(selected)<100: return {'status':'insufficient_confident_reference_cells','selected_representation':'PCA','n_diagnostic_cells':len(selected)}
 lab=pred.reference_coarse.to_numpy()[selected];fine=pred.low_label.astype(str).to_numpy()[selected];batch=library[selected];result={};states=state_scores.iloc[selected]
 for key,rep in [('PCA',pca),('Harmony',harmony)]:
  z=np.asarray(rep[selected],dtype=np.float32);nn=NearestNeighbors(n_neighbors=min(31,len(z)),algorithm='brute',n_jobs=8).fit(z).kneighbors(z,return_distance=False)[:,1:]
  purity=np.mean(lab[nn]==lab[:,None]);fine_purity=np.mean(fine[nn]==fine[:,None]);entropy=[]
  for j,nei in enumerate(nn):
   same=nei[lab[nei]==lab[j]];nb=len(pd.unique(batch[lab==lab[j]]))
   if len(same)<3 or nb<2:continue
   _,n=np.unique(batch[same],return_counts=True);p=n/n.sum();entropy.append(float(-np.sum(p*np.log(p))/np.log(nb)))
  from scipy.stats import spearmanr
  state_correlations={};state_details=[]
  for field in states:
   raw=states[field].to_numpy();smoothed=raw[nn].mean(axis=1);rho=[];weights=[]
   for lineage in pd.unique(lab):
    ii=np.flatnonzero(lab==lineage)
    if len(ii)<100 or np.std(raw[ii])<.02 or np.std(smoothed[ii])<1e-8:continue
    value=float(spearmanr(raw[ii],smoothed[ii]).statistic)
    if not np.isfinite(value):continue
    rho.append(value);weights.append(len(ii));state_details.append({'state':field,'lineage':lineage,'n_cells':len(ii),'neighbor_smoothing_spearman':value})
   if weights:state_correlations[field]=float(np.average(rho,weights=weights))
  result[key]={'RNA_state_neighbor_correlations':state_correlations,'RNA_state_lineage_details':state_details,'reference_coarse_neighbor_purity':float(purity),'reference_fine_neighbor_purity':float(fine_purity),'within_coarse_library_entropy':float(np.mean(entropy))if entropy else 0.0,'n_entropy_cells':len(entropy)}
 gain=result['Harmony']['within_coarse_library_entropy']-result['PCA']['within_coarse_library_entropy'];loss=result['PCA']['reference_coarse_neighbor_purity']-result['Harmony']['reference_coarse_neighbor_purity']
 fine_loss=result['PCA']['reference_fine_neighbor_purity']-result['Harmony']['reference_fine_neighbor_purity']
 state_losses={k:result['PCA']['RNA_state_neighbor_correlations'][k]-result['Harmony']['RNA_state_neighbor_correlations'][k]for k in result['PCA']['RNA_state_neighbor_correlations']if k in result['Harmony']['RNA_state_neighbor_correlations']};max_state_loss=max(state_losses.values(),default=0)
 chosen='Harmony'if gain>.02 and loss<=.02 and fine_loss<=.03 and max_state_loss<=.05 and result['Harmony']['reference_coarse_neighbor_purity']>=.9 else'PCA'
 return {'status':'evaluated','n_diagnostic_cells':len(selected),'sampling':'up to600confident reference-labeled cells/library; fixed seed2024','metrics':result,'within_coarse_entropy_gain':gain,'coarse_purity_loss':loss,'fine_purity_loss':fine_loss,'RNA_state_neighbor_correlation_loss':state_losses,'max_RNA_state_correlation_loss':max_state_loss,'selected_representation':chosen,'selection_rule':'Harmony only if entropy improves>0.02, coarse-neighbor purity decreases<=0.02 and remains>=0.90, fine-reference purity decreases<=0.03, and within-lineage IFN/cytotoxic/early-response/cell-cycle RNA-score neighborhood correlation decreases<=0.05; otherwise retain PCA. These technical/reference checks do not establish preservation of unknown donor/clinical effects.'}

def embed(out):
 import scanpy as sc,harmonypy as hm,torch
 sc.settings.n_jobs=8;torch.set_num_threads(8)
 a=ad.read_h5ad(out/'checkpoints/normalized_RNA.h5ad');pred=pd.read_csv(out/'labels/reference_predictions.tsv.gz',sep='\t',index_col=0);assert list(pred.index)==list(a.obs_names)
 sc.pp.highly_variable_genes(a,n_top_genes=min(3000,a.n_vars-1),flavor='seurat',batch_key='library_id',subset=False)
 nuisance=np.asarray(a.var_names.str.match(r'^(MT-|RPS[0-9]+[A-Z]?$|RPS4[XY][0-9]?$|RPL[0-9]+[A-Z]?$|RPLP[0-9]+$)'));protected=np.asarray(a.var_names.isin(PROTECTED));mask=(a.var.highly_variable.to_numpy()&~nuisance)|protected
 a.var['protected_lineage_marker']=np.asarray(a.var_names.isin(LINEAGE_PROTECTED));a.var['protected_state_marker']=np.asarray(a.var_names.isin([g for p in STATE_PANELS.values()for g in p]));a.var['protected_biological_marker']=protected;a.var['integration_feature']=mask
 a.var.to_csv(out/'checkpoints/HVG_and_feature_selection.tsv',sep='\t')
 state_scores=pd.DataFrame(index=a.obs_names)
 for field,genes in STATE_PANELS.items():
  present=[g for g in genes if g in a.var_names]
  if len(present)>=3:state_scores[field]=np.asarray(a.X[:,a.var_names.get_indexer(present)].mean(axis=1)).ravel()
 state_scores.to_csv(out/'annotation_evidence/diagnostic_RNA_state_mean_log_expression.tsv.gz',sep='\t',compression='gzip')
 write_json(out/'annotation_evidence/RNA_state_diagnostic_methods.json',{'panels':STATE_PANELS,'definition':'Mean log1p(total-normalized RNA) over available panel genes, including zeros; only an embedding-preservation diagnostic, not a disease activity test or a substitute for AUCell/UCell/singscore','sampling':'Up to600coarse-reference-agreement and High-model-probability>=0.5 cells per library; fixed seed2024','state_preservation':'Within-lineage Spearman correlation between each original RNA-score and mean neighbor score, requiring100cells and nontrivial variance; decreases>0.05 block Harmony selection','clinical_inference':False})
 hv=ad.AnnData(a.X[:,mask].copy(),obs=a.obs.copy(),var=a.var.loc[mask].copy());del a;gc.collect()
 sc.pp.scale(hv,zero_center=False,max_value=10);npcs=min(50,hv.n_vars-1,hv.n_obs-1)
 sc.pp.pca(hv,n_comps=npcs,zero_center=True,svd_solver='arpack',random_state=SEED)
 pca=np.asarray(hv.obsm['X_pca'],dtype=np.float32);print('PCA_READY',pca.shape,flush=True)
 if hv.obs.library_id.nunique()>1:
  meta=pd.DataFrame({'library_id':hv.obs.library_id.astype(str)},index=hv.obs_names)
  h=hm.run_harmony(pca,meta,['library_id'],theta=2,lamb=1,nclust=min(100,max(2,round(hv.n_obs/100))),max_iter_harmony=15,max_iter_kmeans=20,random_state=SEED,device='cpu',verbose=True)
  z=h.Z_corr
  if hasattr(z,'detach'):z=z.detach().cpu().numpy()
  z=np.asarray(z,dtype=np.float32);harmony=z.T if z.shape[0]==pca.shape[1]else z;assert harmony.shape==pca.shape and np.isfinite(harmony).all();del h;gc.collect()
  diag=representation_diagnostics(pca,harmony,hv.obs,pred,state_scores,out)
 else:harmony=pca.copy();diag={'status':'single_library_pilot','selected_representation':'PCA','reason':'No cross-library correction required for the one-library technical pilot.'}
 write_json(out/'annotation_evidence/PCA_Harmony_diagnostics.json',diag)
 hv.obsm['X_analysis']=harmony if diag['selected_representation']=='Harmony'else pca
 sc.pp.neighbors(hv,n_neighbors=min(30,hv.n_obs-1),use_rep='X_analysis',method='umap',transformer='pynndescent',random_state=SEED)
 sc.tl.leiden(hv,resolution=.8,key_added='leiden',flavor='igraph',n_iterations=2,directed=False,random_state=SEED)
 sc.tl.umap(hv,min_dist=.3,spread=1,maxiter=300,random_state=SEED)
 np.savez_compressed(out/'checkpoints/embeddings.npz',cell_ids=np.asarray(hv.obs_names,dtype=str),X_pca=pca,X_harmony=harmony,X_umap=np.asarray(hv.obsm['X_umap'],dtype=np.float32),leiden=hv.obs.leiden.astype(str).to_numpy(dtype=str),selected_representation=np.asarray([diag['selected_representation']]))
 sp.save_npz(out/'checkpoints/connectivities.npz',hv.obsp['connectivities']);sp.save_npz(out/'checkpoints/distances.npz',hv.obsp['distances'])
 pd.DataFrame(hv.varm['PCs'],index=hv.var_names).to_csv(out/'checkpoints/PCA_loadings.tsv.gz',sep='\t',compression='gzip')
 write_json(out/'checkpoints/embedding.json',{'status':'validated','n_cells':hv.n_obs,'n_integration_features':hv.n_vars,'n_protected_present':int(protected.sum()),'n_PCs':npcs,'neighbors':30,'leiden_resolution':.8,'n_clusters':int(hv.obs.leiden.nunique()),'UMAP_min_dist':.3,'selected_representation':diag['selected_representation'],'scVI_assessment':'Large size and27technical library partitions made scVI worth considering. The installed scvi import failed due NumPy StringDType incompatibility. A transparent PCA/Harmony workflow was selected without changing global dependencies; both embeddings retained and mixing/reference-lineage preservation measured. No latent model can identify unresolved donor/clinical effects.'})
 print('EMBEDDING_READY',hv.n_obs,hv.obs.leiden.nunique(),diag['selected_representation'],flush=True)

def expression_evidence(a,b,clusters,out):
 levels=pd.Index(pd.unique(clusters.astype(str)));code=levels.get_indexer(clusters.astype(str));indicator=sp.csr_matrix((np.ones(a.n_obs,dtype=np.float32),(code,np.arange(a.n_obs))),shape=(len(levels),a.n_obs));sizes=np.bincount(code,minlength=len(levels))
 rgenes=[g for p in list(PANELS.values())+list(FINE_PANELS.values())for g in p if g in a.var_names];rgenes=list(dict.fromkeys(rgenes));rx=a.X[:,a.var_names.get_indexer(rgenes)].tocsr();rmean=np.asarray((indicator@rx).toarray())/sizes[:,None];rpct=np.asarray((indicator@(rx>0)).toarray())/sizes[:,None]
 adtgenes=[g for p in ADT_PANELS.values()for g in p if g in b.var_names];adtgenes=list(dict.fromkeys(adtgenes));ax=b.X[:,b.var_names.get_indexer(adtgenes)].tocsr();amean=np.asarray((indicator@ax).toarray())/sizes[:,None]
 rows=[]
 for i,c in enumerate(levels):
  for j,g in enumerate(rgenes):rows.append({'cluster':c,'modality':'RNA','marker':g,'mean_log_normalized':rmean[i,j],'pct_expressing':rpct[i,j],'n_cells':sizes[i]})
  for j,g in enumerate(adtgenes):rows.append({'cluster':c,'modality':'ADT','marker':g,'mean_log_normalized':amean[i,j],'pct_expressing':np.nan,'n_cells':sizes[i]})
 table=pd.DataFrame(rows);table.to_csv(out/'annotation_evidence/RNA_ADT_marker_evidence.tsv',sep='\t',index=False)
 means=pd.DataFrame(rmean,index=levels,columns=rgenes);pct=pd.DataFrame(rpct,index=levels,columns=rgenes);support=pd.DataFrame(index=levels)
 for group,genes in PANELS.items():
  actual=[g for g in genes if g in means];support[group]=pct[actual].mean(axis=1)if actual else 0
 support.to_csv(out/'annotation_evidence/cluster_RNA_marker_detection.tsv',sep='\t')
 fine_support=pd.DataFrame(index=levels)
 for group,genes in FINE_PANELS.items():
  actual=[g for g in genes if g in pct]
  fine_support[group]=((pct[actual].mean(axis=1)>=.05)|(pct[actual].max(axis=1)>=.2)) if actual else False
 fine_support.to_csv(out/'annotation_evidence/fine_label_RNA_marker_support.tsv',sep='\t')
 # Within-library discrimination avoids interpreting technical-library shifts as protein support.
 # AUROC describes relative cluster enrichment, not absolute antibody positivity or expert validation.
 from scipy.stats import rankdata
 libraries=b.obs.library_id.astype(str).to_numpy();auc_rows=[]
 for lib in pd.unique(libraries):
  ix=np.flatnonzero(libraries==lib);local=code[ix];n=np.bincount(local,minlength=len(levels));n0=len(ix)-n;eligible=(n>=5)&(n0>=5)
  local_x=b.X[ix][:,b.var_names.get_indexer(adtgenes)].tocsc()
  for j,g in enumerate(adtgenes):
   values=np.asarray(local_x[:,j].toarray()).ravel();ranks=rankdata(values,method='average');rs=np.bincount(local,weights=ranks,minlength=len(levels));u=rs-n*(n+1)/2
   for k in np.flatnonzero(eligible):auc_rows.append({'cluster':levels[k],'library_id':lib,'marker':g,'AUROC_cluster_vs_rest':float(u[k]/(n[k]*n0[k])),'n_cluster_cells':int(n[k]),'n_rest_cells':int(n0[k])})
 auc=pd.DataFrame(auc_rows);auc.to_csv(out/'annotation_evidence/ADT_within_library_cluster_AUROC.tsv.gz',sep='\t',index=False,compression='gzip')
 adt_support=pd.DataFrame(index=levels)
 for group,genes in ADT_PANELS.items():
  values=[]
  for c in levels:
   selected=auc[(auc.cluster==c)&auc.marker.isin(genes)] if len(auc) else auc
   marker_auc=[float(np.average(g.AUROC_cluster_vs_rest,weights=g.n_cluster_cells)) for _,g in selected.groupby('marker')] if len(selected) else []
   values.append(max(marker_auc) if marker_auc else np.nan)
  adt_support[group]=values
 adt_support.to_csv(out/'annotation_evidence/cluster_ADT_lineage_enrichment.tsv',sep='\t')
 write_json(out/'annotation_evidence/ADT_crosscheck_methods.json',{'method':'Mann-Whitney AUROC for each predeclared surface marker: Leiden cluster versus remaining QC-retained cells within each technical library; tied ranks averaged; at least5cells per side; per-marker AUROC averaged using cluster-cell counts, then strongest panel marker retained','interpretation':'Relative cluster enrichment; does not establish absolute antibody positivity, donor identity, disease status or external expert confirmation','flags':'AUROC>=0.65 relative support; 0.55–0.65 limited enrichment; <0.55 weak enrichment; no panel/no eligible comparison unavailable. These exploratory flags do not automatically erase RNA-supported labels.'})
 return support,table,fine_support,adt_support

def markers(a,groups,out,name,cap):
 import scanpy as sc
 rng=np.random.default_rng(SEED);idx=[]
 for group in pd.unique(groups):
  ix=np.flatnonzero(groups==group)
  if len(ix)>=20:idx.extend(rng.choice(ix,min(cap,len(ix)),replace=False))
 idx=np.asarray(sorted(idx),dtype=int);sub=ad.AnnData(a.X[idx].copy(),obs=pd.DataFrame({'group':pd.Categorical(groups[idx])},index=a.obs_names[idx]),var=a.var.copy())
 assert sub.obs.group.nunique()>=2
 sc.tl.rank_genes_groups(sub,'group',method='wilcoxon',use_raw=False,n_genes=sub.n_vars,pts=True,tie_correct=True)
 d=sc.get.rank_genes_groups_df(sub,group=None);d=d.rename(columns={'group':'group1','names':'gene','logfoldchanges':'avg_log2FC','pvals':'p_val','pvals_adj':'p_val_adj','pct_nz_group':'pct.1','pct_nz_reference':'pct.2'})
 for c in ['p_val','p_val_adj','avg_log2FC']:d[c]=pd.to_numeric(d[c],errors='coerce')
 bad=d[~np.isfinite(d.avg_log2FC)].copy();bad.to_csv(out/'annotation_evidence'/(Path(name).name+'.nonfinite_log2FC.tsv.gz'),sep='\t',index=False,compression='gzip')
 d=d[np.isfinite(d.avg_log2FC)&(d.p_val_adj<=.05)&(d.avg_log2FC.abs()>=.5)&(d['pct.1']>=.1)].copy();d['group2']='others'
 tested_groups=sorted(map(str,sub.obs.group.cat.categories));d['test_group']=';'.join(tested_groups);d['test_group_number']=len(tested_groups)
 cols=['p_val','avg_log2FC','pct.1','pct.2','p_val_adj','gene','group1','group2','test_group_number','test_group'];d[cols].to_csv(out/name,sep='\t',index=False)
 sampling=pd.DataFrame({'cell_id':sub.obs_names,'group':sub.obs.group.astype(str),'library_id':a.obs.library_id.iloc[idx].astype(str).to_numpy()});sampling.to_csv(out/'annotation_evidence'/(Path(name).name+'.sampling.tsv.gz'),sep='\t',index=False,compression='gzip')
 write_json(out/'annotation_evidence'/(Path(name).name+'.methods.json'),{'method':'Wilcoxon cell-type/cluster versus remaining sampled cells; Benjamini-Hochberg per group; tie correction','sampling':'fixed seed2024; up to'+str(cap)+'QC-retained cells/group; no donor or disease comparison','n_test_cells':sub.n_obs,'n_genes_tested':sub.n_vars,'group_counts':sub.obs.group.value_counts().to_dict(),'display_rules':'adjustedP<=0.05, abs(log2FC)>=0.5, pct.group>=0.1','n_rows':len(d),'n_nonfinite_log2FC_excluded':len(bad),'nonfinite_log2FC_policy':'Preserve original nonfinite estimates in the separate evidence table; never replace with zero','tested_groups':tested_groups,'p_value_scope':'Cells are observational units for annotation markers, not independent people or clinical case-control testing.'});return d

def finalize(out,pilot=False):
 a=ad.read_h5ad(out/'checkpoints/normalized_RNA.h5ad');b=ad.read_h5ad(out/'ADT_normalized.h5ad');pred=pd.read_csv(out/'labels/reference_predictions.tsv.gz',sep='\t',index_col=0);e=np.load(out/'checkpoints/embeddings.npz')
 assert list(a.obs_names)==list(b.obs_names)==list(pred.index)==list(e['cell_ids'])
 clusters=e['leiden'].astype(str);support,evidence,fine_support,adt_support=expression_evidence(a,b,clusters,out)
 candidate=pred.reference_coarse.to_numpy(dtype=str);confident=(pred.reference_compatible.astype(bool)&(pred.high_probability>=.5)).to_numpy();rna_support=np.asarray([float(support.loc[c,g])if g in support else 0.0 for c,g in zip(clusters,candidate)])
 pct=evidence[evidence.modality=='RNA'].pivot(index='cluster',columns='marker',values='pct_expressing');nonimmune_specific=np.ones(a.n_obs,dtype=bool);nonimmune_rows=[]
 for g in ['Epithelial cells','Fibroblasts','Endothelial cells']:
  genes=[x for x in PANELS[g]if x in pct];n_support=(pct[genes]>=.1).sum(axis=1)if genes else pd.Series(0,index=pct.index);passes=n_support>=2
  ii=np.flatnonzero(candidate==g);nonimmune_specific[ii]=[bool(passes.loc[clusters[j]])for j in ii]
  for c in pct.index:nonimmune_rows.append({'cluster':c,'candidate_lineage':g,'n_panel_markers_detected_in_at_least_10pct_cells':int(n_support.loc[c]),'RNA_set_supported':bool(passes.loc[c]),'panel':';'.join(genes)})
 pd.DataFrame(nonimmune_rows).to_csv(out/'annotation_evidence/rare_nonimmune_RNA_set_evidence.tsv',sep='\t',index=False)
 source_counts=pd.crosstab(pd.Series(clusters,index=a.obs_names,name='leiden'),a.obs.library_id.astype(str));source_fraction=source_counts.max(axis=1)/source_counts.sum(axis=1);source_top=source_counts.idxmax(axis=1);source_counts.to_csv(out/'annotation_evidence/cluster_library_counts.tsv',sep='\t')
 good=confident&(rna_support>=.08)&nonimmune_specific;labels=np.where(good,candidate,'Unassigned')
 confidence=pred.high_probability.to_numpy()*np.minimum(1,rna_support/.3)
 status=np.full(a.n_obs,'evidence_consistent',dtype=object);status[~pred.reference_compatible.astype(bool).to_numpy()]='reference_incompatible';status[pred.reference_compatible.astype(bool).to_numpy()&~confident]='low_reference_probability';status[confident&(rna_support<.08)]='RNA_marker_support_low';status[good&pred.reference_rule.eq('high_ILC_low_NK_compatible').to_numpy()]='ILC_NK_compatible_RNA_supported';status[pred.high_cycling_state.to_numpy()]='high_Cycling_state_lineage_unresolved';status[confident&~nonimmune_specific]='nonimmune_marker_set_insufficient'
 fine=pred.low_label.to_numpy(dtype=str);fine_prob=pred.low_probability.to_numpy();cluster_summary=[];fine_out=np.asarray([g+' (subtype unassigned)'if g!='Unassigned'else'Unassigned'for g in labels],dtype=object)
 unsuitable=['Alveolar','Intestinal','Kidney-resident','Kupffer','Hofbauer','thymocyte','ETP','ELP','pre-B','pro-B','Pro-B','Germinal center','germinal center']
 for c in pd.unique(clusters):
  idx=np.flatnonzero(clusters==c);values=pd.Series(fine[idx]);mode=values.value_counts().index[0];fraction=float(values.eq(mode).mean());coarse_mode=pd.Series(labels[idx]).value_counts().index[0]
  cluster_summary.append({'leiden':c,'n_cells':len(idx),'MajorCellType_mode':coarse_mode,'low_model_mode':mode,'low_model_mode_fraction':fraction,'fraction_evidence_consistent':float(good[idx].mean()),'median_reference_score':float(np.median(confidence[idx])),'source_dominant_library':str(source_top.loc[c]),'source_dominant_fraction':float(source_fraction.loc[c]),'source_dominated_review':bool(source_fraction.loc[c]>=.9)})
  eligible=idx[good[idx]&(fine_prob[idx]>=.7)&(confidence[idx]>=.5)]
  for j in eligible:
   if fraction>=.35 and fine[j]==mode and fine[j] in fine_support and bool(fine_support.loc[c,fine[j]]) and not any(z in fine[j]for z in unsuitable):fine_out[j]=fine[j]
 for col in pred:a.obs['celltypist_'+col]=pred[col].to_numpy()
 a.obs['coarse_reference_confident']=confident;a.obs['annotation_reference_basis']=pd.Categorical(np.where(good,pred.reference_rule.astype(str),'Unassigned'))
 a.obs['leiden']=pd.Categorical(clusters);a.obs['MajorCellType']=pd.Categorical(labels);a.obs['celltype_fine']=pd.Categorical(fine_out);a.obs['candidate_annotation_score']=confidence.astype(np.float32);a.obs['annotation_confidence']=np.where(good,confidence,0).astype(np.float32);a.obs['annotation_uncertainty']=pd.Categorical(status);a.obs['RNA_lineage_marker_support']=rna_support.astype(np.float32)
 adt_auc=np.asarray([float(adt_support.loc[c,g])if g in adt_support else np.nan for c,g in zip(clusters,labels)],dtype=np.float32);a.obs['ADT_cluster_lineage_AUROC']=adt_auc
 a.obs['ADT_crosscheck_status']=pd.Categorical(np.where(~np.isfinite(adt_auc),'unavailable',np.where(adt_auc>=.65,'relative_support',np.where(adt_auc>=.55,'limited_enrichment','weak_enrichment'))))
 a.obs['cluster_source_dominant_library']=pd.Categorical([str(source_top.loc[c])for c in clusters]);a.obs['cluster_source_library_max_fraction']=np.asarray([float(source_fraction.loc[c])for c in clusters],dtype=np.float32);a.obs['source_dominated_cluster_review']=a.obs.cluster_source_library_max_fraction>=.9
 a.obs['fine_label_status']=pd.Categorical(np.where(labels=='Unassigned','lineage_unassigned',np.where(fine_out==fine,'reference_and_marker_supported','subtype_unassigned')))
 a.obs['donor_id']=pd.Categorical(['Unresolved']*a.n_obs);a.obs['clinical_status']='Unresolved';a.obs['clinical_comparison_eligible']=False
 for c in ['Group','source_condition','Disease_detail','Disease_subtype']:a.obs[c]=pd.Categorical(['Unresolved']*a.n_obs)
 a.obsm['X_pca']=e['X_pca'];a.obsm['X_harmony']=e['X_harmony'];a.obsm['X_umap']=e['X_umap'];a.obsp['connectivities']=sp.load_npz(out/'checkpoints/connectivities.npz');a.obsp['distances']=sp.load_npz(out/'checkpoints/distances.npz')
 hv=pd.read_csv(out/'checkpoints/HVG_and_feature_selection.tsv',sep='\t',index_col=0);assert list(hv.index)==list(a.var_names)
 for col in ['highly_variable','protected_lineage_marker','protected_state_marker','protected_biological_marker','integration_feature','highly_variable_nbatches']:
  if col in hv:a.var[col]=hv[col].to_numpy()
 a.uns['scaid']['selected_representation']=str(e['selected_representation'][0])
 a.uns['scaid']['annotation_status']='Two CellTypist references with RNA/ADT evidence cross-check; automated labels, not external expert confirmation'
 a.uns['annotation_methods']={'high_model':'Immune_All_High','low_model':'Immune_All_Low','coarse_reference_rule':'High-model probability>=0.5 and exact coarse agreement, or explicitly marked High-ILC/Low-NK compatibility, plus RNA-marker support. High Cycling state is not independent lineage agreement. Low subtype probability does not veto an otherwise RNA-supported coarse lineage','high_model_coarse_probability_threshold':.5,'RNA_cluster_marker_detection_threshold':.08,'rare_nonimmune_rule':'Epithelial/fibroblast/endothelial calls additionally require at least2predeclared panel genes detected in>=10%of cluster cells; otherwise retain cells as Unassigned','source_dominance_review':'Cluster with>=90%records from one technical library is flagged for review; cells are not removed based on source dominance','fine_probability_threshold':.7,'fine_combined_score_threshold':.5,'fine_cluster_mode_fraction':.35,'fine_marker_rule':'Predeclared subtype panel mean detection>=0.05 or a marker detection>=0.20 within the new Leiden cluster; otherwise retain coarse subtype-unassigned label','confidence_interpretation':'Coarse score=High-model reference probability times bounded RNA marker-support heuristic; low-model subtype probabilities retained separately. Not calibrated clinical probabilities','source_clinical_state':'Unresolved; no AD/HC comparison','ADT_crosscheck':'Within-library, cluster versus rest marker AUROC; source antibodies remain separate from RNA. Relative enrichment flags retained without treating background-limited protein evidence as a hard veto.'}
 pd.DataFrame(cluster_summary).to_csv(out/'annotation_evidence/cluster_annotation_summary.tsv',sep='\t',index=False)
 print('MARKERS_START',flush=True)
 # Balanced, explicitly documented cell sampling for annotation evidence, not case-control inference.
 major=markers(a,labels,out,'DEGs.MajorCellType.txt',1000 if not pilot else 500)
 (out/'MajorCellType_markers.tsv').write_bytes((out/'DEGs.MajorCellType.txt').read_bytes())
 cluster=markers(a,clusters,out,'annotation_evidence/leiden_markers.tsv',200)
 print('MARKERS_READY',len(major),len(cluster),flush=True)
 a.X.sum_duplicates();a.X.sort_indices();a.layers['counts'].sum_duplicates();a.layers['counts'].sort_indices()
 for c in a.obs:
  if a.obs[c].dtype.kind in 'OUS':a.obs[c]=a.obs[c].astype(str).astype('category')
 a.uns=prune_none(dict(a.uns));b.uns=prune_none(dict(b.uns));a.uns['scaid']['H5_compatibility']='Null-valued uns entries omitted; natural logarithm transformation retained'
 b.obs=a.obs.copy();b.obsm['X_umap']=a.obsm['X_umap'].copy();b.write_h5ad(out/'ADT_normalized.h5ad',compression='gzip',compression_opts=1)
 final=out/('pilot_normalized_annotated.h5ad'if pilot else'RNA_normalized_annotated.h5ad');tmp=final.with_suffix('.h5ad.tmp');a.write_h5ad(tmp,compression='gzip',compression_opts=1);os.replace(tmp,final)
 a.obs.to_csv(out/'cell_metadata.tsv.gz',sep='\t',compression='gzip');pd.Series(a.var_names).to_csv(out/'genes.tsv',sep='\t',header=False,index=False)
 make_plots(a,b,evidence,pd.DataFrame(cluster_summary),out)
 # Standalone final readback verifies axes, matrix representation and declared normalization.
 check=ad.read_h5ad(final,backed='r');checkadt=ad.read_h5ad(out/'ADT_normalized.h5ad',backed='r');assert check.obs_names.equals(a.obs_names)and check.var_names.equals(a.var_names)
 assert checkadt.obs_names.equals(check.obs_names)and checkadt.var_names.equals(b.var_names)and check.obs.equals(checkadt.obs)
 assert np.array_equal(check.obsm['X_umap'],checkadt.obsm['X_umap'])and check.layers['counts'].dtype.kind in 'iu'and checkadt.layers['counts'].dtype.kind in 'iu'
 modern={'status':'validated','anndata_version':im.version('anndata'),'RNA_shape':list(check.shape),'ADT_shape':list(checkadt.shape),'axes_equal':True,'metadata_equal':True,'UMAP_equal_finite':True,'counts_integer_dtype':True}
 assert check.obsm['X_umap'].shape==(a.n_obs,2)and np.isfinite(check.obsm['X_umap']).all()
 with h5py.File(final,'r')as f:
  assert f['X'].attrs['encoding-type']=='csr_matrix';assert f['uns/scaid/object_id'].asstr()[:].tolist()==[OBJECT]
 for start in np.unique(np.linspace(0,a.n_obs-1,min(50,a.n_obs),dtype=int)):
  x=check.X[start:start+1].toarray();assert np.isclose(np.expm1(x).sum(),10000,rtol=1e-4,atol=.2)
 check.file.close();checkadt.file.close();write_json(out/'provenance/modern_AnnData_readback.json',modern)
 legacy=validate_legacy_h5ad(final,out/'ADT_normalized.h5ad',out)
 versions={n:im.version(n)for n in ['numpy','scipy','pandas','anndata','scanpy','harmonypy','celltypist','scikit-learn','igraph','leidenalg']}
 ready={'status':'validated','scope':'single-library technical pilot; not a portal release'if pilot else'Full QC-retained mixed AD/HC cohort; clinical mapping unresolved; no case-control analysis','RNA_path':str(final),'RNA_sha256':sha(final),'pipeline_code_sha256':sha(Path(__file__)),'capture_manifest_sha256':sha(ROOT/'capture_manifest.tsv'),'QC_sources':json.loads((out/'checkpoints/assembly.json').read_text())['sources'],'ADT_path':str(out/'ADT_normalized.h5ad'),'ADT_sha256':sha(out/'ADT_normalized.h5ad'),'cell_metadata_path':str(out/'cell_metadata.tsv.gz'),'genes_path':str(out/'genes.tsv'),'marker_table_path':str(out/'MajorCellType_markers.tsv'),'n_cells':a.n_obs,'n_RNA_features':a.n_vars,'n_ADT_features':b.n_vars,'n_libraries':int(a.obs.library_id.nunique()),'cell_axis_sha256':axis_sha(a.obs_names),'feature_axis_sha256':axis_sha(a.var_names),'UMAP_sha256':hashlib.sha256(np.ascontiguousarray(a.obsm['X_umap'],dtype='<f4').data).hexdigest(),'normalized_X_data_sha256':hashlib.sha256(a.X.data.data).hexdigest(),'normalized_X_indices_sha256':hashlib.sha256(a.X.indices.data).hexdigest(),'normalized_X_indptr_sha256':hashlib.sha256(a.X.indptr.data).hexdigest(),'MajorCellType_counts':a.obs.MajorCellType.value_counts().to_dict(),'celltype_fine_counts':a.obs.celltype_fine.value_counts().to_dict(),'MajorCellType_field':'MajorCellType','fine_label_field':'celltype_fine','source_clinical_state':'Unresolved','expert_confirmation_claimed':False,'old_annotations_or_embeddings_reused':False,'software_versions':versions,'legacy_AnnData_readback':legacy,'modern_AnnData_readback':modern,'selected_representation':str(e['selected_representation'][0]),'cell_axis_order':'capture_manifest library order, original retained-cell order within each QC bridge'}
 write_json(out/('pilot_ready.json'if pilot else'analysis_candidate.json'),ready);print('PILOT_READY'if pilot else'ANALYSIS_CANDIDATE',final,a.n_obs,a.n_vars,flush=True)

def make_plots(a,b,evidence,summary,out):
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 from matplotlib.colors import Normalize
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
 xy=a.obsm['X_umap'];rng=np.random.default_rng(SEED);order=rng.permutation(a.n_obs)
 for field,name in [('MajorCellType','UMAP_MajorCellType'),('library_id','UMAP_library'),('pool_id','UMAP_pool'),('leiden','UMAP_leiden')]:
  categories=pd.Categorical(a.obs[field].astype(str));fig,ax=plt.subplots(figsize=(11 if field=='library_id'else 9,6));colors=plt.get_cmap('tab20'if len(categories.categories)<=20 else'turbo',len(categories.categories));ax.scatter(xy[order,0],xy[order,1],c=categories.codes[order],cmap=colors,s=.6,alpha=.65,rasterized=True,linewidths=0);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
  if len(categories.categories)<=20:
   handles=[plt.Line2D([],[],marker='o',linestyle='',color=colors(i),markersize=5,label=c)for i,c in enumerate(categories.categories)];ax.legend(handles=handles,loc='center left',bbox_to_anchor=(1,0.5),frameon=False)
  elif field=='library_id':
   handles=[plt.Line2D([],[],marker='o',linestyle='',color=colors(i),markersize=3,label=c.replace('GSE332623_TotalseqB_',''))for i,c in enumerate(categories.categories)];ax.legend(handles=handles,loc='center left',bbox_to_anchor=(1,0.5),frameon=False,fontsize=6,ncol=2)
  elif field=='leiden':
   for c in categories.categories:
    idx=np.flatnonzero(categories==c);center=np.median(xy[idx],axis=0);ax.text(*center,c,ha='center',va='center',fontsize=6,bbox={'facecolor':'white','alpha':.65,'edgecolor':'none','pad':.6})
  fig.tight_layout();fig.savefig(out/'plots'/(name+'.png'));fig.savefig(out/'plots'/(name+'.pdf'));plt.close(fig)
 # These source-stratified compositions are descriptive retained-record summaries, not donor comparisons.
 composition=pd.crosstab(a.obs.library_id,a.obs.MajorCellType);composition.to_csv(out/'annotation_evidence/MajorCellType_by_library_counts.tsv',sep='\t');fractions=composition.div(composition.sum(axis=1),axis=0);fractions.to_csv(out/'annotation_evidence/MajorCellType_by_library_fractions.tsv',sep='\t')
 fig,ax=plt.subplots(figsize=(11,max(5,len(fractions)*.23)));fractions.rename(index=lambda c:str(c).replace('GSE332623_TotalseqB_','')).plot.barh(stacked=True,ax=ax,width=.85,legend=False,colormap='tab20');ax.set_xlabel('Fraction of QC-retained library-barcode records');ax.set_ylabel('Technical library');ax.legend(loc='center left',bbox_to_anchor=(1,.5),frameon=False,fontsize=7);fig.tight_layout();fig.savefig(out/'plots/library_celltype_composition.png');fig.savefig(out/'plots/library_celltype_composition.pdf');plt.close(fig)
 subset=np.sort(rng.choice(a.n_obs,min(12000,a.n_obs),replace=False));fig,axs=plt.subplots(2,2,figsize=(11,8));allcoarse=pd.Categorical(a.obs.MajorCellType.astype(str));allpool=pd.Categorical(a.obs.pool_id.astype(str))
 for row,rep in enumerate(['X_pca','X_harmony']):
  for col,(field,cat) in enumerate([('MajorCellType',allcoarse),('pool_id',allpool)]):
   ax=axs[row,col];z=np.asarray(a.obsm[rep][subset,:2]);cmap=plt.get_cmap('tab20',len(cat.categories));ax.scatter(z[:,0],z[:,1],c=cat.codes[subset],cmap=cmap,s=1.2,alpha=.6,linewidths=0,rasterized=True);ax.set_xlabel(('PC'if row==0 else'Harmony')+'1');ax.set_ylabel(('PC'if row==0 else'Harmony')+'2');ax.set_title(('PCA'if row==0 else'Harmony')+' — '+('RNA-supported lineage'if col==0 else'Source pool'))
   if row==0:ax.legend(handles=[plt.Line2D([],[],marker='o',linestyle='',color=cmap(i),markersize=3,label=c)for i,c in enumerate(cat.categories)],loc='lower center',bbox_to_anchor=(.5,1.18),frameon=False,fontsize=6,ncol=3)
 fig.tight_layout();fig.savefig(out/'plots/PCA_Harmony_reference_pool.png');fig.savefig(out/'plots/PCA_Harmony_reference_pool.pdf');plt.close(fig)
 # Marker dot plot uses exact means/detection from all retained cells, not the marker-test sample.
 r=evidence[evidence.modality=='RNA'];genes=[g for g in ['CD3D','TRAC','NKG7','GNLY','MS4A1','CD79A','MZB1','JCHAIN','LYZ','LST1','CD1C','LILRA4','IL3RA','PF4','PPBP','HBB']if g in set(r.marker)];groups=a.obs.MajorCellType.astype(str);levels=pd.Index(pd.unique(groups));rows=[]
 for g in levels:
  idx=np.flatnonzero(groups.to_numpy()==g);x=a.X[idx][:,a.var_names.get_indexer(genes)];rows.append((np.asarray(x.mean(axis=0)).ravel(),np.asarray((x>0).mean(axis=0)).ravel()))
 means=np.vstack([x[0]for x in rows]);pct=np.vstack([x[1]for x in rows]);scaled=(means-means.mean(axis=0))/(means.std(axis=0)+1e-6);fig,ax=plt.subplots(figsize=(11,max(4,len(levels)*.35)));xx,yy=np.meshgrid(np.arange(len(genes)),np.arange(len(levels)));s=ax.scatter(xx,yy,s=15+200*pct,c=scaled,cmap='RdBu_r',vmin=-2,vmax=2,linewidth=.25,edgecolors='grey');ax.set_xticks(np.arange(len(genes)),genes,rotation=65,ha='right');ax.set_yticks(np.arange(len(levels)),levels);ax.invert_yaxis();fig.colorbar(s,ax=ax,label='Gene-wise standardized mean log expression',shrink=.65);fig.tight_layout();fig.savefig(out/'plots/RNA_marker_dotplot.png');fig.savefig(out/'plots/RNA_marker_dotplot.pdf');plt.close(fig)
 if len(evidence[evidence.modality=='ADT']):
  t=evidence[evidence.modality=='ADT'].pivot(index='cluster',columns='marker',values='mean_log_normalized');z=(t-t.mean())/(t.std()+1e-6);fig,ax=plt.subplots(figsize=(11,max(4,len(t)*.16)));imx=ax.imshow(z.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2);ax.set_xticks(range(len(t.columns)),t.columns,rotation=70,ha='right',fontsize=6);ax.set_yticks(range(len(t)),t.index,fontsize=6);ax.set_ylabel('Leiden cluster');fig.colorbar(imx,ax=ax,label='Antibody-wise standardized mean CLR',shrink=.6);fig.tight_layout();fig.savefig(out/'plots/ADT_marker_crosscheck.png');fig.savefig(out/'plots/ADT_marker_crosscheck.pdf');plt.close(fig)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['assemble','classify','embed','finalize']);p.add_argument('--out',type=Path,default=ROOT/'integration');p.add_argument('--bridges',nargs='*',type=Path);p.add_argument('--pilot',action='store_true');args=p.parse_args();setup(args.out)
 if args.stage=='assemble':assemble(args.out,args.bridges)
 elif args.stage=='classify':classify(args.out)
 elif args.stage=='embed':embed(args.out)
 elif args.stage=='finalize':finalize(args.out,args.pilot)
if __name__=='__main__':main()

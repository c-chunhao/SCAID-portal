"""Full-library, sparse, frozen-model T1D reannotation pilot (no mass relabel).

All cells are classified independently in bounded 5000-cell blocks. There is no
majority vote, no clustering inferred from candidate labels, and no integrated
or z-scored expression as model input. Legacy/source labels are comparison
columns, never classifier inputs. Outputs are candidate labels for review.
"""
from pathlib import Path
import hashlib, json, time, importlib.metadata
import numpy as np, pandas as pd, h5py
from scipy import sparse
import anndata as ad
import celltypist
from celltypist import models
from released_rna import validate_raw_csr_arrays

OUT=Path(__file__).resolve().parent
H5=None
SOURCE=None
SEED=20261005
MODEL_DIR=Path.home()/'.celltypist/data/models'

MARKERS={
 'T/NK cells':['CD3D','CD3E','TRAC','NKG7','GNLY','KLRD1'],
 'B/Plasma cells':['CD79A','CD79B','MS4A1','CD37','JCHAIN','MZB1'],
 'Mononuclear Phagocytes':['FCN1','S100A8','S100A9','LST1','CTSS','FCGR3A'],
 'DC':['CD1C','FCER1A','CLEC10A','CLEC9A','XCR1'],
 'pDC':['TCF4','IL3RA','CLEC4C','GZMB'],
}
SOURCE_MAP={'CD4_T':'T/NK cells','CD8_T':'T/NK cells','NK':'T/NK cells','MAIT':'T/NK cells','T_reg':'T/NK cells','VD2p':'T/NK cells','B_Naive':'B/Plasma cells','B_SM':'B/Plasma cells','B_Plasma':'B/Plasma cells','C_Monocyte':'Mononuclear Phagocytes','NC_Monocyte':'Mononuclear Phagocytes','I_Monocyte':'Mononuclear Phagocytes','cDCs':'DC','pDCs':'DC'}
FINAL_MAP={'T cells':'T/NK cells','NK cells':'T/NK cells','NKT cells':'T/NK cells','ILC':'T/NK cells','B cells':'B/Plasma cells','Plasma cells':'B/Plasma cells','DCs':'DC','pDCs':'DC','pDC':'DC'}

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while True:
   b=f.read(1024*1024)
   if not b:break
   h.update(b)
 return h.hexdigest()

def hcol(h,c):
 q=h['obs'][c]
 if isinstance(q,h5py.Group):
  codes=q['codes'][:];cats=q['categories'].asstr()[:]
  return np.where(codes<0,'',cats[np.maximum(codes,0)])
 return q.asstr()[:] if q.dtype.kind in ['S','O','U'] else q[:]

def broad(label):
 # Frozen mapping of model vocabulary. Progenitors, cycling-without-lineage,
 # and nonimmune cells remain unknown; no posthoc adaptation to source labels.
 if label in ['DC','DC precursor','DC1','DC2','DC3','Migratory DCs','Cycling DCs','Transitional DC','pDC','pDC precursor']:return 'DC'
 if label in ['MNP','Mono-mac','Monocytes','Monocyte precursor','Classical monocytes','Non-classical monocytes','Cycling monocytes','Hofbauer cells','Kupffer cells'] or 'macrophage' in label.lower():return 'Mononuclear Phagocytes'
 if label in ['B-cell lineage','Plasma cells','Plasmablasts','B cells'] or ' B cells' in label:return 'B/Plasma cells'
 if label in ['ILC','ILC precursor','ILC1','ILC2','ILC3','T cells','NK cells','NKT cells','CD16+ NK cells','CD16- NK cells','CD8a/a','CD8a/b(entry)','Transitional NK','Cycling NK cells','MAIT cells','T(agonist)','Double-negative thymocytes','Double-positive thymocytes'] or ' T cells' in label or label.startswith(('Tcm/','Tem/','Trm ')) or label=='Treg(diff)':return 'T/NK cells'
 return 'Other / Unresolved'

def run():
 start=time.time();np.random.seed(SEED)
 OUT.mkdir(parents=True,exist_ok=True)
 cached={n:MODEL_DIR/(n+'.pkl') for n in ['Immune_All_High','Immune_All_Low']}
 model={n:models.Model.load(str(p)) for n,p in cached.items()}
 results=[];norm_audit=[]
 with h5py.File(H5) as h:
  ids=h['names_obs'].asstr()[:];genes=h['names_var'].asstr()[:]
  assert len(set(ids))==len(ids) and len(set(genes))==len(genes)
  legacy=pd.Series(hcol(h,'MajorCellType')).replace(FINAL_MAP).to_numpy();donor=hcol(h,'Sample')
  source=pd.read_csv(SOURCE,sep='\t',dtype=str,keep_default_na=False).set_index('cell_id').reindex(ids)
  assert source.Sample_ID.notna().all() and np.array_equal(source.Sample_ID.to_numpy(),donor)
  gene_index={g:i for i,g in enumerate(genes)}
  marker_idx={n:[gene_index[g] for g in gs if g in gene_index] for n,gs in MARKERS.items()}
  blocks=h['assay/RNA/layers/rawdata'];offset=0
  for k in sorted(blocks):
   g=blocks[k];shape=tuple(int(v) for v in g.attrs['shape'])
   assert shape[1]==len(genes)
   counts=validate_raw_csr_arrays(g['data'][:],g['indices'][:],g['indptr'][:],shape).astype(np.float32)
   total=np.asarray(counts.sum(axis=1)).ravel();assert (total>0).all()
   b_ids=ids[offset:offset+shape[0]]
   df=pd.DataFrame({'cell_id':b_ids,'donor':donor[offset:offset+shape[0]],'legacy_curated_coarse':legacy[offset:offset+shape[0]],'source_author_label':source.Cluster_Annotation_Merged.iloc[offset:offset+shape[0]].to_numpy()})
   df['source_author_coarse']=df.source_author_label.map(SOURCE_MAP).fillna('Other / Unresolved')
   for n,idx in marker_idx.items():
    df[n+'_marker_detected_genes']=np.asarray((counts[:,idx]>0).sum(axis=1)).ravel()
   # Normalize by ALL retained raw genes before any model/HVG restriction.
   x=counts.multiply((10000/total)[:,None]).tocsr();x.data=np.log1p(x.data)
   sample_sum=np.asarray(np.expm1(x[:min(12,shape[0])].toarray()).sum(axis=1)).ravel()
   assert np.max(np.abs(sample_sum-10000))<.05
   norm_audit.append({'block':k,'cells':shape[0],'full_features':shape[1],'raw_nnz':counts.nnz,'raw_library_min':float(total.min()),'raw_library_max':float(total.max()),'normalized_full_library_sum_min':float(sample_sum.min()),'normalized_full_library_sum_max':float(sample_sum.max())})
   a=ad.AnnData(X=x,obs=pd.DataFrame(index=b_ids),var=pd.DataFrame(index=genes))
   for n,m in model.items():
    pred=celltypist.annotate(a,model=m,majority_voting=False)
    labels=pred.predicted_labels.predicted_labels.astype(str).to_numpy()
    prob=pred.probability_matrix.to_numpy();top=np.sort(prob,axis=1)[:,-2:]
    df[n+'_label']=labels;df[n+'_coarse']=[broad(v) for v in labels]
    df[n+'_max_sigmoid_score']=top[:,-1];df[n+'_top_score_margin']=top[:,-1]-top[:,-2]
    del pred,prob,top
   agreed=df.Immune_All_High_coarse.eq(df.Immune_All_Low_coarse)&df.Immune_All_High_coarse.ne('Other / Unresolved')
   # Agreement supplies a review candidate, not calibrated correctness.
   df['vNext_lineage_candidate']=np.where(agreed,df.Immune_All_High_coarse,'Unresolved')
   df['candidate_matches_legacy']=df.vNext_lineage_candidate.eq(df.legacy_curated_coarse)
   df['candidate_matches_source']=df.vNext_lineage_candidate.eq(df.source_author_coarse)
   adequate=(df.Immune_All_High_max_sigmoid_score.ge(.5)&df.Immune_All_Low_max_sigmoid_score.ge(.5)&df.Immune_All_High_top_score_margin.ge(.1)&df.Immune_All_Low_top_score_margin.ge(.1))
   df['review_status']=np.select([~agreed,agreed&~adequate,agreed&adequate&~df.candidate_matches_legacy],[ 'model_coarse_disagreement_or_unknown','model_agreement_low_margin_requires_review','model_agreement_changes_legacy_requires_marker_review'],default='model_agreement_supports_legacy')
   # DC and mono overlapping markers are displayed as evidence, never used
   # to force all source-labelled DCs or all misclassified cells to one label.
   dc=(df['DC_marker_detected_genes']>=2)|(df['pDC_marker_detected_genes']>=2)
   mono=df['Mononuclear Phagocytes_marker_detected_genes']>=3
   df['DC_Mono_marker_review']=np.select([dc&mono,dc&~mono,~dc&mono],['mixed_marker_evidence','DC_marker_support','monocyte_marker_support'],default='sparse_or_no_specific_marker_support')
   results.append(df);offset+=shape[0]
   print(f'Finished {k}: {offset}/{len(ids)} cells',flush=True)
   del counts,x,a
  assert offset==len(ids)
  gene_coverage={n:{'model_features':len(m.features),'full_raw_features_present':len(set(m.features)&set(genes)),'missing_model_features':sorted(set(m.features)-set(genes))} for n,m in model.items()}
 d=pd.concat(results,ignore_index=True);assert d.cell_id.is_unique and len(d)==len(ids)
 d.to_csv(OUT/'T1D_full_RNA_CellTypist_vNext_candidate_labels.tsv.gz',sep='\t',index=False,compression='gzip')
 d.groupby(['source_author_coarse','legacy_curated_coarse','vNext_lineage_candidate','review_status'],observed=True).size().reset_index(name='cells').to_csv(OUT/'T1D_full_RNA_candidate_decision_counts.tsv',sep='\t',index=False)
 d.loc[d.source_author_coarse.eq('DC')|d.legacy_curated_coarse.eq('DC')|d.vNext_lineage_candidate.eq('DC')].to_csv(OUT/'T1D_DC_Mono_candidate_review.tsv.gz',sep='\t',index=False,compression='gzip')
 per=d.groupby('donor').agg(cells=('cell_id','size'),candidate_agreed_cells=('vNext_lineage_candidate',lambda v:v.ne('Unresolved').sum()),legacy_concordant=('candidate_matches_legacy','sum'),source_concordant=('candidate_matches_source','sum')).reset_index()
 per.to_csv(OUT/'T1D_full_RNA_candidate_per_donor.tsv',sep='\t',index=False)
 cm=pd.crosstab(d.legacy_curated_coarse,d.vNext_lineage_candidate);cm.to_csv(OUT/'T1D_legacy_vs_full_RNA_candidate_counts.tsv',sep='\t')
 comparable=d.source_author_coarse.ne('Other / Unresolved')
 summary={'object':'T1DM_syn53641849_PBMC','cells':len(d),'full_RNA_features':len(genes),'source_donors':d.donor.nunique(),'full_library_normalization':'Raw CSR counts / all-gene cell library * 10000, then log1p; no HVG normalization or integrated expression.','prediction':'Frozen cached Immune_All_High and Immune_All_Low models; majority_voting=False, no source labels used as predictors.','coarse_vocabulary_mapping':{n:{v:broad(v) for v in m.cell_types} for n,m in model.items()},'model_feature_coverage':gene_coverage,'model_coarse_agreement_cells':int(d.vNext_lineage_candidate.ne('Unresolved').sum()),'unresolved_cells':int(d.vNext_lineage_candidate.eq('Unresolved').sum()),'candidate_legacy_changes_among_agreed':int((~d.candidate_matches_legacy&d.vNext_lineage_candidate.ne('Unresolved')).sum()),'review_status_counts':d.review_status.value_counts().to_dict(),'source_DC_candidate_counts':d.loc[d.source_author_coarse.eq('DC'),'vNext_lineage_candidate'].value_counts().to_dict(),'old_source_DC_curated_counts':d.loc[d.source_author_coarse.eq('DC'),'legacy_curated_coarse'].value_counts().to_dict(),'source_comparable_cells':int(comparable.sum()),'cell_weighted_source_concordance_of_legacy':float(d.loc[comparable].legacy_curated_coarse.eq(d.loc[comparable].source_author_coarse).mean()),'cell_weighted_source_concordance_of_candidate_including_unresolved':float(d.loc[comparable].candidate_matches_source.mean()),'limits':['Source-author agreement is not blind annotation accuracy.','Model agreement/probabilities are not independent or calibrated correctness.','No fine-subtype label is accepted merely from model prediction.','No production annotations are replaced; candidate changes require marker review and donor/source context.','Pilot is one T1D PBMC object, not a whole-atlas reannotation.'],'normalization_audit':norm_audit,'seed':SEED,'wall_seconds':time.time()-start,'source_sha256':{str(H5):sha(H5),str(SOURCE):sha(SOURCE),**{str(p):sha(p) for p in cached.values()}},'versions':{n:importlib.metadata.version(n) for n in ['celltypist','anndata','scanpy','numpy','scipy','scikit-learn']}}
 (OUT/'T1D_full_RNA_reannotation_pilot_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
 (OUT/'T1D_frozen_model_annotation_policy.json').write_text(json.dumps({'input':'full-RNA log1p(CPM10000)','coarse_mapping':summary['coarse_vocabulary_mapping'],'candidate_acceptance':'High and Low agree on known coarse lineage; all disagreements are Unresolved. This does not accept a production relabel.','review_flag_sigmoid_threshold':.5,'review_flag_top_score_margin':.1,'probability_warning':'Independent sigmoid classifier scores; not calibrated cell-type posterior probabilities.','source_labels_not_inputs':True,'seed':SEED},indent=2)+'\n')
 print(json.dumps({k:v for k,v in summary.items() if k not in ['coarse_vocabulary_mapping','model_feature_coverage','normalization_audit']},indent=2),flush=True)

if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser()
 parser.add_argument('--h5',type=Path,required=True)
 parser.add_argument('--source-metadata',type=Path,required=True)
 parser.add_argument('--models',type=Path,required=True)
 parser.add_argument('--out',type=Path,required=True)
 args=parser.parse_args();H5=args.h5;SOURCE=args.source_metadata;MODEL_DIR=args.models;OUT=args.out
 run()

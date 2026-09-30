import os
"""Second source-author broad-lineage check: official T1D repository RDS labels.
Read-only source metadata extraction; exact barcode/source-individual joins.
"""
from pathlib import Path
import json,hashlib
import h5py,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix,cohen_kappa_score,precision_recall_fscore_support
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures'
SOURCE=T/'T1D_original_author_metadata.tsv.gz';H5=Path(os.environ['SCAID_T1D_RELEASE_H5'])
ONTOLOGY={'CD4_T':'T/NK cells','CD8_T':'T/NK cells','NK':'T/NK cells','MAIT':'T/NK cells','T_reg':'T/NK cells','VD2p':'T/NK cells','B_Naive':'B/Plasma cells','B_SM':'B/Plasma cells','B_Plasma':'B/Plasma cells','C_Monocyte':'Mononuclear Phagocytes','NC_Monocyte':'Mononuclear Phagocytes','I_Monocyte':'Mononuclear Phagocytes','cDCs':'DC','pDCs':'DC'}
FINAL_ALIASES={'T cells':'T/NK cells','NK cells':'T/NK cells','NKT cells':'T/NK cells','ILC':'T/NK cells','B cells':'B/Plasma cells','Plasma cells':'B/Plasma cells','DCs':'DC','pDCs':'DC','pDC':'DC'}
CLASSES=['T/NK cells','B/Plasma cells','Mononuclear Phagocytes','DC']
(T/'T1D_annotation_frozen_ontology.json').write_text(json.dumps({'source_Cluster_Annotation_Merged_to_curated_coarse':ONTOLOGY,'curated_label_vocabulary_aliases':FINAL_ALIASES,'excluded_categories':['Undetermined','HSCs','Platelet','Erythrocytes'],'source_undetermined_policy':'No source reference lineage assigned; excluded from comparable coverage, not promoted to correct predictions.','unknown_curated_labels':'Retained as disagreements for all comparable reference cells'},indent=2)+'\n')
def col(h,k):
 n=h['obs'][k]
 if isinstance(n,h5py.Group):
  cats=n['categories'].asstr()[:];codes=n['codes'][:];return np.where(codes<0,'',cats[np.maximum(codes,0)])
 return n.asstr()[:]
with h5py.File(H5)as h:ids=h['names_obs'].asstr()[:];final=col(h,'MajorCellType');donor=col(h,'Sample')
m=pd.read_csv(SOURCE,sep='\t',keep_default_na=False).set_index('cell_id');assert m.index.is_unique
m=m.reindex(ids);assert m.Sample_ID.notna().all();assert np.array_equal(m.Sample_ID.astype(str).to_numpy(),donor);assert m.COND.eq('T1D').all()
final=pd.Series(final).replace(FINAL_ALIASES).to_numpy();ref=m.Cluster_Annotation_Merged.map(ONTOLOGY);keep=ref.notna().to_numpy();a=ref[keep].to_numpy();p=np.where(np.isin(final[keep],CLASSES),final[keep],'Other / Unknown');d=donor[keep]
labels=CLASSES+['Other / Unknown'];cm=confusion_matrix(a,p,labels=labels);pd.DataFrame(cm,index=labels,columns=labels).to_csv(T/'T1D_source_author_vs_curated_coarse_confusion_counts.tsv',sep='\t')
pr,re,f1,sup=precision_recall_fscore_support(a,p,labels=CLASSES,zero_division=0);pd.DataFrame({'lineage':CLASSES,'precision':pr,'recall':re,'F1':f1,'source_cells':sup}).to_csv(T/'T1D_source_author_concordance_per_lineage.tsv',sep='\t',index=False)
per=pd.DataFrame({'donor':d,'agree':a==p}).groupby('donor').agg(n_cells=('agree','size'),agree_cells=('agree','sum'),agreement=('agree','mean')).reset_index();per.to_csv(T/'T1D_source_author_concordance_per_donor.tsv',sep='\t',index=False)
v=per.agreement.to_numpy();rng=np.random.default_rng(20260930);boots=np.mean(v[rng.integers(0,len(v),size=(10000,len(v)))],axis=1)
s={'object':'T1DM_syn53641849_PBMC','source_author_object':os.environ['SCAID_T1D_SOURCE_RDS'],'source_metadata_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'published_cells':len(ids),'exact_source_barcodes':len(ids),'exact_source_donor_matches':len(ids),'source_disease_consistent':'All matched source COND=T1D; actual individual IDs are source Sample_ID.','reference_comparable_lineage_cells':int(keep.sum()),'reference_excluded_by_category':m.loc[~keep,'Cluster_Annotation_Merged'].value_counts().to_dict(),'coverage':float(keep.mean()),'agreement_cell_weighted':float((a==p).mean()),'agreement_donor_unweighted_mean':float(v.mean()),'donor_cluster_bootstrap_mean_95CI':np.quantile(boots,[.025,.975]).tolist(),'donors':len(per),'cohen_kappa':float(cohen_kappa_score(a,p,labels=labels)),'macro_F1_four_lineages':float(f1.mean()),'curated_other_unknown_among_comparable':int((p=='Other / Unknown').sum()),'cannot_claim':['Blinded gold-standard annotation accuracy; author annotations could have informed curation.','Fine-subtype validity','Whole-atlas accuracy.','Separate cell-type definitions outside frozen comparable lineages.'],'reference':'Original source-author Cluster_Annotation_Merged in official repository Seurat object; no rerun of the stored CellTypist classifier.'}
(B/'T1D_annotation_reference_summary.json').write_text(json.dumps(s,indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(6.5,4.6));norm=cm[:4,:]/cm[:4,:].sum(axis=1,keepdims=True);im=ax.imshow(norm,aspect='auto',cmap='Blues',vmin=0,vmax=1);short=['T/NK','B/plasma','Mono. phagocytes','DC','Other / Unknown'];ax.set_xticks(range(5),short,rotation=25,ha='right');ax.set_yticks(range(4),short[:4]);ax.set_xlabel('Final curated lineage');ax.set_ylabel('Source-author lineage');ax.set_title('Source-reference agreement: T1D syn53641849')
for i in range(4):
 for j in range(5):ax.text(j,i,f'{norm[i,j]:.1%}',ha='center',va='center',fontsize=8,color='white'if norm[i,j]>.5 else'black')
fig.colorbar(im,ax=ax,label='Fraction of source-author lineage',fraction=.035,pad=.04);fig.text(.5,.02,f'Exact barcode/individual joins; {len(per)} source donors; broad lineage comparison.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.045,1,1]);fig.savefig(F/'Annotation_source_reference_T1D.pdf');fig.savefig(F/'Annotation_source_reference_T1D.png',dpi=600);fig.savefig(F/'Annotation_source_reference_T1D.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig);print(json.dumps(s,indent=2))

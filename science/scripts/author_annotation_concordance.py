import os
"""Compare source-author cg_cov with final curated MajorCellType.
The source barcode/donor check is independent of the stored CellTypist rerun.
Agreement is not blinded annotation accuracy and is not generalized to 48 objects.
"""
from pathlib import Path
import hashlib,json
import h5py,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix,cohen_kappa_score,precision_recall_fscore_support
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures'
SOURCE=Path(os.environ['SCAID_SLE_SOURCE_METADATA'])
H5=Path(os.environ['SCAID_SLE_RELEASE_H5'])
# Frozen ontology: compare broad immune lineages. Proliferation is a cell state;
# progenitor coverage is tiny and does not map to the same granularity.
ONTOLOGY={'T4':'T/NK cells','T8':'T/NK cells','NK':'T/NK cells','B':'B/Plasma cells','PB':'B/Plasma cells','cM':'Mononuclear Phagocytes','ncM':'Mononuclear Phagocytes','cDC':'DC','pDC':'DC'}
CLASSES=['T/NK cells','B/Plasma cells','Mononuclear Phagocytes','DC']
(T/'annotation_frozen_ontology.json').write_text(json.dumps({'source_author_cg_cov_to_curated_coarse':ONTOLOGY,'excluded_reference_categories':{'Prolif':'state not lineage','Progen':'granularity differs'},'unknown_curated_labels':'included as disagreements; not silently dropped'},indent=2)+'\n')
def col(h,k):
 n=h['obs'][k]
 return n.asstr()[:] if not isinstance(n,h5py.Group) else n['categories'].asstr()[:][n['codes'][:]]
with h5py.File(H5)as h:
 ids=h['names_obs'].asstr()[:];final=col(h,'MajorCellType');donor=col(h,'Sample');savedcg=col(h,'cg_cov')
meta=pd.read_csv(SOURCE,usecols=['Unnamed: 0','ind_cov','cg_cov','SLE_status']);meta=meta.rename(columns={meta.columns[0]:'cell_id'}).set_index('cell_id')
assert meta.index.is_unique
m=meta.reindex(ids);assert m.ind_cov.notna().all();assert np.array_equal(m.ind_cov.astype(str).to_numpy(),donor);assert np.array_equal(m.cg_cov.astype(str).to_numpy(),savedcg)
ref=m.cg_cov.map(ONTOLOGY);keep=ref.notna().to_numpy();a=ref.iloc[np.flatnonzero(keep)].to_numpy();p=np.where(np.isin(final[keep],CLASSES),final[keep],'Other / Unknown');d=donor[keep]
labels=CLASSES+['Other / Unknown'];cm=confusion_matrix(a,p,labels=labels)
pd.DataFrame(cm,index=labels,columns=labels).to_csv(T/'SLE_source_author_vs_curated_coarse_confusion_counts.tsv',sep='\t')
pr,re,f1,sup=precision_recall_fscore_support(a,p,labels=CLASSES,zero_division=0)
per=pd.DataFrame({'lineage':CLASSES,'precision':pr,'recall':re,'F1':f1,'source_cells':sup});per.to_csv(T/'SLE_source_author_concordance_per_lineage.tsv',sep='\t',index=False)
perdon=pd.DataFrame({'donor':d,'agree':a==p,'reference':a,'curated':p}).groupby('donor').agg(n_cells=('agree','size'),agree_cells=('agree','sum'),agreement=('agree','mean')).reset_index()
perdon.to_csv(T/'SLE_source_author_concordance_per_donor.tsv',sep='\t',index=False)
# Clustered bootstrap treats donors as units, not individual cells.
rng=np.random.default_rng(20260930);v=perdon.agreement.to_numpy();boots=np.mean(v[rng.integers(0,len(v),size=(10000,len(v)))],axis=1)
s={'object':'SLE_GSE174188_PBMC','source_metadata_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'published_cells':len(ids),'exact_source_barcodes':len(ids),'exact_source_donor_matches':len(ids),'exact_saved_cg_cov_matches':len(ids),'reference_comparable_lineage_cells':int(keep.sum()),'reference_excluded_by_category':m.loc[~keep,'cg_cov'].value_counts().to_dict(),'coverage':float(keep.mean()),'agreement_cell_weighted':float((a==p).mean()),'agreement_donor_unweighted_mean':float(v.mean()),'donor_cluster_bootstrap_mean_95CI':np.quantile(boots,[.025,.975]).tolist(),'donors':len(perdon),'cohen_kappa':float(cohen_kappa_score(a,p,labels=labels)),'macro_F1_four_lineages':float(f1.mean()),'curated_other_unknown_among_comparable':int((p=='Other / Unknown').sum()),'cannot_claim':['Blinded gold-standard accuracy','Fine subtype validity','Whole-atlas 48-object benchmark','Independence of reference from possible manual curation; original-author labels may have been available to curators.'],'reference':'Original study cg_cov from source author metadata, exactly barcode-joined, not a repeat of the stored CellTypist classifier.'}
(B/'annotation_reference_summary.json').write_text(json.dumps(s,indent=2)+'\n')
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42})
fig,ax=plt.subplots(figsize=(6.5,4.6));norm=np.divide(cm[:4,:],cm[:4,:].sum(axis=1,keepdims=True),where=cm[:4,:].sum(axis=1,keepdims=True)>0)
im=ax.imshow(norm,cmap='Blues',vmin=0,vmax=1,aspect='auto');short=['T/NK','B/plasma','Mono. phagocytes','DC','Other / Unknown'];ax.set_xticks(range(5),short,rotation=25,ha='right');ax.set_yticks(range(4),short[:4]);ax.set_ylabel('Source-author lineage');ax.set_xlabel('Final curated lineage');ax.set_title('Source-reference agreement: GSE174188')
for i in range(4):
 for j in range(5):ax.text(j,i,f'{norm[i,j]:.1%}',ha='center',va='center',fontsize=8,color='white'if norm[i,j]>.5 else'black')
fig.colorbar(im,ax=ax,label='Fraction of source-author lineage',fraction=.035,pad=.04)
fig.text(.5,.02,f'Exact barcode joins; {len(perdon)} source donors; broad lineage comparison',ha='center',fontsize=8)
fig.tight_layout(rect=[0,.045,1,1]);fig.savefig(F/'Annotation_source_reference.pdf');fig.savefig(F/'Annotation_source_reference.png',dpi=600);fig.savefig(F/'Annotation_source_reference.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
print(json.dumps(s,indent=2))

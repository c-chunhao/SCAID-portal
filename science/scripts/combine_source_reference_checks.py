import os
from pathlib import Path
import json
import h5py,numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures'
ALIASES={'T cells':'T/NK cells','NK cells':'T/NK cells','NKT cells':'T/NK cells','ILC':'T/NK cells','B cells':'B/Plasma cells','Plasma cells':'B/Plasma cells','DCs':'DC','pDCs':'DC','pDC':'DC'}
def col(h,k):
 n=h['obs'][k]
 if isinstance(n,h5py.Group):
  cats=n['categories'].asstr()[:];co=n['codes'][:];return np.where(co<0,'',cats[np.maximum(co,0)])
 return n.asstr()[:]
with h5py.File(os.environ['SCAID_T1D_RELEASE_H5'])as h:ids=h['names_obs'].asstr()[:];raw=col(h,'MajorCellType');donor=col(h,'Sample')
m=pd.read_csv(T/'T1D_original_author_metadata.tsv.gz',sep='\t',keep_default_na=False).set_index('cell_id').reindex(ids)
final=pd.Series(raw).replace(ALIASES).to_numpy();srcmap=json.loads((T/'T1D_annotation_frozen_ontology.json').read_text())['source_Cluster_Annotation_Merged_to_curated_coarse'];ref=m.Cluster_Annotation_Merged.map(srcmap).fillna('Not comparable').to_numpy()
pairs=pd.DataFrame({'source_author_label':m.Cluster_Annotation_Merged.to_numpy(),'curated_raw_label':raw,'source_coarse_reference':ref,'curated_coarse_alias':final}).groupby(['source_author_label','curated_raw_label','source_coarse_reference','curated_coarse_alias']).size().rename('cells').reset_index();pairs.to_csv(T/'T1D_source_author_raw_and_canonical_label_pairs.tsv',sep='\t',index=False)
bad=(ref=='DC')&(final=='Mononuclear Phagocytes');pd.DataFrame({'source_donor':donor[bad],'source_author_label':m.Cluster_Annotation_Merged.to_numpy()[bad],'curated_label':raw[bad]}).groupby(['source_donor','source_author_label','curated_label']).size().rename('cells').reset_index().to_csv(T/'T1D_source_DC_curated_monocyte_disagreement_by_donor.tsv',sep='\t',index=False)
(B/'T1D_ontology_QA_correction.json').write_text(json.dumps({'initial_invalid_string_matching_result':.2626486703193906,'cause':'The initial comparator classified final T cells/NK cells/DCs/pDCs as Other because its source ontology used combined T/NK and DC broad categories. This was an incompatible-vocabulary implementation error, not a changed annotation outcome.','correction':'Explicit semantic curated vocabulary aliases standardized before canonical lineage comparison. No source labels, final labels, cells or models were reannotated.','final_vocabulary_aliases':ALIASES,'raw_and_canonical_pairs':str(T/'T1D_source_author_raw_and_canonical_label_pairs.tsv'),'remaining_true_DC_discordance':{'source_DC_cells':int((ref=='DC').sum()),'classified_final_Mononuclear_Phagocytes':int(bad.sum()),'fraction':float(bad.sum()/(ref=='DC').sum())},'no_accuracy_claim':'Non-blinded original-author reference concordance; source-author labels could have informed final curation.'},indent=2)+'\n')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(6.5,4.1));short=['T/NK','B/plasma','Mono','DC','Other']
for ax,study,prefix,letter in zip(axes,['SLE GSE174188','T1D syn53641849'],['SLE','T1D'],'AB'):
 cm=pd.read_csv(T/(prefix+'_source_author_vs_curated_coarse_confusion_counts.tsv'),sep='\t',index_col=0).to_numpy()[:4,:];v=cm/cm.sum(axis=1,keepdims=True);im=ax.imshow(v,cmap='Blues',vmin=0,vmax=1,aspect='auto');ax.set_yticks(range(4),short[:4]);ax.set_xticks(range(5),short,rotation=40,ha='right');ax.set_title(letter+'  '+study,loc='left',fontsize=9);ax.set_xlabel('Final curated broad lineage');
 for i in range(4):
  for j in range(5):ax.text(j,i,f'{v[i,j]*100:.1f}',ha='center',va='center',fontsize=8,color='white'if v[i,j]>.5 else'black')
axes[0].set_ylabel('Source-author broad lineage');fig.subplots_adjust(left=.105,right=.975,bottom=.24,top=.84,wspace=.5);cax=fig.add_axes([.30,.91,.4,.025]);cb=fig.colorbar(im,cax=cax,orientation='horizontal');cb.set_ticks([0,.5,1],labels=['0%','50%','100%']);fig.text(.5,.023,'2/48 objects; exact-source-barcode joins; broad-lineage concordance; non-blinded.',ha='center',fontsize=8);fig.savefig(F/'FigureS12_source_reference.pdf');fig.savefig(F/'FigureS12_source_reference.png',dpi=600);fig.savefig(F/'FigureS12_source_reference.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
with pd.ExcelWriter(B/'Supplementary_Table_S4_Composition_Annotation_Validation.xlsx',engine='openpyxl',mode='a',if_sheet_exists='replace')as w:
 for sheet,fn in [('T1D_author_raw_label_pairs','T1D_source_author_raw_and_canonical_label_pairs.tsv'),('T1D_author_confusion','T1D_source_author_vs_curated_coarse_confusion_counts.tsv'),('T1D_author_lineage_stats','T1D_source_author_concordance_per_lineage.tsv'),('T1D_author_donor_stats','T1D_source_author_concordance_per_donor.tsv'),('T1D_DC_Mono_disagreements','T1D_source_DC_curated_monocyte_disagreement_by_donor.tsv')]:pd.read_csv(T/fn,sep='\t').to_excel(w,sheet_name=sheet,index=False)
print('Source-reference combined figure and honest label-pair diagnostics completed')

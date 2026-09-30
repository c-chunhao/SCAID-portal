import os
"""Recompute the actual original integrated annotation layers with complete zeros.
Fixed immune denominator from full released metadata + integrated lineage membership.
All plots are explicitly descriptive. No pooled cross-condition ANOVA is retained.
"""
from pathlib import Path
import json,re,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures'
OLD=Path(os.environ['SCAID_SOURCE_AUDIT_DIR'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42})
immune={'Antigen-presenting cells','B cells','B/Plasma cells','pDCs','DC','DCs','Eosinophils','Granulocytes','HSCs','ILC','Mast cells','Mononuclear Phagocytes','Myeloid cells','NK cells','NKT cells','pDC','Plasma cells','T cells','T/NK cells'}
def save(fig,name):
 fig.savefig(F/(name+'.pdf'));fig.savefig(F/(name+'.png'),dpi=600);fig.savefig(F/(name+'.tiff'),dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
m=pd.read_csv(T/'all_release_original_metadata.tsv.gz',sep='\t',keep_default_na=False)
assert len(m)==2572480 and m.Data_object.nunique()==48
m['Data_object']=m.Data_object.str.replace(r'\d+$','',regex=True)
m.loc[m.Data_object.eq('SV_GSE198616_PBMC'),'Data_object']='SV_GSE198616_PBMC_RNA'
hr=m.Data_object.eq('IBD_HRA000072_colon');m.loc[hr,'Data_object']='IBD_'+m.loc[hr,'Disease']+'_HRA000072_colon'
m['Disease']=m.Disease.str.strip().replace({'AH':'AIH','T1DM':'T1D','SS':'SSc'});m.loc[m.Dataset.eq('E-MTAB-8207'),'Disease']='PsA';m.loc[m.Dataset.eq('GSE198616'),'Disease']='BD'
m['key']=m.OldCells+'|'+m.Sample+'|'+m.Dataset
assert m.key.is_unique,'Duplicated atlas records must be reconciled'
lin={}
for name in ['T_NK','B_Plasma']:
 d=pd.read_csv(T/(name+'_original_integration_metadata.tsv.gz'),sep='\t',keep_default_na=False)
 d['key']=d.OldCells+'|'+d.Sample+'|'+d.Dataset;assert d.key.is_unique
 unmatched=~d.key.isin(m.key)
 assert not unmatched.any(), f'{name}: {int(unmatched.sum())} integration cells not in release'
 for c in ['celltype','celltype1','celltype2']:
  d[c]=d[c].replace({'C19':'CD8_MemTex','C13':'CD8_ISGs'}).str.replace('CD8_MenTex','CD8_MemTex',regex=False).str.replace('CD4_Naive','CD4_Tn',regex=False).str.replace('CD4_Effector','CD4_Te',regex=False).str.replace('CD8_Naive','CD8_Tn',regex=False).str.replace('CD8_Memory','CD8_Tm',regex=False)
 lin[name]=d
assert not set(lin['T_NK'].key)&set(lin['B_Plasma'].key)
member=set(lin['T_NK'].key)|set(lin['B_Plasma'].key)
m['denominator_member']=m.MajorCellType.isin(immune)|m.key.isin(member)
unit=m.groupby(['Data_object','Sample'],observed=True).agg(Dataset=('Dataset','first'),Condition=('Disease','first'),Tissue=('SampleTypeFull','first'),All_cells=('key','size'),Immune_denominator=('denominator_member','sum')).reset_index()
assert len(unit)==625
source=pd.read_csv(OLD/'source_cohorts/sample_source_mapping.tsv',sep='\t',keep_default_na=False)
unit=unit.merge(source[['Data_object','Local_sample_label','Source_donor_IDs','Donor_status','Reference_type','Evidence_url']],left_on=['Data_object','Sample'],right_on=['Data_object','Local_sample_label'],validate='one_to_one',how='left')
assert unit.Donor_status.notna().all()
def eligible(x):
 status=x.Donor_status.lower()
 if 'pool' in status:return 'exclude: pooled cells lack individual assignment'
 if not x.Source_donor_IDs.strip():return 'exclude: clinical/source individual key unresolved'
 if any(t in status for t in ['unverified','provisional','independent clinical key not available','alias','partially recoverable']):return 'exclude: identity/independence not established'
 return 'verified source pseudonym; donor grouping available; clinical/visit design still required'
unit['Donor_analysis_eligibility']=unit.apply(eligible,axis=1);unit['Specimen_id']=unit.Data_object+'|'+unit.Sample
unit.to_csv(T/'composition_specimen_denominators_and_eligibility.tsv',sep='\t',index=False)
objs=unit.groupby('Data_object').agg(Condition=('Condition','first'),Source=('Dataset','first'),Tissue=('Tissue','first'),Specimens=('Sample','size'),Cells=('All_cells','sum'),Immune_denominator=('Immune_denominator','sum')).reset_index().sort_values(['Condition','Tissue','Source','Data_object'])
objs['Object_code']=[f'O{i+1:02d}'for i in range(len(objs))];objs.to_csv(T/'composition_object_key.tsv',sep='\t',index=False)
rows=[]
for name,d in lin.items():
 base=d[['key','celltype','celltype1','celltype2']].merge(m[['key','Data_object','Sample']],on='key',validate='one_to_one')
 for layer,col in [('layer1','celltype'),('layer2','celltype1'),('layer3','celltype2')]:
  labels=sorted(base[col].unique());assert '' not in labels
  ct=base.groupby(['Data_object','Sample',col],observed=True).size().rename('Count').reset_index().rename(columns={col:'Cell_type'})
  grid=unit[['Specimen_id','Data_object','Sample','Condition','Tissue','Dataset','Immune_denominator']].merge(pd.DataFrame({'Cell_type':labels}),how='cross')
  z=grid.merge(ct,on=['Data_object','Sample','Cell_type'],how='left',validate='one_to_one');z.Count=z.Count.fillna(0).astype(int)
  z['Proportion']=np.where(z.Immune_denominator.gt(0),z.Count/z.Immune_denominator,np.nan);z['Lineage']=name;z['Layer']=layer
  assert len(z)==625*len(labels);assert z.Proportion.dropna().between(0,1).all()
  z.to_csv(T/(name+'_'+layer+'_complete_zero_specimen_proportions.tsv.gz'),sep='\t',index=False,compression='gzip');rows.append(z)
allp=pd.concat(rows,ignore_index=True)
allp.groupby(['Lineage','Layer','Cell_type','Condition'],observed=True).agg(Specimens=('Specimen_id','size'),Valid_fraction_specimens=('Proportion','count'),Nonzero_specimens=('Count',lambda s:int(s.gt(0).sum())),Median_fraction=('Proportion','median'),Q25=('Proportion',lambda s:s.quantile(.25)),Q75=('Proportion',lambda s:s.quantile(.75)),Cells=('Count','sum')).reset_index().to_csv(T/'all_layers_condition_descriptive_medians.tsv',sep='\t',index=False)
# E-H: same biological annotation layer as legacy Figure 3; all specimens including
# true zeros. Display catalogue distributions without non-identifiable p-values.
cts=[('T_NK','CD4_T'),('T_NK','CD8_T'),('B_Plasma','B cells'),('B_Plasma','Plasma cells')]
conditions=sorted(unit.Condition.unique())
fig,axes=plt.subplots(2,2,figsize=(6.5,7.7),sharey=True)
for ax,(lineage,ct),letter in zip(axes.ravel(),cts,'EFGH'):
 z=allp[(allp.Lineage==lineage)&(allp.Layer=='layer1')&(allp.Cell_type==ct)];assert len(z)==625,(ct,len(z))
 stats=z.groupby('Condition').Proportion.agg(median='median',q25=lambda s:s.quantile(.25),q75=lambda s:s.quantile(.75)).reindex(conditions)
 yy=np.arange(len(conditions));med=stats['median'].to_numpy();ax.barh(yy,med,color='#587f8a',height=.65);ax.errorbar(med,yy,xerr=np.vstack([med-stats.q25,stats.q75-med]),fmt='none',color='#243f47',linewidth=.65,capsize=1.5)
 ax.set_yticks(yy,conditions);ax.invert_yaxis();ax.set_xlim(0,max(.025,float(stats.q75.max())*1.15));ax.set_title(f'{letter}  {ct}',loc='left',fontsize=10);ax.set_xlabel('Median specimen fraction');ax.spines[['top','right']].set_visible(False)
fig.text(.5,.015,'Bars: medians; whiskers: IQR; descriptive catalogue summaries; no disease-effect tests.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.04,1,1],h_pad=2,w_pad=1.2);save(fig,'Figure3_EH_recomputed')
# Readable, single-page heatmaps retain each original subtype label. A cell
# represents the median specimen fraction, with zeros included; source/tissue
# resolution is separately tabulated and plotted with object identifiers.
def heatmap(z,name,title,height):
 v=z.pivot_table(index='Cell_type',columns='Condition',values='Proportion',aggfunc='median').reindex(columns=conditions).fillna(0)
 fig,ax=plt.subplots(figsize=(6.5,height));im=ax.imshow(np.sqrt(v.to_numpy()),cmap='YlGnBu',aspect='auto',vmin=0,vmax=max(.1,float(np.sqrt(v.max().max()))));ax.set_xticks(range(len(v.columns)),v.columns,rotation=90);ax.set_yticks(range(len(v)),v.index);ax.set_title(title,pad=10);ax.set_xlabel('Source-corrected condition label');ax.set_ylabel('Original integrated subtype')
 cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.025);ticks=np.linspace(0,im.norm.vmax,5);cb.set_ticks(ticks,labels=[f'{100*t*t:.1f}%'for t in ticks]);cb.set_label('Median specimen fraction\n(square-root colour scale)',fontsize=8)
 fig.text(.5,.016,'Complete zero counts; fixed annotation-based immune denominator; descriptive only.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.045,1,1]);save(fig,name)
s3=allp[allp.Layer.eq('layer2')].copy();s3['Cell_type']=s3.Lineage.str.replace('_','/')+': '+s3.Cell_type
heatmap(s3,'FigureS3_recomputed','Intermediate subtype composition',7.6)
heatmap(allp[(allp.Lineage=='T_NK')&(allp.Layer=='layer3')],'FigureS8_recomputed','Fine T/NK subtype composition',8.0)
heatmap(allp[(allp.Lineage=='B_Plasma')&(allp.Layer=='layer3')],'FigureS9_recomputed','Fine B/plasma subtype composition',5.5)
# Full study/tissue-stratified summaries; object codes point to the complete key.
with PdfPages(F/'Composition_source_tissue_strata.pdf')as pp:
 for (lineage,layer),z in allp.groupby(['Lineage','Layer'],observed=True):
  v=z.pivot_table(index='Data_object',columns='Cell_type',values='Proportion',aggfunc='median').reindex(objs.Data_object).fillna(0)
  for start in range(0,v.shape[1],12):
   p=v.iloc[:,start:start+12];fig,ax=plt.subplots(figsize=(6.5,8));im=ax.imshow(np.sqrt(p.to_numpy()),aspect='auto',cmap='YlGnBu');ax.set_yticks(range(48),objs.Object_code+' '+objs.Condition);ax.set_xticks(range(p.shape[1]),p.columns,rotation=90);ax.set_title(f'{lineage.replace("_","/")} {layer}: source/tissue strata',fontsize=10);ax.set_xlabel('Original subtype; object identities in accompanying table');fig.colorbar(im,ax=ax,fraction=.025,pad=.02,label='Square root of median fraction');fig.tight_layout();pp.savefig(fig);plt.close(fig)
summary={'release_cells':len(m),'objects':len(objs),'specimen_labels':len(unit),'conditions':conditions,'n_condition_labels':len(conditions),'T_NK_integrated_cells':len(lin['T_NK']),'B_Plasma_integrated_cells':len(lin['B_Plasma']),'immune_denominator_cells':int(m.denominator_member.sum()),'zero_denominator_specimens':int(unit.Immune_denominator.eq(0).sum()),'zero_denominator_records':unit.loc[unit.Immune_denominator.eq(0),['Data_object','Sample','Condition','All_cells','Immune_denominator']].to_dict(orient='records'),'specimens_with_defined_fraction':int(unit.Immune_denominator.gt(0).sum()),'total_grid_rows':len(allp),'zero_count_rows':int(allp.Count.eq(0).sum()),'source_corrections':{'E-MTAB-8207':'PsA','GSE198616':'BD'},'layers':allp.groupby(['Lineage','Layer']).Cell_type.nunique().to_dict().__str__(),'clinical_mapping_eligibility':unit.Donor_analysis_eligibility.value_counts().to_dict(),'statistics':'Descriptive medians/IQR; old ANOVA/F/P labels fully withdrawn. Conditions confounded with studies/tissues; no pooled cross-condition inferential comparison performed.','denominator':'Union of explicitly annotated immune lineages in full released metadata and exact cells of the integrated T/NK or B/plasma objects. Annotation based; not CD45 protein measurement. One fixed denominator per object-scoped specimen applied to every layer/lineage.','identity_scope':'All 625 specimen labels retained descriptively, including pools/unknown clinical identities; none promoted to independent donor replicates. Within-study inference is separate and uses verified source ind_cov.'}
(B/'composition_recomputation_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)

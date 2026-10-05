"""Recompute descriptive integrated-lineage composition directly from exact joins.

Example:
  python recompute_atlas_composition_vNext.py --input-tables PATH \
    --source-map sample_source_mapping.tsv --output NEW_DIRECTORY \
    --regression-tables frozen_corrected_tables --plot

No cells are donors here. No pooled disease tests or missing-value zero imputation.
The source inputs are read-only; outputs must not be inside the input directory.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd

IMMUNE={'Antigen-presenting cells','B cells','B/Plasma cells','pDCs','DC','DCs','Eosinophils','Granulocytes','HSCs','ILC','Mast cells','Mononuclear Phagocytes','Myeloid cells','NK cells','NKT cells','pDC','Plasma cells','T cells','T/NK cells'}
MAJOR={'T_NK':{'T/NK cells','T cells','NK cells','NKT cells'},'B_Plasma':{'B/Plasma cells','B cells','Plasma cells'}}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for chunk in iter(lambda:f.read(8*1024**2),b''):h.update(chunk)
 return h.hexdigest()
def identity_key(frame):
 # Explicitly reject delimiter collisions instead of accepting a false join.
 fields=frame[['OldCells','Sample','Dataset']]
 if fields.isna().any().any():raise ValueError('Missing composite cell identity')
 fields=fields.astype(str)
 if fields.eq('').any().any() or any(fields[col].str.contains('|',regex=False).any() for col in fields):
  raise ValueError('Empty identity or unsupported pipe delimiter in composite cell identity')
 return fields['OldCells']+'|'+fields['Sample']+'|'+fields['Dataset']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input-tables',type=Path,required=True);ap.add_argument('--source-map',type=Path);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--regression-tables',type=Path);ap.add_argument('--plot',action='store_true');args=ap.parse_args()
 inp=args.input_tables.resolve();out=args.output.resolve()
 if out==inp or inp in out.parents:raise ValueError('Output must be separate from the immutable input-tables directory')
 out.mkdir(parents=True,exist_ok=True);tables=out/'tables';tables.mkdir(exist_ok=True);manifest=[]
 def read(name):
  p=inp/name;manifest.append({'path':str(p),'SHA256':sha(p)});return pd.read_csv(p,sep='\t',keep_default_na=False)
 m=read('all_release_original_metadata.tsv.gz');m['Data_object']=m.Data_object.str.replace(r'\d+$','',regex=True)
 m.loc[m.Data_object.eq('SV_GSE198616_PBMC'),'Data_object']='SV_GSE198616_PBMC_RNA'
 h=m.Data_object.eq('IBD_HRA000072_colon');m.loc[h,'Data_object']='IBD_'+m.loc[h,'Disease']+'_HRA000072_colon'
 m['Disease']=m.Disease.str.strip().replace({'AH':'AIH','T1DM':'T1D','SS':'SSc'});m.loc[m.Dataset.eq('E-MTAB-8207'),'Disease']='PsA';m.loc[m.Dataset.eq('GSE198616'),'Disease']='BD'
 m['key']=identity_key(m);assert m.key.is_unique
 lineage={};members=set()
 for name in ['T_NK','B_Plasma']:
  d=read(name+'_original_integration_metadata.tsv.gz');d['key']=identity_key(d)
  assert d.key.is_unique and d.key.isin(m.key).all()
  assert not set(d.key)&members;members.update(d.key)
  for col in ['celltype','celltype1','celltype2']:
   d[col]=d[col].replace({'C19':'CD8_MemTex','C13':'CD8_ISGs'}).str.replace('CD8_MenTex','CD8_MemTex',regex=False).str.replace('CD4_Naive','CD4_Tn',regex=False).str.replace('CD4_Effector','CD4_Te',regex=False).str.replace('CD8_Naive','CD8_Tn',regex=False).str.replace('CD8_Memory','CD8_Tm',regex=False)
  d=d[['key','celltype','celltype1','celltype2']].merge(m[['key','Data_object','Sample']],on='key',validate='one_to_one');lineage[name]=d
 m['immune']=m.MajorCellType.isin(IMMUNE)|m.key.isin(members)
 unit=m.groupby(['Data_object','Sample'],observed=True).agg(Dataset=('Dataset','first'),Condition=('Disease','first'),Tissue=('SampleTypeFull','first'),All_cells=('key','size'),Immune_denominator=('immune','sum')).reset_index();unit['Specimen_id']=unit.Data_object+'|'+unit.Sample
 assert unit.Specimen_id.is_unique
 if args.source_map:
  source=pd.read_csv(args.source_map,sep='\t',keep_default_na=False);manifest.append({'path':str(args.source_map.resolve()),'SHA256':sha(args.source_map)})
  sample_key='Sample'if 'Sample'in source else 'Local_sample_label'
  if sample_key!='Sample':source=source.rename(columns={sample_key:'Sample'})
  cols=['Data_object','Sample']+[c for c in ['Source_donor_IDs','Donor_status','Reference_type','Evidence_url','Donor_analysis_eligibility']if c in source]
  unit=unit.merge(source[cols],on=['Data_object','Sample'],how='left',validate='one_to_one')
  assert unit.Donor_status.notna().all() if 'Donor_status'in unit else True
 unit.to_csv(tables/'composition_specimen_denominators_and_eligibility.tsv',sep='\t',index=False)
 cov=[];allp=[];reg=[]
 for name,d in lineage.items():
  expected=m[m.MajorCellType.isin(MAJOR[name])].groupby(['Data_object','Sample']).size().rename('Source_major_lineage_cells').reset_index()
  observed=d.groupby(['Data_object','Sample']).size().rename('Integrated_lineage_cells').reset_index()
  c=unit[['Data_object','Sample']].merge(expected,on=['Data_object','Sample'],how='left').merge(observed,on=['Data_object','Sample'],how='left')
  c[['Source_major_lineage_cells','Integrated_lineage_cells']]=c[['Source_major_lineage_cells','Integrated_lineage_cells']].fillna(0).astype(int)
  c['Lineage']=name;c['coverage_status']=np.where(c.Source_major_lineage_cells.gt(0)&c.Integrated_lineage_cells.eq(0),'source lineage present but integration unassessed',np.where(c.Integrated_lineage_cells.gt(0),'integrated lineage assessed','no source major-lineage cells'));cov.append(c)
  for layer,col in [('layer1','celltype'),('layer2','celltype1'),('layer3','celltype2')]:
   labels=sorted(d[col].unique());assert ''not in labels
   counts=d.groupby(['Data_object','Sample',col],observed=True).size().rename('Count').reset_index().rename(columns={col:'Cell_type'})
   z=unit[['Specimen_id','Data_object','Sample','Condition','Tissue','Dataset','Immune_denominator']].merge(pd.DataFrame({'Cell_type':labels}),how='cross').merge(counts,on=['Data_object','Sample','Cell_type'],how='left',validate='one_to_one')
   # This zero applies to the observed count grid; inferential availability is
   # explicitly separate and is applied before proportions/summaries are written.
   z['Count']=z.Count.fillna(0).astype(int);z['Lineage']=name;z['Layer']=layer
   z=z.merge(c,on=['Data_object','Sample','Lineage'],validate='many_to_one')
   z['Original_derived_proportion']=np.where(z.Immune_denominator.gt(0),z.Count/z.Immune_denominator,np.nan)
   z['Proportion']=z.Original_derived_proportion.mask(z.coverage_status.eq('source lineage present but integration unassessed'))
   assert z.Proportion.dropna().between(0,1).all() and len(z)==len(unit)*len(labels)
   z['Fraction_availability']=np.where(z.coverage_status.eq('source lineage present but integration unassessed'),'NA: source lineage positive but integrated lineage unassessed',np.where(z.Proportion.notna(),'assessed annotation-based fraction; true subtype zeros retained','NA: no immune denominator'))
   fname=f'{name}_{layer}_coverage_masked_specimen_proportions.tsv.gz';z.to_csv(tables/fname,sep='\t',index=False,compression='gzip');allp.append(z)
   if args.regression_tables:
    reference=pd.read_csv(args.regression_tables/fname,sep='\t');keys=['Data_object','Sample','Cell_type'];a=z.sort_values(keys).reset_index(drop=True);b=reference.sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(a[keys],b[keys]);assert np.array_equal(a.Count.to_numpy(),b.Count.to_numpy());assert np.array_equal(a.Immune_denominator.to_numpy(),b.Immune_denominator.to_numpy());assert np.allclose(a.Proportion,b.Proportion,equal_nan=True,rtol=0,atol=1e-15)
    reg.append({'lineage':name,'layer':layer,'rows':len(z),'count_denominator_and_masked_fraction_roundtrip':True})
 coverage=pd.concat(cov,ignore_index=True);coverage.to_csv(tables/'Integrated_lineage_source_coverage_QA.tsv',sep='\t',index=False)
 allp=pd.concat(allp,ignore_index=True)
 summary=allp.groupby(['Lineage','Layer','Cell_type','Condition']).agg(All_specimen_labels=('Specimen_id','size'),Assessed_fraction_specimens=('Proportion','count'),Unassessed_source_positive_specimens=('coverage_status',lambda x:x.eq('source lineage present but integration unassessed').sum()),Nonzero_assessed_specimens=('Proportion',lambda x:x.gt(0).sum()),Median_fraction=('Proportion','median'),Q25=('Proportion',lambda x:x.quantile(.25)),Q75=('Proportion',lambda x:x.quantile(.75))).reset_index();summary.to_csv(tables/'all_layers_coverage_masked_condition_descriptive.tsv',sep='\t',index=False)
 strata=allp.groupby(['Lineage','Layer','Data_object','Cell_type','Condition','Dataset','Tissue']).agg(Specimens=('Specimen_id','size'),Assessed=('Proportion','count'),Median_fraction=('Proportion','median'),Q25=('Proportion',lambda x:x.quantile(.25)),Q75=('Proportion',lambda x:x.quantile(.75))).reset_index();strata.to_csv(tables/'composition_coverage_masked_source_tissue_strata.tsv.gz',sep='\t',index=False,compression='gzip')
 if args.plot:
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  plots=out/'figures';plots.mkdir(exist_ok=True);cmap=plt.get_cmap('YlGnBu').copy();cmap.set_bad('#d9d9d9')
  for (line,layer),z in summary.groupby(['Lineage','Layer']):
   # pivot preserves all-NA columns/rows and measured zeros; never fill missing.
   v=z.pivot(index='Cell_type',columns='Condition',values='Median_fraction');fig,ax=plt.subplots(figsize=(9,max(3,len(v)*.2)));im=ax.imshow(np.ma.masked_invalid(np.sqrt(v.to_numpy())),cmap=cmap,aspect='auto');ax.set_xticks(range(len(v.columns)),v.columns,rotation=90);ax.set_yticks(range(len(v)),v.index);ax.set_title(f'{line} {layer}: assessed specimen medians; grey = NA');fig.colorbar(im,ax=ax,label='Square root of assessed median fraction');fig.tight_layout();fig.savefig(plots/f'{line}_{layer}_coverage_masked.pdf');plt.close(fig)
 report={'passed':True,'source_release_cell_records':len(m),'source_objects':m.Data_object.nunique(),'specimen_labels':len(unit),'source_positive_unassessed_by_lineage':coverage[coverage.coverage_status.eq('source lineage present but integration unassessed')].groupby('Lineage').size().to_dict(),'defined_fractions_per_lineage':allp[allp.Layer.eq('layer1')].groupby('Lineage').apply(lambda x:int(x.groupby('Specimen_id').Proportion.first().notna().sum()),include_groups=False).to_dict(),'numeric_roundtrip_regression':reg,'source_inputs':manifest,'immune_categories':sorted(IMMUNE),'source_major_lineage_vocabularies':{k:sorted(v)for k,v in MAJOR.items()},'scope':'Descriptive retained integration membership / full-release annotation immune denominator. Source-positive unassessed integration fractions are NA, never biological zero. No donor independence or pooled disease-effect test implied.'}
 (out/'composition_recomputation_QA.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()

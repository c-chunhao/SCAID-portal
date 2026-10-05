"""Source-conditional marker review screens; never automatic annotation fixes."""
from pathlib import Path
import json,argparse
import numpy as np,pandas as pd
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input-dir',type=Path,required=True,help='Audit output directory; no raw counts are reread.')
R=ap.parse_args().input_dir.resolve()
PANELS={'T_core':['TRAC','CD3D'],'NK_core':['NKG7','GNLY','KLRD1'],'CD8_core':['CD8A','CD8B'],'Naive_memory':['CCR7','IL7R'],'Treg':['FOXP3','IL2RA'],'Th1':['TBX21'],'Th17':['RORC','KLRB1'],'Tfh':['CXCR5','PDCD1'],'Exhaustion':['PDCD1','LAG3'],'MAIT_support':['KLRB1','ZBTB16','TRAV1-2'],'Gamma_delta':['TRDC','TRDV2'],'B_core':['MS4A1','CD79A'],'B_naive':['TCL1A'],'GC_support':['BCL6','AICDA'],'Plasma_support':['JCHAIN','MZB1','SDC1'],'IgA':['IGHA1'],'IgG':['IGHG1'],'Cycling':['MKI67'],'Interferon':['ISG15']}
RULES={'minimum_source_or_specimen_label_cells':20,'review_only_not_accuracy_thresholds':True,'rule_version':2,'method_amendment':'TRDC-only initial screen was removed after primary-method review: NK/ILC can express TRDC strongly. Final gamma-delta RNA review requires CD3D/TRDC or CD3D/TRDV2 same-cell Frechet lower bound >=0.10 and never classifies gamma-delta identity. Raw scan and the predeclared 40-gene panel did not change.','primary_reference':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6576116/','NK_T_marker_visible':'T_core mean detection >=0.25','NKT_T_marker_weak':'T_core mean detection <0.10','ILC_T_marker_visible':'T_core mean detection >=0.25','gamma_delta_RNA_candidate':'max(0, detection(CD3D)+detection(TRDC or TRDV2)-1)>=0.10 for at least one available TRD gene; exploratory only','GC_canonical_support_low':'GC_support mean detection <0.02, provided both genes exist on full RNA axis','GC_MKI67_plus_low_cycle':'MKI67+ label but MKI67 detection <0.05','plasmablast_plasma_program_low':'Plasma_support mean detection <0.05, provided all genes exist','CD8_Tn_cytotoxic_marker':'CD8_Tn label with GZMH detection >=0.30 and CCR7 detection <0.10','interpretation':'Conservative descriptive review flags, not calibrated failure cutoffs. Dropout, activation, tissue/source composition and capture technology affect markers; absence alone cannot disprove a label. MAIT/NKT/gamma-delta identity may require TCR/protein evidence.'}
(R/'predeclared_marker_review_rules.json').write_text(json.dumps({'panels':PANELS,'rules':RULES},indent=2))
def review(data,kind):
 keys=['Data_object','Dataset','Condition','Tissue','Lineage','celltype','celltype1','celltype2']
 if kind=='specimen':keys+=['Sample','Source_donor_IDs','Donor_status','Donor_analysis_eligibility']
 detection=data.pivot(index=keys,columns='Marker',values='Detection_fraction')
 expr=data.pivot(index=keys,columns='Marker',values='Mean_log1p_10k')
 cells=data.groupby(keys,dropna=False,observed=True).Cells.first()
 q=detection.index.to_frame(index=False);q['Cells']=cells.reindex(detection.index).to_numpy();q['Review_min20_cells']=q.Cells.ge(20)
 for name,genes in PANELS.items():
  z=detection.reindex(columns=genes);e=expr.reindex(columns=genes);complete=z.notna().all(axis=1)
  q[name+'_available']=complete.to_numpy();q[name+'_mean_detection']=z.mean(axis=1).where(complete).to_numpy();q[name+'_mean_log1p10k']=e.mean(axis=1).where(complete).to_numpy()
 # Frechet bounds are mathematically valid from gene marginals; they are not
 # a measured joint-detection matrix or proof against doublets.
 T=detection.reindex(columns=PANELS['T_core']);NK=detection.reindex(columns=PANELS['NK_core']);ok=T.notna().all(axis=1)&NK.notna().all(axis=1)
 q['T_NK_same_cell_detection_lower_bound']=np.maximum(0,T.max(axis=1)+NK.max(axis=1)-1).where(ok).to_numpy()
 q['T_NK_same_cell_detection_upper_bound']=np.minimum(np.minimum(1,T.sum(axis=1)),np.minimum(1,NK.sum(axis=1))).where(ok).to_numpy()
 for gene in ['TRDC','TRDV2']:
  cd3=detection['CD3D'];trd=detection[gene];q['CD3D_'+gene+'_same_cell_detection_lower_bound']=np.maximum(0,cd3+trd-1).where(cd3.notna()&trd.notna()).to_numpy()
 flags=[]
 for idx,row in q.iterrows():
  f=[];label=row.celltype2
  if row.Cells<20:f.append('too_few_cells_for_marker_review')
  else:
   if label.startswith('NK_')and row.T_core_mean_detection>=.25:f.append('NK_label_with_visible_T_RNA_markers')
   if label.startswith('NKT_')and row.T_core_mean_detection<.10:f.append('NKT_label_with_weak_T_RNA_markers')
   if label=='ILC'and row.T_core_mean_detection>=.25:f.append('ILC_label_with_visible_T_RNA_markers')
   gamma=[v for v in [row.CD3D_TRDC_same_cell_detection_lower_bound,row.CD3D_TRDV2_same_cell_detection_lower_bound]if np.isfinite(v)]
   if gamma and max(gamma)>=.10:f.append('T_CD3D_TRD_joint_bound_review_candidate')
   if label.startswith('Bgc')and row.GC_support_mean_detection<.02:f.append('GC_label_low_BCL6_AICDA_detection')
   if label=='Bgc_MKI67+'and detection.iloc[idx].get('MKI67',np.nan)<.05:f.append('GC_MKI67_plus_low_cycle_marker')
   if label=='Plasmablast'and row.Plasma_support_mean_detection<.05:f.append('Plasmablast_low_plasma_program')
   if label.startswith('CD8_Tn')and detection.iloc[idx].get('GZMH',np.nan)>=.30 and detection.iloc[idx].get('CCR7',np.nan)<.10:f.append('CD8_naive_label_cytotoxic_marker_review')
  if not row.T_core_available and row.Lineage=='T_NK':f.append('T_core_gene_coverage_incomplete')
  if not row.GC_support_available and label.startswith('Bgc'):f.append('GC_marker_gene_coverage_incomplete')
  if not row.Plasma_support_available and label=='Plasmablast':f.append('plasma_marker_gene_coverage_incomplete')
  flags.append(';'.join(f))
 q['Review_flags']=flags;q.to_csv(R/f'{kind}_fine_label_coherence_review.tsv.gz',sep='\t',index=False,compression='gzip')
 if kind=='source':q[q.Review_flags.ne('')&q.Review_min20_cells].to_csv(R/'source_review_candidates.tsv',sep='\t',index=False)
 return q
if __name__=='__main__':
 source=pd.read_csv(R/'all_source_fine_label_marker_evidence.tsv.gz',sep='\t',keep_default_na=True,low_memory=False);specimen=pd.read_csv(R/'all_specimen_fine_label_marker_evidence.tsv.gz',sep='\t',keep_default_na=False,low_memory=False)
 # Empty clinical identity strings are explicit unresolved labels, not missing
 # grouping keys to drop. Only numeric unavailable gene metrics become NA.
 for col in ['Mean_raw_count','Mean_log1p_10k','Detection_fraction']:
  specimen[col]=pd.to_numeric(specimen[col],errors='coerce')
 a=review(source,'source');b=review(specimen,'specimen')
 special=a[a.celltype2.str.contains('NK_|NKT_|ILC|MAIT|Bgc_|Plasmablast')];special.to_csv(R/'NK_NKT_ILC_MAIT_GC_Plasmablast_source_evidence.tsv',sep='\t',index=False)
 gamma=a[a.Review_flags.str.contains('T_CD3D_TRD')&a.Review_min20_cells];gamma.to_csv(R/'gamma_delta_RNA_candidates_source_evidence.tsv',sep='\t',index=False)
 support=a.groupby(['Lineage','celltype2']).agg(Cells=('Cells','sum'),Source_objects=('Data_object','nunique'),Source_groups_ge20=('Review_min20_cells','sum')).reset_index();support.to_csv(R/'all_63_fine_label_audit_coverage.tsv',sep='\t',index=False)
 summary={'source_label_groups':len(a),'source_label_groups_ge20':int(a.Review_min20_cells.sum()),'specimen_label_groups':len(b),'specimen_label_groups_ge20':int(b.Review_min20_cells.sum()),'63_labels_total_cells':int(support.Cells.sum()),'specialist_labels_cells':special.groupby('celltype2').Cells.sum().astype(int).to_dict(),'source_candidates_review_only':int((a.Review_flags.ne('')&a.Review_min20_cells).sum()),'flags':a[a.Review_min20_cells].Review_flags.str.split(';').explode().value_counts().drop('',errors='ignore').astype(int).to_dict(),'no_automatic_relabels':True,'marginal_joint_detection_bounds_note':'T/NK same-cell fractions are Frechet lower/upper bounds computed from raw-gene marginals, not direct measured joint detection or evidence excluding doublets.','missing_feature_values_remain_NA':True}
 (R/'marker_review_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)

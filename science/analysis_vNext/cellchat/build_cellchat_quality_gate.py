"""Publish exact provenance and defensible descriptive network summaries.
Inference gates fail closed when archived pooled networks lack donor keys,
database identity, matching source/tissue, or tested pathway coverage.
"""
from pathlib import Path
import sys, os
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime_paths import required_path, registered_datasets
import hashlib,json
import numpy as np
import pandas as pd
ROOT=required_path('SCAID_CELLCHAT_AUDIT_DIR');T=ROOT/'tables'
S=required_path('SCAID_SCIENCE_TABLES')
q=pd.read_csv(T/'cellchat_object_QA_complete_unique.tsv',sep='\t');p=pd.read_csv(T/'cellchat_evaluated_pathway_descriptive_complete_unique.tsv',sep='\t')
m=pd.read_csv(S/'composition_specimen_denominators_and_eligibility.tsv',sep='\t',keep_default_na=False)
q['Data_object']=q.object_stem.replace({'SV_GSE198616_PBMC':'SV_GSE198616_PBMC_RNA','IBD_Crohns_Disease_HRA000072_colon':'IBD_CD_HRA000072_colon','IBD_Ulcerative_Colitis_HRA000072_colon':'IBD_UC_HRA000072_colon','IBD_Colitis_inflamed_HRA000072_colon':'IBD_CI_HRA000072_colon'})
o=m.groupby('Data_object').agg(Condition=('Condition','first'),Dataset=('Dataset','first'),Tissue=('Tissue','first'),released_cells=('All_cells','sum'),specimens=('Sample','size'),donor_mapping_eligible_specimens=('Donor_analysis_eligibility',lambda s:s.str.startswith('verified').sum())).reset_index()
q=q.merge(o,on='Data_object',validate='one_to_one');assert len(q)==48 and q.Data_object.nunique()==48
q['cells_difference_network_minus_release']=q.n_cells-q.released_cells
q['descriptive_gate_pass']=q.probability_finite&q.probability_nonnegative&q.pvalues_in_range&q.network_dimensions_match&q.n_lr_missing_DB_map.eq(0)
q['donor_case_control_inference_gate_pass']=False
q['inference_block_reason']='archived network pools specimens; verified donor and case/control strata not retained as independent networks'
q['cross_disease_strength_gate_pass']=False
q['strength_block_reason']='source/tissue/number of groups/LR database differ; raw probability sums cannot be interpreted as disease effects'
q.to_csv(T/'CellChat_48_object_provenance_and_inference_gate.tsv',sep='\t',index=False)
p=p.merge(q[['object_stem','Data_object','Condition','Dataset','Tissue','lr_database_SHA256','lr_database_rows','cells_difference_network_minus_release']],on='object_stem',validate='many_to_one')
assert not p.duplicated(['Data_object','pathway']).any()
p['coverage_status']='evaluated in this object; zero significant_edges means no detected conditional edge, not biological absence'
p['missing_pathway_policy']='pathway absent from this table is unassessed/unmatched, not imputed biological zero'
p['clinical_effect_inference']=False
p.to_csv(T/'CellChat_source_tissue_database_stratified_descriptive.tsv',sep='\t',index=False)
# Coverage matrix keeps missing/unmatched pathways as NA, and forbids a misleading
# pooled disease average across unmatched studies/tissues or database generations.
coverage=p.pivot(index='Data_object',columns='pathway',values='significant_edges').reindex(q.Data_object)
coverage.to_csv(T/'CellChat_evaluated_pathway_edge_counts_NA_not_zero.tsv',sep='\t',na_rep='NA')
families=q.groupby(['lr_database_SHA256','lr_database_rows']).agg(objects=('Data_object','size'),cells=('n_cells','sum')).reset_index()
families.to_csv(T/'CellChat_database_generation_counts.tsv',sep='\t',index=False)
summary={'passed_descriptive_structural_audit':bool(q.descriptive_gate_pass.all()),'objects':len(q),'native_network_cells':int(q.n_cells.sum()),'released_cells':int(q.released_cells.sum()),'original_clinical_donor_keys_saved':int(q.n_verified_individual_keys.notna().sum()),'metadata_source_version_note':'Original serialized CellChat package version was not stored. Installed audit reader version is not substituted for generation version.','database_families':families.to_dict('records'),'objects_network_cell_count_differs_release':q.loc[q.cells_difference_network_minus_release.ne(0),['Data_object','n_cells','released_cells','cells_difference_network_minus_release']].to_dict('records'),'object_filename_collision_fixed':'3 HRA000072 objects share a legacy basename; exact full-path manifest gives separate Crohn, UC and inflamed artifacts. No first-file glob fallback.','actual_parameter_note':'Saved options recorded triMean/raw.use=FALSE/population.size=FALSE/100 permutations/seed=1. Parameters are frozen in original object options, with source package version unavailable.','pvalue_limitation':'Permutation p-values are conditional on pooled cells/groups, not donor disease tests. Old nboot=100 p=0 is Monte Carlo resolution, never proof of true zero; no cross-disease differential p-value published.','inference_gate':'No archived pooled object is promoted to donor-level case/control or cross-disease strength inference. All 48 remain descriptive exploratory networks.','missing_value_policy':'Only evaluated pathways can carry zero. Unassessed/unmatched pathways are NA. All summaries retain source, tissue, exact LR database hash, and original cell-group opportunities.'}
(ROOT/'CellChat_quality_gate.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))

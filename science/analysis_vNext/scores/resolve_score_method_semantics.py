"""Separate measured schema validity from unrecorded legacy generation versions.
Uncentred singscore under minimum-rank ties can legitimately be slightly negative.
"""
from pathlib import Path
import sys, os
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime_paths import required_path, registered_datasets
import json,hashlib
import h5py
import numpy as np
import pandas as pd
R=required_path('SCAID_SCORE_AUDIT_DIR');T=R/'tables'
f=T/'Score_H5_144_full_matrix_QA.tsv';x=pd.read_csv(f,sep='\t');assert len(x)==144
reg=registered_datasets();ng={}
for d in reg:
 with h5py.File(d['h5_paths']['RNA'],'r')as h:ng[Path(d['h5_paths']['RNA']).stem]=len(h['names_var'])
x['RNA_full_gene_universe_size']=x.object.map(ng)
x['implementation_consistency_gate_pass']=x.range_gate_pass.copy();x['range_generation_status']='bounded score consistent with known method'
ix=x.method.eq('singscore');N=x.loc[ix,'RNA_full_gene_universe_size'];lower=-499/(2*(N-500))
x.loc[ix,'expected_range']=[f'uncentred min-rank implementation conservative envelope [{a:.10f},1]; exact legacy generation version unrecorded'for a in lower]
x.loc[ix,'implementation_consistency_gate_pass']=(x.loc[ix,'stored_and_implicit_zero_min'].to_numpy()>=lower.to_numpy()-1e-8)&x.loc[ix,'stored_and_implicit_zero_max'].le(1+1e-8)
x['range_gate_pass']=x.range_gate_pass.astype(object);x.loc[ix,'range_gate_pass']=None
x.loc[ix,'range_generation_status']='unverified original generation version; consistent with inspected irGSEA centerScore=FALSE and singscore tiesMethod=min'
x['frozen_KEGG_membership_recovered']=True;x['frozen_KEGG_feature_coverage_quantified']=True;x['cross_version_generation_gate_pass']=False
x['gene_set_gate_reason']='Exact frozen 20240526 KEGG membership recovered and hashes agree across 4 archives; original scoring software versions and resolved rank thresholds not embedded in old H5'
x.to_csv(f,sep='\t',index=False)
assert x.schema_gate_pass.all() and x.RNA_obs_order_equal.all() and x.pathway_order_equal_across_3_methods.all()
info={'audited_files':144,'schema_gate_pass_files':int(x.schema_gate_pass.sum()),'finite_values_all_files':bool(x.n_nonfinite_stored.eq(0).all()),'RNA_obs_aligned_files':int(x.RNA_obs_order_equal.sum()),'pathway_order_consistent_files':int(x.pathway_order_equal_across_3_methods.sum()),'bounded_AUCell_UCell_files':96,'singscore_original_generation_range_unverified_files':48,'inspected_implementation_consistency_pass_files':int(x.implementation_consistency_gate_pass.sum()),'true_structural_or_nonfinite_failures':[],'singscore_observed_min':float(x.loc[ix,'stored_and_implicit_zero_min'].min()),'singscore_observed_max':float(x.loc[ix,'stored_and_implicit_zero_max'].max()),'singscore_explicit_negative_values':int(x.loc[ix,'n_explicit_negative'].sum()),'singscore_range_interpretation':'Do not apply generic [-0.5,0.5]. Inspected irGSEA calls simpleScore(centerScore=FALSE), rankGenes default tiesMethod=min. Under minimum-rank ties, mean ranks can fall below the unique-rank lower bound; small negative uncentred values are mathematically compatible. Never clip these to zero or treat scores as RNA counts.','tie_lower_bound_derivation':'For N full RNA genes and m signature genes, uncentred score=(mean_rank-(m+1)/2)/(N-m). Minimum mean rank under ties is 1, so conservative lower bound=-(m-1)/(2*(N-m)); for source maxGSSize500 use m=500 envelope. This is an implementation consistency bound, not certification of unknown original versions.','generation_provenance_gate_pass':False,'generation_limitations':'Exact generation-library versions/resolved rank parameters were not saved. Full frozen KEGG membership recovered; all 48*359 feature-universe coverage rows quantified. Source scoring uses full RNA data assay rather than HVG2000; msigdb=FALSE.','installed_implementation_evidence':{'irGSEA':'3.3.4','singscore':'1.22.0','source_original_script':'AID_database/流程代码与数据/单个样本流程代码.R:368-386','official_irGSEA_source':'https://github.com/chuiqin/irGSEA/blob/master/R/irGSEA.score.R','official_singscore_manual':'https://www.bioconductor.org/packages/release/bioc/manuals/singscore/man/singscore.pdf'}}
(R/'score_H5_QA.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))

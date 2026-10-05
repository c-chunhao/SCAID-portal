"""Quantify the actual frozen 2024-05-26 KEGG set membership for all 48 assets.
Feature-universe coverage is not reconstructed expression detection or a re-score.
"""
from pathlib import Path
import sys, os
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime_paths import required_path, registered_datasets
import hashlib,json
import h5py
import pandas as pd
R=required_path('SCAID_SCORE_AUDIT_DIR');T=R/'tables'
ds=registered_datasets()
g=pd.read_csv(T/'KEGG_20240526_frozen_membership.tsv',sep='\t');g['feature_name']=g.pathway.str.replace('_','-',regex=False);sets=g.groupby('feature_name').gene.agg(set).to_dict();assert len(sets)==359
snapshots=[]
snapshot_paths=[Path(x) for x in os.environ.get('SCAID_KEGG_SNAPSHOT_PATHS','').split(os.pathsep) if x]
assert len(snapshot_paths)==4, 'Provide four local archived KEGG snapshots in SCAID_KEGG_SNAPSHOT_PATHS'
for p in snapshot_paths:
 snapshots.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert len({x['sha256'] for x in snapshots})==1
rows=[]
for item in ds:
 with h5py.File(item['h5_paths']['RNA'],'r')as h:
  genes=set(h['names_var'].asstr()[:]);hvg=set(h['var/var/_index'].asstr()[:])
 with h5py.File(item['h5_paths']['AUCell'],'r')as h:observed=set(h['matrix/pathway'].asstr()[:])
 for key,signature in sets.items():
  present=signature&genes;presenthvg=signature&hvg;cov=len(present)/len(signature)
  rows.append({'object':Path(item['h5_paths']['RNA']).stem,'pathway':key,'frozen_set_size':len(signature),'present_in_full_RNA_gene_universe':len(present),'full_RNA_feature_coverage':cov,'present_in_HVG2000':len(presenthvg),'HVG_feature_coverage':len(presenthvg)/len(signature),'score_feature_present':key in observed,'screen_low_coverage':len(present)<10 or cov<.2,'coverage_interpretation':'feature presence only; exact min.cells>=3 filtered scoring universe/ranks not embedded in legacy H5','membership_sha256':snapshots[0]['sha256'],'present_genes':';'.join(sorted(present)),'missing_genes':';'.join(sorted(signature-present))})
x=pd.DataFrame(rows);assert len(x)==48*359
x.to_csv(T/'KEGG_20240526_48_object_full_RNA_gene_coverage.tsv.gz',sep='\t',index=False,compression='gzip')
summary=x.groupby('object').agg(frozen_sets=('pathway','size'),score_features=('score_feature_present','sum'),mean_full_RNA_coverage=('full_RNA_feature_coverage','mean'),median_full_RNA_coverage=('full_RNA_feature_coverage','median'),mean_HVG2000_coverage=('HVG_feature_coverage','mean'),low_coverage_feature_screen=('screen_low_coverage','sum')).reset_index()
summary.to_csv(T/'KEGG_48_object_coverage_QA_summary.tsv',sep='\t',index=False)
(R/'KEGG_frozen_membership_coverage_QA.json').write_text(json.dumps({'passed':True,'objects':48,'frozen_sets':359,'coverage_rows':len(x),'membership_snapshots':snapshots,'original_source_semantics':'irGSEA.score custom=TRUE, geneset=hsa, msigdb=FALSE, assay=RNA, slot=data; full RNA assay input, not VariableFeatures/HVG input; min.cells=3, min.feature=0, seed=2024, method=AUCell/UCell/singscore, MaxRank=NULL (defaults).','missing_feature_policy':'359 frozen sets; legacy score assets generally omit hsa01100 Metabolic pathways. Omitted score features remain unavailable, never score=0.','low_coverage_screen':'<10 present genes or <20% feature coverage is an explicit exploratory warning screen, not a published method cutoff or a claim these scores are wrong.','exact_generation_provenance_gate':'Exact frozen membership recovered and hashes agree across four old release folders. Original scoring-library versions, effective min.cells filtered gene universe and exact resolved rank thresholds were not embedded in H5; cross-version numerical reproducibility requires re-score and a recorded manifest.','version_note':'UCell >=2.7.6 changed normalization; installed current version is not claimed as original generation version.'},indent=2)+'\n')
print(summary.to_string(index=False));print('48-object frozen KEGG coverage completed.')

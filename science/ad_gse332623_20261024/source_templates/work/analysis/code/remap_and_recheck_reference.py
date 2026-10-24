#!/usr/bin/env python3
"""Correct only reference-label interpretation, retaining raw predictions and embeddings."""
from pathlib import Path
import runpy,json,shutil
import numpy as np,pandas as pd,h5py
M=runpy.run_path(str(Path(__file__).with_name('analysis_pipeline.py')),run_name='analysis_functions')
OUT=M['ROOT']/'integration';W=M['WORK'];P=OUT/'labels/reference_predictions.tsv.gz';D=OUT/'annotation_evidence/PCA_Harmony_diagnostics.json';E=OUT/'checkpoints/embeddings.npz'
backup=OUT/'provenance/reference_predictions.pre_mapping_fix.tsv.gz'
if not backup.exists():shutil.copy2(P,backup)
old=pd.read_csv(P,sep='\t',index_col=0);new=M['map_reference_predictions'](old)
raw=['high_label','high_probability','high_margin','low_label','low_probability','low_margin']
assert new.index.equals(old.index) and new[raw].equals(old[raw])
tmp=P.with_name(P.name+'.tmp');new.to_csv(tmp,sep='\t',compression='gzip');tmp.replace(P)
M['write_json'](OUT/'checkpoints/reference_classification.json',M['reference_mapping_receipt'](new))
with h5py.File(OUT/'checkpoints/normalized_RNA.h5ad','r')as f:
 g=f['obs/library_id'];libs=g['categories'].asstr()[:][g['codes'][:]] if isinstance(g,h5py.Group)else g.asstr()[:]
 key=f['obs'].attrs['_index'];cells=f['obs'][key].asstr()[:]
obs=pd.DataFrame({'library_id':libs},index=cells);assert obs.index.equals(new.index)
states=pd.read_csv(OUT/'annotation_evidence/diagnostic_RNA_state_mean_log_expression.tsv.gz',sep='\t',index_col=0);assert states.index.equals(new.index)
old_diag=json.loads(D.read_text());old_embedding=json.loads((OUT/'checkpoints/embedding.json').read_text());e=np.load(E)
new_diag=M['representation_diagnostics'](e['X_pca'],e['X_harmony'],obs,new,states,OUT)
new_diag['reference_mapping_policy_version']=2;new_diag['reference_mapping_correction']='High ILC compatible with Low NK only, explicitly separate from exact coarse agreement; High Cycling state is not independent lineage agreement.'
M['write_json'](OUT/'annotation_evidence/PCA_Harmony_diagnostics.corrected_reference_candidate.json',new_diag)
assert new_diag['selected_representation']==str(e['selected_representation'][0])==old_embedding['selected_representation'], 'Corrected references changed representation choice; rebuild only the graph before finalize.'
if not(OUT/'provenance/PCA_Harmony_diagnostics.pre_mapping_fix.json').exists():shutil.copy2(D,OUT/'provenance/PCA_Harmony_diagnostics.pre_mapping_fix.json')
M['write_json'](D,new_diag)
old_embedding['reference_mapping_policy_version']=2;old_embedding['diagnostics_recomputed_after_reference_mapping_fix']=True;M['write_json'](OUT/'checkpoints/embedding.json',old_embedding)
report={'status':'validated','raw_model_predictions_unchanged':True,'counts_or_feature_axes_changed':False,'embeddings_changed':False,'reference_mapping_policy':M['reference_mapping_receipt'](new),'selected_representation_before':old_diag['selected_representation'],'selected_representation_after':new_diag['selected_representation'],'new_diagnostics':{k:v for k,v in new_diag.items()if k!='metrics'},'predictions_sha256':M['sha'](P),'embedding_file_sha256':M['sha'](E),'pipeline_code_sha256':M['sha'](Path(__file__).with_name('analysis_pipeline.py'))}
M['write_json'](OUT/'provenance/reference_mapping_correction.json',report);M['write_json'](W/'validation/reference_mapping_correction.json',report);print(json.dumps(report,indent=2))

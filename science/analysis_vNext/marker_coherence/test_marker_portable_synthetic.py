"""Three explicit synthetic cells; no project biological data are read.

Example: python test_marker_portable_synthetic.py --marker-code-dir .
  --reader ../integration/released_rna.py --output /tmp/marker_synthetic_check
"""
from pathlib import Path
import argparse,json,subprocess,sys,shutil
import h5py,numpy as np,pandas as pd
from scipy import sparse
ap=argparse.ArgumentParser(description=__doc__)
for n in ['marker-code-dir','reader','output']:ap.add_argument('--'+n,type=Path,required=True)
args=ap.parse_args();code=args.marker_code_dir.resolve();F=args.output.resolve();T=F/'tables';H=F/'h5';T.mkdir(parents=True,exist_ok=True);H.mkdir(exist_ok=True)
sys.path.insert(0,str(code));from marker_input_paths import resolve_h5_path,prepare_specimen_source,load_object_crosswalk,exact_left_join
markers=json.loads((code/'predeclared_marker_panel.json').read_text())['markers'];genes=sorted(set(markers+['CD3E','CD3G','CD247','TRGC1','TRGC2','IGJ'])-{'JCHAIN'});ix={g:i for i,g in enumerate(genes)};x=np.zeros((3,len(genes)),np.int32)
for g in ['CD3D','CD3E','CD3G','TRAC','TRDC','TRGC1','NKG7']:x[0,ix[g]]=2
for g in ['IGJ','MZB1','SDC1','CD79A','MKI67']:x[1,ix[g]]=3
x[2,ix['CD79A']]=1;x=sparse.csr_matrix(x);path=H/'synthetic_marker.h5'
with h5py.File(path,'w')as h:
 h.create_dataset('names_obs',data=np.array(['c1','c2','c3'],dtype=h5py.string_dtype()));h.create_dataset('names_var',data=np.array(genes,dtype=h5py.string_dtype()))
 # Three blocks include targets in different blocks and a final nontarget block.
 for block in range(3):
  row=x[block:block+1];g=h.create_group(f'assay/RNA/layers/rawdata/{block:04d}');g.attrs['shape']=row.shape;g.create_dataset('data',data=row.data);g.create_dataset('indices',data=row.indices);g.create_dataset('indptr',data=row.indptr)
pd.DataFrame({'OldCells':['c1','c2','c3'],'Sample':['anonymous']*3,'Dataset':['synthetic']*3,'Data_object':['synthetic_marker']*3}).to_csv(T/'all_release_original_metadata.tsv.gz',sep='\t',index=False,compression='gzip')
for lineage,cell,label in [('T_NK','c1','NKT_KLRC2'),('B_Plasma','c2','Plasmablast')]:
 pd.DataFrame({'OldCells':[cell],'Sample':['anonymous'],'Dataset':['synthetic'],'celltype':[label],'celltype1':[label],'celltype2':[label]}).to_csv(T/f'{lineage}_original_integration_metadata.tsv.gz',sep='\t',index=False,compression='gzip')
source=pd.DataFrame({'Data_object':['synthetic_marker'],'Sample':['anonymous'],'Dataset':['synthetic'],'Condition':['synthetic'],'Tissue':['synthetic']});source.to_csv(T/'composition_specimen_denominators_and_eligibility.tsv',sep='\t',index=False)
(F/'registry.json').write_text(json.dumps({'datasets':[{'h5_paths':{'RNA':'synthetic_marker.h5'}}]}));pd.DataFrame({'Data_object':['synthetic_marker'],'Corrected_Data_object':['synthetic_marker']}).to_csv(F/'crosswalk.tsv',sep='\t',index=False)
checks={'relative_H5_root_resolution':resolve_h5_path('synthetic_marker.h5',H)==path.resolve(),'Sample_only_preserved':prepare_specimen_source(source,{}).Sample.iloc[0]=='anonymous'}
with h5py.File(path,'r')as h:checks['multi_raw_block_fixture']=len(h['assay/RNA/layers/rawdata'])==3
for value,root,label in [('synthetic_marker.h5',None,'absent_root_rejected'),('../synthetic_marker.h5',H,'parent_traversal_rejected')]:
 try:resolve_h5_path(value,root);checks[label]=False
 except ValueError:checks[label]=True
same=source.copy();same['Local_sample_label']='anonymous';checks['matching_dual_sample_columns_preserved']=prepare_specimen_source(same,{}).Sample.iloc[0]=='anonymous';bad=same.copy();bad.Local_sample_label='disagrees'
try:prepare_specimen_source(bad,{});checks['disagreeing_sample_identity_rejected']=False
except ValueError:checks['disagreeing_sample_identity_rejected']=True
common=['--input-tables',str(T),'--registry',str(F/'registry.json'),'--data-h5-root',str(H),'--reader',str(args.reader.resolve()),'--crosswalk',str(F/'crosswalk.tsv')]
commands=[[sys.executable,str(code/'audit_full_RNA_marker_coherence.py'),*common,'--output',str(F/'full')],[sys.executable,str(code/'build_marker_review_evidence.py'),'--input-dir',str(F/'full')],[sys.executable,str(code/'audit_targeted_actual_joint_markers.py'),*common,'--marker-audit-dir',str(F/'full'),'--output',str(F/'joint')]]
for number,cmd in enumerate(commands,1):
 result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);(F/f'command{number}.log').write_text(result.stdout);checks[f'portable_cli_stage{number}_exit0']=result.returncode==0
 if result.returncode:raise RuntimeError(f'Synthetic stage {number} failed; see command log')
j=pd.read_csv(F/'joint/actual_joint_source_evidence.tsv',sep='\t');checks['actual_T_CD3D_TRDC_joint_detected']=j.loc[j.celltype2.eq('NKT_KLRC2'),'CD3D_TRDC_actual_joint_Fraction'].iloc[0]==1;checks['explicit_IGJ_fallback']=j.effective_JCHAIN_feature.eq('IGJ').all();checks['actual_plasma_joint_detected']=j.loc[j.celltype2.eq('Plasmablast'),'plasma_2of3_joint_Fraction'].iloc[0]==1
# The old pipe-concatenated identity could collide; composite keys cannot.
l=pd.DataFrame({'OldCells':['a|b','a'],'Sample':['c','b|c'],'Dataset':['d','d']})
r=l.copy();r['Data_object']=['first','second']
checks['composite_identity_delimiter_collision_avoided']=exact_left_join(l,r,on=['OldCells','Sample','Dataset'],validate='one_to_one',context='synthetic composite identity').Data_object.tolist()==['first','second']
for label,old,new in [('duplicate','x','y'),('conflicting','x','z'),('empty_original','','y'),('empty_corrected','x','')]:
 table=pd.DataFrame({'Data_object':['x',old],'Corrected_Data_object':['y',new]});p=F/f'crosswalk_{label}.tsv';table.to_csv(p,sep='\t',index=False)
 try:load_object_crosswalk(p);checks[f'crosswalk_{label}_rejected']=False
 except ValueError:checks[f'crosswalk_{label}_rejected']=True
for failure in ['unmatched_integrated_cell','unmatched_specimen_source']:
 bad_tables=F/f'{failure}_tables';shutil.copytree(T,bad_tables)
 if failure=='unmatched_integrated_cell':
  p=bad_tables/'all_release_original_metadata.tsv.gz';release=pd.read_csv(p,sep='\t',keep_default_na=False);release=release[release.OldCells.ne('c1')];release.to_csv(p,sep='\t',index=False,compression='gzip')
 else:
  wrong=source.copy();wrong.Sample='no_matching_specimen';wrong.to_csv(bad_tables/'composition_specimen_denominators_and_eligibility.tsv',sep='\t',index=False)
 for stage,original in [('full',commands[0]),('joint',commands[2])]:
  cmd=original.copy();cmd[cmd.index('--input-tables')+1]=str(bad_tables);cmd[cmd.index('--output')+1]=str(F/f'{failure}_{stage}')
  result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True);(F/f'{failure}_{stage}.log').write_text(result.stdout)
  checks[f'{stage}_{failure}_fail_closed']=result.returncode!=0 and 'incomplete exact join' in result.stdout
checks={k:bool(v)for k,v in checks.items()};assert all(checks.values()),checks
report={'passed':True,'checks':checks,'check_count':len(checks),'fixture':'Three explicitly synthetic cells; no project biological records','marker_code_dir':str(code),'reader':str(args.reader.resolve())};(F/'portable_synthetic_QA.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))

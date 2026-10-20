"""Pure synthetic installation checks; no source-study person/cell data are used."""
from pathlib import Path
import os,gzip,json,subprocess,shutil,sys
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('SCAID_SYNTHETIC_OUTPUT',str(ROOT/'synthetic_example/work'))).resolve();OUT.mkdir(parents=True,exist_ok=True)
if sys.byteorder!='little':raise RuntimeError('The example binary map format requires a little-endian host.')
compiler=shutil.which('c++') or shutil.which('g++')
if not compiler:raise RuntimeError('Install a C++17 compiler and zlib development library.')
exe=OUT/'stream_B_watch'
subprocess.run([compiler,'-O2','-std=c++17',str(ROOT/'scripts/stream_B_watch.cpp'),'-lz','-o',str(exe)],check=True)
entries=[(1,1,2),(2,1,3),(4,1,4),(1,2,1),(3,2,7),(2,3,5),(3,3,6),(4,3,1),(1,4,9),(4,5,2),(3,5,8),(1,5,1)]
mtx=OUT/'pure_synthetic.mtx.gz'
with gzip.open(mtx,'wt') as f:
 f.write('%%MatrixMarket matrix coordinate integer general\n% Entirely synthetic 4-gene, 5-cell example\n4 5 12\n')
 for row in entries:f.write('%d %d %d\n'%row)
np.array([-1,0,-1,1,-1,2],dtype='<i4').tofile(OUT/'columns.int32')
np.array([-1,0,-1,1,-1],dtype='<i4').tofile(OUT/'genes.int32')
pre=OUT/'toy'
subprocess.run([str(exe),str(mtx),str(OUT/'columns.int32'),str(OUT/'genes.int32'),str(pre)],check=True,capture_output=True,text=True)
libs=np.fromfile(str(pre)+'_library.uint64',dtype='<u8');watch=np.fromfile(str(pre)+'_watch.uint64',dtype='<u8').reshape(2,3)
assert np.array_equal(libs,[9,12,11])
assert np.array_equal(watch,[[2,0,1],[0,6,8]])
scan=json.loads(Path(str(pre)+'_scan.json').read_text());assert scan['entries_read']==12 and scan['selected_B_count_sum']==32
# A malformed count in a selected column must fail rather than silently enter a score.
bad=OUT/'invalid_selected_count.mtx.gz'
with gzip.open(bad,'wt') as f:f.write('%%MatrixMarket matrix coordinate integer general\n4 5 1\n1 1 -1\n')
failed=subprocess.run([str(exe),str(bad),str(OUT/'columns.int32'),str(OUT/'genes.int32'),str(OUT/'invalid')],capture_output=True,text=True)
assert failed.returncode!=0
# Independent 200-person toy regression; these are not anonymized study records.
rng=np.random.default_rng(20261020);n=200
frame=pd.DataFrame({'synthetic_id':['synthetic_person_%03d'%i for i in range(n)],'condition':pd.Categorical(['Control']*100+['Case']*100,categories=['Control','Case']),'cohort':np.tile(['2','3','4','2'],50),'Age10':rng.normal(0,1,n),'Sex':rng.choice(['Female','Male'],n),'pop_cov':rng.choice(['Synthetic_A','Synthetic_B'],n)})
case=frame.condition.eq('Case').to_numpy();outcomes=[]
for name,true_delta in [('synthetic_IFN_A',.17),('synthetic_IFN_B',.10),('synthetic_null',0.0)]:
 frame['y']=.2+true_delta*case+.025*frame.Age10+rng.normal(0,.07,n)
 fit=smf.ols('y ~ C(cohort) + Age10 + C(Sex) + C(pop_cov) + C(condition)',data=frame).fit(cov_type='HC3');key='C(condition)[T.Case]';design=fit.model.exog
 assert np.linalg.matrix_rank(design)==design.shape[1] and fit.df_resid>=10 and np.isfinite(fit.bse[key]) and fit.bse[key]>0
 ci=fit.conf_int().loc[key]
 outcomes.append({'outcome':name,'true_synthetic_difference':true_delta,'estimate':float(fit.params[key]),'HC3_95CI_low':float(ci.iloc[0]),'HC3_95CI_high':float(ci.iloc[1]),'p':float(fit.pvalues[key]),'people':n,'independent_rows_are_synthetic_people':True})
assert abs(outcomes[0]['estimate']-.17)<.04
q=multipletests([x['p'] for x in outcomes],method='fdr_bh')[1]
for row,value in zip(outcomes,q):row['BH_q']=float(value)
report={'passed':True,'purely_synthetic':True,'source_metadata_or_counts_used':False,'streamed_toy_counts':{'cells':5,'genes':4,'selected_cells':3,'library_totals':libs.tolist(),'watch_counts':watch.tolist(),'malformed_selected_count_rejected':True},'HC3_200_synthetic_people':outcomes,'not_a_reproduction_of_real_SLE_effects':True}
(OUT/'synthetic_check.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'passed':True,'synthetic_people':n,'tiny_matrix_streaming_passed':True,'output':str(OUT/'synthetic_check.json')}))

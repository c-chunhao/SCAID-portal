"""End-to-end public demo with only invented people, values and source pools.
Usage: python synthetic_pool_demo.py --output /tmp/scaid_pool_synthetic
"""
from pathlib import Path
import argparse,json,os,subprocess,sys
import numpy as np,pandas as pd,statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

def main(out):
 rng=np.random.default_rng(803);N=200
 baseline=out/'synthetic_baseline';bcell=out/'synthetic_bcell';result=out/'synthetic_results'
 for p in [baseline/'tables',bcell/'tables']:p.mkdir(parents=True,exist_ok=True)
 sex=np.tile(['Female','Male'],N//2);condition=np.tile(['HD','HD','SLE','SLE','SLE'],N//5)
 base=pd.DataFrame({'ind_cov':[f'SYNTHETIC_PERSON_{i:03d}' for i in range(N)],'batch_cov':np.repeat([f'SYNTHETIC_POOL_{i:02d}' for i in range(20)],10),'cohort':np.repeat(np.repeat([2,3,4],[6,7,7]),10),'condition':condition,'n_cells':1000,'Age10':rng.normal(size=N),'Sex':sex,'pop_cov':np.tile(['European','East_Asian'],N//2)})
 # Randomize ancestry independently of sex to keep the source-pool FE design full rank.
 base['pop_cov']=rng.choice(['European','East_Asian'],N);base['unit_id']=base.ind_cov+'|'+base.batch_cov;base['Age']=40+base.Age10*10;base['Status']='SYNTHETIC'
 f=rng.dirichlet(np.ones(8)*3,size=N);states=['T4','T8','B','NK','cM','ncM','cDC','pDC']
 for i,state in enumerate(states):base[state]=f[:,i]
 base.to_csv(baseline/'tables/GSE174188_unique_donor_composition_complete_zeros.tsv',sep='\t',index=False);base.to_csv(baseline/'tables/GSE174188_all_donor_batch_inventory.tsv',sep='\t',index=False)
 pool_effect=np.repeat(rng.normal(0,.05,20),10);sle=base.condition.eq('SLE').astype(int)
 score={s:.2+.15*sle+pool_effect+rng.normal(0,.05,N) for s in ['cM','T4','B']}
 scores=pd.concat([base[['ind_cov','cohort','condition']].assign(cg_cov=s,IFN12_log1p10k_mean=v,n_cells=80) for s,v in score.items()],ignore_index=True)
 scores.to_csv(baseline/'tables/GSE174188_target_cell_IFN_by_unique_donor.tsv',sep='\t',index=False)
 def reference(z):
  zz=z.copy();zz['cohort']=zz.cohort.astype(str);fit=smf.ols('y ~ C(cohort) + Age10 + C(Sex) + C(pop_cov) + C(condition)',zz).fit(cov_type='HC3');key='C(condition)[T.SLE]';ci=fit.conf_int().loc[key]
  return fit,dict(adjusted_coefficient_SLE_minus_HD=float(fit.params[key]),HC3_95CI_low=float(ci.iloc[0]),HC3_95CI_high=float(ci.iloc[1]),p_raw=float(fit.pvalues[key]))
 original=[]
 for kind,candidates in [('composition',states),('IFN12',['cM','T4','B'])]:
  for s in candidates:
   fit,ref=reference(base.assign(y=base[s] if kind=='composition' else score[s]));original.append(dict(analysis='primary_unique_batch',kind=kind,cell_type=s,**ref))
 old=pd.DataFrame(original);old['p_BH']=np.nan
 for _,ix in old.groupby('kind').groups.items():old.loc[ix,'p_BH']=multipletests(old.loc[ix,'p_raw'],method='fdr_bh')[1]
 old.to_csv(baseline/'tables/GSE174188_unique_donor_HC3_composition_IFN_tests.tsv',sep='\t',index=False)
 bp=base.copy();fractions=rng.dirichlet([8,5,1,1],N)
 for i,s in enumerate(['B_naive','B_mem','B_atypical']):bp['fraction_within_B__'+s]=fractions[:,i]
 bp.to_csv(bcell/'tables/B_individual_composition_and_IFN.tsv',sep='\t',index=False)
 models=[];rows=[]
 def add(analysis,family,state,endpoint,y,subset=None,extra=False):
  z=base.assign(y=y).copy()
  if subset is not None:z=z.loc[subset]
  zz=z.copy();zz['cohort']=zz.cohort.astype(str);formula='y ~ C(cohort) + Age10 + C(Sex) + C(pop_cov) + C(condition)'
  if extra:
   for s in ['B_naive','B_mem','B_atypical']:zz['fraction_within_B__'+s]=bp.loc[z.index,'fraction_within_B__'+s]
   formula+=''.join(' + Q("fraction_within_B__'+s+'")' for s in ['B_naive','B_mem','B_atypical'])
  if min(z.condition.value_counts())>=10:
   fit=smf.ols(formula,zz).fit(cov_type='HC3');key='C(condition)[T.SLE]';ci=fit.conf_int().loc[key];vals=dict(adjusted_SLE_minus_HD=float(fit.params[key]),HC3_CI_low=float(ci.iloc[0]),HC3_CI_high=float(ci.iloc[1]),BH_q=1.)
  else:vals=dict(adjusted_SLE_minus_HD=np.nan,HC3_CI_low=np.nan,HC3_CI_high=np.nan,BH_q=np.nan)
  models.append(dict(analysis=analysis,family=family,state=state,endpoint=endpoint,**vals));rr=z[['ind_cov','unit_id','condition','cohort','Age10','Sex','pop_cov','y']].assign(analysis=analysis,family=family,state=state,endpoint=endpoint,eligible_cell_count=80);rows.append(rr)
 add('same_eligible_B_population','matched_B_population','B','B_fraction_all_PBMC',base.B)
 add('same_eligible_B_population','matched_B_population','B','B_IFN12',score['B'])
 for s in ['B_naive','B_mem','B_atypical']:
  add('source_B_state_primary','state_composition_within_B',s,'state_fraction_within_B',bp['fraction_within_B__'+s])
  sub=list(base.index[base.condition.eq('SLE')][:30])+list(base.index[base.condition.eq('HD')][:3]) if s=='B_atypical' else None
  add('source_B_state_primary','state_IFN',s,'state_IFN12',score['B']+rng.normal(0,.02,N),subset=sub)
 for s in ['B_naive','B_mem']:
  add('common_state_eligible_population','common_state_population',s,'state_fraction_within_B',bp['fraction_within_B__'+s]);add('common_state_eligible_population','common_state_population',s,'state_IFN12',score['B']+rng.normal(0,.02,N))
 add('assessed_state_population_B_IFN_unadjusted','mixture_adjusted_IFN','B','B_IFN12',score['B'])
 add('state_fraction_adjusted_B_IFN','mixture_adjusted_IFN','B','B_IFN12',score['B'],extra=True)
 for s in ['missing_fine_label','other_labelled','total_unresolved']:add('source_label_coverage_diagnostic','source_label_coverage',s,'coverage_fraction_within_B',rng.uniform(0,.1,N))
 pd.DataFrame(models).to_csv(bcell/'tables/B_state_HC3_models.tsv',sep='\t',index=False);pd.concat(rows,ignore_index=True).to_csv(bcell/'tables/B_state_model_input_rows.tsv',sep='\t',index=False)
 env=os.environ.copy();env.update(SCAID_BASELINE_SCIENCE=str(baseline),SCAID_BCELL_OUTPUT=str(bcell),SCAID_POOL_OUTPUT=str(result),OPENBLAS_NUM_THREADS='1')
 here=Path(__file__).resolve().parent
 with (out/'synthetic_analysis.log').open('w') as handle:subprocess.run([sys.executable,str(here/'analyze_source_pool_sensitivity.py')],env=env,stdout=handle,check=True)
 subprocess.run([sys.executable,str(here/'independent_pool_numeric_QA.py'),'--directory',str(result)],env=env,check=True)
 qa=json.loads((result/'provenance/numeric_QA.json').read_text());assert qa['passed'] and qa['estimated']==54 and qa['n_model_rows']==56
 report={'passed':True,'only_invented_inputs':True,'N_synthetic':N,'synthetic_source_pools':20,'model_rows':56,'estimated':54,'unavailable_planned_tests':2}
 (out/'synthetic_demo_QA.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())

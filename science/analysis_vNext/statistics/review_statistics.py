"""Read-only statistical audit plus independent donor-level sensitivity analyses.
Never overwrite the submitted analysis or treat cells/visits as biological replicates.
"""
from pathlib import Path
import sys, os
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from runtime_paths import required_path, registered_datasets
import hashlib,json,platform
import numpy as np
import pandas as pd
from formulaic import model_matrix
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr

ROOT=required_path('SCAID_STATISTICS_AUDIT_DIR')
SRC=required_path('SCAID_SCIENCE_TABLES')
T=ROOT/'tables'; checks=[];inputs=[];design_rows=[];convergence_rows=[]
def read(name):
 p=SRC/name; inputs.append({'path':str(p),'SHA256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 return pd.read_csv(p,sep='\t')
for mode in ['primary_unique_batch','cohort4_sensitivity']:
 for ct in ['cM','T4','B']:
  prefix=f'GSE174188_{ct}_{mode}'
  c=read(prefix+'_counts.tsv').set_index('unit_id'); m=read(prefix+'_metadata.tsv').set_index('unit_id')
  r=read(prefix+'_DESeq2.tsv');q=read(prefix+'_fit_convergence.tsv');q=q.set_index(q.columns[0])
  assert c.index.equals(m.index) and m.ind_cov.is_unique and c.columns.is_unique
  assert np.isfinite(c.to_numpy()).all() and (c.to_numpy()>=0).all() and np.equal(c.to_numpy(),np.rint(c.to_numpy())).all()
  assert m.target_cells.ge(20).all() and m[['Age','Sex','pop_cov','condition']].notna().all().all()
  m['cohort']=m.cohort.astype(str);m['Age10']=(m.Age-40)/10
  design='~cohort + Age10 + Sex + pop_cov + condition' if m.cohort.nunique()>1 else '~Age10 + Sex + pop_cov + condition'
  dm=model_matrix(design,m);x=dm.to_numpy();assert np.linalg.matrix_rank(x)==x.shape[1]
  for col in ['cohort','Sex','pop_cov']:
   zz=pd.crosstab(m[col],m.condition)
   for level,row in zz.iterrows():design_rows.append({'model':prefix,'covariate':col,'level':level,'HD':int(row.get('HD',0)),'SLE':int(row.get('SLE',0))})
  keep=(c>=10).sum(axis=0)>=10;assert set(r.gene)==set(c.columns[keep])
  finite=r.pvalue.notna()&r.padj.notna();bh=multipletests(r.loc[finite,'pvalue'],method='fdr_bh')[1]
  maxerr=float(np.max(np.abs(bh-r.loc[finite,'padj']))) if finite.any() else 0
  assert maxerr<1e-12 and r.gene.is_unique
  invalid=set(q.index[~q['_LFC_converged'].fillna(False)])
  assert not r.loc[r.gene.isin(invalid),'pvalue'].notna().any()
  for flag in q:
   for gene in q.index[q[flag].eq(False)]:convergence_rows.append({'model':prefix,'gene':gene,'flag':flag,'interpretation':'excluded from inference' if flag=='_LFC_converged' else 'dispersion optimizer fallback to grid search; not automatically an invalid fit'})
  checks.append({'model':prefix,'donors':len(m),'HD':int(m.condition.eq('HD').sum()),'SLE':int(m.condition.eq('SLE').sum()),'rank':int(np.linalg.matrix_rank(x)),'parameters':x.shape[1],'scaled_condition_number':float(np.linalg.cond(x)),'n_genes_tested':len(r),'n_valid_adjusted_p':int(finite.sum()),'BH_max_abs_error':maxerr,'integer_count_matrix':True,'unique_donors':True,'metadata_aligned':True,'min_target_cells':int(m.target_cells.min()),'library_min':int(c.sum(axis=1).min()),'library_max':int(c.sum(axis=1).max()),'n_dispersion_genewise_grid_fallback':int(q['_genewise_converged'].eq(False).sum()),'n_dispersion_MAP_grid_fallback':int(q['_MAP_converged'].eq(False).sum()),'n_LFC_nonconvergent_excluded':len(invalid)})
pd.DataFrame(design_rows).to_csv(T/'donor_model_covariate_overlap.tsv',sep='\t',index=False)
pd.DataFrame(convergence_rows).to_csv(T/'donor_DE_convergence_semantics.tsv',sep='\t',index=False)
comp=read('GSE174188_unique_donor_composition_complete_zeros.tsv'); original=read('GSE174188_unique_donor_HC3_composition_IFN_tests.tsv');scores=read('GSE174188_target_cell_IFN_by_unique_donor.tsv')
allct=['B','NK','PB','Progen','Prolif','T4','T8','cDC','cM','ncM','pDC'];tested=['T4','T8','B','NK','cM','ncM','cDC','pDC']
assert comp.ind_cov.is_unique and np.allclose(comp[allct].sum(axis=1),1)
assert comp[allct].ge(0).all().all()
covform='C(cohort) + Age10 + C(Sex) + C(pop_cov) + C(condition)'
rows=[]
def fit(data,y,name,ct,family):
 data=data.copy();data['y']=np.asarray(y);data['cohort']=data.cohort.astype(str)
 cov=covform if data.cohort.nunique()>1 else 'Age10 + C(Sex) + C(pop_cov) + C(condition)'
 model=smf.ols('y ~ '+cov,data=data).fit(cov_type='HC3');key='C(condition)[T.SLE]'
 assert np.linalg.matrix_rank(model.model.exog)==model.model.exog.shape[1]
 ci=model.conf_int().loc[key]
 rows.append({'analysis':name,'family':family,'cell_type':ct,'n_donors':len(data),'HD':int(data.condition.eq('HD').sum()),'SLE':int(data.condition.eq('SLE').sum()),'coefficient':float(model.params[key]),'HC3_CI_low':float(ci.iloc[0]),'HC3_CI_high':float(ci.iloc[1]),'pvalue':float(model.pvalues[key]),'max_hat_leverage':float(model.get_influence().hat_matrix_diag.max()),'effect_scale':name})
 return model
for ct in tested:
 r=fit(comp,comp[ct],'original_fraction_replication',ct,'composition_8')
 old=original[(original.analysis=='primary_unique_batch')&original.kind.eq('composition')&original.cell_type.eq(ct)].iloc[0]
 assert np.isclose(r.params['C(condition)[T.SLE]'],old.adjusted_coefficient_SLE_minus_HD,rtol=1e-10,atol=1e-12)
 # Stabilized transforms acknowledge bounded/relative fractions; they do not
 # recover absolute immune-cell abundances from compositional sequencing.
 counts=comp[allct].mul(comp.n_cells,axis=0);assert np.allclose(counts,np.rint(counts))
 counts=np.rint(counts)
 stab=(counts+.5).div(comp.n_cells+.5*len(allct),axis=0)
 clr=np.log(stab).sub(np.log(stab).mean(axis=1),axis=0)
 fit(comp,np.arcsin(np.sqrt(comp[ct])),'arcsine_fraction_sensitivity',ct,'composition_8')
 fit(comp,clr[ct],'CLR_relative_abundance_sensitivity',ct,'composition_8')
for threshold in [20,50,100]:
 for ct in ['cM','T4','B']:
  z=scores[scores.cg_cov.eq(ct)&scores.n_cells.ge(threshold)].merge(comp[['ind_cov','Age10','Sex','pop_cov']],on='ind_cov',validate='many_to_one')
  assert z.ind_cov.is_unique
  rr=fit(z,z.IFN12_log1p10k_mean,f'IFN12_min_{threshold}_cells',ct,'IFN12_3')
  if threshold==20:
   old=original[(original.analysis=='primary_unique_batch')&original.kind.eq('IFN12')&original.cell_type.eq(ct)].iloc[0]
   assert np.isclose(rr.params['C(condition)[T.SLE]'],old.adjusted_coefficient_SLE_minus_HD,rtol=1e-10)
for name,sub in [('female_only_design_sensitivity',comp[comp.Sex.eq('Female')]),('ancestry_overlap_only_sensitivity',comp[comp.pop_cov.isin(['Asian','European'])])]:
 for ct in tested:fit(sub,sub[ct],name,ct,'composition_8')
 for ct in ['cM','T4','B']:
  z=scores[scores.cg_cov.eq(ct)&scores.n_cells.ge(20)&scores.ind_cov.isin(sub.ind_cov)].merge(sub[['ind_cov','Age10','Sex','pop_cov']],on='ind_cov',validate='many_to_one')
  fit(z,z.IFN12_log1p10k_mean,name,ct,'IFN12_3')
out=pd.DataFrame(rows);out['p_BH']=np.nan
for ix in out.groupby(['analysis','family']).groups.values():out.loc[ix,'p_BH']=multipletests(out.loc[ix,'pvalue'],method='fdr_bh')[1]
out.to_csv(T/'donor_HC3_independent_replication_and_sensitivities.tsv',sep='\t',index=False)
# Stratified resampling of donors, not cells, conditions preserved within cohorts.
rng=np.random.default_rng(20261005);d=comp.copy();d['cohort']=d.cohort.astype(str)
x=model_matrix('~'+covform,d);idx=x.columns.get_loc('C(condition)[T.SLE]');Y=d[tested].to_numpy();strata=[np.asarray(v)for v in d.groupby(['cohort','condition']).indices.values()]
boots=[]
for b in range(1000):
 ix=np.concatenate([rng.choice(v,size=len(v),replace=True)for v in strata]);xx=x.to_numpy()[ix]
 if np.linalg.matrix_rank(xx)<xx.shape[1]:continue
 boots.append(np.linalg.lstsq(xx,Y[ix],rcond=None)[0][idx])
boot=np.asarray(boots);pd.DataFrame({'cell_type':tested,'n_valid_bootstraps':len(boot),'donor_bootstrap_CI_low':np.quantile(boot,.025,axis=0),'donor_bootstrap_CI_high':np.quantile(boot,.975,axis=0),'stratification':'processing cohort and condition; resampling whole donors'}).to_csv(T/'composition_donor_stratified_bootstrap.tsv',sep='\t',index=False)
# Wild donor bootstrap preserves the covariate design, including rare strata;
# HC2-scaled donor residuals account for leverage without treating cells as rows.
xx=x.to_numpy();beta=np.linalg.lstsq(xx,Y,rcond=None)[0];mu=xx@beta;hat=np.einsum('ij,ji->i',xx,np.linalg.pinv(xx));resid=(Y-mu)/np.sqrt(1-hat[:,None]);pinv=np.linalg.pinv(xx)
wild=[]
for b in range(1000):
 weights=rng.choice([-1.,1.],size=len(comp));wild.append((pinv@(mu+resid*weights[:,None]))[idx])
wild=np.asarray(wild);pd.DataFrame({'cell_type':tested,'n_bootstraps':1000,'wild_donor_bootstrap_CI_low':np.quantile(wild,.025,axis=0),'wild_donor_bootstrap_CI_high':np.quantile(wild,.975,axis=0),'method':'fixed design Rademacher wild bootstrap of HC2-scaled donor residuals'}).to_csv(T/'composition_wild_donor_bootstrap.tsv',sep='\t',index=False)
# Primary and cohort-4 are overlapping sensitivity estimates, not replication.
cons=[]
for ct in ['cM','T4','B']:
 a=read(f'GSE174188_{ct}_primary_unique_batch_DESeq2.tsv');b=read(f'GSE174188_{ct}_cohort4_sensitivity_DESeq2.tsv')
 ab=a.merge(b,on='gene',suffixes=('_primary','_cohort4'));valid=ab.padj_primary.notna()&ab.padj_cohort4.notna();ab=ab[valid]
 cons.append({'cell_type':ct,'common_valid_genes':len(ab),'LFC_spearman':float(spearmanr(ab.log2FoldChange_primary,ab.log2FoldChange_cohort4).statistic),'sign_agreement_all':float(np.mean(np.sign(ab.log2FoldChange_primary)==np.sign(ab.log2FoldChange_cohort4))),'primary_FDR05_and_absLFC1':int((ab.padj_primary.lt(.05)&ab.log2FoldChange_primary.abs().ge(1)).sum()),'scope':'overlapping donors, sensitivity only'})
pd.DataFrame(cons).to_csv(T/'primary_cohort4_DE_sensitivity_concordance.tsv',sep='\t',index=False)
(ROOT/'statistical_QA.json').write_text(json.dumps({'passed':True,'models':checks,'composition_original_8_exactly_replicated':True,'IFN12_original_3_exactly_replicated':True,'composition_closure_complete_11_categories':True,'sensitivity_rows':len(out),'donor_bootstrap_seed':20261005,'n_valid_bootstraps':len(boot),'input_manifest':inputs,'software_python':platform.python_version(),'limitations':['SLE/reference analysis remains within a single source study.','Treatment/activity and residual technical confounding are not resolved by transforms or bootstrap.','CLR effects concern relative composition; no absolute abundance claim.','IFN12 is a fixed illustrative interferon program, not an unbiased search for disease-specific novelty.']},indent=2)+'\n')
print(out.to_string(index=False));print('Statistical audit passed.')

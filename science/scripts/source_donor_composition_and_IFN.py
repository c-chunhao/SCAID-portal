import os
"""One preselected source sequencing batch per source individual; complete zeros.
Within-study covariate-adjusted composition/IFN models with donor HC3 variance.
Excluding source activity-switch donors and cohort4 are sensitivity checks, not
independent replication, because some individuals overlap the primary analysis.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42})
source=Path(os.environ['SCAID_SLE_SOURCE_METADATA'])
m=pd.read_csv(source,usecols=['ind_cov','batch_cov','Processing_Cohort','cg_cov','SLE_status']);m['key']=m.ind_cov+'|'+m.batch_cov
frozen=pd.read_csv(T/'GSE174188_frozen_one_batch_per_donor.tsv',sep='\t');frozen['Age10']=(frozen.Age-40)/10
assert frozen.ind_cov.is_unique
inv=pd.read_csv(T/'GSE174188_all_donor_batch_inventory.tsv',sep='\t');activity=inv.groupby('ind_cov').Status.nunique();switch=set(activity[activity>1].index)
m=m[m.key.isin(frozen.unit_id)].copy();den=m.groupby('key').size();ct=pd.crosstab(m.key,m.cg_cov)
ct=ct.reindex(index=frozen.unit_id,fill_value=0);frac=ct.div(den.reindex(ct.index),axis=0)
comp=frozen.set_index('unit_id').join(frac);comp.to_csv(T/'GSE174188_unique_donor_composition_complete_zeros.tsv',sep='\t')
scores=pd.read_csv(T/'GSE174188_target_cell_IFN_by_unique_donor.tsv',sep='\t')
rows=[]
def fit(data,outcome,mode,kind,cell_type):
 assert data.ind_cov.is_unique
 data=data.copy();data['cohort']=data.cohort.astype(str)
 model='y ~ C(cohort) + Age10 + C(Sex) + C(pop_cov) + C(condition)'if data.cohort.nunique()>1 else'y ~ Age10 + C(Sex) + C(pop_cov) + C(condition)'
 data['y']=data[outcome];r=smf.ols(model,data=data).fit(cov_type='HC3');key='C(condition)[T.SLE]';assert key in r.params
 dm=r.model.exog;assert np.linalg.matrix_rank(dm)==dm.shape[1]
 ci=r.conf_int().loc[key]
 return {'analysis':mode,'kind':kind,'cell_type':cell_type,'n_HD':int(data.condition.eq('HD').sum()),'n_SLE':int(data.condition.eq('SLE').sum()),'median_HD':float(data.loc[data.condition.eq('HD'),'y'].median()),'median_SLE':float(data.loc[data.condition.eq('SLE'),'y'].median()),'adjusted_coefficient_SLE_minus_HD':float(r.params[key]),'HC3_95CI_low':float(ci.iloc[0]),'HC3_95CI_high':float(ci.iloc[1]),'p_raw':float(r.pvalues[key]),'design':model,'disease_parameter_estimable':True,'covariance':'HC3, source individuals are the independent rows'}
for mode,sub in [('primary_unique_batch',comp),('activity_switch_excluded',comp[~comp.ind_cov.isin(switch)]),('selected_batch_cohort4',comp[comp.cohort.eq(4)])]:
 for ct in ['T4','T8','B','NK','cM','ncM','cDC','pDC']:
  rows.append(fit(sub,ct,mode,'composition',ct))
 for ct in ['cM','T4','B']:
  s=scores[scores.cg_cov.eq(ct)&scores.n_cells.ge(20)].merge(sub.reset_index()[['ind_cov','Age10','Sex','pop_cov']],on='ind_cov',validate='many_to_one')
  s=s[s.ind_cov.isin(sub.ind_cov)];rows.append(fit(s,'IFN12_log1p10k_mean',mode,'IFN12',ct))
r=pd.DataFrame(rows);r['p_BH']=np.nan
for ix in r.groupby(['analysis','kind']).groups.values():r.loc[ix,'p_BH']=multipletests(r.loc[ix,'p_raw'],method='fdr_bh')[1]
r.to_csv(T/'GSE174188_unique_donor_HC3_composition_IFN_tests.tsv',sep='\t',index=False)
# An honestly captioned within-study result, with primary and explicitly overlapping
# sensitivity effects. No cells used as independent replicates.
fig,axes=plt.subplots(1,2,figsize=(6.5,4.4),gridspec_kw={'width_ratios':[1.1,1]})
a=r[(r.analysis=='primary_unique_batch')&r.kind.eq('composition')].copy();yy=np.arange(len(a));effect=a.adjusted_coefficient_SLE_minus_HD*100;lo=a.HC3_95CI_low*100;hi=a.HC3_95CI_high*100
axes[0].errorbar(effect,yy,xerr=np.vstack([effect-lo,hi-effect]),fmt='o',color='#587f8a',capsize=2,markersize=4);axes[0].set_yticks(yy,a.cell_type);axes[0].invert_yaxis();axes[0].axvline(0,color='#999999',linewidth=.7);axes[0].set_xlabel('Adjusted SLE − HD\nfraction (percentage points)');axes[0].set_title('A  Source cell composition',loc='left')
b=r[(r.analysis=='primary_unique_batch')&r.kind.eq('IFN12')].copy();yy=np.arange(len(b));effect=b.adjusted_coefficient_SLE_minus_HD;lo=b.HC3_95CI_low;hi=b.HC3_95CI_high
axes[1].errorbar(effect,yy,xerr=np.vstack([effect-lo,hi-effect]),fmt='o',color='#9a503b',capsize=2,markersize=4);axes[1].set_yticks(yy,b.cell_type);axes[1].invert_yaxis();axes[1].axvline(0,color='#999999',linewidth=.7);axes[1].set_xlabel('Adjusted SLE − HD\nmean IFN12 log1p(10k)');axes[1].set_title('B  Cell-type IFN12 score',loc='left')
for ax in axes:ax.spines[['top','right']].set_visible(False)
fig.text(.5,.02,'GSE174188: one source batch per donor; age, sex, ancestry and cohort adjusted; HC3 95% CI.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.055,1,1]);fig.savefig(F/'Verified_donor_GSE174188.pdf');fig.savefig(F/'Verified_donor_GSE174188.png',dpi=600);fig.savefig(F/'Verified_donor_GSE174188.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
s={'source_case_control_conflicting_individuals':0,'source_activity_switch_individuals':sorted(switch),'activity_switch_exclusion_reason':'Disease activity/treatment Status varies across source batches; condition SLE/Healthy never conflicts. Exclusion is a sensitivity analysis.','primary_unique_donors':len(comp),'primary_HD':int(comp.condition.eq('HD').sum()),'primary_SLE':int(comp.condition.eq('SLE').sum()),'composition_denominator':'All source-study processed PBMC cells in the single preselected sequencing batch for each individual; source cg_cov assigned in the same vocabulary to both arms.','multiple_testing':'BH within each analysis and outcome family: 8 composition lineages and 3 IFN12 lineage scores. Source cells are not replicates.','limitations':['Observational case/control association; residual batch, clinical activity and treatment confounding can remain.','Source processing cohorts and activity-switch exclusion overlap the primary donor set and are not independent replication.','Within-study source-reference analysis does not claim a cross-disease disease-specific signal.'],'results':r.to_dict(orient='records')}
(B/'donor_HC3_analysis_summary.json').write_text(json.dumps(s,indent=2)+'\n');print(r[['analysis','kind','cell_type','n_HD','n_SLE','adjusted_coefficient_SLE_minus_HD','p_BH']].to_string(index=False))

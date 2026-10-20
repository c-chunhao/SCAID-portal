from pathlib import Path
import os
import pandas as pd,numpy as np,json,sys
from scipy.spatial.distance import cdist
W=Path(os.environ.get('SCAID_DIAGNOSTIC_DIR',str(Path(__file__).resolve().parent)))
A=sys.argv[1]
m=pd.read_csv(W/f'{A}_diagnostic_subset_metadata.tsv.gz',sep='\t').set_index('integration_cell_id')
s=pd.read_csv(W/f'{A}_diagnostic_sampling.tsv.gz',sep='\t')
r=pd.read_csv(W/f'{A}_source_mixing_diagnostics.tsv',sep='\t')
checks=[]
for rep in ['harmony','pca']:
 z=pd.read_csv(W/f'{A}_diagnostic_subset_{rep}.tsv.gz',sep='\t').set_index('integration_cell_id')
 assert z.index.equals(m.index) and np.isfinite(z.values).all()
 rows=r[r.representation.eq(rep)].sort_values(['sampled_cells','condition','tissue','major'])
 # Two full independent distance recomputations per representation; include a small and a large stratum.
 for n in [0,-1]:
  row=rows.iloc[n];key=json.dumps([row.condition,row.tissue,row.major])
  ids=s.loc[s.review.eq('mixing')&s.group.eq(key),'integration_cell_id']
  mm=m.loc[ids];xx=z.loc[ids].to_numpy(dtype=float);labels=mm.Dataset.to_numpy();N=len(xx)
  assert N==row.sampled_cells and (mm.condition.eq(row.condition)&mm.tissue.eq(row.tissue)&mm.major.eq(row.major)).all()
  assert mm.Dataset.value_counts().min()>=50 and mm.Dataset.value_counts().max()<=1000
  same=0
  for a in range(0,N,128):
   d=cdist(xx[a:a+128],xx,metric='sqeuclidean');j=np.arange(min(128,N-a));d[j,a+j]=np.inf
   ni=np.argsort(d,axis=1,kind='stable')[:,:30]
   same+=(labels[ni]==labels[a:a+128,None]).sum()
  observed=same/(N*30);delta=abs(observed-row.mean_same_source_fraction)
  # Float32 neighbour backend ties can change a boundary neighbour; allow <1 neighbour per thousand cells.
  assert delta<1e-4,(A,rep,key,observed,row.mean_same_source_fraction,delta)
  cc=mm.Dataset.value_counts().to_numpy();expected=float((cc*(cc-1)).sum()/(N*(N-1)))
  assert abs(expected-row.expected_same_source_fraction)<1e-12
  checks.append({'representation':rep,'stratum':[row.condition,row.tissue,row.major],'sampled_cells':N,'independent_cdist_same_source_fraction':observed,'reported_same_source_fraction':row.mean_same_source_fraction,'absolute_difference':delta,'finite_sample_expectation_verified':True})
x={'passed':True,'atlas':A,'independent_backend':'scipy.spatial.distance.cdist plus full per-query sorting; main analysis used sklearn nearest neighbours','checks':checks,'no_source_or_label_changes':True}
(W/f'{A}_independent_diagnostic_QA.json').write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(x,ensure_ascii=False,indent=2))

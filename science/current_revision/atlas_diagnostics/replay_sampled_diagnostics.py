"""Reproduce the reported mixing summaries using only supplied sampled coordinates."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
from pathlib import Path
import json,sys
import numpy as np,pandas as pd
from sklearn.neighbors import NearestNeighbors
W=Path(os.environ.get('SCAID_DIAGNOSTIC_DIR',str(Path(__file__).resolve().parent)))
A=sys.argv[1]
m=pd.read_csv(W/f'{A}_diagnostic_subset_metadata.tsv.gz',sep='\t').set_index('integration_cell_id')
s=pd.read_csv(W/f'{A}_diagnostic_sampling.tsv.gz',sep='\t')
r=pd.read_csv(W/f'{A}_source_mixing_diagnostics.tsv',sep='\t')
records=[]
for rep in ['harmony','pca']:
 z=pd.read_csv(W/f'{A}_diagnostic_subset_{rep}.tsv.gz',sep='\t').set_index('integration_cell_id')
 assert z.index.equals(m.index)
 for key,g in s[s.review.eq('mixing')].groupby('group',sort=True):
  condition,tissue,major=json.loads(key);ids=g.integration_cell_id.to_numpy();labels=m.loc[ids,'Dataset'].to_numpy();x=z.loc[ids].to_numpy(dtype=np.float32);N=len(x)
  ni=NearestNeighbors(n_neighbors=31,algorithm='brute',n_jobs=4).fit(x).kneighbors(return_distance=False)[:,:30]
  same=float((labels[ni]==labels[:,None]).mean());cc=pd.Series(labels).value_counts();expect=float((cc*(cc-1)).sum()/(N*(N-1)))
  old=r[r.representation.eq(rep)&r.condition.eq(condition)&r.tissue.eq(tissue)&r.major.eq(major)].iloc[0]
  assert abs(same-old.mean_same_source_fraction)<1e-7
  assert abs(expect-old.expected_same_source_fraction)<1e-12
  records.append({'representation':rep,'condition':condition,'tissue':tissue,'major':major,'sampled_cells':N,'same_source_fraction':same,'expected_same_source_fraction':expect})
qa={'passed':True,'atlas':A,'replayed_strata':len(records),'only_release_sampled_coordinates_used':True,'all_reported_mixing_values_match':True}
(W/f'{A}_sampled_replay_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
print(json.dumps(qa,indent=2))

import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','4');os.environ.setdefault('OMP_NUM_THREADS','4')
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,sys
from sklearn.neighbors import NearestNeighbors
from collections import Counter
import sklearn,scipy
from openpyxl import load_workbook
W=Path(os.environ.get('SCAID_DIAGNOSTIC_DIR',str(Path(__file__).resolve().parent)))
BASE=Path(os.environ.get('SCAID_REFERENCE_METADATA_DIR','reference_metadata'))
A=sys.argv[1];rng=np.random.default_rng(20261020)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for z in iter(lambda:f.read(16*1024*1024),b''):h.update(z)
 return h.hexdigest()
def writej(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
ins=json.loads((W/f'{A}_final_object_inspection.json').read_text())
m=pd.read_csv(W/f'{A}_final_metadata.tsv.gz',sep='\t',dtype=str)
original=pd.read_csv(BASE/f'{A}_original_integration_metadata.tsv.gz',sep='\t',dtype=str)
assert m.shape==original.shape,(m.shape,original.shape)
assert (m.fillna('')==original.fillna('')).all().all()
assert not m.integration_cell_id.duplicated().any()
m['condition']=m.Disease.replace({'SS':'SSc','T1DM':'T1D','AH':'AIH'})
m.loc[m.Dataset.eq('E-MTAB-8207'),'condition']='PsA'
m.loc[m.Dataset.eq('GSE198616'),'condition']='BD'
# Match exact source-specimen keys to current source tissue names; do not pool ileum, colon and rectum.
s4=Path(os.environ.get('SCAID_SUPPLEMENTARY_S4','Supplementary_Table_S4_Composition_Annotation_Validation.xlsx'))
w4=load_workbook(s4,read_only=True,data_only=True)
spec=list(w4['Specimen_denominators'].values)[1:]
lookup={}
for row in spec:
 key=(str(row[2]),str(row[1])); value=(str(row[3]),str(row[4]),str(row[0]))
 assert key not in lookup or lookup[key]==value
 lookup[key]=value
keys=list(zip(m.Dataset,m.Sample));assert all(key in lookup for key in keys)
m['condition']=[lookup[key][0] for key in keys]
m['tissue']=[lookup[key][1] for key in keys]
m['source_object']=[lookup[key][2] for key in keys]
m['major']=m.celltype
m['intermediate']=m.celltype1;m['fine']=m.celltype2
eligible=[];excluded=[];mix_samples={}
for (condition,tissue,major),g in m.groupby(['condition','tissue','major'],sort=True):
 counts=g.Dataset.value_counts();sources=counts[counts>=50].index.tolist()
 rec={'atlas':A,'condition':condition,'tissue':tissue,'major':major,'total_cells':len(g),'all_sources':len(counts),'eligible_sources_ge50':len(sources),'eligible_cells':int(counts.loc[sources].sum())}
 if len(sources)<2:excluded.append(rec);continue
 idx=[]
 for src in sorted(sources):
  ix=g.index[g.Dataset.eq(src)].to_numpy();idx.extend(rng.choice(ix,min(len(ix),1000),replace=False).tolist())
 idx=np.asarray(idx,dtype=int);rec['sampled_cells']=len(idx)
 key=(condition,tissue,major);mix_samples[key]=idx;eligible.append(rec)
# Biological-label review: query only major groups represented in >=2 eligible sources.
bio_samples={}
for (condition,tissue),g in m.groupby(['condition','tissue'],sort=True):
 cc=g.groupby(['major','Dataset']).size();valid=cc[cc>=50]
 majors=[major for major in sorted(g.major.unique()) if sum(1 for a,s in valid.index if a==major)>=2]
 if len(majors)<2:continue
 idx=[]
 for major in majors:
  for src in sorted(s for a,s in valid.index if a==major):
   ix=g.index[g.major.eq(major)&g.Dataset.eq(src)].to_numpy();idx.extend(rng.choice(ix,min(len(ix),500),replace=False).tolist())
 bio_samples[(condition,tissue)]=np.asarray(idx,dtype=int)

sampling=[]
for key,idx in mix_samples.items():
 for position in idx:sampling.append({'review':'mixing','group':json.dumps(key),'original_index':int(position),'integration_cell_id':m.loc[position,'integration_cell_id']})
for key,idx in bio_samples.items():
 for position in idx:sampling.append({'review':'labels','group':json.dumps(key),'original_index':int(position),'integration_cell_id':m.loc[position,'integration_cell_id']})
pd.DataFrame(sampling).to_csv(W/f'{A}_diagnostic_sampling.tsv.gz',sep='\t',index=False)
selected=np.sort(np.unique([item['original_index'] for item in sampling]))
m.loc[selected].to_csv(W/f'{A}_diagnostic_subset_metadata.tsv.gz',sep='\t',index=False)
mix_records=[];bio_records=[];rep_checks={}
for representation in ['harmony','pca']:
 p=W/f'{A}_final_{representation}.tsv.gz'
 if not p.exists():continue
 # Read only dimensions actually used by the saved final global display.
 first=pd.read_csv(p,sep='\t',nrows=0).columns
 z=pd.read_csv(p,sep='\t',usecols=list(first[:21]))
 assert np.array_equal(z.iloc[:,0].to_numpy(),m.integration_cell_id.to_numpy())
 x=z.iloc[:,1:].to_numpy(dtype=np.float32);assert x.shape==(len(m),20);assert np.isfinite(x).all()
 out=pd.DataFrame(x[selected],columns=[f'{representation}_{i}' for i in range(1,21)]);out.insert(0,'integration_cell_id',m.loc[selected,'integration_cell_id'].to_numpy());out.to_csv(W/f'{A}_diagnostic_subset_{representation}.tsv.gz',sep='\t',index=False)
 rep_checks[representation]={'dimensions_used':20,'saved_dimensions':ins['reductions'][representation]['dimensions'],'cells':len(m),'identity_order_matches_original_metadata':True,'path':str(p),'sha256':sha(p)}
 for rec in eligible:
  key=(rec['condition'],rec['tissue'],rec['major']);idx=mix_samples[key];labels=m.Dataset.to_numpy()[idx];n=len(idx)
  k=min(30,n-1);knn=NearestNeighbors(n_neighbors=k+1,algorithm='brute',n_jobs=4).fit(x[idx]);neigh=knn.kneighbors(return_distance=False)[:,:k]
  # sklearn excludes self when querying the fitted training set with X=None.
  assert not np.any(neigh==np.arange(n)[:,None])
  obs=(labels[neigh]==labels[:,None]).mean(axis=1)
  counts=Counter(labels);expect=np.array([(counts[s]-1)/(n-1) for s in labels])
  entropy=[]
  prop=np.array(list(counts.values()),float)/n;global_H=-np.sum(prop*np.log(prop))
  for row in labels[neigh]:
   pp=np.array(list(Counter(row).values()),float)/k;entropy.append(-np.sum(pp*np.log(pp))/global_H)
  mix_records.append({**rec,'representation':representation,'dimensions':20,'neighbors':k,'mean_same_source_fraction':float(obs.mean()),'expected_same_source_fraction':float(expect.mean()),'same_source_excess':float((obs-expect).mean()),'mean_source_entropy_ratio':float(np.mean(entropy)),'sampled_source_counts':json.dumps(dict(sorted(counts.items())))})
 for (condition,tissue),idx in bio_samples.items():
  mm=m.loc[idx];sources=mm.Dataset.to_numpy()
  for src in sorted(set(sources)):
   qpos=np.flatnonzero(sources==src);tpos=np.flatnonzero(sources!=src)
   k=min(30,len(tpos));nn=NearestNeighbors(n_neighbors=k,algorithm='brute',n_jobs=4).fit(x[idx[tpos]])
   neigh=tpos[nn.kneighbors(x[idx[qpos]],return_distance=False)]
   rec={'atlas':A,'condition':condition,'tissue':tissue,'representation':representation,'query_source':src,'query_cells':len(qpos),'other_source_cells':len(tpos),'other_sources':len(set(sources[tpos])),'dimensions':20,'neighbors':k}
   for level in ['major','intermediate','fine']:
    labels=mm[level].to_numpy();same=(labels[neigh]==labels[qpos,None]).mean(axis=1)
    tc=Counter(labels[tpos]);expected=np.array([tc[label]/len(tpos) for label in labels[qpos]])
    rec[f'mean_cross_source_same_{level}_fraction']=float(same.mean());rec[f'expected_same_{level}_fraction']=float(expected.mean())
   bio_records.append(rec)
 pd.DataFrame(mix_records).to_csv(W/f'{A}_source_mixing_diagnostics.tsv',sep='\t',index=False)
 pd.DataFrame(bio_records).to_csv(W/f'{A}_cross_source_label_diagnostics.tsv',sep='\t',index=False)
 print(A,representation,'done',len(mix_records),len(bio_records),flush=True)
 del x,z
pd.DataFrame(excluded).to_csv(W/f'{A}_source_mixing_ineligible_strata.tsv',sep='\t',index=False)
mr=pd.DataFrame(mix_records);br=pd.DataFrame(bio_records);summary=[]
for rep in mr.representation.unique():
 r=mr[mr.representation.eq(rep)];b=br[br.representation.eq(rep)]
 s={'atlas':A,'representation':rep,'eligible_mixing_strata':len(r),'eligible_final_atlas_cells':int(r.eligible_cells.sum()),'sampled_mixing_cells':int(r.sampled_cells.sum()),'all_strata':len(eligible)+len(excluded),'excluded_strata':len(excluded),'median_same_source_fraction':float(r.mean_same_source_fraction.median()),'median_expected_same_source_fraction':float(r.expected_same_source_fraction.median()),'median_same_source_excess':float(r.same_source_excess.median()),'median_source_entropy_ratio':float(r.mean_source_entropy_ratio.median()),'source_query_groups_for_label_review':len(b),'label_review_sampled_cells':int(b.query_cells.sum())}
 for level in ['major','intermediate','fine']:
  if len(b):s[f'cell_weighted_cross_source_same_{level}_fraction']=float(np.average(b[f'mean_cross_source_same_{level}_fraction'],weights=b.query_cells));s[f'cell_weighted_expected_same_{level}_fraction']=float(np.average(b[f'expected_same_{level}_fraction'],weights=b.query_cells))
 summary.append(s)
writej(W/f'{A}_final_diagnostic_summary.json',{'passed':True,'atlas':A,'saved_final_metadata_identical_to_original_reference':True,'final_object_cells':len(m),'representations':rep_checks,'summary':summary,'versions':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scikit_learn':sklearn.__version__,'scipy':scipy.__version__},'interpretation':'Descriptive conditional neighbourhood diagnostics of the saved final atlas. Stored-label concordance is not independent annotation accuracy; source mixing is not proof of removal of disease/tissue confounding. PCA is an integrated pre-Harmony representation, not an uncorrected RNA control.'})
print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

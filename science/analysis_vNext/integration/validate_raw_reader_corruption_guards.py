from pathlib import Path
import json,numpy as np
from released_rna import validate_raw_csr_arrays,iter_raw_counts
OUT=Path(__file__).resolve().parent

def run(h5=None, out=None):
 checks={}
 valid=validate_raw_csr_arrays(np.array([100000.,1.]),np.array([0,1]),np.array([0,2]),(1,2))
 assert valid.shape==(1,2)
 trials={
  'large_fractional_count_100000_1_rejected':(np.array([100000.1,1.]),np.array([0,1]),np.array([0,2]),(1,2)),
  'fractional_pointer_dtype_rejected':(np.array([1.]),np.array([0]),np.array([0.,1.]),(1,1)),
  'nonzero_first_pointer_rejected':(np.array([1.]),np.array([0]),np.array([1,1]),(1,1)),
  'wrong_last_pointer_rejected':(np.array([1.]),np.array([0]),np.array([0,0]),(1,1)),
  'nonmonotone_pointer_rejected':(np.array([1.,2.]),np.array([0,1]),np.array([0,2,1]),(2,2)),
  'feature_index_out_of_bounds_rejected':(np.array([1.]),np.array([1]),np.array([0,1]),(1,1)),
  'negative_feature_index_rejected':(np.array([1.]),np.array([-1]),np.array([0,1]),(1,1)),
  'fractional_feature_index_dtype_rejected':(np.array([1.]),np.array([0.]),np.array([0,1]),(1,1)),
  'duplicate_feature_coordinate_rejected':(np.array([1.,2.]),np.array([0,0]),np.array([0,2]),(1,1)),
  'nan_count_rejected':(np.array([np.nan]),np.array([0]),np.array([0,1]),(1,1)),
  'negative_count_rejected':(np.array([-1.]),np.array([0]),np.array([0,1]),(1,1)),
 }
 for name,args in trials.items():
  try:validate_raw_csr_arrays(*args);ok=False
  except ValueError:ok=True
  assert ok,name;checks[name]=ok
 fixture={'status':'skipped: no --h5 supplied; synthetic guards still tested'}
 if h5 is not None:
  ids,genes,x=next(iter_raw_counts(h5))
  fixture={'status':'passed','cells':len(ids),'features':len(genes),'nonzero_counts':int(x.nnz)}
 report={'checks':checks,'optional_real_fixture':fixture}
 if out is not None:
  out.parent.mkdir(parents=True,exist_ok=True)
  out.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))

if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--h5',type=Path);parser.add_argument('--out',type=Path)
 args=parser.parse_args();run(args.h5,args.out)

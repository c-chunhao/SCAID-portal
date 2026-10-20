#!/usr/bin/env python3
"""Build bounded-memory CSC indices of validated RNA H5AD X for lazy gene plots.

Source H5AD, cells, genes, UMAP and RNA values are immutable. Conversion performs
no normalization, feature/cell filtering, cell annotation or batch correction.
Every numeric entry and cell assignment receives full memmap write/readback.
"""
from pathlib import Path
import argparse,collections,datetime,hashlib,json,mmap,os,shutil,tempfile,time
import h5py
import numpy as np

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
TARGETS=('Gastroesophageal_Reflux_Disease/GSE295274/Esophagus','Sarcoidosis/GSE319770/Heart','T1D/GSE316337/PBMC')
FIG_ROOT=FIGURE_ROOT
MANIFEST=GENE_MANIFEST


import sys
sys.path.insert(0,str(BACKEND_ROOT))
from tisch_api.gene_index_builder import build_index,digest,gene_filename

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--object-id',choices=TARGETS,action='append');parser.add_argument('--available-only',action='store_true');parser.add_argument('--max-nnz-chunk',type=int,default=2000000);args=parser.parse_args()
 manifest={'schema_version':1,'figure_root':str(FIG_ROOT),'index_root':str(INDEX_ROOT),'objects':{}}
 if MANIFEST.exists():manifest=json.loads(MANIFEST.read_text())
 for entry in manifest['objects'].values():
  entry.setdefault('filenames',[gene_filename(g) for g in entry['genes']])
  entry.setdefault('filename_encoded_features',[{'gene':g,'filename':gene_filename(g)} for g in entry['genes'] if gene_filename(g)!='Umap_'+g+'.png'])
 inventory=json.loads((WORK/'matrices/download_manifest.json').read_text());source={r['object_id']:r for r in inventory['files'] if r['method']=='RNA' and r.get('validation_status')=='passed'}
 requested=args.object_id or TARGETS;completed=[]
 for obj in requested:
  if obj in manifest['objects']:
   entry=manifest['objects'][obj]
   if digest(entry['index_path'])!=entry['index_sha256']:raise ValueError('Existing gene index checksum mismatch')
   print(json.dumps({'object_id':obj,'status':'already_validated'}),flush=True);continue
  if obj not in source:
   if args.available_only:continue
   raise ValueError('RNA H5AD has not finished validation: '+obj)
  record=source[obj];output=INDEX_ROOT/(obj.replace('/','_')+'.h5');entry=build_index(record['path'],output,obj,record['sha256'],max_nnz_chunk=args.max_nnz_chunk)
  report_path=WORK/'gene_cache/reports'/(obj.replace('/','_')+'.json');report_path.parent.mkdir(parents=True,exist_ok=True);report_path.write_text(json.dumps({k:v for k,v in entry.items() if k not in ('genes','filenames')},ensure_ascii=False,indent=2))
  manifest['objects'][obj]=entry;manifest['updated_at_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();MANIFEST.parent.mkdir(exist_ok=True,parents=True)
  temp=MANIFEST.with_suffix('.json.tmp');temp.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');os.replace(temp,MANIFEST)
  print(json.dumps({k:v for k,v in entry.items() if k not in ['genes','filenames']},ensure_ascii=False),flush=True);completed.append(obj)
 print(json.dumps({'completed_this_run':completed,'registered_gene_objects':list(manifest['objects'])},ensure_ascii=False),flush=True)
if __name__=='__main__':main()

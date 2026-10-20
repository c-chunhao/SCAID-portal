#!/usr/bin/env python3
"""Bounded two-object export queue; no database registration or deployment."""
import concurrent.futures, json, os, subprocess, sys, threading, time
from pathlib import Path

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
PLAN = load_plan()
ENV = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
LOCK = threading.Lock()

def refresh_manifest():
    with LOCK:
        objects=[];files=[]
        for i,obj in enumerate(PLAN['objects']):
            stem=obj['object_id'].replace('/','_')
            p=WORK/'matrices/objects'/stem/'export_report.json'
            row={'object_index':i,'object_id':obj['object_id'],'catalog_cell_count':obj['catalog']['cell_count'],'status':'pending','report_path':str(p)}
            if p.is_file():
                r=json.loads(p.read_text());row['status']=r['status']
                row['actual_cell_count']=r.get('actual_cell_count')
                if r['status']=='validated':
                    for method,f in r['matrices'].items():
                        files.append({'object_id':obj['object_id'],'method':method,
                                      'counts_available':f.get('counts_available',method=='RNA' and 'counts' in f.get('matrix_checks',{})),
                                      'raw_count_status':f.get('raw_count_status','integer_count_layer_preserved_from_processed_source' if method=='RNA' else 'not_applicable_pathway_scores'),
                                      **{k:f[k] for k in ('expression_status','content_description','source_processed_counts_semantics') if k in f},
                                      **{k:f[k] for k in ('path','bytes','sha256','format','shape','validated','validation_status')}})
            objects.append(row)
        manifest={'schema_version':1,'updated':time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                  'object_count':len(objects),'validated_object_count':sum(x['status']=='validated' for x in objects),
                  'file_count':len(files),'total_validated_bytes':sum(x['bytes'] for x in files),
                  'objects':objects,'files':files,
                  'scope':'Existing processed Seurat data converted without filtering or relabeling; sample labels are not authenticated donor identifiers.'}
        dst=WORK/'matrices/download_manifest.json';dst.parent.mkdir(parents=True,exist_ok=True);tmp=dst.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(manifest,ensure_ascii=False,indent=2));tmp.replace(dst)
        return manifest

def run(i):
    obj=PLAN['objects'][i];stem=obj['object_id'].replace('/','_')
    p=WORK/'matrices/objects'/stem/'export_report.json';p.parent.mkdir(parents=True,exist_ok=True)
    if p.is_file() and json.loads(p.read_text())['status']=='validated':
        return {'index':i,'status':'already_validated'}
    log=p.parent/'export.log'
    with log.open('a') as stream:
        current=json.loads(p.read_text()) if p.is_file() else {}
        if current.get('status')!='exported_awaiting_python_validation':
            c=subprocess.run([os.environ.get('SCAID_RSCRIPT','Rscript'),str(CODE_ROOT/'export_new_aid_h5ad.R'),str(i),str(WORK)],env=ENV,stdout=stream,stderr=subprocess.STDOUT)
            if c.returncode:
                refresh_manifest();return {'index':i,'status':'R_export_failed','log':str(log)}
        c=subprocess.run([sys.executable,'-u',str(CODE_ROOT/'validate_new_aid_h5ad.py'),str(p)],env=ENV,stdout=stream,stderr=subprocess.STDOUT)
    manifest=refresh_manifest()
    print('DONE',obj['object_id'],'return',c.returncode,'validated_objects',manifest['validated_object_count'],'files',manifest['file_count'],flush=True)
    return {'index':i,'status':'validated' if c.returncode==0 else 'validation_failed','log':str(log)}

if __name__=='__main__':
    refresh_manifest()
    queue=sorted(range(len(PLAN['objects'])),key=lambda i:PLAN['objects'][i]['catalog']['cell_count'])
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(run,queue))
    manifest=refresh_manifest();(WORK/'matrices/batch_result.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
    print('FINAL',manifest['validated_object_count'],manifest['file_count'],manifest['total_validated_bytes'],flush=True)
    if any(r['status'] not in ('validated','already_validated') for r in results):sys.exit(1)

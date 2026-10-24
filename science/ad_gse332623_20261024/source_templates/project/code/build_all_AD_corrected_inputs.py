#!/usr/bin/env python3
"""Two-process, resumable, independent source-library ingestion; no analysis runner."""
import os, json, subprocess, time, concurrent.futures, hashlib
from pathlib import Path
SOURCE=Path('@RAW_INPUT@')
OUTPUT=Path('@PROJECT@/00.corrected_data/GSE332623')
SCRIPT=Path('@PROJECT@/code/GSE332623_build_corrected_input.R')
LOGS=OUTPUT/'logs'

def run(library):
    LOGS.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',RCPP_PARALLEL_NUM_THREADS='1')
    start=time.monotonic()
    with (LOGS/(library+'.log')).open('a') as log:
        result=subprocess.run(['/usr/bin/Rscript',str(SCRIPT),library,str(OUTPUT)],env=env,stdout=log,stderr=subprocess.STDOUT)
    receipt=OUTPUT/'receipts'/(library+'.json')
    return {'library':library,'exit_code':result.returncode,'elapsed_seconds':round(time.monotonic()-start,2),'receipt_exists':receipt.exists()}

if __name__=='__main__':
    libraries=sorted(x.name.removesuffix('_matrix.mtx.gz') for x in SOURCE.glob('GSE332623_TotalseqB_NEA[1-7]_ADT[1-4]_matrix.mtx.gz'))
    if len(libraries)!=27:raise SystemExit('Expected exact 27 existing libraries')
    OUTPUT.mkdir(parents=True,exist_ok=True)
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(run,u)for u in libraries]
        for future in concurrent.futures.as_completed(futures):
            item=future.result();results.append(item);print(json.dumps(item),flush=True)
            (OUTPUT/'run_progress.json').write_text(json.dumps({'workers':2,'results':results},indent=2))
    failed=[x for x in results if x['exit_code']!=0 or not x['receipt_exists']]
    receipts=[json.loads((OUTPUT/'receipts'/(u+'.json')).read_text())for u in libraries if (OUTPUT/'receipts'/(u+'.json')).exists()]
    report={'stage':'Counts-only source input; no QC, disease assignment, normalization, embeddings or annotation',
            'status':'passed'if not failed and len(receipts)==27 and all(x['retained_all_source_barcodes'] and x['Group']=='Unresolved' for x in receipts) and sum(x['n_cells']for x in receipts)==754086 else'failed','workers':2,'n_libraries':len(receipts),
            'total_source_barcode_records':sum(x['n_cells']for x in receipts),'RNA_features_by_library':sorted(set(x['rna_features']for x in receipts)),
            'ADT_features_by_library':sorted(set(x['adt_features']for x in receipts)),'all_clinical_groups_unresolved':all(x['Group']=='Unresolved'for x in receipts),
            'all_source_axes_retained':all(x['retained_all_source_barcodes']for x in receipts),'output_bytes':sum(x['output_size_bytes']for x in receipts),
            'libraries':receipts,'process_results':results,'failures':failed}
    (OUTPUT/'manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
    if report['status']!='passed':raise SystemExit(1)
    print(json.dumps({k:v for k,v in report.items()if k not in ['libraries','process_results']}),flush=True)

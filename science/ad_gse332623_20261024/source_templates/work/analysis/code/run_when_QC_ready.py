#!/usr/bin/env python3
"""Pilot on the first validated QC bridge; then analyze all27completed bridges."""
from pathlib import Path
import os,sys,time,json,csv,subprocess,hashlib
ROOT=Path('@PROJECT@/00.corrected_analysis/GSE332623')
WORK=Path('@WORK@/analysis')
SCRIPT=Path(__file__).with_name('analysis_pipeline.py');OUT=ROOT/'integration'
ENV=os.environ.copy();ENV.update(OPENBLAS_NUM_THREADS='8',OMP_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8',VECLIB_MAXIMUM_THREADS='8')
def record(**info):
 info['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());p=WORK/'progress.json';tmp=p.with_suffix('.json.tmp');tmp.write_text(json.dumps(info,ensure_ascii=False,indent=2));os.replace(tmp,p)
def sources():
 with(ROOT/'capture_manifest.tsv').open()as f:rows=list(csv.DictReader(f,delimiter='\t'))
 ready=[]
 for r in rows:
  d=ROOT/'qc'/r['library_id'];receipt=d/'receipt.json';bridge=d/'counts_bridge.h5'
  if receipt.exists()and bridge.exists():
   try:ok=json.loads(receipt.read_text()).get('status')=='validated'
   except Exception:ok=False
   if ok:ready.append(bridge)
 return rows,ready

def run_analysis(out,bridges,pilot):
 stages=[('assemble','checkpoints/assembly.json'),('classify','checkpoints/reference_classification.json'),('embed','checkpoints/embedding.json'),('finalize','pilot_ready.json'if pilot else'analysis_candidate.json')]
 for stage,checkpoint in stages:
  p=out/checkpoint
  if p.exists():
   try:
    d=json.loads(p.read_text())
    if d.get('status')=='validated':continue
   except Exception:pass
  record(status='running',scope='pilot'if pilot else'full_cohort',stage=stage,n_QC_libraries=len(bridges),out=str(out))
  cmd=[sys.executable,str(SCRIPT),stage,'--out',str(out)]
  if stage=='assemble':cmd+=['--bridges']+list(map(str,bridges))
  if pilot:cmd+=['--pilot']
  log=WORK/'logs'/f'{"pilot"if pilot else"full"}_{stage}.log'
  resources=WORK/'logs'/f'{"pilot"if pilot else"full"}_{stage}.resources.txt'
  timed=['/usr/bin/time','-v','-o',str(resources)]+cmd
  with log.open('a')as f:res=subprocess.run(timed,env=ENV,stdout=f,stderr=subprocess.STDOUT)
  if res.returncode:
   record(status='failed',scope='pilot'if pilot else'full_cohort',stage=stage,exit_code=res.returncode,log=str(log));raise SystemExit(res.returncode)

if __name__=='__main__':
 WORK.mkdir(parents=True,exist_ok=True);(WORK/'logs').mkdir(exist_ok=True);pilot_done=(OUT/'pilot/pilot_ready.json').exists()
 while True:
  rows,ready=sources()
  if ready and not pilot_done:
   run_analysis(OUT/'pilot',[ready[0]],True);pilot_done=True
  if len(ready)==len(rows)==27:
   run_analysis(OUT,ready,False)
   record(status='awaiting_analytical_review',scope='full_cohort',stage='candidate_complete',n_QC_libraries=27,analysis_candidate=str(OUT/'analysis_candidate.json'))
   while not(OUT/'analysis_ready.json').exists():time.sleep(20)
   release=json.loads((OUT/'analysis_ready.json').read_text());assert release['status']=='validated'and release['n_libraries']==27
   record(status='validated',scope='full_cohort',stage='complete',n_QC_libraries=27,analysis_ready=str(OUT/'analysis_ready.json'));break
  record(status='waiting_for_QC',n_QC_libraries=len(ready),expected=27,pilot_completed=pilot_done);time.sleep(20)

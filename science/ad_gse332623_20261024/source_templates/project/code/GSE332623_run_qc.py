#!/usr/bin/env python3
"""Resume 27 independent RNA QC channels with four bounded systemd scopes."""
import csv
from datetime import datetime, timezone, timedelta
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

PROJECT = Path('@PROJECT@')
WORK = Path('@WORK@')
BASE = PROJECT / '00.corrected_analysis/GSE332623'
CODE = Path(__file__).resolve().parent
MAX_JOBS = 4
MEMORY_GIB_PER_JOB = 80
RESERVE_GIB = 192


def now():
    return datetime.now(timezone(timedelta(hours=8))).isoformat()


def file_hash(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(8 * 1024 * 1024): digest.update(block)
    return digest.hexdigest()


def available_gib():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'): return int(line.split()[1]) / (1024 * 1024)
    raise RuntimeError('Cannot inspect available memory')


def atomic_json(path, value):
    temp = path.with_name(path.name + '.tmp-' + str(os.getpid()))
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temp.replace(path)


def main():
    import fcntl
    BASE.joinpath('qc').mkdir(exist_ok=True)
    lock = (BASE / 'qc/queue.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    rows = list(csv.DictReader((BASE / 'capture_manifest.tsv').open(), delimiter='\t'))
    assert len(rows) == 27 and len({row['library_id'] for row in rows}) == 27
    code_hashes = {name: file_hash(CODE / name) for name in ('GSE332623_rna_qc.R', 'GSE332623_finalize_qc.R')}
    stamp = datetime.now().strftime('%Y%m%dT%H%M%S')
    active, pending, results = {}, [], []
    for row in rows:
        path = BASE / 'qc' / row['library_id'] / 'receipt.json'
        if path.exists():
            data = json.loads(path.read_text())
            if data.get('status') != 'validated' or data['source_sha256'] != row['input_sha256']:
                raise ValueError(f'Unverified existing QC receipt: {path}')
            if file_hash(Path(data['bridge_path'])) != data['bridge_sha256']:
                raise ValueError(f'QC bridge changed: {path}')
            results.append({'library_id': row['library_id'], 'status': 'validated_existing', 'n_retained': data['n_retained']})
        else:
            pending.append(row)
    state_path = BASE / 'qc/run_progress.json'
    while pending or active:
        for library, item in list(active.items()):
            process = item['process']
            code = process.poll()
            if code is not None:
                item['stream'].close()
                path = BASE / 'qc' / library / 'receipt.json'
                success = code == 0 and path.exists() and json.loads(path.read_text()).get('status') == 'validated'
                results.append({'library_id': library, 'status': 'validated' if success else 'failed',
                                'exit_code': code, 'log': item['log'], 'unit': item['unit'],
                                'started': item['started'], 'finished': now(), 'maximum_observed_memory_bytes': item['peak']})
                del active[library]
            else:
                current = subprocess.run(['systemctl','--user','show', item['unit'] + '.scope',
                                          '-p','MemoryPeak','--value'], capture_output=True, text=True)
                value = current.stdout.strip()
                if value.isdigit(): item['peak'] = max(item['peak'], int(value))
        while pending and len(active) < MAX_JOBS and available_gib() >= RESERVE_GIB:
            row = pending.pop(0); library = row['library_id']
            unit = 'scaid-ad-qc-' + library.replace('GSE332623_TotalseqB_', '').lower().replace('_','-') + '-' + stamp.lower()
            log = WORK / 'logs' / (library + '_' + stamp + '.log')
            stream = log.open('a')
            command = ['systemd-run','--user','--scope','--quiet','--unit='+unit,
                       '-p',f'MemoryMax={MEMORY_GIB_PER_JOB}G','-p','MemorySwapMax=1G','-p','CPUQuota=400%',
                       'env','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=4','MKL_NUM_THREADS=1','R_MAX_VSIZE=72Gb',
                       'bash','-c','Rscript "$1" "$2" && Rscript "$3" "$2"','scaid-qc',
                       str(CODE / 'GSE332623_rna_qc.R'),library,str(CODE / 'GSE332623_finalize_qc.R')]
            process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            active[library] = {'process':process,'stream':stream,'log':str(log),'unit':unit,
                               'started':now(),'peak':0}
        value = {'status':'running','updated':now(),'n_total':len(rows),'n_finished':len(results),
                 'n_validated':sum(x['status'].startswith('validated') for x in results),
                 'n_failed':sum(x['status']=='failed' for x in results),'n_pending':len(pending),
                 'active':{key:{k:v for k,v in item.items() if k not in ('process','stream')}
                           for key,item in active.items()},'results':results,
                 'resource_policy':{'max_concurrent':MAX_JOBS,'memory_max_GiB_each':MEMORY_GIB_PER_JOB,
                                    'aggregate_QC_hard_budget_GiB':MAX_JOBS*MEMORY_GIB_PER_JOB,
                                    'CPU_quota_cores_each':4,'min_available_GiB_before_dispatch':RESERVE_GIB},
                 'code_sha256':code_hashes,'available_memory_GiB':available_gib()}
        atomic_json(state_path, value)
        if pending or active: time.sleep(10)
    value['status'] = 'complete' if value['n_failed'] == 0 else 'failed'
    atomic_json(state_path, value)
    atomic_json(WORK / 'validation/QC_queue_completion.json', value)
    return 0 if value['status']=='complete' else 1


if __name__ == '__main__':
    raise SystemExit(main())

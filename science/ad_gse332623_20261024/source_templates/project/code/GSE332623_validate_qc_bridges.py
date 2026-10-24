#!/usr/bin/env python3
"""Independent chunked count/axis/flag verification of completed QC bridge files."""
import csv
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import h5py
import numpy as np
import pandas as pd

BASE=Path('@PROJECT@/00.corrected_analysis/GSE332623')
WORK=Path('@WORK@')


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        while b:=stream.read(8*1024*1024):h.update(b)
    return h.hexdigest()


def validate(entry):
    directory=BASE/'qc'/entry['library_id'];receipt_path=directory/'receipt.json'
    receipt=json.loads(receipt_path.read_text());assert receipt['status']=='validated'
    target=Path(receipt['bridge_path']);receipt_hash=digest(receipt_path)
    result_path=directory/'independent_bridge_QA.json'
    if result_path.exists():
        old=json.loads(result_path.read_text())
        if old.get('passed') and old['receipt_sha256']==receipt_hash:
            return old
    frame=pd.read_csv(directory/'per_cell_QC.tsv.gz',sep='\t')
    keep=frame.passes_QC.astype(str).str.lower().eq('true')
    basic=frame.passes_basic_qc.astype(str).str.lower().eq('true')
    assert (keep==(basic & frame.doublet_class.eq('singlet'))).all()
    assert len(frame)==int(entry['source_barcodes'])==receipt['n_input']
    assert int(keep.sum())==receipt['n_retained']
    filtered=frame.loc[keep].reset_index(drop=True)
    assert frame.cell_id.is_unique
    assert frame.cell_id.tolist()==[entry['library_id']+'_'+barcode for barcode in frame.source_barcode]
    assert frame.Group.eq('Unresolved').all() and frame.source_condition.eq('Unresolved').all()
    all_records={}
    with h5py.File(target,'r') as file:
        assert file.attrs['schema']=='scaid-rna-adt-counts-v1'
        assert file['obs/cell_id'].asstr()[:].tolist()==filtered.cell_id.tolist()
        assert file['obs/source_barcode'].asstr()[:].tolist()==filtered.source_barcode.tolist()
        for field in ('library_id','capture_id'):
            assert set(file['obs/'+field].asstr()[:])=={entry[field]}
        assert set(file['obs/Group'].asstr()[:])=={'Unresolved'}
        assert set(file['obs/source_condition'].asstr()[:])=={'Unresolved'}
        for method,nfeatures in [('RNA',33694),('ADT',140)]:
            group=file[method];shape=group['shape'][:].astype(np.int64)
            assert shape.tolist()==[len(filtered),nfeatures]
            p=group['indptr'][:].astype(np.int64)
            assert len(p)==len(filtered)+1 and p[0]==0 and (np.diff(p)>=0).all()
            assert p[-1]==len(group['data'])==len(group['indices'])
            assert len(group['feature_ids'])==len(group['source_symbols'])==nfeatures
            names=group['feature_names'].asstr()[:];assert len(set(names))==nfeatures
            row_sums=np.zeros(len(filtered),dtype=np.float64);detected=np.zeros(len(filtered),dtype=np.int64)
            for lo in range(0,len(filtered),512):
                hi=min(len(filtered),lo+512);start,end=int(p[lo]),int(p[hi])
                values=group['data'][start:end];indices=group['indices'][start:end]
                assert np.isfinite(values).all() and (values>=0).all() and np.array_equal(values,np.floor(values))
                assert (indices>=0).all() and (indices<nfeatures).all()
                offsets=p[lo:hi+1]-start
                sums=np.concatenate(([0.],np.cumsum(values,dtype=np.float64)))
                nnz=np.concatenate(([0],np.cumsum(values>0,dtype=np.int64)))
                row_sums[lo:hi]=sums[offsets[1:]]-sums[offsets[:-1]]
                detected[lo:hi]=nnz[offsets[1:]]-nnz[offsets[:-1]]
            assert np.array_equal(row_sums,filtered['nCount_'+method].to_numpy())
            assert np.array_equal(detected,filtered['nFeature_'+method].to_numpy())
            all_records[method]={'shape':shape.tolist(),'nnz':int(p[-1]),'all_values_finite_nonnegative_integer':True,
                                 'all_row_counts_and_detected_features_match_QC':True,'axis_complete':True}
        assert file['RNA_decontX/shape'][:].tolist()==file['RNA/shape'][:].tolist()
        assert file['RNA_decontX/feature_names'].asstr()[:].tolist()==file['RNA/feature_names'].asstr()[:].tolist()
    assert digest(target)==receipt['bridge_sha256']
    result={'passed':True,'library_id':entry['library_id'],'n_input':len(frame),'n_retained':len(filtered),
            'receipt_sha256':receipt_hash,'bridge_sha256':receipt['bridge_sha256'],'matrices':all_records,
            'all_source_flags_accounted':True,'source_cell_ids_prefix_exact':True,'clinical_assignment':'Unresolved',
            'scope':'Independent HDF/Python all-value raw-count check and all-cell metadata check; source MEX/RDS validation retained in original import receipts.'}
    temp=result_path.with_suffix('.tmp');temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(result_path)
    return result


def main():
    entries=list(csv.DictReader((BASE/'capture_manifest.tsv').open(),delimiter='\t'))
    results=[]
    for entry in entries:
        if (BASE/'qc'/entry['library_id']/'receipt.json').exists():
            results.append(validate(entry));print(entry['library_id'],'independent QA passed',flush=True)
    value={'passed_for_completed':all(x['passed'] for x in results),'n_validated':len(results),'n_expected':27,
           'complete':len(results)==27,'n_input_verified':sum(x['n_input'] for x in results),
           'n_retained_verified':sum(x['n_retained'] for x in results),'results':results}
    (WORK/'validation/QC_independent_bridge_QA.json').write_text(json.dumps(value,indent=2)+'\n')
    return 0 if value['passed_for_completed'] else 1


if __name__=='__main__':raise SystemExit(main())

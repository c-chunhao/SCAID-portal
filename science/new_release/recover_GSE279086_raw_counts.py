#!/usr/bin/env python3
"""Recover integer counts by exact GSM/barcode/feature mapping; retain processed values."""
import json, os, re, time
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import gzip, hashlib
import anndata as ad
from anndata.io import read_elem, sparse_dataset
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
OUT=WORK/'matrices/source_recovery_GSE279086'
REPORT=WORK/'matrices/objects/T1D_GSE279086_Kidney/export_report.json'
SOURCE=RAW_ROOT

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        while b:=f.read(8*1024*1024):h.update(b)
    return h.hexdigest()

def unique_r(names):
    """R make.unique's reserved-original-name rule, with the default '.' separator."""
    reserved=set(names);used=set();result=[];next_suffix={}
    for name in names:
        if name not in used:
            value=name
        else:
            i=next_suffix.get(name,1)
            while f'{name}.{i}' in reserved or f'{name}.{i}' in used:i+=1
            value=f'{name}.{i}';next_suffix[name]=i+1
        used.add(value);result.append(value)
    return result

def one(directory,pattern):
    paths=list(directory.glob(pattern))
    if len(paths)!=1:raise ValueError(f'Expected one {pattern} in {directory}, found {len(paths)}')
    return paths[0]

def lines(path):
    with gzip.open(path,'rt') as f:return [s.rstrip('\r\n') for s in f]

def matrix(path):
    with gzip.open(path,'rt') as f:
        header=f.readline().strip()
    with threadpool_limits(limits=2),gzip.open(path,'rb') as f:m=mmread(f).tocsr()
    return m,header

def mapping(axis,requested,what):
    if len(set(axis))!=len(axis):raise ValueError('Non-unique '+what)
    lookup={s:i for i,s in enumerate(axis)}
    missing=[s for s in requested if s not in lookup]
    if missing:raise ValueError(f'{what} mapping missing {len(missing)} values: {missing[:12]}')
    ix=np.asarray([lookup[s] for s in requested],dtype=np.int64)
    if len(set(ix))!=len(ix):raise ValueError('Mapping is not one-to-one: '+what)
    return ix

def put_array(group,name,values):
    if name in group:raise ValueError('Refusing to replace HDF5 item '+name)
    d=group.create_dataset(name,data=values,compression='gzip',compression_opts=1,
                           chunks=True)
    d.attrs['encoding-type']='array';d.attrs['encoding-version']='0.2.0'
    return d

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    audit_path=OUT/'axis_audit.json'
    if audit_path.exists() and not json.loads(audit_path.read_text()).get('all_26_samples_unique_raw_axis_mapping',False):
        raise ValueError('Full-source raw-axis mapping is incomplete. Use the declared processed-expression-only main file and the independently verified 24-GSM raw subset; do not recover a complete count layer.')
    result=json.loads(REPORT.read_text())
    if result['status']!='exported_awaiting_raw_count_recovery':
        raise ValueError('Expected an exported processed matrix awaiting recovery')
    path=Path(result['matrices']['RNA']['temporary_path'])
    with h5py.File(path,'r') as h:
        obs=read_elem(h['obs']);var=read_elem(h['var'])
        genes=var.index.tolist()
        if 'source_processed_counts' not in h['layers'] or 'counts' in h['layers']:
            raise ValueError('Processed/raw matrix layers do not match recovery policy')
    counts_by_sample=[];positions=[];audits=[];files=[];raw_gene_ids=None
    for gsm in sorted(obs['GSM'].astype(str).unique()):
        print('RECOVER',gsm,flush=True)
        directory=SOURCE/gsm
        paths={k:one(directory,pattern) for k,pattern in {
            'raw_matrix':'*_matrix.mtx.gz','raw_barcodes':'*_barcodes.tsv.gz','raw_features':'*_features.tsv.gz',
            'processed_matrix':'*_matrix_processed.mtx.gz','processed_barcodes':'*_barcodes_processed.tsv.gz',
            'processed_features':'*_features_processed.tsv.gz'}.items()}
        raw_features=[s.split('\t') for s in lines(paths['raw_features'])]
        if any(len(s)<2 for s in raw_features):raise ValueError('Raw features lack identifiers/symbols')
        raw_names=unique_r([s[1] for s in raw_features])
        # Seurat feature-name sanitation is only accepted when no such transformation is needed.
        if any('_' in g or '|' in g for g in genes):raise ValueError('Target names need a feature sanitation policy')
        gi=mapping(raw_names,genes,'raw gene symbols made unique using the source reader rule')
        gene_ids=[raw_features[k][0] for k in gi]
        if raw_gene_ids is None:raw_gene_ids=gene_ids
        elif raw_gene_ids!=gene_ids:raise ValueError('Gene identifier mapping differs between GSMs')
        global_positions=np.flatnonzero(obs['GSM'].astype(str).eq(gsm).to_numpy())
        cell_ids=obs.index[global_positions].tolist()
        if not all(s.startswith(gsm+'_') for s in cell_ids):raise ValueError('Cell GSM prefix differs from metadata')
        processed_cells=[s[len(gsm)+1:] for s in cell_ids]
        raw_cells=[s.rsplit('_',1)[-1] for s in cell_ids]
        if not all(re.fullmatch(r'[ACGTN]+-[0-9]+',s) for s in raw_cells):raise ValueError('Unexpected raw barcode format')
        rb=lines(paths['raw_barcodes']);pb=lines(paths['processed_barcodes'])
        ci=mapping(rb,raw_cells,'raw barcode')
        pci=mapping(pb,processed_cells,'processed barcode')
        raw,header=matrix(paths['raw_matrix'])
        if 'coordinate integer general' not in header:raise ValueError('Raw MatrixMarket is not declared integer')
        if raw.shape!=(len(raw_features),len(rb)):raise ValueError('Raw axes differ from header dimensions')
        if not np.isfinite(raw.data).all() or (raw.data<0).any() or not np.equal(raw.data,np.rint(raw.data)).all():
            raise ValueError('Raw source has invalid counts')
        block=raw[gi,:][:,ci].T.tocsr();block.sort_indices()
        if block.data.max(initial=0)>np.iinfo(np.int32).max:raise ValueError('Counts exceed int32 storage')
        block=block.astype(np.int32)
        del raw
        # Check the processed counterpart against every retained cell/feature of the RDS layer.
        pf=lines(paths['processed_features'])
        if len(set(pf))!=len(pf):raise ValueError('Processed feature names are non-unique')
        proc,pheader=matrix(paths['processed_matrix'])
        if proc.shape!=(len(pf),len(pb)):raise ValueError('Processed axes differ from source dimensions')
        plook={g:k for k,g in enumerate(pf)}
        present=np.asarray([k for k,g in enumerate(genes) if g in plook],dtype=np.int64)
        pi=np.asarray([plook[genes[k]] for k in present],dtype=np.int64)
        small=proc[pi,:][:,pci].T.tocoo()
        processed=sparse.csr_matrix((small.data,(small.row,present[small.col])),shape=(len(ci),len(genes)))
        processed.sort_indices();del proc,small
        with h5py.File(path,'r') as h:source=sparse_dataset(h['layers/source_processed_counts'])[global_positions,:].tocsr()
        source.sort_indices()
        if not np.array_equal(source.indptr,processed.indptr) or not np.array_equal(source.indices,processed.indices):
            raise ValueError('Processed source support differs from retained RDS layer')
        if not np.allclose(source.data,processed.data,rtol=1e-13,atol=1e-14):
            raise ValueError('Processed source values differ from retained RDS layer')
        maximum=float(np.max(np.abs(source.data-processed.data),initial=0))
        audit={'GSM':gsm,'retained_cells':len(ci),'retained_genes':len(gi),'raw_matrix_shape':list((len(raw_features),len(rb))),
               'raw_header':header,'processed_header':pheader,'raw_feature_mapping_unique':True,'raw_barcode_mapping_unique':True,
               'processed_feature_mapping_exact':True,'processed_barcode_mapping_exact':True,
               'processed_full_sparse_support_equal':True,'processed_maximum_absolute_value_difference':maximum,
               'recovered_count_nnz':block.nnz,'source_files':{}}
        for kind,f in paths.items():
            row={'GSM':gsm,'kind':kind,'path':str(f),'bytes':f.stat().st_size,'sha256':sha(f)}
            audit['source_files'][kind]=row;files.append(row)
        np.savez_compressed(OUT/f'{gsm}_exact_axis_mapping.npz',
                            retained_global_cell_positions=global_positions,raw_cell_positions=ci,raw_gene_positions=gi,
                            processed_cell_positions=pci,processed_present_gene_positions=present,processed_gene_positions=pi)
        audits.append(audit);counts_by_sample.append(block);positions.extend(global_positions.tolist())
    merged=sparse.vstack(counts_by_sample,format='csr')
    perm=np.argsort(np.asarray(positions));counts=merged[perm,:].tocsr();del merged,counts_by_sample
    if counts.shape!=(len(obs),len(genes)) or len(set(positions))!=len(obs):raise ValueError('Recovered global matrix axes incomplete')
    counts.sort_indices();counts.eliminate_zeros()
    totals=np.asarray(counts.sum(1)).ravel().astype(np.float64)
    detected=np.diff(counts.indptr).astype(np.float64)
    mt=np.asarray([str(g).upper().startswith('MT-') for g in genes])
    mito=np.divide(np.asarray(counts[:,mt].sum(1)).ravel()*100.,totals,out=np.zeros(len(obs)),where=totals>0)
    with h5py.File(path,'r+') as h:
        group=h['layers'].create_group('counts');group.attrs['encoding-type']='csr_matrix';group.attrs['encoding-version']='0.1.0';group.attrs['shape']=counts.shape
        for name,values in [('data',counts.data),('indices',counts.indices),('indptr',counts.indptr)]:
            put_array(group,name,values)
            for start in range(0,len(values),2_000_000):
                if not np.array_equal(group[name][start:start+2_000_000],values[start:start+2_000_000]):raise ValueError('Raw count full-component readback differs')
        obs_group=h['obs'];columns=list(obs_group.attrs['column-order'])
        columns=[s.decode() if isinstance(s,bytes) else str(s) for s in columns]
        for name,values in [('nCount_source_raw',totals),('nFeature_source_raw',detected),('percent_mt_source_raw',mito)]:
            put_array(obs_group,name,values);columns.append(name)
        del obs_group.attrs['column-order'];obs_group.attrs.create('column-order',np.asarray(columns,dtype=h5py.string_dtype('utf-8')))
        v=h['var'];vc=list(v.attrs['column-order']);vc=[s.decode() if isinstance(s,bytes) else str(s) for s in vc]
        for name,values in [('gene_id_source_raw',raw_gene_ids)]:
            d=v.create_dataset(name,data=np.asarray(values,dtype=h5py.string_dtype('utf-8')))
            d.attrs['encoding-type']='string-array';d.attrs['encoding-version']='0.2.0';vc.append(name)
        del v.attrs['column-order'];v.attrs.create('column-order',np.asarray(vc,dtype=h5py.string_dtype('utf-8')))
        description='Integer raw matrices were matched to the retained cells and features by GSM, exact 10x barcode and source-reader make.unique gene symbols. Original fractional processed counts are preserved as source_processed_counts. X and original QC fields remain derived from the processed source; *_source_raw QC columns describe the recovered raw-count layer.'
        d=h['uns/scaid'].create_dataset('raw_count_recovery',data=np.asarray([description],dtype=h5py.string_dtype('utf-8')))
        d.attrs['encoding-type']='string-array';d.attrs['encoding-version']='0.2.0'
    rr=np.unique(np.round(np.linspace(0,counts.shape[0]-1,19)).astype(int));cc=np.unique(np.round(np.linspace(0,counts.shape[1]-1,17)).astype(int))
    pairs=[(int(r),int(c)) for c in cc for r in rr]
    probes={'cell_index_zero':[r for r,c in pairs],'feature_index_zero':[c for r,c in pairs],
            'values':[int(counts[r,c]) for r,c in pairs]}
    result['matrices']['RNA']['matrix_checks']['counts']={'source_class':'scipy.csr_matrix recovered from declared integer source MatrixMarket',
        'shape':list(counts.shape),'nnz':counts.nnz,'minimum':0,'maximum':int(counts.data.max(initial=0)),
        'finite':True,'full_component_readback':True,'nonnegative_integer_values':True,'source_value_probes':probes}
    result['matrices']['RNA']['additional_derived_metadata']=['nCount_source_raw','nFeature_source_raw','percent_mt_source_raw']
    recovery={'status':'passed','source':'GSE279086 integer raw matrices stored alongside processed matrices',
              'source_GSM_count':len(audits),'retained_cell_count':len(obs),'retained_gene_count':len(genes),
              'source_processed_counts_preserved':True,'normalized_X_unchanged':True,'original_QC_metadata_unchanged':True,
              'raw_counts_nonnegative_integer':True,'full_sparse_component_readback':True,
              'recovery_interpretation':description,'sample_audits':audits,'source_files':files,
              'completed':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (OUT/'raw_count_recovery.json').write_text(json.dumps(recovery,ensure_ascii=False,indent=2))
    pd.DataFrame(files).to_csv(OUT/'source_file_hashes.tsv',sep='\t',index=False)
    result['raw_count_recovery_report']=str(OUT/'raw_count_recovery.json')
    result['status']='exported_awaiting_python_validation'
    temp=REPORT.with_suffix('.json.tmp');temp.write_text(json.dumps(result,ensure_ascii=False,indent=2));temp.replace(REPORT)
    print('RECOVERY_PASSED',len(audits),counts.shape,counts.nnz,flush=True)

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Run checkpointed KEGG scoring; publish three verified standard H5AD files."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import h5py
import numpy as np
import pandas as pd
import anndata as ad
from contextlib import ExitStack

CODE=Path(__file__).resolve().parent
WORK=CODE.parent
GENESETS=Path('@PROJECT@/20240526_kegg_hsa_list.rds')
CORRECTED=Path('@PROJECT@/00.corrected_analysis/GSE332623')
METHODS=('AUCell','UCell','singscore')
PRODUCTION_PYTHON=Path('@WEB_PYTHON@')
SYMBOL_ALIASES=WORK.parent/'analysis/validation/hgnc_alias_review/cellchat_symbol_aliases_verified.tsv'
SYMBOL_ALIAS_SHA='51873d0026a7cb6657298baf266f7de7c4f9b1b7ab132e4b2193c932703b18cd'

def verify_symbol_aliases(input_path,path=SYMBOL_ALIASES,sha=SYMBOL_ALIAS_SHA,count=9):
    assert digest(path)==sha
    aliases=pd.read_csv(path,sep='\t')
    assert len(aliases)==count and aliases.source_symbol.is_unique and aliases.canonical_symbol.is_unique
    confirmed=aliases.alias_confirmed if 'alias_confirmed' in aliases else aliases.source_alias_confirmed
    assert confirmed.all() and not aliases.canonical_present_in_RNA.any()
    if 'ensembl_exact_match' in aliases:assert aliases.ensembl_exact_match.all()
    else:
        assert aliases.accepted.all() and not aliases.identity_mapping.any()
        for row in aliases.itertuples():assert row.source_ensembl_id in row.canonical_reference_Ensembl_ids.split(';')
    with h5py.File(input_path,'r') as h:
        genes=h['var'][h['var'].attrs['_index']].asstr()[...].tolist()
        feature_ids=h['var/source_feature_ids'].asstr()[...].tolist()
        for row in aliases.itertuples():
            assert row.source_symbol in genes and row.canonical_symbol not in genes
            assert row.source_ensembl_id in feature_ids[genes.index(row.source_symbol)].split(';')
    return aliases

def production_reader_check(path,shape):
    program="import anndata as ad,json,sys; a=ad.read_h5ad(sys.argv[1],backed='r'); assert a.shape==tuple(json.loads(sys.argv[2])); assert a.obs_names.is_unique and a.var_names.is_unique; assert 'X_umap' in a.obsm; print(json.dumps({'anndata_version':ad.__version__,'backed_read':True,'shape':list(a.shape)})); a.file.close()"
    env=os.environ.copy();env.update({x:'1' for x in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    result=subprocess.run([str(PRODUCTION_PYTHON),'-c',program,str(path),json.dumps(list(shape))],env=env,check=True,capture_output=True,text=True)
    return json.loads(result.stdout)

def digest(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def axis_digest(axis):return hashlib.sha256('\n'.join(map(str,axis)).encode()).hexdigest()

def save(path,data):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n');os.replace(tmp,path)

def declared_hashes(value):
    if isinstance(value,dict):
        for key,item in value.items():
            if key.lower().endswith('sha256') and isinstance(item,str):yield item
            yield from declared_hashes(item)
    elif isinstance(value,list):
        for item in value:yield from declared_hashes(item)

def boundary():
    os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:12])
    resource.setrlimit(resource.RLIMIT_AS,(92*1024**3,92*1024**3))

def prepare(input_path,output,mode,ready,chunk_size,object_id,workers):
    input_path=Path(input_path).resolve();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    assert 1<=chunk_size<=5000
    assert 1<=workers<=8 and (mode=='formal' or workers<=4)
    if mode=='formal':
        assert output.is_relative_to(CORRECTED) and object_id=='Atopic_Dermatitis/GSE332623/PBMC'
        assert ready and Path(ready).is_file(),'Formal run requires validated integration readiness receipt'
        receipt=json.loads(Path(ready).read_text());assert receipt.get('status')=='validated'
        assert receipt.get('n_libraries')==27 and Path(receipt['RNA_path']).resolve()==input_path
    else:receipt=None
    source_hash=digest(input_path)
    if receipt:assert source_hash in set(declared_hashes(receipt)),'RNA file SHA not present in integration receipt'
    source=ad.read_h5ad(input_path,backed='r')
    assert source.n_obs>0 and source.n_vars>500
    assert source.obs_names.is_unique and source.var_names.is_unique
    assert 'MajorCellType' in source.obs and 'X_umap' in source.obsm
    if 'feature_type' in source.var:
        assert set(source.var.feature_type.astype(str))<= {'Gene Expression','RNA'},'ADT or mixed feature types in RNA input'
    genes=source.var_names.astype(str).tolist();cells=source.obs_names.astype(str).tolist()
    assert not any('\n' in x or '\t' in x for x in genes+cells)
    with h5py.File(input_path,'r') as handle:
        assert handle['X'].attrs['encoding-type']=='csr_matrix'
        assert tuple(handle['X'].attrs['shape'])==(source.n_obs,source.n_vars)
    (output/'genes.tsv').write_text('\n'.join(genes)+'\n')
    metadata=source.obs.copy();metadata.insert(0,'cell_id',cells)
    metadata.to_csv(output/'cellsmeta.tsv',sep='\t',index=False)
    cfg={'input':str(input_path),'output':str(output),'mode':mode,'object_id':object_id,
         'code_dir':str(CODE),'raw_scores':str(output/'scores_all_methods.h5'),
         'genesets':str(GENESETS),'genesets_sha256':digest(GENESETS),'chunk_size':chunk_size,'workers':workers,
         'checkpoints':str(output/'checkpoints'),'resource_pause_file':str(WORK/'formal_pause_new_batches.flag') if mode=='formal' else None,'seed':2026,'AUCell_chunk_seed':'seed + chunk_index - 1; explicit SerialParam RNGseed','RNA_gene_universe':source.n_vars,
         'cell_count':source.n_obs,'source_sha256':source_hash,'cell_axis_sha256':axis_digest(cells),
         'gene_axis_sha256':axis_digest(genes),'ready_receipt':str(ready) if ready else None,
         'algorithm_version':'RNA_full_universe_vectorized_v2',
         'source_stat':{'bytes':input_path.stat().st_size,'mtime_ns':input_path.stat().st_mtime_ns},
         'pipeline_code_sha256':{name:digest(CODE/name) for name in ('run_pathways.py','run_pathways.R','pathway_core.R')},
         'assay_semantics':'Existing log1p normalized RNA only; no ADT features, filtering, relabeling or renormalization.',
         'created_utc':datetime.now(timezone.utc).isoformat()}
    previous=output/'scoring_config.json'
    if previous.exists():
        old=json.loads(previous.read_text())
        for key in ('source_sha256','genesets_sha256','chunk_size','cell_axis_sha256','gene_axis_sha256','object_id','algorithm_version'):assert old.get(key)==cfg[key],'Resume input/config changed: '+key
    save(previous,cfg);source.file.close();return cfg

def finalize(cfg):
    output=Path(cfg['output']);source=ad.read_h5ad(cfg['input'],backed='r')
    assert {'bytes':Path(cfg['input']).stat().st_size,'mtime_ns':Path(cfg['input']).stat().st_mtime_ns}==cfg['source_stat'],'Source RNA changed during scoring'
    registry_path=Path(cfg['AUCell_zero_mask']['raw_scores']).parent/'pathway_gene_overlap.tsv' if cfg.get('AUCell_zero_mask') else output/'pathway_gene_overlap.tsv'
    registry=pd.read_csv(registry_path,sep='\t')
    selected=registry.loc[registry.included].copy().set_index('pathway_id')
    selected['low_RNA_gene_coverage']=selected.RNA_gene_overlap<10
    selected.to_csv(output/'pathway_gene_coverage.tsv',sep='\t')
    result={'status':'validated','object_id':cfg['object_id'],'mode':cfg['mode'],
            'source_RNA_sha256':cfg['source_sha256'],'genesets_sha256':cfg['genesets_sha256'],
            'cell_count':source.n_obs,'RNA_gene_universe':source.n_vars,
            'cell_axis_sha256':cfg['cell_axis_sha256'],'pathway_count':len(selected),
            'excluded_gene_sets':registry.loc[~registry.included].to_dict(orient='records'),'files':[],
            'scoring_config_sha256':digest(output/'scoring_config.json'),'gene_overlap_registry':str(registry_path),'gene_overlap_registry_sha256':digest(registry_path),
            'pipeline_code_sha256':cfg['pipeline_code_sha256'],
            'low_RNA_gene_coverage_pathway_count':int(selected.low_RNA_gene_coverage.sum()),
            'baseline_AUCell_publication':'Not published; contains random zero-gene rank contribution in AUCell 1.24.0 default ranking',
            'AUCell_zero_mask':cfg.get('AUCell_zero_mask',{})}
    if cfg.get('AUCell_zero_mask'):
        result.update({key:cfg['AUCell_zero_mask'][key] for key in ('resource_gene_alias_table_sha256','KEGG_effective_alias_count','corrected_pathway_count','unchanged_US_columns')})
    if cfg['mode']=='formal':assert cfg.get('AUCell_zero_mask',{}).get('status')=='validated','Formal AUCell zero-mask branch is required'
    with ExitStack() as stack:
        raw=stack.enter_context(h5py.File(cfg['raw_scores'],'r'))
        original=stack.enter_context(h5py.File(cfg['input'],'r'))
        masked=stack.enter_context(h5py.File(cfg['AUCell_zero_mask']['raw_scores'],'r')) if cfg.get('AUCell_zero_mask') else None
        assert np.all(raw['completed_chunks'][...]==1)
        raw_cells=raw['cell_ids'].asstr()[...].tolist();raw_pathways=raw['pathway_ids'].asstr()[...].tolist()
        assert raw_cells==source.obs_names.astype(str).tolist() and raw_pathways==selected.index.tolist()
        for method in METHODS:
            used=masked if masked is not None else raw
            assert np.all(used['completed_chunks'][...]==1)
            assert used['cell_ids'].asstr()[...].tolist()==raw_cells and used['pathway_ids'].asstr()[...].tolist()==raw_pathways
            scores=used['scores/'+method];assert scores.shape==(source.n_obs,len(selected))
            target_dir=output/'RNA_scores_zero_mask_and_alias_v1' if cfg.get('AUCell_zero_mask') else output
            target_dir.mkdir(parents=True,exist_ok=True)
            target=target_dir/f'GSE332623_PBMC_{method}_corrected.h5ad' if cfg['mode']=='formal' else target_dir/f'TEST_{method}.h5ad'
            method_receipt=target_dir/f'{method}_H5AD_completion.json'
            if target.exists() and method_receipt.exists():
                existing=json.loads(method_receipt.read_text());assert existing['sha256']==digest(target)
                assert existing['source_RNA_sha256']==cfg['source_sha256']
                result['files'].append(existing);continue
            temporary=target.with_suffix('.h5ad.tmp')
            template=ad.AnnData(X=np.empty((0,len(selected)),dtype=np.float32),obs=source.obs.iloc[:0].copy(),var=selected)
            template.uns['scaid_scoring']={'method':method,'object_id':cfg['object_id'],'source_RNA_sha256':cfg['source_sha256'],
                'KEGG_gene_sets_sha256':cfg['genesets_sha256'],'seed':2026,'RNA_gene_universe':source.n_vars,
                'max_gene_set_overlap':500,'min_gene_set_overlap':1,'singscore_centerScore':False,
                'AUCell_maxRank':int(np.ceil(.05*source.n_vars)),'UCell_maxRank':min(1500,source.n_vars),
                'chunk_size':cfg['chunk_size'],'AUCell_chunk_seed':cfg['AUCell_chunk_seed'],'scope':'Mixed study cohort; pathway scores do not assign disease/control or donor identity.'}
            if masked is not None:
                template.uns['scaid_scoring']['symbol_alias_mapping_sha256']=cfg['AUCell_zero_mask']['symbol_alias_mapping_sha256']
                template.uns['scaid_scoring']['resource_gene_alias_table_sha256']=cfg['AUCell_zero_mask']['resource_gene_alias_table_sha256']
                template.uns['scaid_scoring']['symbol_alias_scope']='Verified one-to-one source-symbol/Ensembl mappings for internal KEGG membership only; original RNA gene axis unchanged'
            if method=='AUCell' and masked is not None:
                template.uns['scaid_scoring']['AUCell_zero_policy']='Explicit RNA>0 weight mask; positive-gene ranks and fixed 5% cap unchanged'
                template.uns['scaid_scoring']['AUCell_zero_mask_software_workaround']='AUCell 1.24.0 keepZeroesAsNA flag indexing bug'
                template.uns['scaid_scoring']['AUCell_zero_mask_code_sha256']=cfg['AUCell_zero_mask']['code_sha256']['aucell_zero_mask_core.R']
            template.write_h5ad(temporary)
            minimum=float('inf');maximum=float('-inf')
            with h5py.File(temporary,'r+') as final:
                for group in ('obs','obsm'):
                    if group in final:del final[group]
                    original.copy(group,final)
                del final['X'];data=final.create_dataset('X',shape=scores.shape,dtype='float32',chunks=(min(1000,source.n_obs),len(selected)),compression='gzip',compression_opts=4)
                data.attrs['encoding-type']='array';data.attrs['encoding-version']='0.2.0'
                for start in range(0,source.n_obs,1000):
                    values=np.asarray(scores[start:start+1000,:],dtype=np.float32)
                    assert np.isfinite(values).all();data[start:start+len(values),:]=values
                    minimum=min(minimum,float(values.min()));maximum=max(maximum,float(values.max()))
                final.flush()
                for start in range(0,source.n_obs,1000):assert np.array_equal(data[start:start+1000,:],scores[start:start+1000,:])
                null_nodes=[]
                final.visititems(lambda name,item:null_nodes.append(name) if item.attrs.get('encoding-type')=='null' else None)
                assert not null_nodes,'Unsupported null encoding in score H5AD: '+str(null_nodes)
            check=ad.read_h5ad(temporary,backed='r')
            assert check.shape==(source.n_obs,len(selected)) and check.obs.equals(source.obs)
            assert check.obs_names.equals(source.obs_names) and check.var_names.tolist()==raw_pathways
            for key in source.obsm:assert np.array_equal(check.obsm[key],source.obsm[key])
            check.file.close();production_read=production_reader_check(temporary,(source.n_obs,len(selected)))
            os.replace(temporary,target)
            record={'method':method,'path':str(target),'shape':[source.n_obs,len(selected)],
                'bytes':target.stat().st_size,'sha256':digest(target),'min':minimum,'max':maximum,
                'source_RNA_sha256':cfg['source_sha256'],
                'full_score_component_readback':True,'all_cells_metadata_embeddings_preserved':True,'production_backed_reader':production_read,'no_null_encoding':True,'validated':True}
            save(method_receipt,record);result['files'].append(record)
    source.file.close();save(output/'pathway_completion.json',result);return result

def ensure_aucell_zero_mask(cfg,env):
    """Recompute masked AUCell and the approved alias-affected UCell/singscore columns."""
    output=Path(cfg['output']);branch=output/'RNA_scores_zero_mask_and_alias_v1';branch.mkdir(parents=True,exist_ok=True)
    approval_path=WORK/'KEGG_alias_mapping_approved.json'
    approval=json.loads(approval_path.read_text());assert approval['status']=='approved'
    official_validation=Path(approval['actual_pipeline_validation'])
    assert digest(official_validation)==approval['actual_pipeline_validation_sha256']
    verified=json.loads(official_validation.read_text())
    assert verified['status']=='validated' and all(x==0 for x in verified['official_float32_max_absolute_errors'].values())
    alias_path=Path(approval['path']);alias_sha=approval['sha256'];alias_count=approval['mapping_count']
    verify_symbol_aliases(cfg['input'],alias_path,alias_sha,alias_count)
    cfg['symbol_alias_path']=str(alias_path);cfg['symbol_alias_sha256']=alias_sha
    cfg['symbol_alias_approval']=str(approval_path);cfg['symbol_alias_approval_sha256']=digest(approval_path)
    validation=WORK/'testing/independent_AUCell_zero_mask_validation.json'
    reference=json.loads(validation.read_text())
    assert reference['status']=='validated' and reference['explicit_RNA_zero_to_NA_official_calcAUC_max_absolute_error']==0
    assert reference['zero_detection_nonzero_scores']==0 and reference['cells']==1000 and reference['pathways']==356
    branch_cfg=branch/'zero_mask_config.json'
    code_hash={name:digest(CODE/name) for name in ('aucell_zero_mask_core.R','run_aucell_zero_mask.R','symbol_alias_view.R','run_pathways.py')}
    branch_config={**cfg,'AUCell_zero_mask_code_sha256':code_hash,'AUCell_zero_mask_version':'RNA_scores_zero_mask_and_alias_v1'}
    if branch_cfg.exists():
        previous=json.loads(branch_cfg.read_text())
        for key in ('source_sha256','cell_axis_sha256','gene_axis_sha256','genesets_sha256','chunk_size','seed'):assert previous[key]==branch_config[key]
        assert previous['AUCell_zero_mask_code_sha256']['aucell_zero_mask_core.R']==code_hash['aucell_zero_mask_core.R']
    else:save(branch_cfg,branch_config)
    with (branch/'zero_mask_run.log').open('a') as log:
        subprocess.run(['Rscript','--vanilla',str(CODE/'run_aucell_zero_mask.R'),str(branch_cfg),str(branch)],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,preexec_fn=boundary)
    answer=json.loads((branch/'r_zero_mask_completion.json').read_text())
    assert answer['status']=='all_chunks_complete' and answer['zero_detection_nonzero_scores']==0
    raw=branch/'scores_three_methods_zero_mask_and_alias.h5'
    affected=pd.read_csv(branch/'alias_recomputed_pathways.tsv',sep='\t').pathway_id.tolist()
    with h5py.File(raw,'r') as fixed,h5py.File(cfg['raw_scores'],'r') as baseline:
        ids=fixed['pathway_ids'].asstr()[...].tolist();unaffected=[i for i,x in enumerate(ids) if x not in affected]
        assert len(affected)==approval['affected_included_pathways'] and len(unaffected)==356-len(affected)
        for method in ('UCell','singscore'):
            for start in range(0,cfg['cell_count'],1000):
                assert np.array_equal(fixed['scores/'+method][start:start+1000,:][:,unaffected],baseline['scores/'+method][start:start+1000,:][:,unaffected]),'Unchanged pathway columns differ'
        official_path=Path(approval['official_reference']);assert digest(official_path)==approval['official_reference_sha256']
        with h5py.File(official_path,'r') as official:
            assert fixed['cell_ids'].asstr()[:1000].tolist()==official['cell_ids'].asstr()[...].tolist()
            assert fixed['pathway_ids'].asstr()[...].tolist()==official['pathway_ids'].asstr()[...].tolist()
            for method in METHODS:assert np.array_equal(fixed['scores/'+method][:1000,:],official[method][...].astype('float32')),'Final first-block official equality failed: '+method
    answer.update({'status':'validated','raw_scores':str(raw),'raw_scores_sha256':digest(raw),'code_sha256':code_hash,
        'baseline_raw_scores':cfg['raw_scores'],'baseline_AUCell_not_published':True,
        'official_equivalence_validation':str(validation),'official_equivalence_validation_sha256':digest(validation),
        'official_equivalence_max_absolute_error':0,'unchanged_UCell_singscore_columns':len(unaffected),'recomputed_UCell_singscore_columns':affected,
        'unaffected_columns_bitwise_equal_baseline':True,'symbol_alias_mapping':str(alias_path),'symbol_alias_mapping_sha256':alias_sha,
        'symbol_alias_mapping_count':alias_count,'symbol_alias_approval_sha256':digest(approval_path),
        'resource_gene_alias_table_sha256':alias_sha,'KEGG_effective_alias_count':approval['KEGG_effective_alias_count'],
        'corrected_pathway_count':len(affected),'unchanged_US_columns':len(unaffected),
        'official_full_1000x356_reference':str(official_path),'official_full_1000x356_reference_sha256':approval['official_reference_sha256'],
        'official_full_first_block_float32_equality':True,'independent_alias_gate_sha256':approval['independent_review_sha256'],
        'source_RNA_sha256':cfg['source_sha256'],'cell_axis_sha256':cfg['cell_axis_sha256']})
    save(branch/'zero_mask_completion.json',answer);cfg['AUCell_zero_mask']=answer
    save(output/'finalization_config.json',cfg)
    return cfg

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',required=True,type=Path);parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--mode',choices=['test','pilot','formal'],required=True);parser.add_argument('--ready',type=Path)
    parser.add_argument('--object-id',default='Atopic_Dermatitis/GSE332623/PBMC');parser.add_argument('--chunk-size',type=int,default=1000);parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    resource_parameters=WORK/'formal_resources_parameters.json'
    if args.mode=='formal' and resource_parameters.exists():
        parameters=json.loads(resource_parameters.read_text());assert parameters['score_workers'] in (4,8)
        args.workers=parameters['score_workers']
    if (args.output/'pathway_completion.json').is_file():
        existing=json.loads((args.output/'pathway_completion.json').read_text())
        assert existing['status']=='validated' and existing['source_RNA_sha256']==digest(args.input)
        assert existing['object_id']==args.object_id and existing['mode']==args.mode
        for file in existing['files']:assert digest(file['path'])==file['sha256']
        print(json.dumps({'status':'already_validated','mode':args.mode,'immutable_H5ADs_unchanged':True}));return
    baseline_config=args.output/'scoring_config.json'
    if args.mode=='formal' and baseline_config.exists() and not (args.output/'baseline_scoring_config.json').exists():
        save(args.output/'baseline_scoring_config.json',json.loads(baseline_config.read_text()))
    cfg=prepare(args.input,args.output,args.mode,args.ready,args.chunk_size,args.object_id,args.workers)
    env=os.environ.copy();env.update({x:'1' for x in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
    baseline_done=args.output/'r_scoring_completion.json'
    if baseline_done.exists():
        assert json.loads(baseline_done.read_text())['status']=='all_chunks_complete'
        with h5py.File(cfg['raw_scores'],'r') as baseline:assert np.all(baseline['completed_chunks'][...]==1)
    else:
        with (args.output/'scoring_run.log').open('a') as log:
            subprocess.run(['Rscript','--vanilla',str(CODE/'run_pathways.R'),str(args.output/'scoring_config.json')],
                           stdout=log,stderr=subprocess.STDOUT,env=env,check=True,preexec_fn=boundary)
    if args.mode=='formal':cfg=ensure_aucell_zero_mask(cfg,env)
    result=finalize(cfg);print(json.dumps({'status':result['status'],'cells':result['cell_count'],'pathways':result['pathway_count'],'files':len(result['files']),'mode':result['mode']}))

if __name__=='__main__':main()

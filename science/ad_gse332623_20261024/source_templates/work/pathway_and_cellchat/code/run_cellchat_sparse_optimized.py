#!/usr/bin/env python3
"""Isolated single-worker sparse-exact CellChat candidate; no active analysis is changed."""
from pathlib import Path
import argparse,json,os,subprocess,resource,gc
import psutil
import anndata as ad
import pandas as pd
from PIL import Image
from run_pathways import CODE,CORRECTED,digest,axis_digest,save,declared_hashes,boundary,SYMBOL_ALIASES,SYMBOL_ALIAS_SHA,verify_symbol_aliases

def finalize(output,cfg):
    output=Path(output)
    manifest=pd.read_csv(output/'cellchat_figure_manifest.tsv',sep='\t')
    assert not manifest.path.duplicated().any()
    for index,row in manifest.iterrows():
        path=Path(row.path)
        with Image.open(path) as image:
            assert image.size==(row.width_px,row.height_px);image.verify()
        manifest.loc[index,'sha256']=digest(path);manifest.loc[index,'bytes']=path.stat().st_size
    manifest['bytes']=manifest['bytes'].astype('int64')
    manifest.to_csv(output/'cellchat_figure_manifest.tsv',sep='\t',index=False)
    files=[]
    for path in sorted(output.iterdir()):
        if path.suffix in ('.rds','.tsv','.gz','.json') and path.name!='cellchat_completion.json':files.append({'path':str(path),'bytes':path.stat().st_size,'sha256':digest(path)})
    computed=json.loads((output/'cellchat_computation.json').read_text());rendered=json.loads((output/'cellchat_rendering.json').read_text())
    assert rendered['errors']==0 and (cfg['mode']!='formal' or len(rendered['families'])==9)
    assert manifest.object_id.eq(cfg['object_id']).all()
    result={'status':'validated','mode':cfg['mode'],'object_id':cfg['object_id'],'source_RNA_sha256':cfg['source_RNA_sha256'],'coarse_review_sha256':cfg.get('coarse_review_sha256'),'symbol_alias_sha256':cfg.get('symbol_alias_sha256'),'computation':computed,'rendering':rendered,'files':files,
            'all_pngs_parsed_and_sha256_recorded':True,'figure_manifest':str(output/'cellchat_figure_manifest.tsv')}
    save(output/'cellchat_completion.json',result);return result

def optimized_boundary():
    os.sched_setaffinity(0,{10})
    _,hard=resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS,(128*1024**3,hard))

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=['test','pilot','formal'],required=True);p.add_argument('--ready',type=Path);p.add_argument('--nboot',type=int,default=100)
    p.add_argument('--object-id',default='Atopic_Dermatitis/GSE332623/PBMC');p.add_argument('--plot-limit',type=int,default=0)
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--cache-receipt',type=Path,required=True)
    p.add_argument('--reference-pid',type=int)
    p.add_argument('--implementation-review',type=Path,required=True)
    p.add_argument('--full-axis-review',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    assert a.workers==1,'Isolated optimized candidate has exactly one worker'
    assert a.output.resolve()==(CORRECTED/'pathway_and_cellchat/cellchat_sparse_exact_candidate').resolve(),'Never overwrite the active official reference directory'
    for proc in psutil.process_iter():
        try:
            cmd=proc.cmdline()
            ours=any(str(CODE/name) in arg for name in ('run_cellchat.R','run_cellchat_optimized.R','run_cellchat_sparse_optimized.R') for arg in cmd)
            if proc.uids().real==os.getuid() and proc.status()!=psutil.STATUS_ZOMBIE and ours:
                binding=json.loads((CODE.parent/'CellChat_reference_coexistence_binding.json').read_text())
                assert a.reference_pid==proc.pid==binding['pid'],'Refuse unapproved concurrent CellChat R'
                assert proc.uids().real==binding['uid'] and proc.create_time()==binding['created']
                assert cmd==binding['cmdline'] and digest(binding['config_path'])==binding['config_sha256']
                reference_cfg=json.loads(Path(binding['config_path']).read_text())
                assert reference_cfg['source_RNA_sha256']==binding['source_RNA_sha256']
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    cache_receipt=json.loads(a.cache_receipt.read_text())
    assert cache_receipt['status']=='validated_internal_cache_not_a_model_run',cache_receipt['status']
    equivalence=CODE.parent.parent/'analysis/validation/sparse_triMean_prototype/pilot100.json'
    eq=json.loads(equivalence.read_text())
    assert eq['bootstrap']==100 and eq['probability_max_absolute_error']==0 and eq['pvalue_max_absolute_error']==0 and eq['all_probability_pvalue_and_dimnames_bitwise_identical']
    implementation_review=json.loads(a.implementation_review.read_text())
    full_axis_review=json.loads(a.full_axis_review.read_text())
    assert implementation_review['status']=='passed_implementation_review_not_formal_network'
    assert full_axis_review['status']=='validated_formal_axis_means_only_not_network'
    full_binding=implementation_review['full_axis_review']
    assert Path(full_binding['path']).resolve()==a.full_axis_review.resolve()
    assert full_binding['sha256']==digest(a.full_axis_review) and full_binding['status']==full_axis_review['status']
    independent_binding=implementation_review['independent_review']
    independent_review=json.loads(Path(independent_binding['path']).read_text())
    assert independent_binding['sha256']==digest(independent_binding['path'])
    assert independent_binding['status']==independent_review['status']=='passed_independent_sparse_helper_review_not_formal_network'
    cache_path=a.cache_receipt.parent/'normalized_DB_expression_cache.h5'
    assert cache_path.is_file()

    sha=digest(a.input);assert a.nboot>=2 and 1<=a.workers<=4
    assert sha==implementation_review['source_RNA_sha256']==full_axis_review['source_RNA_sha256']
    assert full_axis_review['source_cache_sha256']==cache_receipt['sha256']==digest(cache_path)
    assert implementation_review['cache_review']['cache_sha256']==cache_receipt['sha256']
    assert implementation_review['cache_review']['receipt_sha256']==digest(a.cache_receipt)
    assert full_axis_review['cells']==528739 and full_axis_review['genes']==1434 and full_axis_review['groups']==8
    assert full_axis_review['seed']==2026 and full_axis_review['permutations_checked']==2
    assert len(full_axis_review['comparisons'])==3 and all(x['bitwise_all_values_and_attributes_identical'] is True for x in full_axis_review['comparisons'])
    if a.reference_pid:
        binding=json.loads((CODE.parent/'CellChat_reference_coexistence_binding.json').read_text())
        assert binding['pid']==a.reference_pid and binding['source_RNA_sha256']==sha
    assert sha in set(declared_hashes(cache_receipt)),'Cache receipt must bind immutable source RNA'
    assert digest(cache_path) in set(declared_hashes(cache_receipt)),'Cache bytes must match independently validated receipt'
    if (a.output/'cellchat_completion.json').is_file():
        existing=json.loads((a.output/'cellchat_completion.json').read_text())
        assert existing['status']=='validated' and existing['source_RNA_sha256']==sha and existing['mode']==a.mode
        assert existing['object_id']==a.object_id and existing['computation']['nboot']==a.nboot
        if a.mode=='formal':
            assert existing['coarse_review_sha256']==digest(CORRECTED/'integration/annotation_evidence/coarse_annotation_review.tsv')
            assert existing['symbol_alias_sha256']==SYMBOL_ALIAS_SHA
        for file in existing['files']:assert digest(file['path'])==file['sha256']
        print(json.dumps({'status':'already_validated','mode':a.mode,'immutable_results_unchanged':True}));return
    if a.mode=='formal':
        assert a.output.resolve().is_relative_to(CORRECTED) and a.nboot==100 and a.plot_limit==0
        receipt=json.loads(a.ready.read_text());assert receipt['status']=='validated'
        assert receipt.get('n_libraries')==27 and Path(receipt['RNA_path']).resolve()==a.input.resolve()
        assert sha in set(declared_hashes(receipt))
        assert a.object_id=='Atopic_Dermatitis/GSE332623/PBMC'
    source=ad.read_h5ad(a.input,backed='r');assert source.obs_names.is_unique and source.var_names.is_unique
    assert 'MajorCellType' in source.obs
    metadata=source.obs.copy();metadata.insert(0,'cell_id',source.obs_names.astype(str))
    metadata.to_csv(a.output/'cellsmeta.tsv',sep='\t',index=False)
    cfg={'input':str(a.input.resolve()),'output':str(a.output.resolve()),'code_dir':str(CODE),'mode':a.mode,'object_id':a.object_id,
         'source_RNA_sha256':sha,'cell_axis_sha256':axis_digest(source.obs_names),'gene_axis_sha256':axis_digest(source.var_names),
         'pipeline_code_sha256':{name:digest(CODE/name) for name in ('run_cellchat_sparse_optimized.py','run_cellchat_sparse_optimized.R','cellchat_exact_mean.R','render_cellchat.R','pathway_core.R','symbol_alias_view.R')},
         'min_cells':10,'chunk_size':1000,'nboot':a.nboot,'seed':2026,'plot_limit':a.plot_limit,'workers':a.workers,
         'DB_expression_cache':str(cache_path),'DB_expression_cache_sha256':digest(cache_path),
         'DB_expression_cache_receipt_sha256':digest(a.cache_receipt),'implementation_equivalence_sha256':digest(equivalence),
         'implementation_review_sha256':digest(a.implementation_review),'full_axis_aggregation_review_sha256':digest(a.full_axis_review),
         'sparse_helper_R':str(CODE/'sparse_trimean_prototype.R'),
         'sparse_helper_cpp':str(CODE/'sparse_trimean_prototype.cpp'),
         'installed_official_function_source_sha256':digest(CODE.parent/'testing/computeCommunProb_installed_source.R')}
    for name in ('sparse_trimean_prototype.R','sparse_trimean_prototype.cpp'):
        cfg['pipeline_code_sha256'][name]=digest(CODE/name)
    for name,sha_reviewed in implementation_review['reviewed_code_sha256'].items():
        assert cfg['pipeline_code_sha256'][name]==sha_reviewed,'Current sparse implementation must match independent review'
    cfg['independent_implementation_review_sha256']=independent_binding['sha256']
    if a.mode=='formal':
        verify_symbol_aliases(a.input)
        cfg['symbol_alias_path']=str(SYMBOL_ALIASES);cfg['symbol_alias_sha256']=SYMBOL_ALIAS_SHA
        review=CORRECTED/'integration/annotation_evidence/coarse_annotation_review.tsv'
        assert review.is_file()
        cfg['coarse_review_path']=str(review);cfg['coarse_review_sha256']=digest(review)
        review_table=pd.read_csv(review,sep='\t')
        types=review_table.loc[review_table.review_status.str.contains('source_dominated',case=False,regex=True),'MajorCellType'].tolist()
        predicted=source.obs.MajorCellType.astype(str).isin(types)
        save(a.output/'provisional_exclusion_receipt.json',{'status':'verified_scope','RNA_sha256':sha,'coarse_review_sha256':cfg['coarse_review_sha256'],'excluded_labels':types,'excluded_cells':int(predicted.sum()),'all_RNA_cells_retained':source.n_obs,'reason':'source_dominated_provisional_nonimmune_prediction'})
    source.file.close();del source,metadata;gc.collect();save(a.output/'cellchat_config.json',cfg)
    env=os.environ.copy();env.update({x:'1' for x in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')})
    with (a.output/'cellchat_run.log').open('a') as log:
        subprocess.run(['Rscript','--vanilla',str(CODE/'run_cellchat_sparse_optimized.R'),str(a.output/'cellchat_config.json')],stdout=log,stderr=subprocess.STDOUT,env=env,check=True,preexec_fn=optimized_boundary)
    result=finalize(a.output,cfg)
    print(json.dumps({'status':result['status'],'mode':a.mode,'modeled_cells':result['computation']['modeled_cells'],'figures':result['rendering']['figures']}))

if __name__=='__main__':main()

#!/usr/bin/env python3
"""Materialize recorded source, verify inputs, and replay frozen-RNA kernels.

No production service, database, or original data is modified. Scientific R
kernels are the recorded sources; this file supplies portable local bindings.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, os, shutil, subprocess, sys

HERE = Path(__file__).resolve().parent
RNA_SHA = 'b9aa963123cb8f010d3e092e03507b70d02170afd93d32a6cdea2efe2a2c09ec'
OBJECT = 'Atopic_Dermatitis/GSE332623/PBMC'

def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def read(path): return json.loads(Path(path).read_text())
def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')

def materialize(root, raw=None, r_lib=None, web_python=None, models=None):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    replacements={'@PROJECT@':str(root/'project'),'@WORK@':str(root/'work'),
      '@RAW_INPUT@':str(Path(raw).resolve() if raw else root/'raw'),
      '@R_LIB@':str(Path(r_lib).resolve() if r_lib else root/'R-library'),
      '@WEB_PYTHON@':str(Path(web_python).resolve() if web_python else Path(sys.executable)),
      '@ANALYSIS_PYTHON@':sys.executable,
      '@CELLTYPIST_MODELS@':str(Path(models).resolve() if models else root/'celltypist-models'),
      '@DOWNLOADS@':str(root/'downloads'),'@GENE_CACHE@':str(root/'gene_cache'),
      '@FIGURES@':str(root/'figures')}
    files=[]
    for source in sorted((HERE/'source_templates').rglob('*')):
        if not source.is_file(): continue
        target=root/source.relative_to(HERE/'source_templates');target.parent.mkdir(parents=True,exist_ok=True)
        text=source.read_text()
        for token,value in replacements.items():text=text.replace(token,value)
        if target.exists() and target.read_text()!=text:raise ValueError(f'Refuse to replace changed materialized source: {target}')
        target.write_text(text);files.append({'path':str(target.relative_to(root)),'sha256':digest(target)})
    resource=root/'resources';resource.mkdir(exist_ok=True)
    for source in (HERE/'resources').iterdir():
        if not source.is_file():continue
        text=source.read_text()
        for token,value in replacements.items():text=text.replace(token,value)
        (resource/source.name).write_text(text)
    inventory=root/'project/code/v2/00.Reference/SampleInventory.tsv'
    inventory.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(resource/'SampleInventory.GSE332623.tsv',inventory)
    manifest=root/'project/00.corrected_analysis/GSE332623/capture_manifest.tsv'
    manifest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(resource/'capture_manifest.tsv',manifest)
    dispatch_code=root/'project/code/GSE332623_analysis'
    dispatch_code.mkdir(parents=True,exist_ok=True)
    for name in ['analysis_pipeline.py','remap_and_recheck_reference.py','extend_coarse_marker_review.py','render_coarse_marker_summary.py','run_when_QC_ready.py']:
        shutil.copyfile(root/'work/analysis/code'/name,dispatch_code/name)
    # Keep historical evidence separate: never rewrite its original hash bindings.
    save(root/'materialized_sources.json',{'status':'path_bindings_materialized_not_an_analysis_run','files':files,'bindings':replacements})
    return root

def load_module(path, name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def execute(argv):
    print('+ '+ ' '.join(map(str,argv)),flush=True)
    env=os.environ.copy();env.update({name:'1' for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']})
    subprocess.run(list(map(str,argv)),env=env,check=True)

def verify_rna(path):
    actual=digest(path)
    if actual!=RNA_SHA:raise ValueError('Frozen RNA SHA256 mismatch; this replay accepts the published 2026-10-24 RNA only')
    return str(Path(path).resolve())

def verify_R_runtime(kind, output):
    observed=read(HERE/'environment-observed.json')
    packages=observed['R_pathways'] if kind=='scores' else observed['R_CellChat']
    packages={**packages,'Matrix':'1.6.5','hdf5r':'1.3.12'}
    versions=','.join(json.dumps(name)+'='+json.dumps(version) for name,version in packages.items())
    program='stopifnot(as.character(getRversion())=="4.3.3"); expected<-c('+versions+'); for(p in names(expected))stopifnot(as.character(packageVersion(p))==expected[[p]])'
    if kind=='cellchat':
        destination=Path(output)/'installed_CellChat_function_canonical.txt'
        program+='; writeLines(deparse(get("computeCommunProb",asNamespace("CellChat")),width.cutoff=500L),commandArgs(TRUE)[1],useBytes=TRUE)'
        execute(['Rscript','--vanilla','-e',program,destination])
        expected=read(HERE/'resources/CellChat_official_function_fingerprint.json')
        if digest(destination)!=expected['canonical_function_sha256']:
            raise ValueError('Installed CellChat scientific function differs from the recorded run; a new implementation review is required')
    else:execute(['Rscript','--vanilla','-e',program])

def scores(args):
    rna=verify_rna(args.rna);out=Path(args.output).resolve()
    if out.exists():raise ValueError('Replay output must not already exist')
    out.mkdir(parents=True);root=materialize(out/'runtime')
    verify_R_runtime('scores',out)
    code=root/'work/pathway_and_cellchat/code';genesets=Path(args.genesets).resolve()
    original=read(HERE/'evidence/zero_mask_config.json')
    if digest(genesets)!=original['genesets_sha256']:raise ValueError('Exact 359-set KEGG RDS checksum mismatch; consult INPUTS.md')
    (root/'project/20240526_kegg_hsa_list.rds').symlink_to(genesets)
    cfg={'input':rna,'output':str(out),'code_dir':str(code),'raw_scores':str(out/'baseline/scores_all_methods.h5'),
      'genesets':str(genesets),'checkpoints':str(out/'baseline/checkpoints'),'chunk_size':1000,
      'workers':min(args.workers,4),'seed':2026,'resource_pause_file':None}
    save(out/'baseline_config.json',cfg)
    execute(['Rscript','--vanilla',code/'run_pathways.R',out/'baseline_config.json'])
    cfg['workers']=args.workers
    cfg['symbol_alias_path']=str(root/'resources/approved_symbol_renames.tsv')
    save(out/'corrected_config.json',cfg)
    execute(['Rscript','--vanilla',code/'run_aucell_zero_mask.R',out/'corrected_config.json',out/'RNA_scores_zero_mask_and_alias_v1'])
    # Independent official first-block methods, identical to the recorded review.
    execute(['Rscript','--vanilla',root/'work/analysis/code/official_versioned_crosswalk_1000.R',code,rna,
      root/'resources/approved_symbol_renames.tsv',root/'resources/KEGG_set_impact.tsv',out/'official_first1000'])
    import h5py, numpy as np
    with h5py.File(out/'RNA_scores_zero_mask_and_alias_v1/scores_three_methods_zero_mask_and_alias.h5') as fast, h5py.File(out/'official_first1000/official_1000x356_unified_crosswalk.h5') as official:
        assert fast['cell_ids'].asstr()[:1000].tolist()==official['cell_ids'].asstr()[:].tolist()
        assert fast['pathway_ids'].asstr()[:].tolist()==official['pathway_ids'].asstr()[:].tolist()
        errors={}
        for method in ['AUCell','UCell','singscore']:
            a=fast['scores/'+method][:1000,:];b=official[method][:].astype('float32')
            assert np.array_equal(a,b),method+' first-block mismatch'
            errors[method]=float(np.max(np.abs(a-b)))
    save(out/'replay_receipt.json',{'status':'computed_and_first_block_officially_verified','input_RNA_sha256':RNA_SHA,
      'kernels':'recorded baseline plus final explicit-zero-mask/resource-alias branch',
      'official_first1000_float32_max_absolute_errors':errors,'all_cells_full_recomputation_claimed':True,
      'H5AD_export_recreated':False,'byte_identical_original_files_claimed':False,'workers':args.workers})

def cellchat(args):
    """Rebuild source CSC/DB cache, then call the unchanged selected R kernel."""
    rna=verify_rna(args.rna);out=Path(args.output).resolve()
    if out.exists():raise ValueError('Replay output must not already exist')
    out.mkdir(parents=True);root=materialize(out/'runtime')
    verify_R_runtime('cellchat',out)
    work=root/'work';code=work/'pathway_and_cellchat/code';base=root/'project/00.corrected_analysis/GSE332623'
    target=out/'cellchat_sparse_exact_replay';target.mkdir()
    import anndata as ad
    import numpy as np
    source=ad.read_h5ad(rna,backed='r');meta=source.obs.copy();meta.insert(0,'cell_id',source.obs_names.astype(str))
    meta.to_csv(target/'cellsmeta.tsv',sep='\t',index=False);source.file.close()
    cfg=read(HERE/'evidence/cellchat_config.json')
    cfg.update(input=rna,output=str(target),code_dir=str(code),workers=1,
      symbol_alias_path=str(root/'resources/cellchat_symbol_aliases_verified.tsv'),
      coarse_review_path=str(root/'resources/coarse_annotation_review.tsv'),
      sparse_helper_R=str(code/'sparse_trimean_prototype.R'),sparse_helper_cpp=str(code/'sparse_trimean_prototype.cpp'))
    # Prepare the same DB matching view and inclusion mask using the original
    # prefix of the selected R script. No probability model is run in this step.
    original=(code/'run_cellchat_sparse_optimized.R').read_text()
    marker='## Internal performance cache:'
    if original.count(marker)!=1:raise ValueError('Selected R preprocessing boundary changed')
    prep=out/'prepare_DB_view.R';prep.write_text(original.split(marker)[0]+"\ncat('DB_VIEW_PREPARED\\n')\n")
    save(out/'cellchat_config.json',cfg);execute(['Rscript','--vanilla',prep,out/'cellchat_config.json'])
    expected_db='45536ae5df0ff6670a12ac97b3610a7455ef4ccdf8667346435131153a40835e'
    if digest(target/'CellChatDB_human_used.rds')!=expected_db:
        raise ValueError('Installed CellChatDB bytes differ from the recorded run; a new resource/implementation review is required')
    index=load_module(work/'website/gene_index_builder_large_cohort.py','ad_replay_index')
    cache_root=work/'pathway_and_cellchat/performance_cache';cache_root.mkdir(parents=True,exist_ok=True)
    entry=index.build_index(rna,cache_root/'full_RNA_csc.h5',OBJECT,RNA_SHA)
    save(work/'website/prepared/gene_entry.json',entry)
    cache=load_module(cache_root/'build_normalized_DB_expression_cache.py','ad_replay_cache')
    cache.INDEX_SHA=entry['index_sha256'];cache.VIEW=target/'CellChat_RNA_symbol_matching_view.tsv';cache.INCLUSION=target/'cell_inclusion.tsv.gz'
    cache.ALIAS=root/'resources/cellchat_symbol_aliases_verified.tsv'
    cache.main()  # Includes every-column/full-axis readback and direct RNA probes.
    receipt=read(cache.RECEIPT)
    cfg.update(DB_expression_cache=str(cache.OUT),DB_expression_cache_sha256=receipt['sha256'],
      DB_expression_cache_receipt_sha256=digest(cache.RECEIPT))
    # Re-run the full-axis observed/two-permutation equivalence gate before modeling.
    execute(['Rscript','--vanilla',work/'analysis/code/validate_sparse_trimean_formal_axis.R',work,out/'equivalence/formal_axis'])
    save(out/'cellchat_config.json',cfg)
    execute(['Rscript','--vanilla',code/'run_cellchat_sparse_optimized.R',out/'cellchat_config.json'])
    computation=read(target/'cellchat_computation.json');rendering=read(target/'cellchat_rendering.json')
    assert computation['nboot']==100 and computation['modeled_cells']==528739
    assert rendering['errors']==0
    save(out/'replay_receipt.json',{'status':'computed_with_source_cache_and_full_axis_aggregation_checks',
      'input_RNA_sha256':RNA_SHA,'cache_receipt_sha256':digest(cache.RECEIPT),
      'original_scientific_kernel_sha256':digest(code/'run_cellchat_sparse_optimized.R'),
      'historical_pilot100_equivalence':'evidence/CellChat_pilot_100bootstrap_equivalence.json',
      'original_full_formal_network_equivalence_claimed':False,'computation':computation,'rendering':rendering,
      'runtime_adaptations':['local paths','regenerated CSC cache hash','no transient coexisting reference PID'],
      'byte_identical_original_container_files_claimed':False})

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    m=sub.add_parser('materialize');m.add_argument('--root',type=Path,required=True)
    m.add_argument('--raw-input',type=Path);m.add_argument('--r-lib',type=Path);m.add_argument('--web-python',type=Path);m.add_argument('--celltypist-models',type=Path)
    v=sub.add_parser('verify-inputs');v.add_argument('--raw-input',type=Path);v.add_argument('--rna',type=Path)
    sub.add_parser('plan')
    for name in ['scores','cellchat']:
        c=sub.add_parser(name);c.add_argument('--rna',type=Path,required=True);c.add_argument('--output',type=Path,required=True)
        if name=='scores':c.add_argument('--genesets',type=Path,required=True);c.add_argument('--workers',type=int,choices=range(1,9),default=8)
    a=p.parse_args()
    if a.command=='materialize':print(materialize(a.root,a.raw_input,a.r_lib,a.web_python,a.celltypist_models))
    elif a.command=='verify-inputs':
        if not a.raw_input and not a.rna:p.error('Provide --raw-input and/or --rna')
        if a.rna:verify_rna(a.rna);print('Frozen RNA full SHA256 matched')
        if a.raw_input:
            files=read(HERE/'resources/raw_inputs.json')['files']
            for f in files:
                if digest(a.raw_input/f['filename'])!=f['sha256']:raise ValueError('MEX checksum mismatch: '+f['filename'])
            print(f'All {len(files)} raw MEX members matched')
    elif a.command=='plan':
        print('1. Read INPUTS.md and ENVIRONMENT_AND_PARAMETERS.md.\n2. Verify the 81 MEX members or frozen RNA.\n3. Materialize raw/QC/integration sources in a new workspace.\n4. Supply pinned external KEGG resource and CellTypist models.\n5. Run scores/cellchat on the exact frozen RNA, or follow raw-stage commands in REPRODUCE.md.\nHistorical evidence is not a fresh-run approval receipt; new raw runs need fresh source and scientific validation.')
    elif a.command=='scores':scores(a)
    else:cellchat(a)

if __name__=='__main__':main()

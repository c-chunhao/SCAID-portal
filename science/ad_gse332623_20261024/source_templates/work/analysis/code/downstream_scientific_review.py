#!/usr/bin/env python3
"""Small independent read-only scientific checks; no scoring or CC rerun."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, re, subprocess, os
import numpy as np
import pandas as pd
import h5py
from scipy.stats import hypergeom

ROOT = Path('@PROJECT@/00.corrected_analysis/GSE332623')
WORK = Path('@WORK@')
OUT = WORK/'analysis/validation'
DS = ROOT/'pathway_and_cellchat'
DW = WORK/'pathway_and_cellchat'
OUT.mkdir(parents=True, exist_ok=True)

def sha(path):
    with Path(path).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def load(path): return json.loads(Path(path).read_text())
def strings(value):
    if isinstance(value,str):yield value
    elif isinstance(value,dict):
        for item in value.values():yield from strings(item)
    elif isinstance(value,list):
        for item in value:yield from strings(item)
def axis(handle, group):
    node=handle[group]; return node[node.attrs['_index']].asstr()[:]
def bh(p):
    order=np.argsort(p); adj=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    result=np.empty_like(adj);result[order]=np.minimum(adj,1);return result

ready=load(ROOT/'integration/analysis_ready.json')
parser=argparse.ArgumentParser()
parser.add_argument('--cellchat-dir',type=Path)
parser.add_argument('--selection',type=Path,default=DS/'selected_cellchat_result.json')
args=parser.parse_args()
selection=None
if args.cellchat_dir is not None:
    CS=args.cellchat_dir.resolve()
    selection={'status':'explicit_directory_argument','cellchat_dir':str(CS),'source_RNA_sha256':ready['RNA_sha256']}
elif args.selection.is_file():
    selection=load(args.selection)
    assert selection['status']=='selected_for_independent_review'
    assert selection['source_RNA_sha256']==ready['RNA_sha256']
    CS=Path(selection['cellchat_dir']).resolve()
else: CS=DS/'cellchat'
assert CS in {DS/'cellchat',DS/'cellchat_exact_candidate',DS/'cellchat_sparse_exact_candidate'},'Only this study reference or isolated reviewed exact candidate may be reviewed'
report={'status':'pending_formal_products','reviewed_utc':datetime.now(timezone.utc).isoformat(),
 'scope':'Independent read-only code, scientific-scope and small numerical checks; no full-matrix repeat or clinical inference.',
 'review_script_sha256':sha(__file__),
 'source_RNA_sha256':ready['RNA_sha256'],'source_ready_sha256':sha(ROOT/'integration/analysis_ready.json'),
 'cellchat_selection':selection,
 'clinical_scope':'Mixed AD/HC source study; cell-level donor and condition unresolved. No disease/control test or independent cohort replication.',
 'sources':[
  'https://github.com/jinworks/CellChat/blob/main/R/modeling.R',
  'https://github.com/jinworks/CellChat/blob/main/R/database.R',
  'https://bioconductor.org/packages/release/bioc/manuals/AUCell/man/AUCell.pdf',
  'https://davislaboratory.github.io/singscore/reference/rankGenes.html',
  'https://code.bioconductor.org/browse/UCell/',
  'https://yulab-smu.top/biomedical-knowledge-mining-book/universal-enrichment.html']}

cfg=load(DS/'pathways/scoring_config.json');assert cfg['source_sha256']==ready['RNA_sha256']
registry_path=DS/'pathways/RNA_scores_zero_mask_and_alias_v1/pathway_gene_overlap.tsv'
if not registry_path.is_file():registry_path=DS/'pathways/pathway_gene_overlap.tsv'
registry=pd.read_csv(registry_path,sep='\t')
selected=registry.loc[registry.included]
assert len(registry)==359 and len(selected)==356 and not registry.pathway_id.duplicated().any()
core=(DW/'code/pathway_core.R').read_text(); prior=(DW/'testing/code_at_equivalence/pathway_core.R').read_text()
math=lambda text: text.split('score_chunk <- function',1)[1].split('run_pathway_chunks <- function',1)[0].strip()
assert math(core)==math(prior)
depth=pd.read_csv(ROOT/'integration/cell_metadata.tsv.gz',sep='\t',usecols=['n_genes_symbol','MajorCellType'])
report['pathways']={'RNA_genes':33660,'cells':581725,'original_sets':359,'scored_sets':356,
 'excluded_sets':registry.loc[~registry.included].to_dict('records'),
 'small_sets_under_10':selected.loc[selected.RNA_gene_overlap<10,['pathway_id','pathway_name','RNA_gene_overlap']].to_dict('records'),
 'detected_RNA_quantiles':{str(q):float(v) for q,v in depth.n_genes_symbol.quantile([0,.1,.25,.5,.75,.9,1]).items()},
 'records_below_AUCell_top_1683':int((depth.n_genes_symbol<1683).sum()),
 'records_below_UCell_maxRank_1500':int((depth.n_genes_symbol<1500).sum()),
 'baseline_core_function_matches_original_official_equivalence':True,
 'baseline_equivalence_receipt_sha256':sha(DW/'testing/real_AD_1000x356_equivalence_receipt.json'),
 'formal_corrections':'AUCell explicit RNA-zero mask and verified gene-set-to-source aliases; all three final methods independently checked against a new official 1000x356 reference.',
 'algorithm_semantics':['All methods rank the same RNA universe; no ADT feature or DB-only renormalization.',
 'AUCell 1.24.0 seeded ranking with aucMaxRank=1683; explicit original-RNA>0 weight mask eliminates stochastic zero-gene contributions while preserving positive-gene ranks.',
 'UCell 2.6.2: maxRank=1500, average ties; earlier n*maxRank normalization, not the later >=2.7.6 formula.',
 'singscore 1.22.0: minimum ties, centerScore=FALSE. Scores of different methods have different scales.',
 'Individual-cell gene-set scores are exploratory relative-expression summaries, not disease association or activity calls.']}
crosswalk=load(OUT/'resource_gene_crosswalk/crosswalk_receipt.json')
assert crosswalk['accepted_renames']==1172 and crosswalk['KEGG_impact']['affected_scored_sets']==159
report['resource_matching']={'crosswalk_receipt':str(OUT/'resource_gene_crosswalk/crosswalk_receipt.json'),
 'crosswalk_receipt_sha256':sha(OUT/'resource_gene_crosswalk/crosswalk_receipt.json'),
 'approved_renames_sha256':crosswalk['files']['approved_symbol_renames.tsv']['sha256'],
 'RNA_source_axis_unchanged':True,'fixed_annotation_version':'org.Hs.eg.db 3.18.0',
 'accepted_RNA_renames':1172,'KEGG_aliases_in_all_359_source_sets':247,
 'KEGG_aliases_in_356_scored_sets':244,'KEGG_restored_scored_member_entries':799,
 'KEGG_sets_corrected':159,'KEGG_sets_unchanged':197,
 'matching':'Existing source RNA members first; otherwise canonical resource gene to unique original RNA source; deduplicate physical rows.',
 'strict_rejection':'Original source multiple Ensembl IDs, stable-ID multiple primary SYMBOL/NCBI IDs, LOC/readthrough, CCL3L copy ambiguity, canonical already present or nonunique RNA target.',
 'scope':'Versioned resource matching only; no claim of a fresh genome-wide HGNC annotation or source-RNA relabeling.',
 'independent_gate':str(WORK/'validation/HGNC_independent_review/global_resource_crosswalk_review.json')}

enrichment_receipt=load(ROOT/'enrichment/enrichment_receipt.json')
assert ready['RNA_sha256'] in set(strings(enrichment_receipt))
assert report['resource_matching']['approved_renames_sha256'] in set(strings(enrichment_receipt))
report['enrichment']={'receipt':str(ROOT/'enrichment/enrichment_receipt.json'),
 'receipt_sha256':sha(ROOT/'enrichment/enrichment_receipt.json'),
 'candidate_background_genes':29936,'rule':'RNA detected in >=3 QC-retained cells',
 'marker_definition':'Positive log2FC>1, within-coarse-group BH P<0.05, pct.1>=0.1; Unassigned excluded',
 'statistical_unit':'Balanced-cell coarse annotation marker discovery; exploratory marker-list ORA, not a patient-level or AD/HC test',
 'multiple_testing':'BH applied separately within each cell type and database; no joint family-wide FDR claim',
 'provisional_Fibroblasts':'Retained exploratory marker-enrichment row, single-source dominated; not validated AD-specific stroma'}
for method in ('GO','KEGG'):
    frame=pd.read_csv(ROOT/f'enrichment/{method}_all_tested.tsv',sep='\t')
    errors=[];probes=[]
    for group, part in frame.groupby('cell_type'):
        errors.append(float(np.max(np.abs(bh(part.pvalue.to_numpy())-part['p.adjust'].to_numpy()))))
        row=part.sort_values('p.adjust').iloc[0]
        k,n=map(int,row.GeneRatio.split('/'));K,N=map(int,row.BgRatio.split('/'))
        independent=float(hypergeom.sf(k-1,N,K,n))
        assert np.isclose(independent,row.pvalue,rtol=1e-10,atol=0)
        probes.append({'group':group,'term':row.ID,'reported_p':float(row.pvalue),'independent_hypergeom_p':independent})
    assert max(errors)<1e-10
    report['enrichment'][method]={'tested_rows':len(frame),
      'effective_annotated_background':sorted(set(map(int,frame.BgRatio.str.split('/').str[-1]))),
      'independent_BH_max_absolute_error':max(errors),'independent_hypergeometric_probes':probes}

ccfg=load(CS/'cellchat_config.json');assert ccfg['source_RNA_sha256']==ready['RNA_sha256']
assert ccfg['nboot']==100 and ccfg['seed']==2026 and ccfg['min_cells']==10
assert ccfg['coarse_review_sha256']==sha(ccfg['coarse_review_path'])
if CS.name in ('cellchat_exact_candidate','cellchat_sparse_exact_candidate') and selection is not None:
    impl_ref=selection.get('implementation_review')
    if impl_ref is None and args.cellchat_dir is not None:
        impl_path=OUT/('exact_sparse_mean_implementation_review.json' if CS.name=='cellchat_sparse_exact_candidate' else 'exact_mean_implementation_review.json')
        impl_ref={'status':'passed','path':str(impl_path),'sha256':sha(impl_path)}
    assert impl_ref['status']=='passed' and impl_ref['sha256']==sha(impl_ref['path'])
    implementation=load(impl_ref['path'])
    assert implementation['status']=='passed_implementation_review_not_formal_network'
    assert implementation['source_RNA_sha256']==ready['RNA_sha256']
    assert ccfg['workers']==1 and ccfg['nboot']==100 and ccfg['seed']==2026
    for name,value in implementation['reviewed_code_sha256'].items():assert ccfg['pipeline_code_sha256'][name]==value
    assert ccfg['DB_expression_cache_sha256']==implementation['cache_review']['cache_sha256']
    assert ccfg['DB_expression_cache_receipt_sha256']==implementation['cache_review']['receipt_sha256']
    if CS.name=='cellchat_sparse_exact_candidate':
        assert ccfg['implementation_review_sha256']==impl_ref['sha256']
        full_binding=implementation['full_axis_review']
        assert sha(full_binding['path'])==full_binding['sha256']==ccfg['full_axis_aggregation_review_sha256']
        full_proof=load(full_binding['path'])
        assert full_binding['status']==full_proof['status']=='validated_formal_axis_means_only_not_network'
        assert full_proof['source_RNA_sha256']==ready['RNA_sha256'] and full_proof['source_cache_sha256']==ccfg['DB_expression_cache_sha256']
        independent=implementation['independent_review']
        assert independent['status']==load(independent['path'])['status']=='passed_independent_sparse_helper_review_not_formal_network'
        assert sha(independent['path'])==independent['sha256']
        supplement=implementation['independent_full_axis_review']
        extra=load(supplement['path'])
        assert supplement['status']==extra['status']=='passed_independent_sparse_formal_axis_supplement_not_network'
        assert sha(supplement['path'])==supplement['sha256'] and extra['formal_axis_proof_sha256']==full_binding['sha256']
        assert ccfg['implementation_equivalence_sha256']==implementation['pilot_review']['sha256']
    report['implementation_review']={'path':impl_ref['path'],'sha256':impl_ref['sha256'],'runtime_code_cache_and_parameters_bound':True}
assert ccfg['symbol_alias_sha256']==sha(ccfg['symbol_alias_path'])
cc_aliases=pd.read_csv(ccfg['symbol_alias_path'],sep='\t')
global_aliases=pd.read_csv(OUT/'resource_gene_crosswalk/approved_symbol_renames.tsv',sep='\t')
alias_keys=['source_symbol','canonical_symbol','source_ensembl_id']
assert len(cc_aliases)==9
assert set(map(tuple,cc_aliases[alias_keys].to_numpy())).issubset(set(map(tuple,global_aliases[alias_keys].to_numpy())))
report['cellchat']={'min_cells':10,'nboot':100,'seed':2026,'p_value_cutoff':'strict p < 0.05',
 'verified_aliases':{'count':9,'table_sha256':ccfg['symbol_alias_sha256'],'same_pairs_and_stable_ids_as_global_resource_crosswalk':True},
 'parameters':{'triMean':True,'raw_use':True,'population_size':False,'distance_use':False,'PPI_projection':False},
 'input_normalization':'Full RNA log1p(cp10k), database-gene subset not renormalized',
 'expected_modeled_records':528739,'expected_modeled_groups':8,
 'excluded_records':{'Unassigned':52809,'Fibroblasts_provisional':172,'Mast_cells':4,'Neutrophils':1},
 'code_sha256':ccfg['pipeline_code_sha256'],
 'significance_scope':'Nominal cell-label-permutation P in 0.01 increments (100 permutations); not FDR-adjusted, no donor or condition replication.',
 'interaction_scope':'Pooled transcript-based hypotheses. Dissociated PBMC provide no measured physical contact or spatial validation.',
 'population_size_false_scope':'Per-group expression-based strength; no inference about total cohort abundance-weighted signaling.',
 'zero_edge_interpretation':'No prediction under these RNA/triMean/database/threshold settings, not evidence of biological absence.'}
inclusion=CS/'cell_inclusion.tsv.gz'
if inclusion.is_file():
    inc=pd.read_csv(inclusion,sep='\t');kept=inc.loc[inc.included]
    assert len(inc)==581725 and len(kept)==528739 and kept.MajorCellType.nunique()==8
    assert not kept.MajorCellType.isin(['Unassigned','Fibroblasts','Mast cells','Neutrophils']).any()
    report['cellchat']['actual_inclusion']={'records':len(kept),'groups':kept.MajorCellType.value_counts().to_dict(),
        'excluded_reasons':inc.loc[~inc.included].reason.value_counts().to_dict(),'inclusion_sha256':sha(inclusion)}
db_overlap=CS/'CellChatDB_gene_overlap.tsv'
if db_overlap.is_file():
    overlap=pd.read_csv(db_overlap,sep='\t')
    report['cellchat']['human_DB_genes']={'all':len(overlap),'matched_RNA_rows_direct_or_verified_alias':int(overlap.present_in_RNA.sum()),
        'missing':overlap.loc[~overlap.present_in_RNA,'DB_gene'].tolist()}
    internal_cache=DW/'performance_cache/normalized_DB_expression_cache_receipt.json'
    if internal_cache.is_file():
        cache_info=load(internal_cache)
        assert cache_info['source_RNA_sha256']==ready['RNA_sha256']
        if cache_info['status']=='validated_internal_cache_not_a_model_run':
            report['cellchat']['human_DB_genes']['matched_rows_with_zero_RNA_expression_in_all_retained_records']=cache_info['zero_nnz_genes_retained']
            report['cellchat']['human_DB_genes']['matching_is_not_expression_or_activity']=True

sensitivity=OUT/'aucell_zero_sensitivity/AUCell_zero_sensitivity.json'
if sensitivity.is_file(): report['AUCell_zero_sensitivity']={'receipt':str(sensitivity),'sha256':sha(sensitivity)}
masked_sensitivity=OUT/'aucell_zero_sensitivity/explicit_zero_mask_sensitivity.json'
if masked_sensitivity.is_file():report['AUCell_corrected_sensitivity']={'receipt':str(masked_sensitivity),'sha256':sha(masked_sensitivity)}

pending=[]
if selection is None:pending.append('explicit_cellchat_selection')
path_complete=DS/'pathways/pathway_completion.json'
if path_complete.is_file():
    completion=load(path_complete);assert completion['status']=='validated'
    assert completion['source_RNA_sha256']==ready['RNA_sha256'] and completion['pathway_count']==356
    assert completion['resource_gene_alias_table_sha256']==report['resource_matching']['approved_renames_sha256']
    assert completion['corrected_pathway_count']==159 and completion['unchanged_US_columns']==197
    assert completion['KEGG_aliases_in_all_source_sets']==247 and completion['KEGG_aliases_in_scored_sets']==244
    assert completion['KEGG_restored_member_entries']==799
    assert completion['AUCell_zero_mask']['status']=='validated'
    assert completion['cell_count']==581725 and len(completion['files'])==3
    assert {item['method'] for item in completion['files']}=={'AUCell','UCell','singscore'}
    assert completion['AUCell_zero_mask']['zero_detection_nonzero_scores']==0
    assert completion['AUCell_zero_mask']['unaffected_columns_bitwise_equal_baseline']
    probes=[];rows=np.unique(np.linspace(0,581724,41,dtype=int));cols=np.array([0,1,5,20,100,200,355])
    with h5py.File(completion['AUCell_zero_mask']['raw_scores'],'r') as raw,h5py.File(ready['RNA_path'],'r') as source, h5py.File(OUT/'resource_gene_crosswalk/official_1000x356_unified_crosswalk.h5','r') as official:
        original_ids=axis(source,'obs')[rows]
        for item in completion['files']:
            with h5py.File(item['path'],'r') as final:
                assert final['X'].shape==(581725,356)
                assert np.array_equal(axis(final,'obs')[rows],original_ids)
                a=final['X'][rows,:][:,cols];b=raw['scores/'+item['method']][rows,:][:,cols]
                assert np.array_equal(a,b) and np.isfinite(a).all()
                assert np.array_equal(final['obsm/X_umap'][rows,:],source['obsm/X_umap'][rows,:])
                assert np.array_equal(final['X'][:1000,:],np.asarray(official[item['method']][:],dtype=np.float32))
                assert np.array_equal(axis(final,'obs')[:1000],official['cell_ids'].asstr()[:])
                assert final['uns/scaid_scoring/resource_gene_alias_table_sha256'].asstr()[()]==report['resource_matching']['approved_renames_sha256']
                if item['method']=='AUCell':assert 'RNA>0' in final['uns/scaid_scoring/AUCell_zero_policy'].asstr()[()]
                assert np.array_equal(final['var/RNA_gene_overlap'][:],selected.RNA_gene_overlap.to_numpy())
                probes.append({'method':item['method'],'cells_checked':len(rows),'sets_checked':len(cols),
                    'score_raw_vs_H5AD_exact':True,'cell_ids_and_UMAP_preserved':True,
                    'full_first_1000x356_vs_official_reference_float32_equal':True,
                    'corrected_coverage_and_zero_policy_metadata_verified':True,
                    'production_reader':item['production_backed_reader']})
    report['pathways']['formal_completion']={'receipt_sha256':sha(path_complete),'independent_targeted_checks':probes}
    figure_receipt=DS/'pathways/pathway_figures_completion.json'
    if figure_receipt.is_file():
        figures=load(figure_receipt)
        assert figures['status']=='validated' and figures['figures']==1068 and figures['pathways_per_method']==356
        assert figures['source_score_receipt_sha256']==sha(path_complete)
        assert figures['all_cells_plotted'] and figures['no_cell_subsampling']
        report['pathways']['formal_rendering_receipt']={'sha256':sha(figure_receipt),**figures}
else: pending.append('pathway_completion.json')
cc_complete=CS/'cellchat_completion.json'
if cc_complete.is_file() and selection is not None:
    if 'cellchat_completion_sha256' in selection:assert selection['cellchat_completion_sha256']==sha(cc_complete)
    completion=load(cc_complete);assert completion['status']=='validated' and completion['source_RNA_sha256']==ready['RNA_sha256']
    computed=completion['computation'];assert computed['modeled_cells']==528739 and len(computed['modeled_cell_types'])==8
    assert computed['nboot']==100 and computed['seed']==2026 and computed['type']=='triMean'
    assert computed['raw_use'] and not computed['population_size'] and not computed['distance_use']
    assert computed['DB_genes_in_RNA']==1434 and computed['restored_RNA_symbol_rows']==9
    table=pd.read_csv(CS/'ligand_receptor_communications.tsv',sep='\t')
    assert (table.pval<.05).all() and (table.prob>0).all()
    counts=pd.read_csv(CS/'modeled_group_counts.tsv',sep='\t')
    assert counts.cell_count.sum()==528739 and len(counts)==8
    report['cellchat']['formal_completion']={'receipt_sha256':sha(cc_complete),
      'positive_LR_records_strict_p_below_05':len(table),'unique_LR':table.interaction_name.nunique(),
      'inferred_pathways':table.pathway_name.nunique(),'modeled_group_counts':counts.to_dict('records'),
      'rendering':completion['rendering']}
    net_review=OUT/f'independent_saved_cellchat_network_scope_{CS.name}.json'
    old=load(net_review) if net_review.is_file() else {}
    if old.get('cellchat_completion_sha256')!=sha(cc_complete):
        env=os.environ.copy();env.update({name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')})
        subprocess.run(['Rscript','--vanilla',str(WORK/'analysis/code/check_small_cellchat_network.R'),str(CS),str(net_review)],env=env,check=True,capture_output=True,text=True)
        record=load(net_review);record['cellchat_completion_sha256']=sha(cc_complete)
        net_review.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    assert load(net_review)['status']=='validated_saved_network_scope'
    report['cellchat']['independent_saved_network_review']={'path':str(net_review),'sha256':sha(net_review)}
else: pending.append('cellchat_completion.json')
report['pending']=pending
if not pending: report['status']='reviewed_with_interpretation_limits'
target=OUT/'downstream_scientific_review.json';tmp=target.with_suffix('.json.tmp')
tmp.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');tmp.replace(target)
note='''新 AD RNA 分析的通路与通讯结果独立复核

本次只读检查正式代码、输入标识、集合筛选、表达背景及少量结果，不重复全矩阵验证或运行昂贵算法。

评分使用同一33660基因RNA轴，原359集合保留356；三项超500基因上限被列入排除表。AUCell默认零值贡献已通过小样本实测，并采用显式RNA>0屏蔽、保留正表达基因排名和5%上限。资源匹配使用固定org.Hs.eg.db3.18及原始Ensembl核验的一对一键，1172个全RNA别名中247个涉及原359个KEGG集合，其中244个出现在实际评分的356个集合中，恢复799个集合成员条目，修复159列；197列UCell/singscore逐位保留。三法使用新映射的1000×356独立官方参照复验。小基因集合和低检测深度按探索性表达评分解释，分数不表示疾病激活。UCell注明2.6.2版本，三种方法分数不可直接互换。

GO/KEGG用实测29936基因候选背景，并对前景、背景和集合使用相同版本匹配键；实际有效背景随数据库注释交集确定，数值见JSON。每群超几何计算与独立计算一致，BH重算最大绝对误差<1e-10。富集解释细胞群标记，不是疾病差异；单库主导的172个间质样预测仍为待复核标签。

CellChat按RNA表达量、humanDB基因及复合物成员匹配；完整RNA归一化后取数据库基因，不再次归一化。预期528739记录/8群；排除52809个Unassigned、172个暂定Fibro、4个Mast和1个Neut。所有原RNA记录仍保留。通讯推断用triMean、100次细胞标签置换、名义P<0.05；该P没有多重检验校正，也不能视为患者级证据。没有空间信息，RNA预测不证明真实接触或AD特异通讯。

当前状态：STATUS
待核验正式结果：PENDING
完整证据与参数见同目录downstream_scientific_review.json。
'''
(OUT/'downstream_scientific_review_说明.txt').write_text(note.replace('STATUS',report['status']).replace('PENDING',', '.join(pending) or '已完成小量结果核验'))
if not pending:
    final={'status':'passed','signed_utc':datetime.now(timezone.utc).isoformat(),
      'source_RNA_sha256':ready['RNA_sha256'],
      'resource_gene_alias_table_sha256':report['resource_matching']['approved_renames_sha256'],
      'pathway_completion_sha256':sha(path_complete),'cellchat_completion_sha256':sha(cc_complete),
      'cellchat_dir':str(CS),'cellchat_selection':selection,
      'enrichment_receipt_sha256':sha(ROOT/'enrichment/enrichment_receipt.json'),
      'review_evidence':str(target),'review_evidence_sha256':sha(target),
      'review_script_sha256':sha(__file__),
      'network_review_script_sha256':sha(WORK/'analysis/code/check_small_cellchat_network.R'),
      'independent_network_review_sha256':sha(net_review),
      'scope':'Read-only scientific parameter, matching, statistical-scope and actual final-product checks. Small full-first-block official score reference and saved network arrays checked; no expensive full analyses repeated.',
      'interpretation_limits':['Mixed source AD/HC cohort with unresolved cell-level clinical/donor identities; no patient-level or disease/control inference',
      'RNA gene-set expression summaries and pooled predicted transcript communication; no proven pathway activation or physical contact',
      'Nominal cell-label permutation P<0.05 for CellChat, not multiple-testing-adjusted patient evidence',
      'Rare/source-dominated automatic labels remain provisional, with ambiguous/small/provisional groups excluded from communication modeling',
      'Reference-model annotation and RNA/ADT consistency review do not constitute team-expert confirmation or independent-cohort validation',
      'The 27 source libraries were processed as technical channels; pool/library labels do not establish independent donors or independently verified physical captures',
      'Genotype-based cell-to-donor demultiplexing was unavailable; RNA-based doublet prediction does not establish donor-singlet identity',
      'Retained library-barcode records are this RNA-only reanalysis, not a reconstruction of the publication demultiplexed participant-level analytical cohort',
      'Fixed org.Hs.eg.db3.18 stable-ID resource matching; not a live comprehensive HGNC reannotation']}
    final_path=OUT/'final_scientific_scope_review.json'
    final_tmp=final_path.with_suffix('.json.tmp')
    final_tmp.write_text(json.dumps(final,ensure_ascii=False,indent=2)+'\n')
    final_tmp.replace(final_path)
print(json.dumps({'status':report['status'],'pending':pending,'report':str(target)},ensure_ascii=False))

"""Unique-donor cell-type pseudobulk DE; independent primary and sensitivity runs.
This is an observational within-study SLE/reference analysis, not a cross-disease
mechanism or independent validation of a proposed novel discovery.
"""
from pathlib import Path
import json,sys
import numpy as np,pandas as pd
from formulaic import model_matrix
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from statsmodels.stats.multitest import multipletests
import pydeseq2
B=Path(__file__).resolve().parents[1];T=B/'tables'
IFN=['ISG15','IFI6','IFI44L','MX1','OAS1','IFIT1','IFIT3','RSAD2','ISG20','XAF1','LY6E','IFI44']
mode=sys.argv[1]if len(sys.argv)>1 else'primary_unique_batch'
results=[]
for ct in ['cM','T4','B']:
 prefix=f'GSE174188_{ct}_{mode}'
 print('Start',prefix,flush=True)
 count=pd.read_csv(T/(prefix+'_counts.tsv'),sep='\t',index_col=0);m=pd.read_csv(T/(prefix+'_metadata.tsv'),sep='\t',index_col=0)
 assert count.index.equals(m.index)and m.ind_cov.is_unique
 m['cohort']=m.cohort.astype(str);m['Age']=m.Age.astype(float);m['Age10']=(m.Age-40)/10
 assert m[['Age','Sex','pop_cov','condition']].notna().all().all()
 # Do not include individual fixed effects: disease is constant within individuals.
 design='~cohort + Age10 + Sex + pop_cov + condition'if m.cohort.nunique()>1 else'~Age10 + Sex + pop_cov + condition'
 dm=model_matrix(design,m);assert np.linalg.matrix_rank(dm.to_numpy())==dm.shape[1], 'Non-identifiable design'
 keep=(count>=10).sum(axis=0)>=10;count=count.loc[:,keep].astype(np.int64)
 dds=DeseqDataSet(counts=count,metadata=m,design=design,n_cpus=4,refit_cooks=True,quiet=True)
 dds.deseq2();stat=DeseqStats(dds,contrast=['condition','SLE','HD'],n_cpus=4,alpha=0.05,quiet=True);stat.summary()
 r=stat.results_df.rename_axis('gene').reset_index();q=dds.var[['_genewise_converged','_MAP_converged','_LFC_converged']].copy();q.to_csv(T/(prefix+'_fit_convergence.tsv'),sep='\t');bad=q.index[~q['_LFC_converged'].fillna(False)];r.loc[r.gene.isin(bad),['pvalue','padj']]=np.nan;eligible=r.pvalue.notna()&r.padj.notna();r.loc[eligible,'padj']=multipletests(r.loc[eligible,'pvalue'],method='fdr_bh')[1];r=r.sort_values('padj',na_position='last')
 r.to_csv(T/(prefix+'_DESeq2.tsv'),sep='\t',index=False);r[r.gene.isin(IFN+['HSPA1A','HSPA1B','IFI27'])].to_csv(T/(prefix+'_IFN_stress_watch.tsv'),sep='\t',index=False)
 norms=pd.DataFrame(dds.layers['normed_counts'],index=count.index,columns=count.columns);watch=[g for g in IFN if g in norms]
 norms[watch].to_csv(T/(prefix+'_IFN_normalized_pseudobulk.tsv'),sep='\t')
 s={'cell_type':ct,'mode':mode,'unique_donors':len(m),'n_HD':int(m.condition.eq('HD').sum()),'n_SLE':int(m.condition.eq('SLE').sum()),'target_source_cells':int(m.target_cells.sum()),'design':design,'design_full_rank':True,'design_columns':dm.columns.tolist(),'n_LFC_nonconverged_excluded':len(bad),'LFC_nonconverged_genes_excluded':bad.tolist(),'Age_covariate':'(Age-40)/10 for numerical scaling; no effect on disease coefficient','n_genes_prefilter':len(keep),'n_genes_retained':int(keep.sum()),'gene_filter':'at least 10 counts in at least 10 unique donors','n_genes_padj_nonmissing':int(r.padj.notna().sum()),'n_genes_FDR05':int(r.padj.lt(.05).sum()),'IFN12_watch':r[r.gene.isin(IFN)][['gene','log2FoldChange','pvalue','padj']].to_dict(orient='records'),'pydeseq2_version':pydeseq2.__version__,'multiple_testing':'Benjamini-Hochberg recomputed after exclusion of nonconvergent LFC fits with the original DESeq2 independent-filter/Cook mask, within each cell type and analysis; primary and cohort4 sensitivity reported separately, never counted as independent replication.','note':'Archived remapped source processed counts; no raw droplet QC implied. Age/Sex/self-reported pop_cov and processing cohort adjusted; residual within-cohort batch, clinical severity and treatment can remain.'}
 results.append(s);(B/(prefix+'_DE_summary.json')).write_text(json.dumps(s,indent=2)+'\n');print(json.dumps({k:s[k]for k in ['cell_type','mode','unique_donors','n_HD','n_SLE','n_genes_FDR05']}),flush=True)
(B/('DE_all_'+mode+'_summary.json')).write_text(json.dumps(results,indent=2)+'\n')

"""BH postprocessing after exclusion of nonconvergent LFC fits.
Keep the original DESeq2 independent-filter/Cook's mask fixed; recompute BH only
among finite, convergence-valid eligible p-values. No DESeq2 model is refitted.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
B=Path(__file__).resolve().parents[1];T=B/'tables';IFN=['ISG15','IFI6','IFI44L','MX1','OAS1','IFIT1','IFIT3','RSAD2','ISG20','XAF1','LY6E','IFI44'];checks=[]
for mode in ['primary_unique_batch','cohort4_sensitivity']:
 summaries=[]
 for ct in ['cM','T4','B']:
  pre=f'GSE174188_{ct}_{mode}';r=pd.read_csv(T/(pre+'_DESeq2.tsv'),sep='\t');q=pd.read_csv(T/(pre+'_fit_convergence.tsv'),sep='\t',index_col=0);bad=q.index[q['_LFC_converged'].eq(False)];old=r.padj.copy();oldn=int(old.lt(.05).sum());r.loc[r.gene.isin(bad),['pvalue','padj']]=np.nan
  eligible=r.pvalue.notna()&r.padj.notna();r.loc[eligible,'padj']=multipletests(r.loc[eligible,'pvalue'],method='fdr_bh')[1];assert not r.loc[r.gene.isin(bad),'pvalue'].notna().any();r=r.sort_values('padj',na_position='last');r.to_csv(T/(pre+'_DESeq2.tsv'),sep='\t',index=False);r[r.gene.isin(IFN+['HSPA1A','HSPA1B','IFI27'])].to_csv(T/(pre+'_IFN_stress_watch.tsv'),sep='\t',index=False)
  s=json.loads((B/(pre+'_DE_summary.json')).read_text());s['n_genes_FDR05']=int(r.padj.lt(.05).sum());s['n_genes_padj_nonmissing']=int(r.padj.notna().sum());s['IFN12_watch']=r[r.gene.isin(IFN)][['gene','log2FoldChange','pvalue','padj']].to_dict(orient='records');s['n_LFC_nonconverged_excluded']=len(bad);s['LFC_nonconverged_genes_excluded']=bad.tolist();s['multiple_testing']='Benjamini-Hochberg recomputed after exclusion of nonconvergent LFC fits, separately within each cell type and analysis, retaining the original DESeq2 independent-filter/Cook mask. Primary and cohort4 sensitivity overlap and are not independent replication.';(B/(pre+'_DE_summary.json')).write_text(json.dumps(s,indent=2)+'\n');summaries.append(s)
  checks.append({'model':pre,'eligible_valid_pvalues':int(eligible.sum()),'LFC_nonconvergent_genes':bad.tolist(),'FDR05_before_valid_only_BH':oldn,'FDR05_after_valid_only_BH':s['n_genes_FDR05'],'mask_scope':'Original independent-filter/Cooks padj-finite mask is frozen; invalid gene p/padj excluded before new BH; all remaining inferential values adjusted over this valid family.'})
 (B/('DE_all_'+mode+'_summary.json')).write_text(json.dumps(summaries,indent=2)+'\n')
(B/'valid_gene_BH_postprocessing_QA.json').write_text(json.dumps({'passed':True,'models':checks,'DESeq2_model_rerun':False,'reason':'Post-fit convergence failure must not remain in the BH correction family.'},indent=2)+'\n');print(json.dumps(checks,indent=2))

"""Review paired protein evidence without relabelling or removing source cells."""
from pathlib import Path
import argparse, hashlib, json, importlib.metadata
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def review(source, out):
    out.mkdir(parents=True,exist_ok=True);fig=out/'figures';fig.mkdir(exist_ok=True)
    before=sha(source);a=ad.read_h5ad(source)
    assert a.obs_names.is_unique and a.var_names.is_unique
    assert sp.issparse(a.X) and np.isfinite(a.X.data).all()
    assert (a.X.data>=0).all() and np.equal(a.X.data,np.floor(a.X.data)).all()
    totals=np.asarray(a.X.sum(axis=1)).ravel()
    assert np.array_equal(totals,a.obs.total_counts.to_numpy())
    names=np.array(a.uns['ADT_feature_names'],dtype=str)
    counts=a.obsm['ADT_counts'];counts=counts.toarray() if sp.issparse(counts) else np.asarray(counts)
    assert counts.shape==(a.n_obs,len(names)) and len(set(names))==len(names)
    assert np.isfinite(counts).all() and (counts>=0).all() and np.equal(counts,np.floor(counts)).all()
    # Seurat CLR, margin=2, applied to feature-by-cell input: each cell is
    # divided by exp(mean(log1p(ADT counts))) before log1p. Includes zero tags.
    clr=np.log1p(counts/np.exp(np.log1p(counts).mean(axis=1,keepdims=True)))
    panel={'CD19':'CD19','CD27':'CD27','CD24':'CD24','CD38':'CD38','CD21':'CR2',
           'IgM':'IGHM','IgD':'IGHD','CD138':'SDC1','CD5':'CD5','HLADR':'HLA-DRA','CD86':'CD86'}
    symbols=a.var.gene_symbol.astype(str).to_numpy()
    cells=a.obs.copy();cells.index.name='cell_id';passmask=cells.passes_QC.to_numpy(dtype=bool)
    correlations=[]
    for j,protein in enumerate(names):
        cells['ADT_raw_'+protein]=counts[:,j];cells['ADT_CLR_'+protein]=clr[:,j]
        gene=panel[protein];g=symbols==gene
        raw=np.asarray(a.X[:,g].sum(axis=1)).ravel()
        rna=np.log1p(np.divide(raw.astype(float)*10000,totals,out=np.zeros_like(raw,dtype=float),where=totals>0))
        cells['RNA_log1p10k_'+gene]=rna
        for channel,idx in cells.loc[passmask].groupby('capture_channel',observed=True).groups.items():
            pos=cells.index.get_indexer(idx)
            constant = np.ptp(rna[pos]) == 0 or np.ptp(clr[pos,j]) == 0
            rho = np.nan if constant else spearmanr(rna[pos],clr[pos,j]).statistic
            correlations.append({'capture_channel':channel,'protein':protein,'RNA_gene':gene,'cells':len(pos),
                                 'RNA_detection_fraction':float((raw[pos]>0).mean()),'rho_descriptive':float(rho),
                                 'correlation_status':'undefined_constant_measurement' if constant else 'defined',
                                 'use':'paired measurement descriptive association; no cell-level inferential P'})
    cells['review_RNA_B_CD19_lower_reference_decile']=False
    cells['review_Unknown_CD19_at_or_above_RNA_B_median']=False
    cut=[]
    for channel,idx in cells.groupby('capture_channel',observed=True).groups.items():
        sub=cells.loc[idx];b=sub.passes_QC & sub.annotation_coarse.eq('B cells')
        if int(b.sum())<20:continue
        q10,median=sub.loc[b,'ADT_CLR_CD19'].quantile([.1,.5]).to_numpy()
        bflag=sub.passes_QC & sub.annotation_coarse.eq('B cells') & sub.ADT_CLR_CD19.le(q10)
        uflag=sub.passes_QC & sub.annotation_coarse.eq('Unknown') & sub.ADT_CLR_CD19.ge(median)
        cells.loc[idx,'review_RNA_B_CD19_lower_reference_decile']=bflag
        cells.loc[idx,'review_Unknown_CD19_at_or_above_RNA_B_median']=uflag
        cut.append({'capture_channel':channel,'RNA_B_reference_cells':int(b.sum()),'CLR_CD19_q10':q10,'CLR_CD19_median':median,
                    'RNA_B_lower_decile_review_cells':int(bflag.sum()),'Unknown_high_CD19_review_cells':int(uflag.sum()),
                    'policy':'relative within-channel review indicators, not validated positivity gates or hard QC filters'})
    cells['review_below_source_protocol_500_genes']=cells.passes_QC & cells.n_genes_by_counts.lt(500)
    cells.to_csv(out/'GSE262240_paired_RNA_ADT_review_cells.tsv.gz',sep='\t',compression='gzip')
    pd.DataFrame(correlations).to_csv(out/'GSE262240_RNA_ADT_descriptive_correlations.tsv',sep='\t',index=False)
    pd.DataFrame(cut).to_csv(out/'GSE262240_ADT_reference_review_flags.tsv',sep='\t',index=False)
    summaries=[]
    for (channel,label),sub in cells.loc[passmask].groupby(['capture_channel','annotation_coarse'],observed=True):
        for p in names:
            q=sub['ADT_CLR_'+p].quantile([.25,.5,.75]).to_numpy()
            summaries.append({'capture_channel':channel,'RNA_coarse_label':label,'cells':len(sub),'protein':p,'CLR_q25':q[0],'CLR_median':q[1],'CLR_q75':q[2]})
    pd.DataFrame(summaries).to_csv(out/'GSE262240_protein_by_channel_RNA_label.tsv',sep='\t',index=False)
    sensitivity=cells.loc[passmask].groupby(['capture_channel','annotation_coarse'],observed=True).agg(
        original_QC_pass=('passes_QC','size'),below_source_500=('review_below_source_protocol_500_genes','sum'))
    sensitivity['remaining_if_source_500_applied']=sensitivity.original_QC_pass-sensitivity.below_source_500
    sensitivity['use']='review sensitivity only; original passes_QC unchanged'
    sensitivity.to_csv(out/'GSE262240_source_500_gene_sensitivity.tsv',sep='\t')
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    channels=list(cells.capture_channel.unique());f,axs=plt.subplots(4,2,figsize=(10,11),layout='constrained')
    for ax,channel in zip(axs.flat,channels):
        sub=cells.loc[passmask & cells.capture_channel.eq(channel).to_numpy()]
        groups=[sub.loc[sub.annotation_coarse.eq('B cells'),'ADT_CLR_CD19'].to_numpy(),
                sub.loc[sub.annotation_coarse.eq('Unknown'),'ADT_CLR_CD19'].to_numpy(),
                sub.loc[~sub.annotation_coarse.isin(['B cells','Unknown']),'ADT_CLR_CD19'].to_numpy()]
        ax.boxplot(groups,tick_labels=[f'RNA B\nn={len(groups[0])}',f'Unknown\nn={len(groups[1])}',f'Other RNA\nn={len(groups[2])}'],showfliers=False)
        ax.set_title(str(channel));ax.set_ylabel('CD19 protein (CLR)')
    f.suptitle('Paired protein review of RNA annotations\nDescriptive within-channel evidence; no relabelling or additional filtering',fontsize=12)
    f.savefig(fig/'GSE262240_CD19_protein_annotation_review.pdf');f.savefig(fig/'GSE262240_CD19_protein_annotation_review.png',dpi=180);plt.close(f)
    # Small complete numeric fixture enables independent comparison to Seurat.
    pd.DataFrame(counts[:120,:].T,index=names).to_csv(out/'CLR_Seurat_reference_input.tsv',sep='\t')
    pd.DataFrame(clr[:120,:].T,index=names).to_csv(out/'CLR_python_reference_output.tsv',sep='\t')
    assert sha(source)==before
    report={'input':str(source),'input_sha256':before,'script_sha256':sha(Path(__file__)),
            'all_source_cells':a.n_obs,'original_QC_pass_cells':int(passmask.sum()),'paired_proteins':len(names),
            'RNA_B_lower_CD19_review_cells':int(cells.review_RNA_B_CD19_lower_reference_decile.sum()),
            'Unknown_high_CD19_review_cells':int(cells.review_Unknown_CD19_at_or_above_RNA_B_median.sum()),
            'QC_pass_below_source_500_genes':int(cells.review_below_source_protocol_500_genes.sum()),
            'raw_RNA_and_ADT_counts_not_modified':True,'existing_QC_and_annotations_not_modified':True,
            'gene_duplicate_symbols_summed_for_each_marker':True,
            'clinical_limitation':'Eight source labels are not independently verified donors; common protocol versus ED8 source group discordance unresolved.',
            'interpretation':'Paired protein evidence is an independent measured modality, not blinded annotation accuracy. Review indicators are relative internal reference quantiles; no additional hard filters, disease-effect tests or automatic relabels.',
            'software':{p:importlib.metadata.version(p) for p in ['anndata','numpy','scipy','pandas','matplotlib']},
            'reference':'https://satijalab.org/seurat/articles/multimodal_vignette.html'}
    (out/'GSE262240_multimodal_review_manifest.json').write_text(json.dumps(report,ensure_ascii=False,indent=2));print(json.dumps(report,ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);s=p.parse_args();review(s.input,s.out)

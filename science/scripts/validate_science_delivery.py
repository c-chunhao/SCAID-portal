from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,pymupdf as fitz
from PIL import Image
from statsmodels.stats.multitest import multipletests
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures';checks={};fonts=[]
s=json.loads((B/'composition_recomputation_summary.json').read_text());assert s['release_cells']==2572480 and s['objects']==48 and s['specimen_labels']==625 and s['n_condition_labels']==24
u=pd.read_csv(T/'composition_specimen_denominators_and_eligibility.tsv',sep='\t');assert len(u)==625;assert not u.duplicated(['Data_object','Sample']).any()
for ln,expected in [('T_NK',937825),('B_Plasma',253635)]:
 for layer in ['layer1','layer2','layer3']:
  d=pd.read_csv(T/(ln+'_'+layer+'_complete_zero_specimen_proportions.tsv.gz'),sep='\t');classes=d.Cell_type.nunique();assert len(d)==625*classes;assert not d.duplicated(['Data_object','Sample','Cell_type']).any();assert int(d.Count.sum())==expected;assert d.Count.ge(0).all();assert d.Proportion.dropna().between(0,1).all()
  checks[ln+'_'+layer]={'specimens':d.Specimen_id.nunique(),'classes':classes,'rows':len(d),'zeros':int(d.Count.eq(0).sum()),'counts_sum':int(d.Count.sum()),'complete_cartesian_grid':True}
for p in T.glob('GSE174188_*_metadata.tsv'):
 m=pd.read_csv(p,sep='\t',index_col=0);assert m.ind_cov.is_unique;cc=pd.read_csv(str(p).replace('_metadata.tsv','_counts.tsv'),sep='\t',index_col=0);assert m.index.equals(cc.index);assert np.isfinite(cc.to_numpy()).all()and(cc.to_numpy()>=0).all();assert np.array_equal(cc.to_numpy(),np.rint(cc.to_numpy()));checks[p.stem]={'rows':len(m),'unique_donors':m.ind_cov.nunique(),'metadata_count_rows_exact_match':True,'count_nonnegative_integers':True}
for p in sorted(F.glob('*.pdf')):
 if p.name=='Figure3_marker_panels.pdf':continue
 doc=fitz.open(p)
 for i,pg in enumerate(doc):
  spans=[sp for bl in pg.get_text('dict')['blocks']for ln in bl.get('lines',[])for sp in ln['spans']if sp['text'].strip()]
  factor=min(1,468/pg.rect.width,576/pg.rect.height);sizes=[x['size']*factor for x in spans];minimum=min(sizes)if sizes else None
  fonts.append({'file':p.name,'page':i+1,'width_in':pg.rect.width/72,'height_in':pg.rect.height/72,'minimum_font_at_6_5in_by_max8in':minimum,'below_7pt':sum(x<6.999 for x in sizes)})
  assert minimum is None or minimum>=6.999,(p.name,i+1,minimum)
 doc.close()
rasters=[]
for p in sorted(F.glob('*.tiff')):
 with Image.open(p)as im:
  dpi=im.info.get('dpi');assert dpi and abs(dpi[0]-600)<.5
  rasters.append({'file':p.name,'pixels':list(im.size),'dpi':[float(x)for x in dpi],'one_frame':im.n_frames==1})
refs=json.loads((B/'annotation_reference_summary.json').read_text());assert refs['exact_source_barcodes']==664791 and refs['exact_source_donor_matches']==664791
checks['source_annotation_reference']=refs
refs_t1d=json.loads((B/'T1D_annotation_reference_summary.json').read_text());assert refs_t1d['exact_source_barcodes']==75746 and refs_t1d['exact_source_donor_matches']==75746 and refs_t1d['reference_comparable_lineage_cells']==74830 and refs_t1d['curated_other_unknown_among_comparable']==359
checks['T1D_source_annotation_reference']=refs_t1d
v=json.loads((B/'T1D_ontology_QA_correction.json').read_text());assert v['remaining_true_DC_discordance']['source_DC_cells']==1411 and v['remaining_true_DC_discordance']['classified_final_Mononuclear_Phagocytes']==324
checks['T1D_true_DC_disagreement']=v['remaining_true_DC_discordance']
for prefix,ncells,ngenes,nnz in [('BD',20375,16973,21859160),('pSS',29567,21176,42860863)]:
 a=json.loads((B/(prefix+'_full_retained_source_count_audit.json')).read_text());assert a['passed'] and a['retained_cells']==ncells and a['retained_genes']==ngenes and a['retained_sparse_nonzero_values']==nnz and a['remaining_numeric_differences']==0 and a['canonical_sparse_SHA256']==a['published_canonical_sparse_SHA256'] and not a['source_or_release_objects_mutated'] and a['count_values_changed']==0
 checks[prefix+'_full_retained_source_count_audit']=a
for mode in ['primary_unique_batch','cohort4_sensitivity']:
 for ct in ['cM','T4','B']:
  prefix=f'GSE174188_{ct}_{mode}';r=pd.read_csv(T/(prefix+'_DESeq2.tsv'),sep='\t');q=pd.read_csv(T/(prefix+'_fit_convergence.tsv'),sep='\t',index_col=0);bad=q.index[~q['_LFC_converged'].fillna(False)];assert r.loc[r.gene.isin(bad),['pvalue','padj']].isna().all().all();valid=r.padj.notna();assert np.allclose(r.loc[valid,'padj'],multipletests(r.loc[valid,'pvalue'],method='fdr_bh')[1],atol=1e-14,rtol=1e-12)
  checks['valid_gene_BH_'+prefix]={'finite_eligible_pvalues':int(valid.sum()),'nonconvergent_genes_excluded':bad.tolist(),'BH_recalculated_after_invalid_exclusion':True}
with pd.ExcelFile(B/'Supplementary_Table_S4_Composition_Annotation_Validation.xlsx') as x:
 assert set(['BD_full_counts_audit','pSS_full_counts_audit','T1D_author_confusion','T1D_DC_Mono_disagreements']).issubset(x.sheet_names)
 med=pd.read_excel(x,sheet_name='Condition_medians');assert 'Valid_fraction_specimens' in med
 checks['S4_final_scope']={'sheet_names':x.sheet_names,'valid_fraction_specimen_counts_reported':True,'BD_pSS_full_counts_and_T1D_disagreements_included':True}

out={'passed':True,'checks':checks,'figure_fonts':fonts,'tiffs':rasters,'scientific_scope':'Descriptive atlas panels; one-source-batch-per-biological-donor within-study comparisons; source-author reference concordance with disclosed non-blinded boundary.'}
(B/'science_validation.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',len(checks),'data checks;',len(fonts),'PDF page font checks;',len(rasters),'single-page TIFFs')

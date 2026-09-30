from pathlib import Path
import pandas as pd,json
B=Path(__file__).resolve().parents[1];T=B/'tables'
def read(name):return pd.read_csv(T/name,sep='\t')
with pd.ExcelWriter(B/'Supplementary_Table_S4_Composition_Annotation_Validation.xlsx',engine='openpyxl')as w:
 for sheet,fn in [('Specimen_denominators','composition_specimen_denominators_and_eligibility.tsv'),('Object_key','composition_object_key.tsv'),('Condition_medians','all_layers_condition_descriptive_medians.tsv'),('BD_full_counts_audit','BD_full_retained_count_audit_by_specimen.tsv'),('BD_feature_name_map','BD_date_corrupted_feature_symbol_reconciliation.tsv'),('Author_confusion','SLE_source_author_vs_curated_coarse_confusion_counts.tsv'),('Author_lineage_stats','SLE_source_author_concordance_per_lineage.tsv'),('Author_donor_stats','SLE_source_author_concordance_per_donor.tsv'),('One_batch_per_donor','GSE174188_frozen_one_batch_per_donor.tsv'),('Activity_switch_visits','GSE174188_activity_status_discordant_batches.tsv'),('Unique_donor_comp_IFN','GSE174188_unique_donor_HC3_composition_IFN_tests.tsv'),('Source_donor_composition','GSE174188_unique_donor_composition_complete_zeros.tsv')]:read(fn).to_excel(w,sheet_name=sheet,index=False)
 notes=[{'Scope':'All atlas composition figures','Interpretation':'Descriptive specimen distributions. Complete zeros, fixed immune denominator, source PsA/BD corrections. No pooled ANOVA retained.'},{'Scope':'Source-author annotation agreement','Interpretation':'GSE174188 and syn53641849, 2/48 objects, broad lineages. Original authors may have informed manual curation. Agreement is not blinded accuracy or whole-atlas validation; T1D author DC to final MonoPhag disagreement is retained.'},{'Scope':'SLE donor models','Interpretation':'One frozen source sequencing batch per verified ind_cov. Independent rows are donors, not cells/batches. Sensitivity analyses overlap primary and do not count as independent validation.'},{'Scope':'Remaining scientific contribution','Interpretation':'Known SLE IFN association reproduced; does not establish a novel cross-disease mechanism or an independently validated disease-specific discovery.'}];pd.DataFrame(notes).to_excel(w,sheet_name='Readme_scope',index=False)
 for sh in w.sheets.values():
  sh.freeze_panes='A2';sh.auto_filter.ref=sh.dimensions
  for col in sh.columns:sh.column_dimensions[col[0].column_letter].width=min(60,max(12,len(str(col[0].value))+3))
with pd.ExcelWriter(B/'Supplementary_Table_S5_Unique_Donor_DE.xlsx',engine='openpyxl')as w:
 for ct in ['cM','T4','B']:
  for mode,short in [('primary_unique_batch','primary'),('cohort4_sensitivity','cohort4')]:
   prefix=f'GSE174188_{ct}_{mode}';read(prefix+'_DESeq2.tsv').to_excel(w,sheet_name=ct+'_'+short,index=False)
   read(prefix+'_metadata.tsv').to_excel(w,sheet_name=ct+'_'+short+'_units',index=False)
   q=read(prefix+'_fit_convergence.tsv');q[q['_LFC_converged'].eq(False)].to_excel(w,sheet_name=ct+'_'+short+'_fit_excl',index=False)
 for sh in w.sheets.values():sh.freeze_panes='A2';sh.auto_filter.ref=sh.dimensions
print('S4/S5 workbooks written')

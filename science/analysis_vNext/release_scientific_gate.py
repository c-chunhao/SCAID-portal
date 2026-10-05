"""Fail closed before adding a dataset or making a donor disease claim.

This validates evidence completeness, not biological correctness. Existing
48-object exploratory release is distinct from admitting new source records.
"""
import argparse,json,re
from pathlib import Path

def evaluate(record):
    required=['source_accession','source_input_sha256','counts_schema_verified',
              'original_counts_provenance_verified','capture_mapping_verified',
              'clinical_group_mapping_verified','all_source_barcodes_preserved',
              'qc_flags_and_thresholds_recorded','actual_doublet_training_verified',
              'annotation_review_completed','source_duplicate_check_completed']
    boolean_required=[k for k in required if k not in ['source_accession','source_input_sha256']]
    missing=[k for k in boolean_required if record.get(k) is not True]
    if not isinstance(record.get('source_accession'),str) or not record['source_accession'].strip():missing.append('source_accession')
    if not isinstance(record.get('source_input_sha256'),str) or not re.fullmatch(r'[0-9a-fA-F]{64}',record.get('source_input_sha256','')):
        if 'source_input_sha256' not in missing:missing.append('source_input_sha256')
    release=not missing
    disease_requirements=['verified_cell_to_individual_mapping','appropriate_case_control_design',
                         'independent_individuals','full_rank_contrast_design',
                         'covariate_overlap_reviewed','raw_count_pseudobulk_and_valid_BH']
    disease_missing=[k for k in disease_requirements if record.get(k) is not True]
    if record.get('same_donor_overlap_with_validation',False) is not False:disease_missing.append('independent_validation')
    return {'accession':record.get('source_accession'),'new_release_eligible':release,
            'missing_release_evidence':missing,'donor_disease_contrast_eligible':release and not disease_missing,
            'missing_donor_contrast_evidence':disease_missing,
            'ambient_status':record.get('ambient_status','not_assessable_without_unfiltered_background'),
            'evidence_completeness_only':True,'imputed_clinical_or_donor_values_allowed':False}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('manifest',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=evaluate(json.loads(a.manifest.read_text()));a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))

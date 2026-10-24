#!/usr/bin/env python3
"""Versioned, strict RNA stable-ID nomenclature bridge; source data unmodified."""
from pathlib import Path
from collections import defaultdict
import pandas as pd, hashlib, json
from datetime import datetime, timezone

OUT=Path('@WORK@/analysis/validation/resource_gene_crosswalk')
AXIS=OUT/'source_RNA_feature_axis.tsv';ANNOT=OUT/'orgHs_3_18_source_ensembl_annotations.tsv'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
axis=pd.read_csv(AXIS,sep='\t',keep_default_na=False)
annotation=pd.read_csv(ANNOT,sep='\t',dtype=str,keep_default_na=False)
source_symbols=set(axis.source_symbol);annotation=annotation[annotation.SYMBOL.ne('NA')]
by_id={key:frame for key,frame in annotation.groupby('ENSEMBL',sort=False)}
canonical_ids=annotation.groupby('SYMBOL').ENSEMBL.unique().to_dict()
explicit_copy_uncertainty={'CCL3L1','CCL3L3'}
records=[];candidates=[]
for row in axis.itertuples(index=False):
    old=row.source_symbol;ids=row.source_feature_ids.split(';');record={'source_symbol':old,
      'source_ensembl_id':row.source_feature_ids,'resource_symbol':old,'action':'unchanged',
      'rename_accepted':False,'canonical_symbol':'','candidate_reason':'','annotation_version':'org.Hs.eg.db 3.18.0'}
    if len(ids)!=1:
        record['candidate_reason']='multiple_source_Ensembl_features';records.append(record);continue
    if ids[0] not in by_id:
        record['candidate_reason']='source_Ensembl_not_annotated_in_fixed_cache';records.append(record);continue
    rows=by_id[ids[0]];symbols=sorted(set(rows.SYMBOL)-{'NA',''})
    for canonical in symbols:
        part=rows.loc[rows.SYMBOL.eq(canonical)];aliases=set(part.ALIAS)-{'NA',''}
        ncbi_ids=set(part.ENTREZID)-{'NA',''};names=sorted(set(part.GENENAME)-{'NA',''})
        reasons=[]
        if len(symbols)!=1:reasons.append('source_Ensembl_has_multiple_official_symbols')
        if len(ncbi_ids)!=1:reasons.append('source_Ensembl_has_multiple_NCBI_gene_ids')
        if canonical.startswith('LOC'):reasons.append('LOC_unresolved_primary_symbol')
        if any('readthrough' in x.lower() for x in names):reasons.append('readthrough_identity_requires_review')
        if canonical in {'PEDS1-UBE2V1','CENPS-CORT','MICOS10-NBL1'}:reasons.append('compound_identity_requires_review')
        if canonical in explicit_copy_uncertainty or old in explicit_copy_uncertainty:reasons.append('ambiguous_CCL3L_copy_nomenclature')
        if canonical!=old:
            if canonical in source_symbols:reasons.append('canonical_already_present_on_RNA_axis')
            if old not in aliases:reasons.append('source_symbol_not_alias_of_same_stable_gene')
        candidate={'source_symbol':old,'canonical_symbol':canonical,'source_ensembl_id':ids[0],
          'source_alias_confirmed':old in aliases,'NCBI_gene_id':';'.join(sorted(ncbi_ids)),
          'official_gene_name':';'.join(names),'canonical_reference_Ensembl_ids':';'.join(sorted(canonical_ids[canonical])),
          'canonical_present_in_RNA':canonical in source_symbols,
          'identity_mapping':old==canonical,'accepted':not reasons,'rejection_reason':';'.join(reasons),
          'annotation_version':'org.Hs.eg.db 3.18.0'}
        candidates.append(candidate)
    record['candidate_reason']=';'.join(symbols);records.append(record)

c=pd.DataFrame(candidates)
# Enforce unique target from the whole RNA axis, before accepting any rename.
tentative=c.loc[c.accepted]
counts=tentative.groupby('canonical_symbol').source_symbol.nunique()
conflicts=set(counts[counts.gt(1)].index)
bad=c.accepted & c.canonical_symbol.isin(conflicts)
c.loc[bad,'accepted']=False
c.loc[bad,'rejection_reason']='multiple_RNA_source_symbols_map_to_same_canonical'
accepted=c[c.accepted].copy();assert not accepted.source_symbol.duplicated().any()
assert not accepted.canonical_symbol.duplicated().any()
renames=accepted[~accepted.identity_mapping].copy()
assert not renames.canonical_present_in_RNA.any() and renames.source_alias_confirmed.all()
assert renames.source_ensembl_id.str.match('^ENSG[0-9]+$').all()
assert not renames.NCBI_gene_id.str.contains(';',regex=False).any()
accepted_by_source=accepted.set_index('source_symbol')
full=pd.DataFrame(records)
for i,row in full.iterrows():
    old=row.source_symbol
    if old in accepted_by_source.index:
        x=accepted_by_source.loc[old];full.loc[i,'canonical_symbol']=x.canonical_symbol
        full.loc[i,'resource_symbol']=x.canonical_symbol
        full.loc[i,'action']='rename' if not x.identity_mapping else 'identity'
        full.loc[i,'rename_accepted']=not x.identity_mapping
        full.loc[i,'candidate_reason']='same_source_Ensembl_unique_official_gene;accepted'
    else:
        part=c.loc[c.source_symbol.eq(old)]
        if len(part):full.loc[i,'candidate_reason']='|'.join(sorted(set(part.rejection_reason)))
assert len(full)==33660 and not full.source_symbol.duplicated().any()
assert not full.resource_symbol.duplicated().any()
tables={'resource_symbol_crosswalk.tsv':full,'approved_symbol_renames.tsv':renames,
 'accepted_identity_mappings.tsv':accepted[accepted.identity_mapping],
 'all_mapping_candidates_with_reasons.tsv':c,
 'rejected_mapping_candidates.tsv':c[~c.accepted]}
for name,table in tables.items():table.to_csv(OUT/name,sep='\t',index=False)
cc=Path('@WORK@/analysis/validation/hgnc_alias_review/cellchat_symbol_aliases_verified.tsv')
known=pd.read_csv(cc,sep='\t');a=dict(zip(renames.source_symbol,renames.canonical_symbol))
assert all(a.get(r.source_symbol)==r.canonical_symbol for r in known.itertuples(index=False))
provenance=json.loads((OUT/'annotation_cache_provenance.json').read_text())
receipt={'status':'strict_versioned_resource_crosswalk_ready','created_utc':datetime.now(timezone.utc).isoformat(),
 'source_RNA_sha256':'b9aa963123cb8f010d3e092e03507b70d02170afd93d32a6cdea2efe2a2c09ec',
 'source_RNA_axis_unchanged':True,'source_counts_and_normalization_unchanged':True,
 'annotation':'org.Hs.eg.db 3.18.0 fixed official NCBI Gene-derived cache; not a live HGNC nomenclature refresh',
 'strict_rules':['Exactly one original Ensembl feature per RNA source symbol',
 'That source Ensembl resolves to exactly one official SYMBOL and one NCBI Gene ID in the fixed cache',
 'For rename: source old symbol is an ALIAS of that same stable gene, canonical is absent from original RNA axis',
 'Canonical target has exactly one accepted original RNA source; no source-axis merging or expression summing',
 'LOC, readthrough, compound identity and ambiguous CCL3L copy nomenclature are excluded',
 'Other reference Ensembl IDs are retained as evidence; only the exact unique present source Ensembl can support a mapping'],
 'source_RNA_features':len(full),'accepted_renames':len(renames),'accepted_identity_rows':int(accepted.identity_mapping.sum()),
 'all_candidate_rows':len(c),'rejected_candidate_rows':int((~c.accepted).sum()),
 'unchanged_unresolved_rows':int(full.action.eq('unchanged').sum()),
 'rejection_reason_counts':c.loc[~c.accepted,'rejection_reason'].value_counts().to_dict(),
 'existing_9_CellChat_table_preserved_sha256':sha(cc),'all_9_CellChat_mappings_pass_this_gate':True,
 'files':{name:{'path':str(OUT/name),'rows':len(table),'sha256':sha(OUT/name)} for name,table in tables.items()},
 'source_axis_dump_sha256':sha(AXIS),'cache_annotation_dump_sha256':sha(ANNOT),
 'annotation_database_sha256':sha(provenance['database'])}
(OUT/'crosswalk_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('status','accepted_renames','accepted_identity_rows','rejected_candidate_rows','unchanged_unresolved_rows')},ensure_ascii=False))

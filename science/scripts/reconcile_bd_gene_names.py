import os
"""Bounded reconciliation of the five archived BD diagnostic sparse vectors.
The source features have date-converted gene symbols. Stable Ensembl IDs identify
correct names. No source or published count value is changed.
"""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,h5py
B=Path(__file__).resolve().parents[1];T=B/'tables';OLD=Path(os.environ['SCAID_HISTORICAL_QC_AUDIT_DIR'])
fp=Path(os.environ['SCAID_BD_SOURCE_FEATURE_FILE']);f=pd.read_csv(fp,sep='\t',header=None,names=['ensembl_id','source_symbol','feature_type'])
# Read10X uses make.unique on symbols in original feature row order.
seen={};unique=[]
for gene in f.source_symbol:
 n=seen.get(gene,0);seen[gene]=n+1;unique.append(gene if n==0 else gene+'.'+str(n))
f['replayed_unique_symbol']=unique
api=json.loads((T/'BD_date_gene_ensembl_lookup.json').read_text());x=f[f.ensembl_id.isin(api)].copy();x['ensembl_current_symbol']=x.ensembl_id.map(lambda k:(api[k]or{}).get('display_name'));assert x.ensembl_current_symbol.notna().all();x['evidence_url']='https://www.ensembl.org/Homo_sapiens/Gene/Summary?g='+x.ensembl_id
x.to_csv(T/'BD_date_corrupted_feature_symbol_reconciliation.tsv',sep='\t',index=False)
correction=dict(zip(x.replayed_unique_symbol,x.ensembl_current_symbol));e=pd.read_csv(OLD/'SV_GSE198616_PBMC_RNA.historical_source_counts_first5.tsv',sep='\t');e['canonical_gene']=e.gene.replace(correction)
assert not e.duplicated(['cell','canonical_gene']).any(), 'No gene collision silently collapsed'
with h5py.File(os.environ['SCAID_BD_RELEASE_H5'])as h:
 ids=h['names_obs'].asstr()[:];genes=h['names_var'].asstr()[:];blocks=h['assay/RNA/layers/rawdata'];parts=[];offset=0
 for key in sorted(blocks):
  block=blocks[key];ptr=block['indptr'][:];N=len(ptr)-1
  for cell in e.cell.unique():
   ix=np.where(ids==cell)[0];assert len(ix)==1;row=int(ix[0])-offset
   if row<0 or row>=N:continue
   a,z=int(ptr[row]),int(ptr[row+1]);gi=block['indices'][a:z];v=block['data'][a:z];parts.append(pd.DataFrame({'cell':cell,'canonical_gene':genes[gi],'published':v}))
  offset+=N
p=pd.concat(parts,ignore_index=True);d=e.rename(columns={'value':'replayed'}).merge(p,on=['cell','canonical_gene'],how='outer');d[['replayed','published']]=d[['replayed','published']].fillna(0);d['difference']=d.published-d.replayed
bad=d[d.difference.ne(0)];assert len(bad)==0
changed=e[e.gene.ne(e.canonical_gene)];changed.to_csv(T/'BD_first5_expressed_gene_alias_changes.tsv',sep='\t',index=False)
s={'source_object':'SV_GSE198616_PBMC_RNA','source_feature_file':str(fp),'source_feature_SHA256':hashlib.sha256(fp.read_bytes()).hexdigest(),'reference_gene_ids':'Stable Ensembl IDs supplied in original source features.tsv.gz. Official Ensembl REST lookup/id, retrieved 2026-09-30, archived response.','cause':'Source feature symbols had been converted to spreadsheet-style dates (e.g. 1-Mar, 6-Sep). Read10X make.unique added .1 for duplicate dates representing different genuine genes. Published names were corrected; literal name comparison therefore failed although count values matched.','source_feature_date_rows_reconciled':len(x),'sampled_cells':e.cell.nunique(),'sampled_sparse_nonzero_entries':len(e),'expressed_date_alias_entries_changed_for_comparison':len(changed),'sampled_counts_after_stable_ID_gene_name_reconciliation_exact_match':True,'count_values_changed':0,'gene_name_collision_after_mapping':False,'full_remaining_difference_entries':len(bad),'published_or_source_objects_mutated':False,'limits':'Five complete sampled sparse vectors only, not an audit of every BD cell or all 48 source objects. Numeric filtering and all BD barcode/library-size/basic metric checks were already exact. Does not recover original execution logs or original doublet/ambient software versions.'}
(B/'BD_first5_source_count_name_reconciliation.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))

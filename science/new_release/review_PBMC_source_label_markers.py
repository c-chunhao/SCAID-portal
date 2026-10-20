#!/usr/bin/env python3
"""Source-specific descriptive marker evidence, preserving every original label."""
import json,os,time
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem,sparse_dataset
from prepare_new_aid_integration import WORK,PLAN,STUDY_SCOPE,sha,lineage

GENES=['CD3D','CD3E','TRAC','CD8A','MS4A1','CD79A','CD79B','CD37','CD74',
       'NKG7','GNLY','PRF1','JCHAIN','MZB1','LYZ','LST1','S100A8','S100A9','CTSS','FCN1','CD14','FCGR3A']
OUT=WORK/'integration/source_label_marker_review'
CAP=300

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261020);rows=[];inputs=[]
    for obj in PLAN['objects']:
        if obj['catalog']['tissue']!='PBMC':continue
        p=WORK/'matrices/objects'/obj['object_id'].replace('/','_')/'export_report.json'
        report=json.loads(p.read_text())
        if report['status']!='validated':raise ValueError('Source input is not validated')
        item=report['matrices']['RNA'];source=obj['catalog']['dataset_id']
        context=STUDY_SCOPE.get(source,{})
        print('MARKER_REVIEW',obj['object_id'],flush=True)
        with h5py.File(item['path'],'r',rdcc_nbytes=64*1024*1024,rdcc_nslots=10007) as h:
            obs=read_elem(h['obs']);var=read_elem(h['var']);labels=obs['MajorCellType'].astype('string').fillna('[missing]')
            valid=[g for g in GENES if g in var.index];gi=np.asarray([var.index.get_loc(g) for g in valid])
            matrix=sparse_dataset(h['X'])
            for label in sorted(labels.unique()):
                eligible=np.flatnonzero(labels.eq(label).to_numpy())
                idx=np.sort(rng.choice(eligible,size=min(CAP,len(eligible)),replace=False))
                values=matrix[idx,:][:,gi]
                means=np.asarray(values.mean(0)).ravel();fractions=np.asarray((values>0).mean(0)).ravel()
                lookup={g:k for k,g in enumerate(valid)}
                coarse=lineage(label)
                eligible_label=(coarse is not None and not context.get('exclude_from_lymphoid_pilot',False)
                    and (not context.get('B_lineage_only_in_pilot',False) or coarse in ['B cells','Plasma cells']))
                for g in GENES:
                    k=lookup.get(g)
                    rows.append({'object_id':obj['object_id'],'source':source,'condition':obj['catalog']['disease_name'],
                        'original_MajorCellType':str(label),'total_source_label_cells':len(eligible),
                        'review_cells':len(idx),'maximum_review_cells_per_object_label':CAP,
                        'assay_scope':context.get('assay_scope','Source-specific PBMC assay'),
                        'provided_label_is_eligible_for_final_pilot':eligible_label,
                        'gene':g,'feature_available':k is not None,
                        'mean_source_RNA_data':float(means[k]) if k is not None else None,
                        'fraction_source_RNA_data_positive':float(fractions[k]) if k is not None else None})
        inputs.append({'object_id':obj['object_id'],'path':item['path'],'sha256':item['sha256']})
    frame=pd.DataFrame(rows);frame.to_csv(OUT/'source_by_provided_label_marker_summary.tsv',sep='\t',index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    short={'Active_Tattoo_Uveitis':'ATU','Inactive_Tattoo_Uveitis':'ITU','Antiphospholipid_Syndrome':'APS',
           'Sarcoidosis':'Sarcoidosis','Crohn_Disease':'CD','Ulcerative_Colitis':'UC','T1D':'T1D'}
    frame['display_group']=[f'{source} · {short.get(obj.split("/")[0],obj.split("/")[0])} · {label}'
          for source,obj,label in zip(frame.source,frame.object_id,frame.original_MajorCellType)]
    pivot=frame.pivot(index='display_group',columns='gene',values='mean_source_RNA_data').reindex(columns=GENES)
    fig,ax=plt.subplots(figsize=(17,max(6,len(pivot)*.27)),layout='constrained')
    cmap=plt.get_cmap('YlGnBu').copy();cmap.set_bad('#e4e4e4')
    im=ax.imshow(pivot.to_numpy(),aspect='auto',cmap=cmap)
    ax.set_xticks(range(len(GENES)),GENES,rotation=55,ha='right',fontsize=9)
    ax.set_yticks(range(len(pivot)),pivot.index,fontsize=8)
    ax.set_title('Marker expression by source and provided cell label',loc='left',fontsize=12)
    fig.colorbar(im,ax=ax,label='Mean original RNA data-layer expression',shrink=.5)
    fig.savefig(OUT/'source_by_provided_label_marker_heatmap.png',dpi=220)
    fig.savefig(OUT/'source_by_provided_label_marker_heatmap.pdf');plt.close(fig)
    result={'status':'completed_descriptive_source_label_review','source_objects':len(inputs),
         'rows':len(frame),'marker_genes':GENES,'maximum_review_cells_per_object_label':CAP,
         'data_layer':'Original RNA X/data from each validated source H5AD; not the pilot embedding',
         'selection':'Deterministic random sample within each object and original MajorCellType label, capped at 300 cells.',
         'interpretation':'Marker summaries support manual review of provided labels and assay carryover. No labels are changed; absence of a source gene is represented as missing, not zero.',
         'assay_warning':'GSE268626 selected monocytes and is excluded from the final lymphoid pilot. Its provided T/B labels may reflect carryover or require annotation review; enrichment alone does not prove a label wrong.',
         'input_files':inputs,'table':str(OUT/'source_by_provided_label_marker_summary.tsv'),
         'table_sha256':sha(OUT/'source_by_provided_label_marker_summary.tsv'),
         'completed':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (OUT/'source_label_marker_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps({k:result[k] for k in ['status','source_objects','rows','table']}),flush=True)

if __name__=='__main__':main()

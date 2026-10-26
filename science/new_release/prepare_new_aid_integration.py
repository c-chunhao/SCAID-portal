#!/usr/bin/env python3
"""Inventory all new objects and produce a descriptive, bounded PBMC pilot."""
import argparse, collections, fcntl, hashlib, json, os, time
from pathlib import Path
from importlib.metadata import version
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('NUMBA_NUM_THREADS','4')
import anndata as ad
from anndata.io import read_elem, sparse_dataset
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
OUT=WORK/'integration'
PLAN=load_plan()
SEED=20261020
SCRIPT_SIGNATURE=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
STUDY_SCOPE={
    'GSE268626':{'assay_scope':'Peripheral blood monocytes isolated by negative magnetic selection; enriched monocyte assay',
                 'clinical_scope':'Newly diagnosed treatment-naive Crohn disease/ulcerative colitis and source-defined healthy references',
                 'exclude_from_lymphoid_pilot':True,
                 'B_lineage_only_in_pilot':False,
                 'verification':'Official GEO SOFT independently checked; independent official GEO assay-scope review (21 October 2026)'},
    'GSE262240':{'assay_scope':'Magnetically sorted B cells from PBMC; enriched B-cell assay',
                 'clinical_scope':'Triple-positive primary antiphospholipid syndrome and source-defined reference samples',
                 'B_lineage_only_in_pilot':True,
                 'verification':'Official GEO SOFT independently checked; official GEO SOFT for GSE262240'},
    'GSE270142':{'assay_scope':'PBMC-derived PE-binding/islet-antigen-reactive B cells; enriched B-cell assay',
                 'clinical_scope':'ND_FDR, autoantibody-positive at-risk and recent-onset T1D source groups',
                 'B_lineage_only_in_pilot':True,
                 'verification':'Official GEO SOFT independently checked; official GEO SOFT for GSE270142'},
    'GSE271063':{'assay_scope':'Study-defined PBMC assay',
                 'clinical_scope':'TN10 at-risk participants; Teplizumab/Placebo are treatment arms, not disease/healthy-control labels',
                 'B_lineage_only_in_pilot':False,
                 'verification':'Official GEO SOFT independently checked; official GEO SOFT for GSE271063'}
}

def sha(path):
    x=hashlib.sha256()
    with open(path,'rb') as f:
        while b:=f.read(8*1024*1024):x.update(b)
    return x.hexdigest()

def reports():
    result=[]
    for obj in PLAN['objects']:
        p=WORK/'matrices/objects'/obj['object_id'].replace('/','_')/'export_report.json'
        r=json.loads(p.read_text()) if p.exists() else {'status':'pending'}
        result.append((obj,r))
    return result

def inventory():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[];groups=[];labels=[];overlap=[];source_seen={};inputs=[]
    for obj,r in reports():
        cat=obj['catalog'];row={'object_id':obj['object_id'],'condition':cat['disease_name'],
              'source':cat['dataset_id'],'tissue':cat['tissue'],'cell_records':cat['cell_count'],
              'status':r['status'],'rna_counts_verified':False,'sample_label_count':None,
              'sample_labels_are_authenticated_donors':False}
        context=STUDY_SCOPE.get(cat['dataset_id'],{})
        row.update(assay_scope=context.get('assay_scope','Retain source assay context; no unselected-whole-tissue composition assumption'),
                   clinical_scope=context.get('clinical_scope','Original clinical/group labels retained; donor-level design requires confirmation'))
        if r['status']=='validated':
            rr=r['matrices']['RNA'];has_counts=rr.get('counts_available','counts' in rr.get('matrix_checks',{}))
            row['rna_counts_verified']=bool(rr['validated'] and has_counts)
            row['raw_count_status']=rr.get('raw_count_status','integer_count_layer_preserved_from_processed_source' if has_counts else 'unavailable')
            row['gene_count']=rr['feature_count'];row['rna_path']=rr['path'];row['rna_sha256']=rr['sha256']
            with h5py.File(rr['path'],'r') as h:meta=read_elem(h['obs'])
            sample=meta['Sample'] if 'Sample' in meta else pd.Series(index=meta.index,dtype='object')
            sample_text=sample.astype('string').str.strip()
            valid_sample=sample_text.notna()&~sample_text.str.lower().isin(['','na','nan','null','none','unknown'])
            row['sample_label_count']=int(sample_text[valid_sample].nunique())
            row['sample_label_missing_cells']=int((~valid_sample).sum())
            for name in ('Group','GSM','MajorCellType','nCount_RNA','nFeature_RNA','percent.mt','doublet_status','doublet_score','ambient_fraction'):
                row[name+'_present']=name in meta
                row[name+'_missing_cells']=int(meta[name].isna().sum()) if name in meta else len(meta)
            for col in ('Group','Disease_subtype','Disease_detail','MajorCellType','doublet_status'):
                if col in meta:
                    for label,n in meta[col].astype('string').fillna('[missing]').value_counts().items():
                        labels.append({'object_id':obj['object_id'],'metadata_column':col,'original_label':label,'cell_count':int(n)})
            if 'Group' in meta:
                for label,n in meta['Group'].astype('string').fillna('[missing]').value_counts().items():
                    groups.append({'source':cat['dataset_id'],'tissue':cat['tissue'],'object_id':obj['object_id'],
                                   'condition':cat['disease_name'],'original_Group_label':label,'cell_records':int(n)})
            # Deduplicate only exact cell identifiers within the same accession and tissue.
            key=(cat['dataset_id'],cat['tissue'])
            known=source_seen.setdefault(key,{})
            for prev,ids in known.items():
                n=len(set(meta.index).intersection(ids))
                if n:overlap.append({'source':key[0],'tissue':key[1],'object_1':prev,'object_2':obj['object_id'],'shared_source_cell_ids':n})
            known[obj['object_id']]=set(meta.index)
            inputs.append({'object_id':obj['object_id'],'condition':cat['disease_name'],'source':cat['dataset_id'],
                           'tissue':cat['tissue'],'RNA':rr['path'],'sha256':rr['sha256'],
                           'counts_layer':'layers/counts' if has_counts else None,
                           'source_fractional_processed_layer':None if has_counts else 'layers/source_processed_counts',
                           'counts_available':bool(has_counts),'raw_count_status':row['raw_count_status'],
                           'normalized_layer':'X','source_annotation':'obs/MajorCellType',
                           'sample_label':'obs/Sample','donor_id_authenticated':False})
        rows.append(row)
    frame=pd.DataFrame(rows);frame.to_csv(OUT/'integration_input_inventory.tsv',sep='\t',index=False)
    condition_design=[]
    for condition,subset in frame.groupby('condition',sort=True):
        condition_design.append({'condition':condition,'object_count':len(subset),
             'source_count':int(subset.source.nunique()),'tissue_count':int(subset.tissue.nunique()),
             'cell_records':int(subset.cell_records.sum()),'multiple_sources_present':bool(subset.source.nunique()>1),
             'sources':';'.join(sorted(subset.source.unique())),'tissues':';'.join(sorted(subset.tissue.unique()))})
    pd.DataFrame(condition_design).to_csv(OUT/'condition_source_coverage.tsv',sep='\t',index=False)
    matched=[]
    for (source,tissue),subset in frame.groupby(['source','tissue'],sort=True):
        if subset.condition.nunique()>1:
            matched.append({'source':source,'tissue':tissue,'conditions':';'.join(sorted(subset.condition.unique())),
                  'object_count':len(subset),'cell_records':int(subset.cell_records.sum()),
                  'sample_labels_per_object':';'.join(f'{r.object_id}={r.sample_label_count}' for r in subset.itertuples()),
                  'design_status':'Same source and tissue; donor identity, clinical stage, treatment and repeated sampling require confirmation before comparisons.'})
    pd.DataFrame(matched).to_csv(OUT/'within_source_tissue_comparison_candidates.tsv',sep='\t',index=False)
    pd.DataFrame(groups).to_csv(OUT/'source_tissue_group_labels.tsv',sep='\t',index=False)
    pd.DataFrame(labels).to_csv(OUT/'original_annotation_and_qc_labels.tsv',sep='\t',index=False)
    pd.DataFrame(overlap,columns=['source','tissue','object_1','object_2','shared_source_cell_ids']).to_csv(OUT/'same_source_cell_overlap.tsv',sep='\t',index=False)
    pd.crosstab(frame['source'],frame['condition'],values=frame['cell_records'],aggfunc='sum').fillna(0).astype('int64').to_csv(OUT/'source_by_condition_cell_records.tsv',sep='\t')
    pd.crosstab(frame['tissue'],frame['condition'],values=frame['cell_records'],aggfunc='sum').fillna(0).astype('int64').to_csv(OUT/'tissue_by_condition_cell_records.tsv',sep='\t')
    result={'created':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'object_count':len(rows),
            'validated_objects':sum(r['status']=='validated' for _,r in reports()),'cell_records':int(frame.cell_records.sum()),
            'input_objects':inputs,'shared_cell_id_pairs':overlap,
            'interpretation':'Cell records, original sample labels and original Group labels are descriptive metadata, not verified independent donors or disease-control assignments.',
            'integration_policy':{'separate_tissues':True,'preserve_original_counts':True,'do_not_overwrite_published_atlases':True,
                                  'no_case_control_pvalues_without_authenticated_donor_and_design':True}}
    (OUT/'integration_input_manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    return result

def lineage(x):
    s=str(x).strip().lower()
    if s in {'t cell','t cells','t_cell','t_cells'}:return 'T cells'
    if s in {'nk cell','nk cells','nk_cell','nk_cells'}:return 'NK cells'
    if s in {'b cell','b cells','b_cell','b_cells'}:return 'B cells'
    if s in {'plasma cell','plasma cells','plasma_cell','plasma_cells','plasmablast','plasmablasts'}:return 'Plasma cells'
    return None

def neighbor_metrics(embedding,obs):
    from sklearn.neighbors import NearestNeighbors
    ix=NearestNeighbors(n_neighbors=31,n_jobs=2).fit(embedding).kneighbors(embedding,return_distance=False)[:,1:]
    study=obs['scaid_source'].astype(str).to_numpy();celltype=obs['pilot_coarse_cell_type'].astype(str).to_numpy()
    sources=np.unique(study);p=np.stack([(study[ix]==s).mean(1) for s in sources],axis=1)
    ent=-(np.where(p>0,p*np.log(np.maximum(p,1e-15)),0)).sum(1)/max(np.log(len(sources)),1)
    purity=(celltype[ix]==celltype[:,None]).mean(1)
    return {'source_neighbor_entropy_mean':float(ent.mean()),'source_neighbor_entropy_median':float(np.median(ent)),
            'provided_coarse_label_neighbor_agreement_mean':float(purity.mean()),'neighbors':30,
            'interpretation':'Technical diagnostics only. Source mixing does not establish biological validity; provided labels are not an independent reference.'}

def pilot():
    import scanpy as sc
    import harmonypy as hm
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    dest=OUT/'PBMC_pilot';dest.mkdir(parents=True,exist_ok=True)
    selected=[(o,r) for o,r in reports() if o['catalog']['tissue']=='PBMC'
              and not STUDY_SCOPE.get(o['catalog']['dataset_id'],{}).get('exclude_from_lymphoid_pilot',False)]
    missing=[o['object_id'] for o,r in selected if r['status']!='validated']
    if missing:raise RuntimeError('Pilot waits for validated RNA inputs: '+str(missing))
    existing_path=dest/'pilot_result.json'
    if existing_path.exists():
        old=json.loads(existing_path.read_text())
        old_sources={o['object_id']:o['source_sha256'] for o in old.get('provenance',{}).get('source_objects',[])}
        current_sources={o['object_id']:r['matrices']['RNA']['sha256'] for o,r in selected}
        artifact=Path(old.get('output','/nonexistent'))
        if old.get('script_sha256')==SCRIPT_SIGNATURE and old_sources==current_sources and artifact.is_file() and sha(artifact)==old.get('sha256'):
            print('REUSE_VERIFIED_PILOT',artifact,flush=True)
            return old
    rng=np.random.default_rng(SEED);datasets=[];selection=[];seen=set()
    for obj,r in selected:
        path=r['matrices']['RNA']['path']
        with h5py.File(path,'r',rdcc_nbytes=64*1024*1024,rdcc_nslots=10007) as h:
            obs=read_elem(h['obs']);var=read_elem(h['var'])
            coarse=obs['MajorCellType'].map(lineage).astype('object')
            context=STUDY_SCOPE.get(obj['catalog']['dataset_id'],{})
            if context.get('B_lineage_only_in_pilot'):
                coarse=coarse.where(coarse.isin(['B cells','Plasma cells']),None)
            eligible=np.flatnonzero(coarse.notna().to_numpy())
            nonduplicate=[];duplicate_count=0
            for k in eligible:
                source_key=(obj['catalog']['dataset_id'],str(obs.index[k]))
                if source_key in seen:duplicate_count+=1
                else:nonduplicate.append(k);seen.add(source_key)
            eligible=np.asarray(nonduplicate,dtype=int)
            chosen=np.sort(rng.choice(eligible,size=min(2000,len(eligible)),replace=False))
            if len(chosen)<30:
                raise ValueError('Insufficient lymphoid cells in '+obj['object_id'])
            counts=sparse_dataset(h['layers/counts'])[chosen,:].tocsr()
            m=obs.iloc[chosen].copy();m['original_cell_id']=m.index.astype(str)
            m.index=pd.Index([obj['object_id']+'|'+x for x in m.index])
            m['pilot_coarse_cell_type']=coarse.iloc[chosen].to_numpy()
            m['pilot_assay_scope']=context.get('assay_scope','Source PBMC assay; source-specific selection not assumed absent')
            m['pilot_lineage']=np.where(m['pilot_coarse_cell_type'].isin(['T cells','NK cells']),'T/NK','B/plasma')
            a=ad.AnnData(X=counts,obs=m,var=var.copy())
            a.layers['counts']=counts.copy()
            datasets.append(a)
            selection.append({'object_id':obj['object_id'],'source':obj['catalog']['dataset_id'],
                 'source_cells':len(obs),'eligible_lymphoid_cells':len(eligible),
                 'duplicate_same_source_cells_omitted':duplicate_count,'selected_cells':len(chosen),
                 'original_label_counts':{str(k):int(v) for k,v in m['MajorCellType'].value_counts().items()},
                 'source_sha256':r['matrices']['RNA']['sha256']})
            selection[-1]['assay_scope']=context.get('assay_scope','Source PBMC assay; ascertainment remains study-specific')
            selection[-1]['B_lineage_only_by_verified_assay_scope']=context.get('B_lineage_only_in_pilot',False)
    combined=ad.concat(datasets,join='inner',merge='same',uns_merge='same')
    combined.obs_names_make_unique()
    # A shared gene universe avoids padding unavailable genes with presumed biological zero.
    shared_genes=combined.n_vars
    if shared_genes<1000:raise ValueError('Shared gene universe too small for this pilot; review feature identifiers first')
    if np.any(np.asarray(combined.layers['counts'].sum(1)).ravel()<=0):
        raise ValueError('A sampled cell has no counts in the shared gene universe; review selection rather than silently dropping cells')
    combined.X=combined.layers['counts'].astype(np.float32)
    sc.pp.normalize_total(combined,target_sum=1e4)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined,n_top_genes=3000,flavor='seurat',batch_key='scaid_source',inplace=True)
    protected=['CD3D','CD3E','TRAC','CD8A','MS4A1','CD79A','CD37','CD74','NKG7','GNLY','PRF1','JCHAIN','MZB1']
    nuisance=combined.var_names.str.match(r'^(MT-|RPL[0-9]|RPS[0-9]|HBA[12]$|HBB$)')
    use=(combined.var['highly_variable'].to_numpy()&~nuisance)|combined.var_names.isin(protected)
    combined.var['pilot_embedding_feature']=use
    pca=combined[:,use].copy()
    sc.pp.scale(pca,zero_center=False,max_value=10)
    sc.tl.pca(pca,n_comps=50,zero_center=True,svd_solver='arpack',random_state=SEED)
    combined.obsm['X_pca']=pca.obsm['X_pca'].copy()
    harmon=hm.run_harmony(combined.obsm['X_pca'],combined.obs[['scaid_source']].astype(str),
                         ['scaid_source'],theta=2.0,max_iter_harmony=20,random_state=SEED,device='cpu',verbose=True)
    corrected=np.asarray(harmon.Z_corr)
    if corrected.shape==(50,combined.n_obs):corrected=corrected.T
    if corrected.shape!=(combined.n_obs,50) or not np.isfinite(corrected).all():raise ValueError('Invalid corrected embedding')
    combined.obsm['X_harmony']=corrected.astype(np.float32)
    sc.pp.neighbors(combined,n_neighbors=30,use_rep='X_pca',key_added='before',random_state=SEED)
    sc.tl.umap(combined,neighbors_key='before',random_state=SEED)
    combined.obsm['X_umap_before']=combined.obsm['X_umap'].copy()
    sc.pp.neighbors(combined,n_neighbors=30,use_rep='X_harmony',key_added='after',random_state=SEED)
    sc.tl.umap(combined,neighbors_key='after',random_state=SEED)
    combined.obsm['X_umap_harmony']=combined.obsm['X_umap'].copy()
    metrics={'before':neighbor_metrics(combined.obsm['X_pca'],combined.obs),
             'after':neighbor_metrics(combined.obsm['X_harmony'],combined.obs)}
    t1d=combined.obs['scaid_object_id'].astype(str).str.startswith('T1D/').to_numpy()
    if t1d.sum()>100 and combined.obs.loc[t1d,'scaid_source'].nunique()>1:
        metrics['within_T1D_before']=neighbor_metrics(combined.obsm['X_pca'][t1d],combined.obs.loc[t1d])
        metrics['within_T1D_after']=neighbor_metrics(combined.obsm['X_harmony'][t1d],combined.obs.loc[t1d])
    sc.set_figure_params(dpi=150,fontsize=10)
    palette={'T cells':'#3d7ca6','NK cells':'#b879a8','B cells':'#de8a45','Plasma cells':'#65916f'}
    fig,axes=plt.subplots(2,2,figsize=(12,10),layout='constrained')
    for row,(key,title) in enumerate([('X_umap_before','Before study adjustment'),('X_umap_harmony','After study adjustment')]):
        xy=combined.obsm[key]
        for col,(field,label) in enumerate([('scaid_source','Study'),('pilot_coarse_cell_type','Provided broad cell type')]):
            ax=axes[row,col];levels=sorted(combined.obs[field].astype(str).unique());colors=plt.get_cmap('tab10').colors
            for i,lev in enumerate(levels):
                ix=combined.obs[field].astype(str).eq(lev).to_numpy();color=palette.get(lev,colors[i%len(colors)])
                ax.scatter(xy[ix,0],xy[ix,1],s=2.5,alpha=.7,c=[color],linewidths=0,rasterized=True,label=lev)
            ax.set_title(title+' · '+label,loc='left',fontsize=11);ax.set_xticks([]);ax.set_yticks([])
            ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.spines[['top','right']].set_visible(False)
            ax.legend(markerscale=3,frameon=False,fontsize=8,loc='best')
    fig.savefig(dest/'PBMC_study_integration_overview.png',dpi=240)
    fig.savefig(dest/'PBMC_study_integration_overview.pdf');plt.close(fig)
    # Marker evidence is descriptive and uses the same normalized shared-gene matrix.
    genes=[g for g in protected if g in combined.var_names]
    rownames=sorted(combined.obs['pilot_coarse_cell_type'].unique())
    mat=np.stack([np.asarray(combined[combined.obs['pilot_coarse_cell_type'].eq(k),genes].X.mean(0)).ravel() for k in rownames])
    fig,ax=plt.subplots(figsize=(10,3.6),layout='constrained')
    im=ax.imshow(mat,aspect='auto',cmap='YlGnBu');ax.set_xticks(range(len(genes)),genes,rotation=50,ha='right')
    ax.set_yticks(range(len(rownames)),rownames);fig.colorbar(im,ax=ax,label='Mean log1p expression')
    ax.set_title('Canonical markers in the sampled PBMC lymphoid cells',loc='left',fontsize=11)
    fig.savefig(dest/'PBMC_pilot_marker_expression.png',dpi=240);fig.savefig(dest/'PBMC_pilot_marker_expression.pdf');plt.close(fig)
    provenance={'purpose':'Descriptive feasibility pilot for the next atlas release; no disease effects or independent replication are claimed.',
        'seed':SEED,'source_objects':selection,'sampling':'At most 2,000 provided T/NK/B/plasma cells per PBMC object, sampled without replacement; mixed PBMC/SFMC objects excluded. Verified B-enriched GSE262240 and GSE270142 assays contribute only B/plasma labels; minor non-B labels are not treated as a broad-lineage validation.',
        'deduplication':'In the fixed release-plan object order, each source-accession/cell-identifier pair is assigned to its first eligible object before per-object sampling. This is deterministic object assignment, not study-level balancing.',
        'excluded_monocyte_enriched_assay':'GSE268626 CD and UC objects are excluded from the lymphoid pilot because the study selected peripheral monocytes. Original labels and downloads remain unchanged; minor T/B labels require source-specific marker review.',
        'retained_source_metadata':'All original cell metadata retained. Sample labels are not authenticated donor IDs.',
        'new_qc_filtering':'None. Existing processed source objects and their original doublet/ambient/QC fields were preserved.',
        'shared_gene_count':shared_genes,'embedding_gene_count':int(use.sum()),
        'normalization':'Counts restricted to the shared gene universe, normalized to 10,000 counts per cell, log1p transformed.',
        'features':'Study-aware Seurat-flavor HVGs plus canonical lymphoid markers; mitochondrial/ribosomal/hemoglobin genes excluded from the embedding feature set.',
        'integration':'50 PCA dimensions; Harmony on source accession with theta=2, maximum 20 iterations; UMAP with 30 neighbors.',
        'method_choice':'A bounded conventional pilot permits direct inspection of before/after study effects. A learned large-atlas model is deferred until the study design and annotations are reviewed.',
        'source_condition_confounding':'Some conditions occur in only one source. This embedding does not separate disease effects from study effects.',
        'clinical_stage_scope':'Directory condition labels may include disease-associated risk or treatment groups. Original Disease_subtype and Group labels are retained; a T1D catalogue label is not treated as proof of diagnosed T1D or of a healthy-control comparison.',
        'original_atlases':'Published T/NK and B/plasma atlases were not overwritten.'}
    combined.uns['pilot_provenance_json']=json.dumps(provenance,ensure_ascii=False)
    combined.uns['pilot_technical_diagnostics_json']=json.dumps(metrics,ensure_ascii=False)
    output=dest/'PBMC_lymphoid_integration_pilot.h5ad'
    combined.write_h5ad(output,compression='gzip',compression_opts=1)
    # Independent reopening checks exact axes and count components.
    check=ad.read_h5ad(output)
    if not combined.obs_names.equals(check.obs_names) or not combined.var_names.equals(check.var_names):raise ValueError('Pilot axis roundtrip failed')
    for x,y in [(combined.layers['counts'],check.layers['counts']),(combined.X,check.X)]:
        if sparse.issparse(x):
            if (x!=y).nnz:raise ValueError('Pilot sparse roundtrip differs')
        elif not np.array_equal(x,y):raise ValueError('Pilot data roundtrip differs')
    for k in combined.obsm:
        if not np.array_equal(combined.obsm[k],check.obsm[k]):raise ValueError('Pilot embedding roundtrip differs')
    pd.DataFrame(selection).to_csv(dest/'pilot_object_selection.tsv',sep='\t',index=False)
    combined.obs.to_csv(dest/'pilot_cell_metadata.tsv.gz',sep='\t',compression='gzip')
    pd.crosstab(combined.obs['scaid_source'],combined.obs['pilot_coarse_cell_type']).to_csv(dest/'study_by_cell_type.tsv',sep='\t')
    result={'status':'completed_descriptive_pilot','cells':combined.n_obs,'shared_genes':combined.n_vars,
            'script_sha256':SCRIPT_SIGNATURE,
            'objects':len(selection),'studies':int(combined.obs['scaid_source'].nunique()),
            'output':str(output),'sha256':sha(output),'bytes':output.stat().st_size,
            'validation':{'AnnData_reopen':True,'exact_cell_and_gene_axes':True,'exact_counts_and_X':True,'exact_embeddings':True},
            'provenance':provenance,'technical_diagnostics':metrics,
            'versions':{k:version(k) for k in ['scanpy','anndata','harmonypy','numpy','scipy','scikit-learn']}}
    (dest/'pilot_result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    (dest/'README_INTEGRATION_PILOT.txt').write_text(
        'This is a descriptive pilot for the next atlas version; it does not replace the two atlases in the original article.\n'
        f'It includes {combined.n_obs:,} PBMC cells with existing T/NK/B/plasma labels from {len(selection)} objects and {result["studies"]} studies.\n'
        'At most 2,000 cells are sampled per object; original counts, labels, Sample and QC fields are retained.\n'
        'Original sample labels do not establish authenticated donors. These results assess integration feasibility, source mixing and broad-lineage markers; no case/control or between-disease significance tests are performed.\n'
        'Next, review markers and annotations, confirm independent samples/donors and clinical mappings, and expand integration within the same tissue and comparable study designs.\n')
    print(json.dumps({k:result[k] for k in ['status','cells','shared_genes','objects','studies','output','sha256']},ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');a=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    with open(OUT/'.integration_run.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        result=inventory();print('INVENTORY',result['validated_objects'],'/',result['object_count'],flush=True)
        if a.pilot:pilot()

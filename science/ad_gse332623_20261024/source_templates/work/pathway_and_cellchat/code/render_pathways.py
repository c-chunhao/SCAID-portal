#!/usr/bin/env python3
"""Render all corrected score UMAPs with an exact method/pathway/cell binding."""
from pathlib import Path
import argparse,json,os,textwrap,fcntl
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import h5py,numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from run_pathways import digest,save,METHODS

COORD=None;VALUES=None;OUT=None;OBJECT=None;METHOD=None;NAMES=None

def draw(index):
    pathway_id=NAMES.index[index];name=str(NAMES.iloc[index].pathway_name).replace('_',' ')
    title=f'{METHOD} · {pathway_id}\n'+textwrap.fill(name,width=58)
    values=VALUES[:,index];order=np.argsort(values,kind='stable')
    target=OUT/METHOD/f'{pathway_id}.png';target.parent.mkdir(parents=True,exist_ok=True)
    temporary=target.with_suffix('.png.tmp')
    fig,ax=plt.subplots(figsize=(7.2,6.2),dpi=180)
    points=ax.scatter(COORD[order,0],COORD[order,1],c=values[order],s=2 if len(order)<100000 else .7,
                      linewidths=0,cmap='viridis',rasterized=True)
    ax.set_aspect('equal',adjustable='box');ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2')
    ax.set_xticks([]);ax.set_yticks([]);ax.spines[['top','right']].set_visible(False)
    ax.set_title(title,fontsize=11,pad=13)
    bar=fig.colorbar(points,ax=ax,fraction=.045,pad=.035);bar.set_label(f'{METHOD} score',fontsize=10);bar.ax.tick_params(labelsize=9)
    fig.tight_layout(pad=1.4);fig.savefig(temporary,format='png',facecolor='white');plt.close(fig)
    os.replace(temporary,target)
    with Image.open(target) as image:assert image.size==(1296,1116);image.verify()
    # Directly bind a few plotted values to the source matrix without rounding.
    probes=[0,len(values)//2,len(values)-1]
    return {'object_id':OBJECT,'method':METHOD,'pathway_id':pathway_id,'pathway_name':name,'figure_id':f'{METHOD}__{pathway_id}',
            'title':title.replace('\n',' | '),'path':str(target),'width_px':1296,'height_px':1116,'cell_count':len(values),
            'min':float(values.min()),'max':float(values.max()),'sha256':digest(target),
            'probe_indices':','.join(map(str,probes)),'probe_values':','.join(repr(float(values[i])) for i in probes)}

def main():
    global COORD,VALUES,OUT,OBJECT,METHOD,NAMES
    p=argparse.ArgumentParser();p.add_argument('--pathways',type=Path,required=True);p.add_argument('--workers',type=int,default=4);p.add_argument('--limit',type=int,default=0)
    a=p.parse_args();assert 1<=a.workers<=4
    lock=(a.pathways/'pathway_rendering.lock').open('a')
    fcntl.flock(lock.fileno(),fcntl.LOCK_EX)
    receipt=json.loads((a.pathways/'pathway_completion.json').read_text());assert receipt['status']=='validated'
    if (a.pathways/'pathway_figures_completion.json').is_file():
        done=json.loads((a.pathways/'pathway_figures_completion.json').read_text())
        if done['status']=='validated' and done['pathways_per_method']==(a.limit or receipt['pathway_count']):
            assert digest(a.pathways/'pathway_figure_manifest.tsv')==done['manifest_sha256']
            for row in pd.read_csv(a.pathways/'pathway_figure_manifest.tsv',sep='\t').itertuples():assert digest(row.path)==row.sha256
            print(json.dumps({'status':'already_validated','immutable_PNGs_unchanged':True}));return
    config=json.loads((a.pathways/'scoring_config.json').read_text())
    assert a.limit==0 or receipt['mode']!='formal'
    source=ad.read_h5ad(config['input'],backed='r');COORD=np.asarray(source.obsm['X_umap'],dtype=np.float32)
    assert COORD.shape==(source.n_obs,2) and np.isfinite(COORD).all();source.file.close()
    NAMES=pd.read_csv(receipt.get('gene_overlap_registry',str(a.pathways/'pathway_gene_overlap.tsv')),sep='\t');NAMES=NAMES.loc[NAMES.included].set_index('pathway_id')
    OUT=a.pathways/'figures';OBJECT=receipt['object_id'];records=[]
    for file in receipt['files']:
        METHOD=file['method'];assert METHOD in METHODS
        with h5py.File(file['path'],'r') as f:VALUES=f['X'][...]
        assert VALUES.shape==(len(COORD),len(NAMES)) and np.isfinite(VALUES).all()
        indices=range(min(a.limit,len(NAMES))) if a.limit else range(len(NAMES))
        with ProcessPoolExecutor(max_workers=a.workers,mp_context=mp.get_context('fork')) as pool:
            records.extend(pool.map(draw,indices,chunksize=1))
        del VALUES
        print(json.dumps({'method':METHOD,'rendered':len(indices),'cells':len(COORD)}),flush=True)
    manifest=pd.DataFrame(records);assert not manifest.figure_id.duplicated().any()
    manifest.to_csv(a.pathways/'pathway_figure_manifest.tsv',sep='\t',index=False)
    result={'status':'validated','mode':receipt['mode'],'object_id':OBJECT,'all_cells_plotted':True,'no_cell_subsampling':True,'figures':len(records),'methods':list(METHODS),'pathways_per_method':len(NAMES) if not a.limit else min(a.limit,len(NAMES)),
            'source_score_receipt_sha256':digest(a.pathways/'pathway_completion.json'),'manifest_sha256':digest(a.pathways/'pathway_figure_manifest.tsv'),'output_dimensions':[1296,1116]}
    save(a.pathways/'pathway_figures_completion.json',result)
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()

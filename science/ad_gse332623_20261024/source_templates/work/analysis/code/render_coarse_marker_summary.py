#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def render(root):
 root=Path(root);d=pd.read_csv(root/'annotation_evidence/extended_coarse_marker_summary.tsv',sep='\t');genes=list(pd.unique(d.gene));available=set(d.MajorCellType);order=[g for g in ['T cells','NK cells','B cells','Plasma cells','Monocytes','Dendritic cells','pDC','Platelets','Fibroblasts','Mast cells','Neutrophils','Unassigned']if g in available]
 means=d.pivot(index='MajorCellType',columns='gene',values='mean_log1p_cp10000').loc[order,genes].to_numpy();pct=d.pivot(index='MajorCellType',columns='gene',values='fraction_expressing').loc[order,genes].to_numpy();scaled=(means-means.mean(axis=0))/(means.std(axis=0)+1e-6);ns=d.groupby('MajorCellType').n_records.first();labels=[g+' ('+str(ns[g])+'; provisional)'if g in ['Fibroblasts','Mast cells','Neutrophils']else g for g in order]
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
 fig=plt.figure(figsize=(16,6));gs=fig.add_gridspec(2,2,width_ratios=[1,.035],height_ratios=[.65,.35],wspace=.12,hspace=.35);ax=fig.add_subplot(gs[:,0]);cax=fig.add_subplot(gs[0,1]);lax=fig.add_subplot(gs[1,1]);lax.axis('off')
 xx,yy=np.meshgrid(np.arange(len(genes)),np.arange(len(order)));dots=ax.scatter(xx,yy,s=12+185*pct,c=scaled,cmap='RdBu_r',vmin=-2,vmax=2,linewidth=.25,edgecolors='grey');ax.set_xticks(np.arange(len(genes)),genes,rotation=60,ha='right',fontsize=9);ax.set_yticks(np.arange(len(order)),labels,fontsize=9);ax.invert_yaxis();fig.colorbar(dots,cax=cax,label='Mean RNA expression\n(z-score per gene)')
 handles=[ax.scatter([],[],s=12+185*f,c='grey',linewidth=.25,edgecolors='grey',label=str(int(f*100))+'%')for f in [.1,.5,.9]];lax.legend(handles=handles,title='Fraction\nexpressing',loc='center',frameon=False,scatterpoints=1,labelspacing=1.3,handletextpad=.8,fontsize=8,title_fontsize=8)
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(root/'plots'/('RNA_marker_dotplot.'+ext))
 plt.close(fig)

if __name__=='__main__':render('@PROJECT@/00.corrected_analysis/GSE332623/integration')

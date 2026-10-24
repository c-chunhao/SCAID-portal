#!/usr/bin/env python3
"""Extend coarse marker evidence/plot without changing RNA/ADT matrices or metadata."""
from pathlib import Path
import json,hashlib
import h5py,numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('@PROJECT@/00.corrected_analysis/GSE332623/integration')
GENES=['CD3D','TRAC','NKG7','GNLY','MS4A1','CD79A','MZB1','JCHAIN','LYZ','LST1','CD1C','LILRA4','IL3RA','PF4','PPBP','HBB','COL1A1','DCN','COL3A1','KIT','TPSAB1','CSF3R','FCGR3B']
def cat(g):
 if isinstance(g,h5py.Group):
  c=g['codes'][:];assert np.all(c>=0);return g['categories'].asstr()[:][c]
 return g.asstr()[:]
with h5py.File(ROOT/'RNA_normalized_annotated.h5ad','r')as f:
 names=f['var'][f['var'].attrs['_index']].asstr()[:];genes=[x for x in GENES if x in set(names)];wanted=np.full(len(names),-1,dtype=np.int32);wanted[pd.Index(names).get_indexer(genes)]=np.arange(len(genes))
 labels=cat(f['obs/MajorCellType']);libs=cat(f['obs/library_id']);levels=pd.Index(pd.unique(labels));codes=levels.get_indexer(labels);sizes=np.bincount(codes,minlength=len(levels));means=np.zeros((len(levels),len(genes)));nexpr=np.zeros_like(means)
 ptr=f['X/indptr'][:];ng=len(genes)
 for start in range(0,len(labels),10000):
  end=min(len(labels),start+10000);lo=int(ptr[start]);hi=int(ptr[end]);feature=f['X/indices'][lo:hi];value=f['X/data'][lo:hi];col=wanted[feature];keep=(col>=0)&(value>0);row=np.repeat(np.arange(start,end),np.diff(ptr[start:end+1]))[keep];flat=codes[row]*ng+col[keep]
  means+=np.bincount(flat,weights=value[keep],minlength=len(levels)*ng).reshape(len(levels),ng);nexpr+=np.bincount(flat,minlength=len(levels)*ng).reshape(len(levels),ng)
means/=sizes[:,None];pct=nexpr/sizes[:,None]
rows=[{'MajorCellType':g,'n_records':int(sizes[i]),'gene':x,'mean_log1p_cp10000':float(means[i,j]),'fraction_expressing':float(pct[i,j])}for i,g in enumerate(levels)for j,x in enumerate(genes)]
pd.DataFrame(rows).to_csv(ROOT/'annotation_evidence/extended_coarse_marker_summary.tsv',sep='\t',index=False)
order=[x for x in ['T cells','NK cells','B cells','Plasma cells','Monocytes','Dendritic cells','pDC','Platelets','Fibroblasts','Mast cells','Neutrophils','Unassigned']if x in levels];ix=levels.get_indexer(order);mu=means[ix];pr=pct[ix];scaled=(mu-mu.mean(axis=0))/(mu.std(axis=0)+1e-6)
label=[];review=[]
for g in order:
 n=int(sizes[levels.get_loc(g)]);source=pd.Series(libs[labels==g]).value_counts();status='main_PBMC_reference_prediction';note='Computational reference/RNA evidence; no external expert confirmation claimed.'
 if g=='Fibroblasts':status='source_dominated_provisional_prediction';note='All172records from NEA7_ADT2/source-dominated cluster23. Model and collagen/DCN RNA support retained; cell identity/tissue origin requires team review, not a validated new finding.'
 if g in ['Mast cells','Neutrophils']:status='very_rare_provisional_prediction';note='Retain original computational prediction for review; no discrete fine subtype or new-cell-population claim. The single neutrophil prediction lacks FCGR3B support.'
 if g=='Unassigned':status='unresolved_computational_annotation';note='All records retained; raw reference labels and uncertainty remain available.'
 review.append({'MajorCellType':g,'n_records':n,'n_source_libraries':len(source),'dominant_library':source.index[0],'dominant_source_fraction':float(source.iloc[0]/n),'review_status':status,'interpretation':note})
 label.append(g+' ('+str(n)+'; provisional)' if g in ['Fibroblasts','Mast cells','Neutrophils']else g)
pd.DataFrame(review).to_csv(ROOT/'annotation_evidence/coarse_annotation_review.tsv',sep='\t',index=False)
from render_coarse_marker_summary import render
render(ROOT)
report={'status':'completed','source':'Final RNA.X log1p(cp10000); streamed CSR rows10000at a time; only group-by-marker summary dense','n_records':len(labels),'n_genes':len(genes),'coarse_markers_added':['COL1A1','DCN','COL3A1','KIT','TPSAB1','CSF3R','FCGR3B'],'RNA_ADT_or_obs_modified':False,'review_table':str(ROOT/'annotation_evidence/coarse_annotation_review.tsv'),'helper_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(ROOT/'provenance/extended_marker_review.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))

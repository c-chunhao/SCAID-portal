import os
"""Retain original A-D vector marker marks; replace legacy E-H statistics.
E-H are re-laid as four labelled rows of a descriptive heatmap, with complete
zero-corrected source specimen distributions in the accompanying detailed figure.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pymupdf as fitz
from PIL import Image
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures'
old=Path(os.environ['SCAID_ORIGINAL_FIGURE3_PDF'])
# Data derived by the zero-complete script, rather than reading plotted heights.
z=pd.read_csv(T/'all_layers_condition_descriptive_medians.tsv',sep='\t');cts=[('T_NK','CD4_T'),('T_NK','CD8_T'),('B_Plasma','B cells'),('B_Plasma','Plasma cells')]
conditions=sorted(z.Condition.unique());arr=np.array([z[(z.Lineage==lin)&z.Layer.eq('layer1')&z.Cell_type.eq(ct)].set_index('Condition').Median_fraction.reindex(conditions).to_numpy()for lin,ct in cts])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(6.5,2.65));im=ax.imshow(arr,aspect='auto',cmap='YlGnBu',vmin=0,vmax=float(arr.max()));ax.set_yticks(range(4),['E  CD4 T','F  CD8 T','G  B cells','H  Plasma cells']);ax.set_xticks(range(len(conditions)),conditions,rotation=90);ax.set_title('Zero-complete specimen composition (descriptive)',loc='left',fontsize=9);cb=fig.colorbar(im,ax=ax,fraction=.025,pad=.015);ticks=np.linspace(0,float(arr.max()),4);cb.set_ticks(ticks,labels=[f'{100*t:.0f}%'for t in ticks]);cb.set_label('Median fraction',fontsize=8);fig.text(.5,.015,'Fixed immune denominator; source/tissue confounded; no pooled ANOVA.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.06,1,1]);sub=F/'Figure3_EH_summary_heatmap.pdf';fig.savefig(sub);plt.close(fig)
original=fitz.open(old)
# An isolated 3.53-pt residual x from the old indexed layout was left in panel C.
# Remove that text-only artifact; never redact vector paths/dot positions/colours.
removed=[]
for b in original[0].get_text('dict')['blocks']:
 for ln in b.get('lines',[]):
  for sp in ln['spans']:
   if sp['bbox'][1]<287 and sp['size']<6.999 and sp['text'].strip():
    assert sp['text'].strip()=='x',f'Unreviewed scientific text {sp["text"]}'
    original[0].add_redact_annot(fitz.Rect(sp['bbox']),fill=None);removed.append({'text':sp['text'],'bbox':sp['bbox'],'size':sp['size']})
if removed:original[0].apply_redactions(images=0,graphics=0,text=0)
d=fitz.open();pg=d.new_page(width=468,height=487.8);pg.show_pdf_page(fitz.Rect(0,0,468,287),original,0,clip=fitz.Rect(0,0,468,287))
summary=fitz.open(sub);pg.show_pdf_page(fitz.Rect(0,293,468,483.8),summary,0);pg.insert_text((8,484.8),'A-D marker labels and scales: Supplementary Table S3 and Figure3_marker_panels.pdf.',fontsize=7)
p=F/'Figure3_recomputed.pdf';d.save(p,garbage=4,deflate=True);pg.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False).save(F/'Figure3_recomputed.png');im=Image.open(F/'Figure3_recomputed.png');im.save(F/'Figure3_recomputed.tiff',dpi=(600,600),compression='tiff_lzw');summary.close();original.close();d.close()
# Full A-D originals retain all gene/cell labels and original quantitative legends.
source=Path(os.environ['SCAID_ORIGINAL_MARKER_PANEL_DIR']);marker=fitz.open()
for stem in ['Fig3A_TNK_layer1_MarkerPlot','Fig3B_CD4T_MarkerPlot','Fig3C_CD8T_MarkerPlot','Fig3D_BPlasma_layer1_MarkerPlot']:
 x=fitz.open(source/(stem+'.pdf'));marker.insert_pdf(x);x.close()
marker.save(F/'Figure3_marker_panels.pdf',garbage=4,deflate=True);marker.close()
(B/'figure3_composition_provenance.json').write_text(json.dumps({'A_D_source':str(old),'A_D_clip_pt':[0,0,468,287],'scientific_A_D_vector_marks':'Source dot positions, sizes and colour vector paths preserved; no raw expression or marker statistic changed. One isolated residual text x from the previous indexed layout removed without redacting vector paths.','artifact_removed':removed,'E_H':'New native descriptive median-fraction heatmap derived from complete-zero specimen tables. No legacy ANOVA/F/P retained. Full distributions in Figure3_EH_recomputed.pdf.','figure_final_inches':[6.5,487.8/72]},indent=2)+'\n')
print('Figure3 assembled')

from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
B=Path(__file__).resolve().parents[1];T=B/'tables';F=B/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42})
IFN=['ISG15','IFI6','IFI44L','MX1','OAS1','IFIT1','IFIT3','RSAD2','ISG20','XAF1','LY6E','IFI44']
tests=pd.read_csv(T/'GSE174188_unique_donor_HC3_composition_IFN_tests.tsv',sep='\t');comp=tests[tests.analysis.eq('primary_unique_batch')&tests.kind.eq('composition')]
score=pd.read_csv(T/'GSE174188_target_cell_IFN_by_unique_donor.tsv',sep='\t');units=pd.read_csv(T/'GSE174188_frozen_one_batch_per_donor.tsv',sep='\t');assert units.ind_cov.is_unique
fig=plt.figure(figsize=(6.5,7.8));gs=gridspec.GridSpec(3,6,figure=fig,height_ratios=[1.05,1.12,1.25],hspace=.78,wspace=1.8)
ax=fig.add_subplot(gs[0,:3]);ax.axis('off');ax.set_title('A  Verified source individuals',loc='left')
items=[('Source processed PBMCs','1,263,676 cells; 261 individuals'),('Eligible case/reference cohorts','Processing cohorts 2–4'),('Frozen independent observation','1 source batch per individual'),('Primary comparison','71 references / 162 SLE cases')]
for i,(title,detail)in enumerate(items):
 y=.88-i*.245;ax.text(.02,y,title,transform=ax.transAxes,fontsize=8.5,fontweight='bold',va='top');ax.text(.02,y-.105,detail,transform=ax.transAxes,fontsize=8,va='top');
 if i<3:ax.plot([.02,.98],[y-.183,y-.183],transform=ax.transAxes,color='#dedede',linewidth=.5)
ax=fig.add_subplot(gs[0,3:]);yy=np.arange(len(comp));e=comp.adjusted_coefficient_SLE_minus_HD*100;lo=comp.HC3_95CI_low*100;hi=comp.HC3_95CI_high*100
ax.errorbar(e,yy,xerr=np.vstack([e-lo,hi-e]),fmt='o',color='#526e80',capsize=2,markersize=3);ax.axvline(0,color='#aaaaaa',linewidth=.7);ax.set_yticks(yy,comp.cell_type);ax.invert_yaxis();ax.set_xlabel('Adjusted SLE − reference\nfraction (percentage points)');ax.set_title('B  Source cell composition',loc='left');ax.spines[['top','right']].set_visible(False)
ax=fig.add_subplot(gs[1,:]);rng=np.random.default_rng(20260930);positions=[];values=[];colors=[];counts=[]
for j,ct in enumerate(['cM','T4','B']):
 s=score[score.cg_cov.eq(ct)&score.n_cells.ge(20)]
 for arm,pos,c in [('HD',j*3+1,'#526e80'),('SLE',j*3+2,'#a96343')]:
  v=s.loc[s.condition.eq(arm),'IFN12_log1p10k_mean'].to_numpy();assert s[s.condition.eq(arm)].ind_cov.is_unique
  positions.append(pos);values.append(v);colors.append(c);counts.append(len(v));ax.scatter(pos+rng.uniform(-.16,.16,len(v)),v,s=3.5,color=c,alpha=.3,linewidths=0,zorder=1)
 bp=ax.boxplot(values[-2:],positions=positions[-2:],widths=.45,showfliers=False,patch_artist=True,medianprops={'color':'black','linewidth':.8},whiskerprops={'linewidth':.65},capprops={'linewidth':.65},boxprops={'linewidth':.65})
 for box,c in zip(bp['boxes'],colors[-2:]):box.set_facecolor(c);box.set_alpha(.35)
ax.set_xticks(positions,[('Ref.'if k%2==0 else'SLE')+'\nn='+str(n)for k,n in enumerate(counts)]);ax.set_ylabel('Donor mean IFN12 log1p(10k)');ax.set_title('C  Interferon signature within source-author cell types',loc='left',pad=27);ax.spines[['top','right']].set_visible(False)
for j,ct in enumerate(['cM','T4','B']):ax.text(j*3+1.5,1.015,{'cM':'Classical monocytes','T4':'CD4 T cells','B':'B cells'}[ct],ha='center',va='bottom',transform=ax.get_xaxis_transform(),fontsize=8)
watch=[]
for j,(ct,letter)in enumerate(zip(['cM','T4','B'],'DEF')):
 ax=fig.add_subplot(gs[2,2*j:2*j+2]);d=pd.read_csv(T/f'GSE174188_{ct}_primary_unique_batch_DESeq2.tsv',sep='\t').set_index('gene').reindex(IFN);assert d.log2FoldChange.notna().all()and d.padj.lt(.05).all();watch.append(d.reset_index().assign(cell_type=ct))
 e=d.log2FoldChange.to_numpy();se=d.lfcSE.to_numpy();yy=np.arange(len(IFN));ax.errorbar(e,yy,xerr=1.96*se,fmt='o',color='#a96343',capsize=1.5,markersize=2.6,elinewidth=.65);ax.set_yticks(yy,IFN);ax.invert_yaxis();ax.axvline(0,color='#aaaaaa',linewidth=.6);ax.set_xlim(-.1,3.55);ax.set_xticks([0,1,2,3]);ax.set_xlabel('Pseudobulk log₂ fold change');ax.set_title(f'{letter}  '+{'cM':'cM','T4':'CD4 T','B':'B'}[ct],loc='left');ax.spines[['top','right']].set_visible(False)
fig.text(.5,.017,'Age, sex, ancestry and processing cohort adjusted; donor-level observations.',ha='center',fontsize=8)
fig.subplots_adjust(left=.115,right=.985,top=.95,bottom=.10)
fig.savefig(F/'Figure5_verified_donor.pdf');fig.savefig(F/'Figure5_verified_donor.png',dpi=600);fig.savefig(F/'Figure5_verified_donor.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'});plt.close(fig)
pd.concat(watch,ignore_index=True).to_csv(T/'Figure5_IFN12_pseudobulk_display_values.tsv',sep='\t',index=False)
(B/'figure5_native_provenance.json').write_text(json.dumps({'figure_inches':[6.5,7.8],'A':'Frozen source individual/batch selection','B':'Unique-source-individual HC3-adjusted fraction differences, 95% CI','C':'One score per biological individual and target cell type; source processed raw count normalization, not pooled cell replicates','D_F':'Predefined 12 IFN genes, cell-type pseudobulk log2 fold change and Wald 95% CI; all shown genes BH<0.05 within their lineage','source_cells':'Full archived processed source MTX, not a new primary read alignment/QC run','claims':'Within-study observational SLE/reference association; reproduces a known IFN pattern, not novel cross-disease mechanism or external independent validation.'},indent=2)+'\n');print('Figure5 written')

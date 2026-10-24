#!/usr/bin/env python3
"""Collect per-library QC receipts and flags without opening expression matrices."""
import csv
from datetime import datetime, timezone, timedelta
import json
import os
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path('@PROJECT@/00.corrected_analysis/GSE332623')
WORK = Path('@WORK@')


def main():
    output = BASE / 'qc/summary'; output.mkdir(exist_ok=True)
    manifest = list(csv.DictReader((BASE / 'capture_manifest.tsv').open(), delimiter='\t'))
    rows, failures = [], []
    flags = ['qc_low_umi','qc_high_umi','qc_low_genes','qc_high_genes','qc_high_mt','passes_basic_qc','passes_QC','ambient_review']
    for entry in manifest:
        directory = BASE / 'qc' / entry['library_id']; path = directory / 'receipt.json'
        if not path.exists(): continue
        record = json.loads(path.read_text())
        if record.get('status') != 'validated': continue
        progress_path=directory/'progress.json'
        if progress_path.exists():
            progress=json.loads(progress_path.read_text())
            if progress.get('status')!='validated':
                progress.update(status='validated',stage='RNA_QC_bridge_validated',
                                receipt=str(path),clinical_status='Unresolved')
                temp=progress_path.with_suffix('.summary-tmp')
                temp.write_text(json.dumps(progress,ensure_ascii=False,indent=2)+'\n');temp.replace(progress_path)
        frame = pd.read_csv(directory / 'per_cell_QC.tsv.gz', sep='\t',
                            usecols=['cell_id','nCount_RNA','nFeature_RNA','percent_mt','percent_ribo','percent_hb',
                                     'doublet_class','ambient_contamination',*flags])
        for flag in flags:
            frame[flag] = frame[flag].astype(str).str.lower().eq('true')
        assert frame.cell_id.is_unique and len(frame) == int(record['n_input'])
        assert int(frame.passes_QC.sum()) == int(record['n_retained'])
        assert int((frame.doublet_class == 'doublet').sum()) == int(record['n_doublets'])
        assert ((frame.passes_QC) == (frame.passes_basic_qc & frame.doublet_class.eq('singlet'))).all()
        row = {'library_id':entry['library_id'],'pool_id':entry['pool_id'],'rna_gsm':entry['rna_gsm'],
               'adt_gsm':entry['adt_gsm'],'n_source_barcodes':len(frame),
               'n_basic_retained':int(frame.passes_basic_qc.sum()),'n_basic_excluded':int((~frame.passes_basic_qc).sum()),
               'n_scDblFinder_doublets':int((frame.doublet_class == 'doublet').sum()),
               'n_final_retained':int(frame.passes_QC.sum()),'n_ambient_review_retained':int((frame.passes_QC & frame.ambient_review).sum()),
               'ambient_is_filter':False,'doublet_fraction_among_basic':float(frame.doublet_class.eq('doublet').sum()/frame.passes_basic_qc.sum()),
               'RNA_median_UMI':float(frame.nCount_RNA.median()),'RNA_median_genes':float(frame.nFeature_RNA.median()),
               'RNA_median_mt_percent':float(frame.percent_mt.median()),'RNA_median_ribo_percent':float(frame.percent_ribo.median()),
               'RNA_hb_99_percentile':float(frame.percent_hb.quantile(.99)),
               'retained_low_RNA_500genes':int((frame.passes_QC & frame.nFeature_RNA.lt(500)).sum()),
               'clinical_assignment':'Unresolved','source_sha256':entry['input_sha256'],
               'bridge_sha256':record['bridge_sha256']}
        row.update({'threshold_'+key:value for key,value in record['thresholds'].items() if isinstance(value,(int,float))})
        row.update({flag+'_count':int(frame[flag].sum()) for flag in flags[:5]})
        row['review_flags'] = ';'.join([name for name,condition in [
            ('basic_retention_below_80pct',row['n_basic_retained']/len(frame) < .8),
            ('doublet_fraction_above_45pct',row['doublet_fraction_among_basic'] > .45),
            ('ambient_review_above_10pct',row['n_ambient_review_retained']/row['n_final_retained'] > .1)] if condition])
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(output / 'QC_summary_per_library.tsv',sep='\t',index=False)
    report = {'status':'complete' if len(rows)==27 else 'partial','updated':datetime.now(timezone(timedelta(hours=8))).isoformat(),
              'n_libraries_validated':len(rows),'n_libraries_expected':27,
              'n_source_barcode_records':int(summary.n_source_barcodes.sum()) if rows else 0,
              'n_basic_retained':int(summary.n_basic_retained.sum()) if rows else 0,
              'n_predicted_doublets':int(summary.n_scDblFinder_doublets.sum()) if rows else 0,
              'n_final_retained':int(summary.n_final_retained.sum()) if rows else 0,
              'n_ambient_review_retained':int(summary.n_ambient_review_retained.sum()) if rows else 0,
              'ambient_is_filter':False,'clinical_assignment':'Unresolved',
              'source_barcodes_are_donor_counts':False,'per_cell_table_invariants_passed':True,
              'source_raw_counts_preserved':'RNA/ADT integer source counts; decontX separate fractional matrix',
              'libraries_requiring_review':[{key:row[key] for key in ('library_id','review_flags')} for row in rows if row['review_flags']]}
    (output / 'QC_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    if not rows: print(json.dumps(report)); return
    labels = summary.library_id.str.replace('GSE332623_TotalseqB_','',regex=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, ax = plt.subplots(figsize=(11,max(4,len(rows)*.28)))
    y=np.arange(len(rows)); bottom=np.zeros(len(rows))
    for column,color,label in [('n_final_retained','#297f86','Retained'),('n_scDblFinder_doublets','#c57969','Predicted doublets'),('n_basic_excluded','#c3c8ca','RNA metric outliers')]:
        value=summary[column].to_numpy();ax.barh(y,value,left=bottom,color=color,label=label);bottom+=value
    ax.set_yticks(y,labels);ax.invert_yaxis();ax.set_xlabel('Source barcode records');ax.legend(frameon=False,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.1));fig.tight_layout();fig.savefig(output/'QC_retention_by_library.png',dpi=180);plt.close(fig)
    fig, axes=plt.subplots(1,3,figsize=(13,max(4,len(rows)*.25)),sharey=True)
    for ax,metric,lower,upper,label in zip(axes,['RNA_median_UMI','RNA_median_genes','RNA_median_mt_percent'],
        ['threshold_UMI_min','threshold_genes_min',None],['threshold_UMI_max','threshold_genes_max','threshold_mt_max'],
        ['RNA UMIs','Detected RNA genes','Mitochondrial RNA (%)']):
        ax.scatter(summary[metric],y,color='#297f86',s=17,label='Median')
        if lower:ax.scatter(summary[lower],y,marker='|',color='#b44932',s=40,label='Lower limit')
        ax.scatter(summary[upper],y,marker='|',color='#b44932',s=40,label='Upper limit')
        if lower:ax.set_xscale('log')
        ax.set_xlabel(label);ax.set_yticks(y,labels);ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    axes[0].legend(frameon=False,fontsize=8);fig.tight_layout();fig.savefig(output/'QC_thresholds_by_library.png',dpi=180);plt.close(fig)
    print(json.dumps(report,ensure_ascii=False))


if __name__=='__main__': main()

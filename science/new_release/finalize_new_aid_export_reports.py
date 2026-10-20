#!/usr/bin/env python3
"""Add explicit content semantics to completed reports and regenerate the manifest."""
import importlib.util, json
from pathlib import Path

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
def main():
    for p in (WORK/'matrices/objects').glob('*/export_report.json'):
        r=json.loads(p.read_text())
        if r['status']!='validated':continue
        for method,item in r['matrices'].items():
            if method=='RNA':
                if 'counts_available' not in item:
                    if not item['independent_validation']['counts']['all_values_nonnegative_integer']:
                        raise ValueError('Missing integer count verification')
                    item.update(counts_available=True,
                        raw_count_status='integer_count_layer_preserved_from_processed_source',
                        expression_status='source_normalized_data_with_integer_count_layer',
                        content_description='RNA expression with the source integer count layer.',
                        source_processed_counts_semantics='Not separately supplied; the original RNA counts layer contains nonnegative integers.')
                r['counts_available']=item['counts_available'];r['raw_count_status']=item['raw_count_status'];r['expression_status']=item['expression_status']
            else:
                item.update(counts_available=False,raw_count_status='not_applicable_pathway_scores',
                            expression_status='original_unscaled_pathway_score',content_description=method+' KEGG pathway scores from the original data layer.')
        q=p.with_suffix('.json.tmp');q.write_text(json.dumps(r,ensure_ascii=False,indent=2));q.replace(p)
    spec=importlib.util.spec_from_file_location('export_batch',CODE_ROOT/'export_new_aid_batch.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    manifest=module.refresh_manifest()
    print(json.dumps({k:manifest[k] for k in ['validated_object_count','file_count','total_validated_bytes']}))

if __name__=='__main__':main()

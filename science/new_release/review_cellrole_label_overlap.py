from pathlib import Path
import concurrent.futures,json,itertools
import pymupdf as fitz
from runtime_paths import WORK, PLAN_PATH, FIGURE_ROOT, INPUT_ROOT, configured
HERE=WORK/'cellchat'

def inspect(task):
 row,groups=task;p=Path(row['source_pdf']);doc=fitz.open(p);page=doc[0];spans=[]
 for block in page.get_text('dict')['blocks']:
  for line in block.get('lines',[]):
   for span in line['spans']:
    if span['text'].strip() in groups:spans.append({'label':span['text'].strip(),'bbox':span['bbox']})
 overlaps=[]
 for a,b in itertools.combinations(spans,2):
  if a['label']==b['label']:continue
  ar,br=fitz.Rect(a['bbox']),fitz.Rect(b['bbox']);inter=ar&br
  ratio=inter.get_area()/min(ar.get_area(),br.get_area()) if not inter.is_empty else 0
  if ratio>=.15:overlaps.append({'labels':[a['label'],b['label']],'overlap_fraction_of_smaller':round(ratio,3),'bbox_a':a['bbox'],'bbox_b':b['bbox']})
 doc.close()
 return {'dataset_object':row['dataset_object'],'source_pdf':str(p),'path':row['path'],'class':row['class'],'matched_label_count':len(spans),'overlaps':overlaps} if overlaps else None

if __name__=='__main__':
 rows=[json.loads(s) for s in (HERE/'asset_index.runtime.jsonl').read_text().splitlines()]
 names={}
 for r in rows:
  if r['class']=='Single' and Path(r['source_pdf']).stem.startswith('Count-'):names.setdefault(r['dataset_object'],set()).add(Path(r['source_pdf']).stem[6:].replace('_','/'))
 tasks=[]
 for r in rows:
  if r['class']=='SignalingCellRole' or (r['class']=='Overview' and r['name']=='cellRole.png'):tasks.append((r,names.get(r['dataset_object'],set())))
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(inspect,tasks,chunksize=8))
 flagged=[r for r in results if r];out={'reviewed_scatter_pdf_count':len(tasks),'flagged_text_overlap_pdf_count':len(flagged),'definition':'Intersection >=15% of the smaller cell-group label bounding box. Flags locate candidates; biological identities and coordinates are not changed.','visual_confirmation_needed':True,'flagged':flagged}
 (HERE/'cellrole_label_overlap_candidates.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
 print(json.dumps({'reviewed':len(tasks),'flagged':len(flagged),'largest_overlap_examples':sorted(flagged,key=lambda r:max(x['overlap_fraction_of_smaller'] for x in r['overlaps']),reverse=True)[:4]},ensure_ascii=False))

#!/usr/bin/env python3
"""Independent file checks and runtime import index for prepared CellChat assets."""
from pathlib import Path
import collections,hashlib,json,os,concurrent.futures
import pymupdf as fitz
from PIL import Image
from runtime_paths import WORK, PLAN_PATH, FIGURE_ROOT, INPUT_ROOT, configured
HERE=WORK/'cellchat'
RUNTIME=FIGURE_ROOT

def check(row):
 errors=[];source=Path(row['source_pdf']);p=Path(row['path'])
 try:
  if hashlib.sha256(source.read_bytes()).hexdigest()!=row['source_sha256']:errors.append('source_hash_changed')
  if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:errors.append('asset_hash_mismatch')
  with Image.open(p) as im:
   if list(im.size)!=[row['image_width'],row['image_height']]:errors.append('png_dimension_mismatch')
   im.verify()
  if max(row['image_width'],row['image_height'])!=1600:errors.append('render_size_not_uniform')
  if len(row['name'])>255 or len(row['relative_path'])>512:errors.append('db_field_too_long')
  if row['pdf_pages']!=1:errors.append('unexpected_source_page_count')
  # Independently verify that PDF content corrections modify no other text.
  if row['text_corrections']:
   doc=fitz.open(source);page=doc[0];before=page.get_text(clip=fitz.INFINITE_RECT())
   expected=before.replace('Unknonw','Unknown')
   if row['class']=='Overview' and source.stem=='Weight':expected=expected.replace('Number of interactions','Interaction strength')
   for x in page.get_contents():
    data=doc.xref_stream(x).replace(b'(Unknonw)',b'(Unknown)')
    if row['class']=='Overview' and source.stem=='Weight':
     for old in [b'[(Number of inter) 10 (actions)] TJ',b'(Number of interactions) Tj']:data=data.replace(old,b'(Interaction strength) Tj')
    doc.update_stream(x,data)
   page=doc.reload_page(page);after=page.get_text(clip=fitz.INFINITE_RECT())
   if after!=expected:errors.append('text_change_beyond_approved_replacements')
   doc.close()
 except Exception as e:errors.append(f'{type(e).__name__}: {e}')
 return {'source_pdf':str(source),'path':str(p),'errors':errors} if errors else None

def main():
 rows=[json.loads(s) for s in (HERE/'asset_index.jsonl').read_text().splitlines()]
 # Publish spelling in menu filenames too, while retaining the exact source path.
 renames=[]
 for row in rows:
  old=Path(row['path']);name=old.name.replace('Unknonw','Unknown')
  if old.name!=name:
   new=old.with_name(name)
   if new.exists():raise RuntimeError(f'Display spelling rename target already exists: {new}')
   os.rename(old,new);row.update(name=name,path=str(new),image=str(new),source_filename=Path(row['source_pdf']).name)
   renames.append({'from':str(old),'to':str(new),'reason':'Unknonw -> Unknown spelling only'})
 (HERE/'asset_index.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in rows))
 (HERE/'display_filename_corrections.json').write_text(json.dumps(renames,ensure_ascii=False,indent=2))
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:results=list(pool.map(check,rows))
 errors=[r for r in results if r]
 runtime=[]
 for row in rows:
  r={**row,'stage_path':row['path']};target=RUNTIME/row['relative_path']/row['name'];r.update(path=str(target),image=str(target));runtime.append(r)
 (HERE/'asset_index.runtime.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in runtime))
 # Refresh correction paths after menu-name changes.
 (HERE/'display_label_corrections.jsonl').write_text(''.join(json.dumps({k:r[k] for k in ['dataset_object','class','source_pdf','path','source_sha256','sha256','text_corrections']},ensure_ascii=False)+'\n' for r in rows if r['text_corrections']))
 manifest=json.loads((HERE/'manifest.json').read_text())
 manifest.update(runtime_asset_root=str(RUNTIME),runtime_index=str(HERE/'asset_index.runtime.jsonl'),display_filename_correction_count=len(renames))
 (HERE/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
 counts=collections.Counter(r['dataset_object'] for r in rows)
 qc={'passed':not errors,'asset_count':len(rows),'object_count':len(counts),'source_hashes_unchanged':not any('source_hash_changed' in x['errors'] for x in errors),
  'all_png_signatures_crcs_hashes_dimensions_valid':not errors,'pdf_text_unchanged_except_approved_spelling_and_weight_title':not any('text_change_beyond_approved_replacements' in x['errors'] for x in errors),
  'runtime_index_prepared':True,'deployed':False,'error_count':len(errors),'errors':errors}
 (HERE/'independent_QA.json').write_text(json.dumps(qc,ensure_ascii=False,indent=2))
 print(json.dumps(qc,ensure_ascii=False))
 if errors:raise SystemExit(1)

if __name__=="__main__":main()

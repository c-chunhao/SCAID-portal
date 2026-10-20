#!/usr/bin/env python3
"""Prepare website PNG copies of existing CellChat PDFs without modifying sources.
Semantic mapping follows the generating CellChatFig.R directory names. Only two
text corrections are permitted: Unknonw -> Unknown and aggregate Weight title
Number of interactions -> Interaction strength. No network values are changed.
"""
from pathlib import Path
import argparse, collections, concurrent.futures, hashlib, json, os, time
import pymupdf as fitz
from PIL import Image, ImageDraw, ImageFont

from runtime_paths import WORK, PLAN_PATH, FIGURE_ROOT, INPUT_ROOT, configured
HERE=WORK/'cellchat'
SOURCE=configured('SCAID_CELLCHAT_INPUT_ROOT',INPUT_ROOT/'02.CellChat')
PLAN=PLAN_PATH
CLASSES=('Overview','LRcircle','LRcontribution','SignalingCellRole','SignalingChord','SignalingCircle','SignalingHt','SignalingRole','Single')
MAX_SIDE=1600

def qa_font(size):
    try:
        return ImageFont.truetype(os.environ.get('SCAID_QA_FONT', 'DejaVuSans.ttf'), size)
    except OSError:
        return ImageFont.load_default(size=size)

def digest(data):return hashlib.sha256(data).hexdigest()

def prepare(task):
    obj, source, family = task
    f=Path(source)
    output=HERE/'assets'/obj/'CellChat'/family/(f.stem+'.png')
    row={'dataset_object':obj,'class':family,'source_pdf':source,'name':output.name,
         'relative_path':f'{obj}/CellChat/{family}','path':str(output),'image':str(output),
         'page_number':1,'is_pdf':'False','status':'pending'}
    try:
        raw=f.read_bytes();row['source_sha256']=digest(raw);row['source_bytes']=len(raw)
        doc=fitz.open(stream=raw,filetype='pdf');row['pdf_pages']=len(doc)
        if len(doc)!=1:raise ValueError('Expected exactly one PDF page')
        page=doc[0]
        raw_text=page.get_text(clip=fitz.INFINITE_RECT())
        if not raw_text.strip():raise ValueError('PDF has no extractable text; requires separate review')
        source_words=page.get_text('words',clip=fitz.INFINITE_RECT())
        outside=[]
        for w in source_words:
            r=fitz.Rect(w[:4]);q=page.rect
            if r.x0 < -.5 or r.y0 < -.5 or r.x1 > q.width+.5 or r.y1 > q.height+.5:
                outside.append({'text':w[4],'bbox':list(r)})
        corrections=[];replacement_count=0;title_count=0
        for x in page.get_contents():
            stream=doc.xref_stream(x)
            typo_count=stream.count(b'(Unknonw)')
            if typo_count:
                stream=stream.replace(b'(Unknonw)',b'(Unknown)');replacement_count+=typo_count
            if family=='Overview' and f.stem=='Weight':
                for old in [b'[(Number of inter) 10 (actions)] TJ',b'(Number of interactions) Tj']:
                    count=stream.count(old)
                    if count:
                        stream=stream.replace(old,b'(Interaction strength) Tj');title_count+=count
            doc.update_stream(x,stream)
        # Reload the page to force text-cache regeneration after stream editing.
        page=doc.reload_page(page)
        fixed_text=page.get_text(clip=fitz.INFINITE_RECT())
        if replacement_count:
            corrections.append({'from':'Unknonw','to':'Unknown','occurrences':replacement_count,'reason':'Spelling only; cell identities unchanged'})
        if 'Unknonw' in fixed_text:raise ValueError('Uncorrected misspelled Unknown label remains')
        if family=='Overview' and f.stem=='Weight':
            if title_count!=1 or 'Interaction strength' not in fixed_text or 'Number of interactions' in fixed_text:raise ValueError('Aggregate weight title correction failed')
            corrections.append({'from':'Number of interactions','to':'Interaction strength','occurrences':title_count,'reason':'Weight.pdf is generated from cellchat@net$weight'})
        # Source fonts, curves and numeric data remain native; no crop is used.
        # Recover text outside a too-small source media box by expanding the
        # publishing page before rendering, preserving all original content.
        old_rect=list(page.rect)
        if outside:
            rect=fitz.Rect(page.mediabox)
            for x in outside:rect.include_rect(fitz.Rect(x['bbox']))
            rect+=(-8,-8,8,8)
            page.set_mediabox(rect)
            page=doc.reload_page(page)
        scale=MAX_SIDE/max(page.rect.width,page.rect.height)
        pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False)
        data=pix.tobytes('png')
        output.parent.mkdir(parents=True,exist_ok=True)
        temp=output.with_suffix('.png.tmp');temp.write_bytes(data);os.replace(temp,output)
        row.update(status='prepared',sha256=digest(data),bytes=len(data),image_width=pix.width,
                   image_height=pix.height,source_page_rect=old_rect,published_page_rect=list(page.rect),
                   text_corrections=corrections,source_text_outside_page=outside,
                   source_text_characters=len(raw_text),rendered_text_characters=len(fixed_text))
        doc.close()
    except Exception as e:
        row.update(status='failed',error=f'{type(e).__name__}: {e}')
    return row


def contacts(rows):
    font=qa_font(16)
    subfont=qa_font(12)
    out=HERE/'visual_qc';out.mkdir(exist_ok=True)
    object_ids=list(dict.fromkeys(r['dataset_object'] for r in rows))
    # Every object is represented by three overview views, plus two sheets of
    # all nine families from the smallest and largest output directories.
    groups=[]
    for stem in ('Count','Weight','cellRole'):
        groups.append((f'overview_{stem}',[r for r in rows if r['class']=='Overview' and Path(r['name']).stem==stem]))
    count=collections.Counter(r['dataset_object'] for r in rows)
    for label,obj in [('smallest',min(count,key=count.get)),('largest',max(count,key=count.get))]:
        pick=[]
        for c in CLASSES:
            rs=[r for r in rows if r['dataset_object']==obj and r['class']==c]
            if rs:pick.append(rs[len(rs)//2])
        groups.append((label+'_all_families',pick))
    for stem,rs in groups:
        cols=4 if len(rs)>12 else 3;w=420;h=450
        sheet=Image.new('RGB',(cols*w,((len(rs)+cols-1)//cols)*h),'#edf1f4');draw=ImageDraw.Draw(sheet)
        for i,r in enumerate(rs):
            x=(i%cols)*w;y=(i//cols)*h
            with Image.open(r['path']) as img:
                img.thumbnail((w-20,h-76))
                sheet.paste(img,(x+(w-img.width)//2,y+10+(h-76-img.height)//2))
            obj=r['dataset_object'].split('/')
            draw.text((x+8,y+h-62),obj[0][:43],font=subfont,fill='black')
            draw.text((x+8,y+h-45),'/'.join(obj[1:]),font=subfont,fill='black')
            draw.text((x+8,y+h-26),f"{r['class']} / {r['name'][:42]}",font=subfont,fill='black')
        sheet.save(out/(stem+'.jpg'),quality=90)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4);parser.add_argument('--limit',type=int);args=parser.parse_args()
    plan=json.loads(PLAN.read_text());tasks=[]
    for obj in plan['objects']:
        oid=obj['object_id'];directory=SOURCE/oid
        for f in sorted(directory.rglob('*.pdf')):
            rel=f.relative_to(directory)
            family='Overview' if len(rel.parts)==1 else rel.parts[0]
            if family not in CLASSES:raise RuntimeError(f'Unknown CellChat family: {rel}')
            if family=='Overview' and f.stem not in ('Count','Weight','cellRole'):raise RuntimeError(f'Unknown overview PDF: {rel}')
            tasks.append((oid,str(f),family))
    if args.limit:tasks=tasks[:args.limit]
    start=time.time();rows=[];partial=HERE/'asset_index.partial.jsonl'
    with partial.open('w') as log,concurrent.futures.ProcessPoolExecutor(max_workers=min(4,args.workers)) as pool:
        for row in pool.map(prepare,tasks,chunksize=4):
            rows.append(row);log.write(json.dumps(row,ensure_ascii=False)+'\n');log.flush()
            if len(rows)%500==0:print(json.dumps({'processed':len(rows),'total':len(tasks),'seconds':round(time.time()-start),'failed':sum(r['status']=='failed' for r in rows)}),flush=True)
    prepared=[r for r in rows if r['status']=='prepared'];failed=[r for r in rows if r['status']!='prepared']
    final=HERE/'asset_index.jsonl'
    final.write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in prepared))
    objects=[]
    for obj in plan['objects']:
        rs=[r for r in rows if r['dataset_object']==obj['object_id']]
        objects.append({'dataset_object':obj['object_id'],'class_counts':dict(collections.Counter(r['class'] for r in rs if r['status']=='prepared')),
          'source_count':len(rs),'prepared_count':sum(r['status']=='prepared' for r in rs),
          'failed_count':sum(r['status']!='prepared' for r in rs),'bytes':sum(r.get('bytes',0) for r in rs),
          'all_nine_classes_present':set(r['class'] for r in rs if r['status']=='prepared')==set(CLASSES)})
    manifest={'schema_version':1,'source_root':str(SOURCE),'asset_root':str(HERE/'assets'),'index':str(final),
      'object_count':len(objects),'source_pdf_count':len(rows),'prepared_count':len(prepared),'failed_count':len(failed),
      'class_counts':dict(collections.Counter(r['class'] for r in prepared)),
      'published_bytes':sum(r['bytes'] for r in prepared),'text_correction_files':sum(bool(r['text_corrections']) for r in prepared),
      'unknown_spelling_files':sum(any(c['from']=='Unknonw' for c in r['text_corrections']) for r in prepared),
      'weight_title_files':sum(any(c['to']=='Interaction strength' for c in r['text_corrections']) for r in prepared),
      'source_outside_text_files':sum(bool(r['source_text_outside_page']) for r in prepared),
      'worker_count':min(4,args.workers),'longest_side_pixels':MAX_SIDE,'elapsed_seconds':round(time.time()-start,2),
      'scientific_values_recomputed':False,'source_files_modified':False,'objects':objects,'failed':failed}
    (HERE/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    (HERE/'display_label_corrections.jsonl').write_text(''.join(json.dumps({k:r[k] for k in ['dataset_object','class','source_pdf','path','source_sha256','sha256','text_corrections']},ensure_ascii=False)+'\n' for r in prepared if r['text_corrections']))
    contacts(prepared)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ('objects','failed')},ensure_ascii=False),flush=True)
    if failed:raise SystemExit(1)
if __name__=='__main__':main()

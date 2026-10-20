#!/usr/bin/env python3
"""Reposition only text and label leaders in two installed CellChat scatter PNGs.
Original PDF/RDS files, marker shapes, plot axes and all numeric values remain
unchanged. This is deterministic scientific-figure rendering from PDF sources.
"""
from pathlib import Path
import hashlib,json,os,re,datetime
import pymupdf as fitz
from PIL import Image
from runtime_paths import WORK, PLAN_PATH, FIGURE_ROOT, INPUT_ROOT, configured
HERE=WORK/'cellchat'
INDEX=HERE/'asset_index.runtime.jsonl'
SPECS={
 'Atopic_Dermatitis/GSE230575/SkinBiopsy/SignalingCellRole/MHC-II.pdf':{
  'Melanocytes':{'origin':[86.0,350.0],'line':[[111.0,354.0],[57.2,380.94]],'color':[.651,.337,.157]}},
 'Crohn_Disease/GSE284230/SmallIntestine/SignalingCellRole/CD99.pdf':{
  'Smooth muscle cells':{'origin':[174.0,134.0],'line':[[258.0,130.0],[283.86,125.58]],'color':[.106,.620,.467]},
  'Pericytes':{'origin':[212.0,150.0],'line':[[250.5,145.5],[267.78,142.7]],'color':[.329,.690,.894]},
  'Glial cells':{'origin':[230.0,115.0],'line':[[270.5,112.3],[289.15,118.0]],'color':[.949,.580,.012]},
  'NK cells':{'origin':[314.0,110.0],'line':[[312.0,107.5],[308.6,109.0]],'color':[.651,.337,.157]}},
}
BT_RE=re.compile(rb'BT\b.*?\bET',re.S)
TM_RE=re.compile(rb'([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+([-+0-9.]+)\s+Tm')

def sha(data):return hashlib.sha256(data).hexdigest()

def main():
 rows=[json.loads(s) for s in INDEX.read_text().splitlines()];records=[]
 for relative,moves in SPECS.items():
  selected=[r for r in rows if r['source_pdf'].endswith('/'+relative)]
  if len(selected)!=1:raise ValueError('Expected one exact installed source mapping')
  row=selected[0];source=Path(row['source_pdf']);target=Path(row['path'])
  if not target.is_file():raise ValueError('Installed PNG is missing: '+str(target))
  if row.get('post_install_label_adjustment'):raise ValueError('Already adjusted; refusing a second mutation')
  before_png=target.read_bytes()
  if sha(before_png)!=row['sha256']:raise ValueError('Installed PNG hash differs from runtime index')
  source_bytes=source.read_bytes()
  if sha(source_bytes)!=row['source_sha256']:raise ValueError('Original PDF hash changed')
  doc=fitz.open(stream=source_bytes,filetype='pdf');page=doc[0]
  before_text=page.get_text(clip=fitz.INFINITE_RECT());before_drawings=page.get_drawings()
  original_numbers=re.findall(r'(?<![A-Za-z])[0-9]+(?:\.[0-9]+)?',before_text)
  positions={};changed={k:0 for k in moves}
  for x in page.get_contents():
   original_stream=doc.xref_stream(x)
   # Apply earlier approved spelling correction if one occurs in the source.
   stream=original_stream.replace(b'(Unknonw)',b'(Unknown)')
   def move_block(match):
    block=match.group();label=b''.join(re.findall(rb'\(([^()]*)\)',block)).decode('latin1')
    if label not in moves:return block
    matrix=TM_RE.search(block)
    if matrix is None:raise ValueError('Target label has no text matrix')
    old=list(map(float,matrix.groups()));new=moves[label]['origin']
    positions[label]={'from':[old[4],float(page.rect.height)-old[5]],'to':new}
    replacement=b' '.join(matrix.groups()[:4])+f' {new[0]:.2f} {page.rect.height-new[1]:.2f} Tm'.encode()
    changed[label]+=1
    return block[:matrix.start()]+replacement+block[matrix.end():]
   updated=BT_RE.sub(move_block,stream)
   if BT_RE.sub(b'',updated)!=BT_RE.sub(b'',original_stream):raise ValueError('Original non-text drawing operators changed')
   doc.update_stream(x,updated)
  if any(count!=1 for count in changed.values()):raise ValueError('Target labels did not map exactly once')
  page=doc.reload_page(page)
  if page.get_drawings()!=before_drawings:raise ValueError('Original vector geometry changed after text repositioning')
  for label,spec in moves.items():
   page.draw_line(fitz.Point(*spec['line'][0]),fitz.Point(*spec['line'][1]),color=spec['color'],width=.43,stroke_opacity=.5,overlay=True)
  page=doc.reload_page(page)
  after_text=page.get_text(clip=fitz.INFINITE_RECT());after_drawings=page.get_drawings()
  if after_drawings[:len(before_drawings)]!=before_drawings:raise ValueError('Original drawn data changed after adding label leaders')
  if after_text!=before_text.replace('Unknonw','Unknown'):raise ValueError('Figure text content changed')
  if re.findall(r'(?<![A-Za-z])[0-9]+(?:\.[0-9]+)?',after_text)!=original_numbers:raise ValueError('Numeric axis or legend labels changed')
  label_boxes=[]
  for block in page.get_text('dict')['blocks']:
   for line in block.get('lines',[]):
    for span in line['spans']:
     if span['text'].strip() in moves:label_boxes.append((span['text'].strip(),fitz.Rect(span['bbox'])))
  all_spans=[span for block in page.get_text('dict')['blocks'] for line in block.get('lines',[]) for span in line['spans']]
  collisions=[]
  for label,box in label_boxes:
   for span in all_spans:
    if span['text'].strip()==label:continue
    other=fitz.Rect(span['bbox']);intersection=box&other
    if not intersection.is_empty and intersection.get_area()/min(box.get_area(),other.get_area())>=.15:collisions.append([label,span['text']])
  if collisions:raise ValueError('Adjusted label still overlaps text: '+str(collisions))
  pix=page.get_pixmap(matrix=fitz.Matrix(1600/max(page.rect.width,page.rect.height),1600/max(page.rect.width,page.rect.height)),alpha=False)
  rendered=pix.tobytes('png')
  backup=HERE/'installed_label_polish_backup'/row['relative_path']/row['name'];backup.parent.mkdir(parents=True,exist_ok=True)
  if backup.exists():raise ValueError('Original installed backup already exists')
  backup.write_bytes(before_png)
  temporary=target.with_suffix('.png.label-polish.tmp');temporary.write_bytes(rendered);os.replace(temporary,target)
  with Image.open(target) as image:
   assert list(image.size)==[row['image_width'],row['image_height']];image.verify()
  adjustment={'original_published_sha256':row['sha256'],'final_sha256':sha(rendered),'original_bytes':row['bytes'],'final_bytes':len(rendered),
   'backup_path':str(backup),'text_origins':positions,'added_label_leader_count':len(moves),'numeric_text_identical':True,
   'all_original_vector_drawing_geometry_identical':True,'source_pdf_unchanged':sha(source.read_bytes())==row['source_sha256'],'all_adjusted_labels_no_significant_text_overlap':not collisions}
  row.update(sha256=sha(rendered),bytes=len(rendered),post_install_label_adjustment=adjustment)
  records.append({'dataset_object':row['dataset_object'],'class':row['class'],'name':row['name'],'path':str(target),'source_pdf':str(source),'adjustment':adjustment})
  doc.close()
 temp=INDEX.with_suffix('.jsonl.tmp');temp.write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in rows));os.replace(temp,INDEX)
 manifest_path=HERE/'manifest.json';manifest=json.loads(manifest_path.read_text());manifest['published_bytes']=sum(r['bytes'] for r in rows)
 for obj in manifest['objects']:obj['bytes']=sum(r['bytes'] for r in rows if r['dataset_object']==obj['dataset_object'])
 manifest['post_install_label_adjustments']=records
 manifest['stage_index_note']='asset_index.jsonl preserves pre-install staging hashes; asset_index.runtime.jsonl is authoritative after the two installed label-layout adjustments.'
 manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
 receipt={'passed':True,'phase':'two_installed_cellrole_label_adjustments','adjusted_image_count':len(records),'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),
  'source_pdf_files_modified':False,'original_vector_drawings_preserved':True,'all_numeric_text_and_data_coordinates_preserved':True,'checks':'Exact original vector drawing commands and paths before/after equal; only text Tm translation plus appended label leader lines changed. PNG CRC, dimensions and SHA validated.','records':records}
 (HERE/'post_install_label_polish_QA.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2))
 print(json.dumps(receipt,ensure_ascii=False))
if __name__=='__main__':main()

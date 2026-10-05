"""Explicit read-only marker-audit input identity and path resolution."""
from pathlib import Path
import pandas as pd

def load_object_crosswalk(path):
 """Require explicit, nonempty, unique original-object mappings."""
 table=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
 columns=['Data_object','Corrected_Data_object']
 if any(c not in table for c in columns):
  raise ValueError('Crosswalk requires Data_object and Corrected_Data_object')
 for column in columns:
  if table[column].str.strip().eq('').any():
   raise ValueError(f'Crosswalk has empty {column}')
 if table.Data_object.duplicated().any():
  raise ValueError('Crosswalk original Data_object is duplicated or conflicting')
 return dict(zip(table.Data_object,table.Corrected_Data_object))

def exact_left_join(left,right,on,validate,context):
 """Match every original row; an incomplete join must never redefine scope."""
 for label,table in [('left',left),('right',right)]:
  for key in on:
   if key not in table:raise ValueError(f'{context}: {label} lacks {key}')
   if table[key].isna().any() or table[key].astype(str).str.strip().eq('').any():
    raise ValueError(f'{context}: {label} has empty identity field {key}')
 marker='_marker_exact_join_status'
 if marker in left or marker in right:raise ValueError('Reserved join-status column exists')
 joined=left.merge(right,on=on,how='left',validate=validate,indicator=marker,sort=False)
 unmatched=int(joined[marker].ne('both').sum())
 if len(joined)!=len(left) or unmatched:
  raise ValueError(f'{context}: incomplete exact join ({unmatched} unmatched rows); refusing to drop input cells/groups')
 return joined.drop(columns=[marker])

def resolve_h5_path(value,data_h5_root=None):
 path=Path(value)
 if path.is_absolute():result=path.resolve()
 else:
  if '..'in path.parts:raise ValueError('Relative H5 path must not contain ..')
  if not path.parts:raise ValueError('Empty relative H5 path')
  if data_h5_root is None:raise ValueError('Relative registry h5_paths require --data-h5-root')
  root=Path(data_h5_root).resolve();result=(root/path).resolve()
  if not result.is_relative_to(root):raise ValueError('Relative H5 path escapes --data-h5-root')
 if not result.is_file():raise FileNotFoundError(f'Registered H5 file does not exist: {result}')
 return result

def prepare_specimen_source(source,crosswalk):
 source=source.copy();source.Data_object=source.Data_object.replace(crosswalk)
 if'Sample'in source:
  if'Local_sample_label'in source and not source.Sample.astype(str).equals(source.Local_sample_label.astype(str)):
   raise ValueError('Sample and Local_sample_label disagree; provide an explicit reconciled specimen mapping')
 elif'Local_sample_label'in source:source=source.rename(columns={'Local_sample_label':'Sample'})
 else:raise ValueError('Specimen source table requires Sample or Local_sample_label')
 required=['Data_object','Sample','Dataset','Condition','Tissue']
 missing=[k for k in required if k not in source]
 if missing:raise ValueError(f'Specimen source table lacks required fields: {missing}')
 for col,value in [('Source_donor_IDs',''),('Donor_status','unresolved: source donor metadata not supplied'),('Donor_analysis_eligibility','not established from provided specimen table')]:
  if col not in source:source[col]=value
 cols=required+['Source_donor_IDs','Donor_status','Donor_analysis_eligibility']
 if source[['Data_object','Sample']].duplicated().any():raise ValueError('Source specimen mapping is not unique')
 return source[cols]

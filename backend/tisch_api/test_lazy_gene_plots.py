"""Original-value/axis conversion, case distinction and strict lazy-cache scope."""
import hashlib,importlib.util,json,os,threading,unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
from unittest.mock import patch
import h5py
import numpy as np
from PIL import Image

from . import lazy_gene_plots as plots

from . import gene_index_builder as PREP


class LazyGeneTests(unittest.TestCase):
 def setUp(self):
  self.tmp=TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name);self.fig=self.root/'Fig';self.indices=self.root/'indices';self.indices.mkdir()
  self.object='Test/GSE1/PBMC';self.genes=['CD3D','C2orf15','C2ORF15','ZERO','THRA1/BTR','THRA1%2FBTR'];self.dense=np.array([[0.,1.,0.,0.,4.3,0.],[2.5,0.,3.2,0.,0.,6.7],[0.,4.1,0.,0.,6.2,0.],[6.,0.,7.,0.,0.,2.1],[0.,0.,0.,0.,0.,0.]])
  self.coords=np.array([[-2.,1.],[1.,0.],[.5,3.],[4.,2.],[1.1,-2.2]])
  source=self.root/'RNA.h5ad';columns=[];values=[];ptr=[0]
  for row in self.dense:
   nonzero=np.flatnonzero(row);columns.extend(nonzero);values.extend(row[nonzero]);ptr.append(len(values))
  with h5py.File(source,'w') as f:
   g=f.create_group('X');g.attrs.update({'encoding-type':'csr_matrix','shape':self.dense.shape});g.create_dataset('data',data=np.array(values));g.create_dataset('indices',data=np.array(columns,dtype=np.int32));g.create_dataset('indptr',data=np.array(ptr,dtype=np.int32))
   for key,axis in [('obs',[f'cell-{i}' for i in range(5)]),('var',self.genes)]:
    g=f.create_group(key);g.attrs['_index']='_index';g.create_dataset('_index',data=np.array(axis,dtype=object),dtype=h5py.string_dtype())
   f.create_group('obsm').create_dataset('X_umap',data=self.coords);f.create_group('uns').create_group('scaid').create_dataset('object_id',data=np.array([self.object],dtype=object),dtype=h5py.string_dtype())
  self.entry=PREP.build_index(source,self.indices/'gene.h5',self.object,hashlib.sha256(source.read_bytes()).hexdigest(),max_nnz_chunk=3)
  self.manifest=self.root/'manifest.json';self.value={'schema_version':1,'figure_root':str(self.fig),'index_root':str(self.indices),'objects':{self.object:self.entry}};self.write_manifest()

 def write_manifest(self):
  plots._load_manifest.cache_clear();self.manifest.write_text(json.dumps(self.value))

 def record(self,gene='CD3D',obj=None):
  folder=(obj or self.object)+'/GeneUmap';name='Umap_'+gene+'.png';target=self.fig/folder/plots._gene_filename(gene)
  return SimpleNamespace(name=name,relative_path=folder,path=str(target),image=str(target),page_number=1)

 def test_all_original_gene_values_umap_cells_and_case_distinction_preserved(self):
  entry={**self.entry,'object_id':self.object}
  for i,gene in enumerate(self.genes):
   got,coords,values=plots.read_gene(entry,i);self.assertEqual(got,gene);np.testing.assert_array_equal(coords,self.coords);np.testing.assert_array_equal(values,self.dense[:,i])
  self.assertEqual(self.entry['source_case_collisions'],[['C2orf15','C2ORF15']]);self.assertTrue(self.entry['full_source_component_readback']);self.assertTrue(self.entry['final_compressed_component_readback'])

 def test_only_registered_object_gene_and_exact_database_paths_may_render(self):
  self.assertIsNone(plots.resolve_registered_gene(self.record('UNKNOWN'),self.manifest));self.assertIsNone(plots.resolve_registered_gene(self.record(obj='Other/GSE1/PBMC'),self.manifest))
  self.assertIsNone(plots.resolve_registered_gene(self.record('../../etc/passwd'),self.manifest))
  record=self.record();record.path=str(self.fig/'Other/GSE1/PBMC/GeneUmap/Umap_CD3D.png')
  with self.assertRaises(plots.GeneAssetError):plots.resolve_registered_gene(record,self.manifest)
  record=self.record();record.image='/etc/passwd'
  with self.assertRaises(plots.GeneAssetError):plots.resolve_registered_gene(record,self.manifest)

 def test_manifest_outside_index_root_and_symlink_target_are_rejected(self):
  self.value['objects'][self.object]['index_path']='/etc/passwd';self.write_manifest()
  with self.assertRaises(plots.GeneAssetError):plots.resolve_registered_gene(self.record(),self.manifest)
  self.value['objects'][self.object]=self.entry;self.write_manifest();target=Path(self.record().path);target.parent.mkdir(parents=True);target.symlink_to(self.root/'elsewhere.png')
  with self.assertRaises(plots.GeneAssetError):plots.resolve_registered_gene(self.record(),self.manifest)

 def test_changed_index_stat_and_shape_fail_closed(self):
  entry={**self.entry,'object_id':self.object,'index_bytes':self.entry['index_bytes']+1}
  with self.assertRaises(plots.GeneAssetError):plots.read_gene(entry,0)
  entry={**self.entry,'object_id':self.object,'cell_count':6}
  with self.assertRaises(plots.GeneAssetError):plots.read_gene(entry,0)

 def test_render_is_complete_exact_size_cache_and_zero_feature_work(self):
  path=plots.ensure_registered_gene_figure(self.record(),self.manifest)
  with Image.open(path) as image:self.assertEqual(image.size,(1650,1500));image.verify()
  original=path.read_bytes()
  with patch.object(plots,'_render',side_effect=AssertionError('Cache hit must not recompute')):self.assertEqual(plots.ensure_registered_gene_figure(self.record(),self.manifest),path)
  self.assertEqual(path.read_bytes(),original)
  zero=plots.ensure_registered_gene_figure(self.record('ZERO'),self.manifest)
  with Image.open(zero) as image:self.assertEqual(image.size,(1650,1500));image.verify()

 def test_competing_requests_generate_one_gene_once(self):
  record=self.record();calls=[];errors=[]
  def render(entry,column,target):calls.append(column);target.write_bytes(b'generated-once')
  def request():
   try:plots.ensure_registered_gene_figure(record,self.manifest)
   except Exception as error:errors.append(str(error))
  with patch.object(plots,'_render',side_effect=render):
   threads=[threading.Thread(target=request) for _ in range(3)]
   for thread in threads:thread.start()
   for thread in threads:thread.join()
  self.assertFalse(errors);self.assertEqual(calls,[0])

 def test_cross_process_lock_generates_one_png_and_registry_errors_return_none(self):
  import multiprocessing,time
  context=multiprocessing.get_context('fork');record=self.record();counter=self.root/'render-calls.txt';results=context.Queue()
  def render(entry,column,target):
   with counter.open('a') as handle:handle.write(str(os.getpid())+'\n')
   time.sleep(.05);target.write_bytes(b'generated-once')
  def worker():
   try:results.put(str(plots.ensure_registered_gene_figure(record,self.manifest)))
   except Exception as error:results.put('ERROR:'+str(error))
  with patch.object(plots,'_render',side_effect=render):
   processes=[context.Process(target=worker) for _ in range(3)]
   for process in processes:process.start()
   for process in processes:process.join(5);self.assertEqual(process.exitcode,0)
  self.assertEqual(len(counter.read_text().splitlines()),1);self.assertEqual([results.get(timeout=1) for _ in range(3)],[record.path]*3)
  record=self.record();record.path='/etc/passwd'
  self.assertIsNone(plots.ensure_registered_gene_figure(record,self.manifest))

 def test_builder_rejects_an_unverified_source_file(self):
  with self.assertRaises(ValueError):PREP.build_index(self.root/'RNA.h5ad',self.indices/'bad.h5',self.object,'0'*64)
  self.assertFalse((self.indices/'bad.h5').exists())

 def test_cached_figure_fails_closed_when_its_index_is_changed(self):
  path=plots.ensure_registered_gene_figure(self.record(),self.manifest);original=path.read_bytes()
  index=Path(self.entry['index_path']);index.chmod(0o644)
  with index.open('ab') as handle:handle.write(b'changed-index')
  self.assertIsNone(plots.ensure_registered_gene_figure(self.record(),self.manifest))
  self.assertEqual(path.read_bytes(),original)

 def test_new_validated_source_version_rebuilds_old_cached_png_and_sidecar(self):
  record=self.record();path=plots.ensure_registered_gene_figure(record,self.manifest);old_hash=hashlib.sha256(path.read_bytes()).hexdigest()
  source=self.root/'RNA.h5ad'
  with h5py.File(source,'r+') as handle:
   position=int(np.flatnonzero(handle['X/indices'][:]==0)[0]);handle['X/data'][position]=22.5
  replacement=PREP.build_index(source,self.indices/'gene-v2.h5',self.object,hashlib.sha256(source.read_bytes()).hexdigest(),max_nnz_chunk=3)
  self.value['objects'][self.object]=replacement;self.write_manifest()
  original_render=plots._render
  with patch.object(plots,'_render',wraps=original_render) as render:
   self.assertEqual(plots.ensure_registered_gene_figure(record,self.manifest),path);self.assertEqual(render.call_count,1)
  self.assertNotEqual(hashlib.sha256(path.read_bytes()).hexdigest(),old_hash)
  sidecar=json.loads(plots._sidecar(path).read_text());self.assertEqual(sidecar['binding']['source_sha256'],replacement['source_sha256']);self.assertEqual(sidecar['binding']['index_sha256'],replacement['index_sha256'])

 def test_source_gene_with_slash_uses_safe_filename_without_renaming_feature(self):
  record=self.record('THRA1/BTR');self.assertEqual(record.name,'Umap_THRA1/BTR.png')
  path=plots.ensure_registered_gene_figure(record,self.manifest)
  self.assertIsNotNone(path);self.assertEqual(path.name,'Umap_THRA1%2FBTR.png');self.assertEqual(path.parent,self.fig/self.object/'GeneUmap')
  self.assertFalse((path.parent/'Umap_THRA1').exists())
  literal=plots.ensure_registered_gene_figure(self.record('THRA1%2FBTR'),self.manifest)
  self.assertIsNotNone(literal);self.assertEqual(literal.name,'Umap_THRA1%252FBTR.png');self.assertNotEqual(path,literal)
  self.assertNotEqual(plots._gene_filename('THRA1/BTR'),plots._gene_filename('THRA1%2FBTR'))
  long_feature='é'*100;self.assertTrue(plots._safe_name(long_feature));self.assertLessEqual(len(plots._gene_filename(long_feature).encode()),255);self.assertTrue(plots._gene_filename(long_feature).startswith('Umap_feature_'))
  gene,_,values=plots.read_gene({**self.entry,'object_id':self.object},4);self.assertEqual(gene,'THRA1/BTR');np.testing.assert_array_equal(values,self.dense[:,4])

if __name__=='__main__':unittest.main()

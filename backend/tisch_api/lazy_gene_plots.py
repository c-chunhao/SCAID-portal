"""Exact-object gene expression figures from validated column-oriented RNA data.

Only pre-registered PDFImage rows listed by the server's gene manifest can render
an asset. No request string is used as a source path. Original UMAP coordinates
and RNA data values are displayed without fitting, filtering or relabeling.
"""
from contextlib import contextmanager
from functools import lru_cache
import hashlib,json,logging,os,threading
from pathlib import Path
from urllib.parse import quote

LOGGER=logging.getLogger(__name__)
_RENDER_LOCK=threading.Lock()
WIDTH,HEIGHT,DPI=1650,1500,150


class GeneAssetError(RuntimeError):
    """Configured index or requested registered asset failed its contract."""


def _configured_manifest():
    try:
        from django.conf import settings
        configured=getattr(settings,'SCAID_GENE_CACHE_MANIFEST',os.environ.get('SCAID_GENE_CACHE_MANIFEST',str(Path(__file__).with_name('data')/'lazy_gene_plots.json')))
    except (ImportError,RuntimeError):
        configured=os.environ.get('SCAID_GENE_CACHE_MANIFEST',str(Path(__file__).with_name('data')/'lazy_gene_plots.json'))
    return Path(configured) if configured else None


def _gene_filename(gene):
    filename='Umap_'+quote(gene,safe='-_.~')+'.png'
    return filename if len(filename.encode('utf-8'))<=255 else 'Umap_feature_'+hashlib.sha256(gene.encode('utf-8')).hexdigest()+'.png'

def _safe_name(gene):
    return bool(gene) and gene==gene.strip() and gene.isprintable() and len('Umap_'+gene+'.png')<=255 and len(_gene_filename(gene).encode('utf-8'))<=255


@lru_cache(maxsize=2)
def _load_manifest(filename,mtime_ns,size):
    value=json.loads(Path(filename).read_text())
    if value.get('schema_version')!=1 or not isinstance(value.get('objects'),dict):
        raise GeneAssetError('Invalid gene figure manifest')
    figure_root=Path(value['figure_root']).resolve()
    index_root=Path(value['index_root']).resolve()
    folders={}
    for object_id,entry in value['objects'].items():
        folder=entry['relative_path']
        if folder!=object_id+'/GeneUmap' or any(p in {'.','..'} for p in Path(folder).parts) or Path(folder).is_absolute():
            raise GeneAssetError('Invalid gene figure folder')
        path=Path(entry['index_path'])
        if not path.is_absolute() or not path.resolve().is_relative_to(index_root):
            raise GeneAssetError('Index lies outside the configured root')
        genes=entry['genes']
        if len(genes)!=entry['feature_count'] or len(set(genes))!=len(genes) or not all(_safe_name(g) for g in genes):
            raise GeneAssetError('Invalid or ambiguous gene axis')
        filenames=entry.get('filenames',[_gene_filename(g) for g in genes])
        if filenames!=[_gene_filename(g) for g in genes] or len(set(filenames))!=len(genes):raise GeneAssetError('Invalid or ambiguous physical gene filenames')
        folders[folder]={**entry,'filenames':filenames,'object_id':object_id,'by_name':{'Umap_'+g+'.png':i for i,g in enumerate(genes)}}
    return figure_root,index_root,folders


def _row_path(record):
    # Both stored database references must agree with the registered filename.
    return str(getattr(record,'path','')),str(getattr(record,'image',''))


def resolve_registered_gene(record,manifest_path=None):
    """Return a strict registry binding, or None for an unrelated figure row."""
    manifest=Path(manifest_path) if manifest_path is not None else _configured_manifest()
    if manifest is None or not manifest.is_file():return None
    try:
        stat=manifest.stat();root,index_root,folders=_load_manifest(str(manifest.resolve()),stat.st_mtime_ns,stat.st_size)
        folder=str(getattr(record,'relative_path',''));entry=folders.get(folder)
        if entry is None:return None
        name=str(getattr(record,'name',''));column=entry['by_name'].get(name)
        if column is None:return None
        if int(getattr(record,'page_number',0))!=1:raise GeneAssetError('Gene figure must have one page')
        target=root/folder/entry['filenames'][column]
        if _row_path(record)!=(str(target),str(target)):
            raise GeneAssetError('Figure paths do not match the exact registered object and gene')
        if target.is_symlink() or target.parent.resolve()!=root/folder or not target.resolve().is_relative_to(root):
            raise GeneAssetError('Gene cache target contains a symlink or escapes its root')
        return entry,column,target
    except (OSError,ValueError,KeyError,TypeError) as exc:
        raise GeneAssetError('Gene figure registry validation failed') from exc


@contextmanager
def _file_lock(path):
    import fcntl
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.is_symlink() or path.parent.is_symlink():raise GeneAssetError('Unsafe gene lock target')
    descriptor=os.open(path,os.O_CREAT|os.O_RDWR|os.O_NOFOLLOW,0o600)
    try:
        with os.fdopen(descriptor,'a+b') as stream:
            fcntl.flock(stream.fileno(),fcntl.LOCK_EX)
            try:yield
            finally:fcntl.flock(stream.fileno(),fcntl.LOCK_UN)
    except Exception:
        raise


def _validate_index(entry,column):
    import h5py
    index=Path(entry['index_path'])
    if index.is_symlink() or not index.is_file():raise GeneAssetError('Gene index unavailable')
    stat=index.stat()
    if stat.st_size!=entry['index_bytes'] or stat.st_mtime_ns!=entry['index_mtime_ns']:
        raise GeneAssetError('Gene index changed since validation')
    with h5py.File(index,'r') as handle:
        if (int(handle.attrs.get('schema_version',0))!=1 or handle.attrs.get('object_id')!=entry['object_id']
                or handle.attrs.get('source_sha256')!=entry['source_sha256']):
            raise GeneAssetError('Gene index provenance mismatch')
        if (handle.attrs.get('gene_axis_sha256')!=entry['gene_axis_sha256'] or handle.attrs.get('cell_axis_sha256')!=entry['cell_axis_sha256']
                or handle['genes'].shape!=(entry['feature_count'],) or handle['cells'].shape!=(entry['cell_count'],)
                or tuple(handle['expression_csc'].attrs['shape'])!=(entry['cell_count'],entry['feature_count'])):
            raise GeneAssetError('Gene/cell index axes or shape do not match the validated registry')
        if not 0<=column<entry['feature_count']:raise GeneAssetError('Gene column outside axis')
        if handle['genes'].asstr()[column]!=entry['genes'][column]:raise GeneAssetError('Gene axis mismatch')


def read_gene(entry,column):
    """Read one exact CSC column and original UMAP, preserving all cell rows."""
    import h5py
    import numpy as np
    _validate_index(entry,column)
    with h5py.File(entry['index_path'],'r') as handle:
        start,end=map(int,handle['expression_csc/indptr'][column:column+2])
        values=np.asarray(handle['expression_csc/data'][start:end],dtype=np.float64)
        cells=np.asarray(handle['expression_csc/indices'][start:end],dtype=np.int64)
        n=entry['cell_count'];coords=np.asarray(handle['umap'],dtype=np.float64)
        if (coords.shape!=(n,2) or not np.isfinite(coords).all() or len(cells)!=len(values)
                or (len(cells) and (cells.min()<0 or cells.max()>=n or np.any(np.diff(cells)<=0)))
                or not np.isfinite(values).all() or (len(values) and values.min()<0)):
            raise GeneAssetError('Invalid gene values or cell coordinates')
        expression=np.zeros(n,dtype=np.float64);expression[cells]=values
    return entry['genes'][column],coords,expression


def _cache_binding(entry,column):
    return {key:entry[key] for key in ('object_id','source_sha256','index_sha256','index_bytes','index_mtime_ns','gene_axis_sha256','cell_axis_sha256','umap_sha256')}|{
        'gene':entry['genes'][column],'column':column,'renderer_version':1,'width':WIDTH,'height':HEIGHT}


def _sidecar(target):
    directory=target.parent/'.cache_metadata'
    if directory.is_symlink():raise GeneAssetError('Unsafe gene cache metadata directory')
    path=directory/(hashlib.sha256(target.name.encode()).hexdigest()+'.json')
    if path.is_symlink():raise GeneAssetError('Unsafe gene cache metadata file')
    return path


def _cached(target,binding):
    if not target.is_file() or target.is_symlink():return False
    sidecar=_sidecar(target)
    if not sidecar.is_file():return False
    try:
        value=json.loads(sidecar.read_text())
        if value.get('binding')!=binding or value.get('png_bytes')!=target.stat().st_size:return False
        with target.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()==value.get('png_sha256')
    except (OSError,ValueError,TypeError):return False


def _record_cache(target,binding):
    sidecar=_sidecar(target);sidecar.parent.mkdir(parents=True,exist_ok=True)
    with target.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    value={'binding':binding,'png_bytes':target.stat().st_size,'png_sha256':digest}
    temporary=sidecar.with_name('.'+sidecar.name+'.'+str(os.getpid())+'.tmp')
    try:
        temporary.write_text(json.dumps(value,ensure_ascii=False,sort_keys=True)+'\n');os.replace(temporary,sidecar)
    finally:
        if temporary.exists():temporary.unlink()


def _render(entry,column,target):
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import LinearSegmentedColormap,Normalize
    gene,umap,expression=read_gene(entry,column)
    # Matplotlib font caches are process-global; isolate each render in a worker.
    with _RENDER_LOCK:
        fig=Figure(figsize=(WIDTH/DPI,HEIGHT/DPI),dpi=DPI,facecolor='white')
        FigureCanvasAgg(fig);ax=fig.add_axes([.10,.12,.73,.77]);ax.set_facecolor('white')
        order=np.argsort(expression,kind='stable')
        maximum=float(expression.max())
        if maximum>0:
            cmap=LinearSegmentedColormap.from_list('scaid_expression',['#E1E1E1','#7663AA','#3E2672'])
            image=ax.scatter(umap[order,0],umap[order,1],c=expression[order],s=2.4,linewidths=0,
                cmap=cmap,norm=Normalize(vmin=0,vmax=maximum),rasterized=True)
            cax=fig.add_axes([.87,.24,.025,.50]);bar=fig.colorbar(image,cax=cax);bar.set_label('Normalized RNA expression',fontsize=12);bar.ax.tick_params(labelsize=11)
        else:
            ax.scatter(umap[:,0],umap[:,1],c='#D3D3D3',s=2.4,linewidths=0,rasterized=True)
            ax.text(.99,.99,'Expression = 0 in this object',ha='right',va='top',transform=ax.transAxes,fontsize=12,color='#444444')
        ax.set_title(gene,fontsize=21,pad=16,parse_math=False);ax.set_xlabel('UMAP 1',fontsize=14);ax.set_ylabel('UMAP 2',fontsize=14)
        ax.tick_params(labelsize=11);ax.spines[['top','right']].set_visible(False);ax.margins(.025);ax.set_aspect('equal',adjustable='box')
        temporary=target.with_name('.'+target.name+'.'+str(os.getpid())+'.tmp')
        try:
            fig.savefig(temporary,format='png',dpi=DPI,facecolor='white',metadata={'Title':gene,'Description':'Original UMAP coordinates and RNA data-layer expression; all cells retained.'})
            with temporary.open('rb') as stream:os.fsync(stream.fileno())
            os.replace(temporary,target)
        finally:
            fig.clear()
            if temporary.exists():temporary.unlink()


def _ensure_registered_gene_figure(record,manifest_path=None):
    """Materialize a registered cached figure; return None for any other row.

    Invalid configured indices raise GeneAssetError for the caller to fail closed.
    Existing complete cached PNGs are returned without reading expression data.
    """
    binding=resolve_registered_gene(record,manifest_path)
    if binding is None:return None
    entry,column,target=binding
    _validate_index(entry,column)
    cache_binding=_cache_binding(entry,column)
    if _cached(target,cache_binding):return target
    target.parent.mkdir(parents=True,exist_ok=True)
    lock=target.parent/'.locks'/(hashlib.sha256(target.name.encode()).hexdigest()+'.lock')
    with _file_lock(lock):
        if target.is_symlink():raise GeneAssetError('Unsafe gene cache target')
        _validate_index(entry,column)
        if not _cached(target,cache_binding):
            _render(entry,column,target)
            _record_cache(target,cache_binding)
    return target


def ensure_registered_gene_figure(record,manifest_path=None):
    """Return an approved cached PNG or None; malformed registries fail closed."""
    try:
        return _ensure_registered_gene_figure(record,manifest_path)
    except (GeneAssetError,OSError,ValueError,KeyError,TypeError,ImportError) as error:
        LOGGER.warning('Registered gene figure unavailable (%s)',type(error).__name__)
        return None

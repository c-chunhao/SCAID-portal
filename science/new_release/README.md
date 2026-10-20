# Additional module source candidate

中文索引：本目录提供格式导出、独立回读、CellChat 图片准备、marker 表校验、延迟基因图缓存和描述性 PBMC pilot 代码。输入数据与输出文件需在本机提供；既有两条 atlas 未被替换。网站新增功能已于2026年10月21日部署；验收范围见仓库根README与核验说明。

These scripts prepare the 31 additional source–condition–tissue objects. They do
not replace either existing published atlas. The original Seurat object, RNA
expression, source labels, cell order and existing reductions are retained during
format conversion. A source sample label is not an authenticated donor identity.
No input matrix, individual metadata, production registry, generated PNG cache,
maintenance credential or internal operational report is supplied here.

## Dependencies and configuration

Use a separate scientific Python environment and `pip install -r requirements.txt`.
The pins describe the observed preparation runtime; a new installation is not
claimed. R export additionally requires Seurat, SeuratObject, Matrix, rhdf5 and
jsonlite, with support for `LayerData` (Seurat v5). Set `SCAID_RSCRIPT` for a custom Rscript executable and `SCAID_QA_FONT` for
an optional review-contact-sheet font. R and Python must read the same
locally authorized source objects. Keep scientific packages separate from Django.

Run from this directory, or set the variables below to absolute local paths.
`runtime_paths.py` resolves Python defaults relative to this directory and passes
its plan/input/output paths into R subprocesses.

```sh
export SCAID_SCIENCE_WORK_ROOT=/tmp/scaid-new-release
export SCAID_SCIENCE_INPUT_ROOT=/path/to/authorized/local-inputs
export SCAID_RELEASE_PLAN="$PWD/examples/release_plan.example.json"
export SCAID_SEURAT_INPUT_ROOT="$SCAID_SCIENCE_INPUT_ROOT/00.data"
export SCAID_ANNOTATION_INPUT_ROOT="$SCAID_SCIENCE_INPUT_ROOT/01.ann"
export SCAID_CELLCHAT_INPUT_ROOT="$SCAID_SCIENCE_INPUT_ROOT/02.CellChat"
export SCAID_GSE279086_RAW_ROOT="$SCAID_SCIENCE_INPUT_ROOT/GSE279086"
export SCAID_EXPORT_ROOT="$SCAID_SCIENCE_WORK_ROOT/exports"
export SCAID_FIGURE_ROOT="$SCAID_SCIENCE_WORK_ROOT/figures"
export SCAID_GENE_INDEX_ROOT="$SCAID_SCIENCE_WORK_ROOT/gene_cache/indexes"
export SCAID_GENE_CACHE_MANIFEST="$SCAID_SCIENCE_WORK_ROOT/gene_cache/lazy_gene_plots.json"
```

`examples/release_plan.example.json` contains the 31 aggregate object identities
and source catalogue fields used by these scripts. It is not a database import.
Provide each original RDS at
`$SCAID_SEURAT_INPUT_ROOT/<condition>/<accession>/<tissue>/<condition>-<accession>-<tissue>.rds`.
Do not reinterpret the source condition/group labels as verified case/control
assignments. `SCAID_BACKEND_SOURCE` optionally selects a separate backend checkout
for the shared CSC index builder; its default is this repository's backend.

## Conversion and independent checks

```sh
Rscript export_new_aid_h5ad.R 0 "$SCAID_SCIENCE_WORK_ROOT"
python validate_new_aid_h5ad.py "$SCAID_SCIENCE_WORK_ROOT/matrices/objects/Active_Tattoo_Uveitis_GSE232802_PBMC/export_report.json"
```

The R writer exports RNA, AUCell, UCell and singscore as cells × features H5AD.
It checks every written sparse component or dense block. The independent Python
reader checks axes, source metadata, source value probes, original reductions and
integer/nonnegative counts, hashes the file and promotes it atomically only on
success. It refuses to overwrite an existing validated output. Resume validates
existing hashes rather than trusting a status flag alone.

`export_new_aid_batch.py` provides a bounded two-object queue for normal integer
RNA source layers. It may complete 30 objects and report the audited fractional
exception as failed; handle that object explicitly before finalizing the manifest.
A complete set is 31 objects × four methods = **124 main H5AD files**.

The GSE279086 Kidney object has fractional processed values in its source
`RNA/counts` layer. Preserve those values rather than rounding them:

```sh
# Index 26 belongs to this object in the supplied example plan; verify your plan.
Rscript export_new_aid_h5ad.R 26 "$SCAID_SCIENCE_WORK_ROOT" recover_raw_counts
python audit_counts_source_GSE279086.py
python declare_GSE279086_processed_only.py
python validate_new_aid_h5ad.py "$SCAID_SCIENCE_WORK_ROOT/matrices/objects/T1D_GSE279086_Kidney/export_report.json"
python finalize_new_aid_export_reports.py
```

The historical `recover_raw_counts` argument permits preservation of this one
fractional source layer; it does not itself reconstruct raw counts. After the
explicit processed-only declaration, `X` contains the original normalized data,
`layers/source_processed_counts` retains fractional source values, and
`layers/counts` is absent. Original `nCount_RNA`/`nFeature_RNA` describe processed
source values. Thus **30 RNA objects retain verified nonnegative integer source
count layers; one complete RNA object is processed-only**. Original integer
source layers are not independently re-derived from sequencing reads here.

`audit_counts_source_GSE279086.py` independently checks raw/processed gene and
barcode axes. It requires the unpromoted `.h5ad.tmp` and the original metadata
export. A full recovery must fail when any retained barcode lacks an exact
mapping; no cells are imputed or filled with invented zeros. The historical full
recovery implementation remains in `recover_GSE279086_raw_counts.py` for review.
`export_GSE279086_verified_raw_subset.py` additionally exports only the **24 fully
mapped GSMs, 23,979 cells and 36,601 original gene rows**. This supplemental file
is optional and outside the main 124-file manifest. It omits the two unmatched
GSMs and 684 retained processed cells; it does not replace the complete object.

## Prepared views and local registration

`prepare_cellchat_assets.py` renders the existing CellChat PDFs, records hashes,
retains network values, expands media boxes for clipped labels, and records only
the approved Unknown spelling and aggregate Weight-title corrections.
`finalize_cellchat_assets.py` independently checks files and prepares an index;
it does not install or register the assets. The overlap review and two explicit
label-position refinements are included with their exact-source guards.
Generated indices describe local paths and stay outside Git.

`prepare_marker_assets.py` requires a local catalogue snapshot at
`$SCAID_SCIENCE_WORK_ROOT/markers/catalog_snapshot.json`. Each row must include
`id` (or `record_id`), `abbreviation`, `dataset_id` and `tissue`. Matching uses the
three scientific fields; the output registry binds the resulting local ID too.
TSV values, including reported zero p values, are preserved. The contrast is
cell type versus the remaining cells within the same object; it is not a
validated donor-level disease comparison. Set `SCAID_MARKERS_ROOT` in Django to
the generated `markers` directory after reviewing and installing the files.
The example in `backend/docs/new_release_assets` does not automatically match a
clean database.

`prepare_additional_scores.py` checks the three remaining source image families.
`prepare_lazy_gene_index.py` builds the three exact-object RNA/UMAP CSC indices
from a validated download manifest. It retains all cells, source RNA values and
coordinates and fully checks the conversion. The API renders registered gene
PNGs on demand and reuses a hash-bound cache; it does not rerun normalization,
annotation or integration. Build/register the appropriate PDFImage rows locally;
this repository does not install the scientific library automatically.

## Descriptive PBMC pilot

```sh
python prepare_new_aid_integration.py --pilot
python review_PBMC_source_label_markers.py
```

The candidate pilot comprises **7 objects, 5 studies and 14,000 sampled PBMC
lymphoid cells**, capped at 2,000 per object. It excludes the monocyte-enriched
GSE268626 assays from the lymphoid pilot, limits verified B-cell-enriched assays
to the B/plasma lineage, and retains study-specific assay and treatment context.
The code uses shared genes, preserves counts, normalizes only the pilot analysis
copy, runs study Harmony/UMAP, and reports technical diagnostics and descriptive
marker summaries. Source labels remain source labels; the pilot does not
establish calibrated annotation accuracy, authenticated donors, whole-PBMC
composition or disease/control effects. It does not perform disease-comparison
significance tests or replace the two original atlases. Enlargement requires
reviewed annotation, assay scope and independent donor/clinical mapping.

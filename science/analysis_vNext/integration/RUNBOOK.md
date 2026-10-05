# Public module runbook

## 无私有矩阵/模型的guard检查

Python需要numpy/scipy/h5py；R需要SeuratObject/Matrix/jsonlite/digest。当前目录应为本integration目录。

```bash
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
python check_released_rna_public.py
Rscript check_integration_safeguards.R
```

Python临时toy实测fullRNA/denseHVG轴、未知dense failclosed、乱序cell-ID、full-library denominator与fraction100000.1/badptr/badindex/重复coordinates拒绝。R临时toy实测counts getter、embedding/metadata错序重排、missing/dup拒绝、稀疏全RNA normalization及从实际scale.data getter写出独立featureaxis。

## 指定自己的真实输入

```bash
python released_rna.py /path/to/released_RNA.h5 --out /path/to/new_schema.json
python reannotate_t1d_full_rna_pilot.py --help
python reannotate_t1d_full_rna_pilot.py --h5 /path/to/T1D_RNA.h5 \
  --source-metadata /path/to/T1D_original_author_metadata.tsv.gz \
  --models /path/to/cached_CellTypist_models --out /path/to/distinct_candidate_output
Rscript assemble_B_Plasma_vNext.R /path/to/B_Plasma_input_root /path/to/distinct_output.rds
```

pilot input是SCAID multi-feature H5（names_obs/names_var与RNA/rawdata CSR）；source-metadata需要cell_id/Sample_ID/Cluster_Annotation_Merged列，与H5的Sample逐项对应，仅用于审阅，模型目录需Immune_All_High.pkl和Immune_All_Low.pkl。模型以及临床cellmetadata未随代码发布。assembly预期inputroot下历史Mult_01/Mult_03(B、Plasma、Bgc)结构，base显式为Mult_01/scRNA.Mult1.rds，branch逐个只保留metadata；巨型执行前安排资源，不覆盖production。

future writer应在序列化的同一Seurat对象上调用 `export_RNA_layer_axes_checked(object,sidecar_path,source_file)`，metadata必须用match按真实matrix rownames排序。legacydense读取时必须独立验证相应source scale.data行序和export关联；fullRNA raw正常可读。

`validate_raw_reader_corruption_guards.py`另含本地真实T1D5000cells regression，需要指定本机release路径，公开环境先运行self-contained check。

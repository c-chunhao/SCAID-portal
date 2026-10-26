# RNA matrix and feature-axis specification

`names_obs` is the cell axis. `names_var` and `var/rawvar/_index` are the full raw-RNA feature axis. CSR blocks under `assay/RNA/layers/rawdata` contain full raw-RNA counts, with order matching those axes.

The existing dense matrix under `assay/RNA/layers/data` contains only 2000 HVG features. Read its column names from `var/var/_index`, not from the first 2000 genes in `names_var`. The actual matrix shape is cells × HVG, while the legacy `shape` attribute records HVG × cells. Validate reads against the real dataset shape and an explicit feature index, retaining the reversed attribute as a compatibility diagnostic. The name `data` does not establish that the matrix contains normalized RNA for all genes. A matrix with negative values cannot be treated as raw counts or CellTypist log1p input.

The new `SCAID-RNA-sidecar-1.0` JSON explicitly records each block's cell offset, actual cells × features shape, feature index, matrix semantics, encoding, reversed-attribute status and permitted analytical uses. `released_rna.py` provides a reusable reader and refuses to guess when the feature axis is uncertain.

Having 2000 axis names must also be distinguished from independently verifying the gene order in every column. The easySCFr exporter selects `rawvar[rownames(rawvar) %in% VariableFeatures(sce), ]` to create var metadata for this layer, without explicitly reordering it by `rownames(scaleData)`. The actual archived postprocessed T1D RDS (`20241116/00.data/type_1_diabetes_mellitus/syn53641849/PBMC`) was successfully read: all 2000 RNA `scale.data` row names exactly matched the H5 declared axis, element by element. The maximum difference across 16 cells × 2000 scaled values was 5.33e-15, consistent with text round-trip precision. Raw counts for the same 16 cells × 20,690 genes exactly matched the independent official source. The T1D dense axis therefore has independent support from the actual object, with no mismatch found. The link to the historical exporter call was not archived; other objects are not automatically validated by this result. The legacy dense/HVG reader rejects reads by default. T1D reading can be explicitly authorized by supplying the independently verified `T1D_export_candidate_2_scale_feature_names.txt`. Raw counts remain readable normally.

New writers must separately record `raw_feature_names = rownames(raw_counts)` and `scaled_HVG_feature_names = rownames(scaleData)`. Corresponding metadata must be reordered to that axis with `match()` and checked with `identical(rownames(metadata), rownames(matrix))`; membership filtering alone is insufficient. Each layer records whether its values are counts, log-normalized or scaled, rather than naming scaled RNA only `data`.

For CellTypist input, first normalize full raw counts to 10000 using the full RNA library, then apply log1p. If a model-gene subset is needed, select it after normalization. Do not subset HVGs first and then renormalize the HVG library to 10000. Each CSR block can be normalized and classified independently; millions of cells need not be converted to a dense matrix.

This specification addresses the risk of incorrect reads; it does not establish that gene mismatches occurred on the existing website. The current website displays pregenerated figures and does not directly query expression matrices through this dense layer. Legacy files remain unchanged. The new reader and JSON sidecar support reproducible analysis and future exports.

Reproduction:

```bash
python released_rna.py /path/to/released_RNA.h5 --out schema.json
```

The returned counts CSR is strictly checked for nonnegative integer values, row coverage and full gene order. `iter_full_log1p` retains sparsity. `read_dense_hvg` reads at most 10000 cells per call and is disabled by default; it requires `verified_source_feature_names`, the verified `scale.data` row names of the corresponding export source. Do not pass the var index read from the legacy H5 back into the reader as if it were an independent source.

Official exporter source: [easySCFr save.R](https://github.com/xleizi/easySCF/blob/main/r/R/save.R). The historical local call is in `20250706/Fig/SV/GSE198616/PBMC/revise.R:74`. The historical package version was not saved; the current official source is used only to review semantics and implementation.

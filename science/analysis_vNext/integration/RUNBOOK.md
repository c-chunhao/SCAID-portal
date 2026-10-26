# Public module runbook

## Guard checks without private matrices or models

Python requires numpy/scipy/h5py; R requires SeuratObject/Matrix/jsonlite/digest. Run from this integration directory.

```bash
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
python check_released_rna_public.py
Rscript check_integration_safeguards.R
```

Temporary Python toy tests check full-RNA/dense-HVG axes, rejection of an unverified dense axis, reordered cell IDs, the full-library denominator, and rejection of the fractional value 100000.1, invalid pointers, invalid indices and duplicate coordinates. Temporary R toy tests check the counts getter, reordering of embeddings/metadata, rejection of missing/duplicate IDs, sparse full-RNA normalization, and writing an independent feature axis from the actual `scale.data` getter.

## Specify your own real inputs

```bash
python released_rna.py /path/to/released_RNA.h5 --out /path/to/new_schema.json
python reannotate_t1d_full_rna_pilot.py --help
python reannotate_t1d_full_rna_pilot.py --h5 /path/to/T1D_RNA.h5 \
  --source-metadata /path/to/T1D_original_author_metadata.tsv.gz \
  --models /path/to/cached_CellTypist_models --out /path/to/distinct_candidate_output
Rscript assemble_B_Plasma_vNext.R /path/to/B_Plasma_input_root /path/to/distinct_output.rds
```

The pilot input is a SCAID multi-feature H5 (`names_obs`/`names_var` and RNA/rawdata CSR). Source metadata require `cell_id`/`Sample_ID`/`Cluster_Annotation_Merged` columns corresponding element by element to the H5 `Sample`, and are used only for review. The model directory requires `Immune_All_High.pkl` and `Immune_All_Low.pkl`. Models and clinical cell metadata are not distributed with the code. Assembly expects the historical `Mult_01`/`Mult_03` (B, Plasma, Bgc) structure under the input root, with the base explicitly set to `Mult_01/scRNA.Mult1.rds`; only metadata are retained from each branch. Arrange resources before large runs, and do not overwrite production outputs.

Future writers should call `export_RNA_layer_axes_checked(object,sidecar_path,source_file)` on the same Seurat object that is serialized. Metadata must be ordered with `match` against actual matrix row names. Reading legacy dense data requires independent verification of the corresponding source `scale.data` row order and its export linkage. Full raw RNA remains readable normally.

`validate_raw_reader_corruption_guards.py` also contains a local regression check using 5000 real T1D cells and requires an explicit local release path. In a public environment, run the self-contained check first.

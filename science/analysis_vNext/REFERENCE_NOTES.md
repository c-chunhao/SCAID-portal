# SCAID analysis vNext — audit reference notes, 5 October 2026

These supplementary notes record the analysis entry points, completed checks and scientific constraints from the 5 October 2026 audit. See [README.md](README.md) for the public module overview and [ENVIRONMENT_AND_REPRODUCIBILITY.md](ENVIRONMENT_AND_REPRODUCIBILITY.md) for environment and reproducibility details.

This is the new analysis entry point. Historical R scripts, original RDS/base H5 files, and earlier articles and figures are retained for provenance. Historical cell-marker intersections and pooled CellChat aggregates cannot support disease inference.

## Currently executable entry points

- `statistics/recompute_atlas_composition.py`: An explicit CLI recomputes composition from all release metadata and exact joins to the two integrations. Source-positive units with no integrated members receive NA; genuinely observed subtype zeros are retained. Inputs are read-only, and outputs go to a separate directory.
- `integration/released_rna.py`: Raw counts and full-RNA normalization are kept separate from the scaled HVG axis. Sparse blocks are read with value and index checks. Full-library normalization precedes gene selection.
- `integration/integration_safeguards.R`: Embeddings are attached by unique cell ID, metadata joins are strict, normalization uses full counts, and getters support the relevant Seurat versions. Incorrect axes or missing IDs cause rejection.
- `qc/sparse_export_guard.py`: Complete CSR-block export and independent shape checks retain singleton tail rows.
- `release_scientific_gate.py`: Typed evidence, hexadecimal SHA values, and gates for technical capture, clinical grouping, annotation and duplicates. Sample labels do not automatically establish donors.
- `new_data/review_aps_multimodal.py`: Descriptive review of paired protein/RNA, internal reference flags and sensitivity to the source 500-gene threshold. It does not automatically overwrite annotations or QC.
- `new_data/audit_new_data_sources.py`: Checks integer counts, CSR shapes/indices, complete MatrixMarket records, features and barcodes. An integer representation does not establish the provenance of source counts.

For integration, see the [README](integration/README.md) and [runbook](integration/RUNBOOK.md) in that directory. Operators must separately retain scientific evidence, input paths, versions and SHA values locally; patient data are not distributed with the source code. The other CellChat, score and statistical replay scripts replay evidence from frozen audit inputs. Public releases must separately provide environment-root configuration and path resolvers; local absolute paths are not a general installation procedure.

## Completed validation

- All 2,572,480 count rows and 3,646,001,623 nonzero values in 48 effective RNA files were checked.
- All six composition layers were independently rerun from original metadata: 563 T/NK and 537 B/plasma defined fractions; 148 unassessed lineage units were not filled with zero.
- All values, structures and identities in 48 CellChat objects and 144 score H5 files were audited.
- Independent source counts matched exactly for 16 T1D cells × 20,690 RNA genes; 16 × 2000 scaled-gene values agreed with the corresponding historical Seurat object.
- A genuine per-capture scDblFinder PSO pilot, sensitivity to QC ordering before/after doublet calling, and source barcode/count preservation were checked. Paired APS ADT CLR was independently compared with Seurat.
- The sparse exporter was checked at the 0/1/4999/5000/5001/10000/10001 boundaries and against deliberately omitted tail rows. The raw reader fails closed on fractional counts and CSR corruption.

## Required scientific constraints

1. Separate RNA from mixed-feature matrices first; antibody or HTO features must not contribute to RNA library size. Source filtered counts are not unfiltered droplets.
2. Split doublet analysis only by independently verified technical captures; retain all original barcodes/counts and QC flags.
3. A unit with a source major lineage but no integration members must receive NA. Proportions among partially retained integration members do not establish absolute abundance in the complete tissue.
4. Inference uses actual independent source individuals and requires a case/control design, clinical keys, a full-rank design and covariate overlap.
5. Mixing in a 2D UMAP and a method's own cluster labels are not independent biological validation. High/Low predictions from the same reference are not a gold standard.
6. Cell-type markers are not disease DEGs; scores are not RNA counts or p values. Probabilities from different CellChat database snapshots cannot be used to compare disease intensity.
7. The 59 new candidates are not automatically added to the total of 48 objects, and missing historical data are not fabricated. Blinded accuracy validation cannot yet be claimed for all fine labels.

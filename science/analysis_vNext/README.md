# Reviewed analysis modules

These modules cover coverage-aware composition, strict released-RNA readers,
integration safeguards, new-source and multimodal checks, score semantics,
CellChat object audit, sparse-count export checks and scientific evidence gates.
They do not distribute assays, clinical tables, cell/individual records or
downloaded annotation models. Third-party inputs remain subject to their
applicable access and reuse terms.

The current composition entry is:

```bash
python science/scripts/recompute_atlas_composition.py \
  --input-tables LOCAL_METADATA_DIRECTORY --output NEW_OUTPUT_DIRECTORY --plot
```

Supply original release and integrated metadata as documented in the script.
`--source-map` can add externally verified identity fields; it never creates
missing donor identities. Source-major-lineage-positive samples with no exact
integrated members have unavailable fractions, while assessed subtype zeros
are preserved. Denominators use annotated immune lineages plus exact integrated
members. Outputs are descriptive retained integration membership, without
assuming specimen labels are independent individuals. The frozen September 30
scripts are in `science/archive/2026-09-30/` and are deprecated.

The integration module has its own [runbook](integration/RUNBOOK.md). Its
synthetic reader tests and R safeguards are self contained. Real annotation
pilots require explicit local `--h5`, `--source-metadata`, `--models` and `--out`
inputs; source labels are used for evaluation, not prediction.

Replay audit modules use explicit local inputs:

| Modules | Required local input variables |
| --- | --- |
| `cellchat/` | `SCAID_CELLCHAT_AUDIT_DIR`; the quality gate also needs `SCAID_SCIENCE_TABLES` |
| `scores/` | `SCAID_SCORE_AUDIT_DIR`, `SCAID_ASSET_REGISTRY`; relative assets require `SCAID_DATA_H5_ROOT` |
| Frozen KEGG coverage | Also `SCAID_KEGG_SNAPSHOT_PATHS`: four archived files separated by the platform path separator |
| `statistics/review_statistics.py` | `SCAID_STATISTICS_AUDIT_DIR`, `SCAID_SCIENCE_TABLES` |

Replay audits expect their verified inventories and frozen membership tables in
the supplied local directories. They check the reviewed SCAID release scope;
they do not claim arbitrary cohorts share the same provenance. Relative assets
are resolved against the explicit data root, including the corrected Lyme RNA
revision. The registry itself is not modified.

Run checks without any real records:

```bash
python science/analysis_vNext/tests/test_portable_contracts.py
python science/analysis_vNext/tests/test_native_matrixmarket_stream.py
python science/analysis_vNext/qc/verify_sparse_export_boundaries.py
python science/analysis_vNext/integration/check_released_rna_public.py
Rscript science/analysis_vNext/integration/check_integration_safeguards.R
```

The native stream test compiles the C++ source in a temporary directory and
requires `g++` and zlib development headers. Native scanner usage is documented
in [NATIVE_SCANNER.md](new_data/NATIVE_SCANNER.md); no executable is committed.

The JSON files under `schemas/` contain only fictitious examples or column
definitions. Scientific release and donor-inference gates require actual
Boolean evidence fields and hexadecimal SHA256 values. Passing this gate
establishes evidence completeness, not biological validity or independent
validation. Missing historical QC and clinical fields are never imputed.

The [fine-label marker audit](marker_coherence/README.md) has three explicit-input
CLI stages for full-RNA detection summaries, uncalibrated review screens and
actual same-cell joint-marker programs. Marginal detection bounds and actual
co-detection remain distinct; missing features retain NA. Anonymous synthetic
CLI checks are documented separately from the original biological audit.

The [frozen PSO doublet replay](qc/PSO_REPLAY.md) requires six independently
verified physical-capture inputs and pinned package versions. Its adapter
records actual classifier training, rejects training failures and preserves
the source counts and barcodes. It does not recover historical QC or infer
capture/donor identities from specimen names.

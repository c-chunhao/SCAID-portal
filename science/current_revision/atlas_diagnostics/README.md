# Final-atlas diagnostic and marker-review scripts

Use the [independent Python 3.12 environment](../INSTALLATION.md). The
package-specific `tested_versions.json` records the original analysis runtime;
it is separate from the clean installation environment.

These scripts diagnose existing atlas representations and summarize previously measured marker evidence. They do not integrate new data, replace labels, apply new QC filters, or claim independent annotation accuracy.

## Synthetic example

Install `requirements.txt`, then run:

```bash
SCAID_DIAGNOSTIC_DIR=synthetic_example python replay_sampled_diagnostics.py SYNTHETIC
SCAID_DIAGNOSTIC_DIR=synthetic_example python verify_final_diagnostics.py SYNTHETIC
```

Every coordinate, label and identifier in `synthetic_example/` is synthetic. No clinical table, real cell metadata or expression matrix is included in this code package.

## Real-data input and provenance

`export_final_reductions.R ATLAS OUTPUT_DIR SOURCE_RDS` reads a supplied final Seurat object and exports its retained PCA/Harmony reductions, metadata and commands. ATLAS is `T_NK` or `B_Plasma`; expected final member totals are verified. R 4.3.3, SeuratObject and jsonlite were used.

`analyze_final_latent.py ATLAS` checks the exported metadata against the previously extracted final member list, joins source-specific tissues through Supplementary Table S4 and calculates conditional neighbourhood diagnostics. Set these paths:

- `SCAID_DIAGNOSTIC_DIR`: exported representations and output directory.
- `SCAID_REFERENCE_METADATA_DIR`: original final-member metadata files.
- `SCAID_SUPPLEMENTARY_S4`: specimen denominator workbook.

The diagnostic supplement supplies sampled 20-dimensional coordinates and sampling records. `replay_sampled_diagnostics.py` reproduces all reported source-mixing estimates from those sampled coordinates. `verify_final_diagnostics.py` independently checks selected strata using SciPy pairwise distances. Biological-label concordance describes existing labels and is not independent accuracy.

`summarize_existing_annotation.py` reconciles the existing 63-label source-conditioned marker results and targeted same-cell co-detection results. Set `SCAID_MARKER_EVIDENCE_DIR` and `SCAID_INTEGRATION_REVIEW_DIR` to the evidence directories. It retains unavailable-gene states and does not reassign production labels.

The measured counts, sampled cell-level evidence and original expression matrices are supplied or referenced in the scientific supplement; they are intentionally absent from this source-only package. `tested_versions.json` records the environment used for the actual diagnostics. The real-data output provenance records identify source files and hashes.

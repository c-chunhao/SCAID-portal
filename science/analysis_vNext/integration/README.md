# Versioned integration and annotation safeguards

This directory is the public code copy of the local integration audit dated 2026-10-05. It contains executable modules, reproducibility instructions and guard validations without cell identities. Full audit evidence, raw matrices, reference models, clinical tables and per-cell tables are outside this code directory.

`released_rna.py` explicitly separates the full raw-RNA and legacy scaled-HVG feature axes and strictly validates CSR structure and nonnegative exact integers. Genes are selected only after full-RNA library normalization to 10k followed by log1p. Legacy dense-HVG reading rejects unknown column order by default and requires an independently verified source `scale.data` feature axis. Passing the H5 var index back into the reader does not constitute an independent source.

`integration_safeguards.R` attaches embeddings and joins labels by unique cell ID. Its Seurat v5/v4-compatible RNA counts getter rejects reading split layers as a combined layer, fixes the UMAP seed parameter, and exports a source-SHA sidecar from the actual axes of counts and `scale.data` separately. `assemble_B_Plasma_vNext.R` fixes the historical base-variable mix-up and refuses to overwrite an existing output RDS.

`reannotate_t1d_full_rna_pilot.py` provides a CLI for candidate annotation using sparse full RNA in 5000-cell blocks and frozen models. Legacy/source labels are review columns, not model inputs. The related High/Low models cannot serve as independent validation; their agreement and sigmoid scores are not accuracy measures. Candidates do not automatically replace production labels.

Completed local checks include a T1D pilot with 75746 cells; paired HVG/full-RNA sensitivity for 53938 cells; chunk/model/column invariance for 40 cells; zero differences from the independent source in raw values for 16 cells × 20690 genes; and exact agreement between all 2000 T1D scaled-axis names and the archived postprocessed RDS, with differences in 16 × 2000 scaled values ≤5.33e-15. This T1D evidence does not automatically authorize other release objects.

See [RUNBOOK.md](RUNBOOK.md) for execution entry points, [method_protocol.vNext.json](method_protocol.vNext.json) and [schema_sidecar_specification.md](schema_sidecar_specification.md) for methodological boundaries, and [DEPRECATED_ANALYSIS.md](DEPRECATED_ANALYSIS.md) for risks in deprecated historical scripts. UMAP is a display; batch/biology validation requires high-dimensional analysis and independent evidence.

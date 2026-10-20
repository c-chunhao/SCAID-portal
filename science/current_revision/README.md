# October 2026 scientific follow-up

These packages contain the current-revision diagnostic and B-state follow-up
scripts alongside their explicitly synthetic installation checks. Their
dependencies, real-data input requirements and tested versions are documented
within each directory. They complement the existing `science/scripts` and
`science/analysis_vNext` entry points; they do not replace the production atlas
or supply raw-data inputs.

- [atlas_diagnostics](atlas_diagnostics/README.md): conditional source-mixing
  diagnostics of the final PCA/Harmony representations, independent numerical
  checks and reconciliation of 63-label marker/co-detection evidence. These
  diagnose an existing integration and do not measure independent annotation
  accuracy. Real sampled cell-level evidence is excluded from this source
  package.
- [bcell_state](bcell_state/README.md): source-recorded B-state composition and
  within-state IFN12 follow-up in GSE174188, fixed original individual/batch
  selection, 40 recorded models (36 estimable), explicit label-coverage checks,
  group-level outputs and independent HC3/BH numerical verification. The
  installation demonstration uses entirely synthetic people/counts. No real
  individual records, clinical covariates, barcodes, raw counts or cell mappings
  are included.

The supplemental diagnostics and source-recorded states support conditional and
within-source interpretation. They do not constitute independent cohort
replication or causal evidence. Source availability and third-party reuse terms
follow the repository's existing licence notices; no licence is assigned here.

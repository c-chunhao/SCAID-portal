# Environment and reproducibility

The runtime actually used in this round is recorded in `observed_python_environment.json` and the R/QC manifests. That file is evidence of the environment used, not a lockfile guaranteed to resolve again with a single pip command. The local machine already has a shared scientific environment; its software versions are not changed.

- The composition entry point provides explicit `--input-tables`, `--source-map` and `--output` arguments. Inputs are read-only, and outputs require a new directory. Outputs include file SHA values, the immune vocabulary, coverage rules and summaries for all six layers.
- The original release/integration metadata were actually rerun, with numerical regression checks for Count, immune denominator, NA mask and fraction across all six layers; the work was not limited to editing the final Excel file.
- RNA reading uses the independent gene axis for all raw counts. Compute the full library size before normalization/log1p and any optional gene subset; scaled HVG data are not RNA counts.
- CellTypist models must be loaded from frozen local files with verified SHA values. The pilot manifest records the current model snapshot. The two related models are not independent gold standards.
- New QC uses real technical captures, with scDblFinder preceding strict tail QC; all counts and flags are retained. The classifier adapter in the source must perform real training and abort on failure, rather than label every cell a singlet after an error.
- Clinical sample labels within an object, source pseudonyms and capture labels must be stored separately. Without an individual mapping, records cannot enter donor-level disease inference.
- Newly installed versions cannot be used to fabricate the software versions behind historical UMAPs or preprocessing. Current checks also cannot recover historical filtering logs that were never saved.
- Local evidence replay scripts have explicit fixed input directories. Public reuse starts from the portable CLI and runbook; there is no claim that every replay script supports arbitrary data sources.
- All large scientific matrices are read sparsely or in blocks. Limit threads and memory on the shared host; do not convert the entire large integrated T/NK RDS to a dense matrix.

The original primary data, earlier articles and previous figure versions are all retained. The only change to the official download data is restoration of the Lyme versioned tail; the original base H5 and source RDS are not overwritten.

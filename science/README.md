# Scientific recomputation

The scripts correct complete-zero catalogue composition, compare final broad
lineages with source-author labels, and perform cell-type pseudobulk and
individual-level SLE analyses. `INPUTS.json` lists the required input-path
environment variables for each script. Provide the original repository inputs
and the released H5 objects under their applicable access and reuse conditions.
Cell/donor metadata, matrices and clinical tables are not redistributed here.
Create `tables/` and `figures/` output directories before running scripts.

Use `extract_release_metadata.R` and `extract_integrated_metadata.R` before
`recompute_atlas_composition.py`. The source-audit directory must contain the
verified `source_cohorts/sample_source_mapping.tsv`; missing mappings must not
be guessed. Source-author comparisons require exact barcode and donor joins.
SLE processing uses `extract_sle_pseudobulk.py`, `fit_sle_donor_models.py`, then
`source_donor_composition_and_IFN.py`. Run the latter first before rendering
`build_verified_donor_figure5.py`. All original inputs are read only.

The original composition ANOVA is withdrawn. Catalogue medians are descriptive;
625 specimen labels include pools and repeat observations. One specimen has no
immune-cell denominator and contributes undefined fractions rather than false
zero fractions. In the source-study SLE comparison, one sequencing batch is
selected per biological source individual before outcome testing, and models
adjust for age, sex, ancestry and processing cohort. This reproduces a known
interferon pattern; it does not establish a novel cross-disease mechanism.
Source-author reference concordance is limited to coarse lineages in two
objects, is not blinded annotation accuracy, and includes unknown predictions
as disagreements. Third-party inputs are subject to their own licences.

Dependencies include pandas, numpy, scipy, h5py, matplotlib, scikit-learn,
statsmodels, PyDESeq2 0.5.4, PyMuPDF and R/Seurat/Matrix. Web-runtime
requirements do not constitute a locked scientific environment. Full methods,
covariate design and gene-level results accompany the manuscript supplement.

# Additional source-pool sensitivity

This package adds processing-pool sensitivity to the existing individual-level
GSE174188 analyses. The source `batch_cov` variable identifies multiplexed
pools, nested within processing cohorts. The original participant/pool
selection and all original HC3 estimates are retained. The additional models
use CR1 covariance clustered by source pool with t(G-1) inference; one variant
also uses pool fixed effects. The analysis reads derived individual summaries,
not expression/count matrices. It is a within-source sensitivity analysis, not
independent replication or a causal model.

## Install and run the synthetic example

Use the separate Python 3.12 environment and exact tested lock in
[INSTALLATION.md](../INSTALLATION.md). Then, from this directory:

```sh
python scripts/synthetic_pool_demo.py --output work/pure_synthetic_pool_example
```

The example invents 200 people and 20 pools, builds every required input, runs
all 56 planned model rows (54 estimable), and independently checks QR/CR1, t
intervals and BH arithmetic. None of the example inputs or IDs derives from a
real participant.

## Real-input replay

Set the following variables to privately supplied/reconstructed inputs:

- `SCAID_BASELINE_SCIENCE` (default `data/baseline_science`): `tables/` contains
  `GSE174188_unique_donor_composition_complete_zeros.tsv`,
  `GSE174188_target_cell_IFN_by_unique_donor.tsv`,
  `GSE174188_all_donor_batch_inventory.tsv`, and
  `GSE174188_unique_donor_HC3_composition_IFN_tests.tsv`. The existing
  `science/scripts` entry points define these source-derived tables.
- `SCAID_BCELL_OUTPUT` (default `work/bcell_state`): the output of the real
  [B-state workflow](../bcell_state/README.md), with `tables/` containing
  `B_individual_composition_and_IFN.tsv`, `B_state_model_input_rows.tsv` and
  `B_state_HC3_models.tsv`.
- `SCAID_POOL_OUTPUT` (default `work/source_pool_sensitivity`): a new private
  local output directory.

```sh
python scripts/analyze_source_pool_sensitivity.py
python scripts/independent_pool_numeric_QA.py --directory "$SCAID_POOL_OUTPUT"
```

Export the variables before the shell command above; Python defaults are not
shell exports. Real individual rows and pool-level leave-out records are
excluded from this public package. Output files prefixed `PRIVATE_`, the input
records and local output directories must not be committed. Aggregate models,
eligibility, leave-one-pool-out summaries, source-field semantics, the additional
analysis plan and numeric QA are included. Approximately 20 retained source
pools provide limited finite-cluster inference; t/CR1 is an approximation.

The source-publication and author-code DOI in the semantics record identify
third-party source evidence; they are not an archive DOI for SCAID. Project and
third-party rights follow the repository notices; no new licence is assigned.

# Source-recorded B-cell states and individual-level IFN12

This package contains the scripts used for the within-GSE174188 B-state follow-up
and aggregate results. It does not contain clinical tables, individual-level
records, original barcodes/count matrices, cell-ID mappings or scan binaries.
It does not run the real analysis without separately supplied source data and
baseline inputs, and does not establish cell-intrinsic activation or independent
replication. Project and third-party licence status follow the repository notices;
this package does not assign a new licence.

## Synthetic installation check

Install `requirements.txt`; a C++17 compiler and zlib development library are
needed for streaming. Run `python synthetic_example/run_synthetic.py`. It creates
a tiny synthetic MTX and validates selected-cell library totals/watch counts, then
fits HC3 models on 200 entirely synthetic people. It is an installation and
mathematical check, not a miniature or anonymized copy of the real SLE results.

## Real-data inputs

Set these variables explicitly, or use the relative defaults shown:

- `SCAID_GSE174188_SOURCE=data/GSE174188`: contains
  `GSE174188_metadata.csv` and `10x/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz}`.
- `SCAID_BASELINE_SCIENCE=data/baseline_science`: contains the four baseline
  tables listed below under `tables/`.
- `SCAID_BCELL_OUTPUT=work/bcell_state`: private local output; its `tables/` and
  `provenance/` contain individual-level derived files and must not be committed.

Required baseline tables:

1. `GSE174188_frozen_one_batch_per_donor.tsv`: 233 preselected individuals,
   `ind_cov`, `batch_cov`, `unit_id`, `cohort`, `condition`, `n_cells`, `Age`,
   `Sex`, `pop_cov`, `Status`. Selection is maximum total processed source-cell
   recovery in cohorts 2–4 per person, lexical batch tie break, before B/IFN tests.
2. `GSE174188_target_cell_IFN_by_unique_donor.tsv`: original lineage-level
   scores including `ind_cov`, `cg_cov`, `IFN12_log1p10k_mean`, `n_cells`.
3. `GSE174188_all_donor_batch_inventory.tsv`: source-individual/batch inventory
   including `ind_cov` and source `Status`, used for within-source sensitivity.
4. `GSE174188_unique_donor_HC3_composition_IFN_tests.tsv`: original aggregate
   coefficients/CIs/P values, used to confirm the original main models unchanged.

The first three contain source-individual information and are not included.
The repository's `science/scripts/extract_sle_pseudobulk.py` defines the original
selection and score calculations from supplied source metadata/counts;
`source_donor_composition_and_IFN.py` defines original aggregate tests. Reconstruct
or supply the baseline records under the source repository's access conditions;
do not replace them with invented rows or relax identity checks. The fourth
group-level reference is included in `aggregate_results/` for inspection.

## Real-data commands

```sh
python scripts/prepare_bcell.py
c++ -O3 -std=c++17 scripts/stream_B_watch.cpp -lz -o work/stream_B_watch
work/stream_B_watch "$SCAID_GSE174188_SOURCE/10x/matrix.mtx.gz" \
  "$SCAID_BCELL_OUTPUT/provenance/B_column_index.int32" \
  "$SCAID_BCELL_OUTPUT/provenance/watch_gene_index.int32" \
  "$SCAID_BCELL_OUTPUT/provenance/selected_B"
python scripts/analyze_bcell_states.py
python scripts/build_readable_results.py
```

Export the variables before using these shell commands; Python's relative
defaults do not export shell variables. Inputs are the archived processed/remapped
source matrix, not unfiltered droplets. The real-analysis assertions deliberately
require the published dimensions and original model agreement.

## Outcome and coverage definitions

Source `ct_cov` states are used within `cg_cov=B`; they are not learned from IFN
expression. IFN requires >=20 target-state cells per individual and >=10 eligible
people in each arm. Fractions retain all B cells in the denominator, including
unresolved labels. Entirely unavailable named-state labels are not state absence.
Full-rank designs, >=10 residual degrees of freedom and stable positive HC3
variance are required. BH is applied by each recorded outcome family; unavailable
planned tests enter as P=1 and remain unavailable. See `analysis_plan.json` and
`coverage_design_addendum.md` for the exploratory/status/coverage qualifications.

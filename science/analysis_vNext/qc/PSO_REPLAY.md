# Frozen PSO doublet replay

`run_pso_scDblFinder_pilot.R INPUT_ROOT [CALLING_STAGE] [ADAPTER_PATH]`
replays the six preverified GSE162183 physical-capture inputs used in the
October 5 audit. INPUT_ROOT/doublets_input contains each source gene-by-cell
sparse MatrixMarket matrix and its _barcodes.tsv with a cell_id column.
The reviewed channels are Ctrl1_10x, Ctrl2_10x, Ctrl3_10x, Psor1_10x,
Psor2_10x and Psor3_10x. These labels do not establish individual identities.
No actual matrices or cell records are distributed in this code module.

The process-local adapter is the exact source used by the recorded pilot.
This replay explicitly requires scDblFinder 1.16.0 and xgboost 3.2.1.1,
along with Matrix, SingleCellExperiment, BiocParallel, jsonlite and digest.
It changes modern API contracts, records real training and prediction
changes, and makes classifier failures fatal. Installed packages on disk
are never modified. Other versions and source designs need their own
validated implementation, not an untested compatibility assumption.

The primary pilot calls doublets before strict distribution-tail QC;
CALLING_STAGE records the separate order-sensitivity replay. All original
counts and barcodes are retained, with separate QC/doublet flags. The
source partition was established separately and is an input prerequisite;
this program cannot infer capture or donor identity from barcode suffixes.
It does not automatically approve any new cohort for publication.

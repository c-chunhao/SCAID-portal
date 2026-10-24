# Reproduction and verification

Read [INPUTS.md](INPUTS.md) and the recorded environment first. Work in a fresh
directory with the pinned scientific runtime and enough disk/memory. The original
Linux run used R 4.3.3, up to four 80-GiB QC jobs, and eight integration threads.
The final sparse CellChat process peaked near 25 GiB RSS; its original cache and
full-axis verification reserve guards wait for at least 320 GiB system-available
memory. Those guards remain in the recorded code. The web Python environment is
separate from the scientific environment.

Activate the appropriate R library with `R_LIBS_USER` before the frozen-RNA
commands. They reject mismatched recorded R/package versions. CellChat additionally
checks a canonical fingerprint of the installed official function and the saved
human DB bytes; a version label alone is insufficient.

```sh
cd science/ad_gse332623_20261024
python reproduce.py plan
python reproduce.py materialize --root /your/new/ad-replay \
  --raw-input /your/extracted/GSE332623 --r-lib /your/R/4.3-library \
  --web-python /your/legacy-anndata-env/bin/python \
  --celltypist-models /your/celltypist-models
```

`materialize` replaces only the documented root tokens in source templates,
writes local source hashes and refuses to replace changed materialized code.
The published [execution map](SOURCE_EXECUTION_MAP.json) retains original and
public hashes. Historical receipts stay in `evidence`; they are not silently
re-signed or accepted as fresh-run scientific approval.

## Frozen-RNA replay

The portable entrypoint runs the recorded scientific R kernels on the exact
published RNA. This separates scientific logic from the production launcher's
original PID, path and old receipt bindings. An existing output directory is
rejected; input RNA and KEGG checksums must match before computation begins.

```sh
python reproduce.py scores --rna /your/RNA.h5ad \
  --genesets /your/20240526_kegg_hsa_list.rds \
  --output /your/new/score-replay --workers 8
python reproduce.py cellchat --rna /your/RNA.h5ad \
  --output /your/new/cellchat-replay
```

`scores` computes the baseline with at most four workers, then all final corrected
scores with eight workers by default. It uses the full 33,660-gene universe,
1,000-cell chunks, seed `2026 + chunk_index - 1`, the explicit RNA-positive AUCell
mask and the accepted resource mappings. It independently calls official AUCell,
UCell and singscore on the first 1,000 cells and requires float32 equality. Outputs
are the full raw score HDF5 matrices and a new receipt. This entrypoint does not
recreate the original H5AD packaging/metadata bytes; the actual export functions
are retained in `run_pathways.py` for inspection and explicit adaptation.

`cellchat` first uses the unchanged preprocessing prefix of the selected R script
to derive DB gene matching and the same inclusion mask. The recorded large-cohort
index builder reconstructs a full source-preserving RNA CSC index with complete
component readback. The actual DBcache builder then copies every selected column,
validates every component and cell/library/pool axis, retains the 48 zero-nonzero
columns and checks direct RNA probes. Only the regenerated index hash and local
path bindings are supplied as runtime glue. The cache retains 581,725 records;
the unchanged modeling mask selects 528,739 across eight groups.

Before modeling, the recorded independent full-axis script rechecks sparse versus
official triMean for the observed grouping and first two seed-2026 permutations.
The unchanged selected R kernel then computes all 100 permutations and renders
the nine figure families. The historical 100-permutation pilot proof is supplied
as dated evidence, not represented as a newly rerun pilot. New replay receipts
describe their narrower checks. Exact DB bytes and installed-function identity
must be checked against the recorded environment; changed resources require a
new implementation review. A separate full dense formal network is not claimed.

## Raw MEX to corrected RNA/ADT

The raw-stage scientific scripts are supplied for a new analysis. After all
81 source-member checks pass, copy the minimal inventory projection to the
materialized expected path and ingest each library:

```sh
cd /your/new/ad-replay
mkdir -p project/code/v2/00.Reference
cp resources/SampleInventory.GSE332623.tsv project/code/v2/00.Reference/SampleInventory.tsv
Rscript project/code/GSE332623_build_corrected_input.R \
  GSE332623_TotalseqB_NEA1_ADT1 project/00.corrected_data/GSE332623
```

Repeat for all 27 library IDs in `resources/capture_manifest.tsv`. The reader
preserves RNA/ADT axes and integer counts and writes source/RDS receipts. Rebuild
the capture manifest with the newly created RDS paths and hashes; the historical
corrected-RDS hashes are not assumed to remain byte identical. Preserve manifest
order and technical-library/unresolved-donor fields.

Run `GSE332623_qc_preflight.R <library> <output>`, then
`GSE332623_rna_qc.R <library>` and `GSE332623_finalize_qc.R <library>` for each
library. `GSE332623_run_qc.py` is the recorded Linux/systemd bounded dispatcher;
per-library R commands can also be scheduled independently with the documented
resources. Summarize and validate bridges with the supplied QC scripts. decontX
contamination flags never remove cells; its fractional corrected layer is kept
apart from the original integer count layer.

The integration pipeline has four explicit stages:

```sh
python work/analysis/code/analysis_pipeline.py assemble --out /your/integration \
  --bridges /your/bridge1.h5 /your/bridge2.h5
python work/analysis/code/analysis_pipeline.py classify --out /your/integration
python work/analysis/code/analysis_pipeline.py embed --out /your/integration
python work/analysis/code/analysis_pipeline.py finalize --out /your/integration
```

The `assemble` example abbreviates the list: supply **all 27 bridge files in
capture-manifest order**, not only the two shown. Place the verified CellTypist
model files in the materialized model directory. Run the included reference
remapping/recheck, coarse-marker review and renderer before downstream selection.
The authoritative final run selected PCA; fresh data or versions must be evaluated
with the gate, rather than forcing the historical decision.

Raw-stage execution requires fresh receipts, resource gates, source matching and
independent review. The original production pathway/CellChat Python launchers
intentionally bind historical source/cache/approval hashes. They are preserved
for provenance and should not be made to pass by overwriting the dated receipts.
This release did not run a clean raw-to-final reproduction or create a complete
recreated lockfile. The unbundled historical KEGG/model/DB bytes and unrecorded
dependency details remain explicit reproducibility limitations.

# GSE332623 recorded environment, actual parameters and reproducibility limits

Audited 2026-10-24. Work root: `@WORK@`; project root: `@PROJECT@`. This audit changed no original source, input or result. `SOURCE_EXECUTION_MAP.json` supplies original-source and published-file SHA256. The local audit checked 26 current-source/final-manifest bindings, all matching. These are final source revisions; intermediate historical revisions were not all preserved.

## Actual final analysis

- Ingest: 27 published MEX libraries, 754,086 barcode records, 33,694 RNA features plus 140 ADT features in every library. `features.tsv.gz` third column separates Gene Expression and Antibody Capture. The 27 technical libraries are operational QC partitions, not verified individual captures or donors. Cross-library identical barcode strings are retained.
- QC: 725,850 basic-pass records; scDblFinder calls 144,125 doublets; final 581,725 retained records. RNA-based filters and doublet calls exclude records; decontX supplies a review flag only (6,724 retained records flagged), with fractional RNA_decontX kept separately. Raw integer RNA/ADT counts drive downstream analyses.
- Integration: exact original-symbol duplicate aggregation 33,694 -> 33,660 rows (34 duplicate symbols), per-cell UMI conserved; RNA total 10,000 then natural log1p; ADT separate per-cell CLR. Batch-aware 3,000 HVGs plus protected embedded lineage/state panels yields 3,015 selected integration features, 92 protected features present, scale(zero_center=FALSE,max_value=10), then centered 50-PC arpack PCA; mitochondrial/ribosomal nuisance genes are excluded unless present in protected panels.
- Harmony was computed on library_id with theta=2, lamb=1, nclust=100, max_iter_harmony=15, max_iter_kmeans=20, random_state=2024; converged after 8 iterations. Final representation is PCA. Neighborhood30, Leiden0.8 (igraph, n_iterations=2), 24 clusters, UMAP min_dist0.3/spread1/maxiter300, seed2024.
- PCA selection gate: Harmony entropy gain >0.02, coarse purity loss <=0.02 and coarse purity >=0.90, fine purity loss <=0.03, within-lineage state-score neighbor Spearman loss <=0.05. Actual early-response loss 0.1756013708 fails the state gate; both PCA/Harmony are saved. Diagnostics sample up to600 confident reference-labeled cells/library (16,200 total), seed2024.
- CellTypist High32 and Low98 immune models run in10,000-cell chunks. High model SHA256 a715fec36c2c421f7c4e31cf4cb4bea883eedab7fd7b20a4b76f091c22660448; Low model SHA256 290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502. Coarse labels require reference evidence, High probability>=0.5 and RNA-marker support. ILC is compatible only with Low-model NK and is distinct from exact coarse agreement; Cycling is a state. Coarse RNA support is cluster mean panel detection>=0.08; nonimmune epithelial/fibroblast/endothelial candidates also need>=2 panel genes detected in>=10% of cluster cells. Candidate confidence=High probability*min(1,RNA support/0.3). Fine labels need Low probability>=0.7, candidate confidence>=0.5, cluster modal subtype fraction>=0.35, per-cell subtype equal to mode, and subtype panel mean detection>=0.05 or any marker>=0.2; predefined unsuitable tissue/progenitor subtypes are excluded. Source-library dominance review>=0.9. ADT exploratory AUROC thresholds0.65 relative support/0.55 limited enrichment with>=5 cells each side per library. All raw predictions and uncertainty remain recorded.
- Final pathway branch is `RNA_scores_zero_mask_and_alias_v1`. RNA log1p(cp10k) is scored over all33,660 RNA genes, without ADT or re-normalization. 359 KEGG sets from historical20240526 resource; overlap1..500 gives356 scored sets (three broad sets excluded), with16 low-coverage (<10 genes) sets flagged. AUCell cap ceil(0.05*33660)=1683, explicit RNA>0 rank-weight mask; UCell2.6.2 cap1500/average descending ties; singscore centerScore=FALSE with min-tie rankGenes. Chunk1000, 582 chunks, seed2026+chunk_index-1, explicit SerialParam RNGseed. Baseline scoring used4 workers; final zero-mask/alias branch8 workers (not the baseline receipt worker value).
- Resource matching uses1,172 verified collision-free symbol renames, direct source hits first and one physical row counted once.247 aliases occur in all359 KEGG sets (244 in356 scored sets), restoring799 member entries in159 scored sets. UCell/singscore recomputed159 columns and retained197 bitwise baseline columns. Final AUCell has zero scores for all22,838,176 no-detected-member cell/set pairs. Independent1000x356 official reference maximum errors0.
- Final CellChat is `cellchat_sparse_exact_candidate`:528,739 records,8 groups; excludes52,809 Unassigned,172 provisional source-dominated Fibroblasts, and5 cells in groups<10 (4 Mast+1 Neutrophil). All581,725 remain in RNA/pathway outputs. Human DB3233 interactions,1434 matched RNA genes (9 verified internal aliases); all original normalization values preserved and no DB-subset re-normalization.
- CellChat parameters: identifyOverExpressedGenes(do.fast=TRUE,min.cells=10,only.pos=TRUE,thresh.pc=0,thresh.fc=0,thresh.p=0.05); triMean/raw.use=TRUE/population.size=FALSE/distance.use=FALSE/nboot=100/seed.use=2026. PPI projection absent. Strict nominal p<0.05 edges; unfiltered probabilities preserved. One calculation worker and sequential centrality seed2026. Final660 significant communication records,79 unique LR pairs,42 inferred pathways,350 rendered PNGs.
- Sparse CellChat replaces only group aggregation in the installed official computeCommunProb AST; finite nonnegative sparse-positive type7 quartiles use collapse2.1.6 arithmetic and unchanged R mean(Q25,Q50,Q50,Q75). All probability/P/dimname arrays are bitwise identical in a17,536-cell100-bootstrap pilot; formal1434gene/8group aggregation matched observed plus first2 permutations (34,416 values). This is bounded implementation-equivalence evidence, not a separately completed full dense formal reference.
- Marker tests are within-cohort annotation explorations, not donor/AD-HC tests: Wilcoxon with tie correction, BH within group, fixed seed2024 stratified samples documented in per-table methods JSON (up to1000 cells/coarse group,7988 tested cells; up to200 cells/Leiden group,4800 tested cells). Enrichment is optional downstream: positive log2FC>1/BH P<0.05/detection>=0.1, Unassigned excluded; background29,936 RNA symbols detected in>=3 cells, mapped GO/KEGG effective backgrounds16,404/7,803, GO BP, minGS10/maxGS500.

## Recorded environment (not a recreated lockfile)

Python executable for integration: `@ANALYSIS_PYTHON@`; logs show Python3.13 site-packages but do not record a patch version. Legacy H5AD readback used `@WEB_PYTHON@` with anndata0.11.3. Main analysis records:

| Python distribution | Recorded version |
|---|---|
| numpy | 1.26.4 |
| scipy | 1.16.3 |
| pandas | 2.3.3 |
| anndata | 0.12.17 |
| scanpy | 1.11.5 |
| harmonypy | 0.2.0 |
| celltypist | 1.7.1 |
| scikit-learn | 1.8.0 |
| igraph | 1.0.0 |
| leidenalg | 0.11.0 |

QC explicitly prepended R library `@R_LIB@`; downstream Rscript inherited its runtime/library environment. R4.3.3 is recorded. QC dependency record:

| R / package | Recorded version |
|---|---|
| R | 4.3.3 |
| Seurat | 5.4.0 |
| SeuratObject | 5.3.0 |
| Matrix | 1.6.5 |
| SingleCellExperiment | 1.24.0 |
| scater | 1.30.1 |
| scuttle | 1.10.3 |
| scran | 1.28.2 |
| scDblFinder | 1.16.0 |
| celda | 1.18.2 |
| BiocParallel | 1.36.0 |
| ggplot2 | 4.0.2 |
| jsonlite | 2.0.0 |
| rhdf5 | 2.46.1 |
| glmGamPoi | 1.14.3 |
| irlba | 2.3.7 |
| igraph | 2.2.2 |
| bluster | 1.10.0 |
| DropletUtils | not installed / not required by this chain |

Pathway recorded versions: AUCell=1.24.0, UCell=2.6.2, singscore=1.22.0, Matrix=1.6.5, hdf5r=1.3.12, BiocParallel=1.36.0. CellChat recorded versions: CellChat=2.2.0.9001, Matrix=1.6.5, hdf5r=1.3.12, future=1.70.0, presto=1.0.0. Sparse prototype additionally records CellChat=2.2.0.9001, collapse=2.1.6, Rcpp=1.1.1. Enrichment: clusterProfiler4.10.1, org.Hs.eg.db3.18.0.

Unrecorded dependencies used by source include exact Python patch, h5py, torch, numba, pynndescent, psutil, Pillow; and R matrixStats, ComplexHeatmap, digest, jsonlite transitive versions beyond the QC observation. No conda/renv/container lockfile was found in this scoped run. Do not invent pins from today's installation.

Resource controls: ingest2 processes with BLAS/OMP1; QC up to4 systemd user scopes,80GiB each/CPUQuota400%, OPENBLAS1/OMP4/MKL1/R_MAX_VSIZE72Gb; integration BLAS/OMP/MKL/NUMBA8; score workers8 with worker BLAS1; selected CellChat1 worker with BLAS1, CPU affinity10 and RLIMIT_AS soft128GiB. Historical large reference resource budgets describe superseded reference processes, not final sparse computational requirements.

## QC thresholds and stochastic settings

For each library x, lower thresholds = max(floor, floor(10^(median(log10(x+1))-3*MAD)-1)), with UMI floor200 and genes floor100. Upper thresholds = max(ceil(q99.9 raw x),ceil(10^(median(log10(x+1))+4*MAD)-1)). MT max=max(15,min(30,median(percent_mt)+4*MAD)). Ribosomal/hemoglobin/ADT are review signals only.

scDblFinder: seed332623+1-based sorted library index (332624..332650), L'Ecuyer-CMRG/SerialParam seed, clusters=TRUE,nfeatures=1352,dims=20, dbr/dbr.sd use package defaults. If >50,000 basic-pass cells, deterministic partitions of<=25,000 are used (reported library basic-pass counts do not trigger this condition in this run). Default doublet prior is a technical-library model assumption, not demultiplexing truth. decontX: nonzero RNA rows of all basic-pass cells,maxIter500,seedQC+1, no supplied empty-droplet background. ambient cutoff=max(0.5,min(0.9,median(contamination)+4*MAD)); flag if strictly above; never excludes.

| Library suffix | UMI min | UMI max | Genes min | Genes max | MT max % |
|---|---:|---:|---:|---:|---:|
| NEA1_ADT1 | 261 | 39496 | 176 | 16336 | 15 |
| NEA1_ADT2 | 308 | 63138 | 189 | 21002 | 15 |
| NEA1_ADT3 | 223 | 25063 | 133 | 11190 | 15.0472900211825 |
| NEA1_ADT4 | 220 | 23677 | 133 | 10636 | 15 |
| NEA2_ADT1 | 368 | 52645 | 235 | 17461 | 15 |
| NEA2_ADT2 | 321 | 47854 | 220 | 15747 | 15 |
| NEA2_ADT3 | 286 | 55406 | 180 | 19941 | 15 |
| NEA3_ADT1 | 967 | 46607 | 568 | 11775 | 15 |
| NEA3_ADT2 | 885 | 61155 | 505 | 14649 | 15 |
| NEA3_ADT3 | 649 | 60610 | 399 | 15800 | 15 |
| NEA3_ADT4 | 854 | 49441 | 486 | 13141 | 15 |
| NEA4_ADT1 | 894 | 29809 | 575 | 7753 | 15 |
| NEA4_ADT2 | 900 | 35087 | 600 | 8227 | 15 |
| NEA4_ADT3 | 751 | 36211 | 523 | 8518 | 15 |
| NEA4_ADT4 | 978 | 27105 | 614 | 7160 | 15 |
| NEA5_ADT1 | 996 | 47075 | 633 | 9268 | 15 |
| NEA5_ADT2 | 593 | 63725 | 455 | 13146 | 15 |
| NEA5_ADT3 | 615 | 47104 | 437 | 11454 | 15 |
| NEA5_ADT4 | 610 | 50971 | 419 | 12559 | 15 |
| NEA6_ADT1 | 410 | 85574 | 358 | 15583 | 15 |
| NEA6_ADT2 | 751 | 34495 | 507 | 9063 | 15 |
| NEA6_ADT3 | 842 | 42851 | 577 | 9785 | 15 |
| NEA6_ADT4 | 767 | 50474 | 536 | 10887 | 15.6565318159248 |
| NEA7_ADT1 | 361 | 39022 | 273 | 12999 | 15 |
| NEA7_ADT2 | 378 | 49828 | 280 | 15131 | 15 |
| NEA7_ADT3 | 497 | 65468 | 375 | 16660 | 15 |
| NEA7_ADT4 | 369 | 54778 | 291 | 16127 | 15 |

## Exact execution references and minimal source set

The publication package includes the actual minimal ingestion/QC/integration/final-scoring sources plus selected independent-validation scripts; SOURCE_EXECUTION_MAP.json is the authoritative per-file list. The broad old v2 functions file is not an ingest dependency: builder sources only code/v2/GSE332623.multimodal.R. The required inventory projection contains 54 public rows and the three columns gse/gsm/sample_title. Marker/state panels are embedded in the integration scripts; no external marker GMT is consumed.

Actual ingress child command: `/usr/bin/Rscript <project>/code/GSE332623_build_corrected_input.R <library_id> <project>/00.corrected_data/GSE332623`; `build_all_AD_corrected_inputs.py` dispatches27 libraries. Actual QC child command: `Rscript GSE332623_rna_qc.R <library_id> && Rscript GSE332623_finalize_qc.R <library_id>` within recorded systemd scopes; first library was a pilot, then26 dispatched queue jobs.

Integration exact timed argv are in `analysis/logs/full_{assemble,classify,embed,finalize}.resources.txt`. They call `<project>/code/GSE332623_analysis/analysis_pipeline.py <stage> --out <analysis>/integration`; assemble additionally passes all27 `--bridges` in capture-manifest order. Source copies in work analysis/code match the project source directory. `remap_and_recheck_reference.py` recomputes reference diagnostics after mapping-policy fix; `extend_coarse_marker_review.py`/renderer add review evidence.

Pathway launch shape from actual source watcher plus final config: `python run_pathways.py --input <RNA> --mode formal --ready <analysis_ready.json> --output <analysis>/pathway_and_cellchat/pathways --chunk-size1000 --workers8` (options require separate tokens). Initial baseline had4 workers; final config/zero-mask receipt binds8 workers. Inside, actual R scientific kernel argv: `Rscript --vanilla run_pathways.R <scoring_config.json>` and `Rscript --vanilla run_aucell_zero_mask.R <zero_mask_config.json> <RNA_scores_zero_mask_and_alias_v1>`. Final scoring log contains a recovered HDF5 lock failure followed by validated completion; use final receipts, not absence of historical errors, as completion evidence.

Exact saved final CellChat full argv is `CellChat_sparse_candidate_plan.json` command and duplicated in audit manifest. It invokes run_cellchat_sparse_optimized.py --input<...RNA> --output<...cellchat_sparse_exact_candidate> --mode formal --ready<...analysis_ready> --nboot100 --plot-limit0 --workers1 --cache-receipt<...normalized_DB_expression_cache_receipt.json> --implementation-review<...exact_sparse_mean_implementation_review.json> --full-axis-review<...formal_axis_means.json>. The historical --reference-pid1887613 tied a coexisting standard reference process and is not a fresh-run argument. The R scientific kernel is `Rscript --vanilla run_cellchat_sparse_optimized.R <cellchat_config.json>`.

## Run-bound gates / inputs to regenerate for another workspace

1. Source MEX:81 files under `@RAW_INPUT@`, exactly27 matrix/features/barcodes triplets. Per-library ingestion receipts record all81 source SHA256, RNA/ADT feature axes and RDS output SHA. The reader feature-type split and exactly-once library prefixes must be preserved. Small manifest and27 receipts are suitable provenance; do not bundle large counts.
2. Capture manifest binds27 corrected input RDS hashes and GEO/SRA/BioSample identifiers; its SHA is16c1e2d1742c3d81f71c37b4bed2300b39e16abcd4ec9d194ce64131bf91a926. Each QC receipt binds corrected RDS SHA, raw-count bridge SHA, final filtered RDS SHA, seed, thresholds and versions. Recomputed byte hashes may differ even if scientific counts match; regenerate receipts for reruns.
3. Models: two CellTypist immune PKLs are external model assets identified by hashes above and model checkpoints. KEGG historical RDS has86,304 bytes/SHA ccb053e043dc678796520c6e7c0a36ad8f09bd4dd8673e582dcd309526e608de. CellChatDB_human_used.rds is2,087,991 bytes/SHA45536ae5df0ff6670a12ac97b3610a7455ef4ccdf8667346435131153a40835e. Publication of source code must not imply these exact historical resource bytes can be acquired from a current package URL.
4. Pathway approval gate: KEGG_alias_mapping_approved.json binds approved_symbol_renames.tsv SHA1c08833b7f9be7ffec22501a97506c803d87d34c1febfcede02fc8f288dee3fd/count1172 and actual_pipeline_validation.json SHA0bac1efde6e2a70a2c7744af9f6ad96f1255e11dcee47e75c6098053dc2cf12c. Independent AUCell zero-mask validation must have1000 cells/356 pathways/max_error0/no-zero-detection nonzero score. Preserve these historical receipts; adapt fresh path/hash mappings explicitly rather than rewriting them silently.
5. DBcache builder needs website/prepared/gene_entry.json -> full-source CSC gene index SHAa1bf380298c897810fa2a8fdac835251f437f42a10a6890c16999bae7ffdf0fa, source RNA SHAb9aa963123cb8f010d3e092e03507b70d02170afd93d32a6cdea2efe2a2c09ec, standard reference cellchat/CellChat_RNA_symbol_matching_view.tsv and cell_inclusion.tsv.gz,9-alias table SHA51873d0026a7cb6657298baf266f7de7c4f9b1b7ab132e4b2193c932703b18cd. The full cache is379,170,215 bytes/SHA557ebe5f3fcc9b9f744c5ba68455b6a25100d31ec03d15e678d82853144fa363,581725x1434,48 zero columns retained. This source builder is website-index-bound; it is not a standalone generic RNA-to-cache entrypoint. A clean cache can instead be reconstructed from the same RNA values/axes plus exact9-alias mapping and inclusion policy, with new receipts and independently verified equivalence.
6. Sparse launcher checks pilot100.json, full_axis_means.json, exact_sparse_mean_implementation_review.json, independent_sparse_helper_review.json; sourceRNA/cache SHA, exact528739cells/1434genes/8groups/seed2026, all bitwise equivalence flags, proof-path/hash crosslinks, official computeCommunProb source SHA and each reviewed script SHA. Path-only source transforms change script hashes and require an explicit original->public hash mapping or new verification; historical proof receipts must remain unchanged. Fixed candidate-directory restriction and process-coexistence binding are production run guards.
7. Final RNA SHAb9aa963123cb8f010d3e092e03507b70d02170afd93d32a6cdea2efe2a2c09ec; ADT SHA f23a4fb3b3e1812b478b80f635d017fb4a9e111686f000fccfc09fa683dfa7c4. Large data hashes in this audit are original receipt evidence, not a new full-data verification. All small sources/proofs/resources in the manifest were freshly hashed.

## Limits and publication claims

- No full clean raw-to-final replay was performed by this audit. Source templates, recorded command/config values, historical outputs and replay kernels must be described separately. No exact environment lockfile exists in this scoped run.
- Earlier reference-correction receipt records pipeline SHA c15345803afa6bdd4d179188cc2a8dc0a6599d90e1d7f7223bfa41967281ee3e; final analysis_ready/code_manifest/current pipeline SHA d46534d44f2b9d6700c52e9f26237f7a90ec2f64757ec4442630b547e8e4d6d8. Earlier extended-marker helper SHA cb73f320ea2e5d3b1df646aafff41185c29ea3057cd24a43da500a6918fe8d5a differs from final summary_helper SHA3953f642cac36f1d9cbe41ae7b74ce917e53587a5600af7c8d877ec15fffb3d3; both roles are preserved in the final marker review receipt. Final files match final manifests; do not claim final file bytes were the only versions used at every stage.
- Do not publish stale `pathway_and_cellchat/code_release_snapshot` as final code; it omits final zero-mask/sparse changes. Standard `cellchat` and earlier `cellchat_exact_candidate` were superseded before completion. Selected sparse result and final independent QA are authoritative.
- Axis hashes intentionally differ by encoding: integration axis_sha appends a final newline; pathway axis_digest joins without final newline; DBcache records JSON-array axis hashes. Compare values using the same encoding, not digest strings across protocols.
- No cell-to-donor/AD-HC mapping; donor/clinical labels remain Unresolved. Technical-library mixing and reference-neighborhood tests cannot establish preservation of unobserved clinical effects. No case-control, donor-level significance or confirmed AD-specific pathways/communication is supported.
- Rare Fibroblast/Mast/Neutrophil predictions are provisional; no external expert confirmation is claimed. ADT crosschecks are exploratory within-library cluster-vs-rest evidence. CellChat has nominal edge P values and describes pooled modeled populations.
- Small analytical source files were scoped for publication and a credential-assignment pattern scan found no hits. This is not permission to copy operator_access, backup, database snapshots, environment secrets, per-cell matrices, models or compiled caches.

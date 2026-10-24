# GSE332623 actual-run source and reproduction package

Version `0.3.1-ad20261024` records the analysis selected for the 24 October 2026
[mixed AD/healthy-reference PBMC browser record](https://scaid01.com/Dataset?record=176).
It contains actual scientific scripts, small reference mappings, original-source
hashes, path-redacted configurations and dated validation receipts. It contains
no expression matrices, cell tables, database dumps, model weights or credentials.

The RNA reanalysis retains **581,725 library–barcode records from 27 technical
libraries**. Donor identities and AD/control assignments remain **Unresolved**.
The cohort is available for descriptive browsing; it supports no AD-versus-control,
age-group or donor-level comparison. The two existing atlases remain unchanged:
937,825 T/NK cells and 253,635 B/plasma cells.

The selected embedding is **PCA**. Harmony was computed and rejected by the
recorded state-preservation gate. Final pathway scores use the
`RNA_scores_zero_mask_and_alias_v1` branch, including the explicit AUCell zero
mask and 1,172 approved resource-symbol mappings. Final communication inference
uses `cellchat_sparse_exact_candidate`: 528,739 modeled records, 1,434 matched
genes, eight populations, 100 permutations and seed 2026. Earlier baseline score
files and incomplete CellChat candidates are not the published final results.

- [Inputs and acquisition](INPUTS.md): 81 raw MEX checksums, five frozen downloads,
  model/resource hashes and missing-resource limitations.
- [Recorded environment and parameters](ENVIRONMENT_AND_PARAMETERS.md): actual
  QC thresholds, priors, seeds, package versions and scope limitations.
- [Reproduction commands](REPRODUCE.md): portable raw-stage materialization and
  executable frozen-RNA scoring/CellChat replay.
- [Source execution map](SOURCE_EXECUTION_MAP.json): original executed file SHA256
  versus the published, path-redacted bytes. Path-only transformations preserve
  scientific code; intermediate historical source revisions were not all retained.
- [Dated evidence](evidence): original numeric conclusions and original hash
  bindings, with private roots redacted. Redacted receipt bytes have new hashes;
  embedded historical hashes are not signatures of the redacted copy.

The original DB cache is a validated subset of the full RNA CSC index. Every
extracted column and all axes were checked against that index; direct RNA probes
covered every library and all nine restored CellChat aliases. It is not a new
normalization. Sparse triMean aggregation preserves the installed CellChat model:
the 17,536-record, 100-permutation pilot had bitwise-identical probability/P-value
arrays, and the full 528,739-record aggregation matched the observed grouping and
first two seeded permutations (34,416 values). These bounded checks do not claim
that a separate dense full-cohort reference network completed.

The new replay entrypoint has been checked for parsing, source materialization,
command configuration and failure on wrong inputs. A clean full raw-to-final run
and an independent clean installation were not repeated for this source release.
Exact external KEGG/model resources must be supplied. No project licence or
permission for third-party data is assigned by publishing this package; existing
root licensing notices and embedded source notices continue to apply.

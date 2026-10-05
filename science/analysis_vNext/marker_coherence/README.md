# Fine-label full-RNA marker coherence audit

The actual October 5 run covered all 937,825 original integrated T/NK cells and 253,635 B/plasma cells (63 fine labels). All cells were exact-joined to the effective released raw-RNA H5. Counts were streamed from the 46 objects containing integrated members; all 48 object schemas/gene axes were checked. The AIH and RA GSE159117 objects have no integrated members and were not assigned invented fine labels.

The panel contains 40 genes, declared before the raw-count scan. At most 2048 cells and those 40 columns are made dense. Normalization uses the sum of the complete retained raw-RNA gene library, then log1p(count/library×10000); the old dense HVG data layer is never used. Source/specimen groups retain actual counts, detected fractions, means, missing-feature NA and clinical identity status. Specimens are not silently treated as independent donors.

The entrypoints are portable and require explicit local inputs. No actual clinical/cell-level tables are included in the source-code release.

```bash
python audit_full_RNA_marker_coherence.py \
  --input-tables /local/released_science/tables \
  --registry /local/dataset_assets.json \
  --data-h5-root /local/Data.H5 \
  --reader ../integration/released_rna.py \
  --crosswalk /local/legacy_object_to_H5_exact_crosswalk.tsv \
  --output /local/new_marker_evidence

python build_marker_review_evidence.py --input-dir /local/new_marker_evidence

python audit_targeted_actual_joint_markers.py \
  --input-tables /local/released_science/tables \
  --registry /local/dataset_assets.json \
  --data-h5-root /local/Data.H5 \
  --reader ../integration/released_rna.py \
  --crosswalk /local/legacy_object_to_H5_exact_crosswalk.tsv \
  --marker-audit-dir /local/new_marker_evidence \
  --output /local/new_actual_joint_evidence
```

The input directory requires the released-cell and original T/NK/B-plasma integration metadata tables plus `composition_specimen_denominators_and_eligibility.tsv`. The registry must point to effective available raw H5 paths. The crosswalk provides explicit original-object/H5 identity; basename or first-match guessing is forbidden.

For relative registry paths, `--data-h5-root` is required and is joined to the path exactly. Parent traversal, symlink escape and absent files are rejected. Absolute input paths remain explicit. Source metadata keeps `Sample` when present, verifies equal `Local_sample_label` if both are supplied, and accepts Local-only input by explicit rename. Missing clinical donor fields remain unresolved rather than converting specimens to donors.

Review thresholds apply to groups with at least 20 cells and are uncalibrated review screens, not accuracy criteria, disease tests or relabeling rules. TRDC alone does not identify gamma-delta T cells: final rule version 2 requires a CD3D/TRD joint-detection lower bound, derived from gene marginals. These are Frechet bounds, not directly measured cell coexpression, TCR identity or proof against doublets. Production labels and all raw assets remain unchanged.

The targeted third entrypoint measures **actual same-cell co-detection**, separately from those bounds, in all original GC/plasmablast/NKT/ILC/MAIT labels and source fine groups selected by review. Its additional 30-gene panel reads CD3D/E/G together with TRDC/TRGC1/TRGC2, and plasma joint programs; direct CD3D/TRD joint fractions can be checked against the original marginal bounds. It explicitly records the physical `JCHAIN` or legacy `IGJ` feature used, retaining both original coverage entries and never adding the two genes together. Unavailable required components make a program NA. Measured RNA programs remain annotation-review evidence and do not establish invariant TCR identity, calibrated accuracy or exclude doublets.

Frozen canonical-symbol requests keep unavailable genes NA. `JCHAIN` can be absent while legacy `IGJ` is present; separate alias availability is reported rather than silently changing the frozen panel. Limited marker coverage, dropout and absent TCR/protein validation prevent claims that every fine label has been independently authenticated. The full source/specimen evidence and specialist tables belong in the local scientific review package, not the public code repository.

All integrated cells must match the released metadata on the separate `OldCells`, `Sample` and `Dataset` fields. Both scans use validated left joins and reject every unmatched cell or specimen group; they never shrink the expected scope after an incomplete join. Empty or duplicated/conflicting crosswalk identities are rejected.

The reproducible anonymous check creates three synthetic cells and runs all three CLIs, including negative tests for missing identities; it reads no real biological tables:

```bash
python test_marker_portable_synthetic.py \
  --marker-code-dir . \
  --reader ../integration/released_rna.py \
  --output /local/new_synthetic_marker_check
```

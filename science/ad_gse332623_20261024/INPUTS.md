# Inputs, downloads and resource acquisition

The source study is [GEO GSE332623](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE332623),
associated with [DOI 10.1186/s12967-026-08288-7](https://doi.org/10.1186/s12967-026-08288-7).
Use the GEO supplementary download `GSE332623_RAW.tar`, extract its 27
`GSE332623_TotalseqB_NEA*_ADT*` matrix/features/barcodes triplets, and compare all
81 compressed members against [raw_inputs.json](resources/raw_inputs.json).
RNA and ADT rows coexist in each MEX matrix and must be separated using the third
feature-table column. Repeated barcode strings across libraries are not sufficient
evidence for duplicate cells; retain the original library prefixes.

```sh
python reproduce.py verify-inputs --raw-input /your/extracted/GSE332623
```

[capture_manifest.tsv](resources/capture_manifest.tsv) records 27 technical
libraries, public GEO/SRA/BioSample identifiers, original corrected-input hashes
and explicit unresolved clinical status. These are not 27 independently confirmed
donors or physical captures. [SampleInventory.GSE332623.tsv](resources/SampleInventory.GSE332623.tsv)
contains only the three columns needed by the reader (54 public sample rows).

The portal's [download API](https://scaid01.com/api/h5-file/?dataset_record_id=176)
provides five current immutable H5AD files. The observed public links and formats
are recorded in [public_downloads.json](resources/public_downloads.json). Download
the RNA entry to replay downstream kernels and verify its complete SHA256:

```sh
python reproduce.py verify-inputs --rna /your/downloads/Atopic_Dermatitis_GSE332623_PBMC_RNA.h5ad
```

| File method | SHA256 |
|---|---|
| RNA | `b9aa963123cb8f010d3e092e03507b70d02170afd93d32a6cdea2efe2a2c09ec` |
| ADT | `f23a4fb3b3e1812b478b80f635d017fb4a9e111686f000fccfc09fa683dfa7c4` |
| AUCell | `a87c7b41532777232d1584bed1b33515c3ff0f930ef7206254f73dbc2d70083e` |
| UCell | `8096fe3faa572cab30b9832e5d83e93e246ee7301ac431666d33f2781d69ddeb` |
| singscore | `0ed4882734f913748e5be922291ad4c5927928041cea4403ec94c9542ff3d01b` |

RNA `X` is source-count-derived log1p(total 10,000); `layers['counts']` holds
integer source counts aggregated by original RNA symbol. ADT is processed
separately with CLR and retains integer counts. Score files contain scores,
not counts, and keep matched RNA metadata and embeddings. All 581,725 records
remain in these five exports, including records excluded from CellChat modeling.
Study data reuse terms apply separately from this software snapshot.

External resources are not bundled:

| Resource | Recorded identity / acquisition |
|---|---|
| KEGG gene sets | Historical `20240526_kegg_hsa_list.rds`, 359 sets, 86,304 bytes; SHA256 `ccb053e043dc678796520c6e7c0a36ad8f09bd4dd8673e582dcd309526e608de`. Obtain the exact resource from the study maintainers subject to its terms. The inspected run did not preserve a public byte-identical download URL. A current KEGG refresh is a different analysis. |
| CellTypist High32 | `Immune_All_High.pkl`; SHA256 `a715fec36c2c421f7c4e31cf4cb4bea883eedab7fd7b20a4b76f091c22660448`. Obtain from the official CellTypist model repository and verify bytes. |
| CellTypist Low98 | `Immune_All_Low.pkl`; SHA256 `290874d35dac039d4c9218c343fde4aac1077709b72a331ce7266f6828c36502`. Obtain from the official CellTypist model repository and verify bytes. |
| CellChat DB | `CellChatDB.human` from the recorded CellChat 2.2.0.9001 runtime; saved RDS SHA256 `45536ae5df0ff6670a12ac97b3610a7455ef4ccdf8667346435131153a40835e`. The historical package commit was not recorded; version text alone does not guarantee identical DB bytes. |
| Human annotation | `org.Hs.eg.db` 3.18.0 for stable-ID matching and optional GO enrichment. Accepted aggregate mappings are included; this is not a live HGNC reannotation. |

Exact model and KEGG bytes are prerequisites for the corresponding stages.
Source publication does not remove these resource-acquisition limitations.

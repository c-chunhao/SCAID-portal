# Local module registration examples

`marker_manifest.example.json` preserves the field structure and reviewed
condition/accession/tissue bindings from a catalogue snapshot. Its numeric IDs
(including 145) refer to that local catalogue only. It contains no marker TSV
values and is not automatically imported. A newly migrated SQLite database is
empty and will normally assign different IDs. Use
`science/new_release/prepare_marker_assets.py` with a reviewed snapshot of your
own database; never bind by ID alone. The API validates both ID and scientific
metadata before serving a table.

The populated gene-cache manifest is intentionally excluded. The empty
`backend/tisch_api/data/lazy_gene_plots.example.json` documents the registry
container. Build a real registry locally with the shared index builder and set
`SCAID_GENE_CACHE_MANIFEST` to its path. Generated entries include local absolute
paths, exact source/index hashes, mtime/size, cells/genes/UMAP axis hashes and
filename mappings. The source package includes no biological matrix or PNG
cache. Keep the writable cache on an explicitly configured data volume.

The main asset map resolves relative figure paths against `SCAID_FIGURE_ROOT`
and relative download paths against the first `SCAID_H5_DOWNLOAD_ROOTS` root.
Any `h5_file_metadata` entries must use the same relative filename as `h5_paths`.
After moving paths, rebuild validated local registrations and the lazy index
registry; changing a root does not rewrite existing database rows.

Final stage2 deployment and anonymous/browser/download QA are maintained in the
release process. These examples do not assert that a local or production
installation has passed that acceptance.

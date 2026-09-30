# Public catalogue registration and clean installation

The public API is read-only. Use a separate maintenance database identity for
migrations and imports; do not grant write access to the web service identity.

## Installation

Install `requirements.txt`, configure Django settings, then run:

```sh
python manage.py migrate --noinput
python manage.py makemigrations --check --dry-run
```

Migration 0006 fills the previously omitted model state for the catalogue,
H5 registry and disease/site relations. Absent tables are created. Existing
complete tables are preserved; equivalent indexes are reused. Missing PDF
columns are added only on incomplete installations. The corrected 0005 also
creates `relative_path` before indexing it on a clean database; databases that
already applied 0005 retain their existing column.

Back up the database before a production migration. The 0006 reverse database
operation is deliberately a no-op: it cannot drop manually provisioned tables
or delete scientific data. Reverting application source does not require a
reverse schema migration.

## Reviewed public manifest

`import_public_manifest` registers existing prepared PNG/JPEG plots and HDF5
assets, along with aggregate catalogue metadata and gene/pathway vocabulary.
It does not process raw sequencing data, validate biological conclusions,
recompute QC, or import patient-level metadata. It rejects unknown catalogue
fields. `sample_size` is a count of sample labels, not an authenticated count
of biological donors.

```json
{
  "schema_version": 1,
  "release_id": "reviewed-public-release-YYYYMMDD",
  "source_url": "https://example.org/public-source",
  "datasets": [{
    "catalog": {
      "disease_name": "Reviewed condition",
      "abbreviation": "CODE",
      "dataset_source": "GEO",
      "dataset_id": "GSE...",
      "tissue": "PBMC",
      "sample_size": 0,
      "cell_count": 0,
      "batch": "YYYYMMDD"
    },
    "figures": [{"path": "CODE/GSE.../PBMC/GeneUmap/Umap_CD3D.png"}],
    "h5_files": [{"path": "RNA/CODE_GSE..._PBMC.h5", "method": "RNA"}]
  }]
}
```

First inspect the complete dry-run report:

```sh
python manage.py import_public_manifest public-manifest.json \
  --figure-root /absolute/reviewed/figures \
  --h5-root /absolute/reviewed/h5 \
  --report import-dry-run.json
```

Apply the same reviewed file explicitly:

```sh
python manage.py import_public_manifest public-manifest.json \
  --figure-root /absolute/reviewed/figures \
  --h5-root /absolute/reviewed/h5 \
  --apply --report import-applied.json
```

Default execution does not write to the database. Assets must exist, resolve
inside the supplied roots and have valid image/HDF5 headers. HDF5 header checks
are not a full contents or scientific-QC validation. Imports use natural keys,
reject duplicate/ambiguous identities and apply every database change within a
single transaction. A repeat import reports unchanged rows rather than adding
duplicates. Reports contain the manifest SHA-256 and public aggregate identity,
not server paths or patient identifiers. Save the input manifest and report with
the release. Do not run concurrent imports for the same new identities: the
legacy tables lack unique constraints on several natural keys.

New assets also need to be inside the runtime figure and H5 download allowlists.
A new dataset must have its reviewed asset relationship added to
`tisch_api/data/dataset_assets.json` to receive `asset_status: linked`; until
then it remains `unverified`. Use the same reviewed accession, condition and
tissue; directory spelling must not silently change the scientific diagnosis.

## Curated conditions and historical files

Public catalogue/tree `abbreviation` uses the manifest's condition code.
E-MTAB-8207 is `PsA`; `AS` selects the true ankylosing-spondylitis cohort.
GSE198616 is systemic `BD`; `BD(Uveitis)` is the separate aqueous-humour cohort.
`SV` remains a documented legacy search alias for the GSE198616 cohort. Original
`AS/...` and `SV/...` filenames and stored catalogue values are retained for
asset compatibility; no scientific objects are rewritten by this change.

## Figure API

Figure JSON uses explicit public fields (`id`, `name`, `relative_path`,
`page_number`, `image`, dimensions, `is_pdf`, `thumb_url`, and curated
`condition_abbreviation`, `disease_label`, `dataset_label`, `tissue_label`). It omits absolute
source `path` and internal timestamps. `image` is a same-origin stable endpoint:
`/api/pdf-images/{id}/original/`; clients must preserve this API URL without
prepending an old media base path. Originals are served only after a registered
path is checked against configured roots. Optional legacy rendered-page roots
must be configured explicitly using `SCAID_EXTRA_FIGURE_ROOTS`.

PNG previews use per-target process locks and independently named temporary
files, then atomic replacement. A cold-cache stampede produces one complete
JPEG rather than shared-temporary-file failures. Lock files stay in place to
preserve inode locking. The original figure's data and pixels are not changed.

## Synthetic fixture

Create the packaged demonstration in a temporary directory. It contains a
plain image, a real small HDF5 container and aggregate zero counts; no clinical
observations or patient identifiers are included. Never import this demonstration
into the scientific production catalogue.

```sh
python docs/create_demo_manifest.py /tmp/scaid-synthetic-demo
python manage.py import_public_manifest /tmp/scaid-synthetic-demo/public-manifest.json \
  --figure-root /tmp/scaid-synthetic-demo/figures \
  --h5-root /tmp/scaid-synthetic-demo/h5
```

The command above is a dry-run. To test `--apply`, configure a separate empty
database with its own Django settings, migrate it first, then inspect all
resulting catalogue, plot and download API routes.

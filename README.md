# SCAID source snapshot

**Version 0.1.0-rc3 · 30 September 2026**

This repository contains the current Vue/Vite frontend and Django REST backend
of the [SCAID portal](https://scaid01.com/), plus a reviewed public asset mapping,
aggregate catalogue snapshot and a synthetic fixture for local installation.
The portal provides catalogue browsing, precomputed gene/pathway/CellChat views
and registered HDF5 downloads. Its public API is read-only.

The frozen aggregate catalogue included here describes **48 disease–tissue
objects, 2,572,480 cell records, 625 sample labels and 24 condition labels**.
Sample labels are not verified independent donors. The condition labels include
subtypes and comparator conditions; they are not 24 equivalent diagnoses.
The 2024–2026 dataset collection is separate ongoing work and is not counted as
an already integrated release.

This is a source release candidate. It does not establish peer-reviewed publication,
a complete raw-data-to-paper reproduction or a separate archival DOI.
Processed matrices, sequencing reads, cell/donor tables, clinical spreadsheets,
database dumps, production credentials and confidential review correspondence
are excluded. The scientific atlas images used by the current frontend are
included as aggregate figure assets; the complete dataset plot/HDF5 library is
not included. See [LICENSE_STATUS.md](LICENSE_STATUS.md) for licensing status.

## Local backend and synthetic fixture

The checked runtime is Python 3.13.11. `backend/requirements.txt` records the
observed web runtime versions. Use a fresh environment and local paths:

```sh
cd backend
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
export SCAID_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_urlsafe(64))')"
export SCAID_SQLITE_PATH=/tmp/scaid-demo.sqlite3
export SCAID_FIGURE_ROOT=/tmp/scaid-synthetic-demo/figures
export SCAID_H5_DOWNLOAD_ROOTS=/tmp/scaid-synthetic-demo/h5
export SCAID_MEDIA_ROOT=/tmp/scaid-synthetic-demo/media
export SCAID_THUMBNAIL_ROOT=/tmp/scaid-synthetic-demo/thumbs
python manage.py check
python manage.py migrate --noinput
python manage.py makemigrations --check --dry-run
python docs/create_demo_manifest.py /tmp/scaid-synthetic-demo
python manage.py import_public_manifest /tmp/scaid-synthetic-demo/public-manifest.json \
  --figure-root /tmp/scaid-synthetic-demo/figures \
  --h5-root /tmp/scaid-synthetic-demo/h5
```

The last command is a **dry-run**. To populate the separate demonstration DB,
repeat it with `--apply`. The fixture contains no patients or biological dataset;
its zero aggregate counts and small HDF5 illustration are explicitly synthetic.
Do not import it into the scientific production catalogue.

```sh
python manage.py test tisch_api.test_api_fixes tisch_api.test_thumbnails \
  tisch_api.test_catalog tisch_api.test_openapi tisch_api.test_tree_query_count \
  tisch_api.test_release_fixes --settings=tisch_api.test_settings
python manage.py runserver 127.0.0.1:18080
```

The test configuration uses an isolated database and a nonproduction test key.
The synthetic demo allows catalogue, image, preview and download routes to be
checked without redistributing patient data. See
[backend/docs/PUBLIC_IMPORT.md](backend/docs/PUBLIC_IMPORT.md) for the input
schema, natural-key matching, transactions, audit reports and limitations.
A clean migrated database is otherwise empty. The aggregate release snapshot
is documentation, not a patient table or an automatic production import.

## Frontend

The checked local Node version is 24.13.0. Install the locked dependencies:

```sh
cd frontend
npm ci
npm run lint
npm test
npm run build
npm run dev -- --host 127.0.0.1
```

Vite proxies `/api` to `http://127.0.0.1:18080`; override with `SCAID_API_TARGET`.
Original figures use `/api/pdf-images/{id}/original/`, and PNG previews use
`/api/pdf-images/{id}/thumb/`. Scientific figure captions prefer the API's
curated condition/dataset/tissue labels, rather than historical asset folders.
Browser test instructions are in [frontend/README.md](frontend/README.md).

## Deploying reviewed data

Production settings are configured entirely outside Git through environment
variables. `.env.example` contains placeholders only and is not loaded
automatically. Set `SCAID_DEPLOYMENT=production`, explicit allowed hosts and a new
secret. Install `backend/requirements-production.txt` for MySQL. Use a
SELECT-only database user for the public web service and a separate maintenance
identity for migrations/imports. Bind the application/database to local
interfaces and configure a reverse proxy; Django's runserver is for local checks.

The packaged asset manifest uses relative paths. `SCAID_FIGURE_ROOT` and the first
`SCAID_H5_DOWNLOAD_ROOTS` root resolve them; optional rendered-page roots are
explicitly configured with `SCAID_EXTRA_FIGURE_ROOTS`. Register actual reviewed
assets at these locations before expecting real plots/downloads. Moving roots
does not automatically rewrite existing database registrations. New catalogue
items remain unverified until their reviewed relationship is added to the
asset manifest. Current PsA, true AS, systemic BD and BD(Uveitis) scopes are
separate; legacy AS/SV filenames remain compatible.

The hardened HTTPS proxy mode requires a proxy that strips untrusted forwarding
headers before setting `X-Forwarded-Proto`. Enable `SCAID_TRUST_PROXY_SSL=true`
only for that deployment topology. H5/preview nginx acceleration is optional and
requires matching `internal` locations; without it Django streams validated
registered files.

## Source provenance and scope

Current application files are copied from the running project's working source,
with portable release-only changes to settings, path resolution and packaging.
The original teaching/template Git history is excluded. File hashes are in
`SOURCE_SHA256SUMS.txt`; `SOURCE_PROVENANCE.json` describes copy/adaptation origins
without exposing the server's private paths. Source availability does not verify
scientific claims, donor identities or every original study's redistribution
terms. The software and dataset licences are separate.

## Scientific analysis scripts

See [science/README.md](science/README.md) for source-reference concordance,
coverage-aware descriptive composition and verified-individual SLE analysis.
The [reviewed analysis modules](science/analysis_vNext/README.md) also provide
strict released-RNA readers, integration and export guards, source checks,
and scientific evidence gates. Their inputs must be supplied locally.

The [fine-label marker audit](science/analysis_vNext/marker_coherence/README.md)
separates marginal detection bounds from measured same-cell marker programs.
The [frozen PSO replay](science/analysis_vNext/qc/PSO_REPLAY.md) documents the
six preverified capture inputs, pinned doublet-calling runtime and fail-closed
classifier checks. These modules do not establish independent donor identities,
calibrated annotation accuracy or complete historical preprocessing records.
Source cell/donor inputs remain subject to repository access and reuse conditions.
The committed September 30 composition outputs are historical; reproduce the
current coverage-aware outputs with the documented CLI.

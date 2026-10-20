# SCAID source snapshot

**Version 0.2.0-rc1 · 20 October 2026**

This repository contains the current Vue/Vite frontend and Django REST backend
of the [SCAID portal](https://scaid01.com/), plus a reviewed public asset mapping,
aggregate catalogue snapshot and a synthetic fixture for local installation.
The portal provides catalogue browsing, precomputed gene/pathway/CellChat views
and registered HDF5 downloads. Its public API is read-only.

The aggregate catalogue now documents **79 source–condition–tissue objects,
4,590,932 cell records and 1,095 entry-specific sample labels under 36 condition
labels**. The 48-object analytical collection contains 2,572,480 records and
625 labels. The 31 additional browser entries contain 28 expression/score
objects and three overview-only objects, without new CellChat or H5 downloads.
Sample labels are not independent donor counts. Catalogue records and module
flags were checked against the anonymous API on 20 October 2026. The atlas
inputs and catalogue scope must not be added together as independent cells.

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
  tisch_api.test_release_fixes tisch_api.test_module_availability --settings=tisch_api.test_settings
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

## October production-source synchronization

This release refreshes the current frontend, API module-availability logic,
case-sensitive gene lookup, current 79-entry asset catalogue and portable analysis
modules. Business logic comes from the running deployment's working source;
settings and asset roots are adapted for local installation. The source package
does not include production database registrations, figure libraries, H5 matrices,
clinical records, secrets or deployment certificates. Local demo and regression
tests are reproducible without them. Runtime mirroring does not establish a full
raw-data-to-paper reproduction or independent scientific validation.

Set `SCAID_PUBLISHED_FIGURES_BASE` for the relative prepared gene/score family
layout in `published_figure_roots.json`. Import reviewed file registrations using
the documented manifest schema. Absent data families return unavailable results;
installing source alone does not redistribute the production data library.

`LICENSE_STATUS.md` remains authoritative: no new project licence is assigned
by this release. Third-party source retains its original terms and provenance.
A GitHub release is a versioned source reference, not a permanent archive DOI.

# SCAID source snapshot

**Version 0.3.1-ad20261024 · 24 October 2026**

This revision synchronizes the current custom application files and adds the
[actual GSE332623 analysis source, environment, parameters and reproduction guide](science/ad_gse332623_20261024/README.md).
The current catalogue has **80 source–condition–tissue objects, 5,172,657 cell
records and 1,122 entry-specific technical/sample labels**. The added mixed
AD/healthy-reference cohort retains 581,725 records in 27 technical libraries;
donor and clinical groups are unresolved and clinical comparisons are disabled.
The two original integrated atlases remain unchanged. Current API and UI source
concordance is recorded in [the source map](science/ad_gse332623_20261024/SOURCE_EXECUTION_MAP.json).

The default source retains all earlier backend modules. The current deployed
Dataset UI has figure tabs and download entries; the earlier candidate's
MarkerTable UI component is retained as source but has no current Dataset tab.
Its earlier candidate validation statements below describe their original dates
and scope. The 24 October source update does not newly reproduce those historical
checks or a complete clean raw-data-to-final scientific analysis.

## Earlier 21 October module release context

This repository contains the current Vue/Vite frontend and Django REST backend
of the [SCAID portal](https://scaid01.com/), plus a reviewed public asset mapping,
aggregate catalogue snapshot and a synthetic fixture for local installation.
The portal provides catalogue browsing, precomputed gene/pathway/CellChat views
and registered HDF5 downloads. Its public API is read-only.

The aggregate catalogue documents **79 source–condition–tissue objects,
4,590,932 cell records and 1,095 entry-specific sample labels under 36 condition
labels**. The existing 48-object collection has 2,572,480 records and 625 labels.
This candidate adds module support for the 31 additional browser objects:
CellChat, cell-type-versus-rest marker tables, 124 H5AD exports and exact-object
lazy gene maps for three previously missing Gene modules. The two existing
atlases remain unchanged. Sample labels are not authenticated independent donors,
and catalogue/atlas scopes must not be added as independent cells.

中文索引：[新增功能与科学脚本](science/new_release/README.md)、
[本地资源注册示例](backend/docs/new_release_assets/README.md)、
[许可状态](LICENSE_STATUS.md)。新增功能已于2026年10月21日部署并完成下述范围的验收。本包仍为源码候选；数据矩阵、
个人元数据、生产密钥、绝对路径运行清单和图片缓存不随源码分发。

The [new-module source notes](science/new_release/README.md) distinguish 30 RNA
objects retaining integer source counts from the one processed-only RNA object;
the optional 24-GSM raw-count subset is additional to the main 124 exports.
The PBMC pilot (7 objects, 5 studies, 14,000 cells) is descriptive and does not
replace either atlas or establish donor-level disease comparisons. The modules were deployed and verified on 21 October 2026 (Asia/Shanghai).
Local verification passed 81 backend tests, 19 frontend tests, lint and build.
The two production acceptance stages recorded 518 API checks (320 + 198);
stage two made 374 HTTPS requests. All 125 download files passed HEAD and Range
checks, and four small files passed complete-download SHA256 verification.
These checks do not claim a network rehash of the complete approximately 40 GB
download library. Desktop/mobile scenarios and six Gene/KEGG interaction
scenarios also passed. Historical
rc2 installation and website checks below retain their original dates/scope.

This is a source release candidate. It does not establish peer-reviewed publication,
a complete raw-data-to-paper reproduction or a separate archival DOI.
Processed matrices, sequencing reads, cell/donor tables, clinical spreadsheets,
database dumps, production credentials and confidential review correspondence
are excluded. The scientific atlas images used by the current frontend are
included as aggregate figure assets; the complete dataset plot/HDF5 library is
not included. See [LICENSE_STATUS.md](LICENSE_STATUS.md) for licensing status.

## Local backend and synthetic fixture

The production web runtime was Python 3.13.11. A separate clean local web
installation and synthetic SQLite/API demo also passed with Python 3.12.12.
The exact local installation is recorded in
`backend/requirements-tested-python312.txt`; the original observed web versions
remain in `backend/requirements.txt`. Use a fresh environment and local paths:

```sh
cd backend
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-tested-python312.txt
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
  tisch_api.test_release_fixes tisch_api.test_module_availability \
  tisch_api.test_markers tisch_api.test_lazy_gene_plots \
  tisch_api.test_new_release_modules tisch_api.test_dynamic_gene_thumbnails \
  --settings=tisch_api.test_settings
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

Install the current scientific packages in the separate
[Python 3.12 environment](science/current_revision/INSTALLATION.md); the
Python 3.13 runtime above refers to the web service.

The [October scientific follow-up](science/current_revision/README.md) adds
final-representation integration diagnostics, 63-label marker-evidence review
and source-recorded B-state/IFN analysis. Each package includes dependencies,
input schemas and an explicitly synthetic installation check; B-state outputs
are aggregate only. The real clinical/count inputs are supplied separately.

## October production-source synchronization

The preceding rc2 release refreshed the frontend, API module-availability logic,
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

The [additional source-pool sensitivity](science/current_revision/source_pool_sensitivity/README.md)
uses existing individual summaries, CR1 covariance and pool fixed effects, with
aggregate results and a fully synthetic end-to-end example. The versioned
verification evidence distinguishes the original live-site/rc1 build checks
from rc2 independent installation and local proxy checks. The v0.2.0-rc2 source release did not redeploy the website.


## Additional module candidate and verification state

The v0.3.0-rc1 candidate adds source for exact-object CellChat overview families,
marker table pagination/filtering/TSV downloads, H5AD content metadata and
source-bound lazy gene rendering. Install the web dependencies from
`backend/requirements.txt` (or the Python 3.12 candidate lock); matplotlib and
NumPy are required for lazy rendering. The preceding rc2 Python 3.12 lock has
been extended for the new renderer and is not a fresh-installation receipt.
Scientific dependencies are separate in `science/new_release/requirements.txt`.

Set `SCAID_MARKERS_ROOT` and `SCAID_GENE_CACHE_MANIFEST` to locally prepared,
reviewed resources. The source's gene registry is absent by default, so it does
not attempt to load server-specific biological indices. Marker example IDs are
local catalogue bindings, not portable database primary keys. See the
[registration notes](backend/docs/new_release_assets/README.md).

The 31 additional objects are deployed with the module and download availability
recorded in this candidate. Production acceptance completed on 21 October 2026
within the API/browser/download scope summarized above and in
[the verification notes](docs/CANDIDATE_CHECKS.md). Release evidence records the
frozen source tree separately from the production acceptance reports.
The original licence and unresolved frontend-template permission notices remain
unchanged; preparing this candidate grants no additional reuse rights.

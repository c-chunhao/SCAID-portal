# SCAID web client

Vue 3 single-page application for the Single-Cell AutoImmune Disease database
(https://scaid01.com/). The client is read-only: it renders precomputed
figures and catalogue records served by the Django REST API under `/api/` and
stable original figure endpoints under `/api/pdf-images/{id}/original/`. No user accounts exist; `/login` is an
explanatory page only.

## Stack

- Vue 3.5, Vue Router 4 (code-split routes), Pinia
- Element Plus 2.8 (auto-imported components), Sass
- Axios with a same-origin `/api` base URL
- Vite 5 build; Node 18+ recommended

## Pages

| Route | Purpose |
|---|---|
| `/Home` | Disease browser and free-text search (`/api/by-site/`, `/api/cell-data-all/?name=`) |
| `/Dataset` | Filter disease–tissue data objects, view Overview / Marker / GO / KEGG / CellChat figures, download H5 files |
| `/Gene`, `/KEGG` | Gene expression and pathway-score UMAP maps (`/api/pdf-images/`) |
| `/Information` | Step-by-step guide, H5 schema, abbreviation glossary (generated from the catalogue), API reference |
| `/TeamInformation` | Atlas figure gallery |
| Other paths | Page-not-found recovery with the main navigation retained |

## Share or restore a query

Dataset and exploration searches update the address bar when a query is
submitted, a record/analysis is selected or a result page changes. **Copy query
link** copies the same view; if clipboard access is unavailable, it exposes a
selectable URL. Browser Back and direct links restore the submitted view.
Changing a draft filter does not silently run a search.

| Page | Query parameters |
|---|---|
| `/Dataset` | `id` (disease abbreviation, retained for existing Home links), `dataset`, `tissue`, `search`, `record` (catalogue record ID), `analysis` (figure category), `page` |
| `/Gene` | `gene`, `disease`, `dataset`, `page`, `size` |
| `/KEGG` | `kegg`, `disease`, `dataset`, `method` (`AUCell`, `UCell`, `singscore`), `page`, `size` |

Empty filters, default page 1, default page size 10 and Dataset's default
`Overview` analysis are omitted. Valid exploration sizes are 10, 20, 50 and
100. For example:

```text
/Dataset?record=130&analysis=LRcircle&page=2
/Gene?gene=TRIM21&disease=CI&dataset=GSE121380
/KEGG?kegg=hsa00010&method=AUCell&page=2&size=20
```

Record IDs and plot availability belong to the current API catalogue. A link
whose record is unavailable for its filters displays a recovery message.

## Develop

```sh
npm install
SCAID_API_TARGET=http://127.0.0.1:18080 SCAID_MEDIA_TARGET=http://127.0.0.1:10209 npm run dev
```

The dev server proxies `/api` and `/system/media/` to the targets above. In
production Nginx serves the `dist/` build and reverse-proxies those two
prefixes to the API host, so the client needs no environment variables.

## Check and build

```sh
npm run lint      # ESLint, no auto-fix
npm test          # node:test: search requests/races, URL restoration and URL handling
npm run build     # outputs dist/
```

The Playwright scripts use fixture APIs while exercising the real running
client. They are not part of `npm test`:

```sh
SCAID_UI_URL=http://127.0.0.1:5173 node tests/dataset-browser-regression.mjs
SCAID_UI_URL=http://127.0.0.1:5173 SCAID_CHECK_OUTPUT=/tmp/scaid-ui-check node tests/frontend-fixes-regression.mjs
```

Install/provide `playwright-core` and Chromium, and set `PLAYWRIGHT_CORE_PATH`
and `CHROMIUM_PATH` when their locations differ from the script defaults.
Create `SCAID_CHECK_OUTPUT` before the second script; it receives screenshots
and a JSON result. Real API availability and catalogue contents need separate
deployment checks.

## Conventions

- All numbers shown to users come from the API (`/api/cell-data-all/`,
  `/api/cell-data/`); do not hard-code catalogue statistics in templates.
- Download links are only rendered when the API returns a same-origin
  `/api/h5-file/<id>/download/` path.
- Original figure links accept same-origin `/api/pdf-images/<id>/original/`
  responses; thumbnail URLs come from the API. Scientific figures are shown
  as supplied. Preview failures offer retry and the original-file link.
- Scientific wording follows the manuscript: cell records, sample labels and
  precomputed views; the site does not run analyses on request.
- Sample labels are source identifiers and may include pools/repeated
  observations. They must not be presented as verified independent donors.

<template>
  <div class="documentation-page">
    <aside class="sidebar" aria-label="Documentation navigation">
      <div class="sidebar-top">
        <p class="sidebar-heading">SCAID guide</p>
        <button type="button" class="contents-toggle" :aria-expanded="contentsOpen" aria-controls="guide-contents" @click="contentsOpen = !contentsOpen">{{ contentsOpen ? 'Close contents' : 'Contents' }}</button>
      </div>
      <nav id="guide-contents" :class="{ 'contents-open': contentsOpen }" aria-label="Guide contents">
        <div v-for="group in menuGroups" :key="group.name" class="menu-group">
          <p class="menu-group-title">{{ group.name }}</p>
          <ul class="submenu">
            <li v-for="item in group.children" :key="item.id">
              <button
                type="button"
                :class="{ active: selectedMenu === item.id }"
                :aria-current="selectedMenu === item.id ? 'location' : undefined"
                :aria-controls="item.id"
                @click="scrollToSection(item.id)"
              >{{ item.name }}</button>
            </li>
          </ul>
        </div>
      </nav>
    </aside>

    <article class="guide-content" aria-label="SCAID documentation">
      <section id="resource-overview" tabindex="-1" aria-labelledby="overview-title">
        <h1 id="overview-title">Explore SCAID, step by step</h1>
        <p class="review-date">Information &amp; usage · Checked 30 September 2026</p>
        <p class="lead">Find a dataset, inspect its precomputed figures, and download the data with the context needed to interpret it.</p>
        <p>SCAID brings together single-cell data from autoimmune and related immune-mediated conditions. The current catalog contains <strong>48 disease–tissue data objects</strong>, linked to <strong>34 repository identifiers</strong>, with <strong>24 reconciled condition labels</strong> and <strong>14 tissue labels</strong>. E-MTAB-8207 is catalogued as psoriatic arthritis and GSE198616 as Behçet disease; the older object identifiers are kept for traceability.</p>
        <div class="scope-summary" aria-label="Audited resource size">
          <div><strong>2,572,480</strong><span>cell records in the audited objects</span></div>
          <div><strong>625</strong><span>distinct sample labels</span></div>
        </div>
        <p>A sample label is not necessarily an independent donor, and a repository identifier is not necessarily an independent study. One study can contribute several tissues or disease groups. Counts and examples on this page describe the checked snapshot and may change with later releases.</p>
        <p>Use the standard HTTPS address <a href="https://scaid01.com/">scaid01.com</a>. No account is required for public browsing.</p>
        <p>The web application source and scientific analysis scripts are available in the <a href="https://github.com/c-chunhao/SCAID-portal" target="_blank" rel="noopener noreferrer">public GitHub repository</a>, with installation instructions and a synthetic demonstration dataset. Consult the repository for software licensing and original dataset reuse conditions.</p>
      </section>

      <section id="dataset-workflow" tabindex="-1" aria-labelledby="dataset-title">
        <h2 id="dataset-title">1. Start with a disease: IBD</h2>
        <ol>
          <li>On Home, enter <code>IBD</code> in <strong>Find a disease or dataset</strong>, then click <strong>Explore</strong> or press Enter.</li>
          <li>The Dataset page retains the query and returns <strong>10 dataset records</strong> in the checked snapshot. IBD is an umbrella text query covering the catalog labels CI, CD and UC; it is not a separate disease abbreviation option.</li>
          <li>Find <strong>CI · GSE121380 · colon</strong> and select <strong>View details</strong> in its dataset column. Check the selected disease, accession and tissue below the table. The <strong>Overview</strong> section contains four figures for this record.</li>
          <li>Select <strong>Marker</strong>, <strong>GO</strong>, <strong>KEGG</strong>, or a <strong>Cellchat</strong> submenu to view available results. Use figure pagination when a category has more than one page.</li>
          <li>To narrow your search, choose <strong>Disease abbreviation</strong>, then <strong>Dataset ID</strong>, then <strong>Tissue</strong>, and click <strong>Search</strong>. Use <strong>Reset</strong> to return to the full catalog.</li>
        </ol>
        <router-link class="example-link" :to="{ path: '/Dataset', query: { search: 'IBD' } }">Open the IBD example →</router-link>
        <p class="note">Selecting another dataset returns its details to Overview. Starting a new search clears the previous selection; click a result row to open its details. The catalog’s Sample Size field must not be assumed to count independent donors.</p>
      </section>

      <section id="gene-workflow" tabindex="-1" aria-labelledby="gene-title">
        <h2 id="gene-title">2. Search for a gene: CD3D</h2>
        <ol>
          <li>Open <router-link to="/Gene">Gene</router-link>, type <code>CD3D</code> into the Gene field, and <strong>select the matching suggestion</strong>. Typing alone does not complete the selection.</li>
          <li>Leave Disease (abbreviation) and Dataset ID empty, then click <strong>Search plots</strong>.</li>
          <li>The checked snapshot returns <strong>48 plots</strong>. Page 1 shows <strong>1–10</strong>; page 5 shows <strong>41–48</strong>. Use the page numbers, Next, or the page-size selector to explore all results.</li>
          <li>To restrict the results, select a disease abbreviation or accession and search again.</li>
        </ol>
        <p>For a broad gene suggestion search, use <strong>Load more</strong> when available. No matching plot means that a figure is unavailable for the selected query; it does not prove zero expression.</p>
      </section>

      <section id="pathway-workflow" tabindex="-1" aria-labelledby="pathway-title">
        <h2 id="pathway-title">3. Explore a pathway: hsa00010</h2>
        <ol>
          <li>Open <router-link to="/KEGG">KEGG</router-link>, type <code>hsa00010</code>, and select the matching suggestion.</li>
          <li>Choose <strong>AUCell</strong> under Scoring method. Leave disease and dataset filters empty, then click <strong>Search plots</strong>.</li>
          <li>The checked snapshot returns <strong>48 plots</strong>, with 10 initially shown. Navigate the result pages to view the rest.</li>
          <li>Choose <strong>UCell</strong> or <strong>singscore</strong> and search again to view that method’s output. Leaving Scoring method empty returns <strong>144 plots</strong> for this example: three methods across 48 objects.</li>
        </ol>
        <p>The KEGG page shows single-cell pathway-score figures. The Dataset page’s KEGG tab shows precomputed enrichment figures. These are different outputs, and score values from different methods or processing contexts do not automatically share a comparable scale.</p>
      </section>

      <section id="downloads" tabindex="-1" aria-labelledby="downloads-title">
        <h2 id="downloads-title">4. Download and inspect an H5 file</h2>
        <p>On Dataset, select a row and verify its disease, accession and tissue. Once the list has loaded, click <strong>Download (4)</strong> and choose a method-labelled entry.</p>
        <div class="table-scroll" role="region" aria-label="Download formats" tabindex="0">
          <table>
            <thead><tr><th scope="col">Category</th><th scope="col">Data product</th><th scope="col">Before reuse</th></tr></thead>
            <tbody>
              <tr><th scope="row">RNA</th><td>RNA assay and associated object information</td><td>Check expression transformation and count provenance; a layer name alone does not certify raw counts.</td></tr>
              <tr><th scope="row">AUCell</th><td>Precomputed pathway scores</td><td rowspan="3">Use the specified scoring method and keep its name with the downloaded file. These are not gene-count matrices.</td></tr>
              <tr><th scope="row">UCell</th><td>Precomputed pathway scores</td></tr>
              <tr><th scope="row">singscore</th><td>Precomputed pathway scores</td></tr>
            </tbody>
          </table>
        </div>
        <p>Method entries may share a basename, such as <code>IBD_GSE121380_colon.h5</code>. Save each in a separate method folder or add the method to its filename to avoid overwriting another download.</p>
        <h3>HDF5 is a container, not a universal single-cell schema</h3>
        <p>The inspected RNA example contains <code>assay</code>, <code>graphs</code>, <code>names_obs</code>, <code>names_var</code>, <code>obs</code>, <code>reductions</code>, <code>uns</code> and <code>var</code>. The inspected score files contain <code>matrix/cell</code>, <code>matrix/pathway</code>, <code>matrix/value</code>, <code>matrix/x</code> and <code>matrix/y</code>.</p>
        <p><strong>Do not pass these files directly to <code>scanpy.read_h5ad</code> or assume that renaming them to .h5ad converts them.</strong> Inspect the schema first and use a reader or conversion matched to that export format. A header-only inspection with Python is:</p>
        <pre><code>import h5py

with h5py.File("RNA/IBD_GSE121380_colon.h5", "r") as handle:
    print(list(handle.keys()))
    print(dict(handle.attrs))</code></pre>
        <p>Before loading expression values, verify matrix orientation, cell and feature identifiers, sparse indexing, and normalization. The example above inspects the file; it is not a complete reader or conversion routine.</p>
      </section>

      <section id="integration" tabindex="-1" aria-labelledby="integration-title">
        <h2 id="integration-title">5. Use the integrated atlases</h2>
        <p>Open <router-link to="/TeamInformation">Data Integration</router-link> to inspect the precomputed T/NK and B/Plasma atlas figures. These offer another view of annotated immune populations across the collection. The page does not re-integrate a selected cohort or run a new statistical comparison.</p>
        <p>Use cell-type labels, tissue and source information together when interpreting clusters. A well-mixed embedding does not establish that RNA expression is free of study or platform effects. Comparisons across diseases require a documented donor/sample design and appropriate controls.</p>
      </section>

      <section id="interpretation" tabindex="-1" aria-labelledby="interpretation-title">
        <h2 id="interpretation-title">6. Interpret the results with their context</h2>
        <p>SCAID displays previously generated figures. Browser searches do not recalculate differential expression, UMAP coordinates, pathway scores or CellChat networks.</p>
        <ul>
          <li><strong>Processing:</strong> retained command records describe LogNormalize with a scale factor of 10,000 in the 48 audited objects; 44 contain Harmony reductions. This does not establish uniform upstream QC or expression comparability across studies.</li>
          <li><strong>Metadata:</strong> sample labels, donors, tissue, treatment and source study are different variables. Missing metadata and unrecovered QC or doublet records should remain unavailable or unverified.</li>
          <li><strong>Annotation:</strong> final cell types reflect expert curation informed by available automated annotations. Agreement with those labels is not an independent accuracy benchmark.</li>
          <li><strong>Inference:</strong> marker and enrichment figures are descriptive. Predicted ligand–receptor communication does not by itself establish a causal interaction or a therapeutic target.</li>
          <li><strong>Replication:</strong> cells from one donor are not independent biological replicates. Formal disease comparisons need appropriate biological units, control groups, covariates and independent validation.</li>
        </ul>
      </section>

      <section id="troubleshooting" tabindex="-1" aria-labelledby="troubleshooting-title">
        <h2 id="troubleshooting-title">7. Troubleshoot and record a reproducible query</h2>
        <dl class="troubleshooting-list">
          <dt>No matching plots</dt>
          <dd>Select an actual Gene or KEGG suggestion, remove restrictive filters, and click Search again. A missing figure is not a negative biological result.</dd>
          <dt>Filters, figures or downloads could not be loaded</dt>
          <dd>This indicates a request failure. Use Retry where available, or search again. If the problem persists, record the exact visible error.</dd>
          <dt>No Download button or dataset details</dt>
          <dd>Select a result row first, then wait for its download list. A new search clears the previous selection.</dd>
          <dt>The legacy numeric address is blocked</dt>
          <dd>Use <a href="https://scaid01.com/">https://scaid01.com/</a>, which uses the standard HTTPS port.</dd>
        </dl>
        <p>For an issue report or an analysis record, include the date, page URL, query, disease, accession, tissue, category or scoring method, page number, downloaded method and filename, and any error message. Keep the original study accession and cite the source study when reusing its data.</p>
        <p>For help, email <a href="mailto:chunhao.bio@gmail.com">chunhao.bio@gmail.com</a> with those details.</p>
      </section>

      <section id="abbreviations" tabindex="-1" aria-labelledby="abbreviations-title">
        <h2 id="abbreviations-title">8. Disease abbreviations used in the catalogue</h2>
        <p>The Dataset filters and result tables use short codes. This list is generated from the live catalogue, so it always matches the current release. IBD is a text query that covers CI, CD and UC; it is not a separate code.</p>
        <p v-if="glossaryLoading" role="status">Loading abbreviations…</p>
        <el-alert v-else-if="glossaryError" :title="glossaryError" type="error" show-icon :closable="false">
          <el-button link type="primary" @click="loadGlossary">Retry</el-button>
        </el-alert>
        <div v-else class="table-scroll" role="region" aria-label="Disease abbreviations" tabindex="0">
          <table>
            <thead><tr><th scope="col">Code</th><th scope="col">Condition label</th><th scope="col">Data objects</th><th scope="col">Tissues</th></tr></thead>
            <tbody>
              <tr v-for="row in glossary" :key="row.code">
                <td><code>{{ row.code }}</code></td>
                <td>{{ row.labels.join('; ') }}</td>
                <td>{{ row.objects }}</td>
                <td>{{ row.tissues.join(', ') }}</td>
              </tr>
            </tbody>
          </table>
        </div>
      </section>

      <section id="api-access" tabindex="-1" aria-labelledby="api-title">
        <h2 id="api-title">9. Programmatic access</h2>
        <p>Every page reads the same read-only JSON API, which needs no key or account. The machine-readable description is at <a href="/api/schema/">/api/schema/</a> (OpenAPI 3.0) and the endpoint index at <a href="/api/">/api/</a>. Write methods are rejected with 405.</p>
        <div class="table-scroll" role="region" aria-label="API endpoints" tabindex="0">
          <table>
            <thead><tr><th scope="col">Endpoint</th><th scope="col">Returns</th><th scope="col">Typical query</th></tr></thead>
            <tbody>
              <tr><td><code>/api/cell-data-all/</code></td><td>Flat list of the 48 disease–tissue data objects</td><td><code>?name=IBD</code>, <code>?abbreviation=SLE</code>, <code>?dataset_id=GSE174188</code></td></tr>
              <tr><td><code>/api/cell-data/</code></td><td>Same catalogue as a disease → accession → tissue tree</td><td>no parameters</td></tr>
              <tr><td><code>/api/pdf-images/</code></td><td>Precomputed figures for one object (paginated)</td><td><code>?dataset_record_id=130&amp;pdf_name=Overview</code>, <code>?name=TRIM21&amp;GeneUmap=1</code></td></tr>
              <tr><td><code>/api/h5-file/</code></td><td>RNA and pathway-score H5 files for one object, with <code>download_url</code></td><td><code>?dataset_record_id=130</code></td></tr>
              <tr><td><code>/api/gene-kegg-enum/</code></td><td>Indexed gene symbols and KEGG pathway ids</td><td><code>?tp=gene&amp;name=TRIM</code></td></tr>
              <tr><td><code>/api/by-site/</code></td><td>Affected sites with their diseases (home page browser)</td><td>no parameters</td></tr>
            </tbody>
          </table>
        </div>
        <pre><code>curl -s "https://scaid01.com/api/cell-data-all/?name=IBD" | python -m json.tool
curl -s "https://scaid01.com/api/h5-file/?dataset_record_id=130"
curl -OJ "https://scaid01.com/api/h5-file/273/download/"   # RNA matrix, IBD_GSE121380_colon.h5</code></pre>
        <p>Responses describe the released objects only. Cell counts are cell records; sample sizes are sample labels, not verified donors. Cite the original study accession when reusing downloaded data.</p>
        <button type="button" class="back-to-top" @click="scrollToSection('resource-overview')">Back to the beginning ↑</button>
      </section>
    </article>
  </div>
</template>

<script setup>
defineOptions({ name: 'ScaidDocumentation' })
import { onMounted, ref } from 'vue'
import { getRequest } from '@/api/home'

const glossary = ref([])
const glossaryLoading = ref(false)
const glossaryError = ref('')

// Build the abbreviation table from the catalogue itself so it cannot drift from the data.
async function loadGlossary() {
  glossaryLoading.value = true
  glossaryError.value = ''
  try {
    const rows = await getRequest('/cell-data-all/', { format: 'json' })
    if (!Array.isArray(rows)) throw new Error('Invalid catalogue response')
    const byCode = new Map()
    for (const row of rows) {
      const code = row.abbreviation || '—'
      const entry = byCode.get(code) || { code, labels: new Set(), tissues: new Set(), objects: 0 }
      const label = row.disease_label || row.disease_full_name || row.disease_name
      if (label) entry.labels.add(label)
      if (row.tissue) entry.tissues.add(row.tissue_label || row.tissue)
      entry.objects += 1
      byCode.set(code, entry)
    }
    glossary.value = Array.from(byCode.values())
      .map(e => ({ code: e.code, labels: [...e.labels].sort(), tissues: [...e.tissues].sort(), objects: e.objects }))
      .sort((a, b) => a.code.localeCompare(b.code))
  } catch {
    glossaryError.value = 'Unable to load the abbreviation list. Please retry.'
  } finally {
    glossaryLoading.value = false
  }
}
onMounted(loadGlossary)

const menuGroups = [
  { name: 'Introduction', children: [{ id: 'resource-overview', name: 'Resource scope' }] },
  { name: 'Explore', children: [
    { id: 'dataset-workflow', name: 'Find an IBD dataset' },
    { id: 'gene-workflow', name: 'Search for CD3D' },
    { id: 'pathway-workflow', name: 'Explore hsa00010' },
  ] },
  { name: 'Reuse & help', children: [
    { id: 'downloads', name: 'Download and read H5 files' },
    { id: 'integration', name: 'Integrated atlases' },
    { id: 'interpretation', name: 'Interpretation and limits' },
    { id: 'troubleshooting', name: 'Troubleshooting' },
  ] },
  { name: 'Reference', children: [
    { id: 'abbreviations', name: 'Disease abbreviations' },
    { id: 'api-access', name: 'Programmatic access (API)' },
  ] },
]
const selectedMenu = ref('resource-overview')
const contentsOpen = ref(false)

const scrollToSection = (sectionId) => {
  const section = document.getElementById(sectionId)
  if (!section) return
  selectedMenu.value = sectionId
  contentsOpen.value = false
  section.focus({ preventScroll: true })
  section.scrollIntoView({
    behavior: window.matchMedia('(prefers-reduced-motion: reduce)').matches ? 'auto' : 'smooth',
    block: 'start',
  })
}
</script>

<style scoped>
.documentation-page { display: flex; max-width: 1320px; margin: 0 auto; background: #fff; color: var(--scaid-ink); }
.sidebar { position: sticky; top: 88px; align-self: flex-start; max-height: calc(100vh - 104px); width: 270px; flex-shrink: 0; overflow-y: auto; padding: 28px 16px; background: var(--scaid-soft); border-right: 1px solid var(--scaid-line); }
.sidebar-heading { margin: 0 8px 20px; font-size: 20px; font-weight: 650; color: var(--scaid-ink); }
.contents-toggle { display: none; }
.menu-group { margin-bottom: 18px; }
.menu-group-title { margin: 0 8px 6px; font-size: 14px; font-weight: 650; color: var(--scaid-muted); }
.submenu { list-style: none; padding: 0; margin: 0; }
.submenu button { display: block; width: 100%; text-align: left; border: 0; border-radius: 5px; background: transparent; color: #495057; padding: 10px 12px; line-height: 1.4; font-size: 15px; cursor: pointer; }
.submenu button:hover { color: var(--scaid-accent); background: #efe6df; }
.submenu button.active { color: var(--scaid-accent); font-weight: 650; background: #ede0d7; }
button:focus-visible, a:focus-visible, .table-scroll:focus-visible { outline: 2px solid var(--scaid-focus); outline-offset: 3px; }
.guide-content { min-width: 0; flex: 1; padding: 36px clamp(24px, 4vw, 56px) 60px; }
.guide-content section { max-width: 75ch; margin: 0 auto 42px; padding-bottom: 32px; border-bottom: 1px solid var(--scaid-line); scroll-margin-top: 100px; }
.guide-content section:last-child { border-bottom: 0; margin-bottom: 0; }
.guide-content section:focus { outline: none; }
.review-date { font-size: 14px; color: var(--scaid-muted); }
h1 { margin: 0 0 12px; font-size: 2rem; font-weight: 650; line-height: 1.25; color: var(--scaid-ink); }
h2 { margin: 0 0 18px; font-size: 25px; line-height: 1.3; color: #493c36; }
h3 { font-size: 20px; margin: 24px 0 12px; line-height: 1.4; }
p, li, dd { font-size: 16px; line-height: 1.75; }
p { margin: 0 0 16px; }
.lead { font-size: 19px; color: #59646c; }
ol, .guide-content ul { padding-left: 25px; margin: 0 0 20px; }
li { margin-bottom: 10px; }
.submenu li { margin: 0; }
a { color: var(--scaid-accent); text-decoration: underline; text-underline-offset: 3px; overflow-wrap: anywhere; }
.example-link { display: inline-block; margin: 2px 0 20px; font-weight: 600; }
.note { border: 1px solid var(--scaid-line); border-radius: 6px; padding: 12px 16px; background: var(--scaid-soft); }
.scope-summary { display: flex; gap: 24px; margin: 22px 0; padding: 22px; background: #f8f5f1; border-radius: 6px; }
.scope-summary div { flex: 1; }
.scope-summary strong { display: block; font-size: 30px; color: #88523c; line-height: 1.3; }
.scope-summary span { display: block; margin-top: 6px; font-size: 14px; }
.table-scroll { overflow-x: auto; margin: 18px 0; }
table { width: 100%; border-collapse: collapse; font-size: 15px; line-height: 1.6; }
th, td { text-align: left; padding: 12px; border: 1px solid #dde2e6; vertical-align: top; }
thead { background: #f1f4f6; }
code { font-size: .9em; color: #6c4434; background: #f4f4f3; padding: 2px 4px; border-radius: 3px; overflow-wrap: anywhere; }
pre { max-width: 100%; padding: 18px; background: #f4f4f3; border: 1px solid #e4e6e8; border-radius: 5px; overflow-x: auto; }
pre code { padding: 0; background: transparent; color: #33424b; overflow-wrap: normal; }
.troubleshooting-list dt { margin-top: 18px; color: #493c36; font-weight: 700; }
.troubleshooting-list dd { margin: 6px 0 0; }
.back-to-top { border: 1px solid #9b5f47; border-radius: 5px; background: #fff; color: #845038; padding: 10px 16px; cursor: pointer; }
@media (max-width: 991px) {
  .documentation-page { display: block; }
  .sidebar { position: static; max-height: none; width: 100%; border-right: 0; border-bottom: 1px solid var(--scaid-line); padding: 14px 16px; overflow: visible; }
  .sidebar-top { display: flex; justify-content: space-between; gap: 16px; align-items: center; }
  .sidebar-heading { margin: 0; font-size: 18px; }
  .contents-toggle { display: inline-flex; align-items: center; min-height: 44px; border: 1px solid #c7b8ae; border-radius: 6px; padding: 8px 12px; background: #fff; color: var(--scaid-accent); font-size: 14px; font-weight: 600; cursor: pointer; }
  .sidebar nav { display: none; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 8px; padding-top: 20px; }
  .sidebar nav.contents-open { display: grid; }
  .menu-group { margin-bottom: 0; }
  .menu-group-title { font-size: 12px; }
  .submenu button { min-height: 44px; font-size: 14px; padding: 9px 8px; }
  .guide-content { overflow: visible; padding: 26px 20px 40px; }
  .guide-content section { scroll-margin-top: 100px; margin-bottom: 28px; padding-bottom: 22px; }
  h1 { font-size: 1.75rem; }
  h2 { font-size: 23px; }
  .scope-summary { gap: 16px; padding: 16px; }
  .scope-summary strong { font-size: 25px; }
}
@media (max-width: 480px) {
  .sidebar nav { grid-template-columns: 1fr 1fr; }
  .menu-group:last-child { grid-column: 1 / -1; }
  .menu-group:last-child .submenu { display: grid; grid-template-columns: 1fr 1fr; }
  .guide-content { padding-right: 16px; padding-left: 16px; }
  .scope-summary { display: block; }
  .scope-summary div + div { margin-top: 16px; }
  th, td { padding: 8px; font-size: 13px; overflow-wrap: anywhere; }
  th[scope="row"], thead th:first-child { min-width: 76px; white-space: nowrap; overflow-wrap: normal; }
}
</style>

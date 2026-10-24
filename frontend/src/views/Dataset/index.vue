<template>
  <div class="container_box dataset-browser">
    <header v-if="!currentRow" class="dataset-intro">
      <div class="container">
        <h1>Dataset Browser</h1>
        <p class="intro-summary">Filter disease–tissue data objects by disease, accession or tissue.</p>
        <details class="analysis-description">
          <summary>About the available analyses</summary>
          <p>
            SCAID shows metadata, cell-type annotation views, cell-type marker genes
            (within-object contrasts, not donor-level disease-versus-control tests), GO and KEGG
            annotations of those markers, predicted cell–cell communication, and the downloadable
            RNA and pathway-score H5 files where published. Available modules vary by object.
          </p>
        </details>
      </div>
    </header>
    <div class="container dataset-workspace">
      <section v-if="!currentRow" class="dataset-filters" aria-labelledby="dataset-filter-title">
        <div class="section-heading">
          <h2 id="dataset-filter-title">Find a dataset</h2>
          <p>Choose a disease to narrow the dataset and tissue options, or search directly.</p>
        </div>
        <el-alert v-if="optionsError" :title="optionsError" type="error" show-icon :closable="false" />
        <el-form label-position="top" :model="form" @submit.prevent="getList">
          <el-row :gutter="20">
            <el-col :xs="24" :sm="8">
              <el-form-item label="Disease abbreviation" prop="abbreviation">
                <el-select v-model="form.abbreviation" clearable filterable placeholder="All diseases" style="width: 100%" @change="onAbbreviationChange">
                  <el-option v-for="item in options" :key="item.abbreviation" :label="abbreviationLabel(item)" :value="item.abbreviation" />
                </el-select>
              </el-form-item>
            </el-col>
            <el-col :xs="24" :sm="8">
              <el-form-item label="Dataset ID" prop="dataset_id">
                <el-select v-model="form.dataset_id" clearable filterable placeholder="All datasets" :disabled="!form.abbreviation" style="width: 100%" @change="onDatasetIdChange">
                  <el-option v-for="item in datasetIdOption" :key="item.dataset_id" :label="item.dataset_label || item.dataset_id" :value="item.dataset_id" />
                </el-select>
              </el-form-item>
            </el-col>
            <el-col :xs="24" :sm="8">
              <el-form-item label="Tissue" prop="tissue">
                <el-select v-model="form.tissue" clearable filterable placeholder="All tissues" :disabled="!form.dataset_id" style="width: 100%">
                  <el-option v-for="item in tissueOption" :key="item.tissue" :label="item.tissue_label || item.tissue" :value="item.tissue" />
                </el-select>
              </el-form-item>
            </el-col>
          </el-row>
          <div class="search-row">
            <el-form-item label="Dataset / tissue / file search" prop="name">
              <el-input v-model="form.name" clearable placeholder="Enter a disease abbreviation, dataset ID, tissue or file name" />
            </el-form-item>
            <div class="compare_btn">
              <el-button type="primary" native-type="submit" :loading="tableLoading">Search</el-button>
              <el-button @click="resetFilters">Reset</el-button>
            </div>
          </div>
        </el-form>
      </section>

      <section v-if="!currentRow" ref="resultsSection" class="dataset-results" aria-labelledby="dataset-results-title" tabindex="-1">
        <div class="results-heading">
          <h2 id="dataset-results-title" class="result-summary" aria-live="polite">{{ tableLoading ? 'Loading datasets…' : tableError ? 'Datasets unavailable' : `${tableData.length} datasets` }}</h2>
          <p v-if="tableData.length">Select a dataset to view figures and downloads.</p>
        </div>
        <el-alert v-if="tableError" :title="tableError" type="error" show-icon :closable="false" />
        <el-empty v-if="!tableLoading && !tableError && tableData.length === 0" description="No data objects match this query. Clear a filter or try a disease abbreviation, accession or tissue." />
        <template v-if="tableLoading || tableData.length">
          <p id="dataset-scroll-help" class="table-scroll-hint">Scroll across the table to see all columns. With keyboard focus, use the left and right arrow keys.</p>
          <div class="table-scroll-region" role="region" aria-label="Dataset results table" aria-describedby="dataset-scroll-help" tabindex="0">
            <el-table :data="tableData" v-loading="tableLoading" highlight-current-row :row-class-name="datasetRowClass" style="width: 100%; min-width: 940px" max-height="500" @row-click="row => openDataset(row)">
              <el-table-column label="Dataset ID" min-width="185">
                <template #default="scope">
                  <button type="button" class="dataset-open" :aria-label="`View ${scope.row.dataset_label || scope.row.dataset_id}, ${scope.row.tissue_label || scope.row.tissue}`" :aria-pressed="currentRow?.id === scope.row.id" @click.stop="openDataset(scope.row)">
                    <span>{{ scope.row.dataset_label || scope.row.dataset_id }}</span>
                    <span class="dataset-open-label">{{ currentRow?.id === scope.row.id ? 'Selected' : 'View details' }}</span>
                  </button>
                </template>
              </el-table-column>
              <el-table-column prop="abbreviation" label="Abbreviation" min-width="125" />
              <el-table-column label="Tissue" min-width="140">
                <template #default="scope">{{ scope.row.tissue_label || scope.row.tissue }}</template>
              </el-table-column>
              <el-table-column prop="sample_size" label="Sample labels" min-width="135" align="right" />
              <el-table-column label="Disease" min-width="270">
                <template #default="scope">{{ scope.row.disease_label || scope.row.disease_full_name || scope.row.disease_name || scope.row.abbreviation }}</template>
              </el-table-column>
              <el-table-column label="Cell records" min-width="125" align="right">
                <template #default="scope">{{ formatCount(scope.row.cell_count) }}</template>
              </el-table-column>
              <el-table-column label="Published views" min-width="175">
                <template #default="scope"><span class="module-summary">{{ moduleSummary(scope.row) }}</span></template>
              </el-table-column>
            </el-table>
          </div>
        </template>
        <p class="sample-label-note">Sample labels count identifiers in the source data. They may represent pooled samples or repeated observations and do not establish the number of independent donors.</p>
        <ShareQueryLink v-if="!tableLoading && !tableError" :href="shareHref" />
        <el-alert v-if="selectionError" :title="selectionError" type="warning" :closable="false" show-icon />
      </section>

      <section v-if="currentRow" class="dataset-details" aria-label="Selected dataset">
        <div class="selection-heading">
          <div>
            <h1 ref="selectionTitle" class="selection-title" tabindex="-1">{{ currentRow.abbreviation }} · {{ currentRow.dataset_label || currentRow.dataset_id }} · {{ currentRow.tissue_label || currentRow.tissue }}</h1>
            <p v-if="currentRow.disease_label || currentRow.disease_full_name || currentRow.disease_name">{{ currentRow.disease_label || currentRow.disease_full_name || currentRow.disease_name }}</p>
            <p v-if="currentRow.study_context" class="sample-label-note">{{ currentRow.study_context }}</p>
          </div>
          <div class="selection-utilities">
            <button type="button" class="back-to-results" @click="returnToResults">Back to datasets</button>
            <ShareQueryLink :href="shareHref" />
          </div>
        </div>
        <dl class="dataset-facts">
          <div><dt>Source</dt><dd><a v-if="sourceHref" :href="sourceHref" target="_blank" rel="noopener noreferrer">{{ currentRow.dataset_source || 'GEO' }} · {{ currentRow.dataset_label || currentRow.dataset_id }}<span class="sr-only"> (opens in a new tab)</span></a><span v-else>{{ currentRow.dataset_source || 'Not reported' }}</span></dd></div>
          <div><dt>Tissue</dt><dd>{{ currentRow.tissue_label || currentRow.tissue || 'Not reported' }}</dd></div>
          <div><dt>Cell records</dt><dd>{{ formatCount(currentRow.cell_count) }}</dd></div>
          <div><dt>Sample labels</dt><dd>{{ formatCount(currentRow.sample_size) }}</dd></div>
        </dl>
        <div class="detail-actions">
          <nav v-if="hasPublishedModule(currentRow, 'gene_umap') || hasPublishedModule(currentRow, 'kegg_score')" class="exploration-links" aria-label="Explore this dataset">
            <RouterLink v-if="hasPublishedModule(currentRow, 'gene_umap')" :to="explorationTarget('gene')">Gene expression <svg aria-hidden="true" viewBox="0 0 20 20"><path d="M4 10h12m-5-5 5 5-5 5" /></svg></RouterLink>
            <RouterLink v-if="hasPublishedModule(currentRow, 'kegg_score')" :to="explorationTarget('kegg')">Pathway scores <svg aria-hidden="true" viewBox="0 0 20 20"><path d="M4 10h12m-5-5 5 5-5 5" /></svg></RouterLink>
          </nav>
          <p v-else class="module-summary">Main figures only. Gene and pathway views are not published for this object.</p>
        </div>
        <details class="publication-details">
          <summary>Data notes &amp; limitations{{ publicationNotes(currentRow).length ? ` (${publicationNotes(currentRow).length})` : '' }}</summary>
          <p class="sample-label-note">Sample labels count identifiers in the source data. They may represent pooled samples or repeated observations and do not establish the number of independent donors. Cell records describe this processed object and are not a deduplicated atlas total.</p>
          <p v-for="note in publicationNotes(currentRow)" :key="note" class="publication-note">{{ note }}</p>
        </details>
        <div class="dataset-toolbar">
          <el-menu :default-active="activeName" mode="horizontal" :ellipsis="false" @select="handleMenuSelect" class="dataset-menu">
            <el-menu-item v-if="hasPublishedAnalysis(currentRow, 'Overview')" index="Overview">Overview</el-menu-item>
            <el-menu-item v-if="hasPublishedAnalysis(currentRow, 'Main_Marker')" index="Main_Marker">Marker</el-menu-item>
            <el-menu-item v-if="hasPublishedAnalysis(currentRow, 'MajorCellType_Go_Ht')" index="MajorCellType_Go_Ht">GO</el-menu-item>
            <el-menu-item v-if="hasPublishedAnalysis(currentRow, 'MajorCellType_KEGG_Ht')" index="MajorCellType_KEGG_Ht">KEGG</el-menu-item>
            <el-sub-menu v-if="hasPublishedModule(currentRow, 'cellchat')" index="CellChat">
              <template #title>CellChat</template>
              <el-menu-item v-for="item in availableCellchatList" :key="item.key" :index="item.key">{{ item.label }}</el-menu-item>
            </el-sub-menu>
          </el-menu>
          <el-dropdown v-if="hasPublishedModule(currentRow, 'downloads')" trigger="click" :disabled="downloadLoading || downList.length === 0" class="download-menu">
            <el-button :loading="downloadLoading" :disabled="downloadLoading || downList.length === 0">
              {{ downloadLoading ? 'Loading downloads' : downloadError ? 'Downloads unavailable' : downList.length ? `Download (${downList.length})` : 'No downloads available' }}
              <svg v-if="downList.length && !downloadLoading" aria-hidden="true" viewBox="0 0 20 20" class="download-chevron"><path d="m5 7.5 5 5 5-5" /></svg>
            </el-button>
            <template #dropdown>
              <el-dropdown-menu>
                <el-dropdown-item v-for="item in downList" :key="item.id">
                  <a :href="item.download_url" :download="getDownloadName(item)" class="download-link">{{ getDownloadLabel(item) }}</a>
                </el-dropdown-item>
              </el-dropdown-menu>
            </template>
          </el-dropdown>
        </div>
        <p v-if="analysisNotice" class="publication-note" role="status">{{ analysisNotice }}</p>
        <p v-if="selectedAnalysis" class="analysis-help"><strong>{{ selectedAnalysis.label }}.</strong> {{ selectedAnalysis.description }} These are precomputed predictions from expression data.</p>
        <el-alert v-if="downloadError" :title="downloadError" type="error" show-icon :closable="false">
          <el-button link type="primary" @click="loadDownloads(currentRow)">Retry downloads</el-button>
        </el-alert>
        <div class="menu-content" :aria-busy="loading">
          <div v-if="loading" class="loading-container" role="status">
            <p>Loading figures…</p>
            <el-skeleton :rows="4" animated />
          </div>
          <el-alert v-else-if="imageError" :title="imageError" type="error" show-icon :closable="false">
            <el-button link type="primary" @click="loadImages(imagePage)">Retry figures</el-button>
          </el-alert>
          <template v-else>
            <el-alert v-if="imageLoadError" title="Some image files could not be loaded. Please retry or report this dataset." type="error" show-icon :closable="false">
              <el-button link type="primary" @click="loadImages(imagePage)">Retry figures</el-button>
            </el-alert>
            <el-empty v-if="fileList.length === 0" description="No figures available for this selection" />
            <template v-else>
              <div class="figure-summary">
                <p class="result-summary" aria-live="polite">{{ imageCount }} figures · Page {{ imagePage }} of {{ Math.ceil(imageCount / IMAGE_PAGE_SIZE) }}</p>
                <p>Open a figure to inspect it at full resolution.</p>
              </div>
              <div v-for="group in imageGroups" :key="group.name" class="fileBox">
                <div v-for="item in group.files" :key="item.id" class="file_item">
                  <div class="figure-heading">
                    <h2 class="figure-title">{{ figureTitle(item) }}</h2>
                    <p class="figure-source">{{ item.name || item.pdf_name || activeName }}</p>
                  </div>
                  <a :href="getFullImageUrl(item.image)" target="_blank" rel="noopener" title="Open the full-resolution figure in a new tab">
                    <img :src="item.thumb_url || getFullImageUrl(item.image)" :style="getImageStyle(item)" :width="item.image_width || undefined" :height="item.image_height || undefined" :alt="`${currentRow.dataset_id}: ${item.name || activeName}`" loading="lazy" decoding="async" @error="imageLoadError = true" />
                    <span class="figure-link-label">Open full-resolution figure <svg aria-hidden="true" viewBox="0 0 20 20"><path d="M11 4h5v5M16 4l-8 8M8 4H4v12h12v-4" /></svg></span>
                  </a>
                </div>
              </div>
            </template>
            <el-pagination v-if="imageCount > IMAGE_PAGE_SIZE" :current-page="imagePage" :page-size="IMAGE_PAGE_SIZE" :total="imageCount" layout="prev, pager, next" @current-change="loadImages" />
          </template>
        </div>
      </section>
    </div>
  </div>
</template>

<script setup>
defineOptions({ name: 'DatasetBrowser' });
import { computed, nextTick, onBeforeUnmount, onMounted, reactive, ref, watch } from 'vue';
import { useRoute, useRouter } from 'vue-router';
import { getRequest } from '@/api/home';
import { getFullImageUrl } from '@/config';
import { datasetQuery, explorationQuery, querySignature, readDatasetQuery } from '@/utils/queryState';
import { hasPublishedAnalysis, hasPublishedModule, publicationNotes } from '@/utils/moduleAvailability';
import ShareQueryLink from '@/components/ShareQueryLink.vue';

const route = useRoute();
const router = useRouter();
const form = reactive({ abbreviation: '', dataset_id: '', tissue: '', name: '' });
const options = ref([]);
const optionsError = ref('');
const datasetIdOption = computed(() => options.value.find(item => item.abbreviation === form.abbreviation)?.children || []);
const tissueOption = computed(() => datasetIdOption.value.find(item => item.dataset_id === form.dataset_id)?.children || []);
const tableData = ref([]);
const tableLoading = ref(false);
const tableError = ref('');
const selectionError = ref('');
const appliedFilters = ref({ abbreviation: '', dataset_id: '', tissue: '', name: '' });
const currentRow = ref(null);
const resultsSection = ref(null);
const selectionTitle = ref(null);
const activeName = ref('Overview');
const analysisNotice = ref('');
const fileList = ref([]);
const loading = ref(false);
const imageError = ref('');
const imageLoadError = ref(false);
const imageCount = ref(0);
const imagePage = ref(1);
const IMAGE_PAGE_SIZE = 24;
const downList = ref([]);
const downloadLoading = ref(false);
const downloadError = ref('');
const cellchatList = [
  { key: 'CellChatOverview', label: 'Interaction overview', description: 'Overall interaction counts, strengths and cell-group signaling roles.' },
  { key: 'LRcircle', label: 'Ligand–receptor networks', description: 'Circle plots for individual ligand–receptor pairs.' },
  { key: 'LRcontribution', label: 'Ligand–receptor contributions', description: 'Contributions of ligand–receptor pairs to a signaling pathway.' },
  { key: 'SignalingCellRole', label: 'Cell signaling roles', description: 'Cell-group roles for the selected signaling pathway.' },
  { key: 'SignalingChord', label: 'Signaling chord diagrams', description: 'Chord diagrams of signaling between cell groups.' },
  { key: 'SignalingCircle', label: 'Signaling circle diagrams', description: 'Circle diagrams of signaling between cell groups.' },
  { key: 'SignalingHt', label: 'Signaling heatmaps', description: 'Heatmaps of signaling between cell groups.' },
  { key: 'SignalingRole', label: 'Signaling role summaries', description: 'Summary views of cell-group signaling roles.' },
  { key: 'Single', label: 'Cell-group network views', description: 'Network views by cell group, labeled Count or Weight in the source figures.' },
];
const availableCellchatList = computed(() => cellchatList.filter(item => hasPublishedAnalysis(currentRow.value, item.key)));
const analyses = ['Overview', 'Main_Marker', 'MajorCellType_Go_Ht', 'MajorCellType_KEGG_Ht', ...cellchatList.map(item => item.key)];
const selectedAnalysis = computed(() => cellchatList.find(item => item.key === activeName.value));
const sourceHref = computed(() => /^GSE\d+$/i.test(currentRow.value?.dataset_id || '')
  ? `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${encodeURIComponent(currentRow.value.dataset_id)}` : '');
const formatCount = value => Number.isFinite(Number(value)) && value !== null && value !== '' ? Number(value).toLocaleString('en-US') : 'Not reported';
function moduleSummary(row) {
  if (!row.module_availability) return 'Browse analyses';
  return row.module_availability.gene_umap && row.module_availability.kegg_score ? 'Main · Gene · Scores' : 'Main figures only';
}
function explorationTarget(kind) {
  return { path: kind === 'gene' ? '/Gene' : '/KEGG', query: explorationQuery(currentRow.value, kind) };
}
const shareHref = computed(() => router.resolve({ path: '/Dataset', query: datasetQuery(appliedFilters.value, currentRow.value?.id || '', activeName.value, imagePage.value) }).href);
// Downloads must stay on this origin's H5 endpoint; never render an API-supplied absolute or javascript: URL.
const SAFE_DOWNLOAD_URL = /^\/api\/h5-file\/\d+\/download\/$/;
// disease_label is the manifest's English label; disease_name holds the source-language name and is not shown.
const abbreviationLabel = (item) => {
  const name = item.disease_label || item.disease_full_name;
  return name && name !== item.abbreviation ? `${item.abbreviation} — ${name}` : item.abbreviation;
};
let tableRequestId = 0;
let imageRequestId = 0;
let downloadRequestId = 0;
let disposed = false;
let appliedFilterKey = '';
let routeVersion = 0;

const imageGroups = computed(() => {
  const groups = new Map();
  for (const item of fileList.value) {
    const name = activeName.value === 'Overview' ? (item.name || 'Overview') : activeName.value;
    if (!groups.has(name)) groups.set(name, []);
    groups.get(name).push(item);
  }
  return Array.from(groups, ([name, files]) => ({ name, files }));
});

function clearDetails() {
  imageRequestId += 1;
  downloadRequestId += 1;
  currentRow.value = null;
  activeName.value = 'Overview';
  analysisNotice.value = '';
  fileList.value = [];
  downList.value = [];
  imageCount.value = 0;
  imagePage.value = 1;
  loading.value = false;
  downloadLoading.value = false;
  imageError.value = '';
  imageLoadError.value = false;
  downloadError.value = '';
}

function syncQuery() {
  const query = datasetQuery(appliedFilters.value, currentRow.value?.id || '', activeName.value, imagePage.value);
  if (querySignature(query) !== querySignature(route.query)) router.push({ path: '/Dataset', query });
}

async function getList(syncUrl = true) {
  const requestId = ++tableRequestId;
  clearDetails();
  selectionError.value = '';
  appliedFilters.value = { ...form, name: form.name.trim() };
  appliedFilterKey = querySignature(appliedFilters.value);
  if (syncUrl !== false) syncQuery();
  tableData.value = [];
  tableLoading.value = true;
  tableError.value = '';
  try {
    const response = await getRequest('/cell-data-all/', { ...appliedFilters.value, format: 'json' });
    if (!Array.isArray(response)) throw new Error('Invalid dataset response');
    if (requestId !== tableRequestId || disposed) return;
    tableData.value = response;
  } catch {
    if (requestId === tableRequestId && !disposed) tableError.value = 'Unable to load datasets. Please try Search again.';
  } finally {
    if (requestId === tableRequestId && !disposed) tableLoading.value = false;
  }
}

function handleClickRow(row, syncUrl = true, analysis = 'Overview', page = 1) {
  clearDetails();
  currentRow.value = row;
  const available = hasPublishedAnalysis(row, analysis);
  activeName.value = available ? analysis : 'Overview';
  imagePage.value = available ? page : 1;
  if (!available) analysisNotice.value = 'The linked analysis has not been published for this object. Showing Overview.';
  // Reset the tab before starting either request for the new dataset.
  if (syncUrl !== false) syncQuery();
  else if (!available) router.replace({ path: '/Dataset', query: datasetQuery(appliedFilters.value, row.id, activeName.value, imagePage.value) });
  loadImages(imagePage.value, false);
  if (hasPublishedModule(row, 'downloads')) loadDownloads(row);
}

function datasetRowClass({ row }) {
  return currentRow.value?.id === row.id ? 'selected-dataset' : '';
}

async function openDataset(row) {
  handleClickRow(row);
  await nextTick();
  selectionTitle.value?.focus({ preventScroll: true });
  selectionTitle.value?.scrollIntoView({ block: 'start' });
}

async function returnToResults() {
  clearDetails();
  syncQuery();
  await nextTick();
  resultsSection.value?.focus({ preventScroll: true });
  resultsSection.value?.scrollIntoView({ block: 'start' });
}

async function loadImages(page = 1, syncUrl = true) {
  if (!currentRow.value) return;
  const requestId = ++imageRequestId;
  const row = currentRow.value;
  const tab = activeName.value;
  imagePage.value = page;
  if (syncUrl !== false) syncQuery();
  fileList.value = [];
  imageCount.value = 0;
  imageError.value = '';
  imageLoadError.value = false;
  loading.value = true;
  try {
    const response = await getRequest('/pdf-images/', {
      format: 'json', dataset_record_id: row.id, abbreviation: row.abbreviation,
      dataset_id: row.dataset_id, tissue: row.tissue, pdf_name: tab,
      page, page_size: IMAGE_PAGE_SIZE,
    });
    if (response.results?.code !== 200 || !Array.isArray(response.results.data)) throw new Error('Invalid figure response');
    if (requestId !== imageRequestId || disposed) return;
    fileList.value = response.results.data;
    imageCount.value = Number(response.count) || fileList.value.length;
  } catch {
    if (requestId === imageRequestId && !disposed) imageError.value = 'Unable to load figures for this dataset.';
  } finally {
    if (requestId === imageRequestId && !disposed) loading.value = false;
  }
}

function handleMenuSelect(index) {
  if (!analyses.includes(index) || !hasPublishedAnalysis(currentRow.value, index)) return;
  analysisNotice.value = '';
  activeName.value = index;
  loadImages(1);
}

async function loadDownloads(row) {
  const requestId = ++downloadRequestId;
  downList.value = [];
  downloadError.value = '';
  downloadLoading.value = true;
  try {
    const response = await getRequest('/h5-file/', { format: 'json', dataset_record_id: row.id });
    if (!Array.isArray(response) || response.some(item => !SAFE_DOWNLOAD_URL.test(item.download_url || ''))) throw new Error('Invalid download response');
    if (requestId === downloadRequestId && !disposed) downList.value = response;
  } catch {
    if (requestId === downloadRequestId && !disposed) downloadError.value = 'Unable to load downloads for this dataset.';
  } finally {
    if (requestId === downloadRequestId && !disposed) downloadLoading.value = false;
  }
}

function onAbbreviationChange() {
  form.dataset_id = '';
  form.tissue = '';
}
function onDatasetIdChange() { form.tissue = ''; }
function resetFilters() {
  Object.assign(form, { abbreviation: '', dataset_id: '', tissue: '', name: '' });
  getList();
}
function getDownloadName(item) {
  const relativeName = String(item.relative_path || '').split(/[\\/]/).pop();
  return String(item.name || item.file_name || (/\.(h5|h5ad|hdf5)$/i.test(relativeName) ? relativeName : `dataset-${item.id}.h5`)).split(/[\\/]/).pop();
}
function getDownloadLabel(item) {
  const method = item.method || String(item.relative_path || '').split(/[\\/]/).find(segment => /^(RNA|AUCell|UCell|singscore)$/i.test(segment));
  return method ? `${method} — ${getDownloadName(item)}` : getDownloadName(item);
}
function getImageStyle(item) { return { width: item.image_width ? `${item.image_width}px` : '100%', maxWidth: '100%', height: 'auto' }; }
function figureTitle(item) {
  const source = String(item.name || item.pdf_name || activeName.value).replace(/\.(pdf|png|jpe?g)$/i, '');
  const titles = {
    Umap_MajorCellType: 'Major cell types',
    Umap_Sample: 'Sample labels',
    Umap_library: 'Technical library labels',
    Umap_Leiden: 'Cell clusters',
    'Umap_RNA_snn_res.1': 'Cell clusters',
    MajorCellType_CellPro: 'Sample composition',
    Main_Marker: 'Cell-type markers',
    MajorCellType_Go_Ht: 'GO enrichment',
    MajorCellType_GO_Ht: 'GO enrichment',
    MajorCellType_KEGG_Ht: 'KEGG enrichment',
  };
  return titles[source] || selectedAnalysis.value?.label || titles[activeName.value] || source;
}

watch(() => querySignature(route.query), async () => {
  const version = ++routeVersion;
  const restored = readDatasetQuery(route.query, analyses);
  const key = querySignature(restored.filters);
  Object.assign(form, restored.filters);
  if (key !== appliedFilterKey) await getList(false);
  if (version !== routeVersion || disposed || tableError.value) return;
  if (!restored.record) {
    if (currentRow.value) clearDetails();
    selectionError.value = '';
    return;
  }
  const row = tableData.value.find(item => String(item.id) === restored.record);
  if (!row) {
    clearDetails();
    selectionError.value = 'The linked dataset is not available for these filters. Clear the filters or choose another dataset.';
    return;
  }
  selectionError.value = '';
  if (String(currentRow.value?.id) !== restored.record) {
    handleClickRow(row, false, restored.analysis, restored.page);
  } else if (activeName.value !== restored.analysis || imagePage.value !== restored.page) {
    if (hasPublishedAnalysis(row, restored.analysis)) {
      activeName.value = restored.analysis;
      analysisNotice.value = '';
      loadImages(restored.page, false);
    } else {
      handleClickRow(row, false, restored.analysis, restored.page);
    }
  }
}, { immediate: true });

onMounted(async () => {
  try {
    const response = await getRequest('/cell-data/', { format: 'json' });
    if (response.code !== 200 || !Array.isArray(response.data)) throw new Error('Invalid filter response');
    if (!disposed) options.value = response.data;
  } catch {
    if (!disposed) optionsError.value = 'Unable to load filter options. You can still search by dataset, tissue or file name.';
  }
});
onBeforeUnmount(() => { disposed = true; });
</script>

<style lang="scss" scoped>
.dataset-browser {
  color: var(--scaid-ink, #2b2320);
  min-width: 0;
  padding-bottom: 40px;
  --dataset-muted: var(--scaid-muted, #6b5f5a);
  --dataset-accent: var(--scaid-accent, #8a503a);
  --dataset-soft: var(--scaid-soft, #f8f5f2);
  --dataset-line: var(--scaid-line, #e5ddd7);
}
.dataset-intro {
  padding: 30px 0 26px;
  background: var(--dataset-soft);
  border-bottom: 1px solid var(--dataset-line);
  h1 { margin: 0 0 10px; font-size: clamp(2rem, 2.7vw, 2.4rem); font-weight: 600; letter-spacing: -.025em; line-height: 1.2; }
}
.intro-summary { max-width: 75ch; margin: 0; color: var(--dataset-muted); font-size: 1rem; line-height: 1.6; }
.analysis-description {
  margin-top: 14px;
  max-width: 75ch;
  color: var(--dataset-muted);
  font-size: 0.9rem;
  summary { width: fit-content; min-height: 44px; cursor: pointer; padding: 10px 0; color: var(--dataset-accent); text-underline-offset: 3px; }
  summary:hover { text-decoration: underline; }
  p { margin: 8px 0 0; line-height: 1.7; }
}
.dataset-workspace { padding-top: 24px; }
.dataset-filters { padding: 24px; border: 1px solid var(--dataset-line); border-radius: 12px; background: var(--dataset-soft); }
.section-heading {
  margin-bottom: 20px;
  h2 { margin: 0 0 6px; font-size: 1.25rem; font-weight: 600; letter-spacing: -.015em; }
  p { margin: 0; color: var(--dataset-muted); font-size: 0.9rem; line-height: 1.5; }
}
.search-row { display: flex; align-items: flex-end; gap: 16px; }
.search-row :deep(.el-form-item) { flex: 1; min-width: 0; margin-bottom: 0; }
.compare_btn { display: flex; flex-shrink: 0; padding-bottom: 0; }
.dataset-filters :deep(.el-form-item__label) { margin-bottom: 8px; color: var(--scaid-ink); font-size: .875rem; font-weight: 600; }
.dataset-filters :deep(.el-select__wrapper),
.dataset-filters :deep(.el-input__wrapper),
.compare_btn :deep(.el-button) { min-height: 44px; }
.dataset-filters :deep(.el-select__wrapper), .dataset-filters :deep(.el-input__wrapper) { background: #fff; }
.dataset-filters :deep(.el-select__wrapper.is-disabled) { background: var(--dataset-soft); }
.dataset-results { margin-top: 28px; scroll-margin-top: 100px; }
.results-heading { display: flex; align-items: baseline; justify-content: space-between; flex-wrap: wrap; column-gap: 16px; row-gap: 6px; margin-bottom: 14px; }
.result-summary { margin: 0; color: var(--dataset-muted); font-size: 0.9rem; font-variant-numeric: tabular-nums; }
.results-heading h2 { color: var(--scaid-ink); font-size: 1.25rem; font-weight: 600; letter-spacing: -.015em; }
.results-heading p { margin: 0; font-size: 0.85rem; color: var(--dataset-muted); }
.table-scroll-hint { display: none; margin: 0 0 10px; font-size: 0.8rem; line-height: 1.5; color: var(--dataset-muted); }
.sample-label-note { max-width: 75ch; margin: 12px 0; color: var(--dataset-muted); font-size: .85rem; line-height: 1.6; }
.publication-note { max-width: 75ch; margin: 12px 0; color: var(--dataset-muted); font-size: .9rem; line-height: 1.6; }
.analysis-help { max-width: 75ch; margin: 12px 0; color: var(--dataset-muted); font-size: .9rem; line-height: 1.6; }
.figure-title { margin: 0; color: var(--scaid-ink); font-size: 1rem; font-weight: 600; line-height: 1.5; overflow-wrap: anywhere; }
.table-scroll-region { overflow-x: auto; border: 1px solid var(--dataset-line); border-radius: 12px; scrollbar-width: thin; scrollbar-color: var(--dataset-muted) var(--dataset-soft); }
.table-scroll-region :deep(.el-table) {
  --el-table-border-color: var(--dataset-line);
  --el-table-header-bg-color: var(--dataset-soft);
  --el-table-header-text-color: var(--scaid-ink, #2b2320);
  --el-table-text-color: var(--scaid-ink, #2b2320);
  --el-table-row-hover-bg-color: var(--dataset-soft);
  --el-table-current-row-bg-color: #f2e9e2;
  font-size: 0.875rem;
  font-variant-numeric: tabular-nums;
}
.table-scroll-region :deep(.el-table th.el-table__cell) { padding: 12px 0; font-weight: 600; }
.table-scroll-region :deep(.el-table td.el-table__cell) { padding: 10px 0; }
.table-scroll-region :deep(.el-table__row) { cursor: pointer; }
.table-scroll-region :deep(.selected-dataset td.el-table__cell) { background: #f2e9e2; }
.dataset-open { display: flex; flex-direction: column; gap: 3px; width: 100%; min-height: 44px; padding: 2px 0; border: 0; border-radius: 3px; background: none; color: var(--dataset-accent); font: inherit; text-align: left; cursor: pointer; }
.dataset-open > span:first-child { font-weight: 600; text-underline-offset: 3px; }
.dataset-open:hover > span:first-child { text-decoration: underline; text-decoration-thickness: 1px; }
.dataset-open-label { color: var(--dataset-muted); font-size: 0.75rem; line-height: 1.3; }
.dataset-open[aria-pressed='true'] .dataset-open-label { color: var(--dataset-accent); font-weight: 600; }
.dataset-details { min-width: 0; }
.selection-heading { display: flex; align-items: flex-start; justify-content: space-between; gap: 24px; margin-bottom: 24px; }
.selection-heading > div { min-width: 0; }
.selection-utilities { display: flex; align-items: center; flex-wrap: wrap; gap: 4px 18px; flex-shrink: 0; }
.selection-title { margin: 0; font-size: clamp(1.5rem, 2.2vw, 1.9rem); font-weight: 600; letter-spacing: -.02em; line-height: 1.3; overflow-wrap: anywhere; scroll-margin-top: 100px; }
.selection-heading p { margin: 7px 0 0; color: var(--dataset-muted); line-height: 1.5; }
.back-to-results { min-height: 44px; flex-shrink: 0; padding: 8px 0; border: 0; background: none; color: var(--dataset-accent); font: inherit; font-size: 0.85rem; text-decoration: underline; text-underline-offset: 3px; cursor: pointer; }
.dataset-toolbar { display: flex; align-items: center; gap: 16px; border-bottom: 1px solid var(--dataset-line); }
.dataset-menu { flex: 1; min-width: 0; background: transparent; border-bottom: 0; --el-menu-text-color: var(--dataset-muted); --el-menu-active-color: var(--dataset-accent); --el-menu-hover-bg-color: var(--dataset-soft); --el-menu-horizontal-height: 52px; }
.dataset-menu :deep(.el-menu-item), .dataset-menu :deep(.el-sub-menu__title) { padding: 0 18px; font-size: .9rem; font-weight: 500; }
.dataset-menu :deep(.el-menu-item.is-active), .dataset-menu :deep(.el-sub-menu.is-active > .el-sub-menu__title) { font-weight: 600; }
.download-menu { flex-shrink: 0; padding-bottom: 8px; }
.download-menu :deep(.el-button) { min-height: 44px; background: var(--dataset-soft); border-color: var(--dataset-line); color: var(--dataset-accent); }
.download-menu :deep(.el-button:not(:disabled):hover) { background: #fff; border-color: var(--dataset-accent); }
.download-chevron { width: 16px; height: 16px; margin-left: 6px; fill: none; stroke: currentColor; stroke-width: 1.5; }
.download-link { display: block; width: 100%; padding: 10px 0; color: inherit; overflow-wrap: anywhere; max-width: calc(100vw - 64px); }
.menu-content { padding: 20px 0 0; min-height: 200px; margin-bottom: 24px; }
.figure-summary { display: flex; justify-content: space-between; flex-wrap: wrap; gap: 6px 16px; }
.figure-summary > p { margin: 0; color: var(--dataset-muted); font-size: 0.85rem; line-height: 1.5; }
.fileBox { margin: 24px 0; }
.figure-heading { padding: 18px 20px 0; }
.figure-source { margin: 4px 0 0; color: var(--dataset-muted); font-size: .8125rem; line-height: 1.5; overflow-wrap: anywhere; }
.dataset-facts { display: grid; grid-template-columns: minmax(0, 1.6fr) repeat(3, minmax(0, 1fr)); gap: 24px; margin: 0 0 16px; padding: 20px 24px; border-radius: 12px; background: var(--dataset-soft); }
.dataset-facts > div { min-width: 0; }
.dataset-facts dt { margin: 0 0 6px; color: var(--dataset-muted); font-size: .8125rem; font-weight: 400; line-height: 1.4; }
.dataset-facts dd { margin: 0; color: var(--scaid-ink); font-size: 1rem; font-weight: 600; line-height: 1.5; font-variant-numeric: tabular-nums; overflow-wrap: anywhere; }
.dataset-facts a { color: var(--dataset-accent); text-underline-offset: 3px; }
.detail-actions { display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 4px 20px; margin-bottom: 8px; }
.exploration-links { display: flex; flex-wrap: wrap; gap: 0 24px; }
.exploration-links a { display: inline-flex; align-items: center; gap: 8px; min-height: 44px; color: var(--dataset-accent); font-size: .9rem; font-weight: 500; text-underline-offset: 4px; }
.exploration-links svg { width: 16px; height: 16px; flex-shrink: 0; fill: none; stroke: currentColor; stroke-width: 1.5; stroke-linecap: round; stroke-linejoin: round; }
.exploration-links a:hover { text-decoration: underline; text-decoration-thickness: 1px; }
.module-summary { color: var(--dataset-muted); font-size: .85rem; line-height: 1.5; }
.publication-details { margin-bottom: 12px; }
.publication-details summary { width: fit-content; min-height: 44px; padding: 10px 0; color: var(--dataset-accent); font-size: .9rem; cursor: pointer; }
.publication-details summary:hover { text-decoration: underline; text-underline-offset: 3px; }
.publication-details[open] { margin-bottom: 20px; }
.file_item { margin-bottom: 24px; border: 1px solid var(--dataset-line); border-radius: 12px; background: #fff; transition: border-color 160ms ease-out; }
.file_item:hover, .file_item:focus-within { border-color: var(--dataset-accent); }
.file_item > a { display: block; padding: 20px; border-radius: 12px; color: var(--dataset-accent); text-decoration: none; }
.file_item img { display: block; margin: 0 auto; object-fit: contain; }
.figure-link-label { display: flex; justify-content: flex-end; align-items: center; gap: 6px; margin-top: 16px; font-size: .8125rem; line-height: 1.5; }
.figure-link-label svg { flex-shrink: 0; width: 16px; height: 16px; fill: none; stroke: currentColor; stroke-width: 1.5; stroke-linecap: round; stroke-linejoin: round; }
.loading-container { padding: 24px 0; min-height: 230px; color: var(--dataset-muted); }
.loading-container p { margin: 0 0 24px; }
.menu-content :deep(.el-pagination) { justify-content: center; flex-wrap: wrap; gap: 6px; margin-top: 24px; }
.el-alert { margin: 12px 0; }
.dataset-browser :is(button, summary, a, [tabindex]):focus-visible { outline: 2px solid var(--dataset-accent); outline-offset: 3px; }
@media (max-width: 1100px) {
  .table-scroll-hint { display: block; }
  .dataset-menu :deep(.el-menu-item), .dataset-menu :deep(.el-sub-menu__title) { padding-right: 14px; padding-left: 14px; }
}
@media (max-width: 767px) {
  .dataset-intro { padding: 24px 0; }
  .dataset-intro h1 { font-size: 1.75rem; }
  .intro-summary { font-size: 1rem; }
  .dataset-workspace { padding-top: 20px; }
  .dataset-filters { padding: 20px 16px; }
  .search-row { flex-direction: column; align-items: stretch; gap: 16px; }
  .compare_btn :deep(.el-button) { flex: 1; }
  .dataset-results { margin-top: 22px; }
  .results-heading { display: block; }
  .results-heading p { margin-top: 7px; line-height: 1.5; }
  .dataset-facts { grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 18px 20px; padding: 18px 16px; }
  .selection-heading { flex-direction: column-reverse; gap: 12px; margin-bottom: 16px; }
  .selection-title { font-size: 1.5rem; }
  .selection-utilities { width: 100%; justify-content: space-between; }
  .exploration-links { gap: 0 16px; }
  .exploration-links a { gap: 6px; font-size: .875rem; }
  .detail-actions { margin-bottom: 4px; }
  .back-to-results { min-height: 44px; }
  .dataset-toolbar { flex-direction: column; align-items: stretch; gap: 12px; }
  .dataset-menu { display: flex; flex-wrap: wrap; width: 100%; height: auto; }
  .dataset-menu :deep(.el-menu-item), .dataset-menu :deep(.el-sub-menu) { flex: 1 0 auto; height: 46px; }
  .dataset-menu :deep(.el-sub-menu__title) { height: 46px; justify-content: center; }
  .download-menu { padding-bottom: 14px; }
  .download-menu :deep(.el-button) { width: 100%; min-height: 44px; }
  .menu-content { padding-top: 16px; }
  .figure-heading { padding: 16px 16px 0; }
  .file_item > a { padding: 16px; }
  .figure-link-label { justify-content: flex-start; }
}
@media (prefers-reduced-motion: reduce) {
  .loading-container :deep(.el-skeleton__item) { animation: none; }
  .file_item { transition: none; }
}
</style>

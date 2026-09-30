<template>
  <div class="container_box dataset-browser">
    <header class="dataset-intro">
      <div class="container">
        <h1>Dataset Browser</h1>
        <p class="intro-summary">Filter disease–tissue data objects by disease, accession or tissue.</p>
        <details class="analysis-description">
          <summary>About the available analyses</summary>
          <p>
            For each object SCAID shows its metadata, cell-type annotation views, cell-type marker genes
            (within-object contrasts, not donor-level disease-versus-control tests), GO and KEGG
            annotations of those markers, predicted cell–cell communication, and the downloadable
            RNA and pathway-score H5 files.
          </p>
        </details>
      </div>
    </header>
    <div class="container dataset-workspace">
      <section class="dataset-filters" aria-labelledby="dataset-filter-title">
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

      <section ref="resultsSection" class="dataset-results" aria-labelledby="dataset-results-title" tabindex="-1">
        <div class="results-heading">
          <h2 id="dataset-results-title" class="result-summary" aria-live="polite">{{ tableLoading ? 'Loading datasets…' : tableError ? 'Datasets unavailable' : `${tableData.length} datasets` }}</h2>
          <p v-if="tableData.length">Select a dataset to view figures and downloads.</p>
        </div>
        <el-alert v-if="tableError" :title="tableError" type="error" show-icon :closable="false" />
        <el-empty v-if="!tableLoading && !tableError && tableData.length === 0" description="No data objects match this query. Clear a filter or try a disease abbreviation, accession or tissue." />
        <template v-if="tableLoading || tableData.length">
          <p id="dataset-scroll-help" class="table-scroll-hint">Scroll across the table to see all columns. With keyboard focus, use the left and right arrow keys.</p>
          <div class="table-scroll-region" role="region" aria-label="Dataset results table" aria-describedby="dataset-scroll-help" tabindex="0">
            <el-table :data="tableData" v-loading="tableLoading" highlight-current-row :row-class-name="datasetRowClass" style="width: 100%; min-width: 940px" max-height="500" @row-click="handleClickRow">
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
              <el-table-column prop="id" label="ID" width="70" align="right" />
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
            <h2 ref="selectionTitle" class="selection-title" tabindex="-1">{{ currentRow.abbreviation }} · {{ currentRow.dataset_label || currentRow.dataset_id }} · {{ currentRow.tissue_label || currentRow.tissue }}</h2>
            <p v-if="currentRow.disease_label || currentRow.disease_full_name || currentRow.disease_name">{{ currentRow.disease_label || currentRow.disease_full_name || currentRow.disease_name }}</p>
          </div>
          <button type="button" class="back-to-results" @click="returnToResults">Back to datasets</button>
        </div>
        <div class="dataset-toolbar">
          <el-menu :default-active="activeName" mode="horizontal" :ellipsis="false" @select="handleMenuSelect" class="dataset-menu">
            <el-menu-item index="Overview">Overview</el-menu-item>
            <el-menu-item index="Main_Marker">Marker</el-menu-item>
            <el-menu-item index="MajorCellType_Go_Ht">GO</el-menu-item>
            <el-menu-item index="MajorCellType_KEGG_Ht">KEGG</el-menu-item>
            <el-sub-menu index="CellChat">
              <template #title>CellChat</template>
              <el-menu-item v-for="item in cellchatList" :key="item.key" :index="item.key">{{ item.label }}</el-menu-item>
            </el-sub-menu>
          </el-menu>
          <el-dropdown trigger="click" :disabled="downloadLoading || downList.length === 0" class="download-menu">
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
                <h3 v-if="activeName === 'Overview'">{{ group.name.replace(/\.pdf$/i, '') }}</h3>
                <div v-for="item in group.files" :key="item.id" class="file_item">
                  <h4 class="figure-title">{{ figureTitle(item) }}</h4>
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
import { datasetQuery, querySignature, readDatasetQuery } from '@/utils/queryState';
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
  { key: 'LRcircle', label: 'Ligand–receptor networks', description: 'Circle plots for individual ligand–receptor pairs.' },
  { key: 'LRcontribution', label: 'Ligand–receptor contributions', description: 'Contributions of ligand–receptor pairs to a signaling pathway.' },
  { key: 'SignalingCellRole', label: 'Cell signaling roles', description: 'Cell-group roles for the selected signaling pathway.' },
  { key: 'SignalingChord', label: 'Signaling chord diagrams', description: 'Chord diagrams of signaling between cell groups.' },
  { key: 'SignalingCircle', label: 'Signaling circle diagrams', description: 'Circle diagrams of signaling between cell groups.' },
  { key: 'SignalingHt', label: 'Signaling heatmaps', description: 'Heatmaps of signaling between cell groups.' },
  { key: 'SignalingRole', label: 'Signaling role summaries', description: 'Summary views of cell-group signaling roles.' },
  { key: 'Single', label: 'Cell-group network views', description: 'Network views by cell group, labeled Count or Weight in the source figures.' },
];
const analyses = ['Overview', 'Main_Marker', 'MajorCellType_Go_Ht', 'MajorCellType_KEGG_Ht', ...cellchatList.map(item => item.key)];
const selectedAnalysis = computed(() => cellchatList.find(item => item.key === activeName.value));
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
  activeName.value = analysis;
  imagePage.value = page;
  // Reset the tab before starting either request for the new dataset.
  if (syncUrl !== false) syncQuery();
  loadImages(page, false);
  loadDownloads(row);
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

function returnToResults() {
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
  if (!analyses.includes(index)) return;
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
function figureTitle(item) { return String(item.name || item.pdf_name || selectedAnalysis.value?.label || activeName.value).replace(/\.pdf$/i, ''); }

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
    activeName.value = restored.analysis;
    loadImages(restored.page, false);
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
  padding: 36px 0 28px;
  background: var(--dataset-soft);
  border-bottom: 1px solid var(--dataset-line);
  h1 { margin: 0 0 12px; font-size: 2.25rem; font-weight: 600; letter-spacing: -0.025em; line-height: 1.2; }
}
.intro-summary { margin: 0; font-size: 1.05rem; line-height: 1.6; }
.analysis-description {
  margin-top: 14px;
  max-width: 76ch;
  color: var(--dataset-muted);
  font-size: 0.9rem;
  summary { width: fit-content; cursor: pointer; padding: 5px 0; color: var(--dataset-accent); text-underline-offset: 3px; }
  summary:hover { text-decoration: underline; }
  p { margin: 8px 0 0; line-height: 1.7; }
}
.dataset-workspace { padding-top: 30px; }
.dataset-filters { padding-bottom: 26px; border-bottom: 1px solid var(--dataset-line); }
.section-heading {
  margin-bottom: 20px;
  h2 { margin: 0 0 6px; font-size: 1.2rem; font-weight: 600; }
  p { margin: 0; color: var(--dataset-muted); font-size: 0.9rem; line-height: 1.5; }
}
.search-row { display: flex; align-items: flex-end; gap: 16px; }
.search-row :deep(.el-form-item) { flex: 1; min-width: 0; margin-bottom: 0; }
.compare_btn { display: flex; flex-shrink: 0; padding-bottom: 0; }
.dataset-filters :deep(.el-form-item__label) { color: var(--scaid-ink, #2b2320); font-weight: 500; }
.dataset-filters :deep(.el-select__wrapper),
.dataset-filters :deep(.el-input__wrapper),
.compare_btn :deep(.el-button) { min-height: 42px; }
.dataset-filters :deep(.el-select__wrapper.is-disabled) { background: var(--dataset-soft); }
.dataset-results { margin-top: 24px; scroll-margin-top: 100px; }
.results-heading { display: flex; align-items: baseline; justify-content: space-between; flex-wrap: wrap; column-gap: 16px; row-gap: 6px; margin-bottom: 14px; }
.result-summary { margin: 0; color: var(--dataset-muted); font-size: 0.9rem; font-variant-numeric: tabular-nums; }
.results-heading h2 { color: var(--scaid-ink, #2b2320); font-size: 1.1rem; font-weight: 600; }
.results-heading p { margin: 0; font-size: 0.85rem; color: var(--dataset-muted); }
.table-scroll-hint { display: none; margin: 0 0 10px; font-size: 0.8rem; line-height: 1.5; color: var(--dataset-muted); }
.sample-label-note { max-width: 75ch; margin: 12px 0; color: var(--dataset-muted); font-size: .85rem; line-height: 1.6; }
.analysis-help { max-width: 75ch; margin: 12px 0; color: var(--dataset-muted); font-size: .9rem; line-height: 1.6; }
.figure-title { margin: 0 0 10px; font-size: .95rem; font-weight: 600; line-height: 1.5; overflow-wrap: anywhere; }
.table-scroll-region { overflow-x: auto; border: 1px solid var(--dataset-line); border-radius: 8px; }
.table-scroll-region :deep(.el-table) {
  --el-table-border-color: var(--dataset-line);
  --el-table-header-bg-color: var(--dataset-soft);
  --el-table-header-text-color: var(--scaid-ink, #2b2320);
  --el-table-text-color: var(--scaid-ink, #2b2320);
  --el-table-row-hover-bg-color: #faf6f2;
  --el-table-current-row-bg-color: #f2e9e2;
  font-size: 0.875rem;
  font-variant-numeric: tabular-nums;
}
.table-scroll-region :deep(.el-table th.el-table__cell) { padding: 12px 0; font-weight: 600; }
.table-scroll-region :deep(.el-table td.el-table__cell) { padding: 8px 0; }
.table-scroll-region :deep(.el-table__row) { cursor: pointer; }
.table-scroll-region :deep(.selected-dataset td.el-table__cell) { background: #f2e9e2; }
.dataset-open { display: flex; flex-direction: column; gap: 3px; width: 100%; min-height: 44px; padding: 2px 0; border: 0; border-radius: 3px; background: none; color: var(--dataset-accent); font: inherit; text-align: left; cursor: pointer; }
.dataset-open > span:first-child { font-weight: 600; text-decoration: underline; text-underline-offset: 3px; }
.dataset-open:hover > span:first-child { text-decoration-thickness: 2px; }
.dataset-open-label { color: var(--dataset-muted); font-size: 0.75rem; line-height: 1.3; }
.dataset-open[aria-pressed='true'] .dataset-open-label { color: var(--dataset-accent); font-weight: 600; }
.dataset-details { margin-top: 38px; padding-top: 28px; border-top: 1px solid var(--dataset-line); }
.selection-heading { display: flex; align-items: flex-start; justify-content: space-between; gap: 16px; margin-bottom: 20px; }
.selection-heading > div { min-width: 0; }
.selection-title { margin: 0; font-size: 1.3rem; font-weight: 600; line-height: 1.45; overflow-wrap: anywhere; scroll-margin-top: 100px; }
.selection-heading p { margin: 7px 0 0; color: var(--dataset-muted); line-height: 1.5; }
.back-to-results { flex-shrink: 0; padding: 8px 0; border: 0; background: none; color: var(--dataset-accent); font: inherit; font-size: 0.85rem; text-decoration: underline; text-underline-offset: 3px; cursor: pointer; }
.dataset-toolbar { display: flex; align-items: center; gap: 16px; border-bottom: 1px solid var(--dataset-line); }
.dataset-menu { flex: 1; min-width: 0; border-bottom: 0; --el-menu-text-color: var(--dataset-muted); --el-menu-active-color: var(--dataset-accent); --el-menu-hover-bg-color: var(--dataset-soft); --el-menu-horizontal-height: 52px; }
.dataset-menu :deep(.el-menu-item), .dataset-menu :deep(.el-sub-menu__title) { padding: 0 18px; font-weight: 500; }
.download-menu { flex-shrink: 0; padding-bottom: 8px; }
.download-menu :deep(.el-button) { min-height: 40px; }
.download-chevron { width: 16px; height: 16px; margin-left: 6px; fill: none; stroke: currentColor; stroke-width: 1.5; }
.download-link { display: block; width: 100%; padding: 10px 0; color: inherit; overflow-wrap: anywhere; max-width: calc(100vw - 64px); }
.menu-content { padding: 20px 0 0; min-height: 200px; margin-bottom: 24px; }
.figure-summary { display: flex; justify-content: space-between; flex-wrap: wrap; gap: 6px 16px; }
.figure-summary > p { margin: 0; color: var(--dataset-muted); font-size: 0.85rem; line-height: 1.5; }
.fileBox { margin: 26px 0; }
.fileBox h3 { margin: 0 0 12px; font-size: 1rem; font-weight: 600; overflow-wrap: anywhere; }
.file_item { margin-bottom: 16px; border: 1px solid var(--dataset-line); background: #fff; }
.file_item > a { display: block; padding: 16px; color: var(--dataset-accent); }
.file_item img { display: block; margin: 0 auto; }
.figure-link-label { display: flex; justify-content: flex-end; align-items: center; gap: 6px; margin-top: 14px; font-size: 0.8rem; line-height: 1.5; }
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
  .dataset-intro { padding: 28px 0 22px; }
  .dataset-intro h1 { font-size: 1.8rem; }
  .intro-summary { font-size: 1rem; }
  .dataset-workspace { padding-top: 24px; }
  .dataset-filters { padding-bottom: 24px; }
  .search-row { flex-direction: column; align-items: stretch; gap: 16px; }
  .compare_btn :deep(.el-button) { flex: 1; }
  .dataset-results { margin-top: 22px; }
  .results-heading { display: block; }
  .results-heading p { margin-top: 7px; line-height: 1.5; }
  .dataset-details { margin-top: 28px; padding-top: 24px; }
  .selection-heading { flex-direction: column; gap: 4px; }
  .selection-title { font-size: 1.15rem; }
  .back-to-results { min-height: 40px; }
  .dataset-toolbar { flex-direction: column; align-items: stretch; gap: 12px; }
  .dataset-menu { display: flex; flex-wrap: wrap; width: 100%; height: auto; }
  .dataset-menu :deep(.el-menu-item), .dataset-menu :deep(.el-sub-menu) { flex: 1 0 auto; height: 46px; }
  .dataset-menu :deep(.el-sub-menu__title) { height: 46px; justify-content: center; }
  .download-menu { padding-bottom: 14px; }
  .download-menu :deep(.el-button) { width: 100%; min-height: 44px; }
  .menu-content { padding-top: 16px; }
  .file_item > a { padding: 12px; }
  .figure-link-label { justify-content: flex-start; }
}
@media (prefers-reduced-motion: reduce) {
  .loading-container :deep(.el-skeleton__item) { animation: none; }
}
</style>

<template>
  <section class="marker-browser" aria-labelledby="markers-title" :aria-busy="loading">
    <div class="marker-heading">
      <div>
        <h2 id="markers-title">Cell-type markers</h2>
        <p>Each cell type is compared with the remaining cells in this dataset. Positive and negative log2 fold changes are retained.</p>
      </div>
      <a v-if="available" :href="fullDownloadUrl" download class="marker-download">Download full table (TSV)</a>
    </div>

    <el-form v-if="available" :model="form" class="marker-filters" label-position="top" @submit.prevent="applyFilters">
      <el-form-item label="Gene contains">
        <el-input v-model="form.gene" clearable maxlength="128" placeholder="For example, NKG7" />
      </el-form-item>
      <el-form-item label="Cell type">
        <el-select v-model="form.cellType" clearable filterable placeholder="All cell types">
          <el-option v-for="group in cellTypes" :key="group" :label="displayCellType(group)" :value="group" />
        </el-select>
      </el-form-item>
      <div class="marker-filter-actions">
        <el-button type="primary" native-type="submit" :loading="loading">Search markers</el-button>
        <el-button :disabled="loading" @click="resetFilters">Reset</el-button>
      </div>
    </el-form>

    <el-alert v-if="error" :title="error" type="error" show-icon :closable="false">
      <el-button link type="primary" @click="loadMarkers">Retry markers</el-button>
    </el-alert>
    <div v-if="loading && !available" class="marker-loading" role="status">
      <p>Loading marker results…</p>
      <el-skeleton :rows="4" animated />
    </div>
    <el-empty v-else-if="!error && !available" description="A searchable marker table has not been published for this dataset. The Marker figures tab contains the available plots." />
    <template v-if="available && !error">
      <div class="marker-results-heading">
        <p role="status" aria-live="polite">{{ loading ? 'Loading marker results…' : `${count.toLocaleString('en-US')} marker rows` }}</p>
        <a v-if="count && !loading" :href="filteredDownloadUrl" download class="marker-download">Download these results (TSV)</a>
      </div>
      <p id="marker-scroll-help" class="marker-table-help">Select a column heading to sort. Scroll across the table on smaller screens.</p>
      <div class="marker-scroll" tabindex="0" role="region" aria-label="Cell-type marker results" aria-describedby="marker-scroll-help">
        <el-table :data="rows" v-loading="loading" :default-sort="{ prop: 'p_val_adj', order: 'ascending' }" style="width: 100%; min-width: 900px" @sort-change="sortMarkers">
          <el-table-column prop="gene" label="Gene" min-width="140" sortable="custom" />
          <el-table-column prop="group1" label="Cell type" min-width="170" sortable="custom">
            <template #default="scope">{{ displayCellType(scope.row.group1) }}</template>
          </el-table-column>
          <el-table-column prop="avg_log2FC" label="log2 fold change" min-width="160" align="right" sortable="custom">
            <template #default="scope">{{ formatEffect(scope.row.avg_log2FC) }}</template>
          </el-table-column>
          <el-table-column prop="p_val_adj" label="Adjusted P" min-width="150" align="right" sortable="custom">
            <template #default="scope">{{ formatP(scope.row.p_val_adj) }}</template>
          </el-table-column>
          <el-table-column prop="pct.1" label="Expressing: target" min-width="170" align="right" sortable="custom">
            <template #default="scope">{{ formatPercent(scope.row['pct.1']) }}</template>
          </el-table-column>
          <el-table-column prop="pct.2" label="Expressing: rest" min-width="165" align="right" sortable="custom">
            <template #default="scope">{{ formatPercent(scope.row['pct.2']) }}</template>
          </el-table-column>
          <template #empty>No marker rows match these filters. Clear the gene search or choose another cell type.</template>
        </el-table>
      </div>
      <p class="marker-footnote">Expression percentages describe cells in the target group and remaining cells. P values of 0 are displayed as reported in the source output. Downloads retain all original columns and precision.</p>
      <el-pagination v-if="count > PAGE_SIZE" :current-page="page" :page-size="PAGE_SIZE" :total="count" :disabled="loading" layout="prev, pager, next" @current-change="changePage" />
    </template>
  </section>
</template>

<script setup>
import { computed, onBeforeUnmount, onMounted, reactive, ref } from 'vue';
import { getRequest } from '@/api/home';

const props = defineProps({ datasetRecordId: { type: [Number, String], required: true } });
const PAGE_SIZE = 25;
const form = reactive({ gene: '', cellType: '' });
const applied = reactive({ gene: '', cellType: '' });
const rows = ref([]);
const cellTypes = ref([]);
const available = ref(false);
const count = ref(0);
const page = ref(1);
const ordering = ref('p_val_adj');
const loading = ref(false);
const error = ref('');
let controller;
let requestVersion = 0;
let disposed = false;
const fullDownloadUrl = computed(() => `/api/markers/download/?${new URLSearchParams({ dataset_record_id: String(props.datasetRecordId) })}`);
const filteredDownloadUrl = computed(() => {
  const params = new URLSearchParams({ dataset_record_id: String(props.datasetRecordId), ordering: ordering.value });
  if (applied.gene) params.set('gene', applied.gene);
  if (applied.cellType) params.set('cell_type', applied.cellType);
  return `/api/markers/download/?${params}`;
});
const displayCellType = value => value === 'Unknonw' ? 'Unknown' : value;
const formatP = value => value === 0 ? '0' : Number(value).toExponential(2);
const formatEffect = value => Number(value).toFixed(3);
const formatPercent = value => `${(Number(value) * 100).toFixed(1)}%`;

async function loadMarkers() {
  controller?.abort();
  controller = new AbortController();
  const version = ++requestVersion;
  loading.value = true;
  error.value = '';
  rows.value = [];
  try {
    const response = await getRequest('/markers/', {
      dataset_record_id: props.datasetRecordId, gene: applied.gene,
      cell_type: applied.cellType, ordering: ordering.value, page: page.value, page_size: PAGE_SIZE,
    }, { signal: controller.signal });
    if (version !== requestVersion || disposed) return;
    if (response.dataset_record_id !== Number(props.datasetRecordId) || typeof response.available !== 'boolean'
      || !Array.isArray(response.results) || !Array.isArray(response.cell_types)
      || !Number.isSafeInteger(response.count) || response.count < 0) throw new Error('Invalid marker response');
    rows.value = response.results;
    available.value = response.available;
    cellTypes.value = response.cell_types;
    count.value = response.count;
  } catch {
    if (version === requestVersion && !disposed) error.value = 'Unable to load marker results for this dataset. Please retry.';
  } finally {
    if (version === requestVersion && !disposed) loading.value = false;
  }
}
function applyFilters() {
  Object.assign(applied, { gene: form.gene.trim(), cellType: form.cellType });
  page.value = 1;
  loadMarkers();
}
function resetFilters() {
  Object.assign(form, { gene: '', cellType: '' });
  applyFilters();
}
function sortMarkers({ prop, order }) {
  ordering.value = prop && order ? `${order === 'descending' ? '-' : ''}${prop}` : 'p_val_adj';
  page.value = 1;
  loadMarkers();
}
function changePage(value) { page.value = value; loadMarkers(); }
onMounted(loadMarkers);
onBeforeUnmount(() => { disposed = true; requestVersion += 1; controller?.abort(); });
</script>

<style lang="scss" scoped>
.marker-browser { min-width: 0; color: var(--scaid-ink, #2b2320); }
.marker-heading { display: flex; justify-content: space-between; align-items: baseline; flex-wrap: wrap; gap: 12px 24px; margin-bottom: 22px; }
.marker-heading h2 { margin: 0 0 8px; font-size: 1.2rem; font-weight: 600; }
.marker-heading p { max-width: 72ch; margin: 0; color: var(--scaid-muted, #6b5f5a); line-height: 1.6; font-size: .9rem; }
.marker-download { display: inline-flex; align-items: center; min-height: 44px; color: var(--scaid-accent, #8a503a); text-decoration: underline; text-underline-offset: 4px; font-size: .9rem; }
.marker-download:hover { text-decoration-thickness: 2px; }
.marker-filters { display: grid; grid-template-columns: minmax(160px, 1fr) minmax(200px, 1fr) auto; align-items: end; gap: 16px; padding-bottom: 20px; border-bottom: 1px solid var(--scaid-line, #e5ddd7); }
.marker-filters :deep(.el-form-item) { min-width: 0; margin-bottom: 0; }
.marker-filters :deep(.el-form-item__label) { color: var(--scaid-ink, #2b2320); }
.marker-filters :deep(.el-select) { width: 100%; }
.marker-filters :deep(.el-input__wrapper), .marker-filters :deep(.el-select__wrapper), .marker-filter-actions :deep(.el-button) { min-height: 42px; }
.marker-filter-actions { display: flex; }
.marker-results-heading { display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 4px 16px; margin: 14px 0 4px; }
.marker-results-heading p, .marker-table-help { margin: 0; color: var(--scaid-muted, #6b5f5a); font-size: .85rem; line-height: 1.5; }
.marker-table-help { margin-bottom: 12px; }
.marker-scroll { overflow-x: auto; border: 1px solid var(--scaid-line, #e5ddd7); border-radius: 8px; }
.marker-scroll :deep(.el-table) { --el-table-border-color: var(--scaid-line, #e5ddd7); --el-table-header-bg-color: var(--scaid-soft, #f8f5f2); --el-table-header-text-color: var(--scaid-ink, #2b2320); --el-table-text-color: var(--scaid-ink, #2b2320); font-variant-numeric: tabular-nums; font-size: .875rem; }
.marker-scroll :deep(.el-table th.el-table__cell) { padding: 12px 0; }
.marker-footnote { max-width: 75ch; color: var(--scaid-muted, #6b5f5a); font-size: .8rem; line-height: 1.6; margin: 12px 0; }
.marker-loading { padding: 20px 0; min-height: 220px; }
.marker-loading p { margin: 0 0 20px; }
.marker-browser :deep(.el-pagination) { justify-content: center; flex-wrap: wrap; gap: 6px; margin-top: 20px; }
.marker-browser :is(a, button, [tabindex]):focus-visible { outline: 2px solid var(--scaid-accent, #8a503a); outline-offset: 3px; }
@media (max-width: 900px) { .marker-filters { grid-template-columns: 1fr 1fr; } .marker-filter-actions { grid-column: 1 / -1; } }
@media (max-width: 600px) {
  .marker-filters { grid-template-columns: minmax(0, 1fr); }
  .marker-filter-actions :deep(.el-button) { flex: 1; }
  .marker-filters :deep(.el-input__inner), .marker-filters :deep(.el-select__selected-item) { font-size: 16px; }
  .marker-scroll :deep(.el-table) { font-size: .875rem; }
}
</style>

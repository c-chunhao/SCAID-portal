<template>
  <div class="container_box exploration-page">
    <div class="jumbotron jumbotron-fluid exploration-intro">
      <div class="container">
        <h1 class="display-5">{{ title }} Exploration</h1>
        <p><slot /></p>
      </div>
    </div>

    <div class="container exploration-content">
      <el-form class="exploration-form" :model="form" label-position="top" @submit.prevent="submit">
        <div class="primary-query-row">
          <el-form-item class="query-field" :label="kind === 'gene' ? 'Gene' : 'KEGG pathway'" prop="pdf_name" :error="state.validationError" required>
            <el-select
              ref="querySelect"
              v-model="form.pdf_name"
              filterable
              remote
              clearable
              :remote-method="fetchOptionsDebounced"
              :loading="state.optionsLoading"
              :placeholder="kind === 'gene' ? 'Type a gene symbol to find suggestions' : 'Type a pathway ID to find suggestions'"
              :no-data-text="state.optionsError ? 'Suggestions unavailable' : 'No matching suggestions'"
              @clear="search.fetchOptions('')"
            >
              <el-option v-for="item in state.options" :key="item.id || item.name" :label="item.name" :value="item.name" />
              <template v-if="state.optionsTotal" #footer>
                <div class="suggestion-footer">
                  <span>{{ state.options.length }} of {{ state.optionsTotal }} suggestions</span>
                  <el-button
                    v-if="state.options.length < state.optionsTotal"
                    size="small"
                    :loading="state.optionsLoading"
                    @mousedown.prevent
                    @click="search.loadMoreOptions"
                  >Load more</el-button>
                </div>
              </template>
            </el-select>
          </el-form-item>
          <div class="form-actions">
            <el-button type="primary" native-type="submit" :loading="state.loading">Search plots</el-button>
            <p>Leave optional filters empty to see all available results.</p>
          </div>
        </div>
        <div v-if="state.optionsError" class="inline-recovery">
          <el-alert :title="state.optionsError" type="error" :closable="false" show-icon />
          <el-button :loading="state.optionsLoading" @click="search.retryOptions">Retry suggestions</el-button>
        </div>

        <details ref="optionalFilters" class="optional-filters" :open="!state.searched">
          <summary>{{ filterSummary }}</summary>
          <fieldset class="filter-fields">
          <legend class="sr-only">Optional filters</legend>
          <div class="filter-grid" :class="{ 'filter-grid--pathway': kind === 'kegg' }">
            <el-form-item label="Disease (abbreviation)" prop="abbreviation">
              <el-select v-model="form.abbreviation" filterable clearable :loading="filtersLoading" placeholder="All diseases" @change="onAbbreviationChange">
                <el-option v-for="item in diseaseOptions" :key="item.abbreviation" :label="item.disease_label && item.disease_label !== item.abbreviation ? `${item.abbreviation} — ${item.disease_label}` : item.abbreviation" :value="item.abbreviation" />
              </el-select>
            </el-form-item>
            <el-form-item label="Dataset ID" prop="dataset_id">
              <el-select v-model="form.dataset_id" filterable clearable :loading="filtersLoading" placeholder="All datasets" @change="onDatasetChange">
                <el-option v-for="item in datasetOptions" :key="item.dataset_id" :label="item.dataset_label || item.dataset_id" :value="item.dataset_id" />
              </el-select>
            </el-form-item>
            <el-form-item label="Tissue" prop="tissue">
              <el-select v-model="form.tissue" filterable clearable :loading="filtersLoading" placeholder="All tissues">
                <el-option v-for="item in tissueOptions" :key="item.tissue" :label="item.tissue_label || item.tissue" :value="item.tissue" />
              </el-select>
            </el-form-item>
            <el-form-item v-if="kind === 'kegg'" label="Scoring method" prop="kegg_param">
              <el-select v-model="form.kegg_param" clearable placeholder="All methods">
                <el-option v-for="method in methods" :key="method" :label="method" :value="method" />
              </el-select>
            </el-form-item>
          </div>
          <div v-if="filterError" class="inline-recovery">
            <el-alert :title="filterError" type="error" :closable="false" show-icon />
            <el-button :loading="filtersLoading" @click="fetchFilters">Retry filters</el-button>
          </div>
          </fieldset>
        </details>
        <p class="publication-scope">Filters include objects with published {{ kind === 'gene' ? 'gene-expression' : 'pathway-score' }} maps. Objects with only overview and annotation figures remain available in the <router-link to="/Dataset">Dataset Browser</router-link>.</p>
      </el-form>

      <section class="search-results" aria-labelledby="results-heading" :aria-busy="state.loading">
        <div class="results-heading">
          <h2 id="results-heading">Plots</h2>
          <p v-if="state.searched && !state.loading && !state.error && state.total" role="status">
            <strong>{{ state.total }} found</strong><span>Showing {{ firstResult }}–{{ lastResult }}</span>
          </p>
        </div>
        <div v-if="state.loading" class="loading-results">
          <p role="status">Loading {{ title }} plots…</p>
          <div class="file-box" aria-hidden="true">
            <div v-for="index in 3" :key="index" class="plot-skeleton">
              <div class="skeleton-image" />
              <div class="skeleton-caption" />
            </div>
          </div>
        </div>
        <div v-else-if="state.error" class="results-recovery">
          <el-alert :title="state.error" type="error" :closable="false" show-icon />
          <el-button @click="search.changePage(state.page)">Retry search</el-button>
        </div>
        <div v-else-if="state.searched && !state.total" class="results-placeholder" role="status">
          <h3>No matching plots</h3>
          <p>Try another {{ kind === 'gene' ? 'gene' : 'pathway' }}, or clear an optional filter and search again.</p>
        </div>
        <template v-else-if="state.searched">
          <div class="results-toolbar">
            <div class="page-size-control">
              <label for="plots-per-page">Plots per page</label>
              <el-select id="plots-per-page" :model-value="state.pageSize" @update:model-value="changePageSize">
                <el-option v-for="size in [10, 20, 50, 100]" :key="size" :label="String(size)" :value="size" />
              </el-select>
            </div>
            <div class="pagination-control">
              <span class="page-position">Page {{ state.page }} of {{ totalPages }}</span>
              <el-pagination
                :current-page="state.page"
                :page-size="state.pageSize"
                :pager-count="5"
                :total="state.total"
                :layout="compactPagination ? 'prev, next' : 'prev, pager, next'"
                @update:current-page="changePage"
              />
            </div>
            <ShareQueryLink :href="shareHref" :disabled="state.loading" />
          </div>
          <div class="file-box">
            <figure v-for="item in state.files" :key="item.id" class="file-item">
              <FigurePreview :src="item.thumb_url || getFullImageUrl(item.image)" :full-src="getFullImageUrl(item.image)" :alt="`${title} plot for ${form.pdf_name} in ${figureLabel(item)}`" square />
              <figcaption>
                <span class="figure-scope">{{ figureLabel(item) }}</span>
                <span class="figure-hint">Select figure to open full resolution</span>
              </figcaption>
            </figure>
          </div>
          <div v-if="state.total > state.pageSize" class="results-footer">
            <span class="page-position">Page {{ state.page }} of {{ totalPages }}</span>
            <el-pagination
              :current-page="state.page"
              :page-size="state.pageSize"
              :pager-count="5"
              :total="state.total"
              :layout="compactPagination ? 'prev, next' : 'prev, pager, next'"
              @update:current-page="changePage"
            />
          </div>
        </template>
        <div v-else class="results-placeholder">
          <h3>{{ form.pdf_name ? 'Ready to search' : `Choose a ${kind === 'gene' ? 'gene' : 'pathway'} to begin` }}</h3>
          <p>{{ form.pdf_name ? 'Select Search plots to view results for your current selection.' : `Type a ${kind === 'gene' ? 'gene symbol' : 'pathway ID'} above and select a suggestion, then search to view the available plots.` }}</p>
        </div>
      </section>
    </div>
  </div>
</template>

<script setup>
import { computed, onBeforeUnmount, onMounted, reactive, ref, watch } from 'vue'
import { useDebounceFn, useMediaQuery } from '@vueuse/core'
import { useRoute, useRouter } from 'vue-router'
import { getRequest } from '@/api/home'
import { getFullImageUrl } from '@/config'
import { createExplorationSearch, createExplorationState, requestErrorMessage } from '@/utils/explorationSearch'
import { explorationQuery, querySignature, readExplorationQuery } from '@/utils/queryState'
import FigurePreview from '@/components/FigurePreview.vue'
import ShareQueryLink from '@/components/ShareQueryLink.vue'

const props = defineProps({ kind: { type: String, required: true } })
const route = useRoute()
const router = useRouter()
const title = computed(() => props.kind === 'gene' ? 'Gene' : 'KEGG')
const querySelect = ref()
const optionalFilters = ref()
// Native disclosure state can change without a Vue render. Close it explicitly
// when submitting or restoring a feature, including an unchanged query.
const closeOptionalFilters = () => optionalFilters.value?.removeAttribute('open')
const compactPagination = useMediaQuery('(max-width: 640px)')
const form = reactive({ pdf_name: '', abbreviation: '', dataset_id: '', tissue: '', kegg_param: '' })
const state = reactive(createExplorationState())
const filterSummary = computed(() => {
  const selected = [form.abbreviation, form.dataset_id, form.tissue, props.kind === 'kegg' ? form.kegg_param : ''].filter(Boolean)
  return selected.length ? `Optional filters · ${selected.join(' · ')}` : 'Optional filters · All published datasets'
})
const search = createExplorationSearch({ request: getRequest, kind: props.kind, state })
// One suggestion request per pause in typing instead of one per keystroke.
const fetchOptionsDebounced = useDebounceFn((query) => search.fetchOptions(query), 300)
const filters = ref([])
const filtersLoading = ref(false)
const filterError = ref('')
const methods = ['AUCell', 'singscore', 'UCell']
let filterController
let filterVersion = 0
const diseaseOptions = computed(() => filters.value.filter(item => item.abbreviation))
const filteredDatasets = computed(() => {
  const groups = form.abbreviation
    ? filters.value.filter(item => item.abbreviation === form.abbreviation)
    : filters.value
  return groups.flatMap(item => item.children || []).filter(item => item.dataset_id)
})
const datasetOptions = computed(() => [...new Map(filteredDatasets.value.map(item => [item.dataset_id, item])).values()])
const tissueOptions = computed(() => {
  const datasets = form.dataset_id
    ? filteredDatasets.value.filter(item => item.dataset_id === form.dataset_id)
    : filteredDatasets.value
  const tissues = datasets.flatMap(item => item.children || []).filter(item => item.tissue)
  return [...new Map(tissues.map(item => [item.tissue, item])).values()]
})
const firstResult = computed(() => (state.page - 1) * state.pageSize + 1)
const lastResult = computed(() => (state.page - 1) * state.pageSize + state.files.length)
const totalPages = computed(() => Math.max(1, Math.ceil(state.total / state.pageSize)))
const shareHref = computed(() => router.resolve({ path: route.path, query: explorationQuery(form, props.kind, state.page, state.pageSize) }).href)

const fetchFilters = async () => {
  const version = ++filterVersion
  filterController?.abort()
  filterController = new AbortController()
  filtersLoading.value = true
  filterError.value = ''
  try {
    const response = await getRequest('/cell-data/', { format: 'json', module: props.kind === 'gene' ? 'gene_umap' : 'kegg_score' }, { signal: filterController.signal })
    if (version !== filterVersion) return
    if (!Array.isArray(response.data)) throw new Error('Invalid filter response')
    filters.value = response.data
  } catch (error) {
    if (version === filterVersion) filterError.value = requestErrorMessage(error, 'Disease, dataset and tissue filters')
  } finally {
    if (version === filterVersion) filtersLoading.value = false
  }
}

// 'IBD_Colitis_inflamed/GSE121380/colon/GeneUmap' -> 'IBD_Colitis_inflamed · GSE121380 · colon'
const figureLabel = (item) => {
  const parts = String(item.relative_path || '').split('/').filter(Boolean)
  const curated = [item.condition_abbreviation, item.dataset_label || item.dataset_id, item.tissue_label].filter(Boolean)
  const scope = curated.length ? curated.join(' · ') : parts.slice(0, 3).join(' · ')
  const method = parts[3] === 'KEGGScore' && parts[4] ? ` · ${parts[4]}` : ''
  return scope ? scope + method : (item.name || '')
}
const onAbbreviationChange = () => { form.dataset_id = ''; form.tissue = '' }
const onDatasetChange = () => { form.tissue = '' }
let restoringRoute = false
let appliedQuery = ''
const applyRouteQuery = async () => {
  const restored = readExplorationQuery(route.query, props.kind)
  const signature = querySignature(explorationQuery(restored, props.kind, restored.page, restored.pageSize))
  if (signature === appliedQuery) return
  appliedQuery = signature
  restoringRoute = true
  Object.assign(form, {
    pdf_name: restored.pdf_name, abbreviation: restored.abbreviation,
    dataset_id: restored.dataset_id, tissue: restored.tissue, kegg_param: restored.kegg_param
  })
  restoringRoute = false
  state.pageSize = restored.pageSize
  search.resetResults()
  // Changing the URL while this component is reused cancels both result and
  // suggestion requests; version guards prevent late responses replacing it.
  search.fetchOptions(restored.pdf_name)
  if (restored.pdf_name) {
    closeOptionalFilters()
    await search.submit(form, restored.page)
  }
}
const updateQuery = async (page = 1, pageSize = state.pageSize) => {
  const query = explorationQuery(form, props.kind, page, pageSize)
  const signature = querySignature(query)
  if (signature === appliedQuery) {
    state.pageSize = pageSize
    await search.submit(form, page)
  } else {
    await router.push({ path: route.path, query })
  }
}
const submit = async () => {
  if (!String(form.pdf_name || '').trim()) {
    await search.submit(form)
    querySelect.value?.focus()
    return
  }
  closeOptionalFilters()
  await updateQuery()
}
const changePage = (page) => updateQuery(page)
const changePageSize = (size) => updateQuery(1, size)
watch(form, () => { if (!restoringRoute) search.resetResults() }, { flush: 'sync' })
watch(() => querySignature(route.query), applyRouteQuery, { immediate: true })
onMounted(fetchFilters)
onBeforeUnmount(() => {
  filterVersion += 1
  filterController?.abort()
  search.dispose()
})
</script>

<style lang="scss" scoped>
.exploration-page {
  padding-bottom: 3.5rem;
  color: var(--scaid-ink);
}
.exploration-intro {
  margin-bottom: 0;
  padding: 30px 0 26px;
  background: var(--scaid-soft);
  border-bottom: 1px solid var(--scaid-line);
  .display-5 { margin: 0 0 10px; text-align: left; font-size: clamp(2rem, 2.7vw, 2.4rem); font-weight: 600; letter-spacing: -.025em; line-height: 1.2; }
  p { max-width: 75ch; margin: 0; font-size: 1rem; line-height: 1.6; color: var(--scaid-muted); }
  :deep(a) { color: var(--scaid-accent); text-decoration: underline; text-underline-offset: .18em; }
}
.exploration-content { padding-top: 24px; }
.exploration-form { max-width: 100%; padding: 24px; background: var(--scaid-soft); border: 1px solid var(--scaid-line); border-radius: 12px; }
.primary-query-row { display: grid; grid-template-columns: minmax(0, 640px) minmax(180px, 1fr); align-items: start; gap: 24px; }
.query-field { min-width: 0; margin-bottom: 0; }
.el-select { width: 100%; }
.exploration-form :deep(.el-form-item__label) { margin-bottom: 8px; color: var(--scaid-ink); font-size: .875rem; font-weight: 600; line-height: 1.4; }
.exploration-form :deep(.el-select__wrapper) { min-height: 44px; background: #fff; }
.query-field :deep(.el-select__wrapper) { min-height: 48px; }
.optional-filters { min-width: 0; margin: 20px 0 0; border-top: 1px solid var(--scaid-line); }
.optional-filters summary { min-height: 44px; padding: 12px 0 10px; color: var(--scaid-accent); font-size: .875rem; line-height: 1.5; cursor: pointer; overflow-wrap: anywhere; }
.optional-filters summary:hover { text-decoration: underline; text-underline-offset: 3px; }
.optional-filters summary:focus-visible { outline: 2px solid var(--scaid-accent); outline-offset: 3px; }
.filter-fields { min-width: 0; margin: 0; padding: 8px 0 0; border: 0; }
.filter-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 20px; }
.filter-grid--pathway { grid-template-columns: repeat(4, minmax(0, 1fr)); }
.filter-grid .el-form-item { min-width: 0; margin-bottom: 0; }
.form-actions { display: flex; align-items: flex-start; flex-direction: column; gap: 8px; padding-top: 28px; }
.form-actions .el-button { min-height: 48px; padding: 0 24px; font-weight: 600; }
.form-actions p { max-width: 34ch; margin: 0; color: var(--scaid-muted); font-size: .8125rem; line-height: 1.5; }
.publication-scope { max-width: 75ch; margin: 12px 0 0; color: var(--scaid-muted); font-size: .8125rem; line-height: 1.6; }
.publication-scope a { color: var(--scaid-accent); text-decoration: underline; text-underline-offset: .18em; }
.inline-recovery { display: flex; align-items: center; flex-wrap: wrap; gap: .75rem; margin: 1rem 0; }
.inline-recovery .el-alert { flex: 1 1 280px; }
.inline-recovery .el-button, .results-recovery .el-button { min-height: 44px; }
.suggestion-footer { display: flex; align-items: center; flex-wrap: wrap; justify-content: space-between; gap: .5rem; white-space: normal; font-size: .85rem; }
.search-results { margin-top: 28px; }
.results-heading { display: flex; align-items: baseline; flex-wrap: wrap; justify-content: space-between; gap: .5rem 1.5rem; margin-bottom: 16px; }
.results-heading h2 { margin: 0; font-size: 1.25rem; font-weight: 600; letter-spacing: -.015em; }
.results-heading p { display: flex; flex-wrap: wrap; gap: .75rem; margin: 0; font-size: .875rem; color: var(--scaid-muted); }
.results-heading strong { color: var(--scaid-ink); font-weight: 600; }
.results-placeholder { padding: 24px; background: var(--scaid-soft); border-radius: 12px; }
.results-placeholder h3 { margin: 0 0 .5rem; font-size: 1rem; font-weight: 600; }
.results-placeholder p, .loading-results > p { max-width: 72ch; margin: 0; color: var(--scaid-muted); font-size: .95rem; line-height: 1.6; }
.results-recovery { display: grid; justify-items: start; gap: 1rem; }
.results-toolbar { display: flex; flex-wrap: wrap; justify-content: space-between; align-items: center; gap: 1rem; margin-bottom: 16px; }
.page-size-control { display: flex; align-items: center; gap: .75rem; }
.page-size-control label, .page-position { margin: 0; color: var(--scaid-muted); font-size: .8125rem; font-variant-numeric: tabular-nums; }
.page-size-control .el-select { width: 80px; }
.page-size-control :deep(.el-select__wrapper) { min-height: 44px; }
.pagination-control { display: flex; align-items: center; flex-wrap: wrap; gap: .75rem; }
.el-pagination { --el-pagination-button-width: 40px; --el-pagination-button-height: 44px; max-width: 100%; padding: 0; }
.results-footer { display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 1rem; margin-top: 24px; padding-top: 16px; border-top: 1px solid var(--scaid-line); }
.file-box { display: grid; grid-template-columns: repeat(auto-fill, minmax(min(100%, 300px), 1fr)); gap: 24px; }
.file-item {
  min-width: 0;
  margin: 0;
  background: #fff;
  border: 1px solid var(--scaid-line);
  border-radius: 12px;
  transition: border-color 160ms ease-out;
  &:hover, &:focus-within { border-color: var(--scaid-accent); }
  figcaption { padding: 16px 20px 18px; border-top: 1px solid var(--scaid-line); font-size: .9rem; line-height: 1.5; overflow-wrap: anywhere; }
}
.figure-scope { display: block; color: var(--scaid-ink); font-weight: 600; }
.figure-hint { display: block; margin-top: 6px; color: var(--scaid-muted); font-size: .8125rem; }
.loading-results .file-box { margin-top: 20px; }
.plot-skeleton { padding: 12px; border: 1px solid var(--scaid-line); border-radius: 12px; }
.skeleton-image { aspect-ratio: 1; background: var(--scaid-soft); border-radius: 6px; }
.skeleton-caption { width: 65%; height: .85rem; margin: 1rem 0 .5rem; background: var(--scaid-line); border-radius: 3px; }
@media (max-width: 991px) {
  .primary-query-row { grid-template-columns: minmax(0, 1fr) 200px; gap: 20px; }
  .filter-grid, .filter-grid--pathway { grid-template-columns: repeat(2, minmax(0, 1fr)); }
}
@media (max-width: 640px) {
  .exploration-page { padding-bottom: 2.5rem; }
  .exploration-intro { padding: 24px 0; }
  .exploration-intro .display-5 { font-size: 1.75rem; }
  .exploration-intro p { font-size: 1rem; }
  .exploration-content { padding-top: 20px; }
  .exploration-form { padding: 20px 16px; }
  .primary-query-row { display: block; }
  .filter-grid, .filter-grid--pathway { grid-template-columns: minmax(0, 1fr); gap: 16px; }
  .form-actions { align-items: stretch; padding-top: 12px; gap: 8px; }
  .form-actions .el-button { width: 100%; }
  .form-actions p { max-width: none; }
  .search-results { margin-top: 24px; }
  .results-placeholder { padding: 20px; }
  .results-toolbar { align-items: stretch; flex-direction: column; }
  .page-size-control, .pagination-control { justify-content: space-between; }
  .el-pagination { --el-pagination-button-width: 44px; }
  .file-item figcaption { padding: 14px 16px 16px; }
  .plot-skeleton:not(:first-child) { display: none; }
}
@media (prefers-reduced-motion: reduce) {
  .file-item { transition: none; }
}
</style>

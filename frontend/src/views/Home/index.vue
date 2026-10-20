<template>
  <div class="container_box">
    <section class="hero" aria-labelledby="hero-title">
      <div class="container hero-inner">
        <h1 id="hero-title" class="hero-title">Single-cell data for autoimmune and related immune-mediated conditions</h1>
        <p class="catalogue-context">Welcome to SCAID.</p>
        <p class="hero-lead">
          Browse annotated cell populations, precomputed gene and pathway maps and predicted cell–cell
          communication where published for a disease–tissue data object. Available expression matrices
          or pathway scores can be downloaded to analyse under your own design.
        </p>
        <form class="hero-search" role="search" @submit.prevent="handleSearch">
          <label class="search-label" for="home-search">Find a disease or dataset</label>
          <div class="search-controls">
          <el-input
            id="home-search"
            v-model="inputValue"
            size="large"
            clearable
            aria-label="Search a disease, abbreviation or accession"
            placeholder="Disease, abbreviation or accession, e.g. lupus, IBD, GSE174188"
          />
          <el-button type="primary" size="large" native-type="submit">Explore</el-button>
          </div>
        </form>
        <p class="hero-examples">
          Try:
          <router-link v-for="chip in examples" :key="chip.label" class="chip" :to="chip.to">{{ chip.label }}</router-link>
        </p>
      </div>
    </section>

    <section class="catalogue-summary container" aria-labelledby="catalogue-size-title" :aria-busy="statsLoading">
        <h2 id="catalogue-size-title">Database Statistics</h2>
        <p v-if="statsLoading" class="summary-status" role="status">Loading catalogue summary…</p>
        <p v-else-if="statsError" class="summary-status" role="status">Catalogue summary is temporarily unavailable. <button type="button" class="text-button" @click="fetchStats">Retry</button></p>
        <dl v-else-if="stats" class="stats" aria-label="Catalogue size, computed from the live catalogue">
          <div><dt>{{ fmt(stats.objects) }}</dt><dd>disease–tissue data objects</dd></div>
          <div><dt>{{ fmt(stats.sources) }}</dt><dd>source identifiers</dd></div>
          <div><dt>{{ fmt(stats.cells) }}</dt><dd>reported cell records</dd></div>
          <div><dt>{{ fmt(stats.samples) }}</dt><dd>reported sample labels</dd></div>
          <div><dt>{{ fmt(stats.conditions) }}</dt><dd>condition labels</dd></div>
          <div><dt>{{ fmt(stats.tissues) }}</dt><dd>tissue labels</dd></div>
        </dl>
        <p class="stats-note">Cell and sample totals sum available object-level values; cells are not deduplicated across objects. Sample labels are not verified donors, and source identifiers are not independent studies. See <router-link to="/Information">Information</router-link> for provenance and limits.</p>
    </section>

    <section class="browse container" aria-labelledby="browse-title">
      <div class="browse-heading">
        <h2 id="browse-title">Disease and condition coverage</h2>
        <p>Explore autoimmune, related immune-mediated, treatment-associated and other conditions included in SCAID. Select a condition below to browse its published data objects.</p>
      </div>
      <div class="coverage-layout">
        <figure class="condition-overview coverage-figure">
          <a :href="conditionOverview" class="overview-image-link" target="_blank" rel="noopener noreferrer" aria-label="Open the body illustration at full resolution in a new tab">
            <img :src="conditionOverview" width="1122" height="1402" alt="A soft, flat illustration of the human body and its organs. Complete condition names and dataset links are provided as text on this page." aria-describedby="overview-caption" loading="lazy" decoding="async" />
          </a>
          <figcaption id="overview-caption">
            <div class="overview-caption-heading">
              <a :href="conditionOverview" download="SCAID-body-overview-flat-2026-10-16.png">Download illustration PNG</a>
            </div>
          </figcaption>
        </figure>
        <div class="condition-index" :aria-busy="statsLoading">
          <div class="condition-index-heading">
            <h3>Browse all catalogued conditions</h3>
            <p>Each link opens the Dataset Browser filtered to that condition.</p>
            <p>The groups aid navigation. Conditions can affect several systems; a group does not identify a dataset’s sampled tissue or imply an exclusive affected organ. The linked index below follows the live catalogue.</p>
          </div>
          <p v-if="statsLoading" class="coverage-status" role="status">Loading condition coverage…</p>
          <el-alert v-else-if="statsError" title="Unable to load the condition index. Please try again." type="error" show-icon :closable="false">
            <el-button link type="primary" @click="fetchStats">Retry condition index</el-button>
          </el-alert>
          <div v-else-if="conditionGroups.length" class="condition-groups">
            <section v-for="group in conditionGroups" :key="group.key" class="condition-group" :aria-labelledby="`coverage-${group.key}`">
              <h4 :id="`coverage-${group.key}`">{{ group.title }}</h4>
              <ul>
                <li v-for="condition in group.conditions" :key="condition.abbreviation">
                  <router-link :to="{ path: '/Dataset', query: { id: condition.abbreviation } }" class="condition-link">
                    <span class="condition-name">{{ condition.label }}</span>
                    <span v-if="condition.label !== condition.abbreviation" class="condition-abbreviation">{{ condition.abbreviation }}</span>
                  </router-link>
                </li>
              </ul>
            </section>
          </div>
          <p v-else class="coverage-status" role="status">No conditions are currently listed in the published catalogue.</p>
        </div>
      </div>
    </section>
  </div>
</template>

<script setup>
defineOptions({ name: 'ScaidHome' });
import { computed, ref, onMounted, onBeforeUnmount } from 'vue';
import { getRequest } from '@/api/home';
import { useRouter } from 'vue-router';
import { buildConditionCoverage } from '@/utils/conditionCoverage';
import conditionOverview from '@/assets/images/scaid-body-overview-flat-20261016.png';
const router = useRouter();
const catalogueRows = ref([]);
const conditionGroups = computed(() => buildConditionCoverage(catalogueRows.value));
const inputValue = ref('');
const stats = ref(null);
const statsLoading = ref(true);
const statsError = ref(false);
let disposed = false;
const examples = [
  { label: 'SLE', to: { path: '/Dataset', query: { id: 'SLE' } } },
  { label: 'RA', to: { path: '/Dataset', query: { id: 'RA' } } },
  { label: 'IBD', to: { path: '/Dataset', query: { search: 'IBD' } } },
  { label: 'Psoriasis', to: { path: '/Dataset', query: { id: 'PSO' } } },
  { label: 'TRIM21', to: { path: '/Gene', query: { gene: 'TRIM21' } } },
];
const fmt = (n) => (n == null ? '—' : Number(n).toLocaleString('en-US'));

// Headline numbers are derived from the catalogue itself so they cannot drift from the data.
async function fetchStats() {
  statsLoading.value = true;
  statsError.value = false;
  try {
    const rows = await getRequest('/cell-data-all/', { format: 'json' });
    if (!Array.isArray(rows) || !rows.length) throw new Error('Invalid catalogue response');
    if (disposed) return;
    catalogueRows.value = rows;
    const uniq = (values) => new Set(values.filter(Boolean)).size;
    stats.value = {
      objects: rows.length,
      sources: uniq(rows.map(r => r.dataset_label || r.dataset_id)),
      cells: rows.reduce((sum, r) => sum + (Number(r.cell_count) || 0), 0),
      samples: rows.reduce((sum, r) => sum + (Number(r.sample_size) || 0), 0),
      conditions: uniq(rows.map(r => r.disease_label)),
      tissues: uniq(rows.map(r => (r.tissue_label || r.tissue || '').toLowerCase())),
    };
  } catch {
    stats.value = null;
    statsError.value = true;
  } finally {
    statsLoading.value = false;
  }
}

function handleSearch() {
  const search = inputValue.value.trim();
  router.push({ path: '/Dataset', query: search ? { search } : {} });
}
onMounted(fetchStats);
onBeforeUnmount(() => { disposed = true; });
</script>

<style scoped lang="scss">
.hero { background: var(--scaid-soft); border-bottom: 1px solid var(--scaid-line); padding: 48px 0 40px; }
.hero-inner { max-width: 920px; }
.hero-title { margin: 0 0 12px; max-width: 32ch; font-size: 2.25rem; line-height: 1.22; color: var(--scaid-ink); font-weight: 650; letter-spacing: -.025em; }
.catalogue-context { margin: 0 0 18px; color: var(--scaid-accent); font-size: .95rem; font-weight: 600; }
.hero-lead { margin: 0 0 28px; font-size: 1rem; line-height: 1.7; color: var(--scaid-muted); max-width: 73ch; }
.hero-search { max-width: 780px; }
.search-label { display: block; margin-bottom: 8px; color: var(--scaid-ink); font-size: .95rem; font-weight: 600; }
.search-controls { display: flex; gap: 10px; }
.search-controls .el-input { flex: 1; min-width: 0; }
.search-controls .el-button { min-width: 104px; }
.hero-examples { margin: 12px 0 0; font-size: .9rem; color: var(--scaid-muted); display: flex; align-items: center; flex-wrap: wrap; gap: 4px 8px; }
.chip { display: inline-flex; align-items: center; min-height: 36px; padding: 4px 12px; border: 1px solid #d8c6bb; border-radius: 999px; color: var(--scaid-accent); text-decoration: none; background: #fff; }
.chip:hover { background: #eaded5; border-color: #bd9a87; }
.catalogue-summary { padding-top: 28px; padding-bottom: 28px; border-bottom: 1px solid var(--scaid-line); }
.catalogue-summary h2 { margin: 0 0 20px; font-size: 1.125rem; font-weight: 600; }
.stats { display: grid; grid-template-columns: repeat(6, minmax(0, 1fr)); gap: 20px; margin: 0; }
.stats dt { font-size: 1.5rem; font-weight: 650; line-height: 1.3; font-variant-numeric: tabular-nums; color: var(--scaid-ink); }
.stats dd { margin: 5px 0 0; font-size: .875rem; color: var(--scaid-muted); }
.stats-note { margin: 18px 0 0; font-size: .875rem; color: var(--scaid-muted); }
.stats-note a { text-decoration: underline; }
.summary-status { color: var(--scaid-muted); margin: 0; }
.text-button { border: 0; background: transparent; color: var(--scaid-accent); text-decoration: underline; text-underline-offset: 3px; cursor: pointer; min-height: 44px; padding: 0 8px; }
.browse { padding-top: 36px; padding-bottom: 64px; }
.browse-heading h2 { margin: 0 0 8px; font-size: 1.5rem; font-weight: 650; color: var(--scaid-ink); }
.browse-heading p { max-width: 75ch; margin: 0 0 28px; color: var(--scaid-muted); }
/* Keep the illustration beside the live condition index, as in the original
   site. All names remain readable HTML rather than text inside the artwork. */
.coverage-layout { display: grid; grid-template-columns: minmax(0, 2fr) minmax(0, 3fr); gap: 36px; align-items: start; }
.condition-overview { min-width: 0; margin: 0; position: sticky; top: 100px; }
.overview-image-link { display: block; background: #fff; }
.overview-image-link img { display: block; width: 100%; height: auto; max-height: calc(100dvh - 220px); object-fit: contain; }
.overview-image-link:focus-visible { outline: 2px solid var(--scaid-focus); outline-offset: 4px; }
.condition-overview figcaption { padding-top: 14px; font-size: .875rem; color: var(--scaid-muted); }
.overview-caption-heading { display: flex; align-items: baseline; justify-content: space-between; flex-wrap: wrap; gap: 8px 24px; }
.overview-caption-heading a { color: var(--scaid-accent); text-decoration: underline; text-underline-offset: 3px; }
.condition-overview figcaption p { max-width: 75ch; margin: 8px 0 0; line-height: 1.65; }
.condition-index { min-width: 0; margin-top: 0; }
.condition-index-heading h3 { margin: 0 0 8px; font-size: 1.25rem; font-weight: 650; color: var(--scaid-ink); }
.condition-index-heading p { max-width: 75ch; margin: 0 0 12px; color: var(--scaid-muted); }
.condition-index-heading p:last-child { margin-bottom: 24px; }
.condition-groups { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 32px 24px; }
.condition-group { min-width: 0; padding-top: 16px; border-top: 1px solid var(--scaid-line); }
.condition-group h4 { margin: 0 0 10px; font-size: 1rem; font-weight: 650; line-height: 1.5; color: var(--scaid-ink); }
.condition-group ul { list-style: none; margin: 0; padding: 0; }
.condition-link { display: flex; align-items: baseline; justify-content: space-between; gap: 8px 12px; min-height: 44px; padding: 9px 0; color: var(--scaid-accent); text-decoration: none; line-height: 1.55; }
.condition-name { min-width: 0; font-size: .9375rem; text-decoration: underline; text-decoration-color: #dac3b6; text-underline-offset: 3px; overflow-wrap: anywhere; }
.condition-abbreviation { flex-shrink: 0; color: var(--scaid-muted); font-size: .8rem; }
.condition-link:hover .condition-name { color: #673924; text-decoration-color: currentColor; }
.coverage-status { margin: 0; color: var(--scaid-muted); }
@media (max-width: 991px) {
  .hero { padding: 36px 0; }
  .hero-title { font-size: 2rem; }
  .stats { grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 22px; }
  .coverage-layout { grid-template-columns: minmax(0, 1fr); gap: 32px; }
  .condition-overview { position: static; width: 100%; max-width: 500px; margin: 0 auto; }
  .overview-image-link img { max-height: none; }
}
@media (max-width: 600px) {
  .hero-title { font-size: 1.75rem; }
  .hero-lead { margin-bottom: 22px; }
  .hero { padding: 28px 0; }
  .search-controls { flex-direction: column; }
  .search-controls .el-button { width: 100%; }
  .chip { min-height: 44px; }
  .stats { grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 22px 16px; }
  .browse { padding-top: 28px; padding-bottom: 40px; }
  .overview-caption-heading { align-items: flex-start; flex-direction: column; gap: 4px; }
  .overview-caption-heading a { display: inline-flex; align-items: center; min-height: 44px; }
  .condition-groups { grid-template-columns: minmax(0, 1fr); gap: 24px; }
  .condition-link { gap: 12px; }
}
</style>

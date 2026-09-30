<template>
  <div class="container_box">
    <section class="hero" aria-labelledby="hero-title">
      <div class="container hero-inner">
        <h1 id="hero-title" class="hero-title">Public single-cell data for autoimmune and related immune-mediated conditions</h1>
        <p class="catalogue-context">One read-only catalogue. No account required.</p>
        <p class="hero-lead">
          Browse annotated cell populations, precomputed gene and pathway maps and predicted cell–cell
          communication for each disease–tissue data object, then download the expression matrix or
          pathway scores to analyse under your own design.
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
        <h2 id="catalogue-size-title">Catalogue at a glance</h2>
        <p v-if="statsLoading" class="summary-status" role="status">Loading catalogue summary…</p>
        <p v-else-if="statsError" class="summary-status" role="status">Catalogue summary is temporarily unavailable. <button type="button" class="text-button" @click="fetchStats">Retry</button></p>
        <dl v-else-if="stats" class="stats" aria-label="Catalogue size, computed from the live catalogue">
          <div><dt>{{ fmt(stats.objects) }}</dt><dd>disease–tissue data objects</dd></div>
          <div><dt>{{ fmt(stats.sources) }}</dt><dd>source identifiers</dd></div>
          <div><dt>{{ fmt(stats.cells) }}</dt><dd>cell records</dd></div>
          <div><dt>{{ fmt(stats.samples) }}</dt><dd>sample labels</dd></div>
          <div><dt>{{ fmt(stats.conditions) }}</dt><dd>condition labels</dd></div>
          <div><dt>{{ fmt(stats.tissues) }}</dt><dd>tissue labels</dd></div>
        </dl>
        <p class="stats-note">Sample labels are not verified donors, and source identifiers are not independent studies. See <router-link to="/Information">Information</router-link> for provenance and limits.</p>
    </section>

    <section class="browse container" aria-labelledby="browse-title">
      <div class="browse-heading">
        <h2 id="browse-title">Browse by affected site</h2>
        <p>Each link opens the Dataset page filtered to that condition. Diseases listed under several sites have data objects from more than one tissue.</p>
      </div>
      <p v-if="loading" role="status">Loading diseases…</p>
      <el-alert v-if="loadError" :title="loadError" type="error" show-icon :closable="false">
        <el-button link type="primary" @click="fetchData">Retry</el-button>
      </el-alert>
      <div class="browse-grid">
        <figure class="body-map">
          <img src="../../assets/images/aid.jpg" alt="Schematic of the human body marking the tissues represented in SCAID" loading="lazy" decoding="async" />
        </figure>
        <div class="site-cards">
          <article v-for="(item, index) in bloodList" :key="item.id || index" class="tissue-card">
            <h3 class="tissue-title">{{ item.name_en }}</h3>
            <ul class="disease-list">
              <li v-for="(option, optionIndex) in item.diseases || []" :key="option.abbreviation || optionIndex">
                <router-link v-if="option.abbreviation" :to="{ path: '/Dataset', query: { id: option.abbreviation } }" class="disease-name">{{ option.name_en }}</router-link>
                <span v-else class="disease-name">{{ option.name_en }}</span>
              </li>
              <li v-if="!item.diseases?.length" class="empty-tissue">No datasets available</li>
            </ul>
          </article>
        </div>
      </div>
    </section>
  </div>
</template>

<script setup>
defineOptions({ name: 'ScaidHome' });
import { ref, onMounted, onBeforeUnmount } from 'vue';
import { getRequest } from '@/api/home';
import { useRouter } from 'vue-router';
const router = useRouter();
const bloodList = ref([]);
const inputValue = ref('');
const loading = ref(false);
const loadError = ref('');
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

async function fetchData() {
  loading.value = true;
  loadError.value = '';
  try {
    const response = await getRequest('/by-site/', { format: 'json' });
    if (response.code !== 200 || !Array.isArray(response.data)) throw new Error('Invalid disease response');
    if (!disposed) bloodList.value = response.data;
  } catch {
    loadError.value = 'Unable to load disease categories. Please try again.';
  } finally {
    loading.value = false;
  }
}

// Headline numbers are derived from the catalogue itself so they cannot drift from the data.
async function fetchStats() {
  statsLoading.value = true;
  statsError.value = false;
  try {
    const rows = await getRequest('/cell-data-all/', { format: 'json' });
    if (!Array.isArray(rows) || !rows.length) throw new Error('Invalid catalogue response');
    if (disposed) return;
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
onMounted(() => { fetchData(); fetchStats(); });
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
.browse-grid { display: grid; grid-template-columns: minmax(240px, 2fr) 3fr; gap: 36px; align-items: start; }
.body-map { margin: 0; position: sticky; top: 100px; }
.body-map img { width: 100%; height: auto; aspect-ratio: 1100 / 633; display: block; }
.site-cards { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 28px 24px; }
.tissue-card { border-top: 1px solid var(--scaid-line); padding-top: 14px; min-width: 0; }
.tissue-title { margin: 0 0 8px; font-size: 1rem; font-weight: 650; color: var(--scaid-ink); }
.disease-list { list-style: none; margin: 0; padding: 0; }
.disease-list li + li { margin-top: 4px; }
.disease-name { display: block; padding: 6px 0; font-size: .95rem; color: var(--scaid-accent); text-decoration: underline; text-decoration-color: #dac3b6; text-underline-offset: 3px; overflow-wrap: anywhere; }
a.disease-name:hover { color: #673924; text-decoration-color: currentColor; }
.empty-tissue { font-size: .9rem; color: var(--scaid-muted); }
@media (max-width: 991px) {
  .hero { padding: 36px 0; }
  .hero-title { font-size: 2rem; }
  .stats { grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 22px; }
  .browse-grid { grid-template-columns: 1fr; gap: 28px; }
  .body-map { position: static; max-width: 560px; margin: 0 auto; }
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
  .site-cards { gap: 24px 18px; }
  .disease-name { min-height: 44px; }
}
</style>

<template>
  <div class="container_box">
    <section class="hero" aria-labelledby="hero-title">
      <div class="container hero-inner">
        <div class="hero-copy">
          <h1 id="hero-title" class="hero-title">Single-cell data for autoimmune and related immune-mediated conditions</h1>
          <p class="catalogue-context">Welcome to SCAID.</p>
          <p class="hero-lead">
            Browse annotated cell populations, precomputed gene and pathway maps and predicted cell–cell
            communication where published for a disease–tissue data object. Available expression matrices
            or pathway scores can be downloaded to analyse under your own design.
          </p>
        </div>
        <div class="hero-explore">
          <form class="hero-search" role="search" @submit.prevent="handleSearch">
            <label class="search-label" for="home-search">Find a disease or dataset</label>
            <div class="search-controls">
              <el-input
                id="home-search"
                v-model="inputValue"
                size="large"
                clearable
                aria-label="Search a disease, abbreviation or accession"
                placeholder="Disease, abbreviation or accession"
              />
              <el-button type="primary" size="large" native-type="submit">
                Explore
                <svg class="explore-arrow" viewBox="0 0 20 20" fill="none" aria-hidden="true"><path d="M4 10h12m-5-5 5 5-5 5" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" /></svg>
              </el-button>
            </div>
          </form>
          <p class="hero-examples">
            <span>Try:</span>
            <router-link v-for="chip in examples" :key="chip.label" class="chip" :to="chip.to">{{ chip.label }}</router-link>
          </p>
        </div>
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
            <details class="coverage-notes">
              <summary>
                <svg viewBox="0 0 16 16" fill="none" aria-hidden="true"><path d="m6 4 4 4-4 4" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" /></svg>
                About this condition index
              </summary>
              <p>The groups aid navigation. Conditions can affect several systems; a group does not identify a dataset’s sampled tissue or imply an exclusive affected organ. The linked index below follows the live catalogue.</p>
            </details>
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
import { summarizeCatalogue } from '@/utils/catalogueSummary';
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
    stats.value = summarizeCatalogue(rows);
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
.hero {
  padding: 56px 0;
  background: var(--scaid-soft);
  border-bottom: 1px solid var(--scaid-line);
}
.hero-inner {
  display: grid;
  grid-template-columns: minmax(0, 1.4fr) minmax(320px, .9fr);
  gap: 64px;
  align-items: center;
}
.hero-copy, .hero-explore { min-width: 0; }
.hero-title {
  max-width: 29ch;
  margin: 0 0 18px;
  color: var(--scaid-ink);
  font-size: clamp(2rem, 2.65vw, 2.625rem);
  font-weight: 600;
  line-height: 1.15;
  letter-spacing: -.03em;
}
.catalogue-context {
  margin: 0 0 10px;
  color: var(--scaid-accent);
  font-size: .9375rem;
  font-weight: 600;
}
.hero-lead {
  max-width: 67ch;
  margin: 0;
  color: var(--scaid-muted);
  font-size: .9375rem;
  line-height: 1.7;
}
.hero-explore {
  padding-left: 32px;
  border-left: 1px solid var(--scaid-line);
}
.search-label {
  display: block;
  margin-bottom: 16px;
  color: var(--scaid-ink);
  font-size: 1.125rem;
  font-weight: 600;
  line-height: 1.4;
}
.search-controls { display: grid; gap: 12px; }
.search-controls .el-input { min-width: 0; }
.search-controls .el-button { width: 100%; }
.explore-arrow { width: 20px; height: 20px; margin-left: 12px; }
.hero-examples {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 4px 7px;
  margin: 16px 0 0;
  color: var(--scaid-muted);
  font-size: .8125rem;
}
.hero-examples > span { margin-right: 3px; }
.chip {
  display: inline-flex;
  align-items: center;
  min-height: 40px;
  padding: 4px 11px;
  border: 1px solid var(--scaid-line);
  border-radius: 999px;
  background: #fff;
  color: var(--scaid-accent);
  text-decoration: none;
}
.chip:hover { border-color: var(--scaid-accent); background: var(--scaid-soft); }
.catalogue-summary {
  padding-top: 30px;
  padding-bottom: 28px;
  border-bottom: 1px solid var(--scaid-line);
}
.catalogue-summary h2 {
  margin: 0 0 22px;
  color: var(--scaid-ink);
  font-size: 1.125rem;
  font-weight: 600;
}
.stats {
  display: grid;
  grid-template-columns: repeat(6, minmax(0, 1fr));
  margin: 0;
}
.stats > div { min-width: 0; padding: 0 20px; border-left: 1px solid var(--scaid-line); }
.stats > div:first-child { padding-left: 0; border-left: 0; }
.stats > div:last-child { padding-right: 0; }
.stats dt {
  color: var(--scaid-ink);
  font-size: clamp(1.5rem, 2.1vw, 1.875rem);
  font-weight: 500;
  line-height: 1.25;
  letter-spacing: -.025em;
  font-variant-numeric: tabular-nums;
}
.stats dd { margin: 7px 0 0; color: var(--scaid-muted); font-size: .8125rem; line-height: 1.5; }
.stats-note {
  margin: 22px 0 0;
  color: var(--scaid-muted);
  font-size: .8125rem;
  line-height: 1.65;
}
.stats-note a { text-decoration: underline; }
.summary-status { margin: 0; color: var(--scaid-muted); }
.text-button {
  min-height: 44px;
  padding: 0 8px;
  border: 0;
  background: transparent;
  color: var(--scaid-accent);
  text-decoration: underline;
  text-underline-offset: 3px;
  cursor: pointer;
}
.browse { padding-top: 44px; padding-bottom: 64px; }
.browse-heading h2 {
  margin: 0 0 10px;
  color: var(--scaid-ink);
  font-size: 1.625rem;
  font-weight: 600;
  letter-spacing: -.02em;
}
.browse-heading p { max-width: 75ch; margin: 0 0 32px; color: var(--scaid-muted); line-height: 1.7; }
/* Preserve the original left-hand artwork and the linked, live condition index. */
.coverage-layout {
  display: grid;
  grid-template-columns: minmax(0, 2fr) minmax(0, 3fr);
  gap: 48px;
  align-items: start;
}
.condition-overview { position: sticky; top: 100px; min-width: 0; margin: 0; }
.overview-image-link { display: block; background: #fff; }
.overview-image-link img {
  display: block;
  width: 100%;
  height: auto;
  max-height: calc(100dvh - 220px);
  object-fit: contain;
}
.overview-image-link:focus-visible { outline: 2px solid var(--scaid-focus); outline-offset: 4px; }
.condition-overview figcaption { padding-top: 14px; color: var(--scaid-muted); font-size: .8125rem; }
.overview-caption-heading { display: flex; flex-wrap: wrap; align-items: baseline; gap: 8px 24px; }
.overview-caption-heading a {
  display: inline-flex;
  align-items: center;
  min-height: 44px;
  color: var(--scaid-accent);
  text-decoration: underline;
  text-underline-offset: 3px;
}
.condition-index { min-width: 0; }
.condition-index-heading { margin-bottom: 24px; }
.condition-index-heading h3 {
  margin: 0 0 7px;
  color: var(--scaid-ink);
  font-size: 1.125rem;
  font-weight: 600;
}
.condition-index-heading p { max-width: 75ch; margin: 0 0 6px; color: var(--scaid-muted); font-size: .875rem; line-height: 1.65; }
.coverage-notes { color: var(--scaid-muted); font-size: .8125rem; }
.coverage-notes summary {
  display: flex;
  align-items: center;
  gap: 5px;
  width: fit-content;
  min-height: 36px;
  color: var(--scaid-accent);
  list-style: none;
  cursor: pointer;
}
.coverage-notes summary::-webkit-details-marker { display: none; }
.coverage-notes summary svg { flex-shrink: 0; width: 16px; height: 16px; }
.coverage-notes[open] summary svg { transform: rotate(90deg); }
.coverage-notes summary:hover { text-decoration: underline; text-underline-offset: 3px; }
.coverage-notes summary:focus-visible { outline: 2px solid var(--scaid-focus); outline-offset: 3px; }
.coverage-notes p { margin: 4px 0 12px; padding-left: 21px; font-size: .8125rem; }
.condition-groups { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 28px 28px; }
.condition-group { min-width: 0; padding-top: 16px; border-top: 1px solid var(--scaid-line); }
.condition-group h4 { margin: 0 0 8px; color: var(--scaid-ink); font-size: .9375rem; font-weight: 600; line-height: 1.5; }
.condition-group ul { margin: 0; padding: 0; list-style: none; }
.condition-link {
  display: flex;
  justify-content: space-between;
  align-items: baseline;
  gap: 8px 12px;
  min-height: 44px;
  margin: 0 -7px;
  padding: 9px 7px;
  border-radius: 6px;
  color: var(--scaid-accent);
  text-decoration: none;
  line-height: 1.55;
}
.condition-name { min-width: 0; font-size: .875rem; overflow-wrap: anywhere; }
.condition-abbreviation { flex-shrink: 0; color: var(--scaid-muted); font-size: .75rem; }
.condition-link:hover { background: var(--scaid-soft); }
.condition-link:hover .condition-name { text-decoration: underline; text-underline-offset: 3px; }
.coverage-status { margin: 0; color: var(--scaid-muted); }
@media (max-width: 1199px) {
  .hero-inner { grid-template-columns: minmax(0, 1.3fr) minmax(280px, .9fr); gap: 36px; }
  .hero-explore { padding-left: 24px; }
  .stats > div { padding: 0 16px; }
  .coverage-layout { gap: 32px; }
  .condition-groups { column-gap: 24px; }
}
@media (max-width: 991px) {
  .hero { padding: 40px 0; }
  .hero-inner { grid-template-columns: minmax(0, 1fr); gap: 28px; }
  .hero-title { max-width: 31ch; font-size: 2.25rem; }
  .hero-explore { padding: 24px 0 0; border-top: 1px solid var(--scaid-line); border-left: 0; }
  .search-label { margin-bottom: 12px; }
  .search-controls { grid-template-columns: minmax(0, 1fr) auto; }
  .search-controls .el-button { min-width: 130px; }
  .hero-examples { margin-top: 12px; }
  .stats { grid-template-columns: repeat(3, minmax(0, 1fr)); row-gap: 26px; }
  .stats > div:nth-child(3n + 1) { padding-left: 0; border-left: 0; }
  .stats dt { font-size: 1.75rem; }
  .coverage-layout { grid-template-columns: minmax(0, 1fr); gap: 32px; }
  .condition-overview { position: static; width: 100%; max-width: 440px; margin: 0 auto; }
  .overview-image-link img { max-height: none; }
}
@media (max-width: 600px) {
  .hero { padding: 30px 0; }
  .hero-title { margin-bottom: 16px; font-size: 2rem; line-height: 1.17; }
  .hero-lead { font-size: .9375rem; }
  .hero-inner { gap: 24px; }
  .hero-explore { padding-top: 22px; }
  .search-controls { grid-template-columns: minmax(0, 1fr); }
  .search-controls .el-button { min-width: 0; }
  .chip { min-height: 44px; }
  .catalogue-summary { padding-top: 26px; padding-bottom: 26px; }
  .catalogue-summary h2 { margin-bottom: 20px; }
  .stats { grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 22px 0; }
  .stats > div, .stats > div:nth-child(3n + 1) { padding: 0 18px; border-left: 1px solid var(--scaid-line); }
  .stats > div:nth-child(2n + 1) { padding-left: 0; border-left: 0; }
  .stats dt { font-size: 1.75rem; }
  .stats-note { margin-top: 20px; }
  .browse { padding-top: 32px; padding-bottom: 40px; }
  .browse-heading h2 { font-size: 1.5rem; line-height: 1.25; }
  .browse-heading p { margin-bottom: 24px; font-size: .9375rem; }
  .coverage-layout { gap: 24px; }
  .condition-groups { grid-template-columns: minmax(0, 1fr); gap: 26px; }
  .condition-link { gap: 12px; }
  .condition-name { font-size: .9375rem; }
  .condition-abbreviation { font-size: .8125rem; }
  .coverage-notes summary { min-height: 44px; }
}
</style>

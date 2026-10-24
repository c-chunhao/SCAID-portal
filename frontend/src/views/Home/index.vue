<template>
  <div class="container_box">
    <section class="hero" aria-labelledby="hero-title">
      <div class="container hero-inner">
        <div class="hero-copy">
          <h1 id="hero-title" class="hero-title">
            <span class="hero-brand">SCAID</span>
            <span class="hero-subtitle">Single-cell data for autoimmune and related immune-mediated conditions</span>
          </h1>
          <p class="hero-lead">
            Welcome to SCAID. Browse annotated cell populations, precomputed gene and pathway maps and predicted cell–cell
            communication where published for a disease–tissue data object. Available expression matrices
            or pathway scores can be downloaded to analyse under your own design.
          </p>
        </div>
        <div class="hero-explore">
          <router-link to="/TeamInformation" class="atlas-preview" aria-label="Explore the T and NK cell atlas">
            <img :src="atlasPreview" width="1600" height="863" alt="Published T and NK cell atlas showing cell-type annotations across UMAP panels." decoding="async" />
            <span class="atlas-preview-caption">
              <span>T / NK cell atlas</span>
              <span class="atlas-preview-action">Explore atlas <svg viewBox="0 0 20 20" fill="none" aria-hidden="true"><path d="M4 10h12m-5-5 5 5-5 5" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" /></svg></span>
            </span>
          </router-link>
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

    <section class="catalogue-summary" aria-labelledby="catalogue-size-title" :aria-busy="statsLoading">
      <div class="container">
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
      </div>
    </section>

    <section class="browse container" aria-labelledby="browse-title">
      <div class="browse-heading">
        <h2 id="browse-title">Disease and condition coverage</h2>
        <p>Explore autoimmune, related immune-mediated, treatment-associated and other conditions included in SCAID. Select a condition below to browse its published data objects.</p>
      </div>
      <nav v-if="conditionGroups.length" class="condition-jump" aria-label="Jump to a condition group">
        <span class="condition-jump-label">Jump to a group</span>
        <div class="condition-jump-links">
          <a v-for="group in conditionGroups" :key="group.key" :href="`#coverage-${group.key}`" class="condition-jump-link">
            {{ group.title }}<span class="condition-jump-count">{{ group.conditions.length }}</span>
          </a>
        </div>
      </nav>
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
              <h4 :id="`coverage-${group.key}`" tabindex="-1">{{ group.title }}</h4>
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
import atlasPreview from '@/assets/images/dataIntegration_preview/FIG1_T_NK_celltype.jpg';
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
  padding: 48px 0 52px;
  background: var(--scaid-soft);
}
.hero-inner {
  display: grid;
  grid-template-columns: minmax(0, 1.3fr) minmax(340px, .9fr);
  gap: 88px;
  align-items: center;
}
.hero-copy, .hero-explore { min-width: 0; }
.hero-title {
  margin: 0 0 22px;
  color: var(--scaid-ink);
  font-weight: 600;
  letter-spacing: -.03em;
}
.hero-brand {
  display: block;
  margin-bottom: 20px;
  color: var(--scaid-accent);
  font-size: clamp(4rem, 7.1vw, 6rem);
  line-height: .94;
}
.hero-subtitle {
  display: block;
  max-width: 30ch;
  font-size: clamp(1.625rem, 2.35vw, 2rem);
  line-height: 1.2;
}
.hero-lead {
  max-width: 66ch;
  margin: 0;
  color: var(--scaid-muted);
  font-size: .9375rem;
  line-height: 1.7;
}
.hero-explore {
  padding-left: 36px;
  border-left: 1px solid var(--scaid-line);
}
.atlas-preview { display: block; margin-bottom: 24px; color: var(--scaid-ink); text-decoration: none; }
.atlas-preview img { display: block; width: 100%; height: auto; max-height: 180px; object-fit: contain; background: #fff; }
.atlas-preview-caption { display: flex; align-items: center; justify-content: space-between; gap: 12px; padding-top: 10px; font-size: .8125rem; font-weight: 500; }
.atlas-preview-action { display: inline-flex; align-items: center; gap: 6px; color: var(--scaid-accent); white-space: nowrap; }
.atlas-preview-action svg { width: 17px; height: 17px; }
.atlas-preview:hover .atlas-preview-action { color: var(--scaid-accent-hover); text-decoration: underline; text-underline-offset: 3px; }
.search-label {
  display: block;
  margin-bottom: 12px;
  color: var(--scaid-ink);
  font-size: 1.125rem;
  font-weight: 600;
  line-height: 1.4;
}
.search-controls { display: grid; grid-template-columns: minmax(0, 1fr) auto; gap: 8px; }
.search-controls .el-input { min-width: 0; }
.search-controls .el-button { padding-right: 16px; padding-left: 16px; }
.explore-arrow { width: 18px; height: 18px; margin-left: 8px; }
.hero-examples {
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 2px 3px;
  margin: 10px 0 0;
  color: var(--scaid-muted);
  font-size: .8125rem;
}
.hero-examples > span { margin-right: 3px; }
.chip {
  display: inline-flex;
  align-items: center;
  min-height: 44px;
  padding: 4px 6px;
  color: var(--scaid-accent);
  text-decoration: underline;
  text-decoration-color: var(--scaid-line-strong);
  text-underline-offset: 4px;
}
.chip:hover { color: var(--scaid-accent-hover); text-decoration-color: currentColor; }
.catalogue-summary {
  padding-top: 30px;
  padding-bottom: 30px;
  background: var(--scaid-accent-soft);
}
.catalogue-summary h2 {
  margin: 0 0 24px;
  color: var(--scaid-ink);
  font-size: 1.125rem;
  font-weight: 600;
}
.stats {
  display: grid;
  grid-template-columns: repeat(6, minmax(0, 1fr));
  margin: 0;
}
.stats > div { min-width: 0; padding: 0 20px; border-left: 1px solid var(--scaid-line-strong); }
.stats > div:first-child { padding-left: 0; border-left: 0; }
.stats > div:last-child { padding-right: 0; }
.stats dt {
  color: var(--scaid-ink);
  font-size: clamp(1.5rem, 2.5vw, 2rem);
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
.browse { padding-top: 56px; padding-bottom: 72px; }
.browse-heading { display: grid; grid-template-columns: minmax(0, 2fr) minmax(0, 3fr); gap: 48px; align-items: start; }
.browse-heading h2 {
  max-width: 20ch;
  margin: 0;
  color: var(--scaid-ink);
  font-size: 1.875rem;
  font-weight: 600;
  letter-spacing: -.02em;
}
.browse-heading p { max-width: 66ch; margin: 0; color: var(--scaid-muted); line-height: 1.7; }
.condition-jump { margin: 28px 0 32px; padding: 20px 0; border-top: 1px solid var(--scaid-line); border-bottom: 1px solid var(--scaid-line); }
.condition-jump-label { display: block; margin-bottom: 12px; color: var(--scaid-muted); font-size: .8125rem; }
.condition-jump-links { display: flex; flex-wrap: wrap; gap: 8px; }
.condition-jump-link { display: inline-flex; align-items: center; gap: 10px; min-height: 40px; padding: 8px 13px; border: 1px solid var(--scaid-line); border-radius: 999px; color: var(--scaid-accent); font-size: .8125rem; line-height: 1.35; text-decoration: none; white-space: nowrap; }
.condition-jump-count { color: var(--scaid-muted); font-size: .75rem; font-variant-numeric: tabular-nums; }
.condition-jump-link:hover { border-color: var(--scaid-accent); background: var(--scaid-accent-soft); color: var(--scaid-accent-hover); }
/* Preserve the original left-hand artwork and the linked, live condition index. */
.coverage-layout {
  display: grid;
  grid-template-columns: minmax(0, 2fr) minmax(0, 3fr);
  gap: 48px;
  align-items: start;
}
.condition-overview { position: sticky; top: 100px; min-width: 0; margin: 0; padding: 24px 24px 12px; border: 1px solid var(--scaid-line); border-radius: var(--scaid-radius); background: #fff; }
.overview-image-link { display: block; background: #fff; }
.overview-image-link img {
  display: block;
  width: 100%;
  height: auto;
  max-height: calc(100dvh - 232px);
  object-fit: contain;
}
.overview-image-link:focus-visible { outline: 2px solid var(--scaid-focus); outline-offset: 4px; }
.condition-overview figcaption { margin-top: 14px; padding-top: 6px; border-top: 1px solid var(--scaid-line); color: var(--scaid-muted); font-size: .8125rem; }
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
.condition-index-heading { margin-bottom: 20px; }
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
.condition-groups { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 28px; }
.condition-group { min-width: 0; padding-top: 16px; border-top: 1px solid var(--scaid-line); }
.condition-group h4 { margin: 0 0 8px; color: var(--scaid-ink); font-size: .9375rem; font-weight: 600; line-height: 1.5; scroll-margin-top: 24px; }
.condition-group h4:focus-visible { outline: 2px solid var(--scaid-focus); outline-offset: 5px; }
.condition-group h4:target { color: var(--scaid-accent); }
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
  .coverage-layout, .browse-heading { gap: 32px; }
  .condition-groups { column-gap: 24px; }
}
@media (max-width: 991px) {
  .hero { padding: 36px 0; }
  .hero-inner { grid-template-columns: minmax(0, 1fr); gap: 28px; }
  .hero-title { max-width: none; }
  .hero-brand { font-size: 4.5rem; margin-bottom: 16px; }
  .hero-subtitle { max-width: 36ch; font-size: 1.875rem; }
  .atlas-preview { display: none; }
  .hero-explore { padding: 24px 0 0; border-top: 1px solid var(--scaid-line); border-left: 0; }
  .search-label { margin-bottom: 12px; }
  .search-controls .el-button { min-width: 130px; }
  .hero-examples { margin-top: 6px; }
  .stats { grid-template-columns: repeat(3, minmax(0, 1fr)); row-gap: 26px; }
  .stats > div:nth-child(3n + 1) { padding-left: 0; border-left: 0; }
  .stats dt { font-size: 1.75rem; }
  .browse-heading { grid-template-columns: minmax(0, 1fr); gap: 12px; }
  .browse-heading h2 { max-width: none; }
  .condition-jump { margin: 24px 0 28px; padding: 16px 0 10px; }
  .condition-jump-label { margin-bottom: 10px; }
  .condition-jump-links { flex-wrap: nowrap; overflow-x: auto; padding: 4px 4px 10px; scroll-padding-inline: 4px; }
  .condition-jump-link { flex-shrink: 0; min-height: 44px; }
  .coverage-layout { grid-template-columns: minmax(0, 1fr); gap: 32px; }
  .condition-overview { position: static; width: 100%; max-width: 440px; margin: 0 auto; }
  .overview-image-link img { max-height: none; }
}
@media (max-width: 600px) {
  .hero { padding: 26px 0 22px; }
  .hero-title { margin-bottom: 16px; }
  .hero-brand { margin-bottom: 14px; font-size: 3.5rem; }
  .hero-subtitle { max-width: 28ch; font-size: 1.5rem; line-height: 1.2; }
  .hero-lead { font-size: .9375rem; line-height: 1.65; }
  .hero-inner { gap: 20px; }
  .hero-explore { padding-top: 18px; }
  .search-label { font-size: 1rem; margin-bottom: 10px; }
  .search-controls .el-button { min-width: 0; padding-right: 12px; padding-left: 12px; }
  .explore-arrow { width: 16px; height: 16px; margin-left: 6px; }
  .hero-examples { column-gap: 1px; }
  .chip { padding-right: 5px; padding-left: 5px; }
  .catalogue-summary { padding-top: 26px; padding-bottom: 26px; }
  .catalogue-summary h2 { margin-bottom: 20px; }
  .stats { grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 22px 0; }
  .stats > div, .stats > div:nth-child(3n + 1) { padding: 0 18px; border-left: 1px solid var(--scaid-line-strong); }
  .stats > div:nth-child(2n + 1) { padding-left: 0; border-left: 0; }
  .stats dt { font-size: 1.75rem; }
  .stats-note { margin-top: 20px; }
  .browse { padding-top: 36px; padding-bottom: 40px; }
  .browse-heading h2 { font-size: 1.5rem; line-height: 1.25; }
  .browse-heading p { font-size: .9375rem; }
  .condition-jump { margin: 20px 0 24px; }
  .coverage-layout { gap: 24px; }
  .condition-overview { max-width: 300px; padding: 18px 18px 8px; }
  .condition-groups { grid-template-columns: minmax(0, 1fr); gap: 26px; }
  .condition-link { gap: 12px; }
  .condition-name { font-size: .9375rem; }
  .condition-abbreviation { font-size: .8125rem; }
  .coverage-notes summary { min-height: 44px; }
}
</style>

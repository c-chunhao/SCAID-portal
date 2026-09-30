<template>
  <ExplorationSearch kind="kegg">
    Explore <a href="https://www.kegg.jp/kegg/pathway.html" target="_blank" rel="noopener noreferrer">KEGG pathways</a>
    across autoimmune disease datasets. Select a pathway ID (for example,
    hsa00010 for Glycolysis / Gluconeogenesis) and an optional scoring method to
    view the precomputed gene-set score maps (AUCell, UCell, singscore). Scores summarize expression relative to a gene set and do not measure biochemical pathway activity.
  </ExplorationSearch>
</template>

<script setup>
import ExplorationSearch from '@/components/ExplorationSearch.vue'
defineOptions({ name: 'KeggExploration' })
</script>

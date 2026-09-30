<template>
  <ExplorationSearch kind="gene">
    Search for a single gene across autoimmune disease datasets and explore its
    expression in precomputed UMAP plots. Select a disease or dataset to narrow
    the search, or leave those filters empty to compare all available results.
  </ExplorationSearch>
</template>

<script setup>
import ExplorationSearch from '@/components/ExplorationSearch.vue'
defineOptions({ name: 'GeneExploration' })
</script>

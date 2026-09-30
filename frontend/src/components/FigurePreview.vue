<template>
  <div class="figure-preview" :class="{ 'square-preview': square }">
    <a v-if="!failed" :href="fullSrc" target="_blank" rel="noopener" :aria-label="fullLabel || `Open full-resolution ${alt} in a new tab`" title="Open the full-resolution figure in a new tab">
      <img :key="attempt" :src="previewSrc" :alt="alt" :width="width || undefined" :height="height || undefined" loading="lazy" decoding="async" @error="failed = true" />
    </a>
    <div v-else class="preview-recovery" role="status">
      <p>Figure preview could not be loaded.</p>
      <button type="button" @click="retry">Retry figure</button>
      <a :href="fullSrc" target="_blank" rel="noopener">Open full-resolution figure</a>
    </div>
  </div>
</template>

<script setup>
import { computed, ref, watch } from 'vue'
const props = defineProps({
  src: { type: String, required: true }, fullSrc: { type: String, required: true },
  alt: { type: String, required: true }, fullLabel: { type: String, default: '' },
  square: Boolean, width: Number, height: Number
})
const failed = ref(false)
const attempt = ref(0)
const previewSrc = computed(() => {
  if (!attempt.value) return props.src
  const url = new URL(props.src, window.location.origin)
  url.searchParams.set('preview_retry', String(attempt.value))
  return url.href
})
const retry = () => { attempt.value += 1; failed.value = false }
watch(() => props.src, () => { failed.value = false; attempt.value = 0 })
</script>

<style scoped>
.figure-preview { min-width: 0; }
.figure-preview > a { display: block; border-radius: 10px; }
.figure-preview img { display: block; width: 100%; height: auto; max-width: 100%; background: #fff; }
.square-preview { aspect-ratio: 1; }
.square-preview > a { display: flex; height: 100%; align-items: center; justify-content: center; padding: .5rem; }
.square-preview img { height: 100%; object-fit: contain; }
.figure-preview a:focus-visible, .preview-recovery button:focus-visible { outline: 3px solid var(--scaid-focus); outline-offset: 3px; }
.preview-recovery { display: flex; min-height: 180px; height: 100%; flex-direction: column; justify-content: center; align-items: center; gap: .75rem; padding: 1.25rem; text-align: center; }
.preview-recovery p { margin: 0; color: var(--scaid-muted); }
.preview-recovery button { min-height: 44px; padding: .5rem .8rem; border: 1px solid var(--scaid-line); border-radius: 6px; background: #fff; color: var(--scaid-accent); font: inherit; cursor: pointer; }
.preview-recovery a { color: var(--scaid-accent); text-decoration: underline; text-underline-offset: .15em; }
</style>

<template>
  <div class="figure-preview" :class="{ 'square-preview': square }">
    <div class="preview-media" :aria-busy="loading">
      <button v-if="!failed" type="button" class="preview-image" :aria-label="`View larger: ${alt}`" @click="openViewer($event)">
        <img :key="previewSrc" :src="previewSrc" :alt="alt" :width="width || undefined" :height="height || undefined" :style="imageStyle" loading="lazy" decoding="async" @load="loading = false" @error="onError" />
      </button>
      <span v-if="loading && !failed" class="preview-loading" role="status">Loading figure…</span>
      <div v-if="failed" class="preview-recovery" role="status">
        <p>Figure preview could not be loaded.</p>
        <button type="button" @click="retry">Retry figure</button>
      </div>
    </div>
    <div class="preview-toolbar" role="group" :aria-label="`Viewing options: ${alt}`">
      <button type="button" @click="openViewer($event)">View larger</button>
      <button type="button" @click="openViewer($event, true)">Full screen</button>
      <a :href="fullSrc" target="_blank" rel="noopener" :aria-label="fullLabel || `Open original ${alt} in a new tab`">Open original</a>
    </div>
    <FigureViewer ref="viewer" :src="fullSrc" :alt="alt" :full-label="fullLabel" @imageerror="emit('imageerror', $event)" />
  </div>
</template>

<script setup>
import { computed, ref, watch } from 'vue'
import FigureViewer from './FigureViewer.vue'
const props = defineProps({
  src: { type: String, required: true }, fullSrc: { type: String, required: true },
  alt: { type: String, required: true }, fullLabel: { type: String, default: '' },
  square: Boolean, width: Number, height: Number,
  imageStyle: { type: [String, Object, Array], default: undefined }
})
const emit = defineEmits(['imageerror'])
const viewer = ref(null)
const failed = ref(false)
const loading = ref(true)
const attempt = ref(0)
const previewSrc = computed(() => {
  if (!attempt.value) return props.src
  const url = new URL(props.src, window.location.origin)
  url.searchParams.set('preview_retry', String(attempt.value))
  return url.href
})
const retry = () => { attempt.value += 1; failed.value = false; loading.value = true }
const onError = event => {
  failed.value = true
  loading.value = false
  emit('imageerror', { source: 'preview', src: props.src, event })
}
const openViewer = (event, fullscreen = false) => viewer.value?.open({ fullscreen, returnFocus: event.currentTarget })
watch(() => props.src, () => { failed.value = false; loading.value = true; attempt.value = 0 })
</script>

<style scoped>
.figure-preview { min-width: 0; }
.preview-media { position: relative; min-width: 0; background: #fff; }
.preview-image { display: flex; width: 100%; min-height: 144px; align-items: center; justify-content: center; padding: 0; border: 0; border-radius: 12px 12px 0 0; background: #fff; cursor: zoom-in; }
.figure-preview img { display: block; width: 100%; height: auto; max-width: 100%; background: #fff; object-fit: contain; }
.square-preview .preview-media { aspect-ratio: 1; }
.square-preview .preview-image, .square-preview .preview-recovery { position: absolute; inset: 0; }
.square-preview .preview-image { height: 100%; min-height: 0; padding: .5rem; }
.square-preview img { height: 100%; }
.preview-loading { position: absolute; top: .75rem; left: .75rem; padding: .25rem .5rem; background: var(--scaid-soft, #f8f7f4); color: var(--scaid-muted, #636963); font-size: .8125rem; pointer-events: none; }
.preview-toolbar { display: flex; flex-wrap: wrap; align-items: center; gap: 4px 8px; padding: 8px 12px; border-top: 1px solid var(--scaid-line, #e3e5df); background: var(--scaid-soft, #f8f7f4); }
.preview-toolbar button, .preview-toolbar a, .preview-recovery button { display: inline-flex; align-items: center; justify-content: center; min-height: 44px; padding: 8px 10px; border: 1px solid var(--scaid-line-strong, #c3c8bf); border-radius: 8px; background: #fff; color: var(--scaid-accent, #88513c); font: inherit; font-size: .8125rem; line-height: 1.3; text-align: center; cursor: pointer; }
.preview-toolbar a { border-color: transparent; background: transparent; text-decoration: underline; text-underline-offset: .18em; }
.preview-toolbar button:hover, .preview-recovery button:hover { border-color: var(--scaid-accent, #88513c); background: var(--scaid-accent-soft, #f2e8e1); }
.figure-preview a:focus-visible, .figure-preview button:focus-visible { outline: 2px solid var(--scaid-focus, #88513c); outline-offset: 3px; }
.figure-preview .preview-image:focus-visible { outline-offset: -4px; }
.preview-recovery { display: flex; min-height: 180px; height: 100%; flex-direction: column; justify-content: center; align-items: center; gap: .75rem; padding: 1.25rem; text-align: center; }
.preview-recovery p { margin: 0; color: var(--scaid-muted); }
</style>

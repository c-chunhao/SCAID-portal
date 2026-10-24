<template>
  <Teleport to="body">
    <dialog ref="dialog" class="figure-viewer" :aria-labelledby="titleId" :aria-describedby="descriptionId" @cancel.prevent="close" @close="afterClose" @keydown="onDialogKeydown">
      <section v-if="isOpen" ref="shell" class="viewer-shell">
        <header class="viewer-heading">
          <div>
            <h2 :id="titleId">Figure viewer</h2>
            <p :id="descriptionId">{{ alt }}</p>
          </div>
          <button ref="closeButton" type="button" class="viewer-button viewer-close" @click="close">Close</button>
        </header>
        <div class="viewer-toolbar" role="group" aria-label="Figure controls">
          <div class="zoom-controls">
            <button type="button" class="viewer-button" :disabled="!loaded || zoom <= MIN_FIGURE_ZOOM" @click="changeZoom(zoom - .25)">Zoom out</button>
            <output class="zoom-value" aria-live="polite" :aria-label="`${zoomPercent}% of fitted size`">{{ zoomPercent }}%</output>
            <button type="button" class="viewer-button" :disabled="!loaded || zoom >= MAX_FIGURE_ZOOM" @click="changeZoom(zoom + .25)">Zoom in</button>
            <button type="button" class="viewer-button" :disabled="!loaded" @click="reset">Reset</button>
          </div>
          <div class="file-controls">
            <button type="button" class="viewer-button" @click="toggleFullscreen">{{ isFullscreen ? 'Exit full screen' : 'Full screen' }}</button>
            <a class="viewer-original" :href="src" target="_blank" rel="noopener" :aria-label="fullLabel || `Open original ${alt} in a new tab`">Open original</a>
          </div>
        </div>
        <p class="viewer-hint">100% fits the whole figure. Zoom in, then drag or use arrow keys to inspect labels.</p>
        <p v-if="fullscreenNotice" class="fullscreen-notice" role="status">{{ fullscreenNotice }}</p>
        <div ref="viewport" class="viewer-viewport" :class="{ 'can-pan': loaded && zoom > 1, 'is-dragging': dragging }" tabindex="0" role="region" aria-label="Figure image viewport" :aria-busy="!loaded && !failed" @pointerdown="startDrag" @pointermove="moveDrag" @pointerup="endDrag" @pointercancel="endDrag" @lostpointercapture="endDrag" @keydown="onViewportKeydown">
          <div class="viewer-stage" :style="stageStyle">
            <img :key="imageSrc" class="viewer-image" :class="{ 'is-loading': !loaded }" :src="imageSrc" :alt="alt" :style="imageStyle" decoding="async" :draggable="false" @load="onLoad" @error="onError" />
          </div>
          <div v-if="!loaded && !failed" class="viewer-status" role="status">Loading original figure…</div>
          <div v-if="failed" class="viewer-status viewer-recovery" role="status">
            <p>The original figure could not be loaded.</p>
            <div>
              <button type="button" class="viewer-button" @click="retry">Retry figure</button>
              <a class="viewer-original" :href="src" target="_blank" rel="noopener">Open original in a new tab</a>
            </div>
          </div>
        </div>
      </section>
    </dialog>
  </Teleport>
</template>

<script setup>
import { computed, nextTick, onBeforeUnmount, onMounted, ref, useId, watch } from 'vue'
import { clampFigureZoom, figureZoomScroll, fitFigure, MAX_FIGURE_ZOOM, MIN_FIGURE_ZOOM } from '../utils/figureViewport.js'

const props = defineProps({
  src: { type: String, required: true },
  alt: { type: String, required: true },
  fullLabel: { type: String, default: '' }
})
const emit = defineEmits(['imageerror'])
const titleId = useId()
const descriptionId = useId()
const dialog = ref(null)
const shell = ref(null)
const viewport = ref(null)
const closeButton = ref(null)
const isOpen = ref(false)
const loaded = ref(false)
const failed = ref(false)
const attempt = ref(0)
const zoom = ref(1)
const naturalSize = ref({ width: 0, height: 0 })
const viewportSize = ref({ width: 0, height: 0 })
const isFullscreen = ref(false)
const fullscreenNotice = ref('')
const dragging = ref(false)
let drag = null
let returnFocus = null
let resizeObserver = null
let previousOverflow = null

const imageSrc = computed(() => {
  if (!attempt.value) return props.src
  const url = new URL(props.src, window.location.origin)
  url.searchParams.set('figure_retry', String(attempt.value))
  return url.href
})
const fittedSize = computed(() => fitFigure(naturalSize.value.width, naturalSize.value.height, viewportSize.value.width, viewportSize.value.height))
const displayedSize = computed(() => ({ width: fittedSize.value.width * zoom.value, height: fittedSize.value.height * zoom.value }))
const stageStyle = computed(() => ({ width: `${displayedSize.value.width}px`, height: `${displayedSize.value.height}px` }))
const imageStyle = computed(() => ({ width: `${displayedSize.value.width}px`, height: `${displayedSize.value.height}px` }))
const zoomPercent = computed(() => Math.round(zoom.value * 100))

function capturePosition() {
  return {
    scrollLeft: viewport.value?.scrollLeft || 0,
    scrollTop: viewport.value?.scrollTop || 0,
    imageWidth: displayedSize.value.width,
    imageHeight: displayedSize.value.height,
    viewportWidth: viewportSize.value.width,
    viewportHeight: viewportSize.value.height
  }
}

async function restorePosition(previous) {
  const position = figureZoomScroll({
    ...previous,
    nextImageWidth: displayedSize.value.width,
    nextImageHeight: displayedSize.value.height,
    nextViewportWidth: viewportSize.value.width,
    nextViewportHeight: viewportSize.value.height
  })
  await nextTick()
  if (viewport.value) {
    viewport.value.scrollLeft = position.left
    viewport.value.scrollTop = position.top
  }
}

function measure() {
  if (!viewport.value) return
  const width = viewport.value.clientWidth
  const height = viewport.value.clientHeight
  if (viewportSize.value.width === width && viewportSize.value.height === height) return
  const previous = capturePosition()
  viewportSize.value = { width, height }
  restorePosition(previous)
}

async function changeZoom(nextZoom) {
  if (!loaded.value) return
  endDrag()
  const previous = capturePosition()
  zoom.value = clampFigureZoom(nextZoom)
  await restorePosition(previous)
}

async function reset() {
  await changeZoom(1)
  if (viewport.value) {
    viewport.value.scrollLeft = 0
    viewport.value.scrollTop = 0
  }
}

async function open(options = {}) {
  if (isOpen.value) return
  returnFocus = options.returnFocus || document.activeElement
  attempt.value = 0
  loaded.value = false
  failed.value = false
  zoom.value = 1
  naturalSize.value = { width: 0, height: 0 }
  fullscreenNotice.value = ''
  isOpen.value = true
  await nextTick()
  dialog.value.showModal()
  previousOverflow = document.documentElement.style.overflow
  document.documentElement.style.overflow = 'hidden'
  closeButton.value?.focus({ preventScroll: true })
  measure()
  resizeObserver = new ResizeObserver(measure)
  resizeObserver.observe(viewport.value)
  if (options.fullscreen) await toggleFullscreen()
}

function onDialogKeydown(event) {
  if (event.key !== 'Tab' || event.ctrlKey || event.metaKey || event.altKey || !dialog.value?.open) return
  const focusable = [...dialog.value.querySelectorAll('button:not(:disabled), a[href], [tabindex]:not([tabindex="-1"])')]
    .filter(element => element.getClientRects().length && getComputedStyle(element).visibility !== 'hidden')
  if (!focusable.length) return
  const first = focusable[0]
  const last = focusable[focusable.length - 1]
  const active = document.activeElement
  if (event.shiftKey && (active === first || !dialog.value.contains(active))) {
    event.preventDefault()
    last.focus({ preventScroll: true })
  } else if (!event.shiftKey && (active === last || !dialog.value.contains(active))) {
    event.preventDefault()
    first.focus({ preventScroll: true })
  }
}

function exitOwnedFullscreen() {
  if (shell.value && document.fullscreenElement === shell.value && document.exitFullscreen) {
    return document.exitFullscreen().catch(() => {})
  }
}

function close() {
  exitOwnedFullscreen()
  dialog.value?.close()
}

function afterClose() {
  endDrag()
  resizeObserver?.disconnect()
  resizeObserver = null
  if (previousOverflow !== null) {
    document.documentElement.style.overflow = previousOverflow
    previousOverflow = null
  }
  isOpen.value = false
  isFullscreen.value = false
  if (returnFocus?.isConnected) returnFocus.focus({ preventScroll: true })
  returnFocus = null
}

async function toggleFullscreen() {
  fullscreenNotice.value = ''
  if (document.fullscreenElement === shell.value) {
    await exitOwnedFullscreen()
    return
  }
  const target = shell.value
  if (!target?.requestFullscreen || document.fullscreenEnabled === false) {
    fullscreenNotice.value = 'Full screen is unavailable in this browser. You can keep using this viewer.'
    return
  }
  try {
    await target.requestFullscreen()
    // A request may finish after the user has already closed the viewer.
    if (!isOpen.value && document.fullscreenElement === target) await document.exitFullscreen()
  } catch {
    if (isOpen.value) fullscreenNotice.value = 'Full screen could not be opened. You can keep using this viewer.'
  }
}

function onFullscreenChange() {
  isFullscreen.value = !!shell.value && document.fullscreenElement === shell.value
  if (isOpen.value) nextTick(measure)
}

function onLoad(event) {
  if (!isOpen.value) return
  naturalSize.value = { width: event.target.naturalWidth, height: event.target.naturalHeight }
  failed.value = false
  loaded.value = true
  zoom.value = 1
  measure()
  nextTick(reset)
}

function onError(event) {
  if (!isOpen.value) return
  endDrag()
  loaded.value = false
  failed.value = true
  emit('imageerror', { source: 'original', src: props.src, event })
}

function retry() {
  attempt.value += 1
  failed.value = false
  loaded.value = false
  zoom.value = 1
  naturalSize.value = { width: 0, height: 0 }
}

function startDrag(event) {
  if (!loaded.value || zoom.value <= 1 || !event.isPrimary || event.button !== 0) return
  event.preventDefault()
  viewport.value.focus({ preventScroll: true })
  drag = { pointerId: event.pointerId, x: event.clientX, y: event.clientY, left: viewport.value.scrollLeft, top: viewport.value.scrollTop }
  dragging.value = true
  viewport.value.setPointerCapture(event.pointerId)
}

function moveDrag(event) {
  if (!drag || event.pointerId !== drag.pointerId) return
  event.preventDefault()
  viewport.value.scrollLeft = drag.left + drag.x - event.clientX
  viewport.value.scrollTop = drag.top + drag.y - event.clientY
}

function endDrag() {
  const pointerId = drag?.pointerId
  drag = null
  dragging.value = false
  if (pointerId !== undefined && viewport.value?.hasPointerCapture(pointerId)) viewport.value.releasePointerCapture(pointerId)
}

function onViewportKeydown(event) {
  if (event.ctrlKey || event.metaKey || event.altKey || !loaded.value) return
  const direction = { ArrowLeft: [-64, 0], ArrowRight: [64, 0], ArrowUp: [0, -64], ArrowDown: [0, 64] }[event.key]
  if (direction && zoom.value > 1) {
    event.preventDefault()
    viewport.value.scrollLeft += direction[0]
    viewport.value.scrollTop += direction[1]
  } else if (event.key === '+' || event.key === '=') {
    event.preventDefault()
    changeZoom(zoom.value + .25)
  } else if (event.key === '-') {
    event.preventDefault()
    changeZoom(zoom.value - .25)
  } else if (event.key === '0') {
    event.preventDefault()
    reset()
  }
}

watch(() => props.src, () => { if (isOpen.value) close() })
onMounted(() => document.addEventListener('fullscreenchange', onFullscreenChange))
onBeforeUnmount(() => {
  document.removeEventListener('fullscreenchange', onFullscreenChange)
  exitOwnedFullscreen()
  if (dialog.value?.open) dialog.value.close()
  afterClose()
})
defineExpose({ open, close })
</script>

<style scoped>
.figure-viewer { width: min(1240px, calc(100vw - 32px)); height: min(900px, calc(100dvh - 32px)); max-width: none; max-height: none; margin: auto; padding: 0; overflow: hidden; border: 1px solid var(--scaid-line-strong, #c3c8bf); border-radius: 12px; background: var(--scaid-soft, #f8f7f4); color: var(--scaid-ink, #252824); font-family: var(--scaid-font, sans-serif); }
.figure-viewer::backdrop { background: rgb(37 32 28 / 62%); }
.viewer-shell { display: flex; flex-direction: column; width: 100%; height: 100%; min-height: 0; background: var(--scaid-soft, #f8f7f4); }
.viewer-shell:fullscreen { width: 100vw; height: 100dvh; padding: env(safe-area-inset-top) env(safe-area-inset-right) env(safe-area-inset-bottom) env(safe-area-inset-left); }
.viewer-heading { display: flex; flex-shrink: 0; align-items: flex-start; justify-content: space-between; gap: 20px; padding: 16px 20px 12px; }
.viewer-heading > div { min-width: 0; }
.viewer-heading h2 { margin: 0 0 4px; color: var(--scaid-ink, #252824); font-size: 1.25rem; font-weight: 600; line-height: 1.3; }
.viewer-heading p { max-width: 75ch; margin: 0; color: var(--scaid-muted, #636963); font-size: .875rem; line-height: 1.45; overflow-wrap: anywhere; }
.viewer-toolbar { display: flex; flex-shrink: 0; flex-wrap: wrap; justify-content: space-between; gap: 8px 20px; padding: 0 20px 8px; }
.zoom-controls, .file-controls { display: flex; flex-wrap: wrap; align-items: center; gap: 4px; }
.viewer-button, .viewer-original { display: inline-flex; align-items: center; justify-content: center; min-height: 44px; padding: 8px 12px; border: 1px solid var(--scaid-line-strong, #c3c8bf); border-radius: 8px; background: #fff; color: var(--scaid-accent, #88513c); font: inherit; font-size: .875rem; line-height: 1.3; cursor: pointer; }
.viewer-button:hover:not(:disabled) { border-color: var(--scaid-accent, #88513c); background: var(--scaid-accent-soft, #f2e8e1); }
.viewer-button:disabled { color: var(--scaid-muted, #636963); opacity: .65; cursor: default; }
.viewer-original { border-color: transparent; background: transparent; text-decoration: underline; text-underline-offset: .18em; }
.viewer-close { flex-shrink: 0; }
.zoom-value { min-width: 4.25rem; color: var(--scaid-ink, #252824); font-size: .875rem; font-variant-numeric: tabular-nums; text-align: center; }
.viewer-hint { flex-shrink: 0; margin: 0; padding: 0 20px 12px; color: var(--scaid-muted, #636963); font-size: .8125rem; line-height: 1.45; }
.fullscreen-notice { flex-shrink: 0; margin: 0; padding: 0 20px 12px; color: var(--scaid-accent, #88513c); font-size: .8125rem; line-height: 1.45; }
.viewer-viewport { position: relative; flex: 1; min-width: 0; min-height: 0; overflow: auto; overscroll-behavior: contain; border-top: 1px solid var(--scaid-line, #e3e5df); background: #fff; scrollbar-color: var(--scaid-accent, #88513c) var(--scaid-soft, #f8f7f4); }
.viewer-stage { display: grid; place-items: center; min-width: 100%; min-height: 100%; }
.viewer-image { display: block; max-width: none; max-height: none; object-fit: contain; background: #fff; user-select: none; -webkit-user-drag: none; }
.viewer-image.is-loading { visibility: hidden; }
.can-pan { cursor: grab; touch-action: pinch-zoom; }
.is-dragging { cursor: grabbing; user-select: none; }
.viewer-status { position: absolute; inset: 0; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 16px; padding: 24px; background: #fff; color: var(--scaid-muted, #636963); font-size: .9375rem; text-align: center; }
.viewer-recovery p { margin: 0; }
.viewer-recovery > div { display: flex; flex-wrap: wrap; align-items: center; justify-content: center; gap: 8px; }
.viewer-button:focus-visible, .viewer-original:focus-visible, .viewer-viewport:focus-visible { outline: 2px solid var(--scaid-focus, #88513c); outline-offset: 3px; }
.viewer-viewport:focus-visible { outline-offset: -4px; }
@media (max-width: 640px) {
  .figure-viewer { width: calc(100vw - 16px); height: calc(100dvh - 16px); }
  .viewer-heading { gap: 12px; padding: 12px 12px 8px; }
  .viewer-heading h2 { font-size: 1.125rem; }
  .viewer-heading p { font-size: .8125rem; }
  .viewer-toolbar { gap: 4px 12px; padding: 0 12px 8px; }
  .viewer-button, .viewer-original { padding: 8px 10px; font-size: .8125rem; }
  .viewer-hint, .fullscreen-notice { padding: 0 12px 8px; }
  .zoom-value { min-width: 3.5rem; }
}
</style>

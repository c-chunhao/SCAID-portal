<template>
  <div class="share-query">
    <button type="button" class="share-query-button" :disabled="disabled" @click="copyLink">Copy query link</button>
    <span v-if="copied" role="status">Link copied.</span>
    <div v-if="manualCopy" class="manual-copy">
      <label :for="inputId">Copy this query URL</label>
      <input :id="inputId" ref="linkInput" :value="absoluteUrl" readonly @focus="$event.target.select()" />
    </div>
  </div>
</template>

<script setup>
import { computed, nextTick, ref, useId, watch } from 'vue'
const props = defineProps({ href: { type: String, required: true }, disabled: Boolean })
const inputId = `query-link-${useId()}`
const linkInput = ref()
const copied = ref(false)
const manualCopy = ref(false)
const absoluteUrl = computed(() => new URL(props.href, window.location.origin).href)
watch(() => props.href, () => { copied.value = false; manualCopy.value = false })
const copyLink = async () => {
  copied.value = false
  try {
    await navigator.clipboard.writeText(absoluteUrl.value)
    copied.value = true
  } catch {
    manualCopy.value = true
    await nextTick()
    linkInput.value?.focus()
    linkInput.value?.select()
  }
}
</script>

<style scoped>
.share-query { display: flex; align-items: center; flex-wrap: wrap; gap: .6rem; font-size: .85rem; }
.share-query-button { min-height: 44px; padding: .5rem .8rem; border: 1px solid var(--scaid-line); border-radius: 6px; background: #fff; color: var(--scaid-accent); font: inherit; cursor: pointer; }
.share-query-button:hover { background: var(--scaid-soft); }
.share-query-button:focus-visible, .manual-copy input:focus-visible { outline: 2px solid var(--scaid-focus); outline-offset: 3px; }
.share-query-button:disabled { opacity: .6; cursor: default; }
.manual-copy { flex-basis: 100%; min-width: 0; }
.manual-copy label { display: block; margin-bottom: .3rem; }
.manual-copy input { width: 100%; min-height: 44px; padding: .5rem; border: 1px solid var(--scaid-line); border-radius: 4px; color: var(--scaid-ink); background: #fff; font: inherit; }
</style>

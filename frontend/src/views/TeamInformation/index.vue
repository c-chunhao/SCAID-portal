<template>
  <div class="container_box">
    <div class="body2 container">
      <header class="atlas-intro">
        <h1>Integrated lymphocyte atlases</h1>
        <p>
          Two integrated atlases give the collection a common annotation vocabulary: a T/NK atlas
          (937,825 cells; 4 major, 27 intermediate and 43 fine-level labels) and a B/plasma atlas
          (253,635 cells; 3, 6 and 20 labels). The figures below are the published atlas views and
          marker panels; cell positions and labels are not recomputed on this site. A well-mixed
          embedding does not by itself establish that expression values are free of study or platform
          effects, and these atlases do not replace donor-level comparisons with appropriate controls.
        </p>
      </header>
      <div class="image-grid">
        <div 
          v-for="image in pngImages" 
          :key="image.name" 
          class="image-item"
        >
          <figure class="image-container">
            <FigurePreview :src="image.url" :full-src="image.original" :alt="caption(image.name)" />
            <figcaption>{{ caption(image.name) }}</figcaption>
          </figure>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import FigurePreview from '@/components/FigurePreview.vue'
defineOptions({ name: 'TeamInformation' })

// Vite resolves both previews and originals into real build assets. If a
// preview is absent, use the packaged original rather than a /src URL.
const previews = import.meta.glob('/src/assets/images/dataIntegration_preview/*.jpg', { eager: true, query: '?url', import: 'default' })
const originals = import.meta.glob('/src/assets/images/dataIntegration/*.png', { eager: true, query: '?url', import: 'default' })
const previewByStem = Object.fromEntries(Object.entries(previews).map(([path, url]) => [path.split('/').pop().replace(/\.jpg$/i, ''), url]))
const pngImages = Object.entries(originals).map(([path, original]) => {
  const name = path.split('/').pop()
  return { name, original, url: previewByStem[name.replace(/\.png$/i, '')] || original }
}).sort((a, b) => a.name.localeCompare(b.name, undefined, { numeric: true }))

const caption = (fileName) => {
  const stem = String(fileName).replace(/\.(png|jpe?g)$/i, '')
  const match = stem.match(/^FIG(\d+)_(.*)$/i)
  const body = (match ? match[2] : stem)
    .replace(/^T_NK/, 'T/NK').replace(/^B_Plasma/, 'B/plasma')
    .replace(/_celltype$/, ' cell types').replace(/_Marker_layer(\d)/, ' marker panel, layer $1')
    .replace(/_/g, ' ')
  return match ? `Fig. ${match[1]} · ${body}` : body
}
</script>

<style scoped lang="scss">
.atlas-intro { max-width: 75ch; margin: 0 0 32px; h1 { font-size: 2.25rem; font-weight: 600; line-height: 1.2; letter-spacing: -.025em; color: var(--scaid-ink); margin: 0 0 16px; } p { color: var(--scaid-muted); line-height: 1.7; margin: 0; } }
figcaption { font-size: .9rem; color: var(--scaid-muted); padding: 16px 20px; background: #fff; border-top: 1px solid var(--scaid-line); }
.body2 {
  padding-top: 40px;
  padding-bottom: 48px;
  width: 100%;
  
  .title {
    text-align: center;
    margin-bottom: 2rem;
    color: #333;
    font-size: 1.8rem;
  }
  
  .image-grid {
    display: flex;
      flex-direction: column;
    gap: 32px;
    padding: 0;
  }
  
  .image-item {
    border: 1px solid var(--scaid-line);
    border-radius: var(--scaid-radius);
    overflow: hidden;
    background: white;
    width: 100%;
    
    .image-container {
      display: flex;
      flex-direction: column;
      align-items: stretch;
      justify-content: center;
      background: #fff;
      margin: 0;
      
      .png-image {
        width: 100%;
        height: auto;
        max-width: 100%;
        display: block;
      }
    }
    
    .image-info {
      padding: 1rem;
      
      .image-name {
        display: block;
        font-weight: 600;
        color: #333;
        margin-bottom: 0.5rem;
        font-size: 0.9rem;
      }
      
      .image-size {
        color: #666;
        font-size: 0.8rem;
      }
    }
  }
}
@media (max-width: 767px) {
  .body2 { padding-top: 28px; }
  .atlas-intro h1 { font-size: 1.8rem; }
  figcaption { padding: 14px 16px; }
}
</style>

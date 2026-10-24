<template>
  <figure class="home-workflow" aria-label="SCAID exploration workflow">
    <svg class="workflow-connectors" viewBox="0 0 420 220" preserveAspectRatio="none" aria-hidden="true">
      <defs>
        <marker :id="arrowId" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
          <path d="M1 1 6 4 1 7" />
        </marker>
      </defs>
      <path d="M184 26H234" :marker-end="`url(#${arrowId})`" />
      <path d="M315 50V65H210V76" />
      <path d="M70 88V76H350V88" />
      <path d="M210 76V88" />
      <path d="M70 148V162H350V148" />
      <path d="M210 148V177" :marker-end="`url(#${arrowId})`" />
    </svg>
    <ol class="workflow-stages">
      <li class="workflow-start">
        <span class="stage-symbol" aria-hidden="true">
          <svg viewBox="0 0 32 32"><ellipse cx="16" cy="7" rx="10" ry="4" /><path d="M6 7v9c0 2.2 4.5 4 10 4s10-1.8 10-4V7M6 16v9c0 2.2 4.5 4 10 4s10-1.8 10-4v-9" /><path d="M11 13h2m-2 9h2" /></svg>
        </span>
        <span>Find data</span>
      </li>
      <li class="workflow-cells">
        <span class="stage-symbol" aria-hidden="true">
          <svg viewBox="0 0 32 32"><circle cx="10" cy="10" r="7" /><circle cx="23" cy="12" r="5" /><circle cx="17" cy="24" r="6" /><circle cx="10" cy="10" r="2" /><circle cx="23" cy="12" r="1.5" /><circle cx="17" cy="24" r="2" /></svg>
        </span>
        <span>Cell types<span class="sr-only">: view annotations</span></span>
      </li>
      <li class="workflow-analyses">
        <span class="sr-only">Explore available precomputed analyses:</span>
        <ul class="analysis-branches">
          <li>
            <span class="stage-symbol" aria-hidden="true">
              <svg viewBox="0 0 32 32"><path d="M8 3c0 10 16 16 16 26M24 3c0 10-16 16-16 26M9 7h14M12 12h8M12 20h8M9 25h14" /></svg>
            </span>
            <span>Genes<span class="sr-only">: gene expression maps</span></span>
          </li>
          <li>
            <span class="stage-symbol" aria-hidden="true">
              <svg viewBox="0 0 32 32"><path d="m16 8-9 9m9-9 9 9m-18 3 9 7 9-7" /><circle cx="16" cy="5" r="3" /><circle cx="6" cy="19" r="3" /><circle cx="26" cy="19" r="3" /><circle cx="16" cy="27" r="3" /></svg>
            </span>
            <span>Pathways<span class="sr-only">: pathway scores</span></span>
          </li>
          <li>
            <span class="stage-symbol" aria-hidden="true">
              <svg viewBox="0 0 32 32"><circle cx="8" cy="16" r="5" /><circle cx="25" cy="7" r="4" /><circle cx="25" cy="25" r="4" /><path d="m13 13 8-4m-8 10 8 4M8 11V5m0 16v6M21 16h7" /></svg>
            </span>
            <span>CellChat<span class="sr-only">: predicted cell communication</span></span>
          </li>
        </ul>
      </li>
      <li class="workflow-download">
        <span class="stage-symbol" aria-hidden="true">
          <svg viewBox="0 0 32 32"><path d="M16 3v17m-6-6 6 6 6-6M5 21v7h22v-7" /></svg>
        </span>
        <span>Download data<span class="sr-only">: available matrices and scores</span></span>
      </li>
    </ol>
    <figcaption class="sr-only">Find a dataset, view its cell annotations, explore published analyses and download the available data.</figcaption>
  </figure>
</template>

<script setup>
import { useId } from 'vue'
const arrowId = `workflow-arrow-${useId()}`
</script>

<style scoped>
.home-workflow { position: relative; width: 100%; height: 220px; margin: 0 0 24px; color: var(--scaid-ink); }
.workflow-connectors { position: absolute; inset: 0; width: 100%; height: 100%; overflow: visible; pointer-events: none; }
.workflow-connectors path { fill: none; stroke: var(--scaid-line-strong); stroke-width: 1.4; stroke-linecap: round; stroke-linejoin: round; vector-effect: non-scaling-stroke; }
.workflow-connectors marker path { stroke: var(--scaid-accent); stroke-width: 1.5; }
.workflow-stages, .analysis-branches { margin: 0; padding: 0; list-style: none; }
.workflow-stages { display: grid; grid-template-columns: 1fr 1fr; grid-template-rows: 52px 36px 60px 29px 43px; height: 100%; }
.workflow-start, .workflow-cells, .workflow-download { display: flex; align-items: center; justify-content: center; gap: 10px; background: var(--scaid-soft); font-size: .9375rem; font-weight: 500; line-height: 1.4; }
.workflow-start { grid-column: 1; grid-row: 1; margin-right: 22px; }
.workflow-cells { grid-column: 2; grid-row: 1; margin-left: 22px; }
.stage-symbol { display: inline-flex; flex: 0 0 auto; align-items: center; justify-content: center; width: 36px; height: 36px; color: var(--scaid-accent); }
.stage-symbol svg { width: 30px; height: 30px; fill: var(--scaid-accent-soft); stroke: currentColor; stroke-width: 1.5; stroke-linecap: round; stroke-linejoin: round; }
.workflow-start .stage-symbol svg { fill: var(--scaid-accent-soft); }
.workflow-analyses { grid-column: 1 / -1; grid-row: 3; }
.analysis-branches { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); height: 100%; }
.analysis-branches > li { display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 4px; color: var(--scaid-accent); font-size: .875rem; font-weight: 500; line-height: 1.3; }
.analysis-branches .stage-symbol { width: 34px; height: 34px; background: var(--scaid-soft); }
.analysis-branches > li > span:last-child { padding: 0 4px; background: var(--scaid-soft); }
.workflow-download { grid-column: 1 / -1; grid-row: 5; justify-self: center; padding: 0 12px; }
.workflow-download .stage-symbol { width: 30px; height: 30px; }
.workflow-download svg { width: 26px; height: 26px; fill: none; }
@media (max-width: 991px) {
  .home-workflow { max-width: 440px; margin-right: auto; margin-left: auto; }
}
@media (max-width: 600px) {
  .home-workflow { height: 196px; margin-bottom: 20px; }
  .workflow-stages { grid-template-rows: 46px 32px 54px 26px 38px; }
  .workflow-start, .workflow-cells, .workflow-download { gap: 6px; font-size: .875rem; }
  .workflow-start { margin-right: 14px; }
  .workflow-cells { margin-left: 14px; }
  .stage-symbol { width: 30px; height: 30px; }
  .stage-symbol svg { width: 26px; height: 26px; }
  .analysis-branches .stage-symbol { width: 30px; height: 30px; }
  .analysis-branches { font-size: .8125rem; }
}
</style>

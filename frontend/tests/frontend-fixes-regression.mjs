// Focused user-visible regressions for routing, shareable queries and preview recovery.
import assert from 'node:assert/strict'
import fs from 'node:fs/promises'
import process from 'node:process'
const { chromium } = await import(process.env.PLAYWRIGHT_CORE_PATH || '/tmp/scaid-audit-tools/node_modules/playwright-core/index.mjs')
const base = process.env.SCAID_UI_URL || 'http://127.0.0.1:15190'
const output = process.env.SCAID_CHECK_OUTPUT || '/home/data/KWQ/全项目修复-2026-09-30/frontend'
const browser = await chromium.launch({ executablePath: process.env.CHROMIUM_PATH || '/home/bio/.cache/ms-playwright/chromium-1234/chrome-linux64/chrome', headless: true, args: ['--no-sandbox'] })
const context = await browser.newContext({ viewport: { width: 1280, height: 900 }, permissions: ['clipboard-read', 'clipboard-write'] })
const page = await context.newPage()
const errors = [], requests = [], checks = []
page.on('pageerror', error => errors.push(error.message))
let failGenePreview = true, failAtlasPreview = true, delayOldGene = false
const svg = '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"><rect width="200" height="100" fill="steelblue"/></svg>'
const rows = [
  { id: 130, abbreviation: 'CI', dataset_id: 'GSE121380', tissue: 'colon', sample_size: 5, disease_label: 'Colitis, inflamed' },
  { id: 131, abbreviation: 'PsA', dataset_id: 'E-MTAB-8207', tissue: 'PBMC', sample_size: 8, disease_label: 'Psoriatic arthritis' }
]
const pause = ms => new Promise(resolve => setTimeout(resolve, ms))
async function until(check, message) {
  for (let attempt = 0; attempt < 160; attempt += 1) {
    if (await check()) return
    await pause(50)
  }
  throw new Error(message)
}
const navigate = query => page.evaluate(path => document.querySelector('#app').__vue_app__.config.globalProperties.$router.push(path), query)
const row = name => page.locator('.el-table__body-wrapper tr').filter({ hasText: name })
const genePlots = () => page.locator('.file-item img')
const datasetPlots = () => page.locator('.dataset-details .file_item img')
await page.route(/\/(?:src\/assets\/images\/dataIntegration_preview\/|assets\/FIG\d+).*\.jpg(?:\?|$)/, async route => {
  if (route.request().resourceType() !== 'image') return route.continue()
  if (failAtlasPreview) { failAtlasPreview = false; return route.fulfill({ status: 503, body: 'temporary fixture preview error' }) }
  return route.continue()
})
await page.route('**/api/**', async route => {
  const url = new URL(route.request().url()), params = Object.fromEntries(url.searchParams), path = url.pathname
  if (!path.startsWith('/api/')) return route.fallback()
  requests.push({ path, ...params })
  if (/^\/api\/pdf-images\/\d+\/original\/$/.test(path)) return route.fulfill({ contentType: 'image/svg+xml', body: svg })
  if (path.startsWith('/api/fixture-thumbnail/')) {
    if (failGenePreview && params.family === 'gene' && params.name === 'CD3D') { failGenePreview = false; return route.fulfill({ status: 503, body: 'temporary fixture thumbnail error' }) }
    return route.fulfill({ contentType: 'image/svg+xml', body: svg })
  }
  let body
  if (path === '/api/cell-data/') body = { code: 200, data: rows.map(item => ({ abbreviation: item.abbreviation, disease_label: item.disease_label, children: [{ dataset_id: item.dataset_id, children: [{ tissue: item.tissue }] }] })) }
  else if (path === '/api/cell-data-all/') body = rows.filter(item => (!params.abbreviation || item.abbreviation === params.abbreviation) && (!params.dataset_id || item.dataset_id === params.dataset_id) && (!params.tissue || item.tissue === params.tissue) && (!params.name || Object.values(item).join(' ').toLowerCase().includes(params.name.toLowerCase())))
  else if (path === '/api/gene-kegg-enum/') body = { count: 1, results: { code: 200, data: [{ id: 1, name: params.name }] } }
  else if (path === '/api/h5-file/') body = [{ id: 1, name: 'fixture.h5', method: 'RNA', download_url: '/api/h5-file/1/download/' }]
  else if (path === '/api/pdf-images/') {
    const family = params.GeneUmap ? 'gene' : params.KEGGScore ? 'kegg' : 'dataset'
    if (delayOldGene && family === 'gene' && params.pdf_name === 'CD4') await pause(600)
    const count = family === 'dataset' ? 31 : family === 'kegg' && !params.kegg_param ? 144 : 48
    const first = (Number(params.page || 1) - 1) * Number(params.page_size || 10)
    const data = Array.from({ length: Math.min(Number(params.page_size || 10), Math.max(0, count - first)) }, (_, index) => {
      const id = first + index + 1
      return { id, name: `${params.pdf_name}.pdf`, image: `/api/pdf-images/${id}/original/`, thumb_url: `/api/fixture-thumbnail/${id}/?family=${family}&name=${encodeURIComponent(params.pdf_name)}`, relative_path: `CI/GSE121380/colon/${family === 'kegg' ? `KEGGScore/${params.kegg_param || 'AUCell'}` : 'GeneUmap'}` }
    })
    body = { count, results: { code: 200, data } }
  } else return route.fulfill({ status: 404, json: { detail: 'Unknown fixture endpoint' } })
  try { return await route.fulfill({ json: body }) } catch (error) { if (!route.request().failure()) throw error }
})
try {
  await page.goto(`${base}/not-a-real-page`)
  await page.getByRole('heading', { name: 'Page not found' }).waitFor()
  assert.ok(await page.getByRole('navigation', { name: 'Main navigation' }).isVisible())
  await page.getByRole('link', { name: 'Browse datasets', exact: true }).click()
  await row('GSE121380').waitFor()
  checks.push('Unknown routes show an accessible recovery page inside the normal navigation.')

  await page.goto(`${base}/Dataset?id=CI&dataset=GSE121380&tissue=colon&record=130&analysis=LRcircle&page=2`)
  await until(async () => await datasetPlots().count() === 7, 'shared dataset selection/page did not restore')
  assert.match(await page.locator('.selection-title').innerText(), /GSE121380/)
  assert.match(await page.locator('.analysis-help').innerText(), /Ligand–receptor networks/)
  assert.equal(await page.locator('.figure-title').first().innerText(), 'LRcircle')
  assert.equal(await page.locator('.file_item > a').first().getAttribute('href'), '/api/pdf-images/25/original/')
  assert.ok(await page.getByText('Sample labels', { exact: true }).isVisible())
  await page.getByRole('button', { name: 'Copy query link', exact: true }).click()
  const copiedDataset = await page.evaluate(() => navigator.clipboard.readText())
  assert.equal(new URL(copiedDataset).searchParams.get('record'), '130')
  assert.equal(new URL(copiedDataset).searchParams.get('analysis'), 'LRcircle')
  assert.equal(new URL(copiedDataset).searchParams.get('page'), '2')
  await page.getByRole('menuitem', { name: 'Marker', exact: true }).click()
  await until(async () => new URL(page.url()).searchParams.get('analysis') === 'Main_Marker' && await datasetPlots().count() === 24, 'analysis selection did not update URL')
  await page.reload()
  await until(async () => await datasetPlots().count() === 24 && (await page.locator('.dataset-menu .el-menu-item.is-active').innerText()).trim() === 'Marker', 'reload did not restore Marker')
  checks.push('Dataset links preserve filters, selected record, analysis, figure page and same-origin original image URLs; copy and reload restore them.')

  await page.goto(`${base}/Gene?gene=CD3D&disease=CI&dataset=GSE121380&page=2&size=20`)
  await until(async () => await genePlots().count() >= 7 && await page.getByRole('button', { name: 'Retry figure', exact: true }).count() >= 1, 'failed Gene preview had no recovery')
  assert.ok(requests.some(item => item.path === '/api/pdf-images/' && item.GeneUmap && item.pdf_name === 'CD3D' && item.page === '2' && item.page_size === '20' && item.abbreviation === 'CI' && item.dataset_id === 'GSE121380'))
  await page.getByRole('button', { name: 'Retry figure', exact: true }).first().click()
  await until(async () => await page.getByRole('button', { name: 'Retry figure', exact: true }).count() === 0, 'Gene preview retry did not recover')
  await page.getByRole('button', { name: 'Copy query link', exact: true }).click()
  assert.equal(new URL(await page.evaluate(() => navigator.clipboard.readText())).searchParams.get('size'), '20')
  await navigate('/Gene?gene=TRIM21')
  await until(async () => await genePlots().count() === 10 && (await genePlots().first().getAttribute('alt')).includes('TRIM21'), 'same-component Gene query retained the old plots')
  assert.match(await page.locator('.query-field .el-select').innerText(), /TRIM21/)
  delayOldGene = true
  await navigate('/Gene?gene=CD4')
  await until(() => requests.some(item => item.GeneUmap && item.pdf_name === 'CD4'), 'delayed old query was not requested')
  await navigate('/Gene?gene=MS4A1')
  await until(async () => await genePlots().count() === 10 && (await genePlots().first().getAttribute('alt')).includes('MS4A1'), 'new query missing')
  await pause(650)
  assert.ok((await genePlots().first().getAttribute('alt')).includes('MS4A1'), 'late old query replaced current content')
  await page.goBack()
  await until(async () => (await page.locator('.query-field .el-select').innerText()).includes('CD4') && await genePlots().count() === 10, 'back navigation did not restore gene')
  checks.push('Gene query changes and browser Back restore controls/results; late prior requests cannot overwrite them; preview failures are visible and retryable.')

  await page.goto(`${base}/KEGG?kegg=hsa00010&disease=CI&dataset=GSE121380&method=AUCell&page=2&size=20`)
  await until(async () => await genePlots().count() === 20, 'shared KEGG query did not restore')
  assert.ok(requests.some(item => item.KEGGScore && item.pdf_name === 'hsa00010' && item.kegg_param === 'AUCell' && item.page === '2' && item.page_size === '20'))
  await page.locator('.results-toolbar .btn-next').click()
  await until(async () => new URL(page.url()).searchParams.get('page') === '3' && await genePlots().count() === 8, 'KEGG pagination did not sync URL')
  await navigate('/KEGG?kegg=hsa00240&method=UCell')
  await until(async () => await genePlots().count() === 10 && (await genePlots().first().getAttribute('alt')).includes('hsa00240'), 'same-component KEGG query retained old plots')
  assert.match(await page.locator('.optional-filters').innerText(), /UCell/)
  checks.push('KEGG links and pagination preserve pathway, method, filters and size; changing a reused query refreshes results.')

  await page.goto(`${base}/TeamInformation`)
  const atlasRetry = page.getByRole('button', { name: 'Retry figure', exact: true })
  await atlasRetry.first().waitFor()
  assert.equal(await page.locator('.image-item').count(), 11)
  await atlasRetry.first().click()
  await until(async () => await atlasRetry.count() === 0, 'atlas preview retry did not recover')
  const atlasOriginal = await page.locator('.image-container a').first().getAttribute('href')
  assert.equal((await context.request.get(new URL(atlasOriginal, base).href)).status(), 200, 'atlas original asset is unavailable')
  checks.push('All 11 atlas figures remain available, with visible preview recovery and compiled original assets.')

  await page.setViewportSize({ width: 320, height: 780 })
  await page.goto(`${base}/Dataset?record=130&analysis=LRcontribution`)
  await until(async () => await datasetPlots().count() === 24, 'mobile details did not load')
  assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), 'mobile details overflow the document')
  await page.screenshot({ path: `${output}/mobile-cellchat.png`, fullPage: false })
  await page.setViewportSize({ width: 1280, height: 900 })
  await page.goto(`${base}/Gene?gene=TRIM21`)
  await until(async () => await genePlots().count() === 10, 'final desktop Gene results missing')
  await page.screenshot({ path: `${output}/desktop-gene.png`, fullPage: false })
  assert.deepEqual(errors, [])
  checks.push('Desktop and 320 px mobile checks show no body overflow or uncaught errors.')
  await fs.writeFile(`${output}/browser-fixes-results.json`, JSON.stringify({ passed: true, checks, requestCount: requests.length, errors }, null, 2))
  console.log(`PASS: ${checks.length} focused browser groups.`)
} catch (error) {
  await page.screenshot({ path: `${output}/browser-failure.png`, fullPage: false })
  await fs.writeFile(`${output}/browser-fixes-results.json`, JSON.stringify({ passed: false, checks, requestCount: requests.length, errors, failure: String(error) }, null, 2))
  throw error
} finally { await browser.close() }

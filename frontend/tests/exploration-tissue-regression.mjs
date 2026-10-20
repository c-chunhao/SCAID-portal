// Run against local Vite with the real, read-only SCAID API proxy. No API fixtures.
import assert from 'node:assert/strict'
import fs from 'node:fs/promises'
import process from 'node:process'
import { URL } from 'node:url'
import { setTimeout as pause } from 'node:timers/promises'

const { chromium } = await import(process.env.PLAYWRIGHT_CORE_PATH || 'playwright-core')
const base = process.env.SCAID_UI_URL || 'http://127.0.0.1:15216'
const output = process.env.SCAID_CHECK_OUTPUT || '/tmp/scaid-browser-check'
const browser = await chromium.launch({
  executablePath: process.env.CHROMIUM_PATH || chromium.executablePath(),
  headless: true,
  args: ['--no-sandbox']
})
const page = await browser.newPage({ viewport: { width: 1280, height: 900 }, permissions: ['clipboard-read', 'clipboard-write'] })
const checks = [], requests = [], errors = []
page.on('pageerror', error => errors.push(error.message))
page.on('request', request => {
  const url = new URL(request.url())
  if (url.pathname === '/api/pdf-images/') requests.push(Object.fromEntries(url.searchParams))
})
const filterSelect = label => page.locator('.optional-filters .el-form-item').filter({ has: page.locator('label', { hasText: label }) }).locator('.el-select')
const navigate = path => page.evaluate(target => document.querySelector('#app').__vue_app__.config.globalProperties.$router.push(target), path)
async function openOptionalFilters() {
  const filters = page.locator('details.optional-filters')
  await filters.waitFor()
  if (await filters.getAttribute('open') === null) await filters.locator('summary').click()
}
async function until(check, message) {
  for (let attempt = 0; attempt < 150; attempt += 1) {
    if (await check()) return
    await pause(50)
  }
  throw new Error(message)
}
async function plotResponse(action, expected) {
  const responseTask = page.waitForResponse(response => {
    const url = new URL(response.url())
    return url.pathname === '/api/pdf-images/' && Object.entries(expected).every(([key, value]) => url.searchParams.get(key) === String(value))
  })
  await action()
  const response = await responseTask
  assert.equal(response.status(), 200)
  const body = await response.json()
  const data = body.results?.data || []
  assert.ok(data.length, 'real query should have published plots')
  if (expected.tissue) assert.ok(data.every(item => item.tissue_label === expected.tissue), 'response crossed the selected tissue')
  await until(async () => await page.locator('.file-item').count() === data.length, 'response not reflected in visible plots')
  return body
}

try {
  for (const [path, kind, feature] of [['Gene', 'gene', 'CD3D'], ['KEGG', 'kegg', 'hsa00010']]) {
    const suffix = kind === 'kegg' ? '&method=AUCell' : ''
    await plotResponse(() => page.goto(`${base}/${path}?${kind}=${feature}&disease=PsA&dataset=8207&tissue=SFMC${suffix}`), { pdf_name: feature, tissue: 'SFMC', dataset_id: '8207', ...(kind === 'kegg' ? { kegg_param: 'AUCell' } : {}) })
    await openOptionalFilters()
    assert.match(await filterSelect('Tissue').innerText(), /SFMC/)
    await page.getByRole('button', { name: 'Copy query link', exact: true }).click()
    const shared = new URL(await page.evaluate(() => navigator.clipboard.readText()))
    assert.equal(shared.searchParams.get('tissue'), 'SFMC')
    assert.equal(shared.searchParams.get('dataset'), '8207')
    assert.equal(shared.searchParams.get(kind), feature)
    await plotResponse(() => navigate(`/${path}?${kind}=${feature}&disease=PsA&dataset=8207&tissue=PBMC${suffix}`), { pdf_name: feature, tissue: 'PBMC', dataset_id: '8207' })
    await plotResponse(() => page.goBack(), { pdf_name: feature, tissue: 'SFMC', dataset_id: '8207' })
    await openOptionalFilters()
    assert.match(await filterSelect('Tissue').innerText(), /SFMC/)
    checks.push(`${path}: real PsA/8207 SFMC and PBMC remain separate; copied links and browser Back restore the tissue and results.`)
  }

  const first = await plotResponse(() => page.goto(`${base}/Gene?gene=CD3D&tissue=PBMC`), { pdf_name: 'CD3D', tissue: 'PBMC', page: '1' })
  assert.ok(first.count > 10, 'published PBMC query must exercise pagination')
  await plotResponse(() => page.locator('.results-toolbar .btn-next').click(), { pdf_name: 'CD3D', tissue: 'PBMC', page: '2' })
  assert.equal(new URL(page.url()).searchParams.get('tissue'), 'PBMC')
  assert.equal(new URL(page.url()).searchParams.get('page'), '2')
  await plotResponse(() => page.reload(), { pdf_name: 'CD3D', tissue: 'PBMC', page: '2' })
  checks.push('A real CD3D/PBMC query retains tissue in pagination, refresh and URL restoration.')

  await page.goto(`${base}/Gene?disease=PsA&dataset=8207&tissue=SFMC`)
  await filterSelect('Tissue').waitFor()
  await until(async () => !(await filterSelect('Tissue').locator('.is-loading').count()), 'module filters did not load')
  const before = requests.length
  const queryInput = page.locator('.query-field input').first()
  await queryInput.fill('CD3D')
  await page.getByRole('option', { name: 'CD3D', exact: true }).click()
  assert.equal(requests.length, before, 'selecting a gene should wait for Search plots')
  await plotResponse(() => page.getByRole('button', { name: 'Search plots', exact: true }).click(), { pdf_name: 'CD3D', tissue: 'SFMC', dataset_id: '8207' })
  const afterSubmit = requests.length
  await openOptionalFilters()
  const dataset = filterSelect('Dataset ID')
  await dataset.hover()
  await dataset.locator('.el-select__clear').click()
  assert.match(await filterSelect('Tissue').innerText(), /All tissues/)
  assert.equal(requests.length, afterSubmit, 'changing a parent filter must not submit automatically')
  checks.push('A Dataset-context link restores tissue without a feature; gene selection waits for Search; clearing the dataset clears its tissue without a new result request.')

  assert.deepEqual(errors, [])
  await fs.mkdir(output, { recursive: true })
  await fs.writeFile(`${output}/exploration-tissue-results.json`, JSON.stringify({ passed: true, checks, requestCount: requests.length, requests, errors }, null, 2))
  process.stdout.write(`PASS: ${checks.length} tissue and route groups against the real SCAID API.\n`)
} catch (error) {
  await fs.mkdir(output, { recursive: true })
  await fs.writeFile(`${output}/exploration-tissue-results.json`, JSON.stringify({ passed: false, checks, requests, errors, failure: String(error) }, null, 2))
  throw error
} finally {
  await browser.close()
}

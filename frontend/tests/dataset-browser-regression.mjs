// Run against the Vite server with Node 20+ and playwright-core available.
// SCAID_UI_URL, PLAYWRIGHT_CORE_PATH and CHROMIUM_PATH can override local defaults.
import assert from 'node:assert/strict';
import process from 'node:process';
const { chromium } = await import(process.env.PLAYWRIGHT_CORE_PATH || '/tmp/scaid-audit-tools/node_modules/playwright-core/index.mjs');
const browser = await chromium.launch({ executablePath: process.env.CHROMIUM_PATH || '/home/bio/.cache/ms-playwright/chromium-1234/chrome-linux64/chrome', headless: true, args: ['--no-sandbox'] });
const context = await browser.newContext({ viewport: { width: 1280, height: 900 }, acceptDownloads: true });
const page = await context.newPage();
const pageErrors = [];
page.on('pageerror', error => pageErrors.push(error.message));
const rows = [
  { id: 1, abbreviation: 'JIA', dataset_id: 'GSE111', tissue: 'PBMC/SFMC', sample_size: 10, file_name: 'first.h5' },
  { id: 2, abbreviation: 'RA', dataset_id: 'GSE222', tissue: 'PBMC', sample_size: 20, file_name: 'second.h5' },
];
const requests = [];
let tableFailure = false;
let imageFailure = false;
let delayOldDataset = false;
let completedOldImages = 0;
let completedOldDownloads = 0;
const wait = milliseconds => new Promise(resolve => setTimeout(resolve, milliseconds));
async function until(check, message) {
  for (let attempt = 0; attempt < 150; attempt += 1) {
    if (await check()) return;
    await wait(50);
  }
  throw new Error(message);
}
await page.route('**/fixture-images/**', route => route.fulfill({ contentType: 'image/svg+xml', body: '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"><rect width="200" height="100" fill="steelblue"/></svg>' }));
await page.route('**/api/**', async route => {
  const url = new URL(route.request().url());
  const params = Object.fromEntries(url.searchParams);
  const path = url.pathname;
  if (!path.startsWith('/api/')) return route.fallback();
  requests.push({ path, ...params });
  let body;
  if (path === '/api/cell-data/') {
    body = { code: 200, data: rows.map(row => ({ abbreviation: row.abbreviation, children: [{ dataset_id: row.dataset_id, children: [{ tissue: row.tissue }] }] })) };
  } else if (path === '/api/cell-data-all/') {
    if (tableFailure) { tableFailure = false; return route.fulfill({ status: 500, json: { detail: 'fixture failure' } }); }
    body = rows.filter(row => (!params.abbreviation || row.abbreviation === params.abbreviation) && (!params.name || Object.values(row).join(' ').toLowerCase().includes(params.name.toLowerCase())));
  } else if (path === '/api/pdf-images/') {
    assert.ok(params.dataset_record_id, 'image request must bind to a dataset record');
    assert.ok(Number(params.page_size) <= 100, 'respect backend page size cap');
    if (imageFailure) { imageFailure = false; return route.fulfill({ status: 500, json: { detail: 'fixture failure' } }); }
    const old = delayOldDataset && params.dataset_record_id === '1';
    if (old) await wait(600);
    const count = params.pdf_name === 'Overview' ? 31 : 2;
    const pageNumber = Number(params.page);
    const pageSize = Number(params.page_size);
    const data = Array.from({ length: Math.min(pageSize, Math.max(0, count - (pageNumber - 1) * pageSize)) }, (_, index) => ({ id: (pageNumber - 1) * pageSize + index + 1, name: params.pdf_name === 'Overview' ? 'Umap_Sample.pdf' : `${params.pdf_name}.pdf`, image: `/fixture-images/${params.dataset_record_id}-${params.pdf_name}-${index}.svg`, image_width: 200 }));
    body = { count, links: { next: pageNumber * pageSize < count ? 'fixture-next' : null, previous: pageNumber > 1 ? 'fixture-prev' : null }, results: { code: 200, data } };
    await route.fulfill({ json: body });
    if (old) completedOldImages += 1;
    return;
  } else if (path === '/api/h5-file/') {
    assert.ok(params.dataset_record_id, 'downloads must bind to a dataset record');
    assert.equal(params.name, undefined, 'missing filename must never request all downloads');
    const old = delayOldDataset && params.dataset_record_id === '1';
    if (old) await wait(600);
    body = ['RNA', 'AUCell', 'UCell', 'singscore'].map((method, index) => ({
      id: Number(params.dataset_record_id) + index * 100,
      name: `dataset-${params.dataset_record_id}.h5`,
      // Exercise both the new method field and legacy relative-path forms.
      method: index === 0 ? method : undefined,
      relative_path: index === 1 ? method : `folder/${method}/dataset-${params.dataset_record_id}.h5`,
      download_url: `/api/h5-file/${Number(params.dataset_record_id) + index * 100}/download/`,
    }));
    await route.fulfill({ json: body });
    if (old) completedOldDownloads += 1;
    return;
  } else if (path.match(/^\/api\/h5-file\/\d+\/download\/$/)) {
    return route.fulfill({ contentType: 'application/octet-stream', headers: { 'Content-Disposition': 'attachment; filename="dataset-2.h5"' }, body: 'fixture h5 payload' });
  } else if (path === '/api/by-site/') {
    body = { code: 200, data: [{ name_en: 'Joint', diseases: [{ abbreviation: 'Gout', name_en: 'Gout' }, { abbreviation: 'JIA', name_en: 'Juvenile idiopathic arthritis' }] }, { name_en: 'Empty category', diseases: [] }] };
  } else {
    return route.fulfill({ status: 404, json: { detail: `Unexpected fixture endpoint ${path}` } });
  }
  return route.fulfill({ json: body });
});
const row = dataset => page.locator('.el-table__body-wrapper tr').filter({ hasText: dataset });
const figure = () => page.locator('.dataset-details .file_item img');
const activeTab = () => page.locator('.dataset-menu .el-menu-item.is-active');
try {
  await page.goto(`${process.env.SCAID_UI_URL || 'http://127.0.0.1:15173'}/Dataset`);
  await row('GSE111').waitFor();
  await row('GSE111').click();
  await until(async () => await figure().count() === 24, 'first figure page did not load');
  assert.match(await page.locator('.dataset-details .result-summary').innerText(), /31 figures.*Page 1 of 2/);
  await page.locator('.el-pagination .btn-next').click();
  await until(async () => await figure().count() === 7, 'second figure page did not load');
  assert.match(await page.locator('.dataset-details .result-summary').innerText(), /Page 2 of 2/);
  await page.getByRole('menuitem', { name: 'Marker', exact: true }).click();
  await until(async () => await figure().count() === 2, 'marker did not load');
  await row('GSE222').click();
  await until(async () => await figure().count() === 24 && (await figure().first().getAttribute('src')).includes('/2-Overview-'), 'dataset switch did not reset Overview');
  assert.equal((await activeTab().innerText()).trim(), 'Overview');

  delayOldDataset = true;
  await row('GSE111').click();
  await until(() => requests.some(item => item.path === '/api/pdf-images/' && item.dataset_record_id === '1' && item.page === '1'), 'older request missing');
  await row('GSE222').click();
  await until(() => completedOldImages === 1 && completedOldDownloads === 1, 'delayed requests did not complete');
  await until(async () => await figure().count() === 24, 'new dataset results missing');
  assert.ok((await figure().first().getAttribute('src')).includes('/2-Overview-'), 'stale image response overwrote the selected dataset');
  assert.match(await page.locator('.selection-title').innerText(), /GSE222/);
  const figureRequestsBeforeDownload = requests.filter(item => item.path === '/api/pdf-images/').length;
  await page.getByRole('button', { name: 'Download (4)', exact: true }).click();
  const expectedLabels = ['RNA', 'AUCell', 'UCell', 'singscore'].map(method => `${method} — dataset-2.h5`);
  await until(async () => await page.locator('.download-link:visible').count() === 4, 'four download methods did not appear');
  assert.deepEqual(await page.locator('.download-link:visible').allTextContents(), expectedLabels, 'same-named H5 matrices must show distinct method labels');
  const downloadLink = page.getByRole('link', { name: 'RNA — dataset-2.h5', exact: true });
  await downloadLink.waitFor();
  assert.equal(await downloadLink.getAttribute('href'), '/api/h5-file/2/download/');
  const downloadPromise = page.waitForEvent('download');
  await downloadLink.click();
  const download = await downloadPromise;
  assert.equal(download.suggestedFilename(), 'dataset-2.h5');
  assert.equal((await activeTab().innerText()).trim(), 'Overview');
  assert.equal(requests.filter(item => item.path === '/api/pdf-images/').length, figureRequestsBeforeDownload, 'download selection issued a figure query');

  imageFailure = true;
  await page.getByRole('menuitem', { name: 'Marker', exact: true }).click();
  await page.getByText('Unable to load figures for this dataset.').waitFor();
  await page.getByRole('button', { name: 'Retry figures' }).click();
  await until(async () => await figure().count() === 2, 'figure retry did not recover');

  await page.getByPlaceholder('Enter a disease abbreviation, dataset ID, tissue or file name').fill('does-not-exist');
  await page.getByRole('button', { name: 'Search', exact: true }).click();
  await until(async () => await page.locator('.dataset-details').count() === 0 && await row('GSE111').count() === 0, 'search retained old details');
  await page.getByRole('button', { name: 'Reset', exact: true }).click();
  await row('GSE111').waitFor();
  tableFailure = true;
  await page.getByRole('button', { name: 'Search', exact: true }).click();
  await page.getByText('Unable to load datasets. Please try Search again.').waitFor();
  await page.getByRole('button', { name: 'Search', exact: true }).click();
  await row('GSE222').waitFor();

  const diseaseSelect = page.locator('.el-select').nth(0);
  const datasetSelect = page.locator('.el-select').nth(1);
  const tissueSelect = page.locator('.el-select').nth(2);
  await diseaseSelect.click();
  await page.getByRole('option', { name: 'JIA', exact: true }).click();
  await datasetSelect.click();
  await page.getByRole('option', { name: 'GSE111', exact: true }).click();
  await tissueSelect.click();
  await page.getByRole('option', { name: 'PBMC/SFMC', exact: true }).click();
  await datasetSelect.hover();
  await datasetSelect.locator('.el-select__clear').click();
  assert.match(await tissueSelect.innerText(), /All tissues/, 'clearing dataset must clear its tissue');
  await datasetSelect.click();
  await page.getByRole('option', { name: 'GSE111', exact: true }).click();
  await diseaseSelect.hover();
  await diseaseSelect.locator('.el-select__clear').click();
  assert.match(await datasetSelect.innerText(), /All datasets/, 'clearing disease must clear its dataset');
  assert.match(await tissueSelect.innerText(), /All tissues/);

  await page.getByRole('link', { name: 'Home', exact: true }).click();
  await page.getByText('No datasets available').waitFor();
  const diseaseLink = page.getByRole('link', { name: 'Juvenile idiopathic arthritis' });
  assert.equal(await diseaseLink.getAttribute('href'), '/Dataset?id=JIA');
  await diseaseLink.click();
  await row('GSE111').waitFor();
  assert.ok(page.url().endsWith('/Dataset?id=JIA'));
  assert.equal(await row('GSE222').count(), 0);
  await page.getByRole('button', { name: 'Reset', exact: true }).click();
  await row('GSE222').waitFor();
  assert.equal(new URL(page.url()).search, '');
  await page.getByRole('link', { name: 'Home', exact: true }).click();
  await page.getByLabel('Search a disease, abbreviation or accession').fill('GSE222');
  await page.getByRole('button', { name: 'Explore', exact: true }).click();
  await row('GSE222').waitFor();
  assert.equal(await row('GSE111').count(), 0);
  assert.equal(new URL(page.url()).searchParams.get('search'), 'GSE222');

  await page.setViewportSize({ width: 390, height: 844 });
  const menuButton = page.getByRole('button', { name: 'Menu', exact: true });
  await menuButton.waitFor();
  assert.equal(await page.getByRole('link', { name: 'Home', exact: true }).isVisible(), false);
  await menuButton.click();
  await page.getByRole('link', { name: 'Home', exact: true }).waitFor();
  const navBox = await page.getByRole('navigation', { name: 'Main navigation' }).boundingBox();
  const titleBox = await page.getByRole('heading', { name: 'Dataset Browser' }).boundingBox();
  assert.ok(titleBox.y >= navBox.y + navBox.height, 'expanded mobile navigation obscures title');
  assert.ok(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), 'mobile layout overflows viewport');
  await page.getByRole('link', { name: 'Home', exact: true }).click();
  await page.getByRole('button', { name: 'Menu', exact: true }).waitFor();
  assert.deepEqual(pageErrors, [], 'uncaught browser errors');
  console.log('PASS: Dataset pagination, tab reset, stale images/downloads, distinct matrix download labels and real download link, errors/retry, search/reset and filter clearing, JIA navigation, home search and mobile menu.');
} finally {
  await browser.close();
}

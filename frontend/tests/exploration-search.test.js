import test from 'node:test'
import assert from 'node:assert/strict'
import { createExplorationSearch } from '../src/utils/explorationSearch.js'
import { getFullImageUrl, getApiUrl } from '../src/config/index.js'
import http from '../src/utils/http.js'

const pageResponse = (data, count = data.length) => ({
  count,
  links: { next: 'http://47.115.146.132/api/pdf-images/?page=2' },
  results: { code: 200, data }
})
const pending = () => {
  let resolve, reject
  const promise = new Promise((ok, fail) => { resolve = ok; reject = fail })
  return { promise, resolve, reject }
}

test('blank Gene and KEGG searches show validation and never request the whole image collection', async () => {
  for (const kind of ['gene', 'kegg']) {
    let calls = 0
    const search = createExplorationSearch({ kind, request: async () => { calls += 1 } })
    await search.submit({ pdf_name: '   ' })
    assert.equal(calls, 0)
    assert.match(search.state.validationError, /Please select/)
    assert.equal(search.state.searched, false)
  }
})

test('all 48 matching plots are reachable, with relative page requests and preserved filters', async () => {
  const requests = []
  const search = createExplorationSearch({ kind: 'gene', request: async (url, params) => {
    requests.push({ url, params })
    const first = (params.page - 1) * params.page_size
    return pageResponse(Array.from({ length: Math.min(params.page_size, 48 - first) }, (_, i) => ({ id: first + i + 1 })), 48)
  } })
  const filters = { pdf_name: 'CD3D', abbreviation: 'AS', dataset_id: 'GSE216883' }
  await search.submit(filters)
  assert.equal(search.state.total, 48)
  assert.equal(search.state.files.length, 10)
  const seen = new Set(search.state.files.map(item => item.id))
  filters.pdf_name = 'Changed form outside controller'
  for (let page = 2; page <= 5; page += 1) {
    await search.changePage(page)
    search.state.files.forEach(item => seen.add(item.id))
  }
  assert.equal(seen.size, 48)
  assert.equal(search.state.files.length, 8)
  assert.ok(requests.every(({ url, params }) => url === '/pdf-images/' && params.pdf_name === 'CD3D' && params.abbreviation === 'AS' && params.dataset_id === 'GSE216883'))
  await search.changePageSize(999)
  assert.equal(search.state.pageSize, 100)
  assert.equal(search.state.page, 1)
  assert.equal(search.state.files.length, 48)
})

test('KEGG requests select the KEGGScore family and retain the selected scoring method', async () => {
  let params
  const search = createExplorationSearch({ kind: 'kegg', request: async (_, value) => {
    params = value
    return pageResponse([{ id: 1 }])
  } })
  await search.submit({ pdf_name: ' hsa00240 ', kegg_param: 'UCell' })
  assert.equal(params.pdf_name, 'hsa00240')
  assert.equal(params.KEGGScore, 'KEGGScore')
  assert.equal(params.kegg_param, 'UCell')
  assert.equal('KEGGUmap' in params, false)
})

test('restoring a shared KEGG query requests its selected page and page size without fetching page one first', async () => {
  const calls = []
  const search = createExplorationSearch({ kind: 'kegg', request: async (_, params) => {
    calls.push(params)
    return pageResponse([{ id: 21 }], 48)
  } })
  search.state.pageSize = 20
  await search.submit({ pdf_name: 'hsa00010', abbreviation: 'BD', dataset_id: 'GSE198616', kegg_param: 'UCell' }, 2)
  assert.equal(calls.length, 1)
  assert.equal(calls[0].page, 2)
  assert.equal(calls[0].page_size, 20)
  assert.equal(calls[0].abbreviation, 'BD')
  assert.equal(search.state.page, 2)
  assert.equal(search.state.total, 48)
})

test('a slower earlier result cannot replace the newly selected gene or stop its loading state', async () => {
  const old = pending(), latest = pending()
  const signals = []
  const search = createExplorationSearch({ kind: 'gene', request: (_, params, options) => {
    signals.push(options.signal)
    return params.pdf_name === 'CD3D' ? old.promise : latest.promise
  } })
  const oldTask = search.submit({ pdf_name: 'CD3D' })
  const latestTask = search.submit({ pdf_name: 'MS4A1' })
  assert.equal(signals[0].aborted, true)
  old.resolve(pageResponse([{ id: 'wrong-gene' }]))
  await oldTask
  assert.equal(search.state.loading, true)
  assert.deepEqual(search.state.files, [])
  latest.resolve(pageResponse([{ id: 'latest-gene' }]))
  await latestTask
  assert.equal(search.state.loading, false)
  assert.deepEqual(search.state.files, [{ id: 'latest-gene' }])
})

test('late errors do not obscure newer successful results; filter changes invalidate active requests', async () => {
  const old = pending(), afterEdit = pending()
  const search = createExplorationSearch({ kind: 'gene', request: (_, params) => {
    if (params.pdf_name === 'CD3D') return old.promise
    if (params.pdf_name === 'CD4') return afterEdit.promise
    return Promise.resolve(pageResponse([{ id: 'latest' }]))
  } })
  const oldTask = search.submit({ pdf_name: 'CD3D' })
  await search.submit({ pdf_name: 'MS4A1' })
  old.reject(new Error('Older request failed'))
  await oldTask
  assert.equal(search.state.error, '')
  assert.deepEqual(search.state.files, [{ id: 'latest' }])
  const editTask = search.submit({ pdf_name: 'CD4' })
  search.resetResults()
  afterEdit.resolve(pageResponse([{ id: 'stale-after-edit' }]))
  await editTask
  assert.equal(search.state.searched, false)
  assert.deepEqual(search.state.files, [])
})

test('an unsuccessful search clears old plots and distinguishes timeout, failure and an empty result', async () => {
  let response = pageResponse([{ id: 'previous' }])
  const search = createExplorationSearch({ kind: 'gene', request: async () => {
    if (response instanceof Error) throw response
    return response
  } })
  await search.submit({ pdf_name: 'CD3D' })
  response = Object.assign(new Error('Timeout'), { code: 'ECONNABORTED' })
  await search.submit({ pdf_name: 'MS4A1' })
  assert.deepEqual(search.state.files, [])
  assert.match(search.state.error, /timed out/)
  assert.equal(search.state.loading, false)
  response = { results: { code: 500, message: 'Service error' } }
  await search.submit({ pdf_name: 'MS4A1' })
  assert.match(search.state.error, /could not be loaded/)
  response = pageResponse([])
  await search.submit({ pdf_name: 'MS4A1' })
  assert.equal(search.state.error, '')
  assert.equal(search.state.searched, true)
  assert.equal(search.state.total, 0)
})

test('suggestions retain all pages beyond the backend 100-item cap and normalize both types', async () => {
  for (const kind of ['gene', 'kegg']) {
    const search = createExplorationSearch({ kind, request: async (url, params) => {
      assert.equal(url, '/gene-kegg-enum/')
      assert.equal(params.tp, kind)
      assert.equal(params.page_size, 100)
      const first = (params.page - 1) * 100
      return pageResponse(Array.from({ length: Math.min(100, 125 - first) }, (_, i) => ({ id: first + i })), 125)
    } })
    await search.fetchOptions('CD')
    assert.equal(search.state.options.length, 100)
    assert.equal(search.state.optionsTotal, 125)
    await search.loadMoreOptions()
    assert.equal(search.state.options.length, 125)
    assert.equal(new Set(search.state.options.map(item => item.id)).size, 125)
  }
})

test('suggestion failures are visible and retryable; clearing the input discards pending suggestions', async () => {
  const delayed = pending()
  let fail = true
  const search = createExplorationSearch({ kind: 'gene', request: async (_, params) => {
    if (params.name === 'MS4') return delayed.promise
    if (fail) throw new Error('Backend unavailable')
    return pageResponse([{ id: 1, name: 'CD3D' }])
  } })
  await search.fetchOptions('CD3')
  assert.match(search.state.optionsError, /could not be loaded/)
  assert.equal(search.state.optionsLoading, false)
  fail = false
  await search.retryOptions()
  assert.equal(search.state.optionsError, '')
  assert.equal(search.state.options[0].name, 'CD3D')
  const task = search.fetchOptions('MS4')
  assert.equal(search.state.optionsLoading, true)
  await search.fetchOptions('')
  delayed.resolve(pageResponse([{ id: 2, name: 'MS4A1' }]))
  await task
  assert.deepEqual(search.state.options, [])
  assert.equal(search.state.optionsLoading, false)
})

test('image and API URLs preserve the public origin and unrelated external image hosts', () => {
  assert.equal(getFullImageUrl('/media/pdf_images/a.jpg'), '/system/media/pdf_images/a.jpg')
  assert.equal(getFullImageUrl('pdf_images/a.jpg'), '/system/media/pdf_images/a.jpg')
  assert.equal(getFullImageUrl('/system/media/pdf_images/a.jpg'), '/system/media/pdf_images/a.jpg')
  assert.equal(getFullImageUrl('/api/pdf-images/130/original/'), '/api/pdf-images/130/original/')
  assert.equal(getFullImageUrl('http://47.115.146.132:10209/system/media/a.jpg?x=1'), '/system/media/a.jpg?x=1')
  assert.equal(getFullImageUrl('http://47.115.146.132/media/a.jpg'), '/system/media/a.jpg')
  assert.equal(getFullImageUrl('https://example.org/a.jpg'), 'https://example.org/a.jpg')
  assert.equal(getApiUrl('/api/pdf-images/'), '/api/pdf-images/')
  assert.equal(getApiUrl('/pdf-images/'), '/api/pdf-images/')
})

test('HTTP calls use one API prefix, a 30-second timeout and accept both pagination envelopes', async () => {
  for (const url of ['/pdf-images/', '/api/pdf-images/']) {
    for (const body of [pageResponse([{ id: 1 }]), { count: 1, results: [{ id: 1 }] }]) {
      const response = await http({ url, adapter: async (config) => {
        assert.equal(http.getUri(config), '/api/pdf-images/')
        assert.equal(config.timeout, 30000)
        return { data: body, status: 200, statusText: 'OK', headers: {}, config }
      } })
      assert.equal(response, body)
    }
  }
})

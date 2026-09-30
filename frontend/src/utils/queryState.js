// Only strings and positive, bounded integers become public search parameters.
// These helpers keep URL restoration and UI submissions on the same contract.
export const queryText = (value) => typeof value === 'string' ? value.trim() : ''
export const queryPage = (value, fallback = 1) => {
  const number = Number(value)
  return Number.isSafeInteger(number) && number > 0 && number <= 100000 ? number : fallback
}
export const compactQuery = (value) => Object.fromEntries(Object.entries(value)
  .filter(([, item]) => item !== '' && item !== undefined && item !== null)
  .map(([key, item]) => [key, String(item)]))
export const querySignature = (value) => JSON.stringify(Object.entries(value).sort(([a], [b]) => a.localeCompare(b)))

export const readExplorationQuery = (query, kind) => ({
  pdf_name: queryText(query[kind]),
  abbreviation: queryText(query.disease),
  dataset_id: queryText(query.dataset),
  kegg_param: kind === 'kegg' && ['AUCell', 'UCell', 'singscore'].includes(query.method) ? query.method : '',
  page: queryPage(query.page),
  pageSize: [10, 20, 50, 100].includes(Number(query.size)) ? Number(query.size) : 10
})
export const explorationQuery = (filters, kind, page = 1, pageSize = 10) => compactQuery({
  [kind]: queryText(filters.pdf_name),
  disease: queryText(filters.abbreviation),
  dataset: queryText(filters.dataset_id),
  method: kind === 'kegg' ? queryText(filters.kegg_param) : '',
  page: page > 1 ? page : '',
  size: pageSize !== 10 ? pageSize : ''
})

export const readDatasetQuery = (query, analyses) => ({
  filters: {
    abbreviation: queryText(query.id),
    dataset_id: queryText(query.dataset),
    tissue: queryText(query.tissue),
    name: queryText(query.search)
  },
  record: /^\d+$/.test(queryText(query.record)) ? queryText(query.record) : '',
  analysis: analyses.includes(query.analysis) ? query.analysis : 'Overview',
  page: queryPage(query.page)
})
export const datasetQuery = (filters, record = '', analysis = 'Overview', page = 1) => compactQuery({
  id: queryText(filters.abbreviation),
  dataset: queryText(filters.dataset_id),
  tissue: queryText(filters.tissue),
  search: queryText(filters.name),
  record,
  analysis: record && analysis !== 'Overview' ? analysis : '',
  page: record && page > 1 ? page : ''
})

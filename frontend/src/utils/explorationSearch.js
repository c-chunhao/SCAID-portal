export const createExplorationState = () => ({
  files: [],
  total: 0,
  page: 1,
  pageSize: 10,
  loading: false,
  searched: false,
  error: '',
  validationError: '',
  options: [],
  optionsTotal: 0,
  optionsPage: 0,
  optionsLoading: false,
  optionsError: ''
})

export const requestErrorMessage = (error, subject) => {
  if (error?.code === 'ECONNABORTED' || error?.code === 'ETIMEDOUT') {
    return `${subject} timed out. Please retry or narrow your search.`
  }
  return `${subject} could not be loaded. Please try again.`
}

const readPage = (response) => {
  const envelope = response?.results ?? response
  if (envelope?.code !== undefined && Number(envelope.code) !== 200) {
    throw new Error(envelope.message || 'Request failed')
  }
  const data = Array.isArray(envelope) ? envelope : envelope?.data
  if (!Array.isArray(data)) throw new Error('Invalid response')
  return { data, count: Number(response.count ?? data.length) }
}

// Accept a reactive state in Vue, or a plain state in regression tests.
export const createExplorationSearch = ({ request, kind, state = createExplorationState() }) => {
  let resultVersion = 0
  let optionVersion = 0
  let resultController
  let optionController
  let submittedQuery = null
  let optionQuery = ''

  const resetResults = () => {
    resultVersion += 1
    resultController?.abort()
    submittedQuery = null
    Object.assign(state, { files: [], total: 0, page: 1, loading: false, searched: false, error: '', validationError: '' })
  }

  const fetchResults = async (page = 1) => {
    if (!submittedQuery) return
    const version = ++resultVersion
    resultController?.abort()
    resultController = new AbortController()
    Object.assign(state, { files: [], page, loading: true, searched: true, error: '' })
    try {
      const response = await request('/pdf-images/', {
        ...submittedQuery,
        page,
        page_size: state.pageSize
      }, { signal: resultController.signal })
      if (version !== resultVersion) return
      const { data, count } = readPage(response)
      state.files = data
      state.total = count
    } catch (error) {
      if (version === resultVersion) state.error = requestErrorMessage(error, 'Search results')
    } finally {
      if (version === resultVersion) state.loading = false
    }
  }

  const submit = async (filters, page = 1) => {
    resetResults()
    const name = String(filters.pdf_name || '').trim()
    if (!name) {
      state.validationError = `Please select a ${kind === 'gene' ? 'gene' : 'KEGG pathway'} before searching.`
      return
    }
    submittedQuery = {
      format: 'json',
      pdf_name: name,
      abbreviation: filters.abbreviation || '',
      dataset_id: filters.dataset_id || '',
      ...(kind === 'gene' ? { GeneUmap: 'GeneUmap' } : { KEGGScore: 'KEGGScore', kegg_param: filters.kegg_param || '' })
    }
    await fetchResults(page)
  }

  const changePageSize = async (size) => {
    state.pageSize = Math.max(1, Math.min(100, Number(size) || 10))
    await fetchResults(1)
  }

  const fetchOptions = async (query, append = false) => {
    if (append && (state.optionsLoading || state.options.length >= state.optionsTotal)) return
    const version = ++optionVersion
    optionController?.abort()
    optionController = new AbortController()
    optionQuery = String(query || '').trim()
    state.optionsError = ''
    if (!append) {
      state.options = []
      state.optionsTotal = 0
      state.optionsPage = 0
    }
    if (!optionQuery) {
      state.optionsLoading = false
      return
    }
    state.optionsLoading = true
    const page = append ? state.optionsPage + 1 : 1
    try {
      const response = await request('/gene-kegg-enum/', {
        format: 'json', tp: kind, name: optionQuery, page_size: 100, page
      }, { signal: optionController.signal })
      if (version !== optionVersion) return
      const { data, count } = readPage(response)
      state.options = append ? [...state.options, ...data] : data
      state.optionsTotal = count
      state.optionsPage = page
    } catch (error) {
      if (version === optionVersion) state.optionsError = requestErrorMessage(error, 'Search suggestions')
    } finally {
      if (version === optionVersion) state.optionsLoading = false
    }
  }

  const dispose = () => {
    resetResults()
    optionVersion += 1
    optionController?.abort()
  }

  return {
    state, submit, resetResults, changePage: fetchResults, changePageSize,
    fetchOptions, loadMoreOptions: () => fetchOptions(optionQuery, true),
    retryOptions: () => fetchOptions(optionQuery), dispose
  }
}

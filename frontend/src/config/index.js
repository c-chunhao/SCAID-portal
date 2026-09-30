export const API_CONFIG = {
  BASE_URL: '/system',
  DOWN_URL: '',
  API_PREFIX: '/api'
}

// Keep requests on the browser's origin, including its scheme and proxy port.
export const normalizeApiPath = (endpoint) => {
  return `/${String(endpoint || '').replace(/^\/+/, '').replace(/^(?:api\/)+/, '')}`
}

export const getApiUrl = (endpoint) => `${API_CONFIG.API_PREFIX}${normalizeApiPath(endpoint)}`

export const getFullImageUrl = (imageUrl) => {
  if (!imageUrl) return imageUrl
  let path = String(imageUrl)
  if (/^(?:https?:)?\/\//i.test(path)) {
    const url = new URL(path, 'https://localhost')
    // Old records contain this deployment's absolute URL. Preserve unrelated hosts.
    if (url.hostname !== '47.115.146.132') return path
    path = `${url.pathname}${url.search}${url.hash}`
  }
  path = `/${path.replace(/^\/+/, '')}`
  if (/^\/api\/pdf-images\/\d+\/original\/(?:[?#]|$)/.test(path)) return path
  if (path.startsWith('/system/')) return path
  return `${API_CONFIG.BASE_URL}${path.startsWith('/media/') ? path : `/media${path}`}`
}

import axios from 'axios'
import { normalizeApiPath } from '../config/index.js'

const service = axios.create({
  baseURL: '/api',
  timeout: 30000
})

service.interceptors.request.use((config) => {
  // Accept older callers that include /api without requesting /api/api/.
  if (!/^(?:https?:)?\/\//i.test(config.url || '')) {
    config.url = normalizeApiPath(config.url)
  }
  // The public API is read-only and unauthenticated; no token is ever attached.
  return config
})

service.interceptors.response.use((response) => {
  const body = response.data
  const envelope = body?.results && !Array.isArray(body.results) ? body.results : body
  if (envelope?.code !== undefined && Number(envelope.code) !== 200) {
    return Promise.reject(new Error(envelope.message || 'The request failed. Please try again.'))
  }
  return body
})

export default service

import request from '@/utils/http'

export const getRequest = (url, params, options = {}) => request({
  ...options,
  url,
  method: 'get',
  params
})

export const postRequest = (url, data) => request({ url, method: 'post', data })
export const putRequest = (url, data) => request({ url, method: 'put', data })
export const deleteRequest = (url, params) => request({ url, method: 'delete', params })

import test from 'node:test'
import assert from 'node:assert/strict'
import { datasetQuery, explorationQuery, readDatasetQuery, readExplorationQuery } from '../src/utils/queryState.js'

test('a copied dataset query restores record, analysis, tissue with a slash, filters and page', () => {
  const query = datasetQuery({ abbreviation: 'JIA', dataset_id: 'GSE205095', tissue: 'PBMC/SFMC', name: 'T cell' }, 130, 'LRcircle', 2)
  const encoded = new URLSearchParams(query)
  const restored = readDatasetQuery(Object.fromEntries(encoded), ['Overview', 'LRcircle'])
  assert.deepEqual(restored, { filters: { abbreviation: 'JIA', dataset_id: 'GSE205095', tissue: 'PBMC/SFMC', name: 'T cell' }, record: '130', analysis: 'LRcircle', page: 2 })
})

test('a copied KEGG query restores all optional filters, tissue, method, page and size', () => {
  const query = explorationQuery({ pdf_name: 'hsa00010', abbreviation: 'PsA', dataset_id: '8207', tissue: 'SFMC', kegg_param: 'AUCell' }, 'kegg', 3, 20)
  assert.deepEqual(readExplorationQuery(query, 'kegg'), { pdf_name: 'hsa00010', abbreviation: 'PsA', dataset_id: '8207', tissue: 'SFMC', kegg_param: 'AUCell', page: 3, pageSize: 20 })
})

test('Gene links preserve tissue with a slash and old links leave tissue unrestricted', () => {
  const query = explorationQuery({ pdf_name: ' CD3D ', abbreviation: 'AS', dataset_id: 'GSE216883', tissue: 'PBMC/SFMC' }, 'gene', 2, 50)
  const restored = readExplorationQuery(Object.fromEntries(new URLSearchParams(query)), 'gene')
  assert.deepEqual(restored, { pdf_name: 'CD3D', abbreviation: 'AS', dataset_id: 'GSE216883', tissue: 'PBMC/SFMC', kegg_param: '', page: 2, pageSize: 50 })
  assert.deepEqual(explorationQuery({ pdf_name: 'CD3D' }, 'gene'), { gene: 'CD3D' })
  assert.equal(readExplorationQuery({ gene: 'CD3D' }, 'gene').tissue, '')
})

test('malformed or repeated query values fall back to safe defaults', () => {
  assert.deepEqual(readExplorationQuery({ gene: ['CD3D', 'TRIM21'], disease: ['RA'], tissue: ['PBMC', 'SFMC'], page: '-1', size: '999', method: 'invalid' }, 'gene'), { pdf_name: '', abbreviation: '', dataset_id: '', tissue: '', kegg_param: '', page: 1, pageSize: 10 })
  assert.equal(readDatasetQuery({ record: 'a130', analysis: 'bad', page: 'Infinity' }, ['Overview']).record, '')
  assert.equal(readDatasetQuery({ page: '1.5' }, ['Overview']).page, 1)
})
